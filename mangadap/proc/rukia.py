# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Implements a full-spectrum fitting algorithm that fits a
wavelength-dependent dispersion.

----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

import os
import time
import numpy

from astropy.io import fits
import astropy.constants
from matplotlib import pyplot, rc, patches
from matplotlib.ticker import NullFormatter, MultipleLocator, FormatStrFormatter
from scipy import interpolate, optimize, ndimage

from ppxf.ppxf_util import gaussian_filter1d

#from mangadap.util.instrument import convolution_variable_sigma, spectrum_velocity_scale
#from mangadap.util.instrument import spectral_resolution, resample_vector_npix, resample_vector
#from mangadap.util.constants import constants
#from mangadap.util.filter import off_diagonal_identity, BoxcarFilter, interpolate_masked_vector

from mangadap.util.sampling import angstroms_per_pixel, Resample
from mangadap.util.bitmask import BitMask

#from mangadap.proc.templatelibrary import TemplateLibrary, TemplateLibraryDef
#from mangadap.proc.util import optimal_scale

#-----------------------------------------------------------------------------

class PieceWiseLinear:
    def __init__(self, deg, c=None, rng=[-1,1]):
        self.deg = deg
        self.np = self.deg
        self.coeff = numpy.ones(self.np, dtype=float) if c is None else c
        self.rng = rng

    def reset_range(self, rng):
        self.rng = rng

    def __call__(self, x, a=None):
        if a is not None:
            if len(a) != self.np:
                raise ValueError('Incorrect number of parameters')
            self.coeff = a
        if self.deg == 1:
            return self.coeff[0]
        return interpolate.interp1d(numpy.linspace(*self.rng, self.deg), self.coeff,
                                    fill_value='extrapolate', assume_sorted=True)(x)


class LegendrePolynomial:
    def __init__(self, deg, c=None):
        self.deg = deg
        self.np = self.deg+1
        self.coeff = numpy.ones(self.np, dtype=float) if c is None else c

    def __call__(self, x, a=None):
        if a is not None:
            if len(a) != self.np:
                raise ValueError('Incorrect number of parameters')
            self.coeff = a
        return numpy.polynomial.legendre.legval(x, self.coeff)


class ScaledEmpirical:
    def __init__(self, x, y, scale_func=None):
        self.deg = 0 if scale_func is None else scale_func.deg
        self.np = 1 if scale_func is None else scale_func.np
        self.sf = LegendrePolynomial(0) if scale_func is None else scale_func
        self.coeff = numpy.array([1.0]) if scale_func is None else scale_func.coeff
        self.interp = interpolate.interp1d(x, y, fill_value='extrapolate')
        self.has_rng = hasattr(self.sf, 'rng')

    def reset_rng(self, rng):
        if self.has_rng:
            self.sf.reset_rng(rng)

    def __call__(self, x, a=None):
        if a is not None:
            if len(a) != self.np:
                raise ValueError('Incorrect number of parameters')
            self.coeff = a
        return self.sf(x, a=self.coeff)*self.interp(x)
        

class ModelSpectrum:
    """
    The domain for all functions is expected to be [-1, 1].

    If sig_func is None, only a shift is applied.

    """
    def __init__(self, tpl_flux, sig_func=None, add_func=None, mult_func=None, sampling_factor=1):
        if len(tpl_flux.shape) != 1:
            raise ValueError('Use vectors, not arrays.')
        if not isinstance(sampling_factor,int):
            raise TypeError('sampling_factor must be an integer; can be 1')
        if len(tpl_flux) % sampling_factor:
            raise ValueError('The length of the flux vector must be an integer number of samples.')

        # Number of pixels to bin for output
        self.sampling_factor = sampling_factor
        # Total number of output pixels (expects no modulus as tested
        # above).
        self.npix = len(tpl_flux)//self.sampling_factor
        # Rebinned samples
        self.x = numpy.linspace(-1., 1., self.npix)
        # Native samples
        self.s = numpy.linspace(-1., 1., self.npix*self.sampling_factor)

        # Set the functions for sigma, the additive polynomial, and the
        # multiplicative polynomial
        self.sf = sig_func
        self.af = add_func
        self.mf = mult_func
        
        self.sd = 0 if self.sf is None else self.sf.np
        self.ad = 0 if self.af is None else self.af.np
        self.md = 0 if self.mf is None else self.mf.np

        self.flux = tpl_flux.copy()
        self.np = 1+self.sd+self.ad+self.md     # all the functional coeffs and the pixel shift
        self.p = numpy.zeros(self.np)
        self.free = numpy.ones(self.np, dtype=bool)
        self.nfree = numpy.sum(self.free)
        self.pname = self._set_parameter_names()


    def _set_parameter_names(self):
        pname = [ 'SHIFT' ]
        for i in range(self.sd):
            pname += [ 'SIG{0}'.format(i+1) ]
        for i in range(self.ad):
            pname += [ 'ADD{0}'.format(i+1) ]
        for i in range(self.md):
            pname += [ 'MULT{0}'.format(i+1) ]
        return pname


    def __call__(self, a):
        return self.get_model(a)

    def _set_par(self, a):
        if a.shape != self.nfree:
            raise ValueError('Number of parameters is not correct!')
        self.p[self.free] = a

    def _sig_coeff(self):
        return self.p[1:1+self.sd]

    def _add_coeff(self):
        return self.p[1+self.sd:1+self.sd+self.ad]

    def _mult_coeff(self):
        return self.p[1+self.sd+self.ad:]

    def fix_par(self, p, fixed):
        if len(fixed) != len(self.free):
            raise ValueError('Fixed vector has incorrect length.')
        if len(p.shape) != 1 or p.size != self.np:
            raise ValueError('Par vector has incorrect shape.')
        self.p = p.copy()
        self.free = ~fixed
        self.nfree = numpy.sum(self.free)

    def get_model(self, a):
        self._set_par(a)
        model = ndimage.interpolation.shift(self.flux, self.p[0]*self.sampling_factor, order=1)
        if self.sf is not None:
            model = convolution_variable_sigma(model,
                            self.sampling_factor*self.sf(self.s, a=self._sig_coeff()))
        if self.sampling_factor > 1:
            model = numpy.mean(model.reshape(-1,self.sampling_factor), axis=1)
        if self.mf is not None:
            model *= self.mf(self.x, a=self._mult_coeff())
        if self.af is not None:
            model += self.af(self.x, a=self._add_coeff())
        return model

    def get_sigma_function(self, a, sample=False):
        if self.sf is None:
            return None
        self._set_par(a)
        return (self.sampling_factor if sample else 1.0) \
                    * self.sf(self.s if sample else self.x, a=self._sig_coeff())

    def get_multc_function(self, a, sample=False):
        if self.mf is None:
            return None
        self._set_par(a)
        return self.mf(self.s if sample else self.x, a=self._mult_coeff())

    def get_addc_function(self, a, sample=False):
        if self.af is None:
            return None
        self._set_par(a)
        return self.af(self.s if sample else self.x, a=self._add_coeff())

    def bounds(self):
        lb = numpy.full(self.np, -numpy.inf, dtype=float)
        ub = numpy.full(self.np, numpy.inf, dtype=float)

        # +/- 5 pixels for offset
        lb[0] = -10
        ub[0] = 10

        # sigma factor must be positive
        lb[1] = 0.0
        # multiplicative factor must be positive
        lb[1+self.sd+self.ad] = 0.0
        return (lb[self.free], ub[self.free])
        


class FitResiduals:
    def __init__(self, obj_flux, model, goodpixels=None):
        self.obj_flux = obj_flux.copy()
        self.model = model
        self.gpm = numpy.arange(obj_flux.size) if goodpixels is None else goodpixels
    def __call__(self, a):
        return (self.obj_flux-self.model(a))[self.gpm]


class FitChiSquare:
    def __init__(self, obj_flux, obj_ferr, model, goodpixels=None):
        self.obj_flux = obj_flux.copy()
        self.obj_ferr = obj_ferr.copy()
        self.model = model
        self.gpm = numpy.arange(obj_flux.size) if goodpixels is None else goodpixels
    def __call__(self, a):
        return ((self.obj_flux[self.gpm]-self.model(a)[self.gpm])/self.obj_ferr[self.gpm])


def template_library_list():

    return [ TemplateLibraryDef(key='BASS10k',
                                file_search='./bass2000_solar_spectrum_10k.fits',
                                fwhm=None,
                                sres_ext='SPECRES',
                                in_vacuum=False,
                                wave_limit=numpy.array([None, None]),
                                lower_flux_limit=None,
                                log10=False),

             TemplateLibraryDef(key='BASS20k',
                                file_search='./bass2000_solar_spectrum_20k.fits',
                                fwhm=None,
                                sres_ext='SPECRES',
                                in_vacuum=False,
                                wave_limit=numpy.array([None, None]),
                                lower_flux_limit=None,
                                log10=False),

             TemplateLibraryDef(key='Kurucz10k',
                                file_search='./kurucz_solar_spectrum_10k.fits',
                                fwhm=None,
                                sres_ext='SPECRES',
                                in_vacuum=True,
                                wave_limit=numpy.array([None, None]),
                                lower_flux_limit=None,
                                log10=False),

             TemplateLibraryDef(key='Kurucz20k',
                                file_search='./kurucz_solar_spectrum_20k.fits',
                                fwhm=None,
                                sres_ext='SPECRES',
                                in_vacuum=True,
                                wave_limit=numpy.array([None, None]),
                                lower_flux_limit=None,
                                log10=False)
            ]


def resample_twilight_spectra(wave, flux, ferr, sres, mask=None, dlogLam=1e-5,
                              wavelength_range=None):

    if not isinstance(flux, numpy.ma.MaskedArray):
        raise TypeError('flux must be a masked array')
    if not isinstance(ferr, numpy.ma.MaskedArray):
        raise TypeError('ferr must be a masked array')

    _mask = numpy.zeros(flux.shape, dtype=int) if mask is None else mask.astype(int)
    _mask |= numpy.ma.getmaskarray(flux).astype(int)

    nspec = flux.shape[0]

    fullRange = [wave[0], wave[-1]] if wavelength_range is None \
                    else numpy.array(wavelength_range).astype(float)

    # Get the number of pixels needed
    npix, _fullRange = resample_vector_npix(outRange=fullRange, dx=dlogLam, log=True)

    # Now resample the spectra.  First allocate the arrays
    flux_out = numpy.zeros((nspec, npix), dtype=numpy.float64)
    ferr_out = numpy.zeros((nspec, npix), dtype=numpy.float64)
    sres_out = numpy.zeros((nspec, npix), dtype=numpy.float64)
    mask_out = numpy.zeros((nspec, npix), dtype=numpy.float64)
    wave_in = wave.copy()
    for i in range(nspec):
        print('Resampling {0}/{1}'.format(i+1,nspec), end='\r')
        # Rebin the observed wavelength range
        wave_out, flux_out[i, :] = resample_vector(flux.data[i,:], xRange=[wave_in[0], wave_in[-1]],
                                                   inLog=True, newRange=fullRange, newLog=True,
                                                   dx=dlogLam, ext_value=-9999., conserve=False,
                                                   flat=False)

        wave_out, ferr_out[i, :] = resample_vector(numpy.square(ferr.data[i,:]), xRange=[wave_in[0],
                                                   wave_in[-1]], inLog=True, newRange=fullRange,
                                                   newLog=True, dx=dlogLam, ext_value=-9999.,
                                                   conserve=False, flat=False)

        wave_out, mask_out[i, :] = resample_vector(_mask[i,:].astype(float), xRange=[wave_in[0],
                                                   wave_in[-1]], inLog=True, newRange=fullRange,
                                                   newLog=True, dx=dlogLam, ext_value=1.0,
                                                   conserve=False, flat=False)

#        pyplot.step(wave_in, flux[i,:], where='mid')
#        pyplot.step(wave_out, flux_out[i,:], where='mid', color='g')
#        pyplot.show()

        # Find the unobserved pixels, set them to have 0. flux, and
        # flag them as having no data
        indx = numpy.where(flux_out[i,:] < -9900.)
        flux_out[i,indx] = 0.0
        ferr_out[i,indx] = 1.0

        ferr_out[i,:] = numpy.sqrt(ferr_out[i,:])

        # Resample the spectral resolution by simple interpolation.
        # Select the good pixels
        # define the interpolator (uses linear interpolation; k=1)
#        interpolator = InterpolatedUnivariateSpline(self.hdu['WAVE'].data[i,:].ravel(),
#                                                    self.hdu['SPECRES'].data[i,:].ravel(), k=1)
        interpolator = interpolate.interp1d(wave_in, sres[i,:], fill_value='extrapolate')
        # And then interpolate
        sres_out[i,:] = interpolator(wave_out)

#        pyplot.step(wave_in, sres[i,:], where='mid')
#        pyplot.step(wave_out, sres_out[i,:], where='mid', color='g')
#        pyplot.show()

    print('Resampling DONE                     ')
    indx = mask_out > 0.1
    mask_out[indx] = 1.0
    mask_out[~indx] = 0.0

#    print(ferr_out.size)
#    print(numpy.sum(numpy.invert(ferr_out>0)))
#    print(numpy.sum(mask_out))
#    print(numpy.sum(numpy.invert(ferr_out>0) & numpy.invert(mask_out.astype(bool))))

    indx = numpy.invert(ferr_out > 0)
    mask_out[indx] = 1.0

#    print(ferr_out.size)
#    print(numpy.sum(numpy.invert(ferr_out>0)))
#    print(numpy.sum(mask_out))
#    print(numpy.sum(numpy.invert(ferr_out>0) & numpy.invert(mask_out.astype(bool))))
#
#    exit()

    return wave_out, flux_out, ferr_out, sres_out, mask_out.astype(bool)


def calculate_parameter_covariance(result):
    a = numpy.dot(result.jac.T,result.jac)
    try:
        return numpy.linalg.inv(a)
    except:
        return numpy.linalg.pinv(a)


def initialize_model(wave, tpl, sampling_factor=10, sig_degree=None, sig_pwl=False,
                     add_degree=None, mult_degree=None, base_mult=None, fom_is_chisqr=False,
                     base_sig=None, fix_sig_diff=False, init_sig_factor=1.0):

    # Build the fitting functions:
    # Convolution kernel sigma:
    sigma_function = (PieceWiseLinear(sig_degree) if sig_pwl else LegendrePolynomial(sig_degree)) \
                        if base_sig is None else ScaledEmpirical(base_sig)
    # Multiplicative Continuum:
    mult_function = None if mult_degree is None else LegendrePolynomial(mult_degree)
    if base_mult is not None:
        mult_function = ScaledEmpirical(base_mult, scale_func=mult_function)
    # Additive Continuum:
    add_function = None if add_degree is None else LegendrePolynomial(add_degree)

    # Instantiate the class that produces the model spectrum
    model = ModelSpectrum(tpl['FLUX'].data[0,:], sig_func=sigma_function, add_func=add_function,
                          mult_func=mult_function, sampling_factor=sampling_factor)

    # Initialize the kernel to a constant offset based on the spectral
    # resolution vectors
    p = numpy.zeros(model.np)
    if base_sig is None:
        p[1] = 1.0 #init_sig_factor * numpy.ma.sqrt(numpy.ma.mean(numpy.square(base_sig)))
    else:
        p[1] = init_sig_factor
        if fix_sig_diff:
            fixp = numpy.zeros(model.np, dtype=bool)
            fixp[1] = True
            model.fix_par(p, fixp)

    # Initialize the multiplicative polynomial to 1
    if mult_degree is not None:
        p[1+model.sd+model.ad] = 1.0
    # The parameter vector that will be optimized
    a0 = p[model.free]

    return model, p, a0


def fit_spectra(wave, flux, ferr, sres, good_pixels, tpl, sampling_factor=None, sig_degree=None,
                sig_pwl=False, add_degree=None, mult_degree=None, base_mult=None,
                fom_is_chisqr=False, base_sig=None, fix_sig_diff=False, init_sig_factor=1.0,
                plot=False, quiet=False):

    nspec, npix = flux.shape

    # Get zeroth model just for setup
    model, p, a0 = initialize_model(wave, tpl, sampling_factor=sampling_factor,
                                    sig_degree=sig_degree, sig_pwl=sig_pwl, add_degree=add_degree,
                                    mult_degree=mult_degree,
                                    base_mult=None if base_mult is None else base_mult[0,:], 
                                    fom_is_chisqr=fom_is_chisqr,
                                    base_sig=None if base_sig is None else base_sig[0,:], 
                                    fix_sig_diff=fix_sig_diff, init_sig_factor=init_sig_factor)

    # Arrays for the results
    best_fit_model = numpy.empty(flux.shape, dtype=float)
    best_fit_par = numpy.empty((nspec, model.np), dtype=float)
    best_fit_par_err = numpy.zeros((nspec, model.np), dtype=float)
    best_fit_chi = numpy.empty(nspec, dtype=float)
    best_fit_rchi = numpy.empty(nspec, dtype=float)
    sigma_resolution_match = numpy.empty(flux.shape, dtype=float)
    estimated_resolution = numpy.empty(flux.shape, dtype=float)

    velscale = spectrum_velocity_scale(wave)
    tpl_sres = spectral_resolution(tpl['WAVE'].data, tpl['SPECRES'].data[0,:], log10=True)

    if not quiet:
        print('Velocity scale of spectrum to fit: {0}'.format(velscale))
        print('Velocity scale of the template spectrum: {0}'.format(
                                                        spectrum_velocity_scale(tpl['WAVE'].data)))
        print(' ')
        print('#{0:>4}'.format('Nfit'), end='')
        for j in range(model.np):
            if not model.free[j]:
                continue
            print(' {0:>8}'.format(model.pname[j]), end='')
        print(' {0:>8} {1:>8}'.format('FOM', 'FOM/DOF'))

        print(' {0:>4}'.format('init'), end='')
        for j in range(model.np):
            print(' {0:8.1e}'.format(a0[j]), end='')
        print(' {0:>8} {1:>8}'.format('...', '...'))

    # Iterate through all spectra
    for i in range(nspec):
#    for i in range(14,15):

#        if plot:
#            pyplot.plot(tpl['WAVE'].data, base_sig[i,:])
#            pyplot.show()
        model, p, a0 = initialize_model(wave, tpl, sampling_factor=sampling_factor,
                                        sig_degree=sig_degree, sig_pwl=sig_pwl,
                                        add_degree=add_degree, mult_degree=mult_degree,
                                        base_mult=None if base_mult is None else base_mult[i,:], 
                                        fom_is_chisqr=fom_is_chisqr,
                                        base_sig=None if base_sig is None else base_sig[i,:], 
                                        fix_sig_diff=fix_sig_diff, init_sig_factor=init_sig_factor)
    
        # Set up the mask
        gpm = numpy.arange(npix).astype(int)[good_pixels[i,:]]

        # Degrees of freedom for reduced chi-square
        dof = len(gpm) - model.nfree

        # Set the figure-of-merit for the fit, fit, and get the
        # parameter errors from the covariance matrix
        resid = FitChiSquare(flux[i,:], ferr[i,:], model, goodpixels=gpm) if fom_is_chisqr \
                    else FitResiduals(flux[i,:], ferr[i,:], goodpixels=gpm)
        result = optimize.least_squares(resid, a0, bounds=model.bounds(), method='trf')
        x_covar = calculate_parameter_covariance(result)

        # Save the model data
        best_fit_model[i,:] = model.get_model(result.x)
        best_fit_par[i,:] = model.p
        best_fit_par_err[i,model.free] = numpy.ma.sqrt(numpy.diagonal(x_covar)).filled(-1.)
        best_fit_chi[i] = 2*result.cost
        best_fit_rchi[i] = best_fit_chi[i]/dof

        sigma_resolution_match[i,:], estimated_resolution[i,:] \
                = get_resolution_data(best_fit_par[i,:], model, velscale, tpl_sres(wave))

        # Report
        if not quiet:
            print(' {0:4d}'.format(i+1), end='')
            for j in range(model.np):
                if not model.free[j]:
                    continue
                print(' {0:8.1e}'.format(best_fit_par[i,j]), end='')
            print(' {0:8.1e} {1:8.1e}'.format(best_fit_chi[i], best_fit_rchi[i]))

            print(' {0:>4}'.format(' err'), end='')
            for j in range(model.np):
                if not model.free[j]:
                    continue
                print(' {0:8.1e}'.format(best_fit_par_err[i,j]), end='')
            print(' {0:>8} {1:>8}'.format('...', '...'))
        
        if plot:
            _best_fit = numpy.ma.MaskedArray(best_fit_model[i,:],
                                             mask=numpy.invert(good_pixels[i,:]))
            fit_diff = numpy.ma.MaskedArray(flux[i,:] - best_fit_model[i,:],
                                            mask=numpy.invert(good_pixels[i,:]))
            w,h = pyplot.figaspect(1)
            fig = pyplot.figure(figsize=(1.5*w,1.5*h))

            Dy = (numpy.ma.amax(_best_fit)-numpy.ma.amin(_best_fit))*1.1
            ylim = [0,0]
            ylim[0] = (numpy.ma.amax(_best_fit)+numpy.ma.amin(_best_fit) - Dy)/2.
            ylim[1] = ylim[0] + Dy

            Ddy = Dy/3.
            dylim = [0,0]
            dylim[0] = (numpy.ma.amax(fit_diff)+numpy.ma.amin(fit_diff) - Ddy)/2.
            dylim[1] = dylim[0] + Ddy

            ax = fig.add_axes([0.1, 0.3, 0.8, 0.6])
            ax.set_ylim(ylim)
            ax.set_xlim([wave[0], wave[-1]])
            ax.step(wave, flux[i,:], where='mid', color='k', lw=0.5)
            ax.plot(wave, _best_fit, color='C3', lw=0.5)

            ax = fig.add_axes([0.1, 0.1, 0.8, 0.2])
            ax.set_xlim([wave[0], wave[-1]])
            ax.set_ylim(dylim)
            ax.step(wave, fit_diff, where='mid', color='k', lw=0.5)
            ax.step(wave, ferr[i,:], where='mid', color='C0', lw=0.5)
    
            pyplot.show()

    return best_fit_model, best_fit_par, best_fit_par_err, best_fit_chi, best_fit_rchi, \
                sigma_resolution_match, estimated_resolution
    

def get_resolution_data(best_fit_par, model, velscale, tpl_sres):

    sigma_resolution_match = model.get_sigma_function(best_fit_par[model.free])*velscale

    s = constants().sig2fwhm
    c = astropy.constants.c.to('km/s').value
    estimated_resolution = numpy.ma.power(numpy.square(c/s/tpl_sres)
                                          + numpy.square(sigma_resolution_match), -0.5)*c/s

    return sigma_resolution_match, estimated_resolution


def get_polynomial_model(npix, par, degree, s=None):
    if degree is None:
        return None
    nsets = par.shape[0]
    poly = LegendrePolynomial(degree)
    x = numpy.linspace(-1., 1., npix)
    model = numpy.empty((nsets,npix), dtype=float)
    for i in range(nsets):
        model[i,:] = poly(x, a=par[i,:] if s is None else par[i,s:s+poly.np])
    return model

class RukiaBitMask(BitMask):
    """
    Global mask bits used for Rukia.
    """
    def __init__(self):
        super(RukiaBitMask, self).__init__(['LOWSIGMA'])

class Rukia:
    """
    Fit a spectrum by convolving a *single* template spectrum with a
    Gaussian kernel that has a wavelength-dependent width.

    This instantiation defines up to three functional forms used to
    define the fitted model. The objects used to describe each
    functional form must have an ``np`` attribute that provides the
    number of parameters and it must be callable with a form::

        func(x, a=a)

    where ``func`` is the object, ``x`` is the coordinate vector, and
    ``a`` are the model parameters. For example, see
    :class:`LegendrePolynomial`.

    The total number of parameters for the fitted model is the sum of
    the number of parameters for each functional form, plus one
    additional parameter if a Doppler shift is also included in the
    fit (see ``fit_shift``).

    Args:
        sigma_model (object):
            An object that provides the parametrized form of sigma as
            a function of wavelength.
        add_model (object, optional):
            An object that provides the parametrized form of an
            additive function included in the model.
        mul_model (object, optional):
            An object that provides the parametrized form of a
            multiplicative function included in the model.
    """
    bitmask = RukiaBitMask()
    def __init__(self, sigma_model, add_model=None, mul_model=None):
        self.sigma = sigma_model
        self.nsp = self.sigma.np
        self.add = add_model
        self.nap = 0 if self.add is None else self.add.np
        self.mul = mul_model
        self.nmp = 0 if self.mul is None else self.mul.np
        self.np = 1 + self.nsp + self.nap + self.nmp

        # Model parameters
        self.par = numpy.ones(self.np, dtype=float)
        # Initialize all parameters as free
        self.free = numpy.ones(self.np, dtype=bool)

        # Objects with data for the spectrum to fit
        self.wave = None
        self.flux = None
        self.err = None
        self.mask = None
        self.sres = None
        # Objects with data for the template spectrum
        self.tpl_wave = None
        self.tpl_flux = None
        self.tpl_sres = None

        # Bitmask value for modeling results/input
        self.flags = bitmask.minimum_dtype()(0)

    @staticmethod
    def _init_data(wave, flux, err, mask, sres):
        # Set the data to fit
        if flux.dim > 1:
            raise ValueError('Flux vector must be 1D.')
        _mask = numpy.zeros(flux.shape, dtype=bool) if mask is None else mask.copy()
        if _mask.shape != flux.shape:
            raise ValueError('Mask and flux vectors do not have the same shape.')
        if isinstance(flux, numpy.ma.MaskedArray):
            _flux = flux.data.copy()
            _mask |= numpy.ma.getmaskarray(flux)
        else:
            _flux = flux.copy()
        if wave.shape != _flux.shape:
            raise ValueError('Wavelength and flux vectors do not have the same shape.')
        _err = None
        if err is not None:
            if err.shape != _flux.shape:
                raise ValueError('Flux error and flux vectors do not have the same shape.')
            if isinstance(err, numpy.ma.MaskedArray) else err.copy()
                _err = err.data.copy()
                _mask |= numpy.ma.getmaskarray(err)
            else:
                _err = err.copy()
        _sres = None
        if sres is not None:
            if sres.shape != _flux.shape:
                raise ValueError('Spectral resolution and flux vectors do not have the same shape.')
            _sres = interpolate_masked_vector(sres) if isinstance(sres, numpy.ma.MaskedArray) \
                        else sres.copy()
        return _flux, _err, _mask, _sres

    @property
    def has_template(self):
        return self.tpl_wave is not None and self.tpl_flux is not None

    def model(self, par, wave=None, tpl_wave=None, tpl_flux=None, rng=None):
        r"""
        Construct the model spectrum.

        Args:
            par (`numpy.ndarray`_):

                List of model parameters. Shape must be
                :math:`(N_{\rm par},)` or the number of free
                parameters (see :attr:`free`).
            
            wave (`numpy.ndarray`_, optional):
                The wavelength vector on which to construct the
                model. If not provided, the wavelengths for the model
                must have already been defined for this instance,
                e.g., via a fit to a provided spectrum (see
                :func:`fit`).
            tpl_wave(`numpy.ndarray`_, optional):
                Wavelength vector for template spectrum. If provided,
                ``tpl_flux`` must also be provided. If not provided,
                the template spectrum must have already been defined
                for this instance, e.g., via a fit to a provided
                spectrum (see :func:`fit`).
            tpl_flux(`numpy.ndarray`_, optional):
                Flux vector for template spectrum. If provided,
                ``tpl_wave`` must also be provided. If not provided,
                the template spectrum must have already been defined
                for this instance, e.g., via a fit to a provided
                spectrum (see :func:`fit`).
            rng (`numpy.ndarray`_, optional):
                The wavelength range used to scale the functional
                forms used for the model parameters. If None, the
                range is assumed to be already defined or irrelevant
                to the construction of the model.

        Returns:
            `numpy.ndarray`_: The model spectrum.
        """
        if par.dim != 1:
            raise ValueError('Parameter array must be a 1D vector.')
        nfree = numpy.sum(self.free)
        if par.size != self.np and par.size != nfree:
            raise ValueError('Parameter vector must have {0} or {1} elements.'.format(
                                self.np, nfree))
        if (tpl_wave is None or tpl_flux is None) and not self.has_template:
            raise ValueError('No template data available to construct model.')
        if (tpl_wave is None and tpl_flux is not None) \
                or (tpl_wave is not None and tpl_flux is None):
            raise ValueError('If providing tpl_wave or tpl_flux, must provide both.')
        if wave is None and self.wave is None:
            raise ValueError('No wavelength vector available for model.')
        _wave = self.wave if wave is None else wave
        _tpl_wave = self.tpl_wave if tpl_wave is None else tpl_wave
        _dw = angstroms_per_pixel(_tpl_wave, regular=False) if self.dw is None else self.dw
        _tpl_flux = self.tpl_flux if tpl_flux is None else tpl_flux
        self._set_par(par)
        if rng is not None:
            self._reset_range(rng)

        # Get the convolution kernel width in angstroms for at each
        # template wavelength, accounting for the shift in the template
        # wavelengths
        sigma = self.sigma(_tpl_wave/(1+self.par[0]), a=self.par[1:1+self.nsp])
        if numpy.any(sigma <= 0.01):
            # Imposes a minimum sigma for the convolution. This is set
            # to be the same as what is used by
            # ppxf.ppxf_util.gaussian_filter1d
            sigma = numpy.clip(sigma, 0.01, None)
            self.flags = self.bitmask.turn_on(self.flags, 'LOWSIGMA')

        # Convolve the template spectrum, where the sigma is converted
        # to the number of pixels
        # TODO: Should insist that spectra are regularly gridded.
        # Otherwise, the convertion to pixels below is not strictly
        # true if there are significant differences in the width of
        # adjacent pixels, relative to the width of the kernel.
        model_flux = gaussian_filter1d(_tpl_flux, sigma/_dw)

        # Resample the convolved data to the output wavelength grid
        model_flux = Resample(model_flux, x=_tpl_wave/(1+self.par[0]), )

           def __init__(self, y, e=None, mask=None, x=None, xRange=None, xBorders=None, inLog=False,
                 newRange=None, newpix=None, newLog=True, newdx=None, base=10.0, ext_value=0.0,
                 conserve=False, step=True): 




    def _set_par(self, par):
        """
        Set the parameters by accounting for any fixed parameters.
        """
        if par.size == self.np:
            self.par = par
            return
        self.par[self.free] = par
        return self.par

    def _reset_range(self, rng):
        for f in [self.sigma, self.add, self.mul]:
            if f is None:
                continue
            if hasattr(f, 'reset_range'):
                f.reset_range(rng)

    def fit(self, wave, flux, tpl_wave, tpl_flux, err=None, mask=None, sres=None, tpl_sres=None,
            shift=0.0, fit_shift=True, width=None, rejiter=None, rejsig=3.0, rejbox=101):
        r"""
        Fit the spectrum.

        If err is provided:
            - the fit figure-of-merit is chi-square
            - rejections are based on error-normalized residuals
        Otherwise:
            - the fit figure-of-merit is the fit RMS
            - rejections are based on un-normalized residuals.

        Args:
            shift (:obj:`float`, optional):
                The Doppler (multiplicative) shift between the input
                spectrum wavelengths and the template wavelengths;
                i.e., this is :math:`z` is the Doppler shift such
                that the best match in the spectra is achieved when
                the wavelengths of the spectrum to fit are
                
                .. math::
                    
                    \lambda = (1+z) \lambda_{\rm tpl}.

                If ``fit_shift`` is False, this shift is fixed during
                the fit; otherwise, it is included as a fitted model
                parameter and the provided value is used as an
                initial guess.
            fit_shift (:obj:`bool`, optional):
                In addition to the broadening, fit a Doppler (i.e.,
                multiplicative) shift between the template and target
                wavelengths.
        """
        # Setup the spectrum to fit
        self.wave = wave
        self.flux, self.err, self.mask, self.sres \
                = self.init_data(wave, flux, err, mask, sres)
        # Setup the template spectrum
        self.tpl_wave = tpl_wave
        self.tpl_flux, _, _, self.tpl_sres \
                = self.init_data(tpl_wave, tpl_flux, None, None, tpl_sres)
        if self.wave[-1] < self.tpl_wave[0] or self.wave[0] > self.tpl_wave[-1]:
            raise ValueError('No overlapping spectral region between the spectrum to fit and '
                             'the template.')
        # Fit the overlapping spectral range
        self.rng = numpy.array([max(self.wave[0], self.tpl_wave[0]),
                                min(self.wave[-1], self.tpl_wave[-1])])




def fit_tpl_to_twilight_spectrum(
                                 # Input data to fit
                                 wave, flux, ferr, sres,
                                 # Model definition
                                 libkey='BASS10k', sig_degree=None, sig_pwl=False, add_degree=None,
                                 mult_degree=None,
                                 # Fit range, masking, and smoothing
                                 wavelength_range=None, fit_mask_wave=None,
                                 smooth_mask_wave=None, smoothing_boxcar=70,
                                 # Resampling details
                                 dlogLam=1e-4, sampling_factor=10,
                                 # Fitting details
                                 pix_buffer=10, niter=None, fom_is_chisqr=False,
                                 use_sig_diff=False, fix_sig_diff=False, init_sig_factor=1.0,
                                 # Fit plot
                                 plot=False, ofile=None):
    """
    If sig_degree is None and use_sig_diff is None, then only the shift
    is applied, not the broadening.
    """
    # Check input
    _niter = 1 if niter is None or niter < 1 else niter
    if sig_pwl and use_sig_diff:
        raise ValueError('Cannot use piece-wise linear function for sigma and the nominal form '
                         'of the difference in resolution.')
    if fit_mask_wave is not None and fit_mask_wave.shape[-1] != 2:
        raise ValueError('Bad fit_mask_wave')
    _fit_mask_wave = None if fit_mask_wave is None else numpy.atleast_2d(fit_mask_wave)
    if smooth_mask_wave is not None and smooth_mask_wave.shape[-1] != 2:
        raise ValueError('Bad smooth_mask_wave')
    _smooth_mask_wave = None if smooth_mask_wave is None else numpy.atleast_2d(smooth_mask_wave)

    # Check the arrays
    if len(flux.shape) not in [ 1, 2 ]:
        raise ValueError('Data arrays must be 1 or 2 dimensional.')
    if ferr.shape != flux.shape:
        raise ValueError('Error array does not match shape of flux array.')
    if sres.shape != flux.shape:
        raise ValueError('Resolution array does not match shape of flux array.')
    _flux = numpy.ma.atleast_2d(flux.copy() if isinstance(flux, numpy.ma.MaskedArray) \
                                            else numpy.ma.MaskedArray(flux.copy()))
    _ferr = numpy.ma.atleast_2d(ferr.copy() if isinstance(ferr, numpy.ma.MaskedArray) \
                                            else numpy.ma.MaskedArray(ferr.copy()))
    _sres = numpy.atleast_2d(sres.copy())
    nspec = _flux.shape[0]

    # Resample the twilight spectrum
    rwave, rflux, rferr, rsres, rmask = resample_twilight_spectra(wave, _flux, _ferr, _sres,
                                                                  dlogLam=dlogLam,
                                                                  wavelength_range=wavelength_range)
#    for i in range(rflux.shape[0]):
#        pyplot.plot(rwave, rsres[i,:], color='k', lw=0.5)
#    pyplot.show()
#    exit()

    # Set the matching wavelength range for the template (slight
    # adjustments to make sure the number of pixels is correct)
    offset = (1.-1./sampling_factor)*0.5*dlogLam
    doffset = 0.1*dlogLam/sampling_factor
    tpl_wavelength_range = [ numpy.power(10, numpy.log10(rwave[0]) - offset),
                             numpy.power(10, numpy.log10(rwave[-1]) + offset + doffset) ]
    
    # Get the template library spectrum
    tpllibs = template_library_list()
    tpl = TemplateLibrary(libkey, tpllib_list=tpllibs, spectral_step=dlogLam/sampling_factor,
                          log=True, renormalize=False, wavelength_range=tpl_wavelength_range,
                          hardcopy=False)

#    print(wavelength_range)
#    print(tpl_wavelength_range)
#    print(tpl['WAVE'].data[0], tpl['WAVE'].data[-1])
#    print(rwave[0], numpy.power(10, numpy.mean(numpy.log10(tpl['WAVE'].data[0:sampling_factor]))))
#    print(rwave[-1], numpy.power(10, numpy.mean(numpy.log10(tpl['WAVE'].data[-sampling_factor:]))))
#    print(spectrum_velocity_scale(wave))
#    print(spectrum_velocity_scale(rwave))
#    print(spectrum_velocity_scale(tpl['WAVE'].data))
#    print(tpl['WAVE'].data.shape, rwave.shape)
#
#    pyplot.step(tpl['WAVE'].data, tpl['FLUX'].data[0,:], where='mid', color='r', lw=0.5, zorder=1)
#    pyplot.step(rwave, numpy.ma.MaskedArray(rflux[0,:], mask=rmask[0,:]), where='mid', color='k',
#                lw=0.5, zorder=2)
#    pyplot.show()
#    exit()

    # Debugging check
    npix = rwave.size
    if tpl['WAVE'].data.size != npix*sampling_factor:
        raise ValueError('wavelength vector shapes must match to proceed')

    # Get the nominal sigma function needed to match the spectral resolution
    tpl_sres = spectral_resolution(tpl['WAVE'].data, tpl['SPECRES'].data[0,:], log10=True)
    base_sig = numpy.empty((nspec,tpl['WAVE'].data.size), dtype=float)
    match_base_sig = numpy.empty(rflux.shape, dtype=float)
    for i in range(nspec):
        obj_sres = spectral_resolution(rwave, rsres[i,:], log10=True)
        tpl_sres.match(obj_sres)
        base_sig[i,:] = tpl_sres.sig_pd/sampling_factor

        interp = interpolate.interp1d(tpl['WAVE'].data, base_sig[i,:])
        match_base_sig[i,:] = interp(rwave)

    # Get flags for the good pixels in each spectrum
    good_pixels = numpy.invert(rmask)
    good_pixels[:,:pix_buffer] = False
    good_pixels[:,(npix-pix_buffer):] = False
    if _fit_mask_wave is not None:
        indx = numpy.zeros(npix, dtype=bool)
        for mw in _fit_mask_wave:
            indx |= ((rwave > mw[0]) & (rwave < mw[1]))
        good_pixels[:,indx] = False

    # Fit all the spectra
    best_fit_model, best_fit_par, best_fit_par_err, best_fit_chi, best_fit_rchi, \
        sigma_resolution_match, estimated_resolution \
                = fit_spectra(rwave, rflux, rferr, rsres, good_pixels, tpl,
                              sampling_factor=sampling_factor, sig_degree=sig_degree,
                              sig_pwl=sig_pwl, add_degree=add_degree, mult_degree=mult_degree, 
                              fom_is_chisqr=fom_is_chisqr, base_sig=base_sig,
                              fix_sig_diff=fix_sig_diff, init_sig_factor=init_sig_factor,
                              plot=plot)
  
    if _niter <= 1:
        return rwave, rflux, rferr, rsres, rmask, good_pixels, best_fit_model, best_fit_par, \
                best_fit_par_err, best_fit_chi, best_fit_rchi, sigma_resolution_match, \
                estimated_resolution, match_base_sig, None

    bf = BoxcarFilter(smoothing_boxcar)
    smoothing_mask = build_smoothing_mask(rwave, rmask, pix_buffer, mask_wave=_smooth_mask_wave)

    nsig_par = 0 if sig_degree is None and not use_sig_diff \
                    else (sig_degree+1 if sig_degree is not None else 1)
    nadd_par = 0 if add_degree is None else add_degree+1

    for j in range(1,_niter):
        add_cnt = numpy.zeros(rflux.shape, dtype=float) if add_degree is None \
                        else get_polynomial_model(npix, best_fit_par, add_degree, s=1+nsig_par)
        mult_cnt = numpy.ones(rflux.shape, dtype=float) if mult_degree is None \
                        else get_polynomial_model(npix, best_fit_par, mult_degree,
                                                  s=1+nsig_par+nadd_par)

        fit_ratio = numpy.ma.divide(rflux, numpy.ma.divide(best_fit_model-add_cnt,mult_cnt))
        smooth_fit_ratio = bf.smooth(fit_ratio, mask=smoothing_mask, lo=5, hi=5, local_sigma=True)
#        good_pixels[bf.output_mask] = False

        if plot:
            pyplot.step(rwave, fit_ratio[14,:], where='mid', color='k', lw=0.5)
            pyplot.plot(rwave, smooth_fit_ratio[14,:], color='C3', lw=0.5)
            pyplot.show()

        # Fit all the spectra
        best_fit_model, best_fit_par, best_fit_par_err, best_fit_chi, best_fit_rchi, \
                sigma_resolution_match, estimated_resolution \
                    = fit_spectra(rwave, rflux, rferr, rsres, good_pixels, tpl,
                                  sampling_factor=sampling_factor, sig_degree=sig_degree,
                                  sig_pwl=sig_pwl, add_degree=add_degree, mult_degree=mult_degree, 
                                  base_mult=smooth_fit_ratio, fom_is_chisqr=fom_is_chisqr,
                                  base_sig=base_sig, fix_sig_diff=fix_sig_diff,
                                  init_sig_factor=init_sig_factor, plot=plot)

        if j < _niter-1:
            best_fit_model = numpy.ma.divide(best_fit_model, smooth_fit_ratio)

    return rwave, rflux, rferr, rsres, rmask, good_pixels, best_fit_model, best_fit_par, \
                best_fit_par_err, best_fit_chi, best_fit_rchi, sigma_resolution_match, \
                estimated_resolution, match_base_sig, smooth_fit_ratio


def init_ax(fig, win):
    ax = fig.add_axes(win, facecolor='0.95')
    ax.minorticks_on()
    ax.tick_params(which='major', length=10)
    ax.tick_params(which='minor', length=5)
    ax.grid(True, which='major', color='0.8', zorder=2, linestyle='-')
    return ax


def init_ax_twinx(ax):
    axt = ax.twinx()
    axt.minorticks_on()
    axt.tick_params(which='major', length=10)
    axt.tick_params(which='minor', length=5)
    return axt


def get_ylim(y, fac, about_zero=False):
    Dy = (numpy.ma.amax(y) - numpy.ma.amin(y))*fac
    ylim = [ (numpy.ma.amax(y) + numpy.ma.amin(y) - Dy)/2., 0 ]
    ylim[1] = ylim[0] + Dy
    if about_zero:
        max_lim = numpy.amax(numpy.absolute(ylim))
        ylim = [-max_lim, max_lim]
    return ylim


def fit_plot(wave, _flux, good_pixels, _best_fit_model, windows=None, ofile=None):

    if windows is not None and windows.shape != (6,2):
        raise ValueError('Must define three windows.')

    font = { 'size' : 10 }
    rc('font', **font)

    flux = numpy.ma.MaskedArray(_flux, mask=numpy.invert(good_pixels))
    best_fit_model = numpy.ma.MaskedArray(_best_fit_model, mask=numpy.invert(good_pixels))

    xlim = [ wave[good_pixels][0], wave[good_pixels][-1] ]
    t = numpy.ma.append(flux, best_fit_model)
    ylim = get_ylim(t, 1.1)
    resid = flux-best_fit_model
    dylim = get_ylim(resid[good_pixels], 1.1, about_zero=True)

    w,h = pyplot.figaspect(1)
    fig = pyplot.figure(figsize=(1.5*w,1.5*h))

    ax = init_ax(fig, [0.05, 0.51, 0.93, 0.47])
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.xaxis.set_major_formatter(NullFormatter())

    ax.step(wave, flux, where='mid', color='k', lw=0.5, zorder=3)
    ax.plot(wave, best_fit_model, color='r', zorder=4)

#    ax.text(0.2, 0.44, 'MaStar: 8921-3703-57434', horizontalalignment='left',
#            verticalalignment='center', transform=ax.transAxes, fontsize=14)
#    ax.text(0.2, 0.38, fit_lib.library['key'], color='r', horizontalalignment='left',
#            verticalalignment='center', transform=ax.transAxes, fontsize=14)
#    ax.text(0.2, 0.32, 'Resolution fixed relative to DRP function: {0}'.format('False'
#                        if isinstance(model.base_sig, float) else 'True'),
#            horizontalalignment='left', verticalalignment='center', transform=ax.transAxes,
#            fontsize=14)
#    ax.text(0.2, 0.26, r'Resolution difference polynomial order: {0}'.format(model.sd-1),
#            horizontalalignment='left', verticalalignment='center', transform=ax.transAxes,
#            fontsize=14)
#    ax.text(0.2, 0.20, r'Additive polynomial order: {0}'.format(model.ad-1),
#            horizontalalignment='left', verticalalignment='center', transform=ax.transAxes,
#            fontsize=14)
#    ax.text(0.2, 0.14, r'Multiplicative polynomial order: {0}'.format(model.md-1),
#            horizontalalignment='left', verticalalignment='center', transform=ax.transAxes,
#            fontsize=14)
    resid_rms = numpy.sqrt(numpy.mean(numpy.square(resid[good_pixels])))
    ax.text(0.2, 0.08, r'Residual RMS: {0:.4f}'.format(resid_rms), horizontalalignment='left',
            verticalalignment='center', transform=ax.transAxes, fontsize=14)

    if windows is not None:
        for w in windows:
            ax.add_patch(patches.Rectangle((w[0], ylim[0]), numpy.diff(w)[0],
                                           numpy.diff(ylim)[0], color='0.9', zorder=1))

    ax = init_ax(fig, [0.05, 0.43, 0.93, 0.08])
    ax.set_xlim(xlim)
    ax.set_ylim(dylim)

    ax.step(wave, resid, where='mid', color='k', lw=0.5, zorder=3)

    if windows is not None:
        for w in windows:
            ax.add_patch(patches.Rectangle((w[0], dylim[0]), numpy.diff(w)[0], numpy.diff(dylim)[0],
                                           color='0.9', zorder=1))

        ax = [ init_ax(fig, [0.05, 0.23, 0.295, 0.17]),
               init_ax(fig, [0.365, 0.23, 0.295, 0.17]),
               init_ax(fig, [0.68, 0.23, 0.295, 0.17]),
               init_ax(fig, [0.05, 0.03, 0.295, 0.17]),
               init_ax(fig, [0.365, 0.03, 0.295, 0.17]),
               init_ax(fig, [0.68, 0.03, 0.295, 0.17]) ]
        for a in ax:
            a.yaxis.set_major_formatter(NullFormatter())
            a.xaxis.set_major_locator(MultipleLocator(50))
            a.xaxis.set_minor_locator(MultipleLocator(10))

        for i in range(len(ax)):
            ax[i].set_xlim(windows[i,:])
            ax[i].set_ylim([0,1])
            indx = (wave > windows[i,0]) & (wave < windows[i,1])
            if numpy.sum(indx) == 0:
                continue
            t = numpy.ma.append(flux[indx], best_fit_model[indx])
            wylim = get_ylim(t, 1.1)
    
            ax[i].set_ylim(wylim)
    
            ax[i].step(wave[indx], flux[indx], where='mid', color='k', lw=0.5, zorder=3)
            ax[i].plot(wave[indx], best_fit_model[indx], color='r', zorder=4)

    if ofile is None:
        pyplot.show()
    else:
        fig.canvas.print_figure(ofile, bbox_inches='tight')
    fig.clear()
    pyplot.close(fig)
    
    return resid_rms

def colors():
    
    elodie_c = 'DodgerBlue'
    munari_s_c = 'MediumSeaGreen'
    munari_f_c = 'GoldenRod'

    return elodie_c, munari_s_c, munari_f_c

def rms_trend(ofile=None):

    elodie_c, munari_s_c, munari_f_c = colors()

    wavelength_range=[3905, 6800]
    elodie_norm, elodie_rms = fit_tpl_to_mastar_rescale_rms(libkey='ELODIE',
                                                            wavelength_range=wavelength_range,
                                                            add_degree=0, mult_degree=3,
                                                            pix_buffer=100, sampling=[0.5,2.0,60],
                                                            log_sampling=False)#True)
    munari_norm, munari_rms = fit_tpl_to_mastar_rescale_rms(libkey='Munari',
                                                            wavelength_range=wavelength_range,
                                                            add_degree=0, mult_degree=3,
                                                            pix_buffer=100, sampling=[0.5,2.0,60],
                                                            log_sampling=False)#True)
    wavelength_range=[3630, 10310]
    mask_wave = numpy.array([7510,7600])
    munari_full_norm, munari_full_rms = fit_tpl_to_mastar_rescale_rms(libkey='Munari',
                                                            wavelength_range=wavelength_range,
                                                            add_degree=0, mult_degree=6,
                                                            pix_buffer=200, mask_wave=mask_wave,
                                                            sampling=[0.5,2.0,60],
                                                            log_sampling=False)#True)

    font = { 'size' : 16 }
    rc('font', **font)

    w,h = pyplot.figaspect(1)
    fig = pyplot.figure(figsize=(1.5*w,1.5*h))

    ax = init_ax(fig, [0.15, 0.15, 0.7, 0.7])
    ax.set_xlim([0.4, 2.1])
    ax.set_ylim([0.017, 0.038])
    ax.grid(True, which='minor', color='0.8', zorder=0, linestyle=':')
#    ax.set_xscale('log')
#    ax.set_yscale('log')

    ax.plot(elodie_norm, elodie_rms, color=elodie_c, lw=0.5, zorder=3)
    ax.scatter(elodie_norm, elodie_rms, s=100, marker='.', color=elodie_c, lw=0, zorder=4)
    ax.plot(munari_norm, munari_rms, color=munari_s_c, lw=0.5, zorder=3)
    ax.scatter(munari_norm, munari_rms, s=100, marker='.', color=munari_s_c, lw=0, zorder=4)
    ax.plot(munari_full_norm, munari_full_rms, color=munari_f_c, lw=0.5, zorder=3)
    ax.scatter(munari_full_norm, munari_full_rms, s=100, marker='.', color=munari_f_c, lw=0,
               zorder=4)

    ax.text(0.5, 1.03, 'MaStar: 8921-3703-57434', horizontalalignment='center',
            verticalalignment='center', transform=ax.transAxes, fontsize=20)
    ax.text(0.82, 0.77, r'Munari', horizontalalignment='center',
            verticalalignment='center', transform=ax.transAxes, color=munari_s_c, fontsize=20,
            rotation=35)
    ax.text(0.84, 0.68, r'Munari Full', horizontalalignment='center',
            verticalalignment='center', transform=ax.transAxes, color=munari_f_c, fontsize=20,
            rotation=32)
    ax.text(0.82, 0.54, r'ELODIE', horizontalalignment='center',
            verticalalignment='center', transform=ax.transAxes, color=elodie_c, fontsize=20,
            rotation=37)
    ax.text(0.5, -0.08, r'$R_{\rm DRP}/R_{\rm fit}$', horizontalalignment='center',
            verticalalignment='center', transform=ax.transAxes)
    ax.text(-0.13, 0.5, r'Fit Residual RMS', horizontalalignment='center',
            verticalalignment='center', transform=ax.transAxes, rotation='vertical')

    if ofile is None:
        pyplot.show()
    else:
        fig.canvas.print_figure(ofile, bbox_inches='tight')
    fig.clear()
    pyplot.close(fig)


def smooth_twilight_spectra(ifile, ofile, leg_order=11):

    # Open the file
    hdu = fits.open(ifile)
    nspec = hdu['FLUX'].data.shape[0]
    wave = hdu['WAVE'].data
    flux = numpy.ma.MaskedArray(hdu['FLUX'].data, mask=(hdu['MASK'].data > 0)
                                                            | ~(hdu['FLUX'].data > 0))

    # Set the mask for wavelength regions to ignore during the smoothing
    smoothing_wavelength_mask = numpy.array([ [ 3790, 3990],
                                              [ 5790, 5990],
                                              [ 6850, 6950],
                                              [ 7550, 7700] ])
    smoothing_mask = numpy.any( numpy.array([ (wave > w[0]) & (wave < w[1])
                                                for w in smoothing_wavelength_mask ]), axis=0)

    # Mask the flux array
    _flux = flux.copy()
    _flux[ numpy.array([smoothing_mask]*nspec) ] = numpy.ma.masked

    # Get the median of all spectra
    med_flux = numpy.ma.median(_flux, axis=0)
    indx = numpy.ma.getmaskarray(med_flux)
    _med_flux = numpy.ma.median(_flux.data, axis=0)
    _med_flux[numpy.invert(indx)] = numpy.ma.masked
    masked = numpy.ma.getmaskarray(med_flux)

    # Fit a Legendre polynomial, omitting the above specified wavelength
    # regions
    fit = numpy.polynomial.legendre.legval(wave,
                numpy.polynomial.legendre.legfit(wave[~masked], med_flux.compressed(), leg_order))

    plot_fit = False
    if plot_fit:
        w,h = pyplot.figaspect(1)
        fig = pyplot.figure(figsize=(1.5*w,1.5*h))
    
        ax = init_ax(fig, [0.15, 0.2, 0.8, 0.7])
        ax.step(wave, _med_flux, where='mid', color='C0', lw=0.5)
        ax.step(wave, med_flux, where='mid', color='k', lw=0.5)
        ax.plot(wave, fit, color='C3')
    
        ax.text(-0.15, 0.50, r'Flux', horizontalalignment='center', verticalalignment='center',
                transform=ax.transAxes, rotation='vertical', fontsize=12)
        ax.text(0.50, -0.1, r'Wavelength ($\AA$)', horizontalalignment='center',
                verticalalignment='center', transform=ax.transAxes, fontsize=12)
        ax.text(0.7, 0.9, r'Data', horizontalalignment='left', verticalalignment='center',
                transform=ax.transAxes, color='k', fontsize=12)
        ax.text(0.7, 0.85, r'Masked', horizontalalignment='left', verticalalignment='center',
                transform=ax.transAxes, color='C0', fontsize=12)
        ax.text(0.7, 0.8, r'Polynomial Fit', horizontalalignment='left', verticalalignment='center',
                transform=ax.transAxes, color='C3', fontsize=12)
    
    #    pyplot.show()
        fig.canvas.print_figure('legfit.png')
        fig.clear()
        pyplot.close()
        exit()
    
    # Get the smoothed flux array by optimally rescaling the fit to the
    # median spectrum for each fiber spectrum
    smooth_flux = numpy.empty(flux.shape, dtype=float)
    for i in range(nspec):
        indx = ~_flux.mask[i,:]
        smooth_flux[i,:] = optimal_scale(_flux.data[i,:][indx], fit[indx])*fit

    plot_smoothed=False
    if plot_smoothed:
        w,h = pyplot.figaspect(1)
        fig = pyplot.figure(figsize=(1.5*w,1.5*h))

        ax = fig.add_axes([0.12, 0.68, 0.8, 0.3])
        ax.xaxis.set_major_formatter(NullFormatter())
        ax.imshow(numpy.ma.log10(flux), origin='lower', interpolation='nearest', cmap='inferno',
                  aspect='auto')
        ax.text(1.05, 0.5, r'mgCFrame-00177325', horizontalalignment='center',
                verticalalignment='center', transform=ax.transAxes, rotation='vertical',
                fontsize=12)
    
        ax = fig.add_axes([0.12, 0.375, 0.8, 0.3])
        ax.xaxis.set_major_formatter(NullFormatter())
        ax.imshow(numpy.ma.log10(smooth_flux), origin='lower', interpolation='nearest',
                  cmap='inferno', aspect='auto')
        ax.text(-0.12, 0.5, r'Fiber Index', horizontalalignment='center',
                verticalalignment='center', transform=ax.transAxes, rotation='vertical',
                fontsize=12)
        ax.text(1.05, 0.5, r'Legendre Poly', horizontalalignment='center',
                verticalalignment='center', transform=ax.transAxes, rotation='vertical',
                fontsize=12)
    
        ax = fig.add_axes([0.12, 0.07, 0.8, 0.3])
        ax.imshow(numpy.ma.log10(flux/smooth_flux), origin='lower', interpolation='nearest',
                  cmap='inferno', aspect='auto')
        ax.text(0.5, -0.18, r'Wavelength Channel', horizontalalignment='center',
                verticalalignment='center', transform=ax.transAxes, fontsize=12)
        ax.text(1.05, 0.5, r'Normalized Flux', horizontalalignment='center',
                verticalalignment='center', transform=ax.transAxes, rotation='vertical',
                fontsize=12)

#        pyplot.show()
        fig.canvas.print_figure('normalized_spectra.png')
        fig.clear()
        pyplot.close()
        exit()

    # Write the output
    print('writing {0}'.format(ofile))
    fits.HDUList([ hdu['PRIMARY'].copy(),
                   fits.ImageHDU(flux.data/smooth_flux, header=hdu['FLUX'].header.copy()),
                   fits.ImageHDU(hdu['IVAR'].data*numpy.square(smooth_flux),
                                 header=hdu['IVAR'].header.copy()),
                   hdu['MASK'].copy(),
                   hdu['WAVE'].copy(),
                   hdu['DISP'].copy(),
                   hdu['SLITMAP'].copy(),
                   hdu['SKY'].copy(),
                   fits.ImageHDU(smooth_flux, header=hdu['FLUX'].header.copy(), name='SMFLUX')
                  ]).writeto(ofile, clobber=True)


def fix_dispersion(wave, disp):
    nspec = disp.shape[0]
    indx = ~(disp > 0)
    _disp = disp.copy()
    for i in range(nspec):
        if numpy.sum(~indx[i,:]) == 0:
            print('Skipping {0}/{1}'.format(i+1, nspec))
            continue
        interpolator = interpolate.interp1d(wave[~indx[i,:]], disp[i,~indx[i,:]],
                                            fill_value='extrapolate')
        _disp[i,indx[i,:]] = interpolator(wave[indx[i,:]])
    return _disp


def reduce_twilight_spectra(ifile, err_wgt=False, standard_error=True, by_block=False):
    """
    Get mean spectrum within each BLOCKID
    """

    hdu = fits.open(ifile)
    wave = hdu['WAVE'].data.copy()
    flux = numpy.ma.MaskedArray(hdu['FLUX'].data,
                                mask=(hdu['MASK'].data > 0) | ~(hdu['FLUX'].data > 0))
    mask = (hdu['MASK'].data > 0).astype(int)
    ferr = numpy.ma.power(hdu['IVAR'].data, -0.5)
    ferr[numpy.ma.getmaskarray(flux)] = numpy.ma.masked
    ferr.data[numpy.ma.getmaskarray(flux) | numpy.ma.getmaskarray(ferr)] = 1.0
    print(numpy.sum(ferr.mask), numpy.sum(flux.mask))
    disp = numpy.ma.MaskedArray(fix_dispersion(hdu['WAVE'].data, hdu['DISP'].data),
                                mask=numpy.ma.getmaskarray(flux))

    if by_block:
        blocks, reduce_indices = numpy.unique(hdu['SLITMAP'].data['BLOCKID'], return_index=True)
    else:
        reduce_indices = numpy.array([0])

#    mean_flux = numpy.add.reduceat(flux.filled(0.0), reduce_indices, axis=0)
#    test = numpy.add.reduceat(numpy.ones(flux.shape[0]), reduce_indices)
#    print(mean_flux.shape)
#    print(test)
#    exit()

    # Get the mean spectrum and its propagated properties
    if not err_wgt:
        mean_flux = numpy.add.reduceat(flux.filled(0.0), reduce_indices, axis=0)
        nsum_flux = numpy.add.reduceat(numpy.invert(numpy.ma.getmaskarray(flux)), reduce_indices,
                                       axis=0)
        sdev_flux = numpy.add.reduceat(numpy.square(flux.filled(0.0)), reduce_indices, axis=0)
        merr_flux = numpy.add.reduceat(numpy.square(ferr.filled(0.0)), reduce_indices, axis=0)
        mean_disp = numpy.add.reduceat(numpy.square(disp.filled(0.0)), reduce_indices, axis=0)

        dof_correction = numpy.ma.divide(nsum_flux, numpy.minimum(nsum_flux-1,
                                                        numpy.zeros(nsum_flux.shape, dtype=int)))
        mean_flux = numpy.ma.divide(mean_flux, nsum_flux)
        merr_flux = numpy.ma.divide(numpy.ma.sqrt(merr_flux), nsum_flux)
        sdev_flux = numpy.ma.sqrt((numpy.ma.divide(sdev_flux,nsum_flux) - numpy.square(mean_flux))
                                    * dof_correction)
        serr_flux = sdev_flux/numpy.ma.sqrt(nsum_flux)
        mean_disp = numpy.ma.sqrt(mean_disp/nsum_flux)

        mean_sres = numpy.ma.power(constants().sig2fwhm * mean_disp / wave, -1)
        for i in range(mean_sres.shape[0]):
            mean_sres[i,:] = interpolate_masked_vector(mean_sres[i,:])

#        pyplot.step(wave, numpy.ma.log10(sdev_flux), where='mid', color='k')
#        pyplot.step(wave, numpy.ma.log10(serr_flux), where='mid', color='r')
#        pyplot.step(wave, numpy.ma.log10(merr_flux), where='mid', color='b')
#        pyplot.show()

        return wave, mean_flux, (serr_flux if standard_error else merr_flux), mean_sres

    wgt = hdu['IVAR'].data.copy()
    wgt[flux.mask | (hdu['IVAR'].data < 0.)] = 0.

    ew_mean_flux = numpy.add.reduceat((wgt*flux).filled(0.0), reduce_indices, axis=0)
    nsum_flux = numpy.add.reduceat(wgt > 0, reduce_indices, axis=0)
    wsum_flux = numpy.add.reduceat(wgt, reduce_indices, axis=0)
    ew_sdev_flux = numpy.add.reduceat((wgt*numpy.square(flux)).filled(0.0), reduce_indices, axis=0)
    ew_merr_flux = numpy.add.reduceat((wgt*numpy.square(ferr)).filled(0.0), reduce_indices, axis=0)
    ew_mean_disp = numpy.add.reduceat((wgt*numpy.square(disp)).filled(0.0), reduce_indices, axis=0)

    dof_correction = numpy.ma.divide(nsum_flux, numpy.maximum(nsum_flux-1,
                                        numpy.zeros(nsum_flux.shape, dtype=int)))
    ew_mean_flux = numpy.ma.divide(ew_mean_flux, wsum_flux)
    ew_merr_flux = numpy.ma.divide(numpy.ma.sqrt(ew_merr_flux),wsum_flux)
    ew_sdev_flux = numpy.ma.sqrt((numpy.ma.divide(ew_sdev_flux,wsum_flux)
                                    - numpy.square(ew_mean_flux)) * dof_correction)
    ew_serr_flux = numpy.ma.divide(ew_sdev_flux,numpy.ma.sqrt(nsum_flux))
    ew_mean_disp = numpy.ma.sqrt(numpy.ma.divide(ew_mean_disp,wsum_flux))

    ew_mean_sres = numpy.ma.power(constants().sig2fwhm * ew_mean_disp / wave, -1)
    for i in range(ew_mean_sres.shape[0]):
        ew_mean_sres[i,:] = interpolate_masked_vector(ew_mean_sres[i,:])

#    print(ew_serr_flux.size)
#    print(numpy.sum(numpy.invert(ew_serr_flux.data > 0)))
#    print(numpy.sum(numpy.ma.getmaskarray(ew_serr_flux)))
#    print(numpy.sum(numpy.invert(ew_serr_flux.data > 0) & numpy.invert(numpy.ma.getmaskarray(ew_serr_flux))))
#
#    print(ew_merr_flux.size)
#    print(numpy.sum(numpy.invert(ew_merr_flux.data > 0)))
#    print(numpy.sum(numpy.ma.getmaskarray(ew_merr_flux)))
#    print(numpy.sum(numpy.invert(ew_merr_flux.data > 0) & numpy.invert(numpy.ma.getmaskarray(ew_merr_flux))))
#    exit()

#    pyplot.step(wave, numpy.ma.log10(ew_sdev_flux), where='mid', color='k')
#    pyplot.step(wave, numpy.ma.log10(ew_serr_flux), where='mid', color='r')
#    pyplot.step(wave, numpy.ma.log10(ew_merr_flux), where='mid', color='b')
#    pyplot.show()

#    pyplot.step(wave, ew_mean_flux, where='mid', color='k', zorder=3)
#    pyplot.errorbar(wave, ew_mean_flux, yerr=ew_sdev_flux, fmt='none', capsize=0, ecolor='0.6',
#                    lw=0.5, zorder=1)
#    pyplot.errorbar(wave, ew_mean_flux, yerr=ew_serr_flux, fmt='none', capsize=0, ecolor='0.3',
#                    lw=1, zorder=2)
#
#    pyplot.show()

    return wave, ew_mean_flux, (ew_serr_flux if standard_error else ew_merr_flux), ew_mean_sres


def mean_twilight_spectrum(ifile, err_wgt=False, standard_error=True):

    hdu = fits.open(ifile)
    wave = hdu['WAVE'].data.copy()
    flux = numpy.ma.MaskedArray(hdu['FLUX'].data,
                                mask=(hdu['MASK'].data > 0) | ~(hdu['FLUX'].data > 0))
    mask = (hdu['MASK'].data > 0).astype(int)
    ferr = numpy.ma.power(hdu['IVAR'].data, -0.5)
    ferr[numpy.ma.getmaskarray(flux)] = numpy.ma.masked
    ferr.data[numpy.ma.getmaskarray(flux) | numpy.ma.getmaskarray(ferr)] = 1.0
    print(numpy.sum(ferr.mask), numpy.sum(flux.mask))
    disp = numpy.ma.MaskedArray(fix_dispersion(hdu['WAVE'].data, hdu['DISP'].data),
                                mask=numpy.ma.getmaskarray(flux))

    # Get the mean spectrum and its propagated properties
    if not err_wgt:
        mean_flux = numpy.ma.sum(flux, axis=0)
        nsum_flux = numpy.sum(~flux.mask, axis=0)
        sdev_flux = numpy.ma.sum(numpy.square(flux), axis=0)
        merr_flux = numpy.ma.sum(numpy.square(ferr), axis=0)
        mean_disp = numpy.ma.sum(numpy.square(disp), axis=0)

        mean_flux /= nsum_flux
        merr_flux = numpy.ma.sqrt(merr_flux)/nsum_flux
        sdev_flux = numpy.sqrt(sdev_flux/nsum_flux - numpy.square(mean_flux))
        serr_flux = sdev_flux/numpy.ma.sqrt(nsum_flux)
        mean_disp = numpy.ma.sqrt(mean_disp/nsum_flux)

        mean_sres = numpy.ma.power(constants().sig2fwhm * mean_disp / wave, -1)

        pyplot.step(wave, numpy.ma.log10(ew_sdev_flux), where='mid', color='k')
        pyplot.step(wave, numpy.ma.log10(ew_serr_flux), where='mid', color='r')
        pyplot.step(wave, numpy.ma.log10(ew_merr_flux), where='mid', color='b')
        pyplot.show()

        return wave, mean_flux, (serr_flux if standard_error else merr_flux), mean_sres

    wgt = hdu['IVAR'].data.copy()
    wgt[flux.mask | (hdu['IVAR'].data < 0.)] = 0.
    ew_mean_flux = numpy.ma.sum(wgt*flux.data, axis=0)
    nsum_flux = numpy.ma.sum(wgt > 0, axis=0)
    wsum_flux = numpy.ma.sum(wgt, axis=0)
    ew_sdev_flux = numpy.ma.sum(wgt*numpy.square(flux), axis=0)
    ew_merr_flux = numpy.ma.sum(wgt*numpy.square(ferr), axis=0)
    ew_mean_disp = numpy.ma.sum(wgt*numpy.square(disp), axis=0)

    wsum_flux[ ~(wsum_flux > 0) ] = numpy.ma.masked

    ew_mean_flux = numpy.ma.divide(ew_mean_flux, wsum_flux)
    ew_merr_flux = numpy.ma.divide(numpy.ma.sqrt(ew_merr_flux),wsum_flux)
    ew_sdev_flux = numpy.ma.sqrt(numpy.ma.divide(ew_sdev_flux,wsum_flux)-numpy.square(ew_mean_flux))
    ew_serr_flux = numpy.ma.divide(ew_sdev_flux,numpy.ma.sqrt(nsum_flux))
    ew_mean_disp = numpy.ma.sqrt(numpy.ma.divide(ew_mean_disp,wsum_flux))

    ew_mean_sres = numpy.ma.power(constants().sig2fwhm * ew_mean_disp / wave, -1)

#    pyplot.step(wave, numpy.ma.log10(ew_sdev_flux), where='mid', color='k')
#    pyplot.step(wave, numpy.ma.log10(ew_serr_flux), where='mid', color='r')
#    pyplot.step(wave, numpy.ma.log10(ew_merr_flux), where='mid', color='b')
#    pyplot.show()

#    pyplot.step(wave, ew_mean_flux, where='mid', color='k', zorder=3)
#    pyplot.errorbar(wave, ew_mean_flux, yerr=ew_sdev_flux, fmt='none', capsize=0, ecolor='0.6',
#                    lw=0.5, zorder=1)
#    pyplot.errorbar(wave, ew_mean_flux, yerr=ew_serr_flux, fmt='none', capsize=0, ecolor='0.3',
#                    lw=1, zorder=2)
#
#    pyplot.show()

    return wave, ew_mean_flux, (ew_serr_flux if standard_error else ew_merr_flux), ew_mean_sres


def build_smoothing_mask(wave, default_mask, pix_buffer, mask_wave=None):
    mask = default_mask.copy()
    mask[:,:pix_buffer] = True
    mask[:,wave.size-pix_buffer:] = True
    if mask_wave is not None:
        indx = numpy.zeros(wave.size, dtype=bool)
        for mw in mask_wave:
            indx |= ((wave > mw[0]) & (wave < mw[1]))
        mask[:,indx] = True
    return mask


#-----------------------------------------------------------------------------

if __name__ == '__main__':
    t = time.clock()

    create_smoothed_spectra=False
    if create_smoothed_spectra:
        smooth_twilight_spectra('drp/trunk/mgCFrame-00177325-LOG.fits.gz',
                                './mgCFrame-00177325-LOG-trunk-smooth.fits.gz')
        smooth_twilight_spectra('drp/drldec16/mgCFrame-00177325-LOG.fits.gz',
                                './mgCFrame-00177325-LOG-dec16-smooth.fits.gz')
        exit()


    make_resolution_plot = False
    if make_resolution_plot:
        w,h = pyplot.figaspect(1)
        fig = pyplot.figure(figsize=(1.5*w,1.5*h))

        ax_res = init_ax(fig, [0.15, 0.55, 0.8, 0.43])
        ax_res.set_xlim([3600,10400])
        ax_res.set_ylim([1100,2490])
        ax_res.xaxis.set_major_formatter(NullFormatter())

        ax_dres = init_ax(fig, [0.15, 0.12, 0.8, 0.43])
        ax_dres.set_xlim([3600,10400])
        ax_dres.set_ylim([30,240])

        hdu = fits.open('./mgCFrame-00177325-LOG-trunk-smooth.fits.gz')
        disp = fix_dispersion(hdu['WAVE'].data, hdu['DISP'].data)
        sres = numpy.ma.divide(hdu['WAVE'].data[None,:], constants().sig2fwhm * disp)

        ax_res.plot(hdu['WAVE'].data, numpy.ma.mean(sres, axis=0), color='C0')
        ax_dres.plot(hdu['WAVE'].data, numpy.ma.std(sres, axis=0), color='C0')

        hdu = fits.open('./mgCFrame-00177325-LOG-dec16-smooth.fits.gz')
        disp = fix_dispersion(hdu['WAVE'].data, hdu['DISP'].data)
        sres = numpy.ma.divide(hdu['WAVE'].data[None,:], constants().sig2fwhm * disp)

        ax_res.plot(hdu['WAVE'].data, numpy.ma.mean(sres, axis=0), color='C1')
        ax_dres.plot(hdu['WAVE'].data, numpy.ma.std(sres, axis=0), color='C1')

        pyplot.show()
        exit()


    wave_trunk, flux_trunk, ferr_trunk, sres_trunk \
            = reduce_twilight_spectra('./mgCFrame-00177325-LOG-trunk-smooth.fits.gz',
                                       err_wgt=True, standard_error=True, by_block=True)

    wave_dec16, flux_dec16, ferr_dec16, sres_dec16 \
            = reduce_twilight_spectra('./mgCFrame-00177325-LOG-dec16-smooth.fits.gz',
                                      err_wgt=True, standard_error=True, by_block=True)
#    print(numpy.sum(numpy.ma.getmaskarray(sres_dec16)))
#    for i in range(flux_dec16.shape[0]):
#        pyplot.plot(wave_dec16, sres_dec16[i,:], color='k', lw=0.5)
#    pyplot.show()
#    exit()

#    print(flux_dec16.shape)
#    exit()
#
#    wave_dec16, flux_dec16, ferr_dec16, sres_dec16 \
#            = mean_twilight_spectrum('./mgCFrame-00177325-LOG-dec16-smooth.fits.gz',
#                                     err_wgt=True, standard_error=True)

    compare_rdx_plot=False
    if compare_rdx_plot:
        w,h = pyplot.figaspect(1)
        fig = pyplot.figure(figsize=(1.5*w,1.5*h))

        ax = init_ax(fig, [0.15, 0.68, 0.8, 0.3])
        ax.set_xlim([3600,10400])
        ax.set_ylim([0.15, 1.29])
        ax.xaxis.set_major_formatter(NullFormatter())
        ax.step(wave_trunk, flux_trunk, where='mid', color='C0', lw=0.5, zorder=3)
        ax.step(wave_dec16, flux_dec16, where='mid', color='C1', lw=0.5, zorder=3)

        ax = init_ax(fig, [0.15, 0.38, 0.8, 0.3])
        ax.set_xlim([3600,10400])
        ax.set_ylim([0.0000005, 0.5])
        ax.set_yscale('log')
        ax.xaxis.set_major_formatter(NullFormatter())
        ax.step(wave_trunk, ferr_trunk, where='mid', color='C0', lw=0.5, zorder=3)
        ax.step(wave_trunk, ferr_dec16, where='mid', color='C1', lw=0.5, zorder=3)
        ax.step(wave_trunk, numpy.absolute(flux_trunk-flux_dec16), where='mid', color='0.5', lw=0.5,
                zorder=2)

        ax = init_ax(fig, [0.15, 0.08, 0.8, 0.3])
        ax.set_xlim([3600,10400])
        ax.set_ylim([0.0005,  500])
        ax.set_yscale('log')
#        ax.step(wave_trunk, ferr_trunk/flux_trunk, where='mid', color='C0', lw=0.5, zorder=3)
#        ax.step(wave_trunk, ferr_dec16/flux_dec16, where='mid', color='C1', lw=0.5, zorder=3)
        ax.step(wave_trunk, numpy.absolute(flux_trunk-flux_dec16)/ferr_trunk, where='mid',
                color='0.5', lw=0.5, zorder=2)

        pyplot.show()
        exit()

    fit_mask_wave = numpy.array([#[ 5990, 6130],
                                 [ 6850, 6950],
                                 [ 7550, 7700] ])
#    fit_mask_wave = None
    smooth_mask_wave = numpy.array([#[ 5990, 6130],
                                    [ 6860, 6925],
                                    [ 7580, 7680] ])
    wavelength_range = [ 3630, 8950 ] #10000 ]
    dlogLam = 1e-4
    sampling_factor = 10
    libkey='BASS10k'

    smoothing_boxcar = 400

    trunk = True
    if trunk:
        rwave, rflux, rferr, rsres, rmask, gpm, best_fit_model, best_fit_par, \
                best_fit_par_err, best_fit_chi, best_fit_rchi, sigma_resolution_match, \
                estimated_resolution, base_sig, base_mult \
                        = fit_tpl_to_twilight_spectrum(
                                    wave_trunk, flux_trunk, ferr_trunk, sres_trunk,
#                                     wave_trunk, flux_trunk, ferr_trunk, sres_trunk,
                                    libkey=libkey,
#                                     sig_degree=20, sig_pwl=True,
                                    add_degree=0, mult_degree=0,
#                                     sig_degree=1, sig_pwl=False, add_degree=0, mult_degree=2,
                                    wavelength_range=wavelength_range, fit_mask_wave=fit_mask_wave,
                                    smooth_mask_wave=smooth_mask_wave,
                                    smoothing_boxcar=smoothing_boxcar,
                                    dlogLam=dlogLam, sampling_factor=sampling_factor,
                                    pix_buffer=10, niter=2, fom_is_chisqr=True, #)
                                    use_sig_diff=True, fix_sig_diff=False, init_sig_factor=1.0)#,
                                    #plot=True)

        ofile = './mgCFrame-00177325-LOG-trunk-smooth-fit-blocks-f{0}.fits.gz'.format(str(smoothing_boxcar).zfill(3))
        fits.HDUList([ fits.PrimaryHDU(),
                       fits.ImageHDU(data=rwave, name='WAVE'),
                       fits.ImageHDU(data=rflux, name='FLUX'),
                       fits.ImageHDU(data=rferr, name='FERR'),
                       fits.ImageHDU(data=rsres, name='SRES'),
                       fits.ImageHDU(data=numpy.invert(gpm).astype(numpy.uint8), name='MASK'),
                       fits.ImageHDU(data=best_fit_model, name='FIT'),
                       fits.ImageHDU(data=best_fit_par, name='PAR'),
                       fits.ImageHDU(data=best_fit_par_err, name='PAR_ERR'),
                       fits.ImageHDU(data=numpy.array([best_fit_chi,best_fit_rchi]).T, name='CHI'),
                       fits.ImageHDU(data=sigma_resolution_match, name='SIG_MATCH'),
                       fits.ImageHDU(data=estimated_resolution, name='FIT_SRES'),
                       fits.ImageHDU(data=base_sig, name='BASE_SIG_MATCH'),
                       fits.ImageHDU(data=base_mult, name='BASE_MULT_MATCH')
                     ]).writeto(ofile, overwrite=True)

    dec16 = True
    if dec16:
        rwave, rflux, rferr, rsres, rmask, gpm, best_fit_model, best_fit_par, \
                best_fit_par_err, best_fit_chi, best_fit_rchi, sigma_resolution_match, \
                estimated_resolution, base_sig, base_mult \
                        = fit_tpl_to_twilight_spectrum(
                                    wave_dec16, flux_dec16, ferr_dec16, sres_dec16,
#                                     wave_trunk, flux_trunk, ferr_trunk, sres_trunk,
                                    libkey=libkey,
#                                     sig_degree=20, sig_pwl=True,
                                    add_degree=0, mult_degree=0,
#                                     sig_degree=1, sig_pwl=False, add_degree=0, mult_degree=2,
                                    wavelength_range=wavelength_range, fit_mask_wave=fit_mask_wave,
                                    smooth_mask_wave=smooth_mask_wave,
                                    smoothing_boxcar=smoothing_boxcar,
                                    dlogLam=dlogLam, sampling_factor=sampling_factor,
                                    pix_buffer=10, niter=2, fom_is_chisqr=True, #)
                                    use_sig_diff=True, fix_sig_diff=False, init_sig_factor=1.0)#,
                                    #plot=True)

        ofile = './mgCFrame-00177325-LOG-dec16-smooth-fit-blocks-f{0}.fits.gz'.format(str(smoothing_boxcar).zfill(3))
        fits.HDUList([ fits.PrimaryHDU(),
                       fits.ImageHDU(data=rwave, name='WAVE'),
                       fits.ImageHDU(data=rflux, name='FLUX'),
                       fits.ImageHDU(data=rferr, name='FERR'),
                       fits.ImageHDU(data=rsres, name='SRES'),
                       fits.ImageHDU(data=numpy.invert(gpm).astype(numpy.uint8), name='MASK'),
                       fits.ImageHDU(data=best_fit_model, name='FIT'),
                       fits.ImageHDU(data=best_fit_par, name='PAR'),
                       fits.ImageHDU(data=best_fit_par_err, name='PAR_ERR'),
                       fits.ImageHDU(data=numpy.array([best_fit_chi,best_fit_rchi]).T, name='CHI'),
                       fits.ImageHDU(data=sigma_resolution_match, name='SIG_MATCH'),
                       fits.ImageHDU(data=estimated_resolution, name='FIT_SRES'),
                       fits.ImageHDU(data=base_sig, name='BASE_SIG_MATCH'),
                       fits.ImageHDU(data=base_mult, name='BASE_MULT_MATCH')
                     ]).writeto(ofile, overwrite=True)

#    pyplot.plot(rwave, rsres[0,:], color='k')
#    pyplot.plot(rwave, estimated_resolution[0,:], color='r')
#    pyplot.show()
    exit()
#
#    pyplot.step(rwave, rflux[0,:], where='mid', color='k', lw=0.5)
#    pyplot.plot(rwave, best_fit_model[0,:]*data_renorm[0,:], color='r', lw=0.5)
#    pyplot.show()
#    
#    pyplot.step(rwave, rflux[0,:]/data_renorm[0,:], where='mid', color='k', lw=0.5)
#    pyplot.plot(rwave, best_fit_model[0,:], color='r', lw=0.5)
#    pyplot.show()

    spec_smooth_box = 500
    mask = build_smoothing_mask(rwave, rmask[0,:], 10, mask_wave=smooth_mask_wave)
    s = numpy.ma.MaskedArray(rflux[0,:]/data_renorm[0,:], mask=mask)
    ss = boxcar_smooth_vector(s, spec_smooth_box, lo=2., hi=None, niter=10, return_mask=False)
    pyplot.step(rwave, rflux[0,:]/data_renorm[0,:]/ss, where='mid', color='k', lw=0.5)
#    pyplot.plot(rwave, interpolate_masked_vector(smooth_masked_vector(s, spec_smooth_box)),
#                color='r')
    pyplot.show()

    interp_ss = interpolate.interp1d(rwave, ss, fill_value='extrapolate')
    interp_nn = interpolate.interp1d(rwave, data_renorm[0,:], fill_value='extrapolate')

    tpllibs = template_library_list()
    tpl = TemplateLibrary(libkey, tpllib_list=tpllibs, spectral_step=dlogLam/sampling_factor,
                          log=True, directory_path='.',
                          processed_file='{0}_tpllib.fits'.format(libkey), clobber=True,
                          renormalize=False, wavelength_range=wavelength_range)
    tpl_flux = tpl['FLUX'].data[0,:]/interp_ss(tpl['WAVE'].data)
    tpl_flux /= numpy.mean(tpl_flux)

    hdr = fits.Header()
    hdr['CRVAL1'] = (numpy.log10(tpl['WAVE'].data[0]), 'Log of first wavelength')
    hdr['CRPIX1'] = (1, 'Reference pixel of first wavelength')
    hdr['CDELT1'] = (dlogLam, 'Change in log(lambda) per pixel')
    hdr['CD1_1'] = dlogLam

    fits.HDUList([ fits.PrimaryHDU(data=tpl_flux, header=hdr) ]).writeto('tpl_for_xc.fits',
                                                                         clobber=True)


    wave_over, flux_over, ferr_over, sres_over, mask_over \
            = resample_twilight_spectra(wave_trunk, numpy.ma.atleast_2d(flux_trunk.copy()),
                                        numpy.ma.atleast_2d(ferr_trunk.copy()),
                                        numpy.atleast_2d(sres_trunk.copy()),
                                        dlogLam=dlogLam/sampling_factor,
                                        wavelength_range=wavelength_range)
    obj_flux = flux_over[0,:]/interp_nn(wave_over)/interp_ss(wave_over)
    obj_flux /= numpy.mean(obj_flux)

    hdr = fits.Header()
    hdr['CRVAL1'] = (numpy.log10(wave_over[0]), 'Log of first wavelength')
    hdr['CRPIX1'] = (1, 'Reference pixel of first wavelength')
    hdr['CDELT1'] = (dlogLam, 'Change in log(lambda) per pixel')
    hdr['CD1_1'] = dlogLam

    fits.HDUList([ fits.PrimaryHDU(data=obj_flux, header=hdr) ]).writeto('obj_for_xc.fits',
                                                                         clobber=True)

    

    pyplot.plot(tpl['WAVE'].data, tpl_flux, color='k')
    pyplot.plot(wave_over, obj_flux, color='r')
    pyplot.show()
    exit()


    pyplot.step(rwave, s, where='mid', color='k', lw=0.5)
    pyplot.plot(rwave, interpolate_masked_vector(smooth_masked_vector(s, spec_smooth_box)),
                color='r')
    pyplot.show()
    exit()

    rwave_trunk, rflux_trunk, rferr_trunk, rsres_trunk, rmask_trunk \
            = resample_twilight_spectra(wave_trunk, flux_trunk.reshape(1,-1),
                                        ferr_trunk.reshape(1,-1), sres_trunk.reshape(1,-1))

    _rflux_trunk = numpy.ma.MaskedArray(rflux_trunk, mask=rmask_trunk.astype(bool))

    pyplot.step(wave_trunk, flux_trunk, where='mid', color='k')
    pyplot.step(rwave_trunk, _rflux_trunk[0,:], where='mid', color='r', lw=0.5)
    pyplot.show()
    exit()
    

    wave_dec16, flux_dec16, ferr_dec16, sres_dec16 \
            = mean_twilight_spectrum('./mgCFrame-00177325-LOG-dec16-smooth.fits.gz',
                                     err_wgt=False, standard_error=True)

    pyplot.plot(wave_trunk, flux_trunk, color='r')
    pyplot.plot(wave_dec16, flux_dec16, color='b')
    pyplot.show()

    exit()

    wave_over, flux_over, ferr_over, sres_over, mask_over \
                = resample_twilight_spectra(hdu['WAVE'].data, hdu['FLUX'].data, ferr, sres, mask)

    av_flux = numpy.ma.mean(flux_over, axis=0)
    ew_av_flux = numpy.ma.average(flux_over, axis=0, weights=numpy.ma.power(ferr_over, -2))

    ew_av_flux = numpy.ma.average(flux_over, axis=0, weights=numpy.ma.power(ferr_over, -2))

    pyplot.step(wave_over, av_flux, where='mid', color='k', lw=0.5)
    pyplot.step(wave_over, ew_av_flux, where='mid', color='b', lw=0.5)
    pyplot.show()
    exit()


    w,h = pyplot.figaspect(1)
    fig = pyplot.figure(figsize=(1.5*w,1.5*h))

    ax = fig.add_axes([0.1, 0.68, 0.8, 0.3])
    ax.xaxis.set_major_formatter(NullFormatter())
    ax.imshow(numpy.ma.log10(flux), origin='lower', interpolation='nearest', cmap='inferno', aspect='auto')

    ax = fig.add_axes([0.1, 0.37, 0.8, 0.3])
    ax.xaxis.set_major_formatter(NullFormatter())
    ax.imshow(numpy.log10(ferr), origin='lower', interpolation='nearest', cmap='inferno', aspect='auto')

    ax = fig.add_axes([0.1, 0.06, 0.8, 0.3])
    ax.imshow(sres, origin='lower', interpolation='nearest', cmap='inferno', aspect='auto')

    pyplot.show()

    exit()

    wave, flux, ferr, sres, mask = resample_twilight_spectra(hdu['WAVE'].data, hdu['FLUX'].data,
                                                             ferr, sres, mask)


    av_flux = numpy.ma.mean(flux/smooth_flux, axis=0)
    fitting_wavelength_mask = numpy.array([ [ 5790, 5990],
                                            [ 6850, 6950],
                                            [ 7550, 7700] ])
    fitting_mask = numpy.any( numpy.array([ (wave > w[0]) & (wave < w[1])
                                                for w in fitting_wavelength_mask ]), axis=0)
    fitting_mask[:] = False

    pyplot.step(wave, flux[0,:], where='mid', color='k', lw=0.5)
    pyplot.plot(wave, smooth_flux[0,:], color='r')
    pyplot.show()
    exit()


    disp = numpy.ma.MaskedArray(hdu['DISP'].data, mask=~(hdu['DISP'].data > 0))

    max_disp = numpy.ma.amax(disp, axis=0)[1:] / numpy.diff(wave)

    boxcar = 200
    hbox = boxcar//2
    i = numpy.arange(hbox,wave.size-hbox)
    dwave = numpy.array([ wave[i+hbox]-wave[i-hbox] for i in numpy.arange(hbox,wave.size-hbox) ])
    wave = numpy.array([ numpy.power(10., (numpy.log10(wave[i+hbox])+numpy.log10(wave[i-hbox]))/2)
                                        for i in numpy.arange(hbox,wave.size-hbox) ])
    pyplot.plot(wave, dwave)
    pyplot.show()
    exit()
    print(dwave.shape)
    exit()

    pyplot.plot(wave[1:], max_disp)
    pyplot.show()
    pyplot.plot(wave[1:], max_disp)
    exit()

    print(disp.shape)
    print(max_disp.shape)

    exit()


    print('Elapsed time: {0} seconds'.format(time.clock() - t))



