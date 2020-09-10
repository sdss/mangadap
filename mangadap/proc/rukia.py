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
import warnings

from IPython import embed

import numpy
from scipy import interpolate, optimize, ndimage
from matplotlib import pyplot, rc, patches
from matplotlib.ticker import NullFormatter, MultipleLocator, FormatStrFormatter

from astropy.io import fits
import astropy.constants

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

def scaled_coordinates(x, rng):
    return x if rng is None else 2 * (x - rng[0]) / (rng[1]-rng[0]) - 1 


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
    def __init__(self, deg, c=None, rng=None):
        self.deg = deg
        self.np = self.deg+1
        self.coeff = numpy.ones(self.np, dtype=float) if c is None else c
        self.rng = rng

    def reset_range(self, rng):
        self.rng = rng

    def __call__(self, x, a=None):
        if a is not None:
            if len(a) != self.np:
                raise ValueError('Incorrect number of parameters')
            self.coeff = a
        return numpy.polynomial.legendre.legval(x, self.coeff)


# TODO:  assess how rng should be used here
class ScaledEmpirical:
    def __init__(self, x, y, scale_func=None, rng=None, assume_sorted=True):
        self.deg = 0 if scale_func is None else scale_func.deg
        self.np = 1 if scale_func is None else scale_func.np
        self.sf = LegendrePolynomial(0) if scale_func is None else scale_func
        self.coeff = numpy.array([1.0]) if scale_func is None else scale_func.coeff
        self.reset_range(rng)
        self.interp = interpolate.interp1d(scaled_coordinates(x, self.rng), y,
                                           fill_value='extrapolate', assume_sorted=assume_sorted)

    def reset_rng(self, rng):
        if hasattr(self.sf, 'reset_range'):
            self.sf.reset_range(rng)
        self.rng = rng

    def __call__(self, x, a=None):
        if a is not None:
            if len(a) != self.np:
                raise ValueError('Incorrect number of parameters')
            self.coeff = a
        return self.sf(x, a=self.coeff)*self.interp(x)
        

class FitResiduals:
    def __init__(self, obj_flux, model, gpm=None):
        self.obj_flux = obj_flux.copy()
        self.model = model
        self.gpm = numpy.ones(obj_flux.size, dtype=bool) if gpm is None else gpm
    def __call__(self, a):
        resid = (self.obj_flux[self.gpm]-self.model(a)[self.gpm])
        print(a) #, numpy.sqrt(numpy.sum(numpy.square(resid))))
        return resid


class FitChiSquare:
    def __init__(self, obj_flux, obj_ferr, model, gpm=None):
        self.obj_flux = obj_flux.copy()
        self.obj_ferr = obj_ferr.copy()
        self.model = model
        self.gpm = numpy.ones(obj_flux.size, dtype=bool) if gpm is None else gpm
    def __call__(self, a):
        return ((self.obj_flux[self.gpm]-self.model(a)[self.gpm])/self.obj_ferr[self.gpm])


class RukiaBitMask(BitMask):
    """
    Global mask bits used for Rukia.
    """
    def __init__(self):
        super(RukiaBitMask, self).__init__(['LOSIGMA', 'HISIGMA'])


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
    """
    Bitmask for Rukia results.
    """

    def __init__(self, sigma_model, mul_model=None, add_model=None):

        # Saved for convenience
        self._c = astropy.constants.c.to('km/s').value

        self.sigma = sigma_model
        self.nsp = self.sigma.np
        self.mul = mul_model
        self.nmp = 0 if self.mul is None else self.mul.np
        self.add = add_model
        self.nap = 0 if self.add is None else self.add.np
        self.np = 1 + self.nsp + self.nmp + self.nap

        # Model parameters
        self.par = numpy.ones(self.np, dtype=float)
        # Initialize all parameters as free
        self.free = numpy.ones(self.np, dtype=bool)
        self.nfree = numpy.sum(self.free)

        # Objects with data for the spectrum to fit
        self.wave_scaled = None
        self.wave = None
        self.dw = None
        self.flux = None
        self.err = None
        self.gpm = None
        self.sres = None
        # Objects with data for the template spectrum
        self.tpl_wave = None
        self.tpl_dw = None
        self.tpl_flux = None
        self.tpl_sres = None

        # Bitmask value for modeling results/input
        self.flags = self.bitmask.minimum_dtype()(0)

        # Range is initialized to None and then changed by the fit()
        # method
        # TODO: Issue warning if any of the models have rng defined?
        # This will overwrite them.
        self._reset_range(None)

    @staticmethod
    def _init_data(wave, flux, err, mask, sres):
        # Set the data to fit
        if flux.ndim > 1:
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
        _dw = angstroms_per_pixel(wave, regular=False)
        _err = None
        if err is not None:
            if err.shape != _flux.shape:
                raise ValueError('Flux error and flux vectors do not have the same shape.')
            if isinstance(err, numpy.ma.MaskedArray):
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
        return _dw, _flux, _err, numpy.logical_not(_mask), _sres

    @property
    def has_template(self):
        return self.tpl_wave is not None and self.tpl_flux is not None

    def shifted_wave(self, wave, shift=None):
        if shift is None:
            shift = self.par[0]
        return wave/(1+shift/self._c)

    def model(self, par, wave=None, tpl_wave=None, tpl_flux=None, rng=None,
              pixel_sigma_bounds=[0.01,100]):
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
        if (tpl_wave is None or tpl_flux is None) and not self.has_template:
            raise ValueError('No template data available to construct model.')
        if (tpl_wave is None and tpl_flux is not None) \
                or (tpl_wave is not None and tpl_flux is None):
            raise ValueError('If providing tpl_wave or tpl_flux, must provide both.')
        if wave is None and self.wave is None:
            raise ValueError('No wavelength vector available for model.')
        if rng is not None:
            self._reset_range(rng)

        self._set_par(par)

        _wave_scaled = self.wave_scaled if wave is None else scaled_coordinates(wave, self.rng)
        _wave = self.wave if wave is None else wave
        _tpl_wave = self.shifted_wave(self.tpl_wave if tpl_wave is None else tpl_wave)
        _tpl_wave_scaled = scaled_coordinates(_tpl_wave, self.rng)
        _tpl_flux = self.tpl_flux if tpl_flux is None else tpl_flux

        # Sigma in angstroms
        sigma = self.sigma(_tpl_wave_scaled, a=self.sigma_par())
        # Convert it to pixels
        dw = angstroms_per_pixel(_tpl_wave, regular=False)
        sigma /= dw
        # Imposes bounds on the allowed sigma in pixels for the
        # convolution. The minimum is, by default, set to be the same
        # as what is used by ppxf.ppxf_util.gaussian_filter1d
        self.flags = 0
        if numpy.any(sigma <= pixel_sigma_bounds[0]):
            self.flags = self.bitmask.turn_on(self.flags, 'LOSIGMA')
        if numpy.any(sigma >= pixel_sigma_bounds[1]):
            self.flags = self.bitmask.turn_on(self.flags, 'HISIGMA')
        # Clip the sigma to the specified bounds
        sigma = numpy.clip(sigma, *pixel_sigma_bounds)

        # TODO: Should insist that spectra are regularly gridded.
        # Otherwise, the convertion to pixels below is not strictly
        # true if there are significant differences in the width of
        # adjacent pixels, relative to the width of the kernel.

        # Convolve the template spectrum. The
        # gaussian_filter1d algorithm results in the first and last
        # ``p`` pixels being set to 0.
        p = int(numpy.ceil(numpy.max(3*sigma)))
        model_flux = gaussian_filter1d(_tpl_flux, sigma)[p:-p]

        # Resample the convolved data to the output wavelength grid
        model_flux = Resample(model_flux, x=_tpl_wave[p:-p], newx=_wave).outy #, step=False).outy

#        model_flux = interpolate.interp1d(_tpl_wave[p:-p], model_flux)(_wave)

        if numpy.any(numpy.isclose(model_flux, 0.)):
            embed()
            exit()

        # Multiply by any polynomial
        if self.mul is not None:
            model_flux *= self.mul(_wave_scaled, a=self.mul_par())
        # Include an additive polynomial
        if self.add is not None:
            model_flux += self.mul(_wave_scaled, a=self.add_par())
        return model_flux

    def _set_par(self, par):
        """
        Set the parameters by accounting for any fixed parameters.
        """
        if par.ndim != 1:
            raise ValueError('Parameter array must be a 1D vector.')
        if par.size == self.np:
            self.par = par.copy()
            return
        if par.size != self.nfree:
            raise ValueError('Must provide {0} or {1} parameters.'.format(self.np, self.nfree))
        self.par[self.free] = par.copy()

    def _reset_range(self, rng):
        self.rng = rng
        for f in [self.sigma, self.add, self.mul]:
            if f is None:
                continue
            if hasattr(f, 'reset_range'):
                f.reset_range(rng)
        if self.wave is not None:
            self.wave_scaled = scaled_coordinates(self.wave, self.rng)

    def _sigma_par_slice(self):
        return slice(1,1+self.nsp)

    def sigma_par(self, par=None):
        if par is None:
            par = self.par
        return par[self._sigma_par_slice()]

    def _mul_par_slice(self):
        s = 1+self.nsp
        return slice(s,s+self.nmp)

    def mul_par(self, par=None):
        if par is None:
            par = self.par
        return par[self._mul_par_slice()]

    def _add_par_slice(self):
        s = 1+self.nsp+self.nmp
        return slice(s,s+self.nap)

    def add_par(self, par=None):
        if par is None:
            par = self.par
        return par[self._add_par_slice()]

    def _init_par(self, shift=0.0, fit_shift=True, sigma_p0=None, add_p0=None, mul_p0=None):
        # Check the input parameters
        if sigma_p0 is not None and len(sigma_p0) != self.nsp:
            raise ValueError('Incorrect number of parameters for the sigma model.')
        if add_p0 is not None and len(add_p0) != self.nap:
            raise ValueError('Incorrect number of parameters for the sigma model.')
        if mul_p0 is not None and len(mul_p0) != self.nmp:
            raise ValueError('Incorrect number of parameters for the sigma model.')
        # Initialize the parameter vector
        self.par[0] = shift
        self.free[0] = fit_shift
        self.nfree = numpy.sum(self.free)
        if sigma_p0 is not None:
            self.par[self._sigma_par_slice()] = sigma_p0
        elif self.nsp > 0:
            self.par[self._sigma_par_slice()] = numpy.array([1.] + [0.]*(self.nsp-1))
        if mul_p0 is not None:
            self.par[self._mul_par_slice()] = mul_p0
        elif self.nmp > 0:
            self.par[self._mul_par_slice()] = numpy.array([1.] + [0.]*(self.nmp-1))
        if add_p0 is not None:
            self.par[self._add_par_slice()] = add_p0
        elif self.nap > 0:
            self.par[self._add_par_slice()] = numpy.array([1.] + [0.]*(self.nap-1))

    def par_scale(self):
        scale = numpy.ones(self.np, dtype=float)
        scale[0] = 1000.
        return scale

    def fitting_range(self, wave, tpl_wave, shift=None, shift_range=None, rng=None):

        # Determine the spectral range to fit
        if wave[-1] < tpl_wave[0] or wave[0] > tpl_wave[-1]:
            raise ValueError('No overlapping spectral region between the spectrum to fit and '
                             'the template.')

        # 300 here is matched to the maximum sigma allowed of 100 in
        # gaussian_filter1d and how it maskes the edges of the
        # spectrum.
        tpl_wave_lim = [self.shifted_wave(tpl_wave[300], shift=shift if shift_range is None
                                                               else shift+shift_range[0]),
                        self.shifted_wave(tpl_wave[-300], shift=shift if shift_range is None
                                                                else shift+shift_range[1])]

        # Get the overlapping spectral range to fit
        # TODO: Alter this based on the range of allowed shifts
        _rng = numpy.array([max(wave[0], tpl_wave_lim[0]), min(wave[-1], tpl_wave_lim[-1])])

        if rng is None:
            return _rng

        if rng[0] < _rng[0] or rng[0] > _rng[1] or rng[1] > _rng[1] or rng[1] < _rng[0]:
            raise ValueError('Provided range is invalid.  Maximum range is: '
                             '{0}, but requested range is {1}'.format(_rng, rng))
        return rng

#    def _reject(self, rejiter=None, rejsig=3.0, rejbox=101):


    @staticmethod
    def reject_model_outliers(obj_flux, ppxf_result, rescale=False, local_sigma=False, boxcar=None,
                              nsigma=3.0, niter=9, loggers=None, quiet=False):
        if boxcar is None and local_sigma:
            raise ValueError('For local sigma determination, must provide boxcar size.')
        model_flux = PPXFFit.compile_model_flux(obj_flux, ppxf_result, rescale=rescale)

        if local_sigma:
            if not quiet:
                log_output(loggers, 1, logging.INFO,
                           'Rejecting using local sigma (boxcar is {0} pixels).'.format(boxcar))
            residual = obj_flux-model_flux # This should be masked where the data were not fit
            bf = BoxcarFilter(boxcar, lo=nsigma, hi=nsigma, niter=niter, y=residual,
                              local_sigma=True)
            obj_flux[bf.output_mask] = numpy.ma.masked
            return obj_flux

        if not quiet:
            log_output(loggers, 1, logging.INFO, 'Rejecting full spectrum outliers.')
        for i in range(niter):
            residual = obj_flux-model_flux  # This should be masked where the data were not fit
            sigma = numpy.array([numpy.ma.std(residual, axis=1)]*residual.shape[1]).T
            indx = numpy.absolute(residual) > nsigma*sigma
            old_mask = numpy.ma.getmaskarray(obj_flux).copy()
            obj_flux[indx] = numpy.ma.masked
            if numpy.sum(numpy.ma.getmaskarray(obj_flux)) == numpy.sum(old_mask):
                break
        return obj_flux

    def fit(self, wave, flux, tpl_wave, tpl_flux, err=None, mask=None, sres=None, tpl_sres=None,
            shift=0.0, fit_shift=True, shift_range=[-300.,300.], rejiter=None, rejsig=3.0,
            rejbox=101, sigma_p0=None, add_p0=None, mul_p0=None, rng=None):
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
        self._init_par(shift=shift, fit_shift=fit_shift, sigma_p0=sigma_p0, add_p0=add_p0,
                       mul_p0=mul_p0)

        self.wave = wave
        self.tpl_wave = tpl_wave
        self._reset_range(self.fitting_range(wave, tpl_wave, shift, shift_range=shift_range,
                          rng=rng))

        indx = (self.wave > self.rng[0]) & (self.wave < self.rng[1])
        if numpy.sum(indx) != len(indx):
            warnings.warn('Limiting fit to {0} - {1}.'.format(*self.rng))

        # Setup the spectrum to fit
        self.wave_scaled = self.wave_scaled[indx]
        self.wave = self.wave[indx]
        self.dw, self.flux, self.err, self.gpm, self.sres \
                = self._init_data(wave[indx], flux[indx], None if err is None else err[indx],
                                  None if mask is None else mask[indx],
                                  None if sres is None else sres[indx])

        # Setup the template spectrum
        self.tpl_dw, self.tpl_flux, _, _, self.tpl_sres \
                = self._init_data(tpl_wave, tpl_flux, None, None, tpl_sres)

        lb = numpy.full(self.np, -numpy.inf, dtype=float)
        ub = numpy.full(self.np, numpy.inf, dtype=float)
        lb[0], ub[0] = shift + numpy.array(shift_range)

        # Set the figure-of-merit for the fit, perform the fit, and get
        # the parameter errors from the covariance matrix
        resid = FitResiduals(self.flux, self.model, gpm=self.gpm) if self.err is None \
                    else FitChiSquare(self.flux, self.ferr, self.model, gpm=self.gpm)

        diff_step = numpy.full(self.np, numpy.sqrt(numpy.finfo(float).eps), dtype=float)
        diff_step[0] = 0.1
        result = optimize.least_squares(resid, self.par[self.free], method='trf',
                                        x_scale=self.par_scale()[self.free],
                                        diff_step=diff_step[self.free],
                                        bounds=(lb[self.free], ub[self.free]), verbose=2)
        self._set_par(result.x)
        result = optimize.least_squares(resid, self.par[self.free], method='trf',
                                        x_scale=self.par_scale()[self.free],
#                                        diff_step=diff_step[self.free],
                                        bounds=(lb[self.free], ub[self.free]), verbose=2)

        print(result.x)
        embed()
        exit()


def calculate_parameter_covariance(result):
    a = numpy.dot(result.jac.T,result.jac)
    try:
        return numpy.linalg.inv(a)
    except:
        return numpy.linalg.pinv(a)



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



