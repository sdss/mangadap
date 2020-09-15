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
from scipy import optimize

import astropy.constants

from ppxf.ppxf_util import gaussian_filter1d
from ppxf.capfit import cov_err

from ..util.filter import interpolate_masked_vector
from ..util.sampling import angstroms_per_pixel, Resample
from ..util.bitmask import BitMask
from ..util import modeling


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
        self.interp = interpolate.interp1d(modeling.scaled_coordinates(x, rng=self.rng), y,
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

        _wave_scaled = self.wave_scaled if wave is None \
                            else modeling.scaled_coordinates(wave, rng=self.rng)
        _wave = self.wave if wave is None else wave
        _tpl_wave = self.shifted_wave(self.tpl_wave if tpl_wave is None else tpl_wave)
        _tpl_wave_scaled = modeling.scaled_coordinates(_tpl_wave, rng=self.rng)
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
        model_flux = Resample(model_flux, x=_tpl_wave[p:-p], newx=_wave).outy

        # Multiply by any polynomial
        if self.mul is not None:
            model_flux *= self.mul(_wave_scaled, a=self.mul_par())
        # Include an additive polynomial
        if self.add is not None:
            model_flux += self.add(_wave_scaled, a=self.add_par())
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
            self.wave_scaled = modeling.scaled_coordinates(self.wave, rng=self.rng)

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

    def _init_par(self, shift=0.0, fit_shift=True, sigma_p0=None, sigma_f=None, add_p0=None,
                  add_f=None, mul_p0=None, mul_f=None):
        # Check the input parameters
        if sigma_p0 is not None and len(sigma_p0) != self.nsp:
            raise ValueError('Incorrect number of parameters for the sigma model.')
        if sigma_f is not None and len(sigma_f) != self.nsp:
            raise ValueError('Incorrect number of parameter flags for the sigma model.')
        if sigma_p0 is None and sigma_f is not None:
            warnigns.warn('Did not provide sigma model parameters, but did provide sigma model '
                          'parameter fitting flags.  Parameters may get fixed to default guesses.')
        if add_p0 is not None and len(add_p0) != self.nap:
            raise ValueError('Incorrect number of parameters for the additive polynomial.')
        if add_f is not None and len(add_f) != self.nap:
            raise ValueError('Incorrect number of parameter flags for the additive polynomial.')
        if add_p0 is None and add_f is not None:
            warnigns.warn('Did not provide additive polynomial parameters, but did provide '
                          'additive polynomial parameter fitting flags.  Parameters may get '
                          'fixed to default guesses.')
        if mul_p0 is not None and len(mul_p0) != self.nmp:
            raise ValueError('Incorrect number of parameters for the multiplicative polynomial.')
        if mul_f is not None and len(mul_f) != self.nmp:
            raise ValueError('Incorrect number of parameter flags for the multiplicative '
                             'polynomial.')
        if mul_p0 is None and mul_f is not None:
            warnigns.warn('Did not provide multiplicative polynomial parameters, but did provide '
                          'multiplicative polynomial parameter fitting flags.  Parameters may get '
                          'fixed to default guesses.')

        # Initialize the parameter vector
        self.par[0] = shift
        self.free[0] = fit_shift

        # Sigma parameters:
        if self.nsp > 0:
            ps = self._sigma_par_slice()
            self.par[ps] = numpy.array([1.] + [0.]*(self.nsp-1)) if sigma_p0 is None else sigma_p0
            self.free[ps] = True if sigma_f is None else sigma_f

        # Multiplicative polynomial parameters:
        if self.nmp > 0:
            ps = self._mul_par_slice()
            self.par[ps] = numpy.array([1.] + [0.]*(self.nmp-1)) if mul_p0 is None else mul_p0
            self.free[ps] = True if mul_f is None else mul_f

        # Additive polynomial parameters:
        if self.nap > 0:
            ps = self._add_par_slice()
            self.par[ps] = numpy.array([1.] + [0.]*(self.nap-1)) if add_p0 is None else add_p0
            self.free[ps] = True if add_f is None else add_f

        self.nfree = numpy.sum(self.free)

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

    def fit(self, wave, flux, tpl_wave, tpl_flux, err=None, mask=None, sres=None, tpl_sres=None,
            shift=0.0, fit_shift=True, shift_range=[-300.,300.], rejiter=None, rejsig=3.0,
            rejbox=101, sigma_p0=None, sigma_f=None, add_p0=None, add_f=None, mul_p0=None,
            mul_f=None, rng=None):
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
        self._init_par(shift=shift, fit_shift=fit_shift, sigma_p0=sigma_p0, sigma_f=sigma_f,
                       add_p0=add_p0, add_f=add_f, mul_p0=mul_p0, mul_f=mul_f)

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
        resid = modeling.FitResiduals(self.flux, self.model, gpm=self.gpm) if self.err is None \
                    else modeling.FitChiSquare(self.flux, self.err, self.model, gpm=self.gpm)

        # TODO: Change this for the sigma parameters too?
        diff_step = numpy.full(self.np, numpy.sqrt(numpy.finfo(float).eps), dtype=float)
        diff_step[0] = 0.1
        result = optimize.least_squares(resid, self.par[self.free], method='trf',
                                        diff_step=diff_step[self.free],
                                        bounds=(lb[self.free], ub[self.free])) #, verbose=2)

        if rejiter is None or rejiter == 0:
            self._set_par(result.x)
            return

        # Keep the initial mask
        init_gpm = self.gpm.copy()

        # Iteratively reject
        i = 0
        while i < rejiter or rejiter < 0:
            # Reject pixels
            rej = modeling.reject_residuals_1d(resid(self.par), lo=rejsig, hi=rejsig,
                                               boxcar=rejbox)
            if not numpy.any(rej):
                # None were rejected so we're done
                break
            # Mask the rejected values
            self.gpm[numpy.where(self.gpm)[0][rej]] = False
            # Refit, starting at the previous fit position
            result = optimize.least_squares(resid, self.par[self.free], method='trf',
                                            diff_step=diff_step[self.free],
                                            bounds=(lb[self.free], ub[self.free])) #, verbose=2)
            # Reset the best-fit parameters
            self._set_par(result.x)
            # Increment iteration
            i += 1




