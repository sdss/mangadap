# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Implements an emission-line profile fitting class.

*License*:
    Copyright (c) 2015, SDSS-IV/MaNGA Pipeline Group
        Licensed under BSD 3-clause license - see LICENSE.rst

*Source location*:
    $MANGADAP_DIR/python/mangadap/proc/elric.py

*Imports and python version compliance*:
    ::

*Class usage examples*:
        Add examples

*Revision history*:
    | **26 Apr 2016**: Original implementation by K. Westfall (KBW)
    | **13 Jul 2016**: (KBW) Include log_bounds determining whether or
        not a returned parameters is near its boundary.
    | **19 Jul 2016**: (KBW) Changed file name

.. _glob.glob: https://docs.python.org/3.4/library/glob.html
.. _configparser.ConfigParser: https://docs.python.org/3/library/configparser.html#configparser.ConfigParser
.. _scipy.optimize.least_squares: http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html
.. _scipy.optimize.OptimizeResult: http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.OptimizeResult.html

.. _astropy.io.fits.hdu.hdulist.HDUList: http://docs.astropy.org/en/v1.0.2/io/fits/api/hdulists.html
.. _astropy.modeling: http://docs.astropy.org/en/stable/modeling/index.html
.. _astropy.modeling.FittableModel: http://docs.astropy.org/en/stable/api/astropy.modeling.FittableModel.html
.. _astropy.modeling.polynomial.Legendre1D: http://docs.astropy.org/en/stable/api/astropy.modeling.polynomial.Legendre1D.html
.. _astropy.modeling.models.CompoundModel: http://docs.astropy.org/en/stable/modeling/compound-models.html

"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import sys
import warnings
if sys.version > '3':
    long = int

import time
import os
import logging

import numpy
from scipy import interpolate, integrate, optimize
from scipy.special import erf
import astropy.constants
from astropy.modeling import FittableModel, Parameter

from ..par.parset import ParSet
from ..util.fileio import init_record_array
from ..util.instrument import spectrum_velocity_scale, resample_vector
from ..util.log import log_output
from ..util.constants import constants
from .spatiallybinnedspectra import SpatiallyBinnedSpectra
from .stellarcontinuummodel import StellarContinuumModel
from .pixelmask import PixelMask, SpectralPixelMask
from .spectralfitting import EmissionLineFit
from .emissionlinedb import EmissionLineDB
from .util import residual_growth

from matplotlib import pyplot

__author__ = 'Kyle B. Westfall'
# Add strict versioning
# from distutils.version import StrictVersion

from astropy.modeling import FittableModel, Parameter
from astropy.modeling.polynomial import Legendre1D

# BASE PROFILE DEFINITIONS ---------------------------------------------
class GaussianLineProfile(FittableModel):
    r"""
    Define a Gaussian line profile as parameterized by its zeroth
    moment, mean, and standard deviation:

    .. math::

        \mathcal{N}(x|f,\mu,\sigma) = \frac{f}{\sqrt{2\pi}\sigma}
        \exp\left(\frac{-(x-\mu)^2}{2\sigma^2}\right)

    The base class is `astropy.modeling.FittableModel`_, which
    facilitates its use in combining multiple components and other
    models in the `astropy.modeling`_ suite.
    """
    # Inputs the x positions
    inputs = ('x',)
    # Returns the renormalized Gaussian PDF at each x position
    outputs = ('y',)

    # Parameters are the peak value of the profile, the profile mean,
    # and the profile standard deviation (sigma)
    zmom=Parameter(default=1.0)
    mean=Parameter(default=0.0)
    sigma=Parameter(default=1.0)

    
    @staticmethod
    def evaluate(x, zmom, mean, sigma):
        # Discrete samples
        return zmom * numpy.exp(-0.5*numpy.square((x - mean)/sigma)) / numpy.sqrt(2*numpy.pi)/sigma
#        # Integrated over a pixel
#        diff = (x-mean)/pixelscale
#        denom = numpy.sqrt(2.0)*sigma/pixelscale
#        return zmom * (erf((diff+0.5)/denom) - erf((diff-0.5)/denom))/2.


    @staticmethod
    def fit_deriv(x, zmom, mean, sigma):
        y = GaussianLineProfile.evaluate(x,zmom,mean,sigma)
        d_zmom = y/zmom
        d_mean = y*(x-mean)/numpy.square(sigma)
        d_sigma = y*(numpy.square(x-mean)/numpy.power(sigma,3) - 1.0/sigma)
        return [d_zmom, d_mean, d_sigma]


    @staticmethod
    def flux(zmom, mean, sigma):
        return zmom


    @staticmethod
    def flux_err(zmom, mean, sigma, zmome, meane, sigmae):
        return zmome


    @staticmethod
    def moment(order, zmom, mean, sigma):
        if order == 0:
            return zmom
        if order == 1:
            return mean
        if order == 2:
            return sigma


    @staticmethod
    def moment_err(order, zmom, mean, sigma, zmome, meane, sigmae):
        if order == 0:
            return zmome
        if order == 1:
            return meane
        if order == 2:
            return sigmae


    @staticmethod
    def integral(zmom, mean, sigma):
        return integrate.quad( GaussianLineProfile.evaluate, -sigma, sigma, args=(zmom,mean,sigma))

    
    @staticmethod
    def scale_flux(zmom, mean, sigma, fac):
        return [zmom*fac, mean, sigma]


    @staticmethod
    def fix_flux():
        return [ True, False, False ]

    
    @staticmethod
    def mean_indx():
        return 1

    
    @staticmethod
    def shift_mean(zmom, mean, sigma, shift):
        return [zmom, mean+shift, sigma ]

    
    @staticmethod
    def fix_mean():
        return [ False, True, False ]

    
    @staticmethod
    def scale_stddev(zmom, mean, sigma, fac):
        return [zmom, mean, sigma*fac ]


    @staticmethod
    def fix_stddev():
        return [ False, False, True ]

    
#-----------------------------------------------------------------------

   
#-----------------------------------------------------------------------
class NCompLineProfile:
    """
    Construct a single line profile from many components with the same
    profile parameterization.
    """
    def __init__(self, ncomp, par=None, err=None, profile=GaussianLineProfile):
        # Check the type of the profile to use
        if not issubclass(profile, FittableModel):
            raise TypeError('Profile must be a subclass of FittableModel.')
        self.profile_class = profile

        # Ensure that the profile can be defined
        if not ncomp > 0:
            raise ValueError('Number of input components must be 1 or higher!')

        # Calculate and save the number of parameters
        self.ncomp = ncomp
        self.npar_per_component = len(self.profile_class.param_names)
        self.npar = self.ncomp*self.npar_per_component

        # Check the parameter vector, if provided
        if par is not None and len(par.shape) > 2:
            raise ValueError('Input parameter array should only be one or two dimensional.')
        _par = None if par is None else par.ravel()
        if _par is not None and len(_par) != self.npar:
            raise ValueError('Incorrect number of parameters provided')
        self.par = None if _par is None else _par.reshape(self.ncomp,-1).astype(float)

        # Check the parameter vector, if provided
        if err is not None and len(err.shape) > 2:
            raise ValueError('Input parameter error array should only be one or two dimensional.')
        _err = None if err is None else err.ravel()
        if _err is not None and len(_err) != self.npar:
            raise ValueError('Incorrect number of parameter errors provided')
        self.err = None if _err is None else _err.reshape(self.ncomp,-1).astype(float)

        # Build the profile function
        self.profile = self.profile_class() \
                if self.par is None else self.profile_class(*self.par[0])
        for i in range(1,ncomp):
            self.profile += self.profile_class() \
                    if self.par is None else self.profile_class(*self.par[i])


    def __call__(self, x, par=None):
        return self.sample(x,par)


    def _quick_sample(self, x):
        """
        Sample without providing/checking new input parameters.
        """
        return self.profile.evaluate(x, *self.par.ravel())


    def _mom1_integrand(self, x):
        return x*self._quick_sample(x)


    def _mom2_integrand(self, x):
        return x*x*self._quick_sample(x)
    

    def assign_par(self, par):
        _par = par.ravel() if len(par.shape) == 2 else par
        if len(_par) != self.npar:
            raise ValueError('Incorrect number of parameters provided')
        self.par = _par.copy().reshape(self.ncomp,-1)


    def assign_err(self, err):
        _err = err.ravel() if len(err.shape) == 2 else err
        if len(_err) != self.npar:
            raise ValueError('Incorrect number of parameters provided')
        self.err = _err.copy().reshape(self.ncomp,-1)


    def sample(self, x, par=None):
        if par is not None:
            self.assign_par(par)
        return self._quick_sample(x)


    def flux(self, par=None):
        if par is not None:
            self.assign_par(par)
        return numpy.sum( [ self.profile_class.flux(*p) \
                            for p in numpy.take(self.par,numpy.arange(self.ncomp),axis=0) ] )


    def flux_err(self, par=None, err=None):
        if par is not None:
            self.assign_par(par)
        if err is not None:
            self.assign_err(err)
        return numpy.sum([ self.profile_class.flux_err(*p, *e) \
                           for p,e in zip(numpy.take(self.par,numpy.arange(self.ncomp),axis=0),
                                          numpy.take(self.err,numpy.arange(self.ncomp),axis=0)) ] )


    def moment(self, order=0, par=None):
        """
        .. todo::
            impose some reasonable limits for the integrals
        """
        if par is not None:
            self.assign_par(par)

        if self.ncomp == 1:
            return self.profile_class.moment(order, *self.par[0,:])

        flux = self.flux(par=par)
        if order == 0:
            return flux

        mom1, _ = integrate.quad(self._mom1_integrand, -numpy.inf, numpy.inf)
        if order == 1:
            return mom1/flux
        mom2, _ = integrate.quad(self._mom2_integrand, -numpy.inf, numpy.inf)
        return numpy.sqrt(mom2/flux-numpy.square(mom1/flux))

    
    def moment_err(self, order=0, par=None, err=None):
        """
        .. todo::
            impose some reasonable limits for the integrals
        """
        if self.err is None and err is None:
            return 0.0

        if par is not None:
            self.assign_par(par)
        if err is not None:
            self.assign_err(err)

        if self.ncomp == 1:
            return self.profile_class.moment_err(order, *self.par[0,:], *self.err[0,:])
        raise NotImplementedError('Unable to calculate error for multiple components.')

    
    def scale_flux(self, fac):
        for i in range(self.ncomp):
            self.par[i,:] = self.profile_class.scale_flux(*self.par[i,:], fac)


    def set_flux(self, flx):
        self.scale_flux( flx/self.flux() )


    def fix_flux(self):
        return numpy.array(self.profile_class.fix_flux()*self.ncomp)


    def mean_indx(self):
        return [ self.profile_class.mean_indx()+i*self.ncomp for i in range(self.ncomp) ]
            

    def shift_mean(self, shift):
        for i in range(self.ncomp):
            self.par[i,:] = self.profile_class.shift_mean(*self.par[i,:], shift)

    
    def set_mean(self, mean):
        self.shift_mean( mean-self.moment(order=1) )

    
    def fix_mean(self):
        return numpy.array(self.profile_class.fix_mean()*self.ncomp)


    def scale_stddev(self, fac):
        if self.ncomp != 1:
            raise NotImplementedError('Cannot scale the standard deviation of multi-component profiles.')
        self.par[0,:] = self.profile_class.scale_stddev(*self.par[0,:], fac)


    def set_stddev(self, stddev):
        self.scale_stddev( stddev/self.moment(order=2) )


    def fix_stddev(self):
        if self.ncomp != 1:
            raise NotImplementedError('Cannot fix the standard deviation of multi-component profiles.')
        return numpy.array(self.profile_class.fix_stddev()*self.ncomp)

#-----------------------------------------------------------------------


# Main fitting engine --------------------------------------------------
class LineProfileFit:
    r"""
    Simultaneously fit multiple line profiles.  Currently only allows
    one to fit using an NCompLineProfile object.  The fitting algorithm
    used is `scipy.optimize.least_squares`_ with fitting method 'trf' to
    allow for bounds.

    .. todo::
        For multicomponent lines, set the first normalization to be the
        normalization for the sum of all components, then force the
        normalization of the subcomponents to be ordered from highest to
        lowest and bounded from 0 to 1.

    Args:
        x (1D array): Independent variable
        y (1D array): Dependent variable
        profile_list (list of or individual :class:`NCompLineProfile`):
            The profile(s) to fit to the dependent variable.
        base_order (int): (**Optional**) The order of the Legendre
            polynomial to include in the model for the baseline trend in
            y below the fitted line profile(s).
        error (1D array): (**Optional**) Error in the dependent
            variables.  If not provided, no error weighting is peformed
            during the fitting process, and the covariance *will not* be
            constructed.
        mask (1D array): (**Optional**) Boolean array used to ignore
            values in `y` during the fit.
        par (1D array): (**Optional**) Initial guess for model
            parameters.  The number of parameters much match the
            expectation based on the provided list of profiles and the
            order of the baseline polynomial.  If not provided, the
            parameters are initialized to 0.
        fixed (1D array): (**Optional**) Flags used to fix parameters
            during the fit.  The number of parameters much match the
            expectation based on the provided list of profiles and the
            order of the baseline polynomial.  If not provided, all
            parameters are freely fit.
        bounds (2-tuple): (**Optional**) Tuple with two array-like
            elements giving the upper and lower bound for each
            parameter.  The length of each array element must match the
            number of parameters.  For an unbounded problem, set
            ``bounds=None``, or use numpy.inf with an appropriate sign
            to disable bounds on all or some variables.
        run_fit (bool): (**Optional**) Flag to run the fit upon
            instantiation of the object, which defaults to True.  If set
            to False, the object is initialized but the fit is not
            executed, and can be executed later using :func:`fit`.
        construct_covariance (bool): (**Optional**) Flag to construct
            the covariance matrix based on the result object provided by
            `scipy.optimize.least_squares`_, which defaults to True.  If
            set to False, the covariance matrix is set to None.
        verbose (int): (**Optional**) Verbosity level for
            `scipy.optimize.least_squares`_; default is 0.

    Attributes:
        x (numpy.ndarray): Independent variable of length :math:`M`.
        y (numpy.ndarray): Dependent variable to be fit of length
            :math:`M`.
        err (numpy.ndarray): Error, :math:`\sigma`, in the dependent
            variable of length :math:`M`.
        mask (numpy.ndarray): Flag to fit dependent variables of length
            :math:`M`.
        nlines (int): Number of lines being fit.
        base_order (int): The order of the Legendre polynomial include
            in the model; see `astropy.modeling.polynomial.Legendre1D`_.
        model (`astropy.modeling.models.CompoundModel`_): The compound
            model being fit composed of the line profiles and the
            baseline (if a baseline is being fit).  That is, this
            defines the function :math:`f(\mathbf{x}|\mathbf{\theta})`,
            where :math:`\mathbf{x}` is the dependent variable and
            :math:`\mathbf{\theta}` is the list of variables.  The
            fitting algorithm minimizes the (error-weighted) residuals
            to approximate :math:`y=f(\mathbf{x}|\mathbf{\theta})`.
        npar (int): Total number of parameters in the model.  The is the
            number of parameters per line and the number of parameters
            included in the baseline.
        nfitpar (int): Number of *free* parameters, defined as
            :math:`N`.
        par (numpy.ndarray): Full list of model parameters, including
            those parameters that have been fixed.
        fixed (numpy.ndarray): Flags to fit (False) or fix (True) a
            given parameter.
        bounds (2-tuple): Tuple with two array-like
            elements giving the upper and lower bound for each
            parameter.  For an unbounded problem, this is set to
            ``bounds=(-numpy.inf,numpy.inf)``.
        result (`scipy.optimize.OptimizeResult`_): Object with the
            results from `scipy.optimize.least_squares`_.  The
            best-fitting parameters, :math:`\mathbf{\theta}`, is
            returned as ``result.x``.
        cov (numpy.ndarray): The formal covariance matrix for the fit.
            The `scipy.optimize.OptimizeResult`_ object provides the
            Jacobian of the model, an :math:`M \times N` array with
            elements

            .. math::

                J_{ij} = \left.\frac{\partial f_i}{\partial
                \theta_j}\right|_{\mathbf{\theta}}

            at location in parameter space of the best-fitting model.
            This is used to construct the covariance matrix by taking
            the inverse of the curvature matrix:

            .. math::

                \mathbf{\alpha}_{kl} = \left[\frac{1}{\sigma_{i}}
                J_{ik}\right]^{\rm T} \left[\frac{1}{\sigma_{i}}
                J_{il}\right].

            That is, :math:`\mathbf{C} = \mathbf{\alpha}^{-1}`.

    Raises:
        TypeError: Raised if the provided profile objects are not
            instances of :class:`NCompLineProfile`.
        ValueError: Raised if any of the provided parameter arrays (par,
            fixed) are not one-dimensional or the number of parameters
            is not as expected based on the number of profile and
            baseline parameters.

    .. todo::
        - Provide a better initialization for the parameters.
        - Need to provide the best-fitting parameters and errors as full
          vectors, including the fixed parameters.

    """
    def __init__(self, x, y, profile_list, base_order=None, error=None, mask=None, par=None,
                 fixed=None, bounds=None, run_fit=True, construct_covariance=True, verbose=0):

        # Check the input list of profiles
        _profile_list = profile_list if isinstance(profile_list, list) else [profile_list]
        self.nlines = len(_profile_list)
        for i in range(self.nlines):
            if not isinstance(_profile_list[i], NCompLineProfile):
                raise TypeError('LineProfileFit only works with NCompLineProfile objects.')

        # Construct and check the parameter vectors
        self.npar = numpy.sum([ p.npar for p in _profile_list])
        self.base_order = -1 if base_order is None else base_order
        if self.base_order > -1:
            self.npar += (self.base_order + 1)
        if par is not None and len(par.shape) != 1:
            raise ValueError('Input parameter array should only be one-dimensional.')
        if fixed is not None and fixed.shape != par.shape:
            raise ValueError('Input parameter and fixed flags should have the same shape.')
        _par = None if par is None else par.ravel()
        if _par is not None and len(_par) != self.npar:
            raise ValueError('Incorrect number of parameters provided')
        self.par = numpy.zeros(self.npar, dtype=numpy.float) \
                            if _par is None else _par.astype(numpy.float)
        self.fixed = numpy.zeros(self.npar, dtype=numpy.bool) \
                            if fixed is None else fixed.astype(numpy.bool)
        if not isinstance(bounds, tuple):
            raise ValueError('Bounds must have type tuple.')
        if len(bounds) != 2:
            raise ValueError('Bounds must only have two elements.')
        if len(bounds[0]) != self.npar or len(bounds[1]) != self.npar:
            raise ValueError('Lower and/or upper bounds have an incorrect number of parameters.')
        self.bounded = bounds is not None
        self.bounds = (-numpy.inf, numpy.inf) if bounds is None else bounds

        # Get the number of free parameters
        self.nfitpar = numpy.sum(numpy.invert(self.fixed))

        # Construct the model
        self.model = _profile_list[0].profile
        for i in range(1,self.nlines):
            self.model += _profile_list[i].profile
        if self.base_order > -1:
            self.model += Legendre1D(self.base_order)
        # Fix the parameters
        if numpy.sum(self.fixed) > 0:
            self.model.fixed[ self.model.param_names[fixed] ] = True

#        print(self.model)

        # Need to setup how to tie fluxes, velocities, dispersions
        self.x = None
        self.y = None
        self.err = None
        self.mask = None
        self.result = None
        self.cov = None

        if run_fit:
            self.fit(x, y, error=error, mask=mask, construct_covariance=construct_covariance,
                     verbose=verbose)


    def _assign_par(self, par):
        """
        Fixed parameters are ignored!
        """
        if len(par) not in [ self.npar, self.nfitpar ]:
            raise ValueError('Incorrect number of parameters provided')
        self.par[~(self.fixed)] = par[~(self.fixed)] if len(par) == self.npar else par
        # Then tie parameters...

        
    def _chi(self, par):
        return self._resid(par)/self.err


    def _resid(self, par):
        self._assign_par(par)
        ymodel = self._quick_sample(self.x)
        return (self.y-ymodel)


    def _calculate_covariance_matrix(self, verbose=0):
        if self.result is None or self.err is None:
            if self.result is not None:
                warnings.warn('Cannot calculate covariance matrix without error vector!')
            return None
#        invsig = numpy.array([numpy.ma.power(self.err,-1.)]*self.result.jac.shape[1]).T
#        a = (invsig.ravel()*self.result.jac.ravel()).reshape(self.result.jac.shape)
#        a = numpy.dot(a.T,a)
        a = numpy.dot(self.result.jac.T,self.result.jac)
        if verbose > 0:
            print('Curvature:')
            print('    ', a)
        try:
            return numpy.linalg.inv(a)
        except:
            return numpy.linalg.pinv(a)


    def _quick_sample(self, x):
        """
        Sample without providing/checking new input parameters.
        """
        #print(*self.par.ravel())
        return self.model.evaluate(x, *self.par.ravel())


    def sample(self, x, par=None):
        if par is not None:
            self._assign_par(par)
        return self._quick_sample(x)

    
#    def d_chi(self, par):
#        self._assign_par(par)
#        model_flux = self._quick_sample(self.wave[~(self.mask)])
#        d_model_flux = numpy.array(self.model.fit_deriv(self.wave[~(self.mask)],*self.par))
#        return numpy.array([ -dm/self.error[~(self.mask)] \
#                             for dm in numpy.take(d_model_flux,numpy.arange(self.par), axis=0) ])


    def fit(self, x, y, error=None, mask=None, construct_covariance=True, verbose=0):
        """
        Fit the line profiles provided upon initialization to the data
        using `scipy.optimize.least_squares`_.
        """
        # Check that the lengths of all the arrays match
        if len(x.shape) != 1:
            raise ValueError('Can only fit single spectra!')
        if x.shape != y.shape \
                or (error is not None and error.shape != y.shape) \
                or (mask is not None and mask.shape != y.shape):
            raise ValueError('All vector shapes must be identical.')

        # If the input is a masked array, combine the two masks
        _mask = numpy.ma.getmaskarray(y)
        if mask is not None:
            _mask |= mask

        # Set the internal vectors to be numpy.ndarrays that only include the unmasked values
        self.x = numpy.array(x[~_mask])
        self.y = numpy.array(y[~_mask])
        self.err = None if error is None else numpy.array(error[~_mask])

        if self.y.size <= self.npar:
            raise ValueError('Insufficient data points ({0}) to fit with this function ' \
                             '(npar={1})!'.format(self.y.size, self.npar))

#        print('Initial guess parameters inside LineProfileFit:')
#        print(self.par)
#        print(self.bounds)

        try:
            self.result = optimize.least_squares(self._resid if self.err is None else self._chi,
                                                 self.par, bounds=self.bounds, jac='2-point',
                                                 method='trf', loss='linear', verbose=verbose)
        except ValueError as e:
            print('ValueError: {0}'.format(e))
            print(self.par)
            print(self.bounds)
            self.result.success = False
            return

#        try:
#            self.result = optimize.least_squares(self._resid if self.err is None else self._chi,
#                                                 self.par, bounds=self.bounds, jac='2-point',
#                                                 #method='lm', loss='linear', verbose=1)
#                                                 method='trf', loss='linear', verbose=2)
#        except:
#            exc_info = sys.exc_info()
#            warnings.warn('Error during fit: {0}:: {1}'.format(exc_info[0], exc_info[1]))
#            self.result = OptimizeResult(x=self.par, cost=None, fun=None, jac=None, grad=None,
#                                         optimality=None, active_mask=None, nfev=None, njev=None,
#                                         status=None, success=False)
#            return

        self.cov = None
        if construct_covariance:
            self.cov = self._calculate_covariance_matrix(verbose=verbose)

        if verbose > 0:
            print('Result: ')
            print('    ', self.result.x)
            if self.cov is not None:
                print('Covariance: ')
                print('    ', self.cov)


class ElricPar(ParSet):
    def __init__(self, emission_lines, base_order, window_buffer, guess_redshift, guess_dispersion,
                 minimum_snr=None, pixelmask=None, stellar_continuum=None):

        arr_in_fl = [ numpy.ndarray, list, int, float ]
        in_fl = [ int, float ]

        pars =     [ 'emission_lines', 'base_order', 'window_buffer', 'guess_redshift',
                     'guess_dispersion', 'minimum_snr', 'pixelmask', 'stellar_continuum' ]
        values =   [ emission_lines, base_order, window_buffer, guess_redshift, guess_dispersion,
                     minimum_snr, pixelmask, stellar_continuum ]
        defaults = [ None, -1,   25.0, None, None, 0.0, None, None ]
        dtypes =   [ EmissionLineDB, int, in_fl, arr_in_fl, arr_in_fl, in_fl, PixelMask,
                     StellarContinuumModel ]

        ParSet.__init__(self, pars, values=values, defaults=defaults, dtypes=dtypes)


    def toheader(self, hdr):
        hdr['LPBASEO'] = (self['base_order'], 'Baseline Legendre polynomial order')
        hdr['LPWIN'] = (self['window_buffer'], 'Buffer for fitting window (ang)')
        return hdr


    def fromheader(self, hdr):
        self['base_order'] = hdr['LPBASEO']
        self['window_buffer'] = hdr['LPWIN']



class TiedLineProfile:
    def __init__(self, profile, base_profile, relative_flux=None, mean_separation=None,
                 relative_stddev=None):
        self.profile = profile
        self.base_profile = base_profile
        self.relative_flux = relative_flux
        self.mean_separation = mean_separation
        self.relative_stddev = relative_stddev


    def tie_flux(self):
        self.profile.set_flux(self.relative_flux*self.base_profile.flux())


    def tie_mean(self):
        self.profile.set_mean(self.base_profile.moment(order=1)+self.velocity_separation)


    def tie_stddev(self):
        self.profile.set_stddev(self.relative_stddev*self.base_profile.moment(order=2))


class ElricFittingWindow:
    r"""

    A utility class for the fitting windows used by Elric

    The class can be instantiated as fully None.

    Args:
        nlines (int): (**Optional**) The number of lines to be fit
            simultaneously.
        db_indx (numpy.ndarray): (**Optional**) An array with the
            database index (0-based) of each line to be fit.
        line_index (numpy.ndarray): (**Optional**) An array with the
            index numbers, read from the database, of each line to be
            fit.
        restwave (numpy.ndarray): (**Optional**) An array with the rest
            wavelengths of each line to be fit.
        profile_set (numpy.ndarray): (**Optional**) An array with the
            profile objects that define the functional form of each
            line.
        fixed_par (numpy.ndarray): (**Optional**) An array with the
            fixed parameters for *all* parameters in the model to be
            fit.
        bounds (numpy.ndarray): (**Optional**) A two-column array with
            the lower (first column) and upper (second column) bounds
            for the fit parameters.
        log_bounds (numpy.ndarray): (**Optional**) The range of the
            boundary should be considered as logarithmic when testing if
            a parameter is near its boundary.
        output_model (bool): (**Optional**) Include the best-fitting
            model in the composite emission-line model for each
            spectrum.  This is only flagged as true if ALL the emission
            lines in the fitting window are to be included according to
            the emission-line database.
        tied_pairs (numpy.ndarray): (**Optional**) The series of
            tied profiles (:class:`TiedLineProfile` objects) that are
            used to tie parameters of the model.
        tied_funcs (numpy.ndarray): (**Optional**) The member functions
            of the :class:`TiedLineProfile` objects that should be
            called *sequentially* to tie model parameters.
    """
    def __init__(self, nlines=None, db_indx=None, line_index=None, restwave=None, profile_set=None,
                 fixed_par=None, bounds=None, log_bounds=None, output_model=None, tied_pairs=None,
                 tied_funcs=None):

        self.nlines = nlines
        self.db_indx = db_indx
        self.line_index = line_index
        self.restwave = restwave
        self.profile_set = profile_set
        self.init_par = numpy.array([ p.par.copy().ravel() for p in self.profile_set ])
        self.fixed_par = fixed_par
        self.bounds = bounds
        self.log_bounds = log_bounds
        self.output_model = output_model
        self.tied_pairs = tied_pairs
        self.tied_funcs = tied_funcs


    def append(self, db_indx, line_index, restwave, profile, fixed_par, bounds, log_bounds,
               output_model):
        """
        Append another line to the fitting window.  Must be a single line!
        """

        # Check that input is for a single line
        if not isinstance(db_indx, int) or not isinstance(line_index, (int, numpy.int64)):
            raise TypeError('Appended index must be a single integer.')
        if not isinstance(restwave, float):
            raise TypeError('Appended restwavelength must be a single float.')
        if not isinstance(profile, NCompLineProfile):
            raise TypeError('Appended profile must have type NCompLineProfile.')
        if len(fixed_par.shape) != 1 or len(fixed_par) != profile.npar:
            raise TypeError('Appended fixed parameter must be a vector with the correct length.')
        if len(bounds.shape) != 2 or bounds.shape[0] != profile.npar or bounds.shape[1] != 2:
            raise TypeError('Appended bounds must have shape: {0}'.format(tuple(profile.npar,2)))
        if len(log_bounds.shape) != 1 or len(log_bounds) != profile.npar:
            raise TypeError('Appended log_bounds parameter must be a vector of the correct length.')
        if not isinstance(output_model, bool):
            raise TypeError('Flag to output model must be boolean.')
        
        self.nlines += 1
        self.db_indx = numpy.append(self.db_indx, db_indx)
        self.line_index = numpy.append(self.line_index, line_index)
        self.restwave = numpy.append(self.restwave, restwave)
        self.profile_set = numpy.append(self.profile_set, profile)
        self.init_par = numpy.append(self.init_par, profile.par.copy(), axis=0)
        self.fixed_par = numpy.append(self.fixed_par, fixed_par.astype(bool).reshape(1,-1), axis=0)
        self.bounds = numpy.append(self.bounds, bounds.reshape(1,-1,2), axis=0)
        self.log_bounds = numpy.append(self.log_bounds, log_bounds.astype(bool).reshape(1,-1),
                                       axis=0)
        self.output_model &= output_model


    def reset_init_mean(self, indx, newmean):
        self.profile_set[indx].set_mean(newmean)
        self.init_par[indx,:] = self.profile_set[indx].par


    def reinit_profiles(self):
        for i in range(self.nlines):
            self.profile_set[i].assign_par(self.init_par[i,:])


class Elric(EmissionLineFit):
    """
    ELRIC: Emission-Line Regression and Inference Class

    https://en.wikipedia.org/wiki/Edward_Elric

    Use LineProfileFit to fit the emission-line properties in a set of
    spectra.

    .. todo::
        - Implement some scheme to penalize multicomponent fits at low S/N

    """
    def __init__(self, bitmask, wave=None, flux=None, emission_lines=None, error=None, mask=None,
                 stellar_continuum=None, base_order=-1, window_buffer=25, guess_redshift=None,
                 guess_dispersion=None, default_dispersion=20.0, run_fit=False, loggers=None,
                 quiet=False):

        EmissionLineFit.__init__(self, 'elric', bitmask)
        # Attributes kept by SpectralFitting:
        #   fit_type='emission_line', bitmask=bitmask, par=None
        # Attributes kept by EmissionLineFit:
        #   fit_method='elric'

        # Declare the attributes kept for convenience
        self.default_dispersion = default_dispersion

        self.loggers=None
        self.quiet=None

        self.nspec = None
        self.nwave = None
        self.emission_lines = None
        self.neml = None

        self.nwindows = None
        self.fitting_window = None

        self.wave = None
        self.sres = None
        self.flux = None
        self.fluxnc = None
        self.error = None
        self.mask = None
        self.continuum = None
        self.no_continuum = None

        self.base_order = None
        self.window_buffer = None
        self.guess_redshift = None
        self.guess_dispersion = None

        self.bestfit = None

        if not run_fit:
            return

        if wave is None or flux is None or emission_lines is None:
            raise ValueError('Cannot run fit without wavelength and flux vectors, and ' \
                             'emission-line database.')

        self.fit(wave, flux, emission_lines, error=error, mask=mask,
                 stellar_continuum=stellar_continuum, base_order=base_order,
                 window_buffer=window_buffer, guess_redshift=guess_redshift,
                 guess_dispersion=guess_dispersion, loggers=loggers, quiet=quiet)


    @staticmethod
    def _find_tied_index(index, emission_lines):
        i = numpy.where(emission_lines['index'] == index)[0]
        if len(i) > 1:
            raise ValueError('Line indices are not unique!')
        if emission_lines['mode'][i[0]] == 'f':
            return i[0], index
        tied = int(emission_lines['mode'][i[0]][1:])
        j = numpy.where(emission_lines['index'] == tied)[0]
        if len(j) > 1:
            raise ValueError('Line indices are not unique!')
        return (j[0], tied) if emission_lines['mode'][j[0]] == 'f' \
                    else Elric._find_tied_index(tied, emission_lines)

    
    @staticmethod
    def _set_profile_ties(base_profiles, base_restwave, base_fixed_par, tied_profiles,
                          tied_restwave, tied_fixed_par, mode, flux):
        nprof = len(base_profiles)
        if len(base_restwave) != nprof or len(tied_profiles) != nprof \
                or len(tied_restwave) != nprof or len(mode) != nprof or len(flux) != nprof:
            raise ValueError('Must provide same length arrays to _tied_funcs.')
        
        tied_pairs = []
        tied_functions = []
        for i in range(nprof):
            if mode[i][0] == 'w':
                continue 
            tied_pars += [ TiedLineProfile(tied_profiles[i], base_profiles[i]) ]
            tied_fixed_par[i,:] |= base_fixed_par[i,:]
            if mode[i][0] in [ 'a', 'x' ]:
                tied_pairs[i].relative_flux = flux[i]
                tied_functions += [ tied_pairs[i].tie_flux ]
                tied_fixed_par[i,:] |= tied_pairs[i].profile.fix_flux()
            if mode[i][0] in [ 'a', 'k', 'v' ]:
                tied_pairs[i].mean_separation \
                    = (tied_restwave[i]/base_restwave[i]-1)*astropy.constants.c.to('km/s').value
                tied_functions += [ tied_pairs[i].tie_mean ]
                tied_fixed_par[i,:] |= tied_pairs[i].profile.fix_mean()
            if mode[i][0] in [ 'a', 'x', 'k', 's' ]:
                tied_pairs[i].relative_stddev = 1.0
                tied_functions += [ tied_pairs[i].tie_stddev ]
                tied_fixed_par[i,:] |= tied_pairs[i].profile.fix_stddev()

        return numpy.array(tied_pairs), numpy.array(tied_functions), tied_fixed_par


    @staticmethod
    def _line_velocity_offset(restwave, restwave_primary):
        return (restwave/restwave_primary - 1.0) * astropy.constants.c.to('km/s').value


    def _parse_emission_line_models(self):
        r"""

        Parse the input emission-line file into a set of windows --- the
        number of windows is :math:`N_{\rm win}`) --- that are fit for
        each spectrum.
        """

        # Setup the number of windows to fit based on the primary lines
        primary_line = self.emission_lines['mode'] == 'f'
        primary_index = self.emission_lines['index'][primary_line]
        self.nwindows = numpy.sum(primary_line)
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO,
                       'Number of emission lines to fit: {0}'.format(self.neml))
            log_output(self.loggers, 1, logging.INFO,
                       'Number of fitting windows: {0}'.format(self.nwindows))

        # Set up the primary line data using the guess parameters from
        # the database
        self.fitting_window = numpy.empty(self.nwindows, dtype=object)
        _db_indx = numpy.arange(self.neml)[primary_line]
        for i in range(self.nwindows):

            db_indx = numpy.array([_db_indx[i]])
            line_index = numpy.array([primary_index[i]])
            restwave = numpy.array([ self.emission_lines['restwave'][primary_line][i] ])
            profile_set = numpy.array([ NCompLineProfile(
                                        self.emission_lines['ncomp'][primary_line][i],
                                        par=self.emission_lines['par'][primary_line][i],
                                profile=eval(self.emission_lines['profile'][primary_line][i])) ])
            fixed_par = self.emission_lines['fix'][primary_line][i].astype(bool).reshape(1,-1)
            bounds = numpy.array([ [l,u] \
                                    for l,u in zip(self.emission_lines['lobnd'][primary_line][i],
                                                   self.emission_lines['hibnd'][primary_line][i])
                                 ]).reshape(1,-1,2)
            log_bounds = self.emission_lines['log_bnd'][primary_line][i].astype(bool).reshape(1,-1)
            output_model = self.emission_lines['output_model'][primary_line][i]
            self.fitting_window[i] = ElricFittingWindow(nlines=1, db_indx=db_indx,
                                                        line_index=line_index, restwave=restwave,
                                                        profile_set=profile_set,
                                                        fixed_par=fixed_par, bounds=bounds,
                                                        log_bounds=log_bounds,
                                                        output_model=output_model)

        # Return if nothing is tied
        if self.nwindows == self.neml:
            return

        # Append additional profiles
        for i, e in enumerate(self.emission_lines):
            if primary_line[i]:
                continue
            # Basic check
            if e['index'] == int(e['mode'][1:]):
                raise ValueError('Cannot tie line to itself!')

            # Run up the depenency chain to find the primary line
            # associated with this one
            indx, tied_index = self._find_tied_index(e['index'], self.emission_lines)
            base_indx = numpy.where( primary_index == tied_index)[0][0]

            # Append the line to this fitting window
            self.fitting_window[base_indx].append(i, e['index'], e['restwave'],
                                                  NCompLineProfile(e['ncomp'], par=e['par'],
                                                                   profile=eval(e['profile'])),
                                                  e['fix'],
                    numpy.array([ [l,u] for l,u in zip(e['lobnd'], e['hibnd']) ]).reshape(-1,2),
                                                  e['log_bnd'], bool(e['output_model']))
            # Set the guess velocity of the profile based on the
            # velocity relative to the main line
            self.fitting_window[base_indx].reset_init_mean(-1,
                            self._line_velocity_offset(self.fitting_window[base_indx].restwave[-1],
                                                       self.fitting_window[base_indx].restwave[0]))

        # Finalize by setting up any tied parameters
        for i in range(self.nwindows):
            if self.fitting_window[i].nlines == 1:
                continue
            # Basic debug check
            if self.emission_lines['mode'][self.fitting_window[i].db_indx[0]] != 'f':
                raise ValueError('DEBUG')

            # The order is important.  Start with the primary line:
            tied = numpy.zeros(self.fitting_window[i].nlines, dtype=numpy.bool)
            tied[0] = True              # Primary line already "tied" to the set
            # Iterate until all lines are tied
            while numpy.sum(tied) != self.fitting_window[i].nlines:
                # Find the indices of lines that are directly tied to
                # any previously tied lines.  At the beginning, this
                # only selects lines that are tied to the primary line.
                indx = numpy.where([
                            not tied[j] and numpy.any(int(mode[1:]) \
                                            == self.fitting_window[i].line_index[tied]) \
                            for j,mode in enumerate(
                                self.emission_lines['mode'][self.fitting_window[i].db_indx[:]])
                                   ])[0]

                # For these lines, find the indices in the list of lines
                # for this window with their associated lines to tie to
                tied_indx = [ int(mode[1:]) for mode in \
                                self.emission_lines['mode'][self.fitting_window[i].db_indx[indx]]]
                profile_indx = [ numpy.where(t == self.fitting_window[i].line_index[:])[0][0] \
                                    for t in tied_indx ]

                # Run the maing function that sets up the seried of
                # tying classes/functions and updates the parameters
                # that should be "fixed" during the fit.
                _tied_pairs, _tied_funcs, _fixed_par \
                        = self._set_profile_ties(self.fitting_window[i].profile_set[profile_indx],
                                                 self.fitting_window[i].restwave[profile_indx],
                                numpy.take(self.fitting_window[i].fixed_par, profile_indx, axis=0),
                                                 self.fitting_window[i].profile_set[indx],
                                                 self.fitting_window[i].restwave[indx],
                                numpy.take(self.fitting_window[i].fixed_par, indx, axis=0),
                                self.emission_lines['mode'][self.fitting_window[i].db_indx[indx]],
                                self.emission_lines['flux'][self.fitting_window[i].db_indx[indx]])

                # Add these to the main output
                self.fitting_window[i].tied_pairs = numpy.append(self.fitting_window[i].tied_pairs,
                                                                 _tied_pairs)
                self.fitting_window[i].tied_funcs = numpy.append(self.fitting_window[i].tied_funcs,
                                                                 _tied_funcs)
                for j, jj in enumerate(indx):
                    self.fitting_window[i].fixed_par[jj,:] = _fixed_par[j,:]

                # These lines are now tied!
                tied[indx] = True


    @staticmethod
    def _velocity_vectors(wave, fitting_window):
        return numpy.array([ (wave/fw.restwave[0]-1.0)*astropy.constants.c.to('km/s').value \
                            for fw in fitting_window ])


    @staticmethod
    def _fit_masks(wave, fitting_window, redshift, window_buffer):

        # Input redshift must be a single value
        if not isinstance(redshift, float):
            raise TypeError('redshift must be a single float')

        nwindows = len(fitting_window)
        fitting_mask = numpy.array([ numpy.logical_and(
                                        wave>(fw.restwave[0]-window_buffer)*(1+redshift),
                                        wave<(fw.restwave[0]+window_buffer)*(1+redshift)) \
                                     for fw in fitting_window ])
#        print(fitting_mask.shape)

        for i in range(nwindows):
            if len(fitting_window[i].restwave) == 1:
                continue
            fitting_mask[i,:] |= numpy.any( numpy.array([ numpy.logical_and(
                                                        wave>(rw-window_buffer)*(1+redshift),
                                                        wave<(rw+window_buffer)*(1+redshift)) \
                                                for rw in fitting_window[i].restwave ]), axis=0)
        return fitting_mask


    @staticmethod
    def _check_db(emission_lines):
        neml = emission_lines.nsets
        for i in range(neml):
            profile = eval(emission_lines['profile'][i])
            npar = len(profile.param_names)
            if emission_lines['par'][i].size != npar*emission_lines['ncomp'][i]:
                raise ValueError('Provided {0} parameters, but expected {1}.'.format(
                                  emission_lines['par'][i].size, npar*emission_lines['ncomp'][i]))
            if emission_lines['fix'][i].size != npar*emission_lines['ncomp'][i]:
                raise ValueError('Provided {0} fix flags, but expected {1}.'.format(
                                  emission_lines['fix'][i].size, npar*emission_lines['ncomp'][i]))
            if numpy.any([f not in [0, 1] for f in emission_lines['fix'][i] ]):
                warnings.warn('Fix values should only be 0 or 1; non-zero values interpreted as 1.')
            if emission_lines['lobnd'][i].size != npar*emission_lines['ncomp'][i]:
                raise ValueError('Provided {0} lower bounds, but expected {1}.'.format(
                                  emission_lines['lobnd'][i].size,
                                  npar*emission_lines['ncomp'][i]))
            if emission_lines['hibnd'][i].size != npar*emission_lines['ncomp'][i]:
                raise ValueError('Provided {0} upper bounds, but expected {1}.'.format(
                                  emission_lines['hibnd'][i].size,
                                  npar*emission_lines['ncomp'][i]))
            if emission_lines['log_bnd'][i].size != npar*emission_lines['ncomp'][i]:
                raise ValueError('Provided {0} log boundaries designations, but expected '
                                 '{1}.'.format(emission_lines['log_bnd'][i].size,
                                               npar*emission_lines['ncomp'][i]))


    @staticmethod
    def _is_near_bound(par, lbnd, ubnd, logbounded, rtol=1e-2, atol=1e-4):
        """
        Determine of any of the parameters within start and start+npar
        are "near" an imposed boundary.
        """
#        print('par: ', par)
#        print('lbnd: ', lbnd)
        lbounded = ~numpy.isinf(lbnd)
#        print('L bounded: ', lbounded)

#        print('ubnd: ', ubnd)
        ubounded = ~numpy.isinf(ubnd)
#        print('U bounded: ', ubounded)

        if numpy.sum(lbounded) == 0 and numpy.sum(ubounded) == 0:
#            print('No bounds')
            return False

        if numpy.any(lbounded & ~ubounded & (par - lbnd < atol)):
#            print('Near lower bound')
            return True

        if numpy.any(~lbounded & ubounded & (ubnd - par < atol)):
#            print('Near upper bound')
            return True

        ulbounded = lbounded & ubounded
#        print('UL bounded', ulbounded)
        if numpy.sum(ulbounded) == 0:
#            print('None bounded from both sides')
            return False

#        print('log bounded', logbounded)
#        print('UL & log bounded', ulbounded & logbounded)
#        print('UL & not log bounded', ulbounded & ~logbounded)
        Dp = numpy.zeros(ulbounded.size, dtype=numpy.float)
        Dp[ulbounded & logbounded] = numpy.log10(ubnd[ulbounded & logbounded]) \
                                        - numpy.log10(lbnd[ulbounded & logbounded])
        Dp[ulbounded & ~logbounded] = ubnd[ulbounded & ~logbounded] - lbnd[ulbounded & ~logbounded]
        tol = Dp*rtol
#        print(tol)
#        print(tol[ulbounded & logbounded])
#        print((numpy.ma.log10(par) - numpy.ma.log10(lbnd))[ulbounded & logbounded])
#        print((numpy.ma.log10(ubnd) - numpy.ma.log10(par))[ulbounded & logbounded])

        if numpy.any( ulbounded & logbounded & ( (numpy.ma.log10(par) - numpy.ma.log10(lbnd) < tol) 
                                        | (numpy.ma.log10(ubnd) - numpy.ma.log10(par) < tol))):
#            print('Near boundary in log')
            return True
                            
#        print(tol[ulbounded & ~logbounded])
#        print((par - lbnd)[ulbounded & ~logbounded])
#        print((ubnd - par)[ulbounded & ~logbounded])

        if numpy.any( ulbounded & ~logbounded & ( (par - lbnd < tol) | (ubnd - par < tol))):
#            print('Near boundary in linear')
            return True
#        print('Not near boundary')
        return False


#        # Check if the parameters are near a boundary
#        if self.bestfit[i,j].bounded:
#            model_fit_par['LOBND'][i,j,:npar] = self.bestfit[i,j].bounds[0]
#            model_fit_par['UPBND'][i,j,:npar] = self.bestfit[i,j].bounds[1]
#            for k in range(npar):
#                if not numpy.isinf(self.bestfit[i,j].bounds[0][k]) \
#                        and model_fit_par['PAR'][i,j,k] - model_fit_par['ERR'][i,j,k] \
#                                    <= self.bestfit[i,j].bounds[0][k]:
#                    near_bound = True
#                    break
#                if not numpy.isinf(self.bestfit[i,j].bounds[1][k]) \
#                        and model_fit_par['PAR'][i,j,k] + model_fit_par['ERR'][i,j,k] \
#                                    >= self.bestfit[i,j].bounds[1][k]:
#                    near_bound = True
#                    break

#            if self.bestfit[i,j].bounded:
#                if numpy.any(~numpy.isinf(self.bestfit[i,j].bounds[0][pl:pl+nlinepar[k]]) \
#                                & (model_fit_par['PAR'][i,j,pl:pl+nlinepar[k]]
#                                    - model_fit_par['ERR'][i,j,pl:pl+nlinepar[k]]
#                                    <= self.bestfit[i,j].bounds[0][pl:pl+nlinepar[k]])):
#                    model_eml_par['MASK'][i,emlj] \
#                            = self.bitmask.turn_on(model_eml_par['MASK'][i,emlj], 'NEAR_BOUND')
#                    near_bound = True
#                if numpy.any(~numpy.isinf(self.bestfit[i,j].bounds[1][pl:pl+nlinepar[k]]) \
#                                & (model_fit_par['PAR'][i,j,pl:pl+nlinepar[k]]
#                                    + model_fit_par['ERR'][i,j,pl:pl+nlinepar[k]] \
#                                    >= self.bestfit[i,j].bounds[1][pl:pl+nlinepar[k]])):
#                    model_eml_par['MASK'][i,emlj] \
#                            = self.bitmask.turn_on(model_eml_par['MASK'][i,emlj], 'NEAR_BOUND')
#                    near_bound = True


    @staticmethod
    def _correct_subeml_par(restwave_0, restwave_k, par, err=None):
        """
        Correct the parameters fit assuming a rest wavelength of
        restwave_0, when it's actually restwave_k.  Input parameters are
        expected to be flux, velocity, velocity dispersion.
        """
        c = astropy.constants.c.to('km/s').value
        _par = par.copy()
        _par[0] *= restwave_k/restwave_0
        _par[1] = (restwave_0*(par[1]/c + 1)/restwave_k - 1)*c
        _par[2] *= restwave_0/restwave_k
        if err is None:
            return _par

        _err = err.copy()
        _err[0] *= restwave_k/restwave_0
        _err[1] *= restwave_0/restwave_k
        _err[2] *= restwave_0/restwave_k
        return _par, _err
        

    def _assess_and_save_fit(self, i, j, model_fit_par, model_eml_par):
        """
        Assess the result of the LineProfileFit results.
            x Check the success failure
            x Calculate the chi2, rchi2, rms, fractional_rms, residuals,
              fractional_residuals
            x Compare the fit parameters to the bounds
            - Check the chi-square and residuals?
            - Check the velocity offset wrt the input
            - For multiple lines, check the order of the lines matches
              the guess parameters
            - For multiple component lines, check the ordering of the
              subcomponents
        """
        include_model = True
        if not self.bestfit[i,j].result.success:
            model_fit_par['MASK'][i,j] = self.bitmask.turn_on(model_fit_par['MASK'][i,j],
                                                              'FIT_FAILED')

        # Assign the fitting window for each emission line
        model_eml_par['WIN_INDEX'][i,self.fitting_window[j].db_indx] = j
#        print(model_eml_par['WIN_INDEX'][i,:])

        # Get the full list of parameters (including fixed and tied values)
        npar = self.bestfit[i,j].npar
        self.bestfit[i,j]._assign_par(self.bestfit[i,j].result.x)
        par = self.bestfit[i,j].par
        model_fit_par['PAR'][i,j,:npar] = par[:]
        model_fit_par['FIXED'][i,j,:npar] = self.bestfit[i,j].fixed
        model_fit_par['FIXED'][i,j,npar:] = True
        model_fit_par['ERR'][i,j,:npar] = numpy.zeros(npar, dtype=numpy.float)
        if self.error is not None:
            tmperr = numpy.diag(self.bestfit[i,j].cov)
            if ~numpy.any(tmperr < 0):
                model_fit_par['ERR'][i,j,~(model_fit_par['FIXED'][i,j])] \
                        = numpy.sqrt(numpy.diag(self.bestfit[i,j].cov))
            else:
                model_fit_par['MASK'][i,j] = self.bitmask.turn_on(model_fit_par['MASK'][i,j],
                                                                  'UNDEFINED_COVAR')


#        print('PAR: ', model_fit_par['PAR'][i,j,:npar])
#        print('ERR: ', model_fit_par['ERR'][i,j,:npar])

        # If the parameters are bounded, save the bounds
        if self.bestfit[i,j].bounded:
            model_fit_par['LOBND'][i,j,:npar] = self.bestfit[i,j].bounds[0]
            model_fit_par['UPBND'][i,j,:npar] = self.bestfit[i,j].bounds[1]

        # Get the per-line parameters
        near_bound = False
        nlinepar = numpy.array([ p.npar for p in self.fitting_window[j].profile_set ])
        for k in range(self.fitting_window[j].nlines):
            pl = 0 if k == 0 else numpy.sum(nlinepar[:k])
            emlj = self.fitting_window[j].db_indx[k]
            model_eml_par['FLUX'][i,emlj] \
                    = self.fitting_window[j].profile_set[k].moment(order=0,
                                                                   par=par[pl:pl+nlinepar[k]])

            model_eml_par['FLUXERR'][i,emlj] \
                    = self.fitting_window[j].profile_set[k].moment_err(order=0,
                                                                       par=par[pl:pl+nlinepar[k]],
                                                    err=model_fit_par['ERR'][i,j,pl:pl+nlinepar[k]])
            # Flux unit conversion
#            print(astropy.constants.c.to('km/s').value/self.fitting_window[j].restwave[k])
#                    /= astropy.constants.c.to('km/s').value/self.fitting_window[j].restwave[k]
            model_eml_par['FLUX'][i,emlj] \
                    *= self.fitting_window[j].restwave[0]/astropy.constants.c.to('km/s').value
            model_eml_par['FLUXERR'][i,emlj] \
                    *= self.fitting_window[j].restwave[0]/astropy.constants.c.to('km/s').value

            model_eml_par['KIN'][i,emlj,0] \
                    = self.fitting_window[j].profile_set[k].moment(order=1,
                                                                   par=par[pl:pl+nlinepar[k]])
            model_eml_par['KINERR'][i,emlj,0] \
                    = self.fitting_window[j].profile_set[k].moment_err(order=1,
                                                                       par=par[pl:pl+nlinepar[k]],
                                                    err=model_fit_par['ERR'][i,j,pl:pl+nlinepar[k]])

            model_eml_par['KIN'][i,emlj,1] \
                    = self.fitting_window[j].profile_set[k].moment(order=2,
                                                                   par=par[pl:pl+nlinepar[k]])

            model_eml_par['KINERR'][i,emlj,1] \
                    = self.fitting_window[j].profile_set[k].moment_err(order=2,
                                                                       par=par[pl:pl+nlinepar[k]],
                                                    err=model_fit_par['ERR'][i,j,pl:pl+nlinepar[k]])

            # Correct for the fact that the fit was done using a
            # velocity vector defined by the primary line
            # TODO: Need to change how fitting is done so that this is
            # NOT necessary!
            if k > 0:
                line_par, line_parerr \
                        = self._correct_subeml_par(self.fitting_window[j].restwave[0],
                                                   self.fitting_window[j].restwave[k],
                                                   numpy.array([ model_eml_par['FLUX'][i,emlj],
                                                                 model_eml_par['KIN'][i,emlj,0],
                                                                 model_eml_par['KIN'][i,emlj,1] ]),
                                                   err= numpy.array([
                                                            model_eml_par['FLUXERR'][i,emlj],
                                                            model_eml_par['KINERR'][i,emlj,0],
                                                            model_eml_par['KINERR'][i,emlj,1] ]))
                model_eml_par['FLUX'][i,emlj], model_eml_par['KIN'][i,emlj,0], \
                    model_eml_par['KIN'][i,emlj,1] = tuple(line_par)
                model_eml_par['FLUXERR'][i,emlj], model_eml_par['KINERR'][i,emlj,0], \
                    model_eml_par['KINERR'][i,emlj,1] = tuple(line_parerr)

            if self.bestfit[i,j].bounded \
                    & self._is_near_bound(model_fit_par['PAR'][i,j,pl:pl+nlinepar[k]],
                                          self.bestfit[i,j].bounds[0][pl:pl+nlinepar[k]],
                                          self.bestfit[i,j].bounds[1][pl:pl+nlinepar[k]],
                            self.fitting_window[j].log_bounds.ravel().copy()[pl:pl+nlinepar[k]]):
                model_eml_par['MASK'][i,emlj] \
                            = self.bitmask.turn_on(model_eml_par['MASK'][i,emlj], 'NEAR_BOUND')
                near_bound = True

#        print('Any line in window near boundary: ', near_bound)
        if near_bound:
            model_fit_par['MASK'][i,j] = self.bitmask.turn_on(model_fit_par['MASK'][i,j],
                                                              'NEAR_BOUND')

        # Save which parameters were fixed or should be ignored
        model_fit_par['IGNORE'][i,j,npar:] = True
        
        # Set the fit statistics
        model_fit_par['CHI2'][i,j] = 0.0 if self.error is None else \
                                            numpy.sum(numpy.square(self.bestfit[i,j]._chi(par)))
#        print('CHI2:', model_fit_par['CHI2'][i,j])
        model_fit_par['RCHI2'][i,j] = model_fit_par['CHI2'][i,j] \
                                            / (model_fit_par['NPIXFIT'][i,j] - npar)
#        print('RCHI2:', model_fit_par['RCHI2'][i,j])
        resid = self.bestfit[i,j]._resid(par)
        model_fit_par['RMS'][i,j] = numpy.sqrt(numpy.mean(numpy.square(resid)))
#        print('RMS:', model_fit_par['RMS'][i,j])

        model_fit_par['RESID'][i,j] = residual_growth(resid, [0.25, 0.50, 0.75, 0.90, 0.99])
#        print('RESID: ', model_fit_par['RESID'][i,j])

        fit = self.bestfit[i,j].sample(self.bestfit[i,j].x, par)
        indx = numpy.absolute(fit) > 1e-4
        if numpy.sum(indx) > 0:
            frac_resid = resid[indx]/fit[indx]
            model_fit_par['FRAC_RMS'][i,j] = numpy.sqrt(numpy.mean(numpy.square(frac_resid)))
#            print('FRMS:', model_fit_par['FRAC_RMS'][i,j])
            if numpy.sum(indx) > 1:
                model_fit_par['FRAC_RESID'][i,j] = residual_growth(frac_resid,
                                                               [0.25, 0.50, 0.75, 0.90, 0.99])
#                print('FRAC_RESID: ', model_fit_par['FRAC_RESID'][i,j])

        return near_bound


    def _instrumental_dispersion_correction(self, model_eml_par):
        """
        Determine the instrumental velocity dispersion at the centroids
        of the fitted emission lines.
        """
        if self.sres is None:
            return
        interpolator = interpolate.interp1d(self.wave, self.sres, fill_value='extrapolate',
                                            assume_sorted=True)
        restwave = numpy.array([self.emission_lines['restwave']]*self.nspec)
        cnst = constants()

        # Get the instrumental dispersion at each (valid) line center
        indx = ~self.bitmask.flagged(model_eml_par['MASK'],
                                     flag=['INSUFFICIENT_DATA', 'FIT_FAILED', 'NEAR_BOUND',
                                           'UNDEFINED_COVAR' ])
        model_eml_par['SINST'][indx] = astropy.constants.c.to('km/s') \
                                            / interpolator((model_eml_par['KIN'][indx,0] 
                                                        / astropy.constants.c.to('km/s').value+1.0)
                                                                *restwave[indx])/cnst.sig2fwhm

#        z = model_eml_par['KIN'][indx,0] / astropy.constants.c.to('km/s').value
#        pyplot.scatter(restwave[indx]*(1.0+z), model_eml_par['SINST'][indx], marker='.', s=30,
#                       color='k')
#        pyplot.show()
#
#        pyplot.scatter(model_eml_par['SINST'][indx], model_eml_par['KIN'][indx,1], marker='.',
#                       s=30, color='k')
#        pyplot.show()


    def _correct_velocity_dispersion(self, model_eml_par):
        """
        Use the **previously calculated** instrumental velocity
        dispersion to correct the measured velocity dispersions of the
        lines to the astrophysical velocity dispersion.
        """

        # Correct the second moments for the instrumental dispersion
        defined = numpy.zeros(indx.shape, dtype=numpy.bool)
        defined[indx] = model_eml_par['SINST'][indx] < model_eml_par['KIN'][indx,1]
        nonzero = numpy.zeros(indx.shape, dtype=numpy.bool)
        nonzero[indx] = numpy.absolute(model_eml_par['SINST'][indx] 
                                            - model_eml_par['KIN'][indx,1]) > 0

#        orig = model_eml_par['KIN'].copy()

        model_eml_par['KINERR'][nonzero,1] \
                    = 2.0*model_eml_par['KIN'][nonzero,1]*model_eml_par['KINERR'][nonzero,1]
        model_eml_par['KIN'][nonzero,1] = numpy.square(model_eml_par['KIN'][nonzero,1]) \
                                                - numpy.square(model_eml_par['SINST'][nonzero])
        model_eml_par['KIN'][nonzero,1] = model_eml_par['KIN'][nonzero,1] \
                            / numpy.sqrt(numpy.absolute(model_eml_par['KIN'][nonzero,1]))
        model_eml_par['KINERR'][nonzero,1] /= numpy.absolute(2.*model_eml_par['KIN'][nonzero,1])

        # Flag undefined values
        model_eml_par['MASK'][indx & ~defined] \
                = self.bitmask.turn_on(model_eml_par['MASK'][indx & ~defined], 'UNDEFINED_SIGMA')

#        flg = self.bitmask.flagged(model_eml_par['MASK'], flag='UNDEFINED_SIGMA')
#        pyplot.scatter(orig[nonzero & flg,1], model_eml_par['KIN'][nonzero & flg,1], marker='.',
#                       s=30, color='0.8')
#        pyplot.scatter(orig[nonzero & ~flg,1], model_eml_par['KIN'][nonzero & ~flg,1], marker='.',
#                       s=30, color='k')
#        pyplot.show()


    def fit_SpatiallyBinnedSpectra(self, binned_spectra, par=None, loggers=None, quiet=False):
        """

        This is a basic interface that is geared for the DAP that
        interacts with the rest of the, more general, parts of the
        class.

        This should not declare anything to self!

        """
        # Assign the parameters if provided
        if par is None:
            raise ValueError('Required parameters for LineProfileFit have not been defined.')

        # SpatiallyBinnedSpectra object always needed
        if binned_spectra is None:
            raise ValueError('Must provide spectra object for fitting.')
        if not isinstance(binned_spectra, SpatiallyBinnedSpectra):
            raise TypeError('Must provide a valid SpatiallyBinnedSpectra object!')
        if binned_spectra.hdu is None:
            raise ValueError('Provided SpatiallyBinnedSpectra object is undefined!')

        # StellarContinuumModel object only used when accounting for
        # underlying absorption
        if par['stellar_continuum'] is not None:
            if not isinstance(par['stellar_continuum'], StellarContinuumModel):
                raise TypeError('Provided stellar continuum must have StellarContinuumModel type!')
            if par['stellar_continuum'].hdu is None:
                raise ValueError('Provided StellarContinuumModel is undefined!')

        # Get the data arrays to fit
        wave = binned_spectra['WAVE'].data
        sres = binned_spectra['SPECRES'].data
        flux = binned_spectra.copy_to_masked_array(flag=binned_spectra.do_not_fit_flags())
        ivar = binned_spectra.copy_to_masked_array(ext='IVAR',
                                                        flag=binned_spectra.do_not_fit_flags())

        good_spec = (binned_spectra['BINS'].data['SNR'] > par['minimum_snr']) \
                        & ~(numpy.array([ b in binned_spectra.missing_bins \
                                                for b in numpy.arange(binned_spectra.nbins)]))
        _flux = flux[good_spec,:]
        noise = numpy.ma.power(ivar[good_spec,:], -0.5)
        guess_redshift = par['guess_redshift'][good_spec]
        guess_dispersion = par['guess_dispersion'][good_spec]

        if par['stellar_continuum'] is not None:
            # Get the models for the binned spectra
            continuum = par['stellar_continuum'].copy_to_masked_array(
                        flag=par['stellar_continuum'].all_except_emission_flags())[good_spec,:]
        else:
            continuum = None

        model_wave, model_flux, model_base, model_mask, model_fit_par, model_eml_par \
            = self.fit(wave, _flux, par['emission_lines'], error=noise, mask=par['pixelmask'],
                       sres=sres, stellar_continuum=continuum, base_order=par['base_order'],
                       window_buffer=par['window_buffer'], guess_redshift=guess_redshift,
                       guess_dispersion=guess_dispersion, loggers=loggers, quiet=quiet)
        
        _model_flux = numpy.zeros(flux.shape, dtype=numpy.float)
        _model_flux[good_spec,:] = model_flux

        _model_base = numpy.zeros(flux.shape, dtype=numpy.float)
        _model_base[good_spec,:] = model_base

        _model_mask = numpy.zeros(flux.shape, dtype=self.bitmask.minimum_dtype())

        _model_mask[numpy.ma.getmaskarray(flux)] \
                = self.bitmask.turn_on(_model_mask[numpy.ma.getmaskarray(flux)], 'DIDNOTUSE')
        bad_snr = (binned_spectra['BINS'].data['SNR'] < par['minimum_snr']) \
                        & ~(numpy.array([ b in binned_spectra.missing_bins \
                                                for b in numpy.arange(binned_spectra.nbins)]))
        _model_mask[bad_snr,:] = self.bitmask.turn_on(_model_mask[bad_snr,:], 'LOW_SNR')

        _model_mask[good_spec,:] = model_mask

        _model_fit_par = init_record_array(flux.shape[0], model_fit_par.dtype)
        _model_fit_par[good_spec] = model_fit_par
        _model_fit_par['BIN_INDEX'] = numpy.arange(flux.shape[0])

        _model_eml_par = init_record_array(flux.shape[0], model_eml_par.dtype)
        _model_eml_par[good_spec] = model_eml_par
        _model_eml_par['BIN_INDEX'] = numpy.arange(flux.shape[0])

        return model_wave, _model_flux, _model_base, _model_mask, _model_fit_par, _model_eml_par


    def fit(self, wave, flux, emission_lines, error=None, mask=None, sres=None,
            stellar_continuum=None, base_order=-1, window_buffer=25, guess_redshift=None,
            guess_dispersion=None, loggers=None, quiet=False):
        """
        The flux array is expected to have size Nspec x Nwave.

        Raises:
            ValueError: Raised if the length of the spectra, errors, or
                mask does not match the length of the wavelength array;
                raised if the wavelength, redshift, or dispersion arrays
                are not 1D vectors; and raised if the number of
                redshifts or dispersions is not a single value or the
                same as the number of input spectra.
        """
        # Check the input emission-line database
        self._check_db(emission_lines)
        self.emission_lines = emission_lines
        self.neml = self.emission_lines.nsets

        # Initialize the reporting
        if loggers is not None:
            self.loggers = loggers
        self.quiet = quiet
        
        # Check the input
        if len(wave.shape) != 1:
            raise ValueError('Input wavelengths must be a vector; all flux vectors should have' \
                             'the same wavelength solution.')
        self.wave = wave
        self.nwave = len(self.wave)

        self.sres = None
        if sres is not None:
            if len(sres.shape) != 1:
                raise ValueError('Input spectral resolution must be a vector; all flux vectors ' \
                                 'should have the same spectral resolution.')
            if len(sres) != self.nwave:
                raise ValueError('Spectral resolution vector length must match the wavelength ' \
                                 'vector.')
            self.sres = sres

        _flux = flux.reshape(1,-1) if len(flux.shape) != 2 else flux
        if _flux.shape[1] != self.nwave:
            raise ValueError('The spectra have to have the same length as the wavelength vector.')
        if error is not None and error.shape != _flux.shape:
            raise ValueError('The shape of any provided error array must match the flux array.')

        # Set the input to masked arrays, if they aren't already
        self.flux = _flux if isinstance(_flux,numpy.ma.MaskedArray) else numpy.ma.MaskedArray(_flux)
        self.nspec = self.flux.shape[0]
        # Include the mask, if provided
        if mask is not None:
            if isinstance(mask, numpy.ndarray) and numpy.issubdtype(mask.dtype,bool):
                if mask is not None and mask.shape != self.flux.shape:
                    raise ValueError('The shape of the mask array must match the flux array.')
                self.flux[mask] = numpy.ma.masked
            if isinstance(mask, SpectralPixelMask):
                self.flux[mask.boolean(self.wave,nspec=self.nspec)] = numpy.ma.masked
        # Do the same with the error array
        self.error = None if error is None else (error if isinstance(error, numpy.ma.MaskedArray) \
                            else numpy.ma.MaskedArray(error, mask=numpy.ma.getmaskarray(self.flux)))

        # Setup the redshift and dispersion vectors
        if guess_redshift is not None:
            self.guess_redshift = numpy.atleast_1d(guess_redshift)
            if len(self.guess_redshift.shape) != 1:
                raise ValueError('Input guess redshifts must be a single vector.')
            if self.guess_redshift.size not in [1, self.nspec]:
                raise ValueError('Provide single redshift or a redshift for each spectrum.')
            if self.guess_redshift.size == 1:
                self.guess_redshift = numpy.full(self.nspec, guess_redshift)
        else:
            self.guess_redshift = numpy.zeros(self.nspec, dtype=numpy.float)

        if guess_dispersion is not None:
            self.guess_dispersion = numpy.asarray(guess_dispersion)
            if len(self.guess_dispersion.shape) != 1:
                raise ValueError('Input guess dispersions must be a single vector.')
            if self.guess_dispersion.size not in [1, self.nspec]:
                raise ValueError('Provide single dispersion or a dispersion for each spectrum.')
            if self.guess_dispersion.size == 1:
                self.guess_dispersion = numpy.full(self.nspec, guess_dispersion)
        else:
            self.guess_dispersion = numpy.full(self.nspec, self.default_dispersion,
                                               dtype=numpy.float)

        # Subtract the stellar continuum if it has been provided
        if stellar_continuum is not None:
            if stellar_continuum.shape != self.flux.shape:
                raise ValueError('Shape of the stellar continuum array must match the flux array.')
            self.continuum = stellar_continuum
            self.fluxnc = self.flux - self.continuum
            if isinstance(self.continuum, numpy.ma.MaskedArray):
                # Get where the stellar-continuum models are masked but
                # the binned spectra are not
                self.no_continuum = numpy.invert(numpy.ma.getmaskarray(self.flux)) \
                                            & numpy.ma.getmaskarray(self.continuum)
                self.fluxnc.mask[self.no_continuum] = False
            else:
                self.no_continuum = None
        else:
            self.fluxnc = self.flux
            self.continuum = None
            self.no_continuum = None

        # No continuum for any spaxel!
        if self.no_continuum is None:
            self.no_continuum = numpy.ones(self.flux.shape, dtype=numpy.bool)

        # Get the spectra to use for the profile fitting
#        pyplot.step(self.wave, self.flux[0,:], where='mid', color='k', lw=0.5, linestyle='-')
#        pyplot.step(self.wave, self.fluxnc[0,:], where='mid', color='g', lw=0.5, linestyle='-')
#        pyplot.plot(self.wave, self.no_continuum[0,:], color='r', lw=0.5, linestyle='-')
#        pyplot.show()

        # Build the emission-line fitting windows
        self._parse_emission_line_models()
        max_npar = numpy.amax([ fw.fixed_par.size for fw in self.fitting_window ])
        self.base_order = base_order
        if self.base_order > -1:
            max_npar += base_order+1
        self.window_buffer = window_buffer

        # Initialize the output arrays
        #  Model flux:
        model_flux = numpy.zeros(self.flux.shape, dtype=numpy.float)
        model_base = numpy.zeros(self.flux.shape, dtype=numpy.float)
        #  Model mask:
        model_mask = numpy.zeros(self.flux.shape, dtype=self.bitmask.minimum_dtype())
        indx = numpy.ma.getmaskarray(self.flux)
        model_mask[indx] = self.bitmask.turn_on(model_mask[indx], 'DIDNOTUSE')
        model_mask[self.no_continuum] = self.bitmask.turn_on(model_mask[self.no_continuum],
                                                             'NOCONTINUUM')

        # Model parameters and fit quality
        model_fit_par = init_record_array(self.nspec,
                                          self._per_fitting_window_dtype(self.nwindows, max_npar,
                                                                self.bitmask.minimum_dtype()))
        model_fit_par['BIN_INDEX'] = numpy.arange(self.nspec)

        # Emission-line parameters
        model_eml_par = init_record_array(self.nspec,
                                          self._per_emission_line_dtype(self.neml, 2,
                                                                self.bitmask.minimum_dtype()))
        model_eml_par['BIN_INDEX'] = numpy.arange(self.nspec)

        t = time.clock()
        # Do the fit
        velocity = self._velocity_vectors(self.wave, self.fitting_window)
        self.bestfit = numpy.empty((self.nspec,self.nwindows), dtype=object)
        for i in range(self.nspec):

#            for j in range(self.nwindows):
#                print('Window: {0}/{1}'.format(j+1, self.nwindows))
#                for p in self.fitting_window[j].profile_set:
#                    print(p.par)
            if not self.quiet:
                log_output(self.loggers, 1, logging.INFO, 'Fit: {0}/{1}'.format(i+1,self.nspec))

            # Get the fitting mask for each emission-line in each window
            fitting_mask = self._fit_masks(self.wave, self.fitting_window, self.guess_redshift[i],
                                           self.window_buffer)

            # Set the mask for any pixels that are NEVER fit for this
            # spectrum
            indx = ~numpy.any(fitting_mask, axis=0)
            model_mask[i,indx] = self.bitmask.turn_on(model_mask[i,indx], 'OUTSIDE_RANGE')

            model_jump = numpy.array([
                        numpy.sum(fm & self.no_continuum[i,:]) not in [0, numpy.sum(fm)] \
                        for fm in numpy.take(fitting_mask, numpy.arange(self.nwindows), axis=0) ])
            spec_to_fit = numpy.ma.array([ self.flux[i,:] if mj else self.fluxnc[i,:] \
                                            for mj in model_jump ])

            cz = self.guess_redshift[i]*astropy.constants.c.to('km/s').value

#            total_near_bound = 0
            # Fit each window
            for j in range(self.nwindows):
                # Rest the profile to the initial parameters
                self.fitting_window[j].reinit_profiles()
#                print('Fitting: ', self.emission_lines['name'][self.fitting_window[j].db_indx])

                # No observed pixels in fitting window
#                print('Sum of fitting mask: {0}'.format(numpy.sum(fitting_mask[j,:])))
                if numpy.sum(fitting_mask[j,:]) == 0:
                    warnings.warn('No valid data for fit to line(s) {0} ({1}).  '
                                  'Continuuing'.format(self.fitting_window[j].line_index,
                                    self.emission_lines['name'][self.fitting_window[j].db_indx]))
                    model_fit_par['MASK'][i,j] = self.bitmask.turn_on(model_fit_par['MASK'][i,j],
                                                                      'INSUFFICIENT_DATA')
                    continue

#                w,h = pyplot.figaspect(1)
#                fig = pyplot.figure(figsize=(1.5*w,1.5*h))
#                ax = fig.add_axes([0.1, 0.3, 0.8, 0.4])
#                ax.set_xlim(velocity[j,0], velocity[j,-1])
#                ax.step(velocity[j,:], spec_to_fit[j,:], where='mid', color='k', linestyle='-',
#                        lw=0.5)
#                ax.plot(velocity[j,:], fitting_mask[j,:], color='r')
#                axt = ax.twiny()
#                axt.set_xlim(self.wave[0], self.wave[-1])
#                pyplot.show()

                # Shift the guess velocity based on the redshift and
                # adjust the guess flux
                # TODO: Use emission-line moment results for the
                # guesses, if provided
                for k,p in enumerate(self.fitting_window[j].profile_set):
                    p.shift_mean(cz)
                    vi = numpy.argsort(numpy.absolute(velocity[j,:]-p.moment(order=1)))[0]
#                    print('guess flux: ', (spec_to_fit[j,vi] if spec_to_fit[j,vi] > 0.1 else 0.1) \
#                                * numpy.sqrt(2*numpy.pi) * p.moment(order=2))
                    p.set_flux((spec_to_fit[j,vi] if spec_to_fit[j,vi] > 0.1 else 0.1) \
                                * numpy.sqrt(2*numpy.pi) * p.moment(order=2))

                # Setup the guess parameters and bounds
                _guess_par = numpy.array([ p.par.copy().ravel() \
                                            for p in self.fitting_window[j].profile_set ]).ravel()
                _fixed_par = self.fitting_window[j].fixed_par.ravel().copy()
                _bounds = self.fitting_window[j].bounds.copy().reshape(-1,2)
                print(_guess_par)
                print(_fixed_par)
                print(_bounds)

                # For individual lines, bound the line mean to be within
                # the fitting window
                if self.fitting_window[j].nlines == 1:
                    indx = self.fitting_window[j].profile_set[0].mean_indx()
                    _bounds[indx,0] = velocity[j,fitting_mask[j,:]][0]
                    _bounds[indx,1] = velocity[j,fitting_mask[j,:]][-1]
                # For multiple lines, bound the line mean to be within
                # its own region of the window defined by the +/- 2/3 of
                # the velocity separation to the nearest line.  These
                # borders overlap!
                else:
                    indx = numpy.array(self.fitting_window[j].profile_set[0].mean_indx())
                    for k in range(1,self.fitting_window[j].nlines):
                        indx = numpy.append(indx,
                                numpy.array(self.fitting_window[j].profile_set[0].mean_indx()) \
                            + sum([self.fitting_window[j].profile_set[kk].npar for kk in range(k)]))
                    srt = numpy.argsort(_guess_par[indx])
                    borders = []
                    for k in range(len(indx)-1):
                        diff = _guess_par[indx][srt[k+1]]-_guess_par[indx][srt[k]]
                        borders += [ _guess_par[indx][srt[k]] + 2.0*diff/3.0,
                                     _guess_par[indx][srt[k]] + diff/3.0 ]
                    borders = numpy.array([ velocity[j,fitting_mask[j,:]][0], *borders,
                                            velocity[j,fitting_mask[j,:]][-1] ])
                    _bounds[indx,0] = borders[2*srt]
                    _bounds[indx,1] = borders[2*srt+1]

                # Add the paramters for the baseline
                if base_order > -1:
                    _guess_par = numpy.append(_guess_par, numpy.zeros(base_order+1))
                    _fixed_par = numpy.append(self.fitting_window[j].fixed_par.ravel(),
                                              numpy.zeros(base_order+1).astype(bool))
                    _bounds = numpy.append(_bounds,
                          numpy.array([-numpy.inf,numpy.inf]*(base_order+1)).reshape(-1,2),axis=0)

                # Check there is sufficient data to perform the fit
                model_fit_par['NPIXFIT'][i,j] = len(spec_to_fit[j,fitting_mask[j,:]].compressed())
                if model_fit_par['NPIXFIT'][i,j] <= _guess_par.size:
                    warnings.warn('Insufficient data points ({0}) to fit with this function ' \
                                  '(npar={1})!'.format(
                                        len(spec_to_fit[j,fitting_mask[j,:]].compressed()),
                                        _guess_par.size))
                    self.bestfit[i,j] = None
                    # Flag that the fit was not performed
                    model_fit_par['MASK'][i,j] = self.bitmask.turn_on(model_fit_par['MASK'][i,j],
                                                                      'INSUFFICIENT_DATA')
                    continue

                # Make sure that the guess is still in the bounds:
                # TODO: How do I deal with this:
                #   - throw a warning and adjust the guess
                #   - as above and include a mask bit
                #   - don't perform the fit because it will be bogus
                indx = numpy.logical_or( _guess_par < _bounds[:,0], _guess_par > _bounds[:,1])
                if numpy.sum(indx) > 0:
                    warnings.warn('Initial guess outside bounds for line(s) {0} ({1}).  '
                                  'Adjusted to center of bounds.'.format(
                                    self.fitting_window[j].line_index,
                                    self.emission_lines['name'][self.fitting_window[j].db_indx]))
                    _guess_par[indx] = numpy.mean(_bounds[indx,:], axis=1)

                # TODO: Get rid of this debugging issue...
                if ~numpy.all(_bounds[:,0]<_bounds[:,1]):

                    print('input')
                    print(_bounds)
                    print(self.fitting_window[j].bounds.reshape(-1,2))
                    print(numpy.sum(fitting_mask[j,:]))

                    if self.fitting_window[j].nlines == 1:
                        indx = self.fitting_window[j].profile_set[0].mean_indx()
                        _bounds[indx,0] = velocity[j,fitting_mask[j,:]][0]
                        _bounds[indx,1] = velocity[j,fitting_mask[j,:]][-1]
                    else:
                        indx = numpy.array(self.fitting_window[j].profile_set[0].mean_indx())
                        for k in range(1,self.fitting_window[j].nlines):
                            indx = numpy.append(indx,
                                numpy.array(self.fitting_window[j].profile_set[0].mean_indx()) \
                                + sum([self.fitting_window[j].profile_set[kk].npar
                                    for kk in range(k)]))
                        srt = numpy.argsort(_guess_par[indx])
                        borders = []
                        for k in range(len(indx)-1):
                            diff = _guess_par[indx][srt[k+1]]-_guess_par[indx][srt[k]]
                            borders += [ _guess_par[indx][srt[k]] + 2.0*diff/3.0,
                                        _guess_par[indx][srt[k]] + diff/3.0 ]
                        borders = numpy.array([ velocity[j,fitting_mask[j,:]][0], *borders,
                                                velocity[j,fitting_mask[j,:]][-1] ])
                        print('borders')
                        print(borders)
                        _bounds[indx,0] = borders[2*srt]
                        _bounds[indx,1] = borders[2*srt+1]

                    print('output')
                    print(_bounds)
                    raise ValueError('WTF')

#                print(_guess_par)
#                print(_bounds)
#                pyplot.errorbar(velocity[j,fitting_mask[j,:]], spec_to_fit[j,fitting_mask[j,:]],
#                                yerr=_error[i,fitting_mask[j,:]], marker='.', color='k',
#                                linestyle='', capsize=0)
#                pyplot.title('Spec: {0}; Line (set): {1}'.format(i+1,j+1))
#                pyplot.show()
#                print('N valid pixels: {0}'.format(
#                                               len(spec_to_fit[j,fitting_mask[j,:]].compressed())))
                
                # Run the fit
#                print('guess par: {0} {1}'.format(i, j))
#                print(_guess_par)
#                print(' ')
                self.bestfit[i,j] = LineProfileFit(velocity[j,:], spec_to_fit[j,:],
                                                   self.fitting_window[j].profile_set.tolist(),
                                            error=None if self.error is None else self.error[i,:],
                                                   base_order=base_order, mask=~(fitting_mask[j,:]),
                                                   par=_guess_par,
                                                   bounds=(_bounds[:,0], _bounds[:,1]),
                                                   construct_covariance=True)
                                                   #, verbose=2)#, fixed=_fixed_par)
#                print('Spec: {0}; Line (set): {1}; Success: {2}'.format(i+1,j+1,
#                                                                self.bestfit[i,j].result.success))

                near_bound = self._assess_and_save_fit(i, j, model_fit_par, model_eml_par)
                if near_bound:
#                    print('Flagging mask as near bound')
                    model_mask[i,fitting_mask[j,:]] \
                            = self.bitmask.turn_on(model_mask[i,fitting_mask[j,:]], 'NEAR_BOUND')
#                    total_near_bound += 1

                # Add the best-fit to the lines to the composite model
                # for this spectrum
                if self.fitting_window[j].output_model and not self.bestfit[i,j].result.success:
                    model_mask[i,fitting_mask[j,:]] \
                            = self.bitmask.turn_on(model_mask[i,fitting_mask[j,:]], 'FIT_FAILED')
                elif not self.fitting_window[j].output_model or model_fit_par['MASK'][i,j] != 0:
                    model_fit_par['MASK'][i,j] = self.bitmask.turn_on(model_fit_par['MASK'][i,j],
                                                                      'EXCLUDED_FROM_MODEL')
                elif self.fitting_window[j].output_model and model_fit_par['MASK'][i,j] == 0:
                    # Get the full model (eventually will only be the
                    # baseline)
                    p = self.bestfit[i,j].result.x.copy()
                    model_base[i,fitting_mask[j,:]] \
                            += self.bestfit[i,j].sample(velocity[j,fitting_mask[j,:]], par=p)
                    # Remove the baseline and only get the line profile
                    p[-base_order-1:] = 0.0
                    model_flux[i,fitting_mask[j,:]] \
                            += self.bestfit[i,j].sample(velocity[j,fitting_mask[j,:]], par=p)

#                if 44 in self.fitting_window[j].line_index:
#                    pyplot.errorbar(velocity[j,fitting_mask[j,:]], spec_to_fit[j,fitting_mask[j,:]],
#                                    yerr=self.error[i,fitting_mask[j,:]], marker='.', color='k',
#                                    linestyle='', capsize=0)
#                    pyplot.plot(velocity[j,fitting_mask[j,:]],
#                                self.bestfit[i,j].sample(velocity[j,fitting_mask[j,:]],
#                                par=self.bestfit[i,j].result.x), linestyle='-', color='g')
#                    pyplot.plot(velocity[j,fitting_mask[j,:]],
#                                self.bestfit[i,j].sample(velocity[j,fitting_mask[j,:]],
#                                par=_guess_par), linestyle='-', color='r')
#                    print(_guess_par)
#                    print(self.bestfit[i,j].result.x)
#                    pyplot.title('Spec: {0}; Line (set): {1}'.format(i+1,j+1))
#                    pyplot.show()

#            pyplot.plot(wave, model_flux[i,:], linestyle='-', color='g')
#            pyplot.plot(wave, model_base[i,:], linestyle='-', color='r')
#            pyplot.show()

#            print('Windows with parameters near boundary: {0}'.format(total_near_bound))

        # Remove the lines from the full model to just provide the
        # baseline model
        model_base -= model_flux

#        pyplot.plot(wave, model_flux[0,:], linestyle='-', color='g')
#        pyplot.plot(wave, model_base[0,:], linestyle='-', color='r')
#        pyplot.show()

        # Determine the instrumental dispersion at the line centers
        self._instrumental_dispersion_correction(model_eml_par)

        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, 'Fits completed in {0:.4e} min.'.format(
                       (time.clock() - t)/60))

        return self.wave, model_flux, model_base, model_mask, model_fit_par, model_eml_par


