# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Implements an emission-line profile fitting class.

.. warning::

    Although it may still be possible to use :class:`Elric`, it is
    currently not actively used by the survey-level operation of the
    MaNGA DAP. See instead :class:`mangadap.proc.sasuke.Sasuke`.

----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

import os
import time
import logging

import numpy
from scipy import interpolate, integrate, optimize
from scipy.special import erf

import astropy.constants
from astropy.modeling.polynomial import Legendre1D

from ..par.parset import KeywordParSet
from ..par.emissionlinedb import EmissionLineDB
from ..util.fileio import init_record_array
from ..util.log import log_output
from ..util.pixelmask import SpectralPixelMask
from ..util import lineprofiles
from .spatiallybinnedspectra import SpatiallyBinnedSpectra
from .stellarcontinuummodel import StellarContinuumModel
from .spectralfitting import EmissionLineFit
from .util import sample_growth

from matplotlib import pyplot

# Add strict versioning
# from distutils.version import StrictVersion

# Main fitting engine --------------------------------------------------
class LineProfileFit:
    r"""
    Simultaneously fit multiple line profiles. Currently only allows
    one to fit using an
    :class:`~mangadap.util.lineprofiles.NCompLineProfile` object. The
    fitting algorithm used is `scipy.optimize.least_squares`_ with
    fitting method 'trf' to allow for bounds.

    .. todo::
        For multicomponent lines, set the first normalization to be the
        normalization for the sum of all components, then force the
        normalization of the subcomponents to be ordered from highest to
        lowest and bounded from 0 to 1.

    Args:
        x (array-like):
            1D vector with the independent variable
        y (array-like):
            1D vector with the dependent variable
        profile_list (:obj:`list`, :class:`~mangadap.util.lineprofiles.NCompLineProfile`):
            The profile(s) to fit to the dependent variable.
        base_order (int): (**Optional**) The order of the Legendre
            polynomial to include in the model for the baseline trend in
            y below the fitted line profile(s).
        error (array-like, optional):
            1D vector with the error in the dependent variables. If
            not provided, no error weighting is performed during the
            fitting process, and the covariance *will not* be
            constructed.
        mask (array-like, optional):
            1D boolean array used to ignore values in `y` during the
            fit.
        par (array-like, optional):
            1D vector with the initial guess for model parameters.
            The number of parameters much match the expectation based
            on the provided list of profiles and the order of the
            baseline polynomial. If not provided, the parameters are
            initialized to 0.
        fixed (array-like, optional):
            1D vector with flags used to fix parameters during the
            fit. The number of parameters much match the expectation
            based on the provided list of profiles and the order of
            the baseline polynomial. If not provided, all parameters
            are freely fit.
        bounds (:obj:`tuple`): (**Optional**) 2-tuple with two
            array-like elements giving the upper and lower bound for
            each parameter. The length of each array element must
            match the number of parameters. For an unbounded problem,
            set ``bounds=None``, or use numpy.inf with an appropriate
            sign to disable bounds on all or some variables.
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
        x (`numpy.ndarray`_): Independent variable of length :math:`M`.
        y (`numpy.ndarray`_): Dependent variable to be fit of length
            :math:`M`.
        err (`numpy.ndarray`_): Error, :math:`\sigma`, in the dependent
            variable of length :math:`M`.
        mask (`numpy.ndarray`_): Flag to fit dependent variables of length
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
        par (`numpy.ndarray`_): Full list of model parameters, including
            those parameters that have been fixed.
        fixed (`numpy.ndarray`_): Flags to fit (False) or fix (True) a
            given parameter.
        bounds (:obj:`tuple`): 2-tuple with two array-like
            elements giving the upper and lower bound for each
            parameter.  For an unbounded problem, this is set to
            ``bounds=(-numpy.inf,numpy.inf)``.
        result (`scipy.optimize.OptimizeResult`_): Object with the
            results from `scipy.optimize.least_squares`_.  The
            best-fitting parameters, :math:`\mathbf{\theta}`, is
            returned as ``result.x``.
        cov (`numpy.ndarray`_): The formal covariance matrix for the fit.
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
        TypeError:
            Raised if the provided profile objects are not instances
            of :class:`~mangadap.util.lineprofiles.NCompLineProfile`.
        ValueError:
            Raised if any of the provided parameter arrays (par,
            fixed) are not one-dimensional or the number of
            parameters is not as expected based on the number of
            profile and baseline parameters.

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
            if not isinstance(_profile_list[i], lineprofiles.NCompLineProfile):
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
        self.par[numpy.invert(self.fixed)] = par[numpy.invert(self.fixed)] \
                                                    if len(par) == self.npar else par
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

        # If the input is a masked array, combine the two masks; the
        # first line below works for both normal numpy.ndarrays and
        # MaskedArray
        _mask = numpy.ma.getmaskarray(y)
        if mask is not None:
            _mask |= mask

        # Set the internal vectors to be numpy.ndarrays that only include the unmasked values
        indx = numpy.invert(_mask)
        self.x = numpy.array(x[indx])
        self.y = numpy.array(y[indx])
        self.err = None if error is None else numpy.array(error[indx])

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
            self.result = optimize.OptimizeResult(success=False)
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


class ElricPar(KeywordParSet):
    """
    Elric emission-line fitting parameters.

    The defined parameters are:

    .. include:: ../tables/elricpar.rst
    """
    def __init__(self, emission_lines=None, base_order=None, window_buffer=None,
                 guess_redshift=None, guess_dispersion=None, minimum_snr=None, pixelmask=None,
                 stellar_continuum=None):

        arr_in_fl = [ numpy.ndarray, list, int, float ]
        in_fl = [ int, float ]

        pars =     [ 'emission_lines', 'base_order', 'window_buffer', 'guess_redshift',
                     'guess_dispersion', 'minimum_snr', 'pixelmask', 'stellar_continuum' ]
        values =   [ emission_lines, base_order, window_buffer, guess_redshift, guess_dispersion,
                     minimum_snr, pixelmask, stellar_continuum ]
        defaults = [ None, -1,   25.0, None, None, 0.0, None, None ]
        dtypes =   [ EmissionLineDB, int, in_fl, arr_in_fl, arr_in_fl, in_fl, SpectralPixelMask,
                     StellarContinuumModel ]

        super(ElricPar, self).__init__(pars, values=values, defaults=defaults, dtypes=dtypes)


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
        db_indx (`numpy.ndarray`_): (**Optional**) An array with the
            database index (0-based) of each line to be fit.
        line_index (`numpy.ndarray`_): (**Optional**) An array with the
            index numbers, read from the database, of each line to be
            fit.
        restwave (`numpy.ndarray`_): (**Optional**) An array with the rest
            wavelengths of each line to be fit.
        profile_set (`numpy.ndarray`_): (**Optional**) An array with the
            profile objects that define the functional form of each
            line.
        fixed_par (`numpy.ndarray`_): (**Optional**) An array with the
            fixed parameters for *all* parameters in the model to be
            fit.
        bounds (`numpy.ndarray`_): (**Optional**) A two-column array with
            the lower (first column) and upper (second column) bounds
            for the fit parameters.
        log_bounds (`numpy.ndarray`_): (**Optional**) The range of the
            boundary should be considered as logarithmic when testing if
            a parameter is near its boundary.
        output_model (bool): (**Optional**) Include the best-fitting
            model in the composite emission-line model for each
            spectrum.  This is only flagged as true if ALL the emission
            lines in the fitting window are to be included according to
            the emission-line database.
        tied_pairs (`numpy.ndarray`_): (**Optional**) The series of
            tied profiles (:class:`TiedLineProfile` objects) that are
            used to tie parameters of the model.
        tied_funcs (`numpy.ndarray`_): (**Optional**) The member functions
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
        if not isinstance(profile, lineprofiles.NCompLineProfile):
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
        self.redshift = None
        self.dispersion = None

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


    # TODO: Convert this to a DataTable
    @staticmethod
    def _per_fitting_window_dtype(nwin, max_npar, mask_dtype):
        r"""
        Construct the record array data type for the output fits
        extension.
        """

        return [ ('BINID',numpy.int),
                 ('BINID_INDEX',numpy.int),
                 ('MASK', mask_dtype, (nwin,)),
                 ('NPIXFIT',numpy.int,(nwin,)),
                 ('PAR',numpy.float,(nwin,max_npar)),
                 ('ERR',numpy.float,(nwin,max_npar)),
                 ('LOBND',numpy.float,(nwin,max_npar)),
                 ('UPBND',numpy.float,(nwin,max_npar)),
                 ('FIXED',numpy.bool,(nwin,max_npar)),
#                 ('TIED',numpy.int,(nwin,max_npar)),
                 ('IGNORE',numpy.bool,(nwin,max_npar)),
                 ('CHI2',numpy.float,(nwin,)),
                 ('RCHI2',numpy.float,(nwin,)),
                 ('RMS',numpy.float,(nwin,)),
                 ('RESID',numpy.float,(nwin,7)),
                 ('FRAC_RMS',numpy.float,(nwin,)),
                 ('FRAC_RESID',numpy.float,(nwin,7))
               ]


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
            profile_set = numpy.array([ lineprofiles.NCompLineProfile(
                                        self.emission_lines['ncomp'][primary_line][i],
                                        par=self.emission_lines['par'][primary_line][i],
                                        profile=eval('lineprofiles.'
                                            +self.emission_lines['profile'][primary_line][i])) ])
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
                                                  lineprofiles.NCompLineProfile(e['ncomp'],
                                                                                par=e['par'],
                                                                   profile=eval('lineprofiles.'
                                                                                +e['profile'])),
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


#    @staticmethod
#    def _check_db(emission_lines):
#
#        # Check the input type
#        if not isinstance(emission_lines, EmissionLineDB):
#            raise TypeError('Input database must have type EmissionLineDB.')
#
#        # Check the database itself
#        neml = emission_lines.nsets
#        for i in range(neml):
#            profile = eval(emission_lines['profile'][i])
#            npar = len(profile.param_names)
#            if emission_lines['par'][i].size != npar*emission_lines['ncomp'][i]:
#                raise ValueError('Provided {0} parameters, but expected {1}.'.format(
#                                  emission_lines['par'][i].size, npar*emission_lines['ncomp'][i]))
#            if emission_lines['fix'][i].size != npar*emission_lines['ncomp'][i]:
#                raise ValueError('Provided {0} fix flags, but expected {1}.'.format(
#                                  emission_lines['fix'][i].size, npar*emission_lines['ncomp'][i]))
#            if numpy.any([f not in [0, 1] for f in emission_lines['fix'][i] ]):
#                warnings.warn('Fix values should only be 0 or 1; non-zero values interpreted as 1.')
#            if emission_lines['lobnd'][i].size != npar*emission_lines['ncomp'][i]:
#                raise ValueError('Provided {0} lower bounds, but expected {1}.'.format(
#                                  emission_lines['lobnd'][i].size,
#                                  npar*emission_lines['ncomp'][i]))
#            if emission_lines['hibnd'][i].size != npar*emission_lines['ncomp'][i]:
#                raise ValueError('Provided {0} upper bounds, but expected {1}.'.format(
#                                  emission_lines['hibnd'][i].size,
#                                  npar*emission_lines['ncomp'][i]))
#            if emission_lines['log_bnd'][i].size != npar*emission_lines['ncomp'][i]:
#                raise ValueError('Provided {0} log boundaries designations, but expected '
#                                 '{1}.'.format(emission_lines['log_bnd'][i].size,
#                                               npar*emission_lines['ncomp'][i]))


    @staticmethod
    def _is_near_bound(par, lbnd, ubnd, logbounded, rtol=1e-2, atol=1e-4):
        """
        Determine of any of the parameters within start and start+npar
        are "near" an imposed boundary.
        """
#        print('par: ', par)
#        print('lbnd: ', lbnd)
        lbounded = numpy.invert(numpy.isinf(lbnd))
#        print('L bounded: ', lbounded)

#        print('ubnd: ', ubnd)
        ubounded = numpy.invert(numpy.isinf(ubnd))
#        print('U bounded: ', ubounded)

        if numpy.sum(lbounded) == 0 and numpy.sum(ubounded) == 0:
#            print('No bounds')
            return False

        if numpy.any(lbounded & numpy.invert(ubounded) & (par - lbnd < atol)):
#            print('Near lower bound')
            return True

        if numpy.any(numpy.invert(lbounded) & ubounded & (ubnd - par < atol)):
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
        indx = ulbounded & logbounded
        Dp[indx] = numpy.log10(ubnd[indx]) - numpy.log10(lbnd[indx])
        indx = ulbounded & numpy.invert(logbounded)
        Dp[indx] = ubnd[indx] - lbnd[indx]
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

        if numpy.any(ulbounded & numpy.invert(logbounded) & ((par-lbnd < tol) | (ubnd-par < tol))):
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

            - (DONE) Check the success failure
            - (DONE) Calculate the chi2, rchi2, rms, fractional_rms,
              residuals, fractional_residuals
            - (DONE) Compare the fit parameters to the bounds
            - Check the chi-square and residuals?
            - Check the velocity offset wrt the input
            - For multiple lines, check the order of the lines matches
              the guess parameters
            - For multiple component lines, check the ordering of the
              subcomponents

        """
        # Instantiate the returned flag that the parameter are near the
        # imposed boundary; only meaningful if the fit did not fail
        near_bound = False

        # Assign the fitting window
        model_eml_par['FIT_INDEX'][i,self.fitting_window[j].db_indx] = j
#        print(model_eml_par['FIT_INDEX'][i,:])

        # If the parameters are bounded, save the bounds
        npar = self.bestfit[i,j].npar
        if self.bestfit[i,j].bounded:
            model_fit_par['LOBND'][i,j,:npar] = self.bestfit[i,j].bounds[0]
            model_fit_par['UPBND'][i,j,:npar] = self.bestfit[i,j].bounds[1]

        # Save which parameters were fixed or should be ignored
        model_fit_par['IGNORE'][i,j,npar:] = True
        
#        include_model = True
        if not self.bestfit[i,j].result.success:
            model_fit_par['MASK'][i,j] = self.bitmask.turn_on(model_fit_par['MASK'][i,j],
                                                              'FIT_FAILED')
            if not self.quiet:
                log_output(self.loggers, 1, logging.ERROR,
                           'FIT FAILED: Spec: {0}/{1} ({2} : {3})'.format(i+1,self.nspec,
                                                            self.fitting_window[j].line_index,
                                    self.emission_lines['name'][self.fitting_window[j].db_indx]))
            
            return near_bound

        # Get the full list of parameters (including fixed and tied values)
        self.bestfit[i,j]._assign_par(self.bestfit[i,j].result.x)
        par = self.bestfit[i,j].par
        model_fit_par['PAR'][i,j,:npar] = par[:]
        model_fit_par['FIXED'][i,j,:npar] = self.bestfit[i,j].fixed
        model_fit_par['FIXED'][i,j,npar:] = True
        model_fit_par['ERR'][i,j,:npar] = numpy.zeros(npar, dtype=numpy.float)
        if self.error is not None:
            tmperr = numpy.diag(self.bestfit[i,j].cov)
            if not numpy.any(tmperr < 0):
                model_fit_par['ERR'][i,j,numpy.invert(model_fit_par['FIXED'][i,j])] \
                        = numpy.sqrt(numpy.diag(self.bestfit[i,j].cov))
            else:
                model_fit_par['MASK'][i,j] = self.bitmask.turn_on(model_fit_par['MASK'][i,j],
                                                                  'UNDEFINED_COVAR')


#        print('PAR: ', model_fit_par['PAR'][i,j,:npar])
#        print('ERR: ', model_fit_par['ERR'][i,j,:npar])

        # Get the per-line parameters
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

        model_fit_par['RESID'][i,j] = sample_growth(numpy.ma.absolute(resid),
                                                    [0.0, 0.25, 0.50, 0.75, 0.90, 0.99, 1.0])
#        model_fit_par['RESID'][i,j] = residual_growth(resid, [0.25, 0.50, 0.75, 0.90, 0.99])
#        print('RESID: ', model_fit_par['RESID'][i,j])

        fit = self.bestfit[i,j].sample(self.bestfit[i,j].x, par)
        indx = numpy.absolute(fit) > 1e-4
        if numpy.sum(indx) > 0:
            frac_resid = resid[indx]/fit[indx]
            model_fit_par['FRAC_RMS'][i,j] = numpy.sqrt(numpy.mean(numpy.square(frac_resid)))
#            print('FRMS:', model_fit_par['FRAC_RMS'][i,j])
            if numpy.sum(indx) > 1:
                model_fit_par['FRAC_RESID'][i,j] = sample_growth(numpy.ma.absolute(frac_resid),
                                                        [0.0, 0.25, 0.50, 0.75, 0.90, 0.99, 1.0])
#                model_fit_par['FRAC_RESID'][i,j] = residual_growth(frac_resid,
#                                                               [0.25, 0.50, 0.75, 0.90, 0.99])
#                print('FRAC_RESID: ', model_fit_par['FRAC_RESID'][i,j])

        return near_bound


#    def _instrumental_dispersion_correction(self, model_eml_par):
#        """
#        Determine the instrumental velocity dispersion at the centroids
#        of the fitted emission lines.
#        """
#        if self.sres is None:
#            return
#
#        interpolator = interpolate.interp1d(self.wave, self.sres, fill_value='extrapolate',
#                                            assume_sorted=True)
#        restwave = numpy.array([self.emission_lines['restwave']]*self.nspec)
#        cnst = constants()
#
#        # Get the instrumental dispersion at each (valid) line center
#        indx = ~self.bitmask.flagged(model_eml_par['MASK'],
#                                     flag=['INSUFFICIENT_DATA', 'FIT_FAILED', 'NEAR_BOUND',
#                                           'UNDEFINED_COVAR' ])
#        model_eml_par['SINST'][indx] = astropy.constants.c.to('km/s') \
#                                            / interpolator((model_eml_par['KIN'][indx,0] 
#                                                        / astropy.constants.c.to('km/s').value+1.0)
#                                                                *restwave[indx])/cnst.sig2fwhm

#        z = model_eml_par['KIN'][indx,0] / astropy.constants.c.to('km/s').value
#        pyplot.scatter(restwave[indx]*(1.0+z), model_eml_par['SINST'][indx], marker='.', s=30,
#                       color='k')
#        pyplot.show()
#
#        pyplot.scatter(model_eml_par['SINST'][indx], model_eml_par['KIN'][indx,1], marker='.',
#                       s=30, color='k')
#        pyplot.show()


# NEVER TESTED
#    def _correct_velocity_dispersion(self, model_eml_par):
#        """
#        Use the **previously calculated** instrumental velocity
#        dispersion to correct the measured velocity dispersions of the
#        lines to the astrophysical velocity dispersion.
#        """
#
#        # Correct the second moments for the instrumental dispersion
#        defined = numpy.zeros(indx.shape, dtype=numpy.bool)
#        defined[indx] = model_eml_par['SINST'][indx] < model_eml_par['KIN'][indx,1]
#        nonzero = numpy.zeros(indx.shape, dtype=numpy.bool)
#        nonzero[indx] = numpy.absolute(model_eml_par['SINST'][indx] 
#                                            - model_eml_par['KIN'][indx,1]) > 0
#
##        orig = model_eml_par['KIN'].copy()
#
#        model_eml_par['KINERR'][nonzero,1] \
#                    = 2.0*model_eml_par['KIN'][nonzero,1]*model_eml_par['KINERR'][nonzero,1]
#        model_eml_par['KIN'][nonzero,1] = numpy.square(model_eml_par['KIN'][nonzero,1]) \
#                                                - numpy.square(model_eml_par['SINST'][nonzero])
#        model_eml_par['KIN'][nonzero,1] = model_eml_par['KIN'][nonzero,1] \
#                            / numpy.sqrt(numpy.absolute(model_eml_par['KIN'][nonzero,1]))
#        model_eml_par['KINERR'][nonzero,1] /= numpy.absolute(2.*model_eml_par['KIN'][nonzero,1])
#
#        # Flag undefined values
#        model_eml_par['MASK'][indx & ~defined] \
#                = self.bitmask.turn_on(model_eml_par['MASK'][indx & ~defined], 'UNDEFINED_SIGMA')
#
##        flg = self.bitmask.flagged(model_eml_par['MASK'], flag='UNDEFINED_SIGMA')
##        pyplot.scatter(orig[nonzero & flg,1], model_eml_par['KIN'][nonzero & flg,1], marker='.',
##                       s=30, color='0.8')
##        pyplot.scatter(orig[nonzero & ~flg,1], model_eml_par['KIN'][nonzero & ~flg,1], marker='.',
##                       s=30, color='k')
##        pyplot.show()

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

        # Continuum accounts for underlying absorption
        if par['stellar_continuum'] is not None \
                and not isinstance(par['stellar_continuum'], StellarContinuumModel):
            raise TypeError('Must provide a valid StellarContinuumModel object.')
        continuum = None if par['stellar_continuum'] is None \
                        else par['stellar_continuum'].fill_to_match(binned_spectra['BINID'].data,
                                                            missing=binned_spectra.missing_bins)
        if continuum is not None:
            if continuum.shape != binned_spectra['FLUX'].data.shape:
                raise ValueError('Provided continuum does not match shape of the binned spectra.')
            if not isinstance(continuum, numpy.ma.MaskedArray):
                continuum = numpy.ma.MaskedArray(continuum)

        # Get the data arrays to fit
        good_snr = binned_spectra.above_snr_limit(par['minimum_snr'])
        wave, flux, ivar, sres = EmissionLineFit.get_spectra_to_fit(binned_spectra,
                                                                    pixelmask=par['pixelmask'],
                                                                    select=good_snr)

        # Return the fitted data
        model_wave, model_flux, model_base, model_mask, model_fit_par, model_eml_par \
                = self.fit(binned_spectra['WAVE'].data, flux, par['emission_lines'],
                           ivar=ivar, sres=sres, continuum=continuum[good_snr,:],
                           base_order=par['base_order'], window_buffer=par['window_buffer'],
                           guess_redshift=par['guess_redshift'][good_snr],
                           guess_dispersion=par['guess_dispersion'][good_snr], loggers=loggers,
                           quiet=quiet)

        # Save the the bin ID numbers indices based on the spectra
        # selected to be fit
        model_fit_par['BINID'] = binned_spectra['BINS'].data['BINID'][good_snr]
        model_fit_par['BINID_INDEX'] = numpy.arange(binned_spectra.nbins)[good_snr]

        model_eml_par['BINID'] = binned_spectra['BINS'].data['BINID'][good_snr]
        model_eml_par['BINID_INDEX'] = numpy.arange(binned_spectra.nbins)[good_snr]

        # Add the equivalent width data
        EmissionLineFit.measure_equivalent_width(binned_spectra['WAVE'].data, flux,
                                                 par['emission_lines'], model_eml_par,
                                                 redshift=par['guess_redshift'][good_snr],
                                                 bitmask=self.bitmask)

        # Only return model and model parameters for the *fitted*
        # spectra
        return model_flux, model_base, model_mask, model_fit_par, model_eml_par, None


    def fit(self, wave, flux, emission_lines, ivar=None, mask=None, sres=None,
            continuum=None, base_order=-1, window_buffer=25, guess_redshift=None,
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
        # Initialize the reporting
        if loggers is not None:
            self.loggers = loggers
        self.quiet = quiet

        # Check the input emission-line database
        EmissionLineFit.check_emission_line_database(emission_lines)
        self.emission_lines = emission_lines
        self.neml = self.emission_lines.neml

        # Prepare the spectra for fitting
        self.wave = wave
        self.flux, self.error, self.sres, self.continuum, self.redshift, self.dispersion \
                    = EmissionLineFit.check_and_prep_input(self.wave, flux, ivar=ivar, mask=mask,
                                                           sres=sres, continuum=continuum,
                                                           redshift=guess_redshift,
                                                           dispersion=guess_dispersion,
                                                        default_dispersion=self.default_dispersion)
        # Keep the shape
        self.nspec, self.nwave = self.flux.shape
        # Subtract the continuum
        self.fluxnc, self.no_continuum = EmissionLineFit.subtract_continuum(self.flux,
                                                                            self.continuum)
        # No continuum for any spaxel!
        if self.no_continuum is None:
            self.no_continuum = numpy.ones(self.flux.shape, dtype=numpy.bool)

#        pyplot.step(self.wave, self.flux[0,:], where='mid', color='k', lw=0.5, linestyle='-')
#        pyplot.step(self.wave, self.fluxnc[0,:], where='mid', color='g', lw=0.5, linestyle='-')
#        pyplot.plot(self.wave, self.no_continuum[0,:], color='r', lw=0.5, linestyle='-')
#        pyplot.show()
#        exit()

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
                                          Elric._per_fitting_window_dtype(self.nwindows, max_npar,
                                                                self.bitmask.minimum_dtype()))
        model_fit_par['BINID'] = numpy.arange(self.nspec)
        model_fit_par['BINID_INDEX'] = numpy.arange(self.nspec)

        # Emission-line parameters
        model_eml_par = self.emission_line_datatable(self.neml, 2, self.bitmask.minimum_dtype(),
                                                     shape=self.nspec)
#        model_eml_par = init_record_array(self.nspec,
#                                          self._per_emission_line_dtype(self.neml, 2,
#                                                                self.bitmask.minimum_dtype()))
        model_eml_par['CONTMPLY'] = numpy.ones(model_eml_par['CONTMPLY'].shape, dtype=float)
        model_eml_par['BINID'] = numpy.arange(self.nspec)
        model_eml_par['BINID_INDEX'] = numpy.arange(self.nspec)

        t = time.perf_counter()
        # Do the fit
        velocity = self._velocity_vectors(self.wave, self.fitting_window)
        self.bestfit = numpy.empty((self.nspec,self.nwindows), dtype=object)
        for i in range(self.nspec):
#        for i in range(5):

#            for j in range(self.nwindows):
#                print('Window: {0}/{1}'.format(j+1, self.nwindows))
#                for p in self.fitting_window[j].profile_set:
#                    print(p.par)
#            if not self.quiet:
#                log_output(self.loggers, 1, logging.INFO, 'Fit: {0}/{1}'.format(i+1,self.nspec))
            print('Fit: {0}/{1}'.format(i+1,self.nspec), end='\r')

            # Get the fitting mask for each emission-line in each window
            fitting_mask = self._fit_masks(self.wave, self.fitting_window, self.redshift[i],
                                           self.window_buffer)

            # Set the mask for any pixels that are NEVER fit for this
            # spectrum
            indx = numpy.invert(numpy.any(fitting_mask, axis=0))
            model_mask[i,indx] = self.bitmask.turn_on(model_mask[i,indx], 'OUTSIDE_RANGE')

            model_jump = numpy.array([
                        numpy.sum(fm & self.no_continuum[i,:]) not in [0, numpy.sum(fm)] \
                        for fm in numpy.take(fitting_mask, numpy.arange(self.nwindows), axis=0) ])
            spec_to_fit = numpy.ma.array([ self.flux[i,:] if mj else self.fluxnc[i,:] \
                                            for mj in model_jump ])

            cz = self.redshift[i]*astropy.constants.c.to('km/s').value

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
                # TODO: Use self.dispersion to set the initial
                # dispersion guesses?  Currently uses guesses from
                # emission-line database (see
                # _parse_emission_line_models)
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
#                print(_guess_par)
#                print(_fixed_par)
#                print(_bounds)

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

                # Add the parameters for the baseline
                if base_order > -1:
                    _guess_par = numpy.append(_guess_par, numpy.zeros(base_order+1))
#                    _guess_par[-base_order-1] = numpy.median(spec_to_fit[j,fitting_mask[j,:]])
#                    print(0.0, _guess_par[-base_order-1])
#                    pyplot.scatter(velocity[j,fitting_mask[j,:]], spec_to_fit[j,fitting_mask[j,:]],
#                                   marker='.', s=100, color='k')
#                    y = velocity.copy()[j,fitting_mask[j,:]]
#                    y[:] = numpy.median(spec_to_fit[j,fitting_mask[j,:]])
#                    pyplot.plot(velocity[j,fitting_mask[j,:]], y, color='r')
#                    pyplot.show()
#                    exit()
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
                if not numpy.all(_bounds[:,0]<_bounds[:,1]):

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
                                                   base_order=base_order,
                                                   mask=numpy.invert(fitting_mask[j,:]),
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

            # Determine the instrumental dispersion correction at the line centers
            indx = numpy.invert(self.bitmask.flagged(model_eml_par['MASK'][i,:],
                                                     flag=['INSUFFICIENT_DATA', 'FIT_FAILED',
                                                           'NEAR_BOUND', 'UNDEFINED_COVAR' ]))
            model_eml_par['SIGMACORR'][i,indx] \
                    = EmissionLineFit.instrumental_dispersion(self.wave, self.sres[i,:],
                                                              emission_lines['restwave'][indx],
                                                              model_eml_par['KIN'][i,indx,0])

        # With these definitions, the instrumental sigma and the sigma correction are
        # identical, and the "template" sigma are all zero
        model_eml_par['SIGMAINST'] = model_eml_par['SIGMACORR'].copy()

        print('Fit: {0}/{0}'.format(i+1,self.nspec))
        # Remove the lines from the full model to just provide the
        # baseline model
        model_base -= model_flux

#        pyplot.plot(wave, model_flux[0,:], linestyle='-', color='g')
#        pyplot.plot(wave, model_base[0,:], linestyle='-', color='r')
#        pyplot.show()

        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, 'Fits completed in {0:.4e} min.'.format(
                       (time.perf_counter() - t)/60))

        return self.wave, model_flux, model_base, model_mask, model_fit_par, model_eml_par


