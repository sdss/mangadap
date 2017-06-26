# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Implements an emission-line fitting class that largely wraps pPXF.

*License*:
    Copyright (c) 2017, SDSS-IV/MaNGA Pipeline Group
        Licensed under BSD 3-clause license - see LICENSE.rst

*Source location*:
    $MANGADAP_DIR/python/mangadap/proc/sasuke.py

*Imports and python version compliance*:
    ::

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
        from scipy import interpolate, fftpack
        import astropy.constants

        from ..par.parset import ParSet
        from ..par.emissionlinedb import EmissionLineDB
        from ..util.fileio import init_record_array
        from ..util.instrument import spectrum_velocity_scale, resample_vector
        from ..util.log import log_output
        from ..util.pixelmask import SpectralPixelMask
        from .spatiallybinnedspectra import SpatiallyBinnedSpectra
        from .stellarcontinuummodel import StellarContinuumModel
        from .spectralfitting import EmissionLineFit
        from .util import residual_growth

*Class usage examples*:
        Add examples

*Revision history*:
    | **24 May 2017**: Original implementation started by K. Westfall (KBW)
    | **23 Jun 2017**: (KBW) Documentation; fix error in
        :func:`Sasuke._save_results`

.. _numpy.ma.MaskedArray: https://docs.scipy.org/doc/numpy-1.12.0/reference/maskedarray.baseclass.html
.. _numpy.recarray: https://docs.scipy.org/doc/numpy/reference/generated/numpy.recarray.html

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
import itertools

import numpy
from scipy import interpolate, fftpack

import astropy.constants

from ..par.parset import ParSet
from ..par.emissionlinedb import EmissionLineDB
from ..contrib.ppxf import ppxf
from ..util.fileio import init_record_array
from ..util.instrument import spectrum_velocity_scale, spectral_coordinate_step
from ..util.instrument import SpectralResolution
from ..util.log import log_output
from ..util.pixelmask import SpectralPixelMask
from ..util.constants import DAPConstants
from ..util import lineprofiles
from .spatiallybinnedspectra import SpatiallyBinnedSpectra
from .stellarcontinuummodel import StellarContinuumModel
from .spectralfitting import EmissionLineFit
from .bandpassfilter import emission_line_equivalent_width
from .util import residual_growth
from .ppxffit import PPXFFit, PPXFFitResult

# For debugging
from matplotlib import pyplot

class SasukePar(ParSet):
    """
    A class specific to the DAP's use of Sasuke.

    This is the object that gets passed to
    :func:`Sasuke.fit_SpatiallyBinnedSpectra`.  In the DAP, it is
    instantiated by
    :func:`mangadap.proc.emissionlinemodel.available_emission_line_modeling_methods`
    and some of its components are filled by
    :func:`mangdap.proc.emissionlinemodel.EmissionLineModel._fill_method_par`.

    .. todo::
        It should be possible to allow for a replacement set of
        templates but I (KBW) need remember how to do this.

    When instantiated, the :class:`mangadap.par.parset.ParSet` objects
    test that the input objects match the provided dtypes.  See
    documentation for :class:`mangadap.par.parset.ParSet` for the list
    of attributes and exceptions raised.

    Args:
        stellar_continuum
            (:class:`mangadap.proc.stellarcontinuummodel.StellarContinuumModel`):
            The result of the previous fit to the stellar continuum.
        emission_lines
            (:class:`mangadap.par.emissionlinedb.EmissionLineDB`):
            Emission-line database with the details of the lines to be
            fit.
        guess_redshift (array-like): Single or per-spectrum redshift to
            use as the initial velocity guess.
        guess_dispersion (array-like): Single or per-spectrum velocity
            dispersion to use as the initial guess.
        minimum_snr (float): Minimum S/N of spectrum to fit.
        pixelmask (:class:`mangadap.util.pixelmask.SpectralPixelMask`):
            Mask to apply to all spectra being fit.
        reject_boxcar (int): Size of the boxcar to use when rejecting
            fit outliers.
        bias (float): pPXF bias parameter.  (Irrelevant because gas is
            currently always fit with moments=2.)
        moments (int): pPXF moments parameter.  (Irrelevant because gas is
            currently always fit with moments=2.)
        degree (int): pPXF degree parameter setting the degree of the
            additive polynomials to use.
        mdegree (int): pPXF mdegree parameter setting the degree of the
            multiplicative polynomials to use.
    """
    def __init__(self, stellar_continuum, emission_lines, guess_redshift=None,
                 guess_dispersion=None, minimum_snr=None, pixelmask=None,
                 reject_boxcar=None, bias=None, moments=None, degree=None, mdegree=None):

        arr_like = [ numpy.ndarray, list ]
        arr_in_fl = [ numpy.ndarray, list, int, float ]
        in_fl = [ int, float ]

        pars =     [ 'stellar_continuum', 'emission_lines', 'guess_redshift', 'guess_dispersion',
                     'minimum_snr', 'pixelmask', 'reject_boxcar', 'bias', 'moments', 'degree',
                     'mdegree' ]
        values =   [ stellar_continuum, emission_lines, guess_redshift, guess_dispersion,
                     minimum_snr, pixelmask, reject_boxcar, bias, moments, degree, mdegree ]
        defaults = [ None, None, None, None, 0.0, None, None, None, 2, -1, 8 ]
        dtypes =   [ StellarContinuumModel, EmissionLineDB, arr_in_fl, arr_in_fl, in_fl,
                     SpectralPixelMask, int, in_fl, int, int, int ]

        ParSet.__init__(self, pars, values=values, defaults=defaults, dtypes=dtypes)


class Sasuke(EmissionLineFit):
    r"""
    Use ninja skills and pPXF to fit emission lines.

    https://en.wikipedia.org/wiki/Sasuke_Uchiha

    Effectively, **nothing** happens during the instantiation of this
    object.  A typical usage of the class to fit a set of emission lines
    would be::

        # Read the emission-line database
        emldb = EmissionLineDB('ELPMILES')
        # Instantiate the emission-line fitter
        el_fitter = Sasuke(EmissionLineModelBitMask())
        # Fit the spectra
        model_wave, model_flux, model_eml_flux, model_mask, model_fit_par, \
            model_eml_par = el_fitter.fit(...)

    See :func:`fit` for the arguments to the main fitting function.

    The class inherits attributes from
    :class:`mangadap.proc.spectralfitting.EmissionLineFit` (`fit_type`,
    `bitmask`, `par`, `fit_method`).  Other attributes are defined upon
    instantiation and set to None.  This isn't necessary (attributes can
    be defined elsewhere in the class methods), but it provides a
    collation of the class attributes for reference.

    Args:
        bitmask (:class:`BitMask`): BitMask object use to flag fit
            results.  This *must* be provided and should typically be an
            instantiation of :class:`EmissionLineModelBitMask`; however,
            it can be any object with :class:`BitMask` as its base
            class.  The flags set within the main fitting function
            (:func:`Sasuke.fit`) are: DIDNOTUSE, INSUFFICIENT_DATA,
            FIT_FAILED, NEAR_BOUND, NO_FIT,
            :attr:`mangadap.proc.ppxffit.PPXFFit.rej_flag`
            (PPXF_REJECT), MIN_SIGMA, BAD_SIGMA, and MAXITER.  Also the
            DAP-specific calling function
            (:func:`Sasuke.fit_SpatiallyBinnedSpectra`) will also assign
            bits NON_POSITIVE_CONTINUUM during the equivalent width
            measurements (see
            :func:`mangadap.spectralfitting.EmissionLineFit.measure_equivalent_width`).
        loggers (list): (**Optional**) List of `logging.Logger`_ objects
            to log progress; ignored if quiet=True.  Logging is done
            using :func:`mangadap.util.log.log_output`.  Default is no
            logging.  This can be reset in some methods.
        quiet (bool): (**Optional**) Suppress all terminal and logging
            output.  Default is False.

    Attributes:
        loggers (list): List of `logging.Logger`_ objects to log
            progress; ignored if quiet=True.  Logging is done using
            :func:`mangadap.util.log.log_output`.
        quiet (bool): Suppress all terminal and logging output.
        obj_wave (numpy.ndarray): Wavelength vector for object spectra.
            Shape is (:math:`N_{\rm pix}`,).
        obj_flux (`numpy.ma.MaskedArray`_): Object spectra to fit.
            Shape is (:math:`N_{\rm spec},N_{\rm pix}`).
        obj_ferr (`numpy.ma.MaskedArray`_): :math:`1\sigma` errors in
            object spectra.  Shape is (:math:`N_{\rm spec},N_{\rm
            pix}`).
        obj_sres (numpy.ndarray): Spectral resolution array for object
            spectra.  Shape is (:math:`N_{\rm spec},N_{\rm pix}`).
        nobj (int): Number of object spectra (i.e., :math:`N_{\rm
            spec}`).
        npix_obj (int): Number of pixels in each object spectrum (i.e.,
            :math:`N_{\rm pix}`).
        input_obj_mask (numpy.ndarray): A copy of the input mask array
            of the object spectra (boolean array).  Shape is
            (:math:`N_{\rm spec},N_{\rm pix}`).
        obj_to_fit (numpy.ndarray): Flag to fit each object spectrum.
            Instantiating by fully masked spectra in
            :attr:`input_obj_mask`.  Shape is (:math:`N_{\rm spec}`,).
        input_cz (numpy.ndarray): Input redshifts (in km/s) for each
            spectrum.  Shape is (:math:`N_{\rm spec}`,).
        velscale (float): Velocity scale (km/s) of the object spectra.
        tpl_wave (numpy.ndarray): Wavelength vector for template spectra.
            Shape is (:math:`N_{\rm pix,t}`,).
        tpl_flux (numpy.ndarray): Template spectra to use in fit.
            Shape is (:math:`N_{\rm tpl},N_{\rm pix,t}`).
        tpl_sres (:class:`mangadap.util.instrument.SpectralResolution`):
            Spectral resolution of the template spectra.  All templates
            are assumed to have the same spectral resolution.
        tpl_to_use (numpy.ndarray): Set of flags used to select the
            templates used to fit each spectrum.  Shape is
            (:math:`N_{\rm spec},N_{\rm tpl}`).
        nstpl (int): Number of stellar templates.
        ntpl (int): Total number of templates (gas + stars).
        npix_tpl (int): Number of pixels in template spectra (i.e.,
            :math:`N_{\rm pix,t}`).
        tpl_npad (int): Nearest length for FFT, :math:`N_{\rm pad}`
        tpl_rfft (numpy.ndarray): The complex array with the real FFT of
            the template spectra.  Shape is (:math:`N_{\rm tpl}, N_{\rm
            pad}/2 + 1`).
        matched_resolution (bool): The spectral resolution of the
            templates is matched to that of the galaxy data.  WARNING:
            This functionality needs to be checked in relation to the
            gas templates!
        velscale_ratio (int): The **integer** ratio between the velocity
            scale of the pixel in the galaxy data to that of the
            template data.
        emldb (:class:`mangadap.par.emissionlinedb.EmissionLineDB'):
            Emission-line database that is parsed to construct the
            emission-line templates (see
            :class:`EmissionLineTemplates`).
        neml (int): Number of emission lines in the database.
        fit_eml (numpy.ndarray): Boolean array setting which emission
            lines are fit.  Shape is (:math:`N_{\rm eml`,).
        eml_tpli (numpy.ndarray): Integer array with the template that
            includes each emission line.  Shape is (:math:`N_{\rm eml`,).
        eml_compi (numpy.ndarray): Integer array with the kinematic
            component that includes each emission line.  Shape is
            (:math:`N_{\rm eml`,).
        ncomp (int): Total number of kinematic components to fit.
        gas_comp (numpy.ndarray): Boolean array set to True for
            emission-line components.  Shape is (:math:`N_{\rm comp}`).
        tpl_comp (numpy.ndarray): The integer kinematic component
            associated with each template.  Shape is (:math:`N_{\rm
            tpl},). 
        tpl_vgrp (numpy.ndarray): The integer velocity group associated
            with each template.  Shape is (:math:`N_{\rm tpl},). 
        tpl_sgrp (numpy.ndarray): The integer sigma group associated
            with each template.  Shape is (:math:`N_{\rm tpl},). 
        comp_moments (numpy.ndarray): Number of moments for each
            component.  Moments with negative numbers have fixed
            kinematics.  Shape is (:math:`N_{\rm comp},).
        comp_start_kin (numpy.ndarray): Array of lists where each list
            provdes the starting parameters for the kinematics of each
            component.  Shape is (:math:`N_{\rm comp},).
        npar_kin (int): The total number of kinematic parameters, which
            is just the sum of the absolute value of
            :attr:`comp_moments`, :math:`N_{\rm kin}`.
        nfree_kin (int): The total number of *free* kinematic
            parameters.
        tied (list): List of lists setting the tied parameters for the
            fit.  See the TIED parameter in ppxf.  Length is
            :math:`N_{\rm comp}`.
        velocity_limits (numpy.ndarray): The upper and lower velocity
            limits imposed by pPXF.  See
            :func:`mangadap.proc.ppxffit.PPXFFit.losvd_limits`.
        sigma_limits (numpy.ndarray): The upper and lower velocity
            dispersion limits imposed by pPXF.  See
            :func:`mangadap.proc.ppxffit.PPXFFit.losvd_limits`.
        gh_limits (numpy.ndarray): The upper and lower limits on *all*
            higher order Gauss-Hermite moments imposed by pPXF.  See
            :func:`mangadap.proc.ppxffit.PPXFFit.losvd_limits`.
        bias (float): pPXF bias parameter.  (Currently irrelevant
            because gas is currently always fit with moments=2.)
        degree (int): pPXF degree parameter setting the degree of the
            additive polynomials to use, :math:`o_{\rm add}`.
        mdegree (int): pPXF mdegree parameter setting the degree of the
            multiplicative polynomials to use, :math:`o_{\rm mult}`.
        reject_boxcar (int): Size of the boxcar to use when rejecting
            fit outliers.
        spectrum_start (numpy.ndarray): Array with the starting index of
            the pixel in the object spectra to fit (inclusive).  Shape
            is (:math:`N_{\rm spec}`,).
        spectrum_end (numpy.ndarray): Array with the ending index of the
            pixel in the object spectra to fit (exclusive).  Shape is
            (:math:`N_{\rm spec}`,)
        dof (int): Degrees of freedom in the fit.
        base_velocity (numpy.ndarray): The base velocity shift between
            the template and object spectra because of the difference in
            their starting wavelength.  Shape is (:math:`N_{\rm
            spec}`,).

    .. todo::

        - :attr:`velocity_limits`, :attr:`sigma_limits`, and
          :attr:`gh_limits` are **not** passed to ppxf during the fit.
          They are expected to match what's in the pPXF code.  Should
          change the code so that they **are** passed using pPXF's BOUND
          keyword.
    """
    def __init__(self, bitmask, loggers=None, quiet=False):

        EmissionLineFit.__init__(self, 'sasuke', bitmask)
        # Attributes kept by SpectralFitting:
        #   fit_type='emission_line', bitmask=bitmask, par=None
        # Attributes kept by EmissionLineFit:
        #   fit_method='sasuke'

        # Logging and terminal output
        self.loggers = loggers
        self.quiet = quiet

        # Data to fit
        self.obj_wave = None
        self.obj_flux = None
        self.obj_ferr = None
        self.obj_sres = None
        self.nobj = None
        self.npix_obj = None
        self.input_obj_mask = None
        self.obj_to_fit = None
        self.input_cz = None
        self.velscale = None
        self.waverange = None

        # Template data
        self.tpl_wave = None
        self.tpl_flux = None
        self.tpl_sres = None
        self.tpl_to_use = None
        self.nstpl = None
        self.ntpl = None
        self.npix_tpl = None
        self.tpl_npad = None
        self.tpl_rfft = None

        self.matched_resolution = None
        self.velscale_ratio = None

        self.emldb = None
        self.neml = None
        self.fit_eml = None
        self.eml_tpli = None
        self.eml_compi = None

        # Kinematic components and tied parameters
        self.ncomp = None
        self.gas_comp = None
        self.tpl_comp = None
        self.tpl_vgrp = None
        self.tpl_sgrp = None
        self.comp_moments = None
        self.comp_start_kin = None
        self.npar_kin = None
        self.nfree_kin = None
        self.tied = None

        # Fitting parameters
        self.velocity_limits = None
        self.sigma_limits = None
        self.gh_limits = None

        self.bias = None
        self.degree = None
        self.mdegree = None
        self.reject_boxcar = None
#        self.fix_kinematics = False

        self.spectrum_start = None
        self.spectrum_end = None
        self.dof = None
        self.base_velocity = None


    @staticmethod
    def _per_fit_dtype(ntpl, nadd, nmult, nkin, mask_dtype):
        r"""
        Static method (independent of self) that returns the data type
        for the result of each pPXF fit.  The data are as follows:

            - ``BINID``: Spectrum ID number.
            - ``BINID_INDEX``: Index of the spectrum in the list of
              provided spectra.
            - ``MASK``: The bitmask value of the fit.  See
              :func:`_save_results`.
            - ``BEGPIX``: Same as :attr:`spectrum_start`.
            - ``ENDPIX``: Same as :attr:`spectrum_end`.
            - ``NPIXTOT``: Length of spectrum passed to pPXF
              (``ENDPIX-BEGPIX``)
            - ``NPIXFIT``: Number of pixels included in the fit
              (excluding rejected/masked pixels)
            - ``KINCMP``: The kinematic component assigned to each
              template.  Shape is (:math:`N_{\rm tpl}`,)
            - ``VELCMP``: The velocity group assigned to each template.
              Shape is (:math:`N_{\rm tpl}`,)
            - ``SIGCMP``: The sigma group assigned to each template.
              Shape is (:math:`N_{\rm tpl}`,)
            - ``USETPL``: Flag that the template was used during the
              fit.  Shape is (:math:`N_{\rm tpl}`,)
            - ``TPLWGT``: Optimal weight of each template in the fit.
              Shape is (:math:`N_{\rm tpl}`,)
            - ``ADDCOEF``: Additive polynomal coefficients.  Shape is
              (:math:`o_{\rm add}+1`,)
            - ``MULTCOEF``: Additive polynomal coefficients.  Shape is
              (:math:`o_{\rm mult}`,)
            - ``KININP``: Initial guess kinematics.  Shape is
              (:math:`N_{\rm kin}`,)
            - ``KIN``: Best-fitting kinematics.  Shape is (:math:`N_{\rm
              kin}`,)
            - ``KINERR``: Errors in the best-fitting kinematics.  Shape
              is (:math:`N_{\rm kin}`,)
            - ``TIEDKIN``: Index of the kinematic parameter to which
              each parameter is tied.  I.e., TIEDKIN[3] = 1 means that
              parameter 4 is tied to paramter 2.  Shape is
              (:math:`N_{\rm kin}`,)
            - ``CHI2``: Chi-square of the fit
            - ``RCHI2``: Reduced chi-square (recalculated after pPXF).
            - ``ROBUST_RCHI2``: Reduced chi-square (returned by pPXF).
            - ``RMS``: RMS of the fit residuals.
            - ``ABSRESID``: The minimum, 68%, 95%, and 99% growth, and
              maximum absolute residual.  Shape is (5,).
            - ``FRMS``: RMS of the fractional fit residuals, where the
              fractional residuals are (data-model)/model.
            - ``FABSRESID``: The minimum, 68%, 95%, and 99% growth, and
              maximum absolute fractional residual.  Shape is (5,).

        In :func:`fit`, BINID and BINID_INDEX are the same.  They're not
        in :func:`fit_SpatiallyBinnedSpectra`.

        Args:
            ntpl (int): Number of templates, :math:`N_{\rm tpl}`.
            nadd (int): Number of additive polynomial coefficents,
                :math:`o_{\rm add}+1`.
            nmult (int): Number of multiplicative polynomial
                coefficients, :math:`o_{\rm mult}`.
            nkin (int): Number of kinematic parameters, :math:`N_{\rm
                kin}`.
            mask_dtype (dtype): Type for the bitmask variable.  See
                :func:`mangadap.util.bitmask.BitMask.minimum_dtype`.

        Returns:
            list: List of tuples with the name and type of each element
            in the array.

        """
        return [ ('BINID',numpy.int),
                 ('BINID_INDEX',numpy.int),
                 ('MASK', mask_dtype),
                 ('BEGPIX', numpy.int),
                 ('ENDPIX', numpy.int),
                 ('NPIXTOT',numpy.int),
                 ('NPIXFIT',numpy.int),
                 ('KINCMP',numpy.int,(ntpl,)),
                 ('VELCMP',numpy.int,(ntpl,)),
                 ('SIGCMP',numpy.int,(ntpl,)),
                 ('USETPL',numpy.bool,(ntpl,)),
                 ('TPLWGT',numpy.float,(ntpl,)),
                 ('ADDCOEF',numpy.float,(nadd,)) if nadd > 1 else ('ADDCOEF',numpy.float),
                 ('MULTCOEF',numpy.float,(nmult,)) if nmult > 1 else ('MULTCOEF',numpy.float),
                 ('KININP',numpy.float,(nkin,)),
                 ('KIN',numpy.float,(nkin,)),
                 ('KINERR',numpy.float,(nkin,)),
                 ('TIEDKIN',numpy.int,(nkin,)),
                 ('CHI2',numpy.float),
                 ('RCHI2',numpy.float),
                 ('ROBUST_RCHI2',numpy.float),
                 ('RMS',numpy.float),
                 ('ABSRESID',numpy.float,(5,)),
                 ('FRMS',numpy.float),
                 ('FABSRESID',numpy.float,(5,))
               ] 


    def _run_fit_iteration(self, obj_flux, obj_ferr, obj_to_fit, weight_errors=False,
                           component_fits=False, plot=False):
        r"""
        Execute a specific fit iteration.  This is used by
        :func:`_fit_all_spectra` to run through as set of fit iteration.
        All spectra in *obj_flux* are fit with a specific set of ppxf
        parameters.  Additional iteration modes can be added by
        including new arguments to the function.
        
        (FB)  This function fits all the spectra in the cube and returns
        a  PPXFFitResult object with the resulting fit paramters

        Args:
            obj_flux (`numpy.ma.MaskedArray`_): Object spectra to fit
                **for this iteration**.  Different from
                :attr:`obj_flux`.
            obj_ferr (`numpy.ma.MaskedArray`_): Error in spectra to fit
                **for this iteration**.  Different from
                :attr:`obj_flux`.
            obj_to_fit (numpy.ndarray): Boolean array to fit each
                spectrum.
            weight_errors (bool): (**Optional**) Flag to calculate and
                assign errors in the weights; see
                :class:`mangadap.proc.ppxffit.PPXFFitResult`.
            component_fits (bool): (**Optional**) Flag to construct the
                optimal fit for each kinematic component.  This is used
                to construct the best-fitting emission-line model in
                :func:`_emission_line_only_model`; see
                :class:`mangadap.proc.ppxffit.PPXFFitResult`.
            plot (bool): (**Optional**) Passed to
                :class:`mangadap.contrib.ppxf.ppxf` to construct a plot
                showing the result **of each fit**.  Should only used
                when debugging.

        Returns:
            numpy.ndarray : Array with :math:`N_{\rm spec}` instances of
            :class:`mangadap.proc.ppxffit.PPXFFitResult`.

        """
#        linear = fix_kinematics and mdegree < 1
        linear = False

        # Create the object to hold all the fits
        result = numpy.empty(self.nobj, dtype=object)

        # Fit each spectrum
#$$$$$$$$ (FB) this is where the code looks over all the spectra in the cube 
        for i in range(self.nobj):
            print('Running pPXF fit on spectrum: {0}/{1}'.format(i+1,self.nobj), end='\r')
            # Meant to ignore this spectrum
#$$$$$$$$ (FB) NO idea what this is doing, but probalby not important            
            if not obj_to_fit[i]:
                result[i] = None
                continue

            # Get the pixels to fit for this spectrum
            gpm = numpy.where( ~ (obj_flux.mask[i, self.spectrum_start[i]:self.spectrum_end[i]]) )[0]

            # Check if there is sufficient data for the fit
            ntpl_to_use = numpy.sum(self.tpl_to_use[i,:])
            if len(gpm) < self.dof+ntpl_to_use:
#$$$$$$$$ (FB) If insufficient data, then throw a warning 
                if not self.quiet:                                        
                    warnings.warn('Insufficient data points ({0}) to fit spectrum {1}'
                                  '(dof={2}).'.format(len(gpm), i+1, self.dof+ntpl_to_use))
#$$$$$$$$ (FB) If insufficient data, generate the PPXFFitResult object with None and go to next spectrum (continue)
                result[i] = PPXFFitResult(self.degree, self.mdegree, self.spectrum_start[i],
                                          self.spectrum_end[i], self.tpl_to_use[i,:],
                                          None, self.ntpl)
                continue

            # Run ppxf
            if plot:
                pyplot.clf()

#$$$$$$$$ (FB) fit spectrum with ppxf and get ppxf output
            _ppxf=ppxf(self.tpl_flux[self.tpl_to_use[i,:],:].T,
                                 obj_flux.data[i,self.spectrum_start[i]:self.spectrum_end[i]],
                                 obj_ferr.data[i,self.spectrum_start[i]:self.spectrum_end[i]],
                                 self.velscale, self.comp_start_kin[i].tolist(), bias=self.bias, 
                                 component=self.tpl_comp[self.tpl_to_use[i,:]], degree=self.degree,
                                 goodpixels=gpm, linear=linear, mdegree=self.mdegree,
                                 moments=self.comp_moments, plot=plot, quiet=(not plot),
                                 templates_rfft=self.tpl_rfft[self.tpl_to_use[i,:],:].T,
                                 tied=self.tied, velscale_ratio=self.velscale_ratio,
                                 vsyst=-self.base_velocity[i])

#$$$$$$$$ (FB) package PPXF output in DAP-specific results object                                 
            result[i] = PPXFFitResult(self.degree, self.mdegree, self.spectrum_start[i],
                                      self.spectrum_end[i], self.tpl_to_use[i,:],
                                      _ppxf, self.ntpl,
                                      weight_errors=weight_errors, component_fits=component_fits)

            # TODO: check output
#            if result[i].kin[1] < 0:
#                result[i].kin[1] = numpy.absolute(result[i].kin[1]) #self.sigma_limits[0]
#                warnings.warn('pPXF gives negative dispersion! Change -{0:.4f} to {0:.4f}'.format(
#                                    result[i].kin[1]))
                
            if result[i].reached_maxiter() and not self.quiet:
                warnings.warn('pPXF optimizer reached maximum number of iterations for spectrum '
                              '{0}.'.format(i+1))
            if plot:
                pyplot.show()

        print('Running pPXF fit on spectrum: {0}/{0}'.format(self.nobj))
#$$$$$$$$ (FB) this function returns a  PPXFFitResult object         
        return result


    def _fit_all_spectra(self, plot=False): #, plot_file_root=None):
        """
        Fit all spectra contained by the class instance
        (:attr:`obj_flux`) using an optimized set of iterations.  For
        example, see
        :func:`mangadap.proc.ppxffit.PPXFFit._fit_all_spectra()`.
        
        $$$$$$$$ (FB) this is a wrapper to   _run_fit_iteration which allows the 
        user to run a second fit with outlier rejection (DO WE NEED THIS?)
        For testing and simplicity I killed this functionality
        INSTEAD REUSE THIS FUCTION
        TO DO:
            1. run first fit iteration with single stellar template (fixed kin) 
            and all lines tied
            2. depixelise (break bins)
            3. run second iteration with preferred parameters

        Args:
            plot (bool): (**Optional**) Passed to
                :class:`mangadap.contrib.ppxf.ppxf` to construct a plot
                showing the result of **each fit** during **each fit
                iteration**.  Should only used when debugging.
        
        Returns:
            numpy.ndarray : Array with :math:`N_{\rm spec}` instances of
            :class:`mangadap.proc.ppxffit.PPXFFitResult`.
        
        """
#        run_rejection = self.reject_boxcar is not None
#        #---------------------------------------------------------------
#        # Fit the spectra
#        if not self.quiet:
#            log_output(self.loggers, 1, logging.INFO,
#                       'Number of object spectra to fit: {0}/{1}'.format(
#                            numpy.sum(self.obj_to_fit), len(self.obj_to_fit)))
#        result = self._run_fit_iteration(self.obj_flux, self.obj_ferr, self.obj_to_fit,
#                                         weight_errors=(not run_rejection),
#                                         component_fits=(not run_rejection), plot=plot)
#        if not run_rejection:
#            # Only a single fit so return
#            return result
#
#        #---------------------------------------------------------------
#        # Rejection iteration
#
#        # Copy the input as to not overwrite the input masks
#        obj_flux = self.obj_flux.copy()
#        obj_ferr = self.obj_ferr.copy()
#        obj_to_fit = self.obj_to_fit.copy()
#
#        # Save which were not fit successfully
#        obj_to_fit &= numpy.invert(numpy.array([ r is None or r.fit_failed() for r in result ]))
#        if not self.quiet:
#            log_output(self.loggers, 1, logging.INFO,
#                       'Number of object spectra to fit (excluding failed fits): {0}/{1}'.format(
#                            numpy.sum(self.obj_to_fit), len(self.obj_to_fit)))
#
#        # Reject model outliers
#        obj_flux = PPXFFit.reject_model_outliers(obj_flux, result, rescale=False,
#                                                 local_sigma=True, boxcar=self.reject_boxcar,
#                                                 loggers=self.loggers, quiet=self.quiet)
#        obj_ferr[numpy.ma.getmaskarray(obj_flux)] = numpy.ma.masked
#
#        # Return results of refit (only ever do one rejection iteration
#        return self._run_fit_iteration(obj_flux, obj_ferr, obj_to_fit, weight_errors=True,
#                                       component_fits=True, plot=plot)
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO,
                       'Number of object spectra to fit: {0}/{1}'.format(
                            numpy.sum(self.obj_to_fit), len(self.obj_to_fit)))
        result = self._run_fit_iteration(self.obj_flux, self.obj_ferr, self.obj_to_fit,
                                         weight_errors=True,
                                         component_fits=True, plot=plot)
        return result


    def _emission_line_only_model(self, result):
        r"""
        Given the set of :class:`mangadap.proc.ppxf.PPXFFitResult`
        instances resulting from the fit to all object spectra,
        construct the best-fitting models that only include the model
        emission lines.  The :class:`mangadap.proc.ppxf.PPXFFitResult`
        instances must have had the best-fitting component models
        constructed when instantiated.  See :func:`_run_fit_iteration`
        and :func:`_save_results`.

        $$$$$$$$ (FB) I think this function genrates the emission line model from
        the results object. The model is a flux vector made up by a sum of Gaussians
        I imagine result[i].bestfit_compis a Gaussian. Need to check with Kyle

        Args:
            result (numpy.ndarray): An array of
                :class:`mangadap.proc.ppxffit.PPXFFitResult` instances
                with the result of the fits to each object spectrum.

        Returns:
            `numpy.ma.MaskedArray`_: A masked array with shape
            :math:`(N_{\rm spec},N_{\rm pix})` with the best-fitting
            emission-line-only model.  All pixels between
            :attr:`spectrum_start` and :attr:`spectrum_end` are not
            masked.
        """
        # Models originally fully masked
        model_eml_flux = numpy.ma.MaskedArray(numpy.zeros(self.obj_flux.shape, dtype=float),
                                              mask=numpy.ones(self.obj_flux.shape, dtype=bool))
        for i in range(self.nobj):
            if result[i] is None or result[i].fit_failed():
                continue
            s = result[i].start
            e = result[i].end
            # Sum all the emission-line components
            model_eml_flux[i,s:e] = numpy.sum(result[i].bestfit_comp[self.gas_comp,:], axis=0)
        return model_eml_flux


#----------------------------------------------------------------------------
# (FB) THESE FUNCTIONS SET FLAGS ON THE KINEMATICS
#----------------------------------------------------------------------------
    def _is_near_bounds(self, kin, kininp, vel_indx, sig_indx, lbound, ubound, tol_frac=1e-2):
        """
        Check if the fitted kinematics are near the imposed limits.
        
        The definition of "near" is that the velocity and higher moments
        cannot be closer than the provided fraction of the total width
        to the boundary.  For the velocity dispersion, the fraction is
        done in log space.

        (FB) This function sets flags on the kinematics based on vicinty
        to bounds. DOES THIS APPLY TO EMISSION LINES AS WELL? (KBW) Yes.
        Is there any reason it shouldn't?

        Args:
            kin (numpy.ndarray): Best-fitting kinematic parameters.
            kininp (numpy.ndarray): The initial guesses for the
                best-fitting kinematics.  This is needed because the
                velocity limits are set relative to the input guess.
            vel_index (numpy.ndarray): Boolean array setting if the
                parameter is a velocity.  This is needed because the
                velocity limits are set relative to the input guess.
            sig_index (numpy.ndarray): Boolean array setting if the
                parameter is a velocity dispersion.  This is needed
                because the proximity to the bounds is logarithmic for
                the velocity dispersion.
            lbound (numpy.ndarray): Lower bound on each parameter.
            ubound (numpy.ndarray): Upper bound on each parameter.
            tol_frac (float): (**Optional**) The fractional tolerance
                for classifying the returned parameter and near the
                boundary.
        
        Returns:
            numpy.ndarray: Two boolean arrays flagging if the parameter
            is near either boundary and if the parameter is near the
            lower boundary.  The latter is important in case the fitted
            parameter is a velocity dispersion and that it's near the
            lower boundary because it has hit the pixel sampling limit.

        """

        # Offset velocity: bounded by *deviations* from input value
        _lbound = lbound
        _lbound[vel_indx] += kininp[vel_indx]
        _ubound = ubound
        _ubound[vel_indx] += kininp[vel_indx]

        # Set the tolerance
        Db = ubound-lbound
        Db[sig_indx] = numpy.log10(ubound[sig_indx])-numpy.log10(lbound[sig_indx])
        tol = Db*tol_frac

        # Determine if the parameter is near the lower boundary (only
        # relevant to the sigma) ... 
        near_lower_bound = kin - _lbound < tol
        # and whether it's close to either
        near_bound = near_lower_bound | (_ubound - kin < tol)

        # Return the two boundary flags
        return near_bound, near_lower_bound


    def _validate_dispersions(self, model_eml_par, rng=[0,400]):
        """
        Check that the corrected velocity dispersion are in the provided
        range.

        (FB) More flagging of kinematic paramters. DOES THIS APPLY TO
        EMISSION LINES? (KBW): Yes.  The corrected dispsersion must be
        larger than 0 and less than 400 km/s.  It's easy to turn this
        off or change the limits.

        Args:
            model_eml_par (`numpy.recarray`_): Record array with the
                parameters measured for each emission line.  See
                :func:`mangadap.proc.spectralfitting.EmissionLineFit._per_emission_line_dtype`.
            rng (list): (**Optional**) Two-element list with the minimum
                and maximum allowed *corrected* velocity dispersion.
                Measurements outside this range are flagged as
                ``BAD_SIGMA``.

        Returns:
            `numpy.recarray`_: Returns the input record array with any
            additional flags.
        
        """
        _rng = numpy.square(rng)
        _fit_eml = numpy.ones(model_eml_par['SIGMACORR'].shape, dtype=numpy.bool)
        _fit_eml[:,numpy.invert(self.fit_eml)] = False
        sigcor = numpy.square(model_eml_par['KIN'][:,:,1]) \
                        - numpy.square(model_eml_par['SIGMACORR'][:,:])
        indx = ((sigcor < _rng[0]) | (sigcor > _rng[1])) & _fit_eml
        if numpy.sum(indx) == 0:
            return
        model_eml_par['MASK'][indx] = self.bitmask.turn_on(model_eml_par['MASK'][indx], 'BAD_SIGMA')

        return model_eml_par
#----------------------------------------------------------------------------


#----------------------------------------------------------------------------
# (FB) SAVE RESULTS
#----------------------------------------------------------------------------
    def _save_results(self, etpl, result, model_mask, model_fit_par, model_eml_par):
        r"""
        Save and assess the results of the ppxf fits.
        
        The results are saved as direct output from pPXF and parsed into
        the results for each emission line.  The function also assesses
        the data to set any necessary flags.  Much of this is the same
        as :class:`mangadap.proc.ppxffit.PPXFFit._save_results`.

        .. todo::
            - Need to check this.  There could be something wrong with
              parsing the data into the individual emission lines.

        Args:
            etpl (:class:`EmissionLineTemplates`): The object used to
                construct and hold the emission-line templates.
            result (numpy.ndarray): An array of
                :class:`mangadap.proc.ppxffit.PPXFFitResult` instances
                with the result of the fits to each object spectrum.
            model_mask (numpy.ndarray): The array of bitmask values
                associated with spectral fitting flags.  Shape is
                :math:`(N_{\rm spec}, N_{\rm pix})`.
            model_fit_par (`numpy.recarray`_): Record array with the
                dtype as defined by :func:`_per_fit_dtype` that holds
                the ppxf output.
            model_eml_par (`numpy.recarray`_): Record array with the
                individual emission-line data; see
                :func:`mangadap.proc.spectralfitting.EmissionLineFit._per_emission_line_dtype`.

        Returns:
            numpy.ndarray, numpy.recarray: The function returns: (1) the
            best-fitting model spectra, (2) the best-fitting
            emission-line only spectra, (3) the bitmask values, (4) the
            per spectrum ppxf result, and (5) the per spectrum
            emission-line parameters.   The first 3 returned objects are
            of type numpy.ndarray and have shape :math:`(N_{\rm spec},
            N_{\rm pix})`; the last two are numpy.recarray instances
            with shape :math:`(N_{\rm spec},)`.

        """
        #---------------------------------------------------------------
        # Get the model spectra
        model_flux = PPXFFit.compile_model_flux(self.obj_flux, result)
        model_eml_flux = self._emission_line_only_model(result)

        # Save the pixel statistics
        model_fit_par['BEGPIX'] = self.spectrum_start
        model_fit_par['ENDPIX'] = self.spectrum_end
        model_fit_par['NPIXTOT'] = self.spectrum_end - self.spectrum_start

        # Calculate the model residuals, which are masked where the data
        # were not fit
        residual = self.obj_flux - model_flux
        fractional_residual = numpy.ma.divide(self.obj_flux - model_flux, model_flux)
        # Get the chi-square for each spectrum
        model_fit_par['CHI2'] = numpy.sum(numpy.square(residual/self.obj_ferr), axis=1)
        # Get the (fractional) residual RMS for each spectrum
        model_fit_par['RMS'] = numpy.sqrt(numpy.ma.mean(numpy.square(residual), axis=1))
        model_fit_par['FRMS'] = numpy.sqrt(numpy.ma.mean(numpy.square(fractional_residual), axis=1))

        # Flag the pixels that were not used
        model_mask[numpy.ma.getmaskarray(self.obj_flux)] \
                        = self.bitmask.turn_on(model_mask[numpy.ma.getmaskarray(self.obj_flux)],
                                               flag='DIDNOTUSE')

        # Mask any lines that were not fit
        model_eml_par['MASK'][:,numpy.invert(self.fit_eml)] \
                    = self.bitmask.turn_on(model_eml_par['MASK'][:,numpy.invert(self.fit_eml)],
                                           flag='NO_FIT')

        # Generate some convenience data:
        #  - Get the list of indices in the flatted kinematics vectors
        #    with the *unfixed*, *defining* parameters used for each
        #    kinematic measurement.  These are used to set the
        #    kinematics and errors for each emission line.
        #  - Generate vectors with the lower and upper bounds for the
        #    kinematic parameters
        #  - Flag parameters that are velocity and sigma components
        lboundi = [ self.velocity_limits[0], self.sigma_limits[0], self.gh_limits[0],
                    self.gh_limits[0], self.gh_limits[0], self.gh_limits[0] ]
        uboundi = [ self.velocity_limits[1], self.sigma_limits[1], self.gh_limits[1],
                    self.gh_limits[1], self.gh_limits[1], self.gh_limits[1] ]
        lbound = []
        ubound = []
        par_indx = []
        vel_indx = numpy.zeros(self.npar_kin, dtype=bool)
        sig_indx = numpy.zeros(self.npar_kin, dtype=bool)
        for j in range(self.ncomp):
            start = numpy.sum(numpy.absolute(self.comp_moments[:j]))
            nmom = numpy.absolute(self.comp_moments[j])
            par_indx += [ [0]*nmom ]
            for k in range(nmom):
                par_indx[j][k] = start+k if len(self.tied[j][k]) == 0 \
                                        else int(self.tied[j][k].split('[')[1].split(']')[0])
            vel_indx[par_indx[j][0]] = True
            sig_indx[par_indx[j][1]] = True
            lbound += [ lboundi[:nmom] ]
            ubound += [ uboundi[:nmom] ]
        lbound = numpy.concatenate(tuple(lbound))
        ubound = numpy.concatenate(tuple(ubound))

        # The set of gas templates, and the kinematic component,
        # velocity group, and sigma group associated with each
        # template are the same for all fits
        model_fit_par['KINCMP'] = numpy.array([self.tpl_comp]*self.nobj)
        model_fit_par['VELCMP'] = numpy.array([self.tpl_vgrp]*self.nobj)
        model_fit_par['SIGCMP'] = numpy.array([self.tpl_sgrp]*self.nobj)
        model_fit_par['TIEDKIN'] = numpy.array([numpy.concatenate(tuple(par_indx))]*self.nobj)

        #---------------------------------------------------------------
        # Need to iterate over each spectrum
        for i in range(self.nobj):

            #-----------------------------------------------------------
            # Set output flags
            # - No fit was performed
            if result[i] is None:
                model_mask[i,:] = self.bitmask.turn_on(model_mask[i,:], 'NO_FIT')
                continue

            # - No fit attempted because of insufficient data
            if result[i].empty_fit():
                model_mask[i,:] = self.bitmask.turn_on(model_mask[i,:], 'NO_FIT')
                model_fit_par['MASK'][i] = self.bitmask.turn_on(model_fit_par['MASK'][i],
                                                                'INSUFFICIENT_DATA')
                model_eml_par['MASK'][i] = self.bitmask.turn_on(model_eml_par['MASK'][i],
                                                                'INSUFFICIENT_DATA')
                continue

            # - Fit attempted but failed
            if result[i].fit_failed():
                model_mask[i,:] = self.bitmask.turn_on(model_mask[i,:], 'FIT_FAILED')
                model_fit_par['MASK'][i] = self.bitmask.turn_on(model_fit_par['MASK'][i],
                                                                'FIT_FAILED')
                model_eml_par['MASK'][i] = self.bitmask.turn_on(model_eml_par['MASK'][i],
                                                                'FIT_FAILED')

            # - Fit successful but hit maximum iterations.
            if result[i].reached_maxiter():
                model_fit_par['MASK'][i] = self.bitmask.turn_on(model_fit_par['MASK'][i], 'MAXITER')
                model_eml_par['MASK'][i] = self.bitmask.turn_on(model_eml_par['MASK'][i], 'MAXITER')

            # - Mask rejected pixels
            original_gpm = numpy.where(numpy.invert(
                               numpy.ma.getmaskarray(self.obj_flux)[i,self.spectrum_start[i]
                                                                      :self.spectrum_end[i]]))[0]
            rejected_pixels = list(set(original_gpm) - set(result[i].gpm))
            if len(rejected_pixels) > 0:
                model_mask[i,self.spectrum_start[i]:self.spectrum_end[i]][rejected_pixels] \
                        = self.bitmask.turn_on(model_mask[i,self.spectrum_start[i]:
                                                            self.spectrum_end[i]][rejected_pixels],
                                               flag=PPXFFit.rej_flag)

            #-----------------------------------------------------------
            # Save the model parameters and figures of merit
            # - Number of fitted pixels
            model_fit_par['NPIXFIT'][i] = len(result[i].gpm)
            # - Templates used
            model_fit_par['USETPL'][i] = result[i].tpl_to_use
            # - Template weights
            model_fit_par['TPLWGT'][i][result[i].tpl_to_use] = result[i].tplwgt
            # Additive polynomial coefficients
            if self.degree >= 0 and result[i].addcoef is not None:
                model_fit_par['ADDCOEF'][i] = result[i].addcoef
            if self.mdegree > 0 and result[i].multcoef is not None:
                model_fit_par['MULTCOEF'][i] = result[i].multcoef
            # Flattened input kinematics vector
            model_fit_par['KININP'][i] = numpy.concatenate(tuple(self.comp_start_kin[i]))
            # Flattened best-fit kinematics vector
            model_fit_par['KIN'][i] = numpy.concatenate(tuple(result[i].kin))
            # Flattened kinematic error vector
            model_fit_par['KINERR'][i] = numpy.concatenate(tuple(result[i].kinerr))
            # Chi-square
            model_fit_par['RCHI2'][i] = model_fit_par['CHI2'][i] \
                                        / (model_fit_par['NPIXFIT'][i] 
                                            - self.dof - numpy.sum(model_fit_par['TPLWGT'][i] > 0))
            model_fit_par['ROBUST_RCHI2'][i] = result[i].robust_rchi2

            # Get growth statistics for the residuals
            model_fit_par['ABSRESID'][i] = residual_growth((residual[i,:]).compressed(),
                                                       [0.68, 0.95, 0.99])
            model_fit_par['FABSRESID'][i] = residual_growth(fractional_residual[i,:].compressed(),
                                                        [0.68, 0.95, 0.99])

            #-----------------------------------------------------------
            # Test if the kinematics are near the imposed boundaries.
            near_bound, near_lower_bound = self._is_near_bounds(model_fit_par['KIN'][i],
                                                                model_fit_par['KININP'][i],
                                                                vel_indx, sig_indx, lbound, ubound)

            # Add the *global* flag for the fit.
            # TODO: These are probably too general.
            # - If the velocity dispersion has hit the lower limit, ONLY
            #   flag the value as having a MIN_SIGMA.
            if numpy.any(near_lower_bound & sig_indx):
                model_fit_par['MASK'][i] = self.bitmask.turn_on(model_fit_par['MASK'][i],
                                                                'MIN_SIGMA')
            # - Otherwise, flag the full fit as NEAR_BOUND, both the
            #   parameters and the model
            if numpy.any((near_lower_bound & numpy.invert(sig_indx)) 
                                | (near_bound & numpy.invert(near_lower_bound))):
                model_fit_par['MASK'][i] = self.bitmask.turn_on(model_fit_par['MASK'][i],
                                                                'NEAR_BOUND')
                model_mask[i,:] = self.bitmask.turn_on(model_mask[i,:], 'NEAR_BOUND')

            # Convert the velocities from pixel units to cz
            model_fit_par['KININP'][i,vel_indx], _ \
                        = PPXFFit.convert_velocity(model_fit_par['KININP'][i,vel_indx],
                                                   numpy.zeros(numpy.sum(vel_indx)))
            model_fit_par['KIN'][i,vel_indx], model_fit_par['KINERR'][i,vel_indx] \
                        = PPXFFit.convert_velocity(model_fit_par['KIN'][i,vel_indx],
                                                   model_fit_par['KINERR'][i,vel_indx])

            # Divvy up the fitted parameters into the result for each
            # emission line
            for j in range(self.neml):
                if not self.fit_eml[j]:
                    continue

                # The "fit index" is the component of the line
                model_eml_par['FIT_INDEX'][i,j] = self.eml_compi[j]

                # EmissionLineTemplates constructs each line to have the
                # flux provided by the emission-line database
                model_eml_par['FLUX'][i,j] = result[i].tplwgt[self.eml_compi[j]] \
                                                * self.emldb['flux'][j]
                model_eml_par['FLUXERR'][i,j] = result[i].tplwgterr[self.eml_compi[j]] \
                                                * self.emldb['flux'][j]

                # Use the flattened vectors to set the kinematics
                indx = par_indx[self.eml_compi[j]]
                model_eml_par['KIN'][i,j,:] = model_fit_par['KIN'][i,indx]
                model_eml_par['KINERR'][i,j,:] = model_fit_par['KINERR'][i,indx]

                # Get the bound masks specific to this emission-line (set)
                if numpy.any(near_lower_bound[indx] & sig_indx[indx]):
                    model_eml_par['MASK'][i,j] = self.bitmask.turn_on(model_eml_par['MASK'][i,j],
                                                                      'MIN_SIGMA')
                if numpy.any((near_lower_bound[indx] & numpy.invert(sig_indx[indx])) 
                                | (near_bound[indx] & numpy.invert(near_lower_bound[indx]))):
                    model_eml_par['MASK'][i,j] = self.bitmask.turn_on(model_eml_par['MASK'][i,j],
                                                                      'NEAR_BOUND')

            # Get the instrumental dispersion in the galaxy data at the
            # location of the fitted lines
            sigma_inst = EmissionLineFit.instrumental_dispersion(self.obj_wave, self.obj_sres[i,:],
                                                        self.emldb['restwave'][self.fit_eml],
                                                        model_eml_par['KIN'][i,self.fit_eml,0])
#            pyplot.scatter(self.emldb['restwave'][self.fit_eml], sigma_inst, marker='.', s=50)
#            pyplot.scatter(self.emldb['restwave'][self.fit_eml], etpl.eml_sigma_inst[self.fit_eml], marker='.', s=50)
#            pyplot.show()

            # The dispersion correction is the quadrature difference
            # between the instrumental dispersion in the galaxy data to
            # the dispersion used when constructing the emission-line
            # templates
            sigma2corr = numpy.square(sigma_inst) - numpy.square(etpl.eml_sigma_inst[self.fit_eml])
            if numpy.any(sigma2corr < 0):
                print(sigma2corr)
                warnings.warn('Encountered negative sigma corrections!')
            model_eml_par['SIGMACORR'][i,self.fit_eml] = numpy.ma.sqrt(numpy.square(sigma_inst)
                                    - numpy.square(etpl.eml_sigma_inst[self.fit_eml])).filled(0.0)
#            print(model_eml_par['SIGMACORR'][i])

        #---------------------------------------------------------------
        # Test if kinematics are reliable
        model_eml_par = self._validate_dispersions(model_eml_par)

        #---------------------------------------------------------------
        # Return the fitting results
        # - model_flux: full model fit to the spectra
        # - model_eml_flux: emission-line only model
        # - model_mask: Bitmask spectra for the fit
        # - model_fit_par: The saved results from the ppxf fit
        # - model_eml_par: The fit results parsed into data for each
        #   emission line
        return model_flux, model_eml_flux, model_mask, model_fit_par, model_eml_par



    def fit_SpatiallyBinnedSpectra(self, binned_spectra, par=None, loggers=None, quiet=False):
        """
        This DAP-specific function interprets the DAP-specific classes
        and constructs call(s) to the general :func:`fit` function to
        fit the spectra.
       
        The format of the calling sequence is dictated by the DAP
        interface: Any emission-line fitter (e.g.
        :class:`mangadap.proc.elric.Elric`) must have the same function
        interface in order to be called within the
        :class:`mangadap.proc.emissionlinemodel.EmissionLineModel`
        object.  Because this is a DAP-specific function, it should not
        declare anything to self.

        The returned emission-line "baseline" array is set to be the
        difference between the best-fitting stellar continuum (passed to
        this function as par['stellar_continuum']) and the best-fitting
        stellar continuum from the combined emission-line and stellar
        spectra produced by :func:`fit`.  In this way, the best-fitting
        model for each spectrum is::
        
            best_fit_model = par['stellar_continuum']['FLUX'].data \
                                + model_eml_flux + model_eml_base

        where `model_eml_flux` and `model_eml_base` are the 2nd and 3rd
        returned objects, respectively.  These are written to the output
        DAP model LOGCUBE file in the extensions EMLINE and EMLINE_BASE,
        respectively.

        .. todo::
            - Use waverange in pixel mask to restrict wavelength range.
              Add to SasukePar.
            - We probably need to rethink the emission-line "baseline"
              output.

        Args:
            binned_spectra
                (:class:`mangadap.proc.spatiallybinnedspectra.SpatiallyBinnedSpectra`):
                Spectra to fit.
            par (:class:`SasukePar`): Parameters provided from the DAP
                to the general Sasuke fitting algorithm (:func:`fit`).
                Althought technically optional given that it is a
                keyword parameter, the :class:`SasukePar` parameter must
                be provided for proper execution of the function.
            loggers (list): (**Optional**) List of `logging.Logger`_ objects
                to log progress; ignored if quiet=True.  Logging is done
                using :func:`mangadap.util.log.log_output`.  Default is
                no logging.
            quiet (bool): (**Optional**) Suppress all terminal and
                logging output.  Default is False.

        Returns:
            numpy.ndarray, numpy.recarray: The function returns: (1)
            wavelength vector of the fitted models, which should match
            the binned spectra wavelength vector
            (binned_spectra['WAVE'].data), (2) the best-fitting
            emission-line model model spectra, (3) the best-fitting
            emission-line baseline (see the description above), (4) the
            bitmask values, (5) the per spectrum ppxf result, and (6)
            the per spectrum emission-line parameters.   The first
            object is a numpy.ndarray instance with shape :math:`(N_{\rm
            pix},)` , the next 3 objects are numpy.ndarray instances
            with shape :math:`(N_{\rm spec}, N_{\rm pix})`, and the last
            two are numpy.recarray instances with shape :math:`(N_{\rm
            spec},)`.
        """
        # Assign the parameters if provided
        if par is None:
            raise ValueError('Must provide parameters!')
        if not isinstance(par, SasukePar):
            raise TypeError('Input parameters must be an instance of SasukePar.')
        # SasukePar checks the types of the stellar continuum,
        # emission-line database, and pixel mask

        # SpatiallyBinnedSpectra object always needed
        if binned_spectra is None:
            raise ValueError('Must provide spectra object for fitting.')
        if not isinstance(binned_spectra, SpatiallyBinnedSpectra):
            raise TypeError('Must provide a valid SpatiallyBinnedSpectra object!')
        if binned_spectra.hdu is None:
            raise ValueError('Provided SpatiallyBinnedSpectra object is undefined!')

        # Get the data arrays to fit
        # TODO: This could be where we pull out the individual spaxels
        # and renormalize the continuum fit from the binned data.  For
        # now this just pulls out the binned spectra
        flux, ferr = EmissionLineFit.get_spectra_to_fit(binned_spectra, pixelmask=par['pixelmask'],
                                                        error=True)
        
                                                        
#        binned_spectra['FLUX']
        
#        binned_spectra['BINDID'] 2D map
#        binned_spectra['BINDID'] 2D map
        # TODO: Also may want to include pixels rejected during stellar
        # kinematics fit
        nobj = flux.shape[0]
#        print(nobj)

        # Get the stellar templates
        # TODO: This could be where we instead construct the
        # optimal template used for each spectrum.
        stellar_templates = None if par['stellar_continuum'] is None else \
                                par['stellar_continuum'].method['fitpar']['template_library']
        stpl_wave = None if par['stellar_continuum'] is None else stellar_templates['WAVE'].data
        stpl_flux = None if par['stellar_continuum'] is None else stellar_templates['FLUX'].data
        if not quiet:
            warnings.warn('Adopting mean spectral resolution of all templates!')
        stpl_sres = None if par['stellar_continuum'] is None \
                        else numpy.mean(stellar_templates['SPECRES'].data, axis=0).ravel()
        velscale_ratio = 1 if par['stellar_continuum'] is None \
                                else par['stellar_continuum'].method['fitpar']['velscale_ratio']
        matched_resolution = False if par['stellar_continuum'] is None \
                                else par['stellar_continuum'].method['fitpar']['match_resolution']

        # Get the stellar kinematics
        # TODO: This could be where we pull out the valid fits and then
        # interpolate (or some other approach) to spaxels without direct
        # stellar kinematics measurements.  For now this just pulls out
        # the data for the binned spectra and replaces those results for
        # spectra with insufficient S/N or where pPXF failed to the
        # median redshift/dispersion of the fitted spectra.  We'll need
        # to revisit this

        stellar_velocity, stellar_dispersion = (None, None) if par['stellar_continuum'] is None \
                        else par['stellar_continuum'].matched_guess_kinematics(binned_spectra,
                                                                               cz=True)
        stellar_kinematics = None if stellar_velocity is None or stellar_dispersion is None \
                                else numpy.array([ stellar_velocity, stellar_dispersion ]).T
#        print(stellar_redshift)
#        print(stellar_dispersion)
#        print(stellar_redshift.shape)

        # Set which stellar templates to use for each spectrum
        # TODO: Here I'm flagging to only use the stellar templates used
        # in the stellar kinematics fit.  You could select to use all
        # templates by setting stpl_to_use=None below, or we could construct
        # the optimal template for each spectrum and just have a single
        # stellar template used in the emission-line fits.
        stpl_to_use = None if par['stellar_continuum'] is None \
                        else par['stellar_continuum'].matched_template_flags(binned_spectra)
#        stpl_to_use = None

        # Get the spectra that meet the S/N criterion
        # TODO: This could be based on the moment assessment of the
        # emission-line S/N instead; for now just based on continuum
        # S/N.
        good_snr = binned_spectra.above_snr_limit(par['minimum_snr'])

        # Determine which spectra have a valid stellar continuum fit
        good_stellar_continuum_fit = numpy.invert(
                    par['stellar_continuum'].bitmask.flagged(
                                            par['stellar_continuum']['PAR'].data['MASK'],
                                            flag=[ 'NO_FIT', 'INSUFFICIENT_DATA', 'FIT_FAILED']))

        # TODO: At the moment, spectra to fit must have both good S/N
        # and a stellar continuum fit
        spec_to_fit = good_snr & good_stellar_continuum_fit

        # TODO: For now can only fit two moments
        if par['moments'] != 2:
            print(par['moments'])
            raise NotImplementedError('Number of gas moments can only be two.')

        # Return the fitted data
        model_wave, model_flux, model_eml_flux, model_mask, model_fit_par, model_eml_par \
                = self.fit(par['emission_lines'], binned_spectra['WAVE'].data, flux[spec_to_fit,:],
                           obj_ferr=ferr[spec_to_fit,:],
                           obj_sres=binned_spectra['SPECRES'].data.copy()[spec_to_fit,:],
                           guess_redshift=par['guess_redshift'][spec_to_fit],
                           guess_dispersion=par['guess_dispersion'][spec_to_fit],
                           reject_boxcar=par['reject_boxcar'],
                           stpl_wave=stpl_wave, stpl_flux=stpl_flux, stpl_sres=stpl_sres,
                           stpl_to_use=None if stpl_to_use is None else stpl_to_use[spec_to_fit,:],
                           stellar_kinematics=None if stellar_kinematics is None else
                                        stellar_kinematics[spec_to_fit,:],
                           velscale_ratio=velscale_ratio, matched_resolution=matched_resolution,
                           bias=par['bias'], degree=par['degree'], mdegree=par['mdegree'],
                           #moments=par['moments'],
                           loggers=loggers, quiet=quiet)
        # Save the the bin ID numbers indices based on the spectra
        # selected to be fit
        model_fit_par['BINID'] = binned_spectra['BINS'].data['BINID'][spec_to_fit]
        model_fit_par['BINID_INDEX'] = numpy.arange(binned_spectra.nbins)[spec_to_fit]

        model_eml_par['BINID'] = binned_spectra['BINS'].data['BINID'][spec_to_fit]
        model_eml_par['BINID_INDEX'] = numpy.arange(binned_spectra.nbins)[spec_to_fit]

        # Add the equivalent width data
        EmissionLineFit.measure_equivalent_width(binned_spectra['WAVE'].data, flux[spec_to_fit,:],
                                                 par['emission_lines'], model_eml_par,
                                                 redshift=par['guess_redshift'][spec_to_fit],
                                                 bitmask=self.bitmask, checkdb=False)
        eml_continuum = model_flux-model_eml_flux

        # Use the previous fit to the stellar continuum to set the
        # emission-line model "baseline"

        # TODO: This is done as a place holder.  We need a better way of
        # propagating the difference between the stellar-kinematics fit
        # and the combined stellar-continuum + emission-line model to
        # the output datamodel.
        stellar_continuum = numpy.zeros(flux.shape, dtype=float) \
                            if par['stellar_continuum'] is None \
                            else par['stellar_continuum'].fill_to_match(binned_spectra).filled(0.0)

        model_eml_base = model_flux-model_eml_flux-stellar_continuum[spec_to_fit,:]
        model_eml_flux += model_eml_base

#        pyplot.plot(model_wave, flux[spec_to_fit,:][0,:], color='k', lw=1, zorder=1)
#        pyplot.plot(model_wave, stellar_continuum[spec_to_fit,:][0,:], color='C0', lw=1.0,
#                    zorder=2)
#        pyplot.plot(model_wave, model_eml_flux[0,:], color='C1', lw=1.0, zorder=3)
#        pyplot.plot(model_wave, model_eml_base[0,:], color='C2', lw=1.0, zorder=4)
#        pyplot.plot(model_wave, eml_continuum[0,:], color='C3', lw=1.0, zorder=4)
#        pyplot.show()
        
        # Only return model and model parameters for the *fitted*
        # spectra
        return model_wave, model_eml_flux, model_eml_base, model_mask, model_fit_par, model_eml_par


    def fit(self, emission_lines, obj_wave, obj_flux, obj_ferr=None, mask=None, obj_sres=None,
            guess_redshift=None, guess_dispersion=None, reject_boxcar=None, stpl_wave=None,
            stpl_flux=None, stpl_sres=None, stpl_to_use=None, stellar_kinematics=None,
            velscale_ratio=None, matched_resolution=True, waverange=None, bias=None, degree=4,
            mdegree=0, max_velocity_range=400., alias_window=None, dvtol=1e-10, loggers=None,
            quiet=False, plot=False):
            #moments=2,
        r"""
        Fit a set of emission lines using pPXF to all provided spectra.

        The input flux arrays are expected to be ordered with spectra
        along **rows**.  This is opposite to how they're provided to
        ppxf.  I.e., to plot the first spectrum in the array, one would
        do::
            
            from matplotlib import pyplot
            pyplot.plot(obj_wave, obj_flux[0,:])
            pyplot.show()

        The function should fit the spectra with or without any provided
        stellar templates.

        The body of this function mostly deals with checking the input
        and setting up the template and object data for use with pPXF.

        The function is meant to be general, but has been tested only on
        MaNGA spectra so far.

        If the spectral resolution is not matched between the templates
        and the object spectra, the provided stellar_kinematics is
        expected to include any resolution difference; i.e., it is
        **not** the astrophysical velocity dispersion.

        .. todo::
            - guess_redshift and guess_dispersion are not actually
              optional.  The function will fail without them!
            - Allow for moments != 2.
            - Allow for fixed components to be set from emission-line
              database
            - Allow for bounds to be set from emission-line database

        Args:
            emission_lines
                (:class:`mangadap.par.emissionlinedb.EmissionLineDB'):
                Emission-line database that is parsed to construct the
                emission-line templates to fit (see
                :class:`EmissionLineTemplates`).
            obj_wave (numpy.ndarray): Wavelength vector for object
                spectra.  Shape is (:math:`N_{\rm pix}`,).
            obj_flux (numpy.ndarray): Object spectra to fit.  Can be
                provided as a `numpy.ma.MaskedArray`_.  Shape is
                (:math:`N_{\rm spec},N_{\rm pix}`).
            obj_ferr (numpy.ndarray): (**Optional**) :math:`1\sigma`
                errors in object spectra.  Can be provided as a
                `numpy.ma.MaskedArray`_.  Shape is (:math:`N_{\rm
                spec},N_{\rm pix}`).  If None, the quadrature sum of the
                fit residuals is used as the fit metric instead of
                chi-square (all errors are set to unity).
            mask (numpy.ndarray or
                :class:`mangadap.util.pixelmask.SpectralPixelMask`):
                (**Optional**) Boolean array or
                :class:`mangadap.util.pixelmask.SpectralPixelMask`
                instance used to censor regions of the spectra to ignore
                during fitting.
            obj_sres (numpy.ndarray): (**Optional**) The spectral
                resolution of the object data.  Can be a single vector
                for all spectra or one vector per object spectrum.
            guess_redshift (array-like): (**Optional**) The starting
                redshift guess.  Can provide a single value for all
                spectra or one value per spectrum.
            guess_dispersion (array-like): (**Optional**) The starting
                velocity dispersion guess.  Can provide a single value
                for all spectra or one value per spectrum.
            reject_boxcar (int): (**Optional**) Size of the boxcar to
                use when rejecting fit outliers.  If None, no rejection
                iteration is performed.
            stpl_wave (numpy.ndarray): (**Optional**) Wavelength vector
                for stellar template spectra.  Shape is (:math:`N_{{\rm
                pix},\ast}`,).
            stpl_flux (numpy.ndarray): (**Optional**) Stellar template
                spectra to use in fit.  Shape is (:math:`N_{{\rm
                tpl},\ast},N_{{\rm pix},\ast}`).
            stpl_sres (numpy.ndarray): (**Optional**) Spectral
                resolution of the stellar template spectra.  All
                templates are assumed to have the same spectral
                resolution.  Shape is (:math:`N_{{\rm pix},\ast}`,).
                This is also used to set the spectral resolution of the
                emission-line templates.  If not provided, the
                emission-line templates are constructed with an LSF that
                has a disperison of 1 pixel.
            stpl_to_use (numpy.ndarray): (**Optional**) Set of flags
                used to select the stellar templates used to fit each
                spectrum.  Shape is (:math:`N_{\rm spec},N_{{\rm
                tpl},\ast}`).
            stellar_kinematics (numpy.ndarray): (**Optional**) The
                kinematics to use for the stellar component.  These
                kinematics are **fixed** for all calls to ppxf.  Shape
                is (:math:`N_{\rm spec},N_{{\rm kin},\ast}`).  The shape
                of this array determines the number of moments assigned
                to the stellar component.
            velscale_ratio (int): (**Optional**) The **integer** ratio
                between the velocity scale of the pixel in the galaxy
                data to that of the template data.  If None, set to
                unity.
            matched_resolution (bool): (**Optional**) The spectral
                resolution of the templates is matched to that of the
                galaxy data.  WARNING: This functionality needs to be
                checked in relation to the gas templates!
            waverange (array-like): (**Optional**) Lower and upper
                wavelength limits to *include* in the fit.  This can be
                a two-element vector to apply the same limits to all
                spectra, or a N-spec x 2 array with wavelength ranges
                for each spectrum to be fit.  Default is to use as much
                of the spectrum as possible.
            bias (float): (**Optional**) pPXF bias parameter.
                (Currently irrelevant because gas is currently always
                fit with moments=2.)
            degree (int): (**Optional**) pPXF degree parameter setting
                the degree of the additive polynomials to use,
                :math:`o_{\rm add}`.
            mdegree (int): (**Optional**) pPXF mdegree parameter setting
                the degree of the multiplicative polynomials to use,
                :math:`o_{\rm mult}`.
            max_velocity_range (float): (**Optional**) Maximum range
                (+/-) expected for the fitted velocities in km/s.
                Default is 400 km/s.
            alias_window (float) : (**Optional**) The window to mask to
                avoid aliasing near the edges of the spectral range in
                km/s.  Default is six times *max_velocity_range*.
            dvtol (float): (**Optional**) The velocity scale of the
                template spectra and object spectrum must be smaller
                than this tolerance.  Default is 1e-10.
            plot (bool): (**Optional**) Show the automatically generated
                pPXF fit plots for each iterations of each spectrum.
            loggers (list): (**Optional**) List of `logging.Logger`_
                objects to log progress; ignored if quiet=True.  Logging
                is done using :func:`mangadap.util.log.log_output`.
            quiet (bool): (**Optional**) Suppress all terminal and
                logging output.

        Returns:
            numpy.ndarray, numpy.recarray: The function returns: (1) the
            wavelength vector for the model spectra (should be the same
            as obj_wave), (2) the best-fitting model spectra, (3) the
            best-fitting emission-line only spectra, (4) the bitmask
            values, (5) the per spectrum ppxf result, and (6) the per
            spectrum emission-line parameters.  The first object is a
            numpy.ndarray instance with shape :math:`(N_{\rm pix},)`,
            the next 3 objects are numpy.ndarray instances with shape
            :math:`(N_{\rm spec}, N_{\rm pix})`, and the last two are
            numpy.recarray instances with shape :math:`(N_{\rm spec},)`.

        Raises:
            ValueError: Raised if the length of the spectra, errors, or
                mask does not match the length of the wavelength array;
                raised if the wavelength, redshift, or dispersion arrays
                are not 1D vectors; and raised if the number of
                redshifts or dispersions is not a single value or the
                same as the number of input spectra.
        """
        #---------------------------------------------------------------
        # Initialize the reporting
        if loggers is not None:
            self.loggers = loggers
        self.quiet = quiet

        #---------------------------------------------------------------
        # Check the input data
        self.obj_wave, self.obj_flux, self.obj_ferr, self.obj_sres \
                = PPXFFit.check_objects(obj_wave, obj_flux, obj_ferr=obj_ferr, obj_sres=obj_sres)
        self.nobj, self.npix_obj = self.obj_flux.shape
        self.waverange = PPXFFit.set_wavelength_range(self.nobj, self.obj_wave, waverange)
        self.input_obj_mask = numpy.ma.getmaskarray(self.obj_flux).copy()
        self.obj_to_fit = numpy.any(numpy.invert(self.input_obj_mask), axis=1)
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO,
                       'Number of object spectra to fit: {0}/{1}'.format(
                            numpy.sum(self.obj_to_fit), len(self.obj_to_fit)))
        self.input_cz, guess_kin = PPXFFit.check_input_kinematics(self.nobj, guess_redshift,
                                                                  guess_dispersion)

        #---------------------------------------------------------------
        # Compare pixel scales and set template wavelength vector
        if stpl_wave is not None:
            self.velscale, self.velscale_ratio \
                    = PPXFFit.check_pixel_scale(stpl_wave, self.obj_wave,
                                                velscale_ratio=velscale_ratio, dvtol=dvtol)
            self.tpl_wave = stpl_wave
        else:
            self.velscale = spectrum_velocity_scale(self.obj_wave)
            self.velscale_ratio = 1
            self.tpl_wave = self.obj_wave

        #---------------------------------------------------------------
        # Check any input stellar template spectra
        # self.tpl_sres has type
        # mangadap.util.instrument.SpectralResolution!
        if stpl_flux is not None:
            if stpl_wave is None:
                raise ValueError('Must provide wavelengths if providing stellar template fluxes.')
            self.tpl_wave, self.tpl_flux, self.tpl_sres \
                    = PPXFFit.check_templates(stpl_wave, stpl_flux, tpl_sres=stpl_sres,
                                              velscale_ratio=self.velscale_ratio)
            self.nstpl = self.tpl_flux.shape[0]
            # Check or instantiate the fit flags
            self.tpl_to_use = PPXFFit.check_template_usage_flags(self.nobj, self.nstpl, stpl_to_use)
        else:
            self.tpl_flux = None
            self.tpl_sres = None
            self.nstpl = 0
            self.tpl_to_use = None

        #---------------------------------------------------------------
        # Set the template spectral resolution.  This is needed to
        # construct the emission-line templates:
        # - Set to the stellar template resolution of that is provided
        #   (done above)
        # - If the object resolution is not available, set such that the
        #   template lines will be initialized with sigma = 1 pixel.
        # - If the object resolution is available, set to the minimum
        #   object resolution at each wavelength channel.
        convertR = astropy.constants.c.to('km/s').value / DAPConstants.sig2fwhm
        self.matched_resolution = matched_resolution
        if self.tpl_sres is None:
            self.matched_resolution = False
            self.tpl_sres = numpy.full(self.npix_obj, convertR/self.velscale, dtype=float) \
                                if self.obj_sres is None else numpy.amin(self.obj_sres, axis=1)
            self.tpl_sres = SpectralResolution(self.tpl_wave, self.tpl_sres, log10=True)

        #---------------------------------------------------------------
        # If provided, check the shapes of the stellar kinematics
        if stellar_kinematics is not None and stellar_kinematics.shape[0] != self.nobj:
            raise ValueError('Provided kinematics do not match the number of input object spectra.')
        stellar_moments = None if stellar_kinematics is None else stellar_kinematics.shape[1]
        if self.nstpl > 0 and stellar_kinematics is None:
            raise ValueError('Must provide stellar kinematics if refiting stellar templates.')
        # Convert velocities from cz to pPXF pixelized velocities
        if self.nstpl > 0:
            _stellar_kinematics = stellar_kinematics.copy()
            _stellar_kinematics[:,0], _ = PPXFFit.revert_velocity(stellar_kinematics[:,0],
                                                                  numpy.zeros(self.nobj))

        #---------------------------------------------------------------
        # Build the emission-line templates; the EmissionLineTemplates
        # object will check the database
        self.emldb = emission_lines
        self.neml = self.emldb.neml
        etpl = EmissionLineTemplates(self.tpl_wave, convertR/self.tpl_sres.sres(), emldb=self.emldb,
                                     loggers=self.loggers, quiet=self.quiet)
        # Set the component associated with each emission line emission line
        self.fit_eml = self.emldb['action'] == 'f'
        self.eml_tpli = etpl.tpli.copy()
        self.eml_compi = numpy.full(self.neml, -1, dtype=int)
        self.eml_compi[self.fit_eml] = numpy.array([ etpl.comp[i] for i in etpl.tpli[self.fit_eml]])

        #---------------------------------------------------------------
        # Save the basic pPXF parameters
        # TODO: Use the emission-line database to set the number of
        # moments to fit to each emission-line component.  For now it is
        # always moments=2!
        self.velocity_limits, self.sigma_limits, self.gh_limits \
                    = PPXFFit.losvd_limits(self.velscale)
        self.bias = bias
        self.degree = degree
        self.mdegree = mdegree
        moments = 2                #numpy.absolute(moments)
        self.reject_boxcar = reject_boxcar
#        self.fix_kinematics = False     #moments < 0

        #---------------------------------------------------------------
        # Compile the template fluxes, components, and velocity and
        # sigma groups
        if self.nstpl == 0:
            self.tpl_flux = etpl.flux
            self.tpl_to_use = numpy.ones((self.nobj,etpl.ntpl), dtype=numpy.bool)
            self.tpl_comp = etpl.comp
            self.tpl_vgrp = etpl.vgrp
            self.tpl_sgrp = etpl.sgrp
            self.ncomp = numpy.amax(self.tpl_comp)+1
            self.comp_moments = numpy.array([moments]*self.ncomp)
            self.comp_start_kin = numpy.array([[gk.tolist()]*self.ncomp for gk in guess_kin ])
        else:
            self.tpl_flux = numpy.append(self.tpl_flux, etpl.flux, axis=0)
            self.tpl_to_use = numpy.append(self.tpl_to_use,
                                           numpy.ones((self.nobj,etpl.ntpl), dtype=numpy.bool),
                                           axis=1)
            self.tpl_comp = numpy.append(numpy.zeros(self.nstpl, dtype=int), etpl.comp+1)
            self.tpl_vgrp = numpy.append(numpy.zeros(self.nstpl, dtype=int), etpl.vgrp+1)
            self.tpl_sgrp = numpy.append(numpy.zeros(self.nstpl, dtype=int), etpl.sgrp+1)
            self.eml_tpli[self.fit_eml] += 1
            self.eml_compi[self.fit_eml] += 1
            self.ncomp = numpy.amax(self.tpl_comp)+1
            self.comp_moments = numpy.array([-stellar_moments] + [moments]*(self.ncomp-1))
            self.comp_start_kin = numpy.array([ [sk.tolist()] + [gk.tolist()]*(self.ncomp-1) 
                                            for sk,gk in zip(_stellar_kinematics, guess_kin) ])
        self.ntpl, self.npix_tpl = self.tpl_flux.shape
        self.tpl_npad = fftpack.next_fast_len(self.npix_tpl)
        self.tpl_rfft = numpy.fft.rfft(self.tpl_flux, self.tpl_npad, axis=1)
        print(self.comp_moments)
        print(self.comp_moments.shape)
        exit()

        # Set which components are gas components
        self.gas_comp = numpy.ones(self.ncomp, dtype=bool)
        if self.nstpl > 0:
            self.gas_comp[0] = False

        #---------------------------------------------------------------
        # Parse the velocity and sigma groups into tied parameters
        self.npar_kin = numpy.sum(numpy.absolute(self.comp_moments))
        self.tied = numpy.empty(self.npar_kin, dtype=object)
        tpl_index = numpy.arange(self.ntpl)
        for i in range(self.ncomp):
            # Do not allow tying to fixed components?
            if self.comp_moments[i] < 0:
                continue
            # Velocity group of this component
            indx = self.tpl_comp[tpl_index[self.tpl_vgrp == i]]
            if len(indx) > 1:
                parn = [ 0 + numpy.sum(numpy.absolute(self.comp_moments[:i])) for i in indx ]
                self.tied[parn[1:]] = 'p[{0}]'.format(parn[0])
            
            # Sigma group of this component
            indx = self.tpl_comp[tpl_index[self.tpl_sgrp == i]]
            if len(indx) > 1:
                parn = [ 1 + numpy.sum(numpy.absolute(self.comp_moments[:i])) for i in indx ]
                self.tied[parn[1:]] = 'p[{0}]'.format(parn[0])

        self.tied[[t is None for t in self.tied ]] = ''
        self.nfree_kin = numpy.sum([len(t) == 0 for t in self.tied])

        # Get the degrees of freedom (accounting for fixed stellar
        # component)
        self.dof = self.nfree_kin + numpy.sum(self.comp_moments[self.comp_moments < 0]) \
                        + max(self.mdegree, 0)
        if self.degree >= 0:
            self.dof += self.degree+1

        # Check if tying parameters is needed
        if self.nfree_kin == self.npar_kin:
            self.tied = None
        else:
            start = [numpy.sum(numpy.absolute(self.comp_moments[:i])) for i in range(self.ncomp) ]
            self.tied = [ self.tied[start[i] : 
                          start[i]+numpy.absolute(self.comp_moments[i])].tolist()
                            for i in range(self.ncomp) ]

        #---------------------------------------------------------------
        # Report the input checks/results
#        if not self.quiet:
#            log_output(self.loggers, 1, logging.INFO, 'Pixel scale: {0} km/s'.format(self.velscale))
#            log_output(self.loggers, 1, logging.INFO, 'Pixel scale ratio: {0}'.format(
#                                                                            self.velscale_ratio))
#            log_output(self.loggers, 1, logging.INFO, 'Dispersion limits: {0} - {1}'.format(
#                                                                            *self.sigma_limits))
#            log_output(self.loggers, 1, logging.INFO, 'Model degrees of freedom: {0}'.format(
#                                                                            self.dof+self.ntpl))
#            log_output(self.loggers, 1, logging.INFO, 'Number of tied parameters: {0}'.format(
#                            0 if self.tied is None else
#                            self.nfree_kin + numpy.sum(self.comp_moments[self.comp_moments < 0])))

        #---------------------------------------------------------------
        # Initialize the output arrays.  This is done here for many of
        # these objects only in case PPXFFit.initialize_pixels_to_fit()
        # fails and the code returns without running the pPXF fitting.
        #  - Model flux
        model_flux = numpy.zeros(self.obj_flux.shape, dtype=numpy.float)
        model_eml_flux = numpy.zeros(self.obj_flux.shape, dtype=numpy.float)
        #  - Model mask:
        model_mask = numpy.zeros(self.obj_flux.shape, dtype=self.bitmask.minimum_dtype())
        indx = numpy.ma.getmaskarray(self.obj_flux)
        model_mask[indx] = self.bitmask.turn_on(model_mask[indx], 'DIDNOTUSE')
        #  - Model parameters and fit quality
        model_fit_par = init_record_array(self.nobj,
                                          self._per_fit_dtype(self.ntpl, self.degree+1,
                                          self.mdegree, self.npar_kin,
                                          self.bitmask.minimum_dtype()))
        model_fit_par['BINID'] = numpy.arange(self.nobj)
        model_fit_par['BINID_INDEX'] = numpy.arange(self.nobj)
        #  - Model emission-line parameters
        model_eml_par = init_record_array(self.nobj,
                                          self._per_emission_line_dtype(self.neml, 2,
                                                                self.bitmask.minimum_dtype()))
        model_eml_par['BINID'] = numpy.arange(self.nobj)
        model_eml_par['BINID_INDEX'] = numpy.arange(self.nobj)

        #---------------------------------------------------------------
        # Initialize the mask and the spectral range to fit
        model_mask, err, self.spectrum_start, self.spectrum_end \
                = PPXFFit.initialize_pixels_to_fit(self.tpl_wave, self.obj_wave, self.obj_flux,
                                                   self.obj_ferr, self.velscale,
                                                   velscale_ratio=self.velscale_ratio,
                                                   waverange=waverange, mask=mask,
                                                   bitmask=self.bitmask,
                                                   velocity_offset=self.input_cz,
                                                   max_velocity_range=max_velocity_range,
                                                   alias_window=alias_window, ensemble=False,
                                                   loggers=self.loggers, quiet=self.quiet)
        ended_in_error = numpy.array([e is not None for e in err])
        if numpy.any(ended_in_error):
            if not self.quiet:
                warnings.warn('Masking failures in some/all spectra.  Errors are: {0}'.format(
                                numpy.array([(i,e) for i,e in enumerate(err)])[ended_in_error]))
            model_fit_par['MASK'][ended_in_error] \
                    = self.bitmask.turn_on(model_fit_par['MASK'][ended_in_error], 'NO_FIT')
            model_eml_par['MASK'][ended_in_error] \
                    = self.bitmask.turn_on(model_eml_par['MASK'][ended_in_error], 'NO_FIT')
        if numpy.all(ended_in_error):
            return self.obj_wave, model_flux, model_eml_flux, model_mask, model_eml_par

        #---------------------------------------------------------------
        # Get the input pixel shift between the object and template
        # wavelength vectors; interpretted by pPXF as a base velocity
        # shift between the two
        self.base_velocity = numpy.array([PPXFFit.ppxf_tpl_obj_voff(self.tpl_wave,
                                                            self.obj_wave[s:e], self.velscale,
                                                            velscale_ratio=self.velscale_ratio)
                                                for s,e in zip(self.spectrum_start,
                                                               self.spectrum_end)])

        #---------------------------------------------------------------
        # Fit all spectra
        t = time.clock()
#        warnings.warn('debugging!')
#        self.obj_to_fit[ numpy.arange(self.nobj)[self.obj_to_fit][2:] ] = False
        
#&&&&&#&&&&&#&&&&&#&&&&&#&&&&&#&&&&&#&&&&&
#&&&&& (FB) use _fit_all_spectra
#&&&&&#&&&&&#&&&&&#&&&&&#&&&&&#&&&&&#&&&&&
        result = self._fit_all_spectra(plot=plot)#, plot_file_root=plot_file_root)
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, 'Fits completed in {0:.4e} min.'.format(
                       (time.clock() - t)/60))

        #---------------------------------------------------------------
        # Save the results
#&&&&&#&&&&&#&&&&&#&&&&&#&&&&&#&&&&&#&&&&&
#&&&&& (FB) use _save_results function
#&&&&&#&&&&&#&&&&&#&&&&&#&&&&&#&&&&&#&&&&&
        model_flux, model_eml_flux, model_mask, model_fit_par, model_eml_par \
                = self._save_results(etpl, result, model_mask, model_fit_par, model_eml_par)

        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, 'Sasuke finished')

        return self.obj_wave, model_flux, model_eml_flux, model_mask, model_fit_par, model_eml_par


#%%%%

# (FB) to himself: DO NOT MODIFY ANYTHING BELOW THIS
class EmissionLineTemplates:
    r"""
    Construct a set of emission-line templates based on an emission-line
    database. Input to Sasuke

    Args:
        wave (array-like): A single wavelength vector with the
            wavelengths for the template spectra.
        sigma_inst (float,array-like): The single value or value as a
            function of wavelength for the instrumental dispersion to
            use for the template construction.
        log (bool): (**Optional**) Flag that the wavelengths have been
            sampled geometrically.
        emldb (:class:`mangadap.par.emissionlinedb.EmissionLineDB'): (**Optional**)
            Emission-line database that is parsed to setup which lines
            to include in the same template because they are modeled as
            having the same velocity, velocity dispersion and flux
            ratio.  If not provided, no templates are constructed in
            instantiation; to build the templates using an existing
            instantiation, use :func:`build_templates`.

    """
    def __init__(self, wave, sigma_inst, log=True, emldb=None, loggers=None, quiet=False):

        self.loggers=None
        self.quiet=None

        self.wave = numpy.asarray(wave)
        if len(self.wave.shape) != 1:
            raise ValueError('Provided wavelengths must be a single vector.')

        _sinst = numpy.full(self.wave.size, sigma_inst, dtype=float) \
                    if isinstance(sigma_inst, (int, float)) else numpy.asarray(sigma_inst)
        if _sinst.shape != self.wave.shape:
            raise ValueError('Provided sigma_inst must be a single number or a vector with the'
                             'same length as the wavelength vector.')
        self.sigma_inst = interpolate.interp1d(self.wave, _sinst, assume_sorted=True)

        self.dv = numpy.full(self.wave.size, spectrum_velocity_scale(wave), dtype=float) if log \
                    else astropy.constants.c.to('km/s').value*spectral_coordinate_step(wave)/wave

        self.emldb = None           # Original database
        self.ntpl = None            # Number of templates
        self.flux = None            # Template fluxes
        self.tpli = None            # Template associated with each emission line
        self.comp = None            # Kinematic component associated with each template
        self.vgrp = None            # Velocity group associated with each template
        self.sgrp = None            # Sigma group associated with each template
        self.eml_sigma_inst = None  # Instrumental dispersion at the center of each line

#        if database is given then load the emission line database and build the templates
        if emldb is not None:
            self.build_templates(emldb, loggers=loggers, quiet=quiet)


    def _tied_index(self, i, primary=False):
        """
        Return the index of the line to which this one is tied and it's
        primary line, which can be the same.
        """
        db_rows = numpy.arange(self.emldb.neml)
        indx = db_rows[self.emldb['index'] == int(self.emldb['mode'][i][1:])][0]
        if not primary:
            return indx
        max_iter = 100
        j = 0
        while self.emldb['mode'][indx] != 'f' and j < max_iter:
            indx = db_rows[self.emldb['index'] == int(self.emldb['mode'][indx][1:])][0]
            j+=1
        if j == max_iter:
            raise ValueError('Line {0} (index={1}) does not trace back to a primary line!'.format(
                                i, self.emldb['index'][i]))
        return indx


    def _parse_emission_line_database(self):
        r"""
        Parse the input emission-line database; see
        :class:`mangadap.par.emissionlinedb.EmissionLinePar`.

        Only lines with `action=f` are included in any template.  The
        array :attr:`tpli` provides the index of the template with each
        line in the emission-line database.  Lines that are not assigned
        to any template, either because they do not have `action=f` or
        their center lies outside the wavelength range in :attr:`wave`,
        are given an index of -1.

        Only lines with `mode=a` (i.e., tie flux, velocity, and velocity
        dispersion) are included in the same template.

        Lines with tied velocities are assigned the same velocity
        component (:attr:`vcomp`) and lines with the tied velocity
        dispersions are assigned the same sigma component
        (:attr:`scomp`).

        .. warning::
            The construction of templates for use with :class:`Sasuke`
            does *not* allow one to tie fluxes while leaving the
            velocities and/or velocity dispersions as independent.

        """
        # Get the list of lines to ignore
        ignore_line = self.emldb['action'] != 'f'

        # The total number of templates to construct is the number of
        # lines in the database minus the number of lines with mode=aN
        tied_all = numpy.array([m[0] == 'a' for m in self.emldb['mode']])
        self.ntpl = self.emldb.neml - numpy.sum(ignore_line) - numpy.sum(tied_all)

        # Initialize the components
        self.comp = numpy.zeros(self.ntpl, dtype=int)-1
        self.vgrp = numpy.zeros(self.ntpl, dtype=int)-1
        self.sgrp = numpy.zeros(self.ntpl, dtype=int)-1

        # All the primary lines go into individual templates, kinematic
        # components, velocity groups, and sigma groups
        self.tpli = numpy.zeros(self.emldb.neml, dtype=int)-1
        primary_line = (self.emldb['mode'] == 'f') & numpy.invert(ignore_line)
        nprimary = numpy.sum(primary_line)
        self.tpli[primary_line] = numpy.arange(nprimary)
        self.comp[:nprimary] = numpy.arange(nprimary)
        self.vgrp[:nprimary] = numpy.arange(nprimary)
        self.sgrp[:nprimary] = numpy.arange(nprimary)

        finished = primary_line | ignore_line
        while numpy.sum(finished) != self.emldb.neml:
            # Find the indices of lines that are tied to finished lines
            for i in range(self.emldb.neml):
                if finished[i]:
                    continue
                indx = self._tied_index(i)
                if not finished[indx]:
                    continue

                finished[i] = True

                # Mode=a: Line is part of an existing template
                if self.emldb['mode'][i][0] == 'a':
                    self.tpli[i] = self.tpli[indx]
                # Mode=k: Line is part of a different template but an
                # existing kinematic component
                if self.emldb['mode'][i][0] == 'k':
                    self.tpli[i] = numpy.amax(self.tpli)+1
                    self.comp[self.tpli[i]] = self.comp[self.tpli[indx]]
                    self.vgrp[self.tpli[i]] = self.vgrp[self.tpli[indx]]
                    self.sgrp[self.tpli[i]] = self.sgrp[self.tpli[indx]]
                # Mode=v: Line is part of a different template and
                # kinematic component with an untied sigma, but tied to
                # an existing velocity group
                if self.emldb['mode'][i][0] == 'v':
                    self.tpli[i] = numpy.amax(self.tpli)+1
                    self.comp[self.tpli[i]] = numpy.amax(self.comp)+1
                    self.sgrp[self.tpli[i]] = numpy.amax(self.sgrp)+1
                    self.vgrp[self.tpli[i]] = self.vgrp[self.tpli[indx]]
                # Mode=s: Line is part of a different template and
                # kinematic component with an untied velocity, but tied
                # to an existing sigma group
                if self.emldb['mode'][i][0] == 's':
                    self.tpli[i] = numpy.amax(self.tpli)+1
                    self.comp[self.tpli[i]] = numpy.amax(self.comp)+1
                    self.vgrp[self.tpli[i]] = numpy.amax(self.vgrp)+1
                    self.sgrp[self.tpli[i]] = self.sgrp[self.tpli[indx]]

        # Debug:
        if numpy.any(self.comp < 0) or numpy.any(self.vgrp < 0) or numpy.any(self.sgrp < 0):
            raise ValueError('DEBUG: Incorrect parsing of emission-line database.')


    def check_database(self, emldb):
        r"""
        Check that the provided emission-line database can be used with
        the :class:`EmissionLineTemplates` class.  Most checks are
        performed by
        :func:`mangadap.proc.spectralfitting.EmissionLineFit.check_emission_line_database`.

        Additional checks specific to :class:`EmissionLineTemplates`
        are:
            - Any lines with mode `w` is treated as `f` and a warning is
              provided.
            - The :class:`EmissionLineTemplates` object *cannot* be used
              with mode `x`; any lines with this mode will cause a
              ValueError to be raised..

        This function does *not* check if the initial parameters
        provided by the database are consistent with other elements in
        the database because they are not used to construct the
        templates.

        Args:
            emldb (:class:`mangadap.par.emissionlinedb.EmissionLineDB'):
                Emission-line database.

        Raises:
            TypeError: Raised if the provided object is not an instance
                of :class:`mangadap.par.emissionlinedb.EmissionLineDB`.
            ValueError: Raised if any line has a mode of `x` or if the
                database does not provide a valid definition for any
                templates.
            NameError: Raised if a defined profile type is not known.
        """
        EmissionLineFit.check_emission_line_database(emldb, wave=self.wave, check_par=False)

        # Check that no lines only tie the fluxes
        if numpy.any([m[0] == 'x' for m in emldb['mode']]):
            raise ValueError('Cannot tie only fluxes in an EmissionLineTemplates object.')

        # Warn user of any lines with mode=w
        if numpy.any([m[0] == 'w' for m in emldb['mode']]):
            warnings.warn('Any line with mode=w treated the same as mode=f.')


    def build_templates(self, emldb, loggers=None, quiet=False):
        r"""
        Build the set of templates for a given emission-line database.
        The function uses the current values in :attr:`wave` and
        :attr:`sigma_inst`.

        Warn the user if any line is undersampled; i.e., the FWHM of the
        line is less than 2.1 or sigma < 0.9.

        Warn the user if any line grouped in the same template falls
        outside the spectral range.

        Args:
            emldb (:class:`mangadap.par.emissionlinedb.EmissionLineDB'):
                Emission-line database.

        Returns:
            numpy.ndarray: Returns 4 arrays: (1) the set of templates
            with shape :math:`N_{\rm tpl}\times N_{\rm wave}`, (2) the
            kinematic component assignement for each template, (3) the
            velocity group associated with each template, and (4) the
            sigma group assocated with each template.
        """
        #---------------------------------------------------------------
        # Initialize the reporting
        if loggers is not None:
            self.loggers = loggers
        self.quiet = quiet

        #---------------------------------------------------------------
        # Check the database can be used with this class
        self.check_database(emldb)
        # Save a pointer to the database
        self.emldb = emldb
        # Parse the database for construction of the templates
        self._parse_emission_line_database()

        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO,
                       'Number of emission lines to fit: {0}'.format(numpy.sum(self.tpli>-1)))
            log_output(self.loggers, 1, logging.INFO,
                       'Number of emission-line templates: {0}'.format(len(self.comp)))
            log_output(self.loggers, 1, logging.INFO,
                       'Number of emission-line kinematic components: {0}'.format(
                                                                    numpy.amax(self.comp)+1))
            log_output(self.loggers, 1, logging.INFO,
                       'Number of emission-line velocity groups: {0}'.format(
                                                                    numpy.amax(self.vgrp)+1))
            log_output(self.loggers, 1, logging.INFO,
                       'Number of emission-line sigma groups: {0}'.format(
                                                                    numpy.amax(self.sgrp)+1))

        # Get the instrumental dispersion at the center of each line
        self.eml_sigma_inst = self.sigma_inst(self.emldb['restwave'])

        # Constuct the templates
        self.flux = numpy.zeros((self.ntpl,self.wave.size), dtype=float)
        for i in range(self.ntpl):
            # Find all the lines associated with this template:
            index = numpy.arange(self.emldb.neml)[self.tpli == i]
            # Add each line to the template
            for j in index:
                profile = eval('lineprofiles.'+self.emldb['profile'][j])()
                p = profile.parameters_from_moments(self.emldb['flux'][j], 0.0,
                                                    self.eml_sigma_inst[j])
                v = astropy.constants.c.to('km/s').value*(self.wave/self.emldb['restwave'][j]-1)
                srt = numpy.argsort(numpy.absolute(v))
                if self.eml_sigma_inst[j]/self.dv[srt[0]] < 0.9:
                    warnings.warn('{0} line is undersampled!'.format(self.emldb['name'][j]))
                self.flux[i,:] += profile(v, p)

        return self.flux, self.comp, self.vgrp, self.sgrp




