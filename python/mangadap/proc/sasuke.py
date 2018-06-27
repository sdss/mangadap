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
        from ..util.log import log_output
        from ..util.pixelmask import SpectralPixelMask
        from .spatiallybinnedspectra import SpatiallyBinnedSpectra
        from .stellarcontinuummodel import StellarContinuumModel
        from .spectralfitting import EmissionLineFit
        from .util import sample_growth

*Class usage examples*:
        Add examples

*Revision history*:
    | **24 May 2017**: Original implementation started by K. Westfall (KBW)
    | **23 Jun 2017**: (KBW) Documentation; fix error in
        :func:`Sasuke._save_results`
    | **08 Sep 2017**: (KBW) Move emission-line template construction to
        :mod:`mangadap.proc.emissionlinetemplates`. Add
        `deconstruct_bins` flag to parameters.
    | **25 Oct 2017**: (KBW) Allow to use reddening fit from pPXF
    | **02 Feb 2018**: (KBW) Use 
        :func:`mangadap.proc.spectralfitting.EmissionLineFit.select_binned_spectra_to_fit`.
        Adjust for change to
        :func:`mangadap.proc.stellarcontinuummodel.StellarContinuumModel.fill_to_match`.
    | **24 Feb 2018**: (KBW) Allow for a new template library to be
        used, different from the one used during the stellar continuum
        fit.
    | **20 Mar 2018**: (KBW) Corrected flux calculation and inclusion
        of provided pixel mask.
    | **22 May 2018**: (KBW) Change import to ppxf package.
    | **29 May 2018**: (KBW) Change xjmc function import and call

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
import logging

import numpy
from scipy import interpolate, fftpack

import astropy.constants
from ppxf import ppxf

from ..par.parset import ParSet
from ..par.emissionlinedb import EmissionLineDB
from ..util.fitsutil import DAPFitsUtil
from ..util.fileio import init_record_array
from ..util.filter import interpolate_masked_vector
from ..util.instrument import spectrum_velocity_scale, spectral_coordinate_step
from ..util.instrument import SpectralResolution
from ..util.log import log_output
from ..util.pixelmask import SpectralPixelMask
from ..util.constants import DAPConstants
from ..util import lineprofiles
from .templatelibrary import TemplateLibrary
from .spatiallybinnedspectra import SpatiallyBinnedSpectra
from .stellarcontinuummodel import StellarContinuumModel
from .spectralfitting import EmissionLineFit
from .bandpassfilter import emission_line_equivalent_width
from .emissionlinetemplates import EmissionLineTemplates
from .util import sample_growth
from .ppxffit import PPXFModel, PPXFFitResult, PPXFFit
from ..contrib.xjmc import emline_fitter_with_ppxf

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

    The continuum templates can either be the string keyword used to
    construct the template library, or the template library itself.  If
    the former, the fit_SpatiallyBinnedSpectra object will construct the
    template library for later callback.

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
        continuum_templates
            (str, :class:`mangadap.proc.templatelibrary.TemplateLibrary`):
            (**Optional**) The new continuum template library to use
            during the emission-line fit.
        etpl_line_sigma_mode (str): (**Optional**) Mode used to set the
            instrumental dispersion of the emission-line templates.
            Mode options are explated by
            :func:`Sasuke.etpl_line_sigma_options`.  Default is
            'default'.
        etpl_line_sigma_min (int,float): (**Optional**) Impose a minimum
            emission-line sigma by offsetting the nominal trend, in
            quadrature, to have this minimum value.  Default is 0.
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
        reddening (float): pPXF reddening parameter setting the initial
            :math:`E(B-V)` to fit, based on a Calzetti law.
    """
    def __init__(self, stellar_continuum, emission_lines, continuum_templates=None,
                 etpl_line_sigma_mode=None, etpl_line_sigma_min=None, velscale_ratio=None,
                 guess_redshift=None, guess_dispersion=None, minimum_snr=None,
                 deconstruct_bins=None, pixelmask=None, reject_boxcar=None, bias=None,
                 moments=None, degree=None, mdegree=None, reddening=None):

        arr_like = [ numpy.ndarray, list ]
        arr_in_fl = [ numpy.ndarray, list, int, float ]
        in_fl = [ int, float ]
        etpl_mode_options = Sasuke.etpl_line_sigma_options()

        pars =     [ 'stellar_continuum', 'emission_lines', 'continuum_templates',
                     'etpl_line_sigma_mode', 'etpl_line_sigma_min', 'velscale_ratio',
                     'guess_redshift', 'guess_dispersion', 'minimum_snr', 'deconstruct_bins',
                     'pixelmask', 'reject_boxcar', 'bias', 'moments', 'degree', 'mdegree',
                     'reddening' ]
        values =   [ stellar_continuum, emission_lines, continuum_templates, etpl_line_sigma_mode,
                     etpl_line_sigma_min, velscale_ratio, guess_redshift, guess_dispersion,
                     minimum_snr, deconstruct_bins, pixelmask, reject_boxcar, bias, moments,
                     degree, mdegree, reddening ]
        defaults = [ None, None, None, 'default', 0.0, 1, None, None, 0.0, False, None, None, None,
                     2, -1, 0, None ]
        options = [ None, None, None, etpl_mode_options, None, None, None, None, None, None, None,
                    None, None, None, None, None, None ]
        dtypes =   [ StellarContinuumModel, EmissionLineDB, [str, TemplateLibrary], str, in_fl,
                     int, arr_in_fl, arr_in_fl, in_fl, bool, SpectralPixelMask, int, in_fl, int,
                     int, int, in_fl ]

        ParSet.__init__(self, pars, values=values, defaults=defaults, options=options,
                        dtypes=dtypes)

        self._check()


    def _check(self):
        if self['mdegree'] > 0 and self['reddening'] is not None:
            raise ValueError('Cannot fit both multiplicative polynomials and an extinction curve.')


class Sasuke(EmissionLineFit):
    r"""
    Use ninja skills and pPXF to fit emission lines.

    https://en.wikipedia.org/wiki/Sasuke_Uchiha

    Effectively, **nothing** happens during the instantiation of this
    object.  Instead, a typical usage of the class when fitting a set of
    emission lines would be::

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
        reddening (float): pPXF reddening parameter setting the initial
            :math:`E(B-V)` to fit, based on a Calzetti law.
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

# REMOVED
#        tpl_sres (:class:`mangadap.util.instrument.SpectralResolution`):
#            Spectral resolution of the template spectra.  All templates
#            are assumed to have the same spectral resolution.

#        gas_comp (numpy.ndarray): Boolean array set to True for
#            emission-line components.  Shape is (:math:`N_{\rm comp}`).
#        self.gas_comp = None

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
        self.obj_skyx = None
        self.obj_skyy = None
        self.nobj = None
        self.npix_obj = None
        self.input_obj_mask = None
        self.obj_to_fit = None
        self.input_cz = None
        self.velscale = None
        self.waverange = None

        # Data to remap to
        self.remap_flux = None
        self.remap_ferr = None
        self.remap_sres = None
        self.nremap = 0
        self.remap_skyx = None
        self.remap_skyy = None

        # Template data
        self.tpl_wave = None
        self.tpl_flux = None
#        self.tpl_sres = None
        self.tpl_to_use = None
        self.nstpl = None
        self.ntpl = None
        self.npix_tpl = None
        self.tpl_npad = None
        self.tpl_rfft = None
        self.gas_tpl = None

        self.matched_resolution = None
        self.velscale_ratio = None

        self.emldb = None
        self.neml = None
        self.fit_eml = None
        self.eml_tpli = None
        self.eml_compi = None

        # Kinematic components and tied parameters
        self.ncomp = None
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
        self.reddening = None
        self.reject_boxcar = None
#        self.fix_kinematics = False

        self.spectrum_start = None
        self.spectrum_end = None
        self.dof = None
        self.base_velocity = None


    @staticmethod
    def _per_fit_dtype(ntpl, nadd, nmult, nkin, mask_dtype):
        r"""
        Returns the data type for the result of each pPXF fit.
        
        The data are as follows:
            - ``BINID``: Spectrum ID number.
            - ``BINID_INDEX``: Index of the spectrum in the list of
              provided spectra.
            - ``NEAREST_BIN``: ID of the nearest binned spectrum when
              deconstructing bins.
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
            - ``APLYMINMAX``: Minimum and maximum of additive
              polynomial.  Shape is (2,).
            - ``MULTCOEF``: Multiplicative polynomal coefficients.
              Shape is (:math:`o_{\rm mult}`,)
            - ``MPLYMINMAX``: Minimum and maximum of multiplicative
              polynomial.  Shape is (2,).
            - ``EBV``: Fitted E(B-V) from pPXF, if applied.
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

        In :func:`fit`, ``BINID`` and ``BINID_INDEX`` are the same.
        They are not in :func:`fit_SpatiallyBinnedSpectra`.

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
                 ('NEAREST_BIN',numpy.int),
                 ('MASK', mask_dtype),
                 ('BEGPIX', numpy.int),
                 ('ENDPIX', numpy.int),
                 ('NPIXTOT',numpy.int),
                 ('NPIXFIT',numpy.int),
                 ('KINCMP',numpy.int,(ntpl,)),
                 ('VELCMP',numpy.int,(ntpl,)),
                 ('SIGCMP',numpy.int,(ntpl,)),
#                 ('USETPL',numpy.bool,(ntpl,)),
                 ('TPLWGT',numpy.float,(ntpl,)),
                 ('TPLWGTERR',numpy.float,(ntpl,)),
                 ('ADDCOEF',numpy.float,(nadd,)) if nadd > 1 else ('ADDCOEF',numpy.float),
                 ('APLYMINMAX',numpy.float,(2,)) if nadd > 1 else ('APLYMINMAX',numpy.float),
                 ('MULTCOEF',numpy.float,(nmult,)) if nmult > 1 else ('MULTCOEF',numpy.float),
                 ('MPLYMINMAX',numpy.float,(2,)) if nmult > 1 else ('MPLYMINMAX',numpy.float),
                 ('EBV',numpy.float),
                 ('KININP',numpy.float,(nkin,)),
                 ('KIN',numpy.float,(nkin,)),
                 ('KINERR',numpy.float,(nkin,)),
                 ('TIEDKIN',numpy.int,(nkin,)),
                 ('CHI2',numpy.float),
                 ('RCHI2',numpy.float),
#                 ('ROBUST_RCHI2',numpy.float),
                 ('RMS',numpy.float),
                 ('ABSRESID',numpy.float,(5,)),
                 ('FRMS',numpy.float),
                 ('FABSRESID',numpy.float,(5,))
               ] 


    @staticmethod
    def _copy_from_stellar_continuum(stellar_continuum):
        return stellar_continuum.method['fitpar']['template_library'], \
                stellar_continuum.method['fitpar']['match_resolution'], \
                stellar_continuum.method['fitpar']['velscale_ratio']


    @staticmethod
    def etpl_line_sigma_options():
        r"""

        Keyword options for the mode used to set the standard deviation
        of the emission lines the the emission-line templates.  Possible
        modes are:

            ``default``: In order of precedence: 
                - Use the spectral resolution of the observed spectra, or
                - use the spectral resolution of the stellar templates, or
                - adopt a FWHM of 2 pixels and calculate the resolution
                  assuming a Gaussian line profile.

            ``zero``:  The width of all lines is set to 0.  The function
                from :mod:`mangadap.util.lineprofiles` used to construct
                the line *must* be able to produce the line if the
                standard deviation is 0!

            ``offset``: Apply a quadrature offset of the
                wavelength-dependent trend of the instrumental
                dispersion resulting from the ``default`` mode such that
                the minimum instrumental dispersion is 0.  Note that the
                minimum dispersion can be set to something other than 0
                using either the `etpl_line_sigma_min` parameter in
                :class:`SasukePar` or the keyword argument in
                :func:`Sasuke.fit`.  If the minimum dispersion is 0, the
                function from :mod:`mangadap.util.lineprofiles` used to
                construct the line *must* be able to produce a line with
                a standard deviation of 0!

        Returns:
            list: List of allowed options.
        """
        return ['default', 'zero', 'offset' ]


    def _check_remapping_coordinates(self, obj_skyx, obj_skyy, remap_skyx, remap_skyy):
        if any(x is None for x in [obj_skyx, obj_skyy, remap_skyx, remap_skyy]):
            raise ValueError('If remapping, must provide all on-sky positions.')
        if len(obj_skyx.shape) != 1:
            raise ValueError('Object spectra coordinates must be a vector.')
        if obj_skyx.size != self.nobj:
            raise ValueError('Object sky coordinates to not match the number of spectra.')
        if obj_skyx.shape != obj_skyy.shape:
            raise ValueError('Mismatch of object x and y coordinate shapes.')
        if len(remap_skyx.shape) != 1:
            raise ValueError('Remapping coordinates must be a vector.')
        if remap_skyx.size != self.nremap:
            raise ValueError('Remapping sky coordinates to not match the number of spectra.')
        if remap_skyx.shape != remap_skyy.shape:
            raise ValueError('Mismatch of the remapping x and y coordinate shapes.')
        return obj_skyx, obj_skyy, remap_skyx, remap_skyy


#    def _run_fit_iteration(self, obj_flux, obj_ferr, obj_to_fit, weight_errors=False,
#                           component_fits=False, plot=False):
#        r"""
#        Execute a specific fit iteration.  This is used by
#        :func:`_fit_all_spectra` to run through as set of fit
#        iterations.  All spectra in *obj_flux* are fit with a specific
#        set of ppxf parameters.  Additional iteration modes can be added
#        by including new arguments to the function.
#        
#        Args:
#            obj_flux (`numpy.ma.MaskedArray`_): Object spectra to fit
#                **for this iteration**.  Different from
#                :attr:`obj_flux`.
#            obj_ferr (`numpy.ma.MaskedArray`_): Error in spectra to fit
#                **for this iteration**.  Different from
#                :attr:`obj_flux`.
#            obj_to_fit (numpy.ndarray): Boolean array to fit each
#                spectrum.
#            weight_errors (bool): (**Optional**) Flag to calculate and
#                assign errors in the weights; see
#                :class:`mangadap.proc.ppxffit.PPXFFitResult`.
#            component_fits (bool): (**Optional**) Flag to construct the
#                optimal fit for each kinematic component.  This is used
#                to construct the best-fitting emission-line model in
#                :func:`_emission_line_only_model`; see
#                :class:`mangadap.proc.ppxffit.PPXFFitResult`.
#            plot (bool): (**Optional**) Passed to
#                :class:`mangadap.contrib.ppxf.ppxf` to construct a plot
#                showing the result **of each fit**.  Should only used
#                when debugging.
#
#        Returns:
#            numpy.ndarray : Array with :math:`N_{\rm spec}` instances of
#            :class:`mangadap.proc.ppxffit.PPXFFitResult`.
#
#        """
##        linear = fix_kinematics and mdegree < 1
#        linear = False
#
#        # Create the object to hold all the fits
#        result = numpy.empty(self.nobj, dtype=object)
#
#        # Fit each spectrum
#        for i in range(self.nobj):
#            print('Running pPXF fit on spectrum: {0}/{1}'.format(i+1,self.nobj), end='\r')
#            # Meant to ignore this spectrum
#            if not obj_to_fit[i]:
#                result[i] = None
#                continue
#
#            # Get the pixels to fit for this spectrum
#            gpm = numpy.where( ~ (obj_flux.mask[i, self.spectrum_start[i]:self.spectrum_end[i]]) )[0]
#
#            # Check if there is sufficient data for the fit
#            ntpl_to_use = numpy.sum(self.tpl_to_use[i,:])
#            if len(gpm) < self.dof+ntpl_to_use:
#                if not self.quiet:                                        
#                    warnings.warn('Insufficient data points ({0}) to fit spectrum {1}'
#                                  '(dof={2}).'.format(len(gpm), i+1, self.dof+ntpl_to_use))
#                # Declare an empty PPXFFitResult
#                result[i] = PPXFFitResult(self.degree, self.mdegree, self.spectrum_start[i],
#                                          self.spectrum_end[i], self.tpl_to_use[i,:],
#                                          None, self.ntpl)
#                continue
#
#            # Clear any existing plot data
#            if plot:
#                pyplot.clf()
#
#            # Run ppxf
#            ppxf_fit = ppxf(self.tpl_flux[self.tpl_to_use[i,:],:].T,
#                            obj_flux.data[i,self.spectrum_start[i]:self.spectrum_end[i]],
#                            obj_ferr.data[i,self.spectrum_start[i]:self.spectrum_end[i]],
#                            self.velscale, self.comp_start_kin[i].tolist(), bias=self.bias, 
#                            component=self.tpl_comp[self.tpl_to_use[i,:]], degree=self.degree,
#                            goodpixels=gpm, linear=linear, mdegree=self.mdegree,
#                            moments=self.comp_moments, plot=plot, quiet=(not plot),
#                            templates_rfft=self.tpl_rfft[self.tpl_to_use[i,:],:].T,
#                            tied=self.tied, velscale_ratio=self.velscale_ratio,
#                            vsyst=-self.base_velocity[i])
#            # Save result
#            result[i] = PPXFFitResult(self.degree, self.mdegree, self.spectrum_start[i],
#                                      self.spectrum_end[i], self.tpl_to_use[i,:],
#                                      ppxf_fit, self.ntpl, weight_errors=weight_errors,
#                                      component_fits=component_fits)
#
#            # TODO: check output
##            if result[i].kin[1] < 0:
##                result[i].kin[1] = numpy.absolute(result[i].kin[1]) #self.sigma_limits[0]
##                warnings.warn('pPXF gives negative dispersion! Change -{0:.4f} to {0:.4f}'.format(
##                                    result[i].kin[1]))
#                
#            if result[i].reached_maxiter() and not self.quiet:
#                warnings.warn('pPXF optimizer reached maximum number of iterations for spectrum '
#                              '{0}.'.format(i+1))
#            if plot:
#                pyplot.show()
#
#        print('Running pPXF fit on spectrum: {0}/{0}'.format(self.nobj))
#        return result


#    def _fit_all_spectra(self, plot=False): #, plot_file_root=None):
#        """
#        Fit all spectra contained by the class instance
#        (:attr:`obj_flux`) using an optimized set of iterations.  For
#        example, see
#        :func:`mangadap.proc.ppxffit.PPXFFit._fit_all_spectra()`.
#        
#        Args:
#            plot (bool): (**Optional**) Passed to
#                :class:`mangadap.contrib.ppxf.ppxf` to construct a plot
#                showing the result of **each fit** during **each fit
#                iteration**.  Should only used when debugging.
#        
#        Returns:
#            numpy.ndarray : Array with :math:`N_{\rm spec}` instances of
#            :class:`mangadap.proc.ppxffit.PPXFFitResult`.
#        
#        """
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
#        # Return results of refit (only ever do one rejection iteration)
#        return self._run_fit_iteration(obj_flux, obj_ferr, obj_to_fit, weight_errors=True,
#                                       component_fits=True, plot=plot)
#
##        if not self.quiet:
##            log_output(self.loggers, 1, logging.INFO,
##                       'Number of object spectra to fit: {0}/{1}'.format(
##                            numpy.sum(self.obj_to_fit), len(self.obj_to_fit)))
##        result = self._run_fit_iteration(self.obj_flux, self.obj_ferr, self.obj_to_fit,
##                                         weight_errors=True,
##                                         component_fits=True, plot=plot)
##        return result


#    def _emission_line_only_model(self, result):
#        r"""
#        Given the set of :class:`mangadap.proc.ppxf.PPXFFitResult`
#        instances resulting from the fit to all object spectra,
#        construct the best-fitting models that only include the model
#        emission lines.  The :class:`mangadap.proc.ppxf.PPXFFitResult`
#        instances must have had the best-fitting components constructed
#        when instantiated (declared with `component_fits=True`); the
#        'best-fitting components' are the templates weighted and
#        modified by any fitted multiplicative polynomial.  See
#        :func:`_run_fit_iteration` and :func:`_save_results`.
#
#        Args:
#            result (numpy.ndarray): An array of
#                :class:`mangadap.proc.ppxffit.PPXFFitResult` instances
#                with the result of the fits to each object spectrum.
#
#        Returns:
#            `numpy.ma.MaskedArray`_: A masked array with shape
#            :math:`(N_{\rm spec},N_{\rm pix})` with the best-fitting
#            emission-line-only model.  All pixels between
#            :attr:`spectrum_start` and :attr:`spectrum_end` are not
#            masked.
#        """
#        # Models originally fully masked
#        model_eml_flux = numpy.ma.MaskedArray(numpy.zeros(self.obj_flux.shape, dtype=float),
#                                              mask=numpy.ones(self.obj_flux.shape, dtype=bool))
#        for i in range(self.nobj):
#            if result[i] is None or result[i].fit_failed():
#                continue
#            s = result[i].start
#            e = result[i].end
#            # Sum all the emission-line components
#            model_eml_flux[i,s:e] = numpy.sum(result[i].bestfit_comp[self.gas_comp,:], axis=0)
#        return model_eml_flux


    def _is_near_bounds(self, kin, kin_inp, vel_indx, sig_indx, lbound, ubound, tol_frac=1e-2,
                        fill_value=-999.):
        r"""
        Check if the kinematics were fit and near the imposed limits.

        Any value that is close to ``fill_value`` is flagged as not
        having been fit, presumably because too much/all of the data
        near the emission line is masked.
        
        The definition of "near" is that the velocity and higher moments
        cannot be closer than the provided fraction of the total width
        to the boundary.  For the velocity dispersion, the fraction is
        done in log space.

        Args:
            kin (numpy.ndarray): Best-fitting kinematic parameters.
                Shape should be :math:`(N_{\rm spec},N_{\rm kin})`.
            kin_inp (numpy.ndarray): The initial guesses for the
                best-fitting kinematics.  This is needed because the
                velocity limits are set relative to the input guess.
                Shape should be :math:`(N_{\rm spec},N_{\rm kin})`.
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
            numpy.ndarray: Three boolean arrays that flag if (1) the
            parameter is near either boundary, (2) the parameter is near
            the lower boundary, and (3) the parameter is near the
            fill_value.  The second array is important in case the
            fitted parameter is a velocity dispersion and that it's near
            the lower boundary because it has hit the pixel sampling
            limit.

        """
        nspec = kin.shape[0]        # Number of spectra fit

        # Offset velocity: bounded by *deviations* from input value
        _lbound = numpy.array([lbound]*nspec)
        _ubound = numpy.array([ubound]*nspec)
        _lbound[:,vel_indx] += kin_inp[:,vel_indx]
        _ubound[:,vel_indx] += kin_inp[:,vel_indx]

        # Set the tolerance
        Db = _ubound-_lbound
        Db[:,sig_indx] = numpy.log10(_ubound[:,sig_indx])-numpy.log10(_lbound[:,sig_indx])
        tol = Db*tol_frac

        # Determine if the parameter is near the lower boundary (only
        # relevant to the sigma) ... 
        near_lower_bound = kin - _lbound < tol
        # and whether it's close to either
        near_bound = near_lower_bound | (_ubound - kin < tol)

        # Determine if the parameter is near the fill value indicating
        # that it was not included in the fit
        tol = 1e-5                      # !! HARDCODED !!
        no_data = kin - fill_value < tol

        # If the parameter was not fit, parameter not near any boundary
        near_bound[no_data] = False
        near_lower_bound[no_data] = False

#        for k,b,l,f in zip(kin[0], near_bound[0], near_lower_bound[0], no_data[0]):
#            print(k,b,l,f)

        # Return the two boundary flags
        return near_bound, near_lower_bound, no_data


    def _validate_dispersions(self, model_eml_par, rng=[0,400]):
        """
        Check that the corrected velocity dispersion are in the provided
        range.

        (FB) More flagging of kinematic paramters. DOES THIS APPLY TO
        EMISSION LINES? (KBW): Yes.  The corrected dispsersion must be
        larger than 0 and less than 400 km/s.  It's easy to turn this
        off or change the limits.

        (KBW) Currently not called...

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
        _fit_eml = numpy.zeros(model_eml_par['SIGMACORR'].shape, dtype=numpy.bool)
        _fit_eml[:,self.fit_eml] = True
        sigcor = numpy.square(model_eml_par['KIN'][:,:,1]) \
                        - numpy.square(model_eml_par['SIGMACORR'][:,:])
        indx = ((sigcor < _rng[0]) | (sigcor > _rng[1])) & _fit_eml
        if numpy.sum(indx) == 0:
            return
        model_eml_par['MASK'][indx] = self.bitmask.turn_on(model_eml_par['MASK'][indx], 'BAD_SIGMA')

        return model_eml_par


    def _reset_to_fill_value(self, model_eml_par, fill_value=-999.):
        """
        Reset all the emission lines masked as having insufficient data
        to fit to the provided fill value.

        .. note::
            The error values are currently reset to 0. by
            :class:`mangadap.proc.emissionlinemodel.EmissionLineModel.`
        """
        indx = self.bitmask.flagged(model_eml_par['MASK'], flag='INSUFFICIENT_DATA')
        model_eml_par['FLUX'][indx] = fill_value
        model_eml_par['FLUXERR'][indx] = fill_value
        model_eml_par['KIN'][indx,:] = fill_value
        model_eml_par['KINERR'][indx,:] = fill_value
        model_eml_par['SIGMACORR'][indx] = fill_value
        model_eml_par['SIGMAINST'][indx] = fill_value
        model_eml_par['SIGMATPL'][indx] = fill_value
        if self.degree > -1:
            model_eml_par['CONTAPLY'][indx] = fill_value
        if self.mdegree > 0:
            model_eml_par['CONTMPLY'][indx] = fill_value
        if self.reddening is not None:
            model_eml_par['CONTRFIT'][indx] = fill_value
        return model_eml_par


    def _save_results(self, etpl, start, end, flux, ferr, model_flux, model_eml_flux, model_mask,
                      model_wgts, model_wgts_err, model_addcoef, model_multcoef, model_reddening,
                      model_kin_inp, model_kin, model_kin_err, model_fit_par, model_eml_par,
                      fill_value=-999.):
        r"""
        Save and assess the results of the ppxf fits.
        
        The results are saved both as the direct output from pPXF and
        after parsing the data into the results for each emission line.
        The function also assesses the data to set any necessary flags.
        Much of this is the same as
        :class:`mangadap.proc.ppxffit.PPXFFit._save_results`.

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
        # Generate some convenience data:
        #  - Generate vectors with the lower and upper bounds for the
        #    kinematic parameters (lbound, ubound)
        #  - Get the list of indices in the flatted kinematics vectors
        #    with the *unfixed*, *defining* parameters used for each
        #    kinematic measurement.  These are used to set the
        #    kinematics and errors for each emission line. (par_indx)
        #  - Flag parameters that are velocity (vel_indx) and sigma
        #    components (sig_indx)
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
            ii = numpy.sum(numpy.absolute(self.comp_moments[:j]))
            nmom = numpy.absolute(self.comp_moments[j])
            par_indx += [ [0]*nmom ]
            for k in range(nmom):
                par_indx[j][k] = ii+k if len(self.tied[j][k]) == 0 \
                                            else int(self.tied[j][k].split('[')[1].split(']')[0])
            vel_indx[ii+0] = True
            sig_indx[ii+1] = True
            lbound += [ lboundi[:nmom] ]
            ubound += [ uboundi[:nmom] ]
        lbound = numpy.concatenate(tuple(lbound))
        ubound = numpy.concatenate(tuple(ubound))

#        print(self.tpl_comp)
#        print(self.tpl_vgrp)
#        print(self.tpl_sgrp)
#        print(lbound)
#        print(ubound)
#        print(numpy.concatenate(tuple(par_indx)))

        # Set the model data to masked arrays
        model_flux = numpy.ma.MaskedArray(model_flux.data, mask=model_mask > 0)
        model_eml_flux = numpy.ma.MaskedArray(model_eml_flux.data, mask=model_mask > 0)

#        pyplot.imshow(numpy.log10(model_flux), origin='lower', interpolation='nearest', aspect='auto')
#        pyplot.colorbar()
#        pyplot.show()

        # Save the pixel statistics
        model_fit_par['BEGPIX'][:] = start
        model_fit_par['ENDPIX'][:] = end
        model_fit_par['NPIXTOT'][:] = end - start
        model_fit_par['NPIXFIT'] = numpy.sum(numpy.invert(model_flux.mask[:,start:end]), axis=1)

        # Calculate the model residuals, which are masked where the data
        # were not fit
        residual = flux - model_flux
        fractional_residual = numpy.ma.divide(flux - model_flux, model_flux)
        # Get the (reduced) chi-square for each spectrum
        model_fit_par['CHI2'] = numpy.ma.sum(numpy.square(residual/ferr), axis=1)
        # Get the (fractional) residual RMS for each spectrum
        model_fit_par['RMS'] = numpy.sqrt(numpy.ma.mean(numpy.square(residual), axis=1))
        model_fit_par['FRMS'] = numpy.sqrt(numpy.ma.mean(numpy.square(fractional_residual), axis=1))

#        print(model_fit_par['CHI2'])

        # Reduced chi-square
        model_fit_par['RCHI2'] = model_fit_par['CHI2'] / (model_fit_par['NPIXFIT'] 
                                    - self.dof - numpy.sum(model_fit_par['TPLWGT'] > 0, axis=1))
#        print(model_fit_par['RCHI2'])

        # Save the weights and errors
        model_fit_par['TPLWGT'] = model_wgts
        model_fit_par['TPLWGTERR'] = model_wgts_err

        # The set of gas templates, and the kinematic component,
        # velocity group, and sigma group associated with each
        # template are the same for all fits
        nspec = self.nobj if self.nremap == 0 else self.nremap
        model_fit_par['KINCMP'] = numpy.array([self.tpl_comp]*nspec)
        model_fit_par['VELCMP'] = numpy.array([self.tpl_vgrp]*nspec)
        model_fit_par['SIGCMP'] = numpy.array([self.tpl_sgrp]*nspec)
        model_fit_par['TIEDKIN'] = numpy.array([numpy.concatenate(tuple(par_indx))]*nspec)

        # Save the polynomial coefficients
        used_apoly = False
        if self.degree > -1 and model_addcoef is not None:
            model_fit_par['ADDCOEF'] = model_addcoef
            used_apoly = True

        used_mpoly = False
        if self.mdegree > 0 and model_multcoef is not None:
            model_fit_par['MULTCOEF'] = model_multcoef
            used_mpoly = True

        used_ebv = False
        if self.reddening is not None and model_reddening is not None:
            model_fit_par['EBV'] = model_reddening
            used_ebv = True

        # Flattened input and output kinematics
        model_fit_par['KININP'] = model_kin_inp
        model_fit_par['KIN'] = model_kin
        model_fit_par['KINERR'] = model_kin_err

        # Test if the kinematics are near the imposed boundaries.
        near_bound, near_lower_bound, no_data = self._is_near_bounds(model_fit_par['KIN'],
                                                                     model_fit_par['KININP'],
                                                                     vel_indx, sig_indx, lbound,
                                                                     ubound, fill_value=fill_value)

        # Flag the fit *globally*
        # - If the velocity dispersion has hit the lower limit for all
        #   lines, ONLY flag the value as having a MIN_SIGMA.
        indx = numpy.all(near_lower_bound & sig_indx[None,:], axis=1)
        if numpy.sum(indx) > 0:
            model_fit_par['MASK'][indx] = self.bitmask.turn_on(model_fit_par['MASK'][indx],
                                                               'MIN_SIGMA')

        # - Otherwise, flag the full fit (parameters and model) as
        #   NEAR_BOUND if all the parameters are near the boundary but
        #   not the lower sigma boundary
        indx = numpy.all( (near_lower_bound & numpy.invert(sig_indx)[None,:])
                            | (near_bound & numpy.invert(near_lower_bound)), axis=1)
        if numpy.sum(indx) > 0:
            model_fit_par['MASK'][indx] = self.bitmask.turn_on(model_fit_par['MASK'][indx],
                                                               'NEAR_BOUND')
            model_mask[indx,:] = self.bitmask.turn_on(model_mask[indx,:], 'NEAR_BOUND')

        # - If no_data is all true for any one spectrum, the ppxf fitter
        # should have raised an error.

        # Convert the velocities from pixel units to cz
        model_fit_par['KININP'][:,vel_indx], _ \
                    = PPXFFit.convert_velocity(model_fit_par['KININP'][:,vel_indx],
                                               numpy.zeros(numpy.sum(vel_indx)))
        model_fit_par['KIN'][:,vel_indx], model_fit_par['KINERR'][:,vel_indx] \
                    = PPXFFit.convert_velocity(model_fit_par['KIN'][:,vel_indx],
                                               model_fit_par['KINERR'][:,vel_indx])

        # Divvy up the fitted parameters into the result for each
        # emission line
        for j in range(self.neml):
            if not self.fit_eml[j]:
                continue
#            print(self.emldb['name'][j], self.emldb['restwave'][j], self.emldb['flux'][j],
#                  self.eml_tpli[j], result[i].tplwgt[self.eml_tpli[j]])
#            pyplot.plot(self.tpl_wave, self.tpl_flux[self.eml_tpli[j],:])
#            pyplot.show()

            # The "fit index" is the component of the line
            model_eml_par['FIT_INDEX'][:,j] = self.eml_compi[j]

            # Use the flattened vectors to set the kinematics
            indx = par_indx[self.eml_compi[j]]
            model_eml_par['KIN'][:,j,:] = model_fit_par['KIN'][:,indx]
            model_eml_par['KINERR'][:,j,:] = model_fit_par['KINERR'][:,indx]
            z = model_eml_par['KIN'][:,j,0]/astropy.constants.c.to('km/s').value

            # Use the fitted weights to set the gas flux; the
            # EmissionLineTemplates class constructs each line to have
            # the flux provided by the emission-line database in the
            # rest frame.  The ppxf convolution keeps the sum of the
            # template constant, meaning that the total flux in the line
            # increases with redshift.
            model_eml_par['FLUX'][:,j] = model_fit_par['TPLWGT'][:,self.eml_tpli[j]] \
                                            * self.emldb['flux'][j] * (1 + z)
            model_eml_par['FLUXERR'][:,j] = model_fit_par['TPLWGTERR'][:,self.eml_tpli[j]] \
                                                * self.emldb['flux'][j] * (1 + z)

            # Get the bound masks specific to this emission-line (set)
            # - Determine if the emission-line was part of a rejected
            # template
            flg = numpy.any(no_data[:,indx], axis=1)
            if numpy.sum(flg) > 0:
                model_eml_par['MASK'][flg,j] = self.bitmask.turn_on(model_eml_par['MASK'][flg,j],
                                                                    'INSUFFICIENT_DATA')

            # - Determine if the velocity dispersion parameter of this
            #   line has hit the lower limit; if so, ONLY flag the value
            #   as having a MIN_SIGMA.
            flg = numpy.any(near_lower_bound[:,indx] & sig_indx[None,indx], axis=1)
            if numpy.sum(flg) > 0:
                model_eml_par['MASK'][flg,j] = self.bitmask.turn_on(model_eml_par['MASK'][flg,j],
                                                                    'MIN_SIGMA')

            # - Determine if any of the kinematic parameters are near
            #   the bound (excluding the lower velocity dispersion
            #   limit)
            flg = numpy.any((near_lower_bound[:,indx] & numpy.invert(sig_indx)[None,indx])
                                | (near_bound[:,indx] & numpy.invert(near_lower_bound[:,indx])),
                            axis=1)
            if numpy.sum(flg) > 0:
                model_eml_par['MASK'][flg,j] = self.bitmask.turn_on(model_eml_par['MASK'][flg,j],
                                                                    'NEAR_BOUND')

        # Flag the pixels that were not used
        model_mask[flux.mask] = self.bitmask.turn_on(model_mask[flux.mask], flag='DIDNOTUSE')

        # Mask any lines that were not fit
        model_eml_par['MASK'][:,numpy.invert(self.fit_eml)] \
                    = self.bitmask.turn_on(model_eml_par['MASK'][:,numpy.invert(self.fit_eml)],
                                           flag='NO_FIT')

        #---------------------------------------------------------------
        # Iterate over each spectrum
        poly_x = numpy.linspace(-1, 1, end-start)
        sres = self.obj_sres if self.nremap == 0 else self.remap_sres
        for i in range(nspec):

            # Get growth statistics for the residuals
            # TODO: Specifically around the emission lines?
            model_fit_par['ABSRESID'][i,:] = sample_growth(numpy.ma.absolute(residual[i,:]),
                                                           [0.0, 0.68, 0.95, 0.99, 1.0])
            model_fit_par['FABSRESID'][i,:] \
                        = sample_growth(numpy.ma.absolute(fractional_residual[i,:]),
                                        [0.0, 0.68, 0.95, 0.99, 1.0])

            if used_apoly or used_mpoly or used_ebv:
                # Sample the polynomials at the fitted wavelength of
                # each line
                obswave = self.emldb['restwave'][self.fit_eml] \
                                    * (1 + model_eml_par['KIN'][i,self.fit_eml,0]
                                                / astropy.constants.c.to('km/s').value)

                # - Additive polynomial:
                if used_apoly:
                    apoly = numpy.polynomial.legendre.legval(poly_x,
                                                             numpy.append(1, model_addcoef[i,:]))
                    model_fit_par['APLYMINMAX'][i,:] = [numpy.amin(apoly), numpy.amax(apoly)]
                    model_eml_par['CONTAPLY'][i,:] \
                            = interpolate.interp1d(self.obj_wave[start:end], apoly,
                                                   bounds_error=False, fill_value=0.0,
                                                   assume_sorted=True)(obswave)

                # - Multiplicative polynomial:
                if used_mpoly:
                    mpoly = numpy.polynomial.legendre.legval(poly_x,
                                                             numpy.append(1, model_multcoef[i,:]))
                    model_fit_par['MPLYMINMAX'][i,:] = [numpy.amin(mpoly), numpy.amax(mpoly)]
                    model_eml_par['CONTMPLY'][i,:] \
                            = interpolate.interp1d(self.obj_wave[start:end], mpoly,
                                                   bounds_error=False, fill_value=1.0,
                                                   assume_sorted=True)(obswave)
                #-------------------------------------------------------
                # As of version 6.6.5, pPXF no longer applies the
                # multiplicative polynomial to the gas templates
                # factor = interpolate.interp1d(self.obj_wave[start:end], mpoly, bounds_error=False,
                #                              fill_value=1.0, assume_sorted=True)(obswave)
                # model_eml_par['FLUX'][i,self.fit_eml] *= factor
                # model_eml_par['FLUXERR'][i,self.fit_eml] *= factor
                #-------------------------------------------------------

                # - Reddening:
                if used_ebv:
                    extcurve = ppxf.reddening_cal00(self.obj_wave, model_reddening[i])
                    model_eml_par['CONTRFIT'][i,:] \
                            = interpolate.interp1d(self.obj_wave, extcurve,
                                                   bounds_error=False, fill_value=1.0,
                                                   assume_sorted=True)(obswave)

            # Get the instrumental dispersion in the galaxy data at the
            # location of the fitted lines
            # TODO: sres has to be provided!
            model_eml_par['SIGMAINST'][i,self.fit_eml] = \
                        EmissionLineFit.instrumental_dispersion(self.obj_wave, sres[i,:],
                                                        self.emldb['restwave'][self.fit_eml],
                                                        model_eml_par['KIN'][i,self.fit_eml,0])

            # Set the instrumental dispersion of the emission line
            # templates to the output database
            model_eml_par['SIGMATPL'][i,self.fit_eml] = etpl.eml_sigma_inst[self.fit_eml]

            # Add the template dispersion into the fitted dispersion to
            # get the observed dispersion
            model_eml_par['KINERR'][i,self.fit_eml,1] \
                        = model_eml_par['KIN'][i,self.fit_eml,1] \
                            * model_eml_par['KINERR'][i,self.fit_eml,1]
            model_eml_par['KIN'][i,self.fit_eml,1] \
                        = numpy.sqrt( numpy.square(model_eml_par['KIN'][i,self.fit_eml,1])
                                        + numpy.square(model_eml_par['SIGMATPL'][i,self.fit_eml]))
            model_eml_par['KINERR'][i,self.fit_eml,1] /= model_eml_par['KIN'][i,self.fit_eml,1]

            # With these definitions, the sigma correction is the same
            # as the instrumental dispersion; see copy function outside
            # the loop below

#            # The dispersion correction is the quadrature difference
#            # between the instrumental dispersion in the galaxy data to
#            # the dispersion used when constructing the emission-line
#            # templates
#            sigma2corr = numpy.square(sigma_inst) - numpy.square(etpl.eml_sigma_inst[self.fit_eml])
#            if numpy.any(sigma2corr < 0):
#                # TODO: Raise an error instead?
#                print(sigma2corr)
#                warnings.warn('Encountered negative sigma corrections!')
#            model_eml_par['SIGMACORR'][i,self.fit_eml] \
#                            = sigma2corr / numpy.sqrt(numpy.absolute(sigma2corr))

#            print(model_eml_par['SIGMACORR'][i,self.fit_eml])
#            pyplot.scatter(self.emldb['restwave'][self.fit_eml], sigma_inst,
#                           marker='.', lw=0, s=100)
#            pyplot.scatter(self.emldb['restwave'][self.fit_eml], etpl.eml_sigma_inst[self.fit_eml],
#                           marker='.', lw=0, s=100)
#            pyplot.show()

            #-----------------------------------------------------------
            # Set output flags
            # - No fit was performed
#            if result[i] is None:
#                model_mask[i,:] = self.bitmask.turn_on(model_mask[i,:], 'NO_FIT')
#                continue

            # - No fit attempted because of insufficient data
#            if result[i].empty_fit():
#                model_mask[i,:] = self.bitmask.turn_on(model_mask[i,:], 'NO_FIT')
#                model_fit_par['MASK'][i] = self.bitmask.turn_on(model_fit_par['MASK'][i],
#                                                                'INSUFFICIENT_DATA')
#                model_eml_par['MASK'][i] = self.bitmask.turn_on(model_eml_par['MASK'][i],
#                                                                'INSUFFICIENT_DATA')
#                continue

            # - Fit attempted but failed
#            if result[i].fit_failed():
#                model_mask[i,:] = self.bitmask.turn_on(model_mask[i,:], 'FIT_FAILED')
#                model_fit_par['MASK'][i] = self.bitmask.turn_on(model_fit_par['MASK'][i],
#                                                                'FIT_FAILED')
#                model_eml_par['MASK'][i] = self.bitmask.turn_on(model_eml_par['MASK'][i],
#                                                                'FIT_FAILED')

            # - Fit successful but hit maximum iterations.
#            if result[i].reached_maxiter():
#                model_fit_par['MASK'][i] = self.bitmask.turn_on(model_fit_par['MASK'][i], 'MAXITER')
#                model_eml_par['MASK'][i] = self.bitmask.turn_on(model_eml_par['MASK'][i], 'MAXITER')

            #-----------------------------------------------------------
            # Save the model parameters and figures of merit

            # TODO: FIX THIS
#            model_fit_par['USETPL'][i,:] = result[i].tpl_to_use

        # With the above definitions (starting at line 1274), the
        # instrumental sigma and the sigma correction are identical
        model_eml_par['SIGMACORR'] = model_eml_par['SIGMAINST'].copy()

        #---------------------------------------------------------------
        # Reset any parameters based on insufficient data to the fill_value
        model_eml_par = self._reset_to_fill_value(model_eml_par, fill_value=fill_value)

        #---------------------------------------------------------------
        # Test if kinematics are reliable
#        model_eml_par = self._validate_dispersions(model_eml_par)

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
        # TODO: May also want to exclude pixels rejected during stellar
        # kinematics fit
        wave, flux, ferr = EmissionLineFit.get_spectra_to_fit(binned_spectra,
                                                              pixelmask=par['pixelmask'],
                                                              error=True)
        sres = binned_spectra.copy_to_array(ext='SPECRES')
        _sres = numpy.ma.MaskedArray(sres.copy(), mask=numpy.invert(sres > 0))
        for i in range(sres.shape[0]):
            sres[i,:] = interpolate_masked_vector(_sres[i,:])

        # Get the binned spectra that meet the S/N criterion
        bins_to_fit = EmissionLineFit.select_binned_spectra_to_fit(binned_spectra,
                                                                   minimum_snr=par['minimum_snr'],
                                                        stellar_continuum=par['stellar_continuum'])
        nobj = flux.shape[0]

        # Get the stellar templates.  If fitting a different set of
        # templates, the spectral resolution of the new templates must
        # be matched to the galaxy data and the velocity dispersion used
        # by the fit must be astrophysical (corrected for any resolution
        # difference).
        stellar_templates = None
        matched_resolution = False
        velscale_ratio = 1
        if par['continuum_templates'] is not None:
            if isinstance(par['continuum_templates'], TemplateLibrary):
                # The template library has already been instantiated so
                # just copy the object.  Assume this means that the
                # spectral resolution has been matched to the MaNGA
                # data because otherwise par['continuum_templates'] is
                # None.
                stellar_templates = par['continuum_templates']
                # TODO: Instead test that the spectral resolution of the
                # template libary has been matched to the MaNGA data?
                matched_resolution = True
                velscale_ratio = par['velscale_ratio']
            else:
                # Otherwise it must be the keyword of the library that
                # needs to be constructed.

                if par['stellar_continuum'] is not None \
                        and par['continuum_templates'] \
                            == par['stellar_continuum'].method['fitpar']['template_library_key']:
                    # The templates used are the same, so warn the user
                    # and just copy over the existing template library.
                    # This maintains the resolution difference
                    warnings.warn('Request emission-line continuum templates identical to those '
                                  'used during the stellar continuum fitting.')
                    stellar_templates, matched_resolution, velscale_ratio \
                                = Sasuke._copy_from_stellar_continuum(par['stellar_continuum'])
#                    stellar_templates \
#                            = par['stellar_continuum'].method['fitpar']['template_library']
#                    matched_resolution \
#                            = par['stellar_continuum'].method['fitpar']['match_resolution']
#                    velscale_ratio \
#                            = par['stellar_continuum'].method['fitpar']['velscale_ratio']
                else:
                    # The template library needs to be constructed based
                    # on the provided keyword
                    # TODO: The TemplateLibrary object uses the median
                    # spectral resolution vector when performing the
                    # resolution match to the MaNGA data.
                    match_resolution = True
                    velscale_ratio = par['velscale_ratio']
                    velocity_offset = numpy.mean(astropy.constants.c.to('km/s').value*par['guess_redshift'][bins_to_fit,:])
                    # TODO: The dapsrc and analysis_path will return to
                    # the defaults!
                    stellar_templates \
                            = TemplateLibrary(par['continuum_templates'],
                                              velocity_offset=velocity_offset,
                                              drpf=binned_spectra.drpf,
                                              match_to_drp_resolution=True,
                                              velscale_ratio=velscale_ratio, hardcopy=False,
                                              loggers=loggers, quiet=quiet)

        elif par['continuum_templates'] is None and par['stellar_continuum'] is not None:
            stellar_templates, matched_resolution, velscale_ratio \
                        = Sasuke._copy_from_stellar_continuum(par['stellar_continuum'])

        stpl_wave = None if stellar_templates is None else stellar_templates['WAVE'].data
        stpl_flux = None if stellar_templates is None else stellar_templates['FLUX'].data
        if not quiet:
            warnings.warn('Adopting mean spectral resolution of all templates!')
        stpl_sres = None if stellar_templates is None \
                        else numpy.mean(stellar_templates['SPECRES'].data, axis=0).ravel()

        # Get the stellar kinematics
        stellar_velocity, stellar_dispersion = (None, None) if par['stellar_continuum'] is None \
                        else par['stellar_continuum'].matched_kinematics(
                                                binned_spectra['BINID'].data, cz=True,
                                                nearest=True, missing=binned_spectra.missing_bins,
                                                corrected=matched_resolution)
        stellar_kinematics = None if stellar_velocity is None or stellar_dispersion is None \
                                else numpy.array([ stellar_velocity, stellar_dispersion ]).T

        # Set which stellar templates to use for each spectrum
        # TODO: Template modes:
        #   - provide all
        #   - provide non-zero
        #   - optimal template from each bin
#        stpl_to_use = None if par['stellar_continuum'] is None \
#                        else par['stellar_continuum'].matched_template_flags(binned_spectra)
        stpl_to_use = None                  # Use all templates

        # TODO: For now can only fit two moments
        if par['moments'] != 2:
            print(par['moments'])
            raise NotImplementedError('Number of gas moments can only be two.')

        if par['deconstruct_bins']:
            # TODO: Only need to refit "bins" that are made up of more
            # than one spaxel...

            # Get the individual spaxel data
            _, spaxel_flux, spaxel_ferr = EmissionLineFit.get_spectra_to_fit(binned_spectra.drpf,
                                                                        pixelmask=par['pixelmask'],
                                                                          error=True)
            # Apply the Galactic extinction
            spaxel_flux, spaxel_ferr = binned_spectra.galext.apply(spaxel_flux, err=spaxel_ferr,
                                                                   deredden=True)

            # TODO: set minimum_fraction as a keyword.  Set to 0.8 by
            # default
            spaxel_to_fit = binned_spectra.check_fgoodpix()

            # Get the spectral resolution:
            # - stack_sres sets whether or not the spectral resolution
            #   is determined on a per-spaxel basis or with a single
            #   vector
            stack_sres = binned_spectra.method['stackpar']['stack_sres']
            # - prepixel_sres sets if the prepixelized version of the
            #   LSF measurements were used
            prepixel_sres = binned_spectra.method['prepixel_sres']
            # This pulls out the appropiate spectral resolution
            spaxel_sres = binned_spectra.drpf.spectral_resolution(
                                            ext=None if stack_sres else 'SPECRES', toarray=True,
                                            fill=True, pre=prepixel_sres)
            # Spaxel coordinates; shape is (Nspaxel,)
            spaxel_x = binned_spectra.rdxqa['SPECTRUM'].data['SKY_COO'][:,0]
            spaxel_y = binned_spectra.rdxqa['SPECTRUM'].data['SKY_COO'][:,1]
            # Get the binned coordinates
            bin_x = binned_spectra['BINS'].data['SKY_COO'][:,0]
            bin_y = binned_spectra['BINS'].data['SKY_COO'][:,1]

            # Return the fits to the individual spaxel data
            model_wave, model_flux, model_eml_flux, model_mask, model_fit_par, model_eml_par \
                    = self.fit(par['emission_lines'], wave, flux[bins_to_fit,:],
                               obj_ferr=ferr[bins_to_fit,:], obj_sres=sres[bins_to_fit,:],
                               guess_redshift=par['guess_redshift'][bins_to_fit],
                               guess_dispersion=par['guess_dispersion'][bins_to_fit],
                               reject_boxcar=par['reject_boxcar'], stpl_wave=stpl_wave,
                               stpl_flux=stpl_flux, stpl_sres=stpl_sres,
                               stpl_to_use=None if stpl_to_use is None \
                                                else stpl_to_use[bins_to_fit,:],
                               stellar_kinematics=None if stellar_kinematics is None \
                                                       else stellar_kinematics[bins_to_fit,:],
                               etpl_sinst_mode=par['etpl_line_sigma_mode'],
                               etpl_sinst_min=par['etpl_line_sigma_min'],
                               remap_flux=spaxel_flux[spaxel_to_fit,:],
                               remap_ferr=spaxel_ferr[spaxel_to_fit,:],
                               remap_sres=spaxel_sres[spaxel_to_fit,:],
                               remap_skyx=spaxel_x[spaxel_to_fit],
                               remap_skyy=spaxel_y[spaxel_to_fit], obj_skyx=bin_x[bins_to_fit],
                               obj_skyy=bin_y[bins_to_fit], velscale_ratio=velscale_ratio,
                               matched_resolution=matched_resolution, bias=par['bias'],
#                               degree=-1, mdegree=0,
                               degree=par['degree'], mdegree=par['mdegree'],
                               reddening=par['reddening'],
                               #moments=par['moments'],
                               loggers=loggers, quiet=quiet)

            # Convert the index number of the nearest bin to the BIN ID
            # number
            model_fit_par['NEAREST_BIN'] \
                = binned_spectra['BINS'].data['BINID'][bins_to_fit][model_fit_par['NEAREST_BIN']]

            model_binid = numpy.full(binned_spectra.spatial_shape, -1, dtype=int)
            model_binid.ravel()[spaxel_to_fit] = numpy.arange(numpy.sum(spaxel_to_fit))

            # Add the equivalent width data
#            redshift = model_eml_par['KIN'][:,0]/astropy.constants.c.to('km/s').value
            # TODO: Redshift is set to fitted redshift.  Should it use
            # this or the input guess?
            EmissionLineFit.measure_equivalent_width(wave, spaxel_flux[spaxel_to_fit,:],
                                                     par['emission_lines'], model_eml_par,
                                                     #redshift=redshift,
                                                     bitmask=self.bitmask, checkdb=False)
        else:
            # Return the fits to the binned data
            model_wave, model_flux, model_eml_flux, model_mask, model_fit_par, model_eml_par \
                    = self.fit(par['emission_lines'], wave, flux[bins_to_fit,:],
                               obj_ferr=ferr[bins_to_fit,:], obj_sres=sres[bins_to_fit,:],
                               guess_redshift=par['guess_redshift'][bins_to_fit],
                               guess_dispersion=par['guess_dispersion'][bins_to_fit],
                               reject_boxcar=par['reject_boxcar'], stpl_wave=stpl_wave,
                               stpl_flux=stpl_flux, stpl_sres=stpl_sres,
                               stpl_to_use=None if stpl_to_use is None \
                                                else stpl_to_use[bins_to_fit,:],
                               stellar_kinematics=None if stellar_kinematics is None \
                                                       else stellar_kinematics[bins_to_fit,:],
                               etpl_sinst_mode=par['etpl_line_sigma_mode'],
                               etpl_sinst_min=par['etpl_line_sigma_min'],
                               velscale_ratio=velscale_ratio,
                               matched_resolution=matched_resolution, bias=par['bias'],
#                               degree=-1, mdegree=0,
                               degree=par['degree'], mdegree=par['mdegree'],
                               reddening=par['reddening'],
                               #moments=par['moments'],
                               loggers=loggers, quiet=quiet)

            # Save the the bin ID numbers indices based on the spectra
            # selected to be fit
            model_fit_par['BINID'] = binned_spectra['BINS'].data['BINID'][bins_to_fit]
            model_fit_par['BINID_INDEX'] = numpy.arange(binned_spectra.nbins)[bins_to_fit]
            model_fit_par['NEAREST_BIN'] = binned_spectra['BINS'].data['BINID'][bins_to_fit]

            model_eml_par['BINID'] = binned_spectra['BINS'].data['BINID'][bins_to_fit]
            model_eml_par['BINID_INDEX'] = numpy.arange(binned_spectra.nbins)[bins_to_fit]

            # Bin IDs are the same as on input
            model_binid = None

            # Add the equivalent width data
#            redshift=par['guess_redshift'][bins_to_fit],
            # TODO: Redshift is set to fitted redshift.  Should it use
            # this or the input guess?  Should the EWs be measured in
            # EmissionLineModel instead?  (Computation not specific to
            # Sasuke...)
            EmissionLineFit.measure_equivalent_width(wave, flux[bins_to_fit,:],
                                                     par['emission_lines'], model_eml_par,
                                                     #redshift=redshift,
                                                     bitmask=self.bitmask, checkdb=False)

        # Reset the equivalent widths to the fill value.  The error
        # values are currently reset to 0. by EmissionLineModel.
        indx = self.bitmask.flagged(model_eml_par['MASK'], flag='INSUFFICIENT_DATA')
        model_eml_par['BMED'][indx] = -999.
        model_eml_par['RMED'][indx] = -999.
        model_eml_par['EWCONT'][indx] = -999.
        model_eml_par['EW'][indx] = -999.
        model_eml_par['EWERR'][indx] = -999.

        # Calculate the "emission-line baseline" as the difference
        # between the stellar continuum model determined for the
        # kinematics and the one determined by the optimized
        # stellar-continuum + emission-line fit:
        if par['stellar_continuum'] is not None:

            if par['deconstruct_bins']:

                # Get the stellar continuum.  The model extends over
                # over regions masked during the fit, but 0 outside the
                # spectral range of the fit.
                sc_continuum = par['stellar_continuum'].fill_to_match(
                                                    binned_spectra['BINID'].data,
                                                    missing=binned_spectra.missing_bins)

#                print(type(sc_continuum))
#                for _i in range(10):
#                    pyplot.plot(binned_spectra['WAVE'].data, sc_continuum.data[_i,:])
#                    pyplot.plot(binned_spectra['WAVE'].data, sc_continuum[_i,:])
#                    pyplot.show()

                # Construct the full 3D cube for the stellar continuum
                sc_model_flux = DAPFitsUtil.reconstruct_cube(binned_spectra.drpf.shape,
                                                             binned_spectra['BINID'].data.ravel(),
                                                             sc_continuum.filled(0.0))

#                print(type(sc_model_flux))
#                print(sc_model_flux.shape)
#                nx = sc_model_flux.shape[0]
#                for _i in range(10):
##                    pyplot.plot(binned_spectra['WAVE'].data, sc_model_flux.data[_i,_i,:])
#                    pyplot.plot(binned_spectra['WAVE'].data, sc_model_flux[nx//2+_i-5,nx//2+_i-5,:])
#                    pyplot.show()

                # Construct the full 3D cube of the new stellar
                # continuum from the combined stellar-continuum +
                # emission-line fit
                el_continuum = DAPFitsUtil.reconstruct_cube(binned_spectra.drpf.shape,
                                                            model_binid.ravel(),
                                                            model_flux - model_eml_flux)

#                print(type(el_continuum))
#                print(el_continuum.shape)
#                nx = el_continuum.shape[0]
#                for _i in range(10):
##                    pyplot.plot(binned_spectra['WAVE'].data, el_continuum.data[_i,_i,:])
#                    pyplot.plot(binned_spectra['WAVE'].data, el_continuum[nx//2+_i-5,nx//2+_i-5,:])
#                    pyplot.show()

                # Get the difference, restructure it to match the shape
                # of the emission-line models
                model_eml_base = (el_continuum - sc_model_flux).reshape(-1,
                                                                self.npix_obj)[spaxel_to_fit,:]
#                if model_mask is not None:
#                    model_eml_base[model_mask==0] = 0.0

#                print(type(model_eml_base))
#                print(model_eml_base.shape)
#                for _i in range(10):
#                    pyplot.plot(binned_spectra['WAVE'].data, model_eml_base[_i,:])
#                    pyplot.show()

            else:
                el_continuum = model_flux - model_eml_flux
                sc_continuum = par['stellar_continuum'].fill_to_match(binned_spectra['BINID'].data,
                                                            missing=binned_spectra.missing_bins)
                model_eml_base = (el_continuum - sc_continuum[bins_to_fit,:]).filled(0.0)
        else:
            model_eml_base = numpy.zeros(model_flux.shape, dtype=float)

#        pyplot.plot(model_wave, flux[spec_to_fit,:][0,:], color='k', lw=1, zorder=1)
#        pyplot.plot(model_wave, stellar_continuum[spec_to_fit,:][0,:], color='C0', lw=1.0,
#                    zorder=2)
#        pyplot.plot(model_wave, model_eml_flux[0,:], color='C1', lw=1.0, zorder=3)
#        pyplot.plot(model_wave, model_eml_base[0,:], color='C2', lw=1.0, zorder=4)
#        pyplot.plot(model_wave, eml_continuum[0,:], color='C3', lw=1.0, zorder=4)
#        pyplot.show()

#        print(par['emission_lines']['name'])
#        print(model_eml_par['FLUX'][0,:])
#        print(model_eml_par['FLUXERR'][0,:])
        
        # Only return model and model parameters for the *fitted*
        # spectra
#        print(type(model_eml_flux))
#        print(type(model_eml_base))
#        print(type(model_mask))
#        for _i in range(10):
#            pyplot.plot(binned_spectra['WAVE'].data, model_eml_flux[_i,:])
#            pyplot.plot(binned_spectra['WAVE'].data, model_eml_base[_i,:])
#            pyplot.plot(binned_spectra['WAVE'].data, model_mask[_i,:])
#            pyplot.show()
        return model_eml_flux, model_eml_base, model_mask, model_fit_par, model_eml_par,\
                    model_binid


    def fit(self, emission_lines, obj_wave, obj_flux, obj_ferr=None, obj_mask=None, obj_sres=None,
            guess_redshift=None, guess_dispersion=None, reject_boxcar=None, stpl_wave=None,
            stpl_flux=None, stpl_sres=None, stpl_to_use=None, stellar_kinematics=None,
            etpl_sinst_mode='default', etpl_sinst_min=0, remap_flux=None, remap_ferr=None,
            remap_mask=None, remap_sres=None, remap_skyx=None, remap_skyy=None, obj_skyx=None,
            obj_skyy=None, velscale_ratio=None, matched_resolution=True, waverange=None, bias=None,
            degree=-1, mdegree=0, reddening=None, max_velocity_range=400., alias_window=None,
            dvtol=1e-10, loggers=None, quiet=False, plot=False, sigma_rej=3.):
            #moments=2,
        r"""
        Fit a set of emission lines using pPXF to all provided spectra.

        The input flux arrays are expected to be ordered with spectra
        along **rows**. I.e., to plot the first spectrum in the array,
        one would do::
            
            from matplotlib import pyplot
            pyplot.plot(obj_wave, obj_flux[0,:])
            pyplot.show()

        The function will fit the spectra with or without any provided
        stellar templates.

        The body of this function mostly deals with checking the input
        and setting up the template and object data for use with pPXF.

        The function is meant to be general, but has only been tested on
        MaNGA spectra.

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
            - Allow to ignore emission-line database tying and just tie
              all kinematics
            - Should probably put the emission-line construction in a
              separate function.

        Args:
            emission_lines
                (:class:`mangadap.par.emissionlinedb.EmissionLineDB`):
                Emission-line database that is parsed to construct the
                emission-line templates to fit (see
                :class:`mangadap.proc.emissionelinetemplates.EmissionLineTemplates`).
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
            obj_mask (numpy.ndarray or
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
                kinematics to use for the stellar component.  If the
                spectral resolution of the templates is different from
                the galaxy data, the provided stellar velocity
                dispersions *must* include the the (assumed
                wavelength-independent) quadrature difference in the
                template and galaxy instrumental resolutions.  The
                stellar kinematics are **fixed** for all calls made to
                ppxf.  Shape is (:math:`N_{\rm spec},N_{{\rm
                kin},\ast}`).  The shape of this array determines the
                number of moments assigned to the stellar component.
            etpl_sinst_mode (str): (**Optional**) Mode used to set the
                instrumental dispersion of the emission-line templates;
                see :func:`etpl_line_sigma_options` for the options.
                Default mode is `default`, imagine that.
            etpl_sinst_min (str): (**Optional**) Minimum allowed
                instrumental dispersion value.  If the mode of
                constructing the instrumental dispersion of the
                templates results in values that are below this value,
                the whole function is offset in quadrature such that no
                value goes below the minimum.  Default is 0 km/s.
            remap_flux (numpy.ndarray): (**Optional**) Ultimately the
                spectra to be modeled, if provided.  The nominal usage
                is, e.g., if the object spectra (`obj_flux`) are binned
                spaxels, this can be flux arrays of all the spaxels that
                were used to construct the binned data.  The function
                would then use a first fit to the binned spectra and a
                matching between the on-sky positions (see `remap_skyx`,
                `remap_skyy`, `obj_skyx`, `obj_skyy`) to ultimately fit
                the individual spaxel data.  See
                :func:`mangadap.contrib.xjmc.emline_fitter_with_ppxf`.
                These spectra must be sampled identically to the object
                spectra (`velscale`) and have the same wavelength vector
                (`obj_wave`).  Can be provided as a
                `numpy.ma.MaskedArray`_.  Shape is (:math:`N_{\rm
                remap},N_{\rm pix}`).
            remap_ferr (numpy.ndarray): (**Optional**)  The 1-sigma
                error in `remap_flux`.  Can be provided as a
                `numpy.ma.MaskedArray`_.  Shape is (:math:`N_{\rm
                remap},N_{\rm pix}`).
            remap_mask (numpy.ndarray or
                :class:`mangadap.util.pixelmask.SpectralPixelMask`):
                (**Optional**) Boolean array or
                :class:`mangadap.util.pixelmask.SpectralPixelMask`
                instance used to censor regions of the remapping spectra
                to ignore during fitting.
            remap_sres (numpy.ndarray): (**Optional**)  The spectral
                resolution of the remapping spectra.  Can be a single
                vector for all spectra or one vector per object
                spectrum.
            remap_skyx (numpy.ndarray): (**Optional**) On-sky x position
                of the remapping spectra.  Shape is (:math:`N_{\rm
                remap}`,).
            remap_skyy (numpy.ndarray): (**Optional**) On-sky y position
                of the remapping spectra.  Shape is (:math:`N_{\rm
                remap}`,).
            obj_skyx (numpy.ndarray): (**Optional**) On-sky x position
                of the object spectra.  Shape is (:math:`N_{\rm
                spec}`,).
            obj_skyy (numpy.ndarray): (**Optional**) On-sky y position
                of the object spectra.  Shape is (:math:`N_{\rm
                spec}`,).
            velscale_ratio (int): (**Optional**) The **integer** ratio
                between the velocity scale of the pixel in the galaxy
                data to that of the template data.  If None, set to
                unity.
            matched_resolution (bool): (**Optional**) The spectral
                resolution of the templates is matched to that of the
                galaxy data.  WARNING: This functionality is never used!
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
                :math:`o_{\rm add}`.  Default is no polynomial included.
            mdegree (int): (**Optional**) pPXF mdegree parameter setting
                the degree of the multiplicative polynomials to use,
                :math:`o_{\rm mult}`.  Default is no polynomial
                included.
            reddening (float): (**Optional**) pPXF reddening parameter
                used to fit :math:`E(B-V)` assuming a Calzetti law.
                Cannot be fit simultaneously with multiplicative
                polynomial.  Must be larger than 0 to start.  Default is
                not fit.
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

        # The minimum instrumental dispersion of the emission-line
        # templates cannot be below 0.
        _etpl_sinst_min = etpl_sinst_min
        if _etpl_sinst_min < 0:
            warnings.warn('Emission-line template instrumental dispersion must be 0 km/s or '
                          'higher.  Setting etpl_sinst_min=0.')
            _etpl_sinst_min = 0

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
                            numpy.sum(self.obj_to_fit), self.nobj))
        # guess_kin are in "ppxf" velocities
        self.input_cz, guess_kin = PPXFFit.check_input_kinematics(self.nobj, guess_redshift,
                                                                  guess_dispersion)
#        print(self.input_cz[0], guess_kin[0,:])
#        print(numpy.mean(self.input_cz))

        # Check the remapping data, if provided
        if remap_flux is not None:
            _, self.remap_flux, self.remap_ferr, self.remap_sres \
                        = PPXFFit.check_objects(obj_wave, remap_flux, obj_ferr=remap_ferr,
                                                obj_sres=remap_sres)
            self.nremap = self.remap_flux.shape[0]
            self.obj_skyx, self.obj_skyy, self.remap_skyx, self.remap_skyy \
                        = self._check_remapping_coordinates(obj_skyx, obj_skyy,
                                                            remap_skyx, remap_skyy)
        else:
            self.nremap = 0
            self.remap_flux = None
            self.remap_ferr = None
            self.remap_sres = None

            self.obj_skyx = None
            self.obj_skyy = None
            self.remap_skyx = None
            self.remap_skyy = None

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
        R_to_sinst = astropy.constants.c.to('km/s').value / DAPConstants.sig2fwhm
        if stpl_flux is not None:
            if stpl_wave is None:
                raise ValueError('Must provide wavelengths if providing stellar template fluxes.')
            self.tpl_wave, self.tpl_flux, tpl_sres \
                    = PPXFFit.check_templates(stpl_wave, stpl_flux, tpl_sres=stpl_sres,
                                              velscale_ratio=self.velscale_ratio)
            self.nstpl = self.tpl_flux.shape[0]
            # Check or instantiate the fit flags
            self.tpl_to_use = PPXFFit.check_template_usage_flags(self.nobj, self.nstpl, stpl_to_use)

            # etpl_sinst defaults to the instrumental dispersion of the
            # stellar templates
            etpl_sinst = R_to_sinst/tpl_sres.sres(copy=True)
        else:
            self.tpl_flux = None
            etpl_sinst = None
            self.nstpl = 0
            self.tpl_to_use = None

        # TODO: Allow the gas template resolution to be some fraction of
        # the stellar template resolution?

        #---------------------------------------------------------------
        # Set the emission-line spectral resolution.
        # - If the object resolution and template resolution are not
        #   available, use the 2-pixel resolution
        # - If the object resolution is not available, set to the
        #   stellar template resolution (i.e., use what is set above)
        # - If the object resolution is available, set it to the minimum
        #   object resolution at each wavelength channel.
        # - Force the instrumental dispersion to 0 if requested by the
        #   input mode.
        # - Apply a quardrature offset if the dispersion dips below the
        #   imposed minimum or if the user requested to offset to the
        #   minimum.
        self.matched_resolution = matched_resolution
        if etpl_sinst is None and self.obj_sres is None:
            self.matched_resolution = False
            etpl_sinst = numpy.full(self.npix_obj, 2*self.velscale/DAPConstants.sig2fwhm,
                                    dtype=float)
        if self.obj_sres is not None:
            interp = interpolate.interp1d(self.obj_wave,
                                          R_to_sinst/numpy.amin(self.obj_sres, axis=0),
                                          assume_sorted=True, bounds_error=False,
                                          fill_value='extrapolate')
            etpl_sinst = interp(self.tpl_wave * (1 + numpy.median(self.input_cz) 
                                                        / astropy.constants.c.to('km/s').value))
#            pyplot.plot(self.obj_wave, interp.y)
#            pyplot.plot(self.tpl_wave, etpl_sinst)
#            _etpl_sinst = etpl_sinst.copy()

        if etpl_sinst_mode == 'zero':
            etpl_sinst[:] = 0.
        min_sinst = numpy.amin(etpl_sinst)
        if (etpl_sinst_mode == 'offset' and min_sinst > _etpl_sinst_min) \
                or min_sinst < _etpl_sinst_min:
            dsigma_inst = numpy.square(min_sinst)-numpy.square(_etpl_sinst_min)
            etpl_sinst = numpy.sqrt(numpy.square(etpl_sinst) - dsigma_inst)
            min_sinst = numpy.amin(etpl_sinst)
#            print(min_sinst)
#        pyplot.plot(self.tpl_wave, etpl_sinst)
#        pyplot.plot(self.tpl_wave, numpy.sqrt(numpy.square(_etpl_sinst) - numpy.square(etpl_sinst)))
#        pyplot.show()
#        exit()

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
        etpl = EmissionLineTemplates(self.tpl_wave, etpl_sinst, emldb=self.emldb,
                                     loggers=self.loggers, quiet=self.quiet)
#        pyplot.plot(self.tpl_wave, numpy.sum(etpl.flux, axis=0))
#        for i in range(etpl.ntpl):
#            pyplot.plot(self.tpl_wave, etpl.flux[i,:])
#        pyplot.show()

        # Report the resolution mode
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO,
                       'Emission-line resolution mode: {0}'.format(etpl_sinst_mode))
            log_output(self.loggers, 1, logging.INFO,
                       'Imposed minimum instrumental dispersion: {0}'.format(etpl_sinst_min))

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
        self.degree = max(degree,-1)
        self.mdegree = max(mdegree,0)
        self.reddening = reddening
        moments = 2                #numpy.absolute(moments)
        self.reject_boxcar = reject_boxcar
#        self.fix_kinematics = False     #moments < 0

        #---------------------------------------------------------------
        # Compile the template fluxes, components, and velocity and
        # sigma groups; the gas templates are always appended after the
        # stellar templates
        if self.nstpl == 0:
            self.gas_tpl = numpy.ones(etpl.ntpl, dtype=bool)
            self.tpl_flux = etpl.flux
            self.tpl_to_use = numpy.ones((self.nobj,etpl.ntpl), dtype=numpy.bool)
            self.tpl_comp = etpl.comp
            self.tpl_vgrp = etpl.vgrp
            self.tpl_sgrp = etpl.sgrp
            self.ncomp = numpy.amax(self.tpl_comp)+1
            self.comp_moments = numpy.array([moments]*self.ncomp)
            self.comp_start_kin = numpy.array([[gk.tolist()]*self.ncomp for gk in guess_kin ])
        else:
            self.gas_tpl = numpy.append(numpy.zeros(self.tpl_flux.shape[0]),
                                        numpy.ones(etpl.ntpl)).astype(bool)
            self.tpl_flux = numpy.append(self.tpl_flux, etpl.flux, axis=0)
            self.tpl_to_use = numpy.append(self.tpl_to_use,
                                           numpy.ones((self.nobj,etpl.ntpl), dtype=numpy.bool),
                                           axis=1)
            self.tpl_comp = numpy.append(numpy.zeros(self.nstpl, dtype=int), etpl.comp+1)
            self.tpl_vgrp = numpy.append(numpy.zeros(self.nstpl, dtype=int), etpl.vgrp+1)
            self.tpl_sgrp = numpy.append(numpy.zeros(self.nstpl, dtype=int), etpl.sgrp+1)
            self.eml_tpli[self.fit_eml] += self.nstpl
            self.eml_compi[self.fit_eml] += 1
            self.ncomp = numpy.amax(self.tpl_comp)+1
            self.comp_moments = numpy.array([-stellar_moments] + [moments]*(self.ncomp-1))
            self.comp_start_kin = numpy.array([ [sk.tolist()] + [gk.tolist()]*(self.ncomp-1) 
                                            for sk,gk in zip(_stellar_kinematics, guess_kin) ])
        self.ntpl, self.npix_tpl = self.tpl_flux.shape
        self.tpl_npad = fftpack.next_fast_len(self.npix_tpl)
        self.tpl_rfft = numpy.fft.rfft(self.tpl_flux, self.tpl_npad, axis=1)

#        # Set which components are gas components
#        self.gas_comp = numpy.ones(self.ncomp, dtype=bool)
#        if self.nstpl > 0:
#            self.gas_comp[0] = False

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
        if not self.quiet:
            if self.nstpl > 0:
                log_output(self.loggers, 1, logging.INFO,
                           'Number of stellar templates: {0}'.format(self.nstpl))
            log_output(self.loggers, 1, logging.INFO, 'Pixel scale: {0} km/s'.format(self.velscale))
            log_output(self.loggers, 1, logging.INFO, 'Pixel scale ratio: {0}'.format(
                                                                            self.velscale_ratio))
            log_output(self.loggers, 1, logging.INFO, 'Dispersion limits: {0} - {1}'.format(
                                                                            *self.sigma_limits))
            log_output(self.loggers, 1, logging.INFO, 'Model degrees of freedom: {0}'.format(
                                                                            self.dof+self.ntpl))
            log_output(self.loggers, 1, logging.INFO, 'Number of tied parameters: {0}'.format(
                            0 if self.tied is None else
                            self.nfree_kin + numpy.sum(self.comp_moments[self.comp_moments < 0])))
            if self.remap_flux is not None:
                log_output(self.loggers, 1, logging.INFO,
                           'Emission-line fits remapped to a set of {0} spectra.'.format(
                            self.nremap))

        #---------------------------------------------------------------
        # Initialize the output arrays.  This is done here for many of
        # these objects only in case PPXFFit.initialize_pixels_to_fit()
        # fails and the code returns without running the pPXF fitting.
        #  - Set the output shape and initial mask
        output_shape = self.obj_flux.shape if self.nremap == 0 else self.remap_flux.shape
        #  - Model flux and error (begins fully masked)
        model_flux = numpy.ma.masked_all(output_shape, dtype=float)
        model_flux.data[:,:] = 0.0
        model_eml_flux = numpy.ma.masked_all(output_shape, dtype=float)
        model_eml_flux.data[:,:] = 0.0
        #  - Model parameters and fit quality
        model_fit_par = init_record_array(output_shape[0],
                                          self._per_fit_dtype(self.ntpl, self.degree+1,
                                          self.mdegree, self.npar_kin,
                                          self.bitmask.minimum_dtype()))
        #  - Model emission-line parameters
        model_eml_par = init_record_array(output_shape[0],
                                          self._per_emission_line_dtype(self.neml, 2,
                                                                self.bitmask.minimum_dtype()))
        model_eml_par['CONTMPLY'] = numpy.ones(model_eml_par['CONTMPLY'].shape, dtype=float)

        #---------------------------------------------------------------
        # Initialize the mask and the spectral range to fit for the
        # input object spectra
        # TODO: This alters the mask of self.obj_flux!!
        obj_model_mask, err, obj_start, obj_end \
                    = PPXFFit.initialize_pixels_to_fit(self.tpl_wave, self.obj_wave, self.obj_flux,
                                                       self.obj_ferr, self.velscale,
                                                       velscale_ratio=self.velscale_ratio,
                                                       waverange=waverange, mask=obj_mask,
                                                       bitmask=self.bitmask,
                                                       velocity_offset=self.input_cz,
                                                       max_velocity_range=max_velocity_range,
                                                       alias_window=alias_window,
                                                       ensemble=True, #False,
                                                       loggers=self.loggers, quiet=self.quiet)
        ended_in_error = numpy.array([e is not None for e in err])
        if numpy.any(ended_in_error):
            if not self.quiet:
                warnings.warn('Masking failures in some/all spectra.  Errors are: {0}'.format(
                                numpy.array([(i,e) for i,e in enumerate(err)])[ended_in_error]))
            if self.nremap == 0:
                model_fit_par['MASK'][ended_in_error] \
                        = self.bitmask.turn_on(model_fit_par['MASK'][ended_in_error], 'NO_FIT')
                model_eml_par['MASK'][ended_in_error] \
                        = self.bitmask.turn_on(model_eml_par['MASK'][ended_in_error], 'NO_FIT')
        if numpy.all(ended_in_error):
            if self.nremap == 0:
                return self.obj_wave, model_flux, model_eml_flux, obj_model_mask, \
                            model_fit_par, model_eml_par
            else:
                raise ValueError('Could not prepare spectra for fitting!')

        #---------------------------------------------------------------
        # Initialize the mask and the spectral range to fit for the
        # remapping spectra
        # TODO: This alters the mask of self.remap_flux!!
        if self.nremap != 0:
            remap_model_mask, err, remap_start, remap_end \
                        = PPXFFit.initialize_pixels_to_fit(self.tpl_wave, self.obj_wave,
                                                           self.remap_flux, self.remap_ferr,
                                                           self.velscale,
                                                           velscale_ratio=self.velscale_ratio,
                                                           waverange=waverange, # mask=remap_mask,
                                                           bitmask=self.bitmask,
                                                        velocity_offset=numpy.median(self.input_cz),
                                                           max_velocity_range=max_velocity_range,
                                                           alias_window=alias_window,
                                                           ensemble=True, #False,
                                                           loggers=self.loggers, quiet=self.quiet)
            ended_in_error = numpy.array([e is not None for e in err])
            if numpy.any(ended_in_error):
                if not self.quiet:
                    warnings.warn('Masking failures in some/all remapping spectra.  '
                                  'Errors are: {0}'.format(
                                    numpy.array([(i,e) for i,e in enumerate(err)])[ended_in_error]))
                model_fit_par['MASK'][ended_in_error] \
                        = self.bitmask.turn_on(model_fit_par['MASK'][ended_in_error], 'NO_FIT')
                model_eml_par['MASK'][ended_in_error] \
                        = self.bitmask.turn_on(model_eml_par['MASK'][ended_in_error], 'NO_FIT')
            if numpy.all(ended_in_error):
                return self.obj_wave, model_flux, model_eml_flux, remap_model_mask, \
                            model_fit_par, model_eml_par
        else:
            remap_model_mask= None
            remap_start = None
            remap_end = None

        #---------------------------------------------------------------
        # Get the input pixel shift between the object and template
        # wavelength vectors; interpretted by pPXF as a base velocity
        # shift between the two.  This is a single number, used for both
        # the object spectra and the remapping spectra if provided
        start = numpy.amax(obj_start if self.nremap == 0 else numpy.append(obj_start, remap_start))
        end = numpy.amin(obj_end if self.nremap == 0 else numpy.append(obj_end, remap_end))

        model_mask = obj_model_mask if self.nremap == 0 else remap_model_mask
        self.base_velocity = -PPXFFit.ppxf_tpl_obj_voff(self.tpl_wave, self.obj_wave[start:end],
                                                        self.velscale,
                                                        velscale_ratio=self.velscale_ratio)

        #---------------------------------------------------------------
        # Call the emission line fitter contributed by Xihan Ji and
        # Michele Cappellari
        t = time.clock()
#        warnings.warn('debugging!')
#        self.obj_to_fit[ numpy.arange(self.nobj)[self.obj_to_fit][2:] ] = False

        # Prep:
        wave = self.obj_wave[start:end]
        if self.nremap == 0:
            flux = self.obj_flux.data[:,start:end]
            ferr = self.obj_ferr.data[:,start:end]
            mask = numpy.invert(self.obj_flux.mask[:,start:end])
            flux_binned = None
            ferr_binned = None
            mask_binned = None
        else:
            flux = self.remap_flux.data[:,start:end]
            ferr = self.remap_ferr.data[:,start:end]
            mask = numpy.invert(self.remap_flux.mask[:,start:end])
            flux_binned = self.obj_flux.data[:,start:end]
            ferr_binned = self.obj_ferr.data[:,start:end]
            mask_binned = numpy.invert(self.obj_flux.mask[:,start:end])

        # Run the fitter
        model_flux[:,start:end], model_eml_flux[:,start:end], _model_mask, model_wgts, \
                model_wgts_err, model_addcoef, model_multcoef, model_reddening, model_kin_inp, \
                model_kin, model_kin_err, nearest_bin \
                            = emline_fitter_with_ppxf(self.tpl_flux, wave, flux, ferr, mask,
                                                      self.velscale, self.velscale_ratio,
                                                      self.tpl_comp, self.gas_tpl,
                                                      self.comp_moments, self.comp_start_kin,
                                                      tied=self.tied, degree=self.degree,
                                                      mdegree=self.mdegree,
                                                      reddening=self.reddening,
                                                      reject_boxcar=self.reject_boxcar,
                                                      vsyst=self.base_velocity,
                                                      tpl_to_use=self.tpl_to_use,
                                                      flux_binned=flux_binned,
                                                      noise_binned=ferr_binned,
                                                      mask_binned=mask_binned,
                                                      x_binned=self.obj_skyx,
                                                      y_binned=self.obj_skyy, x=self.remap_skyx,
                                                      y=self.remap_skyy, plot=plot, quiet=not plot,
                                                      sigma_rej=sigma_rej)
#                                                           plot=True, quiet=False)

        # Construct the bin ID numbers
        model_fit_par['BINID'] = numpy.arange(output_shape[0])
        model_fit_par['BINID_INDEX'] = numpy.arange(output_shape[0])
        model_fit_par['NEAREST_BIN'] = nearest_bin

        model_eml_par['BINID'] = numpy.arange(output_shape[0])
        model_eml_par['BINID_INDEX'] = numpy.arange(output_shape[0])

#        result = self._fit_all_spectra(plot=plot)#, plot_file_root=plot_file_root)
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, 'Fits completed in {0:.4e} min.'.format(
                       (time.clock() - t)/60))

        # Flag pixels as rejected during fitting; _model_mask is True
        # where the pixels were fit
        indx = mask & numpy.invert(_model_mask)
        model_mask[:,start:end][indx] = self.bitmask.turn_on(model_mask[:,start:end][indx],
                                                             PPXFFit.rej_flag)

        #---------------------------------------------------------------
        # Save the results
#        model_flux, model_eml_flux, model_mask, model_fit_par, model_eml_par \
#                = self._save_results(etpl, result, model_mask, model_fit_par, model_eml_par)

#        pyplot.imshow(numpy.log10(model_mask), origin='lower', interpolation='nearest',
#                      aspect='auto')
#        pyplot.show()

        if self.nremap == 0:
            flux = self.obj_flux
            ferr = self.obj_ferr
        else:
            flux = self.remap_flux
            ferr = self.remap_ferr
        model_flux, model_eml_flux, model_mask, model_fit_par, model_eml_par \
                = self._save_results(etpl, start, end, flux, ferr, model_flux, model_eml_flux,
                                     model_mask, model_wgts, model_wgts_err, model_addcoef,
                                     model_multcoef, model_reddening, model_kin_inp, model_kin,
                                     model_kin_err, model_fit_par, model_eml_par)

        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, 'Sasuke finished')

        return self.obj_wave, model_flux, model_eml_flux, model_mask, model_fit_par, model_eml_par


    @staticmethod
    def construct_continuum_models(emission_lines, stpl_wave, stpl_flux, obj_wave, obj_flux_shape,
                                   model_fit_par, select=None, redshift_only=False,
                                   deredshift=False, dispersion_corrections=None, dvtol=1e-10):
        """
        Construct the continuum models using the provided set of model
        parameters.
        
        This is a wrapper for :class:`PPXFModel`, and very similar to
        :func:`mangadap.proc.ppxffit.PPXFFit.contruct_models`.

        The input velocities are expected to be cz, not "ppxf"
        (pixelized) velocities.

        If redshift_only is true, the provided dispersion is set to 1e-9
        km/s, which is numerically identical to 0 (i.e., just shifting
        the spectrum) in the tested applications.  However, beware that
        this is a HARDCODED number.

        """
        if redshift_only and dispersion_corrections is not None:
            raise ValueError('redshift_only and dispersion_corrections are mutually exclusive.')
        if deredshift:
            raise NotImplementedError('Cannot yet deredshift models.')
        
        # Check the spectral sampling
        velscale_ratio = int(numpy.around(spectrum_velocity_scale(obj_wave)
                                / spectrum_velocity_scale(stpl_wave)))
        _velscale, _velscale_ratio \
                = PPXFFit.check_pixel_scale(stpl_wave, obj_wave, velscale_ratio=velscale_ratio,
                                             dvtol=dvtol)
        # Check the input spectra
        obj_flux = numpy.zeros(obj_flux_shape, dtype=float)
        _obj_wave, _obj_flux, _, _ = PPXFFit.check_objects(obj_wave, obj_flux)
        nobj = _obj_flux.shape[0]
        _stpl_wave, _stpl_flux, _ = PPXFFit.check_templates(stpl_wave, stpl_flux,
                                                          velscale_ratio=_velscale_ratio)
        nstpl = _stpl_flux.shape[0]

        # Construct the emission-line templates, just to determine how
        # many of them there should be
        etpl = EmissionLineTemplates(_stpl_wave, numpy.ones(_stpl_wave.shape, dtype=float),
                                     emldb=emission_lines, quiet=True)
        _ntpl = nstpl+etpl.ntpl
        ngas_comp = numpy.amax(etpl.comp)+1
        
        # Check the shape of the input model parameter database
        if model_fit_par['BINID'].size != nobj:
            raise ValueError('Incorrect number of model-parameter sets.')
        if model_fit_par['TPLWGT'].shape[1] != _ntpl:
            raise ValueError('The number of weights does not match the number of templates.')

        # Get the input pixel shift between the object and template
        # wavelength vectors; interpretted by pPXF as a base velocity
        # shift between the two.  For Sasuke, start and end should be
        # the same, but leave it as general.
        vsyst = numpy.array([ -PPXFFit.ppxf_tpl_obj_voff(_stpl_wave, _obj_wave[s:e], _velscale,
                                                         velscale_ratio=_velscale_ratio)
                                for s,e in zip(model_fit_par['BEGPIX'], model_fit_par['ENDPIX'])])

        # Get the additive and multiplicative degree of the polynomials
        degree = model_fit_par['ADDCOEF'].shape
        degree = -1 if len(degree) == 1 else degree[1]-1
        mdegree = model_fit_par['MULTCOEF'].shape
        mdegree = 0 if len(mdegree) == 1 else mdegree[1]
        moments = model_fit_par['KIN'].shape[1]

        # Determine what to use for the reddening; assume it is fit if
        # the multiplicative polynomial is not fit.  Even if neither are
        # fit, Sasuke sets the reddening to 0 such that the correct
        # model will be returned
        reddening = mdegree == 0

        # Determine the number of stellar kinematic moments.  The number
        # of gas moments is currently HARD-WIRED to be 2, meaning that
        # the number stellar kinematic moments is just:
        smoments = model_fit_par['KIN'].shape[1] - 2*ngas_comp

        # Only produce selected models
        skip = numpy.zeros(nobj, dtype=bool) if select is None else numpy.invert(select)

        # Instantiate the output model array
        models = numpy.ma.zeros(_obj_flux.shape, dtype=numpy.float)
        # Initially mask everything
        models[:,:] = numpy.ma.masked

        # Get the kinematics to use
        kin = model_fit_par['KIN'].copy()[:,:smoments]
        kin[:,0],_ = PPXFFit.revert_velocity(model_fit_par['KIN'][:,0],
                                             model_fit_par['KINERR'][:,0])
        if redshift_only:
            kin[:,1] = 1e-9
        elif dispersion_corrections is not None:
            kin[:,1] = numpy.ma.sqrt(numpy.square(model_fit_par['KIN'][:,1]) 
                                        - numpy.square(dispersion_corrections)).filled(1e-9)
#        if deredshift:
#            kin[:,0] = 0.0

        # Construct the model for each (selected) object spectrum
        for i in range(nobj):
            if skip[i]:
                continue
            ebv = model_fit_par['EBV'][i] if reddening else None

            print('Constructing continuum for spectrum: {0}/{1}'.format(i+1,nobj), end='\r')

            # This is redeclared every iteration to allow for the
            # starting and ending pixels to be different (annoying); as
            # will the velocity offset; this means that the FFT of the
            # templates is recalculated at every step...
            f = PPXFModel(_stpl_flux.T,
                          _obj_flux.data[i,model_fit_par['BEGPIX'][i]:model_fit_par['ENDPIX'][i]],
                          _velscale, velscale_ratio=_velscale_ratio, vsyst=vsyst[i],
                          moments=smoments, degree=degree, mdegree=mdegree, reddening=ebv)

            models[i,model_fit_par['BEGPIX'][i]:model_fit_par['ENDPIX'][i]] \
                        = f(kin[i,:], model_fit_par['TPLWGT'][i,:nstpl],
                            addpoly=None if degree < 0 else model_fit_par['ADDCOEF'][i,:],
                            multpoly=None if mdegree < 1 else model_fit_par['MULTCOEF'][i,:],
                            reddening=ebv)

#            if i == 0 or redshift_only:
#                pyplot.plot(_obj_wave, _obj_flux[i,:])
#                pyplot.plot(_obj_wave, models[i,:])
#                pyplot.show()
#                if redshift_only:
#                    exit()
#
        print('Constructing continuum for spectrum: {0}/{0}'.format(nobj))
        return models

