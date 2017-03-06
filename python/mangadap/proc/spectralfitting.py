# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Implements a few base classes used during spectral fitting procedures.

*License*:
    Copyright (c) 2015, SDSS-IV/MaNGA Pipeline Group
        Licensed under BSD 3-clause license - see LICENSE.rst

*Source location*:
    $MANGADAP_DIR/python/mangadap/proc/spectralfitting.py

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

        import numpy

*Class usage examples*:
        Add examples

*Revision history*:
    | **14 Apr 2016**: Implementation begun by K. Westfall (KBW)
    | **19 Apr 2016**: (KBW) First version
    | **26 Apr 2016**: (KBW) Moved PPXFFit to a separate file (ppxffit.py)
    | **03 Nov 2016**: (KBW) Added USETPL column to stellar kinematics
        output table.

.. _astropy.io.fits.hdu.hdulist.HDUList: http://docs.astropy.org/en/v1.0.2/io/fits/api/hdulists.html
.. _glob.glob: https://docs.python.org/3.4/library/glob.html


"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import sys
import warnings
if sys.version > '3':
    long = int

import numpy
from scipy.interpolate import interp1d
import astropy.constants

from ..util.bitmask import BitMask
from ..util.constants import constants

__author__ = 'Kyle B. Westfall'
# Add strict versioning
# from distutils.version import StrictVersion



# BASE CLASS -----------------------------------------------------------
class SpectralFitting():
    """
    Base class for spectral fitting.
    """
    def __init__(self, fit_type, bitmask, par=None):
        self.fit_type = fit_type
        if not isinstance(bitmask, BitMask):
            raise TypeError('Input bit mask must have type BitMask.')
        self.bitmask = bitmask
        self.par = par


# ----------------------------------------------------------------------
class StellarKinematicsFit(SpectralFitting):
    """
    Base class for fitting stellar kinematics.
    """
    def __init__(self, fit_method, bitmask, par=None):
        SpectralFitting.__init__(self, 'stellar_kinematics', bitmask, par=par)
        self.fit_method = fit_method

    @staticmethod
    def _per_stellar_kinematics_dtype(ntpl, nadd, nmult, nkin, mask_dtype):
        r"""
        Construct the record array data type for the output fits
        extension.
        """
        return [ ('BINID',numpy.int),
                 ('BINID_INDEX',numpy.int),
                 ('MASK',mask_dtype),
                 ('BEGPIX', numpy.int),
                 ('ENDPIX', numpy.int),
                 ('NPIXTOT',numpy.int),
                 ('NPIXFIT',numpy.int),
                 ('TPLWGT',numpy.float,(ntpl,)),
                 ('USETPL',numpy.bool,(ntpl,)),
                 ('ADDCOEF',numpy.float,(nadd,)) if nadd > 1 else ('ADDCOEF',numpy.float),
                 ('MULTCOEF',numpy.float,(nmult,)) if nmult > 1 else ('MULTCOEF',numpy.float),
                 ('KININP',numpy.float,(2,)),
                 ('KIN',numpy.float,(nkin,)),
                 ('KINERR',numpy.float,(nkin,)),
                 ('CHI2',numpy.float),
                 ('RCHI2',numpy.float),
                 ('ROBUST_RCHI2',numpy.float),
                 ('RMS',numpy.float),
                 ('ABSRESID',numpy.float,(5,)),
                 ('FRMS',numpy.float),
                 ('FABSRESID',numpy.float,(5,)),
                 ('SIGMACORR',numpy.float)
               ]

# ----------------------------------------------------------------------
class CompositionFit(SpectralFitting):
    """
    Base class for fitting the spectral composition.
    """
    def __init__(self, fit_method, bitmask, par=None):
        SpectralFitting.__init__(self, 'composition', bitmask, par=par)
        self.fit_method = fit_method


# ----------------------------------------------------------------------
class EmissionLineFit(SpectralFitting):
    """
    Base class for fitting emission lines.
    """
    def __init__(self, fit_method, bitmask, par=None):
        SpectralFitting.__init__(self, 'emission_line', bitmask, par=par)
        self.fit_method = fit_method


    @staticmethod
    def _per_emission_line_dtype(neml, nkin, mask_dtype):
        r"""
        Construct the record array data type for the output fits
        extension.
        """
        return [ ('BINID', numpy.int),
                 ('BINID_INDEX', numpy.int),
                 ('WIN_INDEX', numpy.int, (neml,)),
                 ('MASK', mask_dtype, (neml,)),
                 ('FLUX', numpy.float, (neml,)),
                 ('FLUXERR', numpy.float, (neml,)),
                 ('KIN', numpy.float, (neml,nkin)),
                 ('KINERR', numpy.float, (neml,nkin)),
                 ('SINST', numpy.float, (neml,)),
                 ('BMED', numpy.float, (neml,)),
                 ('RMED', numpy.float, (neml,)),
                 ('EWCONT', numpy.float, (neml,)),
                 ('EW', numpy.float, (neml,)),
                 ('EWERR', numpy.float, (neml,))
               ]

#                 ('EWOF',numpy.float,neml),
#                 ('EWOFERR',numpy.float,neml)
#                 ('EW',numpy.float,neml),
#                 ('EWERR',numpy.float,neml)
#                 ('COMBKIN',numpy.float,nkin),
#                 ('COMBKINERR',numpy.float,nkin),
#                 ('COMBKINSTD',numpy.float,nkin),
#                 ('COMBKINN',numpy.int),

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
    def get_binned_data(binned_spectra, pixelmask=None, select=None):
        # Grab the spectra
        flux = binned_spectra.copy_to_masked_array(flag=binned_spectra.do_not_fit_flags())
        ivar = binned_spectra.copy_to_masked_array(ext='IVAR',
                                                   flag=binned_spectra.do_not_fit_flags())
        # Mask any pixels in the pixel mask
        if pixelmask is not None:
            indx = pixelmask.boolean(binned_spectra['WAVE'].data, nspec=binned_spectra.nbins)
            flux[indx] = numpy.ma.masked
            ivar[indx] = numpy.ma.masked
        
        _select = numpy.ones(binned_spectra.nbins, dtype=bool) if select is None else select
        return flux[select,:], ivar[select,:]


    @staticmethod
    def check_and_prep_input(wave, flux, ivar=None, mask=None, sres=None, continuum=None,
                             redshift=None, dispersion=None, default_dispersion=100.0):
        """
        Check the input used for emission-line measurements.

        inverse variance is converted to 1-sigma error

        mask must be a boolean array.

        sres can be a single vector, but will be returned as an array
        with a size that matches flux.

        output all converted to masked arrays with at least two dimensions

        """
        # Check the input wavelength and flux shapes
        if len(wave.shape) != 1:
            raise ValueError('Input wavelengths must be a single vector.')

        # Check the mask shape
        if mask is not None and mask.shape != flux.shape:
            raise ValueError('Input mask must have the same shape as the flux array.')

        # Convert the input arrays to masked arrays if they aren't
        # already, and compare the array shapes
        _flux = numpy.ma.atleast_2d(flux if isinstance(flux, numpy.ma.MaskedArray) \
                                            else numpy.ma.MaskedArray(flux, mask=mask))
        if len(wave) != _flux.shape[1]:
            raise ValueError('Wavelength vector does not match shape of the flux array.')
        expected_shape = _flux.shape
        nspec = expected_shape[0]

        if ivar is None:
            _err = None
        else:
            _ivar = numpy.ma.atleast_2d(ivar if isinstance(ivar, numpy.ma.MaskedArray) \
                                                else numpy.ma.MaskedArray(ivar, mask=mask))
            if _ivar.shape != expected_shape:
                raise ValueError('Input ivar array must be the same shape as the flux array.')
            _err = numpy.ma.sqrt(1.0 /_ivar)

        if sres is None:
            _sres = None
        else:
            _sres = numpy.array([sres]*nspec) if len(sres.shape) == 1 else sres.copy()
            if _sres.shape != expected_shape:
                raise ValueError('Input sres array must match the flux array, either as a single '
                                 'vector as an 2D array.')

        if continuum is None:
            _continuum = None
        else:
            _continuum = numpy.ma.atleast_2d(continuum \
                                    if isinstance(continuum, numpy.ma.MaskedArray) \
                                    else numpy.ma.MaskedArray(continuum, mask=mask))
            if _continuum.shape != expected_shape:
                raise ValueError('Input continuum array must be the same shape as the flux array.')

        # Check the input redshifts and dispersions
        _redshift = numpy.zeros(nspec, dtype=numpy.float) if redshift is None else redshift
        if len(_redshift) != nspec:
            raise ValueError('Must provide one redshift per input spectrum.')

        _dispersion = numpy.full(nspec, default_dispersion, dtype=numpy.float) \
                            if dispersion is None else dispersion
        if len(_dispersion) != nspec:
            raise ValueError('Must provide one dispersion per input spectrum.')

        # Return all arrays (even if None)
        return _flux, _err, _sres, _continuum, _redshift, _dispersion


    @staticmethod
    def subtract_continuum(flux, continuum):
        """
        Subtract the continuum.  Does not check that shapes match.
        Returns the continuum subtracted flux and a boolean array
        setting where the continuum is not defined.
        """
        # Allow continuum to be None
        if continuum is None:
            return flux, None

        # Get where the continuum is masked but the spectra are not
        no_continuum = numpy.invert(numpy.ma.getmaskarray(flux)) & numpy.ma.getmaskarray(continuum)

        # Subtract the continuum (ensure output is a masked array)
        continuum_subtracted_flux = numpy.ma.subtract(flux, continuum)

        # Unmask regions where only the continuum is masked.
        continuum_subtracted_flux.mask[no_continuum] = False

        return continuum_subtracted_flux, no_continuum


    @staticmethod
    def instrumental_dispersion(wave, sres, restwave, cz):
        """
        Determine the instrumental dispersion for a set of rest
        wavelengths and velocities.

        Args:

            wave (numpy.ndarray): Vector with the wavelengths of the
                spectrum.

            sres (numpy.ndarray): Vector with the spectral resolution as
                a function of wavelength.

            restwave (float, numpy.ndarray): Rest wavelengths for a set
                of measured lines.

            cz (float, numpy.ndarray): Redshifts (in km/s) of each or
                all lines.

        Returns:
            numpy.ndarray : The instrumental dispersions for each
                provided line.
        """
        # Check input
        if len(wave.shape) != 1:
            raise ValueError('Input wavelength must be a 1D vector.')
        if wave.shape != sres.shape:
            raise ValueError('Input wavelength and resolution vectors must have the same shape.')
        
        nwave = wave.size
        _restwave = numpy.atleast_1d(restwave)
        nline = _restwave.size
        _cz = numpy.atleast_1d(cz)
        if _cz.size not in [ 1, nline ]:
            raise ValueError('Must provide single redshift or one redshift per line.')
        _cz = numpy.full(nline, cz, dtype=float) if _cz.size == 1 else _cz

        interpolator = interp1d(wave, sres, fill_value='extrapolate', assume_sorted=True)
        c = astropy.constants.c.to('km/s').value
        return c / interpolator((cz/c + 1.0) * restwave)/constants().sig2fwhm






