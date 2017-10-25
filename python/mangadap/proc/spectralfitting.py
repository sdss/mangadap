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
    | **25 Oct 2017**: (KBW) Added PLY columns to emission-line database

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
from ..util.constants import DAPConstants
from ..util import lineprofiles
from ..par.emissionlinedb import EmissionLineDB
from .bandpassfilter import emission_line_equivalent_width

# For debugging
from matplotlib import pyplot

# BASE CLASS -----------------------------------------------------------
class SpectralFitting():
    """
    Base class for spectral fitting.
    """
    def __init__(self, fit_type, bitmask=None, par=None):
        self.fit_type = fit_type
        if bitmask is not None and not isinstance(bitmask, BitMask):
            raise TypeError('Input bit mask must have type BitMask.')
        self.bitmask = bitmask
        self.par = par


# ----------------------------------------------------------------------
class StellarKinematicsFit(SpectralFitting):
    """
    Base class for fitting stellar kinematics.
    """
    def __init__(self, fit_method, bitmask, par=None):
        SpectralFitting.__init__(self, 'stellar_kinematics', bitmask=bitmask, par=par)
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
                 ('TPLWGTERR',numpy.float,(ntpl,)),
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
                 ('SIGMACORR_SRES',numpy.float),
                 ('SIGMACORR_EMP',numpy.float)
               ]

# ----------------------------------------------------------------------
class CompositionFit(SpectralFitting):
    """
    Base class for fitting the spectral composition.
    """
    def __init__(self, fit_method, bitmask, par=None):
        SpectralFitting.__init__(self, 'composition', bitmask=bitmask, par=par)
        self.fit_method = fit_method


# ----------------------------------------------------------------------
class EmissionLineFit(SpectralFitting):
    """
    Base class for fitting emission lines.
    """
    def __init__(self, fit_method, bitmask=None, par=None):
        SpectralFitting.__init__(self, 'emission_line', bitmask=bitmask, par=par)
        self.fit_method = fit_method


    @staticmethod
    def _per_emission_line_dtype(neml, nkin, mask_dtype):
        r"""
        Construct the record array data type for the output fits
        extension.
        """
        return [ ('BINID', numpy.int),
                 ('BINID_INDEX', numpy.int),
                 ('FIT_INDEX', numpy.int, (neml,)),
                 ('MASK', mask_dtype, (neml,)),
                 ('FLUX', numpy.float, (neml,)),
                 ('FLUXERR', numpy.float, (neml,)),
                 ('KIN', numpy.float, (neml,nkin)),
                 ('KINERR', numpy.float, (neml,nkin)),
                 ('SIGMACORR', numpy.float, (neml,)),
                 ('CONTAPLY', numpy.float, (neml,)),
                 ('CONTMPLY', numpy.float, (neml,)),
                 ('BMED', numpy.float, (neml,)),
                 ('RMED', numpy.float, (neml,)),
                 ('EWCONT', numpy.float, (neml,)),
                 ('EW', numpy.float, (neml,)),
                 ('EWERR', numpy.float, (neml,))
               ]


    @staticmethod
    def get_spectra_to_fit(spectra, pixelmask=None, select=None, error=False):
        """
        Get the spectra to fit during the emission-line fitting.

        Args:
            spectra (:class:`mangadap.drpfits.DRPFits`,
                :class:`mangadap.proc.spatiallybinnedspectra.SpatiallBinnedSpectra`):
                Object with the spectra to fit.  Can be one of the
                provided objects.  This works because both have
                `copy_to_masked_array` and `do_not_fit_flags` methods.
            pixelmask
                (:class:`mangadap.util.pixelmask.SpectralPixelMask`):
                (**Optional**) Pixel mask to apply.
            select (numpy.ndarray): (**Optional**) Select specific
                spectra to return.
            error (bool): (**Optional**) Return 1-sigma errors instead
                of inverse variance.

        Returns:
            numpy.ma.MaskedArray: Two masked arrays: the flux data and
            the uncertainties, either as 1-sigma error or the inverse
            variance.

        """
        # Grab the spectra
        flux = spectra.copy_to_masked_array(flag=spectra.do_not_fit_flags())
        ivar = spectra.copy_to_masked_array(ext='IVAR', flag=spectra.do_not_fit_flags())
        nspec = flux.shape[0]

        # Convert inverse variance to error
        if error:
            ivar = numpy.ma.power(ivar, -0.5)
            flux[numpy.ma.getmaskarray(ivar)] = numpy.ma.masked

        # Mask any pixels in the pixel mask
        if pixelmask is not None:
            indx = pixelmask.boolean(spectra['WAVE'].data, nspec=nspec)
            flux[indx] = numpy.ma.masked
            ivar[indx] = numpy.ma.masked
        
        _select = numpy.ones(nspec, dtype=bool) if select is None else select
        return flux[_select,:], ivar[_select,:]


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

        sinst = c / interpolator((cz/c + 1.0) * restwave)/DAPConstants.sig2fwhm

#        pyplot.plot(wave, c/sres/DAPConstants.sig2fwhm, lw=2)
#        pyplot.plot(wave, 2.5*astropy.constants.c.to('km/s').value/wave/DAPConstants.sig2fwhm,
#                    color='k', lw=1)
#        pyplot.scatter((cz/c + 1.0) * restwave, sinst, color='C1', marker='.', lw=0, s=100)
#        pyplot.show()

        return c / interpolator((cz/c + 1.0) * restwave)/DAPConstants.sig2fwhm


    @staticmethod
    def check_emission_line_database(emldb, wave=None, check_par=True):
        r"""
        Check the emission-line database.  Modes are checked by
        :class:`mangadap.par.emissionlinedb.EmissionLinePar`, and the
        indices are checked to be unique by
        :class:`mangadap.par.emissionlinedb.EmissionLineDB`.

            - The type of the object must be
              :class:`mangadap.par.emissionlinedb.EmissionLineDB`
            - The provided profile type of each line must be a defined
              class.
            - At least one line must have mode=`f`
            - All tied lines must be tied to a line with a correctly
              specified index.
            - Warnings will be provided for any line with a centroid
              that falls outside of the provided wavelength range.
            - The database must provide at least one valid line.

        Args:
            emldb (:class:`mangadap.par.emissionlinedb.EmissionLineDB'):
                Emission-line database.
            wave (array-like)

        Raises:
            TypeError: Raised if the provided object is not an instance
                of :class:`mangadap.par.emissionlinedb.EmissionLineDB`.
            ValueError: Raised if any line has a mode of `x` or if the
                database does not provide a valid definition for any
                templates.
            NameError: Raised if a defined profile type is not known.
        """

        # Check the object type
        if not isinstance(emldb, EmissionLineDB):
            raise TypeError('Emission lines must be defined using an EmissionLineDB object.')

        # Check the profile type
        unique_profiles = numpy.unique(emldb['profile'])
        for u in unique_profiles:
            try:
                eval('lineprofiles.'+u)
            except NameError as e:
                raise NameError('Profile type {0} not defined in'
                                'mangadap.util.lineprofiles!'.format(u))

        # There must be one primary line
        if not numpy.any([m[0] == 'f' for m in emldb['mode']]):
            raise ValueError('At least one line in the database must have mode=f.')

        # Check that there are lines to fit
        lines_to_fit = emldb['action'] == 'f'
        if numpy.sum(lines_to_fit) == 0:
            raise ValueError('No lines to fit in the database!')
        if wave is not None:
            _wave = numpy.asarray(wave)
            if len(_wave.shape) != 1:
                raise ValueError('Provided wavelengths must be a single vector.')
            lines_in_range = numpy.array([rw > _wave[0] and rw < _wave[-1] 
                                                for rw in emldb['restwave']]) 
            if numpy.sum(lines_to_fit & lines_in_range) == 0:
                raise ValueError('No lines to fit in the provided spectral range!')

        # Check that the tied line indices exist in the database
        for m in emldb['mode']:
            if m[0] == 'f':
                continue
            tied_index = int(m[1:])
            if numpy.sum(emldb['index'] == tied_index) == 0:
                raise ValueError('No line with index={0} to tie to!'.format(tied_index))

        # Only check the provided parameters if requested
        if not check_par:
            return

        # Check the provided parameters, fix flags, and bounds
        for i in range(emldb.neml):
            profile = eval('lineprofiles.'+emldb['profile'][i])
            npar = len(profile.param_names)
            if emldb['par'][i].size != npar*emldb['ncomp'][i]:
                raise ValueError('Provided {0} parameters, but expected {1}.'.format(
                                  emldb['par'][i].size, npar*emldb['ncomp'][i]))
            if emldb['fix'][i].size != npar*emldb['ncomp'][i]:
                raise ValueError('Provided {0} fix flags, but expected {1}.'.format(
                                  emldb['fix'][i].size, npar*emldb['ncomp'][i]))
            if numpy.any([f not in [0, 1] for f in emldb['fix'][i] ]):
                warnings.warn('Fix values should only be 0 or 1; non-zero values interpreted as 1.')
            if emldb['lobnd'][i].size != npar*emldb['ncomp'][i]:
                raise ValueError('Provided {0} lower bounds, but expected {1}.'.format(
                                  emldb['lobnd'][i].size, npar*emldb['ncomp'][i]))
            if emldb['hibnd'][i].size != npar*emldb['ncomp'][i]:
                raise ValueError('Provided {0} upper bounds, but expected {1}.'.format(
                                  emldb['hibnd'][i].size, npar*emldb['ncomp'][i]))
            if emldb['log_bnd'][i].size != npar*emldb['ncomp'][i]:
                raise ValueError('Provided {0} log boundaries designations, but expected '
                                 '{1}.'.format(emldb['log_bnd'][i].size, npar*emldb['ncomp'][i]))


    @staticmethod
    def measure_equivalent_width(wave, flux, emission_lines, model_eml_par, mask=None,
                                 redshift=None, bitmask=None, checkdb=True):
        """
        The flux array is expected to have size Nspec x Nwave.

        Provided previous emission-line fits, this function adds the
        equivalent width measurements to the output database.

        Errors currently *do not* include the errors in the continuum
        measurement; only the provided error in the flux.

        Raises:
            ValueError: Raised if the length of the spectra, errors, or
                mask does not match the length of the wavelength array;
                raised if the wavelength, redshift, or dispersion arrays
                are not 1D vectors; and raised if the number of
                redshifts or dispersions is not a single value or the
                same as the number of input spectra.
        """

        # Check the input emission-line database
        if checkdb:
            EmissionLineFit.check_emission_line_database(emission_lines)

        # If the redshift is NOT provided, use the fitted velocity
        _redshift = numpy.mean(model_eml_par['KIN'][:,:,0]/astropy.constants.c.to('km/s').value,
                               axis=1) if redshift is None else redshift

        # Calculate the wavelength at which to measure the continuum,
        # matching what is done by
        # :class:`mangadap.proc.emissionlineMoments.EmissionLineMoments`
        line_center = (1+_redshift)[:,None]*emission_lines['restwave'][None,:]

        # Compute the equivalent widths.  The checking done by
        # EmissionLineFit.check_and_prep_input is *identical* to what is
        # done within emission_line_equivalent_width()
        model_eml_par['BMED'], model_eml_par['RMED'], pos, model_eml_par['EWCONT'], \
                model_eml_par['EW'], model_eml_par['EWERR'] \
                        = emission_line_equivalent_width(wave, flux,
                                                         emission_lines['blueside'],
                                                         emission_lines['redside'], line_center,
                                                         model_eml_par['FLUX'], mask=mask,
                                                         redshift=_redshift,
                                                         line_flux_err=model_eml_par['FLUXERR'])

        # Flag non-positive measurements
        if bitmask is not None:
            model_eml_par['MASK'][numpy.invert(pos)] \
                    = bitmask.turn_on(model_eml_par['MASK'][numpy.invert(pos)],
                                      'NON_POSITIVE_CONTINUUM')



