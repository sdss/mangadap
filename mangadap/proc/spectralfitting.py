# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Implements a few base classes used during spectral fitting procedures.

----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import warnings

from IPython import embed

import numpy
from scipy.interpolate import interp1d
import astropy.constants

from ..util.bitmask import BitMask
from ..util.datatable import DataTable
from ..util.constants import DAPConstants
from ..util import lineprofiles
from ..util.filter import interpolate_masked_vector
from ..par.emissionlinedb import EmissionLineDB
from .bandpassfilter import emission_line_equivalent_width, passband_median

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
# TODO: Minimize this to what any kinematics fit should *absolutely*
# provide. All the rest should go into separate DataTables specific to
# the fitting class.
class StellarKinematicsFitDataTable(DataTable):
    """
    Class defining the data table to hold the stellar kinematics fit.

    Table includes:

    .. include:: ../tables/stellarkinematicsfitdatatable.rst

    Args:
        ntpl (:obj:`int`):
            Number of templates used to model the stellar continuum.
        nadd (:obj:`int`):
            Number of coefficients in any additive polynomial
            included in the fit. Can be 0.
        nmult (:obj:`int`):
            Number of coefficients in any multiplicative polynomial
            included in the fit. Can be 0.
        nkin (:obj:`int`):
            Number of kinematic moments included in the fit. Note
            that the number of moments used as input guesses is
            *always* assumed to be 2.
        mask_dtype (:obj:`type`):
            The data type used for the maskbits (e.g., numpy.int16).
            Typically this would be set by
            :func:`mangadap.util.bitmask.BitMask.minimum_dtype`.
        shape (:obj:`int`, :obj:`tuple`, optional):
            The shape of the initial array. If None, the data array
            will not be instantiated; use :func:`init` to initialize
            the data array after instantiation.
    """
    def __init__(self, ntpl=1, nadd=0, nmult=0, nkin=2, mask_dtype=numpy.int16, shape=None):
        # NOTE: This should require python 3.7 to make sure that this
        # is an "ordered" dictionary.
        datamodel = dict(BINID=dict(typ=int, shape=None, descr='Bin ID number'),
                         BINID_INDEX=dict(typ=int, shape=None, descr='0-indexed number of bin'),
                         MASK=dict(typ=mask_dtype, shape=None, descr='Maskbit value'),
                         BEGPIX=dict(typ=int, shape=None,
                                     descr='Index of the first pixel included in the fit'),
                         ENDPIX=dict(typ=int, shape=None,
                                     descr='Index of the pixel just beyond the last pixel '
                                           'included in fit'),
                         NPIXTOT=dict(typ=int, shape=None, 
                                      descr='Total number of pixels in the spectrum to be fit.'),
                         NPIXFIT=dict(typ=int, shape=None,
                                      descr='Number of pixels used by the fit.'),
                         TPLWGT=dict(typ=float, shape=(ntpl,),
                                     descr='Optimal weight of each template.'),
                         TPLWGTERR=dict(typ=float, shape=(ntpl,),
                                        descr='Nominal error in the weight of each template.'),
                         USETPL=dict(typ=bool, shape=(ntpl,),
                                     descr='Flag that each template was included in the fit.'),
                         ADDCOEF=dict(typ=float, shape=(nadd,) if nadd > 1 else None,
                                      descr='Coefficients of the additive polynomial, if '
                                            'included.'),
                         MULTCOEF=dict(typ=float, shape=(nmult,) if nmult > 1 else None,
                                       descr='Coefficients of the multiplicative polynomial, if '
                                             'included.'),
                         KININP=dict(typ=float, shape=(2,),
                                     descr='Input guesses for the kinematics'),
                         KIN=dict(typ=float, shape=(nkin,),
                                     descr='Best-fitting stellar kinematics'),
                         KINERR=dict(typ=float, shape=(nkin,),
                                     descr='Errors in the best-fitting stellar kinematics'),
                         CHI2=dict(typ=float, shape=None,
                                   descr='Chi-square figure-of-merit for the fit'),
                         RCHI2=dict(typ=float, shape=None,
                                    descr='Reduced chi-square figure-of-merit for the fit'),
                         CHIGRW=dict(typ=float, shape=(5,),
                                    descr='Value of the error-normalized residuals at 0, 68%, '
                                          '95%, 99%, and 100% growth'),
                         RMS=dict(typ=float, shape=None,
                                  descr='Root-mean-square of the fit residuals.'),
                         RMSGRW=dict(typ=float, shape=(5,),
                                    descr='Value of absolute value of the fit residuals at 0, '
                                          '68%, 95%, 99%, and 100% growth'),
                         FRMS=dict(typ=float, shape=None,
                                   descr='Root-mean-square of the fractional residuals '
                                         '(i.e., residuals/model).'),
                         FRMSGRW=dict(typ=float, shape=(5,),
                                      descr='Value of absolute value of the fractional residuals '
                                            'at 0, 68%, 95%, 99%, 100% growth'),
                         SIGMACORR_SRES=dict(typ=float, shape=None,
                                             descr='Quadrature correction for the stellar '
                                                   'velocity dispersion determined by the mean '
                                                   'difference in spectral resolution between '
                                                   'galaxy and template data.'),
                         SIGMACORR_EMP=dict(typ=float, shape=None,
                                            descr='Quadrature correciton for the stellar '
                                                  'velocity dispersion determined by fitting the '
                                                  'optimal template to one resolution matched to '
                                                  'the galaxy data.'))

        keys = list(datamodel.keys())
        super(StellarKinematicsFitDataTable,
                self).__init__(keys, [datamodel[k]['typ'] for k in keys],
                               element_shapes=[datamodel[k]['shape'] for k in keys],
                               descr=[datamodel[k]['descr'] for k in keys],
                               shape=shape)


class StellarKinematicsFit(SpectralFitting):
    """
    Base class for fitting stellar kinematics.
    """
    def __init__(self, fit_method, bitmask, par=None):
        SpectralFitting.__init__(self, 'stellar_kinematics', bitmask=bitmask, par=par)
        self.fit_method = fit_method

    @staticmethod
    def init_datatable(ntpl, nadd, nmult, nkin, mask_dtype, shape=None):
        return StellarKinematicsFitDataTable(ntpl=ntpl, nadd=nadd, nmult=nmult, nkin=nkin,
                                             mask_dtype=mask_dtype, shape=shape)


# ----------------------------------------------------------------------
class CompositionFit(SpectralFitting):
    """
    Base class for fitting the spectral composition.
    """
    def __init__(self, fit_method, bitmask, par=None):
        SpectralFitting.__init__(self, 'composition', bitmask=bitmask, par=par)
        self.fit_method = fit_method


# ----------------------------------------------------------------------
class EmissionLineFitDataTable(DataTable):
    """
    Primary data table with the results of the parameterized emission-line fit.

    Table includes:

    .. include:: ../tables/emissionlinefitdatatable.rst

    Args:
        neml (:obj:`int`):
            Number of emission lines being fit
        nkin (:obj:`int`):
            Number of kinematic parameters (e.g., 2 for V and sigma)
        mask_dtype (:obj:`type`):
            The data type used for the maskbits (e.g., numpy.int16).
            Typically this would be set by
            :func:`mangadap.util.bitmask.BitMask.minimum_dtype`.
        shape (:obj:`int`, :obj:`tuple`, optional):
            The shape of the initial array. If None, the data array
            will not be instantiated; use :func:`init` to initialize
            the data array after instantiation.
    """
    def __init__(self, neml=1, nkin=2, mask_dtype=numpy.int16, shape=None):
        # NOTE: This should require python 3.7 to make sure that this
        # is an "ordered" dictionary.
        datamodel = dict(BINID=dict(typ=int, shape=None, descr='Bin ID number'),
                         BINID_INDEX=dict(typ=int, shape=None, descr='0-indexed number of bin'),
                         FIT_INDEX=dict(typ=int, shape=(neml,),
                                        descr='The index in the fit database associated with '
                                              'each emission line.'),
                         MASK=dict(typ=mask_dtype, shape=(neml,),
                                   descr='Maskbit value for each emission line.'),
                         FLUX=dict(typ=float, shape=(neml,),
                                   descr='The best-fitting flux of the emission line.'),
                         FLUXERR=dict(typ=float, shape=(neml,),
                                      descr='The error in the best-fitting emission-line flux'),
                         KIN=dict(typ=float, shape=(neml,nkin),
                                  descr='The best-fitting kinematics in each emission line'),
                         KINERR=dict(typ=float, shape=(neml,nkin),
                                     descr='The error in the best-fitting emission-line '
                                           'kinematics'),
                         SIGMACORR=dict(typ=float, shape=(neml,),
                                        descr='Quadrature correction in the emission-line '
                                              'velocity dispersion'),
                         SIGMAINST=dict(typ=float, shape=(neml,),
                                        descr='Dispersion of the instrumental line-spread '
                                              'function at the location of each emission line.'),
                         SIGMATPL=dict(typ=float, shape=(neml,),
                                       descr='Dispersion of the instrumental line-spread function '
                                             'of the emission-line templates.'),
                         CONTAPLY=dict(typ=float, shape=(neml,),
                                       descr='The value of any additive polynomial included in '
                                             'the fit at the location of each emission line'),
                         CONTMPLY=dict(typ=float, shape=(neml,),
                                       descr='The value of any multiplicative polynomial included '
                                             'in the fit at the location of each emission line'),
                         CONTRFIT=dict(typ=float, shape=(neml,),
                                       descr='The value of any extinction curve included in the '
                                             'fit at the location of each emission line'),
                         LINE_PIXC=dict(typ=int, shape=(neml,),
                                        descr='The integer pixel nearest the center of each '
                                              'emission line.'),
                         AMP=dict(typ=float, shape=(neml,),
                                  descr='The best-fitting amplitude of the emission line.'),
                         ANR=dict(typ=float, shape=(neml,),
                                  descr='The amplitude-to-noise ratio defined as the model '
                                        'amplitude divided by the median noise in the two (blue '
                                        'and red) sidebands defined for the emission line.'),
                         LINE_NSTAT=dict(typ=int, shape=(neml,),
                                         descr='The number of pixels included in the fit metric '
                                               'calculations (LINE_RMS, LINE_FRMS, LINE_CHI2) '
                                               'near each emission line.'),
                         LINE_RMS=dict(typ=float, shape=(neml,),
                                       descr='The root-mean-square residual of the model fit '
                                             'near each emission line.'),
                         LINE_FRMS=dict(typ=float, shape=(neml,),
                                       descr='The root-mean-square of the fractional residuals '
                                             'of the model fit near each emission line.'),
                         LINE_CHI2=dict(typ=float, shape=(neml,),
                                       descr='The chi-square of the model fit near each emission '
                                             'line.'),
                         BMED=dict(typ=float, shape=(neml,),
                                   descr='The median flux in the blue sideband of each emission '
                                         'line'),
                         RMED=dict(typ=float, shape=(neml,),
                                   descr='The median flux in the red sideband of each emission '
                                         'line'),
                         EWCONT=dict(typ=float, shape=(neml,),
                                     descr='The continuum value interpolated at the '
                                           'emission-line center (in the observed frame) used for '
                                           'the equivalent width measurement.'),
                         EW=dict(typ=float, shape=(neml,),
                                 descr='The equivalent width of each emission line'),
                         EWERR=dict(typ=float, shape=(neml,),
                                    descr='The error in the equivalent width of each emission '
                                          'line'))

        keys = list(datamodel.keys())
        super(EmissionLineFitDataTable,
                self).__init__(keys, [datamodel[k]['typ'] for k in keys],
                               element_shapes=[datamodel[k]['shape'] for k in keys],
                               descr=[datamodel[k]['descr'] for k in keys],
                               shape=shape)


class EmissionLineFit(SpectralFitting):
    """
    Base class for fitting emission lines.
    """
    def __init__(self, fit_method, bitmask=None, par=None):
        SpectralFitting.__init__(self, 'emission_line', bitmask=bitmask, par=par)
        self.fit_method = fit_method

    @staticmethod
    def init_datatable(neml, nkin, mask_dtype, shape=None):
        return EmissionLineFitDataTable(neml=neml, nkin=nkin, mask_dtype=mask_dtype, shape=shape)

    @staticmethod
    def select_binned_spectra_to_fit(binned_spectra, minimum_snr=0.0, stellar_continuum=None,
                                     debug=False):
        """
        Select binned spectra for which to fit emission lines.

        .. todo::
            This could be based on the moment assessment of the
            emission-line S/N instead; for now just based on continuum
            S/N.

        Args:
            binned_spectra
                (:class:`mangadap.proc.spatiallybinnedspectra.SpatiallyBinnedSpectra`):
                Binned spectra to be fit.
            minimum_snr (float): The minimum S/N of the binned spectrum
                to fit; see
                :func:`mangadap.proc.spatiallybinnedspectra.SpatiallyBinnedSpectra.above_snr_limit`.
            stellar_continuum
                (:class:`mangadap.proc.stellarcontinuummodel.StellarContinuumModel`,
                optional): Stellar-continuum models that have been fit
                to the binned spectra, if available.  The current
                function will only return True for spectra that are both
                above the S/N limit and have good stellar-continuum
                models.

        Returns:
            numpy.ndarray: Boolean vector with the spectra in the
            binned_spectra object to fit.
        """

        bins_to_fit = binned_spectra.above_snr_limit(minimum_snr, debug=debug)
        if stellar_continuum is None:
            return bins_to_fit

        # Determine which spectra have a valid stellar continuum fit
        indx = numpy.invert(stellar_continuum.bitmask.flagged(
                                    stellar_continuum['PAR'].data['MASK'],
                                    flag=[ 'NO_FIT', 'INSUFFICIENT_DATA', 'FIT_FAILED']))
        with_good_continuum = numpy.zeros(binned_spectra.nbins, dtype=bool)
        with_good_continuum[stellar_continuum['PAR'].data['BINID_INDEX'][indx]] = True
        bins_to_fit &= with_good_continuum

        if debug:
            warnings.warn('DEBUG!!')
            indx = numpy.arange(len(bins_to_fit))[bins_to_fit]
#            numpy.random.shuffle(indx)
            bins_to_fit[indx[2:]] = False

        return bins_to_fit

    @staticmethod
    def select_spaxels_to_fit(binned_spectra, bins_to_fit=None, debug=False):
        """
        Select spaxels for which to fit emission lines.
        """
        if not debug:
            return binned_spectra.check_fgoodpix()

        warnings.warn('DEBUG!!')

        uniq, indx = numpy.unique(binned_spectra['BINID'].data.ravel(), return_index=True)
        if uniq[0] == -1:
            uniq = uniq[1:]
            indx = indx[1:]
        spaxels_to_fit = numpy.zeros(numpy.prod(binned_spectra['BINID'].data.shape), dtype=bool)
        if bins_to_fit is None:
            spaxels_to_fit[indx] = True
        else:
            spaxels_to_fit[indx[bins_to_fit]] = True
        return spaxels_to_fit

#        spaxels_to_fit = binned_spectra.check_fgoodpix()
#        indx = numpy.arange(len(spaxels_to_fit))[spaxels_to_fit]
#        numpy.random.shuffle(indx)
#        spaxels_to_fit[indx[10:]] = False
#        return spaxels_to_fit

    @staticmethod
    def get_spectra_to_fit(binned_spectra, pixelmask=None, select=None, error=False,
                           original_spaxels=False):
        r"""
        Get the spectra to fit during the emission-line fitting.

        Args:
            binned_spectra
                (:class:`mangadap.proc.spatiallybinnedspectra.SpatiallBinnedSpectra`):
                Object with the spectra to fit.
            pixelmask
                (:class:`mangadap.util.pixelmask.SpectralPixelMask`,
                optional): Pixel mask to apply.
            select (`numpy.ndarray`, optional):
                Select specific spectra to return.  Must have the
                correct shape; cf. `original_spaxels`.
            error (:obj:`bool`, optional):
                Return :math:`1\sigma` errors instead of inverse
                variance.
            original_spaxels (:obj:`bool`, optional):
                Instead of the binned spectra, use the `cube` attribute
                of the `binned_spectra` object to return the original
                spaxels, corrected for Galactic extinction.

        Returns:
            Four objects are returned:
                - (1) The wavelength vector of the spectra,
                - (2) a masked numpy array with the flux data,
                - (3) a masked numpy array with the error data (returned
                  as either inverse variance or :math:`1\sigma`,
                - (4) and an array with the spectral resolution for each
                  spectrum, based on the internal binned spectra
                  parameters.
        """
        # Grab the spectra
        wave = binned_spectra['WAVE'].data.copy()

        if original_spaxels:
            flags = binned_spectra.cube.do_not_fit_flags()
            flux = binned_spectra.cube.copy_to_masked_array(flag=flags)
            ivar = binned_spectra.cube.copy_to_masked_array(attr='ivar', flag=flags)
            flux, ivar = binned_spectra.galext.apply(flux, ivar=ivar, deredden=True)
            sres = binned_spectra.cube.copy_to_array(attr='sres')
        else:
            flags = binned_spectra.do_not_fit_flags()
            flux = binned_spectra.copy_to_masked_array(flag=flags)
            ivar = binned_spectra.copy_to_masked_array(ext='IVAR', flag=flags)
            sres = binned_spectra.copy_to_array(ext='SPECRES')
            sres = numpy.apply_along_axis(interpolate_masked_vector, 1,
                                    numpy.ma.MaskedArray(sres, mask=numpy.invert(sres > 0)))
        nspec = flux.shape[0]

        # Convert inverse variance to error
        if error:
            ivar = numpy.ma.power(ivar, -0.5)
            flux[numpy.ma.getmaskarray(ivar)] = numpy.ma.masked

        # Mask any pixels in the pixel mask
        if pixelmask is not None:
            indx = pixelmask.boolean(binned_spectra['WAVE'].data, nspec=nspec)
            flux[indx] = numpy.ma.masked
            ivar[indx] = numpy.ma.masked
        
        _select = numpy.ones(nspec, dtype=bool) if select is None else select
        return wave, flux[_select,:], ivar[_select,:], sres[_select,:]

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
            - At least one line must have ``mode='f'``
            - All tied lines must be tied to a line with a correctly
              specified index.
            - Warnings will be provided for any line with a centroid
              that falls outside of the provided wavelength range.
            - The database must provide at least one valid line.

        Args:
            emldb (:class:`mangadap.par.emissionlinedb.EmissionLineDB`):
                Emission-line database.
            wave (array-like):
                Wavelength vector.
            check_par (:obj:`bool`, optional):
                Validate the provided parameters.

        Raises:
            TypeError:
                Raised if the provided object is not an instance of
                :class:`mangadap.par.emissionlinedb.EmissionLineDB`.
            ValueError:
                Raised if any line has a mode of `x` or if the database
                does not provide a valid definition for any templates.
            NameError:
                Raised if a defined profile type is not known.

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
        for i in range(emldb.size):
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
    def check_emission_line_database(emldb, wave=None):
        r"""
        Check the emission-line database.  Modes are checked by
        :class:`mangadap.par.emissionlinedb.EmissionLinePar`, and the
        indices are checked to be unique by
        :class:`mangadap.par.emissionlinedb.EmissionLineDB`.

            - The type of the object must be
              :class:`mangadap.par.emissionlinedb.EmissionLineDB`
            - At least one line must be fit independently
            - All tied lines must be tied to a line with a correctly
              specified index.
            - The database must provide at least one valid line to fit.

        Args:
            emldb (:class:`mangadap.par.emissionlinedb.EmissionLineDB`):
                Emission-line database.
            wave (array-like):
                Wavelength vector.

        Raises:
            TypeError:
                Raised if the provided object is not an instance of
                :class:`mangadap.par.emissionlinedb.EmissionLineDB`.
            ValueError:
                Raised if there are no independently fitted lines, if there
                are no lines to fit, if there are no valid lines within the
                provided wavelength range, if the wavelength range is not
                provided as a 1D vector, or if the index to tie lines is not
                valid. (Warning: The latter does not catch lines with valid
                indices but where the tied line falls outside of the provided
                wavelength range.)
        """
        # Check the object type
        if not isinstance(emldb, EmissionLineDB):
            raise TypeError('Emission lines must be defined using an EmissionLineDB object.')

        # There must be one primary line
        if numpy.all(emldb['tie_index'] > 0):
            raise ValueError('At least one line in the database must be independent of all '
                             'other lines.  I.e., the tied index must be None.')

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
        for i in emldb['tie_index']:
            if i <= 0:
                continue
            if numpy.sum(emldb['index'] == i) == 0:
                raise ValueError('No line with index={0} to tie to!'.format(tied_index))

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
            ValueError:
                Raised if the length of the spectra, errors, or mask
                does not match the length of the wavelength array;
                raised if the wavelength, redshift, or dispersion arrays
                are not 1D vectors; and raised if the number of
                redshifts or dispersions is not a single value or the
                same as the number of input spectra.
        """

        # Check the input emission-line database
        if checkdb:
            EmissionLineFit.check_emission_line_database(emission_lines)

        nspec = flux.shape[0]
        nbands = len(emission_lines)
        if redshift is None:
            # If the redshift is NOT provided, use the fitted velocity
            # for each emission line.  Replace masked measurements with
            # mean (unmasked) redshift of that spectrum, or of all
            # spectra.
            _redshift = model_eml_par['KIN'][:,:,0]/astropy.constants.c.to('km/s').value
            _redshift = numpy.ma.MaskedArray(_redshift,
                                             mask=bitmask.flagged(model_eml_par['MASK']))
            mean_redshift = numpy.ma.mean(_redshift, axis=1).filled(numpy.ma.mean(_redshift))
            mean_redshift = numpy.array([mean_redshift]*nbands).T
            _redshift[_redshift.mask] = mean_redshift[_redshift.mask]
            _redshift = numpy.asarray(_redshift)
        else:
            # Use provided data and check its shape
            _redshift = redshift
            if _redshift.ndim == 1:
                if len(_redshift) != nspec:
                    raise ValueError('Must provide at least one redshift per input spectrum.')
                _redshift = numpy.array([_redshift]*nbands).T
            if _redshift.ndim == 2 and _redshift.shape != (nspec,nbands):
                raise ValueError('Provided redshift array does not match spectra and bands.')

        # Calculate the wavelength at which to measure the continuum,
        # matching what is done by
        # :class:`mangadap.proc.emissionlineMoments.EmissionLineMoments`
        line_center = (1+_redshift)*emission_lines['restwave'][None,:]

        # Compute the equivalent widths.  The checking done by
        # EmissionLineFit.check_and_prep_input is *identical* to what is
        # done within emission_line_equivalent_width()
        model_eml_par['BMED'], model_eml_par['RMED'], pos, model_eml_par['EWCONT'], \
                model_eml_par['EW'], model_eml_par['EWERR'] \
                        = emission_line_equivalent_width(wave, flux, emission_lines['blueside'],
                                                         emission_lines['redside'], line_center,
                                                         model_eml_par['FLUX'], mask=mask,
                                                         redshift=_redshift,
                                                         line_flux_err=model_eml_par['FLUXERR'])

        # Flag non-positive measurements
        if bitmask is not None:
            model_eml_par['MASK'][numpy.invert(pos)] \
                    = bitmask.turn_on(model_eml_par['MASK'][numpy.invert(pos)],
                                      'NON_POSITIVE_CONTINUUM')

    @staticmethod
    def line_metrics(emission_lines, wave, flux, ferr, model_flux, model_eml_par,
                     mask=None, model_mask=None, bitmask=None, window=15, fill_redshift=False):
        r"""
        Calculate fit-quality metrics near each emission line.

        .. todo::
            - Allow window to be defined in angstroms?

        Args:
            emission_lines (:class:`mangadap.par.emissionlinedb.EmissionLineDB`):
                Emission-line database use during the fit.
            wave (`numpy.ndarray`):
                Wavelength vector for object spectra.  Shape is
                :math:`(N_{\rm pix},)`.
            flux (`numpy.ndarray`):
                Object spectra that have been fit.  Can be provided as a
                `numpy.ma.MaskedArray`_.  Shape is :math:`(N_{\rm
                spec},N_{\rm pix})`.
            ferr (`numpy.ndarray`)
                :math:`1\sigma` errors in the object spectra.  Can be
                provided as a `numpy.ma.MaskedArray`_.  Shape is
                :math:`(N_{\rm spec},N_{\rm pix})`.
            model_flux (`numpy.ndarray`):
                Best-fitting model spectra. Can be provided as a
                `numpy.ma.MaskedArray`_.  Shape is (:math:`N_{\rm
                spec},N_{\rm pix}`).
            model_eml_par (:class:`EmissionLineFitDataTable`):
                A numpy record array with data type given by
                :func:`_per_emission_line_dtype`.  Uses ``FLUX``,
                ``KIN``, and ``MASK``; and assigns results to ``LINE_*``
                columns.
            mask (`numpy.ndarray`, optional):
                A mask for the object spectra that have been fit.  Added
                to mask attribute of `flux` if it is a
                `numpy.ma.MaskedArray`_.
            model_mask (`numpy.ndarray`, optional):
                A *boolean* numpy array with the mask for the model
                spectra.  Added to mask attribute of `model_flux` if it
                is a `numpy.ma.MaskedArray`_.
            bitmask (:class:`mangadap.util.bitmask.BitMask`, optional):
                The `BitMask` object used to interpret the MASK column
                in the `model_eml_par` object.  If None, the MASK column
                is ignored.
            window (:obj:`int`, optional):
                The width of the window used to compute the metrics
                around each line in number of pixels.
            fill_redshift (:obj:`bool`, optional):
                Fill any masked velocity measurement to the masked
                median of the velocities for the unmasked lines in the
                same spectrum when constructing the redshifted bands.
                If False, the A/N measurement is masked.

        Returns:
            `numpy.recarray`: Return the input `model_eml_par` after
            filling the LINE_PIXC, AMP, ANR, LINE_NSTAT, LINE_CHI2,
            LINE_RMS, and LINE_FRMS columns.

        Raises:
            ValueError:
                Raised if various checks of the input array sizes are
                incorrect.
        """
        # Check shapes and set masks for all arrays
        _wave = numpy.atleast_1d(wave)
        if _wave.ndim > 1:
            raise ValueError('Wavelength must be a 1D vector.')

        _flux = flux if isinstance(flux, numpy.ma.MaskedArray) else numpy.ma.MaskedArray(flux)
        _flux = numpy.ma.atleast_2d(_flux)
        if _flux.ndim > 2:
            raise ValueError('Flux must be 1 or 2D.')
        if _flux.shape[1] != len(_wave):
            raise ValueError('Flux has different spectral channels than the wavelength vector.')

        _ferr = ferr if isinstance(ferr, numpy.ma.MaskedArray) else numpy.ma.MaskedArray(ferr)
        _ferr = numpy.ma.atleast_2d(_ferr)
        if _ferr.shape != _flux.shape:
            raise ValueError('Flux and errors must have the same shape.')

        _mask = numpy.zeros_like(_flux, dtype=bool) if mask is None \
                        else numpy.atleast_2d(mask).astype(bool)
        if _mask.shape != _flux.shape:
            raise ValueError('Flux and mask must have the same shape.')
        _flux.mask |= _mask
        _ferr.mask |= _mask

        _model_flux = model_flux if isinstance(model_flux, numpy.ma.MaskedArray) \
                        else numpy.ma.MaskedArray(model_flux)
        _model_flux = numpy.ma.atleast_2d(_model_flux)
        if _flux.shape != _model_flux.shape:
            raise ValueError('Flux and model must have the same shape.')

        _model_mask = numpy.zeros_like(_flux, dtype=bool) if model_mask is None \
                        else numpy.atleast_2d(model_mask).astype(bool)
        if _model_mask.shape != _flux.shape:
            raise ValueError('Flux and model mask must have the same shape.')
        _model_flux.mask |= _model_mask

        # Spectra for figures of merit
        resid = numpy.square(_flux-_model_flux)
        fresid = numpy.square(numpy.ma.divide(_flux-_model_flux, _model_flux))
        chisqr = numpy.square(numpy.ma.divide(_flux-_model_flux, _ferr))
        spec_mask = resid.mask | fresid.mask | chisqr.mask
        resid[spec_mask] = numpy.ma.masked
        fresid[spec_mask] = numpy.ma.masked
        chisqr[spec_mask] = numpy.ma.masked

        # Get the pixel at which to center the metric calculations
        eml_mask = numpy.zeros(model_eml_par['MASK'].shape, dtype=bool) if bitmask is None else\
                        bitmask.flagged(model_eml_par['MASK'],
                                        flag=['INSUFFICIENT_DATA', 'FIT_FAILED', 'UNDEFINED_COVAR',
                                              'NEAR_BOUND'])
        z = numpy.ma.MaskedArray(model_eml_par['KIN'][:,:,0]/astropy.constants.c.to('km/s').value,
                                 mask=eml_mask)
        sample_wave = emission_lines['restwave'][None,:]*(1+z)
        interp = interp1d(wave, numpy.arange(wave.size), bounds_error=False, fill_value=-1,
                          assume_sorted=True)
        model_eml_par['LINE_PIXC'] = numpy.around(interp(sample_wave.data)).astype(int)
        model_eml_par['LINE_PIXC'][sample_wave.mask] = -1
        eml_mask = (model_eml_par['LINE_PIXC'] < 0) | (model_eml_par['LINE_PIXC'] >= wave.size)

        # Get the fitted line amplitude
        sigma_ang = model_eml_par['KIN'][:,:,1]*sample_wave/astropy.constants.c.to('km/s').value
        sigma_ang = numpy.ma.MaskedArray(sigma_ang, mask=numpy.invert(sigma_ang > 0))
        model_eml_par['AMP'] = numpy.ma.divide(model_eml_par['FLUX'], sigma_ang).filled(0.0) \
                                    / numpy.sqrt(2*numpy.pi)

        # Shift the bands to the appropriate redshift
        nspec = model_eml_par.size
        if numpy.any(z.mask) and fill_redshift:
            z_per_spec = numpy.ma.median(z, axis=1)
            for i in range(nspec):
                z[i,z.mask[i,:]] = z_per_spec[i]
        _bluebands = emission_lines['blueside'][None,:,:]*(1.0+z.data[:,:,None])
        _redbands = emission_lines['redside'][None,:,:]*(1.0+z.data[:,:,None])

        # Get the mean noise in the sidebands to either side of each
        # emission line
        noise = numpy.zeros_like(model_eml_par['AMP'], dtype=float)
        for i in range(nspec):
            _bluenoise = passband_median(_wave, _ferr[i,:], passband=_bluebands[i,:,:])
            _bluenoise = numpy.ma.MaskedArray(_bluenoise, mask=numpy.invert(_bluenoise > 0))
            _rednoise = passband_median(_wave, _ferr[i,:], passband=_redbands[i,:,:])
            _rednoise = numpy.ma.MaskedArray(_rednoise, mask=numpy.invert(_rednoise > 0))
            noise[i,:] = ((_bluenoise + _rednoise)/2.).filled(0.0)
        model_eml_par['ANR'] = numpy.ma.divide(model_eml_par['AMP'], noise).filled(0.0)

        neml = len(emission_lines)
        for i in range(neml):
            print('Getting fit metrics for line: {0}/{1}'.format(i+1, neml), end='\r')
            start = model_eml_par['LINE_PIXC'][:,i] - window//2
            end = start + window
            m = numpy.zeros(flux.shape, dtype=float)
            for j in range(nspec):
                if eml_mask[j,i]:
                    continue
                m[j,start[j]:end[j]] = 1.
            m[spec_mask] = 0.

            model_eml_par['LINE_NSTAT'][:,i] = numpy.sum(m,axis=1)
            model_eml_par['LINE_RMS'][:,i] = numpy.ma.sqrt(numpy.ma.mean(m*resid,axis=1)
                                                           ).filled(0.0)
            model_eml_par['LINE_FRMS'][:,i] = numpy.ma.sqrt(numpy.ma.mean(m*fresid,axis=1)
                                                            ).filled(0.0)
            model_eml_par['LINE_CHI2'][:,i] = numpy.ma.sum(m*chisqr,axis=1).filled(0.0)

        print('Getting fit metrics for line: {0}/{0}'.format(neml))
        return model_eml_par


