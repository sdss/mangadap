"""
A class hierarchy that measures moments of the observed emission lines.

----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import inspect
from pathlib import Path
import logging

import numpy

from astropy.io import fits
import astropy.constants

from ..util.fitsutil import DAPFitsUtil
from ..util.bitmask import BitMask
from ..util.dapbitmask import DAPBitMask
from ..util.datatable import DataTable
from ..util.pixelmask import SpectralPixelMask
from ..util.log import log_output
from ..par.parset import KeywordParSet
from ..par.artifactdb import ArtifactDB
from ..par.emissionmomentsdb import EmissionMomentsDB
from .spatiallybinnedspectra import SpatiallyBinnedSpectra
from .stellarcontinuummodel import StellarContinuumModel
from .bandpassfilter import passband_integral, passband_integrated_width
from .bandpassfilter import passband_weighted_mean, passband_weighted_sdev, pseudocontinuum
from .bandpassfilter import emission_line_equivalent_width
from .spectralfitting import EmissionLineFit


class EmissionLineMomentsDef(KeywordParSet):
    """
    A class that holds the parameters necessary to perform the
    emission-line moment measurements.

    The defined parameters are:

    .. include:: ../tables/emissionlinemomentsdef.rst
    """
    def __init__(self, key='EMOMMPL11', minimum_snr=0.0, pixelmask=None, passbands='ELBMPL9',
                 redo_postmodeling=True, fit_vel_name='Ha-6564', overwrite=False):

        # Use the signature to get the parameters and the default values
        sig = inspect.signature(self.__class__)
        pars = list(sig.parameters.keys())
        defaults = [sig.parameters[key].default for key in pars]

        # Remaining definitions done by hand
        in_fl = [int, float]

        values = [key, minimum_snr, pixelmask, passbands, redo_postmodeling, fit_vel_name,
                  overwrite]
        dtypes = [str, in_fl, SpectralPixelMask, str, bool, str, bool]
        descr = ['Keyword used to distinguish between different emission-line moment databases.',
                 'Minimum S/N of spectrum to analyze',
                 'Object used to mask spectral pixels',
                 'Either a string identifying the emission-line bandpass filter database to ' \
                    'use, or the direct path to the parameter file defining the database',
                 'Redo the moment measurements after the emission-line modeling has been ' \
                    'performed',
                 'The name of the emission line used to set the redshift of each spaxel used ' \
                    'to set the observed wavelength of the bandpasses.',
                 'If the output file already exists, redo all the calculations and overwrite it.']
        super().__init__(pars, values=values, defaults=defaults, dtypes=dtypes, descr=descr)
        self.momdb = None
        self._validate()

    @classmethod
    def from_dict(cls, d):
        """
        Instantiate from a dictionary.
        """
        # Copy over the primary keywords.
        # TODO: Get these from inspection of the instatiation signature!
        _d = {}
        for key in ['key', 'minimum_snr', 'passbands', 'redo_postmodeling', 'fit_vel_name',
                    'overwrite']:
            if key in d.keys():
                _d[key] = d[key]

        artifacts = None if 'artifact_mask' not in d.keys() or d['artifact_mask'] is None \
                        else ArtifactDB.from_key(d['artifact_mask'])
        waverange = None if 'waverange' not in d.keys() or d['waverange'] is None \
                        else d['waverange']

        if not all([a is None for a in [artifacts, waverange]]):
            _d['pixelmask'] = SpectralPixelMask(artdb=artifacts, waverange=waverange)

        # Return the instantiation
        return super().from_dict(_d)

    def _validate(self):
        """
        Validate the object and instantiate the passband database to use for the
        emission-line moments.
        """
        if self['passbands'] is None:
            raise ValueError('Must define a passband database.')
        # Check if the provided string is a file
        parfile = Path(self['passbands']).resolve()
        if parfile.exists():
            self.momdb = EmissionMomentsDB(str(parfile))
            return

        # Otherwise, assume this is a keyword selecting a database distributed
        # with the DAP has been chosen.
        self.momdb = EmissionMomentsDB.from_key(self['passbands'])


class EmissionLineMomentsBitMask(DAPBitMask):
    r"""
    Derived class that specifies the mask bits for the emission-line
    moment measurements.  The maskbits defined are:
    
    .. include:: ../tables/emissionlinemomentsbitmask.rst
    """
    cfg_root = 'emission_line_moments_bits'


class EmissionLineMomentsDataTable(DataTable):
    """
    Primary data table with the results of the emission-line-moment
    calculations.

    Table includes:

    .. include:: ../tables/emissionlinemomentsdatatable.rst

    Args:
        neml (:obj:`int`):
            Number of emission lines being measured.
        bitmask (:class:`~mangadap.util.bitmask.BitMask`, optional):
            Object used to flag mask bits. If None, flags are simply
            boolean.
        shape (:obj:`int`, :obj:`tuple`, optional):
            The shape of the initial array. If None, the data array
            will not be instantiated; use :func:`init` to initialize
            the data array after instantiation.
    """
    def __init__(self, neml=1, bitmask=None, shape=None):
        # NOTE: This should require python 3.7 to make sure that this
        # is an "ordered" dictionary.
        datamodel = dict(BINID=dict(typ=int, shape=None, descr='Spectrum/Bin ID number'),
                         BINID_INDEX=dict(typ=int, shape=None,
                                          descr='Index of the spectrum in the list of '
                                                'provided spectra.'),
                         REDSHIFT=dict(typ=float, shape=None,
                                       descr='Redshift used for shifting the passbands'),
                         MASK=dict(typ=bool if bitmask is None else bitmask.minimum_dtype(),
                                   shape=(neml,),
                                   descr='Bad-value boolean or bit mask value for the moments'),
                         BCEN=dict(typ=float, shape=(neml,), descr='Center of the blue sideband.'),
                         BCONT=dict(typ=float, shape=(neml,),
                                    descr='Pseudo-continuum in the blue sideband'),
                         BCONTERR=dict(typ=float, shape=(neml,),
                                       descr='Error in the blue-sideband pseudo-continuum'),
                         RCEN=dict(typ=float, shape=(neml,), descr='Center of the red sideband.'),
                         RCONT=dict(typ=float, shape=(neml,),
                                    descr='Pseudo-continuum in the red sideband'),
                         RCONTERR=dict(typ=float, shape=(neml,),
                                       descr='Error in the red-sideband pseudo-continuum'),
                         CNTSLOPE=dict(typ=float, shape=(neml,),
                                       descr='Continuum slope used to determine the continuum at '
                                             'the line center.'),
                         FLUX=dict(typ=float, shape=(neml,), descr='Summed flux (0th moment)'),
                         FLUXERR=dict(typ=float, shape=(neml,),
                                      descr='Error in the summed flux'),
                         MOM1=dict(typ=float, shape=(neml,),
                                   descr='Line centroid redshift (:math:`cz`; 1st moment)'),
                         MOM1ERR=dict(typ=float, shape=(neml,),
                                      descr='Error in the line centroid redshift'),
                         MOM2=dict(typ=float, shape=(neml,),
                                   descr='Line standard deviation (2nd moment)'),
                         MOM2ERR=dict(typ=float, shape=(neml,),
                                      descr='Error in the line standard deviation'),
                         SINST=dict(typ=float, shape=(neml,),
                                    descr='Instrumental dispersion at the line centroid'),
                         BMED=dict(typ=float, shape=(neml,),
                                   descr='Median flux in the blue sideband used for EW'),
                         RMED=dict(typ=float, shape=(neml,),
                                   descr='Median flux in the red sideband used for EW'),
                         EWCONT=dict(typ=float, shape=(neml,),
                                     descr='Continuum value used for EW calculation'),
                         EW=dict(typ=float, shape=(neml,),
                                 descr='Equivalent width (FLUX/pseudo continuum)'),
                         EWERR=dict(typ=float, shape=(neml,),
                                    descr='Error in the equivalent width'))

        keys = list(datamodel.keys())
        super(EmissionLineMomentsDataTable,
                self).__init__(keys, [datamodel[k]['typ'] for k in keys],
                               element_shapes=[datamodel[k]['shape'] for k in keys],
                               descr=[datamodel[k]['descr'] for k in keys],
                               shape=shape)


class EmissionLineMoments:
    r"""
    Class that holds the emission-line moment measurements.

    Args:
        loggers (list): (**Optional**) List of `logging.Logger`_ objects
            to log progress; ignored if quiet=True.  Logging is done
            using :func:`mangadap.util.log.log_output`.  Default is no
            logging.
        quiet (bool): (**Optional**) Suppress all terminal and logging
            output.  Default is False.

    Attributes:
        loggers (list): List of `logging.Logger`_ objects to log
            progress; ignored if quiet=True.  Logging is done using
            :func:`mangadap.util.log.log_output`.
        quiet (bool): Suppress all terminal and logging output.

    """
#    @profile
    def __init__(self, database, binned_spectra, stellar_continuum=None,
                 emission_line_model=None, redshift=None, database_list=None, artifact_path=None,
                 bandpass_path=None, output_path=None, output_file=None, hardcopy=True,
                 overwrite=None, checksum=False, loggers=None, quiet=False):

        self.loggers = None
        self.quiet = False

        # Define the database properties
#        self.database = None
#        self.pixelmask = None
#        self.artdb = None
#        self.momdb = None

        self.database = database
        if not isinstance(self.database, EmissionLineMomentsDef):
            raise TypeError('Method must have type EmissionLineMomentsDef.')

#        self._define_databases(database_key, database_list=database_list,
#                               artifact_path=artifact_path, bandpass_path=bandpass_path)

        self.nmom = self.database.momdb.size

        self.binned_spectra = None
        self.stellar_continuum = None
        self.emission_line_model = None
        self.redshift = None

        # Define the output directory and file
        self.directory_path = None      # Set in default_paths
        self.output_file = None
        self.hardcopy = None

        # Initialize the objects used in the assessments
        self.bitmask = EmissionLineMomentsBitMask()

        self.hdu = None
        self.checksum = checksum
        self.spatial_shape = None
        self.nspec = None
        self.spatial_index = None
        self.image_arrays = None
        self._assign_image_arrays()

        self.nbins = None
        self.missing_bins = None

        # Run the assessments of the DRP file
        self.measure(binned_spectra, stellar_continuum=stellar_continuum,
                     emission_line_model=emission_line_model, redshift=redshift, 
                     output_path=output_path, output_file=output_file, hardcopy=hardcopy,
                     overwrite=overwrite, loggers=loggers, quiet=quiet)

    def __getitem__(self, key):
        return self.hdu[key]

    @staticmethod
    def default_paths(cube, method_key, rdxqa_method, binning_method, stelcont_method=None,
                      elmodel_method=None, output_path=None, output_file=None):
        """
        Set the default directory and file name for the output file.

        Args:
            cube (:class:`mangadap.datacube.datacube.DataCube`):
                Datacube to analyze.
            method_key (:obj:`str`):
                Keyword designating the method used for the reduction
                assessments.
            rdxqa_method (:obj:`str`):
                The method key for the basic assessments of the datacube.
            binning_method (:obj:`str`):
                The method key for the spatial binning.
            stelcont_method (:obj:`str`, optional):
                The method key for the stellar-continuum fitting method.  If
                None, not included in output file name.
            elmodel_method (:obj:`str`, optional):
                The method key for the emission-line modeling method.  If None,
                not included in the output file name.
            output_path (:obj:`str`, `Path`_, optional):
                The path for the output file.  If None, the current working
                directory is used.
            output_file (:obj:`str`, optional):
                The name of the output "reference" file. The full path of the
                output file will be :attr:`directory_path`/:attr:`output_file`.
                If None, the default is to combine ``cube.output_root`` and the
                method keys.  The order of the keys is the order of operations
                (rdxqa, binning, stellar continuum, emission-line model).

        Returns:
            :obj:`tuple`: Returns a `Path`_ with the output directory and a
            :obj:`str` with the output file name.
        """
        directory_path = Path('.').resolve() if output_path is None \
                                else Path(output_path).resolve()
        method = f'{rdxqa_method}-{binning_method}'
        if stelcont_method is not None:
            method = f'{method}-{stelcont_method}'
        if elmodel_method is not None:
            method = f'{method}-{elmodel_method}'
        method = f'{method}-{method_key}'
        _output_file = f'{cube.output_root}-{method}.fits.gz' if output_file is None \
                            else output_file
        return directory_path, _output_file


    def _initialize_primary_header(self, hdr=None, measurements_binid=None):
        """
        Construct the primary header for the reference file.

        Args:
            hdr (`astropy.io.fits.Header`_, optional):
                Input base header for added keywords. If None, uses
                the :attr:`cube` header (if there is one) and then
                cleans the header using
                :func:`mangadap.util.fitsutil.DAPFitsUtil.clean_dap_primary_header`.
            measurements_binid (`numpy.ndarray`_, optional):
                Bin IDs for moment measurements. Only use is to check
                if this is None, and the boolean result is saved to
                the header to indicate if the emission-line moment
                measurements are disconnected from the
                stellar-continuum measurements.

        Returns:
            `astropy.io.fits.Header`_: Initialized header object.
        """
        # Copy the from the DRP and clean it
        if hdr is None:
            hdr = self.binned_spectra.cube.prihdr.copy()
            hdr = DAPFitsUtil.clean_dap_primary_header(hdr)
        
        hdr['AUTHOR'] = 'Kyle B. Westfall <westfall@ucolick.org>'
        hdr['ELMKEY'] = (self.database['key'], 'Emission-line moments database keyword')
        hdr['ELMMINSN'] = (self.database['minimum_snr'], 'Minimum S/N of spectrum to include')
#        hdr['ARTDB'] = (self.database['artifacts'], 'Artifact database keyword')
        hdr['MOMDB'] = (self.database['passbands'], 'Emission-line moments database keyword')
        if self.stellar_continuum is not None:
            hdr['SCKEY'] = (self.stellar_continuum.method['key'], 'Stellar-continuum model keyword')
        if self.emission_line_model is not None:
            hdr['ELFKEY'] = (self.emission_line_model.method['key'],
                                'Emission-line modeling method keyword')
        hdr['NBINS'] = (self.nbins, 'Number of unique spectra for measurements')
        hdr['ELMREBIN'] = (measurements_binid is not None, 'Bin IDs disconnected from SC binning')
        return hdr

    def _assign_image_arrays(self):
        """
        Set :attr:`image_arrays`, which contains the list of extensions
        in :attr:`hdu` that are on-sky image data.
        """
        self.image_arrays = [ 'BINID' ]

    def _get_missing_bins(self, unique_bins=None):
        if unique_bins is None:
            good_snr = self.binned_spectra.above_snr_limit(self.database['minimum_snr'])
            return numpy.sort(self.binned_spectra.missing_bins + 
                        self.binned_spectra['BINS'].data['BINID'][numpy.invert(good_snr)].tolist())
        return SpatiallyBinnedSpectra._get_missing_bins(unique_bins)


    def _assign_redshifts(self, redshift, measure_on_unbinned_spaxels, good_snr,
                          default_redshift=None):
        """
        Set the redshift to use for each spectrum for the emission-line
        moments.
        
        In terms of precedence, directly provided redshifts override
        those in any available EmissionLineModel.

        If :attr:`emission_line_model` and redshift are None, the default
        redshift is used (or 0.0 if this is also None).

        To get the emission-line stellar kinematics, the function calls
        :func:`mangadap.proc.emissionlinemodel.EmissionLineModel.matched_kinematics`.
        The method must have a fit_vel_name define for this use!

        In this function, the provided redshift must be a single value
        or None; therefore, the means of any vectors should be provided
        instead of the full vector.

        Args:
            redshift (:obj:`float`, `numpy.ndarray`_):
                Redshifts (:math:`z`) to use for each spectrum.
            measure_on_unbinned_spaxels (:obj:`bool`):
                Flag that method expects to measure moments on unbinned
                spaxels.
            good_snr (`numpy.ndarray`_):
                Boolean array setting which spectra have sufficient S/N
                for the measurements.

            default_redshift (:obj:`float`, optional):

                Only used if there are stellar kinematics available.
                Provides the default redshift to use for spectra without
                stellar measurements; see ``redshift`` in
                :func:`mangadap.proc.stellarcontinuummodel.StellarContinuumModel.matched_kinematics`.
                If None (default), the median of the unmasked stellar
                velocities will be used.

        """

        # Construct the binid matrix if measuring on unbinned spaxels
        if measure_on_unbinned_spaxels:
            binid = numpy.full(self.spatial_shape, -1, dtype=int)
            binid.ravel()[good_snr] = numpy.arange(self.nbins)
            missing = []
            nspec = self.binned_spectra.cube.nspec
        else:
            binid = self.binned_spectra['BINID'].data
            missing = self.binned_spectra.missing_bins
            nspec = self.binned_spectra.nbins

        #---------------------------------------------------------------
        # Get the redshifts measured for the selected emission line and
        # use them if no default value is provided
        if self.emission_line_model is not None and self.database['fit_vel_name'] is not None \
                and redshift is None:
            self.redshift, _ = self.emission_line_model.matched_kinematics(
                                        binid, self.database['fit_vel_name'],
                                        redshift=default_redshift, nearest=True, missing=missing)
            if measure_on_unbinned_spaxels:
                tmp = self.redshift.copy()
                self.redshift = numpy.zeros(nspec, dtype=float)
                self.redshift[good_snr] = tmp
            return

        #---------------------------------------------------------------
        # Use the default value(s)
        _redshift = numpy.atleast_1d(redshift)
        if len(_redshift) not in [ 1, nspec ]:
            raise ValueError('Provided redshift must be either a single value or match the '
                             'number of binned spectra or the number of unbinned spaxels.')
        self.redshift = numpy.full(nspec, redshift, dtype=float) \
                                if len(_redshift) == 1 else _redshift.copy()

    def _flag_good_spectra(self, measure_on_unbinned_spaxels):
        if measure_on_unbinned_spaxels:
            return self.binned_spectra.check_fgoodpix()
        return self.binned_spectra.above_snr_limit(self.database['minimum_snr'])

    @staticmethod
    def sideband_pseudocontinua(wave, spec, sidebands, spec_n=None, err=None, log=True):
        """
        Get the side-band integrals.

        err is a single vector that is the same for all spec!

        Returns:
            :obj:`float`, `numpy.ndarray`_: Return five arrays or
            floats: (1) The center of each passband, (2) the mean
            continuum level, (3) the propagated error in the
            continuum level (will be None if no errors are provided),
            (4) flag that part of the passband was masked, (5) flag
            that the passband was fully masked or empty.
        """

        if spec.ndim == 1:
            return pseudocontinuum(wave, spec, passband=sidebands, err=err, log=log)

        # Calculate the pseudo-continua in the sidebands
        nspec = spec.shape[0]
        nbands = sidebands.shape[0]

        # More than one spectrum entered, use spec_n to decide which
        # spectrum to use for the integrals.  If spec_n is None and
        # there is more than one spectrum (improper use of the code!),
        # allow the code to continue but issue a warning
        if spec_n is None:
            warnings.warn('Multiple spectra entered, but no designation for which to use!')
            _spec_n = numpy.zeros(nbands, dtype=int)
        else:
            _spec_n = spec_n
        
        # TODO: Change this to a warning?
        if numpy.amax(_spec_n) > nspec-1:
            raise ValueError('Selected spectrum not provided!')

        # Instantiate the arrays
        continuum = numpy.zeros(nbands, dtype=float)
        continuum_error = None if err is None else numpy.zeros(nbands, dtype=float)
        band_center = numpy.zeros(nbands, dtype=float)
        band_incomplete = numpy.zeros(nbands, dtype=bool)
        band_empty = numpy.zeros(nbands, dtype=bool)

        for i in range(nspec):
            indx = _spec_n==i
            if numpy.sum(indx) == 0:
                continue
            band_center[indx], continuum[indx], _continuum_error, band_incomplete[indx], \
                    band_empty[indx] = pseudocontinuum(wave, spec[i,:], passband=sidebands[indx,:],
                                                       err=err, log=log)
            if err is not None:
                continuum_error[indx] = _continuum_error
        return band_center, continuum, continuum_error, band_incomplete, band_empty

    @staticmethod
    def single_band_moments(wave, spec, passband, restwave, err=None):
        """
        Measure the moments for a single band.

        If spectrum errors are not provided, moment errors are returned
        as 0.0.

        Args:
            wave (`numpy.ndarray`_):
                Wavelengths in angstroms
            spec (`numpy.ndarray`_):
                Flux in flux density (per angstrom)
            restwave (float):
                The rest wavelength of the line.
            err (`numpy.ndarray`_):
                The 1-sigma errors in the flux density.

        Returns:
            flux
            fluxerr
            mom1
            mom1err
            mom2
            mom2err
            incomplete band
            empty band
            division by zero
            undefined 1st moment
            undefined 2nd moment

        """
        _passband = numpy.atleast_1d(passband)
        if _passband.ndim > 1:
            raise ValueError('Should only pass a single band to single_band_moments.')
        if len(_passband) != 2:
            raise ValueError('Band should be a list of 2 numbers.')
        # Get the fraction of the passband that is unmasked
        interval_frac = passband_integrated_width(wave, spec, passband=passband, log=True) \
                                / numpy.diff(passband)[0]
        incomplete = interval_frac < 1.0
        empty = numpy.invert(interval_frac > 0.0)
        if empty:
            return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, incomplete, empty, False, True, True

        # Get the redshift vector
        cz = astropy.constants.c.to('km/s').value*(wave/restwave-1)

        # Get the integrated flux
        flux = passband_integral(wave, spec, passband=passband, log=True)
        fluxerr = 0.0 if err is None else \
                    passband_integral(wave, err, passband=passband, log=True, quad=True)

        # No flux in the passband, meaning division by zero will occur!
        if not numpy.absolute(flux) > 0.0:
            return flux, fluxerr, 0.0, 0.0, 0.0, 0.0, incomplete, empty, True, True, True

        # Get the first and second moments:
        mom1, mom1err = passband_weighted_mean(wave, spec, cz, yerr=err, passband=passband,
                                               log=True)
        mom2, mom2err = passband_weighted_sdev(wave, spec, cz, yerr=err, passband=passband,
                                               log=True)

        # If no errors provided, set errors to 0
        if err is None:
            mom1err = 0.0
            mom2err = 0.0

        # Set masked values to 0.0
        undefined_mom1 = False
        if isinstance(mom1, numpy.ma.core.MaskedConstant) and mom1.mask:
            mom1 = 0.0
            mom1err = 0.0
            undefined_mom1 = True
        undefined_mom2 = False
        if isinstance(mom2, numpy.ma.core.MaskedConstant) and mom2.mask:
            mom2 = 0.0
            mom2err = 0.0
            undefined_mom2 = True
        if isinstance(mom1err, numpy.ma.core.MaskedConstant) and mom1err.mask:
            mom1err = 0.0
        if isinstance(mom2err, numpy.ma.core.MaskedConstant) and mom2err.mask:
            mom2err = 0.0

        return flux, fluxerr, mom1, mom1err, mom2, mom2err, incomplete, empty, False, \
                    undefined_mom1, undefined_mom2

    @staticmethod
    def continuum_subtracted_moments(wave, spec, mainbands, restwave, bcen, bcont, rcen, rcont,
                                     spec_n=None, err=None, sres=None):
        """
        Calculate the continuum-subtracted moments
        
        err and sres are single vectors that are the same for all spec!

        Returns:
            cntm
            cntb
            flux
            fluxerr
            mom1
            mom1err
            mom2
            mom2err
            sinst
            incomplete
            empty
            divbyzero
            undefined_mom2

        """
        # Get the parameters for the linear continuum across the
        # primary passband
        nbands = mainbands.shape[0]
        flux = numpy.zeros(nbands, dtype=float)
        fluxerr = numpy.zeros(nbands, dtype=float)
        mom1 = numpy.zeros(nbands, dtype=float)
        mom1err = numpy.zeros(nbands, dtype=float)
        mom2 = numpy.zeros(nbands, dtype=float)
        mom2err = numpy.zeros(nbands, dtype=float)
        incomplete = numpy.zeros(nbands, dtype=bool)
        empty = numpy.zeros(nbands, dtype=bool)
        divbyzero = numpy.zeros(nbands, dtype=bool)
        undefined_mom1 = numpy.zeros(nbands, dtype=bool)
        undefined_mom2 = numpy.zeros(nbands, dtype=bool)

        if numpy.any(rcen-bcen == 0):
            print('band centers are the same!')
        cntm = (rcont - bcont) / (rcen - bcen)
        cntb = bcont - bcen*cntm

        # If more than one spectrum entered, use spec_n to decide which
        # spectrum to use for the integrals
        if len(spec.shape) > 1:
            # In this case, spec_n should always be provided, but allow
            # the code to continue with a warning
            if spec_n is None:
                warnings.warn('Multiple spectra entered, but no designation for which to use!')
                _spec_n = numpy.zeros(self.nmom, dtype=int)
            else:
                _spec_n = spec_n

        # Get the moments for all the main passbands
        for i,p in enumerate(mainbands):
            cont = cntb[i] + cntm[i]*wave
            _spec = spec - cont if len(spec.shape) == 1 else spec[_spec_n[i],:] - cont

            flux[i], fluxerr[i], mom1[i], mom1err[i], mom2[i], mom2err[i], incomplete[i], \
                        empty[i], divbyzero[i], undefined_mom1[i], undefined_mom2[i] \
                        = EmissionLineMoments.single_band_moments(wave, _spec, p, restwave[i],
                                                                  err=err)

        # Get the instrumental dispersion at the position of the line
        # center
        sinst = numpy.zeros(nbands, dtype=float) if sres is None \
                        else EmissionLineFit.instrumental_dispersion(wave, sres, restwave, mom1)

        return cntm, cntb, flux, fluxerr, mom1, mom1err, mom2, mom2err, sinst, incomplete, empty, \
                    divbyzero, undefined_mom1, undefined_mom2

    @staticmethod
    def measure_moments(momdb, wave, flux, ivar=None, mask=None, sres=None, continuum=None,
                        redshift=None, bitmask=None):
        """
        Measure the emission-line moments.

        If not input as masked arrays, flux is converted to one.
        """
        # Check the input moment database
        if not isinstance(momdb, EmissionMomentsDB):
            raise TypeError('Input database must have type EmissionMomentsDB.')

        # Check the bitmask if provided
        if bitmask is not None and not isinstance(bitmask, BitMask):
            raise TypeError('Input bitmask must have type BitMask.')

        # Prepare the spectra for the measurements
        _flux, _err, _sres, _continuum, _redshift, _ \
                    = EmissionLineFit.check_and_prep_input(wave, flux, ivar=ivar, mask=mask,
                                                           sres=sres, continuum=continuum,
                                                           redshift=redshift)

        # Subtract the continuum and determine where the continuum is
        # undefined.  If _continuum is None, this just returns _fluxnc =
        # _flux and no_continuum = None.
        _fluxnc, no_continuum = EmissionLineFit.subtract_continuum(_flux, _continuum)

        # Initialize the output data
        nspec = _flux.shape[0]
        nmom = momdb.size
        measurements = EmissionLineMomentsDataTable(neml=nmom, bitmask=bitmask, shape=nspec)

        # Save the redshift used
        measurements['REDSHIFT'] = _redshift

        # Mask dummy lines
        if numpy.any(momdb.dummy):
            measurements['MASK'][:,momdb.dummy] = True if bitmask is None else \
                    bitmask.turn_on(measurements['MASK'][:,momdb.dummy], 'UNDEFINED_BANDS')
        
        # No valid lines
        if numpy.all(momdb.dummy):
            return measurements

        # Common arrays used for each spectrum
        blue_fraction = numpy.zeros(nmom, dtype=float)
        red_fraction = numpy.zeros(nmom, dtype=float)

        incomplete = numpy.zeros(nmom, dtype=bool)
        empty = numpy.zeros(nmom, dtype=bool)
        divbyzero = numpy.zeros(nmom, dtype=bool)
        undefined_mom1 = numpy.zeros(nmom, dtype=bool)
        undefined_mom2 = numpy.zeros(nmom, dtype=bool)

        # Perform the measurements for each spectrum
        for i in range(nspec):
            print('Measuring emission-line moments in spectrum: {0}/{1}'.format(i+1,nspec),
                  end='\r')

            incomplete[:] = False
            empty[:] = False
            divbyzero[:] = False
            undefined_mom1[:] = False
            undefined_mom2[:] = False

            # Shift the sidebands to the appropriate redshift
            _bluebands = momdb['blueside'][numpy.invert(momdb.dummy)]*(1.0+_redshift[i])
            _redbands = momdb['redside'][numpy.invert(momdb.dummy)]*(1.0+_redshift[i])
            __err = None if _err is None else _err[i,:]
            __sres = None if _sres is None else _sres[i,:]

            # Check if the model was subtracted over the full range of
            # the sidebands, or if the model was subtracted in one
            # sideband but not the other
            if no_continuum is None:
                spec = _flux[i,:]
                spec_n = None

                # Flag moment calculations without a stellar-absorption
                # correction
                measurements['MASK'][i,:] = True if bitmask is None \
                                                else bitmask.turn_on(measurements['MASK'][i,:],
                                                                     'NO_ABSORPTION_CORRECTION')
            else:
                blue_fraction[numpy.invert(momdb.dummy)] \
                        = passband_integrated_width(wave, no_continuum[i,:], passband=_bluebands,
                                                log=True) / numpy.diff(_bluebands, axis=1).ravel()
                indx = numpy.logical_and(blue_fraction > 0.0, blue_fraction < 1.0)
                measurements['MASK'][i,indx] = True if bitmask is None \
                            else bitmask.turn_on(measurements['MASK'][i,indx], 'BLUE_JUMP')

                red_fraction[numpy.invert(momdb.dummy)] \
                        = passband_integrated_width(wave, no_continuum[i,:], passband=_redbands,
                                                log=True) / numpy.diff(_redbands, axis=1).ravel()
                indx = numpy.logical_and(red_fraction > 0.0, red_fraction < 1.0)
                measurements['MASK'][i,indx] = True if bitmask is None \
                            else bitmask.turn_on(measurements['MASK'][i,indx], 'RED_JUMP')
                
                indx = (numpy.invert(blue_fraction > 0) & (red_fraction > 0)) \
                            | (numpy.invert(red_fraction > 0) & (blue_fraction > 0))
                measurements['MASK'][i,indx] = True if bitmask is None \
                            else bitmask.turn_on(measurements['MASK'][i,indx],'JUMP_BTWN_SIDEBANDS')

                # Only use model-subtracted measurements if there are no
                # jumps in the bands or between the bands
                jumps = measurements['MASK'][i,:] if bitmask is None \
                            else bitmask.flagged(measurements['MASK'][i,:],
                                                 flag=['BLUE_JUMP', 'RED_JUMP',
                                                       'JUMP_BTWN_SIDEBANDS'])

                spec = numpy.ma.append(_flux[i,:].reshape(1,-1), _fluxnc[i,:].reshape(1,-1), axis=0)
                spec_n = numpy.zeros(nmom, dtype=int)
                spec_n[numpy.invert(jumps)] = 1

                # Flag moment calculations as having not been corrected
                # for stellar absorption
                measurements['MASK'][i,spec_n == 0] = True if bitmask is None \
                            else bitmask.turn_on(measurements['MASK'][i,spec_n == 0],
                                                 'NO_ABSORPTION_CORRECTION')

            # Get the blue pseudo continuum
            measurements['BCEN'][i,numpy.invert(momdb.dummy)], \
                measurements['BCONT'][i,numpy.invert(momdb.dummy)], conterr, \
                incomplete[numpy.invert(momdb.dummy)], empty[numpy.invert(momdb.dummy)] \
                    = EmissionLineMoments.sideband_pseudocontinua(wave, spec, _bluebands,
                                                                  spec_n=None if spec_n is None
                                                            else spec_n[numpy.invert(momdb.dummy)],
                                                                  err=__err)
            if __err is not None:
                measurements['BCONTERR'][i,numpy.invert(momdb.dummy)] = conterr
            if numpy.sum(incomplete) > 0:
                measurements['MASK'][i,incomplete] = True if bitmask is None \
                            else bitmask.turn_on(measurements['MASK'][i,incomplete], 'BLUE_INCOMP')
            if numpy.sum(empty) > 0:
                measurements['MASK'][i,empty] = True if bitmask is None \
                            else bitmask.turn_on(measurements['MASK'][i,empty], 'BLUE_EMPTY')

            # Get the red pseudo continuum
            measurements['RCEN'][i,numpy.invert(momdb.dummy)], \
                measurements['RCONT'][i,numpy.invert(momdb.dummy)], conterr, \
                incomplete[numpy.invert(momdb.dummy)], empty[numpy.invert(momdb.dummy)] \
                    = EmissionLineMoments.sideband_pseudocontinua(wave, spec, _redbands,
                                                                  spec_n=None if spec_n is None
                                                            else spec_n[numpy.invert(momdb.dummy)],
                                                                  err=__err)
            if __err is not None:
                measurements['RCONTERR'][i,numpy.invert(momdb.dummy)] = conterr
            if numpy.sum(incomplete) > 0:
                measurements['MASK'][i,incomplete] = True if bitmask is None \
                            else bitmask.turn_on(measurements['MASK'][i,incomplete], 'RED_INCOMP')
            if numpy.sum(empty) > 0:
                measurements['MASK'][i,empty] = True if bitmask is None \
                            else bitmask.turn_on(measurements['MASK'][i,empty], 'RED_EMPTY')

            # Get the main passband integral after subtracting the
            # continuum
            indx = numpy.invert(momdb.dummy) \
                        & numpy.invert(measurements['MASK'][i,:] if bitmask is None
                                        else bitmask.flagged(measurements['MASK'][i,:],
                                                             flag=['BLUE_EMPTY', 'RED_EMPTY',
                                                                   'BLUE_JUMP', 'RED_JUMP']))

            _mainbands = momdb['primary'][indx]*(1.0+_redshift[i])

            measurements['CNTSLOPE'][i,indx], cntb, measurements['FLUX'][i,indx], \
                    measurements['FLUXERR'][i,indx], measurements['MOM1'][i,indx], \
                    measurements['MOM1ERR'][i,indx], measurements['MOM2'][i,indx], \
                    measurements['MOM2ERR'][i,indx], measurements['SINST'][i,indx], \
                    incomplete[indx], empty[indx], divbyzero[indx], undefined_mom1[indx], \
                    undefined_mom2[indx] = \
                            EmissionLineMoments.continuum_subtracted_moments(wave, spec,
                                                               _mainbands, momdb['restwave'][indx],
                                                               measurements['BCEN'][i,indx],
                                                               measurements['BCONT'][i,indx],
                                                               measurements['RCEN'][i,indx],
                                                               measurements['RCONT'][i,indx],
                                                               spec_n=None if spec_n is None 
                                                                    else spec_n[indx], err=__err,
                                                               sres=__sres) 

            # Turn on any necessary bits
            measurements['MASK'][i,incomplete] = True if bitmask is None \
                        else bitmask.turn_on(measurements['MASK'][i,incomplete], 'MAIN_INCOMP')
            measurements['MASK'][i,empty] = True if bitmask is None \
                        else bitmask.turn_on(measurements['MASK'][i,empty], 'MAIN_EMPTY')
            measurements['MASK'][i,divbyzero] = True if bitmask is None \
                        else bitmask.turn_on(measurements['MASK'][i,divbyzero], 'DIVBYZERO')
            measurements['MASK'][i,undefined_mom1] = True if bitmask is None \
                        else bitmask.turn_on(measurements['MASK'][i,undefined_mom1],
                                             'UNDEFINED_MOM1')
            measurements['MASK'][i,undefined_mom2] = True if bitmask is None \
                        else bitmask.turn_on(measurements['MASK'][i,undefined_mom2],
                                             'UNDEFINED_MOM2')

        print('Measuring emission-line moments in spectrum: {0}/{0}'.format(nspec))
        return measurements
        
    def file_name(self):
        """Return the name of the output file."""
        return self.output_file

    def file_path(self):
        """Return the full path to the output file."""
        if self.directory_path is None or self.output_file is None:
            return None
        return self.directory_path / self.output_file

    def measure(self, binned_spectra, stellar_continuum=None, emission_line_model=None,
                redshift=None, output_path=None, output_file=None, hardcopy=True, overwrite=None,
                loggers=None, quiet=False):
        """
        Measure the emission-line moments using the binned spectra.

        If neither stellar-continuum nor emission-line models are
        provided:

            - Moments are measure on the binned spectra
            - No continuum subtraction is performed

        If a stellar-continuum model is provided without an
        emission-line model:

            - Moments are measured on the binned spectra
            - The best-fitting stellar continuum is subtracted

        If an emission-line model is provided without a
        stellar-continuum model:

            - Moments are measured on the relevant (binned or unbinned)
              spectra
            - If the emission-line model includes data regarding the
              stellar-continuum fit (template spectra and template
              weights), the continuum is subtracted before the
              measurements are made; otherwise, no continuum is
              subtracted

        If both stellar-continuum and emission-line models are provided,
        and if the stellar-continuum and emission-line fits are
        performed on the same spectra:

            - Moments are measured on the relevant (binned or unbinned)
              spectra
            - If the emission-line model includes data regarding the
              stellar-continuum fit (template spectra and template
              weights), the continuum is subtracted before the
              measurements are made; otherwise, the stellar continuum is
              used.

        If both stellar-continuum and emission-line models are provided,
        and if the stellar-continuum and emission-line fits are
        performed on different spectra:

            - The behavior is exactly as if the stellar-continuum model
              was not provided.

        Args:
            output_path (:obj:`str`, `Path`_, optional):
                The path for the output file.  If None, the current working
                directory is used.
            output_file (:obj:`str`, optional):
                The name of the output "reference" file. The full path of the
                output file will be :attr:`directory_path`/:attr:`output_file`.
                If None, the default is to combine ``cube.output_root`` and the
                method keys.  The order of the keys is the order of operations
                (rdxqa, binning, stellar continuum, emission-line model).

        """
        # Initialize the reporting
        if loggers is not None:
            self.loggers = loggers
        self.quiet = quiet
        
        # SpatiallyBinnedSpectra object always needed
        if binned_spectra is None:
            raise ValueError('Must provide spectra object for fitting.')
        if not isinstance(binned_spectra, SpatiallyBinnedSpectra):
            raise TypeError('Must provide a valid SpatiallyBinnedSpectra object!')
        if binned_spectra.hdu is None:
            raise ValueError('Provided SpatiallyBinnedSpectra object is undefined!')
        self.binned_spectra = binned_spectra

        # Check stellar-continuum model object, if provided
        self.stellar_continuum = None
        if stellar_continuum is not None:
            if not isinstance(stellar_continuum, StellarContinuumModel):
                raise TypeError('Provided stellar continuum must have StellarContinuumModel type!')
            if stellar_continuum.hdu is None:
                raise ValueError('Provided StellarContinuumModel is undefined!')
            self.stellar_continuum = stellar_continuum

        # Check emission-line model object, if provided
        self.emission_line_model = None
        if emission_line_model is not None:
#            if not isinstance(emission_line_model, EmissionLineModel):
#                raise TypeError('Provided emission line models must be of type EmissionLineModel.')
            if emission_line_model.hdu is None:
                raise ValueError('Provided EmissionLineModel is undefined!')
            self.emission_line_model = emission_line_model

        # What stellar continuum is available?
        #  - Assume the stellar continuum can always be extracted if a
        #    StellarContinuumModel object is provided
        #  - Check if the EmissionLineModel fitter has the appropriate
        #    function; True for Sasuke, not currently true for Elric
        eml_stellar_continuum_available = self.emission_line_model is not None \
                and callable(self.emission_line_model.method['fitclass'].construct_continuum_models)

        # What spectra to use?
        #  - Assume StellarContinuumModel always fits the binned spectra
        #  - The EmissionLineModel fits the binned spectra or unbinned
        #    spaxels as specified by its deconstruct_bins flag
        measure_on_unbinned_spaxels = self.emission_line_model is not None \
                and self.emission_line_model.method['fitpar']['deconstruct_bins'] != 'ignore'

        self.spatial_shape =self.binned_spectra.spatial_shape
        self.nspec = self.binned_spectra.cube.nspec if measure_on_unbinned_spaxels \
                            else self.binned_spectra.nbins
        self.spatial_index = self.binned_spectra.spatial_index.copy()
        
        #---------------------------------------------------------------
        # Get the good spectra
        good_snr = self._flag_good_spectra(measure_on_unbinned_spaxels)

        # Set the number of bins measured and missing bins
        self.nbins = numpy.sum(good_snr)
        self.missing_bins = [] if measure_on_unbinned_spaxels else self._get_missing_bins()

        # Get the redshifts to apply
        self._assign_redshifts(redshift, measure_on_unbinned_spaxels, good_snr)

        # Report
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, '-'*50)
            log_output(loggers, 1, logging.INFO, f'{"EMISSION-LINE MOMENTS":^50}')
            log_output(self.loggers, 1, logging.INFO, '-'*50)
            log_output(self.loggers, 1, logging.INFO, '-'*50)
            log_output(self.loggers, 1, logging.INFO, 'Measurements for {0}'.format(
                        'unbinned spaxels' if measure_on_unbinned_spaxels else 'binned spectra'))
            if self.stellar_continuum is None and self.emission_line_model is None:
                log_output(self.loggers, 1, logging.INFO, 'No stellar continuum available.')
            elif self.emission_line_model is not None and eml_stellar_continuum_available:
                log_output(self.loggers, 1, logging.INFO,
                           'Using stellar continuum from emission-line model fit.')
            else:
                log_output(self.loggers, 1, logging.INFO,
                           'Using stellar continuum from stellar-continuum model fit.')
            log_output(self.loggers, 1, logging.INFO, f'Number of spectra: {self.nspec}')
            if not measure_on_unbinned_spaxels and len(self.binned_spectra.missing_bins) > 0:
                log_output(self.loggers, 1, logging.INFO,
                           f'Missing bins: {len(self.binned_spectra.missing_bins)}')
            log_output(self.loggers, 1, logging.INFO,
                       f'With good S/N and to measure: {numpy.sum(good_snr)}')

        # Make sure there are good spectra
        if numpy.sum(good_snr) == 0:
            raise ValueError('No good spectra for measurements!')

        #---------------------------------------------------------------
        # (Re)Set the output paths
        self.directory_path, self.output_file = \
                EmissionLineMoments.default_paths(self.binned_spectra.cube, self.database['key'],
                        self.binned_spectra.rdxqa.method['key'], self.binned_spectra.method['key'],
                        stelcont_method=None if self.stellar_continuum is None 
                                        else self.stellar_continuum.method['key'],
                        elmodel_method=None if self.emission_line_model is None 
                                       else self.emission_line_model.method['key'],
                        output_path=output_path, output_file=output_file)
#        self._set_paths(directory_path, dapver, analysis_path, output_file)

        _overwrite = self.database['overwrite'] if overwrite is None else overwrite

        #---------------------------------------------------------------
        # Check that the file path is defined
        ofile = self.file_path()
        if ofile is None:
            raise ValueError('File path for output file is undefined!')

        # Report
        if not self.quiet:
            log_output(self.loggers,1,logging.INFO, f'Output path: {self.directory_path}')
            log_output(self.loggers,1,logging.INFO, f'Output file: {self.output_file}')
        
        #---------------------------------------------------------------
        # If the file already exists, and not overwriting, just read the
        # file
        if ofile.exists() and not _overwrite:
            self.hardcopy = True
            if not self.quiet:
                log_output(self.loggers, 1, logging.INFO, 'Using existing file')
            self.read(checksum=self.checksum)
            if not self.quiet:
                log_output(self.loggers, 1, logging.INFO, '-'*50)
            return

        #---------------------------------------------------------------
        # Get the spectra to use for the measurements
        wave, flux, ivar, sres \
                = EmissionLineFit.get_spectra_to_fit(self.binned_spectra,
                                                     pixelmask=self.database['pixelmask'],
                                                     select=good_snr,
                                                     original_spaxels=measure_on_unbinned_spaxels)

        # Get the continuum if there is any
        if eml_stellar_continuum_available:
            if measure_on_unbinned_spaxels:
                binid = numpy.full(self.binned_spectra.spatial_shape, -1, dtype=int)
                binid.ravel()[good_snr] = numpy.arange(self.nbins)
                continuum = self.emission_line_model.fill_continuum_to_match(binid)
            else:
                continuum = self.emission_line_model.fill_continuum_to_match(
                                                self.binned_spectra['BINID'].data,
                                                missing=self.binned_spectra.missing_bins)
        elif self.stellar_continuum is not None:
            continuum = self.stellar_continuum.fill_to_match(self.binned_spectra['BINID'].data,
                                                        missing=self.binned_spectra.missing_bins)

        #---------------------------------------------------------------
        # Perform the moment measurements
        _redshift = self.redshift[good_snr]
        measurements = self.measure_moments(self.database.momdb, wave, flux, ivar=ivar, sres=sres,
                                            continuum=continuum, redshift=_redshift,
                                            bitmask=self.bitmask)

        #---------------------------------------------------------------
        # Perform the equivalent width measurements
        include_band = numpy.array([numpy.invert(self.database.momdb.dummy)]*self.nbins) \
                        & numpy.invert(self.bitmask.flagged(measurements['MASK'],
                                                            flag=['BLUE_EMPTY', 'RED_EMPTY']))

        # Set the line center at the center of the primary passband
        line_center = (1.0+_redshift)[:,None]*self.database.momdb['restwave'][None,:]
        measurements['BMED'], measurements['RMED'], pos, measurements['EWCONT'], \
                measurements['EW'], measurements['EWERR'] \
                        = emission_line_equivalent_width(wave, flux,
                                                         self.database.momdb['blueside'],
                                                         self.database.momdb['redside'],
                                                         line_center, measurements['FLUX'],
                                                         redshift=_redshift,
                                                         line_flux_err=measurements['FLUXERR'],
                                                         include_band=include_band)

        # Flag non-positive measurements
        measurements['MASK'][include_band & numpy.invert(pos)] \
                    = self.bitmask.turn_on(measurements['MASK'][include_band & numpy.invert(pos)],
                                           'NON_POSITIVE_CONTINUUM')

        #---------------------------------------------------------------
        # Set the number of bins measured, missing bins, and bin IDs
        measurements['BINID'] = numpy.arange(self.nbins) if measure_on_unbinned_spaxels \
                                    else self.binned_spectra['BINS'].data['BINID'][good_snr]
        measurements['BINID_INDEX'] = numpy.arange(self.nbins) if measure_on_unbinned_spaxels \
                                        else numpy.arange(self.binned_spectra.nbins)[good_snr]
        measurements_binid = None
        if measure_on_unbinned_spaxels:
            measurements_binid = numpy.full(self.binned_spectra.spatial_shape, -1, dtype=int)
            measurements_binid.ravel()[good_snr] = numpy.arange(self.nbins)

        #---------------------------------------------------------------
        # Initialize the header
        self.hardcopy = hardcopy
        pri_hdr = self._initialize_primary_header(measurements_binid=measurements_binid)
        map_hdr = DAPFitsUtil.build_map_header(self.binned_spectra.cube.fluxhdr,
                                               'K Westfall <westfall@ucolick.org>')

        # Get the spatial map mask
        map_mask = numpy.zeros(self.spatial_shape, dtype=self.bitmask.minimum_dtype())

        # Account for measurements on individual spaxels
        if measurements_binid is None:
            # Add any spaxel not used because it was flagged by the
            # binning step
            indx = self.binned_spectra['MAPMASK'].data > 0
            map_mask[indx] = self.bitmask.turn_on(map_mask[indx], 'DIDNOTUSE')
            # Isolate any spaxels with foreground stars
            indx = self.binned_spectra.bitmask.flagged(self.binned_spectra['MAPMASK'].data,
                                                       'FORESTAR')
            map_mask[indx] = self.bitmask.turn_on(map_mask[indx], 'FORESTAR')
            # Get the bins that were below the S/N limit
            indx = numpy.invert(DAPFitsUtil.reconstruct_map(self.spatial_shape,
                                                        self.binned_spectra['BINID'].data.ravel(),
                                                        good_snr, dtype='bool')) & (map_mask == 0)
            map_mask[indx] = self.bitmask.turn_on(map_mask[indx], 'LOW_SNR')

            # Get the bin ids with measured indices
            bin_indx = DAPFitsUtil.downselect_bins(self.binned_spectra['BINID'].data.ravel(),
                                                measurements['BINID']).reshape(self.spatial_shape)
        else:
            # Assume any model with a binid less than zero is from a
            # spaxel that was not used
            indx = measurements_binid < 0
            map_mask[indx] = self.bitmask.turn_on(map_mask[indx], 'DIDNOTUSE')

            # The number of valid bins MUST match the number of
            # measurements
            nvalid = numpy.sum(numpy.invert(indx))
            if nvalid != len(measurements):
                raise ValueError('Provided id does not match the number of measurements.')

            # Get the bin ids with fitted models
            bin_indx = measurements_binid

        # Save the data to the hdu attribute
        self.hdu = fits.HDUList([fits.PrimaryHDU(header=pri_hdr),
                                 fits.ImageHDU(data=bin_indx, header=map_hdr, name='BINID'),
                                 fits.ImageHDU(data=map_mask, header=map_hdr, name='MAPMASK'),
                                 self.database.momdb.to_datatable().to_hdu(name='ELMBAND'),
                                 measurements.to_hdu(name='ELMMNTS')])

        #---------------------------------------------------------------
        # Write the data, if requested
        if self.hardcopy:
            if not self.directory_path.exists():
                self.diretory_path.mkdir(parents=True)
            self.write(overwrite=_overwrite)
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, '-'*50)

    def write(self, overwrite=False):
        """
        Write the hdu object to the file.
        """
        DAPFitsUtil.write(self.hdu, str(self.file_path()), overwrite=overwrite, checksum=True,
                          loggers=self.loggers, quiet=self.quiet)

    def read(self, ifile=None, strict=True, checksum=False):
        """

        Read an existing file with a previously binned set of spectra.

        """
        _ifile = self.file_path() if ifile is None else Path(ifile).resolve()
        if not _ifile.exists():
            raise FileNotFoundError('File does not exist!: {0}'.format(ifile))

        if self.hdu is not None:
            self.hdu.close()

        self.hdu = DAPFitsUtil.read(str(_ifile), checksum=checksum)

        # Confirm that the internal method is the same as the method
        # that was used in writing the file
        if self.hdu['PRIMARY'].header['ELMKEY'] != self.database['key']:
            if strict:
                raise ValueError('Keywords in header do not match specified method keywords!')
            else:
                warnings.warn('Keywords in header do not match specified method keywords!')
        # TODO: "strict" should also check other aspects of the file to
        # make sure that the details of the method are also the same,
        # not just the keyword

        self.nbins = self.hdu['PRIMARY'].header['NBINS']
        unique_bins = numpy.unique(self.hdu['BINID'].data.ravel()) \
                            if self.hdu['PRIMARY'].header['ELMREBIN'] else None
        self.missing_bins = self._get_missing_bins(unique_bins=unique_bins)

    # TODO: Use the EmissionMomentsDB.channel_names function!
    def channel_names(self):
        return numpy.array([ '{0}-{1}'.format(n,int(w)) \
                        for n,w in zip(self['ELMBAND'].data['NAME'],
                                       self['ELMBAND'].data['RESTWAVE']) ])

