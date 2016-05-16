# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
A class hierarchy that measures moments of the observed emission lines.

*License*:
    Copyright (c) 2015, SDSS-IV/MaNGA Pipeline Group
        Licensed under BSD 3-clause license - see LICENSE.rst

*Source location*:
    $MANGADAP_DIR/python/mangadap/proc/emissionlinemoments.py

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
            try:
                from configparser import ConfigParser
            except ImportError:
                warnings.warn('Unable to import configparser!  Beware!')
        else:
            try:
                from ConfigParser import ConfigParser
            except ImportError:
                warnings.warn('Unable to import ConfigParser!  Beware!')

        import glob
        import os.path
        import numpy
        from astropy.io import fits
        import astropy.constants

        from ..par.parset import ParSet
        from ..config.defaults import default_dap_source, default_dap_file_name
        from ..config.defaults import default_dap_method, default_dap_method_path
        from ..util.fileio import init_record_array, rec_to_fits_type, write_hdu
        from ..util.bitmask import BitMask
        from .artifactdb import ArtifactDB
        from .emissionmomentsdb import EmissionMomentsDB
        from .pixelmask import SpectralPixelMask
        from .spatiallybinnedspectra import SpatiallyBinnedSpectra
        from .stellarcontinuummodel import StellarContinuumModel
        from .bandpassfilter import passband_integral, passband_integrated_width
        from .bandpassfilter import passband_integrated_mean, passband_weighted_mean
        from .util import _select_proc_method

*Class usage examples*:
    Add examples!

*Revision history*:
    | **25 Apr 2016**: Implementation begun by K. Westfall (KBW)

.. _astropy.io.fits.hdu.hdulist.HDUList: http://docs.astropy.org/en/v1.0.2/io/fits/api/hdulists.html
.. _glob.glob: https://docs.python.org/3.4/library/glob.html
.. _configparser.ConfigParser: https://docs.python.org/3/library/configparser.html#configparser.ConfigParser


"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import sys
import warnings
if sys.version > '3':
    long = int
    try:
        from configparser import ConfigParser
    except ImportError:
        warnings.warn('Unable to import configparser!  Beware!')
else:
    try:
        from ConfigParser import ConfigParser
    except ImportError:
        warnings.warn('Unable to import ConfigParser!  Beware!')

import glob
import os.path
import numpy
from astropy.io import fits
import astropy.constants

from ..par.parset import ParSet
from ..config.defaults import default_dap_source, default_dap_file_name
from ..config.defaults import default_dap_method, default_dap_method_path
from ..util.fileio import init_record_array, rec_to_fits_type, write_hdu
from ..util.bitmask import BitMask
from .artifactdb import ArtifactDB
from .emissionmomentsdb import EmissionMomentsDB
from .pixelmask import SpectralPixelMask
from .spatiallybinnedspectra import SpatiallyBinnedSpectra
from .stellarcontinuummodel import StellarContinuumModel
from .bandpassfilter import passband_integral, passband_integrated_width
from .bandpassfilter import passband_integrated_mean, passband_weighted_mean
from .util import _select_proc_method

from matplotlib import pyplot

__author__ = 'Kyle B. Westfall'
# Add strict versioning
# from distutils.version import StrictVersion


class EmissionLineMomentsDef(ParSet):
    """
    A class that holds the parameters necessary to perform the
    emission-line moment measurements.

    Args:
        key (str): Keyword used to distinguish between different
            emission-line moment databases.
        minimum_snr (bool): Minimum S/N of spectrum to fit
        artifacts (str): String identifying the artifact database to use
        passbands (str): String identifying the emission-line bandpass
            filter database to use
    """
    def __init__(self, key, minimum_snr, artifacts, passbands):
        in_fl = [ int, float ]

        pars =     [ 'key', 'minimum_snr', 'artifacts', 'passbands' ]
        values =   [   key,   minimum_snr,   artifacts,   passbands ]
        dtypes =   [   str,         in_fl,         str,         str ]

        ParSet.__init__(self, pars, values=values, dtypes=dtypes)


def validate_emission_line_moments_config(cnfg):
    """ 
    Validate the `configparser.ConfigParser`_ object that is meant to
    define a set of emission-line moment measurements.

    Args:
        cnfg (`configparser.ConfigParser`_): Object meant to contain
            defining parameters of the emission-line moments as needed
            by :class:`EmissionLineMomentsDef'

    Raises:
        KeyError: Raised if any required keywords do not exist.
        ValueError: Raised if keys have unacceptable values.
        FileNotFoundError: Raised if a file is specified but could not
            be found.
    """
    # Check for required keywords
    if 'key' not in cnfg.options('default'):
        raise KeyError('Keyword \'key\' must be provided.')
    if 'minimum_snr' not in cnfg.options('default') or cnfg['default']['minimum_snr'] is None:
        cnfg['default']['minimum_snr']= '0.0'

    if 'artifact_mask' not in cnfg.options('default') \
            or cnfg['default']['artifact_mask'] is None:
        cnfg['default']['artifact_mask'] = 'None'

    if 'emission_passbands' not in cnfg.options('default') \
            or cnfg['default']['emission_passbands'] is None:
        raise ValueError('Must provide a keyword with the emission-line bandpass filter to use!')


def available_emission_line_moment_databases(dapsrc=None):
    """
    Return the list of available emission-line moment databases.

    Available database combinations:

    .. todo::
        Fill in

    Args:
        dapsrc (str): (**Optional**) Root path to the DAP source
            directory.  If not provided, the default is defined by
            :func:`mangadap.config.defaults.default_dap_source`.

    Returns:
        list: A list of :func:`EmissionLineMomentsDef` objects, each
        defining an emission-line moment database to measure.

    Raises:
        NotADirectoryError: Raised if the provided or default
            *dapsrc* is not a directory.
        OSError/IOError: Raised if no emission-line moment configuration
            files could be found.
        KeyError: Raised if the emission-line moment database keywords
            are not all unique.
        NameError: Raised if ConfigParser is not correctly imported.

    .. todo::
        - Somehow add a python call that reads the databases and
          constructs the table for presentation in sphinx so that the
          text above doesn't have to be edited with changes in the
          available databases.
        
    """
    # Check the source directory exists
    dapsrc = default_dap_source() if dapsrc is None else str(dapsrc)
    if not os.path.isdir(dapsrc):
        raise NotADirectoryError('{0} does not exist!'.format(dapsrc))

    # Check the configuration files exist
    search_dir = os.path.join(dapsrc, 'python/mangadap/config/emission_line_moments')
    ini_files = glob.glob(os.path.join(search_dir, '*.ini'))
    if len(ini_files) == 0:
        raise IOError('Could not find any configuration files in {0} !'.format(search_dir))

    # Build the list of library definitions
    moment_set_list = []
    for f in ini_files:
        # Read the config file
        cnfg = ConfigParser(allow_no_value=True)
        cnfg.read(f)
        # Ensure it has the necessary elements to define the template
        # library
        validate_emission_line_moments_config(cnfg)

        moment_set_list += [ EmissionLineMomentsDef(cnfg['default']['key'],
                                                    eval(cnfg['default']['minimum_snr']), 
                                                None if cnfg['default']['artifact_mask'] == 'None' \
                                                    else cnfg['default']['artifact_mask'],
                                                    cnfg['default']['emission_passbands']) ]

    # Check the keywords of the libraries are all unique
    if len(numpy.unique(numpy.array([moment['key'] for moment in moment_set_list]))) \
                != len(moment_set_list):
        raise KeyError('Emission-line moment database keywords are not all unique!')

    # Return the default list of assessment methods
    return moment_set_list


class EmissionLineMomentsBitMask(BitMask):
    r"""
    Derived class that specifies the mask bits for the emission-line
    moment measurements.  See :class:`mangadap.util.bitmask.BitMask` for
    attributes.

    A list of the bits and meanings are provided by the base class
    function :func:`mangadap.util.bitmask.BitMask.info`; i.e.,::

        from mangadap.proc.emissionlinemoments import EmissionLineMomentsBitMask
        bm = EmissionLineMomentsBitMask()
        bm.info()

    """
    def __init__(self, dapsrc=None):
        dapsrc = default_dap_source() if dapsrc is None else str(dapsrc)
        BitMask.__init__(self, ini_file=os.path.join(dapsrc, 'python', 'mangadap', 'config',
                                                     'bitmasks', 'emission_line_moments_bits.ini'))


class EmissionLineMoments:
    r"""

    Class that holds the emission-line moment measurements.

    """
    def __init__(self, database_key, binned_spectra, redshift=None, stellar_continuum=None,
                 database_list=None, artifact_list=None, bandpass_list=None, dapsrc=None,
                 dapver=None, analysis_path=None, directory_path=None, output_file=None,
                 hardcopy=True, clobber=False, verbose=0, checksum=False):

        self.version = '1.0'
        self.verbose = verbose

        # Define the database properties
        self.database = None
        self.artdb = None
        self.momdb = None
        self.pixelmask = None
        self._define_databases(database_key, database_list=database_list,
                               artifact_list=artifact_list, bandpass_list=bandpass_list)

        self.binned_spectra = None
        self.redshift = None
        self.stellar_continuum = None
        self.nmom = self.momdb.nsets

        # Define the output directory and file
        self.directory_path = None      # Set in _set_paths
        self.output_file = None
        self.hardcopy = None

        # Initialize the objects used in the assessments
        self.bitmask = EmissionLineMomentsBitMask(dapsrc=dapsrc)
        self.hdu = None
        self.checksum = checksum

        self.spatial_shape = None
        self.spatial_index = None

        self.nbins = None
        self.missing_bins = None

        # Run the assessments of the DRP file
        self.measure(binned_spectra, redshift=redshift, stellar_continuum=stellar_continuum,
                     dapsrc=dapsrc, dapver=dapver, analysis_path=analysis_path,
                     directory_path=directory_path, output_file=output_file, hardcopy=hardcopy,
                     clobber=clobber, verbose=verbose)


    def __del__(self):
        """
        Deconstruct the data object by ensuring that the fits file is
        properly closed.
        """
        if self.hdu is None:
            return
        self.hdu.close()
        self.hdu = None


    def __getitem__(self, key):
        return self.hdu[key]


    def _define_databases(self, database_key, database_list=None, artifact_list=None,
                          bandpass_list=None, dapsrc=None):
        r"""

        Select the database of bandpass filters.

        """
        # Grab the specific database
        self.database = _select_proc_method(database_key, EmissionLineMomentsDef,
                                            method_list=database_list,
                                            available_func=available_emission_line_moment_databases,
                                            dapsrc=dapsrc)

        # Instantiate the artifact and bandpass filter database
        self.artdb = None if self.database['artifacts'] is None else \
                ArtifactDB(self.database['artifacts'], artdb_list=artifact_list, dapsrc=dapsrc)
        self.pixelmask = SpectralPixelMask(artdb=self.artdb)

        self.momdb = None if self.database['passbands'] is None else \
                EmissionMomentsDB(self.database['passbands'], emldb_list=bandpass_list,
                                  dapsrc=dapsrc)


    def _set_paths(self, directory_path, dapver, analysis_path, output_file):
        """
        Set the :attr:`directory_path` and :attr:`output_file`.  If not
        provided, the defaults are set using, respectively,
        :func:`mangadap.config.defaults.default_dap_reference_path` and
        :func:`mangadap.config.defaults.default_dap_file_name`.

        Args:
            directory_path (str): The exact path to the DAP
                emission-line moments file.  See :attr:`directory_path`.
            dapver (str): DAP version.
            analysis_path (str): The path to the top-level directory
                containing the DAP output files for a given DRP and DAP
                version.
            output_file (str): The name of the file with emission-line
                moment measurements.  See :func:`measure`.
        """
        # Set the output directory path
        method = default_dap_method(binned_spectra=self.binned_spectra,
                                    stellar_continuum=self.stellar_continuum)
        self.directory_path = default_dap_method_path(method, plate=self.binned_spectra.drpf.plate,
                                                      ifudesign=self.binned_spectra.drpf.ifudesign,
                                                      ref=True,
                                                      drpver=self.binned_spectra.drpf.drpver,
                                                      dapver=dapver, analysis_path=analysis_path) \
                                        if directory_path is None else str(directory_path)

        # Set the output file
        ref_method = '{0}-{1}'.format(self.binned_spectra.rdxqa.method['key'],
                                      self.binned_spectra.method['key'])
        if self.stellar_continuum is not None:
            ref_method = '{0}-{1}'.format(method, self.stellar_continuum.method['key'])
        ref_method = '{0}-{1}'.format(method, self.database['key'])
        self.output_file = default_dap_file_name(self.binned_spectra.drpf.plate,
                                                 self.binned_spectra.drpf.ifudesign, ref_method) \
                                        if output_file is None else str(output_file)


    def _clean_drp_header(self, ext='FLUX'):
        """
        Read and clean the header from the DRP fits file extension for
        use in the output file for this class object.

        .. todo::
            - Currently just returns the existing header.  Need to
              decide if/how to manipulate DRP header.

        """
        return self.binned_spectra.drpf[ext].header.copy()


    def _initialize_header(self, hdr):
        """

        Initialize the header.

        """
        hdr['ELMKEY'] = (self.database['key'], 'Emission-line moments database keyword')
        hdr['ELMMINSN'] = (self.database['minimum_snr'], 'Minimum S/N of spectrum to include')
        hdr['ARTDB'] = (self.database['artifacts'], 'Artifact database keyword')
        hdr['MOMDB'] = (self.database['passbands'], 'Emission-line moments database keyword')
        if self.stellar_continuum is not None:
            hdr['SCTPL'] = (self.stellar_continuum.method['key'], 'Stellar-continuum model keyword')
        hdr['NBINS'] = (self.nbins, 'Number of unique spatial bins')
        if len(self.missing_bins) > 0:
            hdr['EMPTYBIN'] = (str(self.missing_bins), 'List of bins with no data')
        # Anything else?
        # Additional database details?
        return hdr


    def _moments_database_dtype(self, name_len):
        r"""
        Construct the record array data type for the output fits
        extension.
        """
        return [ ('ID',numpy.int),
                 ('NAME','<U{0:d}'.format(name_len)),
                 ('RESTWAVE', numpy.float),
                 ('PASSBAND', numpy.float, 2),
                 ('BLUEBAND', numpy.float, 2),
                 ('REDBAND', numpy.float, 2)
               ]

    
    def _compile_database(self):
        """
        Compile the database with the specifications of each index.
        """
        name_len = 0
        for n in self.momdb['name']:
            if name_len < len(n):
                name_len = len(n)

        # Instatiate the table data that will be saved defining the set
        # of emission-line moments measured
        hdu_database = init_record_array(self.nmom, self._moments_database_dtype(name_len))

        hk = [ 'ID', 'NAME', 'RESTWAVE', 'PASSBAND', 'BLUEBAND', 'REDBAND' ]
        mk = [ 'index', 'name', 'restwave', 'primary', 'blueside', 'redside' ]
        for _hk, _mk in zip(hk,mk):
            hdu_database[_hk] = self.momdb[_mk]

        #print(hdu_database)

        return hdu_database


    def _per_bin_dtype(self):
        r"""
        Construct the record array data type for the output fits
        extension.
        """
        return [ ('BIN_INDEX',numpy.int),
                 ('REDSHIFT', numpy.float),
                 ('MASK', self.bitmask.minimum_uint_dtype(), self.nmom),
                 ('BCEN', numpy.float, self.nmom), 
                 ('BCONT', numpy.float, self.nmom), 
                 ('BCONTERR', numpy.float, self.nmom),
                 ('RCEN', numpy.float, self.nmom), 
                 ('RCONT', numpy.float, self.nmom), 
                 ('RCONTERR', numpy.float, self.nmom), 
                 ('CNTSLOPE', numpy.float, self.nmom),
                 ('FLUX', numpy.float, self.nmom), 
                 ('FLUXERR', numpy.float, self.nmom),
                 ('MOM1', numpy.float, self.nmom), 
                 ('MOM1ERR', numpy.float, self.nmom),
                 ('MOM2', numpy.float, self.nmom), 
                 ('MOM2ERR', numpy.float, self.nmom)
               ]


    def _check_snr(self):
        # binned_spectra['BINS'].data['SNR'] has length
        # binned_spectra.nbins
        return self.binned_spectra['BINS'].data['SNR'] > self.database['minimum_snr']

    
    def _missing_flags(self):
        return numpy.array([ b in self.missing_bins for b in numpy.arange(self.nbins)])


    def _bins_to_measure(self):
        """
        Determine which bins to use for the emission-line moment
        measurements.  They must not be designated as "missing" and they
        must have sufficient S/N.
        """
        return (self._check_snr()) & ~(self._missing_flags())


    def _assign_redshifts(self, redshift):
        """
        Set the redshift to apply to each spectrum.  Directly provided
        redshifts take precedence over those in the
        StellarContinuumModel object, if both are provided.

        The spectrum must not be a "missing bin" and have good S/N.
        """
        
        self.redshift = numpy.zeros(self.nbins, dtype=numpy.float)

        if redshift is not None:
            _redshift = numpy.asarray(redshift)
            if len(_redshift) not in [ 1, self.nbins, self.nspec ]:
                raise ValueError('Provided redshift must be either a single value, match the ' \
                                 'number of binned spectra, or match the total number of DRP ' \
                                 'spectra.')
            if len(_redshift) == 1:
                self.redshift[:] = redshift
            elif len(_redshift) == self.nspec:
                self.redshift = _redshift[ self.unique_bins(index=True) ]
            else:   # Has length nbins
                self.redshift = _redshift
        elif self.stellar_continuum is not None:
            self.redshift = self.stellar_continuum['PAR'].data['KIN'][:,0] \
                                / astropy.constants.c.to('km/s').value
        self.redshift = self.redshift[self._bins_to_measure()]


    def _spectra_for_measurements(self):
        """
        Compile the set of spectra for the emission-line moment
        measurements.  If a stellar continuum model is provided, this
        function will also provide the model-subtracted flux and a
        boolean array flagged with regions where the model was not
        subtracted from the data.
        """
        # Get the data arrays
        wave = self.binned_spectra['WAVE'].data
        flux = self.binned_spectra.copy_to_masked_array(flag=self.binned_spectra.do_not_fit_flags())
        ivar = self.binned_spectra.copy_to_masked_array(ext='IVAR',
                                                        flag=self.binned_spectra.do_not_fit_flags())
        indx = self.pixelmask.boolean(wave,nspec=self.nbins)
        flux[indx] = numpy.ma.masked
        ivar[indx] = numpy.ma.masked

        mask = numpy.zeros((self.nbins,self.binned_spectra.nwave),
                           dtype=self.bitmask.minimum_uint_dtype())
        flux_mask = numpy.ma.getmaskarray(flux)
        mask[flux_mask] = self.bitmask.turn_on(mask[flux_mask], 'DIDNOTUSE')

        good_bins = self._bins_to_measure()
        flux = flux[good_bins,:]
        ivar = ivar[good_bins,:]
        mask = mask[good_bins,:]

        if self.stellar_continuum is None:
            return flux, ivar, mask, None, None

        # Get the models for the binned spectra
        model = self.stellar_continuum.copy_to_masked_array(
                        flag=self.stellar_continuum.all_except_emission_flags())[good_bins,:]
        # Get where the stellar-continuum models are masked but the
        # binned spectra are not
        no_model = numpy.invert(numpy.ma.getmaskarray(flux)) & numpy.ma.getmaskarray(model)

        model_subtracted_flux = flux - model
        model_subtracted_flux.mask[no_model] = False

#        pyplot.step(wave, flux[0,:], where='mid', linestyle='-', color='k', lw=0.5, zorder=2)
#        pyplot.plot(wave, model[0,:], linestyle='-', color='r', lw=1.5, zorder=1, alpha=0.5)
#        pyplot.step(wave, model_subtracted_flux[0,:], where='mid', linestyle='-', color='g',
#                    lw=0.5, zorder=3)
#        pyplot.show()

        return flux, ivar, mask, model_subtracted_flux, no_model


    def _sideband_pseudocontinua(self, spec, sidebands, hdu_measurements, spec_n=None, noise=None):
        """Get the side-band integrals."""

        # Calculate the pseudo-continua in the sidebands
        wave = self.binned_spectra['WAVE'].data

        # If more than one spectrum entered, use spec_n to decide which
        # spectrum to use for the integrals
        if len(spec.shape) > 1:
            # In this case, spec_n should always be provided, but allow
            # the code to continue with a warning
            if spec_n is None:
                warnings.warn('Multiple spectra entered, but no designation for which to use!')
                _spec_n = numpy.zeros(self.nmom*2, dtype=numpy.int)
            else:
                _spec_n = numpy.append(spec_n, spec_n)

            # Instantiate the arrays
            pseudocontinuum = numpy.zeros(self.nmom*2)
            pseudocontinuum_error = None if noise is None else numpy.zeros(self.nmom*2)
            flux_weighted_center = numpy.zeros(self.nmom*2)
            interval_frac = numpy.zeros(self.nmom*2)

            # Use the binned spectrum
            pseudocontinuum[_spec_n==0], _pseudocontinuum_error \
                    = passband_integrated_mean(wave, spec[0,:], passband=sidebands[_spec_n==0,:],
                                               err=noise, log=True)
            if noise is not None:
                pseudocontinuum_error[_spec_n==0] = _pseudocontinuum_error
            flux_weighted_center[_spec_n==0], _fwc_err \
                    = passband_weighted_mean(wave, spec[0,:]+1.0, wave,
                                             passband=sidebands[_spec_n==0,:], log=True)
            interval_frac[_spec_n==0] \
                    = passband_integrated_width(wave, spec[0,:], passband=sidebands[_spec_n==0,:],
                                                log=True) \
                                / numpy.diff(sidebands[_spec_n==0,:], axis=1).ravel()

            # Use the model-subtracted binned spectrum
            pseudocontinuum[_spec_n==1], _pseudocontinuum_error \
                    = passband_integrated_mean(wave, spec[1,:], passband=sidebands[_spec_n==1],
                                               err=noise, log=True)
            if noise is not None:
                pseudocontinuum_error[_spec_n==1] = _pseudocontinuum_error
            flux_weighted_center[_spec_n==1], _fwc_err \
                    = passband_weighted_mean(wave, spec[1,:]+1.0, wave,
                                             passband=sidebands[_spec_n==1,:], log=True)
            interval_frac[_spec_n==1] \
                    = passband_integrated_width(wave, spec[1,:], passband=sidebands[_spec_n==1,:],
                                                log=True) \
                                / numpy.diff(sidebands[_spec_n==1,:], axis=1).ravel()

        else:
            # Use same spectrum for all sidebands
            pseudocontinuum, pseudocontinuum_error \
                    = passband_integrated_mean(wave, spec, passband=sidebands, err=noise, log=True)
            flux_weighted_center, _fwc_err \
                    = passband_weighted_mean(wave, spec+1.0, wave, passband=sidebands, log=True)
            # Calculate the fraction of the band that is covered by
            # unmasked pixels
            interval_frac = passband_integrated_width(wave, spec, passband=sidebands, log=True) \
                                / numpy.diff(sidebands, axis=1).ravel()

        # Save the results to the main output array
        hdu_measurements['BCEN'] = flux_weighted_center.reshape(-1,self.nmom)[0,:] - 1.0
        hdu_measurements['RCEN'] = flux_weighted_center.reshape(-1,self.nmom)[1,:] - 1.0
        hdu_measurements['BCONT'] = pseudocontinuum.reshape(-1,self.nmom)[0,:]
        hdu_measurements['RCONT'] = pseudocontinuum.reshape(-1,self.nmom)[1,:]
        if pseudocontinuum_error is not None:
            hdu_measurements['BCONTERR'] = pseudocontinuum_error.reshape(-1,self.nmom)[0,:]
            hdu_measurements['RCONTERR'] = pseudocontinuum_error.reshape(-1,self.nmom)[1,:]

#        pyplot.step(wave, spec[0], where='mid', color='k', lw=0.5, zorder=1)
#        pyplot.step(wave, spec[1], where='mid', color='g', lw=0.5, zorder=2)
#        pyplot.scatter(hdu_measurements['BCEN'], hdu_measurements['BCONT'], marker='.', s=50,
#                       color='b', lw=0, zorder=3)
#        pyplot.scatter(hdu_measurements['RCEN'], hdu_measurements['RCONT'], marker='.', s=50,
#                       color='r', lw=0, zorder=3)
#        pyplot.show()

        # Passband is not fully covered by valid pixels
        blue_fraction = interval_frac.reshape(-1,self.nmom)[0,:]
#        print(blue_fraction)
        hdu_measurements['MASK'][blue_fraction < 1.0] \
                = self.bitmask.turn_on(hdu_measurements['MASK'][blue_fraction < 1.0], 'BLUE_INCOMP')
        red_fraction = interval_frac.reshape(-1,self.nmom)[1,:]
#        print(red_fraction)
        hdu_measurements['MASK'][red_fraction < 1.0] \
                = self.bitmask.turn_on(hdu_measurements['MASK'][red_fraction < 1.0], 'RED_INCOMP')
        # Passband is empty
        hdu_measurements['MASK'][~(blue_fraction > 0.0)] \
            = self.bitmask.turn_on(hdu_measurements['MASK'][~(blue_fraction > 0.0)], 'BLUE_EMPTY')
        hdu_measurements['MASK'][~(red_fraction > 0.0)] \
            = self.bitmask.turn_on(hdu_measurements['MASK'][~(red_fraction > 0.0)], 'RED_EMPTY')


    def _single_band_moments(self, wave, spec, passband, restwave, mask, noise=None):
        """
        Measure the moments for a single band.
        """
        # Get the fraction of the passband that is unmasked
        interval_frac = passband_integrated_width(wave, spec, passband=passband, log=True) \
                                / numpy.diff(passband).ravel()
        if interval_frac < 1.0:
            mask = self.bitmask.turn_on(mask, 'MAIN_INCOMP')
        if not interval_frac > 0.0:
            mask = self.bitmask.turn_on(mask, 'MAIN_EMPTY')
            # Nothing in the passband to integrate!
            return None, None, None, None, None, None, mask

        # Get the redshift vector
        cz = astropy.constants.c.to('km/s').value*(wave/restwave-1)

        # Get the integrated flux
        flux = passband_integral(wave, spec, passband=passband, log=True)
        fluxerr = None if noise is None else \
                    numpy.ma.sqrt( passband_integral(wave, numpy.square(noise), passband=passband,
                                                     log=True))

        # No flux in the passband
        if not numpy.absolute(flux) > 0.0:
            mask = self.bitmask.turn_on(mask, 'DIVBYZERO')
            return flux, fluxerr, None, None, None, None, mask

        # Get the first and second moments; saved to a temporary
        # variable for use in the error calculations (if possible)
        mom1 = passband_integral(wave, cz*spec, passband=passband, log=True)
        _mom1 = mom1 / flux
        mom2 = passband_integral(wave, numpy.square(cz)*spec, passband=passband, log=True)
        _mom2 = mom2 / flux

        # Perform the error calculations
        if noise is not None:
            mom1err = passband_integral(wave, cz*numpy.square(noise), passband=passband, log=True)
            mom1err = numpy.ma.sqrt(mom1err * numpy.square(_mom1/mom1)
                                        + numpy.square(fluxerr * _mom1/flux))
            mom1 = _mom1
                
            if _mom2 < numpy.square(mom1):
                mask = self.bitmask.turn_on(mask, 'UNDEFINED_MOM2')
                return flux, fluxerr, mom1, mom1err, None, None, mask

            mom2err = passband_integral(wave, numpy.square(cz*noise), passband=passband, log=True)
            mom2err = numpy.ma.sqrt(mom2err * numpy.square(_mom2/mom2)
                                        + numpy.square(fluxerr * _mom2/flux))
            mom2err = numpy.ma.sqrt(numpy.square(mom2err) + numpy.square(2*mom1*mom1err))
            mom2 = numpy.ma.sqrt(_mom2 - numpy.square(mom1))
            mom2err /= (2.0*mom2)
        else:
            mom1 = _mom1
            mom1err = None

            if _mom2 < numpy.square(mom1):
                mask = self.bitmask.turn_on(mask, 'UNDEFINED_MOM2')
                return flux, fluxerr, mom1, mom1err, None, None, mask

            mom2 = numpy.ma.sqrt(_mom2 - numpy.square(mom1))
            mom2err = None

        return flux, fluxerr, mom1, mom1err, mom2, mom2err, mask


    def _continuum_subtracted_moments(self, spec, mainbands, hdu_measurements, spec_n=None,
                                      noise=None):
        """
        Calculate the continuum-subtracted moments
        """

        # Get the parameters for the linear continuum across the
        # primary passband

        indx = numpy.invert(self.bitmask.flagged(hdu_measurements['MASK'],
                                    flag=[ 'BLUE_EMPTY', 'RED_EMPTY', 'BLUE_JUMP', 'RED_JUMP' ]))

        hdu_measurements['CNTSLOPE'][indx] \
                = (hdu_measurements['RCONT'][indx] - hdu_measurements['BCONT'][indx]) \
                        / (hdu_measurements['RCEN'][indx] - hdu_measurements['BCEN'][indx])
        continuum_b = hdu_measurements['BCONT'] \
                                - hdu_measurements['BCEN']*hdu_measurements['CNTSLOPE']

        wave = self.binned_spectra['WAVE'].data

        moment_keys = numpy.array(['FLUX', 'FLUXERR', 'MOM1', 'MOM1ERR', 'MOM2', 'MOM2ERR', 'MASK'])

        # If more than one spectrum entered, use spec_n to decide which
        # spectrum to use for the integrals
        if len(spec.shape) > 1:
            # In this case, spec_n should always be provided, but allow
            # the code to continue with a warning
            if spec_n is None:
                warnings.warn('Multiple spectra entered, but no designation for which to use!')
                _spec_n = numpy.zeros(self.nmom, dtype=numpy.int)
            else:
                _spec_n = spec_n

            # Get the moments for all the main passbands
            for i,p in enumerate(mainbands):
                cont = continuum_b[i] + hdu_measurements['CNTSLOPE'][i]*wave
                _spec = spec[_spec_n[i],:] - cont
                moments = list(self._single_band_moments(wave, _spec, p, self.momdb['restwave'][i],
                                                         hdu_measurements['MASK'][i], noise=noise))
                for j,k in enumerate(moment_keys):
                    if moments[j] is None:
                        continue
                    hdu_measurements[k][i] = moments[j]
        else:
            # Get the moments for all the main passbands
            for i,p in enumerate(mainbands):
                cont = continuum_b[i] + hdu_measurements['CNTSLOPE'][i]*wave
                _spec = spec - cont
                moments = list(self._single_band_moments(wave, _spec, p, self.momdb['restwave'][i],
                                                         hdu_measurements['MASK'][i], noise=noise))
                for j,k in enumerate(moment_keys):
                    if moments[j] is None:
                        continue
                    hdu_measurements[k][i] = moments[j]


    def _measure_moments(self, redshift, flux, ivar=None, mask=None, model_subtracted_flux=None,
                         no_model=None):
        """
        Measure the emission-line moments.

        If not input as masked arrays, flux is converted to one.
        """

        # Instatiate the table data that will be saved with the index
        # measurements
        nspec = flux.shape[0]
        wave = self.binned_spectra['WAVE'].data
        hdu_measurements = init_record_array(nspec, self._per_bin_dtype())
        hdu_measurements['REDSHIFT'] = redshift
        
        # Convert the input arrays to masked arrays
        _flux = flux if isinstance(flux, numpy.ma.MaskedArray) else \
                    numpy.ma.MaskedArray(flux, mask=(None if mask is None else mask > 0))

        if ivar is None:
            noise = None
        else:
            _ivar = ivar if isinstance(ivar, numpy.ma.MaskedArray) else \
                        numpy.ma.MaskedArray(ivar, mask=(None if mask is None else mask > 0))
            noise = numpy.ma.sqrt(1.0 /_ivar)

        if model_subtracted_flux is None:
            _model_subtracted_flux = None
        else:
            _model_subtracted_flux = model_subtracted_flux \
                        if isinstance(model_subtracted_flux, numpy.ma.MaskedArray) else \
                            numpy.ma.MaskedArray(model_subtracted_flux,
                                                 mask=(None if mask is None else mask > 0))

        # Use the no_model flags to determine if a band jumps between
        # regions that have and have not had the stellar-continuum model
        # subtracted from it
        _no_model = None if no_model is None else numpy.ma.MaskedArray(numpy.ones(no_model.shape),
                                                                       mask=no_model)
#        print('no_model exists:', _no_model is not None)
            
        # Measure the pseudo-continuum in the sidebands
        sidebands = numpy.append(self.momdb['blueside'], self.momdb['redside'], axis=0)
        bandpass = numpy.append(self.momdb['primary'], self.momdb['redside'], axis=0)
        for i in range(nspec):

            print('Measuring emission-line moments in spectrum: {0}/{1}'.format(i+1,nspec),
                  end='\r')

            # Shift the sidebands to the appropriate redshift
            _sidebands = sidebands*(1.0+redshift[i])
            _mainbands = self.momdb['primary']*(1.0+redshift[i])
            _noise = None if noise is None else noise[i,:]

            # Check if the model was subtracted over the full range of
            # the sidebands, or if the model was subtracted in one
            # sideband but not the other
            if _no_model is not None:
                interval_frac = passband_integrated_width(wave, _no_model[i,:], passband=_sidebands,
                                        log=True) / numpy.diff(_sidebands, axis=1).ravel()
                blue_fraction = interval_frac.reshape(-1,self.nmom)[0,:]
                indx = numpy.logical_and(blue_fraction > 0.0, blue_fraction < 1.0)
                hdu_measurements['MASK'][i,indx] \
                        = self.bitmask.turn_on(hdu_measurements['MASK'][i,indx], 'BLUE_JUMP')
                red_fraction = interval_frac.reshape(-1,self.nmom)[1,:]
                indx = numpy.logical_and(red_fraction > 0.0, red_fraction < 1.0)
                hdu_measurements['MASK'][i,indx] \
                        = self.bitmask.turn_on(hdu_measurements['MASK'][i,indx], 'RED_JUMP')
                
                indx = (~(blue_fraction > 0) & (red_fraction > 0)) \
                            | (~(red_fraction > 0) & (blue_fraction > 0))
                hdu_measurements['MASK'][i,indx] \
                        = self.bitmask.turn_on(hdu_measurements['MASK'][i,indx],
                                               'JUMP_BTWN_SIDEBANDS')

                # Only use model-subtracted measurements if there are no
                # jumps in the bands or between the bands
                jumps = self.bitmask.flagged(hdu_measurements['MASK'][i,:],
                                             flag=['BLUE_JUMP', 'RED_JUMP', 'JUMP_BTWN_SIDEBANDS'])

                spec = numpy.ma.append(_flux[i,:].reshape(1,-1),
                                       _model_subtracted_flux[i,:].reshape(1,-1), axis=0)
                spec_n = numpy.zeros(self.nmom, dtype=numpy.int)
                spec_n[~(jumps)] = 1
                self._sideband_pseudocontinua(spec, _sidebands, hdu_measurements[i],
                                              spec_n=spec_n, noise=_noise)

                self._continuum_subtracted_moments(spec, _mainbands, hdu_measurements[i],
                                                   spec_n=spec_n, noise=_noise)

                # Flag moment calculations as having not corrected for
                # stellar absorption
                hdu_measurements['MASK'][i,spec_n == 0] \
                            = self.bitmask.turn_on(hdu_measurements['MASK'][i,spec_n == 0],
                                                   'NO_ABSORPTION_CORRECTION')

            else:
                spec = _flux[i,:] if _model_subtracted_flux is None else _model_subtracted_flux[i,:]
                self._sideband_pseudocontinua(spec, _sidebands, hdu_measurements[i], noise=_noise)
                self._continuum_subtracted_moments(spec, _mainbands, hdu_measurements[i],
                                                   noise=_noise)

                # Flag moment calculations as having not corrected for
                # stellar aborption
                if _model_subtracted_flux is None:
                    hdu_measurements['MASK'][i,:] \
                            = self.bitmask.turn_on(hdu_measurements['MASK'][i,:],
                                                   'NO_ABSORPTION_CORRECTION')

        print('Measuring emission-line moments in spectrum: {0}/{0}'.format(nspec))
        return hdu_measurements
        

    def file_name(self):
        """Return the name of the output file."""
        return self.output_file


    def file_path(self):
        """Return the full path to the output file."""
        if self.directory_path is None or self.output_file is None:
            return None
        return os.path.join(self.directory_path, self.output_file)


    def measure(self, binned_spectra, redshift=None, stellar_continuum=None, dapsrc=None,
                dapver=None, analysis_path=None, directory_path=None, output_file=None,
                hardcopy=True, clobber=False, verbose=0):
        """

        Measure the emission-line moments using the binned spectra.

        """
        # SpatiallyBinnedSpectra object always needed
        if binned_spectra is None:
            raise ValueError('Must provide spectra object for fitting.')
        if not isinstance(binned_spectra, SpatiallyBinnedSpectra):
            raise TypeError('Must provide a valid SpatiallyBinnedSpectra object!')
        if binned_spectra.hdu is None:
            raise ValueError('Provided SpatiallyBinnedSpectra object is undefined!')
        self.binned_spectra = binned_spectra

        # StellarContinuumModel object only used when accounting for
        # underlying absorption
        if stellar_continuum is not None:
            if not isinstance(stellar_continuum, StellarContinuumModel):
                raise TypeError('Provided stellar continuum must have StellarContinuumModel type!')
            if stellar_continuum.hdu is None:
                raise ValueError('Provided StellarContinuumModel is undefined!')
            self.stellar_continuum = stellar_continuum

        self.spatial_shape =self.binned_spectra.spatial_shape
        self.spatial_index = self.binned_spectra.spatial_index.copy()
        
        self.nbins = self.binned_spectra.nbins
        self.missing_bins = self.binned_spectra.missing_bins

#        pyplot.scatter(numpy.arange(self.nspec), self.redshift, marker='.', s=50, color='k', lw=0)
#        pyplot.show()
            
        # Get the good spectra
        #   - Must have sufficienct S/N, as defined by the input par
        good_snr = self._check_snr()
        # Report
        print('Total spectra: ', len(good_snr))
        print('With good S/N: ', numpy.sum(good_snr))

        # Get the redshifts to apply
        self._assign_redshifts(redshift)

        # (Re)Set the output paths
        self._set_paths(directory_path, dapver, analysis_path, output_file)

        # Check that the file path is defined
        ofile = self.file_path()
        if ofile is None:
            raise ValueError('File path for output file is undefined!')

        # Report
        print('Output path: ', self.directory_path)
        print('Output file: ', self.output_file)

        # If the file already exists, and not clobbering, just read the
        # file
        if os.path.isfile(ofile) and not clobber:
            self.read()
            return

        # Get the spectra to use for the measurements
        flux, ivar, mask, model_subtracted_flux, no_model = self._spectra_for_measurements()

        # Compile the information on the suite of measured indices
        hdu_database = self._compile_database()

        # Perform the measurements on the galaxy spectra
        hdu_measurements = self._measure_moments(self.redshift, flux, ivar=ivar, mask=mask,
                                                 model_subtracted_flux=model_subtracted_flux,
                                                 no_model=no_model)

        # Fill array for any missing bins in prep for writing to the HDU
        if len(hdu_measurements) < self.nbins:
            # Fill the full array for all bins
            _hdu_measurements = hdu_measurements.copy()
            hdu_measurements = init_record_array(self.nbins, self._per_bin_dtype())
            if len(self.missing_bins) > 0:
                hdu_measurements['MASK'][self.missing_bins] \
                    = self.bitmask.turn_on(hdu_measurements['MASK'][self.missing_bins], 'DIDNOTUSE')
            hdu_measurements[good_snr] = _hdu_measurements

        hdu_measurements['BIN_INDEX'] = numpy.arange(self.nbins)
        hdu_measurements['MASK'][~(good_snr)] \
                = self.bitmask.turn_on(hdu_measurements['MASK'][~(good_snr)], 'LOW_SNR') 

        # Initialize the header keywords
        hdr = self._clean_drp_header(ext='PRIMARY')
        self._initialize_header(hdr)

        # Save the data to the hdu attribute
        self.hdu = fits.HDUList([ fits.PrimaryHDU(header=hdr),
                                  fits.BinTableHDU.from_columns( [ fits.Column(name=n,
                                                        format=rec_to_fits_type(hdu_database[n]),
                                array=hdu_database[n]) for n in hdu_database.dtype.names ],
                                                               name='ELMBAND'),
                                  fits.BinTableHDU.from_columns( [ fits.Column(name=n,
                                                    format=rec_to_fits_type(hdu_measurements[n]),
                                array=hdu_measurements[n]) for n in hdu_measurements.dtype.names ],
                                                               name='ELMMNTS')
                                ])

        # Write the data, if requested
        if not os.path.isdir(self.directory_path):
            os.makedirs(self.directory_path)
        self.hardcopy = hardcopy
        if self.hardcopy:
            self.write(clobber=clobber)


    def write(self, clobber=False):
        """
        Write the hdu object to the file.
        """
        # Get the output file and determine if it should be compressed
        ofile = self.file_path()
        write_hdu(self.hdu, ofile, clobber=clobber, checksum=True)


    def read(self, ifile=None, strict=True, checksum=False):
        """

        Read an existing file with a previously binned set of spectra.

        """
        if ifile is None:
            ifile = self.file_path()
        if not os.path.isfile(ifile):
            raise FileNotFoundError('File does not exist!: {0}'.format(ifile))

        if self.hdu is not None:
            self.hdu.close()

        self.hdu = fits.open(ifile, checksum=checksum)

        # Confirm that the internal method is the same as the method
        # that was used in writing the file
        if self.hdu['PRIMARY'].header['ELMKEY'] != self.database['key']:
            if strict:
                raise ValueError('Keywords in header does not match specified method keyword!')
            else:
                warnings.warn('Keywords in header does not match specified method keyword!')
        # TODO: "strict" should also check other aspects of the file to
        # make sure that the details of the method are also the same,
        # not just the keyword

        self.nbins = self.hdu['PRIMARY'].header['NBINS']
        try:
            self.missing_bins = eval(self.hdu['PRIMARY'].header['EMPTYBIN'])
        except KeyError:
            # Assume if this fails, it's because the keyword doesn't
            # exist
            self.missing_bins = []


    def unique_bins(self, index=False):
        """
        Get the unique bins or the indices of the unique bins in the
        flattened spatial dimension.
        """
        unique_bins, first_occurrence = numpy.unique(self.binned_spectra['BINID'].data.ravel(),
                                                     return_index=True)
        return first_occurrence[1:] if index else unique_bins[1:]

