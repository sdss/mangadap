# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
A class hierarchy that performs the spectral-index measurements.

*License*:
    Copyright (c) 2015, SDSS-IV/MaNGA Pipeline Group
        Licensed under BSD 3-clause license - see LICENSE.rst

*Source location*:
    $MANGADAP_DIR/python/mangadap/proc/spectralindices.py

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
        from ..config.defaults import default_dap_source, default_dap_reference_path
        from ..config.defaults import default_dap_file_name
        from ..util.instrument import spectral_resolution, match_spectral_resolution
        from ..util.instrument import spectral_coordinate_step
        from ..util.fileio import init_record_array, rec_to_fits_type, write_hdu
        from ..util.bitmask import BitMask
        from .artifactdb import ArtifactDB
        from .absorptionindexdb import AbsorptionIndexDB
        from .bandheadindexdb import BandheadIndexDB
        from .pixelmask import StellarContinuumPixelMask
        from .spatiallybinnedspectra import SpatiallyBinnedSpectra
        from .templatelibrary import TemplateLibrary
        from .stellarcontinuummodel import StellarContinuumModel
        from .bandpassfilter import passband_integral, passband_integrated_width, passband_integrated_mean
        from .util import _select_proc_method, flux_to_fnu


*Class usage examples*:
    Add examples!

*Revision history*:
    | **20 Apr 2016**: Implementation begun by K. Westfall (KBW)

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

from ..drpfits import DRPFits
from ..par.parset import ParSet
from ..config.defaults import default_dap_source, default_dap_reference_path
from ..config.defaults import default_dap_file_name
from ..util.instrument import spectral_resolution, match_spectral_resolution
from ..util.instrument import spectral_coordinate_step
from ..util.fileio import init_record_array, rec_to_fits_type, write_hdu
from ..util.bitmask import BitMask
from .artifactdb import ArtifactDB
from .absorptionindexdb import AbsorptionIndexDB
from .bandheadindexdb import BandheadIndexDB
from .pixelmask import StellarContinuumPixelMask
from .spatiallybinnedspectra import SpatiallyBinnedSpectra
from .templatelibrary import TemplateLibrary
from .stellarcontinuummodel import StellarContinuumModel
from .bandpassfilter import passband_integral, passband_integrated_width, passband_integrated_mean
from .util import _select_proc_method, flux_to_fnu

from matplotlib import pyplot

__author__ = 'Kyle B. Westfall'
# Add strict versioning
# from distutils.version import StrictVersion


class SpectralIndicesDef(ParSet):
    """
    A class that holds the parameters necessary to perform the
    spectral-index measurements.

    Args:
        key (str): Keyword used to distinguish between different spatial
            binning schemes.
        minimum_snr (bool): Minimum S/N of spectrum to fit
        fwhm (bool): Resolution FWHM in angstroms at which to make the
            measurements.
        artifacts (str): String identifying the artifact database to use
        absindex (str): String identifying the absorption-index database
            to use
        bandhead (str): String identifying the bandhead-index database
            to use
    """
    def __init__(self, key, minimum_snr, fwhm, artifacts, absindex, bandhead):
        in_fl = [ int, float ]

        pars =     [ 'key', 'minimum_snr', 'fwhm', 'artifacts', 'absindex', 'bandhead' ]
        values =   [   key,   minimum_snr,   fwhm,   artifacts,   absindex,   bandhead ]
        dtypes =   [   str,         in_fl,  in_fl,         str,        str,        str ]

        ParSet.__init__(self, pars, values=values, dtypes=dtypes)


def validate_spectral_indices_config(cnfg):
    """ 
    Validate the `configparser.ConfigParser`_ object that is meant to
    define a set of spectral-index measurements.

    Args:

        cnfg (`configparser.ConfigParser`_): Object meant to contain
            defining parameters of the binning method as needed by
            :class:`mangadap.proc.spectralindices.SpectralIndicesDef`

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
    if 'resolution_fwhm' not in cnfg.options('default') \
            or cnfg['default']['resolution_fwhm'] is None:
        cnfg['default']['resolution_fwhm']= '-1'

    if 'artifact_mask' not in cnfg.options('default') \
            or cnfg['default']['artifact_mask'] is None:
        cnfg['default']['artifact_mask'] = 'None'

    if 'absorption_indices' not in cnfg.options('default') \
            or cnfg['default']['absorption_indices'] is None:
        cnfg['default']['absorption_indices'] = 'None'

    if 'bandhead_indices' not in cnfg.options('default') \
            or cnfg['default']['bandhead_indices'] is None:
        cnfg['default']['bandhead_indices'] = 'None'

    if cnfg['default']['absorption_indices'] == 'None' \
            and cnfg['default']['bandhead_indices'] == 'None':
        raise ValueError('Must provide either an absorption-index database or a bandhead-index ' \
                         'database, or both.')


def available_spectral_index_databases(dapsrc=None):
    """
    Return the list of available spectral index databases

    Available database combinations:

    .. todo::
        Fill in

    Args:
        dapsrc (str): (**Optional**) Root path to the DAP source
            directory.  If not provided, the default is defined by
            :func:`mangadap.config.defaults.default_dap_source`.

    Returns:
        list: A list of
        :func:`mangadap.proc.spectralindices.SpectralIndicesDef`
        objects, each defining a set of spectral index databases.

    Raises:
        NotADirectoryError: Raised if the provided or default
            *dapsrc* is not a directory.
        OSError/IOError: Raised if no binning scheme configuration files
            could be found.
        KeyError: Raised if the binning method keywords are not all
            unique.
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
    search_dir = os.path.join(dapsrc, 'python/mangadap/config/spectral_indices')
    ini_files = glob.glob(os.path.join(search_dir, '*.ini'))
    if len(ini_files) == 0:
        raise IOError('Could not find any configuration files in {0} !'.format(search_dir))

    # Build the list of library definitions
    index_set_list = []
    for f in ini_files:
        # Read the config file
        cnfg = ConfigParser(allow_no_value=True)
        cnfg.read(f)
        # Ensure it has the necessary elements to define the template
        # library
        validate_spectral_indices_config(cnfg)

        index_set_list += [ SpectralIndicesDef(cnfg['default']['key'],
                                               eval(cnfg['default']['minimum_snr']), 
                                               eval(cnfg['default']['resolution_fwhm']),
                                               None if cnfg['default']['artifact_mask'] == 'None' \
                                                    else cnfg['default']['artifact_mask'],
                                        None if cnfg['default']['absorption_indices'] == 'None' \
                                                else cnfg['default']['absorption_indices'],
                                        None if cnfg['default']['bandhead_indices'] == 'None' \
                                                else cnfg['default']['bandhead_indices']) ]

    # Check the keywords of the libraries are all unique
    if len(numpy.unique(numpy.array([index['key'] for index in index_set_list]))) \
                != len(index_set_list):
        raise KeyError('Spectral-index set keywords are not all unique!')

    # Return the default list of assessment methods
    return index_set_list


class SpectralIndicesBitMask(BitMask):
    r"""

    Derived class that specifies the mask bits for the spectral-index
    measurements.  See :class:`mangadap.util.bitmask.BitMask` for
    attributes.

    A list of the bits and meanings are provided by the base class
    function :func:`mangadap.util.bitmask.BitMask.info`; i.e.,::

        from mangadap.proc.spectralindices import SpectralIndicesBitMask
        bm = SpectralIndicesBitMask()
        bm.info()

    """
    def __init__(self, dapsrc=None):
        dapsrc = default_dap_source() if dapsrc is None else str(dapsrc)
        BitMask.__init__(self, ini_file=os.path.join(dapsrc, 'python', 'mangadap', 'config',
                                                     'bitmasks', 'spectral_indices_bits.ini'))


class SpectralIndices:
    r"""

    Class that holds the spectral-index measurements.

    """
    def __init__(self, database_key, binned_spectra, redshift=None, stellar_continuum=None,
                 database_list=None, artifact_list=None, absorption_index_list=None,
                 bandhead_index_list=None, dapsrc=None, dapver=None, analysis_path=None,
                 directory_path=None, output_file=None, hardcopy=True, clobber=False, verbose=0,
                 checksum=False):

        self.version = '1.0'
        self.verbose = verbose

        # Define the database properties
        self.database = None
        self.artdb = None
        self.absdb = None
        self.bhddb = None
        self.pixelmask = None
        self._define_databases(database_key, database_list=database_list,
                               artifact_list=artifact_list,
                               absorption_index_list=absorption_index_list,
                               bandhead_index_list=bandhead_index_list)

        self.binned_spectra = None
        self.redshift = None
        self.stellar_continuum = None
        self.nindx = 0
        if self.absdb is not None:
            self.nindx += self.absdb.nsets
        if self.bhddb is not None:
            self.nindx += self.bhddb.nsets

        # Define the output directory and file
        self.directory_path = None      # Set in _set_paths
        self.output_file = None
        self.hardcopy = None

        # Initialize the objects used in the assessments
        self.bitmask = SpectralIndicesBitMask(dapsrc=dapsrc)
        self.hdu = None
        self.checksum = checksum

        self.shape = None
        self.spatial_shape = None
        self.nspec = None
        self.spatial_index = None
        self._assign_spectral_arrays()
        self.dispaxis = None
        self.nwave = None

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
                          absorption_index_list=None, bandhead_index_list=None, dapsrc=None):
        r"""

        Select the database of indices

        """
        # Grab the specific database
        self.database = _select_proc_method(database_key, SpectralIndicesDef,
                                            method_list=database_list,
                                            available_func=available_spectral_index_databases,
                                            dapsrc=dapsrc)

        # Instantiate the artifact, absorption-index, and bandhead-index
        # databases
        self.artdb = None if self.database['artifacts'] is None else \
                ArtifactDB(self.database['artifacts'], artdb_list=artifact_list, dapsrc=dapsrc)
        # TODO: Generalize the name of this object
        self.pixelmask = StellarContinuumPixelMask(self.artdb, None)

        self.absdb = None if self.database['absindex'] is None else \
                AbsorptionIndexDB(self.database['absindex'], indxdb_list=absorption_index_list,
                                  dapsrc=dapsrc)

        self.bhddb = None if self.database['bandhead'] is None else \
                BandheadIndexDB(self.database['bandhead'], indxdb_list=bandhead_index_list,
                                dapsrc=dapsrc)


    def _set_paths(self, directory_path, dapver, analysis_path, output_file):
        """
        Set the :attr:`directory_path` and :attr:`output_file`.  If not
        provided, the defaults are set using, respectively,
        :func:`mangadap.config.defaults.default_dap_reference_path` and
        :func:`mangadap.config.defaults.default_dap_file_name`.

        Args:
            directory_path (str): The exact path to the DAP reduction
                assessments file.  See :attr:`directory_path`.
            dapver (str): DAP version.
            analysis_path (str): The path to the top-level directory
                containing the DAP output files for a given DRP and DAP
                version.
            output_file (str): The name of the file with the reduction assessments.
                See :func:`compute`.

        """
        # Set the output directory path
        self.directory_path = default_dap_reference_path(plate=self.binned_spectra.drpf.plate,
                                                         drpver=self.binned_spectra.drpf.drpver,
                                                         dapver=dapver,
                                                         analysis_path=analysis_path) \
                                        if directory_path is None else str(directory_path)

        # Set the output file
        method = '{0}-{1}'.format(self.binned_spectra.rdxqa.method['key'],
                                  self.binned_spectra.method['key'])
        if self.stellar_continuum is not None:
            method = '{0}-{1}'.format(method, self.stellar_continuum.method['key'])
        method = '{0}-{1}'.format(method, self.database['key'])

        self.output_file = default_dap_file_name(self.binned_spectra.drpf.plate,
                                                 self.binned_spectra.drpf.ifudesign,
                                                 self.binned_spectra.drpf.mode, method) \
                                        if output_file is None else str(output_file)


    def _clean_drp_header(self, ext='FLUX'):
        """
        Read and clean the header from the DRP fits file extension for
        use in the output file for this class object.

        .. todo::
            - Currently just returns the existing header.  Need to
              decide if/how to manipulate DRP header.

        """
        return self.binned_spectra.drpf[ext].header


    def _initialize_header(self, hdr):
        """

        Initialize the header.

        """
        hdr['SIKEY'] = (self.database['key'], 'Spectral-index database keyword')
        hdr['SIMINSN'] = (self.database['minimum_snr'], 'Minimum S/N of spectrum to include')
        hdr['FWHM'] = (self.database['fwhm'], 'FWHM of index system resolution (ang)')
        hdr['ARTDB'] = (self.database['artifacts'], 'Artifact database keyword')
        hdr['ABSDB'] = (self.database['absindex'], 'Absorption-index database keyword')
        hdr['BHDDB'] = (self.database['bandhead'], 'Bandhead-index database keyword')
        if self.stellar_continuum is not None:
            hdr['SCTPL'] = (self.stellar_continuum.method['key'], 'Stellar-continuum model keyword')
        hdr['NBINS'] = (self.nbins, 'Number of unique spatial bins')
        if len(self.missing_bins) > 0:
            hdr['EMPTYBIN'] = (str(self.missing_bins), 'List of bins with no data')
        # Anything else?
        # Additional database details?
        return hdr


    def _initialize_mask(self, good_snr):
        """

        Initialize the mask be setting the DIDNOTUSE, FORESTAR, and LOW_SNR masks

        """
        # Initialize to all zeros
        mask = numpy.zeros(self.shape, dtype=self.bitmask.minimum_uint_dtype())

        # Turn on the flag stating that the pixel wasn't used
        indx = self.binned_spectra.bitmask.flagged(self.binned_spectra['MASK'].data,
                                                   flag=self.binned_spectra.do_not_fit_flags())
        mask[indx] = self.bitmask.turn_on(mask[indx], 'DIDNOTUSE')

        # Turn on the flag stating that the pixel has a foreground star
        indx = self.binned_spectra.bitmask.flagged(self.binned_spectra['MASK'].data,
                                                   flag='FORESTAR')
        mask[indx] = self.bitmask.turn_on(mask[indx], 'FORESTAR')

        # Turn on the flag stating that the S/N in the spectrum was
        # below the requested limit
        indx = numpy.array([numpy.invert(good_snr).reshape(self.spatial_shape).T]*self.nwave).T
        mask[indx] = self.bitmask.turn_on(mask[indx], flag='LOW_SNR')

        return mask


    def _index_database_dtype(self, name_len):
        r"""
        Construct the record array data type for the output fits
        extension.
        """
        return [ ('TYPE','<U10'),
                 ('ID',numpy.int),
                 ('NAME','<U{0:d}'.format(name_len)),
                 ('PASSBAND', numpy.float, 2),
                 ('BLUEBAND', numpy.float, 2),
                 ('REDBAND', numpy.float, 2),
                 ('UNIT', '<U3'),
                 ('COMPONENT', numpy.uint8),
                 ('INTEGRAND', '<U7'),
                 ('ORDER', '<U3')
               ]

    
    def _compile_database(self):
        """
        Compile the database with the specifications of each index.
        """
        name_len = 0
        for n in self.absdb['name']:
            if name_len < len(n):
                name_len = len(n)

        # Instatiate the table data that will be saved defining the set
        # of indices measured
        hdu_database = init_record_array(self.nindx, self._index_database_dtype(name_len))

        hdu_database['TYPE'][:self.absdb.nsets] = 'absorption'
        hk = [ 'ID', 'NAME', 'PASSBAND', 'BLUEBAND', 'REDBAND', 'UNIT', 'COMPONENT' ]
        ak = [ 'index', 'name', 'primary', 'blueside', 'redside', 'units', 'component' ]
        for _hk, _ak in zip(hk,ak):
            hdu_database[_hk][:self.absdb.nsets] = self.absdb[_ak]

        hdu_database['TYPE'][self.absdb.nsets:] = 'bandhead'
        hk = [ 'ID', 'NAME', 'BLUEBAND', 'REDBAND', 'INTEGRAND', 'ORDER' ]
        ak = [ 'index', 'name', 'blueside', 'redside', 'integrand', 'order' ]
        for _hk, _ak in zip(hk,ak):
            hdu_database[_hk][self.absdb.nsets:] = self.bhddb[_ak]

        return hdu_database


    def _per_bin_dtype(self):
        r"""
        Construct the record array data type for the output fits
        extension.
        """
        return [ ('BIN_INDEX',numpy.int),
                 ('REDSHIFT', numpy.float),
                 ('MASK', self.bitmask.minimum_uint_dtype(), self.nindx),
                 ('BCONT', numpy.float, self.nindx), 
                 ('BCONTERR', numpy.float, self.nindx),
                 ('RCONT', numpy.float, self.nindx), 
                 ('RCONTERR', numpy.float, self.nindx), 
                 ('INDX_DISPCORR', numpy.float, self.nindx), 
                 ('INDX', numpy.float, self.nindx), 
                 ('INDXERR', numpy.float, self.nindx)
               ]


    def _assign_spectral_arrays(self):
        self.spectral_arrays = [ 'FLUX', 'IVAR', 'MASK' ]


    def _check_snr(self):
        # binned_spectra['BINS'].data['SNR'] has length
        # binned_spectra.nbins
        return self.binned_spectra['BINS'].data['SNR'] > self.database['minimum_snr']

    
    def _missing_flags(self):
        return numpy.array([ b in self.missing_bins for b in numpy.arange(self.nbins)])


    def _bins_to_measure(self):
        """
        Determine which bins to use for the spectral-index measurements.
        They must not be designated as "missing" and they must have
        sufficient S/N.
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
                print('single')
                self.redshift[:] = redshift
            elif len(_redshift) == self.nspec:
                print('nspec')
                self.redshift = _redshift[ self.unique_bins(index=True) ]
            else:   # Has length nbins
                print('nbin')
                self.redshift = _redshift
        elif self.stellar_continuum is not None:
            print('from stellar continuum')
            self.redshift = self.stellar_continuum['PAR'].data['KIN'][:,0] \
                                / astropy.constants.c.to('km/s').value

        self.redshift = self.redshift[self._bins_to_measure()]
        print(self.redshift)


    def _spectra_for_measurements(self):
        """

        Compile the set of spectra for the spectral-index measurements.
        If the input fwhm is > 0, this function will match the spectral
        resolution of the data to the spectral-index system, based on
        the provided FWHM.

        """
        # Get the data arrays
        wave = self.binned_spectra['WAVE'].data
        flux = self.binned_spectra.copy_to_masked_array(
                                                        flag=self.binned_spectra.do_not_fit_flags())
        ivar = self.binned_spectra.copy_to_masked_array(ext='IVAR',
                                                        flag=self.binned_spectra.do_not_fit_flags())
        indx = self.pixelmask.boolean(wave,nspec=self.nbins)
        flux[indx] = numpy.ma.masked
        ivar[indx] = numpy.ma.masked

        mask = numpy.zeros((self.nbins,self.nwave), dtype=self.bitmask.minimum_uint_dtype())
        flux_mask = numpy.ma.getmaskarray(flux)
        mask[flux_mask] = self.bitmask.turn_on(mask[flux_mask], 'DIDNOTUSE')

        good_bins = self._bins_to_measure()
        flux = flux[good_bins,:]
        ivar = ivar[good_bins,:]
        mask = mask[good_bins,:]

        if self.database['fwhm'] < 0:
            return flux, ivar, mask

#        pyplot.step(wave, flux[0,:], where='mid', linestyle='-', color='k', lw=0.5)
#        pyplot.show()

        # Revert flux and ivar to unmasked arrays
        flux = numpy.asarray(flux)
        ivar = numpy.asarray(ivar)

        # For masked pixels, interpolate across them using simple linear
        # interpolation.
        # TODO: This approach needs to be vetted.
        for i in range(flux.shape[0]):
            mp = mask[i,:] > 0              # Masked pixels
            up = numpy.invert(mp)           # Unmasked pixels
            if numpy.sum(mp) == 0:          # No masked pixels so continue
                continue
            if numpy.sum(up) == 0:          # No unmasked pixels!
                flux[i,:] = 0.0
                continue
            # Linearly interpolate both the flux and inverse variance
            flux[i,mp] = numpy.interp(wave[mp], wave[up], flux[i,up])
            ivar[i,mp] = numpy.interp(wave[mp], wave[up], ivar[i,up])
#        pyplot.step(wave, flux[0,:], where='mid', linestyle='-', color='k', lw=0.5)
#        pyplot.show()

        # Use :func:`mangadap.instrument.match_spectral_resolution` to
        # match the spectral resolution of the binned spectra to the
        # spectral-index system
        existing_sres = self.binned_spectra.drpf['SPECRES'].data
        new_sres = wave/self.database['fwhm']
        new_flux, sres, sigoff, new_mask, new_ivar \
                = match_spectral_resolution(wave, flux, existing_sres, wave, new_sres, ivar=ivar,
                                            log10=True, new_log10=True)

#        pyplot.step(wave, flux[0,:], where='mid', linestyle='-', color='k', lw=0.5)
#        pyplot.step(wave, new_flux[0,:], where='mid', linestyle='-', color='r', lw=2.5)
#        pyplot.show()
#        pyplot.step(wave, ivar[0,:]/numpy.mean(ivar[0,:]), where='mid', linestyle='-', color='k',
#                    lw=0.5)
#        pyplot.step(wave, new_ivar[0,:]/numpy.mean(new_ivar[0,:]), where='mid', linestyle='-', color='r', lw=2.5)
#        pyplot.show()

        # From instrument.py:  "Any pixel that had a resolution that was
        # lower than the target resolution (up to some tolerance defined
        # by *min_sig_pix*) is returned as masked."  This is exactly
        # what should be masked with the SPECRES_LOW bit.
        if numpy.sum(new_mask) > 0:
            mask[new_mask > 0] = self.bitmask.turn_on(mask[new_mask > 0], 'SPECRES_LOW')
        return numpy.ma.MaskedArray(flux,mask=mask>0), numpy.ma.MaskedArray(ivar,mask=mask>0), mask


    def _bandhead_indices(self, hdu_measurements):
        """
        Calculate the bandhead indices based on the measure
        pseudo-continuua in the blue and red sidbands.
        """
        # Initialize the arrays to hold the numerator and denominator
        n = numpy.ma.zeros(self.bhddb.nsets)
        nerr = numpy.ma.zeros(self.bhddb.nsets)
        d = numpy.ma.zeros(self.bhddb.nsets)
        derr = numpy.ma.zeros(self.bhddb.nsets)

        # Calculate the index components
        for i in range(self.bhddb.nsets):
            nkey, dkey = ('BCONT','RCONT') if self.bhddb['order'][i] == 'b_r' else ('RCONT','BCONT')
            n[i] = hdu_measurements[nkey][i+self.absdb.nsets]
            nerr[i] = hdu_measurements['{0}ERR'.format(nkey)][i+self.absdb.nsets]
            d[i] = hdu_measurements[dkey][i+self.absdb.nsets]
            derr[i] = hdu_measurements['{0}ERR'.format(dkey)][i+self.absdb.nsets]

        # Determine which indices have both a valid index and index
        # error calculation
        # TODO: Do something different for the errors
        indx = (numpy.absolute(d) > 0) & (numpy.absolute(n) > 0)
        if numpy.sum(indx) != self.bhddb.nsets:
            hdu_measurements['MASK'][self.absdb.nsets:][numpy.invert(indx)] \
                = self.bitmask.turn_on(
                    hdu_measurements['MASK'][self.absdb.nsets:][numpy.invert(indx)], 'DIVBYZERO')

        # Calculate the index and flag any with division by zero issues
        hdu_measurements['INDX'][self.absdb.nsets:][indx] = n[indx]/d[indx]
        hdu_measurements['INDXERR'][self.absdb.nsets:][indx] = numpy.sqrt(
                numpy.square(nerr[indx]*hdu_measurements['INDX'][self.absdb.nsets:][indx]/n[indx])
              + numpy.square(derr[indx]*hdu_measurements['INDX'][self.absdb.nsets:][indx]/d[indx]) )



    def _absorption_indices(self, wave, flux, hdu_measurements, err=None):

        # Get the redshifted passband centers
        # TODO: This does NOT account for masked pixels...
        blue_center = numpy.mean((1.0+hdu_measurements['REDSHIFT'])*self.absdb['blueside'], axis=1)
        red_center = numpy.mean((1.0+hdu_measurements['REDSHIFT'])*self.absdb['redside'], axis=1)
        # Get the parameters for the linear continuum across the
        # primary passband
        continuum_m = (hdu_measurements['RCONT'][:self.absdb.nsets] 
                        - hdu_measurements['BCONT'][:self.absdb.nsets]) \
                      / (red_center - blue_center)[:self.absdb.nsets]
        continuum_b = hdu_measurements['BCONT'][:self.absdb.nsets] - blue_center*continuum_m

        # Get the redshifted passband and its width
        primaryband = (1.0+hdu_measurements['REDSHIFT'])*self.absdb['primary']
        primaryband_width = numpy.diff(primaryband, axis=1).reshape(primaryband.shape[0])

        # Compute the continuum normalize indices.  This has to be a for
        # loop because the continuum is index-dependent
        for i,p in enumerate(primaryband):
            # From Worthey et al. 1994, eqns. 2 and 3
            cont = continuum_b[i] + continuum_m[i]*wave
            integrand = 1.0 - flux/cont if self.absdb['units'][i] == 'ang' else flux/cont

            # Calculate the integral over the passband
            hdu_measurements['INDX'][i] = passband_integral(wave, integrand, passband=p, log=True)
            if err is not None:
                hdu_measurements['INDXERR'][i] = passband_integral(wave, numpy.square(err/cont),
                                                                   passband=p, log=True)

            # Get the fraction of the band covered by the spectrum and
            # flag bands that are only partially covered or empty
            interval = passband_integrated_width(wave, integrand, passband=p, log=True)
            interval_frac = interval / primaryband_width[i]
            if interval_frac < 1.0:
                hdu_measurements['MASK'][i] = self.bitmask.turn_on(hdu_measurements['MASK'][i],
                                                                   'MAIN_INCOMP')
            if not interval_frac > 0.0:
                hdu_measurements['MASK'][i] = self.bitmask.turn_on(hdu_measurements['MASK'][i],
                                                                   'MAIN_EMPTY')
            # Convert the index to magnitudes
            if self.absdb['units'][i] == 'mag':
                if not interval > 0.0:
                    hdu_measurements['MASK'][i] = self.bitmask.turn_on(hdu_measurements['MASK'][i],
                                                                       'DIVBYZERO')
                    hdu_measurements['INDX'][i] = 0.0
                    hdu_measurements['INDXERR'][i] = 0.0
                else:
                    # Worthey et al. 1994, eqn. 3: The passband interval
                    # cancels out of the error propagation.  The error
                    # calculation is done first so as to not replace the
                    # linear calculation of the index.
                    hdu_measurements['INDXERR'][i] = numpy.absolute(2.5
                                                    * hdu_measurements['INDXERR'][i]
                                                    / hdu_measurements['INDX'][i]/numpy.log(10.0))
                    hdu_measurements['INDX'][i] = -2.5 * numpy.ma.log10(hdu_measurements['INDX'][i] 
                                                                        / interval)


    def _measure_indices(self, redshift, flux, ivar=None, mask=None):
        """
        Measure the indices.

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
            
        # Get the list of bands that integrate over f_nu, not f_lambda,
        # and create the f_nu spectra
        fnu_indx = numpy.append(numpy.append(self.absdb['integrand'],
                                             self.bhddb['integrand'],axis=0),
                                numpy.append(self.absdb['integrand'],
                                             self.bhddb['integrand'],axis=0), axis=0) == 'fnu'
        fnu_bands = numpy.sum(fnu_indx) > 0
        if fnu_bands:
            flux_fnu = flux_to_fnu(numpy.array([wave]*nspec), _flux)

#        nu = 1e-13*astropy.constants.c.to('nm/s').value/wave    # In THz
#        pyplot.plot(wave, _flux[0,:]*wave, linestyle='-', color='r', lw=1, zorder=2)
#        pyplot.plot(wave, flux_fnu[0,:]*nu, linestyle='-', color='0.5', lw=3, zorder=1)
#        pyplot.show()

        # Measure the pseudo-continuum in the sidebands, working under
        # the assumption that the number of f_nu integrals to be far
        # less than the f_lambda integrals.
        sidebands = numpy.append(numpy.append(self.absdb['blueside'],
                                              self.bhddb['blueside'],axis=0),
                                 numpy.append(self.absdb['redside'],
                                              self.bhddb['redside'],axis=0), axis=0)
        for i in range(nspec):
            # Shift the sidebands to the appropriate redshift
            _sidebands = sidebands*(1.0+redshift[i])
            _noise = None if noise is None else noise[i,:]

            # Calculate the pseudo-continuum
            pseudocontinuum, pseudocontinuum_error \
                    = passband_integrated_mean(wave, _flux[i,:], _sidebands, err=_noise, log=True)
            # Calculate the pseudo-continuum (for integrals over f_nu)
            if fnu_bands:
                pseudocontinuum[fnu_indx], fnu_pseudocontinuum_error \
                        = passband_integrated_mean(wave, flux_fnu[i,:], _sidebands[fnu_indx],
                                                   err=_noise, log=True)
                if noise is not None:
                    pseudocontinuum_error[fnu_indx] = fnu_pseudocontinuum_error
                

            # Save the results to the main output array
            hdu_measurements['BCONT'][i,:] = pseudocontinuum.reshape(-1,self.nindx)[0,:]
            hdu_measurements['RCONT'][i,:] = pseudocontinuum.reshape(-1,self.nindx)[1,:]
            if pseudocontinuum_error is not None:
                hdu_measurements['BCONTERR'][i,:] \
                        = pseudocontinuum_error.reshape(-1,self.nindx)[0,:]
                hdu_measurements['RCONTERR'][i,:] \
                        = pseudocontinuum_error.reshape(-1,self.nindx)[1,:]

            # Check for the existence of valid pixels in the sideband
            interval_frac = passband_integrated_width(wave, _flux[i,:], _sidebands, log=True) \
                                    / numpy.diff(_sidebands, axis=1).ravel()
            # Passband is not fully covered by valid pixels
            blue_fraction = interval_frac.reshape(-1,self.nindx)[0,:]
            hdu_measurements['MASK'][i,blue_fraction < 1.0] \
                    = self.bitmask.turn_on(hdu_measurements['MASK'][i,blue_fraction < 1.0],
                                           'BLUE_INCOMP')
            red_fraction = interval_frac.reshape(-1,self.nindx)[1,:]
            hdu_measurements['MASK'][i,red_fraction < 1.0] \
                    = self.bitmask.turn_on(hdu_measurements['MASK'][i,red_fraction < 1.0],
                                           'RED_INCOMP')
            # Passband is empty
            hdu_measurements['MASK'][i,~(blue_fraction > 0.0)] \
                    = self.bitmask.turn_on(hdu_measurements['MASK'][i,~(blue_fraction > 0.0)],
                                           'BLUE_EMPTY')
            hdu_measurements['MASK'][i,~(red_fraction > 0.0)] \
                    = self.bitmask.turn_on(hdu_measurements['MASK'][i,~(red_fraction > 0.0)],
                                           'RED_EMPTY')

#            pyplot.scatter(_interval, interval/_interval, marker='.', color='k', s=50)
#            pyplot.show()

#            center = numpy.mean(_sidebands, axis=1).reshape(-1,self.nindx)
#            pyplot.step(wave, flux[i,:], where='mid', color='k', lw=0.5, linestyle='-')
#            pyplot.scatter(center[0,:], hdu_measurements['BCONT'][i,:], marker='.',
#                           s=100, color='b', lw=0)
#            pyplot.scatter(center[1,:], hdu_measurements['RCONT'][i,:], marker='.',
#                           s=100, color='r', lw=0)
#            pyplot.show()

            # Perform the bandhead measurements
            self._bandhead_indices(hdu_measurements[i])

            # Perform the absorption-line measurements
            self._absorption_indices(wave, flux[i,:], hdu_measurements[i], err=_noise)

#        x = numpy.append(center, numpy.array([[-100]*self.absdb.nsets]).T,axis=1).ravel()
#        _x = numpy.ma.MaskedArray(x, mask=x==-100)
#        print(_x)
#        y = numpy.append(pseudocontinuum, numpy.array([[-100]*self.absdb.nsets]).T,axis=1).ravel()
#        _y = numpy.ma.MaskedArray(y, mask=y==-100)
#        print(type(_y))
#        print(numpy.sum(_y.mask))
#        pyplot.plot(_x, _y, color='g', lw=2, linestyle='-')
#        pyplot.show()

        return hdu_measurements
        

    def _unit_selection(self):
        
        # Flag the indices as either having magnitude or angstrom units
        angu = numpy.full(self.nindx, False, dtype=numpy.bool)
        angu[:self.absdb.nsets] = self.absdb['units'] == 'ang'
        magu = numpy.invert(angu)
        # Bandhead indices are unitless, but the correction follow the
        # same operation as for the absorption-line indices with
        # angstrom units
        magu[self.absdb.nsets:] = False

        return angu, magu


    def _resolution_matched_template_library(self, dapsrc=None):
        """
        Get a version of the template library that has had its
        resolution matched to that of the spectral-index database.
        """
        wave = self.binned_spectra['WAVE'].data
        # Set the spectral resolution
        sres = spectral_resolution(wave, wave/self.database['fwhm'], log10=True)
        velocity_offset=self.stellar_continuum.method['fitpar']['template_library'].velocity_offset
        # Set the file designation for this template library.
        # Move this to config.default?
        designation = '{0}-FWHM{1:.2f}'.format(
                                self.stellar_continuum.method['fitpar']['template_library_key'],
                                self.database['fwhm'])
        # Set the output file name
        processed_file = default_dap_file_name(self.binned_spectra.drpf.plate,
                                               self.binned_spectra.drpf.ifudesign,
                                               self.binned_spectra.drpf.mode, designation)
        # Return the template library object
        return TemplateLibrary(self.stellar_continuum.method['fitpar']['template_library_key'],
                               sres=sres, velocity_offset=velocity_offset,
                               spectral_step=spectral_coordinate_step(wave, log=True), log=True,
                               dapsrc=dapsrc, directory_path=self.directory_path,
                               processed_file=processed_file)
    

    def _calculate_dispersion_corrections(self, redshift, dapsrc=None):
        if self.stellar_continuum is None:
            return None

        # Get the wavelength and mask arrays
        wave = self.binned_spectra['WAVE'].data
        flux = self.binned_spectra.copy_to_masked_array(
                                                    flag=self.binned_spectra.do_not_fit_flags())
        mask = numpy.ma.getmaskarray(flux) | self.pixelmask.boolean(wave, nspec=self.nbins)

        print(self.database['fwhm'])
        print(self.stellar_continuum.method['fitpar']['template_library_key'])

        # Get the template library
        template_library = None if self.database['fwhm'] < 0  \
                or self.stellar_continuum.method['fitpar']['template_library_key'] is None else \
                    self._resolution_matched_template_library(dapsrc=dapsrc)

        # Get the broadened stellar-continuum model
        good_bins = self._bins_to_measure()
        nspec = numpy.sum(good_bins)
        broadened_models = numpy.ma.MaskedArray(self.stellar_continuum.construct_models(
                                                template_library=template_library), mask=mask)
        broadened_models = broadened_models[good_bins,:]

        # Get the unbroadened version
        unbroadened_models = numpy.ma.MaskedArray(self.stellar_continuum.construct_models(
                                                  template_library=template_library,
                                                  redshift_only=True), mask=flux.mask.copy())
        unbroadened_models = unbroadened_models[good_bins,:]

#        gflux = self.binned_spectra.copy_to_masked_array(
#                                                    flag=self.binned_spectra.do_not_fit_flags())
#        gflux[mask] = numpy.ma.masked
#        gflux = gflux[good_bins,:]
#
#        flux = self.binned_spectra.copy_to_masked_array(
#                                                    flag=self.binned_spectra.do_not_fit_flags())
#        pyplot.step(wave, flux[0,:], where='mid', linestyle='-', color='k', lw=0.5, zorder=3)
#        pyplot.step(wave, broadened_models[0,:], where='mid', linestyle='-', color='r', lw=1.5,
#                    zorder=1, alpha=0.5)
#        pyplot.step(wave, unbroadened_models[0,:], where='mid', linestyle='-', color='g', lw=1.5,
#                    zorder=2, alpha=0.5)
#        pyplot.show()

        # Measure the indices for both sets of spectra
        broadened_model_measurements = self._measure_indices(redshift, broadened_models)
        unbroadened_model_measurements = self._measure_indices(redshift, unbroadened_models)

        # Do not apply the correction if any of the bands were empty or
        # would have had to divide by zero
        bad = self.bitmask.flagged(broadened_model_measurements['MASK'],
                                   flag=[ 'MAIN_EMPTY', 'BLUE_EMPTY', 'RED_EMPTY', 'DIVBYZERO' ])
        bad |= self.bitmask.flagged(unbroadened_model_measurements['MASK'],
                                    flag=[ 'MAIN_EMPTY', 'BLUE_EMPTY', 'RED_EMPTY', 'DIVBYZERO' ])
        angu, magu = self._unit_selection()

        good_ang = ~bad & (numpy.array([angu]*nspec)) \
                        & (numpy.absolute(broadened_model_measurements['INDX']) > 0)
        good_mag = ~bad & (numpy.array([magu]*nspec))

        corrections = numpy.zeros((nspec, self.nindx), dtype=numpy.float)
        corrections[good_ang] = unbroadened_model_measurements['INDX'][good_ang] \
                                    / broadened_model_measurements['INDX'][good_ang]
        corrections[good_mag] = unbroadened_model_measurements['INDX'][good_mag] \
                                    - broadened_model_measurements['INDX'][good_mag]
        return corrections, good_ang, good_mag


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

        Measure the spectral indices using the binned spectra.

        """
        # SpatiallyBinnedSpectra object always needed
        if binned_spectra is None:
            raise ValueError('Must provide spectra object for fitting.')
        if not isinstance(binned_spectra, SpatiallyBinnedSpectra):
            raise TypeError('Must provide a valid SpatiallyBinnedSpectra object!')
        if binned_spectra.hdu is None:
            raise ValueError('Provided SpatiallyBinnedSpectra object is undefined!')
        self.binned_spectra = binned_spectra

        # TODO: Will need the emission-line model data as well...
        
        # StellarContinuumModel object only used when calculating dispersion corrections.
        if stellar_continuum is not None:
            if not isinstance(stellar_continuum, StellarContinuumModel):
                raise TypeError('Provided stellar continuum must have StellarContinuumModel type!')
            if stellar_continuum.hdu is None:
                raise ValueError('Provided StellarContinuumModel is undefined!')
            self.stellar_continuum = stellar_continuum

        self.shape = self.binned_spectra.shape
        self.spatial_shape =self.binned_spectra.spatial_shape
        self.nspec = self.binned_spectra.nspec
        self.spatial_index = self.binned_spectra.spatial_index.copy()
        self.dispaxis = self.binned_spectra.dispaxis
        self.nwave = self.binned_spectra.nwave
        
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
        flux, ivar, mask = self._spectra_for_measurements()

        # Compile the information on the suite of measured indices
        hdu_database = self._compile_database()

        # Perform the measurements on the galaxy spectra
        hdu_measurements = self._measure_indices(self.redshift, flux, ivar=ivar, mask=mask)

        # Get the corrections by performing the measurements on the
        # best-fitting continuum models, with and without the velocity
        # dispersion broadening
        hdu_measurements['INDX_DISPCORR'], good_ang, good_mag \
                = self._calculate_dispersion_corrections(self.redshift, dapsrc=dapsrc)

        # Flag bad corrections
        angu, magu = self._unit_selection()
        bad_ang = (numpy.array([angu]*good_ang.shape[0])) & ~good_ang
        hdu_measurements['MASK'][bad_ang] \
                = self.bitmask.turn_on(hdu_measurements['MASK'][bad_ang],
                                       'NO_DISPERSION_CORRECTION')
        bad_mag = (numpy.array([magu]*good_mag.shape[0])) & ~good_mag
        hdu_measurements['MASK'][bad_mag] \
                = self.bitmask.turn_on(hdu_measurements['MASK'][bad_mag],
                                       'NO_DISPERSION_CORRECTION')

        # Apply the corrections
        hdu_measurements['INDX'][good_ang] *= hdu_measurements['INDX_DISPCORR'][good_ang]
        hdu_measurements['INDXERR'][good_ang] \
                *= numpy.abs(hdu_measurements['INDX_DISPCORR'][good_ang])
        hdu_measurements['INDX'][good_mag] += hdu_measurements['INDX_DISPCORR'][good_mag]

        # Fill area for any missing binsPrep for writing
        if len(hdu_measurements) < self.nbins:
            # Fill the full array for all bins
            _hdu_measurements = hdu_measurements.copy()
            hdu_measurements = init_record_array(self.nbins, self._per_bin_dtype())
            if len(self.missing_bins) > 0:
                hdu_measurements['MASK'][self.missing_bins] \
                    = self.bitmask.turn_on(hdu_measurements['MASK'][self.missing_bins], 'DIDNOTUSE')
            hdu_measurements[good_snr] = _hdu_measurements

            _flux = numpy.ma.MaskedArray(numpy.zeros((self.nbins,self.nwave), dtype=numpy.float))
            _flux[good_snr,:] = flux

            _ivar = numpy.ma.MaskedArray(numpy.zeros((self.nbins,self.nwave), dtype=numpy.float))
            _ivar[good_snr,:] = ivar

            _mask = numpy.zeros((self.nbins,self.nwave), dtype=self.bitmask.minimum_uint_dtype())
            _mask[good_snr,:] = mask
            _mask[~good_snr,:] = self.bitmask.turn_on(_mask[~good_snr,:], 'DIDNOTUSE')
        else:
            _flux = flux.copy()
            _ivar = ivar.copy()
            _mask = mask.copy()

        hdu_measurements['BIN_INDEX'] = numpy.arange(self.nbins)
        hdu_measurements['MASK'][~(good_snr)] \
                = self.bitmask.turn_on(hdu_measurements['MASK'][~(good_snr)], 'LOW_SNR') 

        bin_indx = self.binned_spectra['BINID'].data.ravel()
        unique_bins, reconstruct = numpy.unique(bin_indx, return_inverse=True)
        indx = bin_indx > -1

        flux = numpy.zeros(self.shape, dtype=numpy.float)
        flux.reshape(-1,self.nwave)[indx,:] = _flux[unique_bins[reconstruct[indx]],:]
    
        ivar = numpy.zeros(self.shape, dtype=numpy.float)
        ivar.reshape(-1,self.nwave)[indx,:] = _ivar[unique_bins[reconstruct[indx]],:]
   
        _good_snr = numpy.full(numpy.prod(self.spatial_shape), False, dtype=numpy.bool)
        _good_snr[indx] = good_snr[unique_bins[reconstruct[indx]]]
   
        mask = self._initialize_mask(_good_snr)
        mask.reshape(-1,self.nwave)[indx,:] = _mask[unique_bins[reconstruct[indx]],:]


        # Initialize the header keywords
        hdr = self._clean_drp_header(ext='PRIMARY')
        self._initialize_header(hdr)

        # Save the data to the hdu attribute
        self.hdu = fits.HDUList([ fits.PrimaryHDU(header=hdr),
                                  fits.ImageHDU(data=flux.data,
                                                header=self.binned_spectra['FLUX'].header,
                                                name='FLUX'),
                                  fits.ImageHDU(data=ivar.data,
                                                header=self.binned_spectra['IVAR'].header,
                                                name='IVAR'),
                                  fits.ImageHDU(data=mask,
                                                header=self.binned_spectra['MASK'].header,
                                                name='MASK'),
                                  self.binned_spectra['WAVE'],
                                  fits.BinTableHDU.from_columns( [ fits.Column(name=n,
                                                        format=rec_to_fits_type(hdu_database[n]),
                                array=hdu_database[n]) for n in hdu_database.dtype.names ],
                                                               name='SIPAR'),
                                  fits.BinTableHDU.from_columns( [ fits.Column(name=n,
                                                    format=rec_to_fits_type(hdu_measurements[n]),
                                array=hdu_measurements[n]) for n in hdu_measurements.dtype.names ],
                                                               name='SINDX')
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
        print('restructuring')
        # Restructure the data to match the DRPFits file
        if self.binned_spectra.drpf.mode == 'CUBE':
            DRPFits._restructure_cube(self.hdu, ext=self.spectral_arrays, inverse=True)
        elif self.binned_spectra.drpf.mode == 'RSS':
            DRPFits._restructure_rss(self.hdu, ext=self.spectral_arrays, inverse=True)

        # Get the output file and determine if it should be compressed
        ofile = self.file_path()
        write_hdu(self.hdu, ofile, clobber=clobber, checksum=True)

        # Revert the structure
        print('reverting')
        if self.binned_spectra.drpf.mode == 'CUBE':
            DRPFits._restructure_cube(self.hdu, ext=self.spectral_arrays)
        elif self.binned_spectra.drpf.mode == 'RSS':
            DRPFits._restructure_rss(self.hdu, ext=self.spectral_arrays)


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
        if self.hdu['PRIMARY'].header['SIKEY'] != self.database['key']:
            if strict:
                raise ValueError('Keywords in header does not match specified method keyword!')
            else:
                warnings.warn('Keywords in header does not match specified method keyword!')
        # TODO: "strict" should also check other aspects of the file to
        # make sure that the details of the method are also the same,
        # not just the keyword

        if self.binned_spectra.drpf.mode == 'CUBE':
            DRPFits._restructure_cube(self.hdu, ext=self.spectral_arrays)
        elif self.binned_spectra.drpf.mode == 'RSS':
            DRPFits._restructure_rss(self.hdu, ext=self.spectral_arrays)

        self.nbins = self.hdu['PRIMARY'].header['NBINS']
        try:
            self.missing_bins = eval(self.hdu['PRIMARY'].header['EMPTYBIN'])
        except KeyError:
            # Assume if this fails, it's because the keyword doesn't
            # exist
            self.missing_bins = []

        print(self.missing_bins)










    def unique_bins(self, index=False):
        """
        Get the unique bins or the indices of the unique bins in the
        flattened spatial dimension.
        """
        unique_bins, first_occurrence = numpy.unique(self.binned_spectra['BINID'].data.ravel(),
                                                     return_index=True)
        return first_occurrence[1:] if index else unique_bins[1:]



