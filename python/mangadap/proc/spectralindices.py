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

        import glob
        import os.path
        import numpy
        from astropy.io import fits
        import astropy.constants

        from ..par.parset import ParSet
        from ..par.artifactdb import ArtifactDB
        from ..par.absorptionindexdb import AbsorptionIndexDB
        from ..par.bandheadindexdb import BandheadIndexDB
        from ..config.defaults import default_dap_source, default_dap_file_name
        from ..config.defaults import default_dap_method, default_dap_method_path
        from ..config.defaults import default_dap_common_path
        from ..util.instrument import spectral_resolution, match_spectral_resolution
        from ..util.instrument import spectral_coordinate_step, spectrum_velocity_scale
        from ..util.fitsutil import DAPFitsUtil
        from ..util.fileio import init_record_array, rec_to_fits_type
        from ..util.log import log_output
        from ..util.bitmask import BitMask
        from ..util.pixelmask import SpectralPixelMask
        from ..util.parser import DefaultConfig
        from .spatiallybinnedspectra import SpatiallyBinnedSpectra
        from .templatelibrary import TemplateLibrary
        from .stellarcontinuummodel import StellarContinuumModel
        from .emissionlinemodel import EmissionLineModel
        from .bandpassfilter import passband_integral, passband_integrated_width, passband_integrated_mean
        from .bandpassfilter import passband_weighted_mean
        from .util import select_proc_method, flux_to_fnu


*Class usage examples*:
    Add examples!

*Revision history*:
    | **20 Apr 2016**: Implementation begun by K. Westfall (KBW)
    | **09 May 2016**: (KBW) Add subtraction of emission-line models
    | **11 Jul 2016**: (KBW) Allow to not apply dispersion corrections
        for index measurements
    | **28 Jul 2016**: (KBW) Fixed error in initialization of guess
        redshift when stellar continuum is provided.
    | **23 Feb 2017**: (KBW) Use DAPFitsUtil read and write functions.
    | **27 Feb 2017**: (KBW) Use DefaultConfig

.. todo::

    - Also provide index as measured from best-fitting model without
      dispersion?

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

import glob
import os
import logging

import numpy
from astropy.io import fits
import astropy.constants

from ..par.parset import ParSet
from ..par.artifactdb import ArtifactDB
from ..par.absorptionindexdb import AbsorptionIndexDB
from ..par.bandheadindexdb import BandheadIndexDB
from ..config.defaults import default_dap_source, default_dap_file_name
from ..config.defaults import default_dap_method, default_dap_method_path
from ..config.defaults import default_dap_common_path
from ..util.instrument import spectral_resolution, match_spectral_resolution
from ..util.instrument import spectral_coordinate_step, spectrum_velocity_scale
from ..util.fitsutil import DAPFitsUtil
from ..util.fileio import init_record_array, rec_to_fits_type
from ..util.log import log_output
from ..util.bitmask import BitMask
from ..util.pixelmask import SpectralPixelMask
from ..util.parser import DefaultConfig
from .spatiallybinnedspectra import SpatiallyBinnedSpectra
from .templatelibrary import TemplateLibrary
from .stellarcontinuummodel import StellarContinuumModel
from .emissionlinemodel import EmissionLineModel
from .bandpassfilter import passband_integral, passband_integrated_width, passband_integrated_mean
from .bandpassfilter import passband_weighted_mean, pseudocontinuum
from .util import select_proc_method, flux_to_fnu

from matplotlib import pyplot

# Add strict versioning
# from distutils.version import StrictVersion


class SpectralIndicesDef(ParSet):
    """
    A class that holds the parameters necessary to perform the
    spectral-index measurements.

    Args:
        key (str): Keyword used to distinguish between different
            spectral-index databases.
        minimum_snr (bool): Minimum S/N of spectrum to fit
        fwhm (int, float): Resolution FWHM in angstroms at which to make
            the measurements.
        compute_corrections (bool): Flag to compute dispersion corrections
            to indices.  Dispersion corrections are always calculated!
        artifacts (str): String identifying the artifact database to use
        absindex (str): String identifying the absorption-index database
            to use
        bandhead (str): String identifying the bandhead-index database
            to use
    """
    def __init__(self, key, minimum_snr, fwhm, compute_corrections, artifacts, absindex, bandhead):
        in_fl = [ int, float ]

        pars =     [ 'key', 'minimum_snr', 'fwhm', 'compute_corrections', 'artifacts', 'absindex',
                        'bandhead' ]
        values =   [   key,   minimum_snr,   fwhm,   compute_corrections,   artifacts,   absindex,
                          bandhead ]
        dtypes =   [   str,         in_fl,  in_fl,                  bool,         str,        str,
                               str ]

        ParSet.__init__(self, pars, values=values, dtypes=dtypes)


def validate_spectral_indices_config(cnfg):
    """ 
    Validate :class:`mangadap.util.parser.DefaultConfig` object with
    spectral-index measurement parameters.

    Args:
        cnfg (:class:`mangadap.util.parser.DefaultConfig`): Object with
            parameters to validate.

    Raises:
        KeyError: Raised if any required keywords do not exist.
        ValueError: Raised if keys have unacceptable values.
        FileNotFoundError: Raised if a file is specified but could not
            be found.
    """
    # Check for required keywords
    if not cnfg.keyword_specified('key'):
        raise KeyError('Keyword \'key\' must be provided.')

    if not cnfg.keyword_specified('absorption_indices') \
                and not cnfg.keyword_specified('bandhead_indices'):
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
        list: A list of :func:`SpectralIndicesDef` objects, each
        defining a spectral-index database to measure.

    Raises:
        NotADirectoryError: Raised if the provided or default
            *dapsrc* is not a directory.
        OSError/IOError: Raised if no spectral-index configuration files
            could be found.
        KeyError: Raised if the spectral-index database keywords are not
            all unique.

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
        cnfg = DefaultConfig(f=f)
        # Ensure it has the necessary elements to define the template
        # library
        validate_spectral_indices_config(cnfg)

        index_set_list += [ SpectralIndicesDef(cnfg['key'],
                                               cnfg.getfloat('minimum_snr', default=0.), 
                                               cnfg.getfloat('resolution_fwhm', default=-1),
                                               cnfg.getbool('compute_sigma_correction',
                                                            default=False),
                                               cnfg['artifact_mask'], cnfg['absorption_indices'],
                                               cnfg['bandhead_indices']) ]

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


# TODO: These two should have the same base class
class AbsorptionLineIndices:
    """
    Measure a set of spectral indices.
    """
    def __init__(self, wave, flux, bluebands, redbands, mainbands, err=None, log=True, units=None):

        # Check the shape of the input spectrum
        if len(flux.shape) != 1:
            raise ValueError('Input flux must be a single vector!')
        if wave.shape != flux.shape:
            raise ValueError('Input flux and wavelength vectors must have the same shape.')
        
        if len(bluebands.shape) != 2:
            raise ValueError('Band definitions must be two-dimensional: Nindx x 2.')
        if bluebands.shape != redbands.shape or bluebands.shape != mainbands.shape:
            raise ValueError('Input bands must have identical shape.')

        self.nindx = bluebands.shape[0]
        # Make sure the units exist
        if units is None:
            warnings.warn('Input units not provided.  Assuming angstroms.')
            self.units = numpy.full(self.nindx, 'ang', dtype=object)
        else:
            self.units = units
        self.order = None

        # Get the two pseudocontinua and the flux-weighted band centers
        self.blue_center, self.blue_continuum, self.blue_continuum_err, self.blue_incomplete, \
            self.blue_empty = pseudocontinuum(wave, flux, passband=bluebands, err=err, log=log)

#        print(self.blue_center)
#        print('')
#        print(self.blue_continuum)
#        print('')

        self.red_center, self.red_continuum, self.red_continuum_err, self.red_incomplete, \
            self.red_empty = pseudocontinuum(wave, flux, passband=redbands, err=err, log=log)
#        print(self.red_center)
#        print('')
#        print(self.red_continuum)
#        print('')

        
#       # TEST: Change to center instead of flux-weighted center
#        self.blue_center = numpy.mean(bluebands, axis=1).ravel()
#        self.red_center = numpy.mean(redbands, axis=1).ravel()

        # Get the parameters for the linear continuum across the
        # primary passband
        self.continuum_m = (self.red_continuum - self.blue_continuum) \
                                / (self.red_center - self.blue_center)
        self.continuum_b = self.blue_continuum - self.blue_center * self.continuum_m
#        print(self.continuum_m)
#        print('')
#        print(self.continuum_b)
#        print('')

        # Compute the continuum normalized indices.  This has to be done
        # in a for loop because the continuum is index-dependent
        self.index = numpy.zeros(self.nindx, dtype=numpy.float)
        self.index_err = numpy.zeros(self.nindx, dtype=numpy.float)
        self.main_incomplete = numpy.zeros(self.nindx, dtype=numpy.bool)
        self.main_empty = numpy.zeros(self.nindx, dtype=numpy.bool)
        self.divbyzero = numpy.zeros(self.nindx, dtype=numpy.bool)

        for i,m in enumerate(mainbands):

            # From Worthey et al. 1994, eqns. 2 and 3
            cont = self.continuum_b[i] + self.continuum_m[i]*wave
            integrand = 1.0 - flux/cont if self.units[i] == 'ang' else flux/cont

            # Calculate the integral over the passband
#            print(self.units[i])
#            print(m)
#            dw = numpy.diff(wave[0:2])[0]
#            print(dw)
#            print(self.index[i])
            self.index[i] = passband_integral(wave, integrand, passband=m, log=log)
#            print(self.index[i])
#            tmp = numpy.logical_and(wave > m[0], wave < m[1])
#            print(numpy.sum(integrand[tmp])*dw)
#            tmp = numpy.logical_and(wave > bluebands[i,0], wave < redbands[i,1])
#            pyplot.plot(wave[tmp], integrand[tmp])
#            pyplot.show()
            if err is not None:
                self.index_err[i] = numpy.sqrt(passband_integral(wave, numpy.square(err/cont),
                                               passband=m, log=log))

            # Get the fraction of the band covered by the spectrum and
            # flag bands that are only partially covered or empty
            interval = passband_integrated_width(wave, flux, passband=m, log=log)
            interval_frac = interval / numpy.diff(m)[0]
            self.main_incomplete[i] = interval_frac < 1.0
            self.main_empty[i] = ~(interval_frac > 0.0)

#            # TEST: Change calculation
#            cont = self.continuum_b[i] + self.continuum_m[i]*wave
#            fint = passband_integral(wave, flux, passband=m, log=log)
#            cint = passband_integral(wave, cont, passband=m, log=log)
#
##            print(interval*(1-fint/cint) - self.index[i])
#            self.index[i] = interval*(1.0-fint/cint if self.units[i] == 'ang' else fint/cint)
#            if err is not None:
#                eint = numpy.sqrt(passband_integral(wave, numpy.square(err), passband=m, log=log))
#                self.index_err[i] = abs(interval*eint/cint)

            # Convert the index to magnitudes
            if self.units[i] == 'mag':
                if not interval > 0.0:
                    self.divbyzero[i] = True
                    self.index[i] = 0.0
                    self.index_err[i] = 0.0
                else:
                    # Worthey et al. 1994, eqn. 3: The passband interval
                    # cancels out of the error propagation.  The error
                    # calculation is done first so as to not replace the
                    # linear calculation of the index.
                    self.index_err[i] = numpy.absolute(2.5 * self.index_err[i] / self.index[i] 
                                                        / numpy.log(10.0))
                    self.index[i] = -2.5 * numpy.ma.log10(self.index[i] / interval)
#        print('')


class BandheadIndices:
    """
    Measure a set of bandhead indices.
    """
    def __init__(self, wave, flux, bluebands, redbands, err=None, log=True, order=None):

        # Check the shape of the input spectrum
        if len(flux.shape) != 1:
            raise ValueError('Input flux must be a single vector!')
        if wave.shape != flux.shape:
            raise ValueError('Input flux and wavelength vectors must have the same shape.')
        
        if len(bluebands.shape) != 2:
            raise ValueError('Band definitions must be two-dimensional: Nindx x 2.')
        if bluebands.shape != redbands.shape:
            raise ValueError('Input bands must have identical shape.')

        self.nindx = bluebands.shape[0]
        self.units = None
        # Check input order
        if order is None:
            warnings.warn('Input order not specified.  Assuming r_b.')
            self.order = numpy.full(self.nindx, 'r_b', dtype=object)
        else:
            self.order = order

        # Get the two pseudocontinua and the flux-weighted band centers
        self.blue_center, self.blue_continuum, self.blue_continuum_err, self.blue_incomplete, \
            self.blue_empty = pseudocontinuum(wave, flux, passband=bluebands, err=err, log=log)

        self.red_center, self.red_continuum, self.red_continuum_err, self.red_incomplete, \
            self.red_empty = pseudocontinuum(wave, flux, passband=redbands, err=err, log=log)
        
        # Initialize the arrays to hold the numerator and denominator
        blue_n = order == 'b_r'

        n = numpy.ma.zeros(self.nindx, dtype=numpy.float)
        d = numpy.ma.zeros(self.nindx, dtype=numpy.float)
        n[blue_n] = self.blue_continuum[blue_n]
        d[blue_n] = self.red_continuum[blue_n]
        n[~blue_n] = self.red_continuum[~blue_n]
        d[~blue_n] = self.blue_continuum[~blue_n]

        nerr = numpy.ma.zeros(self.nindx)
        derr = numpy.ma.zeros(self.nindx)
        if self.blue_continuum_err is not None:
            nerr[blue_n] = self.blue_continuum_err[blue_n]
            derr[~blue_n] = self.blue_continuum_err[~blue_n]
        if self.red_continuum_err is not None:
            nerr[~blue_n] = self.red_continuum_err[~blue_n]
            derr[blue_n] = self.red_continuum_err[blue_n]

        # Determine which indices have both a valid index and index
        # error calculation
        self.main_incomplete = None
        self.main_empty = None
        self.divbyzero = ~((numpy.absolute(d) > 0) & (numpy.absolute(n) > 0))

        # Calculate the indices and their nominal errors
        self.index = numpy.zeros(self.nindx, dtype=numpy.float)
        self.index[~self.divbyzero] = n[~self.divbyzero]/d[~self.divbyzero]
        self.index_err = numpy.zeros(self.nindx, dtype=numpy.float)
        self.index_err[~self.divbyzero] = numpy.sqrt(
                            numpy.square(nerr[~self.divbyzero]*self.index[~self.divbyzero]
                                            / n[~self.divbyzero])
                          + numpy.square(derr[~self.divbyzero]*self.index[~self.divbyzero]
                                            / d[~self.divbyzero]) )


class SpectralIndices:
    r"""

    Class that holds the spectral-index measurements.

    """
    def __init__(self, database_key, binned_spectra, redshift=None, stellar_continuum=None,
                 emission_line_model=None, database_list=None, artifact_list=None,
                 absorption_index_list=None, bandhead_index_list=None, dapsrc=None, dapver=None,
                 analysis_path=None, directory_path=None, output_file=None, hardcopy=True,
                 tpl_symlink_dir=None, clobber=False, checksum=False, loggers=None, quiet=False):

        self.loggers = None
        self.quiet = False

        # Define the database properties
        self.database = None
        self.artdb = None
        self.pixelmask = None
        self.absdb = None
        self.bhddb = None
        self._define_databases(database_key, database_list=database_list,
                               artifact_list=artifact_list,
                               absorption_index_list=absorption_index_list,
                               bandhead_index_list=bandhead_index_list)

        self.nabs, self.nbhd = self.count_indices(self.absdb, self.bhddb)
        self.nindx = self.nabs + self.nbhd
        self.correct_indices = self.database['compute_corrections']

        self.binned_spectra = None
        self.redshift = None
        self.stellar_continuum = None
        self.emission_line_model = None

        # Define the output directory and file
        self.directory_path = None      # Set in _set_paths
        self.output_file = None
        self.hardcopy = None

        # Initialize the objects used in the assessments
        self.bitmask = SpectralIndicesBitMask(dapsrc=dapsrc)

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
        self.measure(binned_spectra, redshift=redshift, stellar_continuum=stellar_continuum,
                     emission_line_model=emission_line_model, dapsrc=dapsrc, dapver=dapver,
                     analysis_path=analysis_path, directory_path=directory_path,
                     output_file=output_file, hardcopy=hardcopy, tpl_symlink_dir=tpl_symlink_dir,
                     clobber=clobber, loggers=loggers, quiet=quiet)


#    def __del__(self):
#        """
#        Deconstruct the data object by ensuring that the fits file is
#        properly closed.
#        """
#        if self.hdu is None:
#            return
#        self.hdu.close()
#        self.hdu = None


    def __getitem__(self, key):
        return self.hdu[key]


    def _define_databases(self, database_key, database_list=None, artifact_list=None,
                          absorption_index_list=None, bandhead_index_list=None, dapsrc=None):
        r"""

        Select the database of indices

        """
        # Grab the specific database
        self.database = select_proc_method(database_key, SpectralIndicesDef,
                                           method_list=database_list,
                                           available_func=available_spectral_index_databases,
                                           dapsrc=dapsrc)

        # Instantiate the artifact, absorption-index, and bandhead-index
        # databases
        self.artdb = None if self.database['artifacts'] is None else \
                ArtifactDB(self.database['artifacts'], artdb_list=artifact_list, dapsrc=dapsrc)
        # TODO: Generalize the name of this object
        self.pixelmask = SpectralPixelMask(artdb=self.artdb)

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
        :func:`mangadap.config.defaults.default_dap_common_path` and
        :func:`mangadap.config.defaults.default_dap_file_name`.

        Args:
            directory_path (str): The exact path to the DAP
                spectral-index file.  See :attr:`directory_path`.
            dapver (str): DAP version.
            analysis_path (str): The path to the top-level directory
                containing the DAP output files for a given DRP and DAP
                version.
            output_file (str): The name of the file with spectral-index
                moment measurements.  See :func:`measure`.
        """
        # Set the output directory path
        method = default_dap_method(binned_spectra=self.binned_spectra,
                                    stellar_continuum=self.stellar_continuum)
        self.analysis_path = default_analysis_path if analysis_path is None else str(analysis_path)
        self.directory_path = default_dap_method_path(method, plate=self.binned_spectra.drpf.plate,
                                                      ifudesign=self.binned_spectra.drpf.ifudesign,
                                                      ref=True,
                                                      drpver=self.binned_spectra.drpf.drpver,
                                                      dapver=dapver,
                                                      analysis_path=self.analysis_path) \
                                        if directory_path is None else str(directory_path)

        # Set the output file
        ref_method = '{0}-{1}'.format(self.binned_spectra.rdxqa.method['key'],
                                      self.binned_spectra.method['key'])
        if self.stellar_continuum is not None:
            ref_method = '{0}-{1}'.format(ref_method, self.stellar_continuum.method['key'])
        if self.emission_line_model is not None:
            ref_method = '{0}-{1}'.format(ref_method, self.emission_line_model.method['key'])
        ref_method = '{0}-{1}'.format(ref_method, self.database['key'])

        self.output_file = default_dap_file_name(self.binned_spectra.drpf.plate,
                                                 self.binned_spectra.drpf.ifudesign, ref_method) \
                                        if output_file is None else str(output_file)


    def _initialize_primary_header(self, hdr=None):
        """
        Initialize the header of :attr:`hdu`.

        Returns:
            astropy.io.fits.Header : Edited header object.

        """
        # Copy the from the DRP and clean it
        if hdr is None:
            hdr = self.binned_spectra.drpf.hdu['PRIMARY'].header.copy()
            hdr = DAPFitsUtil.clean_dap_primary_header(hdr)
        
        hdr['AUTHOR'] = 'Kyle B. Westfall <westfall@ucolick.org>'
        hdr['SIKEY'] = (self.database['key'], 'Spectral-index database keyword')
        hdr['SIMINSN'] = (self.database['minimum_snr'], 'Minimum S/N of spectrum to include')
        hdr['SIFWHM'] = (self.database['fwhm'], 'FWHM of index system resolution (ang)')
        hdr['ARTDB'] = (self.database['artifacts'], 'Artifact database keyword')
        hdr['ABSDB'] = (self.database['absindex'], 'Absorption-index database keyword')
        hdr['BHDDB'] = (self.database['bandhead'], 'Bandhead-index database keyword')
        if self.stellar_continuum is not None:
            hdr['SCKEY'] = (self.stellar_continuum.method['key'], 'Stellar-continuum model keyword')
        if self.emission_line_model is not None:
            hdr['ELFKEY'] = (self.emission_line_model.method['key'],
                                'Emission-line modeling method keyword')
        hdr['NBINS'] = (self.nbins, 'Number of unique spatial bins')
#        if len(self.missing_bins) > 0:
#            hdr['EMPTYBIN'] = (str(self.missing_bins), 'List of bins with no data')
#        hdr['SICORR'] = (self.correct_indices, 'Indices corrected for velocity dispersion')
        # Anything else?
        # Additional database details?
        return hdr


#    def _initialize_mask(self, good_snr):
#        """
#
#        Initialize the mask be setting the DIDNOTUSE, FORESTAR, and LOW_SNR masks
#
#        """
#        # Initialize to all zeros
#        mask = numpy.zeros(self.shape, dtype=self.bitmask.minimum_dtype())
#
#        # Turn on the flag stating that the pixel wasn't used
##        print('Mask size:', self.binned_spectra['MASK'].data.size)
#        indx = self.binned_spectra.bitmask.flagged(self.binned_spectra['MASK'].data,
#                                                   flag=self.binned_spectra.do_not_fit_flags())
##        print('Masked as DIDNOTUSE:', numpy.sum(indx))
#        mask[indx] = self.bitmask.turn_on(mask[indx], 'DIDNOTUSE')
#
#        # Turn on the flag stating that the pixel has a foreground star
#        indx = self.binned_spectra.bitmask.flagged(self.binned_spectra['MASK'].data,
#                                                   flag='FORESTAR')
##        print('Masked as FORESTAR: ', numpy.sum(indx))
#        mask[indx] = self.bitmask.turn_on(mask[indx], 'FORESTAR')
#
#        # Turn on the flag stating that the S/N in the spectrum was
#        # below the requested limit
#        indx = numpy.array([numpy.invert(good_snr).reshape(self.spatial_shape).T]*self.nwave).T
#        mask[indx] = self.bitmask.turn_on(mask[indx], flag='LOW_SNR')
#
#        return mask


    def _index_database_dtype(self, name_len):
        r"""
        Construct the record array data type for the output fits
        extension.
        """
        return [ ('TYPE','<U10'),
                 ('ID',numpy.int),
                 ('NAME','<U{0:d}'.format(name_len)),
                 ('PASSBAND', numpy.float, (2,)),
                 ('BLUEBAND', numpy.float, (2,)),
                 ('REDBAND', numpy.float, (2,)),
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
        if self.absdb is not None:
            for n in self.absdb['name']:
                if name_len < len(n):
                    name_len = len(n)
        if self.bhddb is not None:
            for n in self.bhddb['name']:
                if name_len < len(n):
                    name_len = len(n)

        # Instatiate the table data that will be saved defining the set
        # of indices measured
        passband_database = init_record_array(self.nindx, self._index_database_dtype(name_len))

        t = 0 if self.absdb is None else self.absdb.nsets
        if self.absdb is not None:
            passband_database['TYPE'][:t] = 'absorption'
            hk = [ 'ID', 'NAME', 'PASSBAND', 'BLUEBAND', 'REDBAND', 'UNIT', 'COMPONENT' ]
            ak = [ 'index', 'name', 'primary', 'blueside', 'redside', 'units', 'component' ]
            for _hk, _ak in zip(hk,ak):
                passband_database[_hk][:t] = self.absdb[_ak]
        if self.bhddb is not None:
            passband_database['TYPE'][t:] = 'bandhead'
            hk = [ 'ID', 'NAME', 'BLUEBAND', 'REDBAND', 'INTEGRAND', 'ORDER' ]
            ak = [ 'index', 'name', 'blueside', 'redside', 'integrand', 'order' ]
            for _hk, _ak in zip(hk,ak):
                passband_database[_hk][t:] = self.bhddb[_ak]
        return passband_database


#    def _assign_spectral_arrays(self):
#        self.spectral_arrays = [ 'FLUX', 'IVAR', 'MASK' ]


    def _assign_image_arrays(self):
        """
        Set :attr:`image_arrays`, which contains the list of extensions
        in :attr:`hdu` that are on-sky image data.
        """
        self.image_arrays = [ 'BINID' ]


    def _get_missing_bins(self):
        good_snr = self.binned_spectra.above_snr_limit(self.database['minimum_snr'])
        return numpy.sort(self.binned_spectra['BINS'].data['BINID'][~good_snr].tolist()
                                + self.binned_spectra.missing_bins) 


#    def _check_snr(self):
#        # binned_spectra['BINS'].data['SNR'] has length
#        # binned_spectra.nbins
#        return self.binned_spectra['BINS'].data['SNR'] > self.database['minimum_snr']

    
#    def _missing_flags(self):
#        return numpy.array([ b in self.missing_bins for b in numpy.arange(self.nbins)])
#

#    def _bins_to_measure(self):
#        """
#        Determine which bins to use for the spectral-index measurements.
#        They must not be designated as "missing" and they must have
#        sufficient S/N.
#        """
#        return (self._check_snr()) & ~(self._missing_flags())
##        good_bins = (self._check_snr()) & ~(self._missing_flags())
##        good_bins[2] = False
##        return good_bins


    def _assign_redshifts(self, redshift):
        """
        Set the redshift to apply to each spectrum.
        """
        self.redshift = numpy.zeros(self.binned_spectra.nbins, dtype=numpy.float)
        if redshift is None:
            return

        _redshift = numpy.atleast_1d(redshift)
        if len(_redshift) not in [ 1, self.binned_spectra.nbins ]:
            raise ValueError('Provided redshift must be either a single value or match the ' \
                             'number of binned spectra.')
        self.redshift = numpy.full(self.binned_spectra.nbins, redshift, dtype=float) \
                            if len(_redshift) == 1 else _redshift.copy()


    @staticmethod
    def spectra_for_index_measurements(binned_spectra, pixelmask=None, select=None,
                                       resolution_fwhm=None, emission_line_model=None):
        """
        Compile the set of spectra for the spectral-index measurements.
        If the input fwhm is > 0, this function will match the spectral
        resolution of the data to the spectral-index system, based on
        the provided FWHM.

        """
        # Get the main data arrays
        wave = binned_spectra['WAVE'].data
        flux = binned_spectra.copy_to_masked_array(flag=binned_spectra.do_not_fit_flags())
        ivar = binned_spectra.copy_to_masked_array(ext='IVAR',
                                                   flag=binned_spectra.do_not_fit_flags())

        # Mask any pixels in the pixel mask
        if pixelmask is not None:
            indx = pixelmask.boolean(wave, nspec=binned_spectra.nbins)
            flux[indx] = numpy.ma.masked
            ivar[indx] = numpy.ma.masked

        # TODO: Need to understand if baseline effects should be taken
        # out as well...
        if emission_line_model is not None:
            eml_model = emission_line_model.fill_to_match(binned_spectra)
            no_eml = numpy.invert(numpy.ma.getmaskarray(flux)) & numpy.ma.getmaskarray(eml_model)
            flux -= eml_model
            flux.mask[no_eml] = False

        # Make sure ivar mask is identical to flux mask
        ivar.mask = flux.mask.copy()

        # Set the selected spectra
        _select = numpy.ones(binned_spectra.nbins, dtype=bool) if select is None else select

#        pyplot.step(wave, flux[0,:], where='mid', linestyle='-', color='k', lw=1)
#        _flux = binned_spectra.copy_to_masked_array(flag=binned_spectra.do_not_fit_flags())
#        pyplot.step(wave, _flux[0,:], where='mid', linestyle='-', color='b', lw=1)
#        pyplot.show()

        return flux[select,:], ivar[select,:] if resolution_fwhm < 0 \
                else adjust_spectral_resolution(wave, flux[select,:], ivar[select,:],
                                                binned_spectra['SPECRES'].data.copy()[select,:],
                                                resolution_fwhm)

    @staticmethod
    def adjust_spectral_resolution(wave, flux, ivar, sres, resolution_fwhm):
        """
        flux and ivar are expected to be masked arrays
        """
        # Revert flux and ivar to unmasked arrays
        mask = numpy.ma.getmaskarray(flux)
        flux = numpy.asarray(flux)
        ivar = numpy.asarray(ivar)
        nspec = flux.shape[0]

        # For masked pixels, interpolate across them using simple linear
        # interpolation.
        # TODO: This approach needs to be vetted.
        for i in range(nspec):
            mp = mask[i,:]              # Masked pixels
            up = numpy.invert(mp)       # Unmasked pixels
            if numpy.sum(mp) == 0:      # No masked pixels so continue
                continue
            if numpy.sum(up) == 0:      # No unmasked pixels!
                flux[i,:] = 0.0
                continue
            # Linearly interpolate both the flux and inverse variance
            flux[i,mp] = numpy.interp(wave[mp], wave[up], flux[i,up])
            ivar[i,mp] = numpy.interp(wave[mp], wave[up], ivar[i,up])

#        pyplot.step(wave, flux[0,:], where='mid', linestyle='-', color='r', lw=0.5)
#        pyplot.show()

        # Use :func:`mangadap.instrument.match_spectral_resolution` to
        # match the spectral resolution of the binned spectra to the
        # spectral-index system
        new_sres = wave/resolution_fwhm
        
        _flux, _sres, sigoff, _mask, _ivar \
                = match_spectral_resolution(wave, flux, sres, wave, new_sres, ivar=ivar,
                                            log10=True, new_log10=True)

        # FOR DEBUGGING
#        warnings.warn('NOT MATCHING SPECTRAL RESOLUTION!')
#        new_flux = flux.copy()
#        new_mask = mask.copy()

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
        if numpy.sum(_mask) > 0:
            mask[_mask > 0] = True
        return numpy.ma.MaskedArray(flux,mask=mask), numpy.ma.MaskedArray(ivar,mask=mask)


    @staticmethod
    def unit_selection(absdb, bhddb):
        nabs, nbhd = SpectralIndices.count_indices(absdb, bhddb)
        
        # Flag the indices as either having magnitude or angstrom units
        angu = numpy.zeros(nabs + nbhd, dtype=bool)
        if nabs > 0:
            angu[:nabs] = absdb['units'] == 'ang'
        magu = numpy.invert(angu)

        # Bandhead indices are unitless, but the correction follows the
        # same operation as for the absorption-line indices with
        # angstrom units
        magu[nabs:] = False
        angu[nabs:] = True

        return angu, magu


    def _resolution_matched_templates(self, dapsrc=None, dapver=None, analysis_path=None,
                                      tpl_symlink_dir=None):
        """
        Get a version of the template library that has had its
        resolution matched to that of the spectral-index database.
        """
        velocity_offset=self.stellar_continuum.method['fitpar']['template_library'].velocity_offset
        return self.stellar_continuum.get_template_library(dapsrc=dapsrc, dapver=dapver,
                                                           analysis_path=analysis_path,
                                                           tpl_symlink_dir=tpl_symlink_dir,
                                                           velocity_offset=velocity_offset,
                                                           resolution_fwhm=self.database['fwhm'])
    

    @staticmethod
    def calculate_dispersion_corrections(absdb, bhddb, wave, flux, continuum, continuum_dcnvlv,
                                         redshift=None, bitmask=None):
        """
        Calculate the dispersion corrections using the best-fitting
        template models.
        """
#        pyplot.step(wave, flux[0,:], where='mid', linestyle='-', color='k', lw=0.5, zorder=1)
#        pyplot.plot(wave, continuum[0,:], linestyle='-', color='g', lw=1.0, zorder=3, alpha=0.5)
#        pyplot.plot(wave, continuum_dcnvlv[0,:], linestyle='-', color='b', lw=1.0, zorder=3,
#                    alpha=0.5)
#        pyplot.show()

        # Make sure continuum includes the flux masks
        nspec = flux.shape[0]
        _flux = numpy.ma.MaskedArray(flux)
        _continuum = numpy.ma.MaskedArray(continuum)
        _continuum_dcnvlv = numpy.ma.MaskedArray(continuum_dcnvlv)
        _continuum[_flux.mask] = numpy.ma.masked
        _continuum_dcnvlv[_flux.mask] = numpy.ma.masked

        # Measure the indices for both sets of spectra
        indx = SpectralIndices.measure_indices(absdb, bhddb, wave, _continuum, redshift=redshift,
                                               bitmask=bitmask)
        dcnvlv_indx = SpectralIndices.measure_indices(absdb, bhddb, wave, _continuum_dcnvlv,
                                                      redshift=redshift, bitmask=bitmask)

        # Do not apply the correction if any of the bands were empty or
        # would have had to divide by zero
        bad_indx = bitmask.flagged(indx['MASK'], flag=[ 'MAIN_EMPTY', 'BLUE_EMPTY', 'RED_EMPTY',
                                                        'DIVBYZERO' ])
        bad_indx |= bitmask.flagged(dcnvlv_indx['MASK'], flag=[ 'MAIN_EMPTY', 'BLUE_EMPTY',
                                                                'RED_EMPTY', 'DIVBYZERO' ])

        print('Number of indices: ', bad_indx.size)
        print('Bad indices: ', numpy.sum(bad_indx))

        # Determine which indices have good measurements
        angu, magu = SpectralIndices.unit_selection(absdb, bhddb)
        good_ang = ~bad_indx & (numpy.array([angu]*nspec)) & (numpy.absolute(indx['INDX']) > 0)
        good_mag = ~bad_indx & (numpy.array([magu]*nspec))

        print('Good angstrom indices: ', numpy.sum(good_ang))
        print('Good magnitude indices: ', numpy.sum(good_mag))

        # Determine and return the corrections to apply
        corrections = numpy.zeros(indx['INDX'].shape, dtype=numpy.float)
        corrections[good_ang] = dcnvlv_indx['INDX'][good_ang] / indx['INDX'][good_ang]
        corrections[good_mag] = dcnvlv_indx['INDX'][good_mag] - indx['INDX'][good_mag]
        return corrections, good_ang, good_mag


    @staticmethod
    def count_indices(absdb, bhddb):
        return 0 if absdb is None else absdb.nsets, 0 if bhddb is None else bhddb.nsets


    @staticmethod
    def output_dtype(nindx, bitmask=None):
        r"""
        Construct the record array data type for the output fits
        extension.
        """
        return [ ('BINID',numpy.int),
                 ('BINID_INDEX',numpy.int),
                 ('REDSHIFT', numpy.float),
                 ('MASK', numpy.bool if bitmask is None else bitmask.minimum_dtype(), (nindx,)),
                 ('BCEN', numpy.float, (nindx,)), 
                 ('BCONT', numpy.float, (nindx,)), 
                 ('BCONTERR', numpy.float, (nindx,)),
                 ('RCEN', numpy.float, (nindx,)), 
                 ('RCONT', numpy.float, (nindx,)), 
                 ('RCONTERR', numpy.float, (nindx,)), 
                 ('INDX_DISPCORR', numpy.float, (nindx,)), 
                 ('INDX', numpy.float, (nindx,)), 
                 ('INDXERR', numpy.float, (nindx,))
               ]


    @staticmethod
    def check_and_prep_input(wave, flux, ivar=None, mask=None, redshift=None, bitmask=None):

        # Check the bitmask if provided
        if bitmask is not None and not isinstance(bitmask, BitMask):
            raise TypeError('Input bitmask must have type BitMask.')

        # Check the input wavelength and flux shapes
        if len(wave.shape) != 1:
            raise ValueError('Input wavelengths must be a single vector.')
        if len(wave) != flux.shape[1]:
            raise ValueError('Wavelength vector does not match shape of the flux array.')

        # Check the mask shape
        if mask is not None and mask.shape != flux.shape:
            raise ValueError('Input mask must have the same shape as the flux array.')

        # Check the input redshifts
        nspec = flux.shape[0]
        _redshift = numpy.zeros(nspec, dtype=numpy.float) if redshift is None else redshift
        if len(_redshift) != nspec:
            raise ValueError('Must provide one redshift per input spectrum (flux.shape[0]).')

        # Convert the input arrays to masked arrays if they aren't
        # already, and compare the array shapes
        _flux = flux if isinstance(flux, numpy.ma.MaskedArray) else \
                    numpy.ma.MaskedArray(flux, mask=(None if mask is None else mask > 0))

        if ivar is None:
            noise = None
        else:
            if ivar.shape != flux.shape:
                raise ValueError('Input ivar array must be the same shape as the flux array.')
            _ivar = ivar if isinstance(ivar, numpy.ma.MaskedArray) else \
                        numpy.ma.MaskedArray(ivar, mask=(None if mask is None else mask > 0))
            noise = numpy.ma.sqrt(1.0 /_ivar)

        return _flux, noise, _redshift


#    @staticmethod
#    def sideband_pseudocontinua(wave, spec, sidebands, noise=None, log=True):
#        """Get the side-band integrals in a single spectrum."""
#
#        # Calculate the pseudo-continua in the sidebands
#        nbands = sidebands.shape[0]
#        
#        pseudocontinuum, pseudocontinuum_error \
#                    = passband_integrated_mean(wave, spec, passband=sidebands, err=noise, log=log)
#        flux_weighted_center, _fwc_err \
#                    = passband_weighted_mean(wave, spec, wave, passband=sidebands, log=log)
#
#        # Calculate the fraction of the band that is covered by unmasked
#        # pixels
#        interval_frac = passband_integrated_width(wave, spec, passband=sidebands, log=log) \
#                                / numpy.diff(sidebands, axis=1).ravel()
#
##        if len(spec.shape) > 1:
##            pyplot.step(wave, spec[0], where='mid', color='k', lw=0.5, zorder=1)
##            pyplot.step(wave, spec[1], where='mid', color='g', lw=0.5, zorder=2)
##        else:
##            pyplot.step(wave, spec, where='mid', color='k', lw=0.5, zorder=1)
##        pyplot.scatter(flux_weighted_center-1.0, pseudocontinuum, marker='.', s=50,
##                       color='b', lw=0, zorder=3)
##        pyplot.show()
#
#        return flux_weighted_center, pseudocontinuum, pseudocontinuum_error, \
#                            interval_frac < 1.0, ~(interval_frac > 0.0)


    @staticmethod
    def set_masks(measurements, blue_incomplete, blue_empty, red_incomplete, red_empty, divbyzero,
                  main_incomplete=None, main_empty=None, bitmask=None):

#        print('blue incomplete: {0}/{1}'.format(numpy.sum(blue_incomplete), blue_incomplete.size))
#        print('blue empty: {0}/{1}'.format(numpy.sum(blue_empty), blue_empty.size))
#        print('red incomplete: {0}/{1}'.format(numpy.sum(red_incomplete), red_incomplete.size))
#        print('red empty: {0}/{1}'.format(numpy.sum(red_empty), red_empty.size))
#        print('divbyzero: {0}/{1}'.format(numpy.sum(divbyzero), divbyzero.size))
#        if main_incomplete is not None:
#            print('main incomplete: {0}/{1}'.format(numpy.sum(main_incomplete),
#                                                    main_incomplete.size))
#        if main_empty is not None:
#            print('main empty: {0}/{1}'.format(numpy.sum(main_empty), main_empty.size))

        measurements['MASK'][blue_incomplete] = True if bitmask is None else \
                bitmask.turn_on(measurements['MASK'][blue_incomplete], 'BLUE_INCOMP')
        measurements['MASK'][blue_empty] = True if bitmask is None else \
                bitmask.turn_on(measurements['MASK'][blue_empty], 'BLUE_EMPTY')
        measurements['MASK'][red_incomplete] = True if bitmask is None else \
                bitmask.turn_on(measurements['MASK'][red_incomplete], 'RED_INCOMP')
        measurements['MASK'][red_empty] = True if bitmask is None else \
                bitmask.turn_on(measurements['MASK'][red_empty], 'RED_EMPTY')
        measurements['MASK'][divbyzero] = True if bitmask is None else \
                bitmask.turn_on(measurements['MASK'][divbyzero], 'DIVBYZERO')
        if main_incomplete is not None:
            measurements['MASK'][main_incomplete] = True if bitmask is None else \
                    bitmask.turn_on(measurements['MASK'][main_incomplete], 'MAIN_INCOMP')
        if main_empty is not None:
            measurements['MASK'][main_empty] = True if bitmask is None else \
                    bitmask.turn_on(measurements['MASK'][main_empty], 'MAIN_EMPTY')
        return measurements


    @staticmethod
    def save_results(results, measurements, good, err=False, bitmask=None):
        if not isinstance(results, (AbsorptionLineIndices, BandheadIndices)):
            raise TypeError('Input must be of type AbsorptionLineIndices or BandheadIndices')

        # Save the data
        measurements['BCEN'][good] = results.blue_center
        measurements['BCONT'][good] = results.blue_continuum
        if err:
            measurements['BCONTERR'][good] = results.blue_continuum_err
        measurements['RCEN'][good] = results.red_center
        measurements['RCONT'][good] = results.red_continuum
        if err:
            measurements['RCONTERR'][good] = results.red_continuum_err
        measurements['INDX'][good] = results.index
        if err:
            measurements['INDXERR'][good] = results.index_err

        # Reshape the flags
        blue_incomplete = numpy.zeros(len(good), dtype=numpy.bool)
        blue_incomplete[good] = results.blue_incomplete
        blue_empty = numpy.zeros(len(good), dtype=numpy.bool)
        blue_empty[good] = results.blue_empty
        red_incomplete = numpy.zeros(len(good), dtype=numpy.bool)
        red_incomplete[good] = results.red_incomplete
        red_empty = numpy.zeros(len(good), dtype=numpy.bool)
        red_empty[good] = results.red_empty
        main_incomplete = None
        if results.main_incomplete is not None:
            main_incomplete = numpy.zeros(len(good), dtype=numpy.bool)
            main_incomplete[good] = results.main_incomplete
        main_empty = None
        if results.main_empty is not None:
            main_empty = numpy.zeros(len(good),dtype=numpy.bool)
            main_empty[good] = results.main_empty
        divbyzero = numpy.zeros(len(good), dtype=numpy.bool)
        divbyzero[good] = results.divbyzero

        # Set the masks
        return SpectralIndices.set_masks(measurements, blue_incomplete, blue_empty, red_incomplete,
                                         red_empty, divbyzero, main_incomplete=main_incomplete,
                                         main_empty=main_empty, bitmask=bitmask)


    @staticmethod
    def measure_indices(absdb, bhddb, wave, flux, ivar=None, mask=None, redshift=None,
                        bitmask=None):
        """
        Measure the spectral indices in a set of spectra.

        If not input as masked arrays, flux is converted to one.

        """
        # Check the input databases
        if absdb is not None and not isinstance(absdb, AbsorptionIndexDB):
            raise TypeError('Input database must have type AbsorptionIndexDB.')
        if bhddb is not None and not isinstance(bhddb, BandheadIndexDB):
            raise TypeError('Input database must have type BandheadIndexDB.')

        # Get the number of indices
        nabs, nbhd = SpectralIndices.count_indices(absdb, bhddb)
        nindx = nabs+nbhd
        if nindx == 0:
            raise ValueError('No indices to measure!')

        _flux, noise, _redshift = SpectralIndices.check_and_prep_input(wave, flux, ivar=ivar,
                                                                       mask=mask, redshift=redshift,
                                                                       bitmask=bitmask)
        nspec = _flux.shape[0]
#        print(_flux.shape)
#        print(None if noise is None else noise.shape)
#        print(wave.shape)
#        print((numpy.array([wave]*nspec)).shape)

        # Create the f_nu spectra
        # The conversion is multiplicative, meaning the calculation of
        # the error calculation can use exactly the same function 
        flux_fnu = flux_to_fnu(numpy.array([wave]*nspec), _flux)
        noise_fnu = None if noise is None else flux_to_fnu(numpy.array([wave]*nspec), noise)

#        nu = 1e-13*astropy.constants.c.to('nm/s').value/wave    # In THz
#        pyplot.plot(wave, _flux[0,:]*wave, linestyle='-', color='r', lw=1, zorder=2)
#        pyplot.plot(wave, flux_fnu[0,:]*nu, linestyle='-', color='0.5', lw=3, zorder=1)
#        pyplot.show()

        # Get the list of good indices of each type
        abs_fnu = numpy.zeros(nabs, dtype=numpy.bool) if absdb is None \
                        else ~absdb.dummy & (absdb['integrand'] == 'fnu')
        good_abs_fnu = numpy.zeros(nindx, dtype=numpy.bool)
        if nabs > 0:
            good_abs_fnu[:nabs][abs_fnu] = True
        abs_flambda = numpy.zeros(nabs, dtype=numpy.bool) if absdb is None \
                        else ~absdb.dummy & (absdb['integrand'] == 'flambda')
        good_abs_flambda = numpy.zeros(nindx, dtype=numpy.bool)
        if nabs > 0:
            good_abs_flambda[:nabs][abs_flambda] = True

        bhd_fnu = numpy.zeros(nbhd, dtype=numpy.bool) if bhddb is None \
                        else ~bhddb.dummy & (bhddb['integrand'] == 'fnu')
        good_bhd_fnu = numpy.zeros(nindx, dtype=numpy.bool)
        if nbhd > 0:
            good_bhd_fnu[nabs:][bhd_fnu] = True
        bhd_flambda = numpy.zeros(nbhd, dtype=numpy.bool) if bhddb is None \
                        else ~bhddb.dummy & (bhddb['integrand'] == 'flambda')
        good_bhd_flambda = numpy.zeros(nindx, dtype=numpy.bool)
        if nbhd > 0:
            good_bhd_flambda[nabs:][bhd_flambda] = True

        # Initialize the output data
        measurements = init_record_array(nspec, SpectralIndices.output_dtype(nindx,bitmask=bitmask))

        # Perform the measurements on each spectrum
        for i in range(nspec):

            print('Measuring spectral indices in spectrum: {0}/{1}'.format(i+1,nspec), end='\r')

            # -----------------------------------
            # Measure the absorption-line indices
            if nabs > 0:

                # Shift the bands
                _bluebands = absdb['blueside']*(1.0+_redshift[i])
                _redbands = absdb['redside']*(1.0+_redshift[i])
                _mainbands = absdb['primary']*(1.0+_redshift[i])

                # Integrate over F_nu
                if numpy.sum(good_abs_fnu) > 0:
                    # Make the measurements ...
                    _noise = None if noise_fnu is None else noise_fnu[i,:]
                    results = AbsorptionLineIndices(wave, flux_fnu[i,:], _bluebands[abs_fnu],
                                                    _redbands[abs_fnu], _mainbands[abs_fnu],
                                                    err=_noise, units=absdb['units'][abs_fnu])
                    # ... and save them
                    measurements[i] = SpectralIndices.save_results(results, measurements[i],
                                                                   good_abs_fnu,
                                                                   err=_noise is not None,
                                                                   bitmask=bitmask)

                # Integrate over F_lambda
                if numpy.sum(good_abs_flambda) > 0:
                    # Make the measurements ...
                    _noise = None if noise is None else noise[i,:]
                    results = AbsorptionLineIndices(wave, _flux[i,:], _bluebands[abs_flambda],
                                                    _redbands[abs_flambda], _mainbands[abs_flambda],
                                                    err=_noise, units=absdb['units'][abs_flambda])
                    # ... and save them
                    measurements[i] = SpectralIndices.save_results(results, measurements[i],
                                                                   good_abs_flambda,
                                                                   err=_noise is not None,
                                                                   bitmask=bitmask)

            # -----------------------------------
            # Measure the bandhead indices
            if nbhd > 0:

                # Shift the bands
                _bluebands = bhddb['blueside']*(1.0+_redshift[i])
                _redbands = bhddb['redside']*(1.0+_redshift[i])

                # Integrate over F_nu
                if numpy.sum(good_bhd_fnu) > 0:
                    # Make the measurements ...
                    _noise = None if noise_fnu is None else noise_fnu[i,:]
                    results = BandheadIndices(wave, flux_fnu[i,:], _bluebands[bhd_fnu],
                                              _redbands[bhd_fnu], err=_noise,
                                              order=bhddb['order'][bhd_fnu])
                    # ... and save them
                    measurements[i] = SpectralIndices.save_results(results, measurements[i],
                                                                   good_bhd_fnu,
                                                                   err=_noise is not None,
                                                                   bitmask=bitmask)

                # Integrate over F_lambda
                if numpy.sum(good_bhd_flambda) > 0:
                    # Make the measurements ...
                    _noise = None if noise is None else noise[i,:]
                    results = BandheadIndices(wave, _flux[i,:], _bluebands[bhd_flambda],
                                              _redbands[bhd_flambda], err=_noise,
                                              order=bhddb['order'][bhd_flambda])
                    # ... and save them
                    measurements[i] = SpectralIndices.save_results(results, measurements[i],
                                                                   good_bhd_flambda,
                                                                   err=_noise is not None,
                                                                   bitmask=bitmask)

#            pyplot.scatter(_interval, interval/_interval, marker='.', color='k', s=50)
#            pyplot.show()

#            center = numpy.mean(_sidebands, axis=1).reshape(-1,self.nindx)
#            pyplot.step(wave, flux[i,:], where='mid', color='k', lw=0.5, linestyle='-')
#            pyplot.scatter(center[0,:], measurements['BCONT'][i,:], marker='.',
#                           s=100, color='b', lw=0)
#            pyplot.scatter(center[1,:], measurements['RCONT'][i,:], marker='.',
#                           s=100, color='r', lw=0)
#            pyplot.show()

#        x = numpy.append(center, numpy.array([[-100]*self.absdb.nsets]).T,axis=1).ravel()
#        _x = numpy.ma.MaskedArray(x, mask=x==-100)
#        print(_x)
#        y = numpy.append(pseudocontinuum, numpy.array([[-100]*self.absdb.nsets]).T,axis=1).ravel()
#        _y = numpy.ma.MaskedArray(y, mask=y==-100)
#        print(type(_y))
#        print(numpy.sum(_y.mask))
#        pyplot.plot(_x, _y, color='g', lw=2, linestyle='-')
#        pyplot.show()

#        print('Total masked: {0}/{1}'.format(numpy.sum(measurements['MASK'] > 0),
#                                             measurements['MASK'].size))

        # Return the data
        return measurements


    def file_name(self):
        """Return the name of the output file."""
        return self.output_file


    def file_path(self):
        """Return the full path to the output file."""
        if self.directory_path is None or self.output_file is None:
            return None
        return os.path.join(self.directory_path, self.output_file)


    def measure(self, binned_spectra, redshift=None, stellar_continuum=None,
                emission_line_model=None, dapsrc=None, dapver=None, analysis_path=None,
                directory_path=None, output_file=None, hardcopy=True, tpl_symlink_dir=None,
                clobber=False, loggers=None, quiet=False):
        """
        Measure the spectral indices using the binned spectra.
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

        # StellarContinuumModel object only used when calculating dispersion corrections.
        if stellar_continuum is not None:
            if not isinstance(stellar_continuum, StellarContinuumModel):
                raise TypeError('Provided stellar continuum must have StellarContinuumModel type!')
            if stellar_continuum.hdu is None:
                raise ValueError('Provided StellarContinuumModel is undefined!')
            self.stellar_continuum = stellar_continuum
        if self.stellar_continuum is not None and self.database['compute_corrections']:
            warnings.warn('Cannot compute dispersion corrections without StellarContinuumModel.')
            self.correct_indices = False

        # EmissionLineModel object used to subtract emission lines
        if emission_line_model is not None:
            if not isinstance(emission_line_model, EmissionLineModel):
                raise TypeError('Provided emission line models must be of type EmissionLineModel.')
            if emission_line_model.hdu is None:
                raise ValueError('Provided EmissionLineModel is undefined!')
            self.emission_line_model = emission_line_model

        self.spatial_shape =self.binned_spectra.spatial_shape
        self.nspec = self.binned_spectra.nspec
        self.spatial_index = self.binned_spectra.spatial_index.copy()
        
        # Get the redshifts to apply
        self._assign_redshifts(redshift)

        #---------------------------------------------------------------
        # Get the good spectra
        good_snr = self.binned_spectra.above_snr_limit(self.database['minimum_snr'])

        # Report
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, '-'*50)
            log_output(self.loggers, 1, logging.INFO, 'SPECTRAL-INDEX MEASUREMENTS:')
            log_output(self.loggers, 1, logging.INFO, '-'*50)
            log_output(self.loggers, 1, logging.INFO, 'Number of binned spectra: {0}'.format(
                                                            self.binned_spectra.nbins))
            if len(self.binned_spectra.missing_bins) > 0:
                log_output(self.loggers, 1, logging.INFO, 'Missing bins: {0}'.format(
                                                            len(self.binned_spectra.missing_bins)))
            log_output(self.loggers, 1, logging.INFO, 'With good S/N and to measure: {0}'.format(
                                                            numpy.sum(good_snr)))
            
        # Make sure there are good spectra
        if numpy.sum(good_snr) == 0:
            raise ValueError('No good spectra for measurements!')

        #---------------------------------------------------------------
        # (Re)Set the output paths
        self._set_paths(directory_path, dapver, analysis_path, output_file)

        #---------------------------------------------------------------
        # Check that the file path is defined
        ofile = self.file_path()
        if ofile is None:
            raise ValueError('File path for output file is undefined!')

        # Report
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, 'Output path: {0}'.format(
                                                                            self.directory_path))
            log_output(self.loggers, 1, logging.INFO, 'Output file: {0}'.format(
                                                                            self.output_file))
        
        #---------------------------------------------------------------
        # If the file already exists, and not clobbering, just read the
        # file
        if os.path.isfile(ofile) and not clobber:
            self.hardcopy = True
            if not self.quiet:
                log_output(self.loggers, 1, logging.INFO, 'Using existing file')
            self.read(checksum=self.checksum)
            if not self.quiet:
                log_output(self.loggers, 1, logging.INFO, '-'*50)
            return

#        pyplot.scatter(numpy.arange(self.nbins), self.redshift, marker='.', s=50, color='k', lw=0)
#        pyplot.show()

        #---------------------------------------------------------------
        # Set the number of bins measured and missing bins
        self.nbins = numpy.sum(good_snr)
        self.missing_bins = self._get_missing_bins()

        #---------------------------------------------------------------
        # Get the spectra to use for the measurements
        flux, ivar = SpectralIndices.spectra_for_index_measurements(self.binned_spectra,
                                                    pixelmask=self.pixelmask, select=good_snr,
                                                    resolution_fwhm=self.database['fwhm'],
                                                    emission_line_model=self.emission_line_model)

        #---------------------------------------------------------------
        # Perform the measurements on the galaxy spectra
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO,
                       'Measuring spectral indices in observed spectra...')
        measurements = self.measure_indices(self.absdb, self.bhddb,
                                            self.binned_spectra['WAVE'].data, flux, ivar=ivar,
                                            redshift=self.redshift[good_snr], bitmask=self.bitmask)

        #---------------------------------------------------------------
        # Determine the velocity dispersion corrections, if requested.

        # Get the template spectra to use
        replacement_templates = None if self.database['fwhm'] < 0 \
                    else self._resolution_matched_templates(dapsrc=dapsrc, dapver=dapver,
                                                            analysis_path=analysis_path,
                                                            tpl_symlink_dir=tpl_symlink_dir)

        # Get the continuum with and without the LOSVD convolution
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, 'Constructing models with LOSVD')
        continuum = self.stellar_continuum.fill_to_match(self.binned_spectra,
                                                        replacement_templates=replacement_templates)
        
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, 'Constructing models without LOSVD')
        continuum_dcnvlv = self.stellar_continuum.fill_to_match(self.binned_spectra,
                                                        replacement_templates=replacement_templates,
                                                        redshift_only=True)
        

        # Get the corrections by performing the measurements on the
        # best-fitting continuum models, with and without the velocity
        # dispersion broadening
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO,
                       'Calculating dispersion corrections using stellar continuum model...')
        measurements['INDX_DISPCORR'], good_ang, good_mag \
                = SpectralIndices.calculate_dispersion_corrections(self.absdb, self.bhddb,
                                                        self.binned_spectra['WAVE'].data,
                                                        flux, continuum[good_snr,:],
                                                        continuum_dcnvlv[good_snr,:],
                                                        redshift=self.redshift[good_snr],
                                                        bitmask=self.bitmask)

        # Add in the ancillary information and initialize the bins
        measurements['BINID'] = self.binned_spectra['BINS'].data['BINID'][good_snr]
        measurements['BINID_INDEX'] = numpy.arange(self.binned_spectra.nbins)[good_snr]

        # Flag bad corrections
        bad_correction = ~good_ang & ~good_mag
        measurements['MASK'][bad_correction] = self.bitmask.turn_on(
                                                        measurements['MASK'][bad_correction],
                                                        'NO_DISPERSION_CORRECTION')

#        # Apply the corrections
#        if self.correct_indices:
#            measurements['INDX'][good_ang] *= measurements['INDX_DISPCORR'][good_ang]
#            measurements['INDXERR'][good_ang] \
#                    *= numpy.absolute(measurements['INDX_DISPCORR'][good_ang])
#            measurements['INDX'][good_mag] += measurements['INDX_DISPCORR'][good_mag]

        #---------------------------------------------------------------
        # Initialize the header keywords
        self.hardcopy = hardcopy
        pri_hdr = self._initialize_primary_header()
        map_hdr = DAPFitsUtil.build_map_header(self.binned_spectra.drpf,
                                               'K Westfall <westfall@ucolick.org>')

        # Get the bin ids with measurements
        bin_indx = DAPFitsUtil.downselect_bins(self.binned_spectra['BINID'].data.ravel(),
                                               measurements['BINID'])

        # Get the spatial map mask
        map_mask = numpy.zeros(self.spatial_shape, dtype=self.bitmask.minimum_dtype())
        # Add any spaxel not used because it was flagged by the binning
        # step
        indx = self.binned_spectra['MAPMASK'].data > 0
        map_mask[indx] = self.bitmask.turn_on(map_mask[indx], 'DIDNOTUSE')
        # Isolate any spaxels with foreground stars
        indx = self.binned_spectra.bitmask.flagged(self.binned_spectra['MAPMASK'].data, 'FORESTAR')
        map_mask[indx] = self.bitmask.turn_on(map_mask[indx], 'FORESTAR')
        # Get the bins that were blow the S/N limit
        indx = numpy.invert(DAPFitsUtil.reconstruct_map(self.spatial_shape,
                                                        self.binned_spectra['BINID'].data.ravel(),
                                                        good_snr, dtype='bool')) & (map_mask == 0)
        map_mask[indx] = self.bitmask.turn_on(map_mask[indx], 'LOW_SNR')

        # Compile the information on the suite of measured indices
        passband_database = self._compile_database()

        # Save the data to the hdu attribute
        self.hdu = fits.HDUList([ fits.PrimaryHDU(header=pri_hdr),
                                  fits.ImageHDU(data=bin_indx.reshape(self.spatial_shape),
                                                header=map_hdr, name='BINID'),
                                  fits.ImageHDU(data=map_mask, header=map_hdr, name='MAPMASK'),
                                  fits.BinTableHDU.from_columns( [ fits.Column(name=n,
                                                    format=rec_to_fits_type(passband_database[n]),
                            array=passband_database[n]) for n in passband_database.dtype.names ],
                                                               name='SIPAR'),
                                  fits.BinTableHDU.from_columns( [ fits.Column(name=n,
                                                    format=rec_to_fits_type(measurements[n]),
                                array=measurements[n]) for n in measurements.dtype.names ],
                                                               name='SINDX')
                                ])

        # Write the data, if requested
        if self.hardcopy:
            if not os.path.isdir(self.directory_path):
                os.makedirs(self.directory_path)
            self.write(clobber=clobber)
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, '-'*50)


#        # Fill array for any missing bins in prep for writing to the HDU
#        bin_indx = self.binned_spectra['BINID'].data.ravel()
#        unique_bins, reconstruct = numpy.unique(bin_indx, return_inverse=True)
#        indx = bin_indx > -1
#
#        flux = numpy.zeros(self.shape, dtype=numpy.float)
#        flux.reshape(-1,self.nwave)[indx,:] = _flux[unique_bins[reconstruct[indx]],:]
#    
#        ivar = numpy.zeros(self.shape, dtype=numpy.float)
#        ivar.reshape(-1,self.nwave)[indx,:] = _ivar[unique_bins[reconstruct[indx]],:]
#   
#        _good_bins = numpy.zeros(numpy.prod(self.spatial_shape), dtype=numpy.bool)
#        _good_bins[indx] = good_bins[unique_bins[reconstruct[indx]]]
#   
#        mask = self._initialize_mask(_good_bins)
#        mask.reshape(-1,self.nwave)[indx,:] = _mask[unique_bins[reconstruct[indx]],:]


    def write(self, clobber=False):
        """
        Write the hdu object to the file.
        """
        DAPFitsUtil.write(self.hdu, self.file_path(), clobber=clobber, checksum=True,
                          loggers=self.loggers, quiet=self.quiet)
        
#        # Restructure the map
#        DAPFitsUtil.restructure_map(self.hdu, ext=self.image_arrays, inverse=True)
#        # Writeh the HDU
#        write_hdu(self.hdu, self.file_path(), clobber=clobber, checksum=True, loggers=self.loggers,
#                  quiet=self.quiet)
#        # Restructure the map
#        DAPFitsUtil.restructure_map(self.hdu, ext=self.image_arrays)


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

#        self.hdu = fits.open(ifile, checksum=checksum)
        self.hdu = DAPFitsUtil.read(ifile, checksum=checksum)

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

#        if not self.quiet:
#            log_output(self.loggers, 1, logging.INFO, 'Reverting to python-native structure.')
#        DAPFitsUtil.restructure_map(self.hdu, ext=self.image_arrays)

        self.nbins = self.hdu['PRIMARY'].header['NBINS']
        self.missing_bins = self._get_missing_bins()


