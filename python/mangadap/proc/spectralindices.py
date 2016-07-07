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

        from ..mangafits import MaNGAFits
        from ..par.parset import ParSet
        from ..config.defaults import default_dap_source, default_dap_file_name
        from ..config.defaults import default_dap_method, default_dap_method_path
        from ..config.defaults import default_dap_common_path
        from ..util.instrument import spectral_resolution, match_spectral_resolution
        from ..util.instrument import spectral_coordinate_step
        from ..util.fileio import init_record_array, rec_to_fits_type, write_hdu
        from ..util.bitmask import BitMask
        from .artifactdb import ArtifactDB
        from .absorptionindexdb import AbsorptionIndexDB
        from .bandheadindexdb import BandheadIndexDB
        from .pixelmask import SpectralPixelMask
        from .spatiallybinnedspectra import SpatiallyBinnedSpectra
        from .templatelibrary import TemplateLibrary
        from .stellarcontinuummodel import StellarContinuumModel
        from .bandpassfilter import passband_integral, passband_integrated_width, passband_integrated_mean
        from .util import _select_proc_method, flux_to_fnu


*Class usage examples*:
    Add examples!

*Revision history*:
    | **20 Apr 2016**: Implementation begun by K. Westfall (KBW)
    | **09 May 2016**: (KBW) Add subtraction of emission-line models

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
import os
import logging

import numpy
from astropy.io import fits
import astropy.constants

from ..mangafits import MaNGAFits
from ..par.parset import ParSet
from ..config.defaults import default_dap_source, default_dap_file_name
from ..config.defaults import default_dap_method, default_dap_method_path
from ..config.defaults import default_dap_common_path
from ..util.instrument import spectral_resolution, match_spectral_resolution
from ..util.instrument import spectral_coordinate_step, spectrum_velocity_scale
from ..util.fileio import init_record_array, rec_to_fits_type, write_hdu
from ..util.log import log_output
from ..util.bitmask import BitMask
from .artifactdb import ArtifactDB
from .absorptionindexdb import AbsorptionIndexDB
from .bandheadindexdb import BandheadIndexDB
from .pixelmask import SpectralPixelMask
from .spatiallybinnedspectra import SpatiallyBinnedSpectra
from .templatelibrary import TemplateLibrary
from .stellarcontinuummodel import StellarContinuumModel
from .emissionlinemodel import EmissionLineModel
from .bandpassfilter import passband_integral, passband_integrated_width, passband_integrated_mean
from .bandpassfilter import passband_weighted_mean
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
        key (str): Keyword used to distinguish between different
            spectral-index databases.
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
            defining parameters of the spectral-index database as needed
            by :class:`SpectralIndicesDef`

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
        list: A list of :func:`SpectralIndicesDef` objects, each
        defining a spectral-index database to measure.

    Raises:
        NotADirectoryError: Raised if the provided or default
            *dapsrc* is not a directory.
        OSError/IOError: Raised if no spectral-index configuration files
            could be found.
        KeyError: Raised if the spectral-index database keywords are not
            all unique.
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
            self.blue_empty = SpectralIndices.sideband_pseudocontinua(wave, flux, bluebands,
                                                                      noise=err, log=log)

#        print(self.blue_center)
#        print('')
#        print(self.blue_continuum)
#        print('')

        self.red_center, self.red_continuum, self.red_continuum_err, self.red_incomplete, \
            self.red_empty = SpectralIndices.sideband_pseudocontinua(wave, flux, redbands,
                                                                     noise=err, log=log)
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
            self.blue_empty = SpectralIndices.sideband_pseudocontinua(wave, flux, bluebands,
                                                                      noise=err, log=log)

        self.red_center, self.red_continuum, self.red_continuum_err, self.red_incomplete, \
            self.red_empty = SpectralIndices.sideband_pseudocontinua(wave, flux, redbands,
                                                                     noise=err, log=log)
        
        # Initialize the arrays to hold the numerator and denominator
        blue_n = order == 'b_r'

        n = numpy.ma.zeros(self.nindx, dtype=numpy.float)
        d = numpy.ma.zeros(self.nindx, dtype=numpy.float)
        n[blue_n] = self.blue_continuum[blue_n]
        d[~blue_n] = self.blue_continuum[~blue_n]
        n[~blue_n] = self.red_continuum[~blue_n]
        d[blue_n] = self.red_continuum[blue_n]

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

        self.version = '1.0'
        self.loggers = None
        self.quiet = False

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
        self.emission_line_model = None
        self.nindx = self.count_indices(self.absdb, self.bhddb)

        # Define the output directory and file
        self.directory_path = None      # Set in _set_paths
        self.output_file = None
        self.hardcopy = None
        self.tpl_symlink_dir = None

        # Initialize the objects used in the assessments
        self.bitmask = SpectralIndicesBitMask(dapsrc=dapsrc)
        self.hdu = None
        self.checksum = checksum

        self.shape = None
        self.spatial_shape = None
        self.nspec = None
        self.spatial_index = None
        self.spectral_arrays = None
        self._assign_spectral_arrays()
        self.image_arrays = None
        self._assign_image_arrays()
        self.dispaxis = None
        self.nwave = None

        self.nbins = None
        self.missing_bins = None

        # Run the assessments of the DRP file
        self.measure(binned_spectra, redshift=redshift, stellar_continuum=stellar_continuum,
                     emission_line_model=emission_line_model, dapsrc=dapsrc, dapver=dapver,
                     analysis_path=analysis_path, directory_path=directory_path,
                     output_file=output_file, hardcopy=hardcopy, tpl_symlink_dir=tpl_symlink_dir,
                     clobber=clobber, loggers=loggers, quiet=quiet)


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
        self.pixelmask = SpectralPixelMask(artdb=self.artdb)

        self.absdb = None if self.database['absindex'] is None else \
                AbsorptionIndexDB(self.database['absindex'], indxdb_list=absorption_index_list,
                                  dapsrc=dapsrc)

        self.bhddb = None if self.database['bandhead'] is None else \
                BandheadIndexDB(self.database['bandhead'], indxdb_list=bandhead_index_list,
                                dapsrc=dapsrc)

    @staticmethod
    def count_indices(absdb, bhddb):
        nindx = 0
        if absdb is not None:
            nindx += absdb.nsets
        if bhddb is not None:
            nindx += bhddb.nsets
        return nindx


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
        hdr['SIKEY'] = (self.database['key'], 'Spectral-index database keyword')
        hdr['SIMINSN'] = (self.database['minimum_snr'], 'Minimum S/N of spectrum to include')
        hdr['FWHM'] = (self.database['fwhm'], 'FWHM of index system resolution (ang)')
        hdr['ARTDB'] = (self.database['artifacts'], 'Artifact database keyword')
        hdr['ABSDB'] = (self.database['absindex'], 'Absorption-index database keyword')
        hdr['BHDDB'] = (self.database['bandhead'], 'Bandhead-index database keyword')
        if self.stellar_continuum is not None:
            hdr['SCKEY'] = (self.stellar_continuum.method['key'], 'Stellar-continuum model keyword')
        if self.emission_line_model is not None:
            hdr['EMLKEY'] = (self.emission_line_model.method['key'], 'Emission-line model keyword')
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
        mask = numpy.zeros(self.shape, dtype=self.bitmask.minimum_dtype())

        # Turn on the flag stating that the pixel wasn't used
#        print('Mask size:', self.binned_spectra['MASK'].data.size)
        indx = self.binned_spectra.bitmask.flagged(self.binned_spectra['MASK'].data,
                                                   flag=self.binned_spectra.do_not_fit_flags())
#        print('Masked as DIDNOTUSE:', numpy.sum(indx))
        mask[indx] = self.bitmask.turn_on(mask[indx], 'DIDNOTUSE')

        # Turn on the flag stating that the pixel has a foreground star
        indx = self.binned_spectra.bitmask.flagged(self.binned_spectra['MASK'].data,
                                                   flag='FORESTAR')
#        print('Masked as FORESTAR: ', numpy.sum(indx))
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


    @staticmethod
    def output_dtype(nindx, bitmask=None):
        r"""
        Construct the record array data type for the output fits
        extension.
        """
        return [ ('BIN_INDEX',numpy.int),
                 ('REDSHIFT', numpy.float),
                 ('MASK', numpy.bool if bitmask is None else bitmask.minimum_dtype(), nindx),
                 ('BCEN', numpy.float, nindx), 
                 ('BCONT', numpy.float, nindx), 
                 ('BCONTERR', numpy.float, nindx),
                 ('RCEN', numpy.float, nindx), 
                 ('RCONT', numpy.float, nindx), 
                 ('RCONTERR', numpy.float, nindx), 
                 ('INDX_DISPCORR', numpy.float, nindx), 
                 ('INDX', numpy.float, nindx), 
                 ('INDXERR', numpy.float, nindx)
               ]


    def _assign_spectral_arrays(self):
        self.spectral_arrays = [ 'FLUX', 'IVAR', 'MASK' ]


    def _assign_image_arrays(self):
        self.image_arrays = [ 'BINID' ]


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
#        good_bins = (self._check_snr()) & ~(self._missing_flags())
#        good_bins[2] = False
#        return good_bins


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
            # Make this a function in StellarContinuumModel
            self.redshift = self.stellar_continuum['PAR'].data['KIN'][:,0] \
                                / astropy.constants.c.to('km/s').value


    def _spectra_for_measurements(self):
        """

        Compile the set of spectra for the spectral-index measurements.
        If the input fwhm is > 0, this function will match the spectral
        resolution of the data to the spectral-index system, based on
        the provided FWHM.

        """
        # Get the data arrays
        wave = self.binned_spectra['WAVE'].data
        flux = self.binned_spectra.copy_to_masked_array(flag=self.binned_spectra.do_not_fit_flags())
        ivar = self.binned_spectra.copy_to_masked_array(ext='IVAR',
                                                        flag=self.binned_spectra.do_not_fit_flags())
        indx = self.pixelmask.boolean(wave,nspec=self.nbins)
        flux[indx] = numpy.ma.masked
        ivar[indx] = numpy.ma.masked

        # TODO: Need to understand if baseline effects should be taken
        # out as well...
        if self.emission_line_model is not None:
            eml_model = self.emission_line_model.copy_to_array()
            flux -= eml_model

        mask = numpy.zeros((self.nbins,self.nwave), dtype=self.bitmask.minimum_dtype())
        flux_mask = numpy.ma.getmaskarray(flux)
        mask[flux_mask] = self.bitmask.turn_on(mask[flux_mask], 'DIDNOTUSE')

        if self.database['fwhm'] < 0:
            return flux, ivar, mask

#        pyplot.step(wave, flux[0,:], where='mid', linestyle='-', color='k', lw=1.5)
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

#        pyplot.step(wave, flux[0,:], where='mid', linestyle='-', color='r', lw=0.5)
#        pyplot.show()

        # Use :func:`mangadap.instrument.match_spectral_resolution` to
        # match the spectral resolution of the binned spectra to the
        # spectral-index system
        existing_sres = self.binned_spectra.drpf['SPECRES'].data
        new_sres = wave/self.database['fwhm']
        
        new_flux, sres, sigoff, new_mask, new_ivar \
                = match_spectral_resolution(wave, flux, existing_sres, wave, new_sres, ivar=ivar,
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
        if numpy.sum(new_mask) > 0:
            mask[new_mask > 0] = self.bitmask.turn_on(mask[new_mask > 0], 'SPECRES_LOW')
        return numpy.ma.MaskedArray(flux,mask=mask>0), numpy.ma.MaskedArray(ivar,mask=mask>0), mask


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


    @staticmethod
    def sideband_pseudocontinua(wave, spec, sidebands, noise=None, log=True):
        """Get the side-band integrals in a single spectrum."""

        # Calculate the pseudo-continua in the sidebands
        nbands = sidebands.shape[0]
        
        pseudocontinuum, pseudocontinuum_error \
                    = passband_integrated_mean(wave, spec, passband=sidebands, err=noise, log=log)
        flux_weighted_center, _fwc_err \
                    = passband_weighted_mean(wave, spec, wave, passband=sidebands, log=log)

        # Calculate the fraction of the band that is covered by unmasked
        # pixels
        interval_frac = passband_integrated_width(wave, spec, passband=sidebands, log=log) \
                                / numpy.diff(sidebands, axis=1).ravel()

#        if len(spec.shape) > 1:
#            pyplot.step(wave, spec[0], where='mid', color='k', lw=0.5, zorder=1)
#            pyplot.step(wave, spec[1], where='mid', color='g', lw=0.5, zorder=2)
#        else:
#            pyplot.step(wave, spec, where='mid', color='k', lw=0.5, zorder=1)
#        pyplot.scatter(flux_weighted_center-1.0, pseudocontinuum, marker='.', s=50,
#                       color='b', lw=0, zorder=3)
#        pyplot.show()

        return flux_weighted_center, pseudocontinuum, pseudocontinuum_error, \
                            interval_frac < 1.0, ~(interval_frac > 0.0)


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
        if not isinstance(absdb, AbsorptionIndexDB):
            raise TypeError('Input database must have type AbsorptionIndexDB.')
        if not isinstance(bhddb, BandheadIndexDB):
            raise TypeError('Input database must have type BandheadIndexDB.')

        # Get the number of indices
        nindx = SpectralIndices.count_indices(absdb, bhddb)
        if nindx == 0:
            raise ValueError('No indices to measure!')

        _flux, noise, _redshift = SpectralIndices.check_and_prep_input(wave, flux, ivar=ivar,
                                                                       mask=mask, redshift=redshift,
                                                                       bitmask=bitmask)
        nspec = _flux.shape[0]

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
        abs_fnu = ~absdb.dummy & (absdb['integrand'] == 'fnu')
        good_abs_fnu = numpy.zeros(nindx, dtype=numpy.bool)
        good_abs_fnu[:absdb.nsets][abs_fnu] = True
        abs_flambda = ~absdb.dummy & (absdb['integrand'] == 'flambda')
        good_abs_flambda = numpy.zeros(nindx, dtype=numpy.bool)
        good_abs_flambda[:absdb.nsets][abs_flambda] = True

        bhd_fnu = ~bhddb.dummy & (bhddb['integrand'] == 'fnu')
        good_bhd_fnu = numpy.zeros(nindx, dtype=numpy.bool)
        good_bhd_fnu[absdb.nsets:][bhd_fnu] = True
        bhd_flambda = ~bhddb.dummy & (bhddb['integrand'] == 'flambda')
        good_bhd_flambda = numpy.zeros(nindx, dtype=numpy.bool)
        good_bhd_flambda[absdb.nsets:][bhd_flambda] = True

        # Initialize the output data
        measurements = init_record_array(nspec, SpectralIndices.output_dtype(nindx,bitmask=bitmask))

        # Perform the measurements on each spectrum
        for i in range(nspec):

            print('Measuring spectral indices in spectrum: {0}/{1}'.format(i+1,nspec), end='\r')

            # -----------------------------------
            # Measure the absorption-line indices

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
                                                               good_abs_fnu, err=_noise is not None,
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
                                                               good_bhd_fnu, err=_noise is not None,
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
#            pyplot.scatter(center[0,:], hdu_measurements['BCONT'][i,:], marker='.',
#                           s=100, color='b', lw=0)
#            pyplot.scatter(center[1,:], hdu_measurements['RCONT'][i,:], marker='.',
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


    def _unit_selection(self):
        
        # Flag the indices as either having magnitude or angstrom units
        angu = numpy.full(self.nindx, False, dtype=numpy.bool)
        angu[:self.absdb.nsets] = self.absdb['units'] == 'ang'
        magu = numpy.invert(angu)
        # Bandhead indices are unitless, but the correction follows the
        # same operation as for the absorption-line indices with
        # angstrom units
        magu[self.absdb.nsets:] = False
        angu[self.absdb.nsets:] = True

        return angu, magu


    def _resolution_matched_template_library(self, dapver=None, dapsrc=None):
        """
        Get a version of the template library that has had its
        resolution matched to that of the spectral-index database.
        """
#        return self.stellar_continuum.method['fitpar']['template_library']
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
                                               designation)

        directory_path = default_dap_common_path(plate=self.binned_spectra.drpf.plate,
                                                 ifudesign=self.binned_spectra.drpf.ifudesign,
                                                 drpver=self.binned_spectra.drpf.drpver,
                                                 dapver=dapver, analysis_path=self.analysis_path)

        velscale_ratio = 1 if self.stellar_continuum.method['fitpar']['velscale_ratio'] is None \
                            else self.stellar_continuum.method['fitpar']['velscale_ratio']
        spectral_step=spectral_coordinate_step(wave, log=True) / velscale_ratio

        # Return the template library object
        return TemplateLibrary(self.stellar_continuum.method['fitpar']['template_library_key'],
                               sres=sres, velocity_offset=velocity_offset,
                               spectral_step=spectral_step, log=True, dapsrc=dapsrc,
                               directory_path=directory_path, symlink_dir=self.tpl_symlink_dir,
                               processed_file=processed_file, loggers=self.loggers,
                               quiet=self.quiet)
    

    def _calculate_dispersion_corrections(self, good_bins, dapver=None, dapsrc=None):
        """
        Calculate the dispersion corrections using the best-fitting
        template models.
        """
        if self.stellar_continuum is None:
            return None, None, None

        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO,
                       'Calculating dispersion corrections using stellar continuum model...')

        # Get the wavelength and mask arrays
        wave = self.binned_spectra['WAVE'].data
        flux = self.binned_spectra.copy_to_masked_array(
                                                    flag=self.binned_spectra.do_not_fit_flags())
        mask = numpy.ma.getmaskarray(flux) | self.pixelmask.boolean(wave, nspec=self.nbins)

        # Get the template library
        template_library = None if self.database['fwhm'] < 0  \
                or self.stellar_continuum.method['fitpar']['template_library_key'] is None else \
                    self._resolution_matched_template_library(dapver=dapver, dapsrc=dapsrc)

        # Get the broadened stellar-continuum model
        good_models = good_bins & (self.stellar_continuum['PAR'].data['MASK'] == 0)
#        print('good_bins: {0}'.format(numpy.sum(good_bins)))
#        print('good models: {0}'.format(numpy.sum(good_models)))
        broadened_models = numpy.ma.MaskedArray(self.stellar_continuum.construct_models(
                                                template_library=template_library), mask=mask)
#        print('Broadened models shape: {0}'.format(broadened_models.shape))

        # Get the unbroadened version
        unbroadened_models = numpy.ma.MaskedArray(self.stellar_continuum.construct_models(
                                                            template_library=template_library,
                                                            redshift_only=True),
                                                  mask=mask)
#        print('Models without broadening shape: {0}'.format(unbroadened_models.shape))

#        model_flux = self.stellar_continuum.copy_to_masked_array(
#                                    flag=self.stellar_continuum.all_except_emission_flags())
#        pyplot.step(wave, flux[0,:], where='mid', linestyle='-', color='k', lw=0.5, zorder=1)
#        pyplot.plot(wave, model_flux[0,:], linestyle='-', color='r', lw=1.5, zorder=2, alpha=0.5)
#        pyplot.plot(wave, broadened_models.data[0,:], linestyle='-', color='g', lw=1.0, zorder=3,
#                    alpha=0.5)
#        pyplot.plot(wave, unbroadened_models.data[0,:], linestyle='-', color='b', lw=1.0, zorder=3,
#                    alpha=0.5)
#        pyplot.show()
#        exit()

        # Measure the indices for both sets of spectra
        broadened_model_measurements = init_record_array(self.nbins, self.output_dtype(self.nindx,
                                                                        bitmask=self.bitmask))
        broadened_model_measurements[good_models] \
                = self.measure_indices(self.absdb, self.bhddb, wave, broadened_models[good_models],
                                       redshift=self.redshift[good_models], bitmask=self.bitmask)

        unbroadened_model_measurements = init_record_array(self.nbins, self.output_dtype(
                                                                self.nindx, bitmask=self.bitmask))
        unbroadened_model_measurements[good_models] \
                = self.measure_indices(self.absdb, self.bhddb, wave,
                                       unbroadened_models[good_models],
                                       redshift=self.redshift[good_models], bitmask=self.bitmask)

        # Do not apply the correction if any of the bands were empty or
        # would have had to divide by zero
        bad_measurement = numpy.array([~good_models]*self.nindx).T
        bad_measurement |= self.bitmask.flagged(broadened_model_measurements['MASK'],
                                                flag=[ 'MAIN_EMPTY', 'BLUE_EMPTY', 'RED_EMPTY',
                                                       'DIVBYZERO' ])
        bad_measurement |= self.bitmask.flagged(unbroadened_model_measurements['MASK'],
                                                flag=[ 'MAIN_EMPTY', 'BLUE_EMPTY', 'RED_EMPTY',
                                                       'DIVBYZERO' ])

        print('Number of measurements: ', bad_measurement.size)
        print('Bad measurements: ', numpy.sum(bad_measurement))

        # Determine which indices have good measurements
        angu, magu = self._unit_selection()

        good_ang = ~bad_measurement & (numpy.array([angu]*self.nbins)) \
                                    & (numpy.absolute(broadened_model_measurements['INDX']) > 0)
        good_mag = ~bad_measurement & (numpy.array([magu]*self.nbins))

        print('Good angstrom measurements: ', numpy.sum(good_ang))
        print('Good magnitude measurements: ', numpy.sum(good_mag))

        # Determine and return the corrections to apply
        corrections = numpy.zeros((self.nbins, self.nindx), dtype=numpy.float)
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
        self.stellar_continuum = None
        if stellar_continuum is not None:
            if not isinstance(stellar_continuum, StellarContinuumModel):
                raise TypeError('Provided stellar continuum must have StellarContinuumModel type!')
            if stellar_continuum.hdu is None:
                raise ValueError('Provided StellarContinuumModel is undefined!')
            self.stellar_continuum = stellar_continuum

        # EmissionLineModel object used to subtract emission lines
        if emission_line_model is not None:
            if not isinstance(emission_line_model, EmissionLineModel):
                raise TypeError('Provided emission line models must be of type EmissionLineModel.')
            if emission_line_model.hdu is None:
                raise ValueError('Provided EmissionLineModel is undefined!')
            self.emission_line_model = emission_line_model

        self.shape = self.binned_spectra.shape
        self.spatial_shape =self.binned_spectra.spatial_shape
        self.nspec = self.binned_spectra.nspec
        self.spatial_index = self.binned_spectra.spatial_index.copy()
        self.dispaxis = self.binned_spectra.dispaxis
        self.nwave = self.binned_spectra.nwave
        
        self.nbins = self.binned_spectra.nbins
        self.missing_bins = self.binned_spectra.missing_bins

        # Get the redshifts to apply
        self._assign_redshifts(redshift if self.stellar_continuum is None else None)

        # Report
        good_bins = self._bins_to_measure()
        good_snr = self._check_snr()
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, '-'*50)
            log_output(self.loggers, 1, logging.INFO, 'SPECTRAL-INDEX MEASUREMENTS:')
            log_output(self.loggers, 1, logging.INFO, '-'*50)
            log_output(self.loggers, 1, logging.INFO, 'Total bins: {0}'.format(
                                                            self.binned_spectra.nbins))
            log_output(self.loggers, 1, logging.INFO, 'Missing bins: {0}'.format(
                                                            len(self.binned_spectra.missing_bins)))
            log_output(self.loggers, 1, logging.INFO, 'With good S/N: {0}'.format(
                                                                            numpy.sum(good_snr)))
            log_output(self.loggers, 1, logging.INFO, 'Total spectra to use: {0}'.format(
                                                                            numpy.sum(good_bins)))
            log_output(self.loggers, 1, logging.INFO, 'Number of indices to measure: {0}'.format(
                                                                            self.nindx))

        # Make sure there are good spectra
        if numpy.sum(good_bins) == 0:
            raise ValueError('No good spectra for measurements!')

        # (Re)Set the output paths
        self._set_paths(directory_path, dapver, analysis_path, output_file)

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
        
        # If the file already exists, and not clobbering, just read the
        # file
        if os.path.isfile(ofile) and not clobber:
            if not self.quiet:
                log_output(self.loggers, 1, logging.INFO, 'Using existing file')
            self.read()
            if not self.quiet:
                log_output(self.loggers, 1, logging.INFO, '-'*50)
            return

#        pyplot.scatter(numpy.arange(self.nbins), self.redshift, marker='.', s=50, color='k', lw=0)
#        pyplot.show()
            
        # Get the spectra to use for the measurements
        _flux, _ivar, _mask = self._spectra_for_measurements()

        # Compile the information on the suite of measured indices
        hdu_database = self._compile_database()

        # Instatiate the table data that will be saved with the index
        # measurements
        hdu_measurements = init_record_array(self.nbins, self.output_dtype(self.nindx,
                                                                           bitmask=self.bitmask))

        # Perform the measurements on the galaxy spectra
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO,
                       'Measuring spectral indices in observed spectra...')
        hdu_measurements[good_bins] = self.measure_indices(self.absdb, self.bhddb,
                                                           self.binned_spectra['WAVE'].data, 
                                                           _flux[good_bins,:].copy(),
                                                           ivar=_ivar[good_bins,:].copy(),
                                                           mask=_mask[good_bins,:].copy(),
                                                           redshift=self.redshift[good_bins],
                                                           bitmask=self.bitmask)

        # Get the corrections by performing the measurements on the
        # best-fitting continuum models, with and without the velocity
        # dispersion broadening
        self.tpl_symlink_dir = tpl_symlink_dir
        hdu_measurements['INDX_DISPCORR'], good_ang, good_mag \
                    = self._calculate_dispersion_corrections(good_bins, dapver=dapver,
                                                             dapsrc=dapsrc)

        # Add in the ancillary information and initialize the bins
        hdu_measurements['BIN_INDEX'] = numpy.arange(self.nbins)
        hdu_measurements['REDSHIFT'] = self.redshift
        hdu_measurements['MASK'][~good_bins] \
                = self.bitmask.turn_on(hdu_measurements['MASK'][~good_bins], 'NO_MEASUREMENT')

#        print('After correction and no measurements: {0}/{1}'.format(
#                        numpy.sum(hdu_measurements['MASK'] > 0), hdu_measurements['MASK'].size))

        # Flag bad corrections
        bad_correction = (numpy.array([good_bins]*self.nindx).T) & ~good_ang & ~good_mag
        hdu_measurements['MASK'][bad_correction] \
                = self.bitmask.turn_on(hdu_measurements['MASK'][bad_correction],
                                       'NO_DISPERSION_CORRECTION')

#        print('After bad corrections: {0}/{1}'.format(
#                        numpy.sum(hdu_measurements['MASK'] > 0), hdu_measurements['MASK'].size))

        # Apply the corrections
        hdu_measurements['INDX'][good_ang] *= hdu_measurements['INDX_DISPCORR'][good_ang]
        hdu_measurements['INDXERR'][good_ang] \
                *= numpy.abs(hdu_measurements['INDX_DISPCORR'][good_ang])
        hdu_measurements['INDX'][good_mag] += hdu_measurements['INDX_DISPCORR'][good_mag]

        # Fill array for any missing bins in prep for writing to the HDU
        bin_indx = self.binned_spectra['BINID'].data.ravel()
        unique_bins, reconstruct = numpy.unique(bin_indx, return_inverse=True)
        indx = bin_indx > -1

        flux = numpy.zeros(self.shape, dtype=numpy.float)
        flux.reshape(-1,self.nwave)[indx,:] = _flux[unique_bins[reconstruct[indx]],:]
    
        ivar = numpy.zeros(self.shape, dtype=numpy.float)
        ivar.reshape(-1,self.nwave)[indx,:] = _ivar[unique_bins[reconstruct[indx]],:]
   
        _good_bins = numpy.zeros(numpy.prod(self.spatial_shape), dtype=numpy.bool)
        _good_bins[indx] = good_bins[unique_bins[reconstruct[indx]]]
   
        mask = self._initialize_mask(_good_bins)
        mask.reshape(-1,self.nwave)[indx,:] = _mask[unique_bins[reconstruct[indx]],:]

        # Initialize the header keywords
        hdr = self._clean_drp_header(ext='PRIMARY')
        self._initialize_header(hdr)

        # Save the data to the hdu attribute
        self.hdu = fits.HDUList([ fits.PrimaryHDU(header=hdr),
                                  fits.ImageHDU(data=flux.data,
                                                header=self.binned_spectra['FLUX'].header.copy(),
                                                name='FLUX'),
                                  fits.ImageHDU(data=ivar.data,
                                                header=self.binned_spectra['IVAR'].header.copy(),
                                                name='IVAR'),
                                  fits.ImageHDU(data=mask,
                                                header=self.binned_spectra['MASK'].header.copy(),
                                                name='MASK'),
                                  self.binned_spectra['WAVE'].copy(),
                                  self.binned_spectra['BINID'].copy(),
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
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, '-'*50)


    def write(self, clobber=False):
        """
        Write the hdu object to the file.
        """
        # Restructure the data to match the DRPFits file
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, 'Restructuring data to match DRP.')
        if self.binned_spectra.drpf.mode == 'CUBE':
            MaNGAFits.restructure_cube(self.hdu, ext=self.spectral_arrays, inverse=True)
            MaNGAFits.restructure_map(self.hdu, ext=self.image_arrays, inverse=True)
        elif self.binned_spectra.drpf.mode == 'RSS':
            MaNGAFits.restructure_rss(self.hdu, ext=self.spectral_arrays, inverse=True)

        # Get the output file and determine if it should be compressed
        ofile = self.file_path()
        write_hdu(self.hdu, ofile, clobber=clobber, checksum=True, loggers=self.loggers,
                  quiet=self.quiet)

        # Revert the structure
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, 'Reverting to python-native structure.')
        if self.binned_spectra.drpf.mode == 'CUBE':
            MaNGAFits.restructure_cube(self.hdu, ext=self.spectral_arrays)
            MaNGAFits.restructure_map(self.hdu, ext=self.image_arrays)
        elif self.binned_spectra.drpf.mode == 'RSS':
            MaNGAFits.restructure_rss(self.hdu, ext=self.spectral_arrays)


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

        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, 'Reverting to python-native structure.')
        if self.binned_spectra.drpf.mode == 'CUBE':
            MaNGAFits.restructure_cube(self.hdu, ext=self.spectral_arrays)
            MaNGAFits.restructure_map(self.hdu, ext=self.image_arrays)
        elif self.binned_spectra.drpf.mode == 'RSS':
            MaNGAFits.restructure_rss(self.hdu, ext=self.spectral_arrays)

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



