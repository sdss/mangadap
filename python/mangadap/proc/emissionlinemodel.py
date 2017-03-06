# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
A class hierarchy that fits the emission lines.

*License*:
    Copyright (c) 2015, SDSS-IV/MaNGA Pipeline Group
        Licensed under BSD 3-clause license - see LICENSE.rst

*Source location*:
    $MANGADAP_DIR/python/mangadap/proc/emissionlinefits.py

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
        import os
        import logging
        
        import numpy
        from astropy.io import fits
        import astropy.constants
        
        from ..drpfits import DRPFits
        from ..par.parset import ParSet
        from ..par.artifactdb import ArtifactDB
        from ..par.emissionlinedb import EmissionLineDB
        from ..config.defaults import default_dap_source, default_dap_file_name
        from ..config.defaults import default_dap_method, default_dap_method_path
        from ..util.fitsutil import DAPFitsUtil
        from ..util.fileio import init_record_array, rec_to_fits_type, rec_to_fits_col_dim
        from ..util.bitmask import BitMask
        from ..util.pixelmask import SpectralPixelMask
        from ..util.log import log_output
        from .spatiallybinnedspectra import SpatiallyBinnedSpectra
        from .stellarcontinuummodel import StellarContinuumModel
        from .elric import Elric, ElricPar, GaussianLineProfile
        from .util import _select_proc_method

*Class usage examples*:
    Add examples!

*Revision history*:
    | **26 Apr 2016**: Implementation begun by K. Westfall (KBW)
    | **28 Jul 2016**: (KBW) Fixed error in initialization of guess
        redshift when stellar continuum is provided.
    | **23 Feb 2017**: (KBW) Use DAPFitsUtil read and write functions.

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

from ..drpfits import DRPFits
from ..par.parset import ParSet
from ..par.artifactdb import ArtifactDB
from ..par.emissionlinedb import EmissionLineDB
from ..config.defaults import default_dap_source, default_dap_file_name
from ..config.defaults import default_dap_method, default_dap_method_path
from ..util.fitsutil import DAPFitsUtil
from ..util.fileio import init_record_array, rec_to_fits_type, rec_to_fits_col_dim
from ..util.bitmask import BitMask
from ..util.pixelmask import SpectralPixelMask
from ..util.log import log_output
from ..util.parser import DefaultConfig
from .spatiallybinnedspectra import SpatiallyBinnedSpectra
from .bandpassfilter import emission_line_equivalent_width
from .elric import Elric, ElricPar, GaussianLineProfile
from .util import _select_proc_method

from matplotlib import pyplot

__author__ = 'Kyle B. Westfall'
# Add strict versioning
# from distutils.version import StrictVersion


class EmissionLineModelDef(ParSet):
    """
    A class that holds the parameters necessary to perform the
    emission-line profile fits.

    Args:
        key (str): Keyword used to distinguish between different
            emission-line moment databases.
        minimum_snr (bool): Minimum S/N of spectrum to fit
        artifacts (str): String identifying the artifact database to use
        emission_lines (str): String identifying the emission-line
            database to use
    """
    def __init__(self, key, minimum_snr, artifacts, emission_lines, fitpar, fitclass, fitfunc):
        in_fl = [ int, float ]
        par_opt = [ ParSet, dict ]

        pars =     [ 'key', 'minimum_snr', 'artifacts', 'emission_lines', 'fitpar', 'fitclass',
                     'fitfunc' ]
        values =   [   key,   minimum_snr,   artifacts,   emission_lines,   fitpar,   fitclass,
                       fitfunc ]
        dtypes =   [   str,         in_fl,         str,              str,  par_opt,       None,
                          None ]

        can_call = [ False,         False,      False,    False,      False,      True ]


        ParSet.__init__(self, pars, values=values, dtypes=dtypes)


def validate_emission_line_modeling_method_config(cnfg):
    """ 

    Validate the :class:`mangadap.util.parser.DefaultConfig` with the
    emission-line modeling method parameters.

    Args:
        cnfg (:class:`mangadap.util.parser.DefaultConfig`): Object meant
            to contain defining parameters of the emission-line modeling
            method needed by :class:`EmissionLineModelDef'

    Raises:
        KeyError: Raised if any required keywords do not exist.
        ValueError: Raised if keys have unacceptable values.
    """
    # Check for required keywords
    if 'key' not in cnfg:
        raise KeyError('Keyword \'key\' must be provided.')


def available_emission_line_modeling_methods(dapsrc=None):
    """
    Return the list of available emission-line modeling methods.

    Available methods are:

    .. todo::
        Fill in

    Args:
        dapsrc (str): (**Optional**) Root path to the DAP source
            directory.  If not provided, the default is defined by
            :func:`mangadap.config.defaults.default_dap_source`.

    Returns:
        list: A list of :func:`EmissionLineModelDef` objects, each
        defining an emission-line modeling method.

    Raises:
        NotADirectoryError: Raised if the provided or default
            *dapsrc* is not a directory.
        OSError/IOError: Raised if no emission-line moment configuration
            files could be found.
        KeyError: Raised if the emission-line modeling method keywords
            are not all unique.

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
    search_dir = os.path.join(dapsrc, 'python/mangadap/config/emission_line_modeling')
    ini_files = glob.glob(os.path.join(search_dir, '*.ini'))
    if len(ini_files) == 0:
        raise IOError('Could not find any configuration files in {0} !'.format(search_dir))

    # Build the list of library definitions
    method_list = []
    for f in ini_files:
        # Read the config file
        cnfg = DefaultConfig(f)
        # Ensure it has the necessary elements to define the template
        # library
        validate_emission_line_modeling_method_config(cnfg)

        # Currently only implement Elric for fitting the emission lines
        fitpar = ElricPar(None, cnfg.getint('baseline_order'), cnfg.getfloat('window_buffer'),
                          None, None, minimum_snr=cnfg.getfloat('minimum_snr'))
        fitclass = Elric(EmissionLineModelBitMask(dapsrc=dapsrc))
        fitfunc = fitclass.fit_SpatiallyBinnedSpectra

        method_list += [ EmissionLineModelDef(cnfg['key'], cnfg.getfloat('minimum_snr',default=0.0),
                                              cnfg['artifact_mask'], cnfg['emission_lines'], fitpar,
                                              fitclass, fitfunc) ]

    # Check the keywords of the libraries are all unique
    if len(numpy.unique(numpy.array([method['key'] for method in method_list]))) \
                != len(method_list):
        raise KeyError('Emission-line fitting method keywords are not all unique!')

    # Return the default list of fitting methods
    return method_list


class EmissionLineModelBitMask(BitMask):
    r"""

    Derived class that specifies the mask bits for the emission-line
    modeling.  See :class:`mangadap.util.bitmask.BitMask` for
    attributes.

    A list of the bits and meanings are provided by the base class
    function :func:`mangadap.util.bitmask.BitMask.info`; i.e.,::

        from mangadap.proc.emissionlinemodel import EmissionLineModelBitMask
        bm = EmissionLineModelBitMask()
        bm.info()

    """
    def __init__(self, dapsrc=None):
        dapsrc = default_dap_source() if dapsrc is None else str(dapsrc)
        BitMask.__init__(self, ini_file=os.path.join(dapsrc, 'python', 'mangadap', 'config',
                                                     'bitmasks', 'emission_line_model_bits.ini'))


class EmissionLineModel:
    r"""

    Class that holds the emission-line model fits.

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
    def __init__(self, method_key, binned_spectra, redshift=None, dispersion=None, continuum=None,
                 continuum_method=None, method_list=None, artifact_list=None,
                 emission_line_db_list=None, dapsrc=None, dapver=None, analysis_path=None,
                 directory_path=None, output_file=None, hardcopy=True, clobber=False,
                 checksum=False, loggers=None, quiet=False):

        self.loggers = None
        self.quiet = False

        # Define the database properties
        self.method = None
        self.pixelmask = None
        self.artdb = None
        self.emldb = None
        self._define_method(method_key, method_list=method_list, artifact_list=artifact_list,
                            emission_line_db_list=emission_line_db_list, dapsrc=dapsrc)

        self.neml = self.emldb.nsets

        self.binned_spectra = None
        self.redshift = None
        self.dispersion = None
        self.continuum = None
        self.continuum_method = None

        # Define the output directory and file
        self.directory_path = None      # Set in _set_paths
        self.output_file = None
        self.hardcopy = None

        # Initialize the objects used in the assessments
        self.bitmask = EmissionLineModelBitMask(dapsrc=dapsrc)

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

        self.nmodels = None
        self.missing_models = None

        # Fit the binned spectra
        self.fit(binned_spectra, redshift=redshift, dispersion=dispersion, continuum=continuum,
                 continuum_method=continuum_method, dapsrc=dapsrc, dapver=dapver,
                 analysis_path=analysis_path, directory_path=directory_path,
                 output_file=output_file, hardcopy=hardcopy, clobber=clobber, loggers=loggers,
                 quiet=quiet)


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


    def _define_method(self, method_key, method_list=None, artifact_list=None,
                       emission_line_db_list=None, dapsrc=None):
        r"""
        Select the modeling method
        """
        # Grab the specific method
        self.method = _select_proc_method(method_key, EmissionLineModelDef, method_list=method_list,
                                          available_func=available_emission_line_modeling_methods,
                                          dapsrc=dapsrc)

        # Instantiate the artifact and emission-line databases
        self.artdb = None if self.method['artifacts'] is None else \
                    ArtifactDB(self.method['artifacts'], artdb_list=artifact_list, dapsrc=dapsrc)
        self.pixelmask = SpectralPixelMask(artdb=self.artdb)

        self.emldb = None if self.method['emission_lines'] is None else \
                        EmissionLineDB(self.method['emission_lines'],
                                       emldb_list=emission_line_db_list, dapsrc=dapsrc)


    def _set_paths(self, directory_path, dapver, analysis_path, output_file):
        """
        Set the :attr:`directory_path` and :attr:`output_file`.  If not
        provided, the defaults are set using, respectively,
        :func:`mangadap.config.defaults.default_dap_method_path` and
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
                                    continuum_method=self.continuum_method)
        self.directory_path = default_dap_method_path(method, plate=self.binned_spectra.drpf.plate,
                                                      ifudesign=self.binned_spectra.drpf.ifudesign,
                                                      ref=True,
                                                      drpver=self.binned_spectra.drpf.drpver,
                                                      dapver=dapver, analysis_path=analysis_path) \
                                        if directory_path is None else str(directory_path)

        # Set the output file
        ref_method = '{0}-{1}'.format(self.binned_spectra.rdxqa.method['key'],
                                      self.binned_spectra.method['key'])
        if self.continuum_method is not None:
            ref_method = '{0}-{1}'.format(ref_method, self.continuum_method)
        ref_method = '{0}-{1}'.format(ref_method, self.method['key'])
        self.output_file = default_dap_file_name(self.binned_spectra.drpf.plate,
                                                 self.binned_spectra.drpf.ifudesign, ref_method) \
                                        if output_file is None else str(output_file)


    def _initialize_primary_header(self, hdr=None):
        """

        Initialize the header.

        """
        # Copy the from the DRP and clean it
        if hdr is None:
            hdr = self.binned_spectra.drpf.hdu['PRIMARY'].header.copy()
            hdr = DAPFitsUtil.clean_dap_primary_header(hdr)
        
        hdr['AUTHOR'] = 'Kyle B. Westfall <westfall@ucolick.org>'
        hdr['ELFKEY'] = (self.method['key'], 'Emission-line modeling method keyword')
        hdr['ELFMINSN'] = (self.method['minimum_snr'], 'Minimum S/N of spectrum to include')
        hdr['ARTDB'] = (self.method['artifacts'], 'Artifact database keyword')
        hdr['EMLDB'] = (self.method['emission_lines'], 'Emission-line database keyword')
        if self.continuum_method is not None:
            hdr['SCKEY'] = (self.continuum_method, 'Stellar-continuum model keyword')
        hdr['NELMOD'] = (self.nmodels, 'Number of unique emission-line models')
#        if len(self.missing_models) > 0:
#            hdr['EMPTYEL'] = (str(self.missing_models), 'List of bins w/o EL model')
        # Anything else?
        # Additional database details?
        return hdr


    def _emission_line_database_dtype(self, name_len, mode_len, prof_len):
        r"""
        Construct the record array data type for the output fits
        extension.
        """
        return [ ('ID',numpy.int),
                 ('NAME','<U{0:d}'.format(name_len)),
                 ('RESTWAVE', numpy.float),
                 ('ACTION', '<U1'),
                 ('FLUXRATIO', numpy.float),
                 ('MODE','<U{0:d}'.format(mode_len)),
                 ('PROFILE', '<U{0:d}'.format(prof_len)),
                 ('NCOMP', numpy.int)
               ]


    def _compile_database(self):
        """
        Compile the database with the specifications of each index.
        """
        name_len = numpy.amax([ len(n) for n in self.emldb['name'] ])
        mode_len = numpy.amax([ len(m) for m in self.emldb['mode'] ])
        prof_len = numpy.amax([ len(p) for p in self.emldb['profile'] ])

        # Instatiate the table data that will be saved defining the set
        # of emission-line moments measured
        line_database = init_record_array(self.neml,
                                self._emission_line_database_dtype(name_len, mode_len, prof_len))

        hk = [ 'ID', 'NAME', 'RESTWAVE', 'ACTION', 'FLUXRATIO', 'MODE', 'PROFILE', 'NCOMP' ]
        mk = [ 'index', 'name', 'restwave', 'action', 'flux', 'mode', 'profile', 'ncomp' ]
        for _hk, _mk in zip(hk,mk):
            line_database[_hk] = self.emldb[_mk]

        return line_database


    def _assign_input_kinematics(self, redshift, dispersion, default_dispersion=100.0):
        """
        Set the initial redshift and velocity dispersion for each
        spectrum.  Directly provided redshifts take precedence over
        those in the StellarContinuumModel object, if both are provided.
        The dispersions are **never** drawn from the stellar kinematics,
        but set to default_dispersion if they're not provided.
        """

        # Redshift
        self.redshift = numpy.zeros(self.binned_spectra.nbins, dtype=numpy.float)
        if redshift is not None:
            _redshift = numpy.atleast_1d(redshift)
            if len(_redshift) not in [ 1, self.binned_spectra.nbins ]:
                raise ValueError('Provided redshift must be either a single value or match the ' \
                                 'number of binned spectra.')
            self.redshift = numpy.full(self.binned_spectra.nbins, redshift, dtype=float) \
                                if len(_redshift) == 1 else _redshift.copy()

        # Same for dispersion
        self.dispersion = numpy.full(self.binned_spectra.nbins, default_dispersion,
                                     dtype=numpy.float)
        if dispersion is not None:
            _dispersion = numpy.atleast_1d(dispersion)

            if len(_dispersion) not in [ 1, self.binned_spectra.nbins ]:
                raise ValueError('Provided dispersion must be either a single value or match the ' \
                                 'number of binned spectra.')
            self.dispersion = numpy.full(self.binned_spectra.nbins, dispersion, dtype=float) \
                                if len(_redshift) == 1 else _dispersion.copy()


    def _fill_method_par(self, dapsrc=None, analysis_path=None):
        """

        Fill in any remaining modeling parameters.

        """
        # Fill the guess kinematics
        self.method['fitpar']['guess_redshift'] = self.redshift
        self.method['fitpar']['guess_dispersion'] = self.dispersion
        self.method['fitpar']['emission_lines'] = self.emldb
        self.method['fitpar']['pixelmask'] = self.pixelmask
        if self.continuum is not None:
            self.method['fitpar']['stellar_continuum'] = self.continuum


    def _add_method_header(self, hdr):
        """Add fitting method information to the header."""
        if self.method['fitclass'] is not None:
            try:
                hdr['ELTYPE'] = self.method['fitclass'].fit_type
                hdr['ELMETH'] = self.method['fitclass'].fit_method
            except:
                if not self.quiet and self.hardcopy:
                    warnings.warn('Fit class object does not have fit_type and/or fit_method ' \
                                  'attributes.  No parameters written to header.')
        if self.method['fitpar'] is not None:
            try:
                hdr = self.method['fitpar'].toheader(hdr)
            except:
                if not self.quiet and self.hardcopy:
                    warnings.warn('Fit parameter class has no toheader() function.  No ' \
                                  'parameters written to header.')
        return hdr


    def _finalize_cube_mask(self, mask):
        """
        Finalize the mask by setting the DIDNOTUSE, FORESTAR, and LOW_SNR masks

        Returns:
            numpy.ndarray : Bitmask array.
        """
        # Turn on the flag stating that the pixel wasn't used
        indx = self.binned_spectra.bitmask.flagged(self.binned_spectra.drpf['MASK'].data,
                                                   flag=self.binned_spectra.do_not_fit_flags())
        mask[indx] = self.bitmask.turn_on(mask[indx], 'DIDNOTUSE')

        # Turn on the flag stating that the pixel has a foreground star
        indx = self.binned_spectra.bitmask.flagged(self.binned_spectra.drpf['MASK'].data,
                                                   flag='FORESTAR')
        mask[indx] = self.bitmask.turn_on(mask[indx], 'FORESTAR')

        # Turn on the flag stating that the S/N in the spectrum was
        # below the requested limit
        low_snr = numpy.invert(self.binned_spectra['BINID'].data == self.hdu['BINID'].data)
        indx = numpy.array([low_snr]*self.nwave).transpose(1,2,0)
        mask[indx] = self.bitmask.turn_on(mask[indx], flag='LOW_SNR')

        return mask

    
    def _assign_spectral_arrays(self):
        self.spectral_arrays = [ 'FLUX', 'BASE', 'MASK' ]


    def _assign_image_arrays(self):
        """
        Set :attr:`image_arrays`, which contains the list of extensions
        in :attr:`hdu` that are on-sky image data.
        """
        self.image_arrays = [ 'BINID' ]


    def _get_missing_models(self):
        good_snr = self.binned_spectra.above_snr_limit(self.method['minimum_snr'])
        return numpy.sort(self.binned_spectra['BINS'].data['BINID'][~good_snr].tolist()
                                + self.binned_spectra.missing_bins) 


    def _construct_2d_hdu(self, good_snr, model_flux, model_base, model_mask, model_fit_par,
                          model_eml_par):
        """
        Construct :attr:`hdu` that is held in memory for manipulation of
        the object.  See :func:`construct_3d_hdu` if you want to convert
        the object into a DRP-like datacube.
        """
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, 'Constructing hdu ...')

        # Initialize the headers
        pri_hdr = self._initialize_primary_header()
        pri_hdr = self._add_method_header(pri_hdr)
        map_hdr = DAPFitsUtil.build_map_header(self.binned_spectra.drpf,
                                               'K Westfall <westfall@ucolick.org>')

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

        # Get the bin ids with fitted models
        bin_indx = DAPFitsUtil.downselect_bins(self.binned_spectra['BINID'].data.ravel(),
                                               model_eml_par['BINID'])

        # Compile the information on the suite of measured indices
        line_database = self._compile_database()

        # Save the data to the hdu attribute
        self.hdu = fits.HDUList([ fits.PrimaryHDU(header=pri_hdr),
                                  fits.ImageHDU(data=model_flux.data, name='FLUX'),
                                  fits.ImageHDU(data=model_base.data, name='BASE'),
                                  fits.ImageHDU(data=model_mask, name='MASK'),
                                  self.binned_spectra['WAVE'].copy(),
                                  fits.ImageHDU(data=bin_indx.reshape(self.spatial_shape),
                                                header=map_hdr, name='BINID'),
                                  fits.ImageHDU(data=map_mask, header=map_hdr, name='MAPMASK'),
                                  fits.BinTableHDU.from_columns( [ fits.Column(name=n,
                                                        format=rec_to_fits_type(line_database[n]),
                                        array=line_database[n]) for n in line_database.dtype.names],
                                                         name='PAR'),
                                  fits.BinTableHDU.from_columns([fits.Column(name=n,
                                                        format=rec_to_fits_type(model_fit_par[n]),
                                                        dim=rec_to_fits_col_dim(model_fit_par[n]),
                                        array=model_fit_par[n]) for n in model_fit_par.dtype.names],
                                                                    name='FIT'),
                                  fits.BinTableHDU.from_columns([fits.Column(name=n,
                                                        format=rec_to_fits_type(model_eml_par[n]),
                                                        dim=rec_to_fits_col_dim(model_eml_par[n]),
                                        array=model_eml_par[n]) for n in model_eml_par.dtype.names],
                                                                    name='EMLDATA')
                                ])


    def file_name(self):
        """Return the name of the output file."""
        return self.output_file


    def file_path(self):
        """Return the full path to the output file."""
        if self.directory_path is None or self.output_file is None:
            return None
        return os.path.join(self.directory_path, self.output_file)


    def fit(self, binned_spectra, redshift=None, dispersion=None, continuum=None,
            continuum_method=None, dapsrc=None, dapver=None, analysis_path=None,
            directory_path=None, output_file=None, hardcopy=True, clobber=False, loggers=None,
            quiet=False):
        """

        Fit the emission lines.

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

        # Continuum accounts for underlying absorption
        self.continuum = None
        if continuum is not None:
            if continuum.shape != self.binned_spectra['FLUX'].data.shape:
                raise ValueError('Provided continuum does not match shape of the binned spectra.')
            self.continuum = continuum if isinstance(continuum, numpy.ma.MaskedArray) \
                                    else numpy.ma.MaskedArray(continuum)
            self.continuum_method = continuum_method

        self.shape = self.binned_spectra.shape
        self.spatial_shape =self.binned_spectra.spatial_shape
        self.nspec = self.binned_spectra.nspec
        self.spatial_index = self.binned_spectra.spatial_index.copy()
        self.dispaxis = self.binned_spectra.dispaxis
        self.nwave = self.binned_spectra.nwave
        
        # Get the guess kinematics
        self._assign_input_kinematics(redshift, dispersion)

        #---------------------------------------------------------------
        # Get the good spectra
        good_snr = self.binned_spectra.above_snr_limit(self.method['minimum_snr'])

        # Report
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, '-'*50)
            log_output(self.loggers, 1, logging.INFO, 'EMISSION-LINE PROFILE FITTING:')
            log_output(self.loggers, 1, logging.INFO, '-'*50)
            log_output(self.loggers, 1, logging.INFO, 'Number of binned spectra: {0}'.format(
                                                            self.binned_spectra.nbins))
            if len(self.binned_spectra.missing_bins) > 0:
                log_output(self.loggers, 1, logging.INFO, 'Missing bins: {0}'.format(
                                                            len(self.binned_spectra.missing_bins)))
            log_output(self.loggers, 1, logging.INFO, 'With good S/N and to fit: {0}'.format(
                                                            numpy.sum(good_snr)))
            
        if numpy.sum(good_snr) == 0:
            raise ValueError('No good spectra to fit!')

        #---------------------------------------------------------------
        # Fill in any remaining binning parameters
        self._fill_method_par(dapsrc=dapsrc, analysis_path=analysis_path)

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

        #---------------------------------------------------------------
        # Fit the spectra
        # Mask should be fully defined within the fitting function
        model_wave, model_flux, model_base, model_mask, model_fit_par, model_eml_par = \
            self.method['fitfunc'](self.binned_spectra, par=self.method['fitpar'],
                                   loggers=self.loggers, quiet=self.quiet)

#        pyplot.step(model_wave, model_flux[0,:], where='mid')
#        pyplot.show()
       
        # The number of models returned should match the number of
        # "good" spectra

        # TODO: Include failed fits in "missing" models?

        # DEBUG
        if model_flux.shape[0] != numpy.sum(good_snr):
            raise ValueError('Unexpected returned shape of fitted continuum models.')
       
        #---------------------------------------------------------------
        # Set the number of models and the missing models
        self.nmodels = model_flux.shape[0]
        self.missing_models = self._get_missing_models()

        # Report
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, 'Fitted models: {0}'.format(self.nmodels))
            if len(self.missing_models) > 0:
                log_output(self.loggers, 1, logging.INFO, 'Missing models: {0}'.format(
                                                            len(self.missing_models)))

        # Construct the 2d hdu list that only contains the fitted models
        self.hardcopy = hardcopy
        self._construct_2d_hdu(good_snr, model_flux, model_base, model_mask, model_fit_par,
                               model_eml_par)

        #---------------------------------------------------------------
        # Write the data, if requested
        if self.hardcopy:
            if not os.path.isdir(self.directory_path):
                os.makedirs(self.directory_path)
            self.write(clobber=clobber)
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, '-'*50)


    def construct_3d_hdu(self):
        """
        Reformat the model spectra into a cube matching the shape of
        the DRP fits file.
        """
        # Report
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO,
                       'Constructing emission-line model datacube ...')

        bin_indx = self.hdu['BINID'].data.copy().ravel()
        model_flux = self.hdu['FLUX'].data.copy()
        model_base = self.hdu['BASE'].data.copy()
        model_mask = self.hdu['MASK'].data.copy()

        flux, base, mask = DAPFitsUtil.reconstruct_cube(self.shape, bin_indx,
                                                        [model_flux, model_base, model_mask])

        mask = self._finalize_cube_mask(mask)

        # Primary header is identical regardless of the shape of the
        # extensions
        hdr = self.hdu['PRIMARY'].header.copy()
        cube_hdr = DAPFitsUtil.build_cube_header(self.binned_spectra.drpf,
                                                 'K Westfall <westfall@ucolick.org>')

        # Return the converted hdu without altering the internal hdu
        return fits.HDUList([ fits.PrimaryHDU(header=hdr),
                              fits.ImageHDU(data=flux, header=cube_hdr, name='FLUX'),
                              fits.ImageHDU(data=base, header=cube_hdr, name='BASE'),
                              fits.ImageHDU(data=mask, header=cube_hdr, name='MASK'),
                              self.hdu['WAVE'].copy(),
                              self.hdu['BINID'].copy(),
                              self.hdu['PAR'].copy(),
                              self.hdu['FIT'].copy(),
                              self.hdu['EMLDATA'].copy()
                            ])


    # Exact same function as used by SpatiallyBinnedSpectra
    def write(self, match_DRP=False, clobber=False):
        """
        Write the hdu object to the file.
        """
        # Convert the spectral arrays in the HDU to a 3D cube and write
        # it
        if match_DRP:
            hdu = self.construct_3d_hdu()
            DAPFitsUtil.write(hdu, self.file_path(), clobber=clobber, checksum=True,
                              loggers=self.loggers, quiet=self.quiet)
#            DAPFitsUtil.write_3d_hdu(hdu, self.file_path(), self.binned_spectra.drpf.mode,
#                                     self.spectral_arrays, self.image_arrays, clobber=clobber,
#                                     checksum=True, loggers=self.loggers, quiet=self.quiet)
            return

        DAPFitsUtil.write(self.hdu, self.file_path(), clobber=clobber, checksum=True,
                          loggers=self.loggers, quiet=self.quiet) 
#        # Restructure the spectral arrays as if they're RSS data, and
#        # restructure any maps
#        DAPFitsUtil.restructure_rss(self.hdu, ext=self.spectral_arrays, inverse=True)
#        DAPFitsUtil.restructure_map(self.hdu, ext=self.image_arrays, inverse=True)
#        # Write the HDU
#        write_hdu(self.hdu, self.file_path(), clobber=clobber, checksum=True, loggers=self.loggers,
#                  quiet=self.quiet)
#        # Revert back to the python native storage for internal use
#        DAPFitsUtil.restructure_rss(self.hdu, ext=self.spectral_arrays)
#        DAPFitsUtil.restructure_map(self.hdu, ext=self.image_arrays)


    def read(self, ifile=None, strict=True, checksum=False):
        """
        Read an existing file with an existing set of emission-line
        models.
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
        if self.hdu['PRIMARY'].header['ELFKEY'] != self.method['key']:
            if strict:
                raise ValueError('Keywords in header do not match specified method keywords!')
            else:
                warnings.warn('Keywords in header do not match specified method keywords!')
        # TODO: "strict" should also check other aspects of the file to
        # make sure that the details of the method are also the same,
        # not just the keyword

#        if not self.quiet:
#            log_output(self.loggers, 1, logging.INFO, 'Reverting to python-native structure.')
#        DAPFitsUtil.restructure_rss(self.hdu, ext=self.spectral_arrays)
#        DAPFitsUtil.restructure_map(self.hdu, ext=self.image_arrays)

        # Attempt to read the modeling parameters
        if self.method['fitpar'] is not None and callable(self.method['fitpar'].fromheader):
            self.method['fitpar'].fromheader(self.hdu['PRIMARY'].header)

        self.nmodels = self.hdu['FLUX'].shape[0]
        self.missing_models = self._get_missing_models()


    def copy_to_array(self, ext='FLUX', waverange=None, include_missing=False):
        r"""

        Wrapper for :func:`mangadap.util.fitsutil.DAPFitsUtil.copy_to_array`
        specific for :class:`EmissionLineModel`.

        Return a copy of the selected data array.  The array size is
        always :math:`N_{\rm models} \times N_{\rm wavelength}`; i.e., the
        data is always flattened to two dimensions and the unique
        spectra are selected.

        Args:
            ext (str) : (**Optional**) Name of the extension from which
                to draw the data.  Must be allowed for the current
                :attr:`mode`; see :attr:`data_arrays`.  Default is
                ``'FLUX'``.
            waverange (array-like) : (**Optional**) Two-element array
                with the first and last wavelength to include in the
                computation.  Default is to use the full wavelength
                range.
            include_missing (bool) : (**Optional**) Create an array with
                a size that accommodates the missing models.

        Returns:
            numpy.ndarray: A 2D array with a copy of the data from the
            selected extension.

        """
        return DAPFitsUtil.copy_to_array(self.hdu, ext=ext, allowed_ext=self.spectral_arrays,
                                         waverange=waverange,
                                    missing_bins=self.missing_models if include_missing else None,
                                         nbins=self.nmodels,
                                        unique_bins=DAPFitsUtil.unique_bins(self.hdu['BINID'].data))


    def copy_to_masked_array(self, ext='FLUX', flag=None, waverange=None, include_missing=False):
        r"""
        Wrapper for
        :func:`mangadap.util.fitsutil.DAPFitsUtil.copy_to_masked_array`
        specific for :class:`EmissionLineModel`.

        Return a copy of the selected data array as a masked array.
        This is functionally identical to :func:`copy_to_array`,
        except the output format is a `numpy.ma.MaskedArray`_.  The
        pixels that are considered to be masked can be specified using
        `flag`.
        
        Args:
            ext (str) : (**Optional**) Name of the extension from which
                to draw the data.  Must be allowed for the current
                :attr:`mode`; see :attr:`data_arrays`.  Default is
                `'FLUX'`.
            flag (str or list): (**Optional**) (List of) Flag names that
                are considered when deciding if a pixel should be
                masked.  The names *must* be a valid bit name as defined
                by :attr:`bitmask`.  If not provided, *ANY* non-zero
                mask bit is omitted.
            waverange (array-like) : (**Optional**) Two-element array
                with the first and last wavelength to include in the
                computation.  Default is to use the full wavelength
                range.
            include_missing (bool) : (**Optional**) Create an array with
                a size that accommodates the missing models.

        Returns:
            numpy.ndarray: A 2D array with a copy of the data from the
            selected extension.
        """
        return DAPFitsUtil.copy_to_masked_array(self.hdu, ext=ext, mask_ext='MASK', flag=flag,
                                                bitmask=self.bitmask,
                                                allowed_ext=self.spectral_arrays,
                                                waverange=waverange,
                                    missing_bins=self.missing_models if include_missing else None,
                                                nbins=self.nmodels,
                                        unique_bins=DAPFitsUtil.unique_bins(self.hdu['BINID'].data))


    def fill_to_match(self, binned_spectra, include_base=False):
        """
        Get the emission-line models that match the shape of the
        provided binned_spectra.
        """
        if binned_spectra is self.binned_spectra:
            emission_lines = self.copy_to_array()
            if include_base:
                emission_lines += self.copy_to_array(ext='BASE')
            emission_lines = numpy.ma.MaskedArray(emission_lines)

            # Number of models matches the numbers of bins
            if binned_spectra.nbins == self.nmodels:
                return emission_lines
    
            # Fill in bins with no models with masked zeros
            _emission_lines = numpy.ma.zeros(binned_spectra['FLUX'].data.shape, dtype=float)
            _emission_lines[:,:] = numpy.ma.masked
            for i,j in enumerate(self.hdu['FIT'].data['BINID_INDEX']):
                _emission_lines[j,:] = emission_lines[i,:]
            return _emission_lines

        raise NotImplementedError('Can only match to internal binned_spectra.')

    # TODO: Copy matched_guess_kinematics() from StellarContinuumModel,
    # specifying a specific emission line

