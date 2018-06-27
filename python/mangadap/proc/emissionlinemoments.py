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

        import glob
        import os
        import logging
        import numpy

        from astropy.io import fits
        import astropy.constants

        from ..config.defaults import dap_source_dir, default_dap_file_name
        from ..config.defaults import default_dap_method, default_dap_method_path
        from ..util.fileio import init_record_array, rec_to_fits_type
        from ..util.bitmask import BitMask
        from ..util.pixelmask import SpectralPixelMask
        from ..util.log import log_output
        from ..util.parser import DefaultConfig
        from ..par.parset import ParSet
        from ..par.artifactdb import ArtifactDB
        from ..par.emissionmomentsdb import EmissionMomentsDB
        from .spatiallybinnedspectra import SpatiallyBinnedSpectra
        from .stellarcontinuummodel import StellarContinuumModel
        from .bandpassfilter import passband_integral, passband_integrated_width, passband_median
        from .bandpassfilter import passband_integrated_mean, passband_weighted_mean
        from .util import select_proc_method

*Class usage examples*:
    Add examples!

*Revision history*:
    | **25 Apr 2016**: Implementation begun by K. Westfall (KBW)
    | **20 May 2016**: (KBW) Added loggers and quiet keyword arguments
        to :class:`EmissionLineMoments`, removed verbose 
    | **28 Jul 2016**: (KBW) Fixed error in initialization of guess
        redshift when stellar continuum is provided.
    | **09 Jan 2017**: (KBW) Generalized so that it no longer reads in
        StellarContinuumModel object.
    | **23 Feb 2017**: (KBW) Use DAPFitsUtil read and write functions.
    | **27 Feb 2017**: (KBW) Use DefaultConfig for ini files.
    | **31 May 2017**: (KBW) Revert to using
        :class:`mangadap.proc.stellarcontinuummodel.StellarContinuumModel`
        on input
    | **02 Feb 2018**: (KBW) Adjust for change to
        :func:`mangadap.proc.stellarcontinuummodel.StellarContinuumModel.fill_to_match`.
    | **15 Feb 2018**: (KBW) Add parameter that will set whether or not
        the moments should be remeasured after the emission-line
        modeling.  Allow to pass an emission-line model for setting up
        the continuum and the velocities to measure.

.. _astropy.io.fits.hdu.hdulist.HDUList: http://docs.astropy.org/en/v1.0.2/io/fits/api/hdulists.html
.. _glob.glob: https://docs.python.org/3.4/library/glob.html
.. _logging.Logger: https://docs.python.org/3/library/logging.html


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

from ..config.defaults import dap_source_dir, default_dap_file_name
from ..config.defaults import default_dap_method, default_dap_method_path
from ..util.fitsutil import DAPFitsUtil
from ..util.fileio import init_record_array, rec_to_fits_type
from ..util.bitmask import BitMask
from ..util.pixelmask import SpectralPixelMask
from ..util.log import log_output
from ..util.parser import DefaultConfig
from ..par.parset import ParSet
from ..par.artifactdb import ArtifactDB
from ..par.emissionmomentsdb import EmissionMomentsDB
from .spatiallybinnedspectra import SpatiallyBinnedSpectra
from .stellarcontinuummodel import StellarContinuumModel
from .bandpassfilter import passband_integral, passband_integrated_width, passband_median
from .bandpassfilter import passband_integrated_mean, passband_weighted_mean, pseudocontinuum
from .bandpassfilter import emission_line_equivalent_width
from .spectralfitting import EmissionLineFit
from .util import select_proc_method

from matplotlib import pyplot
#from memory_profiler import profile

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
    def __init__(self, key, minimum_snr, artifacts, passbands, redo_postmodeling, fit_vel_name):
        in_fl = [ int, float ]

        pars =     [ 'key', 'minimum_snr', 'artifacts', 'passbands', 'redo_postmodeling',
                     'fit_vel_name' ]
        values =   [ key, minimum_snr, artifacts, passbands, redo_postmodeling, fit_vel_name ]
        dtypes =   [ str, in_fl, str, str, bool, str ]

        ParSet.__init__(self, pars, values=values, dtypes=dtypes)


def validate_emission_line_moments_config(cnfg):
    """ 
    Validate the :class:`mangadap.util.parser.DefaultConfig` with the
    emission-line moment measurement parameters.

    Args:
        cnfg (:class:`mangadap.util.parser.DefaultConfig`): Object meant
            to contain defining parameters of the emission-line moments
            as needed by :class:`EmissionLineMomentsDef`

    Raises:
        KeyError: Raised if any required keywords do not exist.
        ValueError: Raised if keys have unacceptable values.
        FileNotFoundError: Raised if a file is specified but could not
            be found.
    """
    # Check for required keywords
    required_keywords = ['key', 'emission_passbands']
    if not cnfg.all_required(required_keywords):
        raise KeyError('Keywords {0} must all have valid values.'.format(required_keywords))


def available_emission_line_moment_databases(dapsrc=None):
    """
    Return the list of available emission-line moment databases.

    Available database combinations:

    .. todo::
        Fill in

    Args:
        dapsrc (str): (**Optional**) Root path to the DAP source
            directory.  If not provided, the default is defined by
            :func:`mangadap.config.defaults.dap_source_dir`.

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

    .. todo::
        - Somehow add a python call that reads the databases and
          constructs the table for presentation in sphinx so that the
          text above doesn't have to be edited with changes in the
          available databases.
        
    """
    # Check the source directory exists
    dapsrc = dap_source_dir() if dapsrc is None else str(dapsrc)
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
        cnfg = DefaultConfig(f)
        # Ensure it has the necessary elements to define the template
        # library
        validate_emission_line_moments_config(cnfg)

        moment_set_list += [ EmissionLineMomentsDef(cnfg['key'],
                                                    cnfg.getfloat('minimum_snr', default=0.),
                                                    cnfg.get('artifact_mask'),
                                                    cnfg['emission_passbands'],
                                                    cnfg.getbool('redo_postmodeling',
                                                                 default=False),
                                                    cnfg.get('fit_vel_name')) ]

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
        dapsrc = dap_source_dir() if dapsrc is None else str(dapsrc)
        BitMask.__init__(self, ini_file=os.path.join(dapsrc, 'python', 'mangadap', 'config',
                                                     'bitmasks', 'emission_line_moments_bits.ini'))


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
    def __init__(self, database_key, binned_spectra, stellar_continuum=None,
                 emission_line_model=None, redshift=None, database_list=None, artifact_list=None,
                 bandpass_list=None, dapsrc=None, dapver=None, analysis_path=None,
                 directory_path=None, output_file=None, hardcopy=True, clobber=False,
                 checksum=False, loggers=None, quiet=False):

        self.loggers = None
        self.quiet = False

        # Define the database properties
        self.database = None
        self.pixelmask = None
        self.artdb = None
        self.momdb = None
        self._define_databases(database_key, database_list=database_list,
                               artifact_list=artifact_list, bandpass_list=bandpass_list)

        self.nmom = self.momdb.nsets

        self.binned_spectra = None
        self.stellar_continuum = None
        self.emission_line_model = None
        self.redshift = None

        # Define the output directory and file
        self.directory_path = None      # Set in _set_paths
        self.output_file = None
        self.hardcopy = None

        # Initialize the objects used in the assessments
        self.bitmask = EmissionLineMomentsBitMask(dapsrc=dapsrc)

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
                     emission_line_model=emission_line_model, redshift=redshift, dapsrc=dapsrc,
                     dapver=dapver, analysis_path=analysis_path, directory_path=directory_path,
                     output_file=output_file, hardcopy=hardcopy, clobber=clobber, loggers=loggers,
                     quiet=quiet)


    def __getitem__(self, key):
        return self.hdu[key]


    def _define_databases(self, database_key, database_list=None, artifact_list=None,
                          bandpass_list=None, dapsrc=None):
        r"""
        Select the database of bandpass filters.
        """
        # Grab the specific database
        self.database = select_proc_method(database_key, EmissionLineMomentsDef,
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
        :func:`mangadap.config.defaults.default_dap_method_path` and
        :func:`mangadap.config.defaults.default_dap_file_name`.

        ..warning::

            File name is different if an emission-line model is
            provided!

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
        continuum_method = None if self.stellar_continuum is None \
                                else self.stellar_continuum.method['key']
        method = default_dap_method(binned_spectra=self.binned_spectra,
                                    continuum_method=continuum_method)
        self.directory_path = default_dap_method_path(method, plate=self.binned_spectra.drpf.plate,
                                                      ifudesign=self.binned_spectra.drpf.ifudesign,
                                                      ref=True,
                                                      drpver=self.binned_spectra.drpf.drpver,
                                                      dapver=dapver, analysis_path=analysis_path) \
                                        if directory_path is None else str(directory_path)

        # Set the output file
        ref_method = '{0}-{1}'.format(self.binned_spectra.rdxqa.method['key'],
                                      self.binned_spectra.method['key'])
        if continuum_method is not None:
            ref_method = '{0}-{1}'.format(ref_method, continuum_method)
        if self.emission_line_model is not None:
            ref_method = '{0}-{1}'.format(ref_method, self.emission_line_model.method['key'])
        ref_method = '{0}-{1}'.format(ref_method, self.database['key'])
        self.output_file = default_dap_file_name(self.binned_spectra.drpf.plate,
                                                 self.binned_spectra.drpf.ifudesign, ref_method) \
                                        if output_file is None else str(output_file)


    def _initialize_primary_header(self, hdr=None, measurements_binid=None):
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
        hdr['ELMKEY'] = (self.database['key'], 'Emission-line moments database keyword')
        hdr['ELMMINSN'] = (self.database['minimum_snr'], 'Minimum S/N of spectrum to include')
        hdr['ARTDB'] = (self.database['artifacts'], 'Artifact database keyword')
        hdr['MOMDB'] = (self.database['passbands'], 'Emission-line moments database keyword')
        if self.stellar_continuum is not None:
            hdr['SCKEY'] = (self.stellar_continuum.method['key'], 'Stellar-continuum model keyword')
        if self.emission_line_model is not None:
            hdr['ELFKEY'] = (self.emission_line_model.method['key'],
                                'Emission-line modeling method keyword')
        hdr['NBINS'] = (self.nbins, 'Number of unique spectra for measurements')
        hdr['ELMREBIN'] = (measurements_binid is not None, 'Bin IDs disconnected from SC binning')
        return hdr


    def _moments_database_dtype(self, name_len):
        r"""
        Construct the record array data type for the output fits
        extension.
        """
        return [ ('ID',numpy.int),
                 ('NAME','<U{0:d}'.format(name_len)),
                 ('RESTWAVE', numpy.float),
                 ('PASSBAND', numpy.float, (2,)),
                 ('BLUEBAND', numpy.float, (2,)),
                 ('REDBAND', numpy.float, (2,))
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
        passband_database = init_record_array(self.nmom, self._moments_database_dtype(name_len))

        hk = [ 'ID', 'NAME', 'RESTWAVE', 'PASSBAND', 'BLUEBAND', 'REDBAND' ]
        mk = [ 'index', 'name', 'restwave', 'primary', 'blueside', 'redside' ]
        for _hk, _mk in zip(hk,mk):
            passband_database[_hk] = self.momdb[_mk]

        return passband_database


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

        If self.emission_line_model and redshift are None, the default
        redshift is used (or 0.0 if this is also None).

        To get the emission_line stellar kinematics, the function calls
        :func:`mangadap.proc.emissionlinemodel.EmissionLineModel.matched_kinematics`.
        The method must have a fit_vel_name define for this use!

        In this function, the provided redshift must be a single value
        or None; therefore, the means of any vectors should be provided
        instead of the full vector.

        Args:

            redshift (float, numpy.ndarray): Redshifts (:math:`z`) to
                use for each spectrum.  If None, the default

            measure_on_unbinned_spaxels (bool): Flag that method expects
                to measure moments on unbinned spaxels.

            good_snr (numpy.ndarray): Boolean array setting which
                spectra have sufficient S/N for the measurements.

            default_redshift (float): (**Optional**) Only used if there
                are stellar kinematics available.  Provides the default
                redshift to use for spectra without stellar
                measurements; see :arg:`redshift` in
                :func:`mangadap.proc.stellarcontinuummodel.StellarContinuumModel.matched_kinematics`.
                If None (default), the median of the unmasked stellar
                velocities will be used.

        """

        # Construct the binid matrix if measuring on unbinned spaxels
        if measure_on_unbinned_spaxels:
            binid = numpy.full(self.spatial_shape, -1, dtype=int)
            binid.ravel()[good_snr] = numpy.arange(self.nbins)
            missing = []
            nspec = self.binned_spectra.drpf.nspec
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
            fgoodpix = self.binned_spectra.check_fgoodpix()
            good_snr = self.binned_spectra.rdxqa['SPECTRUM'].data['SNR'] \
                                > self.database['minimum_snr']
            return fgoodpix & good_snr
        return self.binned_spectra.above_snr_limit(self.database['minimum_snr'])


    @staticmethod
    def output_dtype(nmom, bitmask=None):
        r"""
        Construct the record array data type for the output fits
        extension.

        Returned columns are:

        +-------------+--------------------------------------------+
        |      Column | Description                                |
        +=============+============================================+
        |       BINID | Bin ID of the fitted spectrum              |
        +-------------+--------------------------------------------+
        | BINID_INDEX | Index of the bin ID of the fitted spectrum |
        +-------------+--------------------------------------------+
        |    REDSHIFT | Redshift used for setting the passbands    |
        +-------------+--------------------------------------------+
        |        MASK | Bit mask value for the measurements        |
        +-------------+--------------------------------------------+
        |        BCEN | Blueband center for moments                |
        +-------------+--------------------------------------------+
        |       BCONT | Blueband pseudocontinuum for moments       |
        +-------------+--------------------------------------------+
        |    BCONTERR | Error in the above                         |
        +-------------+--------------------------------------------+
        |        RCEN | Redband center for moments                 |
        +-------------+--------------------------------------------+
        |       RCONT | Redband pseudocontinuum for moments        |
        +-------------+--------------------------------------------+
        |    RCONTERR | Error in the above                         |
        +-------------+--------------------------------------------+
        |    CNTSLOPE | Continuum slope used for the moments       |
        +-------------+--------------------------------------------+
        |        FLUX | Zeroth moment                              |
        +-------------+--------------------------------------------+
        |     FLUXERR | Zeroth moment error                        |
        +-------------+--------------------------------------------+
        |        MOM1 | First velocity moment                      |
        +-------------+--------------------------------------------+
        |     MOM1ERR | First velocity moment error                |
        +-------------+--------------------------------------------+
        |        MOM2 | Second velocity moment                     |
        +-------------+--------------------------------------------+
        |     MOM2ERR | Second velocity moment error               |
        +-------------+--------------------------------------------+
        |       SINST | Instrumental dispersion at mom1            |
        +-------------+--------------------------------------------+
        |        BMED | Median flux in the blue band used for EW   |
        +-------------+--------------------------------------------+
        |        RMED | Median flux in the red band used for EW    |
        +-------------+--------------------------------------------+
        |      EWCONT | Continuum value used for EW                |
        +-------------+--------------------------------------------+
        |          EW | Equivalenth width (FLUX/pseudo continuum)  |
        +-------------+--------------------------------------------+
        |       EWERR | Error in the above                         |
        +-------------+--------------------------------------------+
        """
        return [ ('BINID',numpy.int),
                 ('BINID_INDEX',numpy.int),
                 ('REDSHIFT', numpy.float),
                 ('MASK', numpy.bool if bitmask is None else bitmask.minimum_dtype(), (nmom,)),
                 ('BCEN', numpy.float, (nmom,)), 
                 ('BCONT', numpy.float, (nmom,)), 
                 ('BCONTERR', numpy.float, (nmom,)),
                 ('RCEN', numpy.float, (nmom,)), 
                 ('RCONT', numpy.float, (nmom,)), 
                 ('RCONTERR', numpy.float, (nmom,)), 
                 ('CNTSLOPE', numpy.float, (nmom,)),
                 ('FLUX', numpy.float, (nmom,)), 
                 ('FLUXERR', numpy.float, (nmom,)),
                 ('MOM1', numpy.float, (nmom,)), 
                 ('MOM1ERR', numpy.float, (nmom,)),
                 ('MOM2', numpy.float, (nmom,)), 
                 ('MOM2ERR', numpy.float, (nmom,)),
                 ('SINST', numpy.float, (nmom,)),
                 ('BMED', numpy.float, (nmom,)), 
                 ('RMED', numpy.float, (nmom,)), 
                 ('EWCONT', numpy.float, (nmom,)), 
                 ('EW', numpy.float, (nmom,)), 
                 ('EWERR', numpy.float, (nmom,))
               ]


    @staticmethod
    def sideband_pseudocontinua(wave, spec, sidebands, spec_n=None, err=None, log=True,
                                weighted_center=True):
        """
        Get the side-band integrals.

        err is a single vector that is the same for all spec!
        """

        if len(spec.shape) == 1:
            return pseudocontinuum(wave, spec, passband=sidebands, err=err, log=log,
                                   weighted_center=weighted_center)

        # Calculate the pseudo-continua in the sidebands
        nspec = spec.shape[0]
        nbands = sidebands.shape[0]

        # More than one spectrum entered, use spec_n to decide which
        # spectrum to use for the integrals.  If spec_n is None and
        # there is more than one spectrum (improper use of the code!),
        # allow the code to continue but issue a warning
        if spec_n is None:
            warnings.warn('Multiple spectra entered, but no designation for which to use!')
            _spec_n = numpy.zeros(nbands, dtype=numpy.int)
        else:
            _spec_n = spec_n
        
        # TODO: Change this to a warning?
        if numpy.amax(_spec_n) > nspec-1:
            raise ValueError('Selected spectrum not provided!')

        # Instantiate the arrays
        continuum = numpy.zeros(nbands, dtype=numpy.float)
        continuum_error = None if err is None else numpy.zeros(nbands, dtype=numpy.float)
        band_center = numpy.zeros(nbands, dtype=numpy.float)
        band_incomplete = numpy.zeros(nbands, dtype=bool)
        band_empty = numpy.zeros(nbands, dtype=bool)

        for i in range(nspec):
            indx = _spec_n==i
            if numpy.sum(indx) == 0:
                continue
            band_center[indx], continuum[indx], _continuum_error, band_incomplete[indx], \
                    band_empty[indx] = pseudocontinuum(wave, spec[i,:], passband=sidebands[indx,:],
                                                       err=err, log=log,
                                                       weighted_center=weighted_center)
            if err is not None:
                continuum_error[indx] = _continuum_error
        return band_center, continuum, continuum_error, band_incomplete, band_empty


    @staticmethod
    def single_band_moments(wave, spec, passband, restwave, err=None):
        """
        Measure the moments for a single band.

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
            undefined 2nd moment
        """
        # Get the fraction of the passband that is unmasked
        interval_frac = passband_integrated_width(wave, spec, passband=passband, log=True) \
                                / numpy.diff(passband).ravel()
        incomplete = interval_frac < 1.0
        empty = numpy.invert(interval_frac > 0.0)
        if empty:
            return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, incomplete, empty, False, False

        # Get the redshift vector
        cz = astropy.constants.c.to('km/s').value*(wave/restwave-1)

        # Get the integrated flux
        flux = passband_integral(wave, spec, passband=passband, log=True)
        fluxerr = None if err is None else \
                    numpy.ma.sqrt( passband_integral(wave, numpy.square(err), passband=passband,
                                                     log=True))

        # No flux in the passband, meaning division by zero will occur!
        if not numpy.absolute(flux) > 0.0:
            return flux, fluxerr, 0.0, 0.0, 0.0, 0.0, incomplete, empty, True, False

        # Get the first and second moments; saved to a temporary
        # variable for use in the error calculations (if possible)
        mom1 = passband_integral(wave, cz*spec, passband=passband, log=True)
        _mom1 = mom1 / flux
        mom2 = passband_integral(wave, numpy.square(cz)*spec, passband=passband, log=True)
        _mom2 = mom2 / flux

        # Perform the error calculations
        if err is None:
            mom1 = _mom1
            mom1err = None

            if _mom2 < numpy.square(mom1):
                return flux, fluxerr, mom1, mom1err, 0.0, 0.0, incomplete, empty, False, True

            mom2 = numpy.ma.sqrt(_mom2 - numpy.square(mom1))
            mom2err = None
        else:
            mom1err = passband_integral(wave, cz*numpy.square(err), passband=passband, log=True)
            mom1err = numpy.ma.sqrt(mom1err * numpy.square(_mom1/mom1)
                                        + numpy.square(fluxerr * _mom1/flux))
            mom1 = _mom1
                
            if _mom2 < numpy.square(mom1):
                return flux, fluxerr, mom1, mom1err, 0.0, 0.0, incomplete, empty, False, True

            mom2err = passband_integral(wave, numpy.square(cz*err), passband=passband, log=True)
            mom2err = numpy.ma.sqrt(mom2err * numpy.square(_mom2/mom2)
                                        + numpy.square(fluxerr * _mom2/flux))
            mom2err = numpy.ma.sqrt(numpy.square(mom2err) + numpy.square(2*mom1*mom1err))
            mom2 = numpy.ma.sqrt(_mom2 - numpy.square(mom1))
            mom2err /= (2.0*mom2)

        return flux, fluxerr, mom1, mom1err, mom2, mom2err, incomplete, empty, False, False


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
        flux = numpy.zeros(nbands, dtype=numpy.float)
        fluxerr = numpy.zeros(nbands, dtype=numpy.float)
        mom1 = numpy.zeros(nbands, dtype=numpy.float)
        mom1err = numpy.zeros(nbands, dtype=numpy.float)
        mom2 = numpy.zeros(nbands, dtype=numpy.float)
        mom2err = numpy.zeros(nbands, dtype=numpy.float)
        incomplete = numpy.zeros(nbands, dtype=numpy.bool)
        empty = numpy.zeros(nbands, dtype=numpy.bool)
        divbyzero = numpy.zeros(nbands, dtype=numpy.bool)
        undefined_mom2 = numpy.zeros(nbands, dtype=numpy.bool)

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
                _spec_n = numpy.zeros(self.nmom, dtype=numpy.int)
            else:
                _spec_n = spec_n

        # Get the moments for all the main passbands
        for i,p in enumerate(mainbands):
            cont = cntb[i] + cntm[i]*wave
            _spec = spec - cont if len(spec.shape) == 1 else spec[_spec_n[i],:] - cont

            flux[i], fluxerr[i], mom1[i], mom1err[i], mom2[i], mom2err[i], incomplete[i], \
                    empty[i], divbyzero[i], undefined_mom2[i] \
                        = EmissionLineMoments.single_band_moments(wave, _spec, p, restwave[i],
                                                                  err=err)

        # Get the instrumental dispersion at the position of the line
        # center
        sinst = numpy.zeros(nbands, dtype=numpy.float) if sres is None \
                        else EmissionLineFit.instrumental_dispersion(wave, sres, restwave, mom1)

        return cntm, cntb, flux, fluxerr, mom1, mom1err, mom2, mom2err, sinst, incomplete, empty, \
                    divbyzero, undefined_mom2


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
        # undefined
        _fluxnc, no_continuum = EmissionLineFit.subtract_continuum(_flux, _continuum)

        # Initialize the output data
        nspec = _flux.shape[0]
        nmom = momdb.nsets
        measurements = init_record_array(nspec, EmissionLineMoments.output_dtype(nmom,
                                                                                 bitmask=bitmask))

        # Mask dummy lines
        if numpy.any(momdb.dummy):
            measurements['MASK'][:,momdb.dummy] = True if bitmask is None else \
                    bitmask.turn_on(measurements['MASK'][:,momdb.dummy], 'UNDEFINED_BANDS')
        
        # No valid lines
        if numpy.all(momdb.dummy):
            return measurements

        # Common arrays used for each spectrum
        blue_fraction = numpy.zeros(nmom, dtype=numpy.float)
        red_fraction = numpy.zeros(nmom, dtype=numpy.float)

        incomplete = numpy.zeros(nmom, dtype=numpy.bool)
        empty = numpy.zeros(nmom, dtype=numpy.bool)
        divbyzero = numpy.zeros(nmom, dtype=numpy.bool)
        undefined_mom2 = numpy.zeros(nmom, dtype=numpy.bool)

        # Perform the measurements for each spectrum
        for i in range(nspec):
            print('Measuring emission-line moments in spectrum: {0}/{1}'.format(i+1,nspec),
                  end='\r')

            incomplete[:] = False
            empty[:] = False
            divbyzero[:] = False
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
                spec = _flux[i,:] if _fluxnc is None else _fluxnc[i,:]
                spec_n = None

                # Flag moment calculations as having not corrected for
                # stellar aborption
                if _fluxnc is None:
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
                spec_n = numpy.zeros(nmom, dtype=numpy.int)
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
                                                        spec_n=spec_n[numpy.invert(momdb.dummy)],
                                                                  err=__err, weighted_center=False)
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
                                                          spec_n=spec_n[numpy.invert(momdb.dummy)],
                                                                  err=__err, weighted_center=False)
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
                    incomplete[indx], empty[indx], divbyzero[indx], undefined_mom2[indx] = \
                            EmissionLineMoments.continuum_subtracted_moments(wave, spec,
                                                               _mainbands, momdb['restwave'][indx],
                                                               measurements['BCEN'][i,indx],
                                                               measurements['BCONT'][i,indx],
                                                               measurements['RCEN'][i,indx],
                                                               measurements['RCONT'][i,indx],
                                                               spec_n=spec_n[indx], err=__err,
                                                               sres=__sres) 

            # Turn on any necessary bits
            measurements['MASK'][i,incomplete] = True if bitmask is None \
                        else bitmask.turn_on(measurements['MASK'][i,incomplete], 'MAIN_INCOMP')
            measurements['MASK'][i,empty] = True if bitmask is None \
                        else bitmask.turn_on(measurements['MASK'][i,empty], 'MAIN_EMPTY')
            measurements['MASK'][i,divbyzero] = True if bitmask is None \
                        else bitmask.turn_on(measurements['MASK'][i,divbyzero], 'DIVBYZERO')
            measurements['MASK'][i,undefined_mom2] = True if bitmask is None \
                        else bitmask.turn_on(measurements['MASK'][i,undefined_mom2],
                                             'UNDEFINED_MOM2')

#        if len(spec.shape) > 1:
#            pyplot.step(wave, spec[0], where='mid', color='k', lw=0.5, zorder=1)
#            pyplot.step(wave, spec[1], where='mid', color='g', lw=0.5, zorder=2)
#        else:
#            pyplot.step(wave, spec, where='mid', color='k', lw=0.5, zorder=1)
#        pyplot.scatter(measurements['BCEN'], measurements['BCONT'], marker='.', s=50,
#                       color='b', lw=0, zorder=3)
#        pyplot.scatter(measurements['RCEN'], measurements['RCONT'], marker='.', s=50,
#                       color='r', lw=0, zorder=3)
#        pyplot.show()

        print('Measuring emission-line moments in spectrum: {0}/{0}'.format(nspec))
        return measurements
        

    def file_name(self):
        """Return the name of the output file."""
        return self.output_file


    def file_path(self):
        """Return the full path to the output file."""
        if self.directory_path is None or self.output_file is None:
            return None
        return os.path.join(self.directory_path, self.output_file)


    def measure(self, binned_spectra, stellar_continuum=None, emission_line_model=None,
                redshift=None, dapsrc=None, dapver=None, analysis_path=None, directory_path=None,
                output_file=None, hardcopy=True, clobber=False, loggers=None, quiet=False):
        """
        Measure the emission-line moments using the binned spectra.

        If neither stellar-continuum nor emission-line models are provided:
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
                and self.emission_line_model.method['deconstruct_bins']

        self.spatial_shape =self.binned_spectra.spatial_shape
        self.nspec = self.binned_spectra.drpf.nspec if measure_on_unbinned_spaxels \
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
            log_output(loggers, 1, logging.INFO, '{0:^50}'.format('EMISSION-LINE MOMENTS'))
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
            log_output(self.loggers, 1, logging.INFO, 'Number of spectra: {0}'.format(self.nspec))
            if not measure_on_unbinned_spaxels and len(self.binned_spectra.missing_bins) > 0:
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
            log_output(self.loggers,1,logging.INFO, 'Output path: {0}'.format(self.directory_path))
            log_output(self.loggers,1,logging.INFO, 'Output file: {0}'.format(self.output_file))
        
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
        # Get the spectra to use for the measurements
        wave, flux, ivar = EmissionLineFit.get_spectra_to_fit(self.binned_spectra.drpf
                                                                if measure_on_unbinned_spaxels
                                                                else self.binned_spectra,
                                                              pixelmask=self.pixelmask,
                                                              select=good_snr)
        # Make sure that the spectra have been corrected for Galactic
        # extinction
        if measure_on_unbinned_spaxels:
            flux, ivar = binned_spectra.galext.apply(flux, ivar=ivar, deredden=True)

        # Get the continuum if there is any
        if eml_stellar_continuum_available:
            if measure_on_unbinned_spaxels:
                binid = numpy.full(self.binned_spectra.spatial_shape, -1, dtype=int)
                binid.ravel()[good_snr] = numpy.arange(self.nbins)
                continuum = self.emission_line_model.fill_continuum_to_match(binid)
                specres_ext='SPECRES' if self.binned_spectra.method['spec_res'] == 'cube' else None
                sres = self.binned_spectra.drpf.spectral_resolution(ext=specres_ext, toarray=True,
                                    pre=self.binned_spectra.method['prepixel_sres'], fill=True)
                sres = sres.data[good_snr,:]    # spectral_resolution returns a masked array
            else:
                continuum = self.emission_line_model.fill_continuum_to_match(
                                                self.binned_spectra['BINID'].data,
                                                missing=self.binned_spectra.missing_bins)
                sres = binned_spectra['SPECRES'].data.copy()[good_snr,:]
        elif self.stellar_continuum is not None:
            continuum = self.stellar_continuum.fill_to_match(
                                self.binned_spectra['BINID'].data,
                                missing=self.binned_spectra.missing_bins)
            sres = binned_spectra['SPECRES'].data.copy()[good_snr,:]
        _redshift = self.redshift[good_snr]

        #---------------------------------------------------------------
        # Perform the moment measurements
        measurements = self.measure_moments(self.momdb, wave, flux, ivar=ivar, sres=sres,
                                            continuum=continuum, redshift=_redshift,
                                            bitmask=self.bitmask)

        #---------------------------------------------------------------
        # Perform the equivalent width measurements
        include_band = numpy.array([numpy.invert(self.momdb.dummy)]*self.nbins) \
                        & numpy.invert(self.bitmask.flagged(measurements['MASK'],
                                                            flag=['BLUE_EMPTY', 'RED_EMPTY']))

        # Set the line center at the center of the primary passband
#        line_center = (1.0+_redshift)[:,None]*numpy.mean(self.momdb['primary'], axis=1)[None,:]
        line_center = (1.0+_redshift)[:,None]*self.momdb['restwave'][None,:]
#        line_center = (1.0+measurements['MOM1']/astropy.constants.c.to('km/s').value) \
#                            * self.momdb['restwave'][None,:]

#        pyplot.scatter(line_center_1, line_center_2, marker='.', s=50, color='k', lw=0, zorder=2)
#        pyplot.plot([0,10000], [0,10000], color='r', zorder=1)
#        pyplot.show()

        measurements['BMED'], measurements['RMED'], pos, measurements['EWCONT'], \
                measurements['EW'], measurements['EWERR'] \
                        = emission_line_equivalent_width(wave, flux, self.momdb['blueside'],
                                                         self.momdb['redside'], line_center,
                                                         measurements['FLUX'], redshift=_redshift,
                                                         line_flux_err=measurements['FLUXERR'],
                                                         include_band=include_band)
        measurements['REDSHIFT'] = _redshift

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

#        # Set the undefined bands
#        measurements['MASK'][:,self.momdb.dummy] \
#                    = self.bitmask.turn_on(measurements['MASK'][:,self.momdb.dummy],
#                                           'UNDEFINED_BANDS')

        #---------------------------------------------------------------
        # Initialize the header
        self.hardcopy = hardcopy
        pri_hdr = self._initialize_primary_header(measurements_binid=measurements_binid)
        map_hdr = DAPFitsUtil.build_map_header(self.binned_spectra.drpf,
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

#        # Add any spaxel not used because it was flagged by the binning
#        # step
#        indx = self.binned_spectra['MAPMASK'].data > 0
#        map_mask[indx] = self.bitmask.turn_on(map_mask[indx], 'DIDNOTUSE')
#        # Isolate any spaxels with foreground stars
#        indx = self.binned_spectra.bitmask.flagged(self.binned_spectra['MAPMASK'].data, 'FORESTAR')
#        map_mask[indx] = self.bitmask.turn_on(map_mask[indx], 'FORESTAR')
#        # Get the bins that were below the S/N limit
#        indx = numpy.invert(DAPFitsUtil.reconstruct_map(self.spatial_shape,
#                                                        self.binned_spectra['BINID'].data.ravel(),
#                                                        good_snr, dtype='bool')) & (map_mask == 0)
#        map_mask[indx] = self.bitmask.turn_on(map_mask[indx], 'LOW_SNR')
#
#        # Get the bin ids with measurements
#        bin_indx = DAPFitsUtil.downselect_bins(self.binned_spectra['BINID'].data.ravel(),
#                                               measurements['BINID'])

        # Compile the saved information on the suite of measured indices
        passband_database = self._compile_database()

        # Save the data to the hdu attribute
        self.hdu = fits.HDUList([ fits.PrimaryHDU(header=pri_hdr),
                                  fits.ImageHDU(data=bin_indx, header=map_hdr, name='BINID'),
                                  fits.ImageHDU(data=map_mask, header=map_hdr, name='MAPMASK'),
                                  fits.BinTableHDU.from_columns(
                                        [ fits.Column(name=n,
                                                      format=rec_to_fits_type(passband_database[n]),
                                                      array=passband_database[n])
                                          for n in passband_database.dtype.names ], name='ELMBAND'),
                                  fits.BinTableHDU.from_columns(
                                        [ fits.Column(name=n,
                                                      format=rec_to_fits_type(measurements[n]),
                                                      array=measurements[n])
                                          for n in measurements.dtype.names ], name='ELMMNTS'),
                                ])

        #---------------------------------------------------------------
        # Write the data, if requested
        if self.hardcopy:
            if not os.path.isdir(self.directory_path):
                os.makedirs(self.directory_path)
            self.write(clobber=clobber)
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, '-'*50)


    def write(self, clobber=False):
        """
        Write the hdu object to the file.
        """
        DAPFitsUtil.write(self.hdu, self.file_path(), clobber=clobber, checksum=True,
                          loggers=self.loggers, quiet=self.quiet)


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
        if self.hdu['PRIMARY'].header['ELMKEY'] != self.database['key']:
            if strict:
                raise ValueError('Keywords in header do not match specified method keywords!')
            else:
                warnings.warn('Keywords in header do not match specified method keywords!')
        # TODO: "strict" should also check other aspects of the file to
        # make sure that the details of the method are also the same,
        # not just the keyword

#        if not self.quiet:
#            log_output(self.loggers, 1, logging.INFO, 'Reverting to python-native structure.')
#        DAPFitsUtil.restructure_map(self.hdu, ext=self.image_arrays)

        self.nbins = self.hdu['PRIMARY'].header['NBINS']
        unique_bins = numpy.unique(self.hdu['BINID'].data.ravel()) \
                            if self.hdu['PRIMARY'].header['ELMREBIN'] else None
        self.missing_bins = self._get_missing_bins(unique_bins=unique_bins)


    # TODO: Use the EmissionMomentsDB.channel_names function!
    def channel_names(self):
        return numpy.array([ '{0}-{1}'.format(n,int(w)) \
                        for n,w in zip(self['ELMBAND'].data['NAME'],
                                       self['ELMBAND'].data['RESTWAVE']) ])

