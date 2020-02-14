# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Class that performs a number of assessments of a DRP file needed for
handling of the data by the DAP.  These assessments need only be done
once per DRP data file.

Revision history
----------------

    | **24 Mar 2016**: Implementation begun by K. Westfall (KBW)
    | **11 May 2016**: (KBW) Switch to using
        `pydl.goddard.astro.airtovac`_ instead of internal function
    | **19 May 2016**: (KBW) Added loggers and quiet keyword arguments
        to :class:`ReductionAssessment`, removed verbose 
    | **23 Feb 2017**: (KBW) Use DAPFitsUtil read and write functions.
    | **27 Feb 2017**: (KBW) Use DefaultConfig.
    | **17 May 2017**: (KBW) Added ability to use a response function
        for the flux statistics.

----

.. include license and copyright
.. include:: ../copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../links.rst
"""
import os
import time
import glob
import warnings
import logging

import numpy

from scipy import sparse
from pydl.goddard.astro import airtovac
from astropy.io import fits

from ..drpfits import DRPFits
from ..par.parset import KeywordParSet
from ..config.defaults import dap_source_dir, default_dap_common_path
from ..config.defaults import default_dap_file_name
from ..util.fitsutil import DAPFitsUtil
from ..util.covariance import Covariance
from ..util.geometry import SemiMajorAxisCoo
from ..util.fileio import init_record_array, rec_to_fits_type, create_symlink
from ..util.log import log_output
from ..util.parser import DefaultConfig
from .util import select_proc_method

# Used for testing/debugging
from matplotlib import pyplot
#from memory_profiler import profile

# Add strict versioning
# from distutils.version import StrictVersion

class ReductionAssessmentDef(KeywordParSet):
    """
    Class with parameters used to define how the reduction assessments
    are performed.  At the moment this is just a set of parameters that
    define how the S/N is calculated.

    See :class:`mangadap.par.parset.ParSet` for attributes.

    .. todo::

        - Allow for different ways of calculating covariance?

    The defined parameters are:

    .. include:: ../tables/reductionassessmentdef.rst

    """
    def __init__(self, key=None, waverange=None, response_func=None, covariance=False):
        # Perform some checks of the input
        ar_like = [ numpy.ndarray, list ]
        #covar_opt = covariance_options()
        
        pars =   [ 'key', 'waverange', 'response_func', 'covariance' ]
        values = [   key,   waverange,   response_func,   covariance ]
        #options = [ None,        None,    covar_opt ]
        dtypes = [   str,     ar_like,         ar_like,         bool ]
        descr = ['Keyword to distinguish the assessment method.',
                 'A two-element vector with the starting and ending wavelength (angstroms ' \
                    'in **vacuum**) within which to calculate the signal-to-noise',
                 'A two-column array with a response function to use for the S/N calculation.  ' \
                    'The columns must br the wavelength and amplitude of the response function, ' \
                    'respectively.',
                 'Type of covariance measurement to produce.']

        super(ReductionAssessmentDef, self).__init__(pars, values=values, dtypes=dtypes,
                                                     descr=descr)


def validate_reduction_assessment_config(cnfg):
    """ 
    Validate the :class:`mangadap.util.parser.DefaultConfig` object that
    provides the reduction-assessment method parameters.

    Args:
        cnfg (:class:`mangadap.util.parser.DefaultConfig`): Object with
            the reduction-assessment method parameters needed by
            :class:`mangadap.proc.reductionassessments.ReductionAssessmentDef`.

    Returns:
        bool: Booleans that specify how the reduction assessment should
        be constructed.  The flags specify to use (1) the wavelength
        range, (2) a bandpass filter parameter file, or (3) a file with
        a filter response function.

    Raises:
        KeyError: Raised if required keyword does not exist.
        ValueError: Raised if keys have unacceptable values.
        FileNotFoundError: Raised if a file is specified but could not
            be found.
    """
    # Check for required keywords
    required_keywords = ['key']
    if not cnfg.all_required(required_keywords):
        raise KeyError('Keywords {0} must all have valid values.'.format(required_keywords))

    def_range = cnfg.keyword_specified('wave_limits')
    def_response = cnfg.keyword_specified('response_function_file')

    if numpy.sum([ def_range, def_response] ) != 1:
        raise ValueError('Method undefined.  Must provide either \'waverange\' or '
                         '\'response_function_file\'.')

    if def_response and not os.path.isfile(cnfg['response_function_file']):
        raise FileNotFoundError('response_function_file does not exist: {0}'.format(
                                cnfg['response_function_file']))

    return def_range, def_response


def available_reduction_assessments(dapsrc=None):
    r"""
    Return the list of available reduction assessment methods.  To get a
    list of default methods provided by the DAP do::

        from mangadap.proc.reductionassessments import available_reduction_assessments
        rdx_methods = available_reduction_assessments()
        print(rdx_methods)

    Each element in the `rdx_methods` list is an instance of
    :class:`ReductionAssessmentDef`, which  is printed using the
    :class:`ParSet` base class representation function.

    New methods can be included by adding ini config files to
    `$MANGADAP_DIR/python/mangadap/config/reduction_assessments`.  See
    an example file at
    `$MANGADAP_DIR/python/mangadap/config/example_ini/example_reduction_assessment_config.ini`.

    Args:
        dapsrc (str): (**Optional**) Root path to the DAP source
            directory (i.e., $MANGADAP_DIR).  If not provided, the
            default is defined by
            :func:`mangadap.config.defaults.dap_source_dir`.

    Returns:
        list: A list of :func:`ReductionAssessmentDef` objects, each
        defining a separate assessment method.

    Raises:
        NotADirectoryError: Raised if the provided or default
            *dapsrc* is not a directory.
        OSError/IOError: Raised if no reduction assessment configuration
            files could be found.
        KeyError: Raised if the assessment method keywords are not all
            unique.
        NameError: Raised if either ConfigParser or
            ExtendedInterpolation are not correctly imported.  The
            latter is a *Python 3 only module*!
    """
    # Check the source directory exists
    dapsrc = dap_source_dir() if dapsrc is None else str(dapsrc)
    if not os.path.isdir(dapsrc):
        raise NotADirectoryError('{0} does not exist!'.format(dapsrc))

    # Check the configuration files exist
    ini_files = glob.glob(dapsrc+'/python/mangadap/config/reduction_assessments/*.ini')
    if len(ini_files) == 0:
        raise IOError('Could not find any configuration files in {0} !'.format(
                      dapsrc+'/python/mangadap/config/reduction_assessments'))

    # Build the list of library definitions
    assessment_methods = []
    for f in ini_files:
        # Read the config file
        cnfg = DefaultConfig(f=f, interpolate=True)
        # Ensure it has the necessary elements to define the template
        # library
#        def_range, def_par, def_response = validate_reduction_assessment_config(cnfg)
        def_range, def_response = validate_reduction_assessment_config(cnfg)
        in_vacuum = cnfg.getbool('in_vacuum', default=False)
        if def_range:
            waverange = cnfg.getlist('wave_limits', evaluate=True)
            if not in_vacuum:
                waverange = airtovac(waverange)
            assessment_methods += [ ReductionAssessmentDef(key=cnfg['key'], waverange=waverange,
                                                           covariance=cnfg.getbool('covariance',
                                                                            default=False)) ]
        elif def_response:
            response = numpy.genfromtxt(cnfg['response_function_file'])[:,:2]
            if not in_vacuum:
                response[:,0] = airtovac(response[:,0])
            assessment_methods += [ ReductionAssessmentDef(key=cnfg['key'], response_func=response,
                                                           covariance=cnfg.getbool('covariance',
                                                                                   default=False)) ]
        else:
            raise ValueError('Must define a wavelength range or a response function.')

    # Check the keywords of the libraries are all unique
    if len(numpy.unique( numpy.array([ method['key'] for method in assessment_methods ]) )) \
            != len(assessment_methods):
        raise KeyError('Reduction assessment method keywords are not all unique!')

    # Return the default list of assessment methods
    return assessment_methods


class ReductionAssessment:
    r"""
    Object used to perform and store assessments of a DRP file needed
    for handling of the data by the DAP.  These assessments need only be
    done once per DRP data file.

    See :func:`compute` for the provided data.

    Args:
        method_key (str): Keyword selecting the assessment method to
            use.
        drpf (:class:`mangadap.drpfits.DRPFits`): DRP file (object) to
            use for the assessments.
        pa (float): (**Optional**) On-sky position angle of the major
            axis used to calculate elliptical, semi-major-axis
            coordinates, defined as the angle from North through East
            and denoted :math:`\phi_0`.  Default is 0.0.
        ell (float): (**Optional**) Ellipticity defined as
            :math:`\varepsilon=1-b/a`, based on the semi-minor to
            semi-major axis ratio (:math:`b/a`) of the isophotal ellipse
            used to calculate elliptical, semi-major-axis coordinates.
            Default is 0.0.
        method_list (list): (**Optional**) List of
            :class:`ReductionAssessmentDef` objects that define the
            parameters required to perform the assessments of the DRP
            file.  The *method_key* must select one of these objects.
            Default is to construct list using
            :func:`available_reduction_assessments`.
        dapver (str): (**Optional**) DAP version, which is used to
            define the default DAP analysis path.  Default is defined by
            :func:`mangadap.config.defaults.default_dap_version`
        analysis_path (str): (**Optional**) The path to the top level
            directory containing the DAP output files for a given DRP
            and DAP version.  Default is defined by
            :func:`mangadap.config.defaults.default_analysis_path`.
        directory_path (str): (**Optional**) The exact path for the
            output file.  Default is defined by
            :func:`mangadap.config.defaults.default_dap_common_path`.
        output_file (str): (**Optional**) The name of the file for the
            computed assessments.  The full path of the output file will
            be :attr:`directory_path`/:attr:`output_file`.  Default is
            defined by
            :func:`mangadap.config.defaults.default_dap_file_name`.
        hardcopy (bool): (**Optional**) Flag to write the data to a fits
            file.  Default is True.
        symlink_dir (str): (**Optional**) Create a symlink to the file
            in this directory.  Default is for no symlink.
        clobber (bool): (**Optional**) If the output file already
            exists, this will force the assessments to be redone and the
            output file to be overwritten.  Default is False.
        checksum (bool): (**Optional**) Compare the checksum when
            reading the data.
        loggers (list): (**Optional**) List of `logging.Logger`_ objects
            to log progress; ignored if quiet=True.  Logging is done
            using :func:`mangadap.util.log.log_output`.  Default is no
            logging.
        quiet (bool): (**Optional**) Suppress all terminal and logging
            output.  Default is False.

    Attributes:
        method (:class:`ReductionAssessmentDef`):  Parameters defining the
            method to use for the reduction assessments.
        drpf (:class:`mangadap.drpfits.DRPFits`): DRP file (object) with
            which the template library is associated for analysis
        pa (float): On-sky position angle of the major axis used to
            calculate elliptical, semi-major-axis coordinates, defined
            as the angle from North through East and denoted
            :math:`\phi_0`.
        ell (float): Ellipticity defined as :math:`\varepsilon=1-b/a`,
            based on the semi-minor to semi-major axis ratio
            (:math:`b/a`) of the isophotal ellipse used to calculate
            elliptical, semi-major-axis coordinates.
        directory_path (str): The exact path for the output file.
            Default is defined by
            :func:`mangadap.config.defaults.default_dap_common_path`.
        output_file (str): The name of the file for the computed
            assessments.  The full path of the output file will be
            :attr:`directory_path`/:attr:`output_file`.  Default is
            defined by :func:`_set_paths`.
        hardcopy (bool): Flag to keep a hardcopy of the data by writing
            the data to a fits file.
        symlink_dir (str): Symlink created to the file in this directory
        hdu (`astropy.io.fits.hdu.hdulist.HDUList`_): HDUList with the
            data, with columns as described above.
        correlation (:class:`mangadap.util.covariance.Covariance`):
            Covariance matrix for the mean flux measurements, if
            calculated.
        loggers (list): List of `logging.Logger`_ objects to log
            progress; ignored if quiet=True.  Logging is done using
            :func:`mangadap.util.log.log_output`.
        quiet (bool): Suppress all terminal and logging output.

    """
#    @profile
    def __init__(self, method_key, drpf, pa=0.0, ell=0.0, method_list=None, dapsrc=None,
                 dapver=None, analysis_path=None, directory_path=None, output_file=None,
                 hardcopy=True, symlink_dir=None, clobber=False, checksum=False, loggers=None,
                 quiet=False):
                 
        self.loggers = None
        self.quiet = False

        # Define the method properties
        self.method = None
        self._define_method(method_key, method_list=method_list, dapsrc=dapsrc)

        # Set in compute via _set_paths
        self.drpf = None

        # Define the output directory and file
        self.directory_path = None      # Set in _set_paths
        self.output_file = None
        self.hardcopy = None
        self.symlink_dir = None

        # Initialize the objects used in the assessments
        self.pa = None
        self.ell = None
        self.hdu = None
        self.checksum = checksum
        self.correlation = None
        self.covar_wave = None
        self.covar_channel = None

        # Run the assessments of the DRP file
        self.compute(drpf, pa=pa, ell=ell, dapver=dapver, analysis_path=analysis_path,
                     directory_path=directory_path, output_file=output_file, hardcopy=hardcopy,
                     symlink_dir=symlink_dir, clobber=clobber, loggers=loggers, quiet=quiet)


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


    def _define_method(self, method_key, method_list=None, dapsrc=None):
        """
        Select the assessment method from the provided list.  Used to
        set :attr:`method`; see
        :func:`mangadap.proc.util.select_proc_method`.

        Args:
            method_key (str): Keyword of the selected method.  Available
                methods are provided by
                :func:`available_reduction_assessments`
            method_list (list): (**Optional**) List of
                :class:`ReductionAssessmentDef` objects that define the
                parameters required to assess the reduced data.
            dapsrc (str): (**Optional**) Root path to the DAP source
                directory.  If not provided, the default is defined by
                :func:`mangadap.config.defaults.dap_source_dir`.
        """
        self.method = select_proc_method(method_key, ReductionAssessmentDef,
                                         method_list=method_list,
                                         available_func=available_reduction_assessments)

    def _set_paths(self, directory_path, dapver, analysis_path, output_file):
        """
        Set the I/O paths.  Used to set :attr:`directory_path` and
        :attr:`output_file`.  If not provided, the defaults are set
        using, respectively,
        :func:`mangadap.config.defaults.default_dap_common_path` and
        :func:`mangadap.config.defaults.default_dap_file_name`.

        Args:
            directory_path (str): The exact path to the DAP reduction
                assessments file.  See :attr:`directory_path`.
            dapver (str): DAP version.
            analysis_path (str): The path to the top-level directory
                containing the DAP output files for a given DRP and DAP
                version.
            output_file (str): The name of the file with the reduction
                assessments.  See :func:`compute`.
        """
        self.directory_path, self.output_file \
                = ReductionAssessment.default_paths(self.drpf.plate, self.drpf.ifudesign,
                                                    self.method['key'],
                                                    directory_path=directory_path,
                                                    drpver=self.drpf.drpver, dapver=dapver,
                                                    analysis_path=analysis_path,
                                                    output_file=output_file)


    def _per_spectrum_dtype(self):
        r"""
        Construct the record array data type for the output fits
        extension.
        """
        return [ ('DRP_INDEX',numpy.int,(numpy.asarray(tuple(self.drpf.spatial_index)).shape[1],)),
                 ('SKY_COO',numpy.float,(2,)),
                 ('ELL_COO',numpy.float,(2,)),
                 ('FGOODPIX',numpy.float),
                 ('MINEQMAX',numpy.uint8),
                 ('SIGNAL',numpy.float),
                 ('VARIANCE',numpy.float),
                 ('SNR',numpy.float)
               ]

    
    def _initialize_primary_header(self, hdr=None):
        # Copy the from the DRP and clean it
        if hdr is None:
            hdr = self.drpf.hdu['PRIMARY'].header.copy()
            hdr = DAPFitsUtil.clean_dap_primary_header(hdr)

        hdr['AUTHOR'] = 'Kyle B. Westfall <westfall@ucolick.org>'
        hdr['RDXQAKEY'] = (self.method['key'], 'Method keyword')
        hdr['ECOOPA'] = (self.pa, 'Position angle for ellip. coo')
        hdr['ECOOELL'] = (self.ell, 'Ellipticity (1-b/a) for ellip. coo')
        if self.method['covariance']:
            hdr['BBWAVE'] = (self.covar_wave, 'Covariance channel wavelength')
            hdr['BBINDEX'] = (self.covar_channel, 'Covariance channel index')
        return hdr


    @staticmethod
    def default_paths(plate, ifudesign, method_key, directory_path=None, drpver=None, dapver=None,
                      analysis_path=None, output_file=None):
        """
        Set the default directory and file name for the output file.

        Args:
            plate (int): Plate number of the observation.
            ifudesign (int): IFU design number of the observation.
            method_key (str): Keyword designating the method used for
                the reduction assessments.
            directory_path (str): (**Optional**) The exact path to the
                DAP reduction assessments file.  Default set by
                :func:`mangadap.config.defaults.default_dap_common_path`.
            drpver (str): (**Optional**) DRP version.  Default set by
                :func:`mangadap.config.defaults.default_drp_version`.
            dapver (str): (**Optional**) DAP version.  Default set by
                :func:`mangadap.config.defaults.default_dap_version`.
            analysis_path (str): (**Optional**) The path to the
                top-level directory containing the DAP output files for
                a given DRP and DAP version. Default set by
                :func:`mangadap.config.defaults.default_analysis_path`.
            output_file (str): (**Optional**) The name of the file with
                the reduction assessments.  Default set by
                :func:`mangadap.config.defaults.default_dap_file_name`.

        Returns:
            str: Two strings with the path for the output file and the
            name of the output file.
        """
        _directory_path = default_dap_common_path(plate=plate, ifudesign=ifudesign,
                                                  drpver=drpver, dapver=dapver,
                                                  analysis_path=analysis_path) \
                                        if directory_path is None else str(directory_path)
        _output_file = default_dap_file_name(plate, ifudesign, method_key) \
                                        if output_file is None else str(output_file)
        return _directory_path, _output_file


    def file_name(self):
        """Return the name of the output file."""
        return self.output_file


    def file_path(self):
        """Return the full path to the output file."""
        if self.directory_path is None or self.output_file is None:
            return None
        return os.path.join(self.directory_path, self.output_file)


    def info(self):
        return self.hdu.info()


    def compute(self, drpf, pa=None, ell=None, dapver=None, analysis_path=None, directory_path=None,
                output_file=None, hardcopy=True, symlink_dir=None, clobber=False, loggers=None,
                quiet=False):
        r"""
        Compute and output the main data products.  The list of HDUs
        are:

            - ``PRIMARY`` : Empty apart from the header information.
            - ``SPECTRUM`` : Extension with the main, per-spectrum
              measurements; see below.
            - ``CORREL`` : The correlation matrix between the ``SIGNAL``
              measurements provided in the ``SPECTRUM`` extension.  The
              format of this extension is identical to the nominal
              output of the :class:`mangadap.util.covariance.Covariance`
              object; see
              :func:`mangadap.util.covariance.Covariance.write`.

        The ``SPECTRUM`` extension contains the following columns:

            - ``DRP_INDEX`` : Array coordinates in the DRP file; see
              :attr:`mangadap.drpfits.DRPFits.spatial_index`.  For RSS
              files, this is a single integer; for CUBE files, it is a
              vector of two integers.
            - ``SKY_COO`` : On-sky X and Y coordinates.  Coordinates are
              sky-right offsets from the object center; i.e., positive X
              is along the direction of positive right ascension.  See
              :class:`mangadap.drpfits.DRPFits.mean_sky_coordinates`.
            - ``ELL_COO`` : Elliptical (semi-major axis) radius and
              azimuth angle from N through East with respect to the
              photometric position angle; based on the provided
              ellipticity parameters.  See
              :class:`mangadap.util.geometry.SemiMajorAxisCoo`.
            - ``FGOODPIX`` : Fraction of good pixels in each spectrum.
            - ``MINEQMAX`` : Flag that min(flux) = max(flux) in the
              spectrum; i.e., the spaxel has no data.
            - ``SIGNAL``, ``VARIANCE``, ``SNR`` : Per pixel means of the
              flux, flux variance, and signal-to-noise.  The
              ``VARIANCE`` and ``SNR`` columns use the inverse variance
              provided by the DRP.  See
              :func:`mangadap.drpfits.DRPFits.flux_stats`.

        Args:
            drpf (:class:`mangadap.drpfits.DRPFits`): DRP file (object)
                to use for the assessments.
            pa (float): (**Optional**) On-sky position angle of the
                major axis used to calculate elliptical, semi-major-axis
                coordinates, defined as the angle from North through
                East and denoted :math:`\phi_0`.  Default is 0.0.
            ell (float): (**Optional**) Ellipticity defined as
                :math:`\varepsilon=1-b/a`, based on the semi-minor to
                semi-major axis ratio (:math:`b/a`) of the isophotal
                ellipse used to calculate elliptical, semi-major-axis
                coordinates.  Default is 0.0.
            dapver (str): (**Optional**) DAP version, which is used to
                define the default DAP analysis path.  Default is
                defined by
                :func:`mangadap.config.defaults.default_dap_version`
            analysis_path (str): (**Optional**) The path to the top
                level directory containing the DAP output files for a
                given DRP and DAP version.  Default is defined by
                :func:`mangadap.config.defaults.default_analysis_path`.
            directory_path (str): (**Optional**) The exact path for the
                output file.  Default is defined by
                :func:`mangadap.config.defaults.default_dap_common_path`.
            output_file (str): (**Optional**) The name of the file for
                the computed assessments.  The full path of the output
                file will be :attr:`directory_path`/:attr:`output_file`.
                Default is defined by
                :func:`mangadap.config.defaults.default_reduction_assessments_file`.
            hardcopy (bool): (**Optional**) Flag to write the data to a
                fits file.  Default is True.
            symlink_dir (str): (**Optional**) Create a symlink to the
                file in this directory.  Default is for no symlink.
            clobber (bool): (**Optional**) If the output file already
                exists, this will force the assessments to be redone and
                the output file to be overwritten.  Default is False.
            loggers (list): (**Optional**) List of `logging.Logger`_
                objects to log progress; ignored if quiet=True.  Default
                is no logging.
            quiet (bool): (**Optional**) Suppress all terminal and
                logging output.  Default is False.

        Raises:
            ValueError: Raise if no DRPFits object is provided or if
                the output file is undefined.
        """

        # Initialize the reporting
        if loggers is not None:
            self.loggers = loggers
        self.quiet = quiet

        # Check the input DRPFits object
        if drpf is None:
            raise ValueError('Must provide DRP file object to compute assessments.')
        if not isinstance(drpf, DRPFits):
            raise TypeError('Must provide a valid DRPFits object!')
        if drpf.hdu is None:
            if not self.quiet:
                warnings.warn('DRP file previously unopened.  Reading now.')
            drpf.open_hdu()
        self.drpf = drpf

        # Test if the RSS file exists; cannot compute covariance if not
        if self.method['covariance'] and not drpf.can_compute_covariance:
            warnings.warn('RSS counterpart not available.  Cannot determine covariance matrix!')
            self.method['covariance'] = False

        # Reset the output paths if necessary
        self._set_paths(directory_path, dapver, analysis_path, output_file)
        # Check that the path for or to the file is defined
        ofile = self.file_path()
        if ofile is None:
            raise ValueError('File path for output file is undefined!')

        # Report
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, '-'*50)
            log_output(self.loggers, 1, logging.INFO, '{0:^50}'.format('REDUCTION ASSESSMENT'
                                                                       ' COMPUTATIONS'))
            log_output(self.loggers, 1, logging.INFO, '-'*50)
            log_output(self.loggers, 1, logging.INFO, 'Output path: {0}'.format(
                                                                            self.directory_path))
            log_output(self.loggers, 1, logging.INFO, 'Output file: {0}'.format(
                                                                            self.output_file))

        # If the file already exists, and not clobbering, just read the
        # file
        self.symlink_dir = symlink_dir
        if os.path.isfile(ofile) and not clobber:
            self.hardcopy = True
            if not self.quiet:
                log_output(self.loggers, 1, logging.INFO, 'Reading exiting file')

            # Read the existing output file
            self.hdu = DAPFitsUtil.read(ofile, checksum=self.checksum)

            # Construct the correlation matrix with the appropriate
            # variance
            if self.method['covariance']:
                self.correlation = Covariance.from_fits(self.hdu, ivar_ext=None, row_major=True,
                                                        correlation=True)
                self.correlation = self.correlation.apply_new_variance(
                                                            self.hdu['SPECTRUM'].data['VARIANCE'])
            else:
                self.correlation = None

            # Read the header data
            self.pa = self.hdu['PRIMARY'].header['ECOOPA']
            if not self.quiet and pa is not None and self.pa != pa:
                warnings.warn('Provided position angle different from available file; set ' \
                              'clobber=True to overwrite.')
            self.ell = self.hdu['PRIMARY'].header['ECOOELL']
            if not self.quiet and ell is not None and self.ell != ell:
                warnings.warn('Provided ellipticity different from available file; set ' \
                              'clobber=True to overwrite.')
            if self.method['covariance']:
                self.covar_wave = self.hdu['PRIMARY'].header['BBWAVE']
                if isinstance(self.covar_wave, str):
                    self.covar_wave = eval(self.covar_wave)
                self.covar_channel = self.hdu['PRIMARY'].header['BBINDEX']
                if isinstance(self.covar_channel, str):
                    self.covar_channel = eval(self.covar_channel)
            else:
                self.covar_wave = None
                self.covar_channel = None
#            print(self.covar_wave, self.covar_channel)

            # Make sure the symlink exists
            if self.symlink_dir is not None:
                create_symlink(ofile, self.symlink_dir, clobber=clobber)
            if not self.quiet:
                log_output(self.loggers, 1, logging.INFO, '-'*50)
            return

        # (Re)Initialize some of the attributes
        if pa is not None:
            self.pa = pa
        if ell is not None:
            self.ell = ell

        # Only mask pixels with the following flags set
        flags = ['DONOTUSE', 'FORESTAR']

        # Initialize the record array for the SPECTRUM extension
        spectrum_data = init_record_array(drpf.nspec, self._per_spectrum_dtype())
        spectrum_data['DRP_INDEX'] = numpy.asarray(tuple(drpf.spatial_index))
        spectrum_data['SKY_COO'][:,0], spectrum_data['SKY_COO'][:,1] \
                = drpf.mean_sky_coordinates(waverange=self.method['waverange'],
                                            response_func=self.method['response_func'],
                                            offset=True, flag=flags)

#        print(spectrum_data['DRP_INDEX'].shape)
#        print(spectrum_data['DRP_INDEX'][0:drpf.nx,0])
#        print(spectrum_data['DRP_INDEX'][0:drpf.nx,1])
#        print(spectrum_data['SKY_COO'].shape)
#        print(spectrum_data['SKY_COO'][0:drpf.nx,0])
#        print(spectrum_data['SKY_COO'][0:drpf.nx,1])

#        pyplot.scatter(spectrum_data['SKY_COO'][:,0], spectrum_data['SKY_COO'][:,1], marker='.',
#                       s=30, color='k', lw=0)
#        pyplot.show()

        coord_conversion = SemiMajorAxisCoo(xc=0, yc=0, rot=0, pa=self.pa, ell=self.ell)
        spectrum_data['ELL_COO'][:,0], spectrum_data['ELL_COO'][:,1] \
            = coord_conversion.polar(spectrum_data['SKY_COO'][:,0], spectrum_data['SKY_COO'][:,1])

#        pyplot.scatter(spectrum_data['ELL_COO'][:,0], spectrum_data['ELL_COO'][:,1], marker='.',
#                       s=30, color='k', lw=0)
#        pyplot.show()
       
        flux = drpf.copy_to_masked_array(flag=flags)
        spectrum_data['FGOODPIX'] = numpy.sum(numpy.invert(numpy.ma.getmaskarray(flux)),axis=1) \
                                            / flux.shape[1]

        # TODO: Is this superfluous now?
        frange = numpy.ma.max(flux, axis=1)-numpy.ma.min(flux, axis=1)
        spectrum_data['MINEQMAX'] = (numpy.invert(numpy.ma.getmaskarray(frange))) \
                                        & (numpy.ma.absolute(frange) < 1e-10)
#        print(spectrum_data['MINEQMAX'])
#        print(drpf.nspec, numpy.sum(spectrum_data['MINEQMAX'] & (spectrum_data['FGOODPIX'] > 0.8)))

#        srt = numpy.argsort(spectrum_data['FGOODPIX'])
#        grw = numpy.arange(len(srt))/len(srt)
#        pyplot.step(spectrum_data['FGOODPIX'][srt], grw, color='k')
#        pyplot.show()

        # Get the wavelength range
        waverange = [ drpf['WAVE'].data[0], drpf['WAVE'].data[-1] ] \
                    if self.method['waverange'] is None and self.method['response_func'] is None \
                    else (self.method['waverange'] if self.method['waverange'] is not None else
                        [ self.method['response_func'][0,0], self.method['response_func'][-1,0] ])
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO,
                       'Wavelength range for S and N calculation: {0:.1f} -- {1:.1f}'.format(
                                *waverange))

        # If calculating covariance, make sure that all the pixels in
        # the channel used and selected for analysis (at minimum based
        # on fraction of valid pixels in the spectrum - FGOODPIX) have
        # covariance measurements.  This is done by first assuming the
        # channel at the center of the wavelength range and then
        # flipping back and forth across this channel until a valid
        # channel is found or you reach the limits of the defined
        # wavelength range.  In the latter case, and exception is
        # raised!
        if self.method['covariance']:
            self.covar_wave = drpf._covariance_wavelength(waverange=self.method['waverange'],
                                                        response_func=self.method['response_func'],
                                                          flag=flags)
            self.covar_channel = numpy.argsort(numpy.absolute(drpf['WAVE'].data-self.covar_wave))[0]

            # TODO: Make this fraction an input parameter!
            goodfrac = (spectrum_data['FGOODPIX'].reshape(drpf.spatial_shape) > 0.8)

            off = 1
            sign = 1
            while self.covar_wave > waverange[0] and self.covar_wave < waverange[1] \
                    and numpy.sum(goodfrac & numpy.invert(drpf['IVAR'].data[:,:,self.covar_channel]
                                                          > 0)) > 0:
                self.covar_channel += sign*off
                self.covar_wave = drpf['WAVE'].data[self.covar_channel]
#                print(self.covar_wave, self.covar_channel, off, sign,
#                      numpy.sum(goodfrac & ~(drpf['IVAR'].data[:,:,self.covar_channel] > 0)))
                sign *= -1
                off += 1

            if numpy.sum(goodfrac &
                                numpy.invert(drpf['IVAR'].data[:,:,self.covar_channel] > 0)) > 0:
                raise ValueError('Unable to find wavelength channel within fully valid data.')
            
            if not self.quiet:
                log_output(self.loggers, 1, logging.INFO,
                           'Covariance measured at wavelength: {0} (channel {1})'.format(
                                self.covar_wave, self.covar_channel))
        else:
            self.covar_wave = None
            self.covar_channel = None

        # Determine the statistics of the flux in the specified
        # wavelength range
        spectrum_data['SIGNAL'], spectrum_data['VARIANCE'], spectrum_data['SNR'], self.correlation \
                = drpf.flux_stats(waverange=self.method['waverange'],
                                  response_func=self.method['response_func'],
                                  flag=['DONOTUSE', 'FORESTAR'],
                                  covar=self.method['covariance'], correlation=True,
                                  covar_wave=self.covar_wave)

        # Force the covariance matrix to use the mean variance instead
        # of the variance in the single, specified channel.
        if self.method['covariance']:
            self.correlation = self.correlation.apply_new_variance(spectrum_data['VARIANCE'])
#            self.correlation.show()

#        t = (spectrum_data['FGOODPIX'].reshape(drpf.spatial_shape) > 0.8) & \
#                        (drpf.bitmask.flagged(drpf['MASK'].data[:,:,self.covar_channel],
#                                             flag=['DONOTUSE', 'FORESTAR']))
#        print(numpy.sum(t))
#
#        t = (spectrum_data['FGOODPIX'].reshape(drpf.spatial_shape) > 0.8) & \
#                        ~(drpf['IVAR'].data[:,:,self.covar_channel] > 0)
#        print(numpy.sum(t))
#
#        t = (spectrum_data['FGOODPIX'].reshape(drpf.spatial_shape) > 0.8) & \
#                        ~(drpf['FLUX'].data[:,:,self.covar_channel] > 0)
#        print(numpy.sum(t))
#
#        goodpix_im = spectrum_data['FGOODPIX'].reshape(drpf.spatial_shape)
#        m = drpf['FLUX'].data[:,:,self.covar_channel]
#        pyplot.imshow(m, origin='lower', interpolation='nearest')
#        pyplot.colorbar()
#        pyplot.show()
#
#        maskedm = numpy.ma.MaskedArray(m, mask=~t)
#        pyplot.imshow(maskedm, origin='lower', interpolation='nearest')
#        pyplot.colorbar()
#        pyplot.show()

#        if correlation is not None:
#            correlation.show()
#
#        pyplot.scatter(spectrum_data['SIGNAL'], spectrum_data['VARIANCE'], marker='.', color='k',
#                       s=30)
#        pyplot.plot([0,2], numpy.square(numpy.array([0,2])/50.), color='r')
#        pyplot.show()

        # Construct header
        hdr = self._initialize_primary_header()

        # Get the covariance hdu
        if self.method['covariance']:
            hdr, ivar_hdu, covar_hdu = self.correlation.output_hdus(reshape=True, hdr=hdr)

        # Get the main extension columns and construct the HDUList
        spectrum_cols = [ fits.Column(name=n, format=rec_to_fits_type(spectrum_data[n]),
                                      array=spectrum_data[n]) for n in spectrum_data.dtype.names ]
        self.hdu = fits.HDUList([ fits.PrimaryHDU(header=hdr),
                                  fits.BinTableHDU.from_columns( spectrum_cols, name='SPECTRUM') ])

        # Add the covariance information
        if self.method['covariance']:
            self.hdu += [ covar_hdu ]

        # Write the file
        if not os.path.isdir(self.directory_path):
            os.makedirs(self.directory_path)
        self.hardcopy = hardcopy
        if self.hardcopy:
            DAPFitsUtil.write(self.hdu, ofile, clobber=clobber, checksum=True,
                              symlink_dir=self.symlink_dir, loggers=self.loggers)
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, '-'*50)




