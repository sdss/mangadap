# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Class that performs a number of assessments of a DRP file needed for
handling of the data by the DAP.  These assessments need only be done
once per DRP data file.

----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import os
import time
import glob
import warnings
import logging

from IPython import embed

import numpy

from scipy import sparse
from matplotlib import pyplot

from astropy.io import fits

from pydl.goddard.astro import airtovac

from ..datacube import DataCube
from ..par.parset import KeywordParSet

from ..config import defaults

from ..util.datatable import DataTable
from ..util.fitsutil import DAPFitsUtil
from ..util.covariance import Covariance
from ..util.geometry import SemiMajorAxisCoo
from ..util.fileio import create_symlink
from ..util.log import log_output
from ..util.parser import DefaultConfig
from .util import select_proc_method

# Add strict versioning
# from distutils.version import StrictVersion


# TODO: Add center coordinates
#            center_coo (:obj:`tuple`, optional):
#                A two-tuple with the coordinates in right-ascension
#                and declination for the coordinate-frame origin. If
#                None, no offset is performed.

# TODO: Fix the fact that if the correlation matrix is taken directly
# from the datacube that the wavelength and channel won't be defined.
class ReductionAssessmentDef(KeywordParSet):
    """
    Define how to assess the reduced data.

    At the moment this is just a set of parameters that define how
    the S/N is calculated.

    See :class:`mangadap.par.parset.ParSet` for attributes.

    .. todo::

        - Allow for different ways of calculating covariance?

    The defined parameters are:

    .. include:: ../tables/reductionassessmentdef.rst

    """
    def __init__(self, key=None, waverange=None, response_func=None, covariance=False,
                 minimum_frac=0.8):
        # Perform some checks of the input
        ar_like = [ numpy.ndarray, list ]
        in_fl = [ int, float ]
        
        pars =   ['key', 'waverange', 'response_func', 'covariance', 'minimum_frac']
        values = [key, waverange, response_func, covariance, minimum_frac]
        dtypes = [str, ar_like, ar_like, bool, in_fl]
        descr = ['Keyword to distinguish the assessment method.',
                 'A two-element vector with the starting and ending wavelength (angstroms ' \
                    'in **vacuum**) within which to calculate the signal-to-noise',
                 'A two-column array with a response function to use for the S/N calculation.  ' \
                    'The columns must br the wavelength and amplitude of the response function, ' \
                    'respectively.',
                 'Provide the fiducial spatial covariance.  If this is False, no spatial ' \
                     'covariance will be available for calculations that use the results of the ' \
                     'reduction assessments.  If True and a covariance matrix is available ' \
                     'directly from the datacube, it will be used.  If True and the datacube ' \
                     'does not already provide a computed covariance matrix, one is calculated ' \
                     'using :func:`mangadap.datacube.datacube.DataCube.covariance_matrix` ' \
                     'method.  **WARNING**: If the latter fails either, the DAP will issue a ' \
                     'warning and continue assuming no spatial covariance.',
                 'Minimum fraction of unmasked pixels in a spectrum required for inclusion in ' \
                     'the spatial covariance calculation.  Note this should match the value ' \
                     'used for the spatial-binning module.']

        super(ReductionAssessmentDef, self).__init__(pars, values=values, dtypes=dtypes,
                                                     descr=descr)


def validate_reduction_assessment_config(cnfg):
    """ 
    Validate the :class:`mangadap.util.parser.DefaultConfig` object that
    provides the reduction-assessment method parameters.

    Args:
        cnfg (:class:`mangadap.util.parser.DefaultConfig`):
            Object with the reduction-assessment method parameters
            needed by
            :class:`mangadap.proc.reductionassessments.ReductionAssessmentDef`.

    Returns:
        :obj:`bool`: Booleans that specify how the reduction
        assessment should be constructed. The flags specify to use
        (1) the wavelength range, (2) a bandpass filter parameter
        file, or (3) a file with a filter response function.

    Raises:
        KeyError:
            Raised if required keyword does not exist.
        ValueError:
            Raised if keys have unacceptable values.
        FileNotFoundError:
            Raised if a file is specified but could not be found.
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


def available_reduction_assessments():
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
    `$MANGADAP_DIR/mangadap/config/reduction_assessments`. See an
    example file at
    `$MANGADAP_DIR/example_ini/example_reduction_assessment_config.ini`.

    Returns:
        :obj:`list`: A list of :class:`ReductionAssessmentDef`
        objects, each defining a separate assessment method.

    Raises:
        IOError:
            Raised if no reduction assessment configuration files
            could be found.
        KeyError:
            Raised if the assessment method keywords are not all
            unique.
        ValueError:
            Raised if a wavelength range or a response function are
            not defined by any of the methods.
    """
    # Check the configuration files exist
    search_dir = os.path.join(defaults.dap_config_root(), 'reduction_assessments')
    ini_files = glob.glob(os.path.join(search_dir, '*.ini'))
    if len(ini_files) == 0:
        raise IOError('Could not find any configuration files in {0} !'.format(search_dir))

    # Build the list of library definitions
    assessment_methods = []
    for f in ini_files:
        # Read the config file
        cnfg = DefaultConfig(f=f, interpolate=True)
        # Ensure it has the necessary elements to define the template
        # library
        def_range, def_response = validate_reduction_assessment_config(cnfg)
        in_vacuum = cnfg.getbool('in_vacuum', default=False)
        if def_range:
            waverange = cnfg.getlist('wave_limits', evaluate=True)
            if not in_vacuum:
                waverange = airtovac(waverange)
            assessment_methods += [ReductionAssessmentDef(key=cnfg['key'], waverange=waverange,
                                        covariance=cnfg.getbool('covariance', default=False),
                                        minimum_frac=cnfg.getfloat('minimum_frac', default=0.8))]
        elif def_response:
            response = numpy.genfromtxt(cnfg['response_function_file'])[:,:2]
            if not in_vacuum:
                response[:,0] = airtovac(response[:,0])
            assessment_methods += [ReductionAssessmentDef(key=cnfg['key'], response_func=response,
                                        covariance=cnfg.getbool('covariance', default=False),
                                        minimum_frac=cnfg.getfloat('minimum_frac', default=0.8))]
        else:
            raise ValueError('Must define a wavelength range or a response function.')

    # Check the keywords of the libraries are all unique
    if len(numpy.unique( numpy.array([ method['key'] for method in assessment_methods ]) )) \
            != len(assessment_methods):
        raise KeyError('Reduction assessment method keywords are not all unique!')

    # Return the default list of assessment methods
    return assessment_methods


class ReductionAssessmentDataTable(DataTable):
    """
    Primary data table with the reduction assessments.

    The data table includes:

    .. include:: ../tables/reductionassessmentdatatable.rst

    Args:
        shape (:obj:`int`, :obj:`tuple`, optional):
            The shape of the initial array. If None, the data array
            will not be instantiated; use :func:`init` to initialize
            the data array after instantiation.
    """
    def __init__(self, shape=None):
        # NOTE: This should require python 3.7 to make sure that this
        # is an "ordered" dictionary.
        datamodel = dict(SPAT_INDX=dict(typ=int, shape=(2,),
                            descr='The 2D spatial array indices associated with each measurement. '
                                'See :attr:`~mangadap.datacube.datacube.DataCube.spatial_index`.'),
                         SKY_COO=dict(typ=float, shape=(2,),
                            descr='On-sky X and Y coordinates.  Coordinates are sky-right '
                                  'offsets from the object center; i.e., positive X is along the '
                                  'direction of positive right ascension. See '
                            ':class:`mangadap.datacube.datacube.DataCube.mean_sky_coordinates`.'),
                         ELL_COO=dict(typ=float, shape=(2,),
                                      descr='Elliptical (semi-major axis) radius and azimuth '
                                            'angle from N through E with respect to the '
                                            'photometric position angle; based on the provided '
                                            'ellipticity parameters.  See '
                                            ':class:`mangadap.util.geometry.SemiMajorAxisCoo`.'),
                         FGOODPIX=dict(typ=float, shape=None,
                                       descr='Fraction of good pixels in each spectrum.'),
#                         MINEQMAX=dict(typ=bool, shape=None,
#                                       descr='Flag that min(flux) = max(flux) in the spectrum; '
#                                             'i.e., the spaxel has no data.'),
                         SIGNAL=dict(typ=float, shape=None, descr='Per pixel mean flux'),
                         VARIANCE=dict(typ=float, shape=None, descr='Per pixel mean variance'),
                         SNR=dict(typ=float, shape=None, descr='Per pixel mean S/N'))

        keys = list(datamodel.keys())
        super(ReductionAssessmentDataTable,
                self).__init__(keys, [datamodel[k]['typ'] for k in keys],
                               element_shapes=[datamodel[k]['shape'] for k in keys],
                               descr=[datamodel[k]['descr'] for k in keys],
                               shape=shape)

class ReductionAssessment:
    r"""
    Perform basic assessments of the datacube.

    These assessments are needed for handling of the data by the DAP
    and only to be done once per datacube.

    See :func:`compute` for the provided data.

    Args:
        method_key (:obj:`str`):
            Keyword selecting the assessment method to use.
        cube (:class:`mangadap.datacube.datacube.DataCube`):
            Datacube to analyze.
        pa (:obj:`float`, optional):
            On-sky position angle of the major axis used to calculate
            elliptical, semi-major-axis coordinates, defined as the
            angle from North through East and denoted :math:`\phi_0`.
        ell (:obj:`float`, optional):
            Ellipticity defined as :math:`\varepsilon=1-b/a`, based
            on the semi-minor to semi-major axis ratio (:math:`b/a`)
            of the isophotal ellipse used to calculate elliptical,
            semi-major-axis coordinates.
        method_list (:obj:`list`, optional):
            List of :class:`ReductionAssessmentDef` objects that
            define the parameters required to perform the assessments
            of the DRP file. The ``method_key`` must select one of
            these objects. Default is to construct the list using
            :func:`available_reduction_assessments`.
        dapver (:obj:`str`, optional):
            DAP version, which is used to define the default DAP
            analysis path. (This does **not** define the version of
            the code to use.) Default is defined by
            :func:`mangadap.config.defaults.dap_version`
        analysis_path (:obj:`str`, optional):
            The path to the top-level directory for the DAP output
            files for a given DRP and DAP version. Default is defined
            by :func:`mangadap.config.defaults.dap_analysis_path`.
        directory_path (:obj:`str`, optional):
            The exact path for the output file. Default is defined by
            :func:`mangadap.config.defaults.dap_common_path`.
        output_file (:obj:`str`, optional):
            The name of the file for the computed assessments. The
            full path of the output file will be
            :attr:`directory_path`/:attr:`output_file`. Default is
            defined by
            :func:`mangadap.config.defaults.dap_file_name`.
        hardcopy (:obj:`bool`, optional):
            Flag to write the data to a fits file.
        symlink_dir (:obj:`str`, optional):
            Create a symlink to the file in this directory.
        clobber (:obj:`bool`, optional):
            If the output file already exists, this will force the
            assessments to be redone and the output file to be
            overwritten.
        checksum (:obj:`bool`, optional):
            Compare the checksum when reading the data.
        loggers (:obj:`list`, optional):
            List of `logging.Logger`_ objects to log progress;
            ignored if quiet=True. Logging is done using
            :func:`mangadap.util.log.log_output`. Default is no
            logging.
        quiet (:obj:`bool`, optional):
            Suppress terminal and logging output.

    Attributes:
        method (:class:`ReductionAssessmentDef`):
            Parameters defining the method to use for the reduction
            assessments.
        cube (:class:`mangadap.datacube.datacube.DataCube`):
            Datacube to analyze.
        pa (:obj:`float`):
            On-sky position angle of the major axis used to calculate
            elliptical, semi-major-axis coordinates, defined as the
            angle from North through East and denoted :math:`\phi_0`.
        ell (:obj:`float`):
            Ellipticity defined as :math:`\varepsilon=1-b/a`, based
            on the semi-minor to semi-major axis ratio (:math:`b/a`)
            of the isophotal ellipse used to calculate elliptical,
            semi-major-axis coordinates.
        directory_path (:obj:`str`):
            The exact path for the output file.
        output_file (:obj:`str`):
            The name of the file for the computed assessments. The
            full path of the output file will be
            :attr:`directory_path`/:attr:`output_file`.
        hardcopy (:obj:`bool`):
            Flag to write the data to a fits file.
        symlink_dir (:obj:`str`):
            Create a symlink to the file in this directory.
        hdu (`astropy.io.fits.HDUList`_):
            HDUList with the data, with columns as described above.
        correlation (:class:`mangadap.util.covariance.Covariance`):
            If requested, provides the spatial correlation matrix for
            the flux at a fiducial wavelength in the datacube. See
            :attr:`covar_wave` and :attr:`covar_channel`.
        covar_wave (:obj:`float`):
            The wavelength in angstroms at which the correlation
            matrix (:attr:`correlation`) was calculated, if
            available. Can be None if the covariance matrix is pulled
            directly from :attr:`cube`.
        covar_channel (:obj:`int`):
            The wavelength channel at which the correlation matrix
            (:attr:`correlation`) was calculated, if available.
            Can be None if the covariance matrix is pulled
            directly from :attr:`cube`.
        loggers (:obj:`list`):
            List of `logging.Logger`_ objects to log progress;
            ignored if quiet=True. Logging is done using
            :func:`mangadap.util.log.log_output`.
        quiet (:obj:`bool`):
            Suppress terminal and logging output.

    """
    def __init__(self, method_key, cube, pa=None, ell=None, method_list=None, dapver=None,
                 analysis_path=None, directory_path=None, output_file=None, hardcopy=True,
                 symlink_dir=None, clobber=False, checksum=False, loggers=None, quiet=False):
                 
        self.loggers = None
        self.quiet = False

        # Define the method properties
        self.method = None
        self._define_method(method_key, method_list=method_list)

        # Set in compute via _set_paths
        self.cube = None

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
        self.compute(cube, pa=pa, ell=ell, dapver=dapver, analysis_path=analysis_path,
                     directory_path=directory_path, output_file=output_file, hardcopy=hardcopy,
                     symlink_dir=symlink_dir, clobber=clobber, loggers=loggers, quiet=quiet)

    def __getitem__(self, key):
        return self.hdu[key]

    def _define_method(self, method_key, method_list=None):
        """
        Select the assessment method from the provided list.

        Used to set :attr:`method`; see
        :func:`mangadap.proc.util.select_proc_method`.

        Args:
            method_key (:obj:`str`):
                Keyword of the selected method. Available methods are
                provided by :func:`available_reduction_assessments`
            method_list (:obj:`list`, optional):
                List of :class:`ReductionAssessmentDef` objects that
                define the parameters required to assess the reduced
                data.
        """
        self.method = select_proc_method(method_key, ReductionAssessmentDef,
                                         method_list=method_list,
                                         available_func=available_reduction_assessments)

    def _set_paths(self, directory_path, dapver, analysis_path, output_file):
        """
        Set the I/O paths.

        Used to set :attr:`directory_path` and :attr:`output_file`.
        If not provided, the defaults are set using, respectively,
        :func:`mangadap.config.defaults.dap_common_path` and
        :func:`mangadap.config.defaults.dap_file_name`.

        Args:
            directory_path (:obj:`str`):
                The exact path to the DAP reduction assessments file.
                See :attr:`directory_path`.
            dapver (:obj:`str`):
                DAP version.
            analysis_path (:obj:`str`):
                The path to the top-level directory containing the
                DAP output files for a given DRP and DAP version.
            output_file (:obj:`str`):
                The name of the file with the reduction assessments.
                See :func:`compute`.
        """
        self.directory_path, self.output_file \
                = ReductionAssessment.default_paths(self.cube.plate, self.cube.ifudesign,
                                                    self.method['key'],
                                                    directory_path=directory_path,
                                                    drpver=self.cube.drpver, dapver=dapver,
                                                    analysis_path=analysis_path,
                                                    output_file=output_file)

    def _initialize_primary_header(self, hdr=None):
        """
        Constuct the primary header for the reference file.

        Args:
            hdr (`astropy.io.fits.Header`_, optional):
                Input base header for added keywords. If None, uses
                the :attr:`cube` header (if there is one) and then
                cleans the header using
                :func:`mangadap.util.fitsutil.DAPFitsUtil.clean_dap_primary_header`.

        Returns:
             `astropy.io.fits.Header`_: Initialized header object.
        """
        # Copy the from the DRP and clean it
        if hdr is None:
            hdr = self.cube.prihdr.copy()
            hdr = DAPFitsUtil.clean_dap_primary_header(hdr)
        hdr['AUTHOR'] = 'Kyle B. Westfall <westfall@ucolick.org>'
        hdr['RDXQAKEY'] = (self.method['key'], 'Method keyword')
        if self.pa is not None:
            hdr['ECOOPA'] = (self.pa, 'Position angle for ellip. coo')
        if self.ell is not None:
            hdr['ECOOELL'] = (self.ell, 'Ellipticity (1-b/a) for ellip. coo')
        if self.method['covariance']:
            hdr['BBWAVE'] = ('None' if self.covar_wave is None else self.covar_wave,
                             'Covariance channel wavelength')
            hdr['BBINDEX'] = ('None' if self.covar_channel is None else self.covar_channel,
                              'Covariance channel index')
        return hdr

    @staticmethod
    def default_paths(plate, ifudesign, method_key, directory_path=None, drpver=None, dapver=None,
                      analysis_path=None, output_file=None):
        """
        Set the default directory and file name for the output file.

        Args:
            plate (:obj:`int`):
                Plate number
            ifudesign (:obj:`int`):
                IFU design number
            method_key (:obj:`str`):
                Keyword designating the method used for the reduction
                assessments.
            directory_path (:obj:`str`, optional):
                The exact path to the DAP reduction assessments file.
                Default set by
                :func:`mangadap.config.defaults.dap_common_path`.
            drpver (:obj:`str`, optional):
                DRP version. Default set by
                :func:`mangadap.config.defaults.drp_version`.
            dapver (:obj:`str`, optional):
                DAP version. Default set by
                :func:`mangadap.config.defaults.dap_version`.
            analysis_path (:obj:`str`, optional):
                The path to the top-level directory containing the
                DAP output files for a given DRP and DAP version.
                Default set by
                :func:`mangadap.config.defaults.dap_analysis_path`.
            output_file (:obj:`str`, optional):
                The name of the file with the reduction assessments.
                Default set by
                :func:`mangadap.config.defaults.dap_file_name`.

        Returns:
            :obj:`tuple`: Returns two :obj:`str` objects with the
            path for the output file and the name of the output file.
        """
        _directory_path = defaults.dap_common_path(plate=plate, ifudesign=ifudesign,
                                                   drpver=drpver, dapver=dapver,
                                                   analysis_path=analysis_path) \
                                        if directory_path is None else str(directory_path)
        _output_file = defaults.dap_file_name(plate, ifudesign, method_key) \
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

    def compute(self, cube, pa=None, ell=None, dapver=None, analysis_path=None, directory_path=None,
                output_file=None, hardcopy=True, symlink_dir=None, clobber=False, loggers=None,
                quiet=False):
        r"""
        Compute and output the main data products.  The list of HDUs
        are:

            - ``PRIMARY`` : Empty apart from the header information.
            - ``SPECTRUM`` : Extension with the main, per-spectrum
              measurements; see :class:`ReductionAssessmentDataTable`.
            - ``CORREL`` : The correlation matrix between the ``SIGNAL``
              measurements provided in the ``SPECTRUM`` extension.  The
              format of this extension is identical to the nominal
              output of the :class:`mangadap.util.covariance.Covariance`
              object; see
              :func:`mangadap.util.covariance.Covariance.write`.

        Args:
            cube (:class:`mangadap.datacube.datacube.DataCube`):
                Datacube to analyze.
            pa (:obj:`float`, optional):
                On-sky position angle of the major axis used to
                calculate elliptical, semi-major-axis coordinates,
                defined as the angle from North through East and
                denoted :math:`\phi_0`.
            ell (:obj:`float`, optional):
                Ellipticity defined as :math:`\varepsilon=1-b/a`,
                based on the semi-minor to semi-major axis ratio
                (:math:`b/a`) of the isophotal ellipse used to
                calculate elliptical, semi-major-axis coordinates.
            method_list (:obj:`list`, optional):
                List of :class:`ReductionAssessmentDef` objects that
                define the parameters required to perform the
                assessments of the DRP file. The ``method_key`` must
                select one of these objects. Default is to construct
                the list using
                :func:`available_reduction_assessments`.
            dapver (:obj:`str`, optional):
                DAP version, which is used to define the default DAP
                analysis path. (This does **not** define the version
                of the code to use.) Default is defined by
                :func:`mangadap.config.defaults.dap_version`
            analysis_path (:obj:`str`, optional):
                The path to the top-level directory for the DAP
                output files for a given DRP and DAP version. Default
                is defined by
                :func:`mangadap.config.defaults.dap_analysis_path`.
            directory_path (:obj:`str`, optional):
                The exact path for the output file. Default is
                defined by
                :func:`mangadap.config.defaults.dap_common_path`.
            output_file (:obj:`str`, optional):
                The name of the file for the computed assessments.
                The full path of the output file will be
                :attr:`directory_path`/:attr:`output_file`. Default
                is defined by
                :func:`mangadap.config.defaults.dap_file_name`.
            hardcopy (:obj:`bool`, optional):
                Flag to write the data to a fits file.
            symlink_dir (:obj:`str`, optional):
                Create a symlink to the file in this directory.
            clobber (:obj:`bool`, optional):
                If the output file already exists, this will force the
                assessments to be redone and the output file to be
                overwritten.
            loggers (:obj:`list`, optional):
                List of `logging.Logger`_ objects to log progress;
                ignored if quiet=True. Logging is done using
                :func:`mangadap.util.log.log_output`. Default is no
                logging.
            quiet (:obj:`bool`, optional):
                Suppress terminal and logging output.

        Raises:
            TypeError:
                Raised if the provided ``cube`` is not a
                :class:`mangadap.datacube.datacube.DataCube`.
            ValueError:
                Raised if the output path cannot be defined.
        """

        # Initialize the reporting
        if loggers is not None:
            self.loggers = loggers
        self.quiet = quiet

        # Check the input cube
        if not isinstance(cube, DataCube):
            raise TypeError('Must provide a valid DataCube object!')
        self.cube = cube

        # Test if the RSS file exists; cannot compute covariance if not
        if self.method['covariance'] and not self.cube.can_compute_covariance:
            warnings.warn('Datacube does not provide a covariance matrix and could not load RSS '
                          'counterpart for the covariance computation.  Continuing by assuming '
                          '**no** covariance calculations.')
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

            # Read the header data
            try:
                self.pa = self.hdu['PRIMARY'].header['ECOOPA']
            except:
                if not self.quiet:
                    warnings.warn('Unable to read position angle from file header!')
                self.pa = None
            if not self.quiet and pa is not None and self.pa != pa:
                warnings.warn('Provided position angle different from available file; set ' \
                              'clobber=True to overwrite.')

            try:
                self.ell = self.hdu['PRIMARY'].header['ECOOELL']
            except:
                if not self.quiet:
                    warnings.warn('Unable to read ellipticity from file header!')
                self.ell = None
            if not self.quiet and ell is not None and self.ell != ell:
                warnings.warn('Provided ellipticity different from available file; set ' \
                              'clobber=True to overwrite.')

            if self.method['covariance']:
                # Construct the correlation matrix with the appropriate
                # variance
                self.correlation = Covariance.from_fits(self.hdu, ivar_ext=None,
                                                        correlation=True).apply_new_variance(
                                                            self.hdu['SPECTRUM'].data['VARIANCE'])
                self.covar_wave = self.hdu['PRIMARY'].header['BBWAVE']
                if isinstance(self.covar_wave, str):
                    self.covar_wave = eval(self.covar_wave)
                self.covar_channel = self.hdu['PRIMARY'].header['BBINDEX']
                if isinstance(self.covar_channel, str):
                    self.covar_channel = eval(self.covar_channel)
            else:
                self.correlation = None
                self.covar_wave = None
                self.covar_channel = None

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

        # Initialize the record array for the SPECTRUM extension
#        spectrum_data = init_record_array(self.cube.nspec, self._per_spectrum_dtype())
        spectrum_data = ReductionAssessmentDataTable(shape=self.cube.nspec)
        spectrum_data['SPAT_INDX'] = numpy.asarray(tuple(self.cube.spatial_index))
        spectrum_data['SKY_COO'][:,0], spectrum_data['SKY_COO'][:,1] \
                = (x.ravel() for x in self.cube.mean_sky_coordinates())

        coord_conversion = SemiMajorAxisCoo(xc=0, yc=0, rot=0, pa=self.pa, ell=self.ell)
        spectrum_data['ELL_COO'][:,0], spectrum_data['ELL_COO'][:,1] \
                = coord_conversion.polar(spectrum_data['SKY_COO'][:,0],
                                         spectrum_data['SKY_COO'][:,1])

        # Get the number of good wavelength channels in each spaxel
        # TODO: Move this into DataCube?
        flags = self.cube.do_not_use_flags()
        flux = self.cube.copy_to_masked_array(flag=flags)
        spectrum_data['FGOODPIX'] = numpy.sum(numpy.logical_not(numpy.ma.getmaskarray(flux)),
                                              axis=1) / flux.shape[1]

        # Flag spaxels with a flux that doesn't vary
        # TODO: Is this superfluous now?
#        frange = numpy.ma.amax(flux, axis=1)-numpy.ma.amin(flux, axis=1)
#        spectrum_data['MINEQMAX'] = (numpy.logical_not(numpy.ma.getmaskarray(frange))) \
#                                        & (numpy.ma.absolute(frange) < 1e-10)

        # Set the wavelength range in the statistics. If a response
        # function is used, only use the wavelength region where the
        # response function is provided.
        waverange = [self.cube.wave[0], self.cube.wave[-1]] \
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
        # covariance measurements. This is done by first assuming the
        # channel at the center of the wavelength range and then
        # flipping back and forth across this channel until a valid
        # channel is found or you reach the limits of the defined
        # wavelength range. In the latter case an exception is raised!
        if self.method['covariance'] and self.cube.covar is not None:
            self.correlation = self.cube.covar.copy()
            self.correlation.to_correlation()
            self.covar_wave = None
            self.covar_channel = None
        elif self.method['covariance'] and self.cube.covar is None:
            self.covar_wave = self.cube.central_wavelength(
                                        waverange=self.method['waverange'],
                                        response_func=self.method['response_func'], flag=flags)
            self.covar_channel = numpy.argsort(numpy.absolute(self.cube.wave-self.covar_wave))[0]
            goodfrac = spectrum_data['FGOODPIX'].reshape(self.cube.spatial_shape) \
                            > self.method['minimum_frac']
            off = 1
            sign = 1
            while self.covar_wave > waverange[0] and self.covar_wave < waverange[1] \
                    and numpy.sum(goodfrac & numpy.logical_not(
                                    self.cube.ivar[:,:,self.covar_channel] > 0)) > 0:
                self.covar_channel += sign*off
                self.covar_wave = self.cube.wave[self.covar_channel]
                sign *= -1
                off += 1

            if numpy.sum(goodfrac 
                         & numpy.logical_not(self.cube.ivar[:,:,self.covar_channel] > 0)) > 0:
                raise ValueError('Unable to find wavelength channel within fully valid data.')
            if not self.quiet:
                log_output(self.loggers, 1, logging.INFO,
                           'Covariance measured at wavelength: {0} (channel {1})'.format(
                                self.covar_wave, self.covar_channel))
            # Calculate the covariance matrix
            self.correlation = self.cube.covariance_matrix(self.covar_channel)
            self.correlation.to_correlation()
        else:
            self.correlation = None
            self.covar_wave = None
            self.covar_channel = None

        # Determine the statistics of the flux in the specified
        # wavelength range
        spectrum_data['SIGNAL'], spectrum_data['VARIANCE'], spectrum_data['SNR'] \
                = (x.ravel() for x in 
                        self.cube.flux_stats(waverange=self.method['waverange'],
                                             response_func=self.method['response_func'],
                                             flag=flags))

        # Force the covariance matrix to use the mean variance instead
        # of the variance in the single, specified channel.
        if self.method['covariance']:
            self.correlation = self.correlation.apply_new_variance(spectrum_data['VARIANCE'])

        # Construct header
        hdr = self._initialize_primary_header()

        # Get the covariance hdu
        if self.method['covariance']:
            hdr, ivar_hdu, covar_hdu = self.correlation.output_hdus(reshape=True, hdr=hdr)

        # Get the main extension columns and construct the HDUList
        self.hdu = fits.HDUList([fits.PrimaryHDU(header=hdr),
                                 spectrum_data.to_hdu(name='SPECTRUM')])

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

