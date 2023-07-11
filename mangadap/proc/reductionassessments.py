"""
Class that performs basic assessments of the reduced data products needed for
handling of the data by the DAP.  These assessments need only be done once per
data file.

----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import inspect
from pathlib import Path
import warnings
import logging

from IPython import embed

import numpy

from astropy.io import fits

from pydl.goddard.astro import airtovac

#from ..datacube import DataCube
from ..par.parset import KeywordParSet

from ..config import defaults

from ..util.datatable import DataTable
from ..util.fitsutil import DAPFitsUtil
from ..util.covariance import Covariance
from ..util.geometry import SemiMajorAxisCoo
from ..util.fileio import create_symlink
from ..util.log import log_output

# TODO: Add center coordinates
#            center_coo (:obj:`tuple`, optional):
#                A two-tuple with the coordinates in right-ascension
#                and declination for the coordinate-frame origin. If
#                None, no offset is performed.

# TODO: Fix the fact that if the correlation matrix is taken directly
# from the datacube that the wavelength and channel won't be defined.
# TODO: Change 'waverange' to 'pixelmask'?
class ReductionAssessmentDef(KeywordParSet):
    """
    Define how to assess the reduced data.

    The defined parameters are:

    .. include:: ../tables/reductionassessmentdef.rst

    """
    def __init__(self, key='SNRG', waverange=None, response_func_file='gunn_2001_g_response.db',
                 in_vacuum=True, covariance=True, minimum_frac=0.8, overwrite=False):

        # Use the signature to get the parameters and the default values
        sig = inspect.signature(self.__class__)
        pars = list(sig.parameters.keys())
        defaults = [sig.parameters[key].default for key in pars]

        # Remaining definitions done by hand
        ar_like = [ numpy.ndarray, list ]
        in_fl = [ int, float ]

        values = [key, waverange, response_func_file, in_vacuum, covariance, minimum_frac,
                  overwrite]
        dtypes = [str, ar_like, str, bool, bool, in_fl, bool]
        descr = ['Keyword to distinguish the assessment method.',
                 'A two-element vector with the starting and ending wavelength within which to ' \
                    'calculate the signal-to-noise.  Mutually exclusive with response_func_file',
                 'The name of a file that defines a response function to use for the S/N ' \
                    'calculation.  Must be a local file or distributed with the DAP source code.' \
                    '  Expected to have two columns, with the wavelength and response ' \
                    'efficiency.  Mutually exclusive with waverange.',
                 'Boolean indicating that the wavelengths provided either using waverange or ' \
                    'response_func_file are in vacuum (not air).',
                 'Provide the fiducial spatial covariance.  If this is False, no spatial ' \
                     'covariance will be available for calculations that use the results of the ' \
                     'reduction assessments.  If True and a covariance matrix is available ' \
                     'directly from the datacube, it will be used.  If True and the datacube ' \
                     'does not already provide a computed covariance matrix, one is calculated ' \
                     'using :func:`~mangadap.datacube.datacube.DataCube.covariance_matrix` ' \
                     'method.  **WARNING**: If either of these fail, the DAP will issue a ' \
                     'warning and continue assuming no spatial covariance.',
                 'Minimum fraction of unmasked pixels in a spectrum required for inclusion in ' \
                     'the spatial covariance calculation.  Note this should match the value ' \
                     'used for the spatial-binning module.',
                 'If the output file already exists, redo all the calculations and overwrite it.']

        super().__init__(pars, values=values, defaults=defaults, dtypes=dtypes, descr=descr)
        self.response = None
        self._validate()

    def _validate(self):
        """
        Validate the parameters.
        """
        if self['response_func_file'] in ['none', 'None']:
            self['response_func_file'] = None
        if self['waverange'] is not None and self['response_func_file'] is not None:
            warnings.warn('You have defined both a wavelength range and a response function for '
                          'the reduction assessments; the latter takes precedence, the former is '
                          'ignored!')
            self['waverange'] = None
        if self['waverange'] is None and self['response_func_file'] is None:
            raise ValueError('Reduction assessment method undefined.  Must provide either '
                             '\'waverange\' or \'response_func_file\'.')

        if self['waverange'] is not None:
            if len(self['waverange']) != 2:
                raise ValueError('Defined wavelength range must be a two-element list/array.')
            # Wavelength range was defined
            if not self['in_vacuum']:
                # Convert it to vacuum
                self['waverange'] = airtovac(self['waverange'])
            return

        response_file = Path(self['response_func_file']).resolve()
        if not response_file.exists():
            # File not found, so try to find it in the configuration directory
            response_file = defaults.dap_data_root() / 'filter_response' \
                                / self['response_func_file']
        if not response_file.exists():
            raise FileNotFoundError(f'Could not find {self["response_func_file"]} '
                                    'in the local directory or the DAP source distribution.')

        self.response = numpy.genfromtxt(str(response_file))[:,:2]
        if not self['in_vacuum']:
            self.response[:,0] = airtovac(self.response[:,0])


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
        method (:class:`~mangadap.proc.reductionassessments.ReductionAssessmentDef`):
            Object that defines the main parameters used for the assessment
            method.
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
        reff (:obj:`float`, optional):
            Effective radius of the galaxy.  If None, set to 1.
        output_path (:obj:`str`, `Path`_, optional):
            The path for the output file.  If None, the current working
            directory is used.
        output_file (:obj:`str`, optional):
            The name of the output "reference" file. The full path of the output
            file will be :attr:`directory_path`/:attr:`output_file`.  If None,
            the default is to combine ``cube.output_root`` and the
            ``method_key``.  See
            :func:`~mangadap.proc.reductionassessments.ReductionAssessments.default_paths`.
        hardcopy (:obj:`bool`, optional):
            Flag to write the data to a fits file.
        symlink_dir (:obj:`str`, optional):
            Create a symlink to the file in this directory.
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
        reff (:obj:`float`):
            Effective radius of the galaxy.
        output_path (:obj:`str`):
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
    def __init__(self, method, cube, pa=None, ell=None, reff=None, output_path=None,
                 output_file=None, hardcopy=True, symlink_dir=None, checksum=False, loggers=None,
                 quiet=False):
                 
        self.loggers = None
        self.quiet = False

        # Define the method properties
#        self.method = None
#        self._define_method(method_key, method_list=method_list)
        self.method = method
        if not isinstance(self.method, ReductionAssessmentDef):
            raise TypeError('Method must have type ReductionAssessmentDef.')

        # Set in compute via _set_paths
        self.cube = None

        # Define the output directory and file.  Set by default_paths() in compute()
        self.directory_path = None
        self.output_file = None
        self.hardcopy = None
        self.symlink_dir = None

        # Initialize the objects used in the assessments
        self.pa = None
        self.ell = None
        self.reff = None
        self.hdu = None
        self.checksum = checksum
        self.correlation = None
        self.covar_wave = None
        self.covar_channel = None

        # Run the assessments of the DRP file
        self.compute(cube, pa=pa, ell=ell, reff=reff, output_path=output_path,
                     output_file=output_file, hardcopy=hardcopy, symlink_dir=symlink_dir,
                     loggers=loggers, quiet=quiet)

    def __getitem__(self, key):
        return self.hdu[key]

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
        if self.reff is not None:
            hdr['ECOOREFF'] = (self.reff, 'Effective radius')
        if self.method['covariance']:
            hdr['BBWAVE'] = ('None' if self.covar_wave is None else self.covar_wave,
                             'Covariance channel wavelength')
            hdr['BBINDEX'] = ('None' if self.covar_channel is None else self.covar_channel,
                              'Covariance channel index')
        return hdr

    @staticmethod
    def default_paths(cube, method_key, output_path=None, output_file=None):
        """
        Set the default directory and file name for the output file.

        Args:
            cube (:class:`mangadap.datacube.datacube.DataCube`):
                Datacube to analyze.
            method_key (:obj:`str`):
                Keyword designating the method used for the reduction
                assessments.
            output_path (:obj:`str`, `Path`_, optional):
                The path for the output file.  If None, the current working
                directory is used.
            output_file (:obj:`str`, optional):
                The name of the output "reference" file. The full path of the output
                file will be :attr:`directory_path`/:attr:`output_file`.  If None,
                the default is to combine ``cube.output_root`` and the
                ``method_key``.

        Returns:
            :obj:`tuple`: Returns a `Path`_ with the output directory and a
            :obj:`str` with the output file name.
        """
        directory_path = Path('.').resolve() if output_path is None \
                                else Path(output_path).resolve()
        _output_file = f'{cube.output_root}-{method_key}.fits.gz' if output_file is None \
                                else output_file
        return directory_path, _output_file

    def file_name(self):
        """Return the name of the output file."""
        return self.output_file

    def file_path(self):
        """Return a `Path`_ object with the full path to the output file."""
        if self.directory_path is None or self.output_file is None:
            return None
        return self.directory_path / self.output_file

    def info(self):
        return self.hdu.info()

    def compute(self, cube, pa=None, ell=None, reff=None, output_path=None, output_file=None, hardcopy=True,
                symlink_dir=None, overwrite=None, loggers=None, quiet=False):
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
            reff (:obj:`float`, optional):
                Effective radius of the galaxy.  If None, set to 1.
            output_path (:obj:`str`, `Path`_, optional):
                The path for the output file.  If None, the current working
                directory is used.
            output_file (:obj:`str`, optional):
                The name of the output "reference" file. The full path of the
                output file will be :attr:`directory_path`/:attr:`output_file`.
                If None, the default is to combine ``cube.output_root`` and the
                ``method_key``.  See
                :func:`~mangadap.proc.reductionassessments.ReductionAssessments.default_paths`.
            hardcopy (:obj:`bool`, optional):
                Flag to write the data to a fits file.
            symlink_dir (:obj:`str`, optional):
                Create a symlink to the file in this directory.
            overwrite (:obj:`bool`, optional):
                If the output file already exists, this will force the
                assessments to be redone and the output file to be
                overwritten.  If None, uses the value in :attr:`method`.
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
#        if not isinstance(cube, DataCube):
#            raise TypeError('Must provide a valid DataCube object!')
        self.cube = cube

        # Test if the RSS file exists; cannot compute covariance if not
        if self.method['covariance'] and not self.cube.can_compute_covariance:
            warnings.warn('Datacube does not provide a covariance matrix and could not load RSS '
                          'counterpart for the covariance computation.  Continuing by assuming '
                          '**no** covariance calculations.')
            self.method['covariance'] = False

        # Reset the output paths if necessary
        self.directory_path, self.output_file \
                = ReductionAssessment.default_paths(cube, self.method['key'],
                                                    output_path=output_path,
                                                    output_file=output_file)

        # Check that the path for or to the file is defined
        ofile = self.file_path()
        if ofile is None:
            raise ValueError('File path for output file is undefined!')

        # Report
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, '-'*50)
            log_output(self.loggers, 1, logging.INFO, f'{"REDUCTION ASSESSMENT COMPUTATIONS":^50}')
            log_output(self.loggers, 1, logging.INFO, '-'*50)
            log_output(self.loggers, 1, logging.INFO, f'Output path: {self.directory_path}')
            log_output(self.loggers, 1, logging.INFO, f'Output file: {self.output_file}')

        _overwrite = self.method['overwrite'] if overwrite is None else overwrite

        # If the file already exists, and not overwriting, just read the file
        self.symlink_dir = symlink_dir
        if ofile.exists() and not _overwrite:
            self.hardcopy = True
            if not self.quiet:
                log_output(self.loggers, 1, logging.INFO, 'Reading exiting file')

            # Read the existing output file
            self.hdu = DAPFitsUtil.read(str(ofile), checksum=self.checksum)

            # Read the header data
            try:
                self.pa = self.hdu['PRIMARY'].header['ECOOPA']
            except:
                if not self.quiet:
                    warnings.warn('Unable to read position angle from file header!')
                self.pa = None
            if not self.quiet and pa is not None and self.pa != pa:
                warnings.warn('Provided position angle different from available file; set ' \
                              'overwrite=True to overwrite.')

            try:
                self.ell = self.hdu['PRIMARY'].header['ECOOELL']
            except:
                if not self.quiet:
                    warnings.warn('Unable to read ellipticity from file header!')
                self.ell = None
            if not self.quiet and ell is not None and self.ell != ell:
                warnings.warn('Provided ellipticity different from available file; set ' \
                              'overwrite=True to overwrite.')

            try:
                self.reff = self.hdu['PRIMARY'].header['ECOOREFF']
            except:
                if not self.quiet:
                    warnings.warn('Unable to read effective radius from file header!')
                self.reff = None
            if not self.quiet and reff is not None and self.reff != reff:
                warnings.warn('Provided effective radius different from available file; set ' \
                              'overwrite=True to overwrite.')

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
                create_symlink(str(ofile), self.symlink_dir, overwrite=_overwrite)
            if not self.quiet:
                log_output(self.loggers, 1, logging.INFO, '-'*50)
            return

        # (Re)Initialize some of the attributes
        if pa is not None:
            self.pa = pa
        if ell is not None:
            self.ell = ell
        if reff is not None:
            self.reff = reff

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
                    if self.method['waverange'] is None and self.method.response is None \
                    else (self.method['waverange'] if self.method['waverange'] is not None else
                        [ self.method.response[0,0], self.method.response[-1,0] ])
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
            self.covar_wave = self.cube.central_wavelength(waverange=self.method['waverange'],
                                                           response_func=self.method.response,
                                                           flag=flags)
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
                                             response_func=self.method.response, flag=flags))

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
        if not self.directory_path.exists():
            self.directory_path.mkdir(parents=True)
        self.hardcopy = hardcopy
        if self.hardcopy:
            DAPFitsUtil.write(self.hdu, str(ofile), overwrite=_overwrite, checksum=True,
                              symlink_dir=self.symlink_dir, loggers=self.loggers)
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, '-'*50)



