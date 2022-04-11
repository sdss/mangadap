r"""
Class defining the I/O configuration for MaNGA files in the DAP.

----        

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import time
import warnings
from pathlib import Path
from configparser import ConfigParser

from IPython import embed

import numpy
from astropy.io import fits

from . import defaults
from . import manga_environ
from .analysisplan import AnalysisPlan
from ..util.parser import DefaultConfig
from ..util.constants import DAPConstants
from ..util.filter import interpolate_masked_vector

def drp_version():
    """
    Return the DRP version defined by the environmental variable
    MANGADRP_VER.
    """
    return manga_environ['MANGADRP_VER']


def drp_redux_path(drpver=None):
    """
    Return the main output path for the DRP products using the
    environmental variable MANGA_SPECTRO_REDUX.

    Args:
        drpver (:obj:`str`, optional):
            DRP version. Default is to use :func:`drp_version`.

    Returns:
        :obj:`str`: Path to reduction directory
    """
    # Make sure the DRP version is set
    if drpver is None:
        drpver = drp_version()
    return Path(manga_environ['MANGA_SPECTRO_REDUX']).resolve() / drpver


def drp_directory_path(plate, drpver=None, redux_path=None):
    """
    Return the exact directory path with the DRP file.

    Args:
        plate (:obj:`int`):
            Plate number
        drpver (:obj:`str`, optional):
            DRP version. Default is to use :func:`drp_version`.
        redux_path (:obj:`str`, optional):
            Path to the root reduction directory. Default is to use
            :func:`drp_redux_path`.

    Returns:
        :obj:`str`: Path to the directory with the 3D products of the
        DRP
    """
    # Make sure the redux path is set
    _redux_path = drp_redux_path(drpver=drpver) \
                        if redux_path is None else Path(redux_path).resolve()
    return _redux_path / str(plate) / 'stack'


def drp_finding_chart_path(plate, drpver=None, redux_path=None):
    """
    Return the exact directory path with the finding charts for a given plate.

    Args:
        plate (:obj:`int`):
            Plate number
        drpver (:obj:`str`, optional):
            DRP version. Default is to use :func:`drp_version`.
        redux_path (:obj:`str`, optional):
            Path to the root reduction directory. Default is to use
            :func:`drp_redux_path`.

    Returns:
        :obj:`str`: Path to the directory with the finding charts.
    """
    # Make sure the redux path is set
    _redux_path = drp_redux_path(drpver=drpver) \
                        if redux_path is None else Path(redux_path).resolve()
    return _redux_path / str(plate) / 'images'


def drpall_file(drpver=None, redux_path=None):
    """
    Return the path to the DRPall file.

    Args:
        drpver (:obj:`str`, optional):
            DRP version.  Default is to use :func:`drp_version`.
        redux_path (:obj:`str`, optional):
            Path to the root reduction directory. Default is to use
            :func:`drp_redux_path`.

    Returns:
        :obj:`str`: Full path to the DRPall fits file.
    """
    _drpver = drp_version() if drpver is None else drpver
    _redux_path = drp_redux_path(drpver=_drpver) \
                        if redux_path is None else Path(redux_path).resolve()
    return _redux_path / f'drpall-{_drpver}.fits'


def dap_version():
    """
    Return the DAP version defined by the environmental variable
    MANGADAP_VER.  If that environmental variable does not exist,
    `mangadap.__version__` is returned.
    """
    return manga_environ['MANGADAP_VER']


def dap_analysis_path(drpver=None, dapver=None):
    """
    Return the main output path for the DAP using the environmental
    variable MANGA_SPECTRO_ANALYSIS.

    Args:
        drpver (:obj:`str`, optional):
            DRP version. Default is to use :func:`drp_version`.
        dapver (:obj:`str`, optional):
            DAP version. Default is to use :func:`dap_version`.

    Returns:
        :obj:`str`: Path to analysis directory
    """
    # Make sure the DRP version is set
    if drpver is None:
        drpver = drp_version()
    # Make sure the DAP version is set
    if dapver is None:
        dapver = dap_version()
    return Path(manga_environ['MANGA_SPECTRO_ANALYSIS']).resolve() / drpver / dapver


def manga_fits_root(plate, ifudesign, mode=None):
    """
    Generate the main root name for the output MaNGA fits files for a
    given plate/ifudesign/mode.

    Args:
        plate (:obj:`int`):
            Plate number
        ifudesign (:obj:`int`):
            IFU design number
        mode (:obj:`str`, optional):
            Mode of the output fits file. Options are: ``'LINCUBE'``,
            ``'LINRSS'``, ``'LOGCUBE'``, ``'LOGRSS'``, or ``'MAPS'``.
            Default is that no mode is included in the name.

    Returns:
        :obj:`str`: Root name for a MaNGA fits file:
        ``manga-[PLATE]-[IFUDESIGN]-[MODE]``

    Raises:
        ValueError:
            Raised if mode is not a valid option.
    """
    if mode not in [ None, 'LINCUBE', 'LINRSS', 'LOGCUBE', 'LOGRSS', 'MAPS' ]:
        raise ValueError('Do not understand mode={0}.'.format(mode))
    return f'manga-{plate}-{ifudesign}' if mode is None else f'manga-{plate}-{ifudesign}-{mode}'


def plate_target_files():
    """
    Return the default plateTarget files in mangacore and their
    associated catalog indices. The catalog indices are determined
    assuming the file names are of the form::

        'plateTargets-{0}.par'.format(catalog_id)

    Returns:
        `numpy.ndarray`_: Two arrays: the first contains the
        identified plateTargets files found using the default search
        string, the second provides the integer catalog index
        determined for each file.
    """
    core_dir = Path(manga_environ['MANGACORE_DIR']).resolve()
    if not core_dir.exists():
        raise ValueError(f'MANGACORE directory does not exist! {core_dir}')

    # Default search string
    search_str =  core_dir / 'platedesign' / 'platetargets'
    file_list = sorted(list(search_path.glob('plateTargets*.par')))
    nfiles = len(file_list)
    if nfiles == 0:
        raise ValueError('Unable to find any plateTargets files!')
    trgid = numpy.zeros(nfiles, dtype=int)                  # Array to hold indices
    for i in range(nfiles):
        suffix = file_list[i].split('-')[1]                 # Strip out the '{i}.par'
        trgid[i] = int(suffix[:suffix.find('.')])           # Strip out '{i}'
    
    return numpy.array(file_list), trgid


def redshift_fix_file():
    """
    Return the path to the default redshift fix file.

    Returns:
        :obj:`str`: Expected path to the redshift-fix parameter file.
    """
    return defaults.dap_data_root() / 'fix' / 'redshift_fix.par'


def photometry_fix_file():
    """
    Return the path to the default photometry fix file.

    Returns:
        :obj:`str`: Expected path to the photometry-fix parameter file.
    """
    return defaults.dap_data_root() / 'fix' / 'photometry_fix.par'


def parse_plate_ifu_from_file(name):
    """
    Parse the plate and ifu numbers from the file name of a MaNGA DRP file.
    """
    return tuple([int(n) for n in name.split('-')[1:3]])


# USED BY:
#   - scripts/rundap.py
def dap_config(plate, ifudesign, drpver=None, dapver=None, analysis_path=None,
               directory_path=None):
    """
    Return the full path to the DAP configuration file.

    The configuration file provides the input data necessary to
    instantiate a :class:`mangadap.datacube.manga.MaNGADataCube`.
    
    Args:
        plate (:obj:`int`):
            Plate number
        ifudesign (:obj:`int`):
            IFU design number
        drpver (:obj:`str`, optional):
            DRP version. Default is to use :func:`drp_version`.
        dapver (:obj:`str`, optional):
            DAP version. Default is to use :func:`dap_version`.
        analysis_path (:obj:`str`, optional): 
            Path to the root analysis directory. Default is to use
            :func:`dap_analysis_path`.
        directory_path (:obj:`str`, optional):
            Path to the directory with the DAP output files. Default
            is to use :func:`dap_common_path`

    Returns:
        :obj:`str`: Full path to the DAP par file
    """
    # Make sure the directory path is defined
    _directory_path = dap_common_path(plate=plate, ifudesign=ifudesign, drpver=drpver,
                                      dapver=dapver, analysis_path=analysis_path) \
                            if directory_path is None else Path(directory_path).resolve()
    # Set the name of the par file; put this in its own function?
    return _directory_path /  f'{dap_file_root(plate, ifudesign, mode="CUBE")}.ini'


# USED BY:
#   - scripts/dapall_qa.py
#   - scripts/find_repeat_observations.py
#   - survey/dapall.py
def dapall_file(drpver=None, dapver=None, analysis_path=None):
    """
    Return the path to the DAPall file.

    Args:
        drpver (:obj:`str`, optional):
            DRP version.  Default is to use :func:`drp_version`.
        dapver (:obj:`str`, optional):
            DAP version.  Default is to use :func:`dap_version`.
        analysis_path (:obj:`str`, optional):
            Path to the root analysis directory.  Default is to use
            :func:`dap_analysis_path`

    Returns:
        :obj:`str`: Full path to the DAPall fits file.
    """
    _drpver = drp_version() if drpver is None else drpver
    _dapver = dap_version() if dapver is None else dapver
    _analysis_path = dap_analysis_path(drpver=_drpver, dapver=_dapver) \
                        if analysis_path is None else Path(analysis_path).resolve()
    return _analysis_path / f'dapall-{_drpver}-{_dapver}.fits'


class MaNGAConfig:

    instrument = 'manga'

    def __init__(self, plate, ifudesign, mode='CUBE', log=True, drpver=None, redux_path=None,
                 chart_path=None, directory_path=None):
        """
        Initialize the I/O configuration using the same file as used to read
        MaNGA data.
        """
        # MaNGA-specific attributes
        self.plate = plate
        self.ifudesign = ifudesign
        MaNGAConfig.check_mode(mode)
        self.mode = mode
        self.log = log
        self.drpver = drp_version() if drpver is None else drpver
        # TODO: Depending on what is provided to the function, redux path, chart
        # path, and directory path can all be inconsistent.
        self.redux_path = drp_redux_path(drpver=self.drpver) \
                            if redux_path is None else Path(redux_path).resolve()
        self.chart_path = drp_finding_chart_path(self.plate, drpver=self.drpver,
                                                          redux_path=self.redux_path) \
                                if chart_path is None else Path(chart_path).resolve()
        self.directory_path = drp_directory_path(self.plate, drpver=self.drpver,
                                                          redux_path=self.redux_path) \
                                if directory_path is None else Path(directory_path).resolve()
        # TODO: This is only relevant the CUBE data
        if self.directory_path != drp_directory_path(self.plate, drpver=self.drpver,
                                                     redux_path=redux_path):
            warnings.warn('Paths do not match MaNGA defaults.  May not be able to find paired '
                          'CUBE & RSS files, if requested.')

    @classmethod
    def from_config(cls, cfgfile, mode='CUBE'):
        """
        """
        # Read the configuration file.  This checks if the file exists and will
        # raise an exception if it does not.
        cfg = DefaultConfig(cfgfile, interpolate=True)

        # Set the attributes, forcing a known type
        plate = cfg.getint('plate')
        ifu = cfg.getint('ifu')
        if plate is None or ifu is None:
            raise ValueError('Configuration file must define the plate and IFU.')
        log = cfg.getbool('log', default=True)
        drpver = cfg.get('drpver')
        redux_path = cfg.get('redux_path')
        directory_path = cfg.get('directory_path')

        return cls(plate, ifu, mode=mode, log=log, drpver=drpver, redux_path=redux_path,
                   directory_path=directory_path)

    @classmethod
    def from_file(cls, datafile):
        """
        Initialize the I/O configuration from the datafile for the DAP to
        process.
        """
        ifile = Path(datafile).resolve()
        if not ifile.exists():
            raise FileNotFoundError(f'File does not exist: {ifile}')

        # Parse the relevant information from the filename
        plate, ifudesign, log = ifile.name.split('-')[1:4]

        with fits.open(str(ifile)) as hdu:
            drpver = hdu[0].header['VERSDRP3']

        return cls(int(plate), int(ifudesign), mode='CUBE' if 'CUBE' in log else 'RSS',
                   log='LOG' in log, drpver=drpver, directory_path=ifile.parent)

    def copy(self):
        """
        Create a deep copy.
        """
        return MaNGAConfig(self.plate, self.ifudesign, log=self.log, drpver=self.drpver,
                           directory_path=self.directory_path)

    @property
    def file_name(self):
        root = manga_fits_root(self.plate, self.ifudesign,
                               mode=f'{"LOG" if self.log else "LIN"}{self.mode}')
        return f'{root}.fits.gz'

    @property
    def file_path(self):
        return self.directory_path / self.file_name

    @property
    def image_path(self):
        return self.chart_path / f'{self.ifudesign}.png'

    @staticmethod
    def check_mode(mode):
        """
        Check that the mode is valid.

        Valide modes are set by :func:`mode_options`.

        Args:
            mode (:obj:`str`):
                Mode value to check.

        Raises:
            ValueError:
                Raised if the mode is undefined.
        """
        options = ['CUBE', 'RSS']
        if mode not in options:
            raise ValueError(f'Unknown mode {mode}.  Must be in: {options}')

    @property
    def output_root(self):
        return f'{self.instrument}-{self.plate}-{self.ifudesign}'

    @property
    def cfg_root(self):
        """
        Generate the root name of the MaNGA DAP configuration and script files.
    
        Returns:
            :obj:`str`: Root name for the DAP file:
            ``mangadap-[PLATE]-[IFUDESIGN]-LOG[MODE]``
        """
        return f'mangadap-{self.plate}-{self.ifudesign}-LOG{self.mode}'

    @staticmethod
    def propagate_flags():
        """Flags that should be propagated from the observed data to the analysis data."""
        return ['FORESTAR']

    @staticmethod
    def do_not_use_flags():
        """Return the maskbit names that should not be used."""
        return ['DONOTUSE', 'FORESTAR']

    @staticmethod
    def do_not_fit_flags():
        """Return the maskbit names that should not be fit."""
        return ['DONOTUSE', 'FORESTAR']

    @staticmethod
    def do_not_stack_flags():
        """Return the maskbit names that should not be stacked."""
        return ['DONOTUSE', 'FORESTAR']

    @staticmethod
    def cube_pixelscale():
        """
        Return the default pixel scale of the DRP CUBE files in arcsec.
        """
        return 0.5

    @staticmethod
    def cube_width_buffer():
        """
        Return the default width buffer in pixels used when regridding the DRP
        RSS spectra into the CUBE format.
        """
        return 10

    @staticmethod
    def cube_recenter():
        """
        Return the default recentering flag used when regridding the DRP RSS
        spectra into the CUBE format.
        """
        return False

    @staticmethod
    def regrid_rlim():
        """
        Return the default limiting radius for the Gaussian kernel used when
        regridding the DRP RSS spectra into the CUBE format.
        """
        return 1.6

    @staticmethod
    def regrid_sigma():
        """
        Return the default standard deviation of the Gaussian kernel used when
        regridding the DRP RSS spectra into the CUBE format.
        """
        return 0.7

    @staticmethod
    def spectral_resolution_extension(hdu, ext=None):
        """
        Determine the spectral resolution channel to use.

        Precedence follows this order: ``LSFPRE``, ``PRESPECRES``,
        ``LSFPOST``, ``SPECRES``.

        Args:
            hdu (`astropy.io.fits.HDUList`):
                The opened MaNGA DRP file.
            ext (:obj:`str`, optional):
                Specify the extension with the spectral estimate to
                use. Should be in None, ``LSFPRE``, ``PRESPECRES``,
                ``LSFPOST``, or ``SPECRES``. The default is None,
                which means it will return the extension found first
                in the order above. None is returned if none of the
                extensions are present.

        Returns:
            :obj:`str`: The name of the preferred extension to use.
        """
        allowed = ['LSFPRE', 'LSFPOST', 'PRESPECRES', 'SPECRES']
        if ext is not None and ext not in allowed:
            warnings.warn(f'{ext} is not a viable spectral resolution extension.')
            return None
        available = [h.name for h in hdu if h.name in allowed]
        _ext = ext
        if ext is None:
            _ext = 'LSFPRE'
            if _ext not in available:
                _ext = 'PRESPECRES'
            if _ext not in available:
                _ext = 'LSFPOST'
            if _ext not in available:
                _ext = 'SPECRES'
        return _ext

    @staticmethod
    def spectral_resolution(hdu, ext=None, fill=False, median=False):
        """
        Return the spectral resolution for all spectra.

        See :func:`spectral_resolution_extension` for a description
        of the precedence used when ``ext`` is None.

        Args:
            hdu (`astropy.io.fits.HDUList`):
                The opened MaNGA DRP file.
            ext (:obj:`str`, optional):
                Specify the extension with the spectral estimate to
                use. See :func:`spectral_resolution_extension`.
            fill (:obj:`bool`, optional):
                Fill masked values by interpolation.  Default is to
                leave masked pixels in returned array.
            median (:obj:`bool`, optional):
                Return a single vector with the median spectral
                resolution instead of a per spectrum array. When using
                the ``SPECRES`` extension, this just returns the vector
                provided by the DRP file; when using either of the
                ``LSF`` extensions, this performs a masked median across
                the array and then interpolates any wavelengths that
                were masked in all vectors.

        Returns:
            :obj:`tuple`: Returns a :obj:`str` with the name of the
            extension used for the spectral-resolution measurements
            and a `numpy.ma.MaskedArray`_ with the spectral
            resolution data. Even if interpolated such that there
            should be no masked values, the function returns a masked
            array. Array contains the spectral resolution (:math:`R =
            \lambda/\Delta\lambda`) pulled from the DRP file.
        """
        # Determine which spectral resolution element to use
        _ext = MaNGAConfig.spectral_resolution_extension(hdu, ext=ext)
        # If no valid extension, raise an exception
        if ext is None and _ext is None:
            raise ValueError('No valid spectral resolution extension.')
        if ext is not None and _ext is None:
            raise ValueError('No extension: {0}'.format(ext))
        warnings.warn('Extension {0} used to define spectral resolution.'.format(_ext))
            
        # Set the mode based on the shape of the flux extension
        mode = 'CUBE' if hdu['FLUX'].data.ndim == 3 else 'RSS'

        # Build the spectral resolution vectors
        if 'SPECRES' in _ext:
            sres = numpy.ma.MaskedArray(hdu[_ext].data, mask=numpy.invert(hdu[_ext].data > 0))
            if fill:
                sres = numpy.ma.MaskedArray(interpolate_masked_vector(sres))
            if not median:
                sres = numpy.tile(sres, (hdu['FLUX'].data.shape[0],1)) if mode == 'RSS' \
                         else numpy.tile(sres, (*hdu['FLUX'].data.shape[1:][::-1],1)).T
            return _ext, sres

        # Otherwise dealing with the DISP data
        sres = numpy.ma.MaskedArray(hdu[_ext].data)
        # Mask any non-positive value
        sres[numpy.invert(sres > 0)] = numpy.ma.masked
        # Convert from sigma in angstroms to spectral resolution (based
        # on FWHM). Make sure to treat array shapes correctly for the
        # different modes.
        if mode == 'RSS':
            sres = numpy.ma.divide(hdu['WAVE'].data[None,:], sres) / DAPConstants.sig2fwhm
        else:
            # Mode is 'CUBE'
            sres = numpy.ma.divide(hdu['WAVE'].data[:,None,None], sres) / DAPConstants.sig2fwhm
        # Interpolate over any masked values
        if fill:
            axis = 1 if mode == 'RSS' else 0
            sres = numpy.ma.MaskedArray(
                        numpy.ma.apply_along_axis(interpolate_masked_vector, axis, sres))
        if median:
            sres = numpy.ma.median(sres, axis=0) if mode == 'RSS' \
                        else numpy.ma.median(sres.reshape(sres.shape[0],-1), axis=1)
        return _ext, sres

    @staticmethod
    def write_config(ofile, plate, ifudesign, log=True, z=None, vdisp=None, ell=None, pa=None,
                     reff=None, sres_ext=None, sres_fill=None, covar_ext=None, drpver=None,
                     redux_path=None, directory_path=None, overwrite=True):
        """
        Write the configuration file that can be used to instantiate
        a MaNGA data object
        (:class:`mangadap.datacube.manga.MaNGADataCube` or
        :class:`mangadap.spectra.manga.MaNGARSS`).
        
        See :func:`from_config`.

        Args:
            ofile (:obj:`str`, optional):
                Name of the configuration file.
            plate (:obj:`int`):
                Plate number
            ifudesign (:obj:`int`):
                IFU design
            log (:obj:`bool`, optional):
                Use the datacube that is logarithmically binned in
                wavelength.
            z (:obj:`float`, optional):
                Estimated bulk redshift. If None, some of the DAP
                analysis modules will fault.
            vdisp (:obj:`float`, optional):
                Estimated velocity dispersion. If None, some of the
                DAP analysis modules will assume an initial guess of
                100 km/s.
            ell (:obj:`float`, optional):
                Characteristic isophotal ellipticity (1-b/a). If
                None, some of the DAP modules will issue a warning
                and continue by assuming ``ell=0``.
            pa (:obj:`float`, optional):
                Characteristic isophotal position angle (through E
                from N). If None, some of the DAP modules will issue
                a warning and continue by assuming ``pa=0``.
            reff (:obj:`float`, optional):
                Effective (half-light) radius in arcsec. If None,
                some of the DAP modules will issue a warning and
                continue by assuming ``reff=1``.
            sres_ext (:obj:`str`, optional):
                The extension to use when constructing the spectral
                resolution vectors. See :func:`spectral_resolution`.
            sres_fill (:obj:`bool`, optional):
                Fill masked values by interpolation. Default is to
                leave masked pixels in returned array.
            covar_ext (:obj:`str`, optional):
                Extension to use as the single spatial correlation
                matrix for all wavelength channels, read from the DRP
                file. For generating the covariance matrix directly
                for an arbitrary wavelength channel using the RSS
                file, see
                :func:`mangadap.datacube.datacube.DataCube.covariance_matrix`.
            drpver (:obj:`str`, optional):
                DRP version, which is used to define the default DRP
                redux path. Default is defined by
                :func:`mangadap.config.manga.drp_version`
            redux_path (:obj:`str`, optional):
                The path to the top level directory containing the
                DRP output files for a given DRP version. Default is
                defined by
                :func:`mangadap.config.manga.drp_redux_path`.
            directory_path (:obj:`str`, optional):
                The exact path to the DRP file. Default is defined by
                :func:`mangadap.config.manga.drp_directory_path`.
                Providing this ignores anything provided for
                ``drpver`` or ``redux_path``.
            overwrite (:obj:`bool`, optional):
                Overwrite any existing parameter file.
        """
        _ofile = Path(ofile).resolve()
        if _ofile.exists() and not overwrite:
            raise FileExistsError(f'{_ofile} already exists; to overwrite, set '
                                  'overwrite=True.')

        if not _ofile.parent.exists():
            _ofile.parent.mkdir(parents=True)

        # Build the configuration data
        cfg = ConfigParser(allow_no_value=True)
        cfg['default'] = {'drpver': drpver,
                          'redux_path': redux_path,
                          'directory_path': directory_path,
                          'plate': str(plate),
                          'ifu': str(ifudesign),
                          'log': str(log),
                          'sres_ext': sres_ext,
                          'sres_fill': sres_fill,
                          'covar_ext': covar_ext,
                          'z': None if z is None else '{0:.7e}'.format(z),
                          'vdisp': None if vdisp is None else '{0:.7e}'.format(vdisp),
                          'ell': None if ell is None else '{0:.7e}'.format(ell),
                          'pa': None if pa is None else '{0:.7e}'.format(pa),
                          'reff': None if reff is None else '{0:.7e}'.format(reff)}

        # Write the configuration file
        with _ofile.open('w') as f:
            f.write('# Auto-generated configuration file\n')
            f.write('# {0}\n'.format(time.strftime("%a %d %b %Y %H:%M:%S",time.localtime())))
            f.write('\n')
            cfg.write(f)


class MaNGAAnalysisPlan(AnalysisPlan):
    """
    Define the analysis plan and output paths for MaNGA data.

    This class only redefines the output paths from the base class.
    """
    def __init__(self, planlist, cube=None, analysis_path=None):
#        if cube is None:
#            raise ValueError('For a MaNGA analysis plan, you must provide the cube because it '
#                             'is used to define the output paths.')
        self.cube = cube
        if self.cube is None:
            warnings.warn('Instantiation of the MaNGA analysis plan did not include the cube.  '
                          'Directories will be *undefined*, meaning the code may fault.')
        if analysis_path is None:
            analysis_path = dap_analysis_path(drpver=cube.drpver)
        super().__init__(planlist, cube=cube, analysis_path=analysis_path)

    def common_path(self):
        """
        Return the path for data common to all plans.

        Args:
            cube (:class:`~mangadap.datacube.datacube.DataCube`, optional):
                Cube being analyzed.  Passed for cube-specific path
                specification.  Not used by this base class.

        Returns:
            `Path`_: Path object for the "common" output
        """
        if self.cube is None:
            raise ValueError('Path undefined because cube being processed was not provided '
                             'when MaNGAAnalysisPlan was instantiated!')
        return self.analysis_path / 'common' / str(self.cube.plate) / str(self.cube.ifudesign)

    def method_path(self, plan_index=0, qa=False, ref=False):
        """
        Return the path for method-specific output.

        Args:
            cube (:class:`~mangadap.datacube.datacube.DataCube`, optional):
                Cube being analyzed.  Passed for cube-specific path
                specification.  Not used by this base class.
            plan_index (:obj:`int`, optional):
                The index of the plan.  This is used to select the 'key' of the
                analysis plan being used, which is used as the subdirectory for
                the output.
            qa (:obj:`bool`, optional):
                Flag that the output is specifically quality assessment plots.
            ref (:obj:`bool`, optional):
                Flag that the output is specifically reference data.

        Returns:
            `Path`_: Path object for the method-specific output

        Raises:
            ValueError:
                Raised if the plan index is invalid or if both qa and ref are true.

        """
        if self.cube is None:
            raise ValueError('Path undefined because cube being processed was not provided '
                             'when MaNGAAnalysisPlan was instantiated!')
        if plan_index < 0 or plan_index >= self.nplans:
            raise ValueError(f'Invalid index ({plan_index}); 0 <= index < {self.nplans}.')
        if qa and ref:
            raise ValueError('Cannot provide path for both qa and ref directory.  Pick one.')

        root = self.analysis_path / self.plan[self.plan_keys[plan_index]]['key'] \
                    / str(self.cube.plate) / str(self.cube.ifudesign)
        if not qa and not ref:
            return root
        if qa:
            return root / 'qa'
        return root / 'ref'

