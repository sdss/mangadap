"""
Class defining the I/O configuration for MaNGA files in the DAP.
"""
import warnings
from pathlib import Path
from astropy.io import fits

from . import defaults
from ..util.parser import DefaultConfig
from ..util.drpbitmask import DRPQuality3DBitMask
from .constants import DAPConstants
from .filter import interpolate_masked_vector

class MaNGAConfig:

    instrument = 'manga'

    def __init__(self, plate, ifudesign, mode='CUBE', log=True, drpver=None, redux_path=None,
                 chart_path=None, directory_path=None, analysis_path=None):
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
        self.drpver = defaults.drp_version() if drpver is None else drpver
        # TODO: Depending on what is provided to the function, redux path, chart
        # path, and directory path can all be inconsistent.
        self.redux_path = defaults.drp_redux_path(drpver=self.drpver) \
                            if redux_path is None else Path(redux_path).resolve()
        self.chart_path = defaults.drp_finding_chart_path(self.plate, drpver=self.drpver,
                                                          redux_path=self.redux_path) \
                                if chart_path is None else Path(chart_path).resolve()
        self.directory_path = defaults.drp_directory_path(self.plate, drpver=self.drpver,
                                                          redux_path=self.redux_path) \
                                if directory_path is None else Path(directory_path).resolve()
        # TODO: This is only relevant the CUBE data
        if self.directory_path != defaults.drp_directory_path(self.plate, drpver=self.drpver,
                                                              redux_path=redux_path):
            warnings.warn('Paths do not match MaNGA defaults.  May not be able to find paired '
                          'CUBE & RSS files, if requested.')

        self.analysis_path = defaults.dap_analysis_path(drpver=self.drpver) \
                                if analysis_path is None else Path(analysis_path).resolve()
        # TODO: Turn these all into meta
        self.ebv_key = 'EBVGAL'
        self.qual_bitmask = DRPQuality3DBitMask()
        self.qual_key = 'DRP3QUAL'
        self.qual_flag = 'CRITICAL'


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
                           directory_path=self.directory_path, analysis_path=self.analysis_path)

    @property
    def file_name(self):
        root = defaults.manga_fits_root(self.plate, self.ifudesign,
                                        mode=f'{"LOG" if self.log else "LIN"}{self.mode}')
        return f'{root}.fits.gz'

    @property
    def file_path(self):
        return self.directory_path / self.file_name

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
    def key(self):
        """
        Unique key used to identify the input data.
        """
        return f'{self.plate}-{self.ifudesign}'

    @property
    def output_root(self):
        return f'{self.instrument}-{self.key}'

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


sort this out


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
                :func:`mangadap.config.defaults.drp_version`
            redux_path (:obj:`str`, optional):
                The path to the top level directory containing the
                DRP output files for a given DRP version. Default is
                defined by
                :func:`mangadap.config.defaults.drp_redux_path`.
            directory_path (:obj:`str`, optional):
                The exact path to the DRP file. Default is defined by
                :func:`mangadap.config.defaults.drp_directory_path`.
                Providing this ignores anything provided for
                ``drpver`` or ``redux_path``.
            overwrite (:obj:`bool`, optional):
                Overwrite any existing parameter file.
        """
        _ofile = Path(ofile).resolve()
        if _ofile.exists() and not overwrite:
            raise FileExistsError(f'{_ofile} already exists; to overwrite, set '
                                  'overwrite=True.')

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


