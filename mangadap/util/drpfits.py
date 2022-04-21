# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
r"""
Defines a class used to interface with files produced in the 3D phase of
the MaNGA Data Reduction Pipeline (DRP).

.. todo::

    - Calculation in :func:`DRPFits._cube_dimensions` will only be
      correct if the WCS coordinates have no rotation.
    - Further optimize calculation of transfer matrix
    - Make DRP file class flexible to linear or log-linear wavelength
      sampling?  Incorporate into MODE?
    - Reconstructed broad-band images and PSFs are *not* restructured in
      the CUBE files!  This is why the are transposed in
      :func:`mangadap.drpfits.DRPFits.gri_composite`.
    - Image reconstruction has transpose sense wrt DRP output!
    - Add logging

    - Need to be clear about which functions use the RSS spectra to
      create CUBE related data, like the covariance matrix and
      instrumental dispersion calculations.

    - Computing the approximate covariance cube is currently not
      possible with only the CUBE on disk.  There's a logic problem that
      needs to be fixed.

----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import time
import os
import warnings

from configparser import ConfigParser

import numpy

from .parser import DefaultConfig
from .bitmask import BitMask
from .constants import DAPConstants
from .filter import interpolate_masked_vector
from ..config import defaults


class DRPFitsBitMask(BitMask):
    r"""
    Structure with the DRP mask bits.

    The defined bits are listed at :ref:`metadatamodel-drp3pixmask`.
    """
    def __init__(self, sdss_maskbits=None, mode='CUBE'):
        DRPFits.check_mode(mode)
        _sdss_maskbits = defaults.sdss_maskbits_file() if sdss_maskbits is None else sdss_maskbits
        tmp = BitMask.from_par_file(_sdss_maskbits, 'MANGA_DRP3PIXMASK' if mode == 'CUBE' else 
                                                    'MANGA_DRP2PIXMASK')
        keys, descr = tmp._init_objs()
        super(DRPFitsBitMask, self).__init__(keys, descr=descr)


class DRPQuality3DBitMask(BitMask):
    r"""
    Structure with the definition of the DRP3QUAL mask bits.

    The defined bits are listed at :ref:`metadatamodel-drp3qual`.
    """
    def __init__(self, sdss_maskbits=None):
        _sdss_maskbits = defaults.sdss_maskbits_file() if sdss_maskbits is None else sdss_maskbits
        tmp = BitMask.from_par_file(_sdss_maskbits, 'MANGA_DRP3QUAL')
        keys, descr = tmp._init_objs()
        super(DRPQuality3DBitMask, self).__init__(keys, descr=descr)


class DRPFits:
    r"""

    Provides methods common to both
    :class:`mangadap.datacube.manga.MaNGADataCube` and
    :class:`mangadap.spectra.manga.MaNGARSS`, and not general enough
    for :class:`mangadap.datacube.datacube.DataCube` and
    :class:`mangadap.spectra.rowstackedspectra.RowStackedSpectra`.
    
    """
    def __init__(self, plate, ifudesign, mode, log=True, drpver=None, redux_path=None,
                 directory_path=None):
        DRPFits.check_mode(mode)
        self.plate = plate
        self.ifudesign = ifudesign
        self.mode = mode
        self.samp = 'LOG' if log else 'LIN'
        self.directory_path, self.file_name \
                = DRPFits.default_paths(self.plate, self.ifudesign, self.mode, log=log,
                                        drpver=drpver, redux_path=redux_path,
                                        directory_path=directory_path)
        self.redux_bitmask = DRPQuality3DBitMask()
        self.redux_qual_key = 'DRP3QUAL'
        self.redux_qual_flag = 'CRITICAL'

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
        options = DRPFits.mode_options()
        if mode not in options:
            raise ValueError('Unknown mode {0}.  Must be in: {1}'.format(mode, options))

    @staticmethod
    def mode_options():
        """
        Return the allowed modes.

        Returns:
            :obj:`list`: List of the allowed DRP fits file modes.
        """
        return ['CUBE', 'RSS']

    @staticmethod
    def sampling_options():
        """
        Return the allowed wavelength sampling modes.

        Returns:
            :obj:`list`: List of the allowed DRP fits wavelength
            sampling modes.
        """
        return ['LIN', 'LOG']

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
            warnings.warn('{0} is not a viable spectral resolution extension.'.format(ext))
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
        _ext = DRPFits.spectral_resolution_extension(hdu, ext=ext)
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
    def build_file_name(plate, ifudesign, mode, log=True):
        """
        Return the name of the DRP row-stacked spectra file.

        Args:
            plate (:obj:`int`):
                Plate number
            ifudesign (:obj:`int`):
                IFU design
            mode (:obj:`str`):
                Format mode (CUBE or RSS)
            log (:obj:`bool`, optional):
                Use the spectra that are logarithmically sampled in
                wavelength. If False, sampling is linear in
                wavelength.

        Returns:
            :obj:`str`: The relevant file name.
        """
        DRPFits.check_mode(mode)
        return '{0}.fits.gz'.format(defaults.manga_fits_root(plate, ifudesign,
                                             '{0}{1}'.format('LOG' if log else 'LIN', mode)))

    @staticmethod
    def default_paths(plate, ifudesign, mode, log=True, drpver=None, redux_path=None,
                      directory_path=None):
        """
        Return the primary directory and file name with the DRP fits
        LOG-binned file.

        Args:
            plate (:obj:`int`):
                Plate number
            ifudesign (:obj:`int`):
                IFU design
            mode (:obj:`str`):
                Format mode (CUBE or RSS)
            log (:obj:`bool`, optional):
                Use the spectra that are logarithmically sampled in
                wavelength. If False, sampling is linear in
                wavelength.
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

        Returns:
            :obj:`tuple`: Two strings with the default path to and
            name of the DRP data file.
        """
        _directory_path = defaults.drp_directory_path(plate, drpver=drpver,
                                                      redux_path=redux_path) \
                                if directory_path is None else directory_path
        return _directory_path, DRPFits.build_file_name(plate, ifudesign, mode, log=log)

    def file_path(self):
        """Return the full path to the DRP file"""
        return os.path.join(self.directory_path, self.file_name)

    def finding_chart_path(self):
        """
        Return the full path to the PNG finding chart for the targetted
        object.

        .. todo::
            - Move this to the defaults file?
        """
        return os.path.join(self.directory_path, 'images', str(self.ifudesign)+'.png')

    @classmethod
    def from_config(cls, cfgfile, drpver=None, redux_path=None, directory_path=None):
        """
        Construct a :class:`mangadap.datacube.manga.MaNGADataCube` or
        :class:`mangadap.spectra.manga.MaNGARSS` object from a
        configuration file.

        Using the data read from the configuration file, the method
        instantiates the class using
        :func:`mangadap.spectra.manga.MaNGARSS.from_plateifu` or
        :func:`mangadap.datacube.manga.MaNGADataCube.from_plateifu`.
        This method will therefore fault for this base class!

        The format of the configuration file is:

        .. todo::

            Fill this in.

        Args:
            cfgfile (:obj:`str`):
                Configuration file
            drpver (:obj:`str`, optional):
                DRP version, which is used to define the default DRP
                redux path. Overrides any value in the configuration
                file.
            redux_path (:obj:`str`, optional):
                The path to the top level directory containing the
                DRP output files for a given DRP version. Overrides
                any value in the configuration file.
            directory_path (:obj:`str`, optional):
                The exact path to the DRP file. Providing this
                ignores anything provided for ``drpver`` or
                ``redux_path``. Overrides any value in the
                configuration file.
        """
        # Read the configuration file
        cfg = DefaultConfig(cfgfile, interpolate=True)

        # Set the attributes, forcing a known type
        plate = cfg.getint('plate')
        ifu = cfg.getint('ifu')
        if plate is None or ifu is None:
            raise ValueError('Configuration file must define the plate and IFU.')
        log = cfg.getbool('log', default=True)

        # Overwrite what's in the file with the method keyword arguments
        _drpver = cfg.get('drpver') if drpver is None else drpver
        _redux_path = cfg.get('redux_path') if redux_path is None else redux_path
        _directory_path = cfg.get('directory_path') if directory_path is None else directory_path

        # Get the other possible keywords
        # TODO: Come up with a better way to do this
        kwargs = {}
        kwargs['sres_ext'] = cfg.get('sres_ext')
        kwargs['sres_fill'] = cfg.getbool('sres_fill', default=True)
        kwargs['covar_ext'] = cfg.get('covar_ext')
        kwargs['z'] = cfg.getfloat('z')
        kwargs['vdisp'] = cfg.getfloat('vdisp')
        kwargs['ell'] = cfg.getfloat('ell')
        kwargs['pa'] = cfg.getfloat('pa')
        kwargs['reff'] = cfg.getfloat('reff')

        return cls.from_plateifu(plate, ifu, log=log, drpver=_drpver, redux_path=_redux_path,
                                 directory_path=_directory_path, **kwargs)

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
        if os.path.exists(ofile) and not overwrite:
            raise FileExistsError('Configuration file already exists; to overwrite, set '
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
        with open(ofile, 'w') as f:
            f.write('# Auto-generated configuration file\n')
            f.write('# {0}\n'.format(time.strftime("%a %d %b %Y %H:%M:%S",time.localtime())))
            f.write('\n')
            cfg.write(f)

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
