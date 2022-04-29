"""
Implement the derived class for MaNGA datacubes.

.. todo::

    - Further optimize calculation of transfer matrix/
    - Make DRP file class flexible to linear or log-linear wavelength
      sampling?  Incorporate into MODE?
    - Reconstructed broad-band images and PSFs are *not* restructured in
      the CUBE files!  
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
from pathlib import Path
import warnings

from IPython import embed

import numpy

from astropy.io import fits
from astropy.wcs import WCS

from ..config.manga import MaNGAConfig
from ..util.drpbitmask import DRPFitsBitMask, DRPQuality3DBitMask
from ..util.parser import DefaultConfig
from ..util.constants import DAPConstants
from ..util.covariance import Covariance
from ..spectra import MaNGARSS
from .datacube import DataCube


# NOTE: Because of inheriting from MaNGAConfig and the specific code in the
# MaNGAConfig.from_file method, beware of adding a `from_file` method to this!
# It shouldn't be a problem though, because the main instantiation method is
# effectively the "from_file" method.
class MaNGADataCube(MaNGAConfig, DataCube):
    r"""
    Container class for a MaNGA datacube.

    For additional description and attributes, see the two base
    classes.

    See :func:`from_plateifu` to instantiate using the plate and ifu
    numbers. See
    :func:`~mangadap.datacube.datacube.DataCube.from_config` to
    instantiate from a configuration file.

    Note that ``z``, ``vdisp``, ``ell``, ``pa``, and ``reff`` are
    saved to the :attr:`~mangadap.datacube.datacube.DataCube.meta`
    dictionary.

    Args:
        ifile (:obj:`str`, `Path`_):
            The file with the MaNGA datacube. The name is expected to
            follow the nominal naming convention, which is used to
            parse the plate, ifudesign, and whether or not the
            datacube has been binned logarithmically in wavelength.
        z (:obj:`float`, optional):
            Estimated bulk redshift. If None, some of the DAP
            analysis modules will fault.
        vdisp (:obj:`float`, optional):
            Estimated velocity dispersion. If None, some of the DAP
            analysis modules will assume an initial guess of 100
            km/s.
        ell (:obj:`float`, optional):
            Characteristic isophotal ellipticity (1-b/a). If None,
            some of the DAP modules will issue a warning and continue
            by assuming ``ell=0``.
        pa (:obj:`float`, optional):
            Characteristic isophotal position angle (through E from
            N). If None, some of the DAP modules will issue a warning
            and continue by assuming ``pa=0``.
        reff (:obj:`float`, optional):
            Effective (half-light) radius in arcsec. If None, some of
            the DAP modules will issue a warning and continue by
            assuming ``reff=1``.
        sres_ext (:obj:`str`, optional):
            The extension to use when constructing the spectral
            resolution vectors. See
            :func:`~mangadap.config.manga.MaNGAConfig.spectral_resolution`.
        sres_fill (:obj:`bool`, optional):
            Fill masked values by interpolation. Default is to leave
            masked pixels in returned array.
        covar_ext (:obj:`str`, optional):
            Extension to use as the single spatial correlation matrix
            for all wavelength channels, read from the DRP file. For
            generating the covariance matrix directly for an
            arbitrary wavelength channel using the RSS file, see
            :func:`mangadap.datacube.datacube.DataCube.covariance_matrix`.
    """
    def __init__(self, ifile, z=None, vdisp=None, ell=None, pa=None, reff=None, sres_ext=None,
                 sres_fill=True, covar_ext=None):

        # Use the configuration class method to construct the configuration
        # based on the input data file
        cfg = MaNGAConfig.from_file(ifile)
        # and then use the parent class init to set it to self.  NOTE: This is
        # round about, but the way I got it to work.  There must be a better way
        # to do this...
        MaNGAConfig.__init__(self, cfg.plate, cfg.ifudesign, mode=cfg.mode, log=cfg.log,
                             drpver=cfg.drpver, redux_path=cfg.redux_path,
                             chart_path=cfg.chart_path, directory_path=cfg.directory_path)

        # Collect the metadata into a dictionary
        meta = {}
        meta['z'] = z
        meta['vdisp'] = vdisp
        meta['ell'] = ell
        meta['pa'] = pa
        meta['reff'] = reff

        # Define the relevant BitMask object
        bitmask = DRPFitsBitMask(mode='CUBE')

        # Open the file and initialize the DataCube base class
        with fits.open(ifile) as hdu:
            # Read covariance first
            covar = None if covar_ext is None \
                        else Covariance.from_fits(hdu, ivar_ext=None, covar_ext=covar_ext,
                                                  impose_triu=True, correlation=True)

            print('Reading MaNGA datacube data ...', end='\r')
            prihdr = hdu[0].header
            fluxhdr = hdu['FLUX'].header
            # NOTE: Need to keep a log of what spectral resolution
            # vectors were used so that the same vectors are read from
            # the RSS file, if/when it is loaded.
            self.sres_ext, sres = MaNGAConfig.spectral_resolution(hdu, ext=sres_ext, fill=sres_fill)
            self.sres_fill = sres_fill
            sres = sres.filled(0.0)

            # Pull the Galactic E(B-V) from the header and add it to the metadata
            meta['ebv'] = prihdr['EBVGAL']
            # Determine if the datacube was flagged as having CRITICAL DRP failures
            meta['drpcrit'] = DRPQuality3DBitMask().flagged(prihdr['DRP3QUAL'], flag='CRITICAL')

            # NOTE: Transposes are done here because of how the data is
            # read from the fits file. The covariance is NOT transposed
            # accordingly because of how the correlation matrices are
            # stored by the DRP. Some care needs to be taken here
            # though because the transpose changes the array memory
            # storage from C to Fortran contiguous. However, everything
            # is expected to be okay because numpy array-flattening
            # always performs a C-like flattening, even if the memory
            # storage is Fortran contiguous.
            DataCube.__init__(self, hdu['FLUX'].data.T, wave=hdu['WAVE'].data,
                              ivar=hdu['IVAR'].data.T, mask=hdu['MASK'].data.T, bitmask=bitmask,
                              sres=sres.T, covar=covar, wcs=WCS(header=fluxhdr, fix=True),
                              pixelscale=0.5, log=self.log, meta=meta, prihdr=prihdr,
                              fluxhdr=fluxhdr, name=f'{self.plate}-{self.ifudesign}')
        print('Reading MaNGA datacube data ... DONE')

    @classmethod
    def from_plateifu(cls, plate, ifudesign, log=True, drpver=None, redux_path=None,
                      chart_path=None, directory_path=None, **kwargs):
        """
        Construct a :class:`mangadap.datacube.manga.MaNGADataCube`
        object based on its plate-ifu designation.

        This is a simple wrapper function that calls
        :func:`default_paths` to construct the file names, and then
        calls the class instantiation method.

        Args:
            plate (:obj:`int`):
                Plate number
            ifudesign (:obj:`int`):
                IFU design
            log (:obj:`bool`, optional):
                Use the datacube that is logarithmically binned in
                wavelength.
            drpver (:obj:`str`, optional):
                DRP version, which is used to define the default DRP
                redux path. Default is defined by
                :func:`mangadap.config.manga.drp_version`
            redux_path (:obj:`str`, optional):
                The path to the top level directory containing the
                DRP output files for a given DRP version. Default is
                defined by
                :func:`mangadap.config.manga.drp_redux_path`.
            chart_path (:obj:`str`, optional):
                The path to the directory containing the finding chart images
                for this plate.  Default is defined by
                :func:`mangadap.config.manga.drp_finding_chart_path`.
            directory_path (:obj:`str`, optional):
                The exact path to the DRP file. Default is defined by
                :func:`mangadap.config.manga.drp_directory_path`.
                Providing this ignores anything provided for
                ``drpver`` or ``redux_path``.
            **kwargs:
                Keyword arguments passed directly to the primary
                instantiation method; see
                :class:`mangadap.datacube.manga.MaNGADataCube` or
                :class:`mangadap.spectra.manga.MaNGARSS`.
        """
        # Use the input to instantiate the configuration
        cfg = MaNGAConfig(plate, ifudesign, log=log, drpver=drpver, redux_path=redux_path,
                          chart_path=chart_path, directory_path=directory_path)
        # Then use the configuration to instantiate the object
        return cls(str(cfg.file_path), **kwargs)

    @classmethod
    def from_config(cls, cfgfile):
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
        # First use the base class method to instantiate the configuration
        cfg = MaNGAConfig.from_config(cfgfile, mode='CUBE')

        # Then read the file (again...) to get the relevant keyword arguments
        full_cfg = DefaultConfig(cfgfile, interpolate=True)
        kwargs = {}
        kwargs['z'] = full_cfg.getfloat('z')
        kwargs['vdisp'] = full_cfg.getfloat('vdisp')
        kwargs['ell'] = full_cfg.getfloat('ell')
        kwargs['pa'] = full_cfg.getfloat('pa')
        kwargs['reff'] = full_cfg.getfloat('reff')
        kwargs['sres_ext'] = full_cfg.get('sres_ext')
        kwargs['sres_fill'] = full_cfg.getbool('sres_fill', default=True)
        kwargs['covar_ext'] = full_cfg.get('covar_ext')
        return cls(str(cfg.file_path), **kwargs)

    def load_rss(self, force=False):
        """
        Try to load the source row-stacked spectra for this datacube.

        If :attr:`~mangadap.datacube.datacube.DataCube.rss` is not
        None, this method does not attempt to reload it, unless
        ``force`` is True.

        The method first looks for the relevant file in the same
        directory with the datacube. If the file is not there, it
        uses the DRP version in the header of the datacube file and
        the paths defined by the keyword arguments to construct the
        nominal path to the RSS file. If the file is also not there,
        the method faults.

        Nothing is returned. If successful, the method initializes
        the row-stacked spectra object
        (:class:`mangadap.spectra.manga.MaNGARSS`) to
        :attr:`~mangadap.datacube.datacube.DataCube.rss`.

        Args:
            force (:obj:`bool`, optional):
                Reload the row-stacked spectra if
                :attr:`~mangadap.datacube.datacube.DataCube.rss` is
                not None.
        """
        if self.rss is not None and not force:
            return

        # Get the RSS counterpart
        rss_cfg = MaNGAConfig(self.plate, self.ifudesign, mode='RSS', log=self.log,
                              directory_path=self.directory_path)
        if not rss_cfg.file_path.exists():
            # Not in this directory.  Check the nominal directory
            warnings.warn(f'Could not find: {rss_file_path}.  Trying the default path.')
            rss_cfg = MaNGAConfig(self.plate, self.ifudesign, mode='RSS', log=self.log)
        if not rss_cfg.file_path.exists():
            raise FileNotFoundError(f'Could not find: {rss_file_path}')

        self.rss = MaNGARSS(str(rss_cfg.file_path), sres_ext=self.sres_ext,
                            sres_fill=self.sres_fill)

    def mean_sky_coordinates(self, center_coo=None, offset='OBJ'):
        """
        Calculate the sky coordinates for each spaxel.

        This is a simple wrapper for
        :func:`mangadap.datacube.datacube.DataCube.mean_sky_coordinates`
        that passes the relevant set of center coordinates; see that
        function for the default behavior. Any provided coordinates
        (see ``center_coo``) take precedence.

        Args:
            center_coo (:obj:`tuple`, optional):
                A two-tuple with the coordinates in right-ascension
                and declination for the coordinate-frame origin. If
                None, no offset is performed.
            offset (:obj:`str`, optional):
                The offset coordinates to use from the primary
                header. Must be either 'obj' or 'ifu'. If None, no
                header coordinates are used.

        Returns:
            :obj:`tuple`: Two `numpy.ndarray`_ objects with the RA
            and declination of each pixel in degrees, or its offset
            from the center in arcseconds. In both cases the shape of
            the returned arrays matches the spatial shape of the
            datacube.
        """
        if center_coo is not None and offset is not None:
            warnings.warn('Center coordinates provided directly.  Ignore other offset requests.')
        if offset is not None and offset.upper() not in ['OBJ', 'IFU']:
            raise ValueError('Offset must be None, obj or ifu.')
        _center_coo = (self.prihdr['{0}RA'.format(offset.upper())],
                       self.prihdr['{0}DEC'.format(offset.upper())]) \
                           if center_coo is None and offset is not None else center_coo
        return super().mean_sky_coordinates(center_coo=_center_coo)

    def approximate_correlation_matrix(self, sigma_rho=1.92, rlim=None, rho_tol=None, redo=False):
        """
        Construct an approximate correlation matrix for the datacube.

        This is a simple wrapper for
        :func:`mangadap.datacube.datacube.DataCube.approximate_correlation_matrix`
        that defaults to the expected values for a MaNGA datacube.
        See that method for a description of the arguments. The
        default ``sigma_rho`` is from Westfall et al. (2019, AJ, 158,
        231). The default ``rlim`` is provided by
        :func:`mangadap.config.manga.MaNGAConfig.regrid_rlim`.
        """
        _rlim = MaNGAConfig.regrid_rlim() if rlim is None else rlim
        return super().approximate_correlation_matrix(sigma_rho, _rlim, rho_tol=rho_tol, redo=redo)

    def approximate_covariance_matrix(self, channel, sigma_rho=1.92, rlim=None, rho_tol=None,
                                      csr=False, quiet=False):
        """
        Construct an approximate covariance matrix.
        
        This is a simple wrapper for
        :func:`mangadap.datacube.datacube.DataCube.approximate_covariance_matrix`
        that defaults to the expected values for a MaNGA datacube.
        See that method for a description of the arguments. The
        default ``sigma_rho`` is from Westfall et al. (2019, AJ, 158,
        231). The default ``rlim`` is provided by
        :func:`~mangadap.config.manga.MaNGAConfig.regrid_rlim`.
        """
        _rlim = MaNGAConfig.regrid_rlim() if rlim is None else rlim
        return super().approximate_covariance_matrix(channel, sigma_rho=sigma_rho, rlim=_rlim,
                                                     rho_tol=rho_tol, csr=csr, quiet=quiet)

    def approximate_covariance_cube(self, channels=None, sigma_rho=1.92, rlim=None, rho_tol=None,
                                    csr=False, quiet=False):
        """
        Construct approximate covariance matrices for multiple channels.
        
        This is a simple wrapper for
        :func:`mangadap.datacube.datacube.DataCube.approximate_covariance_cube`
        that defaults to the expected values for a MaNGA datacube.
        See that method for a description of the arguments. The
        default ``sigma_rho`` is from Westfall et al. (2019, AJ, 158,
        231). The default ``rlim`` is provided by
        :func:`mangadap.config.manga.MaNGAConfig.regrid_rlim`.
        """
        _rlim = MaNGAConfig.regrid_rlim() if rlim is None else rlim
        return super().approximate_covariance_cube(channels=channels, sigma_rho=sigma_rho,
                                                   rlim=_rlim, rho_tol=rho_tol, csr=csr,
                                                   quiet=quiet)

