"""
Implement the derived class for MaNGA datacubes.

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

from IPython import embed

import numpy

from astropy.io import fits
from astropy.wcs import WCS

from ..config import defaults
from ..util.drpfits import DRPFits, DRPFitsBitMask
from ..util.parser import DefaultConfig
from ..util.constants import DAPConstants
from ..util.covariance import Covariance
from ..spectra import MaNGARSS
from .datacube import DataCube

class MaNGADataCube(DRPFits, DataCube):
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
        ifile (:obj:`str`):
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
            :func:`~mangadap.util.drpfits.DRPFits.spectral_resolution`.
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

        if not os.path.isfile(ifile):
            raise FileNotFoundError('File does not exist: {0}'.format(ifile))

        # Parse the relevant information from the filename
        directory_path = os.path.split(os.path.abspath(ifile))[0]
        plate, ifudesign, log = ifile.split('-')[1:4]
        # Instantiate the DRPFits base
        DRPFits.__init__(self, int(plate), int(ifudesign), 'CUBE', log='LOG' in log,
                         directory_path=directory_path)

        # Collect the metadata into a dictionary
        meta = {}
        meta['z'] = z
        meta['vdisp'] = vdisp
        meta['ell'] = ell
        meta['pa'] = pa
        meta['reff'] = reff

        # Try to define the BitMask object
        try:
            bitmask = DRPFitsBitMask(mode='CUBE')
        except:
            warnings.warn('Unable to define bit mask for MaNGA datacube.  Can only distinguish '
                          ' between masked (values greater than 0) and unmasked (values of 0).')
            bitmask = None

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
            self.sres_ext, sres = DRPFits.spectral_resolution(hdu, ext=sres_ext, fill=sres_fill)
            self.sres_fill = sres_fill
            sres = sres.filled(0.0)

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
                              pixelscale=0.5, log=self.samp == 'LOG', meta=meta, prihdr=prihdr,
                              fluxhdr=fluxhdr)
        print('Reading MaNGA datacube data ... DONE')

        # Try to use the header to set the DRP version
        self.drpver = self.prihdr['VERSDRP3'] if 'VERSDRP3' in self.prihdr \
                        else defaults.drp_version()
        # Reduction path is always set to the default. A warning is
        # thrown if the default reduction path is not the same as the
        # path expected for the file.
        self.redux_path = defaults.drp_redux_path(drpver=self.drpver)
        if self.directory_path != defaults.drp_directory_path(self.plate, drpver=self.drpver,
                                                              redux_path=self.redux_path):
            warnings.warn('Default reduction path does not match file path.  May not be able to '
                          'find paired RSS file if requested.')

    # TODO: Include a class method that instantiates from (or wraps a Marvin Cube)

    @staticmethod
    def build_file_name(plate, ifudesign, log=True):
        """
        Return the name of the DRP datacube file.

        This is a simple wrapper for
        :func:`mangadap.util.drpfits.DRPFits.build_file_name`,
        specific to the datacube.

        Args:
            plate (:obj:`int`):
                Plate number
            ifudesign (:obj:`int`):
                IFU design
            log (:obj:`bool`, optional):
                Use the spectra that are logarithmically sampled in
                wavelength. If False, sampling is linear in
                wavelength.

        Returns:
            :obj:`str`: The relevant file name.
        """
        return DRPFits.build_file_name(plate, ifudesign, 'CUBE', log=log)

    @staticmethod
    def default_paths(plate, ifudesign, log=True, drpver=None, redux_path=None,
                      directory_path=None):
        """
        Construct the default path and file name with the MaNGA
        datacube.

        This is a simple wrapper for
        :func:`mangadap.util.drpfits.DRPFits.default_paths`, specific
        to the datacube files.

        Args:
            plate (:obj:`int`):
                Plate number
            ifudesign (:obj:`int`):
                IFU design
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
        return DRPFits.default_paths(plate, ifudesign, 'CUBE', log=log, drpver=drpver,
                                     redux_path=redux_path, directory_path=directory_path)

    @classmethod
    def from_plateifu(cls, plate, ifudesign, log=True, drpver=None, redux_path=None,
                      directory_path=None, **kwargs):
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
            **kwargs:
                Keyword arguments passed directly to the primary
                instantiation method; see
                :class:`mangadap.datacube.manga.MaNGADataCube` or
                :class:`mangadap.spectra.manga.MaNGARSS`.
        """
        path, file_name = cls.default_paths(int(plate), int(ifudesign), log=log, drpver=drpver,
                                            redux_path=redux_path, directory_path=directory_path)
        return cls(os.path.join(path, file_name), **kwargs)

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

        rss_file = DRPFits.build_file_name(self.plate, self.ifudesign, 'RSS', log=self.log)
        rss_file_path = os.path.join(self.directory_path, rss_file)
        if not os.path.isfile(rss_file_path):
            # Not in this directory.  Check the nominal directory
            warnings.warn('Could not find: {0}'.format(rss_file_path))
            rss_file_path = os.path.join(
                    defaults.drp_directory_path(self.plate, drpver=self.drpver,
                                                redux_path=self.redux_path), rss_file)
        if not os.path.isfile(rss_file_path):
            warnings.warn('Could not find: {0}'.format(rss_file_path))
            raise FileNotFoundError('Could not find RSS file.')

        self.rss = MaNGARSS(rss_file_path, sres_ext=self.sres_ext, sres_fill=self.sres_fill)

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
        :func:`mangadap.config.defaults.regrid_rlim`.
        """
        _rlim = defaults.regrid_rlim() if rlim is None else rlim
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
        :func:`mangadap.config.defaults.regrid_rlim`.
        """
        _rlim = defaults.regrid_rlim() if rlim is None else rlim
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
        :func:`mangadap.config.defaults.regrid_rlim`.
        """
        _rlim = defaults.regrid_rlim() if rlim is None else rlim
        return super().approximate_covariance_cube(channels=channels, sigma_rho=sigma_rho,
                                                   rlim=_rlim, rho_tol=rho_tol, csr=csr,
                                                   quiet=quiet)

