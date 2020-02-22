"""
Implement the derived class for MaNGA row-stacked spectra.

----

.. include license and copyright
.. include:: ../copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../links.rst
"""
import os
import warnings

from IPython import embed

import numpy

from astropy.io import fits

from ..config import defaults
from ..drpfits import DRPFitsBitMask
from ..util.constants import DAPConstants
from ..util.filter import interpolate_masked_vector
from .rowstackedspectra import RowStackedSpectra

class MaNGARSS(RowStackedSpectra):
    r"""
    Container class for MaNGA row-stacked spectra.

    For additional description and attributes, see
    :class:`mangadap.rowstackedspectra.rowstackedspectra.RowStackedSpectra`.

    See :func:`from_plateifu` to instantiate using the plate and ifu
    numbers.

    Args:
        ifile (:obj:`str`):
            The file with the MaNGA row-stacked spectra. The name is
            expected to follow the nominal naming convention, which
            is used to parse the plate, ifudesign, and whether or not
            the spectra are binned logarithmically in wavelength.
        sres_ext (:obj:`str`, optional):
            The extension to use when constructing the spectral
            resolution vectors. See :func:`spectral_resolution`.
        sres_fill (:obj:`bool`, optional):
            Fill masked values by interpolation. Default is to leave
            masked pixels in returned array.
    """
    def __init__(self, ifile, sres_ext=None, sres_fill=True):

        if not os.path.isfile(ifile):
            raise FileNotFoundError('File does not exist: {0}'.format(ifile))

        # Parse the relevant information from the filename
        self.directory_path, self.file = os.path.split(os.path.abspath(ifile)) 
        self.plate, self.ifudesign, log = self.file.split('-')[1:4]
        self.plate = int(self.plate)
        self.ifudesign = int(self.ifudesign)
        log = 'LOG' in log

        # Try to define the BitMask object
        try:
            bitmask = DRPFitsBitMask(mode='RSS')
        except:
            warnings.warn('Unable to define bit mask for MaNGA row-stacked spectra.  Can only '
                          'distinguish between masked (values greater than 0) and unmasked '
                          '(values of 0).')
            bitmask = None

        # Open the file and initialize the base class
        with fits.open(ifile) as hdu:
            print('Reading MaNGA row-stacked spectra data ...', end='\r')
            self.prihdr = hdu[0].header
            self.fluxhdr = hdu['FLUX'].header
            self.sres_ext, sres = MaNGARSS.spectral_resolution(hdu, ext=sres_ext, fill=sres_fill)
            self.sres_fill = sres_fill
            sres = sres.filled(0.0)

            # The negative for XPOS is because the data in that
            # extension has negative offsets going toward increasing
            # RA. Opposite of the RowStackedSpectra convention. The
            # area of the fiber aperture is normalized to pi arcsec^2
            # by the DRP.
            super(MaNGARSS, self).__init__(hdu['WAVE'].data, hdu['FLUX'].data,
                                           ivar=hdu['IVAR'].data, mask=hdu['MASK'].data,
                                           bitmask=bitmask, sres=sres,
                                           xpos=-hdu['XPOS'].data, ypos=hdu['YPOS'].data,
                                           area=numpy.pi, log=log)
        print('Reading MaNGA row-stacked spectra data ... DONE')

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

    @staticmethod
    def spectral_resolution_extension(hdu, ext=None):
        """
        Determine the spectral resolution channel to use.

        Precedence follows this order: ``PREDISP``, ``PRESPECRES``,
        ``DISP``, ``SPECRES``.

        Args:
            hdu (`astropy.io.fits.HDUList`):
                The opened MaNGA DRP file.
            ext (:obj:`str`, optional):
                Specify the extension with the spectral estimate to
                use. Should be in None, ``PREDISP``, ``PRESPECRES``,
                ``DISP``, or ``SPECRES``. The default is None, which
                means it will return the extension found first in the
                order above. None is returned if none of the
                extensions are present.

        Returns:
            :obj:`str`: The name of the preferred extension to use.
        """
        available = [h.name for h in hdu if h.name in ['PREDISP', 'DISP', 'PRESPECRES', 'SPECRES']]
        _ext = ext
        if ext is None:
            _ext = 'PREDISP'
            if _ext not in available:
                _ext = 'PRESPECRES'
            if _ext not in available:
                _ext = 'DISP'
            if _ext not in available:
                _ext = 'SPECRES'
        return None if _ext not in available else _ext

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
                resolution instead of a per spectrum array.  When using
                the `SPECRES` extension, this just returns the vector
                provided by the DRP file; when using the `DISP`
                extension, this performs a masked median across the
                array and then interpolates any wavelengths that were
                masked in all vectors.

        Returns:
            :obj:`tuple`: Returns a :obj:`str` with the name of the
            extension used for the spectral resolution measurements
            and a `numpy.ma.MaskedArray`_ with the spectral
            resolution data. Even if interpolated such that there
            should be no masked values, the function returns a masked
            array. Array contains the spectral resolution (:math:`R =
            \lambda/\Delta\lambda`) pulled from the DRP file.
        """
        # Determine which spectral resolution element to use
        _ext = MaNGARSS.spectral_resolution_extension(hdu, ext=ext)
        # If no valid extension, raise an exception
        if ext is None and _ext is None:
            raise ValueError('No valid spectral resolution extension.')
        if ext is not None and _ext is None:
            raise ValueError('No extension: {0}'.format(ext))
            
        # Build the spectral resolution vectors
        if 'SPECRES' in _ext:
            sres = numpy.ma.MaskedArray(hdu[_ext].data, mask=numpy.invert(hdu[_ext].data > 0))
            if fill:
                sres = numpy.ma.MaskedArray(interpolate_masked_vector(sres))
            if not median:
                nspec = hdu['FLUX'].data.shape[0]
                sres = numpy.ma.tile(sres, (nspec,1))
            return _ext, sres

        # Otherwise dealing with the DISP cube
        sres = numpy.ma.MaskedArray(hdu[_ext].data)
        # Mask any non-positive value
        sres[numpy.invert(sres > 0)] = numpy.ma.masked
        # Convert from sigma in angstroms to spectral resolution
        # (based on FWHM)
        sres = numpy.ma.divide(hdu['WAVE'].data[None,:], sres) / DAPConstants.sig2fwhm
        # Interpolate over any masked values
        if fill:
            sres = numpy.ma.MaskedArray(numpy.ma.apply_along_axis(interpolate_masked_vector, 1,
                                                                  sres))
        if median:
            sres = numpy.ma.median(sres, axis=0)
        return _ext, sres

    @classmethod
    def from_plateifu(cls, plate, ifudesign, log=True, drpver=None, redux_path=None,
                      directory_path=None, **kwargs):
        """
        Construct a MaNGA row-stacked spectra object based on its
        plate-ifu designation.

        The provided plate, ifudesign, and boolean for the log
        binning are used to construct the name of the file. The DRP
        version and reduction and directory paths are used to
        construct the directory with the file.

        Args:
            plate (:obj:`int`):
                Plate number
            ifudesign (:obj:`int`):
                IFU design
            log (:obj:`bool`, optional):
                Use the spectra that are logarithmically binned in
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
                instantiation method; see :class:`MaNGARSS`.
        """
        # Set the attributes, forcing a known type
        _plate = int(plate)
        _ifudesign = int(ifudesign)

        # Setup the directory path.
        if directory_path is None:
            directory_path = defaults.drp_directory_path(_plate, drpver=drpver,
                                                         redux_path=redux_path)
        return cls(os.path.join(directory_path, MaNGARSS.build_file_name(_plate, _ifudesign,
                                                                         log=log)),
                   **kwargs)

    @staticmethod
    def build_file_name(plate, ifudesign, log=True):
        """
        Return the name of the DRP row-stacked spectra file.

        Args:
            plate (:obj:`int`):
                Plate number
            ifudesign (:obj:`int`):
                IFU design
            log (:obj:`bool`, optional):
                Use the spectra that are logarithmically binned in
                wavelength.

        Returns:
            :obj:`str`: The relevant file name.
        """
        mode = '{0}RSS'.format('LOG' if log else 'LIN')
        return '{0}.fits.gz'.format(defaults.manga_fits_root(plate, ifudesign, mode))

    def file_path(self):
        """Return the full path to the DRP datacube file."""
        return os.path.join(self.directory_path,
                            MaNGARSS.build_file_name(self.plate, self.ifudesign, log=self.log))

    @staticmethod
    def do_not_fit_flags():
        """Return the maskbit names that should not be fit."""
        return ['DONOTUSE', 'FORESTAR']


    @staticmethod
    def do_not_stack_flags():
        """Return the maskbit names that should not be stacked."""
        return ['DONOTUSE', 'FORESTAR']

    @staticmethod
    def default_paths(plate, ifudesign, log=True, drpver=None, redux_path=None,
                      directory_path=None, output_file=None):
        """
        Return the primary directory and file name with the DRP fits
        LOG-binned file.

        Args:
            plate (:obj:`int`):
                Plate number
            ifudesign (:obj:`int`):
                IFU design
            log (:obj:`bool`, optional):
                Use the spectra that are logarithmically binned in
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
            output_file (:obj:`str`, optional):
                The name of the file with the DRP data. Default set
                by :func:`MaNGARSS.build_file_name`.

        Returns:
            :obj:`str`: Two strings with the path to and name of the DRP
            data file.
        """
        _directory_path = defaults.drp_directory_path(plate, drpver=drpver,
                                                      redux_path=redux_path) \
                                if directory_path is None else directory_path
        _output_file = MaNGADataCube.build_file_name(plate, ifudesign, log=log) \
                            if output_file is None else output_file
        return _directory_path, _output_file

    def pointing_offset(self):
        """
        Return the offsets in RA and DEC between the pointing
        coordinates (``IFURA``, ``IFUDEC``) and the designated object
        center coordinates (``OBJRA``, ``OBJDEC``), drawn from the
        primary header of the DRP fits file.

        Returns:
            :obj:`float`: Sky-right arcsecond offsets in RA and DEC.
        """
        return ((self.prihdr['IFURA'] - self.prihdr['OBJRA']) \
                    * numpy.cos(numpy.radians(self.prihdr['OBJDEC'])) * 3600.), \
               ((self.prihdr['IFUDEC'] - self.prihdr['OBJDEC']) * 3600.)

    def mean_sky_coordinates(self, offset=True, **kwargs):
        r"""
        Compute the weighted or unweighted mean sky coordinates for
        each spectrum.

        This is a simple wrapper for
        :func:`mangadap.spectra.rowstackedspectra.RowStackedSpectra.mean_sky_coordinates`;
        see that method for additional keyword arguments.

        .. warning::

            Flux-weighting the coordinates can produce spurious
            results in low-flux regimes.

        Args:
            offset (:obj:`bool`, optional):
                Offset the coordinates in :attr:`xpos` and
                :attr:`ypos` by the difference between the IFU and
                OBJ coordinates in the header (see
                :func:`pointing_offset`). The coordinates are
                nominally defined with respect to the IFU pointing
                center, such that this offsets the coordinates to be
                centered on the galaxy.
            **kwargs:
                Passed directly to
                :func:`mangadap.spectra.rowstackedspectra.RowStackedSpectra.mean_sky_coordinates`;

        Returns:
            :obj:`tuple`: Two `numpy.ndarray`_ objects with the
            fiducial x and y on-sky positions of each spectrum in
            arcseconds relative to a given center. The shape of the
            two arrays is :math:`(N_{\rm spec},)`.
        """
        # TODO: offset cannot be a keyword in kwargs.
        _offset = self.pointing_offset() if offset else None
        return super(MaNGARSS, self).mean_sky_coordinates(offset=_offset, **kwargs)

    def binned_on_sky_area(self, bin_indx, offset=True, **kwargs):
        r"""
        Compute the on-sky area of a set of binned spectra.

        This is a simple wrapper for
        :func:`mangadap.spectra.rowstackedspectra.RowStackedSpectra.binned_on_sky_area`
        and :func:`mean_sky_coordinates`.

        The arguments are passed directly to
        :func:`mean_sky_coordinates`, and the results are then passed
        to the base class function. The area of the fibers
        "beams" are all renormalized to :math:`\pi {\rm arcsec}^2` by
        the DRP.

        Args:
            bin_indx (array-like):
                An array with size :math:`N_{\rm spec}` that gives
                which spaxels were included in each bin. Valid bins
                have indices of :math:`\geq 0`.
            offset (:obj:`bool`, optional):
                Offset the coordinates in :attr:`xpos` and
                :attr:`ypos` by the difference between the IFU and
                OBJ coordinates in the header (see
                :func:`pointing_offset`). The coordinates are
                nominally defined with respect to the IFU pointing
                center, such that this offsets the coordinates to be
                centered on the galaxy.
            **kwargs:
                Passed directly to :func:`mean_sky_coordinates`.

        Returns:
            :obj:`tuple`: Two `numpy.ndarray`_ objects are returned.
            The first has the unique (non-negative) bin indices, and
            the second provides the on-sky area of that bin.
        """
        # TODO: May need a better way to do this. I.e., coordinates
        # will be need to figure out which spectra to bin, and this
        # recalculation of the fiducial coordinates runs the risk that
        # the binned area is not calculated using the same coordinates
        # used to determine which spectra went in each bin.
        x, y = self.mean_sky_coordinates(offset=offset, **kwargs)
        return super(MaNGARSS, self).binned_on_sky_area(bin_indx, x=x, y=y)

    @staticmethod
    def _parse_rectification_parameters(pixelscale, rlim, sigma, recenter, width_buffer):
        """
        Parse the datacube rectification parameters.

        Values that are None are set to their default methods in
        :mod:`mangadap.config.defaults`.

        Returns:
            :obj:`tuple`: The default or input values of arguments,
            in the same order.

        """
        pixelscale = defaults.cube_pixelscale() if pixelscale is None else pixelscale
        rlim = defaults.regrid_rlim() if rlim is None else rlim
        sigma = defaults.regrid_sigma() if sigma is None else sigma
        recenter = defaults.cube_recenter() if recenter is None else recenter
        width_buffer = defaults.cube_width_buffer() if width_buffer is None else width_buffer
        return pixelscale, rlim, sigma, recenter, width_buffer

    def rectification_transfer_matrix(self, channel, pixelscale=None, rlim=None, sigma=None,
                                      recenter=None, width_buffer=None, quiet=False,
                                      rej_flag='3DREJECT'):
        """
        Construct the rectification transfer matrix.

        A simple wrapper for
        :func:`mangadap.spectra.rowstackedspectra.RowStackedSpectra.rectification_transfer_matrix`
        that defaults to the nominal rectification parameters for
        MaNGA datacubes. See this base class function for the
        argument descriptions.
        """
        _pixelscale, _rlim, _sigma, _recenter, _width_buffer \
                = MaNGARSS._parse_rectification_parameters(pixelscale, rlim, sigma, recenter,
                                                           width_buffer)
        return super(MaNGARSS,
            self).rectification_transfer_matrix(channel, _pixelscale, _rlim, _sigma,
                                                recenter=_recenter, width_buffer=_width_buffer,
                                                quiet=quiet, rej_flag='3DREJECT')

    def rectify_wavelength_plane(self, channel, pixelscale=None, rlim=None, sigma=None,
                                 recenter=None, width_buffer=None, quiet=False,
                                 rej_flag='3DREJECT', return_ivar=False, return_covar=False):
        """
        Rectify a wavelength channel into a monochromatic flux map.

        A simple wrapper for
        :func:`mangadap.spectra.rowstackedspectra.RowStackedSpectra.rectify_wavelength_plane`
        that defaults to the nominal rectification parameters for
        MaNGA datacubes. See this base class function for the
        argument descriptions.
        """
        _pixelscale, _rlim, _sigma, _recenter, _width_buffer \
                = MaNGARSS._parse_rectification_parameters(pixelscale, rlim, sigma, recenter,
                                                           width_buffer)
        return super(MaNGARSS, self).rectify_wavelength_plane(channel, pixelscale=_pixelscale,
                                                              rlim=_rlim, sigma=_sigma,
                                                              recenter=_recenter,
                                                              width_buffer=_width_buffer,
                                                              quiet=quiet, rej_flag='3DREJECT',
                                                              return_ivar=return_ivar,
                                                              return_covar=return_covar)

    def covariance_matrix(self, channel, pixelscale=None, rlim=None, sigma=None, recenter=None,
                          width_buffer=None, csr=False, quiet=False, rej_flag='3DREJECT'):
        """
        Construct the covariance matrix for a given wavelength channel.

        A simple wrapper for
        :func:`mangadap.spectra.rowstackedspectra.RowStackedSpectra.covariance_matrix`
        that defaults to the nominal rectification parameters for
        MaNGA datacubes. See this base class function for the
        argument descriptions.
        """
        _pixelscale, _rlim, _sigma, _recenter, _width_buffer \
                = MaNGARSS._parse_rectification_parameters(pixelscale, rlim, sigma, recenter,
                                                           width_buffer)
        return super(MaNGARSS, self).covariance_matrix(channel, pixelscale=_pixelscale, rlim=_rlim,
                                                       sigma=_sigma, recenter=_recenter,
                                                       width_buffer=_width_buffer, csr=csr,
                                                       quiet=quiet, rej_flag=rej_flag)

    def covariance_cube(self, channels=None, pixelscale=None, rlim=None, sigma=None, recenter=None,
                        width_buffer=None, csr=False, quiet=False, rej_flag='3DREJECT'):
        """
        Construct the covariance matrices for many wavelength channels.

        A simple wrapper for
        :func:`mangadap.spectra.rowstackedspectra.RowStackedSpectra.covariance_cube`
        that defaults to the nominal rectification parameters for
        MaNGA datacubes. See this base class function for the
        argument descriptions.
        """
        _pixelscale, _rlim, _sigma, _recenter, _width_buffer \
                = MaNGARSS._parse_rectification_parameters(pixelscale, rlim, sigma, recenter,
                                                           width_buffer)
        return super(MaNGARSS, self).covariance_cube(channels=channels, pixelscale=_pixelscale,
                                                     rlim=_rlim, sigma=_sigma, recenter=_recenter,
                                                     width_buffer=_width_buffer, csr=csr,
                                                     quiet=quiet, rej_flag=rej_flag)

    def instrumental_dispersion_plane(self, channel, dispersion_factor=None, pixelscale=None,
                                      rlim=None, sigma=None, recenter=None, width_buffer=None,
                                      quiet=False, rej_flag='3DREJECT'):
        """
        Construct the instrumental dispersion for a given wavelength channel.

        A simple wrapper for
        :func:`mangadap.spectra.rowstackedspectra.RowStackedSpectra.instrumental_dispersion_plane`
        that defaults to the nominal rectification parameters for
        MaNGA datacubes. See this base class function for the
        argument descriptions.
        """
        _pixelscale, _rlim, _sigma, _recenter, _width_buffer \
                = MaNGARSS._parse_rectification_parameters(pixelscale, rlim, sigma, recenter,
                                                           width_buffer)
        return super(MaNGARSS,
                    self).instrumental_dispersion_plane(channel,
                                                        dispersion_factor=dispersion_factor,
                                                        pixelscale=_pixelscale, rlim=_rlim,
                                                        sigma=_sigma, recenter=_recenter,
                                                        width_buffer=_width_buffer, quiet=quiet,
                                                        rej_flag=rej_flag)

