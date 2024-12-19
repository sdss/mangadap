"""
Implement the derived class for MaNGA row-stacked spectra.

----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
from IPython import embed

import numpy

from astropy.io import fits

from ..config.manga import MaNGAConfig
from ..util.drpbitmask import DRPFitsBitMask
from ..util.parser import DefaultConfig
from .rowstackedspectra import RowStackedSpectra


class MaNGARSS(MaNGAConfig, RowStackedSpectra):
    r"""
    Container class for MaNGA row-stacked spectra.

    For additional description and attributes, see the two base
    classes.

    See :func:`from_plateifu` to instantiate using the plate and ifu
    numbers. See :func:`from_config` to instantiate from a
    configuration file.

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

        # Use the configuration class method to construct the configuration
        # based on the input data file
        cfg = MaNGAConfig.from_file(ifile)
        # and then use the parent class init to set it to self.  NOTE: This is
        # round about, but the way I got it to work.  There must be a better way
        # to do this...
        MaNGAConfig.__init__(self, cfg.plate, cfg.ifudesign, mode=cfg.mode, log=cfg.log,
                             drpver=cfg.drpver, redux_path=cfg.redux_path,
                             chart_path=cfg.chart_path, directory_path=cfg.directory_path)

        # Define the relevant BitMask object
        bitmask = DRPFitsBitMask(mode='RSS')

        # Open the file and initialize the base class
        with fits.open(ifile) as hdu:
            print('Reading MaNGA row-stacked spectra data ...', end='\r')
            prihdr = hdu[0].header
            fluxhdr = hdu['FLUX'].header
            self.sres_ext, sres = MaNGAConfig.spectral_resolution(hdu, ext=sres_ext, fill=sres_fill)
            self.sres_fill = sres_fill
            sres = sres.filled(0.0)

            # The negative for XPOS is because the data in that
            # extension has negative offsets going toward increasing
            # RA. Opposite of the RowStackedSpectra convention. The
            # area of the fiber aperture is normalized to pi arcsec^2
            # by the DRP.
            RowStackedSpectra.__init__(self, hdu['WAVE'].data, hdu['FLUX'].data,
                                       ivar=hdu['IVAR'].data, mask=hdu['MASK'].data,
                                       bitmask=bitmask, sres=sres, xpos=-hdu['XPOS'].data,
                                       ypos=hdu['YPOS'].data, area=numpy.pi,
                                       log=self.log, prihdr=prihdr, fluxhdr=fluxhdr)
        print('Reading MaNGA row-stacked spectra data ... DONE')

    @classmethod
    def from_plateifu(cls, plate, ifudesign, log=True, drpver=None, redux_path=None,
                      chart_path=None, directory_path=None, **kwargs):
        """
        Construct a :class:`mangadap.spectra.manga.MaNGARSS`
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
                Use the spectra that are logarithmically binned in
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
                :class:`mangadap.spectra.manga.MaNGARSS`.
        """
        # Use the input to instantiate the configuration
        cfg = MaNGAConfig(plate, ifudesign, mode='RSS', log=log, drpver=drpver,
                          redux_path=redux_path, chart_path=chart_path,
                          directory_path=directory_path)
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
        cfg = MaNGAConfig.from_config(cfgfile, mode='RSS')

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

        Values that are None are set to their defaults as set by
        :class:`~mangadap.config.manga.MaNGAConfig`.

        Returns:
            :obj:`tuple`: The default or input values of arguments,
            in the same order.
        """
        pixelscale = MaNGAConfig.cube_pixelscale() if pixelscale is None else pixelscale
        rlim = MaNGAConfig.regrid_rlim() if rlim is None else rlim
        sigma = MaNGAConfig.regrid_sigma() if sigma is None else sigma
        recenter = MaNGAConfig.cube_recenter() if recenter is None else recenter
        width_buffer = MaNGAConfig.cube_width_buffer() if width_buffer is None else width_buffer
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

