"""
Base class for a datacube

----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

# TODO: Pilfer the pypeit.DataContainer for this.

# TODO: List the required metadata somewhere in here?

# TODO: Force the datacube data arrays to be read-only?

import warnings

from IPython import embed

import numpy
import warnings

from scipy import sparse, interpolate

from astropy.io import fits


from ..util.bitmask import BitMask
from ..util.mapping import permute_wcs_axes
from ..util.covariance import Covariance
from ..util.pixelmask import SpectralPixelMask
from ..util.sampling import angstroms_per_pixel
from ..util.geometry import polygon_area

class DataCube:
    r"""
    Base container class for a rectilinear datacube.

    Datacubes have three axes: two spatial coordinates and one spectral. The
    datacubes are expected to be rectilinear; i.e., each spatial position has
    the same spectral range and each wavelength channel has the same spatial
    coordinates.

    On input, the ordering of the three dimensions can be arbitrary; however,
    the axes are re-ordered for the internal attributes such that wavelengths
    are ordered along the last axis, with the first two axes being the
    spatial coordinates. Nominally the spatial coordinates are ordered
    predominantly coincident with right-ascension (E toward smaller pixels
    numbers in axis 0) and declination (N toward larger pixel numbers in axis
    1); however, rotation of the datacube spatial coordinates can be non-zero
    with respect to the celestial coordinates, as given by the provided
    world-coordinate system.

    The wavelength vector applicable to all spatial positions can either be
    provided directly or constructed from the provided WCS. Any directly
    provided ``wave`` vector takes precedence over the WCS.

    Spatial covariance/correlation can be provided via the ``covar`` keyword
    argument. If provided as a single array, the provided array is used to
    define the *correlation* matrix and is assumed to be identical for all
    wavelength channels. More options are available if the input is provided
    as a :class:`~mangadap.util.covariance.Covariance` object. The ordering
    of the covariance/correlation matrix is expected to provide the
    correlation between the row-major flattened version of the spatial
    coordinates. That is, for a datacube with spatial shape
    :math:`(N_x,N_y)`, the covariance/correlation value at location
    :math:`i,j` is the correlation between pixels at 2D locations
    :math:`(i_x,i_y)` and :math:`(j_x,j_y)`, where

    .. math::

        \begin{array}{rcl}
        i_x & = & \lfloor i/N_y \rfloor, \\
        i_y & = & i - i_x N_y, \\
        j_x & = & \lfloor j/N_y \rfloor, {\rm and} \\
        j_y & = & j - j_x N_y.
        \end{array}

    If the values of ``axes`` indicates that the provided flux array should
    have it's spatial axes transposed (i.e. ``axes[0] > axes[1]``), the
    provided covariance/correlation data is appropriately restructured to
    correspond to the wavelength channels in :attr:`flux` (again assuming a
    flattened row-major memory block); see
    :func:`~mangadap.util.covariance.Covariance.transpose_raw_shape`.
    
    Derived classes should, in particular, provide the read methods.
    The critical components of the derived classes are:

    .. todo::

        Fill this in.

    Args:
        flux (`numpy.ndarray`_):
            Data array with the flux as a function of spatial and
            spectral position. Shape should be :math:`(N_x, N_y,
            N_\lambda)`.
        wave (`numpy.ndarray`_, optional):
            The wavelength vector associated with each spatial
            position. Shape is :math:`(N_\lambda,)`. If None, ``wcs``
            must be provided such that the wavelength vector can be
            constructed.
        ivar (`numpy.ndarray`_, optional):
            The inverse variance of the flux array. Shape must be the
            same as ``flux``.  If None, errors are ignored.
        mask (`numpy.ndarray`_, optional):
            Mask array for the flux measurements. Can be a boolean
            mask (False=unmasked,good; True=masked,bad) or an integer
            array associated with the provided
            :class:`mangadap.util.bitmask.BitMask` object (see
            ``bitmask``). Shape must be the same as ``flux``.
        bitmask (:class:`mangadap.util.bitmask.BitMask`, optional):
            Object used to select and toggle masked pixels from
            :attr:`mask`. If None, any provided ``mask`` must be a
            boolean array or can be successfully converted to one.
        sres (`numpy.ndarray`_, optional):
            The spectral resolution, :math:`R = \lambda /
            \Delta\lambda`. Can be a single vector with shape
            :math:`(N_\lambda,)`, if the spectral resolution is
            independent of spatial position, or a 3D datacube with a
            spatially and spectrally dependent resolution with a
            shape that matches ``flux``. If None, spectral resolution
            is ignored.
        covar (array-like, :class:`mangadap.util.covariance.Covariance`, optional):
            The spatial covariance in one or more wavelength
            channels. If None, no covariance is assumed. See usage
            and interpretation in class description.
        axes (:obj:`list`, optional):
            The axes with the :math:`x`, :math:`y`, and
            :math:`\lambda` in the provided flux array. For example,
            ``[2,1,0]`` means the spectra are organized along the
            first axis.
        wcs (`astropy.wcs.WCS`_, optional):
            World-coordinate system for the 3D datacube. If None,
            spatial coordinates are set to be centered with a
            pixelscale in arcsec (see ``pixelscale``).
        pixelscale (:obj:`int`, :obj:`float`, optional):
            Spatial extent in arcsec of each spaxel in the datacube,
            assumed to be square. Superseded by ``wcs``, if the
            latter is provided. If the ``wcs`` is not provided and
            the pixelscale is not provided upon instantiation, it is
            set to unity.
        log (:obj:`bool`, optional):
            Flag that the datacube spectral pixels are binned
            logarithmically in wavelength.
        meta (:obj:`dict`, optional):
            A free-form dictionary used to hold metadata relevant to
            the datacube. Metadata required by analysis modules are
            indicated where relevant. If None, :attr:`meta` is
            instantiated as an empty dictionary.
        prihdr (`astropy.io.fits.Header`_, optional):
            Primary header read from datacube fits file. If None,
            instantiated as an empty `astropy.io.fits.Header`_.
        fluxhdr (`astropy.io.fits.Header`_, optional):
            Header specifically for the flux extension of the
            datacube fits file. If None, set to be a copy of the
            primary header.

    Raises:
        ValueError:
            Raised if the wavelength vector cannot be produced, if
            the WCS has the wrong dimensionality, or if there are any
            shape mismatches between the input data arrays.
        TypeError:
            Raised if the input metadata are not provided as a
            dictionary.

    Attributes:
        original_axes (`numpy.ndarray`_):
            Original axis order.
        shape (:obj:`tuple`):
            Datacube shape.
        spatial_shape (:obj:`tuple`):
            Shape of the datacube spatial axes .
        nwave (:obj:`int`):
            Number of wavelength channels.
        spatial_index (`numpy.ndarray`_):
            Array of tuples with the spatial indices of each spectrum
            in the flattened datacube.
        wcs (`astropy.wcs.WCS`_):
            Datacube world-coordinate system.  Can be None.
        pixelscale (:obj:`float`):
            Spatial extent in arcsec of each spaxel in the datacube,
            assumed to be square. Only used if :attr:`wcs` is None.
        bitmask (:class:`mangadap.util.bitmask.BitMask`):
            Object used to select and toggle masked pixels from
            :attr:`mask`.  Can be None.
        log (:obj:`bool`):
            Flag that the datacube spectral pixels are binned
            logarithmically in wavelength.
        meta (:obj:`dict`):
            A free-form dictionary used to hold meta data relevant to
            the datacube. Metadata required by analysis modules are
            indicated where relevant. If no metadata has been
            defined, :attr:`meta` is instantiated as an empty
            dictionary.
        wave (`numpy.ndarray`_):
            Wavelength vector applicable to all spatial positions.
        flux (`numpy.ndarray`_):
            Datacube flux array.
        ivar (`numpy.ndarray`_):
            Datacube flux inverse variance.  Can be None.
        mask (`numpy.ndarray`_):
            Datacube mask.  Can be None.
        sres (`numpy.ndarray`_):
            Datacube spectral resolution.  Can be None.
        covar (:class:`mangadap.util.covariance.Covariance`):
            Datacube spatial covariance.  Can be None.
        prihdr (`astropy.io.fits.Header`_):
            Primary header for the datacube. If not provided on
            instantiation, set to an empty `astropy.io.fits.Header`_.
        fluxhdr (`astropy.io.fits.Header`_):
            Header specifically for the flux array. If not provided
            on instantiation, set to be a copy of :attr:`prihdr`.
        rss (:class:`mangadap.spectra.rowstackedspectra.RowStackedSpectra`):
            The source row-stacked spectra used to build the
            datacube.
        sigma_rho (:obj:`float`):
            The :math:`\sigma_{\rho}` of the Gaussian function used
            to approximate the trend of the correlation coefficient
            with spaxel separation. Used to construct the approximate
            correlation matrix (see
            :func:`approximate_correlation_matrix`).
        correl_rlim (:obj:`float`):
            The limiting radius of the image reconstruction
            (Gaussian) kernel in arcseconds. Used to construct the
            approximate correlation matrix (see
            :func:`approximate_correlation_matrix`).
        approx_correl (:class:`~mangadap.util.covariance.Covariance`):
            Approximate correlation matrix; see
            :func:`approximate_correlation_matrix`.
    """
    # TODO: Add reconstructed PSF?
    def __init__(self, flux, wave=None, ivar=None, mask=None, bitmask=None, sres=None, covar=None,
                 axes=[0,1,2], wcs=None, pixelscale=None, log=True, meta=None, prihdr=None,
                 fluxhdr=None):

        if wcs is None and wave is None:
            raise ValueError('Must either provide a single wavelength vector or a WCS that can '
                             'be used to construct it.')
        if wcs is not None and wcs.wcs.naxis != 3:
            raise ValueError('Provided WCS object must define the coordinate system for each of '
                             'the three datacube axes.')
        self.wcs = None if wcs is None else permute_wcs_axes(wcs, axes)
        self.pixelscale = self._get_pixelscale() if self.wcs is not None and pixelscale is None \
                                else (1 if pixelscale is None else float(pixelscale))

        self.original_axes = numpy.atleast_1d(axes).copy()
        # Re-order so that axes are x, y, lambda
        self.flux = flux.transpose(self.original_axes)
        self.shape = self.flux.shape
        self.spatial_shape = self.shape[:-1]
        self.nwave = self.shape[-1]

        # TODO: Not sure this is useful
        self.spatial_index = numpy.array([(i,j) for i,j in zip(*numpy.unravel_index(
                                                numpy.arange(self.nspec), self.spatial_shape))])

        self.meta = {} if meta is None else meta
        if not isinstance(self.meta, dict):
            raise TypeError('Metadata must be provided as a dictionary.')

        self.wave = None if wave is None else numpy.atleast_1d(wave)
        if self.wave is None:
            self.wave = self._get_wavelength_vector()
        if self.wave.shape != (self.nwave,):
            raise ValueError('Wavelength vector shape is incorrect ')
        self.log = log
        # TODO: Check the above against the wavelength vector and/or
        # set the value automatically based on the input.
        
        self.ivar = None 
        if ivar is not None:
            self.ivar = ivar.transpose(self.original_axes)
            if self.ivar.shape != self.shape:
                raise ValueError('Inverse variance array has incorrect shape.')

        self.mask = None
        if mask is not None:
            self.mask = mask.transpose(self.original_axes)
            if self.mask.shape != self.shape:
                raise ValueError('Mask array has incorrect shape.')
        self.bitmask = None
        if bitmask is not None:
            if self.mask is None:
                warnings.warn('Bitmask provided but no mask data available.')
            if isinstance(bitmask, BitMask):
                self.bitmask = bitmask
            else:
                warnings.warn('Bitmasks must be BitMask objects; provided bitmask is ignored.')
        if self.bitmask is None:
            # Ensure mask is a boolean array
            self.mask = self.mask.astype(bool)

        self.redux_bitmask = None
        self.redux_qual_key = None
        self.redux_qual_flag = None
        
        self.sres = None
        if sres is not None:
            self.sres = numpy.tile(sres, self.spatial_shape+(1,)) \
                            if sres.ndim == 1 else sres.transpose(self.original_axes)
            if self.sres.shape != self.shape:
                raise ValueError('Spectral resolution array has incorrect shape.')

        self.covar = None
        if covar is not None:
            spatial_transpose = self.original_axes[0] > self.original_axes[1]
            raw_shape = self.spatial_shape[::-1] if spatial_transpose else self.spatial_shape
            self.covar = covar if isinstance(covar, Covariance) \
                            else Covariance.from_array(covar, raw_shape=raw_shape)
            if spatial_transpose:
                self.covar = self.covar.transpose_raw_shape()

        # Allocate attributes for primary and flux array fits headers
        self.prihdr = fits.Header() if prihdr is None else prihdr
        self.fluxhdr = self.prihdr.copy() if fluxhdr is None else fluxhdr

        ## KHRR added this
        if not ('OBJRA' in self.prihdr):
            prihdr['OBJRA'] = self.meta['OBJRA']
            prihdr['OBJDEC'] = self.meta['OBJDEC']


        # Allow for a RowStackedSpectrum counterpart
        self.rss = None

        # For the approximate covariance matrix calculations
        self.sigma_rho = None
        self.correl_rlim = None
        self.approx_correl = None

    @classmethod
    def from_config(cls, cfgfile, **kwargs):
        """
        Construct a datacube object using a configuration file.

        The method is undefined for the base class. If called, an
        exception is raised. Derived classes must override this
        method to allow for a configuration file to be used to
        instantiate the relevant :class:`DataCube` subclass.

        Args:
            cfgfile (:obj:`str`):
                Configuration file. See `configparser.ConfigParser`_.
            **kwargs:
                Any other keyword arguments that invoke optional
                instantiation methods. Note that these arguments will
                *never* be used in a command-line level execution of
                the DAP. They should only be available for custom
                scripts.
        """
        raise NotImplementedError('No from_config method is available for {0}.'.format(
                                  self.__class__.__name__))

    # TODO: write a from_rss classmethod

    def _get_pixelscale(self):
        """
        Measure the pixel scale using the WCS.

        The method uses the coordinates of 4 adjacent pixels to
        compute the pixel area. Assuming the pixels are square and
        that the pixel scale is constant in square arcsec over the
        full image, the square-root of the area is the pixel scale.

        The method will fault if :attr:`wcs` is None.

        Returns:
            :obj:`float`: The estimated pixel scale.
        """
        # TODO: Instead get the mean over the full image?
        coo = numpy.array([[1,1,2,2], [1,2,2,1], [1,1,1,1]]).T
        x, y, _ = self.wcs.all_pix2world(coo, 1).T
        x = (x - x[0])*numpy.cos(numpy.radians(y[0]))
        return numpy.sqrt(polygon_area(x, y))*3600

    def _get_wavelength_vector(self):
        """
        Use the `astropy.wcs.WCS`_ attribute (:attr:`wcs`) to
        generate the datacube wavelength vector.

        :attr:`wcs` cannot be None and the wavelength coordinate
        system must be defined along its third axis.

        Returns:
            `numpy.ndarray`_: Vector with wavelengths along the third
            axis defined by :attr:`wcs`.

        Raises:
            ValueError:
                Raised if :attr:`wcs` or :attr:`nwave` is not defined.
        """
        if self.nwave is None:
            raise ValueError('Length of the spectral axis (nwave) must be defined.')
        if self.wcs is None:
            raise ValueError('World coordinate system required to construct wavelength vector.')
        coo = numpy.array([numpy.ones(self.nwave), numpy.ones(self.nwave),
                           numpy.arange(self.nwave)+1]).T
        return self.wcs.all_pix2world(coo, 1)[:,2]*self.wcs.wcs.cunit[2].to('angstrom')

    @property
    def nspec(self):
        """Number of spectra in the datacube."""
        return numpy.prod(self.spatial_shape)

    @staticmethod
    def propagate_flags():
        """
        Flags that should be propagated from the observed data to the analyzed data.

        The base class returns ``None``.
        """
        return None

    @staticmethod
    def do_not_use_flags():
        """
        Return the maskbit names that should not be used.

        The base class returns ``None``.
        """
        return None

    @staticmethod
    def do_not_fit_flags():
        """
        Return the maskbit names that should not be fit.

        The base class returns ``None``.
        """
        return None

    @staticmethod
    def do_not_stack_flags():
        """
        Return the maskbit names that should not be stacked.
        
        The base class returns ``None``.
        """
        return None

    def metakeys(self):
        """Get a :obj:`list` of the keys in the datacube metadata."""
        return list(self.meta.keys())

    # TODO: Add a getitem method that returns the datacube flux?

    @property
    def can_compute_covariance(self):
        """
        Determine if the object can be used to compute the spatial
        covariance.

        If :attr:`covar` is currently not defined, the method tries
        to load the row-stacked spectra used to build the datacube;
        see :func:`load_rss`. If that is successful or if
        :attr:`covar` is already defined, the method returns True. If
        :attr:`covar` is None and the row-stacked spectra cannot be
        loaded, the method returns False.

        Returns:
            :obj:`bool`: Flag that the object can be used to
            calculate the spatial covariance.
        """
        if self.covar is None:
            self.load_rss()
        if self.covar is None and self.rss is None:
            return False
        return True

    def load_rss(self, **kwargs):
        """
        Try to load the source row-stacked spectra for this datacube.

        This method is undefined in the base class, and simply passes.

        Derived classes should override this method if it's possible
        to load the row-stacked spectra. If the load is successful,
        :attr:`rss` should no longer be None.
        """
        pass

    def copy_to_array(self, attr='flux', waverange=None, nbins=None, select_bins=None,
                      missing_bins=None, unique_bins=None):
        r"""
        Return a copy of the selected data array with a flattened
        spatial axis.

        The array size is always :math:`N_{\rm spec} \times N_{\rm
        wavelength}`. The spatial positions within the original
        datacube for each spectrum are given by tuples in
        :attr:`spatial_index`.

        See :func:`copy_to_masked_array` for argument descriptions.

        .. warning::
            Any masking is ignored in this function call. To
            incorporate the mask use :func:`copy_to_masked_array`.

        Returns:
            `numpy.ndarray`_: A 2D array with a copy of the data from the
            selected attribute.
        """
        masked_data = self.copy_to_masked_array(attr=attr, use_mask=False, waverange=waverange,
                                                nbins=nbins, select_bins=select_bins,
                                                missing_bins=missing_bins, unique_bins=unique_bins)
        # For this approach, the wavelengths masked should be
        # *identical* for all spectra
        nwave = numpy.sum(numpy.invert(numpy.ma.getmaskarray(masked_data)), axis=1)
        # No masking should be present except for the wavelength range
        # meaning that the number of unmasked pixels at all spatial
        # positions should be the same.
        if numpy.any(nwave != nwave[0]):
            raise ValueError('Masking in copy_to_mask should only for the wavelength range.')
        if numpy.all(nwave == 0):
            raise ValueError('Full wavelength range has been masked!')
        # Compression returns a flattened array, so it needs to be
        # reshaped for the new number of (unmasked) wavelength channels
        return masked_data.compressed().reshape(-1,nwave[0])

    def copy_to_masked_array(self, attr='flux', use_mask=True, flag=None, waverange=None,
                             nbins=None, select_bins=None, missing_bins=None, unique_bins=None):
        r"""
        Return a copy of the selected data array as a masked array
        with a flattened spatial axis.

        The array size is always :math:`N_{\rm spec} \times N_{\rm
        wavelength}`. The spatial positions within the original
        datacube for each spectrum are given by tuples in
        :attr:`spatial_index`.

        This is functionally identical to :func:`copy_to_array`,
        except the output format is a `numpy.ma.MaskedArray`_. The
        pixels that are considered to be masked can be specified
        using the `flag` option.
        
        Args:
            attr (:obj:`str`, optional):
                The attribute for the returned array. Can be 'flux',
                'ivar', 'sres', or 'mask'. For 'mask', you're likely
                better off using :func:`copy_to_array`. Strings are
                always set to be lower-case, so capitalization
                shouldn't matter.
            use_mask (:obj:`bool`, optional):
                Use the internal mask to mask the data. This is
                largely here to allow for :func:`copy_to_array` to
                wrap this function while not applying the internal
                mask.
            waverange (array-like, optional):
                Two-element array with the first and last wavelength
                to include in the computation. Default is to use the
                full wavelength range.
            flag (:obj:`str`, :obj:`list`, optional):
                One or more bitmask flags used to select bad pixels.
                The names *must* be a valid bit name as defined by
                :attr:`bitmask` (see
                :class:`mangadap.util.bitmask.BitMask`). If not
                provided, *any* non-zero mask bit is omitted.
            nbins (:obj:`int`, optional):
                The total number of defined bins. Default is to find
                the maximum number in the unique_bins list.
            select_bins (array-like, optional):
                A boolean array selecting spectra in the flattened
                cube (with shape :math:`N_{\rm spec} \times
                N_\lambda`) to return.
            missing_bins (:obj:`list`, optional):
                A list of bin numbers that are missing from the
                output selection and that will be replaced with
                masked spectra. If specified, must also provide
                ``unique_bins``.
            unique_bins (array-like, optional):
                The indices of the bins that have valid spectra. The
                length of this array should match the number of
                selected spectra from ``select_bins``.

        Returns:
            `numpy.ma.MaskedArray`_: A 2D array with a copy of the
            data from the selected extension, masked where
            :attr:`mask` is either True or with the selected bitmask
            flags.

        Raises:
            AttributeError:
                Raised if ``attr`` is not a valid attribute.
            ValueError:
                Raised if the unique bins are not specified, but the
                missing bins are, or if the number of unique bins
                does not match the number of spectra in the datacube.
        """
        # Make sure masks and flags make sense
        if flag is not None and self.bitmask is None:
            warnings.warn('No bitmask defined.  Input flags ignored.')

        nspec = self.nspec

        # Create the wavelength mask (will be all true if
        # waverange=None)
        mask = SpectralPixelMask(waverange=waverange).boolean(self.wave, nspec=nspec)

        # Add in any masked data
        if use_mask:
            mask |= self.mask.reshape(nspec,-1) if self.bitmask is None \
                    else self.bitmask.flagged(self.mask.reshape(nspec,-1), flag=flag)

        # Create the output MaskedArray
        # TODO: Handle when attr is None
        a = numpy.ma.MaskedArray(getattr(self, attr.lower()).reshape(nspec,-1), mask=mask)

        # Apply any bin selection
        if select_bins is not None:
            a = a[select_bins,:]

        # No missing bins are specified, so return the array
        if missing_bins is None or len(missing_bins) == 0:
            return a

        # Missing bins have been specified, so the unique_bins also need
        # to be specified
        if unique_bins is None:
            raise ValueError('Must specify the unique bins if missing bins are specified.')
        if len(unique_bins) != a.shape[0]:
            raise ValueError('Number of unique bins does not match the number of spectra.')

        # If the number of bins hasn't been provided, construct it based
        # on the unique bins specified
        _nbins = max(unique_bins)+1 if nbins is None else nbins

        # Create the array with the right shape and mask everything
        _a = numpy.ma.MaskedArray(numpy.zeros((_nbins, a.shape[-1]), dtype=a.dtype),
                                  mask=numpy.ones((_nbins, a.shape[-1]), dtype=bool))
        # Copy the data to the correct bins, which also unmasks these
        # pixels
        _a[unique_bins,:] = a
        # Return the masked array; "missing bins" are fully masked
        return _a

    def interpolate_to_match(self, func, fill_value=0.0):
        r"""
        Interpolate a function to match the datacube wavelength sampling.

        Args:
            func (`numpy.ndarray`_):
                Function to linear interpolate. Shape is
                :math:`(N_{\rm wave}, 2)`, where the first column has
                the wavelength and the second has the function to
                interpolate. If None, simply returns a unity vector
                of the correct length.
            fill_value (:obj:`float`, optional):
                The value to use for spectral regions not sampled by
                the provided function.

        Returns:
            `numpy.ndarray`_: The function interpolated (**not**
            resampled) to match :attr:`wave`. Any spectral regions
            not within the wavelength range of the provided function
            is set to 0.
        """
        if func is None:
            return numpy.ones(self.nwave, dtype=float)
        return interpolate.interp1d(func[:,0], func[:,1], bounds_error=False,
                                    fill_value=fill_value, assume_sorted=True)(self.wave)

    def mean_sky_coordinates(self, center_coo=None):
        r"""
        Compute the mean sky coordinates for each spectrum.

        If the WCS is available, the coordinates can be returned in
        RA and declination (degrees). If center coordinates are
        provided (see ``center_coo``), however, the coordinates are
        offset set as follows:

        .. math::

            x &= (\alpha - \alpha_0) \cos \delta_0 \\
            y &= (\delta - \delta_0),

        where :math:`(\alpha_0, \delta_0)` are the provided
        coordinates.

        If the WCS is not available, the returned coordinates are in
        arcsec from the center of the image (regardless of the value
        of ``center_coo``) determined using the pixelscale. At least
        in this case, the coordinates are assumed to be relative to
        the pixel center (not, e.g., its edge).

        Args:
            center_coo (:obj:`tuple`, optional):
                A two-tuple with the coordinates in right-ascension
                and declination for the coordinate-frame origin. If
                None, no offset is performed.

        Returns:
            :obj:`tuple`: Two `numpy.ndarray`_ objects with the RA
            and declination of each pixel in degrees, or its offset
            from the center in arcseconds. In both cases the shape of
            the returned arrays matches the spatial shape of the
            datacube.
        """
        i, j = numpy.meshgrid(numpy.arange(self.spatial_shape[0]),
                              numpy.arange(self.spatial_shape[1]), indexing='ij')
        if self.wcs is None:
            return (i-self.spatial_shape[0]//2+(self.spatial_shape[0]%2)*0.5)*self.pixelscale, \
                   (j-self.spatial_shape[1]//2+(self.spatial_shape[1]%2)*0.5)*self.pixelscale

        # Generate pixel coordinates to pass to the WCS; just use the
        # first wavelength channel.
        coo = numpy.array([i.ravel()+1, j.ravel()+1, numpy.ones(self.nspec)]).T
        x, y, _ = [c.reshape(self.spatial_shape) for c in self.wcs.all_pix2world(coo, 1).T]
        if center_coo is None:
            # Not offsetting, so we're done
            return x, y

        # TODO: Use astropy.coordinates.SkyOffset 
        # Offset and return
        if len(center_coo) != 2:
            raise ValueError('Provided coordinates are expected to be a single x,y pair.')
        return (x-center_coo[0]) * numpy.cos(numpy.radians(center_coo[1])) * 3600., \
                    (y-center_coo[1]) * 3600.

    def binned_on_sky_area(self, bin_indx):
        r"""
        Compute the on-sky area of a set of binned spectra.

        For each bin, this is just the number of spaxels in the bin
        multiplied by the spaxel area.

        Args:
            bin_indx (array-like):
                An array with size :math:`N_{\rm spec}` that gives
                which spaxels were included in each bin. Valid bins
                have indices of :math:`\geq 0`.

        Returns:
            :obj:`tuple`: Two `numpy.ndarray`_ objects are returned.
            The first has the unique (non-negative) bin indices, and
            the second provides the on-sky area of that bin.
        """
        unique_bins, bin_count = numpy.unique(bin_indx, return_counts=True)
        indx = unique_bins > -1
        return unique_bins[indx], (bin_count[indx]*numpy.square(self.pixelscale)).astype(float)

    def central_wavelength(self, waverange=None, response_func=None, per_pixel=True, flag=None,
                           fluxwgt=False):
        """
        Determine the mean central wavelength for all spectra under
        various conditions.

        The wavelength channel is set to be the center of the
        bandpass selected by ``waverange``, weighted by the response
        function and the flux (if requested/provided). The mask is
        also incorporated in these calculations. By default (i.e., no
        wavelength limits or weighting) the wavelength channel is
        just the central wavelength of the full spectral range. (Note
        that if the spectra are binned logarithmically, this isn't
        necessarily the central wavelength *channel*.)

        Args:
            waverange (array-like, optional):
                Starting and ending wavelength over which to
                calculate the statistics. Default is to use the full
                wavelength range.
            response_func (array-like, optional):
                A two-column array with the wavelength and
                transmission of a broad-band response function to use
                as a weighting function for the calculation.
            per_pixel (:obj:`bool`, optional):
                When providing a response function, base the
                calculation on per pixel measurements, instead of per
                angstrom. Set to False for a per-angstrom
                calculation.
            flag (:obj:`str`, :obj:`list`, optional):
                One or more flag names that are considered when
                deciding if a pixel should be masked. The names
                *must* be a valid bit name as defined by
                :attr:`bitmask`. If :attr:`bitmask` is None, these
                are ignored.
            fluxwgt (:obj:`bool`, optional):
                Flag to weight by the flux when determining the mean
                coordinates.

        Returns:
            :obj:`float`: The mean central wavelength of all spectra.
        """
        if waverange is None and response_func is None and not fluxwgt:
            return (self.wave[0] + self.wave[-1])/2.

        flux = self.copy_to_masked_array(waverange=waverange, flag=flag)
        dw = numpy.ones(self.nwave, dtype=float) if per_pixel \
                else angstroms_per_pixel(self.wave, log=self.log)
        _response_func = self.interpolate_to_match(response_func)

        if fluxwgt:
            norm = numpy.ma.sum(flux*_response_func[None,:]*dw[None,:], axis=1)
            cen_wave = numpy.ma.sum(flux*self.wave[None,:]*_response_func[None,:]*dw[None,:],
                                    axis=1) / norm
            return float(numpy.ma.mean(cen_wave))

        norm = numpy.ma.sum(numpy.invert(numpy.ma.getmaskarray(flux))
                            * _response_func[None,:] * dw[None,:], axis=1)
        cen_wave = numpy.ma.sum(numpy.invert(numpy.ma.getmaskarray(flux)) * self.wave[None,:]
                                * _response_func[None,:] * dw[None,:], axis=1) / norm
        return float(numpy.ma.mean(cen_wave))

    def flux_stats(self, waverange=None, response_func=None, per_pixel=True, flag=None):
        r"""
        Compute the mean flux, propagated error in the mean flux, and
        mean S/N over the specified wavelength range.

        If the wavelength range is not specified, the quantities are
        calculated over the full spectral range.

        Args:
            waverange (array-like, optional):
                Starting and ending wavelength over which to
                calculate the statistics. Default is to use the full
                wavelength range.
            response_func (array-like, optional):
                A two-column array with the wavelength and
                transmission of a broad-band response function to use
                as a weighting function for the calculation.
            per_pixel (:obj:`bool`, optional):
                When providing a response function, continue to
                calculate the statistics per pixel. Set to False for
                a per-angstrom calculation.
            flag (:obj:`str`, :obj:`list`, optional):
                One or more flag names that are considered when
                deciding if a pixel should be masked. The names
                *must* be a valid bit name as defined by
                :attr:`bitmask`.

        Returns:
            `numpy.ndarray`_: Three objects are returned: the mean
            flux, the propagated variance in the mean flux, and the
            mean S/N. The shape of each is the same as the shape of a
            single wavelength channel in the datacube (i.e.,
            :attr:`spatial_shape`).
    
        Raises:
            ValueError:
                Raised of a provided wavelength range object does not
                have two elements or if a response function is
                provided and has an incorrect shape.
        """
        if waverange is not None and len(waverange) != 2:
            raise ValueError('Provided wavelength range must be a two-element vector.')
        if response_func is not None:
            if len(response_func.shape) != 2:
                raise ValueError('Response function object must be two dimensional.')
            if response_func.shape[1] != 2:
                raise ValueError('Response function object must only have two columns.')

        # Grab the masked arrays
        flux = self.copy_to_masked_array(waverange=waverange, flag=flag)
        ivar = self.copy_to_masked_array(attr='ivar', waverange=waverange, flag=flag)
        snr = flux*numpy.ma.sqrt(ivar)

        # Set the response function
        dw = numpy.ones(self.nwave, dtype=float) if per_pixel \
                else angstroms_per_pixel(self.wave, log=self.log)
        _response_func = self.interpolate_to_match(response_func)

        # Calculate the statistics and return
        response_integral = numpy.sum(numpy.invert(numpy.ma.getmaskarray(flux))
                                        * (_response_func*dw)[None,:], axis=1)
        signal = numpy.ma.divide(numpy.ma.sum(flux*(_response_func*dw)[None,:], axis=1),
                                 response_integral).reshape(self.spatial_shape)
        variance = numpy.ma.divide(numpy.ma.sum(numpy.ma.power(ivar, -1.) \
                                    * (_response_func*dw)[None,:], axis=1),
                                   response_integral).reshape(self.spatial_shape)
        snr = numpy.ma.divide(numpy.ma.sum(snr*(_response_func*dw)[None,:], axis=1),
                              response_integral).reshape(self.spatial_shape)
        return signal, variance, snr

    def covariance_matrix(self, channel, **kwargs):
        """
        Construct the formal covariance matrix for the provided
        channel.

        This is a simple wrapper for a call to
        :func:`mangadap.spectra.rowstackedspectra.RowStackedSpectra.covariance_matrix`
        executed using :attr:`rss`, which cannot be None. See that
        method for the argument description.

        .. warning::

            The provided rectifications parameters should be the same
            as used to construct the datacube. If not, the calculated
            covariance matrix will not be correct.

        """
        if self.rss is None:
            raise ValueError('RowStackedSpectra object not available for this datacube!')
        return self.rss.covariance_matrix(channel, **kwargs)

    def covariance_cube(self, **kwargs):
        """
        Construct the formal covariance matrix for one or more
        channels.

        This is a simple wrapper for a call to
        :func:`mangadap.spectra.rowstackedspectra.RowStackedSpectra.covariance_cube`
        executed using :attr:`rss`, which cannot be None. See that
        method for the argument description.

        .. warning::

            The provided rectifications parameters should be the same
            as used to construct the datacube. If not, the calculated
            covariance matrix will not be correct.

        """
        if self.rss is None:
            raise ValueError('RowStackedSpectra object not available for this datacube!')
        return self.rss.covariance_cube(**kwargs)

    def approximate_correlation_matrix(self, sigma_rho, rlim, rho_tol=None, redo=False):
        r"""
        Construct a correlation matrix with correlation coefficients
        that follow a Gaussian in 2D pixel separation.

        This method constructs a correlation matrix with correlation
        coefficients defined as

        .. math::

            \rho_{ij} = \exp(-D_{ij}^2 / 2\sigma_{\rho}^2)

        where :math:`D_{ij}` is the distance between two spaxels in
        the spatial dimension of the datacube (in the number spaxels,
        not arcsec). Any pixels with :math:`D_{ij} > R_{\rm lim}` is
        set to zero, where :math:`R_{\rm lim}` is the limiting radius
        of the kernel used in the datacube rectification; this is
        provided here as ``rlim`` in arcsec, which is converted to
        spaxels using :attr:`pixelscale`.

        We found in Westfall et al. (2019, AJ, 158, 231) that this is
        a reasonable approximation for the formal covariance that
        results from Shepard's rectification method.

        There is an unknown relation between the dispersion of the
        kernel used by Shepard's method and the value of
        :math:`\sigma_{\rho}`. Tests show that it is tantalizingly
        close to :math:`\sigma_{\rho} = \sqrt{2}\sigma`, where
        :math:`\sigma` is in pixels instead of arcsec; however, a
        formal derivation of this hasn't been done and is complicated
        by the focal-plane sampling of the row-stacked spectra.

        Within the limits of how the focal-plane sampling changes
        with wavelength, we found the value of :math:`\sigma_{\rho}`
        varies little with wavelength in MaNGA. Here, we assume
        :math:`\sigma_{\rho}` is fully wavelength independent.
        Therefore, once this method is run once, it doesn't need to
        be run again for different wavelength channels. To force the
        correlation matrix to be recreated, use ``redo`` or change
        the provided :math:`\sigma_{\rho}`.

        Args:
            sigma_rho (:obj:`float`):
                The :math:`\sigma_{\rho}` of the Gaussian function
                used to approximate the trend of the correlation
                coefficient with spaxel separation.
            rlim (:obj:`float`):
                The limiting radius of the image reconstruction
                (Gaussian) kernel in arcseconds.
            rho_tol (:obj:`float`, optional):
                Any correlation coefficient less than this is assumed
                to be equivalent to (and set to) 0.
            redo (:obj:`bool`, optional):
                Force the recalculation of the cube dimensions if
                they are already defined and :math:`\sigma_{\rho}` has
                not changed.

        Returns:
            :class:`mangadap.util.covariance.Covariance`: Correlation
            matrix
        """
        if self.approx_correl is not None \
                and any([sigma_rho == this for this in [None, self.sigma_rho]]) \
                and any([rlim == this for this in [None, self.correl_rlim]]):
            return self.approx_correl

        if sigma_rho is None or rlim is None:
            raise ValueError('Must provide sigma_rho and rlim for approximate correlation '
                             'matrix calculation.')

        # Get the full covariance grid
        ii, jj = map(lambda x: x.ravel(),
                     numpy.meshgrid(numpy.arange(self.nspec), numpy.arange(self.nspec),
                                    indexing='ij'))

        # Convert covariance pixels to spatial pixels along both dimensions
        i_i, i_j = numpy.unravel_index(ii, self.spatial_shape)
        j_i, j_j = numpy.unravel_index(jj, self.spatial_shape)

        # Get the (square of the) distances from each spaxel to every
        # other spaxel
        dij = numpy.square(j_i-i_i) + numpy.square(j_j-i_j)
        indx = dij <= numpy.square(2*rlim/self.pixelscale)

        # Calculate the correlation coefficient
        ii = ii[indx]
        jj = jj[indx]
        rhoij = numpy.exp(-dij[indx]/numpy.square(sigma_rho)/2)

        if rho_tol is not None:
            indx = rhoij > rho_tol
            ii = ii[indx]
            jj = jj[indx]
            rhoij = rhoij[indx]

        # Construct the Covariance object and save the input
        self.sigma_rho = sigma_rho
        self.correl_rlim = rlim
        self.approx_correl = Covariance(sparse.coo_matrix((rhoij, (ii,jj)),
                                                          shape=(self.nspec,self.nspec)).tocsr(),
                                        impose_triu=True, correlation=True,
                                        raw_shape=self.spatial_shape)
        return self.approx_correl

    def approximate_covariance_matrix(self, channel, sigma_rho=None, rlim=None, rho_tol=None,
                                      csr=False, quiet=False):
        r"""
        Return an approximate calculation of the covariance matrix
        assuming
        
        .. math::

            C_{ij} = \rho_{ij}(V_{i} V_{j})^{1/2}

        where :math:`\rho_{ij}` is approximated by
        :func:`approximate_correlation_matrix` and
        :math:`V_i\equivC_{ii}` are the variances provided by the
        inverse of :attr:`ivar`.

        The method first calculates :math:`\rho_{ij}` if it hasn't
        been yet or the provided ``sigma_rho`` and/or ``rlim`` values
        are different to previous calls, which builds
        :attr:`approx_correl`. If :attr:`approx_correl` is not yet
        constructed, ``sigma_rho`` and ``rlim`` must be provided.

        The returned covariance matrix is the correlation matrix
        rescaled by the variance provided by :attr:`ivar` for the
        specified channel.
        
        Args:
            channel (:obj:`int`):
                Index of the spectral channel for which to calculate
                the covariance matrix.
            sigma_rho (:obj:`float`, optional):
                The :math:`\sigma_{\rho}` of the Gaussian function
                used to approximate the trend of the correlation
                coefficient with spaxel separation.
            rlim (:obj:`float`, optional):
                The limiting radius of the image reconstruction
                (Gaussian) kernel in arcseconds.
            rho_tol (:obj:`float`, optional):
                Any correlation coefficient less than this is assumed
                to be equivalent to (and set to) 0.
            csr (:obj:`bool`, optional):
                Instead of returning a
                :class:`mangadap.util.covariance.Covariance` object,
                return the covariance matrix as a
                `scipy.sparse.csr_matrix`_ object. Primarily used by
                :func:`approximate_covariance_cube` for collating the
                covariance matrix of each wavelength channel before
                combining them into a single
                :class:`mangadap.util.covariance.Covariance` object.
            quiet (:obj:`bool`, optional):
                Suppress terminal output

        Returns:
            :class:`mangadap.util.covariance.Covariance`,
            `scipy.sparse.csr_matrix`_: The approximate covariance
            matrix for the designated wavelength channel. The return
            type depends on `csr`.
        """
        self.approximate_correlation_matrix(sigma_rho, rlim, rho_tol=rho_tol)
        var = numpy.ma.power(self.ivar[...,channel], -1).filled(0.0).ravel()
        covar = self.approx_correl.apply_new_variance(var)
        covar.revert_correlation()
        return covar.full() if csr else covar

    def approximate_covariance_cube(self, channels=None, sigma_rho=None, rlim=None, rho_tol=None,
                                    csr=False, quiet=False):
        r"""
        Return the approximate covariance matrices for many
        wavelength channels.

        This is a simple wrapper for
        :func:`approximate_covariance_matrix` that iteratively builds
        the covariance matrix for each wavelength channel.

        If :attr:`approx_correl` is not yet constructed (see
        :func:`approximate_correlation_matrix`), ``sigma_rho`` and
        ``rlim`` must be provided.
        
        Args:
            channels (:obj:`int`, array-like, optional):
                Indices of the spectral channels for which to
                calculate the covariance matrix. If None, the
                covariance matrix is calculated for *all* channels.
            sigma_rho (:obj:`float`, optional):
                The :math:`\sigma_{\rho}` of the Gaussian function
                used to approximate the trend of the correlation
                coefficient with spaxel separation.
            rlim (:obj:`float`, optional):
                The limiting radius of the image reconstruction
                (Gaussian) kernel in arcseconds.
            rho_tol (:obj:`float`, optional):
                Any correlation coefficient less than this is assumed
                to be equivalent to (and set to) 0.
            csr (:obj:`bool`, optional):
                Instead of returning a
                :class:`mangadap.util.covariance.Covariance` object,
                return the covariance matrix as a
                `scipy.sparse.csr_matrix`_ object. Primarily used by
                :func:`approximate_covariance_cube` for collating the
                covariance matrix of each wavelength channel before
                combining them into a single
                :class:`mangadap.util.covariance.Covariance` object.
            quiet (:obj:`bool`, optional):
                Suppress terminal output

        Returns:
            :class:`mangadap.util.covariance.Covariance`,
            `numpy.ndarray`_: If ``csr`` is True, the returned object
            is an `numpy.ndarray`_ of `scipy.sparse.csr_matrix`_
            types.
        """
        # Set the channels
        _channels = numpy.arange(self.nwave) if channels is None \
                        else numpy.atleast_1d(channels)
        nc = len(_channels)

        # Build and return the covariance matrices
        CovCube = numpy.empty(nc, dtype=sparse.csr.csr_matrix)
        for i in range(nc):
            if not quiet:
                print('Calculating covariance matrix: {0}/{1}'.format(i+1,nc), end='\r')
            CovCube[i] = self.approximate_covariance_matrix(_channels[i], sigma_rho=sigma_rho,
                                                            rlim=rlim, rho_tol=rho_tol, csr=True,
                                                            quiet=True)
        if not quiet:
            print('Calculating covariance matrix: {0}/{0}'.format(nc))
        return Covariance(CovCube, input_indx=_channels) if not csr else CovCube

