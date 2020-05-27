"""
Base class for row-stacked spectra

----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

# TODO: Pilfer the pypeit.DataContainer for this.

# TODO: Add meta, as in DataCube.

import warnings
from IPython import embed

import numpy
from scipy import sparse, interpolate

from astropy.io import fits

try:
    from shapely.ops import cascaded_union
    from shapely.geometry import Point
except ImportError:
    warnings.warn('Could not import shapely!', ImportWarning)

from ..util.bitmask import BitMask
from ..util.pixelmask import SpectralPixelMask
from ..util.sampling import angstroms_per_pixel
from ..util.covariance import Covariance
from ..util.constants import DAPConstants

class RowStackedSpectra:
    r"""
    Base container class for a set of 1D extracted spectra.

    Row-stacked spectra have two axes: The first is just used to
    organize the spectra, whereas the second is the spectral axis.
    That is, one can loop through the spectra by, e.g.::

        for s in spectra:
            plt.plot(s)
            ...
    
    The spectral axis is expected to be the same for all spectra.

    Currently no spectral or "spatial" (e.g. cross-talk) covariance
    is accounted for.

    .. warning::

        The class can be used to construct datacubes from row-stacked
        spectra. In these calculations, the aperture area is used to
        rescale the flux to ensure that it is conserved per unit
        area. While the class allows the aperture area to be
        fiber-specific, it does not yet allow the smoothing kernel
        used in the rectification to also be fiber-specific. Beware
        if you're trying to construct a cube using multiple fiber
        formats!
    
    Derived classes should, in particular, provide the read methods.
    The critical components of the derived classes are:

    .. todo::

        - Fill in these "critical components"
        - How do I abstract how to provide the on-sky aperture?

    Args:
        wave (`numpy.ndarray`_):
            The wavelength vector of *all* spectra. Shape is
            :math:`(N_\lambda,)`.
        flux (`numpy.ndarray`_):
            Data array with the flux for each spectrum. Shape should
            be :math:`(N_{\rm spec}, N_\lambda)`.
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
            independent of spatial position, or a 2D set of spectral
            resolution vectors with a shape that matches ``flux``. If
            None, spectral resolution is ignored.
        xpos (`numpy.ndarray`_, optional):
            The wavelength-dependent coordinate along the
            right-ascension direction of the on-sky aperture for all
            spectra. Expected to be arcseconds from some fiducial
            coordinate. Coordinates should increase from W to E.
            Shape must be the same as ``flux``.
        ypos (`numpy.ndarray`_, optional):
            The wavelength-dependent coordinate along the declination
            direction of the on-sky aperture for all spectra.
            Expected to be arcseconds from some fiducial coordinate.
            Coordinates should increase from S to N. Shape must be
            the same as ``flux``.
        area (:obj:`float`, `numpy.ndarray`_, optional):
            The fiducial area subtended by each spectral aperture,
            assumed to be wavelength independent. The aperture is
            assumed to be circular. Can provide a single value for
            all fibers or an array with shape :math:`(N_{\rm spec},)`
            that provides fiber-specific values. If None, the area is
            set to unity, which has follow-on effects on the results
            provided by some of the class methods, as described where
            relevant.
        log (:obj:`bool`, optional):
            Flag that the datacube spectral pixels are binned
            logarithmically in wavelength.
        prihdr (`astropy.io.fits.Header`_, optional):
            Primary header read from source fits file. If None,
            instantiated as an empty `astropy.io.fits.Header`_.
        fluxhdr (`astropy.io.fits.Header`_, optional):
            Header specifically for the flux extension of the source
            fits file. If None, set to be a copy of the primary
            header.

    Attributes:
        shape (:obj:`tuple`):
            Shape of the main spectral arrays.
        nwave (:obj:`int`):
            Number of wavelength channels.
        bitmask (:class:`mangadap.util.bitmask.BitMask`):
            Object used to select and toggle masked pixels from
            :attr:`mask`.  Can be None.
        log (:obj:`bool`):
            Flag that the spectral pixels are binned logarithmically
            in wavelength.
        wave (`numpy.ndarray`_):
            Wavelength vector applicable to all spectra.
        flux (`numpy.ndarray`_):
            Flux array.
        ivar (`numpy.ndarray`_):
            Flux inverse variance.  Can be None.
        mask (`numpy.ndarray`_):
            Spectral mask.  Can be None.
        sres (`numpy.ndarray`_):
            Spectral-resolution array.  Can be None.
        xpos (`numpy.ndarray`_):
            The wavelength-dependent coordinate along the
            right-ascension direction of the on-sky aperture for all
            spectra. Expected to be arcseconds from some fiducial
            coordinate.
        ypos (`numpy.ndarray`_):
            The wavelength-dependent coordinate along the declination
            direction of the on-sky aperture for all spectra.
            Expected to be arcseconds from some fiducial coordinate.
        area (`numpy.ndarray`_):
            The fiducial area subtended by each spectral aperture.
            The aperture is assumed to be circular.
        prihdr (`astropy.io.fits.Header`_):
            Primary header for the row-stacked spectra. If not
            provided on instantiation, set to an empty
            `astropy.io.fits.Header`_.
        fluxhdr (`astropy.io.fits.Header`_):
            Header specifically for the flux array. If not provided
            on instantiation, set to be a copy of :attr:`prihdr`.
    """
    def __init__(self, wave, flux, ivar=None, mask=None, bitmask=None, sres=None, xpos=None,
                 ypos=None, area=None, log=True, prihdr=None, fluxhdr=None):

        # Re-order so that axes are x, y, lambda
        self.wave = wave
        self.flux = flux
        self.shape = self.flux.shape
        self.nwave = self.shape[-1]
        if self.wave.size != self.nwave:
            raise ValueError('Wavelength vector does not match flux array.')
        self.log = log
        # TODO: Check the above against the wavelength vector and/or
        # set the value automatically based on the input.
        
        self.ivar = None 
        if ivar is not None:
            self.ivar = ivar
            if self.ivar.shape != self.shape:
                raise ValueError('Inverse variance array has incorrect shape.')

        self.mask = None
        if mask is not None:
            self.mask = mask
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

        self.sres = None
        if sres is not None:
            self.sres = numpy.tile(sres, (self.shape[0],1)) if sres.ndim == 1 else sres
            if self.sres.shape != self.shape:
                raise ValueError('Spectral resolution array has incorrect shape.')

        self.xpos = None 
        if xpos is not None:
            self.xpos = xpos
            if self.xpos.shape != self.shape:
                raise ValueError('X-coordinate array has incorrect shape.')

        self.ypos = None 
        if ypos is not None:
            self.ypos = ypos
            if self.ypos.shape != self.shape:
                raise ValueError('Y-coordinate array has incorrect shape.')

        self.area = numpy.ones(self.nspec, dtype=float) if area is None else numpy.atleast_1d(area)
        if self.area.size == 1:
            self.area = numpy.full(self.nspec, area, dtype=float)
        if self.area.size != self.nspec:
            raise ValueError('Could not construct area for each aperture; provide a single value '
                             'or one value per spectrum.')

        # Allocate attributes for primary and flux array fits headers
        self.prihdr = fits.Header() if prihdr is None else prihdr
        self.fluxhdr = self.prihdr.deepcopy() if fluxhdr is None else fluxhdr

        # Datacube rectification parameters
        self.pixelscale = None
        self.recenter = None
        self.width_buffer = None
        self.xs = None
        self.nx = None
        self.ys = None
        self.ny = None
        self.rect_rlim = None
        self.rect_sigma = None
        self.rect_channel = None
        self.rect_T = None

    @property
    def nspec(self):
        """Number of spectra in the datacube."""
        return self.shape[0]

    # TODO: Consolidate these functions between RowStackedSpectra and DataCube?
    def copy_to_array(self, attr='flux', waverange=None, nbins=None, select_bins=None,
                      missing_bins=None, unique_bins=None):
        r"""
        Return a copy of the selected data array with a limited
        spectral range and/or a selected set of spectra.

        The array size is always :math:`N_{\rm spec} \times N_{\rm
        wavelength}`. See :func:`copy_to_masked_array` for argument
        descriptions.

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
        Return a copy of the selected data array with selected
        pixels, wavelength ranges, and/or full spectra masked.

        This is functionally identical to :func:`copy_to_array`,
        except the output format is a `numpy.ma.MaskedArray`_. The
        pixels that are considered to be masked can be specified
        using the `flag` option. The array size is always
        :math:`N_{\rm spec} \times N_{\rm wavelength}`.
        
        Args:
            attr (:obj:`str`, optional):
                The attribute for the returned array. Can be 'flux',
                'ivar', 'sres', 'mask', 'xpos', or 'ypos'. For
                'mask', you're likely better off using
                :func:`copy_to_array`. Strings are always set to be
                lower-case, so capitalization shouldn't matter.
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
                A boolean array selecting the spectra to return.
            missing_bins (:obj:`list`, optional):
                A list of bin numbers that are missing from the
                output selection and that will be replaced with
                masked spectra. If specified, must also provide
                ``unique_bins``.
            unique_bins (array-like, optional):
                The indices of the bins that have valid spectra. The
                length of this array should match the select_bins
                array.

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
            mask |= self.mask if self.bitmask is None \
                                else self.bitmask.flagged(self.mask, flag=flag)

        # Create the output MaskedArray
        a = numpy.ma.MaskedArray(getattr(self, attr.lower()), mask=mask)

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

    def mean_sky_coordinates(self, waverange=None, response_func=None, per_pixel=True, offset=None,
                             flag=None, fluxwgt=False):
        r"""
        Compute the weighted or unweighted mean sky coordinates for
        each spectrum.

        Method cannot be performed if coordinate arrays are not
        available.

        .. warning::

            Flux-weighting the coordinates can produce spurious
            results in low-flux regimes.

        Args:
            waverange (array-like, optional):
                Two-element array with the first and last wavelength
                to include in the computation. Default is to use the
                full wavelength range.
            response_func (array-like, optional):
                A two-column array with the wavelength and
                transmission of a broad-band response function to use
                for the calculation. If None, no wavelength-dependent
                weighting function is used.
            per_pixel (:obj:`bool`, optional):
                When providing a response function, continue to
                calculate the statistics per pixel, as opposed to per
                angstrom.
            offset (:obj:`tuple`, optional):
                A two-tuple with an x,y offset to apply.
            flag (:obj:`str`, :obj:`list`, optional):
                One or more flag names that are considered when
                deciding if a pixel should be masked. The names
                *must* be a valid bit name as defined by
                :attr:`bitmask`.
            fluxwgt (:obj:`bool`, optional):
                Flag to weight by the flux when determining the mean
                coordinates.

        Returns:
            :obj:`tuple`: Two `numpy.ndarray`_ objects with the
            fiducial x and y on-sky positions of each spectrum in
            arcseconds relative to a given center. The shape of the
            two arrays is :math:`(N_{\rm spec},)`.
        """
        _offset = (0,0) if offset is None else offset
        if len(_offset) != 2:
            raise ValueError('Offset should be two numbers for the offset in x and y.')

        xpos = self.copy_to_masked_array(attr='xpos', waverange=waverange, flag=flag) + _offset[0]
        ypos = self.copy_to_masked_array(attr='ypos', waverange=waverange, flag=flag) + _offset[1]
        if fluxwgt:
            flux = self.copy_to_masked_array(waverange=waverange, flag=flag)

        # Set the response function
        dw = numpy.ones(self.nwave, dtype=float) if per_pixel \
                else angstroms_per_pixel(self.wave, log=self.log)
        _response_func = self.interpolate_to_match(response_func)

        # Get the normalization and return the flux- or un-weighted coordinates
        if fluxwgt:
            norm = numpy.ma.sum(flux*_response_func[None,:]*dw[None,:], axis=1)
            return numpy.ma.sum(flux*xpos*_response_func[None,:]*dw[None,:],axis=1)/norm, \
                    numpy.ma.sum(flux*ypos*_response_func[None,:]*dw[None,:],axis=1)/norm

        norm = numpy.sum(numpy.invert(numpy.ma.getmaskarray(xpos))*_response_func[None,:]
                                *dw[None,:], axis=1)
        return numpy.ma.sum(xpos*_response_func[None,:]*dw[None,:],axis=1)/norm, \
                    numpy.ma.sum(ypos*_response_func[None,:]*dw[None,:],axis=1)/norm

    def binned_on_sky_area(self, bin_indx, x=None, y=None, **kwargs):
        r"""
        Compute the on-sky area of a set of binned spectra.  

        The method attempts to calculate the overlapping area
        assuming the spectra were obtained with a set of circular
        fibers with the provided area (see :attr:`area`). If not
        provided in the method call (see ``x``, ``y``), the fiber
        fiducial coordinates used for the calculation are calculated
        using :func:`mean_sky_coordinates`.

        The overlapping area is calculated using the `shapely`_
        python package. If that package is not available, the
        returned area of each bin is simply the sum of the area of
        each fiber (does not account for any overlaps).

        Args:
            bin_indx (array-like):
                An array with size :math:`N_{\rm spec}` that gives
                which spaxels were included in each bin. Valid bins
                have indices of :math:`\geq 0`.
            x (array-like, optional):
                On-sky :math:`x` coordinate. Default is to calculate
                :math:`x` and :math:`y` using
                :func:`mean_sky_coordinates`.
            y (array-like, optional):
                On-sky :math:`y` coordinate. Default is to calculate
                :math:`x` and :math:`y` using
                :func:`mean_sky_coordinates`.
            **kwargs:
                Arguments passed directly to
                :func:`mean_sky_coordinates` for the determination of
                the sky coordinates, if they aren't provided directly.

        Returns:
            :obj:`tuple`: Two `numpy.ndarray`_ objects are returned.
            The first has the unique (non-negative) bin indices, and
            the second provides the on-sky area of that bin.
        """
        unique_bins, bin_count = numpy.unique(bin_indx, return_counts=True)
        indx = unique_bins > -1
        nbin = bin_count[indx]
        good_bins = unique_bins[indx]

        try:
            # Test that shapely was imported
            cascaded_union
        except:
            # It wasn't, so do the stupid calculation
            warnings.warn('Could not use \'shapely\' package to compute overlapping fiber area.' \
                          'Return the total fiber area.', ImportWarning)
            return good_bins, numpy.array([numpy.sum(self.area[bin_indx == b]) for b in good_bins])

        if x is None or y is None:
            x, y = self.mean_sky_coordinates(**kwargs)
        if x.size != self.nspec or y.size != self.nspec:
            raise ValueError('Must provide an x and y coordinate for each spectrum.')

        # Assume apertures are circular
        radius = numpy.sqrt(self.area/numpy.pi)
        bin_area = numpy.empty(nbin.size, dtype=float)
        for i, b in enumerate(good_bins):
            _x = x[bin_indx == b]
            _y = y[bin_indx == b]
            bin_area[i] = cascaded_union([Point(xx,yy).buffer(r,64) \
                                                    for xx, yy, r in zip(_x,_y,radius)]).area
        return good_bins, bin_area

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
            return numpy.mean(cen_wave)

        norm = numpy.ma.sum(numpy.invert(numpy.ma.getmaskarray(flux))
                            * _response_func[None,:] * dw[None,:], axis=1)
        cen_wave = numpy.ma.sum(numpy.invert(numpy.ma.getmaskarray(flux)) * self.wave[None,:]
                                * _response_func[None,:] * dw[None,:], axis=1) / norm
        return numpy.mean(cen_wave)

    # TODO: This is virtually identical to the function in DataCube.
    # Consolidate?
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
            mean S/N.
    
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

        # Get the moments
        response_integral = numpy.sum(numpy.invert(numpy.ma.getmaskarray(flux))
                                        * (_response_func*dw)[None,:], axis=1)
        signal = numpy.ma.divide(numpy.ma.sum(flux*(_response_func*dw)[None,:], axis=1),
                                 response_integral)
        variance = numpy.ma.divide(numpy.ma.sum(numpy.ma.power(ivar, -1.) \
                                    * (_response_func*dw)[None,:], axis=1), response_integral)
        snr = numpy.ma.divide(numpy.ma.sum(snr*(_response_func*dw)[None,:], axis=1),
                              response_integral)
        # Done
        return signal, variance, snr

    # Most of the rest of this deals with covariance!
    def _cube_dimensions(self, pixelscale=None, recenter=None, width_buffer=None, redo=False):
        """
        Determine the on-sky dimensions of the reconstructed image for
        all wavelength channels.

        When first calling this method, the pixel scale, recentering
        approach, and buffer must be provided, and they are saved to
        :attr:`pixelscale`, :attr:`recenter`, :attr:`width_buffer`.
        The spatial dimensions of the cube are calculated using the
        data in :attr:`xpos` and :attr:`ypos`, and the results are
        saved to :attr:`xs`, :attr:`ys`, :attr:`nx`, and :attr:`ny`.
        This calculation only needs to be done once per instance and
        settings for ``pixelscale``, ``recenter``, and
        ``width_buffer``. The calculation is forced to be done again
        if any of these are provided to the method and different from
        the saved attributes or if ``redo`` is True.

        Args:
            pixelscale (:obj:`float`, optional):
                Desired pixel scale in arcsec.  
            recenter (:obj:`bool`, optional):
                Flag to recenter the coordinate system.
            width_buffer (:obj:`int`, optional):
                Number of pixels to use as buffer for the image
                reconstruction.
            redo (:obj:`bool`, optional):
                Force the recalculation of the cube dimensions if
                they are already defined.
        """
        if self.xpos is None or self.ypos is None:
            raise ValueError('Cannot construct cube dimensions without the wavelength dependent '
                             'aperture positions.')
        # Check the input
        if self.pixelscale is None and pixelscale is None:
            raise ValueError('Must provide pixelscale on first use of _cube_dimensions.')
        if self.recenter is None and recenter is None:
            raise ValueError('Must provide recenter on first use of _cube_dimensions.')
        if self.width_buffer is None and width_buffer is None:
            raise ValueError('Must provide width_buffer on first use of _cube_dimensions.')

        # Check if the calculation needs to redone. This will always be
        # true on the first use of the method.
        _redo = redo or (pixelscale is not None and self.pixelscale != pixelscale) \
                        or (recenter is not None and self.recenter != recenter) \
                        or (width_buffer is not None and self.width_buffer != width_buffer) \
                        or self.xs is None or self.nx is None or self.ys is None or self.ny is None
        if not _redo:
            return

        # Save the parameters used to create the dimensions. Only save
        # the ones that are actually provided.
        if pixelscale is not None:
            self.pixelscale = pixelscale
        if recenter is not None:
            self.recenter = recenter
        if width_buffer is not None:
            self.width_buffer = width_buffer

        # Get the size in each dimension
        # TODO: This negative here is just to ensure that the
        # calculation matches what's done by the DRP, which defines
        # xpos as increasing from E to W.
        minx = numpy.amin(-self.xpos)
        maxx = numpy.amax(-self.xpos)
        Dx = numpy.floor(maxx-minx)

        miny = numpy.amin(self.ypos)
        maxy = numpy.amax(self.ypos)
        Dy = numpy.floor(maxy-miny)

        # Force the size to be even and the same in both dimensions
        Dx = Dx if Dx > Dy else Dy
        self.nx = int(numpy.floor(Dx/self.pixelscale)+self.width_buffer)
        if self.nx % 2 != 0:
            self.nx += 1
        self.ny = self.nx

        # Set the starting coordinate
        self.xs = -self.nx*self.pixelscale/2.
        self.ys = -self.ny*self.pixelscale/2.

        # Offset to the center, if requested
        if self.recenter:
            self.xs = self.xs + (minx+maxx)/2.0
            self.ys = self.ys + (miny+maxy)/2.0

    def _reinit_rectification_parameters(self, channel=None, pixelscale=None, rlim=None,
                                         sigma=None, recenter=None, width_buffer=None):
        """
        Determine if the rectification paraemters need to be reinitialized.
        """
        if None in [self.pixelscale, self.recenter, self.width_buffer, self.xs, self.nx, self.ys,
                    self.ny, self.rect_rlim, self.rect_sigma, self.rect_channel]:
            return True
        if channel is not None and self.rect_channel != channel:
            return True
        if pixelscale is not None and self.pixelscale != pixelscale:
            return True
        if rlim is not None and self.rect_rlim != rlim:
            return True
        if sigma is not None and self.rect_sigma != sigma:
            return True
        if recenter is not None and self.recenter != recenter:
            return True
        if width_buffer is not None and self.width_buffer != width_buffer:
            return True
        return False

    def _init_rectification_parameters(self, channel=None, pixelscale=None, rlim=None, sigma=None,
                                       recenter=None, width_buffer=None):
        """
        Perform common preparation of rectification parameters.

        Args:
            channel (:obj:`int`, optional):
                Index of the spectral channel for which to calculate
                the transfer matrix.
            pixelscale (:obj:`float`, optional):
                Desired pixel (spaxel) scale in arcsec.
            rlim (:obj:`float`, optional):
                The limiting radius of the image reconstruction
                (Gaussian) kernel in arcseconds.
            sigma (:obj:`float`, optional):
                The sigma of the image reconstruction (Gaussian)
                kernel in arcseconds.
            recenter (:obj:`bool`, optional):
                Flag to recenter the coordinate system.
            width_buffer (:obj:`int`, optional):
                Number of pixels to use as buffer for the image
                reconstruction
        """
        if not self._reinit_rectification_parameters(channel=channel, pixelscale=pixelscale,
                                                     rlim=rlim, sigma=sigma, recenter=recenter,
                                                     width_buffer=width_buffer):
            # Nothing needs to be reinitialized
            return

        # Get the cube dimensions; may not necessarily match DRP calculation
        self._cube_dimensions(pixelscale=pixelscale, recenter=recenter, width_buffer=width_buffer)

        # Set the kernel parameters if they've changed
        if rlim is not None:
            self.rect_rlim = rlim
        if sigma is not None:
            self.rect_sigma = sigma

        # Set the channel if it's changed
        if channel is not None:
            self.rect_channel = channel

    # TODO: Allow for spectrum-dependent (rlim,sigma), or scale a
    # single (rlim,sigma) by the sqrt of the area to account for
    # different size apertures?
    def rectification_transfer_matrix(self, channel, pixelscale, rlim, sigma, recenter=False,
                                      width_buffer=10, quiet=False, rej_flag=None):
        r"""
        Calculate the transfer matrix used to produce a reconstructed
        image of the fiber data at the specified wavelength channel
        using Shepard's method.

        This is done directly using the available on-sky x and y
        coordinates of each fiber as a function of wavelength taken
        from :attr:`xpos` and :attr:`ypos` attributes.

        .. todo::

            - Give more detail on the pixels at which the radius is
              calculated.
            - Describe the algorithm in more depth.

        Args:
            channel (:obj:`int`):
                Index of the spectral channel for which to calculate
                the transfer matrix.
            pixelscale (:obj:`float`):
                Desired pixel (spaxel) scale in arcsec.
            rlim (:obj:`float`):
                The limiting radius of the image reconstruction
                (Gaussian) kernel in arcseconds.
            sigma (:obj:`float`):
                The sigma of the image reconstruction (Gaussian)
                kernel in arcseconds.
            recenter (:obj:`bool`, optional):
                Flag to recenter the coordinate system.
            width_buffer (:obj:`int`, optional):
                Number of pixels to use as buffer for the image
                reconstruction
            quiet (:obj:`bool`, optional):
                Suppress terminal output
            rej_flag (:obj:`str`, optional):
                Bit name to ignore in rectification. Ignored if
                :attr:`bitmask` is None. If None, bitmasks are
                ignored. Use ``'any'`` to mask pixels with any
                flagged bits.

        Returns:
            `scipy.sparse.csr_matrix`_: Transfer matrix
            :math:`{\mathbf T}`

        """
        # Is the matrix already available?
        if self.rect_T is not None \
                and not self._reinit_rectification_parameters(channel=channel,
                                                              pixelscale=pixelscale, rlim=rlim,
                                                              sigma=sigma, recenter=recenter,
                                                              width_buffer=width_buffer):
            return self.rect_T

        # Create it!  
        self._init_rectification_parameters(channel=channel, pixelscale=pixelscale, rlim=rlim,
                                            sigma=sigma, recenter=recenter,
                                            width_buffer=width_buffer)

        # Dimensions of the sparse matrix are
        nim = self.nx*self.ny                   # The number of image pixels
        # by the number of fiber spectra (self.nspec)

        # Get the list of non-zero pixel values in the transfer matrix
        i = numpy.arange(self.nx)
        j = numpy.arange(self.ny)
        ii,jj = numpy.meshgrid(i, j, indexing='ij')         # Mesh of i,j pixel indices

        sp = numpy.empty((self.nx,self.ny), dtype=float)    # Holds spectrum index
        ij = (ii*self.ny+jj)                                # Holds image pixel index
        r2 = numpy.empty((self.nx,self.ny), dtype=float)    # Holds radii
        tot = numpy.zeros((self.nx,self.ny), dtype=float)   # Holds the sum of the weights

        s2 = numpy.square(self.rect_sigma/self.pixelscale)  # sigma^2 of Gaussian
        rl2 = numpy.square(self.rect_rlim/self.pixelscale)  # radius^2 of Gaussian limit

        non_zero_spc = numpy.empty(0, dtype=int)            # Holds triplet spectrum index
        non_zero_pix = numpy.empty(0, dtype=int)            # Holds triplet image index
        non_zero_wgt = numpy.empty(0, dtype=float)          # Holds triplet weight

        # Do not include any pixels with zero inverse variance or
        # pixels that have been masked (either by the boolean mask
        # attribute or flagged with the provided mask bits)
        mask = numpy.invert(self.ivar[:,self.rect_channel] > 0.0)
        if rej_flag is not None and self.bitmask is not None:
            _rej_flag = rej_flag if isinstance(rej_flag, list) or rej_flag != 'any' else None
            mask |= self.bitmask.flagged(self.mask[:,self.rect_channel], flag=_rej_flag)
        elif self.mask is not None:
            mask |= self.mask[:,self.rect_channel]

        # TODO: Can optimize this further?
        for k in range(self.nspec):
            if mask[k]:
                continue

            # Fill spectrum index
            sp.fill(k)

            # NOTE: Calculating full matrix is actually faster than
            # determining submatrix for calculation

            # NOTE: The negative sign in xpos is again to match the
            # calculation done by the DRP, which defines xpos as
            # increasing from E to W.

            # Calculate the distance.
            # ---- WITH RESPECT TO THE EDGE OF THE FIRST PIXEL ----
            #  - matches DRP, but why?
            r2 = numpy.square((-self.xpos[k,self.rect_channel]-self.xs)/self.pixelscale - ii) \
                 + numpy.square((self.ypos[k,self.rect_channel]-self.ys)/self.pixelscale - jj)
            # ---- WITH RESPECT TO THE CENTER OF THE FIRST PIXEL ----
#           r2 = numpy.square( (self.hdu['XPOS'].data[k,channel]-self.xs)/pixelscale-0.5 - ii) \
#                + numpy.square((self.hdu['YPOS'].data[k,channel]-self.ys)/pixelscale-0.5 - jj)

            # Append new indices and weights within rlim
            # TODO: Can this be done quicker if I'm not appending things?
            non_zero_spc = numpy.append(non_zero_spc, sp[r2 < rl2])
            non_zero_pix = numpy.append(non_zero_pix, ij[r2 < rl2])
            wgt = numpy.exp(-r2[r2 < rl2]/s2/2.0)
            tot[r2<rl2] += self.area[k] * wgt
            non_zero_wgt = numpy.append(non_zero_wgt, wgt)

            if not quiet:
                print('Transfer Matrix {:2.1%}'.format((k+1)/self.nspec), end='\r')
        if not quiet:
            print('Transfer Matrix {:2.1%}'.format(1.))

        # Normalize the result and scale by the pixel size to ensure the
        # output cube is in units of calibrated flux density per pixel
        non_zero_wgt *= numpy.square(self.pixelscale) \
                            / tot[numpy.unravel_index(non_zero_pix.astype(int), (self.nx,self.ny))]

        # Set the transfer matrix to a sparse object
        self.rect_T = sparse.coo_matrix((non_zero_wgt, (non_zero_pix, non_zero_spc)),
                                        shape=(nim,self.nspec)).tocsr()
        # TODO: Does this need to be returned?
        return self.rect_T

    # TODO: Include mask propagation

    def rectify_wavelength_plane(self, channel, pixelscale=None, rlim=None, sigma=None,
                                 recenter=False, width_buffer=10, quiet=False, rej_flag=None,
                                 return_ivar=False, return_covar=False):
        r"""
        Return a rectified image for the specified wavelength channel
        using Shepard's method.

        First, the rectification transfer matrix, :math:`{\mathbf
        T}`, is calculated using
        :func:`rectification_transfer_matrix`. Then the wavelength
        image, :math:`{\mathbf I}` is calculated:

        .. math::
        
            {\mathbf T} \times {\mathbf F} = {\mathbf I}

        where :math:`{\mathbf F}` is the vector of fluxes in the
        selected wavelength channel for all the aperture
        measurements. The calculation of :math:`{\mathbf I}` produces
        a vector, but this is reshaped into a 2D array of shape
        :math:`(N_x, N_y)` on output.

        Args:
            channel (:obj:`int`):
                Index of the spectral channel for which to calculate
                the transfer matrix.
            pixelscale (:obj:`float`, optional):
                Desired pixel (spaxel) scale in arcsec. If None,
                :attr:`pixelscale` will be used, but fault if that
                value is also None.
            rlim (:obj:`float`, optional):
                The limiting radius of the image reconstruction
                (Gaussian) kernel in arcseconds. If None,
                :attr:`rect_rlim` will be used, but fault if that
                value is also None.
            sigma (:obj:`float`, optional):
                The sigma of the image reconstruction (Gaussian)
                kernel in arcseconds. If None, :attr:`rect_sigma`
                will be used, but fault if that value is also None.
            recenter (:obj:`bool`, optional):
                Flag to recenter the coordinate system.
            width_buffer (:obj:`int`, optional):
                Number of pixels to use as buffer for the image
                reconstruction
            quiet (:obj:`bool`, optional):
                Suppress terminal output
            rej_flag (:obj:`str`, optional):
                Bit name to ignore in rectification. Ignored if
                :attr:`bitmask` is None. If None, bitmasks are
                ignored. Use ``'any'`` to mask pixels with any
                flagged bits.
            return_ivar (:obj:`bool`, optional):
                Return the nominal inverse variance image (i.e. the
                diagonal of the covariance matrix).
            return_covar (:obj:`bool`, optional):
                Return the covariance matrix. Note that the
                covariance matrix is calculated if either
                ``return_ivar`` or ``return_covar`` are True. If
                ``return_covar`` is True, ``return_ivar`` is ignored
                (because it's automatically true).

        Returns:
            `numpy.ndarray`_, :obj:`tuple`: If both ``return_ivar``
            and ``return_covar`` are False, only the reconstructed
            image is returned. If ``return_covar`` is True, a
            :class:`mangadap.util.covariance.Covariance` object is
            also returned. If ``return_covar`` is False and
            ``return_ivar`` is True, the second returned object is a
            `numpy.ndarray`_ with the propagated inverse variance
            image.
        """
        # Set the transfer matrix (set to self.rect_T; don't need to
        # keep the returned matrix)
        self.rectification_transfer_matrix(channel, pixelscale, rlim, sigma, recenter=recenter,
                                           width_buffer=width_buffer, quiet=quiet,
                                           rej_flag=rej_flag)

        flux = self.rect_T.dot(self.flux[:,channel]).reshape(self.nx,self.ny)

        # Return the regridded data with the proper shape (nx by ny)
        if not return_covar and not return_ivar:
            return flux

        covar = self._formal_covariance_matrix(quiet=quiet)
        if return_covar:
            # Return the rectified image and the full covariance matrix
            return flux, covar

        # Only return the inverse variance, not the full covariance matrix
        return flux, numpy.ma.power(numpy.diagonal(covar.toarray()),
                                    -1).filled(0.0).reshape(self.nx,self.ny)

    def _formal_covariance_matrix(self, csr=False, quiet=False):
        r"""
        Return the formal covariance matrix.
        
        The formal covariance matrix is:

        .. math::
        
             {\mathbf C} = {\mathbf T} \times {\mathbf \Sigma} \times
             {\mathbf T}^{\rm T},

        where :math:`{\mathbf \Sigma}` is the covariance matrix for the
        row-stacked spectra for the specified wavelength channel.

        This method requires that :attr:`rect_T` not be None, it
        assumes that the current transfer matrix is correct for this
        channel, and it is largely a wrapper for
        :func:`mangadap.util.covariance.Covariance.from_matrix_multiplication`.

        Args:
            csr (:obj:`bool`, optional):
                Instead of returning a
                :class:`mangadap.util.covariance.Covariance` object,
                return the covariance matrix as a
                `scipy.sparse.csr_matrix`_ object. Primarily used by
                :func:`covariance_cube` for collating the covariance
                matrix of each wavelength channel before combining
                them into a single
                :class:`mangadap.util.covariance.Covariance` object.
            quiet (:obj:`bool`, optional):
                Suppress terminal output

        Returns:
            :class:`mangadap.util.covariance.Covariance`,
            `scipy.sparse.csr_matrix`_: The covariance matrix for the
            designated wavelength channel. The return type depends on
            ``csr``.
        """
        if self.rect_T is None:
            raise ValueError('To calculate the formal covariance matrix, must first calculate '
                             'the rectification transfer matrix.')
        var = numpy.ma.power(self.ivar[:,self.rect_channel], -1.0).filled(0.0)
        covar = Covariance.from_matrix_multiplication(self.rect_T, var)
        return covar.full() if csr else covar

    def covariance_matrix(self, channel, pixelscale=None, rlim=None, sigma=None, recenter=False,
                          width_buffer=10, csr=False, quiet=False, rej_flag=None):
        r"""
        Return the covariance matrix for the specified wavelength
        channel.

        For a regrided cube image with :math:`N_x\times N_y` pixels, the
        covariance matrix has a size that is :math:`(N_x N_y)\times (N_x
        N_y)`; however, the majority of these pixels will be zero.
        Therefore, the covariance matrix is stored as a sparse matrix
        and interfaced with using the
        :class:`mangadap.util.covariance.Covariance` object class.

        The value of the covariance matrix at pixel :math:`(i,j)` is the
        covariance between pixels :math:`(n_{x,0},n_{y,0})` and
        :math:`(n_{x,1},n_{y,1})` at the specified wavelength channel of
        the reconstructed datacube image, where

        .. math::

            n_{x,0} &= \lfloor i / N_y \rfloor \\
            n_{y,0} &= i - n_{x,0} N_y \\
            n_{x,1} &= \lfloor j / N_y \rfloor \\
            n_{y,1} &= j - n_{x,1} N_y

        and :math:`\lfloor m\rfloor` is the "floor" of :math:`m`. The
        diagonal of the covariance matrix (:math:`i=j`) provides the
        variance.

        .. warning::

            Because the covariance is calculated for a higher
            dimensional array, it's important to note how the pixels
            in the image map to the relevant covariance array
            locations. See the equations above. Also see
            :func:`mangadap.util.covariance.Covariance.transpose_raw_shape`.

        The covariance matrix provided by this method is always the
        formal covariance matrix defined as

        .. math::
        
             {\mathbf C} = {\mathbf T} \times {\mathbf \Sigma} \times
             {\mathbf T}^{\rm T},

        where :math:`{\mathbf \Sigma}` is the covariance matrix for
        the row-stacked spectra for the specified wavelength channel,
        which is taken to be a diagonal matrix with the flux
        variances along the diagonal.

        Args:
            channel (:obj:`int`):
                Index of the spectral channel for which to calculate
                the covariance matrix.
            pixelscale (:obj:`float`, optional):
                Desired pixel (spaxel) scale in arcsec. If None,
                :attr:`pixelscale` will be used, but fault if that
                value is also None.
            rlim (:obj:`float`, optional):
                The limiting radius of the image reconstruction
                (Gaussian) kernel in arcseconds. If None,
                :attr:`rect_rlim` will be used, but fault if that
                value is also None.
            sigma (:obj:`float`, optional):
                The sigma of the image reconstruction (Gaussian)
                kernel in arcseconds. If None, :attr:`rect_sigma`
                will be used, but fault if that value is also None.
            recenter (:obj:`bool`, optional):
                Flag to recenter the coordinate system.
            width_buffer (:obj:`int`, optional):
                Number of pixels to use as buffer for the image
                reconstruction
            csr (:obj:`bool`, optional):
                Instead of returning a
                :class:`mangadap.util.covariance.Covariance` object,
                return the covariance matrix as a
                `scipy.sparse.csr_matrix`_ object. Primarily used by
                :func:`covariance_cube` for collating the covariance
                matrix of each wavelength channel before combining
                them into a single
                :class:`mangadap.util.covariance.Covariance` object.
            quiet (:obj:`bool`, optional):
                Suppress terminal output
            rej_flag (:obj:`str`, optional):
                Bit name to ignore in rectification. Ignored if
                :attr:`bitmask` is None. If None, bitmasks are
                ignored. Use ``'any'`` to mask pixels with any
                flagged bits.

        Returns:
            :class:`mangadap.util.covariance.Covariance`,
            `scipy.sparse.csr_matrix`_: The covariance matrix for the
            designated wavelength channel. The return type depends on
            `csr`.
        """
        # Set the transfer matrix 
        self.rectification_transfer_matrix(channel, pixelscale, rlim, sigma, recenter=recenter,
                                           width_buffer=width_buffer, quiet=quiet,
                                           rej_flag=rej_flag)
        # Return the formal covariance matrix
        return self._formal_covariance_matrix(csr=csr, quiet=quiet)

    def covariance_cube(self, channels=None, pixelscale=None, rlim=None, sigma=None,
                        recenter=False, width_buffer=10, csr=False, quiet=False, rej_flag=None):
        """
        Return the covariance matrices for many wavelength channels.

        This is primarily a wrapper for :func:`covariance_matrix`
        that iterates through one or more wavelength channels and
        constructs a 3D :class:`mangadap.util.covariance.Covariance`
        object.
        
        Args:
            channels (:obj:`int`, array-like, optional):
                Indices of the spectral channels for which to
                calculate the covariance matrix. If None, the
                covariance matrix is calculated for *all* channels.
            pixelscale (:obj:`float`, optional):
                Desired pixel (spaxel) scale in arcsec. If None,
                :attr:`pixelscale` will be used, but fault if that
                value is also None.
            rlim (:obj:`float`, optional):
                The limiting radius of the image reconstruction
                (Gaussian) kernel in arcseconds. If None,
                :attr:`rect_rlim` will be used, but fault if that
                value is also None.
            sigma (:obj:`float`, optional):
                The sigma of the image reconstruction (Gaussian)
                kernel in arcseconds. If None, :attr:`rect_sigma`
                will be used, but fault if that value is also None.
            recenter (:obj:`bool`, optional):
                Flag to recenter the coordinate system.
            width_buffer (:obj:`int`, optional):
                Number of pixels to use as buffer for the image
                reconstruction
            csr (:obj:`bool`, optional):
                Instead of returning a
                :class:`mangadap.util.covariance.Covariance` object,
                return the covariance matrix as a
                `scipy.sparse.csr_matrix`_ object. Primarily used by
                :func:`covariance_cube` for collating the covariance
                matrix of each wavelength channel before combining
                them into a single
                :class:`mangadap.util.covariance.Covariance` object.
            quiet (:obj:`bool`, optional):
                Suppress terminal output
            rej_flag (:obj:`str`, optional):
                Bit name to ignore in rectification. Ignored if
                :attr:`bitmask` is None. If None, bitmasks are
                ignored. Use ``'any'`` to mask pixels with any
                flagged bits.

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
                print('Covariance Cube {0}/{1}'.format(i+1,nc), end='\r')
            CovCube[i] = self.covariance_matrix(_channels[i], pixelscale=pixelscale, rlim=rlim,
                                                sigma=sigma, recenter=recenter,
                                                width_buffer=width_buffer, csr=True, quiet=True,
                                                rej_flag=rej_flag)
        if not quiet:
            print('Covariance Cube {0}/{0}'.format(nc))
        return Covariance(CovCube, input_indx=_channels) if not csr else CovCube

    def instrumental_dispersion_plane(self, channel, dispersion_factor=None, pixelscale=None,
                                      rlim=None, sigma=None, recenter=None, width_buffer=None,
                                      quiet=False, rej_flag=None):
        r"""
        Return the instrumental dispersion for the rectified
        wavelength plane.

        Method requires that :attr:`sres` is not None.

        The method first calculates the rectification transfer
        matrix, :math:`{\mathbf T}`, using
        :func:`regrid_transfer_matrix`. The transfer matrix is then
        used to construct the rectified wavelength plane image,
        :math:`{\mathbf I}`, by computing :math:`{\mathbf T} \times
        {\mathbf F} = {\mathbf I}`, where :math:`{\mathbf F}` is the
        vector of the fiber fluxes. Under the assumption that the
        line-spread-function (LSF) is Gaussian, we determine the
        instrumental dispersion for the data in the wavelength
        channel of the reconstructed CUBE by calculating the second
        moment of the weighted sum of Gaussians of the appropriate
        dispersion. Assuming all the Gaussians have the normalized
        form:

        .. math::
            
            \mathcal{N}(x|\mu=0,\sigma) = \frac{1}{\sqrt{2\pi}\sigma}
            \exp\left(-\frac{x^2}{2\sigma^2}\right),

        the combined instrumental dispersion becomes

        .. math::
            
            \sigma_{\rm inst,j}^2 = \frac{\sum_i T_{ij}
            \sigma^2_i}{\sum_i T_{ij}},

        where :math:`T_{ij}` are the elements of the transfer matrix.
        In terms of matrix multiplications, this can be written as

        .. math::

            {\mathbf S} = \frac{ {\mathbf T} \times {\mathbf V} }{
            {\mathbf T_c} },

        where :math:`{\mathbf T_c} = {\mathbf T_c} \times {\mathbf 1}`
        is the vector with the sum :math:`\sum_j T_{ij}`,
        :math:`{\mathbf V}` is the instrumental variance for all fibers
        at the designated wavelength plane, and :math:`{\mathbf S}` is
        the variance for all the spaxels in the reconstructed wavelength
        image; the division by :math:`{\mathbf T_c}` is element-wise.

        The returned matrix is the element-wise square-root of
        :math:`{\mathbf S}`, rearranged into a 2D array of size
        :attr:`nx` by :attr:`ny`.

        Args:

            channel (:obj:`int`):
                Index of the spectral channel for which to calculate
                the covariance matrix.
            dispersion_factor (:obj:`float`, optional):
                Artificially change the resolution by multiplying the
                instrumental LSF by this factor before calculating
                the reconstructed dispersion. If None, factor is
                unity.
            pixelscale (:obj:`float`, optional):
                Desired pixel (spaxel) scale in arcsec. If None,
                :attr:`pixelscale` will be used, but fault if that
                value is also None.
            rlim (:obj:`float`, optional):
                The limiting radius of the image reconstruction
                (Gaussian) kernel in arcseconds. If None,
                :attr:`rect_rlim` will be used, but fault if that
                value is also None.
            sigma (:obj:`float`, optional):
                The sigma of the image reconstruction (Gaussian)
                kernel in arcseconds. If None, :attr:`rect_sigma`
                will be used, but fault if that value is also None.
            recenter (:obj:`bool`, optional):
                Flag to recenter the coordinate system.
            width_buffer (:obj:`int`, optional):
                Number of pixels to use as buffer for the image
                reconstruction
            quiet (:obj:`bool`, optional):
                Suppress terminal output
            rej_flag (:obj:`str`, optional):
                Bit name to ignore in rectification. Ignored if
                :attr:`bitmask` is None. If None, bitmasks are
                ignored. Use ``'any'`` to mask pixels with any
                flagged bits.

        Returns:
            `numpy.ndarray`_: The reconstructed LSF dispersion for
            the specified wavelength channel.
        """
        if self.sres is None:
            raise ValueError('Spectral resolution must be available for this method.')

        # Get the instrumental dispersion
        disp = self.wave[channel]/self.sres[:,channel]/DAPConstants.sig2fwhm
        if dispersion_factor is not None:
            disp *= dispersion_factor

        # Set the transfer matrix 
        self.rectification_transfer_matrix(channel, pixelscale, rlim, sigma, recenter=recenter,
                                           width_buffer=width_buffer, quiet=quiet,
                                           rej_flag=rej_flag)

        # Rectify the data
        Tc = self.rect_T.sum(axis=1).flatten()
        Tc[numpy.invert(Tc>0)] = 1.0                # Control for zeros
        # NOTE: This returns a numpy.ndarray instead of a numpy.matrix
        return numpy.asarray(numpy.sqrt(self.rect_T.dot(numpy.square(disp)) 
                                / Tc).reshape(self.nx, self.ny))
