# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Provides a set of functions to handle and alter the instrumental
resolution.

----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import warnings

from IPython import embed

import numpy

from scipy import interpolate
from scipy.special import erf

import astropy.constants
from .constants import DAPConstants

from .sampling import spectrum_velocity_scale

#from matplotlib import pyplot

class VariableGaussianKernel:
    r"""
    Support class for variable sigma convolution.  See
    :func:`convolution_variable_sigma`.

    Modest modification of M. Cappellari's
    :func:`ppxf.ppxf_util.gaussian_filter1d` function.

    Args:
        sigma (numpy.ndarray):
            Coordinate-dependent standard deviation of the Gaussian
            kernel in pixels.
        minsig (:obj:`float`, optional):
            The minimum sigma value allowed.
        nsig (:obj:`float`, optional):
            The extent of the calculation of the kernel for the
            convolution.  By default, the kernel is not calculated
            beyond 3*sigma.
        integral (:obj:`bool`, optional):
            Construct the kernel based on the integral of the Gaussian
            over the pixel width using the error function, instead of
            just sampling the Gaussian function directly.

    Raises:
        ValueError:
            Raised if the input sigma has no length.

    Attributes:
        sigma (numpy.ndarray):
            Coordinate-dependent standard deviation of the Gaussian
            kernel.
        p (int):
            The maximum number of pixels on either side of the kernel
            center needed for the convolution.
        m (int):
            Number of pixels in the kernel (:math:`m = 2p + 1`).
        kernel (numpy.ndarray):
            The calculated kernel to use for the convolution.
    """
    def __init__(self, sigma, minsig=0.01, nsig=3.0, integral=False):
        self.n = sigma.size                                     # Vector length
        if self.n == 0:
            raise ValueError('Input sigma has zero length!')
        self.sigma = sigma.clip(min=minsig)                     # Force sigmas to minimum 
        self.p = int(numpy.ceil(numpy.amax(nsig*self.sigma)))   # Kernel covers up to nsig*sigma
        self.m = 2*self.p + 1                                   # Kernel length
        x = numpy.linspace(-self.p, self.p, self.m)             # X^2 for kernel
        x2 = numpy.square(x)

        # Kernel will have size m x n
        self.kernel = (erf((x[:,None]+0.5)/numpy.sqrt(2)/self.sigma) 
                            - erf((x[:,None]-0.5)/numpy.sqrt(2)/self.sigma))/2. if integral else \
                      numpy.exp(-x2[:, None]/(2*numpy.square(self.sigma))) \
                                / numpy.sqrt(2*numpy.pi) / self.sigma

        self.kernel /= numpy.sum(self.kernel, axis=0)[None, :]       # Normalize kernel


    def _check_shape(self, y, ye=None):
        """
        Make sure that the shapes are appropriate for the defined kernel.
        """
        if len(y.shape) != 1:
            raise ValueError('y must be a 1D array!')
        if y.size != self.n:
            raise ValueError('y and sigma must have the same shape!')
        if ye is not None and (len(ye.shape) != 1 or ye.size != self.n):
            raise ValueError('ye length does not must have the correct shape!')


    def _create_a(self, y):
        a = numpy.zeros(self.kernel.shape)
        for i in range(self.m):
            a[i,self.p:-self.p] = y[i:self.n-self.m+i+1]
        return a


    def convolve(self, y, ye=None):
        """
        Convolve the input vector with the variable-width Gaussian.

        The input `y` must be a 1D vector with the same length as
        :attr:`sigma`.  If errors are provided, they must also adhere to
        these shape limitations and are propagated to the convolved
        vector.

        Args:
            y (numpy.ndarray):
                Vector to convolve
            ye (numpy.ndarray, optional):
                Error in the vector to convolve
    
        Raises:
            ValueError:
                Raised if *y* is not a 1D vector, or if the shape of *y*
                and *sigma* (and *ye* if provided) are different.

        Returns:
            `numpy.ndarray`_: The convolved vector. If the errors are
            provided, the function returns two vectors, the convolved
            vector and its error.
        """
        self._check_shape(y, ye=ye)

        # Create m copies of the shifted input function
        a = self._create_a(y)
        if ye is None:
            return numpy.sum(a*self.kernel, axis=0)

        # Construct the error
        ae = self._create_a(ye)
        return numpy.sum(a*self.kernel, axis=0), \
                    numpy.sqrt(numpy.sum(numpy.square(ae*self.kernel), axis=0))


def convolution_variable_sigma(y, sigma, ye=None, integral=False, large_sigma=10.):
    r"""
    Convolve a discretely sampled function :math:`y(x)` with a Gaussian
    kernel, :math:`g`, where the standard deviation of the kernel is a
    function of :math:`x`, :math:`\sigma(x)`.  Nominal calculations can
    be performed to propagate the error in the result; these
    calculations **do not** include the covariance between the pixels,
    which will mean that the calculations likely have significant error!

    The convolution is defined as:

    .. math::

        (y\ast g)(x) &= \int_{-\infty}^{\infty} y(X)\ g(\sigma,x-X)\ dX \\
                     &= \int_{-\infty}^{\infty} \frac{y(X)}{\sqrt{2\pi}\
                        \sigma(X)}\ \exp\left(-\frac{(x-X)^2}{2\
                        \sigma(X)^2}\right) dX .

    To minimize edge effects and account for the censoring of the data
    (finite range in :math:`x`), the convolution is actually calculated
    as a definite integral and normalized as follows:

    .. math::

        (y\ast g)(x) \sim\frac{
        \int_{x-n_\sigma*\sigma(x)}^{x+n_\sigma*\sigma(x)} y(X)\
        g(\sigma,x-X)\ dX}{
        \int_{x-n_\sigma*\sigma(x)}^{x+n_\sigma*\sigma(x)}
        g(\sigma,x-X)\ dX} .

    The above is identical to getting the weighted mean of :math:`y` at
    each position :math:`x`, where the weights are defined by a Gaussian
    kernel centered at :math:`x` with a variable dispersion.

    Use of this function requires:
        - *y* and *sigma* must be 1D vectors
        - *y* and *sigma* must be uniformly sampled on the same grid
        - *sigma* must be in pixel units.

    .. todo::
        This should be deprecated and/or merged with
        :class:`VariableGaussianKernel`.

    Args:
        y (`numpy.ndarray`_):
            A uniformly sampled function to convolve.
        sigma (`numpy.ndarray`_):
            The standard deviation of the Gaussian kernel sampled at
            the same positions as ``y``. The units of ``sigma``
            *must* be in pixels.
        ye (`numpy.ndarray`_, optional):
            Errors in the function :math:`y(x)`. Must be the same
            shape as ``y``.
        integral (:obj:`bool`, optional):
            Construct the kernel based on the integral of the
            Gaussian over the pixel width using the error function,
            instead of just sampling the Gaussian function directly.
            See :class:`VariableGaussianKernel`.
        large_sigma (:obj:`float`, optional):
            A convenience parameter that causes a warning to be
            issued if any of the values in ``sigma`` are larger than
            this threshold.

    Returns:
        `numpy.ndarray`_: Arrays with the convolved function
        :math:`(y\ast g)(x)` sampled at the same positions as the
        input :math:`x` vector and its error. The second array will
        be returned as None if the error vector is not provided.
    """
    if numpy.any(sigma > large_sigma):
        warnings.warn('Gaussian kernel dispersion larger than {0:.1f} pixels!'.format(large_sigma))
    return VariableGaussianKernel(sigma, integral=integral).convolve(y, ye=ye)


class SpectralResolution:
    r"""
    Container class for the resolution, :math:`R =
    \lambda/\Delta\lambda`, of a spectrum.  The primary functionality is
    to determine the parameters necessary to match the resolution of one
    spectrum to another.  It can also be used as a function to
    interpolate the spectral resolution at a given wavelenth.

    Args:
        wave (numpy.ndarray):
            A 1D vector with the wavelength in angstroms.  The sampling
            may be either in linear steps of wavelength or
            :math:`\log_{10}` steps.
        sres (numpy.ndarray):
            A 1D vector with the spectral resolution, :math:`R`, sampled
            at the positions of the provided wavelength vector.
        log10 (:obj:`bool`, optional):
            Flag that the spectrum has been binned logarithmically (base
            10) in wavelength
        interp_ext (:obj:`int`, :obj:`tuple`, :obj:`str`, optional):
            The value to pass as ``fill_value`` to the interpolator,
            which defines its behavior when attempting to sample the
            spectral resolution beyond where it is defined. See
            `scipy.interpolate.interp1d`_.
        bounds_error (:obj:`bool`, optional):
            Force `scipy.interpolate.interp1d`_ to raise an exception
            if the interpolation is outside the bounds of the defined
            spectral resolution function.

    Raises:
        ValueError: Raised if *wave* is not a 1D vector or if *wave* and
            *sres* do not have the same shape.

    Attributes:
        interpolator (`scipy.interpolate.interp1d`_):
            An object used to interpolate the spectral resolution at any
            given wavelength.  The interpolation is hard-wired to be
            *linear* and its extrapolation behavior is defined by
            ``interp_ext``.  The wavelength and resolution vectors are
            held by this object for later reference if needed.
        log10 (:obj:`bool`):
            Flag that the spectrum has been binned logarithmically (base
            10) in wavelength
        c (:obj:`float`):
            The speed of light; provided by `astropy.constants`_.
        dv (:obj:`float`):
            The velocity step per pixel in km/s.  Defined using
            :func:`spectrum_velocity_scale` if :attr:`log10` is True;
            otherwise set to None.
        dw (:obj:`float`):
            The wavelength step per pixel in angstroms.  Defined as the
            wavelength step between the first two pixels if
            :attr:`log10` is False; otherwise set to None.
        min_sig (:obj:`float`):
            Minimum standard deviation allowed for the kernel used to
            match two spectral resolutions.  See
            :func:`GaussianKernelDifference`.
        sig_pd (`numpy.ndarray`):
            The standard deviation in pixels required to match the
            spectral resolution of this object to the resolution defined
            by a different :class:`SpectralResolution` object.  See
            :func:`GaussianKernelDifference`.
        sig_mask (`numpy.ndarray`):
            A *uint* vector used to identify measurements of
            :attr:`sig_pd` that should **not** be used to match the
            spectral resolution of this object to the resolution defined
            by a different :class:`SpectralResolution` object.  See
            :func:`GaussianKernelDifference`.
        sig_vo (:obj:`float`):
            A constant offset of the kernel standard deviation **in
            km/s** that has been applied to :attr:`sig_pd`.  See
            :func:`GaussianKernelDifference`.
        
    .. todo::

        - Allow it to determine if the binning is linear or geometric,
          then use the *log10* keyword to distinguish between natural
          log and :math:`log_{10}` binning.
        - Allow for more than one type of line-spread function?

    .. warning::

        The default behavior of the interpolator is to extrapolate the
        input spectral resolution vector when trying to sample from
        regions beyond where it is sampled.  Use ``interp_ext`` change
        this; see `scipy.interpolate.interp1d`_.

    """
    def __init__(self, wave, sres, log10=False, interp_ext='extrapolate', bounds_error=False):
        # Check the sizes
        if len(wave.shape) != 1:
            raise ValueError('wave must be a 1D array!')
        if wave.shape != sres.shape:
            raise ValueError('wave and sres must have the same shape!')

        # Use linear interpolation when sampling the spectral resolution
        self.interpolator = interpolate.interp1d(wave, sres, fill_value=interp_ext,
                                                 bounds_error=bounds_error)
        self.log10 = log10
        self.c = astropy.constants.c.to('km/s').value

        self.dv = spectrum_velocity_scale(wave) if log10 else None
        self.dw = None if log10 else wave[1] - wave[0]

        # No resolution matching calculated yet
        self.min_sig = None
        self.sig_pd = None
        self.sig_mask = None
        self.sig_vo = None


    def __call__(self, w):
        """Interpolate the spectral resolution at wavelength *w*."""
        return self.interpolator(w)


    def _finalize_GaussianKernelDifference(self, sig2_pd):
        r"""
        Given the calculated :math:`\sigma^2_{p,d}`, calculate and save
        the attributes :attr:`sig_pd` and :attr:`sig_mask`.  See
        :func:`GaussianKernelDifference`.
        """
#        print('finalize b')
        indx = numpy.isclose(sig2_pd, 0.0)
        nindx = numpy.invert(indx)
        self.sig_pd = sig2_pd.copy()
        self.sig_pd[nindx] = sig2_pd[nindx]/numpy.sqrt(numpy.absolute(sig2_pd[nindx]))
        self.sig_pd[indx] = 0.0
#        self.sig_mask = numpy.array(self.sig_pd < -self.min_sig).astype(numpy.uint)
        self.sig_mask = numpy.array(self.sig_pd < self.min_sig).astype(numpy.uint)
#        print('finalize e')


    def _convert_vd2pd(self, sig2_vd):
        r"""
        Convert from :math:`\sigma^2_{v,d}` to :math:`\sigma^2_{p,d}`.
        See :func:`GaussianKernelDifference`.
        """
#        print('vd2pd')
        return sig2_vd / numpy.square(self.dv) if self.log10 else \
               sig2_vd / numpy.square(self.c*self.dw/self.wave())


    def _convert_pd2vd(self, sig2_pd):
        r"""
        Convert from :math:`\sigma^2_{p,d}` to :math:`\sigma^2_{v,d}`.
        See :func:`GaussianKernelDifference`.
        """
#        print('pd2vd')
        return sig2_pd * numpy.square(self.dv) if self.log10 else \
               sig2_pd * numpy.square(self.c*self.dw/self.wave())


    def wave(self, copy=False):
        """
        Return the wavelength vector; held by :attr:`interpolator`.
        """
#        return self.interpolator._data[0]
        return self.interpolator.x.copy() if copy else self.interpolator.x


    def sres(self, copy=False):
        """
        Return the resolution vector; held by :attr:`interpolator`.
        """
#        return self.interpolator._data[1]
        return self.interpolator.y.copy() if copy else self.interpolator.y


    def instrumental_dispersion(self, w=None):
        r"""

        Return the instrumental dispersion by converting from :math:`R`
        to :math:`\sigma_{v,inst}` according to:
            
        .. math::
            
            \sigma_{v,inst} = \frac{c}{\sqrt{8\ln 2}\ R}.

        If w is None, just convert the internal interpolator values.
        Otherwise, return the values sampled at w.
        """
        if w is None:
            return numpy.ma.divide(self.c, DAPConstants.sig2fwhm*self.interpolator.y).filled(0.0)
        return numpy.ma.divide(self.c, DAPConstants.sig2fwhm*self.interpolator(w)).filled(0.0)


    def match(self, new_sres, no_offset=True, min_sig_pix=0.0):
        """
        Currently only an alias for :func:`GaussianKernelDifference`.
        """
        self.GaussianKernelDifference(new_sres, no_offset=no_offset, min_sig_pix=min_sig_pix)


    def GaussianKernelDifference(self, new_sres, no_offset=True, min_sig_pix=0.0):
        r"""
        Determine the parameters of a Gaussian kernel required to
        convert the resolution of this object to the resolution of a
        different the :class:`SpectralResolution` object, *new_sres*.

        The spectral resolution is defined as :math:`R =
        \lambda/\Delta\lambda`, where :math:`\Delta\lambda` is the FWHM
        of the spectral resolution element.  The standard deviation of
        the resolution element in angstroms is then
   
        .. math::
    
            \sigma_\lambda = \frac{\lambda}{f R}, \ \ {\rm where} \ \  f
            = \frac{{\rm FWHM_\lambda}}{\sigma_\lambda} = \sqrt{8\ln 2}.

        Assuming a Gaussian (in angstroms) line-spread function:

        .. math::

            \sigma^2_{\lambda,2} = \sigma^2_{\lambda,1} +
            \sigma^2_{\lambda,d}

        such that

        .. math::

            \sigma^2_{\lambda,d} = \left(\frac{\lambda}{f}\right)^2
            (R^{-2}_2 - R^{-2}_1)

        is the defining parameter of the Gaussian kernel needed to take
        a spectrum of resolution :math:`R_1` to one with a resolution of
        :math:`R_2`.

        For input to :func:`convolution_variable_sigma`, the
        *wavelength-dependent* parameter of the Gaussian kernel is
        converted to pixels.  This function allows for the input spectra
        to be linearly sampled in angstroms or log10(angstroms).  For
        the former (*log10=False*), 

        .. math::

            \sigma^2_{p,d} = \left(\frac{\lambda}{f\
            \delta\lambda}\right)^2 (R^{-2}_2 - R^{-2}_1)

        where :math:`\delta\lambda` is the size of the pixel in
        angstroms.  If the units are log10(angstrom) (*log10=True*), we
        approximate the velocity width of each pixel to be :math:`\delta
        v = c \ln(10.0) (\log\lambda[1]-\log\lambda[0])`, such that
    
        .. math::

            \sigma^2_{p,d} &= \left(\frac{c}{ \delta v \lambda}\right)^2
            \sigma^2_{\lambda,d} \\ &= \left(\frac{c}{ f\ \delta
            v}\right)^2 (R^{-2}_2 - R^{-2}_1)\ ;

        :math:`c` is the speed of light in km/s.

        The nominal use of this algorithm assumes :math:`R_1 \geq R_2`.
        However, in practice, :func:`convolution_variable_sigma` only
        uses a Gaussian kernel up to some minimum value of
        :math:`\epsilon_\sigma`; below this, the kernel is assumed to be
        a Delta function.  Therefore, as long as
    
        .. math::
    
            \sigma_{p,d} \equiv \sigma^2_{p,d}/\sqrt{|\sigma^2_{p,d}|}
            \geq -\epsilon_\sigma\ ,
        
        the behavior of :func:`convolution_variable_sigma` should not be
        affected.  However, in regions with
    
        .. math::
    
            \sigma_{p,d} \equiv \sigma^2_{p,d}/\sqrt{|\sigma^2_{p,d}|}
            \leq \epsilon_\sigma\ ,

        the behavior of :func:`convolution_variable_sigma` does *not*
        produce an accurate convolution!        
    
        To deal with spectral regions that do not have
        :math:`\sigma_{p,d} \geq \epsilon_\sigma`, there are three
        choices:

            (**Option 1**) trim the spectral range to only those
            spectral regions where the existing resolution is better
            than the target resolution up to this limit,
        
            (**Option 2**) match the existing resolution to the target
            resolution up to some constant offset that must be accounted
            for in subsequent analyses, or

            (**Option 3**) allow for a wavelength dependent difference
            in the spectral resolution that must be accounted for in
            subsequent analyses.

        The choice of either Option 1 or 2 is selected by setting
        *no_offset* to, respectively, True or False; Option 1 is the
        default behavior.  Currently, Option 3 is not allowed.

        For Option 1, pixels with :math:`\sigma_{p,d} <
        \epsilon_\sigma` are masked (*sigma_mask = 1*); however, the
        returned values of :math:`\sigma_{p,d}` are left unchanged.

        For Option 2, we define

        .. math::

            \sigma^2_{v,o} = | {\rm min}(0.0, {\rm min}(\sigma^2_{v,d})
            - {\rm max}(\epsilon_\sigma \delta v)^2) |

        where :math:`\delta v` is constant for the logarithmically
        binned spectrum and is wavelength dependent for the linearly
        binned spectra; in the latter case, the velocity step is
        determined for each pixel::

            _wave = self.wave()
            dv = self.c * (2.0*(_wave[1:] - _wave[0:-1]) / (_wave[1:] + _wave[0:-1]))

        If :math:`\sigma^2_{v,o} > 0.0`, it must be that :math:`{\rm
        min}(\sigma^2_{v,d}) < {\rm max}(\epsilon_\sigma \delta v)^2`,
        such that an offset should be applied.  In that case, the
        returned kernel parameters are

        .. math::

            \sigma^\prime_{v,d} = \sqrt{\sigma^2_{v,d} +
            \sigma^2_{v,o}}\ .

        with the units converted to pixels using the equations above, no
        pixels are masked, and :math:`\sqrt{\sigma^2_{v,o}}` is returned
        for the offset.  Otherwise, the offset is set to 0.

        Args:
            new_sres (:class:`SpectralResolution`):
                Spectral resolution to match to.
            no_offset (:obj:`bool`, optional):
                Force :math:`\sigma^2_{v,o} = 0` by masking regions with
                :math:`\sigma_{p,d} < \epsilon_\sigma`; i.e., the value
                of this arguments selects Option 1 (True) or Option 2
                (False).
            min_sig_pix (:obj:`float`, optional):
                Minimum value of the standard deviation allowed before
                assuming the kernel is a Delta function.

        """
#        print('kernel difference')
        # Save the minimum pixel sigma to allow
        self.min_sig = min_sig_pix

        # Interpolate the new spectral resolution vector at the wavelengths
        # of the input spectral resolution
        _wave = self.wave()
        _sres = self.sres()
        interp_sres = new_sres(_wave)

        # Determine the variance (in angstroms) of Gaussian needed to match
        # input resolution to the new values
#        print('sig2_wd')
        sig2_wd = numpy.square(_wave/DAPConstants.sig2fwhm) \
                  * (1.0/numpy.square(interp_sres) - 1.0/numpy.square(_sres))

#        pyplot.plot(_wave, _sres)
#        pyplot.plot(_wave, interp_sres)
#        pyplot.show()

#        print('sig2_vd')
        # Convert to km/s
        sig2_vd = numpy.square(self.c/_wave) * sig2_wd

        # Option 1:
        if no_offset:
            # Convert the variance to pixel coordinates
            sig2_pd = sig2_vd / numpy.square(self.dv) if self.log10 else \
                      sig2_wd / numpy.square(self.dw)
            # No offset applied
            self.sig_vo = 0.0

        # Option 2:
        else:
            # 1% fudge so pixel at min_sig is not masked!
            fudge = 1.01
            
            # Calculate the velocity step of each pixel
            dv = self.c * (2.0*(_wave[1:] - _wave[0:-1]) / (_wave[1:] + _wave[0:-1]))
            # Get the needed *velocity* offset (this is the square)
            self.sig_vo = fudge*numpy.absolute(min(0.0, numpy.amin(sig2_vd)
                                                  - numpy.square(self.min_sig * numpy.amax(dv))))
            # Apply it if it's larger than 0
            if self.sig_vo > 0:
                sig2_vd += self.sig_vo
                self.sig_vo = numpy.sqrt(self.sig_vo)

            # Convert the variance to pixel coordinates
            sig2_pd = self._convert_vd2pd(sig2_vd)

        self._finalize_GaussianKernelDifference(sig2_pd)

#        print('kernel difference done')

#    def ZeroGaussianKernelDifference(self, min_sig_pix=0.0):
#        self.min_sig = min_sig_pix
#        sig2_pd = numpy.zeros(self.wave().shape, dtype=numpy.float64)
#        self._finalize_GaussianKernelDifference(sig2_pd)


    def offset_GaussianKernelDifference(self, offset):
        r"""
        If the properties required to match the resolution of one
        spectrum to another has already been calculated (see
        :func:`GaussianKernelDifference`), this allows for one to apply
        an additional offset.  The additional offset **must** be in km/s
        (not pixels).

        The offset is applied in quadrature; however, the offset can be
        negative such that one can reduce :attr:`sig_pd`.  Once
        converted to km/s, the offset is applied by calculating:

        .. math::
        
            \sigma^{\prime\ 2}_{v,d} = \sigma^{2}_{v,d} +
            \sigma_{off}|\sigma_{off}|\ .

        :attr:`sig_vo` is adjusted in the same way, and the change in
        :math:`\sigma^{\prime\ 2}_{v,d}` is then propagated to
        :attr:`sig_pd` and :attr:`sig_mask`.
        
        Args:
            offset (float): Value of the standard deviation in km/s to
                add in quadrature to a previously calculated
                :attr:`sig_pd`.
        
        Raises:
            ValueError: Raised if the kernel properties have not yet
                been defined.
        """
        if self.min_sig is None or self.sig_pd is None or self.sig_mask is None \
                or self.sig_vo is None:
#            print('WARNING: No kernel difference yet defined.  Assuming 0.')
#            self.ZeroGaussianKernelDifference()
            raise ValueError('No kernel defined yet.  Run GaussianKernelDifference first.')
        if numpy.isclose(offset,0.0):
            return
        off2 = offset*numpy.absolute(offset)
        sig2_vo = self.sig_vo*numpy.absolute(self.sig_vo) + off2
        self.sig_vo = 0.0 if numpy.isclose(sig2_vo, 0.0) \
                          else sig2_vo/numpy.sqrt(numpy.absolute(sig2_vo))
        sig2_vd = self._convert_pd2vd(self.sig_pd*numpy.absolute(sig_pd)) + off2
        self._finalize_GaussianKernelDifference(self._convert_vd2pd(sig2_vd))


    def adjusted_resolution(self, indx=None):
        r"""

        Return the resolution that should result from applying
        :func:`convolution_variable_sigma` to the spectrum associated
        with this spectral resolution object using :attr:`sigma_pd`.
        I.e., calculate:

        .. math::

            R_2 = \left[ \left(\frac{f}{c}\right)^2 \sigma^2_{v,d} +
            R^{-2}_1\right]^{-1/2}\ . 

        Args:
            indx (tuple): (**Optional**) Selection tuple used to return
                a subset of the full resolution vector.

        Returns:
            numpy.ndarray: The (full or selected) vector with the
            adjusted resolution.

        .. todo::
            Allow to reset the resolution of this object to the adjusted
            resolution and reset the kernel variables to None.

        """
        if indx is None:
            return 1.0/numpy.sqrt( numpy.square(DAPConstants.sig2fwhm/self.c) \
                                   * self._convert_pd2vd(self.sig_pd*numpy.absolute(self.sig_pd)) \
                                   + 1.0/numpy.square(self.sres()) )

        return 1.0/numpy.sqrt( numpy.square(DAPConstants.sig2fwhm/self.c) \
                            * (self._convert_pd2vd(self.sig_pd*numpy.absolute(self.sig_pd)))[indx] \
                               + 1.0/numpy.square(self.sres()[indx]) )


def match_spectral_resolution(wave, flux, sres, new_sres_wave, new_sres, ivar=None, mask=None, 
                              min_sig_pix=0.0, no_offset=True, variable_offset=False, log10=False,
                              new_log10=False, integral=True, quiet=False):
    r"""
    Adjust the existing spectral resolution of a spectrum to a **lower**
    resolution as best as possible.  The primary functionality is in
    :class:`SpectralResolution`, which determines the Gaussian kernel
    parameters needed to match the resolution, and
    :func:`convolve_variable_sigma`, which actually performs the
    convolution to match the resolution.

    In particular, see
    :func:`SpectralResolution.GaussianKernelDifference` for a
    description of how the kernel parameters are determined and how
    regions where the target resolution is **higher** than the existing
    resolution.  In this case, one of the options is to adopt an offset
    of the resolution (in km/s) that could be corrected for in
    subsequent analysis.  In this case, setting *variable_offset* to
    True allows the offset to be different for all input spectra.  If
    one expects to combine the spectra, the default behavior should be
    used, forcing all the spectra to have a constant offset.

    Args:
        wave (numpy.ndarray): A 1D or 2D (:math:`N_{\rm spec}\times
            N_{\rm pix}`) array with the wavelength in angstroms for a
            set of spectra.  The sampling may be either in linear steps
            of wavelength or :math:`\log_{10}` steps, as set using
            *log10*.
        flux (numpy.ndarray): A 1D or 2D (:math:`N_{\rm spec}\times
            N_{\rm pix}`) array with the flux sampled at the provided
            wavelengths.
        sres (numpy.ndarray): A 1D or 2D (:math:`N_{\rm spec}\times
            N_{\rm pix}`) array with the spectral resolution, :math:`R`,
            at the provided wavelengths.
        new_sres_wave (numpy.ndarray): A 1D vector with the wavelength
            in angstroms at which the new resolution of the input
            spectra has been sampled.  The sampling may be either in
            linear steps of wavelength or :math:`\log_{10}` steps, as
            set using *new_log10*.
        new_sres (numpy.ndarray): A 1D vector with the new resolution
            for the input spectra.
        ivar (numpy.ndarray): (**Optional**) A 1D or 2D (:math:`N_{\rm
            spec}\times N_{\rm pix}`) array with the inverse variance of
            the flux sampled at the provided wavelengths.  This vector
            is used to estimate the noise in the resolution-matched
            spectra.

            .. warning::
                The accuracy of the errors still remain untested!
            
        mask (numpy.ndarray): (**Optional**) A 1D or 2D (:math:`N_{\rm
            spec}\times N_{\rm pix}`) array with a *uint* mask for the
            flux sampled at the provided wavelengths.
        no_offset (bool): (**Optional**) Force :math:`\sigma^2_{v,o} =
            0` by masking regions with :math:`\sigma_{p,d} <
            -\epsilon_\sigma`; i.e., the value of this arguments selects
            Option 1 (True) or Option 2 (False).  See
            :func:`SpectralResolution.GaussianKernelDifference`.
        min_sig_pix (float): (**Optional**) Minimum value of the
            standard deviation in pixels allowed before assuming the
            kernel is a Delta function.
        variable_offset (bool): (**Optional**) Flag to allow the offset
            applied to each spectrum (when the input contains more than
            one spectraum) to be tailored to each spectrum.  Otherwise
            (*variable_offset=False*) the offset is forced to be the
            same for all spectra.
        log10 (bool): (**Optional**) Flag that the spectrum has been
            binned logarithmically (base 10) in wavelength
        new_log10 (bool): (**Optional**) Flag that the coordinates of
            the new spectral resolution are  spectrum as been binned
            logarithmically (base 10) in wavelength.

    Returns: 
        numpy.ndarray: Five objects are returned:

            - A 1D or 2D (:math:`N_{\rm spec}\times N_{\rm pix}`) array
              with the resolution-matched flux sampled at the input
              wavelengths.
            - A 1D or 2D (:math:`N_{\rm spec}\times N_{\rm pix}`) array
              with the spectral resolution, :math:`R`, of the
              resolution-matched spectra at the provided wavelengths.
            - A 1D vector with any constant offset in resolution **in
              km/s** between the targetted value and the result.  See
              :func:`SpectralResolution.GaussianKernelDifference`.
            - A 1D or 2D (:math:`N_{\rm spec}\times N_{\rm pix}`) array
              with a *uint* mask for the resolution-matched flux sampled
              at the input wavelengths.  This is returned regardless of
              whether an input mask was provided.  Any pixel that had a
              resolution that was lower than the target resolution (up
              to some tolerance defined by *min_sig_pix*) is returned as
              masked.
            - A 1D or 2D (:math:`N_{\rm spec}\times N_{\rm pix}`) array
              with the inverse variance of the resolution-matched flux
              sampled at the input wavelengths.  If *ivar*
              is not provided, a 'None' returned as the last element

    Raises:
        ValueError: Raised if:

            - the input *wave* array is 2D and the *sres* array is not;
              a 1D wavelength array is allowed for a 2D *sres* array but
              not vice versa

            - the number of spectral pixels in *wave*, *flux*, and
              *sres* is not the same

            - the shape of the *flux*, *mask* (if provided), and *ivar*
              (if provided) are not the same

            - the shape of the *new_sres_wave* and *new_sres* arrays
              are not the same and/or not 1D

    .. todo::

        - Add interp_ext != 'extrapolate' option?
        - Better way to use warnings?

    """
    # Check the dimensionality of wave and sres
    wave_matrix = len(wave.shape) == 2
    sres_matrix = len(sres.shape) == 2
    if wave_matrix and not sres_matrix:
        raise ValueError('If input wavelength array is 2D, the spectral resolution array must' \
                         ' also be 2D')

    # Check the shapes
    if (wave_matrix == sres_matrix and wave.shape != sres.shape) or \
       (not wave_matrix and sres_matrix and wave.shape[0] != sres.shape[1]):
        raise ValueError('Input spectral resolution and coordinate arrays must have the same' \
                         ' number of spectral channels!')
    if (wave_matrix and wave.shape != flux.shape) or \
       (not wave_matrix and len(flux.shape) == 2 and wave.shape[0] != flux.shape[1]) or \
       (not wave_matrix and len(flux.shape) == 1 and wave.shape != flux.shape):
        raise ValueError('Input flux and coordinate arrays must have the same number of' \
                         ' spectral channels!')
    if (mask is not None and mask.shape != flux.shape):
        raise ValueError('Input flux and mask arrays must have the same shape!')
    if (ivar is not None and ivar.shape != flux.shape):
        raise ValueError('Input flux and ivar arrays must have the same shape!')

    if len(sres.shape) > len(flux.shape):
        raise ValueError('Shape of the spectral resolution array must be <= to the flux array.')
        
    if len(new_sres_wave.shape) != 1 or len(new_sres.shape) != 1:
        raise ValueError('New spectral resolution and coordinate arrays must be 1D!')
    if new_sres_wave.shape != new_sres.shape:
        raise ValueError('New spectral resolution and coordinate arrays must have the same shape!')

    # Raise a warning if the new_sres vector will have to be
    # extrapolated for the input wavelengths
    if not quiet and (numpy.amin(wave) < new_sres_wave[0] or numpy.amax(wave) > new_sres_wave[-1]):
        warnings.warn('Mapping to the new spectral resolution will require extrapolating the' \
                      ' provided input vectors!')

    # Initialize some variables
    nspec = 1 if len(flux.shape) == 1 else flux.shape[0]
    nsres = 1 if len(sres.shape) == 1 else sres.shape[0]
    if sres_matrix and nspec != nsres:
        raise ValueError('For 2D matrices, number of spectral resolution vectors must match the ' \
                         'number of spectra.')
    spec_dim = len(flux.shape)
    sres_dim = len(sres.shape)
    sigma_offset = numpy.zeros(nspec, dtype=numpy.float64)
    new_res = SpectralResolution(new_sres_wave, new_sres, log10=new_log10,
                                 interp_ext=(new_sres[0],new_sres[-1]))

#    pyplot.plot(new_sres_wave, new_sres)
#    pyplot.show()
#    exit()

    res = numpy.empty(nspec, dtype=object)

    # Get the kernel parameters necessary to match all spectra to the
    # new resolution
    if nsres == 1 and sres_dim == 1:
        res[0] = SpectralResolution(wave, sres, log10=log10, interp_ext=(sres[0], sres[-1]))
        res[0].match(new_res, no_offset=no_offset, min_sig_pix=min_sig_pix)
        sigma_offset[0] = res[0].sig_vo
        for i in range(1,nspec):
            res[i] = res[0]
#        pyplot.plot(wave, res[0].sig_pd)
#        pyplot.show()
    else:
        for i in range(0,nsres):
            _wave = wave[i,:].ravel() if wave_matrix else wave
            _sres = sres[i,:].ravel() if sres_matrix else sres
            res[i] = SpectralResolution(_wave, _sres, log10=log10,
                                        interp_ext=(sres[i,0],sres[i,-1]))
            res[i].match(new_res, no_offset=no_offset, min_sig_pix=min_sig_pix)
            sigma_offset[i] = res[i].sig_vo
#            pyplot.plot(_wave, res[i].sig_pd)
#            pyplot.plot(_wave, numpy.sqrt(res[i]._convert_pd2vd(numpy.square(res[i].sig_pd))))
#            pyplot.show()

    # Force all the offsets to be the same, if requested
    if not no_offset and not variable_offset:
        common_offset = numpy.max(sigma_offset)
        offset_diff = numpy.sqrt( numpy.square(common_offset) - numpy.square(sigma_offset))
        for r, o in zip(res,offset_diff):
            r.offset_GaussianKernelDifference(o)

    # Perform the convolutions
    out_flux = flux.copy()
    out_ivar = None if ivar is None else numpy.ma.MaskedArray(ivar.copy())
    noise = None if ivar is None else numpy.ma.sqrt(1.0/out_ivar)
    out_sres = sres.copy()
    _mask = numpy.zeros(flux.shape, dtype=numpy.uint) if mask is None else mask
    out_mask = _mask.copy()

#    print('test div by zero')
    if nspec == 1 and spec_dim == 1:
        sig_kernel = res[0].sig_pd.copy()
        sig_kernel[sig_kernel < min_sig_pix] = 0.0
        indx = res[0].sig_pd > min_sig_pix
        try:
            if ivar is None:
                out_flux[indx] = convolution_variable_sigma(flux, sig_kernel,
                                                            integral=integral)[indx]
            else:
                oflux, oivar = convolution_variable_sigma(flux, sig_kernel,
                                                          ye=None if ivar is None else noise[indx],
                                                          integral=integral)

                out_flux[indx] = oflux[indx]
                out_ivar[indx] = oivar[indx]
                del oflux, oivar
            out_sres[indx] = res[0].adjusted_resolution(indx=indx)
            out_mask = numpy.array((res[0].sig_mask == 1) | (_mask == 1)).astype(numpy.uint)
        except ValueError as e:
            warnings.warn('Encountered ValueError: {0} ; continuing but resolution is NOT '
                          'changed and mask is set.'.format(e))
            out_mask = numpy.ones(flux.shape, dtype=numpy.uint)
    else:
        for i in range(0,nspec):
            if not quiet:
                print('Matching resolution: {0}/{1}'.format(i+1,nspec), end='\r')
            try:
#                indx = numpy.where(res[i].sig_pd > min_sig_pix)
                sig_kernel = res[i].sig_pd.copy()
                sig_kernel[sig_kernel < min_sig_pix] = 0.0
                indx = res[i].sig_pd > min_sig_pix
                if ivar is None:
                    out_flux[i,indx] = convolution_variable_sigma(flux[i,:], sig_kernel,
                                                                  integral=integral)[indx]
                else:
                    oflux, oivar = convolution_variable_sigma(flux[i,:], sig_kernel,
                                                              ye=None if ivar is None 
                                                                    else noise[i,:].ravel(),
                                                              integral=integral)
                    out_flux[i,indx] = oflux[indx]
                    out_ivar[i,indx] = oivar[indx]
                    del oflux, oivar
                out_mask[i,:] = numpy.array((res[i].sig_mask == 1) \
                                        | (_mask[i,:] == 1)).astype(numpy.uint)
                if len(out_sres.shape) == 1 and i == 0:
                    out_sres[indx] = res[i].adjusted_resolution(indx=indx)
                else:
                    out_sres[i,indx] = res[i].adjusted_resolution(indx=indx)
            except ValueError as e:
                warnings.warn('Encountered ValueError: {0} ; continuing but resolution is NOT '
                              'changed and mask is set.'.format(e))
                out_mask[i,:] = numpy.ones(flux.shape[1], dtype=numpy.uint)
        if not quiet:
            print('Matching resolution: {0}/{0}'.format(nspec))

    # TODO: Add this functionality from the IDL version?
    #
    # Finally, the code masks a number of pixels at the beginning and
    # end of the spectra to remove regions affected by errors in the
    # convolution due to the censoring of the data.  The number of
    # pixels is the FWHM of the largest Gaussian applied in the
    # convolution: ceil(sig2fwhm*max(diff_sig_w)/dw).  This is currently
    # hard-wired and should be tested.

    if ivar is not None:
        out_ivar = numpy.square(1.0/out_ivar)
        # When returning out_ivar, convert it to a normal array
        return out_flux, out_sres, sigma_offset, out_mask, numpy.asarray(out_ivar)
    return out_flux, out_sres, sigma_offset, out_mask, None
    

