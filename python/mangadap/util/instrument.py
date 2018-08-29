# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Provides a set of functions to handle instrumental sampling and
resolution effects.

*License*:
    Copyright (c) 2015, SDSS-IV/MaNGA Pipeline Group
        Licensed under BSD 3-clause license - see LICENSE.rst

*Source location*:
    $MANGADAP_DIR/python/mangadap/util/instrument.py

*Imports and python version compliance*:
    ::

        from __future__ import division
        from __future__ import print_function
        from __future__ import absolute_import
        from __future__ import unicode_literals

        import sys
        if sys.version > '3':
            long = int

        import warnings
        import numpy
        from scipy import integrate
        from scipy.interpolate import InterpolatedUnivariateSpline
        from scipy.special import erf
        import astropy.constants

        from .constants import constants

*Revision history*:
    | **27 May 2015**: Original implementation by K. Westfall (KBW)
        based on downgrader_MANGA.f provided by D.  Thomas, O. Steele,
        D.  Wilkinson, D. Goddard.
    | **13 Jun 2015**: D.Wilkinson edit to not calculate unimportant
        convolution terms -> runs 5x faster.
    | **03 Feb 2016**: (KBW) Generalized M. Cappellari's
        :func:`log_rebin` to :func:`resample_vector` to allow for linear
        interpolation of flux density across the pixel.
    | **04 Feb 2016**: (KBW) Further fixes to :func:`resample_vector`.
        Added :func:`spectral_coordinate_step`, and propagated change to
        :func:`spectrum_velocity_scale`
    | **21 Apr 2016**: (KBW) It's the Queen's 90th birthday!  Removed
        log10 keyword from :func:`spectrum_velocity_scale`.
    | **20 May 2016**: (KBW) Corrected match between number of pixels
        and output range computed in :func:`resample_vector_npix`; now
        returns an adjusted range to make sure that the sampling and
        range results in an exact integer number of pixels.
    | **05 Jul 2016**: (KBW) To avoid confusion, commented out
        log_rebin.
    | **25 Oct 2016**: (KBW) Modified :func:`spectral_coordinate_step`
        to be a mean over the full spectrum to avoid numerical precision
        errors.
    | **06 Apr 2017**: (KBW) Add :func:`angstroms_per_pixel`.
    | **30 Aug 2017**: (KBW) Add :func:`resample1d`;
        :func:`resample_vector` should be deprecated.
    | **27 Sep 2017**: (KBW) Added `integral` keyword to
        :func:`match_spectral_resolution` so that it can be passed to
        :func:`convolution_variable_sigma`.
    | **19 Jul 2018**: (KBW) Fixed a bug in
        :class:`VariableGaussianKernel` when constructing the kernel
        based on the error function.
"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import sys
if sys.version > '3':
    long = int

import warnings
import numpy
from scipy import integrate, interpolate
#from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.special import erf
import astropy.constants

from .constants import DAPConstants

from matplotlib import pyplot

def spectral_coordinate_step(wave, log=False, base=10.0):
    """
    Return the sampling step for the input wavelength vector.

    If the sampling is logarithmic, return the change in the logarithm
    of the wavelength; otherwise, return the linear step in angstroms.

    Args: 
        wave (numpy.ndarray): Wavelength coordinates of each spectral
            channel in angstroms.
        log (bool): (**Optional**) Input spectrum has been sampled
            geometrically.
        base (float): (**Optional**) If sampled geometrically, the
            sampling is done using a logarithm with this base.  For
            natural logarithm, use numpy.exp(1).

    Returns:
        float: Spectral sampling step in either angstroms (log=False) or
        the step in log(angstroms).
    """
    dw = numpy.diff(numpy.log(wave))/numpy.log(base) if log else numpy.diff(wave)
    if numpy.any( numpy.absolute(numpy.diff(dw)) > 100*numpy.finfo(dw.dtype).eps):
#        pyplot.plot(dw[:-1], numpy.absolute(numpy.diff(dw)))
#        pyplot.show()
        raise ValueError('Wavelength vector is not uniformly sampled to numerical accuracy.')
    return numpy.mean(dw)


def spectrum_velocity_scale(wave):
    """
    Determine the velocity sampling of an input wavelength coordinate
    vector.
    
    .. note::
        The wavelength vector is assumed to be geometrically sampled!
        However, the input units expected to be in angstroms, not, e.g.,
        log(angstrom).

    Args: 
        wave (numpy.ndarray): Wavelength coordinates of each spectral
            channel in angstroms.  It is expected that the spectrum has
            been sampled geometrically

    Returns:
        float: Velocity scale of the spectrum in km/s.

    """
    return astropy.constants.c.to('km/s').value*spectral_coordinate_step(wave, log=True,
                                                                         base=numpy.exp(1.))


def angstroms_per_pixel(wave, log=False, base=10.0, regular=True):
    """
    Return a vector with the angstroms per pixel at each channel.

    When `regular=True`, the function assumes that the wavelengths are
    either sampled linearly or geometrically.  Otherwise, it calculates
    the size of each pixel as the difference between the wavelength
    coordinates.  The first and last pixels are assumed to have a width
    as determined by assuming the coordinate is at it's center.

    Args:
        wave (numpy.ndarray): (Geometric) centers of the spectrum pixels
            in angstroms.
        log (numpy.ndarray): (**Optional**) The vector is geometrically
            sampled.
        base (float): (**Optional**) Base of the logarithm used in the
            geometric sampling.
        regular (bool): (**Optional**) The vector is regularly sampled.
            

    Returns:
        numpy.ndarray: The angstroms per pixel.
    """
    if regular:
        ang_per_pix = spectral_coordinate_step(wave, log=log, base=base)
        if log:
            ang_per_pix *= wave*numpy.log(base)
    else:
        ang_per_pix = numpy.diff([(3*wave[0]-wave[1])/2] 
                                    + ((wave[1:] + wave[:-1])/2).tolist()
                                    + [(3*wave[-1]-wave[-2])/2])
    return ang_per_pix


class convolution_integral_element:
    """
    Support class for variable sigma convolution.  See
    :func:`convolution_variable_sigma`.

    OUT OF DATE; DO NOT USE

    Args:
        y (numpy.ndarray): Vector to convolve
        sigma (numpy.ndarray): Coordinate-dependent standard deviation of the
            Gaussian kernel
        ye (numpy.ndarray): (**Optional**) Error in the vector to
            convolve

    Raises:
        ValueError: Raised if *y* is not a 1D vector, or if the shape of
            *y* and *sigma* (and *ye* if provided) are different.

    Attributes:
        x (numpy.ndarray): Pixel coordinate vector
        y (numpy.ndarray): Vector to convolve
        ye (numpy.ndarray): Error in the vector to convolve
        sigma (numpy.ndarray): Coordinate-dependent standard deviation of the
            Gaussian kernel
        norm (numpy.ndarray): Gaussian normalization; calculated once for
            efficiency

    .. todo::

        - Allow to switch to pixel sampled Gaussian kernel?

    """
    def __init__(self, y, sigma, ye=None):
        if len(y.shape) != 1:
            raise ValueError('y must be a 1D array!')
        if y.shape != sigma.shape:
            raise ValueError('y and sigma must have the same shape!')
        if ye is not None and ye.shape != y.shape:
            raise ValueError('y and ye must have the same shape!')
        self.x = numpy.arange(sigma.size, dtype=numpy.float64)
        self.y = y
        self.ye = ye
        self.sigma = sigma
        self.norm = numpy.sqrt(2.0*numpy.pi)*self.sigma


    def _get_kernel(self, xc):
        """Calculate the kernel vector when centered at *xc*.

        .. todo::

            - Function is about 30% slower when using erf() as opposed
              to exp().  erf() needed when sigma is small, but may be
              efficient to include some decision on when it's safe to
              use the quick way.
        
        """

        d = (self.x-xc)
        gf = numpy.square(d/self.sigma) 
        close_value = numpy.where(gf < 50.0)
#        outkern = numpy.exp(-0.5*gf[close_value])/self.norm[close_value]
        outkern = (erf((d[close_value]+0.5)/numpy.sqrt(2)/self.sigma[close_value])
                    - erf((d[close_value]-0.5)/numpy.sqrt(2)/self.sigma[close_value]))/2.
        return close_value, outkern


#    def _get_kernel(self, xc):
#        """Calculate the kernel vector when centered at *xc*."""
#        return numpy.exp(-0.5*numpy.square((self.x-xc)/self.sigma))/self.norm


    def __call__(self, xc):
        """
        Calculates the weighted mean of :attr:`y`, where the weights are
        defined by a Gaussian with standard deviation :attr:`sigma` and
        centered at xc.

        Args:
            xc (float): Center for the Gaussian weighting function

        Returns:
            float: The weighted mean of :attr:`y`
        """
#        kernel = self._get_kernel(xc)
#        return numpy.sum(self.y*kernel) / numpy.sum(kernel)
        close_array, kernel = self._get_kernel(xc)
        return numpy.sum(self.y[close_array]*kernel) / numpy.sum(kernel)

        # Need to test equivalence and speed of these two implementations
#        return integrate.simps(self.y*kernel) / integrate.simps(kernel)


    def error(self, xc):
        """
        Calculates the error in the weighted mean of :attr:`y` using
        nominal error propagation.  The weights are defined by a
        Gaussian with standard deviation :attr:`sigma` and centered at
        xc.

        Args:
            xc (float): Center for the Gaussian weighting function

        Returns:
            float: The error in the weighted mean of :attr:`y`
        """
        close_array, kernel = self._get_kernel(xc)
        return numpy.sqrt(numpy.sum(numpy.square(self.ye[close_array]*kernel)) / numpy.sum(kernel))
#        kernel = self._get_kernel(xc)
#        return numpy.sqrt(numpy.sum(numpy.square(self.ye*kernel)) / numpy.sum(kernel))


class VariableGaussianKernel:
    """
    Support class for variable sigma convolution.  See
    :func:`convolution_variable_sigma`.

    Stolen from M. Cappellari's gaussian_filter1d function.

    Args:
        y (numpy.ndarray): Vector to convolve
        sigma (numpy.ndarray): Coordinate-dependent standard deviation of the
            Gaussian kernel
        ye (numpy.ndarray): (**Optional**) Error in the vector to
            convolve

    Raises:
        ValueError: Raised if *y* is not a 1D vector, or if the shape of
            *y* and *sigma* (and *ye* if provided) are different.

    Attributes:
        x (numpy.ndarray): Pixel coordinate vector
        y (numpy.ndarray): Vector to convolve
        ye (numpy.ndarray): Error in the vector to convolve
        sigma (numpy.ndarray): Coordinate-dependent standard deviation of the
            Gaussian kernel
        norm (numpy.ndarray): Gaussian normalization; calculated once for
            efficiency

    .. todo::

        - Allow to switch to pixel sampled Gaussian kernel?

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
        self._check_shape(y, ye=ye)

        # Create m copies of the shifted input function
        a = self._create_a(y)
        if ye is None:
            return numpy.sum(a*self.kernel, axis=0)

        # Construct the error
        ae = self._create_a(numpy.square(ye))
        return numpy.sum(a*self.kernel, axis=0), numpy.sqrt(numpy.sum(ae*self.kernel, axis=0))


def convolution_variable_sigma(y, sigma, ye=None, integral=False):
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

    Args:
        y (numpy.ndarray): A uniformly sampled function to convolve.
        sigma (numpy.ndarray): The standard deviation of the Gaussian
            kernel sampled at the same positions as *y*.  The units of
            *sigma* **must** be in pixels.
        ye (numpy.ndarray): (**Optional**) Errors in the function
            :math:`y(x)`.

    Returns:
        numpy.ndarray: Arrays with the convolved function :math:`(y\ast
        g)(x)` sampled at the same positions as the input :math:`x`
        vector and its error.  The second array will be returned as None
        if the error vector is not provided.
    """
#    kernel = convolution_integral_element(y,sigma,ye=ye)
#    conv = numpy.array([kernel(x) for x in kernel.x])
#    if ye is None:
#        return conv
#    return conv, numpy.array([kernel.error(x) for x in kernel.x])

    return VariableGaussianKernel(sigma, integral=integral).convolve(y,ye=ye)


class SpectralResolution:
    r"""
    
    Container class for the resolution, :math:`R =
    \lambda/\Delta\lambda`, of a spectrum.  The primary functionality is
    to determine the parameters necessary to match the resolution of one
    spectrum to another.  It can also be used as a function to
    interpolate the spectral resolution at a given wavelenth.

    Args:
        wave (numpy.ndarray): A 1D vector with the wavelength in
            angstroms.  The sampling may be either in linear steps of
            wavelength or :math:`\log_{10}` steps.
        sres (numpy.ndarray): A 1D vector with the spectral resolution,
            :math:`R`, sampled at the positions of the provided
            wavelength vector.
        log10 (bool): (**Optional**) Flag that the spectrum has been
            binned logarithmically (base 10) in wavelength
        interp_ext (int or str): (**Optional**) The value to pass as
            *ext* to the interpolator, which defines its behavior when
            attempting to sample the spectral resolution beyond where it
            is defined.  See
            `scipy.interpolate.InterpolatedUnivariateSpline`_.  Default
            is to extrapolate.

    Raises:
        ValueError: Raised if *wave* is not a 1D vector or if *wave* and
            *sres* do not have the same shape.

    Attributes:
        interpolator
            (`scipy.interpolate.InterpolatedUnivariateSpline`_): An
            object used to interpolate the spectral resolution at any
            given wavelength.  The interpolation is hard-wired to be
            **linear** and its extrapolation behavior is defined by
            *interp_ext*.  The wavelength and resolution vectors are
            held by this object for later reference if needed.
        log10 (bool): Flag that the spectrum has been binned
            logarithmically (base 10) in wavelength
        c (float): The speed of light; provided by `astropy.constants`_.
        dv (float): The velocity step per pixel in km/s.  Defined using
            :func:`spectrum_velocity_scale` if :attr:`log10` is True;
            otherwise set to None.
        dw (float): The wavelength step per pixel in angstroms.  Defined
            as the wavelength step between the first two pixels if
            :attr:`log10` is False; otherwise set to None.
        min_sig (float): Minimum standard deviation allowed for the
            kernel used to match two spectral resolutions.  See
            :func:`GaussianKernelDifference`.
        sig_pd (numpy.ndarray): The standard deviation in pixels
            required to match the spectral resolution of this object to
            the resolution defined by a different
            :class:`SpectralResolution` object.  See
            :func:`GaussianKernelDifference`.
        sig_mask (numpy.ndarray): A *uint* vector used to identify
            measurements of :attr:`sig_pd` that should **not** be used
            to match the spectral resolution of this object to the
            resolution defined by a different
            :class:`SpectralResolution` object.  See
            :func:`GaussianKernelDifference`.
        sig_vo (float): A constant offset of the kernal standard
            deviation **in km/s** that has been applied to
            :attr:`sig_pd`.  See :func:`GaussianKernelDifference`.
        
    .. todo::

        - Allow it to determine if the binning is linear or geometric,
          then use the *log10* keyword to distinguish between natural
          log and :math:`log_{10}` binning.
        - Allow for more than one type of line-spread function?

    .. warning::

        The default behavior of the interpolator is to extrapolate the
        input spectral resolution vector when trying to sample from
        regions beyond where it is sampled.  Use *interp_ext* change
        this; see `scipy.interpolate.InterpolatedUnivariateSpline`_.

    .. _scipy.interpolate.InterpolatedUnivariateSpline: http://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.InterpolatedUnivariateSpline.html
    .. _astropy.constants: http://docs.astropy.org/en/stable/constants/index.html

    """
    def __init__(self, wave, sres, log10=False, interp_ext='extrapolate'):
        # Check the sizes
        if len(wave.shape) != 1:
            raise ValueError('wave must be a 1D array!')
        if wave.shape != sres.shape:
            raise ValueError('wave and sres must have the same shape!')

        # k=1; always use linear interpolation
#        if sys.version < '3':
#            print('WARNING: InterpolatedUnivariateSpline may have different behavior in python2!')
#            self.interpolator = InterpolatedUnivariateSpline(wave, sres, k=1)
#        else:
#            self.interpolator = InterpolatedUnivariateSpline(wave, sres, k=1, ext=interp_ext)
        self.interpolator = interpolate.interp1d(wave, sres, fill_value=interp_ext)
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
            new_sres (:class:`SpectralResolution`): Spectral resolution
                to match to.
            no_offset (bool): (**Optional**) Force :math:`\sigma^2_{v,o}
                = 0` by masking regions with :math:`\sigma_{p,d} <
                \epsilon_\sigma`; i.e., the value of this arguments
                selects Option 1 (True) or Option 2 (False).
            min_sig_pix (float): (**Optional**) Minimum value of the
                standard deviation allowed before assuming the kernel is
                a Delta function.
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
    new_res = SpectralResolution(new_sres_wave, new_sres, log10=new_log10)

#    pyplot.plot(new_sres_wave, new_sres)
#    pyplot.show()
#    exit()

    res = numpy.empty(nspec, dtype=object)

    # Get the kernel parameters necessary to match all spectra to the
    # new resolution
    if nsres == 1 and sres_dim == 1:
        res[0] = SpectralResolution(wave, sres, log10=log10)
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
            res[i] = SpectralResolution(_wave, _sres, log10=log10)
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
            out_mask = numpy.array((res[0].sig_mask == 1) | (mask == 1)).astype(numpy.uint)
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
    

#def log_rebin(lamRange, spec, oversample=None, velscale=None, flux=False, log10=False,
#              newRange=None, wave_in_ang=False, unobs=0.0):
#    """
#    .. note::
#     
#        Copyright (C) 2001-2014, Michele Cappellari
#        E-mail: cappellari_at_astro.ox.ac.uk
#     
#        This software is provided as is without any warranty whatsoever.
#        Permission to use, for non-commercial purposes is granted.
#        Permission to modify for personal or internal use is granted,
#        provided this copyright and disclaimer are included unchanged at
#        the beginning of the file. All other rights are reserved.
#
#    Logarithmically rebin a spectrum, while rigorously conserving the
#    flux.  Basically the photons in the spectrum are simply
#    ridistributed according to a new grid of pixels, with non-uniform
#    size in the spectral direction.
#    
#    This routine makes the `standard' zero-order assumption that the
#    spectrum is *constant* within each pixels. It is possible to perform
#    log-rebinning by assuming the spectrum is represented by a
#    piece-wise polynomial of higer degree, while still obtaining a
#    uniquely defined linear problem, but this reduces to a deconvolution
#    and amplifies noise.
#
#    .. warning::
#
#        This assumption can be poor for sharp features in the spectrum.
#        Beware if resampling spectra with strong, marginally sampled
#        features!
#    
#    This same routine can be used to compute approximate errors of the
#    log-rebinned spectrum. To do this type the command
#    
#    >>> err2New, logLam, velscale = log_rebin(lamRange, numpy.square(err))
#    
#    and the desired errors will be given by numpy.sqrt(err2New).
#    
#    .. warning::
#    
#        This rebinning of the error-spectrum is very *approximate* as it
#        does not consider the correlation introduced by the rebinning!
#
#    Args:
#
#        lamRange (numpy.ndarray): two elements vector containing the
#            central wavelength of the first and last pixels in the
#            spectrum, which is assumed to have constant wavelength
#            scale! E.g. from the values in the standard FITS keywords:
#            LAMRANGE = CRVAL1 + [0,CDELT1*(NAXIS1-1)].  It must be
#            LAMRANGE[0] < LAMRANGE[1].
#        spec (numpy.ndarray): Input spectrum.
#        oversample (int): (**Optional**) Oversampling can be done, not
#            to loose spectral resolution, especally for extended
#            wavelength ranges and to avoid aliasing.  Default is to
#            provide the same number of output pixels as input.
#        velscale (float): (**Optional**) Velocity scale in km/s per
#            pixels. If this variable is not defined, then it will
#            contain in output the velocity scale.  If this variable is
#            defined by the user it will be used to set the output number
#            of pixels and wavelength scale.
#        flux (bool): (**Optional**) Set this keyword to preserve total
#            flux.  In this case the log rebinning changes the pixels
#            flux in proportion to their dLam so the following command
#            will show large differences beween the spectral shape before
#            and after :func:`log_rebin`::
#     
#                # Plot log-rebinned spectrum
#                pyplot.plot(exp(logLam), specNew)
#                pyplot.plot(numpy.arange(lamRange[0],lamRange[1],spec.size), spec, 'g')
#                pyplot.show()
#     
#            By default, when this keyword is *not* set, the above two
#            lines produce two spectra that almost perfectly overlap each
#            other.
#        log10 (bool): (**Optional**) Flag that the spectrum should be
#            binned in units of base-10 log wavelength, instead of
#            natural log
#        newRange (numpy.ndarray): (**Optional**) Force the spectrum to
#            be sampled to a new spectral range (lamRange is the
#            *existing* spectral range).
#        wave_in_ang (bool): (**Optional**) Return the wavelength
#            coordinates in angstroms, not log(angstroms)
#        unobs (float): (**Optional**) Default value for unobserved
#            spectral regions.
#
#    Returns:
#        numpy.ndarray, float: Returns three variables: logarithmically
#        rebinned spectrum, the log of the wavelength at the geometric
#        center of each pixel, and the velocity scale of each pixel in
#        km/s.
#        
#    Raises:
#        ValueError: Raised if the input spectrum is not a
#            one-dimensional numpy.ndarray.
#        
#    *Modification History*:
#        | **V1.0.0**: Using interpolation. Michele Cappellari, Leiden,
#            22 October 2001
#        | **V2.0.0**: Analytic flux conservation. MC, Potsdam, 15 June
#            2003
#        | **V2.1.0**: Allow a velocity scale to be specified by the
#            user.  MC, Leiden, 2 August 2003
#        | **V2.2.0**: Output the optional logarithmically spaced
#            wavelength at the geometric mean of the wavelength at the
#            border of each pixel.  Thanks to Jesus Falcon-Barroso. MC,
#            Leiden, 5 November 2003
#        | **V2.2.1**: Verify that lamRange[0] < lamRange[1].  MC,
#            Vicenza, 29 December 2004
#        | **V2.2.2**: Modified the documentation after feedback from
#            James Price.  MC, Oxford, 21 October 2010
#        | **V2.3.0**: By default now preserve the shape of the spectrum,
#            not the total flux. This seems what most users expect from
#            the procedure.  Set the keyword /FLUX to preserve flux like
#            in previous version.  MC, Oxford, 30 November 2011
#        | **V3.0.0**: Translated from IDL into Python. MC, Santiago, 23
#            November 2013
#        | **V3.1.0**: Fully vectorized log_rebin. Typical speed up by
#            two orders of magnitude.  MC, Oxford, 4 March 2014
#        | **05 Jun 2015**: (K. Westfall, KBW) Pulled from ppxf_util.py.
#            Conform to mangadap documentation standard.  Transcribe
#            edits made to IDL version that provides for the log10 and
#            newRange arguments.  Add option to return wavelength in
#            angstroms, not log(angstroms).  Break out determination of
#            input and output spectrum pixel coordinates to a new
#            function, :func:`log_rebin_pix`.  Added default value for
#            unobserved pixels.  Default behavior unchanged.
#
#    .. todo::
#
#        - Allow to resample an already geometrically binned spectrum
#    
#    """
#    lamRange = numpy.asarray(lamRange)
#
#    if type(spec) != numpy.ndarray:
#        raise ValueError('Input spectrum must be a numpy.ndarray')
#    s = spec.shape
#    if len(s) != 1:
#        raise ValueError('input spectrum must be a vector')
#    n = s[0]
#
#    # This is broken out into a separate procedure so that it can be
#    # called to determine the size of the rebinned spectra without
#    # actually doing the rebinning
#    dLam, m, logscale, velscale = \
#        log_rebin_pix(lamRange, n, oversample=oversample, velscale=velscale, log10=log10,
#                      newRange=newRange)
#    print(dLam)                        
#    print(m)
#    print(logscale)
#    print(velscale)
#
#    # Get the sampling of the existing spectrum
#    lim = lamRange/dLam + [-0.5, 0.5]           # All in units of dLam
#    borders = numpy.linspace(*lim, num=n+1)     # Linearly sampled pixels
#
#    print(borders)
#    print(dLam)
#
#    # Set limits to a new wavelength range
#    if newRange is not None:
#        lim = numpy.asarray(newRange)/dLam + [-0.5, 0.5]
#
#    # Set the limits to the (base-10 or natural) log of the wavelength
#    logLim = numpy.log(lim) if not log10 else numpy.log10(lim)
#    logLim[1] = logLim[0] + m*logscale      # Set last wavelength, based on integer # of pixels
#
#    # Geometrically spaced pixel borders for the new spectrum
##    newBorders = numpy.logspace(*logLim, num=m+1, base=(10.0 if log10 else numpy.exp(1)))
#    newBorders = numpy.power(10., numpy.linspace(*logLim, num=m+1)) if log10 else \
#                 numpy.exp(numpy.linspace(*logLim, num=m+1))
#
#
#    print(newBorders)
#    print(m)
#    print(logscale)
#
#    # Get the new spectrum by performing an analytic integral
#    k = (newBorders - borders[0]).clip(0, n-1).astype(int)
#    specNew = numpy.add.reduceat(spec, k)[:-1]
#    specNew *= numpy.diff(k) > 0                # fix for design flaw of reduceat()
#    specNew += numpy.diff((newBorders - borders[k])*spec[k])
#
#    # Don't conserve the flux
#    if not flux:
#        specNew /= numpy.diff(newBorders)
#
#    # Output log(wavelength): log of geometric mean
#    LamNew = numpy.sqrt(newBorders[1:]*newBorders[:-1])*dLam
#
#    # Set values for unobserved regions
#    if newRange is not None and (newRange[0] < lamRange[0] or newRange[1] > lamRange[1]):
#            specNew[ (LamNew < lamRange[0]) | (LamNew > lamRange[1]) ] = unobs
#
#    # Return log(wavelength), if requested
#    if not wave_in_ang:
#        LamNew = numpy.log10(LamNew) if log10 else numpy.log(LamNew)
#
#    # Return spectrum, wavelength coordinates, and pixel size in km/s
#    return specNew, LamNew, velscale
#
#
#
#def log_rebin_pix(lamRange, n, oversample=None, velscale=None, log10=False, newRange=None):
#    """
#    Determine the number of new pixels and their coordinate step when
#    rebinning a spectrum in geometrically stepped bins.  The input
#    spectrum must be sampled linearly in wavelength.  This is primarily
#    a support routine for :func:`log_rebin`.
#
#    Although breaking this out from the main :func:`log_rebin` function
#    leads to a few repeat calculations in that function, the use of this
#    function is in determine a common wavelength range for a large
#    number of spectra before resampling the spectra themselves.  See
#    :class:`mangadap.TemplateLibrary` .
#
#    Args:
#        lamRange (numpy.ndarray): two elements vector containing the
#            central wavelength of the first and last pixels in the
#            spectrum, which is assumed to have constant wavelength
#            scale! E.g. from the values in the standard FITS keywords:
#            LAMRANGE = CRVAL1 + [0,CDELT1*(NAXIS1-1)].  It must be
#            LAMRANGE[0] < LAMRANGE[1].
#        n (int): Number of pixels in the original spectrum.
#        oversample (int): (**Optional**) Oversampling can be done, not
#            to loose spectral resolution, especally for extended
#            wavelength ranges and to avoid aliasing.  Default is to
#            provide the same number of output pixels as input.
#        velscale (float): (**Optional**) Velocity scale in km/s per
#            pixels. If this variable is not defined, then it will
#            contain in output the velocity scale.  If this variable is
#            defined by the user it will be used to set the output number
#            of pixels and wavelength scale.
#        log10 (bool): (**Optional**) Flag that the spectrum should be
#            binned in units of base-10 log wavelength, instead of
#            natural log
#        newRange (numpy.ndarray): (**Optional**) Force the spectrum to
#            be sampled to a new spectral range (lamRange is the
#            *existing* spectral range).
#
#    Returns:
#        float, int: Returns
#            
#            1. the linear wavelength step of each pixel in the input
#            spectrum, 
#            2. the number of pixels for the rebinned spectrum, 
#            3. the log-linear wavelength step for each pixel in the new
#            spectrum, and
#            4. the velocity step for each pixel in the new spectrum.
#
#    Raises:
#        ValueError: Raised if the input wavelength range (*lamRange* or
#            *newRange*) does not have two elements or is not sorted.
#    """
#    lamRange = numpy.asarray(lamRange)
#    if len(lamRange) != 2:
#        raise ValueError('lamRange must contain two elements')
#    if lamRange[0] >= lamRange[1]:
#        raise ValueError('It must be lamRange[0] < lamRange[1]')
#
#    # Size of output spectrum
#    m = int(n) if oversample is None else int(n*oversample)
#
#    # Get the sampling of the existing spectrum
#    dLam = numpy.diff(lamRange)/(n - 1.)        # Assume constant dLam
#    lim = lamRange/dLam + [-0.5, 0.5]           # All in units of dLam
#
#    # Get the sampling for the new spectrum, if requested to be
#    # different
#    if newRange is not None:
#        newRange = numpy.asarray(newRange)
#        if len(newRange) != 2:
#            raise ValueError('newRange must contain two elements')
#        if newRange[0] >= newRange[1]:
#            raise ValueError('It must be newRange[0] < newRange[1]')
#        lim = newRange/dLam + [-0.5, 0.5]       # Set limits to a new wavelength range
#
#        # Adjust the length
#        nn = int((lamRange[1]-lamRange[0]-newRange[1]+newRange[0])/dLam)
#        m = m-nn if oversample is None else m-nn*oversample
#
#    # Set the limits to the (base-10 or natural) log of the wavelength
#    logLim = numpy.log(lim) if not log10 else numpy.log10(lim)
#
#    c = astropy.constants.c.to('km/s').value    # Speed of light in km/s (use astropy definition)
#
#    # Set the velocity scale, if velscale not provided; otherwise force
#    # the sampling based in the input velscale
#    if velscale is None:                        # Velocity scale is not set by user
#        velscale = numpy.diff(logLim)[0]/m*c    # Only for output
#        if log10:
#            velscale *= numpy.log(10.)          # Adjust to log base-10
#
#    logscale = velscale/c                       # dlambda/lambda = dln(lambda)
#    if log10:
#        logscale /= numpy.log(10.)              # Convert to dlog10(lambda)
#    m = int(numpy.diff(logLim)/logscale)        # Number of output pixels
#
#    return dLam, m, logscale, velscale


def _pixel_centers(xlim, npix, log=False, base=10.0):
    """
    Determine the centers of pixels in a linearly or geometrically
    sampled vector.

    Args:
        xlim (numpy.ndarray) : (Geometric) Centers of the first and last
            pixel in the vector.
        npix (int) : Number of pixels in the vector.
        log (bool) : (**Optional**) The input range is (to be)
            logarithmically sampled.
        base (float) : (**Optional**) The base of the logarithmic
            sampling.  The default is 10.0; use numpy.exp(1.) for the
            natural logarithm.

    Returns:
        numpy.ndarray, float: A vector with the (npix+1) borders of the
        pixels and the sampling rate.  If logarithmically binned, the
        sampling is the step in :math`\log x`.
    """
    if log:
        logRange = numpy.log(xlim)/numpy.log(base)
        dlogx = numpy.diff(logRange)/(npix-1.)
        centers = numpy.power(base, numpy.linspace(*(logRange/dlogx), num=npix)*dlogx)
        return centers, dlogx
    dx = numpy.diff(xlim)/(npix-1.)
    centers = numpy.linspace(*(xlim/dx), num=npix)*dx
    return centers, dx


def _pixel_borders(xlim, npix, log=False, base=10.0):
    """
    Determine the borders of the pixels in a vector.

    Args:
        xlim (numpy.ndarray) : (Geometric) Centers of the first and last
            pixel in the vector.
        npix (int) : Number of pixels in the vector.
        log (bool) : (**Optional**) The input range is (to be)
            logarithmically sampled.
        base (float) : (**Optional**) The base of the logarithmic
            sampling.  The default is 10.0; use numpy.exp(1.) for the
            natural logarithm.

    Returns:
        numpy.ndarray, float: A vector with the (npix+1) borders of the
        pixels and the sampling rate.  If logarithmically binned, the
        sampling is the step in :math`\log x`.
    """
    if log:
        logRange = numpy.log(xlim)/numpy.log(base)
        dlogx = numpy.diff(logRange)/(npix-1.)
        borders = numpy.power(base, numpy.linspace(*(logRange/dlogx + [-0.5, 0.5]),
                                                   num=npix+1)*dlogx)
        return borders, dlogx
    dx = numpy.diff(xlim)/(npix-1.)
    borders = numpy.linspace(*(xlim/dx + numpy.array([-0.5, 0.5])), num=npix+1)*dx
    return borders, dx


def resample_vector_npix(outRange=None, dx=None, log=False, base=10.0, default=None):
    """
    Determine the number of pixels needed to resample the vector.

    Args:
        outRange (list or numpy.ndarray) : Two-element array with the
            starting and ending x coordinate of the pixel centers to
            divide into pixels of a given width.  If *log* is True, this
            must still be the linear value of the x coordinate, not
            log(x)!.
        dx (float) : Linear or logarithmic pixel width.
        log (bool) : Flag that the range should be logarithmically
            binned.
        base (float) : Base for the logarithm
        default (int) : Default number of pixels to use.  The default is
            returned if either *outRange* or *dx* are not provided.

    Returns:
        int, numpy.ndarray: Returns two objects: The number of pixels to
        cover *outRange* with pixels of width *dx* and the adjusted
        range such that number of pixels of size dx is the exact integer.

    Raises:
        ValueError: Raised if the range is not a two-element vector
    """
    # If the range or sampling are not provided, the number of pixels is
    # already set
    if outRange is None or dx is None:
        return default, outRange
    if len(outRange) != 2:
        raise ValueError('Output range must be a 2-element vector.')

    _outRange = numpy.atleast_1d(outRange).copy()
    npix = int( numpy.diff(numpy.log(_outRange))/numpy.log(base) / dx) + 1 if log else \
                int(numpy.diff(_outRange)/dx) + 1
#    _outRange = outRange
    _outRange[1] = numpy.power(base, numpy.log(_outRange[0])/numpy.log(base) + dx*(npix-1)) \
                            if log else _outRange[0] + dx*(npix-1)
    return npix, _outRange


class Resample:
    r"""
    Resample regularly or irregularly sampled, 1-dimensional data to a
    new grid using integration.
    
    This is basically a generalization of the routine :func:`log_rebin`
    provided by Michele Cappellari in the pPXF package.

    The abscissa coordinates (`x`) for the data (`y`) should be provided
    for irregularly sampled data.  If the input data is linearly or
    geometrically sampled (`inLog=True`), the abscissa coordinates can
    can be generated using the input range for the (geometric) center of
    each grid point.  If neither `x` nor `xRange` are provided, the
    function assumes grid coordinates of `x=numpy.arange(y.size)`.

    The function resamples the data by constructing the borders of the
    output grid using the `new*` keywords and integrating the input
    function between those borders.  The output data will be set to
    `ext_value` for any data beyond the abscissa limits of the input
    data.

    The `conserve` keyword sets how the units of the input data should
    be treated.  If `conserve=False`, the input data are expected to be
    in units such that the integral over :math:`dx` is independent of
    the units of :math:`x` (i.e., flux per unit angstrom, or flux
    density).  If `conserve=True`, the value of the data is assumed to
    have been integrated over the size of each pixel (i.e., units of
    flux).  If `conserve=True`, :math:`y` is converted to units of per
    step in :math:`x` such that the integral before and after the
    resample is the same.  For example, if :math:`y` is a spectrum in
    units of flux, the function first converts the units to flux density
    and then computes the integral over each new pixel to produce the
    new spectra with units of flux.

    Args:
        y (numpy.ndarray):
            Data values to resample.  Shape can be 1 or 2D.  If 1D, the
            shape must be :math:`(N_{\rm pix},)`; otherwise, it must be
            :math:`(N_y,N_{\rm pix})`.  I.e., the length of the last
            axis must match the input coordinates.
        x (numpy.ndarray, optional):
            Abcissa coordinates for the data, which do not need to be
            regularly sampled.  If the pixel borders are not provided,
            they are assumed to be half-way between adjacent pixels, and
            the first and last borders are assumed to be equidistant
            about the provided value.  If these coordinates are not
            provided, they are determined by the input borders, the
            input range, or just assumed to be the indices,
            :math:`0..N_{\rm pix}-1`.
        xRange (array-like, optional):
            A two-element array with the starting and ending value for
            the coordinates of the centers of the first and last pixels
            in y.  Default is :math:`[0,N_{\rm pix}-1]`.
        xBorders (numpy.ndarray, optional):
            An array with the borders of each pixel that must have a
            length of :math:`N_{\rm pix}+1`.
        inLog (:obj:`bool`, optional):
            Flag that the input is logarithmically binned, primarily
            meaning that the coordinates are at the geometric center of
            each pixel and the centers are spaced logarithmically.  If
            false, the sampling is expected to be linear.
        newRange (array-like, optional):
            A two-element array with the (geometric) centers of the
            first and last pixel in the output vector.  If not provided,
            assumed to be the same as the input range.
        newpix (:obj:`int`, optional): 
            Number of pixels for the output vector.  If not provided,
            assumed to be the same as the input vector.
        newLog (:obj:`bool`, optional):
            The output vector should be logarithmically binned.
        newdx (:obj:`float`, optional):
            The sampling step for the output vector.  If `newLog=True`,
            this has to be the change in the logarithm of x for the
            output vector!  If not provided, the sampling is set by the
            output range (see `newRange` above) and number of pixels
            (see `newpix` above).
        base (:obj:`float`, optional):
            The base of the logarithm used for both input and output
            sampling, if specified.  The default is 10; use
            `numpy.exp(1)` for natural logarithm.
        ext_value (:obj:`float`, optional):
            Set extrapolated values to the provided float.  By default,
            extrapolated values are set to 0.  If set to None, values
            are just set to the linear exatrapolation of the data beyond
            the provided limits; use `ext_value=None` with caution!
        conserve (:obj:`bool`, optional):
            Conserve the integral of the input vector.  For example, if
            the input vector is a spectrum in flux units, you should
            conserve the flux in the resampling; if the spectrum is in
            units of flux density, you do not want to conserve the
            integral.
    
    Attributes:
        oy...
    
    Raises:
        ValueError: Raised if *y* is not of type numpy.ndarray, if *y*
            is not one-dimensional, or if *xRange* is not provided and
            the input vector is logarithmically binned (see *inLog*
            above).
    """
    def __init__(self, y, x=None, xRange=None, xBorders=None, inLog=False, newRange=None,
                 newpix=None, newLog=True, newdx=None, base=10.0, ext_value=0.0, conserve=False):

        # Check operation can be performed
        if not isinstance(y, numpy.ndarray):
            raise ValueError('Input vector must be a numpy.ndarray!')
        if len(y.shape) > 2:
            raise ValueError('Input must be a 1D or 2D array!')

        self.y = y.copy()
        self.x = None
        self.xborders = None
        self._input_coordinates(x, xRange, xBorders, inLog, base)

        # If conserving integral, assume input is integrated over pixel
        # width and convert to a density function to be integrated
        if conserve:
            self.y /= numpy.diff(inBorders)

        self._output_coordinates(newRange, newpix, newLog, newdx, base)

        self.resample(ext_value=ext_value, conserve=conserve)

    def _input_coordinates(self, x, xRange, xBorders, inLog, base):
        """
        Determine the centers and pixel borders of the input
        coordinates.
        """
        if (x is not None or xBorders is not None) and xRange is not None:
            warnings.warn('Provided both x or x borders and the x range.  Ignoring range.')
        _xRange = xRange if x is None and xBorders is None else None
        
        if x is not None:
            if x.ndim != 1:
                raise ValueError('Coordinate vector must be 1D.')
            if x.size != self.y.shape[-1]:
                raise ValueError('Coordinate vector must match last dimension of value array.')
        if xBorders is not None:
            if xBorders.ndim != 1:
                raise ValueError('Coordinate borders must be 1D.')
            if xBorders.size != self.y.shape[-1]+1:
                raise ValueError('Coordinate borders must match last dimension of value array.')

        if x is None:
            if xBorders is not None:
                self.x = numpy.sqrt(xBorders[:-1]*xBorders[1:]) if inLog \
                            else (xBorders[:-1]+xBorders[1:])/2.0
            elif xRange is not None:
                self.x = _pixel_centers(xRange, self.y.shape[-1], log=inLog, base=base)
            else:
                self.x = numpy.arange(self.y.shape[-1]) + 0.5
        else:
            self.x = x

        if xBorders is None:
            dx = numpy.diff(numpy.log(self.x)) if inLog else numpy.diff(self.x)
            self.xborders = numpy.exp(numpy.append(numpy.log(self.x[:-1]) - dx/2,
                                        numpy.log(self.x[-1]) + numpy.array([-1,1])*dx[-1]/2)) \
                                if inLog else numpy.append(self.x - dx/2, [self.x[-1] + dx[-1]/2])
        else:
            self.xborders = xBorders

    def _output_coordinates(self, newRange, newpix, newLog, newdx, base):
        """Set the output coordinates."""

        # Set the output range and number of pixels
        outRange = numpy.array([self.x[0], self.x[-1]]) if newRange is None \
                        else numpy.array(newRange)
        m, _outRange = resample_vector_npix(outRange=outRange, log=newLog, base=base, dx=newdx,
                                        default=(self.y.shape[-1] if newpix is None else newpix))
        outRange = outRange if _outRange is None else _outRange

        # Get the output pixel borders
        self.outborders = _pixel_borders(outRange, m, log=newLog, base=base)[0]

        # Get the output coordinate vector
        self.outx = numpy.sqrt(self.outborders[:-1]*self.outborders[1:]) if newLog \
                        else (self.outborders[:-1]+self.outborders[1:])/2.0

    def old_resample(self, ext_value=0.0, conserve=False):
        """Resample the vectors."""

        # Combine the input coordinates and the output borders into a
        # single vector
        combinedX = numpy.append(self.outborders, self.x)
        srt = numpy.argsort(combinedX)
        combinedX = combinedX[srt]
        border = numpy.ones(combinedX.size, dtype=bool)
        border[self.outborders.size:] = False

        # Linearly interpolate the input function at the output border positions
        interp = interpolate.interp1d(self.x, self.y, assume_sorted=True, fill_value='extrapolate')
        combinedY = numpy.append(interp(self.outborders), self.y)[srt]

        # Use reduceat to calculate the integral
        integrand = (combinedY[1:]+combinedY[:-1])*numpy.diff(combinedX)/2.0
        k = numpy.arange(combinedX.size)[border[srt]]
        self.outy = numpy.add.reduceat(integrand, k[:-1]) if k[-1] == combinedY.size-1 \
                            else numpy.add.reduceat(integrand, k)[:-1]

        # Do not conserve the integral over the size of the pixel
        if not conserve:
            self.outy /= numpy.diff(self.outborders)

        # Set values for extrapolated regions
        if ext_value is not None:
            indx = (self.outborders[:-1] < self.xborders[0]) \
                        | (self.outborders[1:] > self.xborders[-1]) 
            if numpy.sum(indx) > 0:
                self.outy[indx] = ext_value

        return self.outx, self.outy


    def resample(self, ext_value=0.0, conserve=False):
        """Resample the vectors."""

        # Convert y to a step function
        _y = numpy.repeat(self.y, 2)
        _x = numpy.repeat(self.xborders, 2)[1:-1]

        # Combine the input coordinates and the output borders into a
        # single vector
        indx = numpy.searchsorted(_x, self.outborders)
        combinedX = numpy.insert(_x, indx, self.outborders)

        # Flag the output borders
        border = numpy.insert(numpy.zeros(_x.size, dtype=bool), indx,
                                numpy.ones(self.outborders.size, dtype=bool))

        # Insert points at the borders of the output function
        interp = interpolate.interp1d(_x, _y, assume_sorted=True, fill_value='extrapolate')
        combinedY = numpy.insert(_y, indx, interp(self.outborders))

        # Use reduceat to calculate the integral
        integrand = combinedY[1:]*numpy.diff(combinedX)
        k = numpy.arange(combinedX.size)[border]
        self.outy = numpy.add.reduceat(integrand, k[:-1]) if k[-1] == combinedY.size-1 \
                            else numpy.add.reduceat(integrand, k)[:-1]

        # Do not conserve the integral over the size of the pixel
        if not conserve:
            self.outy /= numpy.diff(self.outborders)

        # Set values for extrapolated regions
        if ext_value is not None:
            indx = (self.outborders[:-1] < self.xborders[0]) \
                        | (self.outborders[1:] > self.xborders[-1]) 
            if numpy.sum(indx) > 0:
                self.outy[indx] = ext_value

