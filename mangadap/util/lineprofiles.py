# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Implements a set of line profile parameterizations.

----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import numpy
from scipy import special, fftpack
from astropy.modeling import FittableModel, Parameter

#-----------------------------------------------------------------------
class GaussianLineProfile(FittableModel):
    r"""
    Define a Gaussian line profile as parameterized by its zeroth
    moment, mean, and standard deviation:

    .. math::

        \mathcal{N}(x|f,\mu,\sigma) = \frac{f}{\sqrt{2\pi}\sigma}
        \exp\left(\frac{-(x-\mu)^2}{2\sigma^2}\right)

    The base class is `astropy.modeling.FittableModel`_, which
    facilitates its use in combining multiple components and other
    models in the `astropy.modeling`_ suite.
    """
    # Inputs the x positions
    inputs = ('x',)
    # Returns the renormalized Gaussian PDF at each x position
    outputs = ('y',)

    # Parameters are the peak value of the profile, the profile mean,
    # and the profile standard deviation (sigma)
    zmom=Parameter(default=1.0)
    mean=Parameter(default=0.0)
    sigma=Parameter(default=1.0)

    
    @staticmethod
    def evaluate(x, zmom, mean, sigma):
        # Discrete samples
        return zmom * numpy.exp(-0.5*numpy.square((x - mean)/sigma)) / numpy.sqrt(2*numpy.pi)/sigma
#        # Integrated over a pixel
#        diff = (x-mean)/pixelscale
#        denom = numpy.sqrt(2.0)*sigma/pixelscale
#        return zmom * (erf((diff+0.5)/denom) - erf((diff-0.5)/denom))/2.


    @staticmethod
    def fit_deriv(x, zmom, mean, sigma):
        y = GaussianLineProfile.evaluate(x,zmom,mean,sigma)
        d_zmom = y/zmom
        d_mean = y*(x-mean)/numpy.square(sigma)
        d_sigma = y*(numpy.square(x-mean)/numpy.power(sigma,3) - 1.0/sigma)
        return [d_zmom, d_mean, d_sigma]


    @staticmethod
    def flux(zmom, mean, sigma):
        return zmom


    @staticmethod
    def flux_err(zmom, mean, sigma, zmome, meane, sigmae):
        return zmome


    @staticmethod
    def moment(order, zmom, mean, sigma):
        if order == 0:
            return zmom
        if order == 1:
            return mean
        if order == 2:
            return sigma


    @staticmethod
    def moment_err(order, zmom, mean, sigma, zmome, meane, sigmae):
        if order == 0:
            return zmome
        if order == 1:
            return meane
        if order == 2:
            return sigmae


    @staticmethod
    def integral(zmom, mean, sigma):
        return integrate.quad( GaussianLineProfile.evaluate, -sigma, sigma, args=(zmom,mean,sigma))

    
    @staticmethod
    def scale_flux(zmom, mean, sigma, fac):
        return [zmom*fac, mean, sigma]


    @staticmethod
    def fix_flux():
        return [ True, False, False ]

    
    @staticmethod
    def mean_indx():
        return 1

    
    @staticmethod
    def shift_mean(zmom, mean, sigma, shift):
        return [zmom, mean+shift, sigma ]

    
    @staticmethod
    def fix_mean():
        return [ False, True, False ]

    
    @staticmethod
    def scale_stddev(zmom, mean, sigma, fac):
        return [zmom, mean, sigma*fac ]


    @staticmethod
    def fix_stddev():
        return [ False, False, True ]


class GaussianLSF:
    r"""
    Define a Gaussian line profile, *sampled* over the width of the
    sampling step, parameterized by its integral (:math:`F`), center
    (:math:`\mu`), and standard deviation (:math:`\sigma`). I.e:

    .. math::

        \mathcal{N}(x|f,\mu,\sigma) = \frac{f}{\sqrt{2\pi}\sigma}
        \exp\left(\frac{-\Delta^2}{2\sigma^2}\right)

    where :math:`\Delta = x-\mu`.  The coordinate vector :math:`x` does
    not need to be uniformly sampled.

    Args:
        p (array-like): (**Optional**) Input parameters ordered as the
            total integral of the profile, the profile center, and the
            profile standard deviation.  Assumed to be (1.0, 0.0, 1.0)
            by default.

    Attributes:
        p (numpy.ndarray): Most recently used parameters

    Raises:
        ValueError: Raised if the provided parameter vector is not 3
            elements long.
    """
    def __init__(self, p=None):
        self.set_par(p)


    def __call__(self, x, p):
        """Calculate the profile.

        Args:
            x (array-like): Independent variable.
            p (array-like): LSF parameters.
        """
        self.set_par(p)
        return self.sample(x)


    @staticmethod
    def npar():
        return 3


    def set_par(self, p):
        """
        Set the internal parameters to the provided set.

        Args:
            p (array-like): LSF parameters.

        Raises:
            ValueError: Raised if the provided parameter vector is not 3
                elements long.
        """
        if p is None:
            self.p = numpy.array([1.0, 0.0, 1.0])
            return
        if len(p) != GaussianLSF.npar():
            raise ValueError('Must provide {0} parameters.'.format(GaussianLSF.npar()))
        self.p = numpy.asarray(p)


    def sample(self, x):
        """
        Sample the profile.

        Args:
            x (array-like): Independent variable.
        """
        return self.p[0] * numpy.exp(-numpy.square((x-self.p[1])/self.p[2])/2.) \
                    / numpy.sqrt(2.0*numpy.pi) / self.p[2]


    def parameters_from_moments(self, mom0, mom1, mom2):
        """
        Provided the 0th, 1st, and 2nd moments, produce a set of
        parameters for the profile.
        """
        return numpy.array([mom0, mom1, mom2])


class IntegratedGaussianLSF(GaussianLSF):
    r"""
    Define a Gaussian line profile, *integrated* over the width of the
    sampling step, parameterized by its integral (:math:`F`), center
    (:math:`\mu`), and standard deviation (:math:`\sigma`). I.e:

    .. math::

        \mathcal{N}(x|F,\mu,\sigma) = \frac{F}{2} \left[
        {\rm erf}\left(\frac{\Delta+\delta_x/2}{\sqrt{2}\sigma}\right) - 
        {\rm erf}\left(\frac{\Delta-\delta_x/2}{\sqrt{2}\sigma}\right)\right]

    where :math:`{\rm erf}(x)` is the error function, :math:`\Delta =
    x-\mu`, and :math:`\delta_x` is the sampling step.  The sampling
    *must* be uniform in :math:`x`.

    Args:
        p (array-like): (**Optional**) Input parameters ordered as the
            total integral of the profile, the profile center, and the
            profile standard deviation.  Assumed to be (1.0, 0.0, 1.0)
            by default.
        dx (float): (**Optional**) Sampling width.  Default is 1.

    Attributes:
        p (numpy.ndarray): Most recently used parameters
        dx (float): Assumed sampling.

    Raises:
        ValueError: Raised if the provided parameter vector is not 3
            elements long.
    """
    def __init__(self, p=None, dx=None):
        self.set_par(p)
        self.dx = 1.0 if dx is None else dx


    def sample(self, x):
        """
        Sample the profile.

        .. warning::
            Does **not** check if the provided :math:`x` values are
            sampled at :attr:`dx`.

        Args:
            x (array-like): Independent variable.
        """
        n = numpy.sqrt(2.)*self.p[2]
        d = numpy.asarray(x)-self.p[1]
        return self.p[0] * (special.erf((d+self.dx/2.)/n) - special.erf((d-self.dx/2.)/n))/2.


   
class FFTGaussianLSF(GaussianLSF):
    r"""

    Define a Gaussian line profile by first constructing the analytic
    FFT of the profile and then returning the inverse real FFT.  See
    ppxf_util.emline by M. Cappellari.  The sampling *must* be uniform
    in :math:`x`.
    
    Args:
        p (array-like): (**Optional**) Input parameters ordered as the
            total integral of the profile, the profile center, and the
            profile standard deviation.  Assumed to be (1.0, 0.0, 1.0)
            by default.
        dx (float): (**Optional**) Sampling width.  Default is 1.
        pixel (bool): (**Optional**) Flag to produce profile integrated
            over the sampling width.
        
    Attributes:
        p (numpy.ndarray): Most recently used parameters
        dx (float): Assumed sampling.
        pixel (bool): Flag to produce profile integrated over the
            sampling width.

    Raises:
        ValueError: Raised if the provided parameter vector is not 3
            elements long.
    """
    def __init__(self, p=None, dx=None, pixel=True):
        self.set_par(p)
        self.dx = 1.0 if dx is None else dx
        self.pixel = pixel

    def sample(self, x):
        """
        Sample the profile.

        .. warning::
            Does **not** check if the provided :math:`x` values are
            sampled at :attr:`dx`.

        Args:
            x (array-like): Independent variable.
        """
        # Require that the center of the line be within the range of x
        if self.p[1] < x[0] or self.p[1] > x[-1]:
            return numpy.zeros_like(x)
        xsig = self.p[2]/self.dx
        x0 = (self.p[1]-x[0])/self.dx
        npad = 2**int(numpy.ceil(numpy.log2(x.size)))
#        npad = fftpack.next_fast_len(x.size)
        w = numpy.linspace(0,numpy.pi,npad//2+1)
        rfft = self.p[0]*numpy.exp(-0.5*numpy.square(w*xsig) - 1j*w*x0)
        if self.pixel:
            rfft *= numpy.sinc(w/(2*numpy.pi))
        lsf = numpy.fft.irfft(rfft, n=npad)[:x.size]
        return lsf if self.pixel else lsf/self.dx


#-----------------------------------------------------------------------
class NCompLineProfile:
    """
    Construct a single line profile from many components with the same
    profile parameterization.
    """
    def __init__(self, ncomp, par=None, err=None, profile=GaussianLineProfile):
        # Check the type of the profile to use
        if not issubclass(profile, FittableModel):
            raise TypeError('Profile must be a subclass of FittableModel.')
        self.profile_class = profile

        # Ensure that the profile can be defined
        if not ncomp > 0:
            raise ValueError('Number of input components must be 1 or higher!')

        # Calculate and save the number of parameters
        self.ncomp = ncomp
        self.npar_per_component = len(self.profile_class.param_names)
        self.npar = self.ncomp*self.npar_per_component

        # Check the parameter vector, if provided
        if par is not None and len(par.shape) > 2:
            raise ValueError('Input parameter array should only be one or two dimensional.')
        _par = None if par is None else par.ravel()
        if _par is not None and len(_par) != self.npar:
            raise ValueError('Incorrect number of parameters provided')
        self.par = None if _par is None else _par.reshape(self.ncomp,-1).astype(float)

        # Check the parameter vector, if provided
        if err is not None and len(err.shape) > 2:
            raise ValueError('Input parameter error array should only be one or two dimensional.')
        _err = None if err is None else err.ravel()
        if _err is not None and len(_err) != self.npar:
            raise ValueError('Incorrect number of parameter errors provided')
        self.err = None if _err is None else _err.reshape(self.ncomp,-1).astype(float)

        # Build the profile function
        self.profile = self.profile_class() \
                if self.par is None else self.profile_class(*self.par[0])
        for i in range(1,ncomp):
            self.profile += self.profile_class() \
                    if self.par is None else self.profile_class(*self.par[i])


    def __call__(self, x, par=None):
        return self.sample(x,par)


    def _quick_sample(self, x):
        """
        Sample without providing/checking new input parameters.
        """
        return self.profile.evaluate(x, *self.par.ravel())


    def _mom1_integrand(self, x):
        return x*self._quick_sample(x)


    def _mom2_integrand(self, x):
        return x*x*self._quick_sample(x)
    

    def assign_par(self, par):
        _par = par.ravel() if len(par.shape) == 2 else par
        if len(_par) != self.npar:
            raise ValueError('Incorrect number of parameters provided')
        self.par = _par.copy().reshape(self.ncomp,-1)


    def assign_err(self, err):
        _err = err.ravel() if len(err.shape) == 2 else err
        if len(_err) != self.npar:
            raise ValueError('Incorrect number of parameters provided')
        self.err = _err.copy().reshape(self.ncomp,-1)


    def sample(self, x, par=None):
        if par is not None:
            self.assign_par(par)
        return self._quick_sample(x)


    def flux(self, par=None):
        if par is not None:
            self.assign_par(par)
        return numpy.sum( [ self.profile_class.flux(*p) \
                            for p in numpy.take(self.par,numpy.arange(self.ncomp),axis=0) ] )


    def flux_err(self, par=None, err=None):
        if par is not None:
            self.assign_par(par)
        if err is not None:
            self.assign_err(err)
        return numpy.sum([ self.profile_class.flux_err(*p, *e) \
                           for p,e in zip(numpy.take(self.par,numpy.arange(self.ncomp),axis=0),
                                          numpy.take(self.err,numpy.arange(self.ncomp),axis=0)) ] )


    def moment(self, order=0, par=None):
        """
        .. todo::
            impose some reasonable limits for the integrals
        """
        if par is not None:
            self.assign_par(par)

        if self.ncomp == 1:
            return self.profile_class.moment(order, *self.par[0,:])

        flux = self.flux(par=par)
        if order == 0:
            return flux

        mom1, _ = integrate.quad(self._mom1_integrand, -numpy.inf, numpy.inf)
        if order == 1:
            return mom1/flux
        mom2, _ = integrate.quad(self._mom2_integrand, -numpy.inf, numpy.inf)
        return numpy.sqrt(mom2/flux-numpy.square(mom1/flux))

    
    def moment_err(self, order=0, par=None, err=None):
        """
        .. todo::
            impose some reasonable limits for the integrals
        """
        if self.err is None and err is None:
            return 0.0

        if par is not None:
            self.assign_par(par)
        if err is not None:
            self.assign_err(err)

        if self.ncomp == 1:
            return self.profile_class.moment_err(order, *self.par[0,:], *self.err[0,:])
        raise NotImplementedError('Unable to calculate error for multiple components.')

    
    def scale_flux(self, fac):
        for i in range(self.ncomp):
            self.par[i,:] = self.profile_class.scale_flux(*self.par[i,:], fac)


    def set_flux(self, flx):
        self.scale_flux( flx/self.flux() )


    def fix_flux(self):
        return numpy.array(self.profile_class.fix_flux()*self.ncomp)


    def mean_indx(self):
        return [ self.profile_class.mean_indx()+i*self.ncomp for i in range(self.ncomp) ]
            

    def shift_mean(self, shift):
        for i in range(self.ncomp):
            self.par[i,:] = self.profile_class.shift_mean(*self.par[i,:], shift)

    
    def set_mean(self, mean):
        self.shift_mean( mean-self.moment(order=1) )

    
    def fix_mean(self):
        return numpy.array(self.profile_class.fix_mean()*self.ncomp)


    def scale_stddev(self, fac):
        if self.ncomp != 1:
            raise NotImplementedError('Cannot scale the standard deviation of multi-component profiles.')
        self.par[0,:] = self.profile_class.scale_stddev(*self.par[0,:], fac)


    def set_stddev(self, stddev):
        self.scale_stddev( stddev/self.moment(order=2) )


    def fix_stddev(self):
        if self.ncomp != 1:
            raise NotImplementedError('Cannot fix the standard deviation of multi-component profiles.')
        return numpy.array(self.profile_class.fix_stddev()*self.ncomp)


