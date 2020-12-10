# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Implements a wrapper class for pPXF.

.. todo::

    Allow new iteration mode that iteratively fits until the velocity is
    not up against one of the +/- 2000 km/s limits?  Could be useful for
    poor redshift guesses.

----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import time
import warnings
import logging

import numpy
from scipy import interpolate, fftpack
import astropy.constants

from ppxf import ppxf, capfit

from ..par.parset import KeywordParSet
from ..util.pixelmask import PixelMask, SpectralPixelMask
from ..util.filter import BoxcarFilter
from ..util.log import log_output
from ..util.sampling import spectrum_velocity_scale, angstroms_per_pixel, Resample
from ..util.resolution import match_spectral_resolution, SpectralResolution
from ..util.constants import DAPConstants
from .spatiallybinnedspectra import SpatiallyBinnedSpectra
from .templatelibrary import TemplateLibrary
from .spectralfitting import StellarKinematicsFit
from .util import sample_growth, optimal_scale
#from .util import residual_growth, optimal_scale

# For debugging
from matplotlib import pyplot
#from ..contrib import ppxf_util

class PPXFFitPar(KeywordParSet):
    r"""
    Define a parameter set used by the pPXF fitting method.

    .. todo::
        The overlap between this and
        :class:`mangadap.proc.stellarcontinuummodel.StellarContinuumModelDef`
        is not well designed.

    The defined parameters are:

    .. include:: ../tables/ppxffitpar.rst

    """
    def __init__(self, template_library_key=None, template_library=None, guess_redshift=None,
                 guess_dispersion=None, iteration_mode=None, reject_boxcar=None,
                 filter_boxcar=None, filter_operation=None, filter_iterations=None,
                 match_resolution=None, velscale_ratio=None, minimum_snr=None, pixelmask=None,
                 bias=None, degree=None, mdegree=None, filt_degree=None, filt_mdegree=None,
                 moments=None):
    
        arr_in_fl = [ numpy.ndarray, list, int, float ] # guess kinematics
        in_fl = [ int, float ]                          # bias, minimum S/N

        _def = self._keyword_defaults()
        
        iter_opt = PPXFFit.iteration_modes()
        moment_opt = [ 2, 4, 6 ]
        filter_operation_opt = [ 'divide', 'subtract' ]

        pars =     [ 'template_library_key', 'template_library', 'guess_redshift',
                     'guess_dispersion', 'iteration_mode', 'reject_boxcar', 'filter_boxcar',
                     'filter_operation', 'filter_iterations', 'match_resolution', 'velscale_ratio',
                     'minimum_snr', 'pixelmask', 'bias', 'degree', 'mdegree', 'filt_degree',
                     'filt_mdegree', 'moments' ]
        values =   [ template_library_key, template_library, guess_redshift, guess_dispersion,
                     iteration_mode, reject_boxcar, filter_boxcar, filter_operation,
                     filter_iterations, match_resolution, velscale_ratio, minimum_snr, pixelmask,
                     bias, degree, mdegree, filt_degree, filt_mdegree, moments ]
        options =  [ None, None, None, None, iter_opt, None, None, filter_operation_opt, None,
                     None, None, None, None, None, None, None, None, None, moment_opt ]
        defaults = [ None, None, None, None, 'global_template', None, None, 'divide', 0, True,
                     None, None, None, _def['bias'], _def['degree'], _def['mdegree'],
                     _def['filt_degree'], _def['filt_mdegree'], _def['moments'] ]
        dtypes =   [ str, TemplateLibrary, arr_in_fl, arr_in_fl, str, int, int, str, int, bool,
                     int, in_fl, PixelMask, in_fl, int, int, int, int, int ]


        descr = ['Keyword of the library to fit.  See ' \
                    ':func:`mangadap.proc.templatelibrary.available_template_libraries`.',
                 'Object with the spectra in the template library that have been prepared ' \
                    'for analysis of the data.',
                 'Initial guess for the redshift (:math:`cz`) of each binned spectrum.',
                 'Initial guess for the velocity dispersion for each binned spectrum.',
                 'Iteration mode to use; see :func:`PPXFFit.iteration_modes`.',
                 'Number of pixels in the boxcar used to determine the local sigma for ' \
                    'rejecting outliers.',
                 'Size of the boxcar in pixels used in a high-pass filter applied before ' \
                    'fitting the spectra. (**To be deprecated**)',
                 'Operation to use when constructing the filtered spectra.  The boxcar smoothed ' \
                    'version of the spectrum is either subtracted or divided into the original ' \
                    'spectrum if ``filter_operation`` is ``subtract`` or ``divide``, '
                    'respectively. (**To be deprecated**)',
                 'Number of fit-reject-filter iterations. (**To be deprecated**)',
                 'Match the spectral resolution of the template to that of the galaxy data.  ' \
                    'This is used only when constructing the template library.  Default is True.',
                 'The **integer** ratio between the velocity scale of the pixel in the galaxy ' \
                    'data to that of the template data.  This is used only when constructing ' \
                    'the template library.  Default is None, which is the same as assuming ' \
                    'that the velocity scales are identical.',
                 'Minimum S/N ratio to include in the fitting.',
                 'Pixel mask to include during the fitting.',
                 '`ppxf`_ ``bias`` parameter used to penalize low S/N spectra toward a ' \
                    'Gaussian LOSVD.',
                 '`ppxf`_ ``degree`` parameter used to set the order of the additive polynomial ' \
                    'to include in the fit.',
                 '`ppxf`_ ``mdegree`` parameter used to set the order of the multiplicative ' \
                    'polynomial to include in the fit.',
                 'Order of the additive polynomial to include when fitting high-pass filtered ' \
                    'spectra.  (**To be deprecated**)',
                 'Order of the multiplicative polynomial to include when fitting high-pass ' \
                    'filtered spectra.  (**To be deprecated**)',
                 r'`ppxf`_ ``moments`` parameter used to set the number of moments of the ' \
                    r'LOSVD to fit.  The DAP has not been well tested for fits that include ' \
                    r'any more than :math:`V` and :math:`\sigma`.']

        super(PPXFFitPar, self).__init__(pars, values=values, defaults=defaults, options=options,
                                         dtypes=dtypes, descr=descr)
        self._check()


    @staticmethod
    def _keyword_defaults():
        """
        Return the keyword defaults.  Pulled from
        :class:`mangadap.contrib.ppxf.ppxf`.
        """
#        return { 'bias':None, 'clean':False, 'degree':4, 'mdegree':0, 'moments':2, 'regul':0,
#                 'reddening':None, 'component':0, 'reg_dim':None }
        return { 'bias':None, 'degree':8, 'mdegree':0, 'filt_degree':8, 'filt_mdegree':0,
                 'moments':2 }


    def _check(self):
        """
        Perform some preliminary checks on the values of the parameters.
        """
        if self['filter_iterations'] > 0 and self['iteration_mode'] != 'fit_reject_filter':
            warnings.warn('Only the \'fit_reject_filter\' iteration mode includes the filtering '
                          'iterations.')
#        if self['filter_iterations'] > 0 and not self['mdegree'] > 0:
#            raise ValueError('If filtering, the multiplicative polynomial must have a non-zero'
#                             ' order (mdegree > 0)')
        
#        if self['reddening'] is not None:
#            if self['mdegree'] > 0:
#                warnings.warn('Cannot both fit multiplicative polynomial and reddening.' \
#                              'Ignoring mdegree.')
#                self['mdegree'] = 0
        # Other checks (and the one above) done within pPXF


    def toheader(self, hdr):
        hdr['PPXFTPLK'] = (self['template_library_key'], 'Template library key used with pPXF')
        hdr['PPXFMODE'] = (self['iteration_mode'], 'pPXF iteration mode')
        hdr['PPXFBIAS'] = (str(self['bias']) if self['bias'] is None else self['bias'],
                            'pPXF bias value')
        hdr['PPXFAO'] = (self['degree'], 'Additive order in pPXF')
        hdr['PPXFMO'] = (self['mdegree'], 'Multiplicative order in pPXF')
        hdr['PPXFFAO'] = (self['filt_degree'], 'Additive order for filtered spectra')
        hdr['PPXFFMO'] = (self['filt_mdegree'], 'Multiplicative order for filtered spectra')
        hdr['PPXFMOM'] = (self['moments'], 'Number of fitted LOSVD moments in pPXF')
        if self['reject_boxcar'] is not None:
            hdr['PPXFRBOX'] = (self['reject_boxcar'], 'pPXF rejection boxcar')
        if self['filter_boxcar'] is not None:
            hdr['PPXFFBOX'] = (self['filter_boxcar'], 'pPXF filtering boxcar')
        if self['filter_iterations'] is not None:
            hdr['PPXFFILT'] = (self['filter_iterations'], 'pPXF number of filtering iterations')
        return hdr


    def fromheader(self, hdr):
        self['template_library_key'] = hdr['PPXFTPLK']
        self['bias'] = eval(hdr['PPXFBIAS'])
        self['degree'] = hdr['PPXFAO']
        self['mdegree'] = hdr['PPXFMO']
        self['filt_degree'] = hdr['PPXFFAO']
        self['filt_mdegree'] = hdr['PPXFFMO']
        self['moments'] = hdr['PPXFMOM']

        try:
            self['reject_boxcar'] = hdr['PPXFRBOX']
        except KeyError as e:
            warnings.warn('Input header does not specify rejection boxcar.')
            self['reject_boxcar'] = None

        try:
            self['filter_boxcar'] = hdr['PPXFFBOX']
        except KeyError as e:
            warnings.warn('Input header does not specify filtering boxcar.')
            self['filter_boxcar'] = None

        try:
            self['filter_iterations'] = hdr['PPXFFILT']
        except KeyError as e:
            warnings.warn('Input header does not specify number of filtering iterations.')
            self['filter_iterations'] = None
        


class PPXFFitResult:
    """
    A basic utility to save the critical parts of the pPXF model.
    """
    def __init__(self, degree, mdegree, start, end, tpl_to_use, ppxf_fit, ntpl,
                 weight_errors=False, component_fits=False):
        self.start = start
        self.end = end
        self.npixtot = end-start
        self.tpl_to_use = tpl_to_use.copy()
        self.gpm = None if ppxf_fit is None else ppxf_fit.goodpixels.copy()
        self.bestfit = None if ppxf_fit is None else ppxf_fit.bestfit.copy()
        self.tplwgt = None if ppxf_fit is None else ppxf_fit.weights[0:ntpl].copy()
#        if numpy.sum(self.tplwgt) == 0:
#            print(ppxf_fit.weights[0:ntpl])
#            print(self.tplwgt)
#            print(self.status)
#            print(ppxf_fit.sol)
#            print(ppxf_fit.error)
        self.tplwgterr = None

        # Set status
        if ppxf_fit is None:
            self.status = None 
        elif numpy.sum(self.tplwgt) == 0:
            warnings.warn('No non-zero templates!  Setting status as failed.')
            self.status = -2            # All template weights are zero!
        else:
            self.status = ppxf_fit.status
#        self.status = None if ppxf_fit is None else ppxf_fit.status

        if ppxf_fit is not None and weight_errors:
            design_matrix = ppxf_fit.matrix[ppxf_fit.goodpixels,:] \
                                / ppxf_fit.noise[ppxf_fit.goodpixels, None]
            cov, _tplwgterr = capfit.cov_err(design_matrix)
            self.tplwgterr = _tplwgterr[degree+1:]
#            covariance_matrix = numpy.linalg.inv(numpy.matmul(design_matrix.T, design_matrix))
#            self.tplwgterr = numpy.sqrt(numpy.diag(covariance_matrix))[degree+1:]
        self.addcoef = None if ppxf_fit is None or degree < 0 else ppxf_fit.polyweights.copy()
        self.multcoef = None if ppxf_fit is None or mdegree <= 0 else ppxf_fit.mpolyweights.copy()
        self.kin = None if ppxf_fit is None else ppxf_fit.sol.copy()
        self.kinerr = None if ppxf_fit is None else ppxf_fit.error.copy()
#        self.robust_rchi2 = None if ppxf_fit is None else ppxf_fit.chi2

        self.bestfit_comp = None
        if ppxf_fit is None or not component_fits:
            return

        # Construct the model for each component
        component_weights = ppxf_fit.weights[None,:]*numpy.array([
                                                numpy.array(ppxf_fit.component == i).astype(int)
                                                        for i in range(ppxf_fit.ncomp)])
        self.bestfit_comp = numpy.dot(component_weights, ppxf_fit.matrix.T[degree+1:,:])

        # Include the multiplicative polynomial
        if self.multcoef is not None:
            mpoly = numpy.polynomial.legendre.legval(numpy.linspace(-1, 1, end-start),
                                                     numpy.append(1.0,self.multcoef))
            self.bestfit_comp[:,start:end] *= mpoly[None,:]

#        pyplot.plot(ppxf_fit.bestfit, color='k')
#        for i in range(ppxf_fit.ncomp):
#            pyplot.plot(self.bestfit_comp[i,:])
#        pyplot.show()


    def empty_fit(self):
        return self.status is None

    def reached_maxiter(self):
        return self.status == 5

    def fit_failed(self):
        return self.empty_fit() or self.status <= 0


class PPXFModel:
    """
    Class that reconstructs a pPXF model given a set of templates and
    the model parameters.

    This pulls functions from M. Cappellari's ppxf class, version 6.7.0.

    Input common to ppxf() input should have the same format as pPXF
    input.

    The only reason the galaxy spectrum (or spectra) are provided is to
    set if the data are to fit as reflection-symmetric and how many
    spectral pixels there are.

    """
    def __init__(self, templates, galaxy, velscale, velscale_ratio=None, vsyst=None, sigma_diff=0,
                 moments=2, degree=4, mdegree=0, component=0, gas_component=None, lam=None,
                 reddening=None, gas_reddening=None, reddening_func=None, sky=None,
                 templates_rfft=None, trig=False, quiet=False):

        # Suppress terminal output
        self.quiet = quiet

        # Dimensions of spectrum (or spectra) to fit
        self.nspec = galaxy.ndim     # nspec=2 for reflection-symmetric LOSVD
        self.npix = galaxy.shape[0]  # total pixels in the galaxy spectrum
        if not numpy.all(numpy.isfinite(galaxy)):
            raise ValueError('GALAXY data must be finite')
        if self.nspec == 2 and vsyst is None:
            raise ValueError('VSYST must be defined for two-sided fitting')

        # Set templates to use
        self.templates = templates.reshape(templates.shape[0], -1)
        self.npix_temp, self.ntemp = self.templates.shape

        # Pixel sampling
        self.velscale = velscale
        self.vsyst = 0 if vsyst is None else vsyst/velscale
        self.sigma_diff = sigma_diff/velscale
        self.factor = 1   # default value
        if velscale_ratio is not None:
            if not isinstance(velscale_ratio, int):
                raise ValueError('VELSCALE_RATIO must be an integer')
            self.npix_temp -= self.npix_temp % velscale_ratio
            # Make size multiple of velscale_ratio
            self.templates = self.templates[:self.npix_temp, :]
            # This is the size after rebin()
            self.npix_temp //= velscale_ratio
            self.factor = velscale_ratio
        if self.npix_temp < galaxy.shape[0]:
            raise ValueError('TEMPLATES length cannot be smaller than GALAXY')

        # Set the FFT of the templates
        self.npad = 2**int(numpy.ceil(numpy.log2(self.templates.shape[0])))
#        self.npad = fftpack.next_fast_len(self.templates.shape[0])
        self.templates_rfft = numpy.fft.rfft(self.templates, self.npad, axis=0) \
                                if templates_rfft is None else templates_rfft

        # Set sky spectrum
        if sky is None:
            self.sky = None
        else:
            if sky.shape[0] != galaxy.shape[0]:
                raise ValueError('GALAXY and SKY must have the same size')
            self.sky = sky.reshape(sky.shape[0], -1)
        self.nsky = 0 if self.sky is None else self.sky.shape[1]

        # Select polynomials to use
        self.degree = max(degree, -1)
        self.mdegree = max(mdegree, 0)
        if trig:
            if self.degree > -1 and self.degree % 2 != 0:
                raise ValueError('`degree` must be even with `trig=True`')
            if self.mdegree > -1 and self.mdegree % 2 != 0:
                raise ValueError('`mdegree` must be even with `trig=True`')
            self.polyval = ppxf.trigval
            self.polyvander = ppxf.trigvander
        else:
            self.polyval = numpy.polynomial.legendre.legval
            self.polyvander = numpy.polynomial.legendre.legvander

        # Initialize the design matrix.  Additive polynomial components
        # are independent of the input parameters.
        self.npoly = (self.degree + 1)*self.nspec  # Number of additive polynomials in fit
        self.ndmcols = self.npoly + self.nsky*self.nspec + self.ntemp
        self.matrix = numpy.zeros((self.npix*self.nspec, self.ndmcols))
        if self.degree > -1:
            vand = self.polyvander(numpy.linspace(-1, 1, self.npix), self.degree)
            self.matrix[: self.npix, : self.npoly//self.nspec] = vand
            if self.nspec == 2:     # poly for right spectrum
                self.matrix[self.npix :, self.npoly//self.nspec : self.npoly] = vand

        # Kinematic components
        self.component = numpy.atleast_1d(component)
        if self.component.dtype != int:
            raise TypeError('COMPONENT must be integers')
        if self.component.size == 1 and self.ntemp > 1:
            # component is a scalar so all templates have the same LOSVD
            self.component = numpy.zeros(self.ntemp, dtype=int)
        elif self.component.size != self.ntemp:
            raise ValueError('There must be one kinematic COMPONENT per template')
        tmp = numpy.unique(self.component)
        self.ncomp = tmp.size
        if not numpy.array_equal(tmp, numpy.arange(self.ncomp)):
            raise ValueError('COMPONENT must range from 0 to NCOMP-1')

        # Isolate gas components
        if gas_component is None:
            self.gas_component = numpy.zeros(self.component.size, dtype=bool)
            self.gas_any = False
        else:
            self.gas_component = numpy.asarray(gas_component)
            if self.gas_component.dtype != bool:
                raise TypeError('`gas_component` must be boolean')
            if self.gas_component.size != component.size:
                raise ValueError('`gas_component` and `component` must have the same size')
            if not numpy.any(gas_component):
                raise ValueError('At least one `gas_component` must be set to true')
            self.gas_any = True

        # Set reddening
        if lam is None and (reddening is not None or gas_reddening is not None):
            raise ValueError('LAM must be given with REDDENING or GAS_REDDENING keyword')
        if self.mdegree > 0 and reddening is not None:
            raise ValueError('MDEGREE cannot be used with REDDENING keyword')
        if gas_reddening is not None and not self.gas_any:
            raise ValueError('GAS_COMPONENT must be nonzero with GAS_REDDENING keyword')
        if lam is not None and lam.shape != galaxy.shape:
            raise ValueError('GALAXY and LAM must have the same size')
        self.lam = lam
        self.reddening = reddening
        self.gas_reddening = gas_reddening
        if reddening_func is None:
            self.reddening_func = ppxf.reddening_cal00
        else:
            if not callable(reddening_func):
                raise TypeError('`reddening_func` must be callable')
            self.reddening_func = reddening_func

        # Kinematic moments to fit
        self.moments = numpy.atleast_1d(moments)
        if self.moments.size == 1:
            # moments is scalar: all LOSVDs have same number of G-H moments
            self.moments = numpy.full(self.ncomp, self.moments, dtype=int)
        self.fixall = self.moments < 0  # negative moments --> keep entire LOSVD fixed
        self.moments = numpy.abs(self.moments)
        if self.moments.size != self.ncomp:
            raise ValueError('MOMENTS must be an array of length NCOMP')


    def __call__(self, kin, tplwgts, addpoly=None, multpoly=None, reddening=None,
                 gas_reddening=None):
        """Same as :func:`construct`"""
        return self.construct(kin, tplwgts, addpoly=addpoly, multpoly=multpoly,
                              reddening=reddening, gas_reddening=gas_reddening)


    def construct(self, kin, tplwgts, addpoly=None, multpoly=None, reddening=None,
                  gas_reddening=None):
        """
        Construct a pPXF model spectrum provided the input
        parameters.  Mostly a copy of ppxf._linear_fit().

        Args:
            kin (list, numpy.ndarray): Must have the kinematics of each
                component.  Must be a list or array of vectors with a
                shape (NCOMP,).  The length of each vector in the
                list/array must be the same as MOMENTS.
            tplwgts (numpy.ndarray): Weights to apply to each template.
                Shape must be (NTEMP,).
            addpoly (numpy.ndarray): (**Optional**) Coefficients of the
                additive polynomials.  Shape must be
                (NSPEC*(DEGREE+1),).  Exception is raised if
                coefficients are expected (self.degree > -1), but no
                coefficiencts are provided.
            multpoly (numpy.ndarray): (**Optional**) Coefficients of the
                multiplicative polynomials.  Shape must be
                (NSPEC*MDEGREE,).  Exception is raised if coefficients
                are expected (self.mdegree > 0), but no coefficiencts
                are provided.
            reddening (float): (**Optional**) E(B-V) value for the
                continuum fit.  Exception is raised if value is expected
                (self.reddening is not None), but none is provided.
            gas_reddening (float): (**Optional**) E(B-V) value for the
                gas components.  Exception is raised if value is
                expected (self.reddening is not None), but none is
                provided.

        Returns:

            numpy.ndarray:
        
        """

        # Check the input kinematics
        params = [kin] if self.ncomp == 1 else list(kin)
        if len(params) != self.ncomp:
            raise ValueError('There must be one set of kinematic moments per component')
        if not numpy.all([hasattr(p, "__len__") and len(p) == m 
                                for p,m in zip(params, self.moments)]):
            raise ValueError('Input kinematics have the incorrect length')
        # Convert velocity from km/s to pixels
        for j, flx in enumerate(params):
            params[j][:2] = flx[:2]/self.velscale
        params = numpy.concatenate(params)

        # Check the optional input
        if self.degree > -1 and (addpoly is None or len(addpoly) != (self.degree + 1)*self.nspec):
            raise ValueError('addpoly must be provided with length {0}'.format(
                                                                    (self.degree+1)*self.nspec))
        if self.mdegree > 0 and (multpoly is None or len(multpoly) != self.mdegree*self.nspec):
            raise ValueError('multpoly must be provided with length {0}'.format(
                                                                    self.mdegree*self.nspec))
        if self.reddening is not None and reddening is None:
            raise ValueError('reddening must be provided')
        if self.gas_reddening is not None and gas_reddening is None:
            raise ValueError('gas_reddening must be provided')

        # Check the input template weights
        weights = numpy.atleast_1d(tplwgts)
        if weights.size != self.ntemp:
            raise ValueError('Incorrect number of template weights provided.')
        # Include the additive polynomial weights
        if self.degree > -1:
            weights = numpy.append(addpoly, weights)

        # Get the real FFT of the LOSVDs
        losvd_rfft = ppxf.losvd_rfft(params, self.nspec, self.moments,
                                     self.templates_rfft.shape[0], self.ncomp, self.vsyst,
                                     self.factor, self.sigma_diff)

        # Get the multiplicative polynomials
        x = numpy.linspace(-1, 1, self.npix)
        if self.mdegree > 0:
            if self.nspec == 2: # Different multiplicative poly for left/right spectra
                mpoly1 = self.polyval(x, numpy.append(1.0, multpoly[::2]))
                mpoly2 = self.polyval(x, numpy.append(1.0, multpoly[1::2]))
                mpoly = numpy.append(mpoly1, mpoly2).clip(0.1)
            else:
                mpoly = self.polyval(x, numpy.append(1.0, multpoly)).clip(0.1)
        else:
            mpoly = None if reddening is None else self.reddening_func(self.lam, reddening)
        gas_mpoly = None if gas_reddening is None else self.reddening_func(self.lam, gas_reddening)

        # Perform the LOSVD convolution and finish the design matrix
        tmp = numpy.empty((self.nspec, self.npix_temp))
        for j, template_rfft in enumerate(self.templates_rfft.T):  # columns loop
            for k in range(self.nspec):
                pr = template_rfft*losvd_rfft[:, self.component[j], k]
                tt = numpy.fft.irfft(pr, self.npad)
                tmp[k, :] = ppxf.rebin(tt[:self.npix_temp*self.factor], self.factor)
            self.matrix[:, self.npoly + j] = tmp[:, :self.npix].ravel()
            if self.gas_component[j]:
                if gas_mpoly is not None:
                    self.matrix[:, self.npoly + j] *= gas_mpoly
            elif mpoly is not None:
                self.matrix[:, self.npoly + j] *= mpoly

        if self.nsky > 0:
            k = self.npoly + self.ntemp
            self.matrix[: self.npix, k : k + self.nsky] = self.sky
            if nspec == 2:      # Sky for right spectrum
                self.matrix[self.npix :, k + self.nsky : k + 2*self.nsky] = self.sky

        # return the model
        return self.matrix.dot(weights)


class PPXFFit(StellarKinematicsFit):
    """
    Use pPXF to measure the stellar kinematics.  Although it can also
    fit the composition and emission lines, for now just force it to be
    a :class:`StellarKinematicsFit` objec.

    Attributes:
        bitmask (BitMask): Bitmask to use for masking spectral pixels
            and parameter values.  For
            :func:`fit_SpatiallyBinnedSpectra`, must have bit for
            'LOW_SNR'.  For :func:`fit` must have bits for 'TPL_PIXELS',
            'TRUNCATED', 'PPXF_REJECT', 'LARGE_CHI2', 'LARGE_RESID',
            'INSUFFICIENT_DATA', 'FIT_FAILED', 'NEAR_BOUND'.
    """
    # Define bit names as global to the class
    snr_flag = 'LOW_SNR'
    rng_flag = 'OUTSIDE_RANGE'
    tpl_flag = 'TPL_PIXELS'
    trunc_flag = 'TRUNCATED'
    rej_flag = 'PPXF_REJECT'
    def __init__(self, bitmask, par=None):
        StellarKinematicsFit.__init__(self, 'ppxf', bitmask, par=par)
        # Logging and terminal output
        self.loggers = None
        self.quiet = False

        # Imposed fitting boundaries
        self.velocity_limits = None
        self.sigma_limits = None
        self.gh_limits = None

        # Spectral data
        self.tpl_wave = None
        self.tpl_sres = None
        self.tpl_flux = None
        self.ntpl = None
        self.npix_tpl = None

        self.obj_wave = None
        self.obj_sres = None
        self.obj_flux = None
        self.obj_ferr = None
        self.input_obj_mask = None
        self.obj_to_fit = None
        self.nobj = None
        self.npix_obj = None

        # Fitting parameters
        self.velscale = None
        self.velscale_ratio = None
        self.matched_resolution = None
        self.guess_kin = None
        self.spectrum_start = None
        self.spectrum_end = None
        self.dof = None
        self.filt_dof = None
        self.base_velocity = None

        # Fitting options
        self.iteration_mode = None
        self.reject_boxcar = None
        self.filter_boxcar = None
        self.filter_operation = None
        self.filter_iterations = None
        self.bias = None
        self.degree = None
        self.mdegree = None
        self.filt_degree = None
        self.filt_mdegree = None
        self.moments = None

        self.fix_kinematics = None

    @staticmethod
    def iteration_modes():
        r"""
        Possible iteration methods:

          * ``none``: Fit all bins with all templates with a single call
            to pPXF.

          * ``fit_reject_filter``: Perform the following procedure:
                - Fit each spectrum
                - for n iterations:
                    - Reject outliers
                    - Filter the object and template spectra
                    - Fit the filtered spectra
                - Fit the unfiltered spectra with the kinematics fixed
                  to result of the final filtered fit

          * ``no_global_wrej``: Do not fit the global spectrum
            first, but include a rejection iteration.  All templates are
            fit in each step.

          * ``global_template``:  Fit the global spectrum with all
            templates and include a single rejection iteration.  The
            pixel mask for this fit is the base mask for all fits to the
            individual bins.  A single rejection iteration is done for
            each bin.  **Only the global template is used when fitting
            each bin.**

          * ``nonzero_templates``:  Fit the global spectrum with
            all templates and include a single rejection iteration.  The
            pixel mask for this fit is the base mask for all fits to the
            individual bins.  A single rejection iteration is done for
            each bin.  **Only the templates with non-zero weights are
            used when fitting each bin.**

          * ``all_templates``:  Fit the global spectrum with all
            templates and include a single rejection iteration.  The
            pixel mask for this fit is the base mask for all fits to the
            individual bins.  A single rejection iteration is done for
            each bin.  **All templates are used when fitting each bin.**

        Returns:
            list: List of allowed options.
        """
        return [ 'none',
                 'fit_reject_filter',
                 'no_global_wrej',
                 'global_template',
                 'nonzero_templates',
                 'all_templates']

    def _mode_uses_global_spectrum(self):
        return self.iteration_mode in ['global_template', 'nonzero_templates', 'all_templates']

    def _mode_uses_global_template(self):
        return self.iteration_mode in ['global_template']

            
    def _mode_uses_nonzero_templates(self):
        return self.iteration_mode in ['nonzero_templates']


    def _mode_uses_all_templates(self):
        return self.iteration_mode in ['none', 'fit_reject_filter', 'no_global_wrej',
                                       'all_templates']

    def _mode_includes_rejection(self):
        return self.iteration_mode != 'none'

    def _mode_uses_filter(self):
        return self.iteration_mode == 'fit_reject_filter'

    def _check_mode(self, iteration_mode, reject_boxcar, filter_boxcar, filter_operation,
                    filter_iterations, mdegree):
        if iteration_mode not in self.iteration_modes():
            raise ValueError('Do not understand iteration mode \'{0}\''.format(iteration_mode))
        self.iteration_mode = iteration_mode
        self.filter_iterations = 0 if self.iteration_mode != 'fit_reject_filter' \
                                        or filter_iterations is None else filter_iterations
        if self.iteration_mode == 'fit_reject_filter' and filter_boxcar is None:
            warnings.warn('Must provide boxcar for filtering iterations.  Using default (100).')
            self.filter_boxcar = 100
        else:
            self.filter_boxcar = filter_boxcar
        self.reject_boxcar = reject_boxcar
        if self.filter_iterations > 0:
            if filter_operation not in [ 'divide', 'subtract']:
                warnings.warn('Must specify the filtering operation as either \'divide\' or '
                              '\'subtract\'; using divide as default.')
                self.filter_operation = 'divide'
            else:
                self.filter_operation = filter_operation

    @staticmethod
    def check_template_usage_flags(nobj, ntpl, usetpl):
        if usetpl is not None and usetpl.shape != (nobj, ntpl) \
                and usetpl.shape != (ntpl,):
            raise ValueError('Provided template selection object does not have the correct shape!')
        if usetpl is None:
            _usetpl = numpy.ones((nobj,ntpl), dtype=numpy.bool)
        else:
            _usetpl = usetpl.astype(bool)
            if _usetpl.shape == (ntpl,):
                _usetpl = numpy.array([_usetpl]*nobj)
        return _usetpl


    @staticmethod
    def check_input_kinematics(nobj, guess_redshift, guess_dispersion):
        # Get the input guess kinematics
        _guess_redshift = numpy.atleast_1d(guess_redshift)
        if len(_guess_redshift) != 1 and len(_guess_redshift) != nobj:
            raise ValueError('Must provide a single redshift or one per object spectrum.')
        if len(_guess_redshift) == 1:
            _guess_redshift = numpy.full(nobj, _guess_redshift[0], dtype=numpy.float)
        _guess_dispersion = numpy.atleast_1d(guess_dispersion)
        if len(_guess_dispersion) != 1 and len(_guess_dispersion) != nobj:
            raise ValueError('Must provide a single dispersion or one per object spectrum.')
        if len(_guess_dispersion) == 1:
            _guess_dispersion = numpy.full(nobj, _guess_dispersion[0], dtype=numpy.float)

        # Set the input redshifts
        input_cz = _guess_redshift*astropy.constants.c.to('km/s').value

        # Set the input pPXF, pixel-based kinematics (see
        # convert_velocity)
        guess_kin = numpy.array([numpy.log(_guess_redshift+1)*astropy.constants.c.to('km/s').value,
                                    _guess_dispersion]).T
        return input_cz, guess_kin


    @staticmethod
    def check_resolution_match(tpl_sres, obj_sres, matched_resolution):
        # Confirm there is enough information to handle an unmatched
        # resolution
        if not matched_resolution and (obj_sres is None or tpl_sres is None):
            raise ValueError('If the spectral resolution is not matched between the template and '
                             'the object data, you must provide the spectral resolution for both.')
        return matched_resolution


    @staticmethod
    def set_wavelength_range(nobj, obj_wave, waverange=None):
        _waverange = numpy.array([[obj_wave[0]-1, obj_wave[-1]+1]]) \
                                if waverange is None else numpy.atleast_2d(waverange)
        _waverange[numpy.equal(_waverange[:,0], None),0] = obj_wave[0]-1
        _waverange[numpy.equal(_waverange[:,1], None),1] = obj_wave[-1]+1
        if _waverange.shape[0] == 1:
            _waverange = numpy.array([_waverange[0,:]]*nobj)
        if _waverange.shape != (nobj,2):
            raise ValueError('Input wavelength range array does not have the correct shape.')
        return _waverange


    @staticmethod
    def initialize_model_mask(obj_wave, obj_flux, mask=None, bitmask=None, velocity_offset=None):
        model_mask = numpy.zeros(obj_flux.shape,
                                 dtype=bool if bitmask is None else bitmask.minimum_dtype())
        if mask is None:
            return model_mask

        # Include the mask, if provided
        if isinstance(mask, numpy.ndarray):
            if mask is not None and mask.shape != obj_flux.shape:
                raise ValueError('Shape of object mask array must match its flux array.')
            model_mask = mask if bitmask is None else mask.astype(bitmask.minimum_dtype())
        if isinstance(mask, SpectralPixelMask):
            model_mask = mask.boolean(obj_wave, nspec=obj_flux.shape[0],
                                      velocity_offsets=velocity_offset) if bitmask is None \
                            else mask.bits(bitmask, obj_wave, nspec=obj_flux.shape[0],
                                           velocity_offsets=velocity_offset)
        return model_mask


    @staticmethod
    def initialize_pixels_to_fit(tpl_wave, obj_wave, obj_flux, obj_ferr, velscale,
                                 velscale_ratio=None, waverange=None, mask=None, bitmask=None,
                                 velocity_offset=None, max_velocity_range=400., alias_window=None,
                                 ensemble=True, loggers=None, quiet=False):
        """
        obj_flux and obj_ferr *must* be masked arrays.

        .. todo::

            - This function can alter obj_flux and obj_ferr !!  Return
              them instead?

        """
        # Initialize the mask (does not include any mask in the obj_flux
        # object)
        nobj = obj_flux.shape[0]
        err = numpy.empty(nobj, dtype=object)  # Empty error messages
        _waverange = PPXFFit.set_wavelength_range(nobj, obj_wave, waverange=waverange)
        model_mask = PPXFFit.initialize_model_mask(obj_wave, obj_flux, mask=mask, bitmask=bitmask,
                                                   velocity_offset=velocity_offset)

        # Include the masked object pixels
        model_mask[obj_flux.mask] = True if bitmask is None else \
                                        bitmask.turn_on(model_mask[obj_flux.mask], 'DIDNOTUSE')

        # Assess the regions that need to be masked during fitting
        if ensemble:
            # Treat all spectra as part of an ensemble so mask the same
            # region for all spectra
            _waverange = numpy.array([numpy.amax(_waverange[:,0]), numpy.amin(_waverange[:,1])])
            try:
                fit_indx, waverange_mask, npix_mask, alias_mask \
                        = PPXFFit.fitting_mask(tpl_wave, obj_wave, velscale,
                                               velscale_ratio=velscale_ratio,
                                               waverange=_waverange,
                                               velocity_offset=velocity_offset,
                                               max_velocity_range=max_velocity_range,
                                               alias_window=alias_window, loggers=loggers,
                                               quiet=True)#quiet)
            except ValueError as e:
                return model_mask, numpy.array([e]*nobj), numpy.zeros(nobj, dtype=int), \
                            numpy.array([obj_flux.shape[1]]*nobj).astype(int)

            fit_indx = numpy.array([ fit_indx ]*nobj)
            waverange_mask = numpy.array([ waverange_mask ]*nobj)
            npix_mask = numpy.array([ npix_mask ]*nobj)
            alias_mask = numpy.array([ alias_mask ]*nobj)
        else:
            # Treat all spectra independently, such that the masks are
            # independent
            fit_indx = numpy.zeros(obj_flux.shape, dtype=bool)
            waverange_mask = numpy.zeros(obj_flux.shape, dtype=bool)
            npix_mask = numpy.zeros(obj_flux.shape, dtype=bool)
            alias_mask = numpy.zeros(obj_flux.shape, dtype=bool)
            _velocity_offset = numpy.asarray(velocity_offset) \
                                        if isinstance(velocity_offset,(list,numpy.ndarray)) \
                                        else numpy.array([velocity_offset]*nobj)
            if len(_velocity_offset) != nobj:
                raise ValueError('Incorrect number of velocity offsets provided.')
            for i in range(nobj):
                try:
                    fit_indx[i,:], waverange_mask[i,:], npix_mask[i,:], alias_mask[i,:] \
                            = PPXFFit.fitting_mask(tpl_wave, obj_wave, velscale,
                                                   velscale_ratio=velscale_ratio,
                                                   waverange=_waverange[i,:],
                                                   velocity_offset=_velocity_offset[i],
                                                   max_velocity_range=max_velocity_range,
                                                   alias_window=alias_window, loggers=loggers,
                                                   quiet=True)#quiet)
                except ValueError as e:
                    err[i] = e

        # Add these regions to the mask
        if bitmask is None:
            model_mask[waverange_mask] = True
            model_mask[npix_mask] = True
            model_mask[alias_mask] = True
        else:
            model_mask[waverange_mask] = bitmask.turn_on(model_mask[waverange_mask],
                                                         PPXFFit.rng_flag)
            model_mask[npix_mask] = bitmask.turn_on(model_mask[npix_mask], PPXFFit.tpl_flag)
            model_mask[alias_mask] = bitmask.turn_on(model_mask[alias_mask], PPXFFit.trunc_flag)

        # Make sure that the errors are valid
        indx = numpy.invert(obj_ferr.data > 0) | numpy.invert(numpy.isfinite(obj_ferr.data))
        if numpy.sum(indx) > 0:
            model_mask[indx] = True if bitmask is None else\
                                    bitmask.turn_on(model_mask[indx], 'INVALID_ERROR')
        # To avoid having pPXF throw an exception, make sure that all
        # of these bad errors are set to unity
        obj_ferr[indx] = 1.0

        # Update the internal mask of the data
        obj_flux[model_mask > 0] = numpy.ma.masked
        obj_ferr[model_mask > 0] = numpy.ma.masked

        # Determine the starting and ending pixels and return
        # TODO: does this work if all the pixels are masked in a given
        # spectrum?
        pix = numpy.ma.MaskedArray(numpy.array([numpy.arange(obj_flux.shape[1])]*nobj),
                                   mask=numpy.invert(fit_indx))
        return model_mask, err, numpy.ma.amin(pix, axis=1), numpy.ma.amax(pix, axis=1)+1


    def _run_fit_iteration(self, obj_flux, obj_ferr, start, end, base_velocity, tpl_flux,
                           tpl_rfft, guess_kin, fix_kinematics=False, obj_to_fit=None,
                           tpl_to_use=None, degree=None, mdegree=None, dof=None,
                           weight_errors=False, plot=False):
        r"""
        Fit all the object spectra in obj_flux.

        .. todo::
            - Calculate DOF in this function?
            - Explicitly set the bounds to use instead of using the pPXF
              defaults?

        Args:
            obj_flux (numpy.ma.MaskedArray): Size is :math:`N_{\rm
                spec}\times N_{\rm chan}`, object spectra.
            obj_ferr (numpy.ma.MaskedArray): Size is :math:`N_{\rm
                spec}\times N_{\rm chan}`, object errors
            start (array): Size is :math:`N_{\rm spec}`, starting pixel
                for each spectrum
            end (array): Size is :math:`N_{\rm spec}`, ending pixel (+1)
                for each spectrum
            base_velocity (array): Size is :math:`N_{\rm spec}`, base
                velocity offset between each object spectrum and the
                template spectra.
            tpl_flux (array): Size is :math:`N_{\rm tpl}\times N_{\rm
                tpl chan}`, template spectra
            tpl_rfft (array): Size is :math:`N_{\rm tpl}\times N_{\rm
                tpl pad}`, real FFT of the template spectra
            guess_kin (array): Initial guess for kinematics.  Size is
                :math:`N_{\rm spec}\times N_{\rm moments}`.
            fix_kinematics (bool): (**Optional**) Flag to fix the
                kinematics to the input values during the fit.
            obj_to_fit (array): (**Optional**) Size is :math:`N_{\rm
                spec}`, boolean flag to fit object spectrum
            tpl_to_use (array): (**Optional**) Size is :math:`N_{\rm
                spec}\times N_{\rm tpl}`, boolean flag to use a template
                for the fit to each object spectrum
            plot (bool): (**Optional**) Produce the default ppxf fit
                plot.
            degree (int): (**Optional**) Additive polynomial order.
                Default is to use the internal attribute :attr:`degree`.
            mdegree (int): (**Optional**) Multiplicative polynomial
                order.  Default is to use the internal attribute
                :attr:`mdegree`.
            dof (int): (**Optional**) Number of degrees of freedom in
                the fit.  Default is to use the internal attribute
                :attr:`dof`.

        Returns:
            numpy.ndarray : Array with :math:`N_{\rm spec}` instances of
            :class:`PPXFFitResult`.
        """
        # Get the list of templates to use
        nspec = obj_flux.shape[0]
        _obj_to_fit = numpy.ones(nspec, dtype=bool) if obj_to_fit is None else obj_to_fit
        ntpl = tpl_flux.shape[0]
        _tpl_to_use = numpy.ones((nspec,ntpl), dtype=bool) if tpl_to_use is None else tpl_to_use

#        input_kin = self.guess_kin if fixed_kin is None else fixed_kin
        moments = -self.moments if fix_kinematics else self.moments
        degree = self.degree if degree is None else degree
        mdegree = self.mdegree if mdegree is None else mdegree
        dof = self.dof if dof is None else dof

#        linear = fixed_kin is not None and mdegree < 1
        linear = fix_kinematics and mdegree < 1

        # Create the object to hold all the fits
        result = numpy.empty(nspec, dtype=object)

#        plot=True

        # Fit each spectrum
        for i in range(nspec):
            print('Running pPXF fit on spectrum: {0}/{1}'.format(i+1,nspec), end='\r')
            # Meant to ignore this spectrum
            if not _obj_to_fit[i]:
                result[i] = None
                continue

            # Get the pixels to fit for this spectrum
            gpm = numpy.where(numpy.invert(obj_flux.mask[i,start[i]:end[i]]))[0]

            # Check if there is sufficient data for the fit
            ntpl_to_use = numpy.sum(_tpl_to_use[i,:])
            if len(gpm) < dof+ntpl_to_use:
                if not self.quiet:
                    warnings.warn('Insufficient data points ({0}) to fit spectrum {1}'
                                  '(dof={2}).'.format(len(gpm), i+1, dof+ntpl_to_use))
                result[i] = PPXFFitResult(degree, mdegree, start[i], end[i], _tpl_to_use[i,:],
                                          None, ntpl)
                continue

            # Run ppxf
            if plot:
                pyplot.clf()
            result[i] = PPXFFitResult(degree, mdegree, start[i], end[i], _tpl_to_use[i,:],
                            ppxf.ppxf(tpl_flux[tpl_to_use[i,:],:].T,
                                      obj_flux.data[i,start[i]:end[i]],
                                      obj_ferr.data[i,start[i]:end[i]], self.velscale,
                                      guess_kin[i,:], velscale_ratio=self.velscale_ratio,
                                      goodpixels=gpm, bias=self.bias, degree=degree,
                                      mdegree=mdegree, moments=moments, vsyst=-base_velocity[i],
                                      quiet=(not plot), plot=plot, linear=linear,
                                      templates_rfft=tpl_rfft[tpl_to_use[i,:],:].T), ntpl,
                                      weight_errors=weight_errors)
#                                      linear_method='lsqlin'), ntpl,
#                                      weight_errors=weight_errors)

            if result[i].kin[1] < 0:
#                result[i].kin[1] = numpy.absolute(result[i].kin[1]) #self.sigma_limits[0]
#                warnings.warn('pPXF gives negative dispersion! Change -{0:.4f} to {0:.4f}'.format(
#                                    result[i].kin[1]))
                raise ValueError('Dispersion less than 0! {0}/{1} {2}'.format(
                                        i+1,nspec,result[i].kin[1]))
                
            if result[i].reached_maxiter() and not self.quiet:
                warnings.warn('pPXF optimizer reached maximum number of iterations for spectrum '
                              '{0}.'.format(i+1))
            if plot:
                pyplot.show()

#            model = PPXFModel(tpl_flux[tpl_to_use[i,:],:].T, obj_flux.data[i,start[i]:end[i]],
#                              self.velscale, velscale_ratio=self.velscale_ratio, degree=degree,
#                              mdegree=mdegree, moments=moments, vsyst=-base_velocity[i],
#                              quiet=True, templates_rfft=tpl_rfft[tpl_to_use[i,:],:].T)
#           
#            modelfit = model(result[i].kin, result[i].tplwgt, addpoly=result[i].addcoef,
#                             multpoly=result[i].multcoef)
#
#            print(numpy.sum(result[i].bestfit-modelfit))
#
#            pyplot.plot(self.obj_wave[start[i]:end[i]], modelfit, color='C0')
#            pyplot.plot(self.obj_wave[start[i]:end[i]], result[i].bestfit, color='C3')
#            pyplot.plot(self.obj_wave[start[i]:end[i]], result[i].bestfit-modelfit, color='C1')
#            pyplot.show()

        print('Running pPXF fit on spectrum: {0}/{1}'.format(nspec,nspec))
        return result


    def _fit_global_spectrum(self, obj_to_include=None, plot=False):
        """
        Fit the global spectrum.  This:
            - Sets the base-level good pixel mask for the fits to the individual
              bins
            - Gives the template weights for the global template

        .. todo::
            - Only include spectra above a given S/N in global spectrum?
            - Allow for a number of iterations as input.
        """

        _obj_to_include = numpy.ones(self.nobj, dtype=bool) \
                            if obj_to_include is None else obj_to_include
        if len(_obj_to_include) != self.nobj:
            raise ValueError('Incorrect number of object flags.')

        # Create the global spectrum and its error
        # TODO: Does not include covariance!
        global_spectrum = numpy.ma.sum(self.obj_flux[_obj_to_include,:], axis=0).reshape(1,-1)
        global_spectrum_err = numpy.ma.sqrt(numpy.ma.sum(
                                                numpy.square(self.obj_ferr[_obj_to_include,:]),
                                                         axis=0)).reshape(1,-1)
        global_spectrum_err[numpy.ma.getmaskarray(global_spectrum)] = 1.0   # To avoid pPXF error

        # TODO: Because of how it's used, setting start, end, and
        # base_vel this way will mess things up later in the fit()
        # function UNLESS all the spectra have the same start and end;
        # ie., when fitting the global spectrum, the object spectra
        # provided to fit() must be treated as an ensemble.

        # Set the fitting region and base velocity offset
        start = numpy.array([numpy.amax(self.spectrum_start)])
        base_vel = numpy.array([self.base_velocity[numpy.argmax(self.spectrum_start)]])
        end = numpy.array([numpy.amin(self.spectrum_end)])

        # Use any template that is request for any of the individual
        # spectra
        usetpl = numpy.any(self.usetpl, axis=0).reshape(1,-1)

#        plot=True
        # Fit the spectrum
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, 'First fit to global spectrum.')
        result = self._run_fit_iteration(global_spectrum, global_spectrum_err, start, end,
                                         base_vel, self.tpl_flux, self.tpl_rfft, self.guess_kin,
                                         fix_kinematics=self.fix_kinematics, tpl_to_use=usetpl,
                                         plot=plot)

        # Return if pPXF failed
        if result[0].fit_failed():
            if not self.quiet:
                log_output(self.loggers, 1, logging.INFO, 'pPXF failed to fit global spectrum!')
            return None, None

        # Perform a single fit rejection
        global_spectrum = PPXFFit.reject_model_outliers(global_spectrum, result, rescale=True,
                                                        loggers=self.loggers, quiet=self.quiet)

        # refit the spectrum
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, 'Fit to global spectrum after rejection.')
        return self._run_fit_iteration(global_spectrum, global_spectrum_err, start, end,
                                       base_vel, self.tpl_flux, self.tpl_rfft, self.guess_kin,
                                       fix_kinematics=self.fix_kinematics, tpl_to_use=usetpl,
                                       plot=plot)[0]


    def _fill_ppxf_par(self, kin, no_shift=True):

        # Moments for each kinematic component
        ncomp = 1
        moments = numpy.atleast_1d(self.moments)
        _moments = numpy.full(ncomp, numpy.absolute(moments), dtype=int) if moments.size == 1 \
                              else numpy.absolute(moments)

        # Construct the LOSVD parameter vector
        vj = numpy.append(0, numpy.cumsum(_moments)[:-1])
        par = kin.copy()
        par[0:2] /= self.velscale
        if no_shift:
            par[0] = 0.0

        return par, _moments, vj


    def _get_losvd_kernels(self, result, no_shift=True):
        
        nspec = len(result)
        losvd_kernel_rfft = numpy.empty(nspec, dtype=object)
        for i in range(nspec):
            if result[i] is None or result[i].fit_failed():
                continue
            par, _moments, vj = self._fill_ppxf_par(result[i].kin, no_shift=no_shift)
            losvd_kernel_rfft[i] = ppxf.losvd_rfft(par, 1, _moments, self.tpl_rfft.shape[1], 1,
                                        0.0 if no_shift else -self.base_velocity[i]/self.velscale,
                                                   self.velscale_ratio, 0.0)[:,0,0]
        return losvd_kernel_rfft
            

    def _matched_mask_filter(self, bf, obj_mask, tpl_flux, tpl_rfft, tpl_to_use, result):
        """
        self.tpl_npad must exist!
        """
        # Get the LOSVD kernels for each fit
        losvd_kernel_rfft = self._get_losvd_kernels(result)
        npix_tpl_resampled = self.npix_tpl // self.velscale_ratio

        # Instantiate the output arrays
        ntpl_per_obj = numpy.sum(tpl_to_use, axis=1)
        tpl_flux_filt = numpy.ma.MaskedArray(numpy.ma.zeros((numpy.sum(ntpl_per_obj),
                                                             self.npix_tpl), dtype=float))
        tpl_to_use_filt = numpy.zeros((self.nobj, numpy.sum(ntpl_per_obj)), dtype=bool)

        # For each spectrum:
        for i in range(losvd_kernel_rfft.shape[0]):
            if result[i] is None:
                continue
            print('Masking and smoothing templates for object spectrum: {0}/{1}'.format(
                    i+1, losvd_kernel_rfft.shape[0]), end='\r')
#            t = time.perf_counter()

            # Get all the templates convolved by the LOSVD for this fit
            cnvlv_tpl_flux = numpy.fft.irfft(tpl_rfft[tpl_to_use[i],:]
                                                    * losvd_kernel_rfft[i][None,:],
                                             self.tpl_npad, axis=1)[:,:self.npix_tpl]
            if self.velscale_ratio > 1:
                cnvlv_tpl_flux = numpy.mean(cnvlv_tpl_flux.reshape(ntpl_per_obj[i], -1,
                                                                   self.velscale_ratio), axis=2)
#            print('fft: time: {0} seconds'.format(time.perf_counter() - t))

            # Get the object-spectrum mask shifted to the template frame
            _obj_mask = numpy.array([obj_mask[i,:]]*ntpl_per_obj[i], dtype=bool)
            shift = numpy.floor((result[i].kin[0]-self.base_velocity[i])/self.velscale).astype(int)
            if shift != 0:
                _obj_mask = numpy.roll(_obj_mask, shift)
                if shift > 0:
                    _obj_mask[:,:shift] = True
                else:
                    _obj_mask[:,shift:] = True
           
            # Apply the mask to the convolved spectra
            _cnvlv_tpl_flux = numpy.ma.MaskedArray(cnvlv_tpl_flux,
                                                   mask=_obj_mask[:,:npix_tpl_resampled])

            # Smooth the template spectra using the same smoothing
            # function as used for the object data
#            t2 = time.perf_counter()
            sm_cnvlv_tpl_flux = bf.smooth(_cnvlv_tpl_flux)
#            print('smooth: time: {0} seconds'.format(time.perf_counter() - t2))

            # Interpolate the smoothing function to the original pixels
            # of the template spectra
#            t2 = time.perf_counter()
            pixcoo = numpy.arange(self.npix_tpl*self.ntpl)
            interpolator = interpolate.interp1d(numpy.mean(pixcoo.reshape(-1,self.velscale_ratio),
                                                           axis=1), sm_cnvlv_tpl_flux.ravel(),
                                                fill_value='extrapolate', assume_sorted=True)
            sm_tpl_flux = interpolator(pixcoo).reshape(self.ntpl, -1)
#            print('interpolate: time: {0} seconds'.format(time.perf_counter() - t2))

            # Set the filtered templates for this object spectrum
            tpls = numpy.sum(ntpl_per_obj[:i])
            tple = numpy.sum(ntpl_per_obj[:i+1])
            tpl_to_use_filt[i,tpls:tple] = True
            if self.filter_operation == 'divide':
                tpl_flux_filt[tpls:tple,:] = numpy.ma.divide(tpl_flux[tpl_to_use[i],:],
                                                             sm_tpl_flux)
            else:
                tpl_flux_filt[tpls:tple,:] = tpl_flux[tpl_to_use[i],:] - sm_tpl_flux
#            print('obj: {0}, total time: {1} seconds'.format(i+1, time.perf_counter() - t))
        
        print('Masking and smoothing templates for object spectrum:              DONE')
        return tpl_flux_filt, tpl_to_use_filt


    def _fit_all_spectra(self, templates, templates_rfft, tpl_to_use, plot=False,
                         plot_file_root=None):
        """
        Fit all spectra provided.

        - Get an initial fit
        - Reject
        - Mask and smooth templates and objects
        - Fit ratio
        - Mask and smooth
        - Fit ratio
        - Fit unmasked with fixed kinematics to get models (with
          emission lines?)

        """
        #---------------------------------------------------------------
        # Fit the spectra
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO,
                       'Number of object spectra to fit: {0}/{1}'.format(
                            numpy.sum(self.obj_to_fit), len(self.obj_to_fit)))

        result = self._run_fit_iteration(self.obj_flux, self.obj_ferr, self.spectrum_start,
                                         self.spectrum_end, self.base_velocity, templates,
                                         templates_rfft, self.guess_kin,
                                         fix_kinematics=self.fix_kinematics,
                                         obj_to_fit=self.obj_to_fit, tpl_to_use=tpl_to_use,
                                         weight_errors=not self._mode_includes_rejection(),
                                         plot=plot)
        if not self._mode_includes_rejection():
            # Only a single fit so return
            return result

        #---------------------------------------------------------------
        # Copy the input as to not overwrite the input masks
        obj_flux = self.obj_flux.copy()
        obj_ferr = self.obj_ferr.copy()
        obj_to_fit = self.obj_to_fit.copy()
        nspec = obj_flux.shape[0]

        # Save which were not fit successfully
        obj_to_fit &= numpy.invert(numpy.array([ r is None or r.fit_failed() for r in result ]))
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO,
                       'Number of object spectra to fit (excluding failed fits): {0}/{1}'.format(
                            numpy.sum(self.obj_to_fit), len(self.obj_to_fit)))

        #---------------------------------------------------------------
        # Reject model outliers
        # TODO: This will cause an error if boxcar is None!
        obj_flux = PPXFFit.reject_model_outliers(obj_flux, result, rescale=False,
                                                 local_sigma=True, boxcar=self.reject_boxcar,
                                                 loggers=self.loggers, quiet=self.quiet)
        # Copy the new mask to the errors
        obj_ferr[numpy.ma.getmaskarray(obj_flux)] = numpy.ma.masked

        if self.filter_iterations == 0:
            # Refit and return results
            return self._run_fit_iteration(obj_flux, obj_ferr, self.spectrum_start,
                                           self.spectrum_end, self.base_velocity, templates,
                                           templates_rfft, self.guess_kin,
                                           fix_kinematics=self.fix_kinematics,
                                           obj_to_fit=obj_to_fit, tpl_to_use=tpl_to_use,
                                           weight_errors=True, plot=plot)

        #---------------------------------------------------------------
        # Iteratively filter, fit and reject outliers
        print('Filtering. Boxcar size is {0}'.format(self.filter_boxcar))
        bf = BoxcarFilter(self.filter_boxcar)
        for i in range(self.filter_iterations):

            # Get the filtered object spectra and errors
            sm_obj_flux = bf.smooth(obj_flux)
            if self.filter_operation == 'divide':
                obj_flux_filt = numpy.ma.divide(obj_flux, sm_obj_flux)
                obj_ferr_filt = numpy.ma.absolute(numpy.ma.divide(obj_ferr, sm_obj_flux))
                obj_ferr_filt[numpy.ma.getmaskarray(obj_ferr_filt)] = 1.0
            else:
                obj_flux_filt = obj_flux - sm_obj_flux
                obj_ferr_filt = obj_ferr.copy()
                obj_ferr_filt[numpy.ma.getmaskarray(obj_ferr_filt)] = 1.0

#            pyplot.imshow(obj_flux, origin='lower', interpolation='nearest', aspect='auto')
#            pyplot.show()
#            pyplot.imshow(numpy.ma.log10(obj_flux_filt), origin='lower', interpolation='nearest',
#                          aspect='auto')
#            pyplot.show()
        
            # Get the filtered template spectra
            tpl_flux_filt, tpl_to_use_filt \
                        = self._matched_mask_filter(bf, numpy.ma.getmaskarray(obj_flux_filt),
                                                    templates, templates_rfft, tpl_to_use, result)
            if not self.quiet:
                log_output(self.loggers, 1, logging.INFO, 'Constructing FFT of filtered templates')
            tpl_flux_filt_rfft = numpy.fft.rfft(tpl_flux_filt.data, self.tpl_npad, axis=1)

#            pyplot.imshow(templates, origin='lower', interpolation='nearest', aspect='auto')
#            pyplot.show()
#            pyplot.imshow(numpy.ma.log10(tpl_flux_filt), origin='lower', interpolation='nearest',
#                          aspect='auto')
#            pyplot.show()
       
            # Fit the filtered spectra        
            if numpy.sum(numpy.ma.getmaskarray(tpl_flux_filt)) > 0:
                warnings.warn('There are masked template pixels!')
            result = self._run_fit_iteration(obj_flux_filt, obj_ferr_filt, self.spectrum_start,
                                             self.spectrum_end, self.base_velocity,
                                             tpl_flux_filt.data, tpl_flux_filt_rfft, self.guess_kin,
                                             fix_kinematics=self.fix_kinematics,
                                             obj_to_fit=obj_to_fit, tpl_to_use=tpl_to_use_filt,
                                             degree=self.filt_degree, mdegree=self.filt_mdegree,
                                             dof=self.filt_dof, plot=plot,
                                             weight_errors=(i==self.filter_iterations-1))

            if i == self.filter_iterations - 1:
                break

            # Reject outliers
            # TODO: This will cause an error if boxcar is None!
            obj_flux_filt = PPXFFit.reject_model_outliers(obj_flux_filt, result, rescale=False,
                                                          local_sigma=True,
                                                          boxcar=self.reject_boxcar,
                                                          loggers=self.loggers, quiet=self.quiet)
            # Copy the new mask to the unfiltered spectra
            obj_flux[numpy.ma.getmaskarray(obj_flux_filt)] = numpy.ma.masked
            obj_ferr[numpy.ma.getmaskarray(obj_flux)] = numpy.ma.masked

        # Reset template flags back to original, add the smoothed object
        # spectra into the best-fitting models, and recalculate chi-sqr
        for i in range(self.nobj):
            if result[i] is None or result[i].fit_failed():
                continue
            result[i].tpl_to_use = tpl_to_use[i,:]
            if self.filter_operation == 'divide':
                result[i].bestfit *= sm_obj_flux[i,result[i].start:result[i].end]
            else:
                result[i].bestfit += sm_obj_flux[i,result[i].start:result[i].end]
            chi = (obj_flux[i,result[i].start:result[i].end] - result[i].bestfit)[result[i].gpm] \
                        / obj_ferr[i,result[i].start:result[i].end][result[i].gpm]
#            result[i].robust_rchi2 = numpy.sum(numpy.square(chi))/(chi.size - result[i].kin.size -
#                                                                   max(self.filt_mdegree, 0))
#            pyplot.plot(self.obj_wave[result[i].start:result[i].end],
#                        self.obj_flux[i,result[i].start:result[i].end], color='k', lw=0.5)
#            pyplot.plot(self.obj_wave[result[i].start:result[i].end],
#                        result[i].bestfit, color='C3', lw=1.0)
#            pyplot.show()

        # Return the filtered results
        return result


    def _fit_dispersion_correction(self, templates, templates_rfft, result,
                                   baseline_dispersion=None):
        """
        Calculate the dispersion correction:
          - Construct the optimized, redshifted template *without* the
            convolution with the best-fitting LOSVD.
          - Convolve it with the resolution difference in the data.
          - Use pPXF to fit the matched version to the unmatched version;
            this should be the dispersion correction.

        .. todo::
            Decide how to deal with regions below 2-pixel resolution in
            synthetic spectrum being fit.  How does this alter the
            correction?

        """
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO,
                       'Determining velocity dispersion corrections using fit to '
                       'resolution-matched templates.')

        # Construct the model spectra with only the resolution
        # difference
        model_wlosvd = numpy.ma.empty(self.obj_flux.shape, dtype=float)
        model_wlosvd_msres = numpy.ma.empty(self.obj_flux.shape, dtype=float)
        model_template = numpy.empty((self.nobj, templates.shape[1]), dtype=float)
        model_template_rfft = numpy.empty((self.nobj, templates_rfft.shape[1]), dtype=complex)
        start = numpy.empty(self.nobj, dtype=int)
        end = numpy.empty(self.nobj, dtype=int)
        guess_kin = self.guess_kin.copy()
        model_tpl_to_use = numpy.zeros((self.nobj, self.nobj), dtype=bool)
        res_match_offset = numpy.zeros(self.nobj, dtype=float)

        tpl_ang_per_pix = angstroms_per_pixel(self.tpl_wave, log=True, base=10.)
        unity = numpy.ones(self.tpl_wave.size, dtype=float)

        _obj_to_fit = self.obj_to_fit.copy()

        # Get the nominal kinematics
        nominal_redshift = numpy.zeros(self.nobj, dtype=float)
        nominal_dispersion = numpy.zeros(self.nobj, dtype=float)
        for i in range(self.nobj):
            if _obj_to_fit[i] and (result[i] is None or result[i].fit_failed()):
                _obj_to_fit[i] = False
                continue
            if not _obj_to_fit[i] or result[i] is None or result[i].fit_failed():
                continue
            v = PPXFFit.convert_velocity(result[i].kin[0], 0.0)[0]
            nominal_redshift[i] = v/astropy.constants.c.to('km/s').value
            nominal_dispersion[i] = result[i].kin[1] if baseline_dispersion is None \
                                        else baseline_dispersion

        # Create the spectra to use for the correction
        niter = 1
        for j in range(niter):
            for i in range(self.nobj):
                if not _obj_to_fit[i] or result[i] is None or result[i].fit_failed():
                    continue

                # Get the convolution kernel with 0 velocity shift
                nominal_par, _moments, vj = self._fill_ppxf_par(numpy.array([0.0,
                                                                            nominal_dispersion[i]]))
                nominal_losvd_kernel_rfft = ppxf.losvd_rfft(nominal_par, 1, _moments,
                                                            templates_rfft.shape[1], 1, 0.0,
                                                            self.velscale_ratio, 0.0)[:,0,0]

                # Get the composite template
                model_template[i,:] = numpy.dot(result[i].tplwgt,
                                                templates[result[i].tpl_to_use,:])
                model_template_rfft[i,:] = numpy.dot(result[i].tplwgt,
                                                     templates_rfft[result[i].tpl_to_use,:])

                # Convolve the template to the fitted velocity dispersion
                tmp_wlosvd = numpy.fft.irfft(model_template_rfft[i,:]*nominal_losvd_kernel_rfft,
                                             self.tpl_npad)[:self.npix_tpl]
                model_tpl_to_use[i,i] = True

                # Match the spectral resolution
                tmp_wlosvd_msres, sres, res_match_offset[i], mask, ivar = \
                         match_spectral_resolution(self.tpl_wave, tmp_wlosvd,
                                                   self.tpl_sres.sres(),
                                                   self.obj_wave/(1+nominal_redshift[i]),
                                                   self.obj_sres[i,:], min_sig_pix=0.0, log10=True,
                                                   new_log10=True, quiet=True, no_offset=False)

                # Check 2-pixel resolution limit in resolution-matched
                # spectrum
#                pix_per_fwhm = numpy.ma.divide(numpy.ma.divide(self.tpl_wave, sres),
#                                               tpl_ang_per_pix)
#                pix_per_fwhm = interpolate.interp1d(self.tpl_wave, pix_per_fwhm,
#                                                    fill_value='extrapolate')(self.obj_wave)

                # Resample to match the object spectra
                inRange = self.tpl_wave[[0,-1]] * (1.0 + nominal_redshift[i])
                model_wlosvd[i,:] = Resample(tmp_wlosvd,
                                             x=self.tpl_wave*(1 + nominal_redshift[i]),
                                             inLog=True, newRange=self.obj_wave[[0,-1]],
                                             newpix=self.obj_wave.size, newLog=True).outy
                model_wlosvd_msres[i,:] = Resample(tmp_wlosvd_msres,
                                                   x=self.tpl_wave*(1 + nominal_redshift[i]),
                                                   inLog=True, newRange=self.obj_wave[[0,-1]],
                                                   newpix=self.obj_wave.size, newLog=True).outy

                # Check 2-pixel resolution limit in resampled data
#                _, npix = resample_vector(unity, xRange=inRange, inLog=True,
#                                          newRange=self.obj_wave[[0,-1]],
#                                          newpix=self.obj_wave.size, newLog=True, conserve=True,
#                                          flat=False)
#                new_pix_per_fwhm = numpy.ma.divide(pix_per_fwhm, npix)
#                indx = new_pix_per_fwhm < 2
#                print('Number below 2 pixels: {0}'.format(numpy.sum(indx)))

                # Set the pixels to fit
                start[i] = result[i].start
                end[i] = result[i].end

                model_wlosvd[i,:start[i]] = 0.0
                model_wlosvd[i,:start[i]] = numpy.ma.masked
                model_wlosvd[i,end[i]:] = 0.0
                model_wlosvd[i,end[i]:] = numpy.ma.masked

                model_wlosvd_msres[i,:start[i]] = 0.0
                model_wlosvd_msres[i,:start[i]] = numpy.ma.masked
                model_wlosvd_msres[i,end[i]:] = 0.0
                model_wlosvd_msres[i,end[i]:] = numpy.ma.masked

                # Set the guess kinematics
                guess_kin[i,0] = numpy.log(nominal_redshift[i]+1) \
                                    * astropy.constants.c.to('km/s').value
                guess_kin[i,1] = nominal_dispersion[i]

            # Fit the model to the resolution-matched model
            model_ferr = numpy.ma.ones(self.obj_flux.shape, dtype=float)
            model_ferr[ numpy.ma.getmaskarray(model_wlosvd_msres) ] = numpy.ma.masked
            result_wlosvd = self._run_fit_iteration(model_wlosvd, model_ferr, start, end,
                                                    self.base_velocity, model_template,
                                                    model_template_rfft, guess_kin,
                                                    obj_to_fit=_obj_to_fit,
                                                    tpl_to_use=model_tpl_to_use)#,
                                                    #plot=True)

            result_wlosvd_msres = self._run_fit_iteration(model_wlosvd_msres, model_ferr, start,
                                                          end, self.base_velocity, model_template,
                                                          model_template_rfft, guess_kin,
                                                          obj_to_fit=_obj_to_fit,
                                                          tpl_to_use=model_tpl_to_use)#,
                                                          #plot=True)

            dispersion_correction_err = numpy.array([ rcl is None or rcl.fit_failed()
                                                        or rcl.reached_maxiter()
                                                    or rcls is None or rcls.fit_failed()
                                                        or rcls.reached_maxiter() 
                                        for rcl, rcls in zip(result_wlosvd, result_wlosvd_msres) ])
            dispersion_correction_err &= _obj_to_fit
            indx = numpy.invert(dispersion_correction_err) & _obj_to_fit

            disp_wlosvd = numpy.zeros(self.nobj, dtype=float)
            disp_wlosvd[indx] = numpy.array([ rc.kin[1] for rc in result_wlosvd[indx] ])

            disp_wlosvd_msres = numpy.zeros(self.nobj, dtype=float)
            disp_wlosvd_msres[indx] = numpy.array([ numpy.square(rc.kin[1]) - numpy.square(sigoff) 
                                            for rc, sigoff in zip(result_wlosvd_msres[indx],
                                                                  res_match_offset[indx])])
            disp_wlosvd_msres = numpy.ma.sqrt(disp_wlosvd_msres).filled(0.0)

            dispersion_correction = numpy.zeros(self.nobj, dtype=float)
            dispersion_correction[indx] = numpy.ma.sqrt(numpy.square(disp_wlosvd_msres[indx])
                                            -numpy.square(disp_wlosvd[indx])).filled(0.0) 

            # TODO: The looping and this update should probably be
            # removed
            if j < niter-1:
                nominal_dispersion[indx] = numpy.ma.sqrt( numpy.square(nominal_dispersion[indx])
                                            - numpy.square(dispersion_correction[indx])).filled(0.0)
                if not self.quiet:
                    log_output(self.loggers, 1, logging.INFO,
                               'Dispersion corrections larger than measurement: {0}'.format(
                                        numpy.sum(dispersion_correction > nominal_dispersion)))
                nominal_dispersion[dispersion_correction > nominal_dispersion] \
                        = self.sigma_limits[0]

        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO,
                       'Mean dispersion correction: {0}'.format(
                                            numpy.mean(dispersion_correction[indx])))

        return dispersion_correction, dispersion_correction_err


    def _nominal_dispersion_correction(self, obj_sres, gpm, cz):
        """
        Calculate the dispersion corrections as the quadrature
        difference between the spectral resolution of the template and
        object spectra.  Returns a masked array!
        """
        c = astropy.constants.c.to('km/s').value
        fwhm_inst_obj = c/obj_sres[gpm]
        fwhm_inst_tpl = c/self.tpl_sres(self.obj_wave[gpm]/(1+cz/c))
        mean_fwhm_sqr_diff = numpy.mean(numpy.square(fwhm_inst_obj) - numpy.square(fwhm_inst_tpl))
        if mean_fwhm_sqr_diff < 0:
            return 0., True
        return numpy.sqrt(mean_fwhm_sqr_diff)/DAPConstants.sig2fwhm, False


    def _is_near_bounds(self, result, guess_velocity, tol_frac=1e-2):
        """
        Check if the fitted kinematics are near the imposed limits.
        
        The definition of "near" is that the velocity and higher moments
        cannot be closer than the provided fraction of the total width
        to the boundary.  For the velocity dispersion, the fraction is
        done in log space.
        """

        near_bounds = numpy.zeros(self.moments, dtype=numpy.bool)
        near_lower_sigma_bound = False

        # Velocity
        _velocity_limits = self.velocity_limits + guess_velocity
        Dv = numpy.diff(_velocity_limits)[0]
        tol = Dv*tol_frac

        near_bounds[0] = result.kin[0]-_velocity_limits[0] < tol \
                            or _velocity_limits[1]-result.kin[0] < tol

        # Velocity dispersion
        Ds = numpy.diff(numpy.log10(self.sigma_limits))[0]
        tol = Ds*tol_frac

        near_lower_sigma_bound \
                = numpy.log10(result.kin[1]) - numpy.log10(self.sigma_limits[0]) < tol
        near_bounds[1] = near_lower_sigma_bound \
                        or numpy.log10(self.sigma_limits[1]) - numpy.log10(result.kin[1]) < tol

        if self.moments == 2:
            return near_bounds, near_lower_sigma_bound

        # H3 and H4
        Dh = numpy.diff(self.gh_limits)
        tol = Dh*tol_frac
        near_bounds[2] = result.kin[2] - self.gh_limits[0] < tol \
                                or self.gh_limits[1] - result.kin[2] < tol
        near_bounds[3] = result.kin[3] - self.gh_limits[0] < tol \
                                or self.gh_limits[1] - result.kin[3] < tol
        
        if self.moments == 4:
            return near_bounds, near_lower_sigma_bound

        # H5 and H6
        near_bounds[4] = result.kin[4] - self.gh_limits[0] < tol \
                                or self.gh_limits[1] - result.kin[4] < tol
        near_bounds[5] = result.kin[5] - self.gh_limits[0] < tol \
                                or self.gh_limits[1] - result.kin[5] < tol

        return near_bounds, near_lower_sigma_bound


#    def _set_and_report_failed_status(self, model_mask, model_par, message):
#        if not self.quiet:
#            log_output(self.loggers, 1, logging.INFO, message)
#        return self.bitmask.turn_on(model_mask, 'FIT_FAILED'), \
#               self.bitmask.turn_on(model_par, 'FIT_FAILED')


    def _validate_kinematics(self, model_mask, model_par):
        """
        Validate the returned kinematics.

        Checks:
            - corrected velocity dispersion must be in the range 50-400
              km/s
        """
        sigcor = numpy.square(model_par['KIN'][:,1]) - numpy.square(model_par['SIGMACORR_EMP'])
        indx = ((sigcor < 2500.) | (sigcor > 1.6e5)) & self.obj_to_fit
        if numpy.sum(indx) == 0:
            return
#        model_mask[indx,:] = self.bitmask.turn_on(model_mask[indx,:], 'BAD_SIGMA')
        model_par['MASK'][indx] = self.bitmask.turn_on(model_par['MASK'][indx], 'BAD_SIGMA')


    def _save_results(self, global_fit_result, templates, templates_rfft, result, model_mask,
                      model_par):

        #---------------------------------------------------------------
        # Get the model spectra
        model_flux = PPXFFit.compile_model_flux(self.obj_flux, result)

        # Calculate the model residuals.  They still need to be masked
        # in regions that were *not* included in the fit; see the
        # iterations below
        residual = self.obj_flux - model_flux
        fractional_residual = numpy.ma.divide(self.obj_flux - model_flux, model_flux)
        chi2 = numpy.square(numpy.ma.divide(residual, self.obj_ferr))

        # Instantiate a bad pixel mask to be used
        bpm = numpy.zeros_like(self.obj_wave, dtype=bool)

        # TODO: This is done in initialize_mask
#        # Flag the pixels that were not used
#        model_mask[numpy.ma.getmaskarray(self.obj_flux)] \
#                        = self.bitmask.turn_on(model_mask[numpy.ma.getmaskarray(self.obj_flux)],
#                                               flag='DIDNOTUSE')

        #---------------------------------------------------------------
        # Need to iterate over each spectrum
        for i in range(self.nobj):

            #-----------------------------------------------------------
            # Set output flags
            # No fit was performed
            if result[i] is None:
                model_mask[i,:] = self.bitmask.turn_on(model_mask[i,:], 'NO_FIT')
                model_par['MASK'][i] = self.bitmask.turn_on(model_par['MASK'][i], 'NO_FIT')
                continue

            # No fit attempted because of insufficient data
            if result[i].empty_fit():
                model_mask[i,:] = self.bitmask.turn_on(model_mask[i,:], 'NO_FIT')
                model_par['MASK'][i] = self.bitmask.turn_on(model_par['MASK'][i],
                                                            'INSUFFICIENT_DATA')
                continue

            # Fit attempted but failed
            if result[i].fit_failed():
                model_mask[i,:] = self.bitmask.turn_on(model_mask[i,:], 'FIT_FAILED')
                model_par['MASK'][i] = self.bitmask.turn_on(model_par['MASK'][i], 'FIT_FAILED')

            # Fit successful but hit maximum iterations.
            if result[i].reached_maxiter():
                model_par['MASK'][i] = self.bitmask.turn_on(model_par['MASK'][i], 'MAXITER')

            # Test if the kinematics are near the imposed boundaries.
            # The best fitting model is still provided, but masked.
            near_bound, near_lower_sigma_bound = self._is_near_bounds(result[i],
                                                                      self.guess_kin[i,0])

            # If the velocity dispersion has hit the lower limit, ONLY
            # flag the value as having a MIN_SIGMA.
            if near_lower_sigma_bound:
                near_bound[1] = False
                model_par['MASK'][i] = self.bitmask.turn_on(model_par['MASK'][i], 'MIN_SIGMA')
            # Otherwise, flag both the model and parameter set
            if numpy.any(near_bound):
                if not self.quiet:
                    warnings.warn('Returned parameters for spectrum {0} too close to '
                                  'bounds.'.format(i+1))
                model_par['MASK'][i] = self.bitmask.turn_on(model_par['MASK'][i], 'NEAR_BOUND')
#                model_mask[i,:] = self.bitmask.turn_on(model_mask[i,:], 'NEAR_BOUND')

            # Mask rejected pixels
            s = self.spectrum_start[i]
            e = self.spectrum_end[i]
            original_gpm = numpy.where(numpy.invert(self.obj_flux.mask[i,s:e]))[0]
            rejected_pixels = list(set(original_gpm) - set(result[i].gpm))
            if len(rejected_pixels) > 0:
                model_mask[i,s:e][rejected_pixels] \
                        = self.bitmask.turn_on(model_mask[i,s:e][rejected_pixels],
                                               flag=PPXFFit.rej_flag)

            #-----------------------------------------------------------
            # Save the model parameters and figures of merit

            # Number of fitted pixels
            model_par['NPIXFIT'][i] = len(result[i].gpm)
            # Template weights and errors
            if self._mode_uses_global_template():
                model_par['TPLWGT'][i] = result[i].tplwgt[0]*global_fit_result.tplwgt
                if result[i].tplwgterr is not None:
                    model_par['TPLWGTERR'][i] = result[i].tplwgterr[0]*global_fit_result.tplwgt
            else:
                model_par['TPLWGT'][i][result[i].tpl_to_use] = result[i].tplwgt
                if result[i].tplwgterr is not None:
                    model_par['TPLWGTERR'][i][result[i].tpl_to_use] = result[i].tplwgterr
            # Templates used
            model_par['USETPL'][i] = result[i].tpl_to_use
            # Additive polynomial coefficients
            # TODO: Save nominal polynomial errors as well?
            if self.degree >= 0 and result[i].addcoef is not None:
                model_par['ADDCOEF'][i] = result[i].addcoef
            if self.mdegree > 0 and result[i].multcoef is not None:
                model_par['MULTCOEF'][i] = result[i].multcoef
            # Input kinematics
            model_par['KININP'][i] = self.guess_kin[i,:]
            # Best-fit kinematics
            model_par['KIN'][i] = result[i].kin
            # Kinematic errors
            model_par['KINERR'][i] = result[i].kinerr

            # Ensure that only pixels used in the fit are included in
            # the fit metrics
            chi2[i,model_mask[i,:]>0] = numpy.ma.masked
            residual[i,model_mask[i,:]>0] = numpy.ma.masked
            fractional_residual[i,model_mask[i,:]>0] = numpy.ma.masked
#            bpm[:] = True
#            bpm[result[i].start+result[i].gpm] = False
#            chi2[i,bpm] = numpy.ma.masked
#            residual[i,bpm] = numpy.ma.masked
#            fractional_residual[i,bpm] = numpy.ma.masked

#            # Get the chi-square and rms metrics
#            model_par['CHI2'][i] = 0.0 if numpy.all(residual.mask[i] | self.obj_ferr.mask[i]) \
#                                    else numpy.sum(numpy.square(residual[i]/self.obj_ferr[i]))
#            model_par['RMS'][i] = 0.0 if numpy.all(residual.mask[i]) \
#                                    else numpy.sqrt(numpy.ma.mean(numpy.square(residual[i])))
#            model_par['FRMS'][i] = 0.0 if numpy.all(fractional_residual.mask[i]) \
#                            else numpy.sqrt(numpy.ma.mean(numpy.square(fractional_residual[i])))
#
#            # Chi-square
#            model_par['RCHI2'][i] = model_par['CHI2'][i] \
#                                        / (model_par['NPIXFIT'][i] 
#                                            - self.dof - numpy.sum(model_par['TPLWGT'][i] > 0))
#            model_par['ROBUST_RCHI2'][i] = result[i].robust_rchi2

            # Get growth statistics for the three figures of merit
            model_par['CHIGRW'][i] \
                        = sample_growth(numpy.ma.sqrt(chi2[i,:]).compressed(),
                                        [0.0, 0.68, 0.95, 0.99, 1.0], use_interpolate=False)
            model_par['RMSGRW'][i] \
                        = sample_growth(numpy.ma.absolute(residual[i,:]).compressed(),
                                        [0.0, 0.68, 0.95, 0.99, 1.0], use_interpolate=False)
            model_par['FRMSGRW'][i] \
                        = sample_growth(numpy.ma.absolute(fractional_residual[i,:]).compressed(),
                                        [0.0, 0.68, 0.95, 0.99, 1.0], use_interpolate=False)

            # Calculate the dispersion correction if necessary
            if not self.matched_resolution:
                model_par['SIGMACORR_SRES'][i], err \
                        = self._nominal_dispersion_correction(self.obj_sres[i], result[i].gpm,
                                        PPXFFit.convert_velocity(model_par['KIN'][i,0], 0)[0])
                if err:
                    model_par['MASK'][i] = self.bitmask.turn_on(model_par['MASK'][i],
                                                                'BAD_SIGMACORR_SRES')

        # Get the figures-of-merit for each spectrum
        model_par['CHI2'] = numpy.ma.sum(chi2, axis=1).filled(0.0)
        model_par['RMS'] = numpy.sqrt(numpy.ma.mean(numpy.square(residual), axis=1).filled(0.0))
        model_par['FRMS'] = numpy.sqrt(numpy.ma.mean(numpy.square(fractional_residual),
                                                     axis=1).filled(0.0))
        dof = numpy.sum(model_par['TPLWGT'], axis=1) + self.dof 
        model_par['RCHI2'] = model_par['CHI2'] / (model_par['NPIXFIT'] - dof)

        #---------------------------------------------------------------
        # Calculate the dispersion corrections, if necessary
        # TODO: Get the velocity dispersion corrections here and above
        # regardless of the resolution matching?
        if not self.matched_resolution:
            model_par['SIGMACORR_EMP'], err \
                    = self._fit_dispersion_correction(templates, templates_rfft, result,
                                                      baseline_dispersion=100)
            if numpy.sum(err) > 0:
                model_par['MASK'][err] = self.bitmask.turn_on(model_par['MASK'][err],
                                                              'BAD_SIGMACORR_EMP')

        #---------------------------------------------------------------
        # Test if kinematics are reliable
        # TODO: Skipped for now because unreliable flag not well defined
#        self._validate_kinematics(model_mask, model_par)

        #---------------------------------------------------------------
        # Convert the velocities from pixel units to redshift
        model_par['KININP'][:,0], _ = PPXFFit.convert_velocity(model_par['KININP'][:,0],
                                                               numpy.zeros(self.nobj))
        model_par['KIN'][:,0], model_par['KINERR'][:,0] \
                = PPXFFit.convert_velocity(model_par['KIN'][:,0], model_par['KINERR'][:,0])

#        indx = numpy.invert(self.bitmask.flagged(model_par['MASK'],
#                            flag=['NO_FIT', 'INSUFFICIENT_DATA', 'FIT_FAILED']))
#        pyplot.scatter(model_par['SIGMACORR_SRES'][indx], model_par['SIGMACORR_EMP'][indx],
#                       marker='.', s=50)
#        pyplot.show()

        #---------------------------------------------------------------
        return model_flux, model_mask, model_par


    def fit_SpatiallyBinnedSpectra(self, binned_spectra, par=None, loggers=None, quiet=False,
                                   debug=False):
        """

        This is a basic interface that is geared for the DAP that
        interacts with the rest of the, more general, parts of the
        class.

        This should not declare anything to self!

        """
        # Parameter must be provided
        if par is None:
            raise ValueError('Required parameters for PPXFFit have not been defined.')
        # Check the parameters
        _def = PPXFFitPar._keyword_defaults()
#        if par['regul'] != _def['regul'] or par['reddening'] != _def['reddening'] \
#                or par['component'] != _def['component'] or par['reg_dim'] != _def['reg_dim']:
#            raise NotImplementedError('Cannot use regul, reddening, component, or regul_dim yet.')

        # SpatiallyBinnedSpectra object always needed
        if binned_spectra is None or not isinstance(binned_spectra, SpatiallyBinnedSpectra):
            raise TypeError('Must provide a SpatiallyBinnedSpectra object for fitting.')
        if binned_spectra.hdu is None:
            raise ValueError('Provided SpatiallyBinnedSpectra object is undefined!')
        # If specified, the binned_spectra should have spectral
        # resolution measurements based on a pre-pixelized Gaussian
        if 'STCKPRE' in binned_spectra.hdu['PRIMARY'].header \
                and not binned_spectra.hdu['PRIMARY'].header['STCKPRE']:
            raise ValueError('PPXFFit expects LSF measurements based on a pre-pixelized Gaussian.')

        # TemplateLibrary object always needed
        if par['template_library'] is None \
                or not isinstance(par['template_library'], TemplateLibrary):
            raise TypeError('Must provide a TemplateLibrary object for fitting.')
        if par['template_library'].hdu is None:
            raise ValueError('Provided TemplateLibrary object is undefined!')

        # Select the spectra that meet the selection criteria
        # TODO: Link this to the StellarContinuumModel._bins_to_fit()
        # function...
        good_spec = binned_spectra.above_snr_limit(par['minimum_snr'], debug=debug)

        # Get the object data
        obj_wave = binned_spectra['WAVE'].data.copy()
        obj_flux = binned_spectra.copy_to_masked_array(flag=binned_spectra.do_not_fit_flags())
        obj_ferr = numpy.ma.power(binned_spectra.copy_to_masked_array(ext='IVAR',
                                                    flag=binned_spectra.do_not_fit_flags()) , -0.5)
        obj_sres = binned_spectra.copy_to_array(ext='SPECRES')

        # Warn the user that only a single spectral resolution is used
        # for the templates
        if not self.quiet:
            warnings.warn('Adopting mean spectral resolution of all templates!')
        tpl_sres = numpy.mean(par['template_library']['SPECRES'].data, axis=0).ravel()

        # Perform the fit
        # TODO: Alias window is never used...
        model_wave, model_flux, model_mask, model_par \
                = self.fit(par['template_library']['WAVE'].data.copy(),
                           par['template_library']['FLUX'].data.copy(),
                           binned_spectra['WAVE'].data.copy(), obj_flux[good_spec,:],
                           obj_ferr[good_spec,:], par['guess_redshift'][good_spec],
                           par['guess_dispersion'][good_spec], iteration_mode=par['iteration_mode'],
                           reject_boxcar=par['reject_boxcar'], filter_boxcar=par['filter_boxcar'],
                           filter_operation=par['filter_operation'],
                           filter_iterations=par['filter_iterations'], ensemble=True,
                           velscale_ratio=par['velscale_ratio'], mask=par['pixelmask'],
                           matched_resolution=par['match_resolution'],
                           tpl_sres=tpl_sres, obj_sres=obj_sres[good_spec,:],
                           waverange=par['pixelmask'].waverange, bias=par['bias'],
                           degree=par['degree'], mdegree=par['mdegree'],
                           filt_degree=par['filt_degree'], filt_mdegree=par['filt_mdegree'],
                           moments=par['moments'], loggers=loggers, quiet=quiet, dvtol=1e-9)
                           #plot=True)

        # DEBUG
        assert numpy.sum(model_wave - binned_spectra['WAVE'].data) == 0, \
                    'Incorrect wavelength range'

        # Save the the bin ID numbers indices based on the spectra
        # selected to be fit
        model_par['BINID'] = binned_spectra['BINS'].data['BINID'][good_spec]
        model_par['BINID_INDEX'] = numpy.arange(binned_spectra.nbins)[good_spec]

        # Only return model and model parameters for the *fitted*
        # spectra
        return model_wave, model_flux, model_mask, model_par

    
    def fit(self, tpl_wave, tpl_flux, obj_wave, obj_flux, obj_ferr, guess_redshift,
            guess_dispersion, iteration_mode='global_template', reject_boxcar=100, 
            filter_boxcar=100, filter_operation='divide', filter_iterations=0, ensemble=True,
            velscale_ratio=None, mask=None, usetpl=None, matched_resolution=True, tpl_sres=None,
            obj_sres=None, waverange=None, bias=None, degree=4, mdegree=0, filt_degree=4,
            filt_mdegree=0, moments=2, loggers=None, quiet=False, max_velocity_range=400.,
            alias_window=None, dvtol=1e-10, plot=False, plot_file_root=None):
        r"""
        Wrapper for pPXF with some additional convenience functions.
        Limited implementation at the moment.

        Args:
            tpl_wave (numpy.ndarray):
                1D vector of template wavelengths at rest in angstroms.
            tpl_flux (numpy.ndarray):
                N-templates x N-wavelengths array of template spectra to
                fit.
            obj_wave (numpy.ndarray):
                1D vector of object wavelengths in angstroms.  Does NOT
                need to be same as the template wavelengths.
            obj_flux (numpy.ndarray):
                N-spec x N-wavelengths array of object spectra to fit.
                Can be a numpy.ma.MaskedArray.
            obj_ferr (numpy.ndarray):
                N-spec x N-wavelengths array with the errors in the
                object spectra.
            guess_redshift (:obj:`float`, numpy.ndarray):
                Single or spectrum-specific redshift used to set the
                initial guess kinematics.
            guess_dispersion (:obj:`float`, numpy.ndarray):
                Single or spectrum-specific velocity dispersion used to
                set the initial guess kinematics.
            iteration_mode (:obj:`str`, optional):
                Iteration sequence to perform.  See
                :func:`iteration_modes`.
            reject_boxcar (:obj:`int`, optional):
                Size of the boxcar to use during the rejection
                iteration.  Default is 100.  If None, rejection uses the
                entire residual spectrum.
            filter_boxcar (:obj:`int`, optional):
                Size of the boxcar to use when filtering the spectra.
                Cannot be None.
            filter_iterations (:obj:`int`, optional):
                Number of filtering iterations for the
                'fit_reject_filter' iteration mode.
            ensemble (:obj:`bool`, optional):
                Treat the list of input spectra as an ensemble.
                Currently, this only affects how the spectra are masked.
                Default is to treat them as an ensemble.  When not
                treated as an ensemble, each spectrum is masked
                individually according to the input wavelength range and
                velocity offsets.  *It does not make sense to set the
                iteration_mode to something that will include a fit to
                the global spectrum if you're not treating the list of
                object spectra as an ensemble.*
            velscale_ratio (:obj:`int`, optional):
                Ratio of velocity scale per pixel in the object spectra
                to that for the template spectra.  Default is that they
                should be identical.
            mask (numpy.ndarray, :class:`mangadap.util.pixelmask.SpectralPixelMask`, optional):
                A baseline pixel mask to use during the fitting.  Other
                pixels may be masked via the convenience functions, but
                these pixels will always be masked.
            matched_resolution (:obj:`bool`, optional):
                Flag that the object and template spectra have identical
                spectral resolution.  Default is True.
            tpl_sres (numpy.ndarray, optional):
                One-dimensional vector with the spectral resolution
                (:math:`R = \lambda/\Delta\lambda`) at each wavelength
                of the template spectra.  Default is the resolution is
                not provided and assumed to be same as the object
                resolution.
            obj_sres (numpy.ndarray, optional):
                One- or Two-dimensional array with the spectral
                resolution (:math:`R = \lambda/\Delta\lambda`) at each
                wavelength for (each of) the object spectra.  Default is
                the resolution is not provided and assumed to be same as
                the template resolution.
            waverange (array-like, optional):
                Lower and upper wavelength limits to *include* in the
                fit.  This can be a two-element vector to apply the same
                limits to all spectra, or a N-spec x 2 array with
                wavelength ranges for each spectrum to be fit.  Default
                is to use as much of the spectrum as possible.
            bias (:obj:`float`, optional):
                From the pPXF documentation: This parameter biases the
                (h3, h4, ...) measurements towards zero (Gaussian LOSVD)
                unless their inclusion significantly decreases the error
                in the fit.  Set this to BIAS=0.0 not to bias the fit:
                the solution (including [V, sigma]) will be noisier in
                that case. The default BIAS should provide acceptable
                results in most cases, but it would be safe to test it
                with Monte Carlo simulations. This keyword precisely
                corresponds to the parameter \lambda in the Cappellari &
                Emsellem (2004) paper. Note that the penalty depends on
                the *relative* change of the fit residuals, so it is
                insensitive to proper scaling of the NOISE vector. A
                nonzero BIAS can be safely used even without a reliable
                NOISE spectrum, or with equal weighting for all pixels.
            degree (:obj:`int`, optional):
                From the pPXF documentation: degree of the *additive*
                Legendre polynomial used to correct the template
                continuum shape during the fit (default: 4).  Set DEGREE
                = -1 not to include any additive polynomial.
            mdegree (:obj:`int`, optional):
                From the pPXF documentation: degree of the
                *multiplicative* Legendre polynomial (with mean of 1)
                used to correct the continuum shape during the fit
                (default: 0). The zero degree multiplicative polynomial
                is always included in the fit as it corresponds to the
                weights assigned to the templates.  Note that the
                computation time is longer with multiplicative
                polynomials than with the same number of additive
                polynomials.

                .. note::
                
                    **IMPORTANT**: Multiplicative polynomials cannot be
                    used when the REDDENING keyword is set.

            filt_degree (:obj:`int`, optional):
                The order of the additive polynomial to use when fitting
                the filtered spectra.
            filt_mdegree (:obj:`int`, optional):
                The order of the multiplicative polynomial to use when
                fitting the filtered spectra.
            moments (:obj:`int`, optional):
                From the pPXF documentation: Order of the Gauss-Hermite
                moments to fit. Set this keyword to 4 to fit [h3, h4]
                and to 6 to fit [h3, h4, h5, h6]. Note that in all cases
                the G-H moments are fitted (non-linearly) *together*
                with [V, sigma].

                    - If MOMENTS=2 or MOMENTS is not set then only [V,
                      sigma] are fitted and the other parameters are
                      returned as zero.
                    - If MOMENTS is negative then the kinematics of the
                      given COMPONENT are kept fixed to the input
                      values.
                    - EXAMPLE: We want to keep fixed component 0, which
                      has an LOSVD described by [V, sigma, h3, h4] and
                      is modelled with 100 spectral templates; At the
                      same time we fit [V, sigma] for COMPONENT=1, which
                      is described by 5 templates (this situation may
                      arise when fitting stellar templates with
                      pre-determined stellar kinematics, while fitting
                      the gas emission).  We should give in input to
                      ppxf() the following parameters: component =
                      [0]*100 + [1]*5   # --> [0, 0, ..., 0, 1, 1, 1, 1,
                      1] moments = [-4, 2] start = [[V, sigma, h3, h4],
                      [V, sigma]]

            loggers (:obj:`list`, optional):
                List of `logging.Logger`_ objects to log progress;
                ignored if quiet=True.  Logging is done using
                :func:`mangadap.util.log.log_output`.  Default is no
                logging.
            quiet (:obj:`bool`, optional):
                Suppress all terminal and logging output.  Default is
                False.
            max_velocity_range (:obj:`float`, optional):
                Maximum range (+/-) expected for the fitted velocities
                in km/s.
            alias_window (:obj:`float`, optional):
                The window to mask to avoid aliasing near the edges of
                the spectral range in km/s.  Default is six times
                *max_velocity_range*.
            dvtol (:obj:`float`, optional):
                The velocity scale of the template spectra and object
                spectrum must be smaller than this tolerance.
            plot (:obj:`bool`, optional):
                Show the automatically generated pPXF fit plots.
                
        Returns:
            tuple: Returns 4 numpy.ndarray objects:

                1. The wavelengths of the best fitting model spectra.
                   Nominally the same as the wavelengths of the input
                   object spectra (*obj_wave*).

                2. The fluxes of the best-fitting model spectra.

                3. A mask for the best-fitting models spectra, following
                   from the internal bitmask.

                4. A record array with the fitted model parameters; see
                   :class:`spectralfitting.StellarKinematicsFit._per_stellar_kinematics_dtype`.

        Raises:
            ValueError:
                Raised if the input arrays are not of the correct shape
                or if the pixel scale of the template and object spectra
                is greater than the specified tolerance.

        """
        #---------------------------------------------------------------
        # Initialize the reporting
        if loggers is not None:
            self.loggers = loggers
        self.quiet = quiet

        #---------------------------------------------------------------
        # Check the input; the following operations are necessarily
        # sequential in some cases because of attributes assigned along
        # the way

        # - Mode
        self._check_mode(iteration_mode, reject_boxcar, filter_boxcar, filter_operation,
                         filter_iterations, mdegree)
        # - Objects
        self.obj_wave, self.obj_flux, self.obj_ferr, self.obj_sres \
                = PPXFFit.check_objects(obj_wave, obj_flux, obj_ferr=obj_ferr, obj_sres=obj_sres)
        self.nobj, self.npix_obj = self.obj_flux.shape
        self.input_obj_mask = numpy.ma.getmaskarray(self.obj_flux).copy()
        self.obj_to_fit = numpy.any(numpy.invert(self.input_obj_mask), axis=1)
#        print('Objects to fit: {0}/{1}'.format(numpy.sum(self.obj_to_fit), self.nobj))

        # - Input pixel scale
        self.velscale, self.velscale_ratio \
                = PPXFFit.check_pixel_scale(tpl_wave, self.obj_wave,
                                            velscale_ratio=velscale_ratio, dvtol=dvtol)
        # - Templates
        self.tpl_wave, self.tpl_flux, self.tpl_sres \
                = PPXFFit.check_templates(tpl_wave, tpl_flux, tpl_sres=tpl_sres,
                                          velscale_ratio=self.velscale_ratio)
        self.ntpl, self.npix_tpl = self.tpl_flux.shape
        self.tpl_npad = 2**int(numpy.ceil(numpy.log2(self.npix_tpl)))
#        self.tpl_npad = fftpack.next_fast_len(self.npix_tpl)
        self.tpl_rfft = numpy.fft.rfft(self.tpl_flux, self.tpl_npad, axis=1)

        # - Template usage
        self.usetpl = PPXFFit.check_template_usage_flags(self.nobj, self.ntpl, usetpl)
        # - Input kinematics
        self.input_cz, self.guess_kin = PPXFFit.check_input_kinematics(self.nobj, guess_redshift,
                                                                       guess_dispersion)
        # - Spectral resolution
        self.matched_resolution = PPXFFit.check_resolution_match(self.tpl_sres, self.obj_sres,
                                                                 matched_resolution)
        # - Selected wavelength range: always has shape (self.nobj,2)
        # TODO: When is this used?
        self.waverange = PPXFFit.set_wavelength_range(self.nobj, self.obj_wave, waverange=waverange)

        #---------------------------------------------------------------
        # Make sure that any mode that fits the global spectrum treats
        # the individual spectra as part of an ensemble.  This step is
        # important for setting the spectral range and the masking that
        # is assumed throughout the rest of the function!
        if not ensemble and self._mode_uses_global_spectrum():
            warnings.warn('When fitting the global spectrum, the spectra MUST be treated as an '
                          'ensemble.')
            _ensemble = True
        else:
            _ensemble = ensemble

        #---------------------------------------------------------------
        # Initialize the mask and the spectral range to fit; this was
        # moved here (7 June 2017) so that the filtering order could be
        # set below
        model_mask, init_pix_err, self.spectrum_start, self.spectrum_end \
                = PPXFFit.initialize_pixels_to_fit(self.tpl_wave, self.obj_wave, self.obj_flux,
                                                   self.obj_ferr, self.velscale,
                                                   velscale_ratio=self.velscale_ratio,
                                                   waverange=waverange, mask=mask,
                                                   bitmask=self.bitmask,
                                                   velocity_offset=self.input_cz,
                                                   max_velocity_range=max_velocity_range,
                                                   alias_window=alias_window, ensemble=_ensemble,
                                                   loggers=self.loggers, quiet=self.quiet)

        #---------------------------------------------------------------
        # Set the basic pPXF parameters

        # - Polynomials
        self.degree = max(degree,-1)
        self.mdegree = max(mdegree,0)
        if self.filter_iterations > 0:
            # If filtering, use the size of the filter and the fitted
            # spectral range to set an effective polynomial order used
            # during the initial fit in _fit_all_spectra()
            filter_order = int(numpy.round(numpy.mean(self.spectrum_end-self.spectrum_start)
                                                / self.filter_boxcar))-1
            if self.filter_operation == 'divide':
                self.degree = -1
                self.mdegree = max(1, filter_order)
            else:
                self.degree = max(0, filter_order)
                self.mdegree = 0
        self.filt_degree = max(filt_degree,-1)
        self.filt_mdegree = max(filt_mdegree,0)

        # - Kinematics
        self.velocity_limits, self.sigma_limits, self.gh_limits \
                    = PPXFFit.losvd_limits(self.velscale)
        self.bias = bias
        self.moments = numpy.absolute(moments)
        self.fix_kinematics = moments < 0

        # - Degrees of freedom
        self.dof = self.moments + max(self.mdegree, 0)
        if self.degree >= 0:
            self.dof += self.degree+1
        self.filt_dof = self.moments + max(self.filt_mdegree, 0)
        if self.filt_degree >= 0:
            self.filt_dof += self.filt_degree+1
        if self.fix_kinematics:
            self.dof -= self.moments
            self.filt_dof -= self.moments

        #---------------------------------------------------------------
        # Report the input checks/results
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, 'Number of templates: {0}'.format(self.ntpl))
            log_output(self.loggers, 1, logging.INFO, 'Number of object spectra: {0}'.format(
                                                                                        self.nobj))
            log_output(self.loggers, 1, logging.INFO, 'Pixel scale: {0} km/s'.format(self.velscale))
            log_output(self.loggers, 1, logging.INFO, 'Pixel scale ratio: {0}'.format(
                                                                            self.velscale_ratio))
            log_output(self.loggers, 1, logging.INFO, 'Dispersion limits: {0} - {1}'.format(
                                                                            *self.sigma_limits))
            log_output(self.loggers, 1, logging.INFO, 'Kinematics fixed: {0}'.format(
                                                                            self.fix_kinematics))
            log_output(self.loggers, 2, logging.INFO, 'Model degrees of freedom: {0}'.format(
                                                                            self.dof+self.ntpl))
            log_output(self.loggers, 1, logging.INFO, 'Iteration mode: {0}'.format(
                                                                            self.iteration_mode))
            if self._mode_uses_filter():
                log_output(self.loggers, 2, logging.INFO,
                           'Model degrees of freedom for filtered spectra: {0}'.format(
                           self.filt_dof+self.ntpl))
                log_output(self.loggers, 2, logging.INFO, 'Filter boxcar size: {0}'.format(
                                                            self.filter_boxcar))
                log_output(self.loggers, 2, logging.INFO, 'Filtering operation: {0}'.format(
                                                            self.filter_operation))
                log_output(self.loggers, 2, logging.INFO,
                           'Effective polynomial order: {0}'.format(
                            self.degree if self.filter_operation == 'subtract' else self.mdegree))
            else:
                log_output(self.loggers, 2, logging.INFO,
                           'Additive polynomial order: {0}'.format(self.degree 
                                                                if self.degree > -1 else 'None'))
                log_output(self.loggers, 2, logging.INFO,
                           'Multiplicative polynomial order: {0}'.format(self.mdegree 
                                                                if self.mdegree > 0 else 'None'))

        #---------------------------------------------------------------
        # Initialize the output data
        model_flux = numpy.zeros(self.obj_flux.shape, dtype=numpy.float)
        model_par = self.init_datatable(self.ntpl,
                                        0 if self._mode_uses_filter() else self.degree+1,
                                        0 if self._mode_uses_filter() else max(self.mdegree,0),
                                        self.moments, self.bitmask.minimum_dtype(),
                                        shape=self.nobj)

        # Set the bins; here the ID and index are identical
        model_par['BINID'] = numpy.arange(self.nobj)
        model_par['BINID_INDEX'] = numpy.arange(self.nobj)

        # Save the pixel statistics
        model_par['BEGPIX'] = self.spectrum_start
        model_par['ENDPIX'] = self.spectrum_end
        model_par['NPIXTOT'] = self.spectrum_end - self.spectrum_start

        #---------------------------------------------------------------
        # Flag any errors in the fitted spectral range initialization
        ended_in_error = numpy.array([e is not None for e in init_pix_err])
        if numpy.any(ended_in_error):
            if not self.quiet:
                warnings.warn('Masking failures in some/all spectra.  Errors are: {0}'.format(
                                numpy.array([(i,e)
                                            for i,e in enumerate(init_pix_err)])[ended_in_error]))
            model_par['MASK'][ended_in_error] \
                    = self.bitmask.turn_on(model_par['MASK'][ended_in_error], 'NO_FIT')
        if numpy.all(ended_in_error):
            return self.obj_wave, model_flux, model_mask, model_par

        #---------------------------------------------------------------
        # Get the input pixel shift between the object and template
        # wavelength vectors; interpretted by pPXF as a base velocity
        # shift between the two
        self.base_velocity = numpy.array([PPXFFit.ppxf_tpl_obj_voff(self.tpl_wave,
                                                            self.obj_wave[s:e], self.velscale,
                                                            velscale_ratio=self.velscale_ratio)
                                                for s,e in zip(self.spectrum_start,
                                                               self.spectrum_end)])

        #---------------------------------------------------------------
        # Fit the global spectrum if requested by the iteration mode
        global_fit_result = None
        if self._mode_uses_global_spectrum():
            if not self.quiet:
                log_output(self.loggers, 1, logging.INFO, 'Fitting global spectrum.')
            global_fit_result = self._fit_global_spectrum(obj_to_include=self.obj_to_fit, plot=plot)
            if global_fit_result.fit_failed():
                # pPXF failed!  
                model_mask[:] = self.bitmask.turn_on(model_mask[:], 'FIT_FAILED')
                model_par['MASK'][:] = self.bitmask.turn_on(model_par['MASK'][:], 'FIT_FAILED')
                return self.obj_wave, model_flux, model_mask, model_par
            if not self.quiet:
                log_output(self.loggers, 1, logging.INFO,
                           'Number of non-zero templates in fit to global spectrum: {0}'.format(
                                numpy.sum(global_fit_result.tplwgt > 0)))

        #---------------------------------------------------------------
        # Initialize the template set according to the iteration mode
        if self._mode_uses_global_template():
            templates = numpy.dot(global_fit_result.tplwgt, self.tpl_flux).reshape(1,-1)
            tpl_to_use = numpy.ones((self.nobj,1), dtype=numpy.bool)
            templates_rfft = numpy.fft.rfft(templates, self.tpl_npad, axis=1)
        elif self._mode_uses_nonzero_templates():
            templates = self.tpl_flux
            tpl_to_use = numpy.array([global_fit_result.tplwgt > 0]*self.nobj) & self.usetpl
            templates_rfft = self.tpl_rfft
        elif self._mode_uses_all_templates():
            templates = self.tpl_flux
            tpl_to_use = self.usetpl
            templates_rfft = self.tpl_rfft
        else:
            raise ValueError('Unrecognized iteration mode: {0}'.format(self.iteration_mode))

#        #---------------------------------------------------------------
#        # Print heading of fitted kinematics for each bin
#        if not self.quiet:
#            log_output(self.loggers, 2, logging.INFO, '{0:>5}'.format('INDX')
#                       + ' {:>7} {:>7}'.format('SWAVE', 'EWAVE')
#                       + (' {:>9}'*self.moments).format(*(['KIN{0}'.format(i+1)
#                                                                for i in range(self.moments)]))
#                       + ' {0:>5} {1:>9} {2:>5} {3:>4}'.format('NZTPL', 'CHI2', 'RCHI2', 'STAT'))

        # Fit all spectra
        result = self._fit_all_spectra(templates, templates_rfft, tpl_to_use, plot=plot,
                                       plot_file_root=plot_file_root)

        #---------------------------------------------------------------
        # Save the results
        model_flux, model_mask, model_par \
                = self._save_results(global_fit_result, templates, templates_rfft, result,
                                     model_mask, model_par)

        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, 'pPXF finished')

        return self.obj_wave, model_flux, model_mask, model_par


    @staticmethod
    def obj_tpl_pixelmatch(velscale, tpl_wave, velscale_ratio=None, dvtol=1e-10):
        """
        Confirm that the pixel scale of the template and object data are
        the same within some tolerance, accounting for and input ratio
        of the two.
        """
        _velscale_ratio = 1 if velscale_ratio is None else velscale_ratio
        return numpy.absolute(velscale - spectrum_velocity_scale(tpl_wave)*_velscale_ratio) < dvtol


    @staticmethod
    def fitting_mask(tpl_wave, obj_wave, velscale, velscale_ratio=None, waverange=None,
                     velocity_offset=None, max_velocity_range=400., alias_window=None,
                     loggers=None, quiet=False):
        """
        Return a list of pixels in the object spectrum to be fit using pPXF.
 
        Be clear between velocity (ppxf) vs. redshift (cz) !

        The limits applied to the fitted pixels are:
    
            - Apply the provided wavelength range limit (*waverange*).
            - pPXF will only allow a fit when the number of template
              pixels is the same as or exceeds the number of pixels in
              the object spectrum.  The first step toward limiting
              template spectra that are too long is to truncate the blue
              and red edges that likely won't be used given the provided
              velocity offsets (*velocity_offset*) and the expected
              velocity range (*max_velocity_range*).
            - Remove leading and trailing pixels that will cause alias
              problems during the convolution with the LOSVD
              (*alias_window*).
    
        Args:
            obj_wave (array): Wavelength vector of the object spectrum to be
                fit.
            tpl_wave (array): Wavelength vector of the template library to
                fit to the object spectrum.
            velscale (float): Velocity scale of the pixel.
            velscale_ratio (int): (**Optional**) Ratio of the object
                velscale to the template velscale.  Default is 1 (i.e.
                the two have the same pixel scale).
            waverange (array): (**Optional**) Lower and upper wavelength
                limits to *include* in the analysis.  The array can
                either define a single wavelength range -- shape is (2,)
                -- or a set of wavelength ranges -- shape is (n,2).
                Default is to apply no wavelength range limitation.
            velocity_offset (array): (**Optional**) Vector with the
                velocity offset (expected or actual) between the
                template and the object spectrum in km/s.  Used to
                estimate which wavelengths can be removed from the
                template.  This can be a single offset or a set of
                offsets.  If both waverange and velocity_offset are 2D
                arrays, the number of wavelength ranges and offsets must
                be the same.  Default is that there is no velocity
                offset.
            max_velocity_range (float): (**Optional**) Maximum range
                (+/-) expected for the fitted velocities in km/s.
                Default is 400 km/s.
            alias_window (float) : (**Optional**) The window to mask to
                avoid aliasing near the edges of the spectral range in
                km/s.  Default is six times *max_velocity_range*.
            loggers (list): (**Optional**) List of `logging.Logger`_
                objects to log progress; ignored if quiet=True.  Logging
                is done using :func:`mangadap.util.log.log_output`.
                Default is no logging.
            quiet (bool): (**Optional**) Suppress all terminal and
                logging output.  Default is False.
    
        Returns:
            numpy.ndarray: Four boolean vectors are returned:

                - flags for pixels to include in the fit
                - flags for pixels that were excluded because they were
                  outside the designated wavelength range
                - flags for pixels that were excluded to ensure the
                  proper length of the template spectrum wrt the object
                  spectrum.
                - flags for pixels that were truncated to avoid
                  convolution aliasing

        Raises:
            ValueError: Raised if (1) no pixels are valid for the fit,
                (2) the template and object do not have overlapping
                spectral regions given the expected velocity offset
                between the two, or (3) the turncation to deal with
                aliasing removes all remaining pixels.

        """
        # Check the input
        _velocity_offset = numpy.zeros(1, dtype=float) if velocity_offset is None \
                        else numpy.atleast_1d(velocity_offset)
        _alias_window = 6*max_velocity_range if alias_window is None else alias_window

        # 1. Apply the wavelength range limit, if provided
        now=len(obj_wave)                               # Number of object wavelengths
        if not quiet:
            log_output(loggers, 1, logging.INFO,
                                            'Original number of object pixels: {0}'.format(now))
        if waverange is not None and len(waverange) == 2:
            fit_indx = numpy.ones(now, dtype=bool)
            if waverange[0] is not None:
                fit_indx &= (obj_wave > waverange[0])
            if waverange[1] is not None:
                fit_indx &= (obj_wave < waverange[1])
            if numpy.sum(fit_indx) == 0:
                raise ValueError('Selected wavelength range for analysis contains no pixels!')
        else:
            fit_indx = numpy.full(now, True, dtype=numpy.bool)
        waverange_mask = numpy.invert(fit_indx)

        # Minimum and maximum redshift about primary offsets
        c = astropy.constants.c.to('km/s').value
        z_min = (numpy.amin(_velocity_offset) - max_velocity_range)/c
        z_max = (numpy.amax(_velocity_offset) + max_velocity_range)/c
        if not quiet:
            log_output(loggers, 1, logging.INFO,
                       'Minimum/Maximum redshift: {0}/{1}'.format(z_min, z_max))

        # 2. If the number of template pixels is not >= number of fitted galaxy pixels,
        #    further limit the blue and red edges of the galaxy spectra.
        #    **This must account for the relative pixel scale.**
        now=numpy.sum(fit_indx)                 # Number of good object pixels
        _velscale_ratio = 1 if velscale_ratio is None else int(velscale_ratio)
        ntw=len(tpl_wave)//_velscale_ratio       # Number of template pixels
        if not quiet:
            log_output(loggers, 1, logging.INFO,
                       'After selecting the wavelength range to analyze: {0}'.format(now))
            log_output(loggers, 1, logging.INFO,
                'Number of template pixels (in units of the galaxy pixel scale): {0}'.format(ntw))
        if ntw < now:
            # Indices of wavelengths redward of the redshifted template
            # TODO: Change this to the rigorous calculation of the pPXF
            # velocity: see 
            indx = obj_wave > tpl_wave[0]*(1. + z_min)

            if numpy.sum(indx) == 0:
                raise ValueError('No overlapping wavelengths between galaxy and template!')
            if not quiet:
                log_output(loggers, 1, logging.INFO,
                           'Initial wavelength of template: {0}'.format(tpl_wave[0]))
                log_output(loggers, 1, logging.INFO,
                           'Initial wavelength of redshifted template: {0}'.format(
                                                                    tpl_wave[0]*(1. + z_min)))
                log_output(loggers, 1, logging.INFO,
                           'Initial wavelength of spectrum: {0}'.format(obj_wave[0]))
                log_output(loggers, 1, logging.INFO,
                           'Pixels blueward of redshifted template: {0}'.format(
                                                                    len(obj_wave)-numpy.sum(indx)))
            # Merge with current index
            fit_indx &= indx
            if numpy.sum(fit_indx) == 0:
                raise ValueError('No overlapping wavelengths between galaxy and template!')
            now_= numpy.sum(fit_indx)
            if not quiet:
                log_output(loggers, 1, logging.INFO,
                           'After merging with the specified wavelength range: {0}'.format(now_))
        else:
            now_ = now
    
        # New number of good object pixels
        if ntw < now_:
            fit_indx[ numpy.where(fit_indx)[0][ntw:] ] = False      # Truncate the red edge as well
            if not quiet:
                log_output(loggers, 1, logging.INFO, 
                           'Limit to at most the number of template pixels: {0} {1}'.format(
                                                                        numpy.sum(fit_indx), ntw))
        npix_mask = numpy.invert(fit_indx) & numpy.invert(waverange_mask)

        # 3. Limit wavelength range to avoid aliasing problems in the template convolution
        nalias = int(numpy.floor(_alias_window/velscale))            # Number of pixels to mask
        # Mask to the range that should be unaffected by alias errors
        waverange_tpl = [ tpl_wave[nalias*_velscale_ratio]*(1+z_max),
                          tpl_wave[(ntw-nalias)*_velscale_ratio-1]*(1+z_min) ]
        if not quiet:
            log_output(loggers, 1, logging.INFO,
                       'Mask to these wavelengths to avoid convolution aliasing: {0} - {1}'.format(
                                                            waverange_tpl[0], waverange_tpl[1]))
        indx = numpy.logical_and(obj_wave > waverange_tpl[0], obj_wave < waverange_tpl[1])

        # Merge with current index
        fit_indx &= indx
        if numpy.sum(fit_indx) == 0:
            raise ValueError('No wavelengths available in this range!')

        if not quiet:
            log_output(loggers, 1, logging.INFO,
                       'Final wavelength range to fit: {0} {1}'.format(obj_wave[fit_indx][0],
                                                                       obj_wave[fit_indx][-1]))
        alias_mask = numpy.invert(fit_indx) & numpy.invert(npix_mask)

        return fit_indx, waverange_mask, npix_mask, alias_mask


    @staticmethod
    def ppxf_tpl_obj_voff(tpl_wave, obj_wave, velscale, velscale_ratio=None):
        """
        Determine the pseudo offset in velocity between the template and
        object spectra, just due to the difference in the starting
        wavelengths.
    
        This calculation is independent of the base of the logarithm used
        the sampling of the spectra.

        Assumes wavelengths are logarithmically binned.
    
        Args:
            tpl_wave (numpy.ndarray): Wavelength vector for the template
                library to fit to the object spectrum.
            obj_wave (numpy.ndarray): Wavelength vector for the object
                spectrum to be fit.
            velscale (float): Velocity step per pixel in km/s for the
                **object** spectrum.
            velscale_ratio (int): (**Optional**) The **integer** ratio
                between the velocity scale of the pixel in the galaxy
                data to that of the template data.  This is used only
                when constructing the template library.  Default is
                None, which is the same as assuming that the velocity
                scales are identical.
    
        Returns:
            float: Velocity offset in km/s between the initial wavelengths
            of the template and object spectra.

        .. todo::
            - Implement a check that calculates the velocity ratio directly?
        """
        dlogl = numpy.log(obj_wave[0])-numpy.log(tpl_wave[0]) if velscale_ratio is None \
                    else numpy.log(obj_wave[0])-numpy.mean(numpy.log(tpl_wave[0:velscale_ratio]))
        return dlogl*velscale / numpy.diff(numpy.log(obj_wave[0:2]))[0]


    @staticmethod
    def check_templates(tpl_wave, tpl_flux, tpl_sres=None, velscale_ratio=None):
        r"""
        Check that the input template data is valid for use with
        pPXFFit.
        """
        if len(tpl_wave.shape) != 1:
            raise ValueError('Input template wavelengths must be a vector; all spectra should '
                             'have the same wavelength solution.')
        if isinstance(tpl_flux, numpy.ma.MaskedArray):
            raise TypeError('Template spectra cannot be masked arrays!')
        if tpl_flux.shape[1] != len(tpl_wave):
            raise ValueError('The template spectra fluxes must have the same length as the '
                             'wavelength array.')
        if tpl_sres is not None and tpl_sres.shape != tpl_wave.shape:
            raise ValueError('Provided template resolution vector does not have the correct shape.')

        # Force the number of template pixels to be an integer number of
        # object pixels
        if velscale_ratio is not None and velscale_ratio > 1:
            npix_tpl = len(tpl_wave) - len(tpl_wave) % velscale_ratio
            _tpl_wave = tpl_wave[:npix_tpl]
            _tpl_flux = tpl_flux[:,:npix_tpl]
            _tpl_sres = None if tpl_sres is None else tpl_sres[:npix_tpl]
        else:
            _tpl_wave = tpl_wave
            _tpl_flux = tpl_flux
            _tpl_sres = tpl_sres

        _tpl_sres = None if _tpl_sres is None \
                            else SpectralResolution(_tpl_wave, _tpl_sres, log10=True)
        # TODO: Allow spectral resolution to be spectrum dependent?
        return _tpl_wave, _tpl_flux, _tpl_sres


    @staticmethod
    def check_objects(obj_wave, obj_flux, obj_ferr=None, obj_sres=None):
        r"""
        Check that the input object data is valid for use with pPXFFit.

        Args:
            obj_wave (numpy.ndarray): 1D vector of object wavelengths in
                angstroms.  Does NOT need to be same as the template
                wavelengths.
            obj_flux (numpy.ndarray): :math:`N_{\rm spec}\times N_{\rm
                pix}` array of object spectra to fit.  Can be a
                numpy.ma.MaskedArray.
            obj_ferr (numpy.ndarray): (**Optional**) :math:`N_{\rm
                spec}\times N_{\rm pix}` array with the errors in the
                object spectra.  Can be a numpy.ma.MaskedArray.
            obj_sres (numpy.ndarray): (**Optional**) 1D or 2D array with
                the spectral resolution (:math:`R =
                \lambda/\Delta\lambda`) at each wavelength for (each of)
                the object spectra.  Default is the resolution is not
                provided and assumed to be same as the template
                resolution.

        Returns:
            :obj:`tuple`: Returns four arrays: the object wavelength
            vector, the object flux array, the object flux error and
            the object spectral resolution. The latter two objects
            can be None if the relevant object is None on input. All
            arrays have the same shape as the input ``obj_flux``. If
            ``obj_sres`` is input as a 1D vector, the spectral
            resolution is repeated for all input flux vectors. All
            input arrays are **copies** of the input.
        """
        if len(obj_wave.shape) != 1:
            raise ValueError('Input object wavelengths must be a vector; all spectra should '
                             'have the same wavelength solution.')
        if obj_flux.shape[1] != len(obj_wave):
            raise ValueError('The object spectra fluxes must have the same length as the '
                             'wavelength array.')
        _obj_flux = obj_flux.copy() if isinstance(obj_flux, numpy.ma.MaskedArray) \
                                else numpy.ma.MaskedArray(obj_flux, copy=True)
        if obj_ferr is not None and obj_ferr.shape != _obj_flux.shape:
            raise ValueError('The shape of any provided error array must match the flux array.')
        _obj_ferr = numpy.ma.MaskedArray(numpy.ones(_obj_flux.shape, dtype=float),
                                         mask=numpy.ma.getmaskarray(_obj_flux).copy()) \
                        if obj_ferr is None else \
                            (obj_ferr.copy() if isinstance(obj_ferr, numpy.ma.MaskedArray) \
                                else numpy.ma.MaskedArray(obj_ferr.copy()))
        if obj_sres is not None and obj_sres.shape != obj_wave.shape \
                and obj_sres.shape != _obj_flux.shape:
            raise ValueError('Provided object resolution vector does not have the correct shape.')
        _obj_sres = None if obj_sres is None else obj_sres.copy()
        if _obj_sres is not None and _obj_sres.shape != _obj_flux.shape:
            _obj_sres = numpy.array([_obj_sres]*_obj_flux.shape[0])
        return obj_wave, _obj_flux, _obj_ferr, _obj_sres


    @staticmethod
    def check_pixel_scale(tpl_wave, obj_wave, velscale_ratio=None, dvtol=1e-10):
        """
        Confirm that the pixel scale of the template and object spectra
        are identical within a certain tolerance, accounting for an
        input pixel-scale ratio.  Returns the velocity scale of the
        object spectra and the velocity scale ratio wrt the template
        spectra.
        """
        # Get the pixel scale
        velscale = spectrum_velocity_scale(obj_wave)
        _velscale_ratio = 1 if velscale_ratio is None else int(velscale_ratio)
        if not PPXFFit.obj_tpl_pixelmatch(velscale, tpl_wave, velscale_ratio=_velscale_ratio,
                                          dvtol=dvtol):
            raise ValueError('Pixel scale of the object and template spectra must be identical.')
        return velscale, _velscale_ratio


    @staticmethod
    def losvd_limits(velscale):
        r"""
        Return the limits on the LOSVD parameters used by pPXF.

            - Velocity limits are :math:`\pm 2000` km/s
            - Velocity-disperison limits are from 1/10 pixels to 1000 km/s
            - Limits of the higher orders moments are from -0.3 to 0.3
        """
        velocity_limits = numpy.array([-2e3, 2e3])
        sigma_limits = numpy.array([0.01*velscale, 1e3])
        gh_limits = numpy.array([-0.3, 0.3])
        return velocity_limits, sigma_limits, gh_limits


    @staticmethod
    def reject_model_outliers(obj_flux, ppxf_result, rescale=False, local_sigma=False, boxcar=None,
                              nsigma=3.0, niter=9, loggers=None, quiet=False):
        if boxcar is None and local_sigma:
            raise ValueError('For local sigma determination, must provide boxcar size.')
        model_flux = PPXFFit.compile_model_flux(obj_flux, ppxf_result, rescale=rescale)

        if local_sigma:
            if not quiet:
                log_output(loggers, 1, logging.INFO,
                           'Rejecting using local sigma (boxcar is {0} pixels).'.format(boxcar))
            residual = obj_flux-model_flux # This should be masked where the data were not fit
            bf = BoxcarFilter(boxcar, lo=nsigma, hi=nsigma, niter=niter, y=residual,
                              local_sigma=True)
            obj_flux[bf.output_mask] = numpy.ma.masked
            return obj_flux

        if not quiet:
            log_output(loggers, 1, logging.INFO, 'Rejecting full spectrum outliers.')
        for i in range(niter):
            residual = obj_flux-model_flux  # This should be masked where the data were not fit
            sigma = numpy.array([numpy.ma.std(residual, axis=1)]*residual.shape[1]).T
            indx = numpy.absolute(residual) > nsigma*sigma
            old_mask = numpy.ma.getmaskarray(obj_flux).copy()
            obj_flux[indx] = numpy.ma.masked
            if numpy.sum(numpy.ma.getmaskarray(obj_flux)) == numpy.sum(old_mask):
                break
        return obj_flux


    @staticmethod
    def compile_model_flux(obj_flux, ppxf_result, rescale=False):
        """
        Return the model flux but pulling the models out of the ppxf
        results.  The model can be rescaled to the data based on the
        fitted pixels if rescale is True.

        The output array is masked in the spectral regions below and
        above the fitted wavelength range; any intervening pixels are
        *not* masked, even if they're not included in the fit.
        """
        model_flux = numpy.ma.MaskedArray(numpy.zeros(obj_flux.shape, dtype=float),
                                          mask=numpy.ones(obj_flux.shape, dtype=bool))
        for i in range(obj_flux.shape[0]):
            if ppxf_result[i] is None or ppxf_result[i].fit_failed():
                continue
            s = ppxf_result[i].start
            e = ppxf_result[i].end
            g = ppxf_result[i].gpm
            scale = optimal_scale(obj_flux[i,s:e][g], ppxf_result[i].bestfit[g]) if rescale else 1.
            model_flux[i,s:e] = scale*ppxf_result[i].bestfit
        return model_flux


    @staticmethod
    def convert_velocity(v, verr):
        r"""
        Convert kinematics from pPXF from pixel shifts to redshifts.
        pPXF determines the velocity offset by making the approximation
        that every pixel (logarithmically binned in wavelength) is a
        constant change in velocity.  Over large velocity shifts, this
        approximation can become poor.  An e-mail summary from Michele
        Cappellari:

        The velocity scale per pixel is input as
        
        .. math::
            
            \delta v = c \delta\ln\lambda = c (\ln\lambda_1 -
            \ln\lambda_0)

        The velocites output by pPXF are:

        .. math::

            V = \delta v N_{\rm shift} = c \ln(\lambda_{N_{\rm
            shift}}/\lambda_0)

        which implies that the relation between PPXF output velocity and
        redshift is

        .. math::
            
            1 + z = exp(V/c),

        which reduces z~vel/c in the low-redshift limit.  This function
        converts the :math:`V` values provided by pPXF to :math:`cz`
        velocities.

        .. note::

            **IMPORTANT**: After this conversion, you must revert the
            velocities back to the "pixel-based velocities" (using
            :func:`_revert_velocity`) before using the velocities to
            reconstruct the pPXF fitted model.
        """
        c=astropy.constants.c.to('km/s').value
        return (numpy.exp(v/c)-1.0)*c, verr*numpy.absolute(numpy.exp(v/c))


    @staticmethod
    def revert_velocity(v, verr):
        """
        Revert the velocity back to the "pixelized" velocity returned by
        pPXF.

        Order matters here.  The error computation is NOT a true error
        propagation; it's just the inverse of the "convert" operation
        """
        c=astropy.constants.c.to('km/s').value
        _v = c*numpy.log(v/c+1.0)
        return _v, verr/numpy.absolute(numpy.exp(_v/c))


    @staticmethod
    def reconstruct_model(tpl_wave, templates, obj_wave, kin, weights, velscale, polyweights=None,
                          mpolyweights=None, start=None, end=None, redshift_only=False,
                          sigma_corr=0.0, velscale_ratio=None, dvtol=1e-10, revert_velocity=True):
        """
        Construct a pPXF model spectrum based on a set of input spectra
        and parameters.

        **This function is outdated!  Use :func:`construct_models` or
        :class:`PPXFModel` instead.**
        """
        # Make sure the pixel scales match
        _velscale_ratio = 1 if velscale_ratio is None else velscale_ratio
        if not PPXFFit.obj_tpl_pixelmatch(velscale, tpl_wave, velscale_ratio=_velscale_ratio,
                                          dvtol=dvtol):
            raise ValueError('Pixel scale of the object and template spectra must be identical.')

        # Moments for each kinematic component
        ncomp = 1
        moments = numpy.atleast_1d(kin.size)
        _moments = numpy.full(ncomp, numpy.absolute(moments), dtype=int) if moments.size == 1 \
                            else numpy.absolute(moments)

        # Get the redshift to apply
        redshift = kin[0]/astropy.constants.c.to('km/s').value

        # Check that the corrected sigma is defined
        corrected_sigma = numpy.square(kin[1]) - numpy.square(sigma_corr)
        if not redshift_only and not corrected_sigma > 0:
            warnings.warn('Corrected sigma is 0 or not defined.  Redshifting template only.')
            _redshift_only = True
        else:
            _redshift_only = redshift_only

        # Start and end pixel in the object spectrum to fit
        _start = 0 if start is None else start
        _end = obj_wave.size if end is None else end

        # Construct the composite template
        composite_template = numpy.dot(weights, templates)

#        pyplot.step(tpl_wave, composite_template, where='mid', linestyle='-', lw=0.5,
#                    color='k')
#        pyplot.show()

        # Construct the output models
        model = numpy.ma.zeros(obj_wave.size, dtype=numpy.float)
        if _redshift_only:
            # Resample the redshifted template to the wavelength grid of
            # the binned spectra
            model = Resample(composite_template, x=tpl_wave*(1.0 + redshift), inLog=True,
                             newRange=obj_wave[[0,-1]], newpix=obj_wave.size, newLog=True).outy
            model[:_start] = 0.0
            model[:_start] = numpy.ma.masked
            model[_end:] = 0.0
            model[_end:] = numpy.ma.masked
        else:
            # Perform the same operations as pPXF v6.0.0

            # Get the offset velocity just due to the difference in the
            # initial wavelength of the template and object data
            vsyst = -PPXFFit.ppxf_tpl_obj_voff(tpl_wave, obj_wave[_start:_end], velscale,
                                               velscale_ratio=_velscale_ratio)
            # Account for a modulus in the number of object pixels in
            # the template spectra
            if _velscale_ratio > 1:
                npix_tpl = composite_template.size - composite_template.size % _velscale_ratio
                _composite_template = composite_template[:npix_tpl].reshape(1,-1)
            else:
                _composite_template = composite_template.reshape(1,-1)
            npix_tpl = _composite_template.shape[1]

            # Get the FFT of the composite template
            npad = 2**int(numpy.ceil(numpy.log2(npix_tpl)))
#            npad = fftpack.next_fast_len(npix_tpl)
            ctmp_rfft = numpy.fft.rfft(_composite_template, npad, axis=1)

            # Construct the LOSVD parameter vector
            par = kin.copy()
            # Convert the velocity to pixel units
            if revert_velocity:
                par[0], verr = PPXFFit.revert_velocity(par[0], 1.0)
            # Convert the velocity dispersion to ignore the
            # resolution difference
            par[1] = numpy.sqrt(numpy.square(par[1]) - numpy.square(sigma_corr))
            # Convert the kinematics to pixel units
            par[0:2] /= velscale

            # Construct the model spectrum
            kern_rfft = ppxf.losvd_rfft(par, 1, _moments, ctmp_rfft.shape[1], 1, vsyst/velscale,
                                        _velscale_ratio, 0.0)
            _model = numpy.fft.irfft(ctmp_rfft[0,:] * kern_rfft[:,0,0])[:npix_tpl]
            if _velscale_ratio > 1:
                _model = numpy.mean(_model.reshape(-1,_velscale_ratio), axis=1)
            
            # Copy the model to the output vector
            model[_start:_end] = _model[:_end-_start]

#            pyplot.plot(tpl_wave[:npix_tpl], _composite_template[0,:])
#            pyplot.plot(obj_wave, model)
#            pyplot.show()

        # Account for the polynomials
        x = numpy.linspace(-1, 1, _end-_start)
        if mpolyweights is not None:
            model[_start:_end] *= numpy.polynomial.legendre.legval(x,numpy.append(1.0,mpolyweights))
        if polyweights is not None:
            model[_start:_end] += numpy.polynomial.legendre.legval(x,polyweights)

        return model


    @staticmethod
    def construct_models(tpl_wave, tpl_flux, obj_wave, obj_flux_shape, model_par, select=None,
                         redshift_only=False, deredshift=False, corrected_dispersion=False,
                         dvtol=1e-10):
        """
        Construct models using the provided set of model parameters.
        This is a wrapper for :class:`PPXFModel`.

        Only the shape of the object data is needed, not the data
        itself.

        Allows for a replacement template library that must have the
        same shape as :attr:`tpl_flux`.

        The input velocities are expected to be cz, not "ppxf"
        (pixelized) velocities.

        If redshift_only is true, the provided dispersion is set to 1e-9
        km/s, which is numerically identical to 0 (i.e., just shifting
        the spectrum) in the tested applications.  However, beware that
        this is a HARDCODED number.
        
        .. warning::
            This will not work if the parameters are the result of a
            filtered fit! (iteration_mode = 'fit_reject_filter')

        To convolve the model to the corrected dispersion, instead of
        the uncorrected dispersion, set corrected_dispersion=True.
        Correction *always* uses SIGMACORR_EMP data.
        """
        if redshift_only and corrected_dispersion:
            raise ValueError('The redshift_only and corrected_dispersion are mutually exclusive.')
        if deredshift:
            raise NotImplementedError('Cannot yet deredshift models.')

        # Check the spectral sampling
        velscale_ratio = int(numpy.around(spectrum_velocity_scale(obj_wave)
                                        / spectrum_velocity_scale(tpl_wave)))
        _velscale, _velscale_ratio \
                = PPXFFit.check_pixel_scale(tpl_wave, obj_wave, velscale_ratio=velscale_ratio,
                                             dvtol=dvtol)
        # Check the input spectra
        obj_flux = numpy.zeros(obj_flux_shape, dtype=float)
        _obj_wave, _obj_flux, _, _ = PPXFFit.check_objects(obj_wave, obj_flux)
        nobj = _obj_flux.shape[0]
        _tpl_wave, _tpl_flux, _ = PPXFFit.check_templates(tpl_wave, tpl_flux,
                                                          velscale_ratio=_velscale_ratio)
        ntpl = _tpl_flux.shape[0]

        # Check the shape of the input model parameter database
        if model_par['BINID'].size != nobj:
            raise ValueError('Incorrect number of model-parameter sets.')
        if model_par['TPLWGT'].shape[1] != ntpl:
            raise ValueError('The number of weights does not match the number of templates.')

        # Get the input pixel shift between the object and template
        # wavelength vectors; interpretted by pPXF as a base velocity
        # shift between the two
        vsyst = numpy.array([ -PPXFFit.ppxf_tpl_obj_voff(_tpl_wave, _obj_wave[s:e], _velscale,
                                                         velscale_ratio=_velscale_ratio)
                                    for s,e in zip(model_par['BEGPIX'], model_par['ENDPIX'])])

        # Get the additive and multiplicative degree of the polynomials
        degree = model_par['ADDCOEF'].shape
        degree = -1 if len(degree) == 1 else degree[1]-1
        mdegree = model_par['MULTCOEF'].shape
        mdegree = 0 if len(mdegree) == 1 else mdegree[1]
        moments = model_par['KIN'].shape[1]

        # Only produce selected models
        skip = numpy.zeros(nobj, dtype=bool) if select is None else numpy.invert(select)

        # Instantiate the output model array
        models = numpy.ma.zeros(_obj_flux.shape, dtype=numpy.float)
        # Initially mask everything
        models[:,:] = numpy.ma.masked

        # Set the kinematics
        kin = model_par['KIN'].copy()
        kin[:,0],_ = PPXFFit.revert_velocity(model_par['KIN'][:,0], model_par['KINERR'][:,0])
        if redshift_only:
            kin[:,1] = 1e-9
        elif corrected_dispersion:
            kin[:,1] = numpy.ma.sqrt(numpy.square(model_par['KIN'][:,1]) 
                                        - numpy.square(model_par['SIGMACORR_EMP'])).filled(1e-9)
#        if deredshift:
#            kin[:,0] = 0.0

        # Construct the model for each (selected) object spectrum
        for i in range(nobj):
            if skip[i]:
                continue

            print('Constructing model for spectrum: {0}/{1}'.format(i+1,nobj), end='\r')

            # This has to be redeclared every iteration because the
            # starting and ending pixels might be different (annoying);
            # as will the velocity offset; this means that the FFT of
            # the templates is recalculated at every step...
            f = PPXFModel(_tpl_flux.T,
                          _obj_flux.data[i,model_par['BEGPIX'][i]:model_par['ENDPIX'][i]],
                          _velscale, velscale_ratio=_velscale_ratio, vsyst=vsyst[i],
                          moments=moments, degree=degree, mdegree=mdegree)

            models[i,model_par['BEGPIX'][i]:model_par['ENDPIX'][i]] \
                        = f(kin[i,:], model_par['TPLWGT'][i,:],
                            addpoly=None if degree < 0 else model_par['ADDCOEF'][i,:],
                            multpoly=None if mdegree < 1 else model_par['MULTCOEF'][i,:])

        print('Constructing model for spectrum: {0}/{0}'.format(nobj))
        return models

