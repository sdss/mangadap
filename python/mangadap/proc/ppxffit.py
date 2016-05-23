# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Implements a wrapper class for pPXF.

*License*:
    Copyright (c) 2015, SDSS-IV/MaNGA Pipeline Group
        Licensed under BSD 3-clause license - see LICENSE.rst

*Source location*:
    $MANGADAP_DIR/python/mangadap/proc/ppxffit.py

*Imports and python version compliance*:
    ::

        from __future__ import division
        from __future__ import print_function
        from __future__ import absolute_import
        from __future__ import unicode_literals

        import sys
        import warnings
        if sys.version > '3':
            long = int

        import numpy
        import scipy.interpolate
        import scipy.signal
        import astropy.constants

        from ..par.parset import ParSet
        from ..util.bitmask import BitMask
        from ..util.fileio import init_record_array
        from ..util.instrument import spectrum_velocity_scale, resample_vector
        from ..contrib.ppxf import ppxf
        from .spatiallybinnedspectra import SpatiallyBinnedSpectra
        from .templatelibrary import TemplateLibrary
        from .pixelmask import PixelMask
        from .spectralfitting import StellarKinematicsFit

*Class usage examples*:
        Add examples

*Revision history*:
    | **26 Apr 2016**: Moved from spectralfitting.py to its own file by K. Westfall (KBW)

.. _astropy.io.fits.hdu.hdulist.HDUList: http://docs.astropy.org/en/v1.0.2/io/fits/api/hdulists.html
.. _glob.glob: https://docs.python.org/3.4/library/glob.html
.. _configparser.ConfigParser: https://docs.python.org/3/library/configparser.html#configparser.ConfigParser


"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import sys
import warnings
if sys.version > '3':
    long = int

import logging

import numpy
import scipy.interpolate
import scipy.signal
import astropy.constants

from ..par.parset import ParSet
from ..util.bitmask import BitMask
from ..util.fileio import init_record_array
from ..util.log import log_output
from ..util.instrument import spectrum_velocity_scale, resample_vector
from ..contrib.ppxf import ppxf
from .spatiallybinnedspectra import SpatiallyBinnedSpectra
from .templatelibrary import TemplateLibrary
from .pixelmask import PixelMask, SpectralPixelMask
from .spectralfitting import StellarKinematicsFit
from .util import residual_growth

from matplotlib import pyplot

__author__ = 'Kyle B. Westfall'
# Add strict versioning
# from distutils.version import StrictVersion


class PPXFFitPar(ParSet):
    def __init__(self, template_library_key, template_library, guess_redshift, guess_dispersion,
                 minimum_snr=None, pixelmask=None, bias=None, clean=False, degree=None,
                 mdegree=None, moments=None, oversample=None, sky=None, regul=None, reddening=None,
                 component=None, reg_dim=None):
    
        arr = [ numpy.ndarray, list ]                   # sky, p0, lam, component
        arr_in_fl = [ numpy.ndarray, list, int, float ] # component, reg_dim
        in_fl = [ int, float ]                          # Reddening, bias, regul, vsyst

        _def = self._keyword_defaults()
        
        moment_opt = [ 2, 4, 6 ]

        pars =     [ 'template_library_key', 'template_library', 'guess_redshift',
                     'guess_dispersion', 'minimum_snr', 'pixelmask', 'bias', 'clean', 'degree',
                     'mdegree', 'moments', 'oversample', 'sky', 'regul', 'reddening', 'component',
                     'reg_dim' ]
        values =   [ template_library_key, template_library, guess_redshift, guess_dispersion,
                     minimum_snr, pixelmask, bias, clean, degree, mdegree, moments, oversample,
                     sky, regul, reddening, component, reg_dim ]
        options =  [ None, None, None, None, None, None, None, None, None, None, moment_opt, None,
                     None, None, None, None, None ]
        defaults = [ None, None, None, None, None, None, _def['bias'], _def['clean'],
                     _def['degree'], _def['mdegree'], _def['moments'], _def['oversample'], None,
                     _def['regul'], _def['reddening'], _def['component'], _def['reg_dim'] ]
        dtypes =   [ str, TemplateLibrary, arr_in_fl, arr_in_fl, in_fl, PixelMask, in_fl, bool, int,
                     int, int, int, arr, in_fl, in_fl, arr_in_fl, arr_in_fl ]

        ParSet.__init__(self, pars, values=values, defaults=defaults, options=options,
                        dtypes=dtypes)
        self._check()


    @staticmethod
    def _keyword_defaults():
        """
        Return the keyword defaults.  Pulled from
        :class:`mangadap.contrib.ppxf.ppxf`.
        """
        return { 'bias':None, 'clean':False, 'degree':4, 'mdegree':0, 'moments':2,
                 'oversample':None, 'regul':0, 'reddening':None, 'component':0, 'reg_dim':None }


    def _check(self):
        """
        Perform some preliminary checks on the values of the parameters.
        """
        if self['reddening'] is not None:
            if self['mdegree'] > 0:
                warnings.warn('Cannot both fit multiplicative polynomial and reddening.' \
                              'Ignoring mdegree.')
                self['mdegree'] = 0
        # Other checks (and the one above) done within pPXF


    def toheader(self, hdr):
        hdr['PPXFTPLK'] = (self['template_library_key'], 'Template library key used with pPXF')
        hdr['PPXFBIAS'] = (str(self['bias']) if self['bias'] is None else self['bias'],
                            'pPXF bias value')
        hdr['PPXFCLN'] = (str(self['clean']), 'pPXF cleaned spectrum while fitting')
        hdr['PPXFMOM'] = (self['moments'], 'Number of fitted LOSVD moments in pPXF')
        hdr['PPXFDEG'] = (self['degree'], 'Order of additive polynomial in pPXF')
        hdr['PPXFMDEG'] = (self['mdegree'], 'Order of multiplicative polynomial in pPXF')
        hdr['PPXFOVER'] = (str(self['oversample']), 'Templates oversampled in pPXF')


    def fromheader(self, hdr):
        self['template_library_key'] = hdr['PPXFTPLK']
        self['bias'] = eval(hdr['PPXFBIAS'])
        self['clean'] = eval(hdr['PPXFCLN'])
        self['moments'] = hdr['PPXFMOM']
        self['degree'] = hdr['PPXFDEG']
        self['mdegree'] = hdr['PPXFMDEG']
        self['oversample'] = eval(hdr['PPXFOVER'])




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
    def __init__(self, bitmask, par=None):
        StellarKinematicsFit.__init__(self, 'ppxf', bitmask, par=par)
        self.snr_flag = 'LOW_SNR'
        self.rng_flag = 'OUTSIDE_RANGE'
        self.tpl_flag = 'TPL_PIXELS'
        self.trunc_flag = 'TRUNCATED'
        self.rej_flag = 'PPXF_REJECT'
        self.loggers = None
        self.quiet = False

        self.velocity_limits = None
        self.sigma_limits = None
        self.gh_limits = None

        self.tpl_wave = None
        self.tpl_flux = None
        self.ntpl = None
        self.obj_wave = None
        self.obj_flux = None
        self.obj_ferr = None
        self.nobj = None
        self.velScale = None
        self.guess_kin = None
        self.dof = None
        self.base_velocity = None
#        self.vsyst = None


    @staticmethod
    def _losvd_limits(velScale):
        r"""
        Return the limits on the LOSVD parameters used by pPXF.

            - Velocity limits are :math:`\pm 2000` km/s
            - Velocity-disperison limits are from 1/10 pixels to 1000 km/s
            - Limits of the higher orders moments are from -0.3 to 0.3
        """
        velocity_limits = numpy.array([-2e3, 2e3])
        sigma_limits = numpy.array([0.1*velScale, 1e3])
        gh_limits = numpy.array([-0.3, 0.3])
        return velocity_limits, sigma_limits, gh_limits


    @staticmethod
    def _fitting_mask(tpl_wave, obj_wave, velScale, velocity_offset, waverange=None,
                      max_velocity_range=400., alias_window=2400., loggers=None, quiet=False):
        """
        Return a list of pixels in the object spectrum to be fit using pPXF.
    
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
            velScale (float): Velocity scale of the pixel.
            velocity_offset (array) : Vector with the velocity offset
                (expected or actual) between the template and the object
                spectrum in km/s.  Used to estimate which wavelengths can be
                removed from the template.
            waverange (array type): (**Optional**) Lower and upper
                wavelength limits to *include* in the analysis.
            max_velocity_range (float): (**Optional**) Maximum range (+/-)
                expected for the fitted velocities in km/s.
            alias_window (float) : (**Optional**) The window to mask to avoid
                aliasing near the edges of the spectral range in km/s.
                Default is six times the default *max_velocity_range*.
    
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

        # 1. Apply the wavelength range limit, if provided
        now=len(obj_wave)                               # Number of object wavelengths
        if not quiet:
            log_output(loggers, 1, logging.INFO,
                                            'Original number of object pixels: {0}'.format(now))
        if waverange is not None and len(waverange) == 2:
            fit_indx = numpy.logical_and(obj_wave > waverange[0], obj_wave < waverange[1])
            if numpy.sum(fit_indx) == 0:
                raise ValueError('Selected wavelength range for analysis contains no pixels!')
        else:
            fit_indx = numpy.full(now, True, dtype=numpy.bool)
        waverange_mask = ~fit_indx

        # Minimum and maximum redshift about primary offsets
        c = astropy.constants.c.to('km/s').value
        z_min = (numpy.amin(velocity_offset) - max_velocity_range)/c
        z_max = (numpy.amax(velocity_offset) + max_velocity_range)/c

        # 2. If the number of template pixels is not >= number of fitted galaxy pixels,
        #    further limit the blue and red edges of the galaxy spectra
        now=numpy.sum(fit_indx)         # Number of good object pixels
        ntw=len(tpl_wave)               # Number of template pixels
        if not quiet:
            log_output(loggers, 1, logging.INFO,
                       'After selecting the wavelength range to analyze: {0}'.format(now))
            log_output(loggers, 1, logging.INFO, 'Number of template pixels: {0}'.format(ntw))
        if ntw < now:
            # Indices of wavelengths redward of the redshifted template
            indx = obj_wave > (tpl_wave[0]*(1. + z_min))
            if numpy.sum(indx) == 0:
                raise ValueError('No overlapping wavelengths between galaxy and template!')
            if not quiet:
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
    
        # New number of good object pixels
        if ntw < now_:
            fit_indx[ numpy.where(fit_indx)[0][ntw:] ] = False      # Truncate the red edge as well
            if not quiet:
                log_output(loggers, 1, logging.INFO, 
                           'Impose same as number of template pixels: {0} {1}'.format(
                                                                        numpy.sum(fit_indx), ntw))
        npix_mask = ~(fit_indx) & ~(waverange_mask)

        # 3. Limit wavelength range to avoid aliasing problems in the template convolution
        nalias = int(numpy.floor(alias_window/velScale))            # Number of pixels to mask
        # Mask to the range that should be unaffected by alias errors
        waverange_tpl = [ tpl_wave[nalias]*(1+z_max), tpl_wave[ntw-nalias-1]*(1+z_min) ]
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
        alias_mask = ~(fit_indx) & ~(npix_mask)

        return fit_indx, waverange_mask, npix_mask, alias_mask


    @staticmethod
    def ppxf_tpl_obj_voff(tpl_wave, obj_wave, velScale):
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
            velScale (float): Velocity step per pixel in km/s common to
                both the template and object spectrum.
    
        Returns:
            float: Velocity offset in km/s between the initial wavelengths
            of the template and object spectra.
        """
        return (numpy.log(obj_wave[0])-numpy.log(tpl_wave[0]))*velScale \
                    / numpy.diff(numpy.log(obj_wave[0:2]))[0]


    def fit_SpatiallyBinnedSpectra(self, binned_spectra, par=None, loggers=None, quiet=False):
        """

        This is a basic interface that is geared for the DAP that
        interacts with the rest of the, more general, parts of the
        class.

        This should not declare anything to self!

        """
        # Assign and check parameters if provided
        if par is None:
            raise ValueError('Required parameters for PPXFFit have not been defined.')
        # Check the parameters
        _def = PPXFFitPar._keyword_defaults()
        if par['regul'] != _def['regul'] or par['reddening'] != _def['reddening'] \
                or par['component'] != _def['component'] or par['reg_dim'] != _def['reg_dim']:
            raise NotImplementedError('Cannot use regul, reddening, component, or regul_dim yet.')

        # SpatiallyBinnedSpectra object always needed
        if binned_spectra is None:
            raise ValueError('Must provide spectra object for fitting.')
        if not isinstance(binned_spectra, SpatiallyBinnedSpectra):
            raise TypeError('Must provide a SpatiallyBinnedSpectra object for fitting.')
        if binned_spectra.hdu is None:
            raise ValueError('Provided SpatiallyBinnedSpectra object is undefined!')

        # Template library is required!
        if par['template_library'] is None \
                or not isinstance(par['template_library'], TemplateLibrary):
            raise TypeError('Must provide a TemplateLibrary object for fitting.')
        if par['template_library'].hdu is None:
            raise ValueError('Provided TemplateLibrary object is undefined!')

        # Get the object spectra
        obj_wave = binned_spectra['WAVE'].data.copy()
        obj_flux = binned_spectra.copy_to_masked_array(flag=binned_spectra.do_not_fit_flags())
        obj_ferr = numpy.ma.power(binned_spectra.copy_to_masked_array(ext='IVAR',
                                                    flag=binned_spectra.do_not_fit_flags()), -0.5)

#        pyplot.plot(obj_wave, obj_flux[0,:])
#        pyplot.plot(obj_wave, obj_ferr[0,:])
#        pyplot.show()

        # Select the spectra that meet the selection criteria
        # TODO: Link this to the StellarContinuumModel._bins_to_fit()
        # function...
        good_spec = (binned_spectra['BINS'].data['SNR'] > par['minimum_snr']) \
                        & ~(numpy.array([ b in binned_spectra.missing_bins \
                                                for b in numpy.arange(binned_spectra.nbins)]))
#        print(good_spec)
#        good_spec[2:] = False
#        print(good_spec)

        # Perform the fit
        model_wave, model_flux, model_mask, model_par \
                = self.fit(par['template_library']['WAVE'].data.copy(),
                           par['template_library']['FLUX'].data.copy(),
                           binned_spectra['WAVE'].data.copy(), obj_flux[good_spec,:],
                           obj_ferr[good_spec,:], par['guess_redshift'][good_spec],
                           par['guess_dispersion'][good_spec], mask=par['pixelmask'],
                           waverange=par['pixelmask'].waverange, bias=par['bias'],
                           clean=par['clean'], degree=par['degree'], mdegree=par['mdegree'],
                           moments=par['moments'], oversample=par['oversample'], loggers=loggers,
                           quiet=quiet)

        # Reshape the data to include space for binned spectra that were
        # not fit
        _model_flux = numpy.zeros(obj_flux.shape, dtype=numpy.float)
        _model_flux[good_spec,:] = model_flux

        _model_mask = numpy.zeros(obj_flux.shape, dtype=self.bitmask.minimum_dtype())

        flux_mask = numpy.ma.getmaskarray(obj_flux)
        _model_mask[flux_mask] = self.bitmask.turn_on(_model_mask[flux_mask], 'DIDNOTUSE')
        bad_spec = (binned_spectra['BINS'].data['SNR'] < par['minimum_snr']) \
                        & ~(numpy.array([ b in binned_spectra.missing_bins \
                                                for b in numpy.arange(binned_spectra.nbins)]))
#        bad_spec = ~good_spec

        _model_mask[bad_spec,:] = self.bitmask.turn_on(_model_mask[bad_spec,:], 'LOW_SNR')

        _model_mask[good_spec,:] = model_mask

        _model_par = init_record_array(obj_flux.shape[0], model_par.dtype)
        _model_par[good_spec] = model_par
        _model_par['BIN_INDEX'] = numpy.arange(obj_flux.shape[0])
        _model_par['MASK'][bad_spec] = self.bitmask.turn_on(_model_par['MASK'][bad_spec], 'NO_FIT')

        return model_wave, _model_flux, _model_mask, _model_par


    def fit(self, tpl_wave, tpl_flux, obj_wave, obj_flux, obj_ferr, guess_redshift,
            guess_dispersion, mask=None, waverange=None, bias=None, clean=False, degree=4,
            mdegree=0, moments=2, oversample=None, loggers=None, quiet=False, dvtol=1e-10):

        """

        Wrapper for pPXF with some additional convenience functions.
        Limited implementation at the moment.

        Args:
            tpl_wave (numpy.ndarray): 1D vector of template wavelengths
                at rest in angstroms.
            tpl_flux (numpy.ndarray): N-templates x N-wavelengths array
                of template spectra to fit.
            obj_wave (numpy.ndarray): 1D vector of object wavelengths in
                angstroms.  Does NOT need to be same as the template
                wavelengths.
            obj_flux (numpy.ndarray): N-spec x N-wavelengths array
                of object spectra to fit.
            obj_ferr (numpy.ndarray): N-spec x N-wavelengths array
                with the errors in the object spectra.
            guess_redshift (float or numpy.ndarray): Single or
                spectrum-specific redshift used to set the initial guess
                kinematics.
            guess_dispersion (float or numpy.ndarray): Single or
                spectrum-specific velocity dispersion used to set the
                initial guess kinematics.
            mask (numpy.ndarray):  A
                baseline pixel mask to use during the fitting  Other
                pixels may be masked via the convenience functions, but
                these pixels will always be masked.
            bias (float): (**Optional**) Defaults to 0.0. From the pPXF
                documentation: This parameter biases the (h3, h4, ...)
                measurements towards zero (Gaussian LOSVD) unless their
                inclusion significantly decreases the error in the fit.
                Set this to BIAS=0.0 not to bias the fit: the solution
                (including [V, sigma]) will be noisier in that case. The
                default BIAS should provide acceptable results in most
                cases, but it would be safe to test it with Monte Carlo
                simulations. This keyword precisely corresponds to the
                parameter \lambda in the Cappellari & Emsellem (2004)
                paper. Note that the penalty depends on the *relative*
                change of the fit residuals, so it is insensitive to
                proper scaling of the NOISE vector. A nonzero BIAS can
                be safely used even without a reliable NOISE spectrum,
                or with equal weighting for all pixels.
            clean (bool): (**Optional**) Default is False. From the pPXF
                documentation: set this keyword to use the iterative
                sigma clipping method described in Section 2.1 of
                Cappellari et al.  (2002, ApJ, 578, 787).  This is
                useful to remove from the fit unmasked bad pixels,
                residual gas emissions or cosmic rays.

                .. note::

                    **IMPORTANT**: This is recommended *only* if a
                    reliable estimate of the NOISE spectrum is
                    available. See also note below for SOL.

            degree (int): (**Optional**) Default is 4.  From the pPXF
                documentation: degree of the *additive* Legendre
                polynomial used to correct the template continuum shape
                during the fit (default: 4).  Set DEGREE = -1 not to
                include any additive polynomial.
            mdegree (int): (**Optional**) Default is 0.  From the pPXF
                documentation: degree of the *multiplicative* Legendre
                polynomial (with mean of 1) used to correct the
                continuum shape during the fit (default: 0). The zero
                degree multiplicative polynomial is always included in
                the fit as it corresponds to the weights assigned to the
                templates.  Note that the computation time is longer
                with multiplicative polynomials than with the same
                number of additive polynomials.

                .. note::
                
                    **IMPORTANT**: Multiplicative polynomials cannot be
                    used when the REDDENING keyword is set.

            moments (int): (**Optional**) Default is 2.  From the pPXF
                documentation: Order of the Gauss-Hermite moments to
                fit. Set this keyword to 4 to fit [h3, h4] and to 6 to
                fit [h3, h4, h5, h6]. Note that in all cases the G-H
                moments are fitted (non-linearly) *together* with [V,
                sigma].

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

            oversample (int): (**Optional**) Default is None (for no
                oversampling).  From the pPXF documentation: Set this
                keyword to oversample the template by a factor 30x
                before convolving it with a well sampled LOSVD. This can
                be useful to extract proper velocities, even when sigma
                < 0.7*velScale and the dispersion information becomes
                totally unreliable due to undersampling.  IMPORTANT: One
                should sample the spectrum more finely, if possible,
                before resorting to the use of this keyword!
            loggers (list): (**Optional**) List of `logging.Logger`_ objects
                to log progress; ignored if quiet=True.  Logging is done
                using :func:`mangadap.util.log.log_output`.  Default is
                no logging.
            quiet (bool): (**Optional**) Suppress all terminal and
                logging output.  Default is False.
            dvtol (float): (**Optional**) The velocity scale of the
                template spectra and object spectrum must be smaller
                than this tolerance.  Default is 1e-10.
                
        Returns:
            numpy.ndarray: Returns 4 objects:

                1. The wavelengths of the best fitting model spectra.
                Nominally the same as the wavelengths of the input
                object spectra (*obj_wave*).

                2. The fluxes of the best-fitting model spectra.

                3. A mask for the best-fitting models spectra, following
                from the internal bitmask.

                4. A record array with the fitted model parameters; see
                :class:`spectralfitting.StellarKinematicsFit._per_stellar_kinematics_dtype`.

        Raises:
            ValueError: Raised if the input arrays are not of the
                correct shape or if the pixel scale of the template and
                object spectra is greater than the specified tolerance.
        """
        # Initialize the reporting
        if loggers is not None:
            self.loggers = loggers
        self.quiet = quiet

        # Check the input
        # - Templates
        if len(tpl_wave.shape) != 1:
            raise ValueError('Input template wavelengths must be a vector; all spectra should '
                             'have the same wavelength solution.')
        self.tpl_wave = tpl_wave
        if isinstance(tpl_flux, numpy.ma.MaskedArray):
            raise TypeError('Template spectra cannot be masked arrays!')
        if tpl_flux.shape[1] != len(self.tpl_wave):
            raise ValueError('The template spectra fluxes must have the same length as the '
                             'wavelength array.')
        self.tpl_flux = tpl_flux
        self.ntpl = self.tpl_flux.shape[0]

        # - Objects
        if len(obj_wave.shape) != 1:
            raise ValueError('Input object wavelengths must be a vector; all spectra should '
                             'have the same wavelength solution.')
        self.obj_wave = obj_wave
        if obj_flux.shape[1] != len(self.obj_wave):
            raise ValueError('The object spectra fluxes must have the same length as the '
                             'wavelength array.')
        self.obj_flux = obj_flux if isinstance(obj_flux, numpy.ma.MaskedArray) \
                                    else numpy.ma.MaskedArray(obj_flux)
        if obj_ferr is not None and obj_ferr.shape != self.obj_flux.shape:
            raise ValueError('The shape of any provided error array must match the flux array.')
        self.obj_ferr = obj_ferr if isinstance(obj_ferr, numpy.ma.MaskedArray) \
                                    else numpy.ma.MaskedArray(obj_ferr)
        self.nobj = self.obj_flux.shape[0]

        # Get the pixel scale
        self.velScale = spectrum_velocity_scale(self.obj_wave)
        if numpy.absolute(self.velScale - spectrum_velocity_scale(self.tpl_wave)) > dvtol:
            raise ValueError('Pixel scale of the object and template spectra must be identical.')
        self.velocity_limits, self.sigma_limits, self.gh_limits \
                    = PPXFFit._losvd_limits(self.velScale)

        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, 'Number of templates: {0}'.format(self.ntpl))
            log_output(self.loggers, 1, logging.INFO, 'Number of object spectra: {0}'.format(
                                                                                        self.nobj))
            log_output(self.loggers, 1, logging.INFO, 'Pixel scale: {0} km/s'.format(self.velScale))
            log_output(self.loggers, 1, logging.INFO, 'Dispersion limits: {0} - {1}'.format(
                                                                            *self.sigma_limits))

        # Get the input guess kinematics
        _guess_redshift = numpy.atleast_1d(guess_redshift)
        if len(_guess_redshift) != 1 and len(_guess_redshift) != self.nobj:
            raise ValueError('Must provide a single redshift or one per object spectrum.')
        if len(_guess_redshift) == 1:
            _guess_redshift = numpy.full(self.nobj, _guess_redshift[0], dtype=numpy.float)
        _guess_dispersion = numpy.atleast_1d(guess_dispersion)
        if len(_guess_dispersion) != 1 and len(_guess_dispersion) != self.nobj:
            raise ValueError('Must provide a single dispersion or one per object spectrum.')
        if len(_guess_dispersion) == 1:
            _guess_dispersion = numpy.full(self.nobj, _guess_dispersion[0], dtype=numpy.float)
        self.guess_kin = numpy.array([_guess_redshift*astropy.constants.c.to('km/s').value,
                                        _guess_dispersion]).T

        # Initialize the output arrays
        #  Model flux:
        model_flux = numpy.zeros(self.obj_flux.shape, dtype=numpy.float)
        #  Model mask:
        model_mask = numpy.zeros(self.obj_flux.shape, dtype=self.bitmask.minimum_dtype())
        if mask is not None:
            # Include the mask, if provided
            if isinstance(mask, numpy.ndarray):
                if mask is not None and mask.shape != self.obj_flux.shape:
                    raise ValueError('Shape of object mask array must match its flux array.')
                model_mask = mask.astype(self.bitmask.minimum_dtype())
            if isinstance(mask, SpectralPixelMask):
                model_mask = mask.bits(self.bitmask, self.obj_wave, nspec=self.nobj,
                                       velocity_offsets=self.guess_kin[:,0])
        #  Record array with model parameters
        model_par = init_record_array(self.nobj,
                        self._per_stellar_kinematics_dtype(self.ntpl, degree+1, max(mdegree,0),
                                                           moments, self.bitmask.minimum_dtype()))
        model_par['BIN_INDEX'] = numpy.arange(self.nobj)

        # TODO: Need parameter keywords for max_velocity_range and
        # alias_window
        try:
            fit_indx, waverange_mask, npix_mask, alias_mask \
                = self._fitting_mask(self.tpl_wave, self.obj_wave, self.velScale,
                                     self.guess_kin[:,0], waverange=waverange,
                                     loggers=self.loggers, quiet=self.quiet)
        except ValueError as e:
            if not self.quiet:
                warnings.warn('No fitting done because of a masking error: {0}'.format(e))
            model_par['MASK'] = self.bitmask.turn_on(model_par['MASK'], 'NO_FIT')
            return obj_wave, model_flux, model_mask, model_par

        # Add to the mask
        model_mask[:,waverange_mask] = self.bitmask.turn_on(model_mask[:,waverange_mask],
                                                            self.rng_flag)
        model_mask[:,npix_mask] = self.bitmask.turn_on(model_mask[:,npix_mask], self.tpl_flag)
        model_mask[:,alias_mask] = self.bitmask.turn_on(model_mask[:,alias_mask], self.trunc_flag)
        # And update the internal mask of the data
        self.obj_flux[model_mask > 0] = numpy.ma.masked
        self.obj_ferr[model_mask > 0] = numpy.ma.masked

        # Determine the starting and ending pixels
        start, end = numpy.where(fit_indx)[0][ [0,-1] ]
        end += 1
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, 'Fitting pixels {0} --> {1}'.format(
                       start, end))
            log_output(self.loggers, 1, logging.INFO,
                       'Corresponds to wavelengths {0} --> {1}'.format(self.obj_wave[start],
                                                                       self.obj_wave[end-1]))

        # Save the pixel statistics
        model_par['BEGPIX'][:] = start
        model_par['ENDPIX'][:] = end
        model_par['NPIXTOT'][:] = end-start

        # Determine the degrees of freedom of the fit, used for the
        # brute force calculation of the reduced chi-square
        self.dof = moments + max(mdegree, 0)
        if degree >= 0:
            self.dof += degree+1
        if not self.quiet:
            log_output(self.loggers, 2, logging.INFO, 'Total degrees of freedom: {0}'.format(
                                                                            self.dof+self.ntpl))

        # Get the input pixel shift between the object and template
        # wavelength vectors; interpretted by pPXF as a base velocity
        # shift between the two
        self.base_velocity = self.ppxf_tpl_obj_voff(self.tpl_wave, self.obj_wave[start:end],
                                                    self.velScale)
        self.guess_kin[:,0] -= self.base_velocity
#        self.vsyst = 0 #numpy.mean(self.guess_kin[:,0])

        if not self.quiet:
            log_output(self.loggers, 2, logging.INFO,
                       'Initial wavelength yields baseline velocity offset of: {0:.2f}'.format(
                                                                               self.base_velocity))
#            log_output(self.loggers, 2, logging.INFO,
#                       'Systemic velocity offset provided to pPXF: {0:.2f}'.format(self.vsyst))

        _oversample = False if oversample is None else oversample

        if not self.quiet:
            log_output(self.loggers, 2, logging.INFO, '{0:>5}'.format('INDX')
                    + (' {:>14}'*moments).format(*(['KIN{0}'.format(i+1) for i in range(moments)]))
                    + ' {0:>5} {1:>9} {2:>9} {3:>4}'.format('NZTPL', 'CHI2', 'RCHI2', 'STAT'))

        # Fit each binned spectrum:
        for i in range(self.nobj):
            if not self.quiet:
                log_output(self.loggers, 1, logging.INFO, 'Fit: {0}/{1}'.format(i+1,self.nobj))

            # Get the list of good pixels
            gpm = numpy.where( ~(self.obj_flux.mask[i,start:end]))[0]

            # Check if there is sufficient data for the fit
            # TODO: Too conservative?
            if len(gpm) < self.dof+self.ntpl:
                if not self.quiet:
                    warnings.warn('Insufficient data points ({0}) to fit spectrum {1}'
                                  '(dof={2}).'.format(_nfit, i+1, self.dof+self.ntpl))
                # Flag that the fit was not performed
                model_mask[i,:] = self.bitmask.turn_on(model_mask[i,:], 'NO_FIT')
                model_par['MASK'][i] = self.bitmask.turn_on(model_par['MASK'][i],
                                                            'INSUFFICIENT_DATA')
                continue
            
#            pyplot.plot(self.obj_wave[start:end], self.obj_flux[i,start:end])
#            pyplot.plot(self.obj_wave[start:end], self.obj_flux[i,start:end])
#            pyplot.show()

#            pyplot.plot(self.obj_wave[start:end],
#                        self.obj_flux[i,start:end]/numpy.ma.mean(self.obj_flux[i,start:end]))
#            pyplot.plot(self.tpl_wave*(1+(self.vsyst-self.base_velocity)
#                            / astropy.constants.c.to('km/s').value),
#                        self.tpl_flux[31,:]/numpy.mean(self.tpl_flux[31,:]))
#            pyplot.plot(self.tpl_wave*(1+self.base_velocity/astropy.constants.c.to('km/s').value),
#                        self.tpl_flux[31,:]/numpy.mean(self.tpl_flux[31,:]))
#            pyplot.show()
#            print(self.guess_kin[i,:])

            # Run ppxf
            ppxf_fit = ppxf(self.tpl_flux.T, self.obj_flux.data[i,start:end],
                            self.obj_ferr.data[i,start:end], self.velScale, self.guess_kin[i,:],
                            goodpixels=gpm, bias=bias, clean=clean, degree=degree,
                            mdegree=mdegree, moments=moments, oversample=_oversample, quiet=True)
#                            vsyst=self.vsyst, quiet=True, plot=True)

            # TODO: Add the status to the output table?
            # Check the status of the fit
            if not ppxf_fit.mp.status > 0:
                # Fit failed so log it and continue
                if not self.quiet:
                    log_output(self.loggers, 1, logging.INFO,
                               'Failed pPXF status for spectrum {0}; nothing saved.'.format(i+1))
                model_mask[i,:] = self.bitmask.turn_on(model_mask[i,:], 'FIT_FAILED')
                model_par['MASK'][i] = self.bitmask.turn_on(model_par['MASK'][i], 'FIT_FAILED')
                continue
            if ppxf_fit.mp.status == 5:
                # Fit reached the maximum number of iterations.  Keep
                # the data but throw a warning and log the status via
                # the bitmask.
                if not self.quiet:
                    warnings.warn('pPXF optimizer reached maximum number of iterations for '
                                  'spectrum {0}.'.format(i+1))
                model_par['MASK'][i] = self.bitmask.turn_on(model_par['MASK'][i], 'MAXITER')

            # Save the result
            model_flux[i,start:end] = ppxf_fit.bestfit
            if len(gpm) != len(ppxf_fit.goodpixels):
                rejected_pixels = list(set(gpm)-set(ppxf_fit.goodpixels))
                if not self.quiet:
                    log_output(self.loggers, 1, logging.INFO,
                               'Rejected {0} pixels during fit.'.format(len(rejected_pixels)))
                model_mask[i,start:end][rejected_pixels] \
                        = self.bitmask.turn_on(model_mask[i,start:end][rejected_pixels],
                                               flag=self.rej_flag)

            model_par['NPIXFIT'][i] = len(ppxf_fit.goodpixels)
            model_par['TPLWGT'][i] = ppxf_fit.weights[0:self.ntpl]   # TODO: Store sparsely?

            if degree >= 0:
                model_par['ADDCOEF'][i] = ppxf_fit.polyweights
            if mdegree > 0:
                model_par['MULTCOEF'][i] = ppxf_fit.mpolyweights

            model_par['KININP'][i] = self.guess_kin[i,:]
            model_par['KININP'][i][0] += self.base_velocity
            model_par['KIN'][i] = ppxf_fit.sol
            model_par['KIN'][i][0] += self.base_velocity
            model_par['KINERR'][i] = ppxf_fit.error

            resid = self.obj_flux[i,start:end] - ppxf_fit.bestfit
            model_par['CHI2'][i] = numpy.sum( numpy.square(
                                        (resid/self.obj_ferr[i,start:end])[ppxf_fit.goodpixels] ))
            model_par['RCHI2'][i] = model_par['CHI2'][i] \
                            / (model_par['NPIXFIT'][i] 
                                - self.dof - numpy.sum(model_par['TPLWGT'][i] > 0))
            model_par['ROBUST_RCHI2'][i] = ppxf_fit.chi2

            model_par['RMS'][i] = numpy.sqrt(numpy.ma.mean( numpy.square(
                                                                resid[ppxf_fit.goodpixels]) ))
            # Get growth statistics for the residuals
            model_par['RESID'][i] = residual_growth(resid[ppxf_fit.goodpixels],
                                                        [0.25, 0.50, 0.75, 0.90, 0.99])

            indx = numpy.absolute(ppxf_fit.bestfit) > 0
            _goodpixels = numpy.intersect1d(numpy.arange(len(resid))[indx], ppxf_fit.goodpixels,
                                                assume_unique=True)
            if len(_goodpixels) > 0:
                frac_resid = resid[_goodpixels]/ppxf_fit.bestfit[_goodpixels]
                model_par['FRAC_RMS'][i] = numpy.sqrt(numpy.ma.mean(numpy.square(frac_resid)))
                if len(_goodpixels) > 1:
                    model_par['FRAC_RESID'][i] = residual_growth(frac_resid,
                                                                 [0.25, 0.50, 0.75, 0.90, 0.99])

#            pyplot.step(obj_wave, flux[i,:], where='mid', linestyle='-', lw=0.5, color='k',
#                        zorder=3)
#            pyplot.plot(obj_wave, model_flux[i,:], linestyle='-', lw=1.5, color='r',
#                        zorder=1, alpha=0.5)
#            pyplot.show()

            if not self.quiet:
                log_output(self.loggers, 2, logging.INFO, '{0:>5d}'.format(i+1)
                            + (' {:>14.7e}'*moments).format(*model_par['KIN'][i]) 
                    + ' {0:>5d} {1:>9.2e} {2:>9.2e} {3:>4d}'.format(
                                                            numpy.sum(model_par['TPLWGT'][i] > 0),
                                                            model_par['CHI2'][i],
                                                            model_par['RCHI2'][i],
                                                            ppxf_fit.mp.status))

        model_par['KIN'][:,0], model_par['KINERR'][:,0] \
                = self._convert_velocity(model_par['KIN'][:,0], model_par['KINERR'][:,0])

        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, 'pPXF finished')

        return self.obj_wave, model_flux, model_mask, model_par


    def iterative_fit(self, tpl_wave, tpl_flux, obj_wave, obj_flux, obj_ferr, guess_redshift,
                      guess_dispersion, mask=None, waverange=None, bias=None, clean=False,
                      degree=4, mdegree=0, moments=2, oversample=None, loggers=None, quiet=False,
                      dvtol=1e-10):
        """

        - Construct a "global" spectrum from all bins/spaxels within 1 Re
                - pPXF fit to determine the optimal template weights
                    - need radii (x,y,pa,ell) and Re
            - keep polynomial degree low; additive only!
            - initial estimate for masking everything
            - only fit first two moments (V, sig)
            - fit all templates
                - test different libraries
            - produce a single "global" template to be used in what follows

        - Start loop over all Voronoi bins:
            - In all fits use single template from "global" fit
            
            1st pPXF fit:
                - Still use initial mask for gas
                - fit only V, sig
                - measure residuals
                - estimate error spectrum from residuals (rms vs. wavelength)

            2nd pPXF fit:
                - refit spectrum with new error spectrum
                    - first compare error spectrum and residual spectrum
                - again only fit first two moments
                - apply cleaning algorithm to include strongly deviant pixels in mask
                    - chi-square clip? first only clip in emission-line
                      regions (and redward of SII - to catch red sky lines?)
                - This will produce (V_Gauss and sig_Gauss) needed e.g. for Jeans models
                
            3rd pPXF fit:
                - apply mask (goodpixels) from previous fit
                - fit four moments
                - include penalty/bias in h3,h4
                - This will produce (V, sig, h3, h4) needed e.g. for Schwarzschild models

            4th pPXF fit---Monte Carlo errors: Repeat ~100 times
                - Reshuffle the residuals in small wavelength intervals
                - pPXF fit with NO penalty (or very small)
        """

        # Initialize the reporting
        if loggers is not None:
            self.loggers = loggers
        self.quiet = quiet

        # Check the input
        # - Templates
        if len(tpl_wave.shape) != 1:
            raise ValueError('Input template wavelengths must be a vector; all spectra should '
                             'have the same wavelength solution.')
        self.tpl_wave = tpl_wave
        if isinstance(tpl_flux, numpy.ma.MaskedArray):
            raise TypeError('Template spectra cannot be masked arrays!')
        if tpl_flux.shape[1] != len(self.tpl_wave):
            raise ValueError('The template spectra fluxes must have the same length as the '
                             'wavelength array.')
        self.tpl_flux = tpl_flux
        self.ntpl = self.tpl_flux.shape[0]

        # - Objects
        if len(obj_wave.shape) != 1:
            raise ValueError('Input object wavelengths must be a vector; all spectra should '
                             'have the same wavelength solution.')
        self.obj_wave = obj_wave
        if obj_flux.shape[1] != len(self.obj_wave):
            raise ValueError('The object spectra fluxes must have the same length as the '
                             'wavelength array.')
        self.obj_flux = obj_flux if isinstance(obj_flux, numpy.ma.MaskedArray) \
                                    else numpy.ma.MaskedArray(obj_flux)
        if obj_ferr is not None and obj_ferr.shape != self.obj_flux.shape:
            raise ValueError('The shape of any provided error array must match the flux array.')
        self.obj_ferr = obj_ferr if isinstance(obj_ferr, numpy.ma.MaskedArray) \
                                    else numpy.ma.MaskedArray(obj_ferr)
        self.nobj = self.obj_flux.shape[0]

        # Get the pixel scale
        self.velScale = spectrum_velocity_scale(self.obj_wave)
        if numpy.absolute(self.velScale - spectrum_velocity_scale(self.tpl_wave)) > dvtol:
            raise ValueError('Pixel scale of the object and template spectra must be identical.')
        self.velocity_limits, self.sigma_limits, self.gh_limits \
                    = PPXFFit._losvd_limits(self.velScale)

        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, 'Number of templates: {0}'.format(self.ntpl))
            log_output(self.loggers, 1, logging.INFO, 'Number of object spectra: {0}'.format(
                                                                                        self.nobj))
            log_output(self.loggers, 1, logging.INFO, 'Pixel scale: {0} km/s'.format(self.velScale))
            log_output(self.loggers, 1, logging.INFO, 'Dispersion limits: {0} - {1}'.format(
                                                                            *self.sigma_limits))

        # Get the input guess kinematics
        _guess_redshift = numpy.atleast_1d(guess_redshift)
        if len(_guess_redshift) != 1 and len(_guess_redshift) != self.nobj:
            raise ValueError('Must provide a single redshift or one per object spectrum.')
        if len(_guess_redshift) == 1:
            _guess_redshift = numpy.full(self.nobj, _guess_redshift[0], dtype=numpy.float)
        _guess_dispersion = numpy.atleast_1d(guess_dispersion)
        if len(_guess_dispersion) != 1 and len(_guess_dispersion) != self.nobj:
            raise ValueError('Must provide a single dispersion or one per object spectrum.')
        if len(_guess_dispersion) == 1:
            _guess_dispersion = numpy.full(self.nobj, _guess_dispersion[0], dtype=numpy.float)
        self.guess_kin = numpy.array([_guess_redshift*astropy.constants.c.to('km/s').value,
                                        _guess_dispersion]).T

        # Initialize the output arrays
        #  Model flux:
        model_flux = numpy.zeros(self.obj_flux.shape, dtype=numpy.float)
        #  Model mask:
        model_mask = numpy.zeros(self.obj_flux.shape, dtype=self.bitmask.minimum_dtype())
        if mask is not None:
            # Include the mask, if provided
            if isinstance(mask, numpy.ndarray):
                if mask is not None and mask.shape != self.obj_flux.shape:
                    raise ValueError('Shape of object mask array must match its flux array.')
                model_mask = mask.astype(self.bitmask.minimum_dtype())
            if isinstance(mask, SpectralPixelMask):
                model_mask = mask.bits(self.bitmask, self.obj_wave, nspec=self.nobj,
                                       velocity_offsets=self.guess_kin[:,0])
        #  Record array with model parameters
        model_par = init_record_array(self.nobj,
                        self._per_stellar_kinematics_dtype(self.ntpl, degree+1, max(mdegree,0),
                                                           moments, self.bitmask.minimum_dtype()))
        model_par['BIN_INDEX'] = numpy.arange(self.nobj)

        # TODO: Need parameter keywords for max_velocity_range and
        # alias_window
        try:
            fit_indx, waverange_mask, npix_mask, alias_mask \
                = self._fitting_mask(self.tpl_wave, self.obj_wave, self.velScale,
                                     self.guess_kin[:,0], waverange=waverange,
                                     loggers=self.loggers, quiet=self.quiet)
        except ValueError as e:
            if not self.quiet:
                warnings.warn('No fitting done because of a masking error: {0}'.format(e))
            model_par['MASK'] = self.bitmask.turn_on(model_par['MASK'], 'NO_FIT')
            return obj_wave, model_flux, model_mask, model_par

        # Add to the mask
        model_mask[:,waverange_mask] = self.bitmask.turn_on(model_mask[:,waverange_mask],
                                                            self.rng_flag)
        model_mask[:,npix_mask] = self.bitmask.turn_on(model_mask[:,npix_mask], self.tpl_flag)
        model_mask[:,alias_mask] = self.bitmask.turn_on(model_mask[:,alias_mask], self.trunc_flag)
        # And update the internal mask of the data
        self.obj_flux[model_mask > 0] = numpy.ma.masked
        self.obj_ferr[model_mask > 0] = numpy.ma.masked

        # Determine the starting and ending pixels
        start, end = numpy.where(fit_indx)[0][ [0,-1] ]
        end += 1
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, 'Fitting pixels {0} --> {1}'.format(
                       start, end))
            log_output(self.loggers, 1, logging.INFO,
                       'Corresponds to wavelengths {0} --> {1}'.format(self.obj_wave[start],
                                                                       self.obj_wave[end-1]))

        # Save the pixel statistics
        model_par['BEGPIX'][:] = start
        model_par['ENDPIX'][:] = end
        model_par['NPIXTOT'][:] = end-start

        # Determine the degrees of freedom of the fit, used for the
        # brute force calculation of the reduced chi-square
        self.dof = moments + max(mdegree, 0)
        if degree >= 0:
            self.dof += degree+1
        if not self.quiet:
            log_output(self.loggers, 2, logging.INFO, 'Total degrees of freedom: {0}'.format(
                                                                            self.dof+self.ntpl))

        # Get the input pixel shift between the object and template
        # wavelength vectors; interpretted by pPXF as a base velocity
        # shift between the two
        self.base_velocity = self.ppxf_tpl_obj_voff(self.tpl_wave, self.obj_wave[start:end],
                                                    self.velScale)
        self.guess_kin[:,0] -= self.base_velocity
#        self.vsyst = 0 #numpy.mean(self.guess_kin[:,0])

        if not self.quiet:
            log_output(self.loggers, 2, logging.INFO,
                       'Initial wavelength yields baseline velocity offset of: {0:.2f}'.format(
                                                                               self.base_velocity))
#            log_output(self.loggers, 2, logging.INFO,
#                       'Systemic velocity offset provided to pPXF: {0:.2f}'.format(self.vsyst))

        _oversample = False if oversample is None else oversample

        if not self.quiet:
            log_output(self.loggers, 2, logging.INFO, '{0:>5}'.format('INDX')
                    + (' {:>14}'*moments).format(*(['KIN{0}'.format(i+1) for i in range(moments)]))
                    + ' {0:>5} {1:>9} {2:>9} {3:>4}'.format('NZTPL', 'CHI2', 'RCHI2', 'STAT'))

        # Fit each binned spectrum:
        for i in range(self.nobj):
            if not self.quiet:
                log_output(self.loggers, 1, logging.INFO, 'Fit: {0}/{1}'.format(i+1,self.nobj))

            # Get the list of good pixels
            gpm = numpy.where( ~(self.obj_flux.mask[i,start:end]))[0]

            # Check if there is sufficient data for the fit
            # TODO: Too conservative?
            if len(gpm) < self.dof+self.ntpl:
                if not self.quiet:
                    warnings.warn('Insufficient data points ({0}) to fit spectrum {1}'
                                  '(dof={2}).'.format(_nfit, i+1, self.dof+self.ntpl))
                # Flag that the fit was not performed
                model_mask[i,:] = self.bitmask.turn_on(model_mask[i,:], 'NO_FIT')
                model_par['MASK'][i] = self.bitmask.turn_on(model_par['MASK'][i],
                                                            'INSUFFICIENT_DATA')
                continue
            
#            pyplot.plot(self.obj_wave[start:end], self.obj_flux[i,start:end])
#            pyplot.plot(self.obj_wave[start:end], self.obj_flux[i,start:end])
#            pyplot.show()

#            pyplot.plot(self.obj_wave[start:end],
#                        self.obj_flux[i,start:end]/numpy.ma.mean(self.obj_flux[i,start:end]))
#            pyplot.plot(self.tpl_wave*(1+(self.vsyst-self.base_velocity)
#                            / astropy.constants.c.to('km/s').value),
#                        self.tpl_flux[31,:]/numpy.mean(self.tpl_flux[31,:]))
#            pyplot.plot(self.tpl_wave*(1+self.base_velocity/astropy.constants.c.to('km/s').value),
#                        self.tpl_flux[31,:]/numpy.mean(self.tpl_flux[31,:]))
#            pyplot.show()
#            print(self.guess_kin[i,:])

            # Run ppxf
            ppxf_fit = ppxf(self.tpl_flux.T, self.obj_flux.data[i,start:end],
                            self.obj_ferr.data[i,start:end], self.velScale, self.guess_kin[i,:],
                            goodpixels=gpm, bias=bias, clean=clean, degree=degree,
                            mdegree=mdegree, moments=moments, oversample=_oversample, quiet=True)
#                            vsyst=self.vsyst, quiet=True, plot=True)

            # TODO: Add the status to the output table?
            # Check the status of the fit
            if not ppxf_fit.mp.status > 0:
                # Fit failed so log it and continue
                if not self.quiet:
                    log_output(self.loggers, 1, logging.INFO,
                               'Failed pPXF status for spectrum {0}; nothing saved.'.format(i+1))
                model_mask[i,:] = self.bitmask.turn_on(model_mask[i,:], 'FIT_FAILED')
                model_par['MASK'][i] = self.bitmask.turn_on(model_par['MASK'][i], 'FIT_FAILED')
                continue
            if ppxf_fit.mp.status == 5:
                # Fit reached the maximum number of iterations.  Keep
                # the data but throw a warning and log the status via
                # the bitmask.
                if not self.quiet:
                    warnings.warn('pPXF optimizer reached maximum number of iterations for '
                                  'spectrum {0}.'.format(i+1))
                model_par['MASK'][i] = self.bitmask.turn_on(model_par['MASK'][i], 'MAXITER')

            # Save the result
            model_flux[i,start:end] = ppxf_fit.bestfit
            if len(gpm) != len(ppxf_fit.goodpixels):
                rejected_pixels = list(set(gpm)-set(ppxf_fit.goodpixels))
                if not self.quiet:
                    log_output(self.loggers, 1, logging.INFO,
                               'Rejected {0} pixels during fit.'.format(len(rejected_pixels)))
                model_mask[i,start:end][rejected_pixels] \
                        = self.bitmask.turn_on(model_mask[i,start:end][rejected_pixels],
                                               flag=self.rej_flag)

            model_par['NPIXFIT'][i] = len(ppxf_fit.goodpixels)
            model_par['TPLWGT'][i] = ppxf_fit.weights[0:self.ntpl]   # TODO: Store sparsely?

            if degree >= 0:
                model_par['ADDCOEF'][i] = ppxf_fit.polyweights
            if mdegree > 0:
                model_par['MULTCOEF'][i] = ppxf_fit.mpolyweights

            model_par['KININP'][i] = self.guess_kin[i,:]
            model_par['KININP'][i][0] += self.base_velocity
            model_par['KIN'][i] = ppxf_fit.sol
            model_par['KIN'][i][0] += self.base_velocity
            model_par['KINERR'][i] = ppxf_fit.error

            resid = self.obj_flux[i,start:end] - ppxf_fit.bestfit
            model_par['CHI2'][i] = numpy.sum( numpy.square(
                                        (resid/self.obj_ferr[i,start:end])[ppxf_fit.goodpixels] ))
            model_par['RCHI2'][i] = model_par['CHI2'][i] \
                            / (model_par['NPIXFIT'][i] 
                                - self.dof - numpy.sum(model_par['TPLWGT'][i] > 0))
            model_par['ROBUST_RCHI2'][i] = ppxf_fit.chi2

            model_par['RMS'][i] = numpy.sqrt(numpy.ma.mean( numpy.square(
                                                                resid[ppxf_fit.goodpixels]) ))
            # Get growth statistics for the residuals
            model_par['RESID'][i] = residual_growth(resid[ppxf_fit.goodpixels],
                                                        [0.25, 0.50, 0.75, 0.90, 0.99])

            indx = numpy.absolute(ppxf_fit.bestfit) > 0
            _goodpixels = numpy.intersect1d(numpy.arange(len(resid))[indx], ppxf_fit.goodpixels,
                                                assume_unique=True)
            if len(_goodpixels) > 0:
                frac_resid = resid[_goodpixels]/ppxf_fit.bestfit[_goodpixels]
                model_par['FRAC_RMS'][i] = numpy.sqrt(numpy.ma.mean(numpy.square(frac_resid)))
                if len(_goodpixels) > 1:
                    model_par['FRAC_RESID'][i] = residual_growth(frac_resid,
                                                                 [0.25, 0.50, 0.75, 0.90, 0.99])

#            pyplot.step(obj_wave, flux[i,:], where='mid', linestyle='-', lw=0.5, color='k',
#                        zorder=3)
#            pyplot.plot(obj_wave, model_flux[i,:], linestyle='-', lw=1.5, color='r',
#                        zorder=1, alpha=0.5)
#            pyplot.show()

            if not self.quiet:
                log_output(self.loggers, 2, logging.INFO, '{0:>5d}'.format(i+1)
                            + (' {:>14.7e}'*moments).format(*model_par['KIN'][i]) 
                    + ' {0:>5d} {1:>9.2e} {2:>9.2e} {3:>4d}'.format(
                                                            numpy.sum(model_par['TPLWGT'][i] > 0),
                                                            model_par['CHI2'][i],
                                                            model_par['RCHI2'][i],
                                                            ppxf_fit.mp.status))

        model_par['KIN'][:,0], model_par['KINERR'][:,0] \
                = self._convert_velocity(model_par['KIN'][:,0], model_par['KINERR'][:,0])

        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, 'pPXF finished')

        return self.obj_wave, model_flux, model_mask, model_par


    @staticmethod
    def _convert_velocity(v, verr):
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

            V = \delta V N_{\rm shift} = c \ln(\lambda_{N_{\rm
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
    def _revert_velocity(v,verr):
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
    def _losvd_kernel(par, velScale, base_velocity=0.0):#, nwave):
        """
        Return the sampled kernel.

        Meant to match computation of model in pPXF, see
        :func:`mangadap.contrib.ppxf._fitfunc`

        CURRENTLY WILL NOT WORK WITH:
            - multiple kinematic components
            - point-symmetric fits
            - oversampling

        Input expected to be the parameters *after* the pPXF velocities
        had been converted to redshifts.

        """
        _par = par.copy()
        _par[0], verr = PPXFFit._revert_velocity(_par[0], 1.0)
        _par[0] += base_velocity
        #_par[0] = 0.0
        _par[0:2] /= velScale      # Convert to pixels

        dx = int(abs(_par[0])+5.0*par[1])
        nl = 2*dx+1
        x = numpy.linspace(-dx,dx,nl)
        w = (x - _par[0])/_par[1]
        w2 = numpy.square(w)
        gauss = numpy.exp(-0.5*w2)
        losvd = gauss/numpy.sum(gauss)

        # Add the higher order terms
        if len(_par) > 2:
            # h_3 h_4
            poly = 1 + _par[2]/numpy.sqrt(3)*(w*(2*w2-3)) + _par[3]/numpy.sqrt(24)*(w2*(4*w2-12)+3)
            if len(_par) > 4:
                # h_5 h_6
                poly += _par[4]/np.sqrt(60)*(w*(w2*(4*w2-20)+15)) \
                            + _par[5]/np.sqrt(720)*(w2*(w2*(8*w2-60)+90)-15)
            losvd *= poly

        return _par[0]*velScale, losvd

#        npad = int(2**np.ceil(np.log2(nwave + nl/2)))
#        losvd_pad = numpy.zeros(npad)

        

    @staticmethod
    def construct_models(binned_spectra, template_library, model_par, degree, mdegree,
                         redshift_only=False):
        """

        Construct models using the provided set of model parameters and
        template library.  The number of template weights in the
        provided model parameters must match the number of template
        spectra.

        If redshift_only is true, ignore the LOSVD and simply shift the
        model to the correct redshift.

        .. todo::

            This doesn't appear to give exactly the correct
            reconstruction of the pPXF models.

        """
        if not isinstance(binned_spectra, SpatiallyBinnedSpectra):
            raise TypeError('Must provide a SpatiallyBinnedSpectra object.')
        if not isinstance(template_library, TemplateLibrary):
            raise TypeError('Must provide a TemplateLibrary object.')
        if model_par['TPLWGT'].shape[1] != template_library.ntpl:
            raise ValueError('The number of weights does not match the number of templates.')

        velScale = spectrum_velocity_scale(binned_spectra['WAVE'].data)

        tpl_wave = template_library['WAVE'].data
        obj_wave = binned_spectra['WAVE'].data

        outRange = obj_wave[ [0,-1] ]
        outpix = binned_spectra.nwave

        models = numpy.ma.zeros((binned_spectra.nbins, binned_spectra.nwave), dtype=numpy.float)

#        model_valid = numpy.zeros(binned_spectra.nbins, dtype=numpy.bool)

        for i in range(binned_spectra.nbins):
            if i in binned_spectra.missing_bins:
                models[i,:] = numpy.ma.masked
                continue
            if model_par['MASK'][i] > 0:
                models[i,:] = numpy.ma.masked
                warnings.warn('Fit to spectrum {0} masked.'.format(i+1))
                continue

            # Get the composite template
            indx = model_par['TPLWGT'][i] > 0
            composite_template = numpy.dot(model_par['TPLWGT'][i][indx],
                                           template_library['FLUX'].data[indx,:])

#            pyplot.step(tpl_wave, composite_template, where='mid', linestyle='-', lw=0.5,
#                         color='k')
#            pyplot.show()

            redshift = model_par['KIN'][i,0]/astropy.constants.c.to('km/s').value
            start = model_par['BEGPIX'][i]
            end = model_par['ENDPIX'][i]

            # Redshift and broadened using the LOSVD parameters
            if redshift_only:
                # Resample the redshifted template to the wavelength grid of
                # the binned spectra
                inRange = tpl_wave[[0,-1]] * (1.0 + redshift)
                _composite_template = composite_template
            else:
#                print('{0}/{1}'.format(i+1,binned_spectra.nbins))
#                print('Start,End: ', start, end)
#                print(obj_wave[start:end].shape)
#                print(tpl_wave.shape)
                base_velocity = PPXFFit.ppxf_tpl_obj_voff(tpl_wave, obj_wave[start:end], velScale)
                vel, losvd = PPXFFit._losvd_kernel(model_par['KIN'][i,:], velScale,
                                                   base_velocity=base_velocity)
                _composite_template = scipy.signal.fftconvolve(composite_template, losvd,
                                                               mode="same")
                voff, verr = PPXFFit._convert_velocity(-base_velocity, 1.0)
                #inRange = tpl_wave[ [0,-1] ]*(1+redshift)
                inRange = tpl_wave[[0,-1]] * (1.0 + voff/astropy.constants.c.to('km/s').value)
                
            # Resample the template to the galaxy wavelengths
            wave, models[i,:] = resample_vector(_composite_template, xRange=inRange, inLog=True,
                                                newRange=outRange, newpix=outpix, newLog=True,
                                                flat=False)


#            pyplot.plot(tpl_wave*(1+redshift), composite_template, linestyle='-', lw=1.5,
#                        color='r', zorder=1, alpha=0.5)#, where='mid')
#            pyplot.plot(wave, models[i,:], linestyle='-', lw=0.5, color='k', zorder=3)
#                        #, where='mid')
#            pyplot.show()

            x = numpy.linspace(-1, 1, model_par['NPIXTOT'][i])
            if mdegree > 0:
                models[i,start:end] *= numpy.polynomial.legendre.legval(x,
                                                        numpy.append(1.0,model_par['MULTCOEF'][i]))
            if degree > -1:
                model[i,start:end] += numpy.polynomial.legendre.legval(x, model_par['ADDCOEF'][i])

        return models













