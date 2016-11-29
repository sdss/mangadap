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
    | **26 Apr 2016**: Moved from spectralfitting.py to its own file by
        K. Westfall (KBW)
    | **05 Jul 2016**: (KBW) V6.0.0 of pPXF does not use the oversample
        keyword given a better solution; see Cappellari (in prep).  This
        keyword was therefore removed from the parameter set.
    | **06 Jul 2016**: (KBW) Use v6.0.0 pPXF functions to compute models
        using new LOSVD kernel functionality.
    | **10 Oct 2016**: (KBW) Fixed error in calculation of velocity
        offset between template and object spectra to account for
        different size pixels.
    | **31 Oct 2016**: (KBW) Allow the spectral resolution to be a
        vector per spaxel or a vector per set of input spectra in
        :func:`PPXFFit.fit`.
    | **01 Nov 2016**: (KBW) Added new iteration method that does not do
        a first fit to the global spectrum but does include the
        rejection iteration.  Allows users to treat each spectrum
        provided to :func:`PPXFFit.fit` individually.  Fixed goodpixel
        mask when not first fitting the global spectrum.
    | **02 Nov 2016**: (KBW) Added ability to limit which templates are
        fit to each spectrum in :func:`PPXFFit.fit` using usetpl kwarg.

.. todo::
    - Include input to :func:`PPXFFit.fit` that allows the user to
      specify a subset of templates to use for each object spectrum.
    - Additional iteration modes?
        - Allow the fit to continue to iterate until the velocity is not
          up against one of the +/- 2000 km/s limits?  Could be useful
          for poor redshift guesses.
    - Keep track of which templates were non-zero in the fit to the
      global spectrum, if it's fit

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
from ..util.instrument import spectrum_velocity_scale, resample_vector, spectral_resolution
from ..util.constants import constants
from ..contrib.ppxf import ppxf, _templates_rfft, _losvd_rfft
from ..contrib import ppxf_util
from .spatiallybinnedspectra import SpatiallyBinnedSpectra
from .templatelibrary import TemplateLibrary
from .pixelmask import PixelMask, SpectralPixelMask
from .spectralfitting import StellarKinematicsFit
# from .spectralstack import SpectralStack
from .util import residual_growth

from matplotlib import pyplot

__author__ = 'Kyle B. Westfall'
# Add strict versioning
# from distutils.version import StrictVersion


class PPXFFitPar(ParSet):
    r"""

    Define a parameter set used by the pPXF fitting method.

    .. todo::
        The overlap between this and
        :class:`mangadap.proc.stellarcontinuummodel.StellarContinuumModelDef`
        is not well designed.

    Args:
        template_library_key (str): Keyword of the library to fit.  See
            :func:`mangadap.proc.templatelibrary.available_template_libraries`.

        template_library
            (:class:`mangadap.proc.templatelibrary.TemplateLibrary`):
            Object with the spectra in the template library that have
            been prepared for analysis of the data.

        guess_redshift (array-like): Initial guess for the redshift
            (:math:`cz`) of each binned spectrum.

        guess_dispersion (array-like): Initial guess for the velocity
            dispersion for each binned spectrum.

        iteration_mode (str): (**Optional**) Iteration mode to use; see
            :func:`PPXFFit.iteration_modes`.

        match_resolution (bool): (**Optional**) Match the spectral
            resolution of the template to that of the galaxy data.  This
            is used only when constructing the template library.
            Default is True.
        
        velscale_ratio (int): (**Optional**) The **integer** ratio
            between the velocity scale of the pixel in the galaxy data
            to that of the template data.  This is used only when
            constructing the template library.  Default is None, which
            is the same as assuming that the velocity scales are
            identical.
        
        minimum_snr (float): (**Optional**) Minimum S/N ratio to include
            in the fitting.

        pixelmask (:class:`mangadap.proc.pixelmask.PixelMask`):
            (**Optional**) Pixel mask to include during the fitting.
        
        bias, clean, degree, mdegree, moments, sky, regul, reddening,
            component, reg_dim (various): (**Optional**) See
            :class:`mangadap.contrib.ppxf.ppxf` documentation.

    """

    def __init__(self, template_library_key, template_library, guess_redshift, guess_dispersion,
                 iteration_mode='global_template', match_resolution=True, velscale_ratio=None,
                 minimum_snr=None, pixelmask=None, bias=None, clean=False, degree=None,
                 mdegree=None, moments=None, sky=None, regul=None, reddening=None, component=None,
                 reg_dim=None):
    
        arr = [ numpy.ndarray, list ]                   # sky, p0, lam, component
        arr_in_fl = [ numpy.ndarray, list, int, float ] # component, reg_dim
        in_fl = [ int, float ]                          # Reddening, bias, regul

        _def = self._keyword_defaults()
        
        iter_opt = PPXFFit.iteration_modes()
        moment_opt = [ 2, 4, 6 ]

        pars =     [ 'template_library_key', 'template_library', 'guess_redshift',
                     'guess_dispersion', 'iteration_mode', 'match_resolution', 'velscale_ratio',
                     'minimum_snr', 'pixelmask', 'bias', 'clean', 'degree', 'mdegree', 'moments',
                     'sky', 'regul', 'reddening', 'component', 'reg_dim' ]
        values =   [ template_library_key, template_library, guess_redshift, guess_dispersion,
                     iteration_mode, match_resolution, velscale_ratio, minimum_snr, pixelmask,
                     bias, clean, degree, mdegree, moments, sky, regul, reddening, component,
                     reg_dim ]
        options =  [ None, None, None, None, iter_opt, None, None, None, None, None, None, None,
                     None, moment_opt, None, None, None, None, None ]
        defaults = [ None, None, None, None, 'global_template', True, None, None, None,
                     _def['bias'], _def['clean'], _def['degree'], _def['mdegree'], _def['moments'],
                     None, _def['regul'], _def['reddening'], _def['component'], _def['reg_dim'] ]
        dtypes =   [ str, TemplateLibrary, arr_in_fl, arr_in_fl, str, bool, int, in_fl, PixelMask,
                     in_fl, bool, int, int, int, arr, in_fl, in_fl, arr_in_fl, arr_in_fl ]

        ParSet.__init__(self, pars, values=values, defaults=defaults, options=options,
                        dtypes=dtypes)
        self._check()


    @staticmethod
    def _keyword_defaults():
        """
        Return the keyword defaults.  Pulled from
        :class:`mangadap.contrib.ppxf.ppxf`.
        """
        return { 'bias':None, 'clean':False, 'degree':4, 'mdegree':0, 'moments':2, 'regul':0,
                 'reddening':None, 'component':0, 'reg_dim':None }


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
        return hdr


    def fromheader(self, hdr):
        self['template_library_key'] = hdr['PPXFTPLK']
        self['bias'] = eval(hdr['PPXFBIAS'])
        self['clean'] = eval(hdr['PPXFCLN'])
        self['moments'] = hdr['PPXFMOM']
        self['degree'] = hdr['PPXFDEG']
        self['mdegree'] = hdr['PPXFMDEG']




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
        # Specify some mask bit names
        self.snr_flag = 'LOW_SNR'
        self.rng_flag = 'OUTSIDE_RANGE'
        self.tpl_flag = 'TPL_PIXELS'
        self.trunc_flag = 'TRUNCATED'
        self.rej_flag = 'PPXF_REJECT'
        # Logging and terminal output
        self.loggers = None
        self.quiet = False

        # Imposed fitting boundaries
        self.velocity_limits = None
        self.sigma_limits = None
        self.gh_limits = None

        # Fitting parameters
        self.tpl_wave = None
        self.tpl_sres = None
        self.tpl_flux = None
        self.ntpl = None
        self.obj_wave = None
        self.obj_sres = None
        self.obj_flux = None
        self.obj_ferr = None
        self.nobj = None
        self.velscale = None
        self.velscale_ratio = None
        self.matched_resolution = None
        self.guess_kin = None
        self.spectrum_start = None
        self.spectrum_end = None
        self.dof = None
        self.base_velocity = None

        # Fitting options
        self.bias = None
        self.clean = None
        self.degree = None
        self.mdegree = None
        self.moments = None


    @staticmethod
    def iteration_modes():
        r"""
        
        Possible iteration methods:

            ``none``: Fit all bins with all templates with a single call
            to pPXF.

            ``no_global_wrej``: Do not fit the global spectrum
            first, but include a rejection iteration.  All templates are
            fit in each step.

            ``global_template``:  Fit the global spectrum with all
            templates and include a single rejection iteration.  The
            pixel mask for this fit is the base mask for all fits to the
            individual bins.  A single rejection iteration is done for
            each bin.  **Only the global template is used when fitting
            each bin.**

            ``nonzero_templates``:  Fit the global spectrum with
            all templates and include a single rejection iteration.  The
            pixel mask for this fit is the base mask for all fits to the
            individual bins.  A single rejection iteration is done for
            each bin.  **Only the templates with non-zero weights are
            used when fitting each bin.**

            ``all_templates``:  Fit the global spectrum with all
            templates and include a single rejection iteration.  The
            pixel mask for this fit is the base mask for all fits to the
            individual bins.  A single rejection iteration is done for
            each bin.  **All templates are used when fitting each bin.**

            ``global_template_plus_emlines``:  Same as
            ``global_template``, but includes an iteration where the
            emission-lines are simultaneously fit in an effort to fix
            extrapolation errors beneath the each line.

            ``nonzero_templates_plus_emlines``:  Same as
            ``nonzero_templates``, but includes an iteration where the
            emission-lines are simultaneously fit in an effort to fix
            extrapolation errors beneath the each line.

            ``all_templates_plus_emlines``:  Same as ``all_templates``,
            but includes an iteration where the emission-lines are
            simultaneously fit in an effort to fix extrapolation errors
            beneath the each line.

        Returns:
            list: List of allowed options.
        """
        return [ 'none', 'no_global_wrej', 'global_template', 'nonzero_templates', 'all_templates',
                 'global_template_plus_emlines', 'nonzero_templates_plus_emlines',
                 'all_templates_plus_emlines' ]


    @staticmethod
    def _obj_tpl_pixelmatch(velscale, tpl_wave, velscale_ratio=None, dvtol=1e-10):
        _velscale_ratio = 1 if velscale_ratio is None else velscale_ratio
        print(numpy.absolute(velscale - spectrum_velocity_scale(tpl_wave)*_velscale_ratio), dvtol)
        return numpy.absolute(velscale - spectrum_velocity_scale(tpl_wave)*_velscale_ratio) < dvtol


    # TODO: Be clear between velocity (ppxf) vs. redshift (cz)
    @staticmethod
    def _fitting_mask(tpl_wave, obj_wave, velscale, velscale_ratio=None, waverange=None,
                      velocity_offset=None, max_velocity_range=400., alias_window=None,
                      loggers=None, quiet=False):
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
            fit_indx = numpy.logical_and(obj_wave > waverange[0], obj_wave < waverange[1])
            if numpy.sum(fit_indx) == 0:
                raise ValueError('Selected wavelength range for analysis contains no pixels!')
        else:
            fit_indx = numpy.full(now, True, dtype=numpy.bool)
        waverange_mask = ~fit_indx

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
        _velscale_ratio = 1 if velscale_ratio is None else velscale_ratio
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
        npix_mask = ~(fit_indx) & ~(waverange_mask)

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
        alias_mask = ~(fit_indx) & ~(npix_mask)

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


    def _update_rejection(self, ppxf_fit, galaxy, gpm, maxiter=9):
        """
        Reject aberrant pixels from pPXF fit.
        """
        for i in range(maxiter):
#            if not self.quiet:
#                log_output(self.loggers, 1, logging.INFO,
#                           'Rejection iteration, remaining pixels: {0} {1}'.format(i+1, len(gpm)))
            scale = numpy.sum(galaxy[gpm]*ppxf_fit.bestfit[gpm]) \
                            / numpy.sum(numpy.square(ppxf_fit.bestfit[gpm]))
            resid = scale*ppxf_fit.bestfit[gpm] - galaxy[gpm]
            err = numpy.std(resid)
            ok_old = gpm
            gpm = numpy.flatnonzero(numpy.absolute(ppxf_fit.bestfit - galaxy) < 3*err)
            gpm = numpy.intersect1d(ok_old, gpm)
            if numpy.array_equal(gpm, ok_old):
                break 
        return gpm


    def _dispersion_correction(self, obj_sres, gpm):
        """
        Calculate the dispersion corrections as the quadrature
        difference between the spectral resolution of the template and
        object spectra.
        """
        cnst = constants()
        sigma_inst_obj = astropy.constants.c.to('km/s').value/obj_sres[gpm]
        sigma_inst_tpl = astropy.constants.c.to('km/s').value/self.tpl_sres(self.obj_wave[gpm])
        return numpy.sqrt(numpy.mean(numpy.square(sigma_inst_obj)
                                     - numpy.square(sigma_inst_tpl)))/cnst.sig2fwhm


    def _is_near_bounds(self, ppxf_fit, guess_velocity, tol_frac=1e-2):
        """
        Check if the fitted kinematics are near the imposed limits.
        
        The definition of "near" is that the velocity and higher moments
        cannot be closer than the provided fraction of the total width
        to the boundary.  For the velocity dispersion, the fraction is
        done in log space.
        """

        near_bounds = numpy.zeros(self.moments, dtype=numpy.bool)
        near_lower_sigma_bound = False

#        print(ppxf_fit.sol)
#        print(guess_velocity)

        # Velocity
        _velocity_limits = self.velocity_limits + guess_velocity
        Dv = numpy.diff(_velocity_limits)[0]
        tol = Dv*tol_frac

#        print(tol, ppxf_fit.sol[0]-_velocity_limits[0], _velocity_limits[1]-ppxf_fit.sol[0])

        near_bounds[0] = ppxf_fit.sol[0]-_velocity_limits[0] < tol \
                            or _velocity_limits[1]-ppxf_fit.sol[0] < tol

        # Velocity dispersion
        Ds = numpy.diff(numpy.log10(self.sigma_limits))[0]
        tol = Ds*tol_frac

#        print(tol, numpy.log10(ppxf_fit.sol[1]) - numpy.log10(self.sigma_limits[0]),
#              numpy.log10(self.sigma_limits[1]) - numpy.log10(ppxf_fit.sol[1]))

        near_lower_sigma_bound \
                = numpy.log10(ppxf_fit.sol[1]) - numpy.log10(self.sigma_limits[0]) < tol
        near_bounds[1] = near_lower_sigma_bound \
                        or numpy.log10(self.sigma_limits[1]) - numpy.log10(ppxf_fit.sol[1]) < tol

        if self.moments == 2:
            return near_bounds, near_lower_sigma_bound

        # H3 and H4
        Dh = numpy.diff(self.gh_limits)
        tol = Dh*tol_frac
        near_bounds[2] = ppxf_fit.sol[2] - self.gh_limits[0] < tol \
                                or self.gh_limits[1] - ppxf_fit.sol[2] < tol
        near_bounds[3] = ppxf_fit.sol[3] - self.gh_limits[0] < tol \
                                or self.gh_limits[1] - ppxf_fit.sol[3] < tol
        
        if self.moments == 4:
            return near_bounds, near_lower_sigma_bound

        # H5 and H6
        near_bounds[4] = ppxf_fit.sol[4] - self.gh_limits[0] < tol \
                                or self.gh_limits[1] - ppxf_fit.sol[4] < tol
        near_bounds[5] = ppxf_fit.sol[5] - self.gh_limits[0] < tol \
                                or self.gh_limits[1] - ppxf_fit.sol[5] < tol

        return near_bounds, near_lower_sigma_bound
        
        
    def _fit_global_spectrum(self, plot=False):
        """
        Fit the global spectrum.  This:
            - Sets the base-level good pixel mask for the fits to the individual
              bins
            - Gives the template weights for the global template

        .. todo::
            - Only include spectra above a given S/N in global spectrum?
            - Allow for a number of iterations as input.
        """

#        print('fit global')

        # Sum all spectra to create the global spectrum
        global_spectrum = numpy.ma.sum(self.obj_flux, axis=0)
        global_spectrum_err = numpy.ma.sqrt(numpy.ma.sum(numpy.square(self.obj_ferr), axis=0))
        global_spectrum_err[numpy.ma.getmaskarray(global_spectrum)] = 1.0   # To avoid pPXF error

#        stk = SpectralStack()
#        _global_wave, _global_flux, _global_sdev, _global_npix, _global_ivar, _global_covar \
#                = stk.stack(self.obj_wave, self.obj_flux,
#                                      ivar=numpy.ma.power(self.obj_ferr,-2), log=True)
#
#        print(_global_wave.shape)
#        print(type(_global_wave))
#
#        print(_global_flux.shape)
#        print(type(_global_flux))
#
#        print(_global_ivar.shape)
#        print(type(_global_ivar))
#
#        pyplot.step(self.obj_wave, global_spectrum, where='mid', color='k', lw=0.5)
#        pyplot.step(_global_wave, _global_flux[0,:], where='mid', color='r', lw=0.5)
#        pyplot.show()
#        exit()

        # Initialize the good pixel mask based on the mask in the
        # fitting region

        # TODO: Because of how it's used, this will mess things up later
        # in the fit() function UNLESS all the spectra have the same
        # start and end; ie., when fitting the global spectrum, the
        # object spectra provided to fit() must be treated as an
        # ensemble.
        start = numpy.amax(self.spectrum_start)
        base_vel = self.base_velocity[ numpy.argmax(self.spectrum_start) ]
        end = numpy.amin(self.spectrum_end)
        gpm = numpy.where( ~(global_spectrum.mask[start:end]))[0]

        # Use any template that is request for any of the individual
        # spectra
        usetpl = numpy.any(self.usetpl, axis=0)

        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, 'Fitting global spectrum ...')
            log_output(self.loggers, 1, logging.INFO,
                       'Original number of fitted pixels: {0}'.format(len(gpm)))

#        pyplot.step(self.obj_wave, global_spectrum, where='mid', color='k')
#        pyplot.step(self.obj_wave, global_spectrum_err, where='mid', color='k')
#        pyplot.plot(self.obj_wave[numpy.array([self.spectrum_start,self.spectrum_start])], [0,1],
#                    color='r')
#        pyplot.plot(self.obj_wave[numpy.array([self.spectrum_end-1,self.spectrum_end-1])], [0,1],
#                    color='r')
#        pyplot.show()

        # First fit of global spectrum
        ppxf_fit = ppxf(self.tpl_flux[usetpl,:].T, global_spectrum.data[start:end],
                        global_spectrum_err.data[start:end], self.velscale, self.guess_kin[0,:],
                        velscale_ratio=self.velscale_ratio, goodpixels=gpm, bias=self.bias,
                        clean=self.clean, degree=self.degree, mdegree=self.mdegree,
                        moments=self.moments, vsyst=-base_vel, quiet=(not plot), plot=plot)
        if plot:
            pyplot.show()

        # Return if pPXF failed
        if not ppxf_fit.mp.status > 0:
            if not self.quiet:
                log_output(self.loggers, 1, logging.INFO, 'pPXF failed to fit global spectrum!')
            return None, None

#        print(ppxf_fit.sol[0] + self.base_velocity)
#        print(self.guess_kin[0,0] + self.base_velocity)
#        print(ppxf_fit.mp.status)
#        pyplot.show()

#        # Save the non-zero templates
#        nonzero_templates = ppxf_fit.weights[0:self.ntpl] > 0
#        print(numpy.sum(nonzero_templates))
#        print(numpy.where( nonzero_templates )[0])

        # Reject pixels based on their deviation from the model
        gpm = self._update_rejection(ppxf_fit, global_spectrum.data[start:end], gpm)
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO,
                       'Number of pixels to fit after first rejection: {0}'.format(len(gpm)))

        # Refit with the new bad pixel mask and with a starting guess
        # based on the previous fit
        ppxf_fit = ppxf(self.tpl_flux[usetpl,:].T, global_spectrum.data[start:end],
                        global_spectrum_err.data[start:end], self.velscale, ppxf_fit.sol[0:2],
                        velscale_ratio=self.velscale_ratio, goodpixels=gpm, bias=self.bias,
                        clean=self.clean, degree=self.degree, mdegree=self.mdegree,
                        moments=self.moments, vsyst=-base_vel, quiet=(not plot), plot=plot)
        if plot:
            pyplot.show()
#        print(ppxf_fit.sol[0] + self.base_velocity)
#        print(self.guess_kin[0,0] + self.base_velocity)
#        print(ppxf_fit.mp.status)
#        pyplot.show()

        # Return if pPXF failed
        if not ppxf_fit.mp.status > 0:
            if not self.quiet:
                log_output(self.loggers, 1, logging.INFO, 'pPXF failed to fit global spectrum!')
            return None, None

        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO,
                       'Non-zero templates fit to global spectrum: {0}'.format(
                                                    numpy.sum(ppxf_fit.weights[0:self.ntpl] > 0)))

        # Return the good pixel mask and the template weights
        weights = numpy.zeros(self.ntpl, dtype=float)
        weights[usetpl] = ppxf_fit.weights[0:numpy.sum(usetpl)]
        return gpm, weights


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

        # Select the spectra that meet the selection criteria
        # TODO: Link this to the StellarContinuumModel._bins_to_fit()
        # function...
        good_spec = (binned_spectra['BINS'].data['SNR'] > par['minimum_snr']) \
                        & ~(numpy.array([ b in binned_spectra.missing_bins \
                                                for b in numpy.arange(binned_spectra.nbins)]))
#        print(good_spec)
#        good_spec[2:] = False
#        print(good_spec)

#        for i in range(510,515):
#            pyplot.plot(obj_wave, obj_flux[good_spec,:][i,:])
#            pyplot.plot(obj_wave, obj_ferr[good_spec,:][i,:])
#            pyplot.show()

#        print(par['template_library']['SPECRES'].data.shape)
#        for i in range(par['template_library']['SPECRES'].data.shape[0]):
#            pyplot.plot(par['template_library']['WAVE'].data,
#                        par['template_library']['SPECRES'].data[i,:])

        warnings.warn('Adopting mean spectral resolution of all templates!')
        tpl_sres = numpy.mean(par['template_library']['SPECRES'].data, axis=0).ravel()
#        print(tpl_sres.shape)
#        pyplot.plot(par['template_library']['WAVE'].data, tpl_sres)
#        pyplot.show()

#        # Get the spectral resolution of the object data
#        obj_sres = binned_spectra.spectral_resolution_matrix()
#
#HERE !!

        # Perform the fit
        model_wave, model_flux, model_mask, model_par \
                = self.fit(par['template_library']['WAVE'].data.copy(),
                           par['template_library']['FLUX'].data.copy(),
                           binned_spectra['WAVE'].data.copy(), obj_flux[good_spec,:],
                           obj_ferr[good_spec,:], par['guess_redshift'][good_spec],
                           par['guess_dispersion'][good_spec], iteration_mode=par['iteration_mode'],
                           velscale_ratio=par['velscale_ratio'], mask=par['pixelmask'],
                           matched_resolution=par['match_resolution'],
                           tpl_sres=tpl_sres, obj_sres=binned_spectra['SPECRES'].data.copy(),
                           waverange=par['pixelmask'].waverange, bias=par['bias'],
                           clean=par['clean'], degree=par['degree'], mdegree=par['mdegree'],
                           moments=par['moments'], loggers=loggers, quiet=quiet, dvtol=1e-9)
#                           tpl_sres=tpl_sres, obj_sres=obj_sres[good_spec,:],

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

#        t = numpy.ma.MaskedArray(_model_flux, mask=_model_mask > 0)
#        print(_model_flux.shape)
#        print(numpy.sum(numpy.ma.getmaskarray(_model_flux)))
#        pyplot.step(model_wave, obj_flux.data[0,:], where='mid', color='0.5', lw=1)
#        pyplot.step(model_wave, obj_flux[0,:], where='mid', color='k', lw=0.5)
#        pyplot.plot(model_wave, t.data[0,:], color='purple', lw=1)
#        pyplot.plot(model_wave, t[0,:], color='r')
#        pyplot.show()

        return model_wave, _model_flux, _model_mask, _model_par

    
    def _validate_kinematics(self, speci, model_mask, model_par):
        """
        Validate the returned kinematics.

        Checks:
            - corrected velocity dispersion must be in the range 50-400
              km/s
        """
        sigcor = numpy.square(model_par['KIN'][speci,1]) \
                    - numpy.square(model_par['SIGMACORR'][speci])
        if sigcor < 2500. or sigcor > 1.6e5:
            model_mask[speci,:] = self.bitmask.turn_on(model_mask[speci,:], 'BAD_SIGMA')
            model_par['MASK'][speci] = self.bitmask.turn_on(model_par['MASK'][speci], 'BAD_SIGMA')


    def _check_mode(self, iteration_mode):
        if iteration_mode not in self.iteration_modes():
            raise ValueError('Do not understand iteration mode \'{0}\''.format(iteration_mode))
        self.iteration_mode = iteration_mode


    def _check_templates(self, tpl_wave, tpl_flux, tpl_sres):
        """
        Check that the input template data is valid.

        .. todo::
            Make log10 in call to spectral_resolution an argument to
            this function.

        """
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
        if tpl_sres is not None and tpl_sres.shape != self.tpl_wave.shape:
            raise ValueError('Provided template resolution vector does not have the correct shape.')
        self.tpl_sres = None if tpl_sres is None \
                            else spectral_resolution(self.tpl_wave, tpl_sres, log10=True)
        # TODO: Allow spectral resolution to be spectrum dependent?


    def _check_objects(self, obj_wave, obj_flux, obj_ferr, obj_sres):
        """
        Check that the input object data is valid.
        """
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
        self.obj_sres = obj_sres
        if self.obj_sres is not None and self.obj_sres.shape != self.obj_wave.shape \
                and self.obj_sres.shape != self.obj_flux.shape:
            raise ValueError('Provided object resolution vector does not have the correct shape.')
        if self.obj_sres is not None and self.obj_sres.shape != self.obj_flux.shape:
            self.obj_sres = numpy.array([self.obj_sres]*self.nobj)


    def _check_template_usage_flags(self, usetpl):
        if usetpl is not None and usetpl.shape != (self.nobj, self.ntpl) \
                and usetpl.shape != (self.ntpl,):
            raise ValueError('Provided template selection object does not have the correct shape!')
        if usetpl is None:
            self.usetpl = numpy.ones((self.nobj,self.ntpl), dtype=numpy.bool)
        else:
            self.usetpl = usetpl.astype(bool)
            if self.usetpl.shape == (self.ntpl,):
                self.usetpl = numpy.array([self.usetpl]*self.nobj)


    def _check_input_kinematics(self, guess_redshift, guess_dispersion):
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

        # Set the input redshifts
        self.input_cz = _guess_redshift*astropy.constants.c.to('km/s').value
        # Set the input pPXF, pixel-based kinematics
        # (see _convert_velocity)
        self.guess_kin = numpy.array([
                            numpy.log(_guess_redshift+1)*astropy.constants.c.to('km/s').value,
                            _guess_dispersion]).T


    def _check_pixel_scale(self, velscale_ratio, dvtol=1e-10):
        # Get the pixel scale
        self.velscale = spectrum_velocity_scale(self.obj_wave)
        self.velscale_ratio = 1 if velscale_ratio is None else velscale_ratio
        if not self._obj_tpl_pixelmatch(self.velscale, self.tpl_wave,
                                       velscale_ratio=self.velscale_ratio, dvtol=dvtol):
            raise ValueError('Pixel scale of the object and template spectra must be identical.')


    def _check_resolution_match(self, matched_resolution):
        # Confirm there is enough information to handle an unmatched
        # resolution
        if not matched_resolution and (self.obj_sres is None or self.tpl_sres is None):
            raise ValueError('If the spectral resolution is not matched between the template and '
                             'the object data, you must provide the spectral resolution for both.')
        self.matched_resolution = matched_resolution


    def _check_wavelength_range(self, waverange):
        self.waverange = numpy.array([[self.obj_wave[0]-1, self.obj_wave[-1]+1]]) \
                                if waverange is None else numpy.atleast_2d(waverange)
        if self.waverange.shape[0] == 1:
            self.waverange = numpy.array([self.waverange[0,:]]*self.nobj)
        if self.waverange.shape != (self.nobj,2):
            raise ValueError('Input wavelength range array does not have the correct shape.')


    @staticmethod
    def _losvd_limits(velscale):
        r"""
        Return the limits on the LOSVD parameters used by pPXF.

            - Velocity limits are :math:`\pm 2000` km/s
            - Velocity-disperison limits are from 1/10 pixels to 1000 km/s
            - Limits of the higher orders moments are from -0.3 to 0.3
        """
        velocity_limits = numpy.array([-2e3, 2e3])
        sigma_limits = numpy.array([0.1*velscale, 1e3])
        gh_limits = numpy.array([-0.3, 0.3])
        return velocity_limits, sigma_limits, gh_limits


    def _initialize_output(self, mask):
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
                                       velocity_offsets=self.input_cz)

        # Record array with model parameters
        # TODO: Fits using emission lines have hardwired polynomial
        # degrees; Need to fix this!
        model_par = init_record_array(self.nobj,
                        self._per_stellar_kinematics_dtype(self.ntpl, self.degree+1,
                                                           max(self.mdegree,0), self.moments,
                                                           self.bitmask.minimum_dtype())
                            if self.iteration_mode not in [ 'global_template_plus_emlines',
                                                            'nonzero_templates_plus_emlines',
                                                            'all_templates_plus_emlines' ] else
                        self._per_stellar_kinematics_dtype(self.ntpl, 0, 8, self.moments,
                                                           self.bitmask.minimum_dtype()))
        model_par['BIN_INDEX'] = numpy.arange(self.nobj)

        return model_flux, model_mask, model_par


    def _initialize_pixels_to_fit(self, model_mask, ensemble=True, max_velocity_range=400.,
                                  alias_window=None):

        err = numpy.empty(self.nobj, dtype=object)  # Empty error messages

        # Assess the regions that need to be masked during fitting
        if ensemble:
            print('ensemble is true!')
            self.waverange = numpy.array([numpy.amax(self.waverange[:,0]),
                                          numpy.amin(self.waverange[:,1])])
            try:
                fit_indx, waverange_mask, npix_mask, alias_mask \
                        = self._fitting_mask(self.tpl_wave, self.obj_wave, self.velscale,
                                             velscale_ratio=self.velscale_ratio,
                                             waverange=self.waverange,
                                             velocity_offset=self.input_cz,
                                             max_velocity_range=max_velocity_range,
                                             alias_window=alias_window, loggers=self.loggers,
                                             quiet=self.quiet)
            except ValueError as e:
                return model_mask, numpy.array([e]*self.nobj)

            fit_indx = numpy.array([ fit_indx ]*self.nobj)
            waverange_mask = numpy.array([ waverange_mask ]*self.nobj)
            npix_mask = numpy.array([ npix_mask ]*self.nobj)
            alias_mask = numpy.array([ alias_mask ]*self.nobj)
        else:
            fit_indx = numpy.zeros(self.obj_flux.shape, dtype=bool)
            waverange_mask = numpy.zeros(self.obj_flux.shape, dtype=bool)
            npix_mask = numpy.zeros(self.obj_flux.shape, dtype=bool)
            alias_mask = numpy.zeros(self.obj_flux.shape, dtype=bool)
            for i in range(self.nobj):
                try:
                    fit_indx[i,:], waverange_mask[i,:], npix_mask[i,:], alias_mask[i,:] \
                            = self._fitting_mask(self.tpl_wave, self.obj_wave, self.velscale,
                                                 velscale_ratio=self.velscale_ratio,
                                                 waverange=self.waverange[i,:],
                                                 velocity_offset=self.input_cz[i],
                                                 max_velocity_range=max_velocity_range,
                                                 alias_window=alias_window, loggers=self.loggers,
                                                 quiet=self.quiet)

                except ValueError as e:
                    model_par['MASK'][i] = self.bitmask.turn_on(model_par['MASK'], 'NO_FIT')
                    err[i] = e

        # Add these regions to the mask
        model_mask[waverange_mask] = self.bitmask.turn_on(model_mask[waverange_mask], self.rng_flag)
        model_mask[npix_mask] = self.bitmask.turn_on(model_mask[npix_mask], self.tpl_flag)
        model_mask[alias_mask] = self.bitmask.turn_on(model_mask[alias_mask], self.trunc_flag)

        # Make sure that the errors are valid
        indx = ~(self.obj_ferr.data > 0) | ~(numpy.isfinite(self.obj_ferr.data))
        if numpy.sum(indx) > 0:
            model_mask[indx] = self.bitmask.turn_on(model_mask[indx], 'INVALID_ERROR')
        # To avoid having pPXF throw an exception, make sure that all
        # of these bad errors are set to unity
        self.obj_ferr[indx] = 1.0

        # Update the internal mask of the data
        self.obj_flux[model_mask > 0] = numpy.ma.masked
        self.obj_ferr[model_mask > 0] = numpy.ma.masked

        # Determine the starting and ending pixels
        pix = numpy.ma.MaskedArray(numpy.array([numpy.arange(self.obj_flux.shape[1])]*self.nobj),
                                   mask = ~fit_indx)
        self.spectrum_start = numpy.ma.amin(pix, axis=1)
        self.spectrum_end = numpy.ma.amax(pix, axis=1)+1

        return model_mask, err


    def _mode_uses_global_spectrum(self):
        return self.iteration_mode in [ 'global_template', 'global_template_plus_emlines',
                                        'nonzero_templates', 'nonzero_templates_plus_emlines',
                                        'all_templates', 'all_templates_plus_emlines' ]

            
    def _mode_includes_rejection(self):
        return self.iteration_mode != 'none'


    def _set_and_report_failed_status(self, model_mask, model_par, message):
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, message)
        return self.bitmask.turn_on(model_mask, 'FIT_FAILED'), \
               self.bitmask.turn_on(model_par, 'FIT_FAILED')


    def fit(self, tpl_wave, tpl_flux, obj_wave, obj_flux, obj_ferr, guess_redshift,
            guess_dispersion, iteration_mode='global_template', ensemble=True, velscale_ratio=None,
            mask=None, usetpl=None, matched_resolution=True, tpl_sres=None, obj_sres=None,
            waverange=None, bias=None, clean=False, degree=4, mdegree=0, moments=2, loggers=None,
            quiet=False, max_velocity_range=400., alias_window=None, dvtol=1e-10, plot=False):
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
            iteration_mode (str): (**Optional**) Iteration sequence to
                perform.  See :func:`iteration_modes`.
            ensemble (bool): (**Optional**) Treat the list of input
                spectra as an ensemble.  Currently, the only affects how
                the spectra are masked.  Default is to treat them as an
                ensemble.  When not treated as an ensemble, each
                spectrum is masked individually according to the input
                wavelength range and velocity offsets.  *It does not
                make sense to set the iteration_mode to something that
                will include a fit to the global spectrum if you're not
                treating the list of object spectra as an ensemble.*
            velscale_ratio (int): (**Optional**) Ratio of velocity scale
                per pixel in the object spectra to that for the template
                spectra.  Default is that they should be identical.
            mask (numpy.ndarray): (**Optional**) A
                baseline pixel mask to use during the fitting  Other
                pixels may be masked via the convenience functions, but
                these pixels will always be masked.
            matched_resolution (bool): (**Optional**)  Flag that the
                object and template spectra have identical spectral
                resolution.  Default is True.
            tpl_sres (numpy.ndarray): (**Optional**) One-dimensional
                vector with the spectral resolution (:math:`R =
                \lambda/\Delta\lambda`) at each wavelength of the
                template spectra.  Default is the resolution is not
                provided and assumed to be same as the object
                resolution.
            obj_sres (numpy.ndarray): (**Optional**) One- or Two-dimensional
                array with the spectral resolution (:math:`R =
                \lambda/\Delta\lambda`) at each wavelength for (each of)
                the object spectra.  Default is the resolution is not
                provided and assumed to be same as the template
                resolution.
            waverange (array-like): (**Optional**) Lower and upper
                wavelength limits to *include* in the fit.  This can be
                a two-element vector to apply the same limits to all
                spectra, or a N-spec x 2 array with wavelength ranges
                for each spectrum to be fit.  Default is to use as much
                of the spectrum as possible.
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

            loggers (list): (**Optional**) List of `logging.Logger`_ objects
                to log progress; ignored if quiet=True.  Logging is done
                using :func:`mangadap.util.log.log_output`.  Default is
                no logging.
            quiet (bool): (**Optional**) Suppress all terminal and
                logging output.  Default is False.
            max_velocity_range (float): (**Optional**) Maximum range
                (+/-) expected for the fitted velocities in km/s.
                Default is 400 km/s.
            alias_window (float) : (**Optional**) The window to mask to
                avoid aliasing near the edges of the spectral range in
                km/s.  Default is six times *max_velocity_range*.
            dvtol (float): (**Optional**) The velocity scale of the
                template spectra and object spectrum must be smaller
                than this tolerance.  Default is 1e-10.
            plot (bool): (**Optional**) Show the automatically generated
                pPXF fit plots.
                
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

        # Check the input; following operation is necessarily sequential
        # in some cases because of assigned attributes
        # - Mode
        self._check_mode(iteration_mode)
        # - Templates
        self._check_templates(tpl_wave, tpl_flux, tpl_sres)
        # - Objects
        self._check_objects(obj_wave, obj_flux, obj_ferr, obj_sres)
        # - Template usage (needs to know the number of object spectra)
        self._check_template_usage_flags(usetpl)
        # - Input kinematics
        self._check_input_kinematics(guess_redshift, guess_dispersion)
        # - Input pixel scale
        self._check_pixel_scale(velscale_ratio, dvtol=dvtol)
        # - Spectral resolution
        self._check_resolution_match(matched_resolution)
        # - Selected wavelength range: always has shape (self.nobj,2)
        self._check_wavelength_range(waverange)

        # Set the parameter limits
        self.velocity_limits, self.sigma_limits, self.gh_limits \
                    = PPXFFit._losvd_limits(self.velscale)

        # Save the pPXF parameters
        self.bias = bias
        self.clean = clean
        self.degree = degree
        self.mdegree = mdegree
        self.moments = moments

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

        # Initialize the output data
        model_flux, model_mask, model_par = self._initialize_output(mask)

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

        # Initialize the mask and the spectral range to fit
        model_mask, err = self._initialize_pixels_to_fit(model_mask, ensemble=_ensemble,
                                                         max_velocity_range=max_velocity_range,
                                                         alias_window=alias_window)
        ended_in_error = numpy.array([e is not None for e in err])
        if numpy.any(ended_in_error):
            if not self.quiet:
                warnings.warn('Masking failures in some/all spectra.  Errors are: {0}'.format(
                                numpy.array([(i,e) for i,e in enumerate(err)])[ended_in_error]))
            model_par['MASK'][ended_in_error] \
                    = self.bitmask.turn_on(model_par['MASK'][ended_in_error], 'NO_FIT')
        if numpy.all(ended_in_error):
            return self.obj_wave, model_flux, model_mask, model_par

#        if not self.quiet:
#            log_output(self.loggers, 1, logging.INFO, 'Fitting pixels {0} --> {1}'.format(
#                       self.spectrum_start, self.spectrum_end))
#            log_output(self.loggers, 1, logging.INFO,
#                       'Corresponds to wavelengths {0} --> {1}'.format(
#                                                            self.obj_wave[self.spectrum_start],
#                                                            self.obj_wave[self.spectrum_end-1]))

#        _test_mask = self.bitmask.turn_off(model_mask, 'EML_REGION')
#        pyplot.step(self.obj_wave, self.obj_flux[0,:], where='mid', lw=0.5, color='k')
#        pyplot.step(self.obj_wave, model_mask[0,:], where='mid', lw=1., color='r')
#        pyplot.step(self.obj_wave, _test_mask[0,:], where='mid', lw=1., color='g')
#        pyplot.show()

        # Save the pixel statistics
        model_par['BEGPIX'] = self.spectrum_start
        model_par['ENDPIX'] = self.spectrum_end
        model_par['NPIXTOT'] = self.spectrum_end - self.spectrum_start

#        pyplot.step(self.obj_wave, self.obj_flux[0,:]/numpy.median(self.obj_flux[0,:]),
#                    where='mid', color='k', lw=0.5)
#        pyplot.step(self.tpl_wave*(1+guess_redshift[0]),
#                    self.tpl_flux[31,:]/numpy.median(self.tpl_flux[31,:]), where='mid',
#                    color='red', lw=1.5)
#        pyplot.show()

        # Determine the degrees of freedom of the fit, used for the
        # brute force calculation of the reduced chi-square
        self.dof = self.moments + max(self.mdegree, 0)
        if self.degree >= 0:
            self.dof += self.degree+1
        if not self.quiet:
            log_output(self.loggers, 2, logging.INFO, 'Total degrees of freedom: {0}'.format(
                                                                            self.dof+self.ntpl))

        # Get the input pixel shift between the object and template
        # wavelength vectors; interpretted by pPXF as a base velocity
        # shift between the two
        self.base_velocity = numpy.array([self.ppxf_tpl_obj_voff(self.tpl_wave, self.obj_wave[s:e],
                                                                 self.velscale,
                                                                 velscale_ratio=self.velscale_ratio)
                                                for s,e in zip(self.spectrum_start,
                                                               self.spectrum_end)])

        # Report the base velocity
#        if not self.quiet:
#            log_output(self.loggers, 2, logging.INFO,
#                       'Initial wavelength yields baseline velocity offset of: {0:.2f}'.format(
#                                                                               self.base_velocity))

#        pyplot.imshow(self.obj_flux, origin='lower', interpolation='nearest', aspect='auto')
#        pyplot.colorbar()
#        pyplot.show()

        # Initialize the good pixel mask and the template set according
        # to the iteration method
#        plot = False
#        print('Fit to global spectrum')

        # Fit the global spectrum if requested by the iteration mode
        if self._mode_uses_global_spectrum():
            gpm0, global_weights0 = self._fit_global_spectrum(plot=plot)
            if gpm0 is None or global_weights0 is None:
                # pPXF failed!  
                model_mask[:] = self.bitmask.turn_on(model_mask[:], 'FIT_FAILED')
                model_par['MASK'][:] = self.bitmask.turn_on(model_par['MASK'][:], 'FIT_FAILED')
                return self.obj_wave, model_flux, model_mask, model_par

        # TODO: Warn user that selection of templates doesn't work with
        # using the global template

        # Use the global template
        if self.iteration_mode in [ 'global_template', 'global_template_plus_emlines' ]:
            templates = numpy.dot(global_weights0, self.tpl_flux).reshape(1,-1)
            tpl_to_fit = numpy.ones(1, dtype=numpy.bool)
        # Use the non-zero templates
        elif self.iteration_mode in [ 'nonzero_templates', 'nonzero_templates_plus_emlines' ]:
            templates = self.tpl_flux
            tpl_to_fit = numpy.array([global_weights0 > 0]*self.nobj) & self.usetpl
        # Use all templates
        elif self.iteration_mode in [ 'all_templates', 'all_templates_plus_emlines', 'none',
                                      'no_global_wrej' ]:
            templates = self.tpl_flux
            tpl_to_fit = self.usetpl
        else:
            raise ValueError('Unrecognized iteration method: {0}'.format(self.iteration_mode))

        # Print heading of fitted kinematics for each bin
        if not self.quiet:
            log_output(self.loggers, 2, logging.INFO, '{0:>5}'.format('INDX')
                       + ' {:>7} {:>7}'.format('SWAVE', 'EWAVE')
                       + (' {:>9}'*self.moments).format(*(['KIN{0}'.format(i+1)
                                                                for i in range(self.moments)]))
                       + ' {0:>5} {1:>9} {2:>5} {3:>4}'.format('NZTPL', 'CHI2', 'RCHI2', 'STAT'))

#        indx = ~(self.obj_ferr > 0) | ~(numpy.isfinite(self.obj_ferr))
#        print(numpy.sum(indx))

#        plot=True
        # Fit each binned spectrum:
#        for i in range(20,21):
        for i in range(5):
#        for i in range(self.nobj):

            if self.iteration_mode not in [ 'none', 'no_global_wrej']:
                # Reset the good-pixel mask to the default, and the
                # global weights to the original values
                gpm = gpm0.copy()
                global_weights = global_weights0.copy()
            else:
                # Set the good-pixel mask for this spectrum specifically
                gpm = numpy.where(~(self.obj_flux.mask[i,
                                            self.spectrum_start[i]:self.spectrum_end[i]]))[0]
        
            # Check if there is sufficient data for the fit
            # TODO: Make sure this is a reasonable approach
            ntpl_to_fit = numpy.sum(tpl_to_fit[i,:])
            if len(gpm) < self.dof+ntpl_to_fit:
                if not self.quiet:
                    warnings.warn('Insufficient data points ({0}) to fit spectrum {1}'
                                  '(dof={2}).'.format(len(gpm), i+1, self.dof+ntpl_to_fit))
                # Flag that the fit was not performed
                model_mask[i,:] = self.bitmask.turn_on(model_mask[i,:], 'NO_FIT')
                model_par['MASK'][i] = self.bitmask.turn_on(model_par['MASK'][i],
                                                            'INSUFFICIENT_DATA')
                # Move on to the next bin
                continue
            
#            pyplot.plot(self.obj_wave[self.spectrum_start:self.spectrum_end],
#            self.obj_flux[i,self.spectrum_start:self.spectrum_end])
#            pyplot.plot(self.obj_wave[self.spectrum_start:self.spectrum_end],
#            self.obj_flux[i,self.spectrum_start:self.spectrum_end])
#            pyplot.show()

#            pyplot.plot(self.obj_wave[self.spectrum_start:self.spectrum_end],
#                        self.obj_flux[i,self.spectrum_start:self.spectrum_end]/numpy.ma.mean(
#                                                            self.obj_flux[i,self.spectrum_start:self.spectrum_end]))
#            pyplot.plot(self.tpl_wave*(1+(self.vsyst-self.base_velocity)
#                            / astropy.constants.c.to('km/s').value),
#                        self.tpl_flux[31,:]/numpy.mean(self.tpl_flux[31,:]))
#            pyplot.plot(self.tpl_wave*(1+self.base_velocity/astropy.constants.c.to('km/s').value),
#                        self.tpl_flux[31,:]/numpy.mean(self.tpl_flux[31,:]))
#            pyplot.show()
#            print(self.guess_kin[i,:])

            # Run ppxf
#            print('fit individual')
            ppxf_fit = ppxf(templates[tpl_to_fit[i,:],:].T,
                            self.obj_flux.data[i,self.spectrum_start[i]:self.spectrum_end[i]],
                            self.obj_ferr.data[i,self.spectrum_start[i]:self.spectrum_end[i]],
                            self.velscale, self.guess_kin[i,:], velscale_ratio=self.velscale_ratio,
                            goodpixels=gpm, bias=self.bias, clean=self.clean, degree=self.degree,
                            mdegree=self.mdegree, moments=self.moments,
                            vsyst=-self.base_velocity[i], quiet=(not plot), plot=plot)
            if plot:
                pyplot.show()

            # Fit failed
            if not ppxf_fit.mp.status > 0:
                model_mask[i,:], model_par['MASK'][i] = self._set_and_report_failed_status(
                                model_mask[i,:], model_par['MASK'][i],
                                'Initial pPXF status for spectrum {0}; nothing saved.'.format(i+1))
                continue


            # Add the pixels to the DIDNOTUSE mask

# Check the validity of the 'reconstruct_model' function
#            tmpkin = ppxf_fit.sol
#            tmpkin[0], terr = self._convert_velocity(tmpkin[0], 0.0)
#
#            tmpwgts = numpy.zeros(self.ntpl, dtype=numpy.float)
#            if self.iteration_mode == 'global_template':
#                tmpwgts = ppxf_fit.weights[0]*global_weights
#            elif self.iteration_mode == 'nonzero_templates':
#                tmpwgts[nonzero_templates] = ppxf_fit.weights[0:templates.shape[0]]
#            else:
#                tmpwgts = ppxf_fit.weights[0:self.ntpl]
#            
#            rbestfit = self.reconstruct_model(self.tpl_wave, self.tpl_flux, self.obj_wave, tmpkin,
#                                              tmpwgts, self.velscale,
#                                              polyweights=None if degree < 0 
#                                                               else ppxf_fit.polyweights,
#                                              mpolyweights=None if mdegree < 1
#                                                                else ppxf_fit.mpolyweights,
#                                              start=self.spectrum_start, end=self.spectrum_end,
#                                              velscale_ratio=self.velscale_ratio, dvtol=dvtol)
#            pyplot.plot(self.obj_wave[self.spectrum_start:self.spectrum_end], ppxf_fit.bestfit)
#            pyplot.plot(self.obj_wave, rbestfit)
#            pyplot.show()

            # Perform rejection iteration
            if self._mode_includes_rejection():
                # Reject
                gpm_rej = self._update_rejection(ppxf_fit, self.obj_flux.data[i,
                                                    self.spectrum_start[i]:self.spectrum_end[i]],
                                                    gpm)
                # Update input kinematics
                self.guess_kin[i,:] = ppxf_fit.sol[0:2]
                # Refit
                ppxf_fit = ppxf(templates[tpl_to_fit[i,:],:].T,
                                self.obj_flux.data[i,self.spectrum_start[i]:self.spectrum_end[i]],
                                self.obj_ferr.data[i,self.spectrum_start[i]:self.spectrum_end[i]],
                                self.velscale, self.guess_kin[i,:],
                                velscale_ratio=self.velscale_ratio, goodpixels=gpm_rej,
                                bias=self.bias, clean=self.clean, degree=self.degree,
                                mdegree=self.mdegree, moments=self.moments,
                                vsyst=-self.base_velocity[i], quiet=(not plot), plot=plot)
                if plot:
                    pyplot.show()

                # Fit failed
                if not ppxf_fit.mp.status > 0:
                    model_mask[i,:], model_par['MASK'][i] = self._set_and_report_failed_status(
                                    model_mask[i,:], model_par['MASK'][i],
                                   'Rejection iteration pPXF status for spectrum {0}; '
                                   'nothing saved.'.format(i+1))
                    continue

            #-----------------------------------------------------------
            # REFIT WITH EMISSION LINES:
            #   - updated mask; emission-lines unmasked; stellar
            #   kinematics fixed
            # TODO: This section is out of date!
            if self.iteration_mode in [ 'global_template_plus_emlines',
                                        'nonzero_templates_plus_emlines',
                                        'all_templates_plus_emlines' ]:
                raise NotImplementedError('Iterations including emission lines have been '
                                          'deprecated for the time-being.')

                # Construct a set of Gaussian emission line templates
                z = (numpy.exp(ppxf_fit.sol[0]/astropy.constants.c.to('km/s').value)-1.0)
                obj_rest_waverange = numpy.array([self.obj_wave[self.spectrum_start[i]],
                                                  self.obj_wave[self.spectrum_end[i]-1] ]) / (1+z)
                                                    #_guess_redshift[i])
                # TODO: Just adopting FWHM = 2.5 for now
                gas_templates, line_names, line_wave = \
                        ppxf_util.emission_lines(numpy.log(self.tpl_wave), obj_rest_waverange, 2.5)
#                print('GAS')
#                print( gas_templates.shape )
#                print( numpy.sum(gas_templates, axis=0) )
#                for ii in range(gas_templates.shape[1]):
#                    pyplot.plot(numpy.arange(gas_templates.shape[0]), gas_templates[:,ii])
#                pyplot.show()

                # Construct the best-fitting stellar template based on
                # the first fit

                # TODO: Need to check that this works with the change in
                # how the templates are selected for individual spectra
                if self.iteration_mode == 'global_template_plus_emlines':
                    global_weights *= ppxf_fit.weights[0]
#                elif self.iteration_mode == 'nonzero_templates_plus_emlines':
#                    global_weights[tpl_to_fit[i,:]] = ppxf_fit.weights[0:ntpl_to_fit]
                    #templates.shape[0]]
                else:
                    global_weights[tpl_to_fit[i,:]] = ppxf_fit.weights[0:ntpl_to_fit]
#                    global_weights = ppxf_fit.weights[0:self.ntpl]
                global_template = numpy.dot(global_weights, self.tpl_flux).reshape(1,-1)

#                pyplot.plot(self.tpl_wave, global_template[0,:])
#                pyplot.show()

                # Include emission-line templates with others
                _templates = numpy.row_stack([global_template, gas_templates.T])
#                for ii in range(_templates.shape[0]):
#                    pyplot.plot(numpy.arange(_templates.shape[1]), _templates[ii,:])
#                pyplot.show()

                # Set the component values
                _component = [0]*1 + [1]*gas_templates.shape[1]
                # do not fit stars but use previous solution
                _moments = [-self.moments, 2]
                # initialize the gas kinematics based on the stars
                _guess_kin = [ppxf_fit.sol[0:2], [ppxf_fit.sol[0], 50]]
#                _guess_kin = [ppxf_fit.sol[0:2], ppxf_fit.sol[0:2]]
                # unmask the emission lines
                _gpm = numpy.union1d(gpm,
                                     numpy.where(self.bitmask.flagged(
                                        model_mask[i,self.spectrum_start[i]:self.spectrum_end[i]],
                                        flag='EML_REGION'))[0])
                # Fit with emission lines, using only multiplicative polynomials
                _ppxf_fit = ppxf(_templates.T,
                                 self.obj_flux.data[i,self.spectrum_start[i]:self.spectrum_end[i]],
                                 self.obj_ferr.data[i,self.spectrum_start[i]:self.spectrum_end[i]],
                                 self.velscale, _guess_kin, velscale_ratio=self.velscale_ratio,
                                 goodpixels=_gpm, bias=self.bias, clean=self.clean, degree=-1,
                                 mdegree=8, moments=_moments, vsyst=-self.base_velocity[i],
                                 quiet=(not plot), plot=plot, component=_component)
                if plot:
                    pyplot.show()

                # Fit failed
                if not ppxf_fit.mp.status > 0:
                    model_mask[i,:], model_par['MASK'][i] = self._set_and_report_failed_status(
                                    model_mask[i,:], model_par['MASK'][i],
                                   'Emission-line iteration pPXF status for spectrum {0}; '
                                   'nothing saved.'.format(i+1))
                    continue

                # Save the necessary new data to the old fit (without
                # emission lines)
                ppxf_fit.weights = global_weights * _ppxf_fit.weights[0]
                ppxf_fit.polyweights = None
                ppxf_fit.mpolyweights = _ppxf_fit.mpolyweights

#                old_bestfit = ppxf_fit.bestfit
                ppxf_fit.bestfit = self.reconstruct_model(self.tpl_wave, self.tpl_flux,
                                        self.obj_wave, ppxf_fit.sol, ppxf_fit.weights[0:self.ntpl],
                                        self.velscale, polyweights=ppxf_fit.polyweights,
                                        mpolyweights=ppxf_fit.mpolyweights,
                                        start=self.spectrum_start[i], end=self.spectrum_end[i],
                                        velscale_ratio=self.velscale_ratio, dvtol=dvtol,
                                revert_velocity=False)[self.spectrum_start[i]:self.spectrum_end[i]]
#                pyplot.plot(self.obj_wave[self.spectrum_start:self.spectrum_end], ppxf_fit.bestfit)
#                pyplot.plot(self.obj_wave[self.spectrum_start:self.spectrum_end], _ppxf_fit.bestfit)
#                pyplot.plot(self.obj_wave[self.spectrum_start:self.spectrum_end], old_bestfit)
#                pyplot.show()
            #-------------------------------------------------------


            # Check if the fit reached the maximum number of iterations.
            # Keep the data but throw a warning and log the status via
            # the bitmask.
            if ppxf_fit.mp.status == 5:
                if not self.quiet:
                    warnings.warn('pPXF optimizer reached maximum number of iterations for '
                                  'spectrum {0}.'.format(i+1))
                model_par['MASK'][i] = self.bitmask.turn_on(model_par['MASK'][i], 'MAXITER')

            # Test if the kinematics are near the imposed boundaries.
            # The best fitting model is still provided, but masked.
            near_bound, near_lower_sigma_bound = self._is_near_bounds(ppxf_fit, self.guess_kin[i,0])

            # If the velocity dispersion has hit the lower limit, ONLY
            # flag the value as having a MIN_SIGMA.
            if near_lower_sigma_bound:
                near_bound[1] = False
                model_par['MASK'][i] = self.bitmask.turn_on(model_par['MASK'][i], 'MIN_SIGMA')
            # Otherwise, flag both the parameter set and the spectrum as
            # NEAR_BOUND
            if numpy.any(near_bound):
                if not self.quiet:
                    warnings.warn('Returned parameters for spectrum {0} too close to '
                                  'bounds.'.format(i+1))
                model_par['MASK'][i] = self.bitmask.turn_on(model_par['MASK'][i], 'NEAR_BOUND')
                model_mask[i,:] = self.bitmask.turn_on(model_mask[i,:], 'NEAR_BOUND')

            # Save the result
            model_flux[i,self.spectrum_start[i]:self.spectrum_end[i]] = ppxf_fit.bestfit
            if self._mode_includes_rejection() and len(gpm_rej) != len(ppxf_fit.goodpixels):
                rejected_pixels = list(set(gpm_rej)-set(ppxf_fit.goodpixels))
                if not self.quiet:
                    log_output(self.loggers, 1, logging.INFO,
                               'Rejected {0} pixels during fit.'.format(len(rejected_pixels)))
                model_mask[i,self.spectrum_start[i]:self.spectrum_end[i]][rejected_pixels] \
                        = self.bitmask.turn_on(model_mask[i,self.spectrum_start[i]:
                                                            self.spectrum_end[i]][rejected_pixels],
                                               flag=self.rej_flag)

            # Mask the pixels that were never used
            masked_pixels = list(set(numpy.arange(self.spectrum_end[i]-self.spectrum_start[i]))
                                    - set(gpm))
            model_mask[i,self.spectrum_start[i]:self.spectrum_end[i]][masked_pixels] \
                        = self.bitmask.turn_on(model_mask[i,self.spectrum_start[i]:
                                                            self.spectrum_end[i]][masked_pixels],
                                               flag='DIDNOTUSE')

            model_par['NPIXFIT'][i] = len(ppxf_fit.goodpixels)

            # TODO: Store template weights sparsely?
            if self.iteration_mode in [ 'global_template', 'global_template_plus_emlines']:
                model_par['TPLWGT'][i] = ppxf_fit.weights[0]*global_weights
#            elif self.iteration_mode == 'nonzero_templates' \
#                    or self.iteration_mode == 'nonzero_templates_plus_emlines':
#                model_par['TPLWGT'][i][tpl_to_fit[i,:]] = ppxf_fit.weights[0:ntpl_to_fit]
                #templates.shape[0]]
            else:
                model_par['TPLWGT'][i][tpl_to_fit[i,:]] = ppxf_fit.weights[0:ntpl_to_fit]
#                model_par['TPLWGT'][i] = ppxf_fit.weights[0:self.ntpl]

            model_par['USETPL'][i] = tpl_to_fit[i,:]

            if self.degree >= 0 and ppxf_fit.polyweights is not None:
                model_par['ADDCOEF'][i] = ppxf_fit.polyweights
            if self.mdegree > 0 and ppxf_fit.mpolyweights is not None:
                model_par['MULTCOEF'][i] = ppxf_fit.mpolyweights

            model_par['KININP'][i] = self.guess_kin[i,:]
            model_par['KIN'][i] = ppxf_fit.sol
            model_par['KINERR'][i] = ppxf_fit.error

            resid = self.obj_flux[i,self.spectrum_start[i]:self.spectrum_end[i]] - ppxf_fit.bestfit
            model_par['CHI2'][i] = numpy.sum(numpy.square( (resid 
                                        / self.obj_ferr[i,self.spectrum_start[i]:
                                                    self.spectrum_end[i]])[ppxf_fit.goodpixels]))
            model_par['RCHI2'][i] = model_par['CHI2'][i] \
                                        / (model_par['NPIXFIT'][i] 
                                            - self.dof - numpy.sum(model_par['TPLWGT'][i] > 0))
            model_par['ROBUST_RCHI2'][i] = ppxf_fit.chi2

            model_par['RMS'][i] = numpy.sqrt(numpy.ma.mean( numpy.square(
                                                                resid[ppxf_fit.goodpixels]) ))
            # Get growth statistics for the residuals
            model_par['ABSRESID'][i] = residual_growth(resid[ppxf_fit.goodpixels],
                                                       [0.68, 0.95, 0.99])

            indx = numpy.absolute(ppxf_fit.bestfit) > 0
            _goodpixels = numpy.intersect1d(numpy.arange(len(resid))[indx], ppxf_fit.goodpixels,
                                            assume_unique=True)
            if len(_goodpixels) > 0:
                frac_resid = resid[_goodpixels]/ppxf_fit.bestfit[_goodpixels]
                model_par['FRMS'][i] = numpy.sqrt(numpy.ma.mean(numpy.square(frac_resid)))
                if len(_goodpixels) > 1:
                    model_par['FABSRESID'][i] = residual_growth(frac_resid, [0.68, 0.95, 0.99])

#            pyplot.step(obj_wave, flux[i,:], where='mid', linestyle='-', lw=0.5, color='k',
#                        zorder=3)
#            pyplot.plot(obj_wave, model_flux[i,:], linestyle='-', lw=1.5, color='r',
#                        zorder=1, alpha=0.5)
#            pyplot.show()

            # Calculate the dispersion correction if necessary
            if not self.matched_resolution:
                model_par['SIGMACORR'][i] = self._dispersion_correction(self.obj_sres[i],
                                                                        ppxf_fit.goodpixels)

            # Test if kinematics are reliable
            self._validate_kinematics(i, model_mask, model_par)

            if not self.quiet:
                log_output(self.loggers, 2, logging.INFO, '{0:>5d}'.format(i+1)
                           + ' {:>7.2f} {:>7.2f}'.format(self.obj_wave[self.spectrum_start[i]],
                                                         self.obj_wave[self.spectrum_end[i]-1])
                           + (' {:>9.2f}'*self.moments).format(*model_par['KIN'][i]) 
                           + ' {0:>4d} {1:>9.2f} {2:>5.2f} {3:>4d}'.format(
                                                            numpy.sum(model_par['TPLWGT'][i] > 0),
                                                            model_par['CHI2'][i],
                                                            model_par['RCHI2'][i],
                                                            ppxf_fit.mp.status))

        model_par['KININP'][:,0], _ = self._convert_velocity(model_par['KININP'][:,0],
                                                           numpy.zeros(self.nobj))
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



#    @staticmethod
#    def _losvd_kernel(par, velscale, vsyst=0.0):
#        """
#        Return the sampled kernel.
#
#        Meant to match computation of model in pPXF, see
#        :func:`mangadap.contrib.ppxf._fitfunc`
#
#        CURRENTLY WILL NOT WORK WITH:
#            - multiple kinematic components
#            - point-symmetric fits
#            - oversampling
#
#        Input expected to be the parameters *after* the pPXF velocities
#        had been converted to redshifts.
#
#        """
#        _par = par.copy()
#        _par[0], verr = PPXFFit._revert_velocity(_par[0], 1.0)
#        _par[0] += vsyst
#        _par[0:2] /= velscale      # Convert to pixels
#
#        dx = int(abs(_par[0])+5.0*par[1])
#        nl = 2*dx+1
#        x = numpy.linspace(-dx,dx,nl)
#        w = (x - _par[0])/_par[1]
#        w2 = numpy.square(w)
#        gauss = numpy.exp(-0.5*w2)
#        losvd = gauss/numpy.sum(gauss)
#
#        # Add the higher order terms
#        if len(_par) > 2:
#            # h_3 h_4
#            poly = 1 + _par[2]/numpy.sqrt(3)*(w*(2*w2-3)) + _par[3]/numpy.sqrt(24)*(w2*(4*w2-12)+3)
#            if len(_par) > 4:
#                # h_5 h_6
#                poly += _par[4]/numpy.sqrt(60)*(w*(w2*(4*w2-20)+15)) \
#                            + _par[5]/numpy.sqrt(720)*(w2*(w2*(8*w2-60)+90)-15)
#            losvd *= poly
#
##        return losvd
#        return _par[0]*velscale, losvd
#
##        npad = int(2**numpy.ceil(numpy.log2(nwave + nl/2)))
##        losvd_pad = numpy.zeros(npad)

        
    @staticmethod
    def _shift_via_resample(wave, flux, redshift, inLog, outRange, outpix, outLog):
        inRange = wave[[0,-1]] * (1.0 + redshift)
        new_wave, new_flux = resample_vector(flux, xRange=inRange, inLog=inLog, newRange=outRange,
                                             newpix=outpix, newLog=outLog, flat=False)
        return new_flux


    @staticmethod
    def reconstruct_model(tpl_wave, templates, obj_wave, kin, weights, velscale, polyweights=None,
                          mpolyweights=None, start=None, end=None, redshift_only=False,
                          sigma_corr=0.0, velscale_ratio=None, dvtol=1e-10, revert_velocity=True):

        moments = kin.size

        _start = 0 if start is None else start
        _end = obj_wave.size if end is None else end

        _velscale_ratio = 1 if velscale_ratio is None else velscale_ratio

        outRange = obj_wave[ [0,-1] ]
        outpix = obj_wave.size

        model = numpy.ma.zeros(outpix, dtype=numpy.float)

        # Get the composite template
        indx = weights > 0
        redshift = kin[0]/astropy.constants.c.to('km/s').value
        composite_template = numpy.dot(weights, templates)

#        pyplot.step(tpl_wave, composite_template, where='mid', linestyle='-', lw=0.5,
#                    color='k')
#        pyplot.show()

        # Redshift and broadened using the LOSVD parameters
        if redshift_only:
            # Resample the redshifted template to the wavelength grid of
            # the binned spectra
            model[:] = PPXFFit._shift_via_resample(tpl_wave, composite_template, redshift, True,
                                                   outRange, outpix, True)
        else:
#            print('{0}/{1}'.format(i+1,binned_spectra.nbins))
#            print('Start,End: ', start, end)
#            print(obj_wave[start:end].shape)
#            print(tpl_wave.shape)

            # for PPXF v6.0.0
            vsyst = -PPXFFit.ppxf_tpl_obj_voff(tpl_wave, obj_wave[_start:_end], velscale,
                                               velscale_ratio=_velscale_ratio)
            if _velscale_ratio > 1:
                npix_temp = composite_template.size - composite_template.size % _velscale_ratio
                _composite_template = composite_template[:npix_temp].reshape(-1,1)
                npix_temp //= _velscale_ratio
            else:
                _composite_template = composite_template
                npix_temp = composite_template.size

            ctmp_rfft, npad = _templates_rfft(_composite_template)
            _moments = numpy.atleast_1d(moments)
            ncomp = 1
            if _moments.size == 1:
                _moments = numpy.full(ncomp, numpy.absolute(_moments), dtype=int)
            else:
                _moments = numpy.absolute(_moments)
            vj = numpy.append(0, numpy.cumsum(_moments)[:-1])
            par = kin.copy()
            # Convert the velocity to pixel units
            if revert_velocity:
                par[0], verr = PPXFFit._revert_velocity(par[0], 1.0)
            # Convert the velocity dispersion to ignore the
            # resolution difference
            if sigma_corr > 0:
                par[1] = numpy.square(par[1]) - numpy.square(sigma_corr)
                par[1] = numpy.sqrt(par[1]) if par[1] > 0 else 0.0

            # Construct the model
            if par[1] > 0:
                # Velocity dispersion is valid, follow the pPXF
                # method
                par /= velscale
                kern_rfft = _losvd_rfft(par, 1, _moments, vj, npad, 1, vsyst/velscale,
                                        _velscale_ratio, 0.0)
                _model = numpy.fft.irfft(ctmp_rfft[:,0]
                                         * kern_rfft[:,0,0]).ravel()[:npix_temp*_velscale_ratio]
                if _velscale_ratio > 1:
                    _model = numpy.mean(_model.reshape(-1,_velscale_ratio), axis=1).ravel()
#                pyplot.plot(obj_wave[_start:_end], _model[:_end-_start])
#                pyplot.show()
                model[_start:_end] = _model[:_end-_start]
            else:
                # Velocity dispersion is undefined, return a
                # zero-dispersion spectrum
                model[:] = PPXFFit._shift_via_resample(tpl_wave, composite_template, redshift,
                                                       True, outRange, outpix, True)

        # Account for the polynomial
        x = numpy.linspace(-1, 1, _end-_start)
        if mpolyweights is not None:
            model[_start:_end] *= numpy.polynomial.legendre.legval(x,numpy.append(1.0,mpolyweights))
        if polyweights is not None:
            model[_start:_end] += numpy.polynomial.legendre.legval(x,polyweights)

        return model


    @staticmethod
    def construct_models(binned_spectra, template_library, model_par, degree, mdegree, moments,
                         redshift_only=False, velscale_ratio=None, dvtol=1e-10):
        """
        Construct models using the provided set of model parameters and
        template library.  The number of template weights in the
        provided model parameters must match the number of template
        spectra.

        If redshift_only is true, ignore the LOSVD and simply shift the
        model to the correct redshift.

        CURRENTLY WILL NOT WORK WITH:
            - multiple kinematic components
            - point-symmetric fits

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

        tpl_wave = template_library['WAVE'].data.copy()
        obj_wave = binned_spectra['WAVE'].data.copy()

        templates = template_library['FLUX'].data.copy()

        velscale = spectrum_velocity_scale(obj_wave)
        if not PPXFFit._obj_tpl_pixelmatch(velscale, tpl_wave, velscale_ratio=velscale_ratio,
                                           dvtol=dvtol):
            raise ValueError('Pixel scale of the object and template spectra must be identical.')
        _velscale_ratio = 1 if velscale_ratio is None else velscale_ratio

        obj_flux = binned_spectra.copy_to_masked_array(flag=binned_spectra.do_not_fit_flags())

        outRange = obj_wave[ [0,-1] ]
        outpix = binned_spectra.nwave

        models = numpy.ma.zeros((binned_spectra.nbins, binned_spectra.nwave), dtype=numpy.float)

        for i in range(binned_spectra.nbins):
            if i in binned_spectra.missing_bins:
                models[i,:] = numpy.ma.masked
                continue
            if model_par['MASK'][i] > 0:
                models[i,:] = numpy.ma.masked
                warnings.warn('Fit to spectrum {0} masked.'.format(i+1))
                continue

            models[i,:] = PPXFFit.reconstruct_model(tpl_wave, templates, obj_wave,
                                                   model_par['KIN'][i,:], model_par['TPLWGT'][i],
                                                   velscale, polyweights=model_par['ADDCOEF'][i],
                                                   mpolyweights=model_par['MULTCOEF'][i],
                                                   start=model_par['BEGPIX'][i],
                                                   end=model_par['ENDPIX'][i],
                                                   redshift_only=redshift_only,
                                                   velscale_ratio=_velscale_ratio, dvtol=dvtol)

        return models


#    @staticmethod
#    def construct_models(binned_spectra, template_library, model_par, degree, mdegree, moments,
#                         redshift_only=False, velscale_ratio=None, dvtol=1e-10):
#        """
#        Construct models using the provided set of model parameters and
#        template library.  The number of template weights in the
#        provided model parameters must match the number of template
#        spectra.
#
#        If redshift_only is true, ignore the LOSVD and simply shift the
#        model to the correct redshift.
#
#        CURRENTLY WILL NOT WORK WITH:
#            - multiple kinematic components
#            - point-symmetric fits
#
#        .. todo::
#
#            This doesn't appear to give exactly the correct
#            reconstruction of the pPXF models.
#
#        """
#        if not isinstance(binned_spectra, SpatiallyBinnedSpectra):
#            raise TypeError('Must provide a SpatiallyBinnedSpectra object.')
#        if not isinstance(template_library, TemplateLibrary):
#            raise TypeError('Must provide a TemplateLibrary object.')
#        if model_par['TPLWGT'].shape[1] != template_library.ntpl:
#            raise ValueError('The number of weights does not match the number of templates.')
#
#        tpl_wave = template_library['WAVE'].data
#        obj_wave = binned_spectra['WAVE'].data
#
#        velscale = spectrum_velocity_scale(obj_wave)
#        if not PPXFFit._obj_tpl_pixelmatch(velscale, tpl_wave, velscale_ratio=velscale_ratio,
#                                           dvtol=dvtol):
#            raise ValueError('Pixel scale of the object and template spectra must be identical.')
#        _velscale_ratio = 1 if velscale_ratio is None else velscale_ratio
#
#        obj_flux = binned_spectra.copy_to_masked_array(flag=binned_spectra.do_not_fit_flags())
#
#        outRange = obj_wave[ [0,-1] ]
#        outpix = binned_spectra.nwave
#
#        models = numpy.ma.zeros((binned_spectra.nbins, binned_spectra.nwave), dtype=numpy.float)
#
#        for i in range(binned_spectra.nbins):
#            if i in binned_spectra.missing_bins:
#                models[i,:] = numpy.ma.masked
#                continue
#            if model_par['MASK'][i] > 0:
#                models[i,:] = numpy.ma.masked
#                warnings.warn('Fit to spectrum {0} masked.'.format(i+1))
#                continue
#
#            # Get the composite template
#            indx = model_par['TPLWGT'][i] > 0
#            redshift = model_par['KIN'][i,0]/astropy.constants.c.to('km/s').value
#            composite_template = numpy.dot(model_par['TPLWGT'][i][indx],
#                                           template_library['FLUX'].data[indx,:])
#
##            pyplot.step(tpl_wave, composite_template, where='mid', linestyle='-', lw=0.5,
##                         color='k')
##            pyplot.show()
#
#            start = model_par['BEGPIX'][i]
#            end = model_par['ENDPIX'][i]
#
#            # Redshift and broadened using the LOSVD parameters
#            if redshift_only:
#                # Resample the redshifted template to the wavelength grid of
#                # the binned spectra
#                models[i,:] = PPXFFit._shift_via_resample(tpl_wave, composite_template, redshift,
#                                                          True, outRange, outpix, True)
#            else:
##                print('{0}/{1}'.format(i+1,binned_spectra.nbins))
##                print('Start,End: ', start, end)
##                print(obj_wave[start:end].shape)
##                print(tpl_wave.shape)
#
#                # for PPXF v6.0.0
#                vsyst = -PPXFFit.ppxf_tpl_obj_voff(tpl_wave, obj_wave[start:end], velscale)
#                if _velscale_ratio > 1:
#                    npix_temp = composite_template.size - composite_template.size % _velscale_ratio
#                    _composite_template = composite_template[:npix_temp].reshape(-1,1)
#                    npix_temp //= _velscale_ratio
#
#                ctmp_rfft, npad = _templates_rfft(_composite_template)
#                _moments = numpy.atleast_1d(moments)
#                ncomp = 1
#                if _moments.size == 1:
#                    _moments = numpy.full(ncomp, numpy.absolute(_moments), dtype=int)
#                else:
#                    _moments = numpy.absolute(_moments)
#                vj = numpy.append(0, numpy.cumsum(_moments)[:-1])
#                par = model_par['KIN'][i,:].copy()
#                # Convert the velocity to pixel units
#                par[0], verr = PPXFFit._revert_velocity(par[0], 1.0)
#                # Convert the velocity dispersion to ignore the
#                # resolution difference
#                if model_par['SIGMACORR'][i] > 0:
#                    par[1] = numpy.square(par[1]) - numpy.square(model_par['SIGMACORR'][i])
#                    par[1] = numpy.sqrt(par[1]) if par[1] > 0 else 0.0
#
#                # Construct the model
#                if par[1] > 0:
#                    # Velocity dispersion is valid, follow the pPXF
#                    # method
#                    par /= velscale
#                    kern_rfft = _losvd_rfft(par, 1, _moments, vj, npad, 1, vsyst/velscale,
#                                            _velscale_ratio, 0.0)
#                    _model = numpy.fft.irfft(ctmp_rfft[:,0]
#                                             * kern_rfft[:,0,0]).ravel()[:npix_temp*_velscale_ratio]
#                    if _velscale_ratio > 1:
#                        _model = numpy.mean(_model.reshape(-1,_velscale_ratio), axis=1).ravel()
##                    pyplot.plot(obj_wave[start:end], _model[:end-start])
##                    pyplot.show()
#                    models[i,start:end] = _model[:end-start]
#                else:
#                    # Velocity dispersion is undefined, return a
#                    # zero-dispersion spectrum
#                    models[i,:] = PPXFFit._shift_via_resample(tpl_wave, composite_template,
#                                                              redshift, True, outRange, outpix,
#                                                              True)
#
##            pyplot.plot(obj_wave[start:end], obj_flux[i,start:end])
##            pyplot.plot(tpl_wave, composite_template)
##            pyplot.plot(obj_wave[start:end], models[i,start:end])
##            pyplot.show()
##            exit()
#
#            # Account for the polynomial
#            x = numpy.linspace(-1, 1, model_par['NPIXTOT'][i])
#            if mdegree > 0:
#                models[i,start:end] *= numpy.polynomial.legendre.legval(x,
#                                                        numpy.append(1.0,model_par['MULTCOEF'][i]))
#            if degree > -1:
#                models[i,start:end] += numpy.polynomial.legendre.legval(x, model_par['ADDCOEF'][i])
#
##            pyplot.plot(numpy.arange(end-start), obj_flux[i,start:end])
##            pyplot.plot(numpy.arange(end-start), models[i,start:end])
##            pyplot.show()
#
#        return models
#



