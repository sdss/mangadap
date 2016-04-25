# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
A class hierarchy that performs spectral fits.

*License*:
    Copyright (c) 2015, SDSS-IV/MaNGA Pipeline Group
        Licensed under BSD 3-clause license - see LICENSE.rst

*Source location*:
    $MANGADAP_DIR/python/mangadap/proc/spectralfitting.py

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
        import astropy.constants

        from ..par.parset import ParSet
        from ..util.bitmask import BitMask
        from ..util.fileio import init_record_array
        from .templatelibrary import TemplateLibrary
        from .pixelmask import PixelMask
        from ..contrib.ppxf import ppxf
        from ..util.instrument import spectrum_velocity_scale

*Class usage examples*:
        Add examples

*Revision history*:
    | **14 Apr 2016**: Implementation begun by K. Westfall (KBW)
    | **19 Apr 2016**: (KBW) First version

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

from matplotlib import pyplot

__author__ = 'Kyle B. Westfall'
# Add strict versioning
# from distutils.version import StrictVersion



# BASE CLASS -----------------------------------------------------------

class SpectralFitting():
    """
    Base class for spectral fitting.
    """
    def __init__(self, fit_type, bitmask, par=None):
        self.fit_type = fit_type
        if not isinstance(bitmask, BitMask):
            raise TypeError('Input bit mask must have type BitMask.')
        self.bitmask = bitmask
        self.par = par

# ----------------------------------------------------------------------

class StellarKinematicsFit(SpectralFitting):
    """
    Base class for fitting stellar kinematics.
    """
    def __init__(self, fit_method, bitmask, par=None):
        SpectralFitting.__init__(self, 'stellar_kinematics', bitmask, par=par)
        self.fit_method = fit_method

    @staticmethod
    def _per_stellar_kinematics_dtype(ntpl, nadd, nmult, nkin):
        r"""
        Construct the record array data type for the output fits
        extension.
        """
        return [ ('BIN_INDEX',numpy.int),
                 ('BEGPIX', numpy.int),
                 ('ENDPIX', numpy.int),
                 ('NPIXTOT',numpy.int),
                 ('NPIXFIT',numpy.int),
                 ('TPLWGT',numpy.float,ntpl),
                 ('ADDCOEF',numpy.float,nadd) if nadd > 1 else ('ADDCOEF',numpy.float),
                 ('MULTCOEF',numpy.float,nmult) if nmult > 1 else ('MULTCOEF',numpy.float),
                 ('KININP',numpy.float,2),
                 ('KIN',numpy.float,nkin),
                 ('KINERR',numpy.float,nkin),
                 ('CHI2',numpy.float),
                 ('RCHI2',numpy.float),
                 ('ROBUST_RCHI2',numpy.float),
                 ('RMS',numpy.float),
                 ('RESID',numpy.float,7),
                 ('FRAC_RMS',numpy.float),
                 ('FRAC_RESID',numpy.float,7)
               ]

class CompositionFit(SpectralFitting):
    """
    Base class for fitting the spectral composition.
    """
    def __init__(self, fit_method, bitmask, par=None):
        SpectralFitting.__init__(self, 'composition', bitmask, par=par)
        self.fit_method = fit_method


class EmissionLineFit(SpectralFitting):
    """
    Base class for fitting emission lines.
    """
    def __init__(self, fit_method, bitmask, par=None):
        SpectralFitting.__init__(self, 'emission_line', bitmask, par=par)
        self.fit_method = fit_method

# ----------------------------------------------------------------------

class PPXFFitPar(ParSet):
    def __init__(self, template_library_key, template_library, guess_redshift, guess_dispersion,
                 minimum_snr=None, pixelmask=None, bias=None, clean=False, degree=None,
                 mdegree=None, moments=None, oversample=None, sky=None, regul=None, reddening=None,
                 component=None, reg_dim=None):
    
        arr = [ numpy.ndarray, list ] # sky, p0, lam, component
        arr_in_fl = [ numpy.ndarray, list, int, float ] # component, reg_dim
        in_fl = [ int, float ] # Reddening, bias, regul, vsyst

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
                     int, int, bool, arr, in_fl, in_fl, arr_in_fl, arr_in_fl ]

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
                 'oversample':False, 'regul':0, 'reddening':None, 'component':0, 'reg_dim':None }


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
    """
    def __init__(self, bitmask, par=None):
        StellarKinematicsFit.__init__(self, 'ppxf', bitmask, par=par)
        self.snr_flag = 'LOW_SNR'
        self.tpl_flag = 'TPL_PIXELS'
        self.trunc_flag = 'TRUNCATED'
        self.rej_flag = 'PPXF_REJECT'


    def _fitting_mask(self, tpl_wave, obj_wave, velScale, velocity_offset, waverange=None,
                      max_velocity_range=400., alias_window=2400.):
        """
        Return a list of pixels in the object spectrum to be fit using pPXF.
    
        The limits applied to the fitted pixels are:
    
            - Apply the provided wavelength range limit
            (*wave_range_analysis*).
    
            - pPXF will only allow a fit when the number of template pixels
            is the same as or exceeds the number of pixels in the object
            spectra.  The first step toward limiting template spectra that
            are too long is to truncate the blue and red edges that likely
            won't be used given the provided velocity offsets
            (*velocity_offset*) and the expected velocity range
            (*max_velocity_range*).
    
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
            array: Vector with the indices of pixels in the object spectrum
            to fit using pPXF.

        """

        # 1. Apply the wavelength range limit, if provided
        now=len(obj_wave)                               # Number of object wavelengths
        print('Original number of object pixels: {0}'.format(now))
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
        print('After selecting the wavelength range to analyze: {0}'.format(now))
        print('Number of template pixels: {0}'.format(ntw))
        if ntw < now:
            # Indices of wavelengths redward of the redshifted template
            indx = obj_wave > (tpl_wave[0]*(1. + z_min))
            if numpy.sum(indx) == 0:
                raise ValueError('No overlapping wavelengths between galaxy and template!')
            print('Pixels blueward of redshifted template: {0}'.format(
                                                                    len(obj_wave)-numpy.sum(indx)))
            # Merge with current index
            fit_indx &= indx
            if numpy.sum(fit_indx) == 0:
                raise ValueError('No overlapping wavelengths between galaxy and template!')
            now_= numpy.sum(fit_indx)
            print('After merging with the specified wavelength range: {0}'.format(now_))
    
        # New number of good object pixels
        if ntw < now_:
            fit_indx[ numpy.where(fit_indx)[0][ntw:] ] = False      # Truncate the red edge as well
            print('Impose same as number of template pixels: {0} {1}'.format(numpy.sum(fit_indx),
                                                                             ntw))
        npix_mask = ~(fit_indx) & ~(waverange_mask)

        # 3. Limit wavelength range to avoid aliasing problems in the template convolution
        nalias = int(numpy.floor(alias_window/velScale))            # Number of pixels to mask
        # Mask to the range that should be unaffected by alias errors
        waverange_tpl = [ tpl_wave[nalias]*(1+z_max), tpl_wave[ntw-nalias-1]*(1+z_min) ]
        print('Mask to these wavelengths to avoid convolution aliasing: {0} - {1}'.format(
                                                            waverange_tpl[0], waverange_tpl[1]))
        indx = numpy.logical_and(obj_wave > waverange_tpl[0], obj_wave < waverange_tpl[1])

        # Merge with current index
        fit_indx &= indx
        if numpy.sum(fit_indx) == 0:
            raise ValueError('No wavelengths available in this range!')

        print('Final wavelength range to fit: {0} {1}'.format(obj_wave[fit_indx][0],
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
        return (numpy.log(tpl_wave[0])-numpy.log(obj_wave[0]))*velScale \
                    / numpy.diff(numpy.log(obj_wave[0:2]))[0]


    def fit_SpatiallyBinnedSpectra(self, binned_spectra, par=None):
        if par is not None:
            self.par = par
        if self.par is None:
            raise ValueError('Required parameters for PPXFFit have not been defined.')

        if not isinstance(binned_spectra, SpatiallyBinnedSpectra):
            raise TypeError('Must provide a SpatiallyBinnedSpectra object for fitting.')
        if self.par['template_library'] is None \
                or not isinstance(self.par['template_library'], TemplateLibrary):
            raise TypeError('Must provide a TemplateLibrary object for fitting.')

        _def = self.par._keyword_defaults()
        if self.par['regul'] != _def['regul'] \
                or self.par['reddening'] != _def['reddening'] \
                or self.par['component'] != _def['component'] \
                or self.par['reg_dim'] != _def['reg_dim']:
            raise NotImplementedError('Cannot use regul, reddening, component, or regul_dim yet.')

        # Get the spectra
        obj_wave = binned_spectra['WAVE'].data
        velScale = spectrum_velocity_scale(obj_wave)
        flux = binned_spectra.copy_to_masked_array(flag=binned_spectra.do_not_fit_flags())
        noise = numpy.ma.sqrt(numpy.ma.power(binned_spectra.copy_to_masked_array(ext='IVAR',
                                                    flag=binned_spectra.do_not_fit_flags()), -1.0))
        # Get the input guess kinematics
        guess_kin = numpy.empty((binned_spectra.nbins, 2), dtype=numpy.float)
        guess_kin[:,0] = self.par['guess_redshift']*astropy.constants.c.to('km/s').value
        guess_kin[:,1] = self.par['guess_dispersion']

        # Get the template data
        ntpl = self.par['template_library'].ntpl
        tpl_wave = self.par['template_library']['WAVE'].data
        templates = self.par['template_library']['FLUX'].data

        # Initialize the output arrays
        #  Model flux:
        model_flux = numpy.zeros(flux.shape, dtype=numpy.float)
        #  Model mask:
        model_mask = self.par['pixelmask'].bits(self.bitmask, obj_wave, nspec=binned_spectra.nbins,
                                                velocity_offsets=guess_kin[:,0])
        model_mask[numpy.ma.getmaskarray(flux)] \
                = self.bitmask.turn_on(model_mask[numpy.ma.getmaskarray(flux)], 'DIDNOTUSE')
        #  Record array with model parameters
        model_par = init_record_array(binned_spectra.nbins,
                                      self._per_stellar_kinematics_dtype(ntpl, self.par['degree']+1,
                                                                         max(self.par['mdegree'],0),
                                                                         self.par['moments']))
        model_par['BIN_INDEX'] = numpy.arange(binned_spectra.nbins)

        # TODO: Need parameter keywords for max_velocity_range and
        # alias_window
        try:
            fit_indx, waverange_mask, npix_mask, alias_mask \
                = self._fitting_mask(tpl_wave, obj_wave, velScale, guess_kin[:,0],
                                     waverange=self.par['pixelmask'].waverange)
        except ValueError as e:
            warnings.warn('No fitting done because of a masking error: {0}'.format(e))
            return obj_wave, model_flux, model_mask, model_par

        # Add to the mask
        model_mask[:,npix_mask] = self.bitmask.turn_on(model_mask[:,npix_mask], self.tpl_flag)
        model_mask[:,alias_mask] = self.bitmask.turn_on(model_mask[:,alias_mask], self.trunc_flag)
        # And update the internal mask of the data
        flux[model_mask > 0] = numpy.ma.masked
        noise[model_mask > 0] = numpy.ma.masked

        # Determine the starting and ending pixels
        start, end = numpy.where(fit_indx)[0][ [0,-1] ]
        end += 1
        print(start, end, len(obj_wave))

        # Get the input pixel shift between the object and template
        # wavelength vectors; interpretted by pPXF as a base velocity
        # shift between the two
        base_velocity = self.ppxf_tpl_obj_voff(tpl_wave, obj_wave[start:end], velScale)
        guess_kin[:,0] += base_velocity
        vsyst = numpy.mean(guess_kin[:,0])

        # Determine the degrees of freedom of the fit, used for the
        # brute force calculation of the reduced chi-square
        dof = self.par['moments']+max(self.par['mdegree'],0)
        if self.par['degree'] >= 0:
            dof += self.par['degree']+1

        # Fit each binned spectrum:
        for i in range(binned_spectra.nbins):
            print('Fit: {0}/{1}'.format(i+1,binned_spectra.nbins), end='\r')
            # Get the list of good pixels
            gpm = numpy.where( ~(flux.mask[i,start:end]))[0]
            if len(gpm) == 0:
                warnings.warn('No unmasked pixels in bin {0}.'.format(i+1))
                continue

            # Check that it meets the S/N limit
            if binned_spectra['BINS'].data['SNR'][i] < self.par['minimum_snr']:
                model_mask[i,:] = self.bitmask.turn_on(model_mask[i,:], 'LOW_SNR')
                warnings.warn('Bin {0} below S/N limit.' \
                              'S/N={1:.3f}; minimum to fit={2:.1f}.'.format(i+1,
                                binned_spectra['BINS'].data['SNR'][i], self.par['minimum_snr']))
                continue

#            pyplot.plot( obj_wave[start:end], flux[i,start:end]/numpy.ma.mean(flux[i,start:end]))
#            pyplot.plot( tpl_wave*(1+(vsyst-base_velocity)/astropy.constants.c.to('km/s').value),
#                            templates[31,:]/numpy.mean(templates[31,:]))
#            pyplot.show()

            # Run ppxf
            ppxf_fit = ppxf(templates.T, flux.data[i,start:end], noise.data[i,start:end],
                            velScale, guess_kin[i,:], goodpixels=gpm, bias=self.par['bias'],
                            clean=self.par['clean'], degree=self.par['degree'],
                            mdegree=self.par['mdegree'], moments=self.par['moments'],
                            oversample=self.par['oversample'], vsyst=vsyst, quiet=True)#, plot=True)
                            #, regul=par['regul'], lam=obj_wave[fit_indx],
                            #reddening=par['reddening'],
                            #component=par['component'], reg_dim=par['reg_dim'])

            # Save the result
            model_flux[i,start:end] = ppxf_fit.bestfit
            if len(gpm) != len(ppxf_fit.goodpixels):
                rejected_pixels = list(set(gpm)-set(ppxf_fit.goodpixels))
                self.bitmask.turn_on(model_mask[i,start:end][rejected_pixels], flag=self.rej_flag)
            # TODO: Turn on bits that flag as being rejected by pPXF,
            # failed pPXF, etc.
            model_par['BEGPIX'][i] = start
            model_par['ENDPIX'][i] = end
            model_par['NPIXTOT'][i] = end-start
            model_par['NPIXFIT'][i] = len(ppxf_fit.goodpixels)
#            print(model_par['NPIXTOT'][i], model_par['NPIXFIT'][i])
            # TODO: Store sparsely?
            model_par['TPLWGT'][i] = ppxf_fit.weights[0:ntpl]
#            print(model_par['TPLWGT'][i])
            if par['degree'] >= 0:
                model_par['ADDCOEF'][i] = ppxf_fit.polyweights
            if par['mdegree'] > 0:
                model_par['MULTCOEF'][i] = ppxf_fit.mpolyweights
#            print(model_par['MULTCOEF'][i])
            model_par['KININP'][i] = guess_kin[i,:]
            model_par['KININP'][i][0] -= base_velocity
#            print(model_par['KININP'][i])
            model_par['KIN'][i] = ppxf_fit.sol
            model_par['KIN'][i][0] += vsyst - base_velocity
#            print(model_par['KIN'][i])
            model_par['KINERR'][i] = ppxf_fit.error
#            print(model_par['KINERR'][i])
            resid = flux[i,start:end] - ppxf_fit.bestfit
            model_par['CHI2'][i] = numpy.sum( numpy.square(
                                                (resid/noise[i,start:end])[ppxf_fit.goodpixels] ))
#            print(model_par['CHI2'][i])
            model_par['RCHI2'][i] = model_par['CHI2'][i] \
                            / (model_par['NPIXFIT'][i] - dof-numpy.sum(model_par['TPLWGT'][i] > 0))
#            print(model_par['RCHI2'][i])
            model_par['ROBUST_RCHI2'][i] = ppxf_fit.chi2
#            print(model_par['ROBUST_RCHI2'][i])
            model_par['RMS'][i] = numpy.sqrt(numpy.ma.mean( numpy.square(
                                                                resid[ppxf_fit.goodpixels]) ))
#            print(model_par['RMS'][i])
            grw = numpy.arange(model_par['NPIXFIT'][i]).astype(numpy.float)/model_par['NPIXFIT'][i]
            resid_sort = numpy.sort(numpy.absolute(resid[ppxf_fit.goodpixels]))
            interp = scipy.interpolate.interp1d(grw, resid_sort, fill_value='extrapolate')
            model_par['RESID'][i] = numpy.append(numpy.append(resid_sort[0],
                                                            interp([0.25, 0.50, 0.75, 0.90, 0.99])),
                                                              resid_sort[-1])
            indx = numpy.absolute(ppxf_fit.bestfit) > 0
            frac_resid = resid
            frac_resid[indx] /= ppxf_fit.bestfit[indx]
            _goodpixels = numpy.intersect1d(numpy.arange(len(resid))[indx], ppxf_fit.goodpixels,
                                            assume_unique=True)
            model_par['FRAC_RMS'][i] = numpy.sqrt(numpy.ma.mean(numpy.square(
                                                                frac_resid[_goodpixels])))
#            print(model_par['RMS'][i])
            ngood = len(_goodpixels)
            print(ngood, model_par['NPIXFIT'][i])
            grw = numpy.arange(ngood).astype(numpy.float)/ngood
            resid_sort = numpy.sort(numpy.absolute(resid[_goodpixels]))
            interp = scipy.interpolate.interp1d(grw, resid_sort, fill_value='extrapolate')
            model_par['FRAC_RESID'][i] = numpy.append(numpy.append(resid_sort[0],
                                                      interp([0.25, 0.50, 0.75, 0.90, 0.99])),
                                                      resid_sort[-1])
            print(model_par['FRAC_RESID'][i])

#            pyplot.step(obj_wave, flux[i,:], where='mid', linestyle='-', lw=0.5, color='k',
#                        zorder=3)
#            pyplot.plot(obj_wave, model_flux[i,:], linestyle='-', lw=1.5, color='r',
#                        zorder=1, alpha=0.5)
#            pyplot.show()
#            exit()

        model_par['KIN'][:,0], model_par['KINERR'][:,0] \
                = self._convert_velocity(model_par['KIN'][:,0], model_par['KINERR'][:,0])

        print('Fit: {0}/{0}'.format(binned_spectra.nbins))
        return binned_spectra['WAVE'].data, model_flux, model_mask, model_par


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

        


    def construct_models(self, binned_spectra, template_library, model_par, degree, mdegree,
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

        models = numpy.zeros((binned_spectra.nbins, binned_spectra.nwave), dtype=numpy.float)

        for i in range(binned_spectra.nbins):
            if i in binned_spectra.missing_bins:
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
                base_velocity = self.ppxf_tpl_obj_voff(tpl_wave, obj_wave[start:end], velScale)
                vel, losvd = self._losvd_kernel(model_par['KIN'][i,:], velScale,
                                                base_velocity=base_velocity)
                _composite_template = scipy.signal.fftconvolve(composite_template, losvd,
                                                               mode="same")
                voff, verr = self._convert_velocity(-base_velocity, 1.0)
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














