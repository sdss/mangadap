#!/usr/bin/env python3
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from __future__ import unicode_literals

import os
import time
import numpy as np
import pdb
from astropy import constants as const
from astropy import table
import matplotlib.pyplot as plt
from astropy.io import fits
from argparse import ArgumentParser
from scipy import interpolate
import pdb

#from mangadap.drpfits import DRPFits
#from mangadap.util.fitsutil import DAPFitsUtil
from mangadap.util.resolution import SpectralResolution
from mangadap.util.pixelmask import SpectralPixelMask

from mangadap.par.artifactdb import ArtifactDB
#from mangadap.par.emissionmomentsdb import EmissionMomentsDB
from mangadap.par.emissionlinedb import EmissionLineDB

from mangadap.proc.templatelibrary import TemplateLibrary
#from mangadap.proc.emissionlinemoments import EmissionLineMoments
from mangadap.proc.sasuke import Sasuke
from mangadap.proc.ppxffit import PPXFFit

from mangadap.proc.stellarcontinuummodel import StellarContinuumModel, StellarContinuumModelBitMask
from mangadap.proc.emissionlinemodel import EmissionLineModelBitMask
#%%

plt.ioff()


class fit_single_spectrum(object):
    
    def __init__(self, flux, ivar, wave, sres, redshift, sc_tpl_key='MILESHC', 
                 el_tpl_key = 'MILESHC',elmom_key = 'ELBMILES', elfit_key = 'ELPMILES',
                 plot=True):
    
        self.flux=flux
        self.wave=wave
        self.sres=sres
        self.redshift=redshift
        
        self.sc_tpl_key=sc_tpl_key
        self.el_tpl_key=el_tpl_key
        self.elmom_key=elmom_key
        self.elfit_key=elfit_key
        self.plot=plot
        
        spectral_step= np.mean(np.diff(np.log10(self.wave)))

        # Template pixel scale a factor of 4 smaller than galaxy data
        # note velscale MUST be an integer
        velscale_ratio = 4
        print( 'velscale', np.log(self.wave[1]/self.wave[0])*const.c.to('km/s'))
    
        # Get the redshift
        z = [self.redshift]
        print('Redshift: {0}'.format(z[0]))
        dispersion = np.array([100.])
    
        # Fitting functions expect data to be in 2D arrays (for now):
        self.flux = self.flux.reshape(1,-1)
        self.ferr = np.ma.power(ivar, -0.5).reshape(1,-1)
        self.sres = self.sres.reshape(1,-1)
    
        #-------------------------------------------------------------------
        # Fit the stellar continuum
    
        # Mask the 5577 sky line and the emission lines
        sc_pixel_mask = SpectralPixelMask(artdb=ArtifactDB('BADSKY'),
                                          emldb=EmissionLineDB.from_key('ELPFULL'))
    
        # Construct the template library
        sc_tpl = TemplateLibrary(sc_tpl_key, match_to_drp_resolution=False,
                velscale_ratio=velscale_ratio, spectral_step=spectral_step, log=True,
                hardcopy=False)
        sc_tpl_sres = np.mean(sc_tpl['SPECRES'].data, axis=0).ravel()
    
        # Instantiate the fitting class
        ppxf = PPXFFit(StellarContinuumModelBitMask())
   
        # Perform the fit
        cont_wave, cont_flux, cont_mask, cont_par \
            = ppxf.fit(sc_tpl['WAVE'].data.copy(), \
                       sc_tpl['FLUX'].data.copy(), self.wave, self.flux, self.ferr,
                       z, dispersion, iteration_mode='no_global_wrej', reject_boxcar=100,
                       ensemble=False, velscale_ratio=velscale_ratio, mask=sc_pixel_mask,
                       matched_resolution=False, tpl_sres=sc_tpl_sres, 
                       obj_sres=self.sres, degree=8, mdegree=0,
                       moments=2, plot=False)
    
        # Remask the continuum fit
        sc_continuum = StellarContinuumModel.reset_continuum_mask_window(
                            np.ma.MaskedArray(cont_flux, mask=cont_mask>0))
    
        # Show the fit and residual
#        if self.plot ==True:
#            plt.figure()
#            plt.plot(self.wave, self.flux[0,:], label='Data')
#            plt.plot(self.wave, sc_continuum[0,:], label='Model')
#            plt.plot(self.wave, self.flux[0,:] - sc_continuum[0,:], label='Resid')
#            plt.legend()
#            plt.xlabel('Wavelength')
#            plt.ylabel('Flux')
#            plt.show()
        #%%
        #-------------------------------------------------------------------
    
        #-------------------------------------------------------------------
        # Get the emission-line moments using the fitted stellar continuum
    
        # Read the database that define the emission lines and passbands
#        momdb = EmissionMomentsDB(elmom_key)
    
        # Measure the moments
#        elmom = EmissionLineMoments.measure_moments(momdb, self.wave, self.flux, 
#                        continuum=sc_continuum,redshift=z)
        #-------------------------------------------------------------------
    
        #-------------------------------------------------------------------
        # Fit the emission-line model
    
        # Set the emission-line continuum templates if different from those
        # used for the stellar continuum
        if self.sc_tpl_key == self.el_tpl_key:
            # If the keywords are the same, just copy over the previous
            # library ...
            el_tpl = sc_tpl
            el_tpl_sres = sc_tpl_sres
            # ... and the best fitting stellar kinematics
            stellar_kinematics = cont_par['KIN']
        else:
            # If the template sets are different, we need to match the
            # spectral resolution to the galaxy data ...
            _sres = SpectralResolution(self.wave, self.sres[0,:], log10=True)
            el_tpl = TemplateLibrary(self.el_tpl_key, sres=_sres, 
                velscale_ratio=velscale_ratio,spectral_step=spectral_step, 
                log=True, hardcopy=False)
            el_tpl_sres = np.mean(el_tpl['SPECRES'].data, axis=0).ravel()
            # ... and use the corrected velocity dispersions.
            stellar_kinematics = cont_par['KIN']
            stellar_kinematics[:,1] = np.ma.sqrt(np.square(cont_par['KIN'][:,1]) -
                             np.square(cont_par['SIGMACORR_EMP']))
    
        # Mask the 5577 sky line
        el_pixel_mask = SpectralPixelMask(artdb=ArtifactDB('BADSKY'))
    
        # Read the emission line fitting database
        emldb = EmissionLineDB.from_key(self.elfit_key)
        
        # Instantiate the fitting class
        emlfit = Sasuke(EmissionLineModelBitMask())
    
        # Perform the fit
        eml_wave, model_flux, eml_flux, eml_mask, eml_fit_par, eml_eml_par \
            = emlfit.fit(emldb, self.wave, self.flux, obj_ferr=self.ferr, 
            obj_mask=el_pixel_mask, obj_sres=self.sres,
            guess_redshift=z, guess_dispersion=dispersion, reject_boxcar=101,
            stpl_wave=el_tpl['WAVE'].data, stpl_flux=el_tpl['FLUX'].data,
            stpl_sres=el_tpl_sres, stellar_kinematics=stellar_kinematics,
            etpl_sinst_mode='offset', etpl_sinst_min=10.,
            velscale_ratio=velscale_ratio, matched_resolution=False, mdegree=5,
            plot=False, sigma_rej=5.)
#
#        print(emlfit.velscale)
#        pdb.set_trace()
    
        # Get the stellar continuum that was fit for the emission lines
        elcmask = eml_mask.ravel() > 0
        goodpix = np.arange(elcmask.size)[np.invert(elcmask)]
        start, end = goodpix[0], goodpix[-1]+1
        elcmask[start:end] = False
        el_continuum = np.ma.MaskedArray(model_flux - eml_flux,
                                mask=elcmask.reshape(model_flux.shape))
        
        self.model_flux=model_flux[0,:]
        self.el_continuum=el_continuum[0,:]
        self.sc_continuum=sc_continuum[0,:]
        self.emlines_par=table.Table(eml_eml_par)[0]
        self.eml_fit_par=table.Table(eml_fit_par)[0]
        

        # Plot the result
        if self.plot ==True:
            plt.figure()
            plt.plot(self.wave, self.flux[0,:], label='Data')
            plt.plot(self.wave, model_flux[0,:], label='Model')
#            plt.plot(self.wave, el_continuum[0,:], label='EL Cont.')
#            plt.plot(self.wave, sc_continuum[0,:], label='SC Cont.')
#            plt.plot(self.wave, self.flux[0,:] - model_flux[0,:], label='Resid')
            plt.legend()
            plt.xlabel('Wavelength')
            plt.ylabel('Flux')
            plt.title('Emission line fitting and residuals')
            plt.ylim(np.min(self.flux[0,:] - model_flux[0,:]), None)
            
            plt.show()
        
#        pdb.set_trace()
        
        # Remeasure the emission-line moments with the new continuum
#        new_elmom = EmissionLineMoments.measure_moments(momdb, self.wave, self.flux, 
#                    continuum=el_continuum,redshift=z)
    
        # Compare the summed flux and Gaussian-fitted flux for all the
        # fitted lines
#        plt.figure()
#        plt.scatter(emldb['restwave'], (new_elmom['FLUX']-eml_eml_par['FLUX']).ravel(),
#                       c=eml_eml_par['FLUX'].ravel(), cmap='viridis', marker='.', s=60, lw=0, zorder=4)
#        plt.grid()
#        plt.xlabel('Wavelength')
#        plt.ylabel('Summed-Gaussian Difference')
#        plt.show()
    
       

if __name__ == '__main__':
    parser = ArgumentParser()

    parser.add_argument('test', type=str, help='spectrum for testing', default='manga')
    
    t = time.perf_counter()
    mangadap=os.environ['MANGADAP_DIR']
    
    arg = parser.parse_args()
    test=arg.test
    print('fitting galaxy', test)
    
    if test=='manga':
        
        
#        hdulist=fits.open(mangadap+'/data/testdata/manga-7992-1901-LOGCUBE.fits.gz')
        hdulist=fits.open('/Users/francesco/Downloads/manga-8313-1901-LOGCUBE.fits')
        
        wave=hdulist['WAVE'].data
        spectral_step=np.mean(np.diff(np.log10(wave)))
        flux=hdulist['FLUX'].data
        sres=hdulist['SPECRES'].data
        ivar=hdulist['IVAR'].data
        
        sz=flux.shape
        
        flux=flux[:,sz[1]//2,sz[2]//2]
        ivar=ivar[:,sz[1]//2,sz[2]//2]
        
#        redshift=0.0164
        redshift=0.0243
        
        
        out=fit_single_spectrum(flux, ivar, wave, sres, redshift, sc_tpl_key='MILESHC', 
                     el_tpl_key = 'MILESHC', elfit_key = 'ELPMILES',
                     plot=True)
#t0 - Everything free (even doublet fluxes)
#t1 - Tying the groups above in velocity but not dispersion; OIII, OI, NII doublet dispersions and fluxes are tied.
#t2 - Same as t1 but with the dispersion of each group tied.
#t3 - All velocities and dispersions tied; OIII, OI, NII doublet fluxes tied.
#         elfit_key = 'ELPMILES'
        
    elif test=='sdss':
#        hdulist=fits.open(mangadap+'/data/testdata/spec-0352-51694-0030.fits')
        hdulist=fits.open('/Users/francesco/Downloads/spec-1335-52824-0378.fits')
        
        
        coadd = hdulist['COADD'].data
        wave=np.array(10**coadd['loglam'], dtype='float64')
    #    restructure the wavelength vector so it's constant
        spectral_step=np.mean(np.diff(np.log10(wave)))
        wave=10**(np.linspace(np.log10(wave[0]), np.log10(wave[0])+\
                              spectral_step*(len(wave)-1), len(wave)))
        
        pix_ang=np.diff(wave)
        pix_ang = np.append(pix_ang, pix_ang[-1])
        sres=wave/ ( coadd['wdisp']*pix_ang*2.355)
        
        flux=coadd['flux']
        ivar=coadd['ivar']
#        mask=np.where(coadd['and_mask']>0, 1, 0)
        
        redshift=hdulist['SPECOBJ'].data['Z']
        out=fit_single_spectrum(flux, ivar, wave, sres, redshift, sc_tpl_key='MIUSCATTHIN', 
                     el_tpl_key = 'MIUSCATTHIN',elmom_key = 'ELBMILES', elfit_key = 'ELPMILES',
                     plot=True)
        
    elif test=='stack':
        
        hdu=fits.open(mangadap+'/data/testdata/0_2_STACKED_SPECTRA.fits')

        wave=hdu['WAVE'].data
        w = (wave>3700) & ( wave<9600)
        wave=wave[w]
        flux=hdu['STACKED_SPECTRA'].data[w, 2]
        ivar=(flux*0.01)**(-0.5)
        sres= flux*0.0+2000.
        redshift=0.
        
        out=fit_single_spectrum(flux, ivar, wave, sres, redshift, sc_tpl_key='MIUSCATTHIN', 
                     el_tpl_key = 'MIUSCATTHIN', elfit_key = 'ELPT2',
                     plot=True)
    elif test=='joe':
        hdu = fits.open('/Users/francesco/Desktop/spec1678.fits') 
#        sp=hdu['DATA'].data
        wave = hdu[1].data
        flux= hdu[0].data
        sres= flux*0.0+2000
        ivar = 1/(flux*0.01)**0.5
        
        ww  = (wave < 9000)
        wave=wave[ww]
        flux=flux[ww]
        sres=sres[ww]
        ivar=ivar[ww]
        
        redshift=0.5998
        
        
        out=fit_single_spectrum(flux, ivar, wave, sres, redshift, sc_tpl_key='MIUSCATTHIN', 
                     el_tpl_key = 'MIUSCATTHIN', elfit_key = 'ELPMILES',
                     plot=True)
        pdb.set_trace()

#    pdb.set_trace()
    print('Elapsed time: {0} seconds'.format(time.perf_counter() - t))
