#!/usr/bin/env python3
# An excerpt from Kyle's script fit_one_spec.py
import os
import time
import numpy

from astropy.io import fits
from matplotlib import pyplot

from mangadap.drpfits import DRPFits
from mangadap.util.fitsutil import DAPFitsUtil
from mangadap.util.resolution import SpectralResolution
from mangadap.util.pixelmask import SpectralPixelMask

from mangadap.par.artifactdb import ArtifactDB
from mangadap.par.emissionmomentsdb import EmissionMomentsDB
from mangadap.par.emissionlinedb import EmissionLineDB

from mangadap.proc.templatelibrary import TemplateLibrary
from mangadap.proc.emissionlinemoments import EmissionLineMoments
from mangadap.proc.sasuke import Sasuke
from mangadap.proc.ppxffit import PPXFFit
from mangadap.proc.stellarcontinuummodel import StellarContinuumModel, StellarContinuumModelBitMask
from mangadap.proc.emissionlinemodel import EmissionLineModelBitMask

#-----------------------------------------------------------------------------

def get_redshift(plt, ifu): #, drpall_file):
    hdu = fits.open(os.path.join(os.environ['MANGA_SPECTRO_REDUX'], os.environ['MANGADRP_VER'],
                                 'drpall-{0}.fits'.format(os.environ['MANGADRP_VER'])))
#    hdu = fits.open(drpall_file)
    indx = hdu[1].data['PLATEIFU'] == '{0}-{1}'.format(plt, ifu)
    return hdu[1].data['NSA_Z'][indx][0]

#
# def get_spectrum(plt, ifu, x, y, directory_path=None):
#     drpf = DRPFits(plt, ifu, 'CUBE', read=True, directory_path=directory_path)
#     flat_indx = drpf.spatial_shape[1]*x+y
#     # This function always returns as masked array
#     sres = drpf.spectral_resolution(toarray=True, fill=True, pre=True).data
#     flux = drpf.copy_to_masked_array(ext='FLUX', flag=drpf.do_not_fit_flags())
#     ivar = drpf.copy_to_masked_array(ext='IVAR', flag=drpf.do_not_fit_flags())
#     return drpf['WAVE'].data, flux[flat_indx,:], ivar[flat_indx,:], sres[flat_indx,:]


#-----------------------------------------------------------------------------

def fit_one_spect(plt, ifu, flux, ivar, wave):
    if __name__ == '__main__':
        t = time.clock()

        # Show the ppxf plots
    #    fit_plots = True
        fit_plots = False

        # Template keywords
        sc_tpl_key = 'MILES_SSP'
        el_tpl_key = 'MILES_SSP'
    #    el_tpl_key = 'BC03', 'MILESHC'

        # Emission-line database keywords
        elmom_key = 'ELBMILES'
        elfit_key = 'ELPMILES'
    #    elmom_key = 'ELBTT'
    #    elfit_key = 'ELPT4'

        # Template pixel scale a factor of 4 smaller than galaxy data 70
        velscale_ratio = 4 #??
        velscale_ratio = 70
        # Get the redshift
        #drpall_file = os.path.join(os.environ['MANGA_SPECTRO_REDUX'], 'drpall-v2_3_1.fits')
        #z = numpy.array([get_redshift(plt, ifu)]) #, drpall_file)])
        z=0
        #print('Redshift: {0}'.format(z[0]))
        dispersion = numpy.array([300.])
        sres = drpf.spectral_resolution(toarray=True, fill=True, pre=True).data
        sres = sres[1000,:] # Blantant fudge as sres is 4096 x x4593, and needs to be sam shape as flux vector
        # Fitting functions expect data to be in 2D arrays (for now):
        flux = flux.reshape(1,-1)
        ferr = numpy.ma.power(ivar, -0.5).reshape(1,-1)
        #sres = sres.reshape(1,-1)

        #-------------------------------------------------------------------
        # Fit the stellar continuum

        # Mask the 5577 sky line and the emission lines
        sc_pixel_mask = SpectralPixelMask(artdb=ArtifactDB('BADSKY'), emldb=EmissionLineDB('ELPFULL'))

        # Construct the template library
        sc_tpl = TemplateLibrary(sc_tpl_key, match_to_drp_resolution=False,
                                 velscale_ratio=velscale_ratio, spectral_step=1e-4, log=True,
                                 hardcopy=False)
        sc_tpl_sres = numpy.mean(sc_tpl['SPECRES'].data, axis=0).ravel()

        # Instantiate the fitting class
        ppxf = PPXFFit(StellarContinuumModelBitMask())

        # Perform the fit
        cont_wave, cont_flux, cont_mask, cont_par \
            = ppxf.fit(sc_tpl['WAVE'].data.copy(), sc_tpl['FLUX'].data.copy(), wave, flux, ferr,
                       z, dispersion, iteration_mode='no_global_wrej', reject_boxcar=100,
                       ensemble=False, velscale_ratio=velscale_ratio, mask=sc_pixel_mask,
                       matched_resolution=False, tpl_sres=sc_tpl_sres, obj_sres=sres, degree=-1,
                       moments=2, plot=fit_plots) #Degree was 8
        return cont_wave, cont_flux, cont_mask, cont_par
