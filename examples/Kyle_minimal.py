#Minimal working example that reproduces the error I get when trying to add the component and fraction keywords in to the pPXF call
#fit_bulge_spec is based on Kyle's fit_one_spec.py
#/Users/ppzaf/Work/MaNGA/dap/master/examples/Kyle_minimal.py
#All the imports, just in case
import os
import time
import numpy

from astropy.io import fits
from matplotlib import pyplot
import glob

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
from mangadap.proc.ppxffit_frac import PPXFFit_frac
from mangadap.proc.stellarcontinuummodel import StellarContinuumModel, StellarContinuumModelBitMask
from mangadap.proc.emissionlinemodel import EmissionLineModelBitMask

#Values for example spectrum in example galaxy
plt = '7443'
ifu = '9102'
frac = 0.8438
wave = numpy.genfromtxt('wave.txt')
flux = numpy.genfromtxt('flux.txt')
ivar = numpy.genfromtxt('ivar.txt')
sres = numpy.genfromtxt('sres.txt')
z = 0.0916

def fit_bulge_spec(plt, ifu, frac, wave, flux, ivar, sres):
    #Fit each bin using fraction of bulge light in each as 'component' and the single bulge and disk model spectra
    #I constructed
    fit_plots = True
    #flux = bulge_spec
    #ivar = bulge_ivar
    #wave = wave_b
    # Template keywords
    sc_tpl_key = 'MY_SPEC'
    el_tpl_key = 'MY_SPEC'
#    el_tpl_key = 'BC03', 'MILESHC'

    # Emission-line database keywords
    elmom_key = 'ELBMILES'
    elfit_key = 'ELPMILES'
#    elmom_key = 'ELBTT'
#    elfit_key = 'ELPT4'

    # Template pixel scale a factor of 4 smaller than galaxy data 70
    velscale_ratio = 4 #??
    # Get the redshift
    #drpall_file = os.path.join(os.environ['MANGA_SPECTRO_REDUX'], 'drpall-v2_3_1.fits')
    #z = numpy.array([get_redshift(plt, ifu)]) #, drpall_file)])
    #print('Redshift: {0}'.format(z[0]))
    dispersion = numpy.array([300.])
    #sres = drpf.spectral_resolution(toarray=True, fill=True, pre=True).data
    #sres = sres[1000,:] # Blantant fudge as sres is 4096 x x4593, and needs to be same shape as flux vector
    # Fitting functions expect data to be in 2D arrays (for now):
    flux = flux.reshape(1,-1)
    ferr = numpy.ma.power(ivar, -0.5).reshape(1,-1)

    #-------------------------------------------------------------------
    # Fit the stellar continuum

    # Mask the 5577 sky line and the emission lines
    sc_pixel_mask = SpectralPixelMask(artdb=ArtifactDB('BADSKY'), emldb=EmissionLineDB('ELPFULL'))

    # Construct the template library
    sc_tpl = TemplateLibrary(sc_tpl_key, match_to_drp_resolution=True,
                             velscale_ratio=velscale_ratio, spectral_step=1e-4, log=True,
                             hardcopy=False)

    sc_tpl_sres = numpy.mean(sc_tpl['SPECRES'].data, axis=0).ravel()

    # Instantiate the fitting class
    ppxf2 = PPXFFit_frac(StellarContinuumModelBitMask())

    # Perform the fit - this is where the error occurs!! moments has been changed from moments=2 to moments=[2,2], but start value has not been altered to include two components
    cont_wave, cont_flux, cont_mask, cont_par \
        = ppxf2.fit(sc_tpl['WAVE'].data.copy(), sc_tpl['FLUX'].data.copy(), wave, flux, ferr,
                   z, dispersion, frac, component = [0,1], iteration_mode='no_global_wrej', reject_boxcar=100,
                   ensemble=False, velscale_ratio=velscale_ratio, mask=sc_pixel_mask,
                   matched_resolution=False, tpl_sres=sc_tpl_sres, obj_sres=sres, degree=-1,
                   moments=[2,2], plot=fit_plots) #Degree was 8

    return cont_wave, cont_flux, cont_mask, cont_par, sc_tpl
