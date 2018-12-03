#!/usr/bin/env python3
# Script will spectrally decompose a MaNGA datacube
# Python environment is sdss_snakes2
#Working environment is /Users/ppzaf/Work/MaNGA/dap/master/examples

#Imports
import matplotlib, os, pickle
from astropy.io import fits
from astropy.table import Table
import numpy as np
import numpy.ma as ma
from numpy import inf
from numpy.polynomial import legendre
import glob
from os import path
import os
import subprocess
import matplotlib.pyplot as plt
from scipy import ndimage
from sys import stdout
from speclite import redshift
import fileinput
import requests
from requests.auth import HTTPBasicAuth
import time
# MaNGA functions
from mangadap.drpfits import DRPFits
from mangadap.proc.reductionassessments import *
from mangadap.proc.spatiallybinnedspectra import SpatiallyBinnedSpectra
from mangadap.proc.stellarcontinuummodel import StellarContinuumModel
# My functions:
import galfit_utils
from kinematic_utils import *
#import feedme_files
#import make_images

#-----------------------------------------------------------------------------
# Step 1 - TODO: Run Galfit on SDSS r-band image
# - Run Galfit on MaNGA r-band image √
# - Use Galfit result to create disk-only and bulge-only image at MaNGA pixel
# scale with MaNGA PSF √
#-----------------------------------------------------------------------------
# Step 2 - Binning and first fit of kinematics
#-----------------------------------------------------------------------------
# - Bin the DRP cube to a given S:N √
# - Bin the disk-only and bulge-only images from Galfit to the same S:N √
# - Fit stellar kinematics for each bin assuming a single kinematic component
#using SSPs √
#
# DRPFits, ReductionAssesments, SpatiallyBinnedSpectra, StellarContinuumModel
#-----------------------------------------------------------------------------
# Step 3 - Construct disk-only and bulge-only templates
#-----------------------------------------------------------------------------
# - Deredshift all spectra - Speclite √
# - Construct disk-only and bulge-only spectra by weighting the deredshifted
# cube by the binned GalFit images √
# - Fit the bulge and disk spectra with SSPs
# - Use weighted templates to construct "bulge" and "disk" template spectra
#
# SpectralStack (or SpecLite), TemplateLibrary, PPXFFit
#-----------------------------------------------------------------------------
# Step 4 - Fit kinematics of disk and bulge
#-----------------------------------------------------------------------------
# - Fit all bins with the bulge and disk templates set in their own component
# to get the kinematics for each bin
# - Details of fitting the stellar kinematics for two components should be
# worked out: 5x5 grid in velocities (bulge v vs. disk v), or base on the
# original, 1-component kinematics fit.
#
# pPXFFit, NewClass
#-----------------------------------------------------------------------------
# Step 5 - Optimize the bulge-only and disk-only template mix (for each bin)
#-----------------------------------------------------------------------------
# - Fit all bins with all SSP templates but with fixed kinematics
# - Construct a disk-only and bulge-only template spectrum, either one for
# each bin or one for all bins
#
# PPXFFit
#-----------------------------------------------------------------------------
# Step 6 - Iteration
#-----------------------------------------------------------------------------
# - Iterate over 4 and 5 until the kinematics and template mix "converge" or
# after a certain fixed number of iterations
#-----------------------------------------------------------------------------
# Step 7 - Finalize the stellar population properties
#-----------------------------------------------------------------------------
# - Use the weights from the final fit to determine the stellar population
# properties/star-formation history
#-----------------------------------------------------------------------------
# Function definitions
def readspectra(plate, id, vor):
    sdss_verify = open('/Users/ppzaf/Documents/2017/MaNGA/sdss_password.txt')
    up = sdss_verify.readlines()

    cube_top_level_url='https://data.sdss.org/sas/mangawork/manga/spectro/redux/MPL-7/'
    if os.path.isfile('/Users/ppzaf/Work/MaNGA/dap/master/examples/manga-'+str(plate)+'-'+str(id)+'-LOGCUBE.fits.gz'):
        print('Cube file exists!')
        pass
    else:
        r = requests.get(cube_top_level_url+str(plate)+'/stack/manga-'+str(plate)+'-'+str(id)+'-LOGCUBE.fits.gz', auth=HTTPBasicAuth(up[0].rstrip('\n'), up[1].rstrip('\n')), stream=True)
        if r.ok:
            print('##############################################')
            print('Obtaining CUBEFILE...')
            print('##############################################')
            with open('/Users/ppzaf/Work/MaNGA/dap/master/examples/manga-'+str(plate)+'-'+str(id)+'-LOGCUBE.fits.gz', 'wb') as file:
                file.write(r.content)
        else:
            print('No MaNGA data for plate and ID provided: ' +cube_top_level_url+str(plate)+'/stack/manga-'+str(plate)+'-'+str(id)+'-LOGCUBE.fits.gz')
        r = requests.get(cube_top_level_url+str(plate)+'/stack/manga-'+str(plate)+'-'+str(id)+'-LOGRSS.fits.gz', auth=HTTPBasicAuth(up[0].rstrip('\n'), up[1].rstrip('\n')), stream=True)
        if r.ok:
            print('##############################################')
            print('Obtaining RSSFILE...')
            print('##############################################')
            with open('/Users/ppzaf/Work/MaNGA/dap/master/examples/manga-'+str(plate)+'-'+str(id)+'-LOGRSS.fits.gz', 'wb') as file:
                file.write(r.content)
        else:
            print('No MaNGA data for plate and ID provided: ' +cube_top_level_url+str(plate)+'/stack/manga-'+str(plate)+'-'+str(id)+'-LOGRSS.fits.gz')

    hdulist = fits.open('/Users/ppzaf/Work/MaNGA/dap/master/examples/manga-'+str(plate)+'-'+str(id)+'-LOGCUBE.fits.gz') #Extract what we need from the cube file.
    flux = hdulist['FLUX'].data
    ivar = hdulist['IVAR'].data
    mask = hdulist['MASK']. data
    wave = hdulist['WAVE'].data
    rimg = hdulist['RIMG'].data
    rpsf = hdulist['RPSF'].data
    #Reorder Flux and mask arrays from (wavelength, DEC, RA) to (RA, DEC, wavelength)
    flux = np.transpose(flux, axes=(2,1,0))
    ivar = np.transpose(ivar, axes=(2,1,0))
    mask = np.transpose(mask, axes=(2,1,0))
    # MAPS
    maps_top_level_url='https://data.sdss.org/sas/mangawork/manga/spectro/analysis/MPL-7/'
    if os.path.isfile('/Users/ppzaf/Work/MaNGA/dap/master/examples/manga-'+str(plate)+'-'+str(id)+'-MAPS-'+str(vor)+'-GAU-MILESHC.fits.gz'):
        print('Maps file exists!')
        pass
    else:
        r = requests.get(maps_top_level_url+str(vor)+'-GAU-MILESHC/'+str(plate)+'/'+str(id)+'/manga-'+str(plate)+'-'+str(id)+'-MAPS-'+str(vor)+'-GAU-MILESHC.fits.gz', auth=HTTPBasicAuth(up[0].rstrip('\n'), up[1].rstrip('\n')), stream=True)
        if r.ok:
            print('##############################################')
            print('Obtaining MAPSFILE...')
            print('##############################################')
            with open('/Users/ppzaf/Work/MaNGA/dap/master/examples/manga-'+str(plate)+'-'+str(id)+'-MAPS-'+str(vor)+'-GAU-MILESHC.fits.gz', 'wb') as file:
                file.write(r.content)
        else:
            print('No MaNGA data for plate and ID provided: ' +maps_top_level_url+str(vor)+'-GAU-MILESHC/'+str(plate)+'/'+str(id)+'/manga-'+str(plate)+'-'+str(id)+'-MAPS-'+str(vor)+'-GAU-MILESHC.fits.gz')
    hdulist2 = fits.open('/Users/ppzaf/Work/MaNGA/dap/master/examples/manga-'+str(plate)+'-'+str(id)+'-MAPS-'+str(vor)+'-GAU-MILESHC.fits.gz')
    stellar_vel = hdulist2['STELLAR_VEL'].data
    stellar_sig = hdulist2['STELLAR_SIGMA'].data #Must be corrected using stellar_sigmascorr quantity to obtain astrophysical dispersion
    stellar_sig_corr = hdulist2['STELLAR_SIGMACORR'].data
    stellar_mask = hdulist2['STELLAR_VEL_MASK'].data
    header = hdulist2['PRIMARY'].header
    sigma_sig = np.sqrt( np.power(stellar_sig,2) - np.power(stellar_sig_corr,2) )
    #os.remove('/Users/ppzaf/Work/MaNGA/dap/master/examples/manga-'+str(plate)+'-'+str(ifu)+'-LOGCUBE-'+str(vor)+'-GAU-MILESHC.fits.gz')
    #os.remove('/Users/ppzaf/Work/MaNGA/dap/master/examples/manga-'+str(plate)+'-'+str(ifu)+'-MAPS-'+str(vor)+'-GAU-MILESHC.fits.gz')

    return flux, ivar, mask, wave, rimg, rpsf, stellar_vel, stellar_sig, header

def run_galfitm(plate, ifu, feedme_filename):
    #narrowband_feedme_filename = '%s-%s-narrowband-feedme' % (plate, ifu)
    subprocess.call("/usr/local/src/galfitM/galfitm-1.2.1-osx %s" % (feedme_filename), shell=True )
#
#
#-----------------------------------------------------------------------------
# STEP 1
#-----------------------------------------------------------------------------
#TODO: Galfit step on Sloan image. For now, skip straight to MaNGA.
c = 299792.458 # km/s
#NOTE: Change these!!!
Rb = 3.19 #in arcsec
Rd = 7.92
vor='HYB10'
start = time.time()
# 1) Import cube - Both maps and cubefile are flipped up-down c.f. Marvin. You HAVEN'T fixed this yet.

t = Table.read('drpall-v2_4_3.fits')
for i in range (10895, 10896): # 7443-9102 is for testing.
    print('##############################')
    print('Importing MaNGA cubes and maps')
    print('##############################')
    flux, ivar, mask, wave, rimg, rpsf, stellar_vel, stellar_sig, header = readspectra(t['plate'][i], t['ifudsgn'][i], vor)
    #print('###################################')
    #print('de-redshifting datacube')
    #print('###################################')
    #flux_m = np.ma.array(flux, mask = mask)
    #ivar_m = np.ma.array(ivar, mask=mask)
    #waves, fluxes, ivars, header = deredshift(flux_m, wave, stellar_vel, stellar_sig, ivar_m, c, t['nsa_z'][i])
    print('###################################')
    print('Preparing single broadband image')
    print('###################################')
    summed_image = galfit_utils.broadband_images(flux, ivar, mask, header, t['plate'][i], t['ifudsgn'][i], vor)
    #TODO: Find a way of incorporating Rb and Rd into drpall file? Match with Simard.
    galfit_utils.write_feedme(t['plate'][i], t['ifudsgn'][i], summed_image, Rb, Rd, header)
    print('########################################')
    print('Running Galfit on single broadband image')
    print('########################################')
    feedme_filename = '%s-%s-feedme' % (t['plate'][i], t['ifudsgn'][i])
    run_galfitm(t['plate'][i], t['ifudsgn'][i], feedme_filename)

#-----------------------------------------------------------------------------
# STEP 2
#-----------------------------------------------------------------------------
# Bin cube to desired SN (Let's start with 30 and hopefully reduce)
#drpf = DRPFits(obs['plate'], obs['ifudesign'], obs['mode'], drpver=_drpver, redux_path=redux_path, directory_path=directory_path, read=True)
    drpver=os.environ['MANGADRP_VER']
    redux_path = os.path.join(os.environ['MANGA_SPECTRO_REDUX'], drpver)
    directory_path = os.path.join(redux_path, str(t['plate'][i]), 'stack')
    drpf = DRPFits(t['plate'][i], t['ifudsgn'][i], 'CUBE', drpver=drpver, redux_path = redux_path, directory_path=directory_path, read=True)
    #Use defaults.py?
    pa = 51.07
    ell = 0.25 # (1-b/a)
    #rdxqa = ReductionAssessment('SNRR', drpf, pa=obs['pa'], ell=obs['ell'], clobber=False)
    rdxqa = ReductionAssessment('SNRR', drpf, pa=pa, ell=ell, clobber=False)
    reff = 4.956
    binned_spectra = SpatiallyBinnedSpectra('VOR30', drpf, rdxqa, reff=reff, clobber=False)
    # - Bin the disk-only and bulge-only images from Galfit to the same S:N
    analysis_path = os.path.join(os.environ['MANGA_SPECTRO_ANALYSIS'], drpver, os.environ['MANGADAP_VER'], 'common', str(t['plate'][i]), str(t['ifudsgn'][i]))
    hdu1 = fits.open(os.path.join(analysis_path, 'manga-'+ str(t['plate'][i]) + '-' + str(t['ifudsgn'][i]) + '-SNRR-VOR30.fits.gz'))
    binid = hdu1['BINID'].data

    #Load in disk and bulge-only images from subcomps.fits
    hdu2 = fits.open('subcomps.fits')
    bulge = hdu2[1].data #image
    disk = hdu2[2].data #image
    # Sum the flux in each Voronoi bin for bulge and disk regions
    # TODO: My god, find a more efficient way to do this!
    new_bulge_array, new_disk_array = np.empty_like(binid), np.empty_like(binid)
    summed_fluxes_bulge, summed_fluxes_disk = np.zeros(np.max(binid)+1), np.zeros(np.max(binid)+1)
    for k in range(0,np.max(binid)): #for each bin 1-109 in test case:
        for l in range(0, binid.shape[0]):
            for j in range(0, binid.shape[1]):
                if binid[l,j] == k:
                    summed_fluxes_bulge[k]+=(bulge[l,j])
                    summed_fluxes_disk[k]+=(disk[l,j])
    #Now assign the summed fluxes back to a map
    for l in range(0, binid.shape[0]):
        for j in range(0, binid.shape[1]):
            new_bulge_array[l,j] = summed_fluxes_bulge[binid[l,j]]
            new_disk_array[l,j] = summed_fluxes_disk[binid[l,j]]
    bulge_frac = np.divide(new_bulge_array, np.add(new_bulge_array, new_disk_array))
    z = 0.0916
    guess_vel = c*z
    guess_sig = 100
    stellar_continuum = StellarContinuumModel('GAU-MILESSSP', binned_spectra, guess_vel=guess_vel, guess_sig=guess_sig, clobber=False)

    #-----------------------------------------------------------------------------
    # STEP 3
    #-----------------------------------------------------------------------------
    #Grab the kinematics from the single component fit
    hdu4 = fits.open(os.path.join(os.environ['MANGA_SPECTRO_ANALYSIS'], drpver, 'master','common', str(t['plate'][i]) , str(t['ifudsgn'][i]) ,'manga-'+str(t['plate'][i])+'-'+str(t['ifudsgn'][i])+'-SNRR-VOR30.fits.gz')) #For the voronoi bins, binned flux and ivar, distance to bins, and bin_mflux
    hdu3 = fits.open(os.path.join(os.environ['MANGA_SPECTRO_ANALYSIS'], drpver, 'master', 'VOR30-GAU-MILESSSP', str(t['plate'][i]) , str(t['ifudsgn'][i]) , 'ref','manga-'+str(t['plate'][i])+'-'+str(t['ifudsgn'][i])+'-SNRR-VOR30-GAU-MILESSSP.fits.gz')) #For the kinematics only
    par = hdu3['PAR'].data #par.columns
    kin = par['KIN'] # These vals are total vel and sigma one for each bin
    stell_vel = kin[:,0]
    stell_sig = kin[:,1]
    flux_b = hdu4['FLUX'].data
    ivar_b = hdu4['IVAR'].data
    mask_b = hdu4['MASK'].data
    wave_b = hdu4['WAVE'].data
    bins = hdu4['BINS'].data
    lwellcoo = bins['LW_ELL_COO']
    bin_disk = lwellcoo[:,0] #in arcsec. x2 to get spaxel vals
    signal = bins['SIGNAL'] # I'm not 100% sure what this is, but I think it's related to the flux in each bin, so proxy for bin_mflux.
    flux_m = np.ma.array(flux_b, mask = mask_b)
    ivar_m = np.ma.array(ivar_b, mask=mask_b)
    print('###################################')
    print('de-redshifting datacube')
    print('###################################')
    #Deredshift spectra (do I need to broaden dispersion also?) NOTE: Broadening for now, can take out later if req.
    flux_m = np.transpose(flux_m, axes=(1,0))
    ivar_m = np.transpose(ivar_m, axes=(1,0))
    waves, fluxes, ivars = deredshift(flux_m, wave_b, ivar_m, stell_vel, stell_sig, c, z)

    # Construct disk-only and bulge-only spectra by weighting the deredshifted cube by the binned GalFit images
    # Alternatively, create a single bulge and disk spectrum the way I've done it before.
    bulge_spec, disk_spec = create_single_spec(fluxes, bin_disk, signal, Rb)
    # Fit the bulge and disk spectra with SSPs
    
