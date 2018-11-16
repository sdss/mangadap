#!/usr/bin/env python3

# Script will spectrally decompose a MaNGA datacube

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
# My functions:
import galfit_utils
#import feedme_files
#import make_images

#-----------------------------------------------------------------------------
# Step 1 - TODO: Run Galfit on SDSS r-band image
# - Run Galfit on MaNGA r-band image
# - Use Galfit result to create disk-only and bulge-only image at MaNGA pixel
# scale with MaNGA PSF NOTE: Ask Steven how to convolve the model to the orig
# pixel scale
#-----------------------------------------------------------------------------
# Step 2 - Binning and first fit of kinematics
#-----------------------------------------------------------------------------
# - Bin the DRP cube to a given S:N
# - Bin the disk-only and bulge-only images from Galfit to the same S:N
# - Fit stellar kinematics for each bin assuming a single kinematic component
#using SSPs
#
# DRPFits, ReductionAssesments, SpatiallyBinnedSpectra, StellarContinuumModel
#-----------------------------------------------------------------------------
# Step 3 - Construct disk-only and bulge-only templates
#-----------------------------------------------------------------------------
# - Deredshift all spectra - Speclite
# - Construct disk-only and bulge-only spectra by weighting the deredshifted
# cube by the binned GalFit images
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
    if os.path.isfile('/Users/ppzaf/Work/MaNGA/dap/master/temp_files/manga-'+str(plate)+'-'+str(id)+'-LOGCUBE.fits.gz'):
        print('Cube file exists!')
        pass
    else:
        r = requests.get(cube_top_level_url+str(plate)+'/stack/manga-'+str(plate)+'-'+str(id)+'-LOGCUBE.fits.gz', auth=HTTPBasicAuth(up[0].rstrip('\n'), up[1].rstrip('\n')), stream=True)
        if r.ok:
            print('##############################################')
            print('Obtaining CUBEFILE...')
            print('##############################################')
            with open('/Users/ppzaf/Work/MaNGA/dap/master/temp_files/manga-'+str(plate)+'-'+str(id)+'-LOGCUBE.fits.gz', 'wb') as file:
                file.write(r.content)
        else:
            print('No MaNGA data for plate and ID provided: ' +cube_top_level_url+str(plate)+'/stack/manga-'+str(plate)+'-'+str(id)+'-LOGCUBE.fits.gz')
        r = requests.get(cube_top_level_url+str(plate)+'/stack/manga-'+str(plate)+'-'+str(id)+'-LOGRSS.fits.gz', auth=HTTPBasicAuth(up[0].rstrip('\n'), up[1].rstrip('\n')), stream=True)
        if r.ok:
            print('##############################################')
            print('Obtaining RSSFILE...')
            print('##############################################')
            with open('/Users/ppzaf/Work/MaNGA/dap/master/temp_files/manga-'+str(plate)+'-'+str(id)+'-LOGRSS.fits.gz', 'wb') as file:
                file.write(r.content)
        else:
            print('No MaNGA data for plate and ID provided: ' +cube_top_level_url+str(plate)+'/stack/manga-'+str(plate)+'-'+str(id)+'-LOGRSS.fits.gz')

    hdulist = fits.open('/Users/ppzaf/Work/MaNGA/dap/master/temp_files/manga-'+str(plate)+'-'+str(id)+'-LOGCUBE.fits.gz') #Extract what we need from the cube file.
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
    if os.path.isfile('/Users/ppzaf/Work/MaNGA/dap/master/temp_files/manga-'+str(plate)+'-'+str(id)+'-MAPS-'+str(vor)+'-GAU-MILESHC.fits.gz'):
        print('Maps file exists!')
        pass
    else:
        r = requests.get(maps_top_level_url+str(vor)+'-GAU-MILESHC/'+str(plate)+'/'+str(id)+'/manga-'+str(plate)+'-'+str(id)+'-MAPS-'+str(vor)+'-GAU-MILESHC.fits.gz', auth=HTTPBasicAuth(up[0].rstrip('\n'), up[1].rstrip('\n')), stream=True)
        if r.ok:
            print('##############################################')
            print('Obtaining MAPSFILE...')
            print('##############################################')
            with open('/Users/ppzaf/Work/MaNGA/dap/master/temp_files/manga-'+str(plate)+'-'+str(id)+'-MAPS-'+str(vor)+'-GAU-MILESHC.fits.gz', 'wb') as file:
                file.write(r.content)
        else:
            print('No MaNGA data for plate and ID provided: ' +maps_top_level_url+str(vor)+'-GAU-MILESHC/'+str(plate)+'/'+str(id)+'/manga-'+str(plate)+'-'+str(id)+'-MAPS-'+str(vor)+'-GAU-MILESHC.fits.gz')
    hdulist2 = fits.open('/Users/ppzaf/Work/MaNGA/dap/master/temp_files/manga-'+str(plate)+'-'+str(id)+'-MAPS-'+str(vor)+'-GAU-MILESHC.fits.gz')
    stellar_vel = hdulist2['STELLAR_VEL'].data
    stellar_sig = hdulist2['STELLAR_SIGMA'].data #Must be corrected using stellar_sigmascorr quantity to obtain astrophysical dispersion
    stellar_sig_corr = hdulist2['STELLAR_SIGMACORR'].data
    stellar_mask = hdulist2['STELLAR_VEL_MASK'].data
    header = hdulist2['PRIMARY'].header
    sigma_sig = np.sqrt( np.power(stellar_sig,2) - np.power(stellar_sig_corr,2) )
    #os.remove('/Users/ppzaf/Work/MaNGA/dap/master/temp_files/manga-'+str(plate)+'-'+str(ifu)+'-LOGCUBE-'+str(vor)+'-GAU-MILESHC.fits.gz')
    #os.remove('/Users/ppzaf/Work/MaNGA/dap/master/temp_files/manga-'+str(plate)+'-'+str(ifu)+'-MAPS-'+str(vor)+'-GAU-MILESHC.fits.gz')

    return flux, ivar, mask, wave, rimg, rpsf, stellar_vel, stellar_sig, header

def deredshift(flux, wave, stellar_vel, stellar_sig, ivar, c, z):
    #Deredshift AND Dispersion-broaden cube. Output should be cube of deredshifted fluxes and wavelenghts
    #TODO: Does the ivar need to be broadened, too?
    velscale = np.log(wave[1]/wave[0])*c #Sould be 69 km/s
    sig_max = stellar_sig[int(flux.shape[0]/2), int(flux.shape[1]/2)]
    #First, broaden
    FWHM_diff, sigma_corr = [],[]
    flux_smooth = np.empty((flux.shape[0], flux.shape[1], flux.shape[2]))
    for i in range(flux.shape[0]):
        for j in range(flux.shape[1]):
            if stellar_sig[i,j] <= sig_max:
                FWHM_diff = np.sqrt((((2.355*sig_max)/velscale)**2)-((2.355*stellar_sig[i,j])/velscale)**2)
                sigma_corr = FWHM_diff/2.355 #?
                flux_smooth[i,j,:] = ndimage.filters.gaussian_filter1d(flux[i,j,:], sigma_corr)
            else:
                FWHM_diff = 0
                sigma_corr =0
                flux_smooth[i,j,:] = flux[i,j,:]


    z_map = (stellar_vel/c) +z #Mask should come into this somewhere.

    waves, fluxes, ivars = np.empty((flux.shape[0], flux.shape[1], flux.shape[2])),np.empty((flux.shape[0], flux.shape[1], flux.shape[2])), np.empty((flux.shape[0], flux.shape[1], flux.shape[2]))
    for i in range(flux.shape[0]):
        for j in range(flux.shape[1]):
            rules = [dict(name='wave', exponent=+1, array_in = wave), dict(name='flux', exponent=-1, array_in=flux_smooth[i,j]), dict(name='ivar', exponent=2, array_in=ivar[i,j])]#, dict(name='ivar', exponent=+2)]
            result = redshift(z_in = z_map[i,j], z_out=0, rules=rules)
            waves[i,j,:] = result['wave']
            fluxes[i,j,:] = result['flux']
            ivars[i,j,:] = result['ivar']
    return waves, fluxes, ivars, header

#
#
#-----------------------------------------------------------------------------
# STEP 1
#-----------------------------------------------------------------------------
#TODO: Galfit step on Sloan image. For now, skip straight to MaNGA.
c = 299792.458 # km/s
vor='HYB10'
start = time.time()
# 1) Import cube - Both maps and cubefile are flipped up-down c.f. Marvin. You HAVEN'T fixed this yet.

t = Table.read('drpall-v2_4_3.fits')
for i in range (10895, 10896): # 7443-9102 is for testing.
    print('##############################')
    print('Importing MaNGA cubes and maps')
    print('##############################')
    flux, ivar, mask, wave, rimg, rpsf, stellar_vel, stellar_sig, header = readspectra(t['plate'][i], t['ifudsgn'][i], vor)
    print('###################################')
    print('de-redshifting datacube')
    print('###################################')
    flux_m = np.ma.array(flux, mask = mask)
    ivar_m = np.ma.array(ivar, mask=mask)
    waves, fluxes, ivars, header = deredshift(flux_m, wave, stellar_vel, stellar_sig, ivar_m, c, t['nsa_z'][i])
    print('###################################')
    print('Preparing single broadband image')
    print('###################################')
    galfit_utils.broadband_images(flux, ivar, mask, header, t['plate'][i], t['ifudsgn'][i], vor)
