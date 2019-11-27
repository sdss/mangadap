import sys
if sys.version > '3':
    long = int

import os
import gzip
import shutil
import logging
import numpy as np
import pdb

from astropy.io import fits
from astropy.io import ascii

from mangadap.util import sampling
from mangadap.util import fileio

## Need FWHM of resolution element in Ang
## Vacuum or air?


select_ages = np.arange(-1.2,1.2,0.1)
#allZ = np.array([0.001, 0.002, 0.003, 0.004, 0.006, 0.008, 0.010, 0.014, 0.020, 0.030, 0.040])
select_Z = np.array([0.008, 0.02, 0.03])

#bpass_fil = './v2.2.1_imf135_300/spectra-bin-imf135_300.z020.dat'
bpass_root = '/Users/rubin/Research/MaNGA/NaImcmc/Analysis/BPASS/v2.2.1_imf135_300/spectra-bin-imf135_300.z'
out_root = '/Users/rubin/Repo/manga/mangadap/data/spectral_templates/bpass/'

## Set up info common to all files
ncol = 52
col_arr = np.arange(1,ncol+1)
age_arr = 10.0**(6.0 + 0.1*(col_arr-2.0))
logage_arr = np.log10(age_arr / 1.0e9)

#pdb.set_trace()
for Z in select_Z:

    Z_str = "{:.3f}".format(Z)[2:]
    bpass_fil = bpass_root+Z_str+'.dat'

    data = ascii.read(bpass_fil)
    wave = data['col1']
    wvrange = (wave > 3600.0) & (wave < 10400.0)
    dlam = wave[wvrange] - np.roll(wave[wvrange],1)
    npix = len(dlam)

    crval = min(wave[wvrange])
    crpix = 1.0
    cdelt = np.median(dlam)
    #test_wv = (np.arange(1.0,npix+1) - crpix)*cdelt + crval

    for age in select_ages:

        ind = ((logage_arr>age-0.001) & (logage_arr<age+0.001))
        #print(logage_arr[ind])
        if(len(logage_arr[ind])>0):
            colname = 'col'+str(col_arr[ind][0])
            flux = data[colname]
            out_flux = flux[wvrange]

            age_str = "{:.2f}".format(age)
            out_fil = out_root + 'bin-imf135_300-z'+Z_str+'-T'+age_str+'.fits'
            print(out_fil)
   
            ## Maybe can just use fileio.writefits_1dspec
            fileio.writefits_1dspec(out_fil, crval, cdelt, out_flux, clobber=True)

            #step = sampling.spectral_coordinate_step(wave)

            #pdb.set_trace()
