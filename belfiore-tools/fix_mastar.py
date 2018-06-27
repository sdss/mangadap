#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 27 09:33:32 2018

correct the MaSTAR spectra headers so that they can be correctly interpreted by the DAp

@author: francesco
"""

import glob
from astropy.io import fits

files=glob.glob('/Volumes/fbdata/CODE/mangadap/data/spectral_templates/mastar/*.fits')

for file0 in files:
    hdu=fits.open(file0, mode='update')
    aaa=hdu[0].header['CRDELT1']
    del hdu[0].header['CRDELT1']
    hdu[0].header['CDELT1'] =aaa
    hdu.flush()
    hdu.close()