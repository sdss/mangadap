#!/usr/bin/env python3

import os
import time
from IPython import embed
import numpy
from astropy.io import fits
import astropy.constants
from matplotlib import pyplot

from mangadap.util.drpfits import DRPFitsBitMask, DRPFits
from mangadap.dapfits import DAPCubeBitMask
from mangadap.config import defaults

#-----------------------------------------------------------------------------

def collate_data(ofile, tbl, dapver=None, daptype=None):

    # The number of pixels is HARDWIRED!
    nspec = len(tbl)
    npix = 4563

    drpbm = DRPFitsBitMask()
    cube_bm = DAPCubeBitMask()
    flux = numpy.zeros((nspec, npix), dtype=float)
    ivar = numpy.zeros((nspec, npix), dtype=float)
    mask = numpy.zeros((nspec, npix), dtype=drpbm.minimum_dtype())
    model = numpy.zeros((nspec, npix), dtype=float)
    model_mask = numpy.zeros((nspec, npix), dtype=cube_bm.minimum_dtype())
    emline = numpy.zeros((nspec, npix), dtype=float)
    stellar = numpy.zeros((nspec, npix), dtype=float)
    stellar_mask = numpy.zeros((nspec, npix), dtype=cube_bm.minimum_dtype())
    lsf = numpy.zeros((nspec, npix), dtype=float)
    redshift = numpy.zeros(nspec, dtype=float)
    filled = numpy.zeros(nspec, dtype=bool)

    for s in range(nspec):
        print('{0}/{1}'.format(s+1, nspec), end='\r')
        if tbl['NBIN'][s] == 0:
            mask[s,:] = drpbm.turn_on(mask[s,:], 'DONOTUSE')
            filled[s] = True
            continue
        if filled[s]:
            continue

        # Open both the DRP datacube and the DAP model datacube
        drp_hdu = fits.open(os.path.join(*DRPFits.default_paths(tbl['PLT'][s], tbl['IFU'][s],
                                                                'CUBE')))
        dap_hdu = fits.open(os.path.join(
                            defaults.dap_method_path(daptype, plate=tbl['PLT'][s],
                                                     ifudesign=tbl['IFU'][s], dapver=dapver),
                            defaults.dap_file_name(tbl['PLT'][s], tbl['IFU'][s], daptype,
                                                   mode='LOGCUBE')))

        # Find any selected spectra in this datacube
        indx = numpy.where((tbl['PLT'] == tbl['PLT'][s]) & (tbl['IFU'] == tbl['IFU'][s]))[0]
        if len(indx) == 0:
            embed(header='could not find pltifu')
            exit()

        try:
            i, j = numpy.array([list(map(lambda x : x[0],
                                     numpy.where(dap_hdu['BINID'].data[1] == tbl['BIN'][ii]))) 
                            for ii in indx]).T
        except:
            embed(header='bad bin selection')
            exit()

        wave = drp_hdu['WAVE'].data
        flux[indx,:] = drp_hdu['FLUX'].data[:,i,j].T
        ivar[indx,:] = drp_hdu['IVAR'].data[:,i,j].T
        mask[indx,:] = drp_hdu['MASK'].data[:,i,j].T
        model[indx,:] = dap_hdu['MODEL'].data[:,i,j].T
        model_mask[indx,:] = dap_hdu['MODEL_MASK'].data[:,i,j].T
        emline[indx,:] = dap_hdu['EMLINE'].data[:,i,j].T
        stellar[indx,:] = dap_hdu['STELLAR'].data[:,i,j].T
        stellar_mask[indx,:] = dap_hdu['STELLAR_MASK'].data[:,i,j].T
        lsf[indx,:] = dap_hdu['LSF'].data[:,i,j].T
        redshift[indx] = tbl['Z'][indx]
        filled[indx] = True

    print('{0}/{0}'.format(nspec))
    print('All filled: {0}'.format(str(numpy.all(filled))))
    
    fits.HDUList([fits.PrimaryHDU(),
                  fits.ImageHDU(data=wave, name='WAVE'),
                  fits.ImageHDU(data=flux, name='FLUX'),
                  fits.ImageHDU(data=ivar, name='IVAR'),
                  fits.ImageHDU(data=mask, name='MASK'),
                  fits.ImageHDU(data=model, name='MODEL'),
                  fits.ImageHDU(data=model_mask, name='MODEL_MASK'),
                  fits.ImageHDU(data=emline, name='EMLINE'),
                  fits.ImageHDU(data=stellar, name='STELLAR'),
                  fits.ImageHDU(data=stellar_mask, name='STELLAR_MASK'),
                  fits.ImageHDU(data=lsf, name='LSF'),
                  fits.ImageHDU(data=redshift, name='Z')]).writeto(ofile, overwrite=True)


#-----------------------------------------------------------------------------

def main():
    ifile = 'representative_spectra_selection_v2.fits'
    ofile = 'benchmark_spectra_v2.fits'
    dapver = '3.0.1'
    daptype = 'SPX-MILESHC-MASTARHC2'

    # Construct a single file with all the "representative" spectra
    with fits.open(ifile) as hdu:
        collate_data(ofile, hdu['PAR'].data, dapver=dapver, daptype=daptype)
    
if __name__ == '__main__':
    t = time.perf_counter()
    main()
    print('Elapsed time: {0} seconds'.format(time.perf_counter() - t))



