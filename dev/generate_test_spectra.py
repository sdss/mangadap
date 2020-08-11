#!/usr/bin/env python3

import os
import time

from argparse import ArgumentParser

import numpy
from astropy.io import fits
import astropy.constants
from matplotlib import pyplot

from mangadap.drpfits import DRPFitsBitMask, DRPFits
from mangadap.dapfits import DAPCubeBitMask, DAPFits
from mangadap.config.defaults import default_dap_file_name

#-----------------------------------------------------------------------------

def collate_data(tbl, method, minsnr=None, drpver=None, dapver=None, analysis_path=None,
                 ofile=None):

    # Channels with the H-alpha data and the D4000 data
    hac = 18
    d4c = 43

#    ntpl = tbl['TPLWGT'].shape[1]
#    print(ntpl)

    nspec = len(tbl)
    print(nspec)

    # Get the first entry with a valid plate number (rest should also be
    # valid!)
    indx = numpy.where(tbl['PLT'] > 0)[0][0]

    print(tbl['PLT'][indx], tbl['IFU'][indx])
    drpf = DRPFits(tbl['PLT'][indx], tbl['IFU'][indx], 'CUBE', read=True, drpver=drpver)
    print(drpf.file_path())

    npix = drpf['WAVE'].data.size
    print(npix)

    drpbm = DRPFitsBitMask()

    wave = drpf['WAVE'].data
    flux = numpy.zeros((nspec, npix), dtype=numpy.float)
    ivar = numpy.zeros((nspec, npix), dtype=numpy.float)
    mask = numpy.zeros((nspec, npix), dtype=drpbm.minimum_dtype())
    cont = numpy.zeros((nspec, npix), dtype=numpy.float)
    specres = numpy.zeros((nspec, npix), dtype=numpy.float)
    redshift = numpy.zeros(nspec, dtype=numpy.float)
#    tpluse = numpy.zeros( (nspec, ntpl), dtype=numpy.int)

    del drpf

    for s in range(nspec):
#    for s in range(1):
        print('{0}/{1}'.format(s+1, nspec))#, end='\r')
        if tbl['PLT'][s] == 0:
            print('Skipping!')
            mask[s,:] = drpbm.turn_on(mask[s,:], 'DONOTUSE')
            continue
        if minsnr is not None and tbl['SNR'][s] < minsnr:
            print('Skipping!')
            mask[s,:] = drpbm.turn_on(mask[s,:], 'DONOTUSE')
            continue

        # Get the maps file:
        dapf = DAPFits(tbl['PLT'][s], tbl['IFU'][s], method, drpver=drpver, dapver=dapver,
                       analysis_path=analysis_path)
        print(dapf.file_path())

        model_file = os.path.join(dapf.directory_path,
                        default_dap_file_name(tbl['PLT'][s], tbl['IFU'][s], method, mode='LOGCUBE'))
        mod_hdu = fits.open(model_file)
                       
        redshift[s] = dapf.guess_cz()/astropy.constants.c.to('km/s').value
#        print('redshift: ', redshift[s])
#        tpluse[s,:] = (tbl['TPLWGT'][s,:] > 0).astype(int)
#        pyplot.imshow(dapf['BINID'].data[:,:,0].T, origin='lower', interpolation='nearest')
#        pyplot.show()

        indx = dapf['BINID'].data[:,:,0] == tbl['BIN'][s]
        if numpy.sum(indx) == 0:
            raise ValueError('No bin {0} in DAP file.'.format(tbl['BIN'][s]))

        k = numpy.where(indx.ravel())[0][0]
        n = dapf.spatial_shape[0]
        i = k//n
        j = k - i*n

#        print(dapf['SPX_SNR'].data[i,j], tbl['SNR'][s], 0)
#        print(dapf['SPECINDEX'].data[i,j,d4c],tbl['D4000'][s],dapf['SPECINDEX_MASK'].data[i,j,d4c])
#        print(dapf['STELLAR_SIGMA'].data[i,j],tbl['SIGMA'][s],dapf['STELLAR_SIGMA_MASK'].data[i,j])
#        print(dapf['EMLINE_GEW'].data[i,j,hac],tbl['HAEW'][s],dapf['EMLINE_GEW_MASK'].data[i,j,hac])
#        print(dapf['BINID'].data[i,j,:])

        del dapf

 #       pyplot.plot(mod_hdu['WAVE'].data, mod_hdu['FLUX'].data[:,j,i])
 #       pyplot.show()
 #       exit()
        drpf = DRPFits(tbl['PLT'][s], tbl['IFU'][s], 'CUBE', read=True, drpver=drpver)
        flux[s,:] = drpf['FLUX'].data[i,j,:]
        ivar[s,:] = drpf['IVAR'].data[i,j,:]
        mask[s,:] = drpf['MASK'].data[i,j,:]
        specres[s,:] = drpf['SPECRES'].data.copy()
        del drpf
        cont[s,:] = mod_hdu['MODEL'].data[:,j,i] - mod_hdu['EMLINE'].data[:,j,i] \
                            - mod_hdu['EMLINE_BASE'].data[:,j,i]
        mod_hdu.close()
        del mod_hdu

#    pyplot.plot(wave, flux[0,:])
#    pyplot.show()
#    exit()

    print('DONE                      ')
    
    if ofile is not None:
        fits.HDUList([ fits.PrimaryHDU(),
                       fits.ImageHDU(data=wave, name='WAVE'),
                       fits.ImageHDU(data=flux, name='FLUX'),
                       fits.ImageHDU(data=ivar, name='IVAR'),
                       fits.ImageHDU(data=mask, name='MASK'),
                       fits.ImageHDU(data=cont, name='CONT'),
                       fits.ImageHDU(data=specres, name='SPECRES'),
                       fits.ImageHDU(data=redshift, name='Z'),
#                       fits.ImageHDU(data=tpluse, name='USETPL')
                     ]).writeto(ofile, overwrite=True)


#-----------------------------------------------------------------------------

if __name__ == '__main__':
    t = time.clock()

    ifile = 'manga_d4000_sig_haew_snr_selection.fits'
    ofile = 'benchmark_spectra_d4000_sig_haew_snr.fits.gz'

    # Read the table data
    hdu = fits.open(ifile)
    tbl = hdu['PAR'].data.copy()
    hdu.close()

    # Construct a single file with all the spectra with S/N > 20
    method = 'SPX-GAU-MILESHC'
    minsnr=None
    collate_data(tbl, method, minsnr=minsnr, analysis_path=None, ofile=ofile)
    
    print('Elapsed time: {0} seconds'.format(time.clock() - t))



