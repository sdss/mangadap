#!/usr/bin/env python3

import os
import time
import numpy
import glob

from IPython import embed

from astropy.io import fits
import astropy.constants

#import matplotlib
#matplotlib.use('Agg')
#from matplotlib import pyplot

from mangadap.dapfits import DAPQualityBitMask, DAPMapsBitMask
from mangadap.config import defaults
from mangadap.util.datatable import DataTable
from mangadap.util.fileio import channel_dictionary


class ReferencePropertyTable(DataTable):
    def __init__(self, shape=None):
        datamodel = dict(D4000BIN=dict(typ=float, shape=(2,)),
                         SIGMABIN=dict(typ=float, shape=(2,)),
                         HAEWBIN=dict(typ=float, shape=(2,)),
                         SNRBIN=dict(typ=float, shape=(2,)),
                         NBIN=dict(typ=int, shape=None),
                         PLT=dict(typ=int, shape=None),
                         IFU=dict(typ=int, shape=None),
                         BIN=dict(typ=int, shape=None),
                         MNGTARG1=dict(typ=int, shape=None),
                         MNGTARG3=dict(typ=int, shape=None),
                         Z=dict(typ=float, shape=None),
                         SNR=dict(typ=float, shape=None),
                         STZ=dict(typ=float, shape=None),
                         SIGMA=dict(typ=float, shape=None),
                         HAANR=dict(typ=float, shape=None),
                         HAZ=dict(typ=float, shape=None),
                         HASIGMA=dict(typ=float, shape=None),
                         HAEW=dict(typ=float, shape=None),
                         D4000=dict(typ=float, shape=None))

        keys = list(datamodel.keys())
        super(ReferencePropertyTable,
                self).__init__(keys, [datamodel[k]['typ'] for k in keys],
                               element_shapes=[datamodel[k]['shape'] for k in keys],
                               shape=shape)


def get_maps_files(dapver='3.0.1', daptype='SPX-MILESHC-MASTARHC2'):
#    return ['/uufs/chpc.utah.edu/common/home/sdss/mangawork/manga/spectro/analysis/v3_0_1/{0}/{1}/7443/12704/manga-7443-12704-MAPS-{1}.fits.gz'.format(dapver, daptype), '/uufs/chpc.utah.edu/common/home/sdss/mangawork/manga/spectro/analysis/v3_0_1/{0}/{1}/7443/12701/manga-7443-12701-MAPS-{1}.fits.gz'.format(dapver, daptype), '/uufs/chpc.utah.edu/common/home/sdss/mangawork/manga/spectro/analysis/v3_0_1/{0}/{1}/7443/12703/manga-7443-12703-MAPS-{1}.fits.gz'.format(dapver, daptype), '/uufs/chpc.utah.edu/common/home/sdss/mangawork/manga/spectro/analysis/v3_0_1/{0}/{1}/7443/12705/manga-7443-12705-MAPS-{1}.fits.gz'.format(dapver, daptype), '/uufs/chpc.utah.edu/common/home/sdss/mangawork/manga/spectro/analysis/v3_0_1/{0}/{1}/7443/1901/manga-7443-1901-MAPS-{1}.fits.gz'.format(dapver, daptype), '/uufs/chpc.utah.edu/common/home/sdss/mangawork/manga/spectro/analysis/v3_0_1/{0}/{1}/7443/1902/manga-7443-1902-MAPS-{1}.fits.gz'.format(dapver, daptype), '/uufs/chpc.utah.edu/common/home/sdss/mangawork/manga/spectro/analysis/v3_0_1/{0}/{1}/7443/12702/manga-7443-12702-MAPS-{1}.fits.gz'.format(dapver, daptype)]

    root = os.path.join(defaults.dap_analysis_path(dapver=dapver), daptype)
    return glob.glob(os.path.join(root, '*', '*', '*MAPS*fits.gz'))

# This won't work with hybrid binning scheme!
def include_maps_data(tbl, maps_file, exclude_ancillary=False):
    qual_bm = DAPQualityBitMask()
    bm = DAPMapsBitMask()
   
#    # Targeting bits
#    sdssMaskbits = defaults.sdss_maskbits_file()
#    mngtarg1_bm = BitMask.from_par_file(sdssMaskbits, 'MANGA_TARGET1')
#    mngtarg3_bm = BitMask.from_par_file(sdssMaskbits, 'MANGA_TARGET3')

    with fits.open(maps_file) as hdu:
        # Check if the file was CRITICAL
        if qual_bm.flagged(hdu['PRIMARY'].header['DAPQUAL'], flag='DRPCRIT'):
            return 1
        if exclude_ancillary and hdu[0].header['MNGTARG1'] == 0:
            return 1

        d4000_bins = numpy.unique(tbl['D4000BIN'])
        sigma_bins = numpy.unique(tbl['SIGMABIN'])
        haew_bins = numpy.unique(tbl['HAEWBIN'])
        snr_bins = numpy.unique(tbl['SNRBIN'])
        nsigma = len(sigma_bins)-1
        nd4000 = len(d4000_bins)-1
        nhaew = len(haew_bins)-1
        nsnr = len(snr_bins)-1

        el = channel_dictionary(hdu, 'EMLINE_GFLUX')
        si = channel_dictionary(hdu, 'SPECINDEX')

        plt, ifu = map(lambda x : int(x), os.path.split(maps_file)[1].split('-')[1:3])

        # Input redshift
        z = hdu['PRIMARY'].header['SCINPVEL'] / astropy.constants.c.to('km/s').value

        # Stellar-kinematics binning map
        binid, bin_indx = map(lambda x: x[1:], numpy.unique(hdu['BINID'].data[1,:,:].ravel(),
                                                            return_index=True))

        # Get the D4000, H-alpha EW, and sigma values with S/N greater
        # than one and unmasked
        good_snr = hdu['SPX_SNR'].data.ravel() > 1
        good_haew = numpy.logical_not(bm.flagged(hdu['EMLINE_GEW_MASK'].data[el['Ha-6564']].ravel(),
                                                 flag='DONOTUSE'))
        good_d4000 = numpy.logical_not(bm.flagged(hdu['SPECINDEX_MASK'].data[si['D4000']].ravel(),
                                                  flag='DONOTUSE'))
        good_sigma = numpy.logical_not(bm.flagged(hdu['STELLAR_SIGMA_MASK'].data.ravel(),
                                                  flag='DONOTUSE'))

        indx = (good_snr & good_haew & good_d4000 & good_sigma)[bin_indx]
        snr = hdu['SPX_SNR'].data.ravel()[bin_indx][indx]
        stvel = hdu['STELLAR_VEL'].data.ravel()[bin_indx][indx]
        sigma = hdu['STELLAR_SIGMA'].data.ravel()[bin_indx][indx]
        haanr = hdu['EMLINE_GANR'].data[el['Ha-6564']].ravel()[bin_indx][indx]
        havel = hdu['EMLINE_GVEL'].data[el['Ha-6564']].ravel()[bin_indx][indx]
        hasig = hdu['EMLINE_GSIGMA'].data[el['Ha-6564']].ravel()[bin_indx][indx]
        haew = hdu['EMLINE_GEW'].data[el['Ha-6564']].ravel()[bin_indx][indx]
        d4000 = hdu['SPECINDEX'].data[si['D4000']].ravel()[bin_indx][indx]
        binid = binid[indx]

        d4000_i = numpy.digitize(d4000, d4000_bins[1:-1])
        sigma_j = numpy.digitize(sigma, sigma_bins[1:-1])
        haew_k = numpy.digitize(haew, haew_bins[1:-1])
        snr_l = numpy.digitize(snr, snr_bins[1:-1])

        for i in range(nd4000):
            dindx = d4000_i == i
            for j in range(nsigma):
                sindx = sigma_j == j
                for k in range(nhaew):
                    hindx = haew_k == k
                    for l in range(nsnr):
                        nindx = snr_l == l

                        indx = dindx & sindx & hindx & nindx

                        nbin = numpy.sum(indx)
                        if nbin == 0:
                            continue
                       
                        ii = i*nsigma*nhaew*nsnr + j*nhaew*nsnr + k*nsnr + l
                        tbl['NBIN'][ii] += nbin

                        _snr = snr[indx]
                        # Still find the highest S/N spectrum in each
                        # bin
                        m = numpy.argsort(_snr)[-1]
                        if _snr[m] < tbl['SNR'][ii]:
                            continue

                        tbl['PLT'][ii] = plt
                        tbl['IFU'][ii] = ifu
                        tbl['BIN'][ii] = binid[indx][m]
                        tbl['MNGTARG1'][ii] = int(hdu[0].header['MNGTARG1'])
                        tbl['MNGTARG3'][ii] = int(hdu[0].header['MNGTARG3'])
                        tbl['Z'][ii] = z
                        tbl['SNR'][ii] = snr[indx][m]
                        tbl['STZ'][ii] \
                            = stvel[indx][m]*(1+z)/astropy.constants.c.to('km/s').value + z
                        tbl['SIGMA'][ii] = sigma[indx][m]
                        tbl['HAANR'][ii] = haanr[indx][m]
                        tbl['HAZ'][ii] \
                            = havel[indx][m]*(1+z)/astropy.constants.c.to('km/s').value + z
                        tbl['HASIGMA'][ii] = hasig[indx][m]
                        tbl['HAEW'][ii] = haew[indx][m]
                        tbl['D4000'][ii] = d4000[indx][m]

    return 0
    

def main():
    ofile = 'representative_spectra_selection_v2.fits'
    force = False
    if os.path.isfile(ofile) and not force:
        print('File already exists.')
        return

    exclude_ancillary = True
    daptype = 'SPX-MILESHC-MASTARHC2'
    maps_files = get_maps_files()
    nmaps = len(maps_files)
    print('Found {0} MAPS files.'.format(nmaps))

    # Initialize the output database
    # Sigmas are *uncorrected*
    d4000_bins = numpy.array([ 1.0, 1.2, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4 ])
    sigma_bins = numpy.array([0, 25, 50, 75, 100, 150, 200, 250, 2000])
    haew_bins = numpy.array([-1,2,8,16,32,128])
    snr_bins = numpy.array([1,2,5,10,20,40,500])
    nsigma = len(sigma_bins)-1
    nd4000 = len(d4000_bins)-1
    nhaew = len(haew_bins)-1
    nsnr = len(snr_bins)-1

    print(nd4000*nsigma*nhaew*nsnr)
    tbl = ReferencePropertyTable(shape=(nd4000*nsigma*nhaew*nsnr,))

    # Ensure the number in the bin is initialized to 0
    tbl['NBIN'] = 0

    # Set the bin limits
    tbl['D4000BIN'][:,0], tbl['SIGMABIN'][:,0], tbl['HAEWBIN'][:,0], tbl['SNRBIN'][:,0] \
            = tuple( map(lambda x: x.ravel(), numpy.meshgrid(d4000_bins[:-1], sigma_bins[:-1],
                                                             haew_bins[:-1], snr_bins[:-1],
                                                             indexing='ij')))
    tbl['D4000BIN'][:,1], tbl['SIGMABIN'][:,1], tbl['HAEWBIN'][:,1], tbl['SNRBIN'][:,1] \
            = tuple( map(lambda x: x.ravel(), numpy.meshgrid(d4000_bins[1:], sigma_bins[1:],
                                                             haew_bins[1:], snr_bins[1:],
                                                             indexing='ij')))

    # Iterate through the maps files and add the data
    nskipped = 0
    for i in range(nmaps):
        print('{0}/{1}'.format(i+1,nmaps), end='\r')
        nskipped += include_maps_data(tbl, maps_files[i], exclude_ancillary=exclude_ancillary)
    print('{0}/{0}'.format(nmaps))
    print('Number of skipped datacubes: {0}/{1}'.format(nskipped, nmaps))

    # Write the output
    print('Writing: {0}'.format(ofile))
    fits.HDUList([fits.PrimaryHDU(), tbl.to_hdu(name='PAR')]).writeto(ofile, overwrite=True)

#-----------------------------------------------------------------------------

if __name__ == '__main__':
    t = time.perf_counter()
    main()
    print('Elapsed time: {0} seconds'.format(time.perf_counter() - t))



