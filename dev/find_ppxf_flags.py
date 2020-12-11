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


class PPXFFlagTable(DataTable):
    def __init__(self, shape=None):
        datamodel = dict(PLT=dict(typ=int, shape=None),
                         IFU=dict(typ=int, shape=None),
                         CRIT=dict(typ=bool, shape=None),
                         MAIN=dict(typ=bool, shape=None),
                         HAMEAN=dict(typ=float, shape=None),
                         NSPX=dict(typ=int, shape=None),
                         NGOOD=dict(typ=int, shape=None),
                         NFLG=dict(typ=int, shape=None),
                         FLGSNR=dict(typ=float, shape=None),
                         EXSNR=dict(typ=float, shape=None),
                         EXBIN=dict(typ=int, shape=None))
        keys = list(datamodel.keys())
        super(PPXFFlagTable, self).__init__(keys, [datamodel[k]['typ'] for k in keys],
                                            element_shapes=[datamodel[k]['shape'] for k in keys],
                                            shape=shape)


def get_maps_files(dapver='3.0.1', daptype='SPX-MILESHC-MASTARHC2'):
#    return ['/uufs/chpc.utah.edu/common/home/sdss/mangawork/manga/spectro/analysis/v3_0_1/{0}/{1}/7443/12704/manga-7443-12704-MAPS-{1}.fits.gz'.format(dapver, daptype), '/uufs/chpc.utah.edu/common/home/sdss/mangawork/manga/spectro/analysis/v3_0_1/{0}/{1}/7443/12701/manga-7443-12701-MAPS-{1}.fits.gz'.format(dapver, daptype), '/uufs/chpc.utah.edu/common/home/sdss/mangawork/manga/spectro/analysis/v3_0_1/{0}/{1}/7443/12703/manga-7443-12703-MAPS-{1}.fits.gz'.format(dapver, daptype), '/uufs/chpc.utah.edu/common/home/sdss/mangawork/manga/spectro/analysis/v3_0_1/{0}/{1}/7443/12705/manga-7443-12705-MAPS-{1}.fits.gz'.format(dapver, daptype), '/uufs/chpc.utah.edu/common/home/sdss/mangawork/manga/spectro/analysis/v3_0_1/{0}/{1}/7443/1901/manga-7443-1901-MAPS-{1}.fits.gz'.format(dapver, daptype), '/uufs/chpc.utah.edu/common/home/sdss/mangawork/manga/spectro/analysis/v3_0_1/{0}/{1}/7443/1902/manga-7443-1902-MAPS-{1}.fits.gz'.format(dapver, daptype), '/uufs/chpc.utah.edu/common/home/sdss/mangawork/manga/spectro/analysis/v3_0_1/{0}/{1}/7443/12702/manga-7443-12702-MAPS-{1}.fits.gz'.format(dapver, daptype)]

    root = os.path.join(defaults.dap_analysis_path(dapver=dapver), daptype)
    return glob.glob(os.path.join(root, '*', '*', '*MAPS*fits.gz'))

# This won't work with hybrid binning scheme!
def ppxf_flag_info(tbl, maps_file):
    qual_bm = DAPQualityBitMask()
    bm = DAPMapsBitMask()
   
#    # Targeting bits
#    sdssMaskbits = defaults.sdss_maskbits_file()
#    mngtarg1_bm = BitMask.from_par_file(sdssMaskbits, 'MANGA_TARGET1')
#    mngtarg3_bm = BitMask.from_par_file(sdssMaskbits, 'MANGA_TARGET3')

    tbl['PLT'], tbl['IFU'] = map(lambda x : int(x), os.path.split(maps_file)[1].split('-')[1:3])

    with fits.open(maps_file) as hdu:

        # Check if the file was CRITICAL
        tbl['CRIT'] = qual_bm.flagged(hdu['PRIMARY'].header['DAPQUAL'], flag='DRPCRIT')
        tbl['MAIN'] = hdu[0].header['MNGTARG1'] > 0

        el = channel_dictionary(hdu, 'EMLINE_GFLUX')

        haflx = numpy.ma.MaskedArray(hdu['EMLINE_GFLUX'].data[el['Ha-6564']],
                                     mask=hdu['EMLINE_GFLUX_MASK'].data[el['Ha-6564']] > 0)
        tbl['HAMEAN'] = numpy.ma.median(haflx)
        tbl['NSPX'] = haflx.size
        tbl['NGOOD'] = numpy.sum(numpy.logical_not(numpy.ma.getmaskarray(haflx)))
        indx = bm.flagged(hdu['EMLINE_GFLUX_MASK'].data[el['Ha-6564']], 'FITFAILED')
        tbl['NFLG'] = numpy.sum(indx)

        if not numpy.any(indx):
            return

        snr = numpy.ma.MaskedArray(hdu['SPX_SNR'].data,
                                   mask=numpy.logical_not(hdu['SPX_SNR'].data > 0))
        tbl['FLGSNR'] = numpy.ma.median(snr[indx])
        srt = numpy.argsort(snr[indx])
        tbl['EXSNR'] = snr[indx][srt][-1]
        tbl['EXBIN'] = hdu['BINID'].data[3][indx][srt][-1]


def main():
    ofile = 'ppxf_failures_hyb.fits'
    force = False
    if os.path.isfile(ofile) and not force:
        print('File already exists.')
        return

    daptype = 'HYB10-MILESHC-MASTARHC2'
    maps_files = get_maps_files()
    nmaps = len(maps_files)
    print('Found {0} MAPS files.'.format(nmaps))

    tbl = PPXFFlagTable(shape=nmaps)

    # Iterate through the maps files and add the data
    for i in range(nmaps):
        print('{0}/{1}'.format(i+1,nmaps), end='\r')
        ppxf_flag_info(tbl[i], maps_files[i])
    print('{0}/{0}'.format(nmaps))

    # Write the output
    print('Writing: {0}'.format(ofile))
    fits.HDUList([fits.PrimaryHDU(), tbl.to_hdu(name='PAR')]).writeto(ofile, overwrite=True)

#-----------------------------------------------------------------------------

if __name__ == '__main__':
    t = time.perf_counter()
    main()
    print('Elapsed time: {0} seconds'.format(time.perf_counter() - t))



