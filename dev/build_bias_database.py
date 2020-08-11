
import os
import glob

from IPython import embed

import numpy
from matplotlib import pyplot

from astropy.io import fits

from mangadap.dapfits import DAPMapsBitMask, DAPQualityBitMask
from mangadap.config import defaults
from mangadap.util.datatable import DataTable
from mangadap.util.bitmask import BitMask
from mangadap.util.fileio import channel_dictionary
from mangadap.util.fitsutil import DAPFitsUtil

class HalphaBitMask(BitMask):
    def __init__(self):
        super(HalphaBitMask, self).__init__(['CRITCUBE', 'BADFLUX', 'BADEW', 'BADVEL', 'BADSIG'])


class HalphaDataTable(DataTable):
    def __init__(self, shape=None):
        bm = HalphaBitMask()
        _float = numpy.float32
        _int = numpy.int32
        datamodel = dict(PLATE=dict(typ=_int), IFU=dict(typ=_int), X=dict(typ=_int),
                         Y=dict(typ=_int), BINID=dict(typ=_int), R=dict(typ=_float),
                         R_RE=dict(typ=_float), R_KPC=dict(typ=_float), GCNT=dict(typ=_float),
                         GSNR=dict(typ=_float), FLUX=dict(typ=_float), FLUX_ERR=dict(typ=_float),
                         AMP=dict(typ=_float), ANR=dict(typ=_float), RCHI2=dict(typ=_float),
                         CNT=dict(typ=_float), EW=dict(typ=_float), EW_ERR=dict(typ=_float),
                         VEL=dict(typ=_float), VEL_ERR=dict(typ=_float), SIG=dict(typ=_float),
                         SIG_ERR=dict(typ=_float), SIG_INST=dict(typ=_float),
                         MASK=dict(typ=bm.minimum_dtype()))
        keys = list(datamodel.keys())
        super(HalphaDataTable, self).__init__(keys, [datamodel[k]['typ'] for k in keys],
                                              shape=shape)


def get_halpha_data(plt, ifu, hdu):
    el = channel_dictionary(hdu, 'EMLINE_GFLUX')
    uniq_bins, uniq_indx, uniq_cnt = map(lambda x : x[1:],
                                         numpy.unique(hdu['BINID'].data[3].ravel(),
                                                      return_index=True, return_counts=True))
    keep = numpy.logical_not(hdu['EMLINE_GFLUX_MASK'].data.ravel()[uniq_indx] > 0) \
            | numpy.logical_not(hdu['EMLINE_GEW_MASK'].data.ravel()[uniq_indx] > 0) \
            | numpy.logical_not(hdu['EMLINE_GVEL_MASK'].data.ravel()[uniq_indx] > 0) \
            | numpy.logical_not(hdu['EMLINE_GSIGMA_MASK'].data.ravel()[uniq_indx] > 0)
    uniq_bins = uniq_bins[keep]
    uniq_indx = uniq_indx[keep]
    uniq_cnt = uniq_cnt[keep]
    data = HalphaDataTable(shape=uniq_bins.shape)

    data['PLATE'] = plt
    data['IFU'] = ifu
    data['Y'], data['X'] = numpy.unravel_index(uniq_indx, hdu['BINID'].shape[1:])
    data['BINID'] = uniq_bins
    data['R'] = hdu['SPX_ELLCOO'].data[0].ravel()[uniq_indx]
    data['R_RE'] = hdu['SPX_ELLCOO'].data[1].ravel()[uniq_indx]
    data['R_KPC'] = hdu['SPX_ELLCOO'].data[2].ravel()[uniq_indx]
    data['GCNT'] = hdu['SPX_MFLUX'].data.ravel()[uniq_indx]
    data['GSNR'] = hdu['SPX_SNR'].data.ravel()[uniq_indx]
    data['FLUX'] = hdu['EMLINE_GFLUX'].data[el['Ha-6564']].ravel()[uniq_indx]
    data['FLUX_ERR'] = numpy.ma.power(
                            hdu['EMLINE_GFLUX_IVAR'].data[el['Ha-6564']].ravel()[uniq_indx],
                            -0.5).filled(-1)
    data['AMP'] = hdu['EMLINE_GA'].data[el['Ha-6564']].ravel()[uniq_indx]
    data['ANR'] = hdu['EMLINE_GANR'].data[el['Ha-6564']].ravel()[uniq_indx]
    data['RCHI2'] = hdu['EMLINE_LFOM'].data[el['Ha-6564']].ravel()[uniq_indx]
    data['CNT'] = hdu['EMLINE_GEW_CNT'].data[el['Ha-6564']].ravel()[uniq_indx]
    data['EW'] = hdu['EMLINE_GEW'].data[el['Ha-6564']].ravel()[uniq_indx]
    data['EW_ERR'] = numpy.ma.power(hdu['EMLINE_GEW_IVAR'].data[el['Ha-6564']].ravel()[uniq_indx],
                                    -0.5).filled(-1)
    data['VEL'] = hdu['EMLINE_GVEL'].data[el['Ha-6564']].ravel()[uniq_indx]
    data['VEL_ERR'] = numpy.ma.power(hdu['EMLINE_GVEL_IVAR'].data[el['Ha-6564']].ravel()[uniq_indx],
                                    -0.5).filled(-1)
    data['SIG'] = hdu['EMLINE_GSIGMA'].data[el['Ha-6564']].ravel()[uniq_indx]
    data['SIG_ERR'] = numpy.ma.power(
                            hdu['EMLINE_GSIGMA_IVAR'].data[el['Ha-6564']].ravel()[uniq_indx],
                            -0.5).filled(-1)
    data['SIG_INST'] = hdu['EMLINE_INSTSIGMA'].data[el['Ha-6564']].ravel()[uniq_indx]
   
    # Deal with the masking
    qual_bm = DAPQualityBitMask()
    data_bm = HalphaBitMask()
    if qual_bm.flagged(hdu[0].header['DAPQUAL'], flag='CRITICAL'):
        data['MASK'] = data_bm.turn_on(data['MASK'], 'CRITCUBE')
    indx = hdu['EMLINE_GFLUX_MASK'].data[el['Ha-6564']].ravel()[uniq_indx] > 0
    if numpy.any(indx):
        data['MASK'][indx] = data_bm.turn_on(data['MASK'][indx], 'BADFLUX')
    indx = hdu['EMLINE_GEW_MASK'].data[el['Ha-6564']].ravel()[uniq_indx] > 0
    if numpy.any(indx):
        data['MASK'][indx] = data_bm.turn_on(data['MASK'][indx], 'BADEW')
    indx = (hdu['EMLINE_GVEL_MASK'].data.ravel()[uniq_indx] > 0) \
            | numpy.logical_not(hdu['EMLINE_GVEL_IVAR'].data.ravel()[uniq_indx] > 0)
    if numpy.any(indx):
        data['MASK'][indx] = data_bm.turn_on(data['MASK'][indx], 'BADVEL')
    indx = (hdu['EMLINE_GSIGMA_MASK'].data.ravel()[uniq_indx] > 0) \
            | numpy.logical_not(hdu['EMLINE_GSIGMA_IVAR'].data.ravel()[uniq_indx] > 0)
    if numpy.any(indx):
        data['MASK'][indx] = data_bm.turn_on(data['MASK'][indx], 'BADSIG')

    return data


class StellarBitMask(BitMask):
    def __init__(self):
        super(StellarBitMask, self).__init__(['CRITCUBE', 'BADFLUX', 'BADVEL', 'BADSIG'])


class StellarDataTable(DataTable):
    def __init__(self, shape=None):
        bm = StellarBitMask()
        _float = numpy.float32
        _int = numpy.int32
        datamodel = dict(PLATE=dict(typ=_int), IFU=dict(typ=_int), X=dict(typ=_int),
                         Y=dict(typ=_int), BINID=dict(typ=_int), NBIN=dict(typ=_int),
                         R=dict(typ=_float), R_RE=dict(typ=_float), R_KPC=dict(typ=_float),
                         GCNT=dict(typ=_float), GSNR=dict(typ=_float), RCHI2=dict(typ=_float),
                         VEL=dict(typ=_float), VEL_ERR=dict(typ=_float), SIG=dict(typ=_float),
                         SIG_ERR=dict(typ=_float), SIG_CORR=dict(typ=_float),
                         MASK=dict(typ=bm.minimum_dtype()))
        keys = list(datamodel.keys())
        super(StellarDataTable, self).__init__(keys, [datamodel[k]['typ'] for k in keys],
                                               shape=shape)


def get_stellar_data(plt, ifu, hdu):
    uniq_bins, uniq_indx, uniq_cnt = map(lambda x : x[1:],
                                         numpy.unique(hdu['BINID'].data[1].ravel(),
                                                      return_index=True, return_counts=True))
    keep = numpy.logical_not(hdu['BIN_MFLUX_MASK'].data.ravel()[uniq_indx] > 0) \
            | numpy.logical_not(hdu['STELLAR_VEL_MASK'].data.ravel()[uniq_indx] > 0) \
            | numpy.logical_not(hdu['STELLAR_SIGMA_MASK'].data.ravel()[uniq_indx] > 0)
    uniq_bins = uniq_bins[keep]
    uniq_indx = uniq_indx[keep]
    uniq_cnt = uniq_cnt[keep]
    data = StellarDataTable(shape=uniq_bins.shape)

    data['PLATE'] = plt
    data['IFU'] = ifu
    data['Y'], data['X'] = numpy.unravel_index(uniq_indx, hdu['BINID'].shape[1:])
    data['BINID'] = uniq_bins
    data['NBIN'] = uniq_cnt
    data['R'] = hdu['BIN_LWELLCOO'].data[0].ravel()[uniq_indx]
    data['R_RE'] = hdu['BIN_LWELLCOO'].data[1].ravel()[uniq_indx]
    data['R_KPC'] = hdu['BIN_LWELLCOO'].data[2].ravel()[uniq_indx]
    data['GCNT'] = hdu['BIN_MFLUX'].data.ravel()[uniq_indx]
    data['GSNR'] = hdu['BIN_SNR'].data.ravel()[uniq_indx]
    data['RCHI2'] = hdu['STELLAR_FOM'].data[2].ravel()[uniq_indx]
    data['VEL'] = hdu['STELLAR_VEL'].data.ravel()[uniq_indx]
    data['VEL_ERR'] = numpy.ma.power(hdu['STELLAR_VEL_IVAR'].data.ravel()[uniq_indx],
                                     -0.5).filled(-1)
    data['SIG'] = hdu['STELLAR_SIGMA'].data.ravel()[uniq_indx]
    data['SIG_ERR'] = numpy.ma.power(hdu['STELLAR_SIGMA_IVAR'].data.ravel()[uniq_indx],
                                     -0.5).filled(-1)
    data['SIG_CORR'] = hdu['STELLAR_SIGMACORR'].data.ravel()[uniq_indx]
   
    # Deal with the masking
    # maps_bm = DAPMapsBitMask()
    qual_bm = DAPQualityBitMask()
    data_bm = StellarBitMask()
    if qual_bm.flagged(hdu[0].header['DAPQUAL'], flag='CRITICAL'):
        data['MASK'] = data_bm.turn_on(data['MASK'], 'CRITCUBE')
    indx = hdu['BIN_MFLUX_MASK'].data.ravel()[uniq_indx] > 0
    if numpy.any(indx):
        data['MASK'][indx] = data_bm.turn_on(data['MASK'][indx], 'BADFLUX')
    indx = (hdu['STELLAR_VEL_MASK'].data.ravel()[uniq_indx] > 0) \
            | numpy.logical_not(hdu['STELLAR_VEL_IVAR'].data.ravel()[uniq_indx] > 0)
    if numpy.any(indx):
        data['MASK'][indx] = data_bm.turn_on(data['MASK'][indx], 'BADVEL')
    indx = (hdu['STELLAR_SIGMA_MASK'].data.ravel()[uniq_indx] > 0) \
            | numpy.logical_not(hdu['STELLAR_SIGMA_IVAR'].data.ravel()[uniq_indx] > 0)
    if numpy.any(indx):
        data['MASK'][indx] = data_bm.turn_on(data['MASK'][indx], 'BADSIG')

    return data


def get_maps_files(daptype='HYB10-MILESHC-MASTARHC2'):
#    return ['/uufs/chpc.utah.edu/common/home/sdss/mangawork/manga/spectro/analysis/v3_0_1/3.0.1/HYB10-MILESHC-MASTARHC2/7443/12704/manga-7443-12704-MAPS-HYB10-MILESHC-MASTARHC2.fits.gz', '/uufs/chpc.utah.edu/common/home/sdss/mangawork/manga/spectro/analysis/v3_0_1/3.0.1/HYB10-MILESHC-MASTARHC2/7443/12701/manga-7443-12701-MAPS-HYB10-MILESHC-MASTARHC2.fits.gz', '/uufs/chpc.utah.edu/common/home/sdss/mangawork/manga/spectro/analysis/v3_0_1/3.0.1/HYB10-MILESHC-MASTARHC2/7443/12703/manga-7443-12703-MAPS-HYB10-MILESHC-MASTARHC2.fits.gz', '/uufs/chpc.utah.edu/common/home/sdss/mangawork/manga/spectro/analysis/v3_0_1/3.0.1/HYB10-MILESHC-MASTARHC2/7443/12705/manga-7443-12705-MAPS-HYB10-MILESHC-MASTARHC2.fits.gz', '/uufs/chpc.utah.edu/common/home/sdss/mangawork/manga/spectro/analysis/v3_0_1/3.0.1/HYB10-MILESHC-MASTARHC2/7443/1901/manga-7443-1901-MAPS-HYB10-MILESHC-MASTARHC2.fits.gz', '/uufs/chpc.utah.edu/common/home/sdss/mangawork/manga/spectro/analysis/v3_0_1/3.0.1/HYB10-MILESHC-MASTARHC2/7443/1902/manga-7443-1902-MAPS-HYB10-MILESHC-MASTARHC2.fits.gz', '/uufs/chpc.utah.edu/common/home/sdss/mangawork/manga/spectro/analysis/v3_0_1/3.0.1/HYB10-MILESHC-MASTARHC2/7443/12702/manga-7443-12702-MAPS-HYB10-MILESHC-MASTARHC2.fits.gz']

    root = os.path.join(defaults.dap_analysis_path(dapver='3.0.1'), daptype)
    return glob.glob(os.path.join(root, '*', '*', '*MAPS*fits.gz'))
    

def main():
    maps_files = get_maps_files()
    nmaps = len(maps_files)
    nbin = 100

    j = 0
    stellar_data = StellarDataTable()
    halpha_data = HalphaDataTable()
    stellar_root = 'stellar_data_mpl10'
    halpha_root = 'halpha_data_mpl10'
    for i in range(nmaps):
        plt, ifu = map(lambda x : int(x), os.path.split(maps_files[i])[1].split('-')[1:3])
        with fits.open(maps_files[i]) as hdu:
            stellar_data.append(get_stellar_data(plt, ifu, hdu))
            halpha_data.append(get_halpha_data(plt, ifu, hdu))
        print('{0}/{1} : {2} : {3}'.format(i+1, nmaps, stellar_data.size, halpha_data.size),
              end='\r')
        j += 1
        if j % nbin == 0 or i == nmaps-1:
            k = i // nbin
            #print(nmaps, nbin, i, j, k)
            ofile = '{0}_{1}.fits.gz'.format(stellar_root, str(k).zfill(3))
            DAPFitsUtil.write(fits.HDUList([fits.PrimaryHDU(),
                                            stellar_data.to_hdu(name='STELLAR')]),
                              ofile, clobber=True)
            ofile = '{0}_{1}.fits.gz'.format(halpha_root, str(k).zfill(3))
            DAPFitsUtil.write(fits.HDUList([fits.PrimaryHDU(),
                                            halpha_data.to_hdu(name='HALPHA')]),
                              ofile, clobber=True)
            stellar_data = StellarDataTable()
            halpha_data = HalphaDataTable()
    print('{0}/{0} : {1} : {2}'.format(nmaps, stellar_data.size, halpha_data.size))

#    ofile = 'stellar_data_mpl10.fits.gz'
#    DAPFitsUtil.write(fits.HDUList([fits.PrimaryHDU(), stellar_data.to_hdu(name='STELLAR')]),
#                      ofile, clobber=True)
#
#    ofile = 'halpha_data_mpl10.fits.gz'
#    DAPFitsUtil.write(fits.HDUList([fits.PrimaryHDU(), halpha_data.to_hdu(name='HALPHA')]),
#                      ofile, clobber=True)
   
if __name__ == '__main__':
    main()

