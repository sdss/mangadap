from IPython import embed

import numpy as np

from astropy.io import fits
from mangadap.proc import spatialbinning


def test_voronoipar():
    # Single parameter:
    par = spatialbinning.VoronoiBinningPar(target_snr=10., covar_mode='calibrate', covar_par=1.62)

    hdr = par.toheader(fits.Header())
    assert hdr['BINCOV'] == 'calibrate', 'mode is wrong'
    assert eval(hdr['NCALIB']) == 1.62, 'parameter is wrong'

    _par = spatialbinning.VoronoiBinningPar()
    _par.fromheader(hdr)
    assert _par['target_snr'] == par['target_snr'], 'read target S/N incorrectly'
    assert _par['covar_mode'] == par['covar_mode'], 'read mode incorrectly'
    assert _par['covar_par'] == par['covar_par'], 'read parameter incorrectly'

    # Three covar parameters:
    par = spatialbinning.VoronoiBinningPar(target_snr=10., covar_mode='calibrate',
                                           covar_par=[1.62, 1.1, 30])

    hdr = par.toheader(fits.Header())
    assert hdr['BINCOV'] == 'calibrate', 'mode is wrong'
    assert [eval(p) for p in hdr['NCALIB'].split(',')] == par['covar_par'], 'parameters are wrong'

    _par = spatialbinning.VoronoiBinningPar()
    _par.fromheader(hdr)
    assert _par['target_snr'] == par['target_snr'], 'read target S/N incorrectly'
    assert _par['covar_mode'] == par['covar_mode'], 'read mode incorrectly'
    assert _par['covar_par'] == par['covar_par'], 'read parameter incorrectly'


def test_voronoi_bin():

    n = 31
    x, y = np.meshgrid(np.arange(n) - n//2, np.arange(n) - n//2)
    r = np.sqrt(x**2 + y**2)
    img = 100*np.exp(-r/10)
    err = np.sqrt(img)

    binner = spatialbinning.VoronoiBinning()

    no_cov_par = spatialbinning.VoronoiBinningPar(target_snr=7.,
                                                  signal=img.ravel(), noise=err.ravel())
    no_cov_binid = binner.bin_index(x=x.ravel(), y=y.ravel(), par=no_cov_par)
    assert len(np.unique(no_cov_binid)) < x.size, 'Should bin at least a few pixels'

    cov_1p_par = spatialbinning.VoronoiBinningPar(target_snr=7.,
                                                  signal=img.ravel(), noise=err.ravel(),
                                                  covar_mode='calibrate', covar_par=1.62)
    cov_1p_binid = binner.bin_index(x=x.ravel(), y=y.ravel(), par=cov_1p_par)
    assert len(np.unique(cov_1p_binid)) < x.size, 'Should bin at least a few pixels'
    assert len(np.unique(cov_1p_binid)) < len(np.unique(no_cov_binid)), \
            'Covariance should reduce the number of bins'

    cov_3p_par = spatialbinning.VoronoiBinningPar(target_snr=7.,
                                                  signal=img.ravel(), noise=err.ravel(),
                                                  covar_mode='calibrate', covar_par=[1.62,1.1,20])
    cov_3p_binid = binner.bin_index(x=x.ravel(), y=y.ravel(), par=cov_3p_par)
    assert len(np.unique(cov_3p_binid)) < x.size, 'Should bin at least a few pixels'
    assert len(np.unique(cov_3p_binid)) < len(np.unique(no_cov_binid)), \
            'Covariance should reduce the number of bins'
    assert len(np.unique(cov_3p_binid)) < len(np.unique(cov_1p_binid)), \
            '3-parameter covariance calibration should have fewer bins (at least in this case)'
    
