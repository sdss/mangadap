
import pytest
import os

import numpy
from scipy import interpolate
from astropy.io import fits
import astropy.constants

from mangadap.drpfits import DRPFits, DRPFitsBitMask

from mangadap.par.artifactdb import ArtifactDB
from mangadap.par.emissionlinedb import EmissionLineDB

from mangadap.util.pixelmask import SpectralPixelMask

from mangadap.proc.templatelibrary import TemplateLibrary
from mangadap.proc.ppxffit import PPXFFit
from mangadap.proc.stellarcontinuummodel import StellarContinuumModelBitMask

from mangadap.tests.util import data_file

import warnings
warnings.simplefilter("ignore", UserWarning)
warnings.simplefilter("ignore", RuntimeWarning)

def test_ppxffit():
    # Read the data
    specfile = data_file('MaNGA_test_spectra.fits.gz')
    hdu = fits.open(specfile)
    drpbm = DRPFitsBitMask()
    flux = numpy.ma.MaskedArray(hdu['FLUX'].data,
                                mask=drpbm.flagged(hdu['MASK'].data, DRPFits.do_not_fit_flags()))
    ferr = numpy.ma.power(hdu['IVAR'].data, -0.5)
    flux[ferr.mask] = numpy.ma.masked
    ferr[flux.mask] = numpy.ma.masked
    nspec = flux.shape[0]

    # Instantiate the template libary
    velscale_ratio = 4
    tpl = TemplateLibrary('MILESHC', match_to_drp_resolution=False, velscale_ratio=velscale_ratio,
                          spectral_step=1e-4, log=True, hardcopy=False)
    tpl_sres = numpy.mean(tpl['SPECRES'].data, axis=0)

    # Get the pixel mask
    pixelmask = SpectralPixelMask(artdb=ArtifactDB.from_key('BADSKY'),
                                  emldb=EmissionLineDB.from_key('ELPFULL'))

    # Instantiate the fitting class
    ppxf = PPXFFit(StellarContinuumModelBitMask())

    # Perform the fit
    fit_wave, fit_flux, fit_mask, fit_par \
        = ppxf.fit(tpl['WAVE'].data.copy(), tpl['FLUX'].data.copy(), hdu['WAVE'].data, flux, ferr,
                   hdu['Z'].data, numpy.full(nspec, 100.), iteration_mode='no_global_wrej',
                   reject_boxcar=100, ensemble=False, velscale_ratio=velscale_ratio,
                   mask=pixelmask, matched_resolution=False, tpl_sres=tpl_sres,
                   obj_sres=hdu['SRES'].data, degree=8, moments=2)

    # Test the results

    # Rejected pixels
    assert numpy.sum(ppxf.bitmask.flagged(fit_mask, flag='PPXF_REJECT')) == 124, \
                'Different number of rejected pixels'

    # Unable to fit
    assert numpy.array_equal(ppxf.bitmask.flagged_bits(fit_par['MASK'][5]), ['NO_FIT']), \
                'Expected NO_FIT in 6th spectrum'

    # Number of used templates
    assert numpy.array_equal(numpy.sum(fit_par['TPLWGT'] > 0, axis=1),
                             [13, 16, 17, 14, 17,  0,  8, 12]), \
                'Different number of templates with non-zero weights'

    # Number of additive coefficients
    assert fit_par['ADDCOEF'].shape[1] == 9, 'Incorrect number of additive coefficients'

    # No multiplicative coefficients
    assert numpy.all(fit_par['MULTCOEF'] == 0), \
                'No multiplicative coefficients should exist'

    # Kinematics and errors
    assert numpy.allclose(fit_par['KIN'], numpy.array([[14881.7722996 ,   291.85468616],
                                                       [15053.5428018 ,   123.07989289],
                                                       [14788.28367965,   235.8053923 ],
                                                       [ 8293.07541034,   170.86985816],
                                                       [ 9261.89551736,   201.84993174],
                                                       [    0.        ,     0.        ],
                                                       [ 5126.52362297,    66.93778845],
                                                       [ 5458.11430831,    51.92792211]]),
                          rtol=0.0, atol=1e-2), 'Kinematics are different'
    assert numpy.allclose(fit_par['KINERR'], numpy.array([[ 1.9709821 ,  1.85886101],
                                                          [ 1.44396882,  1.65099111],
                                                          [ 2.39500413,  2.34584906],
                                                          [ 2.13102942,  2.26898814],
                                                          [ 1.10522357,  1.11307583],
                                                          [ 0.        ,  0.        ],
                                                          [25.50804271, 30.41783466],
                                                          [ 4.46678658,  7.40489565]]),
                          rtol=0.0, atol=1e-2), 'Kinematic errors are different'

    # Velocity dispersion corrections
    assert numpy.allclose(fit_par['SIGMACORR_SRES'],
                          numpy.array([24.03388208, 11.76975482, 27.80845511, 38.76085538,
                                       23.00081717,  0.        , 63.36255328, 24.45657423]),
                          rtol=0.0, atol=1e-2), 'SRES corrections are different'
    assert numpy.allclose(fit_par['SIGMACORR_EMP'],
                          numpy.array([22.45805838,  0.        , 25.8235458 , 38.03511325,
                                       18.04372944,  0.        , 69.11989189,  0.        ]),
                          rtol=0.0, atol=1e-2), 'EMP corrections are different'

    # Figures of merit
    assert numpy.allclose(fit_par['RCHI2'],
                          numpy.array([ 1.95946775,  1.16851502,  1.48772154,  1.55023733,
                                        2.5465503 ,  0.        ,  1.06056225,  0.87039315]),
                          rtol=0.0, atol=1e-4), 'Reduced chi-square different'

    assert numpy.allclose(fit_par['RMS'],
                          numpy.array([0.0337455 , 0.0189093 , 0.03541922, 0.02358469,
                                       0.04709133, 0.        , 0.01608275, 0.01622584]),
                          rtol=0.0, atol=1e-4), 'RMS different'

    assert numpy.allclose(fit_par['FRMS'],
                          numpy.array([0.01892738, 0.02500785, 0.02596966, 0.03601668,
                                       0.01834503, 0.        , 8.59674419, 0.18296954]),
                          rtol=0.0, atol=1e-4), 'Fractional RMS different'

    assert numpy.allclose(fit_par['RMSGRW'][:,2],

                          numpy.array([0.06759905, 0.03682847, 0.0696026 , 0.04725984,
                                       0.09604793, 0.        , 0.03093494, 0.02879866]),
                          rtol=0.0, atol=1e-4), 'Median absolute residual different'

