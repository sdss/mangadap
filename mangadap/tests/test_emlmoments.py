
import pytest
import os

from IPython import embed

import numpy
from scipy import interpolate
from astropy.io import fits
import astropy.constants

from mangadap.datacube import MaNGADataCube

from mangadap.par.artifactdb import ArtifactDB
from mangadap.par.emissionmomentsdb import EmissionMomentsDB
from mangadap.par.emissionlinedb import EmissionLineDB

from mangadap.util.drpfits import DRPFitsBitMask
from mangadap.util.pixelmask import SpectralPixelMask

from mangadap.proc.templatelibrary import TemplateLibrary
from mangadap.proc.bandpassfilter import emission_line_equivalent_width
from mangadap.proc.ppxffit import PPXFFit
from mangadap.proc.stellarcontinuummodel import StellarContinuumModel, StellarContinuumModelBitMask
from mangadap.proc.emissionlinemoments import EmissionLineMoments, EmissionLineMomentsBitMask

from mangadap.tests.util import data_test_file

import warnings
warnings.simplefilter("ignore", UserWarning)
warnings.simplefilter("ignore", RuntimeWarning)

def test_moments():

    # Read the data
    specfile = data_test_file('MaNGA_test_spectra.fits.gz')
    hdu = fits.open(specfile)
    drpbm = DRPFitsBitMask()
    flux = numpy.ma.MaskedArray(hdu['FLUX'].data, mask=drpbm.flagged(hdu['MASK'].data,
                                                                MaNGADataCube.do_not_fit_flags()))
    ferr = numpy.ma.power(hdu['IVAR'].data, -0.5)
    flux[ferr.mask] = numpy.ma.masked
    ferr[flux.mask] = numpy.ma.masked
    nspec = flux.shape[0]

    # Read the database that define the emission lines and passbands
    momdb = EmissionMomentsDB.from_key('ELBMILES')
    
    # Measure the moments
    elmombm = EmissionLineMomentsBitMask()
    elmom = EmissionLineMoments.measure_moments(momdb, hdu['WAVE'].data, flux,
                                                redshift=hdu['Z'].data, bitmask=elmombm)

    # Measure the EW based on the moments
    include_band = numpy.array([numpy.invert(momdb.dummy)]*nspec) \
                        & numpy.invert(elmombm.flagged(elmom['MASK'],
                                                       flag=['BLUE_EMPTY', 'RED_EMPTY']))
    line_center = (1.0+hdu['Z'].data)[:,None]*momdb['restwave'][None,:]
    elmom['BMED'], elmom['RMED'], pos, elmom['EWCONT'], elmom['EW'], elmom['EWERR'] \
            = emission_line_equivalent_width(hdu['WAVE'].data, flux, momdb['blueside'],
                                             momdb['redside'], line_center, elmom['FLUX'],
                                             redshift=hdu['Z'].data,
                                             line_flux_err=elmom['FLUXERR'],
                                             include_band=include_band)

    # Check the flags
    reference = {'BLUE_INCOMP': 21, 'MAIN_JUMP': 0, 'UNDEFINED_MOM2': 46, 'JUMP_BTWN_SIDEBANDS': 0,
                 'RED_JUMP': 0, 'DIVBYZERO': 0, 'NO_ABSORPTION_CORRECTION': 176, 'RED_EMPTY': 21,
                 'UNDEFINED_BANDS': 8, 'DIDNOTUSE': 0, 'UNDEFINED_MOM1': 0, 'FORESTAR': 0,
                 'NON_POSITIVE_CONTINUUM': 0, 'LOW_SNR': 0, 'MAIN_EMPTY': 21, 'BLUE_JUMP': 0,
                 'RED_INCOMP': 21, 'MAIN_INCOMP': 21, 'BLUE_EMPTY': 21}
    assert numpy.all([ reference[k] == numpy.sum(elmombm.flagged(elmom['MASK'], flag=k))
                            for k in elmombm.keys()]), 'Number of flagged measurements changed'

    # Check that the values are finite
    assert numpy.all([ numpy.all(numpy.isfinite(elmom[n])) for n in elmom.dtype.names]), \
                        'Found non-finite values in output'

    # Check the band definitions
    assert numpy.all(numpy.equal(elmom['REDSHIFT'], hdu['Z'].data)), 'Redshift changed'
    assert numpy.all(numpy.isclose(numpy.mean(momdb['blueside'], axis=1)[None,:],
                                   elmom['BCEN']/(1+hdu['Z'].data[:,None]))
                        | elmombm.flagged(elmom['MASK'], flag='UNDEFINED_BANDS')), \
                'Blue passband center incorrect'
    assert numpy.all(numpy.isclose(numpy.mean(momdb['redside'], axis=1)[None,:],
                                   elmom['RCEN']/(1+hdu['Z'].data[:,None]))
                        | elmombm.flagged(elmom['MASK'], flag='UNDEFINED_BANDS')), \
                'Red passband center incorrect'

    # Check the values
    assert numpy.allclose(elmom['FLUX'][0],
                          numpy.array([ -0.83366296,   0.        ,  -0.7368989 ,  -6.84760392,
                                        -5.8392653 ,  -3.84394899,  -9.63158548, -10.1459227 ,
                                        -1.86639944,   0.19851703,   0.04831539,  -5.58001859,
                                        0.86652478,  -1.3277138 ,   4.48556862,   0.12541773,
                                        -1.37675776,   1.14456948,  -1.41808526,   2.48743805,
                                        -0.31254732,   0.04046428]),
                          rtol=0.0, atol=1e-2), 'Fluxes changed'
    assert numpy.allclose(elmom['MOM1'][0],
                    numpy.array([15403.91870501,     0.        , 13866.58355013, 14816.45834376,
                                 14861.90408263, 14545.21106265, 14929.76054479, 14774.62443577,
                                 14943.56586856, 13010.07824437, 15933.25294444, 14918.25984067,
                                 14425.53398781, 15207.53998774, 14803.71786274, 14160.66542001,
                                 14720.66321017, 14706.89675211, 14880.91017052, 14901.49219165,
                                 14880.79548007, 15615.43369812]),
                          rtol=0.0, atol=1e-1), '1st moments changed'
    assert numpy.allclose(elmom['MOM2'][0],
                          numpy.array([  0.        ,   0.        ,   0.        , 439.76305578,
                                       479.32501708, 325.96571646, 348.71402151, 362.29430475,
                                       128.76827924,   0.        ,   0.        , 322.61461489,
                                       268.26542796,  27.14271982, 259.24977286,   0.        ,
                                       181.94055378, 129.62366078, 147.48288905, 225.76488299,
                                       132.57819153,   0.        ]),
                          rtol=0.0, atol=1e-1), '2nd moments changed'
    assert numpy.allclose(elmom['EW'][0],
                    numpy.array([-0.83148156,  0.        , -0.67854382, -6.65583709, -4.99844209,
                                 -3.06783667, -6.6506484 , -6.86724193, -0.99166185,  0.08843696,
                                 0.01728948, -1.81199184,  0.28592615, -0.46054113,  1.48650809,
                                 0.03822714, -0.40850899,  0.33980593, -0.42043643,  0.73608197,
                                 -0.09406925,  0.01217937]),
                          rtol=0.0, atol=1e-2), 'EW changed'


def test_moments_with_continuum():
    # Read the data
    specfile = data_test_file('MaNGA_test_spectra.fits.gz')
    hdu = fits.open(specfile)
    drpbm = DRPFitsBitMask()
    flux = numpy.ma.MaskedArray(hdu['FLUX'].data, mask=drpbm.flagged(hdu['MASK'].data,
                                                                MaNGADataCube.do_not_fit_flags()))
    ferr = numpy.ma.power(hdu['IVAR'].data, -0.5)
    flux[ferr.mask] = numpy.ma.masked
    ferr[flux.mask] = numpy.ma.masked
    nspec = flux.shape[0]

    # Instantiate the template libary
    velscale_ratio = 4
    tpl = TemplateLibrary('MILESHC', match_resolution=False, velscale_ratio=velscale_ratio,
                          spectral_step=1e-4, log=True, hardcopy=False)
    tpl_sres = numpy.mean(tpl['SPECRES'].data, axis=0)

    # Get the pixel mask
    pixelmask = SpectralPixelMask(artdb=ArtifactDB.from_key('BADSKY'),
                                  emldb=EmissionLineDB.from_key('ELPSCMSK'))

    # Instantiate the fitting class
    ppxf = PPXFFit(StellarContinuumModelBitMask())

    # Perform the fit
    fit_wave, fit_flux, fit_mask, fit_par \
        = ppxf.fit(tpl['WAVE'].data.copy(), tpl['FLUX'].data.copy(), hdu['WAVE'].data, flux, ferr,
                   hdu['Z'].data, numpy.full(nspec, 100.), iteration_mode='no_global_wrej',
                   reject_boxcar=100, ensemble=False, velscale_ratio=velscale_ratio,
                   mask=pixelmask, matched_resolution=False, tpl_sres=tpl_sres,
                   obj_sres=hdu['SRES'].data, degree=8, moments=2)

    # Remask the continuum fit
    sc_continuum = StellarContinuumModel.reset_continuum_mask_window(
                            numpy.ma.MaskedArray(fit_flux, mask=fit_mask>0))

    # Read the database that define the emission lines and passbands
    momdb = EmissionMomentsDB.from_key('ELBMILES')
    
    # Measure the moments
    elmombm = EmissionLineMomentsBitMask()
    elmom = EmissionLineMoments.measure_moments(momdb, hdu['WAVE'].data, flux,
                                                continuum=sc_continuum, redshift=hdu['Z'].data,
                                                bitmask=elmombm)

    # Measure the EW based on the moments
    include_band = numpy.array([numpy.invert(momdb.dummy)]*nspec) \
                        & numpy.invert(elmombm.flagged(elmom['MASK'],
                                                       flag=['BLUE_EMPTY', 'RED_EMPTY']))
    line_center = (1.0+hdu['Z'].data)[:,None]*momdb['restwave'][None,:]
    elmom['BMED'], elmom['RMED'], pos, elmom['EWCONT'], elmom['EW'], elmom['EWERR'] \
            = emission_line_equivalent_width(hdu['WAVE'].data, flux, momdb['blueside'],
                                             momdb['redside'], line_center, elmom['FLUX'],
                                             redshift=hdu['Z'].data,
                                             line_flux_err=elmom['FLUXERR'],
                                             include_band=include_band)

    # Check the flags
    reference = {'BLUE_INCOMP': 21, 'MAIN_JUMP': 0, 'UNDEFINED_MOM2': 42, 'JUMP_BTWN_SIDEBANDS': 0,
                 'RED_JUMP': 0, 'DIVBYZERO': 0, 'NO_ABSORPTION_CORRECTION': 0, 'RED_EMPTY': 21,
                 'UNDEFINED_BANDS': 8, 'DIDNOTUSE': 0, 'UNDEFINED_MOM1': 0, 'FORESTAR': 0,
                 'NON_POSITIVE_CONTINUUM': 0, 'LOW_SNR': 0, 'MAIN_EMPTY': 21, 'BLUE_JUMP': 0,
                 'RED_INCOMP': 21, 'MAIN_INCOMP': 21, 'BLUE_EMPTY': 21}
    assert numpy.all([ reference[k] == numpy.sum(elmombm.flagged(elmom['MASK'], flag=k))
                            for k in elmombm.keys()]), 'Number of flagged measurements changed'

    # Check that the values are finite
    assert numpy.all([ numpy.all(numpy.isfinite(elmom[n])) for n in elmom.dtype.names]), \
                        'Found non-finite values in output'

    # Check the band definitions
    assert numpy.all(numpy.equal(elmom['REDSHIFT'], hdu['Z'].data)), 'Redshift changed'
    assert numpy.all(numpy.isclose(numpy.mean(momdb['blueside'], axis=1)[None,:],
                                   elmom['BCEN']/(1+hdu['Z'].data[:,None]))
                        | elmombm.flagged(elmom['MASK'], flag='UNDEFINED_BANDS')), \
                'Blue passband center incorrect'
    assert numpy.all(numpy.isclose(numpy.mean(momdb['redside'], axis=1)[None,:],
                                   elmom['RCEN']/(1+hdu['Z'].data[:,None]))
                        | elmombm.flagged(elmom['MASK'], flag='UNDEFINED_BANDS')), \
                'Red passband center incorrect'

    # Check the values
    assert numpy.all(numpy.absolute(elmom['FLUX'][0] -
                        numpy.array([ 0.63,  0.00,  0.22, -1.32, -0.88, -0.68, -0.44, -0.13, -1.14,
                                     -0.07, -0.11,  0.01,  0.38,  0.73,  0.71,  0.44,  0.08,  0.74,
                                      1.30,  2.34, 0.55,  0.44])) < 0.01), 'Fluxes too different'

    assert numpy.all(numpy.absolute(elmom['MOM1'][0] -
                        numpy.array([ 14682.6,      0.0, 14843.2, 14865.8, 14890.4, 14404.7,
                                      14208.6,  12376.0, 14662.5, 14148.5, 15804.1, 17948.4,
                                      14874.5,  14774.9, 14840.5, 14746.0, 15093.1, 14857.8,
                                      14839.0,  14840.2, 14876.0, 14859.5])) < 0.1), \
                    '1st moments too different'

    assert numpy.all(numpy.absolute(elmom['MOM2'][0] -
                        numpy.array([322.2,   0.0, 591.4, 436.4, 474.6,   0.0,   0.0,   0.0,
                                     364.6,   0.0,   0.0,   0.0, 289.1, 226.9, 282.6, 283.8,
                                     227.0, 207.7, 207.7, 253.6, 197.0, 212.4])) < 0.1), \
                    '2nd moments too different'

    assert numpy.all(numpy.absolute(elmom['EW'][0] -
                        numpy.array([ 0.63,  0.00,  0.20, -1.28, -0.76, -0.54, -0.30, -0.09,
                                     -0.61, -0.03, -0.04,  0.00,  0.13,  0.25,  0.24,  0.13,
                                      0.02,  0.22,  0.38,  0.69,  0.17,  0.13])) < 0.01), \
                    'EW too different'

