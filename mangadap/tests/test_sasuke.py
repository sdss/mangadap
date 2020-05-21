
import pytest
import os

from IPython import embed

import numpy
from scipy import interpolate
from astropy.io import fits
import astropy.constants

from mangadap.datacube import MaNGADataCube

from mangadap.par.artifactdb import ArtifactDB
from mangadap.par.emissionlinedb import EmissionLineDB

from mangadap.util.drpfits import DRPFitsBitMask
from mangadap.util.pixelmask import SpectralPixelMask

from mangadap.proc.templatelibrary import TemplateLibrary
from mangadap.proc.ppxffit import PPXFFit
from mangadap.proc.stellarcontinuummodel import StellarContinuumModel, StellarContinuumModelBitMask

from mangadap.tests.util import data_test_file

from mangadap.proc.sasuke import Sasuke
from mangadap.proc.emissionlinemodel import EmissionLineModelBitMask

import warnings
warnings.simplefilter("ignore", UserWarning)
warnings.simplefilter("ignore", RuntimeWarning)

def test_sasuke():
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
                                  emldb=EmissionLineDB.from_key('ELPFULL'))

    # Instantiate the fitting class
    ppxf = PPXFFit(StellarContinuumModelBitMask())

    # Perform the fit
    sc_wave, sc_flux, sc_mask, sc_par \
        = ppxf.fit(tpl['WAVE'].data.copy(), tpl['FLUX'].data.copy(), hdu['WAVE'].data, flux, ferr,
                   hdu['Z'].data, numpy.full(nspec, 100.), iteration_mode='no_global_wrej',
                   reject_boxcar=100, ensemble=False, velscale_ratio=velscale_ratio,
                   mask=pixelmask, matched_resolution=False, tpl_sres=tpl_sres,
                   obj_sres=hdu['SRES'].data, degree=8, moments=2)

    # Mask the 5577 sky line
    pixelmask = SpectralPixelMask(artdb=ArtifactDB.from_key('BADSKY'))

    # Read the emission line fitting database
    emldb = EmissionLineDB.from_key('ELPMILES')
    assert emldb['name'][18] == 'Ha', 'Emission-line database names or ordering changed'

    # Instantiate the fitting class
    emlfit = Sasuke(EmissionLineModelBitMask())

    # Perform the fit
    el_wave, model, el_flux, el_mask, el_fit, el_par \
            = emlfit.fit(emldb, hdu['WAVE'].data, flux, obj_ferr=ferr, obj_mask=pixelmask,
                         obj_sres=hdu['SRES'].data, guess_redshift=hdu['Z'].data,
                         guess_dispersion=numpy.full(nspec, 100.), reject_boxcar=101,
                         stpl_wave=tpl['WAVE'].data, stpl_flux=tpl['FLUX'].data,
                         stpl_sres=tpl_sres, stellar_kinematics=sc_par['KIN'],
                         etpl_sinst_mode='offset', etpl_sinst_min=10.,
                         velscale_ratio=velscale_ratio, matched_resolution=False) #, # mdegree=8,
                         #plot=True)

    # Test the results

    # Rejected pixels
    assert numpy.sum(emlfit.bitmask.flagged(el_mask, flag='PPXF_REJECT')) == 277, \
                'Different number of rejected pixels'

    # Unable to fit
    assert numpy.array_equal(emlfit.bitmask.flagged_bits(el_fit['MASK'][5]), ['NO_FIT']), \
                'Expected NO_FIT in 6th spectrum'

    # Number of used templates
    assert numpy.array_equal(numpy.sum(el_fit['TPLWGT'] > 0, axis=1),
                             [51, 50, 52, 48, 36,  0, 42, 48]), \
                'Different number of templates with non-zero weights'

    # No additive coefficients
    assert numpy.all(el_fit['ADDCOEF'] == 0), \
                'No additive coefficients should exist'

    # No multiplicative coefficients
    assert numpy.all(el_fit['MULTCOEF'] == 0), \
                'No multiplicative coefficients should exist'

    # Fit statistics
    assert numpy.allclose(el_fit['RCHI2'],
                          numpy.array([ 2.34441982,  1.21510401,  1.57781064,  1.87956651,
                                        3.19511594,  0.        ,  1.05074234,  0.88393101]),
                          rtol=0.0, atol=1e-4), 'Reduced chi-square different'

    assert numpy.allclose(el_fit['RMS'],
                    numpy.array([ 0.03621057,  0.01914189,  0.0366114 ,  0.02456429,
                                  0.05073385,  0.        ,  0.01283324,  0.01243996]),
                          rtol=0.0, atol=1e-4), 'RMS different'

    assert numpy.allclose(el_fit['FRMS'],
                    numpy.array([ 0.02097277,  0.02650269,  0.02718998,  0.03482793,
                                  0.01856961,  0.        ,  1.12984714,  0.10905279]),
                          rtol=0.0, atol=1e-4), 'Fractional RMS different'

    assert numpy.allclose(el_fit['RMSGRW'][:,2],
                          numpy.array([ 0.07079034,  0.03748393,  0.07066016,  0.04784343,
                                        0.10093365,  0.        ,  0.02707473,  0.02441517]),
                          rtol=0.0, atol=1e-4), 'Median absolute residual different'

    # All lines should have the same velocity
    assert numpy.all(numpy.all(el_par['KIN'][:,:,0] == el_par['KIN'][:,None,0,0], axis=1)), \
                'All velocities should be the same'

    # Test velocity values
    assert numpy.allclose(el_par['KIN'][:,0,0],
                numpy.array([ 14708.51501137,  14885.9407429 ,  14767.12496209,   8159.45965603,
                               9258.7021507 ,      0.        ,   5131.13183361,   5433.87569423]),
                          rtol=0.0, atol=1e-2), 'Velocities are different'

    # H-alpha dispsersions
    # TODO: Need some better examples!
    assert numpy.allclose(el_par['KIN'][:,18,1],
                numpy.array([ 1000.47576421,  1000.47576421,   224.65658093,   114.82127207,
                               170.88868189,     0.        ,    81.28290423,    50.52682452]),
                          rtol=0.0, atol=1e-1), 'H-alpha dispersions are different'
