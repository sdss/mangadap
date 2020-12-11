
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
from mangadap.proc.stellarcontinuummodel import StellarContinuumModelBitMask

from mangadap.tests.util import data_test_file


def test_ppxffit():
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

    # Test the results

    # Rejected pixels
    assert numpy.sum(ppxf.bitmask.flagged(fit_mask, flag='PPXF_REJECT')) == 119, \
                'Different number of rejected pixels'

    # Unable to fit
    assert numpy.array_equal(ppxf.bitmask.flagged_bits(fit_par['MASK'][5]), ['NO_FIT']), \
                'Expected NO_FIT in 6th spectrum'

    # Number of used templates
    assert numpy.array_equal(numpy.sum(numpy.absolute(fit_par['TPLWGT']) > 1e-10, axis=1),
                             [12, 13, 17, 15, 15,  0,  8, 12]), \
                'Different number of templates with non-zero weights'

    # Number of additive coefficients
    assert fit_par['ADDCOEF'].shape[1] == 9, 'Incorrect number of additive coefficients'

    # No multiplicative coefficients
    assert numpy.all(fit_par['MULTCOEF'] == 0), \
                'No multiplicative coefficients should exist'

    # Kinematics and errors
    assert numpy.all(numpy.absolute(fit_par['KIN'] -
                        numpy.array([[ 14880.7, 292.9], [ 15053.4, 123.2],
                                     [ 14787.5, 236.4], [  8291.8, 169.7],
                                     [  9261.4, 202.7], [     0.0,   0.0],
                                     [  5123.5,  63.8], [  5455.6,  51.8]])) < 0.1), \
                'Kinematics are too different'

    assert numpy.all(numpy.absolute(fit_par['KINERR'] -
                        numpy.array([[2.0,1.9], [1.5,1.7], [ 2.4, 2.4], [2.2,2.3],
                                     [1.1,1.1], [0.0,0.0], [26.1,30.8], [4.7,7.5]])) < 0.1), \
                'Kinematic errors are too different'

    # Velocity dispersion corrections
    assert numpy.all(numpy.absolute(fit_par['SIGMACORR_SRES'] -
                        numpy.array([23.5, 10.1, 27.3, 38.7, 22.3,  0.0, 63.8, 23.8])) < 0.1), \
                'SRES corrections are too different'

    assert numpy.all(numpy.absolute(fit_par['SIGMACORR_EMP'] -
                        numpy.array([22.6,  0.0, 26.0, 38.2, 18.0,  0.0, 70.1,  0.0])) < 0.1), \
                'EMP corrections are too different'

    # Figures of merit
    assert numpy.all(numpy.absolute(fit_par['RCHI2'] -
                        numpy.array([ 1.94, 1.18, 1.40, 1.53, 2.50, 0.00, 1.06, 0.86])) < 0.01), \
                'Reduced chi-square too different'

    assert numpy.all(numpy.absolute(fit_par['RMS'] -
                        numpy.array([0.033, 0.019, 0.034, 0.023, 0.046, 0.000, 0.015, 0.015]))
                     < 0.001), 'RMS too different'

    assert numpy.all(numpy.absolute(fit_par['FRMS'] -
                        numpy.array([0.018, 0.023, 0.023, 0.032, 0.018, 0.000, 33.577, 0.148]))
                     < 0.001), 'Fractional RMS too different'

    assert numpy.all(numpy.absolute(fit_par['RMSGRW'][:,2] -
                        numpy.array([0.067, 0.037, 0.068, 0.046, 0.093, 0.000, 0.029, 0.027]))
                     < 0.001), 'Median absolute residual too different'

