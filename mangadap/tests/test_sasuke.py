from IPython import embed

import numpy
from scipy import special
from astropy.io import fits
import astropy.constants

from ppxf import ppxf, ppxf_util

from mangadap.datacube import MaNGADataCube

from mangadap.par.artifactdb import ArtifactDB
from mangadap.par.emissionlinedb import EmissionLineDB

from mangadap.util.drpbitmask import DRPFitsBitMask
from mangadap.util.pixelmask import SpectralPixelMask

from mangadap.proc.templatelibrary import TemplateLibrary
from mangadap.proc.ppxffit import PPXFFit
from mangadap.proc.stellarcontinuummodel import StellarContinuumModelBitMask

from mangadap.tests.util import data_test_file

from mangadap.proc.sasuke import Sasuke
from mangadap.proc.emissionlinemodel import EmissionLineModelBitMask

from mangadap.util.sampling import spectrum_velocity_scale
from mangadap.proc.emissionlinetemplates import EmissionLineTemplates
from mangadap.contrib.xjmc import ppxf_tied_parameters

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
                                  emldb=EmissionLineDB.from_key('ELPSCMSK'))

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
                         velscale_ratio=velscale_ratio) #, matched_resolution=False)

    # Rejected pixels
    # TODO: Commented assertions are too system dependent.  Change these to warnings?
#    assert numpy.sum(emlfit.bitmask.flagged(el_mask, flag='PPXF_REJECT')) == 266, \
#                'Different number of rejected pixels'

    # Unable to fit
    assert numpy.array_equal(emlfit.bitmask.flagged_bits(el_fit['MASK'][5]), ['NO_FIT']), \
                'Expected NO_FIT in 6th spectrum'

    # No *attempted* fits should fail
    assert numpy.sum(emlfit.bitmask.flagged(el_fit['MASK'], flag='FIT_FAILED')) == 0, \
                'Fits should not fail'

    # Number of used templates
#    assert numpy.array_equal(numpy.sum(numpy.absolute(el_fit['TPLWGT']) > 1e-10, axis=1),
#                             [25, 22, 34, 32, 27,  0, 16, 22]), \
#                'Different number of templates with non-zero weights'

    # No additive coefficients
    assert numpy.all(el_fit['ADDCOEF'] == 0), \
                'No additive coefficients should exist'

    # No multiplicative coefficients
    assert numpy.all(el_fit['MULTCOEF'] == 0), \
                'No multiplicative coefficients should exist'

    # Fit statistics
#    assert numpy.all(numpy.absolute(el_fit['RCHI2'] - 
#                                    numpy.array([2.34, 1.22, 1.58, 1.88, 3.20, 0., 1.05, 0.88]))
#                     < 0.02), 'Reduced chi-square are too different'

#    assert numpy.all(numpy.absolute(el_fit['RMS'] -
#                                    numpy.array([0.036, 0.019, 0.036, 0.024, 0.051, 0.000,
#                                                 0.012, 0.012])) < 0.001), 'RMS too different'

#    assert numpy.all(numpy.absolute(el_fit['FRMS'] -
#                                    numpy.array([0.021, 0.025, 0.025, 0.033, 0.018, 0.000,
#                                                 1.052, 0.101])) < 0.001), \
#            'Fractional RMS too different'

#    assert numpy.all(numpy.absolute(el_fit['RMSGRW'][:,2] - 
#                                    numpy.array([0.070, 0.038, 0.071, 0.047, 0.101, 0.000, 0.026,
#                                                 0.024])) < 0.001), \
#            'Median absolute residual too different'

    # All lines should have the same velocity
    assert numpy.all(numpy.all(el_par['KIN'][:,:,0] == el_par['KIN'][:,None,0,0], axis=1)), \
                'All velocities should be the same'

    # Test velocity values
    # TODO: Need some better examples!
#    assert numpy.all(numpy.absolute(el_par['KIN'][:,0,0] -
#                                    numpy.array([14704.9, 14869.3, 14767.1, 8161.9, 9258.7, 0.0,
#                                                  5130.9,  5430.3])) < 0.1), \
#                'Velocities are too different'

#    # H-alpha dispersions
#    assert numpy.all(numpy.absolute(el_par['KIN'][:,18,1] -
#                                    numpy.array([1000.5, 1000.5, 224.7, 124.9, 171.2, 0.0, 81.2,
#                                                   50.0])) < 1e-1), \
#            'H-alpha dispersions are too different'


def test_sasuke_mpl11():
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
    sc_wave, sc_flux, sc_mask, sc_par \
        = ppxf.fit(tpl['WAVE'].data.copy(), tpl['FLUX'].data.copy(), hdu['WAVE'].data, flux, ferr,
                   hdu['Z'].data, numpy.full(nspec, 100.), iteration_mode='no_global_wrej',
                   reject_boxcar=100, ensemble=False, velscale_ratio=velscale_ratio,
                   mask=pixelmask, matched_resolution=False, tpl_sres=tpl_sres,
                   obj_sres=hdu['SRES'].data, degree=8, moments=2)

    # Mask the 5577 sky line
    pixelmask = SpectralPixelMask(artdb=ArtifactDB.from_key('BADSKY'))

    # Read the emission line fitting database
    emldb = EmissionLineDB.from_key('ELPMPL11')
    assert emldb['name'][23] == 'Ha', 'Emission-line database names or ordering changed'

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
                         velscale_ratio=velscale_ratio) #, plot=True) #, matched_resolution=False

    # Rejected pixels
    # TODO: Commented assertions are too system dependent.  Change these to warnings?
#    assert numpy.sum(emlfit.bitmask.flagged(el_mask, flag='PPXF_REJECT')) == 261, \
#                'Different number of rejected pixels'

    # Unable to fit
    assert numpy.array_equal(emlfit.bitmask.flagged_bits(el_fit['MASK'][5]), ['NO_FIT']), \
                'Expected NO_FIT in 6th spectrum'

    # No *attempted* fits should fail
    assert numpy.sum(emlfit.bitmask.flagged(el_fit['MASK'], flag='FIT_FAILED')) == 0, \
                'Fits should not fail'

    # Number of used templates
#    assert numpy.array_equal(numpy.sum(numpy.absolute(el_fit['TPLWGT']) > 1e-10, axis=1),
#                             [24, 22, 36, 35, 31,  0, 17, 23]), \
#                'Different number of templates with non-zero weights'

    # No additive coefficients
    assert numpy.all(el_fit['ADDCOEF'] == 0), \
                'No additive coefficients should exist'

    # No multiplicative coefficients
    assert numpy.all(el_fit['MULTCOEF'] == 0), \
                'No multiplicative coefficients should exist'

    # Fit statistics
#    assert numpy.all(numpy.absolute(el_fit['RCHI2'] - 
#                                    numpy.array([2.33, 1.21, 1.51, 1.88, 3.15, 0., 1.04, 0.87]))
#                     < 0.02), 'Reduced chi-square are too different'

#    assert numpy.all(numpy.absolute(el_fit['RMS'] -
#                                    numpy.array([0.036, 0.019, 0.036, 0.024, 0.050, 0.000,
#                                                 0.012, 0.012])) < 0.001), 'RMS too different'

#    assert numpy.all(numpy.absolute(el_fit['FRMS'] -
#                                    numpy.array([0.020, 0.025, 0.025, 0.033, 0.018, 0.000,
#                                                 1.051, 0.101])) < 0.001), \
#            'Fractional RMS too different'

#    assert numpy.all(numpy.absolute(el_fit['RMSGRW'][:,2] - 
#                                    numpy.array([0.071, 0.037, 0.070, 0.047, 0.099, 0.000, 0.026,
#                                                 0.024])) < 0.001), \
#            'Median absolute residual too different'


    # All (valid) lines should have the same velocity
    masked = emlfit.bitmask.flagged(el_par['MASK'], flag='INSUFFICIENT_DATA')
    assert numpy.all(numpy.all((el_par['KIN'][:,:,0] == el_par['KIN'][:,None,0,0]) | masked,
                               axis=1)), 'All velocities should be the same'

    # Test velocity values
    # TODO: Need some better examples!  This skips over the 4th spectrum because
    # of system-dependent variability in the result.
#    assert numpy.all(numpy.absolute(numpy.append(el_par['KIN'][:3,0,0], el_par['KIN'][4:,0,0]) -
#                                    numpy.array([14655.1, 14390.3, 14768.2, 9259.7, 0.0,
#                                                  5132.6,  5428.7])) < 0.1), \
#                'Velocities are too different'

#    assert numpy.all(numpy.absolute(numpy.array([el_par['KIN'][2,23,0], el_par['KIN'][4,23,0]])
#                                    - numpy.array([14768.2, 9259.7])) < 1.), \
#                'Velocities are too different'

    # H-alpha dispersions
#    assert numpy.all(numpy.absolute(numpy.array([el_par['KIN'][2,23,1], el_par['KIN'][4,23,1]])
#                                    - numpy.array([223.4, 171.2])) < 0.2), \
#                'H-alpha dispersions are too different'


def test_multicomp_basic():
    """ Test a basic multicomponent fit."""

    # Instantiate the template libary
    tpl = TemplateLibrary('MILESHC', match_resolution=False, spectral_step=2e-5)
    tpl_sres = numpy.mean(tpl['SPECRES'].data, axis=0)

    wave = tpl['WAVE'].data.copy()

    velscale = spectrum_velocity_scale(wave)
    emldb = EmissionLineDB(data_test_file('elp_ha_multicomp.par'))

    etpl = EmissionLineTemplates(wave, velscale, emldb=emldb)

    # The broad components are numbered 5 and higher
    stellar_sigma = 150.
    flux = 0.8*ppxf_util.gaussian_filter1d(tpl['FLUX'].data[34], stellar_sigma/velscale)
    narrow_sigma = 100.
    broad_sigma = 500.
#    from matplotlib import pyplot
#    pyplot.plot(wave, flux, color='C0')
    etpl_sigma = []
    for i in range(etpl.ntpl):
        etpl_sigma += [narrow_sigma if etpl.comp[i] < 5 else broad_sigma]
        norm = 20 if etpl.comp[i] < 5 else 80
        line = norm*ppxf_util.gaussian_filter1d(etpl.flux[i], etpl_sigma[-1]/velscale)
        flux += line
#        pyplot.plot(wave, line, color='C3' if etpl.comp[i] < 5 else 'C1')
#    pyplot.plot(wave, flux, color='k')
#    pyplot.show()

    ferr = numpy.full(flux.size, 0.1, dtype=float)
    mask = numpy.zeros(flux.size, dtype=bool)
    sres = tpl_sres.copy()

    templates = numpy.append(tpl['FLUX'].data[34].reshape(1,-1), etpl.flux, axis=0)
    component = numpy.append([0], etpl.comp+1)
    gas_component = numpy.ones(component.size, dtype=bool)
    gas_component[0] = False
    ncomp = numpy.amax(component)+1
    narrow = numpy.ones(ncomp, dtype=bool)
    narrow[0] = False
    narrow[6:] = False
    broad = numpy.logical_not(narrow)
    broad[0] = False
    vgrp = numpy.append([0], etpl.vgrp+1)
    sgrp = numpy.append([0], etpl.sgrp+1)
    moments = numpy.array([-2] + [2]*(ncomp-1))
    gas_comp_sigma = numpy.empty(ncomp-1, dtype=float)
    gas_comp_sigma[etpl.comp] = etpl_sigma
    start_kin = numpy.append([[-50.*velscale, stellar_sigma]],
                             numpy.column_stack(([100.]*(ncomp-1), gas_comp_sigma*1.1)), axis=0)
    A_ineq = numpy.hstack((numpy.zeros((etpl.A_ineq.shape[0], 2), dtype=float),
                           etpl.A_ineq))
    constr_kinem = {'A_ineq': A_ineq, 'b_ineq': etpl.b_ineq}
    tied = ppxf_tied_parameters(component, vgrp, sgrp, moments)

    #from matplotlib import pyplot
    pp = ppxf.ppxf(templates.T, flux[50:-50], ferr[50:-50], velscale, start_kin, moments=moments,
                    degree=-1, mdegree=0, tied=tied, constr_kinem=constr_kinem,
                    component=component, gas_component=gas_component, method='capfit',
                    quiet=True)
    #                quiet=False, plot=True)
    #pyplot.show()

    sol = numpy.array(pp.sol)
    assert numpy.allclose(sol[1,0], sol[2:,0]), 'All gas velocities should be the same'
    assert numpy.all(sol[broad,1] > sol[narrow,1]), 'Constraints not met'
    assert numpy.all(numpy.absolute(sol[1:,0] - start_kin[0,0]) < 1e-2), \
            'Velocities too discrepant from input'
    assert numpy.all(numpy.absolute(sol[narrow,1] - narrow_sigma) < 1e-1), \
            'Narrow dispersion too discrepant from input'
    assert numpy.all(numpy.absolute(sol[broad,1] - broad_sigma) < 1.), \
            'Broad dispersion too discrepant from input'


def test_multicomp_manga():
    # Read the data
    specfile = data_test_file('MaNGA_multicomp_test_spectra.fits.gz')
    hdu = fits.open(specfile)
    drpbm = DRPFitsBitMask()
    flux = numpy.ma.MaskedArray(hdu['FLUX'].data, mask=drpbm.flagged(hdu['MASK'].data,
                                                                MaNGADataCube.do_not_fit_flags()))
#    flux[1:,:] = numpy.ma.masked
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
    sc_wave, sc_flux, sc_mask, sc_par \
        = ppxf.fit(tpl['WAVE'].data.copy(), tpl['FLUX'].data.copy(), hdu['WAVE'].data, flux, ferr,
                   hdu['Z'].data, numpy.full(nspec, 100.), iteration_mode='no_global_wrej',
                   reject_boxcar=100, ensemble=False, velscale_ratio=velscale_ratio,
                   mask=pixelmask, matched_resolution=False, tpl_sres=tpl_sres,
                   obj_sres=hdu['SRES'].data, degree=8, moments=2)

    # Mask the 5577 sky line
    pixelmask = SpectralPixelMask(artdb=ArtifactDB.from_key('BADSKY'))

    # Read the emission line fitting database
    emldb = EmissionLineDB(data_test_file('elp_ha_multicomp_broadv.par'))

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
                         velscale_ratio=velscale_ratio) #, plot=True)

    assert numpy.array_equal(el_fit['MASK'], numpy.zeros(nspec, dtype=el_fit['MASK'].dtype)), \
            'None of the fits should fail!'

    # TODO: Add tests of the results!



