import os
import warnings

import astropy.constants
import numpy
from ppxf import ppxf_util
import pytest

from mangadap.config import defaults
from mangadap.par.emissionlinedb import EmissionLineDB
from mangadap.proc.emissionlinetemplates import EmissionLineTemplates
from mangadap.proc.ppxffit import PPXFFit
from mangadap.proc.templatelibrary import TemplateLibrary
from mangadap.util.constants import DAPConstants
from mangadap.util.sampling import spectrum_velocity_scale


def data_test_file(filename=None):
    root = os.path.join(defaults.dap_data_root(), 'tests')
    return root if filename is None else os.path.join(root, filename)

def remote_data_file(filename=None):
    root = os.path.join(defaults.dap_data_root(), 'remote')
    return root if filename is None else os.path.join(root, filename)

def remote_data_files():
    return ['manga-7815-3702-LINCUBE.fits.gz', 'manga-7815-3702-LINRSS.fits.gz',
            'manga-7815-3702-LOGCUBE.fits.gz', 'manga-7815-3702-LOGRSS.fits.gz']

drp_test_version = 'v3_1_1'
dap_test_version = '3.1.0'
remote_available = all([os.path.isfile(remote_data_file(f)) for f in remote_data_files()])
drpcomplete_available = os.path.isfile(remote_data_file('drpcomplete_{0}.fits'.format(
                                       drp_test_version)))
drpall_available = os.path.isfile(remote_data_file('drpall-{0}.fits'.format(drp_test_version)))

requires_remote = pytest.mark.skipif(not remote_available, reason='Remote data files are missing.')
requires_drpcomplete = pytest.mark.skipif(not drpcomplete_available,
                                          reason='Remote DRPComplete is missing.')
requires_drpall = pytest.mark.skipif(not drpall_available, reason='Remote DRPall is missing.')

if not remote_available:
    warnings.warn('Remote data not available.  Some tests will be skipped.')
if not drpcomplete_available:
    warnings.warn('Remote DRPComplete not available.  Some tests will be skipped.')
if not drpall_available:
    warnings.warn('Remote DRPall not available.  Some tests will be skipped.')


def emline_model_fake_spectrum(emline_match_template_resolution=False):
    """
    Generate a fake spectrum for testing.
    """
    # Instantiate the stellar template libary
    tpl = TemplateLibrary('MILESHC', match_resolution=False, spectral_step=2e-5)
    tpl_sres = numpy.mean(tpl['SPECRES'].data, axis=0)

    wave = tpl['WAVE'].data.copy()
    velscale = spectrum_velocity_scale(wave)

    if emline_match_template_resolution:
        # Set the instrumental resolution to use for the emission lines to be identical to these template spectra
        tpl_wave, tpl_flux, _tpl_sres = PPXFFit.check_templates(
            tpl['WAVE'].data, tpl['FLUX'].data, tpl_sres=tpl_sres, velscale_ratio=1
        )
        R_to_sinst = astropy.constants.c.to('km/s').value / DAPConstants.sig2fwhm
        etpl_sinst = R_to_sinst/_tpl_sres.sres(copy=True)
        min_sinst = numpy.amin(etpl_sinst)
        etpl_sinst_min = 10.
        # Offset such that the minimum over the fitted wavelength
        # range matches the requested minimum
        dsigma_inst = numpy.square(min_sinst)-numpy.square(etpl_sinst_min)
        # Clip to make sure instrumental resolution is real!
        etpl_sinst = numpy.sqrt(numpy.clip(etpl_sinst**2 - dsigma_inst, 0., None))
    else:
        etpl_sinst = velscale

    # Instantiate the gas template library with the same sampling as the stellar
    # library
    emldb = EmissionLineDB.from_key('ELPSTRONG')
    etpl = EmissionLineTemplates(wave, etpl_sinst, emldb=emldb)

    # Construct the fake spectrum
    # - Use one template for the stellar continuum
    stellar_sigma = 150.
    flux = 0.8*ppxf_util.gaussian_filter1d(tpl['FLUX'].data[34], stellar_sigma/velscale)
    # - Add the lines for the gas, all with the same normalization and velocity
    # dispersion
    gas_sigma = 50.
    norm = 20.
    etpl_sigma = []
    for i in range(etpl.ntpl):
        flux += norm*ppxf_util.gaussian_filter1d(etpl.flux[i], gas_sigma/velscale)
    # - Add noise
    ferr = numpy.full(flux.size, 0.05, dtype=float)
    rng = numpy.random.default_rng(seed=99)
    flux += rng.normal(scale=ferr)

    # - Mask the first and last 5-sigma_star pixels
    mask = numpy.zeros(flux.size, dtype=bool)
    nmask = int(5 * stellar_sigma/velscale)
    mask[:nmask] = True
    mask[-nmask:] = True

    return tpl, etpl, stellar_sigma, gas_sigma, velscale, wave, flux, ferr, mask, tpl_sres.copy()


