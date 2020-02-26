import pytest

import numpy

from mangadap.util.resolution import match_spectral_resolution
from mangadap.util.constants import DAPConstants

import warnings
warnings.simplefilter("ignore", UserWarning)
warnings.simplefilter("ignore", RuntimeWarning)

def test_match_spec_res():

    # Build a fake spectrum
    wave = numpy.arange(3650.0, 10501.0, 0.5)

    # Set the resolution to 2.5 angstroms
    fwhm = numpy.zeros(wave.shape, dtype=numpy.float64)
    fwhm += 2.5
    sigma = fwhm/DAPConstants.sig2fwhm
    sres = wave/fwhm

    # Set the flux to a set of uniform emission lines
    flux = numpy.zeros(wave.shape, dtype=numpy.float64)
    gc = numpy.linspace(wave[0], wave[-1], num=20)
    gs = sigma
    for c in gc:
        flux += numpy.exp(-0.5*numpy.square((wave - c)/gs))/numpy.sqrt(2.0*numpy.pi)/gs

    # Set the target spectral resolution
    new_sres = sres/3.0 + 2.0*sres[0]/3.0 + 300
    new_sres = sres-sres+900
    new_fwhm = wave/new_sres
    new_sigma = new_fwhm/DAPConstants.sig2fwhm

    # Set the flux to a set of uniform emission lines
    expected_flux = numpy.zeros(wave.shape, dtype=numpy.float64)
    gc = numpy.linspace(wave[0], wave[-1], num=20)
    for c in gc:
        expected_flux += numpy.exp(-0.5*numpy.square((wave - c)/new_sigma)) \
                                / numpy.sqrt(2.0*numpy.pi)/new_sigma

    # Match the resolution
    new_flux, matched_sres, sigma_offset, new_mask, _ = \
        match_spectral_resolution(wave, flux, sres, wave, new_sres, min_sig_pix=0.0)

    assert numpy.mean(new_flux[50:-50]-expected_flux[50:-50]) < 1e-9, 'Failed resolution match!'


