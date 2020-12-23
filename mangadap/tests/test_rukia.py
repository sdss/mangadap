import pytest
import os

from IPython import embed

import numpy
from scipy import special
from matplotlib import pyplot

from astropy.io import fits
import astropy.constants

from mangadap.util.constants import DAPConstants
from mangadap.tests.util import data_test_file
from mangadap.proc.rukia import Rukia, LegendrePolynomial

#import warnings
#warnings.simplefilter('error', UserWarning)

def pixelated_gaussian(x, c=0.0, s=1.0):
    """
    x must be linearly sampled

    x, c, and s must be in the same units.
    """
    n = numpy.sqrt(2.)*s
    d = numpy.asarray(x)-c
    dx = numpy.mean(numpy.diff(x))
    return (special.erf((d+dx/2.)/n) - special.erf((d-dx/2.)/n))/2.


def read_xsl(ifile, minimum_useful=True):
    if ifile is None:
        return None, None, None
    with fits.open(ifile) as hdu:
        wave = numpy.ma.MaskedArray(hdu['SPECTRUM'].data['WAVE']*10,
                                    mask=numpy.logical_not(numpy.isfinite(
                                                hdu['SPECTRUM'].data['WAVE'])))
        flux = numpy.ma.MaskedArray(hdu['SPECTRUM'].data['FLUX']*1e17,
                                    mask=numpy.logical_not(numpy.isfinite(
                                                hdu['SPECTRUM'].data['FLUX'])))
        ferr = numpy.ma.MaskedArray(hdu['SPECTRUM'].data['ERR']*1e17,
                                    mask=numpy.logical_not(numpy.isfinite(
                                                hdu['SPECTRUM'].data['ERR'])))
        if not minimum_useful:
            return wave, flux, ferr
        indx = wave > float(hdu[0].header['WMIN_UVB'])*10
        return wave[indx], flux[indx], ferr[indx]



def read_miles(ifile):
    with fits.open(ifile) as hdu:
        flux = hdu[0].data
        wave = numpy.arange(flux.size)*hdu[0].header['CDELT1'] + hdu[0].header['CRVAL1']
    return wave, flux


def test_rukia_model():

    sigma = LegendrePolynomial(0)
    mult = LegendrePolynomial(0)
    r = Rukia(sigma, mul_model=mult)

    # Build a fake spectrum
    dw = 1.0
    tpl_wave = numpy.arange(4500.0, 4600.0, dw)
    tpl_flux = numpy.zeros(tpl_wave.shape, dtype=float)
    sigma = 2.5/DAPConstants.sig2fwhm
    for c in numpy.linspace(tpl_wave[0], tpl_wave[-1], num=5)[1:-1]:
        tpl_flux += pixelated_gaussian(tpl_wave, c=c, s=sigma)/dw

    wave = numpy.arange(4500.0, 4600.0, dw) + 3.
    par = numpy.ones(r.np, dtype=float)
    par[0] = 0.
    model = r.model(par, wave=wave, tpl_wave=tpl_wave, tpl_flux=tpl_flux)

    sigma = numpy.sqrt(numpy.square(sigma) + 1)
    flux = numpy.zeros(wave.shape, dtype=float)
    for c in numpy.linspace(tpl_wave[0], tpl_wave[-1], num=5)[1:-1]:
        flux += pixelated_gaussian(wave, c=c, s=sigma)/dw

    assert numpy.amax(numpy.absolute(model-flux)) < 1e-4, 'Bad model construction'


def test_rukia_fit():
    miles_wave, miles_flux = read_miles(data_test_file('MILES_res2.50_star_m0505V.fits'))
    miles_flux /= numpy.median(miles_flux)

    xsl_wave, xsl_flux, _ = read_xsl(data_test_file('xsl_spectrum_X0360_uvb.fits'))
    xsl_flux /= numpy.ma.median(xsl_flux)

    # Setup the fitting functions
    sigma = LegendrePolynomial(1)
    mulp = LegendrePolynomial(5)
    r = Rukia(sigma, mul_model=mulp)

    # Perform the fit
    r.fit(miles_wave, miles_flux, xsl_wave.data, xsl_flux.data, shift=0.0, fit_shift=True,
          rejiter=-1)

    assert numpy.sum(r.gpm) == 2176, 'Change in the number of pixels rejected.'
    assert numpy.sqrt(numpy.mean(numpy.square((r.flux-r.model(r.par))[r.gpm]))) < 0.01, \
            'Fit quality changed.'


