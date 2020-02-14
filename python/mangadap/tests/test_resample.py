
import pytest
import os

import numpy
from scipy import interpolate
from astropy.io import fits

from mangadap.util.sampling import Resample, _pixel_borders
from mangadap.proc.bandpassfilter import passband_integral
from mangadap.tests.util import data_file

import warnings
warnings.simplefilter("ignore", UserWarning)
warnings.simplefilter("ignore", RuntimeWarning)

def test_resample():
    delt = 100
    x = numpy.arange(20)*0.5 + delt
    y = numpy.linspace(0.1, 2.0, 20)
    y = y*y - y*y*y + 0.2*y + 3 + 0.2 * numpy.square(numpy.sin(50*y))
    newRange = numpy.array([2,8]) + delt
    r = Resample(y, xRange=[x[0],x[-1]], newRange=newRange, newpix=40, newLog=True, step=False)
    assert numpy.isclose(r.outx[0], newRange[0]), 'Incorrect starting value'
    assert numpy.mean(interpolate.interp1d(x,y)(r.outx) - r.outy) < 1e-3, \
                'Bad interpolation difference'
    rc = Resample(y, xRange=[x[0],x[-1]], newRange=newRange, newpix=40, newLog=True,
                  conserve=True, step=False)
    assert numpy.isclose(numpy.mean(r.outy/rc.outy), 3.250907334), \
                'Flux conservation test failed'

def test_against_brute_force():
    specfile = data_file('MaNGA_test_spectra.fits.gz')
    hdu = fits.open(specfile)

    # Just test on the first spectrum
    old_wave = hdu['WAVE'].data
    old_flux = numpy.ma.MaskedArray(hdu['FLUX'].data[0,:], mask=hdu['MASK'].data[0,:] > 0)
    old_flux[(old_wave > 5570) & (old_wave < 5586)] = numpy.ma.masked
    old_ferr = numpy.ma.power(hdu['IVAR'].data[0,:], -0.5)
    z = hdu['Z'].data[0]

    # Do the brute force calculation
    borders = _pixel_borders(numpy.array([old_wave[0],old_wave[-1]]), old_wave.size, log=True)[0]
    _p = numpy.repeat(borders, 2)[1:-1].reshape(-1,2)
    new_flux_brute = passband_integral(old_wave/(1+z), old_flux, passband=_p, log=True)
    new_flux_brute /= (_p[:,1]-_p[:,0])

    # Use resample
    r = Resample(old_flux, e=old_ferr, x=old_wave/(1+z), newRange=[old_wave[0], old_wave[-1]],
                 inLog=True, newLog=True)

    # The two shoule be the same
    assert numpy.isclose(numpy.mean(numpy.absolute(new_flux_brute-r.outy)), 0.0), \
                'Resampling and brute force calculation should be identical'

def test_against_brute_force_twod():
    specfile = data_file('MaNGA_test_spectra.fits.gz')
    hdu = fits.open(specfile)

    # Just test on the first spectrum
    old_wave = hdu['WAVE'].data
    old_flux = numpy.ma.MaskedArray(hdu['FLUX'].data[:2,:], mask=hdu['MASK'].data[:2,:] > 0)
    old_flux[:,(old_wave > 5570) & (old_wave < 5586)] = numpy.ma.masked
    old_ferr = numpy.ma.power(hdu['IVAR'].data[:2,:], -0.5)
    z = hdu['Z'].data[0]    # Only use one redshift so that the input coordinates are the same

    # Do the brute force calculation
    borders = _pixel_borders(numpy.array([old_wave[0],old_wave[-1]]), old_wave.size, log=True)[0]
    _p = numpy.repeat(borders, 2)[1:-1].reshape(-1,2)
    new_flux_brute = numpy.expand_dims(passband_integral(old_wave/(1+z), old_flux[0,:],
                                                         passband=_p, log=True), axis=0)
    new_flux_brute = numpy.append(new_flux_brute, numpy.expand_dims(
                                        passband_integral(old_wave/(1+z), old_flux[1,:],
                                                          passband=_p, log=True), axis=0), axis=0) 
    new_flux_brute /= (_p[:,1]-_p[:,0])[None,:]

    # Use resample
    r = Resample(old_flux, e=old_ferr, x=old_wave/(1+z), newRange=[old_wave[0], old_wave[-1]],
                 inLog=True, newLog=True)

    # The two shoule be the same
    assert numpy.isclose(numpy.mean(numpy.absolute(new_flux_brute-r.outy)), 0.0), \
                'Resampling and brute force calculation should be identical'


