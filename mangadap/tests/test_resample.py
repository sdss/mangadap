
import pytest
import os

from IPython import embed

import numpy
from scipy import interpolate
from astropy.io import fits

from mangadap.util import sampling
from mangadap.proc.bandpassfilter import passband_integral
from mangadap.tests.util import data_test_file

import warnings
warnings.simplefilter("ignore", UserWarning)
warnings.simplefilter("ignore", RuntimeWarning)

def test_resample():
    delt = 100
    x = numpy.arange(20)*0.5 + delt
    y = numpy.linspace(0.1, 2.0, 20)
    y = y*y - y*y*y + 0.2*y + 3 + 0.2 * numpy.square(numpy.sin(50*y))
    newRange = numpy.array([2,8]) + delt
    r = sampling.Resample(y, xRange=[x[0],x[-1]], newRange=newRange, newpix=40, newLog=True,
                          step=False)
    assert numpy.isclose(r.outx[0], newRange[0]), 'Incorrect starting value'
    assert numpy.mean(interpolate.interp1d(x,y)(r.outx) - r.outy) < 1e-3, \
                'Bad interpolation difference'

    _r = sampling.Resample(y, xRange=[x[0],x[-1]], newBorders=r.outborders, step=False)
    assert numpy.array_equal(r.outx, _r.outx), 'Coordinates differ'
    assert numpy.array_equal(r.outborders, _r.outborders), 'Coordinates borders differ'
    assert numpy.array_equal(r.outy, _r.outy), 'Resampling differs'

    _r = sampling.Resample(y, xRange=[x[0],x[-1]], newx=r.outx, newLog=True, step=False)
    assert numpy.array_equal(r.outx, _r.outx), 'Coordinates differ'
    assert numpy.allclose(r.outborders, _r.outborders, rtol=0., atol=1e-12), \
                'Coordinates borders differ'
    assert numpy.allclose(r.outy, _r.outy, rtol=0., atol=1e-12), \
                'Coordinates borders differ'
                
    rc = sampling.Resample(y, xRange=[x[0],x[-1]], newRange=newRange, newpix=40, newLog=True,
                           conserve=True, step=False)
    assert numpy.isclose(numpy.mean(r.outy/rc.outy), 3.250907334), \
                'Flux conservation test failed'


def test_resample_covar():
    delt = 100
    x = numpy.arange(20)*0.5 + delt
    y = numpy.linspace(0.1, 2.0, 20)
    y = y*y - y*y*y + 0.2*y + 3 + 0.2 * numpy.square(numpy.sin(50*y))
    newRange = numpy.array([2,8]) + delt
    e = numpy.full_like(y, 2.)
    r = sampling.Resample(y, e=e, xRange=[x[0],x[-1]], newRange=newRange, newpix=40, newLog=True)

    # Test covar without errors
    _r = sampling.Resample(y, xRange=[x[0],x[-1]], newRange=newRange, newpix=40, newLog=True,
                           covar=True)
    assert numpy.array_equal(r.outy, _r.outy), 'Bad resampling with covariance'

    # Test multiple vectors
    _r = sampling.Resample(numpy.row_stack((y,y)), xRange=[x[0],x[-1]], newRange=newRange,
                           newpix=40, newLog=True, covar=True)
    assert numpy.allclose(r.outy, _r.outy[0], rtol=0.0, atol=1e-15), \
            'Bad resampling with covariance'

    # Test covar with errors
    _r = sampling.Resample(y, xRange=[x[0],x[-1]], e=e, newRange=newRange, newpix=40, newLog=True,
                           covar=True)
    assert numpy.array_equal(r.outy, _r.outy), 'Bad resampling with covariance'
    assert numpy.array_equal(r.oute, _r.oute), 'Bad resampling errors with covariance'

    # Test multiple vectors with errors
    _r = sampling.Resample(numpy.row_stack((y,y)), e=numpy.row_stack((e,e)), xRange=[x[0],x[-1]],
                           newRange=newRange, newpix=40, newLog=True, covar=True)
    assert numpy.allclose(r.outy, _r.outy[0], rtol=0.0, atol=1e-15), \
            'Bad resampling with covariance'
    assert numpy.allclose(r.oute, _r.oute[0], rtol=0.0, atol=1e-15), \
            'Bad resampling errors with covariance'

    # Test with conserve
    r = sampling.Resample(y, e=e, xRange=[x[0],x[-1]], newRange=newRange, newpix=40, newLog=True,
                          conserve=True)
    _r = sampling.Resample(y, xRange=[x[0],x[-1]], e=e, newRange=newRange, newpix=40, newLog=True,
                           covar=True, conserve=True)
    assert numpy.allclose(r.outy, _r.outy, rtol=0.0, atol=1e-15), \
            'Bad resampling with covariance'
    assert numpy.allclose(r.oute, _r.oute, rtol=0.0, atol=1e-15), \
            'Bad resampling errors with covariance'

    _r = sampling.Resample(numpy.row_stack((y,y)), e=numpy.row_stack((e,e)), xRange=[x[0],x[-1]],
                           newRange=newRange, newpix=40, newLog=True, covar=True, conserve=True)
    assert numpy.allclose(r.outy, _r.outy[0], rtol=0.0, atol=1e-15), \
            'Bad resampling with covariance'
    assert numpy.allclose(r.oute, _r.oute[0], rtol=0.0, atol=1e-15), \
            'Bad resampling errors with covariance'


def test_against_brute_force():
    specfile = data_test_file('MaNGA_test_spectra.fits.gz')
    hdu = fits.open(specfile)

    # Just test on the first spectrum
    old_wave = hdu['WAVE'].data
    old_flux = numpy.ma.MaskedArray(hdu['FLUX'].data[0,:], mask=hdu['MASK'].data[0,:] > 0)
    old_flux[(old_wave > 5570) & (old_wave < 5586)] = numpy.ma.masked
    old_ferr = numpy.ma.power(hdu['IVAR'].data[0,:], -0.5)
    z = hdu['Z'].data[0]

    # Do the brute force calculation
    borders = sampling.grid_borders(numpy.array([old_wave[0],old_wave[-1]]), old_wave.size,
                                    log=True)[0]
    _p = numpy.repeat(borders, 2)[1:-1].reshape(-1,2)
    new_flux_brute = passband_integral(old_wave/(1+z), old_flux, passband=_p, log=True)
    new_flux_brute /= (_p[:,1]-_p[:,0])

    # Use resample
    r = sampling.Resample(old_flux, e=old_ferr, x=old_wave/(1+z),
                          newRange=[old_wave[0], old_wave[-1]], inLog=True, newLog=True)

    # The two shoule be the same
    assert numpy.isclose(numpy.mean(numpy.absolute(new_flux_brute-r.outy)), 0.0), \
                'Resampling and brute force calculation should be identical'

def test_against_brute_force_twod():
    specfile = data_test_file('MaNGA_test_spectra.fits.gz')
    hdu = fits.open(specfile)

    # Just test on the first spectrum
    old_wave = hdu['WAVE'].data
    old_flux = numpy.ma.MaskedArray(hdu['FLUX'].data[:2,:], mask=hdu['MASK'].data[:2,:] > 0)
    old_flux[:,(old_wave > 5570) & (old_wave < 5586)] = numpy.ma.masked
    old_ferr = numpy.ma.power(hdu['IVAR'].data[:2,:], -0.5)
    z = hdu['Z'].data[0]    # Only use one redshift so that the input coordinates are the same

    # Do the brute force calculation
    borders = sampling.grid_borders(numpy.array([old_wave[0],old_wave[-1]]), old_wave.size,
                                    log=True)[0]
    _p = numpy.repeat(borders, 2)[1:-1].reshape(-1,2)
    new_flux_brute = numpy.expand_dims(passband_integral(old_wave/(1+z), old_flux[0,:],
                                                         passband=_p, log=True), axis=0)
    new_flux_brute = numpy.append(new_flux_brute, numpy.expand_dims(
                                        passband_integral(old_wave/(1+z), old_flux[1,:],
                                                          passband=_p, log=True), axis=0), axis=0) 
    new_flux_brute /= (_p[:,1]-_p[:,0])[None,:]

    # Use resample
    r = sampling.Resample(old_flux, e=old_ferr, x=old_wave/(1+z),
                          newRange=[old_wave[0], old_wave[-1]], inLog=True, newLog=True)

    # The two shoule be the same
    assert numpy.isclose(numpy.mean(numpy.absolute(new_flux_brute-r.outy)), 0.0), \
                'Resampling and brute force calculation should be identical'


def test_grid_npix():
    for s in numpy.random.randint(2, high=1000, size=100):
        x = numpy.linspace(-1.,1.,s)
        nx, rng = sampling.grid_npix(rng=[x[0],x[-1]], dx=x[1]-x[0])
        assert nx == s, 'Bad number of linear grid points'
        assert numpy.allclose(rng, [-1., 1.]), 'Bad linear grid range'

        x = numpy.logspace(-1.,1.,s)
        nx, rng = sampling.grid_npix(rng=[x[0],x[-1]], dx=numpy.log10(x[1])-numpy.log10(x[0]),
                                     log=True)
        assert nx == s, 'Bad number of logarithmic grid points'
        assert numpy.allclose(rng, numpy.power(10., [-1., 1.])), 'Bad logarithmic grid range'


def test_grid_borders():
    x = numpy.linspace(-1,1,100)
    dx = x[1]-x[0]
    nx, rng = sampling.grid_npix(rng=[x[0],x[-1]], dx=dx)
    borders, _dx = sampling.grid_borders(rng, nx)
    assert dx == _dx, 'Bad linear step size'
    assert numpy.allclose(borders, numpy.linspace(-1-dx/2, 1+dx/2, 101)), 'Bad linear borders'


    x = numpy.logspace(-1,1,100)
    dx = numpy.log10(x[1])-numpy.log10(x[0])
    nx, rng = sampling.grid_npix(rng=[x[0],x[-1]], dx=dx, log=True)
    borders, _dx = sampling.grid_borders(rng, nx, log=True)
    assert dx == _dx, 'Bad logarithmic step size'
    assert numpy.allclose(borders, numpy.logspace(-1-dx/2, 1+dx/2, 101)), 'Bad geometric borders'


def test_grid_centers():
    x = numpy.linspace(-1,1,100)
    dx = x[1]-x[0]
    nx, rng = sampling.grid_npix(rng=[x[0],x[-1]], dx=dx)
    centers, _dx = sampling.grid_centers(rng, nx)
    assert dx == _dx, 'Bad linear step size'
    assert numpy.allclose(centers, x), 'Bad linear borders'

    x = numpy.logspace(-1,1,100)
    dx = numpy.log10(x[1])-numpy.log10(x[0])
    nx, rng = sampling.grid_npix(rng=[x[0],x[-1]], dx=dx, log=True)
    centers, _dx = sampling.grid_centers(rng, nx, log=True)
    assert dx == _dx, 'Bad logarithmic step size'
    assert numpy.allclose(centers, x), 'Bad geometric borders'

if __name__ == '__main__':
    test_resample_covar()


