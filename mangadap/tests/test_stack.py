
from IPython import embed

import numpy
from scipy import special
from matplotlib import pyplot

import astropy.constants

from mangadap.proc.spectralstack import SpectralStack

def pixelated_gaussian(x, c=0.0, s=1.0):
    """
    x must be linearly sampled

    x, c, and s must be in the same units.
    """
    n = numpy.sqrt(2.)*s
    d = numpy.asarray(x)-c
    dx = numpy.mean(numpy.diff(x))
    return (special.erf((d+dx/2.)/n) - special.erf((d-dx/2.)/n))/2.

def test_register():

    # Generate some line centers and redshifts
    nlines = 10
    nspec = 20
    rng = numpy.random.default_rng()
    line_flux = rng.uniform(low=10, high=100, size=nlines)
    center = rng.uniform(low=3650., high=10000, size=nlines)
    z = rng.uniform(low=0.03, high=0.1, size=20)
    cz = z*astropy.constants.c.to('km/s').value
    sigma = 1 # In pixels

    # Make the wavelength vector
    wave = numpy.logspace(*numpy.log10([3600., 10300.]).tolist(), 3000)
    # Some convenience things to make sure the Gaussian profiles are sampled in
    # pixels
    dlogw = numpy.mean(numpy.diff(numpy.log10(wave)))
    coff = numpy.log10(wave[0])/dlogw
    x = numpy.arange(wave.size)

    # Construct the spectra
    flux = numpy.zeros((nspec,wave.size), dtype=float)
    for i in range(nspec):
        _center = center*(1+z[i])
        _center = numpy.log10(_center)/dlogw - coff
        for f,c in zip(line_flux,_center):
            flux[i] += f * pixelated_gaussian(x, c=c, s=sigma)

    # Stacking object
    stacker = SpectralStack()
    # Stack the spectra including the offset by the input velocity
    swave, sflux, sfdev, snpix, sivar, ssres, scovar \
            = stacker.stack(wave, flux, cz=cz, log=True, muse_mode=False)

    # Just register the spectra to make sure the registration works;
    # this is a simple check on the above
    rwave, rflux, rivar, rsres = SpectralStack.register(wave, cz, flux, log=True)
    reg_stack_flux = numpy.ma.mean(rflux, axis=0)

    assert numpy.allclose(sflux[0], reg_stack_flux), 'Stack and by-hand check failed.'

    # Construct the deredshifted stack from scratch
    rcoff = numpy.log10(rwave[0])/dlogw
    rx = numpy.arange(rwave.size)
    model_flux = numpy.zeros(reg_stack_flux.size, dtype=float)
    rcenter = numpy.log10(center)/dlogw - rcoff
    for f,c in zip(line_flux,rcenter):
        model_flux += f * pixelated_gaussian(rx, c=c, s=sigma)

    # Compare the stacks against truth
#    pyplot.plot(swave, sflux[0])            # Stacked using stacker.stack (blue)
#    pyplot.plot(rwave, reg_stack_flux)      # Stacked by hand after registration (orange)
#    pyplot.plot(rwave, model_flux)          # Truth (green)
#    pyplot.show()

    xcor = numpy.correlate(sflux[0], model_flux)
    assert numpy.argmax(xcor) == xcor.size//2, 'Should be no lag between the model and stack'
    assert numpy.absolute(numpy.sum(sflux[0])/numpy.sum(model_flux) - 1) < 0.01, \
            'Stack sums should be different by less than 1%'


