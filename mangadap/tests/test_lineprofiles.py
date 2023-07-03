
from IPython import embed

import numpy
from scipy import special
from matplotlib import pyplot

import astropy.constants

from mangadap.util.sampling import spectrum_velocity_scale, spectral_coordinate_step
from mangadap.util.lineprofiles import FFTGaussianLSF
from mangadap.proc import ppxffit


def pixelated_gaussian(x, c=0.0, s=1.0):
    """
    Return a Gaussian integrated over the pixel width.  This is a PDF such that
    the integral is unity.

    Args:
        x (`numpy.ndarray`_):
            X coordinates for the calculation.  This should be linearly sampled.
        c (:obj:`float`, optional):
            The Gaussian center in the same units as the x coordinates.
        s (:obj:`float`, optional):
            The Gaussian dispersion in the same units as the x coordinates.
    """
    n = numpy.sqrt(2.)*s
    d = numpy.asarray(x)-c
    dx = numpy.mean(numpy.diff(x))
    return (special.erf((d+dx/2.)/n) - special.erf((d-dx/2.)/n))/2./dx


def gaussian(x, c=0.0, s=1.0):
    """
    Return a Gaussian sampled at the pixel center.  For a well-sampled profile,
    the integral of the profile is unity.

    Args:
        x (`numpy.ndarray`_):
            X coordinates for the calculation.
        c (:obj:`float`, optional):
            The Gaussian center in the same units as the x coordinates.
        s (:obj:`float`, optional):
            The Gaussian dispersion in the same units as the x coordinates.
    """
    return numpy.exp(-(x - c)**2 / s**2 / 2) / numpy.sqrt(2 * numpy.pi) / s


def test_line_profile():

    # Create an H-alpha line over the MaNGA wavelength range
    restwave = 6564.608
    wave0 = 3621.5959848601933
    wave = 10**(numpy.log10(wave0) + 1e-4*numpy.arange(4563))
    velscale = spectrum_velocity_scale(wave)
    sigma_inst = 2*velscale
    line_flux = 1.

    # Constants
    base = 10.
    c = astropy.constants.c.to('km/s').value

    # Convert from wavelengths to pixel coordinates
    _wave = numpy.log(wave)/numpy.log(base)
    _dw = spectral_coordinate_step(wave, log=True, base=base)
    _restwave = numpy.log(restwave)/numpy.log(base)
    _restwave_pix = (_restwave - _wave[0])/_dw

    # Flux in pixel units
    dl = restwave*(numpy.power(base,_dw/2)-numpy.power(base,-_dw/2))
    _flux = line_flux / dl
    # Dispersion in pixel units
    _sigma = sigma_inst * restwave / c / dl
    # Amplitude
    _amp = _flux / _sigma / numpy.sqrt(2 * numpy.pi)

    # Construct the templates
    pix = numpy.arange(wave.size)

    # --------------------
    # At 0 velocity
    # - FFTGaussianLSF *requires* pixel units
    fft_flux = FFTGaussianLSF()(pix, numpy.array([_flux, _restwave_pix, _sigma]))
    # - And it's best if the pixelated Gaussian is also in pixel units.
    pix_flux = pixelated_gaussian(pix, c=_restwave_pix, s=_sigma) \
                        * _amp * numpy.sqrt(2 * numpy.pi) * _sigma
    # - But it doesn't matter for the direct samples from the Gaussian profile.
    # It can be in pixels ...
    gau_flux = gaussian(pix, c=_restwave_pix, s=_sigma) * _amp * numpy.sqrt(2 * numpy.pi) * _sigma
    #   ... or angstroms
    sigma_ang = sigma_inst * restwave / c
    amp = line_flux / sigma_ang / numpy.sqrt(2 * numpy.pi)
    _gau_flux = gaussian(wave, c=restwave, s=sigma_ang) * amp * numpy.sqrt(2 * numpy.pi) * sigma_ang

    # All the plots should overly each other, except the *gau* values will have
    # slightly higher peaks.
#    pyplot.plot(wave, fft_flux)
#    pyplot.plot(wave, gau_flux)
#    pyplot.plot(wave, pix_flux)
#    pyplot.plot(wave, _gau_flux)
#    pyplot.show()

    # The peak of the direct Gaussian should be larger than the peak after pixel integration
    assert numpy.amax(gau_flux) > numpy.amax(pix_flux), \
            'Pixel integration should lower the peak value'

    # The pixel-integrated profile and the FFT profile should be nearly the same
    # WARNING: This depends on the sigma value!
    assert numpy.allclose(pix_flux, fft_flux), \
            'Pixel-integrated and analytic FFT profiles are too different'
    # --------------------


    # --------------------
    # At z=0.5
    z = 0.5
    cz_vel = c * z

    # Construct the model as done when fitting the spectra
    #   - First make a template.  The template dispersion is 1/sqrt(2) of the
    #     total sigma (just for testing purposes).
    ppxf_vel = ppxffit.PPXFFit.revert_velocity(cz_vel, 0.)[0]
    _tpl_sigma = _sigma / numpy.sqrt(2)
    tpl_flux = FFTGaussianLSF()(pix, numpy.array([_flux, _restwave_pix, _tpl_sigma]))
    wgt = [1.0]

    #   - Convolve with a Gaussian kernel that also has a dispersion that is
    #     1/sqrt(2) of the total sigma.  The quadrature some of the result is
    #     the same dispersion as in the 0 velocity case.
    # NOTE: tpl_flux is used as the bogus galaxy spectrum.  It is unimportant to
    # the construction of the model.
    model = ppxffit.PPXFModel(tpl_flux.reshape(1,-1).T, tpl_flux, velscale, degree=-1)
    shifted_model_flux = model([ppxf_vel, sigma_inst / numpy.sqrt(2)], wgt) / (1+z)

    # Create the same spectrum by putting the line at the redshifted wavelength,
    # but with the expected dispersion (same as in the 0 velocity case).   Just
    # as above, the FFTGaussianLSF and pixelated_gaussian computations requires
    # pixel coordinates...
    obs_wave = restwave * (1 + cz_vel / c)
    _obs_wave = numpy.log(obs_wave)/numpy.log(base)
    _obs_wave_pix = (_obs_wave - _wave[0])/_dw
    _dl = obs_wave*(numpy.power(base,_dw/2)-numpy.power(base,-_dw/2))
    _flux = line_flux / _dl
    shifted_fft_flux = FFTGaussianLSF()(pix, numpy.array([_flux, _obs_wave_pix, _sigma]))
    _amp = _flux / _sigma / numpy.sqrt(2 * numpy.pi)
    shifted_pix_flux = pixelated_gaussian(pix, c=_obs_wave_pix, s=_sigma) \
                        * _amp * numpy.sqrt(2 * numpy.pi) * _sigma
    # ... but the direct samples can be done using pixels ...
    shifted_gau_flux = gaussian(pix, c=_obs_wave_pix, s=_sigma) * _amp * numpy.sqrt(2 * numpy.pi) * _sigma
    # ... or angstroms.
    sigma_ang = sigma_inst * obs_wave / c
    amp = line_flux / sigma_ang / numpy.sqrt(2 * numpy.pi)
    _shifted_gau_flux = gaussian(wave, c=obs_wave, s=sigma_ang) * amp * numpy.sqrt(2 * numpy.pi) * sigma_ang

    # All the plots should overly each other, except the *gau* values will have
    # slightly higher peaks.
#    pyplot.plot(wave, shifted_model_flux)
#    pyplot.plot(wave, shifted_fft_flux)
#    pyplot.plot(wave, shifted_gau_flux)
#    pyplot.plot(wave, shifted_pix_flux)
#    pyplot.plot(wave, _shifted_gau_flux)
#    pyplot.show()

    # The peak of the direct Gaussian should be larger than the peak after pixel integration
    assert numpy.amax(shifted_gau_flux) > numpy.amax(shifted_pix_flux), \
            'Pixel integration should lower the peak value'

    # The model profile, pixel-integrated profile, and FFT profile should be
    # nearly the same
    # WARNING: This depends on the sigma value!
    assert numpy.allclose(shifted_pix_flux, shifted_fft_flux), \
            'Pixel-integrated and analytic FFT profiles are too different'
    assert numpy.allclose(shifted_model_flux, shifted_fft_flux), \
            'Model calculation and direct analytic FFT profiles are too different'


