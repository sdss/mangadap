#!/usr/bin/env python3

import numpy
from os import remove
import os.path
from time import clock

#-----------------------------------------------------------------------------

#from mangadap.util.instrument import spectral_resolution, convolution_variable_sigma
from mangadap.util.instrument import match_spectral_resolution
from mangadap.util.constants import DAPConstants
from matplotlib import pyplot

def run_test():

#    # Define some constants
#    cnst = constants()

    #-------------------------------------------------------------------
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

    #-------------------------------------------------------------------

    #-------------------------------------------------------------------
    # Set the target spectral resolution
    new_sres = sres/3.0 + 2.0*sres[0]/3.0 + 300
    new_sres = sres-sres+900
    new_fwhm = wave/new_sres
    new_sigma = new_fwhm/DAPConstants.sig2fwhm

    # Set the flux to a set of uniform emission lines
    expected_flux = numpy.zeros(wave.shape, dtype=numpy.float64)
    gc = numpy.linspace(wave[0], wave[-1], num=20)
    for c in gc:
        expected_flux += numpy.exp(-0.5*numpy.square((wave - c)/new_sigma))/numpy.sqrt(2.0*numpy.pi)/new_sigma

    # Plot the result
    pyplot.plot(wave, flux)
    pyplot.plot(wave, expected_flux, 'r')
    pyplot.show()

    # Compare old and new
    pyplot.plot(wave, sres)
    pyplot.plot(wave, new_sres, 'g')
    pyplot.show()

    pyplot.plot(wave, fwhm)
    pyplot.plot(wave, new_fwhm, 'g')
    pyplot.show()
    #-------------------------------------------------------------------

#    #-------------------------------------------------------------------
#    # Match the resolution
#    res = []
#    res = res + [spectral_resolution(wave, sres)]
#    res = res + [spectral_resolution(wave, new_sres)]
#
#    min_sig = 0.0
#    res[0].match(res[1], min_sig_pix=min_sig)
#
#    # Show the kernel dispersion as a function of wavelength
#    pyplot.plot(wave, res[0].sig_pd)
#    pyplot.show()
#
##    print(res[0].sig_pd[0])
##    print(res[0].sig_mask[0])
##    print(fwhm[0])
##    print(numpy.sqrt( numpy.square(fwhm[0]/cnst.sig2fwhm) + numpy.square(res[0].sig_pd[0]))*cnst.sig2fwhm)
#
#    # Convolve the spectrum
#    new_flux = numpy.copy(flux)
#    indx = numpy.where(res[0].sig_pd > min_sig)
##    print(indx)
#
#    new_flux[indx] = convolution_variable_sigma(flux[indx], res[0].sig_pd[indx])


    # Match the resolution
    new_flux, matched_sres, sigma_offset, new_mask, _ = \
        match_spectral_resolution(wave, flux, sres, wave, new_sres, min_sig_pix=0.0)

    # Overplot the old and new spectrum
    pyplot.plot(wave, flux)
    pyplot.plot(wave, new_flux, 'g')
    pyplot.show()

    pyplot.plot(wave, flux)
    pyplot.plot(wave, new_flux, 'g')
    pyplot.plot(wave, expected_flux, 'r', lw=0.5)
    pyplot.show()

    pyplot.plot(wave, sres)
    pyplot.plot(wave, matched_sres, 'r')
    pyplot.plot(wave, new_sres, 'g')
    pyplot.show()
    #-------------------------------------------------------------------


    

#-----------------------------------------------------------------------------

if __name__ == '__main__':
    t = clock()

    run_test()

    print('Elapsed time: {0} seconds'.format(clock() - t))



