#!/usr/bin/env python3

import numpy
from os import remove
import os.path
from time import clock

from mangadap.util.instrument import log_rebin
#from mangadap.contrib.ppxf_util import log_rebin
from mangadap.util.constants import constants
from matplotlib import pyplot
from astropy.io import fits

#-----------------------------------------------------------------------------

def run_test():
    # Define some constants
    cnst = constants()

    #-------------------------------------------------------------------
    # Build a fake spectrum
    wave = numpy.arange(3650.0, 10501.0, 0.1)

    # Set the resolution to 2.5 angstroms
    fwhm = numpy.zeros(wave.shape, dtype=numpy.float64)
    fwhm += 2.5
    sigma = fwhm/cnst.sig2fwhm
    sres = wave/fwhm

    # Set the flux to a set of uniform emission lines
    spec = numpy.zeros(wave.shape, dtype=numpy.float64)
    gc = numpy.linspace(wave[0], wave[-1], num=20)
    gs = sigma
    for c in gc:
        spec += numpy.exp(-0.5*numpy.square((wave - c)/gs))/numpy.sqrt(2.0*numpy.pi)/gs

#    hdu = fits.HDUList([fits.PrimaryHDU()])
#    hdu.append(fits.ImageHDU(wave,name='WAVE'))
#    hdu.append(fits.ImageHDU(spec,name='FLUX'))
#    hdu.writeto('fake_spec.fits')

    # Plot the result
    pyplot.plot(wave, spec)
    pyplot.show()
    #-------------------------------------------------------------------

    # Rebin it
    spec_new, wave_new, velscale = log_rebin( [wave[0], wave[-1]], spec, wave_in_ang=True )

    new_range = numpy.array([ 7000., 12000. ])

    spec_10, wave_10, velscale_10 = log_rebin( [wave[0], wave[-1]], spec, newRange=new_range, wave_in_ang=True, unobs=-100.)
#    spec_10, wave_10, velscale_10 = log_rebin( [wave[0], wave[-1]], spec)

    pyplot.plot(wave, spec)
    pyplot.plot(wave_new, spec_new, 'g', linewidth=3)
#    pyplot.plot(numpy.exp(wave_new), spec_new, 'g', linewidth=3)
#    pyplot.plot(numpy.power(10., wave_10), spec_10, 'r')
#    pyplot.plot(numpy.exp(wave_10), spec_10, 'r')
    pyplot.plot(wave_10, spec_10, 'r')
    pyplot.show()


#-----------------------------------------------------------------------------

if __name__ == '__main__':
    t = clock()

    run_test()

    print('Elapsed time: {0} seconds'.format(clock() - t))



