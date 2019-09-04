
import os
import glob

import numpy
from matplotlib import pyplot

from scipy import ndimage

from mangadap.util.fileio import read_template_spectrum

from astropy.io import fits

tpl_files = glob.glob('MaStarHC*.fits')
for tpl in tpl_files:
    wave, flux, ivar, sres = read_template_spectrum(tpl, ivar_ext='IVAR', sres_ext='SRES',
                                                    log10=True)
    pyplot.plot(wave, flux)
    snr = flux*numpy.sqrt(ivar)
    pyplot.plot(wave, snr)
    snr_sm = ndimage.uniform_filter(snr, size=50)
    pyplot.plot(wave, snr_sm)
    pyplot.yscale('log')
    pyplot.title('{0}: {1:.2f} {2:.2f}'.format(tpl, numpy.median(snr), numpy.amin(snr_sm)))
    pyplot.show()

