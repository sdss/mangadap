#!/usr/bin/env python3

import os
import time
import numpy

from matplotlib import pyplot

from astropy.io import fits

from mangadap.util.sampling import Resample, _pixel_borders
from mangadap.proc.bandpassfilter import passband_integral

import spectres

#-----------------------------------------------------------------------------

def resample_test():
    root = os.path.join(os.environ['MANGA_SPECTRO_REDUX'], os.environ['MANGADRP_VER'])
   
    pltifu = '7815-1901'
    hdu = fits.open(os.path.join(root, pltifu.split('-')[0], 'stack',
                                 'manga-{0}-LOGCUBE.fits.gz'.format(pltifu)))

    drpall_file = os.path.join(root, 'drpall-{0}.fits'.format(os.environ['MANGADRP_VER']))
    drpall = fits.open(drpall_file)[1].data

    indx = drpall['PLATEIFU'] == pltifu
    z = drpall['NSA_Z'][indx][0]
    print(z)

    old_wave = hdu['WAVE'].data
    old_flux = numpy.ma.MaskedArray((hdu['FLUX'].data[:,10:12,10]).T,
                                    mask=(hdu['MASK'].data[:,10:12,10]).T > 0)
    old_flux[:,(old_wave > 5570) & (old_wave < 5586)] = numpy.ma.masked
    old_ferr = numpy.ma.power((hdu['IVAR'].data[:,10:12,10]).T, -0.5)
    indx = (old_wave > old_wave[0]/(1+z)) & (old_wave < old_wave[-2]/(1+z))

    t = time.perf_counter()
    new_flux_spectres = numpy.empty((old_flux.shape[0], numpy.sum(indx)), dtype=float)
    new_ferr_spectres = numpy.empty((old_flux.shape[0], numpy.sum(indx)), dtype=float)
    for i in range(old_flux.shape[0]):
        new_flux_spectres[i,:], new_ferr_spectres[i,:] \
                = spectres.spectres(old_wave[indx], old_wave/(1+z), old_flux[i,:].filled(0.0),
                                    spec_errs=old_ferr[i,:].filled(0.0))
    print('SpectRes Time: ', time.perf_counter()-t)

    t = time.perf_counter()
    borders = _pixel_borders(numpy.array([old_wave[0],old_wave[-1]]), old_wave.size, log=True)[0]
    _p = numpy.repeat(borders, 2)[1:-1].reshape(-1,2)
    new_flux_brute = numpy.array([passband_integral(old_wave/(1+z), f, passband=_p, log=True)
                                        for f in old_flux.filled(0.0)])
    new_flux_brute /= (_p[:,1]-_p[:,0])[None,:]
    print('Brute Force Time: ', time.perf_counter()-t)

    t = time.perf_counter()
    r = Resample(old_flux, e=old_ferr, x=old_wave/(1+z), newRange=[old_wave[0], old_wave[-1]],
                 inLog=True, newLog=True)
    print('Resample Time: ', time.perf_counter()-t)

    print('Mean diff:')
    print('    spectres - brute    = {0:.5e}'.format(
            numpy.mean(numpy.absolute(new_flux_spectres-new_flux_brute[:,indx]))))
    print('    spectres - resample = {0:.5e}'.format(
            numpy.mean(numpy.absolute(new_flux_spectres-r.outy[:,indx]))))
    print('    brute - resample    = {0:.5e}'.format(
            numpy.mean(numpy.absolute(new_flux_brute-r.outy))))

    for i in range(old_flux.shape[0]):
        pyplot.plot(old_wave/(1+z), old_flux[i,:])
        pyplot.plot(old_wave[indx], new_flux_spectres[i,:])
        pyplot.plot(old_wave, new_flux_brute[i,:])
        pyplot.plot(r.outx, r.outy[i,:])
        pyplot.plot(r.outx, r.outf[i,:])
        pyplot.show()

#-----------------------------------------------------------------------------

if __name__ == '__main__':
    resample_test()



