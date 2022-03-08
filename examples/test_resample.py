#!/usr/bin/env python3

import os
import time
import numpy

from matplotlib import pyplot

from astropy.io import fits

from mangadap.config import manga, defaults
from mangadap.datacube import MaNGADataCube
from mangadap.util.sampling import Resample, grid_borders
from mangadap.proc.bandpassfilter import passband_integral

try:
    import spectres
except:
    spectres = None

#-----------------------------------------------------------------------------

def get_redshift(plt, ifu, drpall_file=None):
    """
    Get the redshift of a galaxy from the DRPall file.

    Args:
        plt (:obj:`int`):
            Plate number
        ifu (:obj:`int`):
            IFU identifier
        drapall_file (:obj:`str`, optional):
            DRPall file. If None, attempts to use the default path to
            the file using environmental variables.
    
    Returns:
        :obj:`float`: The redshift to the galaxy observed by the
        provided PLATEIFU.
    """
    if drpall_file is None:
        drpall_file = manga.drpall_file()
    if not drpall_file.exists():
        raise FileNotFoundError(f'Could not find DRPall file: {drpall_file}')
    hdu = fits.open(str(drpall_file))
    indx = hdu[1].data['PLATEIFU'] == '{0}-{1}'.format(plt, ifu)
    return hdu[1].data['NSA_Z'][indx][0]

#-----------------------------------------------------------------------------

def resample_test():

    plt = 7815
    ifu = 3702
    drpver = 'v3_1_1'
    directory_path = defaults.dap_source_dir() / 'data' / 'remote'
    cube = MaNGADataCube.from_plateifu(plt, ifu, directory_path=directory_path)
    drpall_file = directory_path / f'drpall-{drpver}.fits'
    z = get_redshift(plt, ifu, drpall_file)

    old_wave = cube.wave
    old_flux = numpy.ma.MaskedArray(cube.flux[10,10:12,:], mask=cube.mask[10,10:12,:] > 0)
    old_flux[:,(old_wave > 5570) & (old_wave < 5586)] = numpy.ma.masked
    old_ferr = numpy.ma.power(cube.ivar[10,10:12,:], -0.5)
    indx = (old_wave > old_wave[0]/(1+z)) & (old_wave < old_wave[-2]/(1+z))

    if spectres is not None:
        t = time.perf_counter()
        new_flux_spectres = numpy.empty((old_flux.shape[0], numpy.sum(indx)), dtype=float)
        new_ferr_spectres = numpy.empty((old_flux.shape[0], numpy.sum(indx)), dtype=float)
        for i in range(old_flux.shape[0]):
            new_flux_spectres[i,:], new_ferr_spectres[i,:] \
                    = spectres.spectres(old_wave[indx], old_wave/(1+z), old_flux[i,:].filled(0.0),
                                        spec_errs=old_ferr[i,:].filled(0.0))
        print('SpectRes Time: ', time.perf_counter()-t)

    t = time.perf_counter()
    borders = grid_borders(numpy.array([old_wave[0],old_wave[-1]]), old_wave.size, log=True)[0]
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
    if spectres is not None:
        print('    spectres - brute    = {0:.5e}'.format(
                numpy.mean(numpy.absolute(new_flux_spectres-new_flux_brute[:,indx]))))
        print('    spectres - resample = {0:.5e}'.format(
                numpy.mean(numpy.absolute(new_flux_spectres-r.outy[:,indx]))))
    print('    brute - resample    = {0:.5e}'.format(
            numpy.mean(numpy.absolute(new_flux_brute-r.outy))))

    for i in range(old_flux.shape[0]):
        pyplot.plot(old_wave/(1+z), old_flux[i,:])
        if spectres is not None:
            pyplot.plot(old_wave[indx], new_flux_spectres[i,:])
        pyplot.plot(old_wave, new_flux_brute[i,:])
        pyplot.plot(r.outx, r.outy[i,:])
        pyplot.plot(r.outx, r.outf[i,:])
        pyplot.show()

#-----------------------------------------------------------------------------

if __name__ == '__main__':
    resample_test()



