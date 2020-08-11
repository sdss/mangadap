#!/usr/bin/env python3

import os
import time
import numpy

#-----------------------------------------------------------------------------

from argparse import ArgumentParser

import numpy
from astropy.io import fits
from matplotlib import pyplot
from matplotlib.ticker import NullFormatter

from mangadap.drpfits import DRPFitsBitMask
from mangadap.dapfits import DAPCubeBitMask
from mangadap.util.fileio import rec_to_fits_type
from mangadap.util.filter import BoxcarFilter
from mangadap.util.fitsutil import DAPFitsUtil

from mangadap.par.artifactdb import ArtifactDB
from mangadap.par.emissionlinedb import EmissionLineDB
from mangadap.proc.templatelibrary import TemplateLibrary
from mangadap.util.pixelmask import SpectralPixelMask

from mangadap.proc.ppxffit import PPXFFit
from mangadap.proc.stellarcontinuummodel import StellarContinuumModelBitMask

#-----------------------------------------------------------------------------
def object_data(f, dispaxis, fit_flag_file):

    drpbm = DRPFitsBitMask()

    hdu = fits.open(f)

    if dispaxis == 1:
        flux = hdu['FLUX'].data.copy()
        ivar = hdu['IVAR'].data.copy()
        mask = hdu['MASK'].data.copy()
        specres = hdu['SPECRES'].data.copy()
    else:
        flux = hdu['FLUX'].data.copy().T
        ivar = hdu['IVAR'].data.copy().T
        mask = hdu['MASK'].data.copy().T
        specres = hdu['SPECRES'].data.copy().T
    wave = hdu['WAVE'].data.copy()
    redshift = hdu['Z'].data.copy()
#    usetpl = hdu['USETPL'].data.copy().astype(bool)
    hdu.close()
    del hdu

    nspec = flux.shape[0]

    fit_flags = numpy.ones(nspec, dtype=bool) if fit_flag_file is None \
                    else (numpy.genfromtxt(fit_flag_file)[:,3]).astype(bool)
    if len(fit_flags) != nspec:
        print(nspec)
        print(len(fit_flags))
        raise ValueError('Incorrect number of fitting flags.')

    ferr = numpy.ma.power(ivar, -0.5)
    bool_mask = drpbm.flagged(mask, flag=['DONOTUSE', 'FORESTAR']) | numpy.ma.getmaskarray(ferr)
    
    return wave, numpy.ma.MaskedArray(flux, mask=bool_mask), \
                numpy.ma.MaskedArray(ferr, mask=bool_mask), \
                numpy.ma.MaskedArray(specres, mask=bool_mask), redshift, fit_flags

def init_ax(fig, pos):
    ax = fig.add_axes(pos) #, facecolor='0.95')
    ax.minorticks_on()
    ax.tick_params(which='major', length=10)
    ax.tick_params(which='minor', length=5)
    ax.tick_params(which='both', direction='out', top='on', right='on')
    return ax

#-----------------------------------------------------------------------------
if __name__ == '__main__':
    t = time.clock()

    parser = ArgumentParser()

    parser.add_argument('inp', type=str, help='input fits file with spectra to fit')
    parser.add_argument('out', type=str, help='output fits file with fitted models')

    parser.add_argument('--output_root', type=str, help='root directory for output', default=None)
    parser.add_argument('--tpl_library', type=str, help='template library to use',
                        default='MILESHC')
    parser.add_argument('--dispaxis', type=int, help='dispersion axis', default=1)
    parser.add_argument('--ppxf_iteration', type=str, help='pPXF iteration mode',
                        default='none')
    parser.add_argument('--filter', type=int, default=0,
                        help='Number of filtering iterations (only valid if the ppxf_iteration '
                             'mode is \'fit_reject_filter\'')
    parser.add_argument('--reject_boxcar', type=int, default=100,
                        help='Boxcar for rejection iterations')
    parser.add_argument('--filter_boxcar', type=int, default=100,
                        help='Boxcar for filtering iterations')
    parser.add_argument('--filter_opt', type=str, default='divide',
                        help='Operation for the smoothed spectra: \'divide\' or \'subtract\'')

    parser.add_argument('--degree', type=int, default=8, help='additive order')
    parser.add_argument('--mdegree', type=int, default=0, help='multiplicative order')
    parser.add_argument('--filt_degree', type=int, default=8, help='filtered additive order')
    parser.add_argument('--filt_mdegree', type=int, default=0, help='filtered multiplicative order')

    parser.add_argument('--spec_flags', type=str, default=None,
                        help='fourth column sets whether (1) or not (0) to fit spectrum')

    arg = parser.parse_args()

    if not os.path.isfile(arg.inp):
        raise FileNotFoundError('No file: {0}'.format(arg.inp))

    pixel_mask = SpectralPixelMask(artdb=ArtifactDB('BADSKY'), emldb=EmissionLineDB('ELPFULL'))
    directory_path = '.' if arg.output_root is None else arg.output_root

    # Open and check the file
    wave, flux, ferr, specres, redshift, fit_spectrum \
            = object_data(arg.inp, arg.dispaxis, arg.spec_flags)

    nspec,npix = flux.shape
    print('Read: {0}'.format(arg.inp))
    print('Contains {0} spectra'.format(nspec))
    print(' each with {0} pixels'.format(npix))
    print('Fitting {0} spectra.'.format(numpy.sum(fit_spectrum)))

#    spec = numpy.arange(flux.shape[0])[fit_spectrum]
#    pyplot.plot(wave, flux[spec[0],:])
#    pyplot.show()

    tpl = TemplateLibrary(arg.tpl_library,
                          match_to_drp_resolution=False,
                          velscale_ratio=4,     # Keep the pixel size at 1/4 the MaNGA step
                          spectral_step=1e-4,
                          log=True,
                          hardcopy=False)
    tpl_sres = numpy.mean(tpl['SPECRES'].data, axis=0).ravel()

#    pyplot.plot(tpl['WAVE'].data, tpl['FLUX'].data[0,:])
#    pyplot.show()

    ppxf = PPXFFit(StellarContinuumModelBitMask())

    guess_dispersion = numpy.full(nspec, 100., dtype=numpy.float)

#    obj_to_fit = numpy.where( numpy.any(numpy.invert(flux.mask), axis=1) & fit_spectrum )[0][:1]
#    print(len(obj_to_fit))
#    flux = flux[obj_to_fit,:]
#    ferr = ferr[obj_to_fit,:]
#    redshift = redshift[obj_to_fit]
#    guess_dispersion = guess_dispersion[obj_to_fit]
#    specres = specres[obj_to_fit,:]

    indx = numpy.any(numpy.invert(flux.mask), axis=1) & fit_spectrum
    flux[~indx,:] = numpy.ma.masked
#    i = numpy.arange(flux.shape[0])[indx][1:]
#    flux[i,:] = numpy.ma.masked

    model_wave, model_flux, model_mask, model_par \
        = ppxf.fit(tpl['WAVE'].data.copy(), tpl['FLUX'].data.copy(), wave, flux, ferr, redshift,
                   guess_dispersion, iteration_mode=arg.ppxf_iteration,
                   reject_boxcar=arg.reject_boxcar, filter_boxcar=arg.filter_boxcar,
                   filter_operation=arg.filter_opt, filter_iterations=arg.filter, ensemble=False,
                   velscale_ratio=4, mask=pixel_mask, matched_resolution=False,
                   tpl_sres=tpl_sres, obj_sres=specres, degree=arg.degree, mdegree=arg.mdegree,
                   filt_degree=arg.filt_degree, filt_mdegree=arg.filt_mdegree, moments=2)
                   #, plot=True)
    
    fits.HDUList([ fits.PrimaryHDU(),
                   fits.ImageHDU(data=model_wave, name='WAVE'),
                   fits.ImageHDU(data=model_flux.data, name='FLUX'),
                   fits.ImageHDU(data=model_mask, name='MASK'),
                   fits.BinTableHDU.from_columns([ fits.Column(name=n,
                                                        format=rec_to_fits_type(model_par[n]),
                                                               array=model_par[n])
                                                      for n in model_par.dtype.names ], name='PAR')
                 ]).writeto(os.path.join(directory_path, arg.out), overwrite=True)
#    DAPFitsUtil.write(hdu, os.path.join(directory_path, arg.out), clobber=True)


    print('Elapsed time: {0} seconds'.format(time.clock() - t))



