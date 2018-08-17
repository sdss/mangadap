#!/usr/bin/env python3

import os
import time
import numpy

from astropy.io import fits
from matplotlib import pyplot
from matplotlib.ticker import NullFormatter

from mangadap.drpfits import DRPFits
from mangadap.util.constants import DAPConstants
from mangadap.util.fitsutil import DAPFitsUtil
from mangadap.util.instrument import SpectralResolution
from mangadap.util.filter import interpolate_masked_vector
from mangadap.util.fileio import read_template_spectrum
from mangadap.util.pixelmask import SpectralPixelMask

from mangadap.par.artifactdb import ArtifactDB
from mangadap.par.emissionmomentsdb import EmissionMomentsDB
from mangadap.par.emissionlinedb import EmissionLineDB

from mangadap.proc.templatelibrary import available_template_libraries
from mangadap.proc.templatelibrary import TemplateLibraryDef, TemplateLibrary
from mangadap.proc.emissionlinemoments import EmissionLineMoments
from mangadap.proc.sasuke import Sasuke
from mangadap.proc.ppxffit import PPXFFit
from mangadap.proc.stellarcontinuummodel import StellarContinuumModel, StellarContinuumModelBitMask
from mangadap.proc.emissionlinemodel import EmissionLineModelBitMask

#-----------------------------------------------------------------------------

def get_redshift(plt, ifu, drpall_file):
#    hdu = fits.open(os.path.join(os.environ['MANGA_SPECTRO_REDUX'], os.environ['MANGADRP_VER'],
#                                 'drpall-{0}.fits'.format(os.environ['MANGADRP_VER'])))
    hdu = fits.open(drpall_file)
    indx = hdu[1].data['PLATEIFU'] == '{0}-{1}'.format(plt, ifu)
    return hdu[1].data['NSA_Z'][indx][0]


def get_spectrum(plt, ifu, x, y):
    drpf = DRPFits(plt, ifu, 'CUBE', read=True)
    flat_indx = drpf.spatial_shape[1]*x+y
    # This function always returns as masked array
    sres = drpf.spectral_resolution(toarray=True, fill=True, pre=True).data
    flux = drpf.copy_to_masked_array(ext='FLUX', flag=drpf.do_not_fit_flags())
    ivar = drpf.copy_to_masked_array(ext='IVAR', flag=drpf.do_not_fit_flags())
    return drpf['WAVE'].data, flux[flat_indx,:], ivar[flat_indx,:], sres[flat_indx,:]


def create_mastar_template_spectrum():
    output_file = 'mastar-8047-9101-tpl.fits'

    hdu = fits.open('mastar-LOG-8047-9101.fits.gz')
    wave = hdu['MASTAR'].data['WAVE'].ravel()
    flux = hdu['MASTAR'].data['FLUX'].ravel()
    sres = interpolate_masked_vector(numpy.ma.power(DAPConstants.sig2fwhm
                                            * hdu['MASTAR'].data['PREDISP'].ravel() / wave, -1))

    hdr = fits.Header()
    hdr['CRPIX1'] = 1
    hdr['CRVAL1'] = numpy.log10(wave[0])
    hdr['CDELT1'] = numpy.diff(numpy.log10(wave[:2]))[0]
    hdr['CD1_1'] = hdr['CDELT1']
    fits.HDUList([ fits.PrimaryHDU(flux, header=hdr),
                   fits.ImageHDU(sres, name='SPECRES')]).writeto(output_file, overwrite=True)
    
        

#-----------------------------------------------------------------------------

if __name__ == '__main__':
    t = time.clock()

    # Get the list of available libraries
    tpllib_list = available_template_libraries()
    print([tpl['key'] for tpl in tpllib_list])
    print(tpllib_list[0].data)

    # ------------------------------------------------------------------
    # Create a library based on a subset of an existing one

    # Add a new library that only has 10 Gyr SSPs from MIUSCAT
    file_search = os.path.join(os.environ['MANGADAP_DIR'], 'data', 'spectral_templates',
                               'miuscat', '*T10.0000.fits')

    # Define the template library
    tpllib_list += [ TemplateLibraryDef('OLD',
                                        file_search=file_search,
                                        fwhm=2.51,
                                        in_vacuum=False,
                                        wave_limit=numpy.array([3480, 9430]),
                                        log10=False) ]

    # ------------------------------------------------------------------
    # Create a new library based on a MaStar spectrum
    create_mastar_template_spectrum()
#    wave, flux, sres = read_template_spectrum('mastar-8047-9101-tpl.fits', sres_ext='SPECRES')
#    pyplot.plot(wave, flux)
#    pyplot.plot(wave, sres)
#    pyplot.show()

    tpllib_list += [ TemplateLibraryDef('MaStar',
                                        file_search='mastar-8047-9101-tpl.fits',
                                        sres_ext='SPECRES',
                                        in_vacuum=True,
                                        log10=True) ]

    # ------------------------------------------------------------------
    # Let's try to use them!

    # Plate-IFU to use
    plt = 7815
    ifu = 6101
    # Spaxel coordinates
    x = 25
    y = 25

    # Template keywords
#    sc_tpl_key = 'MILESHC'
#    sc_tpl_key = 'OLD'
    sc_tpl_key = 'MaStar'

    # Template pixel scale a factor of 4 smaller than galaxy data
    velscale_ratio = 4

    # Get the redshift
    drpall_file = './data/drpall-v2_4_3.fits'
    z = numpy.array([get_redshift(plt, ifu, drpall_file)])
    print('Redshift: {0}'.format(z[0]))
    dispersion = numpy.array([100.])

    # Read a spectrum
    print('reading spectrum')
    wave, flux, ivar, sres = get_spectrum(plt, ifu, x, y)
    # Fitting functions expect data to be in 2D arrays (for now):
    flux = flux.reshape(1,-1)
    ferr = numpy.ma.power(ivar, -0.5).reshape(1,-1)
    sres = sres.reshape(1,-1)

    #-------------------------------------------------------------------
    # Fit the stellar continuum

    # Mask the 5577 sky line and the emission lines
    sc_pixel_mask = SpectralPixelMask(artdb=ArtifactDB('BADSKY'), emldb=EmissionLineDB('ELPFULL'))

    # Construct the template library
    sc_tpl = TemplateLibrary(sc_tpl_key, tpllib_list=tpllib_list, match_to_drp_resolution=False,
                             velscale_ratio=velscale_ratio, spectral_step=1e-4, log=True,
                             hardcopy=False)
    sc_tpl_sres = numpy.mean(sc_tpl['SPECRES'].data, axis=0).ravel()

    # Instantiate the fitting class
    ppxf = PPXFFit(StellarContinuumModelBitMask())

    # Perform the fit
    cont_wave, cont_flux, cont_mask, cont_par \
        = ppxf.fit(sc_tpl['WAVE'].data.copy(), sc_tpl['FLUX'].data.copy(), wave, flux, ferr,
                   z, dispersion, iteration_mode='no_global_wrej', reject_boxcar=100,
                   ensemble=False, velscale_ratio=velscale_ratio, mask=sc_pixel_mask,
                   matched_resolution=False, tpl_sres=sc_tpl_sres, obj_sres=sres, degree=8,
                   moments=2, plot=False)

    # How many templates did it find:
    print(sc_tpl.ntpl)
    # From what files:
    print(sc_tpl.file_list)
    # What were the weights assigned to each template
    print(cont_par['TPLWGT'][0,:])

    # Remask the continuum fit
    sc_continuum = StellarContinuumModel.reset_continuum_mask_window(
                        numpy.ma.MaskedArray(cont_flux, mask=cont_mask>0))

    # Show the fit and residual
    pyplot.plot(wave, flux[0,:], label='Data')
    pyplot.plot(wave, sc_continuum[0,:], label='Model')
    pyplot.plot(wave, flux[0,:] - sc_continuum[0,:], label='Resid')
    pyplot.legend()
    pyplot.xlabel('Wavelength')
    pyplot.ylabel('Flux')
    pyplot.show()
    #-------------------------------------------------------------------

    print('Elapsed time: {0} seconds'.format(time.clock() - t))

