import os
import time

from IPython import embed

import numpy

from astropy.io import fits
from matplotlib import pyplot

from mangadap.datacube import MaNGADataCube
from mangadap.util.fitsutil import DAPFitsUtil
from mangadap.util.resolution import SpectralResolution
from mangadap.util.pixelmask import SpectralPixelMask

from mangadap.par.artifactdb import ArtifactDB
from mangadap.par.emissionmomentsdb import EmissionMomentsDB
from mangadap.par.emissionlinedb import EmissionLineDB

from mangadap.proc.templatelibrary import TemplateLibrary
from mangadap.proc.spectralstack import SpectralStack
from mangadap.proc.emissionlinemoments import EmissionLineMoments
from mangadap.proc.sasuke import Sasuke
from mangadap.proc.ppxffit import PPXFFit
from mangadap.proc.stellarcontinuummodel import StellarContinuumModel, StellarContinuumModelBitMask
from mangadap.proc.emissionlinemodel import EmissionLineModelBitMask
from mangadap.proc.spectralfitting import EmissionLineFit

#-----------------------------------------------------------------------------

def get_redshift(plt, ifu, drpall_file=None):
    hdu = fits.open(os.path.join(os.environ['MANGA_SPECTRO_REDUX'], os.environ['MANGADRP_VER'],
                                 'drpall-{0}.fits'.format(os.environ['MANGADRP_VER']))) \
                if drpall_file is None else fits.open(drpall_file)
    indx = hdu[1].data['PLATEIFU'] == '{0}-{1}'.format(plt, ifu)
    return hdu[1].data['NSA_Z'][indx][0]


def get_spectra(plt, ifu, x, y, directory_path=None):
    cube = MaNGADataCube.from_plateifu(plt, ifu, directory_path=directory_path)
    flat_indx = cube.spatial_shape[1]*x+y
    # This function always returns as masked array
    flux = cube.copy_to_masked_array(attr='flux', flag=cube.do_not_fit_flags())
    ivar = cube.copy_to_masked_array(attr='ivar', flag=cube.do_not_fit_flags())
    sres = cube.copy_to_array(attr='sres')
    return cube.wave, flux[flat_indx,:], ivar[flat_indx,:], sres[flat_indx,:]
        

#-----------------------------------------------------------------------------

if __name__ == '__main__':
    t = time.perf_counter()


    # Plate-IFU to use
    plt = 7815
    ifu = 3702
    # Spaxel coordinates
    x = numpy.array([21,20,19])
    y = numpy.array([21,21,21])
    binid=numpy.array([0,0,1])
    # Number of spectra
    nspec = len(x)

    # Show the ppxf plots
#    fit_plots = True
    fit_plots = False

    # Show summary plots
    usr_plots = True #False

    # Template keywords
    sc_tpl_key = 'MILESHC'
    el_tpl_key = 'MILESHC'
#    el_tpl_key = 'BC03'

    # Emission-line database keywords
    elmom_key = 'ELBMILES'
    elfit_key = 'ELPMILES'
#    elmom_key = 'ELBTT'
#    elfit_key = 'ELPT4'

    # Template pixel scale a factor of 4 smaller than galaxy data
    velscale_ratio = 4

    # Get the redshift
    drpver = 'v2_7_1'
    directory_path = os.path.join(os.environ['MANGADAP_DIR'], 'mangadap', 'data', 'remote')
    drpall_file = os.path.join(directory_path, 'drpall-{0}.fits'.format(drpver))
    z = numpy.array([get_redshift(plt, ifu, drpall_file)]*nspec)
    print('Redshift: {0}'.format(z[0]))
    dispersion = numpy.full_like(z, 100.)

    # Read the spectra
    print('reading spectra')
    directory_path = os.path.join(os.environ['MANGADAP_DIR'], 'mangadap', 'data', 'remote')
    wave, flux, ivar, sres = get_spectra(plt, ifu, x, y, directory_path=directory_path)
    ferr = numpy.ma.power(ivar, -0.5)

    # Fitting functions expect data to be in 2D arrays (for now):
    if len(flux.shape) == 1:
        flux = flux.reshape(1,-1)
        ferr = ferr.reshape(1,-1)
        sres = sres.reshape(1,-1)
        flux_binned = flux.copy()
        ferr_binned = ferr.copy()
        sres_binned = sres.copy()
        x_binned = x.copy()
        y_binned = y.copy()
        z_binned = z.copy()
        dispersion_binned = dispersion.copy()
    else:
        # Stack the spectra
        wave_binned, flux_binned, fsdev_binned, npix_binned, ivar_binned, sres_binned, \
                covar_binned = SpectralStack().stack(wave, flux, binid=binid, ivar=ivar, sres=sres)
        ferr_binned = numpy.ma.power(ivar_binned, -0.5)
        x_binned = numpy.array([numpy.mean(x[binid == i]) for i in numpy.unique(binid)])
        y_binned = numpy.array([numpy.mean(y[binid == i]) for i in numpy.unique(binid)])
        z_binned = numpy.array([numpy.mean(z[binid == i]) for i in numpy.unique(binid)])
        dispersion_binned = numpy.array([numpy.mean(dispersion)])

    if usr_plots:
        for f in flux:
            pyplot.plot(wave, f)
        pyplot.plot(wave_binned, flux_binned[0])
        pyplot.show()

    #-------------------------------------------------------------------
    # Fit the stellar continuum

    # Mask the 5577 sky line and the emission lines
    sc_pixel_mask = SpectralPixelMask(artdb=ArtifactDB.from_key('BADSKY'),
                                      emldb=EmissionLineDB.from_key('ELPSCMSK'))

    # Construct the template library
    sc_tpl = TemplateLibrary(sc_tpl_key, match_resolution=False, velscale_ratio=velscale_ratio,
                             spectral_step=1e-4, log=True, hardcopy=False)
    sc_tpl_sres = numpy.mean(sc_tpl['SPECRES'].data, axis=0).ravel()

    # Instantiate the fitting class
    ppxf = PPXFFit(StellarContinuumModelBitMask())

    # Perform the fit
    cont_wave, cont_flux, cont_mask, cont_par \
        = ppxf.fit(sc_tpl['WAVE'].data.copy(), sc_tpl['FLUX'].data.copy(), wave_binned, 
                   flux_binned, ferr_binned, z_binned, dispersion_binned,
                   iteration_mode='no_global_wrej', reject_boxcar=100,
                   ensemble=False, velscale_ratio=velscale_ratio, mask=sc_pixel_mask,
                   matched_resolution=False, tpl_sres=sc_tpl_sres, obj_sres=sres_binned,
                   degree=8, moments=2, plot=fit_plots)

    # Remask the continuum fit
    sc_continuum = StellarContinuumModel.reset_continuum_mask_window(
                        numpy.ma.MaskedArray(cont_flux, mask=cont_mask>0))

    # Show the fit and residual
    if usr_plots:
        pyplot.plot(wave, flux_binned[0,:], label='Data')
        pyplot.plot(wave, sc_continuum[0,:], label='Model')
        pyplot.plot(wave, flux_binned[0,:] - sc_continuum[0,:], label='Resid')
        pyplot.legend()
        pyplot.xlabel('Wavelength')
        pyplot.ylabel('Flux')
        pyplot.show()
    #-------------------------------------------------------------------

    #-------------------------------------------------------------------
    # Get the emission-line moments using the fitted stellar continuum

    # Read the database that define the emission lines and passbands
    momdb = EmissionMomentsDB.from_key(elmom_key)

    # Measure the moments
    elmom = EmissionLineMoments.measure_moments(momdb, wave, flux_binned, continuum=sc_continuum,
                                                redshift=z_binned)
    #-------------------------------------------------------------------

    #-------------------------------------------------------------------
    # Fit the emission-line model

    # Set the emission-line continuum templates if different from those
    # used for the stellar continuum
    if sc_tpl_key == el_tpl_key:
        # If the keywords are the same, just copy over the previous
        # library ...
        el_tpl = sc_tpl
        el_tpl_sres = sc_tpl_sres
        # ... and the best fitting stellar kinematics
        stellar_kinematics = cont_par['KIN']
    else:
        # If the template sets are different, we need to match the
        # spectral resolution to the galaxy data ...
        _sres = SpectralResolution(wave, sres[0,:], log10=True)
        el_tpl = TemplateLibrary(el_tpl_key, sres=_sres, velscale_ratio=velscale_ratio,
                                 spectral_step=1e-4, log=True, hardcopy=False)
        el_tpl_sres = numpy.mean(el_tpl['SPECRES'].data, axis=0).ravel()
        # ... and use the corrected velocity dispersions.
        stellar_kinematics = cont_par['KIN']
        stellar_kinematics[:,1] = numpy.ma.sqrt(numpy.square(cont_par['KIN'][:,1]) -
                                                    numpy.square(cont_par['SIGMACORR_EMP']))
    
    # Mask the 5577 sky line
    el_pixel_mask = SpectralPixelMask(artdb=ArtifactDB.from_key('BADSKY'))

    # Read the emission line fitting database
    emldb = EmissionLineDB.from_key(elfit_key)

    # Instantiate the fitting class
    emlfit = Sasuke(EmissionLineModelBitMask())

    # Perform the emission-line fit on each spectrum using the stellar
    # kinematics from the stacked spectrum
    eml_wave, model_flux, eml_flux, eml_mask, eml_fit_par, eml_eml_par \
            = emlfit.fit(emldb, wave_binned, flux_binned, obj_ferr=ferr_binned,
                         obj_mask=el_pixel_mask, obj_sres=sres_binned,
                         guess_redshift=z_binned, guess_dispersion=dispersion_binned,
                         reject_boxcar=101, stpl_wave=el_tpl['WAVE'].data,
                         stpl_flux=el_tpl['FLUX'].data, stpl_sres=el_tpl_sres,
                         stellar_kinematics=stellar_kinematics, etpl_sinst_mode='offset',
                         etpl_sinst_min=10., velscale_ratio=velscale_ratio,
                         matched_resolution=False, mdegree=8, plot=fit_plots,
                         remapid=binid, remap_flux=flux, remap_ferr=ferr,
                         remap_mask=el_pixel_mask, remap_sres=sres, remap_skyx=x, remap_skyy=y,
                         obj_skyx=x_binned, obj_skyy=y_binned)

    # Line-fit metrics
    eml_eml_par = EmissionLineFit.line_metrics(emldb, wave, flux, ferr, model_flux, eml_eml_par,
                                               model_mask=eml_mask, bitmask=emlfit.bitmask)

    # Get the stellar continuum that was fit for the emission lines
    elcmask = eml_mask > 0
    edges = numpy.ma.notmasked_edges(numpy.ma.MaskedArray(model_flux, mask=elcmask), axis=1)
    for i,s,e in zip(edges[0][0],edges[0][1],edges[1][1]):
        elcmask[i,s:e+1] = False
    el_continuum = numpy.ma.MaskedArray(model_flux - eml_flux, mask=elcmask)

    # Plot the result
    if usr_plots:
        for i in range(flux.shape[0]):
            pyplot.plot(wave, flux[i,:], label='Data')
            pyplot.plot(wave, model_flux[i,:], label='Model')
            pyplot.plot(wave, el_continuum[i,:], label='EL Cont.')
            pyplot.plot(wave, sc_continuum[binid[i],:], label='SC Cont.')
            pyplot.legend()
            pyplot.xlabel('Wavelength')
            pyplot.ylabel('Flux')
            pyplot.show()

    # Remeasure the emission-line moments with the new continuum
    new_elmom = EmissionLineMoments.measure_moments(momdb, wave, flux, continuum=el_continuum,
                                                    redshift=z)

    # Compare the summed flux and Gaussian-fitted flux for all the
    # fitted lines
    if usr_plots:
        pyplot.scatter(emldb['restwave'], new_elmom['FLUX'][0,:]-eml_eml_par['FLUX'][0,:],
                       c=eml_eml_par['FLUX'][0,:], cmap='viridis', marker='.', s=60, lw=0,
                       zorder=4)
        pyplot.grid()
        pyplot.xlabel('Wavelength')
        pyplot.ylabel('Summed-Gaussian Difference')
        pyplot.show()

    # TODO: Add the spectral index calls...

    print('Elapsed time: {0} seconds'.format(time.perf_counter() - t))

