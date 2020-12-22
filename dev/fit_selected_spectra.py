#!/usr/bin/env python3

import os
import time
import argparse

from IPython import embed

import numpy
from astropy.io import fits
import astropy.constants
from matplotlib import pyplot, ticker

from mangadap.util.constants import DAPConstants
from mangadap.util.filter import interpolate_masked_vector
from mangadap.util.drpfits import DRPFitsBitMask
from mangadap.util.fileio import compress_file
from mangadap.util.resolution import SpectralResolution
from mangadap.util.pixelmask import SpectralPixelMask

from mangadap.par.artifactdb import ArtifactDB
from mangadap.par.emissionmomentsdb import EmissionMomentsDB
from mangadap.par.emissionlinedb import EmissionLineDB

from mangadap.proc.bandpassfilter import emission_line_equivalent_width

from mangadap.proc.templatelibrary import TemplateLibrary
from mangadap.proc.emissionlinemoments import EmissionLineMoments, EmissionLineMomentsBitMask
from mangadap.proc.sasuke import Sasuke
from mangadap.proc.ppxffit import PPXFFit
from mangadap.proc.stellarcontinuummodel import StellarContinuumModel, StellarContinuumModelBitMask
from mangadap.proc.emissionlinemodel import EmissionLineModelBitMask
from mangadap.proc.spectralfitting import EmissionLineFit

#-----------------------------------------------------------------------------
def object_data(f, fit_flag_file):

    with fits.open(f) as hdu:
        flux = hdu['FLUX'].data
        ivar = hdu['IVAR'].data
        mask = hdu['MASK'].data
        sres = hdu['LSF'].data
        wave = hdu['WAVE'].data
        redshift = hdu['Z'].data

    sres = numpy.ma.divide(wave[None,:], sres) / DAPConstants.sig2fwhm
    if numpy.any(numpy.ma.getmaskarray(sres)):
        sres = numpy.ma.apply_along_axis(interpolate_masked_vector, 1, sres)

    bm = DRPFitsBitMask()
    nspec = flux.shape[0]

    fit_flags = numpy.ones(nspec, dtype=bool) if fit_flag_file is None \
                    else (numpy.genfromtxt(fit_flag_file)[:,3]).astype(bool)
    if len(fit_flags) != nspec:
        print(nspec, len(fit_flags))
        raise ValueError('Incorrect number of fitting flags.')

    ferr = numpy.ma.power(ivar, -0.5)
    bool_mask = bm.flagged(mask, flag=['DONOTUSE', 'FORESTAR']) | numpy.ma.getmaskarray(ferr)
    
    return wave, numpy.ma.MaskedArray(flux, mask=bool_mask), \
                numpy.ma.MaskedArray(ferr, mask=bool_mask), sres, redshift, fit_flags

def init_ax(fig, pos):
    ax = fig.add_axes(pos) #, facecolor='0.95')
    ax.minorticks_on()
    ax.tick_params(which='major', length=10)
    ax.tick_params(which='minor', length=5)
    ax.tick_params(which='both', direction='out', top='on', right='on')
    return ax

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('inp', type=str, help='input fits file with spectra to fit')
    parser.add_argument('out', type=str, help='output fits file with fitted models')

    parser.add_argument('--output_root', type=str, help='root directory for output', default=None)
    parser.add_argument('--sc_tpl', type=str, default='MILESHC',
                        help='template library to use for the stellar kinematics')
    parser.add_argument('--sc_vsr', type=int, default=4,
                        help='velocity-scale ratio for templates during stellar kinematics fit')
    parser.add_argument('--sc_deg', type=int, default=8,
                        help='additive order during stellar kinematics fit')
    parser.add_argument('--sc_only', default=False, action='store_true',
                        help='skip the emission-line modeling and only fit the stellar kinematics')

    parser.add_argument('--el_tpl', type=str, default='MILESHC',
                        help='template library to use for the emission-line fitting')
    parser.add_argument('--el_band', type=str, default='ELBMPL9',
                        help='keyword with the emission-line bandpasses for the moment analysis.')
    parser.add_argument('--el_list', type=str, default='ELPMPL9',
                        help='keyword with the emission-line database to fit')
    parser.add_argument('--el_vsr', type=int, default=4,
                        help='velocity-scale ratio for templates during emission-line fit')
    parser.add_argument('--el_deg', type=int, default=8,
                        help='multiplicative order during emission-line fit')

    parser.add_argument('--spec_flags', type=str, default=None,
                        help='fourth column sets whether (1) or not (0) to fit spectrum')

    return parser.parse_args()
    

def write(ofile, wave, cont_flux, cont_mask, cont_par, model_flux=None, model_mask=None,
          eml_flux=None, eml_fit_par=None, eml_eml_par=None, elmom=None):

    compress = False
    if ofile.split('.')[-1] == 'gz':
        _ofile = ofile[:ofile.rfind('.')] 
        compress = True
    else:
        _ofile = ofile

    hdu = fits.HDUList([fits.PrimaryHDU(),
                        fits.ImageHDU(data=wave, name='WAVE'),
                        fits.ImageHDU(data=cont_flux.data, name='STELLAR'),
                        fits.ImageHDU(data=cont_mask, name='STELLAR_MASK'),
                        cont_par.to_hdu(name='STRPAR')])

    if model_flux is not None:
        hdu += [fits.ImageHDU(data=model_flux.data, name='MODEL'),
                fits.ImageHDU(data=model_mask, name='MODEL_MASK'),
                fits.ImageHDU(data=eml_flux.data, name='EMLINE'),
                eml_fit_par.to_hdu(name='EMLFIT'),
                eml_eml_par.to_hdu(name='EMLPAR'),
                elmom.to_hdu(name='EMLMOM')]

    hdu.writeto(_ofile, overwrite=True)
    
    # Compress if desired
    if compress:
        compress_file(_ofile, clobber=True)
        os.remove(_ofile)


def main():
    t = time.perf_counter()
    arg = parse_args()
    if not os.path.isfile(arg.inp):
        raise FileNotFoundError('No file: {0}'.format(arg.inp))
    directory_path = os.getcwd() if arg.output_root is None else os.path.abspath(arg.output_root)
    if not os.path.isdir(directory_path):
        os.makedirs(directory_path)

    data_file = os.path.abspath(arg.inp)
    fit_file = os.path.join(directory_path, arg.out)
    flag_db = None if arg.spec_flags is None else os.path.abspath(arg.spec_flags)

    # Read the data
    spectral_step = 1e-4
    wave, flux, ferr, sres, redshift, fit_spectrum = object_data(data_file, flag_db)
    nspec, npix = flux.shape
    dispersion = numpy.full(nspec, 100., dtype=numpy.float)

#    fit_spectrum[:] = False
#    fit_spectrum[0] = True
#    fit_spectrum[171] = True
#    fit_spectrum[791] = True

    # Mask spectra that should not be fit
    indx = numpy.any(numpy.logical_not(numpy.ma.getmaskarray(flux)), axis=1) & fit_spectrum
    flux[numpy.logical_not(indx),:] = numpy.ma.masked

    print('Read: {0}'.format(arg.inp))
    print('Contains {0} spectra'.format(nspec))
    print(' each with {0} pixels'.format(npix))
    print('Fitting {0} spectra.'.format(numpy.sum(fit_spectrum)))

    #-------------------------------------------------------------------
    #-------------------------------------------------------------------
    # Fit the stellar continuum

    # Construct the template library
    sc_tpl = TemplateLibrary(arg.sc_tpl, match_resolution=False, velscale_ratio=arg.sc_vsr,
                             spectral_step=spectral_step, log=True, hardcopy=False)
    # Set the spectral resolution
    sc_tpl_sres = numpy.mean(sc_tpl['SPECRES'].data, axis=0).ravel()
    # Set the pixel mask
    sc_pixel_mask = SpectralPixelMask(artdb=ArtifactDB.from_key('BADSKY'),
                                      emldb=EmissionLineDB.from_key('ELPMPL8'))
    # Instantiate the fitting class
    ppxf = PPXFFit(StellarContinuumModelBitMask())

    # The following call performs the fit to the spectrum. Specifically
    # note that the code only fits the first two moments, uses an
    # 8th-order additive polynomial, and uses the 'no_global_wrej'
    # iteration mode. See
    # https://sdss-mangadap.readthedocs.io/en/latest/api/mangadap.proc.ppxffit.html#mangadap.proc.ppxffit.PPXFFit.fit
    cont_wave, cont_flux, cont_mask, cont_par \
        = ppxf.fit(sc_tpl['WAVE'].data.copy(), sc_tpl['FLUX'].data.copy(), wave, flux, ferr,
                   redshift, dispersion, iteration_mode='no_global_wrej', reject_boxcar=100,
                   ensemble=False, velscale_ratio=arg.sc_vsr, mask=sc_pixel_mask,
                   matched_resolution=False, tpl_sres=sc_tpl_sres, obj_sres=sres,
                   degree=arg.sc_deg, moments=2) #, plot=True)

    if arg.sc_only:
        write(fit_file, wave, cont_flux, cont_mask, cont_par)
        print('Elapsed time: {0} seconds'.format(time.perf_counter() - t))
        return

#    if numpy.any(cont_par['KIN'][:,1] < 0):
#        embed()
#        exit()
    #-------------------------------------------------------------------


    #-------------------------------------------------------------------
    #-------------------------------------------------------------------
    # Measure the emission-line moments

#    # Remask the continuum fit
#    sc_continuum = StellarContinuumModel.reset_continuum_mask_window(
#                        numpy.ma.MaskedArray(cont_flux, mask=cont_mask>0))
#    # Read the database that define the emission lines and passbands
#    momdb = EmissionMomentsDB.from_key(arg.el_band)
#    # Measure the moments
#    elmom = EmissionLineMoments.measure_moments(momdb, wave, flux, continuum=sc_continuum,
#                                                redshift=redshift)
    #-------------------------------------------------------------------

    #-------------------------------------------------------------------
    # Fit the emission-line model

    # Set the emission-line continuum templates if different from those
    # used for the stellar continuum
    if arg.sc_tpl == arg.el_tpl:
        # If the keywords are the same, just copy over the previous
        # library and the best fitting stellar kinematics
        el_tpl = sc_tpl
        el_tpl_sres = sc_tpl_sres
        stellar_kinematics = cont_par['KIN'].copy()
    else:
        # If the template sets are different, we need to match the
        # spectral resolution to the galaxy data and use the corrected
        # velocity dispersions.
        _sres = SpectralResolution(wave, sres[0,:], log10=True)
        el_tpl = TemplateLibrary(arg.el_tpl, sres=_sres, velscale_ratio=arg.el_vsr,
                                 spectral_step=spectral_step, log=True, hardcopy=False)
        el_tpl_sres = numpy.mean(el_tpl['SPECRES'].data, axis=0).ravel()
        stellar_kinematics = cont_par['KIN'].copy()
        stellar_kinematics[:,1] = numpy.ma.sqrt(numpy.square(cont_par['KIN'][:,1]) -
                                        numpy.square(cont_par['SIGMACORR_SRES'])).filled(0.0)

#    if numpy.any(cont_par['KIN'][:,1] < 0):
#        embed()
#        exit()
#
#    if numpy.any(stellar_kinematics[:,1] < 0):
#        embed()
#        exit()

    # Mask the 5577 sky line
    el_pixel_mask = SpectralPixelMask(artdb=ArtifactDB.from_key('BADSKY'))

    # Read the emission line fitting database
    emldb = EmissionLineDB.from_key(arg.el_list)

    # Instantiate the fitting class
    emlfit = Sasuke(EmissionLineModelBitMask())

    # TODO: Improve the initial velocity guess using the first moment...

    # Perform the fit
    elfit_time = time.perf_counter()
    model_wave, model_flux, eml_flux, model_mask, eml_fit_par, eml_eml_par \
            = emlfit.fit(emldb, wave, flux, obj_ferr=ferr, obj_mask=el_pixel_mask, obj_sres=sres,
                         guess_redshift=redshift, guess_dispersion=dispersion, reject_boxcar=101,
                         stpl_wave=el_tpl['WAVE'].data, stpl_flux=el_tpl['FLUX'].data,
                         stpl_sres=el_tpl_sres, stellar_kinematics=stellar_kinematics,
                         etpl_sinst_mode='offset', etpl_sinst_min=10., velscale_ratio=arg.el_vsr,
                         matched_resolution=False, mdegree=arg.el_deg, ensemble=False)#, plot=True)
    print('EML FIT TIME: ', time.perf_counter() - elfit_time)

    # Line-fit metrics (should this be done in the fit method?)
    eml_eml_par = EmissionLineFit.line_metrics(emldb, wave, flux, ferr, model_flux, eml_eml_par,
                                               model_mask=model_mask, bitmask=emlfit.bitmask)

    # Equivalent widths
    EmissionLineFit.measure_equivalent_width(wave, flux, emldb, eml_eml_par,
                                             bitmask=emlfit.bitmask, checkdb=False)

    # Measure the emission-line moments
    #   - Model continuum
    continuum = StellarContinuumModel.reset_continuum_mask_window(model_flux - eml_flux)
    #   - Updated redshifts
    fit_redshift = eml_eml_par['KIN'][:,numpy.where(emldb['name'] == 'Ha')[0][0],0] \
                        / astropy.constants.c.to('km/s').value
    #   - Set the moment database
    momdb = EmissionMomentsDB.from_key(arg.el_band)
    #   - Set the moment bitmask
    mombm = EmissionLineMomentsBitMask()
    #   - Measure the moments
    elmom = EmissionLineMoments.measure_moments(momdb, wave, flux, ivar=numpy.ma.power(ferr, -2),
                                                continuum=continuum, redshift=fit_redshift,
                                                bitmask=mombm)
    #   - Select the bands that are valid
    include_band = numpy.array([numpy.logical_not(momdb.dummy)]*nspec) \
                        & numpy.logical_not(mombm.flagged(elmom['MASK'],
                                                          flag=['BLUE_EMPTY', 'RED_EMPTY']))
    #   - Set the line center at the center of the primary passband
    line_center = (1.0+fit_redshift)[:,None]*momdb['restwave'][None,:]
    elmom['BMED'], elmom['RMED'], pos, elmom['EWCONT'], elmom['EW'], elmom['EWERR'] \
            = emission_line_equivalent_width(wave, flux, momdb['blueside'], momdb['redside'],
                                             line_center, elmom['FLUX'], redshift=fit_redshift,
                                             line_flux_err=elmom['FLUXERR'],
                                             include_band=include_band)
    #   - Flag non-positive measurements
    indx = include_band & numpy.logical_not(pos)
    elmom['MASK'][indx] = mombm.turn_on(elmom['MASK'][indx], 'NON_POSITIVE_CONTINUUM')
    #   - Set the binids
    elmom['BINID'] = numpy.arange(nspec)
    elmom['BINID_INDEX'] = numpy.arange(nspec)

    write(fit_file, wave, cont_flux, cont_mask, cont_par, model_flux=model_flux,
          model_mask=model_mask, eml_flux=eml_flux, eml_fit_par=eml_fit_par,
          eml_eml_par=eml_eml_par, elmom=elmom)
    print('Elapsed time: {0} seconds'.format(time.perf_counter() - t))


if __name__ == '__main__':
    main()