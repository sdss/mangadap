
# Imports
import os
import time

from IPython import embed

import numpy

from astropy.io import fits
from matplotlib import pyplot

from mangadap.datacube import MaNGADataCube
from mangadap.config import defaults
from mangadap.util.sampling import spectral_coordinate_step

from mangadap.util.fitsutil import DAPFitsUtil
from mangadap.util.resolution import SpectralResolution
from mangadap.util.pixelmask import SpectralPixelMask

from mangadap.par.artifactdb import ArtifactDB
from mangadap.par.emissionmomentsdb import EmissionMomentsDB
from mangadap.par.emissionlinedb import EmissionLineDB

from mangadap.proc.templatelibrary import TemplateLibrary
from mangadap.proc.emissionlinemoments import EmissionLineMoments
from mangadap.proc.sasuke import Sasuke
from mangadap.proc.ppxffit import PPXFFit
from mangadap.proc.stellarcontinuummodel import StellarContinuumModel, StellarContinuumModelBitMask
from mangadap.proc.emissionlinemodel import EmissionLineModelBitMask
from mangadap.proc.spectralfitting import EmissionLineFit

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
        drpall_file = defaults.drpall_file()
    hdu = fits.open(drpall_file)
    indx = hdu[1].data['PLATEIFU'] == '{0}-{1}'.format(plt, ifu)
    return hdu[1].data['NSA_Z'][indx][0]


def get_spectrum(plt, ifu, x, y, directory_path=None):
    """
    Extract a single spectrum from a MaNGA observation.

    Args:
        plt (:obj:`int`):
            Plate number
        ifu (:obj:`int`):
            IFU identifier
        x (:obj:`int`):
            The spaxel coordinate along the RA axis.
        y (:obj:`int`):
            The spaxel coordinate along the DEC axis.
        directory_path (:obj:`str`, optional):
            Directory with the DRP LOGCUBE file. If None, uses the
            default directory path based on the environmental
            variables.

    Returns:
        :obj:`tuple`: Returns 4 numpy vectors: The wavelength, flux,
        flux inverse variance, and spectral resolution extracted from
        the datacube.
    """
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

    # For testing
    drpver = 'v3_0_1'
    directory_path = os.path.join(os.environ['MANGADAP_DIR'], 'mangadap', 'data', 'remote')

    #-------------------------------------------------------------------
    # Read spectra to fit. The following reads a single MaNGA spectrum.
    # This is where you should read in your own spectrum to fit.
    # Plate-IFU to use
    plt = 7815
    ifu = 6101
    # Spaxel coordinates
    x = 25 #30
    y = 25 #37
    # Read a spectrum
    wave, flux, ivar, sres = get_spectrum(plt, ifu, x, y, directory_path=directory_path)
    # In general, the DAP fitting functions expect data to be in 2D
    # arrays with shape (N-spectra,N-wave). So if you only have one
    # spectrum, you need to expand the dimensions:
    flux = flux.reshape(1,-1)
    ferr = numpy.ma.power(ivar, -0.5).reshape(1,-1)
    sres = sres.reshape(1,-1)

    # The majority (if not all) of the DAP methods expect that your
    # spectra are binned logarithmically in wavelength (primarily
    # because this is what pPXF expects). You can either have the DAP
    # function determine this value (commented line below) or set it
    # directly. The value is used to resample the template spectra to
    # match the sampling of the spectra to fit (up to some integer; see
    # velscale_ratio).
#    spectral_step = spectral_coordinate_step(wave, log=True)
    spectral_step = 1e-4

    # Hereafter, the methods expect a wavelength vector, a flux array
    # with the spectra to fit, an ferr array with the 1-sigma errors in
    # the flux, and sres with the wavelength-dependent spectral
    # resolution, R = lambda / Dlambda
    #-------------------------------------------------------------------

    #-------------------------------------------------------------------
    # The DAP needs a reasonable guess of the redshift of the spectrum
    # (within +/- 2000 km/s). In this example, I'm pulling the redshift
    # from the DRPall file. There must be one redshift estimate per
    # spectrum to fit.  Here that means it's a single element array
    drpall_file = os.path.join(directory_path, 'drpall-{0}.fits'.format(drpver))
    z = numpy.array([get_redshift(plt, ifu, drpall_file)])
    print('Redshift: {0}'.format(z[0]))
    # The DAP also requires an initial guess for the velocity
    # dispersion. A guess of 100 km/s is usually robust, but this may
    # depend on your spectral resolution.
    dispersion = numpy.array([100.])
    #-------------------------------------------------------------------

    #-------------------------------------------------------------------
    # The following sets the keyword for the template spectra to use
    # during the fit. You can specify different template sets to use
    # during the stellar-continuum (stellar kinematics) fit and the
    # emission-line modeling.

    # Templates used in the stellar continuum fits
    sc_tpl_key = 'MILESHC'
    # Templates used in the emission-line modeling
    el_tpl_key = 'MASTARSSP'

    # You also need to specify the sampling for the template spectra.
    # The templates must be sampled with the same pixel step as the
    # spectra to be fit, up to an integer factor. The critical thing
    # for the sampling is that you do not want to undersample the
    # spectral resolution element of the template spectra. Here, I set
    # the sampling for the MILES templates to be a factor of 4 smaller
    # than the MaNGA spectrum to be fit (which is a bit of overkill
    # given the resolution difference). I set the sampling of the
    # MaStar templates to be the same as the galaxy data.

    # Template pixel scale a factor of 4 smaller than galaxy data
    sc_velscale_ratio = 4
    # Template sampling is the same as the galaxy data
    el_velscale_ratio = 1

    # You then need to identify the database that defines the
    # emission-line passbands (elmom_key) for the non-parametric
    # emission-line moment calculations, and the emission-line
    # parameters (elfit_key) for the Gaussian emission-line modeling.
    # See
    # https://sdss-mangadap.readthedocs.io/en/latest/emissionlines.html.
    elmom_key = 'ELBMPL9'
    elfit_key = 'ELPMPL9'

    # If you want to also calculate the spectral indices, you can
    # provide a keyword that indicates the database with the passband
    # definitions for both the absorption-line and bandhead/color
    # indices to measure. The script allows these to be None, if you
    # don't want to calculate the spectral indices. See
    # https://sdss-mangadap.readthedocs.io/en/latest/spectralindices.html
    absindx_key = ''
    bhdindx_key = ''

    # Finally, you can set whether or not to show a set of plots.
    #
    # Show the ppxf-generated plots for each fit stage.
    fit_plots = False
    # Show summary plots
    usr_plots = True
    #-------------------------------------------------------------------



    #-------------------------------------------------------------------
    # Fit the stellar continuum

    # First, we construct the template library. The keyword that
    # selects the template library (sc_tpl_key) is defined above. The
    # following call reads in the template library and processes the
    # data to have the appropriate pixel sampling. Note that *no*
    # matching of the spectral resolution to the galaxy spectrum is
    # performed.
    sc_tpl = TemplateLibrary(sc_tpl_key, match_resolution=False, velscale_ratio=sc_velscale_ratio,
                             spectral_step=spectral_step, log=True, hardcopy=False)

    # This calculation of the mean spectral resolution is a kludge. The
    # template library should provide spectra that are *all* at the
    # same spectral resolution. Otherwise, one cannot freely combine
    # the spectra to fit the Doppler broadening of the galaxy spectrum
    # in a robust (constrained) way (without substantially more
    # effort). There should be no difference between what's done below
    # and simply taking the spectral resolution to be that of the first
    # template spectrum (i.e., sc_tpl['SPECRES'].data[0])
    sc_tpl_sres = numpy.mean(sc_tpl['SPECRES'].data, axis=0).ravel()

    # Now we want to construct a pixel mask that excludes regions with
    # known artifacts and emission lines. The 'BADSKY' artifact
    # database only masks the 5577, which can have strong left-over
    # residuals after sky-subtraction. The list of emission lines (set
    # by the ELPMPL8 keyword) can be different from the list of
    # emission lines fit below.
    sc_pixel_mask = SpectralPixelMask(artdb=ArtifactDB.from_key('BADSKY'),
                                      emldb=EmissionLineDB.from_key('ELPMPL8'))

    # Instantiate the fitting class, including the mask that it should
    # use to flag the data. [[This mask should just be default...]]
    ppxf = PPXFFit(StellarContinuumModelBitMask())

    # The following call performs the fit to the spectrum. Specifically
    # note that the code only fits the first two moments, uses an
    # 8th-order additive polynomial, and uses the 'no_global_wrej'
    # iteration mode. See
    # https://sdss-mangadap.readthedocs.io/en/latest/api/mangadap.proc.ppxffit.html#mangadap.proc.ppxffit.PPXFFit.fit
    cont_wave, cont_flux, cont_mask, cont_par \
        = ppxf.fit(sc_tpl['WAVE'].data.copy(), sc_tpl['FLUX'].data.copy(), wave, flux, ferr,
                   z, dispersion, iteration_mode='no_global_wrej', reject_boxcar=100,
                   ensemble=False, velscale_ratio=sc_velscale_ratio, mask=sc_pixel_mask,
                   matched_resolution=False, tpl_sres=sc_tpl_sres, obj_sres=sres, degree=8,
                   moments=2, plot=fit_plots)

    # The returned objects from the fit are the wavelength, model, and
    # mask vectors and the record array with the best-fitting model
    # parameters. The datamodel of the best-fitting model parameters is
    # set by:
    # https://sdss-mangadap.readthedocs.io/en/latest/api/mangadap.proc.spectralfitting.html#mangadap.proc.spectralfitting.StellarKinematicsFit._per_stellar_kinematics_dtype

#    mod_dcnvlv = PPXFFit.construct_models(sc_tpl['WAVE'].data.copy(), sc_tpl['FLUX'].data.copy(),
#                                          wave, flux.shape, cont_par, redshift_only=True)
#
#    pyplot.plot(wave, cont_flux[0])
#    pyplot.plot(wave, mod_dcnvlv[0])
#    pyplot.show()

    # Remask the continuum fit
    sc_continuum = StellarContinuumModel.reset_continuum_mask_window(
                        numpy.ma.MaskedArray(cont_flux, mask=cont_mask>0))

    # Show the fit and residual
    if usr_plots:
        pyplot.plot(wave, flux[0,:], label='Data')
        pyplot.plot(wave, sc_continuum[0,:], label='Model')
        pyplot.plot(wave, flux[0,:] - sc_continuum[0,:], label='Resid')
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
    elmom = EmissionLineMoments.measure_moments(momdb, wave, flux, continuum=sc_continuum,
                                                redshift=z)
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
        el_tpl = TemplateLibrary(el_tpl_key, sres=_sres, velscale_ratio=el_velscale_ratio,
                                 spectral_step=spectral_step, log=True, hardcopy=False)
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

    # Perform the fit
    efit_t = time.perf_counter()
    eml_wave, model_flux, eml_flux, eml_mask, eml_fit_par, eml_eml_par \
            = emlfit.fit(emldb, wave, flux, obj_ferr=ferr, obj_mask=el_pixel_mask, obj_sres=sres,
                         guess_redshift=z, guess_dispersion=dispersion, reject_boxcar=101,
                         stpl_wave=el_tpl['WAVE'].data, stpl_flux=el_tpl['FLUX'].data,
                         stpl_sres=el_tpl_sres, stellar_kinematics=stellar_kinematics,
                         etpl_sinst_mode='offset', etpl_sinst_min=10.,
                         velscale_ratio=el_velscale_ratio, matched_resolution=False, mdegree=8,
                         plot=fit_plots)
    print('TIME: ', time.perf_counter() - efit_t)

    # Line-fit metrics
    eml_eml_par = EmissionLineFit.line_metrics(emldb, wave, flux, ferr, model_flux, eml_eml_par,
                                               model_mask=eml_mask, bitmask=emlfit.bitmask)

    # Get the stellar continuum that was fit for the emission lines
    elcmask = eml_mask.ravel() > 0
    goodpix = numpy.arange(elcmask.size)[numpy.invert(elcmask)]
    start, end = goodpix[0], goodpix[-1]+1
    elcmask[start:end] = False
    el_continuum = numpy.ma.MaskedArray(model_flux - eml_flux,
                                        mask=elcmask.reshape(model_flux.shape))

    # Plot the result
    if usr_plots:
        pyplot.plot(wave, flux[0,:], label='Data')
        pyplot.plot(wave, model_flux[0,:], label='Model')
        pyplot.plot(wave, el_continuum[0,:], label='EL Cont.')
        pyplot.plot(wave, sc_continuum[0,:], label='SC Cont.')
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
        pyplot.scatter(emldb['restwave'], (new_elmom['FLUX']-eml_eml_par['FLUX']).ravel(),
                       c=eml_eml_par['FLUX'].ravel(), cmap='viridis', marker='.', s=60, lw=0,
                       zorder=4)
        pyplot.grid()
        pyplot.xlabel('Wavelength')
        pyplot.ylabel('Summed-Gaussian Difference')
        pyplot.show()

    print('Elapsed time: {0} seconds'.format(time.perf_counter() - t))

