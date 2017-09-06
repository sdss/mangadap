#!/usr/bin/env python3

import time
import warnings
import numpy
import astropy.constants

from astropy.io import fits

from captools import ppxf_util

from mangadap.drpfits import DRPFits
from mangadap.proc.reductionassessments import ReductionAssessment
from mangadap.proc.spatiallybinnedspectra import SpatiallyBinnedSpectra
from mangadap.proc.stellarcontinuummodel import StellarContinuumModel
from mangadap.proc.ppxffit import PPXFFit
from mangadap.util.instrument import spectrum_velocity_scale
from mangadap.util.fitsutil import DAPFitsUtil
from mangadap.util.fileio import init_record_array
from mangadap.proc.emissionlinemoments import EmissionLineMoments
from mangadap.proc.spectralfitting import EmissionLineFit
from mangadap.proc.emissionlinemodel import EmissionLineModelDef, EmissionLineModel
from mangadap.par.emissionlinedb import EmissionLineDB
from mangadap.proc.spectralindices import SpectralIndices
from mangadap.dapfits import construct_maps_file, construct_cube_file

from mangadap.contrib.xjmc import emline_fitter_with_ppxf

#-----------------------------------------------------------------------------

class XJMCEmissionLineFitter(EmissionLineFit):
    def __init__(self, par=None):
        if par is None:
            # Set the default parameter set.  The guess_redshift,
            # stellar_continuum, and emission_lines values can be filled
            # by EmissionLineModel._fill_method_par()
            par = { 'guess_redshift': None,     # The guess redshift for each binned spectrum
                    'stellar_continuum': None,  # The StellarContinuumModel object
                    'emission_lines': None,     # The EmissionLineDB object
                    'degree': -1,                # Additive polynomial order
                    'mdegree': 10 }              # Multiplicative polynomial order
        EmissionLineFit.__init__(self, 'XJMC', None, par=par)


    def line_database(self):
        """
        For use in the output to distinguish the lines being fit.
        """
        neml = 11
        emission_line_name = [ 'Hd', 'Hg', 'Hb', 'Ha', 'OII', 'OII', 'SII', 'SII', 'OIIId', 'OId',
                               'NIId']
        emission_line_restwave = [ 4103, 4342, 4861, 6564, 3727, 3729, 6718, 6732, 4984, 6334,
                                   6567 ]

        name_len = numpy.amax([ len(n) for n in emission_line_name ])
        db = init_record_array(neml, [ ('NAME','<U{0:d}'.format(name_len)), ('RESTWAVE', float) ])
        db['NAME'] = emission_line_name
        db['RESTWAVE'] = emission_line_restwave
        return db
        
    
    def fit(self, binned_spectra, par=None, loggers=None, quiet=False):
        if par is not None:
            self.par = par

        # Check the parameter keys
        required_keys = [ 'guess_redshift', 'stellar_continuum', 'emission_lines', 'degree',
                          'mdegree' ]
        if numpy.any([ reqk not in self.par.keys() for reqk in required_keys ]):
            raise ValueError('Parameter dictionary does not have all the required keys.')

        # Wavelengths are in vacuum
        wave0 = binned_spectra['WAVE'].data.copy()
        # Velocity step per pixel
        velscale = spectrum_velocity_scale(wave0)
        
        # Get the best-fitting stellar kinematics for the binned spectra
        # And correct sigmas with instrumental resolutions (if convolved templates are to be used)
        if self.par['stellar_continuum'] is None \
                or not isinstance(self.par['stellar_continuum'], StellarContinuumModel):
            raise ValueError('Must provide StellarContinuumModel object as the '
                             '\'stellar_continuum\' item in the parameter dictionary')

        stars_vel, stars_sig = self.par['stellar_continuum'].matched_guess_kinematics(
                                            binned_spectra, cz=True, corrected=True, nearest=True)

        # Convert the input stellar velocity from redshift (c*z) to ppxf velocity (c*log(1+z))
        stars_vel = PPXFFit.revert_velocity(stars_vel,0)[0]
        
        # Get the stellar templates and the template resolution;
        # shape is (Nstartpl, Ntplwave)
        stars_templates = self.par['stellar_continuum'].get_template_library(
                            velocity_offset=numpy.median(stars_vel), match_to_drp_resolution=True)
        stars_templates_wave = stars_templates['WAVE'].data.copy()
        template_sres = stars_templates['SPECRES'].data[0,:]
        stars_templates = stars_templates['FLUX'].data.copy()            
        velscale_ratio = self.par['stellar_continuum'].method['fitpar']['velscale_ratio']

        # Set mask for the galaxy spectra according to the templates wave range
        mask = PPXFFit.fitting_mask(tpl_wave=stars_templates_wave, obj_wave=wave0,
                                    velscale=velscale, velscale_ratio=velscale_ratio,
                                    velocity_offset=numpy.median(stars_vel))[0]
        wave = wave0[mask]

        # Calculate the velocity offset between the masked spectra and the tempaltes
        dv = -PPXFFit.ppxf_tpl_obj_voff(stars_templates_wave, wave, velscale,
                                            velscale_ratio=velscale_ratio)
        
        # UNBINNED DATA:
        # Flux and noise masked arrays; shape is (Nspaxel,Nwave) where
        # Nspaxel is Nx*Ny
        # Bin ID from VOR10 reference file used to mask buffer spaxels
        #binid0 = binned_spectra['BINID'].data
        # Create mask for parsing the input data
        #binid = binid0.reshape(-1)
        #mask_spaxel = ~(binid==-1)
        # pPXF would run into problems if dircetly masked arrays were used
        # So both fluxes and their masks are to be provided as input

#        flux00 = binned_spectra.drpf.copy_to_masked_array(flag=['DONOTUSE', 'FORESTAR'])
#        mask_drp0 = ~flux00.mask
#        flux = binned_spectra.drpf.copy_to_array(ext='FLUX')
#        ivar = binned_spectra.drpf.copy_to_array(ext='IVAR')
#        flux0, ivar0 = binned_spectra.galext.apply(flux, ivar=ivar, deredden=True)
#        noise0 = numpy.power(ivar0.data, -0.5)
#        mask_drp = mask_drp0.reshape(-1, mask_drp0.shape[-1])[:,mask]
#        flux = flux0.data[:,mask]
#        noise = noise0[:,mask]

        flux0 = binned_spectra.drpf.copy_to_masked_array(flag=['DONOTUSE', 'FORESTAR'])
        ivar0 = binned_spectra.drpf.copy_to_masked_array(ext='IVAR', flag=['DONOTUSE', 'FORESTAR'])
        flux0, ivar0 = binned_spectra.galext.apply(flux0, ivar=ivar0, deredden=True)
        noise = numpy.ma.power(ivar0, -0.5)
        noise[numpy.invert(noise > 0)] = numpy.ma.masked

        mask_drp = numpy.invert(flux0.mask | noise.mask)[:,mask]
        flux = flux0.data[:,mask]
        noise = noise.filled(0.0)[:,mask]

        # stack_sres sets whether or not the spectral resolution is
        # determined on a per-spaxel basis or with a single vector
        sres = binned_spectra.drpf.spectral_resolution(toarray=True, fill=True) \
                    if binned_spectra.method['stackpar']['stack_sres'] else \
                    binned_spectra.drpf.spectral_resolution(ext='SPECRES', toarray=True, fill=True)
        sres = sres[:,mask]

        # Spaxel coordinates; shape is (Nspaxel,)
        x = binned_spectra.rdxqa['SPECTRUM'].data['SKY_COO'][:,0]
        y = binned_spectra.rdxqa['SPECTRUM'].data['SKY_COO'][:,1]

        # BINNED DATA:
        # Binned flux and binned noise masked arrays; shape is (Nbin,Nwave)

#        flux_binned00 = binned_spectra.copy_to_masked_array(flag=binned_spectra.do_not_fit_flags())
#        mask_binned0 = ~flux_binned00.mask
#        flux_binned0 = binned_spectra.copy_to_array(ext='FLUX')
#        noise_binned0 = numpy.power(binned_spectra.copy_to_array(ext='IVAR'), -0.5)
#        mask_binned = mask_binned0[:,mask]
#        flux_binned = flux_binned0[:,mask]
#        noise_binned = noise_binned0[:,mask]
#        sres_binned0 = binned_spectra.copy_to_array(ext='SPECRES')
#        sres_binned = sres_binned0[:,mask]

        flux_binned = binned_spectra.copy_to_masked_array(flag=binned_spectra.do_not_fit_flags())
        noise_binned = numpy.ma.power(binned_spectra.copy_to_masked_array(ext='IVAR',
                                            flag=binned_spectra.do_not_fit_flags()), -0.5)
        noise_binned[numpy.invert(noise_binned > 0)] = numpy.ma.masked
        mask_binned = numpy.invert(flux_binned.mask | noise_binned.mask)[:,mask]
        flux_binned = flux_binned.data[:,mask]
        noise_binned = noise_binned.filled(0.0)[:,mask]
        sres_binned = binned_spectra.copy_to_array(ext='SPECRES')[:,mask]

        # Bin coordinates; shape is (Nbin,)
        x_binned = binned_spectra['BINS'].data['SKY_COO'][:,0]
        y_binned = binned_spectra['BINS'].data['SKY_COO'][:,1]
        
        # Set initial guesses for the velocity and velocity dispersion
        if self.par['guess_redshift'] is not None:
            # Use guess_redshift if provided
            guess_vel0 = self.par['guess_redshift'] * astropy.constants.c.to('km/s').value
            guess_vel = PPXFFit.revert_velocity(guess_vel0,0)[0]
            # And set default velocity dispersion to 100 km/s
            guess_sig = numpy.full(guess_vel.size, 100, dtype=float)
        elif self.par['stellar_continuum'] is not None:
            # Otherwise use the stellar-continuum result
            guess_vel, guess_sig = stars_vel.copy(), stars_sig.copy()
        else:
            raise ValueError('Cannot set guess kinematics; must provide either \'guess_redshift\' '
                             'or \'stellar_continuum\' in input parameter dictionary.')
        
        # Construct gas templates; shape is (Ngastpl, Ntplwave).
        # Template resolution matched between the stellar and gas
        # templates?
        # Decide whether to use convolved gas templates
        # Set the wavelength of lines to be in vacuum
        FWHM = wave/numpy.max(sres, axis=0)
        FWHM_binned = wave/sres_binned[0,:]
        def fwhm_drp(wave_len0):
            wave_len = wave_len0*(1 + numpy.median(self.par['guess_redshift']))
            index = numpy.argmin(abs(wave-wave_len[:,None]),axis=1) \
                        if numpy.asarray(wave_len) is wave_len \
                        else numpy.argmin(abs(wave-wave_len))
            return FWHM[index]
        def fwhm_binned(wave_len0):
            wave_len = wave_len0*(1 + numpy.median(self.par['guess_redshift']))
            index = numpy.argmin(abs(wave-wave_len[:,None]),axis=1) \
                        if numpy.asarray(wave_len) is wave_len \
                        else numpy.argmin(abs(wave-wave_len))
            return FWHM_binned[index]
        lam_range_gal = numpy.array([numpy.min(wave), numpy.max(wave)]) \
                            / (1 + numpy.median(self.par['guess_redshift']))
        gas_templates, gas_names, gas_wave = \
            ppxf_util.emission_lines(numpy.log(stars_templates_wave), lam_range_gal, fwhm_drp)
        gas_templates_binned, gas_names, gas_wave = \
            ppxf_util.emission_lines(numpy.log(stars_templates_wave), lam_range_gal, fwhm_binned)

        # Default polynomial orders
        degree = -1 if self.par['degree'] is None else self.par['degree']
        mdegree = 10 if self.par['mdegree'] is None else self.par['mdegree']
     
        # --------------------------------------------------------------
        # CALL TO EMLINE_FITTER_WITH_PPXF:
        # Input is:
        #   - wave: wavelength vector; shape is (Nwave,)
        #   - flux: observed, unbinned flux; masked array with shape
        #     (Nspaxel,Nwave)
        #   - noise: error in observed, unbinned flux; masked array with
        #     shape (Nspaxel,Nwave)
        #   - sres: spectral resolution (R=lambda/delta lambda) as a
        #     function of wavelength for each unbinned spectrum; shape
        #     is (Nspaxel,Nwave)
        #   - flux_binned: binned flux; masked array with shape
        #     (Nbin,Nwave)
        #   - noise_binned: noise in binned flux; masked array with
        #     shape (Nbin,Nwave)
        #   - sres_binned: spectral resolution (R=lambda/delta lambda)
        #     as a function of wavelength for each binned spectrum;
        #     shape is (Nbin,Nwave)
        #   - velscale: Velocity step per pixel
        #   - velscale_ratio: Ratio of velocity step per pixel in the
        #     observed data versus in the template data
        #   - dv: Velocity offset between the galaxy and template data
        #     due to the difference in the initial wavelength of the
        #     spectra
        #   - stars_vel: Velocity of the stellar component; shape is
        #     (Nbins,)
        #   - stars_sig: Velocity dispersion of the stellar component;
        #     shape is (Nbins,)
        #   - stars_templates: Stellar templates; shape is (Nstartpl,
        #     Ntplwave)
        #   - guess_vel: Initial guess velocity for the gas components;
        #     shape is (Nbins,)
        #   - guess_sig: Initial guess velocity dispersion for the gas
        #     components; shape is (Nbins,)
        #   - gas_templates: Gas template flux; shape is (Ngastpl,
        #     Ntplwave)
        #   - gas_names: Name of the gas templats; shape is (Ngastpl,)
        #   - template_sres: spectral resolution (R=lambda/delta
        #     lambda) as a function of wavelength for all the templates
        #     templates; shape is (Ntplwave,)
        #   - degree: Additive polynomial order
        #   - mdegree: Multiplicative polynomial order
        #   - x: On-sky spaxel x coordinates; shape is (Nspaxel,)
        #   - y: On-sky spaxel y coordinates; shape is (Nspaxel,)
        #   - x_binned: On-sky bin x coordinate; shape is (Nbin,)
        #   - y_binned: On-sky bin y coordinate; shape is (Nbin,)
        model_flux0, model_eml_flux0, model_mask0, model_binid, eml_flux, eml_fluxerr, \
                eml_kin, eml_kinerr, eml_sigmacorr \
                        = emline_fitter_with_ppxf(wave, flux, noise, sres, flux_binned,
                                                  noise_binned, velscale, velscale_ratio, dv,
                                                  stars_vel, stars_sig, stars_templates,
                                                  guess_vel, guess_sig, gas_templates,
                                                  gas_templates_binned, gas_names, template_sres,
                                                  degree, mdegree, x, y, x_binned, y_binned,
                                                  mask_binned, mask_drp,
                                                  numpy.median(self.par['guess_redshift']))
                                                  #, debug=True)
        # Output is:
        #   - model_flux: stellar-continuum + emission-line model; shape
        #     is (Nmod, Nwave); first axis is ordered by model ID number
        #   - model_eml_flux: model emission-line flux only; shape is
        #     (Nmod, Nwave); first axis is ordered by model ID number
        #   - model_mask: boolean or bit mask for fitted models; shape
        #     is (Nmod, Nwave); first axis is ordered by model ID number
        #   - model_binid: ID numbers assigned to each spaxel with a
        #     fitted model; any spaxel without a model should have
        #     model_binid = -1; the number of >-1 IDs must be Nmod;
        #     shape is (Nx,Ny) which is equivalent to:
        #       flux[:,0].reshape((numpy.sqrt(Nspaxel).astype(int),)*2).shape
        #   - eml_flux: Flux of each emission line; shape is (Nmod,Neml)
        #   - eml_fluxerr: Error in emission-line fluxes; shape is
        #     (Nmod, Neml)
        #   - eml_kin: Kinematics (velocity and velocity dispersion) of
        #     each emission line; shape is (Nmod,Neml,Nkin)
        #   - eml_kinerr: Error in the kinematics of each emission line
        #   - eml_sigmacorr: Quadrature corrections required to obtain
        #     the astrophysical velocity dispersion; shape is
        #     (Nmod,Neml); corrections are expected to be applied as
        #     follows:
        #       sigma = numpy.ma.sqrt( numpy.square(eml_kin[:,:,1])
        #                               - numpy.square(eml_sigmacorr))
        # --------------------------------------------------------------
        
        # Convert the output velocity back from ppxf velocity (c*log(1+z)) to redshift (c*z)
        eml_kin[:,:,0] = PPXFFit.convert_velocity(eml_kin[:,:,0],0)[0]
        
        # Mask output data according to model_binid
        model_binid = numpy.asarray(model_binid, dtype=numpy.int16)
        mask_id = (model_binid > -1)
        model_flux0 = model_flux0[mask_id]
        model_eml_flux0 = model_eml_flux0[mask_id]
        model_mask0 = model_mask0[mask_id]
        eml_flux = eml_flux[mask_id]
        eml_fluxerr = eml_fluxerr[mask_id]
        eml_kin = eml_kin[mask_id]
        eml_kinerr = eml_kinerr[mask_id]
        eml_sigmacorr = eml_sigmacorr[mask_id]
        
        #cube_binid = numpy.full_like(mask_spaxel, -1, dtype=numpy.int16)
        #cube_binid[mask_spaxel] = model_binid
        Nspaxel = x.shape[0]
        model_binid = model_binid.reshape((numpy.sqrt(Nspaxel).astype(int),)*2)
        # The ordered indices in the flatted bin ID map with/for each model
        model_srt = numpy.argsort(model_binid.ravel())[model_binid.ravel() > -1]

        # Construct the output emission-line database.  The data type
        # defined by EmissionLineFit._per_emission_line_dtype(); shape
        # is (Nmod,); parameters must be ordered by model ID number
        nmod = len(model_srt)
        neml = eml_flux.shape[1]
        nkin = eml_kin.shape[-1]
        model_eml_par = init_record_array(nmod,
                                EmissionLineFit._per_emission_line_dtype(neml, nkin, numpy.int16))
        model_eml_par['BINID'] = model_binid.ravel()[model_srt]
        model_eml_par['BINID_INDEX'] = numpy.arange(nmod)
        model_eml_par['MASK'][:,:] = 0
        model_eml_par['FLUX'] = eml_flux
        model_eml_par['FLUXERR'] = eml_fluxerr
        model_eml_par['KIN'] = eml_kin
        model_eml_par['KINERR'] = eml_kinerr
        model_eml_par['SIGMACORR'] = eml_sigmacorr

        # Include the equivalent width measurements
        if self.par['emission_lines'] is not None:
            EmissionLineFit.measure_equivalent_width(wave, flux[model_srt,:],
                                                     par['emission_lines'], model_eml_par)
        
        # Change back the wavelength range of models to match that of the galaxy's
        model_flux = numpy.zeros(flux0[mask_id].shape)
        model_eml_flux = numpy.zeros(flux0[mask_id].shape)
        model_mask = numpy.full_like(flux0[mask_id], 0, dtype=bool)
        model_flux[:,mask] = model_flux0
        model_eml_flux[:,mask] = model_eml_flux0
        model_mask[:,mask] = model_mask0
        
        # Calculate the "emission-line baseline" as the difference
        # between the stellar continuum model determined for the
        # kinematics and the one determined by the optimized
        # stellar-continuum + emission-line fit:
        if self.par['stellar_continuum'] is not None:
            # Construct the full 3D cube for the stellar continuum
            # models
            sc_model_flux, sc_model_mask \
                    = DAPFitsUtil.reconstruct_cube(binned_spectra.drpf.shape,
                                            self.par['stellar_continuum']['BINID'].data.ravel(),
                                            [ self.par['stellar_continuum']['FLUX'].data,
                                              self.par['stellar_continuum']['MASK'].data ])
            # Set any masked pixels to 0
            sc_model_flux[sc_model_mask>0] = 0.0

            # Construct the full 3D cube of the new stellar continuum
            # from the combined stellar-continuum + emission-line fit
            el_continuum = DAPFitsUtil.reconstruct_cube(binned_spectra.drpf.shape,
                                                        model_binid.ravel(),
                                                        model_flux - model_eml_flux)
            # Get the difference, restructure it to match the shape
            # of the emission-line models, and zero any masked pixels
            model_eml_base = (el_continuum - sc_model_flux).reshape(-1,wave0.size)[model_srt,:]
            if model_mask is not None:
                model_eml_base[model_mask==0] = 0.0
        else:
            model_eml_base = numpy.zeros(model_flux.shape, dtype=float)

        # Returned arrays are:
        #   - model_eml_flux: model emission-line flux only; shape is
        #     (Nmod, Nwave); first axis is ordered by model ID number
        #   - model_eml_base: difference between the combined fit and the
        #     stars-only fit; shape is (Nmod, Nwave); first axis is
        #     ordered by model ID number
        #   - model_mask: boolean or bit mask for fitted models; shape
        #     is (Nmod, Nwave); first axis is ordered by model ID number
        #   - model_fit_par: This provides the results of each fit;
        #     TODO: The is set to None.  Provide metrics of the ppxf fit
        #     to each spectrum?
        #   - model_eml_par: output model parameters; data type must be
        #     EmissionLineFit._per_emission_line_dtype(); shape is
        #     (Nmod,); parameters must be ordered by model ID number
        #   - model_binid: ID numbers assigned to each spaxel with a
        #     fitted model; any spaxel with a model should have
        #     model_binid = -1; the number of >-1 IDs must be Nmod;
        #     shape is (Nx,Ny)
        return model_eml_flux, model_eml_base, model_mask, None, model_eml_par, model_binid

    
#-----------------------------------------------------------------------------
if __name__ == '__main__':
    t = time.clock()

    # Set the plate, ifu, and initial velocity/redshift
    plate = 7815
    ifu = 6101
    # KBW testing edits
    vel = 5144.7681
    nsa_redshift = vel/astropy.constants.c.to('km/s').value
#    hdu_cat = fits.open('/Users/mac/Downloads/Monty_Python/Catalog/mpl5_cat.fits')
#    mask_z = [i['PLATEIFU'] == str(plate)+'-'+str(ifu) for i in hdu_cat[1].data]
#    nsa_redshift = hdu_cat[1].data[mask_z][0]['NSA_Z']
#    vel = nsa_redshift*astropy.constants.c.to('km/s').value
    
    # KBW testing edits
#    redux_path = '/Volumes/repo/MaNGA/redux/v2_2_0'
#    analysis_path = '/Users/westfall/Work/MaNGA/testing/mpl6/emissionlines/analysis/v2_2_0/test'
#    redux_path='/Users/mac/Desktop/DRP_MPL5/'
    redux_path='/Users/mac/Desktop/DRP_220/'    # v 2_2_0 redux path
    analysis_path='/Users/mac/Volumes/External/MaNGA/DAP_Output/MPL-5/2.0.2/'
    clobber = False
    
    # Read the DRP LOGCUBE file
    print('Reading DRP fits file.')
    drpf = DRPFits(plate, ifu, 'CUBE', read=True, redux_path=redux_path)

    # Calculate the S/N and coordinates
    rdxqa = ReductionAssessment('SNRG', drpf, clobber=clobber, analysis_path=analysis_path)

    # Peform the Voronoi binning to S/N>~10
    binned_spectra = SpatiallyBinnedSpectra('VOR10', drpf, rdxqa, clobber=clobber,
                                            analysis_path=analysis_path)

    # Fit the stellar kinematics
    stellar_continuum = StellarContinuumModel('GAU-MILESHC', binned_spectra, clobber=clobber,
                                              guess_vel=vel, guess_sig=100.,
                                              analysis_path=analysis_path)

    # Get the emission-line moments
    emission_line_moments = EmissionLineMoments('EMOMF', binned_spectra, clobber=clobber,
                                                stellar_continuum=stellar_continuum,
                                                redshift=nsa_redshift,
                                                analysis_path=analysis_path)

    # Get an estimate of the redshift of each bin using the first moment
    # of the H-alpha emission line:
    el_init_redshift = numpy.full(binned_spectra.nbins, nsa_redshift, dtype=float)
    # HARDCODED FOR A SPECIFIC EMISSION-LINE MOMENT DATABASE
    # TODO: Should
    #   - pass the EmissionLineMoments object to EmissionLineModel
    #   - include the channel used for this in EmissionLineModelDef
    halpha_channel = 7
    halpha_mom1_masked = emission_line_moments['ELMMNTS'].data['MASK'][:,halpha_channel] > 0
    # - Use the 1st moment of the H-alpha line
    el_init_redshift[ emission_line_moments['ELMMNTS'].data['BINID_INDEX'] ] \
                = emission_line_moments['ELMMNTS'].data['MOM1'][:,halpha_channel] \
                                / astropy.constants.c.to('km/s').value

    # - For missing bins in the moment measurements and bad H-alpha
    #   moment measurements, use the value for the nearest good bin
    bad_bins = numpy.append(emission_line_moments.missing_bins,
                    emission_line_moments['ELMMNTS'].data['BINID'][halpha_mom1_masked]).astype(int)
    if len(bad_bins) > 0:
        nearest_good_bin_index = binned_spectra.find_nearest_bin(bad_bins, indices=True)
        bad_bin_index = binned_spectra.get_bin_indices(bad_bins)
        el_init_redshift[bad_bin_index] = el_init_redshift[nearest_good_bin_index]

    # Setup the new emission-line fitter
    fitter = XJMCEmissionLineFitter()

    # Setup the new fitting method
    fit_method = EmissionLineModelDef('XJMC',       # Key for the fitting method
                                      0.0,          # Minimum S/N of the binned spectra to include
                                      None,         # Keyword for an artifact mask
                                      None,         # Keyword for an emission-line database
                                      fitter.par,   # Object with fit parameters
                                      fitter,       # Fitting class instance
                                      fitter.fit)   # Fitting function

    # Fit the emission lines
    emission_line_model = EmissionLineModel('XJMC',
                                            binned_spectra,
                                            stellar_continuum=stellar_continuum,
                                            redshift=el_init_redshift, dispersion=100.0,
                                            method_list=fit_method, 
                                            analysis_path=analysis_path, clobber=True)

    # The rest of this is just a single execution of the remaining
    # analysis steps in
    # $MANGADAP_DIR/python/mangadap/survey/manga_dap.py , with some
    # simplifications
#    spectral_indices = SpectralIndices('INDXEN', binned_spectra, redshift=nsa_redshift,
#                                       stellar_continuum=stellar_continuum,
#                                       emission_line_model=emission_line_model,
#                                       analysis_path=analysis_path,
#                                       clobber=False)
    spectral_indices = None

    construct_maps_file(drpf, rdxqa=rdxqa, binned_spectra=binned_spectra,
                        stellar_continuum=stellar_continuum,
                        emission_line_moments=emission_line_moments,
                        emission_line_model=emission_line_model,
                        spectral_indices=spectral_indices, nsa_redshift=nsa_redshift,
                        analysis_path=analysis_path)

    construct_cube_file(drpf, binned_spectra=binned_spectra, stellar_continuum=stellar_continuum,
                        emission_line_model=emission_line_model,
                        analysis_path=analysis_path)

    print('Elapsed time: {0} seconds'.format(time.clock() - t))

