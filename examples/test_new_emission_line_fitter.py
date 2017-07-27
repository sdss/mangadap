#!/usr/bin/env python3

import time
import warnings
import numpy
import astropy.constants

from mangadap.drpfits import DRPFits
from mangadap.proc.reductionassessments import ReductionAssessment
from mangadap.proc.spatiallybinnedspectra import SpatiallyBinnedSpectra
from mangadap.proc.stellarcontinuummodel import StellarContinuumModel
from mangadap.proc.ppxffit import PPXFFit
from mangadap.util.instrument import spectrum_velocity_scale
from mangadap.util.fitsutil import DAPFitsUtil
from mangadap.util.fileio import init_record_array
from mangadap.proc.emissionlinemoments import EmissionLineMoments
from mangadap.proc.spectralfitting import EmissioneLineFit
from mangadap.proc.emissionlinemodel import EmissionLineModelDef, EmissionLineModel
from mangadap.par.emissionlinedb import EmissionLineDB
from mangadap.proc.spectralindices import SpectralIndices
from mangadap.dapfits import construct_maps_file, construct_cube_file

from mangadap.contrib.xjmc import first_second_iteration

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
                    'degree': 8,                # Additive polynomial order
                    'mdegree': 0 }              # Multiplicative polynomial order
        EmissionLineFit.__init__(self, 'XJMC', None, par=par)


    def fit(self, binned_spectra, par=None, loggers=None, quiet=False):
        if par is not None:
            self.par = par

        # Check the parameter keys
        required_keys = [ 'guess_redshift', 'stellar_continuum', 'emission_lines', 'degree',
                          'mdegree' ]
        if numpy.any([ reqk not in self.par.keys() for reqk in required_keys ]):
            raise ValueError('Parameter dictionary does not have all the required keys.')

        # Wavelengths are in vacuum
        wave = binned_spectra['WAVE'].data.copy()
        # Velocity step per pixel
        velscale = spectrum_velocity_scale(wave)
        # Flux and noise masked arrays; shape is (Nspaxels,Nwave) where
        # Nspaxels is Nx*Ny
        flux = binned_spectra.drpf.copy_to_masked_array(flag=['DONOTUSE', 'FORESTAR'])
        noise = numpy.ma.power(binned_spectra.drpf.copy_to_masked_array(ext='IVAR',
                                    flag=['DONOTUSE', 'FORESTAR']), -0.5)
        # Spaxel coordinates; shape is (nspaxels,)
        x = binned_spectra.rdxqa['SPECTRUM'].data['SKY_COO'][:,0]
        y = binned_spectra.rdxqa['SPECTRUM'].data['SKY_COO'][:,1]
        # Binned flux and binned noise masked arrays; shape is (nbins,nwave)
        flux_binned = binned_spectra.copy_to_masked_array(flag=binned_spectra.do_not_fit_flags())
        noise_binned = numpy.ma.power(binned_spectra.copy_to_masked_array(ext='IVAR',
                                                    flag=binned_spectra.do_not_fit_flags()) , -0.5)
        # Bin coordinates; shape is (nbins,)
        x_binned = binned_spectra['BINS'].data['SKY_COO'][:,0]
        y_binned = binned_spectra['BINS'].data['SKY_COO'][:,1]
        # Set initial guesses for the velocity and velocity dispersion
        if self.par['guess_redshift'] is not None:
            # Use guess_redshift if provided
            vel = self.par['guess_redshift'] * astropy.constants.c.to('km/s').value
            # And set default velocity dispersion to 100 km/s
            sig = numpy.full(vel.size, 100, dtype=float)
        elif self.par['stellar_continuum'] is not None:
            # Otherwise use the stellar-continuum result
            vel, sig = self.par['stellar_continuum'].matched_guess_kinematics(binned_spectra,
                                                                              cz=True)
        else:
            # TODO: Set default guess kinematics
            vel, sig = None, None

        # Get the stellar templates;
        # shape is (Ntemplates, Nwave_templates)
        if self.par['stellar_continuum'] is not None:
            stars_templates = self.par['stellar_continuum'].method['fitpar']['template_library']
            stars_templates_wave = stars_templates['WAVE'].data.copy()
            stars_templates = stars_templates['FLUX'].data.copy()
            velscale_ratio = self.par['stellar_continuum'].method['fitpar']['velscale_ratio']
            dv = -PPXFFit.ppxf_tpl_obj_voff(stars_templates_wave, wave, velscale,
                                            velscale_ratio=velscale_ratio)
        else:
            # TODO: Default construction of stellar templates
            stars_templates = None
            velscale_ratio = None
            dv = None

        # TODO: Construct gas templates
        gas_templates = None
        gas_names = None

        # TODO: Default polynomial orders
        degree = 8 if self.par['degree'] is None else self.par['degree']
        mdegree = 0 if self.par['mdegree'] is None else self.par['mdegree']
        
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
        #       flux[:,0].reshape((numpy.sqrt(Nspaxels).astype(int),)*2).shape
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
        model_flux, model_eml_flux, model_mask, model_binid, eml_flux, eml_fluxerr, \
                eml_kin, eml_kinerr, eml_sigmacorr \
                        = first_second_iteration(wave, flux, noise, flux_binned, noise_binned,
                                                 velscale, velscale_ratio, dv, vel, sig,
                                                 stars_templates, gas_templates, gas_names,
                                                 degree, mdegree, x, y, x_binned, y_binned)
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

        # Calculate the "emission-line baseline" as the difference
        # between the stellar continuum model determined for the
        # kinematics and the one determined by the optimized
        # stellar-continuum + emission-line fit:
        if self.par['stellar_continuum'] is not None:
            # Construct the full 3D cube for the stellar continuum
            # models
            sc_model_flux, sc_model_mask \
                    = DAPFitsUtil.reconstruct_cube(binned_spectra.drpf.shape,
                                                   self.par['stellar_continuum']['BINID'].data,
                                                   [ self.par['stellar_continuum']['FLUX'].data,
                                                     self.par['stellar_continuum']['MASK'].data ])
            # Set any masked pixels to 0
            sc_model_flux[sc_model_mask>0] = 0.0

            # Construct the full 3D cube of the new stellar continuum
            # from the combined stellar-continuum + emission-line fit
            el_continuum = DAPFitsUtil.reconstruct_cube(binned_spectra.drpf.shape, model_binid,
                                                        model_flux - model_eml_flux)
            # Get the difference, restructure it to match the shape
            # of the emission-line models, and zero any masked pixels
            model_eml_base = (el_model_flux - sc_model_flux).reshape(-1,wave.size)[model_srt,:]
            if model_mask is not None:
                model_eml_base[model_mask>0] = 0.0
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
    plate = 7495
    ifu = 12704
    vel = 8675.5
    nsa_redshift = vel/astropy.constants.c.to('km/s').value

    # Read the DRP LOGCUBE file
    drpf = DRPFits(plate, ifu, 'CUBE', read=True)

    # Calculate the S/N and coordinates
    rdxqa = ReductionAssessment('SNRG', drpf)

    # Peform the Voronoi binning to S/N>~10
    binned_spectra = SpatiallyBinnedSpectra('VOR10', drpf, rdxqa)

    # Fit the stellar kinematics
    stellar_continuum = StellarContinuumModel('GAU-MILESHC', binned_spectra, guess_vel=vel,
                                              guess_sig=100.)

    # Get the emission-line moments
    emission_line_moments = EmissionLineMoments('EMOMF', binned_spectra,
                                                stellar_continuum=stellar_continuum,
                                                redshift=nsa_redshift)

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
    # TODO: Include an emission-line database with the lines that
    # matches the lines fit by XJMCEmissionLineFitter
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
                                            method_list=fit_method)

    # The rest of this is just a single execution of the remaining
    # analysis steps in
    # $MANGADAP_DIR/python/mangadap/survey/manga_dap.py , with some
    # simplifications
    spectral_indices = SpectralIndices('INDXEN', binned_spectra, redshift=nsa_redshift,
                                       stellar_continuum=stellar_continuum,
                                       emission_line_model=emission_line_model)

    construct_maps_file(drpf, rdxqa=rdxqa, binned_spectra=binned_spectra,
                        stellar_continuum=stellar_continuum,
                        emission_line_moments=emission_line_moments,
                        emission_line_model=emission_line_model,
                        spectral_indices=spectral_indices, nsa_redshift=nsa_redshift)

    construct_cube_file(drpf, binned_spectra=binned_spectra, stellar_continuum=stellar_continuum,
                        emission_line_model=emission_line_model)

    print('Elapsed time: {0} seconds'.format(time.clock() - t))

