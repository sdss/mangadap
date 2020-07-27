#!/usr/bin/env python3

import os
import time
import warnings
import numpy
import astropy.constants

from mangadap.drpfits import DRPFits
from mangadap.proc.reductionassessments import ReductionAssessment
from mangadap.proc.spectralstack import SpectralStackPar, SpectralStack
from mangadap.proc.spatiallybinnedspectra import SpatiallyBinnedSpectra, SpatiallyBinnedSpectraDef
from mangadap.proc.stellarcontinuummodel import StellarContinuumModel
from mangadap.proc.emissionlinemoments import EmissionLineMoments
from mangadap.proc.emissionlinemodel import EmissionLineModel
from mangadap.proc.spectralindices import SpectralIndices
from mangadap.dapfits import construct_maps_file, construct_cube_file

#-----------------------------------------------------------------------------

class ApertureBinning():
    """
    Perform aperture binning

    Args:
        x (array-like): List of on-sky x coordinates for apertures
        y (array-like): List of on-sky y coordinates for apertures
        r (array-like): Single or list of radii of the apertures

    Attributes:
        n (int): Number of apertures
        x (numpy.ndarray): On-sky x coordinates for apertures
        y (numpy.ndarray): On-sky y coordinates for apertures
        r (numpy.ndarray): Aperture radii
    """
    def __init__(self, x, y, r):
        self.x = numpy.asarray(x)
        if len(self.x.shape) != 1:
            raise ValueError('On-sky coordinates must be one-dimensional.')
        self.n = x.size
        if len(y) != self.n:
            raise ValueError('Input coordinates are of different lengths.')
        self.y = numpy.asarray(y)
        _r = numpy.atleast_1d(r)
        if len(_r) != 1 and len(_r) != self.n:
            raise ValueError('Radii must be common to all apertures or unique to each aperture.')
        self.r = numpy.full(self.n, r, dtype=float) if len(_r) == 1 else _r

    def bin_spaxels(self, x, y, par=None):
        _x = numpy.asarray(x)
        if len(_x.shape) != 1:
            raise ValueError('On-sky coordinates must be one-dimensional.')
        nspaxels = _x.size
        if len(y) != nspaxels:
            raise ValueError('Input coordinates are of different lengths.')
        _y = numpy.asarray(y)

        # Find which spaxels land in each aperture
        indx = numpy.square(_x[:,None]-self.x[None,:]) + numpy.square(_y[:,None]-self.y[None,:]) \
                    < numpy.square(self.r[None,:])
        if numpy.any(numpy.sum(indx, axis=1) > 1):
            warnings.warn('Spaxels found in multiple apertures!')

        # Return the aperture index that each spaxel is within,
        # isolating only one aperture per spaxel; spaxels not in any
        # aperture have a bin ID of -1
        binid = numpy.full((nspaxels, self.n), -1, dtype=int)
        binid[indx] = numpy.array([numpy.arange(self.n)]*nspaxels)[indx]
        return numpy.amax(binid, axis=1)

#-----------------------------------------------------------------------------
if __name__ == '__main__':
    t = time.perf_counter()

    # Set the plate, ifu, and initial velocity/redshift
    plate = 7443
    ifu = 6102
    #vel = 8675.5
    #nsa_redshift = vel/astropy.constants.c.to('km/s').value
    nsa_redshift = 0.0280161
    vel = nsa_redshift * astropy.constants.c.to('km/s').value

    analysis_path = os.path.join(os.getcwd(), 'test_binning_output')

    # Read the DRP LOGCUBE file
    drpf = DRPFits(plate, ifu, 'CUBE', read=True)

    # Calculate the S/N and coordinates
    rdxqa = ReductionAssessment('SNRG', drpf, analysis_path=analysis_path)

    # Setup the aperture binning class
    ax = numpy.array([0.0, 3.0, 6.0])
    ay = numpy.array([0.0, 0.0, 0.0])
    apbin = ApertureBinning(ax, ay, 1.25)

    # Setup the stacking operations
    stackpar = SpectralStackPar(# Operation for stack
                                operation='mean',
                                # Apply a velocity registration
                                vel_register=False,
                                # Velocity offsets for registration
                                vel_offsets=None,
                                # Covariance mode and parameters
                                covar_mode='channels',
                                covar_par=SpectralStack.parse_covariance_parameters('channels', 11),
                                # Propagate the LSF through the stacking
                                stack_sres=True,
                                # Use pre-pixelized LSF (KHRR added this)
                                prepixel_sres=True)
    stacker = SpectralStack()

    # Create a new binning method
    binning_method = SpatiallyBinnedSpectraDef(# Key for binning method
                                               key='Aperture',
                                               # Galactic reddening function to use
                                               galactic_reddening='ODonnell',
                                               # Rv for Galactic reddening
                                               galactic_rv=3.1,             
                                               # Minimum S/N to include
                                               minimum_snr=0.0,
                                               # Binning function
                                               binfunc=apbin.bin_spaxels,
                                               # Object with stacking parameters
                                               stackpar=stackpar,
                                               # Stacking class instance
                                               stackclass=stacker,
                                               # Stacking function
                                               stackfunc=stacker.stack_DRPFits,
                                               # Type of LSF characterization to use
                                               spec_res='spaxel',
                                               # Use pre-pixelized LSF (KHRR added this)
                                               prepixel_sres=True)

    # Bin the spectra using the new binning method
    binned_spectra = SpatiallyBinnedSpectra('Aperture',     # Key for binning method
                                            drpf,           # DRP data to bin
                                            rdxqa,          # Cube coordinates and S/N assessments
                                            method_list=binning_method, # Methods to select from
                                            analysis_path=analysis_path)

    # The rest of this is just a single execution of the remaining
    # analysis steps in
    # $MANGADAP_DIR/python/mangadap/survey/manga_dap.py , with some
    # simplifications

    stellar_continuum = StellarContinuumModel('GAU-MILESHC', binned_spectra, guess_vel=vel,
                                              guess_sig=100., analysis_path=analysis_path)

    emission_line_moments = EmissionLineMoments('EMOMM', binned_spectra,
                                                stellar_continuum=stellar_continuum,
                                                redshift=nsa_redshift, analysis_path=analysis_path)

    emission_line_model = EmissionLineModel('EFITM', binned_spectra,
                                            stellar_continuum=stellar_continuum,
                                            redshift=nsa_redshift, dispersion=100.0,
                                            analysis_path=analysis_path)
        
    spectral_indices = SpectralIndices('INDXEN', binned_spectra, redshift=nsa_redshift,
                                       stellar_continuum=stellar_continuum,
                                       emission_line_model=emission_line_model,
                                       analysis_path=analysis_path)

    construct_maps_file(drpf, rdxqa=rdxqa, binned_spectra=binned_spectra,
                        stellar_continuum=stellar_continuum,
                        emission_line_moments=emission_line_moments,
                        emission_line_model=emission_line_model,
                        spectral_indices=spectral_indices, nsa_redshift=nsa_redshift,
                        analysis_path=analysis_path)

    construct_cube_file(drpf, binned_spectra=binned_spectra, stellar_continuum=stellar_continuum,
                        emission_line_model=emission_line_model, analysis_path=analysis_path)

    print('Elapsed time: {0} seconds'.format(time.perf_counter() - t))

