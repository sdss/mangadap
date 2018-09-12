#!/usr/bin/env python3

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
        if len(r) != 1 and len(r) != self.n:
            raise ValueError('Radii must be common to all apertures or unique to each aperture.')
        self.r = numpy.full(self.n, r, dtype=float) if len(r) == 1 else numpy.asarray(r)

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
    t = time.clock()

    # Set the plate, ifu, and initial velocity/redshift
    plate = 7443
    ifu = 6102
    #vel = 8675.5
    #nsa_redshift = vel/astropy.constants.c.to('km/s').value
    nsa_redshift = 0.0280161
    vel = nsa_redshift * astropy.constants.c.to('km/s').value

    # Read the DRP LOGCUBE file
    drpf = DRPFits(plate, ifu, 'CUBE', read=True)

    # Calculate the S/N and coordinates
    rdxqa = ReductionAssessment('SNRG', drpf)

    # Setup the aperture binning class
    ax = numpy.array([0.0, 3.0, 6.0])
    ay = numpy.array([0.0, 0.0, 0.0])
    apbin = ApertureBinning(ax, ay, numpy.array([2.5]))

    # Setup the stacking operations
    stackpar = SpectralStackPar('mean',         # Operation for stack
                                False,          # Apply a velocity registration
                                None,           # Velocity offsets for registration
                                'channels',     # Covariance mode and parameters
                                SpectralStack.parse_covariance_parameters('channels', 11),
                                True,           # Propagate the LSF through the stacking
                                True)  # Use pre-pixelized LSF (KHRR added this)
    stacker = SpectralStack()

    # Create a new binning method
    binning_method = SpatiallyBinnedSpectraDef('Aperture',      # Key for binning method
                                               'ODonnell',      # Galactic reddening function to use
                                               3.1,             # Rv for Galactic reddening
                                               0.0,             # Minimum S/N to include
                                               None,            # Object with binning parameters
                                               None,            # Binning class instance
                                               apbin.bin_spaxels,   # Binning function
                                               stackpar,        # Object with stacking parameters
                                               stacker,         # Stacking class instance
                                               stacker.stack_DRPFits,   # Stacking function
                                               'spaxel',    # Type of LSF characterization to use
                                               True)  # Use pre-pixelized LSF (KHRR added this)

    # Bin the spectra using the new binning method
    binned_spectra = SpatiallyBinnedSpectra('Aperture',     # Key for binning method
                                            drpf,           # DRP data to bin
                                            rdxqa,          # Cube coordinates and S/N assessments
                                            method_list=binning_method) # Methods to select from

    # The rest of this is just a single execution of the remaining
    # analysis steps in
    # $MANGADAP_DIR/python/mangadap/survey/manga_dap.py , with some
    # simplifications

    stellar_continuum = StellarContinuumModel('GAU-MILESHC', binned_spectra, guess_vel=vel,
                                              guess_sig=100.)

    emission_line_moments = EmissionLineMoments('EMOMF', binned_spectra,
                                                stellar_continuum=stellar_continuum,
                                                redshift=nsa_redshift)

    emission_line_model = EmissionLineModel('EFITF', binned_spectra,
                                            stellar_continuum=stellar_continuum,
                                            redshift=nsa_redshift, dispersion=100.0)
        
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

