#!/usr/bin/env python3

import time
import warnings
import numpy
import os
import astropy.constants
import pdb

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
#Based on KBW's test_new_binning_scheme.py

class SquareBinning():

    """
    Perform binning of full cube in square apertures
    Length of aperture side is given in arcsec with binsz

    Args:
        binsz (array-like): desired bin size

    """
    def __init__(self, binsz):
        self.binsz = binsz


    def bin_spaxels(self, x, y, par=None):
        _x = numpy.asarray(x)
        if len(_x.shape) != 1:
            raise ValueError('On-sky coordinates must be one-dimensional.')
        nspaxels = _x.size
        if len(y) != nspaxels:
            raise ValueError('Input coordinates are of different lengths.')
        _y = numpy.asarray(y)

        binid = numpy.full(nspaxels, -1, dtype=int)
        binsz = self.binsz
        minx = numpy.min(_x)-0.25
        maxx = numpy.max(_x)+0.25
        miny = numpy.min(_y)-0.25
        maxy = numpy.max(_y)+0.25

        ctbin = 0

        # Start from the center, work out to edge of each quadrant

        # Quadrant 1
        nslicex = numpy.ceil((maxx) / binsz)
        nslicey = numpy.ceil((maxy) / binsz)

        x_lo = numpy.linspace(0.0,(nslicex*binsz),nslicex+1.0)
        y_lo = numpy.linspace(0.0,(nslicey*binsz),nslicey+1.0)

        # Find which spaxels land in each aperture
        for xi in x_lo:
            for yj in y_lo:
                indx = (_x >= xi) & (_x < xi+binsz) & (_y >= yj) & (_y < yj+binsz)
                if any(indx):
                    binid[indx] = ctbin
                    ctbin = ctbin + 1


        # Quadrant 2
        nslicex = numpy.ceil((maxx) / binsz)
        nslicey = numpy.floor((miny) / binsz)

        x_lo = numpy.linspace(0.0, (nslicex * binsz), nslicex + 1.0)
        y_lo = numpy.linspace(0.0, (nslicey * binsz), numpy.abs(nslicey) + 1.0)

        # Find which spaxels land in each aperture
        for xi in x_lo:
            for yj in y_lo:
                indx = (_x >= xi) & (_x < xi + binsz) & (_y <= yj) & (_y > yj - binsz)
                if any(indx):
                    binid[indx] = ctbin
                    ctbin = ctbin + 1

        # Quadrant 3
        nslicex = numpy.floor((minx) / binsz)
        nslicey = numpy.ceil((maxy) / binsz)

        x_lo = numpy.linspace(0.0, (nslicex * binsz), numpy.abs(nslicex) + 1.0)
        y_lo = numpy.linspace(0.0, (nslicey * binsz), numpy.abs(nslicey) + 1.0)

        # Find which spaxels land in each aperture
        for xi in x_lo:
            for yj in y_lo:
                indx = (_x <= xi) & (_x > xi - binsz) & (_y >= yj) & (_y < yj + binsz)
                if any(indx):
                    binid[indx] = ctbin
                    ctbin = ctbin + 1

        # Quadrant 4
        nslicex = numpy.floor((minx) / binsz)
        nslicey = numpy.floor((miny) / binsz)

        x_lo = numpy.linspace(0.0, (nslicex * binsz), numpy.abs(nslicex) + 1.0)
        y_lo = numpy.linspace(0.0, (nslicey * binsz), numpy.abs(nslicey) + 1.0)

        # Find which spaxels land in each aperture
        for xi in x_lo:
            for yj in y_lo:
                indx = (_x <= xi) & (_x > xi - binsz) & (_y <= yj) & (_y > yj - binsz)
                if any(indx):
                    binid[indx] = ctbin
                    ctbin = ctbin + 1


        #pdb.set_trace()
        return binid




#-----------------------------------------------------------------------------
if __name__ == '__main__':
    t = time.clock()

    # Set the plate, ifu, and initial velocity/redshift
    binsz = 2.5        # Desired bin size in arcsec
    plate = 7443
    ifu = 1902
    #vel = 8675.5
    #nsa_redshift = vel/astropy.constants.c.to('km/s').value
    #nsa_redshift = 0.0280161   # for 7443-6102
    #nsa_redshift = 0.020478
    nsa_redshift = 0.0189259
    vel = nsa_redshift * astropy.constants.c.to('km/s').value

    analysis_path = os.getcwd()

    # Read the DRP LOGCUBE file
    drpf = DRPFits(plate, ifu, 'CUBE', read=True)


    # Set up binning function
    sqbin = SquareBinning(binsz)


    # Calculate the S/N and coordinates
    rdxqa = ReductionAssessment('SNRG', drpf, analysis_path=analysis_path)
    x = rdxqa.hdu['SPECTRUM'].data['SKY_COO'][:, 0]
    y = rdxqa.hdu['SPECTRUM'].data['SKY_COO'][:, 1]
    fgdpix = rdxqa.hdu['SPECTRUM'].data['FGOODPIX'] > 0.8

    #binid = sqbin.bin_spaxels(x[fgdpix], y[fgdpix], par=None)
    #pdb.set_trace()

    # Setup the stacking operations
    stackpar = SpectralStackPar('mean',  # Operation for stack
                                False,  # Apply a velocity registration
                                None,  # Velocity offsets for registration
                                'channels',  # Covariance mode and parameters
                                SpectralStack.parse_covariance_parameters('channels', 11),
                                True,  # Propagate the LSF through the stacking
                                True)  # Use pre-pixelized LSF (KHRR added this)
    stacker = SpectralStack()

    # Create a new binning method
    binning_method = SpatiallyBinnedSpectraDef('5x5n',  # Key for binning method
                                               'ODonnell',  # Galactic reddening function to use
                                               3.1,  # Rv for Galactic reddening
                                               0.0,  # Minimum S/N to include
                                               None,  # Object with binning parameters
                                               None,  # Binning class instance
                                               sqbin.bin_spaxels,  # Binning function
                                               stackpar,  # Object with stacking parameters
                                               stacker,  # Stacking class instance
                                               stacker.stack_DRPFits,  # Stacking function
                                               'spaxel',  # Type of LSF characterization to use
                                               True)  # Use pre-pixelized LSF (KHRR added this)

    # Bin the spectra using the new binning method
    binned_spectra = SpatiallyBinnedSpectra('5x5n',  # Key for binning method
                                            drpf,  # DRP data to bin
                                            rdxqa,  # Cube coordinates and S/N assessments
                                            method_list=binning_method,  # Methods to select from
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
    #emission_line_moments = None

    emission_line_model = EmissionLineModel('EFITM', binned_spectra,
                                            stellar_continuum=stellar_continuum,
                                            redshift=nsa_redshift, dispersion=100.0,
                                            analysis_path=analysis_path)
    #emission_line_model = None

    #    spectral_indices = SpectralIndices('INDXEN', binned_spectra, redshift=nsa_redshift,
    #                                       stellar_continuum=stellar_continuum,
    #                                       emission_line_model=emission_line_model,
    #                                       analysis_path=analysis_path)
    spectral_indices = None

    construct_maps_file(drpf, rdxqa=rdxqa, binned_spectra=binned_spectra,
                        stellar_continuum=stellar_continuum,
                        emission_line_moments=emission_line_moments,
                        emission_line_model=emission_line_model,
                        spectral_indices=spectral_indices, nsa_redshift=nsa_redshift,
                        analysis_path=analysis_path)

    construct_cube_file(drpf, binned_spectra=binned_spectra, stellar_continuum=stellar_continuum,
                        emission_line_model=emission_line_model, analysis_path=analysis_path)

    print('Elapsed time: {0} seconds'.format(time.clock() - t))



    #pdb.set_trace()