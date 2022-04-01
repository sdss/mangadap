# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Provides the main wrapper function for the MaNGA DAP.

----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import logging
import resource
import time
import os
import warnings

from IPython import embed

import numpy

import astropy.constants

from mangadap import __version__

from ..util.log import init_DAP_logging, module_logging, log_output
from ..datacube import DataCube
from ..config.analysisplan import AnalysisPlan
from ..proc.reductionassessments import ReductionAssessment
from ..proc.spatiallybinnedspectra import SpatiallyBinnedSpectra
from ..proc.stellarcontinuummodel import StellarContinuumModel
from ..proc.emissionlinemoments import EmissionLineMoments
from ..proc.emissionlinemodel import EmissionLineModel
from ..proc.spectralindices import SpectralIndices
from ..dapfits import construct_maps_file, construct_cube_file


def manga_dap(cube, plan, dbg=False, log=None, verbose=0):
    r"""
    Main wrapper function for the MaNGA DAP.

    This function is designed to be called once per datacube. The
    :class:`mangadap.config.analysisplan.AnalysisPlan` instance sets
    the types of analyses to perform on this observation. Each
    analysis plan results in a MAPS and model LOGCUBE file.

    The procedure is as follows.  For each plan in ``plan``:

        - Determine basic assessments of the data, including the S/N
          and spaxel coordinates. See
          :class:`mangadap.proc.reductionassessments.ReductionAssessment`.
        - Bin the spectra. See
          :class:`mangadap.proc.spatiallybinnedspectra.SpatiallyBinnedSpectra`.
        - Fit the stellar continuum for stellar kinematics. See
          :class:`mangadap.proc.stellarcontinuummodel.StellarContinuumModel`.
        - Measure the emission-line moments. See
          :class:`mangadap.proc.emissionlinemoments.EmissionLineMoments`.
        - Fit parameterized line profiles to the emission lines. See
          :class:`mangadap.proc.emissionlinemodel.EmissionLineModel`.
        - Subtract the fitted emission-line models and measure the
          spectral indices. See
          :class:`mangadap.proc.spectralindices.SpectralIndices`.
        - Construct the primary output files based on the plan
          results. See :func:`mangadap.dapfits.construct_maps_file`
          and :func:`mangadap.dapfits.construct_cube_file`.

    Verbose levels (still under development):
        0. Nothing but errors.
        1. Basic status updates. E.g., start and end of each block,
           minor progress within blocks.
        2. Block-level and sub-function updates.
        3. Same as above, with figures.

    Args:
        cube (:class:`mangadap.datacube.datacube.DataCube`):
            Datacube to analyze.
        plan (:class:`mangadap.config.analysisplan.AnalysisPlan`):
            Object that sets the analysis plans to implement and the output
            paths for DAP files.
        dbg (:obj:`bool`, optional):
            Flag to run the DAP in debug mode, see above; default is
            False. Limited use; still under development.
        log (:obj:`str`, optional):
            File name for log output, see above; no log file is
            produced by default.
        verbose (:obj:`int`, optional):
            Verbosity level, see above; default is 0.

    Returns:
        :obj:`int`: Status flag (under development; currently always
        0).
    """
    # Check input
    #  - Check that input cube type
    if not isinstance(cube, DataCube):
        raise TypeError('Input objects must be a subclass of DataCube.')
    # Check that the datacube can be analyzed by the DAP
    if not cube.can_analyze:
        raise ValueError('Cube metadata is insufficient for analysis.  See DataCube.can_analyze.')
    # Check that input plan type
    if not isinstance(plan, AnalysisPlan):
        raise TypeError('plan must be of type AnalysisPlan')

    # Initialize the logging objects
    init_DAP_logging(log)#, simple_warnings=False)

    # Start log
    loggers = module_logging(__name__, verbose)

    # Check the root output path exists, and create it if not.
    if not plan.analysis_path.exists():
        plan.analysis_path.mkdir(parents=True)

    log_output(loggers, 1, logging.INFO, '-'*50)
    log_output(loggers, 1, logging.INFO, '{0:^50}'.format('MaNGA Data Analysis Pipeline'))
    log_output(loggers, 1, logging.INFO, '-'*50)
    log_output(loggers, 1, logging.INFO, '     VERSION: {0}'.format(__version__))
    log_output(loggers, 1, logging.INFO, '       START: {0}'.format(
                                    time.strftime('%a %d %b %Y %H:%M:%S',time.localtime())))
    t = time.perf_counter()
    log_output(loggers, 1, logging.INFO, '  INPUT ROOT: {0}'.format(cube.directory_path))
    log_output(loggers, 1, logging.INFO, '   CUBE FILE: {0}'.format(cube.file_name))
    log_output(loggers, 1, logging.INFO, '     N PLANS: {0}'.format(plan.nplans))
    log_output(loggers, 1, logging.INFO, ' OUTPUT ROOT: {0}'.format(plan.analysis_path))
    log_output(loggers, 1, logging.INFO, '-'*50)

    status = 0

    # Iterate over plans:
    for i, key in enumerate(plan.keys()):
        # Construct some directories to ensure all the reference files
        # are placed in the same path.  The references files placed in
        # the common directory are:
        #   - ReductionAssessments
        #   - SpatiallyBinnedSpectra
        #   - StellarContinuumModel
        # These reference files are symlinked to each DAPTYPE reference
        # directory.  The reference files specific to each DAPTYPE are
        #   - EmissionLineMoments (although this could also go in the
        #   common directory)
        #   - EmissionLineModel
        #   - SpectralIndices
        common_dir = plan.common_path()
        method_dir = plan.method_path(plan_index=i)
        method_ref_dir = plan.method_path(plan_index=i, ref=True)

        log_output(loggers, 1, logging.INFO, '-'*50)
        log_output(loggers, 1, logging.INFO, f'{f"Plan {i+1}: {key}":^50}')
        log_output(loggers, 1, logging.INFO, '-'*50)
        log_output(loggers, 1, logging.INFO, f' METHOD KEY: {plan[key]["key"]}')
        log_output(loggers, 1, logging.INFO, f'COMMON PATH: {common_dir}')
        log_output(loggers, 1, logging.INFO, f'METHOD PATH: {method_dir}')

        #---------------------------------------------------------------
        # S/N Assessments: placed in the common/ directory
        #---------------------------------------------------------------
        rdxqa = None if plan.rdxqa[key] is None else \
                    ReductionAssessment(plan.rdxqa[key], cube, pa=cube.meta['pa'],
                                        ell=cube.meta['ell'], reff=cube.meta['reff'],
                                        output_path=common_dir, symlink_dir=method_ref_dir,
                                        loggers=loggers)

        #---------------------------------------------------------------
        # Spatial Binning: placed in the common/ directory
        #---------------------------------------------------------------
        binned_spectra = None if plan.binning[key] is None else \
                    SpatiallyBinnedSpectra(plan.binning[key], cube, rdxqa, reff=cube.meta['reff'],
                                           ebv=cube.meta['ebv'], output_path=common_dir,
                                           symlink_dir=method_ref_dir, loggers=loggers)

        #---------------------------------------------------------------
        # Stellar Continuum Fit: placed in the common/ directory
        #---------------------------------------------------------------
        # TODO: Allow tpl_hardcopy to be an input parameter?
        stellar_continuum = None if plan.continuum[key] is None else \
                    StellarContinuumModel(plan.continuum[key], binned_spectra,
                                          guess_vel=cube.meta['vel'], guess_sig=cube.meta['vdisp'],
                                          output_path=common_dir, symlink_dir=method_ref_dir,
                                          loggers=loggers)

        #---------------------------------------------------------------
        # For both the emission-line measurements (moments and model):
        # To use the NSA redshift for all moments:
        #   - set redshift=nsa_redshift
        # To use the stellar kinematics to set the redshift:
        #   - set redshift=None
        #---------------------------------------------------------------
        # Emission-line Moment measurements: placed in the DAPTYPE/ref/
        # directory
        #---------------------------------------------------------------
#        warnings.filterwarnings('error', message='Warning: converting a masked element to nan.')
        # TODO: enable artifact_path and bandpass_path?
        emission_line_moments = None if plan.elmom[key] is None else \
                    EmissionLineMoments(plan.elmom[key], binned_spectra,
                                        stellar_continuum=stellar_continuum,
                                        redshift=cube.meta['z'], output_path=method_ref_dir,
                                        loggers=loggers)

        #---------------------------------------------------------------
        # Emission-line Fit: placed in the DAPTYPE/ref/ directory
        #---------------------------------------------------------------
        # To use a fixed velocity dispersion for the initial guess, e.g.:
        #   - set dispersion=100.0
        # To use the stellar kinematics as the initial guess:
        #   - set dispersion=None
        # The provided redshift is only the initial guess for the
        # emission-line model (it's FIXED for the moment measurements
        # above).
        emission_line_model = None if plan.elfit[key] is None else \
                   EmissionLineModel(plan.elfit[key], binned_spectra,
                                     stellar_continuum=stellar_continuum,
                                     emission_line_moments=emission_line_moments, dispersion=100.0,
                                     minimum_error=numpy.finfo(numpy.float32).eps,
                                     output_path=method_ref_dir, loggers=loggers)

        #---------------------------------------------------------------
        # If requested by the emission-line moments method, remeasure
        # the moments after the emission-line modeling.  This will
        # produce a new reference file that will have a different name
        # than the one produced above.  This is placed in the
        # DAPTYPE/ref/ directory.
        #---------------------------------------------------------------
        if emission_line_moments is not None \
                and emission_line_moments.database['redo_postmodeling']:
            emission_line_moments.measure(binned_spectra, stellar_continuum=stellar_continuum,
                                          emission_line_model=emission_line_model,
                                          output_path=method_ref_dir, loggers=loggers)

        #---------------------------------------------------------------
        # Spectral-Index Measurements: placed in the DAPTYPE/ref/
        # directory
        #---------------------------------------------------------------
        spectral_indices = None if plan.sindx[key] is None else \
                    SpectralIndices(plan.sindx[key], binned_spectra, redshift=cube.meta['z'],
                                    stellar_continuum=stellar_continuum,
                                    emission_line_model=emission_line_model,
                                    output_path=method_ref_dir, loggers=loggers)

        #-------------------------------------------------------------------
        # Construct the main output file(s)
        #-------------------------------------------------------------------
        construct_maps_file(cube, method=plan[key]['key'], rdxqa=rdxqa,
                            binned_spectra=binned_spectra, stellar_continuum=stellar_continuum,
                            emission_line_moments=emission_line_moments,
                            emission_line_model=emission_line_model,
                            spectral_indices=spectral_indices, redshift=cube.meta['z'],
                            output_path=method_dir, overwrite=True, loggers=loggers,
                            single_precision=True)

        construct_cube_file(cube, method=plan[key]['key'], rdxqa=rdxqa,
                            binned_spectra=binned_spectra, stellar_continuum=stellar_continuum,
                            emission_line_moments=emission_line_moments,
                            emission_line_model=emission_line_model,
                            spectral_indices=spectral_indices, output_path=method_dir,
                            overwrite=True, loggers=loggers, single_precision=True)

        # Mark for next garbage collection.
        # TODO: Not sure this is useful because these variables may
        # just get immediately defined again in the next iteration.
        del spectral_indices
        del emission_line_model
        del emission_line_moments
        del stellar_continuum
        del binned_spectra
        del rdxqa

    #-------------------------------------------------------------------
    # End log
    #-------------------------------------------------------------------
    log_output(loggers, 1, logging.INFO, '-'*50)
    log_output(loggers, 1, logging.INFO, 'EXECUTION SUMMARY')
    log_output(loggers, 1, logging.INFO, '-'*50)
    log_output(loggers, 1, logging.INFO,
               f'    FINISH: {time.strftime("%a %d %b %Y %H:%M:%S",time.localtime())}')
    if time.perf_counter() - t < 60:
        tstr = f'  DURATION: {time.perf_counter() - t:.5f} sec'
    elif time.perf_counter() - t < 3600:
        tstr = f'  DURATION: {(time.perf_counter() - t)/60.:.5f} min'
    else:
        tstr = f'  DURATION: {(time.perf_counter() - t)/3600.:.5f} hr'
    log_output(loggers, 1, logging.INFO, tstr)

    mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss \
            + resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    if mem < 1024:
        mstr = f'MAX MEMORY: {mem:.4f} bytes'
    elif mem < 1024*1024:
        mstr = f'MAX MEMORY: {mem/1024:.4f} Kbytes'
    elif mem < 1024*1024*1024:
        mstr = f'MAX MEMORY: {mem/1024/1024:.4f} Mbytes'
    else:
        mstr = f'MAX MEMORY: {mem/1024/1024/1024:.4f} Gbytes'
    log_output(loggers, 1, logging.INFO, mstr)
    log_output(loggers, 1, logging.INFO, '-'*50)
    log_output(loggers, 1, logging.INFO, '-'*50)

    return status


