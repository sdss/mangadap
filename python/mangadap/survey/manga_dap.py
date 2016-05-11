# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Provides the main wrapper function for the MaNGA DAP.

*License*:
    Copyright (c) 2015, SDSS-IV/MaNGA Pipeline Group
        Licensed under BSD 3-clause license - see LICENSE.rst

*Source location*:
    $MANGADAP_DIR/python/mangadap/survey/manga_dap.py

*Imports and python version compliance*:
    ::

        from __future__ import division
        from __future__ import print_function
        from __future__ import absolute_import
        from __future__ import unicode_literals
    
        import sys
        if sys.version > '3':
            long = int

        import logging
        import resource
        import time
        import os
        from ..drpfits import DRPFits
        from ..util.log import init_DAP_logging, module_logging, log_output
        from ..proc.reductionassessments import ReductionAssessment

*Usage*:
    Add usage examples

*Revision history*:
    | **21 Mar 2016**: Original implementation by K. Westfall (KBW): started.

"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import sys
if sys.version > '3':
    long = int

import numpy
import logging
import resource
import time
import os

from ..util.log import init_DAP_logging, module_logging, log_output
from ..drpfits import DRPFits
from ..proc.reductionassessments import ReductionAssessment
from ..proc.spatiallybinnedspectra import SpatiallyBinnedSpectra
from ..proc.stellarcontinuummodel import StellarContinuumModel
from ..proc.emissionlinemoments import EmissionLineMoments
from ..proc.emissionlinemodel import EmissionLineModel
from ..proc.spectralindices import SpectralIndices
from ..dapmaps import construct_maps_file
#from ..proc.templatelibrary import TemplateLibrary

from ..util.covariance import Covariance
from matplotlib import pyplot

__author__ = 'Kyle B. Westfall'
__email__ = 'kbwestfall@gmail.com'
__credits__ = [ 'Kyle Westfall', 'Lodovico Coccato', 'Brett Andrews' ]
__copyright__ = 'Copyright 2016, SDSS-IV/MaNGA Pipeline Group'
__license__ = 'BSD3'
__version__ = 2.0
__status__ = 'Development'

def manga_dap(obs, plan, dbg=False, log=None, verbose=0, redux_path=None, dapsrc=None,
              directory_path=None, analysis_path=None):
    """
    Main wrapper function for the MaNGA DAP.

    This should be called once per ${plate}-${ifudesign} (set by obs),
    ${binmode}_${tpl}_${contmode} (set by plan) combination.

    The procedure is as follows:

        - Read the DRP fits file specified by `obs`
        - For each plan in `plan`:

            - Construct the S/N in the DRP cube (plan dependent because
              wavelength range for S/N can change)
            - Bin the spectra
            - Fit the stellar continuum for stellar kinematics
            - Measure the emission-line moments
            - Fit parameterized line profiles to the emission lines
            - Subtract the fitted emission-line models and measure the
              spectral indices

        - Construct the primary output file based on the plan results.

    Verbose levels:

        0. Nothing but errors.
        1. Basic status updates. E.g., start and end of each block,
           minor progress within blocks.
        2. Block-level and sub-function updates.
        3. Same as above, with figures.

    Args:
        obs (:class:`mangadap.par.obsinput.ObsInputPar`): Object with the input
            parameters.
        plan (:class:`mangadap.par.analysisplan.AnalysisPlan`): Object with the
            analysis plan to implement.
        dbg (bool) : (Optional) Flag to run the DAP in debug mode, see
            above; default is False.
        log (str) : (Optional) File name for log output, see above; no
            log file is produced by default.
        verbose (int) : (Optional) Verbosity level, see above; default
            is 0.
        redux_path (str) : (Optional) Top-level directory with the DRP
            products; default is defined by
            :func:`mangadap.config.defaults.default_redux_path`.
        directory_path (str) : (Optional) Direct path to directory
            containing the DRP output file; default is defined by
            :func:`mangadap.config.defaults.default_drp_directory_path`

        analysis_path (str) : (Optional) Top-level directory for the DAP
            output data; default is defined by
            :func:`mangadap.config.defaults.default_analysis_path`.

    Returns:
        int: Status flag
    """

    init_DAP_logging(log, simple_warnings=False)

    # Start log
    loggers = module_logging(__name__, verbose)

    log_output(loggers, 1, logging.INFO, '   DAPVERS: {0}'.format(__version__))
    log_output(loggers, 1, logging.INFO, '     START: {0}'.format(
                                    time.strftime("%a %d %b %Y %H:%M:%S",time.localtime())))
    t = time.clock()
    log_output(loggers, 1, logging.INFO, '     PLATE: {0}'.format(obs['plate']))
    log_output(loggers, 1, logging.INFO, ' IFUDESIGN: {0}'.format(obs['ifudesign']))
    log_output(loggers, 1, logging.INFO, '   3D MODE: {0}\n'.format(obs['mode']))

    log_output(loggers, 1, logging.INFO, '   N PLANS: {0}\n'.format(plan.nplans))
    # TODO: Print plan details

    status = 0

    #-------------------------------------------------------------------
    # Declare the DRP fits file
    #-------------------------------------------------------------------
    drpf = DRPFits(obs['plate'], obs['ifudesign'], obs['mode'], redux_path=redux_path,
#    drpf = DRPFits(obs['plate'], obs['ifudesign'], 'RSS', redux_path=redux_path,
                   directory_path=directory_path)
    log_output(loggers, 1, logging.INFO, ' Read DRP file: {0}\n'.format(drpf.file_path()))

    if not os.path.isdir(analysis_path):
        os.makedirs(analysis_path)

    # Iterate over plans:
    for i in range(plan.nplans):

        #---------------------------------------------------------------
        # S/N Assessments
        #---------------------------------------------------------------
        rdxqa = None if plan['drpqa_key'][i] is None else \
                    ReductionAssessment(plan['drpqa_key'][i], drpf, pa=obs['pa'], ell=obs['ell'],
                                        dapsrc=dapsrc, analysis_path=analysis_path, verbose=verbose,
                                        clobber=plan['drpqa_clobber'][i])
 
        #---------------------------------------------------------------
        # Spatial Binning
        #---------------------------------------------------------------
        binned_spectra = None if plan['bin_key'][i] is None else \
                    SpatiallyBinnedSpectra(plan['bin_key'][i], drpf, rdxqa, reff=obs['reff'],
                                           dapsrc=dapsrc, analysis_path=analysis_path,
                                           verbose=verbose, clobber=plan['bin_clobber'][i])

        #---------------------------------------------------------------
        # Stellar Continuum Fit
        #---------------------------------------------------------------
        stellar_continuum = None if plan['continuum_key'][i] is None else \
                    StellarContinuumModel(plan['continuum_key'][i], drpf, binned_spectra,
                                          guess_vel=obs['vel'], guess_sig=obs['vdisp'],
                                          dapsrc=dapsrc, analysis_path=analysis_path,
                                          verbose=verbose, clobber=plan['continuum_clobber'][i])

        #---------------------------------------------------------------
        # Emission-line Moment measurements
        #---------------------------------------------------------------
        emission_line_moments = None if plan['elmom_key'][i] is None else \
                    EmissionLineMoments(plan['elmom_key'][i], binned_spectra,
                                        stellar_continuum=stellar_continuum, dapsrc=dapsrc,
                                        analysis_path=analysis_path, verbose=verbose,
                                        clobber=plan['elmom_clobber'][i])

        #---------------------------------------------------------------
        # Emission-line Fit
        #---------------------------------------------------------------
        emission_line_model = None if plan['elfit_key'][i] is None else \
                    EmissionLineModel(plan['elfit_key'][i], binned_spectra, guess_vel=obs['vel'],
                                      stellar_continuum=stellar_continuum, dapsrc=dapsrc,
                                      analysis_path=analysis_path, verbose=verbose,
                                      clobber=plan['elfit_clobber'][i])

        # Still need to add equivalent-width measurements

        #---------------------------------------------------------------
        # Spectral-Index Measurements
        #---------------------------------------------------------------
        spectral_indices = None if plan['spindex_key'][i] is None else \
                    SpectralIndices(plan['spindex_key'][i], binned_spectra,
                                    stellar_continuum=stellar_continuum,
                                    emission_line_model=emission_line_model, dapsrc=dapsrc,
                                    analysis_path=analysis_path, verbose=verbose,
                                    clobber=plan['spindex_clobber'][i])

        #-------------------------------------------------------------------
        # Construct the main output file(s)
        #-------------------------------------------------------------------
        construct_maps_file(drpf, rdxqa=rdxqa, binned_spectra=binned_spectra,
                            stellar_continuum=stellar_continuum,
                            emission_line_moments=emission_line_moments,
                            emission_line_model=emission_line_model,
                            spectral_indices=spectral_indices, dapsrc=dapsrc,
                            analysis_path=analysis_path, clobber=True)#clobber)

    #-------------------------------------------------------------------
    # End log
    #-------------------------------------------------------------------
    log_output(loggers, 1, logging.INFO, '    FINISH: {0}'.format(
                                    time.strftime("%a %d %b %Y %H:%M:%S",time.localtime())))
    if time.clock() - t < 60:
        tstr = '  DURATION: {0:.5f} sec'.format((time.clock() - t))
    elif time.clock() - t < 3600:
        tstr = '  DURATION: {0:.5f} min'.format((time.clock() - t)/60.)
    else:
        tstr = '  DURATION: {0:.5f} hr'.format((time.clock() - t)/3600.)
    log_output(loggers, 1, logging.INFO, tstr)


    mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss \
            + resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    if mem < 1024:
        mstr = 'MAX MEMORY: {0:.4f} bytes'.format(mem)
    elif mem < 1024*1024:
        mstr = 'MAX MEMORY: {0:.4f} Kbytes'.format(mem/1024)
    elif mem < 1024*1024*1024:
        mstr = 'MAX MEMORY: {0:.4f} Mbytes'.format(mem/1024/1024)
    else:
        mstr = 'MAX MEMORY: {0:.4f} Gbytes'.format(mem/1024/1024/1024)
    log_output(loggers, 1, logging.INFO, mstr)

    return status

