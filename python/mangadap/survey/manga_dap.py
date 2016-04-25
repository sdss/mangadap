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
from ..proc.spectralindices import SpectralIndices
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

    init_DAP_logging(log)#, simple_warnings=False)

    # Start log
    loggers = module_logging(__name__, verbose)

    log_output(loggers, 1, logging.INFO, '   DAPVERS: {0}'.format(__version__))
    log_output(loggers, 1, logging.INFO, '     START: {0}'.format(
                                    time.strftime("%a %d %b %Y %H:%M:%S",time.localtime())))
    t = time.clock()
    log_output(loggers, 1, logging.INFO, '     PLATE: {0}'.format(obs['plate']))
    log_output(loggers, 1, logging.INFO, ' IFUDESIGN: {0}'.format(obs['ifudesign']))
    log_output(loggers, 1, logging.INFO, '   3D MODE: {0}\n'.format(obs['mode']))

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

#    n_plans = plan.complete_plans
    n_plans = 1

    # Iterate over plans:
    for i in range(n_plans):

        #-------------------------------------------------------------------
        # S/N Assessments
        #-------------------------------------------------------------------
        snr = ReductionAssessment('RFWHMC', drpf, pa=obs['pa'], ell=obs['ell'], dapsrc=dapsrc,
                                  analysis_path=analysis_path, verbose=verbose, clobber=False)
 
        #-------------------------------------------------------------------
        # Spatial Binning
        #-------------------------------------------------------------------
        binned_spectra = SpatiallyBinnedSpectra('LOGR10C', drpf, snr, reff=obs['reff'],
                                                dapsrc=dapsrc, analysis_path=analysis_path,
                                                verbose=verbose, clobber=False)

#        print(binned_spectra['MASK'].data)
#
#        print(binned_spectra.method)
#        flux = binned_spectra.copy_to_array()
#        #flux = binned_spectra.copy_to_masked_array()
#        print(numpy.sum(numpy.ma.getmaskarray(flux)))
#        print(flux.shape)
#        for i in range(binned_spectra.nbins):
#            pyplot.step(binned_spectra['WAVE'].data, flux[i,:], where='mid', linestyle='-')
#        pyplot.show()


        #-------------------------------------------------------------------
        # Stellar Continuum Fit
        #-------------------------------------------------------------------
        stellar_continuum = StellarContinuumModel('GAU01-MILESTHIN', drpf, binned_spectra,
                                                  guess_vel=obs['vel'], guess_sig=obs['vdisp'],
                                                  dapsrc=dapsrc, analysis_path=analysis_path,
                                                  verbose=verbose, clobber=False)

#        flux = binned_spectra.copy_to_masked_array()
#        model = stellar_continuum.copy_to_masked_array()
#        wave = binned_spectra['WAVE'].data
#        flux = binned_spectra.copy_to_masked_array(flag=binned_spectra.do_not_fit_flags())
#        model = stellar_continuum.copy_to_masked_array(
#                                                flag=stellar_continuum.all_except_emission_flags())
#        disp_free = numpy.ma.MaskedArray(stellar_continuum.construct_models(redshift_only=True),
#                                         mask=model.mask.copy())
#        for i in range(binned_spectra.nbins):
#            pyplot.step(wave, flux[i,:], where='mid', linestyle='-', color='k', lw=0.5, zorder=3)
#            pyplot.plot(wave, model[i,:], linestyle='-', color='b', lw=3, zorder=1, alpha=0.5)
#            pyplot.plot(wave, disp_free[i,:], linestyle='-', color='r', lw=1, zorder=1)
#            pyplot.show()

#
#       #-------------------------------------------------------------------
#       # Emission-line Moment measurements
#       #-------------------------------------------------------------------
#       emission_line_moments = EmissionLineMoments(drpf, plan.emission, snr=snr,
#                                                   binned_spectra=binned_spectra,
#                                                   stellar_continuum=stellar_continuum)
#
#       #-------------------------------------------------------------------
#       # Emission-line Fit
#       #-------------------------------------------------------------------
#       emission_lines_fits = EmissionLineFits(drpf, plan.emission, snr=snr,
#                                              binned_spectra=binned_spectra,
#                                              stellar_continuum=stellar_continuum,
#                                              emission_line_moments=emission_line_moments)
#
        #-------------------------------------------------------------------
        # Spectral-Index Measurements
        #-------------------------------------------------------------------
        spectral_indices = SpectralIndices('STANDARD', binned_spectra,
                                           stellar_continuum=stellar_continuum, dapsrc=dapsrc,
                                           analysis_path=analysis_path, verbose=verbose,
                                           clobber=False)

    #-------------------------------------------------------------------
    # Construct the main output file
    #-------------------------------------------------------------------


    # End log
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

