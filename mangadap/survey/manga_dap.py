# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Provides the main wrapper function for the MaNGA DAP.

Revision history
----------------

    | **21 Mar 2016**: Original implementation by K. Westfall (KBW): started.
    | **19 Jul 2016**: (KBW) Always provide the NSA redshift to the modules!
    | **05 May 2017**: (KBW) Clean up documentation; use import to get
        __version__

----

.. include license and copyright
.. include:: ../copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../links.rst
"""
import logging
import resource
import time
import os
import warnings
import numpy

import astropy.constants

from mangadap import __version__

from ..config import defaults
from ..util.log import init_DAP_logging, module_logging, log_output
from ..drpfits import DRPFits
from ..par.obsinput import ObsInputPar
from ..par.analysisplan import AnalysisPlanSet
from ..proc.reductionassessments import ReductionAssessment
from ..proc.spatiallybinnedspectra import SpatiallyBinnedSpectra
from ..proc.stellarcontinuummodel import StellarContinuumModel
from ..proc.emissionlinemoments import EmissionLineMoments
from ..proc.emissionlinemodel import EmissionLineModel
from ..proc.spectralindices import SpectralIndices
from ..dapfits import construct_maps_file, construct_cube_file

# For testing/debugging
from ..util.fitsutil import DAPFitsUtil
from matplotlib import pyplot

def manga_dap(obs, plan, dbg=False, log=None, verbose=0, drpver=None, redux_path=None,
              directory_path=None, dapver=None, dapsrc=None, analysis_path=None):
    r"""
    Main wrapper function for the MaNGA DAP.

    This function is designed to be called once per PLATE-IFUDESIGN as
    set by the provided :class:`mangadap.par.obsinput.ObsInputPar`
    instance.  The :class:`mangadap.par.analysisplan.AnalysisPlanSet`
    instance sets the types of analyses to perform on this observation.
    Each analysis plan results in a MAPS and model LOGCUBE file.

    The procedure is as follows:
        - Read the DRP fits file specified by `obs`
        - For each plan in `plan`:
            - Determine basic assessments of the data, including the S/N
              and spaxel coordinates.  See
              :class:`mangadap.proc.reductionassessments.ReductionAssessment`.
            - Bin the spectra.  See
              :class:`mangadap.proc.spatiallybinnedspectra.SpatiallyBinnedSpectra`.
            - Fit the stellar continuum for stellar kinematics.  See
              :class:`mangadap.proc.stellarcontinuummodel.StellarContinuumModel`.
            - Measure the emission-line moments.  See
              :class:`mangadap.proc.emissionlinemoments.EmissionLineMoments`.
            - Fit parameterized line profiles to the emission lines.
              See
              :class:`mangadap.proc.emissionlinemodel.EmissionLineModel`.
            - Subtract the fitted emission-line models and measure the
              spectral indices.  See
              :class:`mangadap.proc.spectralindices.SpectralIndices`.
            - Construct the primary output files based on the plan
              results.  See :func:`mangadap.dapfits.construct_maps_file`
              and :func:`mangadap.dapfits.construct_cube_file`.

    Verbose levels (still under development):
        0. Nothing but errors.
        1. Basic status updates. E.g., start and end of each block,
           minor progress within blocks.
        2. Block-level and sub-function updates.
        3. Same as above, with figures.

    Args:
        obs (dict, :class:`mangadap.par.obsinput.ObsInputPar`): Object
            with the input parameters.
        plan (:class:`mangadap.par.analysisplan.AnalysisPlanSet`):
            Object with the analysis plans to implement.
        dbg (bool) : (**Optional**) Flag to run the DAP in debug mode,
            see above; default is False.  Limited use; still under
            development.
        log (str) : (**Optional**) File name for log output, see above;
            no log file is produced by default.
        verbose (int) : (**Optional**) Verbosity level, see above;
            default is 0.
        drpver (str): (**Optional**) DRP version.  Default determined by
            :func:`mangadap.config.defaults.default_drp_version`.
        redux_path (str) : (**Optional**) Top-level directory with the
            DRP products; default is defined by
            :func:`mangadap.config.defaults.default_redux_path`.
        directory_path (str) : (**Optional**) Direct path to directory
            containing the DRP output file; default is defined by
            :func:`mangadap.config.defaults.default_drp_directory_path`
        dapver (str): (**Optional**) DAP version.  Default determined by
            :func:`mangadap.config.defaults.default_dap_version`.
        dapsrc (str): (**Optional**) Source directory of the DAP.
            Default determined by
            :func:`mangadap.config.defaults.dap_source_dir`.
        analysis_path (str) : (**Optional**) Top-level directory for the DAP
            output data; default is defined by
            :func:`mangadap.config.defaults.default_analysis_path`.

    Returns:
        int: Status flag (under development; currently always 0)

    """
    # Check input
    if not isinstance(obs, (ObsInputPar, dict)):
        raise TypeError('obs must be of type dict or ObsInputPar')
    if not isinstance(plan, AnalysisPlanSet):
        raise TypeError('plan must be of type AnalysisPlanSet')

    # Initialize the logging objects
    init_DAP_logging(log)#, simple_warnings=False)

    # Start log
    loggers = module_logging(__name__, verbose)

    log_output(loggers, 1, logging.INFO, '-'*50)
    log_output(loggers, 1, logging.INFO, '{0:^50}'.format('MaNGA Data Analysis Pipeline'))
    log_output(loggers, 1, logging.INFO, '-'*50)
    log_output(loggers, 1, logging.INFO, '   VERSION: {0}'.format(__version__))
    log_output(loggers, 1, logging.INFO, '     START: {0}'.format(
                                    time.strftime("%a %d %b %Y %H:%M:%S",time.localtime())))
    t = time.perf_counter()
    log_output(loggers, 1, logging.INFO, '     PLATE: {0}'.format(obs['plate']))
    log_output(loggers, 1, logging.INFO, ' IFUDESIGN: {0}'.format(obs['ifudesign']))
    log_output(loggers, 1, logging.INFO, '   3D MODE: {0}'.format(obs['mode']))
    log_output(loggers, 1, logging.INFO, '   N PLANS: {0}'.format(plan.nplans))
    log_output(loggers, 1, logging.INFO, '-'*50)

    status = 0
    if obs['mode'] != 'CUBE':
        status = 1
        raise NotImplementedError('DAP can currently only analyze CUBE files.')

    #-------------------------------------------------------------------
    # Declare the DRP fits file
    #-------------------------------------------------------------------
    _drpver = defaults.default_drp_version() if drpver is None else drpver
    _dapver = defaults.default_dap_version() if dapver is None else dapver
    drpf = DRPFits(obs['plate'], obs['ifudesign'], obs['mode'], drpver=_drpver,
                   redux_path=redux_path, directory_path=directory_path, read=True)
    log_output(loggers, 1, logging.INFO, ' Opened DRP file: {0}\n'.format(drpf.file_path()))

    # Test if the RSS file exists
    drpf_rss = DRPFits(obs['plate'], obs['ifudesign'], 'RSS', drpver=_drpver,
                       redux_path=redux_path, directory_path=directory_path, read=False)
    if not os.path.isfile(drpf_rss.file_path()):
        warnings.warn('RSS counterpart not available.  Some functionality may be limited!')
    del drpf_rss

    # Set the the analysis path and make sure it exists
    _analysis_path = defaults.default_analysis_path(drpver=drpver, dapver=dapver) \
                            if analysis_path is None else analysis_path
    if not os.path.isdir(_analysis_path):
        os.makedirs(_analysis_path)

    # Set the nsa redshift
    nsa_redshift = obs['vel']/astropy.constants.c.to('km/s').value

    # Iterate over plans:
    for i in range(plan.nplans):

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
        bin_method = SpatiallyBinnedSpectra.define_method(plan['bin_key'][i])
        sc_method = StellarContinuumModel.define_method(plan['continuum_key'][i])
        el_method = None if plan['elfit_key'][i] is None \
                            else EmissionLineModel.define_method(plan['elfit_key'][i])

        method = defaults.default_dap_method(bin_method['key'],
                                             sc_method['fitpar']['template_library_key'],
                                             'None' if el_method is None
                                                else el_method['continuum_tpl_key'])
        method_dir = defaults.default_dap_method_path(method, plate=obs['plate'],
                                                      ifudesign=obs['ifudesign'], drpver=drpver,
                                                      dapver=dapver, analysis_path=_analysis_path)
        method_ref_dir = defaults.default_dap_method_path(method, plate=obs['plate'],
                                                          ifudesign=obs['ifudesign'], ref=True,
                                                          drpver=drpver, dapver=dapver,
                                                          analysis_path=_analysis_path)

        log_output(loggers, 1, logging.INFO, '-'*50)
        log_output(loggers, 1, logging.INFO, '{0:^50}'.format('Plan {0}'.format(i+1)))
        log_output(loggers, 1, logging.INFO, '-'*50)
        log_output(loggers, 1, logging.INFO, '    METHOD: {0}'.format(method))
        log_output(loggers, 1, logging.INFO, '    OUTPUT: {0}'.format(method_dir))
        log_output(loggers, 1, logging.INFO, 'REF OUTPUT: {0}'.format(method_ref_dir))

        #---------------------------------------------------------------
        # S/N Assessments: placed in the common/ directory
        #---------------------------------------------------------------
        rdxqa = None if plan['drpqa_key'][i] is None else \
                    ReductionAssessment(plan['drpqa_key'][i], drpf, pa=obs['pa'], ell=obs['ell'],
                                        dapsrc=dapsrc, analysis_path=_analysis_path,
                                        symlink_dir=method_ref_dir,
                                        clobber=plan['drpqa_clobber'][i], loggers=loggers)

        #---------------------------------------------------------------
        # Spatial Binning: placed in the common/ directory
        #---------------------------------------------------------------
        binned_spectra = None if plan['bin_key'][i] is None else \
                    SpatiallyBinnedSpectra(plan['bin_key'][i], drpf, rdxqa, reff=obs['reff'],
                                           dapsrc=dapsrc, analysis_path=_analysis_path,
                                           symlink_dir=method_ref_dir,
                                           clobber=plan['bin_clobber'][i], loggers=loggers)

        #---------------------------------------------------------------
        # Stellar Continuum Fit: placed in the common/ directory
        #---------------------------------------------------------------
        stellar_continuum = None if plan['continuum_key'][i] is None else \
                    StellarContinuumModel(plan['continuum_key'][i], binned_spectra,
                                          guess_vel=obs['vel'], guess_sig=obs['vdisp'],
                                          dapsrc=dapsrc, analysis_path=_analysis_path,
                                          symlink_dir=method_ref_dir,
                                          tpl_symlink_dir=method_ref_dir,
                                          clobber=plan['continuum_clobber'][i], loggers=loggers)

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
        emission_line_moments = None if plan['elmom_key'][i] is None else \
                    EmissionLineMoments(plan['elmom_key'][i], binned_spectra,
                                        stellar_continuum=stellar_continuum, redshift=nsa_redshift,
                                        dapsrc=dapsrc, analysis_path=_analysis_path,
                                        directory_path=method_ref_dir,
                                        clobber=plan['elmom_clobber'][i], loggers=loggers)

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
        emission_line_model = None if plan['elfit_key'][i] is None else \
                   EmissionLineModel(plan['elfit_key'][i], binned_spectra,
                                     stellar_continuum=stellar_continuum,
                                     emission_line_moments=emission_line_moments, dispersion=100.0,
                                     minimum_error=numpy.finfo(numpy.float32).eps,
                                     dapsrc=dapsrc, analysis_path=_analysis_path,
                                     directory_path=method_ref_dir,
                                     clobber=plan['elfit_clobber'][i], loggers=loggers)

#        model_flux = emission_line_model['FLUX'].data
#        model_base = emission_line_model['BASE'].data
#        model_mask = emission_line_model['MASK'].data
#        model_binid = emission_line_model['BINID'].data
#        eml_par = emission_line_model['EMLDATA'].data
#
#        test = emission_line_model._get_line_fit_metrics(model_flux, model_base, model_mask, eml_par, model_binid)
#

#        cen = drpf.spatial_shape[0]//2
#        indx = emission_line_model['BINID'].data[cen,cen]
#        indx = numpy.arange(len(emission_line_model['EMLDATA'].data['BINID']))[emission_line_model['EMLDATA'].data['BINID'] == indx][0]
#        print(indx)
#        pyplot.plot(emission_line_model['WAVE'].data, emission_line_model['FLUX'].data[indx,:])
#        pyplot.show()
#
#        el_continuum = emission_line_model.fill_continuum_to_match(
#                                            emission_line_model['BINID'].data,
#                                            missing=emission_line_model.missing_models)
#        el_continuum_dcnvlv = emission_line_model.fill_continuum_to_match(
#                                            emission_line_model['BINID'].data,
#                                            missing=emission_line_model.missing_models,
#                                            redshift_only=True)
#        el_continuum_dcnvlv_rest = emission_line_model.fill_continuum_to_match(
#                                            emission_line_model['BINID'].data,
#                                            missing=emission_line_model.missing_models,
#                                            redshift_only=True, deredshift=True)
#
#        pyplot.imshow(el_continuum, origin='lower', interpolation='nearest', aspect='auto')
#        pyplot.colorbar()
#        pyplot.show()
#
#        pyplot.imshow(el_continuum_dcnvlv, origin='lower', interpolation='nearest', aspect='auto')
#        pyplot.colorbar()
#        pyplot.show()
#
#        pyplot.imshow(el_continuum_dcnvlv_rest, origin='lower', interpolation='nearest', aspect='auto')
#        pyplot.colorbar()
#        pyplot.show()
#
#        pyplot.plot(binned_spectra['WAVE'].data, el_continuum[indx,:])
#        pyplot.plot(binned_spectra['WAVE'].data, el_continuum_dcnvlv[indx,:])
#        pyplot.plot(binned_spectra['WAVE'].data, el_continuum_dcnvlv_rest[indx,:])
#        pyplot.show()

#        model_flux = emission_line_model['FLUX'].data
#        model_base = emission_line_model['BASE'].data
#        model_mask = emission_line_model['MASK'].data
#        model_eml_par = emission_line_model['EMLDATA'].data
#        model_binid = emission_line_model['BINID'].data \
#                        if emission_line_model.method['deconstruct_bins'] else None
#        test_eml_par = emission_line_model._get_line_fit_metrics(model_flux, model_base,
#                                                                 model_mask, model_eml_par,
#                                                                 model_binid)

#        el_continuum = emission_line_model.fill_continuum_to_match(binned_spectra['BINID'].data,
#                                                        missing=binned_spectra.missing_bins)
#        sc_continuum = stellar_continuum.fill_to_match(binned_spectra['BINID'].data,
#                                               missing=binned_spectra.missing_bins)
#        pyplot.plot(binned_spectra['WAVE'].data, binned_spectra['FLUX'].data[0,:])
#        pyplot.plot(binned_spectra['WAVE'].data, el_continuum[0,:])
#        pyplot.plot(binned_spectra['WAVE'].data, sc_continuum[0,:])
#        pyplot.show()

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
                                          emission_line_model=emission_line_model, dapsrc=dapsrc,
                                          analysis_path=_analysis_path,
                                          directory_path=method_ref_dir,
                                          clobber=plan['elmom_clobber'][i], loggers=loggers)

#        print(emission_line_moments['ELMMNTS'].data['FLUX'].shape)
#
#        pyplot.scatter(ha_flux,
#                       numpy.ma.divide(emission_line_moments['ELMMNTS'].data['FLUX'][:,18],ha_flux),
#                       marker='.', s=20, lw=0)
#        pyplot.show()
#
#        pyplot.scatter(emission_line_model['EMLDATA'].data['FLUX'][:,18],
#                       numpy.ma.divide(ha_flux,emission_line_model['EMLDATA'].data['FLUX'][:,18]),
#                       marker='.', s=20, lw=0)
#        pyplot.show()
#
#        pyplot.scatter(emission_line_model['EMLDATA'].data['FLUX'][:,18],
#                       numpy.ma.divide(emission_line_moments['ELMMNTS'].data['FLUX'][:,18],
#                                       emission_line_model['EMLDATA'].data['FLUX'][:,18]),
#                       marker='.', s=20, lw=0)
#        pyplot.show()

        #---------------------------------------------------------------
        # Spectral-Index Measurements: placed in the DAPTYPE/ref/
        # directory
        #---------------------------------------------------------------
        spectral_indices = None if plan['spindex_key'][i] is None else \
                    SpectralIndices(plan['spindex_key'][i], binned_spectra, redshift=nsa_redshift,
                                    stellar_continuum=stellar_continuum,
                                    emission_line_model=emission_line_model, dapsrc=dapsrc,
                                    analysis_path=_analysis_path, directory_path=method_ref_dir,
                                    tpl_symlink_dir=method_ref_dir,
                                    clobber=plan['spindex_clobber'][i], loggers=loggers)

        #-------------------------------------------------------------------
        # Construct the main output file(s)
        #-------------------------------------------------------------------
        construct_maps_file(drpf, obs=obs, rdxqa=rdxqa, binned_spectra=binned_spectra,
                            stellar_continuum=stellar_continuum,
                            emission_line_moments=emission_line_moments,
                            emission_line_model=emission_line_model,
                            spectral_indices=spectral_indices, nsa_redshift=nsa_redshift,
                            dapsrc=dapsrc, analysis_path=_analysis_path, clobber=True,
                            loggers=loggers, single_precision=True)

        construct_cube_file(drpf, obs=obs, binned_spectra=binned_spectra,
                            stellar_continuum=stellar_continuum,
                            emission_line_model=emission_line_model,
                            dapsrc=dapsrc, analysis_path=_analysis_path, clobber=True,
                            loggers=loggers, single_precision=True)

        # Force memory to be freed
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
    log_output(loggers, 1, logging.INFO, '    FINISH: {0}'.format(
                                    time.strftime("%a %d %b %Y %H:%M:%S",time.localtime())))
    if time.perf_counter() - t < 60:
        tstr = '  DURATION: {0:.5f} sec'.format((time.perf_counter() - t))
    elif time.perf_counter() - t < 3600:
        tstr = '  DURATION: {0:.5f} min'.format((time.perf_counter() - t)/60.)
    else:
        tstr = '  DURATION: {0:.5f} hr'.format((time.perf_counter() - t)/3600.)
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
    log_output(loggers, 1, logging.INFO, '-'*50)
    log_output(loggers, 1, logging.INFO, '-'*50)

    return status


