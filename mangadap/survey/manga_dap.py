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

from ..config import defaults
from ..util.log import init_DAP_logging, module_logging, log_output
from ..datacube import DataCube
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

def get_manga_dap_meta(cube):
    r"""
    Get the metadata required to run the DAP.

    The metadata is pulled from the provided DataCube and
    checked. The only required metadata keyword is ``z``, which
    sets the initial guess for the bulk redshift of the galaxy.
    If this key is not available or its value doesn't meet the
    criterion below, this function will raise an exception,
    meaning the DAP will fault before it starts processing the
    DataCube.

    The metadata provided by the DataCube must meet the following
    **critical criteria** or the method will fault:

        - Velocity (:math:`cz`) has to be greater than -500.

    For the remainder of the metadata, if the keyword does not
    exist, if the value is None, or if the value is outside the
    accepted range, a default is chosen. The metadata keywords,
    acceptable ranges, and defaults are provided below.

    +-----------+--------------------------------+---------+
    |   Keyword |                          Range | Default |
    +===========+================================+=========+
    | ``vdisp`` |             :math:`\sigma > 0` |     100 |
    +-----------+--------------------------------+---------+
    | ``ell``   | :math:`0 \leq \varepsilon < 1` |    None |
    +-----------+--------------------------------+---------+
    | ``pa``    |    :math:`0 \leq \phi_0 < 360` |    None |
    +-----------+--------------------------------+---------+
    | ``reff``  |        :math:`R_{\rm eff} > 0` |    None |
    +-----------+--------------------------------+---------+

    Returns:
        :obj:`dict`: Returns a dictionary with the following
        keywords:

            - ``z``: The bulk redshift of the galaxy, used to
              calculate :math:`cz`.
            - ``vel``: The initial guess velocity (:math:`cz`) in
              km/s.
            - ``vdisp``: The initial guess velocity dispersion in
              km/s.
            - ``ell``: The isophotal ellipticity (:math:`1-b/a`) to
              use when calculating semi-major axis coordinates.
            - ``pa``: The isophotal position angle in deg from N
              through E, used when calculating semi-major axis
              coordinates.
            - ``reff``: The effective radius in arcsec (DataCube WCS
              coordinates are expected to be in deg), used as a
              normalization of the semi-major axis radius in various
              output data.

    Raises:
        ValueError:
            See **critical criteria** above.
    """
    metadata = {}
    keys = cube.metakeys()
    if 'z' not in keys or cube.meta['z'] is None:
        raise ValueError('DataCube must provide bulk redshift metadata.')
    metadata['z'] = cube.meta['z']
    metadata['vel'] = astropy.constants.c.to('km/s').value * cube.meta['z']
    if metadata['vel'] < -500:
        raise ValueError('Velocity must be > -500 km/s!')

    if 'vdisp' not in keys or cube.meta['vdisp'] is None or not cube.meta['vdisp'] > 0:
        warnings.warn('Velocity dispersion not provided or not greater than 0 km/s.  '
                      'Adopting initial guess of 100 km/s.')
        metadata['vdisp'] = 100.0
    else:
        metadata['vdisp'] = cube.meta['vdisp']

    if 'ell' not in keys or cube.meta['ell'] is None or cube.meta['ell'] < 0 \
            or cube.meta['ell'] > 1:
        warnings.warn('Ellipticity not provided or not in the range 0 <= ell <= 1.  '
                      'Setting to None.')
        metadata['ell'] = None
    else:
        metadata['ell'] = cube.meta['ell']

    if 'pa' not in keys or cube.meta['pa'] is None:
        warnings.warn('Position angle not provided. Setting to None.')
        metadata['pa'] = None
    else:
        metadata['pa'] = cube.meta['pa']
    # Impose expected range
    if metadata['pa'] is not None and (metadata['pa'] < 0 or metadata['pa'] >= 360):
        warnings.warn('Imposing 0-360 range on position angle.')
        while metadata['pa'] < 0:
            metadata['pa'] += 360
        while metadata['pa'] >= 360:
            metadata['pa'] -= 360

    if 'reff' not in keys or cube.meta['reff'] is None or not cube.meta['reff'] > 0:
        warnings.warn('Effective radius not provided or not greater than 0. Setting to None.')
        metadata['reff'] = None
    else:
        metadata['reff'] = cube.meta['reff']

    return metadata


def manga_dap(cube, plan, dbg=False, log=None, verbose=0, drpver=None, redux_path=None,
              directory_path=None, dapver=None, analysis_path=None):
    r"""
    Main wrapper function for the MaNGA DAP.

    This function is designed to be called once per datacube. The
    :class:`mangadap.par.analysisplan.AnalysisPlanSet` instance sets
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
        plan (:class:`mangadap.par.analysisplan.AnalysisPlanSet`):
            Object with the analysis plans to implement.
        dbg (:obj:`bool`, optional):
            Flag to run the DAP in debug mode, see above; default is
            False. Limited use; still under development.
        log (:obj:`str`, optional):
            File name for log output, see above; no log file is
            produced by default.
        verbose (:obj:`int`, optional):
            Verbosity level, see above; default is 0.
        drpver (:obj:`str`, optional):
            DRP version. Default determined by
            :func:`mangadap.config.defaults.drp_version`.
        redux_path (:obj:`str`, optional):
            Top-level directory with the DRP products; default is
            defined by
            :func:`mangadap.config.defaults.drp_redux_path`.
        directory_path (:obj:`str`, optional):
            Direct path to directory containing the DRP output file;
            default is defined by
            :func:`mangadap.config.defaults.drp_directory_path`
        dapver (:obj:`str`, optional):
            DAP version. Default determined by
            :func:`mangadap.config.defaults.dap_version`.
        analysis_path (:obj:`str`, optional):
            Top-level directory for the DAP output data; default is
            defined by
            :func:`mangadap.config.defaults.dap_analysis_path`.

    Returns:
        :obj:`int`: Status flag (under development; currently always
        0).
    """
    # Check input
    if not isinstance(cube, DataCube):
        raise TypeError('Input objects must be a subclass of DataCube.')
    # Get the needed metadata
    metadata = get_manga_dap_meta(cube)
    if not isinstance(plan, AnalysisPlanSet):
        raise TypeError('plan must be of type AnalysisPlanSet')

    # Initialize the logging objects
    init_DAP_logging(log)#, simple_warnings=False)

    # Start log
    loggers = module_logging(__name__, verbose)

    _drpver = defaults.drp_version() if drpver is None else drpver
    _dapver = defaults.dap_version() if dapver is None else dapver
    # Set the the analysis path and make sure it exists
    _analysis_path = defaults.dap_analysis_path(drpver=_drpver, dapver=_dapver) \
                            if analysis_path is None else analysis_path
    if not os.path.isdir(_analysis_path):
        os.makedirs(_analysis_path)

    cube_root, cube_file = os.path.split(cube.file_path())

    log_output(loggers, 1, logging.INFO, '-'*50)
    log_output(loggers, 1, logging.INFO, '{0:^50}'.format('MaNGA Data Analysis Pipeline'))
    log_output(loggers, 1, logging.INFO, '-'*50)
    log_output(loggers, 1, logging.INFO, '     VERSION: {0}'.format(__version__))
    log_output(loggers, 1, logging.INFO, '       START: {0}'.format(
                                    time.strftime('%a %d %b %Y %H:%M:%S',time.localtime())))
    t = time.perf_counter()
    log_output(loggers, 1, logging.INFO, '  INPUT ROOT: {0}'.format(cube_root))
    log_output(loggers, 1, logging.INFO, '   CUBE FILE: {0}'.format(cube_file))
    log_output(loggers, 1, logging.INFO, '     N PLANS: {0}'.format(plan.nplans))
    log_output(loggers, 1, logging.INFO, ' OUTPUT ROOT: {0}'.format(_analysis_path))
    log_output(loggers, 1, logging.INFO, '-'*50)

    status = 0

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

        method = defaults.dap_method(bin_method['key'],
                                     sc_method['fitpar']['template_library_key'],
                                     'None' if el_method is None
                                        else el_method['continuum_tpl_key'])
        method_dir = defaults.dap_method_path(method, plate=cube.plate, ifudesign=cube.ifudesign,
                                              drpver=_drpver, dapver=_dapver,
                                              analysis_path=_analysis_path)
        method_ref_dir = defaults.dap_method_path(method, plate=cube.plate,
                                                  ifudesign=cube.ifudesign, ref=True,
                                                  drpver=_drpver, dapver=_dapver,
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
                    ReductionAssessment(plan['drpqa_key'][i], cube, pa=metadata['pa'],
                                        ell=metadata['ell'], analysis_path=_analysis_path,
                                        symlink_dir=method_ref_dir,
                                        clobber=plan['drpqa_clobber'][i], loggers=loggers)

        #---------------------------------------------------------------
        # Spatial Binning: placed in the common/ directory
        #---------------------------------------------------------------
        binned_spectra = None if plan['bin_key'][i] is None else \
                    SpatiallyBinnedSpectra(plan['bin_key'][i], cube, rdxqa, reff=metadata['reff'],
                                           analysis_path=_analysis_path,
                                           symlink_dir=method_ref_dir,
                                           clobber=plan['bin_clobber'][i], loggers=loggers)

        #---------------------------------------------------------------
        # Stellar Continuum Fit: placed in the common/ directory
        #---------------------------------------------------------------
        stellar_continuum = None if plan['continuum_key'][i] is None else \
                    StellarContinuumModel(plan['continuum_key'][i], binned_spectra,
                                          guess_vel=metadata['vel'], guess_sig=metadata['vdisp'],
                                          analysis_path=_analysis_path, symlink_dir=method_ref_dir,
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
                                        stellar_continuum=stellar_continuum,
                                        redshift=metadata['z'], analysis_path=_analysis_path,
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
                                     analysis_path=_analysis_path, directory_path=method_ref_dir,
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
                                          emission_line_model=emission_line_model,
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
                    SpectralIndices(plan['spindex_key'][i], binned_spectra,
                                    redshift=metadata['z'], stellar_continuum=stellar_continuum,
                                    emission_line_model=emission_line_model,
                                    analysis_path=_analysis_path, directory_path=method_ref_dir,
                                    tpl_symlink_dir=method_ref_dir,
                                    clobber=plan['spindex_clobber'][i], loggers=loggers)

        #-------------------------------------------------------------------
        # Construct the main output file(s)
        #-------------------------------------------------------------------
        construct_maps_file(cube, metadata, rdxqa=rdxqa, binned_spectra=binned_spectra,
                            stellar_continuum=stellar_continuum,
                            emission_line_moments=emission_line_moments,
                            emission_line_model=emission_line_model,
                            spectral_indices=spectral_indices, redshift=metadata['z'],
                            analysis_path=_analysis_path, clobber=True, loggers=loggers,
                            single_precision=True)

        construct_cube_file(cube, metadata, binned_spectra=binned_spectra,
                            stellar_continuum=stellar_continuum,
                            emission_line_model=emission_line_model,
                            analysis_path=_analysis_path, clobber=True, loggers=loggers,
                            single_precision=True)

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


