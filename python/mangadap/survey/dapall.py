# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
r"""
Defines the class that constructs the DAPall summary table.

This is a post processing script that **must** be run after the DRPall
file is created.

*License*:
    Copyright (c) 2016, SDSS-IV/MaNGA Pipeline Group
        Licensed under BSD 3-clause license - see LICENSE.rst

*Source location*:
    $MANGADAP_DIR/python/mangadap/dapall.py

*Imports and python version compliance*:
    ::

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

        from scipy import interpolate

        from astropy.io import fits
        import astropy.constants
        from astropy.cosmology import FlatLambdaCDM
        import astropy.units
        from astropy.stats import sigma_clip

        from ..drpfits import DRPFits, DRPQuality3DBitMask
        from ..dapfits import DAPMapsBitMask, DAPQualityBitMask
        from ..par.analysisplan import AnalysisPlanSet
        from .drpcomplete import DRPComplete
        from ..config import defaults
        from ..util.fileio import init_record_array, rec_to_fits_type, channel_dictionary
        from ..proc.util import select_proc_method
        from ..par.absorptionindexdb import AbsorptionIndexDB
        from ..par.bandheadindexdb import BandheadIndexDB
        from ..par.emissionlinedb import EmissionLineDB
        from ..proc.spectralindices import SpectralIndicesDef, available_spectral_index_databases
        from ..proc.emissionlinemodel import EmissionLineModelDef
        from ..proc.emissionlinemodel import available_emission_line_modeling_methods

*Class usage examples*:
    Add some usage comments here!

*Revision history*:
    | **19 Aug 2016**: Original Implementation by K. Westfall (KBW)

.. _astropy.io.fits.open: http://docs.astropy.org/en/stable/io/fits/api/files.html#astropy.io.fits.open
.. _astropy.io.fits.hdu.hdulist.HDUList: http://docs.astropy.org/en/v1.0.2/io/fits/api/hdulists.html
.. _astropy.cosmology.FlatLambdaCDM: http://docs.astropy.org/en/stable/api/astropy.cosmology.FlatLambdaCDM.html#flatlambdacdm
.. _logging.Logger: https://docs.python.org/3/library/logging.html

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

from scipy import interpolate

from astropy.io import fits
import astropy.constants
from astropy.cosmology import FlatLambdaCDM
import astropy.units
from astropy.stats import sigma_clip

from .drpcomplete import DRPComplete
from ..drpfits import DRPFits, DRPQuality3DBitMask
from ..dapfits import DAPMapsBitMask, DAPQualityBitMask
from ..config import defaults
from ..util.fileio import init_record_array, rec_to_fits_type, channel_dictionary, channel_units
from ..util.log import init_DAP_logging, module_logging, log_output
from ..par.analysisplan import AnalysisPlanSet
from ..par.absorptionindexdb import AbsorptionIndexDB
from ..par.bandheadindexdb import BandheadIndexDB
from ..par.emissionlinedb import EmissionLineDB
from ..proc.util import select_proc_method, sample_growth
from ..proc.spectralindices import SpectralIndicesDef, available_spectral_index_databases
from ..proc.emissionlinemodel import EmissionLineModelDef
from ..proc.emissionlinemodel import available_emission_line_modeling_methods

from matplotlib import pyplot

#def growth(a, growth_fracs, default=-9999.):
#    _a = a.compressed() if isinstance(a, numpy.ma.MaskedArray) else numpy.atleast_1d(a).ravel()
##    print(len(_a))
#    if len(_a) < 2:
#        return tuple([default]*len(growth_fracs))
#    srt = numpy.argsort(_a)
#    grw = (numpy.arange(_a.size,dtype=float)+1)/_a.size
##    print(grw[0:2], '...', grw[-2:])
#    interpolator = interpolate.interp1d(grw, _a[srt], fill_value='extrapolate')
#    return tuple(interpolator(growth_fracs))

class DAPall:
    """
    Construct the summary table information for the DAP.

    Any observation in the DRPComplete file with:

        - MaNGAID != 'NULL'
        - MANGA_TARGET1 > 0 or MANGA_TARGET3 > 0
        - VEL > 0

    are assumed to have been analyzed by the DAP.  The success flag only
    checks that the appropriate maps file exists. 

    Args:
        plan (:class:`mangadap.par.analysisplan.AnalysisPlanSet`): The
            plan object used by the DAP.
        methods (str,list): (**Optional**) Specify a set of methods in
            the DAP plan file to include in the DAPall file.  If not
            provided, all methods in the plan file are included.
        drpver (str): (**Optional**) DRP version.  Default determined by
            :func:`mangadap.config.defaults.default_drp_version`.
        redux_path (str) : (**Optional**) Top-level directory with the
            DRP products; default is defined by
            :func:`mangadap.config.defaults.default_redux_path`.
        dapver (str): (**Optional**) DAP version.  Default determined by
            :func:`mangadap.config.defaults.default_dap_version`.
        dapsrc (str): (**Optional**) Source directory of the DAP.
            Default determined by
            :func:`mangadap.config.defaults.default_dap_source`.
        analysis_path (str) : (**Optional**) Top-level directory for the DAP
            output data; default is defined by
            :func:`mangadap.config.defaults.default_analysis_path`.
        dap_common_path (str): (**Optional**) Path to the DAP `common`
            directory.  Default is set by
            :func:`mangadap.config.defaults.default_dap_common_path`.
        readonly(bool): (**Optional**) If it exists, open any existing
            file and disallow any modifications of the database.
            Default is to automatically check for any need to update the
            file.
        loggers (list): List of `logging.Logger`_ objects used to log
            progress.
        quiet (bool): Suppress all terminal and logging output.

    Raises:
        TypeError: Raised if the input plan is not a
            :class:`mangadap.par.analysisplan.AnalysisPlanSet` object.
        FileNotFoundError: Raise if the DRPall or the DRPComplete file
            are not available.
    
    Attributes:
        drpver (str): DRP version
        redux_path (str): Path to top-level redux directory
        dapsrc (str): Path to DAP source distribution
        dapver (str): DAP version
        analysis_path (str): Path to top-level analysis directory
        dap_common_path (str): Path to DAP common directory
        drpall_file (str): Path to the DRPall file
        plan (:class:`mangadap.par.analysisplan.AnalysisPlanSet`): The
            plan object used by the DAP.
        plan_methods (numpy.ndarray): Array of all methods read from the
            provided :class:`AnalysisPlanSet`
        methods (numpy.ndarray): Array of the methods included in the
            DAPall file.
        h (float): The Hubble parameter (currently always set to unity)
        H0 (float): Hubbles constant (km/s/Mpc)
        cosmo (`astropy.cosmology.FlatLambdaCDM`_): Object used for
            distance calculations
        maps_bm (:class:`mangadap.dapfits.DAPMapsBitMask`): Bitmask used
            to mask map pixels.
        neml (int): Number of emission lines included in each method.
        nindx (int): Number of spectral indices included in each method.
        str_len (dict): Dictionary with the number of characters used
            for each string entry in the database.
        hdu (`astropy.io.fits.hdu.hdulist.HDUList`_): Object with the
            table data.
        ndap (int): Number of rows in the data table; equivalent to the
            number of processed MAPS files.
        plate (numpy.ndarray): Array of the plate numbers
        ifudesign (numpy.ndarray): Array of the ifudesigns
        readonly (bool): Object is prohibited from updating the
            database.
        loggers (list): List of `logging.Logger`_ objects
            to log progress; ignored if quiet=True.  Logging is done
            using :func:`mangadap.util.log.log_output`.  Default is no
            logging.
        quiet (bool): (**Optional**) Suppress all terminal and logging
            output.  Default is False.

    """
    def __init__(self, plan, methods=None, drpver=None, redux_path=None, dapver=None, dapsrc=None,
                 analysis_path=None, dap_common_path=None, readonly=False, loggers=None,
                 quiet=False):

        # Initialize the reporting
        if loggers is not None:
            self.loggers = loggers
        self.quiet = quiet

        # Check the input type
        # TODO: Plan isn't really needed.  Could just troll the
        # analysis_path for the directory names...
        if not isinstance(plan, AnalysisPlanSet):
            raise TypeError('Input plan must have type AnalysisPlanSet.')
        self.plan = plan

        # Path definitions
        self.drpver = defaults.default_drp_version() if drpver is None else str(drpver)
        self.redux_path = defaults.default_redux_path(self.drpver) if redux_path is None \
                                                                  else str(redux_path)

        self.dapsrc = defaults.default_dap_source() if dapsrc is None else str(dapsrc)

        self.dapver = defaults.default_dap_version() if dapver is None else str(dapver)
        self.analysis_path = defaults.default_analysis_path(self.drpver, self.dapver) \
                             if analysis_path is None else str(analysis_path)
        self.dap_common_path = defaults.default_dap_common_path(drpver=self.drpver,
                                            dapver=self.dapver, analysis_path=self.analysis_path) \
                             if dap_common_path is None else str(dap_common_path)

        # Set the name of the DRPall file
        self.drpall_file = os.path.join(self.redux_path, 'drpall-{0}.fits'.format(self.drpver))
        if not os.path.isfile(self.drpall_file):
            raise FileNotFoundError('{0} does not exist!'.format(self.drpall_file))

        # Get the full list of available plan methods
        self.plan_methods = numpy.array([ defaults.default_dap_method(plan=p) for p in self.plan ])

        # Check that the plan methods are unique.  This should never be
        # raised, unless it also causes problems in the DAP itself.
        unique_methods = numpy.unique(self.plan_methods)
        if len(unique_methods) != len(self.plan_methods):
            raise ValueError('All of the plan methods must be unique!')

        # Initialize the cosmology object
        self.h = 1.0
        self.H0 = 100 * self.h * astropy.units.km / astropy.units.s / astropy.units.Mpc
        self.cosmo = FlatLambdaCDM(H0=self.H0, Om0=0.3)

        # Initialize the bitmask
        self.maps_bm = DAPMapsBitMask()

        # Set the length for the string elements of the table
        self.str_len = self._init_string_lengths()

        # Declare the other attributes for use later on
        self.neml = None
        self.nindx = None
        self.hdu = None
        self.ndap = None
        self.methods = None
        self.plate = None
        self.ifudesign = None

        # Give an initial report
        self.readonly = readonly
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, '-'*50)
            log_output(self.loggers, 1, logging.INFO, '{0:^50}'.format('DAPALL DATABASE'))
            log_output(self.loggers, 1, logging.INFO, '-'*50)
            log_output(self.loggers, 1, logging.INFO, 'DRPall file: {0}'.format(self.drpall_file))
            log_output(self.loggers, 1, logging.INFO, 'Plan methods: {0}'.format(self.plan_methods))
            log_output(self.loggers, 1, logging.INFO, 'Output path: {0}'.format(self.analysis_path))
            log_output(self.loggers, 1, logging.INFO, 'Output file: {0}'.format(self.file_name()))
            if self.readonly:
                log_output(self.loggers, 1, logging.INFO, 'READ-ONLY MODE')

        # If the output file already exists, read it
        if os.path.isfile(self.file_path()):
            self._read()

        if self.readonly and self.hdu is None:
            raise ValueError('Selected read-only, but no database can be read.')

        # Automatically attempt to update the database if no database
        # exists or if the database is not meant to be read-only
        if self.hdu is None or not self.readonly:
            self.update(methods=methods)


    def __getitem__(self, key):
        return self.hdu['DAPALL'].data[key]


    @staticmethod
    def _number_of_spectral_indices(key, dapsrc=None):
        if key == 'None':
            return 0
        db = select_proc_method(key, SpectralIndicesDef,
                                available_func=available_spectral_index_databases, dapsrc=dapsrc)
        nabs = 0 if db['absindex'] is None \
                    else AbsorptionIndexDB(db['absindex'], dapsrc=dapsrc).nsets
        nbhd = 0 if db['bandhead'] is None \
                    else BandheadIndexDB(db['bandhead'], dapsrc=dapsrc).nsets
        return nabs+nbhd


    @staticmethod
    def _number_of_emission_lines(key, dapsrc=None):
        if key == 'None':
            return 0
        db = select_proc_method(key, EmissionLineModelDef,
                                available_func=available_emission_line_modeling_methods,
                                dapsrc=dapsrc)
        return EmissionLineDB(db['emission_lines'], dapsrc=dapsrc).nsets


    def _read(self):
        """
        Read the data in the existing file at :func:`file_path`.
        """
        # Close the HDUList if it's already open
        if self.hdu is not None:
            self.hdu.close()
            self.hdu = None

        # Read the file
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, 'Reading exiting file')
        self.hdu = fits.open(self.file_path())

        self.ndap = self.hdu['DAPALL'].header['NAXIS2']

        # Find the methods in the file
        self.methods = numpy.unique(self.hdu['DAPALL'].data['DAPTYPE'])

        # Find the plate-IFU combinations in the file
        pltifu = numpy.unique(self.hdu['DAPALL'].data['PLATEIFU'])
        self.plate, self.ifudesign = map(lambda x: numpy.array(x),
                                         numpy.array([ [int(p),int(f)] 
                                                        for p,f in pltifu.split('-')]).T.tolist())

        # Set the number of emission lines and spectral indices
        self.neml = self.hdu['DAPALL'].data['EMLINE_GFLUX_1RE'].shape[1]
        self.nindx = self.hdu['DAPALL'].data['SPECINDEX_1RE'].shape[1]

        # Report
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, 'Rows: {0}'.format(self.ndap))
            log_output(self.loggers, 1, logging.INFO, 'Available methods: {0}'.format(self.methods))
            log_output(self.loggers, 1, logging.INFO,
                       'Number of observations: {0}'.format(len(self.plate)))
            log_output(self.loggers, 1, logging.INFO,
                       'Maximum number of emission lines: {0}'.format(self.neml))
            log_output(self.loggers, 1, logging.INFO,
                       'Maximum number of spectral indices: {0}'.format(self.nindx))


    def _init_string_lengths(self):
        return { 'PLATEIFU':12, 
                 'MANGAID':10,
                 'MODE':4,
                 'DAPTYPE':20,
                 'VERSDRP2':8,
                 'VERSDRP3':8,
                 'VERSCORE':8,
                 'VERSUTIL':8,
                 'VERSDAP':8,
                 'RDXQAKEY':14,
                 'BINKEY':14,
                 'SCKEY':14,
                 'ELMKEY':14,
                 'ELFKEY':14,
                 'SIKEY':14,
                 'BINTYPE':10,
                 'TPLKEY':15,
                 'DATEDAP':10,
                 'EMLINE_NAME':20,
                 'SPECINDEX_NAME':20
                 'SPECINDEX_UNIT':5
               }

    
    def _table_dtype(self, neml, nindx):
        return [ ('PLATE', numpy.int),
                 ('IFUDESIGN', numpy.int),
                 ('PLATEIFU', '<U{0:d}'.format(self.str_len['PLATEIFU'])),
                 ('MANGAID', '<U{0:d}'.format(self.str_len['MANGAID'])),
                 ('DRPALLINDX', numpy.int),
                 ('MODE', '<U{0:d}'.format(self.str_len['MODE'])),
                 ('DAPTYPE', '<U{0:d}'.format(self.str_len['DAPTYPE'])),
                 ('DAPDONE', numpy.bool),
                 ('OBJRA', numpy.float),
                 ('OBJDEC', numpy.float),
                 ('IFURA', numpy.float),
                 ('IFUDEC', numpy.float),
                 ('MNGTARG1', numpy.int),
                 ('MNGTARG2', numpy.int),
                 ('MNGTARG3', numpy.int),
                 ('Z', numpy.float),
                 ('LDIST_Z', numpy.float),
                 ('ADIST_Z', numpy.float),
                 ('NSA_Z', numpy.float),
                 ('NSA_ZDIST', numpy.float),
                 ('LDIST_NSA_Z', numpy.float),
                 ('ADIST_NSA_Z', numpy.float),
                 ('NSA_ELPETRO_BA', numpy.float),
                 ('NSA_ELPETRO_PHI', numpy.float),
                 ('NSA_ELPETRO_TH50_R', numpy.float),
                 ('NSA_SERSIC_BA', numpy.float),
                 ('NSA_SERSIC_PHI', numpy.float),
                 ('NSA_SERSIC_TH50', numpy.float),
                 ('NSA_SERSIC_N', numpy.float),
                 ('VERSDRP2', '<U{0:d}'.format(self.str_len['VERSDRP2'])),
                 ('VERSDRP3', '<U{0:d}'.format(self.str_len['VERSDRP3'])),
                 ('VERSCORE', '<U{0:d}'.format(self.str_len['VERSCORE'])),
                 ('VERSUTIL', '<U{0:d}'.format(self.str_len['VERSUTIL'])),
                 ('VERSDAP', '<U{0:d}'.format(self.str_len['VERSDAP'])),
                 ('DRP3QUAL', numpy.int),
                 ('DAPQUAL', numpy.int),
                 ('RDXQAKEY', '<U{0:d}'.format(self.str_len['RDXQAKEY'])),
                 ('BINKEY', '<U{0:d}'.format(self.str_len['BINKEY'])),
                 ('SCKEY', '<U{0:d}'.format(self.str_len['SCKEY'])),
                 ('ELMKEY', '<U{0:d}'.format(self.str_len['ELMKEY'])),
                 ('ELFKEY', '<U{0:d}'.format(self.str_len['ELFKEY'])),
                 ('SIKEY', '<U{0:d}'.format(self.str_len['SIKEY'])),
                 ('BINTYPE', '<U{0:d}'.format(self.str_len['BINTYPE'])),
                 ('BINSNR', numpy.float),
                 ('TPLKEY', '<U{0:d}'.format(self.str_len['TPLKEY'])),
                 ('DATEDAP', '<U{0:d}'.format(self.str_len['DATEDAP'])),
                 ('DAPBINS', numpy.int),
                 ('RCOV90', numpy.float),
                 ('SNR_MED', numpy.float, (4,)),
                 ('SNR_RING', numpy.float, (4,)),
                 ('SB_1RE', numpy.float),
                 ('BIN_RMAX', numpy.float),
                 ('BIN_R_N', numpy.float, (3,)),
                 ('BIN_R_SNR', numpy.float, (3,)),
                 ('STELLAR_Z', numpy.float),
                 ('STELLAR_VEL_LO', numpy.float),
                 ('STELLAR_VEL_HI', numpy.float),
                 ('STELLAR_VEL_LO_CLIP', numpy.float),
                 ('STELLAR_VEL_HI_CLIP', numpy.float),
                 ('STELLAR_SIGMA_1RE', numpy.float),
                 ('STELLAR_CONT_RCHI2_1RE', numpy.float),
                 ('HA_Z', numpy.float),
                 ('HA_GVEL_LO', numpy.float),
                 ('HA_GVEL_HI', numpy.float),
                 ('HA_GVEL_LO_CLIP', numpy.float),
                 ('HA_GVEL_HI_CLIP', numpy.float),
                 ('HA_GSIGMA_1RE', numpy.float),
                 ('HA_GSIGMA_HI', numpy.float),
                 ('HA_GSIGMA_HI_CLIP', numpy.float),
                 ('EMLINE_NAME', '<U{0:d}'.format(self.str_len['EMLINE_NAME']*neml)),
                 ('EMLINE_GFLUX_CEN', numpy.float, (neml,)),        # NOT IMPLEMENTED
                 ('EMLINE_GFLUX_1RE', numpy.float, (neml,)),        # NOT IMPLEMENTED
                 ('EMLINE_GFLUX_TOT', numpy.float, (neml,)),
                 ('EMLINE_GSB_1RE', numpy.float, (neml,)),
                 ('EMLINE_GSB_PEAK', numpy.float, (neml,)),
                 ('EMLINE_GEW_1RE', numpy.float, (neml,)),
                 ('EMLINE_GEW_PEAK', numpy.float, (neml,)),         # NOT IMPLEMENTED
                 ('SPECINDEX_NAME', '<U{0:d}'.format(self.str_len['SPECINDEX_NAME']*nindx)),
                 ('SPECINDEX_UNIT', '<U{0:d}'.format(self.str_len['SPECINDEX_UNIT']*nindx)),
                 ('SPECINDEX_LO', numpy.float, (nindx,)),
                 ('SPECINDEX_HI', numpy.float, (nindx,)),
                 ('SPECINDEX_LO_CLIP', numpy.float, (nindx,)),
                 ('SPECINDEX_HI_CLIP', numpy.float, (nindx,)),
                 ('SPECINDEX_1RE', numpy.float, (nindx,)),         # NOT IMPLEMENTED
                 ('SFR_1RE', numpy.float)                       # NOT IMPLEMENTED
               ]


    def _get_completed_observations(self):
        """
        Get the list of DRP completed (and prospectively DAP analyzed)
        observations.

        .. todo:
            Save drpc as an attribute of the class?

        Returns:
            numpy.ndarray: Two arrays listing the plate and ifudesign
            numbers of the available observations that the DAP should
            have analyzed.

        Raises:
            FileNotFoundError: Raised if the DRPComplete file is not
                available.
        """
        # Initialize the DRP complete file (read-only)
        drpc = DRPComplete(readonly=True, drpver=self.drpver, redux_path=self.redux_path,
                           dapver=self.dapver, analysis_path=self.analysis_path)
        if not os.path.isfile(drpc.file_path()):
            raise FileNotFoundError('Could not access DRP complete file: {0}'.format(
                                                                                drpc.file_path()))
    
        # Find the list of file that should have been analyzed
        indx = (drpc['MANGAID'] != 'NULL') \
                & ((drpc['MANGA_TARGET1'] > 0) | (drpc['MANGA_TARGET3'] > 0)) \
                & (drpc['VEL'] > 0)

        return drpc['PLATE'][indx], drpc['IFUDESIGN'][indx]


    def _combine_plateifu_methods(self, plt, ifu):
        """
        Run the combinatorics to create the full list of plates,
        ifudesigns, and methods.
        """
        nmethods = len(self.methods)
        ncomplete = len(plt)
        methods = numpy.array([ numpy.array(self.methods) ]*ncomplete).T.ravel()
        plates = numpy.array([ numpy.array(plt) ]*nmethods).ravel() #.astype(str)
        ifudesigns = numpy.array([ numpy.array(ifu) ]*nmethods).ravel() #.astype(str)

        return methods, plates, ifudesigns


    def _construct_maps_file_list(self, plt, ifu):
        """
        Construct the list of MAPS files that should be in/added to the
        DAPall file.
        """
        # Construct the file names
        methods, plates, ifudesigns = self._combine_plateifu_methods(plt, ifu)
        return numpy.array( ['{0}/{1}/{2}/{3}/manga-{2}-{3}-MAPS-{1}.fits.gz'.format(
                                self.analysis_path, m, p, i) 
                                    for p,i,m in zip(plates, ifudesigns, methods) ])


    def _srd_snr_metric(self, plt, ifu, r_re_map):
        """
        Calculate the SRD S/N metric for a given plate-ifu in each of the griz bands.

        .. todo::
            Is the DRP LOGCUBE file needed for other quantities?
        """
        filter_response_file = [ os.path.join(self.dapsrc, 'data', 'filter_response',
                                               f) for f in [ 'gunn_2001_g_response.db',
                                                             'gunn_2001_r_response.db',
                                                             'gunn_2001_i_response.db',
                                                             'gunn_2001_z_response.db'] ]
        for f in filter_response_file:
            if not os.path.isfile(f):
                raise FileNotFoundError('{0} does not exist!'.format(f))

        # Read the DRPFits file
        drpf = DRPFits(plt, ifu, 'CUBE', drpver=self.drpver, redux_path=self.redux_path, read=True)
#        if not self.quiet:
#            log_output(loggers, 1, logging.INFO, 'Opened DRP file: {0}\n'.format(drpf.file_path()))

        # Put the map in the correct order for use with DRPFits and
        # unravel it
        r = r_re_map.T.ravel()
        # Get the spaxels withint the radius limits
        indx = numpy.arange(r.size)[(r > 1) & (r < 1.5)]
        # Get the appropriate covariance pixels to select
        ci, cj = map( lambda x: x.ravel(), numpy.meshgrid(indx, indx) )

        # Run the S/N calculation; this is the same calculation as done
        # in
        # mangadap.proc.spatialbinning.VoronoiBinning.sn_calculation_covariance_matrix
        nfilter = len(filter_response_file)
        snr_med = numpy.zeros(nfilter, dtype=float)
        snr_ring = numpy.zeros(nfilter, dtype=float)
        for i in range(nfilter):
            response_func = numpy.genfromtext(filter_response_file[i])[:,:2]
            signal, variance, snr, correl = drpf.flux_stats(response_func=response_func,
                                                            flag=['DONOTUSE', 'FORESTAR'],
                                                            covar=True, correlation=True)
            covar = correl.apply_new_variance(variance).toarray()
            snr_med[i] = numpy.median(snr[indx])
            snr_ring[i] = numpy.sum(signal[indx])/numpy.sqrt(numpy.sum(covar[ci,cj]))

        return snr_med, snr_ring

#        t = ((r > 1) & (r < 1.5)).astype(float)
#        sum_signal = numpy.dot(t, signal)
#        sum_variance = t.dot(covar._with_lower_triangle().cov.dot(t))


    @staticmethod
    def _radial_coverage_metric(r, ell):
        step = 1.0
        width = 2.5
        pixelscale = 0.5
        coverage_limit = 0.9

        maxr = numpy.ma.amax(r)
        binr = numpy.arange(width/2, maxr+width/4, width/4)
        coverage = numpy.zeros(len(binr), dtype=float)
        coverage[0] = numpy.sum(r < binr[0])/(numpy.pi*(1-ell)*numpy.square(binr[0]))
        for i in range(1,len(binr)):
            binrl = binr[i] - width/2
            binru = binrl + width
            ellipse_area = numpy.pi*(1-ell)*(numpy.square(binru) - numpy.square(binrl))
            coverage[i] = numpy.sum((r > binrl) & (r < binru))*numpy.square(pixelscale) \
                                / ellipse_area
        indx = numpy.arange(len(binr))[coverage < coverage_limit]
        interp = interpolate.interp1d(coverage[indx[-2]:], binr[indx[-2]:])
        return interp(coverage_limit)


    @staticmethod
    def _binning_metrics(dapmaps):
        """
        Return binning metrics:
            - Maximum radius of any valid bin
            - Number and median S/N of bins between 0-1, 0.5-1.5, and
              1.5-2.5 Re
        """
        # Unique bins
        # TODO: Select the BINID channel using the channel dictionary
        unique, indx = numpy.unique(dapmaps['BINID'].data[0,:,:].ravel(), return_index=True)
        unique = unique[1:]
        unique_indx = indx[1:]

        # Get the radii (in units of Re) of the unmasked, unique
        # bins constructed during the spatial binning
        r_re = dapmaps['BIN_LWELLCOO'].data.copy()[1,:,:].ravel()[unique_indx]

        # Get the S/N metric resulting from the spatial binning
        snr = dapmaps['BIN_SNR'].data.copy().ravel()[unique_indx]

        # Get the number of bins and median S/N between 0-1,
        # 0.5-1.5, and 1.5-2.5 Re
        rlim = numpy.array([ [0,1], [0.5,1.5], [1.5-2.5]])
        bin_r_n = numpy.empty(len(rlim), dtype=float)
        bin_r_snr = numpy.empty(len(rlim), dtype=float)
        for i in range(len(rlim)):
            indx = (r_re > rlim[i,0]) & (r_re < rlim[i,1])
            bin_r_n[i] = numpy.sum(indx)
            bin_r_snr[i] = numpy.median(snr[indx])
        
        return numpy.amax(r_re), bin_r_n, bin_r_snr


    def _stellar_kinematics_metrics(self, dapmaps, redshift):

        # Unique bins
        # TODO: Select the BINID channel using the channel dictionary
        unique, indx = numpy.unique(dapmaps['BINID'].data[1,:,:].ravel(), return_index=True)
        unique = unique[1:]
        unique_indx = indx[1:]

        # Pull the data from the maps file
        r_re = dapmaps['BIN_LWELLCOO'].data.copy()[1,:,:].ravel()[unique_indx]

        bin_flux = numpy.ma.MaskedArray(dapmaps['BIN_MFLUX'].data.copy(),
                                        mask=self.maps_bm.flagged(dapmaps['BIN_MFLUX_MASK'].data,
                                                                 'DONOTUSE')).ravel()[unique_indx]

        svel = numpy.ma.MaskedArray(dapmaps['STELLAR_VEL'].data.copy(),
                                    mask=self.maps_bm.flagged(dapmaps['STELLAR_VEL_MASK'].data,
                                                             'DONOTUSE')).ravel()[unique_indx]

        ssig = numpy.ma.MaskedArray(dapmaps['STELLAR_SIGMA'].data.copy(),
                                    mask=self.maps_bm.flagged(dapmaps['STELLAR_SIGMA_MASK'].data,
                                                              'DONOTUSE')).ravel()[unique_indx]

        ssig2corr = numpy.square(ssig) - numpy.square(
                                    dapmaps['STELLAR_SIGMACORR'].data[1,:,:].ravel()[unique_indx])
        scchi = dapmaps['STELLAR_CONT_RCHI2'].data.copy().ravel()[unique_indx]

        # Convert velocities to redshift
        z = (svel*(1+redshift) + astropy.constants.c.to('km/s').value*redshift) \
                                    / astropy.constants.c.to('km/s').value

        # Get the bins within an on-sky circular aperture of 2.5 arcsec
        d2 = numpy.sum(numpy.square(dapmaps['BIN_LWSKYCOO'].data), axis=0).ravel()[unique_indx]
        center = d2 < 1.25*1.25

        # Get the flux-weighted mean velocity at the center
        wsum = numpy.ma.sum(bin_flux[center])
        stellar_z = -999.0 if numpy.isclose(wsum, 0.) \
                                else numpy.ma.sum(bin_flux[center]*z[center])/wsum

        # Get the growth ranges of the stellar velocity
        stellar_vel_lo, stellar_vel_hi = sample_growth(svel.compressed(), [0.025,0.975])
        stellar_vel_lo_clip, stellar_vel_hi_clip = sample_growth(sigma_clip(svel).compressed()),
                                                                 [0.025,0.975])

        # Get the flux-weighted, corrected sigma within 1 Re
        within_1re = r_re < 1.
        wsum = numpy.ma.sum(bin_flux[within_1re])
        stellar_sigma2_1re = -999. if numpy.isclose(wsum, 0.) \
                                else numpy.ma.sum(bin_flux[within_1re]*ssig2corr[within_1re])/wsum
        stellar_sigma_1re = numpy.sqrt(stellar_sigma2_1re) if stellar_sigma2_1re > 0 else -999.

        # Get the median reduced chi-square
        stellar_cont_rchi2_1re = numpy.ma.median(scchi[within_1re])

        return stellar_z, stellar_vel_lo, stellar_vel_hi, \
                        stellar_vel_lo_clip, stellar_vel_hi_clip, stellar_sigma_1re, \
                        stellar_cont_rchi2_1re


    def _halpha_kinematics_metrics(self, dapmaps, redshift, spx_coo=False):

        # Unique bins
        # TODO: Select the BINID channel using the channel dictionary
        unique, indx = numpy.unique(dapmaps['BINID'].data[3,:,:].ravel(), return_index=True)
        unique = unique[1:]
        unique_indx = indx[1:]

        # Pull the data from the maps file
        cooext = 'SPX_ELLCOO' if spx_coo else 'BIN_LWELLCOO'
        r_re = dapmaps[cooext].data.copy()[1,:,:].ravel()[unique_indx]

        flux = numpy.ma.MaskedArray(dapmaps['EMLINE_GFLUX'].data.copy()[emline['Ha-6564'],:,:],
                                    mask=self.maps_bm.flagged(
                                        dapmaps['EMLINE_GFLUX_MASK'].data[emline['Ha-6564'],:,:],
                                                              'DONOTUSE')).ravel()[unique_indx]

        vel = numpy.ma.MaskedArray(dapmaps['EMLINE_GVEL'].data.copy()[emline['Ha-6564'],:,:],
                                   mask=self.maps_bm.flagged(
                                        dapmaps['EMLINE_GVEL_MASK'].data[emline['Ha-6564'],:,:],
                                                             'DONOTUSE')).ravel()[unique_indx]

        sig = numpy.ma.MaskedArray(dapmaps['EMLINE_GSIGMA'].data.copy()[emline['Ha-6564'],:,:],
                                   mask=self.maps_bm.flagged(
                                        dapmaps['EMLINE_GSIGMA_MASK'].data[emline['Ha-6564'],:,:],
                                                             'DONOTUSE')).ravel()[unique_indx]

        sig2corr = numpy.square(sig) - numpy.square(
                    dapmaps['EMLINE_INSTSIGMA'].data[emline['Ha-6564'],:,:].ravel()[unique_indx])


        # Convert velocities to redshift
        z = (vel*(1+redshift) + astropy.constants.c.to('km/s').value*redshift) \
                                    / astropy.constants.c.to('km/s').value

        # Get the bins within an on-sky circular aperture of 2.5 arcsec
        cooext = 'SPX_SKYCOO' if spx_coo else 'BIN_LWSKYCOO'
        d2 = numpy.sum(numpy.square(dapmaps[cooext].data), axis=0).ravel()[unique_indx]
        center = d2 < 1.25*1.25

        # Get the flux-weighted mean velocity at the center
        wsum = numpy.ma.sum(flux[center])
        halpha_z = -999.0 if numpy.isclose(wsum, 0.) \
                                else numpy.ma.sum(flux[center]*z[center])/wsum

        # Get the growth ranges of the H-alpha velocity
        halpha_gvel_lo, halpha_gvel_hi = sample_growth(vel.compressed(), [0.025,0.975])
        halpha_gvel_lo_clip, halpha_gvel_hi_clip = sample_growth(sigma_clip(vel).compressed()),
                                                                 [0.025,0.975])

        # Get the flux-weighted, corrected sigma within 1 Re
        within_1re = r_re < 1.
        wsum = numpy.ma.sum(flux[within_1re])
        halpha_gsigma2_1re = -999. if numpy.isclose(wsum, 0.) \
                                else numpy.ma.sum(flux[within_1re]*sig2corr[within_1re])/wsum
        halpha_gsigma_1re = numpy.sqrt(halpha_gsigma2_1re) if halpha_gsigma2_1re > 0 else -999.

        # Get the high-growth of the H-alpha velocity dispersion
        halpha_gsigma2_hi = sample_growth(sig2corr.compressed(), [0.975])
        halpha_gsigma2_hi_clip = sample_growth(sigma_clip(sig2corr).compressed()), [0.975])

        halpha_gsigma_hi = numpy.sqrt(halpha_gsigma2_hi) if halpha_gsigma2_hi > 0 else -999.
        halpha_gsigma_hi_clip = numpy.sqrt(halpha_gsigma2_hi_clip) \
                                        if halpha_gsigma2_hi_clip > 0 else -999.

        return halpha_z, halpha_gvel_lo, halpha_gvel_hi, halpha_gvel_lo_clip, halpha_gvel_hi_clip, \
                        halpha_gsigma_1re, halpha_gsigma_hi, halpha_gsigma_hi_clip


    def _emission_line_metrics(self, dapmaps, spx_coo=False):

        # Unique bins
        # TODO: Select the BINID channel using the channel dictionary
        unique, indx = numpy.unique(dapmaps['BINID'].data[3,:,:].ravel(), return_index=True)
        unique = unique[1:]
        unique_indx = indx[1:]

        neml = dapmaps['EMLINE_GFLUX'].shape[0]

        # Pull the data from the maps file
        cooext = 'SPX_ELLCOO' if spx_coo else 'BIN_LWELLCOO'
        r_re = dapmaps[cooext].data.copy()[1,:,:].ravel()[unique_indx]
        flux = numpy.ma.MaskedArray(dapmaps['EMLINE_GFLUX'].data.copy(),
                                    mask=self.maps_bm.flagged(dapmaps['EMLINE_GFLUX_MASK'].data,
                                                        'DONOTUSE')).reshape(neml,-1)[:,unique_indx]
        ew = numpy.ma.MaskedArray(dapmaps['EMLINE_GEW'].data.copy(),
                                  mask=self.maps_bm.flagged(dapmaps['EMLINE_GEW_MASK'].data,
                                                        'DONOTUSE')).reshape(neml,-1)[:,unique_indx]

        # Get the bins within an on-sky circular aperture of 2.5 arcsec
        cooext = 'SPX_SKYCOO' if spx_coo else 'BIN_LWSKYCOO'
        d2 = numpy.sum(numpy.square(dapmaps[cooext].data).reshape(neml,-1)[:,unique_indx], axis=0)
        center = d2 < 1.25*1.25
        within_1re = r_re < 1.

        # Get the total flux within a set of apertures
        emline_gflux_cen = numpy.ma.sum(flux[center], axis=1)
        emline_gflux_1re = numpy.ma.sum(flux[:,within_1re], axis=1)
        emline_gflux_tot = numpy.ma.sum(flux, axis=1)

        # Get the mean surface-brightness within 1 Re and the peak SB
        emline_gsb_1re = numpy.ma.mean(flux[:,within_1re], axis=1)
        emline_gsb_peak = numpy.ma.amax(flux, axis=1)

        # Get the mean equivalent width within 1 Re and the peak EW
        emline_gew_1re = numpy.ma.mean(ew[:,within_1re], axis=1)
        emline_gew_peak = numpy.ma.amax(ew, axis=1)

        return emline_gflux_cen, emline_gflux_1re, emline_gflux_tot, emline_gsb_1re, \
                        emline_gsb_peak, emline_gew_1re, emline_gew_peak


    def _spectral_index_metrics(self, dapmaps, specindex_units):

        # Unique bins
        # TODO: Select the BINID channel using the channel dictionary
        unique, indx = numpy.unique(dapmaps['BINID'].data[4,:,:].ravel(), return_index=True)
        unique = unique[1:]
        unique_indx = indx[1:]

        nindx = dapmaps['SPECINDEX'].shape[0]

        # Pull the data from the maps file
        r_re = dapmaps['BIN_LWELLCOO'].data.copy()[1,:,:].ravel()[unique_indx]
        specindex = numpy.ma.MaskedArray(dapmaps['SPECINDEX'].data.copy(),
                                         mask=self.maps_bm.flagged(dapmaps['SPECINDEX_MASK'].data,
                                                    'DONOTUSE')).reshape(nindx,-1)[:,unique_indx]
        specindex_corr = numpy.ma.MaskedArray(dapmaps['SPECINDEX_CORR'].data.copy(),
                                         mask=self.maps_bm.flagged(dapmaps['SPECINDEX_MASK'].data,
                                                    'DONOTUSE')).reshape(nindx,-1)[:,unique_indx]

        # Get the corrected indices
        ang = specindex_units == 'ang'
        mag = numpy.invert(ang)
        specindex_corr[ang] = specindex[ang]*specindex_corr[ang]
        specindex_corr[mag] = specindex[mag]+specindex_corr[mag]

        # Get the growth limits
        specindex_lo = (numpy.empty(nindx, dtype=float)
        specindex_hi = numpy.empty(nindx, dtype=float)
        specindex_lo_clip = numpy.empty(nindx, dtype=float)
        specindex_hi_clip = numpy.empty(nindx, dtype=float)
        for i in range(nindx): 
            specindex_lo[i], specindex_hi[i] = sample_growth(specindex_corr[i,:].compressed(),
                                                             [0.025,0.975])
            specindex_lo_clip[i] specindex_hi_clip[i] \
                        = sample_growth(sigma_clip(specindex_corr[i,:]).compressed(), [0.025,0.975])

        # Get the median within 1 Re
        specindex_1re = numpy.ma.median(specindex[:,r_re < 1.], axis=1)

        return specindex_lo, specindex_hi, specindex_lo_clip, specindex_hi_clip, specindex_1re


    def file_name(self):
        """Return the name of the DRP complete database."""
        return ('dapall-{0}-{1}.fits'.format(self.drpver, self.dapver))


    def file_path(self):
        """Return the full pat to the DRP complete database."""
        return os.path.join(self.analysis_path, self.file_name())


    def update(self, methods=None):
        """
        Update the DAPall file
        
        If clobber is True, the entire DAPall file is reconstructed from
        scratch.  Otherwise, any additional data is appended to the end
        of the file.

        Args:
            methods (str,list): (**Optional**) Specify a set of methods
                in the DAP plan file to include in the DAPall file.  If
                not provided, all methods in the plan file are included.

        Raises:
            ValueError: Raised if a provided method is not in the provided
                set of analysis plans.
        """

        # TODO: Counting the number of emission lines may need to allow
        # for a function that does not use an EmissionLineDB object.
        # Set the number of emission lines and spectral indices using a
        # dummy maps file for each plan?

        # Use the plans to set the number of emission lines for each
        # method
        self.neml = numpy.array([ self._number_of_emission_lines(p['elfit_key'],
                                                                 dapsrc=self.dapsrc)
                                    for p in self.plan ])

        # Use the plans to set the number of spectral indices
        self.nindx = numpy.array([ self._number_of_spectral_indices(p['spindex_key'],
                                                                    dapsrc=self.dapsrc)
                                    for p in self.plan ])

        # Get the list of methods to look for
        if methods is None:
            self.methods = self.plan_methods
        else:
            self.methods = numpy.atleast_1d(methods)
            if not numpy.all( [m in self.plan_methods for m in self.methods ]):
                raise ValueError('All provided methods must be in provided plan file.')

        # Report
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, '-'*50)
            log_output(self.loggers, 1, logging.INFO, 'RUNNING DAPALL UPDATE')
            log_output(self.loggers, 1, logging.INFO, '-'*50)
            log_output(self.loggers, 1, logging.INFO, 'Plan methods: {0}'.format(self.plan_methods))
            log_output(self.loggers, 1, logging.INFO,
                       'Emission lines in each plan: {0}'.format(self.neml))
            log_output(self.loggers, 1, logging.INFO,
                       'Spectral indices in each plan: {0}'.format(self.nindx))
            log_output(self.loggers, 1, logging.INFO,
                       'Methods for output: {0}'.format(self.methods))

        # Construct the list of files to look for
        self.plate, self.ifudesign = self._get_completed_observations()
        methodlist, platelist, ifulist = self._combine_plateifu_methods(self.plate, self.ifudesign)
        dapfiles = self._construct_maps_file_list(self.plate, self.ifudesign)

        # Report
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO,
                       'Completed observations: {0}'.format(len(self.plate)))
            log_output(self.loggers, 1, logging.INFO,
                       'Total number of expected DAP MAPS files: {0}'.format(len(dapfiles)))

#        # If the table already exists and not clobbering, remove any
#        # successfully completed files from the list
#        if not clobber and self.hdu is not None:
#
#            # Check that there's enough room for the emission lines
#            self['EMLINE_NAME'].shape
#
#            existing_files = self._construct_maps_file_list(self['PLATE'][self['DAPDONE']],
#                                                            self['IFUDESIGN'][self['DAPDONE']])
#            dapfiles = list( set(dapfiles) - set(existing_files) )
#            # TODO: Need to test this
#            emptydb = init_record_array(len(dapfiles),
#                                        self._table_dtype(max(self.neml), max(self.nindx)))

        # If clobbering and the HDUList is already initialized, close it
        # and reinitialize the hdu
        if clobber and self.hdu is not None:
            self.hdu.close()
            self.hdu = None
            self.ndap = None

        # Open the DRPall file
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO,
                       'Reading DRPall: {0}'.format(len(self.drpall_file)))
        hdu_drpall = fits.open(self.drpall_file)
        drpall = hdu_drpall[1].data

        # Initialize the output database
        self.ndap = len(dapfiles)
        max_neml = numpy.amax(self.neml)
        max_nindx = numpy.amax(self.nindx)
        db = init_record_array(self.ndap, self._table_dtype(max_neml, max_nindx))

        # Add the basic information
        db['PLATE'] = platelist
        db['IFUDESIGN'] = ifulist
        db['PLATEIFU'] = numpy.array([ '{0}-{1}'.format(p,i) for p,i in zip(platelist, ifulist) ])
        # For now only analyzing cubes!
        db['MODE'][:] = 'CUBE'
        db['DAPTYPE'] = methodlist

        # Copy from DRPall
        columns_from_drpall = [ 'MANGAID', 'OBJRA', 'OBJDEC', 'IFURA', 'IFUDEC', 'MNGTARG1',
                                'MNGTARG2', 'MNGTARG3', 'NSA_Z', 'NSA_ZDIST',
                                'NSA_ELPETRO_BA', 'NSA_ELPETRO_PHI', 'NSA_ELPETRO_TH50_R',
                                'NSA_SERSIC_BA', 'NSA_SERSIC_PHI', 'NSA_SERSIC_TH50',
                                'NSA_SERSIC_N', 'VERSDRP2', 'VERSDRP3', 'VERSCORE', 'VERSUTIL',
                                'DRP3QUAL' ]

        columns_from_maps_hdr = [ 'VERSDAP', 'DAPQUAL', 'RDXQAKEY', 'BINKEY', 'SCKEY', 'ELMKEY',
                                  'ELFKEY', 'SIKEY', 'BINTYPE' ]

        # Go through each file
        for i, f in enumerate(dapfiles):
            print('Processing {0}/{1}'.format(i+1,self.ndap))#, end='\r')

            # Find the index in the drpall file
            indx = numpy.where(drpall['PLATEIFU'] == db['PLATEIFU'][i])[0]
            # TODO: Set index to -1 instead of throwing an error?
            if len(indx) != 1:
                raise ValueError('Could not find {0} in DRPall file.'.format(db['PLATEIFU'][i]))
            db['DRPALLINDX'][i] = indx[0]

            # Add data from the DRPall file
            for c in columns_from_drpall:
                db[c][i] = drpall[c][db['DRPALLINDX'][i]]

            # Set the cosmological distances based on the NSA redshift
            db['LDIST_NSA_Z'][i] = self.cosmo.luminosity_distance(db['NSA_Z'][i]).value
            db['ADIST_NSA_Z'][i] = self.cosmo.angular_diameter_distance(db['NSA_Z'][i]).value

            # Open the maps file
            if not os.path.isfile(f):
                db['DAPDONE'][i] = 0    # A fault occurred such that the MAPS file was not produced
                continue
            db['DAPDONE'][i] = 1        # DAP successfully produced the file
            dapmaps = fits.open(f)

            # Add information from the DAP MAPS file header
            for c in columns_from_maps_hdr:
                db[c][i] = dapmaps['PRIMARY'].header[c]

            # Additional elements from the header that need to be handled
            # or have different names
            # TODO: REMOVE all of the latter to make sure keywords match
            try:
                db['BINSNR'][i] = dapmaps['PRIMARY'].header['BINSNR']
            except:
                db['BINSNR'][i] = -999.0
            db['TPLKEY'][i] = dapmaps['PRIMARY'].header['PPXFTPLK']
            db['DAPBINS'][i] = dapmaps['PRIMARY'].header['NBINS']
            try:
                db['DATEDAP'][i] = dapmaps['PRIMARY'].header['DATEDAP']
            except:
                db['DATEDAP'][i] = time.strftime('%Y-%m-%d',time.localtime(os.path.getmtime(f)))

            # TODO: Add any checks that the DRPall data matches the
            # header data in the MAPS files?

            # Basic mask for all maps: mask any spaxel not included
            # in any bin
            basic_map_mask = dapmaps['BINID'].data[0,:,:] < 0

            # Set input redshift and cosmological distances based on
            # this redshift
            db['Z'][i] = dapmaps['PRIMARY'].header['SCINPVEL']/astropy.constants.c.to('km/s').value
            db['LDIST_Z'][i] = self.cosmo.luminosity_distance(db['Z'][i]).value
            db['ADIST_Z'][i] = self.cosmo.angular_diameter_distance(db['Z'][i]).value

            # Get the dictionary with the channel names
            emline = channel_dictionary(dapmaps, 'EMLINE_GFLUX')
            neml = len(emline)
            specindex_names = channel_dictionary(dapmaps, 'SPECINDEX')
            specindex_units = channel_units(dapmaps, 'SPECINDEX')
            nindx = len(specindex_names)

            # Set the channel names and units to the table; the sorting is
            # required to make sure the names are put in the correct
            # order.
            srt = numpy.argsort(list(emline.values()))
            db['EMLINE_NAME'][i] = ''.join([ '{0}'.format(n).rjust(self.str_len['EMLINE_NAME'])
                                                for n in numpy.array(list(emline.keys()))[srt] ])
            srt = numpy.argsort(list(specindex_names.values()))
            db['SPECINDEX_NAME'][i] \
                                = ''.join([ '{0}'.format(n).rjust(self.str_len['SPECINDEX_NAME'])
                                        for n in numpy.array(list(specindex_names.keys()))[srt] ])
            db['SPECINDEX_UNIT'][i] \
                                = ''.join([ '{0}'.format(n).rjust(self.str_len['SPECINDEX_UNIT'])
                                            for n in specindex_units ])

            # Determine the coverage of the datacube (independent of the
            # DAP method)
            r = numpy.ma.MaskedArray(dapmaps['SPX_ELLCOO'].data.copy()[0,:,:], mask=basic_map_mask)
            db['RCOV90'][i] = _radial_coverage_metric(r.compressed(), 1-db['NSA_ELPETRO_BA'][i])

            # Get the SRD-based S/N metric (independent of the DAP
            # method)
            r_re = numpy.ma.MaskedArray(dapmaps['SPX_ELLCOO'].data.copy()[1,:,:],
                                        mask=numpy.invert(dapmaps['SPX_SNR'] > 0))
            db['SNR_MED'][i], db['SNR_RING'][i] = _srd_snr_metric(db['PLATE'][i],
                                                                  db['IFUDESIGN'][i], r_re):

            # Get the mean surface brightness within 1 Re used by the
            # stellar-continuum fitting (independent of the DAP method)
            mflux = numpy.ma.MaskedArray(dapmaps['SPX_MFLUX'].data.copy(),
                                        mask=numpy.invert(dapmaps['SPX_SNR'] > 0))
            db['SB_1RE'][i] = numpy.ma.mean(mflux[r_re < 1])

            # Get the spatial binning metrics
            db['BIN_RMAX'][i], db['BIN_R_N'][i], db['BIN_R_SNR'][i] = _binning_metrics(dapmaps)

            # Get the metrics for the stellar kinematics
            db['STELLAR_Z'][i], db['STELLAR_VEL_LO'][i], db['STELLAR_VEL_HI'][i], \
                    db['STELLAR_VEL_LO_CLIP'][i], db['STELLAR_VEL_HI_CLIP'][i], \
                    db['STELLAR_SIGMA_1RE'][i], db['STELLAR_CONT_RCHI2_1RE'] \
                            = self._stellar_kinematics_metrics(dapmaps, db['Z'][i])

            # Get the metrics for the emission-line kinematics
            # TODO: Check for hybrid binning method, then use
            # spx_coo=True if emission lines measured on individual
            # spaxels
            db['HA_Z'][i], db['HA_GVEL_LO'][i], db['HA_GVEL_HI'][i], db['HA_GVEL_LO_CLIP'][i], \
                    db['HA_GVEL_HI_CLIP'][i], db['HA_GSIGMA_1RE'][i], db['HA_GSIGMA_HI'][i], \
                    db['HA_GSIGMA_HI_CLIP'][i] = self._halpha_kinematics_metrics(dapmaps,
                                                                                 db['Z'][i])

            # Get the emission-line metrics
            # TODO: Check for hybrid binning method, then use
            # spx_coo=True if emission lines measured on individual
            # spaxels
            db['EMLINE_GFLUX_CEN'][i], db['EMLINE_GFLUX_1RE'][i], db['EMLINE_GFLUX_TOT'][i], \
                    db['EMLINE_GSB_1RE'][i], db['EMLINE_GSB_PEAK'][i], db['EMLINE_GEW_1RE'][i], \
                    db['EMLINE_GEW_PEAK'][i] = self._emission_line_metrics(dapmaps)

            # Get the spectral-index metrics
            db['SPECINDEX_LO'][i], db['SPECINDEX_HI'][i], db['SPECINDEX_LO_CLIP'][i], \
                    db['SPECINDEX_HI_CLIP'][i], db['SPECINDEX_1RE'][i] \
                            = self._spectral_index_metrics(dapmaps, specindex_units)

            # Estimate the star-formation rate
            log_Mpc_in_cm = numpy.log10(astropy.constants.pc.to('cm').value) + 6
            log_halpha_luminosity_1re = numpy.log10(4*numpy.pi) \
                                            + numpy.log10(db['EMLINE_GFLUX_1RE'][i]) - 17 \
                                            + 2*numpy.log10(db['LDIST_Z'][i]) + 2*log_Mpc_in_cm
            db['SFR_1RE'][i] = numpy.power(10, log_halpha_luminosity_1re - 41.27)

        print('Processing {0}/{0}'.format(self.ndap))

        # Create the primary header
        hdr = fits.Header()
        hdr['DATE'] = (time.strftime('%Y-%m-%d',time.gmtime()), 'UTC date created')
        hdr['VERSDRP3'] = (self.drpver, 'DRP version')
        hdr['VERSDAP'] = (self.dapver, 'DAP version')

        # Set the fits HDU list
        self.hdu = fits.HDUList([ fits.PrimaryHDU(header=hdr),
                                  fits.BinTableHDU.from_columns(
                                    [ fits.Column(name=n, format=rec_to_fits_type(db[n]),
                                      array=db[n]) for n in db.dtype.names ],
                                    name='DAPALL') ])

        # Write the file (this function is probably overkill)
        DAPFitsUtil.write(self.hdu, self.file_path(), clobber=True, checksum=True,
                          loggers=self.loggers)
        

