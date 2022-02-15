# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
r"""
Defines the class that constructs the DAPall summary table.

This is a post processing script that **must** be run after the DRPall
file is created.

.. todo::

    - Include columns with just the absolute luminosity of the
      H-alpha line as the intermediate step between the apparent
      luminosity and the star-formation rate.
    - Provide a rough Balmer decrement measurement?

----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

import logging
import time
import os
import warnings
import numpy

from scipy import interpolate

from astropy.io import fits
import astropy.constants
from astropy.cosmology import FlatLambdaCDM
import astropy.units
from astropy.stats import sigma_clip

from .drpcomplete import DRPComplete
from ..dapfits import DAPMapsBitMask, DAPQualityBitMask
from ..util.drpfits import DRPQuality3DBitMask
from ..config import defaults
from ..util.fileio import init_record_array, rec_to_fits_type, channel_dictionary, channel_units
from ..util.fitsutil import DAPFitsUtil
from ..util.log import init_DAP_logging, module_logging, log_output
from ..par.analysisplan import AnalysisPlanSet
from ..par.absorptionindexdb import AbsorptionIndexDB
from ..par.bandheadindexdb import BandheadIndexDB
from ..par.emissionlinedb import EmissionLineDB
from ..par.emissionmomentsdb import EmissionMomentsDB
from ..proc.util import select_proc_method, sample_growth
from ..proc.spatiallybinnedspectra import SpatiallyBinnedSpectra
from ..proc.stellarcontinuummodel import StellarContinuumModel
from ..proc.emissionlinemoments import EmissionLineMomentsDef
from ..proc.emissionlinemoments import available_emission_line_moment_databases
from ..proc.emissionlinemodel import EmissionLineModel, EmissionLineModelDef
from ..proc.emissionlinemodel import available_emission_line_modeling_methods
from ..proc.spectralindices import SpectralIndicesDef, available_spectral_index_databases

from matplotlib import pyplot

class DAPall:
    """
    Construct the summary table information for the DAP.

    Any observation in the
    :class:`mangadap.survey.drpcomplete.DRPComplete` file that
    satisfies the criteria set by
    :func:`mangadap.survey.drpcomplete.DRPComplete.can_analyze` is
    considered. If the expected output files exist, they are analyzed
    and included in the DAPall database. If the expected files do not
    exist, the construction of the DAPall (and the relevant flag)
    assumes that the DAP analysis of this cube failed.

    Args:
        plan (:class:`mangadap.par.analysisplan.AnalysisPlanSet`):
            The plan object used by the DAP.
        methods (:obj:`str`, :obj:`list`, optional):
            Specify a set of methods in the DAP plan file to include
            in the DAPall file. If not provided, all methods in the
            plan file are included. Each method produces an extension
            in the DAPall file, where the order of the extensions is
            alphanumeric.
        drpver (:obj:`str`, optional):
            DRP version.  Default determined by
            :func:`mangadap.config.defaults.drp_version`.
        redux_path (:obj:`str`, optional):
            Top-level directory with the DRP products; default is
            defined by
            :func:`mangadap.config.defaults.drp_redux_path`.
        dapver (:obj:`str`, optional):
            DAP version.  Default determined by
            :func:`mangadap.config.defaults.dap_version`.
        analysis_path (:obj:`str`, optional):
            Top-level directory for the DAP output data; default is
            defined by
            :func:`mangadap.config.defaults.dap_analysis_path`.
        readonly (:obj:`bool`, optional):
            If it exists, open any existing file and disallow any
            modifications of the database. Default is to
            automatically build the entire database if it doesn't
            exist.
        loggers (:obj:`list`, optional):
            List of `logging.Logger`_ objects used to log progress.
        quiet (:obj:`bool`, optional):
            Suppress all terminal and logging output.

    Raises:
        TypeError: Raised if the input plan is not a
            :class:`mangadap.par.analysisplan.AnalysisPlanSet` object.
        FileNotFoundError: Raise if the DRPall or the DRPComplete file
            are not available.
    
    Attributes:
        drpver (:obj:`str`):
            DRP version
        redux_path (:obj:`str`):
            Path to top-level redux directory
        dapver (:obj:`str`):
            DAP version
        analysis_path (:obj:`str`):
            Path to top-level analysis directory
        drpall_file (:obj:`str`):
            Path to the DRPall file
        h (:obj:`float`):
            The Hubble parameter (currently always set to unity)
        H0 (:obj:`float`):
            Hubbles constant (km/s/Mpc)
        cosmo (`astropy.cosmology.FlatLambdaCDM`_):
            Object used for distance calculations.
        maps_bm (:class:`mangadap.dapfits.DAPMapsBitMask`):
            Bitmask used to mask map pixels.
        neml (:obj:`int`):
            Number of emission lines included in each method.
        nindx (:obj:`int`):
            Number of spectral indices included in each method.
        str_len (:obj:`dict`):
            Dictionary with the number of characters used for each
            string entry in the database.
        hdu (`astropy.io.fits.HDUList`_):
            Object with the table data.
        ndap (:obj:`int`):
            Number of rows in the data table; equivalent to the
            number of processed MAPS files.
        plate (`numpy.ndarray`_):
            Array of the plate numbers
        ifudesign (`numpy.ndarray`_):
            Array of the ifudesigns
        methods (`numpy.ndarray`_):
            Array of the methods included in the DAPall file.
        readonly (:obj:`bool`):
            Object is prohibited from updating the database.
        loggers (:obj:`list`):
            List of `logging.Logger`_ objects to log progress;
            ignored if quiet=True.
        quiet (:obj:`bool`):
            Suppress all terminal and logging output.
        float_dtype (:obj:`str`):
            Sets precision for floating-point numbers.
    """
    def __init__(self, plan, methods=None, drpver=None, redux_path=None, dapver=None,
                 analysis_path=None, readonly=False, loggers=None, quiet=False,
                 single_precision=False):

        # Initialize the reporting
        self.loggers = loggers
        self.quiet = quiet

        self.float_dtype = 'float32' if single_precision else 'float'

        # Path definitions
        self.drpver = defaults.drp_version() if drpver is None else str(drpver)
        self.redux_path = defaults.drp_redux_path(self.drpver) \
                                if redux_path is None else str(redux_path)

        self.dapver = defaults.dap_version() if dapver is None else str(dapver)
        self.analysis_path = defaults.dap_analysis_path(self.drpver, self.dapver) \
                                if analysis_path is None else str(analysis_path)

        # Set the name of the DRPall file
        self.drpall_file = os.path.join(self.redux_path, 'drpall-{0}.fits'.format(self.drpver))
        if not os.path.isfile(self.drpall_file):
            raise FileNotFoundError('{0} does not exist!'.format(self.drpall_file))

        # Initialize the bitmask
        self.maps_bm = DAPMapsBitMask()

        # Set the length for the string elements of the table
        self.str_len = self._init_string_lengths()

        # Initialize the cosmology object
        self.h = 1.0
        self.H0 = 100 * self.h * astropy.units.km / astropy.units.s / astropy.units.Mpc
        self.cosmo = FlatLambdaCDM(H0=self.H0, Om0=0.3)

        # Declare the other attributes for use later on
        self.nmom = None
        self.elmom_channels = None
        self.neml = None
        self.elfit_channels = None
        self.nindx = None
        self.spindx_channels = None
        self.spindx_units = None
        self.hdu = None
        self.ndap = None
        self.plate = None
        self.ifudesign = None
        self.methods = None

        # Give an initial report
        self.readonly = readonly
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, '-'*50)
            log_output(self.loggers, 1, logging.INFO, '{0:^50}'.format('DAPALL DATABASE'))
            log_output(self.loggers, 1, logging.INFO, '-'*50)
            log_output(self.loggers, 1, logging.INFO, 'DRPall file: {0}'.format(self.drpall_file))
            log_output(self.loggers, 1, logging.INFO, 'Output path: {0}'.format(self.analysis_path))
            log_output(self.loggers, 1, logging.INFO, 'Output file: {0}'.format(self.file_name()))
            if self.readonly:
                log_output(self.loggers, 1, logging.INFO, 'READ-ONLY MODE')

        # If the output file already exists, read it
        if os.path.isfile(self.file_path()):
            self._read()

        # Check the combination of read-only and the read data makes
        # sense
        if self.readonly and self.hdu is None:
            raise ValueError('Selected read-only, but no database can be read.')
        if not self.readonly and self.hdu is not None:
            warnings.warn('DAPall file already exists!  Update will overwrite existing data.')

        # Automatically attempt to update the database if no database
        # exists or if the database is not meant to be read-only
        if self.hdu is None or not self.readonly:
            self.update(plan, methods=methods)

#    def __getitem__(self, key):
#        return self.hdu['DAPALL'].data[key]

    @staticmethod
    def _emission_line_moment_db_info(key):
        if key == 'None':
            return 'None', None
        db = select_proc_method(key, EmissionLineMomentsDef,
                                available_func=available_emission_line_moment_databases)
        return db['passbands'], EmissionMomentsDB.from_key(db['passbands']).channel_names()

    @staticmethod
    def _emission_line_db_info(key):
        if key == 'None':
            return 'None', None
        db = select_proc_method(key, EmissionLineModelDef,
                                available_func=available_emission_line_modeling_methods)
        return db['emission_lines'], EmissionLineDB.from_key(db['emission_lines']).channel_names()

    @staticmethod
    def _spectral_index_db_info(key):
        if key == 'None':
            return 0, 'None', 'None', None
        db = select_proc_method(key, SpectralIndicesDef,
                                available_func=available_spectral_index_databases)
        absdb = AbsorptionIndexDB.from_key(db['absindex'])
        bhddb = BandheadIndexDB.from_key(db['bandhead'])
        nindx = absdb.size + bhddb.size
        units = numpy.array(absdb['units'].tolist() + ['']*bhddb.size)
        channels = absdb.channel_names()
        channels.update(bhddb.channel_names(offset=absdb.size))
        if len(channels) != nindx:
            raise ValueError('Spectral index channel names not unique!')
        return db['absindex'], db['bandhead'], channels, units

    # TODO: Fix this!
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

        # Methods are the extension names
        self.methods = [h.name for h in self.hdu][1:]

        # Number of cubes analyzed by the DAP
        self.ndap = self.hdu[self.methods[0]].header['NAXIS2']

        # Find the plate-IFU combinations in the file
        pltifu = numpy.unique(self.hdu[self.methods[0]].data['PLATEIFU'])
        self.plate, self.ifudesign = map(lambda x: numpy.array(x),
                                         numpy.array([ [int(p),int(f)] 
                                                        for p,f in pltifu.split('-')]).T.tolist())

        # Set the number of emission lines and spectral indices
        # TODO: Get this from the number of header keywords?
        self.nmom = self.hdu[self.methods[0]].data['EMLINE_SFLUX_1RE'].shape[1]
        self.neml = self.hdu[self.methods[0]].data['EMLINE_GFLUX_1RE'].shape[1]
        self.nindx = self.hdu[self.methods[0]].data['SPECINDEX_1RE'].shape[1]

        # Report
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, 'Rows per method: {0}'.format(self.ndap))
            log_output(self.loggers, 1, logging.INFO, 'Methods: {0}'.format(self.methods))
            log_output(self.loggers, 1, logging.INFO,
                       'Number of observations: {0}'.format(len(self.plate)))
            log_output(self.loggers, 1, logging.INFO,
                       'Number of emission-line moments: {0}'.format(self.nmom))
            log_output(self.loggers, 1, logging.INFO,
                       'Number of emission lines: {0}'.format(self.neml))
            log_output(self.loggers, 1, logging.INFO,
                       'Number of spectral indices: {0}'.format(self.nindx))

    def _init_string_lengths(self):
        return { 'PLATEIFU':12, 
                 'MANGAID':10,
                 'MODE':4,
                 'DAPTYPE':30,
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
                 'EMLMOM_NAME':20,
                 'EMLINE_NAME':20,
                 'SPECINDEX_NAME':20,
                 'SPECINDEX_UNIT':5
               }

    def _table_dtype(self, nmom, neml, nindx):
        return [ ('PLATE', numpy.int),
                 ('IFUDESIGN', numpy.int),
                 ('PLATEIFU', '<U{0:d}'.format(self.str_len['PLATEIFU'])),
                 ('MANGAID', '<U{0:d}'.format(self.str_len['MANGAID'])),
                 ('DRPALLINDX', numpy.int),
                 ('MODE', '<U{0:d}'.format(self.str_len['MODE'])),
                 ('DAPTYPE', '<U{0:d}'.format(self.str_len['DAPTYPE'])),
                 ('DAPDONE', numpy.bool),
                 ('OBJRA', self.float_dtype),
                 ('OBJDEC', self.float_dtype),
                 ('IFURA', self.float_dtype),
                 ('IFUDEC', self.float_dtype),
                 ('MNGTARG1', numpy.int),
                 ('MNGTARG2', numpy.int),
                 ('MNGTARG3', numpy.int),
                 ('Z', self.float_dtype),
                 ('LDIST_Z', self.float_dtype),
                 ('ADIST_Z', self.float_dtype),
                 ('NSA_Z', self.float_dtype),
                 ('NSA_ZDIST', self.float_dtype),
                 ('LDIST_NSA_Z', self.float_dtype),
                 ('ADIST_NSA_Z', self.float_dtype),
                 ('NSA_ELPETRO_BA', self.float_dtype),
                 ('NSA_ELPETRO_PHI', self.float_dtype),
                 ('NSA_ELPETRO_TH50_R', self.float_dtype),
                 ('NSA_SERSIC_BA', self.float_dtype),
                 ('NSA_SERSIC_PHI', self.float_dtype),
                 ('NSA_SERSIC_TH50', self.float_dtype),
                 ('NSA_SERSIC_N', self.float_dtype),
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
                 ('BINSNR', self.float_dtype),
                 ('TPLKEY', '<U{0:d}'.format(self.str_len['TPLKEY'])),
                 ('DATEDAP', '<U{0:d}'.format(self.str_len['DATEDAP'])),
                 ('DAPBINS', numpy.int),
                 ('RCOV90', self.float_dtype),
                 ('SNR_MED', self.float_dtype, (4,)),
                 ('SNR_RING', self.float_dtype, (4,)),
                 ('SB_1RE', self.float_dtype),
                 ('BIN_RMAX', self.float_dtype),
                 ('BIN_R_N', self.float_dtype, (3,)),
                 ('BIN_R_SNR', self.float_dtype, (3,)),
                 ('STELLAR_Z', self.float_dtype),
                 ('STELLAR_VEL_LO', self.float_dtype),
                 ('STELLAR_VEL_HI', self.float_dtype),
                 ('STELLAR_VEL_LO_CLIP', self.float_dtype),
                 ('STELLAR_VEL_HI_CLIP', self.float_dtype),
                 ('STELLAR_SIGMA_1RE', self.float_dtype),
                 ('STELLAR_RCHI2_1RE', self.float_dtype),
                 ('HA_Z', self.float_dtype),
                 ('HA_GVEL_LO', self.float_dtype),
                 ('HA_GVEL_HI', self.float_dtype),
                 ('HA_GVEL_LO_CLIP', self.float_dtype),
                 ('HA_GVEL_HI_CLIP', self.float_dtype),
                 ('HA_GSIGMA_1RE', self.float_dtype),
                 ('HA_GSIGMA_HI', self.float_dtype),
                 ('HA_GSIGMA_HI_CLIP', self.float_dtype),
                 ('EMLINE_RCHI2_1RE', self.float_dtype),
#                 ('EMLMOM_NAME', '<U{0:d}'.format(self.str_len['EMLMOM_NAME']*nmom)),
                 ('EMLINE_SFLUX_CEN', self.float_dtype, (nmom,)),
                 ('EMLINE_SFLUX_1RE', self.float_dtype, (nmom,)),
                 ('EMLINE_SFLUX_TOT', self.float_dtype, (nmom,)),
                 ('EMLINE_SSB_1RE', self.float_dtype, (nmom,)),
                 ('EMLINE_SSB_PEAK', self.float_dtype, (nmom,)),
                 ('EMLINE_SEW_1RE', self.float_dtype, (nmom,)),
                 ('EMLINE_SEW_PEAK', self.float_dtype, (nmom,)),
#                 ('EMLINE_NAME', '<U{0:d}'.format(self.str_len['EMLINE_NAME']*neml)),
                 ('EMLINE_GFLUX_CEN', self.float_dtype, (neml,)),
                 ('EMLINE_GFLUX_1RE', self.float_dtype, (neml,)),
                 ('EMLINE_GFLUX_TOT', self.float_dtype, (neml,)),
                 ('EMLINE_GSB_1RE', self.float_dtype, (neml,)),
                 ('EMLINE_GSB_PEAK', self.float_dtype, (neml,)),
                 ('EMLINE_GEW_1RE', self.float_dtype, (neml,)),
                 ('EMLINE_GEW_PEAK', self.float_dtype, (neml,)),
#                 ('SPECINDEX_NAME', '<U{0:d}'.format(self.str_len['SPECINDEX_NAME']*nindx)),
#                 ('SPECINDEX_UNIT', '<U{0:d}'.format(self.str_len['SPECINDEX_UNIT']*nindx)),
                 ('SPECINDEX_LO', self.float_dtype, (nindx,)),
                 ('SPECINDEX_HI', self.float_dtype, (nindx,)),
                 ('SPECINDEX_LO_CLIP', self.float_dtype, (nindx,)),
                 ('SPECINDEX_HI_CLIP', self.float_dtype, (nindx,)),
                 ('SPECINDEX_1RE', self.float_dtype, (nindx,)),
                 ('SFR_1RE', self.float_dtype),
                 ('SFR_TOT', self.float_dtype)
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
    
        # Find the list of files that should have been analyzed
        indx = drpc.can_analyze()
        return drpc['PLATE'][indx], drpc['IFUDESIGN'][indx]

#    def _combine_plateifu_methods(self, plt, ifu):
#        """
#        Run the combinatorics to create the full list of plates,
#        ifudesigns, and methods.
#        """
#        nmethods = len(self.methods)
#        ncomplete = len(plt)
#        methods = numpy.array([ numpy.array(self.methods) ]*ncomplete).T.ravel()
#        plates = numpy.array([ numpy.array(plt) ]*nmethods).ravel() #.astype(str)
#        ifudesigns = numpy.array([ numpy.array(ifu) ]*nmethods).ravel() #.astype(str)
#        return methods, plates, ifudesigns

    def _construct_maps_file_list(self, method):
        """
        Construct the list of MAPS files that should be in/added to the
        DAPall file for a given method.
        """
        # Construct the file names
        return numpy.array(['{0}/{1}/{2}/{3}/manga-{2}-{3}-MAPS-{1}.fits.gz'.format(
                                self.analysis_path, method, p, i) 
                                    for p,i in zip(self.plate, self.ifudesign)])
    
    def _add_channels_to_header(self, hdr, channels, prefix, comment, units=None):
        """
        units is supposed to be a list or numpy array
        """
        ch = numpy.array(list(channels.values())).astype(int)
        srt = numpy.argsort(ch)
        ch = ch[srt]
        keys = numpy.array(list(channels.keys()))[srt]
        ndig = int(numpy.log10(ch[-1]))+1
        for c,k in zip(ch,keys):
            hdr[prefix+'{0}'.format(c+1).zfill(ndig)] = (k, comment)
            if units is not None:
                hdr[prefix+'U'+'{0}'.format(c+1).zfill(ndig)] = (units[c], comment+' unit')
        return hdr

    def _srd_snr_metric(self, dapmaps):
        """Grab the SNR metrics from the MAPS headers."""
        filters = [ 'G', 'R', 'I', 'Z' ]
        return numpy.array([dapmaps['PRIMARY'].header['SNR{0}MED'.format(f)] for f in filters]), \
                numpy.array([dapmaps['PRIMARY'].header['SNR{0}RING'.format(f)] for f in filters])

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

        indx = numpy.arange(len(binr))[coverage > coverage_limit]
        if len(indx) == 0:
            warnings.warn('No ring had larger than 90% coverage!')
            return 0.
        interp = interpolate.interp1d(coverage[indx[-1]:], binr[indx[-1]:])
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
        if unique[0] == -1:
            unique = unique[1:]
            unique_indx = indx[1:]

        # Get the radii (in units of Re) of the unmasked, unique
        # bins constructed during the spatial binning
        r_re = dapmaps['BIN_LWELLCOO'].data.copy()[1,:,:].ravel()[unique_indx]

        # Get the S/N metric resulting from the spatial binning
        snr = dapmaps['BIN_SNR'].data.copy().ravel()[unique_indx]

        # Get the number of bins and median S/N between 0-1,
        # 0.5-1.5, and 1.5-2.5 Re
        rlim = numpy.array([ [0,1], [0.5,1.5], [1.5,2.5]])
        bin_r_n = numpy.empty(len(rlim), dtype=float)
        bin_r_snr = numpy.empty(len(rlim), dtype=float)
        for i in range(len(rlim)):
            indx = (r_re > rlim[i,0]) & (r_re < rlim[i,1])
            if numpy.sum(indx) == 0:
                bin_r_n[i] = 0
                bin_r_snr[i] = 0
                continue
            bin_r_n[i] = numpy.sum(indx)
            bin_r_snr[i] = numpy.median(snr[indx])
        
        return numpy.amax(r_re), bin_r_n, bin_r_snr

    def _stellar_kinematics_metrics(self, dapmaps, redshift):

        # Unique bins
        # TODO: Select the BINID channel using the channel dictionary
        unique, indx = numpy.unique(dapmaps['BINID'].data[1,:,:].ravel(), return_index=True)
        if unique[0] == -1:
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
#                                    dapmaps['STELLAR_SIGMACORR'].data.ravel()[unique_indx])
                                    dapmaps['STELLAR_SIGMACORR'].data[0,:,:].ravel()[unique_indx])

        scchi = dapmaps['STELLAR_FOM'].data[2,:,:].copy().ravel()[unique_indx]

        # Convert velocities to redshift
        z = (svel*(1+redshift) + astropy.constants.c.to('km/s').value*redshift) \
                                    / astropy.constants.c.to('km/s').value

        # Get the bins within an on-sky circular aperture of 2.5 arcsec
        d2 = numpy.sum(numpy.square(dapmaps['BIN_LWSKYCOO'].data), axis=0).ravel()[unique_indx]
        center = d2 < 1.25*1.25

        # Get the flux-weighted mean velocity at the center
        wsum = numpy.ma.sum(bin_flux[center])
        stellar_z = -999. if numpy.all(bin_flux.mask[center] | z.mask[center]) \
                                or numpy.isclose(wsum, 0.) \
                            else numpy.ma.sum(bin_flux[center]*z[center])/wsum
        if not numpy.isfinite(stellar_z):
            raise ValueError('Stellar_z is not finite: {0}'.format(stellar_z))

        # Get the growth ranges of the stellar velocity
        stellar_vel_lo, stellar_vel_hi = sample_growth(svel.compressed(), [0.025,0.975])
        stellar_vel_lo_clip, stellar_vel_hi_clip = sample_growth(sigma_clip(svel).compressed(),
                                                                 [0.025,0.975])

        # Get the flux-weighted, corrected sigma within 1 Re
        within_1re = r_re < 1.
        wsum = numpy.ma.sum(bin_flux[within_1re])
        stellar_sigma2_1re = -999. if numpy.all(bin_flux.mask[within_1re]) \
                                            or numpy.isclose(wsum, 0.) \
                                else numpy.ma.sum(bin_flux[within_1re]*ssig2corr[within_1re])/wsum
        stellar_sigma_1re = numpy.sqrt(stellar_sigma2_1re) if stellar_sigma2_1re > 0 else -999.

        # Get the median reduced chi-square
        stellar_rchi2_1re = -999. if numpy.sum(within_1re) == 0 \
                                        else numpy.median(scchi[within_1re])

        return stellar_z, stellar_vel_lo, stellar_vel_hi, \
                        stellar_vel_lo_clip, stellar_vel_hi_clip, stellar_sigma_1re, \
                        stellar_rchi2_1re

    def _halpha_kinematics_metrics(self, dapmaps, deconstruct, redshift):

        # Unique bins
        # TODO: Select the BINID channel using the channel dictionary
        unique, indx = numpy.unique(dapmaps['BINID'].data[3,:,:].ravel(), return_index=True)
        if unique[0] == -1:
            unique = unique[1:]
            unique_indx = indx[1:]

        # Channel dictionary
        emline = channel_dictionary(dapmaps, 'EMLINE_GFLUX')

        # Pull the data from the maps file
        cooext = 'SPX_ELLCOO' if deconstruct else 'BIN_LWELLCOO'
        r_re = dapmaps[cooext].data.copy()[1,:,:].ravel()[unique_indx]

        try:
            flux = numpy.ma.MaskedArray(dapmaps['EMLINE_GFLUX'].data.copy()[emline['Ha-6564'],:,:],
                                        mask=self.maps_bm.flagged(
                                        dapmaps['EMLINE_GFLUX_MASK'].data[emline['Ha-6564'],:,:],
                                                              'DONOTUSE')).ravel()[unique_indx]
        except KeyError as e:
            print(e)
            warnings.warn('No halpha-data will be available.')
            return (-999.,)*8

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

        chi = dapmaps['EMLINE_FOM'].data[2,:,:].copy().ravel()[unique_indx]

        # Convert velocities to redshift
        z = (vel*(1+redshift) + astropy.constants.c.to('km/s').value*redshift) \
                                    / astropy.constants.c.to('km/s').value

        # Get the bins within an on-sky circular aperture of 2.5 arcsec
        cooext = 'SPX_SKYCOO' if deconstruct else 'BIN_LWSKYCOO'
        d2 = numpy.sum(numpy.square(dapmaps[cooext].data), axis=0).ravel()[unique_indx]
        center = d2 < 1.25*1.25

        # Get the flux-weighted mean velocity at the center
        wsum = numpy.ma.sum(flux[center])

        halpha_z = -999. if numpy.all(flux.mask[center] | z.mask[center]) \
                                or numpy.isclose(wsum, 0.) \
                            else numpy.ma.sum(flux[center]*z[center])/wsum
        if not numpy.isfinite(halpha_z):
            raise ValueError('H-alpha z is not finite: {0}'.format(halpha_z))

#        halpha_z = -999. if numpy.all(flux.mask[center]) or numpy.isclose(wsum, 0.) \
#                            else numpy.ma.sum(flux[center]*z[center])/wsum

        # Get the growth ranges of the H-alpha velocity
        if numpy.all(vel.mask):
            halpha_gvel_lo = -999.
            halpha_gvel_hi = -999.
            halpha_gvel_lo_clip = -999.
            halpha_gvel_hi_clip = -999.
        else:
            halpha_gvel_lo, halpha_gvel_hi = sample_growth(vel.compressed(), [0.025,0.975])
            halpha_gvel_lo_clip, halpha_gvel_hi_clip = sample_growth(sigma_clip(vel).compressed(),
                                                                     [0.025,0.975])

        # Get the flux-weighted, corrected sigma within 1 Re
        within_1re = r_re < 1.
        wsum = numpy.ma.sum(flux[within_1re])
        halpha_gsigma2_1re = -999. if numpy.all(flux.mask[within_1re] | sig2corr.mask[within_1re]) \
                                        or numpy.isclose(wsum, 0.) \
                                    else numpy.ma.sum(flux[within_1re]*sig2corr[within_1re])/wsum
        if not numpy.isfinite(halpha_gsigma2_1re):
            raise ValueError('H-alpha gsigma^2 is not finite: {0}'.format(halpha_gsigma2_1re))
        halpha_gsigma_1re = numpy.sqrt(halpha_gsigma2_1re) if halpha_gsigma2_1re > 0 else -999.

        # Get the high-growth of the H-alpha velocity dispersion
        if numpy.all(sig2corr.mask):
            halpha_gsigma2_hi = -999.
            halpha_gsigma2_hi_clip = -999.
        else:
            halpha_gsigma2_hi = sample_growth(sig2corr.compressed(), 0.975)
            halpha_gsigma2_hi_clip = sample_growth(sigma_clip(sig2corr).compressed(), 0.975)

        halpha_gsigma_hi = numpy.sqrt(halpha_gsigma2_hi) if halpha_gsigma2_hi > 0 else -999.
        halpha_gsigma_hi_clip = numpy.sqrt(halpha_gsigma2_hi_clip) \
                                        if halpha_gsigma2_hi_clip > 0 else -999.

        # Get the median reduced chi-square for the emission-line module
        emline_rchi2_1re = -999. if numpy.sum(within_1re) == 0 \
                                        else numpy.median(chi[within_1re])
        

        return halpha_z, halpha_gvel_lo, halpha_gvel_hi, halpha_gvel_lo_clip, halpha_gvel_hi_clip, \
                        halpha_gsigma_1re, halpha_gsigma_hi, halpha_gsigma_hi_clip, \
                        emline_rchi2_1re

    def _emission_line_metrics(self, dapmaps, deconstruct, moment0=False):

        # Unique bins
        # TODO: Select the BINID channel using the channel dictionary
        unique, indx = numpy.unique(dapmaps['BINID'].data[3,:,:].ravel(), return_index=True)
        if unique[0] == -1:
            unique = unique[1:]
            unique_indx = indx[1:]

        flux_ext = 'EMLINE_SFLUX' if moment0 else 'EMLINE_GFLUX'
        neml = dapmaps[flux_ext].shape[0]

        # Pull the data from the maps file
        cooext = 'SPX_ELLCOO' if deconstruct else 'BIN_LWELLCOO'
        r_re = dapmaps[cooext].data.copy()[1,:,:].ravel()[unique_indx]

        flux = numpy.ma.MaskedArray(dapmaps[flux_ext].data.copy(),
                                    mask=self.maps_bm.flagged(dapmaps[flux_ext+'_MASK'].data,
                                                    'DONOTUSE')).reshape(neml,-1)[:,unique_indx]

        ew_ext = 'EMLINE_SEW' if moment0 else 'EMLINE_GEW'
        ew = numpy.ma.MaskedArray(dapmaps[ew_ext].data.copy(),
                                  mask=self.maps_bm.flagged(dapmaps[ew_ext+'_MASK'].data,
                                                    'DONOTUSE')).reshape(neml,-1)[:,unique_indx]

        # Get the total flux at the center
        cooext = 'SPX_SKYCOO' if deconstruct else 'BIN_LWSKYCOO'
        d2 = numpy.sum(numpy.square(dapmaps[cooext].data), axis=0).ravel()[unique_indx]
        center = d2 < 1.25*1.25
        if numpy.sum(center) == 0:
            warnings.warn('No data near center!')
            emline_flux_cen = numpy.full(neml, -999., dtype=numpy.float)
        else:
            emline_flux_cen = numpy.ma.sum(flux[:,center], axis=1).filled(-999.)
            
        # Get the total flux, mean surface-brightness, and mean
        # equivalent width within 1 Re
        within_1re = r_re < 1.
        if numpy.sum(within_1re) == 0:
            warnings.warn('No data within 1 Re!')
            emline_flux_1re = numpy.full(neml, -999., dtype=numpy.float)
            emline_sb_1re = numpy.full(neml, -999., dtype=numpy.float)
            emline_ew_1re = numpy.full(neml, -999., dtype=numpy.float)
        else:
            emline_flux_1re = numpy.ma.sum(flux[:,within_1re], axis=1).filled(-999.)
            emline_sb_1re = numpy.ma.mean(flux[:,within_1re], axis=1).filled(-999.)
            emline_ew_1re = numpy.ma.mean(ew[:,within_1re], axis=1).filled(-999.)

        # Get the total flux within the full field-of-fiew
        emline_flux_tot = numpy.ma.sum(flux, axis=1).filled(-999.)

        # Get the peak surface-brightness and equivalent width within
        # the full field-of-fiew
        emline_sb_peak = numpy.ma.amax(flux, axis=1).filled(-999.)
        emline_ew_peak = numpy.ma.amax(ew, axis=1).filled(-999.)

        return emline_flux_cen, emline_flux_1re, emline_flux_tot, emline_sb_1re, \
                        emline_sb_peak, emline_ew_1re, emline_ew_peak

    def _spectral_index_metrics(self, dapmaps, deconstruct, specindex_units):

        # Unique bins
        # TODO: Select the BINID channel using the channel dictionary
        unique, indx = numpy.unique(dapmaps['BINID'].data[4,:,:].ravel(), return_index=True)
        if unique[0] == -1:
            unique = unique[1:]
            unique_indx = indx[1:]

        nindx = dapmaps['SPECINDEX'].shape[0]

        # Pull the data from the maps file
        cooext = 'SPX_ELLCOO' if deconstruct else 'BIN_LWELLCOO'
        r_re = dapmaps[cooext].data.copy()[1,:,:].ravel()[unique_indx]
        specindex = numpy.ma.MaskedArray(dapmaps['SPECINDEX'].data.copy(),
                                         mask=self.maps_bm.flagged(dapmaps['SPECINDEX_MASK'].data,
                                                    'DONOTUSE')).reshape(nindx,-1)[:,unique_indx]
        specindex_corr = numpy.ma.MaskedArray(
                            dapmaps['SPECINDEX_CORR'].data.copy().reshape(nindx,-1)[:,unique_indx],
                                              mask=specindex.mask.copy())

        # Get the corrected indices; Unitless and ang units are treated
        # the same.
        mag = specindex_units == 'mag'
        ang = numpy.invert(mag)
        specindex_corr[ang] = specindex[ang]*specindex_corr[ang]
        specindex_corr[mag] = specindex[mag]+specindex_corr[mag]

        # Get the growth limits
        specindex_lo = numpy.empty(nindx, dtype=float)
        specindex_hi = numpy.empty(nindx, dtype=float)
        specindex_lo_clip = numpy.empty(nindx, dtype=float)
        specindex_hi_clip = numpy.empty(nindx, dtype=float)
        for i in range(nindx): 
            if numpy.all(specindex_corr.mask[i,:]):
                specindex_lo[i] = -999.
                specindex_hi[i] = -999.
                specindex_lo_clip[i] = -999.
                specindex_hi_clip[i] = -999.
                continue
            specindex_lo[i], specindex_hi[i] = sample_growth(specindex_corr[i,:].compressed(),
                                                             [0.025,0.975])
            specindex_lo_clip[i], specindex_hi_clip[i] \
                        = sample_growth(sigma_clip(specindex_corr[i,:]).compressed(), [0.025,0.975])

        # Get the median within 1 Re
        specindex_1re = numpy.ma.median(specindex[:,r_re < 1.], axis=1).filled(-999.)

        return specindex_lo, specindex_hi, specindex_lo_clip, specindex_hi_clip, specindex_1re

    def file_name(self):
        """Return the name of the DRP complete database."""
        name = defaults.dapall_file(drpver=self.drpver, dapver=self.dapver,
                                    analysis_path=self.analysis_path)
        return os.path.split(name)[1]

    def file_path(self):
        """Return the full pat to the DRP complete database."""
        return defaults.dapall_file(drpver=self.drpver, dapver=self.dapver,
                                    analysis_path=self.analysis_path)

    def method_database(self, method, deconstruct, drpall):
        """
        Construct the DAPall database for the provided DAP method.

        Args:
            method (:obj:`str`):
                Name of the DAP method
            deconstruct (:obj:`bool`):
                Flag that the method deconstructs the bins during the
                emission-line fitting module.
            drpall (`numpy.recarray`_):
                Primary DRPall data table

        Returns:
            `numpy.recarray`_: DAPall database.
        """
        # Construct the list of files to look for
        dapfiles = self._construct_maps_file_list(method)

        # Initialize the output database
        db = init_record_array(self.ndap, self._table_dtype(self.nmom, self.neml, self.nindx))

        # Add the basic information
        db['PLATE'] = self.plate
        db['IFUDESIGN'] = self.ifudesign
        db['PLATEIFU'] = numpy.array(['{0}-{1}'.format(p,i) 
                                            for p,i in zip(self.plate, self.ifudesign)])
        # For now only analyzing cubes!
        db['MODE'][:] = 'CUBE'
        db['DAPTYPE'][:] = method

        # Copy from DRPall
        columns_from_drpall = ['MANGAID', 'OBJRA', 'OBJDEC', 'IFURA', 'IFUDEC', 'MNGTARG1',
                               'MNGTARG2', 'MNGTARG3', 'NSA_Z', 'NSA_ZDIST',
                               'NSA_ELPETRO_BA', 'NSA_ELPETRO_PHI', 'NSA_ELPETRO_TH50_R',
                               'NSA_SERSIC_BA', 'NSA_SERSIC_PHI', 'NSA_SERSIC_TH50',
                               'NSA_SERSIC_N', 'VERSDRP2', 'VERSDRP3', 'VERSCORE', 'VERSUTIL',
                               'DRP3QUAL']

        columns_from_maps_hdr = ['VERSDAP', 'DAPQUAL', 'RDXQAKEY', 'BINKEY', 'SCKEY', 'ELMKEY',
                                 'ELFKEY', 'SIKEY', 'BINTYPE']

        # Itereate through each maps file
        for i in range(self.ndap):
            print('Processing {0}; observation {1}/{2}'.format(method, i+1,self.ndap))#, end='\r')

            # Find the index in the drpall file
            indx = numpy.where(drpall['PLATEIFU'] == db['PLATEIFU'][i])[0]
            if len(indx) > 1:
                warnings.warn('More than one entry in DRPall file: {0}.'.format(db['PLATEIFU'][i]))
            if len(indx) != 1:
                indx = [-1]
                warnings.warn('Could not find {0} in DRPall file!'.format(db['PLATEIFU'][i]))
            db['DRPALLINDX'][i] = indx[0]

            # Add data from the DRPall file
            if db['DRPALLINDX'][i] != -1:
                for c in columns_from_drpall:
                    db[c][i] = drpall[c][db['DRPALLINDX'][i]]

            # Set the cosmological distances based on the NSA redshift
            db['LDIST_NSA_Z'][i] = self.cosmo.luminosity_distance(db['NSA_Z'][i]).value \
                                        if db['NSA_Z'][i] > 0 else -999.
            db['ADIST_NSA_Z'][i] = self.cosmo.angular_diameter_distance(db['NSA_Z'][i]).value \
                                        if db['NSA_Z'][i] > 0 else -999.

            # Open the maps file
            if not os.path.isfile(dapfiles[i]):
                db['DAPDONE'][i] = 0    # A fault occurred such that the MAPS file was not produced
                continue
            db['DAPDONE'][i] = 1        # DAP successfully produced the file
            dapmaps = fits.open(dapfiles[i])

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
                db['DATEDAP'][i] = time.strftime('%Y-%m-%d',
                                                 time.localtime(os.path.getmtime(dapfiles[i])))

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

            # Determine the coverage of the datacube (independent of the
            # DAP method)
            r = numpy.ma.MaskedArray(dapmaps['SPX_ELLCOO'].data.copy()[0,:,:], mask=basic_map_mask)
            db['RCOV90'][i] = self._radial_coverage_metric(r.compressed(),
                                                           1-db['NSA_ELPETRO_BA'][i])

            # Get the SRD-based S/N metric (independent of the DAP
            # method)
            r_re = numpy.ma.MaskedArray(dapmaps['SPX_ELLCOO'].data.copy()[1,:,:],
                                        mask=numpy.invert(dapmaps['SPX_SNR'].data > 0))
            db['SNR_MED'][i], db['SNR_RING'][i] = self._srd_snr_metric(dapmaps)

            # Get the mean surface brightness within 1 Re used by the
            # stellar-continuum fitting (independent of the DAP method)
            mflux = numpy.ma.MaskedArray(dapmaps['SPX_MFLUX'].data.copy(),
                                         mask=numpy.invert(dapmaps['SPX_SNR'].data > 0))
            within_1re = r_re < 1
            if numpy.sum(within_1re) == 0:
                warnings.warn('No data within 1 Re!')
                db['SB_1RE'][i] = -999.
            else:
                db['SB_1RE'][i] = numpy.ma.mean(mflux[within_1re])

            # Get the spatial binning metrics
            db['BIN_RMAX'][i], db['BIN_R_N'][i], db['BIN_R_SNR'][i] \
                            = self._binning_metrics(dapmaps)

            # Get the metrics for the stellar kinematics
            db['STELLAR_Z'][i], db['STELLAR_VEL_LO'][i], db['STELLAR_VEL_HI'][i], \
                    db['STELLAR_VEL_LO_CLIP'][i], db['STELLAR_VEL_HI_CLIP'][i], \
                    db['STELLAR_SIGMA_1RE'][i], db['STELLAR_RCHI2_1RE'][i] \
                            = self._stellar_kinematics_metrics(dapmaps, db['Z'][i])

            # Get the metrics for the emission-line kinematics
            db['HA_Z'][i], db['HA_GVEL_LO'][i], db['HA_GVEL_HI'][i], db['HA_GVEL_LO_CLIP'][i], \
                    db['HA_GVEL_HI_CLIP'][i], db['HA_GSIGMA_1RE'][i], db['HA_GSIGMA_HI'][i], \
                    db['HA_GSIGMA_HI_CLIP'][i], db['EMLINE_RCHI2_1RE'][i] \
                            = self._halpha_kinematics_metrics(dapmaps, deconstruct, db['Z'][i])

            # Get the moment-based emission-line metrics
            db['EMLINE_SFLUX_CEN'][i], db['EMLINE_SFLUX_1RE'][i], db['EMLINE_SFLUX_TOT'][i], \
                    db['EMLINE_SSB_1RE'][i], db['EMLINE_SSB_PEAK'][i], db['EMLINE_SEW_1RE'][i], \
                    db['EMLINE_SEW_PEAK'][i] \
                            = self._emission_line_metrics(dapmaps, deconstruct, moment0=True)

            # Get the Gaussian-based emission-line metrics
            db['EMLINE_GFLUX_CEN'][i], db['EMLINE_GFLUX_1RE'][i], db['EMLINE_GFLUX_TOT'][i], \
                    db['EMLINE_GSB_1RE'][i], db['EMLINE_GSB_PEAK'][i], db['EMLINE_GEW_1RE'][i], \
                    db['EMLINE_GEW_PEAK'][i] = self._emission_line_metrics(dapmaps, deconstruct)

            # Get the spectral-index metrics
            db['SPECINDEX_LO'][i], db['SPECINDEX_HI'][i], db['SPECINDEX_LO_CLIP'][i], \
                    db['SPECINDEX_HI_CLIP'][i], db['SPECINDEX_1RE'][i] \
                            = self._spectral_index_metrics(dapmaps, deconstruct, self.spindx_units)

            # Estimate the star-formation rate
            log_Mpc_in_cm = numpy.log10(astropy.constants.pc.to('cm').value) + 6
            log_halpha_luminosity_1re = numpy.log10(4*numpy.pi) \
                        + numpy.ma.log10(db['EMLINE_GFLUX_1RE'][i,self.elfit_channels['Ha-6564']]) \
                        - 17 + 2*numpy.ma.log10(db['LDIST_Z'][i]) + 2*log_Mpc_in_cm
            db['SFR_1RE'][i] = numpy.ma.power(10, log_halpha_luminosity_1re - 41.27).filled(-999.)
            log_halpha_luminosity_tot = numpy.log10(4*numpy.pi) \
                        + numpy.ma.log10(db['EMLINE_GFLUX_TOT'][i,self.elfit_channels['Ha-6564']]) \
                        - 17 + 2*numpy.ma.log10(db['LDIST_Z'][i]) + 2*log_Mpc_in_cm
            db['SFR_TOT'][i] = numpy.ma.power(10, log_halpha_luminosity_tot - 41.27).filled(-999.)

        print('Processing {0}; observation {1}/{1}'.format(method, self.ndap))

        # Check that all the data is finite.
        found = False
        for name in db.dtype.names:
            if not (numpy.issubdtype(db[name].dtype, numpy.integer) \
                        or numpy.issubdtype(db[name].dtype, numpy.floating)):
                continue
            if not numpy.all(numpy.isfinite(db[name])):
                print('All values not finite in column {0} for method {1}'.format(name, method))
                found = True
        if found:
            raise ValueError('All values not finite in DAPall table.')

        return db

    def update(self, plan, methods=None):
        """
        Update the DAPall file.

        Currently, this constructs the entire DAPall file from
        scratch, regardless of whether or not some fraction of the
        output files have already been analyzed.
        
        Args:
            plan (:class:`mangadap.par.analysisplan.AnalysisPlanSet`):
                The plan object used by the DAP.
            methods (:obj:`str`, :obj:`list`, optional):
                Specify a set of methods in the DAP plan file to
                include in the DAPall file. If not provided, all
                methods in the plan file are included.

        Raises:
            ValueError:
                Raised if a provided method is not in the provided
                set of analysis plans, or if the methods do not all
                have the same number of emission-line bandpasses,
                emission-lines to fit, and number of spectral
                indices.
        """
        if self.readonly:
            raise ValueError('Object is read-only.')

        # Check the input type
        # TODO: Plan isn't really needed.  Could just troll the
        # analysis_path for the directory names...
        if not isinstance(plan, AnalysisPlanSet):
            raise TypeError('Input plan must have type AnalysisPlanSet.')

        # Get the full list of available plan methods
        plan_methods = []
        deconstruct = []
        for p in plan:
            bin_method = SpatiallyBinnedSpectra.define_method(p['bin_key'])
            sc_method = StellarContinuumModel.define_method(p['continuum_key'])
            el_method = EmissionLineModel.define_method(p['elfit_key'])
            plan_methods += [defaults.dap_method(bin_method['key'],
                                                 sc_method['fitpar']['template_library_key'],
                                                 el_method['continuum_tpl_key'])]
            deconstruct += [el_method['deconstruct_bins'] != 'ignore']

        # Check that the plan methods are unique.  This should never be
        # raised, unless it also causes problems in the DAP itself.
        if len(numpy.unique(plan_methods)) != len(plan_methods):
            raise ValueError('All of the plan methods must be unique!')

        # Get the list of methods to look for
        if methods is None:
            self.methods = numpy.array(plan_methods)
            self.deconstruct = numpy.array(deconstruct)
            use_plans = numpy.ones(len(plan_methods), dtype=bool)
        else:
            # Make sure the provided list is unique
            if len(numpy.unique(methods)) != len(methods):
                raise ValueError('Must provide unique list of methods.')
            # All the methods to use must be in the provided plan
            if not numpy.all( [m in plan_methods for m in methods ]):
                raise ValueError('All provided methods must be in provided plan file.')

            # Find the index of the plans to use
            # TODO: Re-write this as a list compression
            self.methods = []
            self.deconstruct = []
            use_plans = numpy.zeros(len(plan_methods), dtype=bool)
            for i in range(len(plan_methods)):
                if plan_method[i] in methods:
                    use_plans[i] = True
                    self.methods += [plan_method[i]]
                    self.deconstruct += [deconstruct[i]]
            self.methods = numpy.array(self.methods)
            self.deconstruct = numpy.array(self.deconstruct)

        # Check that all the plans to use have the same emission-line
        # bandpass channels, ...
        elmom_keys, elmom_channels \
                    = numpy.array([DAPall._emission_line_moment_db_info(p['elmom_key'])
                                        for p in plan[use_plans]], dtype=object).T
        if numpy.any([ len(k) != len(elmom_channels[0]) or len(set(k)-set(elmom_channels[0])) > 0
                        for k in elmom_channels]):
            raise ValueError('All plans must provide the same emission-line bandpass channels.')
        self.elmom_channels = elmom_channels[0].copy()
        del elmom_channels
        self.nmom = len(self.elmom_channels)

        # ..., emission-line fit channels, and...
        elfit_keys, elfit_channels \
                    = numpy.array([DAPall._emission_line_db_info(p['elfit_key'])
                                        for p in plan[use_plans]], dtype=object).T
        if numpy.any([ len(k) != len(elfit_channels[0]) or len(set(k)-set(elfit_channels[0])) > 0
                        for k in elfit_channels]):
            raise ValueError('All plans must provide the same emission-line fit channels.')
        self.elfit_channels = elfit_channels[0].copy()
        del elfit_channels
        self.neml = len(self.elfit_channels)

        # ..., and spectral-index channels
        abs_keys, bhd_keys, spindx_channels, spindx_units \
                    = numpy.array([DAPall._spectral_index_db_info(p['spindex_key'])
                                        for p in plan[use_plans]], dtype=object).T
        if numpy.any([ len(k) != len(spindx_channels[0]) or len(set(k)-set(spindx_channels[0])) > 0
                        for k in spindx_channels]):
            raise ValueError('All plans must provide the same emission-line fit channels.')
        self.spindx_channels = spindx_channels[0].copy()
        self.spindx_units = spindx_units[0].copy()
        del spindx_channels, spindx_units
        self.nindx = len(self.spindx_channels)

        # PlateIFUs that should be analyzed per DAP method
        # TODO: Somehow sort these?
        self.plate, self.ifudesign = self._get_completed_observations()
        # Number of completed observations per analysis method
        self.ndap = len(self.plate)

        # Report
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, '-'*50)
            log_output(self.loggers, 1, logging.INFO, '{0:^50}'.format('RUNNING DAPALL UPDATE'))
            log_output(self.loggers, 1, logging.INFO, '-'*50)
            log_output(self.loggers, 1, logging.INFO,
                       'Completed observations: {0}'.format(self.ndap))
            log_output(self.loggers, 1, logging.INFO,
                       'Available methods: {0}'.format(plan_methods))
            log_output(self.loggers, 1, logging.INFO,
                       'Methods for output: {0}'.format(self.methods))
            log_output(self.loggers, 1, logging.INFO,
                       'Emission line moments: {0}'.format(self.nmom))
            log_output(self.loggers, 1, logging.INFO,
                       'Emission lines fit: {0}'.format(self.neml))
            log_output(self.loggers, 1, logging.INFO,
                       'Spectral indices: {0}'.format(self.nindx))

        # If the HDUList is already initialized, close and reinitialize
        # it
        if self.hdu is not None:
            self.hdu.close()
            self.hdu = None

        # Open the DRPall file
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO,
                       'Reading DRPall: {0}'.format(self.drpall_file))
        hdu_drpall = fits.open(self.drpall_file)
        drpall = hdu_drpall[1].data

        # Create the primary header
        hdr = fits.Header()
        hdr['DATE'] = (time.strftime('%Y-%m-%d',time.gmtime()), 'UTC date created')
        hdr['VERSDRP3'] = (self.drpver, 'DRP version')
        hdr['VERSDAP'] = (self.dapver, 'DAP version')

        # Add the emission-line moment channel header
        hdr = self._add_channels_to_header(hdr, self.elmom_channels, 'ELS',
                                           'Summed emission-line element')
        hdr = self._add_channels_to_header(hdr, self.elfit_channels, 'ELG',
                                           'Gaussian fit emission-line element')
        hdr = self._add_channels_to_header(hdr, self.spindx_channels, 'SPI',
                                           'Spectral-indx element', units=self.spindx_units)

        # Construct a DAPall database for each analysis method and add
        # it as a binary table to the list of HDUs.
        self.hdu = [fits.PrimaryHDU(header=hdr)]
        for method, deconstruct in zip(self.methods, self.deconstruct):
            method_db = self.method_database(method, deconstruct, drpall)
            self.hdu += [fits.BinTableHDU.from_columns(
                                [fits.Column(name=n, format=rec_to_fits_type(method_db[n]),
                                             array=method_db[n]) for n in method_db.dtype.names],
                                                       name=method)]
        self.hdu = fits.HDUList(self.hdu)

        # Write the file (use of this function is probably overkill)
        DAPFitsUtil.write(self.hdu, self.file_path(), clobber=True, checksum=True,
                          loggers=self.loggers)
        

