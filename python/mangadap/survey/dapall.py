# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
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

*Class usage examples*:
    Add some usage comments here!

*Revision history*:
    | **19 Aug 2016**: Original Implementation by K. Westfall (KBW)

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

# For versioning
from scipy import interpolate

from astropy.io import fits
import astropy.constants
from astropy.cosmology import FlatLambdaCDM
import astropy.units
from astropy.stats import sigma_clip

from ..drpfits import DRPFits, DRPQuality3DBitMask
from ..dapqual import DAPQualityBitMask
from ..dapmaps import DAPMapsBitMask
from ..par.analysisplan import AnalysisPlanSet
from .drpcomplete import DRPComplete
from ..config import defaults
from ..util.fileio import init_record_array, rec_to_fits_type, write_hdu
from ..proc.util import _select_proc_method
from ..proc.absorptionindexdb import AbsorptionIndexDB
from ..proc.bandheadindexdb import BandheadIndexDB
from ..proc.emissionlinedb import EmissionLineDB
from ..proc.spectralindices import SpectralIndicesDef, available_spectral_index_databases
from ..proc.emissionlinemodel import EmissionLineModelDef, available_emission_line_modeling_methods

from matplotlib import pyplot

__author__ = 'Kyle Westfall'

def column_dictionary(hdu, ext):
    columndict = {}
    for k, v in hdu[ext].header.items():
        if k[0] == 'C':
            try:
                i = int(k[1:])-1
            except ValueError:
                continue
            columndict[v] = i
    return columndict


def growth(a, growth_fracs, default=-9999.):
    _a = a.compressed() if isinstance(a, numpy.ma.MaskedArray) else numpy.atleast_1d(a).ravel()
#    print(len(_a))
    if len(_a) < 2:
        return tuple([default]*len(growth_fracs))
    srt = numpy.argsort(_a)
    grw = (numpy.arange(_a.size,dtype=float)+1)/_a.size
#    print(grw[0:2], '...', grw[-2:])
    interpolator = interpolate.interp1d(grw, _a[srt], fill_value='extrapolate')
    return tuple(interpolator(growth_fracs))


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
        drpver (str): (**Optional**) DRP version.  Default is to use
            $MANGADRP_VER.
        redux_path (str): (**Optional**) Top-level directory with the
            DRP products; defaults to $MANGA_SPECTRO_REDUX/$MANGADRP_VER
        dapver (str): (**Optional**) DAP version.  Default is to use
            $MANGADAP_VER.
        dapsrc (str): (**Optional**) Path to the DAP source code.
            Default is to use $MANGADAP_DIR.
        analysis_path (str): (**Optional**) Top-level output directory
            for the DAP results; defaults to
            $MANGA_SPECTRO_ANALYSIS/$MANGADRP_VER/$MANGADAP_VER
        update (bool): (**Optional**) Update any existing file with any
            new output that has been found

    Raises:
        TypeError: Raised if the input plan is not a
            :class:`mangadap.par.analysisplan.AnalysisPlanSet` object.
        FileNotFoundError: Raise if the DRPall or the DRPComplete file
            are not available.
    """
    def __init__(self, plan, methods=None, drpver=None, redux_path=None, dapver=None, dapsrc=None,
                 analysis_path=None, dap_common_path=None, clobber=False, update=True,
                 readonly=False):

        # Check the input type
        if not isinstance(plan, AnalysisPlanSet):
            raise TypeError('Input plan must have type AnalysisPlanSet.')

        # Set version
        self.version = '0.1'

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
        self.drpall = os.path.join(self.redux_path, 'drpall-{0}.fits'.format(self.drpver))
        print(self.drpall)

        # Get the full list of available plan methods
        self.plan_methods = numpy.array([ defaults.default_dap_method(plan=p) for p in plan ])

        # TODO: Check that the plan methods are unique?  (Will cause
        # other problems in the DAP itself if they're not!)

        # Initialize the output objects
        # Use the plans to set the number of emission lines for each
        # method
        self.neml = numpy.array([ self._number_of_emission_lines(p['elfit_key'],
                                                                 dapsrc=self.dapsrc)
                                    for p in plan ])

        # Use the plans to set the number of spectral indices
        self.nindx = numpy.array([ self._number_of_spectral_indices(p['spindex_key'],
                                                                    dapsrc=self.dapsrc)
                                    for p in plan ])

        self.str_len = self._init_string_lengths()
        self.hdu = None
        self.ndap = None
        self.methods = None
        self.plate = None
        self.ifudesign = None
        self.h = 1.0
        self.cosmo = None
        self.maps_bm = None

        if os.path.exists(self.file_path()):
            print('READING')
            self._read(mode='readonly' if readonly else 'update')

            # Check that the number of emission lines and indices
            # matches the fits file

        if self.hdu is None or update:
            self.update(methods=methods, clobber=clobber)


    def __getitem__(self, key):
        return self.hdu['DAPALL'].data[key]


    @staticmethod
    def _number_of_spectral_indices(key, dapsrc=None):
        db = _select_proc_method(key, SpectralIndicesDef,
                                 available_func=available_spectral_index_databases, dapsrc=dapsrc)
        nabs = 0 if db['absindex'] is None \
                    else AbsorptionIndexDB(db['absindex'], dapsrc=dapsrc).nsets
        nbhd = 0 if db['bandhead'] is None \
                    else BandheadIndexDB(db['bandhead'], dapsrc=dapsrc).nsets
        return nabs+nbhd


    @staticmethod
    def _number_of_emission_lines(key, dapsrc=None):
        db = _select_proc_method(key, EmissionLineModelDef,
                                 available_func=available_emission_line_modeling_methods,
                                 dapsrc=dapsrc)
        return EmissionLineDB(db['emission_lines'], dapsrc=dapsrc).nsets


    def _read(self, mode='readonly'):
        """Read the data in the existing file at :func:`file_path`."""
        if self.hdu is not None:
            self.hdu.close()
            self.hdu = None
        self.hdu = fits.open(self.file_path(), mode=mode)
        self.ndap = self.hdu['DAPALL'].header['NAXIS2']
        print('Read data: {0} rows'.format(self.ndap))


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
                 'EML_NAME':20,
                 'SPINDX_NAME':20
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
                 ('NSA_Z', numpy.float),
                 ('NSA_ZDIST', numpy.float),
                 ('NSA_ELPETRO_BA', numpy.float),
                 ('NSA_ELPETRO_PHI', numpy.float),
                 ('NSA_ELPETRO_TH50_R', numpy.float),
                 ('NSA_SERSIC_BA', numpy.float),
                 ('NSA_SERSIC_PHI', numpy.float),
                 ('NSA_SERSIC_TH50', numpy.float),
                 ('NSA_SERSIC_N', numpy.float),
                 ('LDIST_NSA_Z', numpy.float),
                 ('ADIST_NSA_Z', numpy.float),
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
                 ('RMAX_RE', numpy.float),
#                 ('RCOV90', numpy.float),
#                 ('SNR', numpy.float),
                 ('SNR_HRE', numpy.float),
                 ('SNR_1RE', numpy.float),
                 ('SNR_2RE', numpy.float),
                 ('STAR_FLUX', numpy.float),
                 ('STAR_Z', numpy.float),
                 ('STAR_VEL_LO', numpy.float),
                 ('STAR_VEL_HI', numpy.float),
                 ('STAR_VEL_LO_CLIP', numpy.float),
                 ('STAR_VEL_HI_CLIP', numpy.float),
                 ('STAR_DISP_1RE', numpy.float),
                 ('SC_RCHI2_1RE', numpy.float),
                 ('HA_FLUX', numpy.float),
                 ('HA_Z', numpy.float),
                 ('HA_VEL_LO', numpy.float),
                 ('HA_VEL_HI', numpy.float),
                 ('HA_VEL_LO_CLIP', numpy.float),
                 ('HA_VEL_HI_CLIP', numpy.float),
                 ('HA_DISP_1RE', numpy.float),
                 ('HA_DISP_HI', numpy.float),
                 ('HA_DISP_HI_CLIP', numpy.float),
                 ('EML_NAME', '<U{0:d}'.format(self.str_len['EML_NAME']*neml)),
                 ('EML_FLUX_TOTAL', numpy.float, neml),
                 ('EML_SB_1RE', numpy.float, neml),
                 ('EML_SB_PEAK', numpy.float, neml),
#                 ('EML_EW_PEAK', numpy.float, neml),
                 ('SPINDX_NAME', '<U{0:d}'.format(self.str_len['SPINDX_NAME']*nindx)),
                 ('SPINDX_LO', numpy.float, nindx),
                 ('SPINDX_MD', numpy.float, nindx),
                 ('SPINDX_HI', numpy.float, nindx),
                 ('SPINDX_LO_CLIP', numpy.float, nindx),
                 ('SPINDX_MD_CLIP', numpy.float, nindx),
                 ('SPINDX_HI_CLIP', numpy.float, nindx)#,
#                 ('SFR', numpy.float)
               ]


    def _get_completed_observations(self):
        """
        Get the list of DRP completed (and prospectively DAP analyzed)
        observations.

        Raises:
            FileNotFoundError: Raised if the DRPComplete file is not
                available.
        .. todo::
            - Create the DRPComplete file if it doesn't exist?
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
        plates = numpy.array([ numpy.array(plt) ]*nmethods).ravel().astype(str)
        ifudesigns = numpy.array([ numpy.array(ifu) ]*nmethods).ravel().astype(str)

        return methods, plates, ifudesigns


    def _construct_maps_file_list(self, plt, ifu):
        """
        Construct the list of MAPS files that should be in/added to the
        DAPall file.
        """
        # Construct the file names
        methods, plates, ifudesigns = self._combine_plateifu_methods(plt, ifu)
        return numpy.array( ['{0}/{1}/{2}/{3}/manga-{2}-{3}-MAPS-{1}.fits.gz'.format(
                                self.analysis_path, p, i, m) 
                                    for p,i,m in zip(methods, plates, ifudesigns) ])


    def file_name(self):
        """Return the name of the DRP complete database."""
        return ('dapall-{0}-{1}.fits'.format(self.drpver, self.dapver))


    def file_path(self):
        """Return the full pat to the DRP complete database."""
        return os.path.join(self.analysis_path, self.file_name())


    def update(self, methods=None, clobber=False):
        """
        Update the DAPall file
        
        If clobber is True, the entire DAPall file is reconstructed from
        scratch.  Otherwise, any additional data is appended to the end
        of the file.

        Args:
            methods (str,list): (**Optional**) Specify a set of methods
                in the DAP plan file to include in the DAPall file.  If
                not provided, all methods in the plan file are included.
            clobber (bool): (**Optional**) Clobber any existing data and
                rewrite the file.  Default is False.

        Raises:
            ValueError: Raised if a provided method is not in the provided
                set of analysis plans.
        """

        # Get the list of methods to look for
        if methods is None:
            self.methods = self.plan_methods
        else:
            self.methods = numpy.atleast_1d(methods)
            if ~numpy.all( [m in self.plan_methods for m in self.methods ]):
                raise ValueError('All provided methods must be in provided plan file.')

        # Construct the list of files to look for
        self.plate, self.ifudesign = self._get_completed_observations()
        methodlist, platelist, ifulist = self._combine_plateifu_methods(self.plate, self.ifudesign)
        dapfiles = self._construct_maps_file_list(self.plate, self.ifudesign)

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
        if not os.path.isfile(self.drpall):
            raise FileNotFoundError('No file: {0}'.format(self.drpall))
        hdu_drpall = fits.open(self.drpall)
        drpall = hdu_drpall[1].data

        # Initialize the cosmology object
        self.H0 = 100 * self.h * astropy.units.km / astropy.units.s / astropy.units.Mpc
        self.cosmo = FlatLambdaCDM(H0=self.H0, Om0=0.3)

        # Initialize the bitmask
        self.maps_bm = DAPMapsBitMask()

        # Initialize the output database
        self.ndap = len(dapfiles)
        db = init_record_array(self.ndap, self._table_dtype(max(self.neml), max(self.nindx)))

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
            indx = numpy.where( drpall['PLATEIFU'] == db['PLATEIFU'][i])[0]
            if len(indx) != 1:
                raise ValueError('Could not find {0} in DRPall file.'.format(db['PLATEIFU'][i]))
            db['DRPALLINDX'][i] = indx[0]

            # Add data from the DRPall file
            for c in columns_from_drpall:
                db[c][i] = drpall[c][db['DRPALLINDX'][i]]

            # Open the maps file
            if not os.path.isfile(f):
                db['DAPDONE'][i] = 0    # DAP successfully produced the file
                print('No file: {0}'.format(f))
                continue
            db['DAPDONE'][i] = 1        # A fault occurred such that the MAPS file was not produced
            dapmaps = fits.open(f)

            # Get the dictionary with the channel names
            emline = column_dictionary(dapmaps, 'EMLINE_GFLUX')
            neml = len(emline)
            sindx = column_dictionary(dapmaps, 'SPECINDEX')
            nindx = len(sindx)
            # And set the channel names to the table
            srt = numpy.argsort(list(emline.values()))
            db['EML_NAME'][i] = ''.join([ '{0}'.format(n).rjust(self.str_len['EML_NAME']) \
                                                for n in numpy.array(list(emline.keys()))[srt] ])
            srt = numpy.argsort(list(sindx.values()))
            db['SPINDX_NAME'][i] = ''.join([ '{0}'.format(n).rjust(self.str_len['SPINDX_NAME']) \
                                                for n in numpy.array(list(sindx.keys()))[srt] ])

            # Add any checks that the DRPall data matches the header
            # data in the MAPS files?

            # Add information from the DAP MAPS file header
            for c in columns_from_maps_hdr:
                db[c][i] = dapmaps['PRIMARY'].header[c]

            # Additional elements from the header the need to be handled
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

            # Use the NSA redshift and astropy to compute luminosity and
            # angular-diameter distances
            db['LDIST_NSA_Z'][i] = self.cosmo.luminosity_distance(db['NSA_Z'][i]).value
            db['ADIST_NSA_Z'][i] = self.cosmo.angular_diameter_distance(db['NSA_Z'][i]).value

            r = dapmaps['SPX_ELLCOO'].data.copy()[0,:,:]
            r_re = r/db['NSA_ELPETRO_TH50_R'][i]
            db['RMAX_RE'][i] = numpy.amax(r_re)

            snr = numpy.ma.MaskedArray(dapmaps['BIN_SNR'].data.copy(),
                                       mask=~(dapmaps['BIN_SNR'].data>0))
            bin_flux = numpy.ma.MaskedArray(dapmaps['BIN_MFLUX'].data.copy(),
                                       mask=self.maps_bm.flagged(dapmaps['BIN_MFLUX_MASK'].data,
                                                                 'DONOTUSE'))
            bin_fvar = numpy.ma.power(numpy.ma.MaskedArray(dapmaps['BIN_MFLUX_IVAR'].data.copy(),
                                       mask=self.maps_bm.flagged(dapmaps['BIN_MFLUX_MASK'].data,
                                                                 'DONOTUSE')), -1.0)
            svel = numpy.ma.MaskedArray(dapmaps['STELLAR_VEL'].data.copy(),
                                        mask=self.maps_bm.flagged(dapmaps['STELLAR_VEL_MASK'].data,
                                                             'DONOTUSE'))
            ssig = numpy.ma.MaskedArray(dapmaps['STELLAR_SIGMA'].data.copy(),
                                    mask=self.maps_bm.flagged(dapmaps['STELLAR_SIGMA_MASK'].data,
                                                              'DONOTUSE'))
            ssig2corr = numpy.square(ssig) - numpy.square(
                                numpy.ma.MaskedArray(dapmaps['STELLAR_SIGMACORR'].data.copy(),
                                    mask=self.maps_bm.flagged(dapmaps['STELLAR_SIGMA_MASK'].data,
                                                              'DONOTUSE')))
            scchi = numpy.ma.MaskedArray(dapmaps['STELLAR_CONT_RCHI2'].data.copy(),
                                         mask=~(dapmaps['BINID'].data > -1))

            havel = numpy.ma.MaskedArray(dapmaps['EMLINE_GVEL'].data.copy()[emline['Ha-6564'],:,:],
                                         mask=self.maps_bm.flagged(
                                            dapmaps['EMLINE_GVEL_MASK'].data[emline['Ha-6564'],:,:],
                                                'DONOTUSE'))
            hasig = numpy.ma.MaskedArray(dapmaps['EMLINE_GVEL'].data.copy()[emline['Ha-6564'],:,:],
                                         mask=self.maps_bm.flagged(
                                            dapmaps['EMLINE_GVEL_MASK'].data[emline['Ha-6564'],:,:],
                                                'DONOTUSE'))
            hasig2corr = numpy.square(hasig) - numpy.square(
                numpy.ma.MaskedArray(dapmaps['EMLINE_INSTSIGMA'].data.copy()[emline['Ha-6564'],:,:],
                                     mask=self.maps_bm.flagged(
                                        dapmaps['EMLINE_GSIGMA_MASK'].data[emline['Ha-6564'],:,:],
                                                'DONOTUSE')))

            line_flux = numpy.ma.MaskedArray(dapmaps['EMLINE_GFLUX'].data.copy(),
                                             mask=self.maps_bm.flagged(
                                                    dapmaps['EMLINE_GFLUX_MASK'].data, 'DONOTUSE'))

            spec_index = numpy.ma.MaskedArray(dapmaps['SPECINDEX'].data.copy(),
                                              mask=self.maps_bm.flagged(
                                                        dapmaps['SPECINDEX_MASK'].data, 'DONOTUSE'))

            skyx = dapmaps['SPX_SKYCOO'].data.copy()[0,:,:]
            skyy = dapmaps['SPX_SKYCOO'].data.copy()[1,:,:]

            # TODO: Need to incorporate covariance
#            indx = (r_re > 1.0) & (r_re < 1.5)
#            db['SNR'][i] = -9999.0 if numpy.sum(indx) == 0 else numpy.ma.median(snr[indx])

            indx = (r_re > 0.4) & (r_re < 0.6)
            db['SNR_HRE'][i] = -9999.0 if numpy.sum(indx) == 0 else numpy.ma.median(snr[indx])
            indx = (r_re > 0.9) & (r_re < 1.1)
            db['SNR_1RE'][i] = -9999.0 if numpy.sum(indx) == 0 else numpy.ma.median(snr[indx])
            indx = (r_re > 1.9) & (r_re < 2.1)
            db['SNR_2RE'][i] = -9999.0 if numpy.sum(indx) == 0 else numpy.ma.median(snr[indx])


            center = skyx*skyx+skyy*skyy < 1.25*1.25
            _svel = svel.copy()
            _svel[center] = numpy.ma.masked
            _havel = havel.copy()
            _havel[center] = numpy.ma.masked

            if numpy.sum(~numpy.ma.getmaskarray(_svel)) > 0:
                svel_clipped = sigma_clip(_svel)
                star_mask = numpy.ma.getmaskarray(svel_clipped)
                star_flux = numpy.ma.MaskedArray(dapmaps['BIN_MFLUX'].data.copy(), mask=star_mask)
#                svel_var = numpy.ma.power(numpy.ma.MaskedArray(
#                                                        dapmaps['STELLAR_VEL_IVAR'].data.copy(),
#                                          mask=star_mask), -1.0)
                db['STAR_FLUX'][i] = numpy.ma.mean(star_flux)
                db['STAR_Z'][i] = (numpy.ma.mean(svel_clipped)*(1+db['NSA_Z'][i]) 
                                        + astropy.constants.c.to('km/s').value*db['NSA_Z'][i]) \
                                    / astropy.constants.c.to('km/s').value
            else:
                db['STAR_FLUX'][i] = -9999.0
                db['STAR_Z'][i] = -9999.0

            db['STAR_VEL_LO'][i], db['STAR_VEL_HI'][i] = growth(svel, [0.025,0.975])
            svel_clipped = sigma_clip(svel)
            db['STAR_VEL_LO_CLIP'][i], db['STAR_VEL_HI_CLIP'][i] \
                        = growth(svel_clipped, [0.025,0.975])

            if numpy.sum(~numpy.ma.getmaskarray(_havel)) > 0:
                havel_clipped = sigma_clip(_havel)
                halp_mask = numpy.ma.getmaskarray(havel_clipped)
                halp_flux = numpy.ma.MaskedArray(
                                        dapmaps['EMLINE_GFLUX'].data.copy()[emline['Ha-6564'],:,:],
                                        mask=halp_mask)
#                havel_var = numpy.ma.power(numpy.ma.MaskedArray(
#                                    dapmaps['EMLINE_GVEL_IVAR'].data.copy()[emline['Ha-6564'],:,:],
#                                           mask=halp_mask), -1.0)
                db['HA_FLUX'][i] = numpy.ma.mean(halp_flux)
                db['HA_Z'][i] = (numpy.ma.mean(havel_clipped)*(1+db['NSA_Z'][i]) 
                                        + astropy.constants.c.to('km/s').value*db['NSA_Z'][i]) \
                                    / astropy.constants.c.to('km/s').value
            else:
                db['HA_FLUX'][i] = -9999.0
                db['HA_Z'][i] = -9999.0

            db['HA_VEL_LO'][i], db['HA_VEL_HI'][i] = growth(havel, [0.025,0.975])
            svel_clipped = sigma_clip(havel)
            db['HA_VEL_LO_CLIP'][i], db['HA_VEL_HI_CLIP'][i] = growth(svel_clipped, [0.025,0.975])

            db['HA_VEL_LO'][i], db['HA_VEL_HI'][i] = growth(havel, [0.025,0.975])


            indx = (r_re < 1.0)
            if numpy.sum(indx) > 0:
                numer = numpy.ma.sum(bin_flux[indx]*ssig2corr[indx])
                db['STAR_DISP_1RE'][i] = -9999.0 if numer < 0 \
                                            else numpy.ma.sqrt(numer/numpy.ma.sum(bin_flux[indx]))
                db['SC_RCHI2_1RE'][i] = numpy.ma.median(scchi[indx])

                numer = numpy.ma.sum(line_flux[emline['Ha-6564'],indx]*hasig2corr[indx])
                db['HA_DISP_1RE'][i] = -9999.0 if numer < 0 \
                                            else numpy.ma.sqrt(numer 
                                                / numpy.ma.sum(line_flux[emline['Ha-6564'],indx]))
                db['EML_SB_1RE'][i,:] = numpy.ma.median(line_flux[:,indx].reshape(neml,-1), axis=1)
            else:
                db['STAR_DISP_1RE'][i] = -9999.0
                db['SC_RCHI2_1RE'][i] = -9999.0
                db['HA_DISP_1RE'][i] = -9999.0
                db['EML_SB_1RE'][i,:] = numpy.full(neml, -9999.0, dtype=float)

            lo, hi = growth(hasig2corr, [0.025,0.975])
            db['HA_DISP_HI'][i] = -9999.0 if hi < 0 else numpy.sqrt(hi)
            hasig2corr_clipped = sigma_clip(hasig2corr)
            lo, hi = growth(hasig2corr_clipped, [0.025,0.975])
            db['HA_DISP_HI_CLIP'][i] = -9999.0 if hi < 0 else numpy.sqrt(hi)

            # Get the total flux and maximum surface brightness of the
            # emission lines
            db['EML_FLUX_TOTAL'][i,:] = numpy.ma.sum(line_flux.reshape(neml,-1), axis=1)
            db['EML_SB_PEAK'][i,:] = numpy.ma.amax(line_flux.reshape(neml,-1), axis=1)

            # Get the spectral index statistics
            spec_index_clipped = sigma_clip(spec_index.reshape(nindx,-1), axis=1)
            for j in range(nindx):
                db['SPINDX_LO'][i,j], db['SPINDX_MD'][i,j], db['SPINDX_HI'][i,j] \
                        = growth(spec_index[j,:], [0.025,0.5,0.975])
                db['SPINDX_LO_CLIP'][i,j], db['SPINDX_MD_CLIP'][i,j], db['SPINDX_HI_CLIP'][i,j] \
                        = growth(spec_index_clipped[j,:], [0.025,0.5,0.975])

        print('Processing: DONE           ')

        # Create the primary header
        hdr = fits.Header()
        hdr['DATE'] = (time.strftime('%Y-%m-%d',time.gmtime()), 'UTC date created')
        hdr['DRPVER'] = (self.drpver, 'DRP version')
        hdr['DAPVER'] = (self.dapver, 'DAP version')
        hdr['DAPALLV'] = (self.version, 'DAPall code version')

#        print(db['EML_NAME'].shape)
#        print(db['EML_NAME'].dtype)
#        print(db['EML_NAME'].dtype.str)
#        print(rec_to_fits_type(db['EML_NAME']))
#        exit()
#        print(rec_to_fits_type(db['EML_NAME']))
#        print('({0},{1})'.format(neml,self.ndap))
#        test = fits.Column(name='EML_NAME', format=rec_to_fits_type(db['EML_NAME']),
#                                      array=db['EML_NAME'], dim='({0},{1})'.format(neml,self.ndap))
#        print(test)
#        exit()

        # Write the file
        self.hdu = fits.HDUList([ fits.PrimaryHDU(header=hdr),
                                  fits.BinTableHDU.from_columns(
                                    [ fits.Column(name=n, format=rec_to_fits_type(db[n]),
                                      array=db[n]) for n in db.dtype.names ],
                                    name='DAPALL') ])
        write_hdu(self.hdu, self.file_path(), clobber=True, checksum=True)
        

