# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
A class hierarchy that fits the emission lines.

*License*:
    Copyright (c) 2015, SDSS-IV/MaNGA Pipeline Group
        Licensed under BSD 3-clause license - see LICENSE.rst

*Source location*:
    $MANGADAP_DIR/python/mangadap/proc/emissionlinefits.py

*Imports and python version compliance*:
    ::

        from __future__ import division
        from __future__ import print_function
        from __future__ import absolute_import
        from __future__ import unicode_literals

        import sys
        import warnings
        if sys.version > '3':
            long = int
        
        import glob
        import os
        import logging
        
        import numpy
        from astropy.io import fits
        
        from ..drpfits import DRPFits
        from ..par.parset import ParSet
        from ..par.artifactdb import ArtifactDB
        from ..par.emissionlinedb import EmissionLineDB
        from ..config.defaults import dap_source_dir, default_dap_file_name
        from ..config.defaults import default_dap_method, default_dap_method_path
        from ..util.fitsutil import DAPFitsUtil
        from ..util.fileio import init_record_array, rec_to_fits_type, rec_to_fits_col_dim
        from ..util.bitmask import BitMask
        from ..util.pixelmask import SpectralPixelMask
        from ..util.log import log_output
        from .spatiallybinnedspectra import SpatiallyBinnedSpectra
        from .stellarcontinuummodel import StellarContinuumModel
        from .elric import Elric, ElricPar, GaussianLineProfile
        from .util import select_proc_method

*Class usage examples*:
    Add examples!

*Revision history*:
    | **26 Apr 2016**: Implementation begun by K. Westfall (KBW)
    | **28 Jul 2016**: (KBW) Fixed error in initialization of guess
        redshift when stellar continuum is provided.
    | **23 Feb 2017**: (KBW) Use DAPFitsUtil read and write functions.
    | **31 May 2017**: (KBW) Revert to using
        :class:`mangadap.proc.stellarcontinuummodel.StellarContinuumModel`
        on input
    | **08 Sep 2017**: (KBW) Add `deconstruct_bins` flag to parameters.
    | **02 Feb 2018**: (KBW) Use
        :func:`mangadap.proc.spectralfitting.EmissionLineFit.select_binned_spectra_to_fit`.
    | **15 Feb 2018**: (KBW) No longer imports
        :class:`mangadap.proc.emissionlinemoments.EmissionLineMoments`
        to avoid circular imports (this should be coded better...)
    | **24 Feb 2018**: (KBW) Added keyword for a new template library
        that can be used instead of the same templates used during the
        stellar-continuum fit.
    
.. _astropy.io.fits.hdu.hdulist.HDUList: http://docs.astropy.org/en/v1.0.2/io/fits/api/hdulists.html
.. _glob.glob: https://docs.python.org/3.4/library/glob.html
.. _logging.Logger: https://docs.python.org/3/library/logging.html
.. _numpy.ma.MaskedArray: http://docs.scipy.org/doc/numpy-1.10.1/reference/maskedarray.baseclass.html#numpy.ma.MaskedArray


"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import sys
import warnings
if sys.version > '3':
    long = int

import glob
import os
import logging

import numpy

from scipy import interpolate

from astropy.io import fits
import astropy.constants

from ..drpfits import DRPFits
from ..par.parset import ParSet
from ..par.artifactdb import ArtifactDB
from ..par.emissionlinedb import EmissionLineDB
from ..config.defaults import dap_source_dir, default_dap_file_name
from ..config.defaults import default_dap_method, default_dap_method_path
from ..util.fitsutil import DAPFitsUtil
from ..util.fileio import init_record_array, rec_to_fits_type, rec_to_fits_col_dim
from ..util.bitmask import BitMask
from ..util.pixelmask import SpectralPixelMask
from ..util.log import log_output
from ..util.parser import DefaultConfig
from .spatiallybinnedspectra import SpatiallyBinnedSpectra
from .stellarcontinuummodel import StellarContinuumModel
#from .emissionlinemoments import EmissionLineMoments
from .bandpassfilter import emission_line_equivalent_width
from .elric import Elric, ElricPar
from .sasuke import Sasuke, SasukePar
from .util import select_proc_method, replace_with_data_from_nearest_coo
from .spectralfitting import EmissionLineFit

from matplotlib import pyplot

class EmissionLineModelDef(ParSet):
    """
    A class that holds the parameters necessary to perform the
    emission-line profile fits.

    .. todo::
        Include waverange?

    Args:
        key (str): Keyword used to distinguish between different
            emission-line moment databases.
        minimum_snr (bool): Minimum S/N of spectrum to fit
        artifacts (str): String identifying the artifact database to use
        emission_lines (str): String identifying the emission-line
            database to use
        continuum_tpl_key (str): String identifying the continuum
            templates to use
    """
    def __init__(self, key, minimum_snr, deconstruct_bins, mom_vel_name, mom_disp_name,
                 artifacts, emission_lines, continuum_tpl_key, fitpar, fitclass, fitfunc):
        in_fl = [ int, float ]
        par_opt = [ ParSet, dict ]

        pars =     [ 'key', 'minimum_snr', 'deconstruct_bins', 'mom_vel_name', 'mom_disp_name',
                     'artifacts', 'emission_lines', 'continuum_tpl_key', 'fitpar', 'fitclass',
                     'fitfunc' ]
        values =   [ key, minimum_snr, deconstruct_bins, mom_vel_name, mom_disp_name, artifacts,
                     emission_lines, continuum_tpl_key, fitpar, fitclass, fitfunc ]
        dtypes =   [ str, in_fl, bool, str, str, str, str, str, par_opt, None, None ]
        can_call = [ False, False, False, False, False, False, False, False, False, False, True ]

        ParSet.__init__(self, pars, values=values, dtypes=dtypes)


def validate_emission_line_modeling_method_config(cnfg):
    """ 
    Validate the :class:`mangadap.util.parser.DefaultConfig` with the
    emission-line modeling method parameters.

    Args:
        cnfg (:class:`mangadap.util.parser.DefaultConfig`): Object meant
            to contain defining parameters of the emission-line modeling
            method needed by :class:`EmissionLineModelDef`

    Raises:
        KeyError: Raised if any required keywords do not exist.
        ValueError: Raised if keys have unacceptable values.
    """
    # Check for required keywords
    if 'key' not in cnfg:
        raise KeyError('Keyword \'key\' must be provided.')


def available_emission_line_modeling_methods(dapsrc=None):
    """
    Return the list of available emission-line modeling methods.

    Available methods are:

    .. todo::
        Fill in

    Args:
        dapsrc (str): (**Optional**) Root path to the DAP source
            directory.  If not provided, the default is defined by
            :func:`mangadap.config.defaults.dap_source_dir`.

    Returns:
        list: A list of :func:`EmissionLineModelDef` objects, each
        defining an emission-line modeling method.

    Raises:
        NotADirectoryError: Raised if the provided or default
            *dapsrc* is not a directory.
        OSError/IOError: Raised if no emission-line moment configuration
            files could be found.
        KeyError: Raised if the emission-line modeling method keywords
            are not all unique.
        NotImplementedError: Raised if the method requests the
            deconstruction of the bins into spaxls for using the Elric
            fitter.

    .. todo::
        Possible to add a python call that reads the databases and
        constructs the table for presentation in sphinx so that the text
        above doesn't have to be edited with changes in the available
        databases?

    """
    # Check the source directory exists
    dapsrc = dap_source_dir() if dapsrc is None else str(dapsrc)
    if not os.path.isdir(dapsrc):
        raise NotADirectoryError('{0} does not exist!'.format(dapsrc))

    # Check the configuration files exist
    search_dir = os.path.join(dapsrc, 'python/mangadap/config/emission_line_modeling')
    ini_files = glob.glob(os.path.join(search_dir, '*.ini'))
    if len(ini_files) == 0:
        raise IOError('Could not find any configuration files in {0} !'.format(search_dir))

    # Build the list of method definitions
    method_list = []
    for f in ini_files:
        # Read the config file
        cnfg = DefaultConfig(f)
        # Ensure it has the necessary elements to define the template
        # library
        validate_emission_line_modeling_method_config(cnfg)
        deconstruct_bins = cnfg.getbool('deconstruct_bins', default=False)
        minimum_snr = cnfg.getfloat('minimum_snr', default=0.0)
        continuum_tpl_key = cnfg.get('continuum_templates')

        if cnfg['fit_method'] == 'elric':
            # Chose to use Elric: Parameter set has defaults to handle
            # missing or None values for baseline_order, window_buffer,
            # and minimum_snr
            if deconstruct_bins:
                raise NotImplementedError('When using Elric, cannot deconstruct bins into spaxels.')
            if continuum_tpl_key is not None:
                raise NotImplementedError('When using Elric, cannot change continuum templates.')

            fitpar = ElricPar(None, cnfg.getint('baseline_order'), cnfg.getfloat('window_buffer'),
                              None, None, minimum_snr=minimum_snr)
            fitclass = Elric(EmissionLineModelBitMask(dapsrc=dapsrc))
            fitfunc = fitclass.fit_SpatiallyBinnedSpectra

        elif cnfg['fit_method'] == 'sasuke':
            # Chose to use Sasuke: Parameter set has defaults to handle
            # missing or None values for reject_boxcar, bias, moments,
            # degree, mdegree; if provided new continuum templates are
            # constructed during the _fill_method_par call.
            fitpar = SasukePar(None, None, continuum_templates=continuum_tpl_key,
                               etpl_line_sigma_mode=cnfg.get('etpl_line_sigma_mode'),
                               etpl_line_sigma_min=cnfg.getfloat('etpl_line_sigma_min'),
                               velscale_ratio=cnfg.getint('velscale_ratio'),
                               guess_redshift=None, guess_dispersion=None, minimum_snr=minimum_snr,
                               deconstruct_bins=deconstruct_bins, pixelmask=None,
                               reject_boxcar=cnfg.getint('reject_boxcar'),
                               bias=cnfg.getfloat('bias'), moments=cnfg.getint('moments'),
                               degree=cnfg.getint('degree'), mdegree=cnfg.getint('mdegree'),
                               reddening=cnfg.getfloat('internal_reddening'))
            fitclass = Sasuke(EmissionLineModelBitMask(dapsrc=dapsrc))
            fitfunc = fitclass.fit_SpatiallyBinnedSpectra

        method_list += [ EmissionLineModelDef(cnfg['key'], minimum_snr, deconstruct_bins,
                                              cnfg.get('mom_vel_name'), cnfg.get('mom_disp_name'),
                                              cnfg.get('artifact_mask'), cnfg['emission_lines'],
                                              continuum_tpl_key, fitpar, fitclass, fitfunc) ]

    # Check the keywords of the libraries are all unique
    if len(numpy.unique(numpy.array([method['key'] for method in method_list]))) \
                != len(method_list):
        raise KeyError('Emission-line fitting method keywords are not all unique!')

    # Return the default list of fitting methods
    return method_list


class EmissionLineModelBitMask(BitMask):
    r"""

    Derived class that specifies the mask bits for the emission-line
    modeling.  See :class:`mangadap.util.bitmask.BitMask` for
    attributes.

    A list of the bits and meanings are provided by the base class
    function :func:`mangadap.util.bitmask.BitMask.info`; i.e.,::

        from mangadap.proc.emissionlinemodel import EmissionLineModelBitMask
        bm = EmissionLineModelBitMask()
        bm.info()

    """
    def __init__(self, dapsrc=None):
        dapsrc = dap_source_dir() if dapsrc is None else str(dapsrc)
        BitMask.__init__(self, ini_file=os.path.join(dapsrc, 'python', 'mangadap', 'config',
                                                     'bitmasks', 'emission_line_model_bits.ini'))


class EmissionLineModel:
    r"""
    Class that holds the emission-line model fits.

    Args:
        loggers (list): (**Optional**) List of `logging.Logger`_ objects
            to log progress; ignored if quiet=True.  Logging is done
            using :func:`mangadap.util.log.log_output`.  Default is no
            logging.
        quiet (bool): (**Optional**) Suppress all terminal and logging
            output.  Default is False.

    Attributes:
        loggers (list): List of `logging.Logger`_ objects to log
            progress; ignored if quiet=True.  Logging is done using
            :func:`mangadap.util.log.log_output`.
        quiet (bool): Suppress all terminal and logging output.

    """
#    @profile
    def __init__(self, method_key, binned_spectra, stellar_continuum=None,
                 emission_line_moments=None, redshift=None, dispersion=None, minimum_error=None,
                 method_list=None, artifact_list=None, emission_line_db_list=None, dapsrc=None,
                 dapver=None, analysis_path=None, directory_path=None, output_file=None,
                 hardcopy=True, clobber=False, checksum=False, loggers=None, quiet=False):

        self.loggers = None
        self.quiet = False

        # Define the database properties
        self.method = None
        self.pixelmask = None
        self.artdb = None
        self.emldb = None
        self._define_method(method_key, method_list=method_list, artifact_list=artifact_list,
                            emission_line_db_list=emission_line_db_list, dapsrc=dapsrc)

        self.neml = None if self.emldb is None else self.emldb.neml

        self.binned_spectra = None
        self.redshift = None
        self.dispersion = None
        self.stellar_continuum = None

        # Define the output directory and file
        self.directory_path = None      # Set in _set_paths
        self.output_file = None
        self.hardcopy = None

        # Initialize the objects used in the assessments
        self.bitmask = EmissionLineModelBitMask(dapsrc=dapsrc)

        self.hdu = None
        self.checksum = checksum
        self.shape = None
        self.spatial_shape = None
#       nspec is never used and it won't be correct if the bins are
#       deconstructed.
#        self.nspec = None
        self.spatial_index = None
        self.spectral_arrays = None
        self._assign_spectral_arrays()
        self.image_arrays = None
        self._assign_image_arrays()
        self.dispaxis = None
        self.nwave = None

        self.nmodels = None
        self.missing_models = None

        # Fit the binned spectra
        self.fit(binned_spectra, stellar_continuum=stellar_continuum,
                 emission_line_moments=emission_line_moments, redshift=redshift,
                 dispersion=dispersion, minimum_error=minimum_error, dapsrc=dapsrc, dapver=dapver,
                 analysis_path=analysis_path, directory_path=directory_path,
                 output_file=output_file, hardcopy=hardcopy, clobber=clobber, loggers=loggers,
                 quiet=quiet)


#    def __del__(self):
#        """
#        Deconstruct the data object by ensuring that the fits file is
#        properly closed.
#        """
#        if self.hdu is None:
#            return
#        self.hdu.close()
#        self.hdu = None


    def __getitem__(self, key):
        return self.hdu[key]


    def _define_method(self, method_key, method_list=None, artifact_list=None,
                       emission_line_db_list=None, dapsrc=None):
        r"""
        Select the modeling method
        """
        # Grab the specific method
        self.method = select_proc_method(method_key, EmissionLineModelDef, method_list=method_list,
                                         available_func=available_emission_line_modeling_methods,
                                         dapsrc=dapsrc)

        # Instantiate the artifact and emission-line databases
        self.artdb = None if self.method['artifacts'] is None else \
                    ArtifactDB(self.method['artifacts'], artdb_list=artifact_list, dapsrc=dapsrc)
        self.pixelmask = SpectralPixelMask(artdb=self.artdb)

        self.emldb = None if self.method['emission_lines'] is None else \
                        EmissionLineDB(self.method['emission_lines'],
                                       emldb_list=emission_line_db_list, dapsrc=dapsrc)


    def _set_paths(self, directory_path, dapver, analysis_path, output_file):
        """
        Set the :attr:`directory_path` and :attr:`output_file`.  If not
        provided, the defaults are set using, respectively,
        :func:`mangadap.config.defaults.default_dap_method_path` and
        :func:`mangadap.config.defaults.default_dap_file_name`.

        Args:
            directory_path (str): The exact path to the DAP
                emission-line moments file.  See :attr:`directory_path`.
            dapver (str): DAP version.
            analysis_path (str): The path to the top-level directory
                containing the DAP output files for a given DRP and DAP
                version.
            output_file (str): The name of the file with emission-line
                moment measurements.  See :func:`measure`.
        """
        # Set the output directory path
        continuum_method = None if self.stellar_continuum is None \
                                else self.stellar_continuum.method['key']
        method = default_dap_method(binned_spectra=self.binned_spectra,
                                    continuum_method=continuum_method)
        self.directory_path = default_dap_method_path(method, plate=self.binned_spectra.drpf.plate,
                                                      ifudesign=self.binned_spectra.drpf.ifudesign,
                                                      ref=True,
                                                      drpver=self.binned_spectra.drpf.drpver,
                                                      dapver=dapver, analysis_path=analysis_path) \
                                        if directory_path is None else str(directory_path)

        # Set the output file
        ref_method = '{0}-{1}'.format(self.binned_spectra.rdxqa.method['key'],
                                      self.binned_spectra.method['key'])
        if continuum_method is not None:
            ref_method = '{0}-{1}'.format(ref_method, continuum_method)
        ref_method = '{0}-{1}'.format(ref_method, self.method['key'])
        self.output_file = default_dap_file_name(self.binned_spectra.drpf.plate,
                                                 self.binned_spectra.drpf.ifudesign, ref_method) \
                                        if output_file is None else str(output_file)


    def _initialize_primary_header(self, hdr=None):
        """

        Initialize the header.

        """
        # Copy the from the DRP and clean it
        if hdr is None:
            hdr = self.binned_spectra.drpf.hdu['PRIMARY'].header.copy()
            hdr = DAPFitsUtil.clean_dap_primary_header(hdr)
        
        hdr['AUTHOR'] = 'Kyle B. Westfall <westfall@ucolick.org>'
        hdr['ELFKEY'] = (self.method['key'], 'Emission-line modeling method keyword')
        hdr['ELFMINSN'] = (self.method['minimum_snr'], 'Minimum S/N of spectrum to include')
        hdr['ARTDB'] = (self.method['artifacts'], 'Artifact database keyword')
        hdr['EMLDB'] = (self.method['emission_lines'], 'Emission-line database keyword')
        if self.stellar_continuum is not None:
            hdr['SCKEY'] = (self.stellar_continuum.method['key'], 'Stellar-continuum model keyword')
        hdr['NELMOD'] = (self.nmodels, 'Number of unique emission-line models')
        return hdr


    def _emission_line_database_dtype(self, name_len, mode_len, prof_len):
        r"""
        Construct the record array data type for the output fits
        extension.
        """
        return [ ('ID', int),
                 ('NAME','<U{0:d}'.format(name_len)),
                 ('RESTWAVE', float),
                 ('ACTION', '<U1'),
                 ('FLUXRATIO', float),
                 ('MODE','<U{0:d}'.format(mode_len)),
                 ('PROFILE', '<U{0:d}'.format(prof_len)),
                 ('NCOMP', int)
               ]


    def _compile_database(self, quiet=False):
        """
        Compile the database with the specifications of each index.
        """
        if self.emldb is None:
            try:
                # Try to use member functions of the fitting class to
                # set at least the NAME and RESTWAVE of the line for use
                # in dapfits.py
                line_database = self.method['fitclass'].line_database()
            except AttributeError:
                if not quiet:
                    warnings.warn('Emission-line fit did not use an EmissionLineDB object, and '
                                  'does not have a line_database() method.  Will not be able to '
                                  'use output to identify which lines were fit.')
                line_database = None
            return line_database

        name_len = numpy.amax([ len(n) for n in self.emldb['name'] ])
        mode_len = numpy.amax([ len(m) for m in self.emldb['mode'] ])
        prof_len = numpy.amax([ len(p) for p in self.emldb['profile'] ])

        # Instatiate the table data that will be saved defining the set
        # of emission-line moments measured
        line_database = init_record_array(self.emldb.neml,
                                self._emission_line_database_dtype(name_len, mode_len, prof_len))

        hk = [ 'ID', 'NAME', 'RESTWAVE', 'ACTION', 'FLUXRATIO', 'MODE', 'PROFILE', 'NCOMP' ]
        mk = [ 'index', 'name', 'restwave', 'action', 'flux', 'mode', 'profile', 'ncomp' ]
        for _hk, _mk in zip(hk,mk):
            line_database[_hk] = self.emldb[_mk]

        return line_database


    def _assign_input_kinematics(self, emission_line_moments, redshift, dispersion,
                                 default_dispersion=100.0):
        """
        Set the initial redshift and velocity dispersion for each
        spectrum.
        
        In terms of precedence, directly provided redshifts and
        dispersions override those in the EmissionLineMoments model,
        which overrides those in the StellarContinuumModel object.  The
        default_dispersion does *not* take precedence over *any*
        provided disperison.

        If emission_line_moments, redshift, and
        :attr:`stellar_continuum` are all None, the redshift is set to
        0.0.  If dispersion is also None, the dispersion is set to
        default_dispersion.

        To get the stellar kinematics, the function calls
        :func:`mangadap.proc.stellarcontinuummodel.StellarContinuumModel.matched_kinematics`.
        In this function, the provided redshift and dispersion must be a
        single value or None; therefore, the means of any vectors
        provided as redshift or disperison are passsed to this function
        instead of the full vector.

        Args:
            emission_line_moments
                (:class:`mangadap.proc.emissionlinemoments.EmissionLineMoments`):
                Object with the results of the emission-line-moment
                measurements
            redshift (float, numpy.ndarray): Redshifts (:math:`z`) to
                use for each spectrum.
            dispersion (float, numpy.ndarray): Velocity dispersion
                (km/s) to use for each spectrum.
            default_dispersion (float, numpy.ndarray): (**Optional**)
                Default velocity dispersion to use (km/s), if relevant.
        """

        #---------------------------------------------------------------
        # Get the redshift and dispersion measured for the stars or
        # based on the moments if the relevant objects are provided
        obj_redshift, obj_dispersion = (None, None)
        if emission_line_moments is not None:
            names = emission_line_moments.channel_names()
            vel_channel = -1 if self.method['mom_vel_name'] is None \
                                else numpy.arange(len(names))[self.method['mom_vel_name']
                                                                == names][0]
            if vel_channel > -1:
                obj_redshift = numpy.full(self.binned_spectra.nbins, 0.0, dtype=float)
            
                mom1_masked = emission_line_moments['ELMMNTS'].data['MASK'][:,vel_channel] > 0
                # - Use the 1st moment of the named line
                obj_redshift[ emission_line_moments['ELMMNTS'].data['BINID_INDEX'] ] \
                        = emission_line_moments['ELMMNTS'].data['MOM1'][:,vel_channel] \
                                / astropy.constants.c.to('km/s').value

                # - For missing bins in the moment measurements and bad
                # line moment measurements, use the value for the
                # nearest good bin
                bad_bins = numpy.append(emission_line_moments.missing_bins,
                                        emission_line_moments['ELMMNTS'].data['BINID']\
                                            [mom1_masked]).astype(int)
                if len(bad_bins) > 0 and len(bad_bins) != self.binned_spectra.nbins:
                    nearest_good_bin_index = self.binned_spectra.find_nearest_bin(bad_bins,
                                                                                  indices=True)
                    bad_bin_index = self.binned_spectra.get_bin_indices(bad_bins)
                    obj_redshift[bad_bin_index] = obj_redshift[nearest_good_bin_index]

            disp_channel = -1 if self.method['mom_disp_name'] is None \
                                else numpy.where(self.method['mom_disp_name'] == names)[0]
            if disp_channel > -1:
                obj_dispersion = numpy.full(self.binned_spectra.nbins, default_dispersion,
                                            dtype=float)
            
                mom2_masked = emission_line_moments['ELMMNTS'].data['MASK'][:,disp_channel] > 0
                # - Use the 2nd moment of the named line
                obj_dispersion[ emission_line_moments['ELMMNTS'].data['BINID_INDEX'] ] \
                                = emission_line_moments['ELMMNTS'].data['MOM2'][:,disp_channel]

                # - For missing bins in the moment measurements and bad
                # line moment measurements, use the value for the
                # nearest good bin
                bad_bins = numpy.append(emission_line_moments.missing_bins,
                                        emission_line_moments['ELMMNTS'].data['BINID']\
                                            [mom2_masked]).astype(int)
                if len(bad_bins) > 0 and len(bad_bins) != self.binned_spectra.nbins:
                    nearest_good_bin_index = self.binned_spectra.find_nearest_bin(bad_bins,
                                                                                  indices=True)
                    bad_bin_index = self.binned_spectra.get_bin_indices(bad_bins)
                    obj_dispersion[bad_bin_index] = obj_dispersion[nearest_good_bin_index]

        if self.stellar_continuum is not None and (obj_redshift is None or obj_dispersion is None):
            sc_redshift, sc_dispersion = \
                    self.stellar_continuum.matched_kinematics(self.binned_spectra['BINID'].data,
                                redshift=None if redshift is None else numpy.ma.mean(redshift),
                                dispersion=default_dispersion if dispersion is None 
                                                                else numpy.ma.mean(dispersion),
                                corrected=True, nearest=True,
                                missing=self.binned_spectra.missing_bins)
            if obj_redshift is None:
                obj_redshift = sc_redshift
            if obj_dispersion is None:
                obj_dispersion = sc_dispersion
        #---------------------------------------------------------------

        # Redshift: use the stellar continuum values if present,
        # otherwise set the default to 0.
        self.redshift = obj_redshift if redshift is None else \
                            numpy.zeros(self.binned_spectra.nbins, dtype=float)
        if redshift is not None:
            _redshift = numpy.atleast_1d(redshift)
            if len(_redshift) not in [ 1, self.binned_spectra.nbins ]:
                raise ValueError('Provided redshift must be either a single value or match the ' \
                                 'number of binned spectra.')
            self.redshift = numpy.full(self.binned_spectra.nbins, redshift, dtype=float) \
                                if len(_redshift) == 1 else _redshift.copy()

        # Same for dispersion
        # Redshift: use the stellar continuum values if present,
        # otherwise set the default to 0.
        self.dispersion = obj_dispersion if dispersion is None else \
                            numpy.full(self.binned_spectra.nbins, default_dispersion, dtype=float)
        if dispersion is not None:
            _dispersion = numpy.atleast_1d(dispersion)

            if len(_dispersion) not in [ 1, self.binned_spectra.nbins ]:
                raise ValueError('Provided dispersion must be either a single value or match the ' \
                                 'number of binned spectra.')
            self.dispersion = numpy.full(self.binned_spectra.nbins, dispersion, dtype=float) \
                                if len(_dispersion) == 1 else _dispersion.copy()


    def _fill_method_par(self, dapsrc=None, analysis_path=None):
        """

        Fill in any remaining modeling parameters.

        .. todo::
            - Construct the replacement template library here instead of
              in Sasuke?
        """
        # Fit parameters not defined so continue
        if self.method['fitpar'] is None:
            return

        # Fill the guess kinematics
        if 'guess_redshift' in self.method['fitpar'].keys():
            self.method['fitpar']['guess_redshift'] = self.redshift
        if 'guess_dispersion' in self.method['fitpar'].keys():
            self.method['fitpar']['guess_dispersion'] = self.dispersion
        if 'emission_lines' in self.method['fitpar'].keys():
            self.method['fitpar']['emission_lines'] = self.emldb
        if 'pixelmask' in self.method['fitpar'].keys():
            self.method['fitpar']['pixelmask'] = self.pixelmask
        if self.stellar_continuum is not None \
                and 'stellar_continuum' in self.method['fitpar'].keys():
            self.method['fitpar']['stellar_continuum'] = self.stellar_continuum


    def _add_method_header(self, hdr, model_binid=None):
        """Add fitting method information to the header."""
        if self.method['fitclass'] is not None:
            try:
                hdr['ELTYPE'] = self.method['fitclass'].fit_type
                hdr['ELMETH'] = self.method['fitclass'].fit_method
            except:
                if not self.quiet and self.hardcopy:
                    warnings.warn('Fit class object does not have fit_type and/or fit_method ' \
                                  'attributes.  No parameters written to header.')
        if self.method['fitpar'] is not None:
            try:
                hdr = self.method['fitpar'].toheader(hdr)
            except:
                if not self.quiet and self.hardcopy:
                    warnings.warn('Fit parameter class has no toheader() function.  No ' \
                                  'parameters written to header.')

        # If the model bin IDs are provided, assume the data was
        # rebinned, meaning that the bin ID of this object is not
        # directly tied to the bin ID for the spatially binned spectra
        # object (self.binned_spectra)
        hdr['ELREBIN'] = (model_binid is not None, 'Bin IDs disconnected from SC binning')

        return hdr


    def _finalize_cube_mask(self, mask):
        """
        Finalize the mask by setting the DIDNOTUSE, FORESTAR, and LOW_SNR masks

        Returns:
            numpy.ndarray : Bitmask array.
        """
        # Turn on the flag stating that the pixel wasn't used
        indx = self.binned_spectra.bitmask.flagged(self.binned_spectra.drpf['MASK'].data,
                                                   flag=self.binned_spectra.do_not_fit_flags())
        mask[indx] = self.bitmask.turn_on(mask[indx], 'DIDNOTUSE')

        # Turn on the flag stating that the pixel has a foreground star
        indx = self.binned_spectra.bitmask.flagged(self.binned_spectra.drpf['MASK'].data,
                                                   flag='FORESTAR')
        mask[indx] = self.bitmask.turn_on(mask[indx], 'FORESTAR')

        # Turn on the flag stating that an individual bin/spaxel was not
        # used, anything without a non-negative binid, in the fit
        indx = self['BINID'].data < 0
        mask[indx,:] = self.bitmask.turn_on(mask[indx,:], 'DIDNOTUSE')

        # Turn on the flag stating that an individual bin/spaxel was not
        # used because the S/N was too low
        indx = (self['BINID'].data < 0) & (self.binned_spectra['BINID'].data > -1)
        mask[indx,:] = self.bitmask.turn_on(mask[indx,:], 'LOW_SNR')

#        # Turn on the flag stating that the S/N in the spectrum was
#        # below the requested limit
#        low_snr = numpy.invert(self.binned_spectra['BINID'].data == self.hdu['BINID'].data)
#        indx = numpy.array([low_snr]*self.nwave).transpose(1,2,0)
#        mask[indx] = self.bitmask.turn_on(mask[indx], flag='LOW_SNR')

        return mask

    
    def _assign_spectral_arrays(self):
        self.spectral_arrays = [ 'FLUX', 'BASE', 'MASK' ]


    def _assign_image_arrays(self):
        """
        Set :attr:`image_arrays`, which contains the list of extensions
        in :attr:`hdu` that are on-sky image data.
        """
        self.image_arrays = [ 'BINID' ]


    def _get_missing_models(self, unique_bins=None):
        if unique_bins is None:
            good_snr = self.binned_spectra.above_snr_limit(self.method['minimum_snr'])
            return numpy.sort(self.binned_spectra.missing_bins + 
                        self.binned_spectra['BINS'].data['BINID'][numpy.invert(good_snr)].tolist())
        return SpatiallyBinnedSpectra._get_missing_bins(unique_bins)


    def _get_line_fit_metrics(self, model_flux, model_base, model_mask, model_eml_par,
                              model_binid, metric_window=15):
        """
        Compute some line fit metrics.

        Fills in the LINE_PIXC, AMP, ANR, LINE_NSTAT, LINE_RMS,
        LINE_FRMS, and LINE_CHI2 elements of model_eml_par.

        metric_window sets the size in pixels to use for the metrics

        """

        # Try to construct the line database
        line_database = self._compile_database()
        if line_database is None:
            warnings.warn('Cannot construct line metrics.')
            return model_eml_par
        neml = len(line_database)

#        print(neml)
#        print(line_database['RESTWAVE'])

        # Construct the bin ID map
        bin_indx = DAPFitsUtil.downselect_bins(self.binned_spectra['BINID'].data.ravel(),
                                            model_eml_par['BINID']).reshape(self.spatial_shape) \
                            if model_binid is None else model_binid

#        pyplot.imshow(bin_indx.T, origin='lower', interpolation='nearest')
#        pyplot.colorbar()
#        pyplot.show()

        wave = self.binned_spectra['WAVE'].data

        # Get the baseline continuum
        continuum = numpy.ma.MaskedArray(numpy.zeros(model_base.shape, dtype=float)) \
                        if self.stellar_continuum is None else \
                    self.stellar_continuum.fill_to_match(bin_indx, missing=self.missing_models)
        indx = numpy.unique(bin_indx.ravel())[1:]
        continuum = continuum[indx,:] + model_base

        # Allow the fitting function to return a boolean model mask
        if model_mask.dtype is numpy.dtype('bool'):
            _model_mask = numpy.zeros(model_flux.shape, dtype=self.bitmask.minimum_dtype())
            _model_mask[model_mask] = self.bitmask.turn_on(_model_mask[model_mask], 'DIDNOTUSE')
        else:
            _model_mask = model_mask

        # Get the fitted spectra and errors
        if self.method['deconstruct_bins']:
            # Get the individual spaxels
            _, flux, ferr = EmissionLineFit.get_spectra_to_fit(self.binned_spectra.drpf,
                                                               pixelmask=self.pixelmask,
                                                               error=True)
            # Apply the Galactic extinction
            flux, ferr = self.binned_spectra.galext.apply(flux, err=ferr, deredden=True)

            spaxel_to_fit = self.binned_spectra.check_fgoodpix()
            flux = flux[spaxel_to_fit,:]
            ferr = ferr[spaxel_to_fit,:]
        else:
            # Get the binned spectra
            _, flux, ferr = EmissionLineFit.get_spectra_to_fit(self.binned_spectra,
                                                               pixelmask=self.pixelmask,
                                                               error=True)
            bins_to_fit = EmissionLineFit.select_binned_spectra_to_fit(self.binned_spectra,
                                                        minimum_snr=self.method['minimum_snr'],
                                                        stellar_continuum=self.stellar_continuum)
            flux = flux[bins_to_fit,:]
            ferr = ferr[bins_to_fit,:]

        # Get the continuum-subtracted flux
        flux_nc = flux - continuum

        # Mask the model flux
        _model_flux = numpy.ma.MaskedArray(model_flux+continuum, mask=_model_mask>0)

        nspec = continuum.shape[0]
#        print(nspec, flux.shape[0])
#        pyplot.plot(wave, flux[nspec//2,:], color='k', lw=1.5)
#        pyplot.plot(wave, _model_flux[nspec//2,:], color='C3')
#        pyplot.plot(wave, flux_nc[nspec//2,:], color='k', lw=1.5)
#        pyplot.plot(wave, _model_flux[nspec//2,:]-continuum[nspec//2,:], color='C0')
#        pyplot.show()

        # Get the data for the metrics
        resid = numpy.square(flux-_model_flux)
        fresid = numpy.square(numpy.ma.divide(flux-_model_flux, _model_flux))
        chisqr = numpy.square(numpy.ma.divide(flux-_model_flux, ferr))
        spec_mask = resid.mask | fresid.mask | chisqr.mask
        resid[spec_mask] = numpy.ma.masked
        fresid[spec_mask] = numpy.ma.masked
        chisqr[spec_mask] = numpy.ma.masked

#        pyplot.plot(wave, flux[nspec//2,:], color='k', lw=1.5)
#        pyplot.plot(wave, resid[nspec//2,:], color='0.5', lw=1.5)
#        pyplot.show()
#        pyplot.plot(wave, flux[nspec//2,:], color='k', lw=1.5)
#        pyplot.plot(wave, fresid[nspec//2,:], color='0.5', lw=1.5)
#        pyplot.show()
#        pyplot.plot(wave, flux[nspec//2,:], color='k', lw=1.5)
#        pyplot.plot(wave, chisqr[nspec//2,:], color='0.5', lw=1.5)
#        pyplot.show()

        # Get the pixel at which to center the metric calculations
        flags = ['INSUFFICIENT_DATA', 'FIT_FAILED', 'UNDEFINED_COVAR', 'NEAR_BOUND' ]
        z = numpy.ma.MaskedArray(model_eml_par['KIN'][:,:,0] / astropy.constants.c.to('km/s').value,
                                 mask=self.bitmask.flagged(model_eml_par['MASK'], flag=flags))
        sample_wave = line_database['RESTWAVE'][None,:]*(1+z)
        interp = interpolate.interp1d(wave, numpy.arange(wave.size), bounds_error=False,
                                      fill_value=-1, assume_sorted=True)
        model_eml_par['LINE_PIXC'] = numpy.around(interp(sample_wave.data)).astype(int)
        model_eml_par['LINE_PIXC'][sample_wave.mask] = -1
        mask = (model_eml_par['LINE_PIXC'] < 0) | (model_eml_par['LINE_PIXC'] >= wave.size)

#        print(sample_wave.shape)
#        print(sample_wave[nspec//2,:])
#        print(interp(sample_wave.data[nspec//2,:]))
#
#        print(len(wave), flux.shape[1])
#        print(model_eml_par['LINE_PIXC'].shape)
#        print(model_eml_par['LINE_PIXC'][nspec//2,:])

#        pyplot.plot(wave, flux[nspec//2,:], color='k', lw=1.5)
#        pyplot.plot(wave, _model_flux[nspec//2,:], color='C3')
#        pyplot.scatter(sample_wave[nspec//2,:], [numpy.amax(flux[nspec//2,:])]*neml, marker='.',
#                       s=50, color='C0')
#        pyplot.show()
#
#        pyplot.plot(numpy.arange(len(wave)), flux[nspec//2,:], color='k', lw=1.5)
#        pyplot.plot(numpy.arange(len(wave)), _model_flux[nspec//2,:], color='C3')
#        pyplot.scatter(model_eml_par['LINE_PIXC'][nspec//2,:],
#                       [numpy.amax(flux[nspec//2,:])]*neml, marker='.', s=50, color='C0')
#        pyplot.show()

        # Iterate over the number of lines
        nspec = flux.shape[0]
        for i in range(neml):
            if not self.quiet:
                print('Getting fit metrics for line: {0}/{1}'.format(i+1, neml), end='\r')
            start = model_eml_par['LINE_PIXC'][:,i] - metric_window//2
            end = start + metric_window
            m = numpy.zeros(flux.shape, dtype=float)
            for j in range(nspec):
                if mask[j,i]:
                    continue
                m[j,start[j]:end[j]] = 1.
                masked_pixel = flux_nc.mask[j,model_eml_par['LINE_PIXC'][j,i]] \
                                | ferr.mask[j,model_eml_par['LINE_PIXC'][j,i]]
                if masked_pixel:
                    model_eml_par['AMP'][j,i] = 0.0 
                    model_eml_par['ANR'][j,i] = 0.0
                else:     
                    model_eml_par['AMP'][j,i] = flux_nc[j,model_eml_par['LINE_PIXC'][j,i]]
                    model_eml_par['ANR'][j,i] = flux_nc[j,model_eml_par['LINE_PIXC'][j,i]] \
                                                    / ferr[j,model_eml_par['LINE_PIXC'][j,i]]
            m[spec_mask] = 0.

            model_eml_par['LINE_NSTAT'][:,i] = numpy.sum(m,axis=1)
            model_eml_par['LINE_RMS'][:,i] = numpy.ma.sqrt(numpy.ma.mean(m*resid,axis=1)
                                                           ).filled(0.0)
            model_eml_par['LINE_FRMS'][:,i] = numpy.ma.sqrt(numpy.ma.mean(m*fresid,axis=1)
                                                            ).filled(0.0)
            model_eml_par['LINE_CHI2'][:,i] = numpy.ma.sum(m*chisqr,axis=1).filled(0.0)

        if not self.quiet:
            print('Getting fit metrics for line: {0}/{0}'.format(neml))
#        print(model_eml_par['LINE_RMS'][nspec//2,:])
        return model_eml_par


    def _construct_2d_hdu(self, good_snr, model_flux, model_base, model_mask, model_eml_par,
                          model_fit_par=None, model_binid=None):
        """
        Construct :attr:`hdu` that is held in memory for manipulation of
        the object.  See :func:`construct_3d_hdu` if you want to convert
        the object into a DRP-like datacube.
        """
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, 'Constructing hdu ...')

        # Initialize the headers
        pri_hdr = self._initialize_primary_header()
        pri_hdr = self._add_method_header(pri_hdr, model_binid=model_binid)
        map_hdr = DAPFitsUtil.build_map_header(self.binned_spectra.drpf,
                                               'K Westfall <westfall@ucolick.org>')

        # Get the spatial map mask
        map_mask = numpy.zeros(self.spatial_shape, dtype=self.bitmask.minimum_dtype())

        # Allow the fitting function to redefine the bin IDs associated
        # with each model
        if model_binid is None:
            # Add any spaxel not used because it was flagged by the
            # binning step
            indx = self.binned_spectra['MAPMASK'].data > 0
            map_mask[indx] = self.bitmask.turn_on(map_mask[indx], 'DIDNOTUSE')
            # Isolate any spaxels with foreground stars
            indx = self.binned_spectra.bitmask.flagged(self.binned_spectra['MAPMASK'].data,
                                                       'FORESTAR')
            map_mask[indx] = self.bitmask.turn_on(map_mask[indx], 'FORESTAR')
            # Get the bins that were below the S/N limit
            indx = numpy.invert(DAPFitsUtil.reconstruct_map(self.spatial_shape,
                                                    self.binned_spectra['BINID'].data.ravel(),
                                                    good_snr, dtype='bool')) & (map_mask == 0)
            map_mask[indx] = self.bitmask.turn_on(map_mask[indx], 'LOW_SNR')

            # Get the bin ids with fitted models
            bin_indx = DAPFitsUtil.downselect_bins(self.binned_spectra['BINID'].data.ravel(),
                                            model_eml_par['BINID']).reshape(self.spatial_shape)
        else:
            # Assume any model with a binid less than zero is from a
            # spaxel that was not used
            indx = model_binid < 0
            map_mask[indx] = self.bitmask.turn_on(map_mask[indx], 'DIDNOTUSE')

            # The number of valid bins MUST match the number of
            # measurements
            nvalid = numpy.sum(numpy.invert(indx))
            if nvalid != len(model_eml_par):
                raise ValueError('Provided model id does not match the number of measurements.')

            # Get the bin ids with fitted models
            bin_indx = model_binid

        # Allow the fitting function to return a boolean model mask
        if model_mask.dtype is numpy.dtype('bool'):
            _model_mask = numpy.zeros(model_flux.shape, dtype=self.bitmask.minimum_dtype())
            _model_mask[model_mask] = self.bitmask.turn_on(_model_mask[model_mask], 'DIDNOTUSE')
        else:
            _model_mask = model_mask

        # Compile the information on the suite of measured indices
        line_database = self._compile_database()

        # Save the data to the hdu attribute
        par_hdu = fits.BinTableHDU(data=None, name='PAR') if line_database is None else \
                        fits.BinTableHDU.from_columns([ fits.Column(name=n,
                                                        format=rec_to_fits_type(line_database[n]),
                                                        array=line_database[n])
                                                            for n in line_database.dtype.names],
                                                      name='PAR')

        fit_hdu = fits.BinTableHDU(data=None, name='FIT') if model_fit_par is None else \
                        fits.BinTableHDU.from_columns([fits.Column(name=n,
                                                        format=rec_to_fits_type(model_fit_par[n]),
                                                        dim=rec_to_fits_col_dim(model_fit_par[n]),
                                                        array=model_fit_par[n])
                                                            for n in model_fit_par.dtype.names],
                                                      name='FIT')

        self.hdu = fits.HDUList([ fits.PrimaryHDU(header=pri_hdr),
                                  fits.ImageHDU(data=model_flux.data, name='FLUX'),
                                  fits.ImageHDU(data=model_base.data, name='BASE'),
                                  fits.ImageHDU(data=_model_mask, name='MASK'),
                                  self.binned_spectra['WAVE'].copy(),
                                  fits.ImageHDU(data=bin_indx, header=map_hdr, name='BINID'),
                                  fits.ImageHDU(data=map_mask, header=map_hdr, name='MAPMASK'),
                                  par_hdu, fit_hdu,
                                  fits.BinTableHDU.from_columns([fits.Column(name=n,
                                                        format=rec_to_fits_type(model_eml_par[n]),
                                                        dim=rec_to_fits_col_dim(model_eml_par[n]),
                                        array=model_eml_par[n]) for n in model_eml_par.dtype.names],
                                                                    name='EMLDATA')
                                ])


    def file_name(self):
        """Return the name of the output file."""
        return self.output_file


    def file_path(self):
        """Return the full path to the output file."""
        if self.directory_path is None or self.output_file is None:
            return None
        return os.path.join(self.directory_path, self.output_file)


    def fit(self, binned_spectra, stellar_continuum=None, emission_line_moments=None,
            redshift=None, dispersion=None, minimum_error=None, dapsrc=None, dapver=None,
            analysis_path=None, directory_path=None, output_file=None, hardcopy=True,
            clobber=False, loggers=None, quiet=False):
        """
        Fit the emission lines.
        """
        # Initialize the reporting
        if loggers is not None:
            self.loggers = loggers
        self.quiet = quiet

        # SpatiallyBinnedSpectra object always needed
        if binned_spectra is None:
            raise ValueError('Must provide spectra object for fitting.')
        if not isinstance(binned_spectra, SpatiallyBinnedSpectra):
            raise TypeError('Must provide a valid SpatiallyBinnedSpectra object!')
        if binned_spectra.hdu is None:
            raise ValueError('Provided SpatiallyBinnedSpectra object is undefined!')
        self.binned_spectra = binned_spectra

        # Continuum accounts for underlying absorption
        if stellar_continuum is not None \
                and not isinstance(stellar_continuum, StellarContinuumModel):
            raise TypeError('Must provide a valid StellarContinuumModel object.')
        self.stellar_continuum = stellar_continuum

        # Emission-line moments can be used to set the initial guess
        # kinematics
        if emission_line_moments is not None:
            # No longer checks for this; errors will be raised
            # just by trying to use the object...
#            if not isinstance(emission_line_moments, EmissionLineMoments):
#                raise TypeError('Must provide a valid EmissionLineMoments object.')
            names = emission_line_moments.channel_names()
            if self.method['mom_vel_name'] is not None \
                        and self.method['mom_vel_name'] not in names:
                raise ValueError('{0}'.format(self.method['mom_vel_name']) +
                                 ' is not a valid name in the EmissionLineMoments object.')
            if self.method['mom_disp_name'] is not None \
                        and self.method['mom_disp_name'] not in names:
                raise ValueError('{0}'.format(self.method['mom_disp_name']) +
                                 ' is not a valid name in the EmissionLineMoments object.')

        self.shape = self.binned_spectra.shape
        self.spatial_shape =self.binned_spectra.spatial_shape
#       nspec is never used and it won't be correct if the bins are
#       deconstructed.
#        self.nspec = self.binned_spectra.nspec
        self.spatial_index = self.binned_spectra.spatial_index.copy()
        self.dispaxis = self.binned_spectra.dispaxis
        self.nwave = self.binned_spectra.nwave
        
        # Get the guess kinematics
        self._assign_input_kinematics(emission_line_moments, redshift, dispersion)

        #---------------------------------------------------------------
        # Get the good spectra
        # TODO: Put this in a function that can be used by Sasuke...
        good_snr = EmissionLineFit.select_binned_spectra_to_fit(self.binned_spectra,
                                                        minimum_snr=self.method['minimum_snr'],
                                                        stellar_continuum=self.stellar_continuum)

        # Report
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, '-'*50)
            log_output(loggers, 1, logging.INFO, '{0:^50}'.format('EMISSION-LINE PROFILE FITTING'))
            log_output(self.loggers, 1, logging.INFO, '-'*50)
            log_output(self.loggers, 1, logging.INFO, 'Number of binned spectra: {0}'.format(
                                                            self.binned_spectra.nbins))
            if len(self.binned_spectra.missing_bins) > 0:
                log_output(self.loggers, 1, logging.INFO, 'Missing bins: {0}'.format(
                                                            len(self.binned_spectra.missing_bins)))
            log_output(self.loggers, 1, logging.INFO, 'With good S/N and to fit: {0}'.format(
                                                            numpy.sum(good_snr)))
            log_output(self.loggers, 1, logging.INFO, 'Bins deconstructed in fitting: {0}'.format(
                                                      str(self.method['deconstruct_bins'])))
            
        if numpy.sum(good_snr) == 0:
            raise ValueError('No good spectra to fit!')

        #---------------------------------------------------------------
        # Fill in any remaining binning parameters
        self._fill_method_par(dapsrc=dapsrc, analysis_path=analysis_path)

        # (Re)Set the output paths
        self._set_paths(directory_path, dapver, analysis_path, output_file)

        #---------------------------------------------------------------
        # Check that the file path is defined
        ofile = self.file_path()
        if ofile is None:
            raise ValueError('File path for output file is undefined!')

        # Report
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO,
                       'Output path: {0}'.format(self.directory_path))
            log_output(self.loggers, 1, logging.INFO,
                       'Output file: {0}'.format(self.output_file))

        #---------------------------------------------------------------
        # If the file already exists, and not clobbering, just read the
        # file
        if os.path.isfile(ofile) and not clobber:
            self.hardcopy = True
            if not self.quiet:
                log_output(self.loggers, 1, logging.INFO, 'Using existing file')
            self.read(checksum=self.checksum)
            if not self.quiet:
                log_output(self.loggers, 1, logging.INFO, '-'*50)
            return

        #---------------------------------------------------------------
        # Fit the spectra
        # Mask should be fully defined within the fitting function
        model_flux, model_base, model_mask, model_fit_par, model_eml_par, model_binid = \
                self.method['fitfunc'](self.binned_spectra, par=self.method['fitpar'],
                                       loggers=self.loggers, quiet=self.quiet)

        # Impose a minimum error because of the conversion to inverse
        # variance when constructing a MAPS file; applied to the error
        # quantities from spectralfitting.EmissionLineFit dtype
        # parameters: FLUXERR, KINERR, EWERR
        # TODO: Add a maskbit as well?  Yes, but in MAPS construction
        indx = model_eml_par['FLUXERR'] < (numpy.finfo(model_eml_par['FLUX'].dtype).eps \
                                            if minimum_error is None else minimum_error)
        model_eml_par['FLUXERR'][indx] = 0.0
        indx = model_eml_par['KINERR'] < (numpy.finfo(model_eml_par['KIN'].dtype).eps \
                                            if minimum_error is None else minimum_error)
        model_eml_par['KINERR'][indx] = 0.0
        indx = model_eml_par['EWERR'] < (numpy.finfo(model_eml_par['EW'].dtype).eps \
                                            if minimum_error is None else minimum_error)
        model_eml_par['EWERR'][indx] = 0.0

        # TODO: Test if, when deconstruct_bins is true, the model_binid
        # spaxels (with fits) all have unique bin IDs.

        # The emission-line fitting can proceed without knowing which
        # lines are fit; this sets the number of lines
        if self.neml is None:
            self.neml = model_eml_par['FLUX'].shape[1]

#        pyplot.step(model_wave, model_flux[0,:], where='mid')
#        pyplot.show()
       
        # The number of models returned should be the number of "good" binned
        # spectra if the fitter does *not* deconstruct the bins.
        # TODO: This may now cause problems with Elric...
        if not self.method['deconstruct_bins'] and model_flux.shape[0] != numpy.sum(good_snr):
            print(model_flux.shape[0])
            print(numpy.sum(good_snr))
            raise ValueError('Unexpected returned shape of fitted emission-line models.')
       
        #---------------------------------------------------------------
        # Set the number of models and the missing models
        self.nmodels = model_flux.shape[0]
        unique_bins = None if model_binid is None else numpy.unique(model_binid.ravel())
        self.missing_models = self._get_missing_models(unique_bins=unique_bins)

        # Report
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, 'Fitted models: {0}'.format(self.nmodels))
            if len(self.missing_models) > 0:
                log_output(self.loggers, 1, logging.INFO, 'Missing models: {0}'.format(
                                                            len(self.missing_models)))
            log_output(self.loggers, 1, logging.INFO, 'Calculating fit metrics')

        # Compute the line-fit metrics
        model_eml_par = self._get_line_fit_metrics(model_flux, model_base, model_mask,
                                                   model_eml_par, model_binid)

        # Construct the 2d hdu list that only contains the fitted models
        self.hardcopy = hardcopy
        self._construct_2d_hdu(good_snr, model_flux, model_base, model_mask, model_eml_par,
                               model_fit_par=model_fit_par, model_binid=model_binid)

        #---------------------------------------------------------------
        # Write the data, if requested
        if self.hardcopy:
            if not os.path.isdir(self.directory_path):
                os.makedirs(self.directory_path)
            self.write(clobber=clobber)
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, '-'*50)


    def construct_3d_hdu(self):
        """
        Reformat the model spectra into a cube matching the shape of
        the DRP fits file.
        """
        # Report
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO,
                       'Constructing emission-line model datacube ...')

        bin_indx = self.hdu['BINID'].data.copy()
        model_flux = self.hdu['FLUX'].data.copy()
        model_base = self.hdu['BASE'].data.copy()
        model_mask = self.hdu['MASK'].data.copy()

        flux, base, mask = DAPFitsUtil.reconstruct_cube(self.shape, bin_indx.ravel(),
                                                        [model_flux, model_base, model_mask])

        mask = self._finalize_cube_mask(mask)

        # Primary header is identical regardless of the shape of the
        # extensions
        hdr = self.hdu['PRIMARY'].header.copy()
        cube_hdr = DAPFitsUtil.build_cube_header(self.binned_spectra.drpf,
                                                 'K Westfall <westfall@ucolick.org>')

        par_hdu = self.hdu['PAR'].copy() if 'data' in self.hdu['PAR'].__dict__ \
                        else fits.BinTableHDU(data=None, name='PAR')
        fit_hdu = self.hdu['FIT'].copy() if 'data' in self.hdu['FIT'].__dict__ \
                        else fits.BinTableHDU(data=None, name='FIT')

        # Return the converted hdu without altering the internal hdu
        return fits.HDUList([ fits.PrimaryHDU(header=hdr),
                              fits.ImageHDU(data=flux, header=cube_hdr, name='FLUX'),
                              fits.ImageHDU(data=base, header=cube_hdr, name='BASE'),
                              fits.ImageHDU(data=mask, header=cube_hdr, name='MASK'),
                              self.hdu['WAVE'].copy(),
                              self.hdu['BINID'].copy(),
                              par_hdu,
                              fit_hdu,
                              self.hdu['EMLDATA'].copy()
                            ])


    # Exact same function as used by SpatiallyBinnedSpectra
    def write(self, match_DRP=False, clobber=False):
        """
        Write the hdu object to the file.
        """
        # Convert the spectral arrays in the HDU to a 3D cube and write
        # it
        if match_DRP:
            hdu = self.construct_3d_hdu()
            DAPFitsUtil.write(hdu, self.file_path(), clobber=clobber, checksum=True,
                              loggers=self.loggers, quiet=self.quiet)
            return
        # Just write the unique (2D) data
        DAPFitsUtil.write(self.hdu, self.file_path(), clobber=clobber, checksum=True,
                          loggers=self.loggers, quiet=self.quiet) 


    def read(self, ifile=None, strict=True, checksum=False):
        """
        Read an existing file with an existing set of emission-line
        models.
        """
        if ifile is None:
            ifile = self.file_path()
        if not os.path.isfile(ifile):
            raise FileNotFoundError('File does not exist!: {0}'.format(ifile))

        if self.hdu is not None:
            self.hdu.close()

#        self.hdu = fits.open(ifile, checksum=checksum)
        self.hdu = DAPFitsUtil.read(ifile, checksum=checksum)

        # Force load the PAR and FIT extensions because it may be
        # checked by self.construct_3d_hdu() or construct_maps_file in
        # dapfits.py
        try:
            if self.hdu['PAR'].data is None:
                pass
        except:
            warnings.warn('Could not load emission-line parameter data.')
        try:
            if self.hdu['FIT'].data is None:
                pass
        except:
            warnings.warn('Could not load emission-line fit data.')

        # Make sure that the number of emission-lines is set
        _neml = self.hdu['EMLDATA'].data['FLUX'].shape[1]
        if self.neml is not None and _neml != self.neml:
            raise ValueError('Num. of emission lines read does not match emission-line database.')
        if self.neml is None:
            self.neml = _neml

        # Confirm that the internal method is the same as the method
        # that was used in writing the file
        if self.hdu['PRIMARY'].header['ELFKEY'] != self.method['key']:
            if strict:
                raise ValueError('Keywords in header do not match specified method keywords!')
            else:
                warnings.warn('Keywords in header do not match specified method keywords!')
        # TODO: "strict" should also check other aspects of the file to
        # make sure that the details of the method are also the same,
        # not just the keyword

#        if not self.quiet:
#            log_output(self.loggers, 1, logging.INFO, 'Reverting to python-native structure.')
#        DAPFitsUtil.restructure_rss(self.hdu, ext=self.spectral_arrays)
#        DAPFitsUtil.restructure_map(self.hdu, ext=self.image_arrays)

        # Attempt to read the modeling parameters
        if self.method['fitpar'] is not None:
            try:
                self.method['fitpar'].fromheader(self.hdu['PRIMARY'].header)
            except:
                if not self.quiet:
                    warnings.warn('Fit parameter object does not have a fromheader ' \
                                  'attribute.  No parameters head from header.')

        self.nmodels = self.hdu['FLUX'].shape[0]
        unique_bins = numpy.unique(self.hdu['BINID'].data.ravel()) \
                            if self.hdu['PRIMARY'].header['ELREBIN'] else None
        self.missing_models = self._get_missing_models(unique_bins=unique_bins)


    def copy_to_array(self, ext='FLUX', waverange=None, include_missing=False):
        r"""

        Wrapper for :func:`mangadap.util.fitsutil.DAPFitsUtil.copy_to_array`
        specific for :class:`EmissionLineModel`.

        Return a copy of the selected data array.  The array size is
        always :math:`N_{\rm models} \times N_{\rm wavelength}`; i.e., the
        data is always flattened to two dimensions and the unique
        spectra are selected.

        Args:
            ext (str) : (**Optional**) Name of the extension from which
                to draw the data.  Must be allowed for the current
                :attr:`mode`; see :attr:`data_arrays`.  Default is
                ``'FLUX'``.
            waverange (array-like) : (**Optional**) Two-element array
                with the first and last wavelength to include in the
                computation.  Default is to use the full wavelength
                range.
            include_missing (bool) : (**Optional**) Create an array with
                a size that accommodates the missing models.

        Returns:
            numpy.ndarray: A 2D array with a copy of the data from the
            selected extension.

        """
        return DAPFitsUtil.copy_to_array(self.hdu, ext=ext, allowed_ext=self.spectral_arrays,
                                         waverange=waverange,
                                    missing_bins=self.missing_models if include_missing else None,
                                         nbins=self.nmodels,
                                        unique_bins=DAPFitsUtil.unique_bins(self.hdu['BINID'].data))


    def copy_to_masked_array(self, ext='FLUX', flag=None, waverange=None, include_missing=False):
        r"""
        Wrapper for
        :func:`mangadap.util.fitsutil.DAPFitsUtil.copy_to_masked_array`
        specific for :class:`EmissionLineModel`.

        Return a copy of the selected data array as a masked array.
        This is functionally identical to :func:`copy_to_array`,
        except the output format is a `numpy.ma.MaskedArray`_.  The
        pixels that are considered to be masked can be specified using
        `flag`.
        
        Args:
            ext (str) : (**Optional**) Name of the extension from which
                to draw the data.  Must be allowed for the current
                :attr:`mode`; see :attr:`data_arrays`.  Default is
                `'FLUX'`.
            flag (str or list): (**Optional**) (List of) Flag names that
                are considered when deciding if a pixel should be
                masked.  The names *must* be a valid bit name as defined
                by :attr:`bitmask`.  If not provided, *ANY* non-zero
                mask bit is omitted.
            waverange (array-like) : (**Optional**) Two-element array
                with the first and last wavelength to include in the
                computation.  Default is to use the full wavelength
                range.
            include_missing (bool) : (**Optional**) Create an array with
                a size that accommodates the missing models.

        Returns:
            numpy.ndarray: A 2D array with a copy of the data from the
            selected extension.
        """
        return DAPFitsUtil.copy_to_masked_array(self.hdu, ext=ext, mask_ext='MASK', flag=flag,
                                                bitmask=self.bitmask,
                                                allowed_ext=self.spectral_arrays,
                                                waverange=waverange,
                                    missing_bins=self.missing_models if include_missing else None,
                                                nbins=self.nmodels,
                                        unique_bins=DAPFitsUtil.unique_bins(self.hdu['BINID'].data))


#    def fill_to_match(self, binned_spectra, include_base=False):
#        """
#        Get the emission-line models that match the shape of the
#        provided binned_spectra.
#        """
#        if binned_spectra is self.binned_spectra:
#            emission_lines = self.copy_to_array()
#            if include_base:
#                emission_lines += self.copy_to_array(ext='BASE')
#            emission_lines = numpy.ma.MaskedArray(emission_lines)
#
#            # Number of models matches the numbers of bins
#            if binned_spectra.nbins == self.nmodels:
#                return emission_lines
#    
#            # Fill in bins with no models with masked zeros
#            _emission_lines = numpy.ma.zeros(binned_spectra['FLUX'].data.shape, dtype=float)
#            _emission_lines[:,:] = numpy.ma.masked
#            for i,j in enumerate(self.hdu['EMLDATA'].data['BINID_INDEX']):
#                _emission_lines[j,:] = emission_lines[i,:]
#            return _emission_lines
#
#        raise NotImplementedError('Can only match to internal binned_spectra.')


    def fill_to_match(self, binid, include_base=False, missing=None):
        """
        Get the emission-line model that matches the input bin ID
        matrix.  The output is a 2D matrix ordered by the bin ID; any
        skipped index numbers in the maximum of the union of the unique
        numbers in the `binid` and `missing` input are masked.

        Use `include_base` to include the baseline in the output
        emission-line models.

        """
        # The input bin id array must have the same shape as self
        if binid.shape != self.spatial_shape:
            raise ValueError('Input bin ID matrix has incorrect shape.')

        # Construct the best-fitting models
        best_fit_model = self.copy_to_array()
        if include_base:
            best_fit_model += self.copy_to_array(ext='BASE')

        # Get the number of output models
        nbins = numpy.amax(binid).astype(int)+1 if missing is None or len(missing) == 0 else \
                    numpy.amax( numpy.append(binid.ravel(), missing) ).astype(int)+1

        # Map the BINID to the spectrum index, assuming bins are sorted
        # and that the BINID map has -1 BINID values
        u, indx, reconstruct = numpy.unique(self['BINID'].data.ravel(), return_index=True,
                                            return_inverse=True)
        u_bin_indx = numpy.arange(len(u))-1
        _bin_indx = u_bin_indx[reconstruct].reshape(self.spatial_shape)

        # Fill in bins with no models with masked zeros
        emission_lines = numpy.ma.zeros((nbins,best_fit_model.shape[1]), dtype=float)
        emission_lines[:,:] = numpy.ma.masked
        for i,j in zip(binid.ravel(), _bin_indx.ravel()):
            if i < 0 or j < 0:
                continue
            emission_lines[i,:] = best_fit_model[j,:]

        return emission_lines


    def matched_kinematics(self, binid, line_name, redshift=None, dispersion=100.0, constant=False,
                           cz=False, corrected=False, min_dispersion=None, nearest=False,
                           missing=None):
        """
        Return the redshift and velocity dispersion for all the
        spectra based on the emission-line model fit to a single line.

        For spectra the were not fit, use the median (unmasked)
        kinematic measurement if no default value is provided as an
        argument.

        Args:
            binid (numpy.ndarray): 2D array with the bin ID associate
                with each spaxel in the datacube.
            line (str): Line to use for the kinematics.  Must match one
                of the lines created by :func:`channel_names`.  Function
                will raise an exception if the line_name is None or if
                the channel names cannot be constructed.
            redshift (float): (**Optional**) The default redshift to use
                for spectra without an emission-line fit.  Default is to
                use the median of the unmasked measurements.
            dispersion (float): (**Optional**) The default velocity
                dispersion to use for spectra without an emission-line
                fit.  Default is 100 km/s.
            constant (bool): (**Optional**) Force the function to return
                a constant redshift and dispersion for each spectrum,
                regardless of any fitted kinematics.
            cz (bool): (**Optional**) Return the redshift as cz in km/s,
                as opposed to unitless z.
            corrected (bool): (**Optional**) By default, the returned
                velocity dispersions are the measured values, which will
                include the instrumental resolution.  Setting corrected
                to True will return the corrected dispersions.
            min_dispersion (float): (**Optional**) Impose a minimum
                dispersion.
            nearest (bool): (**Optional**) Instead of the median of the
                results for the spectra that were not fit, use the value
                from the nearest bin.
            missing (list): (**Optional**) Any bin ID numbers missing
                from the input `binid` image needed for constructing the
                output matched data.

        Returns:
            numpy.ndarray: Returns arrays with a redshift (or cz) and
            dispersion for each spectrum.

        Raises:
            TypeError: Raised if the input redshift or dispersion values
                are not single numbers.
        """
        # Check input
        if redshift is not None and not isinstance(redshift, (float,int)):
            raise TypeError('Redshift must be a single number or None.')
        if dispersion is not None and not isinstance(dispersion, (float,int)):
            raise TypeError('Dispersion must be a single number or None.')
        if binid.shape != self.spatial_shape:
            raise ValueError('Input bin ID map must match the spatial shape of the DRP cube.')

        # Get the channel names
        names = self.channel_names()
        if line_name not in names:
            raise ValueError('Line {0} not present in channel list.'.format(line_name))
        vel_channel = numpy.arange(len(names))[names == line_name][0]

        # Get the number of output kinematics
        nbins = numpy.amax(binid).astype(int)+1 if missing is None or len(missing) == 0 else \
                    numpy.amax( numpy.append(binid.ravel(), missing) ).astype(int)+1

        # Mask bad values
        mask = self.bitmask.flagged(self.hdu['EMLDATA'].data['MASK'][:,vel_channel].copy(),
                                    [ 'NO_FIT', 'FIT_FAILED', 'INSUFFICIENT_DATA', 'NEAR_BOUND' ])

        if numpy.all(mask):
            # Unlike StellarContinuumModel.matched_kinematics, the input
            # guesses here always are based on the number of binned
            # spectra.  So if the class is meant to deconstruct these
            # bins, we can't just use the guess kinematics when all the
            # output is masked.  So we'll just use the defaults.
            warnings.warn('All emission line models have been masked!  Using defaults.')
            eml_z = numpy.ma.MaskedArray(numpy.full(self.nmodels,
                                                    (0.0 if redshift is None else redshift),
                                                    dtype=float))
            eml_d = numpy.ma.MaskedArray(numpy.full(self.nmodels,
                                    (100.0 if dispersion is None else dispersion), dtype=float))
        else:
            # Pull out the data
            eml_z = numpy.ma.MaskedArray(self.hdu['EMLDATA'].data['KIN'][:,vel_channel,0].copy(),
                                            mask=mask) / astropy.constants.c.to('km/s').value
            eml_d = numpy.ma.MaskedArray(self.hdu['EMLDATA'].data['KIN'][:,vel_channel,1].copy(),
                                            mask=mask)

            # Apply the sigma corrections if requested.  Values with the
            # correction larger than the measured value will be masked!
            if corrected:
                sigma_corr = numpy.ma.MaskedArray(
                                self.hdu['EMLDATA'].data['SIGMACORR'][:,vel_channel], mask=mask)
                eml_d = numpy.ma.sqrt(numpy.square(eml_d) - numpy.square(sigma_corr))

        # Set the default values to use when necessary
        _redshift = numpy.ma.median(eml_z) if redshift is None else redshift
        _dispersion = numpy.ma.median(eml_d) if dispersion is None else dispersion
        if min_dispersion is not None and _dispersion < min_dispersion:
            _dispersion = min_dispersion

        # Just return the single value
        if constant:
            eml_z = numpy.full(nbins, _redshift * astropy.constants.c.to('km/s').value
                                        if cz else _redshift, dtype=float)
            eml_d = numpy.full(nbins, _dispersion, dtype=float)
            return eml_z, eml_d

        if nearest:
            # Fill masked values with the nearest datum
            best_fit_kinematics = numpy.ma.append([eml_z], [eml_d], axis=0).T
            if self.method['deconstruct_bins']:
                indx = self['BINID'].data.ravel() > -1
                coo = self.binned_spectra.rdxqa['SPECTRUM'].data['SKY_COO'][indx,:]
            else:
                valid_bins = numpy.unique(self['BINID'].data)[1:]
                coo = self.binned_spectra['BINS'].data['SKY_COO'][valid_bins,:]
            replace = eml_z.mask | eml_d.mask
            kinematics = replace_with_data_from_nearest_coo(coo, best_fit_kinematics, replace)
            eml_z = kinematics[:,0]
            eml_d = kinematics[:,1]
        else:
            # Fill any masked values with the single estimate
            eml_z = eml_z.filled(_redshift)
            eml_d = eml_d.filled(_dispersion)

        # Map the BINID to the spectrum index, assuming bins are sorted
        # and that the BINID map has -1 BINID values
        u, indx, reconstruct = numpy.unique(self['BINID'].data.ravel(), return_index=True,
                                            return_inverse=True)
        u_bin_indx = numpy.arange(len(u))-1
        _bin_indx = u_bin_indx[reconstruct].reshape(self.spatial_shape)

        # Match the kinematics to the output bin ID map
        _eml_z = numpy.ma.masked_all(nbins, dtype=float)
        _eml_d = numpy.ma.masked_all(nbins, dtype=float)
        for i,j in zip(binid.ravel(), _bin_indx.ravel()):
            if i < 0 or j < 0:
                continue
            _eml_z[i] = eml_z[j]
            _eml_d[i] = eml_d[j]

        eml_z = _eml_z.filled(_redshift)
        eml_d = _eml_d.filled(_dispersion)
        # Convert to cz velocities (km/s)
        if cz:
            eml_z *= astropy.constants.c.to('km/s').value
        return eml_z, eml_d


    def fill_continuum_to_match(self, binid, replacement_templates=None, redshift_only=False,
                                deredshift=False, corrected_dispersion=False, missing=None):
        """
        Get the emission-line continuum model, if possible, that matches
        the input bin ID matrix.  The output is a 2D matrix ordered by
        the bin ID; any skipped index numbers in the maximum of the
        union of the unique numbers in the `binid` and `missing` input
        are masked.

        Use `replacement_templates` only if the number of templates is
        identical to the number used during the fitting procedure.  Use
        `redshift_only` to produce the best-fitting model without and
        velocity dispersion.  Use `corrected_dispersion` to produce the
        model using the corrected velocity dispersion.

        """
        if self.stellar_continuum is None:
            raise ValueError('EmissionLineModel object has no associated StellarContinuumModel.')
        if self.method['fitclass'] is None:
            raise ValueError('No class object available for constructing the model!')
        if not callable(self.method['fitclass'].construct_continuum_models):
            raise AttributeError('Provided fit class object has no callable ' \
                                 '\'construct_continuum_models\' attribute!')
        if corrected_dispersion and 'NEAREST_BIN' not in self.hdu['FIT'].columns.names:
            raise KeyError('Nearest stellar continuum bin not set in fit parameters.  Cannot '
                           'apply dispersion corrections.')

        # The input bin id array must have the same shape as self
        if binid.shape != self.spatial_shape:
            raise ValueError('Input bin ID matrix has incorrect shape.')

        # Pull out the dispersion corrections, if requested
        dispersion_corrections = None
        if corrected_dispersion:
            nsc = self.stellar_continuum['PAR'].data['BINID'].size
            indx = [ numpy.arange(nsc)[self.stellar_continuum['PAR'].data['BINID'] == nb][0]
                        for nb in self.hdu['FIT'].data['NEAREST_BIN'] ]
            dispersion_corrections = self.stellar_continuum['PAR'].data['SIGMACORR_EMP'][indx]

        # Get the continuum models
        best_fit_continuum = self.construct_continuum_models(
                                            replacement_templates=replacement_templates,
                                            deredshift=deredshift, redshift_only=redshift_only,
                                            dispersion_corrections=dispersion_corrections)

        # Get the number of output continuua
        nbins = numpy.amax(binid).astype(int)+1 if missing is None or len(missing) == 0 else \
                    numpy.amax( numpy.append(binid.ravel(), missing) ).astype(int)+1

        # Map the BINID to the spectrum index, assuming bins are sorted
        # and that the BINID map has -1 BINID values
        u, indx, reconstruct = numpy.unique(self['BINID'].data.ravel(), return_index=True,
                                            return_inverse=True)
        u_bin_indx = numpy.arange(len(u))-1
        _bin_indx = u_bin_indx[reconstruct].reshape(self.spatial_shape)

        # Fill in bins with no models with masked zeros
        continuum = numpy.ma.zeros((nbins,best_fit_continuum.shape[1]), dtype=float)
        continuum[:,:] = numpy.ma.masked
        for i,j in zip(binid.ravel(), _bin_indx.ravel()):
            if i < 0 or j < 0:
                continue
            continuum[i,:] = best_fit_continuum[j,:]

        return continuum


    def construct_continuum_models(self, replacement_templates=None, redshift_only=False,
                                   deredshift=False, dispersion_corrections=None):

        if self.stellar_continuum is None:
            raise ValueError('EmissionLineModel object has no associated StellarContinuumModel.')

        if self.method['fitclass'] is None:
            raise ValueError('No class object available for constructing the model!')
        if not callable(self.method['fitclass'].construct_continuum_models):
            raise AttributeError('Provided fit class object has no callable ' \
                                 '\'construct_continuum_models\' attribute!')
        if replacement_templates is not None \
                and not isinstance(replacement_templates, TemplateLibrary):
            raise TypeError('Provided template replacements must have type TemplateLibrary.')

        # Only the selected models are constructed, others are masked
        select = numpy.invert(self.bitmask.flagged(self.hdu['FIT'].data['MASK'],
                                            flag=['NO_FIT', 'FIT_FAILED', 'INSUFFICIENT_DATA']))
        templates = self.stellar_continuum.method['fitpar']['template_library'] \
                            if replacement_templates is None else replacement_templates

        return self.method['fitclass'].construct_continuum_models(
                                        self.emldb, templates['WAVE'].data, templates['FLUX'].data,
                                        self.hdu['WAVE'].data, self.hdu['FLUX'].data.shape,
                                        self.hdu['FIT'].data, select=select,
                                        redshift_only=redshift_only, deredshift=deredshift,
                                        dispersion_corrections=dispersion_corrections, dvtol=1e-9)

    # TODO: Use the EmissionMomentsDB.channel_names function!
    def channel_names(self):
        line_database = self._compile_database(quiet=True)
        if line_database is None:
            raise ValueError('Channel names undefined because line database is not available.')
        
        return numpy.array([ '{0}-{1}'.format(n,int(w)) \
                        for n,w in zip(self['PAR'].data['NAME'],
                                       self['PAR'].data['RESTWAVE']) ])


