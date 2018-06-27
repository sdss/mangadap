# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
A class hierarchy that performs the stellar-continuum fitting.

*License*:
    Copyright (c) 2015, SDSS-IV/MaNGA Pipeline Group
        Licensed under BSD 3-clause license - see LICENSE.rst

*Source location*:
    $MANGADAP_DIR/python/mangadap/proc/stellarcontinuummodel.py

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

        import os
        import glob
        import logging
        import numpy

        from astropy.io import fits
        import astropy.constants

        from ..drpfits import DRPFits
        from ..par.parset import ParSet
        from ..par.artifactdb import ArtifactDB
        from ..par.emissionlinedb import EmissionLineDB
        from ..util.log import log_output
        from ..util.fitsutil import DAPFitsUtil
        from ..util.fileio import rec_to_fits_type
        from ..util.instrument import spectrum_velocity_scale
        from ..util.bitmask import BitMask
        from ..util.pixelmask import SpectralPixelMask
        from ..util.parser import DefaultConfig
        from ..config.defaults import dap_source_dir, default_dap_file_name
        from ..config.defaults import default_dap_method, default_dap_method_path
        from .spatiallybinnedspectra import SpatiallyBinnedSpectra
        from .templatelibrary import TemplateLibrary
        from .ppxffit import PPXFFitPar, PPXFFit
        from .util import select_proc_method

*Class usage examples*:
        Add examples

*Revision history*:
    | **14 Apr 2016**: Implementation begun by K. Westfall (KBW)
    | **19 Apr 2016**: (KBW) First version
    | **19 May 2016**: (KBW) Added loggers and quiet keyword arguments
        to :class:`StellarContinuumModel`, removed verbose 
    | **05 Jul 2016**: (KBW) Removed oversample keyword from instantion
        of :class:`mangadap.proc.ppxffit.PPXFFit` objects.
    | **08 Nov 2016**: (KBW) Moved
        :func:`StellarContinuumModel.reset_continuum_mask_window` from
        :class:`Elric` to here.  The function allows one to deal with
        the subtraction of the continuum over the fully viable spectral
        range, ignoring small spectral regions that were ignored during
        the stellar continuum fit.  Also added wrapper function,
        :func:`StellarContinuumModel.emission_line_continuum_model`.
    | **11 Jan 2017**: (KBW) Changed
        StellarContinuumModel.emission_line_continuum_model to
        :func:`StellarContinuumModel.unmasked_continuum_model`.
    | **23 Feb 2017**: (KBW) Use DAPFitsUtil read and write functions.

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

import os
import glob
import logging
import numpy

from astropy.io import fits
import astropy.constants

from ..drpfits import DRPFits
from ..par.parset import ParSet
from ..par.artifactdb import ArtifactDB
from ..par.emissionlinedb import EmissionLineDB
from ..util.log import log_output
from ..util.fitsutil import DAPFitsUtil
from ..util.fileio import rec_to_fits_type
from ..util.instrument import spectral_coordinate_step, spectrum_velocity_scale, SpectralResolution
from ..util.bitmask import BitMask
from ..util.pixelmask import SpectralPixelMask
from ..util.parser import DefaultConfig
from ..config.defaults import dap_source_dir, default_dap_file_name
from ..config.defaults import default_dap_method, default_dap_method_path
from .spatiallybinnedspectra import SpatiallyBinnedSpectra
from .templatelibrary import TemplateLibrary
from .ppxffit import PPXFFitPar, PPXFFit
from .util import select_proc_method, replace_with_data_from_nearest_coo

from matplotlib import pyplot
#from memory_profiler import profile

# Add strict versioning
# from distutils.version import StrictVersion


class StellarContinuumModelDef(ParSet):
    """
    A class that holds the parameters necessary to perform the
    stellar-continuum fitting.

    The function use to fit the stellar continuum model must have the
    form::

        model_wave, model_flux, model_mask, model_par \
                = fitfunc(self, binned_spectra, par=None)

    where `binned_spectra` has type
    :class:`mangadap.proc.spatiallybinnedspetra.SpatiallyBinnedSpectra`,
    and the returned objects are the wavelength vector, the fitted model
    flux, a bitmask for the model, and the set of model parameters.  For
    example, see
    :func:`mangadap.proc.ppxffit.PPXFFit.fit_SpatiallyBinnedSpectra`.

    Args:
        key (str): Keyword used to distinguish between different spatial
            binning schemes.
        minimum_snr (bool): Minimum S/N of spectrum to fit
        fit_type (str): Currently can be anything.  In the DAP, this is
            used to identify the primary purpose of the fit as either
            producing stellar kinematics, composition measurements, or
            emission line fits; see
            :mod:`mangadap.proc.spectralfitting'.  The purpose for this
            type is to isolate the expected format of the binary table
            data; see, e.g.,
            :func:`mangadap.proc.spectralfitting._per_stellar_kinematics_dtype`.
        fitpar (:class:`mangadap.par.parset.ParSet` or dict): Any
            additional parameters, aside from the spectra themselves,
            required by the fitting function.
        fitclass (object): Instance of class object to use for the
            model fitting.  Needed in case fitfunc is a non-static member
            function of a class.
        fitfunc (callable): The function that models the spectra
    """
    def __init__(self, key, minimum_snr, fit_type, fitpar, fitclass, fitfunc):
        in_fl = [ int, float ]
        par_opt = [ ParSet, dict ]

        pars =     [ 'key', 'minimum_snr', 'fit_type', 'fitpar', 'fitclass', 'fitfunc' ]
        values =   [   key,   minimum_snr,   fit_type,   fitpar,   fitclass,   fitfunc ]
        dtypes =   [   str,         in_fl,        str,  par_opt,       None,      None ]
        can_call = [ False,         False,      False,    False,      False,      True ]

        ParSet.__init__(self, pars, values=values, dtypes=dtypes, can_call=can_call)


def validate_stellar_continuum_modeling_method_config(cnfg):
    """ 
    Validate the :class:`mangadap.util.parser.DefaultConfig` object with
    the stellar-continuum-modeling method parameters.

    Args:
        cnfg (:class:`mangadap.util.parser.DefaultConfig`): Object with
            the stellar-continuum-modeling method parameters to
            validate.

    Raises:
        KeyError: Raised if any required keywords do not exist.
        ValueError: Raised if keys have unacceptable values.
        FileNotFoundError: Raised if a file is specified but could not
            be found.
    """
    # Check for required keywords
    required_keywords = ['key', 'fit_type', 'fit_method' ]
    if not cnfg.all_required(required_keywords):
        raise KeyError('Keywords {0} must all have valid values.'.format(required_keywords))

    if cnfg['fit_method'] not in [ 'ppxf' ]:
        raise ValueError('Unknown fitting method: {0}'.format(cnfg['fit_method']))

    # Method specific keywords
    if cnfg['fit_method'] == 'ppxf' and not cnfg.keyword_specified('template_library'):
        raise KeyError('Keyword \'template_library\' must be provided for pPXF fitting.')


def available_stellar_continuum_modeling_methods(dapsrc=None):
    """
    Return the list of available stellar-continuum modeling methods.

    pPXF methods:

    .. todo::
        Fill in

    Args:
        dapsrc (str): (**Optional**) Root path to the DAP source
            directory.  If not provided, the default is defined by
            :func:`mangadap.config.defaults.dap_source_dir`.

    Returns:
        list: A list of
        :func:`mangadap.proc.stellarcontinuummodel.StellarContinuumModelDef`
        objects, each defining a stellar-continuum modeling approach.

    Raises:
        NotADirectoryError: Raised if the provided or default
            *dapsrc* is not a directory.
        OSError/IOError: Raised if no binning scheme configuration files
            could be found.
        KeyError: Raised if the binning method keywords are not all
            unique.

    .. todo::
        - Somehow add a python call that reads the databases and
          constructs the table for presentation in sphinx so that the
          text above doesn't have to be edited with changes in the
          available databases.
        
    """
    # Check the source directory exists
    dapsrc = dap_source_dir() if dapsrc is None else str(dapsrc)
    if not os.path.isdir(dapsrc):
        raise NotADirectoryError('{0} does not exist!'.format(dapsrc))

    # Check the configuration files exist
    search_dir = os.path.join(dapsrc, 'python/mangadap/config/stellar_continuum_modeling')
    ini_files = glob.glob(os.path.join(search_dir, '*.ini'))
    if len(ini_files) == 0:
        raise IOError('Could not find any configuration files in {0} !'.format(search_dir))

    # Build the list of library definitions
    # TODO: Should only read the keywords for the pixel mask so that the
    # artifact and emission-line databases are not read into memory
    modeling_methods = []
    for f in ini_files:
        # Read and validate the config file
        cnfg = DefaultConfig(f=f)
        validate_stellar_continuum_modeling_method_config(cnfg)

        # TODO: Pull this out; shouldn't be instantiating these already
        artifacts = None if cnfg['artifact_mask'] is None else \
                        ArtifactDB(cnfg['artifact_mask'], dapsrc=dapsrc)
        emission_lines = None if cnfg['emission_line_mask'] is None else \
                            EmissionLineDB(cnfg['emission_line_mask'], dapsrc=dapsrc)
        waverange = cnfg.getlist('waverange', evaluate=True)
        minimum_snr = cnfg.getfloat('minimum_snr', default=0.0)

        if cnfg['fit_method'] == 'ppxf':
            # sky is always none for now
            fitpar = PPXFFitPar(cnfg['template_library'], None, None, None,
                                iteration_mode=cnfg.get('fit_iter', default='global_template'),
                                reject_boxcar=cnfg.getint('reject_boxcar'),
                                filter_boxcar=cnfg.getint('filter_boxcar'),
                                filter_operation=cnfg.get('filter_op', default='subtract'),
                                filter_iterations=cnfg.getint('filter_iter'),
                                match_resolution=cnfg.getbool('match_resolution', default=True),
                                velscale_ratio=cnfg.getint('velscale_ratio', default=1),
                                minimum_snr=minimum_snr,
                                pixelmask=SpectralPixelMask(artdb=artifacts, emldb=emission_lines,
                                                            waverange=waverange),
                                bias=cnfg.getfloat('bias'), degree=cnfg.getint('degree'),
                                mdegree=cnfg.getint('mdegree'),
                                filt_degree=cnfg.getint('filter_degree'),
                                filt_mdegree=cnfg.getint('filter_mdegree'),
                                moments=cnfg.getint('moments') )
            fitclass = PPXFFit(StellarContinuumModelBitMask(dapsrc=dapsrc))
            fitfunc = fitclass.fit_SpatiallyBinnedSpectra
        else:
            raise ValueError('Unknown fitting method: {0}'.format(cnfg['default']['fit_method']))

        modeling_methods += [ StellarContinuumModelDef(cnfg['key'], 
                                                       cnfg.getfloat('minimum_snr', default=0.0),
                                                       cnfg['fit_type'],
                                                       fitpar, fitclass, fitfunc) ]

    # Check the keywords of the libraries are all unique
    if len(numpy.unique( numpy.array([ method['key'] for method in modeling_methods ]) )) \
            != len(modeling_methods):
        raise KeyError('Stellar-continuum modeling method keywords are not all unique!')

    # Return the default list of assessment methods
    return modeling_methods


class StellarContinuumModelBitMask(BitMask):
    r"""
    Derived class that specifies the mask bits for the stellar-continuum
    modeling.  See :class:`mangadap.util.bitmask.BitMask` for
    attributes.

    A list of the bits and meanings are provided by the base class
    function :func:`mangadap.util.bitmask.BitMask.info`; i.e.,::

        from mangadap.proc.stellarcontinuummodel import StellarContinuumModelBitMask
        bm = StellarContinuumModelBitMask()
        bm.info()

    """
    def __init__(self, dapsrc=None):
        dapsrc = dap_source_dir() if dapsrc is None else str(dapsrc)
        BitMask.__init__(self, ini_file=os.path.join(dapsrc, 'python', 'mangadap', 'config',
                                                     'bitmasks',
                                                     'stellar_continuum_model_bits.ini'))


class StellarContinuumModel:
    r"""
    Class that holds the stellar-continuum model results.

    Args:
        method_key (str): The keyword that designates which method,
            provided in *method_list*, to use for the continuum-fitting
            procedure.  
        binned_spectra
            (:class:`mangadap.proc.spatiallybinnedspectra.SpatiallyBinnedSpectra`):
            The binned spectra to fit.
        guess_vel (float, numpy.ndarray): A single or spectrum-dependent
            initial estimate of the redshift in km/s (:math:`cz`).
        guess_sig (float, numpy.ndarray): (**Optional**) A single or
            spectrum-dependent initial estimate of the velocity
            dispersion in km/s.  Default is 100 km/s.
        method_list (list): (**Optional**) List of
            :class:`StellarContinuumModelDef` objects that define one or
            more methods to use for the stellar continuum binning.  The
            default list is provided by the config files in the DAP
            source directory and compiled into this list using
            :func:`available_stellar_continuum_modeling_methods`.
        dapsrc (str): (**Optional**) Root path to the DAP source
            directory.  If not provided, the default is defined by
            :func:`mangadap.config.defaults.dap_source_dir`.
        dapver (str): (**Optional**) The DAP version to use for the
            analysis, used to override the default defined by
            :func:`mangadap.config.defaults.default_dap_version`.
        analysis_path (str): (**Optional**) The top-level path for the
            DAP output files, used to override the default defined by
            :func:`mangadap.config.defaults.default_analysis_path`.
        directory_path (str): The exact path to the directory with DAP
            output that is common to number DAP "methods".  See
            :attr:`directory_path`.
        output_file (str): (**Optional**) Exact name for the output
            file.  The default is to use
            :func:`mangadap.config.defaults.default_dap_file_name`.
        hardcopy (bool): (**Optional**) Flag to write the HDUList
            attribute to disk.  Default is True; if False, the HDUList
            is only kept in memory and would have to be reconstructed.
        tpl_symlink_dir (str): (**Optional**) Create a symbolic link to
            the created template library file in the supplied directory.
            Default is to produce no symbolic link.
        clobber (bool): (**Optional**) Overwrite any existing files.
            Default is to use any existing file instead of redoing the
            analysis and overwriting the existing output.
        checksum (bool): (**Optional**) Use the checksum in the fits
            header to confirm that the data has not been corrupted.  The
            checksum is **always** written to the fits header when the
            file is created; this argument does not toggle that
            functionality.
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

    .. todo::
        - Allow for a prior?

    """
#    @profile
    def __init__(self, method_key, binned_spectra, guess_vel, guess_sig=None,
                 method_list=None, dapsrc=None, dapver=None, analysis_path=None,
                 directory_path=None, output_file=None, hardcopy=True, tpl_symlink_dir=None,
                 clobber=False, checksum=False, loggers=None, quiet=False):

        self.loggers = None
        self.quiet = False

        # Define the method properties
        self.method = None
        self._define_method(method_key, method_list=method_list, dapsrc=dapsrc)

        self.binned_spectra = None
        self.guess_vel = None
        self.guess_sig = None

        # Define the output directory and file
        self.directory_path = None      # Set in _set_paths
        self.output_file = None
        self.hardcopy = None

        # Define the bitmask
        self.bitmask = StellarContinuumModelBitMask(dapsrc=dapsrc)

        # Initialize the main class attributes
        self.hdu = None
        self.checksum = checksum
        self.shape = None
        self.spatial_shape = None
        self.nspec = None
        self.spatial_index = None
        self.spectral_arrays = None
        self._assign_spectral_arrays()
        self.image_arrays = None
        self._assign_image_arrays()
        self.dispaxis = None
        self.nwave = None

        self.nmodels = None
        self.missing_models = None

        # TODO: Include covariance between measured properties

        # Run the assessments of the DRP file
        self.fit(binned_spectra, guess_vel, guess_sig=guess_sig, dapsrc=dapsrc,
                 dapver=dapver, analysis_path=analysis_path, directory_path=directory_path,
                 output_file=output_file, hardcopy=hardcopy, tpl_symlink_dir=tpl_symlink_dir,
                 clobber=clobber, loggers=loggers, quiet=quiet)


#    def __del__(self):
#        """
#        Deconstruct the data object by ensuring that the fits file is
#        properly closed.
#        """
#        if self.hdu is None:
#            return
#        self.hdu.close()
#        self.hdu = None
#

    def __getitem__(self, key):
        return self.hdu[key]


    def _define_method(self, method_key, method_list=None, dapsrc=None):
        r"""
        Select the method
        """
        # Grab the specific method
        self.method = select_proc_method(method_key, StellarContinuumModelDef,
                                         method_list=method_list,
                                       available_func=available_stellar_continuum_modeling_methods,
                                         dapsrc=dapsrc)


    def _fill_method_par(self, dapsrc=None, dapver=None, analysis_path=None, tpl_symlink_dir=None):
        """
        Fill in any remaining modeling parameters.

        This creates a full array for the guess redshifts and guess velocity dispersions.

        It also creates/reads the template library.

        .. todo:
            
            HARDCODED now to NEVER save the processed template library.
            Should add this as an option in StellarContinuumModelDef or PPXFFitPar.

        """
        # Report
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, 'Determining initial guess kinematics')

        # Fill the guess kinematics
        c = astropy.constants.c.to('km/s').value
        nbins = self.binned_spectra.nbins
        if isinstance(self.guess_vel, (list, numpy.ndarray)):
            if len(self.guess_vel) > 1 and len(self.guess_vel) != nbins:
                raise ValueError('Incorrect number of guess velocities provided.  Expect one per ' \
                                 'binned spectrum = {0}'.format(nbins))
            self.method['fitpar']['guess_redshift'] = numpy.asarray(self.guess_vel/c).astype(
                                                                                    numpy.float)
        else:
            self.method['fitpar']['guess_redshift'] = numpy.full(nbins, self.guess_vel/c,
                                                                 dtype=numpy.float)

        if isinstance(self.guess_sig, (list, numpy.ndarray)):
            if len(self.guess_sig) > 1 and len(self.guess_sig) != nbins:
                raise ValueError('Incorrect number of guess velocity dispersions provided.  ' \
                                 'Expect one per binned spectrum = {0}'.format(nbins))
            self.method['fitpar']['guess_dispersion'] = numpy.asarray(self.guess_sig).astype(
                                                                                    numpy.float)
        else:
            self.method['fitpar']['guess_dispersion'] = numpy.full(nbins, self.guess_sig,
                                                                   dtype=numpy.float)

        if self.method['fitpar']['template_library_key'] is not None:
            if not self.quiet:
                log_output(self.loggers, 1, logging.INFO, 'Instantiating template library...')
            self.method['fitpar']['template_library'] = self.get_template_library(
                            dapsrc=dapsrc, dapver=dapver, analysis_path=analysis_path,
                            tpl_symlink_dir=tpl_symlink_dir,
                            velocity_offset=numpy.mean(c*self.method['fitpar']['guess_redshift']),
                            match_to_drp_resolution=self.method['fitpar']['match_resolution'],
                            hardcopy=False)


    def get_template_library(self, dapsrc=None, dapver=None, analysis_path=None,
                             tpl_symlink_dir=None, velocity_offset=None,
                             match_to_drp_resolution=False, resolution_fwhm=None,
                             hardcopy=False):

        if resolution_fwhm is None:
            return TemplateLibrary(self.method['fitpar']['template_library_key'],
                                   velocity_offset=velocity_offset, drpf=self.binned_spectra.drpf,
                                   match_to_drp_resolution=match_to_drp_resolution,
                                   velscale_ratio=self.method['fitpar']['velscale_ratio'],
                                   dapsrc=dapsrc, analysis_path=analysis_path,
                                   hardcopy=hardcopy, symlink_dir=tpl_symlink_dir,
                                   loggers=self.loggers, quiet=self.quiet)
        else:
            # Set the spectral resolution
            wave = self.binned_spectra['WAVE'].data
            sres = SpectralResolution(wave, wave/resolution_fwhm, log10=True)

            # Set the file designation for this template library.
            # TODO: Move this to config.default?
            designation = '{0}-FWHM{1:.2f}'.format(self.method['fitpar']['template_library_key'],
                                                   resolution_fwhm)
            # Set the output file name
            processed_file = default_dap_file_name(self.binned_spectra.drpf.plate,
                                                   self.binned_spectra.drpf.ifudesign, designation)

            directory_path = default_dap_common_path(plate=self.binned_spectra.drpf.plate,
                                                     ifudesign=self.binned_spectra.drpf.ifudesign,
                                                     drpver=self.binned_spectra.drpf.drpver,
                                                     dapver=dapver, analysis_path=analysis_path)

            velscale_ratio = 1 if self.method['fitpar']['velscale_ratio'] is None \
                                else self.method['fitpar']['velscale_ratio']
            spectral_step=spectral_coordinate_step(wave, log=True) / velscale_ratio

            return TemplateLibrary(self.method['fitpar']['template_library_key'], sres=sres,
                                   velocity_offset=velocity_offset, spectral_step=spectral_step,
                                   log=True, dapsrc=dapsrc, directory_path=directory_path,
                                   hardcopy=hardcopy, symlink_dir=tpl_symlink_dir,
                                   processed_file=processed_file, loggers=self.loggers,
                                   quiet=self.quiet)


    def _set_paths(self, directory_path, dapver, analysis_path, output_file):
        """
        Set the :attr:`directory_path` and :attr:`output_file`.  If not
        provided, the defaults are set using, respectively,
        :func:`mangadap.config.defaults.default_dap_method_path` and
        :func:`mangadap.config.defaults.default_dap_file_name`.

        Args:
            directory_path (str): The exact path to the DAP reduction
                assessments file.  See :attr:`directory_path`.
            dapver (str): DAP version.
            analysis_path (str): The path to the top-level directory
                containing the DAP output files for a given DRP and DAP
                version.
            output_file (str): The name of the file with the reduction assessments.
                See :func:`compute`.

        """
        self.directory_path, self.output_file \
                = StellarContinuumModel.default_paths(self.binned_spectra.drpf.plate,
                                                      self.binned_spectra.drpf.ifudesign,
                                                      self.binned_spectra.rdxqa.method['key'],
                                                      self.binned_spectra.method['key'],
                                                      self.method['key'],
                                                      directory_path=directory_path, dapver=dapver,
                                                      analysis_path=analysis_path,
                                                      output_file=output_file)


    def _initialize_primary_header(self, hdr=None):
        """
        Initialize the header of :attr:`hdu`.

        Returns:
            astropy.io.fits.Header : Edited header object.

        """
        # Copy the from the DRP and clean it
        if hdr is None:
            hdr = self.binned_spectra.drpf.hdu['PRIMARY'].header.copy()
            hdr = DAPFitsUtil.clean_dap_primary_header(hdr)
        
        # Add keywords specific to this object
        hdr['AUTHOR'] = 'Kyle B. Westfall <westfall@ucolick.org>'
        hdr['VSTEP'] = (spectrum_velocity_scale(self.binned_spectra['WAVE'].data),
                        'Velocity step per spectral channel.')
        hdr['SCKEY'] = (self.method['key'], 'Stellar-continuum modeling method keyword')
        hdr['SCMINSN'] = (self.method['minimum_snr'], 'Minimum S/N of spectrum to include')

        # Guess kinematics
        # TODO: These are currently single numbers, but they may
        # eventually be vectors!  Include guess kinematics as a map?
        if self.guess_vel is not None:
            hdr['SCINPVEL'] = (self.guess_vel, 'Initial guess velocity')
        if self.guess_sig is not None:
            hdr['SCINPSIG'] = (self.guess_sig, 'Initial guess velocity dispersion')

        # Include the number of models
        # TODO: Is it a problem if missing_models is too long?
        hdr['NSCMOD'] = (self.nmodels, 'Number of unique stellar-continuum models')
#        if len(self.missing_models) > 0:
#            hdr['EMPTYSC'] = (str(self.missing_models), 'Bins w/o models')
        return hdr


    def _add_method_header(self, hdr):
        """Add fitting method information to the header."""
        if self.method['fitclass'] is not None:
            try:
                hdr['SCTYPE'] = self.method['fitclass'].fit_type
                hdr['SCMETH'] = self.method['fitclass'].fit_method
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

        # Turn on the flag stating that the S/N in the spectrum was
        # below the requested limit
        low_snr = numpy.invert(self.binned_spectra['BINID'].data == self.hdu['BINID'].data)
        indx = numpy.array([low_snr]*self.nwave).transpose(1,2,0)
        mask[indx] = self.bitmask.turn_on(mask[indx], flag='LOW_SNR')

        return mask


#    def _check_snr(self):
#        """
#        Determine which spectra in :attr:`binned_spectra` have a S/N
#        greater than the minimum set by :attr:`method`.  Only these
#        spectra will be analyzed.
#    
#        Returns:
#            numpy.ndarray : Boolean array for the spectra that satisfy
#            the criterion.
#        """
##        print(self.binned_spectra['BINS'].data['SNR'])
##        print(self.method['minimum_snr'])
#        return self.binned_spectra.above_snr_limit(self.method['minimum_snr'])


#    def _bins_to_fit(self):
#        """Return flags for the bins to fit."""
#        return (self._check_snr()) \
#                    & ~(numpy.array([ b in self.binned_spectra.missing_bins \
#                                        for b in numpy.arange(self.binned_spectra.nbins)]))


    def _assign_spectral_arrays(self):
        """
        Set :attr:`spectral_arrays`, which contains the list of
        extensions in :attr:`hdu` that contain spectral data.
        """
        self.spectral_arrays = [ 'FLUX', 'MASK' ]


    def _assign_image_arrays(self):
        """
        Set :attr:`image_arrays`, which contains the list of extensions
        in :attr:`hdu` that are on-sky image data.
        """
        self.image_arrays = [ 'BINID' ]


    def _get_missing_models(self):
        good_snr = self.binned_spectra.above_snr_limit(self.method['minimum_snr'])
        return numpy.sort(self.binned_spectra['BINS'].data['BINID'][numpy.invert(good_snr)].tolist()
                                + self.binned_spectra.missing_bins) 


    def _construct_2d_hdu(self, good_snr, model_flux, model_mask, model_par):
        """
        Construct :attr:`hdu` that is held in memory for manipulation of
        the object.  See :func:`construct_3d_hdu` if you want to convert
        the object into a DRP-like datacube.
        """
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, 'Constructing hdu ...')

        # Initialize the headers
        pri_hdr = self._initialize_primary_header()
        pri_hdr = self._add_method_header(pri_hdr)
        map_hdr = DAPFitsUtil.build_map_header(self.binned_spectra.drpf,
                                               'K Westfall <westfall@ucolick.org>')

        # Get the spatial map mask
        map_mask = numpy.zeros(self.spatial_shape, dtype=self.bitmask.minimum_dtype())
        # Add any spaxel not used because it was flagged by the binning
        # step
        indx = self.binned_spectra['MAPMASK'].data > 0
        map_mask[indx] = self.bitmask.turn_on(map_mask[indx], 'DIDNOTUSE')
        # Isolate any spaxels with foreground stars
        indx = self.binned_spectra.bitmask.flagged(self.binned_spectra['MAPMASK'].data, 'FORESTAR')
        map_mask[indx] = self.bitmask.turn_on(map_mask[indx], 'FORESTAR')
        # Get the bins that were blow the S/N limit
        indx = numpy.invert(DAPFitsUtil.reconstruct_map(self.spatial_shape,
                                                        self.binned_spectra['BINID'].data.ravel(),
                                                        good_snr, dtype='bool')) & (map_mask == 0)
        map_mask[indx] = self.bitmask.turn_on(map_mask[indx], 'LOW_SNR')

        # Get the bin ids with fitted models
        bin_indx = DAPFitsUtil.downselect_bins(self.binned_spectra['BINID'].data.ravel(),
                                               model_par['BINID']).reshape(self.spatial_shape)
#        vmin = numpy.amin(self.binned_spectra['BINID'].data)
#        vmax = numpy.amax(self.binned_spectra['BINID'].data)
#        pyplot.imshow(self.binned_spectra['BINID'].data, origin='lower', interpolation='nearest',
#                      vmin=vmin, vmax=vmax)
#        pyplot.show()
#        pyplot.imshow(bin_indx.reshape(self.spatial_shape), origin='lower', interpolation='nearest',
#                      vmin=vmin, vmax=vmax)
#        pyplot.show()

        # Save the data to the hdu attribute
        self.hdu = fits.HDUList([ fits.PrimaryHDU(header=pri_hdr),
                                  fits.ImageHDU(data=model_flux.data, name='FLUX'),
                                  fits.ImageHDU(data=model_mask, name='MASK'),
                                  self.binned_spectra['WAVE'].copy(),
                                  fits.ImageHDU(data=bin_indx, header=map_hdr, name='BINID'),
                                  fits.ImageHDU(data=map_mask, header=map_hdr, name='MAPMASK'),
                                  fits.BinTableHDU.from_columns( [ fits.Column(name=n,
                                                             format=rec_to_fits_type(model_par[n]),
                                            array=model_par[n]) for n in model_par.dtype.names ],
                                                               name='PAR')
                                ])

    @staticmethod
    def default_paths(plate, ifudesign, rdxqa_method, binning_method, method_key,
                      directory_path=None, drpver=None, dapver=None, analysis_path=None,
                      output_file=None):
        """
        Set the default directory and file name for the output file.

        Args:
            plate (int): Plate number of the observation.
            ifudesign (int): IFU design number of the observation.
            rdxqa_method (str): Keyword designating the method used to
                construct the reductions assessments.
            bin_method (str): Keyword designating the method use for the
                spatial binning.
            method_key (str): Keyword designating the method used for
                the stellar-continuum fitting.
            directory_path (str): (**Optional**) The exact path to the
                DAP stellar-continuum file.  Default set by
                :func:`mangadap.config.defaults.default_dap_method_path`.
            drpver (str): (**Optional**) DRP version.  Default set by
                :func:`mangadap.config.defaults.default_drp_version`.
            dapver (str): (**Optional**) DAP version.  Default set by
                :func:`mangadap.config.defaults.default_dap_version`.
            analysis_path (str): (**Optional**) The path to the
                top-level directory containing the DAP output files for
                a given DRP and DAP version. Default set by
                :func:`mangadap.config.defaults.default_analysis_path`.
            output_file (str): (**Optional**) The name of the file with
                the reduction assessments.  Default set by
                :func:`mangadap.config.defaults.default_dap_file_name`.

        Returns:
            str: Two strings with the path for the output file and the
            name of the output file.
        """
        # Set the output directory path
        method = default_dap_method(binning_method=binning_method, continuum_method=method_key)
        directory_path = default_dap_method_path(method, plate=plate, ifudesign=ifudesign,
                                                 ref=True, drpver=drpver, dapver=dapver,
                                                 analysis_path=analysis_path) \
                                        if directory_path is None else str(directory_path)

        # Set the output file
        ref_method = '{0}-{1}-{2}'.format(rdxqa_method, binning_method, method_key)
        output_file = default_dap_file_name(plate, ifudesign, ref_method) \
                                        if output_file is None else str(output_file)
        return directory_path, output_file


    def file_name(self):
        """Return the name of the output file."""
        return self.output_file


    def file_path(self):
        """Return the full path to the output file."""
        if self.directory_path is None or self.output_file is None:
            return None
        return os.path.join(self.directory_path, self.output_file)

    
    def info(self):
        return self.hdu.info()


    def all_spectrum_flags(self):
        return ['DIDNOTUSE', 'FORESTAR', 'LOW_SNR', 'ARTIFACT', 'OUTSIDE_RANGE', 'EML_REGION',
                'TPL_PIXELS', 'TRUNCATED', 'PPXF_REJECT', 'INVALID_ERROR', 'FIT_FAILED',
                'NEAR_BOUND' ]


    def all_except_emission_flags(self):
        return ['DIDNOTUSE', 'FORESTAR', 'LOW_SNR', 'ARTIFACT', 'OUTSIDE_RANGE', 'TPL_PIXELS',
                'TRUNCATED', 'PPXF_REJECT', 'INVALID_ERROR', 'FIT_FAILED', 'NEAR_BOUND' ]


    def fit(self, binned_spectra, guess_vel, guess_sig=None, dapsrc=None, dapver=None,
            analysis_path=None, directory_path=None, output_file=None, hardcopy=True,
            tpl_symlink_dir=None, clobber=False, loggers=None, quiet=False):
        """
        Fit the binned spectra using the specified method.

        Args:
            binned_spectra
                (:class:`mangadap.proc.spatiallybinnedspectra.SpatiallyBinnedSpectra`):
                The binned spectra to fit.
            guess_vel (float, numpy.ndarray): A single or spectrum-dependent
                initial estimate of the redshift in km/s (:math:`cz`).
            guess_sig (float, numpy.ndarray): (**Optional**) A single or
                spectrum-dependent initial estimate of the velocity
                dispersion in km/s.  Default is 100 km/s.
            dapsrc (str): (**Optional**) Root path to the DAP source
                directory.  If not provided, the default is defined by
                :func:`mangadap.config.defaults.dap_source_dir`.
            dapver (str): (**Optional**) The DAP version to use for the
                analysis, used to override the default defined by
                :func:`mangadap.config.defaults.default_dap_version`.
            analysis_path (str): (**Optional**) The top-level path for
                the DAP output files, used to override the default
                defined by
                :func:`mangadap.config.defaults.default_analysis_path`.
            directory_path (str): The exact path to the directory with
                DAP output that is common to number DAP "methods".  See
                :attr:`directory_path`.
            output_file (str): (**Optional**) Exact name for the output
                file.  The default is to use
                :func:`mangadap.config.defaults.default_dap_file_name`.
            hardcopy (bool): (**Optional**) Flag to write the HDUList
                attribute to disk.  Default is True; if False, the
                HDUList is only kept in memory and would have to be
                reconstructed.
            tpl_symlink_dir (str): (**Optional**) Create a symbolic link
                to the created template library file in the supplied
                directory.  Default is to produce no symbolic link.
            clobber (bool): (**Optional**) Overwrite any existing files.
                Default is to use any existing file instead of redoing
                the analysis and overwriting the existing output.
            loggers (list): (**Optional**) List of `logging.Logger`_
                objects to log progress; ignored if quiet=True.  Logging
                is done using :func:`mangadap.util.log.log_output`.
                Default is no logging.
            quiet (bool): (**Optional**) Suppress all terminal and
                logging output.  Default is False.

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

        self.shape = self.binned_spectra.shape
        self.spatial_shape =self.binned_spectra.spatial_shape
        self.nspec = self.binned_spectra.nspec
        self.spatial_index = self.binned_spectra.spatial_index.copy()
        self.dispaxis = self.binned_spectra.dispaxis
        self.nwave = self.binned_spectra.nwave
        
        # Get the guess kinematics
        if guess_vel is not None:
            self.guess_vel=guess_vel
        if self.guess_vel is None:
            raise ValueError('Must provide initial guess velocity.')
        self.guess_sig = 100.0 if guess_sig is None else guess_sig

        #---------------------------------------------------------------
        # Get the good spectra
        good_snr = self.binned_spectra.above_snr_limit(self.method['minimum_snr'])

        # Report
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, '-'*50)
            log_output(loggers, 1, logging.INFO, '{0:^50}'.format('STELLAR CONTINUUM FITTING'))
            log_output(self.loggers, 1, logging.INFO, '-'*50)
            log_output(self.loggers, 1, logging.INFO, 'Number of binned spectra: {0}'.format(
                                                            self.binned_spectra.nbins))
            if len(self.binned_spectra.missing_bins) > 0:
                log_output(self.loggers, 1, logging.INFO, 'Missing bins: {0}'.format(
                                                            len(self.binned_spectra.missing_bins)))
            log_output(self.loggers, 1, logging.INFO, 'With good S/N and to fit: {0}'.format(
                                                            numpy.sum(good_snr)))
#            log_output(self.loggers, 1, logging.INFO, 'Total to fit: {0}'.format(
#                                                            numpy.sum(good_bins)))

        # TODO: Is there a better way to handle this?
        if numpy.sum(good_snr) == 0:
            raise ValueError('No good spectra to fit!')

        #---------------------------------------------------------------
        # Fill in any remaining binning parameters
        self._fill_method_par(dapsrc=dapsrc, analysis_path=analysis_path,
                              tpl_symlink_dir=tpl_symlink_dir)

        # (Re)Set the output paths
        self._set_paths(directory_path, dapver, analysis_path, output_file)

        #---------------------------------------------------------------
        # Check that the file path is defined
        ofile = self.file_path()
        if ofile is None:
            raise ValueError('File path for output file is undefined!')

        # Report
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, 'Output path: {0}'.format(
                                                                            self.directory_path))
            log_output(self.loggers, 1, logging.INFO, 'Output file: {0}'.format(
                                                                            self.output_file))
        
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
        model_wave, model_flux, model_mask, model_par \
                = self.method['fitfunc'](self.binned_spectra, par=self.method['fitpar'],
                                         loggers=self.loggers, quiet=self.quiet)

        # The number of models returned should match the number of
        # "good" spectra

        # TODO: Include failed fits in "missing" models?

        # DEBUG
        if model_flux.shape[0] != numpy.sum(good_snr):
            raise ValueError('Unexpected returned shape of fitted continuum models.')

        # Set the number of models and the missing models
        self.nmodels = model_flux.shape[0]
        self.missing_models = self._get_missing_models()

        # Report
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, 'Fitted models: {0}'.format(self.nmodels))
            if len(self.missing_models) > 0:
                log_output(self.loggers, 1, logging.INFO, 'Missing models: {0}'.format(
                                                            len(self.missing_models)))

        # Construct the 2d hdu list that only contains the fitted models
        self.hardcopy = hardcopy
        self._construct_2d_hdu(good_snr, model_flux, model_mask, model_par)

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
            log_output(self.loggers, 1, logging.INFO, 'Constructing stellar continuum datacube ...')

        bin_indx = self.hdu['BINID'].data.copy()
        model_flux = self.hdu['FLUX'].data.copy()
        model_mask = self.hdu['MASK'].data.copy()

        flux, mask = DAPFitsUtil.reconstruct_cube(self.shape, bin_indx.ravel(),
                                                  [model_flux, model_mask])

        mask = self._finalize_cube_mask(mask)

        # Primary header is identical regardless of the shape of the
        # extensions
        hdr = self.hdu['PRIMARY'].header.copy()
        cube_hdr = DAPFitsUtil.build_cube_header(self.binned_spectra.drpf,
                                                 'K Westfall <westfall@ucolick.org>')

        # Return the converted hdu without altering the internal hdu
        return fits.HDUList([ fits.PrimaryHDU(header=hdr),
                              fits.ImageHDU(data=flux, header=cube_hdr, name='FLUX'),
                              fits.ImageHDU(data=mask, header=cube_hdr, name='MASK'),
                              self.hdu['WAVE'].copy(),
                              self.hdu['BINID'].copy(),
                              self.hdu['PAR'].copy()
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
        Read an existing file with an existing set of stellar-continuum
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

        # Confirm that the internal method is the same as the method
        # that was used in writing the file
        if self.hdu['PRIMARY'].header['SCKEY'] != self.method['key']:
            if strict:
                raise ValueError('Keywords in header do not match specified method keywords!')
            elif not self.quiet:
                warnings.warn('Keywords in header do not match specified method keywords!')
        # TODO: "strict" should also check other aspects of the file to
        # make sure that the details of the method are also the same,
        # not just the keyword

#        if not self.quiet:
#            log_output(self.loggers, 1, logging.INFO, 'Reverting to python-native structure.')
#        DAPFitsUtil.restructure_rss(self.hdu, ext=self.spectral_arrays)
#        DAPFitsUtil.restructure_map(self.hdu, ext=self.image_arrays)

        # Attempt to read the input guess velocity and sigma
        try:
            self.guess_vel = self.hdu['PRIMARY'].header['SCINPVEL']
            self.guess_sig = self.hdu['PRIMARY'].header['SCINPSIG']
        except:
            if not self.quiet:
                warnings.warn('Unable to read input guess kinematics from file header!')
            self.guess_vel = None
            self.guess_sig = None

        # Attempt to read the modeling parameters
        if self.method['fitpar'] is not None and callable(self.method['fitpar'].fromheader):
            self.method['fitpar'].fromheader(self.hdu['PRIMARY'].header)

        self.nmodels = self.hdu['FLUX'].shape[0]
        self.missing_models = self._get_missing_models()


    def copy_to_array(self, ext='FLUX', waverange=None, include_missing=False):
        r"""

        Wrapper for :func:`mangadap.util.fitsutil.DAPFitsUtil.copy_to_array`
        specific for :class:`StellarContinuumModel`.

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
                                         nbins=self.nbins,
                                    unique_bins=DAPFitsUtil.unique_bins(self.hdu['BINID'].data))


    def copy_to_masked_array(self, ext='FLUX', flag=None, waverange=None, include_missing=False):
        r"""
        Wrapper for
        :func:`mangadap.util.fitsutil.DAPFitsUtil.copy_to_masked_array`
        specific for :class:`StellarContinuumModel`.

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


    def construct_models(self, replacement_templates=None, redshift_only=False,
                         deredshift=False, corrected_dispersion=False):
        """
        Reconstruct the best fitting models.
        """
        if self.method['fitclass'] is None:
            raise ValueError('No class object available for constructing the model!')
        if not callable(self.method['fitclass'].construct_models):
            raise AttributeError('Provided fit class object has no callable ' \
                                 '\'construct_models\' attribute!')
        if replacement_templates is not None \
                and not isinstance(replacement_templates, TemplateLibrary):
            raise TypeError('Provided template replacements must have type TemplateLibrary.')

        # Only the selected models are constructed, others are masked
        select = numpy.invert(self.bitmask.flagged(self.hdu['PAR'].data['MASK'],
                                            flag=['NO_FIT', 'FIT_FAILED', 'INSUFFICIENT_DATA']))
        templates = self.method['fitpar']['template_library'] if replacement_templates is None \
                                            else replacement_templates

        # The only reason it needs the flux data is to get the shape, so
        # it doesn't matter that I'm passing the model flux
        return self.method['fitclass'].construct_models(templates['WAVE'].data,
                                                        templates['FLUX'].data,
                                                        self.hdu['WAVE'].data,
                                                        self.hdu['FLUX'].data.shape,
                                                        self.hdu['PAR'].data, select=select,
                                                        redshift_only=redshift_only,
                                                        deredshift=deredshift,
                                                        corrected_dispersion=corrected_dispersion,
                                                        dvtol=1e-9)


    @staticmethod
    def reset_continuum_mask_window(continuum, dispaxis=1, quiet=False):
        """
        Reset the mask of the stellar continuum to a continuous window
        from the minimum to maximum valid wavelength.
        """

        spatial_shape = DAPFitsUtil.get_spatial_shape(continuum.shape, dispaxis)
        if len(spatial_shape) != 1:
            raise ValueError('Input array should be two-dimensional!')

        # No pixels are masked!
        if numpy.sum(numpy.ma.getmaskarray(continuum)) == 0:
            return continuum

        _continuum = continuum.copy() if dispaxis == 1 else continuum.copy().T

        nspec,npix = _continuum.shape
        pix = numpy.ma.MaskedArray(numpy.array([ numpy.arange(npix) ]*nspec),
                                   mask=numpy.ma.getmaskarray(_continuum).copy())

        min_good_pix = numpy.ma.amin(pix, axis=dispaxis)
        max_good_pix = numpy.ma.amax(pix, axis=dispaxis)
        for c,s,e in zip(_continuum, min_good_pix,max_good_pix+1):
            if isinstance(s, numpy.ma.core.MaskedConstant) \
                    or isinstance(e, numpy.ma.core.MaskedConstant):
#                if not quiet:
#                    warnings.warn('Full continuum fit is masked.')
                continue
            c.mask[s:e] = False

        return _continuum if dispaxis == 1 else _continuum.T


    def unmasked_continuum_model(self, replacement_templates=None, redshift_only=False,
                                 deredshift=False, corrected_dispersion=False):
        """
        Return the stellar continuum unmasked over the full continuous
        fitting region.  Models are reconstructed based on the model
        parameters if new template fluxes are provided or if the
        returned models should not include the LOSVD convolution
        (redshift_only=True).
        """
        # Get the models for the binned spectra
        reconstruct = replacement_templates is not None or redshift_only or corrected_dispersion
        continuum = self.construct_models(replacement_templates=replacement_templates,
                                          redshift_only=redshift_only, deredshift=deredshift,
                                          corrected_dispersion=corrected_dispersion) \
                        if reconstruct \
                        else self.copy_to_masked_array(flag=self.all_spectrum_flags())
        return StellarContinuumModel.reset_continuum_mask_window(continuum, quiet=self.quiet)


#    def fill_to_match(self, binned_spectra, replacement_templates=None, redshift_only=False):
#        """
#        Get the stellar continuum that matches the shape of the provided binned_spectra.
#        """
#        if binned_spectra is self.binned_spectra:
#            continuum = self.unmasked_continuum_model(replacement_templates=replacement_templates,
#                                                        redshift_only=redshift_only)
#
#            # Number of models matches the numbers of bins
#            if binned_spectra.nbins == self.nmodels:
#                return continuum
#    
#            # Fill in bins with no models with masked zeros
#            _continuum = numpy.ma.zeros(binned_spectra['FLUX'].data.shape, dtype=float)
#            _continuum[:,:] = numpy.ma.masked
#            for i,j in enumerate(self.hdu['PAR'].data['BINID_INDEX']):
#                _continuum[j,:] = continuum[i,:]
#            return _continuum
#
#        raise NotImplementedError('Can only match to internal binned_spectra.')

    def fill_to_match(self, binid, replacement_templates=None, redshift_only=False,
                      deredshift=False, corrected_dispersion=False, missing=None):
        """
        Get the stellar-continuum model from this objects that matches
        the input bin ID matrix.  The output is a 2D matrix ordered by
        the bin ID; any skipped index numbers in the maximum of the
        union of the unique numbers in the `binid` and `missing` input
        are masked.

        All this does is find the stellar continuum from spaxels in this
        model and match them to the input spaxels.  Depending on how
        differently the data data sets have been binned, this says
        nothing about whether or not the returned models are a good fit
        to the data!

        replacement_templates only works if the number of templates is
        the same as those used during the actual fitting process.

        use redshift_only if you want to exclude the best-fitting
        dispersion from the model.

        """
        # The input bin id array must have the same shape as self
        if binid.shape != self.spatial_shape:
            raise ValueError('Input bin ID matrix has incorrect shape.')

        # Construct the best-fitting models
        best_fit_continuum = self.unmasked_continuum_model(
                                    replacement_templates=replacement_templates,
                                    redshift_only=redshift_only, deredshift=deredshift,
                                    corrected_dispersion=corrected_dispersion)

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


    def matched_kinematics(self, binid, redshift=None, dispersion=100.0, constant=False, cz=False,
                           corrected=False, min_dispersion=None, nearest=False, missing=None):
        """
        Return the guess redshift and velocity dispersion for all the
        spectra based on the stellar kinematics.

        For spectra the were not fit, use the median (unmasked)
        kinematic measurement if no default value is provided as an
        argument.

        Args:
            binid (numpy.ndarray): 2D array with the bin ID associate
                with each spaxel in the datacube.  Shape must be the
                same as :attr:`spatial_shape`.
            redshift (float): (**Optional**) The default redshift to use
                for spectra without a stellar-continuum fit.  Default is
                to use the median of the unmasked measurements
            dispersion (float): (**Optional**) The default velocity
                dispersion to use for spectra without a
                stellar-continuum fit.  Default is 100 km/s.
            constant (bool): (**Optional**) Force the function to return
                a constant redshift and dispersion for each spectrum,
                regardless of any fitted kinematics.
            cz (bool): (**Optional**) Return the redshift as cz in km/s,
                as opposed to unitless z.
            corrected (bool): (**Optional**) By default, the returned
                velocity dispersions are the measured values, which will
                include any resolution difference between the templates
                and the galaxy data.  Setting corrected to True will
                return the corrected dispersions.
            min_dispersion (float): (**Optional**) Impose a minimum
                dispersion.
            nearest (bool): (**Optional**) Instead of the median of the
                results for the spectra that were not fit, use the value
                from the nearest bin.
            missing (list): (**Optional**) Any bin ID numbers missing
                from the input `binid` image needed for constructing the
                output matched data.

        Returns:
            numpy.ndarray: Returns (unmasked!) arrays with a redshift
            (or cz) and dispersion for each binned spectrum.

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

        # Get the number of output kinematics
        nbins = numpy.amax(binid).astype(int)+1 if missing is None or len(missing) == 0 else \
                    numpy.amax( numpy.append(binid.ravel(), missing) ).astype(int)+1

        # Mask bad values
        mask = self.bitmask.flagged(self.hdu['PAR'].data['MASK'].copy(),
                                    [ 'NO_FIT', 'FIT_FAILED', 'INSUFFICIENT_DATA', 'NEAR_BOUND' ])

        if numpy.all(mask):
            # Everything is masked so use input guesses
            warnings.warn('All stellar continuum fits have been masked!  Using input guesses.')
            str_z = numpy.ma.MaskedArray(self.method['fitpar']['guess_redshift'].copy())
            str_d = numpy.ma.MaskedArray(self.method['fitpar']['guess_dispersion'].copy())
            # Select one per *stellar-continuum model*
            str_z = str_z[self['PAR'].data['BINID_INDEX']]
            str_d = str_d[self['PAR'].data['BINID_INDEX']]
        else:
            # Pull out the data
            str_z = numpy.ma.MaskedArray(self.hdu['PAR'].data['KIN'][:,0].copy(), mask=mask) \
                            / astropy.constants.c.to('km/s').value
            str_d = numpy.ma.MaskedArray(self.hdu['PAR'].data['KIN'][:,1].copy(), mask=mask)

            # Apply the sigma corrections if requested.  Values with the
            # correction larger than the measured value will be masked!
            if corrected:
                sigma_corr = numpy.ma.MaskedArray(self.hdu['PAR'].data['SIGMACORR_EMP'], mask=mask)
                str_d = numpy.ma.sqrt( numpy.square(str_d) - numpy.square(sigma_corr) )

        # Set the default values to use when necessary
        _redshift = numpy.ma.median(str_z) if redshift is None else redshift
        _dispersion = numpy.ma.median(str_d) if dispersion is None else dispersion
        if min_dispersion is not None and _dispersion < min_dispersion:
            _dispersion = min_dispersion

        # Just return the single value
        if constant:
            str_z = numpy.full(nbins, _redshift * astropy.constants.c.to('km/s').value
                                        if cz else _redshift, dtype=float)
            str_d = numpy.full(nbins, _dispersion, dtype=float)
            return str_z, str_d

        # Fill masked values with the nearest datum
        if nearest:
            # Fill masked values with the nearest datum
            best_fit_kinematics = numpy.ma.append([str_z], [str_d], axis=0).T
            valid_bins = numpy.unique(self['BINID'].data)[1:]
            coo = self.binned_spectra['BINS'].data['SKY_COO'][valid_bins,:]
            replace = str_z.mask | str_d.mask
            kinematics = replace_with_data_from_nearest_coo(coo, best_fit_kinematics, replace)
            str_z = kinematics[:,0]
            str_d = kinematics[:,1]
        else:
            # Fill any masked values with the single estimate
            str_z = str_z.filled(_redshift)
            str_d = str_d.filled(_dispersion)

        # Map the BINID to the spectrum index, assuming bins are sorted
        # and that the BINID map has -1 BINID values
        u, indx, reconstruct = numpy.unique(self['BINID'].data.ravel(), return_index=True,
                                            return_inverse=True)
        u_bin_indx = numpy.arange(len(u))-1
        _bin_indx = u_bin_indx[reconstruct].reshape(self.spatial_shape)

        # Match the kinematics to the output bin ID map
        _str_z = numpy.ma.masked_all(nbins, dtype=float)
        _str_d = numpy.ma.masked_all(nbins, dtype=float)
        for i,j in zip(binid.ravel(), _bin_indx.ravel()):
            if i < 0 or j < 0:
                continue
            _str_z[i] = str_z[j]
            _str_d[i] = str_d[j]

        str_z = _str_z.filled(_redshift)
        str_d = _str_d.filled(_dispersion)
        # Convert to cz velocities (km/s)
        if cz:
            str_z *= astropy.constants.c.to('km/s').value
        return str_z, str_d


    def matched_template_flags(self, binned_spectra):
        """
        TODO: Function out of date!

        Return the templates used during the fit to each spectrum,
        matched to the spectra in the binned_spectra object.  For
        spectra with no models, the flags are all set to true.
        """
        if binned_spectra is self.binned_spectra:
            usetpl = self.hdu['PAR'].data['TPLWGT'] > 0

            # Number of models matches the numbers of bins
            if binned_spectra.nbins == self.nmodels:
#                print('returning usetpl')
                return usetpl
    
            # Fill in bins with no models with masked zeros
#            print('Fill in bins with no models with masked zeros')
            ntpl = self.method['fitpar']['template_library'].ntpl
            _usetpl = numpy.ones((binned_spectra.nbins,ntpl), dtype=bool)
            for i,j in enumerate(self.hdu['PAR'].data['BINID_INDEX']):
                _usetpl[j,:] = usetpl[i,:]
            return _usetpl

        raise NotImplementedError('Can only match to internal binned_spectra.')


