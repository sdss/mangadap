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
            try:
                from configparser import ConfigParser
            except ImportError:
                warnings.warn('Unable to import configparser!  Beware!')
        else:
            try:
                from ConfigParser import ConfigParser
            except ImportError:
                warnings.warn('Unable to import ConfigParser!  Beware!')

        import os
        import glob
        import logging
        import numpy

        from astropy.io import fits
        import astropy.constants

        from ..mangafits import MaNGAFits
        from ..drpfits import DRPFits
        from ..par.parset import ParSet
        from ..util.log import log_output
        from ..util.fileio import rec_to_fits_type, write_hdu
        from ..util.instrument import spectrum_velocity_scale
        from ..util.bitmask import BitMask
        from ..config.defaults import default_dap_source, default_dap_file_name
        from ..config.defaults import default_dap_method, default_dap_method_path
        from .spatiallybinnedspectra import SpatiallyBinnedSpectra
        from .templatelibrary import TemplateLibrary
        from .ppxffit import PPXFFitPar, PPXFFit
        from .artifactdb import ArtifactDB
        from .emissionlinedb import EmissionLineDB
        from .pixelmask import SpectralPixelMask
        from .util import _select_proc_method

*Class usage examples*:
        Add examples

*Revision history*:
    | **14 Apr 2016**: Implementation begun by K. Westfall (KBW)
    | **19 Apr 2016**: (KBW) First version
    | **19 May 2016**: (KBW) Added loggers and quiet keyword arguments
        to :class:`StellarContinuumModel`, removed verbose 
    | **05 Jul 2016**: (KBW) Removed oversample keyword from instantion
        of :class:`mangadap.proc.ppxffit.PPXFFit` objects.

.. _astropy.io.fits.hdu.hdulist.HDUList: http://docs.astropy.org/en/v1.0.2/io/fits/api/hdulists.html
.. _glob.glob: https://docs.python.org/3.4/library/glob.html
.. _configparser.ConfigParser: https://docs.python.org/3/library/configparser.html#configparser.ConfigParser
.. _logging.Logger: https://docs.python.org/3/library/logging.html

"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import sys
import warnings
if sys.version > '3':
    long = int
    try:
        from configparser import ConfigParser
    except ImportError:
        warnings.warn('Unable to import configparser!  Beware!')
else:
    try:
        from ConfigParser import ConfigParser
    except ImportError:
        warnings.warn('Unable to import ConfigParser!  Beware!')

import os
import glob
import logging
import numpy

from astropy.io import fits
import astropy.constants

from ..mangafits import MaNGAFits
from ..drpfits import DRPFits
from ..par.parset import ParSet
from ..util.log import log_output
from ..util.fileio import rec_to_fits_type, write_hdu
from ..util.instrument import spectrum_velocity_scale
from ..util.bitmask import BitMask
from ..config.defaults import default_dap_source, default_dap_file_name
from ..config.defaults import default_dap_method, default_dap_method_path
from .spatiallybinnedspectra import SpatiallyBinnedSpectra
from .templatelibrary import TemplateLibrary
from .ppxffit import PPXFFitPar, PPXFFit
from .artifactdb import ArtifactDB
from .emissionlinedb import EmissionLineDB
from .pixelmask import SpectralPixelMask
from .util import _select_proc_method

from matplotlib import pyplot

__author__ = 'Kyle B. Westfall'
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
    Validate the `configparser.ConfigParser`_ object that is meant to
    define a stellar continuum modeling method.

    Args:
        cnfg (`configparser.ConfigParser`_): Object meant to contain
            defining parameters of the binning method as needed by
            :class:`mangadap.proc.stellarcontinuummodel.StellarContinuumModelDef

    Raises:
        KeyError: Raised if any required keywords do not exist.
        ValueError: Raised if keys have unacceptable values.
        FileNotFoundError: Raised if a file is specified but could not
            be found.
    """
    # Check for required keywords
    if 'key' not in cnfg.options('default'):
        raise KeyError('Keyword \'key\' must be provided.')
    if 'fit_type' not in cnfg.options('default'):
        raise KeyError('Keyword \'fit_type\' must be provided.')
    if 'fit_method' not in cnfg.options('default'):
        raise KeyError('Keyword \'fit_method\' must be provided.')

    if 'minimum_snr' not in cnfg.options('default') or cnfg['default']['minimum_snr'] is None:
        cnfg['default']['minimum_snr']= '0.0'

    if 'waverange' not in cnfg.options('default') or cnfg['default']['waverange'] is None:
        cnfg['default']['waverange']= 'None'

    if 'artifact_mask' not in cnfg.options('default') or cnfg['default']['artifact_mask'] is None:
        cnfg['default']['artifact_mask'] = 'None'

    if 'emission_line_mask' not in cnfg.options('default') \
            or cnfg['default']['emission_line_mask'] is None:
        cnfg['default']['emission_line_mask'] = 'None'

    if cnfg['default']['fit_method'] == 'ppxf':
        if 'template_library' not in cnfg.options('default'):
            raise KeyError('Keyword \'template_library\' must be provided for pPXF fitting.')
        if 'sky_spectra' in cnfg.options('default') and cnfg['default']['sky_spectra'] != 'None':
            raise NotImplementedError('Cannot yet use sky spectra via a config file.')
        if 'fit_iter' not in cnfg.options('default'):
            cnfg['default']['fit_iter'] = 'None'
        if 'match_resolution' not in cnfg.options('default'):
            cnfg['default']['match_resolution'] = 'None'
        if 'velscale_ratio' not in cnfg.options('default'):
            cnfg['default']['velscale_ratio'] = 'None'
        for k in PPXFFitPar._keyword_defaults().keys():
            if k not in cnfg.options('default'):
                cnfg['default'][k] = 'None'
        return

    raise ValueError('{0} is not a recognized fitting method.'.format(
                                                                    cnfg['default']['fit_method']))


def available_stellar_continuum_modeling_methods(dapsrc=None):
    """
    Return the list of available stellar-continuum modeling methods.

    pPXF methods:

    .. todo::
        Fill in

    Args:
        dapsrc (str): (**Optional**) Root path to the DAP source
            directory.  If not provided, the default is defined by
            :func:`mangadap.config.defaults.default_dap_source`.

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
        NameError: Raised if ConfigParser is not correctly imported.

    .. todo::
        - Somehow add a python call that reads the databases and
          constructs the table for presentation in sphinx so that the
          text above doesn't have to be edited with changes in the
          available databases.
        
    """
    # Check the source directory exists
    dapsrc = default_dap_source() if dapsrc is None else str(dapsrc)
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
        # Read the config file
        cnfg = ConfigParser(allow_no_value=True)
        cnfg.read(f)
        # Ensure it has the necessary elements to define the template
        # library
        validate_stellar_continuum_modeling_method_config(cnfg)

        artifacts = ArtifactDB(cnfg['default']['artifact_mask'], dapsrc=dapsrc) \
                        if cnfg['default']['artifact_mask'] != 'None' else None
        emission_lines = EmissionLineDB(cnfg['default']['emission_line_mask'], dapsrc=dapsrc) \
                        if cnfg['default']['emission_line_mask'] != 'None' else None

        if cnfg['default']['fit_method'] == 'ppxf':
            # sky is always none for now
            fitpar = PPXFFitPar(cnfg['default']['template_library'], None, None, None,
                                cnfg['default']['fit_iter'],
                                eval(cnfg['default']['match_resolution']),
                                eval(cnfg['default']['velscale_ratio']),
                                eval(cnfg['default']['minimum_snr']),
                                SpectralPixelMask(artdb=artifacts, emldb=emission_lines,
                                                  waverange=eval(cnfg['default']['waverange'])),
                                eval(cnfg['default']['bias']), eval(cnfg['default']['clean']),
                                eval(cnfg['default']['degree']), eval(cnfg['default']['mdegree']),
                                eval(cnfg['default']['moments']), None,
                                eval(cnfg['default']['regul']), eval(cnfg['default']['reddening']),
                                eval(cnfg['default']['component']),
                                eval(cnfg['default']['reg_dim']) )
            fitclass = PPXFFit(StellarContinuumModelBitMask(dapsrc=dapsrc))
            fitfunc = fitclass.fit_SpatiallyBinnedSpectra
        else:
            raise ValueError('Unknown fitting method: {0}'.format(cnfg['default']['fit_method']))

        modeling_methods += [ StellarContinuumModelDef(cnfg['default']['key'], 
                                                       cnfg['default'].getfloat('minimum_snr'),
                                                       cnfg['default']['fit_type'],
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
        dapsrc = default_dap_source() if dapsrc is None else str(dapsrc)
        BitMask.__init__(self, ini_file=os.path.join(dapsrc, 'python', 'mangadap', 'config',
                                                     'bitmasks',
                                                     'stellar_continuum_model_bits.ini'))


class StellarContinuumModel:
    r"""

    Class that holds the stellar-continuum model results.

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

    .. todo::
        - Allow for a prior?

    """
    def __init__(self, method_key, binned_spectra, guess_vel=None, guess_sig=None,
                 method_list=None, dapsrc=None, dapver=None, analysis_path=None,
                 directory_path=None, output_file=None, hardcopy=True, tpl_symlink_dir=None,
                 clobber=False, checksum=False, loggers=None, quiet=False):

        self.version = '1.0'
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
        self.tpl_symlink_dir = None

        # Initialize the objects used in the assessments
        self.bitmask = StellarContinuumModelBitMask(dapsrc=dapsrc)
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

        # Run the assessments of the DRP file
        self.fit(binned_spectra, guess_vel=guess_vel, guess_sig=guess_sig, dapsrc=dapsrc,
                 dapver=dapver, analysis_path=analysis_path, directory_path=directory_path,
                 output_file=output_file, hardcopy=hardcopy, tpl_symlink_dir=tpl_symlink_dir,
                 clobber=clobber, loggers=loggers, quiet=quiet)


    def __del__(self):
        """
        Deconstruct the data object by ensuring that the fits file is
        properly closed.
        """
        if self.hdu is None:
            return
        self.hdu.close()
        self.hdu = None


    def __getitem__(self, key):
        return self.hdu[key]


    def _define_method(self, method_key, method_list=None, dapsrc=None):
        r"""

        Select the method

        """
        # Grab the specific method
        self.method = _select_proc_method(method_key, StellarContinuumModelDef,
                                          method_list=method_list,
                                        available_func=available_stellar_continuum_modeling_methods,
                                          dapsrc=dapsrc)


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
        # Set the output directory path
        method = default_dap_method(binned_spectra=self.binned_spectra, stellar_continuum=self)
        #method = '{0}-{1}'.format(self.binned_spectra.method['key'], self.method['key'])
        self.directory_path = default_dap_method_path(method, plate=self.binned_spectra.drpf.plate,
                                                      ifudesign=self.binned_spectra.drpf.ifudesign,
                                                      ref=True,
                                                      drpver=self.binned_spectra.drpf.drpver,
                                                      dapver=dapver, analysis_path=analysis_path) \
                                        if directory_path is None else str(directory_path)

        # Set the output file
        ref_method = '{0}-{1}-{2}'.format(self.binned_spectra.rdxqa.method['key'],
                                      self.binned_spectra.method['key'], self.method['key'])
        self.output_file = default_dap_file_name(self.binned_spectra.drpf.plate,
                                                 self.binned_spectra.drpf.ifudesign, ref_method) \
                                        if output_file is None else str(output_file)


    def _fill_method_par(self, dapsrc=None, analysis_path=None):
        """
        Fill in any remaining modeling parameters.
        """
        # Report
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, 'Determining initial guess kinematics')

        # Fill the guess kinematics
        c = astropy.constants.c.to('km/s').value
        if isinstance(self.guess_vel, (list, numpy.ndarray)):
            if len(self.guess_vel) > 1 and len(self.guess_vel) != self.nmodels:
                raise ValueError('Incorrect number of guess velocities provided.  Expect one per ' \
                                 'binned spectrum = {0}'.format(self.nmodels))
            self.method['fitpar']['guess_redshift'] = numpy.asarray(self.guess_vel/c).astype(
                                                                                    numpy.float)
        else:
            self.method['fitpar']['guess_redshift'] = numpy.full(self.nmodels, self.guess_vel/c,
                                                                 dtype=numpy.float)

        if isinstance(self.guess_sig, (list, numpy.ndarray)):
            if len(self.guess_sig) > 1 and len(self.guess_sig) != self.binned_spectra.nmodels:
                raise ValueError('Incorrect number of guess velocity dispersions provided.  ' \
                                 'Expect one per binned spectrum = {0}'.format(self.nmodels))
            self.method['fitpar']['guess_dispersion'] = numpy.asarray(self.guess_sig).astype(
                                                                                    numpy.float)
        else:
            self.method['fitpar']['guess_dispersion'] = numpy.full(self.nmodels, self.guess_sig,
                                                                   dtype=numpy.float)

        if self.method['fitpar']['template_library_key'] is not None: 
            if not self.quiet:
                log_output(self.loggers, 1, logging.INFO, 'Instantiating template library...')
            self.method['fitpar']['template_library'] \
                    = TemplateLibrary(self.method['fitpar']['template_library_key'],
                            velocity_offset=numpy.mean(c*self.method['fitpar']['guess_redshift']),
                                      drpf=self.binned_spectra.drpf,
                                match_to_drp_resolution=self.method['fitpar']['match_resolution'],
                                      velscale_ratio=self.method['fitpar']['velscale_ratio'],
                                      dapsrc=dapsrc, analysis_path=analysis_path,
                                      symlink_dir=self.tpl_symlink_dir, loggers=self.loggers,
                                      quiet=self.quiet)
                                      #, clobber=True)


    def _clean_drp_header(self, ext='FLUX'):
        """
        Read and clean the header from the DRP fits file extension for
        use in the output file for this class object.

        .. todo::
            - Currently just returns the existing header.  Need to
              decide if/how to manipulate DRP header.

        """
        return self.binned_spectra.drpf[ext].header.copy()


    def _initialize_header(self, hdr):
        """

        Initialize the header.

        """
        hdr['VSTEP'] = (spectrum_velocity_scale(self.binned_spectra['WAVE'].data),
                        'Velocity step per spectral channel.')
        hdr['SCKEY'] = (self.method['key'], 'Stellar-continuum modeling method keyword')
        hdr['SCMINSN'] = (self.method['minimum_snr'], 'Minimum S/N of spectrum to include')
        # Guess kinematics may eventually be vectors!
        if self.guess_vel is not None:
            hdr['INPVEL'] = (self.guess_vel, 'Initial guess velocity')
        if self.guess_sig is not None:
            hdr['INPSIG'] = (self.guess_sig, 'Initial guess velocity dispersion')
        hdr['NMOD'] = (self.nmodels, 'Number of unique stellar-continuum models')
        if len(self.missing_models) > 0:
            hdr['EMPTYMOD'] = (str(self.missing_models), 'List of models with no data')
        return hdr


    def _initialize_mask(self, good_snr):
        """

        Initialize the mask be setting the DIDNOTUSE, FORESTAR, and LOW_SNR masks

        """
        # Initialize to all zeros
        mask = numpy.zeros(self.shape, dtype=self.bitmask.minimum_dtype())

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
        indx = numpy.array([numpy.invert(good_snr).reshape(self.spatial_shape).T]*self.nwave).T
        mask[indx] = self.bitmask.turn_on(mask[indx], flag='LOW_SNR')

        return mask


    def _check_snr(self):
        # binned_spectra['BINS'].data['SNR'] has length
        # binned_spectra.nbins
#        print(self.binned_spectra['BINS'].data['SNR'])
#        print(self.method['minimum_snr'])
        return self.binned_spectra['BINS'].data['SNR'] > self.method['minimum_snr']


    def _bins_to_fit(self):
        """Return flags for the bins to fit."""
        return (self._check_snr()) \
                    & ~(numpy.array([ b in self.binned_spectra.missing_bins \
                                        for b in numpy.arange(self.binned_spectra.nbins)]))


    def _assign_spectral_arrays(self):
        self.spectral_arrays = [ 'FLUX', 'MASK' ]


    def _assign_image_arrays(self):
        self.image_arrays = [ 'BINID' ]


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


    def all_except_emission_flags(self):
        return ['DIDNOTUSE', 'FORESTAR', 'LOW_SNR', 'ARTIFACT', 'OUTSIDE_RANGE', 'TPL_PIXELS',
                'TRUNCATED', 'PPXF_REJECT', 'FIT_FAILED' ]


    def fit(self, binned_spectra, guess_vel=None, guess_sig=None, dapsrc=None, dapver=None,
            analysis_path=None, directory_path=None, output_file=None, hardcopy=True,
            tpl_symlink_dir=None, clobber=False, loggers=None, quiet=False):
        """

        Fit the binned spectra

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
        
        self.nmodels = self.binned_spectra.nbins
        self.missing_models = self.binned_spectra.missing_bins

        # Get the guess kinematics
        if guess_vel is not None:
            self.guess_vel=guess_vel
        if self.guess_vel is None:
            raise ValueError('Must provide initial guess velocity.')
        self.guess_sig = 100.0 if guess_sig is None else guess_sig

        # Report
        good_bins = self._bins_to_fit()
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, '-'*50)
            log_output(self.loggers, 1, logging.INFO, 'STELLAR CONTINUUM FITTING:')
            log_output(self.loggers, 1, logging.INFO, '-'*50)
            log_output(self.loggers, 1, logging.INFO, 'Total bins: {0}'.format(
                                                            self.binned_spectra.nbins))
            log_output(self.loggers, 1, logging.INFO, 'Missing bins: {0}'.format(
                                                            len(self.binned_spectra.missing_bins)))
            log_output(self.loggers, 1, logging.INFO, 'With good S/N: {0}'.format(
                                                            numpy.sum(self._check_snr())))
            log_output(self.loggers, 1, logging.INFO, 'Total to fit: {0}'.format(
                                                            numpy.sum(good_bins)))

        if numpy.sum(good_bins) == 0:
            raise ValueError('No good spectra to fit!')

        # Fill in any remaining binning parameters
        self.tpl_symlink_dir = tpl_symlink_dir
        self._fill_method_par(dapsrc=dapsrc, analysis_path=analysis_path)

        # (Re)Set the output paths
        self._set_paths(directory_path, dapver, analysis_path, output_file)

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
        
        # If the file already exists, and not clobbering, just read the
        # file
        if os.path.isfile(ofile) and not clobber:
            if not self.quiet:
                log_output(self.loggers, 1, logging.INFO, 'Using existing file')
            self.read()
            if not self.quiet:
                log_output(self.loggers, 1, logging.INFO, '-'*50)
            return

        # Fit the spectra
        # Mask should be fully defined within the fitting function
        model_wave, model_flux, model_mask, model_par \
                = self.method['fitfunc'](self.binned_spectra, par=self.method['fitpar'],
                                         loggers=self.loggers, quiet=self.quiet)
        
        assert numpy.sum(model_wave - self.binned_spectra['WAVE'].data) == 0, \
                    'Incorrect wavelength range'

        unique_bins, reconstruct = numpy.unique(self.binned_spectra.hdu['BINID'].data.ravel(),
                                                return_inverse=True)

        # TODO: Use missing_models and missing_bins differently such
        # that the former actually reflects the list of missing models!
        # self.nmodels = numpy.amax(unique_bins)+1
        # self.missing_models = list(set(numpy.arange(self.nmodels)) - set(unique_bins))

        # Restructure the output to match the input DRP file
        indx = self.binned_spectra.hdu['BINID'].data.ravel() > -1

        _good_bins = numpy.full(numpy.prod(self.spatial_shape), False, dtype=numpy.bool)
        _good_bins[indx] = good_bins[unique_bins[reconstruct[indx]]]

        flux = numpy.zeros(self.shape, dtype=numpy.float)
        flux.reshape(-1,self.nwave)[indx,:] = model_flux[unique_bins[reconstruct[indx]],:]

        mask = self._initialize_mask(_good_bins)
        mask.reshape(-1,self.nwave)[indx,:] = model_mask[unique_bins[reconstruct[indx]],:]

        # Initialize the header keywords
        hdr = self._clean_drp_header(ext='PRIMARY')
        self._initialize_header(hdr)

        if self.method['fitclass'] is not None:
            try:
                hdr['SCTYPE'] = self.method['fitclass'].fit_type
                hdr['SCMETH'] = self.method['fitclass'].fit_method
            except:
                if not self.quiet and hardcopy:
                    warnings.warn('Fit class object does not have fit_type and/or fit_method ' \
                                  'attributes.  No parameters written to header.')
        if self.method['fitpar'] is not None:
            try:
                self.method['fitpar'].toheader(hdr)
            except:
                if not self.quiet and hardcopy:
                    warnings.warn('Fit parameter class has no toheader() function.  No ' \
                                  'parameters written to header.')


        # Save the data to the hdu attribute
        # TODO: Write the bitmask to the header?
        # TODO: Convert some/all of the binary table data into images
        self.hdu = fits.HDUList([ fits.PrimaryHDU(header=hdr),
                                  fits.ImageHDU(data=flux.data,
                                                header=self.binned_spectra['FLUX'].header.copy(),
                                                name='FLUX'),
                                  fits.ImageHDU(data=mask,
                                                header=self.binned_spectra['MASK'].header.copy(),
                                                name='MASK'),
                                  self.binned_spectra['WAVE'].copy(),
                                  self.binned_spectra['BINID'].copy(),
                                  fits.BinTableHDU.from_columns( [ fits.Column(name=n,
                                                             format=rec_to_fits_type(model_par[n]),
                                            array=model_par[n]) for n in model_par.dtype.names ],
                                                                    name='PAR')
                                ])

        # Write the data, if requested
        if not os.path.isdir(self.directory_path):
            os.makedirs(self.directory_path)
        self.hardcopy = hardcopy
        if self.hardcopy:
            self.write(clobber=clobber)
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, '-'*50)


    # Exact same function as used by SpatiallyBinnedSpectra
    def write(self, clobber=False):
        """
        Write the hdu object to the file.
        """

        # Restructure the data to match the DRPFits file
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, 'Restructuring data to match DRP.')
        if self.binned_spectra.drpf.mode == 'CUBE':
            MaNGAFits.restructure_cube(self.hdu, ext=self.spectral_arrays, inverse=True)
            MaNGAFits.restructure_map(self.hdu, ext=self.image_arrays, inverse=True)
        elif self.binned_spectra.drpf.mode == 'RSS':
            MaNGAFits.restructure_rss(self.hdu, ext=self.spectral_arrays, inverse=True)

        # Get the output file and determine if it should be compressed
        ofile = self.file_path()
        write_hdu(self.hdu, ofile, clobber=clobber, checksum=True, loggers=self.loggers,
                  quiet=self.quiet)

        # Revert the structure
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, 'Reverting to python-native structure.')
        if self.binned_spectra.drpf.mode == 'CUBE':
            MaNGAFits.restructure_cube(self.hdu, ext=self.spectral_arrays)
            MaNGAFits.restructure_map(self.hdu, ext=self.image_arrays)
        elif self.binned_spectra.drpf.mode == 'RSS':
            MaNGAFits.restructure_rss(self.hdu, ext=self.spectral_arrays)


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

        self.hdu = fits.open(ifile, checksum=checksum)

        # Confirm that the internal method is the same as the method
        # that was used in writing the file
        if self.hdu['PRIMARY'].header['SCKEY'] != self.method['key']:
            if strict:
                raise ValueError('Keywords in header does not match specified method keyword!')
            elif not self.quiet:
                warnings.warn('Keywords in header does not match specified method keyword!')
        # TODO: "strict" should also check other aspects of the file to
        # make sure that the details of the method are also the same,
        # not just the keyword

        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, 'Reverting to python-native structure.')
        if self.binned_spectra.drpf.mode == 'CUBE':
            MaNGAFits.restructure_cube(self.hdu, ext=self.spectral_arrays)
            MaNGAFits.restructure_map(self.hdu, ext=self.image_arrays)
        elif self.binned_spectra.drpf.mode == 'RSS':
            MaNGAFits.restructure_rss(self.hdu, ext=self.spectral_arrays)

        # Attempt to read the input guess velocity and sigma
        try:
            self.guess_vel = self.hdu['PRIMARY'].header['INPVEL']
            self.guess_sig = self.hdu['PRIMARY'].header['INPSIG']
        except:
            if not self.quiet:
                warnings.warn('Unable to read input guess kinematics from file header!')
            self.guess_vel = None
            self.guess_sig = None

        # Attempt to read the modeling parameters
        if self.method['fitpar'] is not None and callable(self.method['fitpar'].fromheader):
            self.method['fitpar'].fromheader(self.hdu['PRIMARY'].header)

        self.nmodels = self.hdu['PRIMARY'].header['NMOD']
        try:
            self.missing_models = eval(self.hdu['PRIMARY'].header['EMPTYMOD'])
        except KeyError:
            # Assume if this fails, it's because the keyword doesn't
            # exist
            self.missing_models = []


    def unique_models(self, index=False):
        """
        Get the unique models or the indices of the unique models in the
        flattened spatial dimension.
        """
        unique_models, first_occurrence = numpy.unique(self.hdu['BINID'].data.ravel(),
                                                       return_index=True)
        return first_occurrence[1:] if index else unique_models[1:]


    def copy_to_array(self, waverange=None, ext='FLUX'):
        r"""
        Return a copy of the selected data array.  The array size is
        always :math:`N_{\rm models} \times N_{\rm wavelength}`; i.e., the
        data is always flattened to two dimensions and the unique
        spectra are selected.  See :func:`unique_models` and
        :func:`mangadap.drpfits.DRPFits.copy_to_array`.

        Args:
            waverange (array-like) : (**Optional**) Two-element array
                with the first and last wavelength to include in the
                computation.  Default is to use the full wavelength
                range.
            ext (str) : (**Optional**) Name of the extension from which
                to draw the data.  Must be allowed for the current
                :attr:`mode`; see :attr:`data_arrays`.  Default is
                ``'FLUX'``.

        Returns:
            numpy.ndarray: A 2D array with a copy of the data from the
            selected extension.
        """
        arr = DRPFits.copy_to_array(self, waverange=waverange, ext=ext)[
                                                                self.unique_models(index=True),:]
        if len(self.missing_models) == 0:
            return arr
        _arr = numpy.zeros( (self.nmodels, self.nwave), dtype=self.hdu[ext].data.dtype)
        _arr[self.unique_models(),:] = arr
        return _arr


    def copy_to_masked_array(self, waverange=None, ext='FLUX', mask_ext='MASK', flag=None):
        r"""

        Return a copy of the selected data array as a masked array.
        This is functionally identical to :func:`copy_to_array`,
        except the output format is a `numpy.ma.MaskedArray`_.  The
        pixels that are considered to be masked can be specified using
        `flag`.
        
        Args:
            waverange (array-like) : (**Optional**) Two-element array
                with the first and last wavelength to include in the
                computation.  Default is to use the full wavelength
                range.
            ext (str) : (**Optional**) Name of the extension from which
                to draw the data.  Must be allowed for the current
                :attr:`mode`; see :attr:`data_arrays`.  Default is
                `'FLUX'`.
            mask_ext (str) : (**Optional**) Name of the extension with
                the mask bit data.  Must be allowed for the current
                :attr:`mode`; see :attr:`data_arrays`.  Default is
                `'MASK'`.
            flag (str or list): (**Optional**) (List of) Flag names that
                are considered when deciding if a pixel should be
                masked.  The names *must* be a valid bit name as defined
                by :attr:`bitmask` (see :class:`DRPFitsBitMask`).  If
                not provided, *ANY* non-zero mask bit is omitted.

        Returns:
            numpy.ndarray: A 2D array with a copy of the data from the
            selected extension.
        """
        arr = DRPFits.copy_to_masked_array(self, waverange=waverange, ext=ext, mask_ext=mask_ext,
                                           flag=flag)[self.unique_models(index=True),:]
        if len(self.missing_models) == 0:
            return arr
        _arr = numpy.ma.zeros( (self.nmodels, self.nwave), dtype=self.hdu[ext].data.dtype)
        _arr[self.unique_models(),:] = arr
        _arr[self.missing_models,:] = numpy.ma.masked
        return _arr


    def construct_models(self, template_library=None, redshift_only=False):
        if self.method['fitclass'] is None:
            raise ValueError('No class object available for constructing the model!')
        if not callable(self.method['fitclass'].construct_models):
            raise AttributeError('Provided fit class object has no callable ' \
                                 '\'construct_models\' attribute!')

        _template_library = self.method['fitpar']['template_library'] \
                    if template_library is None else template_library
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO,
                            'Constructing model spectra {0} dispersion'.format(
                                            'not including' if redshift_only else 'including'))

        return self.method['fitclass'].construct_models(self.binned_spectra, _template_library,
                                                        self.hdu['PAR'].data,
                                                        self.method['fitpar']['degree'],
                                                        self.method['fitpar']['mdegree'],
                                                        self.method['fitpar']['moments'],
                                                        redshift_only=redshift_only,
                                        velscale_ratio=self.method['fitpar']['velscale_ratio'],
                                                        dvtol=1e-9)


