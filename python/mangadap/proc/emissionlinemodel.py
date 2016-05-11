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


*Class usage examples*:
    Add examples!

*Revision history*:
    | **26 Apr 2016**: Implementation begun by K. Westfall (KBW)

.. _astropy.io.fits.hdu.hdulist.HDUList: http://docs.astropy.org/en/v1.0.2/io/fits/api/hdulists.html
.. _glob.glob: https://docs.python.org/3.4/library/glob.html
.. _configparser.ConfigParser: https://docs.python.org/3/library/configparser.html#configparser.ConfigParser


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

import glob
import os.path
import numpy
from astropy.io import fits
import astropy.constants

from ..drpfits import DRPFits
from ..par.parset import ParSet
from ..config.defaults import default_dap_source, default_dap_reference_path
from ..config.defaults import default_dap_file_name
from ..util.fileio import init_record_array, rec_to_fits_type, rec_to_fits_col_dim, write_hdu
from ..util.bitmask import BitMask
from .artifactdb import ArtifactDB
from .emissionlinedb import EmissionLineDB
from .pixelmask import SpectralPixelMask
from .spatiallybinnedspectra import SpatiallyBinnedSpectra
from .stellarcontinuummodel import StellarContinuumModel
from .lineprofilefit import Elric, ElricPar, GaussianLineProfile
from .util import _select_proc_method

from matplotlib import pyplot

__author__ = 'Kyle B. Westfall'
# Add strict versioning
# from distutils.version import StrictVersion


class EmissionLineModelDef(ParSet):
    """
    A class that holds the parameters necessary to perform the
    emission-line profile fits.

    Args:
        key (str): Keyword used to distinguish between different
            emission-line moment databases.
        minimum_snr (bool): Minimum S/N of spectrum to fit
        artifacts (str): String identifying the artifact database to use
        emission_lines (str): String identifying the emission-line
            database to use
    """
    def __init__(self, key, minimum_snr, artifacts, emission_lines, fitpar, fitclass, fitfunc):
        in_fl = [ int, float ]
        par_opt = [ ParSet, dict ]

        pars =     [ 'key', 'minimum_snr', 'artifacts', 'emission_lines', 'fitpar', 'fitclass',
                     'fitfunc' ]
        values =   [   key,   minimum_snr,   artifacts,   emission_lines,   fitpar,   fitclass,
                       fitfunc ]
        dtypes =   [   str,         in_fl,         str,              str,  par_opt,       None,
                          None ]

        can_call = [ False,         False,      False,    False,      False,      True ]


        ParSet.__init__(self, pars, values=values, dtypes=dtypes)


def validate_emission_line_modeling_method_config(cnfg):
    """ 
    Validate the `configparser.ConfigParser`_ object that is meant to
    define an emission-line modeling method.

    Args:
        cnfg (`configparser.ConfigParser`_): Object meant to contain
            defining parameters of the emission-line modeling method
            needed by :class:`EmissionLineModelDef'

    Raises:
        KeyError: Raised if any required keywords do not exist.
        ValueError: Raised if keys have unacceptable values.
    """
    # Check for required keywords
    if 'key' not in cnfg.options('default'):
        raise KeyError('Keyword \'key\' must be provided.')
    if 'minimum_snr' not in cnfg.options('default') or cnfg['default']['minimum_snr'] is None:
        cnfg['default']['minimum_snr']= '0.0'

    if 'artifact_mask' not in cnfg.options('default') \
            or cnfg['default']['artifact_mask'] is None:
        cnfg['default']['artifact_mask'] = 'None'

    if 'emission_lines' not in cnfg.options('default') \
            or cnfg['default']['emission_lines'] is None:
        raise ValueError('Must provide a keyword with the emission-line parameters to use!')

#    if 'profile_type' not in cnfg.options('default') or cnfg['default']['emission_lines'] is None:
#        warnings.warn('No profile provided.  Assuming Gaussian.')
#        cnfg['default']['profile_type'] = 'GaussianLineProfile'
#
#    if 'number_of_components' not in cnfg.options('default') \
#            or cnfg['default']['number_of_components'] is None:
#        cnfg['default']['number_of_components'] = '1'

    if 'baseline_order' not in cnfg.options('default') \
            or cnfg['default']['baseline_order'] is None:
        cnfg['default']['baseline_order'] = '-1'

    if 'window_buffer' not in cnfg.options('default') \
            or cnfg['default']['window_buffer'] is None:
        cnfg['default']['window_buffer'] = '25'


def available_emission_line_modeling_methods(dapsrc=None):
    """
    Return the list of available emission-line modeling methods.

    Available methods are:

    .. todo::
        Fill in

    Args:
        dapsrc (str): (**Optional**) Root path to the DAP source
            directory.  If not provided, the default is defined by
            :func:`mangadap.config.defaults.default_dap_source`.

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
    search_dir = os.path.join(dapsrc, 'python/mangadap/config/emission_line_modeling')
    ini_files = glob.glob(os.path.join(search_dir, '*.ini'))
    if len(ini_files) == 0:
        raise IOError('Could not find any configuration files in {0} !'.format(search_dir))

    # Build the list of library definitions
    method_list = []
    for f in ini_files:
        # Read the config file
        cnfg = ConfigParser(allow_no_value=True)
        cnfg.read(f)
        # Ensure it has the necessary elements to define the template
        # library
        validate_emission_line_modeling_method_config(cnfg)

        # Currently only implement Elric for fitting the emission lines
        fitpar = ElricPar(None, cnfg['default'].getint('baseline_order'),
                          cnfg['default'].getfloat('window_buffer'), None, None,
                          minimum_snr=cnfg['default'].getfloat('minimum_snr'))
        fitclass = Elric(EmissionLineModelBitMask(dapsrc=dapsrc))
        fitfunc = fitclass.fit_SpatiallyBinnedSpectra

        method_list += [ EmissionLineModelDef(cnfg['default']['key'],
                                              eval(cnfg['default']['minimum_snr']),
                                              None if cnfg['default']['artifact_mask'] == 'None' \
                                                    else cnfg['default']['artifact_mask'],
                                              cnfg['default']['emission_lines'], fitpar, fitclass,
                                              fitfunc) ]

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
        dapsrc = default_dap_source() if dapsrc is None else str(dapsrc)
        BitMask.__init__(self, ini_file=os.path.join(dapsrc, 'python', 'mangadap', 'config',
                                                     'bitmasks', 'emission_line_model_bits.ini'))


class EmissionLineModel:
    r"""

    Class that holds the emission-line model results.

    """
    def __init__(self, method_key, binned_spectra, guess_vel=None, guess_sig=None,
                 stellar_continuum=None, method_list=None, artifact_list=None,
                 emission_line_db_list=None, dapsrc=None, dapver=None, analysis_path=None,
                 directory_path=None, output_file=None, hardcopy=True, clobber=False, verbose=0,
                 checksum=False):

        self.version = '1.0'
        self.verbose = verbose

        # Define the database properties
        self.method = None
        self.artdb = None
        self.emldb = None
        self.pixelmask = None
        self._define_method(method_key, method_list=method_list, artifact_list=artifact_list,
                            emission_line_db_list=emission_line_db_list)

        self.binned_spectra = None
        self.redshift = None
        self.dispersion = None
        self.stellar_continuum = None
        self.neml = self.emldb.nsets

        # Define the output directory and file
        self.directory_path = None      # Set in _set_paths
        self.output_file = None
        self.hardcopy = None

        # Initialize the objects used in the assessments
        self.bitmask = EmissionLineModelBitMask(dapsrc=dapsrc)
        self.hdu = None
        self.checksum = checksum

        self.spatial_shape = None
        self.spatial_index = None

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

        # Fit the binned spectra
        self.fit(binned_spectra, guess_vel=guess_vel, guess_sig=guess_sig,
                 stellar_continuum=stellar_continuum, dapsrc=dapsrc, dapver=dapver,
                 analysis_path=analysis_path, directory_path=directory_path,
                 output_file=output_file, hardcopy=hardcopy, clobber=clobber, verbose=verbose)


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


    def _define_method(self, method_key, method_list=None, artifact_list=None,
                       emission_line_db_list=None, dapsrc=None):
        r"""
        Select the modeling method
        """
        # Grab the specific method
        self.method = _select_proc_method(method_key, EmissionLineModelDef, method_list=method_list,
                                          available_func=available_emission_line_modeling_methods,
                                          dapsrc=dapsrc)

        # Instantiate the artifact and emission-line databases
        self.artdb = None if self.method['artifacts'] is None else \
                        ArtifactDB(self.method['artifacts'],artdb_list=artifact_list,dapsrc=dapsrc)
        # TODO: Generalize the name of this object
        self.pixelmask = SpectralPixelMask(artdb=self.artdb)

        self.emldb = None if self.method['emission_lines'] is None else \
                        EmissionLineDB(self.method['emission_lines'],
                                       emldb_list=emission_line_db_list, dapsrc=dapsrc)


    def _set_paths(self, directory_path, dapver, analysis_path, output_file):
        """
        Set the :attr:`directory_path` and :attr:`output_file`.  If not
        provided, the defaults are set using, respectively,
        :func:`mangadap.config.defaults.default_dap_reference_path` and
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
        self.directory_path = default_dap_reference_path(plate=self.binned_spectra.drpf.plate,
                                                         drpver=self.binned_spectra.drpf.drpver,
                                                         dapver=dapver,
                                                         analysis_path=analysis_path) \
                                        if directory_path is None else str(directory_path)

        # Set the output file
        method = '{0}-{1}'.format(self.binned_spectra.rdxqa.method['key'],
                                  self.binned_spectra.method['key'])
        if self.stellar_continuum is not None:
            method = '{0}-{1}'.format(method, self.stellar_continuum.method['key'])
        method = '{0}-{1}'.format(method, self.method['key'])

        self.output_file = default_dap_file_name(self.binned_spectra.drpf.plate,
                                                 self.binned_spectra.drpf.ifudesign,
                                                 self.binned_spectra.drpf.mode, method) \
                                        if output_file is None else str(output_file)


    def _emission_line_database_dtype(self, name_len, mode_len, prof_len):
        r"""
        Construct the record array data type for the output fits
        extension.
        """
        return [ ('ID',numpy.int),
                 ('NAME','<U{0:d}'.format(name_len)),
                 ('RESTWAVE', numpy.float),
                 ('ACTION', '<U1'),
                 ('FLUXRATIO', numpy.float),
                 ('MODE','<U{0:d}'.format(mode_len)),
                 ('PROFILE', '<U{0:d}'.format(prof_len)),
                 ('NCOMP', numpy.int)
               ]

    
    def _compile_database(self):
        """
        Compile the database with the specifications of each index.
        """
        name_len = numpy.amax([ len(n) for n in self.emldb['name'] ])
        mode_len = numpy.amax([ len(m) for m in self.emldb['mode'] ])
        prof_len = numpy.amax([ len(p) for p in self.emldb['profile'] ])

        # Instatiate the table data that will be saved defining the set
        # of emission-line moments measured
        hdu_database = init_record_array(self.neml,
                                self._emission_line_database_dtype(name_len, mode_len, prof_len))

        hk = [ 'ID', 'NAME', 'RESTWAVE', 'ACTION', 'FLUXRATIO', 'MODE', 'PROFILE', 'NCOMP' ]
        mk = [ 'index', 'name', 'restwave', 'action', 'flux', 'mode', 'profile', 'ncomp' ]
        for _hk, _mk in zip(hk,mk):
            hdu_database[_hk] = self.emldb[_mk]

        #print(hdu_database)

        return hdu_database


    def _assign_input_kinematics(self, redshift, dispersion, default_dispersion=100.0):
        """
        Set the initial redshift and velocity dispersion for each
        spectrum.  Directly provided redshifts take precedence over
        those in the StellarContinuumModel object, if both are provided.
        The dispersions are **never** drawn from the stellar kinematics,
        but set to default_dispersion if they're not provided.
        """

        self.redshift = numpy.zeros(self.nmodels, dtype=numpy.float)
        if redshift is not None:
            _redshift = numpy.atleast_1d(redshift)
            if len(_redshift) not in [ 1, self.nmodels, self.nspec ]:
                raise ValueError('Provided redshift must be either a single value, match the ' \
                                 'number of model spectra, or match the total number of DRP ' \
                                 'spectra.')
            if len(_redshift) == 1:
                self.redshift[:] = redshift
            elif len(_redshift) == self.nspec:
                self.redshift = _redshift[ self.unique_models(index=True) ]
            else:   # Has length nmodels
                self.redshift = _redshift
        elif self.stellar_continuum is not None:
            self.redshift = self.stellar_continuum['PAR'].data['KIN'][:,0] \
                                / astropy.constants.c.to('km/s').value

        self.dispersion = numpy.full(self.nmodels, default_dispersion, dtype=numpy.float)
        if dispersion is not None:
            _dispersion = numpy.atleast_1d(dispersion)
            if len(_dispersion) not in [ 1, self.nmodels, self.nspec ]:
                raise ValueError('Provided dispersion must be either a single value, match the ' \
                                 'number of model spectra, or match the total number of DRP ' \
                                 'spectra.')
            if len(_disperison) == 1:
                self.dispersion[:] = dispersion
            elif len(_redshift) == self.nspec:
                self.dispersion = _dispersion[ self.unique_models(index=True) ]
            else:   # Has length nmodels
                self.dispersion = _dispersion


    def _fill_method_par(self, dapsrc=None, analysis_path=None):
        """

        Fill in any remaining modeling parameters.

        """
        # Fill the guess kinematics
        self.method['fitpar']['guess_redshift'] = self.redshift
        self.method['fitpar']['guess_dispersion'] = self.dispersion
        self.method['fitpar']['emission_lines'] = self.emldb
        self.method['fitpar']['pixelmask'] = self.pixelmask
        if self.stellar_continuum is not None:
            self.method['fitpar']['stellar_continuum'] = self.stellar_continuum


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
        hdr['ELFKEY'] = (self.method['key'], 'Emission-line modeling method keyword')
        hdr['ELFMINSN'] = (self.method['minimum_snr'], 'Minimum S/N of spectrum to include')
        hdr['ARTDB'] = (self.method['artifacts'], 'Artifact database keyword')
        hdr['EMLDB'] = (self.method['emission_lines'], 'Emission-line database keyword')
        if self.stellar_continuum is not None:
            hdr['SCTPL'] = (self.stellar_continuum.method['key'], 'Stellar-continuum model keyword')
        hdr['NMOD'] = (self.nmodels, 'Number of unique emission-line models')
        if len(self.missing_models) > 0:
            hdr['EMPTYMOD'] = (str(self.missing_models), 'List of models with no data')
        # Anything else?
        # Additional database details?
        return hdr


    def _initialize_mask(self):
        """

        Initialize the mask be setting the DIDNOTUSE, FORESTAR, and LOW_SNR masks

        """
        # Initialize to all zeros
        mask = numpy.zeros(self.shape, dtype=self.bitmask.minimum_uint_dtype())

        # Turn on the flag stating that the pixel wasn't used
        indx = self.binned_spectra.bitmask.flagged(self.binned_spectra['MASK'].data,
                                                   flag=self.binned_spectra.do_not_fit_flags())
        mask[indx] = self.bitmask.turn_on(mask[indx], 'DIDNOTUSE')

        # Turn on the flag stating that the pixel has a foreground star
        indx = self.binned_spectra.bitmask.flagged(self.binned_spectra['MASK'].data,
                                                   flag='FORESTAR')
        mask[indx] = self.bitmask.turn_on(mask[indx], 'FORESTAR')

        return mask


    def _check_snr(self):
        # binned_spectra['BINS'].data['SNR'] has length
        # binned_spectra.nbins
        return self.binned_spectra['BINS'].data['SNR'] > self.method['minimum_snr']


    def _assign_spectral_arrays(self):
        self.spectral_arrays = [ 'FLUX', 'MASK' ]


    def _assign_image_arrays(self):
        self.image_arrays = [ 'BINID' ]


    def _missing_flags(self):
        return numpy.array([ b in self.missing_models for b in numpy.arange(self.nmodels)])


    def file_name(self):
        """Return the name of the output file."""
        return self.output_file


    def file_path(self):
        """Return the full path to the output file."""
        if self.directory_path is None or self.output_file is None:
            return None
        return os.path.join(self.directory_path, self.output_file)

        self.fit(binned_spectra, guess_vel=guess_vel, guess_sig=guess_sig,
                 stellar_continuum=stellar_continuum, dapsrc=dapsrc, dapver=dapver,
                 analysis_path=analysis_path, directory_path=directory_path,
                 output_file=output_file, hardcopy=hardcopy, clobber=clobber, verbose=verbose)

    def fit(self, binned_spectra, guess_vel=None, guess_sig=None, stellar_continuum=None,
            dapsrc=None, dapver=None, analysis_path=None, directory_path=None, output_file=None,
            hardcopy=True, clobber=False, verbose=0):
        """

        Fit the emission lines.

        """

        # SpatiallyBinnedSpectra object always needed
        if binned_spectra is None:
            raise ValueError('Must provide spectra object for fitting.')
        if not isinstance(binned_spectra, SpatiallyBinnedSpectra):
            raise TypeError('Must provide a valid SpatiallyBinnedSpectra object!')
        if binned_spectra.hdu is None:
            raise ValueError('Provided SpatiallyBinnedSpectra object is undefined!')
        self.binned_spectra = binned_spectra

        # StellarContinuumModel object only used when calculating dispersion corrections.
        if stellar_continuum is not None:
            if not isinstance(stellar_continuum, StellarContinuumModel):
                raise TypeError('Provided stellar continuum must have StellarContinuumModel type!')
            if stellar_continuum.hdu is None:
                raise ValueError('Provided StellarContinuumModel is undefined!')
            self.stellar_continuum = stellar_continuum

        self.shape = self.binned_spectra.shape
        self.spatial_shape =self.binned_spectra.spatial_shape
        self.nspec = self.binned_spectra.nspec
        self.spatial_index = self.binned_spectra.spatial_index.copy()
        self.dispaxis = self.binned_spectra.dispaxis
        self.nwave = self.binned_spectra.nwave
        
        self.nmodels = self.binned_spectra.nbins
        self.missing_models = self.binned_spectra.missing_bins

        # Get the guess kinematics
        redshift = None if guess_vel is None else guess_vel / astropy.constants.c.to('km/s').value
        dispersion = None if guess_sig is None else guess_sig
        self._assign_input_kinematics(redshift, dispersion)

        # Get the good spectra
        #   - Must have sufficienct S/N, as defined by the input par
        good_snr = self._check_snr()
        # Report
        print('Total spectra: ', len(good_snr))
        print('With good S/N: ', numpy.sum(good_snr))

        # Fill in any remaining binning parameters
        self._fill_method_par(dapsrc=dapsrc, analysis_path=analysis_path)

        # (Re)Set the output paths
        self._set_paths(directory_path, dapver, analysis_path, output_file)

        # Check that the file path is defined
        ofile = self.file_path()
        if ofile is None:
            raise ValueError('File path for output file is undefined!')

        # Report
        print('Output path: ', self.directory_path)
        print('Output file: ', self.output_file)

        # If the file already exists, and not clobbering, just read the
        # file
        if os.path.isfile(ofile) and not clobber:
            self.read()
            return

        # Fit the spectra
        # Mask should be fully defined within the fitting function
        model_wave, model_flux, model_mask, model_fit_par, model_eml_par = \
            self.method['fitfunc'](self.binned_spectra, par=self.method['fitpar'])

#        pyplot.step(model_wave, model_flux[0,:], where='mid')
#        pyplot.show()
#        exit()
        
        # Compile the information on the suite of measured indices
        hdu_database = self._compile_database()

        # Get the unique bins and how to reconstruct them
        unique_models, reconstruct = numpy.unique(self.binned_spectra.hdu['BINID'].data.ravel(),
                                                return_inverse=True)

        # Restructure the output to match the input DRP file
        indx = self.binned_spectra.hdu['BINID'].data.ravel() > -1

        flux = numpy.zeros(self.shape, dtype=numpy.float)
        flux.reshape(-1,self.nwave)[indx,:] = model_flux[unique_models[reconstruct[indx]],:]

        mask = self._initialize_mask()
        mask.reshape(-1,self.nwave)[indx,:] = model_mask[unique_models[reconstruct[indx]],:]

        # Initialize the header keywords
        hdr = self._clean_drp_header(ext='PRIMARY')
        self._initialize_header(hdr)

        if self.method['fitclass'] is not None:
            try:
                hdr['EMLTYPE'] = self.method['fitclass'].fit_type
                hdr['EMLMETH'] = self.method['fitclass'].fit_method
            except:
                if hardcopy:
                    warnings.warn('Fit class object does not have fit_type and/or fit_method ' \
                                  'attributes.  No parameters written to header.')
        if self.method['fitpar'] is not None:
            try:
                self.method['fitpar'].toheader(hdr)
            except:
                if hardcopy:
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
                                                        format=rec_to_fits_type(hdu_database[n]),
                                        array=hdu_database[n]) for n in hdu_database.dtype.names ],
                                                         name='PAR'),
                                  fits.BinTableHDU.from_columns([fits.Column(name=n,
                                                        format=rec_to_fits_type(model_fit_par[n]),
                                                        dim=rec_to_fits_col_dim(model_fit_par[n]),
                                        array=model_fit_par[n]) for n in model_fit_par.dtype.names],
                                                                    name='FIT'),
                                  fits.BinTableHDU.from_columns([fits.Column(name=n,
                                                        format=rec_to_fits_type(model_eml_par[n]),
                                                        dim=rec_to_fits_col_dim(model_eml_par[n]),
                                        array=model_eml_par[n]) for n in model_eml_par.dtype.names],
                                                                    name='EMLDATA')
                                ])

        # Write the data, if requested
        if not os.path.isdir(self.directory_path):
            os.makedirs(self.directory_path)
        self.hardcopy = hardcopy
        if self.hardcopy:
            self.write(clobber=clobber)


    def write(self, clobber=False):
        """
        Write the hdu object to the file.
        """
        # Restructure the data to match the DRPFits file
        if self.binned_spectra.drpf.mode == 'CUBE':
            DRPFits._restructure_cube(self.hdu, ext=self.spectral_arrays, inverse=True)
            SpatiallyBinnedSpectra._restructure_map(self.hdu, ext=self.image_arrays, inverse=True)
        elif self.binned_spectra.drpf.mode == 'RSS':
            DRPFits._restructure_rss(self.hdu, ext=self.spectral_arrays, inverse=True)
        
        # Get the output file and determine if it should be compressed
        ofile = self.file_path()
        write_hdu(self.hdu, ofile, clobber=clobber, checksum=True)

        # Revert the structure
        if self.binned_spectra.drpf.mode == 'CUBE':
            DRPFits._restructure_cube(self.hdu, ext=self.spectral_arrays)
            SpatiallyBinnedSpectra._restructure_map(self.hdu, ext=self.image_arrays)
        elif self.binned_spectra.drpf.mode == 'RSS':
            DRPFits._restructure_rss(self.hdu, ext=self.spectral_arrays)


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

        self.hdu = fits.open(ifile, checksum=checksum)

        # Confirm that the internal method is the same as the method
        # that was used in writing the file
        if self.hdu['PRIMARY'].header['ELFKEY'] != self.method['key']:
            if strict:
                raise ValueError('ELFKEY in header does not match specified method keyword!')
            else:
                warnings.warn('ELFKEY in header does not match specified method keyword!')
        # TODO: "strict" should also check other aspects of the file to
        # make sure that the details of the method are also the same,
        # not just the keyword

        if self.binned_spectra.drpf.mode == 'CUBE':
            DRPFits._restructure_cube(self.hdu, ext=self.spectral_arrays)
            SpatiallyBinnedSpectra._restructure_map(self.hdu, ext=self.image_arrays)
        elif self.binned_spectra.drpf.mode == 'RSS':
            DRPFits._restructure_rss(self.hdu, ext=self.spectral_arrays)

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




