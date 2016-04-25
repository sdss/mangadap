# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""

A class hierarchy that perform that spatially bins a set of
two-dimensional data.

The base class allows for user-defined definitions of binning
procedures.

*License*:
    Copyright (c) 2015, SDSS-IV/MaNGA Pipeline Group
        Licensed under BSD 3-clause license - see LICENSE.rst

*Source location*:
    $MANGADAP_DIR/python/mangadap/proc/spatialbins.py

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
        
        import glob
        import os.path
        from os import remove, environ
        from scipy import sparse
        from astropy.io import fits
        import astropy.constants
        import time
        import numpy

        from ..par.parset import ParSet
        from ..config.defaults import default_dap_source, default_dap_reference_path
        from ..config.defaults import default_dap_file_name
        from ..util.idlutils import airtovac
        from ..util.geometry import SemiMajorAxisCoo
        from ..util.fileio import init_record_array
        from ..drpfits import DRPFits

*Class usage examples*:

    .. todo::
        Add examples

*Revision history*:
    | **01 Apr 2016**: Implementation begun by K. Westfall (KBW)

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
from os import remove, environ
from scipy import sparse
from astropy.io import fits
import astropy.constants
import time
import numpy

from . import spatialbinning
from .reductionassessments import ReductionAssessment
from .spectralstack import SpectralStackPar, SpectralStack
from ..par.parset import ParSet
from ..config.defaults import default_dap_source, default_dap_reference_path
from ..config.defaults import default_dap_file_name, default_cube_pixelscale
from ..util.fileio import init_record_array, rec_to_fits_type, write_hdu
from ..util.bitmask import BitMask
from ..util.covariance import Covariance
from ..util.geometry import SemiMajorAxisCoo
from ..drpfits import DRPFits
from .util import _select_proc_method, _fill_vector

from matplotlib import pyplot

__author__ = 'Kyle B. Westfall'
# Add strict versioning
# from distutils.version import StrictVersion


class SpatiallyBinnedSpectraDef(ParSet):
    """
    A class that holds the two parameter sets and the key designator for
    the binning scheme.

    binfunc must be able to function properly with the following call:
        
        id = bin(x, y, par=par)

    where ``x`` and ``y`` are the on-sky Cartesian coordiantes provided
    by::
        
        x = ReductionAssessment.hdu['SPECTRUM'].data['SKY_COO'][:,0]
        y = ReductionAssessment.hdu['SPECTRUM'].data['SKY_COO'][:,1]

    id is the index of the bin to place each spectrum.

    stackfunc must be able to function properly with the following
    call::
        
        stack_flux, stack_sdev, stack_npix, stack_ivar, stack_covar = \
                stack(drpf, id, par=par)

    where `drpf` is a :class:`mangadap.drpfits.DRPFits` object.  Note
    that the input and output wavelength ranges are expected to be the
    same!  This is why the wavelengths are not returned!

    As long as they are mutable, the values in par can change, meaning
    that some products of teh bin algorithm can be passed to the stack
    algorithm.  For example, if you want to weight the inclusion of the
    spectrum in the bin, you'd have to provide both the binning and
    stacking routines.  Actually, if that's the case, you're better off
    providing a class object that will both bin and stack the spectra!

    Args:
        key (str): Keyword used to distinguish between different spatial
            binning schemes.
        binpar (:class:`mangadap.par.parset.ParSet` or dict): The
            parameter set defining how to place each spectrum in a bin.
        binclass (object): Instance of class object to use for the
            binning.  Needed in case binfunc is a non-static member
            function of the class.
        binfunc (callable): The function that determines which spectra
            go into each bin.
        stackpar (:class:`mangadap.par.parset.ParSet` or dict): The parameter
            set defining how to stack the spectra in each bin.
        stackclass (object): Instance of class object to used to stack
            the spectra.  Needed in case stackfunc is a non-static
            member function of the class.
        stackfunc (callable): The function that stacks the spectra in a
            given bin.
    """
    def __init__(self, key, minimum_snr, binpar, binclass, binfunc, stackpar, stackclass,
                 stackfunc):
        in_fl = [ int, float ]
#        bincls_opt = [ spatialbinning.SpatialBinning ]
#        stackcls_opt = [ SpectralStack ]
        par_opt = [ ParSet, dict ]

        pars =     [ 'key', 'minimum_snr', 'binpar', 'binclass', 'binfunc', 'stackpar',
                        'stackclass', 'stackfunc' ]
        values =   [   key,   minimum_snr,   binpar,   binclass,   binfunc,   stackpar,
                          stackclass,   stackfunc ]
        dtypes =   [   str,         in_fl,  par_opt,       None,      None,    par_opt,
                                None,        None ]
        can_call = [ False,         False,    False,      False,      True,      False,
                               False,        True ]

        ParSet.__init__(self, pars, values=values, dtypes=dtypes, can_call=can_call)


def validate_spatial_binning_scheme_config(cnfg):
    """ 
    Validate the `configparser.ConfigParser`_ object that is meant to
    define a spatial binning scheme.

    Args:
        cnfg (`configparser.ConfigParser`_): Object meant to contain
            defining parameters of the binning method as needed by
            :class:`mangadap.proc.spatiallybinnedspectra.SpatiallyBinnedSpectraDef`.

    Raises:
        KeyError: Raised if any required keywords do not exist.
        ValueError: Raised if keys have unacceptable values.
        FileNotFoundError: Raised if a file is specified but could not
            be found.
    """
    # Check for required keywords
    if 'key' not in cnfg.options('default'):
        raise KeyError('Keyword \'key\' must be provided.')
    if 'method' not in cnfg.options('default'):
        raise KeyError('Keyword \'method\' must be provided.')

    if 'minimum_snr' not in cnfg.options('default') or cnfg['default']['minimum_snr'] is None:
        cnfg['default']['minimum_snr']= '0.0'

    if 'operation' not in cnfg.options('default') or cnfg['default']['operation'] is None:
        cnfg['default']['operation'] = 'mean'

    if 'velocity_register' not in cnfg.options('default') \
            or cnfg['default']['velocity_register'] is None:
        cnfg['default']['velocity_register'] = 'False'

    if 'stack_covariance_mode' not in cnfg.options('default') \
            or cnfg['default']['stack_covariance_mode'] is None:
        cnfg['default']['stack_covariance_mode'] = 'none'

    covar_par_needed_modes = SpectralStack.covariance_mode_options(par_needed=True)
    if cnfg['default']['stack_covariance_mode'] in covar_par_needed_modes \
            and ('stack_covariance_par' not in cnfg.options('default') \
                    or cnfg['default']['stack_covariance_par'] is None):
        raise ValueError('For covariance mode = {0}, must provide a parameter!'.format(
                         cnfg['default']['stack_covariance_mode']))
        
    if 'stack_covariance_par' not in cnfg.options('default') \
            or cnfg['default']['stack_covariance_par'] is None:
        cnfg['default']['stack_covariance_par'] = 'none'

    if cnfg['default']['method'] in [ 'none', 'global' ]:
        return

    if cnfg['default']['method'] == 'voronoi':
        if 'target_snr' not in cnfg.options('default'):
            raise KeyError('Keyword \'target_snr\' must be provided for Voronoi binning.')
        return

    if cnfg['default']['method'] == 'radial':
        if 'center' not in cnfg.options('default'):
            raise KeyError('Keyword \'center\' must be provided for radial binning.')
        if 'pa' not in cnfg.options('default'):
            raise KeyError('Keyword \'pa\' must be provided for radial binning.')
        if 'ell' not in cnfg.options('default'):
            raise KeyError('Keyword \'ell\' must be provided for radial binning.')
        if 'radii' not in cnfg.options('default'):
            raise KeyError('Keyword \'radii\' must be provided for radial binning.')

        if 'radius_scale' not in cnfg.options('default') or cnfg['default']['radius_scale'] is None:
            cnfg['default']['radius_scale'] = '1.0'
        if 'log_step' not in cnfg.options('default') or cnfg['default']['log_step'] is None:
            cnfg['default']['log_step'] = 'False'
        return

    raise ValueError('{0} is not a recognized binning method.'.format(cnfg['default']['method']))


def available_spatial_binning_methods(dapsrc=None):
    """
    Return the list of available binning schemes.  The following methods
    are predefined by the DAP configuration files.

    Voronoi binning methods:
    +---------+-----+------+-------+-----------+-------+------+
    |         | Min |      |  Vel. |     Covar | Covar | Targ |
    |     Key | S/N | Stat |  Reg. |      Mode |   par |  S/N |
    +=========+=====+======+=======+===========+=======+======+
    |   VOR05 | 0.0 | Mean | False | calibrate |  1.62 |    5 |
    +---------+-----+------+-------+-----------+-------+------+
    |   VOR10 | 0.0 | Mean | False | calibrate |  1.62 |   10 |
    +---------+-----+------+-------+-----------+-------+------+
    |   VOR30 | 0.0 | Mean | False | calibrate |  1.62 |   30 |
    +---------+-----+------+-------+-----------+-------+------+

    Radial binning methods:
    +----------+-----+------+-------+-----------+-------+--------+----+-----+-------+----------+-------+
    |          | Min |      | Vel.  |     Covar | Covar |        |    |     | R     |          |       |
    |      Key | S/N | Stat | Reg.  |      Mode |   par | Center | PA | Ell | Scale |    Radii |   Log |
    +==========+=====+======+=======+===========+=======+========+====+=====+=======+==========+=======+
    |   LOGR10 | 0.0 | Mean | False | calibrate |  1.62 |    0,0 | -1 |  -1 |   1.0 | -1,-1,10 |  True |
    +----------+-----+------+-------+-----------+-------+--------+----+-----+-------+----------+-------+
    |  LOGR10V | 0.0 | Mean |  True | calibrate |  1.62 |    0,0 | -1 |  -1 |   1.0 | -1,-1,10 |  True |
    +----------+-----+------+-------+-----------+-------+--------+----+-----+-------+----------+-------+
    |   RE1025 | 0.0 | Mean | False | calibrate |  1.62 |    0,0 | -1 |  -1 |    -1 | 0,2.5,10 | False |
    +----------+-----+------+-------+-----------+-------+--------+----+-----+-------+----------+-------+
    |  RE1025V | 0.0 | Mean |  True | calibrate |  1.62 |    0,0 | -1 |  -1 |    -1 | 0,2.5,10 | False |
    +----------+-----+------+-------+-----------+-------+--------+----+-----+-------+----------+-------+

    Global binning methods:
    +----------+-----+------+-------+-----------+-------+
    |          | Min |      | Vel.  |     Covar | Covar |
    |      Key | S/N | Stat | Reg.  |      Mode |   par |
    +==========+=====+======+=======+===========+=======+
    |   GLOBAL | 0.0 | Mean | False | calibrate |  1.62 |
    +----------+-----+------+-------+-----------+-------+
    |  GLOBALV | 0.0 | Mean |  True | calibrate |  1.62 |
    +----------+-----+------+-------+-----------+-------+

    Unbinned binning methods:
    +----------+-----+
    |          | Min |
    |      Key | S/N |
    +==========+=====+
    |   SPAXEL | 0.0 |
    +----------+-----+

    Args:
        dapsrc (str): (**Optional**) Root path to the DAP source
            directory.  If not provided, the default is defined by
            :func:`mangadap.config.defaults.default_dap_source`.

    Returns:
        list: A list of
        :func:`mangadap.proc.spatiallybinnedspectra.SpatiallyBinnedSpectraDef`
        objects, each defining a separate binning method.

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
    ini_files = glob.glob(dapsrc+'/python/mangadap/config/spatial_binning/*.ini')
    if len(ini_files) == 0:
        raise IOError('Could not find any configuration files in {0} !'.format(
                      dapsrc+'/python/mangadap/config/spatial_binning'))

    # Build the list of library definitions
    binning_methods = []
    for f in ini_files:
        # Read the config file
        cnfg = ConfigParser(allow_no_value=True)
        cnfg.read(f)
        # Ensure it has the necessary elements to define the template
        # library
        validate_spatial_binning_scheme_config(cnfg)
        if cnfg['default']['method'] == 'global':
            binpar = None
            binclass = spatialbinning.GlobalBinning()
            binfunc = binclass.bin_index
        elif cnfg['default']['method'] == 'radial':
            center = [ float(e.strip()) for e in cnfg['default']['center'].split(',') ]
            radii = [ float(e.strip()) for e in cnfg['default']['radii'].split(',') ]
            radii[2] = int(radii[2])
            binpar = spatialbinning.RadialBinningPar(center, cnfg['default'].getfloat('pa'),
                                                     cnfg['default'].getfloat('ell'),
                                                     cnfg['default'].getfloat('radius_scale'),
                                                     radii, cnfg['default'].getboolean('log_step'))
            binclass = spatialbinning.RadialBinning()
            binfunc = binclass.bin_index
        elif cnfg['default']['method'] == 'voronoi':
            covar = 1.0 if cnfg['default']['noise_calib'] is None \
                        else cnfg['default'].getfloat('noise_calib')
            binpar = spatialbinning.VoronoiBinningPar(cnfg['default'].getfloat('target_snr'),
                                                      None, None, covar)
            binclass = spatialbinning.VoronoiBinning()
            binfunc = binclass.bin_index
        else:   # Do not bin!
            binpar = None
            binclass = None
            binfunc = None
            
        stackpar = SpectralStackPar(cnfg['default']['operation'],
                                    cnfg['default'].getboolean('velocity_register'),
                                    None,
                                    cnfg['default']['stack_covariance_mode'],
                                    SpectralStack.parse_covariance_parameters(
                                            cnfg['default']['stack_covariance_mode'],
                                            cnfg['default']['stack_covariance_par']))
        stackclass = SpectralStack()
        stackfunc = stackclass.stack_DRPFits

        binning_methods += [ SpatiallyBinnedSpectraDef(cnfg['default']['key'],
                                                       cnfg['default'].getfloat('minimum_snr'),
                                                       binpar, binclass, binfunc, stackpar,
                                                       stackclass, stackfunc) ]

    # Check the keywords of the libraries are all unique
    if len(numpy.unique( numpy.array([ method['key'] for method in binning_methods ]) )) \
            != len(binning_methods):
        raise KeyError('Spatial binning method keywords are not all unique!')

    # Return the default list of assessment methods
    return binning_methods


class SpatiallyBinnedSpectraBitMask(BitMask):
    r"""
    Derived class that specifies the mask bits for the spatially binned
    spectra.  See :class:`mangadap.util.bitmask.BitMask` for attributes.

    A list of the bits and meanings are provided by the base class
    function :func:`mangadap.util.bitmask.BitMask.info`; i.e.,::

        from mangadap.proc.spatiallybinnedspectra import SpatiallyBinnedSpectraBitMask
        bm = SpatiallyBinnedSpectraBitMask()
        bm.info()

    """
    def __init__(self, dapsrc=None):
        dapsrc = default_dap_source() if dapsrc is None else str(dapsrc)
        BitMask.__init__(self, ini_file=os.path.join(dapsrc, 'python', 'mangadap', 'config',
                                                     'bitmasks',
                                                     'spatially_binned_spectra_bits.ini'))


class SpatiallyBinnedSpectra:
    r"""

    Class that holds spatially binned spectra.

    .. todo::

        - Add prior with velocity offsets for spaxels.
   
    """
    def __init__(self, method_key, drpf, rdxqa, reff=None, method_list=None, dapsrc=None,
                 dapver=None, analysis_path=None, directory_path=None, output_file=None,
                 hardcopy=True, clobber=False, verbose=0, checksum=False):

        self.version = '1.0'
        self.verbose = verbose

        # Define the method properties
        self.method = None
        self._define_method(method_key, method_list=method_list, dapsrc=dapsrc)

        self.drpf = None
        self.rdxqa = None
        self.reff = None

        # Define the output directory and file
        self.directory_path = None      # Set in _set_paths
        self.output_file = None
        self.hardcopy = None

        # Define the bitmask
        self.bitmask = SpatiallyBinnedSpectraBitMask(dapsrc=dapsrc)

        # Initialize the main class attributes
        self.hdu = None
        self.checksum = checksum
        self.shape = None
        self.spatial_shape = None
        self.nspec = None
        self.spatial_index = None
        self._assign_spectral_arrays()
        self.dispaxis = None
        self.nwave = None

        self.nbins = None
        self.missing_bins = None

        self.covariance = None

        # Bin the spectra
        self.bin_spectra(drpf, rdxqa, reff=reff, dapver=dapver, analysis_path=analysis_path,
                         directory_path=directory_path, output_file=output_file, hardcopy=hardcopy,
                         clobber=clobber, verbose=verbose)


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
        self.method = _select_proc_method(method_key, SpatiallyBinnedSpectraDef,
                                          method_list=method_list,
                                          available_func=available_spatial_binning_methods,
                                          dapsrc=dapsrc)


    def _fill_method_par(self, good_spec):
        """

        Fill in any remaining binning parameters.

        """
        if self.method['binclass'] is None:
            return

        # For the radial binning, fill in the isophotal and scaling
        # parameters, if necessary
        if self.method['binclass'].bintype == 'radial':
            if self.method['binpar']['pa'] < 0:
                self.method['binpar']['pa'] = self.rdxqa.pa
            if self.method['binpar']['ell'] < 0:
                self.method['binpar']['ell'] = self.rdxqa.ell
            if self.method['binpar']['radius_scale'] < 0:
                self.method['binpar']['radius_scale'] = 1.0 if self.reff is None else self.reff

        # For the Voronoi binning type, add the signal and noise (or
        # covariance)
        if self.method['binclass'].bintype == 'voronoi':
            self.method['binpar']['signal'] = self.rdxqa['SPECTRUM'].data['SIGNAL'][good_spec]
            if self.rdxqa.correlation is not None:
                # Overwrite any existing calibration coefficient
                self.rdxqa.correlation.revert_correlation()
                covar = self.rdxqa.correlation.toarray()[good_spec,:][:,good_spec]
                self.rdxqa.correlation.to_correlation()
                i, j = numpy.meshgrid(numpy.arange(covar.shape[0]), numpy.arange(covar.shape[1]))
                self.method['binpar']['covar'] = Covariance(inp=sparse.coo_matrix(
                                                (covar[covar > 0], (i[covar > 0], j[covar > 0])),
                                                shape=covar.shape).tocsr())
            else:
                self.method['binpar']['noise'] \
                        = numpy.sqrt(self.rdxqa['SPECTRUM'].data['VARIANCE'][good_spec])

        # Nothing to add for binning types 'none' or 'global', or
        # user-defined function!



    def _set_paths(self, directory_path, dapver, analysis_path, output_file):
        """
        Set the I/O path to the processed template library.  Used to set
        :attr:`directory_path` and :attr:`output_file`.  If not
        provided, the defaults are set using, respectively,
        :func:`mangadap.config.defaults.default_dap_reference_path` and
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
        self.directory_path = default_dap_reference_path(plate=self.drpf.plate,
                                                         drpver=self.drpf.drpver,
                                                         dapver=dapver,
                                                         analysis_path=analysis_path) \
                                        if directory_path is None else str(directory_path)

        # Set the output file
        method = '{0}-{1}'.format(self.rdxqa.method['key'], self.method['key'])
        self.output_file = default_dap_file_name(self.drpf.plate, self.drpf.ifudesign,
                                                 self.drpf.mode, method) \
                                        if output_file is None else str(output_file)
                                                 #, compressed=False) \


    def _clean_drp_header(self, ext='FLUX'):
        """
        Read and clean the header from the DRP fits file extension for
        use in the output file for this class object.

        .. todo::
            - Currently just returns the existing header.  Need to
              decide if/how to manipulate DRP header.

        """
        return self.drpf[ext].header


    def _initialize_header(self, hdr):
        """

        Initialize the header.

        """
        if self.reff is not None:
            hdr['REFF'] = self.reff
        hdr['BINKEY'] = (self.method['key'], 'Spectal binning method keyword')
        hdr['BINMINSN'] = (self.method['minimum_snr'], 'Minimum S/N of spectrum to include')
        hdr['FSPECCOV'] = (0.8, 'Minimum allowed fraction of good pixels')
        hdr['NBINS'] = (self.nbins, 'Number of unique spatial bins')
        if len(self.missing_bins) > 0:
            hdr['EMPTYBIN'] = (str(self.missing_bins), 'List of bins with no data')
        return hdr


    def _check_fgoodpix(self, minimum_fraction=0.8):
        return self.rdxqa['SPECTRUM'].data['FGOODPIX'] > minimum_fraction


    def _check_snr(self):
        return self.rdxqa['SPECTRUM'].data['SNR'] > self.method['minimum_snr']


    def _initialize_mask(self, good_fgoodpix, good_snr):
        """

        Initialize the mask be setting the DIDNOTUSE, FORESTAR, LOW_SPECCOV, and LOW_SNR masks

        """
        # Initialize to all zeros
        mask = numpy.zeros(self.shape, dtype=self.bitmask.minimum_uint_dtype())

        # Turn on the flag stating that the pixel wasn't used
        indx = self.drpf.bitmask.flagged(self.drpf['MASK'].data,
                                         flag=self.drpf.do_not_stack_flags())
        mask[indx] = self.bitmask.turn_on(mask[indx], 'DIDNOTUSE')

        # Turn on the flag stating that the pixel has a foreground star
        indx = self.drpf.bitmask.flagged(self.drpf['MASK'].data, flag='FORESTAR')
        mask[indx] = self.bitmask.turn_on(mask[indx], 'FORESTAR')

        # Turn on the flag stating that the number of valid channels in
        # the spectrum was below the input fraction.
        indx = numpy.array([numpy.invert(good_fgoodpix).reshape(self.spatial_shape).T]*self.nwave).T
        mask[indx] = self.bitmask.turn_on(mask[indx], flag='LOW_SPECCOV')

        # Turn on the flag stating that the S/N in the spectrum was
        # below the requested limit
        indx = numpy.array([numpy.invert(good_snr).reshape(self.spatial_shape).T]*self.nwave).T
        mask[indx] = self.bitmask.turn_on(mask[indx], flag='LOW_SNR')

        return mask


    def _per_bin_dtype(self):
        r"""
        Construct the record array data type for the output fits
        extension.
        """
        return [ ('BIN_INDEX',numpy.int),
                 ('NBIN',numpy.int),
                 ('SKY_COO',numpy.float,(2,)),
                 ('LW_SKY_COO',numpy.float,(2,)),
                 ('ELL_COO',numpy.float,(2,)),
                 ('LW_ELL_COO',numpy.float,(2,)),
                 ('AREA',numpy.float),
                 ('AREA_FRAC',numpy.float),
                 ('SIGNAL',numpy.float),
                 ('VARIANCE',numpy.float),
                 ('SNR',numpy.float)
               ]


    def _unbinned_data(self, bin_indx):
        """
        Construct the output data table just using the unbinned data.
        """
        good_spec = bin_indx > -1
        self.nbins = numpy.amax(bin_indx)+1
        bin_data = init_record_array(self.nbins, self._per_bin_dtype())

        bin_data['BIN_INDEX'] = numpy.arange(self.nbins)
        bin_data['NBIN'] = numpy.ones(self.nbins, dtype=numpy.int)
        bin_data['SKY_COO'][:,0] = self.rdxqa['SPECTRUM'].data['SKY_COO'][good_spec,0]
        bin_data['SKY_COO'][:,1] = self.rdxqa['SPECTRUM'].data['SKY_COO'][good_spec,1]
        bin_data['LW_SKY_COO'] = bin_data['SKY_COO'].copy()
        bin_data['ELL_COO'][:,0] = self.rdxqa['SPECTRUM'].data['ELL_COO'][good_spec,0]
        bin_data['ELL_COO'][:,1] = self.rdxqa['SPECTRUM'].data['ELL_COO'][good_spec,1]
        bin_data['LW_ELL_COO'] = bin_data['ELL_COO'].copy()

        if self.drpf.pixelscale is None:
            self.drpf.pixelscale = default_cube_pixelscale()

        bin_data['AREA'] = numpy.full(self.nbins, self.drpf.pixelscale, dtype=numpy.float)
        bin_data['AREA_FRAC'] = numpy.ones(self.nbins, dtype=numpy.int)

        bin_data['SIGNAL'] = self.rdxqa['SPECTRUM'].data['SIGNAL'][good_spec]
        bin_data['VARIANCE'] = self.rdxqa['SPECTRUM'].data['VARIANCE'][good_spec]
        bin_data['SNR'] = self.rdxqa['SPECTRUM'].data['SNR'][good_spec]

        return bin_data


    def _fill_covariance(self, stack_covar, bin_indx):
        nspec = len(bin_indx)
        unique_bins, reconstruct = numpy.unique(bin_indx, return_inverse=True)
        indx = bin_indx > -1
        nchan = stack_covar.shape[0]
        self.covariance = numpy.empty(nchan, dtype=sparse.csr.csr_matrix)
        for i in range(nchan):
            j = stack_covar.input_indx[i]
#            stack_covar.show(plane=j)

            # Input spectrum and covariance indices
            ii = unique_bins[reconstruct[indx]]
            ii_i, ii_j = map( lambda x: x.ravel(), numpy.meshgrid(ii, ii) )

            # Output spectrum and covariance indices
            oi = numpy.arange(nspec)[indx]
            oi_i, oi_j = map( lambda x: x.ravel(), numpy.meshgrid(oi, oi) )

            _covar = numpy.zeros((nspec, nspec), dtype=numpy.float)
            _covar[oi_i, oi_j] = stack_covar.toarray(plane=j)[ii_i,ii_j]
#            pyplot.imshow(_covar, origin='lower', interpolation='nearest')
#            pyplot.colorbar()
#            pyplot.show()
            self.covariance[i] = sparse.triu(_covar).tocsr()
        self.covariance = Covariance(inp=self.covariance, input_indx=stack_covar.input_indx)


    def _extract_stack_covariance(self, bin_indx):
        nspec = len(bin_indx)
        unique_bins, unique_indx = numpy.unique(bin_indx, return_index=True)
        self.nbins = numpy.amax(unique_bins)+1
        indx = bin_indx > -1
        nchan = self.covariance.shape[0]
        _covar = numpy.empty(nchan, dtype=sparse.csr.csr_matrix)
        for i in range(nchan):
            j = self.covariance.input_indx[i]
#            self.covariance.show(plane=j)

            # Input spectrum and covariance indices
            ii = unique_indx[1:]
            ii_i, ii_j = map( lambda x: x.ravel(), numpy.meshgrid(ii, ii) )

            # Output spectrum and covariance indices
            oi = unique_bins[1:]
            oi_i, oi_j = map( lambda x: x.ravel(), numpy.meshgrid(oi, oi) )

            _covar = numpy.zeros((self.nbins, self.nbins), dtype=numpy.float)
            _covar[oi_i, oi_j] = self.covariance.toarray(plane=j)[ii_i,ii_j]
#            pyplot.imshow(_covar, origin='lower', interpolation='nearest')
#            pyplot.colorbar()
#            pyplot.show()
            _covar[i] = sparse.triu(_covar).tocsr()
        return Covariance(inp=covar, input_indx=stack_covar.input_indx)


    def _assign_spectral_arrays(self):
        self.spectral_arrays = [ 'FLUX', 'IVAR', 'MASK', 'FLUXD', 'NPIX' ]

    
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


    def do_not_fit_flags(self):
        return ['DIDNOTUSE', 'FORESTAR', 'LOW_SPECCOV', 'LOW_SNR', 'NONE_IN_STACK']


    # TODO: This function needs to be broken up
    def bin_spectra(self, drpf, rdxqa, reff=None, dapver=None, analysis_path=None,
                    directory_path=None, output_file=None, hardcopy=True, clobber=False, verbose=0):
        """

        Bin and stack the spectra.

        .. todo::
            - Allow to weight by S/N?
        """
        # DRPFits object always needed
        if drpf is None:
            raise ValueError('Must provide DRP file object to compute assessments.')
        if not isinstance(drpf, DRPFits):
            raise TypeError('Must provide a valid DRPFits object!')
        if drpf.hdu is None:
            warnings.warn('DRP file previously unopened.  Reading now.')
            drpf.open_hdu()
        self.drpf = drpf
        self.shape = self.drpf.shape
        self.spatial_shape =self.drpf.spatial_shape
        self.nspec = self.drpf.nspec
        self.spatial_index = self.drpf.spatial_index.copy()
        self.dispaxis = self.drpf.dispaxis
        self.nwave = self.drpf.nwave

        # ReductionAssessment object always needed
        if rdxqa is None:
            raise ValueError('Must provide a ReductionAssessment object to bin spectra.')
        if not isinstance(rdxqa, ReductionAssessment):
            raise TypeError('Must provide a valid ReductionAssessment object!')
        self.rdxqa = rdxqa

        # Save the effective radius if provided.  Only used if/when
        # scaling the radii by the effective radius in the radial
        # binning approach
        if reff is not None:
            self.reff = reff

        # Get the good spectra
        #   - Must have valid pixels over more than 80% of the spectral
        #   range 
        good_fgoodpix = self._check_fgoodpix()
        #   - Must have sufficienct S/N, as defined by the input par
        good_snr = self._check_snr()
        # Good spectra have both these criteria
        good_spec = good_fgoodpix & good_snr
        # Report
        print('Total spectra: ', len(good_fgoodpix))
        print('With 80% coverage: ', numpy.sum(good_fgoodpix))
        print('With good S/N: ', numpy.sum(good_snr))
        print('Good spectra: ', numpy.sum(good_spec))

        # Fill in any remaining binning parameters
        self._fill_method_par(good_spec)

        # (Re)Set the output paths
        self._set_paths(directory_path, dapver, analysis_path, output_file)

        # No binning so just fill the hdu with the appropriate data from
        # the DRP file
        if self.method['binpar'] is None and self.method['binclass'] is None \
                and self.method['binfunc'] is None:
        
            # Initialize the header keywords
            hdr = self._clean_drp_header(ext='PRIMARY')
            self._initialize_header(hdr)

            # Initialize the basics of the mask
            mask = self._initialize_mask(good_fgoodpix, good_snr)
            npix = numpy.ones(self.drpf['FLUX'].data.shape, dtype=numpy.int)
            npix[mask > 0] = 0

            # Generate pseudo bin index
            bin_indx = numpy.full(self.drpf.spatial_shape, -1, dtype=numpy.int)
            i = numpy.asarray(tuple(drpf.spatial_index[good_spec]))
            bin_indx[i[:,0],i[:,1]] = numpy.arange(numpy.sum(good_spec))
#            pyplot.imshow(bin_indx.T, origin='lower', interpolation='nearest')
#            pyplot.show()

            # Grab the per-spectrum data
            bin_data = self._unbinned_data(bin_indx.ravel())

            # Set the HDUList just based on the original DRP data
            print('Building HDUList...', end='\r')
            self.hdu = fits.HDUList([ fits.PrimaryHDU(header=hdr), self.drpf['FLUX'],
                                      self.drpf['IVAR'], 
                                      fits.ImageHDU(data=mask, header=self.drpf['MASK'].header,
                                                    name='MASK'),
                                      self.drpf['WAVE'],
                                      fits.ImageHDU(data=numpy.zeros(self.drpf['FLUX'].data.shape),
                                                    header=self.drpf['FLUX'].header, name='FLUXD'),
                                      fits.ImageHDU(data=npix, header=self.drpf['FLUX'].header,
                                                    name='NPIX'),
                                      fits.ImageHDU(data=bin_indx, name='BINID'),
                                      fits.BinTableHDU.from_columns( [ fits.Column(name=n,
                                                             format=rec_to_fits_type(bin_data[n]),
                                                array=bin_data[n]) for n in bin_data.dtype.names ],
                                                               name='BINS')
                                      ])
            print('Building HDUList... done')

            # Try to add the covariance data if requested
            if self.method['stackpar']['covar_mode'] not in ['none', None]:
                try:
                    covar = self.method['stackclass'].build_covariance_data_DRPFits(self.drpf,
                                                        self.method['stackpar']['covar_mode'], 
                                                        self.method['stackpar']['covar_par'])
                except AttributeError as e:
                    warnings.warn('Could not build covariance data:: AttributeError: {0}'.format(e))
                    covar = None

                if isinstance(covar, Covariance):
                    cov_hdr = fits.Header()
                    indx_col, covar_col, var_col, plane_col = covar.binary_columns(hdr=cov_hdr)
                    self.hdu += [ fits.BinTableHDU.from_columns( ([ indx_col, covar_col ] \
                                                                  if var_col is None else \
                                                                [ indx_col, var_col, covar_col ]),
                                                                name='COVAR', header=cov_hdr) ]
                    if plane_col is not None:
                        self.hdu += [ fits.BinTableHDU.from_columns([plane_col],
                                                                    name='COVAR_PLANE') ]

            # Never save a hard copy of the per-spectrum data; already
            # saved in the ReductionAssessment and DRPFits objects
            self.hardcopy = False
            return

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
            self.read(checksum=self.checksum)
            if reff is not None and self.reff != reff:
                warnings.warn('Provided effective radius different from available file; set ' \
                              'clobber=True to overwrite.')
            return

        # Bin the spectra.  To be included in any bin, the spectra must
        # be selected as 'good' according to the selections made above.
        print('Binning spectra ...')
        x = self.rdxqa['SPECTRUM'].data['SKY_COO'][:,0]
        y = self.rdxqa['SPECTRUM'].data['SKY_COO'][:,1]
        bin_indx = numpy.full(self.drpf.nspec, -1, dtype=numpy.int)
        bin_indx[good_spec] = self.method['binfunc'](x[good_spec], y[good_spec],
                                                     par=self.method['binpar'])

        warnings.warn('You\'re forcing bins 2 and 3 to be empty!')
        bin_indx[ (bin_indx == 2) | (bin_indx == 3) ] = 1

        # Stack the spectra in each bin
        print('Stacking spectra ...')
        stack_wave, stack_flux, stack_sdev, stack_npix, stack_ivar, stack_covar = \
                self.method['stackfunc'](self.drpf, bin_indx, par=self.method['stackpar'])

#        print(stack_covar.input_indx)
#        stack_covar.show(plane=0)

#        indx = int(numpy.where(bin_indx == 4)[0][0])
#        pyplot.step(self.drpf['WAVE'].data, self.drpf.select(indx), where='mid',
#                    linestyle='-', color='b')
#        pyplot.step(stack_wave, stack_flux[4,:], where='mid', linestyle='-', color='g')
#        pyplot.show()
#        pyplot.step(self.drpf['WAVE'].data, self.drpf.select(indx, ext='IVAR'), where='mid',
#                    linestyle='-', color='b')
#        pyplot.step(stack_wave, stack_ivar[4,:], where='mid', linestyle='-', color='g')
#        pyplot.show()


#        pyplot.scatter(stack_signal, stack_noise, marker='.', color='k', s=40, lw=0)
#        pyplot.show()

        # Get the unique bins, how to reconstruct the bins from the
        # unique set, and the number of spectra in each bin
        unique_bins, reconstruct, bin_count = numpy.unique(bin_indx, return_inverse=True,
                                                           return_counts=True)

        # Get the number of returned bins, based on their ID number.
        # - The number of returned bins MAY NOT BE THE SAME as the
        #   number of requested bins.  E.g., if radial bins were
        #   requested out to 2.5 Reff, when coverage only goes to 1.5
        #   Reff.
        # - There may also be bins with NO spectra; it is expected that
        #   the stacking function DOES NOT return empty spectra but only
        #   produces a spectrum for the unique bin IDs in the list
        #   provided.  Missing bins interspersed with other bins are
        #   filled in; however, expected bins with indices higher then
        #   self.nbin-1 will not be.
        self.nbins = numpy.amax(unique_bins)+1
        self.missing_bins = list(set(numpy.arange(self.nbins)) - set(unique_bins))

        # Recalculate the mean signal, mean noise, and mean S/N of the
        # binned spectra
        #   - stack_ivar is always expected to exist!
        #   - all returned arrays are expected to be MaskedArrays!
        # !! THIS IS OVER ALL WAVELENGTHS !!
#        print(numpy.sum(stack_ivar.mask))
#        print(numpy.sum(stack_flux.mask))
#        assert numpy.all(stack_ivar.mask == stack_flux.mask)

        _mask = numpy.ma.getmaskarray(stack_flux) | numpy.ma.getmaskarray(stack_ivar) \
                    | numpy.array([
                        self.drpf.wavelength_mask(waverange=self.rdxqa.method['waverange'])
                                  ]*stack_flux.shape[0])
        _stack_flux = numpy.ma.MaskedArray(stack_flux.data, mask=_mask)
        _stack_ivar = numpy.ma.MaskedArray(stack_ivar.data, mask=_mask)
        stack_nsum = numpy.sum(numpy.invert(_mask), axis=1)
        stack_signal = numpy.ma.sum(_stack_flux,axis=1) / stack_nsum
        stack_variance = numpy.ma.sum(numpy.ma.power(_stack_ivar, -1.), axis=1)/stack_nsum
        stack_snr = numpy.ma.sum(_stack_flux*numpy.ma.sqrt(_stack_ivar), axis=1)/stack_nsum
        del _mask, _stack_flux, _stack_ivar

#        pyplot.scatter(stack_snr, stack_signal, marker='.', color='k', s=40, lw=0)
#        pyplot.show()

        # Intialize the data for the binned spectra
        bin_data = init_record_array(self.nbins, self._per_bin_dtype())
        bin_data['BIN_INDEX'] = numpy.arange(self.nbins)

        bin_data['SIGNAL'] = _fill_vector(stack_signal, self.nbins, self.missing_bins, -999) \
                                if len(self.missing_bins) > 0 else stack_signal
        bin_data['VARIANCE'] = _fill_vector(stack_variance, self.nbins, self.missing_bins, -999) \
                                if len(self.missing_bins) > 0 else stack_variance
        bin_data['SNR'] = _fill_vector(stack_snr, self.nbins, self.missing_bins, -999) \
                            if len(self.missing_bins) > 0 else stack_snr

        # Sort the list of bin ids and determine where the spectra jump
        # between bins
        srt = numpy.argsort(bin_indx)
        bin_change = numpy.where(numpy.diff(bin_indx[srt]) > 0)[0] + 1

        # Get the signal-weighted
        wsum = numpy.add.reduceat(self.rdxqa['SPECTRUM'].data['SIGNAL'][srt], bin_change)
        wx = numpy.add.reduceat((x*self.rdxqa['SPECTRUM'].data['SIGNAL'])[srt], bin_change)/wsum
        wy = numpy.add.reduceat((y*self.rdxqa['SPECTRUM'].data['SIGNAL'])[srt], bin_change)/wsum
        wr = numpy.add.reduceat((self.rdxqa['SPECTRUM'].data['ELL_COO'][:,0] \
                                    * self.rdxqa['SPECTRUM'].data['SIGNAL'])[srt], bin_change)/wsum
        wt = numpy.add.reduceat((self.rdxqa['SPECTRUM'].data['ELL_COO'][:,1] \
                                    * self.rdxqa['SPECTRUM'].data['SIGNAL'])[srt], bin_change)/wsum

        # ... and unweighted on-sky coordinates
        bn = bin_count[1:].astype(int)
        bx = numpy.add.reduceat(x[srt], bin_change)/bn
        by = numpy.add.reduceat(y[srt], bin_change)/bn
        br = numpy.add.reduceat(self.rdxqa['SPECTRUM'].data['ELL_COO'][srt,0], bin_change)/bn
        bt = numpy.add.reduceat(self.rdxqa['SPECTRUM'].data['ELL_COO'][srt,1], bin_change)/bn

        # The mean azimuth can go wrong when bins cross the major axis
        _wt = numpy.add.reduceat(((self.rdxqa['SPECTRUM'].data['ELL_COO'][:,1]+180) \
                                    * self.rdxqa['SPECTRUM'].data['SIGNAL'])[srt], bin_change)/wsum
#        pyplot.scatter( wt, _wt-180-wt, marker='.', color='r', s=40, lw=0 )
        indx = numpy.absolute(_wt-180-wt) > 1.0         # HARDWIRED
        wt[indx] = _wt[indx]-180

        _bt = numpy.add.reduceat(self.rdxqa['SPECTRUM'].data['ELL_COO'][srt,1]+180, bin_change)/bn
#        pyplot.scatter( bt, _bt-180-bt, marker='.', color='b', s=40, lw=0 )
#        pyplot.show()
        indx = numpy.absolute(_bt-180-bt) > 1.0         # HARDWIRED
        bt[indx] = _bt[indx]-180

        # Compute the area covered by each bin
        area = self.drpf.binned_on_sky_area(bin_indx, x=x, y=y)

        # Some bins may not have spectra; _fill_vector() inserts default values
        bin_data['NBIN'] = _fill_vector(bn, self.nbins, self.missing_bins, -999) \
                            if len(self.missing_bins) > 0 else bn
        bin_data['SKY_COO'][:,0] = _fill_vector(bx, self.nbins, self.missing_bins, -999) \
                                    if len(self.missing_bins) > 0 else bx
        bin_data['SKY_COO'][:,1] = _fill_vector(by, self.nbins, self.missing_bins, -999) \
                                    if len(self.missing_bins) > 0 else by
        bin_data['ELL_COO'][:,0] = _fill_vector(br, self.nbins, self.missing_bins, -999) \
                                    if len(self.missing_bins) > 0 else br
        bin_data['ELL_COO'][:,1] = _fill_vector(bt, self.nbins, self.missing_bins, -999) \
                                    if len(self.missing_bins) > 0 else bt
        bin_data['LW_SKY_COO'][:,0] = _fill_vector(wx, self.nbins, self.missing_bins, -999) \
                                    if len(self.missing_bins) > 0 else wx
        bin_data['LW_SKY_COO'][:,1] = _fill_vector(wy, self.nbins, self.missing_bins, -999) \
                                    if len(self.missing_bins) > 0 else wy
        bin_data['LW_ELL_COO'][:,0] = _fill_vector(wr, self.nbins, self.missing_bins, -999) \
                                    if len(self.missing_bins) > 0 else wr
        bin_data['LW_ELL_COO'][:,1] = _fill_vector(wt, self.nbins, self.missing_bins, -999) \
                                    if len(self.missing_bins) > 0 else wt
        bin_data['AREA'] = _fill_vector(area, self.nbins, self.missing_bins, -999) \
                            if len(self.missing_bins) > 0 else area

        # Calculate the fractional area of the bin covered by the
        # spectra, if possible; if not, the fractional area is unity
        try:
            bin_area = self.method['binclass'].bin_area()
        except AttributeError as e:
            warnings.warn('Could not calculate nominal bin area:: AttributeError: {0}'.format(e))
            bin_area = bin_data['AREA']
        bin_data['AREA_FRAC'] = bin_data['AREA']/bin_area[:self.nbins]

        # Restructure the output to match the input DRP file
        indx = bin_indx > -1

        flux = numpy.zeros(self.shape, dtype=numpy.float)
        flux.reshape(-1,self.nwave)[indx,:] = stack_flux[unique_bins[reconstruct[indx]],:]

        mask = numpy.zeros(self.shape, dtype=numpy.bool)
        if isinstance(stack_flux, numpy.ma.MaskedArray):
            mask.reshape(-1,self.nwave)[indx,:] \
                = numpy.ma.getmaskarray(stack_flux)[unique_bins[reconstruct[indx]],:]

        sdev = numpy.zeros(self.shape, dtype=numpy.float)
        sdev.reshape(-1,self.nwave)[indx,:] = stack_sdev[unique_bins[reconstruct[indx]],:]

        ivar = numpy.zeros(self.shape, dtype=numpy.float)
        ivar.reshape(-1,self.nwave)[indx,:] = stack_ivar[unique_bins[reconstruct[indx]],:]

        npix = numpy.zeros(self.shape, dtype=numpy.int)
        npix.reshape(-1,self.nwave)[indx,:] = stack_npix[unique_bins[reconstruct[indx]],:]

        self.covariance = None
        if stack_covar is not None:
            self._fill_covariance(stack_covar, bin_indx)

#        print(self.covariance.input_indx)
#        self.covariance.show(plane=self.covariance.input_indx[0])

#        pyplot.imshow(flux.reshape(self.drpf['FLUX'].shape)[:,:,1000].T, origin='lower',
#                      interpolation='nearest')
#        pyplot.show()
#        pyplot.step(stack_wave, flux.reshape(self.drpf['FLUX'].shape)[20,20,:], where='mid',
#                    linestyle='-', color='b')
#        pyplot.show()

        # Initialize the header keywords
        hdr = self._clean_drp_header(ext='PRIMARY')
        self._initialize_header(hdr)
        if self.method['binclass'] is not None:
            try:
                hdr['BINTYPE'] = (self.method['binclass'].bintype, 'Binning method')
            except AttributeError:
                if hardcopy:
                    warnings.warn('Binning parameter class has no attribute bintype.  No type ' \
                                  'written to the header of the output fits file.')

        if self.method['binpar'] is not None:
            if hardcopy and not callable(self.method['binpar'].toheader):
                warnings.warn('Binning parameter class does not have toheader() function.  ' \
                              'No binning parameters written to the header of the output ' \
                              'fits file.')
            elif hardcopy:
                self.method['binpar'].toheader(hdr)

        if self.method['stackpar'] is not None:
            if hardcopy and not callable(self.method['stackpar'].toheader):
                warnings.warn('Stacking parameter class does not have toheader() function.  ' \
                              'No stacking parameters written to the header of the output ' \
                              'fits file.')
            elif hardcopy:
                self.method['stackpar'].toheader(hdr)

        # Initialize the basics of the mask
        mask_bit_values = self._initialize_mask(good_fgoodpix, good_snr)

        # Add flags based on the availability of pixels to stack
        indx = numpy.invert(npix>0)
        mask_bit_values[indx] = self.bitmask.turn_on(mask_bit_values[indx],
                                                     ['NONE_IN_STACK', 'NO_STDDEV'])
        indx = numpy.invert(npix>1)
        mask_bit_values[indx] = self.bitmask.turn_on(mask_bit_values[indx], 'NO_STDDEV')

        # DEBUG
        if numpy.sum( (mask_bit_values == 0) & mask ) > 0:
            raise ValueError('Something should have been masked!')

#        pyplot.imshow(bin_indx.reshape(self.drpf.spatial_shape).T, origin='lower',
#                      interpolation='nearest')
#        pyplot.show()

        # Save the data to the hdu attribute
        self.hdu = fits.HDUList([ fits.PrimaryHDU(header=hdr),
                                  fits.ImageHDU(data=flux.data,
                                                header=self.drpf['FLUX'].header, name='FLUX'),
                                  fits.ImageHDU(data=ivar.data,
                                                header=self.drpf['IVAR'].header, name='IVAR'),
                                  fits.ImageHDU(data=mask_bit_values,
                                                header=self.drpf['MASK'].header, name='MASK'),
                                  self.drpf['WAVE'],
                                  fits.ImageHDU(data=sdev.data,
                                                header=self.drpf['FLUX'].header, name='FLUXD'),
                                  fits.ImageHDU(data=npix.data,
                                                header=self.drpf['FLUX'].header, name='NPIX'),
                                  fits.ImageHDU(data=bin_indx.reshape(self.drpf.spatial_shape),
                                                name='BINID'),
                                  fits.BinTableHDU.from_columns( [ fits.Column(name=n,
                                                             format=rec_to_fits_type(bin_data[n]),
                                                array=bin_data[n]) for n in bin_data.dtype.names ],
                                                               name='BINS')
                                ])

        # Fill the covariance matrix
        if self.covariance is not None:
            cov_hdr = fits.Header()
            indx_col, covar_col, var_col, plane_col = self.covariance.binary_columns(hdr=cov_hdr)
            self.hdu += [ fits.BinTableHDU.from_columns( ([ indx_col, covar_col ] \
                                                            if var_col is None else \
                                                            [ indx_col, var_col, covar_col ]),
                                                        name='COVAR', header=cov_hdr) ]
            if plane_col is not None:
                self.hdu += [ fits.BinTableHDU.from_columns([plane_col], name='COVAR_PLANE') ]

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

        print('restructuring')
        # Restructure the data to match the DRPFits file
        if self.drpf.mode == 'CUBE':
            DRPFits._restructure_cube(self.hdu, ext=self.spectral_arrays, inverse=True)
        elif self.drpf.mode == 'RSS':
            DRPFits._restructure_rss(self.hdu, ext=self.spectral_arrays, inverse=True)

        # Get the output file and determine if it should be compressed
        ofile = self.file_path()
        write_hdu(self.hdu, ofile, clobber=clobber, checksum=True)

        # Revert the structure
        print('reverting')
        if self.drpf.mode == 'CUBE':
            DRPFits._restructure_cube(self.hdu, ext=self.spectral_arrays)
        elif self.drpf.mode == 'RSS':
            DRPFits._restructure_rss(self.hdu, ext=self.spectral_arrays)


    def read(self, ifile=None, strict=True, checksum=False):
        """

        Read an existing file with a previously binned set of spectra.

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
        if self.hdu['PRIMARY'].header['BINKEY'] != self.method['key']:
            if strict:
                raise ValueError('Keywords in header does not match specified method keyword!')
            else:
                warnings.warn('Keywords in header does not match specified method keyword!')
        # TODO: "strict" should also check other aspects of the file to
        # make sure that the details of the method are also the same,
        # not just the keyword

        if self.drpf.mode == 'CUBE':
            DRPFits._restructure_cube(self.hdu, ext=self.spectral_arrays)
        elif self.drpf.mode == 'RSS':
            DRPFits._restructure_rss(self.hdu, ext=self.spectral_arrays)

        # Attempt to read and construct the covariance object
        try:
            self.covariance = Covariance(source=self.hdu, primary_ext='COVAR',
                                         plane_ext='COVAR_PLANE')
        except:
            warnings.warn('Unable to find/read covariance data.')
            self.covariance = None

        # Attempt to read the effective radius
        try:
            self.reff = self.hdu['PRIMARY'].header['REFF']
        except:
            warnings.warn('Unable to read effective radius from file header!')
            self.reff = None

        # Attempt to read binning parameters
        if self.method['binpar'] is not None and callable(self.method['binpar'].fromheader):
            self.method['binpar'].fromheader(self.hdu['PRIMARY'].header)

        # Attempt to stacking parameters
        if self.method['stackpar'] is not None and callable(self.method['stackpar'].fromheader):
            self.method['stackpar'].fromheader(self.hdu['PRIMARY'].header)

        self.nbins = self.hdu['PRIMARY'].header['NBINS']
        try:
            self.missing_bins = eval(self.hdu['PRIMARY'].header['EMPTYBIN'])
        except KeyError:
            # Assume if this fails, it's because the keyword doesn't
            # exist
            self.missing_bins = []

        print(self.missing_bins)


    def unique_bins(self, index=False):
        """
        Get the unique bins or the indices of the unique bins in the
        flattened spatial dimension.
        """
        unique_bins, first_occurrence = numpy.unique(self.hdu['BINID'].data.ravel(),
                                                     return_index=True)
        return first_occurrence[1:] if index else unique_bins[1:]


    def copy_to_array(self, waverange=None, ext='FLUX'):
        r"""
        Return a copy of the selected data array.  The array size is
        always :math:`N_{\rm bins} \times N_{\rm wavelength}`; i.e., the
        data is always flattened to two dimensions and the unique
        spectra are selected.  See :func:`unique_bins` and
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
                                                                self.unique_bins(index=True),:]
        if len(self.missing_bins) == 0:
            return arr
        _arr = numpy.zeros( (self.nbins, self.nwave), dtype=self.hdu[ext].data.dtype)
        _arr[self.unique_bins(),:] = arr
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
                                           flag=flag)[self.unique_bins(index=True),:]
        if len(self.missing_bins) == 0:
            return arr
        _arr = numpy.ma.zeros( (self.nbins, self.nwave), dtype=self.hdu[ext].data.dtype)
        _arr[self.unique_bins(),:] = arr
        _arr[self.missing_bins,:] = numpy.ma.masked
        return _arr



