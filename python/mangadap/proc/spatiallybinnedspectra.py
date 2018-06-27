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
    $MANGADAP_DIR/python/mangadap/proc/spatiallybinnedspectra.py

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
        import os.path
        from os import remove, environ
        from scipy import sparse
        from astropy.io import fits
        import time
        import numpy

        from ..par.parset import ParSet
        from ..config.defaults import dap_source_dir, default_dap_common_path
        from ..config.defaults import default_dap_file_name
        from ..util.geometry import SemiMajorAxisCoo
        from ..util.fileio import init_record_array
        from ..util.parser import DefaultConfig
        from ..drpfits import DRPFits

*Class usage examples*:

    .. todo::
        Add examples

*Revision history*:
    | **01 Apr 2016**: Implementation begun by K. Westfall (KBW)
    | **19 May 2016**: (KBW) Include SPECRES and SPECRESD extensions
        from DRP file in output. Added loggers and quiet keyword
        arguments to :class:`SpatiallyBinnedSpectra`, removed verbose 
    | **06 Jul 2016**: (KBW) Make the application of a reddening
        correction an input parameter.
    | **25 Aug 2016**: (KBW) Fixed error in bin area when calling
        :func:`SpatiallyBinnedSpetra._unbinned_data_table`
    | **01 Dec 2016**: (KBW) Include keyword that describes how to
        handle the spectral resolution.
    | **02 Dec 2016**: (KBW) Incorporate
        :class:`mangadap.util.extinction.GalacticExtinction`.  Revert
        main file structure to be bin-based instead of cube-based;
        include convenience functions that construct the cube for each
        extension as requested.
    | **06 Dec 2016**: (KBW) Significantly restructured.
    | **23 Feb 2017**: (KBW) Use DAPFitsUtil read and write functions.
    | **21 Aug 2017**: (KBW) Use the new PRE-pixelized assessments of
        the LSF.
        
.. _astropy.io.fits.hdu.hdulist.HDUList: http://docs.astropy.org/en/v1.0.2/io/fits/api/hdulists.html
.. _glob.glob: https://docs.python.org/3.4/library/glob.html
.. _logging.Logger: https://docs.python.org/3/library/logging.html
.. _numpy.ma.MaskedArray: http://docs.scipy.org/doc/numpy-1.10.1/reference/maskedarray.baseclass.html#numpy.ma.MaskedArray

.. todo::
    - Check binning an RSS file

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
import time
import logging
import numpy

from scipy import sparse, spatial, interpolate

from astropy.wcs import WCS
from astropy.io import fits

from ..drpfits import DRPFits
from ..par.parset import ParSet
from ..util.fitsutil import DAPFitsUtil
from ..util.fileio import init_record_array, rec_to_fits_type, create_symlink
from ..util.parser import DefaultConfig
from ..util.bitmask import BitMask
from ..util.pixelmask import SpectralPixelMask
from ..util.instrument import spectral_coordinate_step
from ..util.covariance import Covariance
from ..util.geometry import SemiMajorAxisCoo
from ..util.extinction import GalacticExtinction
from ..util.log import log_output
from ..config.defaults import dap_source_dir, default_dap_common_path
from ..config.defaults import default_dap_file_name, default_cube_pixelscale
from . import spatialbinning
from .reductionassessments import ReductionAssessment
from .spectralstack import SpectralStackPar, SpectralStack
from .util import select_proc_method, replace_with_data_from_nearest_coo

from matplotlib import pyplot
#from memory_profiler import profile

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
        
        stack_wave, stack_flux, stack_sdev, stack_npix, stack_ivar, \
                stack_sres, stack_covar = stack(drpf, id, par=par)

    where `drpf` is a :class:`mangadap.drpfits.DRPFits` object.  Note
    that the input and output wavelength ranges are expected to be the
    same!  This is why the wavelengths are not returned!

    As long as they are mutable, the values in par can change, meaning
    that some products of the bin algorithm can be passed to the stack
    algorithm.  For example, if you want to weight the inclusion of the
    spectrum in the bin, you'd have to provide both the binning and
    stacking routines.  Actually, if that's the case, you're better off
    providing a class object that will both bin and stack the spectra!

    .. todo::
        - Impose a set of options of the form for the Galactic
          reddening.

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
        spec_res (str): Keyword defining the treatment of the spectral
            resolution.  See
            :func:`SpatiallyBinnedSpectra.spectral_resolution_options`
            for a list of the options.
        prepixel_sres (bool): Use the prepixelized version of the LSF
            measurements.
    """
    def __init__(self, key, galactic_reddening, galactic_rv, minimum_snr, binpar, binclass,
                 binfunc, stackpar, stackclass, stackfunc, spec_res, prepixel_sres):
        in_fl = [ int, float ]
        res_opt = SpatiallyBinnedSpectra.spectral_resolution_options()
#        bincls_opt = [ spatialbinning.SpatialBinning ]
#        stackcls_opt = [ SpectralStack ]
        par_opt = [ ParSet, dict ]

        pars =     [ 'key', 'galactic_reddening', 'galactic_rv', 'minimum_snr', 'binpar',
                     'binclass', 'binfunc', 'stackpar', 'stackclass', 'stackfunc', 'spec_res',
                     'prepixel_sres' ]
        values =   [ key, galactic_reddening, galactic_rv, minimum_snr, binpar, binclass, binfunc,
                     stackpar, stackclass, stackfunc, spec_res, prepixel_sres ]
        options =  [ None, None, None, None, None, None, None, None, None, None, res_opt, None ]
        dtypes =   [ str, str, in_fl, in_fl, par_opt, None, None, par_opt, None, None, str, bool ]
        can_call = [ False, False, False, False, False, False, True, False, False, True, False,
                     False ]

        ParSet.__init__(self, pars, values=values, options=options, dtypes=dtypes,
                        can_call=can_call)


def validate_spatial_binning_scheme_config(cnfg):
    """ 
    Validate the :class:`mangadap.util.parser.DefaultConfig` object with
    spatial-binning scheme parameters.

    Args:
        cnfg (:class:`mangadap.util.parser.DefaultConfig`): Object with
            the spatial-binning method parameters as needed by
            :class:`mangadap.proc.spatiallybinnedspectra.SpatiallyBinnedSpectraDef`.

    Raises:
        KeyError: Raised if any required keywords do not exist.
        ValueError: Raised if keys have unacceptable values.
        FileNotFoundError: Raised if a file is specified but could not
            be found.
    """
    # Check for required keywords

    required_keywords = [ 'key', 'method', 'spec_res' ]
    if not cnfg.all_required(required_keywords):
        raise KeyError('Keywords {0} must all have valid values.'.format(required_keywords))

    if cnfg['method'] not in [ 'none', 'global', 'voronoi', 'radial' ]:
        raise ValueError('Unknown binning method: {0}'.format(cnfg['method']))

    covar_par_needed_modes = SpectralStack.covariance_mode_options(par_needed=True)
    if cnfg['stack_covariance_mode'] in covar_par_needed_modes \
                and not cnfg.keyword_specified('stack_covariance_par'):
        raise ValueError('must provide a parameter for covariance mode {0}.'.format(
                         cnfg['stack_covariance_mode']))

    if cnfg['method'] in [ 'none', 'global' ]:
        return

    if cnfg['method'] == 'voronoi' and not cnfg.keyword_specified('target_snr'):
        raise KeyError('Keyword \'target_snr\' must be provided for Voronoi binning.')

    required_keywords = [ 'center', 'pa', 'ell', 'radii' ]
    if cnfg['method'] == 'radial' and not cnfg.all_required(required_keywords):
        raise KeyError('Keywords {0} must all have valid values for radial binning.'.format(
                        required_keywords))



def available_spatial_binning_methods(dapsrc=None):
    """
    Return the list of available binning schemes.
    
    The list of available binning schemes is defined by the set of
    configuration files at::

        config_path = os.path.join(dapsrc, 'python', 'mangadap', 'config', 'spatial_binning')
        ini_files = glob.glob(os.path.join(config_path, '/*.ini'))

    Args:
        dapsrc (str): (**Optional**) Root path to the DAP source
            directory.  If not provided, the default is defined by
            :func:`mangadap.config.defaults.dap_source_dir`.

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
    ini_files = glob.glob(dapsrc+'/python/mangadap/config/spatial_binning/*.ini')
    if len(ini_files) == 0:
        raise IOError('Could not find any configuration files in {0} !'.format(
                      dapsrc+'/python/mangadap/config/spatial_binning'))

    # Build the list of library definitions
    binning_methods = []
    for f in ini_files:
        # Read the config file
        cnfg = DefaultConfig(f=f)
        # Ensure it has the necessary elements to define the template
        # library
        validate_spatial_binning_scheme_config(cnfg)
        if cnfg['method'] == 'global':
            binpar = None
            binclass = spatialbinning.GlobalBinning()
            binfunc = binclass.bin_index
        elif cnfg['method'] == 'radial':
            center = cnfg.getlist('center', evaluate=True)
            radii = cnfg.getlist('radii', evaluate=True)
            radii[2] = int(radii[2])
            binpar = spatialbinning.RadialBinningPar(center, cnfg.getfloat('pa', default=0.),
                                                     cnfg.getfloat('ell', default=0.),
                                                     cnfg.getfloat('radius_scale', default=1.),
                                                     radii, cnfg.getbool('log_step', default=False))
            binclass = spatialbinning.RadialBinning()
            binfunc = binclass.bin_index
        elif cnfg['method'] == 'voronoi':
            binpar = spatialbinning.VoronoiBinningPar(cnfg.getfloat('target_snr'), None, None,
                                                      cnfg.getfloat('noise_calib'))
            binclass = spatialbinning.VoronoiBinning()
            binfunc = binclass.bin_index
        else:   # Do not bin!
            binpar = None
            binclass = None
            binfunc = None

        stack_spec_res = cnfg.get('spec_res') == 'spaxel'
        prepixel_sres = cnfg.getbool('prepixel_sres', default=True)

        stackpar = SpectralStackPar(cnfg.get('operation', default='mean'),
                                    cnfg.getbool('velocity_register', default=False), None,
                                    cnfg.get('stack_covariance_mode', default='none'),
                                    SpectralStack.parse_covariance_parameters(
                                            cnfg.get('stack_covariance_mode', default='none'),
                                            cnfg['stack_covariance_par']),
                                    stack_spec_res, prepixel_sres)
        stackclass = SpectralStack()
        stackfunc = stackclass.stack_DRPFits

        binning_methods += [ SpatiallyBinnedSpectraDef(cnfg['key'], cnfg['galactic_reddening'],
                                                       cnfg.getfloat('galactic_rv', default=3.1),
                                                       cnfg.getfloat('minimum_snr', default=0.),
                                                       binpar, binclass, binfunc, stackpar,
                                                       stackclass, stackfunc, cnfg['spec_res'],
                                                       prepixel_sres) ]

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
        dapsrc = dap_source_dir() if dapsrc is None else str(dapsrc)
        BitMask.__init__(self, ini_file=os.path.join(dapsrc, 'python', 'mangadap', 'config',
                                                     'bitmasks',
                                                     'spatially_binned_spectra_bits.ini'))


class SpatiallyBinnedSpectra:
    r"""
    Class that holds spatially binned spectra.

    Args:
        method_key (str): The keyword that designates which method,
            provided in *method_list*, to use for the binning procedure.  
        drpf (:class:`mangadap.drpfits.DRPFits`): The DRP datacube with
            the spectra to bin.
        rdxqa (:class:`mangadap.proc.reductionassessments.ReductionAssessments`):
            The basic assessments of the DRP data that are needed for
            the binning procedures.
        reff (float): (**Optional**) The effective radius of the galaxy.
        method_list (list): (**Optional**) List of
            :class:`SpatiallyBinnedSpectraDef` objects that define one
            or more methods to use for the spatial binning.  Default is
            to use the config files in the DAP source directory to
            construct the available methods using
            :func:`available_spatial_binning_methods`.
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
        symlink_dir (str): (**Optional**) Create a symbolic link to the
            created file in the supplied directory.  Default is to
            produce no symbolic link.
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
        - Allow velocity offsets for registration.
   
    """
#    @profile
    def __init__(self, method_key, drpf, rdxqa, reff=None, method_list=None, dapsrc=None,
                 dapver=None, analysis_path=None, directory_path=None, output_file=None,
                 hardcopy=True, symlink_dir=None, clobber=False, checksum=False, loggers=None,
                 quiet=False):

        self.loggers = None
        self.quiet = False

        # Define the method properties
        self.method = None
        self._define_method(method_key, method_list=method_list, dapsrc=dapsrc)

        self.drpf = None
        self.rdxqa = None
        self.reff = None
        self.galext = None

        # Define the output directory and file
        self.directory_path = None      # Set in _set_paths
        self.output_file = None
        self.hardcopy = None
        self.symlink_dir = None

        # Define the bitmask
        self.bitmask = SpatiallyBinnedSpectraBitMask(dapsrc=dapsrc)

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

        self.nbins = None
        self.missing_bins = None

        self.covariance = None

        # Bin the spectra
        self.bin_spectra(drpf, rdxqa, reff=reff, dapver=dapver, analysis_path=analysis_path,
                         directory_path=directory_path, output_file=output_file, hardcopy=hardcopy,
                         symlink_dir=symlink_dir, clobber=clobber, loggers=loggers, quiet=quiet)


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


    def _define_method(self, method_key, method_list=None, dapsrc=None):
        r"""

        Select the method

        """
        # Grab the specific method
        self.method = select_proc_method(method_key, SpatiallyBinnedSpectraDef,
                                         method_list=method_list,
                                         available_func=available_spatial_binning_methods,
                                         dapsrc=dapsrc)


    def _fill_method_par(self, good_spec):
        """
        Finalize the binning parameters, as needed.

        **For the radial binning**, set the ellipticity and position
        angle.  Set scale radius to the effective radius if desired.

        **For the Voronoi binning**, set the signal and noise.

        Args:
            good_spec (numpy.ndarray):  List of spectra to include in
                the binning.  See :func:`check_fgoodpix` and
                :func:`_check_snr`.
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
#                self.rdxqa.correlation.show()
                covar = self.rdxqa.correlation.toarray()[good_spec,:][:,good_spec]
                self.rdxqa.correlation.to_correlation()
                i, j = numpy.meshgrid(numpy.arange(covar.shape[0]), numpy.arange(covar.shape[1]))
                self.method['binpar']['covar'] = Covariance(
                                inp=sparse.coo_matrix((covar[covar > 0].ravel(),
                                                      (i[covar > 0].ravel(), j[covar > 0].ravel())),
                                                      shape=covar.shape).tocsr())
#                self.method['binpar']['covar'].show()
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
        :func:`mangadap.config.defaults.default_dap_common_path` and
        :func:`mangadap.config.defaults.default_dap_file_name`.

        Args:
            directory_path (str): The exact path to the directory with
                DAP output that is common to number DAP "methods".  See
                :attr:`directory_path`.
            dapver (str): DAP version.
            analysis_path (str): The path to the top-level directory
                containing the DAP output files for a given DRP and DAP
                version.
            output_file (str): The name of the file with the reduction assessments.
                See :func:`compute`.

        """
        # Set the output directory path
        self.directory_path = default_dap_common_path(plate=self.drpf.plate,
                                                      ifudesign=self.drpf.ifudesign,
                                                      drpver=self.drpf.drpver, dapver=dapver,
                                                      analysis_path=analysis_path) \
                                        if directory_path is None else str(directory_path)

        # Set the output file
        method = '{0}-{1}'.format(self.rdxqa.method['key'], self.method['key'])
        self.output_file = default_dap_file_name(self.drpf.plate, self.drpf.ifudesign, method) \
                                        if output_file is None else str(output_file)


    def _initialize_primary_header(self, hdr=None):
        """
        Initialize the header of :attr:`hdu`.

        Returns:
            astropy.io.fits.Header : Edited header object.

        """
        # Copy the from the DRP and clean it
        if hdr is None:
            hdr = self.drpf.hdu['PRIMARY'].header.copy()
            hdr = DAPFitsUtil.clean_dap_primary_header(hdr)
        
        # Add keywords specific to this object
        hdr['AUTHOR'] = 'Kyle B. Westfall <westfall@ucolick.org>'
        if self.reff is not None:
            hdr['REFF'] = self.reff
        hdr['BINKEY'] = (self.method['key'], 'Spectal binning method keyword')
        hdr['BINMINSN'] = (self.method['minimum_snr'], 'Minimum S/N of spectrum to include')
        hdr['FSPCOV'] = (0.8, 'Minimum allowed fraction of good pixels')
        hdr['NBINS'] = (self.nbins, 'Number of unique spatial bins')
#        if len(self.missing_bins) > 0:
#            hdr['EMPTYBIN'] = (str(self.missing_bins), 'List of bins with no data')
        return hdr


    def _add_method_header(self, hdr):
        """
        Add method-specific metadata to the header.
        """
        if self.is_unbinned:
            hdr['BINTYPE'] = ('None', 'Binning method')
        if self.method['binclass'] is not None:
            try:
                hdr['BINTYPE'] = (self.method['binclass'].bintype, 'Binning method')
            except AttributeError:
                if not self.quiet and self.hardcopy:
                    warnings.warn('Binning parameter class has no attribute bintype.  No type ' \
                                  'written to the header of the output fits file.')

        if self.method['binpar'] is not None:
            if not self.quiet and self.hardcopy and not callable(self.method['binpar'].toheader):
                warnings.warn('Binning parameter class does not have toheader() function.  ' \
                              'No binning parameters written to the header of the output ' \
                              'fits file.')
            else:
                hdr = self.method['binpar'].toheader(hdr)

        if self.method['stackpar'] is not None:
            if not self.quiet and self.hardcopy and not callable(self.method['stackpar'].toheader):
                warnings.warn('Stacking parameter class does not have toheader() function.  ' \
                              'No stacking parameters written to the header of the output ' \
                              'fits file.')
            else:
                hdr = self.method['stackpar'].toheader(hdr)
        return hdr


    def _add_reddening_header(self, hdr):
        """
        Set the relevant reddening information to the header.

        Args:
            hdr (astropy.io.fits.Header): Input header object to edit.

        Returns:
            astropy.io.fits.Header : Edited header object.
        """
        hdr['EBVGAL'] = (self.drpf['PRIMARY'].header['EBVGAL'], 'Galactic reddening E(B-V)')
        hdr['GEXTLAW'] = (str(self.galext.form), 'Galactic extinction law')
        hdr['RVGAL'] = (self.galext.rv, 'Ratio of total to selective extinction, R(V)')
        return hdr


    def _initialize_cube_mask(self):
        """
        Initialize the mask, copying the DIDNOTUSE and FORESTAR bits
        from the DRP file, and setting the LOW_SPECCOV and LOW_SNR bits
        based on the input boolean arrays.

        Returns:
            numpy.ndarray : Bitmask array.
        """

        # Get the spaxels with good spectral coverage and S/N
        good_fgoodpix = self.check_fgoodpix()
        good_snr = self._check_snr()

        # Initialize to all zeros
        mask = numpy.zeros(self.shape, dtype=self.bitmask.minimum_dtype())

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


    def _check_snr(self):
        """
        Determine which spectra in :attr:`rdxqa` have a S/N greater than
        the minimum set by :attr:`method`.  Only these spectra will be
        included in the binning.
    
        Returns:
            numpy.ndarray : Boolean array for the spectra that satisfy
            the criterion.
        """
        return self.rdxqa['SPECTRUM'].data['SNR'] > self.method['minimum_snr']


    def _assign_spectral_arrays(self):
        """
        Set :attr:`spectral_arrays`, which contains the list of
        extensions in :attr:`hdu` that contain spectral data.
        """
        self.spectral_arrays = [ 'FLUX', 'IVAR', 'MASK', 'SPECRES', 'FLUXD', 'NPIX' ]


    def _assign_image_arrays(self):
        """
        Set :attr:`image_arrays`, which contains the list of extensions
        in :attr:`hdu` that are on-sky image data.
        """
        self.image_arrays = [ 'BINID' ]


    def _per_bin_dtype(self):
        """
        Construct the record array data type that holds the information
        for the binned spectra.

        .. todo::

            There currently is no mask; however, could add masks for:
                - Insufficient good wavelength channels in the spectrum.
                - Variance in flux in bin is too large.
        """
        return [ ('BINID',numpy.int),
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


    @staticmethod
    def _get_missing_bins(unique_bins):
        return list(set(numpy.arange(numpy.amax(unique_bins)+1)) - set(unique_bins))


    def _unbinned_data_table(self, bin_indx):
        """
        Construct the output data table for the unbinned spectra.

        Args:
            bin_indx (numpy.ndarray): The integer vector with the bin
                associated with each spectrum in the DRP cube.  This is
                the flattened BINID array.

        Returns:
            numpy.recarray : The record array that is put in the BINS
            extension of :attr:`hdu`.
    
        """
        # Find the spectra that have a bin ID
        good_spec = bin_indx > -1
        self.nbins = numpy.amax(bin_indx)+1
        self.missing_bins = []                  # No missing bins

        # Copy the data from the ReductionAssessments object
        bin_data = init_record_array(self.nbins, self._per_bin_dtype())

        bin_data['BINID'] = numpy.arange(self.nbins)
        bin_data['NBIN'] = numpy.ones(self.nbins, dtype=numpy.int)
        bin_data['SKY_COO'][:,0] = self.rdxqa['SPECTRUM'].data['SKY_COO'][good_spec,0]
        bin_data['SKY_COO'][:,1] = self.rdxqa['SPECTRUM'].data['SKY_COO'][good_spec,1]
        bin_data['LW_SKY_COO'] = bin_data['SKY_COO'].copy()
        bin_data['ELL_COO'][:,0] = self.rdxqa['SPECTRUM'].data['ELL_COO'][good_spec,0]
        bin_data['ELL_COO'][:,1] = self.rdxqa['SPECTRUM'].data['ELL_COO'][good_spec,1]
        bin_data['LW_ELL_COO'] = bin_data['ELL_COO'].copy()

        if self.drpf.pixelscale is None:
            self.drpf.pixelscale = default_cube_pixelscale()

        bin_data['AREA'] = numpy.full(self.nbins, numpy.square(self.drpf.pixelscale),
                                      dtype=numpy.float)
        bin_data['AREA_FRAC'] = numpy.ones(self.nbins, dtype=numpy.int)

        bin_data['SIGNAL'] = self.rdxqa['SPECTRUM'].data['SIGNAL'][good_spec]
        bin_data['VARIANCE'] = self.rdxqa['SPECTRUM'].data['VARIANCE'][good_spec]
        bin_data['SNR'] = self.rdxqa['SPECTRUM'].data['SNR'][good_spec]

        return bin_data


    def _interpolated_response_function(self):
        response_func = self.rdxqa.method['response_func']
        if response_func is None:
            # Just return a uniform response function
            return numpy.ones(self.drpf.nwave, dtype=float)
        
        interp = interpolate.interp1d(response_func[:,0], response_func[:,1], bounds_error=False,
                                      fill_value=0.0, assume_sorted=True)
        return interp(self.drpf['WAVE'].data)


    def _binned_data_table(self, bin_indx, stack_flux, stack_ivar, per_pixel=True):
        r"""
        Construct the output data table for the binned spectra.

        Args:
            bin_indx (numpy.ndarray): The integer vector with the bin
                associated with each spectrum in the DRP cube.  This is
                the flattened BINID array.
            stack_flux (numpy.ndarray): The stacked spectra, organized
                as :math:`N_{\rm spec}\times\N_\lambda}`.
            stack_flux (numpy.ndarray): The stacked inverse variance,
                organized as :math:`N_{\rm spec}\times\N_\lambda}`.

        Returns:
            numpy.recarray : The record array that is put in the BINS
            extension of :attr:`hdu`.
    
        """
        # Get the unique bins and the number of spectra in each bin
        unique_bins, bin_count = map(lambda x : x[1:], numpy.unique(bin_indx, return_counts=True))
        unique_bins, unique_indx, bin_count = map(lambda x : x[1:],
                                                  numpy.unique(bin_indx, return_index=True,
                                                               return_counts=True))

        # Get the number of returned bins:
        # - The number of returned bins MAY NOT BE THE SAME as the
        #   number of requested bins.  E.g., if radial bins were
        #   requested out to 2.5 Reff, when coverage only goes to 1.5
        #   Reff.
        # - Bins with no spectra are not included in the data table!
        # - Missing bins are only identified as those without indices in
        #   the range of provided bin numbers.
        self.nbins = len(unique_bins)
        self.missing_bins = SpatiallyBinnedSpectra._get_missing_bins(unique_bins)

        # Intialize the data for the binned spectra
        bin_data = init_record_array(self.nbins, self._per_bin_dtype())
        bin_data['BINID'] = unique_bins.astype(int)
        bin_data['NBIN'] = bin_count.astype(int)

        # Recalculate the mean signal, mean noise, and mean S/N of the
        # binned spectra
#        print(numpy.sum(stack_ivar.mask))
#        print(numpy.sum(stack_flux.mask))
#        assert numpy.all(stack_ivar.mask == stack_flux.mask)
        if stack_ivar is None:
            warnings.warn('No inverse variance for stack.  Errors set to unity.')

        _wavelength_mask = SpectralPixelMask(waverange=
                                self.rdxqa.method['waverange']).boolean(self.drpf['WAVE'].data)
        _mask = numpy.ma.getmaskarray(stack_flux) | _wavelength_mask[None,:]

        if stack_ivar is not None:
            _mask |= numpy.ma.getmaskarray(stack_ivar)

        _stack_flux = numpy.ma.MaskedArray(stack_flux.data, mask=_mask)

        # Same as in mangadap.drpfits.DRPFits.flux_stats() ...
        # Set the response function
        dw = numpy.ones(self.drpf.nwave, dtype=float) if per_pixel else \
                    spectral_coordinate_step(self.drpf['WAVE'].data,
                                             log=True)*numpy.log(10.)*self.drpf['WAVE'].data
        _response_func = self._interpolated_response_function()
        # Get the signal
        response_integral = numpy.sum(numpy.invert(numpy.ma.getmaskarray(_stack_flux))
                                        * (_response_func*dw)[None,:], axis=1)
#        print(bin_data['SIGNAL'].shape)
#        print(_stack_flux.shape)
        bin_data['SIGNAL'] = numpy.ma.divide(numpy.ma.sum(_stack_flux*(_response_func*dw)[None,:],
                                                          axis=1), response_integral).filled(0.0)
        # And the variance and SNR if the inverse variance is available
        if stack_ivar is not None:
            _stack_ivar = numpy.ma.MaskedArray(stack_ivar.data, mask=_mask)
            bin_data['VARIANCE'] = numpy.ma.divide(numpy.ma.sum(numpy.ma.power(_stack_ivar, -1.) \
                                                        * (_response_func*dw)[None,:], axis=1),
                                                            response_integral).filled(0.0)
            bin_data['SNR'] = numpy.ma.divide(numpy.ma.sum(_stack_flux * numpy.ma.sqrt(_stack_ivar)
                                                            * (_response_func*dw)[None,:], axis=1),
                                                            response_integral).filled(0.0)
            del _stack_ivar
        del _mask, _stack_flux

#        indx = bin_indx == 0
#        print(self.rdxqa['SPECTRUM'].data['SIGNAL'][indx])
#        print(bin_data['SIGNAL'][0])
#        print(self.rdxqa['SPECTRUM'].data['VARIANCE'][indx])
#        print(bin_data['VARIANCE'][0])
#        print(self.rdxqa['SPECTRUM'].data['SNR'][indx])
#        print(bin_data['SNR'][0])
#
#        pyplot.scatter(self.rdxqa['SPECTRUM'].data['SNR'][unique_indx], bin_data['SNR'])
#        pyplot.plot([0,100],[0,100],color='k')
#        pyplot.show()
    
#        pyplot.imshow(bin_indx.reshape((int(numpy.sqrt(len(bin_indx))),)*2))
#        pyplot.colorbar()
#        pyplot.show()

#        pyplot.plot(self.drpf['WAVE'].data, _response_func)
#        pyplot.show()

        # TODO: This only works with a limiting wavelength range, not a
        # response function
#        stack_nsum = numpy.sum(numpy.invert(_mask), axis=1)
#        bin_data['SIGNAL'] = numpy.ma.sum(_stack_flux,axis=1) / stack_nsum
#        if stack_ivar is not None:
#            _stack_ivar = numpy.ma.MaskedArray(stack_ivar.data, mask=_mask)
#            bin_data['VARIANCE'] = numpy.ma.sum(numpy.ma.power(_stack_ivar, -1.), axis=1)/stack_nsum
#            bin_data['SNR'] = numpy.ma.sum(_stack_flux * numpy.ma.sqrt(_stack_ivar),
#                                           axis=1)/stack_nsum
#            del _stack_ivar
#        del _mask, _stack_flux


#        pyplot.scatter(bin_data['SNR'], bin_data['SIGNAL'], marker='.', color='k', s=40, lw=0)
#        pyplot.show()

        # Sort the list of bin ids and determine where the spectra jump
        # between bins
        srt = numpy.argsort(bin_indx)
        bin_change = numpy.where(numpy.diff(bin_indx[srt]) > 0)[0] + 1

        # Some convenience reference variables
        x = self.rdxqa['SPECTRUM'].data['SKY_COO'][:,0]
        y = self.rdxqa['SPECTRUM'].data['SKY_COO'][:,1]
        r = self.rdxqa['SPECTRUM'].data['ELL_COO'][:,0]
        t = self.rdxqa['SPECTRUM'].data['ELL_COO'][:,1]
        s = self.rdxqa['SPECTRUM'].data['SIGNAL']

        # Get the signal-weighted on-sky coordinates
        wsum = numpy.add.reduceat(s[srt], bin_change)
        bin_data['LW_SKY_COO'][:,0] = numpy.add.reduceat((x*s)[srt], bin_change)/wsum
        bin_data['LW_SKY_COO'][:,1] = numpy.add.reduceat((y*s)[srt], bin_change)/wsum
        bin_data['LW_ELL_COO'][:,0] = numpy.add.reduceat((r*s)[srt], bin_change)/wsum
        bin_data['LW_ELL_COO'][:,1] = numpy.add.reduceat((t*s)[srt], bin_change)/wsum

        # ... and the unweighted on-sky coordinates
        bin_data['SKY_COO'][:,0] = numpy.add.reduceat(x[srt], bin_change)/bin_data['NBIN']
        bin_data['SKY_COO'][:,1] = numpy.add.reduceat(y[srt], bin_change)/bin_data['NBIN']
        bin_data['ELL_COO'][:,0] = numpy.add.reduceat(r[srt], bin_change)/bin_data['NBIN']
        bin_data['ELL_COO'][:,1] = numpy.add.reduceat(t[srt], bin_change)/bin_data['NBIN']

        # The mean azimuth can go wrong when bins cross the major axis
        _wt = numpy.add.reduceat(((t+180)*s)[srt], bin_change)/wsum
        indx = numpy.absolute(_wt-180-bin_data['LW_ELL_COO'][:,1]) > 1.0        # HARDWIRED
        bin_data['LW_ELL_COO'][indx,1] = _wt[indx]-180

        _bt = numpy.add.reduceat(t[srt]+180, bin_change) / bin_data['NBIN']
        indx = numpy.absolute(_bt-180-bin_data['ELL_COO'][:,1]) > 1.0           # HARDWIRED
        bin_data['ELL_COO'][indx,1] = _bt[indx]-180

#        pyplot.scatter( wt, _wt-180-bin_data['LW_ELL_COO'][:,1], marker='.',color='r',s=40,lw=0 )
#        pyplot.scatter( bt, _bt-180-bin_data['ELL_COO'][:,1], marker='.',color='b',s=40,lw=0 )
#        pyplot.show()

        # Compute the area covered by each bin
        bin_data['AREA'] = self.drpf.binned_on_sky_area(bin_indx, x=x, y=y)

        # Calculate the fractional area of the bin covered by the
        # spectra, if possible; if not, the fractional area is unity
        try:
            bin_total_area = self.method['binclass'].bin_area()[unique_bins]
        except AttributeError as e:
            if not self.quiet:
                warnings.warn('Could not calculate nominal bin area:: '
                              'AttributeError: {0}'.format(e))
            bin_total_area = bin_data['AREA']
        bin_data['AREA_FRAC'] = bin_data['AREA']/bin_total_area

#        pyplot.scatter(bin_data['ELL_COO'][:,0], bin_data['SNR'])
#        pyplot.show()

        return bin_data

   
    def _apply_reddening(self, flux, ivar, sdev, covar, deredden=True):
        if self.galext.form is None:
            return flux, ivar, sdev, covar

        # Apply the reddening correction
        flux, ivar = self.galext.apply(flux, ivar=ivar, deredden=deredden)
        if sdev is not None:
            sdev = self.galext.apply(sdev, deredden=deredden)
        if covar is not None:
            for i,j in enumerate(covar.input_indx):
                covar.cov[i] = covar.cov[i] * numpy.square(self.galext.redcorr[j]) if deredden \
                                    else covar.cov[i] / numpy.square(self.galext.redcorr[j])
        return flux, ivar, sdev, covar


    def _construct_2d_hdu(self, bin_indx, good_fgoodpix, good_snr, bin_data, stack_flux,
                          stack_sdev, stack_ivar, stack_npix, stack_mask, stack_sres, stack_covar):
        """
        Construct :attr:`hdu` that is held in memory for manipulation of
        the object.  See :func:`construct_3d_hdu` if you want to convert
        the object into a DRP-like datacube.

        bin_indx is 2d

        stack_sres is expected to be a MaskedArray

        """
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, 'Constructing hdu ...')
        if stack_sres is None:
            stack_sres = numpy.ma.MaskedArray([self.drpf['SPECRES'].data]*self.nbins)

        self.covariance = None if (stack_covar is None or not isinstance(stack_covar, Covariance)) \
                                    else stack_covar.copy()
#        print(self.covariance.shape)
#        self.covariance.show(plane=self.covariance.input_indx[0])

        # Initialize the headers
        pri_hdr = self._initialize_primary_header()
        pri_hdr = self._add_method_header(pri_hdr)
        pri_hdr = self._add_reddening_header(pri_hdr)
        red_hdr = fits.Header()
        red_hdr = self._add_reddening_header(red_hdr)
        map_hdr = DAPFitsUtil.build_map_header(self.drpf, 'K Westfall <westfall@ucolick.org>')

        # Get the spatial map mask
        # Marginalize the DRP spectra over wavelength
        map_mask = DAPFitsUtil.marginalize_mask(self.drpf['MASK'].data,
                                                [ 'NOCOV', 'LOWCOV', 'DEADFIBER', 'FORESTAR',
                                                  'DONOTUSE' ], self.drpf.bitmask, self.bitmask,
                                                out_flag='DIDNOTUSE')
        map_mask = DAPFitsUtil.marginalize_mask(self.drpf['MASK'].data, ['FORESTAR'],
                                                self.drpf.bitmask, self.bitmask, out_mask=map_mask)
        drp_bad = (map_mask > 0)
        # Add the spectra with low spectral coverage
        indx = numpy.invert(good_fgoodpix.reshape(self.spatial_shape)) & numpy.invert(drp_bad)
        map_mask[indx] = self.bitmask.turn_on(map_mask[indx], 'LOW_SPECCOV')
        # Add the spectra with low S/N
        indx = numpy.invert(good_snr.reshape(self.spatial_shape)) & numpy.invert(drp_bad)
        map_mask[indx] = self.bitmask.turn_on(map_mask[indx], 'LOW_SNR')

        # Fill the covariance HDUs
        if self.covariance is not None:
            pri_hdr, ivar_hdu, covar_hdu = self.covariance.output_hdus(hdr=pri_hdr)

        # Save the data to the hdu attribute
        self.hdu = fits.HDUList([ fits.PrimaryHDU(header=pri_hdr),
                                  fits.ImageHDU(data=stack_flux.data, name='FLUX'),
                                  fits.ImageHDU(data=stack_ivar.data, name='IVAR'),
                                  fits.ImageHDU(data=stack_mask, name='MASK'),
                                  self.drpf['WAVE'].copy(),
                                  fits.ImageHDU(data=stack_sres.data, name='SPECRES'),
                                  fits.ImageHDU(data=self.galext.redcorr.data, header=red_hdr,
                                                name='REDCORR'),
                                  fits.ImageHDU(data=stack_sdev.data, name='FLUXD'),
                                  fits.ImageHDU(data=stack_npix.data, name='NPIX'),
                                  fits.ImageHDU(data=bin_indx, header=map_hdr, name='BINID'),
                                  fits.ImageHDU(data=map_mask, header=map_hdr, name='MAPMASK'),
                                  fits.BinTableHDU.from_columns( [ fits.Column(name=n,
                                                             format=rec_to_fits_type(bin_data[n]),
                                                array=bin_data[n]) for n in bin_data.dtype.names ],
                                                               name='BINS')
                                ])

        # Fill the covariance matrix
        if self.covariance is not None:
            self.hdu += [ covar_hdu ]


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


    @property
    def is_unbinned(self):
        """Determine if the spectra are unbinned."""
        return self.method['binpar'] is None and self.method['binclass'] is None \
                    and self.method['binfunc'] is None


    def do_not_fit_flags(self):
        return ['DIDNOTUSE', 'FORESTAR', 'LOW_SPECCOV', 'LOW_SNR', 'NONE_IN_STACK']


    def check_fgoodpix(self, minimum_fraction=0.8):
        """
        Determine which spectra in :attr:`rdxqa` have a fractional
        spectral coverage of greater than the provided minimum fraction.
        Only these spectra will be included in the binning.
    
        Args:
            minimum_fraction (float): (**Optional**) The minimum
                fraction of the spectrum that must be valid for the
                spectrum to be included in any bin.  Default is 0.8.

        Returns:
            numpy.ndarray : Boolean array for the spectra that satisfy
            the criterion.
        """
        return self.rdxqa['SPECTRUM'].data['FGOODPIX'] > minimum_fraction


    def above_snr_limit(self, sn_limit):
        """
        Flag bins above a provided S/N limit.
        """
#        warnings.warn('You\'re setting all but two spectra as bad!')
#        test = self.hdu['BINS'].data['SNR'] > sn_limit
#        test[:-2] = False
#        return test
        return self.hdu['BINS'].data['SNR'] > sn_limit


    @staticmethod
    def spectral_resolution_options():
        """
        Return the allowed options for treating the spectral resolution.

        Options are:

            'spaxel': If available, use the spectral resolution
            determined for each spaxel.  This is pulled from the 'DISP'
            extension in the DRP file; an exception will be raised if
            this extension does not exist!

            'cube': Only consider the median spectral resolution
            determine for the entire datacube.  This is pulled from the
            'SPECRES' extension in the DRP file; an exception will be
            raised if this extension does not exist!

        Returns:
            list : List of the available method keywords.
        """
        return ['spaxel', 'cube']


    def bin_spectra(self, drpf, rdxqa, reff=None, dapver=None, analysis_path=None,
                    directory_path=None, output_file=None, hardcopy=True, symlink_dir=None,
                    clobber=False, loggers=None, quiet=False):
        """
        Bin and stack the spectra.

        Args:
            drpf (:class:`mangadap.drpfits.DRPFits`): The DRP datacube
                with the spectra to bin.
            rdxqa (:class:`mangadap.proc.reductionassessments.ReductionAssessments`):
                The basic assessments of the DRP data that are needed
                for the binning procedures.
            reff (float): (**Optional**) The effective radius of the
                galaxy.
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
            symlink_dir (str): (**Optional**) Create a symbolic link to
                the created file in the supplied directory.  Default is
                to produce no symbolic link.
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

        # DRPFits object always needed
        if drpf is None:
            raise ValueError('Must provide DRP file object to compute assessments.')
        if not isinstance(drpf, DRPFits):
            raise TypeError('Must provide a valid DRPFits object!')
        if drpf.hdu is None:
            if not self.quiet:
                warnings.warn('DRP file previously unopened.  Reading now.')
            drpf.open_hdu()

        # Test if the RSS file exists; cannot compute covariance if not
        # TODO: this will break for a different stacking class
        if self.method['stackpar']['covar_mode'] is not 'none' and not drpf.can_compute_covariance:
            warnings.warn('Cannot determine covariance matrix!  Continuing without!')
            self.method['stackpar']['covar_mode'] = 'none'

        # TODO: How many of these attributes should I keep, vs. just use
        # drpf attributes?
        self.drpf = drpf
        self.shape = self.drpf.shape
        self.spatial_shape = self.drpf.spatial_shape
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

        # Set the Galactic extinction correction defined by the method
        self.galext = GalacticExtinction(form=self.method['galactic_reddening'],
                                         wave=drpf['WAVE'].data,
                                         ebv=drpf['PRIMARY'].header['EBVGAL'],
                                         rv=self.method['galactic_rv'])

        # Save the effective radius if provided.  Only used if/when
        # scaling the radii by the effective radius in the radial
        # binning approach
        if reff is not None:
            self.reff = reff

        #---------------------------------------------------------------
        # Get the good spectra
        #   - Must have valid pixels over more than 80% of the spectral
        #   range 
        good_fgoodpix = self.check_fgoodpix()
        #   - Must have sufficienct S/N, as defined by the input par
        good_snr = self._check_snr()
        # Good spectra have both these criteria
        good_spec = good_fgoodpix & good_snr

        # Report
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, '-'*50)
            log_output(loggers, 1, logging.INFO, '{0:^50}'.format('SPATIALLY BINNING SPECTRA'))
            log_output(self.loggers, 1, logging.INFO, '-'*50)

            log_output(self.loggers, 1, logging.INFO, 'Total spectra: {0}'.format(
                                                                            len(good_fgoodpix)))
            log_output(self.loggers, 1, logging.INFO, 'With 80% spectral coverage: {0}'.format(
                                                                        numpy.sum(good_fgoodpix)))
            log_output(self.loggers, 1, logging.INFO, 'With good S/N: {0}'.format(
                                                                            numpy.sum(good_snr)))
            log_output(self.loggers, 1, logging.INFO, 'Number of good spectra: {0}'.format(
                                                                            numpy.sum(good_spec)))
            if self.galext.form is None:
                log_output(self.loggers, 1, logging.INFO, 'Galactic dereddening not applied.')
            else:
                log_output(self.loggers, 1, logging.INFO,
                           'Galactic extinction law: {0}'.format(self.galext.form))
                log_output(self.loggers, 1, logging.INFO,
                           'Galactic E(B-V) = {0}'.format(self.galext.ebv))
                log_output(self.loggers, 1, logging.INFO,
                           'Galactic R(V): {0}'.format(self.galext.rv))

        if numpy.sum(good_spec) == 0:
            raise ValueError('No good spectra!')

        #---------------------------------------------------------------
        # Fill in any remaining binning parameters
        self._fill_method_par(good_spec)

        # (Re)Set the output paths
        self._set_paths(directory_path, dapver, analysis_path, output_file)

        #---------------------------------------------------------------
        # No binning so just fill the hdu with the appropriate data from
        # the DRP file
        if self.is_unbinned:

            # TODO: This is just a short cut.  Would be better if I
            # didn't use this.  I.e., I should treat the function calls
            # exactly the same, but just create binning and stacking
            # classes for when the data is unbinned...
        
            # Report data is unbinned
            if not self.quiet:
                log_output(self.loggers, 1, logging.INFO,
                            'No binning requested; analyzing DRP spaxel data directly.')

            # Generate pseudo bin index
            bin_indx = numpy.full(self.drpf.spatial_shape, -1, dtype=numpy.int)
            i = numpy.asarray(tuple(drpf.spatial_index[good_spec]))
            bin_indx[i[:,0],i[:,1]] = numpy.arange(numpy.sum(good_spec))

#            pyplot.imshow(bin_indx.reshape(self.drpf.spatial_shape), origin='lower',
#                          interpolation='nearest')
#            pyplot.show()

            # Build the data arrays directly from the DRP file
            flux = drpf.copy_to_masked_array(flag=drpf.do_not_stack_flags())[good_spec,:]

            sdev = numpy.ma.zeros(flux.shape, dtype=float)
            ivar = drpf.copy_to_masked_array(ext='IVAR',
                                             flag=drpf.do_not_stack_flags())[good_spec,:]
            npix = numpy.ones(flux.shape, dtype=numpy.int)
            npix[numpy.ma.getmaskarray(flux)] = 0
            mask = drpf.copy_to_array(ext='MASK')[good_spec,:]

            # The definition of the extension below:
            # - If the user wants the single resolution vector for the
            #   entire cube (spec_res = cube), then that extension is
            #   defined explicitly
            # - Currently the only other options is spec_res=spaxel.  In
            #   that case, the extension is set to None so that the
            #   DRPFits class can properly decide which data to return
            #   based on whether or not the DISP extension is present.
            # For MPL-5/DR14 data and earlier, the two spectral
            # resolution options should result in identical output!
            # !! Use the new pre-pixelized LSF measurements !!
            specres_ext='SPECRES' if self.method['spec_res'] == 'cube' else None
            sres = self.drpf.spectral_resolution(ext=specres_ext, toarray=True,
                                                 pre=self.method['prepixel_sres'],
                                                 fill=True)[good_spec,:]

            # TODO: Does this work with any covariance mode?  Don't like
            # this back and forth between what is supposed to be a stack
            # only function
            # (SpectralStack.build_covariance_data_DRPFits) and
            # something that is supposed to be indepenent of the details
            # of the stacking implementation (SpatiallyBinnedSpectra)...

            # Try to add the covariance data, if requested
            covar=None
            if self.method['stackpar']['covar_mode'] not in ['none', None]:
                if not self.quiet:
                    log_output(self.loggers, 1, logging.INFO, 'Attempting to compute covariance.')

                try:
                    covar = self.method['stackclass'].build_covariance_data_DRPFits(self.drpf,
                                                        self.method['stackpar']['covar_mode'], 
                                                        self.method['stackpar']['covar_par'])
                except AttributeError as e:
                    if not self.quiet:
                        warnings.warn('Could not build covariance data:: '
                                      'AttributeError: {0}'.format(e))
                    covar = None

                if isinstance(covar, Covariance):
                    covar = covar.spaxel_to_bin_covariance(bin_indx.ravel())

            # Fill the table with the per-spectrum data
            bin_data = self._unbinned_data_table(bin_indx.ravel())

            # Deredden the spectra if the method requests it
            flux, ivar, _, covar = self._apply_reddening(flux, ivar, None, covar)

            # Build the internal HDUList object and covariance attribute
            if (hardcopy or symlink_dir is not None) and not self.quiet:
                warnings.warn('Hardcopies and symlinks of the SpatiallyBinnedSpectra object are '
                              'never kept when analyzing unbinned data.')
            self.hardcopy = False
            self.symlink_dir = None
            self._construct_2d_hdu(bin_indx, good_fgoodpix, good_snr, bin_data, flux, sdev, ivar,
                                   npix, mask, sres, covar)
            # Finish
            if not self.quiet:
                log_output(self.loggers, 1, logging.INFO, '-'*50)
            return

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
        self.symlink_dir = symlink_dir
        if os.path.isfile(ofile) and not clobber:
            self.hardcopy = True
            if not self.quiet:
                log_output(self.loggers, 1, logging.INFO, 'Reading existing file')
            self.read(checksum=self.checksum)
            if not self.quiet and reff is not None and self.reff != reff:
                warnings.warn('Provided effective radius different from available file; set ' \
                              'clobber=True to overwrite.')
            # Make sure the symlink exists
            if self.symlink_dir is not None:
                create_symlink(ofile, self.symlink_dir, clobber=clobber)
            if not self.quiet:
                log_output(self.loggers, 1, logging.INFO, '-'*50)
            return

        #---------------------------------------------------------------
        # Determine how to bin the spectra.  To be included in any bin,
        # the spectra must be selected as 'good' according to the
        # selections made above.
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, 'Binning spectra ...')
        bin_indx = numpy.full(self.nspec, -1, dtype=numpy.int)
        bin_indx[good_spec] \
                = self.method['binfunc'](self.rdxqa['SPECTRUM'].data['SKY_COO'][good_spec,0],
                                         self.rdxqa['SPECTRUM'].data['SKY_COO'][good_spec,1],
                                         par=self.method['binpar'])
        if numpy.sum(bin_indx > -1) == 0:
            raise ValueError('No spectra in ANY bin!')

#        # Done for testing missing bins
#        warnings.warn('You\'re forcing bins 2 and 3 to be empty!')
#        time.sleep(3)
#        warnings.warn('Proceeding...')
#        bin_indx[ (bin_indx == 2) | (bin_indx == 3) ] = 1

        #---------------------------------------------------------------
        # Stack the spectra in each bin
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, 'Stacking spectra ...')
        stack_wave, stack_flux, stack_sdev, stack_npix, stack_ivar, stack_sres, stack_covar = \
                self.method['stackfunc'](self.drpf, bin_indx, par=self.method['stackpar'])
        if stack_covar is not None and stack_covar.is_correlation:
            raise ValueError('Unexpected to have stack covariance be a correlation matrix.')

        # Construct the mask
        stack_mask = numpy.zeros(stack_flux.shape, dtype=self.bitmask.minimum_dtype())
        indx = numpy.invert(stack_npix>0)
        stack_mask[indx] = self.bitmask.turn_on(stack_mask[indx], ['NONE_IN_STACK', 'NO_STDDEV'])
        indx = numpy.invert(stack_npix>1)
        stack_mask[indx] = self.bitmask.turn_on(stack_mask[indx], 'NO_STDDEV')

        #---------------------------------------------------------------
        # Fill the table with the per-bin data (BEFORE applying the
        # reddening)
        bin_data = self._binned_data_table(bin_indx, stack_flux, stack_ivar)

        #---------------------------------------------------------------
        # Deredden the spectra if the method requests it
        stack_flux, stack_ivar, stack_sdev, stack_covar \
                = self._apply_reddening(stack_flux, stack_ivar, stack_sdev, stack_covar)

#        if self.galext.form is not None:
#            stack_flux, stack_ivar = self.galext.apply(stack_flux, ivar=stack_ivar)
#            stack_sdev = self.galext.apply(stack_sdev)
#            if stack_covar is not None:
#                for i,j in enumerate(stack_covar.input_indx):
#                    stack_covar.cov[i] *= numpy.square(self.galext.redcorr[j])

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

        #---------------------------------------------------------------
        # Build the internal HDUList object and covariance attribute
        self.hardcopy = hardcopy
        self._construct_2d_hdu(bin_indx.reshape(self.spatial_shape), good_fgoodpix, good_snr,
                               bin_data, stack_flux, stack_sdev, stack_ivar, stack_npix,
                               stack_mask, stack_sres, stack_covar)

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
        Reformat the binned spectra into a cube matching the shape of
        the DRP fits file from which it was derived.
        """

        # Report
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, 'Constructing binned spectra datacube ...')

        # Get/Copy the necessary data arrays
        bin_indx = self.hdu['BINID'].data.copy()
        stack_flux = self.hdu['FLUX'].data.copy()
        stack_sdev = self.hdu['FLUXD'].data.copy()
        stack_ivar = self.hdu['IVAR'].data.copy()
        stack_npix = self.hdu['NPIX'].data.copy()
        stack_sres = self.hdu['SPECRES'].data.copy()

        # Reconstruct the stacked spectra into a DRP-like datacube
        flux, mask, sdev, ivar, npix, sres \
                = DAPFitsUtil.reconstruct_cube(self.shape, bin_indx.ravel(),
                                               [ stack_flux, numpy.ma.getmaskarray(stack_flux),
                                                 stack_sdev, stack_ivar, stack_npix, stack_sres ])

        # TODO: Move this to the Covariance class?
        if self.covariance is not None:
            covariance = self.covariance.bin_to_spaxel_covariance(bin_indx.ravel())

#        print(self.covariance.input_indx)
#        self.covariance.show(plane=self.covariance.input_indx[0])

#        pyplot.imshow(flux.reshape(self.drpf['FLUX'].shape)[:,:,1000].T, origin='lower',
#                      interpolation='nearest')
#        pyplot.show()
#        pyplot.step(stack_wave, flux.reshape(self.drpf['FLUX'].shape)[20,20,:], where='mid',
#                    linestyle='-', color='b')
#        pyplot.show()

        # Initialize the basics of the mask
        mask_bit_values = self._initialize_cube_mask()

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

        # Primary header is identical regardless of the shape of the
        # extensions
        pri_hdr = self.hdu['PRIMARY'].header.copy()
        cube_hdr = DAPFitsUtil.build_cube_header(self.drpf, 'K Westfall <westfall@ucolick.org>')

        # Fill the covariance HDUs
        if self.covariance is not None:
            pri_hdr, ivar_hdu, covar_hdu = covariance.output_hdus(reshape=True, hdr=pri_hdr)

        # Save the data to the hdu attribute
        # TODO: Strip headers of sub extensions of everything but
        # the WCS and HDUCLAS information.  Key WCS information in
        # BINID.
        hdu = fits.HDUList([ fits.PrimaryHDU(header=pri_hdr),
                             fits.ImageHDU(data=flux, header=cube_hdr, name='FLUX'),
                             fits.ImageHDU(data=ivar, header=cube_hdr, name='IVAR'),
                             fits.ImageHDU(data=mask_bit_values, header=cube_hdr, name='MASK'),
                             self.hdu['WAVE'].copy(),
                             fits.ImageHDU(data=sres, header=cube_hdr, name='SPECRES'),
                             self.hdu['REDCORR'].copy(),
                             fits.ImageHDU(data=sdev, header=cube_hdr, name='FLUXD'),
                             fits.ImageHDU(data=npix, header=cube_hdr, name='NPIX'),
                             self.hdu['BINID'].copy(),
                             self.hdu['BINS'].copy()
                           ])

        # Fill the covariance matrix
        if self.covariance is not None:
            self.hdu += [ covar_hdu ]

        return hdu


    def write(self, match_DRP=False, clobber=False):
        """
        Write the hdu object to the file.
        """
        # Convert the spectral arrays in the HDU to a 3D cube and write
        # it
        if match_DRP:
            hdu = self.construct_3d_hdu()
            DAPFitsUtil.write(hdu, self.file_path(), clobber=clobber, checksum=True,
                              symlink_dir=self.symlink_dir, loggers=self.loggers, quiet=self.quiet)
            return
        # Just write the unique (2D) data
        DAPFitsUtil.write(self.hdu, self.file_path(), clobber=clobber, checksum=True,
                          symlink_dir=self.symlink_dir, loggers=self.loggers, quiet=self.quiet) 


    def read(self, ifile=None, strict=True, checksum=False):
        """

        Read an existing file with a previously binned set of spectra.

        Args:
            ifile (str): (**Optional**) Name of the file with the data.
                Default is to use the name provided by
                :func:`file_path`.
            strict (bool): (**Optional**) Force a strict reading of the
                file to make sure that it adheres to the expected
                format.  At the moment, this only checks to make sure
                the method keyword printed in the header matches the
                expected value in :attr:`method`.
            checksum (bool): (**Optional**) Use the checksum in the
                header to check for corruption of the data.  Default is
                False.

        Raises:
            FileNotFoundError: Raised if the file does not exist.
            ValueError: Raised if strict is true and the header keyword
                'BINKEY' does not match the method keyword.

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
        if self.hdu['PRIMARY'].header['BINKEY'] != self.method['key']:
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

        # Attempt to read and construct the covariance object
        # TODO: Covariance matrices need to have the input_index
        # written to the file
#        self.covariance = Covariance.from_fits(self.hdu, ivar_ext=None, row_major=True,
#                                               correlation=True)
#        var = numpy.ma.power(self.hdu['IVAR'].data[:,:,self.covariance.input_index],
#                             -1).filled(0.0).reshape(-1, self.covariance.shape[-1])
#        self.covariance = self.covariance.apply_new_variance(var)

        try:
            self.covariance = Covariance.from_fits(self.hdu, ivar_ext=None, row_major=True,
                                                   correlation=True)
            var = numpy.ma.power(self.hdu['IVAR'].data[:,:,self.covariance.input_index],
                                 -1).filled(0.0).reshape(-1, self.covariance.shape[-1])
            self.covariance = self.covariance.apply_new_variance(var)
        except Exception as e:
            if not self.quiet:
                print(e)
                warnings.warn('Unable to find/read covariance data.')
            self.covariance = None

        # Attempt to read the effective radius
        try:
            self.reff = self.hdu['PRIMARY'].header['REFF']
        except:
            if not self.quiet:
                warnings.warn('Unable to read effective radius from file header!')
            self.reff = None

        # Attempt to read binning parameters
        if self.method['binpar'] is not None and callable(self.method['binpar'].fromheader):
            self.method['binpar'].fromheader(self.hdu['PRIMARY'].header)

        # Attempt to stacking parameters
        if self.method['stackpar'] is not None and callable(self.method['stackpar'].fromheader):
            self.method['stackpar'].fromheader(self.hdu['PRIMARY'].header)

        self.nbins = self.hdu['PRIMARY'].header['NBINS']
        self.missing_bins = SpatiallyBinnedSpectra._get_missing_bins(self.hdu['BINS'].data['BINID'])
#        try:
#            self.missing_bins = eval(self.hdu['PRIMARY'].header['EMPTYBIN'])
#        except KeyError:
#            # Assume if this fails, it's because the keyword doesn't
#            # exist
#            self.missing_bins = []

        # Galactic extinction data already set by bin_spectra


    def copy_to_array(self, ext='FLUX', waverange=None, include_missing=False):
        r"""
        Wrapper for :func:`mangadap.util.fitsutil.DAPFitsUtil.copy_to_array`
        specific for :class:`SpatiallyBinnedSpectra`.

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
                                        missing_bins=self.missing_bins if include_missing else None,
                                         nbins=self.nbins,
                                        unique_bins=DAPFitsUtil.unique_bins(self.hdu['BINID'].data))


    def copy_to_masked_array(self, ext='FLUX', flag=None, waverange=None, include_missing=False):
        """
        Wrapper for
        :func:`mangadap.util.fitsutil.DAPFitsUtil.copy_to_masked_array`
        specific for :class:`SpatiallyBinnedSpectra`.

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
                                       missing_bins=self.missing_bins if include_missing else None,
                                              nbins=self.nbins,
                                        unique_bins=DAPFitsUtil.unique_bins(self.hdu['BINID'].data))


    def get_bin_indices(self, bins):
        """
        Return the indices of the bins in the BIN table.

        Args:
            bins (array-like): The bin ID numbers to find
        
        Returns:
            numpy.ndarray: Integer array with the index of each bin ID in
            the BINID columns of the BINS extension.
        """
        return numpy.array([numpy.where(self.hdu['BINS'].data['BINID'] == b)[0][0] for b in bins])


    def find_nearest_bin(self, input_bins, weighted=False, indices=False):
        """
        Use a KDTree to find the bins nearest to and excluding a list of
        input bins.

        Args:
            input_bins (array-like): One or more bin ID numbers to use
                to locate the bin nearest to it based on its on-sky
                coordinates.  The list must be a unique set.
            weighted (bool): (**Optional**) Use the weighted coordinates
                (LW_SKY_COO) instead of the unweighted coordinates
                (SKY_COO) when finding the nearest bin.
            indices (bool): (**Optional**) Return the indices of the
                nearest bins instead of their ID number (default).
        Returns:
            numpy.ndarray: The bin IDs, one per input bin, of the bin
            closest to each input bin.  Any bins that are in
            :attr:`missing_bins` have a return value of -1; there are no
            coordinates for these bins.

        Raises:
            ValueError: Raised if the set of input bins is not unique or
                contains bin IDs that are not present.

        """
        # Check the input bin IDs
        _input_bins = numpy.atleast_1d(input_bins)
        single_value = len(_input_bins) == 1 and not isinstance(input_bins, (list, numpy.ndarray))
        if len(numpy.unique(_input_bins)) != len(_input_bins):
            raise ValueError('Must provide a unique set of input bin IDs.')
        maxbinid = numpy.amax(self.hdu['BINS'].data['BINID'])
        if numpy.amin(_input_bins) < 0 or numpy.amax(_input_bins) > maxbinid:
            return ValueError('Input contains invalid bin IDs.')

        # Account for any missing bins
        valid_input_bin = numpy.ones(_input_bins.size, dtype=bool) if len(self.missing_bins) == 0 \
                            else numpy.invert([ ib in self.missing_bins for ib in _input_bins ])
        if not numpy.all(valid_input_bin):
            warnings.warn('Input bin IDs include missing bins.  Returned as -1.')
        _input_bins = _input_bins[valid_input_bin]

        # Select the coordinate column
        col = 'LW_SKY_COO' if weighted else 'SKY_COO'

        # The reference bins are all bins EXCEPT those provided
        reference_bin = numpy.in1d(self.hdu['BINS'].data['BINID'], _input_bins, assume_unique=True,
                                   invert=True)

        # Get the coordinates of the reference grid
        ref_coo = numpy.array([ self.hdu['BINS'].data[col][reference_bin,0],
                                self.hdu['BINS'].data[col][reference_bin,1] ]).T

        # Construct the KDTree
        kd = spatial.KDTree(ref_coo)

        # Get the indices of the input bins in the internal list
        input_bin_indx = self.get_bin_indices(_input_bins)

        # Get the coordinates of the bins to match to the nearest
        # reference bin
        input_coo = numpy.array([ self.hdu['BINS'].data[col][input_bin_indx,0],
                                  self.hdu['BINS'].data[col][input_bin_indx,1] ]).T

        # Get the indices of the nearest bins
        dist, nearest_bin = kd.query(input_coo)

        # Return either the bin indices or the bin ID numbers
        output_bins = numpy.zeros(_input_bins.size, dtype=int) - 1
        output_bins[valid_input_bin] = numpy.arange(self.nbins)[reference_bin][nearest_bin] \
                        if indices else self.hdu['BINS'].data['BINID'][reference_bin][nearest_bin]

        return output_bins[0] if single_value else output_bins


    def replace_with_data_from_nearest_bin(self, data, bad_bins):
        """
        Replace data in the list of provided bad bins with the data from
        the nearest good bin.

        Args:
            data (array-like): Data for each bin.  The length must be
                the same as :attr:`nbins`.
            bad_bins (array-like): The list of indices (must not be a
                boolean array) with bad values to be replaced.
        
        Returns:
            numpy.ndarray: A new array with the bad data filled with the
            data from the nearest bin.

        Raises:
            ValueError: Raised if the input array doesn't have the
                correct shape or if the list of bad bins has numbers
                outside the viable range (0,self.nbins-1).
        """
        if len(data) != self.nbins:
            raise ValueError('Input data must have {0} elements.'.format(self.nbins))
            
        # No bad bins so just return the input
        if len(bad_bins) == 0:
            return data

        # Select the bins to replace
        replace = numpy.zeros(self.nbins, dtype=bool)
        replace[self.get_bin_indices(bad_bins)] = True

        # Get the coordinates of the reference grid
        coo = numpy.array([ self.hdu['BINS'].data['SKY_COO'][:,0],
                            self.hdu['BINS'].data['SKY_COO'][:,1] ]).T

        # Return the replaced data
        return replace_with_data_from_nearest_coo(coo, data, replace)

