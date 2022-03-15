# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
A class hierarchy that perform that spatially bins a set of
two-dimensional data.

The base class allows for user-defined definitions of binning
procedures.

----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import os
import glob
import warnings
import time
import logging

from IPython import embed

import numpy

from scipy import sparse, spatial, interpolate

from astropy.io import fits

from ..datacube import DataCube
from ..datacube import MaNGADataCube, MUSEDataCube
from ..par.parset import KeywordParSet, ParSet
from ..util.datatable import DataTable
from ..util.fitsutil import DAPFitsUtil
from ..util.fileio import create_symlink
from ..util.parser import DefaultConfig
from ..util.dapbitmask import DAPBitMask
from ..util.pixelmask import SpectralPixelMask
from ..util.sampling import angstroms_per_pixel
from ..util.covariance import Covariance
from ..util.geometry import SemiMajorAxisCoo
from ..util.extinction import GalacticExtinction
from ..util.log import log_output
from ..config import defaults
from . import spatialbinning
from .reductionassessments import ReductionAssessment
from .spectralstack import SpectralStackPar, SpectralStack
from .util import select_proc_method, replace_with_data_from_nearest_coo

from matplotlib import pyplot


class SpatiallyBinnedSpectraDef(KeywordParSet):
    """
    A class that holds the two parameter sets and the key designator
    for the binning scheme.

    The provided ``binfunc`` method must have a call with the
    following form::
        
        id = bin(x, y, par=par)

    where ``id`` is the index of the bin to place each spectrum and
    ``x`` and ``y`` are the on-sky Cartesian coordinates;
    e.g.::
        
        x = ReductionAssessment.hdu['SPECTRUM'].data['SKY_COO'][:,0]
        y = ReductionAssessment.hdu['SPECTRUM'].data['SKY_COO'][:,1]

    The provided ``stackfunc`` method must have a call with the
    following form::
        
        stack_wave, stack_flux, stack_sdev, stack_npix, stack_ivar, \
                stack_sres, stack_covar = stack(cube, id, par=par)

    where ``cube`` is a :class:`mangadap.datacube.datacube.DataCube`
    object. Note that the wavelengths are not returned because the
    input and output wavelength ranges are expected to be the same!

    As long as they are mutable, the values in ``par`` can change,
    meaning that some products of the bin algorithm can be passed to
    the stack algorithm. For example, if you want to weight the
    inclusion of the spectrum in the bin, you'd have to provide both
    the binning and stacking routines. Actually, if that's the case,
    you're better off providing a class object that will both bin and
    stack the spectra!

    The defined parameters are:

    .. include:: ../tables/spatiallybinnedspectradef.rst
    """
    def __init__(self, key=None, galactic_reddening=None, galactic_rv=None, minimum_snr=None,
                 minimum_frac=None, binpar=None, binclass=None, binfunc=None, stackpar=None,
                 stackclass=None, stackfunc=None):
        in_fl = [ int, float ]
        par_opt = [ ParSet, dict ]

        pars =     ['key', 'galactic_reddening', 'galactic_rv', 'minimum_snr', 'minimum_frac',
                    'binpar', 'binclass', 'binfunc', 'stackpar', 'stackclass', 'stackfunc']
        values =   [key, galactic_reddening, galactic_rv, minimum_snr, minimum_frac, binpar,
                    binclass, binfunc, stackpar, stackclass, stackfunc]
        options =  [None, None, None, None, None, None, None, None, None, None, None]
        dtypes =   [str, str, in_fl, in_fl, in_fl, par_opt, None, None, par_opt, None, None]
        can_call = [False, False, False, False, False, False, False, True, False, False, True]

        descr = ['Keyword used to distinguish between different spatial binning schemes.',
                 'The string identifier for the Galactic extinction curve to use.  See ' \
                    ':func:`mangadap.util.extinction.GalacticExtinction.valid_forms` for the ' \
                    'available curves.  Default is ``ODonnell``.',
                 'Ratio of V-band extinction to the B-V reddening.  Default is 3.1.',
                 'Minimum S/N of spectra to include in any bin.',
                 'Minimum fraction of unmasked pixels in each spectrum included in any bin.',
                 'The parameter set defining how to place each spectrum in a bin.',
                 'Instance of class object to use for the binning.  Needed in case binfunc ' \
                    'is a non-static member function of the class.',
                 'The function that determines which spectra go into each bin.',
                 'The parameter set defining how to stack the spectra in each bin.',
                 'Instance of class object to used to stack the spectra.  Needed in case ' \
                    'stackfunc is a non-static member function of the class.',
                 'The function that stacks the spectra in a given bin.']

        super(SpatiallyBinnedSpectraDef, self).__init__(pars, values=values, options=options,
                                                        dtypes=dtypes, can_call=can_call,
                                                        descr=descr)


def validate_spatial_binning_scheme_config(cnfg):
    """ 
    Validate the :class:`mangadap.util.parser.DefaultConfig` object with
    spatial-binning scheme parameters.

    Args:
        cnfg (:class:`mangadap.util.parser.DefaultConfig`):
            Object with the spatial-binning method parameters as
            needed by
            :class:`mangadap.proc.spatiallybinnedspectra.SpatiallyBinnedSpectraDef`.

    Raises:
        KeyError:
            Raised if any required keywords do not exist.
        ValueError:
            Raised if keys have unacceptable values.
    """

    # TODO: Can most of this stuff move to the relevant classes
    # themselves? Seems like the binning method should be agnostic
    # about the parameters needed by specific class objects.

    # Check for required keywords
    required_keywords = ['key', 'method']
    if not cnfg.all_required(required_keywords):
        raise KeyError('Keywords {0} must all have valid values.'.format(required_keywords))

    if cnfg['method'] not in [ 'none', 'global', 'voronoi', 'radial', 'square' ]:
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

    if cnfg['method'] == 'square' and not cnfg.keyword_specified('binsz'):
        raise KeyError('Keyword \'binsz\' must be provided for square binning.')



def available_spatial_binning_methods():
    """
    Return the list of available binning schemes.
    
    Returns:
        :obj:`list`: A list of
        :func:`mangadap.proc.spatiallybinnedspectra.SpatiallyBinnedSpectraDef`
        objects, each defining a separate binning method.

    Raises:
        IOError:
            Raised if no binning scheme configuration files could be
            found.
        KeyError:
            Raised if the binning method keywords are not all unique.

    .. todo::
        - Somehow add a python call that reads the databases and
          constructs the table for presentation in sphinx so that the
          text above doesn't have to be edited with changes in the
          available databases.
        
    """
    # Check the configuration files exist
    search_dir = os.path.join(defaults.dap_config_root(), 'spatial_binning')
    ini_files = glob.glob(os.path.join(search_dir, '*.ini'))
    if len(ini_files) == 0:
        raise IOError('Could not find any configuration files in {0} !'.format(search_dir))

    # TODO: Can most of this stuff move to the relevant classes
    # themselves? Seems like the binning method should be agnostic
    # about the parameters needed by specific class objects.

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
            binpar = spatialbinning.RadialBinningPar(center=center,
                                            pa=cnfg.getfloat('pa', default=0.),
                                            ell=cnfg.getfloat('ell', default=0.),
                                            radius_scale=cnfg.getfloat('radius_scale', default=1.),
                                            radii=radii,
                                            log_step=cnfg.getbool('log_step', default=False))
            binclass = spatialbinning.RadialBinning()
            binfunc = binclass.bin_index
        elif cnfg['method'] == 'voronoi':
            binpar = spatialbinning.VoronoiBinningPar(target_snr=cnfg.getfloat('target_snr'),
                                                      covar=cnfg.getfloat('noise_calib'))
            binclass = spatialbinning.VoronoiBinning()
            binfunc = binclass.bin_index
        elif cnfg['method'] == 'square':
            binpar = spatialbinning.SquareBinningPar(binsz=cnfg.getfloat('binsz', default=2.0))
            binclass = spatialbinning.SquareBinning()
            binfunc = binclass.bin_index
        else:   # Do not bin!
            binpar = None
            binclass = None
            binfunc = None

        stackpar = SpectralStackPar(operation=cnfg.get('operation', default='mean'),
                                    register=cnfg.getbool('velocity_register', default=False),
                                    covar_mode=cnfg.get('stack_covariance_mode', default='none'),
                                    covar_par=SpectralStack.parse_covariance_parameters(
                                            cnfg.get('stack_covariance_mode', default='none'),
                                            cnfg['stack_covariance_par']))
        stackclass = SpectralStack()
        stackfunc = stackclass.stack_datacube

        binning_methods += [SpatiallyBinnedSpectraDef(key=cnfg['key'],
                                        galactic_reddening=cnfg['galactic_reddening'],
                                        galactic_rv=cnfg.getfloat('galactic_rv', default=3.1),
                                        minimum_snr=cnfg.getfloat('minimum_snr', default=0.),
                                        minimum_frac=cnfg.getfloat('minimum_frac', default=0.8),
                                        binpar=binpar, binclass=binclass, binfunc=binfunc,
                                        stackpar=stackpar, stackclass=stackclass,
                                        stackfunc=stackfunc)]

    # Check the keywords of the libraries are all unique
    if len(numpy.unique(numpy.array([method['key'] for method in binning_methods]))) \
            != len(binning_methods):
        raise KeyError('Spatial binning method keywords are not all unique!')

    # Return the default list of assessment methods
    return binning_methods


class SpatiallyBinnedSpectraBitMask(DAPBitMask):
    r"""
    Derived class that specifies the mask bits for the spatially binned
    spectra.  The maskbits defined are:
    
    .. include:: ../tables/spatiallybinnedspectrabitmask.rst
    """
    cfg_root = 'spatially_binned_spectra_bits'


class SpatiallyBinnedSpectraDataTable(DataTable):
    """
    Primary data table with the results of the spatial binning.

    Table includes:

    .. include:: ../tables/spatiallybinnedspectradatatable.rst

    .. todo::

        There currently is no mask; however, could add masks for:
            - Insufficient good wavelength channels in the spectrum.
            - Variance in flux in bin is too large.

    Args:
        shape (:obj:`int`, :obj:`tuple`, optional):
            The shape of the initial array. If None, the data array
            will not be instantiated; use :func:`init` to initialize
            the data array after instantiation.
    """
    def __init__(self, shape=None):
        # NOTE: This should require python 3.7 to make sure that this
        # is an "ordered" dictionary.
        datamodel = dict(BINID=dict(typ=int, shape=None, descr='Bin ID number'),
                         NBIN=dict(typ=int, shape=None, descr='Number of spaxels in the bin'),
                         SKY_COO=dict(typ=float, shape=(2,),
                                      descr='The mean on-sky coordinates of the binned spaxels.'),
                         LW_SKY_COO=dict(typ=float, shape=(2,),
                                         descr='The luminosity-weighted mean on-sky coordinates '
                                               'of the binned spaxels.'),
                         ELL_COO=dict(typ=float, shape=(2,),
                                      descr='The mean elliptical coordinates of the binned '
                                            'spaxels.'),
                         LW_ELL_COO=dict(typ=float, shape=(2,),
                                         descr='The luminosity-weighted mean elliptical '
                                               'coordinates of the binned spaxels.'),
                         AREA=dict(typ=float, shape=None, descr='Total on-sky area of the bin'),
                         AREA_FRAC=dict(typ=float, shape=None,
                                        descr='Fraction of the expected area covered by good '
                                              'spaxels (not relevant to all binning schemes)'),
                         SIGNAL=dict(typ=float, shape=None,
                                     descr='Mean flux per pixel in the binned spectrum.'),
                         VARIANCE=dict(typ=float, shape=None,
                                       descr='Mean variance per pixel in the binned spectrum.'),
                         SNR=dict(typ=float, shape=None,
                                       descr='Mean S/N per pixel in the binned spectrum.'))

        keys = list(datamodel.keys())
        super(SpatiallyBinnedSpectraDataTable,
                self).__init__(keys, [datamodel[k]['typ'] for k in keys],
                               element_shapes=[datamodel[k]['shape'] for k in keys],
                               descr=[datamodel[k]['descr'] for k in keys],
                               shape=shape)



class SpatiallyBinnedSpectra:
    r"""
    Class that holds spatially binned spectra.

    Args:
        method_key (:obj:`str`):
            The keyword that designates which method, provided in
            ``method_list``, to use for the binning procedure.
        cube (:class:`mangadap.datacube.datacube.DataCube`):
            The datacube with the spectra to bin.
        rdxqa (:class:`mangadap.proc.reductionassessments.ReductionAssessments`):
            The basic assessments of the datacube that are used for
            the binning procedures.
        reff (:obj:`float`, optional):
            The effective radius of the galaxy in arcsec.
        method_list (:obj:`list`, optional):
            List of :class:`SpatiallyBinnedSpectraDef` objects that
            define one or more methods to use for the spatial
            binning. Default is to use the config files in the DAP
            source directory to construct the available methods using
            :func:`available_spatial_binning_methods`.
        dapver (:obj:`str`, optional):
            The DAP version. Used to construct the output paths,
            overriding the default defined by
            :func:`mangadap.config.defaults.dap_version`. Does
            **not** select the version of the code to use.
        analysis_path (:obj:`str`, optional):
            The top-level path for the DAP output files, used to
            override the default defined by
            :func:`mangadap.config.defaults.dap_analysis_path`.
        directory_path (:obj:`str`, optional):
            The exact path to the directory with DAP output that is
            common to the DAP "methods". Default is defined by
            :func:`mangadap.confgi.defaults.dap_common_path`.
        output_file (:obj:`str`, optional):
            Exact name for the output file. The default is to use
            :func:`mangadap.config.defaults.dap_file_name`.
        hardcopy (:obj:`bool`, optional):
            Flag to write the `astropy.io.fits.HDUList`_
            (:attr:`hdu`) to disk. If False, the object data is only
            kept in memory.
        symlink_dir (:obj:`str`, optional):
            Create a symbolic link to the created file in the
            supplied directory. If None, no symbolic link is created.
        clobber (:obj:`bool`, optional):
            Overwrite any existing files. If False, any existing
            files will be used. If True, the analysis is redone and
            any existing output is overwritten.
        checksum (:obj:`bool`, optional):
            Use the checksum in the fits header to confirm that the
            data has not been corrupted. The checksum is **always**
            written to the fits header when the file is created.
        loggers (:obj:`list`, optional):
            List of `logging.Logger`_ objects to log progress;
            ignored if quiet=True. Logging is done using
            :func:`mangadap.util.log.log_output`. If None, no logging
            is performed and output is just written to ``stdout``.
        quiet (:obj:`bool`, optional):
            Suppress all terminal and logging output.

    Attributes:
        loggers (:obj:`list`):
            List of `logging.Logger`_ objects to log progress;
            ignored if quiet=True. Logging is done using
            :func:`mangadap.util.log.log_output`.
        quiet (:obj:`bool`):
            Suppress all terminal and logging output.

    .. todo::
        - Allow velocity offsets for registration.
        - Fill in attributes.
   
    """
    def __init__(self, method_key, cube, rdxqa, reff=None, method_list=None, dapver=None,
                 analysis_path=None, directory_path=None, output_file=None, hardcopy=True,
                 symlink_dir=None, clobber=False, checksum=False, loggers=None, quiet=False):

        self.loggers = None
        self.quiet = False

        # Define the method properties
        self.method = self.define_method(method_key, method_list=method_list)

        self.cube = None
        self.rdxqa = None
        self.reff = None
        self.galext = None

        # Define the output directory and file
        self.directory_path = None      # Set in _set_paths
        self.output_file = None
        self.hardcopy = None
        self.symlink_dir = None

        # Define the bitmask
        self.bitmask = SpatiallyBinnedSpectraBitMask()

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
        self.nwave = None

        self.nbins = None
        self.missing_bins = None

        self.covariance = None

        # Bin the spectra
        self.bin_spectra(cube, rdxqa, reff=reff, dapver=dapver, analysis_path=analysis_path,
                         directory_path=directory_path, output_file=output_file, hardcopy=hardcopy,
                         symlink_dir=symlink_dir, clobber=clobber, loggers=loggers, quiet=quiet)

    def __getitem__(self, key):
        return self.hdu[key]

    @staticmethod
    def define_method(method_key, method_list=None):
        r"""
        Select the method
        """
        # Grab the specific method
        return select_proc_method(method_key, SpatiallyBinnedSpectraDef, method_list=method_list,
                                  available_func=available_spatial_binning_methods)

    def _fill_method_par(self, good_spec):
        """
        Finalize the binning parameters, as needed.

        **For the radial binning**, set the ellipticity and position
        angle.  Set scale radius to the effective radius if desired.

        **For the Voronoi binning**, set the signal and noise.

        .. todo::
            Abstract this so that this wrapper class doesn't need to
            know about the internal constraints of each binning
            approach.

        Args:
            good_spec (`numpy.ndarray`_):
                List of spectra to include in the binning. See
                :func:`check_fgoodpix` and :func:`_check_snr`.

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
                self.method['binpar']['covar'] = Covariance(
                                inp=sparse.coo_matrix((covar[covar > 0].ravel(),
                                                      (i[covar > 0].ravel(), j[covar > 0].ravel())),
                                                      shape=covar.shape).tocsr())
            else:
                self.method['binpar']['noise'] \
                        = numpy.sqrt(self.rdxqa['SPECTRUM'].data['VARIANCE'][good_spec])

        # Nothing to add for binning types 'none' or 'global', or
        # user-defined function!

    def _set_paths(self, directory_path, dapver, analysis_path, output_file):
        """
        Construct the main output paths.

        This method sets :attr:`directory_path` and
        :attr:`output_file`. If not provided as arguments, the
        defaults are set using, respectively,
        :func:`mangadap.config.defaults.dap_common_path` and
        :func:`mangadap.config.defaults.dap_file_name`.

        Args:
            dapver (:obj:`str`, optional):
                The DAP version. Used to construct the output paths,
                overriding the default defined by
                :func:`mangadap.config.defaults.dap_version`. Does
                **not** select the version of the code to use.
            analysis_path (:obj:`str`, optional):
                The top-level path for the DAP output files, used to
                override the default defined by
                :func:`mangadap.config.defaults.dap_analysis_path`.
            directory_path (:obj:`str`, optional):
                The exact path to the directory with DAP output that
                is common to the DAP "methods". Default is defined by
                :func:`mangadap.confgi.defaults.dap_common_path`.
            output_file (:obj:`str`, optional):
                Exact name for the output file. The default is to use
                :func:`mangadap.config.defaults.dap_file_name`.
        """
        # Set the output directory path
        self.directory_path = defaults.dap_common_path(plate=self.cube.plate,
                                                       ifudesign=self.cube.ifudesign,
                                                       drpver=self.cube.drpver, dapver=dapver,
                                                       analysis_path=analysis_path) \
                                        if directory_path is None else str(directory_path)

        # Set the output file
        method = '{0}-{1}'.format(self.rdxqa.method['key'], self.method['key'])
        self.output_file = defaults.dap_file_name(self.cube.plate, self.cube.ifudesign, method) \
                                        if output_file is None else str(output_file)

    def _initialize_primary_header(self, hdr=None):
        """
        Construct the primary header for the reference file.

        Args:
            hdr (`astropy.io.fits.Header`_, optional):
                Input base header for added keywords. If None, uses
                the :attr:`cube` header (if there is one) and then
                cleans the header using
                :func:`mangadap.util.fitsutil.DAPFitsUtil.clean_dap_primary_header`.

        Returns:
            `astropy.io.fits.Header`_: Initialized header object.
        """
        # Copy the from the DRP and clean it
        if hdr is None:
            hdr = self.cube.prihdr.copy()
            hdr = DAPFitsUtil.clean_dap_primary_header(hdr)
        hdr['AUTHOR'] = 'Kyle B. Westfall <westfall@ucolick.org>'
        if self.reff is not None:
            hdr['REFF'] = self.reff
        hdr['BINKEY'] = (self.method['key'], 'Spectal binning method keyword')
        hdr['BINMINSN'] = (self.method['minimum_snr'], 'Minimum S/N of spectrum to include')
        hdr['FSPCOV'] = (self.method['minimum_frac'], 'Minimum required fraction of good pixels')
        hdr['NBINS'] = (self.nbins, 'Number of unique spatial bins')
        return hdr

    def _add_method_header(self, hdr):
        """
        Add method-specific metadata to the header.

        Note that the header object is both edited in-place and
        returned.

        Args:
            hdr (`astropy.io.fits.Header`_):
                Input base header for added keywords. Expected to
                have been initialized using
                :func:`_initialize_primary_header`.

        Returns:
            `astropy.io.fits.Header`_: Edited header object.
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

        Note that the header object is both edited in-place and
        returned.

        Args:
            hdr (`astropy.io.fits.Header`_):
                Input header object to edit.

        Returns:
            `astropy.io.fits.Header`_: Edited header object.
        """
        hdr['EBVGAL'] = (self.cube.prihdr['EBVGAL'], 'Galactic reddening E(B-V)')
        hdr['GEXTLAW'] = (str(self.galext.form), 'Galactic extinction law')
        hdr['RVGAL'] = (self.galext.rv, 'Ratio of total to selective extinction, R(V)')
        return hdr

    def _finalize_cube_mask(self, mask):
        """
        Finalize the mask after the 2D mask has been reconstructed
        into a 3D cube.

        This mostly handles the masks for regions outside the
        datacube field of view.

        This propagates the bits flagging pixels that should not be
        included in the stacked spectra, and sets the the LOW_SPECCOV
        and LOW_SNR bits based on the results from
        :func:`check_fgoodpix` and :func:`_check_snr`.

        Note that the input mask is both edited in-place and
        returned.

        .. todo::
            This needs to be abstracted for non-DRP datacubes.

        Args:
            mask (`numpy.ndarray`_):
                3D array with the current bitmask data.

        Returns:
            `numpy.ndarray`_: Edited bitmask data.
        """
        # Get the spaxels with good spectral coverage and S/N
        good_fgoodpix = self.check_fgoodpix()
        good_snr = self._check_snr()

        ## KHRR added if statement:
#        if isinstance(self.cube, MaNGADataCube):
#            # Turn on the flag stating that the pixel wasn't used
#            indx = self.cube.bitmask.flagged(self.cube.mask, flag=self.cube.do_not_stack_flags())
#            mask[indx] = self.bitmask.turn_on(mask[indx], 'DIDNOTUSE')
#
#            # Turn on the flag stating that the pixel has a foreground star
#            indx = self.cube.bitmask.flagged(self.cube.mask, flag='FORESTAR')
#            mask[indx] = self.bitmask.turn_on(mask[indx], 'FORESTAR')
        if self.cube.bitmask is not None:

            # Turn on the flag stating that the pixel wasn't used
            indx = self.cube.bitmask.flagged(self.cube.mask, flag=self.cube.do_not_stack_flags())
            mask[indx] = self.bitmask.turn_on(mask[indx], 'DIDNOTUSE')

            # Turn on the flag stating that the pixel has a foreground star
            if self.cube.propagate_flags() is not None:
                flags = self.cube.propagate_flags()
                # TODO: This will barf if flags is a numpy array
                flags = flags if isinstance(flags, list) else [flags]
                for flag in flags:
                    indx = self.cube.bitmask.flagged(self.cube.mask, flag=flag)
                    if numpy.any(indx):
                        if flag in self.bitmask.keys():
                            mask[indx] = self.bitmask.turn_on(mask[indx], 'FORESTAR')
                        warnings.warn(f'{flag} is not a flag in {self.bitmask.__class__.__name__} '
                                      'and cannot be propagated!')


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
            `numpy.ndarray`_: Boolean array for the spectra that
            satisfy the criterion.
        """
        return self.rdxqa['SPECTRUM'].data['SNR'] > self.method['minimum_snr']

    def _assign_spectral_arrays(self):
        """
        Set :attr:`spectral_arrays`, which contains the list of
        extensions in :attr:`hdu` that contain spectral data.
        """
        self.spectral_arrays = ['FLUX', 'IVAR', 'MASK', 'SPECRES', 'FLUXD', 'NPIX']

    def _assign_image_arrays(self):
        """
        Set :attr:`image_arrays`, which contains the list of extensions
        in :attr:`hdu` that are on-sky image data.
        """
        self.image_arrays = ['BINID']

    @staticmethod
    def _get_missing_bins(unique_bins):
        """Return bin IDs omitted from a sequential list."""
        return list(set(numpy.arange(numpy.amax(unique_bins)+1)) - set(unique_bins))

    def _unbinned_data_table(self, bin_indx):
        """
        Construct the output data table for the unbinned spectra.

        Args:
            bin_indx (`numpy.ndarray`_):
                The integer vector with the bin associated with each
                spectrum in the DRP cube. This is the flattened
                ``BINID`` array.

        Returns:
            `numpy.recarray`_: The record array that is put in the
            BINS extension of :attr:`hdu`.
        """
        # Find the spectra that have a bin ID
        good_spec = bin_indx > -1
        self.nbins = numpy.amax(bin_indx)+1
        self.missing_bins = []                  # No missing bins

        # Copy the data from the ReductionAssessments object
        bin_data = SpatiallyBinnedSpectraDataTable(shape=self.nbins)

        bin_data['BINID'] = numpy.arange(self.nbins)
        bin_data['NBIN'] = numpy.ones(self.nbins, dtype=int)
        bin_data['SKY_COO'][:,0] = self.rdxqa['SPECTRUM'].data['SKY_COO'][good_spec,0]
        bin_data['SKY_COO'][:,1] = self.rdxqa['SPECTRUM'].data['SKY_COO'][good_spec,1]
        bin_data['LW_SKY_COO'] = bin_data['SKY_COO'].copy()
        bin_data['ELL_COO'][:,0] = self.rdxqa['SPECTRUM'].data['ELL_COO'][good_spec,0]
        bin_data['ELL_COO'][:,1] = self.rdxqa['SPECTRUM'].data['ELL_COO'][good_spec,1]
        bin_data['LW_ELL_COO'] = bin_data['ELL_COO'].copy()

        if self.cube.pixelscale is None:
            if self.cube.wcs is not None:
                self.cube._get_pixelscale()
            else:
                raise ValueError('Must provide datacube spatial pixelscale, or a WCS so that it '
                                 'can be calculated.')

        bin_data['AREA'] = numpy.full(self.nbins, numpy.square(self.cube.pixelscale), dtype=float)
        bin_data['AREA_FRAC'] = numpy.ones(self.nbins, dtype=int)
        bin_data['SIGNAL'] = self.rdxqa['SPECTRUM'].data['SIGNAL'][good_spec]
        bin_data['VARIANCE'] = self.rdxqa['SPECTRUM'].data['VARIANCE'][good_spec]
        bin_data['SNR'] = self.rdxqa['SPECTRUM'].data['SNR'][good_spec]
        return bin_data

    def _binned_data_table(self, bin_indx, stack_flux, stack_ivar, per_pixel=True):
        r"""
        Construct the output data table for the binned spectra.

        Args:
            bin_indx (`numpy.ndarray`_):
                The integer vector with the bin associated with each
                spectrum in the DRP cube. This is the flattened
                ``BINID`` array.
            stack_flux (`numpy.ndarray`_):
                The stacked spectra with shape :math:`N_{\rm
                spec}\times\N_\lambda}`.
            stack_ivar (`numpy.ndarray`_):
                The stacked inverse variance with shape :math:`N_{\rm
                spec}\times\N_\lambda}`.
            per_pixel (:obj:`bool`, optional):
                Base the flux statistics on per-pixel measurements.
                Set to False for a per-angstrom calculation.

        Returns:
            `numpy.recarray`_: The record array that is put in the
            ``BINS`` extension of :attr:`hdu`.
        """
        # Get the unique bins and the number of spectra in each bin
        # TODO: Below assumes that the first unique value is always -1!!

#        unique_bins, bin_count = map(lambda x : x[1:], numpy.unique(bin_indx, return_counts=True))
        unique_bins, bin_count = numpy.unique(bin_indx, return_counts=True)
        if unique_bins[0] == -1:
            unique_bins = unique_bins[1:]
            bin_count = bin_count[1:]

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
        bin_data = SpatiallyBinnedSpectraDataTable(shape=self.nbins)
        bin_data['BINID'] = unique_bins.astype(int)
        bin_data['NBIN'] = bin_count.astype(int)

        # Recalculate the mean signal, mean noise, and mean S/N of the
        # binned spectra
        if stack_ivar is None:
            warnings.warn('No inverse variance for stack.  Errors set to unity.')

        _wavelength_mask = SpectralPixelMask(waverange=self.rdxqa.method['waverange']
                                             ).boolean(self.cube.wave)
        _mask = numpy.ma.getmaskarray(stack_flux) | _wavelength_mask[None,:]

        if stack_ivar is not None:
            _mask |= numpy.ma.getmaskarray(stack_ivar)

        _stack_flux = numpy.ma.MaskedArray(stack_flux.data, mask=_mask)

        # TODO: Somehow consolidate with
        # mangadap.datacube.datacube.DataCube.flux_stats?
        # Set the response function
        dw = numpy.ones(self.cube.nwave, dtype=float) if per_pixel \
                else angstroms_per_pixel(self.cube.wave, log=self.cube.log)
        _response_func = self.cube.interpolate_to_match(self.rdxqa.method['response_func'])
        # Get the signal
        response_integral = numpy.sum(numpy.invert(numpy.ma.getmaskarray(_stack_flux))
                                        * (_response_func*dw)[None,:], axis=1)
        bin_data['SIGNAL'] = numpy.ma.divide(numpy.ma.sum(_stack_flux*(_response_func*dw)[None,:],
                                                          axis=1), response_integral).filled(0.0)
        # And the variance and SNR, if the inverse variance is available
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

        # Sort the list of bin ids and determine if any spectra jump
        # between bins
        srt = numpy.argsort(bin_indx)

        bin_change = numpy.where(numpy.diff(bin_indx[srt]) > 0)[0] + 1
        if bin_change.size != self.nbins:
            bin_change = numpy.append([0], bin_change)

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

        # Compute the area covered by each bin
        bin_data['AREA'] = self.cube.binned_on_sky_area(bin_indx)[1]

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

        return bin_data

    def _apply_reddening(self, flux, ivar, sdev, covar, deredden=True):
        """
        Correct the spectra for Galactic reddening.

        Largely a wrapper for executing
        :func:`mangadap.util.extinction.GalacticExtinction.apply` on
        the provided arrays.

        If the reddening law is undefined (:attr:`galext.form` is
        None), the method simply returns the input.

        Otherwise, the input arrays are both modified directly and
        returned.

        Args:
            flux (`numpy.ndarray`_):
                Flux array
            ivar (`numpy.ndarray`_):
                Inverse variance array. Can be None.
            sdev (`numpy.ndarray`_):
                The standard deviation of the stacked spectra, if
                relevant to ``flux``.  Can be None.
            covar (:class:`mangadap.util.covariance.Covariance`):
                Spatial covariance in the binned spectra. Assumed to
                be a 3D covariance cube. Can be None.

        Returns:
            :obj:`tuple`: Returns four `numpy.ndarray`_ objects with
            the reddening corrected flux, inverse variance, standard
            deviation in the stack, and spatial covariance. Any of
            the latter three can be None if the corresponding input
            is None.
        """
        if self.galext.form is None:
            return flux, ivar, sdev, covar

        # NOTE: I'm doing some gymnastics here to make sure that the
        # masked values are altered by the reddening, but that they're
        # still masked on output.

        _flux = flux.data if isinstance(flux, numpy.ma.MaskedArray) else flux.copy()
        _ivar = None if ivar is None \
                    else (ivar.data if isinstance(ivar, numpy.ma.MaskedArray) else ivar.copy())
        _sdev = None if sdev is None \
                    else (sdev.data if isinstance(sdev, numpy.ma.MaskedArray) else sdev.copy())
        _covar = None if covar is None else covar.copy()

        # Apply the reddening correction
        _flux, _ivar = self.galext.apply(_flux, ivar=_ivar, deredden=deredden)
        if sdev is not None:
            _sdev = self.galext.apply(_sdev, deredden=deredden)
        if covar is not None:
            for i,j in enumerate(_covar.input_indx):
                _covar.cov[i] = _covar.cov[i] * numpy.square(self.galext.redcorr[j]) if deredden \
                                    else _covar.cov[i] / numpy.square(self.galext.redcorr[j])

        if isinstance(flux, numpy.ma.MaskedArray):
            _flux = numpy.ma.MaskedArray(_flux, mask=numpy.ma.getmaskarray(flux).copy())
        if _ivar is not None and isinstance(ivar, numpy.ma.MaskedArray):
            _ivar = numpy.ma.MaskedArray(_ivar, mask=numpy.ma.getmaskarray(ivar).copy())
        if _sdev is not None and isinstance(_sdev, numpy.ma.MaskedArray):
            _sdev = numpy.ma.MaskedArray(_sdev, mask=numpy.ma.getmaskarray(sdev).copy())

        return _flux, _ivar, _sdev, _covar

    # TODO: Allow both stack_sres and self.cube.sres to be None.
    def _construct_2d_hdu(self, bin_indx, good_fgoodpix, good_snr, bin_data, stack_flux,
                          stack_sdev, stack_ivar, stack_npix, stack_mask, stack_sres, stack_covar):
        r"""
        Construct :attr:`hdu` that is held in memory for manipulation
        of the object. See :func:`construct_3d_hdu` to convert the
        object into a datacube.

        Args:
            bin_indx (`numpy.ndarray`_):
                2D array with the bin associated with each spaxel.
            good_fgoodpix (`numpy.ndarray`_):
                Boolean array selecting spaxels with good spectral
                coverage.
            good_snr (`numpy.ndarray`_):
                Boolean array selecting spaxels with good S/N.
            bin_data (:class:`SpatiallyBinnedSpectraDataTable`):
                Data table with relevant metadata for each binned
                spectrum.
            stack_flux (`numpy.ndarray`_):
                Array with the stacked spectral flux. Shape is
                :math:`(N_{\rm bin}, N_{\rm wave})`.
            stack_sdev (`numpy.ndarray`_):
                Array with the standard deviation in the stacked
                spectral flux. Shape is :math:`(N_{\rm bin}, N_{\rm
                wave})`.
            stack_ivar (`numpy.ndarray`_):
                Array with the inverse variance in the stacked
                spectral flux. Shape is :math:`(N_{\rm bin}, N_{\rm
                wave})`.
            stack_npix (`numpy.ndarray`_):
                Integer array with the number of spectral channels
                that were stacked to construct the binned spectrum.
                Shape is :math:`(N_{\rm bin}, N_{\rm wave})`.
            stack_mask (`numpy.ndarray`_):
                Integer array with the bitmasks associated with pixel
                in the stacked spectra. Shape is :math:`(N_{\rm bin},
                N_{\rm wave})`.
            stack_sres (`numpy.ndarray`_):
                Array with the spectral resolution of the stacked
                spectra. Shape is :math:`(N_{\rm bin}, N_{\rm
                wave})`. If None, the spectral resolution is set to
                be the median of the resolution in the datacube. If
                the datacube also has no spectral resolution data,
                the method faults.
            stack_covar (:class:`mangadap.util.covariance.Covariance`):
                Spatial covariance in the stacked spectra. Can be
                None.
        """
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, 'Constructing hdu ...')
        if stack_sres is None:
            if self.cube.sres is None:
                raise ValueError('Must be able to define spectral resolution.')
            stack_sres = numpy.ma.MaskedArray([numpy.median(self.cube.sres.reshape(-1, 
                                                            self.cube.nwave), axis=0)]*self.nbins)

        self.covariance = None if stack_covar is None or not isinstance(stack_covar, Covariance) \
                                    else stack_covar.copy()

        # Initialize the headers
        pri_hdr = self._initialize_primary_header()
        pri_hdr = self._add_method_header(pri_hdr)
        pri_hdr = self._add_reddening_header(pri_hdr)
        red_hdr = fits.Header()
        red_hdr = self._add_reddening_header(red_hdr)
        map_hdr = DAPFitsUtil.build_map_header(self.cube.fluxhdr,
                                               'K Westfall <westfall@ucolick.org>')

        # Get the spatial map mask
        # Marginalize the DRP spectra over wavelength
        # TODO: This needs to be abstracted for a general datacube, and
        # allow for that datacube to have no defined bitmask.
#        if(self.cube.bitmask==None):
#            map_mask = numpy.zeros(self.spatial_shape, dtype=int)
#            ### KHRR -- need to fix this!!!
#
#        else:
#            map_mask = DAPFitsUtil.marginalize_mask(self.cube.mask, ['NOCOV', 'LOWCOV', 'DEADFIBER',
#                                                                     'FORESTAR', 'DONOTUSE'],
#                                                    self.cube.bitmask, self.bitmask,
#                                                    out_flag='DIDNOTUSE')
#            map_mask = DAPFitsUtil.marginalize_mask(self.cube.mask, ['FORESTAR'], self.cube.bitmask,
#                                                    self.bitmask, out_mask=map_mask)
#        drp_bad = (map_mask > 0)
        # Marginalize the spectral masks over wavelength
        map_mask = DAPFitsUtil.marginalize_mask(self.cube.mask,
                                                inp_flags=self.cube.do_not_use_flags(),
                                                inp_bitmask=self.cube.bitmask, 
                                                out_flag='DIDNOTUSE', out_bitmask=self.bitmask)
        if self.cube.propagate_flags() is not None:
            for flag in self.cube.propagate_flags():
                map_mask = DAPFitsUtil.marginalize_mask(self.cube.mask, inp_flags=flag,
                                                        inp_bitmask=self.cube.bitmask,
                                                        out_flag=flag, out_bitmask=self.bitmask,
                                                        out_mask=map_mask)
        drp_bad = map_mask > 0
        # Add the spectra with low spectral coverage

        ### KHRR added this if
#        if(self.cube.bitmask):
#            indx = numpy.invert(good_fgoodpix.reshape(self.spatial_shape)) & numpy.invert(drp_bad)
#            map_mask[indx] = self.bitmask.turn_on(map_mask[indx], 'LOW_SPECCOV')
#            # Add the spectra with low S/N
#            indx = numpy.invert(good_snr.reshape(self.spatial_shape)) & numpy.invert(drp_bad)
#            map_mask[indx] = self.bitmask.turn_on(map_mask[indx], 'LOW_SNR')

        # TODO: (KBW) This should be okay because LOW_SPECCOV and LOW_SNR are
        # DAP specific flags, not DRP ones.  But need to check!

        # Turn on the flag stating that the number of valid channels in
        # the spectrum was below the input fraction.
        indx = numpy.invert(good_fgoodpix.reshape(self.spatial_shape)) & numpy.invert(drp_bad)
        map_mask[indx] = self.bitmask.turn_on(map_mask[indx], 'LOW_SPECCOV')
        # Add the spectra with low S/N
        indx = numpy.invert(good_snr.reshape(self.spatial_shape)) & numpy.invert(drp_bad)
        map_mask[indx] = self.bitmask.turn_on(map_mask[indx], 'LOW_SNR')

        # Fill the covariance HDUs
        if self.covariance is not None:
            pri_hdr, ivar_hdu, covar_hdu = self.covariance.output_hdus(hdr=pri_hdr)

        # Save the data to the hdu attribute
        self.hdu = fits.HDUList([fits.PrimaryHDU(header=pri_hdr),
                                 fits.ImageHDU(data=stack_flux.data, name='FLUX'),
                                 fits.ImageHDU(data=stack_ivar.data, name='IVAR'),
                                 fits.ImageHDU(data=stack_mask, name='MASK'),
                                 fits.ImageHDU(data=self.cube.wave, name='WAVE'),
                                 fits.ImageHDU(data=stack_sres.data, name='SPECRES'),
                                 fits.ImageHDU(data=self.galext.redcorr.data, header=red_hdr,
                                               name='REDCORR'),
                                 fits.ImageHDU(data=stack_sdev.data, name='FLUXD'),
                                 fits.ImageHDU(data=stack_npix.data, name='NPIX'),
                                 fits.ImageHDU(data=bin_indx, header=map_hdr, name='BINID'),
                                 fits.ImageHDU(data=map_mask, header=map_hdr, name='MAPMASK'),
                                 bin_data.to_hdu(name='BINS')])

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

    @staticmethod
    def do_not_fit_flags():
        return ['DIDNOTUSE', 'FORESTAR', 'LOW_SPECCOV', 'LOW_SNR', 'NONE_IN_STACK', 'IVARINVALID']

    def check_fgoodpix(self, minimum_fraction=None):
        r"""
        Determine which spaxels in :attr:`rdxqa` have a fractional
        spectral coverage of greater than the provided minimum fraction.
        Only these spectra will be included in the binning.
    
        Args:
            minimum_fraction (:obj:`float`, optional):
                The minimum fraction of the spectrum that must be
                valid for the spectrum to be included in any bin. If
                None, ``self.method['minimum_frac']`` is used.

        Returns:
            `numpy.ndarray`_: Boolean array for the spectra that
            satisfy the criterion. Shape is :math:`(N_{\rm
            spaxel},)`.
        """
        if minimum_fraction is None:
            minimum_fraction = self.method['minimum_frac']
        return self.rdxqa['SPECTRUM'].data['FGOODPIX'] > minimum_fraction

    def above_snr_limit(self, sn_limit, debug=False):
        """
        Flag bins above a provided S/N limit.

        Args:
            sn_limit (:obj:`float`):
                S/N threshold.
            debug (:obj:`bool`, optional):
                Run in debug mode. This just selects the first 2
                spectra that meet the S/N criterion so that the rest
                of the analysis is just done using those two bins.

        Returns:
            `numpy.ndarray`_: Boolean array selecting those bins that
            meet the S/N threshold.
        """
        if debug:
            warnings.warn('You are setting all but two spectra as bad!')
            test = self.hdu['BINS'].data['SNR'] > sn_limit
            indx = numpy.arange(len(test))[test]
            test[indx[2:]] = False
            return test
        return self.hdu['BINS'].data['SNR'] > sn_limit

    # TODO: Need to abstract further for non-DRP cubes, and for cubes
    # without a mask.
    def get_mask_array_from_cube(self, select=None):
        r"""
        Convert the datacube mask for individual spaxels to the binned
        spectra mask.

        Any pixel with bit flags in the list returned by
        :func:`mangadap.datacube.datacube.DataCube.do_not_stack_flags`
        for :attr:`cube` are consolidated into the ``DIDNOTUSE``
        flag. Any ``FORESTAR`` flags are propagated.

        Args:
            select (`numpy.ndarray`_, optional):
                Flattened, boolean array with the spaxels in the
                datacube to select. Shape is :math:`(N_{\rm
                spaxel},)`. If None, all spaxels are included.

        Returns:
            `numpy.ndarray`_: Integer array with the consolidated
            mask bits.
        """
        drp_mask = self.cube.copy_to_array(attr='mask')
        if select is not None:
            drp_mask = drp_mask[select,:]

        # Initialize
        mask = numpy.zeros(drp_mask.shape, dtype=self.bitmask.minimum_dtype())

        # Consolidate pixels flagged to be excluded from any
        # stacking into DIDNOTUSE
        flags = self.cube.do_not_stack_flags()
        indx = self.cube.bitmask.flagged(drp_mask, flag=flags)
        mask[indx] = self.bitmask.turn_on(mask[indx], 'DIDNOTUSE')

        # Propagate the FORESTAR flags
        indx = self.cube.bitmask.flagged(drp_mask, flag='FORESTAR')
        mask[indx] = self.bitmask.turn_on(mask[indx], 'FORESTAR')

        return mask

    def bin_spectra(self, cube, rdxqa, reff=None, dapver=None, analysis_path=None,
                    directory_path=None, output_file=None, hardcopy=True, symlink_dir=None,
                    clobber=False, loggers=None, quiet=False):
        """
        Bin and stack the spectra.

        This is the core funtion of this class, constructing its main
        data container, :attr:`hdu`.

        .. todo::
            Describe algorithm.

        Args:
            cube (:class:`mangadap.datacube.datacube.DataCube`):
                The datacube with the spectra to bin.
            rdxqa (:class:`mangadap.proc.reductionassessments.ReductionAssessments`):
                The basic assessments of the datacube that are used
                for the binning procedures.
            reff (:obj:`float`, optional):
                The effective radius of the galaxy in arcsec.
            dapver (:obj:`str`, optional):
                The DAP version. Used to construct the output paths,
                overriding the default defined by
                :func:`mangadap.config.defaults.dap_version`. Does
                **not** select the version of the code to use.
            analysis_path (:obj:`str`, optional):
                The top-level path for the DAP output files, used to
                override the default defined by
                :func:`mangadap.config.defaults.dap_analysis_path`.
            directory_path (:obj:`str`, optional):
                The exact path to the directory with DAP output that
                is common to the DAP "methods". Default is defined by
                :func:`mangadap.confgi.defaults.dap_common_path`.
            output_file (:obj:`str`, optional):
                Exact name for the output file. The default is to use
                :func:`mangadap.config.defaults.dap_file_name`.
            hardcopy (:obj:`bool`, optional):
                Flag to write the `astropy.io.fits.HDUList`_
                (:attr:`hdu`) to disk. If False, the object data is
                only kept in memory.
            symlink_dir (:obj:`str`, optional):
                Create a symbolic link to the created file in the
                supplied directory. If None, no symbolic link is
                created.
            clobber (:obj:`bool`, optional):
                Overwrite any existing files. If False, any existing
                files will be used. If True, the analysis is redone
                and any existing output is overwritten.
            loggers (:obj:`list`, optional):
                List of `logging.Logger`_ objects to log progress;
                ignored if quiet=True. Logging is done using
                :func:`mangadap.util.log.log_output`. If None, no
                logging is performed and output is just written to
                ``stdout``.
            quiet (:obj:`bool`, optional):
                Suppress all terminal and logging output.

        Raises:
            ValueError:
                Raised if ``cube`` or ``rdxqa`` is None, or if no
                spectra in the datacube meet the S/N or spectral
                coverage criteria. If the spectra are actually being
                binned (i.e., the binning type is not 'none'), this
                error is also raised if the output file cannot be
                defined, if no spectra are assigned to any bin, or if
                the stacking function results in a correlation matrix
                instead of a full covariance matrix.
            TypeError:
                Raised if the input ``cube`` is not derived from
                :class:`mangadap.datacube.datacube.DataCube` or if
                ``rdxqa`` is not a
                :class:`mangadap.proc.reductionassessments.ReductionAssessment`.
        """
        # Initialize the reporting
        if loggers is not None:
            self.loggers = loggers
        self.quiet = quiet

        # DataCube always needed
        if cube is None:
            raise ValueError('Must provide datacube object!')
        if not isinstance(cube, DataCube):
            raise TypeError('Must provide a valid DataCube object!')

        # Test if the RSS file exists; cannot compute covariance if not
        # TODO: this will break for a different stacking class
        if self.method['stackpar']['covar_mode'] is not 'none' and not cube.can_compute_covariance:
            warnings.warn('Cannot determine covariance matrix!  Continuing without!')
            self.method['stackpar']['covar_mode'] = 'none'

        # TODO: How many of these attributes should I keep vs. just use
        # the cube attributes?
        self.cube = cube
        self.shape = self.cube.shape
        self.spatial_shape = self.cube.spatial_shape
        self.nspec = self.cube.nspec
        self.spatial_index = self.cube.spatial_index.copy()
        self.nwave = self.cube.nwave

        # ReductionAssessment object always needed
        if rdxqa is None:
            raise ValueError('Must provide a ReductionAssessment object to bin spectra.')
        if not isinstance(rdxqa, ReductionAssessment):
            raise TypeError('Must provide a valid ReductionAssessment object!')
        self.rdxqa = rdxqa

        # Set the Galactic extinction correction defined by the method
        # TODO: Abstract this to allow ebvgal to be input directly or
        # be part of the datacube metadata
        self.galext = GalacticExtinction(form=self.method['galactic_reddening'],
                                         wave=self.cube.wave, ebv=self.cube.prihdr['EBVGAL'],
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

            log_output(self.loggers, 1, logging.INFO,
                       'Total spectra: {0}'.format(len(good_fgoodpix)))
            log_output(self.loggers, 1, logging.INFO,
                       'With 80% spectral coverage: {0}'.format(numpy.sum(good_fgoodpix)))
            log_output(self.loggers, 1, logging.INFO,
                       'With good S/N: {0}'.format(numpy.sum(good_snr)))
            log_output(self.loggers, 1, logging.INFO,
                       'Number of good spectra: {0}'.format(numpy.sum(good_spec)))
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
            bin_indx = numpy.full(self.cube.spatial_shape, -1, dtype=numpy.int)
            i = numpy.asarray(tuple(self.cube.spatial_index[good_spec]))
            bin_indx[i[:,0],i[:,1]] = numpy.arange(numpy.sum(good_spec))

            # Build the data arrays directly from the DRP file
            flags = self.cube.do_not_stack_flags()
            flux = self.cube.copy_to_masked_array(flag=flags)[good_spec,:]
            sdev = numpy.ma.zeros(flux.shape, dtype=float)
            ivar = self.cube.copy_to_masked_array(attr='ivar', flag=flags)[good_spec,:]
            sres = self.cube.copy_to_array(attr='sres')[good_spec,:]
            npix = numpy.ones(flux.shape, dtype=numpy.int)
            npix[numpy.ma.getmaskarray(flux)] = 0

            # Build the mask, converting from the DRP bits to the
            # BinnedSpectra bits
            mask = self.get_mask_array_from_cube(select=good_spec)

            # Mask anything with an invalid ivar
            # TODO: Not sure this is necessary
            indx = numpy.logical_not(ivar > 0)
            if numpy.any(indx):
                flux.mask |= indx
                ivar.mask |= indx
                mask[indx] = self.bitmask.turn_on(mask[indx], ['IVARINVALID', 'DIDNOTUSE'])

            # TODO: Does this work with any covariance mode?  Don't like
            # this back and forth between what is supposed to be a stack
            # only function
            # (SpectralStack.build_covariance_data_DRPFits) and
            # something that is supposed to be independent of the details
            # of the stacking implementation (SpatiallyBinnedSpectra)...

            # Try to add the covariance data, if requested
            covar=None
            if self.method['stackpar']['covar_mode'] not in ['none', None]:
                if not self.quiet:
                    log_output(self.loggers, 1, logging.INFO, 'Attempting to compute covariance.')

                try:
                    covar = self.method['stackclass'].build_covariance_data(self.cube,
                                                        self.method['stackpar']['covar_mode'], 
                                                        self.method['stackpar']['covar_par'])
                except AttributeError as e:
                    if not self.quiet:
                        warnings.warn('Could not build covariance data:: '
                                      'AttributeError: {0}'.format(str(e)))
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
            log_output(self.loggers, 1, logging.INFO,
                       'Output path: {0}'.format(self.directory_path))
            log_output(self.loggers, 1, logging.INFO,
                       'Output file: {0}'.format(self.output_file))

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
                self.method['stackfunc'](self.cube, bin_indx, par=self.method['stackpar'])
        if stack_covar is not None and stack_covar.is_correlation:
            raise ValueError('Unexpected to have stack covariance be a correlation matrix.')

        # Construct the mask
        stack_mask = numpy.zeros(stack_flux.shape, dtype=self.bitmask.minimum_dtype())
        indx = numpy.invert(stack_npix>0)
        stack_mask[indx] = self.bitmask.turn_on(stack_mask[indx], ['NONE_IN_STACK', 'NO_STDDEV'])
        indx = numpy.invert(stack_npix>1)
        stack_mask[indx] = self.bitmask.turn_on(stack_mask[indx], 'NO_STDDEV')
        indx = numpy.invert(stack_ivar>0)
        stack_mask[indx] = self.bitmask.turn_on(stack_mask[indx], ['IVARINVALID', 'DIDNOTUSE'])

        #---------------------------------------------------------------
        # Fill the table with the per-bin data (BEFORE applying the
        # reddening)
        bin_data = self._binned_data_table(bin_indx, stack_flux, stack_ivar)

        #---------------------------------------------------------------
        # Deredden the spectra if the method requests it
        stack_flux, stack_ivar, stack_sdev, stack_covar \
                = self._apply_reddening(stack_flux, stack_ivar, stack_sdev, stack_covar)

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
        the datacube from which it was derived.
        """
        # Report
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, 'Constructing binned spectra datacube ...')

        # Get/Copy the necessary data arrays
        bin_indx = self.hdu['BINID'].data.copy()
        stack_flux = self.hdu['FLUX'].data.copy()
        stack_sdev = self.hdu['FLUXD'].data.copy()
        stack_ivar = self.hdu['IVAR'].data.copy()
        stack_mask = self.hdu['MASK'].data.copy()
        stack_npix = self.hdu['NPIX'].data.copy()
        stack_sres = self.hdu['SPECRES'].data.copy()

        # Reconstruct the stacked spectra into a datacube
        flux, mask, sdev, ivar, npix, sres \
                = DAPFitsUtil.reconstruct_cube(self.shape, bin_indx.ravel(),
                                               [stack_flux, stack_mask, stack_sdev, stack_ivar,
                                                stack_npix, stack_sres ])

        # TODO: Move this to the Covariance class?
        if self.covariance is not None:
            covariance = self.covariance.bin_to_spaxel_covariance(bin_indx.ravel())

        # Finalize the mask
        mask = self._finalize_cube_mask(mask)

        # Primary header is identical regardless of the shape of the
        # extensions
        pri_hdr = self.hdu['PRIMARY'].header.copy()
        cube_hdr = DAPFitsUtil.build_cube_header(self.cube, 'K Westfall <westfall@ucolick.org>')

        # Fill the covariance HDUs
        if self.covariance is not None:
            pri_hdr, ivar_hdu, covar_hdu = covariance.output_hdus(reshape=True, hdr=pri_hdr)

        # Save the data to the hdu attribute
        # TODO: Strip headers of sub extensions of everything but the
        # WCS and HDUCLAS information. Keep WCS information in BINID.
        hdu = fits.HDUList([fits.PrimaryHDU(header=pri_hdr),
                            fits.ImageHDU(data=flux, header=cube_hdr, name='FLUX'),
                            fits.ImageHDU(data=ivar, header=cube_hdr, name='IVAR'),
                            fits.ImageHDU(data=mask, header=cube_hdr, name='MASK'),
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

    def write(self, match_datacube=False, clobber=False):
        """
        Write the hdu object to the file.

        Args:
            match_datacube (:obj:`bool`, optional):
                Match the shape of the data arrays to the input
                datacube. I.e., convert them to 3D and replicate the
                binned spectrum to each spaxel in the bin.
            clobber (:obj:`bool`, optional):
                Overwrite any existing file.
        """
        # Convert the spectral arrays in the HDU to a 3D cube and write
        # it
        if match_datacube:
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
            ifile (:obj:`str`, optional):
                Name of the file with the data. Default is to use the
                name provided by :func:`file_path`.
            strict (:obj:`bool`, optional):
                Force a strict reading of the file to make sure that
                it adheres to the expected format. At the moment,
                this only checks to make sure the method keyword
                printed in the header matches the expected value in
                :attr:`method`.
            checksum (:obj:`bool`, optional):
                Use the checksum in the header to check for
                corruption of the data.

        Raises:
            FileNotFoundError:
                Raised if the file does not exist.
            ValueError:
                Raised if ``strict`` is true and the header keyword
                ``BINKEY`` does not match the method keyword.
        """
        if ifile is None:
            ifile = self.file_path()
        if not os.path.isfile(ifile):
            raise FileNotFoundError('File does not exist!: {0}'.format(ifile))

        if self.hdu is not None:
            self.hdu.close()

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

        # Attempt to read and construct the covariance object
        # TODO: Covariance matrices need to have the input_index
        # written to the file
        try:
            self.covariance = Covariance.from_fits(self.hdu, ivar_ext=None, correlation=True)
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

    def copy_to_array(self, ext='FLUX', waverange=None, include_missing=False):
        r"""
        Wrapper for :func:`mangadap.util.fitsutil.DAPFitsUtil.copy_to_array`
        specific for :class:`SpatiallyBinnedSpectra`.

        Return a copy of the selected data array. The array size is
        always :math:`N_{\rm bins} \times N_{\rm wavelength}`; i.e.,
        the data is always flattened to two dimensions and the unique
        spectra are selected.

        Args:
            ext (:obj:`str`, optional):
                Name of the extension from which to draw the data.
                Must be allowed for the current :attr:`mode`; see
                :attr:`data_arrays`.
            waverange (array-like, optional):
                Two-element array with the first and last wavelength
                to include in the computation. If None, use the full
                wavelength range.
            include_missing (:obj:`bool`, optional):
                Create an array with a size that accommodates the
                missing models.

        Returns:
            `numpy.ndarray`_: A 2D array with a copy of the data from
            the selected extension.
        """
        return DAPFitsUtil.copy_to_array(self.hdu, ext=ext, allowed_ext=self.spectral_arrays,
                                waverange=waverange, missing_bins=self.missing_bins
                                                                    if include_missing else None,
                                nbins=self.nbins,
                                unique_bins=DAPFitsUtil.unique_bins(self.hdu['BINID'].data))

    def copy_to_masked_array(self, ext='FLUX', flag=None, waverange=None, include_missing=False):
        """
        Wrapper for
        :func:`mangadap.util.fitsutil.DAPFitsUtil.copy_to_masked_array`
        specific for :class:`SpatiallyBinnedSpectra`.

        Return a copy of the selected data array as a masked array.
        This is functionally identical to :func:`copy_to_array`,
        except the output format is a `numpy.ma.MaskedArray`_. The
        pixels that are considered to be masked can be specified
        using `flag`.
        
        Args:
            ext (:obj:`str`, optional):
                Name of the extension from which to draw the data.
                Must be allowed for the current :attr:`mode`; see
                :attr:`data_arrays`.
            flag (:obj:`str`, :obj:`list`, optional):
                (List of) Flag names that are considered when
                deciding if a pixel should be masked. The names
                *must* be a valid bit name as defined by
                :attr:`bitmask`. If not provided, *ANY* non-zero mask
                bit is omitted.
            waverange (array-like, optional):
                Two-element array with the first and last wavelength
                to include in the computation. If None, use the full
                wavelength range.
            include_missing (:obj:`bool`, optional):
                Create an array with a size that accommodates the
                missing models.

        Returns:
            `numpy.ma.MaskedArray`_: A 2D array with a copy of the
            data from the selected extension.
        """
        return DAPFitsUtil.copy_to_masked_array(self.hdu, ext=ext, mask_ext='MASK', flag=flag,
                                bitmask=self.bitmask, allowed_ext=self.spectral_arrays,
                                waverange=waverange, missing_bins=self.missing_bins
                                                                    if include_missing else None,
                                nbins=self.nbins,
                                unique_bins=DAPFitsUtil.unique_bins(self.hdu['BINID'].data))

    def get_bin_indices(self, bins):
        """
        Return the indices of the bins in the BIN table.

        Args:
            bins (array-like):
                The bin ID numbers to find
        
        Returns:
            `numpy.ndarray`_: Integer array with the index of each
            bin ID in the BINID columns of the BINS extension.
        """
        return numpy.array([numpy.where(self.hdu['BINS'].data['BINID'] == b)[0][0] for b in bins])

    def find_nearest_bin(self, input_bins, weighted=False, indices=False):
        """
        Use a `scipy.spatial.KDTree`_ to find the bins nearest to and
        excluding a list of input bins.

        Args:
            input_bins (array-like):
                One or more bin ID numbers to use to locate the bin
                nearest to it based on its on-sky coordinates. The
                list must be a unique set.
            weighted (:obj:`bool`, optional):
                Use the weighted coordinates (LW_SKY_COO) instead of
                the unweighted coordinates (SKY_COO) when finding the
                nearest bin.
            indices (:obj:`bool`, optional):
                Return the indices of the nearest bins instead of
                their ID number (default).

        Returns:
            `numpy.ndarray`_: The bin IDs, one per input bin, of the
            bin closest to each input bin. Any bins that are in
            :attr:`missing_bins` have a return value of -1; there are
            no coordinates for these bins.

        Raises:
            ValueError:
                Raised if the set of input bins is not unique or
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
        ref_coo = numpy.array([self.hdu['BINS'].data[col][reference_bin,0],
                               self.hdu['BINS'].data[col][reference_bin,1] ]).T

        # Construct the KDTree
        kd = spatial.KDTree(ref_coo)

        # Get the indices of the input bins in the internal list
        input_bin_indx = self.get_bin_indices(_input_bins)

        # Get the coordinates of the bins to match to the nearest
        # reference bin
        input_coo = numpy.array([self.hdu['BINS'].data[col][input_bin_indx,0],
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
            data (array-like):
                Data for each bin. The length must be the same as
                :attr:`nbins`.
            bad_bins (array-like):
                The list of indices (must not be a boolean array)
                with bad values to be replaced.
        
        Returns:
            `numpy.ndarray`_: A new array with the bad data filled
            with the data from the nearest bin.

        Raises:
            ValueError:
                Raised if the input array doesn't have the correct
                shape or if the list of bad bins has numbers outside
                the viable range (0,self.nbins-1).
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

