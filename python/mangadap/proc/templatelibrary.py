# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
r"""
Class that reads and prepares template libraries for use in fitting the
stellar-continuum of a spectrum.  See
:func:`available_template_libraries` for the list of default template
libraries.

The "raw" template libraries are expected to consist of 1D fits files.
The wavelength sampling can be either linear or logarithmic in
wavelength.  The reference frame of the template wavelengths must be
defined as either vacuum or air.  It is expected that the object spectra
to be fit are calibrated to vacuum wavelengths.  When preparing the
template spectra for analysis, this class will use
`pydl.goddard.astro.airtovac`_ to convert the template wavelengths to
vacuum.  Finally, one can specify that the template library is only
valid within a certain wavelength range and above a certian flux limit;
see :class:`TemplateLibraryDef`.

Preparation of the template library for use in stellar-continuum fitting
consists of a number of steps as described in the documentation of
:func:`TemplateLibrary.process_template_library`.

A template library that has been prepared for analysis is automatically
written to disk for later recovery.

Two support classes are also provided. One is a derived
:class:`mangadap.par.parset.KeywordParSet` instance that provides the
defining parameters of a DAP template library. The second is a
derived :class:`mangadap.util.bitmask.BitMask` instance that defines
the bitmasks for the template library spectra.

.. _templatelibrary-usage:

Class usage examples
--------------------

Assuming you have the default directory structure setup, you can do::

    # Imports
    from mangadap.drpfits import DRPFits
    from mangadap.proc.templatelibrary import TemplateLibrary
    from matplotlib import pyplot

    # Define the DRP file
    drpf = DRPFits(7495, 12703, 'CUBE')

    # Build the template library
    tpl_lib = TemplateLibrary('M11MILES', drpf=drpf, directory_path='.')
    # Writes: ./manga-7495-12703-LOGCUBE_M11-MILES.fits

    # Plot one of the spectra
    pyplot.plot(tpl_lib.hdu['WAVE'].data, tpl_lib.hdu['FLUX'].data[0,:])
    pyplot.show()

As part of the instantiation of the :class:`TemplateLibrary` object
in the above call, the template library is prepared for use in
fitting the specified DRP file.  If the required processing has
already been done, the instantiation of the :class:`TemplateLibrary`
object simply reads the existing file.  If you do not have the
default directory structure setup, you'll need to define the paths
to, e.g., the DRP file; see :class:`mangadap.drpfits.DRPFits`.

If you do not want to process the template (match the spectral
resolution and sampling to that of the DRP data), you can force the
:class:`TemplateLibrary` object to only provide the "raw" spectra::

    # Imports
    from mangadap.proc.templatelibrary import TemplateLibrary
    from matplotlib import pyplot

    # Read the raw template library
    tpl_lib = TemplateLibrary('M11MILES', process=False)
    # Nothing should be written to disk

    # Plot one of the spectra
    pyplot.plot(tpl_lib.hdu['WAVE'].data[0,:], tpl_lib.hdu['FLUX'].data[0,:])
    pyplot.show()

Note that in the previous example, the wavelength data was already
one-dimensional, whereas for the raw library, the wavelength vector
can be spectrum-dependent.

In the above examples, the user has not provided a list of template
libraries, meaning that the default set available to the DAP is
used.  The default set is defined in
:func:`available_template_libraries`.  If you want to use your own
template library, you have to define its parameters using
:class:`TemplateLibraryDef`.  Currently, the template library
spectra are expected to be 1D fits files with WCS header keywords
defining the wavelength solution; see above.  Using an existing DAP
library as an example::

    # Imports
    from mangadap.proc.templatelibrary import TemplateLibraryDef
    from mangadap.proc.templatelibrary import TemplateLibrary
    from mangadap.config.defaults import dap_source_dir
    from matplotlib import pyplot

    # Define the search string for the library
    search_str = dap_source_dir()+'/data/stellar_templates/miles/*.fits'

    # Define the template library parameters
    tpl_par = TemplateLibraryDef(key='MILES',    # Unique keyword for the library
                            file_search=search_str, # Search string
                            fwhm=2.50,              # FWHM of resolution element (assumed const)
                            in_vacuum=False,        # Wavelength in vacuum?
                            wave_limit=numpy.array([ 3575., 7400. ]),   # Range of valid lambda
                            lower_flux_limit=0.0)   # Lower limit for valid flux

    # Read the raw template library
    tpl_lib = TemplateLibrary('MILES', tpllib_list=tpl_par, process=False)
    # Nothing should be written to disk

    # Plot one of the spectra
    pyplot.plot(tpl_lib.hdu['WAVE'].data[0,:], tpl_lib.hdu['FLUX'].data[0,:])
    pyplot.show()

Note that the keyword you use must be provided both to the parameter
set and when instantiating the :class:`TemplateLibrary` object.  In
the example above, I have not processed the library, but you can by
following a similar approach to the first example.

You can also process the spectra to a user-provided resolution and
pixel sampling, toggle on/off the renormalization of the library to
a mean flux of unity and define various paths if you're not working
within the nominal DAP directory structure.  See the optional
instantiation arguments for :class:`TemplateLibrary`.

Revision history
----------------

    | **23 Apr 2015**: Implementation begun by K. Westfall (KBW)
    | **26 May 2015**: (KBW) Added some Sphinx documentation.
    | **17 Jun 2015**: (KBW) Added flexibility in definition of template
        libraries from which to choose using the new
        :class:`TemplateLibraryDef` class.
    | **23 Jun 2015**: (KBW) Allow user to provided non-DRP input
        spectra, meaning they need to provide the velocity scale and the
        wavelength and spectral resolution vectors.  They must also
        directly set the name of the output processed file.
    | **01 Feb 2016**: (KBW) Propagated changes to
        :class:`TemplateLibraryDef` and changed procedure for reading
        raw template library.
    | **04 Feb 2016**: (KBW) Converted from
        :func:`mangadap.util.instrument.log_rebin` to
        :func:`mangadap.util.instrument.resample_vector`.  Allows input
        library to be logarithmically sampled and provides a different
        approach to subsampling spectra.  Implemented other convenience
        operations, such as selecting a certain wavelength range for the
        processed library.
    | **19 May 2016**: (KBW) Added loggers and quiet keyword arguments
        to :class:`TemplateLibrary`
    | **22 Jun 2016**: (KBW) Included MILESHC library in documentation.
        Allow to specify how the resolution and sampling are matched to
        the DRP data.  The :class:`TemplateLibrary` class should be
        generalized to make this more transparent (and unnecessary).
    | **03 Apr 2017**: (KBW) Include arguments for
        :func:`mangadap.util.instrument.match_spectral_resolution` that
        specify the minimum sigma in pixels for the convolution kernel
        and any offset in the resolution required in the call for the
        template library.  Check that the resampling of the spectrum
        does not reach the two-pixel resolution limit; flag the spectrum
        in those regions that do and change the resolution to the two
        pixel limit.
    | **30 Aug 2017**: (KBW) Switch from using
        :func:`mangadap.util.instrument.resample_vector` to
        :func:`mangadap.util.instrument.resample1d`
    | **30 Aug 2018**: (KBW) Switch from using resample1d to
        :class:`mangadap.util.sampling.Resample`.

----

.. include license and copyright
.. include:: ../copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../links.rst
"""

import os
import glob
import time
import warnings
import logging

import numpy
from scipy import sparse, interpolate

from astropy.io import fits
import astropy.constants

from pydl.goddard.astro import airtovac

from ..par.parset import KeywordParSet
from ..config.defaults import dap_source_dir, default_dap_common_path
from ..config.defaults import default_dap_file_name
from ..util.bitmask import BitMask
from ..util.dapbitmask import DAPBitMask
from ..util.log import log_output
from ..util.fileio import readfits_1dspec, read_template_spectrum, writefits_1dspec, create_symlink
from ..util.fitsutil import DAPFitsUtil
from ..util.sampling import Resample, resample_vector_npix, angstroms_per_pixel
from ..util.sampling import spectral_coordinate_step, spectrum_velocity_scale
from ..util.resolution import SpectralResolution, match_spectral_resolution
from ..util.parser import DefaultConfig
from .util import select_proc_method, HDUList_mask_wavelengths

from matplotlib import pyplot

# Add strict versioning
# from distutils.version import StrictVersion

class TemplateLibraryDef(KeywordParSet):
    """
    Class with parameters used to define the template library.
    See :class:`mangadap.par.parset.ParSet` for attributes.

    The defined parameters are:

    .. include:: ../tables/templatelibrarydef.rst
    """
    def __init__(self, key=None, file_search=None, fwhm=None, sres_ext=None, in_vacuum=False,
                 wave_limit=None, lower_flux_limit=None, log10=False): 
        # Perform some checks of the input
        in_fl = [ int, float ]
        
        pars =   [ 'key', 'file_search', 'fwhm', 'sres_ext', 'in_vacuum',  'wave_limit',
                        'lower_flux_limit', 'log10' ]
        values = [   key,   file_search,   fwhm,   sres_ext,   in_vacuum,    wave_limit,
                          lower_flux_limit,   log10 ]
        dtypes = [   str,           str,  in_fl,        str,        bool, numpy.ndarray,
                                     in_fl,    bool ]
        descr = ['Keyword to distinguish the template library.',
                 'Search string used by glob to find the 1D fits spectra to include in the ' \
                    'template library.',
                 'FWHM of the resolution element in angstroms.',
                 'Extension in the fits files with measurements of the spectral resolution as ' \
                    'a function of wavelength.',
                 'Flag that the wavelengths of the spectra are in vacuum, not air.',
                 'Two-element array with the starting and ending wavelengths for the valid ' \
                    'spectral range of the templates.',
                 'Minimum valid flux in the template spectra.',
                 'Flag that the template spectra have been binned logarithmically in wavelength.']

        super(TemplateLibraryDef, self).__init__(pars, values=values, dtypes=dtypes, descr=descr)


def validate_spectral_template_config(cnfg):
    """ 
    Validate the :class:`mangadap.util.parser.DefaultConfig` object with
    the template library parameters.

    Args:
        cnfg (:class:`mangadap.util.parser.DefaultConfig`): Object with
            the template library parameters to validate.

    Raises:
        KeyError: Raised if required keyword does not exist.
        ValueError: Raised if key has unacceptable value.

    """
    # Check for required keywords
    required_keywords = [ 'key', 'file_search' ]
    if not cnfg.all_required(required_keywords):
        raise KeyError('Keywords {0} must all have valid values.'.format(required_keywords))
    
    if not cnfg.keyword_specified('fwhm') and not cnfg.keyword_specified('sres_ext'):
        raise KeyError('Must provide either \'fwhm\' or \'sres_ext\'.')


def available_template_libraries(dapsrc=None):
    r"""
    Return the list of library keys, the searchable string of the 1D
    template library fits files for the template libraries available for
    use by the DAP, the FWHM of the libraries, and whether or not the
    wavelengths are in vacuum.

    The stellar template library files should be a list of 1D fits
    files, and be associated with one of the following library keys:

    Args:
        dapsrc (:obj:`str`, optional):
            Root path to the DAP source directory.  If not provided, the
            default is defined by
            :func:`mangadap.config.defaults.dap_source_dir`.

    Returns:
        list: A list of
        :func:`mangadap.proc.templatelibrary.TemplateLibraryDef`
        objects, each defining a separate template library.

    Raises:
        NotADirectoryError:
            Raised if the provided or default *dapsrc* is not a
            directory.
        OSError/IOError:
            Raised if no template configuration files could be found.
        KeyError:
            Raised if the template-library keywords are not all unique.

    .. todo::

        - Add backup function for Python 2.
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
    ini_files = glob.glob(dapsrc+'/python/mangadap/config/spectral_templates/*.ini')
    if len(ini_files) == 0:
        raise IOError('Could not find any configuration files in {0} !'.format(
                      dapsrc+'/python/mangadap/config/spectral_templates'))

    # Build the list of library definitions
    template_libraries = []
    for f in ini_files:
        # Read and validate the config file
        cnfg = DefaultConfig(f=f, interpolate=True)
        validate_spectral_template_config(cnfg)

        # Convert wave_limit and lower_flux_limit to types acceptable by
        # TemplateLibraryDef
        wave_limit = None if cnfg['wave_limit'] is None \
                        else numpy.array(cnfg.getlist('wave_limit', evaluate=True))
        if wave_limit is not None and len(wave_limit) != 2:
            raise ValueError('Must specify two wavelength limits, can be \'None, None\'.')

        # Append the definition of the template library 
        template_libraries += \
            [ TemplateLibraryDef(key=cnfg['key'], file_search=cnfg['file_search'],
                                 fwhm=cnfg.getfloat('fwhm'), sres_ext=cnfg['sres_ext'],
                                 in_vacuum=cnfg.getbool('in_vacuum', default=False),
                                 wave_limit=wave_limit,
                                 lower_flux_limit=cnfg.getfloat('lower_flux_limit'),
                                 log10=cnfg.getbool('log10', default=False) )
            ]

    # Check the keywords of the libraries are all unique
    if len(numpy.unique( numpy.array([tpl['key'] for tpl in template_libraries]) )) \
            != len(template_libraries):
        raise KeyError('Template-library keywords are not all unique!')

    # Return the default list of template libraries
    return template_libraries


class TemplateLibraryBitMask(DAPBitMask):
    """
    Derived class that specifies the mask bits for the template library
    data.  The maskbits defined are:

    .. include:: ../tables/templatelibrarybitmask.rst
    """
    cfg_root = 'spectral_template_bits'


class TemplateLibrary:
    r"""
    Object used to read, store, and prepare template libraries for use
    in analyzing object spectra.

    The default list of available libraries provided by the MaNGA DAP
    defined in :func:`available_template_libraries`.  The user can
    provide their own library for use with this class provided they are
    contained in 1D fits spectra, sampled linearly in wavelength with
    the wavelength coordinates available via the WCS keywords
    (``CRPIX1``, ``CRVAL1``, ``CDELT1``), and they have an appropriately
    defined spectral resolution (FWHM in angstroms that is constant as a
    function of wavelength).  See :class:`TemplateLibraryDef` and
    :func:`_build_raw_hdu`.

    The class is optimized for use in analyzing MaNGA DRP files;
    however, one can provide the necessary information so that the class
    can be used with a non-DRP spectrum.  In the latter case, the user
    must supply the velocity scale of the pixel for the logarithmically
    resampled template library, and a
    :class:`mangadap.util.resolution.SpectralResolution` object the
    defines the instrumental resolution of the spectrum/spectra to be
    analyzed.

    .. todo::

        - below is out of date.
        - Only works with DRP files that have log-linear wavelength
          binning!
        - Allow to process, where process is just to change the sampling
          or the resolution (not necessarily both).
        - Need to make this more general, removing all dependence on
          DRPFits object.  This would simplify the functionality to
          change how the resolution and sampling matching is specified.

    On initialization, if the DRP file object is not provided (is None),
    the default behavior is to read the raw template library if
    read=True.  If the DRP file is provided, the routine will check for
    the resolution matched fits file; if it doesn't exist and read is
    ``True``, it will prepare the template library for use in analyzing
    the DRP file and write the prepared library file.  If clobber=True,
    the preparation and writing of the template library will be done
    even if the library already exists.

    Args:
        library_key (:obj:`str`):
            Keyword selecting the library to use.
        tpllib_list (:obj:`list`, optional):
            List of :class:`TemplateLibraryDef` objects that define the
            parameters required to read and interpret a template
            library.  The ``library_key`` must select one of the objects
            in this list.
        dapsrc (:obj:`str`, optional):
            Root path to the DAP source directory.  If not provided, the
            default is defined by
            :func:`mangadap.config.defaults.dap_source_dir`.
        drpf (:class:`mangadap.drpfits.DRPFits`, optional):
            DRP file (object) with which the template library is
            associated for analysis
        match_to_drp_resolution (:obj:`bool`, optional):
            Match the spectral resolution of the template library to the
            resolution provided by the :class:`mangadap.drpfits.DRPFits`
            object; the latter must be provided for this argument to
            have any use.
        velscale_ratio (:obj:`int`, optional):
            Resample the template spectra such that the ratio of the
            pixel scale in the provided
            :class:`mangadap.drpfits.DRPFits` object is this many times
            larger than the pixels in the resampled template spectrum.
        sres (:class:`mangadap.util.resolution.SpectralResolution`, optional):
            The object is used simply to access the spectral resolution
            and associated wavelength coordinate vector needed when
            matching the spectral resolution of the template library;
            this is used in place of the attributes in any provided DRP
            file object.
        velocity_offset (:obj:`float`, optional):
            Velocity offset to use when matching the spectral resolution
            between the template library and the galaxy spectra.
        spectral_step (:obj:`float`, optional):
            Target logarithmic (``log=True``) or linear (``log=False``)
            step in wavelength for the template library.
        log (:obj:`bool`, optional):
            Flag to force the library to be logarithmically sampled in
            wavelength.
        wavelength_range (array-like, optional):
            Force the template library to covert this spectral range.
            Unobserved spectral regions will be flagged.
        renormalize (:obj:`bool`, optional):
            After processing, renormalize the flux to unity.
        dapver (:obj:`str`, optional):
            DAP version, which is used to define the default DAP
            analysis path.  Default is defined by
            :func:`mangadap.config.defaults.default_dap_version`
        analysis_path (:obj:`str`, optional):
            The path to the top level directory containing the DAP
            output files for a given DRP and DAP version.  Default is
            defined by
            :func:`mangadap.config.defaults.default_analysis_path`.
        directory_path (:obj:`str`, optional):
            The exact path to the processed template library file.
            Default is defined by
            :func:`mangadap.config.defaults.default_dap_common_path`.
        processed_file (:obj:`str`, optional):
            The name of the file containing the prepared template
            library output file.  The file should be found at
            :attr:`directory_path`/:attr:`processed_file`.  Default is
            defined by
            :func:`mangadap.config.defaults.default_dap_file_name`.
        process (:obj:`bool`, optional):
            If :attr:`drpf` is defined and the prepared template library
            does not exist, this will process the template library in
            preparation for use in fitting the provided DRP file.
        hardcopy (:obj:`bool`, optional):
            Flag to keep a hardcopy of the processed template library.
            Default is True.
        symlink_dir (:obj:`str`, optional):
            Create a symlink to the file in this directory.  Default is
            for no symlink.
        clobber (:obj:`bool`, optional):
            If :attr:`drpf` is define and ``process=True``, this will
            clobber any existing processed template library.

    Attributes:
        bitmask (class:`mangadap.util.bitmask.BitMask`):
            Object used to toggle mask values; see
            :func:`TemplateLibraryBitMask`.
        library (:class:`TemplateLibraryDef`):
            Parameter set required to read and prepare the library.
        file_list (:obj:`list`):
            The list of files found using `glob.glob`_ and
            :attr:`file_search`.
        ntpl (:obj:`int`):
            Number of template spectra in the library
        drpf (:class:`mangadap.drpfits.DRPFits`):
            DRP file (object) with which the template library is
            associated for analysis
        sres (:class:`mangadap.util.resolution.SpectralResolution`):
            The object is used simply to access the spectral resolution
            and associated wavelength coordinate vector needed when
            matching the spectral resolution of the template library;
            this is used in place of the attributes in any provided DRP
            file object.
        velocity_offset (:obj:`float`):
            Velocity offset to use when matching the spectral resolution
            between the template library and the galaxy spectra.
        spectral_step (:obj:`float`):
            Target logarithmic (``log10_sampling=True``) or linear
            (``log10_sampling=False``) step in wavelength for the
            template library.
        log10_sampling (:obj:`bool`):
            Flag that the processed template library is logarithmically
            sampled in wavelength.
        directory_path (:obj:`str`):
            The exact path to the processed template library file.
        processed_file (:obj:`str`):
            The name of the file containing (to contain) the prepared
            template library output file.  The file should be found at
            :attr:`directory_path`/:attr:`processed_file`.
        processed (:obj:`bool`):
            Flag that the template library has been prepared for use in
            the DAP.
        hardcopy (:obj:`bool`):
            Flag to keep a hardcopy of the processed template library.
        symlink_dir (:obj:`str`):
            Symlink created to the file in this directory
        hdu (`astropy.io.fits.hdu.hdulist.HDUList`_):
            HDUList read from the DAP file

    """
    # Class attribute
    supported_libraries = ['BC03', 'BPASS', 'M11ELODIE', 'M11MARCS', 'M11MILES', 'M11STELIB',   
                           'M11STELIBZSOL', 'MASTARHC', 'MILES', 'MILESAVG', 'MILESHC',
                           'MILESTHIN', 'MIUSCAT', 'MIUSCATTHIN', 'STELIB']

    def __init__(self, library_key, tpllib_list=None, dapsrc=None, drpf=None,
                 match_to_drp_resolution=True, velscale_ratio=None, sres=None, velocity_offset=0.0,
                 min_sig_pix=0.0, no_offset=True, spectral_step=None, log=True,
                 wavelength_range=None, renormalize=True, dapver=None, analysis_path=None,
                 directory_path=None, processed_file=None, read=True, process=True, hardcopy=True,
                 symlink_dir=None, clobber=False, checksum=False, loggers=None, quiet=False):

        self.loggers = loggers
        self.quiet = quiet

        # Define the TemplateLibraryBitMask object
        self.bitmask = TemplateLibraryBitMask()

        # Define the properties needed to modify the spectral resolution
        self.sres = None
        self.velocity_offset = None
        self.min_sig_pix = None
        self.no_offset = None

        # Define the target spectral sampling properties
        self.spectral_step = None
        self.log10_sampling = None

#        self.velscale = None
#        velscale (float): The velocity scale of the spectrum that will
#            be analyzed with the library; this is used in place of the
#            attributes in any provided DRP file object.

        # Define the library
        self.library = None
        self.file_list = None
        self.ntpl = None
        self.library = self.define_library(library_key, tpllib_list=tpllib_list, dapsrc=dapsrc)

        # Define the processed file and flag, and the HDUList used to
        # keep the data
        self.directory_path = None
        self.processed_file = None
        self.processed = False
        self.hardcopy = True
        self.symlink_dir = None
        self.hdu = None
        self.checksum = checksum

        # Do not read the library
        if not read:
            if not self.quiet:
                log_output(self.loggers, 1, logging.INFO, 'Nothing read, by request.')
            return
       
        # Do not process the library
        if not process:
            if not self.quiet:
                log_output(self.loggers, 1, logging.INFO, 'Reading raw library without processing.')
            self._read_raw()
            return

        # Read and process the library
        self.process_template_library(drpf=drpf, match_to_drp_resolution=match_to_drp_resolution,
                                      velscale_ratio=velscale_ratio, sres=sres,
                                      velocity_offset=velocity_offset, min_sig_pix=min_sig_pix,
                                      no_offset=no_offset, spectral_step=spectral_step, log=log,
                                      wavelength_range=wavelength_range, renormalize=renormalize,
                                      dapver=dapver, analysis_path=analysis_path,
                                      directory_path=directory_path, processed_file=processed_file,
                                      hardcopy=hardcopy, symlink_dir=symlink_dir, clobber=clobber,
                                      loggers=loggers, quiet=quiet)

    def __getitem__(self, key):
        return self.hdu[key]

    @staticmethod
    def define_library(library_key, tpllib_list=None, dapsrc=None):
        """
        Select the library from the provided list.  Used to set
        :attr:`library`; see
        :func:`mangadap.proc.util.select_proc_method`.

        Args:
            library_key (:obj:`str`):
                Keyword of the selected library.  Available libraries
                are proved by :func:`available__template_libraries`
            tpllib_list (:obj:`list`, optional):
                List of :class:`TemplateLibraryDef` objects that define
                the parameters required to read and interpret a template
                library.
            dapsrc (:obj:`str`, optional):
                Root path to the DAP source directory.  If not provided,
                the default is defined by
                :func:`mangadap.config.defaults.dap_source_dir`.
        """
        # Get the details of the selected template library
        return select_proc_method(library_key, TemplateLibraryDef, method_list=tpllib_list,
                                  available_func=available_template_libraries)

    def _can_set_paths(self, directory_path, drpf, processed_file, quiet=False):
        # Check that the directory_path can be set
        if drpf is None and directory_path is None:
            if not quiet:
                raise ValueError('Cannot define the default directory path without a DRP file; ' \
                                 'provide a directory path.')
            return False

        # Check that the processed file can be set
        if drpf is None and processed_file is None:
            if not quiet:
                raise ValueError('Cannot define the default output file name without a DRP file; ' \
                                 'provide a file name for the processed template library.')
            return False

        return True


    def _set_paths(self, directory_path, dapver, analysis_path, drpf, processed_file):
        """
        Set the I/O path to the processed template library.  Used to set
        :attr:`directory_path` and :attr:`processed_file`.  If not
        provided, the defaults are set using, respectively,
        :func:`mangadap.config.defaults.default_dap_common_path` and
        :func:`mangadap.config.defaults.default_dap_file_name`.

        Args:
            directory_path (str): The exact path to the DAP template
                library file.  See :attr:`directory_path`.
            dapver (str): DAP version.
            analysis_path (str): The path to the top-level directory
                containing the DAP output files for a given DRP and DAP
                version.
            drpf (:class:`mangadap.drpfits.DRPFits`): The container
                object of the DRP file that is used to construct the
                path for the processed template library.
            processed_file (str): The name of the file with the prepared
                template library.  See :attr:`processed_file`.

        """
        # Use this to raise the necessary exceptions
        self._can_set_paths(directory_path, drpf, processed_file)

        # Set the output directory path
        self.directory_path = default_dap_common_path(plate=drpf.plate, ifudesign=drpf.ifudesign,
                                                      drpver=drpf.drpver, dapver=dapver,
                                                      analysis_path=analysis_path) \
                                        if directory_path is None else str(directory_path)

        # Set the output file
        self.processed_file = default_dap_file_name(drpf.plate, drpf.ifudesign,
                                                    self.library['key']) \
                                        if processed_file is None else str(processed_file)

    
    def _get_file_list(self):
        files = glob.glob(self.library['file_search'])
        if len(files) == 0:
            raise ValueError('Library search string did not find any files!')
        return numpy.sort(files)


    def _read_raw(self):
        """
        Read the 'raw' versions of the template library; i.e., the
        library before it has been resolution and sampling matched to a
        DRP file.

        The file list read by the search key is sorted for consistency.
        """
        self.file_list = self._get_file_list()
        self.ntpl = len(self.file_list)
        npix = self._get_nchannels()
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO,
                            'Found {0} {1} templates'.format(self.ntpl, self.library['key']))
            log_output(self.loggers, 1, logging.INFO,
                            'Maximum number of wavelength channels: {0}'.format(npix))
        self._build_raw_hdu(npix)


    def _get_nchannels(self):
        """
        Get the maximum number of wavelength channels needed to store
        all the template spectra in a single array.

        .. todo::
            - What happens if the spectrum has an empty primary
              extension?

        Returns:
            int: Maximum number of pixels used by the listed spectra.

        Raises:
            ValueError: Raised if the input template spectra are not
                one-dimensional.
        """
        max_npix = 0
        for f in self.file_list:
            if fits.getval(f, 'NAXIS') != 1:
                raise ValueError('{0} is not one dimensional!'.format(f))
            npix = fits.getval(f, 'NAXIS1')
            if max_npix < npix:
                max_npix = npix
        return max_npix


    def _build_raw_hdu(self, npix):
        r"""
        Build the "raw" template library arrays.  This simply reads the
        provided list of fits files and puts them into arrays of size
        :math:`N_{\rm tpl} \times N_{\rm pix}`.

        This will *force* reading of the data, even if the :attr:`hdu`
        is already initialized.

        The :attr:`hdu` will contain the appropriate extensions, but it
        is important to note that the wavelength vectors will **not**
        necessarily be the same.  That is, reading of the raw template
        spectra can accommodate spectra that have different wavelength
        coordinates.  Any pixels that have no data are masked using the
        'NO_DATA' bitmask flag; see
        :func:`TemplateLibraryBitMask`.

        The spectral resolution is set using :attr:`fwhm`, and the
        spectral resolution offset is initialized to zero (see
        :func:`mangadap.util.resolution.GaussianKernelDifference`).

        .. warning::

            Currently no errors are saved because none are expected for
            the template libraries.

        Args:
            npix (:obj:`int`):
                Number of spectral channels for the output arrays

        """
        if self.hdu is not None:
            if not self.quiet:
                log_output(self.loggers, 1, logging.INFO, 'Closing existing HDUList.')
            self.hdu.close()
            self.hdu = None

        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, 'Attempting to build raw data ...')
        # Allocate the vectors
        wave = numpy.zeros((self.ntpl, npix), dtype=numpy.float64)
        flux = numpy.zeros((self.ntpl, npix), dtype=numpy.float64)
        sres = numpy.zeros((self.ntpl, npix), dtype=numpy.float64)
        mask = numpy.zeros((self.ntpl, npix), dtype=self.bitmask.minimum_dtype())
        soff = numpy.zeros(self.ntpl, dtype=numpy.float64)
#        ivar = numpy.zeros((self.ntpl, npix), dtype=numpy.float64)
#        ivar[:] = 1.0

        # Read and save each spectrum and mask the unobserved
        # wavelengths
        for i in range(0,self.ntpl):
            if self.library['sres_ext'] is None:
                wave_, flux_ = read_template_spectrum(self.file_list[i],
                                                      log10=self.library['log10'])
                wave[i,:wave_.size] = numpy.copy(wave_)
                flux[i,:wave_.size] = numpy.copy(flux_)
                # Set the spectral resolution
                sres = wave/self.library['fwhm']
            else:
                wave_, flux_, sres_ = read_template_spectrum(self.file_list[i],
                                                             sres_ext=self.library['sres_ext'],
                                                             log10=self.library['log10'])
                wave[i,:wave_.size] = numpy.copy(wave_)
                flux[i,:wave_.size] = numpy.copy(flux_)
                sres[i,:wave_.size] = numpy.copy(sres_)

            if wave_.size != npix:
                mask[i,wave_.size:] = self.bitmask.turn_on(mask[i,wave_.size:],'NO_DATA')
            if self.library['lower_flux_limit'] is not None:
                indx = numpy.invert( flux_ > self.library['lower_flux_limit'] )
                mask[i,indx] = self.bitmask.turn_on(mask[i,indx], 'FLUX_INVALID')
            if self.library['wave_limit'] is not None and self.library['wave_limit'][0] is not None:
                indx = wave[i,:].ravel() < self.library['wave_limit'][0]
                mask[i,indx] = self.bitmask.turn_on(mask[i,indx], 'WAVE_INVALID')
            if self.library['wave_limit'] is not None and self.library['wave_limit'][1] is not None:
                indx = wave[i,:].ravel() > self.library['wave_limit'][1]
                mask[i,indx] = self.bitmask.turn_on(mask[i,indx], 'WAVE_INVALID')

        # (Re)Set the HDUList object
        self._reset_hdu(wave, flux, mask, sres, soff)

        # Add some keywords to the header
        self.hdu[0].header['TPLPROC'] = (0, 'Flag that library has been processed')
        self.processed = False
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, '... done')

    def _reset_hdu(self, wave, flux, mask, sres, soff):
        r"""
        (Re)Set :attr:`hdu` to a new HDUList object using the input
        arrays.  Also sets the header items indicating the version of
        the template reader and the keyword for the library.

        .. warning:
            No checking is done concerning the size of each array!

        Args:
            wave (numpy.ndarray): Array with the wavelengths of each
                pixel.
            flux (numpy.ndarray): Array with the flux in each pixel.
            mask (numpy.ndarray): Bitmask values for each pixel.
            sres (numpy.ndarray): Spectral resolution,
                :math:`R=\lambda/\delta\lambda`, at each pixel.
            soff (numpy.ndarray): The spectral resolution offset for
                each spectrum (see
                :func:`mangadap.util.resolution.GaussianKernelDifference`). 

        """
        if self.hdu is not None:
            self.hdu.close()
            self.hdu = None

        self.hdu = fits.HDUList([ fits.PrimaryHDU(),
                                  fits.ImageHDU(wave, name='WAVE'),
                                  fits.ImageHDU(flux, name='FLUX'),
                                  fits.ImageHDU(mask, name='MASK'),
                                  fits.ImageHDU(sres, name='SPECRES'),
                                  fits.ImageHDU(soff, name='SIGOFF')
                                ])

#        self.hdu['PRIMARY'].header['VDAPTPL'] = (self.version, 'Version of DAP template reader')
        self.hdu['PRIMARY'].header['LIBKEY'] = (self.library['key'], 'Library identifier')


    def _wavelength_range(self, flag=None):
        """
        Return the valid wavelength range for each spectrum based on the
        first and last unmasked pixel; interspersed masked regions are
        not considered.

        Args:
            flag (str or list): (**Optional**) Flags to consider when
                determining the wavelength range; see
                :func:`mangadap.util.bitmasks.BitMask.flagged`.

        Returns:
            numpy.ndarray: Two-element vector with wavelengths of the
            first and last valid pixels.
        """

        indx = numpy.where(numpy.invert(self.bitmask.flagged(self.hdu['MASK'].data, flag=flag)))
        return numpy.array([ numpy.amin(self.hdu['WAVE'].data[indx]),
                             numpy.amax(self.hdu['WAVE'].data[indx])]).astype(float)


    def _minimum_sampling(self):
        """
        Return the minimum sampling of the available wavelength vectors.
        
        Returns:
            float : minimum sampling of all (will just be one if the
            library has been processed) wavelength vectors.
        """
        if self.processed:
            return self.spectral_step

        nspec = self.hdu['WAVE'].data.shape[0]
        minimum_step = spectral_coordinate_step(self.hdu['WAVE'].data[0,:].ravel(),
                                     log=self.library['log10'])
        for i in range(1,nspec):
            step = spectral_coordinate_step(self.hdu['WAVE'].data[0,:].ravel(),
                                            log=self.library['log10'])
            if minimum_step > step:
                minimum_step = step

        return minimum_step


    def _rebin_masked(self, i, flag, fullRange, rmsk_lim=0.5):
        """
        Determine the mask value to adopt for a rebinned spectrum by
        rebinning the mask pixels, setting a masked pixel to unity, and
        an unmasked pixel to zero.  After rebinning, any pixel with a
        value of larger than *rmsk_lim* is masked; otherwise it is left
        unmasked.

        Although the function can be used with multiple flags, its
        intended use is to determine which pixels should be masked with
        a specific flag.

        Args:
            i (int): Index of the spectrum to be rebinned.
            flag (str or list): Flags to consider when determining which
                pixels to mask; see
                :func:`mangadap.util.bitmasks.BitMask.flagged`.
            fullRange (numpy.ndarray): Two-element array with the
                wavelength range for the rebinned spectrum.
            rmsk_lim (float): Limit of the rebinned mask value that is
                allowed before considering the pixel as masked.

        Returns:
            numpy.ndarray: Boolean array indicating which pixels in the
            rebinned spectrum should be masked.
        """
        unity = numpy.ma.MaskedArray(numpy.ones(self.hdu['WAVE'].data.shape[1], dtype=float),
                                     mask=self.bitmask.flagged(self.hdu['MASK'].data[i,:],
                                                               flag=flag))
        r = Resample(unity, x=self.hdu['WAVE'].data[i,:], inLog=self.library['log10'],
                     newRange=fullRange, newLog=self.log10_sampling, newdx=self.spectral_step)
        return r.outf < rmsk_lim

    def _modify_spectral_resolution(self):
        """
        Modify the spectral resolution to match the provided
        :attr:`sres`.

        Returns:
            float: Redshift used when matching the spectral resolution.
        """
        # Calculate to use as an offset of the match to the spectral
        # resolution.  Used to better match the spectral resolution to
        # at the *observed* wavelengths of the object spectrum to which
        # the TemplateLibrary will be fit.
        redshift = self.velocity_offset/astropy.constants.c.to('km/s').value

        # Mask wavelengths where the spectral resolution will have to be
        # extrapolated.
        sres_wave = self.sres.wave()
        wavelim = numpy.array([ sres_wave[0]/(1.+redshift), sres_wave[-1]/(1.+redshift) ])
        self.hdu = HDUList_mask_wavelengths(self.hdu, self.bitmask, 'SPECRES_EXTRAP', wavelim,
                                            invert=True)

#        i = 10
#        oldwave = numpy.copy(self.hdu['WAVE'].data[i,:]).ravel()
#        oldflux = numpy.copy(self.hdu['FLUX'].data[i,:]).ravel()
#        pyplot.plot(self.hdu['WAVE'].data[i,:], self.hdu['FLUX'].data[i,:]) 
#        pyplot.show()

        # Match the resolution of the templates.  ivar is returned, but
        # is always None because ivar is not provided to
        # match_spectral_resolution
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, 'Modifying spectral resolution ... ')
        print('min_sig_pix: ', self.min_sig_pix)
#        print(spectrum_velocity_scale(self.hdu['WAVE'].data))
        self.hdu['FLUX'].data, self.hdu['SPECRES'].data, self.hdu['SIGOFF'].data, res_mask, ivar = \
            match_spectral_resolution(self.hdu['WAVE'].data, self.hdu['FLUX'].data,
                                      self.hdu['SPECRES'].data, sres_wave/(1.+redshift),
                                      self.sres.sres(), min_sig_pix=self.min_sig_pix,
                                      no_offset=self.no_offset, log10=self.library['log10'],
                                      new_log10=True, quiet=self.quiet)
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, '... done')

        # Mask any pixels where the template resolution was too low to
        # match to the galaxy resolution
        self.hdu['MASK'].data[res_mask == 1] = \
            self.bitmask.turn_on(self.hdu['MASK'].data[res_mask == 1], 'SPECRES_LOW')

#        pyplot.plot(oldwave, oldflux)
#        pyplot.plot(self.hdu['WAVE'].data[i,:], self.hdu['FLUX'].data[i,:], 'g')
#        pyplot.show()

        return redshift


    def _process_library(self, wavelength_range=None, renormalize=True):
        """
        Process the template library for use in analyzing object
        spectra.   See :func:`process_template_library`.

        .. todo::

            - Add wavelength coordinate WCS information to the
              appropriate extension headers.
            - Include 'step' as an argument
            - Include limit on covering fraction used to mask data as an
              argument
        """
        # Convert to vacuum wavelengths
#        pyplot.plot(self.hdu['WAVE'].data[0,:], self.hdu['FLUX'].data[0,:]) 
        if not self.library['in_vacuum']:
            self.hdu['WAVE'].data = airtovac(self.hdu['WAVE'].data)

        # TODO: This conversion to vacuum wavelengths is non-linear.
        # This is accounted for in the resampling of the spectra (the
        # pixel size can be non-uniform), but *not* the resolution
        # matching.  This is okay for DR15 because we never resolution
        # match the data, but it will matter when we start to do that by
        # switching templates between the stellar and gas kinematics
        # fit.

#        pyplot.plot(self.hdu['WAVE'].data[0,:], self.hdu['FLUX'].data[0,:], 'g')
#        pyplot.show()

#        oldwave = numpy.copy(self.hdu['WAVE'].data[0,:]).ravel()
#        oldflux = numpy.copy(self.hdu['FLUX'].data[0,:]).ravel()

        #---------------------------------------------------------------
        # Modify the spectral resolution to a target function, if one
        # has been provided.
        redshift = 0.0 if self.sres is None else self._modify_spectral_resolution()

        #---------------------------------------------------------------
        # Resample the templates to a logarithmic binning step.
        #
        # Typically this will be run to match the sampling of the
        # template library to that of a set of galaxy spectra to be fit.
        # However, the sampling can be left alone and this part will
        # force all the spectra to have the same wavelength range.
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, 'Matching sampling ... ')

        # Number of angstroms per pixel
        ang_per_pix = numpy.array([angstroms_per_pixel(w, log=self.library['log10'], base=10.,
                                                       regular=self.library['in_vacuum'])
                                        for w in self.hdu['WAVE'].data])

        # TODO: I think this kludge was necessary because of a bug in
        # angstroms_per_pixel.  If regular was True and log was False,
        # it only returns a single integer, otherwise it returned a
        # vector.  I've changed angstroms_per_pixel, but have left the
        # kludge below for now for testing.

        #KHRR KLUDGE
        if len(ang_per_pix.shape) == 1:
            warnings.warn('Executing KLUDGE!  Bug fix insufficient.')
            tmp_ang_per_pix = numpy.empty(self.hdu['SPECRES'].data.shape)
            for ind in range(len(ang_per_pix)):
                tmp_ang_per_pix[ind,:] = ang_per_pix[ind]
            ang_per_pix = tmp_ang_per_pix

        # Number of pixels per resolution element
        pix_per_fwhm = numpy.ma.divide(self.hdu['WAVE'].data,
                                       self.hdu['SPECRES'].data * ang_per_pix)

        # The raw spectra are allowed to have wavelength ranges that
        # differ.  First, determine the wavelength range that encloses
        # all spectra.  Only ignore pixels that were flagged as having
        # no data.
        fullRange = self._wavelength_range(flag='NO_DATA') if wavelength_range is None \
                            else numpy.array(wavelength_range).astype(float)

        # Get the spectral step if it hasn't been set yet
        if self.spectral_step is None:
            self.spectral_step = self._minimum_sampling()

        # Get the number of pixels needed
        npix, _fullRange = resample_vector_npix(outRange=fullRange, dx=self.spectral_step,
                                                log=self.log10_sampling)

        # Any pixels without data after resampling are given a value
        # that is the minimum flux - 100 so that they can be easily
        # identified afterward.  The minimum flux is:
        min_flux = numpy.amin(self.hdu['FLUX'].data.ravel())
        # the observed pixels are
        no_data = self.bitmask.flagged(self.hdu['MASK'].data, flag='NO_DATA')

        # Now resample the spectra.  First allocate the arrays
        flux = numpy.zeros((self.ntpl, npix), dtype=numpy.float64)
        sres = numpy.zeros((self.ntpl, npix), dtype=numpy.float64)
        mask = numpy.zeros((self.ntpl, npix), dtype=self.bitmask.minimum_dtype())

        for i in range(self.ntpl):
            # Rebin the observed wavelength range
            r = Resample(self.hdu['FLUX'].data[i,:], mask=no_data[i,:],
                         x=self.hdu['WAVE'].data[i,:], inLog=self.library['log10'],
                         newRange=fullRange, newLog=self.log10_sampling,
                         newdx=self.spectral_step, step=False)
            
            # Save the result
            wave = r.outx
            flux[i,:] = r.outy
            indx = r.outf < 0.9
            flux[i,indx] = 0.0
            mask[i,indx] = self.bitmask.turn_on(mask[i,indx], 'NO_DATA')

            # Resample the spectral resolution by simple interpolation.
            sres[i,:] = interpolate.interp1d(self.hdu['WAVE'].data[i,:],
                                             self.hdu['SPECRES'].data[i,:], assume_sorted=True,
                                             fill_value='extrapolate')(wave)
            sres[i,indx] = 0.0

            # Recalculate the number of pixels per fwhm

            # Number of angstroms per pixel
            _ang_per_pix = angstroms_per_pixel(wave, log=self.log10_sampling, base=10.,
                                               regular=True)
            # Number of pixels per resolution element
            _pix_per_fwhm = numpy.ma.divide(wave, sres[i,:] * _ang_per_pix)

#            pyplot.plot(self.hdu['WAVE'].data[i,:-1], numpy.diff(self.hdu['WAVE'].data[i,:]))
#            pyplot.plot(self.hdu['WAVE'].data[i,:], pix_per_fwhm[i,:])
#            pyplot.plot(wave, _pix_per_fwhm)
#            pyplot.show()
#            exit()

            indx = _pix_per_fwhm < 2
            if numpy.sum(indx) > 0:
                warnings.warn('Resampling has cause resolution below the two pixel limit!')
                sres[i,indx] = wave[indx]/2./_ang_per_pix[indx]
                mask[i,indx] = self.bitmask.turn_on(mask[i,indx], 'SPECRES_2PIXEL')
            #-----------------------------------------------------------

            # Rebin the masks, bit-by-bit:
            # Pixels outside the wavelength limits
            indx = self._rebin_masked(i, 'WAVE_INVALID', fullRange, rmsk_lim=0.1)
            mask[i,indx] = self.bitmask.turn_on(mask[i,indx], 'WAVE_INVALID')
            # Pixels below the flux limit
            indx = self._rebin_masked(i, 'FLUX_INVALID', fullRange, rmsk_lim=0.1)
            mask[i,indx] = self.bitmask.turn_on(mask[i,indx], 'FLUX_INVALID')
            # Pixels that required an extrapolation of the spectral
            # resolution
            indx = self._rebin_masked(i, 'SPECRES_EXTRAP', fullRange, rmsk_lim=0.1)
            mask[i,indx] = self.bitmask.turn_on(mask[i,indx], 'SPECRES_EXTRAP')
            # Pixels that had a spectral resolution that was too low to
            # match the galaxy resolution
            indx = self._rebin_masked(i, 'SPECRES_LOW', fullRange, rmsk_lim=0.1)
            mask[i,indx] = self.bitmask.turn_on(mask[i,indx], 'SPECRES_LOW')
            
#            pyplot.plot(self.hdu['WAVE'].data[i,:], self.hdu['FLUX'].data[i,:])
#            pyplot.plot(wave, flux[i,:])
#            pyplot.show()
#            exit()

        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, '... done')
            log_output(self.loggers, 1, logging.INFO,
                       'After resampling (target): {0}'.format(self.spectral_step))
            log_output(self.loggers, 1, logging.INFO,
                       'After resampling (calculated): {0}'.format(spectral_coordinate_step(wave,
                                                                                        log=True)))
#        pyplot.plot(oldwave, oldflux)
#        pyplot.plot(wave, flux[0,:], 'g')
#        pyplot.show()

        # Normalize the templates to the mean flux value after excluding
        # any flagged pixels.
        if renormalize:
            indx = numpy.invert(self.bitmask.flagged(mask, flag=['NO_DATA', 'WAVE_INVALID', 
                                                                 'FLUX_INVALID']))
            if numpy.sum(indx) == 0:
                if not self.quiet:
                    warnings.warn('All pixels masked.  Unable to renormalize TemplateLibrary.')
                flux_norm = 1.0
            else:
                flux_norm = numpy.mean(flux[numpy.where(indx)])
            flux /= flux_norm
        else:
            flux_norm = 1.0
#        print(flux_norm)
#
#        pyplot.plot(oldwave, oldflux/flux_norm)
#        pyplot.plot(wave, flux[0,:], 'g')
#        pyplot.show()

        # Reset the HDUList object
        self._reset_hdu(wave, flux, mask, sres, self.hdu['SIGOFF'].data)

        # Update the header with the redshift offset, the flux
        # normalization, and flag the data as having been prepared for
        # fitting the DRP data
        self.hdu[0].header['ZGUESS'] = (redshift,'Redshift used when matching spectral resolution')
        self.hdu[0].header['FLXNRM'] = (flux_norm,'Flux normalization')
        self.hdu[0].header['TPLPROC'] = 1
        self.processed = True


    def file_name(self):
        """Return the name of the processed file."""
        return self.processed_file


    def file_path(self):
        """Return the full path to the processed file."""
        if self.directory_path is None or self.processed_file is None:
            return None
        return os.path.join(self.directory_path, self.processed_file)


    def read_raw_template_library(self, library_key=None, tpllib_list=None, dapsrc=None):
        """
        Read the identified template library.  If all the arguments are
        the default, the preset attributes from the initialization of
        the object are used.

        Args:
            library_key (str): (**Optional**) Keyword selecting the library
                to use.
            tpllib_list (list): (**Optional**) List of
                :class:`TemplateLibraryDef`
                objects that define the parameters required to read and
                interpret a template library.
            dapsrc (str): (**Optional**) Root path to the DAP source
                directory.  If not provided, the default is defined by
                :func:`mangadap.config.defaults.dap_source_dir`.

        """
        if library_key is not None:
            # Redefine the library
            self.library = self.define_library(library_key, tpllib_list=tpllib_list, dapsrc=dapsrc)
        self._read_raw()


    def process_template_library(self, library_key=None, tpllib_list=None, dapsrc=None, drpf=None,
                                 match_to_drp_resolution=True, velscale_ratio=None, sres=None,
                                 velocity_offset=0.0, min_sig_pix=0.0, no_offset=True,
                                 spectral_step=None, log=True, wavelength_range=None,
                                 renormalize=True, dapver=None, analysis_path=None,
                                 directory_path=None, processed_file=None, hardcopy=True,
                                 symlink_dir=None, clobber=False, loggers=None, quiet=False):
        r"""
        Process the template library for use in analyzing object
        spectra.  Primary steps are to:

            - Read the raw 1D fits files; see :func:`_read_raw`.

            - Convert the wavelengths to vacuum, if necessary; see
              :func:`mangadap.util.idlutils.airtovac`.

            - Mask wavelengths outside the rest wavelength range of the
              DRP spectrum, due to need to extrapolate these values; see
              :func:`mangadap.util.bitmasks.HDUList_mask_wavelengths`.

            - Match the spectral resolution of the template to that of
              the DRP spectra; see
              :func:`mangadap.util.resolution.match_spectral_resolution`.

            - Mask the template pixels where the spectral resolution was
              too low to match to the DRP spectra; see
              :func:`mangadap.util.bitmasks.BitMask.turn_on`.

            - Force a common wavelength range and sampling for all
              templates, where the sampling is forced to match the
              sampling of the DRP spectra; see
              :func:`mangadap.util.sampling.resample_vector_pix`.  The
              masks are appropriately resampled as well; see
              :func:`_rebin_masked`.

        .. warning::

            The routine **does not** check that that an existing
            processed file or the existing object has been processed
            using the same DRPFits, velocity_offset, velscale, or sres
            input.  If unsure, use clobber=True.
        
        Args:
            library_key (:obj:`str`):
                Keyword selecting the library to use.
            tpllib_list (:obj:`list`, optional):
                List of :class:`TemplateLibraryDef` objects that define
                the parameters required to read and interpret a template
                library.  The ``library_key`` must select one of the
                objects in this list.
            dapsrc (:obj:`str`, optional):
                Root path to the DAP source directory.  If not provided,
                the default is defined by
                :func:`mangadap.config.defaults.dap_source_dir`.
            drpf (:class:`mangadap.drpfits.DRPFits`, optional):
                DRP file (object) with which the template library is
                associated for analysis
            match_to_drp_resolution (:obj:`bool`, optional):
                Match the spectral resolution of the template library to
                the resolution provided by the
                :class:`mangadap.drpfits.DRPFits` object; the latter
                must be provided for this argument to have any use.
            velscale_ratio (:obj:`int`, optional):
                Resample the template spectra such that the ratio of the
                pixel scale in the provided
                :class:`mangadap.drpfits.DRPFits` object is this many
                times larger than the pixels in the resampled template
                spectrum.
            sres (:class:`mangadap.util.resolution.SpectralResolution`, optional):
                The object is used simply to access the spectral
                resolution and associated wavelength coordinate vector
                needed when matching the spectral resolution of the
                template library; this is used in place of the
                attributes in any provided DRP file object.
            velocity_offset (:obj:`float`, optional):
                Velocity offset to use when matching the spectral
                resolution between the template library and the galaxy
                spectra.
            spectral_step (:obj:`float`, optional):
                Target logarithmic (``log=True``) or linear
                (``log=False``) step in wavelength for the template
                library.
            log (:obj:`bool`, optional):
                Flag to force the library to be logarithmically sampled
                in wavelength.
            wavelength_range (array-like, optional):
                Force the template library to covert this spectral
                range.  Unobserved spectral regions will be flagged.
            renormalize (:obj:`bool`, optional):
                After processing, renormalize the flux to unity.
            dapver (:obj:`str`, optional):
                DAP version, which is used to define the default DAP
                analysis path.  Default is defined by
                :func:`mangadap.config.defaults.default_dap_version`
            analysis_path (:obj:`str`, optional):
                The path to the top level directory containing the DAP
                output files for a given DRP and DAP version.  Default
                is defined by
                :func:`mangadap.config.defaults.default_analysis_path`.
            directory_path (:obj:`str`, optional):
                The exact path to the processed template library file.
                Default is defined by
                :func:`mangadap.config.defaults.default_dap_common_path`.
            processed_file (:obj:`str`, optional):
                The name of the file containing the prepared template
                library output file.  The file should be found at
                :attr:`directory_path`/:attr:`processed_file`.  Default
                is defined by
                :func:`mangadap.config.defaults.default_dap_file_name`.
            hardcopy (:obj:`bool`, optional):
                Flag to keep a hardcopy of the processed template
                library.  Default is True.
            symlink_dir (:obj:`str`, optional):
                Create a symlink to the file in this directory.  Default
                is for no symlink.
            clobber (:obj:`bool`, optional):
                If :attr:`drpf` is define and ``process=True``, this
                will clobber any existing processed template library.

        Raises:
            ValueError:
                If the velocity scale, spectral resolution, or file name
                for the processed data are not define.

        .. todo::
            - type checking
            - If a DRP file is provided, the processing to a logarithmic
              binning is done by default (log=True).  But linearly
              sampled DRP data are available, so need to have
              :class:`mangadap.drpfits.DRPFits` return the spectral
              sampling type.
            - Documentation needs updating!

        """

        # Initialize the reporting
        if loggers is not None:
            self.loggers = loggers
        self.quiet = quiet

        # Ignore velscale_ratio if it is set to unity
        if velscale_ratio is not None and velscale_ratio == 1:
            velscale_ratio = None

        if drpf is not None:
            # Use the DRP file object to set the spectral resolution and
            # velocity scale to match to.  TODO: Need to update the
            # matching for 'DISP' in new LOGCUBEs?
            if drpf.hdu is None:
                if not self.quiet:
                    warnings.warn('DRP file previously unopened.  Reading now.')
                drpf.open_hdu()

            # If resolution matching is requested, get the spectral
            # resolution vector.
            #   - The procedure requires that the resolution be measured
            #     without the pixel convolution; however, the PRESPECRES
            #     extension is not available with all releases of the
            #     DRP.  So first try to get the prepixelized
            #     measurements, and then warn the user and get the
            #     pixelized measuremeents if that fails.
            #   - The spectral_resolution function returns a full array
            #     with the spectral resolution for each spectrum by
            #     default.  Calling it with median=True will only return
            #     the single vector with the median spectral resolution.
            if match_to_drp_resolution:
                try:
                    self.sres = drpf.spectral_resolution(ext='SPECRES', toarray=True, fill=True,
                                                         pre=True, median=True)
                except:
                    self.sres = drpf.spectral_resolution(ext='SPECRES', toarray=True, fill=True,
                                                         median=True)
                self.sres = SpectralResolution(drpf['WAVE'].data, self.sres, log10=True)
            else:
                self.sres = None

            # Set log sampling by default, but should instead query
            # DRPFits to check if the wavelengths are log binned!
            self.log10_sampling = True
            self.spectral_step = spectral_coordinate_step(drpf.hdu['WAVE'].data, log=True)
        # Use the provided input values
        else:
            self.sres = sres
            self.log10_sampling = log
            self.spectral_step = spectral_step

        self.velocity_offset = velocity_offset
        self.min_sig_pix = min_sig_pix
        self.no_offset = no_offset

        # Adjust for the velocity scale ratio between the template and
        # object data to be fit
        if not self.log10_sampling and velscale_ratio is not None:
            raise ValueError('velscale_ratio only valid with logarithmically sampled spectra.')
        if velscale_ratio is not None:
            self.spectral_step /= velscale_ratio
#        print(self.spectral_step)

        # Set the paths if possible
        directory_path = self.directory_path if directory_path is None else directory_path
        processed_file = self.processed_file if processed_file is None else processed_file
        if self._can_set_paths(directory_path, drpf, processed_file, quiet=True):
#            print('setting paths')
            self._set_paths(directory_path, dapver, analysis_path, drpf, processed_file)
        self.hardcopy = hardcopy
        self.symlink_dir = symlink_dir

        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, '-'*50)
            log_output(self.loggers, 1, logging.INFO, 'Template library output path: {0}'.format(
                                                                            self.directory_path))
            log_output(self.loggers, 1, logging.INFO, 'Template library output file: {0}'.format(
                                                                            self.processed_file))

        # Check that the path for or to the file is defined
        ofile = self.file_path()
        if self.hardcopy and ofile is None:
            raise ValueError('File path for output file is undefined!')

#        if not force and self.velocity_offset is not None \
#           and not numpy.isclose(velocity_offset, self.velocity_offset):
#            print('Forcing processing due to change in velocity offset.')
#            force = True

        # Read and use a pre-existing file
        if self.hardcopy and os.path.isfile(ofile) and not clobber:
            if not self.quiet:
                log_output(self.loggers, 1, logging.INFO, 'Reading existing file')
#            self.hdu = fits.open(ofile, checksum=self.checksum)
            self.hdu = DAPFitsUtil.read(ofile, checksum=self.checksum)
            self.file_list = self._get_file_list()
#            self.file_list = glob.glob(self.library['file_search'])
            self.ntpl = self.hdu['FLUX'].data.shape[0]
            self.processed = True
            # Make sure the symlink exists
            if self.symlink_dir is not None:
                create_symlink(ofile, self.symlink_dir, clobber=clobber)
            if not self.quiet:
                log_output(self.loggers, 1, logging.INFO, '-'*50)
            return

        # Warn the user that the file will be overwritten
        if self.hardcopy and os.path.isfile(ofile):
            if not self.quiet:
                warnings.warn('Overwriting existing file: {0}'.format(self.processed_file))
            os.remove(ofile)

        # Read the raw data
        self._read_raw()
        # Process the library
        self._process_library(wavelength_range=wavelength_range, renormalize=renormalize)
        # Write the fits file
        if self.hardcopy:
            DAPFitsUtil.write(self.hdu, self.file_path(), clobber=clobber, checksum=True,
                              symlink_dir=self.symlink_dir, loggers=self.loggers, quiet=self.quiet)
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, '-'*50)


    def single_spec_to_fits(self, i, ofile, clobber=True):
        """
        Write one of the template spectra to a 1D fits file.

        Output is a one-dimensional fits file with a single extension
        containing the flux measurements and headers keywords needed to
        determine the wavelength of each pixel.

        Args:
            i (int) : Index of the spectrum in the avaiable list to
                output.
            ofile (str) : Name of the file to write
            clobber (bool) : (**Optional**) Flag to clobber any existing
                file of the same name.

        Raises:
            ValueError: Raised if the selected index is not available.

        """

        if i < 0 or i >= self.ntpl:
            raise ValueError('Invalid index ({0})!  Library contains {1} spectra.'.format(i,
                                                                                        self.ntpl))
        wave = self.hdu['WAVE'].data if self.processed else self.hdu['WAVE'].data[i,:]
        log = self.log10_sampling if self.processed else self.library['log10']
        writefits_1dspec(ofile, (numpy.log10(wave[0]) if log else wave[0]), self.spectral_step,
                         self.hdu['FLUX'].data[i,:], clobber=clobber)


