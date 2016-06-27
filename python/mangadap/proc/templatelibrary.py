# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""

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

Two support classes are also provided.  One is a derived
:class:`mangadap.par.parset.ParSet` instance that provides the defining
parameters of a DAP template library.  The second is a derived
:class:`mangadap.util.bitmask.BitMask` instance that defines the
bitmasks for the template library spectra.

*License*:
    Copyright (c) 2015, SDSS-IV/MaNGA Pipeline Group
        Licensed under BSD 3-clause license - see LICENSE.rst

*Source location*:
    $MANGADAP_DIR/python/mangadap/proc/templatelibrary.py

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
            try:
                from configparser import ExtendedInterpolation
            except ImportError:
                warnings.warn('Unable to import ExtendedInterpolation!  Some configurations will fail!')
        else:
            try:
                from ConfigParser import ConfigParser
            except ImportError:
                warnings.warn('Unable to import ConfigParser!  Beware!')
            try:
                from ConfigParser import ExtendedInterpolation
            except ImportError:
                warnings.warn('Unable to import ExtendedInterpolation!  Some configurations will fail!')


        import glob
        import os
        import time
        import logging

        import numpy
        from scipy import sparse
        from scipy.interpolate import InterpolatedUnivariateSpline
        from astropy.io import fits
        import astropy.constants
        from pydl.goddard.astro import airtovac

        from ..util.bitmask import BitMask
        from ..par.parset import ParSet
        from ..config.defaults import default_dap_source, default_dap_common_path
        from ..config.defaults import default_dap_file_name
        from ..util.log import log_output
        from ..util.fileio import readfits_1dspec, read_template_spectrum, writefits_1dspec, write_hdu
        from ..util.instrument import resample_vector, resample_vector_npix, spectral_resolution
        from ..util.instrument import match_spectral_resolution, spectral_coordinate_step
        from .util import _select_proc_method, HDUList_mask_wavelengths

.. warning::

    Because of the use of the ``ExtendedInterpolation`` in
    `configparser.ConfigParser`_, :func:`available_template_libraries`
    is not python 2 compiliant.
    
*Class usage examples*:
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
    to, e.g., the DRP file; see :class:`mangadap.dapfile.dapfile`.

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
        from mangadap.config.defaults import default_dap_source
        from matplotlib import pyplot

        # Define the search string for the library
        search_str = default_dap_source()+'/data/stellar_templates/miles/*.fits'

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

*Revision history*:
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

.. _astropy.io.fits.hdu.hdulist.HDUList: http://docs.astropy.org/en/v1.0.2/io/fits/api/hdulists.html
.. _glob.glob: https://docs.python.org/3.4/library/glob.html
.. _configparser.ConfigParser: https://docs.python.org/3/library/configparser.html#configparser.ConfigParser
.. _pydl.goddard.astro.airtovac: http://pydl.readthedocs.io/en/stable/api/pydl.goddard.astro.airtovac.html#pydl.goddard.astro.airtovac
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
        warnings.warn('Unable to import configparser!  Beware!', ImportWarning)
    try:
        from configparser import ExtendedInterpolation
    except ImportError:
        warnings.warn('Unable to import ExtendedInterpolation!  Some configurations will fail!',
                      ImportWarning)
else:
    try:
        from ConfigParser import ConfigParser
    except ImportError:
        warnings.warn('Unable to import ConfigParser!  Beware!', ImportWarning)
    try:
        from ConfigParser import ExtendedInterpolation
    except ImportError:
        warnings.warn('Unable to import ExtendedInterpolation!  Some configurations will fail!',
                      ImportWarning)

import glob
import os
import time
import logging

import numpy
from scipy import sparse
from scipy.interpolate import InterpolatedUnivariateSpline
from astropy.io import fits
import astropy.constants
from pydl.goddard.astro import airtovac

from ..util.bitmask import BitMask
from ..par.parset import ParSet
from ..config.defaults import default_dap_source, default_dap_common_path
from ..config.defaults import default_dap_file_name
from ..util.log import log_output
from ..util.fileio import readfits_1dspec, read_template_spectrum, writefits_1dspec
from ..util.fileio import write_hdu, create_symlink
from ..util.instrument import resample_vector, resample_vector_npix, spectral_resolution
from ..util.instrument import match_spectral_resolution, spectral_coordinate_step
from .util import _select_proc_method, HDUList_mask_wavelengths

#from matplotlib import pyplot

__author__ = 'Kyle B. Westfall'
# Add strict versioning
# from distutils.version import StrictVersion

class TemplateLibraryDef(ParSet):
    """
    Class with parameters used to define the template library.
    See :class:`mangadap.par.parset.ParSet` for attributes.

    Args:
        key (str): Keyword to distinguish the template library.
        file_search (str): Search string used by glob to find the 1D
            fits spectra to include in the template library.
        fwhm (int or float): FWHM of the resolution element in
            angstroms.
        sres_ext (str): Extension in the fits files with measurements of
            the spectral resolution as a function of wavelength.
        in_vacuum (bool): Flag that the wavelengths of the spectra are
            in vacuum, not air.
        wave_limit (numpy.ndarray): 2-element array with the starting
            and ending wavelengths for the valid spectral range of the
            templates.
        lower_flux_limit (int or float): Minimum valid flux in the
            template spectra.
        log10 (bool): Flag that the template spectra have been binned
            logarithmically in wavelength.

    """
    def __init__(self, key, file_search, fwhm, sres_ext, in_vacuum, wave_limit, lower_flux_limit,
                 log10): 
        # Perform some checks of the input
        in_fl = [ int, float ]
        
        pars =   [ 'key', 'file_search', 'fwhm', 'sres_ext', 'in_vacuum',  'wave_limit',
                        'lower_flux_limit', 'log10' ]
        values = [   key,   file_search,   fwhm,   sres_ext,   in_vacuum,    wave_limit,
                          lower_flux_limit,   log10 ]
        dtypes = [   str,           str,  in_fl,        str,        bool, numpy.ndarray,
                                     in_fl,    bool ]

        ParSet.__init__(self, pars, values=values, dtypes=dtypes)


def validate_spectral_template_config(cnfg):
    """ 
    Validate the `configparser.ConfigParser`_ object that is meant to
    define a template library.

    Args:
        cnfg (`configparser.ConfigParser`_): Object meant to contain
            defining parameters of the template library needed by
            :class:`mangadap.proc.templatelibrary.TemplateLibraryDef`.

    Raises:
        KeyError: Raised if required keyword does not exist.
        ValueError: Raised if key has unacceptable value.

    """
    # Check for required keywords
    if 'key' not in cnfg.options('default'):
        raise KeyError('Keyword \'key\' must be provided.')
    if 'file_search' not in cnfg.options('default'):
        raise KeyError('Keyword \'file_search\' must be provided.')
    if 'fwhm' not in cnfg.options('default') and 'sres_ext' not in cnfg.options('default'):
        raise KeyError('Must provided keyword \'fwhm\' or \'sres_ext\'.')

    # Put in default values
    if 'in_vacuum' not in cnfg.options('default') or cnfg['default']['in_vacuum'] is None:
        cnfg['default']['in_vacuum'] = 'False'
    if 'wave_limit' not in cnfg.options('default') or cnfg['default']['wave_limit'] is None:
        cnfg['default']['wave_limit'] = 'None, None'
    if 'lower_flux_limit' not in cnfg.options('default') \
      or cnfg['default']['lower_flux_limit'] is None:
        cnfg['default']['lower_flux_limit'] = 'None'
    if 'log10' not in cnfg.options('default') or cnfg['default']['log10'] is None:
        cnfg['default']['log10'] = 'False'


def available_template_libraries(dapsrc=None):
    """
    Return the list of library keys, the searchable string of the 1D
    template library fits files for the template libraries available for
    use by the DAP, the FWHM of the libraries, and whether or not the
    wavelengths are in vacuum.

    The stellar template library files should be a list of 1D fits
    files, and be associated with one of the following library keys:

    +-----------------+------------+---------+-------------+-------+
    |                 |   Spectral |         |  Wavelength | Lower |
    |             KEY |  res (ang) | Vacuum? | Range (ang) | Limit |
    +=================+============+=========+=============+=======+
    |        M11MARCS |       2.73 |      No |        full |  None |
    +-----------------+------------+---------+-------------+-------+
    |       M11STELIB |       3.40 |      No |        full |  None |
    +-----------------+------------+---------+-------------+-------+
    |   M11STELIBZSOL |       3.40 |      No |        full |  None |
    +-----------------+------------+---------+-------------+-------+
    |       M11ELODIE |       0.55 |      No |      < 6795 |  None |
    +-----------------+------------+---------+-------------+-------+
    |        M11MILES |       2.54 |      No | 3550 - 7400 |  None |
    +-----------------+------------+---------+-------------+-------+
    |           MILES |       2.50 |      No |      < 7400 |  None |
    +-----------------+------------+---------+-------------+-------+
    |        MILESAVG |       2.50 |      No |      < 7400 |  None |
    +-----------------+------------+---------+-------------+-------+
    |       MILESTHIN |       2.50 |      No |      < 7400 |  None |
    +-----------------+------------+---------+-------------+-------+
    |         MILESHC |       2.50 |      No |      < 7400 |  None |
    +-----------------+------------+---------+-------------+-------+
    |          STELIB |       3.40 |      No |        full |   0.0 |
    +-----------------+------------+---------+-------------+-------+
    |         MIUSCAT |       2.51 |      No | 3480 - 9430 |  None |
    +-----------------+------------+---------+-------------+-------+
    |     MIUSCATTHIN |       2.51 |      No | 3480 - 9430 |  None |
    +-----------------+------------+---------+-------------+-------+

    .. warning::

        Function is currently only valid for Python 3.2 or greater!

    Args:
        dapsrc (str): (**Optional**) Root path to the DAP source
            directory.  If not provided, the default is defined by
            :func:`mangadap.config.defaults.default_dap_source`.

    Returns:
        list: A list of
        :func:`mangadap.proc.templatelibrary.TemplateLibraryDef`
        objects, each defining a separate template library.

    Raises:
        NotADirectoryError: Raised if the provided or default
            *dapsrc* is not a directory.
        OSError/IOError: Raised if no template configuration files could
            be found.
        KeyError: Raised if the template-library keywords are not all
            unique.
        NameError: Raised if either ConfigParser or
            ExtendedInterpolation are not correctly imported.  The
            latter is a *Python 3 only module*!

    .. todo::
        - Add backup function for Python 2.
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
    ini_files = glob.glob(dapsrc+'/python/mangadap/config/spectral_templates/*.ini')
    if len(ini_files) == 0:
        raise IOError('Could not find any configuration files in {0} !'.format(
                      dapsrc+'/python/mangadap/config/spectral_templates'))

    # Build the list of library definitions
    template_libraries = []
    for f in ini_files:
        # Read the config file
        cnfg = ConfigParser(os.environ, allow_no_value=True, interpolation=ExtendedInterpolation())
        cnfg.read(f)
        # Ensure it has the necessary elements to define the template
        # library
        validate_spectral_template_config(cnfg)
        # Convert wave_limit and lower_flux_limit to types acceptable by
        # TemplateLibraryDef
        wave_limit = numpy.array([ None if 'None' in e else float(e.strip()) \
                                        for e in cnfg['default']['wave_limit'].split(',') ])
        lower_flux_limit = None if cnfg['default']['lower_flux_limit'] is 'None' else \
                           cnfg['default'].getfloat('lower_flux_limit')
        # Append the definition of the template library 
        template_libraries += \
            [ TemplateLibraryDef(key=cnfg['default']['key'],
                                 file_search=cnfg['default']['file_search'],
                                 fwhm=cnfg['default'].getfloat('fwhm'),
                                 sres_ext=cnfg['default']['sres_ext'],
                                 in_vacuum=cnfg['default'].getboolean('in_vacuum'),
                                 wave_limit=wave_limit,
                                 lower_flux_limit=lower_flux_limit,
                                 log10=cnfg['default'].getboolean('log10') )
            ]

    # Check the keywords of the libraries are all unique
    if len(numpy.unique( numpy.array([tpl['key'] for tpl in template_libraries]) )) \
            != len(template_libraries):
        raise KeyError('Template-library keywords are not all unique!')

    # Return the default list of template libraries
    return template_libraries


class TemplateLibraryBitMask(BitMask):
    """
    Derived class that specifies the mask bits for the template library
    data.  See :class:`mangadap.util.bitmask.BitMask` for attributes.

    A list of the bits and meanings are provided by the base class
    function :func:`mangadap.util.bitmask.BitMask.info`; i.e.,::

        from mangadap.proc.templatelibrary import TemplateLibraryBitMask
        t = TemplateLibraryBitMask()
        t.info()

    """
    def __init__(self, dapsrc=None):
        dapsrc = default_dap_source() if dapsrc is None else str(dapsrc)
        BitMask.__init__(self, ini_file=os.path.join(dapsrc, 'python', 'mangadap', 'config',
                                                     'bitmasks', 'spectral_template_bits.ini'))
#        t = BitMask.from_ini_file(os.path.join(dapsrc, 'python', 'mangadap', 'config',
#                                                  'bitmasks', 'spectral_template_bits.ini'))
#        BitMask.__init__(self, t.keys(), t.descr)


class TemplateLibrary:
    r"""
    Object used to read, store, and prepare template libraries for use
    in analyzing object spectra.

    The default list of available libraries provided by the MaNGA DAP
    defined in :func:`available_template_libraries`.  The user can
    provide their own library for use with this class provided they are
    contained in 1D fits spectra, sampled linearly in wavelength with
    the wavelength coordinates available via the WCS keywords (CRPIX1,
    CRVAL1, CDELT1), and they have an appropriately defined spectral
    resolution (FWHM in angstroms that is constant as a function of
    wavelength).  See :class:`TemplateLibraryDef` and
    :func:`_build_raw_hdu`.

    The class is optimized for use in analyzing MaNGA DRP files;
    however, one can provide the necessary information so that the class
    can be used with a non-DRP spectrum.  In the latter case, the user
    must supply the velocity scale of the pixel for the logarithmically
    resampled template library, and a
    :class:`mangadap.util.instrument.spectral_resolution` object the
    defines the instrumental resolution of the spectrum/spectra to be
    analyzed.

    .. todo::
        - below is out of date.

    On initialization, if the DRP file object is not provided (is None),
    the default behavior is to read the raw template library if
    read=True.  If the DRP file is provided, the routine will check for
    the resolution matched fits file; if it doesn't exist and read is
    True, it will prepare the template library for use in analyzing the
    DRP file and write the prepared library file.  If clobber=True, the
    preparation and writing of the template library will be done even if
    the library already exists.

    Args:
        library_key (str): Keyword selecting the library to use.
        tpllib_list (list): (**Optional**) List of
            :class:`TemplateLibraryDef` objects that define the
            parameters required to read and interpret a template
            library.  The *library_key* must select one of the objects
            in this list.
        dapsrc (str): (**Optional**) Root path to the DAP source
            directory.  If not provided, the default is defined by
            :func:`mangadap.config.defaults.default_dap_source`.
        drpf (:class:`mangadap.drpfits.DRPFits`): (**Optional**) DRP
            file (object) with which the template library is associated
            for analysis
        match_to_drp_resolution (bool): (**Optional**) Match the
            spectral resolution of the template library to the
            resolution provided by the :class:`mangadap.drpfits.DRPFits`
            object; the latter must be provided for this argument to
            have any use.
        velscale_ratio (int): (**Optional**) Resample the template
            spectra such that the ratio of the pixel scale in the
            provided :class:`mangadap.drpfits.DRPFits` object is this
            many times larger than the pixels in the resampled template
            spectrum.
        sres (:class:`mangadap.util.instrument.spectral_resolution`):
            (**Optional**) The object is used simply to access the
            spectral resolution and associated wavelength coordinate
            vector needed when matching the spectral resolution of the
            template library; this is used in place of the attributes in
            any provided DRP file object.
        velocity_offset (float): (**Optional**) Velocity offset to use
            when matching the spectral resolution between the template
            library and the galaxy spectra.
        spectral_step (float) : (**Optional**) Target logarithmic
            (*log*=True) or linear (*log*=False) step in wavelength for
            the template library.
        log (bool) : (**Optional**) Flag to force the library to be
            logarithmically sampled in wavelength.
        wavelength_range (list or numpy.ndarray) : (**Optional**) Force
            the template library to covert this spectral range.
            Unobserved spectral regions will be flagged.
        renormalize (bool) : (**Optional**) After processing,
            renormalize the flux to unity.
        dapver (str): (**Optional**) DAP version, which is used to
            define the default DAP analysis path.  Default is defined by
            :func:`mangadap.config.defaults.default_dap_version`
        analysis_path (str): (**Optional**) The path to the top level
            directory containing the DAP output files for a given DRP
            and DAP version.  Default is defined by
            :func:`mangadap.config.defaults.default_analysis_path`.
        directory_path (str): (**Optional**) The exact path to the
            processed template library file.  Default is defined by
            :func:`mangadap.config.defaults.default_dap_common_path`.
        processed_file (str): (**Optional**) The name of the file
            containing the prepared template library output file.  The
            file should be found at
            :attr:`directory_path`/:attr:`processed_file`.  Default is
            defined by
            :func:`mangadap.config.defaults.default_dap_file_name`.
        process (bool): (**Optional**) If :attr:`drpf` is defined and
            the prepared template library does not exist, this will
            process the template library in preparation for use in
            fitting the provided DRP file.
        hardcopy (bool): (**Optional**) Flag to keep a hardcopy of the
            processed template library.  Default is True.
        symlink_dir (str): (**Optional**) Create a symlink to the file
            in this directory.  Default is for no symlink.
        clobber (bool): (**Optional**) If :attr:`drpf` is define and
            *process* is True, this will clobber any existing processed
            template library.

    Attributes:
        version (str): Version number
        bitmask (BitMask): A BitMask object used to toggle mask values;
            see :func:`TemplateLibraryBitMask`.
        library (:class:`TemplateLibraryDef`): Parameter set required to
            read and prepare the library.
        file_list (list): The list of files found using `glob.glob`_ and
            :attr:`file_search`.
        ntpl (int): Number of template spectra in the library
        drpf (:class:`mangadap.drpfits.DRPFits`): DRP file (object) with
            which the template library is associated for analysis
        sres (:class:`mangadap.util.instrument.spectral_resolution`):
            The object is used simply to access the spectral resolution
            and associated wavelength coordinate vector needed when
            matching the spectral resolution of the template library;
            this is used in place of the attributes in any provided DRP
            file object.
        velocity_offset (float): Velocity offset to use when matching
            the spectral resolution between the template library and the
            galaxy spectra.
        spectral_step (float) : Target logarithmic
            (:attr:`log10_sampling`=True) or linear
            (:attr:`log10_sampling`=False) step in wavelength for the
            template library.
        log10_sampling (bool): Flag that the processed template library is
            logarithmically sampled in wavelength.
        directory_path (str): The exact path to the processed template
            library file.  Default is defined by
            :func:`mangadap.config.defaults.default_dap_common_path`.
        processed_file (str): The name of the file containing (to
            contain) the prepared template library output file.  The
            file should be found at
            :attr:`directory_path`/:attr:`processed_file`.
        processed (bool): Flag that the template library has been
            prepared for use in the DAP.
        hardcopy (bool): Flag to keep a hardcopy of the processed
            template library.
        symlink_dir (str): Symlink created to the file in this directory
        hdu (`astropy.io.fits.hdu.hdulist.HDUList`_): HDUList read from
            the DAP file

    .. todo::
        - Only works with DRP files that have log-linear wavelength
          binning!
        - Allow to process, where process is just to change the
          sampling or the resolution (not necessarily both).
        - Need to make this more general, removing all dependence on
          DRPFits object.  This would simplify the functionality to
          change how the resolution and sampling matching is specified.

    """
    def __init__(self, library_key, tpllib_list=None, dapsrc=None, drpf=None,
                 match_to_drp_resolution=True, velscale_ratio=None, sres=None, velocity_offset=0.0,
                 spectral_step=None, log=True, wavelength_range=None, renormalize=True, dapver=None,
                 analysis_path=None, directory_path=None, processed_file=None, read=True,
                 process=True, hardcopy=True, symlink_dir=None, clobber=False, checksum=False,
                 loggers=None, quiet=False):

        self.version = '2.1'
        self.loggers = loggers
        self.quiet = quiet

        # Define the TemplateLibraryBitMask object
        self.bitmask = TemplateLibraryBitMask(dapsrc=dapsrc)

        # Define the properties needed to modify the spectral resolution
        self.sres = None
        self.velocity_offset = None

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
        self._define_library(library_key, tpllib_list=tpllib_list, dapsrc=dapsrc)

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
                                      velocity_offset=velocity_offset, spectral_step=spectral_step,
                                      log=log, wavelength_range=wavelength_range,
                                      renormalize=renormalize, dapver=dapver,
                                      analysis_path=analysis_path, directory_path=directory_path,
                                      processed_file=processed_file, hardcopy=hardcopy,
                                      symlink_dir=symlink_dir, clobber=clobber, loggers=loggers,
                                      quiet=quiet)


    def __del__(self):
        """
        Deconstruct the template library object by ensuring that the
        fits file is properly closed.
        """
        if self.hdu is None:
            return
        self.hdu.close()
        self.hdu = None


    def __getitem__(self, key):
        return self.hdu[key]


    def _define_library(self, library_key, tpllib_list=None, dapsrc=None):
        """
        Select the library from the provided list.  Used to set
        :attr:`library`; see
        :func:`mangadap.proc.util._select_proc_method`.

        Args:
            library_key (str): Keyword of the selected library.
                Available libraries are proved by
                :func:`available__template_libraries`
            tpllib_list (list): (**Optional**) List of 
                :class:`TemplateLibraryDef'
                objects that define the parameters required to read and
                interpret a template library.
            dapsrc (str): (**Optional**) Root path to the DAP source
                directory.  If not provided, the default is defined by
                :func:`mangadap.config.defaults.default_dap_source`.
        """
        # Get the details of the selected template library
        self.library = _select_proc_method(library_key, TemplateLibraryDef, method_list=tpllib_list,
                                           available_func=available_template_libraries,
                                           dapsrc=dapsrc)


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


    def _read_raw(self):
        """
        Read the 'raw' versions of the template library; i.e., the
        library before it has been resolution and sampling matched to a
        DRP file.
        """
        self.file_list = glob.glob(self.library['file_search'])
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
        :func:`mangadap.util.instrument.GaussianKernelDifference`).

        .. warning::

            Currently no errors are saved because none are expected for
            the template libraries.

        Args:
            npix (int): Number of spectral channels for the output
                arrays

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
                wave_, flux_ = readfits_1dspec(self.file_list[i], log10=self.library['log10'])
                wave[i,0:wave_.size] = numpy.copy(wave_)
                flux[i,0:wave_.size] = numpy.copy(flux_)
                # Set the spectral resolution
                sres = wave/self.library['fwhm']
            else:
                wave_, flux_, sres_ = read_template_spectrum(self.file_list[i],
                                                             sres_ext=self.library['sres_ext'],
                                                             log10=self.library['log10'])
                wave[i,0:wave_.size] = numpy.copy(wave_)
                flux[i,0:wave_.size] = numpy.copy(flux_)
                sres[i,0:wave_.size] = numpy.copy(sres_)

            if wave_.size != npix:
                mask[i,wave_.size:] = self.bitmask.turn_on(mask[i,wave_.size:],'NO_DATA')
            if self.library['lower_flux_limit'] is not None:
                indx = numpy.invert( flux_ > self.library['lower_flux_limit'] )
                mask[i,indx] = self.bitmask.turn_on(mask[i,indx], 'FLUX_INVALID')
            if self.library['wave_limit'][0] is not None:
                indx = wave[i,:].ravel() < self.library['wave_limit'][0]
                mask[i,indx] = self.bitmask.turn_on(mask[i,indx], 'WAVE_INVALID')
            if self.library['wave_limit'][1] is not None:
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
                :func:`mangadap.util.instrument.GaussianKernelDifference`). 

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

        self.hdu['PRIMARY'].header['VDAPTPL'] = (self.version, 'Version of DAP template reader')
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
                             numpy.amax(self.hdu['WAVE'].data[indx])])


    def _minimum_sampling(self):
        """
        Return the minimum sampling of the available wavelength vetors.
        
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
            tuple: The indices of the pixels in the rebinned spectrum
            that should be masked.
        """

        mask_ex = self.bitmask.flagged(self.hdu['MASK'].data[i,:], flag=flag).astype(numpy.float64)
        wave, mask_ex = resample_vector(mask_ex, xRange=[self.hdu['WAVE'].data[i,0],
                                                         self.hdu['WAVE'].data[i,-1]],
                                        inLog=self.library['log10'], newRange=fullRange,
                                        newLog=self.log10_sampling, dx=self.spectral_step,
                                        conserve=False) # flat=True!

        return numpy.where(mask_ex > rmsk_lim)


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

#        oldwave = numpy.copy(self.hdu['WAVE'].data[0,:]).ravel()
#        oldflux = numpy.copy(self.hdu['FLUX'].data[0,:]).ravel()
#        pyplot.plot(self.hdu['WAVE'].data[0,:], self.hdu['FLUX'].data[0,:]) 
#        pyplot.show()

        # Match the resolution of the templates.  ivar is returned, but
        # is always None because ivar is not provided to
        # match_spectral_resolution
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, 'Modifying spectral resolution ... ')
        self.hdu['FLUX'].data, self.hdu['SPECRES'].data, self.hdu['SIGOFF'].data, res_mask, ivar = \
            match_spectral_resolution(self.hdu['WAVE'].data, self.hdu['FLUX'].data,
                                      self.hdu['SPECRES'].data, sres_wave/(1.+redshift),
                                      self.sres.sres(), min_sig_pix=0.0,
                                      log10=self.library['log10'], new_log10=True)
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, '... done')

#        pyplot.plot(oldwave, oldflux)
#        pyplot.plot(self.hdu['WAVE'].data[0,:], self.hdu['FLUX'].data[0,:], 'g') 
#        pyplot.show()

        # Mask any pixels where the template resolution was too low to
        # match to the galaxy resolution
        self.hdu['MASK'].data[res_mask == 1] = \
            self.bitmask.turn_on(self.hdu['MASK'].data[res_mask == 1], 'SPECRES_LOW')

#        pyplot.plot(oldwave, oldflux)
#        pyplot.plot(self.hdu['WAVE'].data[0,:], self.hdu['FLUX'].data[0,:], 'g')
#        pyplot.show()

        return redshift


    def _process_library(self, wavelength_range=None, renormalize=True):
        """
        Process the template library for use in analyzing object
        spectra.   See :func:`process_template_library`.

        .. todo::

            - Add wavelength coordinate WCS information to the
              appropriate extension headers.

        """
        # Convert to vacuum wavelengths
#        pyplot.plot(self.hdu['WAVE'].data[0,:], self.hdu['FLUX'].data[0,:]) 

        if not self.library['in_vacuum']:
            self.hdu['WAVE'].data \
                    = airtovac(self.hdu['WAVE'].data.ravel()).reshape(self.hdu['WAVE'].data.shape)

#        pyplot.plot(self.hdu['WAVE'].data[0,:], self.hdu['FLUX'].data[0,:], 'g')
#        pyplot.show()

#        oldwave = numpy.copy(self.hdu['WAVE'].data[0,:]).ravel()
#        oldflux = numpy.copy(self.hdu['FLUX'].data[0,:]).ravel()

        ################################################################
        # Modify the spectral resolution to a target function, if one
        # has been provided.
        redshift = 0.0 if self.sres is None else self._modify_spectral_resolution()
#        warnings.warn('Running test that does not match the spectral resolution!')
#        redshift = 0.0

        ################################################################
        # Resample the templates to a logarithmic binning step.
        #
        # Typically this will be run to match the sampling of the
        # template library to that of a set of galaxy spectra to be fit.
        # However, the sampling can be left alone and this part will
        # force all the spectra to have the same wavelength range.

        # The raw spectra are allowed to have wavelength ranges that
        # differ.  First, determine the wavelength range that encloses
        # all spectra.  Only ignore pixels that were flagged as having
        # no data.
        fullRange = self._wavelength_range(flag='NO_DATA') if wavelength_range is None \
                            else numpy.array(wavelength_range)
        # Get the spectral step if it hasn't been set yet
        if self.spectral_step is None:
            self.spectral_step = self._minimum_sampling()

        # Get the number of pixels needed
#        print(fullRange)
        npix, _fullRange = resample_vector_npix(outRange=fullRange, dx=self.spectral_step,
                                               log=self.log10_sampling)
#        print(fullRange, _fullRange, self.spectral_step, npix)

        # Any pixels without data after resampling are given a value
        # that is the minimum flux - 100 so that they can be easily
        # identified afterward.  The minimum flux is:
        min_flux = numpy.amin(self.hdu['FLUX'].data.ravel())
        # the observed pixels are
        observed = numpy.invert(self.bitmask.flagged(self.hdu['MASK'].data, flag='NO_DATA'))
        
        # Now resample the spectra.  First allocate the arrays
        flux = numpy.zeros((self.ntpl, npix), dtype=numpy.float64)
        sres = numpy.zeros((self.ntpl, npix), dtype=numpy.float64)
        mask = numpy.zeros((self.ntpl, npix), dtype=self.bitmask.minimum_dtype())
        if not self.quiet and self.loggers is not None:
            log_output(self.loggers, 1, logging.INFO, 'Matching sampling ... ')
        for i in range(0,self.ntpl):
            # Observed wavelengths
            wave_in = self.hdu['WAVE'].data[i,observed[i,:]].ravel()
            # Rebin the observed wavelength range
            wave, flux[i,:] = resample_vector(self.hdu['FLUX'].data[i,observed[i,:]].ravel(),
                                              xRange=[wave_in[0], wave_in[-1]],
                                              inLog=self.library['log10'], newRange=fullRange,
                                              newLog=self.log10_sampling, dx=self.spectral_step,
                                              ext_value=min_flux-100., conserve=False, flat=False)

#            pyplot.step(oldwave, oldflux, where='mid')
#            pyplot.step(wave, flux[i,:], 'g', where='mid')
#            pyplot.show()

            # Find the unobserved pixels, set them to have 0. flux, and
            # flag them as having no data
            indx = numpy.where(flux[i,:] < min_flux-10.)
            flux[i,indx] = 0.0
            mask[i,indx] = self.bitmask.turn_on(mask[i,indx], 'NO_DATA')

            # Resample the spectral resolution by simple interpolation.
            # Select the good pixels
            indx = numpy.where(numpy.invert(flux[i,:] < min_flux-10.))
            # define the interpolator (uses linear interpolation; k=1)
            interpolator = InterpolatedUnivariateSpline(self.hdu['WAVE'].data[i,:].ravel(),
                                                        self.hdu['SPECRES'].data[i,:].ravel(), k=1)
            # And then interpolate
            sres[i,indx] = interpolator(wave[indx])

            # Finally, rebin the masks:
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
            
#            pyplot.plot(oldwave, oldflux)
#            pyplot.plot(wave, flux[i,:], 'g')
#            pyplot.show()

        if not self.quiet and self.loggers is not None:
            log_output(self.loggers, 1, logging.INFO, '... done')

#        pyplot.plot(oldwave, oldflux)
#        pyplot.plot(wave, flux[0,:], 'g')
#        pyplot.show()

        # Normalize the templates to the mean flux value after excluding
        # any flagged pixels.
        if renormalize:
            indx = numpy.invert(self.bitmask.flagged(mask))
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


#    def dlogLam(self, quiet=True):
#        """
#        Return the logarithmic wavelength sampling of the spectra.
#
#        Always assumes the base of the logarithm is 10!
#
#        Returns:
#            float : Logarithmic wavelength sampling.  Returns None if
#                the velocity scale is not defined.
#
#        """
#        if self.velscale is None:
#            if not quiet:
#                warnings.warn('Velocity scale undefined!')
#            return None
#
#        return self.velscale / numpy.log(10.0) / astropy.constants.c.to('km/s').value


#    def velocity_step(self, quiet=True):
#        """
#        If the library is logarithmically sampled, return the pixel
#        scale in km/s.
#
#        Always assumes the base of the logarithm is 10!
#
#        Args:
#            quiet (bool) : Suppress warning that library is linearly sampled.
#
#        Returns:
#            float: Pixel scale in km/s for logarithmic wavelength
#            sampling.  Returns None if the sampling is linear.
#
#        Raises:
#            ValueError: Raised if library has not been processed.
#
#        .. todo::
#
#            - Allow to provide the velocity step, even if the library
#              has not been processed.  Would then return a vector with
#              the sampling of each spectrum.
#
#        """
#        if not self.processed:
#            raise ValueError('Template has not been processed!')
#        if not self.log10_sampling:
#            if not quiet:
#                warnings.warn('Template library is not logarithmically sampled in wavelength!')
#            return None
#
#        return self.spectral_step * numpy.log(10.0) * astropy.constants.c.to('km/s').value


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
                :func:`mangadap.config.defaults.default_dap_source`.

        """
        if library_key is not None:
            # Redefine the library
            self._define_library(library_key, tpllib_list=tpllib_list, dapsrc=dapsrc)
        self._read_raw()


    def process_template_library(self, library_key=None, tpllib_list=None, dapsrc=None, drpf=None,
                                 match_to_drp_resolution=True, velscale_ratio=None, sres=None,
                                 velocity_offset=0.0, spectral_step=None, log=True,
                                 wavelength_range=None, renormalize=True, dapver=None,
                                 analysis_path=None, directory_path=None, processed_file=None,
                                 hardcopy=True, symlink_dir=None, clobber=False, loggers=None,
                                 quiet=False):
        """

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
              :func:`mangadap.util.instrument.match_spectral_resolution`.

            - Mask the template pixels where the spectral resolution was
              too low to match to the DRP spectra; see
              :func:`mangadap.util.bitmasks.BitMask.turn_on`.

            - Force a common wavelength range and sampling for all
              templates, where the sampling is forced to match the
              sampling of the DRP spectra; see
              :func:`mangadap.util.instrument.resample_vector_pix`.  The masks
              are appropriately resampled as well; see
              :func:`_rebin_masked`.

        .. warning::

            The routine **does not** check that that an existing
            processed file or the existing object has been processed
            using the same DRPFits, velocity_offset, velscale, or sres
            input.  If unsure, use clobber=True.
        
        Args:

            library_key (str): (**Optional**) Keyword selecting the library
                to use; default is to use existing :attr:`library`.
            tpllib_list (list): (**Optional**) List of
                :class:`TemplateLibraryDef`
                objects that define the parameters required to read and
                interpret a template library.  Input ignored if
                *library_key* is None.
            dapsrc (str): (**Optional**) Root path to the DAP source
                directory.  If not provided, the default is defined by
                :func:`mangadap.config.defaults.default_dap_source`.
                Input ignored if *library_key* is None.
            drpf (:class:`mangadap.drpfits.DRPFits`): (**Optional**) DRP
                file (object) with which the template library is
                associated for analysis.  **If not provided**, the user
                must define *velscale* and *sres* such that the library
                can be processed; and the user must provide
                *directory_path* and *processed_file* such that the
                output file can be written.
            match_to_drp_resolution (bool): (**Optional**) Match the
                spectral resolution of the template library to the
                resolution provided by the
                :class:`mangadap.drpfits.DRPFits` object; the latter
                must be provided for this argument to have any use.
            velscale_ratio (int): (**Optional**) Resample the template
                spectra such that the ratio of the pixel scale in the
                provided :class:`mangadap.drpfits.DRPFits` object is
                this many times larger than the pixels in the resampled
                template spectrum.
            sres (:class:`mangadap.util.instrument.spectral_resolution`):
                (**Optional**) The object is used simply to access the
                spectral resolution and associated wavelength coordinate
                vector needed when matching the spectral resolution of
                the template library.  This takes prededence over the
                values provided by the DRP file object.
            velocity_offset (float): (**Optional**) Velocity offset to use
                when matching the spectral resolution between the
                template library and the galaxy spectra.
            spectral_step (float) : (**Optional**) Target logarithmic
                (*log*=True) or linear (*log*=False) step in wavelength for
                the template library.
            log (bool) : (**Optional**) Flag to force the library to be
                logarithmically sampled in wavelength.
            wavelength_range (list or numpy.ndarray) : (**Optional**) Force the
                template library to covert this spectral range.  Unobserved
                spectral regions will be flagged.
            renormalize (bool) : (**Optional**) After processing, renormalize
                the flux to unity.
            dapver (str): (**Optional**) DAP version, which is used to
                define the default DAP analysis path.  Default is
                defined by
                :func:`mangadap.config.defaults.default_dap_version`
            analysis_path (str): (**Optional**) The path to the top level
                directory containing the DAP output files for a given
                DRP and DAP version.  Default is defined by
                :func:`mangadap.config.defaults.default_analysis_path`.
            directory_path (str): (**Optional**) The exact path for the
                processed template library file.  Default is defined by
                :func:`mangadap.config.defaults.default_dap_common_path`.
            processed_file (str): (**Optional**) The name of the file
                containing (to contain) the prepared template library
                output file.  The file should be found at
                :attr:`directory_path`/:attr:`processed_file`.  Default
                is defined by
                :func:`mangadap.config.defaults.default_dap_file_name`.
            hardcopy (bool): (**Optional**) Flag to keep a hardcopy of
                the processed template library.  Default is True.
            symlink_dir (str): (**Optional**) Create a symlink to the
                file in this directory.  Default is for no symlink.
            clobber (bool): (**Optional**) Clobber any existing processed
                library.

        Raises:
            ValueError: If the velocity scale, spectral resolution, or
                file name for the processed data are not define.

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

        # Use the DRP file object to set the spectral resolution and
        # velocity scale to match to.
        if drpf is not None:
            if drpf.hdu is None:
                if not self.quiet:
                    warnings.warn('DRP file previously unopened.  Reading now.')
                drpf.open_hdu()
            self.sres = spectral_resolution(drpf.hdu['WAVE'].data, drpf.hdu['SPECRES'].data,
                                            log10=True) if match_to_drp_resolution else None
#            self.velscale = spectrum_velocity_scale(drpf.hdu['WAVE'].data, log10=True)
            # Set this by default, but need DRPFits to return if the
            # file is logarithmically sampled!
            self.log10_sampling = True
            self.spectral_step = spectral_coordinate_step(drpf.hdu['WAVE'].data, log=True)
            if velscale_ratio is not None:
                self.spectral_step /= velscale_ratio
#            print(self.spectral_step)
        # Use the provided input values
        else:
            self.sres = sres
            self.log10_sampling = log
            self.spectral_step = spectral_step
#            self.velscale = velscale
        self.velocity_offset = velocity_offset

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
        if ofile is None:
            raise ValueError('File path for output file is undefined!')

#        if not force and self.velocity_offset is not None \
#           and not numpy.isclose(velocity_offset, self.velocity_offset):
#            print('Forcing processing due to change in velocity offset.')
#            force = True

        # Read and use a pre-existing file
        if os.path.isfile(ofile) and not clobber:
            if not self.quiet:
                log_output(self.loggers, 1, logging.INFO, 'Reading existing file')
            self.hdu = fits.open(ofile, checksum=self.checksum)
            self.file_list = glob.glob(self.library['file_search'])
            self.ntpl = self.hdu['FLUX'].data.shape[0]
            self.processed = True
            # Make sure the symlink exists
            if self.symlink_dir is not None:
                create_symlink(ofile, self.symlink_dir, clobber=clobber)
            if not self.quiet:
                log_output(self.loggers, 1, logging.INFO, '-'*50)
            return

        # Warn the user that the file will be overwritten
        if os.path.isfile(ofile):
            if not self.quiet:
                warnings.warn('Overwriting existing file: {0}'.format(self.processed_file))
            os.remove(ofile)

        # Read the raw data
        self._read_raw()
        # Process the library
        self._process_library(wavelength_range=wavelength_range, renormalize=renormalize)
        # Write the fits file
        if self.hardcopy:
            write_hdu(self.hdu, self.file_path(), clobber=clobber, checksum=True,
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


