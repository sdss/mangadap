# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
r"""
Provides a set of functions that define and return defaults used by the
MaNGA DAP, such as paths and file names.

*Source location*:
    $MANGADAP_DIR/python/mangadap/util/defaults.py

*Python2/3 compliance*::

    from __future__ import division
    from __future__ import print_function
    from __future__ import absolute_import
    from __future__ import unicode_literals
    
    import sys
    if sys.version > '3':
        long = int

*Imports*::

    import os.path
    from os import environ
    import glob
    import numpy
    from mangadap.par.database_definitions import TemplateLibraryDef

*Revision history*:
    | **23 Apr 2015**: Original implementation by K. Westfall (KBW)
    | **19 May 2015**: (KBW) Documentation and Sphinx tests
    | **20 May 2015**: (KBW) TODO: Change the default names of the par
        and output fits files
    | **15 Jun 2015**: (KBW) Added default functions moved from
        :class:`mangadap.drpfile`
    | **16 Jun 2015**: (KBW) Added wavelength limits and lower flux
        limit to :func:`available_template_libraries`
    | **17 Jun 2015**: (KBW) Moved declarations of template library keys
        to its own function: :func:`default_template_library_keys', and
        edited :func:`available_template_libraries` accordingly
    | **27 Aug 2015**: (KBW) Changed the name of the plan file; added
        :func:`default_dap_file_root` based on file_root() from
        :class:`mangadap.survey.rundap`.
    | **28 Aug 2015**: (KBW) Added
        :func:`mangadap.util.defaults.default_manga_fits_root` and
        :func:`mangadap.util.defaults.default_cube_covariance_file`
    | **07 Oct 2015**: (KBW) Adjusted for changes to naming of the
        template library database definitions.  Added M11-STELIB-ZSOL
        library.
    | **29 Jan 2016**: (KBW) Changed :func:`default_template_libraries`
        to use configparser ini files to define each template library.
    | **03 Feb 2016**: (KBW) Added checks for required environmental
        variables.
    | **17 Feb 2016**: (KBW) Added try/except blocks for importing
        ConfigParser.
"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import sys
if sys.version > '3':
    long = int
    try:
        from configparser import ConfigParser
    except ImportError:
        print('WARNING: Unable to import configparser!  Beware!')
    try:
        from configparser import ExtendedInterpolation
    except ImportError:
        print('WARNING: Unable to import ExtendedInterpolation!  Some configurations will fail!')
else:
    try:
        from ConfigParser import ConfigParser
    except ImportError:
        print('WARNING: Unable to import ConfigParser!  Beware!')
    try:
        from ConfigParser import ExtendedInterpolation
    except ImportError:
        print('WARNING: Unable to import ExtendedInterpolation!  Some configurations will fail!')

import os.path
from os import environ
import glob
import numpy
from mangadap.par.database_definitions import TemplateLibraryDef
from mangadap.config.util import validate_template_config
from mangadap.util.exception_tools import check_environment_variable

__author__ = 'Kyle B. Westfall'

def default_drp_version():
    """
    Return the DRP version defined by the environmental variable
    MANGADRP_VER.

    """
    check_environment_variable('MANGADRP_VER')
    return environ['MANGADRP_VER']


def default_redux_path(drpver=None):
    """
    Return the main output path for the DRP products using the
    environmental variable MANGA_SPECTRO_REDUX.

    Args:
        drpver (str): DRP version

    Returns:
        str: Path to reduction directory
    
    """
    # Make sure the DRP version is set
    if drpver is None:
        drpver = default_drp_version()
    check_environment_variable('MANGA_SPECTRO_REDUX')
    return os.path.join(environ['MANGA_SPECTRO_REDUX'], drpver)


def default_drp_directory_path(drpver, redux_path, plate):
    """
    Return the exact directory path with the DRP file.

    Args:
        drpver (str): DRP version
        redux_path (str): Path to the root reduction directory
        plate (int): Plate number

    Returns:
        str: Path to the directory with the 3D products of the DRP

    """
    # Make sure the redux path is set
    if redux_path is None:
        redux_path = default_redux_path(drpver)

    return os.path.join(redux_path, str(plate), 'stack')


# TODO: Are these values kept in MANGACORE somewhere?
def default_cube_pixelscale():
    """
    Return the default pixel scale of the DRP CUBE files in arcsec.
    """
    return 0.5


def default_cube_width_buffer():
    """
    Return the default width buffer in pixels used when regridding the
    DRP RSS spectra into the CUBE format.
    """
    return 10


def default_cube_recenter():
    """
    Return the default recentering flag used when regridding the DRP RSS
    spectra into the CUBE format.
    """
    return False


def default_regrid_rlim():
    """
    Return the default limiting radius for the Gaussian kernel used when
    regridding the DRP RSS spectra into the CUBE format.
    """
    return 1.6


def default_regrid_sigma():
    """
    Return the default standard deviation of the Gaussian kernel used
    when regridding the DRP RSS spectra into the CUBE format.
    """
    return 0.7


def default_cube_covariance_file(plate, ifudesign):
    """
    Return the default name for the covariance calculated for a 'CUBE'
    fits file.

    Args:
        plate (int): Plate number
        ifudesign (int): IFU design number
    
    Returns:
        str : The default name of the covariance file.
    
    """
    root = default_manga_fits_root(plate, ifudesign, 'CUBE')
    return ('{0}_COVAR.fits'.format(root))


def default_dap_source():
    """
    Return the root path to the DAP source directory using the
    environmental variable MANGADAP_DIR.

    """
    check_environment_variable('MANGADAP_DIR')
    return environ['MANGADAP_DIR']


def default_dap_version():
    """
    Return the DAP version defined by the environmental variable
    MANGADAP_VER.

    """
    check_environment_variable('MANGADAP_VER')
    return environ['MANGADAP_VER']


def default_analysis_path(drpver=None, dapver=None):
    """
    Return the main output path for the DAP using the environmental
    variable MANGA_SPECTRO_ANALYSIS.

    Args:
        drpver (str): DRP version
        dapver (str): DAP version

    Returns:
        str: Path to analysis directory

    """
    # Make sure the DRP version is set
    if drpver is None:
        drpver = default_drp_version()
    # Make sure the DAP version is set
    if dapver is None:
        dapver = default_dap_version()
    check_environment_variable('MANGA_SPECTRO_ANALYSIS')
    return os.path.join(environ['MANGA_SPECTRO_ANALYSIS'], drpver, dapver)


def default_dap_directory_path(drpver, dapver, analysis_path, plate, ifudesign):
    """
    Return the exact directory path with the DAP file.

    Args:
        drpver (str): DRP version
        dapver (str): DAP version
        analysis_path (str): Path to the root analysis directory
        plate (int): Plate number
        ifudesign (int): IFU design number

    Returns:
        str: Path to the directory with DAP output files

    """
    # Make sure the DRP version is set
    if analysis_path is None:
        analysis_path = default_analysis_path(drpver, dapver)

    return os.path.join(analysis_path, str(plate), str(ifudesign))


def default_manga_fits_root(plate, ifudesign, mode):
    """
    Generate the main root name for the output MaNGA fits files for a
    given plate/ifudesign/mode.

    .. todo::
        - Include a "sampling mode" (LIN vs. LOG)?

    Args:
        plate (int): Plate number
        ifudesign (int): IFU design number
        mode (str): Mode of the DRP reduction; either RSS or CUBE

    Returns:
        str : Root name for a MaNGA fits file:
            manga-[PLATE]-[IFUDESIGN]-LOG[MODE]
    """
    return 'manga-{0}-{1}-LOG{2}'.format(plate, ifudesign, mode)


def default_dap_file_root(plate, ifudesign, mode):
    """
    Generate the root name of the MaNGA DAP parameter and script files
    for a given plate/ifudesign/mode.
    
    Args:
        plate (int): Plate number
        ifudesign (int): IFU design number
        mode (str): Mode of the DRP reduction; either RSS or CUBE

    Returns:
        str : Root name for the DAP file:
            mangadap-[PLATE]-[IFUDESIGN]-LOG[MODE]

    """
    return 'mangadap-{0}-{1}-LOG{2}'.format(plate, ifudesign, mode)


def default_dap_par_file(drpver, dapver, analysis_path, directory_path, plate, ifudesign, mode):
    """
    Return the full path to the DAP par file.

    Args:
        drpver (str): DRP version
        dapver (str): DAP version
        analysis_path (str): Path to the root analysis directory
        directory_path (str): Path to the directory with the DAP output
            files
        plate (int): Plate number
        ifudesign (int): IFU design number
        mode (str): Mode of the DRP reduction; either RSS or CUBE

    Returns:
        str: Full path to the DAP par file

    """
    # Make sure the directory path is defined
    if directory_path is None:
        directory_path = default_dap_directory_path(drpver, dapver, analysis_path, plate, ifudesign)
    # Set the name of the par file; put this in its own function?
    par_file = '{0}.par'.format(default_dap_file_root(plate, ifudesign, mode))
#    par_file = 'mangadap-{0}-{1}-LOG{2}.par'.format(plate, ifudesign, mode)
    return os.path.join(directory_path, par_file)

    
def default_dap_plan_file(drpver, dapver, analysis_path, directory_path, plate, ifudesign, mode):
    """
    Return the full path to the DAP plan file.

    Args:
        drpver (str): DRP version
        dapver (str): DAP version
        analysis_path (str): Path to the root analysis directory
        directory_path (str): Path to the directory with the DAP output
            files
        plate (int): Plate number
        ifudesign (int): IFU design number
        mode (str): Mode of the DRP reduction; either RSS or CUBE

    Returns:
        str: Full path to the DAP plan file

    """
    # Make sure the directory path is defined
    if directory_path is None:
        directory_path = default_dap_directory_path(drpver, dapver, analysis_path, plate, ifudesign)
    # Set the name of the plan file; put this in its own function?
#    plan_file = 'manga-{0}-{1}-LOG{2}-dapplan.par'.format(plate, ifudesign, mode)
    plan_file = '{0}-plan.par'.format(default_dap_file_root(plate, ifudesign, mode))
    return os.path.join(directory_path, plan_file)


def default_dap_file_name(plate, ifudesign, mode, bintype, niter):
    """
    Return the name of the DAP output fits file.

    Args:
        plate (int): Plate number
        ifudesign (int): IFU design number
        mode (str): Mode of the DRP reduction; either RSS or CUBE
        bintype (str): Binning type used by the DAP
        niter (int): Iteration number

    Returns:
        str: Name of the DAP output file.

    """
    # Number of spaces provided for iteration number is 3
    siter = str(niter).zfill(3)
    root = default_manga_fits_root(plate, ifudesign, mode)
    # TODO: ? Use default_dap_file_root() or change to:
#   return 'manga-{0}-{1}-LOG{2}_{3}-{4}.fits'.format(plate, ifudesign, mode, bintype, siter)
    return '{0}_BIN-{1}-{2}.fits'.format(root, bintype, siter)


def default_template_library_file(plate, ifudesign, mode, library_key, spindex_key=None):
    """
    Return the name of the processed template file.

    Args:
        plate (int): Plate number
        ifudesign (int): IFU design number
        mode (str): Mode of the DRP reduction; either RSS or CUBE
        library_key (str): Template library keyword
        spindex_key (str): (Optional) Spectral index library keyword

    Returns:
        str: Name of the fits file containing the template library as
            prepared for use by the DAP.

    .. todo:
        - Add a docstring link to the default template libraries

    """
    # TODO: Use default_dap_file_root()?
    root = default_manga_fits_root(plate, ifudesign, mode)
    if spindex_key is not None:
        return '{0}_{1}_{2}.fits'.format(root, library_key, spindex_key)
    return '{0}_{1}.fits'.format(root, library_key)


def default_template_libraries(dapsrc=None):
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
    |       M11-MARCS |       2.73 |      No |        full |  None |
    +-----------------+------------+---------+-------------+-------+
    |      M11-STELIB |       3.40 |      No |        full |  None |
    +-----------------+------------+---------+-------------+-------+
    | M11-STELIB-ZSOL |       3.40 |      No |        full |  None |
    +-----------------+------------+---------+-------------+-------+
    |      M11-ELODIE |       0.55 |      No |      < 6795 |  None |
    +-----------------+------------+---------+-------------+-------+
    |       M11-MILES |       2.54 |      No | 3550 - 7400 |  None |
    +-----------------+------------+---------+-------------+-------+
    |           MILES |       2.50 |      No |      < 7400 |  None |
    +-----------------+------------+---------+-------------+-------+
    |       MILES-AVG |       2.50 |      No |      < 7400 |  None |
    +-----------------+------------+---------+-------------+-------+
    |      MILES-THIN |       2.50 |      No |      < 7400 |  None |
    +-----------------+------------+---------+-------------+-------+
    |          STELIB |       3.40 |      No |        full |   0.0 |
    +-----------------+------------+---------+-------------+-------+
    |         MIUSCAT |       2.51 |      No | 3480 - 9430 |  None |
    +-----------------+------------+---------+-------------+-------+
    |    MIUSCAT-THIN |       2.51 |      No | 3480 - 9430 |  None |
    +-----------------+------------+---------+-------------+-------+

    .. warning::

        Function is currently only valid for Python 3.2 or greater!

    Args:
        dapsrc (str): (Optional) Root path to the DAP source
            directory.  If not provided, the default is defined by
            :func:`mangadap.util.defaults.default_dap_source`.

    Returns:
        list : An list with of :class:`mangadap.par.ParSet.ParSet`
            objects, each of which is created using
            :func:`mangadap.par.database_definitions.TemplateLibraryDef`.

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
    """
    # Check the source directory exists
    dapsrc = default_dap_source() if dapsrc is None else str(dapsrc)
    if not os.path.isdir(dapsrc):
        raise NotADirectoryError('{0} does not exist!'.format(dapsrc))

    # Check the configuration files exist
    ini_files = glob.glob(dapsrc+'/python/mangadap/config/templates/*.ini')
    if len(ini_files) == 0:
        raise IOError('Could not find any template library configuration files in {0} !'.format(
                      dapsrc+'/python/mangadap/config/templates'))

    # Build the list of library definitions
    template_libraries = []
    for f in ini_files:
        # Read the config file
        cnfg = ConfigParser(environ, allow_no_value=True, interpolation=ExtendedInterpolation())
        cnfg.read(f)
        # Ensure it has the necessary elements to define the template
        # library
        validate_template_config(cnfg)
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


def default_plate_target_files():
    """
    Return the default plateTarget files in mangacore and their
    associated catalog indices.  The catalog indices are determined
    assuming the file names are of the form::

        'plateTargets-{0}.par'.format(catalog_id)

    Returns:
        numpy.array : Two arrays: the first contains the identified
            plateTargets files found using the default search
            string, the second provides the integer catalog index
            determined for each file.
    """
    # Default search string
    check_environment_variable('MANGACORE_DIR')
    search_str = os.path.join(environ['MANGACORE_DIR'], 'platedesign', 'platetargets',
                              'plateTargets*.par')
    file_list = glob.glob(search_str)                       # List of files
    nfiles = len(file_list)
    trgid = numpy.zeros(nfiles, dtype=numpy.int)            # Array to hold indices
    for i in range(nfiles):
        suffix = file_list[i].split('-')[1]                 # Strip out the '{i}.par'
        trgid[i] = int(suffix[:suffix.find('.')])           # Strip out '{i}'
    
    return numpy.array(file_list), trgid


