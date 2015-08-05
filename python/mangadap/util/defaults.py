"""

Provides a set of functions that define and return defaults used by the
MaNGA DAP, such as paths and file names.

*Source location*:
    $MANGADAP_DIR/python/mangadap/util/defaults.py

*Python2/3 compliance*::

    from __future__ import division
    from __future__ import print_function
    from __future__ import absolute_import
    
    import sys
    if sys.version > '3':
        long = int

*Imports*::

    import os.path
    from os import environ
    import numpy
    from mangadap.util.par import TemplateLibraryParSet

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

"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import sys
if sys.version > '3':
    long = int

import os.path
from os import environ
import numpy
from mangadap.util.par import TemplateLibraryParSet

__author__ = 'Kyle B. Westfall'

def default_drp_version():
    """
    Return the DRP version defined by the environmental variable
    MANGADRP_VER.

    """
    return environ['MANGADRP_VER']


def default_redux_path(drpver):
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


def default_dap_source():
    """
    Return the root path to the DAP source directory using the
    environmental variable MANGADAP_DIR.

    """
    return environ['MANGADAP_DIR']


def default_dap_version():
    """
    Return the DAP version defined by the environmental variable
    MANGADAP_VER.

    """
    return environ['MANGADAP_VER']


def default_analysis_path(drpver, dapver):
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
    par_file = 'mangadap-{0}-{1}-LOG{2}.par'.format(plate, ifudesign, mode)
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
    # TODO: Change to:
#   plan_file = 'mangadap-{0}-{1}-LOG{2}-plan.par'.format(plate, ifudesign, mode)
    plan_file = 'manga-{0}-{1}-LOG{2}-dapplan.par'.format(plate, ifudesign, mode)
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
    # TODO: Change to:
#   return 'manga-{0}-{1}-LOG{2}_{3}-{4}.fits'.format(plate, ifudesign, mode, bintype, siter)
    return 'manga-{0}-{1}-LOG{2}_BIN-{3}-{4}.fits'.format(plate, ifudesign, mode, bintype, siter)


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
    if spindex_key is not None:
        return 'manga-{0}-{1}-LOG{2}_{3}_{4}.fits'.format(plate, ifudesign, mode, library_key,
                                                          spindex_key)
    return 'manga-{0}-{1}-LOG{2}_{3}.fits'.format(plate, ifudesign, mode, library_key)


def default_template_libraries(dapsrc=None):
    """
    Return the list of library keys, the searchable string of the 1D
    template library fits files for the template libraries available for
    use by the DAP, the FWHM of the libraries, and whether or not the
    wavelengths are in vacuum.

    The stellar template library files should be a list of 1D fits
    files, and be associated with one of the following library keys:

    +--------------+------------+---------+-------------+-------+
    |              |   Spectral |         |  Wavelength | Lower |
    |          KEY |  res (ang) | Vacuum? | Range (ang) | Limit |
    +==============+============+=========+=============+=======+
    |    M11-MARCS |       2.73 |      No |        full |  None |
    +--------------+------------+---------+-------------+-------+
    |   M11-STELIB |       3.40 |      No |        full |  None |
    +--------------+------------+---------+-------------+-------+
    |   M11-ELODIE |       0.55 |      No |      < 6795 |  None |
    +--------------+------------+---------+-------------+-------+
    |    M11-MILES |       2.54 |      No | 3550 - 7400 |  None |
    +--------------+------------+---------+-------------+-------+
    |        MILES |       2.50 |      No |      < 7400 |  None |
    +--------------+------------+---------+-------------+-------+
    |    MILES-AVG |       2.50 |      No |      < 7400 |  None |
    +--------------+------------+---------+-------------+-------+
    |   MILES-THIN |       2.50 |      No |      < 7400 |  None |
    +--------------+------------+---------+-------------+-------+
    |       STELIB |       3.40 |      No |        full |   0.0 |
    +--------------+------------+---------+-------------+-------+
    |      MIUSCAT |       2.51 |      No | 3480 - 9430 |  None |
    +--------------+------------+---------+-------------+-------+
    | MIUSCAT-THIN |       2.51 |      No | 3480 - 9430 |  None |
    +--------------+------------+---------+-------------+-------+

    Args:
        dapsrc (str): (Optional) Root path to the DAP source
            directory.  If not provided, the default is defined by
            :func:`mangadap.util.defaults.default_dap_source`.

    Returns:
        list : An list with of :class:`mangadap.util.par.ParSet`
            objects, each of which is created using
            :func:`mangadap.util.par.TemplateLibraryParSet`.

    Raises:
        NotADirectoryError: Raised if the provided or default
            *dapsrc* is not a directory.

    """
    dapsrc = default_dap_source() if dapsrc is None else str(dapsrc)
    if not os.path.isdir(dapsrc):
        raise NotADirectoryError('{0} does not exist!'.format(dapsrc))

    template_libraries = []

    # M11-MARCS
    template_libraries += \
        [ TemplateLibraryParSet(key='M11-MARCS',
                                file_search=dapsrc+'/external/templates/m11_marcs/*_s.fits',
                                # TODO: This is the resolution in the header of the files, is it
                                # right?
                                fwhm=2.73, in_vacuum=False,
                                wave_limit=numpy.array([ None, None ]),
                                lower_flux_limit=None)
        ]

    # M11-STELIB
    template_libraries += \
        [ TemplateLibraryParSet(key='M11-STELIB',
                                file_search=dapsrc+'/external/templates/m11_stelib/*_s.fits',
                                # Changed in IDL version on 27 Jan 2015 to match Maraston &
                                # Strombach (2011) Table 1
                                fwhm=3.40, in_vacuum=False,
                                wave_limit=numpy.array([ None, None ]),
                                lower_flux_limit=None)
        ]

    # M11-ELODIE
    template_libraries += \
        [ TemplateLibraryParSet(key='M11-ELODIE',
                                file_search=dapsrc+'/external/templates/m11_elodie/*.fits',
                                # Resolution taken from Maraston & Strombach (2011)
                                fwhm=0.55, in_vacuum=False,
                                wave_limit=numpy.array([ None, 6795. ]),
                                lower_flux_limit=None)
        ]

    # M11-MILES
    template_libraries += \
        [ TemplateLibraryParSet(key='M11-MILES',
                                file_search=dapsrc+'/external/templates/m11_miles/*.fits',
                                # TODO: Should this be the same as MILES?
                                fwhm=2.54, in_vacuum=False,
                                wave_limit=numpy.array([ 3550., 7400. ]),
                                lower_flux_limit=None)
        ]

    # MILES
    template_libraries += \
        [ TemplateLibraryParSet(key='MILES',
                                file_search=dapsrc+'/external/templates/miles/*.fits',
                                # TODO: Unknown if this library is in vacuum or in air
                                # Resolution from Falcon-Barroso et al. (2011, A&A, 532, 95)
                                fwhm=2.50, in_vacuum=False,
                                wave_limit=numpy.array([ 3575., 7400. ]),
                                lower_flux_limit=0.0)
        ]

    # MILES-AVG
    template_libraries += \
        [ TemplateLibraryParSet(key='MILES-AVG',
                                file_search=dapsrc+'/external/templates/miles_avg/*.fits',
                                # TODO: Unknown if this library is in vacuum or in air
                                # Resolution from Falcon-Barroso et al. (2011, A&A, 532, 95)
                                fwhm=2.50, in_vacuum=False,
                                wave_limit=numpy.array([ 3575., 7400. ]),
                                lower_flux_limit=None)
        ]

    # MILES-THIN
    template_libraries += \
        [ TemplateLibraryParSet(key='MILES-THIN',
                                file_search=dapsrc+'/external/templates/miles_thin/*.fits',
                                # TODO: Unknown if this library is in vacuum or in air
                                # Resolution from Falcon-Barroso et al. (2011, A&A, 532, 95)
                                fwhm=2.50, in_vacuum=False,
                                wave_limit=numpy.array([ 3575., 7400. ]),
                                lower_flux_limit=None)
        ]

    # STELIB
    template_libraries += \
        [ TemplateLibraryParSet(key='STELIB',
                                file_search=dapsrc+'/external/templates/stelib/*.fits',
                                # TODO: Unknown if this library is in vacuum or in air
                                # Only given as approximate in Le Borgne et al. (2003, A&A, 402,
                                # 433); assume value from Maraston & Strombach (2011, MNRAS,
                                # 418, 2785)
                                fwhm=3.40, in_vacuum=False,
                                wave_limit=numpy.array([ None, None ]),
                                lower_flux_limit=0.0)
        ]

    # MIUSCAT
    template_libraries += \
        [ TemplateLibraryParSet(key='MIUSCAT',
                                file_search=dapsrc+'/external/templates/miuscat/*.fits',
                                # TODO: Unknown if this library is in vacuum or in air
                                # Using value provided by the header in all the fits files
                                fwhm=2.51, in_vacuum=False,
                                wave_limit=numpy.array([ 3480., 9430. ]),
                                lower_flux_limit=None)
        ]

    # MIUSCAT-THIN
    template_libraries += \
        [ TemplateLibraryParSet(key='MIUSCAT-THIN',
                                file_search=dapsrc+'/external/templates/miuscat_thin/*.fits',
                                # TODO: Unknown if this library is in vacuum or in air
                                # Using value provided by the header in all the fits files
                                fwhm=2.51, in_vacuum=False,
                                wave_limit=numpy.array([ 3480., 9430. ]),
                                lower_flux_limit=None)
        ]

    return template_libraries


