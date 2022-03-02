# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
r"""
Provides a set of functions that define and return defaults used by the
MaNGA DAP, such as paths and file names.

----        

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from pathlib import Path
import os
from pkg_resources import resource_filename

import numpy

from . import manga_environ


# TODO: Move this to util/pkg.py?
def dap_source_dir():
    from pkg_resources import resource_filename
    return resource_filename('mangadap', '')


def dap_data_root():
    """Return the root directory with the DAP data."""
    return Path(resource_filename('mangadap', 'data')).resolve()


def dap_config_root():
    """Return the root directory with the DAP config data."""
    return Path(resource_filename('mangadap', 'config')).resolve()


def sdss_maskbits_file():
    """Return the path to the sdss maskbits yanny file."""
    maskbits_file = dap_data_root() / 'sdss' / 'sdssMaskbits.par'
    return maskbits_file if maskbits_file.exists() else None


#def drp_version():
#    """
#    Return the DRP version defined by the environmental variable
#    MANGADRP_VER.
#    """
#    return manga_environ['MANGADRP_VER']
#
#
#def drp_redux_path(drpver=None):
#    """
#    Return the main output path for the DRP products using the
#    environmental variable MANGA_SPECTRO_REDUX.
#
#    Args:
#        drpver (:obj:`str`, optional):
#            DRP version. Default is to use :func:`drp_version`.
#
#    Returns:
#        :obj:`str`: Path to reduction directory
#    """
#    # Make sure the DRP version is set
#    if drpver is None:
#        drpver = drp_version()
#    return Path(manga_environ['MANGA_SPECTRO_REDUX']).resolve() / drpver
#
#
#def drp_directory_path(plate, drpver=None, redux_path=None):
#    """
#    Return the exact directory path with the DRP file.
#
#    Args:
#        plate (:obj:`int`):
#            Plate number
#        drpver (:obj:`str`, optional):
#            DRP version. Default is to use :func:`drp_version`.
#        redux_path (:obj:`str`, optional):
#            Path to the root reduction directory. Default is to use
#            :func:`drp_redux_path`.
#
#    Returns:
#        :obj:`str`: Path to the directory with the 3D products of the
#        DRP
#    """
#    # Make sure the redux path is set
#    _redux_path = drp_redux_path(drpver=drpver) \
#                        if redux_path is None else Path(redux_path).resolve()
#    return _redux_path / str(plate) / 'stack'
#
#
#def drp_finding_chart_path(plate, drpver=None, redux_path=None):
#    """
#    Return the exact directory path with the finding charts for a given plate.
#
#    Args:
#        plate (:obj:`int`):
#            Plate number
#        drpver (:obj:`str`, optional):
#            DRP version. Default is to use :func:`drp_version`.
#        redux_path (:obj:`str`, optional):
#            Path to the root reduction directory. Default is to use
#            :func:`drp_redux_path`.
#
#    Returns:
#        :obj:`str`: Path to the directory with the finding charts.
#    """
#    # Make sure the redux path is set
#    _redux_path = drp_redux_path(drpver=drpver) \
#                        if redux_path is None else Path(redux_path).resolve()
#    return _redux_path / str(plate) / 'images'
#
#
#def drpall_file(drpver=None, redux_path=None):
#    """
#    Return the path to the DRPall file.
#
#    Args:
#        drpver (:obj:`str`, optional):
#            DRP version.  Default is to use :func:`drp_version`.
#        redux_path (:obj:`str`, optional):
#            Path to the root reduction directory. Default is to use
#            :func:`drp_redux_path`.
#
#    Returns:
#        :obj:`str`: Full path to the DRPall fits file.
#    """
#    _drpver = drp_version() if drpver is None else drpver
#    _redux_path = drp_redux_path(drpver=_drpver) \
#                        if redux_path is None else Path(redux_path).resolve()
#    return _redux_path / f'drpall-{_drpver}.fits'


def dapall_file(drpver=None, dapver=None, analysis_path=None):
    """
    Return the path to the DAPall file.

    Args:
        drpver (:obj:`str`, optional):
            DRP version.  Default is to use :func:`drp_version`.
        dapver (:obj:`str`, optional):
            DAP version.  Default is to use :func:`dap_version`.
        analysis_path (:obj:`str`, optional):
            Path to the root analysis directory.  Default is to use
            :func:`dap_analysis_path`

    Returns:
        :obj:`str`: Full path to the DAPall fits file.
    """
    _drpver = drp_version() if drpver is None else drpver
    _dapver = dap_version() if dapver is None else dapver
    _analysis_path = dap_analysis_path(drpver=_drpver, dapver=_dapver) \
                        if analysis_path is None else Path(analysis_path).resolve()
    return _analysis_path / f'dapall-{_drpver}-{_dapver}.fits'


def cube_pixelscale():
    """
    Return the default pixel scale of the DRP CUBE files in arcsec.
    """
    return 0.5


def cube_width_buffer():
    """
    Return the default width buffer in pixels used when regridding the
    DRP RSS spectra into the CUBE format.
    """
    return 10


def cube_recenter():
    """
    Return the default recentering flag used when regridding the DRP RSS
    spectra into the CUBE format.
    """
    return False


def regrid_rlim():
    """
    Return the default limiting radius for the Gaussian kernel used when
    regridding the DRP RSS spectra into the CUBE format.
    """
    return 1.6


def regrid_sigma():
    """
    Return the default standard deviation of the Gaussian kernel used
    when regridding the DRP RSS spectra into the CUBE format.
    """
    return 0.7


def dap_version():
    """
    Return the DAP version defined by the environmental variable
    MANGADAP_VER.  If that environmental variable does not exist,
    `mangadap.__version__` is returned.
    """
    return manga_environ['MANGADAP_VER']


def dap_analysis_path(drpver=None, dapver=None):
    """
    Return the main output path for the DAP using the environmental
    variable MANGA_SPECTRO_ANALYSIS.

    Args:
        drpver (:obj:`str`, optional):
            DRP version. Default is to use :func:`drp_version`.
        dapver (:obj:`str`, optional):
            DAP version. Default is to use :func:`dap_version`.

    Returns:
        :obj:`str`: Path to analysis directory
    """
    # Make sure the DRP version is set
    if drpver is None:
        drpver = drp_version()
    # Make sure the DAP version is set
    if dapver is None:
        dapver = dap_version()
    return Path(manga_environ['MANGA_SPECTRO_ANALYSIS']).resolve() / drpver / dapver








def dap_common_path(plate=None, ifudesign=None, drpver=None, dapver=None, analysis_path=None):
    """
    Return the path to the path to the directory with data common to
    multiple binning schemes.

    Args:
        plate (:obj:`int`, optional):
            Plate number, for reference directory of a specific
            plate.
        ifudesign (:obj:`int`, optional):
            IFU design number.
        drpver (:obj:`str`, optional):
            DRP version. Default is to use :func:`drp_version`.
        dapver (:obj:`str`, optional):
            DAP version.  Default is to use :func:`dap_version`.
        analysis_path (:obj:`str`, optional):
            Path to the root analysis directory. Default is to use
            :func:`dap_analysis_path`.

    Returns:
        :obj:`str`: Path to the directory with DAP reference files

    Raises:
        ValueError:
            Raised if IFU design is provided and plate is not.
    """
    # For hierarchy, if ifudesign is given, plate must also be given
    if plate is None and ifudesign is not None:
        raise ValueError('For IFU design subdirectory, must provide plate number.')

    # Get the main analysis path
    _analysis_path = dap_analysis_path(drpver=drpver, dapver=dapver) \
                        if analysis_path is None else Path(analysis_path).resolve()

    output_path = _analysis_path / 'common'
    if plate is None:
        return output_path
    output_path = output_path / str(plate)
    return output_path if ifudesign is None else output_path / str(ifudesign)


#def dap_method(binning_method, stellar_continuum_templates, emission_line_model_templates):
#    """
#    Construct the ``DAPTYPE`` based on the analysis plan.
#
#    .. include:: ../include/daptype.rst
#
#    Args:
#        binning_method (:obj:`str`):
#            String used to define the binning method.
#        stellar_continuum_templates (:obj:`str`):
#            String defining the template library used in the stellar
#            continuum (stellar kinematics) modeling.
#        emission_line_model_templates (:obj:`str`):
#            String defining the template library used in the
#            emission-line modeling.  Can be None, meaning that the
#            continuum templates are the same as used for the
#            stellar-continuum modeling.
#
#    Returns:
#        :obj:`str`: The string representation of the analysis method.
#
#    """
#    _eltpl = stellar_continuum_templates if emission_line_model_templates is None \
#                                                else emission_line_model_templates
#    return f'{binning_method}-{stellar_continuum_templates}-{_eltpl}'


def dap_method_path(method, plate=None, ifudesign=None, qa=False, ref=False, drpver=None,
                    dapver=None, analysis_path=None):
    """
    Return the path to the designated subdirectory built using the plan
    key identifiers or directly using the provided method.

    The "method" identifying each plan is currently build using the
    keywords for the binning type and the continuum-fitting key.
    Nominally, the latter should include the template set used.

    Args:
        method (:obj:`str`):
            String defining the method identifier for a set of DAP
            output files. These should be built using
            :func:`dap_method`.
        plate (:obj:`int`, optional):
            Plate number.
        ifudesign (:obj:`int`, optional):
            IFU design number.
        qa (:obj:`bool`, optional):
            Give the path to the qa/ subdirectory
        ref (:obj:`bool`, optional):
            Give the path to the ``ref/`` subdirectory.
        drpver (:obj:`str`, optional):
            DRP version. Default is to use :func:`drp_version`.
        dapver (:obj:`str`, optional):
            DAP version. Default is to use :func:`dap_version`.
        analysis_path (:obj:`str`, optional):
            Path to the root analysis directory. Default is to use
            :func:`dap_analysis_path`.

    Returns:
        :obj:`str`: Path to the plan subdirectory

    Raises:
        ValueError:
            Raised if IFU design is provided and plate is not, or if
            either qa or ref are true and one or both of plate and
            IFU design are not provided..
    """
    # For heirarchy, if ifudesign is given, plate must also be given
    if (plate is not None and ifudesign is None) or (plate is None and ifudesign is not None):
        raise ValueError('Must provide both plate and ifudesign, if providing either.')
    if qa and ref:
        raise ValueError('Cannot provide path for both qa and ref directory.  Pick one.')

    # Get the main analysis path
    _analysis_path = dap_analysis_path(drpver=drpver, dapver=dapver) \
                        if analysis_path is None else Path(analysis_path).resolve()

    # Build the plan subirectory
    output_path = _analysis_path / method
    if plate is None:
        return output_path
    output_path = output_path / str(plate) / str(ifudesign)
    if not qa and not ref:
        return output_path
    return output_path / ('qa' if qa else 'ref')


def manga_fits_root(plate, ifudesign, mode=None):
    """
    Generate the main root name for the output MaNGA fits files for a
    given plate/ifudesign/mode.

    Args:
        plate (:obj:`int`):
            Plate number
        ifudesign (:obj:`int`):
            IFU design number
        mode (:obj:`str`, optional):
            Mode of the output fits file. Options are: ``'LINCUBE'``,
            ``'LINRSS'``, ``'LOGCUBE'``, ``'LOGRSS'``, or ``'MAPS'``.
            Default is that no mode is included in the name.

    Returns:
        :obj:`str`: Root name for a MaNGA fits file:
        ``manga-[PLATE]-[IFUDESIGN]-[MODE]``

    Raises:
        ValueError:
            Raised if mode is not a valid option.
    """
    if mode not in [ None, 'LINCUBE', 'LINRSS', 'LOGCUBE', 'LOGRSS', 'MAPS' ]:
        raise ValueError('Do not understand mode={0}.'.format(mode))
    return f'manga-{plate}-{ifudesign}' if mode is None else f'manga-{plate}-{ifudesign}-{mode}'


def dap_file_root(plate, ifudesign, mode=None):
    """
    Generate the root name of the MaNGA DAP parameter and script
    files for a given plate/ifudesign/mode.
    
    Args:
        plate (:obj:`int`):
            Plate number
        ifudesign (:obj:`int`):
            IFU design number
        mode (:obj:`str`, optional):
            Mode of the DRP reduction; either ``RSS`` or ``CUBE``. If
            None, the mode is excluded from the file root.

    Returns:
        :obj:`str`: Root name for the DAP file:
        ``mangadap-[PLATE]-[IFUDESIGN]`` or
        ``mangadap-[PLATE]-[IFUDESIGN]-LOG[MODE]``
    """
    return f'mangadap-{plate}-{ifudesign}' if mode is None else \
                f'mangadap-{plate}-{ifudesign}-LOG{mode}'

    
def dap_config(plate, ifudesign, drpver=None, dapver=None, analysis_path=None,
               directory_path=None):
    """
    Return the full path to the DAP configuration file.

    The configuration file provides the input data necessary to
    instantiate a :class:`mangadap.datacube.manga.MaNGADataCube`.
    
    Args:
        plate (:obj:`int`):
            Plate number
        ifudesign (:obj:`int`):
            IFU design number
        drpver (:obj:`str`, optional):
            DRP version. Default is to use :func:`drp_version`.
        dapver (:obj:`str`, optional):
            DAP version. Default is to use :func:`dap_version`.
        analysis_path (:obj:`str`, optional): 
            Path to the root analysis directory. Default is to use
            :func:`dap_analysis_path`.
        directory_path (:obj:`str`, optional):
            Path to the directory with the DAP output files. Default
            is to use :func:`dap_common_path`

    Returns:
        :obj:`str`: Full path to the DAP par file
    """
    # Make sure the directory path is defined
    _directory_path = dap_common_path(plate=plate, ifudesign=ifudesign, drpver=drpver,
                                      dapver=dapver, analysis_path=analysis_path) \
                            if directory_path is None else Path(directory_path).resolve()
    # Set the name of the par file; put this in its own function?
    return _directory_path /  f'{dap_file_root(plate, ifudesign, "CUBE")}.ini'

    
def dap_plan_file(drpver=None, dapver=None, analysis_path=None):
    """
    Return the full path to the DAP plan file.

    Args:
        drpver (:obj:`str`, optional):
            DRP version. Default is to use :func:`drp_version`.
        dapver (:obj:`str`, optional):
            DAP version. Default is to use :func:`dap_version`.
        analysis_path (:obj:`str`, optional):
            Path to the root analysis directory. Default is to use
            :func:`dap_analysis_path`

    Returns:
        :obj:`str`: Full path to the DAP plan file
    """
    # Get the main analysis path
    if dapver is None:
        dapver = dap_version()

    _analysis_path = dap_analysis_path(drpver=drpver, dapver=dapver) \
                        if analysis_path is None else Path(analysis_path).resolve()
    
    # Set the name of the plan file
    return _analysis_path / f'mangadap-plan-{dapver}.par'


def dap_file_name(plate, ifudesign, output_mode, mode=None, compressed=True):
    """
    Return the name of the DAP output fits file.

    Args:
        plate (:obj:`int`):
            Plate number
        ifudesign (:obj:`int`):
            IFU design number
        output_mode (:obj:`str`):
            Output mode designation
        mode (:obj:`str`, optional):
            Mode of the output fits file. Options are: ``'LINCUBE'``,
            ``'LINRSS'``, ``'LOGCUBE'``, ``'LOGRSS'``, or ``'MAPS'``.
            Default is that no mode is included in the name.
        compressed (:obj:`bool`, optional):
            Append '.gz' to the output file name.

    Returns:
        :obj:`str`: Name of the DAP output file.
    """
    # Number of spaces provided for iteration number is 3
    root = manga_fits_root(plate, ifudesign, mode=mode)
    return ('{0}-{1}.fits.gz' if compressed else '{0}-{1}.fits').format(root, output_mode)


#def plate_target_files():
#    """
#    Return the default plateTarget files in mangacore and their
#    associated catalog indices. The catalog indices are determined
#    assuming the file names are of the form::
#
#        'plateTargets-{0}.par'.format(catalog_id)
#
#    Returns:
#        `numpy.ndarray`_: Two arrays: the first contains the
#        identified plateTargets files found using the default search
#        string, the second provides the integer catalog index
#        determined for each file.
#    """
#    # Default search string
#    search_str = Path(manga_environ['MANGACORE_DIR']).resolve() / 'platedesign' / 'platetargets'
#    file_list = sorted(list(search_path.glob('plateTargets*.par')))
#    nfiles = len(file_list)
#    if nfiles == 0:
#        raise ValueError('Unable to find any plateTargets files!')
#    trgid = numpy.zeros(nfiles, dtype=int)                  # Array to hold indices
#    for i in range(nfiles):
#        suffix = file_list[i].split('-')[1]                 # Strip out the '{i}.par'
#        trgid[i] = int(suffix[:suffix.find('.')])           # Strip out '{i}'
#    
#    return numpy.array(file_list), trgid
#
#
#def redshift_fix_file():
#    """
#    Return the path to the default redshift fix file.
#
#    Returns:
#        :obj:`str`: Expected path to the redshift-fix parameter file.
#    """
#    return dap_data_root() / 'fix' / 'redshift_fix.par'
#
#
#def photometry_fix_file():
#    """
#    Return the path to the default photometry fix file.
#
#    Returns:
#        :obj:`str`: Expected path to the photometry-fix parameter file.
#    """
#    return dap_data_root() / 'fix' / 'photometry_fix.par'

