# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
r"""
Provides a set of functions that define and return defaults used by the
MaNGA DAP, such as paths and file names.

Revision history
----------------

    | **23 Apr 2015**: Original implementation by K. Westfall (KBW)
    | **19 May 2015**: (KBW) Documentation and Sphinx tests
    | **20 May 2015**: (KBW) TODO: Change the default names of the par
        and output fits files
    | **15 Jun 2015**: (KBW) Added default functions moved from
        :class:`mangadap.drpfile`
    | **16 Jun 2015**: (KBW) Added wavelength limits and lower flux
        limit to
        :func:`mangadap.config.inputdata.available_template_libraries`
    | **17 Jun 2015**: (KBW) Moved declarations of template library keys
        to its own function: *default_template_library_keys*, and
        edited :func:`available_template_libraries` accordingly
    | **27 Aug 2015**: (KBW) Changed the name of the plan file; added
        :func:`default_dap_file_root` based on file_root() from
        :class:`mangadap.survey.rundap`.
    | **28 Aug 2015**: (KBW) Added :func:`default_manga_fits_root` and
        :func:`default_cube_covariance_file`
    | **07 Oct 2015**: (KBW) Adjusted for changes to naming of the
        template library database definitions.  Added M11-STELIB-ZSOL
        library.
    | **29 Jan 2016**: (KBW) Changed
        :func:`manga.config.inputdata.available_template_libraries` to
        use configparser ini files to define each template library.
    | **03 Feb 2016**: (KBW) Added checks for required environmental
        variables.
    | **17 Feb 2016**: (KBW) Added try/except blocks for importing
        ConfigParser.
    | **16 Mar 2016**: (KBW) Created :mod:`mangadap.config.inputdata`
        and moved **default_template_libraries** there (and changed it
        to
        :func:`mangadap.config.inputdata.available_template_libraries`.
        No longer need ConfigParser here.
    | **12 Jul 2016**: (KBW) Changed :func:`default_manga_fits_root` to
        accommodate MAPS output name.
    | **05 May 2017**: (KBW) :func:`default_dap_version` returns
        internal version if MANGADAP_VER environmental variable is not
        defined.
    | **01 Mar 2018**: (KBW) Added :func:`default_redshift_fix_file`.

----        

.. include license and copyright
.. include:: ../copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../links.rst
"""

import os
import glob
import numpy

from ..util.exception_tools import check_environment_variable
from mangadap import __version__, dap_source_dir

def dap_source_dir():
    """Return the root path to the DAP source directory."""
    dirlist = os.path.dirname(os.path.abspath(__file__)).split('/')[:-3]
    return os.path.join(os.sep, *dirlist) if dirlist[0] == '' else os.path.join(*dirlist)


def idlutils_dir():
    """Return the default IDLUTILS directory."""
    check_environment_variable('IDLUTILS_DIR')
    return os.environ['IDLUTILS_DIR']


def sdss_maskbits_file():
    """Return the path to the sdss maskbits yanny file."""
    maskbits_file = os.path.join(dap_source_dir(), 'data', 'sdss', 'sdssMaskbits.par')
    if os.path.isfile(maskbits_file):
        return maskbits_file

    try:
        maskbits_file = os.path.join(idlutils_dir(), 'data', 'sdss', 'sdssMaskbits.par')
        if os.path.isfile(maskbits_file):
            return maskbits_file
    except:
        return None


def default_drp_version():
    """
    Return the DRP version defined by the environmental variable
    MANGADRP_VER.
    """
    check_environment_variable('MANGADRP_VER')
    return os.environ['MANGADRP_VER']


def default_redux_path(drpver=None):
    """
    Return the main output path for the DRP products using the
    environmental variable MANGA_SPECTRO_REDUX.

    Args:
        drpver (str): (**Optional**) DRP version.  Default is to use
            :func:`default_drp_version`.

    Returns:
        str: Path to reduction directory
    """
    # Make sure the DRP version is set
    if drpver is None:
        drpver = default_drp_version()
    check_environment_variable('MANGA_SPECTRO_REDUX')
    return os.path.join(os.path.abspath(os.environ['MANGA_SPECTRO_REDUX']), drpver)


def default_drp_directory_path(plate, drpver=None, redux_path=None):
    """
    Return the exact directory path with the DRP file.

    Args:
        plate (int): Plate number
        drpver (str): (**Optional**) DRP version.  Default is to use
            :func:`default_drp_version`.
        redux_path (str): (**Optional**) Path to the root reduction
            directory.  Default is to use :func:`default_redux_path`.

    Returns:
        str: Path to the directory with the 3D products of the DRP
    """
    # Make sure the redux path is set
    _redux_path = default_redux_path(drpver=drpver) \
                        if redux_path is None else os.path.abspath(redux_path)
    return os.path.join(_redux_path, str(plate), 'stack')


def default_drpall_file(drpver=None, redux_path=None):
    """
    Return the path to the DRPall file.

    Args:
        drpver (:obj:`str`, optional):
            DRP version.  Default is to use :func:`default_drp_version`.
        redux_path (:obj:`str`, optional):
            Path to the root reduction directory.  Default is to use
            :func:`default_redux_path`.

    Returns:
        :obj:`str`: Full path to the DRPall fits file.
    """
    _drpver = default_drp_version() if drpver is None else drpver
    _redux_path = default_redux_path(drpver=_drpver) \
                        if redux_path is None else os.path.abspath(redux_path)
    return os.path.join(_redux_path, 'drpall-{0}.fits'.format(_drpver))

def default_dapall_file(drpver=None, dapver=None, analysis_path=None):
    """
    Return the path to the DAPall file.

    Args:
        drpver (:obj:`str`, optional):
            DRP version.  Default is to use :func:`default_drp_version`.
        dapver (:obj:`str`, optional):
            DAP version.  Default is to use :func:`default_dap_version`.
        analysis_path (:obj:`str`, optional):
            Path to the root analysis directory.  Default is to use
            :func:`default_analysis_path`

    Returns:
        :obj:`str`: Full path to the DAPall fits file.
    """
    _drpver = default_drp_version() if drpver is None else drpver
    _dapver = default_dap_version() if dapver is None else dapver
    _analysis_path = default_analysis_path(drpver=_drpver, dapver=_dapver) \
                        if analysis_path is None else os.path.abspath(analysis_path)
    return os.path.join(_analysis_path, 'dapall-{0}-{1}.fits'.format(_drpver, _dapver))

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


def default_dap_version():
    """
    Return the DAP version defined by the environmental variable
    MANGADAP_VER.  If that environmental variable does not exist,
    `mangadap.__version__` is returned.
    """
    no_environ_var=False
    try:
        check_environment_variable('MANGADAP_VER')
    except EnvironmentError as e:
        warnings.warn('$MANGADAP_VER undefined in environment; returning internal version')
        no_environ_var = True
    return __version__ if no_environ_var else os.environ['MANGADAP_VER']


def default_analysis_path(drpver=None, dapver=None):
    """
    Return the main output path for the DAP using the environmental
    variable MANGA_SPECTRO_ANALYSIS.

    Args:
        drpver (str): (**Optional**) DRP version.  Default is to use
            :func:`default_drp_version`.
        dapver (str): (**Optional**) DAP version.  Default is to use
            :func:`default_dap_version`.

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
    return os.path.join(os.path.abspath(os.environ['MANGA_SPECTRO_ANALYSIS']), drpver, dapver)


def default_dap_common_path(plate=None, ifudesign=None, drpver=None, dapver=None,
                            analysis_path=None):
    """
    Return the path to the path to the directory with data common to
    multiple binning schemes.

    Args:
        plate (int): (**Optional**) Plate number, for reference
            directory of a specific plate.
        ifudesign (int): (**Optional**)  IFU design number.
        drpver (str): (**Optional**) DRP version.  Default is to use
            :func:`default_drp_version`.
        dapver (str): (**Optional**) DAP version.  Default is to use
            :func:`default_dap_version`.
        analysis_path (str): (**Optional**) Path to the root analysis
            directory.  Default is to use :func:`default_analysis_path`

    Returns:
        str: Path to the directory with DAP reference files

    Raises:
        ValueError: Raised if IFU design is provided and plate is not.
    """
    # For heirarchy, if ifudesign is given, plate must also be given
    if plate is None and ifudesign is not None:
        raise ValueError('For IFU design subdirectory, must provide plate number.')

    # Get the main analysis path
    _analysis_path = default_analysis_path(drpver=drpver, dapver=dapver) \
                        if analysis_path is None else os.path.abspath(analysis_path)

    output_path = os.path.join(_analysis_path, 'common')
    if plate is None:
        return output_path
    output_path = os.path.join(output_path, str(plate))
    return output_path if ifudesign is None else os.path.join(output_path, str(ifudesign))


#def default_dap_method(plan=None, binned_spectra=None, stellar_continuum=None,
#                       binning_method=None, continuum_method=None):
#    """
#    Return the method for a provided plan.  The name currently uses the
#    keyword identifiers for the binned spectra and the stellar continuum
#    models.
#    """
#    if plan is not None:
#        if not isinstance(plan['bin_key'], str) and len(plan['bin_key']) > 1:
#            raise ValueError('Must only provide a single plan.')
#        if not isinstance(plan['continuum_key'], str) and len(plan['continuum_key']) > 1:
#            raise ValueError('Must only provide a single plan.')
#        _binning_method = str(plan['bin_key'])
#        _continuum_method = str(plan['continuum_key'])
#    else:
##        binning_method = 'None' if binned_spectra is None else binned_spectra.method['key']
#        if binned_spectra is not None:
#            _binning_method = binned_spectra.method['key']
#        elif binning_method is not None:
#            _binning_method = binning_method
#        else:
#            _binning_method = 'None'
#        if stellar_continuum is not None:
#            _continuum_method = stellar_continuum.method['key']
#        elif continuum_method is not None:
#            _continuum_method = continuum_method
#        else:
#            _continuum_method = 'None'
#    return '{0}-{1}'.format(_binning_method, _continuum_method)


def default_dap_method(binning_method, stellar_continuum_templates, emission_line_model_templates):

    """
    Construct the "DAPTYPE" based on the analysis plan.

    The construction is based on directly provided strings for
    `binning_method`, `stellar_continuum_templates`, and
    `emission_line_model_templates`.  The DAPTYPE is constructed as
    these three strings separated by dashes.  E.g., DAPTYPE =
    'HYB10-MILESHC-MILESHC' when the binning method is HYB10 and the
    MILESHC library is used as templates for both the stellar and
    emission-line fitting.

    With a analysis plan file, run the following to get the list of
    daptypes::

        from mangadap.par.analysisplan import AnalysisPlanSet
        from mangadap.proc.spatiallybinnedspectra import SpatiallyBinnedSpectra
        from mangadap.proc.stellarcontinuummodel import StellarContinuumModel
        from mangadap.proc.emissionlinemodel import EmissionLineModel
        from mangadap.config.defaults import default_dap_method

        plans = AnalysisPlanSet.from_par_file(plan_file)
        daptypes = []
        for plan in plans:
            bin_method = SpatiallyBinnedSpectra.define_method(plan['bin_key'])
            sc_method = StellarContinuumModel.define_method(plan['continuum_key'])
            el_method = EmissionLineModel.define_method(plan['elfit_key'])
            daptypes += [default_dap_method(bin_method['key'], sc_method['template_library'],
                                            el_method['continuum_templates'])]

    Args:
        binning_method (:obj:`str`):
            String used to define the binning method.
        stellar_continuum_templates (:obj:`str`):
            String defining the template library used in the stellar
            continuum (stellar kinematics) modeling.

        emission_line_model_templates (:obj:`str`):
            String defining the template library used in the
            emission-line modeling.  Can be None, meaning that the
            continuum templates are the same as used for the
            stellar-continuum modeling.
    Returns:
        str: The string representation of the analysis method.

    """
    _eltpl = stellar_continuum_templates if emission_line_model_templates is None \
                                                else emission_line_model_templates
    return '{0}-{1}-{2}'.format(binning_method, stellar_continuum_templates, _eltpl)


def default_dap_method_path(method, plate=None, ifudesign=None, qa=False, ref=False,
                            drpver=None, dapver=None, analysis_path=None):
    """
    Return the path to the designated subdirectory built using the plan
    key identifiers or directly using the provided method.

    The "method" identifying each plan is currently build using the
    keywords for the binning type and the continuum-fitting key.
    Nominally, the latter should include the template set used.

    Args:
        method (str): String defining the method identifier for a set of
            DAP output files.  These should be built using
            :func:`default_dap_method`.
        plate (int): (**Optional**) Plate number.
        ifudesign (int): (**Optional**)  IFU design number.
        qa (bool): (**Optional**) Give the path to the qa/ subdirectory
        ref (bool): (**Optional**) Give the path to the ref/
            subdirectory
        drpver (str): (**Optional**) DRP version.  Default is to use
            :func:`default_drp_version`.
        dapver (str): (**Optional**) DAP version.  Default is to use
            :func:`default_dap_version`.
        analysis_path (str): (**Optional**) Path to the root analysis
            directory.  Default is to use :func:`default_analysis_path`

    Returns:
        str: Path to the plan subdirectory

    Raises:
        ValueError: Raised if IFU design is provided and plate is not,
            or if either qa or ref are true and one or both of plate and
            IFU design are not provided..
    """
    # For heirarchy, if ifudesign is given, plate must also be given
    if (plate is not None and ifudesign is None) or (plate is None and ifudesign is not None):
        raise ValueError('Must provide both plate and ifudesign, if providing either.')
    if qa and ref:
        raise ValueError('Cannot provide path for both qa and ref directory.  Pick one.')

    # Get the main analysis path
    _analysis_path = default_analysis_path(drpver=drpver, dapver=dapver) \
                        if analysis_path is None else os.path.abspath(analysis_path)

    # Build the plan subirectory
    output_path = os.path.join(_analysis_path, method)
    if plate is None:
        return output_path
    output_path = os.path.join(output_path, str(plate), str(ifudesign))
    if not qa and not ref:
        return output_path
    return os.path.join(output_path, ('qa' if qa else 'ref'))


def default_manga_fits_root(plate, ifudesign, mode=None):
    """
    Generate the main root name for the output MaNGA fits files for a
    given plate/ifudesign/mode.

    Args:
        plate (int): Plate number
        ifudesign (int): IFU design number
        mode (str): (**Optional**) Mode of the output fits file.
            Options are: 'LINCUBE', 'LINRSS', 'LOGCUBE', 'LOGRSS', or
            'MAPS'.  Default is that no mode is included in the name.

    Returns:
        str: Root name for a MaNGA fits file:
        `manga-[PLATE]-[IFUDESIGN]-[MODE]`

    Raises:

        ValueError: Raised if mode is not a valid option.

    """
    if mode not in [ None, 'LINCUBE', 'LINRSS', 'LOGCUBE', 'LOGRSS', 'MAPS' ]:
        raise ValueError('Do not understand mode={0}.'.format(mode))
    return 'manga-{0}-{1}'.format(plate, ifudesign) if mode is None else \
                    'manga-{0}-{1}-{2}'.format(plate, ifudesign, mode)


def default_dap_file_root(plate, ifudesign, mode=None):
    """
    Generate the root name of the MaNGA DAP parameter and script files
    for a given plate/ifudesign/mode.
    
    Args:
        plate (int): Plate number
        ifudesign (int): IFU design number
        mode (str): (**Optional**) Mode of the DRP reduction; either RSS
            or CUBE.  If None (default), the mode is excluded from the
            file root.

    Returns:
        str: Root name for the DAP file: `mangadap-[PLATE]-[IFUDESIGN]`
        or `mangadap-[PLATE]-[IFUDESIGN]-LOG[MODE]`
    """
    return 'mangadap-{0}-{1}'.format(plate, ifudesign) if mode is None else \
                    'mangadap-{0}-{1}-LOG{2}'.format(plate, ifudesign, mode)


def default_dap_par_file(plate, ifudesign, mode, partype='input', drpver=None, dapver=None,
                         analysis_path=None, directory_path=None):
    """
    Return the full path to a par file used by the DAP to analyze the
    specified DRP output file.

    Args:
        plate (int): Plate number
        ifudesign (int): IFU design number
        mode (str): Mode of the DRP reduction; either RSS or CUBE
        partype (str):  An "unregulated" type for the parameter file.
            The default is ``'input'``, signifying the input set of
            observational parameters; see
            :class:`mangadap.par.obsinput.ObsInputPar`.
        drpver (str): (**Optional**) DRP version.  Default is to use
            :func:`default_drp_version`.
        dapver (str): (**Optional**) DAP version.  Default is to use
            :func:`default_dap_version`.
        analysis_path (str): (**Optional**) Path to the root analysis
            directory.  Default is to use :func:`default_analysis_path`
        directory_path (str): (**Optional**) Path to the directory with
            the DAP output files.  Default is to use
            :func:`default_dap_common_path`

    Returns:
        str: Full path to the DAP par file
    """
    # Make sure the directory path is defined
    _directory_path = default_dap_common_path(plate=plate, ifudesign=ifudesign, drpver=drpver,
                                             dapver=dapver, analysis_path=analysis_path) \
                            if directory_path is None else os.path.abspath(directory_path)
    # Set the name of the par file; put this in its own function?
    par_file = '{0}-{1}.par'.format(default_dap_file_root(plate, ifudesign, mode), partype)
    return os.path.join(_directory_path, par_file)

    
def default_dap_plan_file(drpver=None, dapver=None, analysis_path=None):
    """
    Return the full path to the DAP plan file.

    Args:
        drpver (str): (**Optional**) DRP version.  Default is to use
            :func:`default_drp_version`.
        dapver (str): (**Optional**) DAP version.  Default is to use
            :func:`default_dap_version`.
        analysis_path (str): (**Optional**) Path to the root analysis
            directory.  Default is to use :func:`default_analysis_path`

    Returns:
        str: Full path to the DAP plan file
    """
    # Get the main analysis path
    if dapver is None:
        dapver = default_dap_version()

    _analysis_path = default_analysis_path(drpver=drpver, dapver=dapver) \
                        if analysis_path is None else os.path.abspath(analysis_path)
    
    # Set the name of the plan file
    plan_file = 'mangadap-plan-{0}.par'.format(dapver)
    return os.path.join(_analysis_path, plan_file)


def default_dap_file_name(plate, ifudesign, output_mode, mode=None, compressed=True):
    """
    Return the name of the DAP output fits file.

    Args:
        plate (int): Plate number
        ifudesign (int): IFU design number
        output_mode (str) : Output mode designation
        mode (str): (**Optional**) Mode of the output fits file.
            Options are: 'LINCUBE', 'LINRSS', 'LOGCUBE', 'LOGRSS', or
            'MAPS'.  Default is that no mode is included in the name.
        compressed (bool): (**Optional**) Append '.gz' to the output
            file name.  Default is True.

    Returns:
        str: Name of the DAP output file.
    """
    # Number of spaces provided for iteration number is 3
    root = default_manga_fits_root(plate, ifudesign, mode=mode)
    return ('{0}-{1}.fits.gz' if compressed else '{0}-{1}.fits').format(root, output_mode)


def default_plate_target_files():
    """
    Return the default plateTarget files in mangacore and their
    associated catalog indices.  The catalog indices are determined
    assuming the file names are of the form::

        'plateTargets-{0}.par'.format(catalog_id)

    Returns:
        numpy.array: Two arrays: the first contains the identified
        plateTargets files found using the default search string, the
        second provides the integer catalog index determined for each
        file.
    """
    # Default search string
    check_environment_variable('MANGACORE_DIR')
    search_str = os.path.join(os.path.abspath(os.environ['MANGACORE_DIR']), 'platedesign',
                              'platetargets', 'plateTargets*.par')
    file_list = glob.glob(search_str)                       # List of files
    nfiles = len(file_list)
    trgid = numpy.zeros(nfiles, dtype=numpy.int)            # Array to hold indices
    for i in range(nfiles):
        suffix = file_list[i].split('-')[1]                 # Strip out the '{i}.par'
        trgid[i] = int(suffix[:suffix.find('.')])           # Strip out '{i}'
    
    return numpy.array(file_list), trgid


def default_redshift_fix_file():
    """
    Return the path to the default redshift fix file.

    Returns:
        str: Expected path to the redshift-fix parameter file.
    """

    # Default search string
    dapsrc = dap_source_dir()
    return os.path.join(dapsrc, 'data', 'fix', 'redshift_fix.par')

