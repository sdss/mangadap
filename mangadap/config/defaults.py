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
from pkg_resources import resource_filename


def dap_source_dir():
    from pkg_resources import resource_filename
    return Path(resource_filename('mangadap', '')).resolve()


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


# USED BY:
#   - scripts/dapall_qa.py
#   - scripts/find_repeat_observations.py
#   - survey/dapall.py
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


# USED BY:
#   - scripts/rundap.py
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

    
# USED BY:
#   - scripts/rundap.py
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

    

