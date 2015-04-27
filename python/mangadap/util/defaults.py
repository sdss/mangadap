
"""
defaults.py

Functions that define and return defaults used by the DAP, such as paths
and file names.

REVISION HISTORY:
    23 Apr 2015: Original implementation by K. Westfall (KBW)
"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import sys
if sys.version > '3':
    long = int

from os import environ

__author__ = 'Kyle B. Westfall'


def default_drp_version():
    """
    Return the DRP version defined by the environmental variable.
    """
    return environ['MANGADRP_VER']


def default_redux_path(drpver):
    """Return the main output path for the DRP products."""

    # Make sure the DRP version is set
    if drpver is None:
        drpver = default_drp_version()

    return os.path.join(environ['MANGA_SPECTRO_REDUX'], drpver)


def default_drp_directory_path(redux_path, plate):
    """Return the exact directory path with the file."""

    # Make sure the DRP version is set
    if redux_path is None:
        redux_path = default_redux_path()

    return os.path.join(redux_path, str(plate), 'stack')


def default_dap_source():
    """Return the root path to the DAP source."""
    return environ['MANGADAP_DIR']


def default_dap_version():
    """Return the DAP version defined by the environmental variable."""
    return environ['MANGADAP_VER']


def default_analysis_path(drpver, dapver):
    """Return the directory path used by the DAP."""
    # Make sure the DRP version is set
    if drpver is None:
        drpver = default_drp_version()
    # Make sure the DAP version is set
    if dapver is None:
        dapver = default_dap_version()
    return os.path.join(environ['MANGA_SPECTRO_ANALYSIS'], drpver, dapver)


def default_dap_directory_path(drpver, dapver, analysis_path, plate, ifudesign):
    """Return the directory path used by the DAP."""
    # Make sure the DRP version is set
    if analysis_path is None:
        analysis_path = default_analysis_path(drpver, dapver)
    return os.path.join(analysis_path, str(plate), str(ifudesign))


def default_dap_par_file(drpver, dapver, analysis_path, directory_path, plate, ifudesign, mode):
    if directory_path is None:
        directory_path = default_dap_directory_path(drpver, dapver, analysis_path, plate, ifudesign)
    par_file = 'mangadap-{0}-{1}-LOG{2}.par'.format(plate, ifudesign, mode)
    return os.path.join(directory_path, par_file)

    
def default_dap_plan_file(drpver, dapver, analysis_path, directory_path, plate, ifudesign, mode):
    if directory_path is None:
        directory_path = default_dap_directory_path(drpver, dapver, analysis_path, plate, ifudesign)
    plan_file = 'manga-{0}-{1}-LOG{2}-dapplan.par'.format(plate, ifudesign, mode)
    return os.path.join(directory_path, plan_file)

    

