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


