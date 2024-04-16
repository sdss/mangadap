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

from importlib import resources

def dap_source_dir():
    return resources.files('mangadap').parent


def dap_data_root():
    """Return the root directory with the DAP data."""
    return dap_source_dir() / 'mangadap' / 'data'


def dap_config_root():
    """Return the root directory with the DAP config data."""
    return dap_source_dir() / 'mangadap' / 'config'


def sdss_maskbits_file():
    """Return the path to the sdss maskbits yanny file."""
    maskbits_file = dap_data_root() / 'sdss' / 'sdssMaskbits.par'
    return maskbits_file if maskbits_file.is_file() else None


