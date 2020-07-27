# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Implements a BitMask derived class that simplifies the implementation of
the derived classes used by the DAP that all define their bits using the
files in the configuration directory.

----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import os

from .bitmask import BitMask
from ..config import defaults


class DAPBitMask(BitMask):
    r"""
    Class that simplifies the derived classes used by the DAP by always
    sourcing an ini file in the configuration directory.

    This class is an *abstract class* that has no direct use beyond
    being a base class.
    """
    cfg_root = None
    def __init__(self):
        if self.cfg_root is None:
            raise ValueError('Class does not define the configuration file root name!')
        tmp = BitMask.from_ini_file(os.path.join(defaults.dap_config_root(), 'bitmasks',
                                                 '{0}.ini'.format(self.cfg_root)))
        keys, descr = tmp._init_objs()
        super(DAPBitMask, self).__init__(keys, descr=descr)

