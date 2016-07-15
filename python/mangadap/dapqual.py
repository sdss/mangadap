# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Defines the DAP quality bit mask class.

*License*:
    Copyright (c) 2015, SDSS-IV/MaNGA Pipeline Group
        Licensed under BSD 3-clause license - see LICENSE.rst

*Source location*:
    $MANGADAP_DIR/python/mangadap/dapqual.py

*Imports and python version compliance*:
    ::

*Class usage examples*:
    Add some usage comments here!

*Revision history*:
    | **13 Jul 2016**: Original Implementation by K. Westfall (KBW);
        originally in $MANGADAP_DIR/python/mangadap/dapmaps.py

.. todo::
    Do something different than add an empty extension when the data are
    not available?
"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import os

from .util.bitmask import BitMask
from .config.defaults import default_dap_source

__author__ = 'Kyle Westfall'

class DAPQualityBitMask(BitMask):
    """
    .. todo::
        - Force read IDLUTILS version as opposed to internal one?
    """
    def __init__(self, dapsrc=None):
        dapsrc = default_dap_source() if dapsrc is None else str(dapsrc)
        BitMask.__init__(self, ini_file=os.path.join(dapsrc, 'python', 'mangadap', 'config',
                                                     'bitmasks', 'dap_quality_bits.ini'))


