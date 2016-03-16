# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""

Provides a set of functions that define and return options available to
the MaNGA DAP.

*Source location*:
    $MANGADAP_DIR/python/mangadap/proc/options.py

*Python2/3 compliance*::

    from __future__ import division
    from __future__ import print_function
    from __future__ import absolute_import
    from __future__ import unicode_literals
    
    import sys
    if sys.version > '3':
        long = int

*Imports*::

*Revision history*:
    | **18 Jun 2015**: Original implementation by K. Westfall (KBW)
    | **15 Mar 2016**: (KBW) Add DRP output options

"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import sys
if sys.version > '3':
    long = int

__author__ = 'Kyle B. Westfall'

def drp_wave_sampling_options():
    return [ 'LIN', 'LOG' ]

def drp_3dmode_options():
    return [ 'CUBE', 'RSS' ]

def binning_options():
    return [ 'NONE', 'SNR', 'RADIAL', 'ALL' ]

def bin_weighting_options():
    return [ 'Uniform', 'Optimal' ]

def spectral_analysis_options():
    return [ 'STRCNT', 'EMLINE', 'SPINDX' ]


