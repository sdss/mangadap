"""

Provides a set of functions that define and return options available to
the MaNGA DAP.

*Source location*:
    $MANGADAP_DIR/python/mangadap/util/options.py

*Python2/3 compliance*::

    from __future__ import division
    from __future__ import print_function
    from __future__ import absolute_import
    from __future__ import unicode_literals
    
    import sys
    if sys.version > '3':
        long = int

*Imports*::

    import os.path
    from os import environ
    import numpy

*Revision history*:
    | **18 Jun 2015**: Original implementation by K. Westfall (KBW)

"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import sys
if sys.version > '3':
    long = int

import os.path
from os import environ
import numpy

__author__ = 'Kyle B. Westfall'


def binning_options():
    return [ 'NONE', 'SNR', 'RADIAL', 'ALL' ]

def bin_weighting_options():
    return [ 'Uniform', 'Optimal' ]

def spectral_analysis_options():
    return [ 'STRCNT', 'EMLINE', 'SPINDX' ]


