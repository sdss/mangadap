"""

Define some utility classes used to hold parameters.

*Source location*:
    $MANGADAP_DIR/python/mangadap/util/par.py

*Python2/3 compliance*::

    from __future__ import division
    from __future__ import print_function
    from __future__ import absolute_import
    from __future__ import unicode_literals
    
    import sys
    if sys.version > '3':
        long = int

*Imports*::

    import numpy
    from mangadap.util.options import binning_options, bin_weighting_options
    from mangadap.util.options import spectral_analysis_options

*Class usage examples*:

    .. todo::

        Add some usage comments here!

*Revision history*:
    | **16 Jun 2015**: Original implementation by K. Westfall (KBW)

"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import sys
if sys.version > '3':
    long = int

import numpy
from mangadap.par.ParSet import ParSet

__author__ = 'Kyle B. Westfall'


class SNRPar(ParSet):
    """
    Base class with parameters used to calculate S/N.
    """
    def __init__(self, snr_lamrange=None):

        pars =     [ 'snr_lamrange' ]
        values =   [   snr_lamrange ]
        defaults = [ numpy.ndarray([5560.00, 6942.00]) ]
        dtypes =   [  numpy.ndarray ]

        ParSet.__init__(self, pars, values=values, defaults=defaults, dtypes=dtypes)


