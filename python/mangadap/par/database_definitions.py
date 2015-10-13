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


class TemplateLibraryDef(ParSet):
    """
    Class with parameters used to define the template library.
    Options and defaults in ParSet base class are set to None.
    """
    def __init__(self, key, file_search, fwhm, in_vacuum, wave_limit, lower_flux_limit): 
        # Perform some checks of the input
        in_fl = [ int, float ]
        
        pars =   [ 'key', 'file_search', 'fwhm', 'in_vacuum',  'wave_limit', 'lower_flux_limit' ]
        values = [   key,   file_search,   fwhm,   in_vacuum,    wave_limit,   lower_flux_limit ]
        dtypes = [   str,           str,  in_fl,        bool, numpy.ndarray,              in_fl ]

        ParSet.__init__(self, pars, values=values, dtypes=dtypes)


class EmissionLineDef(ParSet):
    """
    Class with parameters used to define an emission-line database.
    Options and defaults in ParSet base class are set to None.
    """
    def __init__(self, key, file_path):

        pars =   [ 'key', 'file_path' ]
        values = [   key,   file_path ]
        dtypes = [   str,         str ]

        ParSet.__init__(self, pars, values=values, dtypes=dtypes)


class SpectralIndexDef(ParSet):
    """
    Class with parameters used to define a spectral-index database.
    Options and defaults in ParSet base class are set to None.
    """
    def __init__(self, key, file_path):

        pars =   [ 'key', 'file_path' ]
        values = [   key,   file_path ]
        dtypes = [   str,         str ]

        ParSet.__init__(self, pars, values=values, dtypes=dtypes)


