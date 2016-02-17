# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""

Define container classes for database parameters.

*Source location*:
    $MANGADAP_DIR/python/mangadap/par/database_definitions.py

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
    from mangadap.par.ParSet import ParSet

*Class usage examples*:

    .. todo::

        Add some usage comments here!

*Revision history*:
    | **16 Jun 2015**: Original implementation by K. Westfall (KBW)
    | **29 Jan 2016**: (KBW) Included sres_ext and log10 keywords in
        TemplateLibraryDef.  Documentation.
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

    Args:
        key (str): Keyword to distinguish the template library.
        file_search (str): Search string used by glob to find the 1D
            fits spectra to include in the template library.
        fwhm (int or float): FWHM of the resolution element in
            angstroms.
        sres_ext (str): Extension in the fits files with measurements of
            the spectral resolution as a function of wavelength.
        in_vacuum (bool): Flag that the wavelengths of the spectra are
            in vacuum, not air.
        wave_limit (numpy.ndarray): 2-element array with the starting
            and ending wavelengths for the valid spectral range of the
            templates.
        lower_flux_limit (int or float): Minimum valid flux in the
            template spectra.
        log10 (bool): Flag that the template spectra have been binned
            logarithmically in wavelength.

    """
    def __init__(self, key, file_search, fwhm, sres_ext, in_vacuum, wave_limit, lower_flux_limit,
                 log10): 
        # Perform some checks of the input
        in_fl = [ int, float ]
        
        pars =   [ 'key', 'file_search', 'fwhm', 'sres_ext', 'in_vacuum',  'wave_limit',
                        'lower_flux_limit', 'log10' ]
        values = [   key,   file_search,   fwhm,   sres_ext,   in_vacuum,    wave_limit,
                          lower_flux_limit,   log10 ]
        dtypes = [   str,           str,  in_fl,        str,        bool, numpy.ndarray,
                                     in_fl,    bool ]

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


