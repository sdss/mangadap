"""

Provides an AnalysisPlan class used to create, read, and edit MaNGA DAP
analysis plans.

This is mostly just an interface with a yanny object.

*Source location*:
    $MANGADAP_DIR/python/mangadap/util/AnalysisPlan.py

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

*Class usage examples*:

    .. todo::

        Add some usage comments here!

*Revision history*:
    | **19 Mar 2016**: Original implementation by K. Westfall (KBW)

"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import sys
if sys.version > '3':
    long = int

import numpy
from ..util.yanny import yanny

__author__ = 'Kyle B. Westfall'

class AnalysisPlan:
    """
    Generic class to handle and manipulate MaNGA DAP analysis plans.

    SignalToNoisePlan
        - wavelength range(s)
        - covariance
        - clobber
    BinningPlan
        - type: All
        - type: None
        - type: Voronoi
            - S/N target
            - S/N calculation
        - type: Radial
            - geometry: cx, cy, pa, boa
            - rscale
            - sampling: rs (in scale units), re (in scale units), nr, rlog
        - velocity register
        - weighting
        - S/N Threshold
        - rejection
        - covariance
        - clobber
    StellarContinuumPlan: can be NULL
        - (binned) spectra
        - S/N Threshold
        - template library
        - pPXFPlan
            - moments
            - degree
            - mdegree
            - reddening
            - oversample
            ...
        - clobber
    EmissionLinePlan: can be NULL
        - (binned) spectra
        - S/N Threshold
        - emission-line moment database
        - emission-line database
        - integrated flux limit
        - allow negative amplitude?
        - clobber
    SpectralIndexPlan: can be NULL
        - (binned) spectra
        - S/N Threshold
        - absorption-line index database
        - bandhead index database
        - pseudo-continuum limit?
        - clobber

    - index
    - LIN vs LOG
    - CUBE vs. RSS
    - prior
    - execute


    open empty, open from file
    add new plan, and iteratively do so
    delete (all) plan(s)
    print plans to screen or log file
    write plans

    """



