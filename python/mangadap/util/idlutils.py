"""

Contains transcriptions of some IDLUTILS functions to python.

*Source location*:
    $MANGADAP_DIR/python/mangadap/util/idlutils.py

*Python2/3 compliance*::

    from __future__ import division
    from __future__ import print_function
    from __future__ import absolute_import
    
    import sys
    if sys.version > '3':
        long = int

*Imports*::

    import numpy

*Revision history*:
    | **23 Apr 2015**: Original implementation by K. Westfall (KBW)
    | **20 May 2015**: (KBW) Documentation and Sphinx tests
"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import sys
if sys.version > '3':
    long = int

import numpy

__author__ = 'Kyle B. Westfall'

#-----------------------------------------------------------------------
def airtovac(wave_air):
    """ 
    Wavelengths are corrected for the index of refraction of air under
    standard conditions.  Wavelength values below 2000 A will not be
    altered.  Uses formula from Ciddor 1996, Applied Optics 62, 958.

    Args:
        wave_air (int or float): Wavelength in Angstroms, scalar or
            vector. If this is the only parameter supplied, it will be
            updated on output to contain double precision vacuum
            wavelength(s). 

    Returns:
        numpy.float64 : The wavelength of the line in vacuum.

    Example:
        If the air wavelength is  W = 6056.125 (a Krypton line), then
        :func:`airtovac` returns vacuum wavelength of W = 6057.8019.
 
    *Revision history*:
        | Written W. Landsman                November 1991
        | Use Ciddor (1996) formula for better accuracy in the infrared 
        |   Added optional output vector, W Landsman Mar 2011
        | Iterate for better precision W.L./D. Schlegel  Mar 2011
        | Transcribed to python, K.B. Westfall Apr 2015

    .. note::
        Take care within 1 A of 2000 A.   Wavelengths below 2000 A *in
        air* are not altered.       

    """

    # Copy the data
    wave_vac = wave_air.astype(numpy.float64) if hasattr(wave_air, "__len__") else float(wave_air)
    g = wave_vac > 2000.0                           # Only modify above 2000 A
    Ng = numpy.sum(g)
    
    if Ng > 0:
        # Handle both arrays and scalars
        if hasattr(wave_air, "__len__"):
            _wave_air = wave_air[g].astype(numpy.float64)
            _wave_vac = wave_vac[g]
        else:
            _wave_air = float(wave_air)
            _wave_vac = float(wave_vac)

        for i in range(0,2):
            sigma2 = numpy.square(1.0e4/_wave_vac)     #Convert to wavenumber squared
            fact = 1.0 + 5.792105e-2/(238.0185 - sigma2) + 1.67917e-3/(57.362 - sigma2)
            _wave_vac = _wave_air*fact

        if hasattr(wave_air, "__len__"):        # Save the result
            wave_vac[g] = _wave_vac
        else:
            wave_vac = _wave_vac

    return wave_vac


#-----------------------------------------------------------------------


