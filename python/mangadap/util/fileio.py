"""

Provides a set of file I/O routines.

*Source location*:
    $MANGADAP_DIR/python/mangadap/util/fileio.py

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
    from astropy.io import fits

*Revision history*:
    | **27 May 2015**: Original implementation by K. Westfall (KBW)

"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import sys
if sys.version > '3':
    long = int

import numpy
from astropy.io import fits

def readfits_1dspec(filename):
    """
    Read a 1D fits spectrum and return two vectors with the wavelength
    and flux.

    Args:
        filename (str): Name of the file to read.

    Returns:
        numpy.ndarray : Two numpy.float64 arrays with the wavelength and
        flux read for the spectrum.

    Raises:
        Exception: Raised if the input fits file has more than one
            extension or its primary extension has more than two
            dimensions.
    """
    hdu = fits.open(filename, mode='readonly')

    if (len(hdu)) != 1:
        raise Exception('{0} has more than one extension.'.format(filename))
    
    if hdu[0].header['NAXIS'] != 1:
        raise Exception('{0} is not one dimensional!'.format(filename))
    
    spec = numpy.copy(hdu[0].data).astype(numpy.float64)

    crval = hdu[0].header['CRVAL1']
    crpix = hdu[0].header['CRPIX1']
    cdelt = hdu[0].header['CDELT1']

    hdu.close()

    npix = spec.size
    wave = ((numpy.arange(1.0,npix+1) - crpix)*cdelt + crval).astype(numpy.float64)

    return wave, spec


