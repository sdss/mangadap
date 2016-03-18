# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""

Provides a set of file I/O routines.

*License*:
    Copyright (c) 2015, Kyle B. Westfall
    Licensed under BSD 3-clause license - see LICENSE.rst

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
    | **29 Jan 2016**: (KBW) Added :func:`read_template_spectrum`
    | **01 Feb 2016**: (KBW) Moved wavelength calculation to a common
        function.
    | **09 Feb 2016**: (KBW) Added :func:`writefits_1dspec`.

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


def wavelength_vector(npix, header, log10=False):
    crval = header['CRVAL1']
    crpix = header['CRPIX1']
    cdelt = header['CDELT1']
    wave = ((numpy.arange(1.0,npix+1) - crpix)*cdelt + crval).astype(numpy.float64)
    return numpy.power(10., wave) if log10 else wave


def readfits_1dspec(filename, log10=False):
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
    wave = wavelength_vector(spec.size, hdu[0].header, log10=log10)
    hdu.close()
    return wave, spec


def writefits_1dspec(ofile, crval1, cdelt1, flux, clobber=True):
    """

    Write a simple one-dimensional spectrum.

    Args:
        ofile (str): Name of the file to write.
        crval1 (float): (Log base 10 of the) Initial wavelength, which
            is included in the header with the keyword 'CRVAL1';
            'CRPIX1' is always set to 1.
        cdelt1 (float): The change in (log base 10) wavelength per
            pixel, which is included in the header with the keywords
            'CDELT1' and 'CD1_1'; 'CRPIX1' is always set to 1.
        flux (array): Vector of the flux values.
        clobber (bool): (Optional) Flag to overwrite any existing file
            of the same name.
    """
    hdr = fits.Header()
    hdr['CRPIX1'] = 1
    hdr['CRVAL1'] = crval1
    hdr['CDELT1'] = cdelt1
    hdr['CD1_1'] = hdr['CDELT1']
    fits.HDUList([ fits.PrimaryHDU(flux, header=hdr) ]).writeto(ofile, clobber=clobber)


def read_template_spectrum(filename, data_ext=None, ivar_ext=None, sres_ext=None, log10=False):
    r"""
    Read a template spectrum.
    
    Template spectra are "raw format" files with template data and are,
    at minimum, expected to have the following components::
    
        hdu[0].header['CRVAL1']
        hdu[0].header['CRPIX1']
        hdu[0].header['CDELT1']
        hdu[data_ext].data

    The latter has the flux data.  If `log10` is true, the wavelength
    solution above is expected to be in log wavelengths.
    
    Args:
        filename (str): Name of the fits file to read.

        data_ext (str): (Optional) Name of the extension with the flux
            data.  If None, default is 0.

        ivar_ext (str): (Optional) Name of the extension with the
            inverse variance data.  If None, no inverse data are
            returned.

        sres_ext (str): (Optional) Name of the extension with the
            specral resolution (:math:R=\lambda/\delta\lambda`)
            measurements.  If None, no spectral resolution data are
            returned.

        log10 (bool): (Optional) Flag the WCS wavelength coordinates as
            being in base-10 log wavelength, instead of linear.  Default
            is to assume linear.

    Returns:
        numpy.ndarray : Up to four numpy.float64 arrays with the
            wavelength, flux, inverse variance (if `ivar_ext` is
            provided), and spectral resolution (if `sres_ext` is
            provided) of the template spectrum.

    Raises:
        ValueError: Raised if fits file is not one-dimensional.

        KeyError: Raised if various header keywords or extension names
            are not available.
    """
    if data_ext is None:
        data_ext = 0

    hdu = fits.open(filename, mode='readonly')
    if len(hdu[data_ext].data.shape) != 1:
        raise ValueError('Spectrum in {0} is not one dimensional!'.format(filename))
    
    spec = numpy.copy(hdu[data_ext].data).astype(numpy.float64)
    wave = wavelength_vector(spec.size, hdu[0].header, log10=log10)

    ret = (wave, spec)

    if ivar_ext is not None:
        ret += (numpy.copy(hdu[ivar_ext].data).astype(numpy.float64), )
    if sres_ext is not None:
        ret += (numpy.copy(hdu[sres_ext].data).astype(numpy.float64), )

    hdu.close()

    return ret
    

