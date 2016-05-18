# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
r"""

Defines a class with common functions for MaNGA fits files.

*License*:
    Copyright (c) 2016, SDSS-IV/MaNGA Pipeline Group
        Licensed under BSD 3-clause license - see LICENSE.rst

*Source location*:
    $MANGADAP_DIR/python/mangadap/mangafits.py

*Imports and python version compliance*:
    ::

        from __future__ import division
        from __future__ import print_function
        from __future__ import absolute_import

        import sys
        if sys.version > '3':
            long = int

        import numpy
        import warnings

        from astropy.io import fits

*Class usage examples*:
    Add some usage comments here!

*Revision history*:
    | **17 May 2015**: Original implementation by K. Westfall (KBW)

.. _astropy.io.fits: http://docs.astropy.org/en/stable/io/fits/index.html
.. _astropy.io.fits.hdu.hdulist.HDUList: http://docs.astropy.org/en/v1.0.2/io/fits/api/hdulists.html
.. _astropy.wcs.wcs.WCS: http://docs.astropy.org/en/v1.0.2/api/astropy.wcs.WCS.html
.. _numpy.meshgrid: http://docs.scipy.org/doc/numpy/reference/generated/numpy.meshgrid.html
.. _scipy.sparse.csr_matrix: http://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csr_matrix.html
.. _numpy.ma.MaskedArray: http://docs.scipy.org/doc/numpy-1.10.1/reference/maskedarray.baseclass.html#numpy.ma.MaskedArray

"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import sys
if sys.version > '3':
    long = int

import numpy
import warnings

from astropy.io import fits

__author__ = 'Kyle B. Westfall'

class MaNGAFits:
    r"""
    An abstract class that only implements a set of static methods for
    interaction with MaNGA fits files.

    .. todo::
        This should eventually be a base class for all the fits files
        used/created by the DAP...

    """
    @staticmethod
    def check_mode(mode):
        """
        Check that the mode is valid.

        Args:

            mode (str): Mode value to check.  Valid modes are `CUBE` and
            `RSS`.
        """
        if mode not in [ 'CUBE', 'RSS']:
            raise ValueError('{0} is not a viable mode.  Must be RSS or CUBE.'.format(mode))
        

    @staticmethod
    def restructure_cube(hdu, ext=['FLUX', 'IVAR', 'MASK'], inverse=False):
        """
        Restructure the cube such that the axes are [x, y, lambda].

        .. warning::
            - Make sure this function remains current with DRP changes.
            - Function will return None if the fits file has not yet
              been read.

        """
        if hdu is None:
            raise ValueError('Input HDUList object is None!')
        if not isinstance(hdu, fits.HDUList):
            raise TypeError('Input must be an astropy.io.fits.HDUList object.')
        _ext = numpy.atleast_1d(ext).ravel()
        for i in range(len(_ext)):
            if len(hdu[_ext[i]].data.shape) != 3:
                raise ValueError('Selected extension is not three-dimensional: {0}'.format(_ext[i]))
        for i in range(len(_ext)):
            hdu[_ext[i]].data = numpy.asarray(hdu[_ext[i]].data.T, order=('F' if inverse else 'C'))
            hdu[_ext[i]].header['NAXIS1'], hdu[_ext[i]].header['NAXIS2'], \
                    hdu[_ext[i]].header['NAXIS3'] = hdu[_ext[i]].data.shape


    @staticmethod
    def restructure_rss(hdu, ext=['FLUX', 'IVAR', 'MASK', 'DISP', 'XPOS', 'YPOS'], inverse=False):
        """
        Restructure the RSS file.  At this point, this only means fixing
        the header.

        .. warning::
            - Make sure this function remains current with DRP changes.
            - Function will return None if the fits file has not yet
              been read.

        """
        if hdu is None:
            raise ValueError('Input HDUList object is None!')
        if not isinstance(hdu, fits.HDUList):
            raise TypeError('Input must be an astropy.io.fits.HDUList object.')
        _ext = numpy.atleast_1d(ext).ravel()
        for i in range(len(_ext)):
            if len(hdu[_ext[i]].data.shape) != 2:
                raise ValueError('Selected extension is not two-dimensional: {0}'.format(_ext[i]))
        for i in range(len(_ext)):
            hdu[_ext[i]].header['NAXIS1'], hdu[_ext[i]].header['NAXIS2'] \
                    = hdu[_ext[i]].data.T.shape if inverse else hdu[_ext[i]].data.shape


    @staticmethod
    def restructure_map(hdu, ext=['BINID'], inverse=False):
        """
        Restructure a mapped quantity such that the axes are [x, y].

        .. warning::
            - Make sure this function remains current with DRP changes.
            - Function will return None if the fits file has not yet
              been read.
        """
        if ext is None:
            raise ValueError('Must provide extensions to restructure.')
        _ext = numpy.atleast_1d(ext).ravel()
        if hdu is None:
            raise ValueError('Input HDUList object is None!')
        if not isinstance(hdu, fits.HDUList):
            raise TypeError('Input must be an astropy.io.fits.HDUList object.')
        for i in range(len(_ext)):
            if len(hdu[_ext[i]].data.shape) != 2:
                raise ValueError('Selected extension is not two-dimensional: {0}'.format(_ext[i]))
        for i in range(len(_ext)):
            hdu[_ext[i]].data = numpy.asarray(hdu[_ext[i]].data.T, order=('F' if inverse else 'C'))
            hdu[_ext[i]].header['NAXIS1'], hdu[_ext[i]].header['NAXIS2'] = hdu[_ext[i]].data.shape


    @staticmethod
    def mode_options():
        """
        Return the allowed modes.

        Returns:
            list: List of the allowed DRP fits file modes.
        """
        return [ 'CUBE', 'RSS' ]


    @staticmethod
    def sampling_options():
        """
        Return the allowed wavelength sampling modes.

        Returns:
            list: List of the allowed DRP fits wavelength sampling
            modes.
        """
        return [ 'LIN', 'LOG' ]


    @staticmethod
    def get_spatial_shape(shape, dispaxis):
        return shape[:dispaxis] + shape[dispaxis+1:]


