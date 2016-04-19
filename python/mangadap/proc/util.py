# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Provides a set of processing utility functions for the MaNGA DAP.

*License*:
    Copyright (c) 2015, SDSS-IV/MaNGA Pipeline Group
        Licensed under BSD 3-clause license - see LICENSE.rst

*Source location*:
    $MANGADAP_DIR/python/mangadap/proc/util.py

*Imports and python version compliance*:
    ::

        from __future__ import division
        from __future__ import print_function
        from __future__ import absolute_import
        from __future__ import unicode_literals

        import sys
        if sys.version > '3':
            long = int

        import numpy
        import astropy.constants

*Revision history*:
    | **01 Feb 2016**: Original implementation by K. Westfall (KBW)

.. _astropy.io.fits.hdu.hdulist.HDUList: http://docs.astropy.org/en/v1.0.2/io/fits/api/hdulists.html

"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import sys
if sys.version > '3':
    long = int

import numpy
import astropy.constants

__author__ = 'Kyle B. Westfall'


def HDUList_mask_wavelengths(hdu, bitmask, bitmask_flag, wave_limits, wave_ext='WAVE', \
                             mask_ext='MASK', invert=False):
    """
    Mask pixels in a specified wavelength range by turning on the bit
    value in the specified extention in a provided HDUList object.

    Args:
        hdu (`astropy.io.fits.hdu.hdulist.HDUList`_): HDUList to alter
        bitmask (class:`BitMask`): Bit mask object used to turn on the
            named bit mask.
        bitmask_flag (str): Name of the bit to turn on.
        wave_limits (list or numpy.ndarray): Two-element array with the
            low and high wavelength limits.
        wave_ext (str): (Optional) Name of the wavelength extension in
            *hdu*.
        mask_ext (str): (Optional) Name of the mask extension in *hdu*.
        invert (bool): (Optional) Invert the sense of the masking.
            Instead of masking pixel in the wavelength interval, mask
            all pixels *outside* it.

    Returns:
        `astropy.io.fits.hdu.hdulist.HDUList`_: The modified HDUList
        object.

    Raises:
        Exception: Raised if *wave_limits* does not have a length of
            two.

    """
    if len(wave_limits) != 2:
        raise Exception('Wavelength limits must be a two-element vector.')

    indx = numpy.where( (hdu[wave_ext].data < wave_limits[0]) \
                        | (hdu[wave_ext].data > wave_limits[1])) if invert else \
           numpy.where( (hdu[wave_ext].data >= wave_limits[0]) \
                        & (hdu[wave_ext].data <= wave_limits[1]))
    hdu[mask_ext].data[indx] = bitmask.turn_on(hdu[mask_ext].data[indx], bitmask_flag)
    return hdu


def _select_proc_method(method_key, method_type, method_list=None, available_func=None,
                        dapsrc=None):
    r"""

    Select a method from a list.  One of method_list or available_func
    must be provided.

    Args:
        method_key (str): Keyword used to select the method.
        method_list (list): (Optional) List of methods from which to
            find the selection keyword.  If none is provided,
            *available_func* **must** be provided.
        available_func (callable): (Optional) Callable function that
            returns a list of default methods in place of *method_list.
            For example, see
            :func:`mangadap.proc.templatelibrary.available_template_libraries`.
        dapsrc (str): (Optional) Root path to the DAP source
            directory.  If not provided, the default is defined by
            :func:`mangadap.config.defaults.default_dap_source`.

    Returns:
        :class:`mangadap.par.ParSet`: A set of parameters used to define
        a method or database.

    Raises:
        KeyError: Raised if the selected keyword is not among the
            provided list or if the provided list has more than one
            identical keyword.
        TypeError: Raised if the input *method_list* object is not a
            list or *method_type*, or if available function is not a
            callable function.
    """
    # Get the default methods if no list provided
    if method_list is None:
        if not callable(available_func):
            raise TypeError('If not providing a list, must provide a callable function to ' \
                            'produce the default list of methods/databases/libraries.')
        method_list = available_func(dapsrc=dapsrc)

#    for m in method_list:
#        print(m['key'])

    # Make sure the methods have the right type
    if not isinstance(method_list, list):
        method_list = [method_list]
    for l in method_list:
        if not isinstance(l, method_type):
            raise TypeError('Input method/database/library must have type {0}')

    # Find the selected method via its keyword
    selected_method = [ l['key'] == method_key for l in method_list ]
    if numpy.sum(selected_method) == 0:
        raise KeyError('{0} is not a valid method/database/library!'.format(method_key))
    if numpy.sum(selected_method) > 1:
        raise KeyError('Keywords are not all unique!')

    # Return the method selected via the input keyword
    indx = numpy.where(selected_method)[0][0]
    return method_list[indx]
    

def _fill_vector(v, length, missing, value):
    if v.size == length:
        v[missing] = value
        return v
    _v = numpy.full(length, value, dtype=v.dtype)
    _v[ list(set( numpy.arange(length) ) - set(missing)) ] = v
    return _v


        # TODO: Requires a spectrum in each bin!
#        if any(numpy.diff(bin_indx[srt]) > 1):
#            rep = numpy.ones(bin_change.size, dtype=numpy.int)
#            i = 1
#            while any(numpy.diff(bin_indx[srt]) > i):
#                rep[ numpy.where(numpy.diff(bin_indx[srt]) > i)[0]+1 == bin_change ] += 1
#                i += 1
#            bin_change = numpy.repeat(bin_change, rep)


def _pixel_scale(wave, units='km/s', base=10.0):
    """
    Calculate the velocity scale per pixel.  Assumes the wavelengths
    are logarithmically sampled.
    """
    if units == 'km/s':
        return numpy.diff(numpy.log(wave[0:2]))[0]*astropy.constants.c.to('km/s').value
    if units == 'logw':
        return numpy.diff(numpy.log(wave[0:2]))[0]/numpy.log(base)
    if units == 'ang':
        return numpy.diff(wave[0:2])[0]

    raise ValueError('Units should be km/s, logw, or ang.  Unknown unit: {0}'.format(units))


