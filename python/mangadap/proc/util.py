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
        from astropy import constants

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
from astropy import constants

__author__ = 'Kyle B. Westfall'


def ppxf_fitting_mask(obj_wave, tpl_wave, velScale, velocity_offset, wave_range_analysis=None,
                      max_velocity_range=400., alias_window=2400.):
    """
    Return a list of pixels in the object spectrum to be fit using pPXF.

    The limits applied to the fitted pixels are:

        - Apply the provided wavelength range limit
          (*wave_range_analysis*).

        - pPXF will only allow a fit when the number of template pixels
          is the same as or exceeds the number of pixels in the object
          spectra.  The first step toward limiting template spectra that
          are too long is to truncate the blue and red edges that likely
          won't be used given the provided velocity offsets
          (*velocity_offset*) and the expected velocity range
          (*max_velocity_range*).

        - Remove leading and trailing pixels that will cause alias
          problems during the convolution with the LOSVD
          (*alias_window*).

    Args:
        obj_wave (array): Wavelength vector of the object spectrum to be
            fit.
        tpl_wave (array): Wavelength vector of the template library to
            fit to the object spectrum.
        velScale (float): Velocity scale of the pixel.
        velocity_offset (array) : Vector with the velocity offset
            (expected or actual) between the template and the object
            spectrum in km/s.  Used to estimate which wavelengths can be
            removed from the template.
        wave_range_analysis (array type): (Optional) Lower and upper
            wavelength limits to *include* in the analysis.
        max_velocity_range (float): (Optional) Maximum range (+/-)
            expected for the fitted velocities in km/s.
        alias_window (float) : (Optional) The window to mask to avoid
            aliasing near the edges of the spectral range in km/s.
            Default is six times the default *max_velocity_range*.

    Returns:
        array: Vector with the indices of pixels in the object spectrum
        to fit using pPXF.

    """

    # 1. Apply the wavelength range limit, if provided
    now=len(obj_wave)                               # Number of object wavelengths
    print('Original number of object pixels: {0}'.format(now))
    if wave_range_analysis is not None and len(wave_range_analysis) == 2:
        fit_indx = numpy.where( (obj_wave > wave_range_analysis[0]) \
                                & (obj_wave < wave_range_analysis[1]) )[0]
        if len(fit_indx) == 0:
            raise ValueError('Selected wavelength range for analysis contains no pixels!')
    else:
        fit_indx = numpy.arange(now)

    # Minimum and maximum redshift about primary offsets
    z_min = (numpy.min(velocity_offset) - max_velocity_range)/constants.c.to('km/s').value
    z_max = (numpy.max(velocity_offset) + max_velocity_range)/constants.c.to('km/s').value

    # 2. If the number of template pixels is not >= number of fitted galaxy pixels,
    #    further limit the blue and red edges of the galaxy spectra
    now=len(fit_indx)               # Number of good object pixels
    ntw=len(tpl_wave)               # Number of template pixels
    print('After selecting the wavelength range to analyze: {0}'.format(now))
    print('Number of template pixels: {0}'.format(ntw))
    if ntw < now:
        # Indices of wavelengths redward of the redshifted template
        indx=numpy.where( obj_wave > (tpl_wave[0]*(1. + z_min)))[0]
        if len(indx) == 0:
            raise ValueError('No overlapping wavelengths between galaxy and template!')
        print('Pixels blueward of redshifted template: {0}'.format(len(obj_wave)-len(indx)))

        # Merge with current index
        fit_indx = numpy.intersect1d(fit_indx, indx, assume_unique=True)
        if len(fit_indx) == 0:
            raise ValueError('No overlapping wavelengths between galaxy and template!')
        now_= len(fit_indx)
        print('Object pixels after merging with the specified wavelength range: {0}'.format(now_))

        # New number of good object pixels
        if ntw < now_:
            fit_indx=fit_indx[:ntw]                      # Truncate the red edge as well
            print('Impose same as number of template pixels: {0} {1}'.format(len(fit_indx), ntw))

    # 3. Limit wavelength range to avoid aliasing problems in the template convolution
    nalias = numpy.floor(alias_window/velScale)   # Number of pixels to mask
    # Mask to the range that should be unaffected by alias errors
    wave_range_tpl_unalias = [ tpl_wave[nalias]*(1+z_max), tpl_wave[ntw-nalias-1]*(1+z_min) ]
    print('Mask to these wavelengths to avoid convolution aliasing: {0} - {1}'.format(
          wave_range_tpl_unalias[0], wave_range_tpl_unalias[1]))
    indx = numpy.where( (obj_wave > wave_range_tpl_unalias[0]) \
                        & (obj_wave < wave_range_tpl_unalias[1]) )[0]
    # Merge with current index
    fit_indx = numpy.intersect1d(fit_indx, indx, assume_unique=True)
    if len(fit_indx) == 0:
        raise ValueError('No intersection between wave_range_tpl_unalias and fit_indx!')

    print('Final wavelength range to fit: {0} {1}'.format(obj_wave[fit_indx[0]],
                                                obj_wave[fit_indx[len(fit_indx)-1]]))
    return fit_indx


def ppxf_tpl_obj_voff(tpl_wave, obj_wave, velocity_scale):
    """
    Determine the pseudo offset in velocity between the template and
    object spectra, just due to the difference in the starting
    wavelengths.

    This calculation is independent of the base of the logarithm used
    the sampling of the spectra.

    Args:
        tpl_wave (numpy.ndarray): Wavelength vector for the template
            library to fit to the object spectrum.
        obj_wave (numpy.ndarray): Wavelength vector for the object
            spectrum to be fit.
        velocity_scale (float): Velocity step per pixel in km/s common
            to both the template and object spectrum.

    Returns:
        float: Velocity offset in km/s between the initial wavelengths
        of the template and object spectra.
    """
    return (numpy.log(tpl_wave[0])-numpy.log(obj_wave[0]))*velocity_scale \
                / numpy.diff(numpy.log(obj_wave[0:2]))


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

    for m in method_list:
        print(m['key'])

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
    _v = numpy.full(length, value, dtype=v.dtype)
    diff = list(set( numpy.arange(length) ) - set(missing))
    _v[diff] = v
    return _v


        # TODO: Requires a spectrum in each bin!
#        if any(numpy.diff(bin_indx[srt]) > 1):
#            rep = numpy.ones(bin_change.size, dtype=numpy.int)
#            i = 1
#            while any(numpy.diff(bin_indx[srt]) > i):
#                rep[ numpy.where(numpy.diff(bin_indx[srt]) > i)[0]+1 == bin_change ] += 1
#                i += 1
#            bin_change = numpy.repeat(bin_change, rep)


