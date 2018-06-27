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

import warnings
import numpy
import astropy.constants
from scipy import interpolate, spatial

from matplotlib import pyplot

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


def select_proc_method(method_key, method_type, method_list=None, available_func=None, dapsrc=None):
    r"""
    Select a method from a list.  One of method_list or available_func
    must be provided.

    Args:
        method_key (str): Keyword used to select the method.
        method_type (object): Object type to check *method_list*
            against.
        method_list (list): (**Optional**) List of methods from which to
            find the selection keyword.  If none is provided,
            *available_func* **must** be provided.
        available_func (callable): (**Optional**) Callable function that
            returns a list of default methods in place of *method_list*.
            For example, see
            :func:`mangadap.proc.templatelibrary.available_template_libraries`.
        dapsrc (str): (**Optional**) Root path to the DAP source
            directory.  If not provided, the default is defined by
            :func:`mangadap.config.defaults.dap_source_dir`.

    Returns:
        object: An object with base class :class:`mangadap.par.ParSet`,
        containing a set of parameters used to define a method or
        database.

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
    

#def _fill_vector(v, length, missing, value):
#    if v.size == length:
#        v[missing] = value
#        return v
#    _v = numpy.full(length, value, dtype=v.dtype)
#    _v[ list(set( numpy.arange(length) ) - set(missing)) ] = v
#    return _v


def flux_to_fnu(wave, flambda, unit_norm=1e-17):
    r"""

    Convert a spectrum with flux per unit wavelength to flux per unit
    frequency; i.e., calculate
    
    .. math::
        
        F_{\nu} = F_{\lambda} \frac{d\lambda}{d\nu} = F_{\lambda}
        \frac{\lambda^2}{c},

    where the first two arguments of the function are :math:`\lambda`
    and :math:`F_{\lambda}`.  The input wavelength units are expected to
    be angstroms, and the input flux units are expected to be :math:`n\
    {\rm erg\ s}^{-1}\ {\rm cm}^{-2}\ {\rm A}^{-1}`, where :math:`n` is the
    value of *unit_norm*.  The output flux units are microjanskys,
    :math:`10^{-29} {\rm erg\ s}^{-1}\ {\rm cm}^{-2}\ {\rm Hz}^{-1}`.

    Args:
        wave (numpy.ndarray, list): The vector with the wavelengths in
            angstroms.
        flambda (numpy.ndarray, list): The vector with the flux per unit
            wavelength (angstroms).
        unit_norm (float): (**Optional**) The unit normalization of the
            flux.  For example, this is :math:`10^{-17}` when the flux
            units are :math:`10^{-17} {\rm erg\ s}^{-1}\ {\rm cm}^{-2}\
            {\rm A}^{-1}`.

    Returns:
        float,numpy.ndarray: The flux in units of microjanskys.

    Raises:
        ValueError: Raised if the arguments do not have the same shape.
    """
    _wave = [wave] if isinstance(wave, float) else wave
    _wave = numpy.array(_wave) if isinstance(_wave, list) else _wave
    _flambda = [flambda] if isinstance(flambda, float) else flambda
    _flambda = numpy.array(_flambda) if isinstance(_flambda, list) else _flambda
    if _wave.shape != _flambda.shape:
        raise ValueError('Wavelength and flux arrays must have the same shape.')
    fnu = _flambda*numpy.square(_wave)*unit_norm*1e30/astropy.constants.c.to('nm/s').value
    return fnu[0] if isinstance(flambda, float) else fnu


        # TODO: Requires a spectrum in each bin!
#        if any(numpy.diff(bin_indx[srt]) > 1):
#            rep = numpy.ones(bin_change.size, dtype=numpy.int)
#            i = 1
#            while any(numpy.diff(bin_indx[srt]) > i):
#                rep[ numpy.where(numpy.diff(bin_indx[srt]) > i)[0]+1 == bin_change ] += 1
#                i += 1
#            bin_change = numpy.repeat(bin_change, rep)


# OLD VERSION
#def residual_growth(resid, growth_samples):
#    """
#    Interpolate the growth curve at distinct fractions, bracketed by the
#    minimum and maximum.
#    """
#    np = resid.size
#    grw = numpy.arange(np).astype(numpy.float)/np
#    resid_sort = numpy.sort(numpy.absolute(resid))
#    interp = interpolate.interp1d(grw, resid_sort, fill_value='extrapolate')
#    return numpy.append(numpy.append(resid_sort[0], interp(growth_samples)), resid_sort[-1])


#def residual_growth(resid, growth_samples):
#    """
#    Sample a set of residuals at specific the growth intervals.  No
#    interpolation is performed.
#    """
#    np = resid.size
#    resid_sort = numpy.sort(numpy.absolute(resid))
#    i = numpy.zeros(len(growth_samples)+2, dtype=float)
#    i[1:-1] = np*numpy.asarray(growth_samples)
#    i[-1] = np-1
#    i[i < 0] = 0
#    i[i >= np] = np-1
#    return resid_sort[i.astype(int)]


def sample_growth(a, samples, default=-9999., use_interpolate=True):
    _samples = numpy.asarray(samples)
    if numpy.any((_samples < 0) | (_samples > 1)):
        raise ValueError('Growth samples must be between 0 and 1.')
    _a = a.compressed() if isinstance(a, numpy.ma.MaskedArray) else numpy.atleast_1d(a).ravel()
    ns = _samples.size
    if len(_a) < 2:
        return [default]*ns if ns > 1 else default
    srt = numpy.argsort(_a)
    n = _a.size
    grw = (numpy.arange(n,dtype=float)+1)/n
    if use_interpolate:
        interpolator = interpolate.interp1d(grw, _a[srt], fill_value='extrapolate')
        g = interpolator(_samples)
        return tuple(g) if ns > 1 else g
    i = (n*_samples).astype(int)
    i[i > n-1] = n-1
    return _a[srt][i]


def growth_lim(a, lim, fac=1.0, midpoint=None, default=[0., 1.]):
    """
    Set the plots limits of an array based on two growth limits.

    Args:
        a (array-like): Array for which to determine limits.
        lim (float): Percentage of the array values to cover.
        fac (float): (**Optional**) Factor to increase the range based
            on the growth limits.  Default is no increase.
        midpoint (float): (**Optional**) Force the midpoint of the range
            to be centered on this value.  Default is to middle of
            growth range.
        default (list): (**Optional**) Default range to return if `a`
            has no data.  Default is 0 to 1.

    Returns:
        list: Lower and upper limits for the range of a plot of the data
        in `a`.
    """
    # Get the values to plot
    _a = a.compressed() if isinstance(a, numpy.ma.MaskedArray) else numpy.asarray(a).ravel()
    if len(_a) == 0:
        # No data so return the default range
        return default

    # Sort the values
    srt = numpy.ma.argsort(_a)

    # Set the starting and ending values based on a fraction of the
    # growth
    _lim = 1.0 if lim > 1.0 else lim
    start = int(len(_a)*(1.0-_lim)/2)
    end = int(len(_a)*(_lim + (1.0-_lim)/2))
    if end == len(_a):
        end -= 1

    # Set the full range and increase it by the provovided factor
    Da = (_a[srt[end]] - _a[srt[start]])*fac

    # Set the midpoint if not provided
    mid = (_a[srt[start]] + _a[srt[end]])/2 if midpoint is None else midpoint

    # Return the range for the plotted data
    return [ mid - Da/2, mid + Da/2 ]


def optimal_scale(dat, mod, wgt=None):
    r"""
    Calculate the optimal scaling of an input model that minimizes the
    weighted root-mean-square difference between a set of data and a
    model.  When defining the weighted RMS as:

    .. math::
        
        {\rm RMS}^2 = \frac{1}{N} \sum_i w_i^2(d_i - f m_i)^2
    
    The optimal renormalization factor that minimizes the RMS is:

    .. math::
        
        f = \frac{\mathbf{d}^\prime \dot
            \mathbf{m}^\prime}{||\mathbf{m}^\prime||^2}

    where :math:`d^\prime_i = w_i d_i` and :math:`m^\prime_i = w_i m_i`.

    Args:
        dat (array-like): Array of data.
        mod (array-like): Model to renormalize.
        wgt (array-like): (**Optional**) Array of weights to apply to
            each residual.

    Returns:
        float : The optimal scaling that minimizes the weighted
        root-mean-square difference betwen the data and the model.

    Raises:
        ValueError: Raised if the array sizes do not match.
    """
    _dat = numpy.atleast_1d(dat)
    _mod = numpy.atleast_1d(mod)
    _wgt = numpy.ones(_dat.shape, dtype=float) if wgt is None else numpy.atleast_1d(wgt)
    if _mod.shape != _dat.shape or (_wgt is not None and _wgt.shape != _dat.shape):
        raise ValueError('Shapes of all input arrays must match.')

#    print(len(_dat))
#    print(len(_mod))
#    print(len(_wgt))
#    pyplot.scatter(_dat, _mod, marker='.', color='k', lw=0)
#    pyplot.scatter(_dat, _wgt, marker='.', color='r', lw=0)
#    pyplot.show()

    dp = _dat*_wgt
    norm_mp = numpy.sum(numpy.square(_wgt*_mod))
    if norm_mp == 0:
        warnings.warn('Scale not determined because model norm is 0.')

    return 1.0 if norm_mp == 0 else numpy.sum(numpy.square(_wgt)*_dat*_mod)/norm_mp


def replace_with_data_from_nearest_coo(coo, data, replace):
    """
    Replace data in array with the spatially nearest neighbor.

    Args:
        coo (numpy.ndarray): A 2D array with the x and y coordinates of
            all the data.  Shape must be (NDATA, 2).

        data (numpy.ndarray): A 1D or 2D array with data to replace.
            The length of the first (or only) axis must be NDATA.

        replace (numpy.ndarray): Boolean array that is True for elements
            that should be replaced on output.  Shape must be (NDATA,).

    Returns:
            numpy.ndarray: The data array with the selected rows
            replaced with the nearest data set.

    Raises:
        ValueError: Raised if the array sizes are inappropriate.

    """
    # Check the input
    _coo = numpy.asarray(coo)
    if len(_coo.shape) != 2:
        raise ValueError('Input coordinate array must be two-dimensional.')
    if _coo.shape[1] != 2:
        raise ValueError('Currently only works with two coordinates.')
    ndata = _coo.shape[0]
    _data = numpy.asarray(data)
    if _data.shape[0] != ndata:
        raise ValueError('Coordinate and data arrays have a mismatched shape.')
    oned = len(_data.shape) == 1
    _replace = numpy.asarray(replace, dtype=bool)
    if len(_replace.shape) != 1 or _replace.shape[0] != ndata:
        raise ValueError('Input replacement selection array has an incorrect shape.')

    # Nothing flagged to replace, so just return a copy of the input
    if numpy.sum(_replace) == 0:
        return data.copy()

    # Use the coordinates to replace to set the KDTree reference grid
    do_not_replace = numpy.invert(_replace)
    kd = spatial.KDTree(_coo[do_not_replace,:])

    # Get the indices of the nearest data points
    dist, nearest_bin = kd.query(_coo[_replace,:])

    # Replace the existing data with the nearest one and return it
    new_data = _data.copy()
    if oned:
        new_data[_replace] = _data[do_not_replace][nearest_bin]
    else:
        new_data[_replace,:] = _data[do_not_replace,:][nearest_bin,:]
    return new_data


