# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""

A catch-all module with miscellaneous utility functions.

*License*:
    Copyright (c) 2015, SDSS-IV/MaNGA Pipeline Group
        Licensed under BSD 3-clause license - see LICENSE.rst

*Source location*:
    $MANGADAP_DIR/python/mangadap/util/misc.py

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

*Revision history*:
    | **2015**: Original implementation by K. Westfall (KBW)
    | **20 May 2015**: (KBW) Documentation and Sphinx tests
    | **04 Jun 2015**: (KBW) Added :func:`where_not`
    | **29 Jul 2016**: (KBW) Change asarray to atleast_1d
    | **17 Nov 2016**: (KBW) Added :func:`high_pass_filter`

.. _numpy.where: http://docs.scipy.org/doc/numpy/reference/generated/numpy.where.html

"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import sys
if sys.version > '3':
    long = int

import numpy

__author__ = 'Kyle Westfall'

def line_coeff(p1, p2):
    # the 'r' here is REQUIRED to get the math below to work
    r"""
    Given two points on a line return the slope and intercept calculated as

    .. math:: 

        m &= \frac{y_1 - y_2}{x_1 - x_2} \\
        b &= y_2 - m\ x_2

    Args:
        p1 (array-like): A two-element :math:`(x_1,y_1)` array with one
            of the two points on the line.
        p2 (array-like): A two-element :math:`(x_2,y_2)` array with one
            of the two points on the line.

    Returns:
        float: The slope (:math:`m`) and intercept (:math:`b`) of the
        line, respectively.

    .. warning:: 
        Performs **no** checks of the input.
    
    """
    m = (p1[1] - p2[1])/(p1[0]-p2[0])
    b = p2[1] - m*p2[0]
    return m, b


def where_not(indx, size):
    """
    Return a tuple with the indices of a vector that were *not* selected
    by a call to `numpy.where`_.  **The function currently only works
    for 1D vectors.**

    Args:
        indx (tuple): Tuple returned by a call to `numpy.where`_ for a
            1D vector.
        size (int): Length of the original vector in the call to
            `numpy.where`_.

    .. warning:: 
        Performs **no** checks of the input.

    """
    return (numpy.setdiff1d(numpy.arange(0,size), indx[0]),)
    

def inverse_with_zeros(v, absolute=True):
    """
    Invert an array with zeros, and handle them by setting the inverse
    to zero.

    Args:
        v (array-like): Array to invert
        absolute (bool): Forces that the absolute value must be larger
            than 0.  Otherwise, any non-positive value is also masked
            with a zero.

    Returns:
        numpy.ndarray: Inverse of the array when the vector has values
        that are greater than 0, otherwise set to 0.0.
    """
    if isinstance(v, numpy.ma.MaskedArray):
        if not absolute:
#            v.mask |= numpy.invert(v > 0)
            v[numpy.invert(v > 0)] = numpy.ma.masked
        return 1.0/v
#    _v = numpy.asarray(v).astype(float)
    _v = numpy.atleast_1d(v).astype(float)
    indx = numpy.absolute(_v) > 0 if absolute else _v > 0
    _v[indx] = 1.0/_v[indx]
    _v[numpy.invert(indx)] = 0.0
    return _v


def high_pass_filter(flux, dw=0, k=None, Dw=None):
    """
    Pulled from FIREFLY's hpf() function and edited on 17 Nov 2016.

    Apply a high-pass filter to the input vector.

    """

    n = flux.size
    _dw = int(dw)

    if k is None and Dw is None:
        Dw = 100
        _k = n//Dw
    elif k is not None and Dw is None:
        _k = int(k)
        Dw = n//_k
    elif k is None and Dw is not None:
        _k = n//Dw
    else:
        raise ValueError('Cannot provide both k and Dw.')

    # PRINT A REPORT

    # Rita's typical inputs for SDSS:
    # w = 10 # 10
    # windowsize = 20 # 20

    # My MaNGA inputs:
    # if w == 0 and windowsize == 0:
    #     print "HPF parameters not specified: choosing default window (stepsize=40)"
    #     w = 40
    #     windowsize = 0

    h           = numpy.fft.fft(flux)
    h_filtered  = numpy.zeros(n, dtype=complex)
    window      = numpy.zeros(n)
    unwindow    = numpy.zeros(n)
    window[0]   = 1                             # keep F0 for normalisation
    unwindow[0] = 1

    for i in range(_dw):
        window[_k+i] = (i+1.0)/_dw
        window[n-1-(_k+_dw-i)] = (_dw-i)/_dw
    window[_k+_dw:n-(_k+_dw)] = 1

    unwindow        = 1 - window
    unwindow[0]     = 1
    
    h_filtered      = h * window
    un_h_filtered   = h * unwindow

    res     = numpy.real(numpy.fft.ifft(h_filtered))
    unres   = numpy.real(numpy.fft.ifft(un_h_filtered)) 
    res_out = (1.0+(res-numpy.median(res))/unres) * numpy.median(res) 

    return res_out, window, res, unres



