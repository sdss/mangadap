"""

A catch-all module with miscellaneous utility functions.

*Source location*:
    $MANGADAP_DIR/python/mangadap/util/misc.py

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
    | **2015**: Original implementation by K. Westfall (KBW)
    | **20 May 2015**: (KBW) Documentation and Sphinx tests
    | **04 Jun 2015**: (KBW) Added :func:`where_not`


.. _numpy.where: http://docs.scipy.org/doc/numpy/reference/generated/numpy.where.html

"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

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
        real, real: Respectively, the slope (:math:`m`) and intercept
            (:math:`b`) of the line.

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
    


