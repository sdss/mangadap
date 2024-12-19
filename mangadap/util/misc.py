# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""

A catch-all module with miscellaneous utility functions.

----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

def line_coeff(p1, p2):
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


def is_number(s):
    """
    Check if the provided object is a number by trying to convert it to a
    floating-point object.

    Args:
        s (:obj:`object`):
            The object to convert

    Returns:
        :obj:`bool`: True if the object can be converted.
    """
    try:
        float(s)
    except (ValueError, TypeError):
        return False
    else:
        return True


