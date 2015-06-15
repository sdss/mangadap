"""

Provides a set of utility functions dealing with computational geometry.

*Source location*:
    $MANGADAP_DIR/python/mangadap/util/geometry.py

*Python2/3 compliance*::

    from __future__ import division
    from __future__ import print_function
    from __future__ import absolute_import
    
    import sys
    if sys.version > '3':
        long = int

*Imports*:
    None

*Revision history*:
    | **2015**: Original implementation by K. Westfall (KBW)
    | **20 May 2015**: (KBW) Documentation and Sphinx tests
"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import sys
if sys.version > '3':
    long = int

__author__ = 'Kyle Westfall'

def polygon_winding_number(polygon, point):
    """
    Determine the winding number of a 2D polygon about a point.  The
    code does **not** check if the polygon is simple (no interesecting
    line segments).  Algorithm taken from Numerical Recipies Section
    21.4.

    Args:
        polygon (numpy.ndarray): An Nx2 array containing the x,y
            coordinates of a polygon.  The points should be ordered
            either counter-clockwise or clockwise.
        point (numpy.ndarray): A 2-element array defining the x,y
            position of the point to use as a reference for the winding
            number.

    Returns:
        int : Winding number of *polygon* w.r.t. *point*

    Raises:
        Exception: Raised if *polygon* is not 2D, if *polygon* does not
            have two columns, or if *point* is not a 2-element array.
    """

    # Check input shape is for 2D only
    if len(polygon.shape) != 2:
        raise Exception('Polygon must be an Nx2 array.')
    if polygon.shape[1] != 2:
        raise Exception('Polygon must be in two dimensions.')
    if point.size != 2:
        raise Exception('Point must contain two elements.')

    # Get the winding number
    np=polygon.shape[0]
    x0 = polygon[np-1,0]
    y0 = polygon[np-1,1]
    wind = 0
    for i in range(0,np):
        if (y0 > point[1]):
            if polygon[i,1] <= point[1] and \
               (x0-point[0])*(polygon[i,1]-point[1]) - (y0-point[1])*(polygon[i,0]-point[0]) < 0:
                wind -= 1
        else:
            if polygon[i,1] > point[1] and \
               (x0-point[0])*(polygon[i,1]-point[1]) - (y0-point[1])*(polygon[i,0]-point[0]) > 0:
                wind += 1
        x0 = polygon[i,0]
        y0 = polygon[i,1]

    return wind


def point_inside_polygon(polygon, point):
    """

    Determine if a point is inside a polygon using the winding number.

    Args:
        polygon (numpy.ndarray): An Nx2 array containing the x,y
            coordinates of a polygon.  The points should be ordered
            either counter-clockwise or clockwise.
        point (numpy.ndarray): A 2-element array defining the x,y
            position of the point to use as a reference for the winding
            number.

    Returns:
        bool : True if the point is inside the polygon.

    .. warning:: 
        If the point is **on** the polygon (or very close to it w.r.t.
        the machine precision), the returned value is false.

    """
    return (abs(polygon_winding_number(polygon, point)) == 1)
        

