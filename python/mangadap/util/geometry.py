from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import sys
if sys.version > '3':
    long = int

import inspect

__author__ = 'Kyle Westfall'

def polygon_winding_number(polygon, point):
    """
    Determine the winding number of a polygon in 2D about a point.
    Taken from Numerical Recipies Section 21.4.
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
    Determine if a point is inside a polygon.

    If the point is *on* the polygon (or very close to it wrt the
    machine precision), the returned value is false.
    """
    return (abs(polygon_winding_number(polygon, point)) == 1)
        

