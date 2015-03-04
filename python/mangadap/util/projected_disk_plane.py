from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import sys
if sys.version > '3':
    long = int

import os.path
from os import environ
import time
import numpy
from scipy import linalg


__author__ = 'Kyle B. Westfall'

class projected_disk_plane:
    """
    Calculate projected galaxy coordinates given a set of input parameters.

    REVISION HISTORY:
        01 Mar 2015: Original Implementation by K. Westfall
    """

    def __init__(self, xc=None, yc=None, rot=None, pa=None, inc=None):
        """
        ARGUMENTS:
            xc - X position of galaxy center. Default is 0.0
            yc - Y position of galaxy center. Default is 0.0
            rot - On-sky rotation of the IFU. Default is 0.0
            pa - On-sky position angle. Default is 0.0
            inc - Plane inclination (0=face-on). Default is 0.0
        """

        self.xc = 0.0 if xc is None else xc
        self.yc = 0.0 if yc is None else yc
        self.rot = 0.0 if rot is None else rot
        self.pa = 0.0 if pa is None else pa
        self.inc = 0.0 if inc is None else inc

        self.A = None
        self.Alu = None
        self.Apiv = None

        self.B = None

        self._setA()

    def _defined(self):
        if self.A is None:
            return False
        if self.Alu is None:
            return False
        if self.Apiv is None:
            return False
        return True
    

    def _setA(self):
        cosr = numpy.cos( numpy.radians(self.rot) )
        sinr = numpy.sin( numpy.radians(self.rot) )
        cosp = numpy.cos( numpy.radians(self.pa) )
        sinp = numpy.sin( numpy.radians(self.pa) )
        cosi = numpy.cos( numpy.radians(self.inc) )

        self.A = numpy.array([ [   1.0,  0.0,   0.0,  0.0,  0.0,   0.0 ],
                               [   0.0,  1.0,   0.0,  0.0,  0.0,   0.0 ],
                               [  cosr, sinr,  -1.0,  0.0,  0.0,   0.0 ],
                               [ -sinr, cosr,   0.0, -1.0,  0.0,   0.0 ],
                               [   0.0,  0.0,  sinp, cosp, -1.0,   0.0 ],
                               [   0.0,  0.0, -cosp, sinp,  0.0, -cosi ] ])

        self.Alu, self.Apiv = linalg.lu_factor(self.A)
        

    def _setB(self, x, y):
        self.B = numpy.array([ x, y, -self.xc, -self.yc, 0.0, 0.0 ])


    def _solve(self, x, y):
        self._setB(x,y)
        return linalg.lu_solve((self.Alu, self.Apiv), self.B)


    def _calculate_polar(self, x, y):

        R = numpy.sqrt( x*x + y*y)
        theta = numpy.degrees( numpy.arctan2(-y, x) )
        # Returned range in theta is -pi,pi: convert to 0,2pi
        if (theta < 0):
            theta += 360

        return R, theta
    

    # TODO: Allow these functions to accept array-like x and y?
    def coo(self, x, y):
        if not self._defined():
            raise Exception('projected_coo object not fully defined!')

        c = self._solve(x, y)
        R, theta = self._calculate_polar(c[4], c[5])
        
        return c[4], c[5], R, theta


    def polar(self, x, y):
        if not self._defined():
            raise Exception('projected_coo object not fully defined!')

        c = self._solve(x, y)
        return self._calculate_polar(c[4], c[5])

    
    def cartesian(self, x, y):
        if not self._defined():
            raise Exception('projected_coo object not fully defined!')

        c = self._solve(x, y)
        return c[4], c[5]

    
