# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
r"""

Defines a class to calculate and convert between on-sky and a projected
(disk) plane.

*Source location*:
    $MANGADAP_DIR/python/mangadap/util/projected_disk_plane.py

*Python2/3 compliance*::

    from __future__ import division
    from __future__ import print_function
    from __future__ import absolute_import
    from __future__ import unicode_literals
    
    import sys
    if sys.version > '3':
        long = int

*Imports*::

    import os.path
    import time
    import numpy
    from scipy import linalg

*Class usage examples*:

    Declare a disk with a position angle of 45 degrees, an inclination
    of 60 degrees, and an on-sky coordinate system that is centered at
    :math:`(x_0,y_0) = (-1.0,0.5)`::

        from mangadap.util.projected_disk_plane import projected_disk_plane
        disk = projected_disk_plane(xc=-1.0, yc=0.5, pa=45, inc=60)

    Determine the disk-plane radius and azimuth for a measurement at a
    focal-plane position of :math:`(x_f,y_f) = (1.0,-0.5)`::

        r, theta = disk.polar(1.0, -0.5)
        print(r, theta)

    Calculate the focal-plane positions for an input
    :math:`(R,\theta)`::

        xf,yf = disk.polar_invert(r, theta)
        print(xf, yf)

*Revision history*:
    | **01 Mar 2015**: Original Implementation by K. Westfall (KBW)
    | **20 May 2015**: (KBW) Documentation and Sphinx tests.  Changed
        :func:`_solve` from hidden to visible.  Moved the check that
        :math:`{\mathbf A}` is define to only occur in :func:`solve`
        since it is called by all the other coordinate calculation
        methods.
    | **04 Jul 2015**: (KBW) Allow for inverse calculations, i.e.
        calculations of the focal-plane cartesian coordinates given the
        disk-plane (Cartesian or polar) coordinates.

.. _scipy.linalg.lu_solve: http://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.lu_solve.html#scipy.linalg.lu_solve
"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import sys
if sys.version > '3':
    long = int

import os.path
import time
import numpy
from scipy import linalg

__author__ = 'Kyle B. Westfall'

class projected_disk_plane:
    r"""

    Calculate projected galaxy coordinates given a set of input
    parameters by calculating :math:`{\mathbf x} = {\mathbf A}^{-1}\
    {\mathbf b}`, where

    .. math::

        {\mathbf A} = \left[
        \begin{array}{rrrrrr}
            1 & 0 & 0 & 0 & 0 & 0 \\
            0 & 1 & 0 & 0 & 0 & 0 \\
            \cos\psi & \sin\psi & -1 & 0 & 0 & 0 \\
            -\sin\psi & \cos\psi & 0 & -1 & 0 &  0 \\
            0 & 0 & \sin\phi_0 & \cos\phi_0 & -1 & 0 \\
            0 & 0 & -\cos\phi_0 & \sin\phi_0 & 0 & -\cos i
        \end{array}
        \right]

        {\mathbf b} = \left[
        \begin{array}{r}
            x_f \\
            y_f \\
            -x_0 \\
            -y_0 \\
            0 \\
            0 
        \end{array}
        \right]

    such that

    .. math::

        {\mathbf x} = \left[
        \begin{array}{r}
            x_f \\
            y_f \\
            x_s \\
            y_s \\
            x_d \\
            y_d 
        \end{array}
        \right]

    and:
        - :math:`\psi` is the Cartesian rotation of the focal-plane
          relative to the sky-plane (+x toward East; +y toward North),
        - :math:`\phi_0` is the on-sky position angle of the major axis
          of the ellipse traced on the sky by a circle in the projected
          plane, defined as the angle from North through East
        - :math:`i` is the inclination of the plane, the angle between
          the plane normal and the normal of the sky plane (:math:`i=0`
          is a face-on plane).
        - :math:`(x_f,y_f)` is the sky-right, focal-plane position
          relative to a reference on-sky position :math:`(x_0,y_0)`
          relative to the center of the projected plane (galaxy center),
        - :math:`(x_s,y_s)` is the on-sky position of :math:`(x_f,y_f)`
          relative to the center of the projected plane, and
        - :math:`(x_d,y_d)` is the Cartesian position of
          :math:`(x_f,y_f)` in the projected plane.

    This form is used such that :math:`{\mathbf A}` need only be defined
    once per class instance.

    The class also allows for inverse calculations, i.e., calculating
    the focal-plane positions provide the disk-plane positions.  In this
    case,

    .. math::

        {\mathbf C} = \left[
        \begin{array}{rrrr}
            \cos\psi & \sin\psi & -1 & 0 \\
            -\sin\psi & \cos\psi & 0 & -1  \\
            0 & 0 & \sin\phi_0 & \cos\phi_0 \\
            0 & 0 & -\cos\phi_0 & \sin\phi_0
        \end{array}
        \right]

        {\mathbf d} = \left[
        \begin{array}{r}
            -x_0 \\
            -y_0 \\
            x_d \\
            y_d \cos i
        \end{array}
        \right]

    such that

    .. math::

        {\mathbf f} = \left[
        \begin{array}{r}
            x_f \\
            y_f \\
            x_s \\
            y_s
        \end{array}
        \right]


    and :math:`{\mathbf f} = {\mathbf C}^{-1}\ {\mathbf d}`.

    Args:
        xc (float): Same as :math:`x_0`, defined above
        yc (float): Same as :math:`y_0`, defined above
        rot (float): Same as :math:`\psi`, defined above
        pa (float): Same as :math:`\phi_0`, defined above
        inc (float): Same as :math:`i`, defined above

    Attributes:
        xc,yc (float,float): a reference on-sky position relative to the
            center of the projected plane (galaxy center); same as
            :math:`(x_0,y_0)` defined above
        rot (float): Cartesian rotation of the focal-plane relative to
            the sky-plane (+x toward East; +y toward North); same as
            :math:`\psi` defined above
        pa (float): On-sky position angle of the major axis of the
            ellipse traced on the sky by a circle in the projected
            plane, defined as the angle from North through East and is
            the same as :math:`\phi_0` defined above
        inc (float): Inclination of the plane, angle between the disk
            normal and the normal of the sky plane (:math:`i=0` is a
            face-on plane).
        A (numpy.ndarray): The coordinate transformation matrix
        Alu (numpy.ndarray): The **lu** array returned by
            `scipy.linalg.lu_factor`_, which is used to calculate the LU
            decomposition of :math:`{\mathbf A}` 
        Apiv (numpy.ndarray): The **piv** array returned by
            `scipy.linalg.lu_factor`_, which is used to calculate the LU
            decomposition of :math:`{\mathbf A}` 
        B (numpy.ndarray): The vector :math:`{\mathbf b}`, as defined
            above, used to calculate :math:`{\mathbf x} = {\mathbf
            A}^{-1}\ {\mathbf b}`
        C (numpy.ndarray): The coordinate transformation matrix use for
            the inverse operations
        Clu (numpy.ndarray): The **lu** array returned by
            `scipy.linalg.lu_factor`_, which is used to calculate the LU
            decomposition of :math:`{\mathbf C}` 
        Cpiv (numpy.ndarray): The **piv** array returned by
            `scipy.linalg.lu_factor`_, which is used to calculate the LU
            decomposition of :math:`{\mathbf C}` 
        D (numpy.ndarray): The vector :math:`{\mathbf d}`, as defined
            above, used to calculate :math:`{\mathbf f} = {\mathbf
            C}^{-1}\ {\mathbf d}`
        
    .. _scipy.linalg.lu_factor: http://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.lu_factor.html#scipy.linalg.lu_factor

    """
    def __init__(self, xc=None, yc=None, rot=None, pa=None, inc=None):

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

        self.C = None
        self.Clu = None
        self.Cpiv = None
        self.D = None

        self._setC()


    def _defined(self):
        """
        Determine if the object is defined such that its methods can be
        used to convert between coordinate systems.
        
        """
        if self.A is None:
            return False
        if self.Alu is None:
            return False
        if self.Apiv is None:
            return False
        if self.C is None:
            return False
        if self.Clu is None:
            return False
        if self.Cpiv is None:
            return False
        return True
    

    def _setA(self):
        """
        Set the transformation matrix and calculate its LU
        decomposition for forward operations.
        """
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
        """Set the on-sky coordinate vector for forward operations."""
        self.B = numpy.array([ x, y, -self.xc, -self.yc, 0.0, 0.0 ])


    def _setC(self):
        """
        Set the transformation matrix and calculate its LU
        decomposition for inverse operations.
        """
        cosr = numpy.cos( numpy.radians(self.rot) )
        sinr = numpy.sin( numpy.radians(self.rot) )
        cosp = numpy.cos( numpy.radians(self.pa) )
        sinp = numpy.sin( numpy.radians(self.pa) )

        self.C = numpy.array([ [  cosr, sinr,  -1.0,  0.0 ],
                               [ -sinr, cosr,   0.0, -1.0 ],
                               [   0.0,  0.0,  sinp, cosp ],
                               [   0.0,  0.0, -cosp, sinp ] ])

        self.Clu, self.Cpiv = linalg.lu_factor(self.C)
        

    def _setD(self, x, y):
        """
        Set the disk-plane coordinate vector for inverse operations.
        """
        cosi = numpy.cos( numpy.radians(self.inc) )
        self.D = numpy.array([ -self.xc, -self.yc, x, cosi*y ])


    def _calculate_polar(self, x, y):
        r"""
        Calculate the polar coordinates (radius and azimuth) provided
        the Cartesian disk-plane coordinates :math:`(x_d,y_d)` using
        
        .. math::

            R &= \sqrt{x_d^2 + y_d^2} \\
            \theta &= \tan^{-1}\left(\frac{-y_d}{x_d}\right)

        Args:
            x,y (float,float): The disk-plane Cartesian coordinates
                :math:`(x_d,y_d)`.

        Returns:
            float, float: The disk-plane polar coordinates: :math:`R,
                \theta`.
        
        """
        R = numpy.sqrt( x*x + y*y)
        theta = numpy.degrees( numpy.arctan2(-y, x) )
        # Returned range in theta is -pi,pi: convert to 0,2pi
        if theta < 0:
            theta += 360

        return R, theta
    

    def _calculate_cartesian(self, r, theta):
        r"""
        Invert the calculation of the projected polar coordinates to
        calculate the Cartesian disk-plane coordinates :math:`(x_d,y_d)`
        using
        
        .. math::

            x_d &= \pm R / \sqrt{1 + \tan^2\theta}\\
            y_d &= -x_d\ \tan\theta

        where :math:`x_d` is negative when :math:`\pi/2 \leq \theta < 3\pi/2`.

        Args:
            r,theta (float,float): The disk-plane polar coordinates
                :math:`(R,\theta)`.

        Returns:
            float, float: The disk-plane Cartesian coordinates:
                :math:`x_d, y_d`.
        
        """
        tant = numpy.tan(numpy.radians(theta))
        xd = r/numpy.sqrt(1.0 + numpy.square(tant))
        if theta > 90 and theta <= 270:
            xd *= -1
        yd = -xd*tant
        return xd, yd
    

    def solve(self, x, y):
        r"""
        Use `scipy.linalg.lu_solve`_ to solve :math:`{\mathbf x} =
        {\mathbf A}^{-1}\ {\mathbf b}`.

        Args:
            x,y (float,float): The coordinate :math:`(x_f,y_f)`, which
                are the sky-right, focal-plane Cartesian coordinates
                relative to a reference on-sky position
                :math:`(x_0,y_0)`, which is relative to the center of
                the projected plane (galaxy center).

        Returns:
            numpy.ndarray : The :math:`{\mathbf x}` vector as defined by
                the solution to :math:`{\mathbf A}^{-1}\ {\mathbf b}`

        Raises:
            Exception: Raised if object was not properly defined.
        """
        if not self._defined():
            raise Exception('projected_coo object not fully defined!')

        self._setB(x,y)
        return linalg.lu_solve((self.Alu, self.Apiv), self.B)


    def solve_inverse(self, x, y):
        r"""
        Use `scipy.linalg.lu_solve`_ to solve :math:`{\mathbf f} =
        {\mathbf C}^{-1}\ {\mathbf d}`.

        Args:
            x,y (float,float): The disk-plane Cartesian coordinates
                :math:`(x_d,y_d)`.

        Returns:
            numpy.ndarray : The :math:`{\mathbf f}` vector as defined by
                the solution to :math:`{\mathbf C}^{-1}\ {\mathbf d}`

        Raises:
            Exception: Raised if object was not properly defined.
        """
        if not self._defined():
            raise Exception('projected_coo object not fully defined!')

        self._setD(x,y)
        return linalg.lu_solve((self.Clu, self.Cpiv), self.D)


    # TODO: Allow these functions to accept array-like x and y?
    def coo(self, x, y):
        r"""
        Calculate :math:`{\mathbf x}` using :func:`solve` for the
        provided :math:`(x_f,y_f)` and return the projected Cartesian
        and polar coordinates, :math:`(x_d,y_d)` and :math:`(R,\theta)`.
        This combines the functionality of :func:`cartesian` and
        :func:`polar`, and so is more efficient than using these both
        separately.

        Args:
            x,y (float,float): The coordinate :math:`(x_f,y_f)`, which
                is the sky-right, focal-plane position relative to a
                reference on-sky position :math:`(x_0,y_0)` relative to
                the center of the projected plane (galaxy center),

        Returns:
            float, float, float, float: The projected Cartesian and
                polar coordinates:  :math:`x_d, y_d, R, \theta`.
        """
        coo = self.solve(x, y)
        R, theta = self._calculate_polar(coo[4], coo[5])
        
        return coo[4], coo[5], R, theta


    def polar(self, x, y):
        r"""
        Calculate :math:`{\mathbf x}` using :func:`solve` for the
        provided :math:`(x_f,y_f)` and return the projected polar
        coordinates, :math:`(R,\theta)`, where
        
        .. math::

            R &= \sqrt{x_d^2 + y_d^2} \\
            \theta &= \tan^{-1}\left(\frac{-y_d}{x_d}\right)

        Args:
            x,y (float,float): The coordinate :math:`(x_f,y_f)`, which
                is the sky-right, focal-plane position relative to a
                reference on-sky position :math:`(x_0,y_0)` relative to
                the center of the projected plane (galaxy center),

        Returns:
            float, float: The projected polar coordinates: :math:`R,
                \theta`.

        """
        coo = self.solve(x, y)
        return self._calculate_polar(coo[4], coo[5])

    
    def polar_invert(self, r, theta):
        r"""
        Calculate :math:`{\mathbf f}` using :func:`solve` for the
        provided :math:`(R,\theta)` and return focal-plane cartesian
        coordinates :math:`(x_f,y_f)`.

        Args:
            r,theta (float,float): The disk-plane polar coordinates
                :math:`(R,\theta)`.

        Returns:
            float, float: The focal-plane Cartesian coordinates
                :math:`(x_f,y_f)`.

        """
        xd, yd = self._calculate_cartesian(r, theta)
        coo = self.solve_inverse(xd, yd)
        return coo[0], coo[1]


    def cartesian(self, x, y):
        r"""
        Calculate :math:`{\mathbf x}` using :func:`solve` for the
        provided :math:`(x_f,y_f)` and return the projected Cartesian
        and coordinates, :math:`(x_d,y_d)`.

        Args:
            x,y (float,float): The coordinate :math:`(x_f,y_f)`, which
                is the sky-right, focal-plane position relative to a
                reference on-sky position :math:`(x_0,y_0)` relative to
                the center of the projected plane (galaxy center),

        Returns:
            float, float: The projected Cartesian coordinates,
                :math:`x_d, y_d`.

        """
        coo = self.solve(x, y)
        return coo[4], coo[5]


    def cartesian_invert(self, x, y):
        r"""
        Calculate :math:`{\mathbf f}` using :func:`solve` for the
        provided :math:`(x_d,y_d)` and return focal-plane cartesian
        coordinates :math:`(x_f,y_f)`.

        Args:
            x,y (float,float): The disk-plane Cartesian coordinates
                :math:`(x_d,y_d)`.

        Returns:
            float, float: The focal-plane Cartesian coordinates
                :math:`(x_f,y_f)`.

        """
        coo = self.solve_inverse(x, y)
        return coo[0], coo[1]


    
