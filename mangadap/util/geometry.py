# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""

Provides a set of utility functions dealing with computational geometry.

----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

import numpy
from scipy import linalg

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
        int: Winding number of `polygon` w.r.t. `point`

    Raises:
        ValueError: Raised if `polygon` is not 2D, if `polygon` does not
            have two columns, or if `point` is not a 2-element array.
    """

    # Check input shape is for 2D only
    if len(polygon.shape) != 2:
        raise ValueError('Polygon must be an Nx2 array.')
    if polygon.shape[1] != 2:
        raise ValueError('Polygon must be in two dimensions.')
    if point.size != 2:
        raise ValueError('Point must contain two elements.')

    # Get the winding number
    np=polygon.shape[0]
    x0 = polygon[np-1,0]
    y0 = polygon[np-1,1]
    wind = 0
    for i in range(np):
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
        bool: True if the point is inside the polygon.

    .. warning:: 
        If the point is **on** the polygon (or very close to it w.r.t.
        the machine precision), the returned value is `False`.

    """
    _point = numpy.atleast_2d(point)
    if _point.shape[-1] != 2:
        raise ValueError('Provided point must have two elements in last dimension.')
    if _point.shape[0] == 1:
        return (abs(polygon_winding_number(polygon, point)) == 1)
    return numpy.array([ abs(polygon_winding_number(polygon, p)) == 1 for p in _point ])


def polygon_area(x, y):
    """
    Compute the area of a polygon using the `Shoelace formula
    <https://en.wikipedia.org/wiki/Shoelace_formula>`_.

    Inspired by `this discussion
    <https://stackoverflow.com/questions/24467972/calculate-area-of-polygon-given-x-y-coordinates>`_.

    Args:
        x (`numpy.ndarray`_):
            Vector with the Cartesian x-coordinates of the polygon
            vertices.
        y (`numpy.ndarray`_):
            Vector with the Cartesian y-coordinates of the polygon
            vertices.
            
    Returns:
        :obj:`float`: Polygon area
    """
    return 0.5 * numpy.abs(numpy.dot(x, numpy.roll(y, 1)) - numpy.dot(y, numpy.roll(x, 1)))


class SemiMajorAxisCoo:
    r"""
    Calculate the semi-major axis coordinates given a set of input
    parameters following :math:`{\mathbf x} = {\mathbf A}^{-1}\ {\mathbf
    b}`, where

    .. math::

        {\mathbf A} = \left[
        \begin{array}{rrrrrr}
            1 & 0 & 0 & 0 & 0 & 0 \\
            0 & 1 & 0 & 0 & 0 & 0 \\
            \cos\psi & \sin\psi & -1 & 0 & 0 & 0 \\
            -\sin\psi & \cos\psi & 0 & -1 & 0 &  0 \\
            0 & 0 & \sin\phi_0 & \cos\phi_0 & -1 & 0 \\
            0 & 0 & -\cos\phi_0 & \sin\phi_0 & 0 & \varepsilon-1
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
            x_a \\
            y_a 
        \end{array}
        \right]

    and:
        - :math:`\psi` is the Cartesian rotation of the focal-plane
          relative to the sky-plane (+x toward East; +y toward North),
        - :math:`\phi_0` is the on-sky position angle of the major axis
          of the ellipse, defined as the angle from North through East
        - :math:`\varepsilon=1-b/a` is the ellipticity based on the the
          semi-minor to semi-major axis ratio (:math:`b/a`).
        - :math:`(x_f,y_f)` is the sky-right, focal-plane position
          relative to a reference on-sky position :math:`(x_0,y_0)`
          relative to the center of the ellipse (galaxy center),
        - :math:`(x_s,y_s)` is the on-sky position of :math:`(x_f,y_f)`
          relative to the center of the ellipse, and
        - :math:`(x_a,y_a)` is the Cartesian position of
          :math:`(x_f,y_f)` in units of the semi-major axis.

    This form is used such that :math:`{\mathbf A}` need only be defined
    once per class instance.

    The class also allows for inverse calculations, i.e., calculating
    the focal-plane positions provide the semi-major axis coordinates.
    In this case,

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
            x_a \\
            y_a (1-\varepsilon)
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
        ell (float): Same as :math:`\varepsilon`, defined above

    Attributes:
        xc,yc (float,float): a reference on-sky position relative to the
            center of the ellipse (galaxy center); same as
            :math:`(x_0,y_0)` defined above
        rot (float): Cartesian rotation of the focal-plane relative to
            the sky-plane (+x toward East; +y toward North); same as
            :math:`\psi` defined above
        pa (float): On-sky position angle of the major axis of the
            ellipse, defined as the angle from North through East and is
            the same as :math:`\phi_0` defined above
        ell (float):  Ellipticity define as :math:`\varepsilon=1-b/a`,
            based on the semi-minor to semi-major axis ratio
            (:math:`b/a`) of the ellipse.
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
        
    """
    def __init__(self, xc=None, yc=None, rot=None, pa=None, ell=None):

        self.xc = 0.0 if xc is None else xc
        self.yc = 0.0 if yc is None else yc
        self.rot = 0.0 if rot is None else rot
        self.pa = 0.0 if pa is None else pa
        self.ell = 0.0 if ell is None else ell

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
        #cosi = numpy.cos( numpy.radians(self.inc) )

        self.A = numpy.array([ [   1.0,  0.0,   0.0,  0.0,  0.0,   0.0 ],
                               [   0.0,  1.0,   0.0,  0.0,  0.0,   0.0 ],
                               [  cosr, sinr,  -1.0,  0.0,  0.0,   0.0 ],
                               [ -sinr, cosr,   0.0, -1.0,  0.0,   0.0 ],
                               [   0.0,  0.0,  sinp, cosp, -1.0,   0.0 ],
                               [   0.0,  0.0, -cosp, sinp,  0.0, self.ell-1. ] ])

        self.Alu, self.Apiv = linalg.lu_factor(self.A)
        

    def _setB(self, x, y):
        """
        Set the on-sky coordinate vector for forward operations.

        Args:
            x, y (float) : Single values for use in calculating the
                semi-major-axis coordinates.
        """
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
        Set the semi-major-axis coordinate vector for inverse operations.

        Args:
            x, y (float) : Single values for use in calculating the
                on-sky focal plane coordinates.
        """
        self.D = numpy.array([ -self.xc, -self.yc, x, (1-self.ell)*y ])


    def _calculate_polar(self, x, y):
        r"""
        Calculate the polar coordinates (radius and azimuth) provided
        the Cartesian semi-major-axis coordinates :math:`(x_a,y_a)`
        using
        
        .. math::

            R &= \sqrt{x_a^2 + y_a^2} \\
            \theta &= \tan^{-1}\left(\frac{-y_a}{x_a}\right)

        Args:
            x,y (array-like): The semi-major-axis Cartesian coordinates
                :math:`(x_a,y_a)`.

        Returns:
            numpy.ndarray: The semi-major-axis polar coordinates:
            :math:`R, \theta`.
        """
#        _x = numpy.asarray(x)
#        _y = numpy.asarray(y)
        _x = numpy.atleast_1d(x)
        _y = numpy.atleast_1d(y)
        if _x.size != _y.size:
            raise ValueError('X and Y arrays must have the same size')
        R = numpy.sqrt( _x*_x + _y*_y)
        theta = numpy.degrees( numpy.arctan2(-_y, _x) )
        if hasattr(theta, '__len__'):
            theta[theta < 0] += 360.
        elif theta < 0:
            theta += 360.
        # Returned range in theta is -pi,pi: convert to 0,2pi
        return R, theta
    

    def _calculate_cartesian(self, r, theta):
        r"""
        Invert the calculation of the semi-major-axis polar coordinates
        to calculate the semi-major-axis Cartesian coordinates
        :math:`(x_a,y_a)` using
        
        .. math::

            x_a &= \pm R / \sqrt{1 + \tan^2\theta}\\
            y_a &= -x_a\ \tan\theta

        where :math:`x_a` is negative when :math:`\pi/2 \leq \theta <
        3\pi/2`.

        Args:

            r,theta (array-like): The semi-major-axis polar coordinates
                :math:`(R,\theta)`.

        Returns:
            numpy.ndarray: The semi-major-axis Cartesian coordinates:
            :math:`x_a, y_a`.
        """
#        _r = numpy.asarray(r)
#        _theta = numpy.asarray(theta)
        _r = numpy.atleast_1d(r)
        _theta = numpy.atleast_1d(theta)
        if _r.size != _theta.size:
            raise ValueError('R and THETA arrays must have the same size')
        tant = numpy.tan(numpy.radians(_theta))
        xd = _r/numpy.sqrt(1.0 + numpy.square(tant))
        if hasattr(xd, '__len__'):
            xd[(_theta > 90) & (_theta <= 270)] *= -1.
        elif _theta > 90 and _theta <= 270:
            xd *= -1.
        yd = -xd*tant
        return xd, yd
    

    def solve(self, x, y):
        r"""
        Use `scipy.linalg.lu_solve`_ to solve :math:`{\mathbf x} =
        {\mathbf A}^{-1}\ {\mathbf b}`.

        Args:
            x,y (array-like): The coordinates :math:`(x_f,y_f)`, which
                are the sky-right, focal-plane Cartesian coordinates
                relative to a reference on-sky position
                :math:`(x_0,y_0)`, which is relative to the center of
                the ellipse (galaxy center).

        Returns:
            numpy.ndarray: The :math:`{\mathbf x}` vectors (separated by
            rows) as defined by the solution to :math:`{\mathbf A}^{-1}\
            {\mathbf b}`

        Raises:
            ValueError: Raised if object was not properly defined or if
                the X and Y arrays do not have the same size.
        """
        if not self._defined():
            raise ValueError('SemiMajorAxisCoo object not fully defined!')
#        _x = numpy.asarray(x)
#        _y = numpy.asarray(y)
        _x = numpy.atleast_1d(x)
        _y = numpy.atleast_1d(y)
        if _x.size != _y.size:
            raise ValueError('X and Y arrays must have the same size')
        out = numpy.empty((_x.size, 6), dtype=float)
        for i, (xx, yy) in enumerate(zip(_x.ravel(),_y.ravel())):
            self._setB(xx,yy)
            out[i,:] = linalg.lu_solve((self.Alu, self.Apiv), self.B)
        return out.reshape(_x.shape + (6,))


    def solve_inverse(self, x, y):
        r"""
        Use `scipy.linalg.lu_solve`_ to solve :math:`{\mathbf f} =
        {\mathbf C}^{-1}\ {\mathbf d}`.

        Args:
            x,y (array-like): The semi-major-axis Cartesian coordinates
                :math:`(x_a,y_a)`.

        Returns:
            numpy.ndarray: The :math:`{\mathbf f}` vector as defined by
            the solution to :math:`{\mathbf C}^{-1}\ {\mathbf d}`

        Raises:
            ValueError: Raised if object was not properly defined or if
                the X and Y arrays do not have the same size.
        """
        if not self._defined():
            raise ValueError('SemiMajorAxisCoo object not fully defined!')
#        _x = numpy.asarray(x)
#        _y = numpy.asarray(y)
        _x = numpy.atleast_1d(x)
        _y = numpy.atleast_1d(y)
        if _x.size != _y.size:
            raise ValueError('X and Y arrays must have the same size')
        out = numpy.empty((_x.size, 4), dtype=float)
        for i, (xx, yy) in enumerate(zip(_x.ravel(),_y.ravel())):
            self._setD(xx,yy)
            out[i,:] = linalg.lu_solve((self.Clu, self.Cpiv), self.D)
        return out.reshape(_x.shape + (4,))


    def coo(self, x, y):
        r"""
        Calculate :math:`{\mathbf x}` using :func:`solve` for the
        provided :math:`(x_f,y_f)` and return the semi-major-axis
        Cartesian and polar coordinates, :math:`(x_a,y_a)` and
        :math:`(R,\theta)`.  This combines the functionality of
        :func:`cartesian` and :func:`polar`, and so is more efficient
        than using these both separately.

        Args:
            x,y (array-like): The coordinates :math:`(x_f,y_f)`, which
                are the sky-right, focal-plane position relative to a
                reference on-sky position :math:`(x_0,y_0)` relative to
                the center of the ellipse (galaxy center),

        Returns:
            numpy.ndarray: Four arrays with the semi-major-axis
            Cartesian and polar coordinates:  :math:`x_a, y_a, R,
            \theta`.

        """
        coo = self.solve(x, y)
        return (coo[...,4], coo[...,5]) + self._calculate_polar(coo[...,4], coo[...,5])


    def polar(self, x, y):
        r"""
        Calculate :math:`{\mathbf x}` using :func:`solve` for the
        provided :math:`(x_f,y_f)` and return the semi-major-axis polar
        coordinates, :math:`(R,\theta)`, where
        
        .. math::

            R &= \sqrt{x_a^2 + y_a^2} \\
            \theta &= \tan^{-1}\left(\frac{-y_a}{x_a}\right)

        Args:
            x,y (array-like): The coordinate :math:`(x_f,y_f)`, which
                is the sky-right, focal-plane position relative to a
                reference on-sky position :math:`(x_0,y_0)` relative to
                the center of the ellipse (galaxy center),

        Returns:
            numpy.ndarray: Two arrays with the semi-major-axis polar
            coordinates: :math:`R, \theta`.
        """
        coo = self.solve(x, y)
        return self._calculate_polar(coo[...,4], coo[...,5])

    
    def polar_invert(self, r, theta):
        r"""
        Calculate :math:`{\mathbf f}` using :func:`solve` for the
        provided :math:`(R,\theta)` and return focal-plane cartesian
        coordinates :math:`(x_f,y_f)`.

        Args:
            r,theta (array-like): The semi-major-axis polar coordinates
                :math:`(R,\theta)`.

        Returns:
            numpy.ndarray: Two arrays with the focal-plane Cartesian
            coordinates :math:`(x_f,y_f)`.
        """
        xd, yd = self._calculate_cartesian(r, theta)
        coo = self.solve_inverse(xd, yd)
        return coo[...,0], coo[...,1]


    def cartesian(self, x, y):
        r"""
        Calculate :math:`{\mathbf x}` using :func:`solve` for the
        provided :math:`(x_f,y_f)` and return the semi-major-axis
        Cartesian and coordinates, :math:`(x_a,y_a)`.

        Args:
            x,y (array-like): The coordinate :math:`(x_f,y_f)`, which
                is the sky-right, focal-plane position relative to a
                reference on-sky position :math:`(x_0,y_0)` relative to
                the center of the ellipse (galaxy center),

        Returns:
            numpy.ndarray: Two arrays with the semi-major-axis Cartesian
            coordinates, :math:`x_a, y_a`.
        """
        coo = self.solve(x, y)
        return coo[...,4], coo[...,5]


    def cartesian_invert(self, x, y):
        r"""
        Calculate :math:`{\mathbf f}` using :func:`solve` for the
        provided :math:`(x_a,y_a)` and return focal-plane cartesian
        coordinates :math:`(x_f,y_f)`.

        Args:
            x,y (array-like): The semi-major-axis Cartesian coordinates
                :math:`(x_a,y_a)`.

        Returns:
            numpy.ndarray: The focal-plane Cartesian coordinates
            :math:`(x_f,y_f)`.
        """
        coo = self.solve_inverse(x, y)
        return coo[...,0], coo[...,1]


    



