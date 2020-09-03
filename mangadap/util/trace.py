# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Provides a class to handle an SDSS trace set; cf. pydl TraceSet.

----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import numpy

class TraceSet:
    """
    Interact with a trace set written by IDLUTILS.

    Args:
        func (:obj:`str`):
        xmin (scalar-like):
        xmax (scalar-like):
        coeff (array-like):

    """
    def __init__(self, func, xmin, xmax, coeff):
        if func not in self.allowed_functions():
            raise ValueError('{0} is not a valid function for use with TraceSet.'.format(func))
        self.func = func
        self.val, self.vander = self._get_val_and_vander()

        self.coeff = numpy.atleast_2d(coeff)
        if self.coeff.ndim != 2:
            raise ValueError('Coefficients array can have no more than two dimensions.')

        # The number of traces and the order of each trace
        self.size, self.order = self.coeff.shape
        self.order -= 1

        self.xmin = xmin
        self.xmax = xmax
        self.Dx = xmax-xmin
        self.xmid = xmin + self.Dx/2

    @classmethod
    def from_fits_table(cls, data):
        return cls(data['FUNC'][0], data['XMIN'][0], data['XMAX'][0], data['COEFF'][0])

    @staticmethod
    def allowed_functions():
        return ['poly', 'legendre', 'chebyshev']

    def _get_val_and_vander(self):
        if self.func == 'poly':
            return numpy.polynomial.polynomial.polyval, numpy.polynomial.polynomial.polyvander
        if self.func == 'legendre':
            return numpy.polynomial.legendre.legval, numpy.polynomial.legendre.legvander
        if self.func == 'chebyshev':
            return numpy.polynomial.chebyshev.chebval, numpy.polynomial.chebyshev.chebvander

        raise ValueError('{0} is not a valid function for use with TraceSet.'.format(self.func))

    def _normalized_coordinates(self, x):
        return 2 * (x-self.xmid)/self.Dx

    def sample(self, x, index=None):
        """
        Sample the trace at the provided coordinates.

        Args:
            x (array-like):
                One-dimensional array with coordinates to sample.

            index (:obj:`int`, optional):
                Sample a single trace with the provided index.  Default
                is to return an array sampling all traces.

        Returns:
            numpy.ndarray: A 1 or 2D array with the samples of one or
            more traces.
        """
        if index is not None:
            return self.val(self._normalized_coordinates(x), self.coeff[index,:])

        return numpy.dot(coeff, self.vander(self._normalized_coordinates(x), self.order))

