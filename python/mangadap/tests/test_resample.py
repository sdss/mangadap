
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import pytest

import numpy
from scipy import interpolate

from mangadap.util.instrument import resample1d

def test_resample():
    delt = 100
    x = numpy.arange(20)*0.5 + delt
    y = numpy.linspace(0.1, 2.0, 20)
    y = y*y - y*y*y + 0.2*y + 3 + 0.2 * numpy.square(numpy.sin(50*y))
    newRange = numpy.array([2,8]) + delt
    Xf, Yf = resample1d(y, xRange=[x[0],x[-1]], newRange=newRange, newpix=40, newLog=True)
    assert numpy.isclose(Xf[0], newRange[0]), 'Incorrect starting value'
    assert numpy.mean(interpolate.interp1d(x,y)(Xf) - Yf) < 1e-3, 'Bad interpolation difference'
    X, Y = resample1d(y, xRange=[x[0],x[-1]], newRange=newRange, newpix=40, newLog=True,
                      conserve=True)
    assert numpy.isclose(numpy.mean(Yf/Y), 3.250907334)

