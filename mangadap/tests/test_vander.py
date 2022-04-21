
import pytest
import os

from IPython import embed

import numpy

from matplotlib import pyplot

from mangadap.util.vander import Legendre1D


def test_1d_fit():
    rng = numpy.random.default_rng(seed=8001)

    x = numpy.linspace(-1, 1, 1000)
    y = 3 + 0.2*x + x*x - x*x*x 
    y += rng.normal(scale=0.01, size=x.size)

    c = numpy.polynomial.legendre.legfit(x,y,4)
    leg = Legendre1D(x,4)
    _c = leg.fit(y)

    assert numpy.allclose(c, _c), 'Legendre1D should match numpy.'


def test_2d_fit():
    rng = numpy.random.default_rng(seed=8001)

    x = numpy.linspace(-1, 1, 1000)
    c = numpy.array([[3, 0.2,   1,   -1],
                     [4, 0.4, 0.8, -0.8]])

    y = numpy.sum([_c[:,None] * numpy.power(x,i)[None,:] for i,_c in enumerate(c.T)], axis=0)
    y += rng.normal(scale=0.01, size=y.size).reshape(y.shape)

    c = numpy.polynomial.legendre.legfit(x,y.T,4)

    leg = Legendre1D(x,4)
    _c = leg.fit(y)

    assert numpy.allclose(c, _c.T), 'Legendre1D should match numpy.'


def test_w1d_1d_fit():

    rng = numpy.random.default_rng(seed=8001)

    x = numpy.linspace(-1, 1, 1000)
    y = 3 + 0.2*x + x*x - x*x*x 
    # Add noise
    y += rng.normal(scale=0.01, size=x.size)
    # Add in deviates
    indx = numpy.unique(rng.integers(1000, size=10))
    y[indx] += rng.normal(scale=10, size=len(indx))
    w = numpy.ones(y.shape, dtype=float)
    w[indx] = 0.

    c_nw = numpy.polynomial.legendre.legfit(x,y,4)
    c = numpy.polynomial.legendre.legfit(x,y,4,w=w)

    leg = Legendre1D(x,4)
    _c_nw = leg.fit(y)
    _c = leg.fit(y,w=w)

    assert numpy.allclose(c_nw, _c_nw), 'Legendre1D fits should match numpy.'
    assert numpy.allclose(c, _c), 'Weighted Legendre1D fits should match numpy.'


def test_w1d_2d_fit():

    rng = numpy.random.default_rng(seed=8001)

    x = numpy.linspace(-1, 1, 1000)
    c = numpy.array([[3, 0.2,   1,   -1],
                     [4, 0.4, 0.8, -0.8]])

    y = numpy.sum([_c[:,None] * numpy.power(x,i)[None,:] for i,_c in enumerate(c.T)], axis=0)
    # Add noise
    y += rng.normal(scale=0.01, size=y.size).reshape(y.shape)
    # Add in deviates
    indx = numpy.unique(rng.integers(1000, size=10))
    y[0,indx] += rng.normal(scale=10, size=len(indx))
    y[1,indx] += rng.normal(scale=10, size=len(indx))
    # Weight all y vectors the same
    w = numpy.ones(x.shape, dtype=float)
    w[indx] = 0.

    c_nw = numpy.polynomial.legendre.legfit(x,y.T,4)
    c = numpy.polynomial.legendre.legfit(x,y.T,4,w=w)

    leg = Legendre1D(x,4)
    _c_nw = leg.fit(y)
    _c = leg.fit(y,w=w)

    assert numpy.allclose(c_nw, _c_nw.T), 'Legendre1D fits should match numpy.'
    assert numpy.allclose(c, _c.T), 'Weighted Legendre1D fits should match numpy.'


def test_w2d_2d_fit():

    rng = numpy.random.default_rng(seed=8001)

    x = numpy.linspace(-1, 1, 1000)
    c = numpy.array([[3, 0.2,   1,   -1],
                     [4, 0.4, 0.8, -0.8]])

    y = numpy.sum([_c[:,None] * numpy.power(x,i)[None,:] for i,_c in enumerate(c.T)], axis=0)
    # Add noise
    y += rng.normal(scale=0.01, size=y.size).reshape(y.shape)
    # Add in deviates
    indx = numpy.unravel_index(numpy.unique(rng.integers(y.size, size=20)), y.shape)
    y[indx] += rng.normal(scale=10, size=len(indx[0]))
    # Weight y vectors independently
    w = numpy.ones(y.shape, dtype=float)
    w[indx] = 0.

    c = numpy.zeros((y.shape[0],5), dtype=float)
    for i in range(y.shape[0]):
        c[i] = numpy.polynomial.legendre.legfit(x,y[i],4,w=w[i])

    leg = Legendre1D(x,4)
    _c = leg.fit(y,w=w)

    assert numpy.allclose(c, _c), 'Weighted Legendre1D fits should match numpy.'

    # Tests both the construction of the model and that the residuals
    # are close to the expectation from the noise added to the data
    # (i.e. is the fit good).
    resid = numpy.ma.std(numpy.ma.MaskedArray(y, mask=numpy.logical_not(w>0)) - leg.model(_c))
    assert numpy.isclose(numpy.round(resid, decimals=2), 0.01), 'Residuals too large'


def test_glbrej_1d_fit():

    rng = numpy.random.default_rng(seed=8001)

    x = numpy.linspace(-1, 1, 1000)
    y = 3 + 0.2*x + x*x - x*x*x 
    # Add noise
    sigma = 0.01
    y += rng.normal(scale=sigma, size=x.size)
    # Add in deviates
    outlier = numpy.unique(rng.integers(1000, size=10))
    offset = rng.normal(scale=10, size=len(outlier))
    indx = numpy.absolute(offset) > 10*sigma
    outlier = outlier[indx]
    y[outlier] += offset[indx] 

    leg = Legendre1D(x,4)
    c, rej = leg.fit(y, rej_iter=-1)

    # Rejected all the outliers?
    assert numpy.all(rej[outlier]), 'Should have rejected all the outliers'


def test_locrej_1d_fit():

    rng = numpy.random.default_rng(seed=8001)

    x = numpy.linspace(-1, 1, 1000)
    y = 3 + 0.2*x + x*x - x*x*x 
    # Add noise
    sigma = 0.01
    y += rng.normal(scale=sigma, size=x.size) * numpy.absolute(x) * 5
    # Add in deviates
    outlier = numpy.unique(rng.integers(1000, size=10))
    offset = rng.normal(scale=10, size=len(outlier))
    indx = numpy.absolute(offset) > 10*sigma
    outlier = outlier[indx]
    y[outlier] += offset[indx] 

    leg = Legendre1D(x,4)
    c, rej = leg.fit(y, rej_iter=-1, rej_win=51)

    # Rejected all the outliers?
    assert numpy.all(rej[outlier]), 'Should have rejected all the outliers'

    # Don't use a local sigma
    _c, _rej = leg.fit(y, rej_iter=-1)

    # Should reject more when not considering local sigma
    assert numpy.sum(_rej) > numpy.sum(rej), 'Should have rejected more when using global sigma'

    # Don't use a local sigma and only one iteration
    _c, __rej = leg.fit(y, rej_iter=1)

    # Should reject more when allowed to continue until no more rejections
    assert numpy.sum(_rej) > numpy.sum(__rej), 'Should have rejected fewer in one iteration'


def test_glbrej_2d_fit():

    rng = numpy.random.default_rng(seed=8001)

    x = numpy.linspace(-1, 1, 1000)
    c = numpy.array([[3, 0.2,   1,   -1],
                     [4, 0.4, 0.8, -0.8]])

    y = numpy.sum([_c[:,None] * numpy.power(x,i)[None,:] for i,_c in enumerate(c.T)], axis=0)
    # Add noise
    sigma = 0.01
    y += rng.normal(scale=sigma, size=y.size).reshape(y.shape)
    # Add in deviates
    outlier = numpy.unravel_index(numpy.unique(rng.integers(y.size, size=20)), y.shape)
    offset = rng.normal(scale=10, size=len(outlier[0]))
    indx = numpy.absolute(offset) > 10*sigma
    outlier = (outlier[0][indx], outlier[1][indx])
    y[outlier] += offset[indx] 

    leg = Legendre1D(x,4)
    c, rej = leg.fit(y, rej_iter=-1)

    # Rejected all the outliers?
    assert numpy.all(rej[outlier]), 'Should have rejected all the outliers'


def test_locrej_2d_fit():

    rng = numpy.random.default_rng(seed=8001)

    x = numpy.linspace(-1, 1, 1000)
    c = numpy.array([[3, 0.2,   1,   -1],
                     [4, 0.4, 0.8, -0.8]])

    y = numpy.sum([_c[:,None] * numpy.power(x,i)[None,:] for i,_c in enumerate(c.T)], axis=0)
    # Add noise
    sigma = 0.01
    y += rng.normal(scale=sigma, size=y.size).reshape(y.shape) \
            * numpy.absolute(x)[None,:] * 5
    # Add in deviates
    outlier = numpy.unravel_index(numpy.unique(rng.integers(y.size, size=20)), y.shape)
    offset = rng.normal(scale=10, size=len(outlier[0]))
    indx = numpy.absolute(offset) > 10*sigma
    outlier = (outlier[0][indx], outlier[1][indx])
    y[outlier] += offset[indx] 

    leg = Legendre1D(x,4)
    c, rej = leg.fit(y, rej_iter=-1, rej_win=51)

    # Rejected all the outliers?
    assert numpy.all(rej[outlier]), 'Should have rejected all the outliers'

    # Don't use a local sigma
    _c, _rej = leg.fit(y, rej_iter=-1)

    # Should reject more when not considering local sigma
    assert numpy.sum(_rej) > numpy.sum(rej), 'Should have rejected more when using global sigma'


def test_ew_glbrej_1d_fit():

    rng = numpy.random.default_rng(seed=8001)

    x = numpy.linspace(-1, 1, 1000)
    y = 3 + 0.2*x + x*x - x*x*x 

    # Add noise
    sigma = 0.01
    noise = rng.normal(scale=sigma, size=x.size)
    err = numpy.full(y.shape, sigma, dtype=float)

    outlier = numpy.unique(rng.integers(1000, size=10))
    offset = rng.normal(scale=10., size=len(outlier))
    noise[outlier] = offset
    err[outlier] = 10.

    # Add noise including larger error
    y += noise

    leg = Legendre1D(x,4)
    c, rej = leg.fit(y, err=err, rej_iter=-1)
    _c, _rej = leg.fit(y, rej_iter=-1)

    # Fewer points should be rejected when accounting for errors
    assert numpy.sum(_rej) > numpy.sum(rej), \
            'Should not have rejected as many points when accounting for error'
