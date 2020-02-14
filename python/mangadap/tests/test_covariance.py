import pytest

import os
import warnings
warnings.simplefilter("ignore", UserWarning)
warnings.simplefilter("ignore", RuntimeWarning)

import numpy

from astropy.io import fits

from matplotlib import pyplot

from mangadap.util.covariance import Covariance

#-----------------------------------------------------------------------------

def test_samples():
    
    # Build a bogus covariance matrix
    m = numpy.zeros(10, dtype=float)
    c = numpy.diag(numpy.full(10-2, 0.2, dtype=float), k=-2) \
            + numpy.diag(numpy.full(10-1, 0.5, dtype=float), k=-1) \
            + numpy.diag(numpy.full(10, 1.0, dtype=float), k=0) \
            + numpy.diag(numpy.full(10-1, 0.5, dtype=float), k=1) \
            + numpy.diag(numpy.full(10-2, 0.2, dtype=float), k=2)

    # Draw samples
    s = numpy.random.multivariate_normal(m, c, size=100000)

    # Instantiate
    covar = Covariance.from_samples(s.T, cov_tol=0.1)

    # Check the values are very nearly the same as the input
    assert numpy.all(numpy.absolute(c - covar.toarray()) < 0.02), 'Covariances are too different'

    # Check that `find` returns the same indices
    # NOTE: For some reason the ordering of numpy.where and
    # scipy.sparse.find are different.  So I run lexsort before checking
    # if the arrays are the same.

    coo = numpy.array([ [i,j] for i,j in zip(*numpy.where(c > 0))])
    srt = numpy.lexsort((coo[:,0], coo[:,1]))
    coo = coo[srt,:]

    _coo = numpy.array([ [i,j] for i,j,k in zip(*covar.find())])
    srt = numpy.lexsort((_coo[:,0], _coo[:,1]))
    _coo = _coo[srt,:]

    assert numpy.array_equal(coo, _coo), 'Did not find the same indices'


def test_mult():

    # Build a bogus covariance matrix
    c = numpy.diag(numpy.full(10-2, 0.2, dtype=float), k=-2) \
            + numpy.diag(numpy.full(10-1, 0.5, dtype=float), k=-1) \
            + numpy.diag(numpy.full(10, 1.0, dtype=float), k=0) \
            + numpy.diag(numpy.full(10-1, 0.5, dtype=float), k=1) \
            + numpy.diag(numpy.full(10-2, 0.2, dtype=float), k=2)

    x = numpy.ones(10, dtype=float)

    # Uncorrelated
    t = numpy.zeros((3,10), dtype=float)
    t[0,0] = 1.0
    t[1,3] = 1.0
    t[2,6] = 1.0

    y = numpy.dot(t, x)

    covar = Covariance.from_matrix_multiplication(t, c)
    assert numpy.array_equal(covar.toarray(), numpy.identity(3)), 'Result should be uncorrelated.'

    # Correlated by 0.2
    t = numpy.zeros((3,10), dtype=float)
    t[0,0] = 1.0
    t[1,2] = 1.0
    t[2,4] = 1.0
    _c = numpy.diag(numpy.full(3-1, 0.2, dtype=float), k=-1) \
            + numpy.diag(numpy.full(3, 1.0, dtype=float), k=0) \
            + numpy.diag(numpy.full(3-1, 0.2, dtype=float), k=1)
    y = numpy.dot(t, x)
    covar = Covariance.from_matrix_multiplication(t, c)
    assert numpy.array_equal(covar.toarray(), _c), 'Result should have off-diagonals = 0.2'

    # Correlated by 0.5 and 0.2
    t = numpy.zeros((3,10), dtype=float)
    t[0,0] = 1.0
    t[1,1] = 1.0
    t[2,2] = 1.0
    _c = numpy.diag(numpy.full(3-2, 0.2, dtype=float), k=-2) \
            + numpy.diag(numpy.full(3-1, 0.5, dtype=float), k=-1) \
            + numpy.diag(numpy.full(3, 1.0, dtype=float), k=0) \
            + numpy.diag(numpy.full(3-1, 0.5, dtype=float), k=1) \
            + numpy.diag(numpy.full(3-2, 0.2, dtype=float), k=2)
    y = numpy.dot(t, x)
    covar = Covariance.from_matrix_multiplication(t, c)
    assert numpy.array_equal(covar.toarray(), _c), 'Result should have off-diagonals = 0.5,0.2'


def test_var():
    var = numpy.ones(3, dtype=float)
    covar = Covariance.from_variance(var)
    assert numpy.array_equal(covar.toarray(), numpy.identity(3)), \
            'Result should be an identity matrix'

def test_io():

    # Clean up in case of a failure
    ofile = 'test_covar_io.fits'
    if os.path.isfile(ofile):
        os.remove(ofile)

    # Build a bogus covariance matrix
    m = numpy.zeros(10, dtype=float)
    c = numpy.diag(numpy.full(10-2, 0.2, dtype=float), k=-2) \
            + numpy.diag(numpy.full(10-1, 0.5, dtype=float), k=-1) \
            + numpy.diag(numpy.full(10, 1.0, dtype=float), k=0) \
            + numpy.diag(numpy.full(10-1, 0.5, dtype=float), k=1) \
            + numpy.diag(numpy.full(10-2, 0.2, dtype=float), k=2)

    # Draw samples
    s = numpy.random.multivariate_normal(m, c, size=100000)

    # Instantiate
    covar = Covariance.from_samples(s.T, cov_tol=0.1)
    # Write
    covar.write(ofile)
    # Read
    _covar = Covariance.from_fits(ofile)
    # Should be the same
    assert numpy.allclose(covar.toarray(), _covar.toarray()), 'Bad I/O'
    # Clean-up 
    os.remove(ofile)

if __name__ == '__main__':
    test_io()


    

