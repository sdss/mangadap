import pytest

import os
import warnings
warnings.simplefilter("ignore", UserWarning)
warnings.simplefilter("ignore", RuntimeWarning)

from IPython import embed

import numpy

from astropy.io import fits

from matplotlib import pyplot

from mangadap.datacube import MaNGADataCube
from mangadap.util.covariance import Covariance
from mangadap.util.constants import DAPConstants
from mangadap.tests.util import requires_remote, remote_data_file

#-----------------------------------------------------------------------------

def test_samples():

    rng = numpy.random.default_rng(seed=8001)
    
    # Build a bogus covariance matrix
    m = numpy.zeros(10, dtype=float)
    c = numpy.diag(numpy.full(10-2, 0.2, dtype=float), k=-2) \
            + numpy.diag(numpy.full(10-1, 0.5, dtype=float), k=-1) \
            + numpy.diag(numpy.full(10, 1.0, dtype=float), k=0) \
            + numpy.diag(numpy.full(10-1, 0.5, dtype=float), k=1) \
            + numpy.diag(numpy.full(10-2, 0.2, dtype=float), k=2)

    # Draw samples
    s = rng.multivariate_normal(m, c, size=100000)

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


def test_array():
    # Construct the Covariance matrix from a pre-calculated array
    c = numpy.diag(numpy.full(10-2, 0.2, dtype=float), k=-2) \
            + numpy.diag(numpy.full(10-1, 0.5, dtype=float), k=-1) \
            + numpy.diag(numpy.full(10, 1.0, dtype=float), k=0) \
            + numpy.diag(numpy.full(10-1, 0.5, dtype=float), k=1) \
            + numpy.diag(numpy.full(10-2, 0.2, dtype=float), k=2)
    covar = Covariance.from_array(c)
    # Should be the same as the identity matrix.
    assert numpy.array_equal(covar.toarray(), c), 'Arrays should be identical'


def test_io():

    rng = numpy.random.default_rng(seed=8001)
    
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
    s = rng.multivariate_normal(m, c, size=100000)

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


@requires_remote
def test_read_drp():
    drpfile = os.path.join(remote_data_file(), MaNGADataCube.build_file_name(7815, 3702))
    
    assert os.path.isfile(drpfile), 'Did not find file'

    with fits.open(drpfile) as hdu:
        covar = Covariance.from_fits(hdu, ivar_ext=None, covar_ext='GCORREL', impose_triu=True,
                                     correlation=True)
        var = numpy.ma.power(hdu['IVAR'].data[hdu['GCORREL'].header['BBINDEX']].T.ravel(),
                             -1).filled(0.0)

    covar = covar.apply_new_variance(var)
    covar.revert_correlation()

    assert numpy.array_equal(var, numpy.diag(covar.toarray())), 'New variance not applied'


@requires_remote
def test_rectification_recovery():
    cube = MaNGADataCube.from_plateifu(7815, 3702, directory_path=remote_data_file(),
                                       covar_ext='GCORREL')
    cube.load_rss()

    hdu = fits.open(cube.file_path())
    channel = hdu['GCORREL'].header['BBINDEX']

    gcorrel = numpy.zeros(eval(hdu['GCORREL'].header['COVSHAPE']), dtype=float)
    i = numpy.ravel_multi_index((hdu['GCORREL'].data['INDXI_C1'],hdu['GCORREL'].data['INDXI_C2']),
                                cube.spatial_shape)

    j = numpy.ravel_multi_index((hdu['GCORREL'].data['INDXJ_C1'],hdu['GCORREL'].data['INDXJ_C2']),
                                cube.spatial_shape)
    gcorrel[i,j] = hdu['GCORREL'].data['RHOIJ']
    gcorrel[j,i] = hdu['GCORREL'].data['RHOIJ']

    assert numpy.allclose(cube.covar.toarray(), gcorrel), 'Bad covariance read'

    flux, C = cube.rss.rectify_wavelength_plane(channel, return_covar=True)
    assert numpy.allclose(cube.flux[...,channel], flux), 'Bad flux rectification'

    ivar = numpy.ma.power(C.variance().reshape(cube.spatial_shape), -1).filled(0.0)
    assert numpy.allclose(cube.ivar[...,channel], ivar), 'Bad inverse variance rectification'

    C.to_correlation()
    assert numpy.allclose(C.toarray(), gcorrel), 'Bad covariance calculation'

    sres = numpy.ma.divide(cube.rss.wave[channel],
                           cube.rss.instrumental_dispersion_plane(channel).ravel()) \
                / DAPConstants.sig2fwhm

    # WARNING: The computations done by the DRP and DAP are different
    # in detail, but (at least for this test cube) the results are
    # virtually identical except for notable outliers.
    assert numpy.ma.median(cube.sres[...,channel].ravel() - sres) < 0.1, \
            'Bad spectral resolution rectification'

#    zoom = 12
#    xs = int(C.shape[0]/2 - C.shape[0]/2/zoom)
#    xe = xs + int(C.shape[0]/zoom) + 1
#    ys = int(C.shape[1]/2 - C.shape[1]/2/zoom)
#    ye = ys + int(C.shape[1]/zoom) + 1
#
#    fig = pyplot.figure(figsize=(pyplot.figaspect(1)))
#    ax = fig.add_axes([0.1, 0.3, 0.35, 0.35])
#    ax.imshow(C.toarray()[xs:xe,ys:ye], origin='lower', interpolation='nearest',
#              cmap='inferno', aspect='auto')
#    ax.set_title('DAP')
#    ax = fig.add_axes([0.55, 0.3, 0.35, 0.35])
#    cax = fig.add_axes([0.91, 0.3, 0.01, 0.35])
#    ax.set_title('DRP')
#    cs = ax.imshow(cube.covar.toarray()[xs:xe,ys:ye], origin='lower',
#                   interpolation='nearest', aspect='auto', cmap='inferno')
#    pyplot.colorbar(cs, cax=cax)
#    pyplot.show()

