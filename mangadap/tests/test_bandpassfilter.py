import pytest

import numpy
from scipy import stats

from mangadap.proc import bandpassfilter

import warnings
warnings.simplefilter("ignore", UserWarning)
warnings.simplefilter("ignore", RuntimeWarning)

#-----------------------------------------------------------------------------

def test_passband_integral():
    x = numpy.arange(3)
    y = numpy.ones(2)
    Y = bandpassfilter.passband_integral(x,y,borders=True)
    E = bandpassfilter.passband_integral(x,y,borders=True, quad=True)
    assert numpy.all(numpy.isclose([Y,E], [2.0, numpy.sqrt(2)])), 'Bad flat integral'
    Y = bandpassfilter.passband_integral(x,y,borders=True,passband=[0.5,1.5])
    E = bandpassfilter.passband_integral(x,y,borders=True, quad=True,passband=[0.5,1.5])
    assert numpy.all(numpy.isclose([Y,E], [1.0, 1/numpy.sqrt(2)])), 'Bad passband integral'


def test_passband_width():
    x = numpy.arange(5)
    y = numpy.ma.MaskedArray(numpy.ones(4))
    y[2] = numpy.ma.masked
    assert numpy.all(numpy.isclose(bandpassfilter.passband_integrated_width(x,y,borders=True,
                                   passband=numpy.array([[0,3],[0.5,2.5],[3,3.2]])),
                                   numpy.array([2,1.5,0.2]))), 'Incorrect widths'


def test_passband_mean():
    x = numpy.arange(5)
    y = numpy.ma.MaskedArray(numpy.ones(4))
    y[2] = numpy.ma.masked
    assert numpy.all(numpy.isclose(bandpassfilter.passband_integrated_mean(x,y,borders=True,
                                   passband=numpy.array([[0,3],[0.5,2.5],[3,3.2]]))[0],
                                   numpy.ones(3))), 'Incorrect means'


def test_gaussian_integral():
    x = numpy.linspace(-6,6,100)
    g = stats.norm.pdf(x)
    passbands = numpy.array([[-6,6], [-5,5], [-4,4], [-3,3], [-2,2], [-1,1]])
    assert numpy.all(numpy.isclose(bandpassfilter.passband_integral(x,g,passband=passbands),
                            [1., 0.99999942, 0.99993731, 0.99727786, 0.95423538, 0.68231932])), \
            'Incorrect Gaussian integrals'


def test_passband_weighted_stats():
    x = numpy.linspace(-15,15,500)
    g = stats.norm.pdf(x, loc=-9, scale=0.3) + stats.norm.pdf(x, loc=-3, scale=0.6) \
                + stats.norm.pdf(x, loc=8, scale=1.2)
    p = numpy.array([[-11,-7], [-7,1], [1,15]])
    assert numpy.all(numpy.isclose(bandpassfilter.passband_weighted_mean(x,g,x,passband=p)[0],
                                   [-9, -3, 8])), 'Incorrect weighted mean'
    assert numpy.all(numpy.isclose(bandpassfilter.passband_weighted_sdev(x,g,x,passband=p)[0],
                                   [0.3, 0.6, 1.2])), 'Incorrect weighted standard deviation'


#def test_passband_weighted_stats_errors():
#    # Errors are not really robust enough for a useful test.
#    x = numpy.linspace(-15,15,500)
#    g = stats.norm.pdf(x, loc=-9, scale=0.3) + stats.norm.pdf(x, loc=-3, scale=0.6) \
#                + stats.norm.pdf(x, loc=8, scale=1.2)
#    p = numpy.array([[-11,-7], [-7,1], [1,15]])
#
#    nsim = 1000
#    mom0 = numpy.ma.zeros((nsim,3), dtype=float)
#    mom0e = numpy.ma.zeros((nsim,3), dtype=float)
#    mean = numpy.ma.zeros((nsim,3), dtype=float)
#    meane = numpy.ma.zeros((nsim,3), dtype=float)
#    sdev = numpy.ma.zeros((nsim,3), dtype=float)
#    sdeve = numpy.ma.zeros((nsim,3), dtype=float)
#    err = 0.01
#    e = err*numpy.ones_like(x)
#    for i in range(nsim):
#        n = stats.norm.rvs(scale=0.05, size=x.size)
#        mom0[i,:] = bandpassfilter.passband_integral(x,g+n,passband=p)
#        mom0e[i,:] = bandpassfilter.passband_integral(x,g+n,passband=p,quad=True)
#        mean[i,:], meane[i,:] = bandpassfilter.passband_weighted_mean(x,g+n,x,yerr=e,passband=p)
#        sdev[i,:], sdeve[i,:] = bandpassfilter.passband_weighted_sdev(x,g+n,x,yerr=e,passband=p)
#
#    print('Mom0:')
#    print(numpy.mean(mom0, axis=0))
#    print(numpy.std(mom0, axis=0))
#    print(numpy.mean(mom0e, axis=0))
#    print(numpy.std(mom0e, axis=0))
#
#    print('Mean:')
#    print(numpy.mean(mean, axis=0))
#    print(numpy.std(mean, axis=0))
#    print(numpy.mean(meane, axis=0))
#    print(numpy.std(meane, axis=0))
#
#    print('SDEV:')
#    print(numpy.mean(sdev, axis=0))
#    print(numpy.std(sdev, axis=0))
#    print(numpy.mean(sdeve, axis=0))
#    print(numpy.std(sdeve, axis=0))
#
#    from matplotlib import pyplot
#
#    pyplot.hist(mean[:,0].compressed(), range=[-15,15], bins=500)
#    pyplot.hist(mean[:,1].compressed(), range=[-15,15], bins=500)
#    pyplot.hist(mean[:,2].compressed(), range=[-15,15], bins=500)
#    pyplot.show()
#
#    pyplot.hist(sdev[:,0].compressed(), range=[-15,15], bins=500)
#    pyplot.hist(sdev[:,1].compressed(), range=[-15,15], bins=500)
#    pyplot.hist(sdev[:,2].compressed(), range=[-15,15], bins=500)
#    pyplot.show()
    

#if __name__ == '__main__':
#    test_passband_weighted_stats_errors()
#    test_moments()
    
