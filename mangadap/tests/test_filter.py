from IPython import embed

import numpy
from mangadap.util.filter import BoxcarFilter

def test_boxcar():

    rng = numpy.random.default_rng(1001)
    wave = numpy.arange(3650.0, 10501.0, 0.5)
    x = (wave - numpy.mean(wave))/(wave[-1]-wave[0])
    sig = 1/(0.5 + 0.5*(x+0.5)**2)
    flux = 3 + 0.1*x + 0.5*x**2 + 3*(x+1)**4
    dflux = rng.normal(scale=sig)
    indx = rng.integers(x.size, size=100)
    dflux[indx] += rng.normal(scale=10*sig)[indx]

    indx = ((wave > 4567) & (wave < 4587)) \
                | ((wave > 5567) & (wave < 5587)) \
                | ((wave > 9067) & (wave < 9087))
    sim = numpy.ma.MaskedArray(flux + dflux, mask=indx)

    bf = BoxcarFilter(100, lo=3., hi=3., niter=None, y=sim, local_sigma=True)
    assert numpy.sum(bf.output_mask[0]) == 212, 'Number of rejections changed.'

#    from matplotlib import pyplot
#    pyplot.plot(wave, sim, zorder=1)
#    pyplot.plot(wave, bf.smoothed_y[0], zorder=2)
#    pyplot.scatter(wave[bf.output_mask[0]], sim[bf.output_mask[0]], marker='.', color='C3',
#                   s=100, lw=0, zorder=3)
#    pyplot.xlabel('x')
#    pyplot.ylabel('y')
#    pyplot.show()


