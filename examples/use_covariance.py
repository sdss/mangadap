#!/usr/bin/env python3

import os
import time
import numpy

from astropy.io import fits
from matplotlib import pyplot

from mangadap.drpfits import DRPFits
from mangadap.util.covariance import Covariance

#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------

if __name__ == '__main__':
    t = time.perf_counter()

    # Read the DRP file
    print('reading')
    drpf = DRPFits(8138, 1901, 'CUBE', directory_path='.', read=True)

    # Get the variance at the wavelength used for the g-band correlation
    # matrix
    bbindex = int(drpf['GCORREL'].header['BBINDEX'])
    var = numpy.ma.power(drpf['IVAR'].data[:,:,bbindex].ravel(), -1).filled(0.0)

    # Use the Covariance class to setup the correlation matrix using the
    # DRP GCORREL extension
    C = Covariance.from_fits(drpf.hdu, ivar_ext=None, covar_ext='GCORREL', correlation=True,
                             impose_triu=True, row_major=True)
    C.show()
    C = C.apply_new_variance(var)
    C.show()
    C.revert_correlation()
    C.show()

    # Recalculate the covariance from scratch
    recalc_C = drpf.covariance_matrix(bbindex)
    recalc_C.show()

    # Calculate the difference between the recalculation at the DRP
    # output
    diff = numpy.ma.MaskedArray(recalc_C.toarray() - C.toarray())
    diff[numpy.invert(numpy.absolute(diff) > 0)] = numpy.ma.masked

    # Absolute difference
    print(numpy.ma.mean(diff))
    print(numpy.ma.std(diff))

    # Fractional difference
    frac = numpy.ma.divide(diff, C.toarray())
    print(numpy.ma.mean(frac))
    print(numpy.ma.std(frac))

#    # Show the difference
#    pyplot.imshow(diff.filled(0.0), origin='lower', interpolation='nearest')
#    pyplot.colorbar()
#    pyplot.show()

    # Flip between the correlation and covariance matrices, and show the
    # difference
    recalc_C.to_correlation()
    recalc_C.show()
    print(recalc_C.is_correlation)
    recalc_C.revert_correlation()
    recalc_C.show()

    print('Elapsed time: {0} seconds'.format(time.perf_counter() - t))



