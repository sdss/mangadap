
import pytest

import numpy
from astropy.io import fits
import astropy.constants

from mangadap.drpfits import DRPFits, DRPFitsBitMask
from mangadap.tests.util import data_file

import warnings
warnings.simplefilter("ignore", UserWarning)
warnings.simplefilter("ignore", RuntimeWarning)

def test_drpfitsbitmask():
    # Read the data
    specfile = data_file('MaNGA_test_spectra.fits.gz')
    hdu = fits.open(specfile)
    drpbm = DRPFitsBitMask()
    assert numpy.sum(drpbm.flagged(hdu['MASK'].data, DRPFits.do_not_fit_flags())) == 4601, \
                'Flags changed'

