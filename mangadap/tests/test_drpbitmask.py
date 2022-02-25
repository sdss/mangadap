import numpy
from astropy.io import fits

from mangadap.datacube import MaNGADataCube
from mangadap.util.drpbitmask import DRPFitsBitMask
from mangadap.tests.util import data_test_file

def test_drpbitmask():
    # Read the data
    specfile = data_test_file('MaNGA_test_spectra.fits.gz')
    hdu = fits.open(specfile)
    drpbm = DRPFitsBitMask()
    assert numpy.sum(drpbm.flagged(hdu['MASK'].data, MaNGADataCube.do_not_fit_flags())) == 4601, \
                'Flags changed'


