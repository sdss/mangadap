
import pytest
import os

import numpy
from astropy.io import fits

from mangadap.par.artifactdb import ArtifactDB
from mangadap.par.emissionlinedb import EmissionLineDB
from mangadap.util.pixelmask import SpectralPixelMask
from mangadap.tests.util import data_file

import warnings
warnings.simplefilter("ignore", UserWarning)
warnings.simplefilter("ignore", RuntimeWarning)

def test_pixelmask():
    specfile = data_file('MaNGA_test_spectra.fits.gz')
    hdu = fits.open(specfile)
    pixelmask = SpectralPixelMask(artdb=ArtifactDB.from_key('BADSKY'),
                                  emldb=EmissionLineDB.from_key('ELPFULL'))
    assert numpy.sum(pixelmask.boolean(hdu['WAVE'].data, nspec=1)) == 487, \
                'Incorrect number of masked pixels'

