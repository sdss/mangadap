
import pytest
import os

from IPython import embed

import numpy
from astropy.io import fits

from mangadap.par.artifactdb import ArtifactDB
from mangadap.par.emissionlinedb import EmissionLineDB
from mangadap.util.pixelmask import SpectralPixelMask
from mangadap.tests.util import data_test_file

def test_pixelmask():
    specfile = data_test_file('MaNGA_test_spectra.fits.gz')
    hdu = fits.open(specfile)
    pixelmask = SpectralPixelMask(artdb=ArtifactDB.from_key('BADSKY'),
                                  emldb=EmissionLineDB.from_key('ELPSCMSK'))
    assert numpy.sum(pixelmask.boolean(hdu['WAVE'].data, nspec=1)) == 489, \
                'Incorrect number of masked pixels'
