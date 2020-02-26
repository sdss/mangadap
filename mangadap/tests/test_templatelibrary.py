
import pytest
import os

import numpy
from scipy import interpolate
from astropy.io import fits
import astropy.constants

from mangadap.proc.templatelibrary import TemplateLibrary, available_template_libraries
from mangadap.util.sampling import spectrum_velocity_scale, spectral_coordinate_step

import warnings
warnings.simplefilter("ignore", UserWarning)
warnings.simplefilter("ignore", RuntimeWarning)

def test_read():

    tpl_list = available_template_libraries()

    for key in TemplateLibrary.supported_libraries:
        tpl = TemplateLibrary(key, tpllib_list=tpl_list, match_to_drp_resolution=False,
                              velscale_ratio=1, spectral_step=1e-4, log=True, hardcopy=False)


def test_mileshc():
    tpl = TemplateLibrary('MILESHC', match_to_drp_resolution=False, velscale_ratio=4,
                          spectral_step=1e-4, log=True, hardcopy=False)
    # Shape
    assert tpl.ntpl == 42, 'Incorrect number of templates'
    assert tpl['FLUX'].data.shape == (42, 12830), 'Incorrect shape of the flux array'

    # Wavelength and velocity coordinates
    assert numpy.isclose(spectral_coordinate_step(tpl['WAVE'].data, log=True), 1e-4/4), \
                'Coordinate step is wrong'
    assert numpy.isclose(spectrum_velocity_scale(tpl['WAVE'].data),
                         astropy.constants.c.to('km/s').value*1e-4*numpy.log(10.)/4), \
                'Velocity step is wrong'

    # Mask
    nmasked = {}
    for k in tpl.bitmask.keys():
        nmasked[k] = numpy.sum(tpl.bitmask.flagged(tpl['MASK'].data, flag=k))
    assert numpy.sum(list(nmasked.values())) == 8278, 'Masking has changed'
    assert nmasked['NO_DATA'] == 0, 'There should be no NO_DATA flags'

    # Flux values
    flux = numpy.ma.MaskedArray(tpl['FLUX'].data, mask=tpl['MASK'].data > 0)
    assert numpy.isclose(numpy.ma.mean(flux), 1.0), 'Templates should be normalized'

