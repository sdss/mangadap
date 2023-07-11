import os
from IPython import embed

import numpy
import astropy.constants

from mangadap.datacube import MaNGADataCube
from mangadap.proc.templatelibrary import TemplateLibrary, TemplateLibraryDef
from mangadap.util.sampling import spectrum_velocity_scale, spectral_coordinate_step
from mangadap.tests.util import requires_remote, remote_data_file, data_test_file
from mangadap.config import defaults


def test_read():
    for key in TemplateLibrary.supported_libraries:
        tpl = TemplateLibrary(key, match_resolution=False, velscale_ratio=1, spectral_step=1e-4,
                              log=True, hardcopy=False)

def test_mileshc():
    tpl = TemplateLibrary('MILESHC', match_resolution=False, velscale_ratio=4, spectral_step=1e-4,
                          log=True, hardcopy=False)

    # Shape
    assert tpl.ntpl == 42, 'Incorrect number of templates'
    assert tpl['FLUX'].data.shape == (42, 12634), 'Incorrect shape of the flux array'

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

    assert numpy.sum(list(nmasked.values())) == 192, 'Masking has changed'
    assert nmasked['NO_DATA'] == 0, 'There should be no NO_DATA flags'

    # Flux values
    flux = numpy.ma.MaskedArray(tpl['FLUX'].data, mask=tpl['MASK'].data > 0)
    assert numpy.isclose(numpy.ma.mean(flux), 1.0), 'Templates should be normalized'


def test_write():
    directory, ofile = TemplateLibrary.default_paths('MILESHC', output_path=data_test_file())
    file_name = directory / ofile
    if file_name.exists():
        os.remove(str(file_name))
    tpl = TemplateLibrary('MILESHC', match_resolution=False, velscale_ratio=4, spectral_step=1e-4,
                          log=True, output_path=data_test_file(), hardcopy=True)
    assert file_name.exists(), 'File not written'
    os.remove(str(file_name))


@requires_remote
def test_match_resolution():
    cube = MaNGADataCube.from_plateifu(7815, 3702, directory_path=remote_data_file())
    tpl = TemplateLibrary('MILESHC', cube=cube, match_resolution=True, velscale_ratio=4,
                          hardcopy=False, output_path=remote_data_file())

    # Resolution should be virtually identical in unmasked regions
    indx = tpl['MASK'].data == 0
    assert numpy.std(tpl.sres(tpl['WAVE'].data[indx[0]]) - tpl['SPECRES'].data[0,indx[0]]) < 0.1, \
                'Spectral resolution difference is above tolerance.'

    # Check the file that would have been written has the expected path
    assert cube.directory_path == tpl.directory_path, 'Cube and TPL paths should match.'
    assert tpl.file_name().startswith(cube.output_root), 'TPL file should start with the cube root'


def test_new_library():
    file_search = str(defaults.dap_data_root() / 'spectral_templates' / 'miles'
                        / 'MILES_res2.50_star_m00*.fits')
    tpllib = TemplateLibraryDef('TestLib',
                                file_search=file_search,
                                fwhm=2.5,
                                in_vacuum=False,
                                wave_limit=[3575.,7400.],
                                lower_flux_limit=0.0,
                                log10=False) 
    tpl = TemplateLibrary(tpllib, match_resolution=False, velscale_ratio=4, spectral_step=1e-4,
                          log=True, hardcopy=False)

    assert tpl.ntpl == 99, 'Wrong number of templates'
    assert len(tpl['WAVE'].data) == 12639, 'Wrong spectrum length'

