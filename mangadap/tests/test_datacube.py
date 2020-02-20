
import pytest

from IPython import embed

import numpy
from astropy.io import fits

from mangadap.util.covariance import Covariance
from mangadap.datacube import MaNGADataCube
from mangadap.tests.util import remote_data_file, requires_remote

import warnings
warnings.simplefilter("ignore", UserWarning)
warnings.simplefilter("ignore", RuntimeWarning)

@requires_remote
def test_sres_ext():
    file = remote_data_file(filename=MaNGADataCube.build_file_name(7815, 3702, log=True))
    hdu = fits.open(file)
    assert MaNGADataCube.spectral_resolution_extension(hdu) == 'DISP', \
                'Bad spectral resolution extension selection'
    assert MaNGADataCube.spectral_resolution_extension(hdu, pre=True) == 'PREDISP', \
                'Bad spectral resolution extension selection'
    assert MaNGADataCube.spectral_resolution_extension(hdu, ext='SPECRES', pre=True) \
                == 'PRESPECRES', 'Bad spectral resolution extension selection'
    assert MaNGADataCube.spectral_resolution_extension(hdu, ext='junk') is None, \
                'Should return None for a bad extension name.'

@requires_remote
def test_read():
    cube = MaNGADataCube.from_plateifu(7815, 3702, directory_path=remote_data_file())

    assert cube.file == MaNGADataCube.build_file_name(cube.plate, cube.ifudesign, log=cube.log), \
            'Name mismatch'
    assert cube.log, 'Should read the log-binned version by default.'
    assert cube.wcs is not None, 'WCS should be defined.'
    assert cube.shape[:2] == cube.spatial_shape, 'Spatial shape should be first two axes.'
    assert cube.nspec == numpy.prod(cube.spatial_shape), 'Definition of number of spectra changed.'
    assert cube.sres is not None, 'Spectral resolution data was not constructed.'
    assert cube.sres_ext == 'DISP' and cube.sres_pre, 'Should default to PREDISP extension.'
    assert abs(cube.pixelscale - cube._get_pixelscale()) < 1e-6, 'Bad match in pixel scale.'
    # NOTE: This is worse than it should be because of how the WCS in MaNGA is defined.
    assert numpy.all(numpy.absolute(cube.wave - cube._get_wavelength_vector(cube.nwave)) < 2e-4), \
            'Bad calculation of wavelength vector.'
    assert cube.covar is None, 'Covariance should not have been read'

@requires_remote
def test_read_correl():
    cube = MaNGADataCube.from_plateifu(7815, 3702, directory_path=remote_data_file(),
                                       covar_ext='GCORREL')
    assert isinstance(cube.covar, Covariance), 'Incorrect type for covariance.'
    assert cube.covar.shape == (cube.nspec,cube.nspec), 'Covariance has incorrect shape.'
    assert cube.covar.is_correlation, 'Covariance object should be in a correlation mode.'
    
    # Check that the variances are all unity (or close to it when it's defined)
    unique_var = numpy.unique(cube.covar.var)
    assert numpy.allclose(unique_var[unique_var>0], 1.), 'Bad variance values'

@requires_remote
def test_wcs():
    cube = MaNGADataCube.from_plateifu(7815, 3702, directory_path=remote_data_file())
    x, y = cube.mean_sky_coordinates()
    assert x[0,0] > x[-1,0], 'RA should increase from large to small indices'
    assert y[0,0] < y[0,-1], 'DEC should increase from small to small indices'
    assert numpy.unravel_index(numpy.argmin( numpy.square(x - cube.prihdr['OBJRA']) 
                                            + numpy.square(y - cube.prihdr['OBJDEC'])), x.shape) \
                == (21,21), 'Object should be at cube center.'
    x, y = cube.mean_sky_coordinates(center_coo=(x[0,0], y[0,0]))
    assert numpy.isclose(x[0,0], 0.0) and numpy.isclose(y[0,0], 0.0), 'Offset incorrect'
    x, y = cube.mean_sky_coordinates(offset='obj')
    assert abs(x[21,21]) < 1e-2 and abs(y[21,21]) < 1e-2, 'Offset incorrect'

@requires_remote
def test_read_lin():
    cube = MaNGADataCube.from_plateifu(7815, 3702, directory_path=remote_data_file(), log=False)
    assert not cube.log, 'Wavelength sampling should be linear'
    assert numpy.isclose(numpy.std(numpy.diff(cube.wave)), 0.), \
                'Wavelength sampling should be linear'

#@requires_remote
def test_load_rss():
    cube = MaNGADataCube.from_plateifu(7815, 3702, directory_path=remote_data_file())
    cube.load_rss()

    embed()
    exit()

if __name__ == '__main__':
    test_load_rss()
