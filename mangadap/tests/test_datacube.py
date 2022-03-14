
import pytest

from IPython import embed

import numpy
from astropy.io import fits

from mangadap.proc.reductionassessments import available_reduction_assessments
from mangadap.util.covariance import Covariance
from mangadap.datacube import MaNGADataCube
from mangadap.datacube import MUSEDataCube
from mangadap.tests.util import data_test_file, remote_data_file, requires_remote

import warnings
warnings.simplefilter("ignore", UserWarning)
warnings.simplefilter("ignore", RuntimeWarning)


@requires_remote
def test_sres_ext():
    file = remote_data_file(filename=MaNGADataCube.build_file_name(7815, 3702, log=True))
    hdu = fits.open(file)
    assert MaNGADataCube.spectral_resolution_extension(hdu) == 'LSFPRE', \
                'Bad spectral resolution extension selection'
    assert MaNGADataCube.spectral_resolution_extension(hdu, ext='SPECRES') == 'SPECRES', \
                'Bad spectral resolution extension selection'
    assert MaNGADataCube.spectral_resolution_extension(hdu, ext='junk') is None, \
                'Should return None for a bad extension name.'


@requires_remote
def test_read():
    cube = MaNGADataCube.from_plateifu(7815, 3702, directory_path=remote_data_file())

    assert cube.file_name == MaNGADataCube.build_file_name(cube.plate, cube.ifudesign,
                    log=cube.log), 'Name mismatch'
    assert cube.log, 'Should read the log-binned version by default.'
    assert cube.wcs is not None, 'WCS should be defined.'
    assert cube.shape[:2] == cube.spatial_shape, 'Spatial shape should be first two axes.'
    assert cube.nspec == numpy.prod(cube.spatial_shape), 'Definition of number of spectra changed.'
    assert cube.sres is not None, 'Spectral resolution data was not constructed.'
    assert cube.sres_ext == 'LSFPRE', 'Should default to LSFPRE extension.'
    assert abs(cube.pixelscale - cube._get_pixelscale()) < 1e-6, 'Bad match in pixel scale.'
    # NOTE: This is worse than it should be because of how the WCS in MaNGA is defined.
    assert numpy.all(numpy.absolute(cube.wave - cube._get_wavelength_vector()) < 2e-4), \
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
    x, y = cube.mean_sky_coordinates(offset=None)
    assert x[0,0] > x[-1,0], 'RA should increase from large to small indices'
    assert y[0,0] < y[0,-1], 'DEC should increase from small to small indices'
    assert numpy.unravel_index(numpy.argmin( numpy.square(x - cube.prihdr['OBJRA']) 
                                            + numpy.square(y - cube.prihdr['OBJDEC'])), x.shape) \
                == (21,21), 'Object should be at cube center.'
    x, y = cube.mean_sky_coordinates(center_coo=(x[0,0], y[0,0]))
    assert numpy.isclose(x[0,0], 0.0) and numpy.isclose(y[0,0], 0.0), 'Offset incorrect'
    x, y = cube.mean_sky_coordinates()
    assert abs(x[21,21]) < 1e-2 and abs(y[21,21]) < 1e-2, 'Offset incorrect'


@requires_remote
def test_copyto():
    cube = MaNGADataCube.from_plateifu(7815, 3702, directory_path=remote_data_file())
    flux = cube.copy_to_array()
    assert not isinstance(flux, numpy.ma.MaskedArray), 'Should output normal array'
    assert flux.shape[0] == cube.nspec, 'Should be flattened into a 2D array.'
    assert flux.shape[1] == cube.nwave, 'Should be flattened into a 2D array.'

    # Apply a wavelength mask
    waverange = [5000, 7000]
    flux = cube.copy_to_array(waverange=waverange)
    indx = (cube.wave > waverange[0]) & (cube.wave < waverange[1])
    assert flux.shape[1] == numpy.sum(indx), 'Wavelength range masking failed'

    # Find the spaxels with non-zero signal
    methods = available_reduction_assessments()
    i = numpy.where([m['key'] == 'SNRG' for m in methods])[0]
    assert len(i) == 1, 'Could not find correct reduction assessment definition.'
    sig, var, snr = cube.flux_stats(response_func=methods[i[0]]['response_func'])
    indx = ((sig > 0) & numpy.invert(numpy.ma.getmaskarray(sig))).data.ravel()
    ngood = numpy.sum(indx)

    # Select the spaxels with non-zero signal
    flux = cube.copy_to_array(waverange=waverange, select_bins=indx)
    assert flux.shape[0] == ngood, 'Bin selection failed'

    # Get the masked array
    flux = cube.copy_to_masked_array()
    assert isinstance(flux, numpy.ma.MaskedArray), 'Should output a masked array'
    assert flux.shape[0] == cube.nspec, 'Should be flattened into a 2D array.'
    assert flux.shape[1] == cube.nwave, 'Should be flattened into a 2D array.'

    # Select the spaxels with non-zero signal
    flux = cube.copy_to_masked_array(select_bins=indx)
    assert flux.shape[0] == ngood, 'Bin selection failed'

    # Try to get the inverse variance
    i = cube.nspec//2 + cube.spatial_shape[1]//2
    ivar = cube.copy_to_masked_array(attr='ivar')
    assert ivar.shape == (cube.nspec, cube.nwave), 'Bad ivar shape'
    assert numpy.array_equal(cube.ivar[numpy.unravel_index(i, cube.spatial_shape)],
                             ivar[i].data), 'Did not pull ivar data.'

    # Try to get the spectral resolution
    sres = cube.copy_to_masked_array(attr='sres')
    assert sres.shape == (cube.nspec, cube.nwave), 'Bad sres shape'
    assert numpy.array_equal(cube.sres[numpy.unravel_index(i, cube.spatial_shape)],
                             sres[i].data), 'Did not pull sres data.'


@requires_remote
def test_stats():
    cube = MaNGADataCube.from_plateifu(7815, 3702, directory_path=remote_data_file())

    # Create a fake bin map
    bin_indx = numpy.arange(cube.nspec/4, dtype=int).reshape(cube.spatial_shape[0]//2,
                                                             cube.spatial_shape[0]//2)
    bin_indx = numpy.repeat(bin_indx, 2, axis=0)
    bin_indx = numpy.repeat(bin_indx, 2, axis=1)

    # Get the bin area
    bins, area = cube.binned_on_sky_area(bin_indx)

    assert numpy.array_equal(bins, numpy.arange(cube.nspec/4)), 'Bad bin list'
    assert numpy.allclose(area, 1.), 'Bad area calculation'

    methods = available_reduction_assessments()
    i = numpy.where([m['key'] == 'SNRG' for m in methods])[0]
    assert len(i) == 1, 'Could not find correct reduction assessment definition.'

    cen_wave = cube.central_wavelength(response_func=methods[i[0]]['response_func'],
                                       flag=cube.do_not_use_flags())
    assert numpy.isclose(cen_wave, 4638.0), 'Central wavelength changed.'

    cen_wave = cube.central_wavelength(waverange=[4000,8000], flag=cube.do_not_use_flags(),
                                       fluxwgt=True)
    assert numpy.isclose(cen_wave, 5895.7), 'Central wavelength changed.'

    cen_wave = cube.central_wavelength(waverange=[4000,8000], flag=cube.do_not_use_flags(),
                                       per_pixel=False)
    assert numpy.isclose(cen_wave, 6044.9), 'Central wavelength changed.'

    sig, var, snr = cube.flux_stats(response_func=methods[i[0]]['response_func'])
    assert sig.shape == cube.spatial_shape, 'Should be shaped as a map.'
    assert isinstance(sig, numpy.ma.MaskedArray), 'Expected masked arrays'
    assert numpy.ma.amax(snr) > 60, 'S/N changed'

    # Try it with the linear cube
    cube = MaNGADataCube.from_plateifu(7815, 3702, directory_path=remote_data_file(), log=False)
    _sig, _var, _snr = cube.flux_stats(response_func=methods[i[0]]['response_func'])
    # TODO: Not sure why these are not closer.
    assert numpy.absolute(numpy.ma.median((sig-_sig)/_sig)) < 0.01, \
            'Signal should be the same to better than 1%.'
    assert numpy.absolute(numpy.ma.median((var-_var)/_var)) < 0.03, \
            'Variance should be the same to better than 3%.'
    assert numpy.absolute(numpy.ma.median((snr-_snr)/_snr)) < 0.02, \
            'S/N should be the same to better than 2%.'


@requires_remote
def test_read_lin():
    cube = MaNGADataCube.from_plateifu(7815, 3702, directory_path=remote_data_file(), log=False)
    assert not cube.log, 'Wavelength sampling should be linear'
    assert numpy.isclose(numpy.std(numpy.diff(cube.wave)), 0.), \
                'Wavelength sampling should be linear'


@requires_remote
def test_from_config():
    cube = MaNGADataCube.from_config(data_test_file('datacube.ini'))
    assert cube.meta['z'] == 0.0293823, 'Bad config file read'
    assert cube.meta['ell'] == 0.110844, 'Bad config file read'


@requires_remote
def test_load_rss():
    cube = MaNGADataCube.from_plateifu(7815, 3702, directory_path=remote_data_file())
    cube.load_rss()


@requires_remote
def test_covariance():
    cube = MaNGADataCube.from_plateifu(7815, 3702, directory_path=remote_data_file())

    with pytest.raises(ValueError):
        # Have to load the RSS first
        cube.covariance_matrix(1000)

    # Load the RSS
    cube.load_rss()

    # Construct a covariance matrix
    C = cube.covariance_matrix(1000)
    assert C.shape == (1764, 1764), 'Bad covariance shape'

    # Make it a correlation matrix and check it
    C.to_correlation()

    # Check that the variances are all unity (or close to it when it's defined)
    unique_var = numpy.unique(numpy.diag(C.toarray()))
    assert numpy.allclose(unique_var[unique_var>0], 1.), 'Bad correlation diagonal'

    # Try multiple channels
    C = cube.covariance_cube(channels=[1000,2000])
    assert numpy.array_equal(C.input_indx, [1000,2000]), 'Bad matrix indices'
    assert C.shape == (1764, 1764, 2), 'Bad covariance shape'

    # Try to convert multiple channels
    C.to_correlation()
    # And reverting it
    C.revert_correlation()

    # Try to generate an approximate correlation matrix, covariance
    # matrix, and covariance cube
    approxC = cube.approximate_correlation_matrix()
    approxC = cube.approximate_covariance_matrix(1000)
    approxC = cube.approximate_covariance_cube(channels=[1000,2000])

    # Variance should be the same for direct and approximate calculations
    assert numpy.allclose(approxC.variance(), C.variance()), 'Variances should be the same.'

