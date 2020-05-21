
import pytest

from IPython import embed

import numpy
from astropy.io import fits

from mangadap.proc.reductionassessments import available_reduction_assessments
from mangadap.util.covariance import Covariance
from mangadap.datacube import MaNGADataCube
from mangadap.spectra import MaNGARSS
from mangadap.tests.util import remote_data_file, requires_remote

import warnings
warnings.simplefilter("ignore", UserWarning)
warnings.simplefilter("ignore", RuntimeWarning)


@requires_remote
def test_sres_ext():
    file = remote_data_file(filename=MaNGARSS.build_file_name(7815, 3702, log=True))
    hdu = fits.open(file)
    assert MaNGARSS.spectral_resolution_extension(hdu) == 'LSFPRE', \
                'Bad spectral resolution extension selection'
    assert MaNGARSS.spectral_resolution_extension(hdu, ext='SPECRES') == 'SPECRES', \
                'Bad spectral resolution extension selection'
    assert MaNGARSS.spectral_resolution_extension(hdu, ext='junk') is None, \
                'Should return None for a bad extension name.'


@requires_remote
def test_read():
    rss = MaNGARSS.from_plateifu(7815, 3702, directory_path=remote_data_file())

    assert rss.file_name == MaNGARSS.build_file_name(rss.plate, rss.ifudesign, log=rss.log), \
            'Name mismatch'
    assert rss.log, 'Should read the log-binned version by default.'
    assert len(rss.shape) == 2, 'Row-stacked spectra are 2D'
    assert rss.shape == (rss.nspec, rss.nwave), 'Shape mismatch'
    assert rss.sres is not None, 'Spectral resolution data was not constructed.'
    assert rss.sres_ext == 'LSFPRE', 'Should default to LSFPRE extension.'
    assert rss.xpos.shape == rss.shape, 'On-sky coordinates are wavelength-dependent'
    assert numpy.all(rss.area == numpy.pi), 'Area is pi square arcsec'
    assert rss.area.shape == (rss.nspec,), 'Area is wavelength-independent'

    assert numpy.all(numpy.absolute(numpy.asarray(rss.pointing_offset())) < 0.01), \
                'Pointing offset for this observation should be small.'


@requires_remote
def test_copyto():
    rss = MaNGARSS.from_plateifu(7815, 3702, directory_path=remote_data_file())
    flux = rss.copy_to_array()
    assert not isinstance(flux, numpy.ma.MaskedArray), 'Should output normal array'
    assert flux.shape == rss.shape, 'Both should be 2D arrays.'

    # Apply a wavelength mask
    waverange = [5000, 7000]
    flux = rss.copy_to_array(waverange=waverange)
    indx = (rss.wave > waverange[0]) & (rss.wave < waverange[1])
    assert flux.shape[1] == numpy.sum(indx), 'Wavelength range masking failed'

    # Find the spaxels with non-zero signal
    methods = available_reduction_assessments()
    i = numpy.where([m['key'] == 'SNRG' for m in methods])[0]
    assert len(i) == 1, 'Could not find correct reduction assessment definition.'
    sig, var, snr = rss.flux_stats(response_func=methods[i[0]]['response_func'])
    indx = ((sig > 0) & numpy.invert(numpy.ma.getmaskarray(sig))).data.ravel()
    ngood = numpy.sum(indx)

    # Select the spaxels with non-zero signal
    flux = rss.copy_to_array(waverange=waverange, select_bins=indx)
    assert flux.shape[0] == ngood, 'Bin selection failed'

    # Get the masked array
    flux = rss.copy_to_masked_array()
    assert isinstance(flux, numpy.ma.MaskedArray), 'Should output a masked array'
    assert flux.shape == rss.shape, 'Both should be 2D arrays.'

    # Select the spaxels with non-zero signal
    flux = rss.copy_to_masked_array(select_bins=indx)
    assert flux.shape[0] == ngood, 'Bin selection failed'

    # Try to get the inverse variance
    i = rss.nspec//2
    ivar = rss.copy_to_masked_array(attr='ivar')
    assert ivar.shape == rss.shape, 'Bad ivar shape'
    assert numpy.array_equal(rss.ivar[i], ivar[i].data), 'Did not pull ivar data.'

    # Try to get the spectral resolution
    sres = rss.copy_to_masked_array(attr='sres')
    assert sres.shape == rss.shape, 'Bad sres shape'
    assert numpy.array_equal(rss.sres[i], sres[i].data), 'Did not pull sres data.'


@requires_remote
def test_wcs():
    rss = MaNGARSS.from_plateifu(7815, 3702, directory_path=remote_data_file())
    # Unrestricted
    x, y = rss.mean_sky_coordinates()

    methods = available_reduction_assessments()
    i = numpy.where([m['key'] == 'SNRG' for m in methods])[0]
    assert len(i) == 1, 'Could not find correct reduction assessment definition.'
    # Weighted by the g-band
    _x, _y = rss.mean_sky_coordinates(response_func=methods[i[0]]['response_func'])
    assert numpy.ma.amax(x-_x) - numpy.ma.amin(x-_x) > 0, 'Should be different'
    assert numpy.ma.amax(y-_y) - numpy.ma.amin(y-_y) > 0, 'Should be different'

    # Find a cluster of dithers
    d = numpy.ma.sqrt(numpy.square(_x) + numpy.square(_y))
    srt = numpy.argsort(d)
    theta = numpy.arctan2(-_y, _x)
    indx = theta[srt[:10]] < -2
    assert numpy.sum(indx) == 5, 'Should find close fiber positions'
    indx = srt[:10][indx]
    bin_indx = numpy.full(rss.nspec, -1, dtype=int)
    bin_indx[indx] = 0

    bins, area = rss.binned_on_sky_area(bin_indx, response_func=methods[i[0]]['response_func'])
    assert numpy.array_equal(bins, [0]), 'Should only be one bin'
    try:
        import shapely
    except:
        # Could not import shapely for proper calcultion, so test the stupid calculation
        assert area[0] == numpy.sum(rss.area[bin_indx > -1]), \
                'Stupid calculation is just the sum of the fiber area'
    else:
        # Check the proper calculation
        assert area[0] < numpy.sum(rss.area[bin_indx > -1]), \
                'Area should be substantially smaller than the stupid calculation yields.'

@requires_remote
def test_stats():
    rss = MaNGARSS.from_plateifu(7815, 3702, directory_path=remote_data_file())

    methods = available_reduction_assessments()
    i = numpy.where([m['key'] == 'SNRG' for m in methods])[0]
    assert len(i) == 1, 'Could not find correct reduction assessment definition.'

    cenwave = rss.central_wavelength(response_func=methods[i[0]]['response_func'])
    assert numpy.isclose(cenwave, 4686.2), 'Central wavelength calculation changed'

    sig, var, snr = rss.flux_stats(response_func=methods[i[0]]['response_func'])

    assert sig.shape == (rss.nspec,), 'Should be one measurement per spectrum.'
    assert isinstance(sig, numpy.ma.MaskedArray), 'Expected masked arrays'
    assert numpy.ma.amax(snr) > 15, 'S/N changed'
    assert numpy.ma.median(snr) < 3, 'S/N changed'

    # Try it with the linear rss
    rss = MaNGARSS.from_plateifu(7815, 3702, directory_path=remote_data_file(), log=False)
    _sig, _var, _snr = rss.flux_stats(response_func=methods[i[0]]['response_func'])
    # TODO: Not sure why these are not closer.
    assert numpy.absolute(numpy.ma.median((sig-_sig)/_sig)) < 0.01, \
            'Signal should be the same to better than 1%.'
    assert numpy.absolute(numpy.ma.median((var-_var)/_var)) < 0.03, \
            'Variance should be the same to better than 3%.'
    assert numpy.absolute(numpy.ma.median((snr-_snr)/_snr)) < 0.02, \
            'S/N should be the same to better than 2%.'


@requires_remote
def test_read_lin():
    rss = MaNGARSS.from_plateifu(7815, 3702, directory_path=remote_data_file(), log=False)
    assert not rss.log, 'Wavelength sampling should be linear'
    assert numpy.isclose(numpy.std(numpy.diff(rss.wave)), 0.), \
                'Wavelength sampling should be linear'


@requires_remote
def test_rectification_shape():
    # Load the datacube and the row-stacked spectra
    cube = MaNGADataCube.from_plateifu(7815, 3702, directory_path=remote_data_file())
    cube.load_rss()

    # Get the recitification parameters
    pixelscale, rlim, sigma, recenter, width_buffer \
            = MaNGARSS._parse_rectification_parameters(None, None, None, None, None)
    # Get the cube dimensions
    cube.rss._cube_dimensions(pixelscale=pixelscale, recenter=recenter, width_buffer=width_buffer)
    # Make sure they match what the DRP produced
    assert cube.spatial_shape == (cube.rss.nx, cube.rss.ny), 'Mismatched cube spatial dimensions'


@requires_remote
def test_covariance():
    rss = MaNGARSS.from_plateifu(7815, 3702, directory_path=remote_data_file())

    # Construct a covariance matrix
    C = rss.covariance_matrix(1000)
    assert C.shape == (1764, 1764), 'Bad covariance shape'

    # Make it a correlation matrix and check it
    C.to_correlation()

    # Check that the variances are all unity (or close to it when it's defined)
    unique_var = numpy.unique(numpy.diag(C.toarray()))
    assert numpy.allclose(unique_var[unique_var>0], 1.), 'Bad correlation diagonal'

    # Try multiple channels
    C = rss.covariance_cube(channels=[1000,2000])
    assert numpy.array_equal(C.input_indx, [1000,2000]), 'Bad matrix indices'
    assert C.shape == (1764, 1764, 2), 'Bad covariance shape'

    # Try to convert multiple channels
    C.to_correlation()
    # And reverting it
    C.revert_correlation()

