
.. include:: include/links.rst

.. _spatialcovariance:

Spatial Covariance
==================

Spatial covariance is an important consideration when 
propagating the error in aggregated spaxel data in an individual
datacube. These issues are discussed at length in Section 9.3 of `Law
et al. (2016, AJ, 152, 83)`_ and Section 6 of `Westfall et al. (2019,
AJ, 158, 231)`_. See the latter for the specific use of spatial
correlation matrices in the survey-level execution of the DAP.

Spatial covariance matrices for specific wavelength channels have
been calculated by the DRP and are provided in the primary datacube
files. The DAP provides methods that can be used to calculate the
spatial covariance for any wavelength channel or for all of them. The
examples below show how the DAP can be used to read and work with
correlation matrices produced by the DRP and how to produce new
spatial covariance matrices.

Using the DRP-provided correlation matrices
-------------------------------------------

The DRP provides a single correlation matrix at a fiducial wavelength
channel for each of the SDSS *griz* bands. You can use the DAP to
read these correlation matrices as follows:

.. code-block:: python

    from mangadap.config import defaults
    from mangadap.util.covariance import Covariance
    from mangadap.datacube import MaNGADataCube

    plate = 7815
    ifu = 3702
    correl = 'GCORREL'
    path = defaults.drp_directory_path()

    drpfile = os.path.join(path(), MaNGADataCube.build_file_name(plate, ifu))
    
    assert os.path.isfile(drpfile), 'Did not find file'

    with fits.open(drpfile) as hdu:
        covar = Covariance.from_fits(hdu, ivar_ext=None, covar_ext='GCORREL', impose_triu=True,
                                     correlation=True)
        var = numpy.ma.power(hdu['IVAR'].data[hdu['GCORREL'].header['BBINDEX']].T.ravel(),
                             -1).filled(0.0)
        covar = covar.apply_new_variance(var)
        covar.revert_correlation()




    cube = MaNGADataCube.from_plateifu(7815, 3702, directory_path=remote_data_file(),
                                       covar_ext='GCORREL')
    cube.load_rss()

    hdu = fits.open(cube.file_path())
    channel = hdu['GCORREL'].header['BBINDEX']

    gcorrel = numpy.zeros(eval(hdu['GCORREL'].header['COVSHAPE']), dtype=float)
    i = numpy.ravel_multi_index((hdu['GCORREL'].data['INDXI_C1'],hdu['GCORREL'].data['INDXI_C2']),
                                cube.spatial_shape)

    j = numpy.ravel_multi_index((hdu['GCORREL'].data['INDXJ_C1'],hdu['GCORREL'].data['INDXJ_C2']),
                                cube.spatial_shape)
    gcorrel[i,j] = hdu['GCORREL'].data['RHOIJ']
    gcorrel[j,i] = hdu['GCORREL'].data['RHOIJ']

    assert numpy.allclose(cube.covar.toarray(), gcorrel), 'Bad covariance read'

    flux, C = cube.rss.rectify_wavelength_plane(channel, return_covar=True)
    assert numpy.allclose(cube.flux[...,channel], flux), 'Bad flux rectification'

    ivar = numpy.ma.power(C.variance().reshape(cube.spatial_shape), -1).filled(0.0)
    assert numpy.allclose(cube.ivar[...,channel], ivar), 'Bad inverse variance rectification'

    C.to_correlation()
    assert numpy.allclose(C.toarray(), gcorrel), 'Bad covariance calculation'

    sres = numpy.ma.divide(cube.rss.wave[channel],
                           cube.rss.instrumental_dispersion_plane(channel).ravel()) \
                / DAPConstants.sig2fwhm

    # WARNING: The computations done by the DRP and DAP are different
    # in detail, but (at least for this test cube) the results are
    # virtually identical except for notable outliers.
    assert numpy.ma.median(cube.sres[...,channel].ravel() - sres) < 0.1, \
            'Bad spectral resolution rectification'

, and


