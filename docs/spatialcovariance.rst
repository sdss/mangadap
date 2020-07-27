
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

    # Imports
    from mangadap.config import defaults
    from mangadap.util.covariance import Covariance
    from mangadap.datacube import MaNGADataCube

    # Set the observation and correlation matrix extension
    plate = 7815
    ifu = 3702
    correl = 'GCORREL'
    path = defaults.drp_directory_path()

    # Set the file and open it
    drpfile = os.path.join(path, MaNGADataCube.build_file_name(plate, ifu))
    with fits.open(drpfile) as hdu:
        # Read the correlation data; setting ivar_ext to None here is
        # important.
        covar = Covariance.from_fits(hdu, ivar_ext=None, covar_ext=correl, impose_triu=True,
                                     correlation=True)
        # Calculate the variance in the wavelength channel used to
        # construct the correlation matrix.
        var = numpy.ma.power(hdu['IVAR'].data[hdu[correl].header['BBINDEX']].T.ravel(),
                             -1).filled(0.0)

    # Apply the variance to the correlation matrix and revert back to a
    # covariance matrix
    covar = covar.apply_new_variance(var)
    covar.revert_correlation()

    # Plot the covariance matrix
    covar.show()

You can also read the DRP-produce correlation matrix when
instantiating the MaNGA datacube. However, note that this assumes the
correlation matrix is applicable to *all* wavelength channels:

.. code-block:: python

    # Imports
    from mangadap.config import defaults
    from mangadap.datacube import MaNGADataCube

    # Read the MaNGA datacube
    path = defaults.drp_directory_path()
    cube = MaNGADataCube.from_plateifu(7815, 3702, directory_path=path, covar_ext='GCORREL')

    # Show the correlation matrix
    cube.covar.show()

You can also check the covariance matrix produced by the DRP against
what would be produced by the DAP for the same wavelength channel:

.. code-block:: python

    # Imports
    from mangadap.config import defaults
    from mangadap.datacube import MaNGADataCube

    # Read the MaNGA datacube
    path = defaults.drp_directory_path()
    cube = MaNGADataCube.from_plateifu(7815, 3702, directory_path=path, covar_ext='GCORREL')

    # Load the RSS spectra. This will only work if the RSS file is
    # either in the same directory as the datacube (i.e.
    # cube.directory_path) or at the nominal path (these are the same
    # in this example).
    cube.load_rss()

    # Check the data read by the MaNGADataCube instantiation vs a brute
    # force read of the correlation matrix
    hdu = fits.open(cube.file_path())
    channel = hdu['GCORREL'].header['BBINDEX']
    gcorrel = numpy.zeros(eval(hdu['GCORREL'].header['COVSHAPE']), dtype=float)
    i = numpy.ravel_multi_index((hdu['GCORREL'].data['INDXI_C1'],hdu['GCORREL'].data['INDXI_C2']),
                                cube.spatial_shape)
    j = numpy.ravel_multi_index((hdu['GCORREL'].data['INDXJ_C1'],hdu['GCORREL'].data['INDXJ_C2']),
                                cube.spatial_shape)
    gcorrel[i,j] = hdu['GCORREL'].data['RHOIJ']
    gcorrel[j,i] = hdu['GCORREL'].data['RHOIJ']

    # This should pass if the two methods of reading the correlation
    # matrix produce the same result.
    assert numpy.allclose(cube.covar.toarray(), gcorrel), 'Bad covariance read'

    # Check the from-scratch calculation of the rectified wavelength
    # channel and the associated correlation matrix against the data
    # provided by the DRP
    flux, C = cube.rss.rectify_wavelength_plane(channel, return_covar=True)

    # Should pass if the rectified flux is the same as produced by the
    # DRP and DAP.
    assert numpy.allclose(cube.flux[...,channel], flux), 'Bad flux rectification'

    # Should pass if the rectified inverse variance is the same as
    # produced by the DRP and DAP.
    ivar = numpy.ma.power(C.variance().reshape(cube.spatial_shape), -1).filled(0.0)
    assert numpy.allclose(cube.ivar[...,channel], ivar), 'Bad inverse variance rectification'

    # Should pass if the correlation matrix is the same as produced by
    # the DRP and DAP.
    C.to_correlation()
    assert numpy.allclose(C.toarray(), gcorrel), 'Bad covariance calculation'

Calculating additional covariance matrices
------------------------------------------

You can use the DAP to calculate the covariance matrix for any
channel in the DRP-provided datacubes.

To calculate the **formal covariance matrix**, use
:func:`~mangadap.datacube.datacube.DataCube.covariance_matrix` or
:func:`~mangadap.datacube.datacube.DataCube.covariance_cube`:

.. code-block:: python

    # Read the datacube (assumes the default paths)
    from mangadap.datacube import MaNGADataCube
    cube = MaNGADataCube.from_plateifu(7815, 3702)

    # Load the source row-stacked spectra
    cube.load_rss()

    # Construct a covariance matrix for the (single) wavelength channel
    # with index 1000
    C = cube.covariance_matrix(1000)

    # Try multiple channels
    C = cube.covariance_cube(channels=[1000,2000])
    assert numpy.array_equal(C.input_indx, [1000,2000]), 'Bad matrix indices'
    assert C.shape == (1764, 1764, 2), 'Bad covariance shape'

Note that in the above, you can calculate the covariance matrices for
*all* wavelength channels if you do not specify the channels in the
call to :func:`~mangadap.datacube.datacube.DataCube.covariance_cube`.

In addition to the formal covariance matrices, we found that the
correlation matrix is well approximated by adopting a Gaussian in
pixel separation for the form of the correlation coefficients
(:math:`{\mathbf \rho}`; see
:func:`~mangadap.datacube.manga.MaNGADataCube.approximate_covariance_matrix`).
To calculate the **approximate correlation matrix**:

.. code-block:: python

    # Read the datacube (assumes the default paths)
    from mangadap.datacube import MaNGADataCube
    cube = MaNGADataCube.from_plateifu(7815, 3702)

    # Try to generate an approximate correlation matrix
    approxC = cube.approximate_correlation_matrix()

    # To get the covariance matrix for a specific wavelength channel,
    # renormalize the approximate correlation matrix by the variance in
    # the desired channel.
    var = numpy.ma.power(cube.ivar[:,:,1000].ravel(), -1).filled(0.0)
    approxC = approxC.apply_new_variance(var)
    approxC.revert_correlation()
    approxC.show()

