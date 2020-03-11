
The DAP :class:`~mangadap.util.covariance.Covariance` object is a
general utility for constructing, visualizing, and storing two- and
three-dimensional covariance matrices.

.. _covariance-construction:

Construction
------------

Beyond the nominal instantiation method, there are numerous methods
used to construct a :class:`~mangadap.util.covariance.Covariance`
object.

 #. The simplest approaches are if you have a variance vector or a
    pre-calculated covariance matrix as a dense array and you want to
    construct a :class:`~mangadap.util.covariance.Covariance` object
    for further calculations:

    .. code-block:: python

        import numpy
        from mangadap.util.covariance import Covariance

        # Construct the Covariance object just from a variance vector
        var = numpy.ones(3, dtype=float)
        covar = Covariance.from_variance(var)
        # Should be the same as the identity matrix.
        assert numpy.array_equal(covar.toarray(), numpy.identity(3))

        # Construct the Covariance object from a pre-calculated array
        c = numpy.diag(numpy.full(10-2, 0.2, dtype=float), k=-2) \
                + numpy.diag(numpy.full(10-1, 0.5, dtype=float), k=-1) \
                + numpy.diag(numpy.full(10, 1.0, dtype=float), k=0) \
                + numpy.diag(numpy.full(10-1, 0.5, dtype=float), k=1) \
                + numpy.diag(numpy.full(10-2, 0.2, dtype=float), k=2)
        covar = Covariance.from_array(c)
        # Should be the same as the identity matrix.
        assert numpy.array_equal(covar.toarray(), c)

    Note the use of the
    :func:`~mangadap.util.covariance.Covariance.toarray` to access the
    array; see :ref:`covariance-access`.

 #. You can construct a covariance matrix based on samples from a
    distribution:

    .. code-block:: python

        import numpy
        from mangadap.util.covariance import Covariance

        # Build a bogus covariance matrix
        m = numpy.zeros(10, dtype=float)
        c = numpy.diag(numpy.full(10-2, 0.2, dtype=float), k=-2) \
                + numpy.diag(numpy.full(10-1, 0.5, dtype=float), k=-1) \
                + numpy.diag(numpy.full(10, 1.0, dtype=float), k=0) \
                + numpy.diag(numpy.full(10-1, 0.5, dtype=float), k=1) \
                + numpy.diag(numpy.full(10-2, 0.2, dtype=float), k=2)

        # Draw samples
        s = numpy.random.multivariate_normal(m, c, size=100000)

        # Instantiate
        covar = Covariance.from_samples(s.T, cov_tol=0.1)

        # Check the values are very nearly the same as the input
        assert numpy.all(numpy.absolute(c - covar.toarray()) < 0.02)

    Here, we've drawn samples from a known multivariate normal
    distribution with a given covariance between its 10 axes and
    checked the reconstruction of the covariance matrix from those
    samples using
    :func:`~mangadap.util.covariance.Covariance.from_samples`. This
    construction method is a simple wrapper for `numpy.cov`_ and
    :func:`~mangadap.util.covariance.Covariance.from_array`.

 #. You can construct the covariance matrix that results from a matrix
    multiplication:

    .. code-block:: python

        import numpy
        from mangadap.util.covariance import Covariance

        # Build a bogus covariance matrix
        c = numpy.diag(numpy.full(10-2, 0.2, dtype=float), k=-2) \
                + numpy.diag(numpy.full(10-1, 0.5, dtype=float), k=-1) \
                + numpy.diag(numpy.full(10, 1.0, dtype=float), k=0) \
                + numpy.diag(numpy.full(10-1, 0.5, dtype=float), k=1) \
                + numpy.diag(numpy.full(10-2, 0.2, dtype=float), k=2)

        # Vector and matrix to multiply
        x = numpy.ones(10, dtype=float)
        t = numpy.zeros((3,10), dtype=float)
        t[0,0] = 1.0
        t[1,2] = 1.0
        t[2,4] = 1.0
        # Apply the matrix multiplication (not used; just for illustration)
        y = numpy.dot(t, x)
        # Expected covariance from the result
        _c = numpy.diag(numpy.full(3-1, 0.2, dtype=float), k=-1) \
                + numpy.diag(numpy.full(3, 1.0, dtype=float), k=0) \
                + numpy.diag(numpy.full(3-1, 0.2, dtype=float), k=1)
        covar = Covariance.from_matrix_multiplication(t, c)
        assert numpy.array_equal(covar.toarray(), _c)

 #. Finally, you can construct the covariance matrix from a previous
    instance that was saved to a fits file using the
    :ref:`covariance-fitsio`.

.. _covariance-access:

Accessing the covariance data
-----------------------------

The :class:`~mangadap.util.covariance.Covariance` object is primarily
a storage and IO utility. Internally, the object only keeps the upper
triangle of the matrix, which means that use of the :attr:`cov`
attribute is not recommended unless you know what you're doing.

Also note that the object can be either "2D" or "3D". When the object
is 3D (``covar.dim == 3``), this is just a convenience for storing
multiple covariance matrices that have the same shape.

There are two ways to access the full covariance matrix, the
:func:`~mangadap.util.covariance.Covariance.full` and
:func:`~mangadap.util.covariance.Covariance.toarray` methods
depending on whether you want a sparse or dense matrix, respectively.
Also note that the use of the
:func:`~mangadap.util.covariance.Covariance.full` method requires you
to specify a single channel for 3D objects. The output of these two
methods can be used as you would use any `scipy.sparse.csr_matrix`_
or `numpy.ndarray`_ object.

To show the covariance matrix, you can use its
:func:`~mangadap.util.covariance.Covariance.show` method to quickly
produce a plot, which is a simple wrapper for the
:func:`~mangadap.util.covariance.Covariance.toarray` method and
`pyplot.imshow`_.

.. _covariance-correl:

Toggling between covariance and correlation matrices
----------------------------------------------------

The :class:`~mangadap.util.covariance.Covariance` object allows you to
toggle between the full covariance matrix, :math:`{\mathbf C}` and a
correlation matrix, :math:`{\mathbf \rho}`, where:

.. math::

    \rho_{ij} = \frac{C_{ij}}{(V_i V_j)^{1/2}}

where :math:`{\mathbf V}` is the variance vector (the diagonal
elements of :math:`{\mathbf C}`). To convert a
:class:`~mangadap.util.covariance.Covariance` object to a correlation
matrix (or ensure that it already is one), use
:func:`~mangadap.util.covariance.Covariance.to_correlation`. To revert
back to a covariance matrix, use
:func:`~mangadap.util.covariance.Covariance.revert_correlation`.

.. _covariance-fitsio:

Fits file I/O methods
---------------------

:class:`~mangadap.util.covariance.Covariance` objects can be saved as
a binary table in a fits file using the
:func:`~mangadap.util.covariance.Covariance.write` method. To reload
the covariance matrix, use the
:func:`~mangadap.util.covariance.Covariance.from_fits` instantiation
method:

.. code-block:: python

    import numpy
    from mangadap.util.covariance import Covariance

    ofile = 'test_covar_io.fits'

    # Build a bogus covariance matrix
    m = numpy.zeros(10, dtype=float)
    c = numpy.diag(numpy.full(10-2, 0.2, dtype=float), k=-2) \
            + numpy.diag(numpy.full(10-1, 0.5, dtype=float), k=-1) \
            + numpy.diag(numpy.full(10, 1.0, dtype=float), k=0) \
            + numpy.diag(numpy.full(10-1, 0.5, dtype=float), k=1) \
            + numpy.diag(numpy.full(10-2, 0.2, dtype=float), k=2)

    # Draw samples
    s = numpy.random.multivariate_normal(m, c, size=100000)

    # Instantiate
    covar = Covariance.from_samples(s.T, cov_tol=0.1)
    # Write
    covar.write(ofile)
    # Read
    _covar = Covariance.from_fits(ofile)
    # Should be the same
    assert numpy.allclose(covar.toarray(), _covar.toarray())

The details of how the covariance data are stored are described by
the :func:`~mangadap.util.covariance.Covariance.write` method.

