
The DAP :class:`mangadap.util.covariance.Covariance` object is a
general utility for constructing, visualizing, and storing two- and
three-dimensional covariance matrices.

.. _covariance-construction:

Construction
------------

Beyond the nominal instantiation method, there are numerous methods
used to construct a :class:`mangadap.util.covariance.Covariance`
object.

 #. The simplest approaches are if you have a variance vector or a
    pre-calculated covariance matrix as a dense array and you want to
    construct a :class:`mangadap.util.covariance.Covariance` object
    for further calculations::

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
    :func:`mangadap.util.covariance.Covariance.toarray` to access the
    array; see :ref:`covariance-access`.

 #. You can also construct a covariance matrix provided a set of
    samples from a distribution::

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
    :func:`mangadap.util.covariance.Covariance.from_samples`. This
    construction method is a simple wrapper for `numpy.cov`_ and
    :func:`mangadap.util.covariance.Covariance.from_array`.

 #. You can also construct the covariance matrix that results from
    a matrix multiplication::

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

Here's how you access the data.

.. _covariance-fitsio:

Fits file I/O methods
---------------------
