# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
r"""
Defines a class used to store and interface with covariance matrices.

.. todo::
    - Allow for calculation of the inverse of the covariance matrix.
    - Instead of 3D covariance cubes being an array of sparse objects,
      make the whole thing a sparse array.

Usage examples
--------------

You can calculate the covariance matrix for a given wavelength
channel in a :class:`mangadap.drpfits.DRPFits` object::

    # Access the DRP RSS file
    from mangadap.drpfits import DRPFits
    drpf = DRPFits(7495, 12703, 'RSS', read=True)

    # Calculate a single covariance matrix
    C = drpf.covariance_matrix(2281)

    # Show the result in an image
    C.show()

    # Access specific elements
    print(C[0,0])

    # Convert to a 'dense' array
    dense_C = C.toarray()

    # Get the triplets of the non-zero elements
    i, j, v = C.find()

    # Write it to disk (clobber existing file)
    C.write('test_covariance.fits', clobber=True)

The covariance matrix is stored in "coordinate" format in a fits
binary table. Since the covariance matrix is symmetric by definition,
only those non-zero elements in the upper triangle (:math:`C_{ij}`
where :math:`i\leq j`) are saved (in memory or on disk). You can read
an existing covariance matrix fits file::

    from mangadap.util.covariance import Covariance
    C = Covariance(ifile='test_covariance.fits')

You can calculate a set of covariance matrices or the full covariance
cube::

    # Access the DRP RSS file
    from mangadap.drpfits import DRPFits
    drpf = DRPFits(7495, 12703, 'RSS', read=True)

    # Calculate the full cube
    CC = drpf.covariance_cube()             # BEWARE: This may require a LOT of memory

    # Calculate fewer but many covariance matrices
    channels = [ 0, 1000, 2000, 3000, 4000 ]
    C = drpf.covariance_cube(channels=channels)

    # Access individual elements
    print(C[0,0,0])
    print(C[0,0,1000])

    # Write the full cube, or set of channels
    CC.write('full_covariance_cube.fits')   # BEWARE: This will be a BIG file
    C.write('covariance_channel_set.fits')

Although you can access the data in the covariance matrix as
explained above, this is generally inefficient because the
:func:`Covariance.__getitem__` function currently cannot handle
slices. If you need to perform operations with the covariance matrix,
you're better off working with the :attr:`Covariance.cov` attribute
directly. To do this, you won't be able to use the aliasing of the
channel indices for 3D covariance matrices.

The :class:`Covariance` class also allows you to toggle between
accessing the matrix as a true covariance matrix, or by splitting the
matrix into its variance and correlation components. For example,::

    # Get the covariance matrix for a single wavelength channel
    from mangadap.drpfits import DRPFits
    drpf = DRPFits(7495, 12703, 'RSS', read=True)
    C = drpf.covariance_matrix(2281)

    # Show the covariance matrix before and after changing it to a
    # correlation matrix
    C.show()
    C.to_correlation()
    C.show()
    print(C.is_correlation)
    C.revert_correlation()

Covariance matrices that have been converted to correlation matrices
can be written and read in without issue. See
:func:`Covariance.write` and :func:`Covariance.read`. For example::

    # Get the covariance matrix for a single wavelength channel
    import numpy
    from mangadap.util.covariance import Covariance
    from mangadap.drpfits import DRPFits

    drpf = DRPFits(7495, 12703, 'RSS', read=True)
    channels = [ 0, 1000, 2000, 3000, 4000 ]
    Cov = drpf.covariance_cube(channels=channels)

    Cov.to_correlation()
    Cov.show(channel=2000)
    Cov.write('correlation_matrix.fits', clobber=True)
    Cov.revert_correlation()

    Corr = Covariance(ifile='correlation_matrix.fits')
    Corr.revert_correlation()

    assert not (numpy.abs(numpy.sum(Cov.toarray(channel=2000) 
                        - Corr.toarray(channel=2000))) > 0.0)

Revision history
----------------

    | **23 Feb 2015**: Original Implementation by K. Westfall (KBW)
    | **04 Aug 2015**: (KBW) Sphinx documentation and minor edits.
    | **29 Mar 2016**: (KBW) Allow the object to flip between a true
        covariance matrix or a correlation matrix with saved variances.
        Allow for a function that only returns the
        `astropy.io.fits.BinTableHDU`_ object with the coordinate-format
        covariance data.
    | **06 Apr 2016**: (KBW) Allow :func:`Covariance.read` to read from
        a file or an `astropy.io.fits.hdu.hdulist.HDUList`_ object, and
        allow the specification of the extensions to read the header,
        covariance, and plane data.
    | **23 Feb 2017**: (KBW) Major revisions: Use DAPFitsUtil to
        read/write the HDUList.  Output is always a correlation matrix.
        Rearrange ordering of 3D covariance matrices to reflect the
        ordering in the DRPFits cube and the DAP Maps channels.  Change
        *plane* keyword to *channel*.  Change the format of the
        covariance matrices when written/read from a file to match what
        is done by the DRP.  Change read method to :func:`from_fits`
        class method.  Include :func:`from_samples` class method.

----

.. include license and copyright
.. include:: ../copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../links.rst
"""
import os
import numpy
import warnings

from scipy import sparse
from astropy.io import fits
from matplotlib import pyplot

from .fitsutil import DAPFitsUtil

class Covariance:
    r"""
    A general utility for storing, manipulating, and file I/O sparse
    covariance matrices.

    .. todo::
        - Can create empty object.  Does this make sense?

    Args:
        inp (`scipy.sparse.csr_matrix`_, `numpy.ndarray`_):
            Covariance matrix to store. Input **must** be covariance
            data, not correlation data. Data type can be either a
            single `scipy.sparse.csr_matrix`_ object or a 1-D array
            of them. Assumes all sparse matrices in the input ndarray
            have the same size. If None, the covariance object is
            instantiated empty.
        input_indx (`numpy.ndarray`_, optional):
            If *inp* is an array of `scipy.sparse.csr_matrix`_
            objects, this is an integer array specifying a pseudo-set
            of indices to use instead of the direct index in the
            array. I.e., if *inp* is an array with 5 elements, one
            can provide a 5-element array for :attr:`input_indx` that
            are used instead as the designated index for the
            covariance matrices. See: :func:`__getitem__`,
            :func:`_grab_true_index`.
        impose_triu (:obj:`bool`, optional):
            Flag to force the `scipy.sparse.csr_matrix`_ object to
            only be the upper triangle of the covariance matrix. The
            covariance matrix is symmetric such that C_ij = C_ji, so
            it's not necessary to keep both values. This flag will
            force a call to `scipy.sparse.triu`_ when setting the
            covariance matrix. Otherwise, the input matrix is
            **assumed** to only have the upper triangle of numbers.
        correlation (:obj:`bool`, optional):
            Convert the input to a correlation matix. The input
            **must** be the covariance matrix.

    Raises:
        TypeError:
            Raised if the input array is not one-dimensional or the
            input covariance matrix are not
            `scipy.sparse.csr_matrix`_ objects.
        ValueError:
            Raised if the :attr:`input_indx` either does not have the
            correct size or is not required due to the covariance
            object being a single matrix.

    Attributes:
        cov (`scipy.sparse.csr_matrix`_, `numpy.ndarray`_):
            The covariance matrix stored in sparse format.
        shape (:obj:`tuple`):
            Shape of the full array.
        dim (:obj:`int`):
            The number of dimensions in the covariance matrix. This
            can either be 2 (a single covariance matrix) or 3
            (multiple covariance matrices).
        nnz (:obj:`int`):
            The number of non-zero covariance matrix elements.
        input_indx (`numpy.ndarray`_):
            If :attr:`cov` is an array of `scipy.sparse.csr_matrix`_
            objects, this is an integer array specifying a pseudo-set
            of indices to use instead of the direct index in the
            array. I.e., if :attr:`cov` is an array with 5 elements,
            one can provide a 5-element array for :attr:`input_indx`
            that are used instead as the designated index for the
            covariance matrices. See: :func:`__getitem__`,
            :func:`_grab_true_index`.
        inv (`scipy.sparse.csr_matrix`_):
            The inverse of the covariance matrix. **This is not
            currently calculated!**
        var (`numpy.ndarray`):
            Array with the variance provided by the diagonal of
            the/each covariance matrix. This is only populated if
            necessary, either by being requested (:func:`variance`)
            or if needed to convert between covariance and
            correlation matrices.
        is_correlation (:obj:`bool`):
            Flag that the covariance matrix has been saved as a
            variance vector and a correlation matrix.

    """
    def __init__(self, inp, input_indx=None, impose_triu=False, correlation=False):

        self.cov = inp
        self.shape = None
        self.dim = None
        self.nnz = None
        self.input_indx = None
        self.inv = None
        self.var = None
        self.is_correlation = False

        # Return empty object
        if self.cov is None:
            return

        # Set the dimensionality, check that each element of the array
        # has the correct type, count the number of non-zero covariance
        # values, and (if necessary) initialize the indices for each
        # covariance matrix
        if isinstance(self.cov, numpy.ndarray):
            if len(self.cov.shape) > 1:
                raise TypeError('Input ndarray can only be one-dimensional')
            self.dim = 3
            self.nnz = 0
            for cov in self.cov:
                if not sparse.isspmatrix_csr(cov):
                    raise TypeError('Input covariance matrix (or elements) must be csr_matrices.')
                self.nnz += cov.nnz

            if input_indx is not None:
                if input_indx.shape != self.cov.shape:
                    raise ValueError('Input array and input index array must be the same size.')
                if numpy.unique(input_indx).size != input_indx.size:
                    raise ValueError('Input list of indices must all be unique!')
                self.input_indx = input_indx
            else:
                self.input_indx = numpy.arange(len(self.cov))
        else:
            self.dim = 2
            if not sparse.isspmatrix_csr(self.cov):
                raise TypeError('Input covariance matrix (or elements) must be csr matrices.')
            self.nnz = self.cov.nnz
            if input_indx is not None:
                raise ValueError('Input indices only allowed when allocating multiple matrices.')

        # Set the shape of the full matrix/cube
        self._set_shape()

        # If requested, impose that the input matrix only have values in
        # its upper triangle.
        if impose_triu:
            self._impose_upper_triangle()

        # Set the inverse of the covariance matrix
        # **NOT IMPLEMENTED YET**

        # Set the variance array and the correlation matrix flag
        if correlation:
            self.to_correlation()
    
    def __getitem__(self, *args):
        """
        Return the covariance value at a provided 2D or 3D position.

        Args:
            *args (tuple):
                2 or 3 integers designating the covariance value to
                return. Number of values must match the
                dimensionality of the object.

        Returns:
            :obj:`float`: The value of the covariance matrix at the
            designated index.

        Raises:
            ValueError:
                Raised if the number of arguments does not match the
                dimensionality of the object.
        """
        indx = tuple(*args)
        if len(indx) != self.dim:
            raise ValueError('Incorrect number of dimensions!')
        if self.dim == 2:
            return self.cov[tuple(sorted(indx))]
        return self.cov[self._grab_true_index(indx[2])][tuple(sorted(indx[:2]))]

    @classmethod
    def from_samples(cls, samples, cov_tol=None, rho_tol=None):
        r"""
        Define a covariance object using descrete samples.

        The covariance is generated using `numpy.cov`_ for a set of
        discretely sampled data for an :math:`N`-dimensional
        parameter space.

        Args:
            samples (`numpy.ndarray`_):
                Array with samples drawn from an
                :math:`N`-dimensional parameter space. The shape of
                the input array must be :math:`N_{\rm par}\times
                N_{\rm samples}`.
            cov_tol (:obj:`float`, optional):
                Any covariance value less than this is assumed to be
                equivalent to (and set to) 0.
            rho_tol (:obj:`float`, optional):
                Any correlation coefficient less than this is assumed
                to be equivalent to (and set to) 0.

        Returns:
            :class:`Covariance`: An :math:`N_{\rm par}\times N_{\rm
            par}` covariance matrix built using the provided samples.

        """
        if len(samples.shape) != 2:
            raise ValueError('Input samples for covariance matrix must be a 2D array!')

        npar = samples.shape[0]

        covariance = numpy.cov(samples)
        if rho_tol is not None:
            variance = numpy.diag(covariance)
            correlation = covariance / numpy.ma.sqrt(variance[:,None]*variance[None,:])
            correlation[ numpy.ma.absolute(correlation) < rho_tol ] = 0.0
            covariance = correlation.filled(0.0) \
                            * numpy.ma.sqrt(variance[:,None]*variance[None,:]).filled(0.0)
        if cov_tol is not None:
            covariance[ covariance < cov_tol ] = 0.0

        indx = covariance > 0.0
        i, j = numpy.meshgrid(numpy.arange(npar), numpy.arange(npar), indexing='ij')

        return cls(sparse.coo_matrix((covariance[indx].ravel(),
                                        (i[indx].ravel(), j[indx].ravel())),
                                     shape=(npar,npar)).tocsr(), impose_triu=True)

    @classmethod
    def from_fits(cls, source, ivar_ext='IVAR', covar_ext='CORREL', row_major=False,
                  impose_triu=False, correlation=False, quiet=False):
        r"""
        Read an existing covariance object previously written to disk;
        see :func:`write`.  The class can read covariance data written
        by other programs *as long as they have a commensurate format*.
        
        If the extension names and column names are correct,
        :class:`Covariance` can read fits files that were not
        necessarily produced by the class itself.  This is useful for
        MaNGA DAP products that include the covariance data in among
        other output products.  Because of this, :func:`output_hdus` and
        :func:`output_fits_data` are provided to produce the data that
        can be placed in any fits file.

        The function determines if the output data were reordered by
        checking the number of binary columns.

        Args:
            source (:obj:`str`, `astropy.io.fits.hdu.hdulist.HDUList`_):
                Initialize the object using an
                `astropy.io.fits.hdu.hdulist.HDUList`_ object or path
                to a fits file.
            ivar_ext (:obj:`str`, optional):
                If reading the data from ``source``, this is the name
                of the extension with the inverse variance data.
                Default is ``'IVAR'``. If None, the variance is taken
                as unity.
            covar_ext (:obj:`str`, optional):
                If reading the data from ``source``, this is the name
                of the extension with covariance data. Default is
                ``'CORREL'``.
            row_major (:obj:`bool`, optional):
                If reading the data from an
                `astropy.io.fits.hdu.hdulist.HDUList`_, this sets if
                the data arrays have been rearranged into a
                python-native (row-major) structure. See
                :func:`mangadap.util.fitsutil.DAPFitsUtil.transpose_image_data`.
                Default is False.
            impose_triu (:obj:`bool`, optional):
                Flag to force the `scipy.sparse.csr_matrix`_ object
                to only be the upper triangle of the covariance
                matrix. The covariance matrix is symmetric such that
                :math:`C_{ij} = C_{ji}`, so it's not necessary to
                keep both values. This flag will force a call to
                `scipy.sparse.triu`_ when setting the covariance
                matrix. Otherwise, the input matrix is *assumed* to
                only have the upper triangle of numbers.
            correlation (:obj:`bool`, optional):
                Return the matrix as a correlation matrix. Default
                (False) is to use the data (always saved in
                correlation format; see :func:`write`) to construct
                the covariance matrix.
            quiet (:obj:`bool`, optional):
                Suppress terminal output.
        """
        if isinstance(source, fits.HDUList):
            hdu = source if row_major else DAPFitsUtil.transpose_image_data(source)
        else:
            try:
                hdu = DAPFitsUtil.read(source)
            except Exception as e:
                print(e)
                warnings.warn('Problem reading covariance from file.')
                return

        # Read a coordinate data
        shape = eval(hdu[covar_ext].header['COVSHAPE'])
        dim = len(shape)
        reshape = len(hdu[covar_ext].columns.names) in [ 5, 6 ]
        if reshape:
            i_c1, i_c2, j_c1, j_c2, rhoij = [ hdu[covar_ext].data[ext].copy() for ext in
                                                [ 'INDXI_C1', 'INDXI_C2', 'INDXJ_C1', 'INDXJ_C2',
                                                  'RHOIJ' ] ]
            n = Covariance.reshape_size(shape[0])
            i = Covariance.ravel_indices(n, i_c1, i_c2)
            j = Covariance.ravel_indices(n, j_c1, j_c2)
        else:
            i, j, rhoij = [ hdu[covar_ext].data[ext] for ext in ['INDXI', 'INDXJ', 'RHOIJ'] ]

        # Number of non-zero elements
        nnz = len(rhoij)

        # Set correlation data
        if dim == 2:
            var = numpy.ones(shape[1:], dtype=float) if ivar_ext is None \
                    else numpy.ma.power(hdu[ivar_ext].data.ravel(), -1).filled(0.0)
            cij = rhoij * numpy.sqrt( var[i]*var[j] )
            cov = sparse.coo_matrix( (cij, (i, j)), shape=shape).tocsr()
            input_indx = None
        else:
            k = hdu[covar_ext].data['INDXK'].copy()
            input_indx = numpy.unique(k)
            if len(input_indx) > shape[-1]:
                raise ValueError('Number of unique channel indices is greater than provided shape.')
            if len(input_indx) < shape[-1]:
                warnings.warn('Fewer unique channels than all available.')

            var = numpy.ones(shape[1:], dtype=float) if ivar_ext is None \
                    else numpy.ma.power(hdu[ivar_ext].data.reshape(-1, shape[-1]), -1).filled(0.0)
            cov = numpy.empty(shape[-1], dtype=sparse.csr.csr_matrix)
            for ii, uk in enumerate(input_indx):
                indx = k == uk
                cij = rhoij[indx] * numpy.sqrt( var[i[indx],ii]*var[j[indx],ii] )
                cov[ii] = sparse.coo_matrix((cij, (i[indx], j[indx])), shape=shape[:-1]).tocsr()

        # Report
        # TODO: Convert report to use logging
        if not quiet:
            print('Read covariance cube:')
            print('       output type: {0}'.format('Correlation' \
                                                    if correlation else 'Covariance'))
            print('        dimensions: {0}'.format(dim))
            print('             shape: {0}'.format(shape))
            if dim == 3:
                print('    pseudo-indices: ', input_indx)
            print('   non-zero values: {0}'.format(nnz))

        return cls(cov, input_indx=input_indx, correlation=correlation)

    @classmethod
    def from_matrix_multiplication(cls, T, Sigma):
        r"""
        Construct the covariance matrix that results from a matrix
        multiplication.
        
        The matrix multiplication should be of the form:

        .. math::
        
            {\mathbf T} \times {\mathbf X} = {\mathbf Y}

        where :math:`{\mathbf T}` is a transfer matrix of size
        :math:`N_y\times N_x`, :math:`{\mathbf X}` is a vector of
        size :math:`N_x`, and :math:`{\mathbf Y}` is the vector of
        length :math:`{N_y}` that results from the multiplication.

        The covariance matrix is then
        
        .. math::

             {\mathbf C} = {\mathbf T} \times {\mathbf \Sigma} \times
             {\mathbf T}^{\rm T},

        where :math:`{\mathbf \Sigma}` is the covariance matrix for
        the elements of :math:`{\mathbf X}`. If `Sigma` is provided
        as a vector of length :math:`N_x`, it is assumed that the
        elements of :math:`{\mathbf X}` are independent and the
        provided vector gives the variance in each element; i.e., the
        provided data represent the diagonal of :math:`{\mathbf
        \Sigma}`.

        Args:
            T (`scipy.sparse.csr_matrix`_, `numpy.ndarray`_):
                Transfer matrix.  See above.
            Sigma (`scipy.sparse.csr_matrix`_, `numpy.ndarray`_):
                Covariance matrix.  See above.
        """
        if len(T.shape) != 2:
            raise ValueError('Input transfer matrix must be two-dimensional.')
        nx = T.shape[1]
        if Sigma.shape != (nx,nx) and Sigma.shape != (nx,):
            raise ValueError('Shape of input variance matrix must be either '
                             '({0},{0}) or ({0},).'.format(nx))
        # If it isn't already, convert T to a csr_matrix
        _T = T if isinstance(T, sparse.csr.csr_matrix) else sparse.csr_matrix(T)
        # Set the covariance matrix in X
        _Sigma = sparse.coo_matrix( (Sigma, (numpy.arange(nx),numpy.arange(nx))),
                                   shape=(nx,nx)).tocsr() if len(Sigma.shape) == 1 \
                    else (Sigma if isinstance(Sigma, sparse.csr.csr_matrix) \
                                else sparse.csr_matrix(Sigma))
        # Construct the covariance matrix
        return cls( sparse.triu(_T.dot(_Sigma.dot(_T.transpose()))).tocsr() )


    @classmethod
    def from_variance(cls, variance, correlation=False):
        r"""
        Construct a diagonal covariance matrix using the provided variance.

        Args:
            variance (`numpy.ndarray`_):
                The variance vector.
            correlation (:obj:`bool`, optional):
                Upon instantiation, convert the :class:`Covariance`
                object to a correlation matrix.
        """
        return cls(sparse.csr.csr_matrix(numpy.diagflat(variance)), correlation=correlation)

    def _grab_true_index(self, inp):
        """
        In the case of a 3D array, return the true array-element
        index given the pseudo index.

        .. todo::
            Search operation is inefficient.  Apparently a better option
            is a feature request for numpy 2.0

        Args:
            inp (:obj:`int`):
                Requested index to convert to the real index in the
                stored array of covariance matrices.

        Returns:
            :obj:`int`: The index in the stored array of the
            requested covariance matrix.
        """
        if self.input_indx is not None:
            i = numpy.where(self.input_indx == inp)[0]
            if i.size != 1:
                raise IndexError('Invalid index: {0}!'.format(inp))
            return i[0]
        return inp

    def _set_shape(self):
        """Set the 2D or 3D shape of the covariance matrix."""
        self.shape = self.cov.shape
        if self.dim == 2:
            return
        self.shape = self.cov.ravel()[0].shape + self.shape

    def _impose_upper_triangle(self):
        """
        Force :attr:`cov` to only contain non-zero elements in its
        upper triangle.
        """
        if self.dim == 2:
            self.cov = sparse.triu(self.cov).tocsr()
            self.nnz = self.cov.nnz
            return

        self.nnz = 0
        for i in range(self.shape[0]):
            self.cov[i] = sparse.triu(self.cov[i]).tocsr()
            self.nnz += self.cov[i].nnz


    def _with_lower_triangle(self, channel=None):
        r"""
        Return a `scipy.sparse.csr_matrix`_ object with both its
        upper and lower triangle filled, ensuring that they are
        symmetric.

        Args:
            channel (:obj:`int`, optional):
                The pseudo-index of the covariance matrix to return.

        Returns:
            `scipy.sparse.csr_matrix`_: The sparse matrix with both
            the upper and lower triangles filled (with symmetric
            information).

        Raises:
            ValueError:
                Raised if the object is 3D and *channel* is not
                provided.
        """
        if self.dim == 2:
            a = self.cov
        else:
            if channel is None:
                raise ValueError('Must define channel!  Use channel=...')
            a = self.cov[self._grab_true_index(channel)]

        return (sparse.triu(a) + sparse.triu(a,1).T)

    @staticmethod
    def reshape_size(size):
        """
        Determine the shape of a 2D square array resulting from
        reshaping a vector.
        
        Args:
            size (:obj:`int`):
                Size of the vector.

        Returns:
            :obj:`int`: Length of one axis of the square array.
        """
        # Get the size of the square image (and make sure it's square)
        n = numpy.floor(numpy.sqrt(size)).astype(int)
        if n*n != size:
            raise ValueError('{0} is not the square of an integer!'.format(size))
        return n

    @staticmethod
    def reshape_indices(size, i):
        """
        Return the indices in the 2D array that results from reshaping a
        vector of length ``size`` into a square 2D array.  The indices of
        interest in the vector are give by i.

        Args:
            size (:obj:`int`):
                Size of the vector.
            i (:obj:`int`):
                Index in the vector.

        Returns:
            tuple: The row and column indices in the square array
            corresponsding to index ``i`` in the vector.
        """
        n = Covariance.reshape_size(size)
        ii = i//n
        return ii, i-ii*n

    @staticmethod
    def ravel_indices(size, i, j):
        """
        Return the indices in the 1D vector that results from
        flattening a ``size``-by-``size`` square array.

        Args:
            size (:obj:`int`):
                Length along one axis of the square 2D array.
            i (:obj:`int`):
                Row (first axis) index.
            j (:obj:`int`):
                Column (second axis) index.

        Returns:
            :obj:`int`: The index in the flattened vector with the
            relevant value.
        """
        return i*size + j

    def apply_new_variance(self, var):
        """
        Using the same correlation coefficients, return a new
        :class:`Covariance` object with the provided variance.

        Args:
            var (`numpy.ndarray`_):
                Variance vector. Must have a length that matches the
                shape of this :class:`Covariance` instance.

        Returns:
            :class:`Covariance`: A covariance matrix with the same
            shape and correlation coefficients and this object, but
            with different variance.

        Raises:
            ValueError:
                Raised if the length of the variance vector is
                incorrect.
        """
        if var.shape != self.shape[1:]:
            raise ValueError('Provided variance has incorrect shape.')

        # Convert to a correlation matrix, if needed
        is_correlation = self.is_correlation
        if not is_correlation:
            self.to_correlation()

        if self.dim == 2:
            i, j, c = sparse.find(self.cov)
            new_cov = sparse.coo_matrix( (c*numpy.sqrt(var[i]*var[j]), (i,j)),
                                           shape=self.shape).tocsr()
        else:
            new_cov = numpy.empty(self.shape[-1], dtype=sparse.csr.csr_matrix)
            for p in range(self.shape[-1]):
                i, j, c = sparse.find(self.cov[p])
                new_cov[p] = sparse.coo_matrix((c*numpy.sqrt(var[i,p]*var[j,p]), (i,j)),
                                               shape=self.shape[:-1]).tocsr()

        # Revert to covariance matrix, if needed
        if not is_correlation:
            self.revert_correlation()

        # Return a new covariance matrix
        return Covariance(new_cov, input_indx=self.input_indx, correlation=is_correlation)

    def copy(self):
        """
        Return a copy of this Covariance object.
        """
        # If the data is saved as a correlation matrix, first revert to
        # a covariance matrix
        is_correlation=self.is_correlation
        if self.is_correlation:
            self.revert_correlation()

        # Create the new Covariance instance with a copy of the data
        cp = Covariance(self.cov.copy(), input_indx=self.input_indx.copy())

        # If necessary, convert the data to a correlation matrix
        if is_correlation:
            self.to_correlation()
            cp.to_correlation()
        return cp

    def toarray(self, channel=None):
        """
        Convert the sparse covariance matrix to a dense array, filled
        with zeros where appropriate.

        Args:
            channel (:obj:`int`, optional):
                The pseudo-index of the covariance matrix to plot.
                Required if the covariance object is 3D.

        Returns:
            `numpy.ndarray`_: Dense array with the full covariance
            matrix.
        """
        if self.dim == 2 or channel is not None:
            return (self._with_lower_triangle(channel=channel)).toarray()
        
        arr = numpy.empty(self.shape, dtype=numpy.float)
        for k in range(self.shape[-1]):
            indx = k if self.input_indx is None else self.input_indx[k]
            arr[:,:,k] = (self._with_lower_triangle(channel=indx)).toarray()
        return arr

    def show(self, channel=None, zoom=None, ofile=None, log10=False):
        """
        Show a covariance/correlation matrix data.

        This converts the (selected) covariance matrix to a filled
        array and plots the array using `matplotlib.pyplot.imshow`_.
        If an output file is provided, the image is redirected to the
        designated output file; otherwise, the image is plotted to
        the screen.

        Args:
            channel (:obj:`int`, optional):
                The pseudo-index of the covariance matrix to plot.
                Required if the covariance object is 3D.
            zoom (:obj:`float`, optional):
                Factor by which to zoom in on the center of the image
                by *removing the other regions of the array*. E.g.
                *zoom=2* will show only the central quarter of the
                covariance matrix.
            ofile (:obj:`str`, optional):
                If provided, the array is output to this file instead
                of being plotted to the screen.
        """
        # Convert the covariance matrix to an array
        a = self.toarray(channel)

        # Remove some fraction of the array to effectively zoom in on
        # the center of the covariance matrix
        if zoom is not None:
            xs = int(self.shape[self.dim-2]/2 - self.shape[self.dim-2]/2/zoom)
            xe = xs + int(self.shape[self.dim-2]/zoom) + 1
            ys = int(self.shape[self.dim-1]/2 - self.shape[self.dim-1]/2/zoom)
            ye = ys + int(self.shape[self.dim-1]/zoom) + 1
            a = a[xs:xe,ys:ye]

        # Print the plot to the screen if no output file is provided.
        if ofile is None:
            im = pyplot.imshow(numpy.ma.log10(a) if log10 else a,
                               interpolation='nearest', origin='lower')
            pyplot.colorbar()
            pyplot.show()
            return

        # Print the plot to the designated output file
        fig = pyplot.figure(1)
        im = pyplot.imshow(numpy.ma.log10(a) if log10 else a,
                           interpolation='nearest', origin='lower')
        pyplot.colorbar()
        fig.canvas.print_figure(ofile)
        fig.clear()
        
    def find(self, channel=None):
        """
        Find the non-zero values in the **full** covariance matrix (not
        just the upper triangle).

        Args:
            channel (:obj:`int`, optional):
                The pseudo-index of the covariance matrix to plot.
                Required if the covariance object is 3D.

        Returns:
            tuple: A tuple of arrays ``i``, ``j``, and ``c``. The
            arrays ``i`` and ``j`` contain the index coordinates of
            the non-zero values, and ``c`` contains the values
            themselves.
        """
        return sparse.find(self._with_lower_triangle(channel=channel))


#    def issingular(self):
#        if self.dim > 2:
#            nsingular = 0
#            for cov in self.cov.ravel():
#                print(numpy.linalg.cond(sparse.triu(cov).toarray()+sparse.triu(cov,1).T.toarray()))
#                print(1/sys.float_info.epsilon)
#        else:
#            print(numpy.linalg.cond(sparse.triu(self.cov).toarray()+sparse.triu(self.cov,1).T.toarray()))
#            print(1/sys.float_info.epsilon)

#   def inverse(self, redo=False):

#       if self.inv is not None and not redo:
#           return self.inv

#       if self.dim > 2:
#           self.inv = numpy.empty(self.shape[:self.dim-2], dtype=sparse.csr.csr_matrix)
#           print(self.inv.shape)
#           print(type(self.inv[0]))
#           for cov, inv in zip(self.cov.ravel(), self.inv.ravel()):
#               try:
#                   inv = sparse.triu(sparse.csr_matrix( \
#                                                   scipy.linalg.inv(sparse.triu(cov).toarray() \
#                                                   + sparse.triu(cov,1).T.toarray()))).tocsr()
#               except numpy.linalg.linalg.LinAlgError as e:
#                   print('Caught error: {0}'.format(e))
#           
#       else:
#           print(self.cov.nnz)
#           try:
#               self.inv = sparse.triu(sparse.csr_matrix( \
#                                               scipy.linalg.inv(sparse.triu(self.cov).toarray() \
#                                               + sparse.triu(self.cov,1).T.toarray()))).tocsr()
#           except numpy.linalg.linalg.LinAlgError as e:
#               print('Caught error: {0}'.format(e))
#           
#       return self.inv

    def output_fits_data(self, reshape=False):
        r"""
        Construct the data arrays to write to a fits file, providing the
        covariance matrix/matrices in coordinate format.
        
        Regardless of whether or not the current internal data is a
        covariance matrix or a correlation matrix, the data is always
        returned as a correlation matrix.

        If the covariance matrix is 2-dimensional, four columns are
        output:

            - :math:`i,j`: The row and column indices, respectively,
              of the covariance matrix.

            - :math:`rho_{ij}`: The correlation coefficient between
              pixels :math:`i` and :math:`j`.

            - :math:`V_i`: The variance in each pixel; i.e., the
              value of :math:`C_{ii} \forall i`. If reshape is True,
              this will be output as a two-dimenional array.

        If the covariance matrix is 3-dimensional, one additional column
        is output:

            - :math:`k`: The indices of the channels associated with
              each covariance matrix. This is either just the index
              of the covariance matrix or the provided pseudo-indices
              of each channel.

        If using the reshape option, the :math:`i,j` indices are
        converted to two columns each providing the indices in the
        associated reshaped array with coordinates :math:`c_1,c_2`.
            
        Args:
            reshape (:obj:`bool`, optional):
                Reshape the output in the case when :math:`i,j` are
                actually pixels in a two-dimensional image. The shape
                of the image is expected to be square, such that the
                shape of the covariance matrix is :math:`N_x\times
                N_y`, where :math:`N_x = N_y`. Each of the ``I`` and
                ``J`` output columns are then split into two columns
                according to associate coordinates, such that there
                are 8 output columns.

        Returns:
            tuple: Either 6 or 8 `numpy.ndarray`_ objects depending
            on if the data has been reshaped into an image.
        """
        # Only ever print correlation matrices
        is_correlation = self.is_correlation
        if not is_correlation:
            self.to_correlation()

        # Get the square size if reshaping
        reshape_n = Covariance.reshape_size(self.shape[0]) if reshape else None

        # Only one covariance matrix
        if self.dim == 2:
            # Get the data
            i, j, rhoij = sparse.find(self.cov)

            # If object was originally a covariance matrix, revert it back
            if not is_correlation:
                self.revert_correlation()

            # Return the data
            if not reshape:
                # Returns four arrays
                return i, j, rhoij, self.var.copy()

            # Returns six arrays
            i_c1, i_c2 = Covariance.reshape_indices(self.shape[0], i)
            j_c1, j_c2 = Covariance.reshape_indices(self.shape[0], j)
            return i_c1, i_c2, j_c1, j_c2, rhoij, self.var.reshape(reshape_n, reshape_n).copy()

        # More than one covariance matrix
        i = numpy.empty(self.nnz, dtype=int)
        j = numpy.empty(self.nnz, dtype=int)
        k = numpy.zeros(self.nnz, dtype=int)
        rhoij = numpy.empty(self.nnz, dtype=float)
        ii = 0
        for kk in range(self.shape[-1]):
            nnz = self.cov[kk].nnz
            k[ii:ii+nnz] = kk if self.input_indx is None else self.input_indx[kk]
            i[ii:ii+nnz], j[ii:ii+nnz], rhoij[ii:ii+nnz] = sparse.find(self.cov[kk])
            ii += nnz
        
        # If object was originally a covariance matrix, revert it back
        if not is_correlation:
            self.revert_correlation()

        # Return the data
        if not reshape:
            # Returns five arrays
            return i, j, k, rhoij, self.var.copy()

        # Returns seven arrays
        i_c1, i_c2 = Covariance.reshape_indices(self.shape[0], i)
        j_c1, j_c2 = Covariance.reshape_indices(self.shape[0], j)
        return i_c1, i_c2, j_c1, j_c2, k, rhoij, self.var.reshape(reshape_n, reshape_n, -1).copy()

    def output_hdus(self, reshape=False, hdr=None):
        r"""
        Construct the output HDUs and header that contain the covariance
        data.

        Args:
            reshape (:obj:`bool`, optional):
                Reshape the output in the case when :math:`i,j` are
                actually pixels in a two-dimensional image. The shape
                of the image is expected to be square, such that the
                shape of the covariance matrix is :math:`N_x\times
                N_y`, where :math:`N_x = N_y`. Each of the ``I`` and
                ``J`` output columns are then split into two columns
                according to associate coordinates, such that there
                are 8 output columns.
            hdr (`astropy.io.fits.Header`_, optional):
                `astropy.io.fits.Header`_ instance to which to add
                covariance keywords. If None, a new
                `astropy.io.fits.Header`_ instance is returned.

        Returns:
            tuple: Returns three objects:

                - The header for the primary HDU.
                - An `astropy.io.fits.ImageHDU`_ object with the
                  variance vector.
                - An `astropy.io.fits.BinTableHDU`_ object with the
                  correlation coefficients.
        """
        # Use input header or create a minimal one
        _hdr = fits.Header() if hdr is None else hdr
        # Ensure the input header has the correct type
        if not isinstance(_hdr, fits.Header):
            raise TypeError('Input header must have type astropy.io.fits.Header.')
        _hdr['COVSHAPE'] = (str(self.shape), 'Shape of the correlation matrix')

        tbl_hdr = fits.Header()
        tbl_hdr['COVSHAPE'] = (str(self.shape), 'Shape of the correlation matrix')

        # Construct the correlation binary table HDU
        if self.dim == 2:
            if reshape:
                i_c1, i_c2, j_c1, j_c2, rhoij, var \
                            = self.output_fits_data(reshape=reshape)
                covar_hdu = fits.BinTableHDU.from_columns([
                                fits.Column(name='INDXI_C1', format='1J', array=i_c1),
                                fits.Column(name='INDXI_C2', format='1J', array=i_c2),
                                fits.Column(name='INDXJ_C1', format='1J', array=j_c1),
                                fits.Column(name='INDXJ_C2', format='1J', array=j_c2),
                                fits.Column(name='RHOIJ', format='1D', array=rhoij)
                                                          ], name='CORREL', header=tbl_hdr)
            else:
                i, j, rhoij, var = self.output_fits_data(reshape=reshape)
                covar_hdu = fits.BinTableHDU.from_columns([
                                fits.Column(name='INDXI', format='1J', array=i),
                                fits.Column(name='INDXJ', format='1J', array=j),
                                fits.Column(name='RHOIJ', format='1D', array=rhoij)
                                                          ], name='CORREL', header=tbl_hdr)
        else:
            if reshape:
                i_c1, i_c2, j_c1, j_c2, k, rhoij, var \
                            = self.output_fits_data(reshape=reshape)
                covar_hdu = fits.BinTableHDU.from_columns([
                                fits.Column(name='INDXI_C1', format='1J', array=i_c1),
                                fits.Column(name='INDXI_C2', format='1J', array=i_c2),
                                fits.Column(name='INDXJ_C1', format='1J', array=j_c1),
                                fits.Column(name='INDXJ_C2', format='1J', array=j_c2),
                                fits.Column(name='INDXK', format='1J', array=k),
                                fits.Column(name='RHOIJ', format='1D', array=rhoij)
                                                          ], name='CORREL', header=tbl_hdr)
            else:
                i, j, k, rhoij, var = self.output_fits_data(reshape=reshape)
                covar_hdu = fits.BinTableHDU.from_columns([
                                fits.Column(name='INDXI', format='1J', array=i),
                                fits.Column(name='INDXJ', format='1J', array=j),
                                fits.Column(name='INDXK', format='1J', array=k),
                                fits.Column(name='RHOIJ', format='1D', array=rhoij)
                                                          ], name='CORREL', header=tbl_hdr)

        ivar_hdu = fits.ImageHDU(data=numpy.ma.power(var,-1.).filled(0.0), name='IVAR')

        return _hdr, ivar_hdu, covar_hdu

    def write(self, ofile, reshape=False, hdr=None, clobber=False):
        r"""
        Write the covariance object to a fits file.

        Objects written using this function can be reinstantiated
        using :func:`from_fits`.

        The covariance matrix (matrices) are stored in "coordinate"
        format using fits binary tables; see
        `scipy.sparse.coo_matrix`_. The matrix is *always* stored as
        a correlation matix, even if the object is currently in the
        state holding the covariance data.

        Independent of the dimensionality of the covariance matrix, the
        written file has a ``PRIMARY`` extension with the keyword
        ``COVSHAPE`` that specifies the original dimensions of the
        covariance matrix; see :attr:`shape`.
            
        The correlation data are written to the ``CORREL`` extension.
        The number of columns in this extension depends on the provided
        keywords; see :func:`output_fits_data`.  The column names are:
            
            - ``INDXI``, ``INDXJ``, ``INDXK``: indices in the covariance
              matrix.  The last index is provided only if the object is
              3D.  ``INDXI`` and ``INDXJ`` are separated into two
              columns if the output is reshaped; these columns are
              ``INDXI_C1``, ``INDXI_C2``, ``INDXJ_C1``, ``INDXJ_C2``.

            - ``RHOIJ``: The non-zero correlation coefficients located
              the specified coordinates

        The inverse of the variance along the diagonal of the covariance
        matrix is output in an ImageHDU in extension ``IVAR``.

        For 3D matrices, if pseudo-indices have been provided, these are
        used in the ``INDXK`` column; however, the ``COVSHAPE`` header
        keyword only gives the shape of the unique indices of the
        covariance channels.
        
        Args:
            ofile (:obj:`str`):
                File name for the output.
            reshape (:obj:`bool`, optional):
                Reshape the output in the case when :math:`i,j` are
                actually pixels in a two-dimensional image. The shape
                of the image is expected to be square, such that the
                shape of the covariance matrix is :math:`N_x\times
                N_y`, where :math:`N_x = N_y`. Each of the ``I`` and
                ``J`` output columns are then split into two columns
                according to associate coordinates, such that there
                are 8 output columns.
            hdr (`astropy.io.fits.Header`_, optional):
                A header object to include in the PRIMARY extension.
                The SHAPE keyword will be added/overwritten.
            clobber (:obj:`bool`, optional):
                Overwrite any existing file.

        Raises:
            FileExistsError:
                Raised if the output file already exists and clobber is False.
            TypeError:
                Raised if the input ``hdr`` does not have the correct
                type.
        """
        if os.path.isfile(ofile) and not clobber:
            raise FileExistsError('{0} exists!  Use \'clobber=True\' to overwrite.'.format(ofile))

        # Construct HDUList and write the fits file
        _hdr, ivar_hdu, covar_hdu = self.output_hdus(reshape=reshape, hdr=hdr)
        DAPFitsUtil.write(fits.HDUList([ fits.PrimaryHDU(header=_hdr), ivar_hdu, covar_hdu ]),
                          ofile, clobber=clobber, checksum=True)

    def variance(self, copy=True):
        """
        Return the variance vector(s) of the covariance matrix.

        Args:
            copy (:obj:`bool`, optional):
                Return a copy instead of a reference.
        """
        if self.var is not None:
            return self.var.copy() if copy else self.var

        if self.dim == 2:
            self.var = numpy.diag(self.cov.toarray()).copy()
            return self.var

        self.var = numpy.empty(self.shape[1:], dtype=numpy.float)
        for p in range(self.shape[-1]):
            self.var[:,p] = numpy.diag(self.cov[p].toarray()).copy()

        return self.var.copy() if copy else self.var

    def to_correlation(self):
        r"""
        Convert the covariance matrix into a correlation matrix by
        dividing each element by the variances.
        
        If the matrix is a correlation matrix already (see
        :attr:`is_correlation`), no operations are performed.
        Otherwise, the variance vectors are computed, if necessary, and
        used to normalize the covariance values.

        A :class:`Covariance` object can be reverted from a correlation
        matrix using :func:`revert_correlation`.
        """
        # Object is already a correlation matrix
        if self.is_correlation:
            return

        # Ensure that the variance has been calculated
        self.variance()

        self.is_correlation = True
        if self.dim == 2:
            i, j, c = sparse.find(self.cov)
            self.cov = sparse.coo_matrix( (c/numpy.sqrt(self.var[i]*self.var[j]), (i,j)),
                                          shape=self.shape).tocsr()
            return
       
        for p in range(self.shape[-1]):
            i, j, c = sparse.find(self.cov[p])
            self.cov[p] = sparse.coo_matrix( (c/numpy.sqrt(self.var[i,p]*self.var[j,p]), (i,j)),
                                             shape=self.shape[:-1]).tocsr()

    def revert_correlation(self):
        r"""
        Revert the object from a correlation matrix back to a full
        covariance matrix.

        This function does nothing if the correlation flag has not been
        flipped.  The variances must have already been calculated!
        """
        if not self.is_correlation:
            return

        if self.dim == 2:
            i, j, c = sparse.find(self.cov)
            self.cov = sparse.coo_matrix( (c*numpy.sqrt(self.var[i]*self.var[j]), (i,j)),
                                           shape=self.shape).tocsr()
            self.is_correlation = False
            return
       
        for p in range(self.shape[-1]):
            i, j, c = sparse.find(self.cov[p])
            self.cov[p] = sparse.coo_matrix( (c*numpy.sqrt(self.var[i,p]*self.var[j,p]), (i,j)),
                                             shape=self.shape[:-1]).tocsr()
        self.is_correlation = False

    ####################################################################
    # TODO: Generalize these to fill and thin or something
    def bin_to_spaxel_covariance(self, bin_indx):
        r"""
        Propagate the covariance matrix data for the stacked spectra
        into the full cube.

        This (self) should be the covariance/correlation matrix for the
        binned spectra.

        .. todo::
            - Does this still need to be tested?
            - Generalize the nomenclature.
        
        Args:
            bin_indx (`numpy.ndarray`_):
                The integer vector with the bin associated with each
                spectrum in the DRP cube. This is the flattened BINID
                array.

        Returns:
            :class:`Covariance`: Covariance/Correlation matrix for
            the spaxelized binned data.
        """
        # Total number of spectra
        nspec = len(bin_indx)

        # Get the unique bins and how to reconstruct the bins from the
        # unique set
        unique_bins, reconstruct = numpy.unique(bin_indx, return_inverse=True)

        # Need to handle missing bins
        if unique_bins.size != unique_bins[-1]+2 \
                or numpy.any((unique_bins - numpy.arange(-1,unique_bins.size-1)) != 0):
            warnings.warn('Bin numbers and indices do not match.  Spectra are expected '
                          'to be sorted by their bin number.')
            unique_bins = numpy.arange(-1,unique_bins.size-1)

        # Get the valid bins
        indx = bin_indx > -1

        # Expand the covariance matrix by repeating the full matrix
        # elements for repeated bin values in different spaxels
        nchan = self.shape[-1]
        spaxel_covar = numpy.empty(nchan, dtype=sparse.csr.csr_matrix)
        for i in range(nchan):
            j = self.input_indx[i]

            # Input bin and covariance indices
            ii = unique_bins[reconstruct[indx]]
            ii_i, ii_j = map( lambda x: x.ravel(), numpy.meshgrid(ii, ii) )

            # Output spaxel and covariance indices
            oi = numpy.arange(nspec)[indx]
            oi_i, oi_j = map( lambda x: x.ravel(), numpy.meshgrid(oi, oi) )

            _covar = numpy.zeros((nspec, nspec), dtype=numpy.float)
            _covar[oi_i, oi_j] = self.toarray(channel=j)[ii_i,ii_j]
            spaxel_covar[i] = sparse.triu(_covar).tocsr()

        return Covariance(spaxel_covar, input_indx=self.input_indx)

    def spaxel_to_bin_covariance(self, bin_indx):
        r"""
        Opposite of :func:`bin_to_spaxel_covariance`: Revert the covariance
        matrix to the covariance between the unique binned spectra.

        This (self) should be the covariance/correlation matrix for the
        binned spectra redistributed to the size of the original spaxels
        map.

        .. todo::
            - Does this still need to be tested?

        .. warning::

            **This does NOT propagate the covariance between the spaxels
            into the covariance in the binned data.** That operation is
            done by, e.g.,
            :func:`mangadap.proc.spectralstack.SpectralStack._stack_with_covariance`.
            This **only** performs the inverse operation of
            :func:`bin_to_spaxel_covariance`.

        Args:
            bin_indx (numpy.ndarray): The integer vector with the bin
                associated with each spectrum in the DRP cube.  This is
                the flattened BINID array.

        Returns:
            :class:`Covariance`: Covariance/Correlation matrix for
            the stacked spectra.
        """
        # Get the unique bins and their first occurrence in the bin list
        unique_bins, unique_indx = numpy.unique(bin_indx, return_index=True)

        # Need to handle missing bins
        if unique_bins.size != unique_bins[-1]+2 \
                or numpy.any((unique_bins - numpy.arange(-1,unique_bins.size-1)) != 0):
            warnings.warn('Bin numbers and indices do not match.  Spectra are expected '
                          'to be sorted by their bin number.')
            unique_bins = numpy.arange(-1,unique_bins.size-1)

        # Total number of bins
        nbins = len(unique_bins)-1

        nchan = self.shape[-1]
        bin_covar = numpy.empty(nchan, dtype=sparse.csr.csr_matrix)
        for i in range(nchan):
            j = self.input_indx[i]

            # Input spectrum and covariance indices
            ii = unique_indx[1:]
            ii_i, ii_j = map( lambda x: x.ravel(), numpy.meshgrid(ii, ii) )

            # Output spectrum and covariance indices
            oi = unique_bins[1:]
            oi_i, oi_j = map( lambda x: x.ravel(), numpy.meshgrid(oi, oi) )

            _covar = numpy.zeros((nbins, nbins), dtype=numpy.float)
            _covar[oi_i, oi_j] = self.toarray(channel=j)[ii_i,ii_j]
            bin_covar[i] = sparse.triu(_covar).tocsr()
        return Covariance(bin_covar, input_indx=self.input_indx)

