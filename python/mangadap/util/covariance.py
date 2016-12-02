# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
r"""
Defines a class used to store and interface with covariance matrices.

*License*:
    Copyright (c) 2015, SDSS-IV/MaNGA Pipeline Group
        Licensed under BSD 3-clause license - see LICENSE.rst

*Source location*:
    $MANGADAP_DIR/python/mangadap/util/covariance.py

*Imports and python version compliance*:
    ::

        from __future__ import division
        from __future__ import print_function
        from __future__ import absolute_import
        from __future__ import unicode_literals
    
        import sys
        if sys.version > '3':
            long = int
    
        import numpy
        from scipy import sparse
        from astropy.io import fits
        from matplotlib import pyplot

*Usage examples*:
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
    binary table.  Since the covariance matrix is symmetric by
    definition, only those non-zero elements in the upper triangle
    (:math:`C_{ij}` where :math:`i\leq j`) are saved (in memory or on
    disk).  You can read an existing covariance matrix fits file::

        from mangadap.util.covariance import Covariance
        C = Covariance(ifile='test_covariance.fits')

    You can calculate a set of covariance matrices or the full
    covariance cube::

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
        print(C[1000,0,0])

        # Write the full cube, or set of channels
        CC.write('full_covariance_cube.fits')   # BEWARE: This will be a BIG file
        C.write('covariance_channel_set.fits')

    Although you can access the data in the covariance matrix as
    explained above, this is generally inefficient.  You're likely
    better of working with the :attr:`Covariance.cov` attribute
    directly.  To do this, you won't be able to use the aliasing of the
    channel indices for 3D covariance matrices.

    The :class:`Covariance` class also allows you to toggle between
    accessing the matrix as a true covariance matrix, or by splitting
    the matrix into its variance and correlation components.  For
    example,::

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
    can be written and read in without issue.  See
    :func:`Covariance.write` and :func:`Covariance.read`.  For example::

        # Get the covariance matrix for a single wavelength channel
        import numpy
        from mangadap.util.covariance import Covariance
        from mangadap.drpfits import DRPFits
        drpf = DRPFits(7495, 12703, 'RSS', read=True)
        channels = [ 0, 1000, 2000, 3000, 4000 ]
        Cov = drpf.covariance_cube(channels=channels)

        Cov.to_correlation()
        Cov.show(plane=2000)
        Cov.write('correlation_matrix.fits', clobber=True)
        Cov.revert_correlation()

        Corr = Covariance(ifile='correlation_matrix.fits')
        Corr.revert_correlation()

        assert ~(numpy.abs(numpy.sum(Cov.toarray(plane=2000) - Corr.toarray(plane=2000))) > 0.0)


*Revision history*:
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

.. todo::
    - Allow for calculation of the inverse of the covariance matrix.
    - Instead of 3D covariance cubes being an array of sparse objects,
      make the whole thing a sparse array.

.. _scipy.sparse.csr_matrix: http://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csr_matrix.html
.. _scipy.sparse.coo_matrix: http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.sparse.coo_matrix.html
.. _scipy.sparse.triu: http://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.triu.html
.. _matplotlib.pyplot.imshow: http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.imshow
.. _astropy.io.fits.Header: http://docs.astropy.org/en/stable/io/fits/api/headers.html#header
.. _astropy.io.fits.BinTableHDU: http://docs.astropy.org/en/stable/io/fits/api/tables.html#astropy.io.fits.BinTableHDU
.. _astropy.io.fits.Column: http://docs.astropy.org/en/stable/io/fits/api/tables.html#astropy.io.fits.Column
.. _astropy.io.fits.hdu.hdulist.HDUList: http://docs.astropy.org/en/v1.0.2/io/fits/api/hdulists.html


"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import sys
if sys.version > '3':
    long = int

import numpy
import warnings

from scipy import sparse
from astropy.io import fits
from matplotlib import pyplot

__author__ = 'Kyle B. Westfall'

class Covariance:
    r"""
    A general utility for storing, manipulating, and file I/O of large
    but sparse covariance matrices.

    Works under the assumption that covariance matrices are symmetric by
    definition.

    Assumes all sparse matrices in the input ndarray have the same size.

    Args:
        inp (`scipy.sparse.csr_matrix`_ or numpy.ndarray): (**Optional**)
           Covariance matrix to store.  Data type can be either a single
           `scipy.sparse.csr_matrix`_ object or a 1-D array of them.
           If not provided, the covariance object is instantiated empty.
        input_indx (numpy.ndarray): (**Optional**) If *inp* is an array of
            `scipy.sparse.csr_matrix`_ objects, this is an integer
            array specifying a pseudo-set of indices to use instead of
            the direct index in the array.  I.e., if *inp* is an array
            with 5 elements, one can provide a 5-element array for
            :attr:`input_indx` that are used instead as the designated
            index for the covariance matrices.  See:
            :func:`__getitem__`, :func:`_grab_true_index`.
        impose_triu (bool): (**Optional**) Flag to force the
            `scipy.sparse.csr_matrix`_ object to only be the upper
            triangle of the covariance matrix.  The covariance matrix is
            symmetric such that C_ij = C_ji, so it's not necessary to
            keep both values.  This flag will force a call to
            `scipy.sparse.triu`_ when setting the covariance matrix.
            Otherwise, the input matrix is **assumed** to only have the
            upper triangle of numbers.
        source (str or `astropy.io.fits.hdu.hdulist.HDUList`_):
            (**Optional**) Initialize the object using an
            `astropy.io.fits.hdu.hdulist.HDUList`_ object or path to a
            fits file.  See :func:`read`.
        primary_ext (str): (**Optional**) If reading the data from
            *source*, this is the name of the extension with the header
            information needed to construct the :class:`Covariance`
            object.  Default is 'PRIMARY'.  See :func:`read`.
        covar_ext (str): (**Optional**) If reading the data from
            *source*, this is the name of the extension with covariance
            data.  Default is 'COVAR'.  See :func:`read`.
        plane_ext (str): (**Optional**) If reading the data from
            *source*, this is the name of the extension with the
            covariance plane indices, if necessary.  Default is 'PLANE'.
            See :func:`read`.

    Raises:
        TypeError: Raised if the input array is not one-dimensional or
            the input covariance matrix are not
            `scipy.sparse.csr_matrix`_ objects.
        Exception: Raised if the :attr:`input_indx` either does not have
            the correct size or is not required due to the covariance
            object being a single matrix.

    Attributes:
        cov (`scipy.sparse.csr_matrix`_ or numpy.ndarray): The
            covariance matrix stored in sparse format.
        dim (int): The number of dimensions in the covariance matrix.
            This can either be 2 (a single covariance matrix) or 3
            (multiple covariance matrices).
        nnz (int): The number of non-zero covariance matrix elements.
        input_indx (numpy.ndarray): If :attr:`cov` is an array of
            `scipy.sparse.csr_matrix`_ objects, this is an integer
            array specifying a pseudo-set of indices to use instead of
            the direct index in the array.  I.e., if :attr:`cov` is an
            array with 5 elements, one can provide a 5-element array for
            :attr:`input_indx` that are used instead as the designated
            index for the covariance matrices.  See:
            :func:`__getitem__`, :func:`_grab_true_index`.
        inv (`scipy.sparse.csr_matrix`_ or numpy.ndarray): The inverse
            of the covariance matrix.  **This is not currently
            calculated!**
        var (numpy.ndarray): Array with the variance provided by the
            diagonal of the/each covariance matrix.  This is only
            populated if necessary, either by being requested
            (:func:`variance`) or if needed to convert between
            covariance and correlation matrices.
        is_correlation (bool): Flag that the covariance matrix has been
            saved as a variance vector and a correlation matrix.

    """
    def __init__(self, inp=None, input_indx=None, impose_triu=False, source=None,
                 primary_ext='PRIMARY', covar_ext='COVAR', plane_ext='PLANE'):

        # If a file is provided with the covariance matrix, read it and
        # return
        if source is not None:
            self.read(source, impose_triu=impose_triu, primary_ext=primary_ext,
                      covar_ext=covar_ext, plane_ext=plane_ext)
            return

        # If no input is provided, free any existing data and return
        if inp is None:
            self._free()
            return
        # Otherwise, save the input covariance matrix
        else:
            self.cov = inp

        # Set the dimensionality, check that each element of the array
        # has the correct type, count the number of non-zero covariance
        # values, and (if necessary) initialize the indices for each
        # covariance matrix
        self.input_indx = None
        if isinstance(self.cov, numpy.ndarray):
            if len(self.cov.shape) > 1:
                raise TypeError('Input ndarray can only be one-dimensional')
            self.dim = 3
            self.nnz = 0
            for cov in self.cov:
                if not sparse.isspmatrix_csr(cov):
                    self._free()
                    raise TypeError('Input covariance matrix (or elements) must by csr_matrices.')
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
                self._free()
                raise TypeError('Input covariance matrix (or elements) must by csr_matrices.')
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
        # **NOT YET IMPLEMENTED**
        self.inv = None

        # Set the variance array and the correlation matrix flag
        self.var = None
        self.is_correlation = False

    
    def __getitem__(self, *args):
        """
        Return the covariance value at a provided 2D or 3D position.

        Args:
            *args (pointer): 2 or 3 integers designating the covariance
                value to return.  Number of values must match the
                dimensionality of the object.

        Returns:
            float: The value of the covariance matrix at the designated
            index.

        Raises:
            IndexError: Raised if the number of arguments does not match
                the dimensionality of the object.
        """
        indx = tuple(*args)
        if len(indx) != self.dim:
            raise IndexError('Incorrect number of dimensions!')
        if self.dim == 2:
            return self.cov[tuple(sorted(indx))]
        return self.cov[self._grab_true_index(indx[0])][tuple(sorted(indx[1:]))]


    def __del__(self):
        """Deconstructor, ensuring memory storage freed."""
        self._free()


    def _grab_true_index(self, inp):
        """
        In the case of a 3D array, return the true array-element index
        given the pseudo index.

        Args:
            inp (int): Requested index to convert to the real index in
                the stored array of covariance matrices.

        Returns:
            int: The index in the stored array of the requested
            covariance matrix.

        .. todo::
            Search operation is inefficient.  Apparently a better option
            is a feature request for numpy 2.0
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
        self.shape += self.cov.ravel()[0].shape


    def _impose_upper_triangle(self):
        """
        Force :attr:`cov` to only contain non-zero elements in its upper
        triangle.
        """
        if self.dim == 2:
            self.cov = sparse.triu(self.cov).tocsr()
            self.nnz = self.cov.nnz
            return

        self.nnz = 0
        for i in range(self.shape[0]):
            self.cov[i] = sparse.triu(self.cov[i]).tocsr()
            self.nnz += self.cov[i].nnz


    def _with_lower_triangle(self, plane=None):
        """
        Return a `scipy.sparse.csr_matrix`_ object with both its
        upper and lower triangle filled, ensuring that they are
        symmetric.

        Args:
            plane (int): (**Optional**) The pseudo-index of the
                covariance matrix to return.

        Returns:
            `scipy.sparse.csr_matrix`_: The sparse matrix with both the
            upper and lower triangles filled (with symmetric
            information).

        Raises:
            ValueError: Raised if the object is 3D and *plane* is not
                provided.
        """
        if self.dim == 2:
            a = self.cov
        else:
            if plane is None:
                raise ValueError('Must define plane!  Use plane=...')
            a = self.cov[self._grab_true_index(plane)]

        return (sparse.triu(a) + sparse.triu(a,1).T)


    def _free(self):
        """Free the memory allocated to the object."""
        self.dim=0
        self.shape=None

        try:
            del self.cov
            self.cov = None
        except (NameError,AttributeError):
            self.cov = None

        try:
            del self.inv
            self.inv = None
        except (NameError,AttributeError):
            self.inv = None


    def show(self, plane=None, zoom=None, ofile=None):
        """
        Convert the (selected) covariance matrix to a filled array and
        plot the array using `matplotlib.pyplot.imshow`_.  If an output
        file is provided, the image is redirected to the designated
        output file; otherwise, the image is plotted to the screen.

        Args:
            plane (int): (**Optional**) The pseudo-index of the
                covariance matrix to plot.  Required if the covariance
                object is 3D.
            zoom (float): (**Optional**) Factor by which to zoom in on
                the center of the image by *removing the other regions
                of the array*.  E.g. *zoom=2* will show only the central
                quarter of the covariance matrix.
            ofile (str): (**Optional**) If provided, the array is output
                to this file instead of being plotted to the screen.
        """
        # Convert the covariance matrix to an array
        a = self.toarray(plane)

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
            im = pyplot.imshow(a, interpolation='nearest', origin='lower')
            pyplot.colorbar()
            pyplot.show()
            return

        # Print the plot to the designated output file
        fig = pyplot.figure(1)
        im = pyplot.imshow(a, interpolation='nearest', origin='lower')
        pyplot.colorbar()
        fig.canvas.print_figure(ofile)
        pyplot.show()
        

    def toarray(self, plane=None):
        """
        Convert the covariance to a full array, filled with zeros when
        appropriate.

        Args:
            plane (int): (**Optional**) The pseudo-index of the
                covariance matrix to plot.  Required if the covariance
                object is 3D.

        Returns:
            numpy.ndarray: Dense array with the full covariance matrix.
        """
        if self.dim == 2 or plane is not None:
            return (self._with_lower_triangle(plane=plane)).toarray()
        
        arr = numpy.empty(self.shape, dtype=numpy.float)
        for i in range(self.shape[0]):
            indx = i if self.input_indx is None else self.input_indx[i]
            arr[i,:,:] = (self._with_lower_triangle(plane=indx)).toarray()
        return arr
        

    def find(self, plane=None):
        """
        Find the non-zero values in the **full** covariance matrix (not
        just the upper triangle).

        Args:
            plane (int): (**Optional**) The pseudo-index of the
                covariance matrix to plot.  Required if the covariance
                object is 3D.

        Returns:
            tuple: A tuple of arrays *i*, *j*, and *c*.  The arrays *i*
            and *j* contain the index coordinates of the non-zero
            values, and *c* contains the values themselves.
        """
        return sparse.find(self._with_lower_triangle(plane=plane))


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


    def binary_columns(self, hdr=None):
        r"""
        Construct the binary columns for the output fits file.

        The four columns are ``INDX``, ``COVAR``, ``VARIANCE``, and
        ``INP_INDX``; see :func:`write`.

        Args:
            hdr (`astropy.io.fits.Header`_) : (**Optional**) A header
                object that, if provided, will have the keywords
                ``COVSHAPE`` and ``COVTYPE`` added based on,
                respectively, the values of :attr:`shape` and
                :attr:`is_correlation`.

        Returns:
            `astropy.io.fits.Column`_: Up to four
            `astropy.io.fits.Column`_ objects with the column data, in
            the sequence listed above.  If the data is not in the form
            of a correlation matrix, the ``VARIANCE`` column is returned
            as None.  If the covariance matrix is only two-dimensional,
            the ``INP_INDX`` column is returned as None.
        """
        # Add the shape to the header
        if hdr is not None:
            hdr['COVSHAPE'] = (str(self.shape), 'Shape of the covariance matrix')
            hdr['COVTYPE'] = ('Correlation' if self.is_correlation else 'Covariance',
                              'Type of covariance data storage')

        # Create the binary table data
        coo = numpy.empty((self.nnz, self.dim), dtype=numpy.int32)
        coo_form = str(self.dim)+'J'
        if self.dim == 2:
            ii, jj, covar_value = sparse.find(self.cov)
            coo[:,0] = ii
            coo[:,1] = jj
            if self.is_correlation:
                var_value = numpy.zeros(self.nnz, dtype=numpy.float64)
                indx = ii==jj
                var_value[indx] = self.var[ii[indx]]
        else:
            covar_value = numpy.empty(self.nnz, dtype=numpy.float64)
            var_value = numpy.zeros(self.nnz, dtype=numpy.float64)
            i = 0
            j = 0
            for cov in self.cov:
                kk = numpy.arange(cov.nnz)+j
                coo[kk,0] = i
                ii, jj, vv = sparse.find(cov)
                coo[kk,1] = ii
                coo[kk,2] = jj
                covar_value[kk] = vv
                if self.is_correlation:
                    indx = ii==jj
                    var_value[kk[indx]] = self.var[i,ii[indx]]
                j += cov.nnz
                i += 1

        return fits.Column(name='INDX', format=coo_form, array=coo), \
               fits.Column(name='COVAR', format='1D', array=covar_value), \
               (fits.Column(name='VARIANCE', format='1D', array=var_value) \
                    if self.is_correlation else None), \
               (None if self.input_indx is None else \
                    fits.Column(name='INP_INDX', format='1J', array=self.input_indx))

        
    def write(self, ofile, hdr=None, clobber=False):
        r"""
        Write the covariance object to a fits file such that it can be
        read for later use; see :func:`read`.  The covariance matrix
        (matrices) are stored in "coordinate" format using fits binary
        tables; see `scipy.sparse.coo_matrix`_.

        Independent of the dimensionality of the covariance matrix, the
        written file has a ``PRIMARY`` extension with two keywords:

            - ``COVSHAPE``: Specifies the original dimensions of the
              covariance matrix; see :attr:`shape`.
            - ``COVTYPE``: Designates if the covariance matrix is of
              type ``'Covariance'`` or ``'Correlation'``; see
              :attr:`is_correlation`.
            
        The covariance data itself is written to the ``COVAR``
        extension, which has three columns:

            - ``INDX``: The indices of the non-zero elements of the
              covariance matrix.  For 2D matrices, this means each
              column element is a two-element vector, :math:`(i,j)`; for
              3D matrices, it is a three-element vector,
              :math:`(i,j,k)`.
            - ``COVAR``: The value of the covariance at the pixel
              coordinate provided by ``INDX``; e.g., :math:`C_{ij}` or
              :math:`C_{ijk}`.  If the :class:`Covariance` object being
              written is a correlation matrix, these are the values of
              :math:`\rho_{ij}` or :math:`\rho_{ijk}`; along the
              diagonal of matrices, the values of :math:`\rho=1` are
              also stored.
            - ``VARIANCE``: If the :class:`Covariance` object being
              written is a correlation matrix, this column provides the
              variances at pixel coordinate :math:`(i,j)`,
              :math:`(i,j,k)`.  To match the number of rows in the other
              columns, the off-diagonal elements are set to 0.0.

        For 3D matrices, a third extension, ``PLANE``, is written
        containing a binary table with the list of pseudo indices for
        each of the provided covariance matrices. these indices are in
        the single column in this extension, ``INP_INDX``.

        Args:
            ofile (str): File name for the output.
            hdr (`astropy.io.fits.Header`_): (**Optional**) A header
                object to include in the PRIMARY extension.  The SHAPE
                keyword will be added/overwritten.
            clobber (bool): (**Optional**) Overwrite any existing file.

        Raises:
            TypeError: Raise if the input *hdr* does not have the
                correct type.
        """
        # Use input header or create a minimal one
        if hdr is None:
            hdr = fits.Header()
        # Ensure the input header has the correct type
        elif not isinstance(hdr, fits.Header):
            raise TypeError('Input header must have type astropy.io.fits.Header.')

        # Get the fits columns
        indx_col, covar_col, var_col, plane_col = self.binary_columns(hdr=hdr)

        # Construct the main binary table HDU
        covar_hdu = fits.BinTableHDU.from_columns([ indx_col, covar_col ], name='COVAR') \
                        if var_col is None else \
                    fits.BinTableHDU.from_columns([ indx_col, var_col, covar_col ], name='COVAR')

        # Construct the HDUList
        hdu = fits.HDUList([ fits.PrimaryHDU(header=hdr), covar_hdu ])

        # Add the index column, if necessary
        if plane_col is not None:
            hdu.append(fits.BinTableHDU.from_columns([ plane_col ], name='PLANE'))

        # Write the file
        hdu.writeto(ofile, clobber=clobber, checksum=True)


    def read(self, source, primary_ext='PRIMARY', covar_ext='COVAR', plane_ext='PLANE',
             impose_triu=False, to_correlation=None, return_hdr=False, quiet=False):
        r"""
        Read an existing covariance object previously written to disk;
        see :func:`write`.  The class can read covariance data written
        by other programs *as long as they have a commensurate format*.
        
        If the extension names and column names are correct,
        :class:`Covariance` can read fits files that were not
        necessarily produced by the class itself.  This is useful for
        MaNGA DAP products that include the covariance data in among
        other output products.  Because of this,
        :class:`Covariance.binary_columns` is included to provide the
        columns and header data that can be placed in any fits file.

        Args:
            source (str or `astropy.io.fits.hdu.hdulist.HDUList`_):
                `astropy.io.fits.hdu.hdulist.HDUList`_ object or path to
                fits file with the covariance data.
            primary_ext (str): (**Optional**) Name of the extension with
                the header information needed to construct the
                :class:`Covariance` object.  Default is 'PRIMARY'.
            covar_ext (str): (**Optional**) Name of the extension with
                covariance data.  Default is 'COVAR'.
            plane_ext (str): (**Optional**) Name of the extension with
                the covariance plane indices, if necessary.  Default is
                'PLANE'.
            impose_triu (bool): (**Optional**) Flag to force the
                `scipy.sparse.csr_matrix`_ object(s) to only contain
                elements in their upper triangle.  The covariance matrix
                is symmetric such that C_ij = C_ji, so it's not
                necessary to keep both values.  This flag will force a
                call to `scipy.sparse.triu`_ when setting the covariance
                matrix.  Otherwise, the input matrix is **assumed** to
                only have the upper triangle of numbers.
            to_correlation (bool): (**Optional**) Force the operation to
                return a correlation matrix regardless of the type of
                data in the file.
            return_hdr (bool): (**Optional**) Return the
                `astropy.io.fits.Header`_ object read from the primary
                extension of the file.
            quiet (bool): (**Optional**) Suppress screen output

        Returns:
            `astropy.io.fits.Header`_: The header from the primary
            extension, if requested.
        """

        if isinstance(source, fits.HDUList):
            hdu = source
        else:
            try:
                hdu = fits.open(source)
            except Exception as e:
                print(e)
                warnings.warn('Problem reading covariance from file.')
                return

        # Erase anything that might already exist
        self._free()

        # Get the shape and type
        self.shape = eval(hdu[primary_ext].header['COVSHAPE'])
        self.dim = len(self.shape)
        self.nnz = hdu[covar_ext].header['NAXIS2']
        self.is_correlation = hdu[primary_ext].header['COVTYPE'] == 'Correlation'
        self.var = None

        # Read a 2D matrix
        if self.dim == 2:
            self.cov = sparse.coo_matrix( (hdu[covar_ext].data['COVAR'], \
                                          (hdu[covar_ext].data['INDX'][:,0], \
                                          hdu[covar_ext].data['INDX'][:,1])), \
                                          shape=self.shape).tocsr()
            # Read the variances if written as a correlation matrix
            if self.is_correlation:
                self.var = numpy.zeros(self.shape[0], dtype=numpy.float)
                indx = hdu[covar_ext].data['INDX'][:,0] == hdu[covar_ext].data['INDX'][:,1]
                self.var[hdu[covar_ext].data['INDX'][indx,0]] = \
                        hdu[covar_ext].data['VARIANCE'][indx]
            self.input_indx = None
        # Read a 3D matrix
        else:
            # Initialize the arrays
            self.cov = numpy.empty(self.shape[0], dtype=sparse.csr.csr_matrix)
            if self.is_correlation:
                self.var = numpy.zeros((self.shape[0], self.shape[1]), dtype=numpy.float)
            # Read each plane
            for i in range(self.shape[0]):
                tt = numpy.where( hdu[covar_ext].data['INDX'][:,0].flatten() == i )[0]

                ii = hdu[covar_ext].data['INDX'][tt,1].flatten()
                jj = hdu[covar_ext].data['INDX'][tt,2].flatten()
                vv = hdu[covar_ext].data['COVAR'][tt]

                self.cov[i] = sparse.coo_matrix( (vv, (ii, jj)), shape=self.shape[1:]).tocsr()

                # Read the variances if written as a correlation matrix
                if self.is_correlation:
                    indx = ii == jj
                    self.var[i,ii[indx]] = hdu[covar_ext].data['VARIANCE'][tt[indx]]

            # Try to assign the pseudo-indices
            try:
                self.input_indx = hdu[plane_ext].data['INP_INDX']
            except KeyError:
                self.input_indx = numpy.arange(self.shape[0])

        # Force the covariance matrix to only have the upper triangle
        if impose_triu:
            self._impose_upper_triangle()

        self.inv = None

        # Convert to a correlation matrix if requested
        if to_correlation is not None:
            if to_correlation:
                self.to_correlation()
            else:
                self.revert_correlation()

        # Report
        # TODO: Convert report to use logging
        if not quiet:
            print('Read covariance cube:')
            print('        input type: {0}'.format(hdu[primary_ext].header['COVTYPE']))
            print('       output type: {0}'.format('Correlation' \
                                                    if self.is_correlation else 'Covariance'))
            print('        dimensions: {0}'.format(self.dim))
            print('             shape: {0}'.format(self.shape))
            if self.dim == 3:
                print('    pseudo-indices: ', self.input_indx)
            print('   non-zero values: {0}'.format(self.nnz))

        # Return the header object if requested
        if return_hdr:
            return hdu[primary_ext].header

    
    def variance(self):
        """
        If not already done, grab the variances along the diagonal of
        the covariance matrix.  The function returns the variance for
        all planes if more than one exists.
        """
        if self.var is not None:
            return self.var

        if self.dim == 2:
            self.var = numpy.diag(self.cov.toarray()).copy()
            return self.var

        self.var = numpy.empty((self.shape[0], self.shape[1]), dtype=numpy.float)
        for p in range(self.shape[0]):
            self.var[p,:] = numpy.diag(self.cov[p].toarray()).copy()

        return self.var


    def to_correlation(self):
        r"""
        Convert the covariance matrix into a correlation matrix by
        dividing each element by the variances.
        
        If the matrix is a correlation matrix already (see
        :attr:`is_correlation`), no operations are performed.
        Otherwise, the variance vectors are computed, if necessary, and
        used to normalize the covariance values.

        A :class:`Covariance` object can be reverted from a correlation
        matrix; see :func:`revert_correlation`.
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
       
        for p in range(self.shape[0]):
            i, j, c = sparse.find(self.cov[p])
            self.cov[p] = sparse.coo_matrix( (c/numpy.sqrt(self.var[p,i]*self.var[p,j]), (i,j)),
                                             shape=self.shape[1:]).tocsr()


    def revert_correlation(self):
        r"""
        Revert the object from a correlation matrix back to a full
        covariance matrix.

        .. todo::
            Allow this function to take in a new set of variances to
            use.

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
       
        for p in range(self.shape[0]):
            i, j, c = sparse.find(self.cov[p])
            self.cov[p] = sparse.coo_matrix( (c*numpy.sqrt(self.var[p,i]*self.var[p,j]), (i,j)),
                                             shape=self.shape[1:]).tocsr()
        self.is_correlation = False

        




