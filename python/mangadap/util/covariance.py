# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Defines a class used to store and interface with covariance matrices.

*License*:
    Copyright (c) 2015, Kyle B. Westfall
    Licensed under BSD 3-clause license - see LICENSE.rst

*Source location*:
    $MANGADAP_DIR/python/mangadap/util/covariance.py

*Python2/3 compliance*::

    from __future__ import division
    from __future__ import print_function
    from __future__ import absolute_import
    from __future__ import unicode_literals
    
    import sys
    if sys.version > '3':
        long = int

*Imports*::

    import os.path
    from scipy import sparse
    import scipy.linalg
    from astropy.io import fits
    import time
    import numpy
    from matplotlib import pyplot

*Class usage examples*:

    .. todo::

        Add some usage comments here!

*Revision history*:
    | **23 Feb 2015**: Original Implementation by K. Westfall (KBW)
    | **04 Aug 2015**: (KBW) Sphinx documentation and minor edits.

.. todo::
    - Allow for calculation of the inverse of the covariance matrix.

.. _scipy.sparse.csr_matrix: http://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csr_matrix.html
.. _scipy.sparse.triu: http://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.triu.html
.. _matplotlib.pyplot.imshow: http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.imshow
.. _astropy.io.fits.Header: http://docs.astropy.org/en/stable/io/fits/api/headers.html#header
"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import sys
if sys.version > '3':
    long = int

import os.path
from scipy import sparse
import scipy.linalg
from astropy.io import fits
import time
import numpy
from matplotlib import pyplot

__author__ = 'Kyle B. Westfall'

class Covariance:
    """
    Meant to be a general utility for storing, manipulating, and file
    I/O of large but sparse covariance matrices.

    Works under the assumption that covariance matrices are symmetric by
    definition.

    Assumes all sparse matrices in the input ndarray have the same size.

    Args:
        inp (`scipy.sparse.csr_matrix`_ or numpy.array): (Optional)
           Covariance matrix to store.  Data type can be either a single
           `scipy.sparse.csr_matrix`_ object or a 1-D array of them.
           If not provided, the covariance object is instantiated empty.
        input_indx (numpy.array): (Optional) If *inp* is an array of
            `scipy.sparse.csr_matrix`_ objects, this is an integer
            array specifying a pseudo-set of indices to use instead of
            the direct index in the array.  I.e., if *inp* is an array
            with 5 elements, one can provide a 5-element array for
            :attr:`input_indx` that are used instead as the designated
            index for the covariance matrices.  See:
            :func:`__getitem__`, :func:`_grab_true_index`.
        impose_triu (bool): (Optional) Flag to force the
            `scipy.sparse.csr_matrix`_ object to only be the upper
            triangle of the covariance matrix.  The covariance matrix is
            symmetric such that C_ij = C_ji, so it's not necessary to
            keep both values.  This flag will force a call to
            `scipy.sparse.triu`_ when setting the covariance matrix.
            Otherwise, the input matrix is **assumed** to only have the
            upper triangle of numbers.
        ifile (str): File from which to read the covariance matrix.
            See: :func:`write`.

    Raises:
        TypeError: Raised if the input array is not one-dimensional or
            the input covariance matrix are not
            `scipy.sparse.csr_matrix`_ objects.
        Exception: Raised if the :attr:`input_indx` either does not have
            the correct size or is not required due to the covariance
            object being a single matrix.

    Attributes:
        cov (`scipy.sparse.csr_matrix`_ or numpy.array): The
            covariance matrix stored in sparse format.
        dim (int): The number of dimensions in the covariance matrix.
            This can either be 2 (a single covariance matrix) or 3
            (multiple covariance matrices).
        nnz (int): The number of non-zero covariance matrix elements.
        input_indx (numpy.array): If :attr:`cov` is an array of
            `scipy.sparse.csr_matrix`_ objects, this is an integer
            array specifying a pseudo-set of indices to use instead of
            the direct index in the array.  I.e., if :attr:`cov` is an
            array with 5 elements, one can provide a 5-element array for
            :attr:`input_indx` that are used instead as the designated
            index for the covariance matrices.  See:
            :func:`__getitem__`, :func:`_grab_true_index`.
        inv (`scipy.sparse.csr_matrix`_ or numpy.array): The inverse
            of the covariance matrix.  **This is not currently
            calculated!**
    """
    def __init__(self, inp=None, input_indx=None, impose_triu=False, ifile=None):

        # If a file is provided with the covariance matrix, read it and
        # return
        if ifile is not None:
            self.read(ifile, impose_triu=impose_triu)
            return

        # If no input is provided, free any existing data and return
        if inp is None:
            self._free()
            return
        # Otherwise, save the input covariance matrix
        else:
            self.cov = inp

        # Set the dimensionality, check that each element of the array
        # has the correct type, and count the number of non-zero
        # covariance values
        if type(self.cov) == numpy.ndarray:
            if len(self.cov.shape) > 1:
                raise TypeError('Input ndarray can only be one-dimensional')
            self.dim = 3
            self.nnz = 0
            for cov in self.cov:
                if not sparse.isspmatrix_csr(cov):
                    self._free()
                    raise TypeError('Input covariance matrix (or elements) must by csr_matrices.')
                self.nnz += cov.nnz
        else:
            self.dim = 2
            if not sparse.isspmatrix_csr(self.cov):
                self._free()
                raise TypeError('Input covariance matrix (or elements) must by csr_matrices.')
            self.nnz = self.cov.nnz

        # Set the shape of the full matrix/cube
        self._set_shape()

        # Initialize the indices for each covariance matrix, if provided
        # and viable
        self.input_indx = None
        if input_indx is not None:
            if self.dim == 2:
                raise Exception('Input indices only allowed when allocating multiple matrices.')
            if input_indx.shape != self.cov.shape:
                raise Exception('Input array and input index array must be the same size.')
            self.input_indx = input_indx
#           print(self.input_indx)

        # If requested, impose that the input matrix only have values in
        # its upper triangle.
        if impose_triu:
            self._impose_upper_triangle()

        # Set the inverse of the covariance matrix
        # **NOT YET IMPLEMENTED**
        self.inv = None

    
    def __getitem__(self, *args):
        """
        Return the covariance value at a provided 2D or 3D position.

        Args:
            *args (pointer): 2 or 3 integers designating the covariance
                value to return.  Number of values must match the
                dimensionality of the object.

        Returns:
            float : The value of the covariance matrix at the designated
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
        return self.cov[self._grab_true_index(indx[:1])][tuple(sorted(indx[1:]))]


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
            int : The index in the stored array of the requested
                covariance matrix.
        """
        if self.input_indx is not None:
            # Find the matching index
            # TODO: Better way to do this?
            i=0
            for coo in self.input_indx:
                if coo == inp:
                    return i
                i += 1
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
        for i in range(0,self.shape[0]):
            self.cov[i] = sparse.triu(self.cov[i]).tocsr()
            self.nnz += self.cov[i].nnz


    def _with_lower_triangle(self, plane=None):
        """
        Return a `scipy.sparse.csr_matrix`_ object with both its
        upper and lower triangle filled, ensuring that they are
        symmetric.

        Args:
            plane (int): (Optional) The pseudo-index of the covariance
                matrix to return.

        Returns:
            `scipy.sparse.csr_matrix`_ : The sparse matrix with both
                the upper and lower triangles filled (with symmetric
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
        plot the array using `matplotlib.pyplot.imshow`.  If an output
        file is provided, the image is redirected to the designated
        output file; otherwise, the image is plotted to the screen.

        Args:
            plane (int): (Optional) The pseudo-index of the covariance
                matrix to plot.  Required if the covariance object is
                3D.
            zoom (float): (Optional) Factor by which to zoom in on the
                center of the image by *removing the other regions of
                the array*.  E.g. *zoom=2* will show only the central
                quarter of the covariance matrix.
            ofile (str): (Optional) If provided, the array is output to
                this file instead of being plotted to the screen.
        """
        # Convert the covariance matrix to an array
        a = self.toarray(plane)

        # Remove some fraction of the array to effectively zoom in on
        # the center of the covariance matrix
        if zoom is not None:
            xs = self.shape[self.dim-2]/2 - self.shape[self.dim-2]/2/zoom
            xe = xs + self.shape[self.dim-2]/zoom + 1
            ys = self.shape[self.dim-1]/2 - self.shape[self.dim-1]/2/zoom
            ye = ys + self.shape[self.dim-1]/zoom + 1
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
            plane (int): (Optional) The pseudo-index of the covariance
                matrix to plot.  Required if the covariance object is
                3D.

        Returns:
            numpy.array : Dense array with the full covariance matrix.
        """
        return (self._with_lower_triangle(plane)).toarray()
        

    def find(self, plane=None):
        """
        Find the non-zero values in the **full** covariance matrix (not
        just the upper triangle).

        Args:
            plane (int): (Optional) The pseudo-index of the covariance
                matrix to plot.  Required if the covariance object is
                3D.

        Returns:
            tuple : A tuple of arrays *i*, *j*, and *c*.  The arrays *i*
                and *j* contain the index coordinates of the non-zero
                values, and *c* contains the values themselves.
        """
        return sparse.find(self._with_lower_triangle(plane))


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

        
    def write(self, ofile, hdr=None, clobber=False):
        r"""
        Write the covariance object to disk such that it can be read in
        later; see :func:`read`.  The covariance matrix (matrices) are
        stored in fits binary tables.  Only the non-zero values are
        stored, and they are stored in "coordinate" format.
        
        For 2D matrices:  The output fits file has two extensions,
        PRIMARY and COVAR.  The PRIMARY extension is empty apart from
        the header, which contains a keyword SHAPE specifying the
        original dimensions of the covariance matrix.  The COVAR
        extension then gives the covariance data in two columns, INDX
        and COVAR.  The INDX column provides the indices, :math:`(i,j)`,
        and the COVAR column provides the covariance values,
        :math:`C_{ij}`.

        For 3D matrices:  The output fits file has the same extensions,
        and the INDX column of the COVAR extension provides three
        indices, :math:`(i,j,k)`, for the covariance values,
        :math:`C_{ijk}`.  If :attr:`input_indx` is not None, a third
        extension, PLANE, is written containing a binary table with the
        list of pseudo indices for each of the provided covariance
        matrices; these indices are in the single column in this
        extension, INP_INDX.
        
        Args:
            ofile (str): File name for the output.
            hdr (`astropy.io.fits.Header`_): (Optional) A header object
                to include in the PRIMARY extension.  The SHAPE keyword
                will be added/overwritten.
            clobber (bool): (Optional) Overwrite any existing file.

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

        hdr['SHAPE'] = str(self.shape)                  # Add the shape to the header
        primary_hdu = fits.PrimaryHDU(header=hdr)       # Set the primary extension

        # Create the binary table data
        coo = numpy.empty((self.nnz, self.dim), dtype=numpy.int32)
        val = numpy.empty(self.nnz, dtype=numpy.float64)
        coo_form = str(self.dim)+'J'
        if self.dim == 2:
            ii, jj, vv = sparse.find(self.cov)
            coo[0:self.nnz,0] = ii
            coo[0:self.nnz,1] = jj
            val[0:self.nnz] = vv
            col1 = fits.Column(name='INDX', format=coo_form, array=coo) 
        else:
            i = 0
            j = 0
            for cov in self.cov:
                coo[j:j+cov.nnz,0] = i
                ii, jj, vv = sparse.find(cov)
                coo[j:j+cov.nnz,1] = ii
                coo[j:j+cov.nnz,2] = jj
                val[j:j+cov.nnz] = vv
                j += cov.nnz
                i += 1
            col1 = fits.Column(name='INDX', format=coo_form, array=coo)
            
        col2 = fits.Column(name='COVAR', format='1D', array=val)

        table_hdu = fits.BinTableHDU.from_columns([col1, col2], name='COVAR')
        hdu = fits.HDUList([ primary_hdu, table_hdu ])

        # Add the index column, if necessary
        if self.input_indx is not None:
            col1 = fits.Column(name='INP_INDX', format='1J', array=self.input_indx)

            indx_hdu = fits.BinTableHDU.from_columns([col1], name='PLANE')
            hdu.append(indx_hdu)

        # Write the file
        hdu.writeto(ofile, clobber=clobber)


    def read(self, ifile, impose_triu=False, return_hdr=False, quiet=False):
        r"""
        Read an existing covariance object previously written to disk;
        see :func:`write`.  A covariance matrix (or set of matrices) are
        written to disk using this object class are stored in fits
        binary tables.  Only the non-zero values are stored, and they
        are stored in "coordinate" format.  The class can read
        covariance data written by other programs *as long as they have
        a commensurate format*.
        
        For 2D matrices:  The output fits file has two extensions,
        PRIMARY and COVAR.  The PRIMARY extension is empty apart from
        the header, which contains a keyword SHAPE specifying the
        original dimensions of the covariance matrix.  The COVAR
        extension then gives the covariance data in two columns, INDX
        and COVAR.  The INDX column provides the indices, :math:`(i,j)`,
        and the COVAR column provides the covariance values,
        :math:`C_{ij}`.

        For 3D matrices:  The output fits file has the same extensions,
        and the INDX column of the COVAR extension provides three
        indices, :math:`(i,j,k)`, for the covariance values,
        :math:`C_{ijk}`.  If :attr:`input_indx` is not None, a third
        extension, PLANE, contains a binary table with the list of
        pseudo indices for each of the provided covariance matrices;
        these indices are in the single column in this extension,
        INP_INDX.
        
        Args:
            ifile (str): File name with the covariance data.

            impose_triu (bool): (Optional) Flag to force the
                `scipy.sparse.csr_matrix`_ object(s) to only contain
                elements in their upper triangle.  The covariance matrix
                is symmetric such that C_ij = C_ji, so it's not
                necessary to keep both values.  This flag will force a
                call to `scipy.sparse.triu`_ when setting the covariance
                matrix.  Otherwise, the input matrix is **assumed** to
                only have the upper triangle of numbers.

            return_hdr (bool): (Optional) Return the
                `astropy.io.fits.Header`_ object read from the primary
                extension of the file.

            quiet (bool): (Optional) Suppress screen output

        Returns:
            `astropy.io.fits.Header`_ : The header from the primary
                extension, if requested.

        """

        try:
            hdu = fits.open(ifile)
        except Exception as e:
            print('Problem reading input file; covariance object unchanged.')
            return

        # Erase anything that might already exist
        self._free()

        self.shape = eval(hdu['PRIMARY'].header['SHAPE'])
        self.dim = len(self.shape)
        self.nnz = hdu['COVAR'].header['NAXIS2']

        if self.dim == 2:
            self.cov = sparse.coo_matrix( (hdu['COVAR'].data['COVAR'], \
                                          (hdu['COVAR'].data['INDX'][:,0], \
                                          hdu['COVAR'].data['INDX'][:,1])), \
                                          shape=self.shape).tocsr()
            self.input_indx = None
        else:
            self.cov = numpy.empty(self.shape[0], dtype=sparse.csr.csr_matrix)
            for i in range(0,self.shape[0]):
                tt = numpy.where( hdu['COVAR'].data['INDX'][:,0].flatten() == i )

                ii = hdu['COVAR'].data['INDX'][tt,self.dim-2].flatten()
                jj = hdu['COVAR'].data['INDX'][tt,self.dim-1].flatten()
                vv = hdu['COVAR'].data['COVAR'][tt]

                self.cov[i] = sparse.coo_matrix( (vv, (ii, jj)), shape=self.shape[1:]).tocsr()
            try:
                self.input_indx = hdu['PLANE'].data['INP_INDX']
            except KeyError:
                self.input_indx = None

        if impose_triu:
            self._impose_upper_triangle()

        self.inv = None

        if not quiet:
            print('Read covariance cube with:')
            print('        dimensions: {0}'.format(self.dim))
            print('             shape: {0}'.format(self.shape))
            print('   non-zero values: {0}'.format(self.nnz))

        if return_hdr:
            return hdu['PRIMARY'].header

        





        

