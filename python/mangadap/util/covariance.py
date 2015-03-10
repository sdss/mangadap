from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import sys
if sys.version > '3':
    long = int

import os.path
from os import environ
from scipy import sparse
import scipy.linalg
from astropy.io import fits
import time
import numpy
from matplotlib import pyplot


__author__ = 'Kyle B. Westfall'

class covariance:
    """
    Meant to be a general utility for storing, manipulating, and file
    I/O of large but sparse covariance matrices.

    Works under the assumption that covariance matrices are symmetric by
    definition.

    Assumes all sparse matrices in the input ndarray have the same size.

    REVISION HISTORY:
        23 Feb 2015: Original Implementation by K. Westfall
    """

    def __init__(self, inp=None, input_indx=None, impose_triu=False, ifile=None):
        """
        ARGUMENTS:
            inp 
                Input matrix.  Can be either a single
                scipy.sparse.csr.csr_matrix or a 1-D array of them.

            impose_triu
                Flag to force the sparse matrix to only be the upper
                triangle of the covariance matrix.  Covariance matrix is
                symmetric such that C_ij = C_ji, so it's not necessary
                to keep both values.  This flag will force a call to
                sparse.triu when setting the covariance matrix.
                Otherwise, the input matrix is ASSUMED to only have the
                upper triangle of numbers.
        """

        if ifile is not None:
            self.read(ifile, impose_triu=impose_triu)
            return

        if inp is None:
            self._free()
            return
        else:
            self.cov = inp

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

        self._set_shape()

        self.input_indx = None

        if input_indx is not None:
            if self.dim == 2:
                raise Exception('Input indices only allowed when allocating multiple matrices.')
            if input_indx.shape != self.cov.shape:
                raise Exception('Input array and input index array must be the same size.')
            self.input_indx = input_indx
#           print(self.input_indx)

        if impose_triu:
            self._impose_upper_triangle()

        self.inv = None

    
    def __getitem__(self, *args):
        indx = tuple(*args)
        if len(indx) != self.dim:
            raise IndexError('Incorrect number of dimensions!')
        if self.dim == 2:
            return self.cov[tuple(sorted(indx))]
        return self.cov[self._grab_true_index(indx[:1])][tuple(sorted(indx[1:]))]


    def __del__(self):
        self._free()


    def _grab_true_index(self, inp):
        if self.input_indx is not None:
            i=0
            for coo in self.input_indx:
                if coo == inp:
                    return i
                i += 1
        return inp


    def _set_shape(self):
        self.shape = self.cov.shape
        if self.dim == 2:
            return
        self.shape += self.cov.ravel()[0].shape


    def _impose_upper_triangle(self):
        if self.dim == 2:
            self.cov = sparse.triu(self.cov).tocsr()
            self.nnz = self.cov.nnz
            return

        self.nnz = 0
        for i in range(0,self.shape[0]):
            self.cov[i] = sparse.triu(self.cov[i]).tocsr()
            self.nnz += self.cov[i].nnz

    def _with_lower_triangle(self, plane=None):
        if self.dim == 2:
            a = self.cov
        else:
            if plane is None:
                raise ValueError('Must define plane!  Use plane=...')
            a = self.cov[self._grab_true_index(plane)]

        return (sparse.triu(a) + sparse.triu(a,1).T)


    def _free(self):
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

        a = self.toarray(plane)

        if zoom is not None:
            xs = self.shape[self.dim-2]/2 - self.shape[self.dim-2]/2/zoom
            xe = xs + self.shape[self.dim-2]/zoom + 1
            ys = self.shape[self.dim-1]/2 - self.shape[self.dim-1]/2/zoom
            ye = ys + self.shape[self.dim-1]/zoom + 1
            a = a[xs:xe,ys:ye]

        if ofile is None:
            im = pyplot.imshow(a, interpolation='nearest', origin='lower')
            pyplot.colorbar()
            pyplot.show()
            return

        fig = pyplot.figure(1)
        im = pyplot.imshow(a, interpolation='nearest', origin='lower')
        pyplot.colorbar()
        fig.canvas.print_figure(ofile)
        pyplot.show()
        

    def toarray(self, plane=None):
        return (self._with_lower_triangle(plane)).toarray()
        

    def find(self, plane=None):
        return sparse.find(self._with_lower_triangle(plane))

#       if self.dim == 2:
#           return sparse.find(self.cov)
#       else:
#           if plane is None:
#               raise ValueError('Must define plane!  Use plane=...')
#           return sparse.find(self.cov[self._grab_true_index(plane)])

        

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

        # Use input header or create a minimal one
        # TODO: Ensure the input header can be used as a primary header!
        if hdr is None:
            hdr = fits.Header()
        hdr['SHAPE'] = str(self.shape)
        primary_hdu = fits.PrimaryHDU(header=hdr)

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

        # Index column
        if self.input_indx is not None:
            col1 = fits.Column(name='INP_INDX', format='1J', array=self.input_indx)

            indx_hdu = fits.BinTableHDU.from_columns([col1], name='PLANE')
            hdu.append(indx_hdu)

        hdu.writeto(ofile,clobber=clobber)


    def read(self, ifile, impose_triu=False, return_hdr=False, quiet=False):

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

        





        

