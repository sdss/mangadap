from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import sys
if sys.version > '3':
    long = int

import os.path
from os import environ
from scipy import sparse
from astropy.io import fits
import time
import numpy

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

    def __init__(self, inp=None, impose_triu=False):
        """
        ARGUMENTS:
            inp 
                Input matrix.  Can be either a single
                scipy.sparse.csr.csr_matrix or an array of them.

            impose_triu
                Flag to force the sparse matrix to only be the upper
                triangle of the covariance matrix.  Covariance matrix is
                symmetric such that C_ij = C_ji, so it's not necessary
                to keep both values.  This flag will force a call to
                sparse.triu when setting the covariance matrix.
                Otherwise, the input matrix is ASSUMED to only have the
                upper triangle of numbers.
        """

        if inp is None:
            self._free()
        else:
            self.cov = inp

        if type(self.cov) == numpy.ndarray:
            self.dim = len(self.cov.shape)+2
            for cov in self.cov.ravel():
                if not sparse.isspmatrix_csr(cov):
                    print(type(cov))
                    self._free()
                    raise TypeError('Incorrect object type!')
        else:
            self.dim = 2
            if not sparse.isspmatrix_csr(self.cov):
                print(type(self.cov))
                self._free()
                raise TypeError('Incorrect object type!')

        self.shape = []
        if self.dim > 2:
            for nn in self.cov.shape:
                self.shape.append(nn)
            nx, ny = self.cov.ravel()[0].shape
        else:
            nx, ny = self.cov.shape
        self.shape.append(nx)
        self.shape.append(ny)
        self.shape = tuple(self.shape)

        if impose_triu:
            if self.dim != 2:
                for cov in self.cov.ravel():
                    cov = sparse.triu(cov).tocsr()
            else:
                self.cov = sparse.triu(self.cov).tocsr()

    
    def __getitem__(self, *args):
        indx = tuple(*args)
        if len(indx) != self.dim:
            raise IndexError('Incorrect number of dimensions!')

        if self.dim == 2:
            return self.cov[tuple(sorted(indx))]
        return self.cov[indx[:self.dim-2]][tuple(sorted(indx[self.dim-2:]))]


    def __del__(self):
        self._free()


    def _free(self):
        self.dim=0
        self.shape=None

        try:
            del self.cov
            self.cov = None
        except NameError as e:
            self.cov = None

        
#   def write(self, ofile):

#       fits.open




        

