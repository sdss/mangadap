# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Light-weight class for handling data tables using numpy record
arrays.

----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
from IPython import embed

import numpy

from astropy.io import fits

from mangadap.util import fileio

# TODO: Put the general ParSet utilities used here (and by BitMask) in
# a more general 'rst' (?) utility module
from mangadap.par.parset import ParSet

# TODO: Can I just subclass from recarray?
# TODO: Add units

class DataTable:
    """
    A light-weight class that wraps a numpy record array.

    Args:
        keys (:obj:`str`, :obj:`list`):
            Table column keys. Must be a single string or a list of
            strings.
        types (:obj:`type`, :obj:`list`):
            Column data type. Must be a single type or a list of
            types. Supported types are limited to integers and
            floating point values. If a single type and more than one
            key is provided, all columns will have the same type.
        element_shapes (:obj:`tuple`, :obj:`list`, optional):
            Shape of each column element. If None, all elements are
            assumed to be single values. If a list is provided, the
            elements of the list must be either None (for a single
            value) or a tuple with the shape of each element. If a
            single shape and more than one key is provided, all
            columns will have the same element shape.
        descr (:obj:`str`, :obj:`list`, optional):
            Human-readable description of each column. If None, no
            descriptions are available. If a list, each element must
            be a string or None. If not None, must provide a
            description or a None element for each key.
        shape (:obj:`int`, :obj:`tuple`):
            Shape of the data table to instantiate. If None,
            :attr:`data` will be None and instantiation of the array
            will require a call to :func:`init`.

    Raises:
        TypeError:
            Raised if any of the provided arguments or list item have
            an incorrect type.
        ValueError:
            Raised if any of the provided arguments do not provide
            the same number of columns as determined by the column
            keys.

    Attributes:
        keys (:obj:`list`):
            See argument list.
        types (:obj:`list`):
            See argument list.
        element_shapes (:obj:`list`):
            See argument list.
        descr (:obj:`list`):
            See argument list.
        data (`numpy.recarray`_):
            Array with the table data.
        _dtype (:obj:`list`):
            List of tuples with the data type of the record array.
    """
    def __init__(self, keys, types, element_shapes=None, descr=None, shape=None):

        self.keys = [keys] if isinstance(keys, str) else list(keys)
        if not isinstance(self.keys, list):
            raise TypeError('Input keys must be a single string or a list.')
        if not all([isinstance(key, str) for key in self.keys]):
            raise TypeError('Input keys must be strings.')

        self.ncols = len(self.keys)

        self.types = [types] if isinstance(keys, type) else list(types)
        if not isinstance(self.types, list):
            raise TypeError('Input types must be a single type or a list.')
        if not all([isinstance(numpy.dtype(t).type(0),
                               (int, numpy.integer, float, numpy.floating, str, numpy.str_,
                                bool, numpy.bool_)) 
                        for t in self.types]):
            raise TypeError('Input types must cast to an integer or floating point object.')
        if len(self.types) == 1:
            self.types = [self.types[0]]*self.ncols
        if len(self.types) != self.ncols:
            raise ValueError('Number of types does not match the number of columns.')

        self.element_shapes = [None]*self.ncols if element_shapes is None \
                                else ([element_shapes] if isinstance(element_shapes, tuple) 
                                      else element_shapes)
        if not isinstance(self.element_shapes, list):
            raise TypeError('Input element_shapes must be None, a single tuple, or a list.')
        if not all([e is None or isinstance(e, tuple) for e in self.element_shapes]):
            raise TypeError('Input element_shapes must be None or a tuple.')
        if len(self.element_shapes) == 1:
            self.element_shapes = [self.element_shapes[0]]*self.ncols
        if len(self.element_shapes) != self.ncols:
            raise ValueError('Number of element shapes does not match the number of columns.')
        
        self.descr = [None]*self.ncols if descr is None \
                                else ([descr] if isinstance(descr, str) else descr)
        if not isinstance(self.descr, list):
            raise TypeError('Input descriptions must be None, a string, or a list.')
        if not all([d is None or isinstance(d, str) for d in self.descr]):
            raise TypeError('Input descriptions must be None or a string.')
        if len(self.descr) != self.ncols:
            raise ValueError('Number of descriptions does not match the number of columns.')

        # Build the record array dtype
        self._dtype = [(k, t) if s is None else (k,t,s) 
                            for k,t,s in zip(self.keys, self.types, self.element_shapes)]

        # Instantiate the array if the shape is defined
        if shape is None:
            self.data = None
        else:
            self.init(shape)

    def __getitem__(self, key):
        if self.data is None:
            raise ValueError('{0} is empty.  First instantiate using init().'.format(
                             self.__class__.__name__))        
        return self.data[key]

    def __setitem__(self, key, value):
        if self.data is None:
            raise ValueError('{0} is empty.  First instantiate using init().'.format(
                             self.__class__.__name__))        
        self.data[key] = value

    def __len__(self):
        return self.size

    def init(self, shape):
        """
        (Re)Initialize :attr:`data`.

        Note that this *always* re-initializes the array, regardless
        of whether or not there was a pre-existing data array.

        Args:
            shape (:obj:`tuple`):
                The shape of the array to construct.
        """
        self.data = fileio.init_record_array(shape, self._dtype)

    def append(self, rhs):
        """
        Append the data in the provided datatable to this one.

        Both objects (``self`` and ``rhs``) must have the same
        derived class.  This instance is directly modified.

        Args:
            rhs (:class:`DataTable`):

                The instance with the data to append. The derived
                class of this object **must** be the same as the
                object to which you want to append the data.  I.e.::

                    assert isinstance(rhs, self.__class__)

        """
        if not isinstance(rhs, self.__class__):
            #if rhs.__class__.__name__ != self.__class__.__name__:
            raise TypeError('Object to append must be of type {0}!'.format(self.__class__.__name__))
        self.data = rhs.data.copy() if self.data is None \
                        else numpy.append(self.data, rhs.data)

    @property
    def shape(self):
        """
        Return the shape of the data array.

        Returns:
            :obj:`tuple`: Returns the shape of the data array. If the
            :attr:`data` has not been instantiated, returns None.
        """
        return None if self.data is None else self.data.shape

    @property
    def size(self):
        """
        Return the size of the data array.

        Returns:
            :obj:`int`: Returns the number of elements in the array. If the
            :attr:`data` has not been instantiated, returns 0.
        """
        return 0 if self.data is None else self.data.size

    @property
    def dtype(self):
        """
        Return the data type of the record array.

        If :attr:`data` has been initialized, this is identically::

            self.data.dtype

        Returns:
            `numpy.dtype`_: Table data type.
        """
        return numpy.dtype(self._dtype) if self.data is None else self.data.dtype

    def to_rst_table(self, header=True, class_link=True):
        """
        Construct a reStructuredText table describing the data table.

        Args:
            header (:obj:`bool`, optional):
                Include a section header
            class_link (:obj:`bool`, optional):
                Include an rst-style link to the class instantiation
                documentation.

        Returns:
            :obj:`list`: Returns a list of lines that can be written
            to an ``*.rst`` file.
        """
        data_table = numpy.empty((self.ncols+1, 3), dtype=object)
        data_table[0,:] = ['Key', 'Type', 'Description']
        for i in range(self.ncols):
            data_table[i+1,0] = self.keys[i]
            data_table[i+1,1] = '``{0}``'.format(numpy.dtype(self.types[i]).type.__name__)
            data_table[i+1,2] = '..' if self.descr[i] is None else self.descr[i]

        output = []
        if header:
            output += [ '{0} Data Table'.format(type(self).__name__) ]
            output += [ '-'*len(output[0]) ]
            output += [ '' ]
        if class_link:
            output += ['Class Instantiation: ' + ParSet._rst_class_name(self)]
            output += ['']
        output += [ParSet._data_table_string(data_table, delimiter='rst')]
        output += ['']
        return output

    def to_hdu(self, name=None, add_dim=False):
        """
        Use the data table to build an `astropy.io.fits.BinTableHDU`_.

        .. todo::

            - ``add_dim`` shouldn't be optional; the code should just
              do this always or know when it needs to.

        Args:
            name (:obj:`str`, optional):
                HDU extension name.
            add_dim(:obj:`bool`, optional):

                Include the dimensionality of the column in the
                construction of the `astropy.io.fits.Column`_.

        Returns:
            `astropy.io.fits.BinTableHDU`_: FITS binary table with
            the data.
        """
        return fits.BinTableHDU.from_columns(
                    [fits.Column(name=n, format=fileio.rec_to_fits_type(self.data[n]),
                                 dim=fileio.rec_to_fits_col_dim(self.data[n]) if add_dim else None,
                                 array=self.data[n]) for n in self.data.dtype.names],
                    name=name)

    # TODO: Add a from_hdu method.
