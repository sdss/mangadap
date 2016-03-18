# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Define a utility base class used to hold parameters.

*License*:
    Copyright (c) 2015, Kyle B. Westfall
    Licensed under BSD 3-clause license - see LICENSE.rst

*Source location*:
    $MANGADAP_DIR/python/mangadap/par/parset.py

*Python2/3 compliance*::

    from __future__ import division
    from __future__ import print_function
    from __future__ import absolute_import
    from __future__ import unicode_literals
    
    import sys
    if sys.version > '3':
        long = int

*Imports*::

    import numpy

..todo::

    - Add range and length parameters allowing one to define the range
      allowed for the parameter values and number of elements required
      (if the parameter is an array)

*Revision history*:
    | **16 Jun 2015**: Original implementation by K. Westfall (KBW)
    | **18 Mar 2016**: (KBW) Change dtype checking

"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import sys
if sys.version > '3':
    long = int

import numpy

__author__ = 'Kyle B. Westfall'

class ParSet:
    """
    Generic base class to handle and manipulate a list of operational
    parameters.
    """
    def __init__(self, pars, values=None, defaults=None, options=None, dtypes=None):
        # Check that the list of input parameters is a list of strings
        if type(pars) != list:
            raise Exception('Input parameter keys must be provided as a list.')
        for key in pars:
            if type(key) != str:
                raise Exception('Input parameter keys must be strings.')
        
        # Get the length of the parameter list and make sure the list
        # has unique values
        self.npar = len(pars)
        if len(numpy.unique(numpy.array(pars))) != self.npar:
            raise Exception('All input parameter keys must be unique.')

        # Check that the other lists, if provided, have the correct type
        # and length
        if values is not None and (type(values) != list or len(values) != self.npar):
            raise Exception('Values must be a list with the same length as the number of keys.')
        if defaults is not None and (type(defaults) != list or len(defaults) != self.npar):
            raise Exception('Defaults must be a list with the same length as the number of keys.')
        if options is not None and (type(options) != list or len(options) != self.npar):
            raise Exception('Options must be a list with the same length as the number of keys.')
        if dtypes is not None and (type(dtypes) != list or len(dtypes) != self.npar):
            raise Exception('Data types must be a list with the same length as the number of keys.')

        # Set up dummy lists for no input
        if values is None:
            values = [None]*self.npar
        if defaults is None:
            defaults = [None]*self.npar
        if options is None:
            options = [None]*self.npar
        if dtypes is None:
            dtypes = [None]*self.npar

        # Set the valid options
        self.options = dict([ (p, [o]) if o is not None and type(o) != list else (p, o) \
                                       for p, o in zip(pars, options) ])
        # Set the valid types
        self.dtype = dict([ (p, [t]) if t is not None and type(t) != list else (p, t) \
                                     for p, t in zip(pars, dtypes) ])

        # Set the data dictionary using the internal functions
        self.data = {}
        for p, d, v in zip(pars, defaults, values):
            if v is None:
                self.__setitem__(p, d)
                continue
            self.__setitem__(p, v)


    def __getitem__(self, key):
        return self.data[key]


    def __setitem__(self, key, value):
        if value is None:
            self.data[key] = value
            return

        if self.options[key] is not None and value not in self.options[key]:
            raise ValueError('Input value invalid: {0}.\nOptions are: {1}'.format(value,
                                                                                self.options[key]))
        if self.dtype[key] is not None \
                and not any([ isinstance(value, d) for d in self.dtype[key]]):
            raise ValueError('Input value incorrect type: {0}.\nOptions are: {1}'.format(value,
                                                                                self.dtype[key]))

        self.data[key] = value
        

    def __iter__(self):
        return iter(self.data.values())

    
    def add(self, key, value, options=None, dtype=None):
        self.options[key] = [options] if options is not None and type(options) != list else options
        self.dtype[key] = [dtype] if dtype is not None and type(dtype) != list else dtype
        self.__setitem__(key, value)


