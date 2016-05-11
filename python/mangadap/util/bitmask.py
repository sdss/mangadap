# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Base class for handling bit masks by the DAP.

*License*:
    Copyright (c) 2015, SDSS-IV/MaNGA Pipeline Group
        Licensed under BSD 3-clause license - see LICENSE.rst

*Source location*:
    $MANGADAP_DIR/python/mangadap/util/bitmask.py

*Imports and python version compliance*:
    ::

        from __future__ import division
        from __future__ import print_function
        from __future__ import absolute_import
        from __future__ import unicode_literals

        import sys
        if sys.version > '3':
            long = int
            try:
                from configparser import ConfigParser
            except ImportError:
                print('WARNING: Unable to import configparser!  Beware!')
        else:
            try:
                from ConfigParser import ConfigParser
            except ImportError:
                print('WARNING: Unable to import ConfigParser!  Beware!')

        import numpy
        import os
        import textwrap
        from pydl.pydlutils.yanny import yanny

*Class usage examples*:
    Here is an example of reading/creating a TemplateLibrary, then
    creating a plot that shows the full spectrum (blue) and the unmasked
    spectrum (green)::

        # Imports
        import numpy
        from mangadap.drpfits import DRPFits
        from mangadap.proc.TemplateLibrary import TemplateLibrary, TemplateLibraryBitMask
        from matplotlib import pyplot

        # Define the DRP file
        drpf = DRPFits(7495, 12703, 'CUBE')

        # Build the template library
        tpl_lib = TemplateLibrary('M11-MILES', drpf=drpf, directory_path='.')
        # Writes: ./manga-7495-12703-LOGCUBE_M11-MILES.fits

        # Initialize the mask object
        tplbm = TemplateLibraryBitMask()

        # Refactor the mask for the first template spectrum using the bitmask
        unmasked = numpy.invert(tplbm.flagged(tpl_lib.hdu['MASK'].data[0,:]))

        # Plot the full spectrum (blue)
        pyplot.plot(tpl_lib.hdu['WAVE'].data, tpl_lib.hdu['FLUX'].data[0,:])
        # Plot the unmasked pixels (green; points are connected even if there are gaps
        pyplot.plot(tpl_lib.hdu['WAVE'].data[unmasked], tpl_lib.hdu['FLUX'].data[0,unmasked])
        # Show the plot
        pyplot.show()


*Revision history*:
    | **01 Jun 2015**: Original implementation by K. Westfall (KBW)
    | **07 Oct 2015**: (KBW) Added a usage case
    | **29 Jan 2016**: (KBW) Changed attributes of :class:`BitMask` and
        added functionality to print a description of the bits.  Convert
        :class:`mangadap.proc.templatelibrary.TemplateLibraryBitMask` to
        new format where the bits are read from a configuration file.
    | **17 Feb 2016**: (KBW) Minor edit to documentation
    | **16 Mar 2016**: (KBW) Moved TemplateLibraryBitMask to
        :class:`mangadap.proc.templatelibrary.TemplateLibraryBitMask`;
        moved HDUList_mask_wavelengths to
        :func:`mangadap.proc.util.HDUList_mask_wavelengths`.
    | **27 Mar 2016**: (KBW) Added :func:`BitMask.from_ini_file` and
        :func:`BitMask.from_par_file` class methods.  Added
        :func:`BitMask._fill_sequence` static method.  This allows for
        :class:`BitMask` objects to be declared directly from the files,
        and allows the bit values to take on any number.
    | **05 Apr 2016**: (KBW) Added parameters to initialization of
        :class:`BitMask` objects to clean up some of the derived class
        initialization.
    | **11 May 2016**: (KBW) Switch to using `pydl.pydlutils.yanny`_
        instead of internal yanny reader

.. _pydl.pydlutils.yanny: http://pydl.readthedocs.io/en/stable/api/pydl.pydlutils.yanny.yanny.html
.. _SDSS-style parameter file: http://www.sdss.org/dr12/software/par/

"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import sys
if sys.version > '3':
    long = int
    try:
        from configparser import ConfigParser
    except ImportError:
        print('WARNING: Unable to import configparser!  Beware!')
else:
    try:
        from ConfigParser import ConfigParser
    except ImportError:
        print('WARNING: Unable to import ConfigParser!  Beware!')

import numpy
import os
import textwrap
from pydl.pydlutils.yanny import yanny
#from .yanny import yanny

__author__ = 'Kyle B. Westfall'

class BitMask:
    """
    Generic class to handle and manipulate bitmasks.  The input list of
    bit names (keys) must be unique, except that values of 'NULL' are
    ignored.  The index in the input keys determines the bit value;
    'NULL' keys are included in the count.  For example::

        >>> from mangadap.util.bitmask import BitMask
        >>> keys = [ 'key1', 'key2', 'NULL', 'NULL', 'key3' ]
        >>> bm = BitMask(keys)
        >>> bm.info()
                 Bit: key1 = 0

                 Bit: key2 = 1

                 Bit: key3 = 4

    Args:

        keys (str, list, numpy.ndarray): (**Optional**) List of keys (or single key) to
            use as the bit name.  Each key is given a bit number ranging
            from 0..N-1.  If **not** provided, one of either *ini_file*
            or *par_file* must be provided.

        descr (str, list, numpy.ndarray): (Optional) List of descriptions
            (or single discription) provided by :func:`info` for each bit.

    Raises:
        ValueError: Raised if more than 64 bits are provided.
        TypeError: Raised if the provided `keys` do not have the correct
            type.

    Attributes:
        nbits (int): Number of bits
        bits (dict): A dictionary with the bit name and value
        descr (numpy.ndarray) : List of bit descriptions
        max_value (int): The maximum valid bitmask value given the
            number of bits.

    """
    def __init__(self, keys=None, descr=None, ini_file=None, par_file=None, par_grp=None):

        if keys is None and ini_file is None and par_file is None:
            raise ValueError('BitMask undefined!')

        if par_file is not None and par_grp is None:
            raise ValueError('If reading from par file, must provide bit group flag.')

        # Allow for file-based initialization
        _self = None
        if ini_file is not None:
            _self = BitMask.from_ini_file(ini_file)
        if par_file is not None :
            _self = BitMask.from_par_file(par_file, par_grp)
        if _self is not None:
            self.nbits = _self.nbits
            self.bits = _self.bits
            self.descr = _self.descr
            self.max_value = _self.max_value
            return

        if not isinstance(keys, (str, list, numpy.ndarray)):
            raise TypeError('Input argument \'keys\' incorrect type.')
        if not all([isinstance(k, str) for k in keys]):
            raise TypeError('Input keys must have string type.')

        if isinstance(keys, list):
            keys = numpy.array(keys)
        if isinstance(keys, str):
            keys = numpy.array([keys])

        # Do not allow for more that 64 bits
        if len(keys) > 64:
            raise ValueError('Can only define up to 64 bits!')

        # Allow for multiple NULL keys; but check the rest for
        # uniqueness
        diff = set(keys) - set(['NULL'])
        if len(diff) != numpy.unique(keys[ keys != 'NULL' ]).size:
            raise ValueError('All input keys must be unique.')

        self.nbits = len(keys)
        self.bits = dict([ (k,i) for i,k in enumerate(keys) ])
        self.max_value = (1 << self.nbits)-1

        self.descr = None
        if descr is None:
            return

        if not isinstance(descr, (str, list, numpy.ndarray)):
            raise TypeError('Input argument \'descr\' incorrect type.')
        if not all([isinstance(d, str) for d in descr]):
            raise TypeError('Input keys must have string type.')

        if isinstance(descr, list):
            descr = numpy.array(descr)
        if isinstance(descr, str):
            descr = numpy.array([descr])

        if len(descr) != self.nbits:
            raise ValueError('Number of listed descriptions not the same as number of keys.')

        self.descr = descr.copy()


    @classmethod
    def from_ini_file(cls, f):
        r"""
        Define the object using a ini configuration file.  The sections
        of the configuration file define the keys, and each section is
        expected to have `value` and `descr` components that define the
        bit value and provided a description of the bit.  An example ini
        file might look like this::

            [NO_DATA]
             value = 0
             descr = Pixel has no data

            [INVALID_DATA]
             value = 1
             descr = Pixel value is invalid
    
        See
        :class:`mangadap.proc.templatelibrary.TemplateLibraryBitMask`
        for an example that uses this function.

        Args:
            f (str) : File name to use for defining the
                :class:`BitMask`.
        
        Returns:
            :class:`BitMask`: Object with bitmasks defined by the ini
            file.

        Raises:
            FileNotFoundError: Raised if the input file does not exist.
        """
        # Check if the file exists
        if not os.path.isfile(f):
            raise FileNotFoundError('Could not find ini file: {0}'.format(f))
        # Define the parser and read the file
        cnfg = ConfigParser()
        cnfg.read(f)

        # Read the keys, values, and descriptions
        keys = numpy.array(cnfg.sections())
        vals = numpy.zeros( keys.size, dtype=numpy.int )
        descr = numpy.zeros( keys.size, dtype=object )
        for i,k in enumerate(keys):
            vals[i] = cnfg[k]['value']
            descr[i] = cnfg[k]['descr']

        # Slot in NULLs where necessary and return the object instance
        keys, vals, descr = cls._fill_sequence(keys, vals, descr)
        srt = numpy.argsort(vals)
        return cls(keys[srt], descr=descr[srt])


    @classmethod
    def from_par_file(cls, f, name):
        r"""
        Define the object using an `SDSS-style parameter file`_.  This
        has been tailored to work with the sdssMaskbits.par file in
        IDLUTILS; however, it can work with other like files.
        
        See :class:`mangadap.drpfits.DRPFitsBitMask` for an example that
        uses this function.

        Args:
            f (str) : File name to use for defining the
                :class:`BitMask`.
            name (str) : The designation of the bits to assign.  For
                example, in :class:`mangadap.drpfits.DRPFitsBitMask`
                this is `'MANGA_DRP3PIXMASK'`.
        
        Returns:
            :class:`BitMask`: Object with bitmasks defined by the
            parameter file.

        Raises:
            FileNotFoundError: Raised if the input file does not exist.
        """
        # Check the file exists
        if not os.path.isfile(f):
            raise FileNotFoundError('Could not find ini file: {0}'.format(f))

        # Read the full yanny file and only select the maskbits typedef
#        bits = yanny(f)['MASKBITS']
        bits = yanny(filename=f, raw=True)['MASKBITS']

        # Find the bits with the correct designation
        indx = numpy.array(bits['flag']) == name
        keys = numpy.array(bits['label'])[indx]
        vals = numpy.array(bits['bit'])[indx]
        descr = numpy.array(bits['description'])[indx]

        # Slot in NULLs where necessary and return the object instance
        keys, vals, descr = cls._fill_sequence(keys, vals, descr)
        srt = numpy.argsort(vals)
        return cls(keys[srt], descr=descr[srt])


    @staticmethod
    def _fill_sequence(keys, vals, descr):
        r"""
        The instantiation of :class:`BitMask` does not include the value
        of the bit, it just assumes that the bits should be in sequence
        such that the first key has a value of 0, and the last key has a
        value of N-1.  Unfortunately, not all the bits the DAP needs to
        use follow a sequence.  This function finds the range of the
        bits and then slots in NULL keywords and empty descriptions
        where necessary to fill in the full complement of bits.  NULL
        keywords are ignored by the :class:`BitMask` object.

        This is a static method because it doesn't depend on any of the
        attributes of :class:`BitMask`.

        Args:
            keys (list or str): Bit names
            vals (list or int): Bit values
            descr (list or str): Description of each bit

        Returns:
            numpy.ndarray: Three arrays with the filled keys, values,
            and descriptions.

        Raises:
            ValueError: Raised if a bit value is less than 0.
        """
        # Make sure the input is treated as an array
        keys = numpy.asarray(keys)
        vals = numpy.asarray(vals)
        descr = numpy.asarray(descr)

        if numpy.amin(vals) < 0:
            raise ValueError('No bit cannot be less than 0!')
        minv = numpy.amin(vals)
        maxv = numpy.amax(vals)

        if minv != 0 or maxv != len(vals)-1:
            diff = list(set(numpy.arange(maxv)) - set(vals))
            vals = numpy.append(vals, diff)
            keys = numpy.append(keys, numpy.array(['NULL']*len(diff)))
            descr = numpy.append(descr, numpy.array(['']*len(diff)))

        return keys, vals, descr

#        if numpy.amin(vals) != 0:
#            raise ValueError('Bit values must start from 0.')
#        if numpy.amax(vals) != keys.size-1:
#            raise ValueError('Bit values must be sequential starting from 0.')

        
    def keys(self):
        """
        Return a list of the bits; 'NULL' keywords are ignored.
        
        Returns:
            list: List of bit keywords.
        """
        return list(set(self.bits.keys())-set(['NULL']))


#    def nbits(self):
#        """
#        Return the number of bits.
#
#        Returns:
#            int : Number of bits
#        """
#        return len(self.bits)


    def info(self):
        """
        Print the list of bits and, if available, their descriptions.
        """
        try:
            tr, tcols = numpy.array(os.popen('stty size', 'r').read().split()).astype(int)
            tcols -= int(tcols*0.1)
        except:
            tr = None
            tcols = None

        for k,v in sorted(self.bits.items(), key=lambda x:(x[1],x[0])):
            if k == 'NULL':
                continue
            print('         Bit: {0} = {1}'.format(k,v))
            if self.descr is not None:
                if tcols is not None:
                    print(textwrap.fill(' Description: {0}'.format(self.descr[v]), tcols))
                else:
                    print(' Description: {0}'.format(self.descr[v]))
            print(' ')


    def minimum_uint_dtype(self):
        """
        Return the smallest uint datatype that is needed to contain all
        the bits in the mask.
        """
        if self.nbits < 8:
            return numpy.uint8
        if self.nbits < 16:
            return numpy.uint16
        if self.nbits < 32:
            return numpy.uint32
        return numpy.uint64


    def flagged(self, value, flag=None):
        """
        Determine if a bit is on in the provided bitmask value.  The
        function can be used to determine if any individual bit is on or
        any one of many bits is on.

        Args:
            value (uint or array): Bitmask value.  It should be less
                than or equal to :attr:`max_value`; however, that is not
                checked.
            flag (string or array): (Optional) Bit names to check.  If
                None, then it checks if any bit is on.
        
        Returns:
            bool: Boolean flags that the provided flags (or any flag)
            is on for the provided bitmask value.

        Raises:
            KeyError: Raised by the dict data type if the input *flag*
                is not one of the valid :attr:`flags`.
            TypeError: Raised if the provided *flag* does not contain
                one or more strings.
        """
        if numpy.any(numpy.array(flag) == 'NULL'):
            raise ValueError('Flag name NULL is not allowed.')

        if flag is None:
            flag = self.keys()

        if isinstance(flag, str):
            return value & (1 << self.bits[flag]) != 0

        if hasattr(flag, "__iter__"):
            if not all([ isinstance(f, str) for f in flag ]):
                raise TypeError('Provided bit names must be strings!')
            out = value & (1 << self.bits[flag[0]]) != 0
            nn = len(flag)
            for i in range(1,nn):
                out |= (value & (1 << self.bits[flag[i]]) != 0)
            return out

        raise Exception('Provided bit name must be a string or an iterable of strings!')


    def flagged_bits(self, value):
        """
        Return the list of flagged bit names for a single bit value.

        Args:
            value (uint): Bitmask value.  It should be less than or
                equal to :attr:`max_value`; however, that is not
                checked.
        
        Returns:
            list: List of flagged bit value keywords.

        Raises:
            KeyError: Raised by the dict data type if the input *flag*
                is not one of the valid :attr:`flags`.
            TypeError: Raised if the provided *flag* does not contain
                one or more strings.
        """
        if not isinstance(value, int):
            raise TypeError('Input must be a single integer.')
        if value <= 0:
            return []
        keys = numpy.array(self.keys())
        indx = numpy.array([ 1<<self.bits[k] & value != 0 for k in keys])
        return list(keys[indx])


    def toggle(self, value, flag):
        """
        Toggle a bit in the provided bitmask value.

        Args:
            value (uint or array): Bitmask value.  It should be less
                than or equal to :attr:`max_value`; however, that is not
                checked.
            flag (list, numpy.ndarray, or str): Bit name(s) to toggle.
        
        Returns:
            uint: New bitmask value after toggling the selected bit.

        Raises:
            KeyError: Raised by the dict data type if the input *flag*
                is not one of the valid :attr:`flags`.
            Exception: Raised if the provided *flag* is not a string.
        """
        if flag == 'NULL':
            raise ValueError('Flag name NULL is invalid!')
        if isinstance(flag, list) or isinstance(flag, numpy.ndarray):
            _value = value
            for f in flag:
                if not isinstance(f, str):
                    raise TypeError('Provided bit name must be a string!')
                _value ^= (1 << self.bits[f])
            return _value
        if not isinstance(flag, str):
            raise TypeError('Provided bit name must be a list or string!')
        return value ^ (1 << self.bits[flag])


    def turn_on(self, value, flag):
        """
        Ensure that a bit is turned on in the provided bitmask value.

        Args:
            value (uint or array): Bitmask value.  It should be less
                than or equal to :attr:`max_value`; however, that is not
                checked.
            flag (list, numpy.ndarray, or str): Bit name(s) to turn on.
        
        Returns:
            uint: New bitmask value after turning on the selected bit.

        Raises:
            KeyError: Raised by the dict data type if the input *flag*
                is not one of the valid :attr:`flags`.
            Exception: Raised if the provided *flag* is not a string.
        """
        if flag == 'NULL':
            raise ValueError('Flag name NULL is invalid!')
        if isinstance(flag, list) or isinstance(flag, numpy.ndarray):
            _value = value
            for f in flag:
                if not isinstance(f, str):
                    raise TypeError('Provided bit name must be a string!')
                _value |= (1 << self.bits[f])
            return _value
        if not isinstance(flag, str):
            raise TypeError('Provided bit name must be a string!')
        return value | (1 << self.bits[flag])

#        if flag == 'NULL':
#            raise ValueError('Flag name NULL is invalid!')
#        if not isinstance(flag, str):
#            raise Exception('Provided bit name must be a string!')
#        return value | (1 << self.bits[flag])


    def turn_off(self, value, flag):
        """
        Ensure that a bit is turned off in the provided bitmask value.

        Args:
            value (uint or array): Bitmask value.  It should be less
                than or equal to :attr:`max_value`; however, that is not
                checked.
            flag (list, numpy.ndarray, or str): Bit name(s) to turn off.
        
        Returns:
            uint: New bitmask value after turning off the selected bit.

        Raises:
            KeyError: Raised by the dict data type if the input *flag*
                is not one of the valid :attr:`flags`.
            Exception: Raised if the provided *flag* is not a string.
        """
        if flag == 'NULL':
            raise ValueError('Flag name NULL is invalid!')
        if isinstance(flag, list) or isinstance(flag, numpy.ndarray):
            _value = value
            for f in flag:
                if not isinstance(f, str):
                    raise TypeError('Provided bit name must be a string!')
                _value &= ~(1 << self.bits[f])
            return _value
        if not isinstance(flag, str):
            raise TypeError('Provided bit name must be a string!')
        return value & ~(1 << self.bits[flag])

#        if flag == 'NULL':
#            raise ValueError('Flag name NULL is invalid!')
#        if isinstance(flag, str):
#            raise Exception('Provided bit name must be a string!')
#        return value & ~(1 << self.bits[flag])



