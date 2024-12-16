# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Base class for handling bit masks.

Class usage examples
--------------------

.. include:: ../include/bitmask_usage.rst

----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import os
import textwrap
from configparser import ConfigParser

import numpy

from pydl.pydlutils.yanny import yanny

from mangadap.par.parset import ParSet

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

    .. todo::
        - Have the class keep the mask values internally instead of
          having it only operate on the mask array...

    Args:
        keys (:obj:`str`, :obj:`list`):
            List of keys (or single key) to use as the bit name.  Each
            key is given a bit number ranging from 0..N-1.
        descr (:obj:`str`, :obj:`list`, optional):
            List of descriptions (or single description) provided by
            :func:`info` for each bit.  No descriptions by default.

    Raises:
        ValueError:
            Raised if more than 64 bits are provided.
        TypeError:
            Raised if the provided `keys` do not have the correct type.

    Attributes:
        nbits (int):
            Number of bits
        bits (dict):
            A dictionary with the bit name and value
        descr (numpy.ndarray):
            List of bit descriptions
        max_value (int):
            The maximum valid bitmask value given the number of bits.
    """
    prefix = 'BIT'
    def __init__(self, keys, descr=None):

        _keys = keys if hasattr(keys, '__iter__') else [keys]
        _keys = numpy.atleast_1d(_keys).ravel()
        _descr = None if descr is None else numpy.atleast_1d(descr).ravel()
        if _descr is not None:
            for i in range(len(_descr)):
                _descr[i] = _descr[i].strip()

        if _descr is not None:
            if not all([isinstance(d, str) for d in _descr]):
                raise TypeError('Input descriptions must have string type.')
            if len(_descr) != len(_keys):
                raise ValueError('Number of listed descriptions not the same as number of keys.')

        # Do not allow for more that 64 bits
        if len(_keys) > 64:
            raise ValueError('Can only define up to 64 bits!')

        # Allow for multiple NULL keys; but check the rest for
        # uniqueness
        diff = set(_keys) - set(['NULL'])
        if len(diff) != numpy.unique(_keys[[k != 'NULL' for k in _keys]]).size:
            raise ValueError('All input keys must be unique.')

        # Initialize the attributes
        self.nbits = len(_keys)
        self.bits = { k:i for i,k in enumerate(_keys) }
        self.max_value = (1 << self.nbits)-1
        self.descr = _descr

    # TODO: Add a to_ini_file method
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
        vals = numpy.zeros(keys.size, dtype=int)
        descr = numpy.zeros(keys.size, dtype=object)
        for i,k in enumerate(keys):
            vals[i] = cnfg[k]['value']
            descr[i] = cnfg[k]['descr']

        # Slot in NULLs where necessary and return the object instance
        keys, vals, descr = cls._fill_sequence(keys, vals, descr)
        srt = numpy.argsort(vals)
        return cls(keys[srt], descr=descr[srt])

    # TODO: Add a to_par_file method
    @classmethod
    def from_par_file(cls, f, name):
        r"""
        Define the object using an `SDSS-style parameter file`_.  This
        has been tailored to work with the sdssMaskbits.par file in
        IDLUTILS; however, it can work with similar files.

        See :class:`mangadap.util.drpbitmask.DRPFitsBitMask` for an
        example that uses this function.

        Args:
            f (:obj:`str`):
                File name to use for defining the :class:`BitMask`.
            name (:obj:`str`):
                The designation of the bits to assign. For example,
                in :class:`mangadap.util.drpbitmask.DRPFitsBitMask` this
                is `'MANGA_DRP3PIXMASK'`.
        
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
        bits = yanny(filename=f, raw=True)['MASKBITS']

        # Find the bits with the correct designation
        indx = numpy.array(bits['flag']) == name
        keys = numpy.array(bits['label'])[indx]
        vals = numpy.array(bits['bit'])[indx]
        descr = numpy.array(bits['description'])[indx]

        # Slot in NULLs where necessary and return the object instance
        keys, vals, descr = cls._fill_sequence(keys, vals, descr)
        srt = numpy.argsort(vals)
        return cls(keys[srt], descr=None if descr is None else descr[srt])

    def _prep_flags(self, flag):
        """Prep the flags for use."""
        # Flags must be a numpy array
        _flag = numpy.array(self.keys()) if flag is None else numpy.atleast_1d(flag).ravel()
        # NULL flags not allowed
        if numpy.any([f == 'NULL' for f in _flag]):
            raise ValueError('Flag name NULL is not allowed.')
        # Flags should be among the bitmask keys
        if numpy.any([f not in self.keys() for f in _flag]):
            raise ValueError('Some bit names not recognized.')
        return _flag

    def _init_objs(self):
        """
        Return the objects needed to instantate another BitMask object
        that's identical to self.
        """
        keys = self.keys()
        vals = [self.bits[k] for k in keys]
        descr = None if self.descr is None else [self.descr[self.bits[k]] for k in keys]
        keys, vals, descr = BitMask._fill_sequence(keys, vals, descr)
        srt = numpy.argsort(vals)
        return keys[srt], None if descr is None else descr[srt]

    @staticmethod
    def _fill_sequence(keys, vals, descr):
        r"""
        Fill bit sequence with NULL keys if bit values are not
        sequential.
        
        The instantiation of :class:`BitMask` does not include the value
        of the bit, it just assumes that the bits should be in sequence
        such that the first key has a value of 0, and the last key has a
        value of N-1.  This is a convenience function that finds the
        range of the bits and then slots in NULL keywords and empty
        descriptions where necessary to fill in the full complement of
        bits.  NULL keywords are ignored by the :class:`BitMask` object.

        Args:
            keys (:obj:`list`, :obj:`str`):
                Bit names
            vals (:obj:`list`, :obj:`int`):
                Bit values
            descr (:obj:`list`, :obj:`str`, optional):
                The description of each bit. If None, no bit
                descriptions are defined.

        Returns:
            `numpy.ndarray`_: Three arrays with the filled keys, values,
            and descriptions.

        Raises:
            ValueError: Raised if a bit value is less than 0.
        """
        _keys = numpy.atleast_1d(keys).ravel()
        _vals = numpy.atleast_1d(vals).ravel()
        _descr = None if descr is None else numpy.atleast_1d(descr).ravel()
        
        if numpy.amin(_vals) < 0:
            raise ValueError('No bit cannot be less than 0!')
        minv = numpy.amin(_vals)
        maxv = numpy.amax(_vals)

        if minv != 0 or maxv != len(_vals)-1:
            diff = list(set(numpy.arange(maxv)) - set(_vals))
            _vals = numpy.append(_vals, diff)
            _keys = numpy.append(_keys, numpy.array(['NULL']*len(diff)))
            if _descr is not None:
                _descr = numpy.append(_descr, numpy.array(['']*len(diff)))

        return _keys, _vals, _descr

    def keys(self):
        """
        Return a list of the bit keywords.

        Keywords are sorted by their bit value and 'NULL' keywords are
        ignored.
        
        Returns:
            list: List of bit keywords.
        """
        k = numpy.array(list(self.bits.keys()))
        return k[[_k != 'NULL' for _k in k]].tolist()

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

    def minimum_dtype(self, asuint=False):
        """
        Return the smallest int datatype that is needed to contain all
        the bits in the mask.  Output as an unsigned int if requested.

        Args:
            asuint (:obj:`bool`, optional):
                Return an unsigned integer type.  Signed types are
                returned by default.

        .. warning::
            uses int16 if the number of bits is less than 8 and
            asuint=False because of issue astropy.io.fits has writing
            int8 values.
        """
        if self.nbits < 8:
            return numpy.uint8 if asuint else numpy.int16
        if self.nbits < 16:
            return numpy.uint16 if asuint else numpy.int16
        if self.nbits < 32:
            return numpy.uint32 if asuint else numpy.int32
        return numpy.uint64 if asuint else numpy.int64

    def flagged(self, value, flag=None):
        """
        Determine if a bit is on in the provided bitmask value.  The
        function can be used to determine if any individual bit is on or
        any one of many bits is on.

        Args:
            value (int, array-like):
                Bitmask value.  It should be less than or equal to
                :attr:`max_value`; however, that is not checked.
            flag (str, array-like, optional):
                One or more bit names to check.  If None, then it checks
                if *any* bit is on.
        
        Returns:
            bool: Boolean flags that the provided flags (or any flag) is
            on for the provided bitmask value.  Shape is the same as
            `value`.

        Raises:
            KeyError: Raised by the dict data type if the input *flag*
                is not one of the valid :attr:`flags`.
            TypeError: Raised if the provided *flag* does not contain
                one or more strings.
        """
        _flag = self._prep_flags(flag)

        out = value & (1 << self.bits[_flag[0]]) != 0
        if len(_flag) == 1:
            return out

        nn = len(_flag)
        for i in range(1,nn):
            out |= (value & (1 << self.bits[_flag[i]]) != 0)
        return out

    def flagged_bits(self, value):
        """
        Return the list of flagged bit names for a single bit value.

        Args:
            value (int):
                Bitmask value.  It should be less than or equal to
                :attr:`max_value`; however, that is not checked.
        
        Returns:
            list: List of flagged bit value keywords.

        Raises:
            KeyError:
                Raised by the dict data type if the input *flag* is not
                one of the valid :attr:`flags`.
            TypeError:
                Raised if the provided *flag* does not contain one or
                more strings.
        """
        if not numpy.issubdtype(type(value), numpy.integer):
            raise TypeError('Input must be a single integer.')
        if value <= 0:
            return []
        keys = numpy.array(self.keys())
        indx = numpy.array([1<<self.bits[k] & value != 0 for k in keys])
        return list(keys[indx])

    def toggle(self, value, flag):
        """
        Toggle a bit in the provided bitmask value.

        Args:
            value (int, array-like):
                Bitmask value.  It should be less than or equal to
                :attr:`max_value`; however, that is not checked.
            flag (str, array-like):
                Bit name(s) to toggle.

        Returns:
            array-like: New bitmask value after toggling the selected
            bit.

        Raises:
            ValueError:
                Raised if the provided flag is None.
        """
        if flag is None:
            raise ValueError('Provided bit name cannot be None.')

        _flag = self._prep_flags(flag)

        out = value ^ (1 << self.bits[_flag[0]])
        if len(_flag) == 1:
            return out

        nn = len(_flag)
        for i in range(1,nn):
            out ^= (1 << self.bits[_flag[i]])
        return out

    def turn_on(self, value, flag):
        """
        Ensure that a bit is turned on in the provided bitmask value.

        Args:
            value (int, array-like):
                Bitmask value.  It should be less than or equal to
                :attr:`max_value`; however, that is not checked.
            flag (str, array-like):
                Bit name(s) to turn on.
        
        Returns:
            uint: New bitmask value after turning on the selected bit.

        Raises:
            ValueError:
                Raised if the provided flag is None.
        """
        if flag is None:
            raise ValueError('Provided bit name cannot be None.')

        _flag = self._prep_flags(flag)

        out = value | (1 << self.bits[_flag[0]])
        if len(_flag) == 1:
            return out

        nn = len(_flag)
        for i in range(1,nn):
            out |= (1 << self.bits[_flag[i]])
        return out

    def turn_off(self, value, flag):
        """
        Ensure that a bit is turned off in the provided bitmask value.

        Args:
            value (int, array-like):
                Bitmask value.  It should be less than or equal to
                :attr:`max_value`; however, that is not checked.
            flag (str, array-like):
                Bit name(s) to turn off.
        
        Returns:
            uint: New bitmask value after turning off the selected bit.

        Raises:
            ValueError:
                Raised if the provided flag is None.
        """
        if flag is None:
            raise ValueError('Provided bit name cannot be None.')

        _flag = self._prep_flags(flag)

        out = value & ~(1 << self.bits[_flag[0]])
        if len(_flag) == 1:
            return out

        nn = len(_flag)
        for i in range(1,nn):
            out &= ~(1 << self.bits[_flag[i]])
        return out

    def consolidate(self, value, flag_set, consolidated_flag):
        """
        Consolidate a set of flags into a single flag.
        """
        indx = self.flagged(value, flag=flag_set)
        value[indx] = self.turn_on(value[indx], consolidated_flag)
        return value

    def unpack(self, value, flag=None):
        """
        Construct boolean arrays with the selected bits flagged.

        Args:
            value (`numpy.ndarray`_):
                The bitmask values to unpack.
            flag (:obj:`str`, :obj:`list`, optional):
                The specific bits to unpack.  If None, all values are
                unpacked.
        Returns:
            tuple: A tuple of boolean numpy.ndarrays flagged according
            to each bit.
        """
        _flag = self._prep_flags(flag)
        return tuple([self.flagged(value, flag=f) for f in _flag])

    def to_header(self, hdr, prefix=None, quiet=False):
        """
        Write the bits to a fits header.

        .. todo::
            - This is very similar to the function in ParSet.  Abstract
              to a general routine?
            - The comment might have a limited length and be truncated.

        Args:
            hdr (`astropy.io.fits.Header`):
                Header object for the parameters. Modified in-place.
            prefix (:obj:`str`, optional):
                Prefix to use for the header keywords, which
                overwrites the string defined for the class. If None,
                uses the default for the class.
            quiet (:obj:`bool`, optional):
                Suppress print statements.
        """
        if prefix is None:
            prefix = self.prefix
        maxbit = max(list(self.bits.values()))
        ndig = int(numpy.log10(maxbit))+1 
        for key, value in sorted(self.bits.items(), key=lambda x:(x[1],x[0])):
            if key == 'NULL':
                continue
            hdr['{0}{1}'.format(prefix, str(value).zfill(ndig))] = (key, self.descr[value])

    @classmethod
    def from_header(cls, hdr, prefix=None):
        """
        Instantiate the BitMask using data parsed from a fits header.

        .. todo::
            - This is very similar to the function in ParSet.  Abstract
              to a general routine?
            - If comments are truncated by the comment line length,
              they'll be different than a direct instantiation.

        Args:
            hdr (`astropy.io.fits.Header`):
                Header object with the bits.
            prefix (:obj:`str`, optional):
                Prefix of the relevant header keywords, which
                overwrites the string defined for the class. If None,
                uses the default for the class.
        """
        if prefix is None:
            prefix = cls.prefix
        # Parse the bits from the header
        keys, values, descr = cls.parse_bits_from_hdr(hdr, prefix)
        # Fill in any missing bits
        keys, values, descr = cls._fill_sequence(keys, values, descr=descr)
        # Make sure the bits are sorted
        srt = numpy.argsort(values)
        # Instantiate the BitMask
        return cls(keys[srt], descr=descr[srt])

    @staticmethod
    def parse_bits_from_hdr(hdr, prefix):
        """
        Parse bit names, values, and descriptions from a fits header.

        .. todo::
            - This is very similar to the function in ParSet.  Abstract
              to a general routine?

        Args:
            hdr (`astropy.io.fits.Header`):
                Header object with the bits.
            prefix (:obj:`str`):
                The prefix used for the header keywords.
        
        Returns:
            Three lists are returned providing the bit names, values,
            and descriptions.
        """
        keys = []
        values = []
        descr = []
        for k, v in hdr.items():
            # Check if this header keyword starts with the required
            # prefix
            if k[:len(prefix)] == prefix:
                try:
                    # Try to convert the keyword without the prefix
                    # into an integer. Bits are 0 indexed and written
                    # to the header that way.
                    i = int(k[len(prefix):])
                except ValueError:
                    # Assume the value is some other random keyword that
                    # starts with the prefix but isn't a parameter
                    continue

                # Assume we've found a bit entry. Parse the bit name
                # and description and add to the compiled list
                keys += [v]
                values += [i]
                descr += [hdr.comments[k]]
        return keys, values, descr

    def to_rst_table(self, header=True, class_link=True):
        """
        Construct a reStructuredText table describing the bitmask.

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
        keys = self.keys()
        data_table = numpy.empty((len(keys)+1, 3), dtype=object)
        data_table[0,:] = ['Key', 'Bit', 'Description']
        for i,k in enumerate(keys):
            data_table[i+1,0] = k
            data_table[i+1,1] = str(self.bits[k])
            data_table[i+1,2] = '..' if self.descr is None else self.descr[self.bits[k]]

        output = []
        if header:
            output += [ '{0} Bits'.format(type(self).__name__) ]
            output += [ '-'*len(output[0]) ]
            output += [ '' ]
        if class_link:
            output += ['Class Instantiation: ' + ParSet._rst_class_name(self)]
            output += ['']
        output += [ParSet._data_table_string(data_table, delimiter='rst')]
        output += ['']
        return output

