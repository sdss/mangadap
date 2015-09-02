"""

Bitmask class

*Source location*:
    $MANGADAP_DIR/python/mangadap/bitmask.py

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

*Class usage examples*:

    .. todo::

        Add some usage comments here!

*Revision history*:
    | **01 Jun 2015**: Original implementation by K. Westfall (KBW)

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

__author__ = 'Kyle B. Westfall'

class BitMask:
    """
    Generic class to handle and manipulate bitmasks.

    Args:
        keys (str, list, numpy.ndarray): List of keys (or single key) to
            use as the bit name.  Each key is given a bit number ranging
            from 0..N-1.

    Raises:
        Exception: Raised if the provided *keys* do not have the right
            type or if the provided keys are not all unique.

    Attributes:
        keys (numpy.ndarray): List of keys (or single key) to use as the
            bit name.  Each key is given a bit number ranging from
            0..N-1.
        nbits (int): The number of bits (i.e. the length of the
            :attr:`keys` vector).
        flags (dict): A dictionary with the bit name and value
        max_value (int): The maximum valid bitmask value given the
            number of bits.

    """
    def __init__(self, keys):
        if type(keys) not in [str, list, numpy.ndarray]:
            raise Exception('Input bit names do not have the right type.')

        if type(keys) == list:
            keys = numpy.array(keys)
        if type(keys) == str:
            keys = numpy.array([keys])

        if keys.size != numpy.unique(keys).size:
            raise Exception('All input keys must be unique.')

        self.keys = keys
        self.nbits = self.keys.size
        self.flags = dict([ (self.keys[i],i) for i in range(0,self.nbits) ])
        self.max_value = (1 << self.nbits)-1


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
            bool : Boolean flags that the provided flags (or any flag)
            is on for the provided bitmask value.

        Raises:
            KeyError: Raised by the dict data type if the input *flag*
                is not one of the valid :attr:`flags`.
            Exception: Raised if the provided *flag* does not contain
                one or more strings.

        """
        if flag is None:
            flag = self.keys

        if type(flag) == str:
            return value & (1 << self.flags[flag]) != 0

        if hasattr(flag, "__iter__"):
            if any([ (type(f) != str and type(f) != numpy.str_) for f in flag ]):
                raise Exception('Provided bit names must be strings!')
            out = value & (1 << self.flags[flag[0]]) != 0
#            print(out)
            nn = len(flag)
#            print(nn)
            for i in range(1,nn):
                out |= (value & (1 << self.flags[flag[i]]) != 0)
#                print(i, out)
            return out

        raise Exception('Provided bit name must be a string or an iterable of strings!')


    def toggle(self, value, flag):
        """
        Toggle a bit in the provided bitmask value.

        Args:
            value (uint or array): Bitmask value.  It should be less
                than or equal to :attr:`max_value`; however, that is not
                checked.
            flag (string): Bit name to toggle.
        
        Returns:
            uint : New bitmask value after toggling the selected bit.

        Raises:
            KeyError: Raised by the dict data type if the input *flag*
                is not one of the valid :attr:`flags`.
            Exception: Raised if the provided *flag* is not a string.

        """
        if type(flag) != str:
            raise Exception('Provided bit name must be a string!')
        return value ^ (1 << self.flags[flag])


    def turn_on(self, value, flag):
        """
        Ensure that a bit is turned on in the provided bitmask value.

        Args:
            value (uint or array): Bitmask value.  It should be less
                than or equal to :attr:`max_value`; however, that is not
                checked.
            flag (string): Bit name to turn on.
        
        Returns:
            uint : New bitmask value after turning on the selected bit.

        Raises:
            KeyError: Raised by the dict data type if the input *flag*
                is not one of the valid :attr:`flags`.
            Exception: Raised if the provided *flag* is not a string.

        """
        if type(flag) != str:
            raise Exception('Provided bit name must be a string!')
        return value | (1 << self.flags[flag])


    def turn_off(self, value, flag):
        """
        Ensure that a bit is turned off in the provided bitmask value.

        Args:
            value (uint or array): Bitmask value.  It should be less
                than or equal to :attr:`max_value`; however, that is not
                checked.
            flag (string): Bit name to turn off.
        
        Returns:
            uint : New bitmask value after turning off the selected bit.

        Raises:
            KeyError: Raised by the dict data type if the input *flag*
                is not one of the valid :attr:`flags`.
            Exception: Raised if the provided *flag* is not a string.

        """
        if type(flag) != str:
            raise Exception('Provided bit name must be a string!')
        return value & ~(1 << self.flags[flag])


class TemplateLibraryBitMask(BitMask):
    """
    Derived class that specifies a BitMask for the template library
    data.

    The bit names and meanings are:
    
        - 'NO_DATA': Pixel has no data 

        - 'WAVE_INVALID': Used to designate pixels in the 1D spectra
          that are outside the valid wavelength range defined by
          :func:`mangadap.util.defaults.default_template_libraries`

        - 'FLUX_INVALID': Used to designate pixels in the 1D spectra
          that are below the valid flux limit defined by
          :func:`mangadap.util.defaults.default_template_libraries`

        - 'SPECRES_EXTRAP': The spectral resolution has been matched to
          a value that was an extrapolation of the target spectral
          resolution samples.
          
        - 'SPECRES_LOW': The spectral resolution was *not* matched to
          the target value because the target value was *higher* than
          the existing spectral resolution.

    """
    def __init__(self):
        BitMask.__init__(self, numpy.array([
                'NO_DATA',          # Pixel not observed
                'WAVE_INVALID',     # Designated as invalid wavelength region
                'FLUX_INVALID',     # Designated as invalid flux value
                'SPECRES_EXTRAP',   # Spectral resolution match used extrapolated resolution value
                'SPECRES_LOW'       # No spectral resolution match because target higher than input
                         ]))


def HDUList_mask_wavelengths(hdu, bitmask, bitmask_flag, wave_limits, wave_ext='WAVE', \
                             mask_ext='MASK', invert=False):
    """
    Mask pixels in a specified wavelength range by turning on the bit
    value in the specified extention in a provided HDUList object.

    Args:
        hdu (`astropy.io.fits.hdu.hdulist.HDUList`_): HDUList to alter
        bitmask (class:`BitMask`): Bit mask object used to turn on the
            named bit mask.
        bitmask_flag (str): Name of the bit to turn on.
        wave_limits (list or numpy.ndarray): Two-element array with the
            low and high wavelength limits.
        wave_ext (str): (Optional) Name of the wavelength extension in
            *hdu*.
        mask_ext (str): (Optional) Name of the mask extension in *hdu*.
        invert (bool): (Optional) Invert the sense of the masking.
            Instead of masking pixel in the wavelength interval, mask
            all pixels *outside* it.

    Returns:
        `astropy.io.fits.hdu.hdulist.HDUList`_ : The modified HDUList
            object.

    Raises:
        Exception: Raised if *wave_limits* does not have a length of
            two.

    """
    if len(wave_limits) != 2:
        raise Exception('Wavelength limits must be a two-element vector.')

    indx = numpy.where( (hdu[wave_ext].data < wave_limits[0]) \
                        | (hdu[wave_ext].data > wave_limits[1])) if invert else \
           numpy.where( (hdu[wave_ext].data >= wave_limits[0]) \
                        & (hdu[wave_ext].data <= wave_limits[1]))
    hdu[mask_ext].data[indx] = bitmask.turn_on(hdu[mask_ext].data[indx], bitmask_flag)
    return hdu



