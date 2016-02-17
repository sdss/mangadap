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

    Here is an example of reading/creating a TemplateLibrary, then
    creating a plot that shows the full spectrum (blue) and the unmasked
    spectrum (green)::

        # Imports
        import numpy
        from mangadap.drpfile import drpfile
        from mangadap.proc.TemplateLibrary import TemplateLibrary
        from mangadap.util.bitmasks import TemplateLibraryBitMask
        from matplotlib import pyplot

        # Define the DRP file
        drpf = drpfile(7495, 12703, 'CUBE')

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
        :class:`TemplateLibraryBitMask` to new format where the bits are
        read from a configuration file.

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
import os
import textwrap
from mangadap.config.util import _read_dap_mask_bits
from mangadap.util.defaults import default_dap_source

__author__ = 'Kyle B. Westfall'

class BitMask:
    """
    Generic class to handle and manipulate bitmasks.

    Args:
        keys (str, list, numpy.ndarray): List of keys (or single key) to
            use as the bit name.  Each key is given a bit number ranging
            from 0..N-1.
        descr (str, list, numpy.ndarray): (Optional) List of descriptions
            (or single discription) provided by :func:`info` for each bit.

    Raises:
        ValueError: Raised if the provided `keys` are not all unique or
            if the number of descriptions provided is not the same as
            the number of keys.
        TypeError: Raised if the provided `keys` do not have the correct
            type.

    Attributes:
        bits (dict): A dictionary with the bit name and value
        descr (numpy.ndarray) : List of bit descriptions
        max_value (int): The maximum valid bitmask value given the
            number of bits.

    """
    def __init__(self, keys, descr=None):
        if type(keys) not in [str, list, numpy.ndarray]:
            raise TypeError('Input bit names do not have the right type.')

        if type(keys) == list:
            keys = numpy.array(keys)
        if type(keys) == str:
            keys = numpy.array([keys])

        if keys.size != numpy.unique(keys).size:
            raise ValueError('All input keys must be unique.')

        self.bits = dict([ (k,i) for i,k in enumerate(keys) ])
        self.max_value = (1 << self.nbits())-1

        self.descr = None
        if descr is None:
            return

        if type(descr) == list:
            descr = numpy.array(descr)
        if type(descr) == str:
            descr = numpy.array([descr])

        if descr.size != self.nbits():
            raise ValueError('Number of listed descriptions not the same as number of keys.')

        self.descr = descr.copy()


    def keys(self):
        """
        Return a list of the bits.
        
        Returns:
            list : List of bit keywords.
        """
        return list(self.bits.keys())


    def nbits(self):
        """
        Return the number of bits.

        Returns:
            int : Number of bits
        """
        return len(self.bits)


    def info(self):
        """
        Print the list of bits and, if available, their descriptions.
        """
        tr, tcols = numpy.array(os.popen('stty size', 'r').read().split()).astype(int)
        tcols -= int(tcols*0.1)
        for k,v in sorted(self.bits.items(), key=lambda x:(x[1],x[0])):
            print('Bit: {0} = {1}'.format(k,v))
            if self.descr is not None:
                print(textwrap.fill('Description: {0}'.format(self.descr[v]), tcols))
            print(' ')


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
            flag = self.keys()

        if type(flag) == str:
            return value & (1 << self.bits[flag]) != 0

        if hasattr(flag, "__iter__"):
            if any([ (type(f) != str and type(f) != numpy.str_) for f in flag ]):
                raise Exception('Provided bit names must be strings!')
            out = value & (1 << self.bits[flag[0]]) != 0
#            print(out)
            nn = len(flag)
#            print(nn)
            for i in range(1,nn):
                out |= (value & (1 << self.bits[flag[i]]) != 0)
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
        return value ^ (1 << self.bits[flag])


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
        return value | (1 << self.bits[flag])


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
        return value & ~(1 << self.bits[flag])


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
    def __init__(self, dapsrc=None):
        keys, descr = _read_dap_mask_bits('template_bits.ini',
                                          default_dap_source() if dapsrc is None else str(dapsrc))
        BitMask.__init__(self, keys, descr)


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



