# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
r"""
Defines a class with common functions for MaNGA fits files.  This class
needs to be as core as possible (little to no dependencies are
higher-level classes).

----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import sys
import os
import numpy
import logging
import warnings

from IPython import embed

from scipy import sparse

from astropy.wcs import WCS
from astropy.io import fits
import astropy.constants

from mangadap import __version__

from .bitmask import BitMask
from .log import log_output
from .fileio import compress_file, create_symlink
from .pixelmask import SpectralPixelMask

class DAPFitsUtil:
    r"""
    An abstract class that only implements a set of static methods for
    interaction with MaNGA fits files.

    .. todo::
        - Make this a base class for all the fits files used/created by
          the DAP?
        - Use HDUList as a base class for this object?

    """
    @staticmethod
    def storage_mode_options():
        """
        Return the allowed storage modes.

        Returns:
            list: List of the allowed storage mode for DAP fits files.
        """
        return [ '2d', '3d' ]


    @staticmethod
    def get_spatial_shape(shape, dispaxis):
        """
        
        Return the shape of the spatial axes given the specified
        dispersion axis.

        Args:
            shape (tuple): Shape of the full numpy.ndarray object.
            dispaxis (int): Axis with the spectra.

        Returns:
            tuple : Shape of the remaining axes.
        """
        return shape[:dispaxis] + shape[dispaxis+1:]


    @staticmethod
    def clean_dap_primary_header(hdr):
        # Remove some keys that are incorrect for DAP data
        hdr.remove('BSCALE', ignore_missing=True)
        hdr.remove('BZERO', ignore_missing=True)
        hdr.remove('BUNIT', ignore_missing=True)
        hdr.remove('MASKNAME', ignore_missing=True)
        return hdr


    @staticmethod
    def initialize_dap_primary_header(cube, maskname=None):
        # Copy the from the DRP and clean it
        hdr = cube.prihdr.copy()
        hdr = DAPFitsUtil.clean_dap_primary_header(hdr)

        # Change MASKNAME
        if maskname is not None:
            hdr['MASKNAME'] = maskname

        # Add versioning
        hdr['VERSPY'] = ('.'.join([ str(v) for v in sys.version_info[:3]]), 'Python version')
        hdr['VERSNP'] = (numpy.__version__, 'Numpy version')
        import scipy
        hdr['VERSSCI'] = (scipy.__version__, 'Scipy version')
        import astropy
        hdr['VERSAST'] = (astropy.__version__, 'Astropy version')
        import pydl
        hdr['VERSPYDL'] = (pydl.__version__, 'pydl version')
        import ppxf
        hdr['VERSPPXF'] = (ppxf.__version__, 'pPXF version')
        try:
            import vorbin
            hdr['VERSVBIN'] = (vorbin.__version__, 'VorBin version (if used)')
        except:
            pass
        hdr['VERSDAP'] = (__version__, 'MaNGA DAP version')

        return hdr


    @staticmethod
    def finalize_dap_header(hdr, ext, bunit=None, hduclas2='DATA', err=False, qual=False,
                            multichannel=False, bit_type=None, prepend=True, channel_names=None,
                            channel_units=None):

        _hdr = hdr.copy() if channel_names is None else \
                    DAPFitsUtil.add_channel_names(hdr, channel_names, units=channel_units)

        if bunit is not None:
            _hdr['BUNIT'] = (bunit, 'Unit of pixel value')

        # Add the common HDUCLASS keys
        _hdr['HDUCLASS'] = ('SDSS', 'SDSS format class')
        _hdr['HDUCLAS1'] = ('CUBE' if multichannel else 'IMAGE', 'Data format')
        if hduclas2 == 'DATA':
            _hdr['HDUCLAS2'] = 'DATA'
            if err:
                _hdr['ERRDATA'] = (ext+'_IVAR' if prepend else 'IVAR',
                                    'Associated inv. variance extension')
            if qual:
                _hdr['QUALDATA'] = (ext+'_MASK' if prepend else 'MASK',
                                    'Associated quality extension')
            return _hdr

        if hduclas2 == 'ERROR':
            _hdr['HDUCLAS2'] = 'ERROR'
            _hdr['HDUCLAS3'] = ('INVMSE', 'Value is inverse mean-square error')
            _hdr['SCIDATA'] = (ext, 'Associated data extension')
            if qual:
                _hdr['QUALDATA'] = (ext+'_MASK' if prepend else 'MASK',
                                    'Associated quality extension')
            return _hdr

        if hduclas2 == 'QUALITY':
            _hdr['HDUCLAS2'] = 'QUALITY'
            _hdr['HDUCLAS3'] = DAPFitsUtil._mask_data_type(bit_type)
            _hdr['SCIDATA'] = (ext, 'Associated data extension')
            if err:
                _hdr['ERRDATA'] = (ext+'_IVAR' if prepend else 'IVAR',
                                    'Associated inv. variance extension')
            return _hdr
            
        raise ValueError('HDUCLAS2 must be DATA, ERROR, or QUALITY.')


    @staticmethod
    def clean_map_header(hdr, multichannel=False):

        # KHRR made some changes here
        # Change header keywords to the default values for the third axis
        if multichannel:
            hdr['NAXIS'] = 3
            hdr.remove('CTYPE3')
            hdr.remove('CUNIT3')
            hdr['CTYPE3'] = ' '
            hdr['CUNIT3'] = ' '
            hdr['CRPIX3'] = 1
            hdr['CRVAL3'] = 1.
            hdr['CD3_3']  = 1.
            w = WCS(header=hdr)

        else:
            #hdr['NAXIS'] = 2
            #hdr.remove('NAXIS3')
            #hdr.remove('CTYPE3')
            #hdr.remove('CUNIT3')
            #hdr.remove('CRPIX3')
            #hdr.remove('CRVAL3')
            #hdr.remove('CD3_3')
            w_tmp = WCS(header=hdr)
            w = w_tmp.dropaxis(2)


        # Remove everything but the WCS information
        #w = WCS(header=hdr)
        hdr = w.to_header().copy()


        # KHRR - the DATE-OBS keyword is not in MUSE
        if('DATE-OBS' in hdr):
            hdr.comments['DATE-OBS'] = 'Date of median exposure'
            hdr.comments['MJD-OBS'] = '[d] MJD for DATE-OBS'

        # Add back in the BSCALE and BZERO values
        hdr['BSCALE'] = 1.
        hdr['BZERO'] = 0.

        # Add back in the default CTYPE3, CUNIT3
        if multichannel:
            hdr['CTYPE3'] = (' ', 'Undefined type')
            hdr['CUNIT3'] = (' ', 'Undefined units')

        return hdr


    @staticmethod
    def build_map_header(hdr, author, multichannel=False, maskname=None):
        hdr = hdr.copy()
        hdr = DAPFitsUtil.clean_map_header(hdr, multichannel=multichannel)
        hdr['AUTHOR'] = author
        if maskname is not None:
            hdr['MASKNAME'] = maskname
        return hdr


    @staticmethod
    def add_channel_names(hdr, names, units=None):
        _hdr = hdr.copy()
        nchannels = len(names)
        ndig = int(numpy.log10(nchannels))+1
        if units is not None and len(units) != nchannels:
            raise ValueError('Length of column names and units are not the same.')

        for i in range(nchannels):
            _hdr['C'+'{0}'.format(i+1).zfill(ndig)] = (names[i], 'Data in channel {0}'.format(i+1))
            if units is not None:
                _hdr['U'+'{0}'.format(i+1).zfill(ndig)] \
                            = (units[i], 'Units of data in channel {0}'.format(i+1))
        return _hdr


    @staticmethod
    def clean_cube_header(hdr):

        # Remove everything but the WCS information
        w = WCS(header=hdr)
        hdr = w.to_header().copy()

        # Fix the DATE-OBS keyword:
        # KHRR added if
        if ('DATE-OBS' in hdr):
            hdr.comments['DATE-OBS'] = 'Date of median exposure'

        if('MJD-OBS' in hdr):
            hdr.comments['MJD-OBS'] = '[d] MJD for DATE-OBS'

        # Add back in the BSCALE and BZERO values; BUNIT added during
        # "finalize"
        hdr['BSCALE'] = 1.
        hdr['BZERO'] = 0.

        return hdr


#    @staticmethod
#    def build_cube_header(drpf, author, maskname=None):
#        hdr = drpf.hdu['FLUX'].header.copy()
#        hdr = DAPFitsUtil.clean_cube_header(hdr)
#        hdr['AUTHOR'] = author
#        if maskname is not None:
#            hdr['MASKNAME'] = maskname
#        return hdr


    @staticmethod
    def build_cube_header(cube, author, maskname=None):
        hdr = cube.fluxhdr.copy()
        hdr = DAPFitsUtil.clean_cube_header(hdr)
        hdr['AUTHOR'] = author
        if maskname is not None:
            hdr['MASKNAME'] = maskname
        return hdr


    @staticmethod
    def downselect_bins(bin_indx, good_bins):
        """
        Replace values in bin_indx with -1 if they are not in good_bins.
        """
        unique_bins, reconstruct = numpy.unique(bin_indx, return_inverse=True)
        missing_bins = list(set(unique_bins) - set(good_bins))
        unique_bins = numpy.array([ -1 if x in missing_bins else x for x in unique_bins])
        return unique_bins[reconstruct]


    @staticmethod
    def match_image_naxis(hdu, ext=None):
        """
        Match the NAXISn header keywords to the shape of the image data.
        That is, reset the header keywrods such that the shape of the
        data array is (NAXIS1, NAXIS2, NAXIS3, ...). 

        Note that any reordering of the NAXISn keywords will *not*
        change the WCS information.  You will have to account for if you
        need the WCS information.

        .. todo::
            If a WCS is detected, use WCS.swapaxes() to transpose the
            WCS system!
        """
        _ext = [ h.name for h in hdu ] if ext is None else ext
        for e in _ext:
            if isinstance(hdu[e], (fits.hdu.image.PrimaryHDU, fits.hdu.image.ImageHDU)) \
                    and hdu[e].data is not None:
                shape = hdu[e].data.shape
                for i in range(len(shape)):
                    hdu[e].header['NAXIS{0}'.format(i+1)] = shape[i]
        return hdu


    @staticmethod
    def transpose_image_data(hdu, ext=None):
        """
        Transpose all image data in an HDUList.

        The primary use of this function is to ensure the HDUList data
        arrays use C-contiguous memory but with the natural array
        ordering of the axes that follow the NAXIS keywords in the
        header.  I.e., data arrays provided by astropy.io.fits have
        shapes like (..., NAXIS3, NAXIS2, NAXIS1).  This function
        reorders the data arrays to have shapes like (NAXIS1, NAXIS2,
        NAXIS3, ...), and ensures that the NAXISn header keywords
        reflect this.  Typically, one would do::

            hdu = transpose_image_data(fits.open('test.fits'))

        To prepare the image data for writing to a new fits file, one
        would do::
        
            transpose_image_data(hdu).writeto('new.fits')

        If you just want to ensure that the NAXISn keywords refect the
        python native ordering, then use :func:`match_image_naxis`.

        Note that any WCS information will **not** be transposed, which
        you will have to account for if you need the WCS information.

        .. todo::
            If a WCS is detected, use WCS.swapaxes() to transpose the
            WCS system!

        """
        _ext = [ h.name for h in hdu ] if ext is None else ext
        for e in _ext:
            if isinstance(hdu[e], (fits.hdu.image.PrimaryHDU, fits.hdu.image.ImageHDU)) \
                    and hdu[e].data is not None:
                hdu[e].data = numpy.asarray(hdu[e].data.T, order='C')
                shape = hdu[e].data.shape
                for i in range(len(shape)):
                    hdu[e].header['NAXIS{0}'.format(i+1)] = shape[i]
        return hdu


    @staticmethod
    def reconstruct_map(map_shape, bin_indx, arr, dtype=None, quiet=False):
        r"""
        Reconstruct a set of maps with the specified shape based on a
        set of input 1D arrays.  The map shape is expected to be
        :math:`(N_x,N_y)` and the array shapes are expected to be
        :math:`(N_{\rm pix},)`, where :math:`N_{\rm pix}` does *not*
        have to match the number of unique bins.  The bin_indx array is
        expected to be a vector with length :math:`N_x\timesN_y`.

        The data type of the output arrays is set by the input type if
        the provided dtype is None; otherwise, the code attempts to case
        each output map according to each element of dtype.

        Returns a tuple if the number of arrays is larger than one so
        that they can be caught by individual variables.
        """
        # Check the input types
        if not isinstance(map_shape, tuple):
            raise TypeError('Map shape must be a tuple.')
        if not isinstance(bin_indx, numpy.ndarray):
            raise TypeError('Input bin_indx must be a numpy.ndarray.')
        if not isinstance(arr, (list, tuple, numpy.ndarray)):
            raise TypeError('Input arr must be a list or tuple of numpy.ndarray objects or a '
                            'single numpy.ndarray.')
        _arr = [arr] if isinstance(arr, numpy.ndarray) else list(arr)
        if dtype is None:
            _dtype = []
            for a in _arr:
                _dtype += [a.dtype.name]
        else:
            _dtype = [dtype] if isinstance(dtype, str) else list(dtype)

        # Check the input sizes
        if len(bin_indx.shape) != 1:
            raise ValueError('Input bin_indx must be a 1D vector.')
        if bin_indx.size != numpy.prod(map_shape):
            raise ValueError('Input bin_indx length does not match cube shape.')
        narr = len(_arr)
        if len(_dtype) != narr:
            raise ValueError('Number of dtypes must match the number of arrays.')
        for i in range(1,narr):
            if _arr[i].shape != _arr[i-1].shape:
                raise ValueError('All arrays must have the same shape.')

        # Get the unique bins and how to reconstruct the bins from the
        # unique set
        unique_bins, reconstruct = numpy.unique(bin_indx, return_inverse=True)
#        print(unique_bins)

        # Handle missing bins by changing from the bin id to the bin index.
        # Detecting whether or not there is a need for handling missing bins is
        # done by comparing the unique array to a full string of all indices up
        # to the maximum index in the unique array.  **This assumes bins can
        # only either be -1, to indicate that the spaxel/bin was not analyzed,
        # or a non-negative index number** (i.e., all bin IDs must be >= -1).
        # The code handles cases both with and without any ignored bins.
        warn = False
        if numpy.any(unique_bins < 0):
            # Includes ignored bins/spaxels
            if unique_bins.size != unique_bins[-1]+2 \
                    or numpy.any((unique_bins - numpy.arange(-1,unique_bins.size-1)) != 0):
                warn = True
                unique_bins = numpy.arange(-1,unique_bins.size-1)
        else:
            # All bins/spaxels have valid bin numbers
            if unique_bins.size != unique_bins[-1]+1 \
                    or numpy.any((unique_bins - numpy.arange(unique_bins.size)) != 0):
                warn = True
                unique_bins = numpy.arange(unique_bins.size)
        if warn and not quiet:
            warnings.warn('Bin numbers and indices do not match.  Map values are expected to be '
                          'sorted by their bin number.')

#        if unique_bins.size != unique_bins[-1]+2 \
#                or numpy.any((unique_bins - numpy.arange(-1,unique_bins.size-1)) != 0):
#            if not quiet:
#                warnings.warn('Bin numbers and indices do not match.  Map values are expected '
#                              'to be sorted by their bin number.')
#            unique_bins = numpy.arange(-1,unique_bins.size-1)

        # Get the valid bins
        indx = bin_indx > -1
#        print(map_shape, bin_indx.shape, numpy.sum(indx), numpy.prod(bin_indx.shape))

        # Restructure all the arrays to match the cube shape
        new_arr = numpy.empty(narr, dtype=object)
#        print('narr:', narr)
        for i in range(narr):
            new_arr[i] = numpy.zeros(map_shape, dtype=_dtype[i])
            new_arr[i].ravel()[indx] = _arr[i][unique_bins[reconstruct[indx]]].copy()
            s = numpy.sum(numpy.invert(numpy.isfinite(new_arr[i])))
            if s > 0:
                print(i+1, narr)
                print('dtype: ', _dtype[i])
                a = numpy.sum(numpy.invert(numpy.isfinite(_arr[i][unique_bins[reconstruct[indx]]])))
                print('not finite input:', a)
                print('not finite output ravel:',
                      numpy.sum(numpy.invert(numpy.isfinite(new_arr[i].ravel()[indx]))))
                print('not finite output:', s)
                nf = numpy.invert(numpy.isfinite(new_arr[i].ravel()[indx]))
                print('inp: ', _arr[i][unique_bins[reconstruct[indx]]][nf])
                print('out: ', new_arr[i].ravel()[indx][nf])
                new_arr[i].ravel()[indx][nf] = 0.0
                raise ValueError('NaNs in data!')

        return tuple([ a for a in new_arr]) if narr > 1 else new_arr[0]


    @staticmethod
    def reconstruct_cube(cube_shape, bin_indx, arr, dtype=None):
        r"""
        Reconstruct a set of cubes with the specified shape based on a
        set of input 2D arrays.  The cube shape is expected to be
        :math:`(N_x,N_y,N_\lamba)` and the array shapes are expected to
        be :math:`(N_{\rm spec},N_\lambda)`, where :math:`N_{\rm spec}`
        does *not* have to match the number of unique bins.  The
        bin_indx array is expected to be a vector with length
        :math:`N_x\timesN_y`.

        The data type of the output arrays is set by the input type if
        the provided dtype is None; otherwise, the code attempts to case
        each output cube accordint to each element of dtype.

        Returns a tuple if the number of arrays is larger than one so
        that they can be caught by individual variables.

        """

        # Check the input types
        if not isinstance(cube_shape, tuple):
            raise TypeError('Cube shape must be a tuple.')
        if not isinstance(bin_indx, numpy.ndarray):
            raise TypeError('Input bin_indx must be a numpy.ndarray.')
        if not isinstance(arr, (list, tuple, numpy.ndarray)):
            raise TypeError('Input arr must be a list or tuple of numpy.ndarray objects or a '
                            'single numpy.ndarray.')
        _arr = [arr] if isinstance(arr, numpy.ndarray) else list(arr)
        if dtype is None:
            _dtype = []
            for a in _arr:
                _dtype += [a.dtype.name]
        else:
            _dtype = [dtype] if isinstance(dtype, str) else list(dtype)

        # Check the input sizes
        if len(bin_indx.shape) != 1:
            raise ValueError('Input bin_indx must be a 1D vector.')
        if bin_indx.size != numpy.prod(cube_shape[0:2]):
            raise ValueError('Input bin_indx length does not match cube shape.')
        narr = len(_arr)
        if len(_dtype) != narr:
            raise ValueError('Number of dtypes must match the number of arrays.')
        nwave = _arr[0].shape[-1]
        if cube_shape[-1] != nwave:
            raise ValueError('Number of pixels in arrays must match last dimension of cube.')
        for i in range(1,narr):
            if _arr[i].shape != _arr[i-1].shape:
                raise ValueError('All arrays must have the same shape.')

        # Get the unique bins and how to reconstruct the bins from the
        # unique set
        unique_bins, reconstruct = numpy.unique(bin_indx, return_inverse=True)

        # Handle missing bins by changing from the bin id to the bin index.
        # Detecting whether or not there is a need for handling missing bins is
        # done by comparing the unique array to a full string of all indices up
        # to the maximum index in the unique array.  **This assumes bins can
        # only either be -1, to indicate that the spaxel/bin was not analyzed,
        # or a non-negative index number** (i.e., all bin IDs must be >= -1).
        # The code handles cases both with and without any ignored bins.
        warn = False
        if numpy.any(unique_bins < 0):
            # Includes ignored bins/spaxels
            if unique_bins.size != unique_bins[-1]+2 \
                    or numpy.any((unique_bins - numpy.arange(-1,unique_bins.size-1)) != 0):
                warn = True
                unique_bins = numpy.arange(-1,unique_bins.size-1)
        else:
            # All bins/spaxels have valid bin numbers
            if unique_bins.size != unique_bins[-1]+1 \
                    or numpy.any((unique_bins - numpy.arange(unique_bins.size)) != 0):
                warn = True
                unique_bins = numpy.arange(unique_bins.size)
        if warn and not quiet:
            warnings.warn('Bin numbers and indices do not match.  Map values are expected to be '
                          'sorted by their bin number.')

#        if unique_bins.size != unique_bins[-1]+2 \
#                or numpy.any((unique_bins - numpy.arange(-1,unique_bins.size-1)) != 0):
#            warnings.warn('Bin numbers and indices do not match.  Spectra are expected '
#                          'to be sorted by their bin number.')
#            unique_bins = numpy.arange(-1,unique_bins.size-1)

        # Get the valid bins
        indx = bin_indx > -1

        # Restructure all the arrays to match the cube shape
        new_arr = numpy.empty(narr, dtype=object)
        for i in range(narr):
            new_arr[i] = numpy.zeros(cube_shape, dtype=_dtype[i])
            new_arr[i].reshape(-1,nwave)[indx,:] = _arr[i][unique_bins[reconstruct[indx]],:]

        return tuple([ a for a in new_arr]) if narr > 1 else new_arr[0]


    @staticmethod
    def read(ofile, permissions='readonly', checksum=False):
        return DAPFitsUtil.transpose_image_data(fits.open(ofile, mode=permissions,
                                                          checksum=checksum))


    @staticmethod
    def write(hdu, ofile, clobber=False, checksum=False, symlink_dir=None, relative_symlink=True,
              loggers=None, quiet=False):
        """
        Write an HDUList to an output file.  It is expected that the hdu
        was either read by the :func:`read` function, or that the hdu
        was passed through :func:`transpose_image_data` such that
        HDUList data arrays use C-contiguous memory but with the natural
        array ordering of the axes that follow the NAXIS keywords in the
        header.
        """
        # Get the output file and determine if it should be compressed
        compress = False
        if ofile.split('.')[-1] == 'gz':
            _ofile = ofile[:ofile.rfind('.')] 
            compress = True
        else:
            _ofile = ofile
    
        # Transpose the order
        hdu = DAPFitsUtil.transpose_image_data(hdu)

        # Write the data
        if not quiet:
            log_output(loggers, 1, logging.INFO, 'Writing: {0}'.format(_ofile))
        hdu.writeto(_ofile, overwrite=clobber, checksum=checksum)

        # Transpose it back
        hdu = DAPFitsUtil.transpose_image_data(hdu)

        # Compress if desired
        if compress:
            if not quiet:
                log_output(loggers, 1, logging.INFO, 'Compressing: {0}'.format(ofile))
            # And compress it
            compress_file(_ofile, clobber=clobber)
            os.remove(_ofile)
    
        # Create the symlink if requested
        if symlink_dir is not None:
            create_symlink(ofile, symlink_dir, relative_symlink=relative_symlink, loggers=loggers,
                           quiet=quiet)


    @staticmethod
    def unique_bins(bin_indx, return_index=False):
        """
        Get the unique bins and the indices of the unique bins in the
        flattened spatial dimension, ignoring the bins with indices of
        -1.
        """
        uniq, indx = numpy.unique(bin_indx.ravel(), return_index=True)
        if uniq[0] == -1:
            uniq = uniq[1:]
            indx = indx[1:]
        return (uniq, indx) if return_index else uniq
        
#        return tuple(map(lambda x : x[1:], numpy.unique(bin_indx.ravel(), return_index=True))) \
#                    if return_index else numpy.unique(bin_indx.ravel())[1:]


    # TODO: Consolidate these with what's in datacube and rowstackedspectra.
    @staticmethod
    def copy_to_array(hdu, ext='FLUX', allowed_ext=None, waverange=None, select_bins=None,
                      missing_bins=None, nbins=None, unique_bins=None):
        r"""
        Return a copy of the selected data array.  The array size is
        always :math:`N_{\rm spec} \times N_{\rm wavelength}`; i.e., the
        CUBE files are flattened to two dimensions, matching the
        dimensionality of the RSS files.  The spatial positions within
        the original DRP file for each spectrum are given by tuples in
        :attr:`spatial_index`.  See :func:`copy_to_masked_array` for
        arguments.

        Returns:
            numpy.ndarray: A 2D array with a copy of the data from the
            selected extension.
        """
        masked_data = DAPFitsUtil.copy_to_masked_array(hdu, ext=ext, allowed_ext=allowed_ext,
                                                       waverange=waverange, select_bins=select_bins,
                                                       missing_bins=missing_bins, nbins=nbins,
                                                       unique_bins=unique_bins)
        # For this approach, the wavelengths masked should be
        # *identical* for all spectra
        nwave = numpy.sum(numpy.invert(numpy.ma.getmaskarray(masked_data))[0,:])
        if nwave == 0:
            raise ValueError('Full wavelength range has been masked!')
        # Compression returns a flattened array, so it needs to be
        # reshaped for the new number of (unmasked) wavelength channels
        return masked_data.compressed().reshape(-1,nwave)


    @staticmethod
    def copy_to_masked_array(hdu, ext='FLUX', mask_ext=None, flag=None, bitmask=None,
                             allowed_ext=None, waverange=None, select_bins=None, missing_bins=None,
                             nbins=None, unique_bins=None):
        r"""
        Return a copy of the selected data array as a masked array.
        This is functionally identical to :func:`DRPFits.copy_to_array`,
        except the output format is a `numpy.ma.MaskedArray`_.  The
        pixels that are considered to be masked can be specified using
        the `flag` option.
        
        Args:
            waverange (array-like) : (**Optional**) Two-element array
                with the first and last wavelength to include in the
                computation.  Default is to use the full wavelength
                range.
            ext (str) : (**Optional**) Name of the extension from which
                to draw the data.  Must be allowed for the current
                :attr:`mode`; see :attr:`spectral_arrays`.  Default is
                `'FLUX'`.
            mask_ext (str) : (**Optional**) Name of the extension with
                the mask bit data.  Must be allowed for the current
                :attr:`mode`; see :attr:`spectral_arrays`.  Default is
                `'MASK'`.
            flag (str or list): (**Optional**) (List of) Flag names that
                are considered when deciding if a pixel should be
                masked.  The names *must* be a valid bit name as defined
                by the provided bitmask (see
                :class:`mangadap.util.bitmask.BitMask`).  If not
                provided, *ANY* non-zero mask bit is omitted.
            bitmask (:class:`mangadap.util.bitmask.BitMask`):
                (**Optional**) BitMask object used to determine if pixel
                is flagged as being masked by a specific maskbit.
            allowed_ext (array-like): (**Optional**) List of allowed
                extenstions that can be copied to an array; the mask
                extension **must** be one of these.
            waverange (array-like): (**Optional**) A two-element
                array with the valid wavelength range.  Anything outside
                of this range is masked.
            select_bins (array-like): (**Optional**) A list of spectra
                in the flattened cube (with shape :math:`N_{\rm spec}
                \times N_\lambda`) to return.
            missing_bins (list): (**Optional**) A list of bin numbers
                that are missing from the output selection and that will
                be replaced with masked spectra.
            nbins (int): (**Optional**) The total number of defined
                bins.  Default is to find the maximum number in the
                unique_bins list.
            unique_bins (array-like): (**Optional**) The indices of the
                bins that have valid spectra.  The length of this array
                should match the selected_bins array.

        Returns:
            numpy.ndarray: A 2D array with a copy of the data from the
            selected extension.

        Raises:
            KeyError : Raised if `ext` is not a valid extension.
        """

        # Make sure extensions are not lists
        _ext = ext
        if isinstance(ext, list):
            if len(ext) > 1:
                raise ValueError('Must provide only one extension.')
            _ext = ext[0]
        _mask_ext = mask_ext
        if isinstance(_mask_ext, list):
            if len(_mask_ext) > 1:
                raise ValueError('Must provide only one mask extension.')
            _mask_ext = mask_ext[0]

        # Make sure operation is allowed
        if allowed_ext is not None and _ext not in allowed_ext:
            raise KeyError('Cannot copy data to array for extension {0}.'.format(_ext))
        if allowed_ext is not None and _mask_ext is not None and _mask_ext not in allowed_ext:
            raise KeyError('Cannot use extension {0} as mask.'.format(_mask_ext))

        # Make sure masks and flags make sense
        if flag is not None and bitmask is None:
            raise ValueError('If flags are provided, must also provided associated BitMask.')
        if bitmask is not None and not isinstance(bitmask, BitMask):
            raise TypeError('Provided object must be a BitMask.')

        # Total number of available spectra
        nspec = numpy.prod(hdu[_ext].data.shape[0:-1])

        # Create the wavelength mask (will be all true if
        # waverange=None)
        mask = SpectralPixelMask(waverange=waverange).boolean(hdu['WAVE'].data, nspec=nspec)

        # Add in any bitmasks
        if _mask_ext is not None and bitmask is not None:
            mask |= bitmask.flagged(hdu[_mask_ext].data.reshape(nspec,-1), flag=flag)
        if _mask_ext is not None and bitmask is None:
            mask |= (hdu[_mask_ext].data > 0)

        # Create the output MaskedArray
        a = numpy.ma.MaskedArray(hdu[_ext].data.copy().reshape(nspec,-1), mask=mask)

        # Apply any bin selection
        if select_bins is not None:
            a = a[select_bins,:]

        # No missing bins are specified, so return the array
        if missing_bins is None or len(missing_bins) == 0:
            return a

        # Missing bins have been specified, so the unique_bins also need
        # to be specified
        if unique_bins is None:
            raise ValueError('Must specify the unique bins if missing bins are specified.')
        if len(unique_bins) != a.shape[0]:
            raise ValueError('Number of unique bins does not match the number of spectra.')

        # If the number of bins hasn't been provided, construct it based
        # on the unique bins specified
        _nbins = max(unique_bins)+1 if nbins is None else nbins

        # Create the array with the right shape and mask everything
        _a = numpy.ma.MaskedArray(numpy.zeros((_nbins, a.shape[-1]), dtype=hdu[_ext].data.dtype),
                                  mask=numpy.ones((_nbins, a.shape[-1]), dtype=bool))
        # Copy the data to the correct bins, which also unmasks these
        # pixels
        _a[unique_bins,:] = a
        # Return the masked array; "missing bins" are fully masked
        return _a

    @staticmethod
    def marginalize_mask(mask, inp_flags=None, inp_bitmask=None, out_flag=None, out_bitmask=None,
                         out_mask=None, dispaxis=2):
        """
        Marginalize a mask over one of its spectral dimension.

        The output mask selects pixels that are masked by *all* pixels along the
        spectral axis.

        If no optional arguments are provided, the input mask must be a boolean
        array and the result is identically ``numpy.all(mask, axis=dispaxis)``.

        Args:
            mask (`numpy.ndarray`_):
                The mask with more than 1 dimension.  The spectral axis is
                defined by ``dispaxis``.  The array can be either a bitmask or a
                boolean mask; the latter is assumed if ``inp_bitmask`` is not
                provided, and vice versa.
            inp_flags (:obj:`str`, :obj:`list`, optional):
                The list of input flags used to select masked values in bitmask
                arrays.  If None and ``inp_bitmask`` is provided, all flags are
                check individually.
            inp_bitmask (:class:`~mangadap.util.bitmask.BitMask`, optional):
                Object used to identify bits in the input mask.
            out_flag (:obj:`str`, optional):
                If returning a bitmask, this sets the flag to use to identify
                the masked pixels.  If None and ``inp_bitmask`` is provided, the
                output flags are identical to the input flags.
            out_bitmask (:class:`~mangadap.util.bitmask.BitMask`, optional):
                Object used to identify bits in the output mask.
            out_mask (`numpy.ndarray`_, optional):
                If provided, any masking is applied directly to this array;
                i.e., this is both modified and returned by the function.  If
                provided, the array must have the correct shape (same as the
                input array wihtout the spectral axis).
            dispaxis (:obj:`int`, optional):
                The spectral axis of the array.

        Returns:
            `numpy.ndarray`_: A mask array with the marginalized mask.  The
            shape is the same as the input array, but without the spectral axis.

        Raises:
            ValueError:
                Raised if the input arrays have incorrect shapes.
        """
        # Check the input
        if not mask.ndim > 1:
            raise ValueError('Input mask array must have more than one dimension.')

        itype = bool if inp_bitmask is None else inp_bitmask.minimum_dtype()
        # TODO: Need to use 'same_kind' here instead of 'equiv' because MaNGA
        # DRP doesn't necessarily use the minimum dtype for the masks.
        if not numpy.can_cast(mask.dtype, numpy.dtype(itype), casting='same_kind'):
            raise TypeError(f'Cannot cast input mask with type {mask.dtype} to expected type '
                            f'{numpy.dtype(itype)}.')

        spatial_shape = DAPFitsUtil.get_spatial_shape(mask.shape, dispaxis)

        if out_mask is not None and out_mask.shape != spatial_shape:
            raise ValueError('Provided output mask has incorrect shape; expected '
                             f'{spatial_shape}, found {out_mask.shape}.')

        if inp_bitmask is None and inp_flags is not None:
            raise ValueError('Cannot interpret input flags without inp_bitmask.')
        if out_bitmask is None and out_flag is not None:
            raise ValueError('Cannot interpret output flags without out_bitmask.')
        if out_bitmask is not None and out_flag is None:
            raise ValueError('If providing the output bitmask, must also provide the output flag.')

        otype = bool if out_bitmask is None else out_bitmask.minimum_dtype()
        # TODO: May need to use 'same_kind' instead of 'equiv'
        if out_mask is not None \
                and not numpy.can_cast(out_mask.dtype, numpy.dtype(otype), casting='equiv'):
            raise TypeError('Provided output mask has incorrect type; expected '
                            f'{numpy.dtype(otype)}, found {out_mask.dtype}.')
        _out_mask = numpy.zeros(spatial_shape, dtype=otype) if out_mask is None else out_mask

        if inp_bitmask is None:
            # Handle the case when the input mask is a boolean array
            indx = numpy.all(mask, axis=dispaxis)
            if out_bitmask is None:
                _out_mask |= indx
                return _out_mask
            _out_mask[indx] = out_bitmask.turn_on(_out_mask[indx], out_flag)
            return _out_mask

        # Handle the case when the input mask is a bit mask
        if inp_flags is None:
            inp_flags = inp_bitmask.keys()
        _inp_flags = inp_flags if isinstance(inp_flags, list) else [inp_flags]

        for flag in _inp_flags:
            indx = numpy.all(inp_bitmask.flagged(mask, flag=flag), axis=dispaxis)
            if numpy.sum(indx) == 0:
                continue
            if out_bitmask is None:
                _out_mask |= indx
                continue
            _out_mask[indx] = out_bitmask.turn_on(_out_mask[indx],
                                                  (flag if out_flag is None else out_flag))
        return _out_mask


    @staticmethod
    def _mask_data_type(bit_type):
        if bit_type in [numpy.uint64, numpy.int64]:
            return ('FLAG64BIT', '64-bit mask')
        if bit_type in [numpy.uint32, numpy.int32]:
            return ('FLAG32BIT', '32-bit mask')
        if bit_type in [numpy.uint16, numpy.int16]:
            return ('FLAG16BIT', '16-bit mask')
        if bit_type in [numpy.uint8, numpy.int8]:
            return ('FLAG8BIT', '8-bit mask')
        if bit_type == numpy.bool:
            return ('MASKZERO', 'Binary mask; zero values are good/unmasked')
        raise ValueError('Invalid bit_type: {0}!'.format(str(bit_type)))


    @staticmethod
    def empty_hdus(ext):
        hdus = []
        for e in ext:
            hdus += [ fits.ImageHDU(data=None, name=e) ]
        return hdus


    @staticmethod
    def list_of_image_hdus(data, hdr, ext):
        hdus = []
        for d,h,e in zip(data,hdr,ext):
            hdus += [ fits.ImageHDU(data=d, header=h, name=e) ]
        return hdus


    @staticmethod
    def redshift_to_Newtonian_velocity(v, redshift, ivar=False):
        if ivar:
            return v*numpy.square(1.0+redshift)
        return (v-astropy.constants.c.to('km/s').value*redshift)/(1.0+redshift)


