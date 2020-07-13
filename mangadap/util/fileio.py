# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Provides a set of file I/O routines.

----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

import os
import warnings
import gzip
import shutil
import logging
import numpy

from astropy.io import fits
from .log import log_output

def wavelength_vector(npix, header, log10=False):
    """
    Return a vector with wavelength coordinates drawn from the WCS
    coordinates in the header.  The function uses CRVAL1, CRPIX1,
    CDELT1.

    Args:
        npix (:obj:`int`):
            Length of the vector in pixels.
        header (`astropy.io.fits.Header`_):
            Header with the WCS data
        log10 (:obj:`bool`, optional):
            Flag that the wavelengths are logarithmically sampled.

    Returns:
        `numpy.ndarray`_: Wavelengths of a spectrum.

    """
    crval = float(header['CRVAL1'])
    crpix = int(header['CRPIX1'])
    cdelt = float(header['CDELT1'])
    wave = (numpy.arange(1.0,npix+1) - crpix)*cdelt + crval
    return numpy.power(10., wave) if log10 else wave


def readfits_1dspec(filename, log10=False):
    """
    Read a 1D fits spectrum and return two vectors with the wavelength
    and flux.

    Args:
        filename (:obj:`str`):
            Name of the file to read.

    Returns:
        `numpy.ndarray`_: Two numpy.float64 arrays with the
        wavelength and flux read for the spectrum.

    Raises:
        ValueError:
            Raised if the input fits file has more than one extension
            or its primary extension has more than two dimensions.
    """
    hdu = fits.open(filename, mode='readonly')

    if (len(hdu)) != 1:
        warnings.warn('{0} has more than one extension.'.format(filename))
    
    if hdu[0].header['NAXIS'] != 1:
        raise ValueError('{0} is not one dimensional!'.format(filename))
    
    spec = numpy.copy(hdu[0].data).astype(numpy.float64)
    wave = wavelength_vector(spec.size, hdu[0].header, log10=log10)
    hdu.close()
    return wave, spec


def writefits_1dspec(ofile, crval1, cdelt1, flux, hdr=None, clobber=False):
    """

    Write a simple one-dimensional spectrum.

    Args:
        ofile (str): Name of the file to write.
        crval1 (float): (Log base 10 of the) Initial wavelength, which
            is included in the header with the keyword 'CRVAL1';
            'CRPIX1' is always set to 1.
        cdelt1 (float): The change in (log base 10) wavelength per
            pixel, which is included in the header with the keywords
            'CDELT1' and 'CD1_1'; 'CRPIX1' is always set to 1.
        flux (array): Vector of the flux values.
        clobber (bool): (**Optional**) Flag to overwrite any existing
            file of the same name.
    """
    _hdr = fits.Header() if hdr is None else hdr
    _hdr['CRPIX1'] = 1
    _hdr['CRVAL1'] = crval1
    _hdr['CDELT1'] = cdelt1
    _hdr['CD1_1'] = _hdr['CDELT1']
    fits.HDUList([ fits.PrimaryHDU(flux, header=_hdr) ]).writeto(ofile, overwrite=clobber)
                 #clobber=clobber)


def read_template_spectrum(filename, data_ext=0, ivar_ext=None, sres_ext=None, log10=False):
    r"""
    Read a template spectrum.
    
    Template spectra are "raw format" files with template data and are,
    at minimum, expected to have the following components::
    
        hdu[0].header['CRVAL1']
        hdu[0].header['CRPIX1']
        hdu[0].header['CDELT1']
        hdu[data_ext].data

    The latter has the flux data. If ``log10`` is true, the
    wavelength solution above is expected to be in log wavelengths.
    
    Args:
        filename (:obj:`str`):
            Name of the fits file to read.
        data_ext (:obj:`str`, optional):
            Name of the extension with the flux data.
        ivar_ext (:obj:`str`, optional):
            Name of the extension with the inverse variance data. If
            None, no inverse data are returned.
        sres_ext (:obj:`str`, optional):
            Name of the extension with the spectral resolution
            (:math:R=\lambda/\delta\lambda`) measurements. If None,
            no spectral resolution data are returned.
        log10 (:obj:`bool`, optional):
            Flag the WCS wavelength coordinates as being in base-10
            log wavelength, instead of linear. Default is to assume
            linear.

    Returns:
        `numpy.ndarray`_: Up to four vectors with the wavelength,
        flux, inverse variance (if `ivar_ext` is provided), and
        spectral resolution (if `sres_ext` is provided) of the
        template spectrum.

    Raises:
        ValueError:
            Raised if fits file is not one-dimensional.
        KeyError:
            Raised if various header keywords or extension names are
            not available.
    """
    with fits.open(filename, mode='readonly') as hdu:
        if len(hdu[data_ext].data.shape) != 1:
            raise ValueError('Spectrum in {0} is not one dimensional!'.format(filename))
    
        spec = numpy.copy(hdu[data_ext].data).astype(numpy.float64)
        wave = wavelength_vector(spec.size, hdu[0].header, log10=log10)

        # Collect objects to return
        ret = (wave, spec)
        if ivar_ext is not None:
            ret += (numpy.copy(hdu[ivar_ext].data).astype(numpy.float64), )
        if sres_ext is not None:
            ret += (numpy.copy(hdu[sres_ext].data).astype(numpy.float64), )
        return ret
    

def init_record_array(shape, dtype):
    r"""
    Utility function that initializes a record array using a provided
    input data type.  For example::

        dtype = [ ('INDX', numpy.int, (2,) ),
                  ('VALUE', numpy.float) ]

    Defines two columns, one named `INDEX` with two integers per row and
    the one named `VALUE` with a single float element per row.  See
    `numpy.recarray`_.
    
    Args:
        shape (:obj:`int`, :obj:`tuple`):
            Shape of the output array.
        dtype (:obj:`list`):
            List of the tuples that define each element in the record
            array.

    Returns:
        `numpy.recarray`_: Zeroed record array
    """
    data = numpy.zeros(shape, dtype=dtype)
    return data.view(numpy.recarray)


def rec_to_fits_type(rec_element):
    """
    Return the string representation of a fits binary table data type
    based on the provided record array element.
    """
    n = 1 if len(rec_element[0].shape) == 0 else rec_element[0].size
    if rec_element.dtype == numpy.bool:
        return '{0}L'.format(n)
    if rec_element.dtype == numpy.uint8:
        return '{0}B'.format(n)
    if rec_element.dtype == numpy.int16 or rec_element.dtype == numpy.uint16:
        return '{0}I'.format(n)
    if rec_element.dtype == numpy.int32 or rec_element.dtype == numpy.uint32:
        return '{0}J'.format(n)
    if rec_element.dtype == numpy.int64 or rec_element.dtype == numpy.uint64:
        return '{0}K'.format(n)
    if rec_element.dtype == numpy.float32:
        return '{0}E'.format(n)
    if rec_element.dtype == numpy.float64:
        return '{0}D'.format(n)
    
    # If it makes it here, assume its a string
    l = int(rec_element.dtype.str[rec_element.dtype.str.find('U')+1:])
#    return '{0}A'.format(l) if n==1 else '{0}A{1}'.format(l*n,l)
    return '{0}A'.format(l*n)


def rec_to_fits_col_dim(rec_element):
    """
    Return the string representation of the dimensions for the fits
    table column based on the provided record array element.

    The shape is inverted because the first element is supposed to be
    the most rapidly varying; i.e. the shape is supposed to be written
    as row-major, as opposed to the native column-major order in python.
    """
    return None if len(rec_element[0].shape) <= 1 else str(rec_element[0].shape[::-1])


def channel_dictionary(hdu, ext, prefix='C'):
    """
    Construct a dictionary of the channels in a MAPS file.
    """
    channel_dict = {}
    for k, v in hdu[ext].header.items():
        if k[:len(prefix)] == prefix:
            try:
                i = int(k[len(prefix):])-1
            except ValueError:
                continue
            channel_dict[v] = i
    return channel_dict


def channel_units(hdu, ext, prefix='U'):
    """
    Construct an array with the channel units.
    """
    cu = {}
    for k, v in hdu[ext].header.items():
        if k[:len(prefix)] == prefix:
            try:
                i = int(k[len(prefix):])-1
            except ValueError:
                continue
            cu[i] = v.strip()
    channel_units = numpy.empty(max(cu.keys())+1, dtype=object)
    for k, v in cu.items():
        channel_units[k] = v
    return channel_units.astype(str)


def compress_file(ifile, clobber=False):
    """
    Compress a file using gzip.  The output file has the same name as
    the input file with '.gz' appended.

    Any existing file will be overwritten if clobber is true.

    An error is raised if the input file name already has '.gz' appended
    to the end.
    """
    if ifile.split('.')[-1] == 'gz':
        raise ValueError('File appears to already have been compressed! {0}'.format(ifile))

    ofile = '{0}.gz'.format(ifile)
    if os.path.isfile(ofile) and not clobber:
        raise FileExistsError('File already exists: {0}.\nTo overwrite, set clobber=True.'.format(
                                ofile))

    with open(ifile, 'rb') as f_in:
        with gzip.open(ofile, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)


def create_symlink(ofile, symlink_dir, relative_symlink=True, clobber=False, loggers=None,
                   quiet=False):
    """
    Create a symlink to the input file in the provided directory.  If
    relative_symlink is True (default), the path to the file is relative
    to the directory with the symlink.
    """
    # Check if the file already exists
    olink_dest = os.path.join(symlink_dir, ofile.split('/')[-1])
    if os.path.isfile(olink_dest) or os.path.islink(olink_dest):
        if clobber:
            os.remove(olink_dest)
        else:
            return

    # Make sure the symlink directory exists
    if not os.path.isdir(symlink_dir):
        os.makedirs(symlink_dir)

    # Set the relative path for the symlink, if requested
    olink_src = os.path.relpath(ofile, start=os.path.dirname(olink_dest)) \
                    if relative_symlink else ofile
    if not quiet:
        log_output(loggers, 1, logging.INFO, 'Creating symlink: {0}'.format(olink_dest))

    # Create the symlink
    os.symlink(olink_src, olink_dest)

