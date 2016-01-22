#!/usr/bin/env python3

from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from __future__ import unicode_literals

import sys
if sys.version > '3':
    long = int

import numpy
from os import remove
import os.path
from time import clock

from argparse import ArgumentParser
from astropy.io import fits
from astropy.wcs import WCS
from mangadap.dapfile import dapfile
from mangadap.util.exception_tools import print_frame

#-----------------------------------------------------------------------------

# TODO: Temporary.  Should fix for MPL-5
def clean_primary_header(hdr):
    # Remove some keys that are incorrect for DAP data
    hdr.remove('BUNIT')
    hdr.remove('MASKNAME')
    # Remove WCS information
    for i in range(3):
        hdr.remove('CTYPE{0}'.format(i+1))
        hdr.remove('CUNIT{0}'.format(i+1))
        hdr.remove('CRVAL{0}'.format(i+1))
        hdr.remove('CRPIX{0}'.format(i+1))
        hdr.remove('CD{0}_{0}'.format(i+1))
    # Remove HDUCLASS keys
    hdr.remove('HDUCLASS')
    hdr.remove('HDUCLAS1')
    hdr.remove('HDUCLAS2')
    hdr.remove('ERRDATA')
    hdr.remove('QUALDATA')


def clean_data_header(hdr, nchannels=1):

#    # Remove the 3rd axis
#    hdr.remove('CTYPE3')
#    hdr.remove('CUNIT3')
#    hdr.remove('CRPIX3')
#    hdr.remove('CRVAL3')
#    hdr.remove('CD3_3')

    out_hdr = hdr.copy()
    
    # Change header keywords to the default values for the third axis
    if nchannels > 1:
        out_hdr['NAXIS'] = 3
        out_hdr.remove('CTYPE3')
        out_hdr.remove('CUNIT3')
        out_hdr['CTYPE3'] = ' '
        out_hdr['CUNIT3'] = ' '
        out_hdr['CRPIX3'] = 1
        out_hdr['CRVAL3'] = 1.
        out_hdr['CD3_3']  = 1.
    else:
        out_hdr['NAXIS'] = 2
        out_hdr.remove('CTYPE3')
        out_hdr.remove('CUNIT3')
        out_hdr.remove('CRPIX3')
        out_hdr.remove('CRVAL3')
        out_hdr.remove('CD3_3')

#    print(hdr['CTYPE3'])
#    print(hdr['CUNIT3'])
#    print(hdr['CRPIX3'])
#    print(hdr['CRVAL3'])
#    print(hdr['CD3_3'])

    # Remove everything but the WCS information
    w = WCS(header=out_hdr)
    out_hdr = w.to_header().copy()

    # Fix the DATE-OBS keyword:
    out_hdr.comments['DATE-OBS'] = 'Date of median exposure'
    out_hdr.comments['MJD-OBS'] = '[d] MJD for DATE-OBS'

    # Add back in the BSCALE and BZERO values
    out_hdr['BSCALE'] = 1.
    out_hdr['BZERO'] = 0.

    # Add back in the default CTYPE3, CUNIT3
    if nchannels > 1:
        out_hdr['CTYPE3'] = (' ', 'Undefined type')
        out_hdr['CUNIT3'] = (' ', 'Undefined units')

    return out_hdr


def finalize_data_header(data_hdr, ext, bunit=None, err=True, qual=True, nchannels=None):
    if bunit is not None:
        data_hdr['BUNIT'] = (bunit, 'Unit of pixel value')
    # Add the common HDUCLASS keys
    data_hdr['HDUCLASS'] = ('SDSS', 'SDSS format class')
    if nchannels is not None and nchannels == 1:
        data_hdr['HDUCLAS1'] = 'IMAGE'
    else:
        data_hdr['HDUCLAS1'] = 'CUBE'
    data_hdr['HDUCLAS2'] = 'DATA'
#    if err:
#        data_hdr['ERRDATA'] = (ext+'_ERR', 'Associated error extension')
    if err:
        data_hdr['ERRDATA'] = (ext+'_IVAR', 'Associated inv. variance extension')
    if qual:
        data_hdr['QUALDATA'] = (ext+'_MASK', 'Associated quality extension')


def finalize_ivar_header(ivar_hdr, ext, bunit=None, nchannels=None):
    if bunit is not None:
        ivar_hdr['BUNIT'] = (bunit, 'Unit of pixel value')
    # Add the common HDUCLASS keys
    ivar_hdr['HDUCLASS'] = ('SDSS', 'SDSS format class')
    if nchannels is not None and nchannels == 1:
        ivar_hdr['HDUCLAS1'] = 'IMAGE'
    else:
        ivar_hdr['HDUCLAS1'] = 'CUBE'
    ivar_hdr['HDUCLAS2'] = 'ERROR'
#    rmse_hdr['HDUCLAS3'] = ('RMSE', 'Value is root-mean-square error')
    ivar_hdr['HDUCLAS3'] = ('INVMSE', 'Value is inverse mean-square error')
    ivar_hdr['SCIDATA'] = (ext, 'Associated data extension')
    ivar_hdr['QUALDATA'] = (ext+'_MASK', 'Associated quality extension')


def finalize_mask_header(mask_hdr, ext, nchannels=None):
    # Add the common HDUCLASS keys
    mask_hdr['HDUCLASS'] = ('SDSS', 'SDSS format class')
    if nchannels is not None and nchannels == 1:
        mask_hdr['HDUCLAS1'] = 'IMAGE'
    else:
        mask_hdr['HDUCLAS1'] = 'CUBE'
    mask_hdr['HDUCLAS2'] = 'QUALITY'
    mask_hdr['HDUCLAS3'] = ('MASKZERO', 'Zero values are good')
    mask_hdr['SCIDATA'] = (ext, 'Associated data extension')
#    mask_hdr['ERRDATA'] = (ext+'_ERR', 'Associated error extension')
    mask_hdr['ERRDATA'] = (ext+'_IVAR', 'Associated inv. variance extension')


def finalize_bin_header(bin_hdr):
    bin_hdr['HDUCLASS'] = ('SDSS', 'SDSS format class')
#    bin_hdr['HDUCLAS1'] = 'CUBE'
    bin_hdr['HDUCLAS1'] = 'IMAGE'
    bin_hdr['HDUCLAS2'] = 'DATA'


def primary_hdu(dapf):
    hdr = dapf.hdu['PRIMARY'].header.copy()
    clean_primary_header(hdr)
    # Add Authors
    hdr['AUTHOR'] = 'K Westfall & B Andrews <kyle.westfall@port.co.uk, andrewsb@pitt.edu>'
    return fits.PrimaryHDU(header=hdr)


def place_values_in_map(x, y, z, xs, dx, ys, dy, channel):
    """
    Given a grid with the bottom corner of the first pixel at (xs, ys)
    and pixels of size dxdy, place the z values at locations x, y in
    channel, where channel is expected to have the correct size upon
    input.  The x and y coordinates are expected to be sampled at the
    same rate as the grid!
    """
    
    nx = channel.shape[0]
    ny = channel.shape[1]

    i = numpy.floor( (x-xs)/dx )
    j = numpy.floor( (y-ys)/dy )

    for ii,jj,zz in zip(i,j,z):
        if ii < 0 or ii >= nx:
            print('Could not place value in map!')
            continue
        if jj < 0 or jj >= ny:
            print('Could not place value in map!')
            continue
        channel[ii,jj] = zz


def place_values_in_bins(binid, zz, default_value=0., dtype=numpy.float64):
    """
    Given a set of bin values zz and a set of bin assignments, return a
    vector with the bin values assigned to the appropriate element of
    the prebinned data.

    - Allow for an input default (non-zero) value?
    """
    z = numpy.full(binid.size, default_value, dtype=dtype)
    valid = numpy.invert(binid < 0)
    z[valid] = zz[binid[valid]]
    return z


def number_of_channels(dapf, ext, data_col, vec_i=None):
    
    # Number of channels with different species
    dim = dapf.hdu[ext].data[data_col].shape
    ndim = len(dim)
    if ndim == 1:
        nchannels = 1
    elif ndim == 2:
        nchannels = dim[-1] if vec_i is None else 1
    elif ndim == 3 and vec_i is None:
        raise ValueError('Select column is three dimensional.  Must specify vector element.')
    else:
        nchannels = dim[-1]

    return dim, nchannels


def binid_hdu(dapf, output_ext, pixelscale=0.5, mpl=None):#base_hdr, 

    # Number of on-sky X,Y pixels in the output map
    npix = int(numpy.sqrt(dapf.hdu['DRPS'].header['NAXIS2']))

    # Add a third axis to the header keywords if necessary
    bin_hdr = clean_data_header(dapf.hdu['PRIMARY'].header.copy(), nchannels=1)

    # Set the pixel coordinates
    dx = pixelscale
    xs = -dx*npix/2-dx/2.
    ys = xs
    dy = dx

    # Save the x and y coordinates
    X = dapf.hdu['DRPS'].data['XPOS']
    Y = dapf.hdu['DRPS'].data['YPOS']
        
    # Have to adjust the x and y values given the error in DAP v1_0_0
    if mpl == 3:
        X += 0.5
        Y -= 0.5

    # Initialize the output data cube and an individual channel
    data_type = dapf.hdu['DRPS'].data['BINID'].dtype
    output_data = numpy.zeros((npix, npix), dtype=data_type)

    z = dapf.hdu['DRPS'].data['BINID'].ravel()
    place_values_in_map(-X, Y, z, xs, dx, ys, dy, output_data)

    # Finalize the headers
    finalize_bin_header(bin_hdr)

    # Return the hdu as a list
    return [ fits.ImageHDU(data=output_data.transpose(), header=bin_hdr, name=output_ext) ]


def single_quantity_hdu(dapf, ext, col, output_ext, pixelscale=0.5, name_ext=None, #base_hdr, 
                        name_col=None, mpl=None, bunit=None, convert_to_per_spaxel=False,
                        channels=None):

    if (name_ext is not None and name_col is None) or (name_ext is None and name_col is not None):
        raise ValueError('If reading names from file, must provide both extension and column.')

    # Number of on-sky X,Y pixels in the output map
    npix = int(numpy.sqrt(dapf.hdu['DRPS'].header['NAXIS2']))

    # Number of channels with different species
    dim, nchannels = number_of_channels(dapf, ext, col)
    ndim = len(dim)

    # If a channel selection is not provided, output them all
    if channels is None:
        channels = numpy.arange(nchannels)

    # Check that selected channels are available
    if any(numpy.array(channels) >= nchannels):
        raise valueerror('selected channels unavailable!')

    # Get the number of channels to include in the output
    out_nchannels = len(channels)

    # Copy the base header
    data_hdr = clean_data_header(dapf.hdu['PRIMARY'].header.copy(), nchannels=out_nchannels)

    # Initialize the output data cube and an individual channel
    data_type = dapf.hdu[ext].data[col].dtype
    output_data = numpy.zeros((out_nchannels, npix, npix), dtype=data_type)

    channel = numpy.zeros((npix, npix), dtype=data_type)

    # Get the number of digits to use for the channel desciptions
    ndig = int(numpy.log10(out_nchannels))+1

    # Set the pixel coordinates
    dx = pixelscale
    xs = -dx*npix/2-dx/2.
    ys = xs
    dy = dx

    # Convert units from per arcsec^2 to per spaxel if requested
    unit_conversion = dx*dx if convert_to_per_spaxel else 1.0

    # Save the BINID, and x and y coordinates
    binid = dapf.hdu['DRPS'].data['BINID']
    X = dapf.hdu['DRPS'].data['XPOS']
    Y = dapf.hdu['DRPS'].data['YPOS']
        
    # Have to adjust the x and y values given the error in DAP v1_0_0
    if mpl == 3:
        X += 0.5
        Y -= 0.5

    # Get the data
    data = dapf.hdu[ext].data[col]

    if data.shape[-1] != nchannels:
        print(data.shape)
        raise ValueError('Unable to extract data with proper size.')

    # Build the maps for each channel
    j = 0   # Output channel iterable
    for i in range(nchannels):
        if i not in channels:
            continue

        z_data = data[:,i].ravel()

        # Add the values to the maps
        place_values_in_map(-X, Y, place_values_in_bins(binid, z_data, dtype=data_type), xs, dx,
                            ys, dy, channel)
        output_data[j,:,:] = channel[:,:]*unit_conversion

        # Set the header keyword for this channel
        if name_ext is not None and name_col is not None:
            data_hdr['C'+'{0}'.format(j+1).zfill(ndig)] = (dapf.hdu[name_ext].data[name_col][i],
                                                        'Species in cube channel {0}'.format(j+1))
        else:
            data_hdr['C'+'{0}'.format(j+1).zfill(ndig)] = (ext+':'+data_col,
                                                           'DAP-full extension:column')
            if vec_i is not None:
                data_hdr['C'+'{0}'.format(j+1).zfill(ndig)] += ':' + str(vec_i)

        # Increment output channel
        j += 1

    # Finalize the headers
    finalize_data_header(data_hdr, output_ext, bunit=bunit, err=False, qual=False,
                         nchannels=out_nchannels)

    if out_nchannels > 1:
        output_data = numpy.transpose(output_data, (0,2,1))
    else:
        output_data = numpy.transpose(output_data[0,:,:])

    # Return the list of hdu instances
    return [ fits.ImageHDU(data=output_data, header=data_hdr, name=output_ext) ]


def measured_quantity_hdus(dapf, ext, data_col, err_col, output_ext, mask_col=None, #base_hdr, 
                           vec_i=None, pixelscale=0.5, name_ext=None, name_col=None,
                           unit_ext=None, unit_col=None, channel_units=None, mpl=None, bunit=None,
                           convert_to_per_spaxel=False, offset=None, channels=None):

    # Perform some input checks
    if channel_units is not None and unit_ext is not None:
        raise ValueError('Cannot provide units both directly and via an extension.')

    if (unit_ext is not None and unit_col is None) or (unit_ext is None and unit_col is not None):
        raise ValueError('If reading units from file, must provide both extension and column.')

    if (name_ext is not None and name_col is None) or (name_ext is None and name_col is not None):
        raise ValueError('If reading names from file, must provide both extension and column.')

    # Number of on-sky X,Y pixels in the output map
    npix = int(numpy.sqrt(dapf.hdu['DRPS'].header['NAXIS2']))

    # Number of channels with different species
    dim, nchannels = number_of_channels(dapf, ext, data_col, vec_i=vec_i)
    ndim = len(dim)

    # If a channel selection is not provided, output them all
    if channels is None:
        channels = numpy.arange(nchannels)

    # Check that selected channels are available
    if any(numpy.array(channels) >= nchannels):
        raise ValueError('Selected channels unavailable!')

    # Get the number of channels to include in the output
    out_nchannels = len(channels)

    # Check that the appropriate number of units was provided
    if channel_units is not None and len(channel_units) != out_nchannels:
        raise ValueError('Provided channel units of does not match number of channels.')

    # Copy the base header
    data_hdr = clean_data_header(dapf.hdu['PRIMARY'].header.copy(), nchannels=out_nchannels)

    # Initialize the output data cube and an individual channel
    data_type = dapf.hdu[ext].data[data_col].dtype
    output_data = numpy.zeros((out_nchannels, npix, npix), dtype=data_type)
#    output_rmse = numpy.zeros((out_nchannels, npix, npix), dtype=data_type)
    output_ivar = numpy.zeros((out_nchannels, npix, npix), dtype=data_type)
    output_mask = numpy.zeros((out_nchannels, npix, npix), dtype=numpy.uint8)

    channel_flt = numpy.zeros((npix, npix), dtype=data_type)
    channel_int = numpy.zeros((npix, npix), dtype=numpy.uint8)

    # Get the number of digits to use for the channel desciptions
    ndig = int(numpy.log10(out_nchannels))+1

    # Set the pixel coordinates
    dx = pixelscale
    xs = -dx*npix/2-dx/2.
    ys = xs
    dy = dx

    # Convert units from per arcsec^2 to per spaxel if requested
    unit_conversion = dx*dx if convert_to_per_spaxel else 1.0

    # Save the BINID, and x and y coordinates
    binid = dapf.hdu['DRPS'].data['BINID']
    X = dapf.hdu['DRPS'].data['XPOS']
    Y = dapf.hdu['DRPS'].data['YPOS']
        
    # Have to adjust the x and y values given the error in DAP v1_0_0
    if mpl == 3:
        X += 0.5
        Y -= 0.5

    # Get the data
    if vec_i is None:
        data = dapf.hdu[ext].data[data_col]
        rmse = dapf.hdu[ext].data[err_col]
        mask = numpy.full(data.shape, 0, dtype=numpy.uint8) if mask_col is None else \
                    dapf.hdu[ext].data[mask_col].astype(numpy.uint8)
    else:
        if ndim == 2:   # nchannel should be 1
            data = numpy.array([dapf.hdu[ext].data[data_col][:,vec_i].ravel()]).transpose()
            rmse = numpy.array([dapf.hdu[ext].data[err_col][:,vec_i].ravel()]).transpose()
        if ndim == 3:
            data = dapf.hdu[ext].data[data_col][:,vec_i,:].reshape(dim[0],dim[2])
            rmse = dapf.hdu[ext].data[err_col][:,vec_i,:].reshape(dim[0],dim[2])
        # NOTE: vec_i does not apply to the mask column
        mask = numpy.full(data.shape, 0, dtype=numpy.uint8) if mask_col is None else \
                    dapf.hdu[ext].data[mask_col].astype(numpy.uint8)

#    print(data.shape)
#    print(numpy.min(data[numpy.where(data > 0)]), numpy.max(data))
#    print(numpy.min(data[data[:,9]>0,9]), numpy.max(data[:,9]))
#    print(rmse.shape)
#    print(mask.shape)

    if data.shape[-1] != nchannels or rmse.shape[-1] != nchannels or mask.shape[-1] != nchannels:
        print(data.shape)
        print(rmse.shape)
        print(mask.shape)
        raise ValueError('Unable to extract data with proper size.')

    if offset is not None:
        data -= offset

    # Build the maps for each channel
    j = 0   # Output channel iterable
    for i in range(nchannels):
        if i not in channels:
            continue

        z_data = data[:,i].ravel()
#        z_rmse = rmse[:,i].ravel()
        z_ivar = rmse[:,i].ravel()
        z_mask = mask[:,i].ravel()

#        print(binid.shape)
#        print(z_data.shape)
#        print(numpy.min(z_data[z_data > 0]), numpy.max(z_data))

        # Calculate inverse variance and mask anything with error <= 0
        indx = z_ivar > 0
        z_ivar[indx] = 1.0/numpy.square(z_ivar[indx])
        z_ivar[numpy.invert(indx)] = 0.
        z_mask[numpy.invert(indx)] = numpy.uint8(1)

        # Add the values to the maps
        place_values_in_map(-X, Y, place_values_in_bins(binid, z_data, dtype=data_type), xs, dx,
                            ys, dy, channel_flt)
        output_data[j,:,:] = channel_flt[:,:]*unit_conversion
#        print(numpy.min(channel_flt[numpy.where(channel_flt > 0)]), numpy.max(channel_flt))

#        place_values_in_map(-X, Y, place_values_in_bins(binid, z_rmse, dtype=data_type), xs, dx,
#                            ys, dy, channel_flt)
#        output_rmse[j,:,:] = channel_flt[:,:]*unit_conversion

        place_values_in_map(-X, Y, place_values_in_bins(binid, z_ivar, dtype=data_type), xs, dx,
                            ys, dy, channel_flt)
        # Fixed bug: 9 Dec 2015:
        # output_ivar[j,:,:] = channel_flt[:,:]*unit_conversion
        output_ivar[j,:,:] = channel_flt[:,:]/numpy.square(unit_conversion)

        place_values_in_map(-X, Y,
                            place_values_in_bins(binid, z_mask, default_value=1, dtype=numpy.uint8),
                            xs, dx, ys, dy, channel_int)
        output_mask[j,:,:] = channel_int[:,:]

        # Set the header keyword for this channel
        if name_ext is not None and name_col is not None:
            data_hdr['C'+'{0}'.format(j+1).zfill(ndig)] = (dapf.hdu[name_ext].data[name_col][i],
                                                        'Species in cube channel {0}'.format(j+1))
        else:
            data_hdr['C'+'{0}'.format(j+1).zfill(ndig)] = (ext+':'+data_col,
                                                           'DAP-full extension:column')
            if vec_i is not None:
                data_hdr['C'+'{0}'.format(j+1).zfill(ndig)] += ':' + str(vec_i)

        # Set the header keyword for this channel
        if unit_ext is not None and unit_col is not None:
            data_hdr['U'+'{0}'.format(j+1).zfill(ndig)] = (dapf.hdu[unit_ext].data[unit_col][i],
                                                'Unit of species in cube channel {0}'.format(j+1))
        if channel_units is not None:
            data_hdr['U'+'{0}'.format(j+1).zfill(ndig)] = (channel_units[j],
                                                'Unit of species in cube channel {0}'.format(j+1))

        # Increment output channel
        j += 1

    # Finalize the headers
    bunit_ivar = None if bunit is None else '('+bunit+')^{-2}'
    ivar_hdr = data_hdr.copy()
#    rmse_hdr = data_hdr.copy()
    mask_hdr = data_hdr.copy()
    finalize_data_header(data_hdr, output_ext, bunit=bunit, nchannels=out_nchannels)
    finalize_ivar_header(ivar_hdr, output_ext, bunit=bunit_ivar, nchannels=out_nchannels)
#    finalize_rmse_header(rmse_hdr, output_ext, bunit=bunit)
    finalize_mask_header(mask_hdr, output_ext, nchannels=out_nchannels)

    if out_nchannels > 1:
        output_data = numpy.transpose(output_data, (0,2,1))
#        output_rmse = numpy.transpose(output_rmse, (0,2,1))
        output_ivar = numpy.transpose(output_ivar, (0,2,1))
        output_mask = numpy.transpose(output_mask, (0,2,1))
    else:
        output_data = numpy.transpose(output_data[0,:,:])
#        output_rmse = numpy.transpose(output_rmse[0,:,:])
        output_ivar = numpy.transpose(output_ivar[0,:,:])
        output_mask = numpy.transpose(output_mask[0,:,:])

    # Return the list of hdu instances
    return [ fits.ImageHDU(data=output_data, header=data_hdr, name=output_ext), \
#             fits.ImageHDU(data=output_rmse, header=rmse_hdr, name=output_ext+'_ERR'), \
             fits.ImageHDU(data=output_ivar, header=ivar_hdr, name=output_ext+'_IVAR'), \
             fits.ImageHDU(data=output_mask, header=mask_hdr, name=output_ext+'_MASK', uint=True) ]


def create_dap_map_file(plate, ifudesign, bintype, output_file, mpl=3, plan_index=None,
                         directory_path=None):

    if mpl == 4 and plan_index is None:
        raise Exception('Must provide plan index for MPL-4 data!')

    if bintype != 'STON' and bintype != 'NONE' and bintype != 'RADIAL':
        raise KeyError('Binning type must be STON, NONE, or RADIAL!')

    if mpl == 3 and bintype == 'STON':
        plan_index = 1
    if mpl == 3 and bintype == 'NONE':
        plan_index = 2
    
    print('Attempting to read DAP output file ...')
    dapf = dapfile(plate, ifudesign, 'CUBE', bintype, plan_index, directory_path=directory_path,
                   checksum=False)
    print('Found: {0}'.format(dapf.file_path()))

    #-------------------------------------------------------------------
    # Build the primary header
    hdulist = [ primary_hdu(dapf) ]

    #-------------------------------------------------------------------
    # Get a common set of header keywords for all channels
#    base_hdr = clean_data_header(dapf.hdu['PRIMARY'].header.copy())
#    print(base_hdr)
#    exit()

    #-------------------------------------------------------------------
    # Build a set of three extensions for each fitted quantity, and add
    # the list of three extensions to the total
    n_quantities = 0

    # Emission-line Gaussian flux
    ext = 'EMLINE_GFLUX'
    print('Adding extension: {0:>10s}'.format(ext))#, end='\r')
    hdulist += measured_quantity_hdus(dapf, 'ELOFIT', 'FLUX_EW', 'FLUXERR_EW', ext, #base_hdr,
                                      mask_col='ELOMIT_EW', name_ext='ELOPAR', name_col='ELNAME',
                                      mpl=mpl, bunit='1E-17 erg/s/cm^2/spaxel',
                                      convert_to_per_spaxel=(mpl==4),
                                      channels=[2,3,4,5,6,7,8,9,10,11,12])
    n_quantities += 1

    # Emission-line Gaussian velocity
    ext = 'EMLINE_GVEL'
    print('Adding extension: {0:>10s}'.format(ext))#, end='\r')
    hdulist += measured_quantity_hdus(dapf, 'ELOFIT', 'IKIN_EW', 'IKINERR_EW', ext, #base_hdr,
                                      mask_col='ELOMIT_EW', vec_i=0, name_ext='ELOPAR',
                                      name_col='ELNAME', mpl=mpl, bunit='km/s',
                                      offset=dapf.guess_cz(), channels=[2,3,4,5,6,7,8,9,10,11,12])
    hdulist[-3].header['VSYS_NSA'] = (dapf.guess_cz(), 'Systemic velocity (km/s)')
    n_quantities += 1

    # Emission-line Gaussian velocity dispersion
    ext = 'EMLINE_GSIGMA'
    print('Adding extension: {0:>10s}'.format(ext))#, end='\r')
    hdulist += measured_quantity_hdus(dapf, 'ELOFIT', 'IKIN_EW', 'IKINERR_EW', ext, #base_hdr,
                                      mask_col='ELOMIT_EW', vec_i=1, name_ext='ELOPAR',
                                      name_col='ELNAME', mpl=mpl, bunit='km/s',
                                      channels=[2,3,4,5,6,7,8,9,10,11,12])
    n_quantities += 1

    ext = 'EMLINE_INSTSIGMA'
    print('Adding extension: {0:>10s}'.format(ext))#, end='\r')
    hdulist += single_quantity_hdu(dapf, 'ELOFIT', 'SINST_EW', ext, #base_hdr,
                                   name_ext='ELOPAR', name_col='ELNAME', mpl=4, bunit='km/s',
                                   channels=[2,3,4,5,6,7,8,9,10,11,12])
    n_quantities += 1

    # Emission-line Gaussian EW
    ext = 'EMLINE_EW'
    print('Adding extension: {0:>10s}'.format(ext))#, end='\r')
    hdulist += measured_quantity_hdus(dapf, 'ELOFIT', 'EW_EW', 'EWERR_EW', ext, #base_hdr,
                                      mask_col='ELOMIT_EW', name_ext='ELOPAR', name_col='ELNAME',
                                      mpl=mpl, bunit='Ang', channels=[2,3,4,5,6,7,8,9,10,11,12])
    n_quantities += 1

    # Emission-line non-parametric line flux
    ext = 'EMLINE_SFLUX'
    print('Adding extension: {0:>10s}'.format(ext))#, end='\r')
    hdulist += measured_quantity_hdus(dapf, 'ELMMNT', 'FLUX', 'FLUXERR', ext, #base_hdr,
                                      mask_col='OMIT', name_ext='ELBAND', name_col='ELNAME',
                                      mpl=mpl, bunit='1E-17 erg/s/cm^2/spaxel',
                                      convert_to_per_spaxel=(mpl==4))

    hdulist[-3].header['C01'] = ('OIId---3728', 'Species in cube channel 1')
    hdulist[-2].header['C01'] = ('OIId---3728', 'Species in cube channel 1')
    hdulist[-1].header['C01'] = ('OIId---3728', 'Species in cube channel 1')
    n_quantities += 1

    # Stellar velocity
    ext = 'STELLAR_VEL'
    print('Adding extension: {0:>10s}'.format(ext))#, end='\r')
    hdulist += measured_quantity_hdus(dapf, 'STFIT', 'KIN', 'KINERR', ext, #base_hdr,
                                      vec_i=0, mpl=mpl, bunit='km/s', offset=dapf.guess_cz())
    hdulist[-3].header['VSYS_NSA'] = (dapf.guess_cz(), 'Systemic velocity (km/s)')
    n_quantities += 1

    # Stellar velocity dispersion
    ext = 'STELLAR_SIGMA'
    print('Adding extension: {0:>10s}'.format(ext))#, end='\r')
    hdulist += measured_quantity_hdus(dapf, 'STFIT', 'KIN', 'KINERR', ext, #base_hdr,
                                      vec_i=1, mpl=mpl, bunit='km/s')
    n_quantities += 1

    # Selected spectral indices are:
    #   D4000, CaII0p39, HDeltaA, CN1, CN2, Ca4227, HGammaA, Fe4668, Hb,
    #   Mgb, Fe5270, Fe5335, Fe5406, NaD, TiO1, TiO2, NaI0p82,
    #   CaII0p86A, CaII0p86B, CaII0p86C, MgI0p88, TiO0p89, FeH0p99
    ext = 'SPECINDEX'
    sindx_channels = [0,1,2,4,5,6,8,13,14,18,19,20,21,24,25,26,27,28,29,30,31,32,33]
    sindx_units = [ 'Angstrom', 'Angstrom', 'Angstrom', 'mag', 'mag', 'Angstrom', 'Angstrom',
                    'Angstrom', 'Angstrom', 'Angstrom', 'Angstrom', 'Angstrom', 'Angstrom',
                    'Angstrom', 'mag', 'mag', 'Angstrom', 'Angstrom', 'Angstrom', 'Angstrom',
                    'Angstrom', 'Angstrom', 'Angstrom' ]
    print('Adding extension: {0:>10s}'.format(ext))#, end='\r')
    hdulist += measured_quantity_hdus(dapf, 'SINDX', 'INDX', 'INDXERR', ext, #base_hdr,
                                      mask_col='SIOMIT', name_ext='SIPAR', name_col='SINAME',
                                      channel_units=sindx_units, mpl=mpl, channels=sindx_channels)
                                      #, bunit='Ang or mag', convert_to_per_spaxel=(mpl==4))
    n_quantities += 1

    # Bin data
    ext = 'BINID'
    print('Adding extension: {0:>10s}'.format(ext))#, end='\r')
    hdulist += binid_hdu(dapf, ext, mpl=mpl)#base_hdr, 
    n_quantities += 1

    print('Finished adding extensions for all {0} quantities.            '.format(n_quantities))

    # Write the datacube
    print('Writing: {0}'.format(output_file))
    fits.HDUList( hdulist ).writeto(output_file, clobber=os.path.isfile(output_file),
                                    checksum=True)


#-----------------------------------------------------------------------------

if __name__ == '__main__':
    t = clock()

    parser = ArgumentParser()

    parser.add_argument('plate', type=int, help='plate ID to process')
    parser.add_argument('ifudesign', type=int, help='IFU design to process')
    parser.add_argument('mpl', type=int, help='MPL number (3 or 4)')
    parser.add_argument('output_file', type=str, help='Name for output file')
    parser.add_argument('-b', '--bintype', type=str,
                        help='Binning type to process: NONE(def), STON, or RADIAL', default='NONE')
    parser.add_argument('-i', '--plan_index', type=int,
                        help='Index of plan used by DAP (used in the DAP output file name)',
                        default=None)
    parser.add_argument('-d', '--directory_path', type=str, help='Path to DAP output file',
                        default=None)

    arg = parser.parse_args()

    if os.path.isfile(arg.output_file):
        print('WARNING: Overwriting existing file {0}!'.format(arg.output_file))

    create_dap_map_file(arg.plate, arg.ifudesign, arg.bintype, arg.output_file, mpl=arg.mpl,
                         plan_index=arg.plan_index, directory_path=arg.directory_path)
 
    print('Elapsed time: {0} seconds'.format(clock() - t))



