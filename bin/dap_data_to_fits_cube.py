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
from mangadap.dapfile import dapfile
from mangadap.util.exception_tools import print_frame

#-----------------------------------------------------------------------------

def place_values_in_map(x, y, z, xs, dx, ys, dy, plane):
    """
    Given a grid with the bottom corner of the first pixel at (xs, ys)
    and pixels of size dxdy, place the z values at locations x, y in
    plane, where plane is expected to have the correct size upon input.
    The x and y coordinates are expected to be sampled at the same rate
    as the grid!
    """
    
    nx = plane.shape[0]
    ny = plane.shape[1]

    i = numpy.floor( (x-xs)/dx )
    j = numpy.floor( (y-ys)/dy )

    for ii,jj,zz in zip(i,j,z):
        if ii < 0 or ii >= nx:
            print('Could not place value in map!')
            continue
        if jj < 0 or jj >= ny:
            print('Could not place value in map!')
            continue
        plane[ii,jj] = zz


def place_values_in_bins(binid, zz):
    """
    Given a set of bin values zz and a set of bin assignments, return a
    vector with the bin values assigned to the appropriate element of
    the prebinned data.

    - Allow for an input default (non-zero) value?
    """
    z = numpy.zeros(binid.size, dtype=numpy.float64)
    valid = numpy.invert(binid < 0)
    z[valid] = zz[binid[valid]]
    return z
#    nspax = binid.size
#    z = numpy.zeros(nspax, dtype=numpy.float64)
#    for i in range(nspax):
#        if binid[i] < 0:
#            continue
#        z[i] = zz[binid[i]]
#    return z


def create_dap_fits_cube(plate, ifudesign, bintype, output_file, mpl=3, plan_index=None,
                         directory_path=None):

    if mpl == 4 and plan_index is None:
        raise Exception('Must provide plan index for MPL-4 data!')

    if bintype != 'STON' and bintype != 'NONE':
        raise KeyError('Binning type must be STON or NONE!')

    if mpl == 3 and bintype == 'STON':
        plan_index = 1
    if mpl == 3 and bintype == 'NONE':
        plan_index = 2
    
    print('Attempting to read DAP output file ...')
    dapf = dapfile(plate, ifudesign, 'CUBE', bintype, plan_index, directory_path=directory_path)
    print('Found: {0}'.format(dapf.file_path()))

    #-------------------------------------------------------------------
    # Initialize the header
    newhdr = dapf.hdu['PRIMARY'].header

    # Remove the wavelength axis and replace it with just the pixel coordinates
    newhdr.remove('CTYPE3')
    newhdr.remove('CUNIT3')
    newhdr['CRPIX3'] = 1
    newhdr['CRVAL3'] = 1.
    newhdr['CD3_3'] = 1.

    #-------------------------------------------------------------------
    # Determine the number of maps to produce
    nstkin = dapf.hdu['STFIT'].data['KIN'].shape[1]
    nelowgt = 1 if mpl == 3 else 4
    nelokin = int(dapf.hdu['ELOFIT'].data['KIN_EW'].shape[1]/nelowgt)
    neml = dapf.hdu['ELOFIT'].data['ELOMIT_EW'].shape[1]
    if mpl == 4:
        nelmmnt = dapf.hdu['ELMMNT'].data['OMIT'].shape[1]
    else:
        nelmmnt = 0
    nsindx = dapf.hdu['SINDX'].data['INDX'].shape[1]

    #-------------------------------------------------------------------
    # List the columns to print
    drps_col = [ 'FGOODPIX', 'SIGNAL', 'NOISE', 'BINID', 'BINW' ]
    bins_col = [ 'BINSN' ]

    stfit_vec = [ 'KIN', 'KINERR' ]
    stfit_col = [ 'RCHI2' ]

    if mpl == 4:
        elmmnt_vec = [ 'OMIT', 'FLUX_RAW', 'FLUXERR_RAW', 'MOM1_RAW', 'MOM1ERR_RAW', 'MOM2_RAW',
                       'MOM2ERR_RAW', 'BLUF_RAW', 'BLUFERR_RAW',  'REDF_RAW', 'REDFERR_RAW',
                       'BFCEN', 'BCONT', 'BCONTERR', 'RFCEN', 'RCONT', 'RCONTERR', 'FLUX',
                       'FLUXERR', 'MOM1', 'MOM1ERR', 'MOM2', 'MOM2ERR' ]
    else:
        elmmnt_vec = []

    elofit_kin_vec = [ 'KIN', 'KINERR' ]
    elofit_col = [ 'NKIN' ]
    if mpl == 4:
        elofit_kin_vec += [ 'KINSTDE' ]

    elofit_ind_vec = [ 'ELOMIT', 'AMPL', 'AMPLERR', 'IKIN', 'IKINERR', 'SINST', 'FLUX', 'FLUXERR',
                        'EW', 'EWERR' ]
    elofit_ind_vec_nel = [ 1, 1, 1, nelokin, nelokin, 1, 1, 1, 1, 1 ]

    if mpl == 3:
        sindx_vec = [ 'SIOMIT', 'INDX', 'INDXERR', 'INDX_OTPL', 'INDX_BOTPL' ]
    if mpl == 4:
        sindx_vec = [ 'SIOMIT', 'BCONT', 'BCONTERR', 'RCONT', 'RCONTERR', 'INDX_RAW',
                      'INDXERR_RAW', 'INDXPERR_RAW', 'BCONT_OTPL', 'RCONT_OTPL', 'INDX_OTPL',
                      'BCONT_BOTPL', 'RCONT_BOTPL', 'INDX_BOTPL', 'INDX', 'INDXERR', 'INDXPERR' ]

    nmaps = len(drps_col) + len(bins_col) + nstkin*len(stfit_vec) + len(stfit_col) \
            + nelmmnt*len(elmmnt_vec) + 2*(nelokin*nelowgt*len(elofit_kin_vec) + len(elofit_col) \
            + neml*sum(elofit_ind_vec_nel)) + nsindx*len(sindx_vec)

#    nmaps = len(drps_col) + len(bins_col) + nstkin*len(stfit_vec) + len(stfit_col) \
#            + nelmmnt*len(elmmnt_vec) + nelokin*nelowgt*len(elofit_kin_vec) + len(elofit_col) \
#            + neml*sum(elofit_ind_vec_nel) + nsindx*len(sindx_vec)

    ndig = int(numpy.log10(nmaps))+1

    # DRP CUBE file is always assumed to be square in the number of
    # on-sky pixels!
    npix = int(numpy.sqrt(dapf.hdu['DRPS'].header['NAXIS2']))

    print('Number of stellar kinematic moments:                 {0}'.format(nstkin))
    print('Number of emission-line kinematic moments:           {0}'.format(nelokin))
    print('Number of emission-line kinematic weighting schemes: {0}'.format(nelowgt))
    print('Number of fitted emission lines:                     {0}'.format(neml))
    print('Number of emission moment bands:                     {0}'.format(nelmmnt))
    print('Number of spectral indices:                          {0}'.format(nsindx))
    print('Total number of maps to produce:                     {0}'.format(nmaps))
    print('Size of each map:                                    {0} pix x {0} pix'.format(npix))

    output_data = numpy.zeros((nmaps, npix, npix), dtype=numpy.float64)
    plane = numpy.zeros((npix, npix), dtype=numpy.float64)

    xs = -npix/4.0-0.25
    dx = 0.5
    ys = xs
    dy = dx

    # Save the BINID, and x and y coordinates
    binid = dapf.hdu['DRPS'].data['BINID']
    X = dapf.hdu['DRPS'].data['XPOS']
    Y = dapf.hdu['DRPS'].data['YPOS']
        
    # Have to adjust the x and y values given the error in DAP v1_0_0
    if mpl == 3:
        X += 0.5
        Y -= 0.5

    # Create and save the data for all the maps
    p = 0

    #-------------------------------------------------------------------
    # DRPS
    for c in drps_col:
        z = dapf.hdu['DRPS'].data[c].ravel()
        place_values_in_map(-X, Y, z, xs, dx, ys, dy, plane)
        output_data[p,:,:] = plane[:,:]
        p += 1
        print('Finished map {0}/{1}'.format(p,nmaps), end="\r")
        newhdr['P'+'{0}'.format(p).zfill(ndig)] = 'DRPS:{0}'.format(c)
    #-------------------------------------------------------------------

    #-------------------------------------------------------------------
    # BINS
    for c in bins_col:
        z = place_values_in_bins(binid, dapf.hdu['BINS'].data[c])
        place_values_in_map(-X, Y, z, xs, dx, ys, dy, plane)
        output_data[p,:,:] = plane[:,:]
        p += 1
        print('Finished map {0}/{1}'.format(p,nmaps), end="\r")
        newhdr['P'+'{0}'.format(p).zfill(ndig)] = 'BINS:{0}'.format(c)
    #-------------------------------------------------------------------

    #-------------------------------------------------------------------
    # STFIT
    for i in range(nstkin):
        for c in stfit_vec:
            z = place_values_in_bins(binid, dapf.hdu['STFIT'].data[c][:,i])
            place_values_in_map(-X, Y, z, xs, dx, ys, dy, plane)
            output_data[p,:,:] = plane[:,:]
            p += 1
            print('Finished map {0}/{1}'.format(p,nmaps), end="\r")
            newhdr['P'+'{0}'.format(p).zfill(ndig)] = 'STFIT:{0}:{1}'.format(c,i)
    for c in stfit_col:
        z = place_values_in_bins(binid, dapf.hdu['STFIT'].data[c])
        place_values_in_map(-X, Y, z, xs, dx, ys, dy, plane)
        output_data[p,:,:] = plane[:,:]
        p += 1
        print('Finished map {0}/{1}'.format(p,nmaps), end="\r")
        newhdr['P'+'{0}'.format(p).zfill(ndig)] = 'STFIT:{0}'.format(c)
    #-------------------------------------------------------------------

    #-------------------------------------------------------------------
    # ELMMNT
    for i in range(nelmmnt):
        for c in elmmnt_vec:
            z = place_values_in_bins(binid, dapf.hdu['ELMMNT'].data[c][:,i])
            place_values_in_map(-X, Y, z, xs, dx, ys, dy, plane)
            output_data[p,:,:] = plane[:,:]
            p += 1
            print('Finished map {0}/{1}'.format(p,nmaps), end="\r")
            newhdr['P'+'{0}'.format(p).zfill(ndig)] = \
                'ELMMNT:{0}:{1}'.format(c, dapf.hdu['ELBAND'].data['ELNAME'][i])
    #-------------------------------------------------------------------

    #-------------------------------------------------------------------
    # ELOFIT: EW
    for i in range(nelowgt*nelokin):
        for c in elofit_kin_vec:
            c += '_EW'
            z = place_values_in_bins(binid, dapf.hdu['ELOFIT'].data[c][:,i])
            place_values_in_map(-X, Y, z, xs, dx, ys, dy, plane)
            output_data[p,:,:] = plane[:,:]
            p += 1
            print('Finished map {0}/{1}'.format(p,nmaps), end="\r")
            newhdr['P'+'{0}'.format(p).zfill(ndig)] = 'ELOFIT:{0}:{1}'.format(c, i)

    for c in elofit_col:
        c += '_EW'
        z = place_values_in_bins(binid, dapf.hdu['ELOFIT'].data[c])
        place_values_in_map(-X, Y, z, xs, dx, ys, dy, plane)
        output_data[p,:,:] = plane[:,:]
        p += 1
        print('Finished map {0}/{1}'.format(p,nmaps), end="\r")
        newhdr['P'+'{0}'.format(p).zfill(ndig)] = 'ELOFIT:{0}'.format(c)

    for i in range(neml):
        for c, nc in zip(elofit_ind_vec, elofit_ind_vec_nel):
            c += '_EW'
            if nc == 1:
                z = place_values_in_bins(binid, dapf.hdu['ELOFIT'].data[c][:,i])
                place_values_in_map(-X, Y, z, xs, dx, ys, dy, plane)
                output_data[p,:,:] = plane[:,:]
                p += 1
                print('Finished map {0}/{1}'.format(p,nmaps), end="\r")
                newhdr['P'+'{0}'.format(p).zfill(ndig)] = \
                    'ELOFIT:{0}:{1}'.format(c, dapf.hdu['ELOPAR'].data['ELNAME'][i])
            else:
                for j in range(nc):
                    z = place_values_in_bins(binid, dapf.hdu['ELOFIT'].data[c][:,j,i])
                    place_values_in_map(-X, Y, z, xs, dx, ys, dy, plane)
                    output_data[p,:,:] = plane[:,:]
                    p += 1
                    print('Finished map {0}/{1}'.format(p,nmaps), end="\r")
                    newhdr['P'+'{0}'.format(p).zfill(ndig)] = \
                        'ELOFIT:{0}:{1}:{2}'.format(c, j, dapf.hdu['ELOPAR'].data['ELNAME'][i])
    #-------------------------------------------------------------------

    #-------------------------------------------------------------------
    # ELOFIT: FB
    for i in range(nelowgt*nelokin):
        for c in elofit_kin_vec:
            c += '_FB'
            z = place_values_in_bins(binid, dapf.hdu['ELOFIT'].data[c][:,i])
            place_values_in_map(-X, Y, z, xs, dx, ys, dy, plane)
            output_data[p,:,:] = plane[:,:]
            p += 1
            print('Finished map {0}/{1}'.format(p,nmaps), end="\r")
            newhdr['P'+'{0}'.format(p).zfill(ndig)] = 'ELOFIT:{0}:{1}'.format(c, i)

    for c in elofit_col:
        c += '_FB'
        z = place_values_in_bins(binid, dapf.hdu['ELOFIT'].data[c])
        place_values_in_map(-X, Y, z, xs, dx, ys, dy, plane)
        output_data[p,:,:] = plane[:,:]
        p += 1
        print('Finished map {0}/{1}'.format(p,nmaps), end="\r")
        newhdr['P'+'{0}'.format(p).zfill(ndig)] = 'ELOFIT:{0}'.format(c)

    for i in range(neml):
        for c, nc in zip(elofit_ind_vec, elofit_ind_vec_nel):
            c += '_FB'
            if nc == 1:
                z = place_values_in_bins(binid, dapf.hdu['ELOFIT'].data[c][:,i])
                place_values_in_map(-X, Y, z, xs, dx, ys, dy, plane)
                output_data[p,:,:] = plane[:,:]
                p += 1
                print('Finished map {0}/{1}'.format(p,nmaps), end="\r")
                newhdr['P'+'{0}'.format(p).zfill(ndig)] = \
                    'ELOFIT:{0}:{1}'.format(c, dapf.hdu['ELOPAR'].data['ELNAME'][i])
            else:
                for j in range(nc):
                    z = place_values_in_bins(binid, dapf.hdu['ELOFIT'].data[c][:,j,i])
                    place_values_in_map(-X, Y, z, xs, dx, ys, dy, plane)
                    output_data[p,:,:] = plane[:,:]
                    p += 1
                    print('Finished map {0}/{1}'.format(p,nmaps), end="\r")
                    newhdr['P'+'{0}'.format(p).zfill(ndig)] = \
                        'ELOFIT:{0}:{1}:{2}'.format(c, j, dapf.hdu['ELOPAR'].data['ELNAME'][i])
    #-------------------------------------------------------------------

    #-------------------------------------------------------------------
    # SINDX
    for i in range(nsindx):
        for c in sindx_vec:
            z = place_values_in_bins(binid, dapf.hdu['SINDX'].data[c][:,i])
            place_values_in_map(-X, Y, z, xs, dx, ys, dy, plane)
            output_data[p,:,:] = plane[:,:]
            p += 1
            print('Finished map {0}/{1}'.format(p,nmaps), end="\r")
            newhdr['P'+'{0}'.format(p).zfill(ndig)] = \
                'SINDX:{0}:{1}'.format(c, dapf.hdu['SIPAR'].data['SINAME'][i])
    #-------------------------------------------------------------------

    print('Finished map: DONE                  ')

    # Write the datacube
    print('Writing: {0}'.format(output_file))
    hdu = fits.PrimaryHDU(numpy.transpose(output_data, (0,2,1)), header=newhdr)
#    if os.path.isfile(output_file):
#        remove(output_file)
    hdu.writeto(output_file, clobber=os.path.isfile(output_file))


#-----------------------------------------------------------------------------

if __name__ == '__main__':
    t = clock()

    parser = ArgumentParser()

    parser.add_argument('plate', type=int, help='plate ID to process')
    parser.add_argument('ifudesign', type=int, help='IFU design to process')
    parser.add_argument('mpl', type=int, help='MPL number (3 or 4)')
    parser.add_argument('output_file', type=str, help='Name for output file')
    parser.add_argument('-b', '--bintype', type=str,
                        help='Binning type to process: NONE(def) or STON', default='NONE')
    parser.add_argument('-i', '--plan_index', type=int,
                        help='Index of plan used by DAP (used in the DAP output file name)',
                        default=None)
    parser.add_argument('-d', '--directory_path', type=str, help='Path to DAP output file',
                        default=None)

    arg = parser.parse_args()

    if os.path.isfile(arg.output_file):
        print('WARNING: Overwriting existing file {0}!'.format(arg.output_file))

    create_dap_fits_cube(arg.plate, arg.ifudesign, arg.bintype, arg.output_file, mpl=arg.mpl,
                         plan_index=arg.plan_index, directory_path=arg.directory_path)
#    try:
#        create_dap_fits_cube(arg.plate, arg.ifudesign, arg.bintype, arg.output_file, mpl=arg.mpl,
#                             plan_index=arg.plan_index, directory_path=arg.directory_path)
#    except Exception as e:
#        print_frame('Exception')
#        print(e)
#        parser.print_usage()
 
    print('Elapsed time: {0} seconds'.format(clock() - t))



