#!/usr/bin/env python3

from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import sys
if sys.version > '3':
    long = int

import numpy
from os import remove
import os.path
from time import clock

from astropy.io import fits
from mangadap.dapfile import dapfile

#-----------------------------------------------------------------------------

def place_values_in_map(x, y, z, xs, dx, ys, dy, plane):
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
    nspax = binid.size
    z = numpy.zeros(nspax, dtype=numpy.float64)
    for i in range(nspax):
        if binid[i] < 0:
            continue
        z[i] = zz[binid[i]]
    return z


def create_dap_fits_cube(plate, ifudesign, bintype, ofits, directory_path):

    # Check the binning type and set the iteration number
    # !! MPL-3 SPECIFIC !!
    if bintype == 'STON':
        niter = 1
    elif bintype == 'NONE':
        niter = 2
    else:
        raise KeyError('Binning type must be STON or NONE.')

    print('Attempting to read DAP output file ...')
    dapf = dapfile(plate, ifudesign, 'CUBE', bintype, niter, directory_path=directory_path)
    print('Found: {0}'.format(dapf.file_path()))

    # Initialize the header
    newhdr = dapf.hdu['PRIMARY'].header

    # Remove the wavelength axis and replace it with just the pixel coordinates
    newhdr.remove('CTYPE3')
    newhdr.remove('CUNIT3')
    newhdr['CRPIX3'] = 1
    newhdr['CRVAL3'] = 1.
    newhdr['CD3_3'] = 1.

    # Determine the number of maps to produce
    nstkin = dapf.hdu['STFIT'].data['KIN'].shape[1]
    nelokin = dapf.hdu['ELOFIT'].data['KIN_EW'].shape[1]
    neml = dapf.hdu['ELOFIT'].data['ELOMIT_EW'].shape[1]
    nsindx = dapf.hdu['SINDX'].data['INDX'].shape[1]
    nmaps = 6 + nstkin*2 + 1 + 2*(nelokin*2 + neml*(8 + nelokin*2)) + nsindx*5

    ndig = int(numpy.log10(nmaps))+1

    # DRP CUBE file is always assumed to be square in the number of
    # on-sky pixels!
    npix = int(numpy.sqrt(dapf.hdu['DRPS'].header['NAXIS2']))

    print('Number of stellar kinematic moments:         {0}'.format(nstkin))
    print('Number of emission-line kinematic moments:   {0}'.format(nelokin))
    print('Number of fitted emission lines:             {0}'.format(neml))
    print('Number of spectral indices:                  {0}'.format(nsindx))
    print('Total number of maps to produce:             {0}'.format(nmaps))
    print('Size of each map:                            {0} pix x {0} pix'.format(npix))

    output_data = numpy.zeros((nmaps, npix, npix), dtype=numpy.float64)
    plane = numpy.zeros((npix, npix), dtype=numpy.float64)

    xs = -npix/4.0-0.25
    dx = 0.5
    ys = xs
    dy = dx

    # Have to adjust the x and y values given the error in DAP v1_0_0
    X = dapf.hdu['DRPS'].data['XPOS'] + 0.5
    Y = dapf.hdu['DRPS'].data['YPOS'] - 0.5

    # Create and save the data for all the maps
    p = 0
    z = dapf.hdu['DRPS'].data['FGOODPIX'].ravel()
    place_values_in_map(-X, Y, z, xs, dx, ys, dy, plane)
    output_data[p,:,:] = plane[:,:]
    p += 1
    print('Finished map {0}/{1}'.format(p,nmaps), end="\r")
    newhdr['P'+'{0}'.format(p).zfill(ndig)] = 'DRPS:FGOODPIX'

    z = dapf.hdu['DRPS'].data['SIGNAL']
    place_values_in_map(-X, Y, z, xs, dx, ys, dy, plane)
    output_data[p,:,:] = plane[:,:]
    p += 1
    print('Finished map {0}/{1}'.format(p,nmaps), end="\r")
    newhdr['P'+'{0}'.format(p).zfill(ndig)] = 'DRPS:SIGNAL'

    z = dapf.hdu['DRPS'].data['NOISE']
    place_values_in_map(-X, Y, z, xs, dx, ys, dy, plane)
    output_data[p,:,:] = plane[:,:]
    p += 1
    print('Finished map {0}/{1}'.format(p,nmaps), end="\r")
    newhdr['P'+'{0}'.format(p).zfill(ndig)] = 'DRPS:NOISE'

    binid = dapf.hdu['DRPS'].data['BINID']
    z = dapf.hdu['DRPS'].data['BINID']
    place_values_in_map(-X, Y, z, xs, dx, ys, dy, plane)
    output_data[p,:,:] = plane[:,:]
    p += 1
    print('Finished map {0}/{1}'.format(p,nmaps), end="\r")
    newhdr['P'+'{0}'.format(p).zfill(ndig)] = 'DRPS:BINID'

    z = dapf.hdu['DRPS'].data['BINW']
    place_values_in_map(-X, Y, z, xs, dx, ys, dy, plane)
    output_data[p,:,:] = plane[:,:]
    p += 1
    print('Finished map {0}/{1}'.format(p,nmaps), end="\r")
    newhdr['P'+'{0}'.format(p).zfill(ndig)] = 'DRPS:BINW'

    z = place_values_in_bins(binid, dapf.hdu['BINS'].data['BINSN'])
    place_values_in_map(-X, Y, z, xs, dx, ys, dy, plane)
    output_data[p,:,:] = plane[:,:]
    p += 1
    print('Finished map {0}/{1}'.format(p,nmaps), end="\r")
    newhdr['P'+'{0}'.format(p).zfill(ndig)] = 'BINS:BINSN'

    for i in range(nstkin):
        z = place_values_in_bins(binid, dapf.hdu['STFIT'].data['KIN'][:,i])
        place_values_in_map(-X, Y, z, xs, dx, ys, dy, plane)
        output_data[p,:,:] = plane[:,:]
        p += 1
        print('Finished map {0}/{1}'.format(p,nmaps), end="\r")
        newhdr['P'+'{0}'.format(p).zfill(ndig)] = 'STFIT:KIN:{0}'.format(i)
        
        z = place_values_in_bins(binid, dapf.hdu['STFIT'].data['KINERR'][:,i])
        place_values_in_map(-X, Y, z, xs, dx, ys, dy, plane)
        output_data[p,:,:] = plane[:,:]
        p += 1
        print('Finished map {0}/{1}'.format(p,nmaps), end="\r")
        newhdr['P'+'{0}'.format(p).zfill(ndig)] = 'STFIT:KINERR:{0}'.format(i)

    z = place_values_in_bins(binid, dapf.hdu['STFIT'].data['RCHI2'])
    place_values_in_map(-X, Y, z, xs, dx, ys, dy, plane)
    output_data[p,:,:] = plane[:,:]
    p += 1
    print('Finished map {0}/{1}'.format(p,nmaps), end="\r")
    newhdr['P'+'{0}'.format(p).zfill(ndig)] = 'STFIT:RCHI2'

    for i in range(nelokin):
        z = place_values_in_bins(binid, dapf.hdu['ELOFIT'].data['KIN_EW'][:,i])
        place_values_in_map(-X, Y, z, xs, dx, ys, dy, plane)
        output_data[p,:,:] = plane[:,:]
        p += 1
        print('Finished map {0}/{1}'.format(p,nmaps), end="\r")
        newhdr['P'+'{0}'.format(p).zfill(ndig)] = 'ELOFIT:KIN_EW:{0}'.format(i)

        z = place_values_in_bins(binid, dapf.hdu['ELOFIT'].data['KINERR_EW'][:,i])
        place_values_in_map(-X, Y, z, xs, dx, ys, dy, plane)
        output_data[p,:,:] = plane[:,:]
        p += 1
        print('Finished map {0}/{1}'.format(p,nmaps), end="\r")
        newhdr['P'+'{0}'.format(p).zfill(ndig)] = 'ELOFIT:KINERR_EW:{0}'.format(i)

    for i in range(neml):
        z = place_values_in_bins(binid, dapf.hdu['ELOFIT'].data['ELOMIT_EW'][:,i])
        place_values_in_map(-X, Y, z, xs, dx, ys, dy, plane)
        output_data[p,:,:] = plane[:,:]
        p += 1
        print('Finished map {0}/{1}'.format(p,nmaps), end="\r")
        newhdr['P'+'{0}'.format(p).zfill(ndig)] = \
                'ELOFIT:ELOMIT_EW:{0}'.format(dapf.hdu['ELOPAR'].data['ELNAME'][i])

        z = place_values_in_bins(binid, dapf.hdu['ELOFIT'].data['AMPL_EW'][:,i])
        place_values_in_map(-X, Y, z, xs, dx, ys, dy, plane)
        output_data[p,:,:] = plane[:,:]
        p += 1
        print('Finished map {0}/{1}'.format(p,nmaps), end="\r")
        newhdr['P'+'{0}'.format(p).zfill(ndig)] = \
                'ELOFIT:AMPL_EW:{0}'.format(dapf.hdu['ELOPAR'].data['ELNAME'][i])

        z = place_values_in_bins(binid, dapf.hdu['ELOFIT'].data['AMPLERR_EW'][:,i])
        place_values_in_map(-X, Y, z, xs, dx, ys, dy, plane)
        output_data[p,:,:] = plane[:,:]
        p += 1
        print('Finished map {0}/{1}'.format(p,nmaps), end="\r")
        newhdr['P'+'{0}'.format(p).zfill(ndig)] = \
                'ELOFIT:AMPLERR_EW:{0}'.format(dapf.hdu['ELOPAR'].data['ELNAME'][i])

        for j in range(nelokin):
            z = place_values_in_bins(binid, dapf.hdu['ELOFIT'].data['IKIN_EW'][:,j,i])
            place_values_in_map(-X, Y, z, xs, dx, ys, dy, plane)
            output_data[p,:,:] = plane[:,:]
            p += 1
            print('Finished map {0}/{1}'.format(p,nmaps), end="\r")
            newhdr['P'+'{0}'.format(p).zfill(ndig)] = \
                    'ELOFIT:IKIN_EW:{0}:{1}'.format(j, dapf.hdu['ELOPAR'].data['ELNAME'][i])

            z = place_values_in_bins(binid, dapf.hdu['ELOFIT'].data['IKINERR_EW'][:,j,i])
            place_values_in_map(-X, Y, z, xs, dx, ys, dy, plane)
            output_data[p,:,:] = plane[:,:]
            p += 1
            print('Finished map {0}/{1}'.format(p,nmaps), end="\r")
            newhdr['P'+'{0}'.format(p).zfill(ndig)] = \
                    'ELOFIT:IKINERR_EW:{0}:{1}'.format(j, dapf.hdu['ELOPAR'].data['ELNAME'][i])

        z = place_values_in_bins(binid, dapf.hdu['ELOFIT'].data['SINST_EW'][:,i])
        place_values_in_map(-X, Y, z, xs, dx, ys, dy, plane)
        output_data[p,:,:] = plane[:,:]
        p += 1
        print('Finished map {0}/{1}'.format(p,nmaps), end="\r")
        newhdr['P'+'{0}'.format(p).zfill(ndig)] = \
                'ELOFIT:SINST_EW:{0}'.format(dapf.hdu['ELOPAR'].data['ELNAME'][i])

        z = place_values_in_bins(binid, dapf.hdu['ELOFIT'].data['FLUX_EW'][:,i])
        place_values_in_map(-X, Y, z, xs, dx, ys, dy, plane)
        output_data[p,:,:] = plane[:,:]
        p += 1
        print('Finished map {0}/{1}'.format(p,nmaps), end="\r")
        newhdr['P'+'{0}'.format(p).zfill(ndig)] = \
                'ELOFIT:FLUX_EW:{0}'.format(dapf.hdu['ELOPAR'].data['ELNAME'][i])

        z = place_values_in_bins(binid, dapf.hdu['ELOFIT'].data['FLUXERR_EW'][:,i])
        place_values_in_map(-X, Y, z, xs, dx, ys, dy, plane)
        output_data[p,:,:] = plane[:,:]
        p += 1
        print('Finished map {0}/{1}'.format(p,nmaps), end="\r")
        newhdr['P'+'{0}'.format(p).zfill(ndig)] = \
                'ELOFIT:FLUXERR_EW:{0}'.format(dapf.hdu['ELOPAR'].data['ELNAME'][i])

        z = place_values_in_bins(binid, dapf.hdu['ELOFIT'].data['EW_EW'][:,i])
        place_values_in_map(-X, Y, z, xs, dx, ys, dy, plane)
        output_data[p,:,:] = plane[:,:]
        p += 1
        print('Finished map {0}/{1}'.format(p,nmaps), end="\r")
        newhdr['P'+'{0}'.format(p).zfill(ndig)] = \
                'ELOFIT:EW_EW:{0}'.format(dapf.hdu['ELOPAR'].data['ELNAME'][i])

        z = place_values_in_bins(binid, dapf.hdu['ELOFIT'].data['EWERR_EW'][:,i])
        place_values_in_map(-X, Y, z, xs, dx, ys, dy, plane)
        output_data[p,:,:] = plane[:,:]
        p += 1
        print('Finished map {0}/{1}'.format(p,nmaps), end="\r")
        newhdr['P'+'{0}'.format(p).zfill(ndig)] = \
                'ELOFIT:EWERR_EW:{0}'.format(dapf.hdu['ELOPAR'].data['ELNAME'][i])

    for i in range(nelokin):
        z = place_values_in_bins(binid, dapf.hdu['ELOFIT'].data['KIN_FB'][:,i])
        place_values_in_map(-X, Y, z, xs, dx, ys, dy, plane)
        output_data[p,:,:] = plane[:,:]
        p += 1
        print('Finished map {0}/{1}'.format(p,nmaps), end="\r")
        newhdr['P'+'{0}'.format(p).zfill(ndig)] = 'ELOFIT:KIN_FB:{0}'.format(i)

        z = place_values_in_bins(binid, dapf.hdu['ELOFIT'].data['KINERR_FB'][:,i])
        place_values_in_map(-X, Y, z, xs, dx, ys, dy, plane)
        output_data[p,:,:] = plane[:,:]
        p += 1
        print('Finished map {0}/{1}'.format(p,nmaps), end="\r")
        newhdr['P'+'{0}'.format(p).zfill(ndig)] = 'ELOFIT:KINERR_FB:{0}'.format(i)

    for i in range(neml):
        z = place_values_in_bins(binid, dapf.hdu['ELOFIT'].data['ELOMIT_FB'][:,i])
        place_values_in_map(-X, Y, z, xs, dx, ys, dy, plane)
        output_data[p,:,:] = plane[:,:]
        p += 1
        print('Finished map {0}/{1}'.format(p,nmaps), end="\r")
        newhdr['P'+'{0}'.format(p).zfill(ndig)] = \
                'ELOFIT:ELOMIT_FB:{0}'.format(dapf.hdu['ELOPAR'].data['ELNAME'][i])

        z = place_values_in_bins(binid, dapf.hdu['ELOFIT'].data['AMPL_FB'][:,i])
        place_values_in_map(-X, Y, z, xs, dx, ys, dy, plane)
        output_data[p,:,:] = plane[:,:]
        p += 1
        print('Finished map {0}/{1}'.format(p,nmaps), end="\r")
        newhdr['P'+'{0}'.format(p).zfill(ndig)] = \
                'ELOFIT:AMPL_FB:{0}'.format(dapf.hdu['ELOPAR'].data['ELNAME'][i])

        z = place_values_in_bins(binid, dapf.hdu['ELOFIT'].data['AMPLERR_FB'][:,i])
        place_values_in_map(-X, Y, z, xs, dx, ys, dy, plane)
        output_data[p,:,:] = plane[:,:]
        p += 1
        print('Finished map {0}/{1}'.format(p,nmaps), end="\r")
        newhdr['P'+'{0}'.format(p).zfill(ndig)] = \
                'ELOFIT:AMPLERR_FB:{0}'.format(dapf.hdu['ELOPAR'].data['ELNAME'][i])

        for j in range(nelokin):
            z = place_values_in_bins(binid, dapf.hdu['ELOFIT'].data['IKIN_FB'][:,j,i])
            place_values_in_map(-X, Y, z, xs, dx, ys, dy, plane)
            output_data[p,:,:] = plane[:,:]
            p += 1
            print('Finished map {0}/{1}'.format(p,nmaps), end="\r")
            newhdr['P'+'{0}'.format(p).zfill(ndig)] = \
                    'ELOFIT:IKIN_FB:{0}:{1}'.format(j, dapf.hdu['ELOPAR'].data['ELNAME'][i])

            z = place_values_in_bins(binid, dapf.hdu['ELOFIT'].data['IKINERR_FB'][:,j,i])
            place_values_in_map(-X, Y, z, xs, dx, ys, dy, plane)
            output_data[p,:,:] = plane[:,:]
            p += 1
            print('Finished map {0}/{1}'.format(p,nmaps), end="\r")
            newhdr['P'+'{0}'.format(p).zfill(ndig)] = \
                    'ELOFIT:IKINERR_FB:{0}:{1}'.format(j, dapf.hdu['ELOPAR'].data['ELNAME'][i])

        z = place_values_in_bins(binid, dapf.hdu['ELOFIT'].data['SINST_FB'][:,i])
        place_values_in_map(-X, Y, z, xs, dx, ys, dy, plane)
        output_data[p,:,:] = plane[:,:]
        p += 1
        print('Finished map {0}/{1}'.format(p,nmaps), end="\r")
        newhdr['P'+'{0}'.format(p).zfill(ndig)] = \
                'ELOFIT:SINST_FB:{0}'.format(dapf.hdu['ELOPAR'].data['ELNAME'][i])

        z = place_values_in_bins(binid, dapf.hdu['ELOFIT'].data['FLUX_FB'][:,i])
        place_values_in_map(-X, Y, z, xs, dx, ys, dy, plane)
        output_data[p,:,:] = plane[:,:]
        p += 1
        print('Finished map {0}/{1}'.format(p,nmaps), end="\r")
        newhdr['P'+'{0}'.format(p).zfill(ndig)] = \
                'ELOFIT:FLUX_FB:{0}'.format(dapf.hdu['ELOPAR'].data['ELNAME'][i])

        z = place_values_in_bins(binid, dapf.hdu['ELOFIT'].data['FLUXERR_FB'][:,i])
        place_values_in_map(-X, Y, z, xs, dx, ys, dy, plane)
        output_data[p,:,:] = plane[:,:]
        p += 1
        print('Finished map {0}/{1}'.format(p,nmaps), end="\r")
        newhdr['P'+'{0}'.format(p).zfill(ndig)] = \
                'ELOFIT:FLUXERR_FB:{0}'.format(dapf.hdu['ELOPAR'].data['ELNAME'][i])

        z = place_values_in_bins(binid, dapf.hdu['ELOFIT'].data['EW_FB'][:,i])
        place_values_in_map(-X, Y, z, xs, dx, ys, dy, plane)
        output_data[p,:,:] = plane[:,:]
        p += 1
        print('Finished map {0}/{1}'.format(p,nmaps), end="\r")
        newhdr['P'+'{0}'.format(p).zfill(ndig)] = \
                'ELOFIT:EW_FB:{0}'.format(dapf.hdu['ELOPAR'].data['ELNAME'][i])

        z = place_values_in_bins(binid, dapf.hdu['ELOFIT'].data['EWERR_FB'][:,i])
        place_values_in_map(-X, Y, z, xs, dx, ys, dy, plane)
        output_data[p,:,:] = plane[:,:]
        p += 1
        print('Finished map {0}/{1}'.format(p,nmaps), end="\r")
        newhdr['P'+'{0}'.format(p).zfill(ndig)] = \
                'ELOFIT:EWERR_FB:{0}'.format(dapf.hdu['ELOPAR'].data['ELNAME'][i])

    for i in range(nsindx):
        z = place_values_in_bins(binid, dapf.hdu['SINDX'].data['SIOMIT'][:,i])
        place_values_in_map(-X, Y, z, xs, dx, ys, dy, plane)
        output_data[p,:,:] = plane[:,:]
        p += 1
        print('Finished map {0}/{1}'.format(p,nmaps), end="\r")
        newhdr['P'+'{0}'.format(p).zfill(ndig)] = \
                'SINDX:SIOMIT:{0}'.format(dapf.hdu['SIPAR'].data['SINAME'][i])

        z = place_values_in_bins(binid, dapf.hdu['SINDX'].data['INDX'][:,i])
        place_values_in_map(-X, Y, z, xs, dx, ys, dy, plane)
        output_data[p,:,:] = plane[:,:]
        p += 1
        print('Finished map {0}/{1}'.format(p,nmaps), end="\r")
        newhdr['P'+'{0}'.format(p).zfill(ndig)] = \
                'SINDX:INDX:{0}'.format(dapf.hdu['SIPAR'].data['SINAME'][i])

        z = place_values_in_bins(binid, dapf.hdu['SINDX'].data['INDXERR'][:,i])
        place_values_in_map(-X, Y, z, xs, dx, ys, dy, plane)
        output_data[p,:,:] = plane[:,:]
        p += 1
        print('Finished map {0}/{1}'.format(p,nmaps), end="\r")
        newhdr['P'+'{0}'.format(p).zfill(ndig)] = \
                'SINDX:INDXERR:{0}'.format(dapf.hdu['SIPAR'].data['SINAME'][i])

        z = place_values_in_bins(binid, dapf.hdu['SINDX'].data['INDX_OTPL'][:,i])
        place_values_in_map(-X, Y, z, xs, dx, ys, dy, plane)
        output_data[p,:,:] = plane[:,:]
        p += 1
        print('Finished map {0}/{1}'.format(p,nmaps), end="\r")
        newhdr['P'+'{0}'.format(p).zfill(ndig)] = \
                'SINDX:INDX_OTPL:{0}'.format(dapf.hdu['SIPAR'].data['SINAME'][i])

        z = place_values_in_bins(binid, dapf.hdu['SINDX'].data['INDX_BOTPL'][:,i])
        place_values_in_map(-X, Y, z, xs, dx, ys, dy, plane)
        output_data[p,:,:] = plane[:,:]
        p += 1
        print('Finished map {0}/{1}'.format(p,nmaps), end="\r")
        newhdr['P'+'{0}'.format(p).zfill(ndig)] = \
                'SINDX:INDX_BOTPL:{0}'.format(dapf.hdu['SIPAR'].data['SINAME'][i])

    print('Finished map: DONE                  ')

    # Write the datacube
    print('Writing: {0}'.format(ofits))
    hdu = fits.PrimaryHDU(numpy.transpose(output_data, (0,2,1)), header=newhdr)
    if os.path.isfile(ofits):
        remove(ofits)
    hdu.writeto(ofits, clobber=True)


#-----------------------------------------------------------------------------

if __name__ == '__main__':
    t = clock()

    narg = len(sys.argv)
    if narg != 5 and narg != 6 :
        print('Usage: dap_data_to_fits_cube.py <plate> <ifudesign> <STON or NONE> <output file>'
              ' <file directory>')
        print('  or')
        print('       dap_data_to_fits_cube.py <plate> <ifudesign> <STON or NONE> <output file>')

        exit(1)

    ofile = sys.argv[4]
    if os.path.isfile(ofile):
        print('WARNING: Overwriting existing file {0}!'.format(ofile))

    directory_path = None if narg == 5 else str(sys.argv[5])

    create_dap_fits_cube(int(sys.argv[1]), int(sys.argv[2]), str(sys.argv[3]), ofile,
                         directory_path)

    print('Elapsed time: {0} seconds'.format(clock() - t))



