#!/usr/bin/env python3

import numpy
from os import remove
import os.path
from time import clock
import sys

from mangadap.drpfile import drpfile

#-----------------------------------------------------------------------------

def calculate_covariance_cube(plate, ifudesign, nchannels, ofile, directory_path):

    print('     PLATE: {0}'.format(plate))
    print(' IFUDESIGN: {0}'.format(ifudesign))

    # Access the DRP RSS file
    print('Attempting to open RSS file:')
    drpf = drpfile(plate, ifudesign, 'RSS', read=True, directory_path=directory_path)
    print('     FOUND: {0}'.format(drpf.file_path()))

    nw = drpf.hdu['FLUX'].data.shape[1]                     # Number of wavelength channels

    if nchannels >= nw or nchannels == 0:
        print('Calculating full covariance cube ...')
        C = drpf.covariance_cube()
    elif nchannels == 1:
        print('Only one channel selected.  Calculating covariance at central channel...')
        C = drpf.covariance_matrix(nw/2)
    else:
        print('Calculating covariance in {0} evenly spaced wavelength channels...'.format(nchannels))
#        planes = numpy.linspace(0, nw-1, num=nchannels, dtype=numpy.int)
        planes = numpy.linspace(0, nw-1, num=nchannels).astype(numpy.int)
        print(planes)
        C = drpf.covariance_cube(channels=planes)
    print('... done.')

    print('Writing data to {0}.'.format(ofile))
    C.write(ofile, clobber=True)            # Write the data

#-----------------------------------------------------------------------------

if __name__ == '__main__':
    t = clock()

    narg = len(sys.argv)
    if narg != 5 and narg != 6 :
        print('Usage: calculate_covariance.py <plate> <ifudesign> <num channels> <output file> '
              '<file directory>')
        print('  or')
        print('       calculate_covariance.py <plate> <ifudesign> <num channels> <output file>')
        exit(1)

    ofile = sys.argv[4]
    if os.path.isfile(ofile):
        print('WARNING: Overwriting existing file {0}!'.format(ofile))

    directory_path = None if narg == 5 else str(sys.argv[5])

    calculate_covariance_cube(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]),
                              str(sys.argv[4]), directory_path)

    print('Elapsed time: {0} seconds'.format(clock() - t))



