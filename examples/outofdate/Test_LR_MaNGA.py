#!/usr/bin/env python3
""" LambdaR for MaNGA
 Version 0.0.2  - Eric, after editing and adding the relevant functions
 Version 0.0.1  - Kyle 
"""

import sys

from os import remove
import os.path

import numpy
from time import clock
from mangadap.util.exception_tools import print_frame
from mangadap.dapfile import dapfile

from mangadap.contrib.LambdaR_2D_forMaNGA import Derive_LR_VS_Profiles, derive_unbinned_field

from matplotlib import pyplot

# Pull the necessary data out of the DAP file
def dap_data(dapf):
    ell = dapf.nsa_ellipticity()
    pa = dapf.guess_position_angle()
    xpos = dapf.hdu['DRPS'].data['XPOS']
    ypos = dapf.hdu['DRPS'].data['YPOS']
    signal = dapf.hdu['DRPS'].data['SIGNAL']
    binid = dapf.hdu['DRPS'].data['BINID']
    binx = dapf.hdu['BINS'].data['BINXRL']
    biny = dapf.hdu['BINS'].data['BINYRU']
    binf = dapf.hdu['BINS'].data['BINF']
    stellar_v = dapf.hdu['STFIT'].data['KIN'][:,0]
    stellar_sigma = dapf.hdu['STFIT'].data['KIN'][:,1]
    return ell, pa, xpos, ypos, signal, binid, binx, biny, stellar_v, stellar_sigma

#-----------------------------------------------------------------------------

if __name__ == '__main__':
    # Start the clock
    t = clock()

    # Check number of command-line arguments is correct
    narg = 6
    if len(sys.argv) != narg:
        print('Usage: Test_LR_MaNGA <plate> <ifudesign> <mode (CUBE or RSS)> <bintype (e.g.' \
              ' STON)> <iter>')
        raise Exception('Incorrect number of arguments!')

    # Try to open the DAP output file
    try:
        dapf = dapfile(plate=sys.argv[1], ifudesign=sys.argv[2], mode=sys.argv[3],
                       bintype=sys.argv[4], niter=sys.argv[5], read=True)
    except Exception as e:
        print(e)
        print_frame('Exception')
        exit(1)

    # Try to open the associated parameter file
    try: 
        dapf.read_par()
    except Exception as e:
        print(e)
        print_frame('Exception')
        exit(1)

    # Get the necessary data
    ell, pa, xpos, ypos, signal, binid, binx, biny, stellar_v, stellar_sigma = dap_data(dapf)

    # Run Eric's functions to produce output data
    sel_spaxel = (binid >= 0)
    xunb, yunb, V_unb = derive_unbinned_field(binx, biny, stellar_v, 
                                              xpos[sel_spaxel], ypos[sel_spaxel])
    xunb, yunb, S_unb = derive_unbinned_field(binx, biny, stellar_sigma, 
                                              xpos[sel_spaxel], ypos[sel_spaxel])

    # Deriving the maximum Radius and creating an array to sample the radii for LambdaR
    Max_Radius = numpy.max(numpy.sqrt(xpos**2 + ypos**2))
    Min_Radius = 1.0  # Minimum Radius to consider [arcsec]
    N_Radius = 30 # number of points to consider 
    print('Maximum radius is {0}'.format(Max_Radius))

    Result = Derive_LR_VS_Profiles(xunb, yunb, signal[sel_spaxel], V_unb, S_unb, 
                                   Eps=ell, PA=pa, verbose=False, Min_Radius=Min_Radius,
                                   Max_Radius=Max_Radius, N_Radius=N_Radius)

    print('Maximum radius is {0}'.format(Result.Max_Radius))

    print(Result.R_EpsPA.shape)
    print(Result.SemiMajorRadius.shape)
    print(Result.LambdaR.shape)
    print(Result.VS.shape)

    pyplot.plot(Result.SemiMajorRadius, Result.LambdaR)
    pyplot.plot(Result.SemiMajorRadius, Result.VS)
    pyplot.show()

    # Write the output file(s)

    # Finish
    print('Elapsed time: {0} seconds'.format(clock() - t))

# dap/trunk/python/mangadap/contrib/dapcontrib_database.py 8143 12701 'CUBE' 'STON' 1
# dapf = dapfile(plate=8143, ifudesign=12701, mode='CUBE', bintype='STON', niter=1, analysis_path='/home/misc/ManGA/Data/manga/spectro/analysis/MPL-3/')
