from __future__ import division
from __future__ import print_function

import os.path
from os import environ
from astropy.io import fits
import time
import numpy
from drpfile import arginp_to_list

__author__ = 'Kyle Westfall'

def parse_dap_file_name(name):
    """
    Parse the name of a DAP file and return the plate, ifudesign, mode,
    binning type, and iteration number.
    """

    if (type(name) != str):
        raise Exception("File name must be a string.")
    if (str.find(name, 'manga-') == -1 or str.find(name, '-LOG') == -1
        or str.find(name, 'BIN-') == -1 or str.find(name, '.fits') == -1):
        raise Exception("String does not look like a DAP fits file name.")

    plate_start = str.find(name, '-')+1
    plate_end = str.find(name,'-',plate_start)
    plate = int(name[plate_start:plate_end])

    ifudesign_end = str.find(name,'-',plate_end+1)
    ifudesign = int(name[plate_end+1:ifudesign_end])

    mode_start = str.find(name,'LOG')+3
    mode_end = str.find(name,'_BIN-')
    mode = name[mode_start:mode_end]

    bintype_end = str.find(name,'-',mode_end+5)
    bintype = name[mode_end+5:bintype_end]

    niter_end = str.find(name,'.fits',bintype_end)
    niter = int(name[bintype_end:niter_end])

    return plate, ifudesign, mode, bintype, niter


class dapfile:
    """
    Object used to hold properties of and read data from a DAP-produced file.

    REVISION HISTORY:
        16 Dec 2014: (KBW) Original Implementation
    """

    # TODO: Allow default of niter=None, then search for the first file
    # with the appropriate name and set niter

    def __init__(self, dapver, plate, ifudesign, mode, bintype, niter, read=True):
        """
        ARGUMENTS:
            dapver - DAP version to use
            plate - plate designation
            ifudesign - ifu design on the plate
            mode - mode of the produced file (CUBE or RSS)
            bintype - binning type used
            niter - iteration number
        """

        # Set the attributes, forcing a known type
        self.dapver = str(dapver)
        self.plate = int(plate)
        self.ifudesign = int(ifudesign)
        self.mode = str(mode)
        self.bintype = str(bintype)
        self.niter = int(niter)

        self.hdulist = None

        if read:
            self.read()

    def __del__(self):
        """
        Destroy the dapfile object, ensuring that the fits file is
        properly closed.
        """

        if self.hdulist is not None:
            self.hdulist.close()

    
    def _twod_image_data(self, exten, indx):
        try:
            data = self.hdulist[exten].data[:,indx]
        except IndexError as e:
            print('{0}'.format(e))
            return None
        else:
            return data
    


#    def _gather_data(self):
#        """
#        Internal support routine used to gather header information from
#        the file.
#        """
#        if (self.manga_id is not None and self.object_ra is not None and
#            self.object_dec is not None):
#            return
#
#        inp = self.file_path()
#        if (not os.path.exists(inp)):
#            raise Exception('Cannot open file: {0}'.format(inp))
#
#        # TODO: This is very slow!  memmap in fits.open() does not help.
#        # A way to read header without setting up the whole hdulist
#        # object?  hdulist = fits.open(inp,memmap=True)
#        hdr = fits.getheader(inp,0)
#        self.manga_id = hdr['MANGAID']
#        self.object_ra = hdr['OBJRA']
#        self.object_dec = hdr['OBJDEC']
#
##        hdulist = fits.open(inp)
##        self.manga_id = hdulist[0].header['MANGAID']
##        self.object_ra = hdulist[0].header['OBJRA']
##        self.object_dec = hdulist[0].header['OBJDEC']
##        hdulist.close()


    def file_name(self):
        """Return the name of the DAP file"""

        # Number of spaces provided for iteration number is 3
        siter = str(self.niter).zfill(3)
        return 'manga-{0}-{1}-LOG{2}_BIN-{3}-{4}.fits'.format(self.plate, self.ifudesign,
                                                              self.mode, self.bintype, siter)


    def file_path(self):
        """Return the full path to the DAP file"""
        return os.path.join(environ['MANGA_SPECTRO_ANALYSIS'], self.dapver, str(self.plate),
                            str(self.ifudesign), self.file_name())


    def read(self):
        """
        Read all the fits file data and keep it in memory
        """

        inp = self.file_path()
        if (not os.path.exists(inp)):
            raise Exception('Cannot open file: {0}'.format(inp))

        if self.hdulist is not None:
            self.hdulist.close()

        self.hdulist = fits.open(inp)

    def close(self):
        self.hdulist.close()
        self.hdulist = None


    def header(self):
        """Return the primary header"""
        if self.hdulist is None:
            self.read()

        return self.hdulist[0].header


    def sn_stats(self):
        if self.hdulist is None:
            self.read()
        return self.hdulist[0].header['SNWAVE1'], self.hdulist[0].header['SNWAVE2'], \
               self.hdulist[0].header['MINSNBIN']


    def spaxel_size(self):
        if self.hdulist is None:
            self.read()
        return self.hdulist[0].header['SPAXDX'], self.hdulist[0].header['SPAXDY']


    def bin_ston(self):
        if self.hdulist is None:
            self.read()
        if self.bintype is not 'STON':
            raise Exception('Binning type was not STON!')
        return self.hdulist[0].header['BINSN']


    def analysis_stats(self):
        if self.hdulist is None:
            self.read()
        return self.hdulist[0].header['FITWAVE1'], self.hdulist[0].header['FITWAVE2'], \
               self.hdulist[0].header['MINSNFIT']


    def table_column(self, exten, col):
        if self.hdulist is None:
            self.read()
        try:
            data = self.hdulist[exten].data[col]
        except KeyError as e:
            print('Could not find extension ({0}) and/or column ({1})!'.format(exten, col))
            return None
        else:
            return data

    
    def bin_wavelength(self):
        return self.hdulist['WAVE'].data


    def bin_spectrum(self, indx):
        return self._twod_image_data('FLUX', indx)


    def bin_inverse_variance(self, indx):
        return self._twod_image_data('IVAR', indx)
            

    def bin_mask(self, indx):
        return self._twod_image_data('MASK', indx)


    def fit_stars_minus_gas(self, indx):
        bestfit = self._twod_image_data('SGMOD', indx)
        emlmodel = self._twod_image_data('ELMOD', indx)
        if bestfit is not None and emlmodel is not None:
            return (bestfit - emlmodel)
        else:
            return None

            
    def fit_stars_and_gas(self, indx):
        return self._twod_image_data('SGMOD', indx)

    
    def fit_eml_only_ew(self, indx):
        return self._twod_image_data('ELOMEW', indx)
            

    def fit_eml_only_fb(self, indx):
        return self._twod_image_data('ELOMFB', indx)
            


