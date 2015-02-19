from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import sys
if sys.version > '3':
    long = int

import os.path
from os import environ
from astropy.io import fits
from astropy.wcs import WCS
import time
import numpy

from mangadap.util.parser import arginp_to_list
from mangadap.util.exception_tools import print_frame

__author__ = 'Kyle Westfall'

def parse_drp_file_name(name):
    """
    Parse the name of a DRP file and return the plate, ifudesign, and
    mode.
    """

    if (type(name) != str):
        raise Exception("File name must be a string.")
    if (str.find(name, 'manga-') == -1 or str.find(name, '-LOG') == -1 or 
        str.find(name, '.fits.gz') == -1):
        raise Exception("String does not look like a DRP fits file name.")

    plate_start = str.find(name, '-')+1
    plate_end = str.find(name,'-',plate_start)
    plate = long(name[plate_start:plate_end])

    ifudesign_end = str.find(name,'-',plate_end+1)
    ifudesign = long(name[plate_end+1:ifudesign_end])

    mode_start = str.find(name,'LOG')+3
    mode = name[mode_start:str.find(name,'.fits.gz')]

    return plate, ifudesign, mode

# Moved to mangadap.util
#def arginp_to_list(inp, evaluate=False, quiet=True):
#    """
#    Convert an input argument to a list of separated values.
#
#    If evaluate is True, evaluate each element in the list before
#    returning it.
#    """
#
#    out = inp
#
#    # Simply return a None value
#    if out is None:
#        return out
#
#    # If the input value is a string, convert it to a list of strings
#    if type(out) == str:
#        out = out.replace("[", "")
#        out = out.replace("]", "")
#        out = out.replace(" ", "")
#        out = out.split(',')
#
#    # If the input is still not a list, make it one
#    if type(out) != list:
#        out = [out]
#
#    # Evaluate each string, if requested
#    if evaluate:
#        n = len(out)
#        for i in range(0,n):
#            try:
#                tmp = eval(out[i])
#            except (NameError, TypeError):
#                if not quiet:
#                    print_frame('NameError/TypeError')
#                    print('Could not evaluate value {0}. Skipping.'.format(out[i]))
#            else:
#                out[i] = tmp
#
#    return out


def drpfile_list(drpver, platelist, ifudesignlist, modelist, combinatorics=False):
    """
    Provided a list of plates, ifu designs, and modes, return a list of
    DRP files to be analyzed by the DAP.

    ARGUMENTS:
        drpver              - The DRP version, which MUST be a single
                              string value used for all DRP files.
        platelist           - List of plates to use.
        ifudesignlist       - List of IFU designs to use.
        modelist            - List of DRP output modes ('CUBE' or 'RSS')
        combinatorics       - Based on the input, create a list with all
                              possible combinations.

    If the number of elements in each list is the same, the matching is
    assumed to be finished unless combinatorics is True.  If the number
    of elements is not the same, or cominatorics is True, the matched
    list is all the combinations of the provided elements.

    REVISION HISTORY:
        20 Nov 2014: Original implementation by K. Westfall
        12 Feb 2014: (KBW) Added directory_path()
    """

    if (type(drpver) != str):
        raise Exception('DRP drpver must be a string')

    if platelist is None or ifudesignlist is None or modelist is None:
        raise ValueError('Must provide platelist, ifudesignlist, and modelist!')

    platelist_ = arginp_to_list(platelist, evaluate=True)
    ifudesignlist_ = arginp_to_list(ifudesignlist, evaluate=True)
    modelist_ = arginp_to_list(modelist)

    n_plates = len(platelist_)
    n_designs = len(ifudesignlist_)
    n_modes = len(modelist_)

    # Perform the combinatorics
    if (n_plates != n_designs or n_plates != n_modes or combinatorics is True):

        # Force elements in platelist, ifudesignlist, and modelist
        # to be unique before performing combinatorics
        platelist_ = list(set(platelist_))
        ifudesignlist_ = list(set(ifudesignlist_))
        modelist_ = list(set(modelist_))

        n_plates = len(platelist_)
        n_designs = len(ifudesignlist_)
        n_modes = len(modelist_)

        # Each element of plate list is repeated n_designs*n_modes times
        for i in range(0,n_plates):
            for j in range (1,n_designs*n_modes):
                platelist_.insert(i*n_designs*n_modes, platelist_[i*n_designs*n_modes])

        # First repeat each element of ifudesign list n_modes times
        for i in range(0,n_designs):
            for j in range (1,n_modes):
                ifudesignlist_.insert(i*n_modes, ifudesignlist_[i*n_modes])
        # Then repeat the entire list n_plates times
        ifudesignlist__ = ifudesignlist_[:]
        for i in range(1,n_plates):
            ifudesignlist_ += ifudesignlist__

        # Mode list iterated most quickly, so repeat the full list
        # n_plates*n_designs times
        modelist__ = modelist_[:]
        for i in range(1,n_plates*n_designs):
            modelist_ += modelist__

        nn = n_plates*n_designs*n_modes
    else:
        nn = n_plates

    # Create and return the list of DRP files
    return [drpfile(drpver, platelist_[i], ifudesignlist_[i], modelist_[i]) for i in range(0,nn)]



class drpfile:
    """
    Object to hold the properties of a DRP-produced file to be analyzed by the DAP.

    REVISION HISTORY:
        20 Nov 2014: (KBW) Original Implementation
    """


    def __init__(self, drpver, plate, ifudesign, mode, pixelscale=None):
        """
        ARGUMENTS:
            drpver - DRP version to use
            plate - plate designation
            ifudesign - ifu design on the plate
            mode - mode of the produced file (CUBE or RSS)
        """

        # Set the attributes, forcing a known type
        self.drpver = str(drpver)
        self.plate = long(plate)
        self.ifudesign = long(ifudesign)
        self.mode = str(mode)

        # Set the image dimensions, needed for RSS files
        self.pixelscale = pixelscale if pixelscale is not None else self._default_pixel_scale()
        self.xs = None
        self.nx = None
        self.ys = None
        self.ny = None

#        self.manga_id = None
#        self.object_ra = None
#        self.object_dec = None

        self.hdulist = None             # Do not automatically read the data
        self.w = None                   # WCS structure

    def __del__(self):
        """
        Destroy the dapfile object, ensuring that the fits file is
        properly closed.
        """

        if self.hdulist is not None:
            self.hdulist.close()


    def _open_hdulist(self):
        """
        Internal support routine used to gather header information from
        the file.
        """
        if self.hdulist is not None:
            return

        inp = self.file_path()
        if (not os.path.exists(inp)):
            raise Exception('Cannot open file: {0}'.format(inp))

        self.hdulist = fits.open(inp)


    def _default_pixel_scale(self):
        """Return the default pixel scale in arcsec."""
        return 0.5


    def _cube_dimensions(self):
        self._open_hdulist()

        # Get the size in each dimension
        minx = min(self.hdulist['XPOS'].data)
        maxx = max(self.hdulist['XPOS'].data)
        Dx = maxx-minx

        miny = min(self.hdulist['YPOS'].data)
        maxy = max(self.hdulist['YPOS'].data)
        Dy = maxy-miny

        # Force the size to be the same in both dimensions
        Dx = Dx if Dx > Dy else Dy
        self.nx = int(Dx/self.pixelscale)+10
        if self.nx % 2 != 0:
            self.nx += 1
        self.ny = self.nx

        # Set the starting coordinate
        self.xs = -self.nx*self.pixelscale/2. + (minx+maxx)/2.0
        self.ys = -self.ny*self.pixelscale/2. + (miny+maxy)/2.0


    def _fix_header(self):
        self._open_hdulist()
        self.hdulist['FLUX'].header['CUNIT1'] = 'deg'
        self.hdulist['FLUX'].header['CUNIT2'] = 'deg'


    def file_name(self):
        """Return the name of the DRP file"""
        return 'manga-{0}-{1}-LOG{2}.fits.gz'.format(self.plate, self.ifudesign, self.mode)


    def directory_path(self):
        """Return the path to the directory with the DRP file."""
        return os.path.join(environ['MANGA_SPECTRO_REDUX'], self.drpver, str(self.plate), 'stack')


    def file_path(self):
        """Return the full path to the DRP file"""
        return os.path.join(self.directory_path(), self.file_name())


    def finding_chart_path(self):
        """
        Return the full path to the PNG finding chart for the targetted
        object.
        """
        return os.path.join(self.directory_path(), 'images', str(self.ifudesign)+'.png')


    def object_data(self):
        """Return some selected object data from the fits file"""
        self._open_hdulist()
        return self.hdulist[0].header['MANGAID'], self.hdulist[0].header['OBJRA'], \
               self.hdulist[0].header['OBJDEC']


    def object_coo(self):
        """Return the object coordinates read from the fits header"""
        self._open_hdulist()
        return self.hdulist[0].header['OBJRA'], self.hdulist[0].header['OBJDEC']


    def created_today(self):
        """Return True if the file was created today."""

        # Get the current time
        current_time = time.localtime()

        # Last time the DRP file was created (or modified?)
        created_time = time.gmtime(os.path.getctime(self.file_path()))

        return (current_time.tm_year == created_time.tm_year and
                current_time.tm_mon == created_time.tm_mon and 
                current_time.tm_mday == created_time.tm_mday)


    def pix_mesh(self):
        if self.mode is 'RSS':
            self._cube_dimensions()
        else:
            self._open_hdulist()
            self.nx = self.hdulist['FLUX'].header['NAXIS1']
            self.ny = self.hdulist['FLUX'].header['NAXIS2']
        
        return numpy.meshgrid(numpy.arange(0.0,self.nx)+1, numpy.arange(0.0,self.ny)+1)


    def world_mesh(self):
        if self.w is None:
            self._fix_header()
            self.w = WCS(header=self.hdulist['FLUX'].header,fix=False,naxis=(1,2))

        x,y = self.pix_mesh()
        xy = numpy.array([x.reshape(self.nx*self.ny),y.reshape(self.nx*self.ny)]).transpose()
        XY = self.w.all_pix2world(xy, 1)
        return XY[:,0].reshape(self.nx, self.ny), XY[:,1].reshape(self.nx, self.ny)


    def gri_composite(self):
        """
        Return the X and Y coordinate grid, and flux in the
        reconstructed gri image data (in nanomaggies).  The shape of the
        Z array is (NX,NY,3), with the g, r, and i image data in Z[:,:,0],
        Z[:,:,1], and Z[:,:,2], respectively.

        The reconstructed imaging data are not part of the RSS files, so
        this function returns 3 'None's if this drpfile is an RSS file.
        """

        if self.mode is 'RSS':
            return None, None, None

        X, Y = self.world_mesh() 

        Z = numpy.array([ self.hdulist['GIMG'].data, self.hdulist['RIMG'].data, \
              self.hdulist['IIMG'].data ] )

        return X, Y, Z


    def white_light(self, mask_list=None):
        """
        Return the X and Y coordinate grid, and a "white-light" flux density (in
        10^{-17} erg/s/cm^2/angstrom) integrated over the full spectral range,
        ignoring masked pixels:

                                \int f dl
                           wl = ---------    .
                                 \int dl

        The pixels to be masked are those that are flagged with any of
        the maskbits provided (mask_list).

        For RSS files, this takes into account the changing XPOS and
        YPOS with wavelength.  For CUBE files, this is a more
        straight-forward integration of the cube.

        ARGUMENTS:
            - mask_list: numpy string array with the list of
              'MANGA_DRPPIXFLAG' values to ignore in the integration

              If None, ignore any non-zero pixel in the MASK extension.

        """

        if mask_list is not None:
            print("Cannot yet interpret mask_list. Ignoring.")

        self.dlogl = 0.0001

#   TODO: Add log sampling to header of RSS files        

        mask = numpy.zeros(self.hdulist['MASK'].data.shape, dtype=numpy.float64)
        mask[ self.hdulist['MASK'].data.shape > 0. ] = 1.0

        if mode is 'RSS':
            X, Y = self.pix_mesh()
            X = (X-1)*self.pixelscale+self.xs
            X = (Y-1)*self.pixelscale+self.ys
            ipos = (self.hdulist['XPOS'].data - self.xs)/self.pixelscale
            jpos = (self.hdulist['YPOS'].data - self.ys)/self.pixelscale

            Z = numpy.zeros(X.shape, dtype=numpy.float64)

            wgt = numpy.zeros(X.shape, dtype=numpy.float64)
            nwave = (self.hdulist['FLUX'].header['NAXIS1'])
            nspec = (self.hdulist['FLUX'].header['NAXIS2'])

            for i in range(0,nspec):
                for j in range(0,nwave):
                    tmp = mask[i,j] * dlogl * self.hdulist['WAVE'].data[i,j]
                    Z[ipos[i,j], jpos[i,j]] += tmp * self.hdulist['FLUX'].data[i,j]
                    wgt[ipos[i,j], jpos[i,j]] += tmp

            Z[wgt > 0.] /= wgt[wgt > 0.]

            return X, Y, Z

        else:
            X, Y = self.world_mesh()

            nwave = (self.hdulist['FLUX'].header['NAXIS3'])
            Z = self.hdulist['FLUX'].data
            for i in range(0,nwave):
                mask[:,:,i]

            ipos = numpy.zeros(self.hdulist['FLUX'].data.shape, dtype=numpy.int32)
            jpos = numpy.zeros(self.hdulist['FLUX'].data.shape, dtype=numpy.int32)
            ipos[:,:,0], jpos[:,:,0] = numpy.meshgrid( numpy.arange(0, self.nx), \
                                                       numpy.arange(0,self.ny), indexing='ij')
            ipos[:,:,1:-1] = ipos[:,:,0]
            
        Z = numpy.zeros(X.shape)
            
   
        # x_grid will ensure that the hdulist is read
        x = self.x_grid()
        y = self.y_grid()
        X, Y = numpy.meshgrid(x, y)



        

    
        




