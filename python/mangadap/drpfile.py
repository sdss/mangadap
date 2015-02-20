from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import sys
if sys.version > '3':
    long = int

import os.path
from os import environ
from scipy import sparse
from astropy.io import fits
from astropy.wcs import WCS
import time
import numpy

from mangadap.util.parser import arginp_to_list

__author__ = 'Kyle B. Westfall'

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


    def __init__(self, plate, ifudesign, mode, drpver=None, directory_path=None, pixelscale=None,
                 read=False):
        """
        ARGUMENTS:
            plate integer
                    plate designation

            ifudesign integer
            
                    ifu design on the plate

            mode string
                    mode of the produced file (CUBE or RSS)

        OPTIONAL:
            drpver string
                    DRP version to use.  Default is to use the
                    $MANGADAP_VER environmental variable

            directory_path string
                    Directory with the DRP file.  Default is to use:
                    
                    os.path.join(environ['MANGA_SPECTRO_REDUX'],
                                 self.drpver, str(self.plate), 'stack')

            pixelscale double
                    size of the pixel in arcsec.  Default is provided by _default_pixel_scale

        """

        # Set the attributes, forcing a known type
        self.plate = long(plate)
        self.ifudesign = long(ifudesign)
        self.mode = str(mode)

        if (mode is not 'CUBE') & (mode is not 'RSS'):
            raise Exception('{0} is not a viable mode.  Must be RSS or CUBE.'.format(mode))

        self.drpver = self._default_drp_version() if drpver is None else str(drpver)
        if directory_path is None:
            self.directory_path = self._default_directory_path()
        else:
            self.directory_path = str(directory_path)

        # Set the image dimensions, needed for RSS files
        self.pixelscale = pixelscale if pixelscale is not None else self._default_pixel_scale()
        self.xs = None
        self.nx = None
        self.ys = None
        self.ny = None

        self.hdu = None                 # Do not automatically read the data
        self.w = None                   # WCS structure

        # Provide internal variables to keep the currently available
        # regridding transfer matrix for computing a plane of the data
        # cube, the index of the wavelength channel defining that plane,
        # and the properties used to create the matrix.
        self.regrid_T = None            
        self.regrid_channel = None
        self.regrid_rlim = None
        self.regrid_sigma = None

        if read:
            self._open_hdu()
        

    def __del__(self):
        """
        Destroy the drpfile object by ensuring that the fits file is
        properly closed.
        """
        if self.hdu is not None:
            self.hdu.close()
            self.hdu = None


    def _open_hdu(self):
        """
        Internal support routine used to gather header information from
        the file.
        """
        if self.hdu is not None:
            return

        inp = self.file_path()
        if (not os.path.exists(inp)):
            raise Exception('Cannot open file: {0}'.format(inp))

        # Open the fits file, but do NOT allow the file to be
        # overwritten
        self.hdu = fits.open(inp, mode='readonly')


    def _default_pixel_scale(self):
        """Return the default pixel scale in arcsec."""
        return 0.5


    def _default_drp_version(self):
        """
        Return the DRP version defined by the environmental variable.
        """
        return environ['MANGADRP_VER']


    def _default_directory_path(self):
        """Return the directory path used by the DRP."""

        # Make sure the DRP version is set
        if self.drpver is None:
            self.drpver = self._default_drp_version()

        return os.path.join(environ['MANGA_SPECTRO_REDUX'], self.drpver, str(self.plate), 'stack')

    def _dimensions_undefined(self):
        if self.xs is None:
            return True
        if self.nx is None:
            return True
        if self.ys is None:
            return True
        if self.ny is None:
            return True
        return False


    def _regrid_transfer_undefined(self):
        if self.regrid_T is None:
            return True
        if self.regrid_channel is None:
            return True
        if self.regrid_rlim is None:
            return True
        if self.regrid_sigma is None:
            return True
        return False


    def _regrid_T_correct(self, channel, rlim, sigma):
        if self.regrid_channel != channel:
            return False
        if self.regrid_rlim != rlim:
            return False
        if self.regrid_sigma != sigma:
            return False
        return True


    def _cube_dimensions(self, recenter, width_buffer, redo=False):
        """
        Determine the on-sky dimensions of a CUBE file that will include
        all the RSS spectra using the data in the 'XPOS' and 'YPOS'
        extensions.  This is only used if mode='RSS' (TODO: true?).
        
        If recenter is false, the coordinates the XPOS,YPOS coordinates
        are assumed to be centered at 0,0.  If recenter is true, the
        XPOS,YPOS can have any center, and the center of the cube is set
        to be ~ the center of the range in XPOS,YPOS.

        width_buffer is the number of pixels to add to the width of the
        cube in addition to the range needed to cover XPOS and YPOS.
        """

        # Only continue if one of the defining dimensions is None,
        # unless the cube_dimensions should be redetermined (redo=True)
        if not redo and not self._dimensions_undefined():
            return

        # Make sure that the fits file is ready for reading
        self._open_hdu()

        # Get the size in each dimension
        minx = numpy.amin(self.hdu['XPOS'].data)
        maxx = numpy.amax(self.hdu['XPOS'].data)
        Dx = maxx-minx

        miny = numpy.amin(self.hdu['YPOS'].data)
        maxy = numpy.amax(self.hdu['YPOS'].data)
        Dy = maxy-miny

        # Force the size to be even and the same in both dimensions
        Dx = Dx if Dx > Dy else Dy
        self.nx = int(Dx/self.pixelscale)+width_buffer
        if self.nx % 2 != 0:
            self.nx += 1
        self.ny = self.nx

        # Set the starting coordinate
        self.xs = -self.nx*self.pixelscale/2.
        self.ys = -self.ny*self.pixelscale/2.

        # Offset to the center, if requested
        if recenter:
            self.xs = self.xs + (minx+maxx)/2.0
            self.ys = self.ys + (miny+maxy)/2.0


    def _fix_header(self):
        self._open_hdu()

        # So that astropy 0.4.4 wcs can read the first two dimensions of
        # the WCS header values
        self.hdu['FLUX'].header['CUNIT1'] = 'deg'
        self.hdu['FLUX'].header['CUNIT2'] = 'deg'

        # TODO: Add the ability remove the SIMPLE keywords from the
        # fits extensions


    def file_name(self):
        """Return the name of the DRP file"""
        return 'manga-{0}-{1}-LOG{2}.fits.gz'.format(self.plate, self.ifudesign, self.mode)


#    def directory_path(self):
#        """Return the path to the directory with the DRP file."""
#        return os.path.join(environ['MANGA_SPECTRO_REDUX'], self.drpver, str(self.plate), 'stack')


    def file_path(self):
        """Return the full path to the DRP file"""
        return os.path.join(self.directory_path, self.file_name())


    def finding_chart_path(self):
        """
        Return the full path to the PNG finding chart for the targetted
        object.
        """
        return os.path.join(self.directory_path, 'images', str(self.ifudesign)+'.png')


    def object_data(self):
        """
        Return some selected object data from the fits file header.
        """
        self._open_hdu()
        return self.hdu[0].header['MANGAID'], self.hdu[0].header['OBJRA'], \
               self.hdu[0].header['OBJDEC']


    def object_coo(self):
        """Return the object coordinates read from the fits header"""
        self._open_hdu()
        return self.hdu[0].header['OBJRA'], self.hdu[0].header['OBJDEC']


    def created_today(self):
        """Return True if the file was created today."""

        # Get the current time
        current_time = time.localtime()

        # Last time the DRP file was created (or modified?)
        created_time = time.gmtime(os.path.getctime(self.file_path()))

        return (current_time.tm_year == created_time.tm_year and
                current_time.tm_mon == created_time.tm_mon and 
                current_time.tm_mday == created_time.tm_mday)


    def pix_mesh(self, recenter=False, width_buffer=10):
        """
        Return the I and J pixel coordinates for an nx*ny mesh.  For the
        'RSS' files, this is determined following the same procedure as
        the DRP.  For the 'CUBE' files, this is pulled directly from the
        header of the 'FLUX' extension.
        """
        if self.mode is 'RSS':
            self._cube_dimensions(recenter, width_buffer)
        else:
            self._open_hdu()
            self.nx = self.hdu['FLUX'].header['NAXIS1']
            self.ny = self.hdu['FLUX'].header['NAXIS2']
        
        return numpy.meshgrid(numpy.arange(0.0,self.nx)+1, numpy.arange(0.0,self.ny)+1)


    def world_mesh(self):
        """
        Return the world X and Y coordinate for the self.nx*self.ny
        mesh.  This can only be done for the 'CUBE' files because it
        requires WCS parameters in the header of the 'FLUX' extension.
        """
        if self.mode is 'RSS':
            raise Exception('Cannot determine world coordinates for RSS files.')

        if self.w is None:
            self._fix_header()
            self.w = WCS(header=self.hdu['FLUX'].header,fix=False,naxis=(1,2))

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

        Z = numpy.transpose(numpy.array([ self.hdu['GIMG'].data, self.hdu['RIMG'].data, \
                            self.hdu['IIMG'].data ] ), axes=(1,2,0))

        return X, Y, Z


    def regrid_transfer_matrix(self, channel, rlim=1.6, sigma=0.7, quiet=False):
        """
        Calculate the sparse matrix used to regrid the RSS data into a
        data cube.  This can only be done using RSS files.

        ARGUMENTS:
                channel integer
                        Index of the spectral channel for which to
                        calculate the transfer matrix.

                rlim double
                        Limiting radius of the Gaussian interpolate in
                        arcsec for each fiber.  Default is 1.6 arcsec.

                sigma double
                        Sigma of the interpolating Gaussian in arcsec.
                        Default is 0.6 arcsec.

        RETURNS:
                The transfer matrix with type scipy.sparse.csr_matrix.
        """

        # Check if the calculation is necessary
        if not self._regrid_transfer_undefined() and self._regrid_T_correct(channel, rlim, sigma):
            return self.regrid_T

        # For CUBE files, try to use the RSS counterpart to get the
        # transfer matrix
        if self.mode is 'CUBE':
            print('Attempting to use RSS counter-part for calculation.')
            drpf = drpfile(self.plate, self.ifudesign, 'RSS', self.drpver, self.directory_path,
                           self.pixelscale)
            self.regrid_T = drpf.regrid_transfer_matrix(channel, rlim, sigma)
            self.regrid_channel = drpf.regrid_channel
            self.regrid_rlim = drpf.regrid_rlim
            self.regrid_sigma = drpf.regrid_sigma

            # Also need to make sure dimensions are set
            self.nx = self.hdu['FLUX'].header['NAXIS1']
            self.ny = self.hdu['FLUX'].header['NAXIS2']
        
            return self.regrid_T

        # Get the cube dimensions; making sure they match the
        # calculation used by the DRP; will open the hdu if it isn't
        # already
        self._cube_dimensions(recenter=False, width_buffer=10, redo=True)

        # Dimensions of the sparse matrix are
        ns = self.hdu['FLUX'].data.shape[0]     # The number of spectra
        nim = self.nx*self.ny                   # The number of image pixels

        # Get the list of non-zero pixel values in the transfer matrix
        i = numpy.arange(0,self.nx)
        j = numpy.arange(0,self.ny)
        ii,jj = numpy.meshgrid(i,j)                         # Mesh of i,j pixel indices

        sp = numpy.ndarray(nim)                             # Holds spectrum index
        ij = (ii*self.ny+jj).reshape(nim)                   # Holds image pixel index
        r2 = numpy.ndarray(nim)                             # Holds radii

        s2 = numpy.square(sigma/self.pixelscale)            # sigma^2 of Gaussian
        rl2 = numpy.square(rlim/self.pixelscale)            # radius^2 of Gaussian limit

        non_zero_spc = numpy.ndarray(0, dtype=numpy.int64)  # Holds triplet spectrum index
        non_zero_pix = numpy.ndarray(0, dtype=numpy.int64)  # Holds triplet image index
        non_zero_wgt = numpy.ndarray(0, dtype=numpy.float64)# Holds triplet weight

        # TODO: Can optimize this further?
        for k in range(0,ns):

            # Fill spectrum index
            sp.fill(k)

            # Calculate radii
            r2[:] = numpy.square( (self.hdu['XPOS'].data[k,channel]-self.xs)/self.pixelscale \
                                  - ii.reshape(nim)) \
                    + numpy.square((self.hdu['YPOS'].data[k,channel]-self.ys)/self.pixelscale \
                                  - jj.reshape(nim))
           
            # Append new indices and weights within rlim
            non_zero_spc = numpy.append(non_zero_spc, sp[r2 < rl2])
            non_zero_pix = numpy.append(non_zero_pix, ij[r2 < rl2])
            wgt = numpy.exp(-r2[r2 < rl2]/s2/2.0)
            wgt /= numpy.sum(wgt)
            non_zero_wgt = numpy.append(non_zero_wgt, wgt)

            if not quiet:
                print('Transfer Matrix {:2.1%}'.format((k+1)/ns), end="\r")

        if not quiet:
            print('Transfer Matrix Done                     ')
        # Save the regridding input and result
        self.regrid_rlim = rlim
        self.regrid_sigma = sigma
        self.regrid_channel = channel
        self.regrid_T = sparse.coo_matrix( (non_zero_wgt, (non_zero_pix, non_zero_spc)), \
                                           shape=(nim,ns) ).tocsr()

        return self.regrid_T


    def regrid_wavelength_plane(self, channel, rlim=1.6, sigma=0.7, quiet=False):
        """
        Return the regridded wavelength plane.  If this is a CUBE file,
        it first checks that the requested rlim and sigma are the same
        as the default DRP values.  If not, it tries to use the RSS
        counterpart to create the image.  CUBE files should not be used
        in general because calculation of each chennel requires
        re-reading the RSS file EACH TIME.
        """
        if self.mode is 'CUBE':
            if rlim == 1.6 and sigma == 0.7:
                self._open_hdu()
                return self.hdu['FLUX'].data[channel,:,:]
            print('Attempting to use RSS counter-part for calculation.')
            drpf = drpfile(self.plate, self.ifudesign, 'RSS', self.drpver, self.directory_path,
                           self.pixelscale)
            return drpf.regrid_wavelength_plane(channel, rlim, sigma, quiet)

        # Set the transfer matrix (set to self.regrid_T; don't need to
        # keep the returned matrix)
        self.regrid_transfer_matrix(channel, rlim, sigma, quiet)

        return self.regrid_T.dot(self.hdu['FLUX'].data[:,channel]).reshape(self.nx,self.ny)


    def covariance_matrix(self, channel, rlim=1.6, sigma=0.7, quiet=False):
        """
        Return the covariance matrix for a specified wavelength channel.
        If this is a CUBE file, it tries to use the RSS counterpart to
        create the covariance matrix.  CUBE files should not be used in
        general because calculation of each chennel requires re-reading
        the RSS file EACH TIME.
        """

        if self.mode is 'CUBE':
            print('Attempting to use RSS counter-part for calculation.')
            drpf = drpfile(self.plate, self.ifudesign, 'RSS', self.drpver, self.directory_path,
                           self.pixelscale)
            return drpf.covariance_matrix(channel, rlim, sigma, quiet)

        # Set the transfer matrix (set to self.regrid_T; don't need to
        # keep the returned matrix)
        self.regrid_transfer_matrix(channel, rlim, sigma, quiet)

        ns = self.hdu['FLUX'].data.shape[0]     # The number of spectra
        var = numpy.zeros(ns, dtype=numpy.float64)
        indx = numpy.where(self.hdu['IVAR'].data[:,channel] > 0.)
        var[indx] = 1.0/self.hdu['IVAR'].data[indx,channel]

        # Set the IVAR values along the diagonal of a sparse matrix
        ns = self.hdu['FLUX'].data.shape[0]     # The number of spectra
        Sigma = sparse.coo_matrix( (var, (numpy.arange(0,ns),numpy.arange(0,ns))), \
                                   shape=(ns,ns)).tocsr()

        return self.regrid_T.dot(Sigma.dot(self.regrid_T.transpose()))


#    def covariance_cube(self, rlim=1.6, sigma=0.7, quiet=False):
#
#        if self.mode is 'CUBE':
#            print('Attempting to use RSS counter-part for calculation.')
#            drpf = drpfile(self.plate, self.ifudesign, 'RSS', self.drpver, self.directory_path,
#                           self.pixelscale)
#            return drpf.covariance_cube(rlim, sigma, quiet)
#
#        ns = self.hdu['FLUX'].data.shape[0]     # The number of spectra
#        nc = self.hdu['FLUX'].data.shape[1]     # The number of wavelength channels
#
#        non_zero_pix = numpy.ndarray(0, dtype=numpy.int64)  # Holds triplet image index
#        non_zero_chn = numpy.ndarray(0, dtype=numpy.int64)  # Holds triplet image*channel index
#        non_zero_cov = numpy.ndarray(0, dtype=numpy.float64)# Holds triplet covariance value
#
#        # The covariance "cube" will be a sparse matrix that has
#        # dimensions (nx*ny, nx*ny*nc)
#        for i in range(0,nc):
#            C = self.covariance_matrix(i,rlim,sigma,quiet)
#
#            non_zero_pix = numpy.append(non_zero_pix, (sparse.find(C)
#
#
#
#        # Set the transfer matrix (set to self.regrid_T; don't need to
#        # keep the returned matrix)
#        self.regrid_transfer_matrix(channel, rlim, sigma, quiet)
#
#        var = numpy.zeros(ns, dtype=numpy.float64)
#        indx = numpy.where(self.hdu['IVAR'].data[:,channel] > 0.)
#        var[indx] = 1.0/self.hdu['IVAR'].data[indx,channel]
#
#        # Set the IVAR values along the diagonal of a sparse matrix
#        ns = self.hdu['FLUX'].data.shape[0]     # The number of spectra
#        Sigma = sparse.coo_matrix( (var, (numpy.arange(0,ns),numpy.arange(0,ns))), \
#                                   shape=(ns,ns)).tocsr()
#
#        return self.regrid_T.dot(Sigma.dot(self.regrid_T.transpose()))



#   def white_light(self, mask_list=None):
#       """
#       Return the X and Y coordinate grid, and a "white-light" flux density (in
#       10^{-17} erg/s/cm^2/angstrom) integrated over the full spectral range,
#       ignoring masked pixels:

#                               \int f dl
#                          wl = ---------    .
#                                \int dl

#       The pixels to be masked are those that are flagged with any of
#       the maskbits provided (mask_list).

#       For RSS files, this takes into account the changing XPOS and
#       YPOS with wavelength.  For CUBE files, this is a more
#       straight-forward integration of the cube.

#       ARGUMENTS:
#           - mask_list: numpy string array with the list of
#             'MANGA_DRPPIXFLAG' values to ignore in the integration

#             If None, ignore any non-zero pixel in the MASK extension.

#       """

#       if mask_list is not None:
#           print("Cannot yet interpret mask_list. Ignoring.")

#       self.dlogl = 0.0001

#   TODO: Add log sampling to header of RSS files        

#       mask = numpy.zeros(self.hdu['MASK'].data.shape, dtype=numpy.float64)
#       mask[ self.hdu['MASK'].data.shape > 0. ] = 1.0

#       if mode is 'RSS':
#           X, Y = self.pix_mesh()
#           X = (X-1)*self.pixelscale+self.xs
#           X = (Y-1)*self.pixelscale+self.ys
#           ipos = (self.hdu['XPOS'].data - self.xs)/self.pixelscale
#           jpos = (self.hdu['YPOS'].data - self.ys)/self.pixelscale

#           Z = numpy.zeros(X.shape, dtype=numpy.float64)

#           wgt = numpy.zeros(X.shape, dtype=numpy.float64)
#           nwave = (self.hdu['FLUX'].header['NAXIS1'])
#           nspec = (self.hdu['FLUX'].header['NAXIS2'])

#           for i in range(0,nspec):
#               for j in range(0,nwave):
#                   tmp = mask[i,j] * dlogl * self.hdu['WAVE'].data[i,j]
#                   Z[ipos[i,j], jpos[i,j]] += tmp * self.hdu['FLUX'].data[i,j]
#                   wgt[ipos[i,j], jpos[i,j]] += tmp

#           Z[wgt > 0.] /= wgt[wgt > 0.]

#           return X, Y, Z

#       else:
#           X, Y = self.world_mesh()

#           nwave = (self.hdu['FLUX'].header['NAXIS3'])
#           Z = self.hdu['FLUX'].data
#           for i in range(0,nwave):
#               mask[:,:,i]

#           ipos = numpy.zeros(self.hdu['FLUX'].data.shape, dtype=numpy.int32)
#           jpos = numpy.zeros(self.hdu['FLUX'].data.shape, dtype=numpy.int32)
#           ipos[:,:,0], jpos[:,:,0] = numpy.meshgrid( numpy.arange(0, self.nx), \
#                                                      numpy.arange(0,self.ny), indexing='ij')
#           ipos[:,:,1:-1] = ipos[:,:,0]
#           
#       Z = numpy.zeros(X.shape)
#           
#  
#       # x_grid will ensure that the hdu is read
#       x = self.x_grid()
#       y = self.y_grid()
#       X, Y = numpy.meshgrid(x, y)



        

    
        




