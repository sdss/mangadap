from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import sys
if sys.version > '3':
    long = int

import os.path
from os import environ
from astropy.io import fits
import time
import numpy

from mangadap.util.exception_tools import print_frame
from mangadap.util.yanny import yanny
from mangadap.util.projected_disk_plane import projected_disk_plane


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


def default_dap_version():
    """Return the DAP version defined by the environmental variable."""
    return environ['MANGADAP_VER']


def default_analysis_path(dapver):
    """Return the directory path used by the DAP."""
    # Make sure the DRP version is set
    if dapver is None:
        dapver = default_dap_version()
    return os.path.join(environ['MANGA_SPECTRO_ANALYSIS'], dapver)


def default_dap_directory_path(dapver, analysis_path, plate, ifudesign):
    """Return the directory path used by the DAP."""
    # Make sure the DRP version is set
    if analysis_path is None:
        analysis_path = default_analysis_path(dapver)
    return os.path.join(analysis_path, str(plate), str(ifudesign))


def default_dap_par_file(dapver, analysis_path, directory_path, plate, ifudesign, mode):
    if directory_path is None:
        directory_path = default_dap_directory_path(dapver, analysis_path, plate, ifudesign)
    par_file = 'mangadap-{0}-{1}-LOG{2}.par'.format(plate, ifudesign, mode)
    return os.path.join(directory_path, par_file)

    
class dapfile:
    """
    Object used to hold properties of and read data from a DAP-produced file.

    REVISION HISTORY:
        16 Dec 2014: Original Implementation by K. Westfall (KBW)
        19 Mar 2015: (KBW) Added analysis_path
    """

    # TODO: Allow default of niter=None, then search for the first file
    # with the appropriate name and set niter

    def __init__(self, plate, ifudesign, mode, bintype, niter, dapver=None, analysis_path=None,
                 directory_path=None, par_file=None, read=True):
        """
        ARGUMENTS:
            plate - plate designation
            ifudesign - ifu design on the plate
            mode - mode of the produced file (CUBE or RSS)
            bintype - binning type used
            niter - iteration number

        OPTIONAL:
            dapver - DAP version to use; ONLY USED TO DEFINE THE
                     DIRECTORY PATH. Defaults to $MANGADAP_VER.

            analysis_path - Main output path for the DAP.  If not
                            provided, the default is:
                                os.path.join(environ['MANGA_SPECTRO_ANALYSIS'],
                                 self.dapver)

            directory_path - 
            
                             Exact path the DAP file.  If not provided,
                             the default is:
                               os.path.join(self.analysis_path, str(self.plate),
                                            str(self.ifudesign))

            par_file - SDSS parameter file used to provide input
                       parameters for the DAP.  If not provided, the
                       default is:
                            os.path.join(self.directory_path,
                              'mangadap-{0}-{1}-LOG{2}.par'.format(self.plate,
                              self.ifudesign, self.mode)
        """

        # Set the attributes, forcing a known type
        self.plate = int(plate)
        self.ifudesign = int(ifudesign)
        self.mode = str(mode)
        self.bintype = str(bintype)
        self.niter = int(niter)

        self.dapver = default_dap_version() if dapver is None else str(dapver)
        self.analysis_path = default_analysis_path() if analysis_path is None \
                                                           else str(analysis_path)
        self.directory_path = self._default_directory_path() if directory_path is None \
                                                             else str(directory_path)

        self.par_file = self._default_par_file() if par_file is None else str(par_file)

        self.hdu = None
        self.par = None

        if read:
            self._open_hdu()


    def __del__(self):
        """
        Destroy the dapfile object, ensuring that the fits file is
        properly closed.
        """
        if self.hdu is None:
            return
        self.hdu.close()
        self.hdu = None


    def _open_hdu(self, permissions='readonly'):
        """
        Internal support routine used to gather header information from
        the file.
        """
        if self.hdu is not None:
            return

        inp = self.file_path()
        if not os.path.exists(inp):
            raise Exception('Cannot open file: {0}'.format(inp))

        # Open the fits file with the requested read/write permission
        self.hdu = fits.open(inp, mode=permissions)

    
    def _read_par(self):
        if self.par is not None:
            return

        if not os.path.exists(self.par_file):
            raise Exception('Cannot open file: {0}'.format(inp))
            
        self.par = yanny(filename=self.par_file, np=True)


    def _default_dap_version(self):
        """
        Return the DAP version defined by the environmental variable.
        """
        return environ['MANGADAP_VER']


    def _default_analysis_path(self):
        """Return the directory path used by the DAP."""

        # Make sure the DRP version is set
        if self.dapver is None:
            self.dapver = self._default_dap_version()

        return os.path.join(environ['MANGA_SPECTRO_ANALYSIS'], self.dapver)


    def _default_directory_path(self):
        """Return the directory path used by the DAP."""

        # Make sure the DRP version is set
        if self.analysis_path is None:
            self.analysis_path = self._default_analysis_path()

        return os.path.join(self.analysis_path, str(self.plate), str(self.ifudesign))


    def _default_par_file(self):
        if self.directory_path is None:
            self.directory_path = self._default_directory_path()
        par_file = 'mangadap-{0}-{1}-LOG{2}.par'.format(self.plate, self.ifudesign, self.mode)
        return os.path.join(self.directory_path, par_file)

    
    def _twod_image_data(self, exten, indx):
        self._open_hdu()
        try:
            data = self.hdu[exten].data[:,indx]
        except IndexError as e:
            print_frame('IndexError')
            print('{0}'.format(e))
            return None
        else:
            return data
    
    
    def file_name(self):
        """Return the name of the DAP file"""

        # Number of spaces provided for iteration number is 3
        siter = str(self.niter).zfill(3)
        return 'manga-{0}-{1}-LOG{2}_BIN-{3}-{4}.fits'.format(self.plate, self.ifudesign,
                                                              self.mode, self.bintype, siter)


    def file_path(self):
        """Return the full path to the DAP file"""
        return os.path.join(self.directory_path, self.file_name())


    def header(self):
        """Return the primary header"""
        self._open_hdu()
        return self.hdu[0].header


    def sn_stats(self):
        self._open_hdu()
        return self.hdu[0].header['SNWAVE1'], self.hdu[0].header['SNWAVE2'], \
               self.hdu[0].header['MINSNBIN']


    def spaxel_size(self):
        self._open_hdu()
        return self.hdu[0].header['SPAXDX'], self.hdu[0].header['SPAXDY']


    def bin_ston(self):
        if self.bintype is not 'STON':
            raise Exception('Binning type was not STON!')
        self._open_hdu()
        return self.hdu[0].header['BINSN']


    def analysis_stats(self):
        self._open_hdu()
        return self.hdu[0].header['FITWAVE1'], self.hdu[0].header['FITWAVE2'], \
               self.hdu[0].header['MINSNFIT']


    def table_column(self, exten, col):
        self._open_hdu()
        try:
            data = self.hdu[exten].data[col]
        except KeyError as e:
            print_frame('KeyError')
            print('Could not find extension ({0}) and/or column ({1})!'.format(exten, col))
            return None
        else:
            return data

    
    def bin_wavelength(self):
        self._open_hdu()
        return self.hdu['WAVE'].data


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


    def bin_disk_polar_coo(self, xc=None, yc=None, rot=None, pa=None, inc=None):

        # Position angle and inclination
        if pa is None or inc is None:
            try:
                self._read_par()                # Try to read the parameter file
            except:
                print('WARNING: Using default pa and/or inclination.')

            if self.par is not None:
                if pa is None:
                    pa = self.guess_position_angle()
                if inc is None:
                    inc = self.guess_inclination()

        # Declare the galaxy coo object
        disk_plane = projected_disk_plane(xc, yc, rot, pa, inc)

        # Calculate the in-plane radius and azimuth for each bin
        self._open_hdu()
        nbin = self.hdu['BINS'].data['BINXRL'].size
        radius = numpy.zeros(nbin, dtype=numpy.float64)
        azimuth = numpy.zeros(nbin, dtype=numpy.float64)
        for i in range(0,nbin):
            radius[i], azimuth[i] = disk_plane.polar(self.hdu['BINS'].data['BINXRL'][i], \
                                                     self.hdu['BINS'].data['BINYRU'][i])
            
        return radius, azimuth

    
    def guess_cz(self):
        self._read_par()
        return self.par['DAPPAR']['vel'][0]


    def guess_inclination(self):
        self._read_par()
        return numpy.degrees( numpy.arccos(1.0 - self.par['DAPPAR']['ell'][0]) )


    def guess_position_angle(self):
        self._read_par()
        return self.par['DAPPAR']['pa'][0]
        



