"""

Defines a class used to interface with files produced by the MaNGA Data
Analysis Pipeline (DAP).

*Source location*:
    $MANGADAP_DIR/python/mangadap/dapfile.py

*Python2/3 compliance*::

    from __future__ import division
    from __future__ import print_function
    from __future__ import absolute_import
    
    import sys
    if sys.version > '3':
        long = int

*Imports*::

    import os.path
    from os import environ
    from astropy.io import fits
    import time
    import numpy
    from mangadap.util.exception_tools import print_frame
    from mangadap.util.yanny import yanny
    from mangadap.util.projected_disk_plane import projected_disk_plane
    from mangadap.util.defaults import default_drp_version, default_dap_version
    from mangadap.util.defaults import default_analysis_path, default_dap_directory_path
    from mangadap.util.defaults import default_dap_par_file, default_dap_plan_file
    from mangadap.util.defaults import default_dap_file_name

*Class usage examples*:

    .. todo::

        Add some usage comments here!

*Revision history*:
    | **16 Dec 2014**: Original Implementation by K. Westfall (KBW)
    | **19 Mar 2015**: (KBW) Added analysis_path to
        :class:`mangadap.dapfile.dapfile`
    | **23 Mar 2015**: (KBW) Fixed bug in parsing of niter in
        :func:`mangadap.dapfile.parse_dap_file_name`
    | **26 Apr 2015**: (KBW) Include DRP version in output path
    | **22 May 2015**: (KBW) Sphinx documentation tests.
    | **26 May 2015**: (KBW) Stripped out some of the basic slicing
        functions.  Added checksum=True when opening the DRP file.
    | **04 Jun 2015**: (KBW) Moved parse_dap_file_name to
        :func:`parser.parse_dap_file_name`
    | **04 Jul 2015**: (KBW) Allow the guess inclination to include an
        intrinsic oblateness.
    | **05 Aug 2015**: (KBW) Changed how directory path is set;
        previously required drpver, dapver, and analysis_path defaults,
        even if directory_path was provided directly.  **May need to add
        checks to other code to make sure drpver, dapver, and
        analysis_path are not None when directory_path has been directly
        defined.**

.. todo::
    - check that *bintype* is valid 
    - Allow default of niter=None, then search for the first file
          with the appropriate name and set niter

.. _astropy.io.fits.hdu.hdulist.HDUList: http://docs.astropy.org/en/v1.0.2/io/fits/api/hdulists.html

"""

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
from mangadap.util.defaults import default_drp_version, default_dap_version
from mangadap.util.defaults import default_analysis_path, default_dap_directory_path
from mangadap.util.defaults import default_dap_par_file, default_dap_plan_file
from mangadap.util.defaults import default_dap_file_name

__author__ = 'Kyle Westfall'

class dapfile:
    """
    Object used to hold properties of and read data from a DAP-produced file.

    Args:
        plate (int): Plate number
        ifudesign (int): IFU design
        mode (str): 3D mode of the DRP file; must be either 'RSS' or
            'CUBE'
        bintype (str): Binning type used by the DAP
        niter (int): Iteration number, typically the plan number in the
            plan parameter file.
        drpver (str): (Optional) DRP version, which is used to define
            the default DAP analysis path.  Default is defined by
            :func:`mangadap.util.defaults.default_drp_version`
        dapver (str): (Optional) DAP version, which is used to define
            the default DAP analysis path.  Default is defined by
            :func:`mangadap.util.defaults.default_dap_version`
        analysis_path (str): (Optional) The path to the top level
            directory containing the DAP output files for a given DRP
            and DAP version.  Default is defined by
            :func:`mangadap.util.defaults.default_analysis_path`.
        directory_path (str): (Optional) The exact path to the DAP file.
            Default is defined by
            :func:`mangadap.util.defaults.default_dap_directory_path`.
        par_file (str): SDSS parameter file used to provide input
            parameters for the DAP.  Default is defined by
            :func:`mangadap.util.defaults.default_dap_par_file`.
        read (bool) : (Optional) Read the DAP file upon instantiation of
            the object.

    Attributes:

        plate (int): Plate number
        ifudesign (int): IFU design
        mode (str): 3D mode of the DRP file; must be either 'RSS' or
            'CUBE'
        bintype (str): Binning type used by the DAP
        niter (int): Iteration number, typically the plan number in the
            plan parameter file.
        drpver (str): DRP version
        dapver (str): DAP version
        analysis_path (str): The path to the top level
            directory containing the DAP output files for a given DRP
            and DAP version.
        directory_path (str): The exact path to the DAP file.
        par_file (str): SDSS parameter file used to provide input
            parameters for the DAP.
        hdu (`astropy.io.fits.hdu.hdulist.HDUList`_): HDUList read from
            the DRP file
        par (:class:`mangadap.util.yanny.yanny`): List of parameters
            used by the DAP.

    """
    def __init__(self, plate, ifudesign, mode, bintype, niter, drpver=None, dapver=None,
                 analysis_path=None, directory_path=None, par_file=None, read=True):
        # Set the attributes, forcing a known type
        self.plate = int(plate)
        self.ifudesign = int(ifudesign)
        self.mode = str(mode)
        self.bintype = str(bintype)
        self.niter = int(niter)

        if directory_path is None:
            self.drpver = default_drp_version() if drpver is None else str(drpver)
            self.dapver = default_dap_version() if dapver is None else str(dapver)
            self.analysis_path = default_analysis_path(self.drpver, self.dapver) \
                                 if analysis_path is None else str(analysis_path)
            self.directory_path = default_dap_directory_path(self.drpver, self.dapver,
                                                             self.analysis_path, self.plate,
                                                             self.ifudesign)
        else:
            self.drpver = None
            self.dapver = None
            self.analysis_path = None
            self.directory_path = str(directory_path)

        self.par_file = default_dap_par_file(self.drpver, self.dapver, self.analysis_path, \
                                             self.directory_path, self.plate, self.ifudesign, \
                                             self.mode) \
                                             if par_file is None else str(par_file)

        self.hdu = None
        self.par = None

        if read:
            self.open_hdu()


    def __del__(self):
        """
        Destroy the dapfile object, ensuring that the fits file is
        properly closed.
        """
        if self.hdu is None:
            return
        self.hdu.close()
        self.hdu = None


#   def _twod_image_data(self, exten, indx):
#       self.open_hdu()
#       try:
#           data = self.hdu[exten].data[:,indx]
#       except IndexError as e:
#           print_frame('IndexError')
#           print('{0}'.format(e))
#           return None
#       else:
#           return data
    
    
    def open_hdu(self, permissions='readonly', quiet=True):
        """
        Open the fits file and save it to :attr:`hdu`; if :attr:`hdu` is
        not None, the function returns without re-reading the data.

        Args:
            permissions (string): (Optional) Open the fits file with
                these read/write permissions.
            quiet (bool): Suppress terminal output

        Raises:
            Exception: Raised if the DAP file doesn't exist.

        """
        if self.hdu is not None:
            if not quiet:
                print('DAP file already open!')
            return

        inp = self.file_path()
        if not os.path.exists(inp):
            raise Exception('Cannot open file: {0}'.format(inp))

        # Open the fits file with the requested read/write permission
        self.hdu = fits.open(inp, mode=permissions, checksum=True)


    def read_par(self, quiet=True):
        """
        Open the parameter file and save it to :attr:`par`; if
        :attr:`par` is not None, the function returns without re-reading
        the data.

        Args:
            quiet (bool): Suppress terminal output

        Raises:
            Exception: Raised if the parameter file doesn't exist.

        """
        if self.par is not None:
            if not quiet:
                print('Parameter file already read!')
            return

        if not os.path.exists(self.par_file):
            raise Exception('Cannot open file: {0}'.format(inp))
            
        self.par = yanny(filename=self.par_file, np=True)

    
    def file_name(self):
        """Return the name of the DAP file"""
        return default_dap_file_name(self.plate, self.ifudesign, self.mode,
                                     self.bintype, self.niter)


    def file_path(self):
        """Return the full path to the DAP file"""
        return os.path.join(self.directory_path, self.file_name())


    def header(self):
        """Return the primary header"""
        self.open_hdu()
        return self.hdu[0].header


#    def sn_stats(self):
#        self.open_hdu()
#        return self.hdu[0].header['SNWAVE1'], self.hdu[0].header['SNWAVE2'], \
#               self.hdu[0].header['MINSNBIN']
#
#
#    def spaxel_size(self):
#        self.open_hdu()
#        return self.hdu[0].header['SPAXDX'], self.hdu[0].header['SPAXDY']
#
#
#    def bin_ston(self):
#        if self.bintype is not 'STON':
#            raise Exception('Binning type was not STON!')
#        self.open_hdu()
#        return self.hdu[0].header['BINSN']
#
#
#    def analysis_stats(self):
#        self.open_hdu()
#        return self.hdu[0].header['FITWAVE1'], self.hdu[0].header['FITWAVE2'], \
#               self.hdu[0].header['MINSNFIT']
#
#
#    def table_column(self, exten, col):
#        self.open_hdu()
#        try:
#            data = self.hdu[exten].data[col]
#        except KeyError as e:
#            print_frame('KeyError')
#            print('Could not find extension ({0}) and/or column ({1})!'.format(exten, col))
#            return None
#        else:
#            return data
#
#    
#    def bin_wavelength(self):
#        self.open_hdu()
#        return self.hdu['WAVE'].data
#
#
#    def bin_spectrum(self, indx):
#        return self._twod_image_data('FLUX', indx)
#
#
#    def bin_inverse_variance(self, indx):
#        return self._twod_image_data('IVAR', indx)
#            
#
#    def bin_mask(self, indx):
#        return self._twod_image_data('MASK', indx)
#
#
#    def fit_stars_minus_gas(self, indx):
#        bestfit = self._twod_image_data('SGMOD', indx)
#        emlmodel = self._twod_image_data('ELMOD', indx)
#        if bestfit is not None and emlmodel is not None:
#            return (bestfit - emlmodel)
#        else:
#            return None
#
#            
#    def fit_stars_and_gas(self, indx):
#        return self._twod_image_data('SGMOD', indx)
#
#    
#    def fit_eml_only_ew(self, indx):
#        return self._twod_image_data('ELOMEW', indx)
#            
#
#    def fit_eml_only_fb(self, indx):
#        return self._twod_image_data('ELOMFB', indx)


    def bin_disk_polar_coo(self, xc=None, yc=None, rot=None, pa=None, inc=None, flip=False):
        r"""

        Calculate the disk-plane coordinates for all the binned spectra
        using
        :class:`mangadap.util.projected_disk_plane.projected_disk_plane`.
        See the documentation their for further information regarding
        the calculations.

        If not provided, the default position angle and inclination will
        be taken from the input parameter file, if available.  If the
        parameter file cannot be read, the default values for the
        position angle and inclination in
        :class:`mangadap.util.projected_disk_plane.projected_disk_plane`
        are used.

        Args:
            xc,yc (float,float): (Optional) a reference on-sky position
                relative to the center of the projected plane (galaxy
                center).
            rot (float): (Optional) Cartesian rotation of the
                focal-plane relative to the sky-plane (+x toward East;
                +y toward North).
            pa (float): (Optional) On-sky position angle of the major
                axis of the ellipse traced on the sky by a circle in the
                projected plane, defined as the angle from North through
                East.
            inc (float): (Optional) Inclination of the plane, angle
                between the disk normal and the normal of the sky plane
                (inc=0 is a face-on plane).
            flip (bool): (Optional) Offset the provided position angle
                by 180 deg; e.g., to set the position angle will be
                defined along the receding major axis.

        Returns:
            numpy.ndarray: Two numpy arrays with the projected polar
                coordinates: :math:`R, \theta`.
        """

        # Position angle and inclination
        if pa is None or inc is None:
            try:
                self.read_par()                # Try to read the parameter file
            except:
                print('WARNING: Using default pa and/or inclination.')

            if self.par is not None:
                if pa is None:
                    pa = self.guess_position_angle(flip)
                if inc is None:
                    inc = self.guess_inclination()

        # Declare the galaxy coo object
        disk_plane = projected_disk_plane(xc, yc, rot, pa, inc)

        # Calculate the in-plane radius and azimuth for each bin
        self.open_hdu()
        nbin = self.hdu['BINS'].data['BINXRL'].size
        radius = numpy.zeros(nbin, dtype=numpy.float64)
        azimuth = numpy.zeros(nbin, dtype=numpy.float64)
        for i in range(0,nbin):
            radius[i], azimuth[i] = disk_plane.polar(self.hdu['BINS'].data['BINXRL'][i], \
                                                     self.hdu['BINS'].data['BINYRU'][i])
            
        return radius, azimuth

    
    def guess_cz(self):
        """Return the guess redshift (cz) from the parameter file."""
        self.read_par()
        return self.par['DAPPAR']['vel'][0]


    def guess_inclination(self, q0=None):
        r"""

        Return a guess inclination based on the isophotal ellipticity
        from the input parameter file using the equation:

        .. math::
        
            \cos^2 i = \frac{ q^2 - q_0^2 }{ 1 - q_0^2 },

        where :math:`q` is the observed oblateness (the semi-minor to
        semi-major axis ratio, :math:`q = 1-\epsilon = b/a`) and
        :math:`q_0` is the intrinsic oblateness.
        
        If :math:`q < q_0`, the function returns a 90 degree
        inclination.
        
        Args:
            q0 (float): (Optional) The intrinsic oblateness of the
                system.  If not provided, the system is assumed to be
                infinitely thin.

        Raises:
            Exception: Raised if the input intrinsic oblatenss is not
                less than one and greater than 0.
        """
        self.read_par()
        if q0 is None or q0 == 0.0:
            return numpy.degrees( numpy.arccos(1.0 - self.par['DAPPAR']['ell'][0]) )

        if not (q0 < 1.0) or (q0 < 0.0):
            raise Exception('Intrinsic oblateness must be 0.0 <= q0 < 1.0!')

        q = 1.0 - self.par['DAPPAR']['ell'][0]
        if q < q0:
            return 90.0

        q02 = q0*q0
        return numpy.degrees(numpy.arccos(numpy.sqrt( (q*q - q02) / (1.0 - q02) )))


    def guess_position_angle(self, flip=False):
        """
        Return the guess position angle from the parameter file.

        Args:
            flip (bool): (Optional) Offset the provided position angle
                by 180 deg.
        """
        self.read_par()
        pa = self.par['DAPPAR']['pa'][0] if not flip else self.par['DAPPAR']['pa'][0] + 180
        return pa if pa < 360 else pa-360
       

    def effective_radius(self):
        """
        Return the effective radius from the parameter file.
        """
        self.read_par()
        return self.par['DAPPAR']['reff'][0]


    def nsa_ellipticity(self):
        """
        Return the ellipticity from the parameter file.
        """
        self.read_par()
        return self.par['DAPPAR']['ell'][0]



