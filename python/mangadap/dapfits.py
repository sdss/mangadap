# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Defines a class used to interface with the final maps files produced by
the MaNGA Data Analysis Pipeline (DAP).

*License*:
    Copyright (c) 2015, SDSS-IV/MaNGA Pipeline Group
        Licensed under BSD 3-clause license - see LICENSE.rst

*Source location*:
    $MANGADAP_DIR/python/mangadap/dapfits.py

*Imports and python version compliance*:
    ::

*Class usage examples*:
    Add some usage comments here!

*Revision history*:
    | **08 Jun 2016**: Original Implementation by K. Westfall (KBW)
    | **25 Aug 2016**: (KBW) Added :func:`DAPFits.unique_indices`,
        :func:`DAPFits.unique_mask`, and :func:`DAPFits.channel_map`

.. todo::
    - Merge this with dapmaps.py and dapcube.py?
    - Allow DAPFits to read/access both the MAPS and the LOGCUBE files,
      not just the former.

.. _astropy.io.fits.hdu.hdulist.HDUList: http://docs.astropy.org/en/v1.0.2/io/fits/api/hdulists.html
.. _astropy.io.fits.open: http://docs.astropy.org/en/stable/io/fits/api/files.html#astropy.io.fits.open
"""

from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from __future__ import unicode_literals

import sys
if sys.version > '3':
    long = int

import time
import os
import numpy
import warnings

from astropy.io import fits

from .util.fileio import channel_dictionary
from .util.exception_tools import print_frame
from .util.geometry import SemiMajorAxisCoo
from .config.defaults import default_drp_version, default_dap_version, default_analysis_path
from .config.defaults import default_dap_method_path, default_dap_file_name, default_dap_par_file
from .par.obsinput import ObsInputPar
from .dapmaps import DAPQualityBitMask, DAPMapsBitMask
from .mangafits import MaNGAFits

__author__ = 'Kyle Westfall'

class DAPFits:
    """
    Object used to hold properties of and read data from a DAP-produced file.

    Args:
        plate (int): Plate number
        ifudesign (int): IFU design
        method (str): Method type (e.g., ``'SPX-GAU-MILESTHIN'``)
        drpver (str): (**Optional**) DRP version, which is used to
            define the default DAP analysis path.  Default is defined by
            :func:`mangadap.config.defaults.default_drp_version`
        dapver (str): (**Optional**) DAP version, which is used to define
            the default DAP analysis path.  Default is defined by
            :func:`mangadap.config.defaults.default_dap_version`
        analysis_path (str): (**Optional**) The path to the top level
            directory containing the DAP output files for a given DRP
            and DAP version.  Default is defined by
            :func:`mangadap.config.defaults.default_analysis_path`.
        directory_path (str): (**Optional**) The exact path to the DAP
            file.  Default is defined by
            :func:`mangadap.config.defaults.default_dap_method_path`.
        reference_path (str): (**Optional**) The exact path of the DAP
            reference directory.  Default is defined by
            :func:`mangadap.config.defaults.default_dap_method_path`.
        par_file (str): (**Optional**) SDSS parameter file used to
            provide input parameters for the DAP.  Default is defined by
            :func:`mangadap.config.defaults.default_dap_par_file`.
        read (bool) : (**Optional**) Read the DAP file upon
            instantiation of the object.
        checksum (bool): (**Optional**) Flag for `astropy.io.fits.open`_
            checksum operation.

    Attributes:
        plate (int): Plate number
        ifudesign (int): IFU design
        method (str): Method type
        drpver (str): DRP version
        dapver (str): DAP version
        analysis_path (str): The path to the top level
            directory containing the DAP output files for a given DRP
            and DAP version.
        directory_path (str): The exact path to the DAP file.
        reference_path (str): The exact path of the DAP reference directory.
        par_file (str): SDSS parameter file used to provide input
            parameters for the DAP.
        dapqual (:class:`mangadap.dapmaps.DAPQualityBitMask`): Global
            quality bit
        bitmask (:class:`mangadap.dapmaps.DAPMapsBitMask`): Map mask
            bits
        hdu (`astropy.io.fits.hdu.hdulist.HDUList`_): HDUList read from
            the DRP file
        par (:class:`mangadap.par.ObsInputPar`): List of parameters used
            by the DAP.
        checksum (bool): Boolean to use for checksum in
            `astropy.io.fits.open`_.

    """
    def __init__(self, plate, ifudesign, method, drpver=None, dapver=None, analysis_path=None,
                 directory_path=None, reference_path=None, par_file=None, read=True,
                 checksum=False):

        # Set the attributes, forcing a known type
        self.plate = int(plate)
        self.ifudesign = int(ifudesign)
        self.method = method

        # TODO: Do we need drpver, dapver, analysis_path to be kept
        # by self?
        if directory_path is None:
            self.drpver = default_drp_version() if drpver is None else str(drpver)
            self.dapver = default_dap_version() if dapver is None else str(dapver)
            self.analysis_path = default_analysis_path(self.drpver, self.dapver) \
                                 if analysis_path is None else str(analysis_path)
            self.directory_path = default_dap_method_path(method, plate=self.plate,
                                                          ifudesign=self.ifudesign,
                                                          drpver=self.drpver, dapver=self.dapver,
                                                          analysis_path=self.analysis_path)
        else:
            self.drpver = None
            self.dapver = None
            self.analysis_path = None
            self.directory_path = str(directory_path)

        # Set the reference directory path
        self.reference_path = default_dap_method_path(method, plate=self.plate,
                                                      ifudesign=self.ifudesign, ref=True,
                                                      drpver=self.drpver, dapver=self.dapver,
                                                      analysis_path=self.analysis_path) \
                                    if reference_path is None else reference_path

        # drpver, dapver, and analysis_path can be None and par_file
        # will still only use the above defined directory_path
        self.par_file = default_dap_par_file(self.plate, self.ifudesign, 'CUBE',
                                             drpver=self.drpver, dapver=self.dapver,
                                             analysis_path=self.analysis_path,
                                             directory_path=self.reference_path) \
                                             if par_file is None else str(par_file)

        # Set the bitmasks
        self.dapqual = DAPQualityBitMask()
        self.bitmask = DAPMapsBitMask()

        # Read in the data
        self.hdu = None
        self.spatial_shape = None
        self.smap_ext = None      # Extensions with single maps
        self.mmap_ext = None      # Extensions with multiple maps
        self.par = None
        self.channel_dict = None
        self.checksum = checksum
        if read:
            self.open_hdu(checksum=self.checksum)
            self.read_par()


    def __del__(self):
        """
        Destroy the dapfile object, ensuring that the fits file is
        properly closed.
        """
        if self.hdu is None:
            return
        self.hdu.close()
        self.hdu = None


    def __getitem__(self, key):
        """Access elements of the hdu."""
        if self.hdu is None:
            return None
        return self.hdu[key]


    def open_hdu(self, permissions='readonly', checksum=False, quiet=True, restructure=True):
        """
        Open the fits file and save it to :attr:`hdu`; if :attr:`hdu` is
        not None, the function returns without re-reading the data.

        Args:
            permissions (string): (**Optional**) Open the fits file with
                these read/write permissions.
            checksum (bool): (**Optional**) Flag for
                `astropy.io.fits.open`_ checksum operation.  This
                overrides the internal :attr:`checksum` attribute **for
                the current operation only**.
            quiet (bool): (**Optional**) Suppress terminal output
            restructure (bool): (**Optional**) Restructure the data such
                that the ordering is (x,y,value); the fits data are
                stored as (value,y,x).  Set to False for quicker
                execution but beware of the index order!

        Raises:
            FileNotFoundError: Raised if the DAP file doesn't exist.

        """
        if self.hdu is not None:
            if not quiet:
                print('DAP file already open!')
            return

        inp = self.file_path()
        if not os.path.exists(inp):
            raise FileNotFoundError('Cannot open file: {0}'.format(inp))

        # Open the fits file with the requested read/write permission
        #check = self.checksum if checksum is None else checksum
        self.hdu = fits.open(inp, mode=permissions, checksum=checksum)

        self.spatial_shape = self.hdu['BINID'].data.shape
        print(self.spatial_shape)

        self.smap_ext = []      # Extensions with single maps
        self.mmap_ext = []      # Extensions with multiple maps
        for h in self.hdu:
            if h.data is not None:
                if len(h.shape) == 2:
                    self.smap_ext += [ h.name ]
                if len(h.shape) == 3:
                    self.mmap_ext += [ h.name ]

        if len(self.mmap_ext) > 0:
            self.channel_dict = { self.mmap_ext[0]: channel_dictionary(self.hdu, self.mmap_ext[0]) }
            for i in range(1,len(self.mmap_ext)):
                self.channel_dict[self.mmap_ext[i]] = channel_dictionary(self.hdu, self.mmap_ext[i])

        if not restructure:
            return

        if len(self.smap_ext) > 0:
            if not quiet:
                print('Restructuring single map extensions.')
            MaNGAFits.restructure_map(self.hdu, ext=self.smap_ext)
        if len(self.mmap_ext) > 0:
            if not quiet:
                print('Restructuring multi-map extensions.')
            MaNGAFits.restructure_cube(self.hdu, ext=self.mmap_ext)


    def read_par(self, quiet=True):
        """
        Open the parameter file and save it to :attr:`par`; if
        :attr:`par` is not None, the function returns without re-reading
        the data.

        Args:
            quiet (bool): Suppress terminal output

        """
        if self.par is not None:
            if not quiet:
                print('Parameter file already read!')
            return

        if not os.path.exists(self.par_file):
            warnings.warn('Input parameters unavailable; cannot open {0}'.format(self.par_file))
            return
            
        # Set the input parameters
        # TODO: Is all of this also in the file header?
        self.par = ObsInputPar.from_par_file(self.par_file)

    
    def file_name(self):
        """Return the name of the DAP file"""
        return default_dap_file_name(self.plate, self.ifudesign, self.method, mode='MAPS')


    def file_path(self):
        """Return the full path to the DAP file"""
        return os.path.join(self.directory_path, self.file_name())

    
    def list_channels(self, ext):
        """
        Provide a list of the channels in a given extension.
        """
        if self.hdu is None:
            self.open_hdu(checksum=self.checksum)
        if ext in self.smap_ext:
            print('{0} contains a single channel.'.format(ext))
            return
        print('Channels in extension {0}'.format(ext))
        print('#  {0:>3} {1}'.format('CH','KEY'))
        srt = numpy.argsort(list(self.channel_dict[ext].values()))
        for k, v in zip(numpy.array(list(self.channel_dict[ext].keys()))[srt],
                        numpy.array(list(self.channel_dict[ext].values()))[srt]):
            print('   {0:>3} {1}'.format(v,k))


    def thin_disk_polar_coo(self, xc=None, yc=None, rot=None, pa=None, inc=None, flip=False):
        r"""

        Calculate the disk-plane coordinates for all the binned spectra
        using :class:`mangadap.util.geometry.SemiMajorAxisCoo`.  See the
        documentation their for further information regarding the
        calculations.

        If not provided, the default position angle and inclination will
        be taken from the input parameter file, if available.  If the
        parameter file cannot be read, the default values for the
        position angle and ellipticity in
        :class:`mangadap.util.geometry.SemiMajorAxisCoo` are used.

        Args:
            xc,yc (float,float): (**Optional**) a reference on-sky
                position relative to the center of the projected plane
                (galaxy center).
            rot (float): (**Optional**) Cartesian rotation of the
                focal-plane relative to the sky-plane (+x toward East;
                +y toward North).
            pa (float): (**Optional**) On-sky position angle of the
                major axis of the ellipse traced on the sky by a circle
                in the projected plane, defined as the angle from North
                through East.
            inc (float): (**Optional**) Inclination of the plane, angle
                between the disk normal and the normal of the sky plane
                (inc=0 is a face-on plane).
            flip (bool): (**Optional**) Offset the provided position
                angle by 180 deg; e.g., to set the position angle will
                be defined along the receding major axis.

        Returns:
            numpy.ndarray: Two numpy arrays with the projected polar
                coordinates: :math:`R, \theta`.
        """

        # Flip if requested
        if pa is not None and flip:
            pa += 180
            if pa >= 360:
                pa -= 360

        # Position angle and inclination
        if pa is None or inc is None:
            try:
                self.read_par()                # Try to read the parameter file
            except:
                warnings.warn('Using default pa and/or inclination.')

            if self.par is not None:
                if pa is None:
                    pa = self.guess_position_angle(flip)
                if inc is None:
                    inc = self.guess_inclination()

        # Declare the galaxy coo object
        disk_plane = SemiMajorAxisCoo(xc=xc, yc=yc, rot=rot, pa=pa,
                                      ell=numpy.cos(numpy.radians(inc)))

        # Calculate the in-plane radius and azimuth for each bin
        if self.hdu is None:
            self.open_hdu(checksum=self.checksum)
        return disk_plane.polar(self.hdu['SPX_SKYCOO'].data[:,:,0],
                                self.hdu['SPX_SKYCOO'].data[:,:,1])

    
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
            q0 (float): (**Optional**) The intrinsic oblateness of the
                system.  If not provided, the system is assumed to be
                infinitely thin.

        Raises:
            Exception: Raised if the input intrinsic oblatenss is not
                less than one and greater than 0.
        """
        self.read_par()
        if q0 is None or q0 == 0.0:
            return numpy.degrees( numpy.arccos(1.0 - self.par['ell']) )

        if not (q0 < 1.0) or (q0 < 0.0):
            raise Exception('Intrinsic oblateness must be 0.0 <= q0 < 1.0!')

        q = 1.0 - self.par['ell']
        if q < q0:
            return 90.0

        q02 = q0*q0
        return numpy.degrees(numpy.arccos(numpy.sqrt( (q*q - q02) / (1.0 - q02) )))


    def guess_position_angle(self, flip=False):
        """
        Return the guess position angle from the parameter file.

        Args:
            flip (bool): (**Optional**) Offset the provided position
                angle by 180 deg.
        """
        self.read_par()
        pa = self.par['pa'] if not flip else self.par['pa'] + 180
        return pa if pa < 360 else pa-360
       

    def effective_radius(self):
        """
        Return the effective radius from the parameter file.
        """
        self.read_par()
        return self.par['reff']


#    def nsa_ellipticity(self):
#        """
#        Return the ellipticity from the parameter file.
#        """
#        self.read_par()
#        return self.par['ell']
#
#
    def guess_cz(self):
        """Return the guess redshift (cz) from the parameter file."""
        self.read_par()
        return self.par['vel']


    def unique_indices(self):
        """
        Use the BINID extension to select the spaxels with unique data.
        
        Returns:
            numpy.ndarray: Returns the list of indices **in the
            flattened array of a map** that are unique.

        """
        if self.hdu is None:
            self.open_hdu(checksum=self.checksum)

        unique_bins, indx = numpy.unique(self.hdu['BINID'].data.ravel(), return_index=True)
        return indx[unique_bins > -1]


    def unique_mask(self):
        """
        Construct a boolean mask for the unique bins.
        
        The returned map is True for spaxels that are part of a unique
        bin, Fals otherwise.
        """
        indx = self.unique_indices()
        msk = numpy.zeros(self.spatial_shape, dtype=numpy.bool)
        msk.ravel()[indx] = True
        return msk


    def channel_map(self, ext, channel=None, flag=None):
        """
        Return a 2D array with the map for the specified channel.
        """
        if channel is None:
            # Channel must be satisfied
            if len(self.hdu[ext].data.shape) > 2:
                raise ValueError('Must specify channel for extension {0}'.format(ext))
            # Return unmasked array
            if flag is None:
                return numpy.ma.MaskedArray(self.hdu[ext].data.copy())
            # Attempt to get the mask extension from the header
            try:
                mask_ext = self.hdu[ext].header['QUALDATA']
            except:
                raise ValueError('Mask extension not specified in header.  Unable to mask array.')
            # Return a masked array, masking the specified flags
            return numpy.ma.MaskedArray(self.hdu[ext].data.copy(),
                                        mask=self.bitmask.flagged(self.hdu[mask_ext].data,
                                                                  flag=flag))

        if channel not in list(self.channel_dict[ext].keys()):
            raise ValueError('No channel {0} in extension {1}'.format(ext, channel))

        # Return unmasked array
        if flag is None:
            return numpy.ma.MaskedArray(self.hdu[ext].data[:,:,self.channel_dict[ext][channel]])
        # Attempt to get the mask extension from the header
        try:
            mask_ext = self.hdu[ext].header['QUALDATA']
        except:
            raise ValueError('Mask extension not specified in header.  Unable to mask array.')
        
        return numpy.ma.MaskedArray(self.hdu[ext].data[:,:,self.channel_dict[ext][channel]],
                                    mask=self.bitmask.flagged(
                                      self.hdu[mask_ext].data[:,:,self.channel_dict[ext][channel]],
                                      flag=flag))



