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

        from __future__ import print_function
        from __future__ import division
        from __future__ import absolute_import
        from __future__ import unicode_literals

        import sys
        if sys.version > '3':
            long = int

        import logging
        import time
        import os
        import numpy
        import warnings

        from astropy.wcs import WCS
        from astropy.io import fits
        import astropy.constants

        from .drpfits import DRPFits, DRPQuality3DBitMask
        from .util.fitsutil import DAPFitsUtil
        from .util.bitmask import BitMask
        from .util.log import log_output
        from .util.fileio import channel_dictionary
        from .util.exception_tools import print_frame
        from .par.obsinput import ObsInputPar
        from .config.defaults import default_drp_version, dap_source_dir, default_dap_version
        from .config.defaults import default_dap_par_file, default_analysis_path
        from .config.defaults import default_dap_method, default_dap_method_path
        from .config.defaults import default_dap_file_name
        from .proc.reductionassessments import ReductionAssessment
        from .proc.spatiallybinnedspectra import SpatiallyBinnedSpectra
        from .proc.stellarcontinuummodel import StellarContinuumModel
        from .proc.emissionlinemoments import EmissionLineMoments
        from .proc.emissionlinemodel import EmissionLineModel
        from .proc.spectralindices import SpectralIndices

*Class usage examples*:
    Add some usage comments here!

*Revision history*:
    | **08 Jun 2016**: Original Implementation by K. Westfall (KBW)
    | **25 Aug 2016**: (KBW) Added :func:`DAPFits.unique_indices`,
        :func:`DAPFits.unique_mask`, and :func:`DAPFits.channel_map`
    | **Feb 2017**: (KBW) Merged construct maps and cube functions here
        and removed dapmaps.py and dapcube.py.
    | ** 3 May 2017**: (KBW) Allow output maps and cube types to be
        single-precision (float32) floats, and changed BINID from int
        (defaults to 64-bit) to int32.
    | **24 Aug 2017**: (KBW) Added BADGEOM DAPQUAL bit; signifies 
        invalid photometric geometry parameters used when instantiating
        the :class:`mangadap.par.obsinput.ObsInputPar` object.
    | **29 Aug 2017**: (KBW) Add S/N metrics to the MAPS primary header
    | **30 Aug 2017**: (KBW) Convert some changes from Xihan Ji.
    | **19 Mar 2018**: (KBW) Mask spectral indices with incomplete band
        coverage as unreliable.

.. todo::
    - Allow DAPFits to read/access both the MAPS and the LOGCUBE files,
      not just the former.

.. _astropy.io.fits.hdu.hdulist.HDUList: http://docs.astropy.org/en/v1.0.2/io/fits/api/hdulists.html
.. _astropy.io.fits.open: http://docs.astropy.org/en/stable/io/fits/api/files.html#astropy.io.fits.open
.. _astropy.io.fits.Header: http://docs.astropy.org/en/stable/io/fits/api/headers.html

"""

from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from __future__ import unicode_literals

import sys
if sys.version > '3':
    long = int

import logging
import time
import os
import numpy
import warnings

from astropy.wcs import WCS
from astropy.io import fits
import astropy.constants

from .drpfits import DRPFits, DRPQuality3DBitMask
from .util.fitsutil import DAPFitsUtil
from .util.bitmask import BitMask
from .util.log import log_output
from .util.fileio import channel_dictionary
from .util.exception_tools import print_frame
from .util.geometry import SemiMajorAxisCoo
from .util.covariance import Covariance
from .par.obsinput import ObsInputPar
from .config.defaults import default_drp_version, dap_source_dir, default_dap_version
from .config.defaults import default_dap_par_file, default_analysis_path
from .config.defaults import default_dap_method, default_dap_method_path
from .config.defaults import default_dap_file_name
from .proc.reductionassessments import ReductionAssessment
from .proc.spatiallybinnedspectra import SpatiallyBinnedSpectra
from .proc.stellarcontinuummodel import StellarContinuumModel
from .proc.emissionlinemoments import EmissionLineMoments
from .proc.emissionlinemodel import EmissionLineModel
from .proc.spectralindices import SpectralIndices

from matplotlib import pyplot

#-----------------------------------------------------------------------
class DAPQualityBitMask(BitMask):
    """
    .. todo::
        - Force read IDLUTILS version as opposed to internal one?
    """
    def __init__(self, dapsrc=None):
        _dapsrc = dap_source_dir() if dapsrc is None else str(dapsrc)
        BitMask.__init__(self, ini_file=os.path.join(_dapsrc, 'python', 'mangadap', 'config',
                                                     'bitmasks', 'dap_quality_bits.ini'))


class DAPMapsBitMask(BitMask):
    """
    .. todo::
        - Force read IDLUTILS version as opposed to internal one?
    """
    def __init__(self, dapsrc=None):
        _dapsrc = dap_source_dir() if dapsrc is None else str(dapsrc)
        BitMask.__init__(self, ini_file=os.path.join(_dapsrc, 'python', 'mangadap', 'config',
                                                     'bitmasks', 'dap_maps_bits.ini'))


class DAPCubeBitMask(BitMask):
    """
    .. todo::
        - Force read IDLUTILS version as opposed to internal one?
    """
    def __init__(self, dapsrc=None):
        _dapsrc = dap_source_dir() if dapsrc is None else str(dapsrc)
        BitMask.__init__(self, ini_file=os.path.join(_dapsrc, 'python', 'mangadap', 'config',
                                                     'bitmasks', 'dap_cube_bits.ini'))


#-----------------------------------------------------------------------
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


#    def __del__(self):
#        """
#        Destroy the dapfile object, ensuring that the fits file is
#        properly closed.
#        """
#        if self.hdu is None:
#            return
#        self.hdu.close()
#        self.hdu = None


    def __getitem__(self, key):
        """Access elements of the hdu."""
        if self.hdu is None:
            return None
        return self.hdu[key]


    def open_hdu(self, permissions='readonly', checksum=False, quiet=True):
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
#        self.hdu = fits.open(inp, mode=permissions, checksum=checksum)
        self.hdu = DAPFitsUtil.read(inp, permissions=permissions, checksum=checksum)

        self.spatial_shape = self.hdu['BINID'].data.shape
#        print(self.spatial_shape)

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

#        if not restructure:
#            return
#
#        if len(self.smap_ext) > 0:
#            if not quiet:
#                print('Restructuring single map extensions.')
#            DAPFitsUtil.restructure_map(self.hdu, ext=self.smap_ext)
#        if len(self.mmap_ext) > 0:
#            if not quiet:
#                print('Restructuring multi-map extensions.')
#            DAPFitsUtil.restructure_cube(self.hdu, ext=self.mmap_ext)


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


#-----------------------------------------------------------------------
class construct_maps_file:
    """
    Construct a DAP MAPS file.

    Set as a class for coherency reasons, but should not be used as an
    object!

    Should force all intermediate objects to be provided.

    """
    def __init__(self, drpf, obs=None, rdxqa=None, binned_spectra=None, stellar_continuum=None,
                 emission_line_moments=None, emission_line_model=None, spectral_indices=None,
                 nsa_redshift=None, dapsrc=None, dapver=None, analysis_path=None,
                 directory_path=None, output_file=None, clobber=True, loggers=None, quiet=False,
                 single_precision=False):

        #---------------------------------------------------------------
        # Initialize the reporting
        self.loggers = None if loggers is None else loggers
        self.quiet = quiet

        self.float_dtype = 'float32' if single_precision else 'float'

        #---------------------------------------------------------------
        # Check input types
        confirm_dap_types(drpf, obs, rdxqa, binned_spectra, stellar_continuum,
                          emission_line_moments, emission_line_model, spectral_indices)

#        print('inside construct')
#        bin_indx = stellar_continuum['BINID'].data.copy().ravel()
#        pyplot.imshow(DAPFitsUtil.reconstruct_map(drpf.spatial_shape, bin_indx, #vel),
#                                                  stellar_continuum['PAR'].data['KIN'][:,0]),
#                      origin='lower', interpolation='nearest')
#        pyplot.show()
#
#        pyplot.imshow(DAPFitsUtil.reconstruct_map(drpf.spatial_shape, bin_indx, #vel),
#                                                  stellar_continuum['PAR'].data['KIN'][:,1]),
#                      origin='lower', interpolation='nearest')
#        pyplot.show()

        #---------------------------------------------------------------
        # Set the output paths
        self.drpf = drpf
        self.method = None
        self.directory_path = None
        self.output_file = None
        self._set_paths(directory_path, dapver, analysis_path, output_file, binned_spectra,
                        stellar_continuum)

        # Save input for reference
        self.spatial_shape = self.drpf.spatial_shape
        self.nsa_redshift = nsa_redshift
#        self.multichannel_arrays = None
#        self._set_multichannel_arrays(emission_line_moments, emission_line_model, spectral_indices)
#        self.singlechannel_arrays = None
#        self._set_singlechannel_arrays(emission_line_moments, emission_line_model, spectral_indices)

        # Report
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, '-'*50)
            log_output(loggers, 1, logging.INFO, '{0:^50}'.format('CONSTRUCTING OUTPUT MAPS'))
            log_output(self.loggers, 1, logging.INFO, '-'*50)
            log_output(self.loggers, 1, logging.INFO, 'Output path: {0}'.format(
                                                                            self.directory_path))
            log_output(self.loggers, 1, logging.INFO, 'Output file: {0}'.format(
                                                                            self.output_file))
            log_output(self.loggers, 1, logging.INFO, 'Output maps have shape {0}'.format(
                                                                            self.spatial_shape))
            log_output(self.loggers, 1, logging.INFO, 'NSA redshift: {0}'.format(self.nsa_redshift))

        #---------------------------------------------------------------
        # Check if the file already exists
        ofile = os.path.join(self.directory_path, self.output_file)
        if os.path.isfile(ofile) and not clobber:
            # TODO: Perform some checks to make sure the existing file
            # has the correct content?
            warnings.warn('Output file exists!  Set clobber=True to overwrite.')
            return

        #---------------------------------------------------------------
        # Initialize the primary header
        prihdr = DAPFitsUtil.initialize_dap_primary_header(self.drpf, maskname='MANGA_DAPPIXMASK')
        # Add the DAP method
        prihdr['DAPTYPE'] = (default_dap_method(binned_spectra=binned_spectra,
                                                stellar_continuum=stellar_continuum),
                             'DAP analysis method')
        # Add the format of this file
        prihdr['DAPFRMT'] = ('MAPS', 'DAP data file format')

        # Get the base map headers
        self.multichannel_maphdr \
                = DAPFitsUtil.build_map_header(self.drpf,
                                'K Westfall & B Andrews <westfall@ucolick.org, andrewsb@pitt.edu>',
                                               multichannel=True, maskname='MANGA_DAPPIXMASK')
        self.singlechannel_maphdr \
                = DAPFitsUtil.build_map_header(self.drpf,
                                'K Westfall & B Andrews <westfall@ucolick.org, andrewsb@pitt.edu>',
                                               maskname='MANGA_DAPPIXMASK')

        #---------------------------------------------------------------
        # Initialize the pixel mask
        self.bitmask = DAPMapsBitMask(dapsrc=dapsrc)

        # Get the mask for the binned spectra
        self.bin_mask = self._build_binning_mask(binned_spectra)

        #---------------------------------------------------------------
        # Construct the hdu list for each input object.
        # Reduction assessments:
        rdxqalist = self.reduction_assessment_maps(prihdr, obs, rdxqa)
        # Construct the BINID extension
        binidlist = combine_binid_extensions(self.drpf, binned_spectra, stellar_continuum,
                                             emission_line_moments, emission_line_model,
                                             spectral_indices, dtype='int32')
        # Binned spectra:
        bspeclist = self.binned_spectra_maps(prihdr, obs, binned_spectra)
        # Stellar-continuum fits:
        contlist = self.stellar_continuum_maps(prihdr, stellar_continuum)
        # Emission-line moments:
        emlmomlist = self.emission_line_moment_maps(prihdr, emission_line_moments)
        # Emission-line models:
        elmodlist = self.emission_line_model_maps(prihdr, emission_line_model)

#        print(elmodlist[10].name)
#        print(elmodlist[10].data.shape)
#        print(numpy.sum(numpy.invert(numpy.isfinite(elmodlist[10].data))))

        # Spectral indices:
        sindxlist = self.spectral_index_maps(prihdr, spectral_indices)

        # Save the data to the hdu attribute
        prihdr = add_snr_metrics_to_header(prihdr, self.drpf, rdxqalist[1].data[:,:,1].ravel(),
                                           dapsrc=dapsrc)
        
        prihdr = finalize_dap_primary_header(prihdr, self.drpf, obs, binned_spectra,
                                             stellar_continuum, dapsrc=dapsrc,
                                             loggers=self.loggers, quiet=self.quiet)
        self.hdu = fits.HDUList([ fits.PrimaryHDU(header=prihdr),
                                  *rdxqalist,
                                  *binidlist,
                                  *bspeclist,
                                  *contlist,
                                  *emlmomlist,
                                  *elmodlist,
                                  *sindxlist
                                ])

#        print('maps')
#        pyplot.imshow(self.hdu['STELLAR_VEL'].data, origin='lower', interpolation='nearest')
#        pyplot.show()
#
#        pyplot.imshow(self.hdu['STELLAR_SIGMA'].data, origin='lower', interpolation='nearest')
#        pyplot.show()

        #---------------------------------------------------------------
        # TEMPORARY FLAGS:
        # Flag the Gaussian-fitted flux as unreliable if the summed flux
        # is not within the factor below
#        if emission_line_moments is not None and emission_line_model is not None \
#                    and self.hdu['EMLINE_GFLUX'].data.shape != self.hdu['EMLINE_SFLUX'].data.shape:
#            warnings.warn('Cannot compare emission-line moment and model fluxes!')
#        elif emission_line_moments is not None and emission_line_model is not None \
#                    and self.hdu['EMLINE_GFLUX'].data.shape == self.hdu['EMLINE_SFLUX'].data.shape:
#            factor = 5.0
#            indx = (self.hdu['EMLINE_GFLUX_MASK'].data == 0) \
#                    & (self.hdu['EMLINE_SFLUX_MASK'].data == 0) \
#                    & ( (self.hdu['EMLINE_GFLUX'].data < self.hdu['EMLINE_SFLUX'].data/factor)
#                        | (self.hdu['EMLINE_GFLUX'].data > self.hdu['EMLINE_SFLUX'].data*factor) ) \
#                    & ( (emission_line_moments['BINID'].data > -1)
#                        & (emission_line_model['BINID'].data > -1) )[:,:,None]
#            print('unreliable Gaussian flux compared to summed flux: ', numpy.sum(indx))
#            if numpy.sum(indx) > 0:
#                self.hdu['EMLINE_GFLUX_MASK'].data[indx] \
#                    = self.bitmask.turn_on(self.hdu['EMLINE_GFLUX_MASK'].data[indx], 'UNRELIABLE')
#                self.hdu['EMLINE_GVEL_MASK'].data[indx] \
#                    = self.bitmask.turn_on(self.hdu['EMLINE_GVEL_MASK'].data[indx], 'UNRELIABLE')
#                self.hdu['EMLINE_GSIGMA_MASK'].data[indx] \
#                    = self.bitmask.turn_on(self.hdu['EMLINE_GSIGMA_MASK'].data[indx], 'UNRELIABLE')

            # TODO: Add if EMLINE_GSIGMA < EMLINE_INSTSIGMA !!

        # Flag the stellar velocity dispersions measured in spectra with
        # S/N<10 as unreliable
#        indx = (binned_spectra['BINID'].data > -1) & (self.hdu['BIN_SNR'].data < 10) \
#                    & (stellar_continuum['BINID'].data > -1) \
#                    & (self.hdu['STELLAR_SIGMA_MASK'].data == 0)
#        print('unreliable sigma because of low S/N: ', numpy.sum(indx))
#        if numpy.sum(indx):
#            self.hdu['STELLAR_SIGMA_MASK'].data[indx] \
#                    = self.bitmask.turn_on(self.hdu['STELLAR_SIGMA_MASK'].data[indx], 'UNRELIABLE')

        # Flag any inverse variances that are not positive as DONOTUSE
        # and MATHERROR
        ext = [ 'BIN_MFLUX', 'STELLAR_VEL', 'STELLAR_SIGMA', 'EMLINE_SFLUX', 'EMLINE_SEW',
                'EMLINE_GFLUX', 'EMLINE_GEW', 'EMLINE_GVEL', 'EMLINE_GSIGMA', 'SPECINDEX' ]
        for e in ext:
            indx = numpy.invert(self.bitmask.flagged(self.hdu['{0}_MASK'.format(e)].data)) \
                            & numpy.invert(self.hdu['{0}_IVAR'.format(e)].data > 0)
            if numpy.sum(indx) > 0:
                self.hdu['{0}_MASK'.format(e)].data[indx] \
                    = self.bitmask.turn_on(self.hdu['{0}_MASK'.format(e)].data[indx], 'MATHERROR')
                self.hdu['{0}_MASK'.format(e)].data[indx] \
                    = self.bitmask.turn_on(self.hdu['{0}_MASK'.format(e)].data[indx], 'DONOTUSE')

        #---------------------------------------------------------------

        # Check that the path exists
        if not os.path.isdir(self.directory_path):
            os.makedirs(self.directory_path)
        # Write the maps file
        DAPFitsUtil.write(self.hdu, ofile, clobber=clobber, checksum=True, loggers=self.loggers,
                          quiet=self.quiet)
        # End
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, '-'*50)


    def _set_paths(self, directory_path, dapver, analysis_path, output_file, binned_spectra,
                   stellar_continuum):
        """
        Set the paths relevant to the map file.
        """
        # The output method directory is, for now, the combination of
        # the binned_spectrum and stellar_continuum method keys
        if directory_path is None and (binned_spectra is None or stellar_continuum is None):
            raise ValueError('Could not define output directory path.')

        # Set the output directory path
        # TODO: Get DAP version from __version__ string
        self.method = default_dap_method(binned_spectra=binned_spectra,
                                         stellar_continuum=stellar_continuum) \
                                if directory_path is None else None
        self.directory_path = default_dap_method_path(self.method, plate=self.drpf.plate,
                                                      ifudesign=self.drpf.ifudesign,
                                                      drpver=self.drpf.drpver, dapver=dapver,
                                                      analysis_path=analysis_path) \
                                                if directory_path is None else str(directory_path)

        # Set the output file
        self.output_file = default_dap_file_name(self.drpf.plate, self.drpf.ifudesign,
                                                 self.method, mode='MAPS') \
                                    if output_file is None else str(output_file)


    def _consolidate_donotuse(self, mask):
        return self.bitmask.consolidate(mask, [ 'NOCOV', 'LOWCOV', 'DEADFIBER', 'FORESTAR',
                                                'NOVALUE', 'MATHERROR', 'FITFAILED', 'NEARBOUND' ],
                                        'DONOTUSE')


    def _build_binning_mask(self, binned_spectra):

        # Marginalize the DRP mask across wavelengths
        mask = DAPFitsUtil.marginalize_mask(self.drpf['MASK'].data,
                                            [ 'NOCOV', 'LOWCOV', 'DEADFIBER', 'FORESTAR',
                                              'DONOTUSE' ], self.drpf.bitmask, self.bitmask)

        # Add any bits not included in the binning algorithm
        indx = binned_spectra.bitmask.flagged(binned_spectra['MAPMASK'].data,
                                              flag=['LOW_SPECCOV', 'LOW_SNR'])
        mask[indx] = self.bitmask.turn_on(mask[indx], 'NOVALUE')

        # Marginalize the binned spectra mask just for NONE_IN_STACK and
        # convert it to NOVALUE
        bin_mask = DAPFitsUtil.marginalize_mask(binned_spectra['MASK'].data, 'NONE_IN_STACK',
                                                binned_spectra.bitmask, self.bitmask,
                                                out_flag='NOVALUE', dispaxis=1)
        # Reconstruct the full mask with just this flag
        bin_mask = DAPFitsUtil.reconstruct_map(self.spatial_shape,
                                               binned_spectra['BINID'].data.ravel(), bin_mask)
        # Add these bits to the full mask
        indx = self.bitmask.flagged(bin_mask, 'NOVALUE')
        mask[indx] = self.bitmask.turn_on(mask[indx], 'NOVALUE')
        
        # Return the mask after consolidating flags into DONOTUSE
        return self._consolidate_donotuse(mask)


    def _stellar_continuum_mask_to_map_mask(self, stellar_continuum, sc_mask, for_dispersion=False):
        """

        sc_mask constains the reconstructed map of the masks in each
        stellar continuum fit (using the stellar continuum bitmask).

        Propagate and consolidate the stellar-continuum masks to the map
        pixel mask.

        DIDNOTUSE, FORESTAR propagated from already existing mask
            (self.bin_mask)

        LOW_SNR, NO_FIT, INSUFFICIENT_DATA propagated to NOVALUE

        FIT_FAILED, NEAR_BOUND propagated to FITFAILED and NEARBOUND

        NEGATIVE_WEIGHTS consolidated into UNRELIABLE; if
            dispersion=True, include BAD_SIGMA in this
        """

        # Copy the binning mask
        mask = self.bin_mask.copy()

        # Use the stellar continuum map mask to flag low S/N
        # bins/spaxels as NOVALUE
        indx = stellar_continuum.bitmask.flagged(stellar_continuum['MAPMASK'].data, flag='LOW_SNR')
        mask[indx] = self.bitmask.turn_on(mask[indx], 'NOVALUE')

        # Consolidate to NOVALUE
        flgd = stellar_continuum.bitmask.flagged(sc_mask, flag=['NO_FIT', 'INSUFFICIENT_DATA' ])
        mask[flgd] = self.bitmask.turn_on(mask[flgd], 'NOVALUE')

        # Copy FIT_FAILED
        flgd = stellar_continuum.bitmask.flagged(sc_mask, flag='FIT_FAILED')
        mask[flgd] = self.bitmask.turn_on(mask[flgd], 'FITFAILED')

        # Copy NEAR_BOUND
        flgd = stellar_continuum.bitmask.flagged(sc_mask, flag='NEAR_BOUND')
        mask[flgd] = self.bitmask.turn_on(mask[flgd], 'NEARBOUND')

        # Consolidate to UNRELIABLE
        flgd = stellar_continuum.bitmask.flagged(sc_mask, flag='NEGATIVE_WEIGHTS')
        mask[flgd] = self.bitmask.turn_on(mask[flgd], 'UNRELIABLE')

        # Convert MIN_SIGMA into a NEARBOUND for the dispersion
        if for_dispersion:
            flgd = stellar_continuum.bitmask.flagged(sc_mask, flag='BAD_SIGMA')
            mask[flgd] = self.bitmask.turn_on(mask[flgd], 'UNRELIABLE')

            flgd = stellar_continuum.bitmask.flagged(sc_mask, flag='MIN_SIGMA')
            mask[flgd] = self.bitmask.turn_on(mask[flgd], 'NEARBOUND')

        return self._consolidate_donotuse(mask)


    def _emission_line_moment_mask_to_map_mask(self, emission_line_moments, elm_mask):
        """
        Propagate and consolidate the emission-line moment masks to the map
        pixel mask.

        DIDNOTUSE, FORESTAR propagated from already existing mask (self.bin_mask)

        MAIN_EMPTY, BLUE_EMPTY, RED_EMPTY, UNDEFINED_BANDS propagated to NOVALUE

        MAIN_JUMP, BLUE_JUMP, RED_JUMP, JUMP_BTWN_SIDEBANDS propagated to FITFAILED

        MAIN_INCOMP, BLUE_INCOMP, RED_INCOMP propagated to UNRELIABLE

        NO_ABSORPTION_CORRECTION propaged to NOCORRECTION

        DIVBYZERO propagated to MATHERROR

        Second moments not provided so no need to include UNDEFINED_MOM2
        """

        # Copy the binning mask
        mask = self.bin_mask.copy()

        # Use the emission-line moments map mask to flag low S/N
        # bins/spaxels as NOVALUE
        indx = emission_line_moments.bitmask.flagged(emission_line_moments['MAPMASK'].data,
                                                     flag='LOW_SNR')
        mask[indx] = self.bitmask.turn_on(mask[indx], 'NOVALUE')

        # Reshape the mask to include all emission lines
        mask = numpy.array([mask]*emission_line_moments.nmom).transpose(1,2,0)

        # Consolidate to NOVALUE
        flgd = emission_line_moments.bitmask.flagged(elm_mask,
                                                     flag=['MAIN_EMPTY', 'BLUE_EMPTY', 'RED_EMPTY',
                                                           'UNDEFINED_BANDS' ])
        mask[flgd] = self.bitmask.turn_on(mask[flgd], 'NOVALUE')

        # Consolidate to FITFAILED
        flgd = emission_line_moments.bitmask.flagged(elm_mask, flag=['MAIN_JUMP', 'BLUE_JUMP',
                                                                 'RED_JUMP', 'JUMP_BTWN_SIDEBANDS'])
        mask[flgd] = self.bitmask.turn_on(mask[flgd], 'FITFAILED')

        # Consolidate to UNRELIABLE
        flgd = emission_line_moments.bitmask.flagged(elm_mask, flag=['MAIN_INCOMP', 'BLUE_INCOMP',
                                                                     'RED_INCOMP' ])
        mask[flgd] = self.bitmask.turn_on(mask[flgd], 'UNRELIABLE')

        # Consolidate to NOCORRECTION
        flgd = emission_line_moments.bitmask.flagged(elm_mask, flag=['NO_ABSORPTION_CORRECTION'])
        mask[flgd] = self.bitmask.turn_on(mask[flgd], 'NOCORRECTION')

        # Consolidate to MATHERROR
        flgd = emission_line_moments.bitmask.flagged(elm_mask, flag=['DIVBYZERO'])
        mask[flgd] = self.bitmask.turn_on(mask[flgd], 'MATHERROR')

        return self._consolidate_donotuse(mask)


    def _emission_line_model_mask_to_map_mask(self, emission_line_model, elf_mask,
                                              for_dispersion=False):
        """
        Propagate and consolidate the emission-line moment masks to the map
        pixel mask.

        INSUFFICIENT_DATA propagated to NOVALUE

        FIT_FAILED, UNDEFINED_COVAR propagated to FITFAILED

        NEAR_BOUND copied to NEARBOUND

        UNDEFINED_SIGMA propagated to UNRELIABLE
        """

        # Copy the common mask
        mask = self.bin_mask.copy()

        # Use the emission-line model map mask to flag low S/N
        # bins/spaxels as NOVALUE
        indx = emission_line_model.bitmask.flagged(emission_line_model['MAPMASK'].data,
                                                   flag='LOW_SNR')
        mask[indx] = self.bitmask.turn_on(mask[indx], 'NOVALUE')

        # Reshape the mask to include all emission lines
        mask = numpy.array([mask]*emission_line_model.neml).transpose(1,2,0)

        # Consolidate to NOVALUE
        flgd = emission_line_model.bitmask.flagged(elf_mask, flag=['INSUFFICIENT_DATA'])
        mask[flgd] = self.bitmask.turn_on(mask[flgd], 'NOVALUE')

        # Consolidate to FITFAILED
        flgd = emission_line_model.bitmask.flagged(elf_mask, flag=['FIT_FAILED', 'UNDEFINED_COVAR'])
        mask[flgd] = self.bitmask.turn_on(mask[flgd], 'FITFAILED')

        # Copy NEAR_BOUND
        flgd = emission_line_model.bitmask.flagged(elf_mask, flag='NEAR_BOUND')
        mask[flgd] = self.bitmask.turn_on(mask[flgd], 'NEARBOUND')

        # Consolidate dispersion flags to NEARBOUND and UNRELIABLE
        if for_dispersion:
            flgd = emission_line_model.bitmask.flagged(elf_mask, flag='MIN_SIGMA')
            mask[flgd] = self.bitmask.turn_on(mask[flgd], 'NEARBOUND')

            flgd = emission_line_model.bitmask.flagged(elf_mask, flag=['UNDEFINED_SIGMA',
                                                                       'BAD_SIGMA'])
            mask[flgd] = self.bitmask.turn_on(mask[flgd], 'UNRELIABLE')

        return self._consolidate_donotuse(mask)


    def _spectral_index_mask_to_map_mask(self, spectral_indices, si_mask):
        """
        Propagate and consolidate the spectral index masks to the map
        pixel mask.

        DIDNOTUSE, FORESTAR propagated from already existing mask (self.bin_mask)

        MAIN_EMPTY, BLUE_EMPTY, RED_EMPTY, UNDEFINED_BANDS propagated to NOVALUE

        MAIN_INCOMP, BLUE_INCOMP, RED_INCOMP propagated to UNRELIABLE

        NO_DISPERSION_CORRECTION propagated to NOCORRECTION

        DIVBYZERO propagated to MATHERROR

        Need to further assess \*_INCOMP to see if these lead to bad values (BADVALUE).
        """

        # Copy the common mask
        mask = self.bin_mask.copy()

        # Use the spectral index map mask to flag low S/N bins/spaxels
        # as NOVALUE
        indx = spectral_indices.bitmask.flagged(spectral_indices['MAPMASK'].data, flag='LOW_SNR')
        mask[indx] = self.bitmask.turn_on(mask[indx], 'NOVALUE')

        # Reshape the mask to include all the spectral indices
        mask = numpy.array([mask]*spectral_indices.nindx).transpose(1,2,0)

        # Consolidate to NOVALUE
        flgd = spectral_indices.bitmask.flagged(si_mask, flag=['MAIN_EMPTY', 'BLUE_EMPTY',
                                                               'RED_EMPTY', 'UNDEFINED_BANDS' ])
        mask[flgd] = self.bitmask.turn_on(mask[flgd], 'NOVALUE')

        # Consolidate to UNRELIABLE
        flgd = spectral_indices.bitmask.flagged(si_mask, flag=['MAIN_INCOMP', 'BLUE_INCOMP',
                                                               'RED_INCOMP' ])
        mask[flgd] = self.bitmask.turn_on(mask[flgd], 'UNRELIABLE')

        # Consolidate to NOCORRECTION
        # TODO: What do I do about this!
        flgd = spectral_indices.bitmask.flagged(si_mask, flag=['NO_DISPERSION_CORRECTION'])
        mask[flgd] = self.bitmask.turn_on(mask[flgd], 'NOCORRECTION')

        # Consolidate to MATHERROR
        flgd = spectral_indices.bitmask.flagged(si_mask, flag=['DIVBYZERO'])
        mask[flgd] = self.bitmask.turn_on(mask[flgd], 'MATHERROR')

        return self._consolidate_donotuse(mask)


    def reduction_assessment_maps(self, prihdr, obs, rdxqa):
        """
        Constructs the 'SPX_SKYCOO', 'SPX_ELLCOO', 'SPX_MFLUX',
        'SPX_MFLUX_IVAR', and 'SPX_SNR' map extensions.

        .. todo::
            Add the spaxel correlation data.
        """
        #---------------------------------------------------------------
        # Extensions filled from the reduction assessments
        ext = ['SPX_SKYCOO', 'SPX_ELLCOO', 'SPX_MFLUX', 'SPX_MFLUX_IVAR', 'SPX_SNR']

        if rdxqa is None:
            # Construct and return the empty hdus
            return DAPFitsUtil.empty_hdus(ext)

        #---------------------------------------------------------------
        # Add data to the primary header
        prihdr = rdxqa._initialize_primary_header(hdr=prihdr)

        #---------------------------------------------------------------
        # Get the extension headers
        hdr = [ DAPFitsUtil.finalize_dap_header(self.multichannel_maphdr, 'SPX_SKYCOO',
                                                bunit='arcsec', multichannel=True,
                                                channel_names=['On-sky X', 'On-sky Y']),
                DAPFitsUtil.finalize_dap_header(self.multichannel_maphdr, 'SPX_ELLCOO',
                                                multichannel=True,
                                                channel_names=['Elliptical radius', 'R/Re',
                                                               'Elliptical azimuth'],
                                                channel_units=['arcsec', '', 'degrees']),
                DAPFitsUtil.finalize_dap_header(self.singlechannel_maphdr, 'SPX_MFLUX',
                                                bunit='1E-17 erg/s/cm^2/ang/spaxel', err=True),
                DAPFitsUtil.finalize_dap_header(self.singlechannel_maphdr, 'SPX_MFLUX',
                                                hduclas2='ERROR',
                                                bunit='(1E-17 erg/s/cm^2/ang/spaxel)^{-2}'),
                DAPFitsUtil.finalize_dap_header(self.singlechannel_maphdr, 'SPX_SNR')
              ]

        #---------------------------------------------------------------
        # Get the data arrays
        # On-sky coordinates
        minimum_value = numpy.finfo(self.float_dtype).eps
        spx_skycoo = rdxqa['SPECTRUM'].data['SKY_COO'].copy().reshape(*self.spatial_shape, -1)
        spx_skycoo[numpy.absolute(spx_skycoo) < minimum_value] = 0.0
        # Elliptical coordinates
        spx_ellcoo = rdxqa['SPECTRUM'].data['ELL_COO'].copy().reshape(*self.spatial_shape, -1)
        spx_ellcoo = numpy.repeat(spx_ellcoo, [2,1], axis=2)
        if obs is not None:
            spx_ellcoo[:,:,1] /= obs['reff']
        spx_ellcoo[numpy.absolute(spx_ellcoo) < minimum_value] = 0.0

        # Bin signal
        signal = rdxqa['SPECTRUM'].data['SIGNAL'].copy().reshape(self.spatial_shape)
        signal[numpy.absolute(signal) < minimum_value] = 0.0
        # Bin inverse variance
        ivar = rdxqa['SPECTRUM'].data['VARIANCE'].copy().reshape(self.spatial_shape)
        ivar = numpy.ma.power(ivar, -1).filled(0.0)
        ivar[numpy.absolute(ivar) < minimum_value] = 0.0
        # Bin S/N
        snr = rdxqa['SPECTRUM'].data['SNR'].copy().reshape(self.spatial_shape)
        snr[numpy.absolute(snr) < minimum_value] = 0.0

        # Organize the extension data
        data = [ spx_skycoo.astype(self.float_dtype), spx_ellcoo.astype(self.float_dtype),
                 signal.astype(self.float_dtype), ivar.astype(self.float_dtype),
                 snr.astype(self.float_dtype) ]

        #---------------------------------------------------------------
        # Return the map hdus
        return DAPFitsUtil.list_of_image_hdus(data, hdr, ext)

    
    def binned_spectra_maps(self, prihdr, obs, binned_spectra):
        """
        Constructs the 'BIN_LWSKYCOO', 'BIN_LWELLCOO', 'BIN_AREA',
        'BIN_FAREA', 'BIN_MFLUX', 'BIN_MFLUX_IVAR', 'BIN_MFLUX_MASK',
        and 'BIN_SNR' map extensions.
        """
        #---------------------------------------------------------------
        ext = ['BIN_LWSKYCOO', 'BIN_LWELLCOO', 'BIN_AREA', 'BIN_FAREA', 'BIN_MFLUX',
               'BIN_MFLUX_IVAR', 'BIN_MFLUX_MASK', 'BIN_SNR']

        if binned_spectra is None:
            # Construct and return the empty hdus
            return DAPFitsUtil.empty_hdus(ext)

        #---------------------------------------------------------------
        # Add data to the primary header
        # TODO: Apply this to each extension instead of the primary
        # header to allow for extension specific binning?
        prihdr = binned_spectra._initialize_primary_header(hdr=prihdr)
        if binned_spectra.is_unbinned:
            prihdr['BINTYPE'] = ('None', 'Binning method')
        else:
            prihdr = binned_spectra._add_method_header(prihdr)
        prihdr = binned_spectra._add_reddening_header(prihdr)

        #---------------------------------------------------------------
        # Get the extension headers
        hdr = [ DAPFitsUtil.finalize_dap_header(self.multichannel_maphdr, 'BIN_LWSKYCOO',
                                                bunit='arcsec', multichannel=True,
                                                channel_names=['Lum. weighted on-sky X',
                                                               'Lum. weighted on-sky Y']),
                DAPFitsUtil.finalize_dap_header(self.multichannel_maphdr, 'BIN_LWELLCOO',
                                                multichannel=True,
                                                channel_names= ['Lum. weighted elliptical radius',
                                                                'R/Re',
                                                                'Lum. weighted elliptical azimuth'],
                                                channel_units=['arcsec', '', 'degrees']),
                DAPFitsUtil.finalize_dap_header(self.singlechannel_maphdr, 'BIN_AREA',
                                                bunit='arcsec^2'),
                DAPFitsUtil.finalize_dap_header(self.singlechannel_maphdr, 'BIN_FAREA'),
                DAPFitsUtil.finalize_dap_header(self.singlechannel_maphdr, 'BIN_MFLUX',
                                                bunit='1E-17 erg/s/cm^2/ang/spaxel', err=True,
                                                qual=True),
                DAPFitsUtil.finalize_dap_header(self.singlechannel_maphdr, 'BIN_MFLUX',
                                                hduclas2='ERROR', qual=True,
                                                bunit='(1E-17 erg/s/cm^2/ang/spaxel)^{-2}'),
                DAPFitsUtil.finalize_dap_header(self.singlechannel_maphdr, 'BIN_MFLUX',
                                                hduclas2='QUALITY', err=True,
                                                bit_type=self.bitmask.minimum_dtype()),
                DAPFitsUtil.finalize_dap_header(self.singlechannel_maphdr, 'BIN_SNR')
               ]

        #---------------------------------------------------------------
        # Get the data arrays
        arr = [ binned_spectra['BINS'].data['LW_SKY_COO'][:,0],
                binned_spectra['BINS'].data['LW_SKY_COO'][:,1],
                binned_spectra['BINS'].data['LW_ELL_COO'][:,0],
                binned_spectra['BINS'].data['LW_ELL_COO'][:,0],
                binned_spectra['BINS'].data['LW_ELL_COO'][:,1],
                binned_spectra['BINS'].data['AREA'],
                binned_spectra['BINS'].data['AREA_FRAC'],
                binned_spectra['BINS'].data['SIGNAL'],
                numpy.ma.power(binned_spectra['BINS'].data['VARIANCE'].copy(), -1).filled(0.0),
                binned_spectra['BINS'].data['SNR']
              ]

        dtypes = [self.float_dtype]*len(arr)

        # Bin index
        bin_indx = binned_spectra['BINID'].data.copy().ravel()

        # Remap the data to the DRP spatial shape
        arr = list(DAPFitsUtil.reconstruct_map(self.spatial_shape, bin_indx, arr, dtype=dtypes))

        # Get the normalized radius
        if obs is not None:
            arr[3] = (arr[3]/obs['reff']).astype(self.float_dtype)

        # Organize the extension data
        data = [ numpy.array(arr[0:2]).transpose(1,2,0), numpy.array(arr[2:5]).transpose(1,2,0) ] \
                    + arr[5:-1] + [ self.bin_mask.copy(), arr[-1] ]

        #---------------------------------------------------------------
        # Return the map hdus
        return DAPFitsUtil.list_of_image_hdus(data, hdr, ext)


    def stellar_continuum_maps(self, prihdr, stellar_continuum):
        """
        Construct the 'STELLAR_VEL', 'STELLAR_VEL_IVAR',
        'STELLAR_VEL_MASK', 'STELLAR_SIGMA', 'STELLAR_SIGMA_IVAR',
        'STELLAR_SIGMA_MASK', 'STELLAR_SIGMACORR',
        'STELLAR_CONT_FRESID', and 'STELLAR_CONT_RCHI2' maps extensions.
        """
        #---------------------------------------------------------------
        ext = [ 'STELLAR_VEL', 'STELLAR_VEL_IVAR', 'STELLAR_VEL_MASK', 'STELLAR_SIGMA',
                'STELLAR_SIGMA_IVAR', 'STELLAR_SIGMA_MASK', 'STELLAR_SIGMACORR',
                'STELLAR_CONT_FRESID', 'STELLAR_CONT_RCHI2' ]

        if stellar_continuum is None:
            # Construct and return the empty hdus
            return DAPFitsUtil.empty_hdus(ext)

        #---------------------------------------------------------------
        # Add data to the primary header
        prihdr = stellar_continuum._initialize_primary_header(hdr=prihdr)
        prihdr = stellar_continuum._add_method_header(prihdr)

        #---------------------------------------------------------------
        # Get the extension headers
        hdr = [ DAPFitsUtil.finalize_dap_header(self.singlechannel_maphdr, 'STELLAR_VEL',
                                                bunit='km/s', err=True, qual=True),
                DAPFitsUtil.finalize_dap_header(self.singlechannel_maphdr, 'STELLAR_VEL',
                                                hduclas2='ERROR', bunit='(km/s)^{-2}', qual=True),
                DAPFitsUtil.finalize_dap_header(self.singlechannel_maphdr, 'STELLAR_VEL',
                                                hduclas2='QUALITY', err=True,
                                                bit_type=self.bitmask.minimum_dtype()),
                DAPFitsUtil.finalize_dap_header(self.singlechannel_maphdr, 'STELLAR_SIGMA',
                                                bunit='km/s', err=True, qual=True),
                DAPFitsUtil.finalize_dap_header(self.singlechannel_maphdr, 'STELLAR_SIGMA',
                                                hduclas2='ERROR', bunit='(km/s)^{-2}', qual=True),
                DAPFitsUtil.finalize_dap_header(self.singlechannel_maphdr, 'STELLAR_SIGMA',
                                                hduclas2='QUALITY', err=True,
                                                bit_type=self.bitmask.minimum_dtype()),
                DAPFitsUtil.finalize_dap_header(self.singlechannel_maphdr, 'STELLAR_SIGMACORR',
                                                bunit='km/s'),
#                DAPFitsUtil.finalize_dap_header(self.multichannel_maphdr, 'STELLAR_SIGMACORR',
#                                                multichannel=True, bunit='km/s',
#                                                channel_names=['resolution difference', 'fit']),
                DAPFitsUtil.finalize_dap_header(self.multichannel_maphdr, 'STELLAR_CONT_FRESID',
                                                multichannel=True,
                                                channel_names=['68th percentile',
                                                               '99th percentile']),
                DAPFitsUtil.finalize_dap_header(self.singlechannel_maphdr, 'STELLAR_CONT_RCHI2')
              ]


        #---------------------------------------------------------------
        # Get the data arrays
        arr = [ DAPFitsUtil.redshift_to_Newtonian_velocity(
                                    stellar_continuum['PAR'].data['KIN'][:,0], self.nsa_redshift),
                DAPFitsUtil.redshift_to_Newtonian_velocity(numpy.ma.power(
                                                    stellar_continuum['PAR'].data['KINERR'][:,0],
                                                    -2).filled(0.0), self.nsa_redshift, ivar=True),
                stellar_continuum['PAR'].data['KIN'][:,1],
                numpy.ma.power(stellar_continuum['PAR'].data['KINERR'][:,1], -2).filled(0.0),
# Only provide the empirical correction
#                stellar_continuum['PAR'].data['SIGMACORR_SRES'],
                stellar_continuum['PAR'].data['SIGMACORR_EMP'],
                stellar_continuum['PAR'].data['FABSRESID'][:,1],
                stellar_continuum['PAR'].data['FABSRESID'][:,3],
                stellar_continuum['PAR'].data['RCHI2'],
                stellar_continuum['PAR'].data['MASK']
              ]

        dtypes = [self.float_dtype]*(len(arr)-1) + [arr[-1].dtype.name]

        # Bin index
        bin_indx = stellar_continuum['BINID'].data.copy().ravel()
##        pyplot.imshow(bin_indx.reshape(self.drpf.spatial_shape), origin='lower', interpolation='nearest')
##        pyplot.show()
#
#        pyplot.imshow(DAPFitsUtil.reconstruct_map(self.spatial_shape, bin_indx, stellar_continuum['PAR'].data['KIN'][:,0]), origin='lower', interpolation='nearest')
#        pyplot.show()
#
#        print(stellar_continuum['PAR'].data['KIN'][:,0])
#
#        print('nbin: {0}'.format(numpy.sum(bin_indx > -1)))    
#        exit()

        # Remap the data to the DRP spatial shape
        arr = list(DAPFitsUtil.reconstruct_map(self.spatial_shape, bin_indx, arr, dtype=dtypes))
        
        # Get the masks
        vel_mask = self._stellar_continuum_mask_to_map_mask(stellar_continuum, arr[-1].copy())
        sig_mask = self._stellar_continuum_mask_to_map_mask(stellar_continuum, arr[-1].copy(),
                                                            for_dispersion=True)

        # Organize the extension data
        data = arr[0:2] + [ vel_mask ] + arr[2:4] + [ sig_mask ] \
                    + [ arr[4], numpy.array(arr[5:7]).transpose(1,2,0), arr[7] ]
#        data = arr[0:2] + [ vel_mask ] + arr[2:4] + [ sig_mask ] \
#                    + [ numpy.array(arr[4:6]).transpose(1,2,0),
#                        numpy.array(arr[6:8]).transpose(1,2,0), arr[8] ]

#        for i,d in enumerate(data):
#            if len(d.shape) == 2:
#                print(i+1)
#                pyplot.imshow(d, origin='lower', interpolation='nearest')
#                pyplot.show()
#        exit()

#        print(data, hdr, ext)
#        print(len(data), len(hdr), len(ext))
#        exit()

        #---------------------------------------------------------------
        # Return the map hdus
        return DAPFitsUtil.list_of_image_hdus(data, hdr, ext)


    def emission_line_moment_maps(self, prihdr, emission_line_moments):
        """
        Construct the 'EMLINE_SFLUX', 'EMLINE_SFLUX_IVAR',
        'EMLINE_SFLUX_MASK', 'EMLINE_SEW', 'EMLINE_SEW_IVAR', and
        'EMLINE_SEW_MASK' maps extensions.
        """
        #---------------------------------------------------------------
        ext = [ 'EMLINE_SFLUX', 'EMLINE_SFLUX_IVAR', 'EMLINE_SFLUX_MASK', 'EMLINE_SEW',
                'EMLINE_SEW_IVAR', 'EMLINE_SEW_MASK' ]

        if emission_line_moments is None:
            # Construct and return the empty hdus
            return DAPFitsUtil.empty_hdus(ext)

        #---------------------------------------------------------------
        # Add data to the primary header
        prihdr = emission_line_moments._initialize_primary_header(hdr=prihdr)

        #---------------------------------------------------------------
        # Get the extension headers

        # Need multichannel map header if more than one moment
        multichannel = emission_line_moments.nmom > 1
        # Build the channel names
#        names = [ '{0}-{1}'.format(n,int(w)) \
#                        for n,w in zip(emission_line_moments['ELMBAND'].data['NAME'],
#                                       emission_line_moments['ELMBAND'].data['RESTWAVE']) ]
        names = emission_line_moments.channel_names()
        # Create the basic header for all extensions
        base_hdr = DAPFitsUtil.add_channel_names(self.multichannel_maphdr if multichannel
                                                 else self.singlechannel_maphdr, names)

        hdr = [ DAPFitsUtil.finalize_dap_header(base_hdr, 'EMLINE_SFLUX', err=True, qual=True,
                                                bunit='1E-17 erg/s/cm^2/spaxel',
                                                multichannel=multichannel),
                DAPFitsUtil.finalize_dap_header(base_hdr, 'EMLINE_SFLUX', hduclas2='ERROR',
                                                qual=True, bunit='(1E-17 erg/s/cm^2/spaxel)^{-2}',
                                                multichannel=multichannel),
                DAPFitsUtil.finalize_dap_header(base_hdr, 'EMLINE_SFLUX', hduclas2='QUALITY',
                                                err=True, bit_type=self.bitmask.minimum_dtype(),
                                                multichannel=multichannel),
                DAPFitsUtil.finalize_dap_header(base_hdr, 'EMLINE_SEW', err=True, qual=True,
                                                bunit='ang', multichannel=multichannel),
                DAPFitsUtil.finalize_dap_header(base_hdr, 'EMLINE_SEW', hduclas2='ERROR', qual=True,
                                                bunit='(ang)^{-2}', multichannel=multichannel),
                DAPFitsUtil.finalize_dap_header(base_hdr, 'EMLINE_SEW', hduclas2='QUALITY',
                                                err=True, bit_type=self.bitmask.minimum_dtype(),
                                                multichannel=multichannel)
              ]

        #---------------------------------------------------------------
        # Get the data arrays
        arr = [ emission_line_moments['ELMMNTS'].data['FLUX'][:,m]
                    for m in range(emission_line_moments.nmom) ]
        arr += [ numpy.ma.power(emission_line_moments['ELMMNTS'].data['FLUXERR'][:,m],
                                -2).filled(0.0) for m in range(emission_line_moments.nmom) ]
        arr += [ emission_line_moments['ELMMNTS'].data['EW'][:,m]
                    for m in range(emission_line_moments.nmom) ]
        arr += [ numpy.ma.power(emission_line_moments['ELMMNTS'].data['EWERR'][:,m],
                                -2).filled(0.0) for m in range(emission_line_moments.nmom) ]
        arr += [ emission_line_moments['ELMMNTS'].data['MASK'][:,m]
                    for m in range(emission_line_moments.nmom) ]

        dtypes = [self.float_dtype]*(len(arr)-emission_line_moments.nmom) \
                        + [a.dtype.name for a in arr[-emission_line_moments.nmom:]]

        # Bin index
        bin_indx = emission_line_moments['BINID'].data.copy().ravel()

        # Remap the data to the DRP spatial shape
        arr = list(DAPFitsUtil.reconstruct_map(self.spatial_shape, bin_indx, arr, dtype=dtypes))

        data = [ numpy.array(arr[emission_line_moments.nmom*i:
                                 emission_line_moments.nmom*(i+1)]).transpose(1,2,0) 
                        for i in range(5) ]
        
        # Get the mask
        elm_mask = self._emission_line_moment_mask_to_map_mask(emission_line_moments,
                                                               data[-1].copy())

        # Organize the extension data
        data = data[:2] + [elm_mask] + data[2:-1] + [elm_mask]

        #---------------------------------------------------------------
        # Return the map hdus
        return DAPFitsUtil.list_of_image_hdus(data, hdr, ext)


    def emission_line_model_maps(self, prihdr, emission_line_model):
        """
        Construct the 'EMLINE_GFLUX', 'EMLINE_GFLUX_IVAR',
        'EMLINE_GFLUX_MASK', 'EMLINE_GEW', 'EMLINE_GEW_IVAR',
        'EMLINE_GEW_MASK', 'EMLINE_GVEL', 'EMLINE_GVEL_IVAR',
        'EMLINE_GVEL_MASK', 'EMLINE_GSIGMA', 'EMLINE_GSIGMA_IVAR',
        'EMLINE_GSIGMA_MASK', 'EMLINE_INSTSIGMA', 'EMLINE_TPLSIGMA' map extensions.
        """
        #---------------------------------------------------------------
        ext = [ 'EMLINE_GFLUX', 'EMLINE_GFLUX_IVAR', 'EMLINE_GFLUX_MASK', 'EMLINE_GEW',
                'EMLINE_GEW_IVAR', 'EMLINE_GEW_MASK', 'EMLINE_GVEL', 'EMLINE_GVEL_IVAR',
                'EMLINE_GVEL_MASK', 'EMLINE_GSIGMA', 'EMLINE_GSIGMA_IVAR', 'EMLINE_GSIGMA_MASK',
                'EMLINE_INSTSIGMA', 'EMLINE_TPLSIGMA' ]

        if emission_line_model is None:
            # Construct and return the empty hdus
            return DAPFitsUtil.empty_hdus(ext)

        #---------------------------------------------------------------
        # Add data to the primary header
        prihdr = emission_line_model._initialize_primary_header(hdr=prihdr)
        prihdr = emission_line_model._add_method_header(prihdr)

        #---------------------------------------------------------------
        # Get the extension headers

        # Need multichannel map header if more than one moment
        multichannel = emission_line_model.neml > 1

        # Create the basic header for all extensions
        if 'data' in emission_line_model['PAR'].__dict__:
            names= [ '{0}-{1}'.format(n,int(w)) \
                        for n,w in zip(emission_line_model['PAR'].data['NAME'],
                                       emission_line_model['PAR'].data['RESTWAVE']) ]
            base_hdr = DAPFitsUtil.add_channel_names(self.multichannel_maphdr if multichannel
                                                     else self.singlechannel_maphdr, names)
        else:
            # TODO: This should throw a ValueError, not just a warning
            warnings.warn('Emission-line model does not include channel names!')
            base_hdr = self.multichannel_maphdr if multichannel else self.singlechannel_maphdr
       
        hdr = [ DAPFitsUtil.finalize_dap_header(base_hdr, 'EMLINE_GFLUX', err=True, qual=True,
                                                bunit='1E-17 erg/s/cm^2/spaxel',
                                                multichannel=multichannel),
                DAPFitsUtil.finalize_dap_header(base_hdr, 'EMLINE_GFLUX', hduclas2='ERROR',
                                                qual=True, bunit='(1E-17 erg/s/cm^2/spaxel)^{-2}',
                                                multichannel=multichannel),
                DAPFitsUtil.finalize_dap_header(base_hdr, 'EMLINE_GFLUX', hduclas2='QUALITY',
                                                err=True, bit_type=self.bitmask.minimum_dtype(),
                                                multichannel=multichannel),
                DAPFitsUtil.finalize_dap_header(base_hdr, 'EMLINE_GEW', err=True, qual=True,
                                                bunit='ang', multichannel=multichannel),
                DAPFitsUtil.finalize_dap_header(base_hdr, 'EMLINE_GEW', hduclas2='ERROR', qual=True,
                                                bunit='(ang)^{-2}', multichannel=multichannel),
                DAPFitsUtil.finalize_dap_header(base_hdr, 'EMLINE_GEW', hduclas2='QUALITY',
                                                err=True, bit_type=self.bitmask.minimum_dtype(),
                                                multichannel=multichannel),
                DAPFitsUtil.finalize_dap_header(base_hdr, 'EMLINE_GVEL', err=True, qual=True,
                                                bunit='km/s', multichannel=multichannel),
                DAPFitsUtil.finalize_dap_header(base_hdr, 'EMLINE_GVEL', hduclas2='ERROR',
                                                qual=True, bunit='(km/s)^{-2}',
                                                multichannel=multichannel),
                DAPFitsUtil.finalize_dap_header(base_hdr, 'EMLINE_GVEL', hduclas2='QUALITY',
                                                err=True, bit_type=self.bitmask.minimum_dtype(),
                                                multichannel=multichannel),
                DAPFitsUtil.finalize_dap_header(base_hdr, 'EMLINE_GSIGMA', err=True, qual=True,
                                                bunit='km/s', multichannel=multichannel),
                DAPFitsUtil.finalize_dap_header(base_hdr, 'EMLINE_GSIGMA', hduclas2='ERROR',
                                                qual=True, bunit='(km/s)^{-2}',
                                                multichannel=multichannel),
                DAPFitsUtil.finalize_dap_header(base_hdr, 'EMLINE_GSIGMA', hduclas2='QUALITY',
                                                err=True, bit_type=self.bitmask.minimum_dtype(),
                                                multichannel=multichannel),
                DAPFitsUtil.finalize_dap_header(base_hdr, 'EMLINE_INSTSIGMA', bunit='km/s',
                                                multichannel=multichannel),
                DAPFitsUtil.finalize_dap_header(base_hdr, 'EMLINE_TPLSIGMA', bunit='km/s',
                                                multichannel=multichannel)
              ]

        #---------------------------------------------------------------
        # Get the data arrays
        narr = 11
        arr = [ emission_line_model['EMLDATA'].data['FLUX'][:,m]
                    for m in range(emission_line_model.neml) ]
        arr += [ numpy.ma.power(emission_line_model['EMLDATA'].data['FLUXERR'][:,m],
                                -2).filled(0.0) for m in range(emission_line_model.neml) ]
        arr += [ emission_line_model['EMLDATA'].data['EW'][:,m]
                    for m in range(emission_line_model.neml) ]
        arr += [ numpy.ma.power(emission_line_model['EMLDATA'].data['EWERR'][:,m],
                                -2).filled(0.0) for m in range(emission_line_model.neml) ]
        arr += [ DAPFitsUtil.redshift_to_Newtonian_velocity(
                        emission_line_model['EMLDATA'].data['KIN'][:,m,0], self.nsa_redshift)
                    for m in range(emission_line_model.neml) ]
        arr += [ DAPFitsUtil.redshift_to_Newtonian_velocity(
                        numpy.ma.power(emission_line_model['EMLDATA'].data['KINERR'][:,m,0],
                                       -2).filled(0.0), self.nsa_redshift, ivar=True) \
                    for m in range(emission_line_model.neml) ]
        arr += [ emission_line_model['EMLDATA'].data['KIN'][:,m,1]
                    for m in range(emission_line_model.neml) ]
        arr += [ numpy.ma.power(emission_line_model['EMLDATA'].data['KINERR'][:,m,1],
                                -2).filled(0.0) for m in range(emission_line_model.neml) ]
        arr += [ emission_line_model['EMLDATA'].data['SIGMAINST'][:,m]
                    for m in range(emission_line_model.neml) ]
        arr += [ emission_line_model['EMLDATA'].data['SIGMATPL'][:,m]
                    for m in range(emission_line_model.neml) ]
        arr += [ emission_line_model['EMLDATA'].data['MASK'][:,m]
                    for m in range(emission_line_model.neml) ]

        dtypes = [self.float_dtype]*(len(arr)-emission_line_model.neml) \
                        + [a.dtype.name for a in arr[-emission_line_model.neml:]]

##        t = [ numpy.ma.power(emission_line_model['EMLDATA'].data['KINERR'][:,m,1],
##                                -2).filled(0.0) for m in range(emission_line_model.neml) ]
#        t = [ emission_line_model['EMLDATA'].data['KINERR'][:,m,1]
#                                 for m in range(emission_line_model.neml) ]
#        for k,_t in enumerate(t):
##            print(k, numpy.sum(_t > 1e20))
#            print(k, numpy.amin(_t))
#
##        print('len(arr): ', len(arr))
##        print('not finite arr: ', [numpy.sum(numpy.invert(numpy.isfinite(a))) for a in arr])

        # Bin index
        bin_indx = emission_line_model['BINID'].data.copy().ravel()

        # Remap the data to the DRP spatial shape
        arr = list(DAPFitsUtil.reconstruct_map(self.spatial_shape, bin_indx, arr, dtype=dtypes))

        data = [ numpy.array(arr[emission_line_model.neml*i:
                                 emission_line_model.neml*(i+1)]).transpose(1,2,0) \
                        for i in range(narr) ]

#        print(len(data))
#        print(data[7].shape)
#        print([numpy.sum(numpy.invert(numpy.isfinite(d))) for d in data])

        # Get the masks
        base_mask = self._emission_line_model_mask_to_map_mask(emission_line_model, data[-1].copy())
        sig_mask = self._emission_line_model_mask_to_map_mask(emission_line_model, data[-1].copy(),
                                                              for_dispersion=True)

        # Organize the extension data
        data = data[:2] + [ base_mask ] + data[2:4] + [ base_mask ] + data[4:6] + [ base_mask ] \
                + data[6:8] + [ sig_mask ] + data[8:10]

        #---------------------------------------------------------------
        # Return the map hdus
        return DAPFitsUtil.list_of_image_hdus(data, hdr, ext)


    def spectral_index_maps(self, prihdr, spectral_indices):
        """
        Construct the 'SPECINDEX', 'SPECINDEX_IVAR', 'SPECINDEX_MASK',
        and 'SPECINDEX_CORR'.
        """
        #---------------------------------------------------------------
        ext = [ 'SPECINDEX', 'SPECINDEX_IVAR', 'SPECINDEX_MASK', 'SPECINDEX_CORR' ]

        if spectral_indices is None:
            # Construct and return the empty hdus
            return DAPFitsUtil.empty_hdus(ext)

        #---------------------------------------------------------------
        # Add data to the primary header
        prihdr = spectral_indices._initialize_primary_header(hdr=prihdr)

        #---------------------------------------------------------------
        # Get the extension headers

        # Need multichannel map header if more than one index
        multichannel = spectral_indices.nindx > 1
        # Create the basic header for all extensions
        base_hdr = self.multichannel_maphdr if multichannel else self.singlechannel_maphdr

        errunits = [ '({0})'.format(u)+'^{-2}' if len(u) > 0 else ''
                            for u in spectral_indices['SIPAR'].data['UNIT'] ]

        hdr = [ DAPFitsUtil.finalize_dap_header(base_hdr, 'SPECINDEX', err=True, qual=True,
                                                multichannel=multichannel,
                                              channel_names=spectral_indices['SIPAR'].data['NAME'],
                                              channel_units=spectral_indices['SIPAR'].data['UNIT']),
                DAPFitsUtil.finalize_dap_header(base_hdr, 'SPECINDEX', hduclas2='ERROR', qual=True,
                                                multichannel=multichannel,
                                              channel_names=spectral_indices['SIPAR'].data['NAME'],
                                                channel_units=errunits),
                DAPFitsUtil.finalize_dap_header(base_hdr, 'SPECINDEX', hduclas2='QUALITY', err=True,
                                                bit_type=self.bitmask.minimum_dtype(),
                                                multichannel=multichannel,
                                              channel_names=spectral_indices['SIPAR'].data['NAME']),
                DAPFitsUtil.finalize_dap_header(base_hdr, 'SPECINDEX_CORR',
                                                multichannel=multichannel,
                                               channel_names=spectral_indices['SIPAR'].data['NAME'])
              ]

        #---------------------------------------------------------------
        # Get the data arrays
        arr = [ spectral_indices['SINDX'].data['INDX'][:,m]
                    for m in range(spectral_indices.nindx) ]
        arr += [ numpy.ma.power(spectral_indices['SINDX'].data['INDXERR'][:,m], -2.).filled(0.0)
                    for m in range(spectral_indices.nindx) ]
        arr += [ spectral_indices['SINDX'].data['MASK'][:,m]
                    for m in range(spectral_indices.nindx) ]
        arr += [ spectral_indices['SINDX'].data['INDX_DISPCORR'][:,m]
                    for m in range(spectral_indices.nindx) ]

        dtypes = [self.float_dtype]*(len(arr)-2*spectral_indices.nindx) \
                        + [a.dtype.name 
                            for a in arr[-2*spectral_indices.nindx:-spectral_indices.nindx]] \
                        + [self.float_dtype]*spectral_indices.nindx

        # Bin index
        bin_indx = spectral_indices['BINID'].data.copy().ravel()

        # Remap the data to the DRP spatial shape
        arr = list(DAPFitsUtil.reconstruct_map(self.spatial_shape, bin_indx, arr, dtype=dtypes))

        data = [ numpy.array(arr[spectral_indices.nindx*i:
                                 spectral_indices.nindx*(i+1)]).transpose(1,2,0) \
                        for i in range(4) ]
        
        # Get the mask
        data[2] = self._spectral_index_mask_to_map_mask(spectral_indices, data[2].copy())

        #---------------------------------------------------------------
        # Return the map hdus
        return DAPFitsUtil.list_of_image_hdus(data, hdr, ext)




    

#-----------------------------------------------------------------------
class construct_cube_file:
    """
    Construct a DAP cube file based on the input data.

    Set as a class for coherency reasons, but should not be used as an
    object!

    Should force all intermediate objects to be provided.

    """
    def __init__(self, drpf, obs=None, binned_spectra=None, stellar_continuum=None,
                 emission_line_model=None, dapsrc=None, dapver=None, analysis_path=None,
                 directory_path=None, output_file=None, clobber=True, loggers=None, quiet=False,
                 single_precision=False):

        #---------------------------------------------------------------
        # Initialize the reporting
        self.loggers = None if loggers is None else loggers
        self.quiet = quiet

        self.float_dtype = 'float32' if single_precision else 'float'

        #---------------------------------------------------------------
        # Check input types
        confirm_dap_types(drpf, None, None, binned_spectra, stellar_continuum, None,
                          emission_line_model, None)

        #---------------------------------------------------------------
        # Set the output paths
        self.drpf = drpf
        self.method = None
        self.directory_path = None
        self.output_file = None
        self._set_paths(directory_path, dapver, analysis_path, output_file, binned_spectra,
                        stellar_continuum)

        # Save input for reference
        self.shape = self.drpf.shape
        self.spatial_shape = self.drpf.spatial_shape
#        self.cube_arrays = None
#        self._assign_cube_arrays()

        # Report
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, '-'*50)
            log_output(loggers, 1, logging.INFO, '{0:^50}'.format('CONSTRUCTING OUTPUT MODEL CUBE'))
            log_output(self.loggers, 1, logging.INFO, '-'*50)
            log_output(self.loggers, 1, logging.INFO, 'Output path: {0}'.format(
                                                                            self.directory_path))
            log_output(self.loggers, 1, logging.INFO, 'Output file: {0}'.format(
                                                                            self.output_file))
            log_output(self.loggers, 1, logging.INFO, 'Output cubes have shape {0}'.format(
                                                                            self.shape))

        #---------------------------------------------------------------
        # Check if the file already exists
        ofile = os.path.join(self.directory_path, self.output_file)
        if os.path.isfile(ofile) and not clobber:
            # TODO: Perform some checks to make sure the existing file
            # has the correct content?
            warnings.warn('Output file exists!  Set clobber=True to overwrite.')
            return

        #---------------------------------------------------------------
        # Construct the 3D data arrays
        binned_spectra_3d_hdu = None if binned_spectra is None \
                                        else binned_spectra.construct_3d_hdu()
        stellar_continuum_3d_hdu = None if stellar_continuum is None \
                                        else stellar_continuum.construct_3d_hdu()
        emission_line_model_3d_hdu = None if emission_line_model is None \
                                        else emission_line_model.construct_3d_hdu()

        #---------------------------------------------------------------
        # Initialize the primary header
        prihdr = DAPFitsUtil.initialize_dap_primary_header(self.drpf, maskname='MANGA_DAPSPECMASK')
        # Add the DAP method
        prihdr['DAPTYPE'] = (default_dap_method(binned_spectra=binned_spectra,
                                                stellar_continuum=stellar_continuum),
                             'DAP analysis method')
        # Add the format of this file
        prihdr['DAPFRMT'] = ('LOGCUBE', 'DAP data file format')

        # Get the base map header
        self.base_cubehdr = DAPFitsUtil.build_cube_header(self.drpf,
                                'K Westfall & B Andrews <westfall@ucolick.org, andrewsb@pitt.edu>',
                                                        maskname='MANGA_DAPSPECMASK')

        #---------------------------------------------------------------
        # Initialize the pixel mask
        self.bitmask = DAPCubeBitMask(dapsrc=dapsrc)

        self.bin_mask = self._initialize_mask(binned_spectra, binned_spectra_3d_hdu)

        #---------------------------------------------------------------
        # Construct the hdu list for each input object.
        # Binned and model spectra:
        prihdr, speclist = self.model_and_data_cube(prihdr, binned_spectra, binned_spectra_3d_hdu,
                                                    stellar_continuum, stellar_continuum_3d_hdu,
                                                    emission_line_model, emission_line_model_3d_hdu)
        # Emission-line only models:
        prihdr, elmodlist = self.emission_line_model_cube(prihdr, binned_spectra,
                                                          emission_line_model,
                                                          emission_line_model_3d_hdu)

        # Get the BINIDs
        binidlist = combine_binid_extensions(self.drpf, binned_spectra, stellar_continuum, None,
                                             emission_line_model, None, dtype='int32')

        #---------------------------------------------------------------
        # Save the data to the hdu attribute
        prihdr = finalize_dap_primary_header(prihdr, self.drpf, obs, binned_spectra,
                                             stellar_continuum, dapsrc=dapsrc,
                                             loggers=self.loggers, quiet=self.quiet)

        self.hdu = fits.HDUList([ fits.PrimaryHDU(header=prihdr),
                                  *speclist,
                                  *elmodlist,
                                  *binidlist
                                ])

        #---------------------------------------------------------------

        # Check that the path exists
        if not os.path.isdir(self.directory_path):
            os.makedirs(self.directory_path)
        # Write the model cube file
        DAPFitsUtil.write(self.hdu, ofile, clobber=clobber, checksum=True, loggers=self.loggers,
                          quiet=self.quiet)
        # End
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, '-'*50)

#        # Restructure the cubes before writing
#        DAPFitsUtil.restructure_cube(self.hdu, ext=self.cube_arrays, inverse=True)
#        # Write the file
#        write_hdu(self.hdu, ofile, clobber=clobber, checksum=True, loggers=self.loggers,
#                  quiet=self.quiet)


    def _set_paths(self, directory_path, dapver, analysis_path, output_file, binned_spectra,
                   stellar_continuum):
        """
        Set the paths relevant to the map file.
        """
        # The output method directory is, for now, the combination of
        # the binned_spectrum and stellar_continuum method keys
        if directory_path is None and (binned_spectra is None or stellar_continuum is None):
            raise ValueError('Could not define output directory path.')

        # Set the output directory path
        # TODO: Get DAP version from __version__ string
        self.method = default_dap_method(binned_spectra=binned_spectra,
                                         stellar_continuum=stellar_continuum) \
                                if directory_path is None else None
        self.directory_path = default_dap_method_path(self.method, plate=self.drpf.plate,
                                                      ifudesign=self.drpf.ifudesign,
                                                      drpver=self.drpf.drpver, dapver=dapver,
                                                      analysis_path=analysis_path) \
                                                if directory_path is None else str(directory_path)

        # Set the output file
        self.output_file = default_dap_file_name(self.drpf.plate, self.drpf.ifudesign,
                                                 self.method, mode='LOGCUBE') \
                                    if output_file is None else str(output_file)


    def _initialize_mask(self, binned_spectra, binned_spectra_3d_hdu):
        """
        Initialize the mask based on the DRP cube mask.

        IVAR_INVALID:  Test of invalid inverse variance based on the DRP
        file.  THIS SHOULD BE A TEST THAT IS INCLUDED IN THE
        BINNED_SPECTRA OBJECT.  This should match the flags from the
        stellar continuum model object.

        Copy the FORESTAR flags from the DRP file to the this one
        """
        mask = numpy.zeros(self.shape, dtype=self.bitmask.minimum_dtype())
        if binned_spectra is None:
            return mask

        indx = numpy.invert(binned_spectra_3d_hdu['IVAR'].data > 0) \
                    | numpy.invert(numpy.isfinite(binned_spectra_3d_hdu['IVAR'].data))
        mask[indx] = self.bitmask.turn_on(mask[indx], 'IVARINVALID')

        indx = binned_spectra.bitmask.flagged(binned_spectra_3d_hdu['MASK'].data, flag='FORESTAR')
        mask[indx] = self.bitmask.turn_on(mask[indx], 'FORESTAR')

        return mask


#    def _include_stellar_continuum_mask(self, mask, stellar_continuum):
#        """
#        From the stellar continuum model object:
#            - copy ARTIFACT
#            - copy INVALID_ERROR into IVAR_INVALID
#            - consolidate DIDNOTUSE, LOW_SNR, OUTSIDE_RANGE, EML_REGION,
#              TPL_PIXELS, TRUNCATED, PPXF_REJECT, INVALID_ERROR into SC_IGNORED
#            - copy FIT_FAILED, NEAR_BOUND into SC_FAILED
#        """
#        indx = self.bitmask.flagged(self.drpf['MASK'].data, flag='ARTIFACT')
#        mask[indx] = self.bitmask.turn_on(mask[indx], 'ARTIFACT')
#
#        indx = self.bitmask.flagged(self.drpf['MASK'].data, flag='INVALID_ERROR')
#        mask[indx] = self.bitmask.turn_on(mask[indx], 'IVAR_INVALID')
#
#        flags = [ 'DIDNOTUSE', 'LOW_SNR', 'OUTSIDE_RANGE', 'EML_REGION', 'TPL_PIXELS', 'TRUNCATED',
#                  'PPXF_REJECT', 'INVALID_ERROR' ]
#        indx = self.bitmask.flagged(self.drpf['MASK'].data, flag=flags)
#        mask[indx] = self.bitmask.turn_on(mask[indx], 'SC_IGNORED')
#
#        indx = self.bitmask.flagged(self.drpf['MASK'].data, flag='FIT_FAILED')
#        mask[indx] = self.bitmask.turn_on(mask[indx], 'SC_FAILED')
#
#        return mask


    def _assign_cube_arrays(self):
        """
        Set :attr:`cube_arrays`, which contains the list of extensions in
        :attr:`hdu` that contain spectral data.
        """
        self.cube_arrays = [ 'FLUX', 'IVAR', 'MASK', 'MODEL', 'EMLINE', 'EMLINE_BASE',
                             'EMLINE_MASK', 'BINID' ]


    def _get_model_and_data_mask(self, binned_spectra, binned_spectra_3d_hdu, stellar_continuum,
                                 stellar_continuum_3d_hdu, emission_line_model,
                                 emission_line_model_3d_hdu):
        """
        For the binned spectra:
            - consolidate DIDNOTUSE, LOW_SPECCOV, LOW_SNR from
              binned_spectra into IGNORED
            - consolidate NONE_IN_STACK from binned_spectra into
              FLUXINVALID

        For the model spectra:
            - copy INVALID_ERROR from stellar_continuum into
              IVARINVALID

            - copy ARTIFACTs from both stellar_continuum and
              emission_line_model

            - from stellar_continuum, consolidate DIDNOTUSE, LOW_SNR,
              OUTSIDE_RANGE, EML_REGION, TPL_PIXELS, TRUNCATED,
              PPXF_REJECT, INVALID_ERROR into a list of pixels ignored
              by the stellar-continuum fit
            - from emission_line_model, consolidate DIDNOTUSE, LOW_SNR,
              OUTSIDE_RANGE into a list of pixels ignored by the
              emission-line fit
            - flag pixels as FITIGNORED if the pixel is ignored by
              **both** the stellar-continuum and emission-line fits

            - from stellar_continuum, consolidate FIT_FAILED and
              NEAR_BOUND into a list of failed pixels from the
              stellar-continuum fit
            - do the same for the emission-line fitting mask
            - flag pixels as FIT_FAILED if the pixel failed in
              **either** the stellar-continuum or emission-line fits
        """
        # Copy the base-level mask
        mask = self.bin_mask.copy()

        # Draw from binned_spectra
        flags = [ 'DIDNOTUSE', 'LOW_SPECCOV', 'LOW_SNR' ]
        indx = binned_spectra.bitmask.flagged(binned_spectra_3d_hdu['MASK'].data, flag=flags)
        mask[indx] = self.bitmask.turn_on(mask[indx], 'IGNORED')

        indx = binned_spectra.bitmask.flagged(binned_spectra_3d_hdu['MASK'].data,
                                              flag='NONE_IN_STACK')
        mask[indx] = self.bitmask.turn_on(mask[indx], 'FLUXINVALID')

        # Draw from stellar_continuum
        indx = stellar_continuum.bitmask.flagged(stellar_continuum_3d_hdu['MASK'].data,
                                                 flag='INVALID_ERROR')
        mask[indx] = self.bitmask.turn_on(mask[indx], 'IVARINVALID')

        # Draw from stellar_continuum and emission_line_model
        indx = stellar_continuum.bitmask.flagged(stellar_continuum_3d_hdu['MASK'].data,
                                                 flag='ARTIFACT') \
                    | emission_line_model.bitmask.flagged(emission_line_model_3d_hdu['MASK'].data,
                                                          flag='ARTIFACT')
        mask[indx] = self.bitmask.turn_on(mask[indx], 'ARTIFACT')

        flags = [ 'DIDNOTUSE', 'LOW_SNR', 'OUTSIDE_RANGE', 'EML_REGION', 'TPL_PIXELS',
                  'TRUNCATED', 'PPXF_REJECT', 'INVALID_ERROR' ]
        sc_indx = stellar_continuum.bitmask.flagged(stellar_continuum_3d_hdu['MASK'].data,
                                                    flag=flags)
        flags = [ 'DIDNOTUSE', 'LOW_SNR', 'OUTSIDE_RANGE', 'TPL_PIXELS',
                  'TRUNCATED', 'PPXF_REJECT', 'INVALID_ERROR' ]
        el_indx = emission_line_model.bitmask.flagged(emission_line_model_3d_hdu['MASK'].data,
                                                      flag=flags)

        mask[sc_indx & el_indx] = self.bitmask.turn_on(mask[sc_indx & el_indx], 'FITIGNORED')

#        pyplot.imshow(numpy.log10(mask[binned_spectra['BINID'].data > -1,:].reshape(-1,4563)),
#                        interpolation='nearest', aspect='auto')
#        pyplot.show()

        # TODO: What do I do with BAD_SIGMA?
        sc_indx = stellar_continuum.bitmask.flagged(stellar_continuum_3d_hdu['MASK'].data,
                                                    flag=['FIT_FAILED', 'NEAR_BOUND'])
        el_indx = emission_line_model.bitmask.flagged(emission_line_model_3d_hdu['MASK'].data,
                                                      flag=['FIT_FAILED', 'NEAR_BOUND'])
        mask[sc_indx | el_indx] = self.bitmask.turn_on(mask[sc_indx | el_indx], 'FITFAILED')

#        pyplot.imshow(numpy.log10(mask[binned_spectra['BINID'].data > -1,:].reshape(-1,4563)),
#                        interpolation='nearest', aspect='auto')
#        pyplot.show()

        return mask


    def model_and_data_cube(self, prihdr, binned_spectra, binned_spectra_3d_hdu, stellar_continuum,
                            stellar_continuum_3d_hdu, emission_line_model,
                            emission_line_model_3d_hdu):
        """
        Constructs the 'FLUX', 'IVAR', 'MASK', 'WAVE', 'REDCORR', and
        'MODEL' model cube extensions.
        """
        #---------------------------------------------------------------
        ext = [ 'FLUX', 'IVAR', 'MASK', 'WAVE', 'REDCORR', 'MODEL' ]

        if binned_spectra is None or stellar_continuum is None or emission_line_model is None:
            # Construct and return the empty hdus
            return prihdr, DAPFitsUtil.empty_hdus(ext)

        #---------------------------------------------------------------
        # Add data to the primary header
        # TODO: Just copy header from maps file?
        prihdr = binned_spectra.rdxqa._initialize_primary_header(hdr=prihdr)
        prihdr = binned_spectra._initialize_primary_header(hdr=prihdr)
        prihdr = binned_spectra._add_method_header(prihdr)
        prihdr = binned_spectra._add_reddening_header(prihdr)
        prihdr = stellar_continuum._initialize_primary_header(hdr=prihdr)
        prihdr = stellar_continuum._add_method_header(prihdr)
        prihdr = emission_line_model._initialize_primary_header(hdr=prihdr)
        prihdr = emission_line_model._add_method_header(prihdr)

        #---------------------------------------------------------------
        # Get the extension headers
        hdr = [ DAPFitsUtil.finalize_dap_header(self.base_cubehdr, 'FLUX',
                                                bunit='1E-17 erg/s/cm^2/ang/spaxel', err=True,
                                                qual=True),
                DAPFitsUtil.finalize_dap_header(self.base_cubehdr, 'FLUX', hduclas2='ERROR',
                                                bunit='(1E-17 erg/s/cm^2/ang/spaxel)^{-2}',
                                                qual=True, prepend=False),
                DAPFitsUtil.finalize_dap_header(self.base_cubehdr, 'FLUX', hduclas2='QUALITY',
                                                 err=True, bit_type=self.bitmask.minimum_dtype(),
                                                 prepend=False),
                None,
                None,
                DAPFitsUtil.finalize_dap_header(self.base_cubehdr, 'MODEL',
                                                 bunit='1E-17 erg/s/cm^2/ang/spaxel', qual=True,
                                                 prepend=False)
              ]

        #---------------------------------------------------------------
        # Get the data arrays
        # Reddened flux data
        flux, ivar = binned_spectra.galext.apply(binned_spectra_3d_hdu['FLUX'].data, deredden=False,
                                                 ivar=binned_spectra_3d_hdu['IVAR'].data)
        flux = flux.astype(self.float_dtype)
        ivar = ivar.astype(self.float_dtype)
        # Best-fitting composite model
        model = stellar_continuum_3d_hdu['FLUX'].data + emission_line_model_3d_hdu['FLUX'].data \
                    + emission_line_model_3d_hdu['BASE'].data
        # with reddening
        model = binned_spectra.galext.apply(model, deredden=False).astype(self.float_dtype)

        # Spectral mask
        mask = self._get_model_and_data_mask(binned_spectra, binned_spectra_3d_hdu,
                                             stellar_continuum, stellar_continuum_3d_hdu,
                                             emission_line_model, emission_line_model_3d_hdu)

        data = [ flux, ivar, mask, binned_spectra['WAVE'].data.copy().astype(self.float_dtype),
                 binned_spectra['REDCORR'].data.copy().astype(self.float_dtype), model ]

        # Return the primary header and the list of HDUs
        return prihdr, DAPFitsUtil.list_of_image_hdus(data, hdr, ext)


    def _get_emission_line_model_mask(self, emission_line_model, emission_line_model_3d_hdu):
        """
        From the emission-line model object:
            - copy ARTIFACT
            - consolidate DIDNOTUSE, LOW_SNR, OUTSIDE_RANGE into IGNORED_EL
            - copy FIT_FAILED, NEAR_BOUND into ELFAILED
        """
        # Copy the base-level mask
        mask = self.bin_mask.copy()

        # Include the emission-line masks
        indx = emission_line_model.bitmask.flagged(emission_line_model_3d_hdu['MASK'].data,
                                                   flag='ARTIFACT')
        mask[indx] = self.bitmask.turn_on(mask[indx], 'ARTIFACT')

        indx = emission_line_model.bitmask.flagged(emission_line_model_3d_hdu['MASK'].data,
                                                   flag=['DIDNOTUSE', 'LOW_SNR', 'OUTSIDE_RANGE',
                                                         'TPL_PIXELS', 'TRUNCATED', 'PPXF_REJECT',
                                                         'INVALID_ERROR'])
        mask[indx] = self.bitmask.turn_on(mask[indx], 'ELIGNORED')

        indx = emission_line_model.bitmask.flagged(emission_line_model_3d_hdu['MASK'].data,
                                                   flag=['FIT_FAILED', 'NEAR_BOUND'])
        mask[indx] = self.bitmask.turn_on(mask[indx], 'ELFAILED')

        return mask


    def emission_line_model_cube(self, prihdr, binned_spectra, emission_line_model,
                                 emission_line_model_3d_hdu):
        """
        Constructs the 'EMLINE', 'EMLINE_BASE', and 'EMLINE_MASK' model
        cube extensions.
        """
        #---------------------------------------------------------------
        ext = [ 'EMLINE', 'EMLINE_BASE', 'EMLINE_MASK' ]

        if emission_line_model is None:
            return prihdr, DAPFitsUtil.empty_hdus(ext)

        # Add data to the primary header (Repeated from
        # model_and_data_cube
#        prihdr = emission_line_model._initialize_primary_header(hdr=prihdr)
#        prihdr = emission_line_model._add_method_header(prihdr)

        #---------------------------------------------------------------
        # Get the extension headers
        hdr = [ DAPFitsUtil.finalize_dap_header(self.base_cubehdr, 'EMLINE',
                                                 bunit='1E-17 erg/s/cm^2/ang/spaxel', qual=True),
                DAPFitsUtil.finalize_dap_header(self.base_cubehdr, 'EMLINE',
                                                 bunit='1E-17 erg/s/cm^2/ang/spaxel', qual=True),
                DAPFitsUtil.finalize_dap_header(self.base_cubehdr, 'EMLINE', hduclas2='QUALITY',
                                                 bit_type=self.bitmask.minimum_dtype()),
              ]

        #---------------------------------------------------------------
        # Get the data arrays
        # Get the reddened line model and baselines
        model = binned_spectra.galext.apply(emission_line_model_3d_hdu['FLUX'].data, deredden=False)
        base = binned_spectra.galext.apply(emission_line_model_3d_hdu['BASE'].data, deredden=False)

        # Get the mask data
        mask = self._get_emission_line_model_mask(emission_line_model, emission_line_model_3d_hdu)

        data = [ model.astype(self.float_dtype), base.astype(self.float_dtype), mask ]

        return prihdr, DAPFitsUtil.list_of_image_hdus(data, hdr, ext)


#-----------------------------------------------------------------------
def combine_binid_extensions(drpf, binned_spectra, stellar_continuum, emission_line_moments,
                             emission_line_model, spectral_indices, dtype=None):

    """
    Combine the bin IDs from the different analysis steps into a single
    multichannel extension.
    """
    if drpf is None:
        raise ValueError('DRP file must be provided.')

    hdr = DAPFitsUtil.finalize_dap_header( DAPFitsUtil.build_map_header(drpf,
                                'K Westfall & B Andrews <westfall@ucolick.org, andrewsb@pitt.edu>',
                                                    multichannel=True, maskname='MANGA_DAPPIXMASK'),
                                          'BINID', multichannel=True,
                                          channel_names=[ 'Binned spectra', 'Stellar continua',
                                                          'Em. line moments', 'Em. line models',
                                                          'Spectral indices' ])

    # Build the ID data for each analysis product
    binid = numpy.zeros(drpf.spatial_shape+(5,), dtype=int)
    if binned_spectra is not None:
        binid[:,:,0] = binned_spectra['BINID'].data
    if stellar_continuum is not None:
        binid[:,:,1] = stellar_continuum['BINID'].data
    if emission_line_moments is not None:
        binid[:,:,2] = emission_line_moments['BINID'].data
    if emission_line_model is not None:
        binid[:,:,3] = emission_line_model['BINID'].data
    if spectral_indices is not None:
        binid[:,:,4] = spectral_indices['BINID'].data

    _dtype = 'int' if dtype is None else dtype

    return DAPFitsUtil.list_of_image_hdus([ binid.astype(_dtype) ], [ hdr ], [ 'BINID' ])


def confirm_dap_types(drpf, obs, rdxqa, binned_spectra, stellar_continuum, emission_line_moments,
                      emission_line_model, spectral_indices):

    if not isinstance(drpf, DRPFits):
        raise TypeError('Input must have type DRPFits.')
    if obs is not None and not isinstance(obs, ObsInputPar):
        raise TypeError('Input must have type ObsInputPar.')
    if rdxqa is not None and not isinstance(rdxqa, ReductionAssessment):
        raise TypeError('Input must have type ReductionAssessment.')
    if binned_spectra is not None and not isinstance(binned_spectra, SpatiallyBinnedSpectra):
        raise TypeError('Input must have type SpatiallyBinnedSpectra.')
    if stellar_continuum is not None and not isinstance(stellar_continuum,
                                                        StellarContinuumModel):
        raise TypeError('Input must have type StellarContinuumModel.')
    if emission_line_moments is not None and not isinstance(emission_line_moments,
                                                            EmissionLineMoments):
        raise TypeError('Input must have type EmissionLineMoments.')
    if emission_line_model is not None and not isinstance(emission_line_model,
                                                          EmissionLineModel):
        raise TypeError('Input must have type EmissionLineModel.')
    if spectral_indices is not None and not isinstance(spectral_indices, SpectralIndices):
        raise TypeError('Input must have type SpectralIndices.')



def finalize_dap_primary_header(prihdr, drpf, obs, binned_spectra, stellar_continuum, dapsrc=None,
                                loggers=None, quiet=False):

    # Initialize the DAP quality flag
    dapqualbm = DAPQualityBitMask(dapsrc=dapsrc)
    drp3qualbm = DRPQuality3DBitMask()
    dapqual = dapqualbm.minimum_dtype()(0)          # Type casting original flag to 0
    if drp3qualbm.flagged(drpf['PRIMARY'].header['DRP3QUAL'], flag='CRITICAL'):
        if not quiet:
            log_output(loggers, 1, logging.INFO, 'DRP File is flagged CRITICAL!')
        dapqual = dapqualbm.turn_on(dapqual, 'CRITICAL')
        dapqual = dapqualbm.turn_on(dapqual, 'DRPCRIT')

    # Flag the file as CRITICAL if the stellar continuum fits are bad
    # for all spectra
    if stellar_continuum is not None:
        mask = stellar_continuum.bitmask.flagged(stellar_continuum['PAR'].data['MASK'],
                                                 ['NO_FIT', 'FIT_FAILED', 'INSUFFICIENT_DATA',
                                                  'NEAR_BOUND'])
        if numpy.all(mask):
            dapqual = dapqualbm.turn_on(dapqual, 'CRITICAL')

    # Signify that the Voronoi binning resulted in a single bin
    if binned_spectra is not None and binned_spectra.method['binclass'] is not None \
            and binned_spectra.method['binclass'].bintype == 'voronoi' \
            and binned_spectra.nbins == 1:
        dapqual = dapqualbm.turn_on(dapqual, 'SINGLEBIN')

    # Input photometric geometry and scale were invalid
    if obs is not None and not (obs.valid_ell and obs.valid_pa and obs.valid_reff):
        dapqual = dapqualbm.turn_on(dapqual, 'BADGEOM')

    # Determine if there's a foreground star
    if numpy.sum(drpf.bitmask.flagged(drpf['MASK'].data, flag='FORESTAR')) > 0:
        dapqual = dapqualbm.turn_on(dapqual, 'FORESTAR')

    # Commit the quality flag to the header
    if dapqualbm.flagged(dapqual, 'CRITICAL'):
        warnings.warn('DAP files are flagged CRITICAL.')
    prihdr['DAPQUAL'] = (dapqual, 'DAP quality bitmask')

    # Finalize authors
    prihdr['AUTHOR'] = 'K Westfall, B Andrews <westfall@ucolick.org, andrewsb@pitt.edu>'

    return prihdr


def add_snr_metrics_to_header(hdr, drpf, r_re, dapsrc=None):
    """
    For all valid spaxels within 1 Re < R < 1.5 Re calculate the median
    S/N and combined S/N (including covariance) in the griz bands.

    Adds the following keywords (8) to the header:
        - SNR?MED: Median S/N within 1-1.5 Re in each of the griz bands.
        - SNR?RING: The combined S/N of all spaxels within 1-1.5 Re in
          the griz bands, incuding covariance

    The inclusion of covariance is based on the exact calculation of the
    correlation matrix at the flux-weighted center of each band, and
    assuming no change in the correlation between the response-weighted
    correlation within the broad band.

    Args:
        hdr (`astropy.io.fits.Header`_): Header object to edit
        drpf (:class:`mangadap.drpfits.DRPFits`): DRP fits cube
        r_re (numpy.ndarray): *Flattened* array with the semi-major axis
            radius in units of the effective radius for all spaxels in
            the datacube.  Shape should be (Nx*Ny,), where the shape of
            the datacube is (Nx,Ny,Nwave).

    Returns:
        `astropy.io.fits.Header`_: The edited header object.

    Raises:
        FileNotFoundError: Raised if any of the response function files
            cannot be found.
    """
    _dapsrc = dap_source_dir() if dapsrc is None else str(dapsrc)
    filter_response_file = [ os.path.join(_dapsrc, 'data', 'filter_response',
                                           f) for f in [ 'gunn_2001_g_response.db',
                                                         'gunn_2001_r_response.db',
                                                         'gunn_2001_i_response.db',
                                                         'gunn_2001_z_response.db'] ]
    for f in filter_response_file:
        if not os.path.isfile(f):
            raise FileNotFoundError('{0} does not exist!'.format(f))

    # Set the header keywords
    key_med = ['SNRGMED', 'SNRRMED', 'SNRIMED', 'SNRZMED' ]
    com_med = [ 'Median g-band SNR from 1-1.5 Re', 'Median r-band SNR from 1-1.5 Re',
                'Median i-band SNR from 1-1.5 Re', 'Median z-band SNR from 1-1.5 Re' ]
                    
    key_ring = ['SNRGRING', 'SNRRRING', 'SNRIRING', 'SNRZRING' ]
    com_ring = [ 'g-band SNR in 1-1.5 Re bin', 'r-band SNR in 1-1.5 Re bin',
                 'i-band SNR in 1-1.5 Re bin', 'z-band SNR in 1-1.5 Re bin' ]

    # Run the S/N calculation; this is the same calculation as done in
    # mangadap.proc.spatialbinning.VoronoiBinning.sn_calculation_covariance_matrix
    nfilter = len(filter_response_file)
    for i in range(nfilter):
        response_func = numpy.genfromtxt(filter_response_file[i])[:,:2]
        signal, variance, snr, covar = drpf.flux_stats(response_func=response_func,
                                                       flag=['DONOTUSE', 'FORESTAR'],
                                                       covar=True) #, correlation=True)
        # Get the spaxels within the radius limits
        indx = numpy.arange(r_re.size)[(r_re > 1) & (r_re < 1.5) & numpy.invert(signal.mask)]
        if len(indx) == 0:
            hdr[key_med[i]] = (0., com_med[i])
            hdr[key_ring[i]] = (0., com_ring[i])
            continue

        if covar is None:
            warnings.warn('Covariance not available!  Continuing without it.')
            covar = Covariance.from_variance(variance, correlation=True)
            
        # Get the appropriate covariance pixels to select
        ci, cj = map( lambda x: x.ravel(), numpy.meshgrid(indx, indx) )

        # Use the covariance matrix from the single wavelength channel
        # calculation, but renormalize it to the mean variance over the
        # response function
        covar = covar.apply_new_variance(variance).toarray()
        # Set the median S/N ...
        hdr[key_med[i]] = (numpy.ma.median(snr[indx]), com_med[i])
        # ... and the combined S/N
        hdr[key_ring[i]] = (numpy.ma.sum(signal[indx])/numpy.sqrt(numpy.sum(covar[ci,cj])),
                            com_ring[i])

    return hdr






