# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Defines a class used to interface with the final maps files produced by
the MaNGA Data Analysis Pipeline (DAP).

.. todo::
    - Allow DAPFits to read/access both the MAPS and the LOGCUBE files,
      not just the former.

----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import logging
import time
import os
import warnings

from IPython import embed

import numpy

from astropy.io import fits
import astropy.constants
import astropy.units
from astropy.cosmology import FlatLambdaCDM

from .datacube import DataCube
#from .datacube import MaNGADataCube, MUSEDataCube
from .util.constants import DAPConstants
from .util.drpfits import DRPQuality3DBitMask
from .util.fitsutil import DAPFitsUtil
from .util.dapbitmask import DAPBitMask
from .util.log import log_output
from .util.fileio import channel_dictionary
from .util.geometry import SemiMajorAxisCoo
from .util.covariance import Covariance
from .config import defaults
from .proc.reductionassessments import ReductionAssessment
from .proc.spatiallybinnedspectra import SpatiallyBinnedSpectra
from .proc.stellarcontinuummodel import StellarContinuumModel
from .proc.emissionlinemoments import EmissionLineMoments
from .proc.emissionlinemodel import EmissionLineModel
from .proc.spectralindices import SpectralIndices

from matplotlib import pyplot

#-----------------------------------------------------------------------
class DAPQualityBitMask(DAPBitMask):
    """
    Class used to flag/interpret global quality bit.

    The defined bits are listed at :ref:`metadatamodel-dapqual`.

    .. todo::
        - Force read IDLUTILS version as opposed to internal one?
    """
    cfg_root = 'dap_quality_bits'


class DAPMapsBitMask(DAPBitMask):
    """
    Class used to flag/interpret bits in DAP MAPS files.

    The defined bits are listed at :ref:`metadatamodel-dappixmask`.

    .. todo::
        - Force read IDLUTILS version as opposed to internal one?
    """
    cfg_root = 'dap_maps_bits'


class DAPCubeBitMask(DAPBitMask):
    """
    Class used to flag/interpret bits in DAP model LOGCUBE files.

    The defined bits are listed at :ref:`metadatamodel-dapspecmask`.

    .. todo::
        - Force read IDLUTILS version as opposed to internal one?
    """
    cfg_root = 'dap_cube_bits'


##-----------------------------------------------------------------------
## TODO: This class is out of date!
#class DAPFits:
#    """
#    Object used to hold properties of and read data from a DAP-produced file.
#
#    Args:
#        plate (:obj:`int`):
#            Plate number
#        ifudesign (:obj:`int`):
#            IFU design
#        method (:obj:`str`):
#            Method type (e.g., ``'SPX-GAU-MILESTHIN'``)
#        drpver (:obj:`str`, optional):
#            DRP version, which is used to define the default DAP
#            analysis path. Default is defined by
#            :func:`mangadap.config.defaults.drp_version`.
#        dapver (:obj:`str`, optional):
#            DAP version, which is used to define the default DAP
#            analysis path. Default is defined by
#            :func:`mangadap.config.defaults.dap_version`.
#        analysis_path (:obj:`str`, optional):
#            The path to the top level directory containing the DAP
#            output files for a given DRP and DAP version. Default is
#            defined by
#            :func:`mangadap.config.defaults.dap_analysis_path`.
#        directory_path (:obj:`str`, optional):
#            The exact path to the DAP file. Default is defined by
#            :func:`mangadap.config.defaults.dap_method_path`.
#        reference_path (:obj:`str`, optional):
#            The exact path of the DAP reference directory. Default is
#            defined by
#            :func:`mangadap.config.defaults.dap_method_path`.
#        par_file (:obj:`str`, optional):
#            SDSS parameter file used to provide input parameters for
#            the DAP. Default is defined by
#            :func:`mangadap.config.defaults.dap_par_file`.
#        read (:obj:`bool`, optional):
#            Read the DAP file upon instantiation of the object.
#        checksum (:obj:`bool`, optional):
#            Flag for `astropy.io.fits.open`_ checksum operation.
#
#    Attributes:
#        plate (:obj:`int`):
#            Plate number
#        ifudesign (:obj:`int`):
#            IFU design
#        method (:obj:`str`):
#            Method type
#        drpver (:obj:`str`):
#            DRP version
#        dapver (:obj:`str`):
#            DAP version
#        analysis_path (:obj:`str`):
#            The path to the top level directory containing the DAP
#            output files for a given DRP and DAP version.
#        directory_path (:obj:`str`):
#            The exact path to the DAP file.
#        reference_path (:obj:`str`):
#            The exact path of the DAP reference directory.
#        par_file (:obj:`str`):
#            SDSS parameter file used to provide input parameters for
#            the DAP.
#        dapqual (:class:`mangadap.dapmaps.DAPQualityBitMask`):
#            Global quality bit
#        bitmask (:class:`mangadap.dapmaps.DAPMapsBitMask`):
#            Map mask bits
#        hdu (`astropy.io.fits.HDUList`_):
#            HDUList read from the DRP file
#        par (:class:`mangadap.par.ObsInputPar`):
#            List of parameters used by the DAP.
#        checksum (:obj:`bool`):
#            Boolean to use for checksum in `astropy.io.fits.open`_.
#    """
#    def __init__(self, plate, ifudesign, method, drpver=None, dapver=None, analysis_path=None,
#                 directory_path=None, reference_path=None, par_file=None, read=True,
#                 checksum=False):
#
#        # Set the attributes, forcing a known type
#        self.plate = int(plate)
#        self.ifudesign = int(ifudesign)
#        self.method = method
#
#        # TODO: Do we need drpver, dapver, analysis_path to be kept
#        # by self?
#        if directory_path is None:
#            self.drpver = defaults.drp_version() if drpver is None else str(drpver)
#            self.dapver = defaults.dap_version() if dapver is None else str(dapver)
#            self.analysis_path = defaults.dap_analysis_path(self.drpver, self.dapver) \
#                                 if analysis_path is None else str(analysis_path)
#            self.directory_path = defaults.dap_method_path(method, plate=self.plate,
#                                                           ifudesign=self.ifudesign,
#                                                           drpver=self.drpver, dapver=self.dapver,
#                                                           analysis_path=self.analysis_path)
#        else:
#            self.drpver = None
#            self.dapver = None
#            self.analysis_path = None
#            self.directory_path = str(directory_path)
#
#        # Set the reference directory path
#        self.reference_path = defaults.dap_method_path(method, plate=self.plate,
#                                                       ifudesign=self.ifudesign, ref=True,
#                                                       drpver=self.drpver, dapver=self.dapver,
#                                                       analysis_path=self.analysis_path) \
#                                    if reference_path is None else reference_path
#
#        # drpver, dapver, and analysis_path can be None and par_file
#        # will still only use the above defined directory_path
#        self.par_file = defaults.dap_par_file(self.plate, self.ifudesign, 'CUBE',
#                                              drpver=self.drpver, dapver=self.dapver,
#                                              analysis_path=self.analysis_path,
#                                              directory_path=self.reference_path) \
#                                    if par_file is None else str(par_file)
#
#        # Set the bitmasks
#        self.bitmask = DAPMapsBitMask()
#
#        # Read in the data
#        self.hdu = None
#        self.spatial_shape = None
#        self.smap_ext = None      # Extensions with single maps
#        self.mmap_ext = None      # Extensions with multiple maps
#        self.par = None
#        self.channel_dict = None
#        self.checksum = checksum
#        if read:
#            self.open_hdu(checksum=self.checksum)
#            self.read_par()
#
#    def __getitem__(self, key):
#        """Access elements of the hdu."""
#        if self.hdu is None:
#            return None
#        return self.hdu[key]
#
#    def open_hdu(self, permissions='readonly', checksum=False, quiet=True):
#        """
#        Open the fits file and save it to :attr:`hdu`; if :attr:`hdu`
#        is not None, the function returns without re-reading the
#        data.
#
#        Args:
#            permissions (:obj:`str`, optional):
#                Open the fits file with these read/write permissions.
#            checksum (:obj:`bool`, optional):
#                Flag for `astropy.io.fits.open`_ checksum operation.
#                This overrides the internal :attr:`checksum`
#                attribute **for the current operation only**.
#            quiet (:obj:`bool`, optional):
#                Suppress terminal output
#
#        Raises:
#            FileNotFoundError: Raised if the DAP file doesn't exist.
#
#        """
#        if self.hdu is not None:
#            if not quiet:
#                print('DAP file already open!')
#            return
#
#        inp = self.file_path()
#        if not os.path.exists(inp):
#            raise FileNotFoundError('Cannot open file: {0}'.format(inp))
#
#        # Open the fits file with the requested read/write permission
#        self.hdu = DAPFitsUtil.read(inp, permissions=permissions, checksum=checksum)
#
#        self.spatial_shape = self.hdu['BINID'].data.shape
#
#        self.smap_ext = []      # Extensions with single maps
#        self.mmap_ext = []      # Extensions with multiple maps
#        for h in self.hdu:
#            if h.data is not None:
#                if len(h.shape) == 2:
#                    self.smap_ext += [ h.name ]
#                if len(h.shape) == 3:
#                    self.mmap_ext += [ h.name ]
#
#        if len(self.mmap_ext) > 0:
#            self.channel_dict = { self.mmap_ext[0]: channel_dictionary(self.hdu, self.mmap_ext[0]) }
#            for i in range(1,len(self.mmap_ext)):
#                self.channel_dict[self.mmap_ext[i]] = channel_dictionary(self.hdu, self.mmap_ext[i])
#
##    def read_par(self, quiet=True):
##        """
##        Open the parameter file and save it to :attr:`par`; if
##        :attr:`par` is not None, the function returns without re-reading
##        the data.
##
##        Args:
##            quiet (:obj:`bool`, optional):
##                Suppress terminal output
##        """
##        if self.par is not None:
##            if not quiet:
##                print('Parameter file already read!')
##            return
##
##        if not os.path.exists(self.par_file):
##            warnings.warn('Input parameters unavailable; cannot open {0}'.format(self.par_file))
##            return
##            
##        # Set the input parameters
##        # TODO: Is all of this also in the file header?
##        self.par = ObsInputPar.from_par_file(self.par_file)
#
#    def file_name(self):
#        """Return the name of the DAP file"""
#        return defaults.dap_file_name(self.plate, self.ifudesign, self.method, mode='MAPS')
#
#    def file_path(self):
#        """Return the full path to the DAP file"""
#        return os.path.join(self.directory_path, self.file_name())
#    
#    def list_channels(self, ext):
#        """
#        Provide a list of the channels in a given extension.
#        """
#        if self.hdu is None:
#            self.open_hdu(checksum=self.checksum)
#        if ext in self.smap_ext:
#            print('{0} contains a single channel.'.format(ext))
#            return
#        print('Channels in extension {0}'.format(ext))
#        print('#  {0:>3} {1}'.format('CH','KEY'))
#        srt = numpy.argsort(list(self.channel_dict[ext].values()))
#        for k, v in zip(numpy.array(list(self.channel_dict[ext].keys()))[srt],
#                        numpy.array(list(self.channel_dict[ext].values()))[srt]):
#            print('   {0:>3} {1}'.format(v,k))
#
#    def thin_disk_polar_coo(self, xc=None, yc=None, rot=None, pa=None, inc=None, flip=False):
#        r"""
#        Calculate the disk-plane coordinates for all the binned spectra
#        using :class:`mangadap.util.geometry.SemiMajorAxisCoo`.  See the
#        documentation their for further information regarding the
#        calculations.
#
#        If not provided, the default position angle and inclination will
#        be taken from the input parameter file, if available.  If the
#        parameter file cannot be read, the default values for the
#        position angle and ellipticity in
#        :class:`mangadap.util.geometry.SemiMajorAxisCoo` are used.
#
#        Args:
#            xc,yc (float,float): (**Optional**) a reference on-sky
#                position relative to the center of the projected plane
#                (galaxy center).
#            rot (float): (**Optional**) Cartesian rotation of the
#                focal-plane relative to the sky-plane (+x toward East;
#                +y toward North).
#            pa (float): (**Optional**) On-sky position angle of the
#                major axis of the ellipse traced on the sky by a circle
#                in the projected plane, defined as the angle from North
#                through East.
#            inc (float): (**Optional**) Inclination of the plane, angle
#                between the disk normal and the normal of the sky plane
#                (inc=0 is a face-on plane).
#            flip (bool): (**Optional**) Offset the provided position
#                angle by 180 deg; e.g., to set the position angle will
#                be defined along the receding major axis.
#
#        Returns:
#            numpy.ndarray: Two numpy arrays with the projected polar
#                coordinates: :math:`R, \theta`.
#        """
#
#        # Flip if requested
#        if pa is not None and flip:
#            pa += 180
#            if pa >= 360:
#                pa -= 360
#
#        # Position angle and inclination
#        if pa is None or inc is None:
#            try:
#                self.read_par()                # Try to read the parameter file
#            except:
#                warnings.warn('Using default pa and/or inclination.')
#
#            if self.par is not None:
#                if pa is None:
#                    pa = self.guess_position_angle(flip)
#                if inc is None:
#                    inc = self.guess_inclination()
#
#        # Declare the galaxy coo object
#        disk_plane = SemiMajorAxisCoo(xc=xc, yc=yc, rot=rot, pa=pa,
#                                      ell=numpy.cos(numpy.radians(inc)))
#
#        # Calculate the in-plane radius and azimuth for each bin
#        if self.hdu is None:
#            self.open_hdu(checksum=self.checksum)
#        return disk_plane.polar(self.hdu['SPX_SKYCOO'].data[:,:,0],
#                                self.hdu['SPX_SKYCOO'].data[:,:,1])
#
#    
#    def guess_inclination(self, q0=None):
#        r"""
#
#        Return a guess inclination based on the isophotal ellipticity
#        from the input parameter file using the equation:
#
#        .. math::
#        
#            \cos^2 i = \frac{ q^2 - q_0^2 }{ 1 - q_0^2 },
#
#        where :math:`q` is the observed oblateness (the semi-minor to
#        semi-major axis ratio, :math:`q = 1-\epsilon = b/a`) and
#        :math:`q_0` is the intrinsic oblateness.
#        
#        If :math:`q < q_0`, the function returns a 90 degree
#        inclination.
#        
#        Args:
#            q0 (float): (**Optional**) The intrinsic oblateness of the
#                system.  If not provided, the system is assumed to be
#                infinitely thin.
#
#        Raises:
#            Exception: Raised if the input intrinsic oblatenss is not
#                less than one and greater than 0.
#        """
#        self.read_par()
#        if q0 is None or q0 == 0.0:
#            return numpy.degrees( numpy.arccos(1.0 - self.par['ell']) )
#
#        if not (q0 < 1.0) or (q0 < 0.0):
#            raise Exception('Intrinsic oblateness must be 0.0 <= q0 < 1.0!')
#
#        q = 1.0 - self.par['ell']
#        if q < q0:
#            return 90.0
#
#        q02 = q0*q0
#        return numpy.degrees(numpy.arccos(numpy.sqrt( (q*q - q02) / (1.0 - q02) )))
#
#
#    def guess_position_angle(self, flip=False):
#        """
#        Return the guess position angle from the parameter file.
#
#        Args:
#            flip (bool): (**Optional**) Offset the provided position
#                angle by 180 deg.
#        """
#        self.read_par()
#        pa = self.par['pa'] if not flip else self.par['pa'] + 180
#        return pa if pa < 360 else pa-360
#       
#
#    def effective_radius(self):
#        """
#        Return the effective radius from the parameter file.
#        """
#        self.read_par()
#        return self.par['reff']
#
#
##    def nsa_ellipticity(self):
##        """
##        Return the ellipticity from the parameter file.
##        """
##        self.read_par()
##        return self.par['ell']
##
##
#    def guess_cz(self):
#        """Return the guess redshift (cz) from the parameter file."""
#        self.read_par()
#        return self.par['vel']
#
#
#    def unique_indices(self):
#        """
#        Use the BINID extension to select the spaxels with unique data.
#        
#        Returns:
#            numpy.ndarray: Returns the list of indices **in the
#            flattened array of a map** that are unique.
#
#        """
#        if self.hdu is None:
#            self.open_hdu(checksum=self.checksum)
#
#        unique_bins, indx = numpy.unique(self.hdu['BINID'].data.ravel(), return_index=True)
#        return indx[unique_bins > -1]
#
#
#    def unique_mask(self):
#        """
#        Construct a boolean mask for the unique bins.
#        
#        The returned map is True for spaxels that are part of a unique
#        bin, Fals otherwise.
#        """
#        indx = self.unique_indices()
#        msk = numpy.zeros(self.spatial_shape, dtype=numpy.bool)
#        msk.ravel()[indx] = True
#        return msk
#
#
#    def channel_map(self, ext, channel=None, flag=None):
#        """
#        Return a 2D array with the map for the specified channel.
#        """
#        if channel is None:
#            # Channel must be satisfied
#            if len(self.hdu[ext].data.shape) > 2:
#                raise ValueError('Must specify channel for extension {0}'.format(ext))
#            # Return unmasked array
#            if flag is None:
#                return numpy.ma.MaskedArray(self.hdu[ext].data.copy())
#            # Attempt to get the mask extension from the header
#            try:
#                mask_ext = self.hdu[ext].header['QUALDATA']
#            except:
#                raise ValueError('Mask extension not specified in header.  Unable to mask array.')
#            # Return a masked array, masking the specified flags
#            return numpy.ma.MaskedArray(self.hdu[ext].data.copy(),
#                                        mask=self.bitmask.flagged(self.hdu[mask_ext].data,
#                                                                  flag=flag))
#
#        if channel not in list(self.channel_dict[ext].keys()):
#            raise ValueError('No channel {0} in extension {1}'.format(ext, channel))
#
#        # Return unmasked array
#        if flag is None:
#            return numpy.ma.MaskedArray(self.hdu[ext].data[:,:,self.channel_dict[ext][channel]])
#        # Attempt to get the mask extension from the header
#        try:
#            mask_ext = self.hdu[ext].header['QUALDATA']
#        except:
#            raise ValueError('Mask extension not specified in header.  Unable to mask array.')
#        
#        return numpy.ma.MaskedArray(self.hdu[ext].data[:,:,self.channel_dict[ext][channel]],
#                                    mask=self.bitmask.flagged(
#                                      self.hdu[mask_ext].data[:,:,self.channel_dict[ext][channel]],
#                                      flag=flag))


#-----------------------------------------------------------------------
class construct_maps_file:
    """
    Construct a DAP MAPS file.

    Set as a class to enforce a namespace. Should not be used as an
    object!

    Should force all intermediate objects to be provided.

    """
    def __init__(self, cube, metadata, rdxqa=None, binned_spectra=None, stellar_continuum=None,
                 emission_line_moments=None, emission_line_model=None, spectral_indices=None,
                 redshift=None, dapver=None, analysis_path=None, directory_path=None,
                 output_file=None, clobber=True, loggers=None, quiet=False,
                 single_precision=False):

        #---------------------------------------------------------------
        # Initialize the reporting
        self.loggers = None if loggers is None else loggers
        self.quiet = quiet

        self.float_dtype = 'float32' if single_precision else 'float'

        #---------------------------------------------------------------
        # Check input types
        confirm_dap_types(cube, rdxqa, binned_spectra, stellar_continuum, emission_line_moments,
                          emission_line_model, spectral_indices)

        #---------------------------------------------------------------
        # Set the output paths
        self.cube = cube
        self.meta = metadata
        self.method = None
        self.directory_path = None
        self.output_file = None
        self._set_paths(directory_path, dapver, analysis_path, output_file, binned_spectra,
                        stellar_continuum, emission_line_model)

        # Save input for reference
        self.spatial_shape = self.cube.spatial_shape
        self.redshift = redshift

        # Report
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, '-'*50)
            log_output(loggers, 1, logging.INFO, '{0:^50}'.format('CONSTRUCTING OUTPUT MAPS'))
            log_output(self.loggers, 1, logging.INFO, '-'*50)
            log_output(self.loggers, 1, logging.INFO,
                       'Output path: {0}'.format(self.directory_path))
            log_output(self.loggers, 1, logging.INFO,
                       'Output file: {0}'.format(self.output_file))
            log_output(self.loggers, 1, logging.INFO,
                       'Output maps have shape {0}'.format(self.spatial_shape))
            log_output(self.loggers, 1, logging.INFO, 'Redshift: {0}'.format(self.redshift))

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
        prihdr = DAPFitsUtil.initialize_dap_primary_header(self.cube, maskname='MANGA_DAPPIXMASK')

        # Add the DAP method
        prihdr['DAPTYPE'] = (defaults.dap_method(binned_spectra.method['key'],
                                    stellar_continuum.method['fitpar']['template_library_key'],
                                    'None' if emission_line_model is None
                                        else emission_line_model.method['continuum_tpl_key']),
                             'DAP analysis method')
        # Add the format of this file
        prihdr['DAPFRMT'] = ('MAPS', 'DAP data file format')

        # Get the base map headers
        self.multichannel_maphdr \
                = DAPFitsUtil.build_map_header(self.cube.fluxhdr,
                                'K Westfall <westfall@ucolick.org> & SDSS-IV Data Group',
                                               multichannel=True, maskname='MANGA_DAPPIXMASK')

        self.singlechannel_maphdr \
                = DAPFitsUtil.build_map_header(self.cube.fluxhdr,
                                'K Westfall <westfall@ucolick.org> & SDSS-IV Data Group',
                                               maskname='MANGA_DAPPIXMASK')

        #---------------------------------------------------------------
        # Initialize the pixel mask
        self.bitmask = DAPMapsBitMask()

        # Get the mask for the binned spectra
        self.bin_mask = self._build_binning_mask(binned_spectra)

        #---------------------------------------------------------------
        # Construct the hdu list for each input object.
        # Reduction assessments:
        rdxqalist = self.reduction_assessment_maps(prihdr, rdxqa)
        # Construct the BINID extension
        binidlist = combine_binid_extensions(self.cube, binned_spectra, stellar_continuum,
                                             emission_line_moments, emission_line_model,
                                             spectral_indices, dtype='int32')
        # Binned spectra:
        bspeclist = self.binned_spectra_maps(prihdr, binned_spectra)
        # Stellar-continuum fits:
        contlist = self.stellar_continuum_maps(prihdr, stellar_continuum)
        # Emission-line moments:
        emlmomlist = self.emission_line_moment_maps(prihdr, emission_line_moments)
        # Emission-line models:
        elmodlist = self.emission_line_model_maps(prihdr, emission_line_model)

        # Spectral indices:
        sindxlist = self.spectral_index_maps(prihdr, spectral_indices)

        # Save the data to the hdu attribute
        prihdr = add_snr_metrics_to_header(prihdr, self.cube, rdxqalist[1].data[:,:,1].ravel())
        
        prihdr = finalize_dap_primary_header(prihdr, self.cube, self.meta, binned_spectra,
                                             stellar_continuum, loggers=self.loggers,
                                             quiet=self.quiet)
        self.hdu = fits.HDUList([ fits.PrimaryHDU(header=prihdr),
                                  *rdxqalist,
                                  *binidlist,
                                  *bspeclist,
                                  *contlist,
                                  *emlmomlist,
                                  *elmodlist,
                                  *sindxlist
                                ])

        extensions = [ h.name for h in self.hdu ]


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
        ext = ['BIN_MFLUX', 'STELLAR_VEL', 'STELLAR_SIGMA', 'EMLINE_SFLUX', 'EMLINE_SEW',
               'EMLINE_GFLUX', 'EMLINE_GEW', 'EMLINE_GVEL', 'EMLINE_GSIGMA', 'SPECINDEX']
        for e in ext:
            if '{0}_MASK'.format(e) not in extensions \
                    or '{0}_IVAR'.format(e) not in extensions \
                    or self.hdu['{0}_MASK'.format(e)].data is None \
                    or self.hdu['{0}_IVAR'.format(e)].data is None:
                continue
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
                   stellar_continuum, emission_line_model):
        """
        Set the paths relevant to the map file.
        """
        # The output method directory is, for now, the combination of
        # the binned_spectrum and stellar_continuum method keys
        if directory_path is None and (binned_spectra is None or stellar_continuum is None):
            raise ValueError('Could not define output directory path.')

        # Set the output directory path
        # TODO: Get DAP version from __version__ string

        self.method = defaults.dap_method(binned_spectra.method['key'],
                                    stellar_continuum.method['fitpar']['template_library_key'],
                                    'None' if emission_line_model is None 
                                        else emission_line_model.method['continuum_tpl_key']) \
                                if directory_path is None else None
        self.directory_path = defaults.dap_method_path(self.method, plate=self.cube.plate,
                                                       ifudesign=self.cube.ifudesign,
                                                       drpver=self.cube.drpver, dapver=dapver,
                                                       analysis_path=analysis_path) \
                                                if directory_path is None else str(directory_path)

        # Set the output file
        self.output_file = defaults.dap_file_name(self.cube.plate, self.cube.ifudesign,
                                                  self.method, mode='MAPS') \
                                    if output_file is None else str(output_file)

    def _consolidate_donotuse(self, mask):
        return self.bitmask.consolidate(mask, ['NOCOV', 'LOWCOV', 'DEADFIBER', 'FORESTAR',
                                               'NOVALUE', 'MATHERROR', 'FITFAILED', 'NEARBOUND'],
                                        'DONOTUSE')

    def _build_binning_mask(self, binned_spectra):
        """
        Consolidate the map mask values into NOVALUE.
        """
        # Consolidate everything in the map mask into NOVALUE
        indx = binned_spectra.bitmask.flagged(binned_spectra['MAPMASK'].data)
        mask = numpy.zeros(indx.shape, dtype=self.bitmask.minimum_dtype())
        mask[indx] = self.bitmask.turn_on(mask[indx], 'NOVALUE')

        # Isolate FORESTAR
        indx = binned_spectra.bitmask.flagged(binned_spectra['MAPMASK'].data, flag='FORESTAR')
        mask[indx] = self.bitmask.turn_on(mask[indx], 'FORESTAR')

        # Marginalize the binned spectra mask just for NONE_IN_STACK and
        # convert it to NOVALUE
        bin_mask = DAPFitsUtil.marginalize_mask(binned_spectra['MASK'].data,
                                                inp_flags='NONE_IN_STACK',
                                                inp_bitmask=binned_spectra.bitmask, 
                                                out_flag='NOVALUE', out_bitmask=self.bitmask,
                                                dispaxis=1)

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
        # Consolidate everything in the map mask into NOVALUE
        indx = stellar_continuum.bitmask.flagged(stellar_continuum['MAPMASK'].data)
        mask = numpy.zeros(indx.shape, dtype=self.bitmask.minimum_dtype())
        mask[indx] = self.bitmask.turn_on(mask[indx], 'NOVALUE')

        # Isolate FORESTAR
        indx = stellar_continuum.bitmask.flagged(stellar_continuum['MAPMASK'].data,
                                                 flag='FORESTAR')
        mask[indx] = self.bitmask.turn_on(mask[indx], 'FORESTAR')
        
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

        NO_ABSORPTION_CORRECTION propagated to NOCORRECTION

        DIVBYZERO propagated to MATHERROR

        Second moments not provided so no need to include UNDEFINED_MOM2
        """
        # Consolidate everything in the map mask into NOVALUE
        indx = emission_line_moments.bitmask.flagged(emission_line_moments['MAPMASK'].data)
        mask = numpy.zeros(indx.shape, dtype=self.bitmask.minimum_dtype())
        mask[indx] = self.bitmask.turn_on(mask[indx], 'NOVALUE')

        # Isolate FORESTAR
        indx = emission_line_moments.bitmask.flagged(emission_line_moments['MAPMASK'].data,
                                                     flag='FORESTAR')
        mask[indx] = self.bitmask.turn_on(mask[indx], 'FORESTAR')

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
        # Consolidate everything in the map mask into NOVALUE
        indx = emission_line_model.bitmask.flagged(emission_line_model['MAPMASK'].data)
        mask = numpy.zeros(indx.shape, dtype=self.bitmask.minimum_dtype())
        mask[indx] = self.bitmask.turn_on(mask[indx], 'NOVALUE')

        # Isolate FORESTAR
        indx = emission_line_model.bitmask.flagged(emission_line_model['MAPMASK'].data,
                                                   flag='FORESTAR')
        mask[indx] = self.bitmask.turn_on(mask[indx], 'FORESTAR')

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
        # Consolidate everything in the map mask into NOVALUE
        indx = spectral_indices.bitmask.flagged(spectral_indices['MAPMASK'].data)
        mask = numpy.zeros(indx.shape, dtype=self.bitmask.minimum_dtype())
        mask[indx] = self.bitmask.turn_on(mask[indx], 'NOVALUE')

        # Isolate FORESTAR
        indx = spectral_indices.bitmask.flagged(spectral_indices['MAPMASK'].data, flag='FORESTAR')
        mask[indx] = self.bitmask.turn_on(mask[indx], 'FORESTAR')

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

#        if abs_only:
#            # Get the number of absorption-line indices. Indices *must*
#            # be ordered as [0:nabs], [nabs:nabs+nbhd]
#            nabs, nbhd = SpectralIndices.count_indices(spectral_indices.absdb,
#                                                       spectral_indices.bhddb)
#            mask[:,:,nabs:] = self.bitmask.turn_on(mask[:,:,nabs:], 'NOVALUE')

        return self._consolidate_donotuse(mask)

    @staticmethod
    def _get_kpc_per_arcsec(z):
        if z <= 0:
            warnings.warn('Systemic velocity <=0; radii in h^-1 kpc not provided.')
            return -1.

        H0 = 100 * astropy.units.km / astropy.units.s / astropy.units.Mpc
        cosmo = FlatLambdaCDM(H0=H0, Om0=0.3)
        return numpy.radians(1/3600) * 1e3 * cosmo.angular_diameter_distance(z).value

    def reduction_assessment_maps(self, prihdr, rdxqa):
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
                                                               'R h/kpc', 'Elliptical azimuth'],
                                                channel_units=['arcsec', '', 'kpc/h',
                                                               'degrees']),
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
        spx_ellcoo = numpy.repeat(spx_ellcoo, (3,1), axis=2)

        # Calculate the radius normalized by the effective radius
        if 'reff' in self.meta.keys() and self.meta['reff'] is not None:
            spx_ellcoo[:,:,1] /= self.meta['reff']

        # Calculate the radius in units of h^-1 kpc
        hkpc_per_arcsec = self.__class__._get_kpc_per_arcsec(self.meta['z'])
        if hkpc_per_arcsec > 0:
            spx_ellcoo[:,:,2] *= hkpc_per_arcsec
        else:
            spx_ellcoo[:,:,2] = -1
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

    
    def binned_spectra_maps(self, prihdr, binned_spectra):
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
                                                                'R/Re', 'R h/kpc',
                                                                'Lum. weighted elliptical azimuth'],
                                                channel_units=['arcsec', '', 'kpc/h', 'degrees']),
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

        # Calculate the radius normalized by the effective radius
        if 'reff' in self.meta.keys() and self.meta['reff'] is not None:
            arr[3] = (arr[3]/self.meta['reff']).astype(self.float_dtype)
        # Calculate the radius in units of h^-1 kpc
        hkpc_per_arcsec = self.__class__._get_kpc_per_arcsec(self.meta['z'])
        arr[4] = (arr[4]*hkpc_per_arcsec).astype(self.float_dtype) if hkpc_per_arcsec > 0 \
                        else (0*arr[4]-1).astype(self.float_dtype)

        # Organize the extension data
        data = [ numpy.array(arr[0:2]).transpose(1,2,0), numpy.array(arr[2:6]).transpose(1,2,0) ] \
                    + arr[6:-1] + [ self.bin_mask.copy(), arr[-1] ]

        #---------------------------------------------------------------
        # Return the map hdus
        return DAPFitsUtil.list_of_image_hdus(data, hdr, ext)


    @staticmethod
    def _join_maps(a):
        return numpy.array(a).transpose(1,2,0)

    def stellar_continuum_maps(self, prihdr, stellar_continuum):
        """
        Construct the 'STELLAR_VEL', 'STELLAR_VEL_IVAR',
        'STELLAR_VEL_MASK', 'STELLAR_SIGMA', 'STELLAR_SIGMA_IVAR',
        'STELLAR_SIGMA_MASK', 'STELLAR_SIGMACORR',
        'STELLAR_FOM' maps extensions.
        """
        #---------------------------------------------------------------
        ext = [ 'STELLAR_VEL', 'STELLAR_VEL_IVAR', 'STELLAR_VEL_MASK', 'STELLAR_SIGMA',
                'STELLAR_SIGMA_IVAR', 'STELLAR_SIGMA_MASK', 'STELLAR_SIGMACORR', 'STELLAR_FOM' ]

        if stellar_continuum is None:
            # Construct and return the empty hdus
            return DAPFitsUtil.empty_hdus(ext)

        #---------------------------------------------------------------
        # Add data to the primary header
        prihdr = stellar_continuum._initialize_primary_header(hdr=prihdr)
        prihdr = stellar_continuum._add_method_header(prihdr)

        # Figure-of-merit units
        fomunits = ['']*9
        fomunits[0] = '1E-17 erg/s/cm^2/ang/spaxel'

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
#                DAPFitsUtil.finalize_dap_header(self.singlechannel_maphdr, 'STELLAR_SIGMACORR',
#                                                bunit='km/s'),
                DAPFitsUtil.finalize_dap_header(self.multichannel_maphdr, 'STELLAR_SIGMACORR',
                                                multichannel=True, bunit='km/s',
                                                channel_names=['resolution difference', 'fit']),
                DAPFitsUtil.finalize_dap_header(self.multichannel_maphdr, 'STELLAR_FOM',
                                                multichannel=True,
                                                channel_names=['rms', 'frms', 'rchi2',
                                                               '68th perc frac resid',
                                                               '99th perc frac resid',
                                                               'max frac resid',
                                                               '68th perc per pix chi',
                                                               '99th perc per pix chi',
                                                               'max per pix chi'],
                                                channel_units=fomunits)
              ]


        #---------------------------------------------------------------
        # Get the data arrays
        arr = [ DAPFitsUtil.redshift_to_Newtonian_velocity(
                                    stellar_continuum['PAR'].data['KIN'][:,0], self.redshift),
                DAPFitsUtil.redshift_to_Newtonian_velocity(numpy.ma.power(
                                                    stellar_continuum['PAR'].data['KINERR'][:,0],
                                                    -2).filled(0.0), self.redshift, ivar=True),
                stellar_continuum['PAR'].data['KIN'][:,1],
                numpy.ma.power(stellar_continuum['PAR'].data['KINERR'][:,1], -2).filled(0.0),
                stellar_continuum['PAR'].data['SIGMACORR_EMP'],
                stellar_continuum['PAR'].data['SIGMACORR_SRES'],
                stellar_continuum['PAR'].data['RMS'],
                stellar_continuum['PAR'].data['FRMS'],
                stellar_continuum['PAR'].data['RCHI2'],
                stellar_continuum['PAR'].data['FRMSGRW'][:,1],
                stellar_continuum['PAR'].data['FRMSGRW'][:,3],
                stellar_continuum['PAR'].data['FRMSGRW'][:,4],
                stellar_continuum['PAR'].data['CHIGRW'][:,1],
                stellar_continuum['PAR'].data['CHIGRW'][:,3],
                stellar_continuum['PAR'].data['CHIGRW'][:,4],
                stellar_continuum['PAR'].data['MASK']
              ]

        # Data types
        dtypes = [self.float_dtype]*(len(arr)-1) + [arr[-1].dtype.name]

        # Bin index
        bin_indx = stellar_continuum['BINID'].data.copy().ravel()

        # Remap the data to the DRP spatial shape
        arr = list(DAPFitsUtil.reconstruct_map(self.spatial_shape, bin_indx, arr, dtype=dtypes))
        
        # Get the masks
        vel_mask = self._stellar_continuum_mask_to_map_mask(stellar_continuum, arr[-1].copy())
        sig_mask = self._stellar_continuum_mask_to_map_mask(stellar_continuum, arr[-1].copy(),
                                                            for_dispersion=True)

        # Organize the extension data
        data = arr[0:2] + [ vel_mask ] + arr[2:4] + [ sig_mask ] \
                    + [ self._join_maps(arr[4:6]), self._join_maps(arr[6:-1]) ]

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
                'EMLINE_SEW_CNT', 'EMLINE_SEW_IVAR', 'EMLINE_SEW_MASK']

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
                DAPFitsUtil.finalize_dap_header(base_hdr, 'EMLINE_SEW_CNT',
                                                bunit='1E-17 erg/s/cm^2/ang/spaxel',
                                                multichannel=multichannel),
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
        arr += [ emission_line_moments['ELMMNTS'].data['EWCONT'][:,m]
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
                        for i in range(6) ]
        
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
        'EMLINE_GFLUX_MASK', 'EMLINE_GEW', 'EMLINE_GEW_CNT',
        'EMLINE_GEW_IVAR', 'EMLINE_GEW_MASK', 'EMLINE_GVEL',
        'EMLINE_GVEL_IVAR', 'EMLINE_GVEL_MASK', 'EMLINE_GSIGMA',
        'EMLINE_GSIGMA_IVAR', 'EMLINE_GSIGMA_MASK',
        'EMLINE_INSTSIGMA', 'EMLINE_TPLSIGMA', 'EMLINE_GA',
        'EMLINE_GANR', 'EMLINE_FOM', 'EMLINE_LFOM' map extensions.
        """
        #---------------------------------------------------------------
        ext = [ 'EMLINE_GFLUX', 'EMLINE_GFLUX_IVAR', 'EMLINE_GFLUX_MASK', 'EMLINE_GEW',
                'EMLINE_GEW_CNT', 'EMLINE_GEW_IVAR', 'EMLINE_GEW_MASK', 'EMLINE_GVEL',
                'EMLINE_GVEL_IVAR', 'EMLINE_GVEL_MASK', 'EMLINE_GSIGMA', 'EMLINE_GSIGMA_IVAR',
                'EMLINE_GSIGMA_MASK', 'EMLINE_INSTSIGMA', 'EMLINE_TPLSIGMA', 'EMLINE_GA',
                'EMLINE_GANR', 'EMLINE_FOM', 'EMLINE_LFOM' ]

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
            names = [ '{0}-{1}'.format(n,int(w)) \
                        for n,w in zip(emission_line_model['PAR'].data['NAME'],
                                       emission_line_model['PAR'].data['RESTWAVE']) ]
            base_hdr = DAPFitsUtil.add_channel_names(self.multichannel_maphdr if multichannel
                                                     else self.singlechannel_maphdr, names)
        else:
            # TODO: This should throw a ValueError, not just a warning
            warnings.warn('Emission-line model does not include channel names!')
            names = None
            base_hdr = self.multichannel_maphdr if multichannel else self.singlechannel_maphdr

        # Get the figure of merit data to output
        fom_names, fom_units, fom_data = emission_line_model.fit_figures_of_merit()
       
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
                DAPFitsUtil.finalize_dap_header(base_hdr, 'EMLINE_GEW_CNT',
                                                bunit='1E-17 erg/s/cm^2/ang/spaxel',
                                                multichannel=multichannel),
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
                                                multichannel=multichannel),
                DAPFitsUtil.finalize_dap_header(base_hdr, 'EMLINE_GA',
                                                bunit='1E-17 erg/s/cm^2/ang/spaxel',
                                                multichannel=multichannel),
                DAPFitsUtil.finalize_dap_header(base_hdr, 'EMLINE_GANR',
                                                multichannel=multichannel),
                DAPFitsUtil.finalize_dap_header(self.multichannel_maphdr, 'EMLINE_FOM',
                                                multichannel=True, channel_names=fom_names,
                                                channel_units=fom_units),
                DAPFitsUtil.finalize_dap_header(base_hdr, 'EMLINE_LFOM', multichannel=multichannel)
              ]

        #---------------------------------------------------------------
        # Get the data arrays
        n_arr_with_eml_channels = 0
        # Fluxes
        arr = [ emission_line_model['EMLDATA'].data['FLUX'][:,m]
                    for m in range(emission_line_model.neml) ]
        n_arr_with_eml_channels += 1
        # Flux errors
        arr += [ numpy.ma.power(emission_line_model['EMLDATA'].data['FLUXERR'][:,m],
                                -2).filled(0.0) for m in range(emission_line_model.neml) ]
        n_arr_with_eml_channels += 1
        # Equivalent width
        arr += [ emission_line_model['EMLDATA'].data['EW'][:,m]
                    for m in range(emission_line_model.neml) ]
        n_arr_with_eml_channels += 1
        # Equivalent width continuum
        arr += [ emission_line_model['EMLDATA'].data['EWCONT'][:,m]
                    for m in range(emission_line_model.neml) ]
        n_arr_with_eml_channels += 1
        # Equivalent width errors
        arr += [ numpy.ma.power(emission_line_model['EMLDATA'].data['EWERR'][:,m],
                                -2).filled(0.0) for m in range(emission_line_model.neml) ]
        n_arr_with_eml_channels += 1
        # Velocities
        arr += [ DAPFitsUtil.redshift_to_Newtonian_velocity(
                        emission_line_model['EMLDATA'].data['KIN'][:,m,0], self.redshift)
                    for m in range(emission_line_model.neml) ]
        n_arr_with_eml_channels += 1
        # Velocity errors
        arr += [ DAPFitsUtil.redshift_to_Newtonian_velocity(
                        numpy.ma.power(emission_line_model['EMLDATA'].data['KINERR'][:,m,0],
                                       -2).filled(0.0), self.redshift, ivar=True) \
                    for m in range(emission_line_model.neml) ]
        n_arr_with_eml_channels += 1
        # Velocity dispersions
        arr += [ emission_line_model['EMLDATA'].data['KIN'][:,m,1]
                    for m in range(emission_line_model.neml) ]
        n_arr_with_eml_channels += 1
        # Velocity dispersion errors
        arr += [ numpy.ma.power(emission_line_model['EMLDATA'].data['KINERR'][:,m,1],
                                -2).filled(0.0) for m in range(emission_line_model.neml) ]
        n_arr_with_eml_channels += 1
        # Instrumental dispersions
        arr += [ emission_line_model['EMLDATA'].data['SIGMAINST'][:,m]
                    for m in range(emission_line_model.neml) ]
        n_arr_with_eml_channels += 1
        # Template velocity dispersions
        arr += [ emission_line_model['EMLDATA'].data['SIGMATPL'][:,m]
                    for m in range(emission_line_model.neml) ]
        n_arr_with_eml_channels += 1
        # Amplitude
        arr += [ emission_line_model['EMLDATA'].data['AMP'][:,m]
                    for m in range(emission_line_model.neml) ]
        n_arr_with_eml_channels += 1
        # A/N
        arr += [ emission_line_model['EMLDATA'].data['ANR'][:,m]
                    for m in range(emission_line_model.neml) ]
        n_arr_with_eml_channels += 1
        # Line reduced chi-square
        mp = 3  # Number of model parameters, TODO: get this from emission_line_model...
        reduced_chi2 = numpy.ma.divide(emission_line_model['EMLDATA'].data['LINE_CHI2'],
                                       emission_line_model['EMLDATA'].data['LINE_NSTAT'] 
                                            - mp).filled(-999.)
        arr += [ reduced_chi2[:,m] for m in range(emission_line_model.neml) ]
        n_arr_with_eml_channels += 1
        # Full spectrum fit-quality metrics (if available)
        nfom = len(fom_names)
        arr += [ fom_data[:,m] for m in range(nfom) ]
        # Mask data
        arr += [ emission_line_model['EMLDATA'].data['MASK'][:,m]
                    for m in range(emission_line_model.neml) ]

        # Data types
        dtypes = [self.float_dtype]*(n_arr_with_eml_channels*emission_line_model.neml + nfom) \
                        + [a.dtype.name for a in arr[-emission_line_model.neml:]]

        # Bin index
        bin_indx = emission_line_model['BINID'].data.copy().ravel()

        # Remap the data to the DRP spatial shape
        arr = list(DAPFitsUtil.reconstruct_map(self.spatial_shape, bin_indx, arr, dtype=dtypes))

        # Join the maps for each emission line in the set of quantities 
        data = [ self._join_maps(arr[emission_line_model.neml*i:emission_line_model.neml*(i+1)])
                        for i in range(n_arr_with_eml_channels) ]
        data += [ self._join_maps(arr[-emission_line_model.neml-nfom:-emission_line_model.neml]) ]
        data += [ self._join_maps(arr[-emission_line_model.neml:]) ]

        # Get the masks
        base_mask = self._emission_line_model_mask_to_map_mask(emission_line_model, data[-1].copy())
        sig_mask = self._emission_line_model_mask_to_map_mask(emission_line_model, data[-1].copy(),
                                                              for_dispersion=True)

        # Organize the extension data
        data = data[:2] + [ base_mask ] + data[2:5] + [ base_mask ] + data[5:7] + [ base_mask ] \
                + data[7:9] + [ sig_mask ] + data[9:-3] + [data[-2], data[-3]]

        #---------------------------------------------------------------
        # Return the map hdus
        return DAPFitsUtil.list_of_image_hdus(data, hdr, ext)


    def spectral_index_maps(self, prihdr, spectral_indices):
        """
        Construct the 'SPECINDEX', 'SPECINDEX_IVAR', 'SPECINDEX_MASK',
        and 'SPECINDEX_CORR'.
        """
        #---------------------------------------------------------------
#        ext = [ 'SPECINDEX', 'SPECINDEX_IVAR', 'SPECINDEX_MASK', 'SPECINDEX_CORR',
#                'SPECINDEX_BCEN', 'SPECINDEX_BCNT', 'SPECINDEX_RCEN', 'SPECINDEX_RCNT',
#                'SPECINDEX_MODEL' ]
        ext = [ 'SPECINDEX', 'SPECINDEX_IVAR', 'SPECINDEX_MASK', 'SPECINDEX_CORR',
                'SPECINDEX_MODEL', 'SPECINDEX_BF', 'SPECINDEX_BF_IVAR', 'SPECINDEX_BF_MASK',
                'SPECINDEX_BF_CORR', 'SPECINDEX_BF_MODEL', 'SPECINDEX_WGT', 'SPECINDEX_WGT_IVAR',
                'SPECINDEX_WGT_MASK', 'SPECINDEX_WGT_CORR', 'SPECINDEX_WGT_MODEL']

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

        hdr = [# Worthey/Trager indices
               DAPFitsUtil.finalize_dap_header(base_hdr, 'SPECINDEX', err=True, qual=True,
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
                                              channel_names=spectral_indices['SIPAR'].data['NAME']),
               DAPFitsUtil.finalize_dap_header(base_hdr, 'SPECINDEX_MODEL',
                                               multichannel=multichannel,
                                              channel_names=spectral_indices['SIPAR'].data['NAME']),
               # Burstein/Faber indices
               DAPFitsUtil.finalize_dap_header(base_hdr, 'SPECINDEX_BF', err=True, qual=True,
                                               multichannel=multichannel,
                                              channel_names=spectral_indices['SIPAR'].data['NAME'],
                                              channel_units=spectral_indices['SIPAR'].data['UNIT']),
               DAPFitsUtil.finalize_dap_header(base_hdr, 'SPECINDEX_BF', hduclas2='ERROR', qual=True,
                                               multichannel=multichannel,
                                              channel_names=spectral_indices['SIPAR'].data['NAME'],
                                               channel_units=errunits),
               DAPFitsUtil.finalize_dap_header(base_hdr, 'SPECINDEX_BF', hduclas2='QUALITY', err=True,
                                               bit_type=self.bitmask.minimum_dtype(),
                                               multichannel=multichannel,
                                              channel_names=spectral_indices['SIPAR'].data['NAME']),
               DAPFitsUtil.finalize_dap_header(base_hdr, 'SPECINDEX_BF_CORR',
                                               multichannel=multichannel,
                                              channel_names=spectral_indices['SIPAR'].data['NAME']),
               DAPFitsUtil.finalize_dap_header(base_hdr, 'SPECINDEX_BF_MODEL',
                                               multichannel=multichannel,
                                              channel_names=spectral_indices['SIPAR'].data['NAME']),
               # Continuum in the main passband
               DAPFitsUtil.finalize_dap_header(base_hdr, 'SPECINDEX_WGT', err=True, qual=True,
                                               multichannel=multichannel,
                                              channel_names=spectral_indices['SIPAR'].data['NAME'],
                                              channel_units=spectral_indices['SIPAR'].data['UNIT']),
               DAPFitsUtil.finalize_dap_header(base_hdr, 'SPECINDEX_WGT', hduclas2='ERROR', qual=True,
                                               multichannel=multichannel,
                                              channel_names=spectral_indices['SIPAR'].data['NAME'],
                                               channel_units=errunits),
               DAPFitsUtil.finalize_dap_header(base_hdr, 'SPECINDEX_WGT', hduclas2='QUALITY', err=True,
                                               bit_type=self.bitmask.minimum_dtype(),
                                               multichannel=multichannel,
                                              channel_names=spectral_indices['SIPAR'].data['NAME']),
               DAPFitsUtil.finalize_dap_header(base_hdr, 'SPECINDEX_WGT_CORR',
                                               multichannel=multichannel,
                                              channel_names=spectral_indices['SIPAR'].data['NAME']),
               DAPFitsUtil.finalize_dap_header(base_hdr, 'SPECINDEX_WGT_MODEL',
                                               multichannel=multichannel,
                                              channel_names=spectral_indices['SIPAR'].data['NAME'])
              ]

        #---------------------------------------------------------------
        # Get the data arrays
        arr = [ spectral_indices['SINDX'].data['INDX'][:,m]
                    for m in range(spectral_indices.nindx) ]
        arr += [ numpy.ma.power(spectral_indices['SINDX'].data['INDX_ERR'][:,m], -2.).filled(0.0)
                    for m in range(spectral_indices.nindx) ]
        arr += [ spectral_indices['SINDX'].data['INDX_CORR'][:,m]
                    for m in range(spectral_indices.nindx) ]
        arr += [ spectral_indices['SINDX'].data['INDX_MOD'][:,m]
                    for m in range(spectral_indices.nindx) ]

        arr += [ spectral_indices['SINDX'].data['INDX_BF'][:,m]
                    for m in range(spectral_indices.nindx) ]
        arr += [ numpy.ma.power(spectral_indices['SINDX'].data['INDX_BF_ERR'][:,m], -2.).filled(0.0)
                    for m in range(spectral_indices.nindx) ]
        arr += [ spectral_indices['SINDX'].data['INDX_BF_CORR'][:,m]
                    for m in range(spectral_indices.nindx) ]
        arr += [ spectral_indices['SINDX'].data['INDX_BF_MOD'][:,m]
                    for m in range(spectral_indices.nindx) ]

        arr += [ spectral_indices['SINDX'].data['AWGT'][:,m]
                    for m in range(spectral_indices.nindx) ]
        arr += [ numpy.ma.power(spectral_indices['SINDX'].data['AWGT_ERR'][:,m], -2.).filled(0.0)
                    for m in range(spectral_indices.nindx) ]
        arr += [ spectral_indices['SINDX'].data['AWGT_CORR'][:,m]
                    for m in range(spectral_indices.nindx) ]
        arr += [ spectral_indices['SINDX'].data['AWGT_MOD'][:,m]
                    for m in range(spectral_indices.nindx) ]

        arr += [ spectral_indices['SINDX'].data['MASK'][:,m]
                    for m in range(spectral_indices.nindx) ]

        dtypes = [self.float_dtype]*(len(arr)-spectral_indices.nindx) \
                        + [a.dtype.name for a in arr[-spectral_indices.nindx:]]

        # Bin index
        bin_indx = spectral_indices['BINID'].data.copy().ravel()

        # Remap the data to the DRP spatial shape
        arr = list(DAPFitsUtil.reconstruct_map(self.spatial_shape, bin_indx, arr, dtype=dtypes))

        data = [ numpy.array(arr[spectral_indices.nindx*i:
                                 spectral_indices.nindx*(i+1)]).transpose(1,2,0) \
                        for i in range(13) ]

        # Get the mask
        si_mask = self._spectral_index_mask_to_map_mask(spectral_indices, data[-1].copy())
#        bf_si_mask = self._spectral_index_mask_to_map_mask(spectral_indices, data[-1].copy(),
#                                                           abs_only=True)

        # Organize the extension data
        data = data[:2] + [si_mask] + data[2:6] + [si_mask] + data[6:10] + [si_mask] + data[10:-1]

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
    def __init__(self, cube, metadata, binned_spectra=None, stellar_continuum=None,
                 emission_line_model=None, dapver=None, analysis_path=None, directory_path=None,
                 output_file=None, clobber=True, loggers=None, quiet=False,
                 single_precision=False):

        #---------------------------------------------------------------
        # Initialize the reporting
        self.loggers = None if loggers is None else loggers
        self.quiet = quiet

        self.float_dtype = 'float32' if single_precision else 'float'

        #---------------------------------------------------------------
        # Check input types
        confirm_dap_types(cube, None, binned_spectra, stellar_continuum, None,
                          emission_line_model, None)

        #---------------------------------------------------------------
        # Set the output paths
        self.cube = cube
        self.meta = metadata
        self.method = None
        self.directory_path = None
        self.output_file = None
        self._set_paths(directory_path, dapver, analysis_path, output_file, binned_spectra,
                        stellar_continuum, emission_line_model)

        # Save input for reference
        self.shape = self.cube.shape
        self.spatial_shape = self.cube.spatial_shape

        # Report
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, '-'*50)
            log_output(loggers, 1, logging.INFO, '{0:^50}'.format('CONSTRUCTING OUTPUT MODEL CUBE'))
            log_output(self.loggers, 1, logging.INFO, '-'*50)
            log_output(self.loggers, 1, logging.INFO,
                       'Output path: {0}'.format(self.directory_path))
            log_output(self.loggers, 1, logging.INFO,
                       'Output file: {0}'.format(self.output_file))
            log_output(self.loggers, 1, logging.INFO,
                       'Output cubes have shape {0}'.format(self.shape))

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
        prihdr = DAPFitsUtil.initialize_dap_primary_header(self.cube, maskname='MANGA_DAPSPECMASK')
        # Add the DAP method
        prihdr['DAPTYPE'] = (defaults.dap_method(binned_spectra.method['key'],
                                    stellar_continuum.method['fitpar']['template_library_key'],
                                    'None' if emission_line_model is None
                                        else emission_line_model.method['continuum_tpl_key']),
                             'DAP analysis method')
        # Add the format of this file
        prihdr['DAPFRMT'] = ('LOGCUBE', 'DAP data file format')

        # Get the base map header
        self.base_cubehdr = DAPFitsUtil.build_cube_header(self.cube,
                                'K Westfall <westfall@ucolick.org> & SDSS-IV Data Group',
                                                        maskname='MANGA_DAPSPECMASK')

        #---------------------------------------------------------------
        # Initialize the pixel mask
        self.bitmask = DAPCubeBitMask()

        #---------------------------------------------------------------
        # Construct the hdu list for each input object.
        # This adds the binning-specific flags to the mask
        # Binned spectra: ['FLUX', 'IVAR', 'MASK', 'WAVE', 'REDCORR']
        prihdr, binlist = self.binned_data_cube(prihdr, binned_spectra, binned_spectra_3d_hdu)

        # Model spectra: [ 'MODEL', 'MODEL_MASK', 'EMLINE', 'STELLAR', 'STELLAR_MASK' ]
        prihdr, modlist = self.model_cubes(prihdr, binned_spectra, stellar_continuum,
                                           stellar_continuum_3d_hdu, emission_line_model,
                                           emission_line_model_3d_hdu)

        # Get the BINIDs
        binidlist = combine_binid_extensions(self.cube, binned_spectra, stellar_continuum, None,
                                             emission_line_model, None, dtype='int32')

        #---------------------------------------------------------------
        # Save the data to the hdu attribute
        prihdr = finalize_dap_primary_header(prihdr, self.cube, self.meta, binned_spectra,
                                             stellar_continuum, loggers=self.loggers,
                                             quiet=self.quiet)

        # Ensure extensions are in the correct order
        self.hdu = fits.HDUList([ fits.PrimaryHDU(header=prihdr),
                                  *binlist,   # FLUX, IVAR, MASK, LSF, WAVE, REDCORR
                                  *modlist,   # MODEL, MODEL_MASK, EMLINE, STELLAR, STELLAR_MASK
                                  *binidlist  # BINID
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

    def _set_paths(self, directory_path, dapver, analysis_path, output_file, binned_spectra,
                   stellar_continuum, emission_line_model):
        """
        Set the paths relevant to the map file.
        """
        # The output method directory is, for now, the combination of
        # the binned_spectrum and stellar_continuum method keys
        if directory_path is None and (binned_spectra is None or stellar_continuum is None):
            raise ValueError('Could not define output directory path.')

        # Set the output directory path
        self.method = defaults.dap_method(binned_spectra.method['key'],
                                    stellar_continuum.method['fitpar']['template_library_key'],
                                    'None' if emission_line_model is None
                                        else emission_line_model.method['continuum_tpl_key']) \
                                if directory_path is None else None
        self.directory_path = defaults.dap_method_path(self.method, plate=self.cube.plate,
                                                       ifudesign=self.cube.ifudesign,
                                                       drpver=self.cube.drpver, dapver=dapver,
                                                       analysis_path=analysis_path) \
                                                if directory_path is None else str(directory_path)

        # Set the output file
        self.output_file = defaults.dap_file_name(self.cube.plate, self.cube.ifudesign,
                                                          self.method, mode='LOGCUBE') \
                                    if output_file is None else str(output_file)

    def _get_data_mask(self, binned_spectra, binned_spectra_3d_hdu):
        """
        For the binned spectra:
            - consolidate ANY flag except NO_STDDEV from binned_spectra
              into IGNORED; done use SpatiallyBinnedSpectra.do_not_fit_flags
            - propagate NONE_IN_STACK from binned_spectra into
              FLUXINVALID
            - propagate IVARINVALID and FORESTAR
        """
        # Initialize the mask
        mask = numpy.zeros(binned_spectra_3d_hdu['MASK'].data.shape,
                           dtype=self.bitmask.minimum_dtype())

        # Consolidate DIDNOTUSE, FORESTAR, LOW_SPECCOV, LOW_SNR from
        # binned_spectra into IGNORED
        flags = binned_spectra.do_not_fit_flags()
        indx = binned_spectra.bitmask.flagged(binned_spectra_3d_hdu['MASK'].data, flag=flags)
        mask[indx] = self.bitmask.turn_on(mask[indx], 'IGNORED')

        # Propagate NONE_IN_STACK from binned_spectra into FLUXINVALID
        indx = binned_spectra.bitmask.flagged(binned_spectra_3d_hdu['MASK'].data,
                                              flag='NONE_IN_STACK')
        mask[indx] = self.bitmask.turn_on(mask[indx], 'FLUXINVALID')

        # Propagate IVARINVALID
        indx = binned_spectra.bitmask.flagged(binned_spectra_3d_hdu['MASK'].data,
                                              flag='IVARINVALID')
        mask[indx] = self.bitmask.turn_on(mask[indx], 'IVARINVALID')

        # Propagate FORESTAR
        indx = binned_spectra.bitmask.flagged(binned_spectra_3d_hdu['MASK'].data, flag='FORESTAR')
        mask[indx] = self.bitmask.turn_on(mask[indx], 'FORESTAR')

        return mask

    def _get_stellar_continuum_mask(self, stellar_continuum, stellar_continuum_3d_hdu):
        """
        For the stellar continuum models:
            - propagate FORESTAR
            - propagate INVALID_ERROR into IVARINVALID
            - propagate ARTIFACTs
            - consolidate DIDNOTUSE, FORESTAR, LOW_SNR, ARTIFACT,
              OUTSIDE_RANGE, EML_REGION, TPL_PIXELS, TRUNCATED,
              PPXF_REJECT, INVALID_ERROR, NO_FIT into FITIGNORED
            - consolidate FIT_FAILED into FITFAILED
        """
        # Initialize
        mask = numpy.zeros(stellar_continuum_3d_hdu['MASK'].data.shape,
                           dtype=self.bitmask.minimum_dtype())

        # Propagate FORESTAR
        indx = stellar_continuum.bitmask.flagged(stellar_continuum_3d_hdu['MASK'].data,
                                                 flag='FORESTAR')
        mask[indx] = self.bitmask.turn_on(mask[indx], 'FORESTAR')

        # Propagate invalid errors
        indx = stellar_continuum.bitmask.flagged(stellar_continuum_3d_hdu['MASK'].data,
                                                 flag='INVALID_ERROR')
        mask[indx] = self.bitmask.turn_on(mask[indx], 'IVARINVALID')

        # Propagate artifacts
        indx = stellar_continuum.bitmask.flagged(stellar_continuum_3d_hdu['MASK'].data,
                                                 flag='ARTIFACT')
        mask[indx] = self.bitmask.turn_on(mask[indx], 'ARTIFACT')

        # Set which pixel were ignored in the fit
        flags = ['DIDNOTUSE', 'FORESTAR', 'LOW_SNR', 'ARTIFACT', 'OUTSIDE_RANGE', 'EML_REGION',
                 'TPL_PIXELS', 'TRUNCATED', 'PPXF_REJECT', 'INVALID_ERROR', 'NO_FIT' ]
        indx = stellar_continuum.bitmask.flagged(stellar_continuum_3d_hdu['MASK'].data, flag=flags)
        mask[indx] = self.bitmask.turn_on(mask[indx], 'FITIGNORED')

        # Set which pixels failed during the fit
        indx = stellar_continuum.bitmask.flagged(stellar_continuum_3d_hdu['MASK'].data,
                                                 flag='FIT_FAILED')
        mask[indx] = self.bitmask.turn_on(mask[indx], 'FITFAILED')

        return self._include_no_model_mask(mask)

    def _get_emission_line_model_mask(self, emission_line_model, emission_line_model_3d_hdu):
        """
        For the emission-line models:
            - propagate FORESTAR
            - propagate ARTIFACT
            - consolidate DIDNOTUSE, FORESTAR, LOW_SNR, ARTIFACT,
              OUTSIDE_RANGE, TPL_PIXELS, TRUNCATED, PPXF_REJECT,
              INVALID_ERROR into FITIGNORED
        """
        if emission_line_model is None:
            return None

        # Initialize
        mask = numpy.zeros(emission_line_model_3d_hdu['MASK'].data.shape,
                           dtype=self.bitmask.minimum_dtype())

        # Flag FORESTAR
        indx = emission_line_model.bitmask.flagged(emission_line_model_3d_hdu['MASK'].data,
                                                   flag='FORESTAR')
        mask[indx] = self.bitmask.turn_on(mask[indx], 'FORESTAR')

        # Flag artifacts
        indx = emission_line_model.bitmask.flagged(emission_line_model_3d_hdu['MASK'].data,
                                                   flag='ARTIFACT')
        mask[indx] = self.bitmask.turn_on(mask[indx], 'ARTIFACT')

        # Set which pixel were ignored in the fit
        flags = ['DIDNOTUSE', 'FORESTAR', 'LOW_SNR', 'ARTIFACT', 'OUTSIDE_RANGE', 'TPL_PIXELS',
                 'TRUNCATED', 'PPXF_REJECT', 'INVALID_ERROR' ]
        indx = emission_line_model.bitmask.flagged(emission_line_model_3d_hdu['MASK'].data,
                                                   flag=flags)
        mask[indx] = self.bitmask.turn_on(mask[indx], 'FITIGNORED')

        return self._include_no_model_mask(mask)

    def binned_data_cube(self, prihdr, binned_spectra, binned_spectra_3d_hdu):
        """
        Constructs the 'FLUX', 'IVAR', 'MASK', 'LSF', 'WAVE', and
        'REDCORR' model cube extensions, and begins the construction
        of the 'MASK'.

        Returns the primary header, a list of ImageHDUs, and the mask
        array.
        """
        #---------------------------------------------------------------
        ext = ['FLUX', 'IVAR', 'MASK', 'LSF', 'WAVE', 'REDCORR']

        if binned_spectra is None:
            # Construct and return the empty hdus
            return prihdr, DAPFitsUtil.empty_hdus(ext)

        #---------------------------------------------------------------
        # Add data to the primary header
        # TODO: Just copy header from maps file?
        prihdr = binned_spectra.rdxqa._initialize_primary_header(hdr=prihdr)
        prihdr = binned_spectra._initialize_primary_header(hdr=prihdr)
        prihdr = binned_spectra._add_method_header(prihdr)
        prihdr = binned_spectra._add_reddening_header(prihdr)

        # Save the original extension from the DRP file used to
        # construct the datacube spectral resolution data.
        lsfhdr = self.base_cubehdr.copy()

        ## KHRR added if statement:
#        if ('DRPEXT' in lsfhdr):
        if hasattr(binned_spectra.cube, 'sres_ext') and binned_spectra.cube.sres_ext is not None:
            lsfhdr['DRPEXT'] = (binned_spectra.cube.sres_ext, 'Source ext from DRP file')

        #---------------------------------------------------------------
        # Get the extension headers; The WAVE and REDCORR extensions
        # have no header.
        hdr = [DAPFitsUtil.finalize_dap_header(self.base_cubehdr, 'FLUX',
                                               bunit='1E-17 erg/s/cm^2/ang/spaxel', err=True,
                                               qual=True),
               DAPFitsUtil.finalize_dap_header(self.base_cubehdr, 'FLUX', hduclas2='ERROR',
                                               bunit='(1E-17 erg/s/cm^2/ang/spaxel)^{-2}',
                                               qual=True, prepend=False),
               DAPFitsUtil.finalize_dap_header(self.base_cubehdr, 'FLUX', hduclas2='QUALITY',
                                               err=True, bit_type=self.bitmask.minimum_dtype(),
                                               prepend=False),
               DAPFitsUtil.finalize_dap_header(lsfhdr, 'LSF', bunit='angstrom'),
               None, None]

        #---------------------------------------------------------------
        # Get the data arrays
        # Reddened flux data
        flux, ivar = binned_spectra.galext.apply(binned_spectra_3d_hdu['FLUX'].data,
                                                 deredden=False,
                                                 ivar=binned_spectra_3d_hdu['IVAR'].data)

        # Convert spectral resolution to dispersion of line-spread
        # function in angstroms
        lsf = numpy.ma.divide(binned_spectra['WAVE'].data[None,None,:],
                              binned_spectra_3d_hdu['SPECRES'].data).filled(0.0) \
                                / DAPConstants.sig2fwhm

        # Return the primary header, the list of ImageHDUs, and the mask
        data = [flux.astype(self.float_dtype), ivar.astype(self.float_dtype),
                self._get_data_mask(binned_spectra, binned_spectra_3d_hdu),
                lsf.astype(self.float_dtype),
                binned_spectra['WAVE'].data.copy().astype(self.float_dtype),
                binned_spectra['REDCORR'].data.copy().astype(self.float_dtype)]
        return prihdr, DAPFitsUtil.list_of_image_hdus(data, hdr, ext)

    def model_cubes(self, prihdr, binned_spectra, stellar_continuum, stellar_continuum_3d_hdu,
                    emission_line_model, emission_line_model_3d_hdu):
        """
        Constructs the 'MODEL', 'MODEL_MASK', 'EMLINE', 'STELLAR', and
        'STELLAR_MASK' model cube extensions, adds the model information
        to the header, and appends data to the mask.
        """
        #---------------------------------------------------------------
        ext = ['MODEL', 'MODEL_MASK', 'EMLINE', 'STELLAR', 'STELLAR_MASK']

        if binned_spectra is None or stellar_continuum is None or emission_line_model is None:
            # Construct and return empty hdus
            return prihdr, DAPFitsUtil.empty_hdus(ext)

        #---------------------------------------------------------------
        # Add the model data to the primary header
        prihdr = stellar_continuum._initialize_primary_header(hdr=prihdr)
        prihdr = stellar_continuum._add_method_header(prihdr)
        prihdr = emission_line_model._initialize_primary_header(hdr=prihdr)
        prihdr = emission_line_model._add_method_header(prihdr)

        #---------------------------------------------------------------
        # Get the extension headers
        hdr = [ DAPFitsUtil.finalize_dap_header(self.base_cubehdr, 'MODEL',
                                                bunit='1E-17 erg/s/cm^2/ang/spaxel', qual=True),
                DAPFitsUtil.finalize_dap_header(self.base_cubehdr, 'MODEL', hduclas2='QUALITY',
                                                bit_type=self.bitmask.minimum_dtype()),
                DAPFitsUtil.finalize_dap_header(self.base_cubehdr, 'EMLINE',
                                                bunit='1E-17 erg/s/cm^2/ang/spaxel'),
                DAPFitsUtil.finalize_dap_header(self.base_cubehdr, 'STELLAR',
                                                bunit='1E-17 erg/s/cm^2/ang/spaxel', qual=True),
                DAPFitsUtil.finalize_dap_header(self.base_cubehdr, 'STELLAR', hduclas2='QUALITY',
                                                bit_type=self.bitmask.minimum_dtype())
              ]

        #---------------------------------------------------------------
        # Get the data arrays
        # Best-fitting composite model
        model = stellar_continuum_3d_hdu['FLUX'].data.copy()
        if emission_line_model is not None:
            model += (emission_line_model_3d_hdu['FLUX'].data
                            + emission_line_model_3d_hdu['BASE'].data)
        # with reddening
        model = binned_spectra.galext.apply(model, deredden=False)
        if emission_line_model is None:
            line = None
        else:
            line = binned_spectra.galext.apply(emission_line_model_3d_hdu['FLUX'].data,
                                               deredden=False)
            line = line.astype(self.float_dtype)
        scnt = binned_spectra.galext.apply(stellar_continuum_3d_hdu['FLUX'].data, deredden=False)

        # Compile the data arrays
        data = [model.astype(self.float_dtype),
                self._get_emission_line_model_mask(emission_line_model,
                                                   emission_line_model_3d_hdu),
                line, scnt.astype(self.float_dtype),
                self._get_stellar_continuum_mask(stellar_continuum, stellar_continuum_3d_hdu)]

        # Return the primary header and the list of HDUs
        return prihdr, DAPFitsUtil.list_of_image_hdus(data, hdr, ext)

    def _set_no_model(self, mask):
        indx = numpy.arange(len(mask))[numpy.invert(mask > 0)]
        if len(indx) == 0:
            return mask
        mask[:indx[0]] = self.bitmask.turn_on(mask[:indx[0]], 'NOMODEL')
        mask[indx[-1]+1:] = self.bitmask.turn_on(mask[indx[-1]+1:], 'NOMODEL')
        return mask

    def _include_no_model_mask(self, mask):
        """
        Assign the NOMODEL bit to spectral regions outside of the fitted
        spectral range.

        Modeled off of StellarContiuumModel.reset_continuum_mask_window.
        """
        return mask if numpy.sum(mask > 0) == 0 else \
                    numpy.apply_along_axis(self._set_no_model, 2, mask)


#-----------------------------------------------------------------------
def combine_binid_extensions(cube, binned_spectra, stellar_continuum, emission_line_moments,
                             emission_line_model, spectral_indices, dtype=None):

    """
    Combine the bin IDs from the different analysis steps into a single
    multichannel extension.
    """
    if cube is None:
        raise ValueError('DataCube must be provided.')

    hdr = DAPFitsUtil.finalize_dap_header( DAPFitsUtil.build_map_header(cube.fluxhdr,
                                        'K Westfall <westfall@ucolick.org> & SDSS-IV Data Group',
                                                    multichannel=True, maskname='MANGA_DAPPIXMASK'),
                                          'BINID', multichannel=True,
                                          channel_names=[ 'Binned spectra', 'Stellar continua',
                                                          'Em. line moments', 'Em. line models',
                                                          'Spectral indices' ])

    # Build the ID data for each analysis product
    binid = numpy.zeros(cube.spatial_shape+(5,), dtype=int)
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


def confirm_dap_types(cube, rdxqa, binned_spectra, stellar_continuum, emission_line_moments,
                      emission_line_model, spectral_indices):

    if not isinstance(cube, DataCube):
        raise TypeError('Input must have type DataCube.')
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


# TODO: Need to abstract this for non-DRP datacubes.
def finalize_dap_primary_header(prihdr, cube, metadata, binned_spectra, stellar_continuum,
                                loggers=None, quiet=False):

    # Initialize the DAP quality flag
    dapqualbm = DAPQualityBitMask()
    dapqual = dapqualbm.minimum_dtype()(0)          # Type casting original flag to 0

    ## KHRR added if statement:
#    if isinstance(cube, MaNGADataCube):
#        if drp3qualbm.flagged(cube.prihdr['DRP3QUAL'], flag='CRITICAL'):
#            if not quiet:
#                log_output(loggers, 1, logging.INFO, 'DRP File is flagged CRITICAL!')
#            dapqual = dapqualbm.turn_on(dapqual, ['CRITICAL', 'DRPCRIT'])
    # TODO: Move this to the datacube base class?
    if cube.redux_bitmask is not None and cube.redux_qual_key in cube.prihder \
            and cube.redux_bitmask.flagged(cube.prihdr[cube.redux_qual_key],
                                           flag=self.redux_qual_flag):
        if not quiet:
            log_output(loggers, 1, logging.INFO, 'DRP File is flagged CRITICAL!')
        dapqual = dapqualbm.turn_on(dapqual, ['CRITICAL', 'DRPCRIT'])

    # Flag the file as CRITICAL if the stellar continuum fits are bad
    # for all spectra
    if stellar_continuum is not None:
        mask = stellar_continuum.bitmask.flagged(stellar_continuum['PAR'].data['MASK'],
                                                 ['NO_FIT', 'FIT_FAILED', 'INSUFFICIENT_DATA',
                                                  'NEAR_BOUND'])
        if numpy.all(mask):
            dapqual = dapqualbm.turn_on(dapqual, ['CRITICAL', 'DAPCRIT'])

    # TODO: Do similar checks for the other major modules

    # Signify that the Voronoi binning resulted in a single bin
    if binned_spectra is not None and binned_spectra.method['binclass'] is not None \
            and binned_spectra.method['binclass'].bintype == 'voronoi' \
            and binned_spectra.nbins == 1:
        dapqual = dapqualbm.turn_on(dapqual, 'SINGLEBIN')

    # Input photometric geometry and scale were invalid
    if 'ell' not in metadata.keys() or metadata['ell'] is None or \
            'pa' not in metadata.keys() or metadata['pa'] is None or \
            'reff' not in metadata.keys() or metadata['reff'] is None:
        dapqual = dapqualbm.turn_on(dapqual, 'BADGEOM')

    # Propagate some flags
    # TODO: This is a bit of a hack.  I should be more careful about how
    # "propagated" flags are defined.

    # Determine if there's a foreground star
    ## KHRR added if statement:
#    if isinstance(cube, MaNGADataCube):
#        if numpy.sum(cube.bitmask.flagged(cube.mask, flag='FORESTAR')) > 0:
#            dapqual = dapqualbm.turn_on(dapqual, 'FORESTAR')
    if cube.propagate_flags() is not None:
        flags = cube.propagate_flags()
        # TODO: This will barf if flags is a numpy array
        flags = flags if isinstance(flags, list) else [flags]
        for flag in flags:
            if numpy.any(cube.bitmask.flagged(cube.mask, flag=flag)):
                if flag in dapqual.keys():
                    dapqual = dapqualbm.turn_on(dapqual, flag)
                warnings.warn(f'{flag} is not a flag in {dapqualbm.__class__.__name__} and '
                              'cannot be propagated!')

    # Commit the quality flag to the header
    if dapqualbm.flagged(dapqual, 'CRITICAL'):
        warnings.warn('DAP files are flagged CRITICAL.')
    prihdr['DAPQUAL'] = (dapqual, 'DAP quality bitmask')

    # Finalize authors
    prihdr['AUTHOR'] = 'K Westfall <westfall@ucolick.org>'
    return prihdr


def add_snr_metrics_to_header(hdr, cube, r_re):
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
        hdr (`astropy.io.fits.Header`_):
            Header object to edit
        cube (:class:`mangadap.datacube.datacube.DataCube`):
            Datacube
        r_re (`numpy.ndarray`_):
            *Flattened* array with the semi-major axis
            radius in units of the effective radius for all spaxels
            in the datacube. Shape should be (Nx*Ny,), where the
            shape of the datacube is (Nx,Ny,Nwave).

    Returns:
        `astropy.io.fits.Header`_: The edited header object.

    Raises:
        FileNotFoundError: Raised if any of the response function files
            cannot be found.
    """
    filter_response_file = [os.path.join(defaults.dap_data_root(), 'filter_response', f) 
                                for f in ['gunn_2001_g_response.db', 'gunn_2001_r_response.db',
                                          'gunn_2001_i_response.db', 'gunn_2001_z_response.db']]
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
    flags = cube.do_not_use_flags()
#    flags = ['DONOTUSE', 'FORESTAR']
    for i in range(nfilter):
        response_func = numpy.genfromtxt(filter_response_file[i])[:,:2]
        signal, variance, snr = [a.ravel() for a in 
                                 cube.flux_stats(response_func=response_func, flag=flags)]
        covar = None
        if cube.can_compute_covariance:
            covar_wave = cube.central_wavelength(response_func=response_func, flag=flags)
            covar_channel = numpy.argsort(numpy.absolute(cube.wave-covar_wave))[0]
            covar = cube.covariance_matrix(covar_channel)

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
        covar = numpy.identity(len(variance), dtype=float) if covar is None \
                        else covar.apply_new_variance(variance).toarray()
        # Set the median S/N ...
        hdr[key_med[i]] = (numpy.ma.median(snr[indx]), com_med[i])
        # ... and the combined S/N
        hdr[key_ring[i]] = (numpy.ma.sum(signal[indx])/numpy.sqrt(numpy.sum(covar[ci,cj])),
                            com_ring[i])

    return hdr






