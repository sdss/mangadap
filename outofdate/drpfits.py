# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
r"""
Defines a class used to interface with files produced in the 3D phase of
the MaNGA Data Reduction Pipeline (DRP).

.. todo::

    - Calculation in :func:`DRPFits._cube_dimensions` will only be
      correct if the WCS coordinates have no rotation.
    - Further optimize calculation of transfer matrix
    - Make DRP file class flexible to linear or log-linear wavelength
      sampling?  Incorporate into MODE?
    - Reconstructed broad-band images and PSFs are *not* restructured in
      the CUBE files!  This is why the are transposed in
      :func:`mangadap.drpfits.DRPFits.gri_composite`.
    - Image reconstruction has transpose sense wrt DRP output!
    - Add logging

    - Need to be clear about which functions use the RSS spectra to
      create CUBE related data, like the covariance matrix and
      instrumental dispersion calculations.

    - Computing the approximate covariance cube is currently not
      possible with only the CUBE on disk.  There's a logic problem that
      needs to be fixed.

Revision history
----------------

    | **20 Nov 2014**: Original implementation by K. Westfall (KBW)
    | **12 Feb 2014**: (KBW) Added :func:`DRPFits.directory_path`
    | **20 Feb 2015**: (KBW) Add covariance calculation to :class:`DRPFits`
    | **19 Mar 2015**: (KBW) Added redux_path to :class:`DRPFits`.
        Re-arranged arguments in :func:`drpfits_list`, made drpver
        optional, and added redux_path
    | **22 May 2015**: (KBW) Sphinx documentation.  Changed DRPFits.w to
        :attr:`DRPFits.wcs`.
    | **26 May 2015**: (KBW) Added checksum=True when opening the DRP
        file.
    | **04 Jun 2015**: (KBW) Moved parse_drp_file_name to
        :func:`mangadap.util.parser.parse_drp_file_name`
    | **15 Jun 2015**: (KBW) Moved functions that return default values
        (like :func:`DRPFits._default_pixelscale`) to
        :mod:`mangadap.config.defaults`
    | **05 Aug 2015**: (KBW) Changed mode testing to be more robust.
        Added directory_path keyword to :func:`drpfits_list`.  Changed
        how directory path is set; previously required drpver and
        redux_path defaults, even if directory_path was provided
        directly.  **May need to add checks to other code to make sure
        drpver and redux_path are not None when directory_path has been
        directly defined.**
    | **28 Aug 2015**: (KBW) Added usage of
        :func:`mangadap.config.defaults.manga_fits_root`
    | **15 Feb 2016**: (KBW) Added :func:`DRPFits.__getitem__`
        function
    | **17 Feb 2016**: (KBW) Converted drpfile class name to DRPFits
    | **23 Mar 2016**: (KBW) Added functionality that abstracts the
        difference between the RSS and CUBE file formats at the user
        level.  CUBE files are now restructured to matched the
        *intended* orientation provided by the DRP; i.e., [x,y,lambda].
        RSS files are left the same, which is [fiber, lambda].
        Documentation.  Testing, particularly of x,y order.
    | **10 Nov 2016**: (KBW) Included 'DISP' in spectral arrays.
    | **30 Nov 2016**: (KBW) Include
        :func:`DRPFits.spectral_resolution`, which returns spectral
        resolution cube or vector independent of MPL
    | **01 Dec 2016**: (KBW) Added
        :func:`DRPFits.spectral_resolution_header`.
    | **06 Dec 2016**: (KBW) Removed wavelength_mask function, now uses
        :class:`mangadap.util.pixelmask.SpectralPixelMask`.  Moved the
        main functionality of :func:`DRPFits.copy_to_array` and
        :func:`DRPFits.copy_to_masked_array` to
        :class:`mangadap.util.fitsutil.DAPFitsUtil`, what's left are
        wrapper functions for the more general functions in
        :class:`mangadap.util.fitsutil.DAPFitsUtil`.
    | **17 Feb 2017**: (KBW) Return nominal inverse variance in
        :func:`DRPFits.regrid_wavelength_plane` if requested.
    | **17 May 2017**: (KBW) Include a response function in
        :func:`DRPFits.flux_stats` and
        :func:`DAPFits.mean_sky_coordinates`.
    | **21 Aug 2017**: (KBW) In spectral resolution function, select
        pre- vs. post-pixelized Gaussian respresentation.

----

.. include license and copyright
.. include:: ../copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../links.rst
"""
import time
import os.path
import warnings
import shutil

import numpy

from matplotlib import pyplot

from scipy import sparse
from scipy import interpolate
from astropy.io import fits
from astropy.wcs import WCS

try:
    from shapely.ops import cascaded_union
    from shapely.geometry import Point
except ImportError:
    warnings.warn('Could not import shapely!', ImportWarning)

from .util.fitsutil import DAPFitsUtil
from .util.bitmask import BitMask
from .util.constants import DAPConstants
from .util.parser import arginp_to_list
from .util.covariance import Covariance
from .util.pixelmask import SpectralPixelMask
from .util.filter import interpolate_masked_vector
from .util.sampling import spectral_coordinate_step
from .config import defaults

def drpfits_list(platelist, ifudesignlist, modelist, combinatorics=False, drpver=None, 
                 redux_path=None, directory_path=None):
    """
    Provided a list of plates, ifudesigns, and modes, return a list of
    DRP files to be analyzed by the DAP.

    If the number of elements in each list is the same, the matching is
    assumed to be finished unless combinatorics is True.  If the number
    of elements is not the same, or cominatorics is True, the matched
    list is all the combinations of the provided elements.

    Args:
        platelist (str or list): List of plates to use.
        ifudesignlist (str or list): List of IFU designs to use.
        modelist (str or list): List of DRP output modes ('CUBE' or 'RSS')
        combinatorics (bool): (**Optional**) Based on the input
            *platelist* and *ifudesignlist*, create a list with all
            possible combinations.
        drpver (str): (**Optional**) The DRP version, which **must** be
            a single string value used for all DRP files.
        redux_path (str): (**Optional**) The path to the top level
            directory containing the DRP output files; this is the same
            as the *redux_path* in the :class:`DRPFits` class.
        directory_path (str): (**Optional**) The exact path to the DRP
            file.  Default is defined by
            :func:`mangadap.config.defaults.drp_directory_path`.

    Returns:
        list: A list of DRP file objects

    Raises:
        Exception: Raised if *drpver* or *redux_path* are not strings.
        ValueError: Raised if the *platelist*, *ifudesignlist*, or
            *modelist* are None.

    """

    if drpver is not None and not isinstance(drpver, str):
        raise Exception('drpver must be a string')
    if redux_path is not None and not isinstance(redux_path, str):
        raise Exception('redux_path must be a string')
    if directory_path is not None and isinstance(directory_path, str):
        raise Exception('directory_path must be a string')

    if platelist is None or ifudesignlist is None or modelist is None:
        raise ValueError('Must provide platelist, ifudesignlist, and modelist!')

    platelist_ = arginp_to_list(platelist, evaluate=True)
#   print(platelist)
#   print(platelist_)
    ifudesignlist_ = arginp_to_list(ifudesignlist, evaluate=True)
#   print(ifudesignlist)
#   print(ifudesignlist_)
    modelist_ = arginp_to_list(modelist)

    n_plates = len(platelist_)
    n_designs = len(ifudesignlist_)
    n_modes = len(modelist_)

#   print(n_plates)
#   print(n_designs)
#   print(n_modes)

    # Perform the combinatorics
    if n_plates != n_designs or n_plates != n_modes or combinatorics:

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

    # Set the directory path based on the provided main path
    return [DRPFits(platelist_[i], ifudesignlist_[i], modelist_[i], drpver=drpver,
                    redux_path=redux_path, directory_path=directory_path, read=False,
                    checksum=False) \
            for i in range(0,nn)]


class DRPFitsBitMask(BitMask):
    r"""
    Structure with the DRP mask bits.A

    The defined bits are listed at :ref:`metadatamodel-drp3pixmask`.
    """
    def __init__(self, sdss_maskbits=None, mode='CUBE'):
        DRPFits.check_mode(mode)
        _sdss_maskbits = defaults.sdss_maskbits_file() if sdss_maskbits is None else sdss_maskbits
        tmp = BitMask.from_par_file(_sdss_maskbits, 'MANGA_DRP3PIXMASK' if mode == 'CUBE' else 
                                                    'MANGA_DRP2PIXMASK')
        keys, descr = tmp._init_objs()
        super(DRPFitsBitMask, self).__init__(keys, descr=descr)


class DRPQuality3DBitMask(BitMask):
    r"""
    Structure with the definition of the DRP3QUAL mask bits.

    The defined bits are listed at :ref:`metadatamodel-drp3qual`.
    """
    def __init__(self, sdss_maskbits=None):
        _sdss_maskbits = defaults.sdss_maskbits_file() if sdss_maskbits is None else sdss_maskbits
        tmp = BitMask.from_par_file(_sdss_maskbits, 'MANGA_DRP3QUAL')
        keys, descr = tmp._init_objs()
        super(DRPQuality3DBitMask, self).__init__(keys, descr=descr)


class DRPFits:
    r"""
    A general purpose class used to interface with a MaNGA DRP file.

    Args:
        plate (:obj:`int`):
            Plate number
        ifudesign (:obj:`int`):
            IFU design
        mode (:obj:`str`):
            3D mode of the DRP file; must be either 'RSS' or 'CUBE'
        drpver (:obj:`str`, optional):
            DRP version, which is used to define the default DRP
            redux path. Default is defined by
            :func:`mangadap.config.defaults.drp_version`
        redux_path (:obj:`str`, optional):
            The path to the top level directory containing the DRP
            output files for a given DRP version. Default is defined
            by :func:`mangadap.config.defaults.drp_redux_path`.
        directory_path (:obj:`str`, optional):
            The exact path to the DRP file. Default is defined by
            :func:`mangadap.config.defaults.drp_directory_path`.
        read (:obj:`bool`, optional):
            Read the DRP file upon instantiation of the object.
        checksum (:obj:`bool`, optional):
            Check for file corruption.

    Raises:
        ValueError:
            Raised if *mode* is not 'RSS' or 'CUBE'.

    Attributes:
        plate (:obj:`int`):
            Plate number
        ifudesign (:obj:`int`):
            IFU designation
        mode (:obj:`str`):
            3D mode of the DRP file, see above
        drpver (:obj:`str`):
            DRP version, which is used to define the default DRP
            redux path, see above.
        redux_path (:obj:`str`):
            The path to the top level directory containing the DRP
            output files for a given DRP version, see above.
        directory_path (:obj:`str`):
            The exact path to the DRP file, see above.
        pixelscale (:obj:`float`):
            Pixel scale used during the CUBE reconstruction.
        recenter (:obj:`bool`):
            If False, the coordinates in the XPOS and YPOS extensions
            of the DRP file are assumed to be centered at 0,0. If
            True, the XPOS and YPOS can have any center, and the
            center of the CUBE is set to be approximately the center
            of the range in XPOS,YPOS.
        width_buffer (:obj:`float`):
            The number of pixels to add to the width of the cube in
            addition to the range needed to cover XPOS and YPOS.
        xs (:obj:`float`):
            The starting on-sky coordinate (x) of the reconstructed
            image defined by the bottom corner of the first pixel,
            not its center!
        ys (:obj:`float`):
            The starting on-sky coordinate (y) of the reconstructed
            image defined by the bottom corner of the first pixel,
            not its center!
        nx (:obj:`int`):
            The size (number of pixels in x) of the reconstructed
            image.
        ny (:obj:`int`):
            The size (number of pixel in y) of the reconstructed
            image.
        
        regrid_T (`scipy.sparse.csr_matrix`_):
            Transfer matrix :math:`{\mathbf T}` such that:

            .. math::

                {\mathbf T} \times {\mathbf F} = {\mathbf I}

            where :math:`{\mathbf F}` is the vector of fluxes in a
            single wavelength channel for all the fiber measurements
            in the field of view and :math:`{\mathbf I}` is the
            pre-formatted reconstructed image of that wavelength
            channel. Saved such that repeat calls to create T for a
            given wavelength channel do not result in repeat
            calculations.
        regrid_channel (:obj:`int`):
            Wavelength channel for which :attr:`regrid_T` has been
            calculated.
        regrid_rlim (:obj:`float`):
            The limiting radius of the Gaussian interpolation kernel
            used during image construction in arcseconds.
        regrid_sigma (:obj:`float`):
            The sigma of the Gaussian interpolation kernel used
            during image construction in arcseconds.
        sigma_rho (:obj:`float`):
            The sigma, :math:`\sigma_\rho`, of the Gaussian function
            used to approximate the trend of the correlation
            coefficient :math:`\rho` with pixel separation as stored
            in :attr:`cov_rho`. That is:

            .. math::

                \rho_{ij} = \exp\left(\frac{-d^2_{ij}}{2
                \sigma^2_\rho}\right)

            where :math:`d_{ij}` is the distance between pixels
            :math:`i` and :math:`j`.
        cov_rho (`scipy.sparse.csr_matrix`_):
            The matrix :math:`{\mathbf R}` containing the correlation
            coefficents, :math:`\rho`, between elements of the
            covariance matrix, :math:`{\mathbf C}`, as approximated
            using the parameterization of :math:`\rho` with pixel
            separation. This matrix will be
            *independent of wavelength*. In general,

            .. math::

                {\mathbf C} = \left[
                \begin{array}{rrrr}
                    ... & ... & ... & ... \\
                    ... & \sigma^2_i & \rho_{ij}\sigma_i\sigma_j & ... \\
                    ... & \rho_{ji}\sigma_j\sigma_i & \sigma^2_j & ... \\
                    ... & ... & ... & ...
                \end{array}
                \right]

            such that

            .. math::

                {\mathbf R} = \left[
                \begin{array}{rrrr}
                    ... & ... & ... & ... \\
                    ... & 1.0 & \rho_{ij} & ... \\
                    ... & \rho_{ji} & 1.0 & ... \\
                    ... & ... & ... & ...
                \end{array}
                \right].

        bitmask (:class:`DRPFitsBitMask`):
            Object used to interpret the DRP bit mask values in the
            ``MASK`` extension.
        hdu (`astropy.io.fits.HDUList`_):
            HDUList read from the DRP file
        ext (:obj:`list`):
            List of fits extensions in the file
        checksum (:obj:`bool`):
            Flag to check for file corruption when opening the HDU.
        wcs (`astropy.wcs.WCS`_):
            WCS object based on WCS keywords in the header of the
            FLUX extension.
        shape (:obj:`tuple`):
            Shape of the main data arrays
        spatial_shape (:obj:`tuple`):
            Shape of the spatial axes only. For RSS files, this is a
            single element with the number of fibers; for CUBE files,
            this has the x and y dimensions of the data cube. *These
            are transposed w.r.t. the read-in DRP file!*
        nspec (:obj:`int`):
            Number of spectra in the DRP file; this is just::
                
                self.nspec = numpy.prod(self.spatial_shape)
        
        spatial_index (`numpy.ndarray`_):
            Array with tuples used to select spectra at specific
            locations within the data array. This is mainly useful in
            ``CUBE`` mode, where this provides the indices in the
            spatial coordinates. The order is :math:`(x,y)`; i.e.
            this is *different* that what you get if you read the DRP
            CUBE fits file directly using `astropy.io.fits`_. In
            ``RSS`` mode, this is just the index of the spectrum in
            the 2D array. See: :func:`select`.
        spectral_arrays (:obj:`list`):
            List of viable keywords for the data arrays in the DRP
            file. For CUBE files, these are 'FLUX', 'IVAR', and
            'MASK'; for RSS files, this also includes 'DISP', 'XPOS',
            and 'YPOS'.
        dispaxis (:obj:`int`):
            Index of the axis with the spectral channels. The
            internal data structure always has :attr:`dispaxis` as
            the *last* axis in the array. So :attr:`dispaxis` is 2
            for CUBE files and 1 for RSS files. This means that the
            internal data array restructures the input fits data for
            the CUBE files.
        nwave (:obj:`int`):
            The number of wavelength channels; this is just::

                self.nwave = self.shape[self.dispaxis]

    """
    def __init__(self, plate, ifudesign, mode, drpver=None, redux_path=None, directory_path=None,
                 read=False, checksum=False):

        # Set the attributes, forcing a known type
        self.plate = int(plate)
        self.ifudesign = int(ifudesign)
        self.check_mode(mode)
        self.mode = mode

        # Setup the directory path.
        if directory_path is None:
            self.drpver = defaults.drp_version() if drpver is None else str(drpver)
            self.redux_path = defaults.drp_redux_path(drpver=self.drpver) \
                              if redux_path is None else str(redux_path)
            self.directory_path = defaults.drp_directory_path(self.plate, drpver=self.drpver,
                                                              redux_path=self.redux_path)
        else:
            # THE DRP VERSION WILL NOT BE SPECIFIED!
            self.drpver = None
            self.redux_path = None
            self.directory_path = str(directory_path)

        # Initialize the image dimensions and their parameters
        self.pixelscale = None
        self.recenter = None
        self.width_buffer = None
        self.xs = None
        self.nx = None
        self.ys = None
        self.ny = None

        # Provide internal variables to keep the currently available
        # regridding transfer matrix for computing a plane of the data
        # cube, the index of the wavelength channel defining that plane,
        # and the properties used to create the matrix.
        self.regrid_T = None            
        self.regrid_channel = None
        self.regrid_rlim = None
        self.regrid_sigma = None

        # Provide internal variables to keep the rho matrix
        self.sigma_rho = None
        self.cov_rho = None

        # Try to define the BitMask object
        try:
            self.bitmask = DRPFitsBitMask(mode=self.mode)
        except:
            warnings.warn('Unable to define bit mask for DRP file.  Can only distinguish between'
                          'masked (values greater than 0) and unmasked (values of 0).')
            self.bitmask = None

        # Setup the variables for the internal data structure
        self.hdu = None                 # Do not automatically read the data
        self.ext = None                 # Extensions
        self.sres_ext = None            # Spectral resolution extensions
        self.checksum = checksum        # Check the file for corruption
        self.wcs = None                 # WCS structure
        self.shape = None               # Shape of the data array
        self.spatial_shape = None       # Shape of the spatial axes
        self.nspec = None               # Total number of spectra (good or bad)
        self.spatial_index = None       # Abstracted spatial indexing
        self.spectral_arrays = None     # Arrays with spectral data
        self.dispaxis = 2 if self.mode == 'CUBE' else 1
        self.nwave = None

        # Read the file, if requested
        if read:
            self.open_hdu(checksum=self.checksum)


    @classmethod
    def from_file(cls, input_file, plate=0, ifudesign=0, read=False):
        """
        Function to read a DRP-like data cube.

        .. warning::
            Work in progress.  Currently does not work.
        """
        if not os.path.isdir('./mimic_manga'):
            os.makedirs('./mimic_manga')
        dest = './mimic_manga/manga-{0}-{1}-LOGCUBE.fits.gz'.format(plate, ifudesign)
        if os.path.islink(dest):
            print('removing')
            os.remove(dest)
        os.symlink(input_file, dest)
        return cls(plate, ifudesign, 'CUBE', directory_path='./mimic_manga', read=read)


#    def __del__(self):
#        """
#        Deconstruct the DRPFits object by ensuring that the fits file is
#        properly closed.
#        """
#        if self.hdu is None:
#            return
#        self.hdu.close()
#        self.hdu = None
#        self.wcs = None
#        self.shape = None
#        self.spatial_shape = None
#        self.nspec = None
#        self.spatial_index = None


    def __getitem__(self, key):
        """Access elements of the hdu."""
        self.open_hdu(checksum=self.checksum)
        return self.hdu[key]


    def _cube_dimensions_undefined(self):
        """Return True if any of the cube dimensions are None."""
        if self.xs is None:
            return True
        if self.nx is None:
            return True
        if self.ys is None:
            return True
        if self.ny is None:
            return True
        return False


    def _cube_dimensions_correct(self, pixelscale, recenter, width_buffer):
        """
        Check that the saved parameters that define the cube dimensions
        are the same as the desired values.

        Args:
            pixelscale (float): Desired pixel scale in arcsec
            recenter (bool): Flag to recenter the coordinate system
            width_buffer (int): Number of pixels to use as buffer for
                the image reconstruction

        Returns:
            bool: True if the saved and desired values are the same.
        """
        if self.pixelscale != pixelscale:
            return False
        if self.recenter != recenter:
            return False
        if self.width_buffer != width_buffer:
            return False
        return True


    def _regrid_transfer_undefined(self):
        """Return True if :attr:`regrid_T` is None."""
        if self.regrid_T is None:
            return True
        return False


    def _regrid_kernel_correct(self, pixelscale, rlim, sigma):
        """
        Check that the saved parameters used to define the image
        reconstruction kernel are the same as the desired values.
        
        Args:
            pixelscale (float): Desired pixel scale in arcsec
            rlim (float): The limiting radius of the image
                reconstruction kernel in arcseconds.
            sigma (float): The sigma of the image reconstruction kernel
                in arcseconds.
    
        Returns:
            bool: True if the saved and desired values are the same.
        """
        if self.pixelscale != pixelscale:
            return False
        if self.regrid_rlim != rlim:
            return False
        if self.regrid_sigma != sigma:
            return False
        return True


    def _regrid_transfer_correct(self, channel, pixelscale, rlim, sigma):
        r"""
        Check that the saved parameters used to construct the transfer
        matrix, :math:`{\mathbf T}`, are the same as the desired values.

        Args:
            channel (int): Index of the wavelength channel to
                reconstruct
            pixelscale (float): Desired pixel scale in arcsec
            rlim (float): The limiting radius of the image
                reconstruction kernel in arcseconds.
            sigma (float): The sigma of the image reconstruction kernel
                in arcseconds.

        Returns:
            bool: True if the saved and desired values are the same.
        """
        if self.regrid_channel != channel:
            return False
        return self._regrid_kernel_correct(pixelscale, rlim, sigma)


    def _regrid_defaults(self, pixelscale, recenter, width_buffer, rlim, sigma):
        """
        Check that the saved regridding parameters are the same as
        the defaults. See
        :func:`mangadap.config.defaults.cube_pixelscale`,
        :func:`mangadap.config.defaults.cube_recenter`,
        :func:`mangadap.config.defaults.cube_width_buffer`,
        :func:`mangadap.config.defaults.regrid_rlim`, and
        :func:`mangadap.config.defaults.regrid_sigma`.

        Args:
            pixelscale (:obj:`float`):
                Desired pixel scale in arcsec
            recenter (:obj:`bool`):
                Flag to recenter the coordinate system
            width_buffer (:obj:`int`):
                Number of pixels to use as buffer for the image
                reconstruction
            rlim (:obj:`float`):
                The limiting radius of the image reconstruction
                kernel in arcseconds.
            sigma (:obj:`float`):
                The sigma of the image reconstruction kernel in
                arcseconds.

        Returns:
            :obj:`bool`: True if the saved and default values are the
            same.
        """
        if pixelscale != defaults.cube_pixelscale():
            return False
        if recenter != defaults.cube_recenter():
            return False
        if width_buffer != defaults.cube_width_buffer():
            return False
        if rlim != defaults.regrid_rlim():
            return False
        if sigma != defaults.regrid_sigma():
            return False
        return True

    def _init_regrid_pars(self, pixelscale, recenter, width_buffer, rlim, sigma):
        """
        Return the regridding parameters. If provided on input, the
        same value is returned. Otherwise, the returned values are
        the defaults. See the defaults defined in
        :func:`mangadap.config.defaults.cube_pixelscale`,
        :func:`mangadap.config.defaults.cube_recenter`,
        :func:`mangadap.config.defaults.cube_width_buffer`,
        :func:`mangadap.config.defaults.regrid_rlim`, and
        :func:`mangadap.config.defaults.regrid_sigma`.

        Args:
            pixelscale (:obj:`float`):
                Desired pixel scale in arcsec
            recenter (:obj:`bool`):
                Flag to recenter the coordinate system
            width_buffer (:obj:`int`):
                Number of pixels to use as buffer for the image
                reconstruction
            rlim (:obj:`float`):
                The limiting radius of the image reconstruction
                kernel in arcseconds.
            sigma (:obj:`float`):
                The sigma of the image reconstruction kernel in
                arcseconds.

        Returns:
            :obj:`tuple`: The validated regridding parameters in the
            same order as listed by the function arguments.
        """
        pixelscale = defaults.cube_pixelscale() if pixelscale is None else pixelscale
        recenter = defaults.cube_recenter() if recenter is None else recenter
        width_buffer = defaults.cube_width_buffer() if width_buffer is None else width_buffer
        rlim = defaults.regrid_rlim() if rlim is None else rlim
        sigma = defaults.regrid_sigma() if sigma is None else sigma
        return pixelscale, recenter, width_buffer, rlim, sigma

    def _variance_correlation_undefined(self):
        """Return True if :attr:`cov_rho` is None."""
        if self.cov_rho is None:
            return True
        return False

    def _variance_correlation_correct(self, sigma_rho, pixelscale, rlim, sigma):
        r"""
        Check that the saved parameters used to construct the
        correlation coefficient matrix, :math:`{\mathbf R}`, are the
        same as the desired values.
        
        Args:
            sigma_rho (:obj:`float`):
                The sigma of the Gaussian function used to
                approximate the trend of the correlation coefficient
                with pixel separation.
            pixelscale (:obj:`float`):
                Desired pixel scale in arcsec
            rlim (:obj:`float`):
                The limiting radius of the image reconstruction
                kernel in arcseconds.
            sigma (:obj:`float`):
                The sigma of the image reconstruction kernel in
                arcseconds.

        Returns:
            :obj:`bool`: True if the saved and desired values are the
            same.
        """
        if self.sigma_rho != sigma_rho:
            return False
        return self._regrid_kernel_correct(pixelscale, rlim, sigma)


    def _cube_dimensions(self, pixelscale=None, recenter=None, width_buffer=None, redo=False):
        """
        Determine the on-sky dimensions of the reconstructed image for
        all wavelength channels and save them in :attr:`xs`, :attr:`ys`,
        :attr:`nx`, and :attr:`ny`.

        For CUBE files, these dimensions are drawn directly from the WCS
        keywords in the header of the FLUX extension of the DRP fits
        file.  In this case, **any entered parameters are ignored and
        the class attributes are set to the default values used by the
        DRP**.

        .. warning::
            The calculation for the CUBE files is only valid if the WCS
            coordinate system has no rotation.

        For RSS files, the dimensions are determined using the data in
        the 'XPOS' and 'YPOS' extensions and the same algorithm used by
        the DRP; however, it is possible to provide different parameters
        that will alter the dimensions.

        See: :attr:`pixelscale`, :attr:`recenter`, :attr:`width_buffer`.

        Args:
            pixelscale (:obj:`float`, optional):
                Desired pixel scale in arcsec.
            recenter (:obj:`bool`, optional):
                Flag to recenter the coordinate system.
            width_buffer (:obj:`int`, optional):
                Number of pixels to use as buffer for the image
                reconstruction.
            redo (:obj:`bool`, optional):
                Force the recalculation of the cube dimensions if
                they are already defined.
        """

        # Make sure that the fits file is ready for reading
        self.open_hdu(checksum=self.checksum)

        # This will only be correct if the WCS coordinates have no rotation
        if self.mode == 'CUBE':
            self.pixelscale = defaults.cube_pixelscale()
            self.recenter = defaults.cube_recenter()
            self.width_buffer = defaults.cube_width_buffer()
            header = self.hdu['FLUX'].header
            # RA of first pixel edge
            self.xs = header['CRVAL1'] - header['CD1_1']*(header['CRPIX1']-1.5)
            # Offset of first pixel edge
            self.xs = (self.xs - header['OBJRA'])*numpy.cos(numpy.radians(header['OBJDEC']))*3600.
            self.nx = header['NAXIS1']

            # DEC of first pixel edge
            self.ys = header['CRVAL2'] - header['CD2_2']*(header['CRPIX2']-1.5)
            # Offset of first pixel edge
            self.ys = (self.ys - header['OBJDEC']) * 3600.
            self.ny = header['NAXIS2']
            return

        # Set the default values for the input
        if pixelscale is None:
            pixelscale = defaults.cube_pixelscale()
        if recenter is None:
            recenter = defaults.cube_recenter()
        if width_buffer is None:
            width_buffer = defaults.cube_width_buffer()

        # Check if the cube_dimensions already exist and were determined
        # using the correct parameters
        if not redo and not self._cube_dimensions_undefined() \
           and self._cube_dimensions_correct(pixelscale, recenter, width_buffer):
            return

        # Save the parameters used to create the dimensions
        self.pixelscale = pixelscale
        self.recenter = recenter
        self.width_buffer = width_buffer

        # Get the size in each dimension
        minx = numpy.amin(self.hdu['XPOS'].data)
        maxx = numpy.amax(self.hdu['XPOS'].data)
        Dx = numpy.floor(maxx-minx)

        miny = numpy.amin(self.hdu['YPOS'].data)
        maxy = numpy.amax(self.hdu['YPOS'].data)
        Dy = numpy.floor(maxy-miny)

        # Force the size to be even and the same in both dimensions
        Dx = Dx if Dx > Dy else Dy
        self.nx = int(numpy.floor(Dx/self.pixelscale)+width_buffer)
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

    def _set_variance_correlation(self, sigma_rho, pixelscale=None, recenter=None, 
                                  width_buffer=None, rlim=None, sigma=None, redo=False):
        """
        Produce :attr:`cov_rho` based on the provided *sigma_rho*.

        By default, the details of the cube dimensions should be the
        same as the DRP used to produce the 'CUBE' files from the 'RSS'
        spectra; however, these can be changed.

        The resulting :attr:`cov_rho` is independent of wavelength and
        can be used in combination with the inverse-variance produced by
        the DRP to yield a wavelength-dependent covariance matrix that
        is close to the formal calculation.

        See: :attr:`cov_rho`, attr:`sigma_rho`, :attr:`pixelscale`,
        :attr:`recenter`, :attr:`width_buffer`, :attr:`rlim`,
        :attr:`sigma`.

        .. warning::
            *sigma* is not actually used by this function.

        .. todo::
            Test the relation between *sigma*, *rlim*, and *sigma_rho*.
            It may be that sigma_rho should/can be determined by sigma.

        Args:
            sigma_rho (float): The sigma of the Gaussian function used
                to approximate the trend of the correlation coefficient
                with pixel separation.
            pixelscale (float): (**Optional**) Desired pixel scale in
                arcsec
            recenter (bool): (**Optional**) Flag to recenter the
                coordinate system
            width_buffer (int): (**Optional**) Number of pixels to use
                as buffer for the image reconstruction
            rlim (float): (**Optional**) The limiting radius of the
                image reconstruction kernel in arcseconds.
            sigma (float): (**Optional**) The sigma of the image
                reconstruction kernel in arcseconds.
            redo (bool): (**Optional**) Force the recalculation of the
                cube dimensions if they are already defined.
        """
        # Set the default values for the input
        pixelscale, recenter, width_buffer, rlim, sigma = \
            self._init_regrid_pars(pixelscale, recenter, width_buffer, rlim, sigma)

        # Check if the variance correlation coefficients already exist
        # and were determined using the correct parameters
        if not redo and not self._cube_dimensions_undefined() \
           and self._cube_dimensions_correct(pixelscale, recenter, width_buffer) \
           and not self._variance_correlation_undefined() \
           and self._variance_correlation_correct(sigma_rho, pixelscale, rlim, sigma):
            return self.cov_rho

        # Get the cube dimensions; may not necessarily match DRP calculation
        self._cube_dimensions(pixelscale=pixelscale, recenter=recenter, width_buffer=width_buffer)

        # Get the full covariance grid
        nim = self.nx*self.ny
        ii, jj = numpy.meshgrid(numpy.arange(0,nim), numpy.arange(0,nim), indexing='ij')

        # Convert covariance pixel to two image pixels for the upper
        # triangle (including the diagonal):
        i_i = numpy.zeros( (nim,nim), dtype=numpy.float64)
        i_i[jj>=ii] = ii[jj >= ii]//self.ny

        i_j = numpy.zeros( (nim,nim), dtype=numpy.float64)
        i_j[jj>=ii] = ii[jj >= ii]-i_i[jj >= ii]*self.ny

        j_i = numpy.zeros( (nim,nim), dtype=numpy.float64)
        j_i[jj>=ii] = jj[jj >= ii]//self.ny

        j_j = numpy.zeros( (nim,nim), dtype=numpy.float64)
        j_j[jj>=ii] = jj[jj >= ii]-j_i[jj >= ii]*self.ny

        #Get the distances
        dij = numpy.zeros( (nim,nim), dtype=numpy.float64)
        dij[jj>=ii] = numpy.sqrt( numpy.square(j_i[jj>=ii]-i_i[jj>=ii]) 
                                    + numpy.square(j_j[jj>=ii]-i_j[jj>=ii]) )
        dij[dij > 2*rlim/pixelscale] = 0.0

        #Set rho_ij
        rho_ij = numpy.zeros( (nim,nim), dtype=numpy.float64)
        rho_ij[dij > 0] = numpy.exp(-0.5*numpy.square(dij[dij > 0]/sigma_rho))
        rho_ij[ii==jj] = 1.0

        # Set the sparse rho matrix and save the parameters used to
        # generate it
        self.cov_rho = sparse.csr_matrix(rho_ij)
        self.sigma_rho = sigma_rho
        self.regrid_rlim = rlim
        self.regrid_sigma = sigma


    def _fix_header(self):
        """
        Use of 'degrees' in early versions of the DRP did not adhere to
        the fits standard causing the `astropy.wcs.WCS`_ to fail
        when initialized; see, e.g., :func:`world_mesh`. This changes
        the units to be 'deg' instead.

        .. note::
            This function is *obsolete* as of v1_5_1 of the DRP and is
            not actively called in this implementation; see source for
            :func:`world_mesh`.

        """
        self.open_hdu(checksum=self.checksum)
        self.hdu['FLUX'].header['CUNIT1'] = 'deg'
        self.hdu['FLUX'].header['CUNIT2'] = 'deg'


    def _set_spectral_arrays(self):
        """
        Set the list of extensions with spectral data
        """
        self.spectral_arrays = [ 'FLUX', 'IVAR', 'MASK' ]
        if self.mode == 'RSS' or (self.mode == 'CUBE' and 'LSFPOST' in self.ext):
            self.spectral_arrays += [ 'LSFPOST' ]
        if self.mode == 'RSS' or (self.mode == 'CUBE' and 'LSFPRE' in self.ext):
            self.spectral_arrays += [ 'LSFPRE' ]
        if self.mode == 'RSS':
            self.spectral_arrays += [ 'XPOS', 'YPOS' ]


    def _generate_spatial_index(self):
        """
        Generate the tuples with the list of original indices in the
        input DRP file; see :attr:`spatial_index`.
        """
        self.spatial_index = numpy.empty(self.nspec, dtype=object)
        if len(self.spatial_shape) == 1:
            self.spatial_index[:] = [ (a,) for a in numpy.arange(self.nspec) ]
            return
        i = numpy.arange(self.nspec)//self.spatial_shape[1]
        j = numpy.arange(self.nspec) - i*self.spatial_shape[1]
        # Set the coordinates to tuples in the ORIGINAL DRP fits file
        # (i.e., the DRP provides [lambda, y, x] whereas this class
        # transposes this order)
        self.spatial_index[:] = [ (ii,jj) for ii, jj in zip(i,j) ]


    def _spectral_resolution_extension(self, ext=None, pre=False):
        """
        Determine the spectral resolution channel to use.
        """
        _ext = ext
        if ext is None:
#            _ext = 'PREDISP' if pre else 'DISP'
            _ext = 'LSFPRE' if pre else 'LSFPOST'
            if _ext not in self.sres_ext:
                _ext = 'PRESPECRES' if pre else 'SPECRES'
        elif ext == 'SPECRES':
            _ext = 'PRESPECRES' if pre else 'SPECRES'
        return None if _ext not in self.sres_ext else _ext


    @staticmethod
    def mode_options():
        """
        Return the allowed modes.

        Returns:
            list: List of the allowed DRP fits file modes.
        """
        return [ 'CUBE', 'RSS' ]


    @staticmethod
    def check_mode(mode):
        """
        Check that the mode is valid.

        Args:

            mode (str): Mode value to check.  Valid modes are `CUBE` and
            `RSS`.
        """
        options = DRPFits.mode_options()
        if mode not in options:
            raise ValueError('Unknown mode {0}.  Must be in: {1}'.format(mode, options))


    @staticmethod
    def sampling_options():
        """
        Return the allowed wavelength sampling modes.

        Returns:
            list: List of the allowed DRP fits wavelength sampling
            modes.
        """
        return [ 'LIN', 'LOG' ]


    @staticmethod
    def default_paths(plate, ifudesign, mode, drpver=None, redux_path=None, directory_path=None,
                      output_file=None):
        """
        Return the primary directory and file name with the DRP fits
        LOG-binned file.

        Args:
            plate (:obj:`int`):
                Plate number of the observation.
            ifudesign (:obj:`int`):
                IFU design number of the observation.
            mode (:obj:`str`):
                3D mode of the DRP file; must be either 'RSS' or 'CUBE'
            drpver (:obj:`str`, optional):
                DRP version. Default set by
                :func:`mangadap.config.defaults.drp_version`.
            redux_path (:obj:`str`, optional):
                The path to the top level directory containing the DRP
                output files for a given DRP version.  Default is
                defined by
                :func:`mangadap.config.defaults.drp_redux_path`.
            directory_path (:obj:`str`, optional):
                The exact path to the DAP reduction assessments file.
                Default set by
                :func:`mangadap.config.defaults.dap_common_path`.
            output_file (:obj:`str`, optional):
                The name of the file with the DRP data.  Default set by
                :func:`mangadap.config.defaults.manga_fits_root`.

        Returns:
            :obj:`str`: Two strings with the path to and name of the DRP
            data file.
        """    
        _directory_path = defaults.drp_directory_path(plate, drpver=drpver,
                                                      redux_path=redux_path) \
                                if directory_path is None else directory_path
        _output_file = '{0}.fits.gz'.format(defaults.manga_fits_root(plate, ifudesign,
                                                                     'LOG{0}'.format(mode))) \
                            if output_file is None else output_file
        return _directory_path, _output_file


    def file_name(self):
        """Return the name of the DRP file"""
        root = defaults.manga_fits_root(self.plate, self.ifudesign, 'LOG{0}'.format(self.mode))
        return '{0}.fits.gz'.format(root)


    def file_path(self):
        """Return the full path to the DRP file"""
        return os.path.join(self.directory_path, self.file_name())


    def finding_chart_path(self):
        """
        Return the full path to the PNG finding chart for the targetted
        object.

        .. todo::
            - Move this to the defaults file?
        """
        return os.path.join(self.directory_path, 'images', str(self.ifudesign)+'.png')


    def open_hdu(self, checksum=False):
        """
        Read the fits file.
        
        For security purposes, the fits file is *always* opened as read
        only.  If :attr:`hdu` is not None, this function assumes the
        data has already been read and returns to the calling process.
        
        Args:
            checksum (bool): (**Optional**) Check for file corruption.

        Raises:
            FileNotFoundError: Raised if the DRP file does not exist.
        """
        if self.hdu is not None:
            return

        inp = self.file_path()
        if not os.path.exists(inp):
            raise FileNotFoundError('Cannot open file: {0}'.format(inp))

        # Open the fits file, but do NOT allow the file to be
        # overwritten.
        # TODO: This takes a while because it's restructuring ALL of the
        # image arrays.  Can I expedite this somehow?
        self.hdu = DAPFitsUtil.read(inp, permissions='readonly', checksum=checksum)
        self.ext = [ h.name for h in self.hdu ]
        self.sres_ext = [ h.name for h in self.hdu
                            if h.name in [ 'LSFPRE', 'LSFPOST', 'PRESPECRES', 'SPECRES' ] ]
        self._set_spectral_arrays()

#        # Reformat and initialize properties of the data
#        if self.mode == 'CUBE':
#            DAPFitsUtil.restructure_cube(self.hdu, ext=self.spectral_arrays)
#        elif self.mode == 'RSS':
#            DAPFitsUtil.restructure_rss(self.hdu, ext=self.spectral_arrays)

        # If RSS data, transpose the hdus so that the spectra are
        # organized along rows.
        if self.mode == 'RSS':
            self.hdu = DAPFitsUtil.transpose_image_data(self.hdu)

        self.shape = self.hdu['FLUX'].data.shape
        self.spatial_shape = DAPFitsUtil.get_spatial_shape(self.shape, self.dispaxis)
        self.nspec = numpy.prod(self.spatial_shape)
        self.nwave = self.shape[self.dispaxis]
        self._generate_spatial_index()

   
    def info(self):
        """Print the HDU info page."""
        self.open_hdu(checksum=self.checksum)
        warnings.warn('Shapes returned in \'Dimensions\' may be wrong!')
        return self.hdu.info()


    @staticmethod
    def do_not_fit_flags():
        """Return the maskbit names that should not be fit."""
        return ['DONOTUSE', 'FORESTAR']


    @staticmethod
    def do_not_stack_flags():
        """Return the maskbit names that should not be stacked."""
        return ['DONOTUSE', 'FORESTAR']


    def select(self, t, ext='FLUX', order='xy'):
        r"""
        Select a specific vector from the fits data arrays.  The spatial
        position within the original file accessed by index ``i`` is
        given by :attr:`spatial_index`.  For CUBE files this is a tuple
        with the relevant pixel indices from the original DRP file.
        This function is setup such that, e.g.::

            from mangadap.drpfits import DRPFits
            from astropy.io import fits
            import numpy
            drpf = DRPFits(7495, 1901, 'CUBE', read=True)
            hdu = fits.open(drpf.file_path())

            # should be no difference
            flux_direct = hdu['FLUX'].data[:,15,16]
            flux_class = drpf.select( (15,16), order='yx' )
            assert not numpy.any( numpy.absolute(flux_direct - flux_class) > 0 ), \
                'Selection error!'

            # Or to use more natural x y order
            flux_class = drpf.select( (16,15) )
            assert not numpy.any( numpy.absolute(flux_direct - flux_class) > 0 ), \
                'Selection error!'

            # Or select the spectrum based on its index
            for i,t in enumerate(drpf.spatial_index):
                if t == (15,16):
                    break
            flux_class = drpf.select(i)
            assert not numpy.any( numpy.absolute(flux_direct - flux_class) > 0 ), \
                'Selection error!'

            # Or use index and the provided spatial_index tuple
            i = 550
            t = drpf.spatial_index[i]
            flux_direct = hdu['FLUX'].data[:,t[0],t[1]]
            flux_class = drpf.select(i)
            assert not numpy.any( numpy.absolute(flux_direct - flux_class) > 0 ), \
                'Selection error!'
            
            flux_class = drpf.select(t)
            assert not numpy.any( numpy.absolute(flux_direct - flux_class) > 0 ), \
                'Selection error!'

        Args:
            t (tuple or int) : If an integer, this is the flattened
                index of the spectrum to return.  If a tuple, this is
                the position *within the original DRP file* for the
                non-spectral axes.  See the examples above.
            ext (str) : (**Optional**) Name of the extension from which
                to draw the data.  Must be allowed for the current
                :attr:`mode`; see :attr:`spectral_arrays`.  Default is
                ``'FLUX'``.

        Returns:
            numpy.ndarray: Vector with the data values as a function of
            wavelength.

        Raises:
            KeyError : Raised if the selected extension cannot be used.
            ValueError : Raised if the input selection object is not an
                int or a tuple.

        .. todo::
            Add select_near function that uses on-sky coordinates.

        """
        if ext not in self.spectral_arrays:
            raise KeyError('Cannot access {0} extension.'.format(ext))
        if isinstance(t, tuple):
            n = (self.shape[self.dispaxis-1],1) if self.mode == 'CUBE' else \
                    (self.shape[self.dispaxis-1],)
            _t = numpy.dot(t[::-1],n) if order == 'yx' else numpy.dot(t,n)
        elif isinstance(t, int):
            _t = t
        else:
            raise ValueError('Input must be int or tuple.')
        return self.hdu[ext].data.reshape(-1,self.nwave)[_t,:]
        

#    def wavelength_mask(self, waverange=None, toarray=False):
#        r"""
#        Return a mask boolean array with flags for pixels that fall
#        within the selected wavelength range.  If no wavelength range is
#        provided, the function just returns a fully False array with the
#        correct shape.
#
#        .. todo::
#        
#            This method is used by other similar objects; see, e.g.,
#            :class:`mangadap.proc.spatiallybinnedspectra.SpatiallyBinnedSpectra`.
#
#            These and other such methods should be pulled out into a
#            base class that these common classes are derived from.  E.g.
#            :class:`DAPFitsUtil`?
#
#        Args:
#            waverange (array-like) : (**Optional**) A two-element array
#                with the minimum and maximum wavelength to include.
#            toarray (bool) : (**Optional**) Return an array with same
#                shape as the flux array.  Default is to provide only the
#                single vector.
#
#        Returns:
#            numpy.ndarray: Boolean mask array.
#
#        Raises:
#            ValueError: Raised if the input wavelength range does not
#                have two elements.
#        """
#        outshape = self.shape if toarray else (self.nwave,)
#        if waverange is None:
#            return numpy.full(outshape, False, dtype=numpy.bool)
#        if len(waverange) != 2:
#            raise ValueError('Input wavelength range must have 2 and only 2 elements!')
#
#        if self.dispaxis is None:
#            self.dispaxis = len(outshape)-1
#        _waverange = waverange if waverange[0] < waverange[1] else [waverange[1], waverange[0]]
#        selected = (self.hdu['WAVE'].data < _waverange[0]) | (self.hdu['WAVE'].data > _waverange[1])
#        return numpy.array([selected]*self.nspec).reshape(outshape) if toarray else selected

    
    def copy_to_array(self, ext='FLUX', waverange=None):
        """
        Wrapper for :func:`mangadap.util.fitsutil.DAPFitsUtil.copy_to_array`
        specific for :class:`DRPFits`.
        """
        return DAPFitsUtil.copy_to_array(self.hdu, ext=ext, allowed_ext=self.spectral_arrays,
                                         waverange=waverange)


    def copy_to_masked_array(self, ext='FLUX', flag=None, waverange=None):
        """
        Wrapper for
        :func:`mangadap.util.fitsutil.DAPFitsUtil.copy_to_masked_array`
        specific for :class:`DRPFits`.
        """
        return DAPFitsUtil.copy_to_masked_array(self.hdu, ext=ext, mask_ext='MASK', flag=flag,
                                                bitmask=self.bitmask,
                                                allowed_ext=self.spectral_arrays,
                                                waverange=waverange)


    def spectral_resolution(self, ext=None, toarray=False, fill=False, pre=False,
                            median=False):
        """
        Return the spectral resolution at each spatial and spectral
        position.

        Select the extension 'DISP' or 'SPECRES'.  To get the
        pre-pixelized versions, set pre=True.  If you set the extension
        to 'PREDISP' and pre=True, it will try to find the extension
        'PREPREDISP' and fault.

        Args:
            ext (:obj:`str`, optional):
                Specify the extension with the spectral estimate to use.
                Should be in [ None, 'DISP', 'SPECRES'].  The default is
                None, which means it will return, in order of
                precedence, the data in 'DISP', 'SPECRES', or a None
                value if neither are present.
            toarray (:obj:`bool`, optiional):
                Return the spectral resolution as a 2D array: Nspec x
                Nwave, even if the DRP file is a CUBE object, similar to
                :func:`DRPFits.copy_to_array`.  Default is to return an
                object with the same shape as the flux array.
            fill (:obj:`bool`, optional):
                Fill masked values by interpolation.  Default is to
                leave masked pixels in returned array.
            pre (:obj:`bool`, optional):
                Read the pre-pixelized version of the spectral
                resolution, instead of the post-pixelized version.  This
                prepends 'PRE' to the extension name.
            median (:obj:`bool`, optional):
                Return a single vector with the median spectral
                resolution instead of a per spectrum array.  When using
                the `SPECRES` extension, this just returns the vector
                provided by the DRP file; when using the `DISP`
                extension, this performs a masked median across the
                array and then interpolates any wavelengths that were
                masked in all vectors.

        Returns:
            `numpy.ma.MaskedArray`_ : Even if interpolated such that
            there should be not masked values, the function returns a
            masked array.  Array contains the spectral resolution
            (:math:`R = \lambda/\Delta\lambda`) pulled from the DRP
            file.
        """
        # Make sure the fits file has been opened
        self.open_hdu(checksum=self.checksum)

        # Determine which spectral resolution element to use
        _ext = self._spectral_resolution_extension(ext=ext, pre=pre)

        # If no valid extension, raise an exception
        if ext is None and _ext is None:
            raise ValueError('No valid spectral resolution extension.')
        if ext is not None and _ext is None:
            raise ValueError('No extension: {0}'.format(ext))
            
#        # Check the selected base extension exists
#        if ext in ['DISP','SPECRES'] and ext not in self.ext:
#            raise ValueError('No extension: {0}'.format(ext))
#
#        # Set the base extension
#        _ext = ('DISP' if 'DISP' in self.ext else 'SPECRES') if ext is None else ext
#        # Add the 'PRE' qualifier if requested and check that it exists
#        if pre:
#            if 'PRE'+_ext not in self.ext:
#                raise ValueError('No {0} extension in DRP file.'.format('PRE'+_ext))
#            _ext = 'PRE'+_ext

        print('Using extension {0} to define the spectral resolution.'.format(_ext))

        # Build the spectral resolution vectors
        sres = None
        if 'LSF' in _ext:
            disp = numpy.ma.MaskedArray(self.copy_to_array(ext=_ext))
            # Mask any non-positive value
            disp[numpy.invert(disp > 0)] = numpy.ma.masked
            # Convert from sigma in angstroms to spectral resolution
            # (based on FWHM)
            sres = numpy.ma.power(DAPConstants.sig2fwhm * disp / self.hdu['WAVE'].data[None,:], -1)
        elif 'SPECRES' in _ext:
            sres = numpy.ma.MaskedArray(self.hdu[_ext].data if median
                                        else numpy.array([self.hdu[_ext].data] 
                                                            * numpy.prod(self.spatial_shape)))
            sres[numpy.invert(sres > 0)] = numpy.ma.masked
        else:
            raise ValueError('Extension {0} invalid.'.format(_ext))

        # Interpolate over any masked values
        if fill:
            outshape = sres.shape
            sres = numpy.ma.MaskedArray(
                            numpy.apply_along_axis(interpolate_masked_vector, 1,
                                                   sres.reshape(1,-1) if sres.ndim == 1
                                                        else sres.reshape(outshape[0], -1))
                                        ).reshape(outshape)

        if median and sres.ndim > 1 and sres.shape[0] > 1:
            # Determine the median over all spectra if requested and
            # necessary
            sres = numpy.ma.median(sres, axis=0)
        elif not toarray:
            # Convert back to datacube format if array format not
            # requested
            sres = sres.reshape(*self.spatial_shape,self.nwave)

        return sres


    def spectral_resolution_header(self, ext=None, pre=False):
        """
        Return a fits header for the spectral resolution array.  Copies
        the basic header from the relevant extension in the DRP file.

        Args:
            ext (str): (**Optional**) Specify the extension with the
                spectral estimate to use.  Should be in [ None, 'DISP',
                'SPECRES'].  The default is None, which means it will
                return, in order of precedence, the header for 'DISP',
                'SPECRES', or an empty header if neither are present.
            pre (bool): (**Optional**) Read the pre-pixelized version of
                the spectral resolution, instead of the post-pixelized
                version.  This prepends 'PRE' to the extension name.
        """
        self.open_hdu(checksum=self.checksum)
        if ext in ['LSFPOST','SPECRES'] and ext not in self.ext:
            raise ValueError('No extension: {0}'.format(ext))

        # Set the extension
        _ext = ('LSFPOST' if 'LSFPOST' in self.ext else 'SPECRES') if ext is None else ext

        if pre:
            if _ext == 'LSFPOST':
                _ext = 'LSFPRE'
            if _ext == 'SPECRES':
                _ext = 'PRE'+_ext
            if _ext not in self.ext:
                raise ValueError('No {0} extension in DRP file.'.format(_ext))

        return self.hdu[_ext].header.copy()
        
#        if ext in ['DISP','SPECRES'] and ext not in self.ext:
#            raise ValueError('No extension: {0}'.format(ext))
#
#        if ext == 'DISP' or (ext is None and 'DISP' in self.ext):
#            return self.hdu['DISP'].header.copy()
#        elif ext == 'SPECRES' or (ext is None and 'SPECRES' in self.ext):
#            return self.hdu['SPECRES'].header.copy()
#        return fits.Header()


    def object_data(self):
        """
        Return the MaNGA ID ('MANGAID'), object right ascension
        ('OBJRA'), and object declination ('OBJDEC') taken from the DRP
        file *primary* header.  The HDUList is opened if hasn't been
        already.

        Returns:
            str, float, float: MaNGA ID and object RA and DEC.
        """
        self.open_hdu(checksum=self.checksum)
        return (self.hdu[0].header['MANGAID'],) + self.object_coo()


    def object_coo(self):
        """
        Return the object right ascension ('OBJRA') and declination
        ('OBJDEC') taken from the DRP file *primary* header.  The
        HDUList is opened if hasn't been already.

        Returns:
            float: Object RA and DEC in degrees.
        """
        self.open_hdu(checksum=self.checksum)
        return self.hdu[0].header['OBJRA'], self.hdu[0].header['OBJDEC']


    def created_today(self):
        """
        Return True if the file was created today based on its time
        stamp.

        .. warning::
            Intended for use at the survey-level for a daily run of the
            DAP.  Likely not to be used though.
        
        """
        # Get the current time
        current_time = time.localtime()

        # Last time the DRP file was created (or modified?)
        created_time = time.gmtime(os.path.getctime(self.file_path()))

        return (current_time.tm_year == created_time.tm_year and
                current_time.tm_mon == created_time.tm_mon and 
                current_time.tm_mday == created_time.tm_mday)


    def pix_mesh(self, pixelscale=None, recenter=None, width_buffer=None, extent=False):
        r"""
        Uses `numpy.meshgrid`_ to create and return the I and J *pixel*
        coordinates for an nx*ny mesh.

        For CUBE files, these dimensions are drawn directly from the WCS
        keywords in the header of the FLUX extension of the DRP fits
        file.  In this case, **any entered parameters are ignored and
        the class attributes are set to the default values used by the
        DRP**.

        .. warning::
            The calculation for the CUBE files is only valid if the WCS
            coordinate system has no rotation.

        For RSS files, the dimensions are determined using the data in
        the 'XPOS' and 'YPOS' extensions and the same algorithm used by
        the DRP; however, it is possible to provide different parameters
        that will alter the dimensions.

        If *extent* is True, the returned mesh is of the pixel *edges*.
        That is, for a grid of size :math:`N`, there are :math:`N` pixel
        centers that run from 1 to :math:`N`; however, there are
        :math:`N+1` pixel edges that run from from :math:`0.5` to
        :math:`N+0.5`.

        See: :attr:`pixelscale`, :attr:`recenter`, :attr:`width_buffer`.

        Args:
            pixelscale (float): (**Optional**) Desired pixel scale in
                arcsec
            recenter (bool): (**Optional**) Flag to recenter the
                coordinate system
            width_buffer (int): (**Optional**) Number of pixels to use
                as buffer for the image reconstruction
            extent (bool): (**Optional**) Return a grid of the pixel
                *edges* instead of the coordinates of the pixel centers.

        Returns:
            numpy.ndarray: Two arrays with the pixel indices;
            `numpy.meshgrid`_ is run using::

                    indexing='ij'

        """
        self._cube_dimensions(pixelscale=pixelscale, recenter=recenter, width_buffer=width_buffer)
        if extent:
            return numpy.meshgrid(numpy.arange(0.0,self.nx+1)+0.5, \
                                  numpy.arange(0.0,self.ny+1)+0.5, indexing='ij')
        return numpy.meshgrid(numpy.arange(0.0,self.nx)+1, numpy.arange(0.0,self.ny)+1, \
                              indexing='ij')


    def pix_mesh_range(self, pixelscale=None, recenter=None, width_buffer=None):
        r"""
        Return the range in x and y of the reconstructed image pixels,
        including the size of the pixel.  Coordinate :math:`(1,1)` is
        the center of the first pixel, so its bottom corner is at
        :math:`(0.5,0.5)`.

        For CUBE files, these dimensions are drawn directly from the WCS
        keywords in the header of the FLUX extension of the DRP fits
        file.  In this case, **any entered parameters are ignored and
        the class attributes are set to the default values used by the
        DRP**.

        .. warning::
            The calculation for the CUBE files is only valid if the WCS
            coordinate system has no rotation.

        For RSS files, the dimensions are determined using the data in
        the 'XPOS' and 'YPOS' extensions and the same algorithm used by
        the DRP; however, it is possible to provide different parameters
        that will alter the dimensions.

        See: :attr:`pixelscale`, :attr:`recenter`, :attr:`width_buffer`.

        Args:
            pixelscale (float): (**Optional**) Desired pixel scale in
                arcsec
            recenter (bool): (**Optional**) Flag to recenter the
                coordinate system
            width_buffer (int): (**Optional**) Number of pixels to use
                as buffer for the image reconstruction
        
        Returns:
            numpy.ndarray: Two arrays with, respectively, the lower and
            upper x range and the lower and upper y range.

        """
        self._cube_dimensions(pixelscale=pixelscale, recenter=recenter, width_buffer=width_buffer)
        return numpy.array([0.5, self.nx+0.5]), numpy.array([0.5, self.ny+0.5])


    def world_mesh(self, skyright=True):
        r"""
        Return the world X and Y coordinate for the :attr:`nx` by
        :attr:`ny` mesh of the reconstructed image.  The pixel
        coordinates are first determined using :func:`pix_mesh` and then
        converted to the world coordinate system.

        For 'CUBE' files, this is done using the WCS information in the
        header of the 'FLUX' extension.  See :attr:`wcs`.

        .. note::
            Prior to v1_5_1 of the DRP, this function required a
            correction to the DRP header because of its use of
            'degrees', which does not adhere to the fits standard
            causing the `astropy.wcs.WCS`_ to fail when initialized.
            This call is not longer made.

        .. todo::
            - Make check of DRP version 'VERSDRP' to determine if header
              fix is necessary.

        For an 'RSS' file, the returned mesh is based on the offset in
        arcseconds as determined by the dimensions of the reconstructed
        image.  The *skyright* option will force the mesh to have a
        *decreasing* x coordinate as a function of increasing index.
        I.e.; the front edge of the first pixel is set to::

            x0 = self.xs+self.nx*self.pixelscale if skyright else self.xs

        and the offset per pixel is set to::

            dx = -self.pixelscale if skyright else self.pixelscale
        
        Recall that :attr:`xs` and :attr:`ys` are defined at the
        bottom-left edge of the first image pixel.

        Args:
            skyright (bool): If True, return the mesh with x coordinates
                *decreasing* as a function of increase index.

        Returns:
            numpy.ndarray: Two arrays with the world x and y
            coordinates.
        """
        if self.mode == 'RSS':
            # x and y are at the center of the pixel
            x, y = self.pix_mesh()
            # x0 is the front edge of the first pixel if not skyright an
            # the back edge of the last pixel if skyright
            x0 = self.xs+self.nx*self.pixelscale if skyright else self.xs
            dx = -self.pixelscale if skyright else self.pixelscale
            return (x-0.5)*dx+x0, (y-0.5)*self.pixelscale+self.ys

        if self.wcs is None:
            # TODO: Allow this to fix the header if MPL requires it?
#            self._fix_header()
#            print(self.hdu['FLUX'].header)
            self.wcs = WCS(header=self.hdu['FLUX'].header,fix=False,naxis=(1,2))

        x,y = self.pix_mesh()
        xy = numpy.array([x.reshape(self.nx*self.ny),y.reshape(self.nx*self.ny)]).transpose()
        XY = self.wcs.all_pix2world(xy, 1)
        return XY[:,0].reshape(self.nx, self.ny), XY[:,1].reshape(self.nx, self.ny)
#        print(self.wcs)
#        print(xy[0:self.nx,0])
#        print(xy[0:self.nx,1])
#        print(XY[0:self.nx,0])
#        print(XY[0:self.nx,1])



    def world_mesh_range(self, skyright=True):
        """
        Return the range in the world X and Y coordinates including the
        size of the pixel.
        
        For 'RSS' files, converts the result of :func:`pix_mesh_range`
        to world coordinates using the same calculation as in
        :func:`world_mesh`.

        For 'CUBE' files, converts the result of :func:`pix_mesh` with
        *extent=True* to world coordinates using the same calculation as
        in :func:`world_mesh`.  The returned arrays give the range in X
        and Y respectively.

        .. warning:: 
            May not be accurate if the reconstructed image has an on-sky
            rotation.
        
        Args:
            skyright (bool): Return the range in X such that the order
                is [max(X), min(X)] if True; otherwise, values are
                returned as [min(X), max(X)].

        Returns:
            numpy.ndarray: Two arrays with, respectively, the range in X
            and Y.

        .. note::
            Prior to v1_5_1 of the DRP, this function required a
            correction to the DRP header because of its use of
            'degrees', which does not adhere to the fits standard
            causing the `astropy.wcs.WCS`_ to fail when initialized.
            This call is not longer made.

        .. todo::
            - Make check of DRP version 'VERSDRP' to determine if header
              fix is necessary.

        """
        if self.mode == 'RSS':
            # x and y are at the edges of the pixel
            x, y = self.pix_mesh_range()
#            print(x, y)
            # x0 is the front edge of the first pixel if not skyright an
            # the back edge of the last pixel if skyright
            x0 = self.xs+self.nx*self.pixelscale if skyright else self.xs
            dx = -self.pixelscale if skyright else self.pixelscale
#            print(x0, dx)
#            print((x-0.5)*dx+x0, (y-0.5)*self.pixelscale+self.ys)
            return (x-0.5)*dx+x0, (y-0.5)*self.pixelscale+self.ys

        if self.wcs is None:
#            self._fix_header()
            self.wcs = WCS(header=self.hdu['FLUX'].header,fix=False,naxis=(1,2))

        x,y = self.pix_mesh(extent=True)
        ncoo = (self.nx+1)*(self.ny+1)
        xy = numpy.array([x.reshape(ncoo),y.reshape(ncoo)]).transpose()
        XY = self.wcs.all_pix2world(xy, 1)
        if skyright:
            return [numpy.amax(XY[:,0]),numpy.amin(XY[:,0])], \
                   [numpy.amin(XY[:,1]),numpy.amax(XY[:,1])] 

        return [numpy.amin(XY[:,0]),numpy.amax(XY[:,0])], [numpy.amin(XY[:,1]),numpy.amax(XY[:,1])] 


#    def gri_composite(self):
#        """
#        Return the world coordinates (see :func:`world_mesh`) and flux
#        in the reconstructed gri image data (in nanomaggies).  The shape
#        of the Z array is (NX,NY,3), with the g, r, and i image data in
#        Z[:,:,0], Z[:,:,1], and Z[:,:,2], respectively.
#
#        .. warning::
#            The reconstructed imaging data are not provided as part of
#            the 'RSS' files; three 'None's are returned for 'RSS' files.
#
#        Returns:
#            numpy.ndarray: Three arrays with, respectively, the world X,
#            world Y, and Z.
#        """
#
#        if self.mode == 'RSS':
#            return None, None, None
#
#        X, Y = self.world_mesh() 
#        #X,Y = self.pix_mesh()
#
#        Z = numpy.transpose(numpy.array([ self.hdu['GIMG'].data.T, self.hdu['RIMG'].data.T, \
#                            self.hdu['IIMG'].data.T ] ), axes=(1,2,0))
#
#        return X, Y, Z


    def regrid_transfer_matrix(self, channel, pixelscale=None, recenter=None, width_buffer=None,
                               rlim=None, sigma=None, quiet=False, rej_flag='3DREJECT'):
        r"""
        Calculate the transfer matrix used to produce a reconstructed
        image of the fiber data at the specified wavelength channel.
        See: :attr:`regrid_T`.

        For 'RSS' files, this is done directly using the available
        on-sky x and y coordinates of each fiber as a function of
        wavelength taken from the 'XPOS' and 'YPOS' extensions.

        For 'CUBE' files, the function will attempt to use the 'RSS'
        counterpart of the file to produce transfer matrix.  In this
        case, the input parameters *must* be the defaults.

        The default behavior is to use the same parameters as used by
        the DRP. For the defaults, see
        :func:`mangadap.config.defaults.cube_pixelscale`,
        :func:`mangadap.config.defaults.cube_recenter`,
        :func:`mangadap.config.defaults.cube_width_buffer`,
        :func:`mangadap.config.defaults.regrid_rlim`, and
        :func:`mangadap.config.defaults.regrid_sigma`. However, the
        code allows the parameters to be freely chosen by the user.

        See: :attr:`pixelscale`, :attr:`recenter`, :attr:`width_buffer`,
        :attr:`regrid_rlim`, :attr:`regrid_sigma`, :attr:`regrid_T`.

        .. todo::

            - Give more detail on the pixels at which the radius is
              calculated.
            - The details of this calculation need to be checked against
              what is done by the DRP:

                - Are the positions in XPOS sky-right in the sense that
                  positive XPOS is toward positive RA or not?
                - Does the DRP rebin the data using a kernel that
                  defines the coordinate of a pixel at the pixel center
                  or at the pixel edge?

        .. warning::
            The internal data structure of the RSS data is not different
            from simply reading the fits file using `astropy.io.fits`_.
            This means you can use this function and apply the result to
            a separate read of the RSS file (not that you would ever
            want to do that).  But after doing that, the output will be
            transposed wrt the CUBE file spatial orientation!  See
            :func:`regrid_wavelength_plane`.

        Args:
            channel (:obj:`int`):
                Index of the spectral channel for which to calculate
                the transfer matrix.
            pixelscale (:obj:`float`, optional):
                Desired pixel scale in arcsec.
            recenter (:obj:`bool`, optional):
                Flag to recenter the coordinate system
            width_buffer (:obj:`int`, optional):
                Number of pixels to use as buffer for the image
                reconstruction
            rlim (:obj:`float`, optional):
                The limiting radius of the image reconstruction
                kernel in arcseconds.
            sigma (:obj:`float`, optional):
                The sigma of the image reconstruction kernel in
                arcseconds.
            quiet (:obj:`bool`, optional):
                Suppress terminal output

        Returns:
            `scipy.sparse.csr_matrix`_ : Transfer matrix
            :math:`{\mathbf T}`

        Raises:
            ValueError:
                Raised for 'CUBE' files if the input parameters are
                not the defaults.
        """
        # Set the default values for the input
        pixelscale, recenter, width_buffer, rlim, sigma = \
            self._init_regrid_pars(pixelscale, recenter, width_buffer, rlim, sigma)

        # Check if the calculation is necessary
        if not self._regrid_transfer_undefined() \
           and self._cube_dimensions_correct(pixelscale, recenter, width_buffer) \
           and self._regrid_transfer_correct(channel, pixelscale, rlim, sigma):
            return self.regrid_T

        # Allow for calculation for CUBE files under certain conditions
        if self.mode == 'CUBE':

            # Do not perform the calculation if the parameters are not
            # the default used by the DRP to create the CUBE file.
            if not self._regrid_defaults(pixelscale, recenter, width_buffer, rlim, sigma):
                raise ValueError('Must use default pixel scale, rlim, and sigma to get transfer '
                                 + 'matrix for DRP-produced CUBE files.')

            warnings.warn('Attempting to use RSS counter-part for calculation.')
            # Get the RSS counterpart
            drpf = DRPFits(self.plate, self.ifudesign, 'RSS', drpver=self.drpver, \
                           redux_path=self.redux_path, directory_path=self.directory_path)
            # Get the transfer matrix
            self.regrid_T = drpf.regrid_transfer_matrix(channel)
            self.regrid_channel = drpf.regrid_channel
            # Save the parameters (should be the defaults)
            self.pixelscale = drpf.pixelscale
            self.regrid_rlim = drpf.regrid_rlim
            self.regrid_sigma = drpf.regrid_sigma
            # Also need to make sure dimensions are set
            self._cube_dimensions()
        
            return self.regrid_T

        # Get the cube dimensions; may not necessarily match DRP calculation
        self._cube_dimensions(pixelscale=pixelscale, recenter=recenter, width_buffer=width_buffer)

        # Dimensions of the sparse matrix are
        nim = self.nx*self.ny                   # The number of image pixels
        # by the number of fiber spectra (self.nspec)

        # Get the list of non-zero pixel values in the transfer matrix
        i = numpy.arange(self.nx)
        j = numpy.arange(self.ny)
        ii,jj = numpy.meshgrid(i, j, indexing='ij')         # Mesh of i,j pixel indices

        sp = numpy.empty((self.nx,self.ny), dtype=numpy.float64)    # Holds spectrum index
        ij = (ii*self.ny+jj)                                        # Holds image pixel index
        r2 = numpy.empty((self.nx,self.ny), dtype=numpy.float64)    # Holds radii
        tot = numpy.zeros((self.nx,self.ny), dtype=numpy.float64)   # Holds the sum of the weights

        s2 = numpy.square(sigma/pixelscale)                 # sigma^2 of Gaussian
        rl2 = numpy.square(rlim/pixelscale)                 # radius^2 of Gaussian limit

        non_zero_spc = numpy.empty(0, dtype=numpy.int64)    # Holds triplet spectrum index
        non_zero_pix = numpy.empty(0, dtype=numpy.int64)    # Holds triplet image index
        non_zero_wgt = numpy.empty(0, dtype=numpy.float64)  # Holds triplet weight

        # Do not include any pixels with zero inverse variance or pixels
        # that have been flagged with the provided mask bits
        mask = numpy.invert(self.hdu['IVAR'].data[:,channel] > 0.0)
        if rej_flag is not None:
            _rej_flag = rej_flag if isinstance(rej_flag, list) or rej_flag != 'any' else None
            mask |= self.bitmask.flagged(self.hdu['MASK'].data[:,channel], flag=_rej_flag)

        # TODO: Can optimize this further
        for k in range(self.nspec):
            if mask[k]:
                continue

            # Fill spectrum index
            sp.fill(k)

            # NOTE: Calculating full matrix is actually faster than
            # determining submatrix for calculation

            # Calcuate the distance
            # ---- WITH RESPECT TO THE EDGE OF THE FIRST PIXEL ----
            #  - matches DRP, but why?!?!
            r2 = numpy.square( (self.hdu['XPOS'].data[k,channel]-self.xs)/pixelscale - ii) \
                 + numpy.square((self.hdu['YPOS'].data[k,channel]-self.ys)/pixelscale - jj)
            # ---- WITH RESPECT TO THE CENTER OF THE FIRST PIXEL ----
#           r2 = numpy.square( (self.hdu['XPOS'].data[k,channel]-self.xs)/pixelscale-0.5 - ii) \
#                + numpy.square((self.hdu['YPOS'].data[k,channel]-self.ys)/pixelscale-0.5 - jj)

            # Append new indices and weights within rlim
            # TODO: Can this be done quicker if I'm not appending things?
            non_zero_spc = numpy.append(non_zero_spc, sp[r2 < rl2])
            non_zero_pix = numpy.append(non_zero_pix, ij[r2 < rl2])
            wgt = numpy.exp(-r2[r2 < rl2]/s2/2.0)
            tot[r2<rl2] += wgt
            non_zero_wgt = numpy.append(non_zero_wgt, wgt)

            if not quiet:
                print('Transfer Matrix {:2.1%}'.format((k+1)/self.nspec), end="\r")

        if not quiet:
            print('Transfer Matrix Done                     ')

        # Save the regridding input
        self.regrid_rlim = rlim
        self.regrid_sigma = sigma
        self.regrid_channel = channel
        self.pixelscale = pixelscale

        # Normalize the result and scale by the pixel size to ensure the
        # output cube is in units of calibrated flux density per pixel
        scale = pixelscale*pixelscale/numpy.pi
        non_zero_wgt *= scale/tot[numpy.unravel_index(non_zero_pix.astype(int), (self.nx,self.ny))]

        # Set the transfer matrix to a sparse object and return it
        self.regrid_T = sparse.coo_matrix( (non_zero_wgt, (non_zero_pix, non_zero_spc)), \
                                                   shape=(nim,self.nspec) ).tocsr()
        return self.regrid_T


    def regrid_wavelength_plane(self, channel, pixelscale=None, recenter=None, width_buffer=None,
                                rlim=None, sigma=None, quiet=False, return_ivar=False,
                                return_covar=False):
        r"""
        Return the reconstructed image for the specified wavelength
        channel.
        
        For 'CUBE' files, the input parameters *must* be the same as the
        defaults.  If they are, the function simply returns the selected
        2D array from the 'FLUX' extension.

        For 'RSS' files, the transfer matrix is first calculated using
        :func:`regrid_transfer_matrix` and then used to calculate:

        .. math::
        
            {\mathbf T} \times {\mathbf F} = {\mathbf I}

        where :math:`{\mathbf F}` is the vector of fluxes in the
        selected wavelength channel for all the fiber measurements in
        the field of view and :math:`{\mathbf I}` is the pre-formatted
        (flattened) reconstructed image of that wavelength channel.

        **On output, :math:`{\mathbf I}` is rearranged into a 2D array
        of size :attr:`nx` by :attr:`ny`.**

        .. warning::
            Because the internal data structure is transposed with
            respect to the DRP file, the result is transposed on output
            so that it will match the similar operation done directly on
            a DRP file.
        
        Args:
            channel (int): Index of the spectral channel for which
                to calculate the transfer matrix.
            pixelscale (float): (**Optional**) Desired pixel scale in
                arcsec
            recenter (bool): (**Optional**) Flag to recenter the
                coordinate system
            width_buffer (int): (**Optional**) Number of pixels to use
                as buffer for the image reconstruction
            rlim (float): (**Optional**) The limiting radius of the
                image reconstruction kernel in arcseconds.
            sigma (float): (**Optional**) The sigma of the image
                reconstruction kernel in arcseconds.
            quiet (bool): (**Optional**) Suppress terminal output
            return_ivar (bool): (**Optional**) Return the nominal
                inverse variance image (i.e. the diagonal of the
                covariance matrix).

        Returns:
            numpy.ndarray: The reconstructed image for the specified
            wavelength channel.

        Raises:
            Exception: Raised for 'CUBE' files when the input parameters
                are not the same as the defaults.
        
        """
        # Set the default values for the input
        pixelscale, recenter, width_buffer, rlim, sigma = \
            self._init_regrid_pars(pixelscale, recenter, width_buffer, rlim, sigma)

        # Allow CUBE output under certain conditions
        if self.mode == 'CUBE':

            # Do not perform the calculation if the parameters are not
            # the default used by the DRP to create the CUBE file.
            if not self._regrid_defaults(pixelscale, recenter, width_buffer, rlim, sigma):
                raise ValueError('Must use default pixel scale, rlim, and sigma to get '
                                 + 'wavelength-channel image for DRP-produced CUBE files.')

            self.open_hdu(checksum=self.checksum)
            return self.hdu['FLUX'].data[:,:,channel].copy()

        # Set the transfer matrix (set to self.regrid_T; don't need to
        # keep the returned matrix)
        self.regrid_transfer_matrix(channel, pixelscale, recenter, width_buffer, rlim, sigma, quiet)

        flux = self.regrid_T.dot(self.hdu['FLUX'].data[:,channel]).reshape(self.nx,self.ny)
        # Return the regridded data with the proper shape (nx by ny)
        if not return_ivar:
            return flux

        covar = self._formal_covariance_matrix(channel, pixelscale, recenter, width_buffer, rlim,
                                               sigma)
        if return_covar:
            return flux, covar
            
        return flux, numpy.ma.power(numpy.diagonal(covar.toarray()), \
                                    -1).filled(0.0).reshape(self.nx,self.ny)


    def _formal_covariance_matrix(self, channel, pixelscale, recenter, width_buffer, rlim, sigma,
                                  csr=False, quiet=False):
        r"""
        Return the formal covariance matrix as defined by

        .. math::
        
             {\mathbf C} = {\mathbf T} \times {\mathbf \Sigma} \times
             {\mathbf T}^{\rm T},

        where :math:`{\mathbf \Sigma}` is the covariance matrix for the
        'RSS' spectra for the specified wavelength channel.  For 'CUBE'
        files, the function will attempt to use the 'RSS' counterpart of
        the file to produce the transfer matrix, :math:`{\mathbf T}`.
        In this case, the input parameters *must* be the defaults.

        .. note::
        
            The current DRP does not produce spectral covariance
            matrices for the 'RSS' spectra.  Here, it is assumed that
            the spectral covariance matrix is zero everywhere except
            along the diagonal, which contains the inverse of the values
            in the 'IVAR' extension.

        The use of this function should **not** be used with the 'CUBE'
        files for efficiency sake.

        Args:
            channel (int): Index of the spectral channel for which to
                calculate the transfer matrix.
            pixelscale (float): Desired pixel scale in arcsec
            recenter (bool): Flag to recenter the coordinate system
            width_buffer (int): Number of pixels to use as buffer for
                the image reconstruction
            rlim (float): The limiting radius of the image
                reconstruction kernel in arcseconds.
            sigma (float): The sigma of the image reconstruction kernel
                in arcseconds.
            csr (bool): (**Optional**) Instead of reaturning a
                :class:`mangadap.util.covariance.Covariance` object,
                return the covariance matrix as a
                `scipy.sparse.csr_matrix`_ object.  Primarily used by
                :func:`covariance_cube` for collating the covariance
                matrix of each wavelength channel before combining them
                into a single
                :class:`mangadap.util.covariance.Covariance` object
            quiet (bool): (**Optional**) Suppress terminal output

        Returns:
            :class:`mangadap.util.covariance.Covariance` or
            `scipy.sparse.csr_matrix`_: The covariance matrix for the
            designated wavelength channel.  The return type depends on
            *csr*.

        Raises:
            ValueError: Raised for 'CUBE' files when the input
                parameters are not the same as the defaults.
        """
        # Allow CUBE output under certain conditions
        if self.mode == 'CUBE':

            # Do not perform the calculation if the parameters are not
            # the default used by the DRP to create the CUBE file.
            if not self._regrid_defaults(pixelscale, recenter, width_buffer, rlim, sigma):
                raise ValueError('Must use default pixel scale, rlim, and sigma to get '
                                 + 'covariance matrices for DRP-produced CUBE files.')

            warnings.warn('Attempting to use RSS counter-part for calculation.')
            drpf = DRPFits(self.plate, self.ifudesign, 'RSS', drpver=self.drpver, \
                           redux_path=self.redux_path, directory_path=self.directory_path)
            return drpf.covariance_matrix(channel, pixelscale=pixelscale, recenter=recenter,
                                          width_buffer=width_buffer, rlim=rlim, sigma=sigma,
                                          csr=csr, quiet=quiet)

        # Set the transfer matrix (set to self.regrid_T; don't need
        # to keep the returned matrix)
        self.regrid_transfer_matrix(channel, pixelscale, recenter, width_buffer, rlim, sigma,
                                    quiet)

        # Get the variance values, ignoring those that are <= 0
#        var = numpy.zeros(self.nspec, dtype=numpy.float64)
#        indx = numpy.where(self.hdu['IVAR'].data[:,channel] > 0.)
#        var[indx] = 1.0/self.hdu['IVAR'].data[indx,channel]

        # TODO: Change to using Covariance.from_matrix_multiplication

        # Set the covariance matrix of the spectra to be a diagonal
        # matrix with the provided variances
        Sigma = sparse.coo_matrix((numpy.ma.power(self.hdu['IVAR'].data[:,channel],
                                                  -1.0).filled(0.0),
                                    (numpy.arange(0,self.nspec),numpy.arange(0,self.nspec))), \
                                  shape=(self.nspec, self.nspec)).tocsr()

        # Return the covariance matrix from the spatial rebinning
        C = sparse.triu(self.regrid_T.dot(Sigma.dot(self.regrid_T.transpose()))).tocsr()

        return Covariance(C) if not csr else C


    def _approximate_covariance_matrix(self, channel, pixelscale, recenter, width_buffer, rlim,
                                       sigma, sigma_rho, csr=False, quiet=False):
        r"""
        Return an approximate calculation of the covariance matrix
        assuming
        
        .. math::

            C_{ij} = \frac{\rho_{ij}}{(V^{-1}_{ii} V^{-1}_{jj})^{1/2}}
            
        where :math:`\rho_{ij}` is approximated by a Gaussian with a
        standard deviation defined by the provided *sigma_rho*.  See the
        description of attributes :attr:`cov_rho` and :attr:`sigma_rho`.
        For this to work, the variance matrix of the reconstructed image
        must have been already calculated, meaning that this approach is
        most appropriately used with the 'CUBE' files.  For 'RSS' files,
        the covariance matrix is determined in this case by calling
        :func:`covariance_matrix` on the 'CUBE' file, which will raise
        an exception if the parameters defining the dimensions of the
        reconstructed image and kernel are not the defaults.

        In general, this function should **not** be used with 'RSS'
        files for efficiency sake.
        
        Args:
            channel (int): Index of the spectral channel for which to
                calculate the transfer matrix.
            pixelscale (float): Desired pixel scale in arcsec
            recenter (bool): Flag to recenter the coordinate system
            width_buffer (int): Number of pixels to use as buffer for
                the image reconstruction
            rlim (float): The limiting radius of the image
                reconstruction kernel in arcseconds.
            sigma (float): The sigma of the image reconstruction kernel
                in arcseconds.
            sigma_rho (float): The sigma of the Gaussian function used
                to approximate the trend of the correlation coefficient
                with pixel separation.
            csr (bool): (**Optional**) Instead of reaturning a
                :class:`mangadap.util.covariance.Covariance` object,
                return the covariance matrix as a
                `scipy.sparse.csr_matrix`_ object.  Primarily used by
                :func:`covariance_cube` for collating the covariance
                matrix of each wavelength channel before combining them
                into a single
                :class:`mangadap.util.covariance.Covariance` object
            quiet (bool): (**Optional**) Suppress terminal output

        Returns:
            :class:`mangadap.util.covariance.Covariance` or
            `scipy.sparse.csr_matrix`_: The covariance matrix for the
            designated wavelength channel.  The return type depends on
            `csr`.

        Raises:
            ValueError: Raised for 'CUBE' files when the input
                parameters are not the same as the defaults.
        """

        # Allow to process RSS if necessary, but warn user
        if self.mode == 'RSS':
            warnings.warn('Attempting to use CUBE counter-part for calculation.')
            drpf = DRPFits(self.plate, self.ifudesign, 'CUBE', drpver=self.drpver, \
                           redux_path=self.redux_path, directory_path=self.directory_path)
            return drpf.covariance_matrix(channel, pixelscale, recenter, width_buffer, rlim,
                                          sigma, sigma_rho, csr, quiet)

        # Get the variance correlation (rho) matrix (returns
        # existing matrix if available for same parameters)
        self._set_variance_correlation(sigma_rho, pixelscale=pixelscale, recenter=recenter,
                                       width_buffer=width_buffer, rlim=rlim, sigma=sigma)

        # Get the non-zero elements
        ci, cj, rho = sparse.find(self.cov_rho)

        # Get the cube pixels
        i_i = ci//self.ny
        i_j = ci-i_i*self.ny
        j_i = cj//self.ny
        j_j = cj-j_i*self.ny

        # Use the available inverse variance cube to approximately
        # calculate the full covariance matrix
        cov = numpy.sqrt(self.hdu['IVAR'].data[i_i,i_j,channel]
                         * self.hdu['IVAR'].data[j_i,j_j,channel])
        cov[cov>0] = rho[cov>0]/cov[cov>0]

        C = sparse.coo_matrix((cov[cov>0], (ci[cov>0], cj[cov>0])), 
                              shape=(self.nx*self.ny,self.nx*self.ny)).tocsr()
        return Covariance(C) if not csr else C


    def _interpolated_response_function(self, response_func):
        if response_func is None:
            return numpy.ones(self.nwave, dtype=float)
        
        interp = interpolate.interp1d(response_func[:,0], response_func[:,1], bounds_error=False,
                                      fill_value=0.0, assume_sorted=True)
        return interp(self['WAVE'].data)


    def _covariance_wavelength(self, waverange=None, response_func=None, per_pixel=True, flag=None):
        """
        Determine the wavelength at which to calculate the covariance
        matrix.

        Args:
            waverange (array-like): (**Optional**) Starting and ending
                wavelength over which to calculate the statistics.
                Default is to use the full wavelength range.
            response_func (array-like): (**Optional**) A two-column
                array with the wavelength and transmission of a
                broad-band response function to use for the calculation.
            per_pixel (bool): (**Optional**) When providing a response
                function, continue to calculate the statistics per
                pixel, as opposed to per angstrom.  Default is to
                compute the statistics on a per pixel basis.
            flag (str or list): (**Optional**) (List of) Flag names that
                are considered when deciding if a pixel should be
                masked.  The names *must* be a valid bit name as defined
                by :attr:`bitmask` (see :class:`DRPFitsBitMask`).

        Returns:
            float: The response-weighted center of the wavelength region
            used to calculate the S/N, which will be where the
            covariance matrix is calculated.

        """
        if waverange is None and response_func is None:
            return (self['WAVE'].data[0]+self['WAVE'].data[-1])/2.

        wave = self['WAVE'].data
        flux = self.copy_to_masked_array(waverange=waverange, flag=flag)

        dw = numpy.ones(self.nwave, dtype=float) if per_pixel else \
                    spectral_coordinate_step(wave, log=True)*numpy.log(10.)*wave

        _response_func = self._interpolated_response_function(response_func)
        response_integral = numpy.sum(numpy.invert(numpy.ma.getmaskarray(flux))
                                        *_response_func[None,:]*dw[None,:], axis=1)
        return numpy.ma.mean(numpy.ma.divide(numpy.sum(numpy.invert(numpy.ma.getmaskarray(flux))
                                                *(self['WAVE'].data* _response_func*dw)[None,:],
                                                       axis=1), response_integral))


    def covariance_matrix(self, channel, pixelscale=None, recenter=None, width_buffer=None, 
                          rlim=None, sigma=None, sigma_rho=None, csr=False, quiet=False):
        r"""
        Return the covariance matrix for the specified wavelength
        channel.

        For a regrided cube image with :math:`N_x\times N_y` pixels, the
        covariance matrix has a size that is :math:`(N_x N_y)\times (N_x
        N_y)`; however, the majority of these pixels will be zero.
        Therefore, the covariance matrix is stored as a sparse matrix
        and interfaced with using the
        :class:`mangadap.util.covariance.Covariance` object class.

        The value of the covariance matrix at pixel :math:`(i,j)` is the
        covariance between pixels :math:`(n_{x,0},n_{y,0})` and
        :math:`(n_{x,1},n_{y,1})` at the specified wavelength channel of
        the reconstructed CUBE image, where

        .. math::

            n_{x,0} &= \lfloor i / N_y \rfloor \\
            n_{y,0} &= i - n_{x,0} N_y \\
            n_{x,1} &= \lfloor j / N_y \rfloor \\
            n_{y,1} &= j - n_{x,1} N_y

        and :math:`\lfloor m\rfloor` is the "floor" of :math:`m`.  The
        diagonal of the covariance matrix (:math:`i=j`) should directly
        provide the inverse of the IVAR values provided by the DRP.

        .. warning::
            THIS IS IMPORTANT.  Because the sense of the pixel
            coordinates, :math:`(n_x,n_y)`, is flipped with respect to a
            direct read of the DRP fits file using `astropy.io.fits`_,
            it is important that you appreciate this ordering in
            handling the coordinates in the output covariance matrix
            using the equation above.  I.e., if you want to get the
            covariance between two pixels in the DRP file, make sure you
            understand which is the :math:`x` pixel and which is the
            :math:`y` pixel!

        You can compare the CUBE provide inverse variance matrix and the
        inverse variance matrix provided along the diagonal of a
        covariance matrix calculated using this function as follows::

            from matplotlib import pyplot
            from mangadap.drpfits import DRPFits
            from astropy.io import fits

            drpf = DRPFits(7495, 3701, 'CUBE', read=False)
            hdu = fits.open(drpf.file_path())
            drpf = DRPFits(7495, 3701, 'RSS', read=True)

            channel = 2000
            C = drpf.covariance_matrix(channel)

            ivar = numpy.diag(C.toarray()).reshape(-1, \
                                    numpy.sqrt(C.shape[0])).T.copy()
            indx = ivar > 0
            ivar[indx] = 1./ivar[indx]
            ivar[numpy.invert(indx)] = 0.

            fig = pyplot.figure(figsize=(pyplot.figaspect(1)))
            ax = fig.add_axes([0.05, 0.3, 0.35, 0.35])
            ax.imshow(ivar, origin='lower', interpolation='nearest',
                      cmap='inferno', aspect='auto')
            ax.set_title('regrid')
            ax = fig.add_axes([0.5, 0.3, 0.35, 0.35])
            cax = fig.add_axes([0.86, 0.3, 0.01, 0.35])
            ax.set_title('DRP')
            cs = ax.imshow(hdu['IVAR'].data[channel,:,:], origin='lower',
                           interpolation='nearest', cmap='inferno')
            pyplot.colorbar(cs, cax=cax)
            pyplot.show()

        Differences between the inverse variance values should be small
        (:math:`\approx 10^{-3}`, which is about 1 part in
        :math:`10^7`).

        The covariance matrix can be generated in one of two ways:

        1. If *sigma_rho* is None, the returned matrix is the formal
        covariance matrix defined as

        .. math::
        
             {\mathbf C} = {\mathbf T} \times {\mathbf \Sigma} \times
             {\mathbf T}^{\rm T},

        where :math:`{\mathbf \Sigma}` is the covariance matrix for the
        'RSS' spectra for the specified wavelength channel.  For 'CUBE'
        files, the function will attempt to use the 'RSS' counterpart of
        the file to produce the transfer matrix, :math:`{\mathbf T}`.
        In this case, the input parameters *must* be the defaults.

        .. note::
        
            The current DRP does not produce spectral covariance
            matrices for the 'RSS' spectra.  Here, it is assumed that
            the spectral covariance matrix is zero everywhere except
            along the diagonal, which contains the inverse of the values
            in the 'IVAR' extension.

        2. If *sigma_rho* is not None, the returned matrix is an
        approximation of the covariance matrix determined by
        *sigma_rho*.  See the description of :attr:`cov_rho` and
        :attr:`sigma_rho`.  For this to work, the variance matrix of the
        reconstructed image must have been already calculated, meaning
        that this approach is most appropriately used with the 'CUBE'
        files.  For 'RSS' files, the covariance matrix is determined in
        this case by calling :func:`covariance_matrix` on the 'CUBE'
        file, which will raise an exception if the parameters defining
        the dimensions of the reconstructed image and kernel are not the
        defaults.

        In general, use of *sigma_rho* should **not** be used with 'RSS'
        files and the use of the formal calculation should **not** be
        used with the 'CUBE' files for efficiency.

        Args:
            channel (int): Index of the spectral channel for which
                to calculate the transfer matrix.
            pixelscale (float): (**Optional**) Desired pixel scale in
                arcsec
            recenter (bool): (**Optional**) Flag to recenter the
                coordinate system
            width_buffer (int): (**Optional**) Number of pixels to use
                as buffer for the image reconstruction
            rlim (float): (**Optional**) The limiting radius of the
                image reconstruction kernel in arcseconds.
            sigma (float): (**Optional**) The sigma of the image
                reconstruction kernel in arcseconds.
            sigma_rho (float): (**Optional**) The sigma of the Gaussian
                function used to approximate the trend of the
                correlation coefficient with pixel separation.
            csr (bool): (**Optional**) Instead of reaturning a
                :class:`mangadap.util.covariance.Covariance` object,
                return the covariance matrix as a
                `scipy.sparse.csr_matrix`_ object.  Primarily used by
                :func:`covariance_cube` for collating the covariance
                matrix of each wavelength channel before combining them
                into a single
                :class:`mangadap.util.covariance.Covariance` object
            quiet (bool): (**Optional**) Suppress terminal output

        Returns:
            :class:`mangadap.util.covariance.Covariance` or
            `scipy.sparse.csr_matrix`_: The covariance matrix for the
            designated wavelength channel.  The return type depends on
            `csr`.

        Raises:
            ValueError: Raised for 'CUBE' files when the input
                parameters are not the same as the defaults.

        .. todo::
            Need to make sure that the correct masks are being used for
            the RSS files.  Should be 3DREJECT but nothing else?
        
        """

        # Set the default values for the input
        pixelscale, recenter, width_buffer, rlim, sigma = \
            self._init_regrid_pars(pixelscale, recenter, width_buffer, rlim, sigma)

        return self._formal_covariance_matrix(channel, pixelscale, recenter, width_buffer,
                                               rlim, sigma, csr=csr, quiet=quiet) \
               if sigma_rho is None else \
               self._approximate_covariance_matrix(channel, pixelscale, recenter, width_buffer,
                                                    rlim, sigma, sigma_rho, csr=csr, quiet=quiet)


    def covariance_cube(self, channels=None, pixelscale=None, recenter=None, width_buffer=None, 
                        rlim=None, sigma=None, sigma_rho=None, csr=False, quiet=False):
        """
        Return the covariance matrices for all wavelength channels; see
        :func:`covariance_matrix`.
        
        Args:
            pixelscale (float): (**Optional**) Desired pixel scale in
                arcsec
            recenter (bool): (**Optional**) Flag to recenter the
                coordinate system
            width_buffer (int): (**Optional**) Number of pixels to use
                as buffer for the image reconstruction
            rlim (float): (**Optional**) The limiting radius of the
                image reconstruction kernel in arcseconds.
            sigma (float): (**Optional**) The sigma of the image
                reconstruction kernel in arcseconds.
            sigma_rho (float): (**Optional**) The sigma of the Gaussian
                function used to approximate the trend of the
                correlation coefficient with pixel separation.
            csr (bool): (**Optional**) Instead of reaturning a
                :class:`mangadap.util.covariance.Covariance` object,
                return a numpy.ndarray of the covariance matrices for
                each channel, which are `scipy.sparse.csr_matrix`_
                objects.  Primarily used by :func:`covariance_cube` for
                collating the covariance matrix of each wavelength
                channel before combining them into a single
                :class:`mangadap.util.covariance.Covariance` object
            quiet (bool): (**Optional**) Suppress terminal output

        Returns:
            :class:`mangadap.util.covariance.Covariance` or
            numpy.ndarray: The return type depends on `csr`: if True,
            the returned object is an ndarray of
            `scipy.sparse.csr_matrix`_ types.

        Raises:
            Exception: Raised for 'CUBE' files when the input parameters
                are not the same as the defaults.

        """

        # Set the default values for the input
        pixelscale, recenter, width_buffer, rlim, sigma = \
            self._init_regrid_pars(pixelscale, recenter, width_buffer, rlim, sigma)

        # Allow CUBE output under certain conditions
        if self.mode == 'CUBE':

            # Do not perform the calculation if the parameters are not
            # the default used by the DRP to create the CUBE file.
            if not self._regrid_defaults(pixelscale, recenter, width_buffer, rlim, sigma):
                raise ValueError('Must use default pixel scale, rlim, and sigma to get '
                                + 'covariance matrices for DRP-produced CUBE files.')

            warnings.warn('Attempting to use RSS counter-part for calculation.')
            drpf = DRPFits(self.plate, self.ifudesign, 'RSS', drpver=self.drpver, \
                           redux_path=self.redux_path, directory_path=self.directory_path)
            return drpf.covariance_cube(channels=channels, sigma_rho=sigma_rho, csr=csr,
                                        quiet=quiet)

        self.open_hdu(checksum=self.checksum)

        # The number of wavelength channels
        _channels = numpy.arange(self.nwave) if channels is None \
                        else numpy.atleast_1d(channels)

        nc = len(_channels)
        CovCube = numpy.empty(nc, dtype=sparse.csr.csr_matrix)   # Empty ndarray

        for i in range(nc):
            if not quiet:
                print('Covariance Cube {0}/{1}'.format(i+1,nc), end="\r")
            CovCube[i] = self.covariance_matrix(_channels[i], pixelscale, recenter, width_buffer,
                                                rlim, sigma, sigma_rho, csr=True, quiet=True)

        if not quiet:
            print('Covariance Cube Done                     ')

        # Don't provide input indices if the full cube is calculated
        return Covariance(CovCube, input_indx=_channels) if not csr else CovCube


    def instrumental_dispersion_plane(self, channel, dispersion_factor=None, pixelscale=None,
                                      recenter=None, width_buffer=None, rlim=None, sigma=None,
                                      pre=False, quiet=False):
        r"""
        Return the instrumental dispersion for the reconstructed 'CUBE'
        wavelength plane.
        
        For 'CUBE' files, the input parameters *must* be the same as the
        defaults.  If they are, the function must be able to find and
        access the 'RSS' file to construct the instrumental dispersion
        map because the necessary information is not in the 'CUBE'
        files (before MPL-6).

        .. todo::
            - Need to implement something that will recognize that the
              'DISP' extension exists in the MPL-6 and later DRP CUBE
              files.
        
        For 'RSS' files, the transfer matrix, :math:`{\mathbf T}`, is
        first calculated using :func:`regrid_transfer_matrix`.  The
        transfer matrix is used to construct the 'CUBE' wavelength plane
        image, :math:`{\mathbf I}`, by computing :math:`{\mathbf T}
        \times {\mathbf F} = {\mathbf I}`, where :math:`{\mathbf F}` is
        the vector of the fiber fluxes.  Under the assumption that the
        line-spread-function (LSF) is Gaussian, we determine the
        instrumental dispersion for the data in the wavelength channel
        of the reconstructed CUBE by calculating the second moment of
        the weighted sum of Gaussians of the appropriate dispersion.
        Assuming all the Gaussians have the normalized form:

        .. math::
            
            \mathcal{N}(x|\mu=0,\sigma) = \frac{1}{\sqrt{2\pi}\sigma}
            \exp\left(-\frac{x^2}{2\sigma^2}\right),

        the combined instrumental dispersion becomes

        .. math::
            
            \sigma_{\rm inst,j}^2 = \frac{\sum_i T_{ij}
            \sigma^2_i}{\sum_i T_{ij}},

        where :math:`T_{ij}` are the elements of the transfer matrix.
        In terms of matrix multiplications, this can be written as

        .. math::

            {\mathbf S} = \frac{ {\mathbf T} \times {\mathbf V} }{
            {\mathbf T_c} },

        where :math:`{\mathbf T_c} = {\mathbf T_c} \times {\mathbf 1}`
        is the vector with the sum :math:`\sum_j T_{ij}`,
        :math:`{\mathbf V}` is the instrumental variance for all fibers
        at the designated wavelength plane, and :math:`{\mathbf S}` is
        the variance for all the spaxels in the reconstructed wavelength
        image; the division by :math:`{\mathbf T_c}` is element-wise.

        The returned matrix is the element-wise square-root of
        :math:`{\mathbf S}`, rearranged into a 2D array of size
        :attr:`ny` by :attr:`nx`.

        .. warning::
            Because the internal data structure is transposed with
            respect to the DRP file, the result is transposed on output
            so that it will match the pixel coordinates in the
            reconstructed image pulled directly from a DRP file.  See
            :func:`regrid_wavelength_channel`.
        
        Args:
            channel (int): Index of the spectral channel for which to
                calculate the transfer matrix.
            dispersion_factor (float): (**Optional**) Artificially
                multiply the dispersion measurements by this factor
                before calculating the reconstructed dispersion.
            pixelscale (float): (**Optional**) Desired pixel scale in
                arcsec
            recenter (bool): (**Optional**) Flag to recenter the
                coordinate system
            width_buffer (int): (**Optional**) Number of pixels to use
                as buffer for the image reconstruction
            rlim (float): (**Optional**) The limiting radius of the
                image reconstruction kernel in arcseconds.
            sigma (float): (**Optional**) The sigma of the image
                reconstruction kernel in arcseconds.
            pre (bool): (**Optional**) Read the pre-pixelized version of
                the spectral resolution, instead of the post-pixelized
                version.  This prepends 'PRE' to the extension name.
            quiet (bool): (**Optional**) Suppress terminal output

        Returns:
            numpy.ndarray: The reconstructed image for the specified
            wavelength channel.

        Raises:
            ValueError: Raised for 'CUBE' files when the input parameters
                are not the same as the defaults.
        """

        # Set the default values for the input
        pixelscale, recenter, width_buffer, rlim, sigma = \
            self._init_regrid_pars(pixelscale, recenter, width_buffer, rlim, sigma)

        # Allow CUBE output under certain conditions
        if self.mode == 'CUBE':

            # Do not perform the calculation if the parameters are not
            # the default used by the DRP to create the CUBE file.
            if not self._regrid_defaults(pixelscale, recenter, width_buffer, rlim, sigma):
                raise ValueError('Must use default pixel scale, rlim, and sigma to get '
                                 + 'wavelength-channel image for DRP-produced CUBE files.')

            drpf = DRPFits(self.plate, self.ifudesign, 'RSS', drpver=self.drpver, \
                           redux_path=self.redux_path, directory_path=self.directory_path)
            return drpf.instrumental_dispersion_plane(channel, dispersion_factor=dispersion_factor,
                                                      quiet=quiet)

        # Set the transfer matrix (set to self.regrid_T; don't need to
        # keep the returned matrix)
        self.regrid_transfer_matrix(channel, pixelscale, recenter, width_buffer, rlim, sigma, quiet)

        # Get the dispersion factor
        _df = 1.0 if dispersion_factor is None else dispersion_factor

        # Return the regridded data with the proper shape (nx by ny)
        Tc = self.regrid_T.sum(axis=1).flatten()
        Tc[numpy.invert(Tc>0)] = 1.0                # Control for zeros
#        ext = 'PREDISP' if pre else 'DISP'
        ext = 'LSFPRE' if pre else 'LSFPOST'
        return numpy.sqrt( self.regrid_T.dot(numpy.square(_df*self.hdu[ext].data[:,channel]))
                                                / Tc ).reshape(self.nx, self.ny)

    def pointing_offset(self):
        """
        Return the offsets in RA and DEC between the pointing
        coordinates (IFURA, IFUDEC) and the designated object center
        coordinates (OBJRA, OBJDEC), drawn from the primary header of
        the DRP fits file.

        Returns:
            float: Sky-right arcsecond offsets in RA and DEC.
        """
        return ((self.hdu[0].header['IFURA'] - self.hdu[0].header['OBJRA']) \
                    * numpy.cos(numpy.radians(self.hdu[0].header['OBJDEC'])) * 3600.), \
               ((self.hdu[0].header['IFUDEC'] - self.hdu[0].header['OBJDEC']) * 3600.)

    def mean_sky_coordinates(self, waverange=None, response_func=None, per_pixel=True, offset=True,
                             flag=None, fluxwgt=False):
        r"""
        Compute the mean sky coordinates for each spectrum.
        
        For CUBE files, this just returns the spaxel coordinates in
        arcseconds relative to the object center, where the object
        center is :math:`(\alpha_0,\delta_0) = ({\rm OBJRA},{\rm
        OBJDEC)` as provided in the primary header, and

        .. math::
            
            x &= (\alpha - \alpha_0) \cos \delta_0 \\
            y &= (\delta - \delta_0)

        The coordinate grid, :math:`(\alpha, \delta)` is based on the
        WCS coordinates in the header as returned by :func:`world_mesh`.

        For RSS files, this returns either the unweighted or
        flux-weighted XPOS and YPOS values.  For observations where the
        pointing center is different from the object center, the
        returned coordinates are relative to the object center if
        offset=True (default) and relative to the pointing center if
        offset=False.

        .. warning::
            
            Flux-weighting the coordinates can produce spurious results
            in low-flux regimes.

        Args:
            waverange (array-like): (**Optional**) Two-element array
                with the first and last wavelength to include in the
                computation.  Default is to use the full wavelength
                range.
            response_func (array-like): (**Optional**) A two-column
                array with the wavelength and transmission of a
                broad-band response function to use for the calculation.
            per_pixel (bool): (**Optional**) When providing a response
                function, continue to calculate the statistics per
                pixel, as opposed to per angstrom.  Default is to
            offset (bool) : Offset the coordinates to the object
                coordinates.
            flag (str or list): (**Optional**) (List of) Flag names that
                are considered when deciding if a pixel should be
                masked.  The names *must* be a valid bit name as defined
                by :attr:`bitmask` (see :class:`DRPFitsBitMask`).
            fluxwgt (bool) : (**Optional**) Flag to weight by the flux
                when determining the mean coordinates for the RSS spectra.

        Returns:
            numpy.ndarray : A 1D array with the coordinates of each
            spectrum.  The index of the vector matches the index in
            :attr:`spatial_index` in case you want to recreate the a map
            for the CUBE files.

        """
        if self.mode == 'CUBE':
            x, y = self.world_mesh()
            x = (x-self.hdu[0].header['OBJRA']) \
                    * numpy.cos(numpy.radians(self.hdu[0].header['OBJDEC']))
            y = (y-self.hdu[0].header['OBJDEC'])
            return x.ravel()*3600., y.ravel()*3600.

        # The negative here is because the XPOS extension has negative
        # offsets going toward increasing RA.
        xoff, yoff = self.pointing_offset() if offset else 0.0, 0.0
        xpos = xoff - self.copy_to_masked_array(ext='XPOS', waverange=waverange, flag=flag)
        ypos = yoff + self.copy_to_masked_array(ext='YPOS', waverange=waverange, flag=flag)
        if fluxwgt:
            flux = self.copy_to_masked_array(ext='FLUX', waverange=waverange, flag=flag)

        # Set the response function
        wave = self['WAVE'].data
        dw = numpy.ones(self.nwave, dtype=float) if per_pixel else \
                    spectral_coordinate_step(wave, log=True)*numpy.log(10.)*wave
        _response_func = self._interpolated_response_function(response_func)

        # Get the normalization and return the flux- or un-weighted coordinates
        if fluxwgt:
            norm = numpy.ma.sum(flux*_response_func[None,:]*dw[None,:], axis=1)
            return numpy.ma.sum(flux*xpos*_response_func[None,:]*dw[None,:],axis=1)/norm, \
                    numpy.ma.sum(flux*ypos*_response_func[None,:]*dw[None,:],axis=1)/norm

        norm = numpy.sum(numpy.invert(numpy.ma.getmaskarray(xpos))*_response_func[None,:]
                                *dw[None,:], axis=1)
        return numpy.ma.sum(xpos*_response_func[None,:]*dw[None,:],axis=1)/norm, \
                    numpy.ma.sum(ypos*_response_func[None,:]*dw[None,:],axis=1)/norm


    def flux_stats(self, waverange=None, response_func=None, per_pixel=True, flag=None,
                   covar=False, correlation=False, covar_wave=None):
        r"""
        Compute the mean flux, propagated error in the mean flux, and
        mean S/N over the specified wavelength range; if the wavelength
        range is not specified, the quantities are calculated over the
        full spectral range.

        If an RSS file, no covariance is calculated.

        If a CUBE file and covar is True, the code will calculate the
        covariance matrix at the specified wavelength (see
        :func:`covariance_matrix`).  If covar_wave is not provided, the
        covariance is calculated at the center wavelength of the
        provided wavelength range.

        Args:
            waverange (array-like): (**Optional**) Starting and ending
                wavelength over which to calculate the statistics.
                Default is to use the full wavelength range.
            response_func (array-like): (**Optional**) A two-column
                array with the wavelength and transmission of a
                broad-band response function to use for the calculation.
            per_pixel (bool): (**Optional**) When providing a response
                function, continue to calculate the statistics per
                pixel, as opposed to per angstrom.  Default is to
                compute the statistics on a per pixel basis.
            flag (str or list): (**Optional**) (List of) Flag names that
                are considered when deciding if a pixel should be
                masked.  The names *must* be a valid bit name as defined
                by :attr:`bitmask` (see :class:`DRPFitsBitMask`).
            covar (bool): (**Optional**) Flag to calculate covariance
                matrix.
            correlation (bool): (**Optional**) Flag to convert the
                covariance matrix to a correlation matrix on ouput.
            covar_wave (double): (**Optional**) Wavelength to use for
                the covariance calculation.

        Returns:
            numpy.ndarray: Four objects are returned: the mean flux, the
            propagated variance in the mean flux, the mean S/N, and the
            covariance/correlation matrix for a single wavelength
            channel.  If the object is an RSS file or no covariance
            calculation is requested, the last returned object is None.
    
        Raises:
            ValueError: Raised of a provided wavelength range object
                does not have two elements.

        """
        if waverange is not None and len(waverange) != 2:
            raise ValueError('Provided wavelength range must be a two-element vector.')
        if response_func is not None:
            if len(response_func.shape) != 2:
                raise ValueError('Response function object must be two dimensional.')
            if response_func.shape[1] != 2:
                raise ValueError('Response function object must only have two columns.')

        # Grab the masked arrays
        wave = self['WAVE'].data
        flux = self.copy_to_masked_array(waverange=waverange, flag=flag)
        ivar = self.copy_to_masked_array(ext='IVAR', waverange=waverange, flag=flag)
        snr = flux*numpy.ma.sqrt(ivar)

        # Set the response function
        dw = numpy.ones(self.nwave, dtype=float) if per_pixel else \
                    spectral_coordinate_step(wave, log=True)*numpy.log(10.)*wave
        _response_func = self._interpolated_response_function(response_func)

#        print(flux.shape)
#        n = numpy.ma.sum(numpy.invert(numpy.ma.getmaskarray(flux)), axis=1)
#        print(n)
#        i = numpy.argsort(n)[-1]
#        print(i)
#
#        pyplot.plot(self['WAVE'].data, dw)
#        pyplot.plot(self['WAVE'].data, flux[i,:])
#        pyplot.plot(self['WAVE'].data, ivar[i,:])
#        pyplot.plot(self['WAVE'].data, snr[i,:])
#        pyplot.plot(self['WAVE'].data, _response_func)
#        pyplot.show()
#        exit()

        # Get the moments
        response_integral = numpy.sum(numpy.invert(numpy.ma.getmaskarray(flux))
                                        *(_response_func*dw)[None,:], axis=1)
        signal = numpy.ma.divide(numpy.ma.sum(flux*(_response_func*dw)[None,:], axis=1),
                                 response_integral)
        variance = numpy.ma.divide(numpy.ma.sum(numpy.ma.power(ivar, -1.) \
                                    * (_response_func*dw)[None,:], axis=1), response_integral)
        snr = numpy.ma.divide(numpy.ma.sum(snr*(_response_func*dw)[None,:], axis=1),
                              response_integral)

#        pyplot.imshow(response_integral.reshape(self.spatial_shape).T, origin='lower',
#                      interpolation='nearest')
#        pyplot.show()
#
#        pyplot.imshow(signal.reshape(self.spatial_shape).T, origin='lower', interpolation='nearest')
#        pyplot.show()
#
#        pyplot.imshow(variance.reshape(self.spatial_shape).T, origin='lower', interpolation='nearest')
#        pyplot.show()
#
#        pyplot.imshow(snr.reshape(self.spatial_shape).T, origin='lower', interpolation='nearest')
#        pyplot.show()

        # No computation of spatial covariance needed or requested
        if self.mode == 'RSS' or not covar:
            return signal, variance, snr, None

        # Test if the RSS file exists
        if self.mode == 'CUBE' and covar:
            rss = DRPFits(self.plate, self.ifudesign, 'RSS', drpver=self.drpver,
                          redux_path=self.redux_path, directory_path=self.directory_path,
                          read=False)
            if not os.path.isfile(rss.file_path()):
                warnings.warn('RSS counterpart not available.  Cannot determine covariance matrix!')
                return signal, variance, snr, None

        # Only calculate the covariance at the central, or input, wavelength
        _covar_wave = covar_wave if covar_wave is not None \
                        else self._covariance_wavelength(waverange=waverange,
                                                         response_func=response_func,
                                                         per_pixel=per_pixel, flag=flag)

        channel = numpy.argsort( numpy.absolute(wave - _covar_wave) )[0]

#        t = self.bitmask.flagged(self.hdu['MASK'].data[:,:,channel],
#                                 flag=['NOCOV', 'LOWCOV', 'DEADFIBER'])
#        tt = self.bitmask.flagged(self.hdu['MASK'].data[:,:,channel],
#                                  flag=['FORESTAR', 'DONOTUSE'])
#        print(numpy.sum(t))
#        print(numpy.sum(tt))
#        print(numpy.sum(t & ~tt))

        # Return the nominal stats and the covariance matrix
        C = self.covariance_matrix(channel)
        if correlation:
            C.to_correlation()

        return signal, variance, snr, C


    def binned_on_sky_area(self, bin_indx, x=None, y=None):
        r"""
        Compute the on-sky area of a set of binned spectra.  For CUBE
        files, this is just the number of spaxels in the bin times the
        spaxel area (as given by :attr:`pixelscale`).

        For RSS files, this will try to calculate the overlapping area
        of the fibers using the `shapely`_ python package:

            - The fibers "beams" are all renormalized to have an area of
              pi arcsec^2 by the DRP, so it's radius 1 arcsec
            - This function will provide the *total* area, not the
              integration-weighted effective area.

        Args:
            bin_indx (array-like): A vector with size :math:`N_{\rm
                spec}` the gives which spaxels or fibers were included
                in each bin.  Valid bins have indices of :math:`\geq 0`.
            x (array-like): (**Optional**) On-sky :math:`x` coordinate.
                Default is to calculate :math:`x` and :math:`y` using
                :func:`mean_sky_coordinates` with no arguments.
            y (array-like): (**Optional**) On-sky :math:`y` coordinate.
                Default is to calculate :math:`x` and :math:`y` using
                :func:`mean_sky_coordinates` with no arguments.

        Returns:
            numpy.ndarray : The on-sky area of each bin.

        """
        unique_bins, bin_count = numpy.unique(bin_indx, return_counts=True)
        indx = numpy.invert(unique_bins < 0)
        nbin = bin_count[indx]
        if self.mode == 'CUBE':
            if self.pixelscale is None:
                self.pixelscale = defaults.cube_pixelscale()
            return (nbin*numpy.square(self.pixelscale)).astype(float)

        try:
            if x is None or y is None:
                x, y = self.mean_sky_coordinates()
            good_bins = unique_bins[indx]

            area = numpy.empty(nbin.shape, dtype=numpy.float)
            for i, b in enumerate(good_bins):
                _x = x[bin_indx == b]
                _y = y[bin_indx == b]

                area[i] = cascaded_union([ Point(xx,yy).buffer(1.0,64) \
                                                        for xx, yy in zip(_x,_y) ]).area
            return area
        except:
            warnings.warn('Could not use \'shapely\' package to compute overlapping fiber area.' \
                          'Return the total fiber area.', ImportWarning)
            return (nbin*numpy.pi).astype(float)

    @property
    def can_compute_covariance(self):
        if self.mode == 'RSS':
            return True
       
        # Try to find the RSS file
        rss = DRPFits(self.plate, self.ifudesign, 'RSS', drpver=self.drpver,
                          redux_path=self.redux_path, directory_path=self.directory_path,
                          read=False)
        if not os.path.isfile(rss.file_path()):
            return False
        return True
        

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



        

    
        




