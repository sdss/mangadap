"""

Class that reads and prepares template libraries for use in fitting the
stellar-continuum of a spectrum.  See
:func:`mangadap.util.defaults.default_template_libraries` for the list
of default template libraries.

The "raw" template libraries are expected to consist of 1D fits files
that are linearly sampled in wavelength.  The fits files must have
CRVAL1, CRPIX1, and CDELT1 keywords used to define the wavelength
coordinates of each pixel::
 
    wave = (numpy.arange(1,npixels+1) - CRPIX1) * CDELT1 + CRVAL1
 
The reference frame of the template wavelengths must also be defined as
either vacuum or air.  It is expected that the spectra to be fit are
calibrated to vacuum wavelengths.  When preparing the template spectra
for analysis, this class will use
:func:`mangadap.util.idlutils.airtovac` to convert the template
wavelengths to vacuum.  Finally, one can specify that the template
library is only valid within a certain wavelength range and above a
certian flux limit; see
:class:`mangadap.util.par.TemplateLibraryParSet`.

Preparation of the template library for use in stellar-continuum fitting
consists of a number of steps.  See
:func:`TemplateLibrary.process_template_library`.

A template library that has been prepared for analysis is written to
disk for later recovery.

*Source location*:
    $MANGADAP_DIR/python/mangadap/proc/TemplateLibrary.py

*Python2/3 compliance*::

    from __future__ import division
    from __future__ import print_function
    from __future__ import absolute_import
    from __future__ import unicode_literals
    
    import sys
    if sys.version > '3':
        long = int

*Imports*::

    import glob
    import os.path
    from os import environ, remove
    from scipy import sparse
    from astropy.io import fits
    from astropy import constants
    import time
    import numpy
    from scipy.interpolate import InterpolatedUnivariateSpline
    from mangadap.util.defaults import default_dap_source, default_drp_version, default_dap_version
    from mangadap.util.defaults import default_analysis_path, default_dap_directory_path
    from mangadap.util.defaults import default_template_libraries, default_template_library_file
    from mangadap.util.idlutils import airtovac
    from mangadap.util.fileio import readfits_1dspec
    from mangadap.util.bitmasks import TemplateLibraryBitMask, HDUList_mask_wavelengths
    from mangadap.util.instrument import log_rebin, log_rebin_pix
    from mangadap.util.instrument import match_spectral_resolution, spectrum_velocity_scale

*Class usage examples*:

    Assuming you have the default directory structure setup, you can do::

        # Imports
        from mangadap.drpfile import drpfile
        from mangadap.TemplateLibrary import TemplateLibrary
        from matplotlib import pyplot

        # Define the DRP file
        drpf = drpfile(7495, 12703, 'CUBE')

        # Build the template library
        tpl_lib = TemplateLibrary('M11-MILES', drpf=drpf, directory_path='.')
        # Writes: ./manga-7495-12703-LOGCUBE_M11-MILES.fits

        # Plot one of the spectra
        pyplot.plot(tpl_lib.hdu['WAVE'].data, tpl_lib.hdu['FLUX'].data[0,:])
        pyplot.show()

    As part of the instantiation of the :class:`TemplateLibrary` object
    in the above call, the template library is prepared for use in
    fitting the specified DRP file.  If the required processing has
    already been done, the instantiation of the :class:`TemplateLibrary`
    object simply reads the existing file.  If you do not have the
    default directory structure setup, you'll need to define the paths
    to, e.g., the DRP file; see :class:`mangadap.dapfile.dapfile`.

    If you do not want to process the template (match the spectral
    resolution and sampling to that of the DRP data), you can force the
    :class:`TemplateLibrary` object to only provide the "raw" spectra::

        # Imports
        from mangadap.TemplateLibrary import TemplateLibrary
        from matplotlib import pyplot

        # Read the raw template library
        tpl_lib = TemplateLibrary('M11-MILES', process=False)
        # Nothing should be written to disk

        # Plot one of the spectra
        pyplot.plot(tpl_lib.hdu['WAVE'].data[0,:], tpl_lib.hdu['FLUX'].data[0,:])
        pyplot.show()

    Note that in the previous example, the wavelength data was already
    one dimensional, whereas for the raw library, the wavelength vector
    can be spectrum-dependent.

    In the above examples, the user has not provided a list of template
    libraries, meaning that the default set available to the DAP is
    used.  The default set is defined in
    :func:`mangadap.util.defaults.default_template_libraries`.  If you
    want to use your own template library, you have to define its
    parameters using :class:`mangadap.util.par.TemplateLibraryParSet`.
    Currently, the template library spectra are expected to be 1D fits
    files with WCS header keywords defining the wavelength solution; see
    above.  Using an existing DAP library as an example::

        # Imports
        from mangadap.TemplateLibrary import TemplateLibrary
        from mangadap.util.par import TemplateLibraryParSet
        from mangadap.util.defaults import default_dap_source
        from matplotlib import pyplot

        # Define the search string for the library
        search_str = default_dap_source()+'/external/templates/miles/*.fits'

        # Define the template library parameters
        tpl_par = TemplateLibraryParSet(key='MILES',    # Unique keyword for the library
                                file_search=search_str, # Search string
                                fwhm=2.50,              # FWHM of resolution element (assumed const)
                                in_vacuum=False,        # Wavelength in vacuum?
                                wave_limit=numpy.array([ 3575., 7400. ]),   # Range of valid lambda
                                lower_flux_limit=0.0)   # Lower limit for valid flux

        # Read the raw template library
        tpl_lib = TemplateLibrary('MILES', tpllib_list=tpl_par, process=False)
        # Nothing should be written to disk

        # Plot one of the spectra
        pyplot.plot(tpl_lib.hdu['WAVE'].data[0,:], tpl_lib.hdu['FLUX'].data[0,:])
        pyplot.show()

    Note that the keyword you use must be provided both to the parameter
    set and when instantiating the :class:`TemplateLibrary` object.  In
    the example above, I have not processed the library, but you can by
    following a similar approach to the first example.

*Revision history*:
    | **23 Apr 2015**: Implementation begun by K. Westfall (KBW)
    | **26 May 2015**: (KBW) Added some Sphinx documentation.
    | **17 Jun 2015**: (KBW) Added flexibility in definition of template
        libraries from which to choose using the new
        :class:`mangadap.util.par.TemplateLibraryParSet` class.
    | **23 Jun 2015**: (KBW) Allow user to provided non-DRP input
        spectra, meaning they need to provide the velocity scale and the
        wavelength and spectral resolution vectors.  They must also
        directly set the name of the output processed file.

.. _astropy.io.fits.hdu.hdulist.HDUList: http://docs.astropy.org/en/v1.0.2/io/fits/api/hdulists.html
.. _glob.glob: https://docs.python.org/3.4/library/glob.html

"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import sys
if sys.version > '3':
    long = int

import glob
import os.path
from os import environ, remove
from scipy import sparse
from astropy.io import fits
from astropy import constants
import time
import numpy

from scipy.interpolate import InterpolatedUnivariateSpline

from mangadap.util.defaults import default_dap_source, default_drp_version, default_dap_version
from mangadap.util.defaults import default_analysis_path, default_dap_directory_path
from mangadap.util.defaults import default_template_libraries, default_template_library_file
from mangadap.util.idlutils import airtovac
from mangadap.util.fileio import readfits_1dspec
from mangadap.util.bitmasks import TemplateLibraryBitMask, HDUList_mask_wavelengths
from mangadap.util.instrument import log_rebin, log_rebin_pix
from mangadap.util.instrument import match_spectral_resolution, spectrum_velocity_scale
from mangadap.util.par import TemplateLibraryParSet

from matplotlib import pyplot

__author__ = 'Kyle B. Westfall'

class TemplateLibrary:
    """
    Object used to read, store, and prepare template libraries for use
    in analyzing object spectra.

    The default list of available libraries provided by the MaNGA DAP
    defined in
    :func:`mangadap.util.defaults.default_template_libraries`.  The user
    can provide their own library for use with this class provided they
    are contained in 1D fits spectra, sampled linearly in wavelength
    with the wavelength coordinates available via the WCS keywords
    (CRPIX1, CRVAL1, CDELT1), and they have an appropriately defined
    spectral resolution (FWHM in angstroms that is constant as a
    function of wavelength).  See
    :class:`mangadap.util.par.TemplateLibraryParSet` and
    :func:`_build_raw`.

    The class is optimized for use in analyzing MaNGA DRP files;
    however, one can provide the necessary information so that the class
    can be used with a non-DRP spectrum.  In the latter case, the user
    must supply the velocity scale of the pixel for the logarithmically
    resampled template library, and a
    :class:`mangadap.util.instrument.spectral_resolution` object the
    defines the instrumental resolution of the spectrum/spectra to be
    analyzed.


    .. todo::
        - below is out of date.

    On initialization, if the DRP file object is not provided (is None),
    the default behavior is to read the raw template library if
    read=True.  If the DRP file is provided, the routine will check for
    the resolution matched fits file; if it doesn't exist and read is
    True, it will prepare the template library for use in analyzing the
    DRP file and write the prepared library file.  If force=True, the
    preparation and writing of the template library will be done even if
    the library already exists.

    Args:
        library_key (str): Keyword selecting the library to use.
        tpllib_list (list): (Optional) List of
            :class:`mangadap.util.par.TemplateLibraryParSet` objects
            that define the parameters required to read and interpret a
            template library.
        dapsrc (str): (Optional) Root path to the DAP source directory.
            If not provided, the default is defined by
            :func:`mangadap.util.defaults.default_dap_source`.
        drpf (:class:`mangadap.drpfile.drpfile`): (Optional) DRP file
            (object) with which the template library is associated for
            analysis
        velscale (float): (Optional) The velocity scale of the spectrum
            that will be analyzed with the library; this is used in
            place of the attributes in any provided DRP file object.
        sres (:class:`mangadap.util.instrument.spectral_resolution`):
            (Optional) The object is used simply to access the spectral
            resolution and associated wavelength coordinate vector
            needed when matching the spectral resolution of the template
            library; this is used in place of the attributes in any
            provided DRP file object.
        velocity_offset (float): (Optional) Velocity offset to use when
            matching the spectral resolution between the template
            library and the galaxy spectra.
        dapver (str): (Optional) DAP version, which is used to define
            the default DAP analysis path.  Default is defined by
            :func:`mangadap.util.defaults.default_dap_version`
        analysis_path (str): (Optional) The path to the top level
            directory containing the DAP output files for a given DRP
            and DAP version.  Default is defined by
            :func:`mangadap.util.defaults.default_analysis_path`.
        directory_path (str): (Optional) The exact path to the processed
            template library file.  Default is defined by
            :func:`mangadap.util.defaults.default_dap_directory_path`.
        processed_file (str): (Optional) The name of the file containing
            the prepared template library output file.  The file should
            be found at :attr:`directory_path`/:attr:`processed_file`.
            Default is defined by
            :func:`mangadap.util.defaults.default_template_library_file`.
        process (bool): (Optional) If :attr:`drpf` is defined and the
            prepared template library does not exist, this will process
            the template library in preparation for use in fitting the
            provided DRP file.
        force (bool): (Optional) If :attr:`drpf` is define and *process*
            is True, this will force the template library to be
            processed, even if the prepared template library already
            exists.

    Attributes:
        version (str): Version number
        tplbm (BitMask): A BitMask object used to toggle mask values;
            see :func:`mangadap.util.bitmasks.TemplateLibraryBitMask`.
        library (str): Keyword of the selected library to use.
        libparset (:class:`mangadap.util.par.TemplateLibraryParSet`):
            Parameter set required to read and prepare the library.
        file_list (list): The list of files found using `glob.glob`_ and
            :attr:`file_search`.
        ntpl (int): Number of template spectra in the library
        drpf (:class:`mangadap.drpfile.drpfile`): DRP file (object) with
            which the template library is associated for analysis
        velscale (float): The velocity scale of the spectrum that will
            be analyzed with the library; this is used in place of the
            attributes in any provided DRP file object.
        sres (:class:`mangadap.util.instrument.spectral_resolution`):
            The object is used simply to access the spectral resolution
            and associated wavelength coordinate vector needed when
            matching the spectral resolution of the template library;
            this is used in place of the attributes in any provided DRP
            file object.
        velocity_offset (float): Velocity offset to use when matching
            the spectral resolution between the template library and the
            galaxy spectra.
        directory_path (str): The exact path to the processed template
            library file.  Default is defined by
            :func:`mangadap.util.defaults.default_dap_directory_path`.
        processed_file (str): The name of the file containing (to
            contain) the prepared template library output file.  The
            file should be found at
            :attr:`directory_path`/:attr:`processed_file`.
        processed (bool): Flag that the template library has been
            prepared for use in the DAP.
        hdu (`astropy.io.fits.hdu.hdulist.HDUList`_): HDUList read from
            the DAP file

    .. todo::

        - Add USER library key and needed functionality.
        - Allow for input of sres vector and sampling to process the
          file instead of the DRP file

    """
    def __init__(self, library_key, tpllib_list=None, dapsrc=None, drpf=None, velscale=None,
                 sres=None, velocity_offset=0.0, dapver=None, analysis_path=None,
                 directory_path=None, processed_file=None, read=True, process=True, force=False):

        self.version = '2.1'

        # Define the TemplateLibraryBitMask object
        self.tplbm = TemplateLibraryBitMask()

        # Define the library properties
        self.library = None
        self.libparset = None
        self.file_list = None
        self.ntpl = None

        # Define the properties needed to prepare the library for
        # analysis
        self.drpf = drpf
        self.velscale = None
        self.sres = None
        self.velocity_offset = None

        # Define the processed file and flag, and the HDUList used to
        # keep the data
        self.directory_path = None
        self.processed_file = None
        self.processed = False
        self.hdu = None

        self._define_library(library_key, tpllib_list=tpllib_list, dapsrc=dapsrc)

        # Do not read the library
        if not read:
            print('Do not read.')
            return
       
        # Do not process the library
        if not process:
            print('Reading raw library without processing.')
            self._read_raw()
            return

        # Read and process the library
        self.process_template_library(velscale=velscale, sres=sres, velocity_offset=velocity_offset,
                                      dapver=dapver, analysis_path=analysis_path,
                                      directory_path=directory_path, processed_file=processed_file,
                                      force=force)


    def __del__(self):
        """
        Deconstruct the template library object by ensuring that the
        fits file is properly closed.
        """
        if self.hdu is None:
            return
        self.hdu.close()
        self.hdu = None


    def _define_library(self, library_key, tpllib_list=None, dapsrc=None):
        """
        Select the library from the provided list.  Used to set
        :attr:`library` and :attr:`libparset`.

        Args:
            library_key (str): Keyword of the selected library.
                Available libraries are proved by
                :func:`mangadap.util.defaults.default_template_libraries`
            tpllib_list (list): (Optional) List of 
                :class:`mangadap.util.par.TemplateLibraryParSet' objects
                that define the parameters required to read and
                interpret a template library.
            dapsrc (str): (Optional) Root path to the DAP source
                directory.  If not provided, the default is defined by
                :func:`mangadap.util.defaults.default_dap_source`.

        Raises:
            KeyError: Raised if the selected keyword is not among the
                provided list or if the provided list has more than one
                identical keyword.
            TypeError: Raised if the input *tpllib_list* object is not a
                list or a
                :class:`mangadap.util.par.TemplateLibraryParSet`.

        """
        # Get the default libraries if no list provided
        if tpllib_list is None:
            tpllib_list = default_template_libraries(dapsrc=dapsrc)

        # Make sure the input tpllib_list has the right type, and force
        # it to be a list
        if type(tpllib_list) != list and not isinstance(tpllib_list, TemplateLibraryParSet):
            raise TypeError('Input tpllib_list must be a list or a TemplateLibraryParSet!')
        if type(tpllib_list) != list:
            tpllib_list = [tpllib_list]

        # Find the selected library via its keyword
        selected_lib = [ l['key'] == library_key for l in tpllib_list ]
        if numpy.sum(selected_lib) == 0:
            raise KeyError('{0} is not a valid template library!'.format(library_key))
        if numpy.sum(selected_lib) > 1:
            raise KeyError('More than one library has the same key!')

        # Save the parameters for this template library
        indx = numpy.where(selected_lib)[0][0]
        self.library = library_key
        self.libparset = tpllib_list[indx]

#        self.file_search = tpllib_list[indx]['file_search']
#        self.fwhm = tpllib_list[indx]['fwhm']
#        self.invac = tpllib_list[indx]['in_vacuum']
#        self.wave_limit = tpllib_list[indx]['wave_limit']
#        self.flux_limit = tpllib_list[indx]['lower_flux_limit']


    def _can_set_paths(self, directory_path, processed_file, quiet=False):
        # Check that the directory_path can be set
        if self.drpf is None and directory_path is None:
            if not quiet:
                raise Exception('Cannot define the default directory path without a DRP file.')
            return False

        # Check that the directory_path can be set
        if self.drpf is None and processed_file is None:
            if not quiet:
                raise Exception('Cannot define the default output file name without a DRP file.')
            return False

        return True


    def _set_paths(self, dapver, analysis_path, directory_path, processed_file):
        """
        Set the I/O path to the processed template library.  Used to set
        :attr:`directory_path` and :attr:`processed_file`.  If not
        provided, the defaults are set using, respectively,
        :func:`mangadap.util.defaults.default_dap_directory_path` and
        :func:`mangadap.util.defaults.default_template_library_file`.

        .. warning::

            :attr:`drpf` must have been set beforehand!

        Args:
            dapver (str): DAP version.
            analysis_path (str): The path to the top-level directory
                containing the DAP output files for a given DRP and DAP
                version.
            directory_path (str): The exact path to the DAP template
                library file.  See :attr:`directory_path`.
            processed_file (str): The name of the file with the prepared
                template library.  See
                :attr:`processed_file`.

        """
        # Use this to raise the necessary exceptions
        self._can_set_paths(directory_path, processed_file)

        # Set the output directory path
        if directory_path is None:
            dapver = default_dap_version() if dapver is None else str(dapver)
            analysis_path = default_analysis_path(self.drpf.drpver, dapver) \
                            if analysis_path is None else str(analysis_path)
            self.directory_path = default_dap_directory_path(self.drpf.drpver, dapver,
                                                             analysis_path, self.drpf.plate,
                                                             self.drpf.ifudesign)
        else:
            self.directory_path = str(directory_path)

        # Set the output file
        self.processed_file = default_template_library_file(self.drpf.plate, self.drpf.ifudesign, \
                                                            self.drpf.mode, self.library) \
                              if processed_file is None else str(processed_file)


    def _read_raw(self):
        """
        Read the 'raw' versions of the template library; i.e., the
        library before it has been resolution and sampling matched to a
        DRP file.
        """
        self.file_list = glob.glob(self.libparset['file_search'])
        self.ntpl = len(self.file_list)
        print('Found {0} {1} templates'.format(self.ntpl, self.library))

        npix = self._get_nchannels()
        print('Maximum number of wavelength channels: {0}'.format(npix))

        self._build_raw_hdu(npix)


    def _get_nchannels(self):
        """
        Get the maximum number of wavelength channels needed to store
        all the template spectra in a single array.

        Returns:
            int : Maximum number of pixels used by the listed spectra.

        """
        max_npix = 0
        for f in self.file_list:
            if fits.getval(f, 'NAXIS') != 1:
                raise Exception('{0} is not one dimensional!'.format(f))
            npix = fits.getval(f, 'NAXIS1')
            if max_npix < npix:
                max_npix = npix
        return max_npix


    def _build_raw_hdu(self, npix):
        r"""
        Build the "raw" template library arrays.  This simply reads the
        provided list of fits files and puts them into arrays of size
        :math:`N_{\rm tpl} \times N_{\rm pix}`.

        This will *force* reading of the data, even if the :attr:`hdu`
        is already initialized.

        The :attr:`hdu` will contain the appropriate extensions, but it
        is important to note that the wavelength vectors will **not**
        necessarily be the same.  That is, reading of the raw template
        spectra can accommodate spectra that have different wavelength
        coordinates.  Any pixels that have no data are masked using the
        'NO_DATA' bitmask flag; see
        :func:`mangadap.util.bitmasks.TemplateLibraryBitMask`.

        The spectral resolution is set using :attr:`fwhm`, and the
        spectral resolution offset is initialized to zero (see
        :func:`mangadap.util.instrument.GaussianKernelDifference`).

        .. warning::

            Currently no errors are saved because none are expected for
            the template libraries.

        Args:
            npix (int): Number of spectral channels for the output
                arrays

        """
        if self.hdu is not None:
            print('Closing existing HDUList.')
            self.hdu.close()
            self.hdu = None

        print('Attempting to build raw data ...')
        # Add some header cards
        wave = numpy.zeros((self.ntpl, npix), dtype=numpy.float64)
        flux = numpy.zeros((self.ntpl, npix), dtype=numpy.float64)
        mask = numpy.zeros((self.ntpl, npix), dtype=numpy.uint8)
        soff = numpy.zeros(self.ntpl, dtype=numpy.float64)
#        ivar = numpy.zeros((self.ntpl, npix), dtype=numpy.float64)
#        ivar[:] = 1.0

        # Read and save each spectrum and mask the unobserved
        # wavelengths
        for i in range(0,self.ntpl):
            wave_, flux_ = readfits_1dspec(self.file_list[i])
            wave[i,0:wave_.size] = numpy.copy(wave_)
            flux[i,0:wave_.size] = numpy.copy(flux_)
            if wave_.size != npix:
                mask[i,wave_.size:] = self.tplbm.turn_on(mask[i,wave_.size:],'NO_DATA')
            if self.libparset['lower_flux_limit'] is not None:
                indx = numpy.invert( flux_ > self.libparset['lower_flux_limit'] )
                mask[i,indx] = self.tplbm.turn_on(mask[i,indx], 'FLUX_INVALID')
            if self.libparset['wave_limit'][0] is not None:
                indx = wave[i,:].ravel() < self.libparset['wave_limit'][0]
                mask[i,indx] = self.tplbm.turn_on(mask[i,indx], 'WAVE_INVALID')
            if self.libparset['wave_limit'][1] is not None:
                indx = wave[i,:].ravel() > self.libparset['wave_limit'][1]
                mask[i,indx] = self.tplbm.turn_on(mask[i,indx], 'WAVE_INVALID')

        # Set the spectral resolution
        sres = wave/self.libparset['fwhm']

        # (Re)Set the HDUList object
        self._reset_hdu(wave, flux, mask, sres, soff)

        # Add some keywords to the header
        self.hdu[0].header['TPLPROC'] = (0, 'Flag that library has been processed for DAP use')
        print('... done')


    def _reset_hdu(self, wave, flux, mask, sres, soff):
        r"""
        (Re)Set :attr:`hdu` to a new HDUList object using the input
        arrays.  Also sets the header items indicating the version of
        the template reader and the keyword for the library.

        .. warning:
            No checking is done concerning the size of each array!

        Args:
            wave (numpy.ndarray): Array with the wavelengths of each
                pixel.
            flux (numpy.ndarray): Array with the flux in each pixel.
            mask (numpy.ndarray): Bitmask values for each pixel.
            sres (numpy.ndarray): Spectral resolution,
                :math:`R=\lambda/\delta\lambda`, at each pixel.
            soff (numpy.ndarray): The spectral resolution offset for
                each spectrum (see
                :func:`mangadap.util.instrument.GaussianKernelDifference`). 

        """
        if self.hdu is not None:
            self.hdu.close()
            self.hdu = None

        self.hdu = fits.HDUList([fits.PrimaryHDU()])
        self.hdu.append(fits.ImageHDU(wave, name='WAVE'))
        self.hdu.append(fits.ImageHDU(flux, name='FLUX'))
#        self.hdu.append(fits.ImageHDU(ivar, name='IVAR'))
        self.hdu.append(fits.ImageHDU(mask, name='MASK'))
        self.hdu.append(fits.ImageHDU(sres, name='SPECRES'))
        self.hdu.append(fits.ImageHDU(soff, name='SIGOFF'))
        self.hdu[0].header['VDAPTPL'] = (self.version, 'Version of DAP template reader')
        self.hdu[0].header['LIBKEY'] = (self.library, 'Library identifier')


    def _wavelength_range(self, flag=None):
        """
        Return the valid wavelength range for each spectrum based on the
        first and last unmasked pixel; interspersed masked regions are
        not considered.

        Args:
            flag (str or list): (Optional) Flags to consider when
                determining the wavelength range; see
                :func:`mangadap.util.bitmasks.BitMask.flagged`.

        Returns:
            numpy.ndarray : Two-element vector with wavelengths of the
                first and last valid pixels.
        
        """

        indx = numpy.where(numpy.invert(self.tplbm.flagged(self.hdu['MASK'].data, flag=flag)))
        return numpy.array([ numpy.amin(self.hdu['WAVE'].data[indx]),
                             numpy.amax(self.hdu['WAVE'].data[indx])])


    def _rebin_masked(self, i, flag, fullRange, rmsk_lim=0.5):
        """
        Determine the mask value to adopt for a rebinned spectrum by
        rebinning the mask pixels, setting a masked pixel to unity, and
        an unmasked pixel to zero.  After rebinning, any pixel with a
        value of larger than *extrap_lim* should be mask; otherwise it
        is left unmasked.

        Although the function can be used with multiple flags, its
        intended use is to determine which pixels should be masked with
        a specific flag.

        .. todo::
            - Allow the rebinned pixel value distinguishing masked and
              unmasked to be an optional parameter.

        Args:
            i (int): Index of the spectrum to be rebinned.
            flag (str or list): Flags to consider when determining which
                pixels to mask; see
                :func:`mangadap.util.bitmasks.BitMask.flagged`.
            fullRange (numpy.ndarray): Two-element array with the
                wavelength range for the rebinned spectrum.
            rmsk_lim (float): Limit of the rebinned mask value that is
                allowed before considering the pixel as masked.

        Returns:
            tuple : The indices of the pixels in the rebinned spectrum
                that should be masked.

        """
        mask_ex = self.tplbm.flagged(self.hdu['MASK'].data[i,:], flag=flag).astype(numpy.float64)
        mask_ex, wave, dv = log_rebin([self.hdu['WAVE'].data[i,0], self.hdu['WAVE'].data[i,-1]],
                                      mask_ex, velscale=self.velscale, log10=True,
                                      newRange=fullRange, wave_in_ang=True)
        return numpy.where(mask_ex > rmsk_lim)


    def _process_library(self):
        """
        Process the template library for use in analyzing object
        spectra.   See :func:`process_template_library`.

        .. todo::

            - Make the resampling more efficient?
        
        """
        # Convert to vacuum wavelengths
#        pyplot.plot(self.hdu['WAVE'].data[0,:], self.hdu['FLUX'].data[0,:]) 
        if not self.libparset['in_vacuum']:
            self.hdu['WAVE'].data = airtovac(self.hdu['WAVE'].data)
#        pyplot.plot(self.hdu['WAVE'].data[0,:], self.hdu['FLUX'].data[0,:], 'g') 
#        pyplot.show()

        # Calculate the redshift to use to set the DRP observed
        # wavelength to the rest wavelength; important for matching the
        # spectral resolution of the template and galaxy spectra in the
        # same reference frame
        redshift = self.velocity_offset/constants.c.to('km/s').value

        # Get the velocity scale and spectral resolution elements
        if self.velscale is not None and self.sres is not None:
            sres_wave = self.sres.wave()
            sres_val = self.sres.sres()
        else:
            # Make sure the DRP file is open
            if self.drpf.hdu is None:
                print('Opening DRP file ... ')
                self.drpf.open_hdu()
                print('... done.')
            sres_wave = self.drpf.hdu['WAVE'].data
            sres_val = self.drpf.hdu['SPECRES'].data
            self.velscale = spectrum_velocity_scale(self.drpf.hdu['WAVE'].data, log10=True)
            # Close the DRP file?

        # Mask wavelengths where the spectral resolution will have to be
        # extrapolated for cross-matching with the template spectra,
        # accounting for the redshift of the galaxy.  
        print('Masking extrapolation wavelengths ... ')
        wavelim = numpy.array([ sres_wave[0]/(1.+redshift), sres_wave[-1]/(1.+redshift) ])
        self.hdu = HDUList_mask_wavelengths(self.hdu, self.tplbm, 'SPECRES_EXTRAP', wavelim,
                                            invert=True)
        print('... done.')

#        oldwave = numpy.copy(self.hdu['WAVE'].data[0,:]).ravel()
#        oldflux = numpy.copy(self.hdu['FLUX'].data[0,:]).ravel()
#        pyplot.plot(self.hdu['WAVE'].data[0,:], self.hdu['SPECRES'].data[0,:]) 
#        pyplot.plot(self.drpf.hdu['WAVE'].data, self.drpf.hdu['SPECRES'].data, 'r') 

        # Match the resolution of the templates.
        print('Matching spectral resolution ... ')
        self.hdu['FLUX'].data, self.hdu['SPECRES'].data, self.hdu['SIGOFF'].data, res_mask = \
            match_spectral_resolution(self.hdu['WAVE'].data, self.hdu['FLUX'].data,
                                      self.hdu['SPECRES'].data, sres_wave/(1.+redshift), sres_val,
                                      min_sig_pix=0.0, new_log10=True)
        print('... done')

#        pyplot.plot(self.hdu['WAVE'].data[0,:], self.hdu['SPECRES'].data[0,:], 'g') 
#        pyplot.show()

        # Mask any pixels where the template resolution was too low to
        # match to the galaxy resolution
        print('Masking low spectral resolution ... ')
        self.hdu['MASK'].data[res_mask == 1] = \
            self.tplbm.turn_on(self.hdu['MASK'].data[res_mask == 1], 'SPECRES_LOW')
        print('... done')

#        pyplot.plot(oldwave, oldflux)
#        pyplot.plot(self.hdu['WAVE'].data[0,:], self.hdu['FLUX'].data[0,:], 'g')
#        pyplot.show()

        ################################################################
        # Resample the templates to match the object spectra.
        print('Matching sampling ... ')
        # The raw spectra are allowed to have wavelength ranges that
        # differ.  First, determine the wavelength range that encloses
        # all spectra.  Only ignore pixels that were flagged as having
        # no data.
        fullRange = self._wavelength_range(flag='NO_DATA')

        # Get the maximum number of pixels needed to cover this full
        # range for each spectrum, imposing the correct velocity
        # sampling
        npix = self.hdu['FLUX'].data.shape[1]
        dw, m, dlogw, dv = log_rebin_pix([self.hdu['WAVE'].data[0,0], self.hdu['WAVE'].data[0,-1]],
                                         npix, velscale=self.velscale, log10=True,
                                         newRange=fullRange)
        max_m = m
        for i in range(1,self.ntpl):
            dw, m, dlogw, dv = log_rebin_pix([self.hdu['WAVE'].data[i,0],
                                              self.hdu['WAVE'].data[i,-1]], npix,
                                              velscale=self.velscale, log10=True,
                                              newRange=fullRange)
            if max_m < m:
                max_m = m

        print(fullRange, self.velscale, max_m, self.ntpl, npix)

        # Any pixels without data after resampling are given a value
        # that is the minimum flux - 100 so that they can be easily
        # identified afterward.  The minimum flux is:
        min_flux = numpy.amin(self.hdu['FLUX'].data.ravel())
        # the observed pixels are
        observed = numpy.invert(self.tplbm.flagged(self.hdu['MASK'].data, flag='NO_DATA'))
        
        # Now resample the spectra.  First allocate the arrays
        flux = numpy.zeros((self.ntpl, max_m), dtype=numpy.float64)
        sres = numpy.zeros((self.ntpl, max_m), dtype=numpy.float64)
        mask = numpy.zeros((self.ntpl, max_m), dtype=numpy.uint8)
        for i in range(0,self.ntpl):
            # Observed wavelengths
            wave_in = self.hdu['WAVE'].data[i,observed[i,:]].ravel()
            # Rebin the observed wavelength range
            flux[i,:], wave, dv = log_rebin([wave_in[0], wave_in[-1]],
                                            self.hdu['FLUX'].data[i,observed[i,:]].ravel(),
                                            velscale=self.velscale, log10=True, newRange=fullRange,
                                            wave_in_ang=True, unobs=min_flux-100.)
            # Find the unobserved pixels, set them to have 0. flux, and
            # flag them as having no data
            indx = numpy.where(flux[i,:] < min_flux-10.)
            flux[i,indx] = 0.0
            mask[i,indx] = self.tplbm.turn_on(mask[i,indx], 'NO_DATA')

            # Resample the spectral resolution by simple interpolation.
            # Select the good pixels
            indx = numpy.where(numpy.invert(flux[i,:] < min_flux-10.))
            # define the interpolator (uses linear interpolation; k=1)
            interpolator = InterpolatedUnivariateSpline(self.hdu['WAVE'].data[i,:].ravel(),
                                                        self.hdu['SPECRES'].data[i,:].ravel(), k=1)
            # And then interpolate
            sres[i,indx] = interpolator(wave[indx])

            # Finally, rebin the masks:
            # Pixels outside the wavelength limits
            indx = self._rebin_masked(i, 'WAVE_INVALID', fullRange, rmsk_lim=0.1)
            mask[i,indx] = self.tplbm.turn_on(mask[i,indx], 'WAVE_INVALID')
            # Pixels below the flux limit
            indx = self._rebin_masked(i, 'FLUX_INVALID', fullRange, rmsk_lim=0.1)
            mask[i,indx] = self.tplbm.turn_on(mask[i,indx], 'FLUX_INVALID')
            # Pixels that required an extrapolation of the spectral
            # resolution
            indx = self._rebin_masked(i, 'SPECRES_EXTRAP', fullRange, rmsk_lim=0.1)
            mask[i,indx] = self.tplbm.turn_on(mask[i,indx], 'SPECRES_EXTRAP')
            # Pixels that had a spectral resolution that was too low to
            # match the galaxy resolution
            indx = self._rebin_masked(i, 'SPECRES_LOW', fullRange, rmsk_lim=0.1)
            mask[i,indx] = self.tplbm.turn_on(mask[i,indx], 'SPECRES_LOW')

        print('... done')

        # On input, the wavelength is now common to all spectra, such
        # that this HDU is a vector, not an array.

#        pyplot.plot(oldwave, oldflux)
#        pyplot.plot(wave, flux[0,:], 'g')
#        pyplot.show()

        # Normalize the templates to the mean flux value after excluding
        # any flagged pixels.
        indx = numpy.where( numpy.invert(self.tplbm.flagged(mask)) )
        flux_norm = numpy.mean(flux[indx])
        flux /= flux_norm

        # Reset the HDUList object
        self._reset_hdu(wave, flux, mask, sres, self.hdu['SIGOFF'].data)
        self.processed = True

        # Update the header with the redshift offset, the flux
        # normalization, and flag the data as having been prepared for
        # fitting the DRP data
        self.hdu[0].header['ZGUESS'] = (redshift, 'Guess redshift used')
        self.hdu[0].header['FLXNRM'] = (flux_norm, 'Flux normalization')
        self.hdu[0].header['TPLPROC'] = 1


    def _write_hdu(self):
        """Write the HDUList to the output file."""
        print('Writing: {0}'.format(self.file_path()))
        self.hdu.writeto(self.file_path())


    def file_name(self):
        """Return the name of the processed file."""
        return self.processed_file


    def file_path(self):
        """Return the full path to the processed file."""
        if self.directory_path is None or self.processed_file is None:
            return None
        return os.path.join(self.directory_path, self.processed_file)


    def read_raw_template_library(self, library_key=None, tpllib_list=None, dapsrc=None):
        """
        Read the identified template library.  If all the arguments are
        the default, the preset attributes from the initialization of
        the object are used.

        Args:
            library_key (str): (Optional) Keyword selecting the library
                to use.
            tpllib_list (list): (Optional) List of
                :class:`mangadap.util.par.TemplateLibraryParSet` objects
                that define the parameters required to read and
                interpret a template library.
            dapsrc (str): (Optional) Root path to the DAP source
                directory.  If not provided, the default is defined by
                :func:`mangadap.util.defaults.default_dap_source`.

        Raises:
            Exception: Raised if the template library is to be processed
                and the name of the prepared template library is not
                provided or cannot be created using the default
                algorithm.

        """
        if library_key is not None:
            # Redefine the library
            self._define_library(library_key, tpllib_list=tpllib_list, dapsrc=dapsrc)
        self._read_raw()


    def process_template_library(self, library_key=None, tpllib_list=None, dapsrc=None, drpf=None,
                                 velscale=None, sres=None, velocity_offset=0.0, dapver=None,
                                 analysis_path=None, directory_path=None, processed_file=None,
                                 force=False):
        """
        Process the template library for use in analyzing object
        spectra.  Primary steps are to:

            - Read the raw 1D fits files; see :func:`_read_raw`.

            - Convert the wavelengths to vacuum, if necessary; see
              :func:`mangadap.util.idlutils.airtovac`.

            - Mask wavelengths outside the rest wavelength range of the
              DRP spectrum, due to need to extrapolate these values; see
              :func:`mangadap.util.bitmasks.HDUList_mask_wavelengths`.

            - Match the spectral resolution of the template to that of
              the DRP spectra; see
              :func:`mangadap.util.instrument.match_spectral_resolution`.

            - Mask the template pixels where the spectral resolution was
              too low to match to the DRP spectra; see
              :func:`mangadap.util.bitmasks.BitMask.turn_on`.

            - Force a common wavelength range and sampling for all
              templates, where the sampling is forced to match the
              sampling of the DRP spectra; see
              :func:`mangadap.util.instrument.log_rebin_pix`.  The masks
              are appropriately resampled as well; see
              :func:`_rebin_masked`.

        .. warning::

            The routine **does not** check that that an existing
            processed file or the existing object has been processed
            using the same drpfile, velocity_offset, velscale, or sres
            input.  If unsure, use force=True.
        
        Args:

            library_key (str): (Optional) Keyword selecting the library
                to use; default is to use existing :attr:`library` and
                :attr:`libparset`.
            tpllib_list (list): (Optional) List of
                :class:`mangadap.util.par.TemplateLibraryParSet` objects
                that define the parameters required to read and
                interpret a template library.  Input ignored if
                *library_key* is None.
            dapsrc (str): (Optional) Root path to the DAP source
                directory.  If not provided, the default is defined by
                :func:`mangadap.util.defaults.default_dap_source`.
                Input ignored if *library_key* is None.
            drpf (:class:`mangadap.drpfile.drpfile`): (Optional) DRP
                file (object) with which the template library is
                associated for analysis.  **If not provided**, the user
                must define *velscale* and *sres* such that the library
                can be processed; and the user must provide
                *directory_path* and *processed_file* such that the
                output file can be written.
            velscale (float): (Optional) The velocity scale of the
                spectrum that will be analyzed with the library.  This
                value takes precedence over the value provided by the
                DRP file object.
            sres (:class:`mangadap.util.instrument.spectral_resolution`):
                (Optional) The object is used simply to access the
                spectral resolution and associated wavelength coordinate
                vector needed when matching the spectral resolution of
                the template library.  This takes prededence over the
                values provided by the DRP file object.
            velocity_offset (float): (Optional) Velocity offset to use
                when matching the spectral resolution between the
                template library and the galaxy spectra.
            dapver (str): (Optional) DAP version, which is used to
                define the default DAP analysis path.  Default is
                defined by
                :func:`mangadap.util.defaults.default_dap_version`
            analysis_path (str): (Optional) The path to the top level
                directory containing the DAP output files for a given
                DRP and DAP version.  Default is defined by
                :func:`mangadap.util.defaults.default_analysis_path`.
            directory_path (str): (Optional) The exact path for the
                processed template library file.  Default is defined by
                :func:`mangadap.util.defaults.default_dap_directory_path`.
            processed_file (str): (Optional) The name of the file
                containing (to contain) the prepared template library
                output file.  The file should be found at
                :attr:`directory_path`/:attr:`processed_file`.  Default
                is defined by
                :func:`mangadap.util.defaults.default_template_library_file`.
            force (bool): (Optional) Force the template library to be
                processed, even if the prepared template library already
                exists.

        Raises:
            Exception: Raised if the output file name is not defined.

        """
        # Reset the processing attributes
        if drpf is not None:
            self.drpf = drpf
        if velscale is not None:
            self.velscale = velscale
        if sres is not None:
            self.sres = sres
        if velocity_offset is not None:
            self.velocity_offset = velocity_offset

        # Check that the library can be processed
        if self.drpf is None and (self.velscale is None or self.sres is None):
            raise Exception('Cannot process template library without providing a DRP file object '\
                            ' or a new velocity scale and spectral resolution.')

        # Set the paths if possible
        directory_path = self.directory_path if directory_path is None else directory_path
        processed_file = self.processed_file if processed_file is None else processed_file
        if self._can_set_paths(directory_path, processed_file, quiet=True):
            self._set_paths(dapver, analysis_path, directory_path, processed_file)

        # Check that the path for or to the file is defined
        ofile = self.file_path()
        if ofile is None:
            raise Exception('File path for output file is undefined!')

#        if not force and self.velocity_offset is not None \
#           and not numpy.isclose(velocity_offset, self.velocity_offset):
#            print('Forcing processing due to change in velocity offset.')
#            force = True

        # Read and use a pre-existing file
        if os.path.isfile(ofile) and not force:
            print('Using existing file: {0}'.format(self.processed_file))
            self.hdu = fits.open(ofile)
            self.file_list = glob.glob(self.libparset['file_search'])
            self.ntpl = self.hdu['FLUX'].data.shape[0]
            self.processed = True
            return

        # Warn the user that the file will be overwritten
        if os.path.isfile(ofile):
            print('WARNING: Overwriting existing file: {0}'.format(self.processed_file))
            remove(ofile)

        self._read_raw()                        # Read the raw data
        self._process_library()                 # Process the library
        self._write_hdu()                       # Write the fits file


