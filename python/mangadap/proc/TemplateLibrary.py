"""

Class that reads and prepares template libraries for use with the MaNGA
DAP.  See :func:`mangadap.util.defaults.default_template_libraries`
for the list of default template libraries.

The "raw" template libraries are expected to consist of 1D fits files
that are linearly sampled in wavelength.  The fits files must have
CRVAL1, CRPIX1, and CDELT1 keywords used to define the wavelength
coordinates of each pixel::
 
    wave = (numpy.arange(1,npixels+1) - CRPIX1) * CDELT1 + CRVAL1
 
The reference frame of the template wavelengths must also be defined as
either vacuum or air.  It is expected that the DRP spectra are in vacuum
wavelengths.  The DAP will use :func:`mangadap.util.idlutils.airtovac`
to convert the template wavelengths to vacuum.  Finally, one can specify
that the template library is only valid within a certain wavelength
range and above a certian flux limit; see
:class:`mangadap.util.par.TemplateLibraryParSet`.

Preparation of the template library for use in fitting the DRP spectra
consists of the following steps:

    - Read the raw 1D fits files; see :func:`_read_raw`.

    - Convert the wavelengths to vacuum, if necessary; see
      :func:`mangadap.util.idlutils.airtovac`.

    - Mask wavelengths outside the rest wavelength range of the DRP
      spectrum, due to need to extrapolate these values; see
      :class:`mangadap.util.bitmasks.TemplateLibraryBitMask` and
      :func:`mangadap.util.bitmasks.HDUList_mask_wavelengths`.

    - Match the spectral resolution of the template to that of the DRP
      spectra; see
      :func:`mangadap.util.instrument.match_spectral_resolution`.  The
      detailed match of the spectral resolution should account for any
      redshift between the template library and the DRP spectra.  Note
      the usage examples below assume the default velocity offset of the
      DRP spectra, which is 0 km/s.

    - Mask the template pixels where the spectral resolution was too low
      to match to the DRP spectra; see
      :func:`mangadap.util.bitmasks.BitMask.turn_on`.

    - Force a common wavelength range and sampling for all templates,
      where the sampling is forced to match the sampling of the DRP
      spectra; see :func:`mangadap.util.instrument.log_rebin`.  The
      masks are appropriately resampled as well.
        
.. todo::

    - Allow preparation of a template library independent of a DRP file.

*Source location*:
    $MANGADAP_DIR/python/mangadap/proc/TemplateLibrary.py

*Python2/3 compliance*::

    from __future__ import division
    from __future__ import print_function
    from __future__ import absolute_import
    
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
        :class:`mangadap.util.par.ParSet` class.

.. _astropy.io.fits.hdu.hdulist.HDUList: http://docs.astropy.org/en/v1.0.2/io/fits/api/hdulists.html
.. _glob.glob: https://docs.python.org/3.4/library/glob.html

"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

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
    Object used to read, prepare, and store template libraries used with
    the MaNGA DAP.

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
        tpllib_list (list): List of :class:`mangadap.util.par.ParSet`
            objects that define the parameters required to read and
            interpret a template library.  These
            :class:`mangadap.util.par.ParSet` objects should have been
            created using
            :class:`mangadap.util.par.TemplateLibraryParSet`.
        dapsrc (str): (Optional) Root path to the DAP source directory.
            If not provided, the default is defined by
            :func:`mangadap.util.defaults.default_dap_source`.
        drpf (:class:`mangadap.drpfile.drpfile`): (Optional) DRP file
            (object) with which the template library is associated for
            analysis
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
        processed_file (str): (Optional) The name of the file containing
            the prepared template library output file.  The file should
            be found at :attr:`directory_path`/:attr:`processed_file`.
            Default is defined by
            :func:`mangadap.util.defaults.default_template_library_file`.
        read (bool): (Optional) Flag to read the template library if it
            exists.  If :attr:`drpf` is defined and the template library
            file exists in the directory path, this template library
            file is read.  Otherwise, the template library data is read
            from the repository.
        process (bool): (Optional) If :attr:`drpf` is defined and the
            prepared template library does not exist, this will process
            the template library in preparation for use in fitting the
            provided DRP file.
        velocity_offset (float): (Optional) Velocity offset to use when
            matching the spectral resolution between the template
            library and the galaxy spectra.
        force (bool): (Optional) If :attr:`drpf` is define and *process*
            is True, this will force the template library to be
            processed, even if the prepared template library already
            exists.

    Attributes:
        version (str): Version number
        tplbm (BitMask): A BitMask object used to toggle mask values;
            see :func:`mangadap.util.bitmasks.TemplateLibraryBitMask`.
        library (str): Keyword of the selected library to use.
        file_search (str): The searchable string used to grab the list of
            template library files.  E.g.::
        
            /Users/westfall/Work/MaNGA/dap/trunk/external/templates/m11_miles/*.fits

        file_list (list): The list of files found using `glob.glob`_ and
            :attr:`file_search`.
        ntpl (int): Number of template spectra in the library
        fwhm (float): The FWHM of the resolution element in angstroms,
            expected to be constant as a function of wavelength.
        invac (bool): Flag that the wavelengths of the template library
            are in vacuum.
        drpf (:class:`mangadap.drpfile.drpfile`): DRP file (object) with
            which the template library is associated for analysis
        drpver (str): DRP version
        dapver (str): DAP version
        analysis_path (str): The path to the top level directory
            containing the DAP output files for a given DRP and DAP
            version.
        directory_path (str): The exact path to the DAP file.
        processed_file (str): The name of the file containing the
            prepared template library output file.  The file should be
            found at :attr:`directory_path`/:attr:`processed_file`.
        processed (bool): Flag that the template library has been
            prepared for use in the DAP.
        velocity_offset (float): Velocity offset to use when matching
            the spectral resolution between the template library and the
            galaxy spectra.
        hdu (`astropy.io.fits.hdu.hdulist.HDUList`_): HDUList read from
            the DAP file

    .. todo::

        - Add USER library key and needed functionality.
        - Allow for input of sres vector and sampling to process the
          file instead of the DRP file
        - Continue to allow drpver to be different from drpf.drpver?

    """
#    def __init__(self, library_key, drpf=None, drpver=None, dapver=None,
#                 analysis_path=None, directory_path=None, processed_file=None, read=True,
#                 process=True, velocity_offset=0.0, force=False):
    def __init__(self, library_key, tpllib_list=None, dapsrc=None, drpf=None, drpver=None,
                 dapver=None, analysis_path=None, directory_path=None, processed_file=None,
                 read=True, process=True, velocity_offset=0.0, force=False):

        self.version = '2.0'

        self.tplbm = TemplateLibraryBitMask()

        self.library = None
        self.file_search = None
        self.file_list = None
        self.ntpl = None
        self.fwhm = None
        self.invac = None
        self.wave_limit = numpy.array([None, None])
        self.flux_limit = None

        self.drpf = None
        self.drpver = None
        self.dapver = None
        self.analysis_path = None
        self.directory_path = None
        self.processed_file = None
        self.processed = False
        self.velocity_offset = None
        self.hdu = None

        self._define_library(library_key, tpllib_list=tpllib_list, dapsrc=dapsrc)

        if not read:
            print('Do not read.')
            return

        # Read the library
#        print('Reading: calling read_template_library.')
#        self.read_template_library(drpf=drpf, dapsrc=dapsrc, drpver=drpver, dapver=dapver,
        self.read_template_library(drpf=drpf, drpver=drpver, dapver=dapver,
                                   analysis_path=analysis_path, directory_path=directory_path,
                                   processed_file=processed_file, process=process,
                                   velocity_offset=velocity_offset, force=force)


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
        :attr:`library`, :attr:`file_search`, :attr:`fwhm`, and :attr:`invac`.

        Args:
            library_key (str): Keyword of the selected library.
                Available libraries are proved by
                :func:`mangadap.util.defaults.default_template_libraries`
            tpllib_list (list): List of
                :class:`mangadap.util.par.ParSet' objects that define
                the parameters required to read and interpret a template
                library.  These :class:`mangadap.util.par.ParSet'
                objects should have been created using
                :class:`mangadap.util.par.TemplateLibraryParSet`.
            dapsrc (str): (Optional) Root path to the DAP source
                directory.  If not provided, the default is defined by
                :func:`mangadap.util.defaults.default_dap_source`.

        Raises:
            KeyError: Raised if the selected keyword is not among the
                provided list or if the provided list has more than one
                identical keyword.

        """
        # Get the default libraries if no list provided
        if tpllib_list is None:
            tpllib_list = default_template_libraries(dapsrc=dapsrc)

        # Make sure the input tpllib_list has the right type, and force
        # it to be a list
        print(type(tpllib_list))
        print(type(tpllib_list) is 'mangadap.util.par.TemplateLibraryParSet')
        print(isinstance(tpllib_list, TemplateLibraryParSet))
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
        self.file_search = tpllib_list[indx]['file_search']
        self.fwhm = tpllib_list[indx]['fwhm']
        self.invac = tpllib_list[indx]['in_vacuum']
        self.wave_limit = tpllib_list[indx]['wave_limit']
        self.flux_limit = tpllib_list[indx]['lower_flux_limit']


    def _set_paths(self, drpf, drpver, dapver, analysis_path, directory_path, processed_file):
        """
        Set the I/O path to the processed template library.  Used to set
        :attr:`drpf`, :attr:`drpver`, :attr:`dapver`,
        :attr:`analysis_path`, :attr:`directory_path`,
        :attr:`processed_file`.  If not provided, the default for the
        latter is set using
        :func:`mangadap.util.defaults.default_template_library_file`.

        Args:
            drpf (:class:`mangadap.drpfile.drpfile`): DRP file (object)
                with which the template library is associated for
                analysis.  See :attr:`drpf`.
            drpver (str): DRP version. See :attr:`drpver`.
            dapver (str): DAP version. See :attr:`dapver`.
            analysis_path (str): The path to the top-level directory
                containing the DAP output files for a given DRP and DAP
                version.  See :attr:`analysis_path`.
            directory_path (str): The exact path to the DAP template
                library file.  See :attr:`directory_path`.
            processed_file (str): The name of the file with the prepared
                template library.  See
                :attr:`processed_filedirectory_path`.

        """
        # Keep the drpfile
        self.drpf = drpf

        # Setup the DAP version and paths (taken from dapfile.__init__())
        self.drpver = default_drp_version() if drpver is None else str(drpver)
        self.dapver = default_dap_version() if dapver is None else str(dapver)
        self.analysis_path = default_analysis_path(self.drpver, self.dapver) \
                             if analysis_path is None else str(analysis_path)
        self.directory_path = default_dap_directory_path(self.drpver, self.dapver, \
                                                         self.analysis_path, self.drpf.plate, \
                                                         self.drpf.ifudesign) \
                              if directory_path is None else str(directory_path)
        self.processed_file = default_template_library_file(self.drpf.plate, self.drpf.ifudesign, \
                                                            self.drpf.mode, self.library) \
                              if processed_file is None else str(processed_file)


    def _read_raw(self):
        """
        Read the 'raw' versions of the template library; i.e., the
        library before it has been resolution and sampling matched to a
        DRP file.
        """
        self.file_list = glob.glob(self.file_search)
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
        :math:`N_{\rm spec} \times N_{\rm pix}`.

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
            if self.flux_limit is not None:
                indx = numpy.invert( flux_ > self.flux_limit )
                mask[i,indx] = self.tplbm.turn_on(mask[i,indx], 'FLUX_INVALID')
            if self.wave_limit[0] is not None:
                indx = wave[i,:].ravel() < self.wave_limit[0]
                mask[i,indx] = self.tplbm.turn_on(mask[i,indx], 'WAVE_INVALID')
            if self.wave_limit[1] is not None:
                indx = wave[i,:].ravel() > self.wave_limit[1]
                mask[i,indx] = self.tplbm.turn_on(mask[i,indx], 'WAVE_INVALID')

        # Set the spectral resolution
        sres = wave/self.fwhm

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


    def _rebin_masked(self, i, flag, drpf_dv, fullRange, rmsk_lim=0.5):
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
            drpf_dv (float): The velocity scale of each pixel in the
                logarithmically binned DRP spectra in km/s.
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
                                      mask_ex, velscale=drpf_dv, log10=True, newRange=fullRange,
                                      wave_in_ang=True)
        return numpy.where(mask_ex > rmsk_lim)


    def _process_library(self):
        """
        Process the template library for use in fitting spectra in the
        provided DRP file.  Primary steps are to:

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

        .. todo::

            - Make the resampling more efficient?
        
        """
        # Read the raw data
        self._read_raw()

        # Convert to vacuum wavelengths
#        pyplot.plot(self.hdu['WAVE'].data[0,:], self.hdu['FLUX'].data[0,:]) 
        if not self.invac:
            self.hdu['WAVE'].data = airtovac(self.hdu['WAVE'].data)
#        pyplot.plot(self.hdu['WAVE'].data[0,:], self.hdu['FLUX'].data[0,:], 'g') 
#        pyplot.show()

        # Calculate the redshift to use to set the DRP observed
        # wavelength to the rest wavelength; important for matching the
        # spectral resolution of the template and galaxy spectra in the
        # same reference frame
        redshift = self.velocity_offset/constants.c.to('km/s').value

        # Make sure the DRP file is open
        if self.drpf.hdu is None:
            print('Opening DRP file ... ')
            self.drpf.open_hdu()
            print('... done.')

        # Mask wavelengths where the spectral resolution of the DRP
        # spectra will have to be extrapolated for cross-matching with
        # the template spectra, accounting for the redshift of the
        # galaxy.  drpf.hdu['WAVE'].data must be one-dimensional (true
        # as of DRP v1_4_0)
        print('Masking extrapolation wavelengths ... ')
        wavelim = numpy.array([ self.drpf.hdu['WAVE'].data[0]/(1.+redshift), \
                                self.drpf.hdu['WAVE'].data[-1]/(1.+redshift) ])
        self.hdu = HDUList_mask_wavelengths(self.hdu, self.tplbm, 'SPECRES_EXTRAP', wavelim,
                                            invert=True)
        print('... done.')

#        oldwave = numpy.copy(self.hdu['WAVE'].data[0,:]).ravel()
#        oldflux = numpy.copy(self.hdu['FLUX'].data[0,:]).ravel()
#        pyplot.plot(self.hdu['WAVE'].data[0,:], self.hdu['SPECRES'].data[0,:]) 
#        pyplot.plot(self.drpf.hdu['WAVE'].data, self.drpf.hdu['SPECRES'].data, 'r') 

        # Match the resolution of the templates to the galaxy data
        # accounting for the redshift of the galaxy.
        print('Matching spectral resolution ... ')
        self.hdu['FLUX'].data, self.hdu['SPECRES'].data, self.hdu['SIGOFF'].data, res_mask = \
            match_spectral_resolution(self.hdu['WAVE'].data, self.hdu['FLUX'].data,
                                      self.hdu['SPECRES'].data,
                                      self.drpf.hdu['WAVE'].data/(1.+redshift),
                                      self.drpf.hdu['SPECRES'].data, min_sig_pix=0.0,
                                      new_log10=True)
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
        # range for each spectrum, imposing the velocity sampling of the
        # DRP file.
        drpf_dv = spectrum_velocity_scale(self.drpf.hdu['WAVE'].data, log10=True)

        npix = self.hdu['FLUX'].data.shape[1]
        dw, m, dlogw, dv = log_rebin_pix([self.hdu['WAVE'].data[0,0], self.hdu['WAVE'].data[0,-1]],
                                         npix, velscale=drpf_dv, log10=True, newRange=fullRange)
        max_m = m
        for i in range(1,self.ntpl):
            dw, m, dlogw, dv = log_rebin_pix([self.hdu['WAVE'].data[i,0],
                                              self.hdu['WAVE'].data[i,-1]], npix, velscale=drpf_dv,
                                              log10=True, newRange=fullRange)
            if max_m < m:
                max_m = m

        print(fullRange, drpf_dv, max_m, self.ntpl, npix)

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
                                            velscale=drpf_dv, log10=True, newRange=fullRange,
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
            indx = self._rebin_masked(i, 'WAVE_INVALID', drpf_dv, fullRange, rmsk_lim=0.1)
            mask[i,indx] = self.tplbm.turn_on(mask[i,indx], 'WAVE_INVALID')
            # Pixels below the flux limit
            indx = self._rebin_masked(i, 'FLUX_INVALID', drpf_dv, fullRange, rmsk_lim=0.1)
            mask[i,indx] = self.tplbm.turn_on(mask[i,indx], 'FLUX_INVALID')
            # Pixels that required an extrapolation of the spectral
            # resolution
            indx = self._rebin_masked(i, 'SPECRES_EXTRAP', drpf_dv, fullRange, rmsk_lim=0.1)
            mask[i,indx] = self.tplbm.turn_on(mask[i,indx], 'SPECRES_EXTRAP')
            # Pixels that had a spectral resolution that was too low to
            # match the galaxy resolution
            indx = self._rebin_masked(i, 'SPECRES_LOW', drpf_dv, fullRange, rmsk_lim=0.1)
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


    def read_template_library(self, library_key=None, tpllib_list=None, dapsrc=None, drpf=None,
                              drpver=None, dapver=None, analysis_path=None, directory_path=None,
                              processed_file=None, process=True, velocity_offset=0.0, force=False):
        """
        Read the identified template library.  If all the arguments are
        the default, the preset attributes are used.  This is
        essentially the same call as the initialization of the class
        instance; see :class:`TemplateLibrary`.

        Args:
            library_key (str): (Optional) Keyword selecting the library
                to use.
            tpllib_list (list): List of
                :class:`mangadap.util.par.ParSet` objects that define
                the parameters required to read and interpret a template
                library.  These :class:`mangadap.util.par.ParSet`
                objects should have been created using
                :class:`mangadap.util.par.TemplateLibraryParSet`.
            dapsrc (str): (Optional) Root path to the DAP source
                directory.  If not provided, the default is defined by
                :func:`mangadap.util.defaults.default_dap_source`.
            drpf (:class:`mangadap.drpfile.drpfile`): (Optional) DRP
                file (object) with which the template library is
                associated for analysis
            drpver (str): (Optional) DRP version, which is used to
                define the default DAP analysis path.  Default is
                defined by
                :func:`mangadap.util.defaults.default_drp_version`
            dapver (str): (Optional) DAP version, which is used to
                define the default DAP analysis path.  Default is
                defined by
                :func:`mangadap.util.defaults.default_dap_version`
            analysis_path (str): (Optional) The path to the top level
                directory containing the DAP output files for a given
                DRP and DAP version.  Default is defined by
                :func:`mangadap.util.defaults.default_analysis_path`.
            directory_path (str): (Optional) The exact path to the DAP
                file.  Default is defined by
                :func:`mangadap.util.defaults.default_dap_directory_path`.
            processed_file (str): (Optional) The name of the file
                containing the prepared template library output file.
                The file should be found at
                :attr:`directory_path`/:attr:`processed_file`.  Default
                is defined by
                :func:`mangadap.util.defaults.default_template_library_file`.
            process (bool): (Optional) If :attr:`drpf` is defined and
                the prepared template library does not exist, this will
                process the template library in preparation for use in
                fitting the provided DRP file.
            velocity_offset (float): (Optional) Velocity offset to use
                when matching the spectral resolution between the
                template library and the galaxy spectra.
            force (bool): (Optional) If :attr:`drpf` is define and
                *process* is True, this will force the template library
                to be processed, even if the prepared template library
                already exists.

        Raises:
            Exception: Raised if the template library is to be processed
                and the name of the prepared template library is not
                provided or cannot be created using the default
                algorithm.

        """
        if library_key is not None:
            # Redefine the library
            self._define_library(library_key, tpllib_list=tpllib_list, dapsrc=dapsrc)

        if not process:
            # Just read the existing library
            self._read_raw()
            return

        # Use the existing attributes if nothing is input
        drpf = self.drpf if drpf is None else drpf
        drpver = self.drpver if drpver is None else drpver
        dapver = self.dapver if dapver is None else dapver
        analysis_path = self.analysis_path if analysis_path is None else analysis_path
        directory_path = self.directory_path if directory_path is None else directory_path
        processed_file = self.processed_file if processed_file is None else processed_file

        # Check that the name is either provided or can be created using the default
        if drpf is None and (processed_file is None or directory_path is None):
            raise Exception('To use default output file and output file path, must provide DRP' \
                            ' file object, which defines the plate and ifudesign.')

        # Set the paths
        print('Reading: set paths.')
        self._set_paths(drpf, drpver, dapver, analysis_path, directory_path, processed_file)
        print(self.directory_path)

        self.process_template_library(velocity_offset=velocity_offset, force=force)


    def process_template_library(self, velocity_offset=0.0, force=False):
        """
        Process the template library for use in fitting spectra in the
        provided DRP file.  Primary steps are to:

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
        
        Args:
            velocity_offset (float): (Optional) Velocity offset to use
                when matching the spectral resolution between the
                template library and the galaxy spectra.
            force (bool): (Optional) If :attr:`drpf` is define and
                *process* is True, this will force the template library
                to be processed, even if the prepared template library
                already exists.

        Raises:
            Exception: Raised if the output file is not defined.

        """

        # Read and use a pre-existing file
        ofile = self.file_path()

        if ofile is None:
            raise Exception('File path for output file is undefined!')

        if not force and self.velocity_offset is not None \
           and not numpy.isclose(velocity_offset, self.velocity_offset):
            print('Forcing processing due to change in velocity offset.')
            force = True

        if os.path.isfile(ofile) and not force:
            self.hdu = fits.open(ofile)
            self.processed = True
            return

        if os.path.isfile(ofile):
            print('WARNING: Overwriting existing file: {0}'.format(self.processed_file))
            remove(ofile)

        print('Processing...')
        self.velocity_offset = velocity_offset
        self._process_library()
        self._write_hdu()


