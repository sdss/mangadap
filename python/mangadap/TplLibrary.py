"""

Class that reads and manipulates template libraries for use with the
MaNGA DAP.  See
:func:`mangadap.util.defaults.available_template_libraries` for the list
of available template libraries.

The fits files must have CRVAL1, CRPIX1, and CDELT1 keywords used to
define the wavelength coordinates of each pixel::
 
        pix = numpy.arange(1,npixels+1)
        wave = (pix - CRPIX1) * CDELT1 + CRVAL1
 
The reference frame of the template wavelengths must also be defined as
either vacuum or air.  It is expected that the DRP spectra are in vacuum
wavelengths.  The DAP will use :func:`mangadap.util.idlutils.airtovac`
to convert the template wavelengths to vacuum.


        The detailed match of the spectral resolution should account for
        any redshift between the existing and target wavelength frame.
        E.g. if the target frame is for a set of galaxy spectra and the
        existing frame is for a set of templates.  This is not accounted
        for here, assuming that the wavelength vectors of both the
        existing and target resolutions are *in the same reference
        frame*!




*Source location*:
    $MANGADAP_DIR/python/mangadap/TplLibrary.py

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
    from os import environ
    from scipy import sparse
    from astropy.io import fits
    import time
    import numpy
    from mangadap.util.defaults import default_dap_source, default_drp_version, default_dap_version
    from mangadap.util.defaults import default_analysis_path, default_dap_directory_path
    from mangadap.util.defaults import available_template_libraries
    from mangadap.util.fileio import readfits_1dspec

*Class usage examples*:

    .. todo::

        Add some usage comments here!

*Revision history*:
    | **23 Apr 2015**: Implementation begun by K. Westfall (KBW)
    | **26 May 2015**: (KBW) Added some Sphinx documentation.

.. _astropy.io.fits.hdu.hdulist.HDUList: http://docs.astropy.org/en/v1.0.2/io/fits/api/hdulists.html

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
from mangadap.util.defaults import available_template_libraries, default_template_library_file
from mangadap.util.idlutils import airtovac
from mangadap.util.fileio import readfits_1dspec
from mangadap.util.bitmasks import TplLibraryBitMask, HDUList_mask_wavelengths
from mangadap.util.instrument import log_rebin, log_rebin_pix
from mangadap.util.instrument import match_spectral_resolution, spectrum_velocity_scale

#from matplotlib import pyplot

__author__ = 'Kyle B. Westfall'

class TplLibrary:
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
        drpf (:class:`mangadap.drpfile.drpfile`): (Optional) DRP file
            (object) with which the template library is associated for
            analysis
        dapsrc (str): (Optional) Root path to the DAP source directory.
            If not provided, the default is defined by
            :func:`mangadap.util.defaults.default_dap_source`.
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
        force (bool): (Optional) If :attr:`drpf` is define and *process*
            is True, this will force the template library to be
            processed, even if the prepared template library already
            exists.

    Attributes:
        library (str): Keyword of the selected library to use.
        files (str): The searchable string used to grab the list of
            template library files.  E.g.::
        
            /Users/westfall/Work/MaNGA/dap/trunk/external/templates/m11_miles/*.fits

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
        processed (bool): Flag that the template library has been
            prepared for use in the DAP.
        processed_file (str): The name of the file containing the
            prepared template library output file.  The file should be
            found at :attr:`directory_path`/:attr:`processed_file`.
        hdu (`astropy.io.fits.hdu.hdulist.HDUList`_): HDUList read from
            the DAP file

    .. todo::

        - Add USER library key and needed functionality.
        - Allow for input of sres vector and sampling to process the
          file instead of the DRP file

    """

    def __init__(self, library_key, drpf=None, dapsrc=None, drpver=None, dapver=None,
                 analysis_path=None, directory_path=None, processed_file=None, read=True,
                 process=True, velocity_offset=0.0, force=False):

        self.version = '2.0'

        self.tplbm = TplLibraryBitMask()

        self.library = None
        self.files = None
        self.fwhm = None
        self.invac = None

        self._define_library(library_key)

        self.drpf = None
        self.drpver = None
        self.dapver = None
        self.analysis_path = None
        self.directory_path = None
        self.processed_file = None
        self.processed = False
        self.velocity_offset = None
        self.hdu = None

        if not read:
            print('Do not read.')
            return

        # Read the library
        print('Reading: calling read_template_library.')
        self.read_template_library(drpf=drpf, dapsrc=dapsrc, drpver=drpver, dapver=dapver,
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


    def _define_library(self, library_key, dapsrc=None):
        """
        Select the library from the available list.  Used to set
        :attr:`library`, :attr:`files`, :attr:`fwhm`, and :attr:`invac`.

        Args:
            library_key (str): Keyword of the selected library.
                Available libraries are proved by
                :func:`mangadap.util.defaults.available_template_libraries`

        Raises:
            KeyError: Raised if the selected keyword is not among the
                available list.
        """
        keys, paths, lsf_fwhm, vac = available_template_libraries(dapsrc)

        if library_key not in keys:
            raise KeyError('{0} is not a valid template library!'.format(library_key))

        self.library = library_key
        indx = keys == self.library
        self.files = paths[indx][0]
        self.fwhm = lsf_fwhm[indx][0]
        self.invac = vac[indx][0]


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


#    def _need_to_reset_paths(self, drpf, drpver, dapver, analysis_path, directory_path,
#                             processed_file):
#        """
#        Determine if the paths need to be reset.
#
#        Args:
#            drpf (:class:`mangadap.drpfile.drpfile`): DRP file (object)
#                with which the template library is associated for
#                analysis.
#            drpver (str): DRP version
#            dapver (str): DAP version
#            analysis_path (str): The path to the top-level directory
#                containing the DAP output files for a given DRP and DAP
#                version.
#            directory_path (str): The exact path to the prepared
#                template library files.
#            processed_file (str): The name of the file containing
#                the prepared template library output file.  The file
#                should be found at
#                :attr:`directory_path`/:attr:`processed_file`.
#
#        Returns:
#            bool : Flag that paths need to be reset because one or more
#                of the arguments is not none.
# 
#        """
#        itms = [drpf, drpver, dapver, analysis_path, directory_path, processed_file]
#        return any(i is not None for i in itms)


    def _read_raw(self):
        """
        Read the 'raw' versions of the template library; i.e., the
        library before it has been resolution and sampling matched to a
        DRP file.
        """
        file_list = glob.glob(self.files)
        nf = len(file_list)
        print('Found {0} {1} templates'.format(nf, self.library))

        npix = self._get_nchannels(file_list)
        print('Maximum number of wavelength channels: {0}'.format(npix))

        self._build_raw_hdu(file_list, nf, npix)


    def _get_nchannels(self, file_list):
        max_npix = 0
        for f in file_list:
            if fits.getval(f, 'NAXIS') != 1:
                raise Exception('{0} is not one dimensional!'.format(f))
            npix = fits.getval(f, 'NAXIS1')
            if max_npix < npix:
                max_npix = npix
        return max_npix


    def _build_raw_hdu(self, file_list, nf, npix):
        """
        On output, the wavelength vectors are not necessarily all the same!
        """
        if self.hdu is not None:
            self.hdu.close()
            self.hdu = None

        print('Attempting to build raw list ...')
        # Add some header cards
        wave = numpy.zeros((nf, npix), dtype=numpy.float64)
        flux = numpy.zeros((nf, npix), dtype=numpy.float64)
        mask = numpy.zeros((nf, npix), dtype=numpy.uint8)
        soff = numpy.zeros(nf, dtype=numpy.float64)
#        ivar = numpy.zeros((nf, npix), dtype=numpy.float64)
#        ivar[:] = 1.0

        for i in range(0,nf):
            wave_, flux_ = readfits_1dspec(file_list[i])
            wave[i,0:wave_.size] = numpy.copy(wave_)
            flux[i,0:wave_.size] = numpy.copy(flux_)
            if wave_.size != npix:
                print('flagging')
                mask[i,wave_.size:] = self.tplbm.turn_on(mask[i,wave_.size:],'NODATA')

        sres = wave/self.fwhm

        # (Re)Set the header
        self._reset_hdu(wave, flux, mask, sres, soff)

        # Add some keywords to the header
        self.hdu[0].header['TPLPROC'] = (0, 'Flag that library has been processed for DAP use')
        print('... done')


    def _reset_hdu(self, wave, flux, mask, sres, soff):
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
        indx = numpy.where(numpy.invert(self.tplbm.flagged(self.hdu['MASK'].data, flag=flag)))
        return numpy.array([ numpy.amin(self.hdu['WAVE'].data[indx]),
                             numpy.amax(self.hdu['WAVE'].data[indx])])


    def _rebin_masked(self, i, flag, drpf_dv, fullRange):
        mask_ex = self.tplbm.flagged(self.hdu['MASK'].data[i,:], flag=flag).astype(numpy.float64)
        mask_ex, wave, dv = log_rebin([self.hdu['WAVE'].data[i,0], self.hdu['WAVE'].data[i,-1]],
                                      mask_ex, velscale=drpf_dv, log10=True, newRange=fullRange,
                                      wave_in_ang=True)
        return numpy.where(mask_ex > 0.5)


    def _process_library(self):

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

        oldflux = numpy.copy(self.hdu['FLUX'].data[0,:]).ravel()
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
#        pyplot.plot(self.hdu['WAVE'].data[0,:], self.hdu['SPECRES'].data[0,:], 'g') 
#        pyplot.show()
        print('... done')

        # Mask any pixels where the template resolution was too low to
        # match to the galaxy resolution
        print('Masking low spectral resolution ... ')
        self.hdu['MASK'].data[res_mask == 1] = \
            self.tplbm.turn_on(self.hdu['MASK'].data[res_mask == 1], 'SPECRES_LOW')
        print('... done')

#        pyplot.plot(self.hdu['WAVE'].data[0,:], oldflux)
#        pyplot.plot(self.hdu['WAVE'].data[0,:], self.hdu['FLUX'].data[0,:], 'g')
#        pyplot.show()

        # Determine the wavelength range that encloses all spectra
        print('Matching sampling ... ')
        fullRange = self._wavelength_range(flag='NODATA')

        # Get the maximum number of pixels needed to cover this full
        # range for each spectrum, imposing the velocity sampling of the
        # DRP file.
        drpf_dv = spectrum_velocity_scale(self.drpf.hdu['WAVE'].data, log10=True)

        nspec = self.hdu['FLUX'].data.shape[0]
        npix = self.hdu['FLUX'].data.shape[1]
        dw, m, dlogw, dv = log_rebin_pix([self.hdu['WAVE'].data[0,0], self.hdu['WAVE'].data[0,-1]],
                                         npix, velscale=drpf_dv, log10=True, newRange=fullRange)
        max_m = m
        for i in range(1,nspec):
            dw, m, dlogw, dv = log_rebin_pix([self.hdu['WAVE'].data[i,0],
                                              self.hdu['WAVE'].data[i,-1]], npix, velscale=drpf_dv,
                                              log10=True, newRange=fullRange)
            if max_m < m:
                max_m = m

        print(fullRange, drpf_dv, max_m, nspec, npix)

        # Resample the templates to match the object spectra: On input,
        # tpl_wave expected to have the same size as tpl_flux
        # (wavelength defined separately for each spectrum); on output,
        # tpl_wave and tpl_sres are single vectors common to ALL
        # template spectra.  Sampling is forced to match galaxy
        # sampling, via velscale
        min_flux = numpy.amin(self.hdu['FLUX'].data.ravel())

        observed = numpy.invert(self.tplbm.flagged(self.hdu['MASK'].data, flag='NODATA'))

        flux = numpy.zeros((nspec, max_m), dtype=numpy.float64)
        sres = numpy.zeros((nspec, max_m), dtype=numpy.float64)
        mask = numpy.zeros((nspec, max_m), dtype=numpy.uint8)
        for i in range(0,nspec):
            wave_in = self.hdu['WAVE'].data[i,observed[i,:]].ravel()
            flux[i,:], wave, dv = log_rebin([wave_in[0], wave_in[-1]],
                                            self.hdu['FLUX'].data[i,observed[i,:]].ravel(),
                                            velscale=drpf_dv, log10=True, newRange=fullRange,
                                            wave_in_ang=True, unobs=min_flux-100.)
            indx = numpy.where(flux[i,:] < min_flux-10.)
            flux[i,indx] = 0.0
            mask[i,indx] = self.tplbm.turn_on(mask[i,indx], 'NODATA')

            indx = numpy.where(numpy.invert(flux[i,:] < min_flux-10.))
            interpolator = InterpolatedUnivariateSpline(self.hdu['WAVE'].data[i,:].ravel(),
                                                        self.hdu['SPECRES'].data[i,:].ravel(), k=1)
            sres[i,indx] = interpolator(wave[indx])

            indx = self._rebin_masked(i, 'SPECRES_EXTRAP', drpf_dv, fullRange)
            mask[i,indx] = self.tplbm.turn_on(mask[i,indx], 'SPECRES_EXTRAP')

            indx = self._rebin_masked(i, 'SPECRES_LOW', drpf_dv, fullRange)
            mask[i,indx] = self.tplbm.turn_on(mask[i,indx], 'SPECRES_LOW')
        print('... done')

        # Normalize the templates and save the normalization
        indx = numpy.where( numpy.invert(self.tplbm.flagged(mask)) )
        flux_norm = numpy.mean(flux[indx])
        flux[indx] /= flux_norm

        # Reset the HDUList object
        self._reset_hdu(wave, flux, mask, sres, self.hdu['SIGOFF'].data)

        # Update the header
        self.hdu[0].header['ZGUESS'] = (redshift, 'Guess redshift used')
        self.hdu[0].header['FLXNRM'] = (flux_norm, 'Flux normalization')
        self.hdu[0].header['TPLPROC'] = 1

        # Write the template data to a fits file
        self._write_hdu()


    def _write_hdu(self):
        print('Writing: {0}'.format(self.file_path()))
        self.hdu.writeto(self.file_path())


    def file_name(self):
        return self.processed_file


    def file_path(self):
        if self.directory_path is None or self.processed_file is None:
            return None
        return os.path.join(self.directory_path, self.processed_file)


    def read_template_library(self, library_key=None, drpf=None, dapsrc=None, drpver=None,
                              dapver=None, analysis_path=None, directory_path=None,
                              processed_file=None, process=True, velocity_offset=0.0, force=False):
        """
        Read the identified template library.  If all the arguments are
        the default, the preset attributes are used.

        Args:
            library_key (str): (Optional) Keyword selecting the library
                to use.
            drpf (:class:`mangadap.drpfile.drpfile`): (Optional) DRP
                file (object) with which the template library is
                associated for analysis
            dapsrc (str): (Optional) Root path to the DAP source
                directory.  If not provided, the default is defined by
                :func:`mangadap.util.defaults.default_dap_source`.
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
            self._define_library(library_key, dapsrc)

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
#        self._write_hdu()


