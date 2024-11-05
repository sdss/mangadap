# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
r"""
Class that reads and prepares template libraries for use in fitting the
stellar-continuum of a spectrum.  See
:func:`available_template_libraries` for the list of default template
libraries.

The "raw" template libraries are expected to consist of 1D fits files.
The wavelength sampling can be either linear or logarithmic in
wavelength.  The reference frame of the template wavelengths must be
defined as either vacuum or air.  It is expected that the object spectra
to be fit are calibrated to vacuum wavelengths.  When preparing the
template spectra for analysis, this class will use
`pydl.goddard.astro.airtovac`_ to convert the template wavelengths to
vacuum.  Finally, one can specify that the template library is only
valid within a certain wavelength range and above a certian flux limit;
see :class:`TemplateLibraryDef`.

Preparation of the template library for use in stellar-continuum fitting
consists of a number of steps as described in the documentation of
:func:`TemplateLibrary.process_template_library`.

A template library that has been prepared for analysis is automatically
written to disk for later recovery.

Two support classes are also provided. One is a derived
:class:`mangadap.par.parset.KeywordParSet` instance that provides the
defining parameters of a DAP template library. The second is a
derived :class:`mangadap.util.bitmask.BitMask` instance that defines
the bitmasks for the template library spectra.

.. _templatelibrary-usage:

Class usage examples
--------------------

TODO: Update this for new datacube usage

Assuming you have the default directory structure setup, you can do::

    # Imports
    from mangadap.drpfits import DRPFits
    from mangadap.proc.templatelibrary import TemplateLibrary
    from matplotlib import pyplot

    # Define the DRP file
    drpf = DRPFits(7495, 12703, 'CUBE')

    # Build the template library
    tpl_lib = TemplateLibrary('M11MILES', drpf=drpf, directory_path='.')
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
to, e.g., the DRP file; see :class:`mangadap.util.drpfits.DRPFits`.

If you do not want to process the template (match the spectral
resolution and sampling to that of the DRP data), you can force the
:class:`TemplateLibrary` object to only provide the "raw" spectra::

    # Imports
    from mangadap.proc.templatelibrary import TemplateLibrary
    from matplotlib import pyplot

    # Read the raw template library
    tpl_lib = TemplateLibrary('M11MILES', process=False)
    # Nothing should be written to disk

    # Plot one of the spectra
    pyplot.plot(tpl_lib.hdu['WAVE'].data[0,:], tpl_lib.hdu['FLUX'].data[0,:])
    pyplot.show()

Note that in the previous example, the wavelength data was already
one-dimensional, whereas for the raw library, the wavelength vector
can be spectrum-dependent.

In the above examples, the user has not provided a list of template
libraries, meaning that the default set available to the DAP is
used.  The default set is defined in
:func:`available_template_libraries`.  If you want to use your own
template library, you have to define its parameters using
:class:`TemplateLibraryDef`.  Currently, the template library
spectra are expected to be 1D fits files with WCS header keywords
defining the wavelength solution; see above.  Using an existing DAP
library as an example::

    # Imports
    import os
    from mangadap.proc.templatelibrary import TemplateLibraryDef
    from mangadap.proc.templatelibrary import TemplateLibrary
    from mangadap.config.defaults import dap_data_root
    from matplotlib import pyplot

    # Define the search string for the library
    search_str = os.path.join(dap_data_root(), 'stellar_templates/miles/*.fits')

    # Define the template library parameters
    tpl_par = TemplateLibraryDef(key='MILES',    # Unique keyword for the library
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

You can also process the spectra to a user-provided resolution and
pixel sampling, toggle on/off the renormalization of the library to
a mean flux of unity and define various paths if you're not working
within the nominal DAP directory structure.  See the optional
instantiation arguments for :class:`TemplateLibrary`.

----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import inspect
from pathlib import Path
import warnings
import logging

from IPython import embed

import numpy
from scipy import interpolate

from astropy.io import fits
import astropy.constants

from pydl.goddard.astro import airtovac

from ..par.parset import KeywordParSet
from ..config import defaults
from ..util.dapbitmask import DAPBitMask
from ..util.log import log_output
from ..util.fileio import read_template_spectrum, writefits_1dspec, create_symlink
from ..util.fitsutil import DAPFitsUtil
from ..util import sampling
from ..util.resolution import SpectralResolution, match_spectral_resolution
from ..util.parser import DefaultConfig
from .util import HDUList_mask_wavelengths


class TemplateLibraryDef(KeywordParSet):
    """
    Class with parameters used to define the template library.
    See :class:`mangadap.par.parset.ParSet` for attributes.

    The defined parameters are:

    .. include:: ../tables/templatelibrarydef.rst
    """
    def __init__(self, key='MILESHC', file_search='miles_cluster/*.fits', fwhm=2.50, sres_ext=None,
                 in_vacuum=False, wave_limit=None, lower_flux_limit=None, log10=False): 

        # Use the signature to get the parameters and the default values
        sig = inspect.signature(self.__class__)
        pars = list(sig.parameters.keys())
        defaults = [sig.parameters[key].default for key in pars]

        # Remaining definitions done by hand
        in_fl = [int, float]
        arr_like = [numpy.ndarray, list]
        
        values = [key, file_search, fwhm, sres_ext, in_vacuum, wave_limit, lower_flux_limit, log10]
        dtypes = [str, str, in_fl, str, bool, arr_like, in_fl, bool]
        descr = ['Keyword to distinguish the template library.',
                 'Search pattern used to find the 1D fits spectra to include in the ' \
                    'template library.  The search string must either include the full path ' \
                    'or be within the template directory of the DAP source distribution.',
                 'FWHM of the resolution element in angstroms.',
                 'Extension in the fits files with measurements of the spectral resolution as ' \
                    'a function of wavelength.  If present, this supercedes any provided FWHM.',
                 'Flag that the wavelengths of the spectra are in vacuum, not air.',
                 'Two-element array with the starting and ending wavelengths for the valid ' \
                    'spectral range of the templates.',
                 'Minimum valid flux in the template spectra.',
                 'Flag that the template spectra have been binned logarithmically in wavelength.']

        super().__init__(pars, values=values, defaults=defaults, dtypes=dtypes, descr=descr)
        self._validate()

    @classmethod
    def from_key(cls, key):
        """
        Instantiate the template parameters using a keyword only.
        """
       # Check the configuration files exist
        search_dir = defaults.dap_config_root() / 'spectral_templates'
        ini_files = sorted(list(search_dir.glob('*.ini')))
        if len(ini_files) == 0:
            raise IOError(f'Could not find any configuration files in {search_dir} !') 
        for i in range(len(ini_files)):
            cnfg = DefaultConfig(f=ini_files[i], interpolate=True)
            if key != cnfg['key']:
                continue
            wave_limit = None if cnfg['wave_limit'] is None \
                            else numpy.array(cnfg.getlist('wave_limit', evaluate=True))
            return cls(key=cnfg['key'], file_search=cnfg['file_search'],
                       fwhm=cnfg.getfloat('fwhm'), sres_ext=cnfg['sres_ext'],
                       in_vacuum=cnfg.getbool('in_vacuum'), wave_limit=wave_limit,
                       lower_flux_limit=cnfg.getfloat('lower_flux_limit'),
                       log10=cnfg.getbool('log10', default=False))

        # If the code makes it here, the keyword did not select a known template library
        raise ValueError(f'Unknown template library: {key}')

    @classmethod
    def from_dict(cls, d):
        """
        Instantiate from a dictionary.

        If the dictionary only contains the keyword for the template library,
        the function assumes the user is trying to define the library based on
        that keyword.  Otherwise, the function calls the base-class
        ``from_dict`` method.
        """
        # First make sure that the dictionary defines the library keyword
        if 'key' not in d.keys():
            raise ValueError('Template library instantiation must include a keyword identifier.')
        if list(d.keys()) == ['key']:
            return cls.from_key(d['key'])
        return super().from_dict(d)

    def _validate(self):
        """ 
        Validate the template library parameters.
        """
        if self['key'] is None:
            raise ValueError('Must define the template library keyword identifier.')
        if self['file_search'] is None:
            raise ValueError('Must define the search pattern for template library files.')
        if self['fwhm'] is None and self['sres_ext'] is None:
            raise ValueError('Template library must provide either \'fwhm\' or \'sres_ext\'.')

    def get_file_list(self):
        """
        Use the search string to find the template library fits files.

        The search string must be via the default location in the DAP
        distribution or the full path to a local library.

        The file list read by the search key is sorted for
        consistency.

        Returns:
            :obj:`list`: The sorted list of files.

        Raises:
            ValueError:
                Raised if no files are found.
        """
        # Try the DAP directory first
        root = defaults.dap_data_root() / 'spectral_templates' / self['file_search']
        if root.parent.exists():
            files = sorted(list(root.parent.glob(root.name)))
        if not root.parent.exists() or len(files) == 0:
            # Then try the provided path directly
            root = Path(self['file_search']).resolve()
            if root.parent.exists():
                files = sorted(list(root.parent.glob(root.name)))
        if not root.parent.exists() or len(files) == 0:
            raise FileNotFoundError('Unable to find template library directory or any valid files '
                                    'as given or in the DAP source distribution.')
        return files


class TemplateLibraryBitMask(DAPBitMask):
    """
    Derived class that specifies the mask bits for the template library
    data.  The maskbits defined are:

    .. include:: ../tables/templatelibrarybitmask.rst
    """
    cfg_root = 'spectral_template_bits'


class TemplateLibrary:
    r"""
    TODO: Revisit these docs.

    Object used to read, store, and prepare template libraries for
    use in analyzing object spectra.

    The default list of available libraries provided by the MaNGA DAP
    defined in :func:`available_template_libraries`. The user can
    provide their own library for use with this class provided they
    are contained in 1D fits spectra, sampled linearly in wavelength
    with the wavelength coordinates available via the WCS keywords
    (``CRPIX1``, ``CRVAL1``, ``CDELT1``), and they have an
    appropriately defined spectral resolution. See
    :class:`TemplateLibraryDef` and :func:`_build_raw_hdu`.

    The class is optimized for use in analyzing MaNGA DRP files;
    however, one can provide the necessary information so that the class
    can be used with a non-DRP spectrum.  In the latter case, the user
    must supply the velocity scale of the pixel for the logarithmically
    resampled template library, and a
    :class:`mangadap.util.resolution.SpectralResolution` object the
    defines the instrumental resolution of the spectrum/spectra to be
    analyzed.

    .. todo::

        - Docs are out of date.
        - Update the list of attributes.
        - Allow for velocity scale input instead of spectral_step?
        - Only works with DRP files that have log-linear wavelength
          binning!
        - Allow to process, where process is just to change the sampling
          or the resolution (not necessarily both).
        - Need to make this more general, removing all dependence on
          DRPFits object.  This would simplify the functionality to
          change how the resolution and sampling matching is specified.

    On initialization, if the datacube object is not provided (is
    None), the default behavior is to read the raw template library
    if read=True. If the datacube is provided, the routine will check
    for the resolution matched fits file; if it doesn't exist and
    read is ``True``, it will prepare the template library for use in
    analyzing the datacube and write the prepared library file (if
    hardcopy is True). If overwrite=True, the preparation and writing
    of the template library will be done even if the library already
    exists.

    .. todo::

        Change `sres` to be a standard input object?

    Args:
        library (:obj:`str`, :class:`TemplateLibraryDef`):
            A keyword selecting one of the template libraries distributed with
            the DAP or the definition of a user-provided library.
        cube (:class:`mangadap.datacube.datacube.DataCube`, optional):
            The datacube to be analyzed using the template library.
        match_resolution (:obj:`bool`, optional):
            Match the spectral resolution of the template library to
            the resolution provided by the ``cube``; the latter must
            be provided for this argument to have any use.
        velscale_ratio (:obj:`int`, optional):
            Resample the template spectra such that the ratio of the
            pixel scale in the provided ``cube`` is this many times
            larger than the pixels in the resampled template
            spectrum.
        sres (:class:`mangadap.util.resolution.SpectralResolution`, optional):
            The object is used simply to access the spectral
            resolution and associated wavelength coordinate vector
            needed when matching the spectral resolution of the
            template library; this is used in lieu of a datacube.
        velocity_offset (:obj:`float`, optional):
            Velocity offset to use when matching the spectral
            resolution between the template library and the galaxy
            spectra.
        min_sig_pix (:obj:`float`, optional):
            Minimum value of the LSF standard deviation in pixels
            allowed before assuming the kernel is a Delta function.
            See
            :func:`mangadap.util.resolution.match_spectral_resolution`.
        no_offset (:obj:`bool`, optional):
            Prohibit the spectral resolution matching to introduce an
            offset in spectral resolution required to match the
            relative shape of the resolution vectors. See
            :func:`mangadap.util.resolution.match_spectral_resolution`
            and
            :func:`mangadap.util.resolution.SpectralResolution.GaussianKernelDifference`.
        spectral_step (:obj:`float`, optional):
            Logarithmic (``log=True``) or linear (``log=False``) step
            in wavelength for the template library.
        log (:obj:`bool`, optional):
            Flag to force the library to be logarithmically sampled in
            wavelength.
        wavelength_range (array-like, optional):
            Force the template library to cover this spectral range.
            Unobserved spectral regions will be masked.
        renormalize (:obj:`bool`, optional):
            After processing, renormalize the flux of the full
            library to unity. The relative fluxes of the templates
            are maintained.
        output_path (:obj:`str`, `Path`_, optional):
            The path for the output file.  If None, the current working
            directory is used.
        output_file (:obj:`str`, optional):
            The name of the output file. The full path of the
            output file will be :attr:`directory_path`/:attr:`output_file`.
            If None, the default is to combine ``cube.output_root`` and the
            library key; if ``cube`` is None, it is simply the library key.
        read (:obj:`bool`, optional):
            Flag to read the template library data.
        process (:obj:`bool`, optional):
            Process (spectral resolution and sampling operations) the
            template library.  See :func:`process_template_library`.
        hardcopy (:obj:`bool`, optional):
            Flag to keep a hardcopy of the processed template library.
        symlink_dir (:obj:`str`, optional):
            Create a symlink to the file in this directory. If None,
            no symbolic link is created.
        overwrite (:obj:`bool`, optional):
            Overwrite any saved, processed template library file.
        checksum (:obj:`bool`, optional):
            Use the checksum in the fits header to confirm that the
            data has not been corrupted. The checksum is **always**
            written to the fits header when the file is created.
        loggers (:obj:`list`, optional):
            List of `logging.Logger`_ objects to log progress;
            ignored if quiet=True. Logging is done using
            :func:`mangadap.util.log.log_output`. If None, no logging
            is performed and output is just written to ``stdout``.
        quiet (:obj:`bool`, optional):
            Suppress all terminal and logging output.

    Attributes:
        bitmask (class:`mangadap.util.bitmask.BitMask`):
            Object used to toggle mask values; see
            :func:`TemplateLibraryBitMask`.
        library (:class:`TemplateLibraryDef`):
            Parameter set required to read and prepare the library.
        file_list (:obj:`list`):
            The list of files found using `glob.glob`_ and
            :attr:`file_search`.
        ntpl (:obj:`int`):
            Number of template spectra in the library
        cube (:class:`mangadap.datacube.datacube.DataCube`):
            The datacube to be analyzed using the template library.
        sres (:class:`mangadap.util.resolution.SpectralResolution`):
            The object is used simply to access the spectral resolution
            and associated wavelength coordinate vector needed when
            matching the spectral resolution of the template library;
            this is used in place of the attributes in any provided DRP
            file object.
        velocity_offset (:obj:`float`):
            Velocity offset to use when matching the spectral resolution
            between the template library and the galaxy spectra.
        spectral_step (:obj:`float`):
            Target logarithmic (``log10_sampling=True``) or linear
            (``log10_sampling=False``) step in wavelength for the
            template library.
        log10_sampling (:obj:`bool`):
            Flag that the processed template library is logarithmically
            sampled in wavelength.
        directory_path (:obj:`str`):
            The exact path to the processed template library file.
        output_file (:obj:`str`):
            The name of the file containing (to contain) the prepared
            template library output file.  The file should be found at
            :attr:`directory_path`/:attr:`output_file`.
        processed (:obj:`bool`):
            Flag that the template library has been prepared for use in
            the DAP.
        hardcopy (:obj:`bool`):
            Flag to keep a hardcopy of the processed template library.
        symlink_dir (:obj:`str`):
            Symlink created to the file in this directory
        hdu (`astropy.io.fits.HDUList`_):
            HDUList read from the DAP file

    """
    # Class attribute
    supported_libraries = ['BC03', 'BPASS', 'M11ELODIE', 'M11MARCS', 'M11MILES', 'M11STELIB',   
                           'M11STELIBZSOL', 'MASTARHC', 'MASTARHC2', 'MASTARHC2-NOISM',
                           'MASTARSSP', 'MASTARSSP-NOISM', 'MILES',
                           'MILESAVG', 'MILESHC', 'MILESTHIN', 'MIUSCAT', 'MIUSCATTHIN', 'STELIB']
    """Provides the keywords of the supported libraries."""

    def __init__(self, library, cube=None, match_resolution=True, velscale_ratio=None, sres=None,
                 velocity_offset=0.0, min_sig_pix=0.0, no_offset=True, spectral_step=None,
                 log=True, wavelength_range=None, renormalize=True, output_path=None,
                 output_file=None, read=True, process=True, hardcopy=False, symlink_dir=None,
                 overwrite=False, checksum=False, loggers=None, quiet=False):

        self.loggers = loggers
        self.quiet = quiet

        # Define the TemplateLibraryBitMask object
        self.bitmask = TemplateLibraryBitMask()

        # Define the properties needed to modify the spectral resolution
        self.sres = None
        self.velocity_offset = None
        self.min_sig_pix = None
        self.no_offset = None

        # Define the target spectral sampling properties
        self.spectral_step = None
        self.log10_sampling = None

        # Define the library
        self.library = library if isinstance(library, TemplateLibraryDef) \
                            else TemplateLibraryDef.from_key(library)
        self.file_list = None
        self.ntpl = None

        # Define the processed file and flag, and the HDUList used to
        # keep the data
        self.directory_path = None
        self.output_file = None
        self.processed = False
        self.hardcopy = hardcopy
        self.symlink_dir = None
        self.hdu = None
        self.checksum = checksum

        # Do not read the library
        if not read:
            if not self.quiet:
                log_output(self.loggers, 1, logging.INFO, 'Nothing read, by request.')
            return
       
        # Do not process the library
        if not process:
            if not self.quiet:
                log_output(self.loggers, 1, logging.INFO, 'Reading raw library without processing.')
            self._read_raw()
            return

        # Read and process the library
        self.process_template_library(cube=cube, match_resolution=match_resolution,
                                      velscale_ratio=velscale_ratio, sres=sres,
                                      velocity_offset=velocity_offset, min_sig_pix=min_sig_pix,
                                      no_offset=no_offset, spectral_step=spectral_step, log=log,
                                      wavelength_range=wavelength_range, renormalize=renormalize,
                                      output_path=output_path, output_file=output_file,
                                      hardcopy=hardcopy, symlink_dir=symlink_dir,
                                      overwrite=overwrite, loggers=loggers, quiet=quiet)

    def __getitem__(self, key):
        return self.hdu[key]

    @staticmethod
    def default_paths(library_key, cube=None, output_path=None, output_file=None):
        """
        Set the default directory and file name for the output file.

        Args:
            library_key (:obj:`str`):
                Unique library keyword.
            cube (:class:`mangadap.datacube.datacube.DataCube`, optional):
                Datacube to analyze.  The output is specific to the cube if the
                spectral resolution and spectral sampling are matched to the
                cube data.  If None, the output file name is simply the library
                key.
            output_path (:obj:`str`, `Path`_, optional):
                The path for the output file.  If None, the current working
                directory is used.
            output_file (:obj:`str`, optional):
                The name of the output file. The full path of the
                output file will be :attr:`directory_path`/:attr:`output_file`.
                If None, the default is to combine ``cube.output_root`` and the
                library key.

        Returns:
            :obj:`tuple`: Returns a `Path`_ with the output directory and a
            :obj:`str` with the output file name.
        """
        directory_path = Path('.').resolve() if output_path is None \
                                else Path(output_path).resolve()
        root = f'{library_key}' if cube is None else f'{cube.output_root}-{library_key}'
        _output_file = f'{root}.fits.gz' if output_file is None else output_file
        return directory_path, _output_file
    
    def _get_nchannels(self):
        """
        Get the maximum number of wavelength channels needed to store
        all the template spectra in a single array.

        .. todo::
            - What happens if the spectrum has an empty primary
              extension?

        Returns:
            :obj:`int`: Maximum number of pixels used by the listed
            spectra.

        Raises:
            ValueError:
                Raised if the input template spectra are not
                one-dimensional.
        """
        max_npix = 0
        for f in self.file_list:
            if fits.getval(str(f), 'NAXIS') != 1:
                raise ValueError('{0} is not one dimensional!'.format(f))
            npix = fits.getval(f, 'NAXIS1')
            if max_npix < npix:
                max_npix = npix
        return max_npix

    def _build_raw_hdu(self, npix):
        r"""
        Build the "raw" template library arrays.

        This simply reads the provided list of fits files and puts
        them into arrays of size :math:`N_{\rm tpl} \times N_{\rm
        pix}`.

        This will *force* reading of the data, even if the :attr:`hdu`
        is already initialized.

        The :attr:`hdu` will contain the appropriate extensions, but
        it is important to note that the wavelength vectors will
        **not** necessarily be the same. That is, reading of the raw
        template spectra can accommodate spectra that have different
        wavelength coordinates. Any pixels that have no data are
        masked using the 'NO_DATA' bitmask flag; see
        :func:`TemplateLibraryBitMask`.

        The spectral resolution is set using the ``fwhm`` or
        ``sres_ext`` values in the library parameters.

        .. warning::

            Currently no errors are saved because none are expected for
            the template libraries.

        Args:
            npix (:obj:`int`):
                Number of spectral channels for the output arrays
        """
        if self.hdu is not None:
            if not self.quiet:
                log_output(self.loggers, 1, logging.INFO, 'Closing existing HDUList.')
            self.hdu.close()
            self.hdu = None

        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, 'Attempting to build raw data ...')
        # Allocate the vectors
        wave = numpy.zeros((self.ntpl, npix), dtype=float)
        flux = numpy.zeros((self.ntpl, npix), dtype=float)
        sres = numpy.zeros((self.ntpl, npix), dtype=float)
        mask = numpy.zeros((self.ntpl, npix), dtype=self.bitmask.minimum_dtype())
        soff = numpy.zeros(self.ntpl, dtype=float)

        # Read and save each spectrum and mask the unobserved
        # wavelengths
        for i in range(0,self.ntpl):
            if self.library['sres_ext'] is None:
                wave_, flux_ = read_template_spectrum(str(self.file_list[i]),
                                                      log10=self.library['log10'])
                # TODO: Are the copies here needed?
                wave[i,:wave_.size] = numpy.copy(wave_)
                flux[i,:wave_.size] = numpy.copy(flux_)
                sres = wave/self.library['fwhm']
            else:
                wave_, flux_, sres_ = read_template_spectrum(str(self.file_list[i]),
                                                             sres_ext=self.library['sres_ext'],
                                                             log10=self.library['log10'])
                wave[i,:wave_.size] = numpy.copy(wave_)
                flux[i,:wave_.size] = numpy.copy(flux_)
                sres[i,:wave_.size] = numpy.copy(sres_)

            if wave_.size != npix:
                mask[i,wave_.size:] = self.bitmask.turn_on(mask[i,wave_.size:],'NO_DATA')
            if self.library['lower_flux_limit'] is not None:
                indx = numpy.invert( flux_ > self.library['lower_flux_limit'] )
                mask[i,indx] = self.bitmask.turn_on(mask[i,indx], 'FLUX_INVALID')
            if self.library['wave_limit'] is not None \
                    and self.library['wave_limit'][0] is not None:
                indx = wave[i,:].ravel() < self.library['wave_limit'][0]
                mask[i,indx] = self.bitmask.turn_on(mask[i,indx], 'WAVE_INVALID')
            if self.library['wave_limit'] is not None \
                    and self.library['wave_limit'][1] is not None:
                indx = wave[i,:].ravel() > self.library['wave_limit'][1]
                mask[i,indx] = self.bitmask.turn_on(mask[i,indx], 'WAVE_INVALID')

        # (Re)Set the HDUList object
        self._reset_hdu(wave, flux, mask, sres, soff)

        # Add some keywords to the header
        self.hdu[0].header['TPLPROC'] = (0, 'Flag that library has been processed')
        self.processed = False
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, '... done')

    def _read_raw(self):
        """
        Read the 'raw' versions of the template library.

        The "raw" library is the data before its sampling or
        resolution have been altered.

        This is a simple wrapper for
        :func:`~mangadap.proc.templatelibrary.TemplateLibraryDef._get_file_list`,
        :func:`_get_channels` and :func:`_build_raw_hdu`.
        """
        self.file_list = self.library.get_file_list()
        self.ntpl = len(self.file_list)
        npix = self._get_nchannels()
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO,
                            'Found {0} {1} templates'.format(self.ntpl, self.library['key']))
            log_output(self.loggers, 1, logging.INFO,
                            'Maximum number of wavelength channels: {0}'.format(npix))
        self._build_raw_hdu(npix)

    def _reset_hdu(self, wave, flux, mask, sres, soff):
        r"""
        (Re)Set :attr:`hdu` to a new HDUList object using the input
        arrays.  Also sets the header items indicating the version of
        the template reader and the keyword for the library.

        .. warning::

            No checking is done concerning the size of each array!

        Args:
            wave (`numpy.ndarray`_):
                Array with the wavelengths of each pixel.
            flux (`numpy.ndarray`_):
                Array with the flux in each pixel.
            mask (`numpy.ndarray`_):
                Bitmask values for each pixel.
            sres (`numpy.ndarray`_):
                Spectral resolution, :math:`R=\lambda/\delta\lambda`,
                at each pixel.
            soff (`numpy.ndarray`_):
                The spectral resolution offset for each spectrum (see
                :func:`mangadap.util.resolution.GaussianKernelDifference`).
        """
        if self.hdu is not None:
            self.hdu.close()
            self.hdu = None
        self.hdu = fits.HDUList([ fits.PrimaryHDU(),
                                  fits.ImageHDU(wave, name='WAVE'),
                                  fits.ImageHDU(flux, name='FLUX'),
                                  fits.ImageHDU(mask, name='MASK'),
                                  fits.ImageHDU(sres, name='SPECRES'),
                                  fits.ImageHDU(soff, name='SIGOFF')])
        self.hdu['PRIMARY'].header['LIBKEY'] = (self.library['key'], 'Library identifier')

    def _wavelength_range(self, flag=None):
        """
        Return the valid wavelength range for each spectrum based on the
        first and last unmasked pixel; interspersed masked regions are
        not considered.

        Args:
            flag (:obj:`str`, :obj:`list`, optional):
                Flags to consider when determining the wavelength
                range; see
                :func:`mangadap.util.bitmasks.BitMask.flagged`.

        Returns:
            `numpy.ndarray`_: Two-element vector with wavelengths of
            the first and last valid pixels.
        """
        indx = numpy.where(numpy.invert(self.bitmask.flagged(self.hdu['MASK'].data, flag=flag)))
        return numpy.array([numpy.amin(self.hdu['WAVE'].data[indx]),
                            numpy.amax(self.hdu['WAVE'].data[indx])]).astype(float)

    def _minimum_sampling(self):
        """
        Return the minimum sampling of the available wavelength vectors.
        
        Returns:
            :obj:`float`: Minimum sampling of all (will just be one
            if the library has been processed) wavelength vectors.
        """
        if self.processed:
            return self.spectral_step

        nspec = self.hdu['WAVE'].data.shape[0]
        minimum_step = sampling.spectral_coordinate_step(self.hdu['WAVE'].data[0,:].ravel(),
                                                         log=self.library['log10'])
        for i in range(1,nspec):
            step = sampling.spectral_coordinate_step(self.hdu['WAVE'].data[0,:].ravel(),
                                                     log=self.library['log10'])
            if minimum_step > step:
                minimum_step = step
        return minimum_step

    def _rebin_masked(self, i, flag, fullRange, rmsk_lim=0.5):
        """
        Determine the mask value to adopt for a rebinned spectrum by
        rebinning the mask pixels, setting a masked pixel to unity,
        and an unmasked pixel to zero. After rebinning, any pixel
        with a value of larger than ``rmsk_lim`` is masked; otherwise
        it is left unmasked.

        Although the function can be used with multiple flags, its
        intended use is to determine which pixels should be masked
        with a specific flag.

        Args:
            i (:obj:`int`):
                Index of the spectrum to be rebinned.
            flag (:obj:`str`, :obj:`list`):
                Flags to consider when determining which pixels to
                mask; see
                :func:`mangadap.util.bitmasks.BitMask.flagged`.
            fullRange (`numpy.ndarray`_):
                Two-element array with the wavelength range for the
                rebinned spectrum.
            rmsk_lim (:obj:`float`):
                Limit of the rebinned mask value that is allowed
                before considering the pixel as masked.

        Returns:
            `numpy.ndarray`_: Boolean array indicating which pixels
            in the rebinned spectrum should be masked.
        """
        unity = numpy.ma.MaskedArray(numpy.ones(self.hdu['WAVE'].data.shape[1], dtype=float),
                                     mask=self.bitmask.flagged(self.hdu['MASK'].data[i,:],
                                                               flag=flag))
        r = sampling.Resample(unity, x=self.hdu['WAVE'].data[i,:], inLog=self.library['log10'],
                              newRange=fullRange, newLog=self.log10_sampling,
                              newdx=self.spectral_step)
        return r.outf < rmsk_lim

    # TODO: Having this return the redshift is a bit weird.
    def _modify_spectral_resolution(self):
        """
        Modify the spectral resolution to match the provided
        :attr:`sres`.

        Returns:
            :obj:`float`:
                Redshift used when matching the spectral resolution.
        """
        # Calculate to use as an offset of the match to the spectral
        # resolution. Used to better match the spectral resolution to
        # at the *observed* wavelengths of the object spectrum to which
        # the TemplateLibrary will be fit.
        redshift = self.velocity_offset/astropy.constants.c.to('km/s').value

        # Mask wavelengths where the spectral resolution will have to be
        # extrapolated.
        sres_wave = self.sres.wave()
        wavelim = numpy.array([sres_wave[0]/(1.+redshift), sres_wave[-1]/(1.+redshift)])
        self.hdu = HDUList_mask_wavelengths(self.hdu, self.bitmask, 'SPECRES_EXTRAP', wavelim,
                                            invert=True)

        # Match the resolution of the templates.  ivar is returned, but
        # is always None because ivar is not provided to
        # match_spectral_resolution
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, 'Modifying spectral resolution ... ')

#        flux_inp = self.hdu['FLUX'].data.copy()

        self.hdu['FLUX'].data, self.hdu['SPECRES'].data, self.hdu['SIGOFF'].data, res_mask, \
                ivar = match_spectral_resolution(self.hdu['WAVE'].data, self.hdu['FLUX'].data,
                                                 self.hdu['SPECRES'].data, sres_wave/(1.+redshift),
                                                 self.sres.sres(), min_sig_pix=self.min_sig_pix,
                                                 no_offset=self.no_offset,
                                                 log10=self.library['log10'], new_log10=True,
                                                 quiet=self.quiet)
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, '... done')

        # Mask any pixels where the template resolution was too low to
        # match to the galaxy resolution
        self.hdu['MASK'].data[res_mask == 1] = \
            self.bitmask.turn_on(self.hdu['MASK'].data[res_mask == 1], 'SPECRES_LOW')

        # Mask any pixels where the template flux is 0 or below
        indx = numpy.invert(self.hdu['FLUX'].data > 0)
        if numpy.any(indx):
            self.hdu['MASK'].data[indx] = self.bitmask.turn_on(self.hdu['MASK'].data[indx],
                                                               'SPECRES_NOFLUX')

        # Function returns the redshift used for the spectral
        # resolution matching. TODO: This should be input instead!
        return redshift


    def _process_library(self, wavelength_range=None, renormalize=True):
        """
        Process the template library for use in analyzing object
        spectra.   See :func:`process_template_library`.

        .. todo::

            - Add wavelength coordinate WCS information to the
              appropriate extension headers.
            - Include 'step' as an argument
            - Include limit on covering fraction used to mask data as an
              argument

        Args:
            wavelength_range (array-like, optional):
                A two-element object with the lower and upper limits
                for the wavelength range of the library.
            renormalize (:obj:`bool`, optional):
                Renormalize the flux of the full library to unity.
                Relative fluxes of spectra in the library are
                maintained.
        """
        # Convert to vacuum wavelengths
        # TODO: This conversion to vacuum wavelengths is non-linear.
        # This is accounted for in the resampling of the spectra (the
        # pixel size can be non-uniform), but *not* the resolution
        # matching.  This is okay for DR15 because we never resolution
        # match the data, but it will matter when we start to do that by
        # switching templates between the stellar and gas kinematics
        # fit.
        if not self.library['in_vacuum']:
            self.hdu['WAVE'].data = airtovac(self.hdu['WAVE'].data)

        #---------------------------------------------------------------
        # Modify the spectral resolution to a target function, if one
        # has been provided.
        redshift = 0.0 if self.sres is None else self._modify_spectral_resolution()

        #---------------------------------------------------------------
        # Resample the templates to a logarithmic binning step.
        #
        # Typically this will be run to match the sampling of the
        # template library to that of a set of galaxy spectra to be fit.
        # However, the sampling can be left alone and this part will
        # force all the spectra to have the same wavelength range.
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, 'Matching sampling ... ')

        # Number of angstroms per pixel
        ang_per_pix = numpy.array([sampling.angstroms_per_pixel(w, log=self.library['log10'],
                                                                base=10.,
                                                                regular=self.library['in_vacuum'])
                                        for w in self.hdu['WAVE'].data])

        # TODO: I think this kludge was necessary because of a bug in
        # angstroms_per_pixel.  If regular was True and log was False,
        # it only returns a single integer, otherwise it returned a
        # vector.  I've changed angstroms_per_pixel, but have left the
        # kludge below for now for testing.

        #KHRR KLUDGE
        if len(ang_per_pix.shape) == 1:
            warnings.warn('Executing KLUDGE!  Bug fix insufficient.')
            tmp_ang_per_pix = numpy.empty(self.hdu['SPECRES'].data.shape)
            for ind in range(len(ang_per_pix)):
                tmp_ang_per_pix[ind,:] = ang_per_pix[ind]
            ang_per_pix = tmp_ang_per_pix

        # Number of pixels per resolution element
        pix_per_fwhm = numpy.ma.divide(self.hdu['WAVE'].data,
                                       self.hdu['SPECRES'].data * ang_per_pix)

        # TODO: Put in a check of whether or not it needs to do the
        # resampling!

        # The raw spectra are allowed to have wavelength ranges that
        # differ.  First, determine the wavelength range that encloses
        # all spectra.  Only ignore pixels that were flagged as having
        # no data.
        no_data_flags = ['NO_DATA', 'SPECRES_NOFLUX']
        fullRange = self._wavelength_range(flag=no_data_flags) if wavelength_range is None \
                            else numpy.array(wavelength_range).astype(float)

        # Get the spectral step if it hasn't been set yet
        if self.spectral_step is None:
            self.spectral_step = self._minimum_sampling()

        # Get the number of pixels needed
        npix, _fullRange = sampling.grid_npix(rng=fullRange, dx=self.spectral_step,
                                              log=self.log10_sampling)

        # Any pixels without data after resampling are given a value
        # that is the minimum flux - 100 so that they can be easily
        # identified afterward.  The minimum flux is:
        min_flux = numpy.amin(self.hdu['FLUX'].data.ravel())
        # the observed pixels are
        no_data = self.bitmask.flagged(self.hdu['MASK'].data, flag='NO_DATA')

        # Now resample the spectra.  First allocate the arrays
        flux = numpy.zeros((self.ntpl, npix), dtype=float)
        sres = numpy.zeros((self.ntpl, npix), dtype=float)
        mask = numpy.zeros((self.ntpl, npix), dtype=self.bitmask.minimum_dtype())

        for i in range(self.ntpl):
            # Rebin the observed wavelength range
            r = sampling.Resample(self.hdu['FLUX'].data[i,:], mask=no_data[i,:],
                                  x=self.hdu['WAVE'].data[i,:], inLog=self.library['log10'],
                                  newRange=fullRange, newLog=self.log10_sampling,
                                  newdx=self.spectral_step, step=False)

            # Save the result
            wave = r.outx
            flux[i,:] = r.outy
            indx = r.outf < 0.9
            # TODO: Set the flux to 0?
            flux[i,indx] = 0.0
            mask[i,indx] = self.bitmask.turn_on(mask[i,indx], 'NO_DATA')

            # Resample the spectral resolution by simple interpolation.
            sres[i,:] = interpolate.interp1d(self.hdu['WAVE'].data[i,:],
                                             self.hdu['SPECRES'].data[i,:], assume_sorted=True,
                                             fill_value='extrapolate')(wave)
            sres[i,indx] = 0.0

            # Recalculate the number of pixels per fwhm
            # - Number of angstroms per pixel
            _ang_per_pix = sampling.angstroms_per_pixel(wave, log=self.log10_sampling, base=10.,
                                                        regular=True)
            # - Number of pixels per resolution element
            _pix_per_fwhm = numpy.ma.divide(wave, sres[i,:] * _ang_per_pix)

            indx = _pix_per_fwhm < 2
            if numpy.sum(indx) > 0:
                warnings.warn('Resampling results in resolution below the two pixel limit!')
                sres[i,indx] = wave[indx]/2./_ang_per_pix[indx]
                mask[i,indx] = self.bitmask.turn_on(mask[i,indx], 'SPECRES_2PIXEL')
            #-----------------------------------------------------------

            # Rebin the masks, bit-by-bit:
            # Pixels outside the wavelength limits
            indx = self._rebin_masked(i, 'WAVE_INVALID', fullRange, rmsk_lim=0.1)
            mask[i,indx] = self.bitmask.turn_on(mask[i,indx], 'WAVE_INVALID')
            # Pixels below the flux limit
            indx = self._rebin_masked(i, 'FLUX_INVALID', fullRange, rmsk_lim=0.1)
            mask[i,indx] = self.bitmask.turn_on(mask[i,indx], 'FLUX_INVALID')
            # Pixels that required an extrapolation of the spectral
            # resolution
            indx = self._rebin_masked(i, 'SPECRES_EXTRAP', fullRange, rmsk_lim=0.1)
            mask[i,indx] = self.bitmask.turn_on(mask[i,indx], 'SPECRES_EXTRAP')
            # Pixels that had a spectral resolution that was too low to
            # match the galaxy resolution
            indx = self._rebin_masked(i, 'SPECRES_LOW', fullRange, rmsk_lim=0.1)
            mask[i,indx] = self.bitmask.turn_on(mask[i,indx], 'SPECRES_LOW')
            # Pixels that no flux because of the wavelength limits set
            # by the resolution matching.
            indx = self._rebin_masked(i, 'SPECRES_NOFLUX', fullRange, rmsk_lim=0.1)
            mask[i,indx] = self.bitmask.turn_on(mask[i,indx], 'SPECRES_NOFLUX')
            
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, '... done')
            log_output(self.loggers, 1, logging.INFO,
                       'After resampling (target): {0}'.format(self.spectral_step))
            log_output(self.loggers, 1, logging.INFO,
                       'After resampling (calculated): {0}'.format(
                            sampling.spectral_coordinate_step(wave, log=True)))

        # Finally, we need to (again) trim the wavelength range to
        # exclude pixels with invalid or empty flux.
        edges = numpy.ma.notmasked_edges(numpy.ma.MaskedArray(numpy.ones(mask.shape),
                            mask=self.bitmask.flagged(mask,
                            flag=['NO_DATA', 'WAVE_INVALID', 'FLUX_INVALID', 'SPECRES_NOFLUX'])),
                                         axis=1)
        s = numpy.amax(edges[0][1])
        e = numpy.amin(edges[1][1])+1
        wave = wave[s:e]
        flux = flux[:,s:e]
        mask = mask[:,s:e]
        sres = sres[:,s:e]

        # Normalize the templates to the mean flux value after excluding
        # any flagged pixels.
        if renormalize:
            indx = numpy.invert(self.bitmask.flagged(mask, flag=['NO_DATA', 'WAVE_INVALID', 
                                                                 'FLUX_INVALID']))
            if numpy.sum(indx) == 0:
                if not self.quiet:
                    warnings.warn('All pixels masked.  Unable to renormalize TemplateLibrary.')
                flux_norm = 1.0
            else:
                flux_norm = numpy.mean(flux[numpy.where(indx)])
            flux /= flux_norm
        else:
            flux_norm = 1.0

        # Reset the HDUList object
        self._reset_hdu(wave, flux, mask, sres, self.hdu['SIGOFF'].data)

        # Update the header with the redshift offset, the flux
        # normalization, and flag the data as having been prepared for
        # fitting the DRP data
        self.hdu[0].header['ZGUESS'] = (redshift,'Redshift used when matching spectral resolution')
        self.hdu[0].header['FLXNRM'] = (flux_norm,'Flux normalization')
        self.hdu[0].header['TPLPROC'] = 1
        self.processed = True

    def file_name(self):
        """Return the name of the processed file."""
        return self.output_file

    def file_path(self):
        """Return the full path to the processed file."""
        if self.directory_path is None or self.output_file is None:
            return None
        return self.directory_path / self.output_file

    def read_raw_template_library(self, library=None):
        """
        Read the identified template library.  If all the arguments are
        the default, the preset attributes from the initialization of
        the object are used.

        Args:
            library (:obj:`str`, :class:`TemplateLibraryDef`, optional):
                Keyword selecting the library to use or the full set of defining
                parameters.
        """
        if library is not None:
            self.library = library if isinstance(library, TemplateLibraryDef) \
                                else TemplateLibraryDef.from_key(library)
        self._read_raw()

    def process_template_library(self, library=None, cube=None, match_resolution=True,
                                 velscale_ratio=None, sres=None, velocity_offset=0.0,
                                 min_sig_pix=0.0, no_offset=True, spectral_step=None, log=True,
                                 wavelength_range=None, renormalize=True, output_path=None,
                                 output_file=None, hardcopy=False, symlink_dir=None,
                                 overwrite=False, loggers=None, quiet=False):
        r"""
        Process the template library for use in analyzing spectra.

        .. todo::
            - Docs need to be updated
            - type checking
            - If a DRP file is provided, the processing to a logarithmic
              binning is done by default (log=True).  But linearly
              sampled DRP data are available, so need to have
              :class:`mangadap.util.drpfits.DRPFits` return the spectral
              sampling type.
            - Documentation needs updating!
        
        Primary steps are to:

            - Read the raw 1D fits files; see :func:`_read_raw`.

            - Convert the wavelengths to vacuum, if necessary; see
              `pydl.goddard.astro.airtovac`_.

            - Mask wavelengths outside the rest wavelength range of the
              datacube spectra, due to need to extrapolate these values; see
              :func:`mangadap.util.bitmasks.HDUList_mask_wavelengths`.

            - Match the spectral resolution of the template to that of
              the datacube (if ``match_resolution`` is True); see
              :func:`mangadap.util.resolution.match_spectral_resolution`.

            - Mask the template pixels where the spectral resolution was
              too low to match to datacube; see
              :func:`mangadap.util.bitmasks.BitMask.turn_on`.

            - Force a common wavelength range and sampling for all
              templates, where the sampling is forced to match the
              sampling of the DRP spectra; see
              :func:`mangadap.util.sampling.grid_pix`. The masks are
              appropriately resampled as well; see
              :func:`_rebin_masked`.

        .. warning::

            The routine **does not** check that that an existing
            processed file or the existing object has been processed
            using the same datacube, velocity_offset, velscale, or
            sres input. If unsure, use overwrite=True.
        
        Args:
            library (:obj:`str`, :class:`TemplateLibraryDef`, optional):
                Keyword selecting the library to use or the full set of defining
                parameters.
            cube (:class:`mangadap.datacube.datacube.DataCube`, optional):
                The datacube to be analyzed using the template library.
            match_resolution (:obj:`bool`, optional):
                Match the spectral resolution of the template library to
                the resolution provided by the ``cube``; the latter must
                be provided for this argument to have any use.
            velscale_ratio (:obj:`int`, optional):
                Resample the template spectra such that the ratio of the
                pixel scale in the provided ``cube`` is this many times
                larger than the pixels in the resampled template
                spectrum.
            sres (:class:`mangadap.util.resolution.SpectralResolution`, optional):
                The object is used simply to access the spectral
                resolution and associated wavelength coordinate vector
                needed when matching the spectral resolution of the
                template library; this is used in lieu of a datacube.
            velocity_offset (:obj:`float`, optional):
                Velocity offset to use when matching the spectral
                resolution between the template library and the galaxy
                spectra.
            min_sig_pix (:obj:`float`, optional):
                Minimum value of the LSF standard deviation in pixels
                allowed before assuming the kernel is a Delta function.
                See
                :func:`mangadap.util.resolution.match_spectral_resolution`.
            no_offset (:obj:`bool`, optional):
                Prohibit the spectral resolution matching to introduce an
                offset in spectral resolution required to match the
                relative shape of the resolution vectors. See
                :func:`mangadap.util.resolution.match_spectral_resolution`
                and
                :func:`mangadap.util.resolution.SpectralResolution.GaussianKernelDifference`.
            spectral_step (:obj:`float`, optional):
                Logarithmic (``log=True``) or linear (``log=False``) step
                in wavelength for the template library.
            log (:obj:`bool`, optional):
                Flag to force the library to be logarithmically sampled in
                wavelength.
            wavelength_range (array-like, optional):
                Force the template library to cover this spectral range.
                Unobserved spectral regions will be masked.
            renormalize (:obj:`bool`, optional):
                After processing, renormalize the flux of the full
                library to unity. The relative fluxes of the templates
                are maintained.
            output_path (:obj:`str`, `Path`_, optional):
                The path for the output file.  If None, the current working
                directory is used.
            output_file (:obj:`str`, optional):
                The name of the output file. The full path of the
                output file will be :attr:`directory_path`/:attr:`output_file`.
                If None, the default is to combine ``cube.output_root`` and the
                library key; if ``cube`` is None, it is simply the library key.
            hardcopy (:obj:`bool`, optional):
                Flag to keep a hardcopy of the processed template library.
            symlink_dir (:obj:`str`, optional):
                Create a symlink to the file in this directory. If None,
                no symbolic link is created.
            overwrite (:obj:`bool`, optional):
                Overwrite any saved, processed template library file.
            loggers (:obj:`list`, optional):
                List of `logging.Logger`_ objects to log progress;
                ignored if quiet=True. Logging is done using
                :func:`mangadap.util.log.log_output`. If None, no logging
                is performed and output is just written to ``stdout``.
            quiet (:obj:`bool`, optional):
                Suppress all terminal and logging output.

        Raises:
            ValueError:
                If the sampling is undefined, or if a hardcopy is
                requested but the output file could not be defined.
        """
        # Initialize the reporting
        if loggers is not None:
            self.loggers = loggers
        self.quiet = quiet

        # TODO: Not sure we should allow the library key to change in this
        # way...
        if library is not None:
            self.library = library if isinstance(library, TemplateLibraryDef) \
                                else TemplateLibraryDef.from_key(library)

        # Ignore velscale_ratio if it is set to unity
        if velscale_ratio is not None and velscale_ratio == 1:
            velscale_ratio = None

        if cube is not None:
            # Use the provided datacube to set the sampling
            self.sres = None
            self.log10_sampling = cube.log
            self.spectral_step = sampling.spectral_coordinate_step(cube.wave, log=cube.log)

            # Use the datacube to set the spectral resolution, if
            # resolution matching was requested. The fiducial spectral
            # resolution is set to be the median spectral resolution of
            # all the spectra that are not fully masked.
            if match_resolution:
                # Copy the spectral resolution to a masked array
                self.sres = cube.copy_to_masked_array(attr='sres')
                # Only use the spectra that are not fully masked
                indx = numpy.sum(numpy.invert(numpy.ma.getmaskarray(self.sres)), axis=1) > 0
                # Get the median resolution
                self.sres = numpy.median(self.sres.data[indx], axis=0)
                # Instantiate the resolution object
                self.sres = SpectralResolution(cube.wave, self.sres, log10=self.log10_sampling)
        else:
            # Use the provided input values
            self.sres = sres
            self.log10_sampling = log
            self.spectral_step = spectral_step

        self.velocity_offset = velocity_offset
        self.min_sig_pix = min_sig_pix
        self.no_offset = no_offset

        # Adjust for the velocity scale ratio between the template and
        # object data to be fit
        if not self.log10_sampling and velscale_ratio is not None:
            raise ValueError('velscale_ratio only valid with logarithmically sampled spectra.')
        if velscale_ratio is not None:
            self.spectral_step /= velscale_ratio

        # (Re)Set the output paths
        self.directory_path, self.output_file \
                = TemplateLibrary.default_paths(self.library['key'], cube=cube,
                                                output_path=output_path, output_file=output_file)
        self.hardcopy = hardcopy
        self.symlink_dir = symlink_dir

        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, '-'*50)
            log_output(self.loggers, 1, logging.INFO,
                       f'Template library output path: {self.directory_path}')
            log_output(self.loggers, 1, logging.INFO,
                       f'Template library output file: {self.output_file}')

        # Check that the path for or to the file is defined
        ofile = self.file_path()
        if self.hardcopy and ofile is None:
            raise ValueError('File path for output file is undefined!')

        # Read and use a pre-existing file
        if self.hardcopy and ofile.exists() and not overwrite:
            if not self.quiet:
                log_output(self.loggers, 1, logging.INFO, 'Reading existing file')
            self.hdu = DAPFitsUtil.read(ofile, checksum=self.checksum)
            self.file_list = self.library.get_file_list()
            self.ntpl = self.hdu['FLUX'].data.shape[0]
            self.processed = True
            # Make sure the symlink exists
            if self.symlink_dir is not None:
                create_symlink(str(ofile), self.symlink_dir, overwrite=overwrite)
            if not self.quiet:
                log_output(self.loggers, 1, logging.INFO, '-'*50)
            return

        # Warn the user that the file will be overwritten
        if self.hardcopy and ofile.exists():
            if not self.quiet:
                warnings.warn(f'Overwriting existing file: {ofile}')

        # Read the raw data
        self._read_raw()

        # Process the library
        self._process_library(wavelength_range=wavelength_range, renormalize=renormalize)

        # Write the fits file
        if self.hardcopy:
            DAPFitsUtil.write(self.hdu, str(ofile), overwrite=overwrite, checksum=True,
                              symlink_dir=self.symlink_dir, loggers=self.loggers, quiet=self.quiet)
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, '-'*50)

    def single_spec_to_fits(self, i, ofile, overwrite=True):
        """
        Write one of the template spectra to a 1D fits file.

        Output is a one-dimensional fits file with a single extension
        containing the flux measurements and header keywords needed
        to determine the wavelength of each pixel.

        Args:
            i (:obj:`int`):
                Index of the spectrum in the avaiable list to output.
            ofile (:obj:`str`):
                Name of the file to write
            overwrite (:obj:`bool`, optional):
                Overwrite any existing file.

        Raises:
            ValueError:
                Raised if the selected index is not available.
        """
        if i < 0 or i >= self.ntpl:
            raise ValueError('Invalid index ({0})!  Library contains {1} spectra.'.format(i,
                                                                                        self.ntpl))
        wave = self.hdu['WAVE'].data if self.processed else self.hdu['WAVE'].data[i,:]
        log = self.log10_sampling if self.processed else self.library['log10']
        writefits_1dspec(ofile, numpy.log10(wave[0]) if log else wave[0], self.spectral_step,
                         self.hdu['FLUX'].data[i,:], overwrite=overwrite)


