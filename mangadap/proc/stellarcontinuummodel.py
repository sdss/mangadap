# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
A class hierarchy that performs the stellar-continuum fitting.

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

from astropy.io import fits
import astropy.constants

from ..par.parset import KeywordParSet, ParSet
from ..util.log import log_output
from ..util.fitsutil import DAPFitsUtil
from ..util.fileio import create_symlink
from ..util.sampling import spectral_coordinate_step, spectrum_velocity_scale
from ..util.resolution import SpectralResolution
from ..util.dapbitmask import DAPBitMask
from .spatiallybinnedspectra import SpatiallyBinnedSpectra
from .templatelibrary import TemplateLibrary
from .ppxffit import PPXFFitPar, PPXFFit
from .util import replace_with_data_from_nearest_coo


class StellarContinuumModelDef(KeywordParSet):
    """
    A class that holds the parameters necessary to perform the
    stellar-continuum fitting.

    The provided ``fitfunc`` must have the form::

        model_wave, model_flux, model_mask, model_par \
                = fitfunc(self, binned_spectra, par=None)

    where ``binned_spectra`` has type
    :class:`mangadap.proc.spatiallybinnedspetra.SpatiallyBinnedSpectra`,
    and the returned objects are the wavelength vector, the fitted
    model flux, a bitmask for the model, and the set of model
    parameters. For example, see
    :func:`mangadap.proc.ppxffit.PPXFFit.fit_SpatiallyBinnedSpectra`.

    The defined parameters are:

    .. include:: ../tables/stellarcontinuummodeldef.rst
    """
    def __init__(self, key='MILESHC', minimum_snr=1.0, fitpar=None, fitclass=None, fitfunc=None,
                 overwrite=False):

        # Use the signature to get the parameters and the default values
        sig = inspect.signature(self.__class__)
        pars = list(sig.parameters.keys())
        defaults = [sig.parameters[key].default for key in pars]

        # Remaining definitions done by hand
        in_fl = [int, float]
        par_opt = [ParSet, dict]

        values = [key, minimum_snr, fitpar, fitclass, fitfunc, overwrite ]
        dtypes = [str, in_fl, par_opt, None, None, bool]
        can_call = [False, False, False, False, True, False]
        descr = ['Keyword used to distinguish between different spatial binning schemes.',
                 'Minimum S/N of spectrum to fit',
                 'Any additional parameters, aside from the spectra themselves, required by ' \
                    'the fitting function.  If None, uses a default pPXF fit.',
                 'Instance of class object to use for the model fitting.  Needed in case ' \
                    'fitfunc is a non-static member function of a class.  If None, uses a ' \
                    'default pPXF fit.',
                 'The function that models the spectra.  If None, uses a default pPXF fit and ' \
                    'anything provided for fitpar and fitclass are ignored!',
                 'If the output file already exists, redo all the calculations and overwrite it.']

        super().__init__(pars, values=values, defaults=defaults, dtypes=dtypes, can_call=can_call,
                         descr=descr)
        self._validate()

    def _validate(self):
        """
        Validate the modeling method.
        """
        if self['fitfunc'] is None:
            warnings.warn('No fitting function defined.  Defaulting to use mangadap.proc.ppxffit.')
            self['fitpar'] = PPXFFitPar()
            self['fitclass'] = PPXFFit(StellarContinuumModelBitMask())
            self['fitfunc'] = self['fitclass'].fit_SpatiallyBinnedSpectra

    @classmethod
    def from_dict(cls, d):
        """
        Instantiate the object from a nested dictionary.

        The dictionary is expected to have keys for all the class keywords; see
        the class instantiation method.  It should also include the
        ``'fit_method'`` keyword used to select the type of fit to perform, and
        a nested dictionary selected by the keyword ``'fit'`` with the 
        parameters required by the fitting method.  If these are not provided,
        the default fitter, :class:`~mangadap.proc.ppxffit.PPXFFit` and
        parameters will be used.

        .. warning::

            Currently, the only accepted value for ``'fit_method'`` is
            ``'ppxf'``, and the ``'fit'`` dictionary must provide the parameters
            defined by :class:`~mangadap.proc.ppxffit.PPXFFitPar`.

        Args:
            d (:obj:`dict`):
                Dictionary used to instantiate the class.
        """
        # Copy over the primary keywords
        _d = {}
        for key in ['key', 'minimum_snr', 'overwrite']:
            if key in d.keys():
                _d[key] = d[key]

        # Set the spatial binning method
        if 'fit_method' not in d or d['fit_method'] in ['none', None]:
            # Use the default fitter and default fitting parameters!
            _d['fitpar'] = None
            _d['fitclass'] = None
            _d['fitfunc'] = None
        elif d['fit_method'] == 'ppxf':
            _d['fitpar'] = PPXFFitPar.from_dict(d['fit']) if 'fit' in d else PPXFFitPar()
            _d['fitclass'] = PPXFFit(StellarContinuumModelBitMask())
            _d['fitfunc'] = _d['fitclass'].fit_SpatiallyBinnedSpectra
        else:
            raise ValueError('Unrecognized fitting method.  Requires direct instantiation or '
                             'code alterations.')
        # Return the instantiation
        return super().from_dict(_d)


class StellarContinuumModelBitMask(DAPBitMask):
    r"""
    Derived class that specifies the mask bits for the stellar-continuum
    modeling. The maskbits defined are:

    .. include:: ../tables/stellarcontinuummodelbitmask.rst
    """
    cfg_root = 'stellar_continuum_model_bits'


class StellarContinuumModel:
    r"""
    Class that holds the stellar-continuum model results.

    .. todo::

        - Update attributes

    Args:
        method (:class:`~mangadap.proc.stellarcontinuummodel.StellarContinuumModelDef`):
            Object that defines the main parameters used for the stellar
            continuum fitting.
        binned_spectra (:class:`mangadap.proc.spatiallybinnedspectra.SpatiallyBinnedSpectra`):
            The binned spectra to fit. Must be able to access parent
            :class:`mangadap.datacube.datacube.DataCube` attribute
            (``cube``).
        guess_vel (:obj:`float`, `numpy.ndarray`_):
            A single or spectrum-dependent initial estimate of the
            redshift in km/s (:math:`cz`).
        guess_sig (:obj:`float`, `numpy.ndarray`_, optional):
            A single or spectrum-dependent initial estimate of the
            velocity dispersion in km/s. If None, default set to 100
            km/s.
        method_list (:obj:`list`, optional):
            List of :class:`StellarContinuumModelDef` objects that
            define one or more methods to use for the stellar
            continuum fitting. The default list is provided by the
            config files in the DAP source directory and compiled
            into this list using
            :func:`available_stellar_continuum_modeling_methods`.
        output_path (:obj:`str`, `Path`_, optional):
            The path for the output file.  If None, the current working
            directory is used.
        output_file (:obj:`str`, optional):
            The name of the output "reference" file. The full path of the output
            file will be :attr:`directory_path`/:attr:`output_file`.  If None,
            the default is to combine ``cube.output_root`` and the method keys.
            The order of the keys is the order of operations (rdxqa, binning,
            stellar continuum).  See
            :func:`~mangadap.proc.stellarcontinuummodel.StellarContinuumModel.default_paths`.
        hardcopy (:obj:`bool`, optional):
            Flag to write the HDUList attribute to disk. If False,
            the core `astropy.io.fits.HDUList`_ attribute
            (:attr:`hdu`) is only kept in memory and would have to be
            reconstructed.
        symlink_dir (:obj:`str`, optional):
            Create a symbolic link to the created file in the supplied
            directory.  If None, no symbolic link is created.
        tpl_hardcopy (:obj:`bool`, optional):
            Save the processed template library used during the fit.  If True,
            the output path and symlink directory for the template library are
            identical to the main reference file.
        overwrite (:obj:`bool`, optional):
            Overwrite any existing files. Default is to use any
            existing file instead of redoing the analysis and
            overwriting the existing output.
        checksum (:obj:`bool`, optional):
            Use the checksum in the fits header to confirm that the
            data has not been corrupted. The checksum is **always**
            written to the fits header when the file is created.
        loggers (:obj:`list`, optional):
            List of `logging.Logger`_ objects to log progress;
            ignored if quiet=True. Logging is done using
            :func:`mangadap.util.log.log_output`. Default is no
            logging.
        quiet (:obj:`bool`, optional):
            Suppress all terminal and logging output.

    Attributes:
        loggers (:obj:`list`):
            List of `logging.Logger`_ objects to log progress.
        quiet (:obj:`bool`):
            Suppress all terminal and logging output.

    """
    def __init__(self, method, binned_spectra, guess_vel, guess_sig=None, method_list=None,
                 output_path=None, output_file=None, hardcopy=True, symlink_dir=None,
                 tpl_hardcopy=False, overwrite=False, checksum=False, loggers=None, quiet=False):

        self.loggers = None
        self.quiet = False

        # Define the method properties
        self.method = method
        if not isinstance(self.method, StellarContinuumModelDef):
            raise TypeError('Method must have type StellarContinuumModelDef.')

        self.binned_spectra = None
        self.guess_vel = None
        self.guess_sig = None

        # Define the output directory and file
        self.directory_path = None      # Set in default_paths
        self.output_file = None
        self.hardcopy = None
        self.tpl_hardcopy = None
        self.symlink_dir = None

        # Define the bitmask
        self.bitmask = StellarContinuumModelBitMask()

        # Initialize the main class attributes
        self.hdu = None
        self.checksum = checksum
        self.shape = None
        self.spatial_shape = None
        self.nspec = None
        self.spatial_index = None
        self.spectral_arrays = None
        self._assign_spectral_arrays()
        self.image_arrays = None
        self._assign_image_arrays()
        self.nwave = None

        self.nmodels = None
        self.missing_models = None

        # TODO: Include covariance between measured properties

        # Run the assessments of the DRP file
        self.fit(binned_spectra, guess_vel, guess_sig=guess_sig, output_path=output_path,
                 output_file=output_file, hardcopy=hardcopy, symlink_dir=symlink_dir,
                 tpl_hardcopy=tpl_hardcopy, overwrite=overwrite, loggers=loggers, quiet=quiet)

    def __getitem__(self, key):
        return self.hdu[key]

    def get_template_library(self, output_path=None, output_file=None, hardcopy=None,
                             symlink_dir=None, velocity_offset=None, match_resolution=False,
                             resolution_fwhm=None):
        """
        Construct the template library for the continuum fitting.

        If the FWHM (in angstroms; see ``resolution_fwhm``) is not
        provided, this method is a simple wrapper that returns a
        :class:`mangadap.proc.templatelibrary.TemplateLibrary`.
        Otherwise, the provided FWHM is used to construct the
        spectral resolution for the template library.

        In either case, the sampling is set to match the spectra in
        :attr:`binned_spectra` (up to the velocity-scale ratio set by
        :attr:`method`).

        Args:
            output_path (:obj:`str`, `Path`_, optional):
                The path for the output file.  If None, set to 
                :attr:`directory_path`.  An output file is only produced if
                ``hardcopy`` is True.
            output_file (:obj:`str`, optional):
                The name of the output template library file.  If None and
                ``resolution_fwhm`` is None, uses the default template library
                file; see
                :func:`~mangadap.proc.templatelibrary.TemplateLibrary.default_paths`.
                If None and ``resolution_fwhm`` is provided, the default
                template file is also used, but the library key is changed to
                indicate the new FHWM of the output resolution.  An output file
                is only produced if ``hardcopy`` is True.
            hardcopy (:obj:`bool`, optional):
                Flag to keep a hardcopy of the processed template library.  If
                None, the flag is set by :attr:`tpl_hardcopy`.  If no hardcopy
                is to be produced, ``output_path`` and ``output_file`` are
                ignored.
            symlink_dir (:obj:`str`, optional):
                Create a symbolic link to the created template library file in
                this directory. If None, :attr:`symlink_dir` is used.  A symlink
                is produced only if ``hardcopy`` is True.
            velocity_offset (:obj:`float`, optional):
                Velocity offset to use when matching the spectral
                resolution between the template library and the galaxy
                spectra.
            match_resolution (:obj:`bool`, optional):
                Match the spectral resolution of the template library to
                the resolution provided by the ``cube``; the latter must
                be provided for this argument to have any use.
            resolution_fwhm (:obj:`float`, optional):
                The target resolution for the template library set by
                a constant FWHM of the resolution element in
                angstroms over all wavelengths. If None, the
                resolution is either matched to the binned spectra
                (see ``match_resolution``) or kept at the native
                resolution of the library.

        Returns:
            :class:`mangadap.proc.templatelibrary.TemplateLibrary`:
            The template library prepared for fitting the binned
            spectra.
        """

        # Allow output paths and files to be over-ridden
        _output_path = self.directory_path if output_path is None else output_path
        _hardcopy = self.tpl_hardcopy if hardcopy is None else hardcopy
        _symlink_dir = self.symlink_dir if symlink_dir is None else symlink_dir

        if resolution_fwhm is None:
            return TemplateLibrary(self.method['fitpar']['template_library_key'],
                                   velocity_offset=velocity_offset, cube=self.binned_spectra.cube,
                                   match_resolution=match_resolution,
                                   velscale_ratio=self.method['fitpar']['velscale_ratio'],
                                   output_path=_output_path, output_file=output_file,
                                   hardcopy=_hardcopy, symlink_dir=_symlink_dir,
                                   loggers=self.loggers, quiet=self.quiet)

        if output_file is None:
            _key = f'{self.method["fitpar"]["template_library_key"]}-FWHM{resolution_fwhm:.2f}'
            output_file = TemplateLibrary.default_paths(_key, cube=self.binned_spectra.cube,
                                                        output_path=_output_path)[1]

        # Set the spectral resolution
        wave = self.binned_spectra['WAVE'].data
        sres = SpectralResolution(wave, wave/resolution_fwhm, log10=True)

        velscale_ratio = 1 if self.method['fitpar']['velscale_ratio'] is None \
                            else self.method['fitpar']['velscale_ratio']
        spectral_step = spectral_coordinate_step(wave, log=True) / velscale_ratio

        return TemplateLibrary(self.method['fitpar']['template_library_key'], sres=sres,
                                velocity_offset=velocity_offset, spectral_step=spectral_step,
                                log=True, output_path=_output_path, output_file=output_file,
                                hardcopy=_hardcopy, symlink_dir=_symlink_dir,
                                loggers=self.loggers, quiet=self.quiet)

    def _initialize_primary_header(self, hdr=None):
        """
        Construct the primary header for the reference file.

        Args:
            hdr (`astropy.io.fits.Header`_, optional):
                Input base header for added keywords. If None, uses
                the datacube header from :attr:`binned_spectra` (if
                there is one) and then cleans the header using
                :func:`mangadap.util.fitsutil.DAPFitsUtil.clean_dap_primary_header`.

        Returns:
            `astropy.io.fits.Header`_: Initialized header object.
        """
        # Copy the from the DRP and clean it
        if hdr is None:
            hdr = self.binned_spectra.cube.prihdr.copy()
            hdr = DAPFitsUtil.clean_dap_primary_header(hdr)
        
        # Add keywords specific to this object
        hdr['AUTHOR'] = 'Kyle B. Westfall <westfall@ucolick.org>'
        hdr['VSTEP'] = (spectrum_velocity_scale(self.binned_spectra['WAVE'].data),
                        'Velocity step per spectral channel.')
        hdr['SCKEY'] = (self.method['key'], 'Stellar-continuum modeling method keyword')
        hdr['SCMINSN'] = (self.method['minimum_snr'], 'Minimum S/N of spectrum to include')

        # Guess kinematics
        # TODO: These are currently single numbers, but they may
        # eventually be vectors!  Include guess kinematics as a map?
        if self.guess_vel is not None:
            hdr['SCINPVEL'] = (self.guess_vel, 'Initial guess velocity')
        if self.guess_sig is not None:
            hdr['SCINPSIG'] = (self.guess_sig, 'Initial guess velocity dispersion')

        # Include the number of models
        hdr['NSCMOD'] = (self.nmodels, 'Number of unique stellar-continuum models')
        return hdr

    def _add_method_header(self, hdr):
        """
        Add method-specific metadata to the header.

        Note that the header object is both edited in-place and
        returned.

        Args:
            hdr (`astropy.io.fits.Header`_):
                Input base header for added keywords. Expected to
                have been initialized using
                :func:`_initialize_primary_header`.

        Returns:
            `astropy.io.fits.Header`_: Edited header object.
        """
        if self.method['fitclass'] is not None:
            try:
                hdr['SCTYPE'] = self.method['fitclass'].fit_type
                hdr['SCMETH'] = self.method['fitclass'].fit_method
            except:
                if not self.quiet and self.hardcopy:
                    warnings.warn('Fit class object does not have fit_type and/or fit_method ' \
                                  'attributes.  No parameters written to header.')
        if self.method['fitpar'] is not None:
            try:
                hdr = self.method['fitpar'].toheader(hdr)
            except:
                if not self.quiet and self.hardcopy:
                    warnings.warn('Fit parameter class has no toheader() function.  No ' \
                                  'parameters written to header.')
        return hdr


    def _finalize_cube_mask(self, mask):
        """
        Finalize the mask after the 2D mask has been reconstructed
        into a 3D cube.

        This mostly handles the masks for regions outside the
        datacube field of view.

        Note that the input mask is both edited in-place and
        returned.

        .. todo::

            - This needs to be abstracted for non-DRP datacubes.
            - Describe MAPMASK usage

        Args:
            mask (`numpy.ndarray`_):
                3D array with the current bitmask data.

        Returns:
            `numpy.ndarray`_: Edited bitmask data.
        """
        # TODO: I'm not sure this is correct for non-MaNGA data
        # Turn on the flag stating that the pixel wasn't used
        indx = self.binned_spectra.bitmask.flagged(self.binned_spectra.cube.mask,
                                                   flag=self.binned_spectra.do_not_fit_flags())
        mask[indx] = self.bitmask.turn_on(mask[indx], 'DIDNOTUSE')

        # Turn on the flag stating that the pixel has a foreground star
        indx = self.binned_spectra.bitmask.flagged(self.binned_spectra.cube.mask, flag='FORESTAR')
        mask[indx] = self.bitmask.turn_on(mask[indx], 'FORESTAR')

        # Propagate the MAPMASK to the full cube
        for f in ['DIDNOTUSE', 'FORESTAR', 'LOW_SNR']:
            indx = self.bitmask.flagged(self['MAPMASK'].data, f)
            mask[indx,:] = self.bitmask.turn_on(mask[indx,:], f)

        return mask

    def _assign_spectral_arrays(self):
        """
        Set :attr:`spectral_arrays`, which contains the list of
        extensions in :attr:`hdu` that contain spectral data.
        """
        self.spectral_arrays = ['FLUX', 'MASK']

    def _assign_image_arrays(self):
        """
        Set :attr:`image_arrays`, which contains the list of extensions
        in :attr:`hdu` that are on-sky image data.
        """
        self.image_arrays = ['BINID']

    def _get_missing_models(self, debug=False):
        """
        Find the spectra with missing models, either because there is
        no model or the bin is also missing.

        Args:
            debug (:obj:`bool`, optional):
                When checking the S/N limits, run
                :func:`mangadap.proc.spatiallybinnedspectra.SpatiallyBinnedSpectra.above_snr_limit`
                in debug mode.

        Returns:
            `numpy.ndarray`_: Sorted list of missing models.
        """
        good_snr = self.binned_spectra.above_snr_limit(self.method['minimum_snr'], debug=debug)
        return numpy.sort(self.binned_spectra['BINS'].data['BINID'][numpy.invert(good_snr)].tolist()
                                + self.binned_spectra.missing_bins) 

    def _construct_2d_hdu(self, good_snr, model_flux, model_mask, model_par):
        r"""
        Construct :attr:`hdu` that is held in memory for manipulation
        of the object. See :func:`construct_3d_hdu` to convert the
        object into a datacube.

        Args:
            good_snr (`numpy.ndarray`_):
                Boolean array selecting the bins with good S/N.
            model_flux (`numpy.ndarray`_):
                Best-fitting model spectra for each bin.
            model_mask (`numpy.ndarray`_):
                Bitmask flags for each model pixel.
            model_par (:class`~mangadap.proc.spectralfitting.StellarKinematicsFitDataTable`):
                Data table with the best-fitting model parameters.
        """
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, 'Constructing hdu ...')

        # Initialize the headers
        pri_hdr = self._initialize_primary_header()
        pri_hdr = self._add_method_header(pri_hdr)
        map_hdr = DAPFitsUtil.build_map_header(self.binned_spectra.cube.fluxhdr,
                                               'K Westfall <westfall@ucolick.org>')

        # Get the spatial map mask
        map_mask = numpy.zeros(self.spatial_shape, dtype=self.bitmask.minimum_dtype())
        # Add any spaxel not used because it was flagged by the binning
        # step
        indx = self.binned_spectra['MAPMASK'].data > 0
        map_mask[indx] = self.bitmask.turn_on(map_mask[indx], 'DIDNOTUSE')
        # Isolate any spaxels with foreground stars
        indx = self.binned_spectra.bitmask.flagged(self.binned_spectra['MAPMASK'].data, 'FORESTAR')
        map_mask[indx] = self.bitmask.turn_on(map_mask[indx], 'FORESTAR')
        # Get the bins that were blow the S/N limit
        indx = numpy.invert(DAPFitsUtil.reconstruct_map(self.spatial_shape,
                                                        self.binned_spectra['BINID'].data.ravel(),
                                                        good_snr, dtype='bool')) & (map_mask == 0)
        map_mask[indx] = self.bitmask.turn_on(map_mask[indx], 'LOW_SNR')

        # Get the bin ids with fitted models
        bin_indx = DAPFitsUtil.downselect_bins(self.binned_spectra['BINID'].data.ravel(),
                                               model_par['BINID']).reshape(self.spatial_shape)

        # Save the data to the hdu attribute
        self.hdu = fits.HDUList([fits.PrimaryHDU(header=pri_hdr),
                                 fits.ImageHDU(data=model_flux.data, name='FLUX'),
                                 fits.ImageHDU(data=model_mask, name='MASK'),
                                 self.binned_spectra['WAVE'].copy(),
                                 fits.ImageHDU(data=bin_indx, header=map_hdr, name='BINID'),
                                 fits.ImageHDU(data=map_mask, header=map_hdr, name='MAPMASK'),
                                 model_par.to_hdu(name='PAR')])

    @staticmethod
    def default_paths(cube, method_key, rdxqa_method, binning_method, output_path=None,
                      output_file=None):
        """
        Set the default directory and file name for the output file.

        Args:
            cube (:class:`mangadap.datacube.datacube.DataCube`):
                Datacube to analyze.
            method_key (:obj:`str`):
                Keyword designating the method used for the reduction
                assessments.
            rdxqa_method (:obj:`str`):
                The method key for the basic assessments of the datacube.
            binning_method (:obj:`str`):
                The method key for the spatial binning.
            output_path (:obj:`str`, `Path`_, optional):
                The path for the output file.  If None, the current working
                directory is used.
            output_file (:obj:`str`, optional):
                The name of the output "reference" file. The full path of the
                output file will be :attr:`directory_path`/:attr:`output_file`.
                If None, the default is to combine ``cube.output_root`` and the
                method keys.  The order of the keys is the order of operations
                (rdxqa, binning, stellar continuum).

        Returns:
            :obj:`tuple`: Returns a `Path`_ with the output directory and a
            :obj:`str` with the output file name.
        """
        directory_path = Path('.').resolve() if output_path is None \
                                else Path(output_path).resolve()
        method = f'{rdxqa_method}-{binning_method}-{method_key}'
        _output_file = f'{cube.output_root}-{method}.fits.gz' if output_file is None \
                            else output_file
        return directory_path, _output_file

    def file_name(self):
        """Return the name of the output file."""
        return self.output_file

    def file_path(self):
        """Return the full path to the output file."""
        if self.directory_path is None or self.output_file is None:
            return None
        return self.directory_path / self.output_file
    
    def info(self):
        return self.hdu.info()

    def all_spectrum_flags(self):
        return ['DIDNOTUSE', 'FORESTAR', 'LOW_SNR', 'ARTIFACT', 'OUTSIDE_RANGE', 'EML_REGION',
                'TPL_PIXELS', 'TRUNCATED', 'PPXF_REJECT', 'INVALID_ERROR', 'FIT_FAILED',
                'NEAR_BOUND']

    def all_except_emission_flags(self):
        return ['DIDNOTUSE', 'FORESTAR', 'LOW_SNR', 'ARTIFACT', 'OUTSIDE_RANGE', 'TPL_PIXELS',
                'TRUNCATED', 'PPXF_REJECT', 'INVALID_ERROR', 'FIT_FAILED', 'NEAR_BOUND']

    def fit(self, binned_spectra, guess_vel, guess_sig=None, output_path=None, output_file=None,
            hardcopy=True, symlink_dir=None, tpl_hardcopy=None, overwrite=None, loggers=None,
            quiet=False):
        """
        Fit the binned spectra using the specified method.

        Args:
            binned_spectra (:class:`mangadap.proc.spatiallybinnedspectra.SpatiallyBinnedSpectra`):
                The binned spectra to fit. Must be able to access
                parent :class:`mangadap.datacube.datacube.DataCube`
                attribute (``cube``).
            guess_vel (:obj:`float`, `numpy.ndarray`):
                A single or spectrum-dependent initial estimate of the
                redshift in km/s (:math:`cz`).
            guess_sig (:obj:`float`, `numpy.ndarray`, optional):
                A single or spectrum-dependent initial estimate of the
                velocity dispersion in km/s.
            output_path (:obj:`str`, `Path`_, optional):
                The path for the output file.  If None, the current working
                directory is used.
            output_file (:obj:`str`, optional):
                The name of the output "reference" file. The full path of the output
                file will be :attr:`directory_path`/:attr:`output_file`.  If None,
                the default is to combine ``cube.output_root`` and the method keys.
                The order of the keys is the order of operations (rdxqa, binning,
                stellar continuum).  See
                :func:`~mangadap.proc.stellarcontinuummodel.StellarContinuumModel.default_paths`.
            hardcopy (:obj:`bool`, optional):
                Flag to write the HDUList attribute to disk.  Default is
                True; if False, the HDUList is only kept in memory and
                would have to be reconstructed.
            symlink_dir (:obj:`str`, optional):
                Create a symbolic link to the created file in the
                supplied directory.  Default is to produce no symbolic
                link. 
            tpl_hardcopy (:obj:`bool`, optional):
                Save the processed template library used during the fit.  If
                True, the output path and symlink directory for the template
                library are identical to the main reference file.
            overwrite (:obj:`bool`, optional):
                Overwrite any existing files.  Default is to use any
                existing file instead of redoing the analysis and
                overwriting the existing output.
            loggers (:obj:`list`, optional):
                List of `logging.Logger`_ objects to log progress;
                ignored if quiet=True.  Logging is done using
                :func:`mangadap.util.log.log_output`.  Default is no
                logging.
            quiet (:obj:`bool`, optional):
                Suppress all terminal and logging output. 
        """
        # Initialize the reporting
        if loggers is not None:
            self.loggers = loggers
        self.quiet = quiet

        # SpatiallyBinnedSpectra object always needed
        if binned_spectra is None:
            raise ValueError('Must provide spectra object for fitting.')
        if not isinstance(binned_spectra, SpatiallyBinnedSpectra):
            raise TypeError('Must provide a valid SpatiallyBinnedSpectra object!')
        if binned_spectra.hdu is None:
            raise ValueError('Provided SpatiallyBinnedSpectra object is undefined!')
        self.binned_spectra = binned_spectra

        self.shape = self.binned_spectra.shape
        self.spatial_shape =self.binned_spectra.spatial_shape
        self.nspec = self.binned_spectra.nspec
        self.spatial_index = self.binned_spectra.spatial_index.copy()
        self.nwave = self.binned_spectra.nwave
        
        # Get the guess kinematics
        if guess_vel is not None:
            self.guess_vel=guess_vel
        if self.guess_vel is None:
            raise ValueError('Must provide initial guess velocity.')
        self.guess_sig = 100.0 if guess_sig is None else guess_sig

        #---------------------------------------------------------------
        # Get the good spectra
#        debug = True
        debug = False
        good_snr = self.binned_spectra.above_snr_limit(self.method['minimum_snr'], debug=debug)

        # Report
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, '-'*50)
            log_output(loggers, 1, logging.INFO, f'{"STELLAR CONTINUUM FITTING":^50}')
            log_output(self.loggers, 1, logging.INFO, '-'*50)
            log_output(self.loggers, 1, logging.INFO,
                       f'Number of binned spectra: {self.binned_spectra.nbins}')
            if len(self.binned_spectra.missing_bins) > 0:
                log_output(self.loggers, 1, logging.INFO,
                           f'Missing bins: {len(self.binned_spectra.missing_bins)}')
            log_output(self.loggers, 1, logging.INFO,
                       f'With good S/N and to fit: {numpy.sum(good_snr)}')

        # TODO: Is there a better way to handle this?
        if numpy.sum(good_snr) == 0:
            raise ValueError('No good spectra to fit!')

        #---------------------------------------------------------------
        # Fill in any remaining binning parameters
        self.tpl_hardcopy = tpl_hardcopy
#        self._fill_method_par()
        if self.method['fitpar'] is not None and hasattr(self.method['fitpar'], 'fill'):
            self.method['fitpar'].fill(self.binned_spectra, guess_vel=self.guess_vel,
                                       guess_sig=self.guess_sig, output_path=self.directory_path,
                                       hardcopy=self.tpl_hardcopy, symlink_dir=self.symlink_dir,
                                       loggers=self.loggers, quiet=self.quiet)

        _overwrite = self.method['overwrite'] if overwrite is None else overwrite

        # (Re)Set the output paths
        self.directory_path, self.output_file \
                = StellarContinuumModel.default_paths(self.binned_spectra.cube,
                                                      self.method['key'],
                                                      self.binned_spectra.rdxqa.method['key'],
                                                      self.binned_spectra.method['key'],
                                                      output_path=output_path,
                                                      output_file=output_file)

        #---------------------------------------------------------------
        # Check that the file path is defined
        ofile = self.file_path()
        if ofile is None:
            raise ValueError('File path for output file is undefined!')

        # Report
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, f'Output path: {self.directory_path}')
            log_output(self.loggers, 1, logging.INFO, f'Output file: {self.output_file}')
        
        #---------------------------------------------------------------
        # If the file already exists, and not ovewriting, just read the
        # file
        self.symlink_dir = symlink_dir
        if ofile.exists() and not overwrite:
            self.hardcopy = True
            if not self.quiet:
                log_output(self.loggers, 1, logging.INFO, 'Using existing file')
            self.read(checksum=self.checksum, debug=debug)
            # Make sure the symlink exists
            if self.symlink_dir is not None:
                create_symlink(str(ofile), self.symlink_dir, overwrite=overwrite)
            if not self.quiet:
                log_output(self.loggers, 1, logging.INFO, '-'*50)
            return

        #---------------------------------------------------------------
        # Fit the spectra
        # Mask should be fully defined within the fitting function
        model_wave, model_flux, model_mask, model_par \
                = self.method['fitfunc'](self.binned_spectra, gpm=good_snr,
                                         par=self.method['fitpar'], loggers=self.loggers,
                                         quiet=self.quiet, debug=debug)

        # The number of models returned should match the number of
        # "good" spectra

        # TODO: Include failed fits in "missing" models?

        # DEBUG
        if model_flux.shape[0] != numpy.sum(good_snr):
            raise ValueError('Unexpected returned shape of fitted continuum models.')

        # Set the number of models and the missing models
        self.nmodels = model_flux.shape[0]
        self.missing_models = self._get_missing_models(debug=debug)

        # Report
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, f'Fitted models: {self.nmodels}')
            if len(self.missing_models) > 0:
                log_output(self.loggers, 1, logging.INFO,
                           f'Missing models: {len(self.missing_models)}')

        # Construct the 2d hdu list that only contains the fitted models
        self.hardcopy = hardcopy
        self._construct_2d_hdu(good_snr, model_flux, model_mask, model_par)

        #---------------------------------------------------------------
        # Write the data, if requested
        if self.hardcopy:
            if not self.directory_path.exists():
                self.director_path.mkdir(parents=True)
            self.write(overwrite=overwrite)
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, '-'*50)

    def construct_3d_hdu(self):
        """
        Reformat the model spectra into a cube matching the shape of
        the DRP fits file.
        """
        # Report
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, 'Constructing stellar continuum datacube ...')

        bin_indx = self.hdu['BINID'].data.copy()
        model_flux = self.hdu['FLUX'].data.copy()
        model_mask = self.hdu['MASK'].data.copy()

        flux, mask = DAPFitsUtil.reconstruct_cube(self.shape, bin_indx.ravel(),
                                                  [model_flux, model_mask])

        mask = self._finalize_cube_mask(mask)

        # Primary header is identical regardless of the shape of the
        # extensions
        hdr = self.hdu['PRIMARY'].header.copy()
        cube_hdr = DAPFitsUtil.build_cube_header(self.binned_spectra.cube,
                                                 'K Westfall <westfall@ucolick.org>')

        # Return the converted hdu without altering the internal hdu
        return fits.HDUList([ fits.PrimaryHDU(header=hdr),
                              fits.ImageHDU(data=flux, header=cube_hdr, name='FLUX'),
                              fits.ImageHDU(data=mask, header=cube_hdr, name='MASK'),
                              self.hdu['WAVE'].copy(),
                              self.hdu['BINID'].copy(),
                              self.hdu['PAR'].copy()
                            ])

    # Exact same function as used by SpatiallyBinnedSpectra
    def write(self, match_datacube=False, overwrite=False):
        """
        Write the hdu object to the file.

        Args:
            match_datacube (:obj:`bool`, optional):
                Match the shape of the data arrays to the input
                datacube. I.e., convert them to 3D and replicate the
                binned spectrum to each spaxel in the bin.
            overwrite (:obj:`bool`, optional):
                Overwrite any existing file.
        """
        # Convert the spectral arrays in the HDU to a 3D cube and write
        # it
        if match_datacube:
            hdu = self.construct_3d_hdu()
            DAPFitsUtil.write(hdu, str(self.file_path()), overwrite=overwrite, checksum=True,
                              symlink_dir=self.symlink_dir, loggers=self.loggers, quiet=self.quiet)
            return
        # Just write the unique (2D) data
        DAPFitsUtil.write(self.hdu, str(self.file_path()), overwrite=overwrite, checksum=True,
                          symlink_dir=self.symlink_dir, loggers=self.loggers, quiet=self.quiet) 

    def read(self, ifile=None, strict=True, checksum=False, debug=False):
        """
        Read an existing file with a previously fit set of continuum
        models.

        Args:
            ifile (:obj:`str`, optional):
                Name of the file with the data. Default is to use the
                name provided by :func:`file_path`.
            strict (:obj:`bool`, optional):
                Force a strict reading of the file to make sure that
                it adheres to the expected format. At the moment,
                this only checks to make sure the method keyword
                printed in the header matches the expected value in
                :attr:`method`.
            checksum (:obj:`bool`, optional):
                Use the checksum in the header to check for
                corruption of the data.
            debug (:obj:`bool`, optional):
                Run check for missing models in debugging mode.

        Raises:
            FileNotFoundError:
                Raised if the file does not exist.
            ValueError:
                Raised if ``strict`` is true and the header keyword
                ``SCKEY`` does not match the method keyword.
        """
        if ifile is None:
            ifile = self.file_path()
        if not ifile.exists():
            raise FileNotFoundError('File does not exist!: {0}'.format(ifile))

        if self.hdu is not None:
            self.hdu.close()

#        self.hdu = fits.open(ifile, checksum=checksum)
        self.hdu = DAPFitsUtil.read(str(ifile), checksum=checksum)

        # Confirm that the internal method is the same as the method
        # that was used in writing the file
        if self.hdu['PRIMARY'].header['SCKEY'] != self.method['key']:
            if strict:
                raise ValueError('Keywords in header do not match specified method keywords!')
            elif not self.quiet:
                warnings.warn('Keywords in header do not match specified method keywords!')
        # TODO: "strict" should also check other aspects of the file to
        # make sure that the details of the method are also the same,
        # not just the keyword

        # Attempt to read the input guess velocity and sigma
        try:
            self.guess_vel = self.hdu['PRIMARY'].header['SCINPVEL']
            self.guess_sig = self.hdu['PRIMARY'].header['SCINPSIG']
        except:
            if not self.quiet:
                warnings.warn('Unable to read input guess kinematics from file header!')
            self.guess_vel = None
            self.guess_sig = None

        # Attempt to read the modeling parameters
        if self.method['fitpar'] is not None and callable(self.method['fitpar'].fromheader):
            self.method['fitpar'].fromheader(self.hdu['PRIMARY'].header)

        self.nmodels = self.hdu['FLUX'].shape[0]
        self.missing_models = self._get_missing_models(debug=debug)

    def copy_to_array(self, ext='FLUX', waverange=None, include_missing=False):
        r"""
        Wrapper for :func:`mangadap.util.fitsutil.DAPFitsUtil.copy_to_array`
        specific for :class:`StellarContinuumModel`.

        Return a copy of the selected data array. The array size is
        always :math:`N_{\rm models} \times N_{\rm wavelength}`; i.e.,
        the data is always flattened to two dimensions and the unique
        spectra are selected.

        Args:
            ext (:obj:`str`, optional):
                Name of the extension from which to draw the data.
                Must be allowed for the current :attr:`mode`; see
                :attr:`data_arrays`.
            waverange (array-like, optional):
                Two-element array with the first and last wavelength
                to include in the computation. If None, use the full
                wavelength range.
            include_missing (:obj:`bool`, optional):
                Create an array with a size that accommodates the
                missing models.

        Returns:
            `numpy.ndarray`_: A 2D array with a copy of the data from
            the selected extension.
        """
        return DAPFitsUtil.copy_to_array(self.hdu, ext=ext, allowed_ext=self.spectral_arrays,
                                         waverange=waverange,
                                    missing_bins=self.missing_models if include_missing else None,
                                         nbins=self.nbins,
                                    unique_bins=DAPFitsUtil.unique_bins(self.hdu['BINID'].data))

    def copy_to_masked_array(self, ext='FLUX', flag=None, waverange=None, include_missing=False):
        """
        Wrapper for
        :func:`mangadap.util.fitsutil.DAPFitsUtil.copy_to_masked_array`
        specific for :class:`StellarContinuumModel`.

        Return a copy of the selected data array as a masked array.
        This is functionally identical to :func:`copy_to_array`,
        except the output format is a `numpy.ma.MaskedArray`_. The
        pixels that are considered to be masked can be specified
        using `flag`.
        
        Args:
            ext (:obj:`str`, optional):
                Name of the extension from which to draw the data.
                Must be allowed for the current :attr:`mode`; see
                :attr:`data_arrays`.
            flag (:obj:`str`, :obj:`list`, optional):
                (List of) Flag names that are considered when
                deciding if a pixel should be masked. The names
                *must* be a valid bit name as defined by
                :attr:`bitmask`. If not provided, *ANY* non-zero mask
                bit is omitted.
            waverange (array-like, optional):
                Two-element array with the first and last wavelength
                to include in the computation. If None, use the full
                wavelength range.
            include_missing (:obj:`bool`, optional):
                Create an array with a size that accommodates the
                missing models.

        Returns:
            `numpy.ma.MaskedArray`_: A 2D array with a copy of the
            data from the selected extension.
        """
        return DAPFitsUtil.copy_to_masked_array(self.hdu, ext=ext, mask_ext='MASK', flag=flag,
                                                bitmask=self.bitmask,
                                                allowed_ext=self.spectral_arrays,
                                                waverange=waverange,
                                    missing_bins=self.missing_models if include_missing else None,
                                                nbins=self.nmodels,
                                        unique_bins=DAPFitsUtil.unique_bins(self.hdu['BINID'].data))

    def construct_models(self, replacement_templates=None, redshift_only=False,
                         deredshift=False, corrected_dispersion=False):
        """
        Reconstruct the best fitting models.

        This is largely a wrapper for a :func:`construct_models`
        method provided by the :attr:`method` fitting class. As such,
        this method will raise an exception if the fitting class has
        no such method.

        Args:
            replacement_templates (:class:`~mangadap.proc.templatelibary.TemplateLibrary`, optional):
                Instead of using the template library used to fit the
                data, use this library instead. The number of spectra
                in the replacement library **must** match the
                original number of templates. This is only really
                useful when replacing the original set of templates
                with ones at a different spectral resolution.
            redshift_only (:obj:`bool`, optional):
                When constructing the models, only redshift the
                spectra, effectively providing the best-fitting model
                with a velocity dispersion of 0 km/s.
            deredshift (:obj:`bool`, optional):
                Return the best-fitting models in the rest frame.
            corrected_dispersion (:obj:`bool`, optional):
                Return the best-fitting models with a LOSVD that has
                been corrected for any resolution difference between
                the fitted templates and the binned spectra.

        Returns:
            `numpy.ndarray`_: The best-fitting models.

        Raises:
            ValueError: 
                Raised if there is no defined fitting class.
            AttributeError:
                Raised if the fitting class has no method called
                :func:`construct_models`.
            TypeError:
                Raised if the provided ``replacement_templates`` are
                not an instance of
                :class:`mangadap.proc.templatelibary.TemplateLibrary`.
        """
        if self.method['fitclass'] is None:
            raise ValueError('No class object available for constructing the model!')
        if not callable(self.method['fitclass'].construct_models):
            raise AttributeError('Provided fit class object has no callable ' \
                                 '\'construct_models\' attribute!')
        if replacement_templates is not None \
                and not isinstance(replacement_templates, TemplateLibrary):
            raise TypeError('Provided template replacements must have type TemplateLibrary.')

        # Only the selected models are constructed, others are masked
        select = numpy.invert(self.bitmask.flagged(self.hdu['PAR'].data['MASK'],
                                            flag=['NO_FIT', 'FIT_FAILED', 'INSUFFICIENT_DATA']))
        templates = self.method['fitpar']['template_library'] if replacement_templates is None \
                                            else replacement_templates

        # The only reason it needs the flux data is to get the shape, so
        # it doesn't matter that I'm passing the model flux
        return self.method['fitclass'].construct_models(templates['WAVE'].data,
                                                        templates['FLUX'].data,
                                                        self.hdu['WAVE'].data,
                                                        self.hdu['FLUX'].data.shape,
                                                        self.hdu['PAR'].data, select=select,
                                                        redshift_only=redshift_only,
                                                        deredshift=deredshift,
                                                        corrected_dispersion=corrected_dispersion,
                                                        dvtol=1e-9)

    # TODO: Get rid of dispaxis?
    # TODO: Make this a general routine.
    @staticmethod
    def reset_continuum_mask_window(continuum, dispaxis=1, quiet=False):
        """
        Reset the mask of the stellar continuum to a continuous window
        from the minimum to maximum valid wavelength.

        .. todo::
            - Allow continuum to be n-dimensional and use
              `numpy.apply_along_axis`.

        Args:   
            continuum (`numpy.ma.MaskedArray`_):
                Stellar continuum models.
            dispaxis (:obj:`int`, optional):
                Wavelength axis in provided ``continuum`` object.
            quiet (:obj:`bool`, optional):
                Suppress terminal output.

        Returns:
            `numpy.ma.MaskedArray`_: Continuum object but unmasked
            over the full continuous window between the minimum and
            maximum valid wavelength.
        """
        spatial_shape = DAPFitsUtil.get_spatial_shape(continuum.shape, dispaxis)
        if len(spatial_shape) != 1:
            raise ValueError('Input array should be two-dimensional!')

        # No pixels are masked!
        if numpy.sum(numpy.ma.getmaskarray(continuum)) == 0:
            return continuum

        _continuum = continuum.copy() if dispaxis == 1 else continuum.copy().T

        nspec,npix = _continuum.shape
        pix = numpy.ma.MaskedArray(numpy.array([ numpy.arange(npix) ]*nspec),
                                   mask=numpy.ma.getmaskarray(_continuum).copy())

        min_good_pix = numpy.ma.amin(pix, axis=dispaxis)
        max_good_pix = numpy.ma.amax(pix, axis=dispaxis)
        for c,s,e in zip(_continuum, min_good_pix,max_good_pix+1):
            if isinstance(s, numpy.ma.core.MaskedConstant) \
                    or isinstance(e, numpy.ma.core.MaskedConstant):
                continue
            numpy.ma.getmaskarray(c)[s:e] = False

        return _continuum if dispaxis == 1 else _continuum.T


    def unmasked_continuum_model(self, replacement_templates=None, redshift_only=False,
                                 deredshift=False, corrected_dispersion=False):
        """
        Return the stellar continuum unmasked over the full continuous
        fitting region.

        Models are reconstructed based on the model parameters if new
        template fluxes are provided or if the returned models should
        not include the LOSVD convolution (redshift_only=True). As
        such, this is largely a wrapper for :func:`construct_models`
        and :func:`reset_continuum_mask_window`.

        Args:
            replacement_templates (:class:`mangadap.proc.templatelibary.TemplateLibrary`, optional):
                Instead of using the template library used to fit the
                data, use this library instead. The number of spectra
                in the replacement library **must** match the
                original number of templates. This is only really
                useful when replacing the original set of templates
                with ones at a different spectral resolution.
            redshift_only (:obj:`bool`, optional):
                When constructing the models, only redshift the
                spectra, effectively providing the best-fitting model
                with a velocity dispersion of 0 km/s.
            deredshift (:obj:`bool`, optional):
                Return the best-fitting models in the rest frame.
            corrected_dispersion (:obj:`bool`, optional):
                Return the best-fitting models with a LOSVD that has
                been corrected for any resolution difference between
                the fitted templates and the binned spectra.

        Returns:
            `numpy.ndarray`_: The best-fitting models.
        """
        # Get the models for the binned spectra
        reconstruct = replacement_templates is not None or redshift_only or corrected_dispersion
        continuum = self.construct_models(replacement_templates=replacement_templates,
                                          redshift_only=redshift_only, deredshift=deredshift,
                                          corrected_dispersion=corrected_dispersion) \
                        if reconstruct \
                        else self.copy_to_masked_array(flag=self.all_spectrum_flags())
        return StellarContinuumModel.reset_continuum_mask_window(continuum, quiet=self.quiet)

    def fill_to_match(self, binid, replacement_templates=None, redshift_only=False,
                      deredshift=False, corrected_dispersion=False, missing=None):
        """
        Get the stellar-continuum model from this objects that matches
        the input bin ID matrix.

        The output is a 2D matrix ordered by the bin ID; any skipped
        index numbers in the maximum of the union of the unique
        numbers in the ``binid`` and ``missing`` input are masked.

        All this does is find the stellar continuum from spaxels in
        this model and match them to the input spaxels. Depending on
        how differently the data sets have been binned, this says
        nothing about whether or not the returned models are a good
        fit to the data!

        Args:
            binid (`numpy.ndarray`_):
                A map of the bin IDs with shape that matches the
                spatial shape of the parent datacube; see
                :attr:`spatial_shape`. This provides the mapping
                between the current stellar continuum models and the
                output data array.
            replacement_templates (:class:`mangadap.proc.templatelibary.TemplateLibrary`, optional):
                Instead of using the template library used to fit the
                data, use this library instead. The number of spectra
                in the replacement library **must** match the
                original number of templates. This is only really
                useful when replacing the original set of templates
                with ones at a different spectral resolution.
            redshift_only (:obj:`bool`, optional):
                When constructing the models, only redshift the
                spectra, effectively providing the best-fitting model
                with a velocity dispersion of 0 km/s.
            deredshift (:obj:`bool`, optional):
                Return the best-fitting models in the rest frame.
            corrected_dispersion (:obj:`bool`, optional):
                Return the best-fitting models with a LOSVD that has
                been corrected for any resolution difference between
                the fitted templates and the binned spectra.
            missing (array-like, optional):
                Include empty spectra for binids with these numbers.
                Ignored if None or an empty array.

        Returns:
            `numpy.ndarray`_: The best-fitting models matched to the
            input bin ID array.

        Raises:
            ValueError:
                Raised if the shape of ``binid`` does not match
                :attr:`spatial_shape`.
        """
        # The input bin id array must have the same shape as self
        if binid.shape != self.spatial_shape:
            raise ValueError('Input bin ID matrix has incorrect shape.')

        # Construct the best-fitting models
        best_fit_continuum = self.unmasked_continuum_model(
                                    replacement_templates=replacement_templates,
                                    redshift_only=redshift_only, deredshift=deredshift,
                                    corrected_dispersion=corrected_dispersion)

        # Get the number of output continuua
        nbins = numpy.amax(binid).astype(int)+1 if missing is None or len(missing) == 0 else \
                    numpy.amax( numpy.append(binid.ravel(), missing) ).astype(int)+1

        # Map the BINID to the spectrum index, assuming bins are sorted
        u, indx, reconstruct = numpy.unique(self['BINID'].data.ravel(), return_index=True,
                                            return_inverse=True)
#        u_bin_indx = numpy.arange(len(u))-1
        u_bin_indx = numpy.arange(u.size)
        if u[0] == -1:
            u_bin_indx -= 1
        _bin_indx = u_bin_indx[reconstruct].reshape(self.spatial_shape)

        # Fill in bins with no models with masked zeros
        continuum = numpy.ma.zeros((nbins,best_fit_continuum.shape[1]), dtype=float)
        continuum[:,:] = numpy.ma.masked
        for i,j in zip(binid.ravel(), _bin_indx.ravel()):
            if i < 0 or j < 0:
                continue
            continuum[i,:] = best_fit_continuum[j,:]

        return continuum

    def matched_kinematics(self, binid, redshift=None, dispersion=100.0, constant=False, cz=False,
                           corrected=False, min_dispersion=None, nearest=False, missing=None):
        r"""
        Get the stellar kinematics from this objects that matches the
        input bin ID matrix.

        This method has a similar intent to :func:`fill_to_match`,
        but instead of matching the model spectra, its intent is to
        match the stellar kinematics.

        For spectra the were not fit but have a bin ID (see
        ``binid``), the method uses the median (unmasked) kinematic
        measurement if no default values are provided (see
        ``redshift`` and ``dispersion``).

        Args:
            binid (`numpy.ndarray`_):
                2D array with the bin ID associate with each spaxel
                in the datacube. Shape must be the same as
                :attr:`spatial_shape`.
            redshift (:obj:`float`, optional):
                The default redshift to use for spectra without a
                stellar-continuum fit. Default is to use the median
                of the unmasked measurements
            dispersion (:obj:`float`, optional):
                The default velocity dispersion to use for spectra
                without a stellar-continuum fit. Default is 100 km/s.
            constant (:obj:`bool`, optional):
                Force the function to return a constant redshift and
                dispersion for each spectrum, regardless of any
                fitted kinematics.
            cz (:obj:`bool`, optional):
                Return the redshift as cz in km/s, as opposed to
                unitless z.
            corrected (:obj:`bool`, optional):
                Return the velocity dispersions corrected for any
                resolution difference between the templates and the
                galaxy data. If False, return the uncorrected data.
            min_dispersion (:obj:`float`, optional):
                Impose a minimum dispersion.
            nearest (:obj:`bool`, optional):
                Instead of the median of the results for the spectra
                that were not fit, use the value from the nearest
                bin.
            missing (array-like, optional):
                Any bin ID numbers missing from the input ``binid``
                image needed for constructing the output matched
                data.

        Returns:
            :obj:`tuple`: Returns `numpy.ndarray`_ objects with a
            redshift (or cz) and dispersion for each binned spectrum.
            Shape is :math:`(N_{\rm bin},)`.

        Raises:
            TypeError:
                Raised if the input redshift or dispersion values are
                not single numbers.
            ValueError:
                Raised if the shape of ``binid`` does not match
                :attr:`spatial_shape`.
        """
        # Check input
        if redshift is not None and not isinstance(redshift, (float,int)):
            raise TypeError('Redshift must be a single number or None.')
        if dispersion is not None and not isinstance(dispersion, (float,int)):
            raise TypeError('Dispersion must be a single number or None.')
        if binid.shape != self.spatial_shape:
            raise ValueError('Input bin ID map must match the spatial shape of the DRP cube.')

        # Get the number of output kinematics
        nbins = numpy.amax(binid).astype(int)+1 if missing is None or len(missing) == 0 else \
                    numpy.amax( numpy.append(binid.ravel(), missing) ).astype(int)+1

        # Mask bad values
        mask = self.bitmask.flagged(self.hdu['PAR'].data['MASK'].copy(),
                                    [ 'NO_FIT', 'FIT_FAILED', 'INSUFFICIENT_DATA', 'NEAR_BOUND' ])

        if numpy.all(mask):
            # Everything is masked so use input guesses
            warnings.warn('All stellar continuum fits have been masked!  Using input guesses.')
            str_z = numpy.ma.MaskedArray(self.method['fitpar']['guess_redshift'].copy())
            str_d = numpy.ma.MaskedArray(self.method['fitpar']['guess_dispersion'].copy())
            # Select one per *stellar-continuum model*
            str_z = str_z[self['PAR'].data['BINID_INDEX']]
            str_d = str_d[self['PAR'].data['BINID_INDEX']]
        else:
            # Pull out the data
            str_z = numpy.ma.MaskedArray(self.hdu['PAR'].data['KIN'][:,0].copy(), mask=mask) \
                            / astropy.constants.c.to('km/s').value
            str_d = numpy.ma.MaskedArray(self.hdu['PAR'].data['KIN'][:,1].copy(), mask=mask)

            # Apply the sigma corrections if requested.  Values with the
            # correction larger than the measured value will be masked!
            if corrected:
                sigma_corr = numpy.ma.MaskedArray(self.hdu['PAR'].data['SIGMACORR_EMP'], mask=mask)
                str_d = numpy.ma.sqrt( numpy.square(str_d) - numpy.square(sigma_corr) )

        # Set the default values to use when necessary
        _redshift = numpy.ma.median(str_z) if redshift is None else redshift
        _dispersion = numpy.ma.median(str_d) if dispersion is None else dispersion
        if _dispersion is numpy.ma.masked:
            _dispersion = numpy.median(self.method['fitpar']['guess_dispersion'])
        if min_dispersion is not None and _dispersion < min_dispersion:
            _dispersion = min_dispersion

        # Just return the single value
        if constant:
            str_z = numpy.full(nbins, _redshift * astropy.constants.c.to('km/s').value
                                        if cz else _redshift, dtype=float)
            str_d = numpy.full(nbins, _dispersion, dtype=float)
            return str_z, str_d

        # Fill masked values with the nearest datum
        if nearest:
            # Fill masked values with the nearest datum
            replace = numpy.ma.getmaskarray(str_z) | numpy.ma.getmaskarray(str_d)

            if numpy.all(replace):
                # All are bad so replace with _redshift and _dispersion
                str_z = numpy.full(str_z.shape, _redshift, dtype=float)
                str_d = numpy.full(str_d.shape, _dispersion, dtype=float)
            elif numpy.any(replace):
                best_fit_kinematics = numpy.ma.append([str_z], [str_d], axis=0).T
                valid_bins = numpy.unique(self['BINID'].data)
                if valid_bins[0] == -1:
                    valid_bins = valid_bins[1:]
                coo = self.binned_spectra['BINS'].data['SKY_COO'][valid_bins,:]
                kinematics = replace_with_data_from_nearest_coo(coo, best_fit_kinematics, replace)
                str_z = kinematics[:,0]
                str_d = kinematics[:,1]
        else:
            # Fill any masked values with the single estimate
            str_z = str_z.filled(_redshift)
            str_d = str_d.filled(_dispersion)

        # Map the BINID to the spectrum index, assuming bins are sorted
        # and that the BINID map has -1 BINID values
        u, indx, reconstruct = numpy.unique(self['BINID'].data.ravel(), return_index=True,
                                            return_inverse=True)
#        u_bin_indx = numpy.arange(len(u))-1
        u_bin_indx = numpy.arange(u.size)
        if u[0] == -1:
            u_bin_indx -= 1
        _bin_indx = u_bin_indx[reconstruct].reshape(self.spatial_shape)

        # Match the kinematics to the output bin ID map
        _str_z = numpy.ma.masked_all(nbins, dtype=float)
        _str_d = numpy.ma.masked_all(nbins, dtype=float)
        for i,j in zip(binid.ravel(), _bin_indx.ravel()):
            if i < 0 or j < 0:
                continue
            _str_z[i] = str_z[j]
            _str_d[i] = str_d[j]

        str_z = _str_z.filled(_redshift)
        str_d = _str_d.filled(_dispersion)
        # Convert to cz velocities (km/s)
        if cz:
            str_z *= astropy.constants.c.to('km/s').value
        return str_z, str_d

# TODO: Function out of date!  Keep it around for now though...
#    def matched_template_flags(self, binned_spectra):
#        """
#        Return the templates used during the fit to each spectrum,
#        matched to the spectra in the binned_spectra object.  For
#        spectra with no models, the flags are all set to true.
#        """
#        if binned_spectra is self.binned_spectra:
#            usetpl = self.hdu['PAR'].data['TPLWGT'] > 0
#
#            # Number of models matches the numbers of bins
#            if binned_spectra.nbins == self.nmodels:
##                print('returning usetpl')
#                return usetpl
#    
#            # Fill in bins with no models with masked zeros
##            print('Fill in bins with no models with masked zeros')
#            ntpl = self.method['fitpar']['template_library'].ntpl
#            _usetpl = numpy.ones((binned_spectra.nbins,ntpl), dtype=bool)
#            for i,j in enumerate(self.hdu['PAR'].data['BINID_INDEX']):
#                _usetpl[j,:] = usetpl[i,:]
#            return _usetpl
#
#        raise NotImplementedError('Can only match to internal binned_spectra.')

