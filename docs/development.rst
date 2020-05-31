
.. include:: include/links.rst

.. _development:

Development Guidelines
======================

The development strategy for the DAP has been to allow for high-level
interaction with the software to make subtle changes to its execution
parameters.  For help with the high-level interaction with the DAP
parameters, see the description of its :ref:`execution`.

We have also recently abstracted the input to the DAP to allow for
custom :class:`~mangadap.datacube.datacube.DataCube` objects to be
provided by the user, enabling the DAP to analyze non-MaNGA
data; see :ref:`datacube`.

Additionally, we have constructed the low-level, core algorithms in a
way that is largely agnostic to the data being analyzed. This should
allow users to use DAP methods/classes in new scripts tailored to
their own needs and analysis of non-MaNGA data. For example, we
demonstrate :ref:`fitonespec` using the DAP fitting modules.

We are very interested in having the community building on the existing
DAP algorithms.  As a means of introducing some of the nuances of its
structure, the following is list of some of the guidelines we tried to
follow in the development of the code:

 * Minimize code repetition by isolating repeated functionality into a
   flexible function or class.
 * Provide utility level classes that ease interaction with the input
   and output data.
 * Provide utility level classes that ease inclusion of new spectral
   template libraries, emission-line parameters, and spectral-index
   parameters for use in the main analysis classes.
 * Compartmentalize analysis into low-level functions that are
   independent of the specifics of the MaNGA data.
 * Generalize the interfaces of these low-level functions such that they
   can be incorporated into user-level scripts.
 * Where possible and appropriate, ensure the low-level functions
   operate on arrays of spectra that can effectively be considered
   independent of one another.
 * Minimize hard-coded behavior in favor of parameter or
   configuration files.
 * Provide an infrastructure that includes hooks for user-provided
   functions to override the base-level DAP analysis functions while
   still maintaining the output DAP datamodel.
 * Use an `astropy.io.fits.HDUList`_ object as the primary interface to
   the high-level classes to generalize interaction and allow them to be
   directly written to their reference fits files with a well-formed
   datamodel.

.. _templatelibrary-dev-example:

TemplateLibrary example
-----------------------

As an example of what some of this looks like in practice, consider the
:ref:`templatelibrary-usage` for a
:class:`~mangadap.proc.templatelibrary.TemplateLibrary`.

The :class:`~mangadap.proc.templatelibrary.TemplateLibrary` object is
setup to read and process the template library used to model the
stellar continuum. The survey-level execution of the DAP currently
uses the ``MILESHC`` library for the stellar kinematics and the
``MASTARHC`` library to model the stellar continuum in the
emission-line module. The way that the DAP knows where to find the
templates and how to process them is via the configuration files
located in ``$MANGADAP_DIR/mangadap/config/spectral_templates``. The
one specific to the ``MILESHC`` library is
``$MANGADAP_DIR/mangadap/config/spectral_templates/miles_hc.ini`` and
looks like this:

.. code-block:: ini

    [Path]
    dapsrc           = ${MANGADAP_DIR}

    [default]
    key              = MILESHC
    file_search      = ${Path:dapsrc}/mangadap/data/spectral_templates/miles_cluster/*.fits
    fwhm             = 2.50
    sres_ext
    in_vacuum        = False
    wave_limit       = 3575, 7400
    lower_flux_limit = 0.0
    log10            = False

This configuration (ini) file allows one to construct the ``MILESHC``
library as follows:

.. code-block:: python

    tpl = TemplateLibrary('MILESHC',
                          velscale_ratio=4,     # Set the pixel size to 1/4 the MaNGA step
                          spectral_step=1e-4,   # The MaNGA step is dlogLambda = 1e-4
                          log=True,             # Sample the templates logarithmically
                          hardcopy=False)       # Don't save a hardcopy of the library to disk

This is possible because the DAP reads all the ini files in the
``$MANGADAP_DIR/mangadap/config/spectral_templates`` directory
and selects the file with ``key = MILESHC``.  If you have a template
library that you want to use that's not provided by the DAP, you can
write a configuration file and put it in the appropriate directory such
that it can be selected by the key.

The template library key for use with the stellar-continuum fitting
is defined in a separate configuration file specific to the
stellar-continuum fitting class; these are kept in
``$MANGADAP_DIR/mangadap/config/stellar_continuum_modeling``. For
example, the ``MILESHCMPL9`` configuration file looks like this:

.. code-block:: ini

    [default]
    key                    = MILESHCMPL9
    fit_type               = stellar_kinematics
    fit_method             = ppxf
    fit_iter               = nonzero_templates
    reject_boxcar          = 100
    filter_boxcar
    filter_op
    filter_iter
    filter_degree
    filter_mdegree
    minimum_snr            = 1.0
    waverange
    artifact_mask          = BADSKY
    emission_line_mask     = ELPMPL8
    template_library       = MILESHC
    match_resolution       = False
    velscale_ratio         = 4
    moments                = 2
    degree                 = 8
    mdegree                = -1
    bias

You can see that the file defines ``template_library = MILESHC``.  To
execute the full DAP using a new template library is a matter of setting
up these configuration files.

However, you can also write scripts that incorporate the DAP
functionality without the need to add configuration files. Assume you
have a script that uses a
:class:`~mangadap.proc.templatelibary.TemplateLibrary` object, you can
define a new template library in the code itself using the
:class:`~mangadap.proc.templatelibrary.TemplateLibraryDef` object. The
:class:`~mangadap.proc.templatelibrary.TemplateLibraryDef` object is
actually the product of the parsed configuration file within the main
DAP code. For example:

.. code-block:: python

    # Imports
    from mangadap.proc.templatelibrary import TemplateLibraryDef, TemplateLibrary

    # Define the search string for the library
    search_str = '/path/to/library/*.fits'
    search_sres_str = '/path/to/library/with/sres/*.fits'

    # Define the template library parameters
    new_tpl_lst = TemplateLibraryDef(key='MYLIB',            # Unique keyword for the library
                                     file_search=search_str, # Search string
                                     fwhm=2.50,              # FWHM of resolution element
                                     in_vacuum=False,        # Wavelength in vacuum?
                                     wave_limit=numpy.array([ 3575., 7400. ]),   # Valid Range
                                     lower_flux_limit=0.0,   # Lower limit for valid flux
                                     log10=False)            # Log binned?

    # Or if you there is an extension SPECRES in *all* the files with
    # the spectral resolution:
    new_tpl_list = [new_tpl_list,
                    TemplateLibraryDef(key='MYLIB_SRES',       # Unique library keyword
                                       file_search=search_sres_str, # Search string
                                       sres_ext='SPECRES',     # Spectral Resolution Extension
                                       in_vacuum=False,        # Wavelength in vacuum?
                                       wave_limit=numpy.array([ 3575., 7400. ]),   # Valid range
                                       lower_flux_limit=0.0,   # Lower limit for valid flux
                                       log10=False)            # Log binned?
                   ]

    # Read and process template library
    tpl = TemplateLibrary('MYLIB',
                          tpllib_list=new_tpl_lst,  # Available list of template libraries
                          velscale_ratio=4,     # Set the pixel size to 1/4 the MaNGA step
                          spectral_step=1e-4,   # The MaNGA step is dlogLambda = 1e-4
                          log=True,             # Sample the templates logarithmically
                          hardcopy=False)       # Don't save a hardcopy of the library to disk


Adding new functionality
------------------------

The pairing of the defining parameters of a specific analysis method
and the instantiation of the method itself, like what is shown above,
is ubiquitous in the DAP, following from one of the main design
principles. Apart from allowing one to alter the details of how the
DAP proceeds via changing or adding configuration files (instead of
changing the code itself), this also facilitates incorporating new
algorithms within the existing infrastructure.

The implementation of hooks for including new algorithms into the DAP
infrastructure exists; however, it is minimal in some respects and
hasn't been well tested. That means that, for anyone that tries this,
there are sure to be some growing pains in making sure it works
properly. At the moment, a primary limitation to incorporating new
algorithms is that the procedure is not as simple as adding a new
file with code in the ``$MANGADAP_DIR/mangadap/contrib`` directory
and a new configuration file. Instead, one has to alter a number of
bits of code in the DAP. There are a few ways to do this, but testing
of these methods has been limited. The following are two sketched out
examples of including new functionality or algorithms.

1. Adding a new binning scheme
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The class that constructs the spatially binned spectra is
:class:`~mangadap.proc.spatiallybinnedspectra.SpatiallyBinnedSpectra`.
The method used to construct a class instance is defined using
:class:`~mangadap.proc.spatiallybinnedspectra.SpatiallyBinnedSpectraDef`,
which has components that define a set of binning parameters, a binning
class instance, and/or a binning function.  The existing binning schemes
are:

 * :class:`~mangadap.proc.spatialbinning.GlobalBinning`
 * :class:`~mangadap.proc.spatialbinning.RadialBinning`
 * :class:`~mangadap.proc.spatialbinning.VoronoiBinning`
 * :class:`~mangadap.proc.spatialbinning.SquareBinning`

All of these classes provide a common interface that
:class:`~mangadap.proc.spatiallybinnedspectra.SpatiallyBinnedSpectra`
calls to determine which spaxels are assigned to each bin.  The format
of this function must be:

.. code-block:: python

    def binning_function(x, y, par=None):
        # Bin the data
        ...
        return bin_id

That is, the function must take in the on-sky x and y positions of each
spaxel, accept some set of parameters provided by the ``par`` dictionary
and return a bin ID number associated with each x and y position.
 
So let's say that you wanted to bin all spectra in a set of apertures.
You could define and implement a function that performs this binning,
and then execute this binning approach within
:class:`~mangadap.proc.spatiallybinnedspectra.SpatiallyBinnedSpectra` as
follows (the code is untested!):

.. code-block:: python

    #!/usr/bin/env python3

    import time
    import warnings
    import numpy
    import astropy.constants

    from mangadap.datacube import MaNGADataCube
    from mangadap.proc.reductionassessments import ReductionAssessment
    from mangadap.proc.spectralstack import SpectralStackPar, SpectralStack
    from mangadap.proc.spatiallybinnedspectra import SpatiallyBinnedSpectra, SpatiallyBinnedSpectraDef
    from mangadap.proc.stellarcontinuummodel import StellarContinuumModel
    from mangadap.proc.emissionlinemoments import EmissionLineMoments
    from mangadap.proc.emissionlinemodel import EmissionLineModel
    from mangadap.proc.spectralindices import SpectralIndices
    from mangadap.dapfits import construct_maps_file, construct_cube_file

    #-----------------------------------------------------------------------------

    class ApertureBinning():
        """
        Perform aperture binning

        Args:
            x (array-like):
                List of on-sky x coordinates for apertures
            y (array-like):
                List of on-sky y coordinates for apertures
            r (array-like):
                Single or list of radii of the apertures

        Attributes:
            n (:obj:`int`):
                Number of apertures
            x (`numpy.ndarray`_):
                On-sky x coordinates for apertures
            y (`numpy.ndarray`_):
                On-sky y coordinates for apertures
            r (`numpy.ndarray`_):
                Aperture radii
        """
        def __init__(self, x, y, r):
            self.x = numpy.asarray(x)
            if len(self.x.shape) != 1:
                raise ValueError('On-sky coordinates must be one-dimensional.')
            self.n = x.size
            if len(y) != self.n:
                raise ValueError('Input coordinates are of different lengths.')
            self.y = numpy.asarray(y)
            if len(r) != 1 and len(r) != self.n:
                raise ValueError('Radii must be common to all apertures or unique to each aperture.')
            self.r = numpy.full(self.n, r, dtype=float) if len(r) == 1 else numpy.asarray(r)

        def bin_spaxels(self, x, y, par=None):
            _x = numpy.asarray(x)
            if len(_x.shape) != 1:
                raise ValueError('On-sky coordinates must be one-dimensional.')
            nspaxels = _x.size
            if len(y) != nspaxels:
                raise ValueError('Input coordinates are of different lengths.')
            _y = numpy.asarray(y)

            # Find which spaxels land in each aperture
            indx = numpy.square(_x[:,None]-self.x[None,:]) + numpy.square(_y[:,None]-self.y[None,:]) \
                        < numpy.square(self.r[None,:])
            if numpy.any(numpy.sum(indx, axis=1) > 1):
                warnings.warn('Spaxels found in multiple apertures!')

            # Return the aperture index that each spaxel is within,
            # isolating only one aperture per spaxel; spaxels not in any
            # aperture have a bin ID of -1
            binid = numpy.full((nspaxels, self.n), -1, dtype=int)
            binid[indx] = numpy.array([numpy.arange(self.n)]*nspaxels)[indx]
            return numpy.amax(binid, axis=1)

    #-----------------------------------------------------------------------------
    if __name__ == '__main__':
        t = time.perf_counter()

        # Set the plate, ifu, and initial velocity/redshift
        plate = 7495
        ifu = 12704
        vel = 8675.5
        nsa_redshift = vel/astropy.constants.c.to('km/s').value

        # Read the DRP LOGCUBE file
        cube = MaNGADataCube.from_plateifu(plate, ifu)

        # Calculate the S/N and coordinates
        rdxqa = ReductionAssessment('SNRG', cube)

        # Setup the aperture binning class
        ax = numpy.array([0.0, 3.0, 6.0])
        ay = numpy.array([0.0, 0.0, 0.0])
        apbin = ApertureBinning(ax, ay, 2.5)

        # Setup the stacking operations
        stackpar = SpectralStackPar('mean',         # Operation for stack
                                    False,          # Apply a velocity registration
                                    None,           # Velocity offsets for registration
                                    'channels',     # Covariance mode and parameters
                                    SpectralStack.parse_covariance_parameters('channels', 11))
        stacker = SpectralStack()

        # Create a new binning method
        binning_method = SpatiallyBinnedSpectraDef('Aperture',      # Key for binning method
                                                   'ODonnell',      # Galactic reddening function
                                                   3.1,             # Rv for Galactic reddening
                                                   0.0,             # Minimum S/N to include
                                                   None,            # Object with binning pars
                                                   None,            # Binning class instance
                                                   apbin.bin_spaxels,   # Binning function
                                                   stackpar,        # Object with stacking pars
                                                   stacker,         # Stacking class instance
                                                   stacker.stack_datacube)   # Stacking function

        # Bin the spectra using the new binning method
        binned_spectra = SpatiallyBinnedSpectra('Aperture',     # Key for binning method
                                                cube,           # DRP data to bin
                                                rdxqa,          # Cube coordinates and S/N
                                                method_list=binning_method) # Binning methods

        # The rest of this is just a single execution of the remaining
        # analysis steps in
        # $MANGADAP_DIR/python/mangadap/survey/manga_dap.py , with some
        # simplifications
        stellar_continuum = StellarContinuumModel('GAU-MILESHC', binned_spectra, guess_vel=vel,
                                                  guess_sig=100.)

        emission_line_moments = EmissionLineMoments('EMOMF', binned_spectra,
                                                    stellar_continuum=stellar_continuum,
                                                    redshift=nsa_redshift)

        emission_line_model = EmissionLineModel('EFITF', binned_spectra,
                                                stellar_continuum=stellar_continuum,
                                                redshift=nsa_redshift, dispersion=100.0)
        
        spectral_indices = SpectralIndices('INDXEN', binned_spectra, redshift=nsa_redshift,
                                           stellar_continuum=stellar_continuum,
                                           emission_line_model=emission_line_model)

        construct_maps_file(cube, rdxqa=rdxqa, binned_spectra=binned_spectra,
                            stellar_continuum=stellar_continuum,
                            emission_line_moments=emission_line_moments,
                            emission_line_model=emission_line_model,
                            spectral_indices=spectral_indices, nsa_redshift=nsa_redshift)

        construct_cube_file(cube, binned_spectra=binned_spectra,
                            stellar_continuum=stellar_continuum,
                            emission_line_model=emission_line_model)

        print('Elapsed time: {0} seconds'.format(time.perf_counter() - t))

You'll notice that there are some limitations in what one can implement.
In this example, the limitation is that each spaxel must be assigned to
a single unique bin, not multiple bins if the apertures are overlapping.

2. Adding a new emission-line fitter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Adding a new binning scheme is relatively straight-forward because all
that's required is to provide a new binning function that adheres to the
specified form.  Things become much more complicated when you want to
replace a core algorithm that provides much of the content of the output
data model.  Still, it can be done, it's just that one has to follow
more requirements that can be more stringent.

Let's say you want to add a new emission-line fitter (as we did for
MPL-6 in changing from :class:`~mangadap.proc.elric.Elric` to
:class:`~mangadap.proc.sasuke.Sasuke`).  The class that constructs the
parameterized emission-line models is
:class:`~mangadap.proc.emissionlinemodel.EmissionLineModel`.  The method
used to construct a class instance is defined using
:class:`~mangadap.proc.emissionlinemodel.EmissionLineModelDef`, which has
components that define a set of model-fitting parameters, a
model-fitting class instance, and/or a model-fitting function.  The
common function call that any emission-line fitter must provide looks
like:

.. code-block:: python

    def fit(binned_spectra, par=None, loggers=None, quiet=False):
        # Fit the spectra
        ...
        return model_eml_flux, model_eml_base, model_eml_mask, model_fit_par, \
                model_eml_par, model_binid

where ``binned_spectra`` is a
:class:`~mangadap.proc.spatiallybinnedspectra.SpatiallyBinnedSpectra`
object and the returned arrays are:

 - ``model_eml_flux``: Model emission-line flux only; shape is
   :math:`(N_{\rm mod}, N_{\rm wave})`.  The first axis is ordered by model
   ID number.

 - ``model_eml_base``: Any baseline resulting from the emission-line fit
   such that the model fit to each spectrum is: ``stellar_continuum +
   model_eml_flux + model_eml_base``; shape is :math:`(N_{\rm mod},
   N_{\rm wave})`.  The first axis is ordered by model ID number.

 - ``model_eml_mask``: Boolean or bit-mask array for fitted models;
   shape is :math:`(N_{\rm mod}, N_{\rm wave})`.  The first axis is
   ordered by model ID number.

 - ``model_fit_par``: A ``numpy`` record array that provides the results
   of each fit.  This can be ``None`` and the output maps file will
   still be successfully written

 - ``model_eml_par``: A ``numpy`` record array that provides the output
   model parameters.  The data type must be as returned by
   :func:`mangadap.proc.spectralfitting.EmissionLineFit._per_emission_line_dtype`
   and the shape must be :math:`(N_{\rm mod},)` with the parameters
   ordered by the model ID number.

 - ``model_binid``: A 2D map of the ID numbers assigned to each spaxel
   with a fitted model.  Any spaxel without an emission-line model
   should have ``model_binid = -1``, and the number of IDs that are
   greater than -1 must be :math:`N_{\rm mod}`.  The shape is
   :math:`(N_x,N_y)`; i.e., it must match the spatial dimensions of the
   fitted DRP data cube.  This can be {{{None}}}, which indicates that
   the bin IDS are the same as the bin IDs set in the ``binned_spectra``
   object.

In addition to the high-level data model interfaces, like
:class:`~mangadap.proc.emissionlinemodel.EmissionLineModel`, the DAP
attempts to provide a set of base classes that provided functionality
common to a set of abstracted spectral-fitting routines.  For the
emission-line fitting, this is the
:class:`~mangadap.proc.spectralfitting.EmissionLineFit` class.  This
object provides the data table description that should be common to all
emission-line model output for the construction of the output data model
by :class:`~mangadap.proc.emissionlinemodel.EmissionLineModel` (see
``model_eml_par`` above).

.. warning::

    The code discussed/provided below was an initial go at sketching out
    the code for the emission-line module used in MPL-6.  The code will
    not work out of the box and is meant to illustrate the solution to
    the problem.  For the actual solution, see the main DAP interface
    class :class:`~mangadap.proc.sasuke.Sasuke` and the primary fitting
    function written by Xihan Ji and Michele Cappellari (with some
    significant edits by Kyle Westfall), :mod:`mangadap.contrib.xjmc`.

So, let's say you have a function ``first_second_iteration`` that
applies the new fitting approach that we want for the emission lines.
The following is untested code that could be used to implement this
function as the emission-line model fitter while still providing the
same DAP output data model:

.. code-block:: python

    #!/usr/bin/env python3

    import time
    import warnings
    import numpy
    import astropy.constants

    from mangadap.datacube import MaNGADataCube
    from mangadap.proc.reductionassessments import ReductionAssessment
    from mangadap.proc.spatiallybinnedspectra import SpatiallyBinnedSpectra
    from mangadap.proc.stellarcontinuummodel import StellarContinuumModel
    from mangadap.proc.ppxffit import PPXFFit
    from mangadap.util.instrument import spectrum_velocity_scale
    from mangadap.util.fitsutil import DAPFitsUtil
    from mangadap.util.fileio import init_record_array
    from mangadap.proc.emissionlinemoments import EmissionLineMoments
    from mangadap.proc.spectralfitting import EmissioneLineFit
    from mangadap.proc.emissionlinemodel import EmissionLineModelDef, EmissionLineModel
    from mangadap.par.emissionlinedb import EmissionLineDB
    from mangadap.proc.spectralindices import SpectralIndices
    from mangadap.dapfits import construct_maps_file, construct_cube_file

    from mangadap.contrib.xjmc import first_second_iteration

    #-----------------------------------------------------------------------------

    class XJMCEmissionLineFitter(EmissionLineFit):
        def __init__(self, par=None):
            if par is None:
                # Set the default parameter set.  The guess_redshift,
                # stellar_continuum, and emission_lines values can be
                # filled by EmissionLineModel._fill_method_par()
                par = { 'guess_redshift': None,     # The guess redshift for each binned spectrum
                        'stellar_continuum': None,  # The StellarContinuumModel object
                        'emission_lines': None,     # The EmissionLineDB object
                        'degree': 8,                # Additive polynomial order
                        'mdegree': 0 }              # Multiplicative polynomial order
            EmissionLineFit.__init__(self, 'XJMC', None, par=par)


        def fit(self, binned_spectra, par=None, loggers=None, quiet=False):
            if par is not None:
                self.par = par

            # Check the parameter keys
            required_keys = [ 'guess_redshift', 'stellar_continuum', 'emission_lines', 'degree',
                              'mdegree' ]
            if numpy.any([ reqk not in self.par.keys() for reqk in required_keys ]):
                raise ValueError('Parameter dictionary does not have all the required keys.')

            # Wavelengths are in vacuum
            wave = binned_spectra['WAVE'].data.copy()
            # Velocity step per pixel
            velscale = spectrum_velocity_scale(wave)
            # Flux and noise masked arrays; shape is (Nspaxels,Nwave)
            # where Nspaxels is Nx*Ny
            flux = binned_spectra.cube.copy_to_masked_array(flag=['DONOTUSE', 'FORESTAR'])
            noise = numpy.ma.power(binned_spectra.cube.copy_to_masked_array(
                                        attr='ivar', flag=['DONOTUSE', 'FORESTAR']), -0.5)
            # Spaxel coordinates; shape is (nspaxels,)
            x = binned_spectra.rdxqa['SPECTRUM'].data['SKY_COO'][:,0]
            y = binned_spectra.rdxqa['SPECTRUM'].data['SKY_COO'][:,1]
            # Binned flux and binned noise masked arrays; shape is (nbins,nwave)
            flux_binned = binned_spectra.copy_to_masked_array(
                                                        flag=binned_spectra.do_not_fit_flags())
            noise_binned = numpy.ma.power(binned_spectra.copy_to_masked_array(ext='IVAR',
                                                    flag=binned_spectra.do_not_fit_flags()) , -0.5)
            # Bin coordinates; shape is (nbins,)
            x_binned = binned_spectra['BINS'].data['SKY_COO'][:,0]
            y_binned = binned_spectra['BINS'].data['SKY_COO'][:,1]
            # Set initial guesses for the velocity and velocity
            # dispersion
            if self.par['guess_redshift'] is not None:
                # Use guess_redshift if provided
                vel = self.par['guess_redshift'] * astropy.constants.c.to('km/s').value
                # And set default velocity dispersion to 100 km/s
                sig = numpy.full(vel.size, 100, dtype=float)
            elif self.par['stellar_continuum'] is not None:
                # Otherwise use the stellar-continuum result
                vel, sig = self.par['stellar_continuum'].matched_guess_kinematics(binned_spectra,
                                                                                  cz=True)
            else:
                # TODO: Set default guess kinematics
                vel, sig = None, None

            # Get the stellar templates;
            # shape is (Ntemplates, Nwave_templates)
            if self.par['stellar_continuum'] is not None:
                stars_templates = self.par['stellar_continuum'].method['fitpar']['template_library']
                stars_templates_wave = stars_templates['WAVE'].data.copy()
                stars_templates = stars_templates['FLUX'].data.copy()
                velscale_ratio = self.par['stellar_continuum'].method['fitpar']['velscale_ratio']
                dv = -PPXFFit.ppxf_tpl_obj_voff(stars_templates_wave, wave, velscale,
                                                velscale_ratio=velscale_ratio)
            else:
                # TODO: Default construction of stellar templates
                stars_templates = None
                velscale_ratio = None
                dv = None

            # TODO: Construct gas templates
            gas_templates = None
            gas_names = None

            # TODO: Default polynomial orders
            degree = 8 if self.par['degree'] is None else self.par['degree']
            mdegree = 0 if self.par['mdegree'] is None else self.par['mdegree']
        
            # Output is:
            #   - model_flux: stellar-continuum + emission-line model;
            #     shape is (Nmod, Nwave); first axis is ordered by model
            #     ID number
            #   - model_eml_flux: model emission-line flux only; shape
            #     is (Nmod, Nwave); first axis is ordered by model ID
            #     number
            #   - model_mask: boolean or bit mask for fitted models;
            #     shape is (Nmod, Nwave); first axis is ordered by model
            #     ID number
            #   - model_binid: ID numbers assigned to each spaxel with a
            #     fitted model; any spaxel without a model should have
            #     model_binid = -1; the number of >-1 IDs must be Nmod;
            #     shape is (Nx,Ny) which is equivalent to:
            #       flux[:,0].reshape((numpy.sqrt(Nspaxels).astype(int),)*2).shape
            #   - eml_flux: Flux of each emission line; shape is
            #     (Nmod,Neml)
            #   - eml_fluxerr: Error in emission-line fluxes; shape is
            #     (Nmod, Neml)
            #   - eml_kin: Kinematics (velocity and velocity dispersion)
            #     of each emission line; shape is (Nmod,Neml,Nkin)
            #   - eml_kinerr: Error in the kinematics of each emission
            #     line
            #   - eml_sigmacorr: Quadrature corrections required to
            #     obtain the astrophysical velocity dispersion; shape is
            #     (Nmod,Neml); corrections are expected to be applied as
            #     follows:
            #       sigma = numpy.ma.sqrt( numpy.square(eml_kin[:,:,1])
            #                               - numpy.square(eml_sigmacorr))
            model_flux, model_eml_flux, model_mask, model_binid, eml_flux, eml_fluxerr, \
                    eml_kin, eml_kinerr, eml_sigmacorr \
                            = first_second_iteration(wave, flux, noise, flux_binned, noise_binned,
                                                     velscale, velscale_ratio, dv, vel, sig,
                                                     stars_templates, gas_templates, gas_names,
                                                     degree, mdegree, x, y, x_binned, y_binned)

            # The ordered indices in the flatted bin ID map with/for
            # each model
            model_srt = numpy.argsort(model_binid.ravel())[model_binid.ravel() > -1]

            # Construct the output emission-line database.  The data
            # type defined by
            # EmissionLineFit._per_emission_line_dtype(); shape is
            # (Nmod,); parameters must be ordered by model ID number
            nmod = len(model_srt)
            neml = eml_flux.shape[1]
            nkin = eml_kin.shape[-1]
            model_eml_par = init_record_array(nmod,
                                EmissionLineFit._per_emission_line_dtype(neml, nkin, numpy.int16))
            model_eml_par['BINID'] = model_binid.ravel()[model_srt]
            model_eml_par['BINID_INDEX'] = numpy.arange(nmod)
            model_eml_par['MASK'][:,:] = 0
            model_eml_par['FLUX'] = eml_flux
            model_eml_par['FLUXERR'] = eml_fluxerr
            model_eml_par['KIN'] = eml_kin
            model_eml_par['KINERR'] = eml_kinerr
            model_eml_par['SIGMACORR'] = eml_sigmacorr

            # Include the equivalent width measurements
            if self.par['emission_lines'] is not None:
                EmissionLineFit.measure_equivalent_width(wave, flux[model_srt,:],
                                                         par['emission_lines'], model_eml_par)

            # Calculate the "emission-line baseline" as the difference
            # between the stellar continuum model determined for the
            # kinematics and the one determined by the optimized
            # stellar-continuum + emission-line fit:
            if self.par['stellar_continuum'] is not None:
                # Construct the full 3D cube for the stellar continuum
                # models
                sc_model_flux, sc_model_mask \
                        = DAPFitsUtil.reconstruct_cube(binned_spectra.shape,
                                                       self.par['stellar_continuum']['BINID'].data,
                                                   [ self.par['stellar_continuum']['FLUX'].data,
                                                     self.par['stellar_continuum']['MASK'].data ])
                # Set any masked pixels to 0
                sc_model_flux[sc_model_mask>0] = 0.0

                # Construct the full 3D cube of the new stellar
                # continuum from the combined stellar-continuum +
                # emission-line fit
                el_continuum = DAPFitsUtil.reconstruct_cube(binned_spectra.shape, model_binid,
                                                            model_flux - model_eml_flux)
                # Get the difference, restructure it to match the shape
                # of the emission-line models, and zero any masked
                # pixels
                model_eml_base = (el_model_flux - sc_model_flux).reshape(-1,wave.size)[model_srt,:]
                if model_mask is not None:
                    model_eml_base[model_mask>0] = 0.0
            else:
                model_eml_base = numpy.zeros(model_flux.shape, dtype=float)

            # Returned arrays are:
            #   - model_eml_flux: model emission-line flux only; shape
            #     is (Nmod, Nwave); first axis is ordered by model ID
            #     number
            #   - model_eml_base: difference between the combined fit
            #     and the stars-only fit; shape is (Nmod, Nwave); first
            #     axis is ordered by model ID number
            #   - model_mask: boolean or bit mask for fitted models;
            #     shape is (Nmod, Nwave); first axis is ordered by model
            #     ID number
            #   - model_fit_par: This provides the results of each fit;
            #     TODO: The is set to None.  Provide metrics of the ppxf
            #     fit to each spectrum?
            #   - model_eml_par: output model parameters; data type must
            #     be EmissionLineFit._per_emission_line_dtype(); shape
            #     is (Nmod,); parameters must be ordered by model ID
            #     number
            #   - model_binid: ID numbers assigned to each spaxel with a
            #     fitted model; any spaxel with a model should have
            #     model_binid = -1; the number of >-1 IDs must be Nmod;
            #     shape is (Nx,Ny)
            return model_eml_flux, model_eml_base, model_mask, None, model_eml_par, model_binid

    
    #-----------------------------------------------------------------------------
    if __name__ == '__main__':
        t = time.perf_counter()

        # Set the plate, ifu, and initial velocity/redshift
        plate = 7495
        ifu = 12704
        vel = 8675.5
        nsa_redshift = vel/astropy.constants.c.to('km/s').value

        # Read the DRP LOGCUBE file
        cube = MaNGADataCube.from_plateifu(plate, ifu)

        # Calculate the S/N and coordinates
        rdxqa = ReductionAssessment('SNRG', cube)

        # Peform the Voronoi binning to S/N>~10
        binned_spectra = SpatiallyBinnedSpectra('VOR10', cube, rdxqa)

        # Fit the stellar kinematics
        stellar_continuum = StellarContinuumModel('GAU-MILESHC', binned_spectra, guess_vel=vel,
                                                  guess_sig=100.)

        # Get the emission-line moments
        emission_line_moments = EmissionLineMoments('EMOMF', binned_spectra,
                                                    stellar_continuum=stellar_continuum,
                                                    redshift=nsa_redshift)

        # Get an estimate of the redshift of each bin using the first
        # moment of the H-alpha emission line:
        el_init_redshift = numpy.full(binned_spectra.nbins, nsa_redshift, dtype=float)
        # HARDCODED FOR A SPECIFIC EMISSION-LINE MOMENT DATABASE
        halpha_channel = 7
        halpha_mom1_masked = emission_line_moments['ELMMNTS'].data['MASK'][:,halpha_channel] > 0
        # - Use the 1st moment of the H-alpha line
        el_init_redshift[ emission_line_moments['ELMMNTS'].data['BINID_INDEX'] ] \
                    = emission_line_moments['ELMMNTS'].data['MOM1'][:,halpha_channel] \
                                    / astropy.constants.c.to('km/s').value
        # - For missing bins in the moment measurements and bad H-alpha
        #   moment measurements, use the value for the nearest good bin
        bad_bins = numpy.append(emission_line_moments.missing_bins,
                    emission_line_moments['ELMMNTS'].data['BINID'][halpha_mom1_masked]).astype(int)
        if len(bad_bins) > 0:
            nearest_good_bin_index = binned_spectra.find_nearest_bin(bad_bins, indices=True)
            bad_bin_index = binned_spectra.get_bin_indices(bad_bins)
            el_init_redshift[bad_bin_index] = el_init_redshift[nearest_good_bin_index]

        # Setup the new emission-line fitter
        fitter = XJMCEmissionLineFitter()

        # Setup the new fitting method
        fit_method = EmissionLineModelDef('XJMC',       # Key for the fitting method
                                          0.0,          # Minimum S/N of the binned spectra
                                          None,         # Keyword for an artifact mask
                                          None,         # Keyword for an emission-line database
                                          fitter.par,   # Object with fit parameters
                                          fitter,       # Fitting class instance
                                          fitter.fit)   # Fitting function

        # Fit the emission lines
        emission_line_model = EmissionLineModel('XJMC',
                                                binned_spectra,
                                                stellar_continuum=stellar_continuum,
                                                redshift=el_init_redshift, dispersion=100.0,
                                                method_list=fit_method)

        # The rest of this is just a single execution of the remaining
        # analysis steps in
        # $MANGADAP_DIR/python/mangadap/survey/manga_dap.py , with some
        # simplifications
        spectral_indices = SpectralIndices('INDXEN', binned_spectra, redshift=nsa_redshift,
                                           stellar_continuum=stellar_continuum,
                                           emission_line_model=emission_line_model)

        construct_maps_file(cube, rdxqa=rdxqa, binned_spectra=binned_spectra,
                            stellar_continuum=stellar_continuum,
                            emission_line_moments=emission_line_moments,
                            emission_line_model=emission_line_model,
                            spectral_indices=spectral_indices, nsa_redshift=nsa_redshift)

        construct_cube_file(cube, binned_spectra=binned_spectra,
                            stellar_continuum=stellar_continuum,
                            emission_line_model=emission_line_model)

        print('Elapsed time: {0} seconds'.format(time.perf_counter() - t))


