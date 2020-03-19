
.. _emission-line-modeling:

Emission-line Modeling
======================

*Analysis class*: :class:`~mangadap.proc.emissionlinemodel.EmissionLineModel`

*Method definition*: :class:`~mangadap.proc.emissionlinemodel.EmissionLineModelDef`

*Reference root*: ``$MANGA_SPECTRO_ANALYSIS/$MANGADRP_VER/$MANGADAP_VER/[PLATE]/[IFUDESIGN]/ref``

*Reference file*: ``manga-[PLATE]-[IFUDESIGN]-[DRPQA_KEY]-[BIN_KEY]-[CONTINUUM_KEY]-[ELFIT_KEY].fits.gz``

*Config files*: ``$MANGADAP_DIR/mangadap/config/emission_line_modeling``

*Example config*: ``efitmpl9.ini``

.. code-block:: ini

    [default]
    key                  = EFITMPL9
    fit_method           = sasuke
    minimum_snr          = 0.0
    deconstruct_bins     = ignore
    mom_vel_name         = Ha-6564
    mom_disp_name
    waverange
    artifact_mask        = BADSKY
    emission_lines       = ELPMPL9
    etpl_line_sigma_mode = offset
    etpl_line_sigma_min  = 10
    baseline_order
    window_buffer
    reject_boxcar        = 101
    continuum_templates  = MASTARHC
    velscale_ratio       = 1
    bias
    moments
    degree
    mdegree              = 14
    internal_reddening

*Important class dependencies*:

 * :class:`~mangadap.util.pixelmask.SpectralPixelMask`: Used to mask
   spectral regions.
 * :class:`~mangadap.par.emissionlinedb.EmissionLineDB`: Generalized
   class that provides the detailed parameters for a set of emission
   lines to fit
 * :class:`~mangadap.util.lineprofiles`: Generalized class that provides
   the detailed profile shapes for each line.
 * :class:`~mangadap.proc.emissionlinetemplates.EmissionLineTemplates`:
   Used to construct the emission-line templates to use with pPXF.
 * :class:`~mangadap.proc.spectralfitting.EmissionLineFit`: Provides
   functions common to both the moment and Gaussian-fit calculations.
 * :class:`~mangadap.proc.sasuke.Sasuke`: The equivalent of
   :class:`~mangadap.proc.ppxffit.PPXFFit` for the emission-line fits
   that is primarily a wrapper for pPXF.  Critically, this uses the
   emission-line fitter contributed by Xihan Ji and Michele Cappellari,
   :func:`~mangadap.contrib.xjmc.emline_fitter_with_ppxf`.
 * :class:`~mangadap.proc.templatelibrary.TemplateLibrary`: If using
   :class:`~mangadap.proc.sasuke.Sasuke` to perform the fit, this
   generalized interface provides the spectral-template library for use
   in modeling the stellar continuum.
 * :class:`~mangadap.proc.elric.Elric`:  Provides the main fitting
   functions when just fitting Gaussians to continuum-subtracted
   spectra.  **BEWARE**: This class has not been used or tested
   regularly since MPL-5.

*Algorithm*:

 * Setup the fitting method:

    * Instantiate the
      :class:`~mangadap.util.pixelmask.SpectralPixelMask` using the
      ``artifact_mask`` and ``waverange`` from config.

        * The ``BADSKY`` artifact mask is read and used to build an
          :class:`~mangadap.par.artifactdb.ArtifactDB` instance that
          masks the typical residuals around the strong sky line at
          5577 angstroms.
        * The ``waverange`` config parameter can be used to limit the
          fitted spectral range; will fit as much as possible if no
          range is provided.

    * If ``fit_method = elric``, ``baseline_order`` sets the Legendre
      function used to set the baseline in each fitting window and
      ``window_buffer`` sets the +/- window in angstroms around each
      line to use during the fit.
    * If ``fit_method = sasuke``:

        * ``etpl_line_sigma_mode`` and ``etpl_line_sigma_min``
          determines the method used to set the emission-line
          template instrumental dispersion; the available options are
          set by
          :func:`~mangadap.proc.sasuke.Sasuke.etpl_line_sigma_options`.
        * ``reject_boxcar`` sets the size of the boxcar (pixels) to
          use for the rejection iterations.
        * The templates used to fit the stellar continuum during the
          emission-line modeling can be different than those used
          during the stellar kinematics fit. Use
          ``continuum_templates`` and ``velscale_ratio`` to select
          the new templates and set their sampling. If
          ``continuum_templates`` is not given ``velscale_ratio`` is
          ignored and the templates are identical between the two
          modules. If the templates are switched, a new
          :class:`~mangadap.proc.templatelibrary.TemplateLibrary`
          object is instantiated. When switching template libraries,
          the templates **must** have their resolution matched to the
          MaNGA data so that the corrected stellar kinematics from
          the existing
          :class:`~mangadap.proc.stellarcontinuummodel.StellarContinuumModel`
          instance can be held fixed during the fitting.
        * ``bias``, ``degree``, ``mdegree`` are passed directly to
          pPXF (at the moment ``moments`` is ignored and always 2!)
        * ``internal_reddening = True`` forces use of the Calzetti
          (2000) attenuation law; will override any non-zero
          ``mdegree``.

    * If ``deconstruct_bins = True``, method will fit the emission
      lines on an individual spaxel basis instead of the binned
      spectra; currently this can only be used if ``fit_method =
      sasuke``.
    * If ``mom_vel_name`` or ``mom_disp_name`` is defined, the DAP
      will use the corresponding moment measurements from the
      :class:`~mangadap.proc.emissionlinemoments.EmissionLineMoments`
      object to set the initial guess for the velocity (default is
      the datacube metadata redshift) and/or velocity dispersion
      (default is 100 km/s) for the fit.

 * If requested, call the "Elric" fitter:

    * **WARNING**: The Elric fitter has not been used in the DAP for
      some time. It should generally not be selected; if it is, one
      may need to spend some time debugging... For this reason, the
      method is not well documented here. Its main DAP wrapper
      fitting function is
      :func:`~mangadap.proc.elric.Elric.fit_SpatiallyBinnedSpectra`
      and its generalized fitting function is
      :func:`~mangadap.proc.elric.Elric.fit`.

 * Or, call "Sasuke" fitter:

    * The main DAP wrapper function is
      :func:`~mangadap.proc.sasuke.Sasuke.fit_SpatiallyBinnedSpectra`
      and does the following:

        * Get the binned spectra from the
          :class:`~mangadap.proc.spatiallybinnedspectra.SpatiallyBinnedSpectra`
          object
        * Either get the stellar templates from the
          :class:`~mangadap.proc.templatelibrary.TemplateLibrary`
          object pointed to by the
          :class:`~mangadap.proc.stellarcontinuummodel.StellarContinuumModel`
          object or, if new templates were selected, build the new
          :class:`~mangadap.proc.templatelibrary.TemplateLibrary`
          instance (which **must** have its resolution matched to the
          MaNGA data).

        * Get the fitted stellar kinematics from the
          :class:`~mangadap.proc.stellarcontinuummodel.StellarContinuumModel`
          object using
          :func:`~mangadap.proc.stellarcontinuummodel.StellarContinuumModel.matched_kinematics`.

        * Determine which binned spectra have the ``minimum_snr``
          from config, and have a good continuum model (cannot be
          flagged as ``NOVALUE`` or ``FITFAILED``).
          
        * If deconstructing bins:

            * Get the individual spaxel spectra from the
              :class:`~mangadap.datacube.datacube.DataCube` object
            * Apply the reddening defined in the
              :class:`~mangadap.proc.spatiallybinnedspectra.SpatiallyBinnedSpectra`
              object
            * Get the individual on-sky spaxel coordinates from the
              :class:`~mangadap.proc.reductionassessments.ReductionAssessment`
              object and the unweighted on-sky binned-spectra
              coordinates from the
              :class:`~mangadap.proc.spatiallybinnedspectra.SpatiallyBinnedSpectra`
              object.
            * Run the generalized fitting function (see description
              below), providing the spectra to which the
              stellar-continuum results are "remapped" to for fitting
              the emission lines.
            * Measure the equivalent widths for the individual spaxels
              using
              :func:`~mangadap.proc.spectralfitting.EmissionLineFit.measure_equivalent_width`.

        * Otherwise:

            * Run the generalized fitting function (see description
              below), only providing the binned spectra.
            * Measure the equivalent widths for the binned spectra
              using
              :func:`~mangadap.proc.spectralfitting.EmissionLineFit.measure_equivalent_width`.

        * Construct the "emission-line baseline" as the difference
          between continuum+emission-line optimized fit from Sasuke
          and the stellar continuum fit for the stellar kinematics
          (from the
          :class:`~mangadap.proc.stellarcontinuummodel.StellarContinuumModel`
          object)

    * The generalized fitter is
      :func:`mangadap.proc.sasuke.Sasuke.fit` and initially proceeds
      very similarly to :func:`mangadap.proc.ppxffit.PPXFFit.fit`:

        * Check the input spectra to fit, guess kinematics, and
          remapping coordinates if provided
        * Check and set stellar templates and stellar kinematics if
          provided
        * Determine the spectral resolution to use for the emission-line
          templates; available options are set by
          :func:`~mangadap.proc.sasuke.Sasuke.etpl_line_sigma_options`.
        * Construct and add emission-line templates using
          :class:`~mangadap.proc.emissionlinetemplates.EmissionLineTemplates`.

            * Parse the
              :class:`~mangadap.par.emissionlinedb.EmissionLineDB`
              object:  Determine which lines to fit and how to group
              lines into the same template (flux ratio fixed and same
              kinematics) kinematic components (same velocity and
              velocity dispersion), velocity groups, and sigma groups
              using
              :func:`~mangadap.proc.emissionlinetemplates.EmissionLineTemplates._parse_emission_line_database`.
            * Sample the desired spectral resolution at each input line
              center.
            * Convert the profile parameters into pixel coordinates.
            * Construct each template with the specified line profile
              using classes/methods in
              :mod:`mangadap.util.lineprofiles`.  Lines with a fixed
              flux ratio are placed in the same template (this means
              they'll also have tied velocities and velocity
              dispersions).

        * Parse the velocity and sigma groups into tied parameters to
          provide to pPXF.
        * Given the template and object spectral range, determine the
          maximum viable fitting range for pPXF using
          :func:`~mangadap.proc.ppxffit.PPXFFit.fitting_mask`.
        * Run fit iterations using
          :func:`~mangadap.contrib.xjmc.emline_fitter_with_ppxf`.

            * If deconstructing the bins (for the
              :ref:`datamodel-hybrid-binning`):

                * (a.) First fit the binned spectra (with a 3-sigma
                  rejection iteration) forcing all the gas components
                  into a single kinematic component (all velocities and
                  velocity dispersions are tied).
                * (b.) Deconstruct the binned spectrum into its
                  individual spaxels.
                * (c.) Use the result of the first fit to create a
                  single, optimal stellar template and to set the
                  starting kinematics for first fit to each spaxel.
                * (d.) Fit each spaxel (with a 3-sigma rejection
                  iteration) with the optimized template and all gas
                  components in a single kinematic component
                * (e.) Reset the starting guesses and refit each spaxel
                  (*without* a rejection iteration) with the gas
                  components in the appropriate velocity and sigma
                  groups.

            * Otherwise:

                * Perform steps a, c, and e above, but just using the
                  provided spectra (whether or not they're bins or
                  individual spaxels).

        * Parse the results of the fit iterations into the output using
          :func:`~mangadap.proc.sasuke.Sasuke._save_results`.

