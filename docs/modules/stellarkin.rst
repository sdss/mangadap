
.. _stellar-kinematics:

Stellar Kinematics
==================

*Analysis class*: :class:`~mangadap.proc.stellarcontinuummodel.StellarContinuumModel`

*Reference root*: see :class:`~mangadap.config.analysisplan.AnalysisPlan.common_path`;
``$MANGA_SPECTRO_ANALYSIS/$MANGADRP_VER/$MANGADAP_VER/common/[PLATE]/[IFUDESIGN]`` for MaNGA

*Reference file*: see :class:`~mangadap.proc.stellarcontinuummodel.StellarContinuumModel.default_paths`
``manga-[PLATE]-[IFUDESIGN]-[RDXQA]-[BINNING]-[CONTINUUM].fits.gz`` for MaNGA

*Optional Parameters*: see :ref:`plan`.  The table below lists the parameters
defined by :class:`~mangadap.proc.stellarcontinuummodel.StellarContinuumModelDef`

.. include:: ../tables/reductionassessmentdef.rst

*Important class dependencies*:

 * :class:`~mangadap.proc.spectralfitting.StellarKinematicsFit`:
   Provides the abstracted base class for stellar kinematics fits.
 * :class:`~mangadap.util.pixelmask.SpectralPixelMask`: Used to mask the
   fitting regions for the fit.
 * :class:`~mangadap.proc.templatelibrary.TemplateLibrary`: Generalized
   interface that provides the spectral-template library for use in
   spectral fitting.
 * :class:`~mangadap.proc.ppxffit.PPXFFit`: Provides the core pPXF
   wrapper.

*Algorithm*:

 * Establish the fitting method, which includes setting up the
   :class:`~mangadap.util.pixelmask.SpectralPixelMask` using the
   ``artifact_mask``, ``emission_line_mask``, and ``waverange`` from
   config.

    * The ``BADSKY`` artifact mask is read and used to build an
      :class:`~mangadap.par.artifactdb.ArtifactDB` instance that masks
      the typical residuals around the strong sky line at 5577
      angstroms.
    * The ``ELPFULL`` emission-line mask is read and used to build an
      :class:`~mangadap.par.emissionlinedb.EmissionLineDB` instance that
      is used to mask the emission lines in the associated parameter
      file.  Each mask is 1500 km/s centered on the line at the input
      guess redshift.
    * The ``waverange`` config parameter can be used to limit the fitted
      spectral range; will fit as much as possible if no range is
      provided.

 * Using the
   :class:`~mangadap.proc.spatiallybinnedspectra.SpatiallyBinnedSpectra`
   object, select all binned spectra with S/N greater than
   ``minimum_snr`` in config.
 * The DAP nominally provides the stellar-continuum fit with the
   velocity and velocity dispersion from :ref:`execution-config` as
   its initial guess redshift and velocity dispersion.
 * Instantiate the
   :class:`~mangadap.proc.templatelibrary.TemplateLibrary` objects as
   selected by the ``template_library`` config parameter.

    * If matching the spectral resolution to the galaxy data
      (``match_resolution`` in config), the resolution is matched at the
      redshifted wavelengths of the galaxy data, adopting the input
      guess velocity as the redshift.
    * The template wavelength channel width is set to a fraction
      (1/``velscale_ratio``) of the galaxy data.

 * Execute the ``fit_method`` selected in config.  Currently, this can
   only be ``ppxf``.
 * In :func:`~mangadap.proc.ppxffit.PPXFFit.fit_SpatiallyBinnedSpectra`:

    * Mask binned spectra, ignoring pixels masked as ``DONOTUSE``,
      ``IGNORED``, ``FLUXINVALID``, or ``FORESTAR`` in DAP
      ``LOGCUBE`` file.
    * Call :func:`~mangadap.proc.ppxffit.PPXFFit.fit` with the data from
      the :class:`~mangadap.proc.templatelibrary.TemplateLibrary` and
      :class:`~mangadap.proc.spatiallybinnedspectra.SpatiallyBinnedSpectra`
      objects.

        * If rejecting, the size of the boxcar (pixels) is set by
          ``reject_boxcar``.
        * All ``filter_*`` config options are only used with
          ``fit_iter=fit_reject_filter``. **Do not use these
          options!**
        * ``moments``, ``degree``, ``mdegree``, and ``bias`` are passed
          directly to pPXF.

    * Given the template and object spectral range, determine the
      maximum viable fitting range for pPXF using
      :func:`~mangadap.proc.ppxffit.PPXFFit.fitting_mask`.
    * Run through the specified iteration procedure, as selected by
      ``fit_iter`` in config; available options are set by
      :func:`~mangadap.proc.ppxffit.PPXFFit.iteration_modes`.
    * Parse the pPXF results into the data table saved to the reference
      file.

        * Spectra without a fit are flagged as either ``NOFIT`` or
          ``FITFAILED``.
        * Check if returned kinematics are near the imposed boundaries:
          :math:`v \pm 2000` km/s from the input redshift and
          :math:`{\rm d}v/100 < \sigma < 1000` km/s, where :math:`{\rm
          d}v` is the size of the pixel (:math:`\sim 70` km/s).  Leads
          to :ref:`metadatamodel-nearbound` in the ``MAPS`` file.
        * Flag pixels rejected by the sigma-clipping iteration.

    * Calculate the dispersion corrections:

        * First construct three spectra: (1) the optimized template; (2)
          the optimized template redshifted to the best-fitting velocity
          and with a velocity dispersion of 100 km/s; (3) the same as
          spectrum 2 but also convolved to the nominal object spectrum
          resolution.
        * Use pPXF to fit spectra 2 and 3 with spectrum 1.
        * The quadrature difference of the fitted dispersion returned
          for the fit to spectrum 3 and spectrum 2 is provided as the
          correction (``STELLAR_SIGMACORR`` in the ``MAPS`` file)

    * Convert the pPXF velocities and velocity errors to :math:`cz`
      velocities in km/s using
      :func:`~mangadap.proc.ppxffit.PPXFFit.convert_velocity`.

 * Construct stellar-continuum ``BINID`` map. Bin IDs are the same as
   for the binned spectra except that any bin that does not meet the
   S/N limit are given a stellar-continuum bin ID of -1.

