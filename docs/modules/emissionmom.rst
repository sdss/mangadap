
.. _emission-line-moments:

Emission-line Moments
=====================

*Analysis class*: :class:`~mangadap.proc.emissionlinemoments.EmissionLineMoments`

*Reference root*: see :class:`~mangadap.config.analysisplan.AnalysisPlan.method_path`;
``$MANGA_SPECTRO_ANALYSIS/$MANGADRP_VER/$MANGADAP_VER/[METHOD]/[PLATE]/[IFUDESIGN]/ref`` for MaNGA

*Reference file*: see :class:`~mangadap.proc.emissionlinemoments.EmissionLineMoments.default_paths`;
for MaNGA the files are:

    * Before Gaussian emission-line modeling:

        * ``manga-[PLATE]-[IFUDESIGN]-[RDXQA]-[BINNING]-[CONTINUUM]-[ELMOM].fits.gz``

    * After Gaussian emission-line modeling:

        * ``manga-[PLATE]-[IFUDESIGN]-[RDXQA]-[BINNING]-[CONTINUUM]-[ELFIT]-[ELMOM].fits.gz``

*Optional Parameters*: see :ref:`plan`.  The table below lists the parameters
defined by :class:`~mangadap.proc.emissionlinemoments.EmissionLineMomentsDef`

.. include:: ../tables/emissionlinemomentsdef.rst

*Important class dependencies*:

 * :mod:`mangadap.proc.bandpassfilter`: Provides the core functions that
   perform the bandpass integrals.
 * :class:`~mangadap.par.emissionmomentsdb.EmissionMomentsDB`:
   Generalized class that provides the detailed parameters for a set of
   emission-line windows used to perform non-parametric moments.
 * :class:`~mangadap.proc.spectralfitting.EmissionLineFit`: Provides
   functions common to both the moment and Gaussian-fit calculations.
 * :class:`~mangadap.util.pixelmask.SpectralPixelMask`: Used to mask
   spectral regions.

*Algorithm*:

 * Read the artifact database to setup the
   :class:`~mangadap.util.pixelmask.SpectralPixelMask` object based on
   the ``artifact_mask`` config.
 * Set up the
   :class:`~mangadap.par.emissionmomentsdb.EmissionMomentsDB` using
   the ``emission_passbands`` config.
 * Determine the binned spectra above the S/N limit set by the
   ``minimum_snr`` config.
 * Use the
   :class:`~mangadap.proc.stellarcontinuummodel.StellarContinuumModel`
   object to construct the stellar continuum for each binned
   spectrum.
 * Subtract the continuum using
   :func:`~mangadap.proc.spectralfitting.EmissionLineFit.subtract_continuum`.
   WARNING: If a binned spectrum does not have a fitted stellar
   continuum, the moment analysis is performed on the binned spectrum
   without any continuum subtraction.
 * Measure the moments using
   :func:`~mangadap.proc.emissionlinemoments.EmissionLineMoments.measure_moments`.

    * Redshift the emission-line passbands based on the provided
      redshift.
    * Determine the pseudo-continuum in the red and blue bands using
      :func:`~mangadap.proc.bandpassfilter.pseudocontinuum`.
    * Set the slope and intercept of a linear continuum extrapolation
      between the two sidebands for all emission-lines.
    * For each emission line, measure the first 3 moments of the
      pseudo-continuum-subtracted spectra using
      :func:`~mangadap.proc.emissionlinemoments.EmissionLineMoments.single_band_moments`:
      (0) integrated flux; (1) intensity weighted redshift
      (:math:`cz`); and (2) intensity weighted :math:`(cz)^2`.
    * Determine the instrumental dispersion at the 1st moment
      locations of each line using
      :func:`~mangadap.proc.spectralfitting.EmissionLineFit.instrumental_dispersion`.
    * Flag any measurement without a continuum spectrum as
      ``NOCORRECTION``.
    * If any of the passbands (blue, red, main) are incomplete (or
      empty) due to masked pixels or straddle the jump between where
      there is and is not a viable continuum subtracted, or if that
      jump occurs between the blue and red passbands, flag the
      moments as ``FITFAILED`` in the ``MAPS`` file.
    * Mask any "dummy" bands. Dummy bands are used to ensure that the
      emission-line moment channels match the emission-line
      Gaussian-fit channels in the output ``MAPS`` file.

 * Using the 0th moment (integrated flux) and the binned spectra
   (''without'' continuum subtraction), measure the emission-line
   equivalent widths using
   :func:`~mangadap.proc.bandpassfilter.emission_line_equivalent_width`.
 * Construct emission-line-moments ``BINID`` map. Bin IDs are the
   same as for the binned spectra except that any bin that does not
   meet the S/N limit are given a emission-line-moment bin ID of -1.


