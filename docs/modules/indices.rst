
.. _spectral-index-measurements:

Spectral-Index Measurements
===========================

*Analysis class*: :class:`~mangadap.proc.spectralindices.SpectralIndices`

*Method definition*: :class:`~mangadap.proc.spectralindices.SpectralIndicesDef`

*Reference root*: ``$MANGA_SPECTRO_ANALYSIS/$MANGADRP_VER/$MANGADAP_VER/[PLATE]/[IFUDESIGN]/ref``

*Reference file*: ``manga-[PLATE]-[IFUDESIGN]-[DRPQA_KEY]-[BIN_KEY]-[CONTINUUM_KEY]-[ELFIT_KEY]-[SPINDEX_KEY].fits.gz``

*Config files*: ``$MANGADAP_DIR/mangadap/config/spectral_indices``

*Example config*: ``indxen.ini``

.. code-block:: ini

    [default]
    key                      = INDXEN
    minimum_snr              = 0.0
    resolution_fwhm          = -1
    compute_sigma_correction = True
    artifact_mask            = BADSKY
    absorption_indices       = EXTINDX
    bandhead_indices         = BHBASIC

*Important class dependencies*:

 * :class:`~mangadap.par.absorptionindexdb.AbsorptionIndexDB`:
   Generalized class that provides the detailed parameters for a set of
   absorption-line spectral indices.
 * :class:`~mangadap.par.bandheadindexdb.BandheadIndexDB`: Generalized
   class that provides the detailed parameters for a set of bandhead (or
   "color") spectral indices.
 * :mod:`mangadap.proc.bandpassfilter`: Provides the core functions that
   perform the bandpass integrals.

*Algorithm*:

 * Read the artifact database to setup the
   :class:`~mangadap.util.pixelmask.SpectralPixelMask` object based on
   the ``artifact_mask`` config.
 * Setup the :class:`~mangadap.par.absorptionindexdb.AbsorptionIndexDB`
   (using ``absorption_indices`` config) and
   :class:`~mangadap.par.bandheadindexdb.BandheadIndexDB` (using
   ``bandhead_indices`` config) databases with the indices to measure.
 * Determine the binned spectra above the S/N limit set by the
   ``minimum_snr`` config.
 * Mask binned spectra, ignoring pixels masked as ``DONOTUSE``,
   ``IGNORED``, ``FLUXINVALID``, or ``FORESTAR`` in DAP ``LOGCUBE``
   file.
 * Get the best-fitting emission-line models from the
   :class:`~mangadap.proc.emissionlinemodel.EmissionLineModel` object and
   subtract it from them from the data; keep track of where an
   emission-line model is and is not defined.
 * Measure the indices using
   :func:`~mangadap.proc.spectralindices.SpectralIndices.measure_indices`:

    * Compute flux per frequency, needed for some indices; i.e.,
      convert spectra from :math:`F_\lambda` to :math:`F_\nu`.
    * Isolate which indices use each definition (:math:`F_\lambda` vs.
      :math:`F_\nu`)
    * Mask any "dummy" indices.
    * For each spectrum, redshift the band definition, measure the
      absorption-line indices using
      :class:`~mangadap.proc.spectralindices.AbsorptionLineIndices`, and
      the bandhead indices using
      :class:`~mangadap.proc.spectralindices.BandheadIndices`, and save
      the results using
      :func:`~mangadap.proc.spectralindices.SpectralIndices.save_results`.

        * Part of saving the results is to determine which indices
          were successfully measured. Only bands that are completely
          masked (or empty) are flagged as ``NOVALUE`` in the output
          maps. I also keep track of which bands are incomplete (only
          partially masked).

 * Compute the velocity-dispersion corrections:

    * Get the best-fitting continuum model from the
      :class:`~mangadap.proc.stellarcontinuummodel.StellarContinuumModel`,
      both with (``continuum``) and without (``continuum_dcnvlv``)
      the convolution with the best-fitting line-of-sight velocity
      distribution function (LOSVD) using
      :func:`~mangadap.proc.stellarcontinuummodel.StellarContinuumModel.fill_to_match`.
    * Remeasure the indices on these two models (``indx`` and
      ``dcnvlv_indx``, respectively) and the correction based on the
      result using
      :func:`~mangadap.proc.spectralindices.SpectralIndices.calculate_dispersion_corrections`

        * For ``mag`` unit indices, the correction is
          ``dcnvlv_indx-indx``
        * For ``ang`` unit indices, the correction is
          ``dcnvlv_indx/indx``

    * Any index with a bad correction is flagged as ``NOCORRECTION``.

 * Construct spectral-index ``BINID`` map. Bin IDs are the same as
   for the binned spectra except that any bin that does not meet the
   S/N limit are given a spectral-index bin ID of -1.


