
.. _maps-construction:

DAP MAPS File Construction
==========================

*Analysis class*: :class:`~mangadap.dapfits.construct_maps_file`

*File Root*: ``$MANGA_SPECTRO_ANALYSIS/$MANGADRP_VER/$MANGADAP_VER/[PLATE]/[IFUDESIGN]``

*File Template*: ``manga-[PLATE]-[IFUDESIGN]-MAPS-[DAPTYPE].fits.gz``

*Important class dependencies*:

 * All of the main DAP classes:

    * :class:`~mangadap.proc.reductionassessments.ReductionAssessment`
    * :class:`~mangadap.proc.spatiallybinnedspectra.SpatiallyBinnedSpectra`
    * :class:`~mangadap.proc.stellarcontinuummodel.StellarContinuumModel`
    * :class:`~mangadap.proc.emissionlinemoments.EmissionLineMoments`
    * :class:`~mangadap.proc.emissionlinemodel.EmissionLineModel`
    * :class:`~mangadap.proc.spectralindices.SpectralIndices`

 * :class:`~mangadap.util.fitsutil.DAPFitsUtil`: Contains many
   convenience methods that, e.g., reconstruct a set of maps or
   datacubes based on input data sorted by bin ID and a map with the
   location of each bin ID.

*Algorithm*:

 * Check the input types
 * Initialize the primary header by copying over the DRP header and
   adding the DAP version information.
 * Construct the ``SPX_SKYCOO``, ``SPX_ELLCOO``, ``SPX_MFLUX``,
   ``SPX_MFLUX_IVAR``, and ``SPX_SNR`` maps from the
   :class:`~mangadap.proc.reductionassessments.ReductionAssessment`
   object using
   :func:`~mangadap.dapfits.construct_maps_file.reduction_assessment_maps`.
 * Combine the ``BINID`` extensions from the
   :class:`~mangadap.proc.spatiallybinnedspectra.SpatiallyBinnedSpectra`,
   :class:`~mangadap.proc.stellarcontinuummodel.StellarContinuumModel`,
   :class:`~mangadap.proc.emissionlinemoments.EmissionLineMoments`,
   :class:`~mangadap.proc.emissionlinemodel.EmissionLineModel`, and
   :class:`~mangadap.proc.spectralindices.SpectralIndices` objects
   using :func:`~mangadap.dapfits.combine_binid_extensions`.
 * Construct the ``BIN_LWSKYCOO``, ``BIN_LWELLCOO``, ``BIN_AREA``,
   ``BIN_FAREA``, ``BIN_MFLUX``, ``BIN_MFLUX_IVAR``,
   ``BIN_MFLUX_MASK``, and ``BIN_SNR`` extensions from the
   :class:`~mangadap.proc.spatiallybinnedspectra.SpatiallyBinnedSpectra`
   object using
   :func:`~mangadap.dapfits.construct_maps_file.binned_spectra_maps`.
 * Construct the ``STELLAR_VEL``, ``STELLAR_VEL_IVAR``,
   ``STELLAR_VEL_MASK``, ``STELLAR_SIGMA``, ``STELLAR_SIGMA_IVAR``,
   ``STELLAR_SIGMA_MASK``, ``STELLAR_SIGMACORR``, and ``STELLAR_FOM``
   extensions using the
   :class:`~mangadap.proc.stellarcontinuummodel.StellarContinuumModel`
   object using
   :func:`~mangadap.dapfits.construct_maps_file.stellar_continuum_maps`.
 * Construct the ``EMLINE_SFLUX``, ``EMLINE_SFLUX_IVAR``,
   ``EMLINE_SFLUX_MASK``, ``EMLINE_SEW``, ``EMLINE_SEW_CNT``,
   ``EMLINE_SEW_IVAR``, and ``EMLINE_SEW_MASK`` extensions using the
   :class:`~mangadap.proc.emissionlinemoments.EmissionLineMoments`
   object using
   :func:`~mangadap.dapfits.construct_maps_file.emission_line_moment_maps`.
 * Construct the ``EMLINE_GFLUX``, ``EMLINE_GFLUX_IVAR``,
   ``EMLINE_GFLUX_MASK``, ``EMLINE_GEW``, ``EMLINE_GEW_CNT``,
   ``EMLINE_GEW_IVAR, ``EMLINE_GEW_MASK``, ``EMLINE_GVEL``,
   ``EMLINE_GVEL_IVAR``, ``EMLINE_GVEL_MASK``, ``EMLINE_GSIGMA``,
   ``EMLINE_GSIGMA_IVAR``, ``EMLINE_GSIGMA_MASK``,
   ``EMLINE_INSTSIGMA``, ``EMLINE_TPLSIGMA``, ``EMLINE_GA``,
   ``EMLINE_GANR``, ``EMLINE_FOM``, and ``EMLINE_LFOM`` using the
   :class:`~mangadap.proc.emissionlinemodel.EmissionLineModel` object
   using
   :func:`~mangadap.dapfits.construct_maps_file.emission_line_model_maps`.
 * Construct the ``SPECINDEX``, ``SPECINDEX_IVAR``,
   ``SPECINDEX_MASK``, ``SPECINDEX_CORR``, ``SPECINDEX_MODEL``
   ``SPECINDEX_BF``, ``SPECINDEX_BF_IVAR``, ``SPECINDEX_BF_MASK``,
   ``SPECINDEX_BF_CORR``, ``SPECINDEX_BF_MODEL`` ``SPECINDEX_WGT``,
   ``SPECINDEX_WGT_IVAR``, ``SPECINDEX_WGT_MASK``,
   ``SPECINDEX_WGT_CORR``, ``SPECINDEX_WGT_MODEL`` extensions using
   the :class:`~mangadap.proc.spectralindices.SpectralIndices` object
   using
   :func:`~mangadap.dapfits.construct_maps_file.spectral_index_maps`.
 * Compute the *griz* S/N metrics to include in the header and which
   then get propagated to the DAPall file using
   :func:`~mangadap.dapfits.add_snr_metrics_to_header`.
 * Finalize the DAP primary header, which primarily constructs the
   :ref:`metadatamodel-dapqual` bit using
   :func:`~mangadap.dapfits.finalize_dap_primary_header`.
 * Flag any map data that is not already flagged and does not have a
   positive inverse variance with both the ``MATHERROR`` and
   ``DONOTUSE`` bits.

