
.. _cube-construction:

DAP LOGCUBE File Construction
=============================

*Analysis class*: :class:`~mangadap.dapfits.construct_cube_file`

*File Root*: ``$MANGA_SPECTRO_ANALYSIS/$MANGADRP_VER/$MANGADAP_VER/[PLATE]/[IFUDESIGN]``

*File Template*: ``manga-[PLATE]-[IFUDESIGN]-LOGCUBE-[DAPTYPE].fits.gz``

*Important class dependencies*:

 * Main full-spectrum modeling classes:

    * :class:`~mangadap.proc.spatiallybinnedspectra.SpatiallyBinnedSpectra`
    * :class:`~mangadap.proc.stellarcontinuummodel.StellarContinuumModel`
    * :class:`~mangadap.proc.emissionlinemodel.EmissionLineModel`

 * :class:`~mangadap.util.fitsutil.DAPFitsUtil`: Contains many
   convenience methods that, e.g., reconstruct a set of maps or
   datacubes based on input data sorted by bin ID and a map with the
   location of each bin ID.

*Algorithm*:

 * Check the input types
 * Construct the 3D cubes for the binned spectra using
   :func:`mangadap.proc.spatiallybinnedspectra.SpatiallyBinnedSpectra.construct_3d_hdu`,
   for the stellar continuum models using
   :func:`mangadap.proc.stellarcontinuummodel.StellarContinuumModel.construct_3d_hdu`,
   and for the emission-line models using
   :func:`mangadap.proc.emissionlinemodel.EmissionLineModel.construct_3d_hdu`
   and the provided objects.
 * Initialize the primary header by copying over the DRP header and
   adding the DAP version information.
 * Construct the ``FLUX``, ``IVAR``, ``MASK``, ``LSF``, ``WAVE``, and
   ``REDCORR`` extensions using the 3D binned spectra cube and
   :func:`~mangadap.dapfits.construct_cube_file.binned_data_cube`.
   The reddening vector is used to *remove* the reddening correction,
   such that the returned spectra match the observed flux in the DRP
   files.
 * Construct the ``MODEL``, ``MODEL_MASK``, ``EMLINE``, ``STELLAR``,
   and ``STELLAR_MASK`` extensions using the model data and
   :func:`~mangadap.dapfits.construct_cube_file.model_cubes`.
 * Combine the ``BINID`` extensions from the
   :class:`~mangadap.proc.spatiallybinnedspectra.SpatiallyBinnedSpectra`,
   :class:`~mangadap.proc.stellarcontinuummodel.StellarContinuumModel`,
   and :class:`~mangadap.proc.emissionlinemodel.EmissionLineModel`
   objects and :func:`~mangadap.dapfits.combine_binid_extensions`.
   This does *not* include the emission-line moment and
   spectral-index bin IDs; those channels in the ``BINID`` extension
   are blank unlike the ``MAPS`` file.
 * Finalize the DAP primary header, which primarily constructs the
   :ref:`metadatamodel-dapqual` bit using
   :func:`~mangadap.dapfits.finalize_dap_primary_header`.

