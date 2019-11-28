
.. _spatial-binning:

Spatial Binning
===============

*Analysis class:* :class:`mangadap.proc.spatiallybinnedspectra.SpatiallyBinnedSpectra`

*Reference root:* ``$MANGA_SPECTRO_ANALYSIS/$MANGADRP_VER/$MANGADAP_VER/common``

*Reference file:* ``manga-[PLATE]-[IFUDESIGN]-[DRPQA_KEY]-[BIN_KEY].fits.gz``

*Config files:* ``$MANGADAP_DIR/python/mangadap/config/spatial_binning``

*Example config*: ``vor10.ini``

.. code-block:: ini

    [default]
     key                    = VOR10
     galactic_reddening     = ODonnell
     galactic_rv            = 3.1
     method                 = voronoi
     minimum_snr            = 1.0
     operation              = mean
     velocity_register      = False
     stack_covariance_mode  = channels
     stack_covariance_par   = 11
     target_snr             = 10
     spec_res               = spaxel
     prepixel_sres          = True

*Important class dependencies:*

 - :class:`mangadap.util.extinction.GalacticExtinction`: Provides the
   Galactic extinction curves the can be selected and applied to the
   data.
 - :mod:`mangadap.proc.spatialbinning`: Provides the classes that
   perform the binning (e.g.,
   :class:`mangadap.proc.spatialbinning.VoronoiBinning`).

 - :class:`mangadap.proc.spectralstack.SpectralStack`: Provides the core
   functions that perform the spectral stacking

*Algorithm:*

 - Ignore any pixels masked with DONOTUSE or FORESTAR bits by DRP.
 - Calculate the Galactic extinction (REDCORR in the DAP LOGCUBE file)
 - Using
   :class:`mangadap.proc.reductionassessments.ReductionAssessments`
   object, find spaxels with >0.8 fractional spectral coverage and above
   the ``minimum_snr`` in the configuration file.  Only those spaxels
   satisfying both criteria are included in any bin.
 - Determine which spaxels to put in each bin following the ``method``
   specified in the config file:

    - ``none`` (``SPX`` binning type): No binning performed.  Every
      selected spaxels given a unique bin ID.
    - ``global`` (``ALL`` binning type): Bin all valid spaxels into a
      single spectrum.
    - ``radial`` (e.g., ``NRE`` binning type): Use the elliptical
      coordinates from the
      :class:`mangadap.proc.reductionassessments.ReductionAssessments`
      object to assign each spaxel to a unique radial bin.  The binning
      annuli are defined using the ``center``, ``pa``, ``ell``,
      ``radius_scale``, ``radii``, and ``log_step`` config values; see
      ``$MANGADAP_DIR/python/mangadap/config/spatial_binning/nre.ini``
      for the ``NRE`` binning case.  If ``pa``, ``ell``, or
      ``radius_scale`` are -1, they are replaced by ``pa``, ``ell``, and
      ``reff``, respectively, from the :ref:`execution-obs-input`.
    - ``voronoi`` (e.g., ``VOR10`` binning type): Use the Voronoi
      tessellation binning algorithm (written by M. Cappellari; see
      `here <https://pypi.org/project/vorbin/>`__) to continually accrete
      adjacent spaxels to reach a minimum S/N (set by ``target_snr`` in
      config), accounting for covariance if available, using the signal
      and noise measurements from the
      :class:`mangadap.proc.reductionassessments.ReductionAssessments`
      object.

 - Stack all spectra assigned to a single bin:

    - Spectra are combined following the specified ``operation`` in
      config.  Available options are set by
      :func:`mangadap.proc.spectralstack.SpectralStack.operation_options`.
    - Account for covariance according to ``stack_covariance_mode`` and
      ``stack_covariance_par`` in config.  Available options are set by
      :func:`mangadap.proc.spectralstack.SpectralStack.covariance_mode_options`.
    - In the stacked spectra, construct the spectral resolution vectors
      following ``spec_res`` in config.  Available options are set by
      :func:`mangadap.proc.spatiallybinnedspectra.SpatiallyBinnedSpectra.spectral_resolution_options`.
      Use the prepixelized LSF measurements (PREDISP,PRESPECRES in the
      DRP file) if ``prepixel_sres=True``, otherwise use the
      post-pixelized LSF measurements (DISP,SPECRES in the DRP file).
    - Mask any wavelength channels in each spaxel with no unmasked
      pixels from the stack (maskbit set to FLUXINVALID in DAP LOGCUBE
      file).

 - Construct the map with the bin ID of each spaxel (BINID in MAPS file)
 - Calculate the mean signal (BIN_MFLUX in MAPS file), variance (inverse
   of BIN_MFLUX_IVAR in MAPS file) and S/N (BIN_SNR in MAPS file) of the
   stacked spectra.  This is done over the same band/wavelength range as
   done for the individual spaxel data for the
   :class:`mangadap.proc.reductionassessments.ReductionAssessments`
   object.
 - Using the mean signal from the
   :class:`mangadap.proc.reductionassessments.ReductionAssessments`
   object, calculate the luminosity-weighted on-sky (BIN_LWSKYCOO in
   MAPS file) and elliptical (BIN_LWELLCOO in MAPS file) coordinates.
   Also calculate the unweighted coordinates; these are *not* provided
   in the output MAPS file.
 - Calculate the area of each bin (BIN_AREA in MAPS file), and the ratio
   of that area to the expected area (BIN_FAREA in MAPS file) of the
   binning procedure.  The latter is only relevant to the radial
   binning, where the expected area is the area of the bin annulus.
 - Apply the Galactic reddening correction to the binned spectra, where
   the reddening law is defined by the ``galactic_reddening`` and
   ``galactic_rv`` parameters, and the E(B-V) value is taken from the
   DRP header keyword EBVGAL; see
   :class:`mangadap.util.extinction.GalacticExtinction`.  The valid
   reddening laws are:

    - ``ODonnell``: see
      :func:`mangadap.util.extinction.reddening_vector_ccm`.
    - ``CCM``: see
      :func:`mangadap.util.extinction.reddening_vector_ccm`.
    - ``FM``: see :func:`mangadap.util.extinction.reddening_vector_fm`.
    - ``Calzetti``: see
      :func:`mangadap.util.extinction.reddening_vector_calzetti`.

.. note::

    Internally, the DAP performs all spectral fitting on the binned
    spectra (termed as such even if a bin only contains a single spaxel)
    *after* they have been corrected for Galactic extinction.
    Therefore, the output emission-line fluxes have been corrected for
    Galactic extinction.  However, the models and binned spectra in the
    output ``LOGCUBE`` file are reverted to their reddened values for
    direct comparison with the DRP ``LOGCUBE`` file.
