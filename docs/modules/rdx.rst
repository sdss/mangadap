
.. _drp-redux-assessments:

Basic Reduction Assessments
===========================

*Analysis class*: :class:`~mangadap.proc.reductionassessments.ReductionAssessment`

*Reference root*: see :class:`~mangadap.config.analysisplan.AnalysisPlan.common_path`;
``$MANGA_SPECTRO_ANALYSIS/$MANGADRP_VER/$MANGADAP_VER/common/[PLATE]/[IFUDESIGN]`` for MaNGA

*Reference file*: see :class:`~mangadap.proc.reductionassessments.ReductionAssessment.default_paths`;
``manga-[PLATE]-[IFUDESIGN]-[RDXQA].fits.gz`` for MaNGA

*Optional Parameters*: see :ref:`plan`.  The table below lists the parameters
defined by :class:`~mangadap.proc.reductionassessments.ReductionAssessmentDef`

.. include:: ../tables/reductionassessmentdef.rst

*Important class dependencies*:

 * :class:`~mangadap.datacube.datacube.DataCube`: Base class that
   provides the datacube to be assessed.

*Algorithm*:

 * Ignore any pixels that are either masked by the boolean mask or
   flagged with the flags returned by
   :func:`~mangadap.datacube.datacube.DataCube.do_not_use_flags`.
 * Compute sky coordinates of each spaxel using
   :func:`~mangadap.datacube.datacube.DataCube.mean_sky_coordinates`.
   (``SPX_SKYCOO`` in ``MAPS`` file)
 * Use input ellipticity and position angle parameters to compute
   semi-major axis radii using
   :class:`~mangadap.util.geometry.SemiMajorAxisCoo` (``SPX_ELLCOO``
   in ``MAPS`` file)
 * Determine the fraction of unmasked wavelength channels for each
   spaxel
 * If config specifies covariance calculation, compute the covariance
   using the ``LOGRSS`` file at:

   * the center of the wavelength range if wavelength limits specified
     or,
   * the broad-band weighted center of the response function

 * Compute the (band-weighted) mean signal in each spaxel
   (``SPX_MFLUX`` in ``MAPS`` file), (band-weighted) mean variance in
   each spaxel (``SPX_MFLUX_IVAR`` in ``MAPS`` file), and the
   (band-weighted) mean S/N in each spaxel (``SPX_SNR`` in ``MAPS``
   file), using
   :func:`~mangadap.datacube.datacube.DataCube.flux_stats`.

