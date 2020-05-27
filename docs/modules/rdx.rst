
.. _drp-redux-assessments:

DRP assessments
===============

*Analysis class*: :class:`~mangadap.proc.reductionassessments.ReductionAssessment`

*Method definition*: :class:`~mangadap.proc.reductionassessments.ReductionAssessmentDef`

*Reference root*: ``$MANGA_SPECTRO_ANALYSIS/$MANGADRP_VER/$MANGADAP_VER/common``

*Reference file*: ``manga-[PLATE]-[IFUDESIGN]-[DRPQA_KEY].fits.gz``

*Config files*: ``$MANGADAP_DIR/mangadap/config/reduction_assessments``

*Example config*: ``snrg.ini``

.. code-block:: ini

    # The response function used is that provided by Jim Gunn in June 2001,
    # found here:
    #
    # http://www.sdss3.org/binaries/filter_curves.fits
    #
    # and parsed into text files.
    [Path]
    dapsrc                 = ${MANGADAP_DIR}

    [default]
    key                    = SNRG
    wave_limits
    response_function_file = ${Path:dapsrc}/data/filter_response/gunn_2001_g_response.db
    in_vacuum              = True
    covariance             = True

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

