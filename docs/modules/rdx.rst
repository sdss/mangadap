
.. _drp-redux-assessments:

DRP assessments
===============

*Analysis class:* :class:`mangadap.proc.reductionassessments.ReductionAssessments`

*Reference root:* ``$MANGA_SPECTRO_ANALYSIS/$MANGADRP_VER/$MANGADAP_VER/common``

*Reference file:* ``manga-[PLATE]-[IFUDESIGN]-[DRPQA_KEY].fits.gz``

*Config files:* ``$MANGADAP_DIR/python/mangadap/config/reduction_assessments``

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

*Important class dependencies:*

 - :class:`mangadap.drpfits.DRPFits`: Advanced class interface to the
   DRP output files that provides most of the assessments

*Algorithm:*

 - Ignore any pixels masked with DONOTUSE or FORESTAR bits by DRP.
 - Compute sky coordinates of each spaxel using the WCS coordinates
   (SPX_SKYCOO in MAPS file)
 - Use input ellipticity and position angle parameters to compute
   semi-major axis radii (SPX_ELLCOO in MAPS file)
 - Determine the fraction of unmasked wavelength channels for each
   spaxel
 - If config specifies covariance calculation, compute the covariance
   using the ``LOGRSS`` file at:

   - the center of the wavelength range if wavelength limits specified
     or,
   - the broad-band weighted center of the response function

 - Compute the (band-weighted) mean signal in each spaxel (SPX_MFLUX in
   MAPS file), (band-weighted) mean variance in each spaxel
   (SPX_MFLUX_IVAR in MAPS file), and the (band-weighted) mean S/N in
   each spaxel (SPX_SNR in MAPS file).

