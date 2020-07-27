
.. include:: ../include/links.rst

.. _dapall-construction:

DAPall File Construction
========================

*Analysis class*: :class:`~mangadap.survey.dapall.DAPall`

*File Root*: ``$MANGA_SPECTRO_ANALYSIS/$MANGADRP_VER/$MANGADAP_VER``

*File Template*: ``dapall-$MANGADRP_VER-$MANGADAP_VER.fits``

*Important class dependencies*:

 * :class:`~mangadap.survey.drpcomplete.DRPComplete`: Provides the
   database used to construct which ``PLATEIFU`` observations were
   analyzed.
 * :class:`~mangadap.par.analysisplan.AnalysisPlan`: Identifies the
   unique analysis approaches used to analyze each datacube.
 * :class:`~mangadap.dapfits.DAPMapsBitMask`: Interprets the masks in
   the ``MAPS`` files.
 * Many of the other core classes are needed but only to define the
   methods used by the analysis approaches selected.
  
*Algorithm*:

 * Instantiate the class used to perform cosmology calculations,
   `astropy.cosmology.FlatLambdaCDM`_, where we set :math:`h=1`,
   :math:`\Omega_M = 0.3`, and :math:`\Omega_\Lambda = 0.7`.
 * Parse the plan file
 * Ensure that all plans to add compute the same emission-line
   moments, emission-line models, and spectral indices.
 * Create list of possibly complete observations and analysis
   products.
 * For each ``MAPS`` file that *should* exist:

   * Find the associated row in the DRPall file, and copy some of those
     data.
   * Calculate the luminosity and angular diameter distance based on
     the NSA redshift
   * Check that the file exists, and if so continue
   * Grab information from the ``MAPS`` file header
   * Calculate the luminosity and angular diameter distance based on
     the input guess redshift (usually the same as the NSA redshift)
   * Calculate the radial coverage metric using
     :func:`~mangadap.survey.dapall.DAPall._radial_coverage_metric`.
   * Pull the S/N metrics from the ``MAPS`` header
   * Get the mean *g*-band surface brightness within 1 :math:`R_e`.
   * Get the binning metrics using
     :func:`~mangadap.survey.dapall.DAPall._binning_metrics`.
   * Get the stellar kinematics metrics using
     :func:`~mangadap.survey.dapall.DAPall._stellar_kinematics_metrics`.
   * Get the :math:`{\rm H}\alpha` kinematics metrics using
     :func:`~mangadap.survey.dapall.DAPall._halpha_kinematics_metrics`.
   * Get the emission-line metrics using
     :func:`~mangadap.survey.dapall.DAPall._emission_line_metrics`.
   * Get the spectral-index metrics using
     :func:`~mangadap.survey.dapall.DAPall._spectral_index_metrics`.
   * Calculate the star-formation rates based on the :math:`{\rm H}\alpha`
     flux within 1 :math:`R_e` and over the full FOV. E.g.,

     .. code-block:: python

        log_Mpc_in_cm = numpy.log10(astropy.constants.pc.to('cm').value) + 6
        log_halpha_luminosity_1re = numpy.log10(4*numpy.pi) \
                    + numpy.log10(db['EMLINE_GFLUX_1RE'][i,self.elfit_channels['Ha-6564']]) \
                    - 17 + 2*numpy.log10(db['LDIST_Z'][i]) + 2*log_Mpc_in_cm
        db['SFR_1RE'][i] = numpy.power(10, log_halpha_luminosity_1re - 41.27)

 * Add the channel names to the header

