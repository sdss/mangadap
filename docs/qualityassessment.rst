
Quality Assessment Plots
========================

The DAP produces a number of quality assessment plots, ranging from
those specific to the results of a given analysis module for a single
observation to plots of global properties for all observations
analyzed for a given MPL. The QA plots are limited in many ways, but
the following is a description of the current QA plots produced.

----

Observation (``PLATEIFU``) specific
-----------------------------------

Spot-check MAPS
~~~~~~~~~~~~~~~

*Script*: ``$MANGADAP_DIR/bin/spotcheck_dap_maps``

*Output root*: ``$MANGA_SPECTRO_ANALYSIS/$MANGADRP_VER/$MANGADAP_VER/[DAPTYPE]/[PLATE]/[IFU]/qa``

*Output file*: ``manga-[PLATE]-[IFUDESIGN]-MAPS-[DAPTYPE]-spotcheck.png``

.. figure:: figures/manga-8138-12704-MAPS-HYB10-MILESHC-MASTARHC-spotcheck.png
   :width: 60 %

   A selection of maps from the ``MAPS`` file for quick assessment of
   the quality of the analysis. The ``PLATEIFU`` and MaNGA ID are
   given at the top of the figure. The panels show from top-to-bottom
   and left-to-right: :math:`g`-band-weighted mean flux
   (``SPX_MFLUX``), :math:`g`-band-weighted spaxel S/N (``SPX_SNR``),
   bin ID number (``BINID``), :math:`g`-band-weighted bin S/N
   (``BIN_SNR``), 68% growth of the fractional residual of the
   stellar-continuum (kinematics) fit (``STELLAR_FOM``; channel
   ``68th perc frac resid``), :math:`\chi^2_{\nu}` of the
   stellar-continuum (kinematics) fit (``STELLAR_FOM``; channel
   ``rchi2``), stellar velocity (``STELLAR_VEL``), (uncorrected)
   stellar velocity dispersion (``STELLAR_SIGMA``), non-parametric
   :math:`{\rm H}\alpha` flux (``EMLINE_SFLUX``; channel
   ``Ha-6564``), Gaussian-fit :math:`{\rm H}\alpha` flux
   (``EMLINE_GFLUX``; channel ``Ha-6564``), ionized-gas velocity
   (``EMLINE_GVEL``; channel ``Ha-6564``), (uncorrected) :math:`{\rm
   H}\alpha` emission-line velocity dispersion (``EMLINE_GSIGMA``;
   channel ``Ha-6564``), non-parametric :math:`{\rm H}\beta` flux
   (``EMLINE_SFLUX``; channel ``Hb-4862``), Gaussian-fit :math:`{\rm
   H}\beta` flux (``EMLINE_GFLUX``; channel ``Hb-4862``), D4000 index
   (``SPECINDEX``; channel ``D4000``), and the Dn4000 index
   (``SPECINDEX``; channel ``Dn4000``).


pPXF Results
~~~~~~~~~~~~

*Script*: ``$MANGADAP_DIR/bin/dap_ppxffit_qa``

*Output root*: ``$MANGA_SPECTRO_ANALYSIS/$MANGADRP_VER/$MANGADAP_VER/[DAPTYPE]/[PLATE]/[IFU]/qa``

*Output file*: ``manga-[PLATE]-[IFUDESIGN]-MAPS-[DAPTYPE]-ppxffit.png``

.. figure:: figures/manga-8138-12704-MAPS-HYB10-MILESHC-MASTARHC-ppxffit.png
   :width: 60 %

   A number of useful (and not so useful) metrics constructed from
   the results of the fit to the stellar continuum performed by pPXF
   for the purpose of measuring the stellar kinematics (see
   :ref:`stellar-kinematics`).


Full-spectrum fit residuals
~~~~~~~~~~~~~~~~~~~~~~~~~~~

*Script*: ``$MANGADAP_DIR/bin/dap_fit_residuals``

*Output root*: ``$MANGA_SPECTRO_ANALYSIS/$MANGADRP_VER/$MANGADAP_VER/[DAPTYPE]/[PLATE]/[IFU]/qa``

*Output files*:

    * ``manga-[PLATE]-[IFUDESIGN]-LOGCUBE-[DAPTYPE]-sc-fitqa-maps.png`` (left)
    * ``manga-[PLATE]-[IFUDESIGN]-LOGCUBE-[DAPTYPE]-el-fitqa-maps.png`` (right)

.. figure:: figures/manga-8138-12704-LOGCUBE-HYB10-MILESHC-MASTARHC-fitqa-maps.png
   :width: 100 %

----

*Output files*:

    * ``manga-[PLATE]-[IFUDESIGN]-LOGCUBE-[DAPTYPE]-sc-fitqa-lambda.png`` (left)
    * ``manga-[PLATE]-[IFUDESIGN]-LOGCUBE-[DAPTYPE]-el-fitqa-lambda.png`` (right)

.. figure:: figures/manga-8138-12704-LOGCUBE-HYB10-MILESHC-MASTARHC-fitqa-lambda.png
   :width: 100 %

----

*Output files*:

    * ``manga-[PLATE]-[IFUDESIGN]-LOGCUBE-[DAPTYPE]-sc-fitqa-growth.png`` (left)
    * ``manga-[PLATE]-[IFUDESIGN]-LOGCUBE-[DAPTYPE]-el-fitqa-growth.png`` (right)

.. figure:: figures/manga-8138-12704-LOGCUBE-HYB10-MILESHC-MASTARHC-fitqa-growth.png
   :width: 100 %

----

Aggregated per plate
--------------------

*Script*: ``$MANGADAP_DIR/bin/dap_plate_fit_qa``

*Output root*: ``$MANGA_SPECTRO_ANALYSIS/$MANGADRP_VER/$MANGADAP_VER/[DAPTYPE]/[PLATE]/qa``

*Output file*: ``[PLATE]-fitqa.png``

.. figure:: figures/8138-fitqa.png
   :width: 60 %

----

Aggregated per DAPTYPE
----------------------

*Script*: ``$MANGADAP_DIR/bin/dapall_qa``

*Output root*: ``$MANGA_SPECTRO_ANALYSIS/$MANGADRP_VER/$MANGADAP_VER/[DAPTYPE]/qa``

*Output file*: ``dapall_radialcoverage.png``

.. figure:: figures/dapall_radialcoverage.png
   :width: 40%

----

*Output file*: ``dapall_redshift_dist.png``

.. figure:: figures/dapall_redshift_dist.png
   :width: 40%

----

*Output file*: ``dapall_mass_vel.png``

.. figure:: figures/dapall_mass_vel.png
   :width: 40%

----

*Output file*: ``dapall_mass_sigma.png``

.. figure:: figures/dapall_mass_sigma.png
   :width: 40%

----

*Output file*: ``dapall_mass_lha.png``

.. figure:: figures/dapall_mass_lha.png
   :width: 40%

----

*Output file*: ``dapall_ew_d4000.png``

.. figure:: figures/dapall_ew_d4000.png
   :width: 40%

----

*Output file*: ``dapall_mgfe_hbeta.png``

.. figure:: figures/dapall_mgfe_hbeta.png
   :width: 40%



