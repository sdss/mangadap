
.. include:: include/links.rst

*********************
What's New in the DAP
*********************

MPL-11 (3.1.0)
==============

High-level changes
------------------

 * Most analysis approaches (``DAPTYPE``) now use a set of MaStar Simple
   Stellar Population models as the continuum for the emission-line fitting
   module. The models have been published by `Maraston et al. (2020, MNRAS,
   496, 2962)`_; see also https://www.icg.port.ac.uk/MaStar/. The models used
   are a subset of the full model library, which spans the full range of
   parameter space of the full library but enables a faster execution time.
   The weighting of the different models is **not** constrained to be smooth
   or to match any parametric star-formation/chemical-enrichment history. The
   models included in the fit have the following grid of physical
   parameters::

        Age [Gyr]:    0.00316228, 0.01, 0.03162278, 0.100, 0.300, 1.000, 3.000, 9.000, 14.000
        log(Z/Z_sun): -1.35, -1., -0.7, -0.33, 0, 0.35
        IMF Slope:    2.35 (Salpeter)

 * The MPL-11 results include the usual 3 binning schemes used with the
   MaStar SSP templates, as well as a fourth ``DAPTYPE`` that uses the same
   hierarchically clustered MaStar templates from the MPL-10 release. This
   fourth analysis approach provides a useful baseline for assessing the
   effects of the change in the continuum on the emission-line results; i.e.,
   all differences in the analysis results between the
   ``HYB10-MILESHC-MASTARSSP`` and ``HYB10-MILESHC-MASTARHC2`` results should
   be driven by the difference in the continuum templates used in the
   emission-line fitting module.

User-level changes/bug fixes
----------------------------

 * A bug in the reconstruction of the model spectra, **after** the model had
   been fit, meant that some regions that were masked were incorrectly
   reconstructed. This did not affect the parameters in the MAPS files or the
   emission-line-only or stellar-continuum-only spectra. This has now been
   fixed.
 * Significant revision to the format of the input files used to provide the
   list of emission-lines to fit. This change primarily enables the DAP to
   automatically construct the matrices used by `ppxf`_ to constrain the
   kinematics of different kinematic components. This functionality was used
   to test the results of tying the dispersion of emission lines to, say, all
   be within a factor of 2 of the H-alpha velocity dispersion. In the end,
   these constraints **were not** implemented in this MPL; the tying strategy
   for MPL-11 is identical to MPL-10.
 * The rest wavelengths of a number of emission lines have been updated to
   better account for the affects of fine-structure blending. For example,
   the rest wavelength (in vacuum) of the H-alpha line was 6564.632 angstrom
   in MPL-10, but it is 6564.608 in MPL-11. This is representative of the
   wavelength changes; i.e., the differences are typically a few hundredths
   of an angstrom. This change should have a very limited effect on the
   emission-line measurements.

Under-the-hood changes/bug fixes
--------------------------------

 * Added the :class:`~mangadap.util.datatable.DataTable` class to handle the
   binary tables in the reference files, mostly to enable a uniform data model
   and enable automated documentation.
 * Enabled fitting a restricted wavelength range in the emission-line
   module.
 * Add parameter to change the required fraction of spectral coverage to
   perform a fit.
 * Detailed handling of masks for binning, MW extinction correction, and
   emission-line continuum construction.
 * Enable the output grid for :class:`~mangadap.util.sampling.Resample`
   to be defined directly.
 * Include covariance calculation in
   :class:`~mangadap.util.sampling.Resample`.
 * Add a class that reads trace sets produced by IDLUTILS,
   :class:`~mangadap.util.trace.TraceSet`, but not well tested yet.
 * Fix untested velocity registration and stacking code, and added
   registration test.
 * Upgraded package requirements
 * Line profile in
   :class:`~mangadap.proc.emissionlinetempltes.EmissionLineTemplates` now
   hardwired to be
   :class:`~mangadap.util.lineprofiles.FFTGaussianLSF`.

MPL-10 (3.0.1)
==============

High-level changes
------------------

 * Incorporates the new LSF measurements from the DRP.
 * The emission-line module now uses a new set of hierarchically
   cluster MaStar spectra to model the stellar continuum; the library
   is called ``MASTARHC2`` (for now) to distinguish it from the
   templates used in MPL-9.
 * Significant additions to the spectral index calculations. See
   :class:`~mangadap.proc.spectralindices.AbsorptionLineIndices`,
   :class:`~mangadap.proc.spectralindices.BandheadIndices`, and
   changes to the :ref:`datamodel-maps`. Specifically, with respect
   to the change in the data model, the following channels were
   **removed**: ``SPECINDEX_BCEN``, ``SPECINDEX_BCNT``,
   ``SPECINDEX_RCEN``, and ``SPECINDEX_RCNT``. The following channels
   were **added**:

   * ``SPECINDEX_BF``: A calculation of the index with a modified
     definition from Burstein et al. (1984) and Faber et al. (1985),
     instead of the more typical definitions of Worthey et al. (1994)
     and Trager et al. (1998). See
     :class:`~mangadap.proc.spectralindices.AbsorptionLineIndices`.
     The bandhead or color indices in this extension are identical to
     the values in the ``SPECINDEX`` extension.
   * ``SPECINDEX_BF_IVAR``, ``SPECINDEX_BF_MASK``: Inverse variance and
     mask for the "BF" spectral indices.
   * ``SPECINDEX_BF_CORR``: Stellar velocity dispersion corrections for
     the "BF" indices.
   * ``SPECINDEX_BF_MODEL``: BF indices measured using the best-fitting
     model spectrum.
   * ``SPECINDEX_WGT``: Weights to use when calculating an index by
     aggregating index measurements from many spaxels/bins. See
     :class:`~mangadap.proc.spectralindices.AbsorptionLineIndices`
     and :class:`~mangadap.proc.spectralindices.BandheadIndices`.
   * ``SPECINDEX_WGT_IVAR``, ``SPECINDEX_WGT_MASK``: Inverse variance
     and mask for index weights.
   * ``SPECINDEX_WGT_CORR``: Stellar velocity dispersion corrections
     for the index weights.
   * ``SPECINDEX_WGT_MODEL``: Index weight based on the best-fitting
     model spectrum.

 * The estimated dispersion in the line-spread function of the binned
   spectra is now included in the :ref:`datamodel-cube` in the new
   ``LSF`` extension. In MPL-10, the ``LSFPRE`` extension from the
   DRP LOGCUBE is used for the spectral resolution of each spaxel.
 * DAPall file is now divided into one extension for each analysis
   method; see the desciption of the :ref:`metadatamodel-dapall`.

User-level changes/bug fixes
----------------------------

 * Major revision of the front-end to ease incorporation of non-MaNGA
   datacubes.
 * Major restructuring of directories to ease `pypi`_ distribution and
   provide a more formal structure of the executables in ``bin`` and
   their associated scripts.
 * Added ``bin/dap_status`` to check the status of a batch DAP run.
 * For instantiating a :class:`~mangadap.datacube.manga.MaNGADataCube`,
   changed from using a yanny par file to a configuration (ini) file.
   Included code that can write these files using data from the
   DRPall or DRPComplete files.
 * Spectral resolution is no longer chosen by the spectral binning
   module; instead the spectral resolution is selected when reading
   the datacube (makes much more sense!). Led to some clean-up of the
   binning config files. To select different spectral resolution
   extensions, use command-line arguments in ``rundap`` or
   ``write_dap_config``.
 * Include a default analysis plan, so that executions of ``manga_dap``
   don't require an analysis plan yanny parameter file.

Under-the-hood algorithmic changes
----------------------------------

 * The effective wavelength used for the two sidebands when calculating
   absorption-line indices now uses the center of the band, instead
   of the flux-weighted center of the band.

Under-the-hood changes/bug fixes
--------------------------------

 * Remove "default" prefix from methods in
   :mod:`mangadap.config.defaults`.
 * Remove obsolete ``dapsrc`` keyword arguments
 * Some general clean-up of commented code and docs.
 * Import clean-up, including removal of any ``from __future__``
   imports.
 * ``__credits__`` now include all authors of the DAP papers
 * Added new :class:`~mangadap.datacube.datacube.DataCube` and
   :class:`~mangadap.spectra.rowstackedspectra.RowStackedSpectra`
   classes, beginning the generalization of the interface with the
   input data for use with non-MaNGA data.
 * Integrated use of :class:`~mangadap.datacube.manga.MaNGADataCube`
   instead of :class:`~mangadap.util.drpfits.DRPFits` when executing
   the pipeline. Old :class:`~mangadap.util.drpfits.DRPFits` is now
   deprecated; :class:`~mangadap.util.drpfits.DRPFits` now repurposed
   to provide functionality common to both
   :class:`~mangadap.datacube.manga.MaNGADataCube` and
   :class:`~mangadap.spectra.manga.MaNGARSS`.
 * Included a script that will download data into a new
   ``mangadap/data/remote`` directory for testing.  The directory is not
   included in the repo and has been added to ``.gitignore`` to prevent
   it from being accidentally added.
 * Added a number of tests that use the new remote data,
   including a nominal run of ``manga_dap`` using the ``ALL`` binning
   scheme. These will be skipped if the remote data is not available.
 * Significant improvements and testing of
   :mod:`mangadap.util.covariance`. Ensuring that the correlation
   matrices provided by the DRP can be read by
   :class:`~mangadap.datacube.manga.MaNGADataCube` and are
   effectively identical to the calculation performed by
   :class:`~mangadap.spectra.manga.MaNGARSS`.
 * Moved all core script code from ``bin`` to ``mangadap/scripts``.
   Code in ``bin`` now make simple calls to these scripts.  Moved
   ``rundap.py`` from :mod:`mangadap.survey` to :mod:`mangadap.scripts`.
 * Usage of :class:`ObsInputPar` is now obsolete and deprecated.
 * Docstring updates for modules up through
   :class:`~mangadap.proc.stellarcontinuummodel.StellarContinuumModel`,
   but still lots to do.
 * Fixed a :math:`(1+z)` bug in the calculation of the spectral index
   errors.
 * Package version requirements were updated in general.
 * Now uses ``ppxf`` version 7.2.0, but still uses the default
   ``linear_method``. Use of ``linear_method='lsqlin'`` requires more
   testing.


MPL-9 (2.4.1)
=============

High-level changes
------------------

 * Extensions ``SPX_ELLCOO`` and ``BIN_LWELLCOO`` have a new
   channel.  The 4 channels are now radius in arcsec, radius in
   :math:`R_{\rm eff}`, radius in :math:`h^{-1}{\rm kpc}`, and azimuth; the
   third channel is the new one.  So the channels are::

        C1      = 'Elliptical radius'  / Data in channel 1
        U1      = 'arcsec  '           / Units of data in channel 1
        C2      = 'R/Re    '           / Data in channel 2
        U2      = '' / Units of data in channel 2
        C3      = 'R h/kpc '           / Data in channel 3
        U3      = 'kpc/h   '           / Units of data in channel 3
        C4      = 'Elliptical azimuth' / Data in channel 4
        U4      = 'degrees '           / Units of data in channel 4

 * The following extensions were added to the `MAPS` files:

   * ``EMLINE_SEW_CNT``: Continuum used for summed-flux equivalent-width
     measurement
   * ``EMLINE_GEW_CNT``: Continuum used for Gaussian-fit equivalent-width
     measurement
   * ``SPECINDEX_BCEN``: Flux-weighted center of the blue sideband in the
     spectral-index measurement
   * ``SPECINDEX_BCNT``: Continuum level in the blue sideband used in the
     spectral-index measurement
   * ``SPECINDEX_RCEN``: Flux-weighted center of the red sideband in the
     spectral-index measurement
   * ``SPECINDEX_RCNT``: Continuum level in the red sideband used in the
     spectral-index measurement
   * ``SPECINDEX_MODEL``: Index measurement made using the best-fitting
     stellar continuum model.

 * All emission-line extensions now contain 35 lines, and channel
   ordering of the line list has changed compared to MPL-8.  Current
   ordering is::

        C01     = 'OII-3727'           / Data in channel 1
        C02     = 'OII-3729'           / Data in channel 2
        C03     = 'H12-3751'           / Data in channel 3
        C04     = 'H11-3771'           / Data in channel 4
        C05     = 'Hthe-3798'          / Data in channel 5
        C06     = 'Heta-3836'          / Data in channel 6
        C07     = 'NeIII-3869'         / Data in channel 7
        C08     = 'HeI-3889'           / Data in channel 8
        C09     = 'Hzet-3890'          / Data in channel 9
        C10     = 'NeIII-3968'         / Data in channel 10
        C11     = 'Heps-3971'          / Data in channel 11
        C12     = 'Hdel-4102'          / Data in channel 12
        C13     = 'Hgam-4341'          / Data in channel 13
        C14     = 'HeII-4687'          / Data in channel 14
        C15     = 'Hb-4862 '           / Data in channel 15
        C16     = 'OIII-4960'          / Data in channel 16
        C17     = 'OIII-5008'          / Data in channel 17
        C18     = 'NI-5199 '           / Data in channel 18
        C19     = 'NI-5201 '           / Data in channel 19
        C20     = 'HeI-5877'           / Data in channel 20
        C21     = 'OI-6302 '           / Data in channel 21
        C22     = 'OI-6365 '           / Data in channel 22
        C23     = 'NII-6549'           / Data in channel 23
        C24     = 'Ha-6564 '           / Data in channel 24
        C25     = 'NII-6585'           / Data in channel 25
        C26     = 'SII-6718'           / Data in channel 26
        C27     = 'SII-6732'           / Data in channel 27
        C28     = 'HeI-7067'           / Data in channel 28
        C29     = 'ArIII-7137'         / Data in channel 29
        C30     = 'ArIII-7753'         / Data in channel 30
        C31     = 'Peta-9017'          / Data in channel 31
        C32     = 'SIII-9071'          / Data in channel 32
        C33     = 'Pzet-9231'          / Data in channel 33
        C34     = 'SIII-9533'          / Data in channel 34
        C35     = 'Peps-9548'          / Data in channel 35

 * Fixed flux ratios of doublets were improved using the following
   calculation: :math:`\lambda_2/\lambda_1\ \cdot\
   (M1_1+E2_1)/(M1_2+E2_2)`, where :math:`M1` and :math:`E2` are the
   magnetic dipole and electric quadrapole Einstein :math:`A`
   coefficients and :math:`\lambda` is the line wavelength.  The
   Einstein :math:`A` coefficients are all taken from NIST.  In all but
   the [O II] 7320,7330 complex (not included in MPL-9 fit), the
   magnetic dipole term dominates such that the electric quadrapole
   terms are effectively irrelevant.  The fixed flux ratios are:

   * NeIII-3968/NeIII-3968 = 0.3
   * OIII-4960/OIII-5008 = 0.35
   * OI-6365/OI-6302 = 0.32
   * NII-6549/NII-6585 = 0.34
   * SIII-9071/SIII-9533 = 0.41

 * Lines with tied velocity dispersions (some because their fluxes are
   also tied) are:

   * OII-3727 + OII-3729
   * HeI-3889 + Hzet-3890
   * NeIII-3968 + NeIII-3968
   * OIII-4960 + OIII-5008
   * NI-5199 + NI-5201
   * OI-6302 + OI-6365
   * NII-6549 + NII-6585
   * SIII-9071 + SIII-9533

User-level changes/bug fixes
----------------------------

 * Allow emission-line moment, modeling, and spectral-index keys to be
   ``None``, meaning that the DAP only fits the stellar kinematics but
   still produces the MAPS and LOGCUBE files.
 * Added beta version of hierarchically clustered MaStar templates
 * Added a script that computes the weights used during the
   ``ppxffit_qa`` plot construction.
 
Minor changes/bug fixes
-----------------------

 * Hotfix to accommodate change in padding computation in
   ``ppxf>=6.7.15``
 * Fixed units bug in :func:`mangadap.proc.util.flux_to_fnu`.
 * Fixed bug in templates when using ``iteration_mode =
   'global_template'``
 * Change from ``time.clock()`` to ``time.perf_counter()``
 * Bug fix in record array dimensionality when writing to binary table
 * Minor plotting changes for Overview paper plots
 * Fixed masking of emission-line chi-square measurements
 * Remove ``MASKNAME`` header keyword inherited from DRP files

MPL-8 (2.3.0)
=============

High-level changes
------------------

 * Change to ``DAPTYPE`` construction.  ``DAPTYPE`` is now ``binning`` -
   ``stellar templates`` - ``emission-line templates``.  This was done
   because in future releases we plan to switch the templates used for
   the stellar kinematics (likely to remain ``MILESHC``) to a different
   template set for the emission-line modeling with a longer spectral
   range to fit the full MaNGA spectral range.

 * Three additional emission-lines are fit: He I at 3889 angstroms, and
   the [N I] doublet at 5200 angstroms.  The He I line has its
   dispersion tied to :math:`{\rm H}\zeta` at 3890, and the dispersions of the
   [N I] doublet are tied.

 * The ``MAPS`` file extensions were modified:

   * ``STELLAR_SIGMACORR`` now has two channels.  The first provides the
     correction constructed using the same methodology as in
     MPL-7/DR15; the second provides a correction that we are currently
     testing for robustness as a replacement correction.
   * The data in the ``STELLAR_CONT_FRESID`` and ``STELLAR_CONT_RCHI2``
     extensions has been consolidated into a new single extension,
     ``STELLAR_FOM`` (FOM=figure-of-merit).  Channel 3 of ``STELLAR_FOM`` is
     the same as ``STELLAR_CONT_RCHI2``, and channels 4 and 5 are the same
     as the two channels in ``STELLAR_CONT_FRESID``.  See the data model
     for the full description of the additional channels in this
     extension.
   * 4 additional channels are provided related to the emission lines:

     * ``EMLINE_GA``: The amplitude of the fitted Gaussians
     * ``EMLINE_GANR``: The amplitude over noise of the fitted Gaussians
     * ``EMLINE_FOM``: Full-spectrum figures-of-merit for the
       emission-line module; this now has exactly the same format as the
       STELLAR_FOM extension
     * ``EMLINE_LFOM``: The reduced chi-square in 15 pixel windows around
       each fitted emission line.

 * To make the files easier to use, the ``LOGCUBE`` extensions were
   modified:

   * The ``EMLINE_BASE`` and ``EMLINE_MASK`` extensions have been removed.
   * The following extensions have been added:

     * ``MODEL_MASK`` is the mask to use with the ``MODEL`` extension (the
       result of the emission-line+continuum fit)
     * ``STELLAR`` is best-fitting stellar continuum from the stellar
       continuum fit
     * ``STELLAR_MASK`` is the mask for the stellar-continuum fit. 

 * Numerous QA plots have been added; see the data model description.

Under-the-hood algorithmic changes
----------------------------------

 * Allow the emission-line fitter to use the bin ID numbers directly
   instead of matching the spaxels to bins by coordinate proximity
 * Construction of the parameter tying object in the emission-line
   fitter is now done just before each spectrum is fit by ppxf (not
   globally in Sasuke) to better handle when components are omitted.
 * When deconstructing bins into spaxels for the emission-line modeling
   (hybrid binning scheme), the second fit iteration only fits spaxels
   that are components of binned spectra and does not refit spectra that
   constitute an entire bin themselves. I.e. this removes some largely
   redundant fitting. 
 * :class:`~mangadap.proc.ppxffit.PPXFFit` and
   :class:`~mangadap.proc.sasuke.Sasuke` include calculations of the
   chi-square growth; and changed the names of the growth columns in
   the reference files.
 * Changed definitions of :math:`A` to be the model amplitude;
   :math:`A/N` is the model amplitude divided by the median noise in the
   two sidebands defined for the emission-line EW calculation.
 * Major changes to survey-level execution of the DAP, including that
   the input data is now pulled from the ``DRPall`` file instead of the
   ``plateTargets`` files.
 * Ignore minimum :math:`S/N` limitation in emission-line moments and
   spectral indices for hybrid scheme as a stop-gap to minimize
   differences in moments, models, and indices ``BINID``.

User-level changes/bug fixes
----------------------------

 * Fixed an error in setting up the tying structure passed to pPXF that
   led to the [O II] dispersion maps being fully masked.
 * Fixed an error in the propagation of the error in the passband
   integral calculation, which affects the non-parameteric emission-line
   measurement errors (summed data) and the spectral index errors.  The
   calculation is now formally correct, but these propagated errors are
   still underestimated with respect to more robust calculations via an
   MC.
 * Fixed a bug that omitted the ``FORESTAR`` flag from getting propagated
   in the hybrid binning case.
 * Fixed :math:`\chi^2` calculations reported in ``MAPS`` files (does not
   affect chi-square used during fit optimization) for both the
   stellar-continuum fit and the emission-line modeling.
 * Velocity-dispersion corrections are now applied to the spectral
   indices summary data provided in the ``DAPall`` file.
 * Fixed minor issue in propagating masks from the reference files to
   the maps files; primarily an issue for the hybrid binning scheme.
 
Minor changes/bug fixes
-----------------------

 * Fixed the bug that led to the error in the sigma corrections for
   MPL-7, what were replaced before distributing these data via DR15.
 * Fixed bug that was causing multiple instances of "Warning: converting
   a masked element to nan" during the emission-line moment
   measurements.
 * Significant changes to the pixel resampling code, but has a minor
   effect on the results.
 

