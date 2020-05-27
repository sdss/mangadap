
.. include:: include/links.rst

*********************
What's New in the DAP
*********************

MPL-10 (2.5.3dev)
=================

High-level changes
------------------

 * Incorporates the new LSF measurements from the DRP.
 * DAPall file is now divided into one extension for each analysis method.
 * New set of MaSTAR HC spectra

User-level changes/bug fixes
----------------------------

 * Major revision of the front-end to ease incorporation of non-MaNGA
   datacubes.
 * Major restructuring of directories to ease `pypi`_ distribution and
   provide a more formal structure of the executables in ``bin`` and
   their associated scripts.
 * Added ``bin/dap_status`` to check the status of a batch DAP run.
 * For instantiating a MaNGADataCube, changed from using a yanny par
   file to a configuration (ini) file.  Included code that can write
   these files using data from the DRPall or DRPComplete files.
 * Spectral resolution is no longer chosen by the spectral binning
   module; instead the spectral resolution is selected when reading the
   datacube (makes much more sense!).  Led to some clean-up of the
   binning config files.
 * To select different spectral resolution extensions, use command-line
   arguments in ``rundap`` or ``write_dap_config``.
 * Include a default analysis plan, so that executions of ``manga_dap``
   don't require an analysis plan yanny parameter file.

Under-the-hood changes/bug fixes
--------------------------------

 * Remove "default" prefix from methods in
   :mod:`mangadap.config.defaults`.
 * Remove obsolete ``dapsrc`` keyword arguments
 * Some general clean-up of commented code and docs.
 * Import clean-up, including removal of any ``from __future__``
   imports.
 * ``__credits__`` now include all authors of the DAP papers
 * Added new DataCube and RowStackedSpectra classes, beginning the
   generalization of the interface with the input data for use with
   non-MaNGA data.
 * Integrated use of :class:`~mangadap.datacube.manga.MaNGADataCube`
   instead of :class:`~mangadap.util.drpfits.DRPFits` when executing the
   pipeline. Old :class:`DRPFits` is now deprecated; :class:`DRPFits`
   now repurposed to provide functionality common to both
   :class:`MaNGADataCube` and :class:`~mangadap.spectra.manga.MaNGARSS`.
 * Included a script that will download data into a new
   ``mangadap/data/remote`` directory for testing.  The directory is not
   included in the repo and has been added to ``.gitignore`` to prevent
   it from being accidentally added.
 * Included a number of tests that use the new remote data.  These will
   be skipped if the remote data is not available.
 * Significant improvements and testing of
   :mod:`mangadap.util.covariance`.  Ensuring that the correlation
   matrices provided by the DRP can be read by :class:`MaNGADataCube`
   and are effectively identical to the calculation performed by
   :class:`MaNGARSS`.
 * Moved all core script code from ``bin`` to ``mangadap/scripts``.
   Code in ``bin`` now make simple calls to these scripts.  Moved
   ``rundap.py`` from :mod:`mangadap.survey` to :mod:`mangadap.scripts`.
 * Tests now include a nominal run of ``manga_dap`` using the ``ALL``
   binning scheme.
 * Usage of :class:`ObsInputPar` is now obsolete and deprecated.
 * Docstring updates for modules up through
   :class:`StellarContinuumModel`, but still lots to do.


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
   dispersion tied to H:math:`\zeta` at 3890, and the dispersions of the
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
 * :class:`mangadap.proc.ppxffit.PPXFFit` and :class:`mangadap.proc.sasuke.Sasuke`
   include calculations of the chi-square growth; and changed the names
   of the growth columns in the reference files.
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
 

