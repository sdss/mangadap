
.. include:: include/links.rst

.. |ang|   unicode:: U+212B

.. _datamodel:

Data Model
==========

The DAP data model consists of a number of input files and two main
output files, the :ref:`datamodel-maps` and :ref:`datamodel-cube` files,
one for each unique observation (``PLATEIFU``) and analysis approach.

The ``MAPS`` file contains the 2D maps with the DAP-derived properties,
whereas the model ``LOGCUBE`` file contains the best-fitting model
spectra.

.. _datamodel-daptype:

DAP Analysis Approach
---------------------

Each analysis approach used by the DAP is signified by a unique
string called the ``DAPTYPE``, constructed by
:func:`~mangadap.config.defaults.dap_method`. Note that the
construction of the ``DAPTYPE`` has changed since DR15.

.. include:: include/daptype.rst

The table below provides relevant ``DAPTYPE`` keywords:

+---------+-----------------------------------------------------------------+
| Release | ``DAPTYPE``                                                     |
+=========+=================================================================+
| DR-15   | ``VOR10-GAU-MILESHC``, ``HYB10-GAU-MILESHC``                    |
+---------+-----------------------------------------------------------------+
| MPL-11  | ``SPX-MILESHC-MASTARSSP``, ``VOR10-MILESHC-MASTARSSP``,         |
|         | ``HYB10-MILESHC-MASTARSSP``, ``HYB10-MILESHC-MASTARHC2``        |
+---------+-----------------------------------------------------------------+

.. _datamodel-directory-structure:

DAP Directory Structure
-----------------------

The root output directory is codified in the environmental variable
``$MANGA_SPECTRO_ANALYSIS`` (`MANGA_SPECTRO_ANALYSIS`_; internal).

The results of each run of the DAP are tied to the DRP version used
to produce the analyzed datacubes and the version of the DAP used to
do the analysis, such that the top-level directory for the DAP output
is: ``$MANGA_SPECTRO_ANALYSIS/$MANGADRP_VER/$MANGADAP_VER``.

----

The top level within a given DAP version contains the following subdirectories:

 * ``[DAPTYPE]``: User-level directory containing the results from a
   primary analysis method or ``DAPTYPE``.
 * ``log``: Survey-level log files for how the DAP was executed
 * ``common``: Survey-level directory containing files common to
   multiple ``DAPTYPE`` methods.

Most users will only interact with the ``[DAPTYPE]`` directories.  For
the MPL-11, these are:

 * ``SPX-MILESHC-MASTARSSP/``: Analysis of each individual spaxel;
   spaxels must have a valid continuum fit for an emission-line model to
   be fit.
 * ``VOR10-MILESHC-MASTARSSP/``: Analysis of spectra binned to
   :math:`{\rm S/N}\sim 10` using the Voronoi binning algorithm
   (Cappellari & Copin 2003).
 * ``HYB10-MILESHC-MASTARSSP/``: Stellar-continuum analysis of spectra
   binned to :math:`{\rm S/N}\sim 10` for the stellar kinematics (same
   as ``VOR10`` approach); however, the emission-line measurements are
   performed on the individual spaxels.  See a description of the
   :ref:`datamodel-hybrid-binning`.
 * ``HYB10-MILESHC-MASTARHC2/``: Same as the above except the hierarchically
   clustered MaStar stellar spectra are used to fit the stellar continuum in
   the emission-line module.

See the advice for :ref:`gettingstarted-daptype` appropriate for your
science.

----

Within each ``[DAPTYPE]`` directory, you'll find the following:

 * ``[PLATE]``: Top-level directory for analysis of the observations on
   each plate.
 * ``qa``: Quality assessment plots based on data from the DAPall file
   specifically for this ``[DAPTYPE]``.

----

Within each ``[PLATE]`` directory, you'll find the following:

 * ``[IFUDESIGN]``: Subdirectory with the main analysis products for
   each ``PLATEIFU`` combination for the relevant ``[DAPTYPE]``.
   **This is the directory that contains** the DAP ``MAPS`` and model
   ``LOGCUBE`` output files.

 * ``qa``: Quality assessment plots consolidated for all the
   observations on this plate.

----

Within each ``[IFUDESIGN]`` directory, you'll find the following subdirectories:

 * ``qa``: Contains PNG plots for quick QA of the DAP analysis of this
   observation.
 * ``ref``: A reference directory with intermediate files written during
   the analysis.

.. _datamodel-hduclass:

HDUCLASS
--------

As with the `DRP HDUCLASS`_, the DAP output files follow the
``HDUCLASS`` header group convention; see also `here
<ftp://ftp.eso.org/pub/dfs/pipelines/doc/VLT-SPE-ESO-19500-5667_DataFormat.pdf>`_.

All headers specify ``HDUCLASS=SDSS``.  The ``HDUCLASS`` header block
in, e.g., the ``EMLINE_GFLUX`` extension in the ``MAPS`` files looks
like this:

.. code-block:: fortran

    HDUCLASS= 'SDSS    '           / SDSS format class
    HDUCLAS1= 'CUBE    '
    HDUCLAS2= 'DATA    '
    ERRDATA = 'EMLINE_GFLUX_IVAR'  / Associated inv. variance extension
    QUALDATA= 'EMLINE_GFLUX_MASK'  / Associated quality extension

The ``ERRDATA`` and ``QUALDATA`` keywords provide the names of the
extensions with the associated uncertainty and quality flags,
respectively, for the data.  In the headers of the ``EMLINE_GFLUX_IVAR``
and ``EMLINE_GFLUX_MASK`` extensions, the ``SCIDATA`` keyword provides
the extension with the associated science data (``EMLINE_GFLUX`` in this
case).

The ``QUALDATA`` extensions provide bit masks of each property value.
It is **important** that you use these :ref:`metadatamodel-maskbits`
when using the data.

.. _datamodel-maps:

DAP MAPS file
-------------

*File root*: ``$MANGA_SPECTRO_ANALYSIS/$MANGADRP_VER/$MANGADAP_VER/[DAPTYPE]/[PLATE]/[IFUDESIGN]``

*File name*: ``manga-[PLATE]-[IFUDESIGN]-MAPS-[DAPTYPE].fits.gz``

The ``MAPS`` files are the primary output file from the DAP and provide
2D "maps" (i.e., images) of DAP measured properties.  The shape and WCS
of these images identically match that of a single wavelength channel in
the corresponding DRP ``LOGCUBE`` file.  Most properties are provided in
groups of three fits extensions:

  #. ``[property]``: the measurement value,
  #. ``[property]_IVAR``: the measurement uncertainty stored as the
     inverse variance, and
  #. ``[property]_MASK``: a corresponding bit mask for each spaxel.

Extensions can either be a single 2D image (``HDUCLAS1= 'IMAGE'``) or
they can have a series of images that are organized along the third
dimension (``HDUCLAS1= 'CUBE'``).  For the latter, each image is said to
be in a specific "channel".  For example, each Gaussian-fitted
emission-line flux is provided in a single channel in the
``EMLINE_GFLUX`` extension.  The header of extensions with multiple
channels provide the names of the quantities in each channel using
header keyword ``C[n]``, where ``[n]`` is the 1-indexed number of the
channel.

It's best to select the extension and channel based on its *name*, *not*
its extension or channel number; see our
:ref:`gettingstarted-maps-example`.  The ordering of, e.g., the emission
lines in the relevant extensions has changed between different DRs/MPLs
and may change again.

.. note::

    Internally, the DAP performs all spectral fitting on the binned
    spectra (termed as such even if a bin only contains a single spaxel)
    *after* they have been corrected for Galactic extinction.
    Therefore, the output emission-line fluxes have been corrected for
    Galactic extinction.  However, the models and binned spectra in the
    output DAP model ``LOGCUBE`` file are reverted to their reddened
    values for direct comparison with the DRP ``LOGCUBE`` file.

The ``MAPS`` files contain the following extensions:

+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
| HDU |                Name | Channels |                                                Units | Description                                                        |
+=====+=====================+==========+======================================================+====================================================================+
|   0 | PRIMARY             |        0 |                                                      | Empty extension with primary header information.                   |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
| **Coordinate and binning extensions**                                                                                                                            |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|   1 | SPX_SKYCOO          |        2 |                                               arcsec | Sky-right offsets -- +x toward +RA and +y toward +DEC -- of each   |
|     |                     |          |                                                      | spaxel from the galaxy center                                      |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|   2 | SPX_ELLCOO          |        4 |      arcsec, unitless, :math:`h^{-1} {\rm kpc}`, deg | Elliptical polar coordinates of each spaxel from the galaxy        |
|     |                     |          |                                                      | center; :math:`R` in arcsec, :math:`R/R_e`, :math:`R` in           |
|     |                     |          |                                                      | :math:`h^{-1} {\rm kpc}`, and azimuthal angle :math:`\theta`.  In  |
|     |                     |          |                                                      | the limit of tilted thin disk, these are the in-plane disk radius  |
|     |                     |          |                                                      | and azimuth.                                                       |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|   3 | SPX_MFLUX           |        1 |       :math:`10^{-17} {\rm erg/s/cm}^2`/|ang|/spaxel | g-band-weighted mean flux, *not* corrected for Galactic extinction |
|     |                     |          |                                                      | or internal attenuation.                                           |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|   4 | SPX_MFLUX_IVAR      |        1 |                                                      | Inverse variance of g-band-weighted mean flux.                     |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|   5 | SPX_SNR             |        1 |                                                      | Mean g-band weighted signal-to-noise ratio per pixel.              |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|   6 | BINID               |        5 |                                                      | Numerical ID for spatial bins for the binned spectra,              |
|     |                     |          |                                                      | stellar-continuum results, emission-line moment results,           |
|     |                     |          |                                                      | emission-line model results, and spectral-index results;           |
|     |                     |          |                                                      | see :ref:`datamodel-binid-usage`.                                  |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|   7 | BIN_LWSKYCOO        |        2 |                                               arcsec | Light-weighted sky-right offsets -- +x toward +RA and +y toward    |
|     |                     |          |                                                      | +DEC -- of each bin from the galaxy center.                        |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|   8 | BIN_LWELLCOO        |        4 |       arcsec, unitless, :math:`h^{-1} {\rm kpc}`,deg | Light-weighted elliptical polar coordinates of each bin from the   |
|     |                     |          |                                                      | galaxy center; :math:`R` in arcsec, :math:`R/R_e`, :math:`R` in    |
|     |                     |          |                                                      | :math:`h^{-1} {\rm kpc}`, and azimuthal angle :math:`\theta`.  In  |
|     |                     |          |                                                      | the limit of tilted thin disk, these are the in-plane disk radius  |
|     |                     |          |                                                      | and azimuth.                                                       |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|   9 | BIN_AREA            |        1 |                               :math:`{\rm arcsec}^2` | Area of each bin.                                                  |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|  10 | BIN_FAREA           |        1 |                                                      | Fractional area that the bin covers for the expected bin shape     |
|     |                     |          |                                                      | (only relevant for radial binning).                                |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|  11 | BIN_MFLUX           |        1 |       :math:`10^{-17} {\rm erg/s/cm}^2`/|ang|/spaxel | g-band-weighted mean flux for the binned spectra, *not* corrected  |
|     |                     |          |                                                      | for Galactic extinction or internal attenuation.                   |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|  12 | BIN_MFLUX_IVAR      |        1 |                                                      | Inverse variance of g-band-weighted mean flux for the binned       |
|     |                     |          |                                                      | spectra.                                                           |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|  13 | BIN_MFLUX_MASK      |        1 |                                                      | Bit mask for the g-band-weighted mean flux per bin.                |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|  14 | BIN_SNR             |        1 |                                                      | Mean g-band-weighted signal-to-noise ratio per pixel in the binned |
|     |                     |          |                                                      | spectra.                                                           |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
| **Stellar (absorption-line) kinematics**                                                                                                                         |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|  15 | STELLAR_VEL         |        1 |                                                 km/s | Line-of-sight stellar velocity, relative to the input guess        |
|     |                     |          |                                                      | redshift (given as :math:`cz` by the keyword ``SCINPVEL`` in the   |
|     |                     |          |                                                      | header of the ``PRIMARY`` extension, and most often identical to   |
|     |                     |          |                                                      | the NSA redshift).                                                 |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|  16 | STELLAR_VEL_IVAR    |        1 |                                                      | Inverse variance of stellar velocity measurements.                 |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|  17 | STELLAR_VEL_MASK    |        1 |                                                      | Data quality mask for stellar velocity measurements.               |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|  18 | STELLAR_SIGMA       |        1 |                                                 km/s | Raw line-of-sight stellar velocity dispersion; see                 |
|     |                     |          |                                                      | :ref:`corrections` for how to use the ``STELLAR_SIGMACORR`` to     |
|     |                     |          |                                                      | obtain the *astrophysical* stellar velocity dispersion.            |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|  19 | STELLAR_SIGMA_IVAR  |        1 |                                                      | Inverse variance of raw stellar velocity dispersion.               |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|  20 | STELLAR_SIGMA_MASK  |        1 |                                                      | Data quality mask for stellar velocity dispersion.                 |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|  21 | STELLAR_SIGMACORR   |        2 |                                                 km/s | Quadrature correction for STELLAR_SIGMA to obtain the              |
|     |                     |          |                                                      | astrophysical velocity dispersion; see :ref:`corrections` for how  |
|     |                     |          |                                                      | to use this extension with the ``STELLAR_SIGMA`` extension to      |
|     |                     |          |                                                      | obtain the *astrophysical* stellar velocity dispersion.            |
|     |                     |          |                                                      | **Use the first channel.**                                         |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|  22 | STELLAR_FOM         |        9 |                                                      | Figures-of-merit for the stellar-continuum fit in 9 channels: (1)  |
|     |                     |          |                                                      | RMS of residuals (in                                               |
|     |                     |          |                                                      | :math:`10^{-17} {\rm erg/s/cm}^2`/|ang|/spaxel), (2) RMS of        |
|     |                     |          |                                                      | fractional residuals, (3) reduced :math:`\chi^2`, (4-6) 68th and   |
|     |                     |          |                                                      | 99th percentile and maximum value of fractional residuals, and     |
|     |                     |          |                                                      | (7-9) 68th and 99th percentile and maximum value of                |
|     |                     |          |                                                      | error-normalized residual (:math:`\chi`).                          |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
| **Emission-line measurements**                                                                                                                                   |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|  23 | EMLINE_SFLUX        |       35 |       :math:`10^{-17} {\rm erg/s/cm}^2{\rm /spaxel}` | Non-parametric summed flux *after subtracting the*                 |
|     |                     |          |                                                      | *stellar-continuum model*.  The emission-line fluxes account for   |
|     |                     |          |                                                      | Galactic reddening using the E(B-V) value (copied to the DAP       |
|     |                     |          |                                                      | primary headers, see the ``EBVGAL`` header keyword) provided by    |
|     |                     |          |                                                      | the DRP header and assuming an O’Donnell (1994, ApJ, 422, 158)     |
|     |                     |          |                                                      | reddening law; however, no attenuation correction is applied due   |
|     |                     |          |                                                      | to dust internal to the galaxy.                                    |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|  24 | EMLINE_SFLUX_IVAR   |       35 |                                                      | Inverse variance for summed flux measurements.                     |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|  25 | EMLINE_SFLUX_MASK   |       35 |                                                      | Data quality mask for summed flux measurements.                    |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|  26 | EMLINE_SEW          |       35 |                                                |ang| | Non-parametric equivalent widths measurements (based on            |
|     |                     |          |                                                      | the non-parametric fluxes in ``EMLINE_SFLUX``).                    |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|  27 | EMLINE_SEW_CNT      |       35 |       :math:`10^{-17} {\rm erg/s/cm}^2`/|ang|/spaxel | Continuum value used to compute the emission-line equivalent width |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|  28 | EMLINE_SEW_IVAR     |       35 |                                                      | Inverse variance for non-parametric equivalent width measurements. |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|  29 | EMLINE_SEW_MASK     |       35 |                                                      | Data quality mask for non-parametric equivalent width measurements |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|  30 | EMLINE_GFLUX        |       35 |       :math:`10^{-17} {\rm erg/s/cm}^2{\rm /spaxel}` | Gaussian profile integrated flux *from a combined*                 |
|     |                     |          |                                                      | *continuum+emission-line fit*.  The flux ratio of the [OIII],      |
|     |                     |          |                                                      | [OI], and [NII] lines are fixed and cannot be treated as           |
|     |                     |          |                                                      | independent measurements.  The emission-line fluxes account for    |
|     |                     |          |                                                      | Galactic reddening using the E(B-V) (copied to the DAP primary     |
|     |                     |          |                                                      | headers, see the ``EBVGAL`` header keyword) value provided by the  |
|     |                     |          |                                                      | DRP header and assuming an O’Donnell (1994, ApJ, 422, 158)         |
|     |                     |          |                                                      | reddening law; however, no attenuation correction is applied due   |
|     |                     |          |                                                      | to dust internal to the galaxy.                                    |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|  31 | EMLINE_GFLUX_IVAR   |       35 |                                                      | Inverse variance for Gaussian flux measurements                    |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|  32 | EMLINE_GFLUX_MASK   |       35 |                                                      | Data quality mask for Gaussian flux measurements                   |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|  33 | EMLINE_GEW          |       35 |                                                |ang| | Gaussian-fitted equivalent widths measurements (based on the       |
|     |                     |          |                                                      | parametric fluxes in ``EMLINE_GFLUX``).                            |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|  34 | EMLINE_GEW_CNT      |       35 |       :math:`10^{-17} {\rm erg/s/cm}^2`/|ang|/spaxel | Continuum value used to compute the emission-line equivalent width |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|  35 | EMLINE_GEW_IVAR     |       35 |                                                      | Inverse variance of the above.                                     |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|  36 | EMLINE_GEW_MASK     |       35 |                                                      | Data quality mask of the above.                                    |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|  37 | EMLINE_GVEL         |       35 |                                                 km/s | Line-of-sight emission-line velocity, relative to the input guess  |
|     |                     |          |                                                      | redshift (given as :math:`cz` by the keyword ``SCINPVEL`` in the   |
|     |                     |          |                                                      | header of the ``PRIMARY`` extension, and most often identical to   |
|     |                     |          |                                                      | the NSA redshift).  A velocity is provided for each line,          |
|     |                     |          |                                                      | **but the velocities are identical for all lines** because the     |
|     |                     |          |                                                      | parameters are tied during the fitting process.                    |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|  38 | EMLINE_GVEL_IVAR    |       35 |                                                      | Inverse variance for Gaussian-fitted velocity measurements, which  |
|     |                     |          |                                                      | are **the same for all lines and should not be combined as if**    |
|     |                     |          |                                                      | **independent measurements**.                                      |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|  39 | EMLINE_GVEL_MASK    |       35 |                                                      | Data quality mask for Gaussian-fitted velocity measurements.       |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|  40 | EMLINE_GSIGMA       |       35 |                                                 km/s | Gaussian profile velocity dispersion as would be measured from a   |
|     |                     |          |                                                      | direct Gaussian fit; see :ref:`corrections` for how                |
|     |                     |          |                                                      | to use the ``EMLINE_INSTSIGMA`` extension with these data to       |
|     |                     |          |                                                      | obtain the *astrophysical* gas velocity dispersion.  Tied velocity | 
|     |                     |          |                                                      | dispersions ([OII], [OIII], [OI], [NII], [NI] and H-zeta+HeI 3889) |
|     |                     |          |                                                      | cannot be treated as independent measurements.                     |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|  41 | EMLINE_GSIGMA_IVAR  |       35 |                                                      | Inverse variance for Gaussian profile velocity dispersion.         |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|  42 | EMLINE_GSIGMA_MASK  |       35 |                                                      | Data quality mask for Gaussian profile velocity dispersion.        |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|  43 | EMLINE_INSTSIGMA    |       35 |                                                 km/s | The instrumental dispersion at the fitted center of each emission  |
|     |                     |          |                                                      | line.                                                              |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|  44 | EMLINE_TPLSIGMA     |       35 |                                                 km/s | The dispersion of each emission line used in the template spectra; |
|     |                     |          |                                                      | see :ref:`datamodel-eml-tpl-resolution`.                           |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|  45 | EMLINE_GA           |       35 |       :math:`10^{-17} {\rm erg/s/cm}^2`/|ang|/spaxel | The amplitude of the model Gaussian fit to each emission line.     |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|  46 | EMLINE_GANR         |       35 |                                                      | The amplitude of the model Gaussian fit relative to the median     |
|     |                     |          |                                                      | noise in two sidebands near the line; the sidebands are identical  |
|     |                     |          |                                                      | to those used in the equivalent width measurement.                 |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|  47 | EMLINE_FOM          |        9 |                                                      | Figures-of-merit for the continuum+emission-line model fit in 9    |
|     |                     |          |                                                      | channels: (1) RMS of residuals (in                                 |
|     |                     |          |                                                      | :math:`10^{-17} {\rm erg/s/cm}^2`/|ang|/spaxel), (2) RMS of        |
|     |                     |          |                                                      | fractional residuals, (3) reduced :math:`\chi^2`, (4-6) 68th and   |
|     |                     |          |                                                      | 99th percentile and maximum value of fractional residuals, and     |
|     |                     |          |                                                      | (7-9) 68th and 99th percentile and maximum value of                |
|     |                     |          |                                                      | error-normalized residual (:math:`\chi`).                          |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|  48 | EMLINE_LFOM         |       35 |                                                      | The reduced :math:`\chi^2` of the fit to each line calculated in   |
|     |                     |          |                                                      | 15-pixel windows centered on each line.                            |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
| **Spectral index measurements**                                                                                                                                  |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|  49 | SPECINDEX           |       46 |                                           |ang|, mag | Spectral-index measurements.  Indices follow the definition from   |
|     |                     |          |                                                      | Worthey et al. (1994) and Trager et al. (1998), as used in all     |
|     |                     |          |                                                      | previous releases.  See :ref:`spectralindices`.                    |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|  50 | SPECINDEX_IVAR      |       46 |                                                      | Inverse variance for spectral index maps.                          |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|  51 | SPECINDEX_MASK      |       46 |                                                      | Data quality mask for spectral index maps.                         |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|  52 | SPECINDEX_CORR      |       46 |                                                  mag | Corrections to apply to account for the velocity dispersion and    |
|     |                     |          |                                                      | effectively determine the index without Doppler broadening;        |
|     |                     |          |                                                      | see :ref:`corrections`.                                            | 
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|  53 | SPECINDEX_MODEL     |       46 |                                           |ang|, mag | Spectral-index measurements for the best-fitting model spectrum.   |
|     |                     |          |                                                      | Note the extension number is different from MPL-9.                 |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|  54 | SPECINDEX_BF        |       46 |                                           |ang|, mag | Spectral-index measurements calculated using a                     |
|     |                     |          |                                                      | definition similar to Burstein et al. (1984) and                   |
|     |                     |          |                                                      | Faber et al. (1985).  See :ref:`spectralindices`.                  |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|  55 | SPECINDEX_BF_IVAR   |       46 |                                                      | Inverse variance in the BF spectral indices.                       |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|  56 | SPECINDEX_BF_MASK   |       46 |                                                      | Data quality mask for the BF spectral indices.                     |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|  57 | SPECINDEX_BF_CORR   |       46 |                                                  mag | Corrections to apply to account for the                            |
|     |                     |          |                                                      | velocity dispersion and effectively determine the index without    |
|     |                     |          |                                                      | Doppler broadening; see :ref:`corrections`.                        | 
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|  58 | SPECINDEX_BF_MODEL  |       46 |                                           |ang|, mag | Spectral indices with the BF definition                            |
|     |                     |          |                                                      | measured using the best-fitting model spectrum.                    |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|  59 | SPECINDEX_WGT       |       46 |       :math:`10^{-17} {\rm erg/s/cm}^2`/|ang|/spaxel | Weights to use when aggregating spectral-index                     |
|     |                     |          |                                                      | measurements.  See :ref:`spectralindices` and :ref:`aggregation`.  |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|  60 | SPECINDEX_WGT_IVAR  |       46 |                                                      | Inverse variance in the spectral index weights.                    |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|  61 | SPECINDEX_WGT_MASK  |       46 |                                                      | Data quality mask for spectral index weights.                      |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|  62 | SPECINDEX_WGT_CORR  |       46 |                                                      | Corrections to apply to account for the                            |
|     |                     |          |                                                      | effect of the velocity dispersion on the calculation of the index  |
|     |                     |          |                                                      | weight; see :ref:`corrections`.                                    | 
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+
|  63 | SPECINDEX_WGT_MODEL |       46 |       :math:`10^{-17} {\rm erg/s/cm}^2`/|ang|/spaxel | Spectral-index weights determined by the                           |
|     |                     |          |                                                      | best-fitting models.                                               |
+-----+---------------------+----------+------------------------------------------------------+--------------------------------------------------------------------+

.. _datamodel-emission-line-channels:

The emission-line measurements for MPL-11 are (identical to MPL-10):

.. code-block:: fortran

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

.. note::

    For the emission-line moments (``SFLUX``, ``SEW``):

      * Channels 2 ('OII-3729'), 8 ('HeI-3889'), 10 ('NeIII-3968'),
        and 19 ('NI-5201') are empty because the line falls in the
        passband of another line: 'OII-3729' in 'OIId-3728',
        'HeI-3889' in 'Hzet-3890', 'NeIII-3968' in 'Heps-3971', and
        'NI-5201' in 'NI-5199'. To compare these fluxes with the
        Gaussian-fitted values, you should sum the Gaussian-fitted
        fluxes first.
      * OIId is contaminated by H14 and H13
      * Hzet is contaminated by HeI
      * Heps is contaminated by NeIII
      * Red sideband of Hbeta is contaminated by HeI
      * Unknown line at 4990 and may contaminate red sideband of OIII
        4960 and the blue sideband of OIII 5008
      * OIII 5008 contaminated by HeI 5017

----

.. _datamodel-spectral-index-channels:

The spectral-index measurements for MPL-11 are below (identical to
MPL-10). Because the spectral-index measurements can be either
angstroms, magnitudes, or unitless, the header of the spectral-index
extensions also include the units using header keywords ``U[n]``. The
indices and relevant units as included in the relevant extension
header are:

.. code-block:: fortran

    C01     = 'CN1     '           / Data in channel 1
    U01     = 'mag     '           / Units of data in channel 1
    C02     = 'CN2     '           / Data in channel 2
    U02     = 'mag     '           / Units of data in channel 2
    C03     = 'Ca4227  '           / Data in channel 3
    U03     = 'ang     '           / Units of data in channel 3
    C04     = 'G4300   '           / Data in channel 4
    U04     = 'ang     '           / Units of data in channel 4
    C05     = 'Fe4383  '           / Data in channel 5
    U05     = 'ang     '           / Units of data in channel 5
    C06     = 'Ca4455  '           / Data in channel 6
    U06     = 'ang     '           / Units of data in channel 6
    C07     = 'Fe4531  '           / Data in channel 7
    U07     = 'ang     '           / Units of data in channel 7
    C08     = 'C24668  '           / Data in channel 8
    U08     = 'ang     '           / Units of data in channel 8
    C09     = 'Hb      '           / Data in channel 9
    U09     = 'ang     '           / Units of data in channel 9
    C10     = 'Fe5015  '           / Data in channel 10
    U10     = 'ang     '           / Units of data in channel 10
    C11     = 'Mg1     '           / Data in channel 11
    U11     = 'mag     '           / Units of data in channel 11
    C12     = 'Mg2     '           / Data in channel 12
    U12     = 'mag     '           / Units of data in channel 12
    C13     = 'Mgb     '           / Data in channel 13
    U13     = 'ang     '           / Units of data in channel 13
    C14     = 'Fe5270  '           / Data in channel 14
    U14     = 'ang     '           / Units of data in channel 14
    C15     = 'Fe5335  '           / Data in channel 15
    U15     = 'ang     '           / Units of data in channel 15
    C16     = 'Fe5406  '           / Data in channel 16
    U16     = 'ang     '           / Units of data in channel 16
    C17     = 'Fe5709  '           / Data in channel 17
    U17     = 'ang     '           / Units of data in channel 17
    C18     = 'Fe5782  '           / Data in channel 18
    U18     = 'ang     '           / Units of data in channel 18
    C19     = 'NaD     '           / Data in channel 19
    U19     = 'ang     '           / Units of data in channel 19
    C20     = 'TiO1    '           / Data in channel 20
    U20     = 'mag     '           / Units of data in channel 20
    C21     = 'TiO2    '           / Data in channel 21
    U21     = 'mag     '           / Units of data in channel 21
    C22     = 'HDeltaA '           / Data in channel 22
    U22     = 'ang     '           / Units of data in channel 22
    C23     = 'HGammaA '           / Data in channel 23
    U23     = 'ang     '           / Units of data in channel 23
    C24     = 'HDeltaF '           / Data in channel 24
    U24     = 'ang     '           / Units of data in channel 24
    C25     = 'HGammaF '           / Data in channel 25
    U25     = 'ang     '           / Units of data in channel 25
    C26     = 'CaHK    '           / Data in channel 26
    U26     = 'ang     '           / Units of data in channel 26
    C27     = 'CaII1   '           / Data in channel 27
    U27     = 'ang     '           / Units of data in channel 27
    C28     = 'CaII2   '           / Data in channel 28
    U28     = 'ang     '           / Units of data in channel 28
    C29     = 'CaII3   '           / Data in channel 29
    U29     = 'ang     '           / Units of data in channel 29
    C30     = 'Pa17    '           / Data in channel 30
    U30     = 'ang     '           / Units of data in channel 30
    C31     = 'Pa14    '           / Data in channel 31
    U31     = 'ang     '           / Units of data in channel 31
    C32     = 'Pa12    '           / Data in channel 32
    U32     = 'ang     '           / Units of data in channel 32
    C33     = 'MgICvD  '           / Data in channel 33
    U33     = 'ang     '           / Units of data in channel 33
    C34     = 'NaICvD  '           / Data in channel 34
    U34     = 'ang     '           / Units of data in channel 34
    C35     = 'MgIIR   '           / Data in channel 35
    U35     = 'ang     '           / Units of data in channel 35
    C36     = 'FeHCvD  '           / Data in channel 36
    U36     = 'ang     '           / Units of data in channel 36
    C37     = 'NaI     '           / Data in channel 37
    U37     = 'ang     '           / Units of data in channel 37
    C38     = 'bTiO    '           / Data in channel 38
    U38     = 'mag     '           / Units of data in channel 38
    C39     = 'aTiO    '           / Data in channel 39
    U39     = 'mag     '           / Units of data in channel 39
    C40     = 'CaH1    '           / Data in channel 40
    U40     = 'mag     '           / Units of data in channel 40
    C41     = 'CaH2    '           / Data in channel 41
    U41     = 'mag     '           / Units of data in channel 41
    C42     = 'NaISDSS '           / Data in channel 42
    U42     = 'ang     '           / Units of data in channel 42
    C43     = 'TiO2SDSS'           / Data in channel 43
    U43     = 'mag     '           / Units of data in channel 43
    C44     = 'D4000   '           / Data in channel 44
    U44     = '' / Units of data in channel 44
    C45     = 'Dn4000  '           / Data in channel 45
    U45     = '' / Units of data in channel 45
    C46     = 'TiOCvD  '           / Data in channel 46
    U46     = '' / Units of data in channel 46

----

.. _datamodel-cube:

DAP Model LOGCUBE file
----------------------

*File root*: ``$MANGA_SPECTRO_ANALYSIS/$MANGADRP_VER/$MANGADAP_VER/[DAPTYPE]/[PLATE]/[IFUDESIGN]``

*File name*: ``manga-[PLATE]-[IFUDESIGN]-LOGCUBE-[DAPTYPE].fits.gz``

The ``LOGCUBE`` files provide the binned spectra and the best-fitting
model spectrum for each spectrum that was successfully fit. These
files are useful for detailed assessments of the model parameters
because they allow you to return to the spectra and compare the model
against the data. As described by `Westfall et al. (2019, AJ, 158,
231)`_, the DAP fits the spectra in two stages, one to get the
stellar kinematics and the second to determine the emission-line
properties. The emission-line module (used for all binning schemes)
fits both the stellar continuum and the emission lines at the same
time, where the stellar kinematics are fixed by the first fit. The
stellar-continuum models from the first fit are provided in the
``STELLAR`` extension; to get the stellar continuum determined during
the emission-line modeling, you have to subtract the emission-line
model (in the ``EMLINE`` extension) from the full model (in the
``MODEL`` extension). Our :ref:`gettingstarted-cube-example` shows
how to plot the model LOGCUBE data.

.. warning::

    In the ``HYB`` binning case the binned spectra provided in the
    ``LOGCUBE`` files are from the Voronoi binning step.  However, the
    emission-line models are fit to the *individual spaxels*.  So:

        - The stellar-continuum fits from the first iteration, in the
          ``STELLAR`` extension, should be compared to the Voronoi
          binned spectra in the file, but
        - the best-fitting model spectra in the ``MODEL`` extension
          should be compared to the individual spectra from the DRP
          ``LOGCUBE`` file!

.. note::

    Internally, the DAP performs all spectral fitting on the binned
    spectra (termed as such even if a bin only contains a single spaxel)
    *after* they have been corrected for Galactic extinction.
    Therefore, the output emission-line fluxes have been corrected for
    Galactic extinction.  However, the models and binned spectra in the
    output DAP model ``LOGCUBE`` file are reverted to their reddened
    values for direct comparison with the DRP ``LOGCUBE`` file.

The ``LOGCUBE`` files contain the following extensions:

+-----+--------------------+------------------------------------------------------+-----------------------------------------------------------------------+
| HDU |               Name |                                                Units | Description                                                           |
+=====+====================+======================================================+=======================================================================+
|   0 |            PRIMARY |                                                      | Empty extension with primary header information.                      |
+-----+--------------------+------------------------------------------------------+-----------------------------------------------------------------------+
|   1 |               FLUX |       :math:`10^{-17} {\rm erg/s/cm}^2`/|ang|/spaxel | Flux of the ''binned'' spectra                                        |
+-----+--------------------+------------------------------------------------------+-----------------------------------------------------------------------+
|   2 |               IVAR |                                                      | Inverse variance in the binned spectra                                |
+-----+--------------------+------------------------------------------------------+-----------------------------------------------------------------------+
|   3 |               MASK |                                                      | Bitmask for the binned spectra.  Note that this mask only applies to  |
|     |                    |                                                      | the binned spectra.                                                   |
+-----+--------------------+------------------------------------------------------+-----------------------------------------------------------------------+
|   4 |                LSF |                                                |ang| | The dispersion (:math:`\sigma`) of the Gaussian                       |
|     |                    |                                                      | line-spread function of the binned spectra.  For MPL-11, the source   |
|     |                    |                                                      | DRP extension is ``LSFPRE``.                                          |
+-----+--------------------+------------------------------------------------------+-----------------------------------------------------------------------+
|   5 |               WAVE |                                                |ang| | Vacuum-wavelength vector                                              |
+-----+--------------------+------------------------------------------------------+-----------------------------------------------------------------------+
|   6 |            REDCORR |                                                      | Reddening correction applied during the fitting procedures.           |
+-----+--------------------+------------------------------------------------------+-----------------------------------------------------------------------+
|   7 |              MODEL |       :math:`10^{-17} {\rm erg/s/cm}^2`/|ang|/spaxel | The best-fitting model spectra (sum of the fitted continuum and       |
|     |                    |                                                      | emission-line models)                                                 |
+-----+--------------------+------------------------------------------------------+-----------------------------------------------------------------------+
|   8 |         MODEL_MASK |                                                      | The mask from the combined continuum+emission-line model fit          |
+-----+--------------------+------------------------------------------------------+-----------------------------------------------------------------------+
|   9 |             EMLINE |       :math:`10^{-17} {\rm erg/s/cm}^2`/|ang|/spaxel | The model spectrum with *only* the emission lines                     |
+-----+--------------------+------------------------------------------------------+-----------------------------------------------------------------------+
|  10 |            STELLAR |       :math:`10^{-17} {\rm erg/s/cm}^2`/|ang|/spaxel | The best-fitting model spectra fit from the stellar-continuum-only    |
|     |                    |                                                      | fit (used to model the stellar kinematics)                            |
+-----+--------------------+------------------------------------------------------+-----------------------------------------------------------------------+
|  11 |       STELLAR_MASK |                                                      | The mask for the best-fitting model spectra fit from the              |
|     |                    |                                                      | stellar-continuum-only fit (used to model the stellar kinematics)     |
+-----+--------------------+------------------------------------------------------+-----------------------------------------------------------------------+
|  12 |              BINID |                                                      | Numerical ID for spatial bins in 5 channels: (1) binned spectra,      |
|     |                    |                                                      | (2) stellar-continuum results, (3) empty, (4) emission-line model     |
|     |                    |                                                      | results, and (5) empty; i.e., channels 1, 2, and 4 are the same as    |
|     |                    |                                                      | the BINID extension in the ``MAPS`` files and channels 3 and 5 are    |
|     |                    |                                                      | empty.                                                                |
+-----+--------------------+------------------------------------------------------+-----------------------------------------------------------------------+

.. note::

    * The shape and WCS of all extensions with datacubes identically
      match that of the corresponding DRP ``LOGCUBE`` file.

    * To calculate the dereddened flux::

        dereddened_flux = FLUX * REDCORR

----

Special considerations
----------------------

Importantly, please consult the DAP papers (see :ref:`citation`) for
usage guidelines and limitations of the data.

.. _datamodel-binid-usage:

DAP BINIDs and usage
~~~~~~~~~~~~~~~~~~~~

It's important to understand that, for all but the ``SPX`` binning type,
not all of the data in the ``MAPS`` and model ``LOGCUBE`` files are
independent.  Putting aside the issue of :ref:`spatialcovariance`, we
*repeat* measurements for a given binned spectrum in all the spaxels
associated with that bin for consistency between the DAP and DRP data
formats.  Therefore, if you are, e.g., fitting a model to the ``MAPS``
data or calculating azimuthal averages, you should pull out the binned
quantities that are *unique* before proceeding.  In addition to any
associated mask values, you should use the ``BINID`` extension (and,
indeed, its main purpose is) to extract the unique (but still
correlated) data to use in such an analysis.

The ``BINID`` extension has one channel for each of the five main
processing steps: binning, stellar-continuum and -kinematics fitting,
emission-line moment measurements, emission-line Gaussian modeling, and
spectral indices.

Keep in mind the following:

 * ``BINID == -1`` means that the spaxel was *not* included in the
   analysis. For example, ``BINID`` values of -1 in the first
   ``BINID`` channel means that either the spaxel had insufficiently
   good/unmasked pixels or too low S/N to be included in the binning
   procedure. Any spaxel with ``BINID == -1`` should also be masked
   as ``DONOTUSE`` in the respective property map.
 * A ``BINID`` may be :math:`> -1` in one channel and :math:`= -1` in a
   different channel.  For example, a spaxel in the binning ``BINID`` map
   may be :math:`> -1` but -1 in the stellar-continuum ``BINID``.  This
   likely means that the spaxels were successfully binned, but the bin
   had :math:`{\rm S/N} < 1` meaning it was not analyzed by the
   stellar-continuum fitting module.
 * Currently, the only difference in bin IDs is the -1 vs.
   non-negative distinction described in the last point, *except for
   the hybrid binning scheme*. For the ``HYB`` binning case, the
   emission-line moments, emission-line modeling, and spectral-index
   measurements are done on a spaxel-by-spaxel basis, meaning that
   the bin IDs are redetermined and is just a running number (not,
   e.g., ordered by S/N) for the spaxels that were analyzed.

See :ref:`gettingstarted-binid` for usage examples that extracts both
the unique and unmasked data from a ``MAPS`` file to produce the g-band
and :math:`{\rm H}\alpha` surface-brightness profiles.

.. _datamodel-hybrid-binning:

Hybrid (HYB) binning scheme
~~~~~~~~~~~~~~~~~~~~~~~~~~~

In all cases except the ``HYB`` binning approach, each analysis
module only works with the "binned" spectra after the binning is
performed. (I've put "binned" in quotes here because all spectra are
treated the same after the binning step, even if the "bin" only
includes a single spaxel.) In the ``HYB`` case, the emission-line
modeling is done by first fitting the continuum+emission-line data
simultaneously, distributing those results as a starting point for
fitting the spaxels within the bin, and then redoing the simultaneous
fit for each spaxel. By fitting the data as a hybrid between the
``VOR10`` and ``SPX`` binning schemes, there are a few things to keep
in mind:

 * Because the stellar kinematics are held fixed to the binned
   results during the spaxel-by-spaxel continuum+emission-line fit,
   there will be (subtle) spatial covariance issues between spaxels
   associated with a single bin, beyond the :ref:`spatialcovariance`
   from the datacube construction alone.

 * The binned spectra provided in the ``HYB`` model ``LOGCUBE`` files
   are from the Voronoi binning step; however, the emission-line
   models are fit to the *individual spaxels*. When using the model
   ``LOGCUBE`` files for this binning scheme:
   
    * The stellar-continuum fits (in the ``STELLAR`` extension) should
      be compared to the Voronoi binned spectra in the file;
    * **however**, the best-fitting model spectra (stellar continuum +
      gas emission) in the ``MODEL`` extension should be compared to the
      individual spectra from the *DRP LOGCUBE* file!
      
 * Because the emission-line modeling is done on the individual spaxels,
   the emission-line moments are recalculated after the emission-line
   modeling to ensure the stellar continuum used for both the Gaussian
   model and the moment calculation is identical.  In the ``HYB`` case,
   this means the emission-line moments are also provided for the
   individual spaxels.

 * The spectral indices are measured on the individual spaxels because
   the emission-line model is first subtracted from the data before the
   index measurements.

Usage Guidelines: Stellar velocity dispersions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Measurement of stellar (and gas!) velocity dispersions in MaNGA is
complicated by the spectral resolution, particularly at low S/N and
low :math:`\sigma`. Please tread carefully! In particular, please
consult Section 7.7 of `Westfall et al. (2019, AJ, 158, 231)`_ for a
detailed discussion of best practices for the stellar velocity
dispersion data.

In summary, there is no hard and fast rule along the lines of, "Only
use measurements when the S/N is above X". (In fact, having
measurements at the lower S/N level is useful for understanding the
affects of the error distribution.) However, here are some rough
guidelines to consider when handling the velocity dispersion data:

 * Kinematics should smoothly vary between adjacent spaxels
 * All velocities are statistically well behaved, except possibly at
   :math:`{\rm S/N} < 5` for :math:`\sigma \sim \sigma_{\rm inst}/2`
 * Be aware of the *distribution* of :math:`\sigma` at a given radius
   or surface brightness when assessing the data.
 * Don’t trust single :math:`\sigma` measurements at :math:`{\rm
   S/N}<5`, only use them to understand the error distribution.
 * Systematic errors in individual :math:`\sigma` become appreciable at:

    * :math:`{\rm S/N} < 20` for :math:`\sigma \sim \sigma_{\rm inst}/2`
      (:math:`\sim 35` km/s)
    * :math:`{\rm S/N} < 10` for :math:`\sigma \sim \sigma_{\rm inst}`
      (:math:`\sim 70` km/s)

.. _datamodel-eml-tpl-resolution:

Usage Guidelines: Emission-line template resolution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When using the recommended emission-line module
(:class:`~mangadap.proc.sasuke.Sasuke`), the emission lines are fit
in a very similar way to the stellar continuum using a set of
emission-line templates. Given the varying spectral resolution of the
MaNGA data, we setup these templates to have a non-zero "instrumental
dispersion" that is the same as the MaNGA data up to some quadrature
offset. The value of the "template instrumental dispersion" at the
location of each emission line is provided in the ``EMLINE_TPLSIGMA``
extension of the ``MAPS`` files. The velocity dispersion actually
measured by this emission-line module (using pPXF) is the quadrature
difference between the template dispersion and the directly observed
sigma of the emission-line (as fit by a Gaussian).

To keep things consistent between MPLs and provide what people expect,
the ``EMLINE_GSIGMA`` data provide the sigma of the line as it would be if
measured by a direct fit of a Gaussian to the line; i.e., we add back
the template instrumental dispersion in quadrature to the pPXF-fitted
sigma and propagate the error as follows:

    - :math:`\sigma^2 = \sigma_{\rm ppxf}^2 + \sigma_{\rm tpl}^2`
    - :math:`\epsilon[\sigma] = \sigma_{\rm ppxf} \epsilon[\sigma_{\rm
      ppxf}]/\sigma`

The ``EMLINE_TPLSIGMA`` (:math:`\sigma_{\rm tpl}`) extension is provided
so that one can recover the exact output from pPXF following the
equations above, where :math:`\sigma` and
:math:`(\epsilon[\sigma])^{-2}` are provided in ``EMLINE_GSIGMA`` and
``EMLINE_GSIGMA_IVAR``, respectively.  One does *not* need to consider
``EMLINE_TPLSIGMA`` when calculating the astrophysical Doppler
broadening of each line; see :ref:`corrections`.

----

DAP global header data
----------------------

The first extension of each of the main DAP output files (the
``MAPS`` and model ``LOGCUBE``) is empty apart from the header data.
The header data is an exact copy of the primary header for the `DRP
LOGCUBE`_ file except that the ``BSCALE``, ``BZERO``, and ``BUNIT``
keywords are removed and the ``AUTHOR`` and ``MASKNAME`` keywords are
changed.

The following keywords are also added, any keyword enclose in
() are only written under certain conditions:

+------------+--------------------------------------------------------------------------------------+
|    Keyword | Description                                                                          |
+============+======================================================================================+
| VERSPY     | `Python <https://www.python.org/>`_ version                                          |
+------------+--------------------------------------------------------------------------------------+
| VERSNP     | `Numpy <https://www.numpy.org/>`_  version                                           |
+------------+--------------------------------------------------------------------------------------+
| VERSSCI    | `Scipy <https://www.scipy.org/>`_ version                                            |
+------------+--------------------------------------------------------------------------------------+
| VERSAST    | `Astropy <https://www.astropy.org/>`_ version                                        |
+------------+--------------------------------------------------------------------------------------+
| VERSPYDL   | `pydl <https://pydl.readthedocs.io/>`_ version                                       |
+------------+--------------------------------------------------------------------------------------+
| VERSDAP    | MaNGA DAP version                                                                    |
+------------+--------------------------------------------------------------------------------------+
| DAPTYPE    | The analysis method identifier for the DAP analysis                                  |
|            | (e.g., ``HYB10-MILESHC-MASTARHC2``)                                                  |
+------------+--------------------------------------------------------------------------------------+
| DAPFRMT    | The format of this output file, either ``MAPS`` or ``LOGCUBE``                       |
+------------+--------------------------------------------------------------------------------------+
| RDXQAKEY   | Configuration keyword for the method used to assess the reduced data                 |
+------------+--------------------------------------------------------------------------------------+
| ECOOPA     | Position angle used for the semi-major axis polar coordinate calculations            |
+------------+--------------------------------------------------------------------------------------+
| ECOOELL    | Ellipticity (1-b/a) used for the semi-major axis polar coordinate calculations       |
+------------+--------------------------------------------------------------------------------------+
| BBWAVE     | Wavelength of the ``LOGCUBE`` channel used for calculating the covariance used in    |
|            | the per spaxel S/N calculation                                                       |
+------------+--------------------------------------------------------------------------------------+
| BBINDX     | Index of the channel                                                                 |
+------------+--------------------------------------------------------------------------------------+
| REFF       | Effective radius                                                                     |
+------------+--------------------------------------------------------------------------------------+
| BINKEY     | Configuration keyword for the spatial binning method                                 |
+------------+--------------------------------------------------------------------------------------+
| BINMINSN   | Minimum S/N of spectrum to include in the binning                                    |
+------------+--------------------------------------------------------------------------------------+
| FSPCOV     | Minimum allowed fraction of good pixels across the full spectral range               |
+------------+--------------------------------------------------------------------------------------+
| NBINS      | Number of unique spatial bins                                                        |
+------------+--------------------------------------------------------------------------------------+
| (EMPTYBIN) | List of empty bins, if any exist                                                     |
+------------+--------------------------------------------------------------------------------------+
| BINTYPE    | Spatial binning method                                                               |
+------------+--------------------------------------------------------------------------------------+
| (BINCX)    | If radial binning, on-sky X center for all bins                                      |
+------------+--------------------------------------------------------------------------------------+
| (BINCY)    | If radial binning, on-sky Y center for all bins                                      |
+------------+--------------------------------------------------------------------------------------+
| (BINPA)    | If radial binning, position angle used for all bins                                  |
+------------+--------------------------------------------------------------------------------------+
| (BINELL)   | If radial binning, ellipticity (1-b/a) used for all bins                             |
+------------+--------------------------------------------------------------------------------------+
| (BINSCL)   | If radial binning, the radius has been scaled by this value (arcsec)                 |
+------------+--------------------------------------------------------------------------------------+
| (BINRAD)   | If radial binning, provides the start, end, and number of radial bins                |
+------------+--------------------------------------------------------------------------------------+
| (BINLGR)   | If radial binning, the geometric step used to set the radial bins                    |
+------------+--------------------------------------------------------------------------------------+
| (BINSNR)   | If Voronoi binning, the target S/N for each bin                                      |
+------------+--------------------------------------------------------------------------------------+
| (BINCOV)   | If Voronoi binning, the method used to incorporate covariance into the S/N           |
|            | calculation                                                                          |
+------------+--------------------------------------------------------------------------------------+
| (NCALIB)   | If Voronoi binning and using a calibration of the noise vector that incorporates     |
|            | covariance, the noise calibration coefficient                                        |
+------------+--------------------------------------------------------------------------------------+
| (STCKOP)   | If binning spectra, the operation used for stacking spectra                          |
+------------+--------------------------------------------------------------------------------------+
| (STCKVREG) | If binning spectra, a boolean flag that the spectra were shifted in velocity before  |
|            | stacked                                                                              |
+------------+--------------------------------------------------------------------------------------+
| (STCKCRMD) | If binning spectra, the approach used to account for covariance in the resulting     |
|            | inverse variance of the binned spectrum                                              |
+------------+--------------------------------------------------------------------------------------+
| (STCKCRPR) | If binning spectra, the method-specific parameters used to incorporate covariance in |
|            | the stacking procedure                                                               |
+------------+--------------------------------------------------------------------------------------+
| (STCKRES)  | Stacking operation performs a stack of the individual spaxel resolution vectors      |
|            | (DISP) as opposed to the single median vector (SPECRES)                              |
+------------+--------------------------------------------------------------------------------------+
| (STCKPRE)  | Stacking operation uses the pre-pixelized spectral resolution instead of the         |
|            | post-pixelized version                                                               |
+------------+--------------------------------------------------------------------------------------+
| GEXTLAW    | Galactic extinction law used to deredden the data                                    |
+------------+--------------------------------------------------------------------------------------+
| RVGAL      | Ratio of total to selective extinction, :math:`R_V`                                  |
+------------+--------------------------------------------------------------------------------------+
| VSTEP      | Velocity step per spectral channel                                                   |
+------------+--------------------------------------------------------------------------------------+
| SCKEY      | Configuration keyword for the method used to model the stellar-continuum             |
+------------+--------------------------------------------------------------------------------------+
| SCMINSN    | Minimum S/N of spectrum to include in stellar-continuum fits                         |
+------------+--------------------------------------------------------------------------------------+
| SCINPVEL   | Initial guess stellar velocity                                                       |
+------------+--------------------------------------------------------------------------------------+
| SCINPSIG   | Initial guess stellar velocity dispersion                                            |
+------------+--------------------------------------------------------------------------------------+
| NSCMOD     | Number of unique stellar-continuum models                                            |
+------------+--------------------------------------------------------------------------------------+
| (EMPTYSC)  | List of bins without a stellar-continuum model, if any exist                         |
+------------+--------------------------------------------------------------------------------------+
| SCTYPE     | Type of spectral fitting method used for the stellar-continuum fits                  |
+------------+--------------------------------------------------------------------------------------+
| SCMETH     | Algorithm used for the stellar-continuum fits                                        |
+------------+--------------------------------------------------------------------------------------+
| PPXFTPLK   | Configuration keyword for the template library key used with pPXF                    |
+------------+--------------------------------------------------------------------------------------+
| PPXFBIAS   | pPXF bias value                                                                      |
+------------+--------------------------------------------------------------------------------------+
| PPXFMOM    | Number of fitted LOSVD moments in pPXF                                               |
+------------+--------------------------------------------------------------------------------------+
| PPXFAO     | Order of additive polynomial in pPXF                                                 |
+------------+--------------------------------------------------------------------------------------+
| PPXFMO     | Order of multiplicative polynomial in pPXF                                           |
+------------+--------------------------------------------------------------------------------------+
| PPXFRBOX   | Size of the boxcar filter used during rejection iterations                           |
+------------+--------------------------------------------------------------------------------------+
| ELMKEY     | Configuration keyword that defines the emission-line moment measurement method       |
+------------+--------------------------------------------------------------------------------------+
| ELMMINSN   | Minimum S/N of spectrum to include in emission-line moment measurements              |
+------------+--------------------------------------------------------------------------------------+
| ARTDB      | Artifact database keyword                                                            |
+------------+--------------------------------------------------------------------------------------+
| MOMDB      | Emission-line moments database keyword                                               |
+------------+--------------------------------------------------------------------------------------+
| ELFKEY     | Configuration keyword that defines the emission-line modeling method                 |
+------------+--------------------------------------------------------------------------------------+
| ELFMINSN   | Minimum S/N of spectrum to include in emission-line modeling                         |
+------------+--------------------------------------------------------------------------------------+
| EMLDB      | Emission-line database keyword                                                       |
+------------+--------------------------------------------------------------------------------------+
| NELMOD     | Number of unique emission-line models                                                |
+------------+--------------------------------------------------------------------------------------+
| ELTYPE     | Type of spectral fitting method used for the emission-line fits                      |
+------------+--------------------------------------------------------------------------------------+
| ELMETH     | Algorithm used for the emission-line modeling                                        |
+------------+--------------------------------------------------------------------------------------+
| SIKEY      | Configuration keyword that defines the spectral-index measurement method             |
+------------+--------------------------------------------------------------------------------------+
| SIMINSN    | Minimum S/N of spectrum to include in spectral-index measurements                    |
+------------+--------------------------------------------------------------------------------------+
| SIFWHM     | FWHM of index system resolution (ang) to which the galaxy spectra were matched       |
+------------+--------------------------------------------------------------------------------------+
| ABSDB      | Absorption-line index database keyword                                               |
+------------+--------------------------------------------------------------------------------------+
| BHDDB      | Bandhead-index database keyword                                                      |
+------------+--------------------------------------------------------------------------------------+
| SICORR     | Flag that indices have been corrected for velocity dispersion                        |
+------------+--------------------------------------------------------------------------------------+
| SNRGMED    | Median g-band signal-to-noise of spaxels within 1-1.5 :math:`R_e`                    |
+------------+--------------------------------------------------------------------------------------+
| SNRGRING   | Total g-band signal-to-noise of a binned spectrum using spaxels within 1-1.5         |
|            | :math:`R_e` bin                                                                      |
+------------+--------------------------------------------------------------------------------------+
| SNRRMED    | Median r-band signal-to-noise of spaxels within 1-1.5  :math:`R_e`                   |
+------------+--------------------------------------------------------------------------------------+
| SNRRRING   | Total r-band signal-to-noise of a binned spectrum using spaxels within 1-1.5         |
|            | :math:`R_e` bin                                                                      |
+------------+--------------------------------------------------------------------------------------+
| SNRIMED    | Median i-band signal-to-noise of spaxels within 1-1.5  :math:`R_e`                   |
+------------+--------------------------------------------------------------------------------------+
| SNRIRING   | Total i-band signal-to-noise of a binned spectrum using spaxels within 1-1.5         |
|            | :math:`R_e` bin                                                                      |
+------------+--------------------------------------------------------------------------------------+
| SNRZMED    | Median z-band signal-to-noise of spaxels within 1-1.5  :math:`R_e`                   |
+------------+--------------------------------------------------------------------------------------+
| SNRZRING   | Total z-band signal-to-noise of a binned spectrum using spaxels within 1-1.5         |
|            | :math:`R_e` bin                                                                      |
+------------+--------------------------------------------------------------------------------------+
| DAPQUAL    | Global DAP quality bit mask: :ref:`metadatamodel-dapqual`                            |
+------------+--------------------------------------------------------------------------------------+

The headers of the data extensions are more minimal. They include:

 * the WCS information,
 * the :ref:`datamodel-hduclass` keyword block,
 * the channel description for the :ref:`datamodel-maps` files,
 * the units for any single image or datacube extensions (``BUNIT``),
   and
 * the ``DATASUM`` and ``CHECKSUM`` values.

----

Reference Files
---------------

For storage of many more fitting products (so far not deemed useful for
the ``MAPS`` files) and rerunning the code, intermediate reference files
are written after each main analysis step.  The naming convention is
essentially to append the necessary analysis keyword to the file name.
These are identically the keys used in the
:ref:`execution-analysis-plan` file: ``drpqa_key``, ``bin_key``,
``continuum_key``, ``elmom_key``, ``elfit_key``, ``spindex_key``.

The reference files are primarily for developer use, but may contain
information that you want.  A bare-bones description of the content of
these files is forthcoming.  If you're interested in using something in
these files, it's probably best to `Submit an issue`_.

