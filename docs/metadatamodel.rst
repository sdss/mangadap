
.. _metadatamodel:

DAP Metadata Model
==================

.. _metadatamodel-maskbits:

Maskbits
--------

The mask bits for the DAP output files are pulled from IDLUTILS. A
default file is provided with the ``mangadap`` distribution; however,
the DAP ``setup.py`` file will also attempt to download the most
relevant file from the public SDSS SVN repository using the python
`requests library <https://pypi.org/project/requests/>`_.

If you're unfamiliar with maskbits, see the `SDSS Primer
<http://www.sdss.org/dr15/algorithms/bitmasks/>`_

See :ref:`bitmasks` for examples of how to use the mask bits.

.. note::

    Some mask bits are *informational*, not necessarily indicating that
    the pixel should be ignored.  Please read and understand the flags
    listed below to determine if it's reasonable for your science case
    to simply ignore pixels with a non-zero mask value.

    Some mask bits are still *placeholders*.  They're notional bits that
    are actually never set in the current version of the DAP.

.. _metadatamodel-dapqual:

MANGA_DAPQUAL
~~~~~~~~~~~~~

There is a single summary maskbit ``MANGA_DAPQUAL`` included in the
headers of both the MAPS and LOGCUBE files describing the overall
quality of the data:

+-----+-----------+--------------------------------------------------------------+
| Bit | Name      | Description                                                  |
+=====+===========+==============================================================+
|   0 | FORESTAR  | There is a foreground star within the data cube              |
+-----+-----------+--------------------------------------------------------------+
|   1 | BADZ      | NSA redshift does not match derived redshift (*placeholder*) |
+-----+-----------+--------------------------------------------------------------+
|   2 | LINELESS  | No emission lines in data cube (*placeholder*)               |
+-----+-----------+--------------------------------------------------------------+
|   3 | PPXFFAIL  | PPXF fails to fit this object (*placeholder*)                |
+-----+-----------+--------------------------------------------------------------+
|   4 | SINGLEBIN | Voronoi binning resulted in all spectra in a single bin      |
+-----+-----------+--------------------------------------------------------------+
|  28 | DRPCRIT   | Critical failure in DRP                                      |
+-----+-----------+--------------------------------------------------------------+
|  29 | DAPCRIT   | Critical failure in DAP                                      |
+-----+-----------+--------------------------------------------------------------+
|  30 | CRITICAL  | Critical failure in DRP or DAP                               |
+-----+-----------+--------------------------------------------------------------+

Anything with the CRITICAL bit set in MANGA_DAPQUAL should generally not
be used for scientific purposes.

.. _metadatamodel-dappixmask:

MANGA_DAPPIXMASK
~~~~~~~~~~~~~~~~

``MANGA_DAPPIXMASK`` is the 2D image bitmap used to describe the quality
of individual pixel measurements in the DAP ``MAPS`` file.  It can take
values:

+-----+--------------+---------------------------------------------------------------+
| Bit | Name         | Description                                                   |
+=====+==============+===============================================================+
|  0  | NOCOV        | No coverage in this spaxel                                    |
+-----+--------------+---------------------------------------------------------------+
|  1  | LOWCOV       | Low coverage in this spaxel                                   |
+-----+--------------+---------------------------------------------------------------+
|  2  | DEADFIBER    | Major contributing fiber is dead                              |
+-----+--------------+---------------------------------------------------------------+
|  3  | FORESTAR     | Foreground star                                               |
+-----+--------------+---------------------------------------------------------------+
|  4  | NOVALUE      | Spaxel was not fit because it did not meet selection criteria |
+-----+--------------+---------------------------------------------------------------+
|  5  | UNRELIABLE   | Value is deemed unreliable; see definition below              |
+-----+--------------+---------------------------------------------------------------+
|  6  | MATHERROR    | Mathematical error in computing value                         |
+-----+--------------+---------------------------------------------------------------+
|  7  | FITFAILED    | Attempted fit for property failed                             |
+-----+--------------+---------------------------------------------------------------+
|  8  | NEARBOUND    | Fitted value is too near an imposed boundary                  |
+-----+--------------+---------------------------------------------------------------+
|  9  | NOCORRECTION | Appropriate correction not available                          |
+-----+--------------+---------------------------------------------------------------+
| 10  | MULTICOMP    | Multi-component velocity features present (*placeholder*)     |
+-----+--------------+---------------------------------------------------------------+
| 30  | DONOTUSE     | Do not use this spaxel for science                            |
+-----+--------------+---------------------------------------------------------------+

The most important flags are incorporated into the ``DONOTUSE`` bit,
which indicates that a given pixel should not be used for science.

.. _metadatamodel-nearbound:

The NEARBOUND flag
++++++++++++++++++

The ``NEARBOUND`` flag is used to signify that a returned parameter is
likely biased by an imposed boundary on the allowed parameter space.
These are specifically relevant to the pPXF kinematics from
:class:`mangadap.proc.ppxffit.PPXFFit` (stellar) and
:class:`mangdap.proc.sasuke.Sasuke` (gas).  We use pPXF with a
:math:`\pm 2000` km/s limit relative to the input guess velocity (set by
:math:`cz` in ``SCINPVEL`` header keyword in the ``PRIMARY`` extension
and most often identical to the NSA redshift) on the returned velocity
and :math:`{\rm d}v/100 < \sigma < 1000` limit on the velocity
dispersion, where :math:`{\rm d}v` is the size of the MaNGA LOGCUBE
wavelength channel (:math:`\sim 70` km/s; given by the ``VSTEP`` header
keyword in the ``PRIMARY`` extension.  The returned velocity is
determined to be ``NEARBOUND`` if the "best-fit" value is within 1/100
of the width of the allowed range of either boundary; i.e.,
``NEARBOUND`` is triggered if the velocity is :math:`-2000<v<-1980` or
:math:`1980<v<2000`.  For the velocity dispersion, ``NEARBOUND`` is
triggered by the same criterion but geometrically; i.e., if the velocity
dispersion is :math:`0.69 < \sigma < 0.74` or :math:`929.8 < \sigma <
1000`.

.. _metadatamodel-unreliable:

The UNRELIABLE flag
+++++++++++++++++++

The ``UNRELIABLE`` flag is *not* incorporated into the ``DONOTUSE``
flag.  This flag is tripped based on various judgement calls made by the
data products committee (DPC) and the pipeline development teams.  You
are **strongly** encouraged to understand the implications of this flag
on the data and how to properly make the distinction between the
``DONOTUSE`` and ``UNRELIABLE`` flags for your science application.  The
definition of the ``UNRELIABLE`` flag can change with time, in the hope
that we eventually converge to a refined set of criteria that allow
users to determine when measurements can be trusted carte blanche and
when the data should be treated more skeptically.  Only spaxels where
analysis has been attempted (with non-zero bin IDs) are flagged as
``UNRELIABLE`` if they meet the necessary criteria.  *Please let us know
if you find a set of automated criteria that would be useful to the
development team in terms of what you would like to see marked as*
``UNRELIABLE``.

Currently, the use of the ``UNRELIABLE`` flag is still rather limited.
This is not to say that all measurements are reliable, but reflects our
hesitance to set (robust) criteria for isolating unreliable
measurements, either because we don't think we're able to or because we
haven't had sufficient time to do so.  Below, we list the condition
under which ``UNRELIABLE`` flags are tripped, and the affected masks in
the ``MAPS`` file.

+-------------------+---------------------------------------------------+
| Affected mask     | Criteria                                          |
+===================+===================================================+
| EMLINE_SFLUX_MASK | If there are '''any''' masked pixels in the three |
|                   | passbands (blue, main, red) used to construct the |
|                   | measurement.                                      |
+-------------------+---------------------------------------------------+
| EMLINE_SEW_MASK   | If there are '''any''' masked pixels in the three |
|                   | passbands (blue, main, red) used to construct the |
|                   | measurement.                                      |
+-------------------+---------------------------------------------------+
| SPECINDEX_MASK    | If there are '''any''' masked pixels in the three |
|                   | passbands (blue, main, red) used to construct the |
|                   | measurement.                                      |
+-------------------+---------------------------------------------------+

.. _metadatamodel-dapspecmask:

MANGA_DAPSPECMASK
~~~~~~~~~~~~~~~~~

``MANGA_DAPPIXMASK`` is the 3D model cube bitmask used to describe the
quality of individual spaxel fits in the DAP model data cube file.  It
can take values:

+-----+--------------+--------------------------------------------------------+
| Bit | Name         | Description                                            |
+=====+==============+========================================================+
|  0  | IGNORED      | Ignored because of DRP flags or stacking parameters    |
+-----+--------------+--------------------------------------------------------+
|  1  | FORESTAR     | There is a foreground star within the data cube        |
+-----+--------------+--------------------------------------------------------+
|  2  | FLUXINVALID  | Ignored because (stacked) flux not valid               |
+-----+--------------+--------------------------------------------------------+
|  3  | IVARINVALID  | Ignored because inverse variance not valid             |
+-----+--------------+--------------------------------------------------------+
|  4  | ARTIFACT     | Contains a designated artifact                         |
+-----+--------------+--------------------------------------------------------+
|  5  | FITIGNORED   | Ignored during fit                                     |
+-----+--------------+--------------------------------------------------------+
|  6  | FITFAILED    | Fit failed in this region                              |
+-----+--------------+--------------------------------------------------------+
|  7  | ELIGNORED    | Ignored during emission-line fit (**deprecated**)      |
+-----+--------------+--------------------------------------------------------+
|  8  | ELFAILED     | Emission-line fit failed (**deprecated**)              |
+-----+--------------+--------------------------------------------------------+
|  9  | NOMODEL      | Identifies pixels outside of the fitted spectral range |
+-----+--------------+--------------------------------------------------------+

DAP Execution Files
-------------------

The DAP is configured using an input execution plan file created by the
user.  There are additional intermediary script files created by the DAP
to allow for event handling and cluster coordination.

See :ref:`execution` for more general information about execution of
the DAP.  What follows is specifically for the survey level execution of
the DAP.

AnalysisPlan file
~~~~~~~~~~~~~~~~~

For a general description the ``AnalysisPlan`` file, see
:ref:`execution-analysis-plan`.

*File template*:
``$MANGA_SPECTRO_ANALYSIS/$MANGADRP_VER/$MANGADAP_VER/log/[timestamp]/mpl[MPL]_plan.par``

In the file template, ``[timestamp]`` is the time when the ``rundap``
script was executed in the format, e.g., ``01Nov2019T16.58.40UTC`` and
``[MPL]`` is the MPL number (e.g., 9).  This is a single file that lists
the ways in which each DRP ``LOGCUBE`` file is to be analyzed for each
MPL.  This file is created once by the person executing the DAP.

.. _metadatamodel-drpcomplete:

DRPComplete database
~~~~~~~~~~~~~~~~~~~~

*File template*:
``$MANGA_SPECTRO_ANALYSIS/$MANGADRP_VER/$MANGADAP_VER/common/drpcomplete_$MANGADRP_VER.fits``

The :class:`mangadap.survey.drpcomplete.DRPComplete` file is primarily
created for the survey-level execution of the DAP.  It collates
information used to create the input parameter files for each completed
DRP file.  It is created/updated at the beginning of each
:ref:`execution-rundap`.

The :class:`mangadap.survey.drpcomplete.DRPComplete` database is written
to a fits file with a primary extension and a binary-table extension;
the table extension has the following columns:

+-------------------+----------------------------------------------------+
| Column            | Description                                        |
+===================+====================================================+
|         ``PLATE`` | The plate number of the observation                |
+-------------------+----------------------------------------------------+
|     ``IFUDESIGN`` | The IFU used to observe the target                 |
+-------------------+----------------------------------------------------+
|         ``MODES`` | Specifies which DRP files are available on disk:   |
|                   | (1) ``CUBE`` files only; (2) Both ``CUBE`` and     |
|                   | ``RSS`` files.                                     |
+-------------------+----------------------------------------------------+
|       ``MANGAID`` | String representation of the MaNGA ID              |
+-------------------+----------------------------------------------------+
|         ``OBJRA`` | Nominal right ascension of the object center       |
+-------------------+----------------------------------------------------+
|        ``OBJDEC`` | Nominal declination of the object center           |
+-------------------+----------------------------------------------------+
|         ``CATID`` | ID number of the parent catalog                    |
+-------------------+----------------------------------------------------+
|       ``CATINDX`` | 0-based index of the row within that parent        |
|                   | catalog with the target information                |
+-------------------+----------------------------------------------------+
|   ``TRG_VERSION`` | The version of the parent catalog                  |
+-------------------+----------------------------------------------------+
|        ``TRG_ID`` | The ID number of the object in the parent catalog. |
+-------------------+----------------------------------------------------+
| ``MANGA_TARGET1`` | Targeting bits for main survey galaxies            |
+-------------------+----------------------------------------------------+
| ``MANGA_TARGET3`` | Targeting bits for ancillary programs              | 
+-------------------+----------------------------------------------------+
|           ``VEL`` | Redshift (:math:`cz`) of the object used as an     |
|                   | initial guess redshift.                            |
+-------------------+----------------------------------------------------+
|         ``VDISP`` | Characteristic velocity dispersion of the object   |
+-------------------+----------------------------------------------------+
|           ``ELL`` | Photometric ellipticity                            |
+-------------------+----------------------------------------------------+
|            ``PA`` | Photometric position angle                         |
+-------------------+----------------------------------------------------+
|          ``REFF`` | Effective (half-light) radius                      |
+-------------------+----------------------------------------------------+

.. note::

    - The DAP currently only works with the ``LOG`` format, and does not
      search for or analyze the ``LIN`` format.

    - ``OBJRA`` and ``OBJDEC`` are not necessarily located at the center
      of the IFU field of view.  The IFU center coordinates are provided
      in DRPall file (**link**) as ``IFURA`` and ``IFUDEC``.

    - The MaNGA ID is defined as ``[CATID]-[CATINDX]`` (**link**)

    - For the main survey galaxies, ``TRG_VERSION`` and ``TRG_ID`` are
      drawn from the NASA-Sloan atlas and are identical to
      'nsa_nsa_version' and 'nsa_nsaid' in the DRPall file (**link**).

    - The targeting bits are defined (**link**).

    - The redshifts from the NSA and ancillary-program catalogs are
      consolidated into the 'z' column in the DRPall file.  See
      discussion of the "redshift fix file" below.

    - The characteristic velocity dispersion is virtually always not
      available and set to -9999.  In this case, the DAP will default to
      a dispersion of 100 km/s.

    - For main survey galaxies, photometric measurements are taken from
      'nsa_ba', 'nsa_phi' and 'nsa_petro_th50_el' in the DRPall file.
      If any of these values do not exist or are 'nan', they are set to
      -9999.0.  Importantly, *these placeholder values are replaced by
      ``ELL=0; PA=0; REFF=1`` when processed by the DAP.*


.. _metadatamodel-redshift-fix:

Redshift Fix File
~~~~~~~~~~~~~~~~~

*File template*:
``$MANGADAP_DIR/data/fix/redshift_fix.par``

The redshift-fix file is an `SDSS parameter file
<https://www.sdss.org/dr15/software/par/>`_ used to replace any redshift
(:math:`cz`) read from the DRPall or plateTargets files.  It has a
simple format that identifies the plate, ifudesign, and replacement
redshift:

.. code-block:: c

    typedef struct {
        int plate;
        int ifudesign;
        double z;
    } DAPZCORR;
    
    DAPZCORR  9677  6102 0.0
    DAPZCORR  9677  6103 0.0
    ...

This files serves to both provide redshifts for objects that don't have
them and replace incorrect redshifts from, e.g., the NASA-Sloan Atlas.
The redshift-fix file is updated for each version of the DAP.

Execution Script
~~~~~~~~~~~~~~~~

*File template*:
``$MANGA_SPECTRO_ANALYSIS/$MANGADRP_VER/$MANGADAP_VER/log/[timestamp]/[PLATE]/[IFUDESIGN]/mangadap-[PLATE]-[IFUDESIGN]``

In the file template, ``[timestamp]`` is the time when the ``rundap``
script was executed in the format, e.g., ``01Nov2019T16.58.40UTC``,
``[PLATE]`` is the plate number and ``[IFUDESIGN]`` is the IFU number.  These
are the script files that are submitted to the Utah CHPC cluster to
execute the DAP, as created by the ``rundap`` script

See :ref:`execution-rundap`.

Observational parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For a general description the ``ObsInputPar`` file, see
:ref:`execution-obs-input`.

*File template*: ``$MANGA_SPECTRO_ANALYSIS/$MANGADRP_VER/$MANGADAP_VER/common/[PLATE]/[IFUDESIGN]/mangadap-[PLATE]-[IFUDESIGN]-LOG[MODE].par``

*Symlink*: ``$MANGA_SPECTRO_ANALYSIS/$MANGADRP_VER/$MANGADAP_VER/[DAPTYPE]/[PLATE]/[IFUDESIGN]/ref/mangadap-[PLATE]-[IFUDESIGN]-LOG[MODE].par``

In the file templates, ``[PLATE]`` is the plate number, ``[IFUDESIGN]``
is the IFU number, ``[MODE]`` is the data format (always ``CUBE`` for
now), and ``[DAPTYPE]`` is the keyword for the analysis approach.  These
files provide input observational parameters to the DAP, and are almost
entirely from the NASA-Sloan Atlas.

.. _metadatamodel-dapall:

DAPall database
---------------

*File template*:
``$MANGA_SPECTRO_ANALYSIS/$MANGADRP_VER/$MANGADAP_VER/dapall-$MANGADRP_VER-$MANGADAP_VER.fits``

The DAPall file has two extensions:
    
    1. ``PRIMARY``: Empty apart from the header information
    2. ``DAPALL``: Binary table data

The DAPall file contains one row per DAP ``MAPS`` file, such that the
total number of rows is :math:`N_{\rm cube}*N_{\rm daptype}`.

Header data
~~~~~~~~~~~

The ``PRIMARY`` extension is empty apart from the following header
keywords:

+--------------+-------------------------------------------------------+
| Key          | Comment                                               |
+==============+=======================================================+
| ``VERSDRP3`` | DRP version                                           |
+--------------+-------------------------------------------------------+
|  ``VERSDAP`` | DAP version                                           |
+--------------+-------------------------------------------------------+
|   ``ELS[n]`` | Line name for non-parametric (summed) emission-line   |
|              | measurement at vector position ``n``-1 in relevant    |
|              | columns of the database                               |
+--------------+-------------------------------------------------------+
|   ``ELG[n]`` | Line name for Gaussian emission-line measurement at   |
|              | vector position ``n``-1 in relevant columns of the    |
|              | database                                              |
+--------------+-------------------------------------------------------+
|   ``SPI[n]`` | Name of spectral index measurement at vector position |
|              | ``n``-1 in relevant columns of the database           |
+--------------+-------------------------------------------------------+
|  ``SPIU[n]`` | Unit of the spectral index measurement at vector      |
|              | position ``n``-1 in relevant columns of the database  |
+--------------+-------------------------------------------------------+
| ``CHECKSUM`` | Used for checking data fidelity                       |
+--------------+-------------------------------------------------------+
|  ``DATASUM`` | Used for checking data fidelity                       |
+--------------+-------------------------------------------------------+

Binary table data
~~~~~~~~~~~~~~~~~

The binary table in the ``DAPALL`` extension has the following columns:

+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|                        Key |            Type |                                              Units | Comment                                                                       |
+============================+=================+====================================================+===============================================================================+
| **Basic designation and NSA information**                                                                                                                                         |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|                  ``PLATE`` |             int |                                                    | Plate number                                                                  |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|              ``IFUDESIGN`` |             int |                                                    | IFU design number                                                             |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|               ``PLATEIFU`` |             str |                                                    | String combination of ``[PLATE]-[IFUDESIGN]`` to ease searching               |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|                ``MANGAID`` |             str |                                                    | MaNGA ID string                                                               |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|             ``DRPALLINDX`` |             int |                                                    | Row index of the observation in the DRPall file                               |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|                   ``MODE`` |             str |                                                    | 3D mode of the DRP file (``CUBE`` or ``RSS``)                                 |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|                ``DAPTYPE`` |             str |                                                    | Keyword of the analysis approach used (e.g., ``SPX-MILESHC-MASTARHC``)        |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|                ``DAPDONE`` |            bool |                                                    | Flag that MAPS file successfully produced                                     |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|                  ``OBJRA`` |          double |                                                deg | RA of the galaxy center                                                       |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|                 ``OBJDEC`` |          double |                                                deg | Declination of the galaxy center                                              |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|                  ``IFURA`` |          double |                                                deg | RA of the IFU pointing center (generally the same as  ``OBJRA``)              |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|                 ``IFUDEC`` |          double |                                                deg | Declination of the IFU pointing center (generally the same as  ``OBJDEC``)    |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|               ``MNGTARG1`` |             int |                                                    | Main survey targeting bit (**link**)                                          |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|               ``MNGTARG2`` |             int |                                                    | Non-galaxy targeting bit (**link**)                                           |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|               ``MNGTARG3`` |             int |                                                    | Ancillary targeting bit (**link**)                                            |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|                      ``Z`` |          double |                                                    | Redshift used for initial guess velocity (typically identical to ``NSA_Z``)   |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|                ``LDIST_Z`` |          double |                           :math:`h^{-1} {\rm Mpc}` | Luminosity distance based on  ``Z`` and a standard cosmology                  |
|                            |                 |                                                    | (:math:`h=1; \Omega_M=0.3; \Omega_\Lambda=0.7`)                               |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|                ``ADIST_Z`` |          double |                           :math:`h^{-1} {\rm Mpc}` | Angular-diameter distance based on  ``Z`` and a standard cosmology            |
|                            |                 |                                                    | (:math:`h=1; \Omega_M=0.3; \Omega_\Lambda=0.7`)                               |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|                  ``NSA_Z`` |          double |                                                    | Redshift from the NASA-Sloan Atlas                                            |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|              ``NSA_ZDIST`` |          double |                                                    | NSA distance estimate using pecular velocity model of Willick et al. (1997);  |
|                            |                 |                                                    | multiply by :math:`c/H_0` for Mpc.                                            |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|            ``LDIST_NSA_Z`` |          double |                           :math:`h^{-1} {\rm Mpc}` | Luminosity distance based on  ``NSA_Z`` and a standard cosmology              |
|                            |                 |                                                    | (:math:`h=1; \Omega_M=0.3; \Omega_\Lambda=0.7`)                               |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|            ``ADIST_NSA_Z`` |          double |                           :math:`h^{-1} {\rm Mpc}` | Angular-diameter distance based on  ``NSA_Z`` and a standard cosmology        |
|                            |                 |                                                    | (:math:`h=1; \Omega_M=0.3; \Omega_\Lambda=0.7`)                               |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|         ``NSA_ELPETRO_BA`` |          double |                                                    | NSA isophotal axial ratio from elliptical Petrosian analysis                  |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|        ``NSA_ELPETRO_PHI`` |          double |                                                deg | NSA isophotal position angle from elliptical Petrosian analysis               |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|     ``NSA_ELPETRO_TH50_R`` |          double |                                             arcsec | NSA elliptical Petrosian effective radius in the r-band; the is the same as   |
|                            |                 |                                                    | :math:`R_e` below.                                                            |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|          ``NSA_SERSIC_BA`` |          double |                                                    | NSA isophotal axial ratio from Sersic fit                                     |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|         ``NSA_SERSIC_PHI`` |          double |                                                deg | NSA isophotal position angle from Sersic fit                                  |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|        ``NSA_SERSIC_TH50`` |          double |                                             arcsec | NSA effective radius from the Sersic fit                                      |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|           ``NSA_SERSIC_N`` |          double |                                                    | NSA Sersic index                                                              |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
| **Version dependency and quality information**                                                                                                                                    |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|               ``VERSDRP2`` |             str |                                                    | Version of DRP used for 2d reductions                                         |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|               ``VERSDRP3`` |             str |                                                    | Version of DRP used for 3d reductions                                         |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|               ``VERSCORE`` |             str |                                                    | Version of mangacore used by the DAP                                          |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|               ``VERSUTIL`` |             str |                                                    | Version of idlutils used by the DAP                                           |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|                ``VERSDAP`` |             str |                                                    | Version of mangadap                                                           |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|               ``DRP3QUAL`` |             str |                                                    | DRP 3D quality bit (**link**)                                                 |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|                ``DAPQUAL`` |             str |                                                    | DAP quality bit (**link**)                                                    |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|  **DAP analysis flow information**                                                                                                                                                |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|               ``RDXQAKEY`` |             str |                                                    | Configuration keyword for the method used to assess the reduced data          |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|                 ``BINKEY`` |             str |                                                    | Configuration keyword for the spatial binning method                          |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|                  ``SCKEY`` |             str |                                                    | Configuration keyword for the method used to model the stellar-continuum      |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|                 ``ELMKEY`` |             str |                                                    | Configuration keyword that defines the emission-line moment measurement       |
|                            |                 |                                                    | method                                                                        |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|                 ``ELFKEY`` |             str |                                                    | Configuration keyword that defines the emission-line modeling method          |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|                  ``SIKEY`` |             str |                                                    | Configuration keyword that defines the spectral-index measurement method      |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|                ``BINTYPE`` |             str |                                                    | Type of binning used                                                          |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|                 ``BINSNR`` |             int |                                                    | Target for bin S/N, if Voronoi binning                                        |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|                 ``TPLKEY`` |             str |                                                    | The identifier of the template library, e.g., ``MILES``.                      |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|  **Additional info pulled from DAP fits file headers**                                                                                                                            |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|                ``DATEDAP`` |             str |                                                    | Date the DAP file was created and/or last modified.                           |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|                ``DAPBINS`` |             int |                                                    | The number of "binned" spectra analyzed by the DAP.                           |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|  **Data assessments provided specifically in the DAPall file**                                                                                                                    |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|                 ``RCOV90`` |          double |                                             arcsec | Semi-major axis radius (:math:`R`) below which spaxels cover at least 90% of  |
|                            |                 |                                                    | elliptical annuli with width :math:`R\pm 2.5` arcsec.  This should be         |
|                            |                 |                                                    | independent of the ``DAPTYPE``.                                               |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|                ``SNR_MED`` | double (vector) |                                                    | Median S/N per pixel in the ''griz'' bands within 1.0-1.5 :math:`R_e`.  This  |
|                            |                 |                                                    | should be independent of the ``DAPTYPE``.                                     |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|               ``SNR_RING`` | double (vector) |                                                    | S/N in the ''griz'' bands when binning all spaxels within 1.0-1.5             |
|                            |                 |                                                    | :math:`R_e`.  This should be independent of the ``DAPTYPE``.                  |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|                 ``SB_1RE`` |          double | :math:`10^{-17} {\rm erg/s/cm}^2{\rm /\AA/spaxel}` | Mean g-band surface brightness of valid spaxels within 1 :math:`R_e`.  This   |
|                            |                 |                                                    | should be independent of the ``DAPTYPE``.                                     |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|               ``BIN_RMAX`` |          double |                                        :math:`R_e` | Maximum g-band luminosity-weighted semi-major radius of any "valid" binned    |
|                            |                 |                                                    | spectrum.                                                                     |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|                ``BIN_R_N`` | double (vector) |                                                    | Number of binned spectra with g-band luminosity-weighted centers within 0-1,  |
|                            |                 |                                                    | 0.5-1.5, and 1.5-2.5 :math:`R_e`.                                             |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|              ``BIN_R_SNR`` | double (vector) |                                                    | Median g-band S/N of all binned spectra with luminosity-weighted centers      |
|                            |                 |                                                    | within 0-1, 0.5-1.5, and 1.5-2.5 :math:`R_e`.                                 |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|              ``STELLAR_Z`` |          double |                                                    | Flux-weighted mean redshift of the stellar component within a 2.5 arcsec      |
|                            |                 |                                                    | aperture at the galaxy center.                                                |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|         ``STELLAR_VEL_LO`` |          double |                                               km/s | Stellar velocity at 2.5% growth of all valid spaxels.                         |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|         ``STELLAR_VEL_HI`` |          double |                                               km/s | Stellar velocity at 97.5% growth of all valid spaxels.                        |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|    ``STELLAR_VEL_LO_CLIP`` |          double |                                               km/s | Stellar velocity at 2.5% growth after iteratively clipping :math:`3\sigma`    |
|                            |                 |                                                    | outliers.                                                                     |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|    ``STELLAR_VEL_HI_CLIP`` |          double |                                               km/s | Stellar velocity at 97.5% growth after iteratively clipping :math:`3\sigma`   |
|                            |                 |                                                    | outliers.                                                                     |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|      ``STELLAR_SIGMA_1RE`` |          double |                                               km/s | Flux-weighted mean stellar velocity dispersion of all spaxels within 1        |
|                            |                 |                                                    | :math:`R_e`.                                                                  |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
| ``STELLAR_CONT_RCHI2_1RE`` |          double |                                                    | Median reduced :math:`chi^2` of the stellar-continuum fit within 1            |
|                            |                 |                                                    | :math:`R_e`.                                                                  |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|                   ``HA_Z`` |          double |                                                    | Flux-weighted mean redshift of the Hα line within a 2.5 arcsec aperture at    |
|                            |                 |                                                    | the galaxy center.                                                            |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|             ``HA_GVEL_LO`` |          double |                                               km/s | Gaussian-fitted velocity of the H:math:`\alpha` line at 2.5% growth of all    |
|                            |                 |                                                    | valid spaxels.                                                                |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|             ``HA_GVEL_HI`` |          double |                                               km/s | Gaussian-fitted velocity of the H:math:`\alpha` line at 97.5% growth of all   |
|                            |                 |                                                    | valid spaxels.                                                                |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|        ``HA_GVEL_LO_CLIP`` |          double |                                               km/s | Gaussian-fitted velocity of the H:math:`\alpha` line at 2.5% growth after     |
|                            |                 |                                                    | iteratively clipping :math:`3\sigma` outliers.                                |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|        ``HA_GVEL_HI_CLIP`` |          double |                                               km/s | Gaussian-fitted velocity of the H:math:`\alpha` line at 97.5% growth after    |
|                            |                 |                                                    | iteratively clipping :math:`3\sigma` outliers.                                |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|          ``HA_GSIGMA_1RE`` |          double |                                               km/s | Flux-weighted H:math:`\alpha` velocity dispersion (from Gaussian fit) of all  |
|                            |                 |                                                    | spaxels within 1 :math:`R_e`.                                                 |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|           ``HA_GSIGMA_HI`` |          double |                                               km/s | H:math:`\alpha` velocity dispersion (from Gaussian fit) at 97.5% growth of    |
|                            |                 |                                                    | all valid spaxels.                                                            |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|      ``HA_GSIGMA_HI_CLIP`` |          double |                                               km/s | H:math:`\alpha` velocity dispersion (from Gaussian fit) at 97.5% growth after |
|                            |                 |                                                    | iteratively clipping :math:`3\sigma` outliers.                                |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|       ``EMLINE_RCHI2_1RE`` |          double |                                                    | Median reduced :math:`\chi^2` of the continuum+emission-line fit within 1     |
|                            |                 |                                                    | :math:`R_e`.                                                                  |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|       ``EMLINE_SFLUX_CEN`` | double (vector) |                  :math:`10^{-17} {\rm erg/s/cm}^2` | Summed emission-line flux integrated within a 2.5 arcsec aperture at the      |
|                            |                 |                                                    | galaxy center.                                                                |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|       ``EMLINE_SFLUX_1RE`` | double (vector) |                  :math:`10^{-17} {\rm erg/s/cm}^2` | Summed emission-line flux integrated within a 1-:math:`R_e` aperture at the   |
|                            |                 |                                                    | galaxy.                                                                       |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|       ``EMLINE_SFLUX_TOT`` | double (vector) |                  :math:`10^{-17} {\rm erg/s/cm}^2` | Total integrated flux of each summed emission measurement within the full     |
|                            |                 |                                                    | MaNGA field-of-view.                                                          |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|         ``EMLINE_SSB_1RE`` | double (vector) |     :math:`10^{-17} {\rm erg/s/cm}^2{\rm /spaxel}` | Mean emission-line surface-brightness from the summed flux measurements       |
|                            |                 |                                                    | within 1 :math:`R_e`.                                                         |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|        ``EMLINE_SSB_PEAK`` | double (vector) |     :math:`10^{-17} {\rm erg/s/cm}^2{\rm /spaxel}` | Peak summed-flux emission-line surface brightness.                            |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|         ``EMLINE_SEW_1RE`` | double (vector) |                                                ang | Mean emission-line equivalent width from the summed flux measurements within  |
|                            |                 |                                                    | 1 :math:`R_e`.                                                                |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|        ``EMLINE_SEW_PEAK`` | double (vector) |                                                ang | Peak emission-line equivalent width from the summed flux measurements.        |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|       ``EMLINE_GFLUX_CEN`` | double (vector) |                  :math:`10^{-17} {\rm erg/s/cm}^2` | Gaussian-fitted emission-line flux integrated within a 2.5 arcsec aperture    |
|                            |                 |                                                    | at the galaxy center.                                                         |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|       ``EMLINE_GFLUX_1RE`` | double (vector) |                  :math:`10^{-17} {\rm erg/s/cm}^2` | Gaussian-fitted emission-line flux integrated within a 1-:math:`R_e` aperture |
|                            |                 |                                                    | at the galaxy.                                                                |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|       ``EMLINE_GFLUX_TOT`` | double (vector) |                  :math:`10^{-17} {\rm erg/s/cm}^2` | Total integrated flux of the Gaussian fit to each emission line within the    |
|                            |                 |                                                    | full MaNGA field-of-view.                                                     |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|         ``EMLINE_GSB_1RE`` | double (vector) |     :math:`10^{-17} {\rm erg/s/cm}^2{\rm /spaxel}` | Mean emission-line surface-brightness from the Gaussian-fitted flux           |
|                            |                 |                                                    | measurements within 1 :math:`R_e`.                                            |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|        ``EMLINE_GSB_PEAK`` | double (vector) |     :math:`10^{-17} {\rm erg/s/cm}^2{\rm /spaxel}` | Peak Gaussian-fitted emission-line surface brightness.                        |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|         ``EMLINE_GEW_1RE`` | double (vector) |                                                ang | Mean emission-line equivalent width from the Gaussian-fitted flux             |
|                            |                 |                                                    | measurements within 1 :math:`R_e`.                                            |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|        ``EMLINE_GEW_PEAK`` | double (vector) |                                                ang | Peak emission-line equivalent width from the Gaussian-fitted flux             |
|                            |                 |                                                    | measurements.                                                                 |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|           ``SPECINDEX_LO`` | double (vector) |                                            ang,mag | Spectral index at 2.5% growth of all valid spaxels.                           |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|           ``SPECINDEX_HI`` | double (vector) |                                            ang,mag | Spectral index at 97.5% growth of all valid spaxels.                          |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|      ``SPECINDEX_LO_CLIP`` | double (vector) |                                            ang,mag | Spectral index at 2.5% growth after iteratively clipping :math:`3\sigma`      |
|                            |                 |                                                    | outliers.                                                                     |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|      ``SPECINDEX_HI_CLIP`` | double (vector) |                                            ang,mag | Spectral index at 97.5% growth after iteratively clipping :math:`3\sigma`     |
|                            |                 |                                                    | outliers.                                                                     |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|          ``SPECINDEX_1RE`` | double (vector) |                                            ang,mag | Median spectral index within 1 :math:`R_e`.                                   |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|                ``SFR_1RE`` |          double |                    :math:`h^{-2} {\rm M}_\odot/yr` | Simple estimate of the star-formation rate within 1 :math:`R_e` based on the  |
|                            |                 |                                                    | Gaussian-fitted H:math:`\alpha` flux;                                         |
|                            |                 |                                                    | :math:`\log {\rm SFR} = \log L_{{\rm H}\alpha} - 41.27` (Kennicutt & Evans    |
|                            |                 |                                                    | [2012, ARAA, 50, 531], citing Murphy et al. [2011, ApJ, 737, 67] and Hao et   |
|                            |                 |                                                    | al. [2011, ApJ, 741, 124]; Kroupa IMF), where :math:`L_{{\rm H}\alpha}`       |
|                            |                 |                                                    | = 4:math:`\pi` EML_FLUX_1RE (LDIST_Z):math:`^2` and *no* attenuation          |
|                            |                 |                                                    | correction has been applied.                                                  |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+
|                ``SFR_TOT`` |          double |                    :math:`h^{-2} {\rm M}_\odot/yr` | Simple estimate of the star-formation rate within the IFU field-of-view based |
|                            |                 |                                                    | on the Gaussian-fitted H:math:`\alpha` flux;                                  |
|                            |                 |                                                    | :math:`\log {\rm SFR} = \log L_{{\rm H}\alpha} - 41.27` (Kennicutt & Evans    |
|                            |                 |                                                    | [2012, ARAA, 50, 531], citing Murphy et al. [2011, ApJ, 737, 67] and Hao et   |
|                            |                 |                                                    | al. [2011, ApJ, 741, 124]; Kroupa IMF), where :math:`L_{{\rm H}\alpha}`       |
|                            |                 |                                                    | = 4:math:`\pi` EML_FLUX_1RE (LDIST_Z):math:`^2` and *no* attenuation          |
|                            |                 |                                                    | correction has been applied.                                                  |
+----------------------------+-----------------+----------------------------------------------------+-------------------------------------------------------------------------------+

.. note::

 * Distance estimates do not include an estimate of the peculiar
   motions.
 * Volume weights are included in the DRPall file.
 * ``RCOV90`` is calculated for the ``CUBE`` files; however, a more
   sophisticated calculation would use the ``RSS`` files to account
   for the significant overlap of the fiber "beams."
 * All radially averaged or summed properties are calculated within
   ''elliptical'' apertures defined using the NSA ellipticity and
   position angle.
 * Possible Future developments: (1) Provide default set of cross
   matching: SDSS I/II, Galaxy Zoo? (2) Include initial radial profiles
   of the emission-line, spectral-index, and other derived properties?

