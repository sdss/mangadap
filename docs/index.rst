The MaNGA Data Analysis Pipeline
================================

.. image:: https://github.com/sdss/mangadap/actions/workflows/ci_tests.yml/badge.svg
   :target: https://github.com/sdss/mangadap/actions
   :alt: Build Status

.. image:: https://codecov.io/gh/sdss/mangadap/branch/master/graph/badge.svg?token=S4veEPJwS1
   :target: https://codecov.io/gh/sdss/mangadap
   :alt: Coverage Status

.. image:: https://readthedocs.org/projects/sdss-mangadap/badge/?version=latest
  :target: https://sdss-mangadap.readthedocs.io/en/latest/
  :alt: Documentation Status

.. image:: https://img.shields.io/github/license/sdss/mangadap
   :target: https://github.com/sdss/mangadap/blob/master/LICENSE.md
   :alt: License

.. image:: http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)
   :target: http://www.astropy.org/
   :alt: Powered by Astropy

The MaNGA data-analysis pipeline (MaNGA DAP) is the survey-led software package
that has analyzed all galaxy data produced by the MaNGA data-reduction pipeline
(MaNGA DRP).  Its goal is to produce high-level, science-ready data products
derived from MaNGA spectra.  The products currently provided are:

 * Spatially stacked spectra
 * Stellar kinematics
 * Nebular emission-line properties: fluxes, equivalent widths, and
   kinematics
 * Spectral indices: absorption-line (e.g., H-delta) and bandhead/color
   (e.g., TiO, D4000) measurements

That is, the DAP currently focuses on "model-independent" properties.
Higher-level, model-dependent properties, such as stellar-population parameters,
are outside of the scope of the current pipeline.

.. note::
    
    This documentation is for the most recent stable release of the DAP.  Use
    the readthedocs selection menu at the lower-left of your web browser to
    select previous version relevant to specific SDSS data releases:

        * `SDSS-IV/MaNGA DR15 <https://www.sdss.org/dr15/manga/>`__ is based on
          version 2.2.1

        * `SDSS-IV/MaNGA DR17 <https://www.sdss.org/dr17/manga/>`__ is based on
          version 3.1.2.  Note that the version of the code that produced the
          DR17 was version 3.1.0; version 3.1.1 and 3.1.2 only include
          documentation updates.

.. _citation:

Citation
========

If you use the DAP software and/or its output products, please cite the following two papers:

 - *Overview*: `Westfall et al. (2019, AJ, 158, 231)`_
 - *Emission-line Modeling*: `Belfiore et al. (2019, AJ, 158, 160)`_

Additionally, if you use SDSS-IV/MaNGA data, please see:

 * `How to Cite SDSS <https://www.sdss.org/collaboration/citing-sdss/>`__
 * `SDSS Technical Publications <https://www.sdss.org/science/technical_publications/>`__

----

.. toctree::
   :caption: Welcome!
   :maxdepth: 2

   whatsnew
   knownissues

**New in version 4.x**: We have made the MaNGA DAP *much* easier to use with
non-MaNGA datacubes.  These changes are reflected throughout our documentation,
but specifically see :ref:`fitdatacube`.

----

.. toctree::
   :caption: Usage
   :maxdepth: 1

   installation
   execution
   plan
   development

----

.. toctree::
   :caption: Methods
   :maxdepth: 1

   workflow
   templatelibraries
   emissionlines
   spectralindices
   spatialcovariance

----

.. toctree::
   :caption: Output
   :maxdepth: 1

   gettingstarted
   metadatamodel
   datamodel
   qualityassessment
   corrections
   aggregation

----

.. toctree::
   :caption: Beyond MaNGA
   :maxdepth: 1

   fitonespec
   fitdatacube

The MaNGA DAP software was written to be more broadly useful. The
pages above provide some examples of how the code can be used more
generally, beyond its application to MaNGA.

----

.. toctree::
   :caption: Utilities
   :maxdepth: 1

   covariance
   resampling
   resolution
   smoothing
   bitmasks
   parameters

----

.. toctree::
   :caption: API
   :maxdepth: 2

   MaNGA DAP API <api/mangadap>
   MaNGA DAP Modules <api/modules>

----

Search
======

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


.. include links

.. include:: include/links.rst
