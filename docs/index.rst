The MaNGA Data Analysis Pipeline
================================

.. image:: http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)
   :target: http://www.astropy.org/
   :alt: Powered by Astropy

.. image:: https://travis-ci.org/sdss/mangadap.svg?branch=master
   :target: https://travis-ci.org/sdss/mangadap
   :alt: Build Status

.. image:: https://coveralls.io/repos/github/sdss/mangadap/badge.svg?branch=master
   :target: https://coveralls.io/github/sdss/mangadap?branch=master
   :alt: Coverage Status

.. image:: https://readthedocs.org/projects/sdss-mangadap/badge/?version=latest
  :target: https://sdss-mangadap.readthedocs.io/en/latest/
  :alt: Documentation Status

.. image:: https://img.shields.io/github/license/sdss/mangadap
   :target: https://github.com/sdss/mangadap/blob/master/LICENSE.md
   :alt: License

The MaNGA data-analysis pipeline (MaNGA DAP) is the survey-led software
package that analyzes the data produced by the MaNGA data-reduction
pipeline (MaNGA DRP) to produced physical properties derived from the
MaNGA spectroscopy.

These quantities currently include:

 * Spatially stacked spectra
 * Stellar kinematics
 * Nebular emission-line properties: fluxes, equivalent widths, and
   kinematics
 * Spectral indices: absorption-line (e.g., H-delta) and bandhead/color
   (e.g., TiO, D4000) measurements

All survey-led properties are derived from the datacubes,
specifically the LOGCUBE files. However, the core functions are
developed to consider each spectrum largely independently. The DAP
currently focuses on "model-independent" properties. Higher-level,
model-dependent properties, such as stellar-population parameters,
are outside of the scope of the current pipeline.

.. warning::
    
    This documentation is for the most recent tag of the DAP.  If you're
    looking for DR15 (DAP version 2.2.1) documentation, see the `DR15
    MaNGA Overview <https://www.sdss.org/dr15/manga/>`_.

.. _citation:

Citation
========

If you use the DAP software and/or its output products, please cite the following two papers:

 - *Overview*: `Westfall et al. (2019, AJ, 158, 231)`_
 - *Emission-line Modeling*: `Belfiore et al. (2019, AJ, 158, 160)`_

----

.. toctree::
   :caption: Welcome!
   :maxdepth: 2

   whatsnew
   knownissues

----

.. toctree::
   :caption: Usage
   :maxdepth: 1

   installation
   execution
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
   fitonecube

The MaNGA DAP software was written to be more broadly useful. The
pages above provide some examples of how the code can be used more
generally, beyond its application to MaNGA.

----

.. toctree::
   :caption: Utilities
   :maxdepth: 1

   parameters
   bitmasks
   datacube
   covariance
   resampling
   resolution
   smoothing

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
