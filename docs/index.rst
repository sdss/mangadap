Welcome to the MaNGA Data Analysis Pipeline
===========================================

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

All survey-led properties are derived from the datacubes, specifically
the LOGCUBE files.  However, the core functions are developed to
consider each spectrum largely independently.  The DAP currently focuses
on "model-independent" properties.  Higher-level, model-dependent
properties, such as stellar-population parameters, will be included in
future releases on a best-effort basis.

.. _citation:

Citation
========

If you use the DAP software and/or its output products, please cite the following two papers:

 - *Overview*: `Westfall et al. (2019)
   <https://ui.adsabs.harvard.edu/abs/2019AJ....158..231W/abstract>`_
 - *Emission-line Modeling*: `Belfiore et al. (2019)
   <https://ui.adsabs.harvard.edu/abs/2019AJ....158..160B/abstract>`_

----

.. toctree::
   :maxdepth: 2

   whatsnew
   knownissues

Usage
=====

.. toctree::
   :maxdepth: 1

   installation
   execution
   development

Output
======

.. toctree::
   :maxdepth: 1

   gettingstarted
   metadatamodel
   datamodel
   corrections
   qualityassessment

Beyond MaNGA
============

The MaNGA DAP software was written to be more broadely useful.  What
follows are some examples of how the code can be used more generally
beyond its application to MaNGA (see also generality discussion included
in our :ref:`development`).

.. toctree::
   :maxdepth: 1

   fit_one_spec

API Reference
=============

.. toctree::
   :maxdepth: 2

   MaNGA DAP API <api/mangadap>

Search
======

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

