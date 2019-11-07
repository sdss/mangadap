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

.. toctree::
   :maxdepth: 2

   whatsnew
   knownissues

How Tos
=======

.. toctree::
   :maxdepth: 1

   installation
   execution
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

