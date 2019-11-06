Welcome to the MaNGA Data Analysis Pipeline
===========================================

The MaNGA data-analysis pipeline (MaNGA DAP) is the survey-led software
package that analyzes the data produced by the MaNGA data-reduction
pipeline (MaNGA DRP) to produced physical properties derived from the
MaNGA spectroscopy.

These quantities currently include:

 - Spatially stacked spectra
 - Stellar kinematics (V and sigma)
 - Nebular emission-line properties: fluxes, equivalent widths, and kinematics (V and sigma)
 - Spectral Indices: absorption-line (e.g., H-delta) and bandhead (e.g., D4000) measurements

All survey-led properties are derived from the datacubes, specifically
the LOGCUBE files. However, the core functions are developed to consider
each spectrum largely independently. The DAP currently focuses on
"model-independent" properties. Higher-level, model-dependent
properties, such as stellar-population parameters, will be included in
future releases on a best-effort basis.

Getting Started
===============

.. toctree::
   :maxdepth: 2

   intro

How Tos
=======

.. toctree::
   :maxdepth: 2

   fit_one_spec

API Reference
=============

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   MaNGA DAP API <api/mangadap>


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

