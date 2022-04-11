# The MaNGA Data Analysis Pipeline

[![Build Status](https://github.com/sdss/mangadap/actions/workflows/ci_tests.yml/badge.svg)](https://github.com/sdss/mangadap/actions)
[![Coverage Status](https://codecov.io/gh/sdss/mangadap/branch/master/graph/badge.svg?token=S4veEPJwS1)](https://codecov.io/gh/sdss/mangadap)
[![Doc Status](https://readthedocs.org/projects/sdss-mangadap/badge/?version=latest)](https://sdss-mangadap.readthedocs.io/en/latest/)
[![License](https://img.shields.io/github/license/sdss/mangadap)](https://github.com/sdss/mangadap/blob/master/LICENSE.md)
[![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)

The MaNGA data-analysis pipeline (MaNGA DAP) is the survey-led software package
that has analyzed all galaxy data produced by the MaNGA data-reduction pipeline
(MaNGA DRP).  Its goal is to produce high-level, science-ready data products
derived from MaNGA spectra.  The products currently provided are:

 - Spatially stacked spectra
 - Stellar kinematics
 - Nebular emission-line properties: fluxes, equivalent widths, and
   kinematics
 - Spectral indices: absorption-line (e.g., H-delta) and bandhead/color
   (e.g., TiO, D4000) measurements

That is, the DAP currently focuses on "model-independent" properties.
Higher-level, model-dependent properties, such as stellar-population parameters,
are outside of the scope of the current pipeline.

**See our full documentation at [sdss-mangadap.readthedocs.io](https://sdss-mangadap.readthedocs.io/en/latest/).**


## SDSS Data Release Versions

 - [SDSS-IV/MaNGA DR15](https://www.sdss.org/dr15/manga/) is based on version 2.2.1

 - [SDSS-IV/MaNGA DR17](https://www.sdss.org/dr17/manga/) is based on version
   3.1.2.  Note that the version of the code that produced the DR17 was version
   3.1.0; version 3.1.1 and 3.1.2 only include documentation updates.


## Citation

If you use the DAP software and/or its output products, please cite the
following two papers:

 - *Overview*: [Westfall et al. (2019, AJ, 158, 231)](https://ui.adsabs.harvard.edu/abs/2019AJ....158..231W/abstract)
 - *Emission-line Modeling*: [Belfiore et al. (2019, AJ, 158, 160)](https://ui.adsabs.harvard.edu/abs/2019AJ....158..160B/abstract)

Additionally, if you use SDSS-IV/MaNGA data, please see:

 - [How to Cite SDSS](https://www.sdss.org/collaboration/citing-sdss/)
 - [SDSS Technical Publications](https://www.sdss.org/science/technical_publications/)


## Installation

`pip install sdss-mangadap`

See https://sdss-mangadap.readthedocs.io/en/latest/installation.html


