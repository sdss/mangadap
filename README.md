# The MaNGA Data Analysis Pipeline

The MaNGA data-analysis pipeline (MaNGA DAP) is the survey-led software
package that analyzes the data produced by the MaNGA data-reduction
pipeline (MaNGA DRP) to produced physical properties derived from the
MaNGA spectroscopy.

For DR15 (version 2.2.2), these quantities are:

 - Spatially stacked spectra
 - Stellar kinematics (V and sigma)
 - Nebular emission-line properties: fluxes, equivalent widths, and
   kinematics (V and sigma)
 - Spectral Indices: absorption-line (e.g., H-delta) and bandhead (e.g.,
   D4000) measurements 

All survey-led properties are derived from the datacubes, specifically
the LOGCUBE files. However, the core functions are developed to consider
each spectrum largely independently. The DAP currently focuses on
"model-independent" properties. Higher-level, model-dependent
properties, such as stellar-population parameters, will be included in
future releases on a best-effort basis. 

## Installation

To install, export the code from git:

`git clone https://github.com/sdss/mangadap.git`

This **will not include the documentation submodules**; see [Cloning the
full repo](#cloning-the-full-repo) below)

----

Install the DAP by running the setup script (for development use `develop` instead of `install`):

`python3 setup.py install`

----

The tests are rather lacking still, but you can test the installation using

`python3 setup.py test`

----

The DAP uses environmental variable to define the paths to specific data
and other repositories.  If these are not defined, warnings will be
issued everytime the DAP is installed or imported.  The relevant
environmental variables, their default, and their usage are provided
below.

|                 Variable |                           Default |                                       Comments |
|:------------------------ |:--------------------------------- |:---------------------------------------------- |
| `MANGADRP_VER`           | `v2_4_3`                          | Version of the DRP, used for path construction |
| `MANGA_SPECTRO_REDUX`    | `$HOME/MaNGA/redux`               | Root path for the reduced data                 |
| `MANGADAP_VER`           | `mangadap.__version__`            | Version of the DAP, used for path construction |
| `MANGA_SPECTRO_ANALYSIS` | `$HOME/MaNGA/analysis`            | Root path for the analysis data                |
| `MANGACORE_VER`          | `v1_6_2`                          | Version of MaNGA core (survey-level meta data) |
| `MANGACORE_DIR`          | `$HOME/MaNGA/core/$MANGACORE_VER` | Root path with the MaNGA core repository       |

Notes:

    - `$MANGACORE_VER` and `$MANGACORE_DIR` are only needed to perform
      the survey-level execution of the pipeline.
    - The DAP expects to find the DRP `LOGCUBE` files in
      `$MANGA_SPECTRO_REDUX/$MANGADRP_VER/[plate]/stack`, where
      `[plate]` is the desired plate number.
    - The DAP expects to find/write data to
      `$MANGA_SPECTRO_ANALYSIS/$MANGADRP_VER/$MANGADAP_VER`.
    - `$MANGADAP_VER` is only used to set the path names, not to select
      the specific version of the pipeline to use

## Usage

## Citation

## Cloning the full repo

This repo has submodules linked to Overleaf; see docs/README.md.  Of
course, this documentation is not needed to run the code.  However, to
pull across all the submodule files as well as the main repo files, use
the `--recursive` option:

`git clone --recursive https://github.com/sdss/mangadap.git`



