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

To install, first export the code from git:

`git clone https://github.com/sdss/mangadap.git`

This **will not include the documentation submodules**; see [Cloning the
full repo](#cloning-the-full-repo) below)

----

There are a few ways to install (always from within the top-level directory of the repo):

 - To install the DAP, run: `python3 setup.py install`

   On MacOSX, you may need to add `CC=clang`, i.e.: `CC=clang python3 setup.py install`

 - To install the DAP in a way that makes it easier to develop, run: `python3 setup.py develop`

 - To install the DAP and update its dependencies as necessary, run:
   `pip3 install -e .`

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

**Notes:**
 - `$MANGACORE_VER` and `$MANGACORE_DIR` are only needed to perform the
   survey-level execution of the pipeline.
 - The DAP expects to find the DRP `LOGCUBE` files in
   `$MANGA_SPECTRO_REDUX/$MANGADRP_VER/[plate]/stack`, where `[plate]`
   is the desired plate number.
 - The DAP expects to find/write data to
   `$MANGA_SPECTRO_ANALYSIS/$MANGADRP_VER/$MANGADAP_VER`.
 - `$MANGADAP_VER` is only used to set the path names, not to select the
   specific version of the pipeline to use

These environmental variables can be added to your `.bash_profile` file in your home directory or be included in a script that is sourced when you want to run the DAP.  The added lines to your `.bash_profile` file could look something like this:

```
export MANGA_SPECTRO_REDUX=/Volumes/MaNGA/redux
export MANGA_SPECTRO_ANALYSIS=/Volumes/MaNGA/analysis

export MANGADRP_VER=v2_4_3

export MANGADAP_VER=2.2.1
export MANGADAP_DIR=$HOME/Work/MaNGA/dap/repo/mangadap

export MANGACORE_VER=trunk
export MANGACORE_DIR=$HOME/Work/MaNGA/core/$MANGACORE_VER
```

**Note**: Some of these variables are also defined by Marvin; see
(here)[https://sdss-marvin.readthedocs.io/en/stable/installation.html].
It's possible to have both Marvin and the DAP point to the same
directory, but beware that this may mean that some of the files get
overwritten.

## Usage

## Citation

## Cloning the full repo

This repo has submodules linked to Overleaf; see docs/README.md.  Of
course, this documentation is not needed to run the code.  However, to
pull across all the submodule files as well as the main repo files, use
the `--recursive` option:

`git clone --recursive https://github.com/sdss/mangadap.git`



