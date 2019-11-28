# The MaNGA Data Analysis Pipeline

[![Build Status](https://travis-ci.org/sdss/mangadap.svg?branch=master)](https://travis-ci.org/sdss/mangadap)
[![Coverage Status](https://coveralls.io/repos/github/sdss/mangadap/badge.svg?branch=master)](https://coveralls.io/github/sdss/mangadap?branch=master)
[![Doc Status](https://readthedocs.org/projects/sdss-mangadap/badge/?version=latest)](https://sdss-mangadap.readthedocs.io/en/latest/)
[![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)
[![License](https://img.shields.io/github/license/sdss/mangadap)](https://github.com/sdss/mangadap/blob/master/LICENSE.md)

The MaNGA data-analysis pipeline (MaNGA DAP) is the survey-led software
package that analyzes the data produced by the MaNGA data-reduction
pipeline (MaNGA DRP) to produced physical properties derived from the
MaNGA spectroscopy.

For full documentation, see: https://sdss-mangadap.readthedocs.io/en/latest/

**Note that the version of the DAP used for DR15 is [2.2.1](https://github.com/sdss/mangadap/releases/tag/2.2.1).**

## Citation

If you use the DAP software and/or its output products, please cite the following two papers:

 - *Overview*: [Westfall et al. (2019, AJ, 158, 231)](https://ui.adsabs.harvard.edu/abs/2019AJ....158..231W/abstract)
 - *Emission-line Modeling*: [Belfiore et al. (2019, AJ, 158, 160)](https://ui.adsabs.harvard.edu/abs/2019AJ....158..160B/abstract)

## Installation

To install, first clone the repository:

`git clone https://github.com/sdss/mangadap.git`

We recommend using the most recent tag:

```
cd mangadap
./checkout_current_tag
```

----

To install, run:

`python3 setup.py install`

On MacOSX, you may need to add `CC=clang`, i.e.:
   
`CC=clang python3 setup.py install`

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

These environmental variables can be added to your `.bash_profile` file
in your home directory or be included in a script that is sourced when
you want to run the DAP.  The added lines to your `.bash_profile` file
could look something like this:

```
export MANGA_SPECTRO_REDUX=/Volumes/MaNGA/redux
export MANGA_SPECTRO_ANALYSIS=/Volumes/MaNGA/analysis

export MANGADRP_VER=v2_4_3

export MANGADAP_VER=2.2.1

export MANGACORE_VER=v1_6_2
export MANGACORE_DIR=$HOME/MaNGA/core/$MANGACORE_VER
```

**Note**: Some of these variables are also defined by Marvin; see
[here](https://sdss-marvin.readthedocs.io/en/stable/installation.html).
It's possible to have both Marvin and the DAP point to the same
directory, but beware that this may mean that some of the files get
overwritten.




