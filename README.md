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

To install, export the code from git (although note that this **will not include the documentation submodules**; see [Cloning the full repo](#cloning-the-full-repo) below):

`git clone https://github.com/sdss/mangadap.git`

Then run the `setup.py` file to install it.

`python3 setup.py install`

Or

`python3 setup.py develop`

for development purposes.

You can also test the installation using

`python3 setup.py test`

Although the tests are severely lacking so far.

## Usage

## Citation

## Cloning the full repo

This repo has submodules linked to Overleaf; see docs/README.md.  Of
course, this documentation is not needed to run the code.  However, to
pull across all the submodule files as well as the main repo files, use
the `--recursive` option:

`git clone --recursive https://github.com/sdss/mangadap.git`



