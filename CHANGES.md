
3.1.3dev
--------

 - Update to custom datacube class implementation.  Previous version
   required a file within the mangadap/ directory structure; new code
   allows for imports from a local user file.
 - Deal with cubes where all the spaxels have valid data.
 - Updates to `manga_dap_inspector` to use a cursor on the map selector
   and a range slider for the vmin,vmax in the map.
 - Refactor installation details inspired by `pypeit` package.
    - Installation dependencies and setup moved into `setup.cfg`
    - Versioning now done automatically
    - Change scripting protocol
    - tox testing CI    
    - Updated MANIFEST


3.1.2 (29 Jul 2021)
-------------------

 - Additional updates to documentation for DR17 release.
 - Includes includes of MaNGA IDs used in the construction of the MaStar
   HC version 2 library.


3.1.1 (1 Mar 2021)
------------------

 - Updated documentation for MPL-11 release.


3.1.0 (23 Dec 2020)
-------------------

 - Added DataTable, mostly for automated documentation purposes.
 - Removed revision histories
 - Updated docs for MPL-10
 - Enabled fitting a restricted wavelength range in the emission-line
   module.
 - Add parameter to change the required fraction of spectral coverage to
   perform a fit.
 - Detailed handling of masks for binning, MW extinction correction, and
   emission-line continuum construction.
 - Enable the output grid for :class:`mangadap.util.sampling.Resample`
   to be defined directly.
 - Include covariance calculation in
   :class:`mangadap.util.sampling.Resample`.
 - Add a class that reads trace sets produced by IDLUTILS,
   :class:`mangadap.util.trace.TraceSet`, but not well tested yet.
 - Fix untested velocity registration and stacking code, and added
   registration test
 - Upgraded package requirements
 - Line profile in `EmissionLineTemplates` now hardwired to be
   `FFTGaussianLSF`.
 - Changed the format of the emission-line database files:
    - Elric will *really* no longer work
    - Format specifies tying in a 4 element vector for each line
    - Adds ability to specify inequality constraints on velocity and
      dispersion

3.0.1 (31 May 2020)
-------------------

 - Fixed bug in association of spaxels with bins in emission-line module
 - Catch failures in pPXF in emission-line module, spaxel/bin flagged as
   `FIT_FAILED`, execution keeps going.
 - 7443-3703 redshift fix
 - Added LSF extension with LSF of binned spectra to output model cube
   file
 - Fixed corner cases in `StellarContinuumModel.matched_kinematics` and
   `EmissionLineModel.matched_kinematics` where all bins are masked.
 - Fixed bug in use of `ColorbarBase` in some plotting scripts.

3.0.0 (27 May 2020)
-------------------

 - Directory restructuring:  Removed "python" top-level directory and
   moved "data" into "mangadap".  This should make package distribution
   via pip easier.
 - Remove "default" prefix from methods in `mangadap.config.defaults`.
 - Remove obsolete `dapsrc` keyword arguments
 - Some general clean-up of commented code and docs.
 - Import clean-up, including removal of any `from __future__` imports.
 - Added `bin/dap_status` to check the status of a batch DAP run.
 - `__credits__` now include all authors of the DAP papers
 - Added new DataCube and RowStackedSpectra classes, beginning the
   generalization of the interface with the input data for use with
   non-MaNGA data.
 - For instantiating a MaNGADataCube, changed from using a yanny par
   file to a configuration (ini) file.  Included code that can write
   these files using data from the DRPall or DRPComplete files.
 - Included a script that will download data into a new
   `mangadap/data/remote` directory for testing.  The directory is not
   included in the repo and has been added to .gitignore to prevent it
   from being accidentally added.
 - Included a number of tests that use the new remote data.  These will
   be skipped if the remote data is not available.
 - Significant improvements and testing of `mangadap.util.covariance`.
   Ensuring that the correlation matrices provided by the DRP can be
   read by `MaNGADataCube` and are effectively identical to the
   calculation performed by `MaNGARSS`.
 - Moved all core script code from `bin` to `mangadap/scripts`.  Code in
   `bin` now make simple calls to these scripts.  Moved `rundap.py` from
   `mangadap.survey` to `mangadap.scripts`.
 - General clean-up of the executables in `bin`.
 - Integrated use of `MaNGADataCube` instead of `DRPFits` when executing
   the pipeline.  Usage of `DRPFits` is now largely obsolete.  Not
   removed yet.  Interface is still not abstracted enough for a general
   `DataCube` (because of naming conventions), but getting there.
 - Tests now include a nominal run of `manga_dap` using the `ALL`
   binning scheme.
 - Usage of `ObsInputPar` is now obsolete.
 - Spectral resolution is no longer chosen by the spectral binning
   module; instead the spectral resolution is selected when reading the
   datacube (makes much more sense!).  Led to some clean-up of the
   binning config files.
 - To select different spectral resolution extensions, use command-line
   arguments in `rundap` or `write_dap_config`.
 - Started doc edits.
 - Example scripts now use `MaNGADataCube` instead of `DRPFits`, but
   further clean-up of the `examples` directory is needed.
 - Include a default analysis plan, so that executions of `manga_dap`
   don't require an analysis plan yanny parameter file.
 - Docstring updates for modules up through `StellarContinuumModel`, but
   still lots to do.
 - Consolidates the common functionality in `MaNGADataCube` and
   `MaNGARSS` into a new, repurposed `DRPFits` base class.
 - Global data from unique DAPTYPEs now put in separate extensions of
   the DAPall file.
 - Incorporate change from 'PREDISP' to 'LSFPRE' and 'DISP' to 'LSFPOST'
   in the spectral resolution extensions produced by the DRP.
 - Added DataCube metadata.
 - Added Legendre polynomial fitter
 - Added MaStar HC V2 library
 - Include additional spectral index measurements and weights to use
   when combining indices; leads to changes in MAPS data model.
 - Updated versions of required packages, including incrementing version
   dependency to ppxf==7.2.0; code still uses default `linear_method`
   pending testing of `linear_method='lsqlin'`.
 - Added photometry fix file; only adjusts 8083-12702 pending further
   assessments.
 - Added redshift fixes for IC342, Coma, and globular cluster cubes.
 - Fixed minor bug in coordinates used in HYB binning cases for some
   DAPall calculations.


2.5.2 (14 Feb 2020)
-------------------

 - Had to retag for pip release after increase in space limit to 150MB.
   Version 2.5.1 on pip (code-base only) removed.


2.5.1 (13 Feb 2020)
-------------------

 - Jumped the gun on tagging 2.5.0.  Fixed Travis testing and package
   data installation issues.


2.5.0 (13 Feb 2020)
-------------------

 - Added BPASS theoretical spectral template library
 - Added numerous config files / line lists to allow for masking of NaD
   in both stellar continuum modeling and emission line modeling steps
 - Adjusted nsigma parameter in `pixelmask.py` so that it can be set
   externally. (This is needed to avoid masking HeI 5877 emission while
   simultaneously masking NaD.)
 - Fixed a bug in `angstroms_per_pixel` that caused it to return a
   single number instead of a vector.
 - Significant doc improvements by migrating most of the information on
   the internal wiki into the docs/ directory.
 - Significant reworking of I/O interface for "spectral feature
   databases".  These were databases where the config files basically
   just pointed to the file in the data/ directory.  The new interface
   circumvents the need for these config files.
 - Added tests to ensure available configuration files are valid.

2.4.1 (01 Nov 2019)
-------------------
 - Bug fix in construction of `STELLAR` extension in model LOGCUBE file;
   construction used a pointer instead of a copy that cause the
   `STELLAR` and `MODEL` extensions to be identical.

2.4.0 (31 Oct 2019)
-------------------

 - Hotfix to accommodate change in padding computation in ppxf>=6.7.15
 - Fixed units bug in `flux_to_fnu`
 - Allow emission-line moment, modeling, and spectral-index keys to be
   None, meaning that the DAP only fits the stellar kinematics but still
   produces the MAPS and LOGCUBE files.
 - Fixed bug in templates when using iteration_mode = 'global_template'
 - Change from `time.clock()` to `time.perf_counter()`
 - Bug fix in record array dimensionality when writing to binary table
 - Minor plotting changes for Overview paper plots
 - Added a channel with R in h^-1 kpc to both the spaxel and bin
   elliptical-coordinate extensions
 - Remove MASKNAME header keyword inherited from DRP files
 - Added continuum measurements for both emission-line moments and
   spectral indices to MAPS files
 - Add model spectral indices to the MAPS files
 - Added new emission lines to parameter and bandpass database files
 - Added beta version of hierarchically clustered MaStar templates
 - Fixed masking of emission-line chi-square measurements
 - Adjusted fixed flux ratio of doublets based on improved calculation:
   `(M1_1+E2_1)/(M1_2+E2_2) * L_2/L_1`, where `M1` and `E2` are the
   magnetic dipole and electric quadrapole Einstein A coefficients and
   `L` is the line wavelength.  The Einstein A coefficients are all
   taken from NIST, except for the OII 7320,7330 complext (two doublets)
   which use the values from Zeippen (1987), A&A, 173, 410.  In all but
   this OII 7320,7330 complex, the magnetic dipole term dominates such
   that the electric quadrapole terms are effectively irrelevant.
 - Added a script that computes the weights used during the `ppxffit_qa`
   plot construction.

2.3.0 (31 Jan 2019)
-------------------

 - Put in default paths when environmental variables not defined and
   warn user, both when running setup.py and when importing mangadap
 - Start testing suite (still very limited)
 - Fix error in constructing parameter tying vector passed to ppxf
 - Allow change of templates used during stellar-continuum and
   emission-ine fit.
 - Fixed bug in post-fit chi-square calculation in PPXFFit
 - Updated requirements
 - Fixed bug in `passband_integral` error calculation, but MC tests show
   that it's a poor substitute for an MC simulation.
 - Significant changes to the resampling code; replaced `resample1d`
   with `Resample` class.
 - Added test data to `data/tests/`
 - Added tests for more core routines, like PPXFFit.fit() and
   Sasuke.fit().  Tests take about 1 min to execute.
 - Fix FORESTAR bug in spectral indices.
 - Fix some bugs in the bandpass integral functions introduced during
   recent testing.
 - Add spectral indices test
 - Change model-cube data model to use `STCNT` and `STCNT_MASK`
   extensions and remove the `EMLINE_BASE` and `EMLINE_MASK` extensions.
 - Allow the emission-line fitter to use the bin ID numbers directly
   instead of matching the spaxels to bins by coordinate proximity
 - Construction of the parameter tying object in the emission-line
   fitter is now done just before each spectrum is fit by ppxf (not
   globally in Sasuke) to better handle when components are omitted.
 - When deconstructing bins into spaxels, the second fit iteration
   *only* fits spaxels that are components of binned spectra and does
   not refit spectra that constitute an entire bin themselves.  I.e.
   this removes some largely redundant fitting.
 - PPXFFit and Sasuke include calculations of the chi-square growth; and
   changed the names of the growth columns in the reference files.
 - Changed definitions of A to be the model amplitude; A/N is the model
   amplitude divided by the median noise in the two sidebands defined
   for the emission-line EW calculation.
 - Consolidated figures-of-merit for the stellar-continuum fitting into
   a single MAPS extension and added additional quantities.
 - Added MAPS with the amplitude, A/N, and local reduced chi-square for
   each line fit in the emission-line modeling.
 - Added MAPS with the full-spectrum figures-of-merit for the
   emission-line module, similar to the stellar-continuum output from
   ppxf.
 - Major changes to survey-level code: drpcomplete.py and rundap.py.
 - Fixed minor issue in EmissionLineMoments.single_band_moments that
   caused many "UserWarning: Warning: converting a masked element to
   nan." messages.
 - Added MODEL_MASK to model cube file to differentiate between model
   mask and binned spectra mask
 - Added fit QA plots for stellar and emission line modules for each
   PLTIFU, and a QA plot combining all observations on a given plate.
 - Added calls to additional QA plots in rundap.py
 - Changed daptype to binning-stellartemplates-emissionlinetemplates
   using default_dap_method.  StellarContinuumModel reference file moved
   to the common directory.
 - Add databases and configuration files that include the HeI 3889 line
   and the NI ~5200 doublet.
 - Added post-processing qa plots
 - Fixed minor issue in propagating masks from the reference files to
   the maps files; primarily an issue for the hybrid binning scheme.
 - Ignore minimum S/N limitation in emission-line moments and spectral
   indices for hybrid scheme as a stop-gap to minimize differences in
   moments, models, and indices BINID.

2.2.2 (Not released)
--------------------
 - Code migrated from SVN to GitHub

2.2.1 (28 Mar 2018): MPL-7/DR15 Candidate tag
---------------------------------------------
 - Add DAPFRMT and DAPTYPE header keywords to MAPS and LOGCUBE files.

2.2   (23 Mar 2018): Initial tag for DR15
-----------------------------------------
 - Fixed NaN problems in DAPall file
 - Emission-lines moments are remeasured after the emission-line model
   fitting to account for the adjustments made to the continuum during
   that process.
 - During the emission-line modeling, fit metrics are now computed in
   15-pixel wide windows around each fitted line.  The metrics are
   amplitude, amplitude-to-noise, chi-square, root-mean-square
   difference with the model, and fractional RMS wrt the model.  These
   data are currently only in the reference files, not the output maps.
 - Begun coding for switching the templates between what is used for the
   stellar kinematics fitting and what's used for the emission-line
   modeling.  Coding incomplete, but does not effect runs that use the
   same template sets for both analysis steps (as done in DR15).
 - MILES-HC templates changed to reflect re-run of the analysis by M.
   Cappellari so that we have documented which MILES template is
   incorporated into each MILES-HC template.  The total number of
   templates was again 49; however, 7 of these have been removed because
   they had artifacts and/or emission-lines in their spectra.  See
   `./data/spectral_templates/miles_cluster/README` .
 - Included redshift fix file that can be used to replace the redshift
   used for any PLATEIFU either because it's not provided by the
   relevant ancillary programs or the NSA redshift is incorrect.
 - Corrected error in the emission-line fluxes for the fact that the
   input template flux increases with the fitted redshift.  So: `F_true
   = (1+z)*F_MPL6`.
 - Corrected the emission-line equivalent width measurements from
   observed to rest frame.
 - Corrected the absorption-line spectral-index equivalent widths
   (angstrom units) from observed to rest frame.
 - Emission-line moments and spectral-indices (anything depending on the
   integration over a set of bandpasses) are now flagged as UNRELIABLE
   in the output maps if *any* pixels are masked in *any* their
   respective band passes.  This is true for both the main passband and
   any sidebands.
 - Spectral index measurements are now provided for the HYB binning
   scheme.
 - Spectral-index velocity-dispersion corrections are done adopting the
   NSA redshift for all spectra.  This means that there will be a
   *velocity* dependence of the correction that can be more significant
   that the correction for the velocity dispersion itself.  The approach
   to the velocity dispersion corrections should be improved in future
   releases.
    
2.1.3 (21 Nov 2017)
-------------------
 - Minor bug fix

2.1.2 (21 Nov 2017)
-------------------
 - Calculation of binned flux statistics now uses response function
   provided to the ReductionAssessments class
 - Binned flux statistics now determined before applying the reddening
   correction; now matches individual spaxel measurements.
 - Fixed some roundoff errors in the `SPX_*` extensions due to
   conversion from double to single-precision floats in output MAPS
 - Removed resolution-vector-based stellar velocity dispersion
   corrections from output MAPS files.  Only empirical (more robust)
   correction now provided, but resolution-vector calculation still in
   reference files.
 - Propagated this change though to the ppxf and spotcheck qa plots
 - Added `EMLINE_TPLSIGMA` to output MAPS files, providing the line
   width of the templates used during Sasuke.  Also added this into the
   measurements returned by pPXF such that `EMLINE_GSIGMA` and
   `EMLINE_INSTSIGMA` have the same meaning as in MPL-5: `EMLINE_GSIGMA`
   is the effective observed sigma of the line and `EMLINE_INSTSIGMA` is
   the instrumental dispersion at the line centroid.

2.1.1 (13 Nov 2017)
-------------------
 - Added check that pPXF returns some non-zero weights in
   stellar-continuum fitting.  If all weights are zero, the pPXF status
   is set to -2.
 - Check if all stellar-continuum fits are bad.  If so, emission-line
   model set fixed stellar kinematics to input redshift and dispersion =
   100 km/s, and the DAP quality flag set to CRITICAL.
 - Fixed fitting window determination in Sasuke
 - Allow for fit to proceed if input velocity is > -500 km/s (previously
   > 0 km/s).
 - Dummy indices are flagged as NOVALUE in output MAPS file
 - Any IVAR values that are 0 in 
 - Changed QA plotting scripts to account for bug in numpy/matplotlib

2.1   (27 Oct 2017)
-------------------
 - Included code for DAPall file
 - Fixed error in bin area calculation in unbinned data.
 - Added convenience functions to DAPFits for use with the MAPS files.
 - Fixed problem with velocity offset between template and galaxy data
   due to difference in pixel sampling.
 - Updated mangampl.py to include matplotlib version; added MPL-6
   versions
 - Fixed how the stellar-continuum-fitting masks are saved; still need
   to make sure users can reconstruct these exactly from the output
   model cube files.
 - Changed how spectral step is calculated to avoid numerical precision
   issues.
 - Incorporated M. Cappellari's `gaussian_filter1d` into the
   `convolution_variable_sigma` function in
   python/mangadap/util/instrument.py.  It's 2 orders of magnitude
   faster now!
 - Include spectral resolution handling hooks into spatial binning and
   stellar continuum fitting.  The code is compliant with data that both
   does and does not have the 'DISP' extension in the DRP cube file for
   backwards compatibility.
 - "Database" parameter classes (e.g., EmissionLineDB) have been moved
   from python/mangadap/par to python/mangadap/proc/.
 - Incorporated spaxel-by-spaxel instrumental dispersion corrections.
 - Significantly restructuring of most of the core classes to minimize
   the size of the reference files and to use common functions where
   possible (e.g., new DAPFitsUtil class).
 - Added equivalent width measurements for the Gaussian fitted fluxes;
   core function identical for the summed and Gaussian-fit fluxes.
 - Moved construction of main output files to dapfits.py
 - Include a map mask in reference files
 - Allow rundap to read the plate, ifu, and mode combinations to analyze
   from a file
 - Include filtering approach in stellar kinematics fitting
 - Incorporate ppxf v6.0.6
 - Reading/writing of fits data is now done using DAPFitsUtil.  When
   reading with this utility, the fits image data are automatically
   restructured into a row-major ordering to ease the confusion between
   the way astropy.io.fits provides the data and the ordering of the
   NAXIS keywords.  At the moment, this restructuring does not mess up
   the WCS coordinates in the header; however, depending on any
   development of astropy.io.fits in the future, this may become a
   problem.  This feels like a big mess at the moment, but there's no
   easy way around it.
 - Restructured how the covariance data are read and written to adhere
   to the format already adopted by the DRP.  The Covariance class
   should now be able to seemlessly read the broad-band correlation
   matrices provided in the LOGCUBE files.
 - Added DefaultConfig class to help manage config files.  Edited the
   StellarContinuumModelDef parameters to accommodate the new filtering
   approach.  Removed all but one StellarContinuumModel method.
 - Include an estimate of the stellar velocity dispersion correction
   based on fitting the native resolution optimized template to the same
   spectrum that has been rsolution matched to object spectrum being
   fit.
 - Changed velocity dispersion correction calculation for stellar
   kinematics: Old version used mean resolution difference between
   galaxy and template spectra.  New version uses a fit of the
   best-fitting model template to a resolution-matched version of
   itself. 
 - MAPS and LOGCUBE files set to single-precision float for survey
   output.
 - BINID in MAPS and LOGCUBE files set to int32 type for survey output.
 - Added NeIII and higher-order Balmer series lines (up to H10).
 - General improvements to the documentation
 - New emission-line module Sasuke that uses pPXF to fit the lines.
 - Adjusted pPXF filtering approach: remove final, fixed-kinematics fit
   and set polynomial order of initial fit to be the effective order of
   the filter (full length/smoothing box-1)
 - Calculate and include nominal errors in template weights in output
   reference file for stellar continuum model.
 - HARDCODED: Use first moment of H-alpha as the initial velocity guess
   for the emission lines.
 - Include pre-pixelized Gaussian extension from DRP.
 - Use 'z' column in plateTargets files to construct input parameter
   files instead of the `nsa_z` column; this should catch the
   ancillaries not in the NSA but with redshifts in the targeting
   catalog.
 - In rundap, wait for cluster processes to finish then also submit a
   script to construct the DAPall file
 - Updated versions of ppxf and `voronoi_2d_binning` and moved these
   files out of the mangadap package and into a new captools directory
 - Edited Sasuke to use iterations tested by Xihan Ji and Michele
   Cappellari.

2.0.2 (04 Aug 2016)
-------------------
 - Bug in setup for emission-line and spectral-index measurements when
   stellar continuum not fit.
 - Fixed error in return type when an emission-line fit fails, and
   report which line and spectrum failed.
 - Necessary conversions from asarray() to `atleast_1d()`
 - Include spot-check plotting code for maps
 - NEARBOUND is treated differently and separately applied to the
   stellar velocities and stellar velocity dispersions.  pPXF fits where
   the velocity dispersions at the lower bound do not *necessarily* have
   a bad velocity measurement.
 - Minor adjustments to TemplateLibrary and DRPComplete classes
 - Added python versions to mangampl.py and version checking to
   rundap.py; strict version compliance can be ignored if necessary
 - Added version information to dapcube.py and dapmaps.py
 - Corrected HDUCLAS2 for some MAPS extensions

2.0.1 (26 Jul 2016)
-------------------
 - Removed old DAP IDL code, except for `pro/mangadap_version.pro`
 - Included correction for Galactic extinction; current default is to
   use the O'Donnell (1994) coefficents for the CCM extinction law.
 - Added MILES cluster template group.
 - Allow for choice to *not* match the spectral resolution of the
   template to that of the galaxy data, and supply templates at higher
   spectral sampling.
 - Include stellar velocity dispersion quadrature correction.
 - Included v6.0.0 of pPXF; allows analytic calculation of LOSVD FFT for
   better performance at low velocity dispersion
 - Fixed error in spectral-index masks
 - Dispersion correction separated from spectral-index measurements in
   maps file
 - Spectral-index dispersion corrections not automatically applied
 - Added config to ignore resolution matching in spectral-index
   measurements (for testing)
 - Changed the lower limit of the emission-line dispersion to be 30
   km/s, which means the line has a FWHM ~1 pixel at minimum; upper
   bound changed to 800 km/s.
 - Definition of "near" bounds for bounded kinematic fits:
     - Velocities and higher moments: Cannot be nearer to a bound than
       1% of the full width of the allowed parameter space.
     - Velocity dispersions: Cannot be nearer to a bound than 1% of the
       LOG of the full width of the allowed parameter space.  The
       velocity limits for the emission-lines are dynamically defined by
       the lines in the window or to be +/- 25 angstroms about the line
       center, and no higher moments are fit at the moment.  So the only
       fixed limits are for the dispersion; any values within the ranges
       30.:31, 774.16:800 are considered "near" the boundary.
     - For the stellar kinematics, this means:

       | PAR   |   `NEAR_BOUND` ranges   |
       |:-----:| -----------------------:|
       |    V  | -2000:-1960, 1960:2000  |
       | sigma | 6.903:7.255, 951.46:1e3 |
       |  hN   | -0.3:-0.294, 0.294:0.3  |

 - Stellar-continuum model mask now includes `NEAR_BOUND` bits, if
   necessary.
 - Constructs an output model data cube with the stellar continuum and
   emission-line model as a complement to the maps file.
 - Allows for an iteration that modifies the stellar-continuum fit by
   unmasking and fitting the emission lines, while not changing the
   stellar kinematics.
 - Fixed problem with emission-line parameters returned for non-primary
   lines fit in each window.
 - Fixed error in construction of the Hessian (inverse covariance)
   matrix from the emission-line fits.
 - Removed `SPXL_SPCOV` extension; changed extension names from 'SPXL'
   to 'SPX' to match binning scheme name.
 - Bookkeeping error in emission-line moment measurements
 - Incorrect/incomplete checking for sufficient data to fit emission
   lines

2.0   (08 Jun 2016)
-------------------
 - Initial tag after major refactor of DAP to be solely in python.  Some
   directory cleanup performed; however, the pro/ directory is still
   included in this tag.  All further tags will have the pro/ directory
   removed.
 - This tag was used to run initial benchmark tests for the candidate
   MPL-5 data.

1.1.2 (24 Feb 2016)
-------------------
 - Tag stable version of trunk before starting major revisions to python
   directory for MPL-5.

1.1.1 (25 Jan 2016)
-------------------
 - Tagged version run for MPL-4.
 - Tag was created much later than MPL-4 was released.  DAP MPL-4
   release located in the mangawork/sandbox directory, and not
   maintained in the nominal `MANGA_SPECTRO_ANALYSIS` directory
   structure because run on ICG, Portsmouth, cluster.  Differences
   between tagged version and version that produced MPL-4 are minor and
   largely contained within the python/mangadap/plot directory.  Edits
   to plotting routines have been used for plots of MPL-4 data shown in
   Marvin.

v1_0_0 (28 Apr 2015)
--------------------
 - Tagged version run on MPL-3 DRP fits files. 
 - Altered python directory structure in anticipation of a python
   rewrite of DAP
 - Stable version v0.94 merged from `kyle_testing` to trunk.
 - Basic integration at Utah, with some minor edits.
 - Spectral index measurements added, most known bugs address.  QA
   ongoing!
 - Added spectral index measurements; known bugs. QA ongoing!  
 - Completed new structure and data model
 - Completed edits through Block 3
 - Completed edits through Block 2
 - Check in using svn in working directory
 - Database structure
 - New environment scripts: `mdap_setup.sh`, `mdap_environment.sh`
 - Basic compilation errors in `manga_dap.pro`
 - New code: `mdap_match_obs_nsa.pro`, `mdap_create_input_table.pro`
 - Formatting of `manga_drp.pro`
 - Copied from `v0_8` by L. Coccato
 - Basic formatting of `map_read_datacube.pro` and `manga_drp.pro`
 - Adjustments to deal with DRP `v1_0_0` headers (different from
   prototype headers)
 - Begin KBW Development version 

