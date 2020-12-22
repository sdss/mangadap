
# Bad fits to OIII lines biased by broad H-alpha?
7991-1901-109
8713-9102-1171
11944-12704-1211    - illustrative
9000-1901-204       - illustrative
9500-1901-158       - illustrative
220 7963-3702-326   - meh
221 7963-3702-352   - bad
448 9090-9101-1066  - bad



# Traceback 1: ERROR: ValueError: not enough values to unpack (expected 2, got 0) [mangadap.contrib.xjmc]
11020  1902
 8612  1901

239 9487-9102-1032  - poor fits to both OIII and H-beta; independent of tying scheme





MPL-11 TODO
-----------

 - MILES-HC:
    - Reconstruct using same procedure used for MaStar spectra
    * Direct spectral resolution assessments

 - Stellar kinematics:
    - Get rid of the filter keyword from the output headers...
    - Try fitting with some high-resolution templates?

 - Emission lines:
    - Check construction of the continuum and emission-line models.
      Model - emline looks odd in some cases.  E.g. SPX-MILESHC-MASTARHC2/9036-3703-41
    - Datamodel for indicating which emission-line parameters are tied
    - FOM flag that identifies poorly fit lines
    * Test use of MaStar SSP models
    * Test use of reddening instead of multiplicative polynomials
    * Balmer-line (particularly H-beta) flux changes
    - Reference files:
        - Make more independent of stellar-continnuum module reference
          file
        - Change `BASE` extension to `CONTINUUM` extension.
    - Wavelength dependent sigma rejection during fit?

 - Spectral Indices:
    - Enable composite indices?

 - Model cube files
    - Make BINID extension identical to BINID in MAPS file

 - DAPall:
    - provide absolute H-alpha luminosity as a separate column
    - attenuation corrections to H-alpha luminosity, SFR
    - fraction of spectra with `sigma < sigma_corr`
    - change 1 Re computation?
    - errors

Tasks
-----

 - Get rid of ppxf_faults in sasuke.py

 - For the test galaxies:
    - 7815-6101: low mass blue
        - Check D4000
    - 8138-12704: high mass blue
        - Check D4000
        - Check pixels in upper left corner
    - 8131-3702: low mass red ; no emission lines
        - Check emission-line properties and errors
    - 8131-6102: high mass red ; no emission lines
        - Check emission-line properties and errors
    - 8329-1901: low mass red; some emission lines
    - 8258-6102: high mass red; some emission lines
    - Repeat observations?
    - Subsample of galaxies that contain the "representative" spectra
      from above?
 
    ~ Compare old and new measurements

 - Identify reason for bad H-alpha velocities in 8138-12704

 - Add flag for when *observed* gas dispersion is less than 0.7 pixels.
   This adds side-lobe peaks that are ~1% of the model peak, and the RMS
   difference wrt an pixelized Gaussian is about 0.3%.

 - Distribution of RMS difference between stellar continuum fit and
   emission-line continuum model over the common fitting region?

 - H-beta analysis when fitting a common wavelength range for MILESHC,
   MASTARHC, and MASTARHC2
    - (DONE) Repeat analysis of new data
    - (DONE) Test imposing limits on the H-beta dispersion relative to
      the H-alpha dispersion?

 - Check offset velocity used for gas vs. stellar velocity fields
    - Make sure they're the same and the former isn't referenced to the
      input H-alpha velocity moment.

 - Check stellar velocity dispersion corrections between the two
   different methods.

 - Analysis of MPL-10
    - H-beta masking
    X Spaxels with failed pPXF fits from emission-line module

 - Spin-up on determination of resolution MILES-HC spectra
    - Test accuracy of wavelength-dependent convolution in the
      undersampled limit -> implement analytic FFT prescription?
    - First fit with pPXF to get best template mix
    - Then fit with own code to get resolution

 - Analysis of repeat observations
    - What parameter space is covered by the repeat observations?
    X Build table of MPL-10 repeats (@manga:kbw/mpl10_repeats.db)
    - Construct "repeat files" for value comparison
    - Repeat plots for stellar kinematics from Overview
    - Construct similar plots for emission-line flux and kinematics

 - Test use of Maraston+ MaStar SSP models
    X Decide on down-selection of templates
    X Compare with MaStar-HC2, BC03, MILES-SSP, M11-MILES
    X Best-fitting polynomial
        - Compare multiplicative polynomial of order ?-25 and
          attenuation law
    - (?) Check resolution of MaStar SSP templates?
    - (?) Detailed list of what stars went into each SSP

 - Test fitting a non-MaNGA cube
    - SAMI galaxies with MaNGA data?

 - Check all flags are correctly in the sdssmaskbits.par file.

TODO
----

 - Change default to fix bad pixels in spectral resolution?
 - Add DAP versioning to reference files
 - Fix minor issue in moment and index spaxels fit in hybrid case,
   compared to emission-line modeling data.
 - Make velocity dispersion correction for spectral indices use template
   mix from fit to the bin instead of the individual spaxels?
 - Documentation of config examples are out of date!
 - Fix ELMREBIN in 2nd round of moment calculations (not sure this is
   still a problem...)
 - Make the package pip installable
 - Enable composite indices
 - Allow construction of DAP MAPS headers to indicate masks that don't
   have the same prepended extension name... (e.g., quality mask of
   `SPECINDEX_CORR` should be `SPECINDEX_MASK`).
 - Differentiate between mask provided for BF index and mask for main
   continuum value? (e.g., latter only cares about the side bands)
 - Deprecate `mangadap.survey.mangampl`

