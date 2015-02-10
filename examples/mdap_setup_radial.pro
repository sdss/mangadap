;===============================================================================
; SPATIAL BINNING
;===============================================================================
;
;       There are currently four types of spatial binning allowed by the
;       DAP.  This is a list of the options and the parameters in the
;       BinPar structure that MUST accompany any usage of each binning
;       type.  In practice, if bin_par.type is not 'ALL', 'STON', or
;       'RADIAL' it is assumed to be 'NONE'.
;
;           bin_par.type='NONE'
;               No binning is performed, all spectra are analyzed
;               (assuming they meet the analysis criteria).  No
;               parameters are required for this binning type.
;
;           bin_par.type='ALL'
;               All spectra are combined into a single bin for analysis.
;               The parameters required for this binning type are:
;
;                   bin_par.v_register
;                       Use a prior fit of the kinematics to velocity
;                       register the spectra to the median redshift
;                       before binning
;
;                   TODO: Allow for a bin_par.v_tweak keyword that will
;                   use a function to tweak the velocity registration.
;
;                   bin_par.optimal_weighting
;                       Weight spectra by S/N^2 during their combination
;                       (1-yes,0-no)
;
;           bin_par.type = 'STON'
;               Spectra are binned, using the Voronoi binning scheme, to
;               a minimum S/N level. In addition to setting
;               bin_par.v_register and bin_par.optimal_weighting (see
;               above), this bin type also requires:
;
;                   bin_par.ston
;                       Minimum S/N level.
;
;           bin_par.type = 'RADIAL'
;               Spectra are binned radially according to a provided
;               planar projection.  In addition to setting
;               bin_par.v_register and bin_par.optimal_weighting (see
;               above under bin_par.type='ALL'), this bin type also
;               requires:
;
;                   bin_par.cx, bin_par.cy:
;                       On-sky X,Y position to use for R=0.  (0,0 by
;                       default)
;
;                   bin_par.pa:
;                       Position angle from North (0 deg) through east
;                       (90 deg) of the MAJOR axis of a constant radius
;                       ellipse. (0 by default)
;
;                   bin_par.ell:
;                       Ellipticity (1-b/a) of a constant radius
;                       ellipse. (0 by default)
;
;                   bin_par.rs:
;                       Starting radius of the first bin (0 by default)
;
;                   bin_par.re:
;                       Ending radius of the last bin (largest R by
;                       default)
;
;                   bin_par.nr:
;                       Number of radial bins (5 by default)
;
;                   bin_par.rlog:
;                       Flag to logarithmically step the radial bins
;                       (1-true;0-false -- 0 by default)
;
;                   bin_par.rscale:
;                       Value by which to scale the radius in arcsec
;                       (1.0 by default)
;
;               The values of bin_par.pa and bin_par.ell WILL BE SET
;               AUTOMATICALLY using the ell and pa parameters read by
;               MANGA_DAP.  bin_par.rscale can be set to Reff
;               automatically by setting bin_par.rscale = -1.0.
;               
;               TODO: Add more options for rscale...
;
;       TODO: Other binning options maybe added later...
;
;===============================================================================
;===============================================================================
;
;===============================================================================
; STELLAR TEMPLATE LIBRARIES
;===============================================================================
;       The available stellar library keywords are:
;
;                                Spectral
;                 KEY    resolution (ang)                                   DATA
;       -------------    ----------------   ------------------------------------
;           M11-MARCS                2.73    dapsrc/external/templates/m11_marcs
;          M11-STELIB                3.40   dapsrc/external/templates/m11_stelib
;          M11-ELODIE                0.55   dapsrc/external/templates/m11_elodie
;           M11-MILES                2.54    dapsrc/external/templates/m11_miles
;               MILES                2.50        dapsrc/external/templates/miles
;              STELIB                3.40       dapsrc/external/templates/stelib
;
;       See dapsrc/pro/fileio/mdap_read_template_library.pro.  See
;       dapsrs/pro/usr/mdap_define_available_template_libraries.pro for
;       the definitions of the available libraries.
;===============================================================================
;===============================================================================
;
;===============================================================================
; EMISSION-LINE PARAMETER FILES
;===============================================================================
;
;       The available emission-line parameter files are:
;
;                 KEY     FILE
;       -------------    -------------------------------------------------------
;            STANDARD       dapsrc/external/manga_emission_line_list_nominal.par
;
;===============================================================================
;===============================================================================
;
;===============================================================================
; SPECTRAL-INDEX PARAMETER FILES
;===============================================================================
;
;       The available spectral-index parameter files are:
;
;                 KEY     FILE
;       -------------    -------------------------------------------------------
;                LICK     dapsrc/external/absorption_line_indices_definition.dat
;
;===============================================================================
;===============================================================================
;
;===============================================================================
; ANALYSIS STEPS
;===============================================================================
;       The currently available analysis steps are:
;
;           'stellar-cont'
;               Fit the stellar continuum using pPXF.  Requires the
;               template library.  Will mask emission lines, if a set of
;               GANDALF-style emission-line parameters are provided.
;               Determines the optimal template mix and the stellar
;               kinematics.
;
;           'star+gas'
;               This step uses GANDALF to fit the spectrum, which will
;               re-optimize the template mix, but does not adjust
;               the stellar kinematics.  Determines the emission-line
;               kinematics, intensities, fluxes, and equivalent widths.
;
;           'emission-line'
;               Fits an emission-line-only spectrum assuming a zero
;               baseline.  If either 'stellar-cont' or 'star+gas' have
;               been applied, the best-fitting stellar continuum (where
;               the GANDALF results take preference) will be subtracted
;               from the input spectra before fitting the
;               emission-line-only model.  If this has not been
;               performed, **no continuum is subtracted from the
;               spectra!**  This fit, therefore, does not require a
;               template library to be provided.  Currently, the lines
;               fit are hard-coded to be:
;
;               #   Line    Rest Wave (air)
;                   OII     3726.03
;                   Hb      4861.32
;                   OIII    4958.83
;                   OIII    5006.77
;                   NII     6547.96
;                   Ha      6562.80
;                   NII     6583.34
;                   SII     6716.31
;                   SII     6730.68
;
;               The fit is done twice using two different contributed
;               codes, one from Enci Wang and another from Francesco
;               Belfiore.  For these lines only, this step determines
;               the emission-line kinematics, intensities, fluxes, and
;               equivalent widths.
;
;           'abs-indices'
;               Calculates the absorption-line indices.  Requires an
;               absorption-line parameter set and a fit to the stellar
;               continuum.  If an emission-line model has been fit (with
;               the one from GANDALF taking precedence), it will be
;               subtracted from the galaxy spectra before performing the
;               measurements.  Other steps applied before performing the
;               index measurements is to replace aberrant pixels and to
;               match the spectral resolution to that of the index
;               system.  Index measurements performed on the
;               best-fitting template and its broadened version are used
;               to correct the measurements made for the data in a
;               differential way.  This step provides index measurements
;               and their errors.
;
; NOTE: The current list of fits must occur sequentially, except the
; emission-line only fit, which can be done without any continuum fit,
; although this is not recommended!  Based on the requested analysis
; steps, the code has logic in place that will determine what other
; steps need to be performed, and will do so after warning the user.
;===============================================================================
;===============================================================================
;
;===============================================================================
; ANALYSIS PARAMETERS
;===============================================================================
;
;       Some of the analysis steps have parameters that can be set to
;       specify certain aspects of their execution.  These parameters
;       are collected into a single structure called AnalysisPar and
;       defined in mdap_define_analysis_par.pro.  The currently
;       available parameters are:
;
;           AnalysisPar.moments
;               Number of *stellar* kinematics to fit.  The number of
;               emission-line kinematic moments is currently fixed at 2
;               (v, sigma).
;
;           AnalysisPar.degree
;               Degree of the *additive* polynomial used by pPXF when
;               fitting the stellar continuum.
;
;           AnalysisPar.mdegree
;               Degree of the *multiplicative* polynomial used by pPXF
;               and GANDALF when fitting the stellar continuum and full
;               spectrum, respectively.
;
;           AnalysisPar.reddening_order
;               The order of the reddening fit by GANDALF:  0 means no
;               reddening; 1 means fit a single dust-screen model; 2
;               means fit a dust-screen model plus a nebular only
;               component.  NOTE, if the reddening is fit, the value of
;               MDEGREE is ignored because fitting them both is
;               degenerate.
;
;           AnalysisPar.reddening
;               A two-element array used to set the initial guess for
;               the reddening coefficients.  The provided values must
;               match the input reddening_order, but the array must
;               always contain two elements.
;
;           AnalysisPar.zero_instr_disp
;               A flag to force GANDALF to ignore the instrumental
;               dispersion (set it to 0).  Flag is 0-false,1-true;
;               default is false.
;
;===============================================================================
;===============================================================================
;
;===============================================================================
; ANALYSIS PRIORS
;===============================================================================
;
;       When executing a certain execution plan, you may want to use
;       results from a previous run of the DAP as priors for another
;       execution plan.  You can do this through the definition of the
;       analysis_prior.  This is a string array that contains either the
;       file name of the DAP fits file to use (must use the full path or
;       the path from the execution directory) or the index number (in
;       string format) of the execution plan in the current list to use
;       as the prior.
;        
;       For example, for two plans where the first has no prior and the
;       second uses the first as a prior, you would define:
;
;           analysis_prior[0] = ''
;           analysis_prior[1] = '0'
;
;===============================================================================
;===============================================================================
;
; Setup some necessary execution variables for the MANGA DAP

; Signifier is what will be reported to the log file as the configuration file
; for a run of MaNGA_DAP

; dapsrc is an optional input to define the DAP source path instead of
; using environmental varaibles.

PRO MDAP_EXECUTION_SETUP, $
        signifier, bin_par, w_range_sn, threshold_ston_bin, w_range_analysis, $
        threshold_ston_analysis, analysis, tpl_lib_analysis, ems_par_analysis, abs_par_analysis, $
        analysis_par, analysis_prior, overwrite_flag, dapsrc=dapsrc, $
        save_intermediate_steps=save_intermediate_steps, $
        remove_null_templates=remove_null_templates, external_library=external_library

        ; Define the DAP source path
        if n_elements(dapsrc) eq 0 then $
            dapsrc = getenv('MANGADAP_DIR')

        ;-----------------------------------------------------------------------
        ; Flag to save intermediate steps.
        ; TODO: This is no longer used!

        ; save_intermediate_steps = 0 

        ;-----------------------------------------------------------------------
        ; Remove templates with zero weights in one fit from use in another fit.
        ; TODO: Currently not implemented.  Will include this as an option when
        ; applying priors.

        ; remove_null_templates = 1

        ;-----------------------------------------------------------------------
        ; Path to a library of fortran or C codes to be used.  If commented,
        ; internal IDL procedures are used.

        ; external_library=getenv('MANGADAP_DIR')+'/external/F90_32/'
        ; external_library=getenv('MANGADAP_DIR')+'/external/F90_64/'
        external_library=dapsrc+'/external/F90_64/'

        ;=======================================================================
        ; DEFINITION OF EXECUTION PROCEDURES

        ; Define a string used to signify this file in the header of the DAP output file(s)
;       cd, current=directory
;       signifier = directory+'/mdap_setup.pro'
        signifier = dapsrc+'/pro/usr/mdap_execution_setup.pro'

        ;-----------------------------------------------------------------------
        ; Define the number of execution iterations and setup the needed vectors
        ; and allocate the necessary arrays.

        niter = 2                                       ; Number of ExecutionPlans to produce

        bin_par_def = MDAP_DEFINE_BIN_PAR()             ; Define the BinPar structure
        bin_par = replicate( bin_par_def, niter)        ; Create the array of BinPar structures

        w_range_sn = dblarr(niter, 2)                   ; Wavelength range for S/N calculation
        threshold_ston_bin = dblarr(niter)              ; Threshold S/N to include spectrum in bin

        w_range_analysis = dblarr(niter, 2)             ; Wavelength range for the analysis
        threshold_ston_analysis = dblarr(niter)         ; Threshold S/N to analyze spectrum

        max_analysis_blocks = 4                         ; Maximum number of analysis blocks
        analysis = strarr(niter, max_analysis_blocks)   ; Analysis steps to apply

        tpl_lib_analysis = strarr(niter)                ; INDEX of template library to use
        ems_par_analysis = strarr(niter)                ; INDEX of emission-line parameter file
        abs_par_analysis = strarr(niter)                ; INDEX of absorption-line parameter file

        analysis_par_def = MDAP_DEFINE_ANALYSIS_PAR()   ; Define the AnalysisPar structure
        analysis_par = replicate( analysis_par_def, niter)  ; Create array of AnalysisPar structs

        analysis_prior = strarr(niter)                  ; Prior information used for analysis

        overwrite_flag = intarr(niter)                  ; Flag to overwrite any existing output file

        ;-----------------------------------------------------------------------
        ; For each iteration:
        bin_par[0].type = 'STON'
        bin_par[0].optimal_weighting = 1        ; Otherwise uniform weighting
        bin_par[0].ston = 40.0d
        ;   leave everything else as default (no velocity registration)

        ; Try RADIAL using the results from the first ExecutionPlan to
        ; velocity register the data -> set v_register to true here and
        ; add the prior below.
        bin_par[1].type = 'RADIAL'
        bin_par[1].v_register = 1
        bin_par[1].optimal_weighting = 1
        bin_par[1].nr = 10
        bin_par[1].rlog = 1
        ;   leave everything else as default

        ; Define the wavelength range over which to calculate the mean S/N per pixel
        w_range_sn[0,*] = [5560.00, 6942.00]
        w_range_sn[1,*] = [5560.00, 6942.00]

        ; Define the S/N threshold to include spectrum in any bin
        threshold_ston_bin[*] = -300.0d

        ; Define the wavelength range over which to perform ALL analyses
        w_range_analysis[0,*] = [3650.,10300.] 
        w_range_analysis[1,*] = [3650.,10300.] 

        ; Define the S/N threshold to perform analysis
        threshold_ston_analysis[*] = 0.0d

        ; Set the list of analyses to perform.  The available analysis steps are
        ; listed above.

        analysis[*,0] = 'stellar-cont'
        analysis[*,1] = 'star+gas'
        analysis[*,2] = 'emission-line'
        analysis[*,3] = 'abs-indices'

        ; Set the template library to use for the analysis
        ; TODO: Change this to use the library key?
        tpl_lib_analysis[*] = 'M11-STELIB'

        ; Set the emission-line parameter set to use
        ems_par_analysis[*] = 'STANDARD'

        ; Set the absorption-line parameter set to use
        abs_par_analysis[*] = 'LICK'

        ; Set additional parameters needed by the analysis modules
        ; The reddening order can be 0, 1, or 2
        ; TODO: Allow for oversample?
        ; IF NOT SET HERE, the default values are:
        ;       moments=2, degree=-1, mdegree=-1, reddening_order=0
        analysis_par[*].moments = 4
        analysis_par[*].degree = -1
        analysis_par[*].mdegree = 6
        analysis_par[*].reddening_order = 0
        ; analysis_par[0].reddening[*] = [0.01,0.01]

        ; Analysis priors, see description above.
        analysis_prior[0] = ''      ; No prior for the first plan
        analysis_prior[1] = '0'     ; Use the results from the first plan as a prior on the second

        ; Set a flag to overwrite existing output: 1-yes, 0-no
        overwrite_flag[0] = 1
;       overwrite_flag[0] = 0
        overwrite_flag[1] = 1

        ;=======================================================================
        ;=======================================================================

END

