;+
; NAME:
;       MDAP_SURVEY_EXECUTION_SETUP
;
; PURPOSE:
;       Create the execution plan elements used at the survey level.
;
; CALLING SEQUENCE:
;       MDAP_SURVEY_EXECUTION_SETUP, mode, bin_par, w_range_sn, threshold_ston_bin, $
;                                    w_range_analysis, threshold_ston_analysis, analysis, $
;                                    tpl_lib_analysis, ems_par_analysis, abs_par_analysis, $
;                                    analysis_par, analysis_prior, overwrite_flag, execute_flag
;
; INPUTS:
;       mode string
;               Output mode of the DRP.  Must be either 'CUBE' or 'RSS'.
;
;       bin_par BinPar[P]
;               Structure containing the parameters used to perform the
;               spatial binning.  See MDAP_DEFINE_BIN_PAR() for a
;               description of the structure.
;
;       w_range_sn dblarr[P][2]
;               Wavelength range used to calculate the S/N of each spectrum,
;               performed using each individual spectrum.  This is always
;               calculated, regardless of the binning type.
;
;       threshold_ston_bin dblarr[P]
;               The S/N threshold for the inclusion of any DRP spectrum in the
;               binning process.
;
;       w_range_analysis dblarr[P][2]
;               The wavelength range to use during the analysis.
;
;       threshold_ston_analysis dblarr[P]
;               The S/N threshold for the inclusion of a spectrum in the
;               analysis.
;
;       analysis strarr[P][]
;               List of analyses to perform.  Valid values are:
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
;       tpl_lib_analysis intarr[P]
;               Index of template library to use for each of the P execution
;               plans.  -1 if no template library is to be used during the
;               analysis.
;
;       ems_par_analysis intarr[P]
;               Index of emission-line parameter set to use for each of the P
;               execution plans.  -1 if no emission-line parameter set is to be
;               used during the analysis.
;
;       abs_par_analysis intarr[P]
;               Index of absorption-line parameter set to use for each of the P
;               execution plans.  -1 if no absorption-line parameter set is to
;               be used during the analysis.
;
;       analysis_par AnalaysisPar[P]
;               An array of structures that keeps the parameters used by
;               the analysis steps.  See MDAP_DEFINE_ANALYSIS_PAR().
;
;       analysis_prior strarr[P]
;               A string list of DAP files to use as a prior during each
;               execution plan.  On input, these are expected to contain
;               either the name of the DAP file directly (with either
;               the full path or the path within the calling directory)
;               or a string representation of an index of the plan to
;               use for the prior.  (**THIS FUNCTION**) replaces the
;               index numbered priors with the full file name, which is
;               what is included in the ExecutionPlan structure.
;
;       overwrite_flag intarr[P]
;               Flag to overwrite any existing output file with the exact same
;               execution plan.  TODO: Should this actually just be a flag that
;               forces analyses to be redone?
;
;       execute_flag intarr[P]
;               Flag to execute the plan (0-no;1-yes)
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;
; OUTPUT:
;
; OPTIONAL OUTPUT:
;       version string
;               Module version. If requested, the module is not executed and only
;               the version flag is returned
;
; COMMENTS:
;
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
;                   bin_par.noise_calib
;                       Apply a calibration of the noise vector in the
;                       binned spectrum.  Options are currently:
;                           
;                           0 - (default) Do not apply
;                           1 - N_calib = N_nominal * (1 + 1.6*log10(N_bin))
;
;                       where N_bin is the number of binned spectra.
;                       NOTE: Option 1 cannot be combined with optimal
;                       weighting!
;
;           bin_par.type = 'STON'
;               Spectra are binned, using the Voronoi binning scheme, to
;               a minimum S/N level. In addition to setting
;               bin_par.v_register, bin_par.optimal_weighting, and
;               bin_par.noise_calib (see above), this bin type also
;               requires:
;
;                   bin_par.ston
;                       Minimum S/N level.
;
;           bin_par.type = 'RADIAL'
;               Spectra are binned radially according to a provided
;               planar projection.  In addition to setting
;               bin_par.v_register, bin_par.optimal_weighting,
;               bin_par.noise_calib (see above under
;               bin_par.type='ALL'), this bin type also requires:
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
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;       MDAP_ALLOCATE_EXECUTION_PLAN_VARIABLES
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;       17 Mar 2015: Original implementation by K. Westfall (KBW)
;       30 Oct 2015: (KBW) Add the "execute_flag" keyword to allow plans
;                          to be skipped, while keeping the output file
;                          names dependent on the index of the input
;                          plan.  **THIS FILE IS VERY OUT OF DATE!!**
;-
;-----------------------------------------------------------------------


PRO MDAP_SURVEY_EXECUTION_SETUP, $
        mode, bin_par, w_range_sn, threshold_ston_bin, w_range_analysis, threshold_ston_analysis, $
        analysis, tpl_lib_analysis, ems_par_analysis, abs_par_analysis, analysis_par, $
        analysis_prior, overwrite_flag, execute_flag, version=version

        version_module = '0.1'                          ; Version number
        if n_elements(version) ne 0 then begin          ; set version and return
            version = version_module
            return
        endif

        ; Check the mode
        if mode ne 'RSS' && mode ne 'CUBE' then $
            message, 'Unknown MODE type: '+mode

        ;-----------------------------------------------------------------------
        ; Define the number of execution iterations and setup the needed vectors
        ; and allocate the necessary arrays.

        niter = 4                                       ; Number of ExecutionPlans to produce

        MDAP_ALLOCATE_EXECUTION_PLAN_VARIABLES, niter, bin_par, w_range_sn, threshold_ston_bin, $
                                                w_range_analysis, threshold_ston_analysis, $
                                                analysis, tpl_lib_analysis, ems_par_analysis, $
                                                abs_par_analysis, analysis_par, analysis_prior, $
                                                overwrite_flag, execute_flag

        ;-----------------------------------------------------------------------
        ; Nominal pipeline run at Utah

        if mode eq 'CUBE' then $
            bin_par[*].noise_calib = 1

        ; Run the pipeline
        bin_par[0].type = 'NONE'        ; ... without binning

        bin_par[1].type = 'ALL'         ; ... by binning all spectra

        bin_par[2].type = 'STON'        ; ... using a S/N=20 Voronoi binning
        bin_par[2].ston = 20.0d

        bin_par[3].type = 'RADIAL'      ; ... and with radial binning
        bin_par[3].v_register = 1
        bin_par[3].nr = 10
        bin_par[3].rlog = 1

        ; Calculate the S/N in the same region for all binning types
        w_range_sn[0,*] = [5560.00, 6942.00]
        w_range_sn[1,*] = [5560.00, 6942.00]
        w_range_sn[2,*] = [5560.00, 6942.00]
        w_range_sn[3,*] = [5560.00, 6942.00]
        threshold_ston_bin[*] = -300.0d             ; Include all data in S/N calculation

        ; Analyze the same wavelength region for all binning types
        w_range_analysis[0,*] = [3650.,10300.] 
        w_range_analysis[1,*] = [3650.,10300.] 
        w_range_analysis[2,*] = [3650.,10300.] 
        w_range_analysis[3,*] = [3650.,10300.] 
        threshold_ston_analysis[*] = 0.0d           ; Include only positive S/N in analysis

        ; Perform all 4 main analysis steps (GANDALF turned off)
        analysis[*,0] = 'stellar-cont'
;       analysis[*,1] = 'star+gas'
        analysis[*,2] = 'emission-line'
        analysis[*,3] = 'abs-indices'

        tpl_lib_analysis[*] = 'M11-STELIB'                 ; Always use the M11-STELIB library
        ems_par_analysis[*] = 'STANDARD'
        abs_par_analysis[*] = 'LICK'

        ; Always use the same analysis parameters
        analysis_par[*].moments = 4
        analysis_par[*].degree = -1
        analysis_par[*].mdegree = 6
        analysis_par[*].reddening_order = 0
        analysis_par[*].zero_instr_disp = 1     ; Do not use instr dispersion in GANDALF

        analysis_prior[0:2] = ''                ; No priors
        analysis_prior[3] = '2'                 ; ... for all but the RADIAL binning type

        overwrite_flag[*] = 1                   ; Always overwrite the existing data

        execute_flag[*] = 1                     ; Execute all plans

;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------

END

