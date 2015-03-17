;+
; NAME:
;       MDAP_ALLOCATE_EXECUTION_PLAN_VARIABLES
;
; PURPOSE:
;       Allocate the necessary space for the execution plan variables.
;
; CALLING SEQUENCE:
;       MDAP_ALLOCATE_EXECUTION_PLAN_VARIABLES, niter, bin_par, w_range_sn, threshold_ston_bin, $
;                                               w_range_analysis, threshold_ston_analysis, $
;                                               analysis, tpl_lib_analysis, ems_par_analysis, $
;                                               abs_par_analysis, analysis_par, analysis_prior, $
;                                               overwrite_flag
;
; INPUTS:
;       niter int
;               Number of plans to allocate.
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;
; OUTPUT:
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
; OPTIONAL OUTPUT:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;       MDAP_DEFINE_BIN_PAR()
;       MDAP_DEFINE_ANALYSIS_PAR()
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;       17 Mar 2015: Original implementation by K. Westfall (KBW)
;-
;-----------------------------------------------------------------------

PRO MDAP_ALLOCATE_EXECUTION_PLAN_VARIABLES, $
                niter,  bin_par, w_range_sn, threshold_ston_bin, w_range_analysis, $
                threshold_ston_analysis, analysis, tpl_lib_analysis, ems_par_analysis, $
                abs_par_analysis, analysis_par, analysis_prior, overwrite_flag

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
END

