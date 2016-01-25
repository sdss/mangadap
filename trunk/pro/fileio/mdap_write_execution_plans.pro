;+
; NAME:
;   MDAP_WRITE_EXECUTION_PLANS
;
; PURPOSE:
;       Write the data necessary to a yanny parameter file that sets the
;       execution plans for the MaNGA DAP.
;
; CALLING SEQUENCE:
;       MDAP_WRITE_EXECUTION_PLANS, ofile, bin_par, w_range_sn, threshold_ston_bin, $
;                                   w_range_analysis, threshold_ston_analysis, analysis, $
;                                   tpl_lib_analysis, ems_par_analysis, abs_par_analysis, $
;                                   analysis_par, analysis_prior, overwrite_flag, execute_flag, $
;                                   version=version, overwrite=overwrite
;
; INPUTS:
;       ofile string
;               Name of the output file
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
;       /overwrite
;           Overwrite any existing file of the same name.
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
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;       MDAP_WRITE_EXECUTION_PLANS_STRUCT
;       MDAP_WRITE_EXECUTION_PLANS_INSTANCE
;
; REVISION HISTORY:
;       17 Mar 2015: Original implementation by K. Westfall (KBW)
;       30 Oct 2015: (KBW) Add the "execute_flag" keyword to allow plans
;                          to be skipped, while keeping the output file
;                          names dependent on the index of the input
;                          plan.  Expand the length of the template
;                          library, emission-line database, and
;                          spectral-index database names.
;-
;-----------------------------------------------------------------------

PRO MDAP_WRITE_EXECUTION_PLANS_STRUCT, $
                unit

        printf, unit, ''
        printf, unit, ''
        printf, unit, 'typedef struct {'
        printf, unit, '    char bin_type[6];'
        printf, unit, '    int bin_v_register;'
        printf, unit, '    int bin_optimal_weighting;'
        printf, unit, '    int bin_noise_calib;'
        printf, unit, '    double bin_ston;'
        printf, unit, '    double bin_cx;'
        printf, unit, '    double bin_cy;'
        printf, unit, '    double bin_pa;'
        printf, unit, '    double bin_ell;'
        printf, unit, '    double bin_rs;'
        printf, unit, '    double bin_re;'
        printf, unit, '    int bin_nr;'
        printf, unit, '    int bin_rlog;'
        printf, unit, '    double bin_rscale;'
        printf, unit, '    double w_range_sn[2];'
        printf, unit, '    double threshold_ston_bin;'
        printf, unit, '    double w_range_analysis[2];'
        printf, unit, '    double threshold_ston_analysis;'
        printf, unit, '    char analysis[4][20];'
        printf, unit, '    char tpl_lib_analysis[18];'
        printf, unit, '    char ems_par_analysis[18];'
        printf, unit, '    char abs_par_analysis[18];'
        printf, unit, '    int analysis_moments;'
        printf, unit, '    int analysis_degree;'
        printf, unit, '    int analysis_mdegree;'
        printf, unit, '    double analysis_reddening[2];'
        printf, unit, '    int analysis_reddening_order;'
        printf, unit, '    int analysis_zero_instr_disp;'
        printf, unit, '    char analysis_prior[300];'
        printf, unit, '    int overwrite_flag;'
        printf, unit, '    int execute_flag;'
        printf, unit, '} DAPPLAN;'
        printf, unit, ''
        printf, unit, ''

END

PRO MDAP_WRITE_EXECUTION_PLANS_INSTANCE, $
                unit, bin_par, w_range_sn, threshold_ston_bin, w_range_analysis, $
                threshold_ston_analysis, analysis, tpl_lib_analysis, ems_par_analysis, $
                abs_par_analysis, analysis_par, analysis_prior, overwrite_flag, execute_flag

        ; Perform some checks
        if strlen(tpl_lib_analysis) gt 18 then begin
            print, strlen(tpl_lib_analysis)
            message, 'Cannot accommodate Tpl Lib key with length greater than 18!'
        endif
        if strlen(ems_par_analysis) gt 18 then begin
            print, strlen(ems_par_analysis)
            message, 'Cannot accommodate emission-line key with length greater than 18!'
        endif
        if strlen(abs_par_analysis) gt 18 then begin
            print, strlen(abs_par_analysis)
            message, 'Cannot accommodate spectral-index key with length greater than 18!'
        endif
        if strlen(analysis_prior) gt 50 then begin
            print, strlen(analysis_prior)
            message, 'Cannot accommodate prior with length greater than 50!'
        endif
        if n_elements(w_range_sn) ne 2 || n_elements(w_range_analysis) ne 2 $
           || n_elements(analysis) ne 4 then begin
            print, n_elements(w_range_sn), n_elements(w_range_analysis), n_elements(analysis)
            message, 'Incorrect length of input vectors.'
        endif

        ; Set analysis_prior to NULL if it has no length
        prior = analysis_prior
        if strlen(analysis_prior) eq 0 then $
            prior = 'NULL'

        ; Set selected analysis to NULL if has no length
        analysis_name = analysis
        for i=0,3 do $
            if strlen(analysis[i]) eq 0 then $
                analysis_name[i] = 'NULL'

        ; Set selected libraries/parameter files to NULL if they have no length
        tpl_lib_name = tpl_lib_analysis
        if strlen(tpl_lib_analysis) eq 0 then $
            tpl_lib_name = 'NULL'
        ems_par_name = ems_par_analysis
        if strlen(ems_par_analysis) eq 0 then $
            ems_par_name = 'NULL'
        abs_par_name = abs_par_analysis
        if strlen(abs_par_analysis) eq 0 then $
            abs_par_name = 'NULL'

        ; Print the plan
        printf, unit, bin_par.type, bin_par.v_register, bin_par.optimal_weighting, $
                bin_par.noise_calib, bin_par.ston, bin_par.cx, bin_par.cy, bin_par.pa, $
                bin_par.ell, bin_par.rs, bin_par.re, bin_par.nr, bin_par.rlog, bin_par.rscale, $
                w_range_sn[0], w_range_sn[1], threshold_ston_bin, w_range_analysis[0], $
                w_range_analysis[1], threshold_ston_analysis, analysis_name[0], analysis_name[1], $
                analysis_name[2], analysis_name[3], tpl_lib_name, ems_par_name, abs_par_name, $
                analysis_par.moments, analysis_par.degree, analysis_par.mdegree, $
                analysis_par.reddening[0], analysis_par.reddening[1], $
                analysis_par.reddening_order, analysis_par.zero_instr_disp, prior, $
                overwrite_flag, execute_flag, $
                format='("DAPPLAN ", A6, 3(I2), 7(E14.6), I4, I2, E14.6, 2(" { ", 2(E14.6), ' $
                       + '" } ", E14.6), " { ", 4(A20), " } ", 3(A18), 3(I3), " { ", 2(E14.6), ' $
                       + '" } ", 2(I2), A301, 2(I2))'
END


PRO MDAP_WRITE_EXECUTION_PLANS, $
                ofile, bin_par, w_range_sn, threshold_ston_bin, w_range_analysis, $
                threshold_ston_analysis, analysis, tpl_lib_analysis, ems_par_analysis, $
                abs_par_analysis, analysis_par, analysis_prior, overwrite_flag, execute_flag, $
                version=version, overwrite=overwrite

        version_module = '0.1'                          ; Version number
        if n_elements(version) ne 0 then begin          ; set version and return
            version = version_module
            return
        endif

        ; Check if it exists and can be overwritten
        if file_test(ofile) eq 1 && keyword_set(overwrite) then begin
            print, 'WARNING: Overwriting existing ' + ofile + '!'
            file_delete, ofile
        endif else if file_test(ofile) eq 1 && ~keyword_set(overwrite) then $
            message, 'WARNING: ' + ofile + ' already exists!  Flag overwrite or provide plan file.'

        ; Open the file
        OPENW, unit, ofile, /get_lun

        ; Write the structure of the parameter file to the header
        MDAP_WRITE_EXECUTION_PLANS_STRUCT, unit

        ; Write each plan instance
        n_plans = n_elements(bin_par)                   ; Number of plans to execute
        for i=0,n_plans-1 do $
            MDAP_WRITE_EXECUTION_PLANS_INSTANCE, unit, bin_par[i], reform(w_range_sn[i,*]), $
                                                 threshold_ston_bin[i], $
                                                 reform(w_range_analysis[i,*]), $
                                                 threshold_ston_analysis[i], $
                                                 reform(analysis[i,*]), $
                                                 tpl_lib_analysis[i], ems_par_analysis[i], $
                                                 abs_par_analysis[i], analysis_par[i], $
                                                 analysis_prior[i], overwrite_flag[i], $
                                                 execute_flag[i]

        ; Close the file
        printf, unit, ''
        printf, unit, ''
        CLOSE, unit
        FREE_LUN, unit
END

