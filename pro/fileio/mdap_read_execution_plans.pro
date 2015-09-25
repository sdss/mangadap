;+
; NAME:
;       MDAP_READ_EXECUTION_PLANS
;
; PURPOSE:
;       Read the data necessary for a set of DAP execution plans from a
;       provided yanny parameter file.
;
; CALLING SEQUENCE:
;       MDAP_READ_EXECUTION_PLANS, ifile, bin_par, w_range_sn, threshold_ston_bin, $
;                                  w_range_analysis, threshold_ston_analysis, analysis, $
;                                  tpl_lib_analysis, ems_par_analysis, abs_par_analysis, $
;                                  analysis_par, analysis_prior, overwrite_flag
;
; INPUTS:
;       ifile string
;               Name of the input yanny file
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
;       MDAP_ALLOCATE_EXECUTION_PLAN_VARIABLES
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;       17 Mar 2015: Original implementation by K. Westfall (KBW)
;-
;-----------------------------------------------------------------------

PRO MDAP_READ_EXECUTION_PLANS, $
                ifile, bin_par, w_range_sn, threshold_ston_bin, w_range_analysis, $
                threshold_ston_analysis, analysis, tpl_lib_analysis, ems_par_analysis, $
                abs_par_analysis, analysis_par, analysis_prior, overwrite_flag

        ; Check the file exists
        if file_test(ifile) eq 0 then $
            message, 'Input plan file unavailable!'

        data = yanny_readone(ifile, 'DAPPLAN')      ; Read the data
        n_plans = n_elements(data)                  ; Number of plans

        ; Allocate the needed arrays
        MDAP_ALLOCATE_EXECUTION_PLAN_VARIABLES, n_plans, bin_par, w_range_sn, threshold_ston_bin, $
                                                w_range_analysis, threshold_ston_analysis, $
                                                analysis, tpl_lib_analysis, ems_par_analysis, $
                                                abs_par_analysis, analysis_par, analysis_prior, $
                                                overwrite_flag

        ; Copy the parameter data to the output variables
        for i=0,n_plans-1 do begin

            bin_par[i].type = data[i].bin_type

            if bin_par[i].type eq 'NULL' then $
                message, 'Bin type cannot be NULL'

            bin_par[i].v_register = data[i].bin_v_register
            bin_par[i].optimal_weighting = data[i].bin_optimal_weighting
            bin_par[i].noise_calib = data[i].bin_noise_calib
            bin_par[i].ston = data[i].bin_ston
            bin_par[i].cx = data[i].bin_cx
            bin_par[i].cy = data[i].bin_cy
            bin_par[i].pa = data[i].bin_pa
            bin_par[i].ell = data[i].bin_ell
            bin_par[i].rs = data[i].bin_rs
            bin_par[i].re = data[i].bin_re
            bin_par[i].nr = data[i].bin_nr
            bin_par[i].rlog = data[i].bin_rlog
            bin_par[i].rscale = data[i].bin_rscale

            w_range_sn[i,*] = data[i].w_range_sn
            threshold_ston_bin[i] = data[i].threshold_ston_bin

            w_range_analysis[i,*] = data[i].w_range_analysis
            threshold_ston_analysis[i] = data[i].threshold_ston_analysis

            analysis[i,*] = data[i].analysis[*]

            for j=0,3 do $
                if analysis[i,j] eq 'NULL' then $
                    analysis[i,j] = ''

            tpl_lib_analysis[i] = data[i].tpl_lib_analysis
            if tpl_lib_analysis[i] eq 'NULL' then $
                tpl_lib_analysis[i] = ''
            ems_par_analysis[i] = data[i].ems_par_analysis
            if ems_par_analysis[i] eq 'NULL' then $
                ems_par_analysis[i] = ''
            abs_par_analysis[i] = data[i].abs_par_analysis
            if abs_par_analysis[i] eq 'NULL' then $
                abs_par_analysis[i] = ''

            analysis_par[i].moments = data[i].analysis_moments
            analysis_par[i].degree = data[i].analysis_degree
            analysis_par[i].mdegree = data[i].analysis_mdegree
            analysis_par[i].reddening = data[i].analysis_reddening
            analysis_par[i].reddening_order = data[i].analysis_reddening_order
            analysis_par[i].zero_instr_disp = data[i].analysis_zero_instr_disp

            analysis_prior[i] = data[i].analysis_prior
            if analysis_prior[i] eq 'NULL' then $
                analysis_prior[i] = ''
            overwrite_flag[i] = data[i].overwrite_flag

        endfor

END


