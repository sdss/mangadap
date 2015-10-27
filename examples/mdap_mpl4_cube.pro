;===============================================================================
; MPL-4 CUBE FILE PLANS
;===============================================================================
; Focus on a single template set:
;       - MIUSCAT-THIN to start
; All plans perform
;       - stellar-cont
;       - emission-line
;       - abs-indices
; Multiplicative polynomial is always 6th-order
; Only the S/N=30 plan (plan 1) fits 4 moments for the stellar
; kinematics; all others fit two.
;
;   1: STON:
;           - Bin to S/N=30; only include S/N>5 spectra (r-band continuum)
;           - stellar-cont (4 moments)
;   2: STON:
;           - Bin to S/N=5; including all spectra with S/N>0 (r-band continuum)
;   3: NONE:
;           - Fit all spaxels with a continuum S/N>0, no binning
;   4: RADIAL:
;           - 10 logarithmically spaced radial bins; S/N>0
;           - **No velocity registration**
;   5: RADIAL:
;           - 10 logarithmically spaced radial bins; S/N>0
;           - Velocity register based on velocities from iter 2
;   6: ALL:
;           - Bin all "valid" spaxels; S/N>0
;           - **No velocity registration**
;   7: ALL:
;           - Bin all "valid" spaxels; S/N>0
;           - Velocity register based on velocities from iter 2

PRO CREATE_MANGA_DAP_EXECUTION_PLAN, $
                ofile, overwrite=overwrite

        ;---------------------------------------------------------------
        ; Define the number of execution iterations and setup the needed vectors
        ; and allocate the necessary arrays.
        niter = 7                                       ; Number of ExecutionPlans to produce

        MDAP_ALLOCATE_EXECUTION_PLAN_VARIABLES, niter, bin_par, w_range_sn, threshold_ston_bin, $
                                                w_range_analysis, threshold_ston_analysis, $
                                                analysis, tpl_lib_analysis, ems_par_analysis, $
                                                abs_par_analysis, analysis_par, analysis_prior, $
                                                overwrite_flag

        ;---------------------------------------------------------------
        ; Global
        w_range_sn[*,0] = 5560.00                   ; Calculate S/N in r-band
        w_range_sn[*,1] = 6942.00
        threshold_ston_bin[i] = 0.0d                ; Only the first plan sets a binning threshold
        bin_par[i].noise_calib = 1      ; Always perform the noise calibration (ignored for NONE)
        analysis[*,0] = 'stellar-cont'              ; Perform all analyses except for GANDALF
        analysis[*,1] = 'emission-line'
        analysis[*,2] = 'abs-indices'
        w_range_analysis[*,0] = 3650.0              ; Analyze the full spectrum
        w_range_analysis[*,1] = 10300.0
        threshold_ston_analysis[*] = 0.0d           ; Analyze all spectra with positive S/N
        ems_par_analysis[*] = 'STANDARD'            ; Always use the standard emission-line mask
        abs_par_analysis[*] = 'LICK'                ; Always ust the standard index definitions
        analysis_par[*].moments = 2                 ; All but one plan uses two moments
        analysis_par[*].degree = -1                 ; Never fit an additive polynomial
        analysis_par[*].mdegree = 6                 ; Always use 6th order mult. polynomial
        analysis_par[*].reddening_order = 0         ; Never fit a reddening function
        analysis_prior[*] = ''                      ; All but two plans do not use a prior
        overwrite_flag[*] = 1                       ; Always overwrite any existing file

        ; Template selection!
        tpl_lib_analysis[*] = 'MIUSCAT-THIN'

        i = 0

        ;---------------------------------------------------------------
        ; Plan 1
        bin_par[i].type = 'STON'
        bin_par[i].ston = 30.0d
        threshold_ston_bin[i] = 5.0d
        analysis_par[i].moments = 4
        i++

        ;---------------------------------------------------------------
        ; Plan 2
        bin_par[i].type = 'STON'
        bin_par[i].ston = 5.0d
        i++

        ;---------------------------------------------------------------
        ; Plan 3
        bin_par[i].type = 'NONE'
        i++

        ;---------------------------------------------------------------
        ; Plan 4
        bin_par[i].type = 'RADIAL'
        bin_par[i].nr = 10
        bin_par[i].rlog = 1
        i++

        ;---------------------------------------------------------------
        ; Plan 5
        bin_par[i].type = 'RADIAL'
        bin_par[i].v_register = 1
        bin_par[i].nr = 10
        bin_par[i].rlog = 1
        analysis_prior[i] = '1'
        i++

        ;---------------------------------------------------------------
        ; Plan 6
        bin_par[i].type = 'ALL'
        i++

        ;---------------------------------------------------------------
        ; Plan 7
        bin_par[i].type = 'ALL'
        bin_par[i].v_register = 1
        analysis_prior[i] = '1'
        i++

        ;---------------------------------------------------------------
        ;---------------------------------------------------------------

        ; Write the parameter file
        MDAP_WRITE_EXECUTION_PLANS, ofile, bin_par, w_range_sn, threshold_ston_bin, $
                                    w_range_analysis, threshold_ston_analysis, analysis, $
                                    tpl_lib_analysis, ems_par_analysis, abs_par_analysis, $
                                    analysis_par, analysis_prior, overwrite_flag, $
                                    overwrite=overwrite
END

