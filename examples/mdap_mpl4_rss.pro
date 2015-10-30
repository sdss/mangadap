;===============================================================================
; MPL-4 RSS FILE PLANS
;===============================================================================
; All plans perform
;       - stellar-cont
;       - emission-line
;       - abs-indices
; Multiplicative polynomial is always 6th-order
; Only the S/N=30 plan (plan 1) fits 4 moments for the stellar
; kinematics; all others fit two.
;
; Dummy plans are inserted such that the indices match the CUBE plans
;
; USING M11-STELIB-ZSOL
;   1: DUMMY
;   2: DUMMY
;   3: DUMMY
;
;   4: RADIAL:
;           - 10 logarithmically spaced radial bins; S/N>0
;           - **No velocity registration**
;   5: RADIAL:
;           - 10 logarithmically spaced radial bins; S/N>0
;           - Velocity register based on velocities from iter 2
;   6: RADIAL:
;           - 2 bins: 0-1 Re, and 1-2 Re; all spectra with S/N>0 included
;           - **No velocity registration**
;   7: RADIAL:
;           - 2 bins: 0-1 Re, and 1-2 Re; all spectra with S/N>0 included
;           - Velocity register based on velocities from iter 2
;   8: ALL:
;           - Bin all "valid" spaxels; S/N>0
;           - **No velocity registration**
;   9: ALL:
;           - Bin all "valid" spaxels; S/N>0
;           - Velocity register based on velocities from iter 2
;
;  10: DUMMY
;
; USING MIUSCAT-THIN
;  11: DUMMY
;  12: DUMMY
;  13: DUMMY
;
;  14: RADIAL:
;           - 10 logarithmically spaced radial bins; S/N>0
;           - **No velocity registration**
;  15: RADIAL:
;           - 10 logarithmically spaced radial bins; S/N>0
;           - Velocity register based on velocities from iter 2
;  16: RADIAL:
;           - 2 bins: 0-1 Re, and 1-2 Re; all spectra with S/N>0 included
;           - **No velocity registration**
;  17: RADIAL:
;           - 2 bins: 0-1 Re, and 1-2 Re; all spectra with S/N>0 included
;           - Velocity register based on velocities from iter 2
;  18: ALL:
;           - Bin all "valid" spaxels; S/N>0
;           - **No velocity registration**
;  19: ALL:
;           - Bin all "valid" spaxels; S/N>0
;           - Velocity register based on velocities from iter 2
;
;  20: DUMMY
;
; USING MILES-THIN
;  21: DUMMY
;  22: DUMMY
;  23: DUMMY
;
;  24: RADIAL:
;           - 10 logarithmically spaced radial bins; S/N>0
;           - **No velocity registration**
;  25: RADIAL:
;           - 10 logarithmically spaced radial bins; S/N>0
;           - Velocity register based on velocities from iter 2
;  26: RADIAL:
;           - 2 bins: 0-1 Re, and 1-2 Re; all spectra with S/N>0 included
;           - **No velocity registration**
;  27: RADIAL:
;           - 2 bins: 0-1 Re, and 1-2 Re; all spectra with S/N>0 included
;           - Velocity register based on velocities from iter 2
;  28: ALL:
;           - Bin all "valid" spaxels; S/N>0
;           - **No velocity registration**
;  29: ALL:
;           - Bin all "valid" spaxels; S/N>0
;           - Velocity register based on velocities from iter 2

PRO CREATE_MANGA_DAP_EXECUTION_PLAN, $
                ofile, overwrite=overwrite

        ;---------------------------------------------------------------
        ; Define the number of execution iterations and setup the needed vectors
        ; and allocate the necessary arrays.
        niter = 29                                       ; Number of ExecutionPlans to produce

        MDAP_ALLOCATE_EXECUTION_PLAN_VARIABLES, niter, bin_par, w_range_sn, threshold_ston_bin, $
                                                w_range_analysis, threshold_ston_analysis, $
                                                analysis, tpl_lib_analysis, ems_par_analysis, $
                                                abs_par_analysis, analysis_par, analysis_prior, $
                                                overwrite_flag, execute_flag

        ;---------------------------------------------------------------
        ; Global
        w_range_sn[*,0] = 5560.00                   ; Calculate S/N in r-band
        w_range_sn[*,1] = 6942.00
        threshold_ston_bin[*] = 0.0d                ; Only the first plan sets a binning threshold
        bin_par[*].noise_calib = 1      ; Always perform the noise calibration (ignored for NONE)
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
;       execute_flag[*] = 1                         ; Always execute the plan

        ; Template selection!
        tpl_lib_analysis[0:9] = 'M11-STELIB-ZSOL'
        tpl_lib_analysis[10:19] = 'MIUSCAT-THIN'
        tpl_lib_analysis[20:28] = 'MILES-THIN'

        ; DO NOT EXECUTE THE DUMMY PLANS!!
        execute_flag[0:2] = 0
        execute_flag[9] = 0
        execute_flag[10:12] = 0
        execute_flag[19] = 0
        execute_flag[20:22] = 0
        ; ASSIGN DUMMY BINNING TYPES JUST TO FILL THE COLUMN
        bin_par[0:2].type = 'NULL'
        bin_par[9].type = 'NULL'
        bin_par[10:12].type = 'NULL'
        bin_par[19].type = 'NULL'
        bin_par[20:22].type = 'NULL'

        ; Turn on the desired template executions
        execute_flag[3:8] = 0           ; M11-STELIB-ZSOL
        execute_flag[13:18] = 0         ; MIUSCAT-THIN
        execute_flag[23:28] = 1         ; MILES-THIN

        i = 0

        ;---------------------------------------------------------------
        ; Plans 1, 11, 21 are DUMMY
        i++

        ;---------------------------------------------------------------
        ; Plan 2, 12, 22 are DUMMY
        i++

        ;---------------------------------------------------------------
        ; Plan 3, 13, 23 are DUMMY
        i++

        ;---------------------------------------------------------------
        ; Plan 4, 14, 24
        bin_par[i].type = 'RADIAL'
        bin_par[i].nr = 10
        bin_par[i].rlog = 1
;       execute_flag[i] = 0
        bin_par[i+10].type = 'RADIAL'
        bin_par[i+10].nr = 10
        bin_par[i+10].rlog = 1
;       execute_flag[i+10] = 0
        bin_par[i+20].type = 'RADIAL'
        bin_par[i+20].nr = 10
        bin_par[i+20].rlog = 1
;       execute_flag[i+20] = 0
        i++

        ;---------------------------------------------------------------
        ; Plan 5, 15, 25
        bin_par[i].type = 'RADIAL'
        bin_par[i].v_register = 1
        bin_par[i].nr = 10
        bin_par[i].rlog = 1
        analysis_prior[i] = '1'
;       execute_flag[i] = 0
        bin_par[i+10].type = 'RADIAL'
        bin_par[i+10].v_register = 1
        bin_par[i+10].nr = 10
        bin_par[i+10].rlog = 1
        analysis_prior[i+10] = '11'
;       execute_flag[i+10] = 0
        bin_par[i+20].type = 'RADIAL'
        bin_par[i+20].v_register = 1
        bin_par[i+20].nr = 10
        bin_par[i+20].rlog = 1
        analysis_prior[i+20] = '21'
;       execute_flag[i+20] = 0
        i++

        ;---------------------------------------------------------------
        ; Plan 6, 16, 26
        bin_par[i].type = 'RADIAL'
        bin_par[i].nr = 2
        bin_par[i].rs = 0.0             ; Start at R=0
        bin_par[i].re = 2.0             ; End at R = 2 Re
        bin_par[i].rscale = -1.         ; Use the effective radius from the parameter file!
;       execute_flag[i] = 0
        bin_par[i+10].type = 'RADIAL'
        bin_par[i+10].nr = 2
        bin_par[i+10].rs = 0.0             ; Start at R=0
        bin_par[i+10].re = 2.0             ; End at R = 2 Re
        bin_par[i+10].rscale = -1.         ; Use the effective radius from the parameter file!
;       execute_flag[i+10] = 1
        bin_par[i+20].type = 'RADIAL'
        bin_par[i+20].nr = 2
        bin_par[i+20].rs = 0.0             ; Start at R=0
        bin_par[i+20].re = 2.0             ; End at R = 2 Re
        bin_par[i+20].rscale = -1.         ; Use the effective radius from the parameter file!
;       execute_flag[i+20] = 0
        i++

        ;---------------------------------------------------------------
        ; Plan 7, 17, 27
        bin_par[i].type = 'RADIAL'
        bin_par[i].v_register = 1
        bin_par[i].nr = 2
        bin_par[i].rs = 0.0             ; Start at R=0
        bin_par[i].re = 2.0             ; End at R = 2 Re
        bin_par[i].rscale = -1.         ; Use the effective radius from the parameter file!
        analysis_prior[i] = '1'
;       execute_flag[i] = 0
        bin_par[i+10].type = 'RADIAL'
        bin_par[i+10].v_register = 1
        bin_par[i+10].nr = 2
        bin_par[i+10].rs = 0.0             ; Start at R=0
        bin_par[i+10].re = 2.0             ; End at R = 2 Re
        bin_par[i+10].rscale = -1.         ; Use the effective radius from the parameter file!
        analysis_prior[i+10] = '11'
;       execute_flag[i+10] = 0
        bin_par[i+20].type = 'RADIAL'
        bin_par[i+20].v_register = 1
        bin_par[i+20].nr = 2
        bin_par[i+20].rs = 0.0             ; Start at R=0
        bin_par[i+20].re = 2.0             ; End at R = 2 Re
        bin_par[i+20].rscale = -1.         ; Use the effective radius from the parameter file!
        analysis_prior[i+20] = '21'
;       execute_flag[i+20] = 0
        i++

        ;---------------------------------------------------------------
        ; Plan 8, 18, 28
        bin_par[i].type = 'ALL'
;       execute_flag[i] = 0
        bin_par[i+10].type = 'ALL'
;       execute_flag[i+10] = 0
        bin_par[i+20].type = 'ALL'
;       execute_flag[i+20] = 0
        i++

        ;---------------------------------------------------------------
        ; Plan 9, 19, 29
        bin_par[i].type = 'ALL'
        bin_par[i].v_register = 1
        analysis_prior[i] = '1'
;       execute_flag[i] = 0
        bin_par[i+10].type = 'ALL'
        bin_par[i+10].v_register = 1
        analysis_prior[i+10] = '11'
;       execute_flag[i+10] = 0
        bin_par[i+20].type = 'ALL'
        bin_par[i+20].v_register = 1
        analysis_prior[i+20] = '21'
;       execute_flag[i+20] = 0
        i++

        ;---------------------------------------------------------------
        ;---------------------------------------------------------------

        ; Write the parameter file
        MDAP_WRITE_EXECUTION_PLANS, ofile, bin_par, w_range_sn, threshold_ston_bin, $
                                    w_range_analysis, threshold_ston_analysis, analysis, $
                                    tpl_lib_analysis, ems_par_analysis, abs_par_analysis, $
                                    analysis_par, analysis_prior, overwrite_flag, execute_flag, $
                                    overwrite=overwrite
        

END

