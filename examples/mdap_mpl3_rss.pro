;===============================================================================
;===============================================================================
; For an explanation of how to setup this procedure for your desired
; execution plan(s), please see
; $MANGADAP_DIR/examples/README.execution_plans
;===============================================================================
;===============================================================================

;For RSS files:
;   1: 10 logarithmically spaced radial bins; S/N>0
;       **No velocity registration**
;       MIUSCAT-THIN (72 templates: Salpeter IMF; every other age <14 Gyr; three metallicities)
;           stellar-cont (2 moments)
;           emission-line
;           abs-indices
;   2: 10 logarithmically spaced radial bins; S/N>0
;       Use velocity registration (based on velocities from CUBE iter 2)
;       MIUSCAT-THIN (72 templates: Salpeter IMF; every other age <14 Gyr; three metallicities)
;           stellar-cont (2 moments)
;           emission-line
;           abs-indices
;   3: Bin all spaxels; S/N>0
;       Use velocity registration (based on velocities from CUBE iter 2)
;       MILES-THIN (219 stars drawn from a 12x12x12 grid in logT, logg, and [Fe/H])
;           stellar-cont (2 moments)
;           emission-line
;           abs-indices
;   4: Bin all spaxels; S/N>0
;       Use velocity registration (based on velocities from CUBE iter 2)
;       *Full MILES*
;           stellar-cont (2 moments)
;           emission-line
;           abs-indices
;   5: Bin all spaxels; S/N>0
;       Use velocity registration (based on velocities from CUBE iter 2)
;       MIUSCAT-THIN (72 templates: Salpeter IMF; every other age <14 Gyr; three metallicities)
;           stellar-cont (2 moments)
;           emission-line
;           abs-indices

PRO CREATE_MANGA_DAP_EXECUTION_PLAN, $
                ofile, overwrite=overwrite

        ;---------------------------------------------------------------
        ; Define the number of execution iterations and setup the needed vectors
        ; and allocate the necessary arrays.
        niter = 5                                       ; Number of ExecutionPlans to produce

        MDAP_ALLOCATE_EXECUTION_PLAN_VARIABLES, niter, bin_par, w_range_sn, threshold_ston_bin, $
                                                w_range_analysis, threshold_ston_analysis, $
                                                analysis, tpl_lib_analysis, ems_par_analysis, $
                                                abs_par_analysis, analysis_par, analysis_prior, $
                                                overwrite_flag

        ;---------------------------------------------------------------
        ; Global
        bin_par[*].noise_calib = 1
        threshold_ston_bin[*] = 0.0d
        threshold_ston_analysis[*] = 0.0d
        analysis[*,0] = 'stellar-cont'
        analysis[*,1] = 'emission-line'
        analysis[*,2] = 'abs-indices'
        ems_par_analysis[*] = 'STANDARD'
        abs_par_analysis[*] = 'LICK'
        analysis_par[*].moments = 2
        analysis_par[*].degree = -1
        analysis_par[*].mdegree = 6
        analysis_par[*].reddening_order = 0
        overwrite_flag[*] = 1

        ;---------------------------------------------------------------
        ; Plan 1
        bin_par[0].type = 'RADIAL'
        bin_par[0].nr = 10
        bin_par[0].rlog = 1
        w_range_sn[0,*] = [5560.00, 6942.00]
        w_range_analysis[0,*] = [3650.,10300.] 
        tpl_lib_analysis[0] = 'MIUSCAT-THIN'
        analysis_prior[0] = ''

        ;---------------------------------------------------------------
        ; Plan 2
        bin_par[1].type = 'RADIAL'
        bin_par[1].v_register = 1
        bin_par[1].nr = 10
        bin_par[1].rlog = 1
        w_range_sn[1,*] = [5560.00, 6942.00]
        w_range_analysis[1,*] = [3650.,10300.] 
        tpl_lib_analysis[1] = 'MIUSCAT-THIN'
        analysis_prior[1] = '1'

        ;---------------------------------------------------------------
        ; Plan 3
        bin_par[2].type = 'ALL'
        bin_par[2].v_register = 1
        w_range_sn[2,*] = [5560.00, 6942.00]
        w_range_analysis[2,*] = [3650.,10300.] 
        tpl_lib_analysis[2] = 'MILES-THIN'
        analysis_prior[2] = '1'

        ;---------------------------------------------------------------
        ; Plan 4
        bin_par[3].type = 'ALL'
        bin_par[3].v_register = 1
        w_range_sn[3,*] = [5560.00, 6942.00]
        w_range_analysis[3,*] = [3650.,10300.] 
        tpl_lib_analysis[3] = 'MILES'
        analysis_prior[3] = '1'

        ;---------------------------------------------------------------
        ; Plan 5
        bin_par[4].type = 'ALL'
        bin_par[4].v_register = 1
        w_range_sn[4,*] = [5560.00, 6942.00]
        w_range_analysis[4,*] = [3650.,10300.] 
        tpl_lib_analysis[4] = 'MIUSCAT-THIN'
        analysis_prior[4] = '1'

        ;---------------------------------------------------------------
        ;---------------------------------------------------------------

        ; Write the parameter file
        MDAP_WRITE_EXECUTION_PLANS, ofile, bin_par, w_range_sn, threshold_ston_bin, $
                                    w_range_analysis, threshold_ston_analysis, analysis, $
                                    tpl_lib_analysis, ems_par_analysis, abs_par_analysis, $
                                    analysis_par, analysis_prior, overwrite_flag, $
                                    overwrite=overwrite
END

