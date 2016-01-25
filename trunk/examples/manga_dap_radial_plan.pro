;===============================================================================
;===============================================================================
; For an explanation of how to setup this procedure for your desired
; execution plan(s), please see
; $MANGADAP_DIR/examples/README.execution_plans
;===============================================================================
;===============================================================================

PRO CREATE_MANGA_DAP_EXECUTION_PLAN, $
                ofile, overwrite=overwrite

        ;---------------------------------------------------------------
        ; Define the number of execution iterations and setup the needed vectors
        ; and allocate the necessary arrays.
        niter = 2                                       ; Number of ExecutionPlans to produce

        MDAP_ALLOCATE_EXECUTION_PLAN_VARIABLES, niter, bin_par, w_range_sn, threshold_ston_bin, $
                                                w_range_analysis, threshold_ston_analysis, $
                                                analysis, tpl_lib_analysis, ems_par_analysis, $
                                                abs_par_analysis, analysis_par, analysis_prior, $
                                                overwrite_flag

        ;---------------------------------------------------------------
        ;---------------------------------------------------------------
        bin_par[0].type = 'STON'
        bin_par[0].ston = 30.0d

        bin_par[1].type = 'RADIAL'
        bin_par[1].v_register = 1
        bin_par[1].nr = 10
        bin_par[1].rlog = 1

        w_range_sn[0,*] = [5560.00, 6942.00]
        w_range_sn[1,*] = [5560.00, 6942.00]

        threshold_ston_bin[*] = -300.0d

        w_range_analysis[0,*] = [3650.,10300.] 
        w_range_analysis[1,*] = [3650.,10300.] 

        threshold_ston_analysis[*] = 0.0d

        analysis[*,0] = 'stellar-cont'
        analysis[*,1] = 'emission-line'
        analysis[*,2] = 'abs-indices'

        tpl_lib_analysis[*] = 'M11-STELIB'
        ems_par_analysis[*] = 'STANDARD'
        abs_par_analysis[*] = 'LICK'

        analysis_par[*].moments = 4
        analysis_par[*].degree = -1
        analysis_par[*].mdegree = 6
        analysis_par[*].reddening_order = 0

        analysis_prior[0] = ''
        analysis_prior[1] = '0'

        overwrite_flag[*] = 1
        ;---------------------------------------------------------------
        ;---------------------------------------------------------------

        ; Write the parameter file
        MDAP_WRITE_EXECUTION_PLANS, ofile, bin_par, w_range_sn, threshold_ston_bin, $
                                    w_range_analysis, threshold_ston_analysis, analysis, $
                                    tpl_lib_analysis, ems_par_analysis, abs_par_analysis, $
                                    analysis_par, analysis_prior, overwrite_flag, $
                                    overwrite=overwrite
END

