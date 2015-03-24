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
        niter = 3                                       ; Number of ExecutionPlans to produce

        MDAP_ALLOCATE_EXECUTION_PLAN_VARIABLES, niter, bin_par, w_range_sn, threshold_ston_bin, $
                                                w_range_analysis, threshold_ston_analysis, $
                                                analysis, tpl_lib_analysis, ems_par_analysis, $
                                                abs_par_analysis, analysis_par, analysis_prior, $
                                                overwrite_flag

        ;---------------------------------------------------------------
        ; First plan:
        ;---------------------------------------------------------------
        bin_par[0].type = 'NONE'

        w_range_sn[0,*] = [5560.00, 6942.00]
        ; Limits the number of included spectra to those with S/N>5
        threshold_ston_bin[0] = 5.0d

        w_range_analysis[0,*] = [3650.,10300.] 
        ; because of threshold_ston_bin = 5, only spectra with S/N>5
        ; will be analyzed, can set this threshold to 0
        threshold_ston_analysis[0] = 0.0d

        analysis[0,0] = 'stellar-cont'
        analysis[0,1] = 'emission-line'

        tpl_lib_analysis[0] = 'MIUSCAT'
        ems_par_analysis[0] = 'STANDARD'

        analysis_par[0].moments = 2
        analysis_par[0].degree = -1
        analysis_par[0].mdegree = 6
        analysis_par[0].reddening_order = 0

        analysis_prior[0] = ''
        overwrite_flag[0] = 1
        ;---------------------------------------------------------------
        ;---------------------------------------------------------------

        ;---------------------------------------------------------------
        ; Second plan:
        ;---------------------------------------------------------------
        bin_par[1].type = 'RADIAL'
        bin_par[1].v_register = 1
        bin_par[1].noise_calib = 1
        bin_par[1].nr = 10
        bin_par[1].rlog = 1

        w_range_sn[1,*] = [5560.00, 6942.00]
        threshold_ston_bin[1] = 0.0d

        w_range_analysis[1,*] = [3650.,10300.] 
        threshold_ston_analysis[1] = 0.0d

        analysis[1,0] = 'stellar-cont'
        analysis[1,1] = 'emission-line'
        analysis[1,2] = 'abs-indices'

        tpl_lib_analysis[1] = 'MIUSCAT'
        ems_par_analysis[1] = 'STANDARD'
        abs_par_analysis[1] = 'LICK'

        analysis_par[1].moments = 2
        analysis_par[1].degree = -1
        analysis_par[1].mdegree = 6
        analysis_par[1].reddening_order = 0

        analysis_prior[1] = '0'
        overwrite_flag[1] = 1
        ;---------------------------------------------------------------
        ;---------------------------------------------------------------

        ;---------------------------------------------------------------
        ; Third plan:
        ;---------------------------------------------------------------
        bin_par[2].type = 'ALL'
        bin_par[2].v_register = 1
        bin_par[2].noise_calib = 1

        w_range_sn[2,*] = [5560.00, 6942.00]
        threshold_ston_bin[2] = 0.0d

        w_range_analysis[2,*] = [3650.,10300.] 
        threshold_ston_analysis[2] = 0.0d

        analysis[2,0] = 'stellar-cont'
        analysis[2,1] = 'emission-line'
        analysis[2,2] = 'abs-indices'

        tpl_lib_analysis[2] = 'MIUSCAT'
        ems_par_analysis[2] = 'STANDARD'
        abs_par_analysis[2] = 'LICK'

        analysis_par[2].moments = 2
        analysis_par[2].degree = -1
        analysis_par[2].mdegree = 6
        analysis_par[2].reddening_order = 0

        analysis_prior[2] = '0'
        overwrite_flag[2] = 1
        ;---------------------------------------------------------------
        ;---------------------------------------------------------------

        ; Write the parameter file
        MDAP_WRITE_EXECUTION_PLANS, ofile, bin_par, w_range_sn, threshold_ston_bin, $
                                    w_range_analysis, threshold_ston_analysis, analysis, $
                                    tpl_lib_analysis, ems_par_analysis, abs_par_analysis, $
                                    analysis_par, analysis_prior, overwrite_flag, $
                                    overwrite=overwrite
END

