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
        ; First plan:
        ;---------------------------------------------------------------
        ; Use S/N binning with uniform weighing and applying the noise
        ; calibration
        bin_par[0].type = 'STON'
        bin_par[0].noise_calib = 1
        bin_par[0].ston = 30.0d
        ;   leave everything else as default (no velocity registration)
        ; Define the wavelength range over which to calculate the S/N
        w_range_sn[0,*] = [5560.00, 6942.00]

        ; Define the S/N threshold to include spectrum in any bin
        threshold_ston_bin[0] = -300.0d

        ; Define the wavelength range over which to perform ALL analyses
        w_range_analysis[0,*] = [3650.,10300.] 

        ; Define the S/N threshold for spectra to analyze
        threshold_ston_analysis[0] = 0.0d

        ; Define the analyses to apply
        analysis[0,0] = 'stellar-cont'
;        analysis[0,1] = 'star+gas'
        analysis[0,1] = 'emission-line'
;        analysis[0,2] = 'abs-indices'

        ; Define the template library to use
;       tpl_lib_analysis[0] = 'M11-STELIB'
        tpl_lib_analysis[0] = 'M11-MILES'
        ; Define the emission-line list to use
        ems_par_analysis[0] = 'STANDARD'
        ; Define the spectral-index parameters to use
        abs_par_analysis[0] = 'LICK'

        ; Set additional parameters needed by the analysis modules
        ; The reddening order can be 0, 1, or 2
        ; IF NOT SET HERE, the default values are:
        ;       moments=2, degree=-1, mdegree=0, reddening_order=0
        analysis_par[0].moments = 4
        analysis_par[0].degree = -1
        analysis_par[0].mdegree = 6
        analysis_par[0].reddening_order = 0

        ; For this first iteration, apply no analysis prior
        analysis_prior[0] = ''      ; No prior for the first plan

        ; Set a flag to overwrite existing output: 1-yes, 0-no
        overwrite_flag[0] = 1
        ;---------------------------------------------------------------
        ;---------------------------------------------------------------


        ;---------------------------------------------------------------
        ; Second plan:
        ;---------------------------------------------------------------
        ; Try RADIAL using the results from the first ExecutionPlan to
        ; velocity register the data -> set v_register to true here and
        ; add the prior below.
        bin_par[1].type = 'RADIAL'
        bin_par[1].v_register = 1
        bin_par[1].noise_calib = 1
        bin_par[1].nr = 10
        bin_par[1].rlog = 1
        ;   leave everything else as default

        ; Define the wavelength range over which to calculate the S/N
        w_range_sn[1,*] = [5560.00, 6942.00]
        ; Define the S/N threshold to include spectrum in any bin
        threshold_ston_bin[1] = -300.0d
        ; Define the wavelength range over which to perform ALL analyses
        w_range_analysis[1,*] = [3650.,10300.] 
        ; Define the S/N threshold to perform analysis
        threshold_ston_analysis[1] = 0.0d

        ; Set the list of analyses to perform.  The available analysis steps are
        ; listed above.
        ; Define the analyses to apply
        analysis[1,0] = 'stellar-cont'
        analysis[1,1] = 'emission-line'

        ; Define the template library to use
;        tpl_lib_analysis[1] = 'M11-STELIB'
        tpl_lib_analysis[1] = 'M11-MILES'
        ; Define the emission-line list to use
        ems_par_analysis[1] = 'STANDARD'

        ; Set additional parameters needed by the analysis modules
        ; The reddening order can be 0, 1, or 2
        ; IF NOT SET HERE, the default values are:
        ;       moments=2, degree=-1, mdegree=-1, reddening_order=0
        analysis_par[1].moments = 4
        analysis_par[1].degree = -1
        analysis_par[1].mdegree = 6
        analysis_par[1].reddening_order = 0

        ; Use the result from the first plan as a prior for this one
        analysis_prior[1] = '0'

        ; Set a flag to overwrite existing output: 1-yes, 0-no
        overwrite_flag[1] = 1
        ;---------------------------------------------------------------
        ;---------------------------------------------------------------

        ; Write the parameter file
        MDAP_WRITE_EXECUTION_PLANS, ofile, bin_par, w_range_sn, threshold_ston_bin, $
                                    w_range_analysis, threshold_ston_analysis, analysis, $
                                    tpl_lib_analysis, ems_par_analysis, abs_par_analysis, $
                                    analysis_par, analysis_prior, overwrite_flag, $
                                    overwrite=overwrite
END

