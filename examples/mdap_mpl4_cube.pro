;===============================================================================
; MPL-4 CUBE FILE PLANS
;===============================================================================
;   1: STON: Bin to S/N=30; only include S/N>5 spectra (r-band continuum)
;       MILES-THIN (219 stars drawn from a 12x12x12 grid in logT, logg, and [Fe/H])
;           stellar-cont (4 moments)
;           emission-line
;           abs-indices
;   2: NONE: Fit all spaxels with a continuum S/N>0, no binning
;       MIUSCAT-THIN (72 templates: Salpeter IMF; every other age <14 Gyr; three metallicities)
;           stellar-cont (2 moments)
;           emission-line
;           abs-indices
;   3: RADIAL: 10 logarithmically spaced radial bins; S/N>0
;       **No velocity registration**
;       MIUSCAT-THIN (72 templates: Salpeter IMF; every other age <14 Gyr; three metallicities)
;           stellar-cont (2 moments)
;           emission-line
;           abs-indices
;   4: RADIAL: 10 logarithmically spaced radial bins; S/N>0
;       Use velocity registration (based on velocities from iter 2)
;       MIUSCAT-THIN (72 templates: Salpeter IMF; every other age <14 Gyr; three metallicities)
;           stellar-cont (2 moments)
;           emission-line
;           abs-indices
;   5: ALL: Bin all "valid" spaxels; S/N>0
;       Use velocity registration (based on velocities from iter 2)
;       MILES-THIN (219 stars drawn from a 12x12x12 grid in logT, logg, and [Fe/H])
;           stellar-cont (2 moments)
;           emission-line
;           abs-indices
;   6: ALL: Bin all "valid" spaxels; S/N>0
;       Use velocity registration (based on velocities from iter 2)
;       *Full MILES*
;           stellar-cont (2 moments)
;           emission-line
;           abs-indices
;   7: ALL: Bin all "valid" spaxels; S/N>0
;       Use velocity registration (based on velocities from iter 2)
;       MIUSCAT-THIN (72 templates: Salpeter IMF; every other age <14 Gyr; three metallicities)
;           stellar-cont (2 moments)
;           emission-line
;           abs-indices
;   8: ALL: Bin all "valid" spaxels; S/N>0
;       **No velocity registration**
;       MIUSCAT-THIN (72 templates: Salpeter IMF; every other age <14 Gyr; three metallicities)
;           stellar-cont (2 moments)
;           emission-line
;           abs-indices
;   9: STON: Bin to S/N=5; only include S/N>0 spectra (r-band continuum)
;       MIUSCAT-THIN (72 templates: Salpeter IMF; every other age <14 Gyr; three metallicities)
;           stellar-cont (2 moments); **NO multiplicative polynomial**
;           emission-line
;           abs-indices

PRO CREATE_MANGA_DAP_EXECUTION_PLAN, $
                ofile, overwrite=overwrite

        ;---------------------------------------------------------------
        ; Define the number of execution iterations and setup the needed vectors
        ; and allocate the necessary arrays.
        niter = 9                                       ; Number of ExecutionPlans to produce

        MDAP_ALLOCATE_EXECUTION_PLAN_VARIABLES, niter, bin_par, w_range_sn, threshold_ston_bin, $
                                                w_range_analysis, threshold_ston_analysis, $
                                                analysis, tpl_lib_analysis, ems_par_analysis, $
                                                abs_par_analysis, analysis_par, analysis_prior, $
                                                overwrite_flag

        ;---------------------------------------------------------------
        ; Global
        threshold_ston_analysis[*] = 0.0d
        analysis[*,0] = 'stellar-cont'
        analysis[*,1] = 'emission-line'
        analysis[*,2] = 'abs-indices'
        ems_par_analysis[*] = 'STANDARD'
        abs_par_analysis[*] = 'LICK'
        analysis_par[*].degree = -1
        analysis_par[*].mdegree = 6
        analysis_par[*].reddening_order = 0
        overwrite_flag[*] = 1

        ;---------------------------------------------------------------
        ; Plan 1
        bin_par[0].type = 'STON'
        bin_par[0].ston = 30.0d
        bin_par[0].noise_calib = 1
        w_range_sn[0,*] = [5560.00, 6942.00]
        threshold_ston_bin[0] = 5.0d
        w_range_analysis[0,*] = [3650.,10300.] 
        tpl_lib_analysis[0] = 'MILES-THIN'
        analysis_par[0].moments = 4
        analysis_prior[0] = ''

        ;---------------------------------------------------------------
        ; Plan 2
        bin_par[1].type = 'NONE'
        w_range_sn[1,*] = [5560.00, 6942.00]
        threshold_ston_bin[1] = 0.0d
        w_range_analysis[1,*] = [3650.,10300.] 
        tpl_lib_analysis[1] = 'MIUSCAT-THIN'
        analysis_par[1].moments = 2
        analysis_prior[1] = ''

        ;---------------------------------------------------------------
        ; Plan 3
        bin_par[2].type = 'RADIAL'
        bin_par[2].noise_calib = 1
        bin_par[2].nr = 10
        bin_par[2].rlog = 1
        w_range_sn[2,*] = [5560.00, 6942.00]
        threshold_ston_bin[2] = 0.0d
        w_range_analysis[2,*] = [3650.,10300.] 
        tpl_lib_analysis[2] = 'MIUSCAT-THIN'
        analysis_par[2].moments = 2
        analysis_prior[2] = ''

        ;---------------------------------------------------------------
        ; Plan 4
        bin_par[3].type = 'RADIAL'
        bin_par[3].v_register = 1
        bin_par[3].noise_calib = 1
        bin_par[3].nr = 10
        bin_par[3].rlog = 1
        w_range_sn[3,*] = [5560.00, 6942.00]
        threshold_ston_bin[3] = 0.0d
        w_range_analysis[3,*] = [3650.,10300.] 
        tpl_lib_analysis[3] = 'MIUSCAT-THIN'
        analysis_par[3].moments = 2
        analysis_prior[3] = '1'

        ;---------------------------------------------------------------
        ; Plan 5
        bin_par[4].type = 'ALL'
        bin_par[4].v_register = 1
        bin_par[4].noise_calib = 1
        w_range_sn[4,*] = [5560.00, 6942.00]
        threshold_ston_bin[4] = 0.0d
        w_range_analysis[4,*] = [3650.,10300.] 
        tpl_lib_analysis[4] = 'MILES-THIN'
        analysis_par[4].moments = 2
        analysis_prior[4] = '1'

        ;---------------------------------------------------------------
        ; Plan 6
        bin_par[5].type = 'ALL'
        bin_par[5].v_register = 1
        bin_par[5].noise_calib = 1
        w_range_sn[5,*] = [5560.00, 6942.00]
        threshold_ston_bin[5] = 0.0d
        w_range_analysis[5,*] = [3650.,10300.] 
        tpl_lib_analysis[5] = 'MILES'
        analysis_par[5].moments = 2
        analysis_prior[5] = '1'

        ;---------------------------------------------------------------
        ; Plan 7
        bin_par[6].type = 'ALL'
        bin_par[6].v_register = 1
        bin_par[6].noise_calib = 1
        w_range_sn[6,*] = [5560.00, 6942.00]
        threshold_ston_bin[6] = 0.0d
        w_range_analysis[6,*] = [3650.,10300.] 
        tpl_lib_analysis[6] = 'MIUSCAT-THIN'
        analysis_par[6].moments = 2
        analysis_prior[6] = '1'

        ;---------------------------------------------------------------
        ; Plan 8
        bin_par[7].type = 'ALL'
        bin_par[7].noise_calib = 1
        w_range_sn[7,*] = [5560.00, 6942.00]
        threshold_ston_bin[7] = 0.0d
        w_range_analysis[7,*] = [3650.,10300.] 
        tpl_lib_analysis[7] = 'MIUSCAT-THIN'
        analysis_par[7].moments = 2
        analysis_prior[7] = ''

        ;---------------------------------------------------------------
        ; Plan 9
        bin_par[8].type = 'STON'
        bin_par[8].ston = 5.0d
        bin_par[8].noise_calib = 1
        w_range_sn[8,*] = [5560.00, 6942.00]
        threshold_ston_bin[8] = 0.0d
        w_range_analysis[8,*] = [3650.,10300.] 
        tpl_lib_analysis[8] = 'MIUSCAT-THIN'
        analysis_par[8].moments = 2
        analysis_par[8].mdegree = 0
        analysis_prior[8] = ''

        ;---------------------------------------------------------------
        ;---------------------------------------------------------------

        ; Write the parameter file
        MDAP_WRITE_EXECUTION_PLANS, ofile, bin_par, w_range_sn, threshold_ston_bin, $
                                    w_range_analysis, threshold_ston_analysis, analysis, $
                                    tpl_lib_analysis, ems_par_analysis, abs_par_analysis, $
                                    analysis_par, analysis_prior, overwrite_flag, $
                                    overwrite=overwrite
END

