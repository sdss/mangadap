
; TODO: Create procedures that read and write execution plans to the headers of
; fits files

PRO MDAP_GENERATE_OUTPUT_FILE_NAMES, $
                file_root, execution_plan

        n_plans = n_elements(execution_plan)
        for i=0,n_plans-1 do begin
            file=file_root+'BIN-'+execution_plan[i].bin_type+'-'+string(i+1,format='(I03)')+'.fits'
            execution_plan[i].ofile = file
        endfor
END

FUNCTION MDAP_SET_TPL_LIB_OUTPUT_FILE, $
                file_root, library_key
        return, file_root+library_key+'.fits'
END

PRO MDAP_RANDOMIZE, $
                inp, seed=seed
        sz=size(inp)
        inp=randomn(seed, sz[1:sz[0]], /double)
END

PRO MDAP_RANDOMIZE_FLAG, $
                inp, seed=seed
        sz=size(inp)
        out=randomn(seed, sz[1:sz[0]])
        inp = intarr(sz[1:sz[0]])
        indx = where(out gt 0, count)
        if count ne 0 then $
            inp[indx] = 1

END


pro test_write

        configure_file='mdap_march_config_v3.dat'
        READCOL, configure_file, command_line, comment='#', delimiter='%', /silent, format='A'
        for i=0, n_elements(command_line)-1 do $
            d=execute(command_line[i])

        READCOL, total_filelist, root_name_vector, velocity_initial_guess_vector,$
                 velocity_dispersion_initial_guess_vector, ellipticity_vector,$
                 position_angle_vector, fibre_number_vector, Reff_vector, mode_vector,$
                 /silent, format='A,F,F,F,F,F,F,A', comment='#'

        input_number=0
        root_name = root_name_vector[input_number]
        mode = mode_vector[input_number]
        velocity_initial_guess = velocity_initial_guess_vector[input_number]
        velocity_dispersion_initial_guess = velocity_dispersion_initial_guess_vector[input_number]
        ell=ellipticity_vector[input_number]
        pa=position_angle_vector[input_number]
        number_of_fibres=fibre_number_vector[input_number]
        Reff=Reff_vector[input_number]

        MDAP_SETUP_IO, root_name, output_root_dir, datacube_name, file_root, output_dir, $
                       output_file_root

        MDAP_DRP_CHECK_FILETYPE, datacube_name, mode

        ; TODO: TEMPORARY --------------------------------------------------------------

        tpl_library_keys = [ 'M11-MARCS' ]
        template_libraries = [ stellar_library_spatial_binning_1 ]
        emission_line_parameters = [ emission_line_file_spatial_binnin_1 ]
        absorption_line_parameters = [ absorption_line_indices ]

        bin_type = [ 'STON', 'STON', 'STON' ]

        if n_elements(w_range_for_sn_computation_for_gas) eq 0 then begin
            w_range_sn = transpose([ [w_range_for_sn_computation], [w_range_for_sn_computation], $
                                     [w_range_for_sn_computation] ])
        endif else begin
            w_range_sn = transpose([ [w_range_for_sn_computation], [w_range_for_sn_computation], $
                                     [w_range_for_sn_computation_for_gas] ])
        endelse

        bin_par = mode eq 'RSS' ? [ sn1_rss, sn2_rss, sn3_rss ] : $
                                   [ sn1_datacubes, sn2_datacubes, sn3_datacubes ]

        threshold_ston_bin = mode eq 'RSS' ? [ sn_thr_tpl_rss, sn_thr_str_rss, sn_thr_ems_rss ] : $
                                             [ sn_thr_tpl_datacubes, sn_thr_str_datacubes, $
                                               sn_thr_ems_datacubes ]

        if n_elements(weight_for_sn) eq 1 then begin
            bin_weight_by_sn2 = [ 1, 1, 1 ]
        endif else $
            bin_weight_by_sn2 = [ 0, 0, 0 ]

        w_range_analysis = transpose([ [trim_wav_range_spatial_binning_1], $
                                       [trim_wav_range_spatial_binning_2], $
                                       [trim_wav_range_spatial_binning_3] ])


        threshold_ston_analysis = [ 0.0d, 0.0d, 0.0d ]

        analysis_per_bin=[ 'stellar-cont', 'emission-line', 'abs-indices' ]
        analysis = transpose([ [analysis_per_bin], [analysis_per_bin], [analysis_per_bin] ])

        max_el = 0
        if (max_el lt n_elements(spectra_fittin_parameters_patial_binning_1)) then $
            max_el = n_elements(spectra_fittin_parameters_patial_binning_1)
        if (max_el lt n_elements(spectra_fittin_parameters_patial_binning_2)) then $
            max_el = n_elements(spectra_fittin_parameters_patial_binning_2)
        if (max_el lt n_elements(spectra_fittin_parameters_patial_binning_3)) then $
            max_el = n_elements(spectra_fittin_parameters_patial_binning_3)

        analysis_extra = strarr(3, max_el)
        analysis_extra[0,0:n_elements(spectra_fittin_parameters_patial_binning_1)-1] = $
                spectra_fittin_parameters_patial_binning_1
        analysis_extra[1,0:n_elements(spectra_fittin_parameters_patial_binning_2)-1] = $
                spectra_fittin_parameters_patial_binning_2
        analysis_extra[2,0:n_elements(spectra_fittin_parameters_patial_binning_3)-1] = $
                spectra_fittin_parameters_patial_binning_3

        tpl_lib_analysis = [ 0, 0, 0 ]
        ems_par_analysis = [ 0, 0, 0 ]
        abs_par_analysis = [ 0, 0, 0 ]

        overwrite_flag = [ 1, 1, 1 ]

        ;-------------------------------------------------------------------------------

        n_tpl_lib = n_elements(template_libraries)      ; Number of template libraries to use
        n_ems_par = n_elements(emission_line_parameters); Number of emission line parameter files
        n_abs_par = n_elements(absorption_line_parameters);Number of absorption line parameter files

        MDAP_BUILD_EXECUTION_PLANS, n_tpl_lib, n_ems_par, n_abs_par, bin_type, w_range_sn, $
                                    bin_par, threshold_ston_bin, bin_weight_by_sn2, $
                                    w_range_analysis, threshold_ston_analysis, analysis, $
                                    analysis_extra, tpl_lib_analysis, ems_par_analysis, $
                                    abs_par_analysis, overwrite_flag, execution_plan

        MDAP_GENERATE_OUTPUT_FILE_NAMES, output_file_root, execution_plan

        mdap_read_drp_fits_version='0'
        mdap_calculate_sn_version='0'
        mdap_spatial_binning_version='0'

        MDAP_READ_DRP_FITS, version=mdap_read_drp_fits_version
        mdap_fiducial_bin_xy_version = '0.1'
        mdap_tpl_lib_setup_version = '0.1'
        MDAP_CALCULATE_SN, version=mdap_calculate_sn_version
        MDAP_SPATIAL_BINNING, version=mdap_spatial_binning_version

        ;-------------------------------------------------------------------------------
        MDAP_READ_DRP_FITS, datacube_name, header, flux, ivar, mask, wave, sres, skyx, skyy, $
                            type=mode, unit=unit
        
        MDAP_GET_SPAXEL_SIZE, header, spaxel_dx, spaxel_dy, type=mode, unit=unit

        MDAP_FIDUCIAL_BIN_XY, skyx, skyy, bskyx, bskyy

        SXADDPAR, header, 'VDAPREAD', mdap_read_drp_fits_version, 'mdap_read_drp_fits version'
        SXADDPAR, header, 'VDAPFXY', mdap_fiducial_bin_xy_version, 'mdap_fiducial_bin_xy version'
        
        ;-------------------------------------------------------------------------------

        i = 0

        MDAP_SELECT_GOOD_SPECTRA, flux, ivar, mask, gflag, gindx, count=gcount

        ; Select the pixels to use in the S/N calculation
        MDAP_SELECT_WAVE, wave, execution_plan[i].wave_range_sn, lam_sn

        ; Calculate the S/N per pixel over some wavelength range
        MDAP_CALCULATE_SN, flux, ivar, mask, wave, lam_sn, signal, noise, gflag=gflag
        
        SXADDPAR, header, 'VDAPSTON', mdap_calculate_sn_version, 'mdap_calculate_sn version'

        print, 'BASIC DATA'
        MDAP_WRITE_OUTPUT, execution_plan[i].ofile, header=header, dx=spaxel_dx, dy=spaxel_dy, $
                           w_range_sn=execution_plan[i].wave_range_sn, xpos=bskyx, ypos=bskyy, $
                           signal=signal, noise=noise

        ;-------------------------------------------------------------------------------
        MDAP_SPATIAL_BINNING, flux, ivar, mask, signal, noise, gflag, bskyx, bskyy, spaxel_dx, $
                              spaxel_dy, execution_plan[i].bin_type, execution_plan[i].bin_par, $
                              execution_plan[i].threshold_ston_bin, $
                              execution_plan[i].bin_weight_by_sn2, bin_wgts, bin_indx, bin_flux, $
                              bin_ivar, bin_mask, xbin, ybin, bin_area, bin_ston, nbin, $
                              sn_calibration=sn_calibration;, /plot

        SXADDPAR, header, 'VDAPBIN', mdap_spatial_binning_version, 'mdap_spatial_binning version'

        ; Write the binning results
        print, 'BINNING RESULTS'
        MDAP_WRITE_OUTPUT, execution_plan[i].ofile, header=header, $
                           bin_type=execution_plan[i].bin_type, bin_par=execution_plan[i].bin_par, $
                           threshold_ston_bin=execution_plan[i].threshold_ston_bin, $
                           bin_indx=bin_indx, bin_weights=bin_wgts, wave=wave, sres=sres, $
                           bin_flux=bin_flux, bin_ivar=bin_ivar, bin_mask=bin_mask, xbin=xbin, $
                           ybin=ybin, bin_area=bin_area, bin_ston=bin_ston, bin_n=nbin, /read_header
        
        ;-------------------------------------------------------------------------------
        library_key = tpl_library_keys[execution_plan[i].tpl_lib]; 'M11-MARCS'
        tpl_out_fits = MDAP_SET_TPL_LIB_OUTPUT_FILE(output_file_root, library_key)
        
        tpl_flux=READFITS(tpl_out_fits, exten_no=0)
        tpl_wave=READFITS(tpl_out_fits, exten_no=1)
        tpl_ivar=READFITS(tpl_out_fits, exten_no=2)
        tpl_mask=READFITS(tpl_out_fits, exten_no=3)
        tpl_sres=READFITS(tpl_out_fits, exten_no=4)
        tpl_soff=READFITS(tpl_out_fits, exten_no=5)

        eml_par = MDAP_READ_EMISSION_LINE_PARAMETERS(emission_line_parameters[ $
                                                             execution_plan[i].ems_par ])

        star_kin_starting_guesses = dblarr(n_elements(xbin),4)
        star_kin_starting_guesses[*,0] = velocity_initial_guess[0]              ; stellar velocity
        star_kin_starting_guesses[*,1]=velocity_dispersion_initial_guess[0]     ; stellar sigma
        gas_kin_starting_guesses = star_kin_starting_guesses[*,0:1]             ; gas velocity
        gas_kin_starting_guesses[*,1]=50.                                       ; gas sigma
        
        ;-------------------------------------------------------------------------------

        ppxf_only=0

        MDAP_SPECTRAL_FITTING, wave, bin_flux, bin_ivar, bin_mask, sres, tpl_wave, tpl_flux, $
                               tpl_ivar, tpl_mask, wavelength_output, obj_fit_mask_ppxf, $
                               weights_ppxf, add_poly_coeff_ppxf, mult_poly_coeff_ppxf, $
                               bestfit_ppxf, chi2_ppxf, obj_fit_mask_gndf, weights_gndf, $
                               mult_poly_coeff_gndf, bestfit_gndf, chi2_gndf, eml_model, $
                               best_template, best_template_losvd_conv, stellar_kinematics, $
                               stellar_kinematics_err, emission_line_kinematics, $
                               emission_line_kinematics_err, emission_line_omitted, $
                               emission_line_kinematics_individual, $
                               emission_line_kinematics_individual_err, $
                               emission_line_intens, emission_line_intens_err, $
                               emission_line_fluxes, emission_line_fluxes_err, $
                               emission_line_EW, emission_line_EW_err, reddening_output, $
                               reddening_output_err, $
                               extra_inputs=execution_plan[i].analysis_extra, $
                               star_kin_starting_guesses=star_kin_starting_guesses, $
                               gas_kin_starting_guesses=gas_kin_starting_guesses, $
                               eml_par=eml_par, external_library=external_library, $
                               wave_range_analysis=execution_plan[i].wave_range_analysis, $
                               ppxf_only=ppxf_only ;, $ /quiet

        MDAP_RANDOMIZE, obj_fit_mask_ppxf, seed=seed
        MDAP_RANDOMIZE, weights_ppxf, seed=seed
        if n_elements(degree) ne 0 then begin
            if degree gt 0 then $
                MDAP_RANDOMIZE, add_poly_coeff_ppxf, seed=seed
        endif
        if n_elements(mdegree) ne 0 then begin
            if mdegree gt 0 then $
                MDAP_RANDOMIZE, mult_poly_coeff_ppxf, seed=seed
        endif
        MDAP_RANDOMIZE, bestfit_ppxf, seed=seed
        MDAP_RANDOMIZE, chi2_ppxf, seed=seed
        MDAP_RANDOMIZE, obj_fit_mask_gndf, seed=seed
        MDAP_RANDOMIZE, weights_gndf, seed=seed
;       if n_elements(reddening) eq 0 and n_elements(mdegree) ne 0 then begin
        if n_elements(reddening) eq 0 && n_elements(mdegree) ne 0 then begin
            if mdegree gt 0 then $
                MDAP_RANDOMIZE, mult_poly_coeff_gndf, seed=seed
        endif
        MDAP_RANDOMIZE, bestfit_gndf, seed=seed
        MDAP_RANDOMIZE, chi2_gndf, seed=seed
        MDAP_RANDOMIZE, eml_model, seed=seed
        MDAP_RANDOMIZE, best_template, seed=seed
        MDAP_RANDOMIZE, best_template_losvd_conv, seed=seed
        MDAP_RANDOMIZE, stellar_kinematics, seed=seed
        MDAP_RANDOMIZE, stellar_kinematics_err, seed=seed
        MDAP_RANDOMIZE, emission_line_kinematics, seed=seed
        MDAP_RANDOMIZE, emission_line_kinematics_err, seed=seed
        MDAP_RANDOMIZE_FLAG, emission_line_omitted, seed=seed
        MDAP_RANDOMIZE, emission_line_kinematics_individual, seed=seed
        MDAP_RANDOMIZE, emission_line_kinematics_individual_err, seed=seed
        MDAP_RANDOMIZE, emission_line_intens, seed=seed
        MDAP_RANDOMIZE, emission_line_intens_err, seed=seed
        MDAP_RANDOMIZE, emission_line_fluxes, seed=seed
        MDAP_RANDOMIZE, emission_line_fluxes_err, seed=seed
        MDAP_RANDOMIZE, emission_line_EW, seed=seed
        MDAP_RANDOMIZE, emission_line_EW_err, seed=seed
        MDAP_RANDOMIZE, reddening_output, seed=seed
        MDAP_RANDOMIZE, reddening_output_err, seed=seed

        print, 'ANALYSIS WAVE RANGE'
        MDAP_WRITE_OUTPUT, execution_plan[i].ofile, header=header, $
                           w_range_analysis=execution_plan[i].wave_range_analysis, /read_header

;       if (ppxf_only eq 1 and n_elements(eml_par) ne 0 ) or ppxf_only eq 0 then begin
        if (ppxf_only eq 1 && n_elements(eml_par) ne 0 ) || ppxf_only eq 0 then begin
            print, 'EMISSION LINE PARAMETERS'
            MDAP_WRITE_OUTPUT, execution_plan[i].ofile, eml_par=eml_par
            ;stop
        endif

        print, 'STELLAR KINEMATICS'
        MDAP_WRITE_OUTPUT, execution_plan[i].ofile, header=header, $
                           obj_fit_mask_ppxf=obj_fit_mask_ppxf, weights_ppxf=weights_ppxf, $
                           add_poly_coeff_ppxf=add_poly_coeff_ppxf, $
                           mult_poly_coeff_ppxf=mult_poly_coeff_ppxf, $
                           bestfit_ppxf=bestfit_ppxf, chi2_ppxf=chi2_ppxf, $
                           stellar_kinematics_fit=stellar_kinematics, $
                           stellar_kinematics_err=stellar_kinematics_err, $
                           extra_inputs=execution_plan[i].analysis_extra, /read_header
        
        if ppxf_only eq 0 then begin
            print, 'EMISSION-LINE FITS'
            MDAP_WRITE_OUTPUT, execution_plan[i].ofile, eml_par=eml_par, $
                               obj_fit_mask_gndf=obj_fit_mask_gndf, weights_gndf=weights_gndf, $
                               mult_poly_coeff_gndf=mult_poly_coeff_gndf, $
                               emission_line_kinematics_avg=emission_line_kinematics, $
                               emission_line_kinematics_aer=emission_line_kinematics_err, $
                               chi2_gndf=chi2_gndf, $
                               emission_line_kinematics_ind=emission_line_kinematics_individual, $
                               emission_line_kinematics_ier=emission_line_kinematics_individual_err, $
                               emission_line_omitted=emission_line_omitted, $
                               emission_line_intens=emission_line_intens, $
                               emission_line_interr=emission_line_intens_err, $
                               emission_line_fluxes=emission_line_fluxes, $
                               emission_line_flxerr=emission_line_fluxes_err, $
                               emission_line_EWidth=emission_line_EW, $
                               emission_line_EW_err=emission_line_EW_err, $
                               reddening_val=reddening, reddening_err=reddening_err, $
                               bestfit_gndf=bestfit_gndf, eml_model=eml_model
        endif

        print, 'OPTIMAL TEMPLATES'
        MDAP_WRITE_OUTPUT, execution_plan[i].ofile, optimal_template=best_template, $
                           losvd_optimal_template=best_template_losvd_conv 

END



