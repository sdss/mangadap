
PRO TEST_EML

        input_number = 0
        configure_file = 'mdap_march_config_v2.dat'
        
        ; Use configuration file to set some user-defined variables
        READCOL, configure_file, command_line, comment='#', delimiter='%', /silent, format='A'
        for i=0, n_elements(command_line)-1 do $
            d=execute(command_line[i])

        READCOL, total_filelist, root_name_vector, velocity_initial_guess_vector,$
                 velocity_dispersion_initial_guess_vector, ellipticity_vector,$
                 position_angle_vector, fibre_number_vector, Reff_vector, mode_vector,$
                 /silent, format='A,F,F,F,F,F,F,A', comment='#'

        ; Save the one requested by the user
        root_name = root_name_vector[input_number]
        velocity_initial_guess = velocity_initial_guess_vector[input_number]

        MDAP_SETUP_IO, root_name, output_root_dir, datacube_name, file_root, output_dir, $
                           output_file_root, output_filefits, output_idlsession


        MDAP_READ_DRP_FITS, datacube_name, header, flux, ivar, mask, wave, sres, skyx, skyy, $
                                type=mode, unit=unit

        eml_par_binning_1 = MDAP_READ_EMISSION_LINE_PARAMETERS(emission_line_file_spatial_binnin_1)
        MDAP_CHECK_EMISSION_LINES, eml_par_binning_1, wave, velocity=velocity_initial_guess

END

