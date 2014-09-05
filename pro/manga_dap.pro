;+
; NAME:
;	MANGA_DAP
;
; PURPOSE:
;	Main wrapping routine for MaNGA data analysis pipeline
;
; CALLING SEQUENCE:
;	MANGA_DAP, input_number, configure_file
;
; INPUTS:
;    input_number int
;		Index number (line number minus one) of the data to analyze from
;		total_filelist file defined in the configure_file.  
;
;    configure_file string
;		"System file" that defines a number of directories and files
;		used by the analysis pipeline.
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;
; OUTPUT:
;
; OPTIONAL OUTPUT:
;
; TODO:
;	- Use environmental variables instead of configure file?
;	- Need a better description of the configure_file; here or elsewhere?
;	- Change default names
;	- Figure out how restoring an old IDL session works
;		- need to think about how to restart these things
;	- Include MW extinction curve
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;
;	setup_dir_file
;
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;	01 Sept 2014: Copied from v0_8 version by L. Coccato.
;	02 Sept 2014: (KBW) Formating and minor edits
;	03 Sept 2014: (KBW) Basic compilation errors
;-
;------------------------------------------------------------------------------

pro MANGA_DAP,$
		input_number, configure_file

	manga_dap_version = '0.9'	; set the version number
	;stop

	;check_version=check_version,$
	;       dont_remove_null_templates=dont_remove_null_templates, $
	;	datacubes_or_rss=datacubes_or_rss

	c=299792.458d			; TODO: there isn't some idlutils file that defines this?
	t0=systime(/seconds)/60./60.	; start time

	;*****SETTING THE CONFIGURATION PARAMETERS DEFINED IN THE CONFIGURATION FILE
	READCOL, configure_file, command_line, comment='#', delimiter='%', /silent, format='A'

;	for i=0,n_elements(command_line)-1 do $
;	    print, command_line[i]
	
	for i=0, n_elements(command_line)-1 do $
	    d=execute(command_line[i])

	print, total_filelist
	print, output_root_dir
	print, stellar_library_spatial_binning_1
	print, stellar_library_spatial_binning_2
	print, stellar_library_spatial_binning_3
	print, emission_line_file_spatial_binnin_1
	print, emission_line_file_spatial_binnin_2
	print, emission_line_file_spatial_binnin_3
	print, emission_line_file_radial_binning
	print, absorption_line_indices
	print, external_library
	print, instrumental_fwhm_file

	; Restore an old IDL session if it exists (don't really know how this works)
	; TODO: find out what restore does and whether or not I can compress these lines
	; TODO: from Lodo: restore session can overwrite some variables so they
	;	need to be read again
	if keyword_set(check_version) then begin
	    print, 'Reading previous IDL session.'
	    READCOL, total_filelist, root_name_vector, velocity_initial_guess_vector,$
		     velocity_dispersion_initial_guess_vector, ellipticity_vector,$
		     position_angle_vector, fibre_number_vector, Reff_vector, mode_vector,$
		     /silent, format='A,F,F,F,F,F,F,A', comment='#'
	    root_name = root_name_vector[input_number]
	    mode = mode_vector[input_number]

	    MDAP_SETUP_IO, root_name, output_root_dir, datacube_name, file_root, output_dir, $
			   output_file_root, output_filefits, output_idlsession

	    output_idlsession=output_dir+'mdap_session.idl'
	    check_previous_session = file_test(output_idlsession)
	    if check_previous_session eq 1 then restore,output_idlsession
	endif

	; Read the table with the list of fits files and initial guess values
	print, 'Reading input table.'
	READCOL, total_filelist, root_name_vector, velocity_initial_guess_vector,$
	         velocity_dispersion_initial_guess_vector, ellipticity_vector,$
		 position_angle_vector, fibre_number_vector, Reff_vector, mode_vector,$
		 /silent, format='A,F,F,F,F,F,F,A', comment='#'

	; Save the one requested by the user
	root_name = root_name_vector[input_number]
	mode = mode_vector[input_number]
;	datacube_dir = datacube_root_dir+mode[0]+'/'
	velocity_initial_guess = velocity_initial_guess_vector[input_number]
	velocity_dispersion_initial_guess = velocity_dispersion_initial_guess_vector[input_number]
	ell=ellipticity_vector[input_number]
	pa=position_angle_vector[input_number]
	number_of_fibres=fibre_number_vector[input_number]
	Reff=Reff_vector[input_number]

	print, 'ROOT: '+root_name
	print, 'MODE: '+mode
;	print, 'DATACUBE_DIR: '+datacube_dir
	print, 'Guess V', velocity_initial_guess
	print, 'Guess VDISP', velocity_dispersion_initial_guess
	print, 'Guess ELL', ell
	print, 'Guess PA', pa
	print, 'IFU SIZE', number_of_fibres
	print, 'Half-light', Reff

	MDAP_SETUP_IO, root_name, output_root_dir, datacube_name, file_root, output_dir, $
		       output_file_root, output_filefits, output_idlsession

	; Set input file names and output directory
	;     make the directory if it doesn't exist
	;     mode must be either 'datacubes' or 'rss'

	; TODO: mode is input but mode is also determined in mdap_read_datacube, right?
	;       mode should also be apparent from the name returned by DRP
	;	... mode is used above to set the output directory
;	if mode eq 'datacubes' then begin		; mode must be 'datacubes' ...
;	    datacube_name=root_name+'.fits'
;	    res = file_test(output_root_dir+'results_datacubes/'+root_name+'/',/directory)
;	    if res eq 0 then spawn,'mkdir '+output_root_dir+'results_datacubes/'+root_name
;	endif else if mode eq 'rss' then begin		; ... or 'rss'
;	    datacube_name=root_name+'_rss.fits'
;	    res = file_test(output_root_dir+'results_rss/'+root_name+'/',/directory)
;	    if res eq 0 then spawn,'mkdir '+output_root_dir+'results_rss/'+root_name
;	endif else begin
;	    print,'Unrecognized data type'		; mode is non-sense
;	    return
;	endelse

	print, save_intermediate_steps
	print, remove_null_templates
	print, instrumental_fwhm_file

	;*****RE-SETTING THE CONFIGURATION PARAMETERS DEFINED IN THE CONFIGURATION FILE
	; why do these need to be re-read?
	READCOL, configure_file, command_line, comment='#', delimiter='%', /silent, format='A'
	for i=0, n_elements(command_line)-1 do $
	    d=execute(command_line[i])

	print, "RE-EXECUTED"
	print, save_intermediate_steps
	print, remove_null_templates
	print, instrumental_fwhm_file

	; TODO: This statement does nothing, right?
	if n_elements(save_intermediate_steps) eq 0 then save_intermediate_steps=0

	; This means that templates are ALWAY removed; so the value in configure_file does nothing
	if n_elements(remove_null_templates) eq 0 then remove_null_templates = 1

	; This file provides the spectral resolution
	; Columns are:
	;	1. Wavelength in angstroms
	;	2. Resolution (lamda/delta lambda)
	;	3. delta lamga (FWHM) in angstroms
	;	4. delta lamga (FWHM) in km/s
	if n_elements(instrumental_fwhm_file) ne 0 then $
	    READCOL, instrumental_fwhm_file, ww, r, fwhm_ang, fwhm_kms, /silent

	; Tell the user what's happening
	print, ''
	print, '# WORKING ON '+root_name+' ('+mode+')'
	print, ''
	; TODO: REOPEN
;	openw,1,output_dir+'mdap.log'

	print, 'BLOCK 0 ... '
	; BLOCK 0 ----------------------------------------------------------------------
	;	Description?
	;-------------------------------------------------------------------------------

	;*** GET MODULES VERSION AND CHECK PREVIOUS ANALYSIS ******************************
	check_previous_analysis = file_test(output_filefits)	; Check for output fits file
	check_previous_session = file_test(output_idlsession)	; Check for idl session file

	print, check_previous_analysis
	print, check_previous_session

	; if BOTH exist, read the fits file and determine the progress using the header
	if check_previous_analysis eq 1 and check_previous_session eq 1 $
	   and keyword_set(check_version) then begin
	    tmp_header = HEADFITS(output_filefits)		; Grab the header
	    manga_dap_version_previous = strcompress(sxpar(tmp_header,'DAPVER'),/remove_all)

	    mdap_read_datacube_version=strcompress(sxpar(tmp_header,'BLOCK1'),/remove_all)
	    mdap_read_datacube_version_previous = mdap_read_datacube_version

	    mdap_spatial_binning_version=strcompress(sxpar(tmp_header,'BLOCK2'),/remove_all)
	    mdap_spatial_binning_version_previous = mdap_spatial_binning_version

	    mdap_log_rebin_version=strcompress(sxpar(tmp_header,'BLOCK3'),/remove_all)
	    mdap_log_rebin_version_previous=mdap_log_rebin_version

	    mdap_spectral_fitting_version=strcompress(sxpar(tmp_header,'BLOCK4'),/remove_all)
	    mdap_spectral_fitting_version_previous=mdap_spectral_fitting_version

	    mdap_measure_indices_version=strcompress(sxpar(tmp_header,'BLOCK5'),/remove_all)
	    mdap_measure_indices_version_previous=mdap_measure_indices_version

	    mdap_spatial_radial_binning_version=strcompress(sxpar(tmp_header,'BLOCK6'),/remove_all)
	    mdap_spatial_radial_binning_version_previous=mdap_spatial_radial_binning_version
	endif else begin					; Set all versions to '0'
	    manga_dap_version_previous = '0'

	    mdap_read_datacube_version='0'
	    mdap_read_datacube_version_previous='0'

	    mdap_spatial_binning_version='0'
	    mdap_spatial_binning_version_previous='0'

	    mdap_log_rebin_version='0'
	    mdap_log_rebin_version_previous='0'

	    mdap_spectral_fitting_version='0'
	    mdap_spectral_fitting_version_previous='0'

	    mdap_measure_indices_version='0'
	    mdap_measure_indices_version_previous='0'

	    mdap_spatial_radial_binning_version='0'
	    mdap_spatial_radial_binning_version_previous='0'
	endelse

	print, manga_dap_version_previous

	print, mdap_read_datacube_version
	print, mdap_read_datacube_version_previous

	print, mdap_spatial_binning_version
	print, mdap_spatial_binning_version_previous

	print, mdap_log_rebin_version
	print, mdap_log_rebin_version_previous

	print, mdap_spectral_fitting_version
	print, mdap_spectral_fitting_version_previous

	print, mdap_measure_indices_version
	print, mdap_measure_indices_version_previous

	print, mdap_spatial_radial_binning_version
	print, mdap_spatial_radial_binning_version_previous

	; update the versions of each procedure to check if these procedures should be redone
	MDAP_READ_DATACUBE,version=mdap_read_datacube_version
	MDAP_SPATIAL_BINNING,version=mdap_spatial_binning_version
	MDAP_LOG_REBIN,version=mdap_log_rebin_version
	MDAP_SPECTRAL_FITTING,version=mdap_spectral_fitting_version
	MDAP_MEASURE_INDICES,version=mdap_measure_indices_version
	MDAP_SPATIAL_RADIAL_BINNING,version=mdap_spatial_radial_binning_version

	print, mdap_read_datacube_version
	print, mdap_spatial_binning_version
	print, mdap_log_rebin_version
	print, mdap_spectral_fitting_version
	print, mdap_measure_indices_version
	print, mdap_spatial_radial_binning_version

	; END BLOCK 0 ? ----------------------------------------------------------------


;; KINEMETRY CHECK, remove it after the tests
; restore, output_idlsession
; goto, kinemetry_step
;;

	; If the version of THIS procedure is more recent, redo all procedures
	execute_all_modules=0
	if manga_dap_version gt manga_dap_version_previous then execute_all_modules=1

	print, manga_dap_version
	print, manga_dap_version_previous

	print, execute_all_modules

	; Bin the data
	; Rebin the data to logarithmic binning
	; Perform the spectral fitting (emission lines?)
	; Measure the absorption-line indices
	; Radially bin the data

;goto, inizia_da_qui

	; BLOCK 1 ----------------------------------------------------------------------
	;	Description?
	;-------------------------------------------------------------------------------

	; Read the datacube (will convert from log to linear bin?)
	if mdap_read_datacube_version gt mdap_read_datacube_version_previous $
	   or execute_all_modules eq 1 then begin

	    MDAP_READ_DATACUBE, datacube_name, data, error, wavelength, x2d, y2d, signal, noise, $
				cdelt1, cdelt2, header_2d, lrange=w_range_for_sn_computation, $
				/keep_original_step, x2d_reconstructed=x2d_reconstructed, $
				y2d_reconstructed=y2d_reconstructed, $
				signal2d_reconstructed=signal2d_reconstructed, $
				noise2d_reconstructed=noise2d_reconstructed, $
				number_of_fibres=number_of_fibres

return
	    ;same as signal if datacube is read, not RSS.
	    if n_elements(signal2d_reconstructed) eq 0 then signal2d_reconstructed = signal

	    ;same as x2d if datacube is read, not RSS.
	    if n_elements(x2d_reconstructed) eq 0 then x2d_reconstructed = x2d
	    
	    ;same as y2d if datacube is read, not RSS.
	    if n_elements(y2d_reconstructed) eq 0 then y2d_reconstructed = y2d
	    
	    ;same as noise if datacube is read, not RSS.
	    if n_elements(noise2d_reconstructed) eq 0 then noise2d_reconstructed = noise

	    ; Add the version number for block 1
	    SXADDPAR, header_2d, 'BLOCK1', mdap_read_datacube_version, 'mdap_read_datacube version'
	    
	    execute_all_modules=1

	    ;-- I compute the signal and the noise over a different wavelength range (redshifted)
	    ; The mdap_read_datacube is run again, but only the signal and noise output are stored

	    ; Recalculate the S/N for for the gas; wasteful to rerun mdap_read_datacube?
	    if n_elements(w_range_for_sn_computation_for_gas) ne 0 then begin
		w_range_for_sn_computation_for_gas_redshifted = $
			w_range_for_sn_computation_for_gas*(1.+velocity_initial_guess[0]/c)
;			w_range_for_sn_computation_for_gas*(1.+velocity_initial_guess[0]/299792.5)

		MDAP_READ_DATACUBE, datacube_name, junk, junk, junk, junk, junk, signal_for_gas, $
				    noise_for_gas, junk, junk, header_junk, $
				    lrange=w_range_for_sn_computation_for_gas_redshifted, $
				    /keep_original_step, /use_total
		junk = 0.
	    endif else begin
		signal_for_gas = signal
		noise_for_gas = noise
		w_range_for_sn_computation_for_gas = w_range_for_sn_computation
	    endelse
	endif

	; Add the version number for manga_dap to the header
	SXADDPAR, header_2d, 'DAPVER', manga_dap_version, 'manga_dap version', BEFORE='BLOCK1'

	; Print progress
	printf, 1, '[INFO] mdap_read_datacube ver ' $
	            +max([mdap_read_datacube_version, mdap_read_datacube_version_previous])
	sz=size(data)
	if sz[0] eq 3 then begin
	    printf, 1, '[INFO] datacube '+root_name+' read; size: '+MDAP_STC(sz[1],/integer) $
		       +'x'+MDAP_STC(sz[2],/integer)
	endif else if sz[0] eq 2 then begin
	    printf, 1, '[INFO] RSS file '+root_name+' read; Nfibres: '+MDAP_STC(sz[1],/integer)
	endif else begin
	    printf, 1, '[INFO] Unrecognized file size'
	    return
	endelse

	; Initialize the extinction curve
	MW_extinction = 0.			; TODO: replaced with appropriate extension in the
						;       input file, or total_filelist

	; Initialize the instrumental resolution.  TODO: use fits extension
	if n_elements(instrumental_fwhm_file) ne 0 then begin	; use input file
	    fwhm_instr=interpol(fwhm_ang,ww,wavelength)
	    fwhm_instr_kmsec_matrix = fltarr(2,n_elements(ww))
	    fwhm_instr_kmsec_matrix[0,*]=ww
	    fwhm_instr_kmsec_matrix[1,*]=fwhm_kms
	endif else begin					; set to a constant
	    fwhm_instr=wavelength*0.+2.54			; TODO: Don't like this practice
	endelse

	; Set mask.  TODO: use fits extension
	mask_range=[5570.,5590.,5800.,5850.]

	; END BLOCK 1 ------------------------------------------------------------------

	; TODO: Why not SXADDPAR here?
	; Add the initial guesses to the header
	header_2d = [ header_2d[0:n_elements(header_2d)-2], $
		      'Reff = '+MDAP_STC(Reff)+'/ effective radius (arcsec)', $
		      'PA = '+MDAP_STC(PA)+'/ postion angle (0=north; 90 = east)', $
		      'ellipticity = '+MDAP_STC(ell) +'/ ellipticity (1-b/a)' ]
	; ... and the command-line input from the configuration file
	for i = 0,n_elements(command_line)-1 do $
	    header_2d=[header_2d,command_line[i]]

	; BLOCK 2 ----------------------------------------------------------------------
	;	Description?
	;-------------------------------------------------------------------------------

	; TODO: Don't really understand what's going on here
	;	- Can I make the different SN levels part of an array and then
	;	  loop over the array?
	;	- What's the difference between sn? and sn_thr_??? 
	;	- In the configuration file, the sn_thr_??? values are all the
	;	  same but different for rss and datacubes

	; Spatially rebin the data at 6 different S/N levels.  These levels are
	; defined in the configuration_file

	if mode eq 'rss' then begin
	    sn1=sn1_rss
	    sn2=sn2_rss
	    sn3=sn3_rss
	    sn_thr_tpl=sn_thr_tpl_rss
	    sn_thr_str=sn_thr_str_rss
	    sn_thr_ems=sn_thr_ems_rss
	    if n_elements(sn_calibration_rss) ne 0 then sn_calibration = sn_calibration_rss
	endif else if mode eq 'datacubes' then begin
	    sn1=sn1_datacubes
	    sn2=sn2_datacubes
	    sn3=sn3_datacubes
	    sn_thr_tpl=sn_thr_tpl_datacubes
	    sn_thr_str=sn_thr_str_datacubes
	    sn_thr_ems=sn_thr_ems_datacubes
	    if n_elements(sn_calibration_datacubes) ne 0 then begin
		sn_calibration = sn_calibration_datacubes
	    endif
	endif else begin
	     print,'Unrecognized data type.'
	     return
	endelse

	; Only execute the binning if required
	if mdap_spatial_binning_version gt mdap_spatial_binning_version_previous $
	   or execute_all_modules eq 1 then begin

	    ; Perform the spatial binning
	    MDAP_SPATIAL_BINNING, data, error, signal, noise, sn1, x2d, y2d, cdelt1, cdelt2, $
				  spatial_binning_tpl, spectra_tpl, errors_tpl, xbin_tpl, $
				  ybin_tpl, area_bins_tpl, bin_sn_tpl, sn_thr=sn_thr_tpl, $
				  x2d_reconstructed=x2d_reconstructed, $
				  y2d_reconstructed=y2d_reconstructed, $
				  nelements_within_bin=nelements_within_bin_tpl, $
				  sn_calibration=sn_calibration, $
				  user_bin_map=user_bin_map_spatial_binning_1, $
				  weight_for_sn=weight_for_sn

	    ; Add the version of the spatial binning procedure to the header
	    SXADDPAR, header_2d, 'BLOCK2', mdap_spatial_binning_version, $
		      'mdap_spatial_binning version'

	    ; TODO: This will force the binning at other S/N levels to be redone
	    ;       as well?  Does one depend on the other?

	    ; Re-execute all the remaining modules, if that wasn't already set
	    execute_all_modules=1
	endif

	; Print progress
	printf, 1, '[INFO] mdap_spatial_binning ver ' $
		   +max([mdap_spatial_binning_version, mdap_spatial_binning_version_previous])
	printf, 1, '[INFO] datacube '+root_name+' spatial binning 1. SN= ' $
		   +MDAP_STC(sn1,/integer)+' Nbins: '+MDAP_STC(n_elements(xbin_tpl), /integer)

	; TODO: Should be able to loop three times, instead of having 3 repeats
	;       of essentially the same code

	; Run the 2nd S/N level
	if mdap_spatial_binning_version gt mdap_spatial_binning_version_previous $
	   or execute_all_modules eq 1 then begin

	    ; Perform the spatial binning
	    MDAP_SPATIAL_BINNING, data, error, signal, noise, sn2, x2d, y2d, cdelt1, cdelt2, $
				  spatial_binning_str, spectra_str,errors_str, xbin_str, $
				  ybin_str, area_bins_str, bin_sn_str, sn_thr=sn_thr_str, $
				  x2d_reconstructed=x2d_reconstructed, $
				  y2d_reconstructed=y2d_reconstructed, $
				  nelements_within_bin=nelements_within_bin_str, $
				  sn_calibration=sn_calibration, $
				  user_bin_map=user_bin_map_spatial_binning_2, $
				  weight_for_sn=weight_for_sn

	    ; TODO: why no information added to the header?

	    ; Re-execute all the remaining modules, if that wasn't already set
	    execute_all_modules=1
	endif

	; Print progress
	printf, 1, '[INFO] datacube '+root_name+' spatial binning 2. SN= ' $
		   + MDAP_STC(sn2,/integer)+' Nbins: '+MDAP_STC(n_elements(xbin_str), /integer)

	; Run the final S/N level
	if mdap_spatial_binning_version gt mdap_spatial_binning_version_previous $
	   or execute_all_modules eq 1 then begin

	    ; Perform the spatial binning
	    MDAP_SPATIAL_BINNING, data, error, signal_for_gas, noise_for_gas, sn3, x2d, y2d, $
				  cdelt1, cdelt2, spatial_binning_ems, spectra_ems, errors_ems, $
				  xbin_ems, ybin_ems, area_bins_ems, bin_sn_ems, $
				  sn_thr=sn_thr_ems, x2d_reconstructed=x2d_reconstructed, $
				  y2d_reconstructed=y2d_reconstructed, $
				  nelements_within_bin=nelements_within_bin_ems, $
				  sn_calibration=sn_calibration, $
				  user_bin_map=user_bin_map_spatial_binning_3, $
				  weight_for_sn=weight_for_sn

	    ; Re-execute all the remaining modules, if that wasn't already set
   	    execute_all_modules=1

	endif

	; Print the progress
	printf, 1, '[INFO] datacube '+root_name+' spatial binning 3. SN= ' $
	           +MDAP_STC(sn3,/integer)+' Nbins: '+MDAP_STC(n_elements(xbin_ems), /integer)

	; END BLOCK 2 ------------------------------------------------------------------


	; BLOCK 3 ----------------------------------------------------------------------
	;	Description?
	;	
	;	- Determine the difference in instrumental resolution beween the
	;	  stellar templates and the data
	;	- Rebin the spatially binned spectra to a linear step in
	;	  ln(wavelength)
	;	- Select the wavelength range (TODO: for what, and is this
	;	  actually done?)
	;-------------------------------------------------------------------------------

	fwhm_stars = wavelength*0.+2.54			; Spectral resolution of MARCS templates

	fwhm_diff = sqrt(fwhm_instr^2 - fwhm_stars^2)	; Quadrature difference between spectral
							; resolution of input and of stellar
							; templates
							; TODO: Assumes Gaussian line-spread func

	; Check where the template instrumental resolution is larger than the data resolution
	indici = where(finite(fwhm_diff) ne 1)
	if indici[0] ne -1 then fwhm_diff[indici] = 0	; TODO: Forces them to be the same

	; Only execute the rebinning if required; 1st spatial binning scheme
	if mdap_log_rebin_version gt mdap_log_rebin_version_previous $
	   or execute_all_modules eq 1 then begin

	    ; Rebin the spectra to a linear step in ln(wavelength)
	    MDAP_LOG_REBIN, spectra_tpl, errors_tpl, wavelength, $
			    stellar_library_spatial_binning_1, fwhm_diff, log_spc_tpl, $
			    log_err_tpl, log_wav_tpl, library_log_tpl, log_wav_library_tpl, $
			    input_velscale=velscale, wave_range=trim_wav_range_spatial_binning_1, $
			    /gal_wavelength_log_step, /quiet

	    ; Add the version to the header
	    SXADDPAR, header_2d, 'BLOCK3', mdap_log_rebin_version, 'mdap_log_rebin version'

	    ; Re-execute all the remaining modules, if that wasn't already set
	    execute_all_modules=1
	endif

	; Print progress
	printf, 1, '[INFO] mdap_log_rebin ver ' $
		   +max([mdap_log_rebin_version, mdap_log_rebin_version_previous])
	printf, 1, '[INFO] datacube '+root_name+'log_rebin 1'

	; TODO: Should be able to do this with a loop?
	
	if mdap_log_rebin_version gt mdap_log_rebin_version_previous $
	   or execute_all_modules eq 1 then begin

	    ; Rebin the spectra to a linear step in ln(wavelength)
	    MDAP_LOG_REBIN, spectra_str, errors_str, wavelength, $
			    stellar_library_spatial_binning_2, fwhm_diff, log_spc_str, $
			    log_err_str, log_wav_str, library_log_str, log_wav_library_str, $
			    input_velscale=velscale, wave_range=trim_wav_range_spatial_binning_2, $
			    /gal_wavelength_log_step, /quiet

	    ; Re-execute all the remaining modules, if that wasn't already set
	    execute_all_modules=1
	endif

	; Print progress
	printf, 1, '[INFO] datacube '+root_name+'log_rebin 2'

	if mdap_log_rebin_version gt mdap_log_rebin_version_previous $
	   or execute_all_modules eq 1 then begin

	    ; Rebin the spectra to a linear step in ln(wavelength)
	    MDAP_LOG_REBIN, spectra_ems, errors_ems, wavelength, $
			    stellar_library_spatial_binning_3, fwhm_diff, log_spc_ems, $
			    log_err_ems, log_wav_ems, library_log_ems, log_wav_library_ems, $
			    input_velscale=velscale, wave_range=trim_wav_range_spatial_binning_3, $
			    /gal_wavelength_log_step, /quiet

	    ; Re-execute all the remaining modules, if that wasn't already set
	    execute_all_modules=1
	endif

	; Print progress
	printf,1,'[INFO] datacube '+root_name+'log_rebin 3'

	; TODO: Does the rebinning change to a constant step in log10(wave) or log(wave)?
	;       If the former, this equation is wrong.
	; TODO: DPR LOGCUBE are in log10(wave)!
	; Set the velocity scale of each pixel
	if n_elements(velscale) eq 0 then begin
	    velscale = c*(log_wav_tpl[1]-log_wav_tpl[0])
	endif

	; Set to a very small number
	; TODO: Why?
	library_log_tpl = library_log_tpl/10.^25.
	library_log_str = library_log_str/10.^25.
	library_log_ems = library_log_ems/10.^25.

	; TODO: First two blocks are not saved, but this one is.  Why?
	if save_intermediate_steps eq 1 then begin
	    save, filename=root_name+mode+'_block3.idl', /variables
	endif

	; END BLOCK 3 ------------------------------------------------------------------


	; ##############################################################################
	; BEGIN computation of "model-independent" data products #######################


	; BLOCK 4 ----------------------------------------------------------------------
	;	Description?
	;	- Run the spectral fitting procedures
	;-------------------------------------------------------------------------------

	; Run a full spectral fit including stellar population, gas and stellar
	; kinematics

	; Set the starting guesses
	; TODO: get the guesses from a cross-correlation with a guess template?
	star_kin_starting_guesses = fltarr(n_elements(xbin_tpl),4)
	star_kin_starting_guesses[*,0] = velocity_initial_guess[0]		; stellar velocity
	star_kin_starting_guesses[*,1]=velocity_dispersion_initial_guess[0]	; stellar sigma
	gas_kin_starting_guesses = star_kin_starting_guesses[*,0:1]		; gas velocity
	gas_kin_starting_guesses[*,1]=50.					; gas sigma

	; Only execute the fit if required
	if mdap_spectral_fitting_version gt mdap_spectral_fitting_version_previous $
	   or execute_all_modules eq 1 then begin

	    ; Perform the spectral fit

	    ; TODO: Be specific about what this execution of MDAP_SPECTRAL_FITTING does;
	    ;       I assume it is run multiple times using the different S/N limits for
	    ;       different purposes
	    MDAP_SPECTRAL_FITTING, log_spc_tpl, log_err_tpl, log_wav_tpl, library_log_tpl, $
				   log_wav_library_tpl, velscale, stellar_kinematics_tpl, $
				   stellar_kinematics_tpl_err, stellar_weights_tpl, $
				   emission_line_kinematics_tpl, $
				   emission_line_kinematics_tpl_err, $
				   emission_line_kinematics_tpl_individual, $
				   emission_line_kinematics_tpl_individual_err,$
				   emission_line_intens_tpl, emission_line_intens_tpl_err, $
				   emission_line_fluxes_tpl, emission_line_fluxes_tpl_err, $
				   emission_line_equivW_tpl, emission_line_equivW_tpl_err, $
				   wavelength_input=exp(log_wav_library_tpl), $
				   wavelength_output_tpl, best_fit_model_tpl, $
				   galaxy_minus_ems_fit_model_tpl, best_template_tpl, $
				   best_template_LOSVD_conv_tpl, reddening_tpl, $
				   reddening_tpl_err, residuals_tpl, $
				   star_kin_starting_guesses=star_kin_starting_guesses, $
				   gas_kin_starting_guesses=gas_kin_starting_guesses, $
				   emission_line_file=emission_line_file_binning_1, $
				   fwhm_instr_kmsec_matrix=fwhm_instr_kmsec_matrix, $
				   extra_inputs=spectral_fit_par_bin1, mask_range=mask_range, $
				   external_library=external_library, /quiet

	   ; need to save all flux maps (warning: sgandalf computes intensities, not fluxes)

	   ; Convert intensities and fluxes to surface brightness
	   junk = size(emission_line_fluxes_tpl)
	   for i = 0, junk[1]-1 do begin
		emission_line_intens_tpl[i,*] = emission_line_intens_tpl[i,*] / area_bins_tpl[i]
		emission_line_intens_tpl_err[i,*] = emission_line_intens_tpl_err[i,*] $
						    / sqrt(area_bins_tpl[i])
		emission_line_fluxes_tpl[i,*] = emission_line_fluxes_tpl[i,*] / area_bins_tpl[i]
		emission_line_fluxes_tpl_err[i,*] = emission_line_fluxes_tpl_err[i,*] $
						    / sqrt(area_bins_tpl[i])
	   endfor

;   for i = 0, junk[1]-1 do emission_line_intens_tpl[i,*]=emission_line_intens_tpl[i,*]/area_bins_tpl[i]
;   for i = 0, junk[1]-1 do emission_line_intens_tpl_err[i,*]=emission_line_intens_tpl_err[i,*]/sqrt(area_bins_tpl[i])
;   for i = 0, junk[1]-1 do emission_line_fluxes_tpl[i,*]=emission_line_fluxes_tpl[i,*]/area_bins_tpl[i]
;   for i = 0, junk[1]-1 do emission_line_fluxes_tpl_err[i,*]=emission_line_fluxes_tpl_err[i,*]/sqrt(area_bins_tpl[i])

	    ; Only keep templates with positive weights for the next fits
;	    if remove_null_templates[0] eq 1 then begin
	    if remove_null_templates eq 1 then begin
		wp = total(stellar_weights_tpl,1)
		library_log_tpl = temporary(library_log_tpl[*,where(wp gt 0)])
		library_log_str = temporary(library_log_str[*,where(wp gt 0)])
		library_log_ems = temporary(library_log_ems[*,where(wp gt 0)])
		stellar_weights_tpl = temporary(stellar_weights_tpl[*,where(wp gt 0)])
	    endif

	    ;LC: calculate the real S/N of binned spectra
	    ; TODO: Is this really the "real S/N"?
	    bin_sn_tpl_real = bin_sn_tpl
	    indxx = where(wavelength_output_tpl ge w_range_for_sn_computation[0] $
			  and wavelength_output_tpl le w_range_for_sn_computation[1])

	    for i = 0, n_elements(bin_sn_tpl_real) -1 do begin
		MDAP_CALCULATE_SPECTRUM_SN, best_fit_model_tpl[i,indxx], residuals_tpl[i,indxx], $
					    wavelength_output_tpl[indxx], sn_per_angstorm, /rms

		bin_sn_tpl_real[i] = sn_per_angstorm[0]
	    endfor

	    ; Add the spectral fitting version to the header
	    SXADDPAR, header_2d, 'BLOCK4', mdap_spectral_fitting_version, $
		      'mdap_spectral_fitting_version'

	    ; Re-execute all the remaining modules, if that wasn't already set
	    execute_all_modules=1
	endif

	; Print progress
	printf, 1, '[INFO] mdap_spectral_fitting ver '+ $
                   max([mdap_spectral_fitting_version, mdap_spectral_fitting_version_previous])
	printf,1,'[INFO] datacube '+root_name+' spectral fitting 1'

	; Save the progress
	if save_intermediate_steps eq 1 then begin
	    save, filename=root_name+mode+'_block4a.idl', /variables 
	endif

	; Use MDAP_CREATE_STARTING_GUESSES to update the guesses for the stellar kinematics
	; Interpolates the stellar kinematics results over ems grid to get starting guesses based
	;   on the previous results
	MDAP_CREATE_STARTING_GUESSES, stellar_kinematics_tpl, xbin_tpl, ybin_tpl, x2d, y2d, $
				      xbin_str, ybin_str, star_kin_starting_guesses, $
				      velocity_initial_guess[0], $
				      velocity_dispersion_initial_guess[0], 0., 0., /h3h4

	; ... and for gas kinematics
	MDAP_CREATE_STARTING_GUESSES, emission_line_kinematics_tpl, xbin_tpl, ybin_tpl, x2d, y2d, $
				      xbin_str, ybin_str, gas_kin_starting_guesses, $
				      velocity_initial_guess[0], 50., 0., 0.


	; Run the spectral fitting again
	if mdap_spectral_fitting_version gt mdap_spectral_fitting_version_previous $
	   or execute_all_modules eq 1 then begin

	    ; Perform an extra fit that only fits V and sigma, with h3=h4=0, if requested
	    if n_elements(spectra_fittin_parameters_patial_binning_2_two_moments) ne 0 then begin
	    
		MDAP_SPECTRAL_FITTING, log_spc_str, log_err_str, log_wav_str, library_log_str, $
				       log_wav_library_str, velscale, $
				       stellar_kinematics_str_two_moments, $
				       stellar_kinematics_str_two_moments_err, $
				       stellar_weights_str, emission_line_kinematics_str, $
				       emission_line_kinematics_str_err, $
				       emission_line_kinematics_str_individual, $
				       emission_line_kinematics_str_individual_err, $
				       emission_line_intens_str, emission_line_intens_str_err, $
				       emission_line_fluxes_str, emission_line_fluxes_str_err, $
				       emission_line_equivW_str, emission_line_equivW_str_err, $
				       wavelength_input=exp(log_wav_library_str), $
				       wavelength_output_str, best_fit_model_str, $
				       galaxy_minus_ems_fit_model_str, best_template_str, $
				       best_template_LOSVD_conv_str, reddening_str, $
				       reddening_str_err, residuals_str, $
				       star_kin_starting_guesses=star_kin_starting_guesses, $
				       gas_kin_starting_guesses=gas_kin_starting_guesses, $
				       fwhm_instr_kmsec_matrix=fwhm_instr_kmsec_matrix, $
				       emission_line_file=emission_line_file_binning_2, $
				       extra_inputs=spectral_fit_par_bin2_2mns, $
				       mask_range=mask_range, $
				       external_library=external_library, /quiet
	    endif

	    ; Fit the stellar kinematics with h3 and h4 as free variables
	    MDAP_SPECTRAL_FITTING, log_spc_str, log_err_str, log_wav_str, library_log_str, $
				   log_wav_library_str, velscale, stellar_kinematics_str, $
				   stellar_kinematics_str_err, stellar_weights_str, $
				   emission_line_kinematics_str, emission_line_kinematics_str_err, $
				   emission_line_kinematics_str_individual, $
				   emission_line_kinematics_str_individual_err, $
				   emission_line_intens_str, emission_line_intens_str_err, $
				   emission_line_fluxes_str, emission_line_fluxes_str_err, $
				   emission_line_equivW_str, emission_line_equivW_str_err, $
				   wavelength_input=exp(log_wav_library_str), $
				   wavelength_output_str, best_fit_model_str, $
				   galaxy_minus_ems_fit_model_str, best_template_str, $
				   best_template_LOSVD_conv_str, reddening_str, $
				   reddening_str_err, residuals_str, $
				   star_kin_starting_guesses=star_kin_starting_guesses, $
				   gas_kin_starting_guesses=gas_kin_starting_guesses, $
				   fwhm_instr_kmsec_matrix=fwhm_instr_kmsec_matrix, $
				   emission_line_file=emission_line_file_binning_2, $
				   extra_inputs=spectral_fit_par_bin2, $
				   mask_range=mask_range, external_library=external_library, /quiet

	    ; save all flux maps, converting to surface brightness
	    ; TODO: Is this really surface brightness?
	    junk = size(emission_line_fluxes_str)
	    for i = 0, junk[1]-1 do begin
		emission_line_intens_str[i,*] = emission_line_intens_str[i,*] / area_bins_str[i]
		emission_line_intens_str_err[i,*] = emission_line_intens_str_err[i,*] $
						    / sqrt(area_bins_str[i])
		emission_line_fluxes_str[i,*] = emission_line_fluxes_str[i,*] / area_bins_str[i]
		emission_line_fluxes_str_err[i,*] = emission_line_fluxes_str_err[i,*] $
						    / sqrt(area_bins_str[i])
	    endfor

	    ;LC: calculate the real S/N of binned spectra
	    ; TODO: Is this really the "real S/N"?
	    bin_sn_str_real = bin_sn_str
	    indxx = where(wavelength_output_str ge w_range_for_sn_computation[0] $
			  and wavelength_output_str le w_range_for_sn_computation[1] )
	    for i = 0, n_elements(bin_sn_str_real) -1 do begin
		MDAP_CALCULATE_SPECTRUM_SN, best_fit_model_str[i,indxx], residuals_str[i,indxx], $
					    wavelength_output_str[indxx], sn_per_angstorm, /rms
		bin_sn_str_real[i] = sn_per_angstorm[0]
	    endfor

	    ; Re-execute all the remaining modules, if that wasn't already set
	    execute_all_modules = 1

	endif

	; Print progress
	printf,1,'[INFO] datacube '+root_name+' spectral fitting 2'

	if save_intermediate_steps eq 1 then begin
	    save, filename=root_name+mode+'_block4b.idl', /variables 
	endif

	; Use MDAP_CREATE_STARTING_GUESSES to update the guesses for the stellar kinematics
	; Interpolates the stellar kinematics results over ems grid to get starting guesses based
	;   on the previous results
	MDAP_CREATE_STARTING_GUESSES, stellar_kinematics_str, xbin_str, ybin_str, x2d, y2d, $
				      xbin_ems, ybin_ems, star_kin_starting_guesses, $
				      velocity_initial_guess[0], $
				      velocity_dispersion_initial_guess[0], 0., 0., /h3h4

	; ... and for gas kinematics
	MDAP_CREATE_STARTING_GUESSES, emission_line_kinematics_str, xbin_str, ybin_str, x2d, y2d, $
				      xbin_ems, ybin_ems, gas_kin_starting_guesses, $
				      velocity_initial_guess[0], 50., 0., 0.

	; Convert from log wavelength to linear values
	; TODO: So log_wav_ems is in natural log, assume that jives with above?
	wavelength_input=exp(log_wav_ems)

	; Now run the spectral fitting again to get the as kinematics
	if mdap_spectral_fitting_version gt mdap_spectral_fitting_version_previous $
	   or execute_all_modules eq 1 then begin
	   
	   MDAP_SPECTRAL_FITTING, log_spc_ems, log_err_ems, log_wav_ems, library_log_ems, $
				   log_wav_library_ems, velscale, stellar_kinematics_ems, $
				   stellar_kinematics_ems_err, stellar_weights_ems, $
				   emission_line_kinematics_ems, emission_line_kinematics_ems_err, $
				   emission_line_kinematics_ems_individual, $
				   emission_line_kinematics_ems_individual_err, $
				   emission_line_intens_ems, emission_line_intens_ems_err, $
				   emission_line_fluxes_ems, emission_line_fluxes_ems_err, $
				   emission_line_equivW_ems, emission_line_equivW_ems_err, $
				   wavelength_input=wavelength_input, $
				   wavelength_output_rest_frame_log, best_fit_model_ems, $
				   galaxy_minus_ems_fit_model_ems, best_template_ems, $
				   best_template_LOSVD_conv_ems, reddening_ems, $
				   reddening_ems_err, residuals_ems, $
				   star_kin_starting_guesses=star_kin_starting_guesses, $
				   gas_kin_starting_guesses=gas_kin_starting_guesses, $
				   emission_line_file=emission_line_file_binning_3, $
				   extra_inputs=spectral_fit_par_bin3, /rest_frame_log, $
				   fwhm_instr_kmsec_matrix=fwhm_instr_kmsec_matrix, $
				   mask_range=mask_range, external_library=external_library, /quiet

	    ; save all flux maps, converting to surface brightness
	    junk = size(emission_line_fluxes_ems)
	    for i = 0, junk[1]-1 do begin
		emission_line_intens_ems[i,*] = emission_line_intens_ems[i,*] / area_bins_ems[i]
		emission_line_intens_ems_err[i,*] = emission_line_intens_ems_err[i,*] $
						    / sqrt(area_bins_ems[i])
		emission_line_fluxes_ems[i,*] = emission_line_fluxes_ems[i,*] / area_bins_ems[i]
		emission_line_fluxes_ems_err[i,*] = emission_line_fluxes_ems_err[i,*] $
						    / sqrt(area_bins_ems[i])
	    endfor

	    ; Re-execute all the remaining modules, if that wasn't already set
	    execute_all_modules = 1

	    ;LC: calculate the real S/N of binned spectra
	    ; TODO: Is this really the "real S/N"?
	    bin_sn_ems_real = bin_sn_ems
	    indxx = where(wavelength_output_rest_frame_log ge $
			  w_range_for_sn_computation_for_gas[0] $
			  and wavelength_output_rest_frame_log le $
			  w_range_for_sn_computation_for_gas[1] )

	    for i = 0, n_elements(bin_sn_ems_real)-1 do begin
		MDAP_CALCULATE_SPECTRUM_SN, best_fit_model_ems[i,indxx], residuals_ems[i,indxx], $
					    wavelength_output_rest_frame_log[indxx], $
					    sn_per_angstorm, /rms
		bin_sn_ems_real[i] = sn_per_angstorm[0]
	    endfor

	    ;junk = temporary(sn_per_angstorm)
	endif

	; Print progress
	printf,1,'[INFO] datacube '+root_name+' spectral fitting 3'

	; Save the progress
	if save_intermediate_steps eq 1 then begin
	    save, filename=root_name+mode+'_block4c.idl', /variables
	endif
	; END BLOCK 4 ------------------------------------------------------------------


	; BLOCK 5 ----------------------------------------------------------------------
	;	Description?
	;	- Measure of the equivalent width of absorption line indices
	;-------------------------------------------------------------------------------

	; Initialize the instrumental resolution.  TODO: use fits extension
	if n_elements(instrumental_fwhm_file) ne 0 then begin	; use input file
	    fwhm_instr=interpol(fwhm_ang,ww,wavelength_output_tpl) 
	endif else begin					; set to a constant
	    fwhm_instr=wavelength_output_tpl*0.+2.54		; TODO: Don't like this practice
	endelse

	;resolution_x_lick=[4000.,4400.,4900.,5400.,6000.]   ;CONTROLLARE
	;lick_fwhm_y=[11.5,9.2,8.4,8.4,9.8]                     ;CONTROLLARE
	;lick_resolution_tmp=interpol(lick_fwhm_y,resolution_x_lick,wavelength_output_tpl)
	;rrr=poly_fit(wavelength_output_tpl,lick_resolution_tmp,4,yfit=lick_resolution)

	;miles_resolution = fwhm_instr*0.+2.54

	; Determine the difference in instrumental resolution
	lick_resolution = fwhm_instr*0.+8.4		; Spectral resolution of Lick indices
	fwhm_diff_indices = sqrt(lick_resolution^2 - fwhm_instr^2)

	; Check where the Lick instrumental resolution is larger than the data resolution
	indici = where(finite(fwhm_diff_indices) ne 1)
	if indici[0] ne -1 then fwhm_diff_indices[indici] = 0.	; TODO: Forces them to be the same

	; Only measure the indices if required
	if mdap_measure_indices_version gt mdap_measure_indices_version_previous $
	   or execute_all_modules eq 1 then begin
	  
	    ; Warn user
	    print, 'measuring indices on ' + MDAP_STC(n_elements(best_template_tpl[*,0]),/integer) $
		   +' spectra'

	    ; Measure the indices
	    MDAP_MEASURE_INDICES, absorption_line_indices, wavelength_output_tpl, $
				  galaxy_minus_ems_fit_model_tpl, best_template_tpl, $
				  best_template_LOSVD_conv_tpl, stellar_kinematics_tpl[*,0], $
				  residuals_tpl, fwhm_diff_indices, abs_line_indices, $
				  abs_line_indices_errors, abs_line_indices_template, $
				  abs_line_indices_template_losvd, dir=output_dir, /noplot

	    ; Add the version to the header
	    SXADDPAR, header_2d, 'BLOCK5', mdap_measure_indices_version, $
		      'mdap_measure_indices version'

	    ; Re-execute all the remaining modules, if that wasn't already set
	    execute_all_modules=1
	endif

	; Print progress
	printf, 1, '[INFO] datacube '+root_name+' indices measured on ' $
		+ MDAP_STC(n_elements(best_template_tpl[*,0]), /integer)+' spectra'

	; Save the progress
	if save_intermediate_steps eq 1 then begin
	    save, filename=root_name+mode+'_block5.idl', /variables
	endif
	; END BLOCK 5 ------------------------------------------------------------------


	; BLOCK 6 ----------------------------------------------------------------------
	;	Description?
	;	- Determine azimuthally averaged radial profiles
	;   TODO: Isn't this essentially a combination of the other blocks but for a
	;         4th binning scheme?
	;-------------------------------------------------------------------------------

	; Create error array using fit residuals, if required
	if mdap_spatial_radial_binning_version gt mdap_spatial_radial_binning_version_previous $
	   or execute_all_modules eq 1 then begin
	   
	    ; TODO: Why is mdap_get_erro_from_residual not used? Function obsolete?
	    ;mdap_get_error_from_residual,residuals_ems,galaxy_minus_ems_fit_model_ems,input_errors

	    sz_ems = size(log_err_ems)
	    nl = n_elements(wavelength_output_rest_frame_log)
	    input_errors = fltarr(sz_ems[2], nl)
	    log_step_gal = LOG_WAV_EMS[1]-LOG_WAV_EMS[0]
	    for i = 0, sz_ems[2]-1 do begin
		rf_gal_lam = exp(log_wav_ems-stellar_kinematics_ems[i,0] / velscale*(log_step_gal))
		input_errors[i,*] = interpol(log_err_ems[*,i], rf_gal_lam, $
					     wavelength_output_rest_frame_log)
	    endfor

	    ; Perform the radial binning
	    MDAP_SPATIAL_RADIAL_BINNING, bin_sn_ems_real, x2d_reconstructed, y2d_reconstructed, $
					 spatial_binning_ems, xbin_ems, ybin_ems, ell, pa, $
					 galaxy_minus_ems_fit_model_ems, input_errors, $
					 wavelength_output_rest_frame_log, spatial_binning_rad, $
					 r_bin, r_bin_lo, r_bin_up, r2d_bin, r2d_bin_lo, $
					 r2d_bin_up, radially_binned_spectra, $
					 radially_binned_errors, $
					 output_lrange=trim_wav_range_radial_binning, $
					 output_wav=output_wav, $
					 n_elements_bin=nelements_within_bin_radial, $
					 low_radial_bins_user_inputs=low_radial_bins_user_inputs, $
					 upper_radial_bins_user_inputs=upper_radial_bins_user_inputs, $
					 Reff_=Reff, PSFsize_=PSFsize, $
					 add_default_bins=add_default_bins

	    ; Print progress
	    printf, 1, '[INFO] datacube '+root_name+' radial binning: ' $
		       +MDAP_STC(n_elements(r_bin),/integer)+' bins'


	    ; Perform the spectral fit of the radially binned data
	    star_kin_starting_guesses_rbin = fltarr(n_elements(r_bin),4)
	    gas_kin_starting_guesses_rbin  = fltarr(n_elements(r_bin),2)
	    star_kin_starting_guesses_rbin[*,1] = velocity_dispersion_initial_guess[0]
	    gas_kin_starting_guesses_rbin[*,1]=50.
	    
	    loglam_gal_rbin=alog(output_wav)
	    
	    MDAP_SPECTRAL_FITTING, radially_binned_spectra, radially_binned_errors, $
				   loglam_gal_rbin, library_log_ems, log_wav_library_ems, $
				   velscale, stellar_kinematics_rbin, stellar_kinematics_rbin_err, $
				   stellar_weights_rbin, emission_line_kinematics_rbin, $
				   emission_line_kinematics_rbin_err, $
				   emission_line_kinematics_rbin_individual, $
				   emission_line_kinematics_rbin_individual_err, $
				   emission_line_intens_rbin, emission_line_intens_rbin_err, $
				   emission_line_fluxes_rbin, emission_line_fluxes_rbin_err, $
				   emission_line_equivW_rbin, emission_line_equivW_rbin_err, $
				   wavelength_input=output_wav, wavelength_output_rbin, $
				   best_fit_model_rbin, galaxy_minus_ems_fit_model_rbin, $
				   best_template_rbin, best_template_LOSVD_conv_rbin, $
				   reddening_rbin, reddening_rbin_err, residuals_rbin, $
				   star_kin_starting_guesses=star_kin_starting_guesses_rbin, $
				   gas_kin_starting_guesses=gas_kin_starting_guesses_rbin, $
				   emission_line_file=emission_line_file_rad_binning, $
				   extra_inputs=spectral_fit_par_rbin, $
				   fwhm_instr_kmsec_matrix=fwhm_instr_kmsec_matrix/3., $
				   range_v_star=[-50.,50.], range_v_gas=[-50.,50.], $
				   mask_range=mask_range, external_library=external_library, /quiet

	    ; Print progress
	    printf,1,'[INFO] datacube '+root_name+' radial binning: spectral fitting' 

	    ;LC: calculate the real S/N of binned spectra
	    ; TODO: Is this really the "real S/N"?
	    bin_sn_rad_real = r_bin
	    indxx = where(wavelength_output_rbin ge w_range_for_sn_computation[0] $
			  and wavelength_output_rbin le w_range_for_sn_computation[1] )

	    for i = 0, n_elements(bin_sn_rad_real) -1 do begin
		MDAP_CALCULATE_SPECTRUM_SN, best_fit_model_rbin[i,indxx], residuals_rbin[i,indxx], $
					    wavelength_output_rbin[indxx], sn_per_angstorm, /rms
		bin_sn_rad_real[i] = sn_per_angstorm[0]
	    endfor

	    ;junk = temporary(sn_per_angstorm)

	    ; Measure the absorption-line indices using the radially binned spectra

	    ;lick_resolution_tmp=interpol(lick_fwhm_y,resolution_x_lick,wavelength_output_rbin)
	    ;rrr=poly_fit(wavelength_output_rbin,lick_resolution_tmp,4,yfit=lick_resolution)
	    ;fwhm_diff_indices=sqrt(double(lick_resolution)^2.-double(fwhm_instr)^2.)*0.

	    ; Initialize the instrumental resolution.  TODO: use fits extension
	    if n_elements(instrumental_fwhm_file) ne 0 then begin		; use input file
		fwhm_instr=interpol(fwhm_ang,ww,wavelength_output_rbin) 
	    endif else begin						; set to a constant
		fwhm_instr=wavelength_output_rbin*0.+2.54		; TODO: bad practice?
	    endelse

	    ;miles_resolution = fwhm_instr*0.+2.54
	    ;fwhm_diff_indices=sqrt(miles_resolution^2-fwhm_instr^2) ;fwhm in angstrom

	    ; Determine the difference in instrumental resolution
	    lick_resolution = fwhm_instr*0.+8.4
	    fwhm_diff_indices=sqrt(lick_resolution^2 - fwhm_instr^2)

	    ; Check where the Lick instrumental resolution is larger than the data resolution
	    indici = where(finite(fwhm_diff_indices) ne 1)
	    if indici[0] ne -1 then fwhm_diff_indices[indici] = 0.  ; TODO: Forced to be the same

	    ; Measure the absorption-line indices
	    MDAP_MEASURE_INDICES, absorption_line_indices, wavelength_output_rbin, $
				  galaxy_minus_ems_fit_model_rbin, best_template_rbin, $
				  best_template_LOSVD_conv_rbin, stellar_kinematics_rbin[*,0], $
				  residuals_rbin, fwhm_diff_indices, abs_line_indices_rbin, $
				  abs_line_indices_errors_rbin, abs_line_indices_template_rbin, $
				  abs_line_indices_template_losvd_rbin, dir=output_dir+'rbin_', $
				  /noplot

	    ; Print progress
	    printf,1,'[INFO] datacube '+root_name+' radial binning: measured indices'

	    ; Add the version to the header
	    SXADDPAR, header_2d, 'BLOCK6', mdap_spatial_radial_binning_version, $
		      'mdap_spatial_radial_binning_version'

	    ; Re-execute all the remaining modules, if that wasn't already set
	    execute_all_modules=1
	endif

	; Save progress
	if save_intermediate_steps eq 1 then begin
	    save, filename=root_name+mode+'_block6a.idl', /variables
	endif

	; Determine the stellar rotation curves using the kinemetry package
	MDAP_DO_KINEMETRY, signal2d_reconstructed, x2d_reconstructed, y2d_reconstructed, xbin_str, $
			   ybin_str, stellar_kinematics_str[*,0], stellar_kinematics_str_err[*,0], $
			   PA_kin_str, PA_kin_std_str, q_kin_str, q_kin_std_str, Vsyst_str, $
			   Vsyst_std_str, Rad_kin_str, Vrot_str, Vrot_err_str, Vexp_str, $
			   Vexp_err_str

	; Determine the gas rotation curves using the kinemetry package
	MDAP_DO_KINEMETRY, signal2d_reconstructed, x2d_reconstructed, y2d_reconstructed, xbin_ems, $
			   ybin_ems, emission_line_kinematics_ems[*,0], $
			   emission_line_kinematics_EMS_err[*,0], PA_kin_ems, PA_kin_std_ems, $
			   q_kin_ems, q_kin_std_ems, Vsyst_ems, Vsyst_std_ems, Rad_kin_ems, $
			   Vrot_ems, Vrot_err_ems, Vexp_ems, Vexp_err_ems

	; Determine other radial profiles
	;	- lambda
	;	- V/sigma
	;	- sigma
	;  TODO: requires testing (LC)

	if n_elements(A_for_rprofiles) eq 0 then begin
	    A_rprofiles=Rad_kin_str
	endif else begin
	    A_rprofiles=A_for_rprofiles
	endelse

	MDAP_DO_K_RPROFILES, A_rprofiles, stellar_kinematics_str[*,0]-Vsyst_str[0], $
			     stellar_kinematics_str_err[*,0], stellar_kinematics_str[*,1], $
			     stellar_kinematics_str_err[*,1], xbin_str, ybin_str, $
			     total(log_spc_str,1), ell, pa, lambda_profile, lambda_profile_err, $
			     vsigma_profile, vsigma_profile_err, sigma_profile, sigma_profile_err

	radii_rprofiles = sqrt(A_rprofiles^2 * (1.-ell))

	; END BLOCK 6 ------------------------------------------------------------------


	; END computation of "model-independent" data products #########################
	; ##############################################################################


	; BEGIN writing "model-independent" products -----------------------------------
	; TODO: Needs to be optimized (LC)
	; TODO: Why is this not done after every block?

	; Index setup
	NINDICES = n_elements(ABS_LINE_INDICES[0,*])			; # of indices measured
	MDAP_READ_INDICES_DEFINITIONS, absorption_line_indices, indices=indices
	READCOL, emission_line_file_binning_3, cnt, ln_name, ln_wav, comment='#', $
	         format='I, A, A', /silent

	NLINES = n_elements(cnt)					; # of lines

	; TODO: Change names of extensions?
	;	Better the comments once I figure out what "signal" is

	; Signal extension
	WRITEFITS, output_filefits, signal2d_reconstructed, header_2d
	k=1

	; Noise extension
	MDAP_ADD_FITS_LAYER, output_filefits, noise2d_reconstructed, k, 'EXTNAME', 'noise'
	k=k+1

	; Binning extension
	MDAP_ADD_FITS_LAYER, output_filefits, spatial_binning_tpl, k, 'EXTNAME','Binning_map_1'
	k=k+1

	; Set the table column names
	stringa = [ "X", "Y", "area_bin", "StoN", "Nelements" ]
	for i = 0, NINDICES-1 do $
	    stringa = [ stringa, indices[i].name, indices[i].name+"_err" ]

	stringa2="{"
	for i = 0, n_elements(stringa)-2 do $
	    stringa2 = stringa2+stringa[i]+':0.,'
 	stringa2 = stringa2+stringa[i]+':0.}'

	d = execute('str='+stringa2)
	p1=replicate(str,n_elements(xbin_tpl))
	p1.x=xbin_tpl
	p1.y=ybin_tpl
	p1.area_bin = AREA_BINS_TPL
	; p1.reddening = reddening_tpl
	p1.StoN = bin_sn_tpl_real
	p1.Nelements=nelements_within_bin_tpl

	for i = 0, NINDICES-1 do begin
	    d=execute('p1.'+indices[i].name+'=ABS_LINE_INDICES[*,i]')
	    d=execute('p1.'+indices[i].name+'_err=ABS_LINE_INDICES_ERRORS[*,i]')
	endfor
	
	MWRFITS, p1, output_filefits, /silent			;?
	h1=HEADFITS(output_filefits,EXTEN = k)			;?
	SXADDPAR, h1, 'EXTNAME', 'Binning_1_data'		;Add the extension name to header
	MODFITS, output_filefits, 0, h1, exten_no=k		;?

	; Add extension with stellar kinematics
	k=k+1
	MDAP_ADD_FITS_LAYER, output_filefits, spatial_binning_str, k, 'EXTNAME', 'Binning_map_2'
	k=k+1

	stringa = [ "X", "Y", "area_bin", "StoN", "Nelements", "Vel", "Vel_err", "Disp", $
		    "Disp_err", "H3", "H3_err", "H4", "H4_err", "Chi2_DOF", "Vel_2moms", $
		    "Vel_2moms_err", "Disp_2moms", "Disp_2moms_err" ]
	stringa2="{"
	for i = 0, n_elements(stringa)-2 do $
	    stringa2 = stringa2+stringa[i]+':0.,'
	stringa2 = stringa2+stringa[i]+':0.}'

	d = execute('str='+stringa2)
	p2=replicate(str,n_elements(xbin_str))
	p2.x=xbin_str
	p2.y=ybin_str
	p2.area_bin = AREA_BINS_STR
	p2.StoN = bin_sn_str_real
	p2.Nelements=nelements_within_bin_str
	p2.Vel=STELLAR_KINEMATICS_STR[*,0]
	p2.Vel_err=STELLAR_KINEMATICS_STR_ERR[*,0]
	p2.Disp=STELLAR_KINEMATICS_STR[*,1]
	p2.Disp_err=STELLAR_KINEMATICS_STR_ERR[*,1]
	p2.H3=STELLAR_KINEMATICS_STR[*,2]
	p2.H3_err=STELLAR_KINEMATICS_STR_ERR[*,2]
	p2.H4=STELLAR_KINEMATICS_STR[*,3]
	p2.H4_err=STELLAR_KINEMATICS_STR_ERR[*,3]
	p2.Chi2_DOF=STELLAR_KINEMATICS_STR[*,4]

	; If fit with h3=h4=0 not performed, just copy V and sig from free fit
	if n_elements(spectra_fittin_parameters_patial_binning_2_two_moments) eq 0 then begin
	    p2.Vel_2moms=STELLAR_KINEMATICS_STR[*,0]
	    p2.Vel_2moms_err=STELLAR_KINEMATICS_STR_ERR[*,0]
	    p2.Disp_2moms=STELLAR_KINEMATICS_STR[*,1]
	    p2.Disp_2moms_err=STELLAR_KINEMATICS_STR_ERR[*,1]
	endif else begin
	; ... otherwise, save results from fit with h3=h4=0
	    p2.Vel_2moms=STELLAR_KINEMATICS_STR_TWO_MOMENTS[*,0]
	    p2.Vel_2moms_err=STELLAR_KINEMATICS_STR_TWO_MOMENTS_ERR[*,0]
	    p2.Disp_2moms=STELLAR_KINEMATICS_STR_TWO_MOMENTS[*,1]
	    p2.Disp_2moms_err=STELLAR_KINEMATICS_STR_TWO_MOMENTS_ERR[*,1]
	endelse
	
	MWRFITS, p2, output_filefits, /silent
	h1=HEADFITS(output_filefits, EXTEN=k)
	SXADDPAR, h1, 'EXTNAME', 'Binning_2_data'
	MODFITS, output_filefits, 0, h1, exten_no=k

	; Store emission-line kinematics, fluxes, and equivalent widths
	k=k+1
	MDAP_ADD_FITS_LAYER, output_filefits, spatial_binning_ems, k, 'EXTNAME', 'Binning_map_3'
	k=k+1

	; Set column names
	stringa = [ "X", "Y", "area_bin", "StoN", "Nelements", "reddening_star", $
		    "reddening_star_err", "reddening_gas", "reddening_gas_err", "Chi2_DOF" ]
		    
	for i = 0, NLINES-1 do begin
	    stringa = [ stringa, ln_name[i]+'_'+MDAP_STC(round(float(ln_wav[i])))+'_amplitude', $
		        ln_name[i]+'_'+MDAP_STC(round(float(ln_wav[i])))+'_amplitude_err' ]
	endfor
	for i = 0, NLINES-1 do begin
	    stringa = [ stringa, ln_name[i]+'_'+MDAP_STC(round(float(ln_wav[i])))+'_flux', $
			ln_name[i]+'_'+MDAP_STC(round(float(ln_wav[i])))+'_flux_err' ]
	endfor
	for i = 0, NLINES-1 do begin
	    stringa = [ stringa, ln_name[i]+'_'+MDAP_STC(round(float(ln_wav[i])))+'_EW', $
			ln_name[i]+'_'+MDAP_STC(round(float(ln_wav[i])))+'_EW_err' ]
	endfor 

	stringa2="{"
	for i = 0, n_elements(stringa)-2 do $
	    stringa2 = stringa2+stringa[i]+':0.,'
	stringa2 = stringa2+stringa[i]+':0.}'

	d = execute('str='+stringa2)
	p3=replicate(str,n_elements(xbin_ems))
	p3.x=xbin_ems
	p3.y=ybin_ems
	p3.area_bin = AREA_BINS_EMS
	p3.StoN = bin_sn_ems_real
	p3.Nelements=nelements_within_bin_ems
	; p3.Vel=EMISSION_LINE_KINEMATICS_EMS[*,0]
	; p3.Vel_err=EMISSION_LINE_KINEMATICS_EMS_ERR[*,0]
	; p3.Disp=EMISSION_LINE_KINEMATICS_EMS[*,1]
	; p3.Disp_err=EMISSION_LINE_KINEMATICS_EMS_ERR[*,1]
	p3.reddening_star = reddening_ems[*,0]
	p3.reddening_star_err = reddening_ems_err[*,0]
	p3.reddening_gas = reddening_ems[*,1]
	p3.reddening_gas_err = reddening_ems_err[*,1]
	p3.Chi2_DOF=STELLAR_KINEMATICS_EMS[*,4]

	for i = 0, NLINES-1 do begin
	    d=execute('p3.'+ln_name[i]+'_'+mdap_stc(round(float(ln_wav[i]))) $
		      +'_amplitude=EMISSION_LINE_INTENS_EMS[*,i]')
	    d=execute('p3.'+ln_name[i]+'_'+mdap_stc(round(float(ln_wav[i]))) $
		       +'_amplitude_err=EMISSION_LINE_INTENS_EMS_ERR[*,i]')
	endfor
	for i = 0, NLINES-1 do begin
	    d=execute('p3.'+ln_name[i]+'_'+mdap_stc(round(float(ln_wav[i]))) $
		      +'_flux=EMISSION_LINE_FLUXES_EMS[*,i]')
	    d=execute('p3.'+ln_name[i]+'_'+mdap_stc(round(float(ln_wav[i]))) $
		      +'_flux_err=EMISSION_LINE_FLUXES_EMS_ERR[*,i]')
	endfor
	for i = 0, NLINES-1 do begin
	    d=execute('p3.'+ln_name[i]+'_'+mdap_stc(round(float(ln_wav[i]))) $
		      +'_EW=EMISSION_LINE_EQUIVW_EMS[*,i]')
	    d=execute('p3.'+ln_name[i]+'_'+mdap_stc(round(float(ln_wav[i]))) $
		      +'_EW_err=EMISSION_LINE_EQUIVW_EMS_ERR[*,i]')
	endfor 

; for i = 0, NLINES-1 do begin 
;    d=execute('p3.'+ln_name[i]+'_'+mdap_stc(round(float(ln_wav[i])))+'_Vel=emission_line_kinematics_ems_individual[*,i,0]')
;    d=execute('p3.'+ln_name[i]+'_'+mdap_stc(round(float(ln_wav[i])))+'_VEl_err=emission_line_kinematics_ems_individual_err[*,i,0]')
; endfor 
; for i = 0, NLINES-1 do begin 
;    d=execute('p3.'+ln_name[i]+'_'+mdap_stc(round(float(ln_wav[i])))+'_Sig=emission_line_kinematics_ems_individual[*,i,1]')
;    d=execute('p3.'+ln_name[i]+'_'+mdap_stc(round(float(ln_wav[i])))+'_Sig_err=emission_line_kinematics_ems_individual_err[*,i,1]')
; endfor 

	MWRFITS, p3, output_filefits, /silent
	h1=HEADFITS(output_filefits, EXTEN=k)
	SXADDPAR, h1, 'EXTNAME', 'Binning_3_lines'
	MODFITS, output_filefits, 0, h1, exten_no=k

	; storing  gas kinematics
	k=k+1
	stringa = [ "X", "Y", "mean_VEL", "err_mean_VEL", "mean_DISP", "err_mean_DISP" ]
	for i = 0, NLINES-1 do begin
	    stringa = [ stringa, ln_name[i]+'_'+mdap_stc(round(float(ln_wav[i])))+'_Vel', $
			ln_name[i]+'_'+mdap_stc(round(float(ln_wav[i])))+'_Vel_err' ]
	endfor
	for i = 0, NLINES-1 do begin
	    stringa = [ stringa, ln_name[i]+'_'+mdap_stc(round(float(ln_wav[i])))+'_Sig', $
			ln_name[i]+'_'+mdap_stc(round(float(ln_wav[i])))+'_Sig_err' ]
	endfor 

	stringa2="{"
	for i = 0, n_elements(stringa)-2 do $
	    stringa2 = stringa2+stringa[i]+':0.,'
	stringa2 = stringa2+stringa[i]+':0.}'

	d = execute('str='+stringa2)
	p3k=replicate(str,n_elements(xbin_ems))
	p3k.x=xbin_ems
	p3k.y=ybin_ems
	p3k.mean_VEL=EMISSION_LINE_KINEMATICS_EMS[*,0]
	p3k.err_mean_VEL=EMISSION_LINE_KINEMATICS_EMS_ERR[*,0]
	p3k.mean_DISP=EMISSION_LINE_KINEMATICS_EMS[*,1]
	p3k.err_mean_DISP=EMISSION_LINE_KINEMATICS_EMS_ERR[*,1]

	for i = 0, NLINES-1 do begin
	    d=execute('p3k.'+ln_name[i]+'_'+mdap_stc(round(float(ln_wav[i]))) $
		      +'_Vel=emission_line_kinematics_ems_individual[*,i,0]')
	    d=execute('p3k.'+ln_name[i]+'_'+mdap_stc(round(float(ln_wav[i]))) $
		      +'_VEl_err=emission_line_kinematics_ems_individual_err[*,i,0]')
	endfor 

	for i = 0, NLINES-1 do begin
	    d=execute('p3k.'+ln_name[i]+'_'+mdap_stc(round(float(ln_wav[i]))) $
		      +'_Sig=emission_line_kinematics_ems_individual[*,i,1]')
	    d=execute('p3k.'+ln_name[i]+'_'+mdap_stc(round(float(ln_wav[i]))) $
		      +'_Sig_err=emission_line_kinematics_ems_individual_err[*,i,1]')
	endfor 

	MWRFITS, p3k, output_filefits, /silent
	h1=HEADFITS(output_filefits, EXTEN=k)
	SXADDPAR, h1, 'EXTNAME', 'Binning_3_kinematics'
	MODFITS, output_filefits, 0, h1, exten_no=k

	; Store azimuthally averaged data
	k=k+1
	MDAP_ADD_FITS_LAYER, output_filefits, spatial_binning_rad, k, 'EXTNAME', $
			     'Binning_map_radial'
	k=k+1

	stringa = [ "amaj", "amaj_lo", "amaj_up", "StoN", "Nelements", "Disp", "Disp_err", $
		    "Chi2_DOF" ]

	for i = 0, NINDICES-1 do $
	    stringa = [stringa,indices[i].name+"_rad",indices[i].name+"_rad_err"]
	    
	stringa2="{"
	for i = 0, n_elements(stringa)-2 do $
	    stringa2 = stringa2+stringa[i]+':0.,'
	stringa2 = stringa2+stringa[i]+':0.}'

	d = execute('str='+stringa2)
	p4=replicate(str,n_elements(r_bin))
	p4.amaj=r_bin
	p4.amaj_lo=r_bin_lo
	p4.amaj_up=r_bin_up
	p4.StoN=bin_sn_rad_real
	p4.Nelements=nelements_within_bin_radial
	p4.Disp=stellar_kinematics_rbin[*,1]
	p4.Disp_err=stellar_kinematics_rbin_err[*,1]
	p4.Chi2_DOF=stellar_kinematics_rbin[*,4]

	for i = 0, NINDICES-1 do begin
	    d=execute('p4.'+indices[i].name+'_rad=ABS_LINE_INDICES_RBIN[*,i]')
	    d=execute('p4.'+indices[i].name+'_rad_err=ABS_LINE_INDICES_ERRORS_RBIN[*,i]')
	endfor 

	MWRFITS, p4, output_filefits, /silent
	h1=HEADFITS(output_filefits, EXTEN=k)
	SXADDPAR, h1, 'EXTNAME', 'Binning_radial_data'
	modfits, output_filefits, 0, h1, exten_no=k

	; Store the stellar rotation curve
	k=k+1
	stringa = [ "a_rot_curve", "PA_kin_stars", "PA_kin_std_stars", "q_kin_stars", $
		    "q_kin_std_stars","Vsyst_stars","Vsyst_std_stars", "Vrot_stars", $
		    "Vrot_err_stars","Vexp_stars", "Vexp_err_stars" ]

	stringa2="{"
	for i = 0, n_elements(stringa)-2 do $
	    stringa2 = stringa2+stringa[i]+':0.,'
	stringa2 = stringa2+stringa[i]+':0.}'

	d = execute('str='+stringa2)
	p5=replicate(str,n_elements(Rad_kin_str))
	p5.a_rot_curve=Rad_kin_str
	p5.PA_kin_stars=replicate(PA_kin_str,n_elements(Rad_kin_str))
	p5.PA_kin_std_stars=replicate(PA_kin_std_str,n_elements(Rad_kin_str))
	p5.q_kin_stars=replicate(q_kin_str,n_elements(Rad_kin_str))
	p5.q_kin_std_stars=replicate(q_kin_std_str,n_elements(Rad_kin_str))
	p5.Vsyst_stars=Vsyst_str
	p5.Vsyst_std_stars=Vsyst_std_str
	p5.Vrot_stars= Vrot_str
	p5.Vrot_err_stars=Vrot_err_str
	p5.Vexp_stars=Vexp_str
	p5.Vexp_err_stars=Vexp_err_str

	;p5.x0=replicate(Xcenter_used_for_stellar_rot_curve,n_elements(Rad_kin_str)) 
	;p5.y0=replicate(Ycenter_used_for_stellar_rot_curve,n_elements(Rad_kin_str)) 

	MWRFITS, p5, output_filefits, /silent
	h1=HEADFITS(output_filefits, EXTEN=k)
	SXADDPAR, h1, 'EXTNAME', 'Stars Rotation'
	MODFITS, output_filefits, 0, h1, exten_no=k

	; Store the gas rotation curve
	k=k+1
	stringa = [ "a_rot_curve", "PA_kin_gas", "PA_kin_std_gas", "q_kin_gas", "q_kin_std_gas", $
		    "Vsyst_gas", "Vsyst_std_gas", "Vrot_gas", "Vrot_err_gas", "Vexp_gas", $
		    "Vexp_err_gas" ]

	stringa2="{"
	for i = 0, n_elements(stringa)-2 do $
	    stringa2 = stringa2+stringa[i]+':0.,'
	stringa2 = stringa2+stringa[i]+':0.}'

	d = execute('str='+stringa2)
	p6=replicate(str,n_elements(Rad_kin_ems))
	p6.a_rot_curve=Rad_kin_ems
	p6.PA_kin_gas=replicate(PA_kin_ems,n_elements(Rad_kin_ems))
	p6.PA_kin_std_gas=replicate(PA_kin_std_ems,n_elements(Rad_kin_ems))
	p6.q_kin_gas=replicate(q_kin_ems,n_elements(Rad_kin_ems))
	p6.q_kin_std_gas=replicate(q_kin_std_ems,n_elements(Rad_kin_ems))
	p6.Vsyst_gas=Vsyst_ems
	p6.Vsyst_std_gas=Vsyst_std_ems
	p6.Vrot_gas= Vrot_ems
	p6.Vrot_err_gas=Vrot_err_ems
	p6.Vexp_gas=Vexp_ems
	p6.Vexp_err_gas=Vexp_err_ems

	;p6.x0=replicate(Xcenter_used_for_gas_rot_curve,n_elements(Rad_kin_ems))
	;p6.y0=replicate(Ycenter_used_for_gas_rot_curve,n_elements(Rad_kin_ems)) 

	MWRFITS, p6, output_filefits, /silent
	h1=HEADFITS(output_filefits, EXTEN=k)                                         
	SXADDPAR, h1, 'EXTNAME', 'GAS rotation'
	MODFITS, output_filefits, 0, h1, exten_no=k

	; Store the lambda, V/sigma, and sigma profiles
	k = k+1
	stringa = [ "radii", "Amaj", "lambda_profile", "lambda_profile_err", "vsigma_profile", $
		    "vsigma_profile_err", "sigma_profile", "sigma_profile_err" ]

	stringa2="{"
	for i =0, n_elements(stringa)-2 do $
	    stringa2=stringa2+stringa[i]+':0.,'
	stringa2 = stringa2+stringa[i]+':0.}'

	d = execute('str='+stringa2)
	p7=replicate(str,n_elements(A_rprofiles))
	p7.radii=radii_rprofiles
	p7.Amaj=A_rprofiles
	p7.lambda_profile=lambda_profile
	p7.lambda_profile_err=lambda_profile_err
	p7.vsigma_profile=vsigma_profile
	p7.vsigma_profile_err=vsigma_profile_err
	p7.sigma_profile=sigma_profile
	p7.sigma_profile_err=sigma_profile_err

	MWRFITS, p7, output_filefits, /silent
	h1=HEADFITS(output_filefits, EXTEN=k)
	SXADDPAR, h1, "EXTNAME", "lambda, V/S, S"
	MODFITS, output_filefits, 0, h1, exten_no=k
	; END writing "model-independent" products -------------------------------------

	; Print progress
	printf, 1, '[INFO] '+root_name+' output saved: '+output_filefits
	save, /variables, filename=output_idlsession
	printf, 1, '[INFO] '+root_name+' output saved: '+output_idlsession

	; ##############################################################################
	; BEGIN computation of "model-dependent" data products #########################

	; END computation of "model-dependent" data products ###########################
	; ##############################################################################

	; Print execution summary
	printf, 1, root_name+'_'+mode+' execution time '+mdap_stc(systime(/seconds)/60./60.-t0) $
		   +' hours'
	close, 1

	; Print execution file
	openw, 1, root_name+'_'+mode+'.done'
	printf, 1, root_name+'_'+mode+' execution time '+mdap_stc(systime(/seconds)/60./60.-t0) $
	           +' hours'
	close, 1

end

; ;uncomment these lines when performing parallelization on SCIAMA.
; PREF_SET,'IDL_CPU_TPOOL_NTHREADS',1,/COMMIT
; manga_dap,!!num!!
; EXIT
; END



