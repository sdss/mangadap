;+
; NAME:
;	MANGA_DAP
;
; PURPOSE:
;	Main wrapping routine for MaNGA data analysis pipeline.
;
;	For now, see documentation here:
;		https://trac.sdss.org/wiki/MANGA/Data/DAPDevelopment/DAP_v0_9_summary
;
; CALLING SEQUENCE:
;	MANGA_DAP, input_number, /nolog, /quiet
;
; INPUTS:
;    input_number int
;		Index number (line number minus one) of the data to analyze from
;		total_filelist file defined by the MDAP_EXECUTION_SETUP procedure.
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;	/nolog
;		Do not produce the log file.
;
;	/quiet
;		Limit (eliminate?) the output printed to the screen.
;
; OUTPUT:
;
; OPTIONAL OUTPUT:
;
; TODO:
;	- Include MW extinction curve
;	- Put all run logic and file IO checks up front!
;	- Somehow allow to check for overlapping objects (multiple redshifts in
;	  a single IFU)
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;	MDAP_EXECUTION_SETUP
;	READCOL
;	MDAP_SETUP_IO
;	MDAP_DRP_CHECK_FILETYPE
;	MDAP_CHECK_FILE()
;	MDAP_BUILD_EXECUTION_PLANS
;	MDAP_PRINT_EXECUTION_PLAN
;	SXADDPAR
;	MDAP_STC()
;	MDAP_READ_DRP_FITS
;	MDAP_GET_SPAXEL_SIZE
;	MDAP_FIDUCIAL_BIN_XY
;	MDAP_SET_TPL_LIB_OUTPUT_FILE()
;	MDAP_READ_TEMPLATE_LIBRARY
;	MDAP_MATCH_RESOLUTION_TPL2OBJ
;	MDAP_TRIM_SPECTRAL_RANGE
;	MDAP_VELOCITY_SCALE()
;	MDAP_RESAMPLE_TEMPLATES
;	MDAP_NORMALIZE_TEMPLATES
;	WRITEFITS
;	MDAP_ADD_FITS_LAYER
;	MDAP_SELECT_GOOD_SPECTRA
;	MDAP_SELECT_WAVE
;	MDAP_CALCULATE_SN
;	MDAP_WRITE_OUTPUT
;	MDAP_SPATIAL_BINNING
;	MDAP_READ_EMISSION_LINE_PARAMETERS
;	MDAP_ERASEVAR
;	MDAP_SPECTRAL_FITTING
;
; INTERNAL SUPPORT ROUTINES:
;	MDAP_FILE_EXTENSION_DATE()
;	MDAP_LOG_INSTANTIATE
;	MDAP_LOG_CLOSE
;	MDAP_GENERATE_OUTPUT_FILE_NAMES
;	MDAP_SET_TPL_LIB_OUTPUT_FILE
;
; REVISION HISTORY:
;	01 Sep 2014: Copied from v0_8 version by L. Coccato.
;	02 Sep 2014: (KBW) Formating and minor edits
;	03 Sep 2014: (KBW) Basic compilation errors
;	26 Oct 2014: (KBW) v0.9 (see documentation linked to above)
;-
;------------------------------------------------------------------------------

FUNCTION MDAP_FILE_EXTENSION_DATE

	vec = fix(date_conv( systime(/julian,/utc), 'V'))
	ext = MDAP_STC(vec[0], /integer)+'_'+$ 
	      MDAP_STC(vec[1], /integer)+'_'+$
	      MDAP_STC(vec[2], /integer)+'_'+$
	      MDAP_STC(vec[3], /integer)+'_'+$
	      MDAP_STC(vec[4], /integer)

	return, ext
END

PRO MDAP_LOG_INSTANTIATE, $
		output_dir, manga_dap_version, signifier, total_filelist, input_number, $
		datacube_name, mode, file_root, guess_vel, guess_sig, ell, pa, Reff, $
		file_unit

	ofile = output_dir+'/mdap_'+MDAP_FILE_EXTENSION_DATE()+'.log'
	openw, file_unit, ofile, /get_lun
	
	cd, current=dir
	printf, file_unit, 'DAP called from: ', dir
	printf, file_unit, 'Start systime(): ', systime()
	printf, file_unit, ''
	printf, file_unit, 'DAP Version: ', manga_dap_version
	printf, file_unit, ''
	printf, file_unit, 'MDAP_EXECUTION_SETUP signifier: ', signifier
	printf, file_unit, 'Input table: ', total_filelist
	printf, file_unit, 'Line to process: ', (input_number+1)
	printf, file_unit, ''
	printf, file_unit, 'Fits file to process: ', datacube_name
	printf, file_unit, 'File type: ', mode
	printf, file_unit, 'Output directory: ', output_dir
	printf, file_unit, 'Root name for output files: ', file_root
	printf, file_unit, ''
	printf, file_unit, 'Input parameters (or guesses):'
	printf, file_unit, '            Velocity: ', guess_vel
	printf, file_unit, '     Vel. dispersion: ', guess_sig
	printf, file_unit, '         Ellipticity: ', ell
	printf, file_unit, '      Position Angle: ', pa
	printf, file_unit, '    Effective radius: ', Reff
	printf, file_unit, ''

END

;output_filefits, 
;	printf, file_unit, 'Data products:'
;	printf, file_unit, '    Primary fits file: ', output_filefits
;	printf, file_unit, ''

PRO MDAP_LOG_CLOSE, $
		file_unit

	printf, file_unit, 'End systime(): ', systime()
	close, file_unit
	free_lun, file_unit
END

; TODO: Eventually this will generate the files names AND check them against
; existing files.  Any execution plan that matches IN DETAIL the execution plan
; of an existing file (apart from the analyses performed) will be used to input
; already completed analyses, unless the overwrite flag is flipped.   This check
; will *include* checks of the module versions used to generate the existing
; analysis.

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

pro MANGA_DAP, $
	input_number, nolog=nolog, quiet=quiet

	manga_dap_version = '0.9'	; set the version number
	t0=systime(/seconds)/60./60.	; start time

	; Execute the script the sets the user-level execution parameters
	MDAP_EXECUTION_SETUP, total_filelist, output_root_dir, tpl_library_keys, $
			      template_libraries, tpl_vacuum_wave, emission_line_parameters, $
			      ems_vacuum_wave, absorption_line_parameters, abs_vacuum_wave, $
			      signifier, bin_type, bin_par, w_range_sn, threshold_ston_bin, $
			      bin_weight_by_sn2, w_range_analysis, threshold_ston_analysis, $
			      analysis, tpl_lib_analysis, ems_par_analysis, abs_par_analysis, $
			      overwrite_flag, analysis_par, $
			      save_intermediate_steps=save_intermediate_steps, $
			      remove_null_templates=remove_null_templates, $
			      external_library=external_library

	; Check that certain things are setup
	if file_test(total_filelist) eq 0 then $
	    message, total_filelist+' does not exist!'

	; Output directory must have some non-zero size
	if strlen(output_root_dir) eq 0 then $
	    message, 'Output directory must have a non-zero length.  Use ./ for current directory'

	; Initialization ---------------------------------------------------------------
	; TODO: Put all this initialization stuff in a wrapper

	; Use configuration file to set some user-defined variables
;	print, 'Reading configuration file:'
;	READCOL, configure_file, command_line, comment='#', delimiter='%', /silent, format='A'
;	for i=0, n_elements(command_line)-1 do $
;	    d=execute(command_line[i])

	; Read the table with the list of fits files and initial guess values
	print, 'Reading input table:'
	READCOL, total_filelist, root_name_vector, velocity_initial_guess_vector,$
	         velocity_dispersion_initial_guess_vector, ellipticity_vector,$
		 position_angle_vector, fibre_number_vector, Reff_vector, mode_vector,$
		 /silent, format='A,F,F,F,F,F,F,A', comment='#'

	; Save the one requested by the user
	root_name = root_name_vector[input_number]
	mode = mode_vector[input_number]
	velocity_initial_guess = velocity_initial_guess_vector[input_number]
	velocity_dispersion_initial_guess = velocity_dispersion_initial_guess_vector[input_number]
	ell=ellipticity_vector[input_number]
	pa=position_angle_vector[input_number]
	number_of_fibres=fibre_number_vector[input_number]
	Reff=Reff_vector[input_number]

;	print, 'ROOT: '+root_name
;	print, 'MODE: '+mode
;	print, 'DATACUBE_DIR: '+datacube_dir
;	print, 'Guess V', velocity_initial_guess
;	print, 'Guess VDISP', velocity_dispersion_initial_guess
;	print, 'Guess ELL', ell
;	print, 'Guess PA', pa
;	print, 'IFU SIZE', number_of_fibres
;	print, 'Half-light', Reff

	; Set input file names and output directory
	;     make the directory if it doesn't exist
	; TODO: Allow for a single spectrum mode
	MDAP_SETUP_IO, root_name, output_root_dir, datacube_name, file_root, output_dir, $
		       output_file_root

	; TODO: These are not really used yet
	; Check if save_intermediate_steps exits, if not set to not save intermediate steps
	if n_elements(save_intermediate_steps) eq 0 then $
	    save_intermediate_steps = 0

	; Check if remove_null_elemens exits, if not set to remove null templates
	; TODO: Where are the removed from?
	if n_elements(remove_null_templates) eq 0 then $
	    remove_null_templates = 0

;	; Tell the user what's happening
;	print, ''
;	print, '# WORKING ON '+root_name+' ('+mode+')'
;	print, ''

	; TODO: Need to pass log_file_unit to other subroutines?

	; Instantiate the log file
	if ~keyword_set(nolog) then begin
	    MDAP_LOG_INSTANTIATE, output_dir, manga_dap_version, signifier, total_filelist, $
				  input_number, datacube_name, mode, file_root, $
				  velocity_initial_guess, velocity_dispersion_initial_guess, ell, $
				  pa, Reff, log_file_unit
	endif

	MDAP_DRP_CHECK_FILETYPE, datacube_name, mode

	; BLOCK 0 ----------------------------------------------------------------------
	;	Build execution plan
	;-------------------------------------------------------------------------------
	if ~keyword_set(quiet) then $
	    print, 'BLOCK 0 ... '
	mdap_execute_plan_version = 0.1

	; TODO: Same as ntpl_libraries, neml_files, and nabs_files in mdap_execution_setup

	n_tpl_lib = n_elements(template_libraries)	; Number of template libraries to use
	n_ems_par = n_elements(emission_line_parameters); Number of emission line parameter files
	n_abs_par = n_elements(absorption_line_parameters);Number of absorption line parameter files

	success = MDAP_CHECK_FILE(template_libraries, /search)	; Aborts of fails
	if success ne 0 then $
	    message, 'Template libraries not correctly defined.'
		
	success = MDAP_CHECK_FILE(emission_line_parameters)
	if success ne 0 then $
	    message, 'Emission-line parameter files not correctly defined.'

	success = MDAP_CHECK_FILE(absorption_line_parameters)
	if success ne 0 then $
	    message, 'Absorption-line parameter files not correctly defined.'

	MDAP_BUILD_EXECUTION_PLANS, n_tpl_lib, n_ems_par, n_abs_par, bin_type, w_range_sn, $
				    bin_par, threshold_ston_bin, bin_weight_by_sn2, $
				    w_range_analysis, threshold_ston_analysis, analysis, $
				    analysis_par, tpl_lib_analysis, ems_par_analysis, $
				    abs_par_analysis, overwrite_flag, execution_plan

	MDAP_GENERATE_OUTPUT_FILE_NAMES, output_file_root, execution_plan

	n_plans = n_elements(execution_plan)
	if ~keyword_set(quiet) then begin
	    for i=0,n_plans-1 do begin
		print, 'EXECUTION PLAN: ', i+1
		MDAP_PRINT_EXECUTION_PLAN, template_libraries, emission_line_parameters, $
					   absorption_line_parameters, execution_plan[i]
	    endfor
	endif

	SXADDPAR, header, 'VDAPPLAN', mdap_execute_plan_version, $
			  'mdap_build_execution_plans version'

	if ~keyword_set(nolog) then begin
	    printf, log_file_unit, '[INFO] Execution plans built using version: ', $
		    mdap_execute_plan_version
	    printf, log_file_unit, '[INFO] Number of plans to execute: ', n_plans
	endif

	; BLOCK 0.1 --------------------------------------------------------------------
	;	VERSION CONTROL
	;-------------------------------------------------------------------------------
	; TODO: Version control will eventually happen in
	; MDAP_GENERATE_OUTPUT_FILE_NAMES, or preferrably in a more well-named
	; procedure
	
	manga_dap_version_previous=0

	mdap_read_drp_fits_version='0'
	mdap_calculate_sn_version='0'
	mdap_spatial_binning_version='0'

	MDAP_READ_DRP_FITS, version= mdap_read_drp_fits_version
	mdap_fiducial_bin_xy_version = '0.1'
	mdap_tpl_lib_setup_version = '0.1'
	MDAP_CALCULATE_SN, version=mdap_calculate_sn_version
	MDAP_SPATIAL_BINNING, version=mdap_spatial_binning_version

	mdap_spatial_binning_version_previous=0

	mdap_spectral_fitting_version=0
	mdap_spectral_fitting_version_previous=0

	mdap_measure_indices_version=0
	mdap_measure_indices_version_previous=0

	mdap_spatial_radial_binning_version=0
	mdap_spatial_radial_binning_version_previous=0

	execute_all_modules = 1				; TODO: For now, just execute all modules!

	if ~keyword_set(quiet) then $
	    print, 'BLOCK 0 ... DONE.'
	; END BLOCK 0 ------------------------------------------------------------------

	if ~keyword_set(quiet) then $
	    print, 'BLOCK 1 ... '
	; BLOCK 1 ----------------------------------------------------------------------
	;	Read the DRP fits file
	;-------------------------------------------------------------------------------

;	if mdap_read_datacube_version gt mdap_read_datacube_version_previous $
;	   or execute_all_modules eq 1 then begin

	; TODO: Version required for this should depend on the version of DRP
	; output being read

	; TODO: Put in the conditions for when this block is executed:
	;	- If any of the plans have the overwrite flag flipped
	;	- If any of the plans do not have a previously created output
	;	  file with the binned spectra

	; TODO: These operations are independent of the plan ---------------------------
	    ; TODO: Allow to read a single spectrum and skip binning step
	    MDAP_READ_DRP_FITS, datacube_name, header, flux, ivar, mask, wave, sres, skyx, skyy, $
				type=mode, unit=unit
	    mask[*,*] = 0.			; TODO: Unmask everything for now

	    ; Get the spaxel size
	    MDAP_GET_SPAXEL_SIZE, header, spaxel_dx, spaxel_dy, type=mode, unit=unit
;	    print, spaxel_dx, spaxel_dy

	    ; Set a fiducial set of coordinates to use for each spectrum
	    MDAP_FIDUCIAL_BIN_XY, skyx, skyy, bskyx, bskyy

	    ; Add/update the version number for DAP read version
	    SXADDPAR, header, 'VDAPREAD', mdap_read_drp_fits_version, 'mdap_read_drp_fits version'
	    SXADDPAR, header, 'VDAPFXY', mdap_fiducial_bin_xy_version, $
		      'mdap_fiducial_bin_xy version'

	    ; Print information to the log file
	    if ~keyword_set(nolog) then begin
		printf, log_file_unit, '[INFO] Read DRP data using version: ' + $
			mdap_read_drp_fits_version
		printf, log_file_unit, '[INFO] Input data type: '+mode
		sz = size(flux)
		printf, log_file_unit, '[INFO] Total number of spectra: '+MDAP_STC(sz[1])
		printf, log_file_unit, '[INFO] Number of spectral channels: '+MDAP_STC(sz[2])
		printf, log_file_unit, '[INFO] Fiducial bin centers calcualted using version: ' + $
			mdap_fiducial_bin_xy_version
	    endif

	    ;-------------------------------------------------------------------
	    ; TODO: This block of code is essentially obsolete.  We expect to
	    ; always read the spectral resolution from the header of the
	    ; DRP-produced fits file.

	    ; Change the spectral resolution according to an input file
	    ; Columns must be:
	    ;	1. Wavelength in angstroms
	    ;	2. Resolution (lamda/delta lambda)
	    ;	3. delta lamga (FWHM) in angstroms
	    ;	4. delta lamga (FWHM) in km/s
	    if n_elements(instrumental_fwhm_file) ne 0 then begin
		if file_test(instrumental_fwhm_file) eq 0 then begin
		    print, 'Unable to read '+instrumental_fwhm_file+'!  Using fits extention or ' $
			   + 'default.'
		endif
		print, 'Reading '+instrumental_fwhm_file
		READCOL, instrumental_fwhm_file, ww, r, fwhm_ang, fwhm_kms, /silent
	    endif
	    
	    ; Initialize the instrumental resolution.
	    if n_elements(fwhm_ang) ne 0 then begin		; use data from input file
		if ~keyword_set(nolog) then begin
		    printf, log_file_unit, '[INFO] Spectral resolution changed using data in: '+ $
			    instrumental_fwhm_file
		endif
		sres=interpol(fwhm_ang,r,wave)		; Interpolate to the object wavelengths
	    endif else if n_elements(sres) eq 0 then begin	; sres not available from fits file
		if ~keyword_set(nolog) then begin
		    printf, log_file_unit, '[INFO] Spectral resolution unavailable, assume R=2000'
		endif
		sres = make_array(n_elements(wave), /double, value=2000.0d)	; Set R=2000.
	    endif
	    ;-------------------------------------------------------------------

;	    sz=size(mask)
;	    print, sz
;	    print, mask[0,0:100]

	;-------------------------------------------------------------------------------
	if ~keyword_set(quiet) then $
	    print, 'BLOCK 1 ... DONE.'
	; END BLOCK 1 ------------------------------------------------------------------

	; TODO: Once the DRP fits file has been read, the required input to
	;	appropriatly manipulate ALL the template libraries is available

	; BLOCK 2 ----------------------------------------------------------------------
	;	- Read the template library
	;	- Match the resolution of the templates to that of the galaxy
	;	  spectra (as best as possible).
	;	- Resample the templates to match the object sampling
	;-------------------------------------------------------------------------------

	if ~keyword_set(quiet) then $
	    print, 'BLOCK 2 ... '

	; Determine which template libraries are needed
	use_tpl_lib = intarr(n_tpl_lib)
	for i=0,n_plans-1 do begin
	    if execution_plan[i].tpl_lib ne -1 then $
		use_tpl_lib[execution_plan[i].tpl_lib] = 1
	endfor

	; Cycle through all the libraries, generating the template libraries to
	; use for the analysis, if they don't already exist!
	indx=where(use_tpl_lib eq 1)
	if indx[0] ne -1 then begin

	    for i=0,n_tpl_lib-1 do begin
		if use_tpl_lib[i] ne 1 then $
		    continue

		; TODO: Placeholder!  Need to replace with function that grabs
		;	the library keyword.  Not sure how to do this yet.
		;	Maybe have a predefined link between library_key and the
		;	files to search for?
		;library_key = MDAP_GET_TPL_LIB_KEY(template_libraries[i])
		library_key = tpl_library_keys[i]; 'M11-MARCS'
		tpl_out_fits = MDAP_SET_TPL_LIB_OUTPUT_FILE(output_file_root, library_key)

		; TODO: What are the flux units for the templates?  Does it
		; matter because I'm normalizing them?

		; TODO: Should the convolution conserve the integral of the template spectra?

		; The file exists, so just continue
		if file_test(tpl_out_fits) eq 1 then $
		    continue

		; Read the template library
		MDAP_READ_TEMPLATE_LIBRARY, library_key, template_libraries[i], tpl_flux, $
					    tpl_ivar, tpl_mask, tpl_wave, tpl_sres

		; Convert the wavelengths to vacuum if necessary
		; TODO: Does not account for the effect on the spectral resolution
		if tpl_vacuum_wave[i] eq 0 then $
		    AIRTOVAC, tpl_wave

		; TODO: Need to limit template wavelength range before running
		; the resolution matching procedure to avoid extrapolation
		; errors.  As with the matching, this trimming accounts for the
		; guess redshift of the galaxy.

		; Limit the spectral range of the templates to be the same as the galaxy spectra
		c=299792.458d				; Speed of light in km/s
		z_guess = velocity_initial_guess / c	; Guess redshift
		if ~keyword_set(quiet) then begin
		    print, 'Trimming templates to ' + MDAP_STC(min(wave/(1+z_guess))) + $
			   'A to ' + MDAP_STC(max(wave/(1+z_guess))) + 'A'
		endif

		MDAP_TRIM_SPECTRAL_RANGE, tpl_flux, tpl_ivar, tpl_mask, tpl_sres, tpl_wave, $
					  ([min(wave/(1+z_guess)), max(wave/(1+z_guess))])

		; Match the resolution of the templates to the galaxy data
		; accounting for the guess redshift of the galaxy.
		MDAP_MATCH_RESOLUTION_TPL2OBJ, tpl_flux, tpl_ivar, tpl_mask, tpl_wave, $
					       tpl_sres, wave/(1+z_guess), sres, tpl_soff, $
					       /no_offset

		; Get velocity scale of the galaxy data
		velScale=MDAP_VELOCITY_SCALE(wave, /log10)

		; Resample the templates to match the object spectra: On input,
		; tpl_wave expected to have the same size as tpl_flux
		; (wavelength defined separately for each spectrum); on output,
		; tpl_wave and tpl_sres are single vectors common to ALL
		; template spectra.  Sampling is forced to match galaxy
		; sampling, via velscale
		MDAP_RESAMPLE_TEMPLATES, tpl_flux, tpl_ivar, tpl_mask, tpl_wave, tpl_sres, $
					 velScale, /reform_sres

		; Normalize the templates and save the normalization
		MDAP_NORMALIZE_TEMPLATES, tpl_flux, tpl_ivar, tpl_mask, tpl_flux_norm

		; Write the template data to a fits file
		MDAP_WRITE_RESAMPLED_TEMPLATES, tpl_out_fits, library_key, z_guess, tpl_flux_norm, $
						mdap_tpl_lib_setup_version, tpl_wave, tpl_flux, $
						tpl_ivar, tpl_mask, tpl_sres, tpl_soff

		if ~keyword_set(nolog) then begin
		    printf, log_file_unit, '[INFO] Read template library: ' + library_key
		    printf, log_file_unit, '[INFO] Resolution and sampling matched to: ' + $
			    datacube_name
		    printf, log_file_unit,'[INFO] Template library setup done using version: ' + $
			    mdap_tpl_lib_setup_version
		    printf, log_file_unit, '[INFO] Results written to: ' + tpl_out_fits
		endif
	    endfor	; Loop over template libraries
	endif	; If any template libraries are used

	if ~keyword_set(quiet) then $
	    print, 'BLOCK 2 ... DONE.'
	; END BLOCK 2 ------------------------------------------------------------------

	;-------------------------------------------------------------------------------
	; ANALYSIS PREP ----------------------------------------------------------------

	; Initialize the extinction curve
	MW_extinction = 0.			; TODO: replaced with appropriate extension in the
						;       input file, or total_filelist

	; From Oliver Steele (04 Sep 2014)
	; 	Requires galRA, galDEC, equinox, lambda_obs, spectrumdata
	;GLACTC,galRA,galDEC,equinox,gl,gb,1,/DEGREE	; Converts RA/DEC to galactic coords 
	;mw_ebmv = DUST_GETVAL(gl,gb,/interp)		; Grabs E(B-V) from Schlegel dust maps
	;FM_UNRED, lambda_obs, spectrumdata, mw_ebmv	; De-reddens using ?? extinction curve

;	print, star_kin_starting_guesses		; h3 and h4 initialized to 0
;	print, gas_kin_starting_guesses

	;-------------------------------------------------------------------------------
	; Now can begin cycling through execution plans

	for i=0,n_plans-1 do begin

	    if ~keyword_set(quiet) then $
		print, 'Beginning ExecutionPlans: ', i+1

	    ; TODO: If previous file exists, read or overwrite it?

	    ;-------------------------------------------------------------------------------
	    ; S/N CALCULATION
	    ;-------------------------------------------------------------------------------

	    ; Determine the 'good' spectra based on a set of criteria defined by this procedure
	    ;	TODO: Currently gindx is not used.
	    ;	      gflag is used by MDAP_CALCULATE_SN and MDAP_SPATIAL_BINNING
	    MDAP_SELECT_GOOD_SPECTRA, flux, ivar, mask, gflag, gindx

	    ; Select the pixels to use in the S/N calculation
	    MDAP_SELECT_WAVE, wave, execution_plan[i].wave_range_sn, lam_sn

	    ; Calculate the S/N per pixel over some wavelength range
	    MDAP_CALCULATE_SN, flux, ivar, mask, wave, lam_sn, signal, noise, gflag=gflag

	    ; Write the version of the S/N calculate ot the header
	    SXADDPAR, header, 'VDAPSTON', mdap_calculate_sn_version, 'mdap_calculate_sn version'

	    ; Add some information to the log
	    if ~keyword_set(nolog) then begin
		printf, log_file_unit, '[INFO] S/N calculation over range: ', $
			execution_plan[i].wave_range_sn
		printf, log_file_unit, '[INFO] S/N calculation done using version: ' + $
			mdap_calculate_sn_version
	    endif

	    ; Write basics of the DRP fits files
	    ; TODO: Can free memory for variables if they won't be used again.  Do this?
	    MDAP_WRITE_OUTPUT, execution_plan[i].ofile, header=header, dx=spaxel_dx, dy=spaxel_dy, $
			       w_range_sn=execution_plan[i].wave_range_sn, xpos=bskyx, ypos=bskyy, $
			       signal=signal, noise=noise

	    ;-------------------------------------------------------------------------------
	    ;-------------------------------------------------------------------------------


	    ;-------------------------------------------------------------------------------
	    ; Spatial binning
	    ;-------------------------------------------------------------------------------

	    ; TODO: If previous file exists, read or overwrite it?

	    MDAP_SPATIAL_BINNING, flux, ivar, mask, signal, noise, gflag, bskyx, bskyy, spaxel_dx, $
				  spaxel_dy, execution_plan[i].bin_type, $
				  execution_plan[i].bin_par, execution_plan[i].threshold_ston_bin, $
				  execution_plan[i].bin_weight_by_sn2, bin_wgts, bin_indx, $
				  bin_flux, bin_ivar, bin_mask, xbin, ybin, bin_area, bin_ston, $
				  nbin, sn_calibration=sn_calibration;, /plot

	    ; Write the version of the spatial binning to the header
	    SXADDPAR, header, 'VDAPBIN', mdap_spatial_binning_version, $
		      'mdap_spatial_binning version'

	    ; Add some information to the log
	    if ~keyword_set(nolog) then begin
		printf, log_file_unit, '[INFO] Type of spatial binning: ', $
			execution_plan[i].bin_type
		printf, log_file_unit, '[INFO] Spatial bining version: ' + $
			mdap_spatial_binning_version
		if execution_plan[i].bin_type ne 'NONE' then begin
		    if n_elements(execution_plan[i].bin_par) then $
			printf, log_file_unit, '[INFO] Bin parameter: ', execution_plan[i].bin_par
		    printf, log_file_unit, '[INFO] Bin threshold S/N: ', $
			    execution_plan[i].threshold_ston_bin
		    printf, log_file_unit, '[INFO] Weight by S/N^2: ', $
			    execution_plan[i].bin_weight_by_sn2
		endif
		printf, log_file_unit, '[INFO] Number of bins: ', n_elements(nbin)
	    endif

	    ; Write the binning results
	    MDAP_WRITE_OUTPUT, execution_plan[i].ofile, header=header, $
			       bin_type=execution_plan[i].bin_type, $
			       bin_par=execution_plan[i].bin_par, $
			       threshold_ston_bin=execution_plan[i].threshold_ston_bin, $
			       bin_indx=bin_indx, bin_weights=bin_wgts, wave=wave, sres=sres, $
			       bin_flux=bin_flux, bin_ivar=bin_ivar, bin_mask=bin_mask, $
			       xbin=xbin, ybin=ybin, bin_area=bin_area, bin_ston=bin_ston, $
			       bin_n=nbin, /read_header

	    ;-------------------------------------------------------------------------------
	    ;-------------------------------------------------------------------------------

	    ;-------------------------------------------------------------------------------
	    ; Check if any analysis was requested, if not continue with the loop
	    ; TODO: UNCOMMENT WHEN LOOP IS INSTATED
;	    indx = where(execution_plan[i].analysis eq 1)
;	    if indx[0] eq -1 then $		; No analysis to perform
;	    	continue

	    ;-------------------------------------------------------------------------------
	    ; Read the pre-existing template library
	    ;	The template library is required for ANY analysis
	    ; TODO: Should not be able to reach here without first creating tpl_out_fits!
	    library_key = tpl_library_keys[execution_plan[i].tpl_lib]; 'M11-MARCS'
	    tpl_out_fits = MDAP_SET_TPL_LIB_OUTPUT_FILE(output_file_root, library_key)
	   
	    MDAP_READ_RESAMPLED_TEMPLATES, tpl_out_fits, tpl_wave, tpl_flux, tpl_ivar, tpl_mask, $
					   tpl_sres, tpl_soff

	    ;-------------------------------------------------------------------------------
	    ; Read the emission line file; this is done at every iteration
	    ;   because it's fairly quick.
	    ; TODO: Check that this I/O doesn't slow things down too much!
	    if execution_plan[i].ems_par ne -1 then begin
		eml_par = MDAP_READ_EMISSION_LINE_PARAMETERS(emission_line_parameters[ $
							     execution_plan[i].ems_par ])
		if ems_vacuum_wave[execution_plan[i].ems_par] eq 0 then begin
		    neml = n_elements(eml_par)
		    for j=0,neml-1 do begin
			AIRTOVAC, eml_par[j].lambda, lam_vac
			eml_par[j].lambda = lam_vac
		    endfor
		endif
	    endif else $
		MDAP_ERASEVAR, eml_par	; Make sure it doesn't exist if no ems_par is defined

	    ; ##############################################################################
	    ; BEGIN computation of "model-independent" data products #######################

	    ; BLOCK 4 ----------------------------------------------------------------------
	    ; Run a full spectral fit including stellar population, gas and stellar
	    ; kinematics
	    ;-------------------------------------------------------------------------------

	    ; Set the starting guesses for the kinematics
	    ; TODO: get the guesses from a cross-correlation with a guess template?

	    ; TODO: convert this to a procedure that could also read in guesses
	    ; from a previous DAP-produced file
	    star_kin_starting_guesses = dblarr(n_elements(xbin),4)
	    star_kin_starting_guesses[*,0] = velocity_initial_guess		; stellar velocity
	    star_kin_starting_guesses[*,1] = velocity_dispersion_initial_guess	; stellar sigma
	    gas_kin_starting_guesses = star_kin_starting_guesses[*,0:1]		; gas velocity
	    gas_kin_starting_guesses[*,1]=50.					; gas sigma

	    ; Ignore any lines outside of the wavelength range of the galaxy data
	    ; TODO: Do this for each spectrum to allow for slightly different velocities?
	    ; TODO: Now done IN MDAP_SPECTRAL_BINNING
;	    MDAP_CHECK_EMISSION_LINES, eml_par, wave, velocity=velocity_initial_guess[0]

	    ; Perform the spectral fit
;	    print, spectra_fittin_parameters_patial_binning_1

	    ppxf_only=0
	    if execution_plan[i].analysis[0] eq 1 and execution_plan[i].analysis[1] eq 0 then $
		ppxf_only=1		; Only execute the PPXF step

	    ; TODO: Be specific about what this execution of MDAP_SPECTRAL_FITTING does;
	    ;       I assume it is run multiple times using the different S/N limits for
	    ;       different purposes
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
				   analysis_par=execution_plan[i].analysis_par, $
				   star_kin_starting_guesses=star_kin_starting_guesses, $
				   gas_kin_starting_guesses=gas_kin_starting_guesses, $
				   eml_par=eml_par, external_library=external_library, $
				   wave_range_analysis=execution_plan[i].wave_range_analysis, $
				   ppxf_only=ppxf_only, quiet=quiet

	    ; TODO: Add the spectral fitting version to the header of the output file
	    ; TODO: Add information to the log file

	    ; Write the analysis wavelength range to the header
	    MDAP_WRITE_OUTPUT, execution_plan[i].ofile, header=header, $
			       w_range_analysis=execution_plan[i].wave_range_analysis,
			       threshold_ston_analysis=execution_plan[i].threshold_ston_analysis,
			       /read_header
	    
	    ; Write the emission line parameters
	    if (execution_plan[i].analysis[0] eq 1 and n_elements(eml_par) ne 0 ) or $
		execution_plan[i].analysis[1] eq 1 then begin
		MDAP_WRITE_OUTPUT, execution_plan[i].ofile, eml_par=eml_par
	    endif

	    ; Write the stellar kinematics results; will always be done if
	    ; program makes it here
	    MDAP_WRITE_OUTPUT, execution_plan[i].ofile, header=header, $
			       obj_fit_mask_ppxf=obj_fit_mask_ppxf, weights_ppxf=weights_ppxf, $
			       add_poly_coeff_ppxf=add_poly_coeff_ppxf, $
			       mult_poly_coeff_ppxf=mult_poly_coeff_ppxf, $
			       bestfit_ppxf=bestfit_ppxf, chi2_ppxf=chi2_ppxf, $
			       stellar_kinematics_fit=stellar_kinematics, $
			       stellar_kinematics_err=stellar_kinematics_err, $
			       analysis_par=execution_plan[i].analysis_par, /read_header

	    ; Write the results of the emission-line fitting
	    if execution_plan[i].analysis[1] eq 1 then begin
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


	    ; Write the optimal templates and optimal template convolved with
	    ; the LOSVD; see MDAP_SPECTRAL_FITTING on how these are generated.
	    MDAP_WRITE_OUTPUT, execution_plan[i].ofile, optimal_template=best_template, $
			       losvd_optimal_template=best_template_losvd_conv

	if ~keyword_set(nolog) then $
	    MDAP_LOG_CLOSE, log_file_unit
	return
	endfor


	   ; need to save all flux maps (warning: sgandalf computes intensities, not fluxes)

	   ; Convert intensities and fluxes to surface brightness
	   junk = size(emission_line_fluxes_tpl)
	   for i = 0, junk[1]-1 do begin
		emission_line_intens_tpl[i,*] = emission_line_intens_tpl[i,*] / area_tpl[i]
		emission_line_intens_tpl_err[i,*] = emission_line_intens_tpl_err[i,*] $
						    / sqrt(area_tpl[i])
		emission_line_fluxes_tpl[i,*] = emission_line_fluxes_tpl[i,*] / area_tpl[i]
		emission_line_fluxes_tpl_err[i,*] = emission_line_fluxes_tpl_err[i,*] $
						    / sqrt(area_tpl[i])
	   endfor

;   for i = 0, junk[1]-1 do emission_line_intens_tpl[i,*]=emission_line_intens_tpl[i,*]/area_tpl[i]
;   for i = 0, junk[1]-1 do emission_line_intens_tpl_err[i,*]=emission_line_intens_tpl_err[i,*]/sqrt(area_tpl[i])
;   for i = 0, junk[1]-1 do emission_line_fluxes_tpl[i,*]=emission_line_fluxes_tpl[i,*]/area_tpl[i]
;   for i = 0, junk[1]-1 do emission_line_fluxes_tpl_err[i,*]=emission_line_fluxes_tpl_err[i,*]/sqrt(area_tpl[i])

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
	    bin_sn_tpl_real = ston_tpl
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
				       emission_line_file=emission_line_file_binning_2, $
				       extra_inputs=spectral_fit_par_bin2_2mns, $
				       mask_range=mask_range, $
				       external_library=external_library, quiet=quiet
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
				   emission_line_file=emission_line_file_binning_2, $
				   extra_inputs=spectral_fit_par_bin2, $
				   mask_range=mask_range, external_library=external_library, /quiet

	    ; save all flux maps, converting to surface brightness
	    ; TODO: Is this really surface brightness?
	    junk = size(emission_line_fluxes_str)
	    for i = 0, junk[1]-1 do begin
		emission_line_intens_str[i,*] = emission_line_intens_str[i,*] / area_str[i]
		emission_line_intens_str_err[i,*] = emission_line_intens_str_err[i,*] $
						    / sqrt(area_str[i])
		emission_line_fluxes_str[i,*] = emission_line_fluxes_str[i,*] / area_str[i]
		emission_line_fluxes_str_err[i,*] = emission_line_fluxes_str_err[i,*] $
						    / sqrt(area_str[i])
	    endfor

	    ;LC: calculate the real S/N of binned spectra
	    ; TODO: Is this really the "real S/N"?
	    bin_sn_str_real = ston_str
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
				   mask_range=mask_range, external_library=external_library, /quiet

	    ; save all flux maps, converting to surface brightness
	    junk = size(emission_line_fluxes_ems)
	    for i = 0, junk[1]-1 do begin
		emission_line_intens_ems[i,*] = emission_line_intens_ems[i,*] / area_ems[i]
		emission_line_intens_ems_err[i,*] = emission_line_intens_ems_err[i,*] $
						    / sqrt(area_ems[i])
		emission_line_fluxes_ems[i,*] = emission_line_fluxes_ems[i,*] / area_ems[i]
		emission_line_fluxes_ems_err[i,*] = emission_line_fluxes_ems_err[i,*] $
						    / sqrt(area_ems[i])
	    endfor

	    ; Re-execute all the remaining modules, if that wasn't already set
	    execute_all_modules = 1

	    ;LC: calculate the real S/N of binned spectra
	    ; TODO: Is this really the "real S/N"?
	    bin_sn_ems_real = ston_ems
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

;	; Initialize the instrumental resolution.  TODO: use fits extension
;	if n_elements(fwhm_ang) ne 0 then begin			; use input file
;	    fwhm_instr=interpol(fwhm_ang,ww,wavelength_output_tpl) 
;	endif else begin					; set to a constant
;	    fwhm_instr=wavelength_output_tpl*0.+2.54		; TODO: Don't like this practice
;	endelse

	;resolution_x_lick=[4000.,4400.,4900.,5400.,6000.]   ;CONTROLLARE
	;lick_fwhm_y=[11.5,9.2,8.4,8.4,9.8]                     ;CONTROLLARE
	;lick_resolution_tmp=interpol(lick_fwhm_y,resolution_x_lick,wavelength_output_tpl)
	;rrr=poly_fit(wavelength_output_tpl,lick_resolution_tmp,4,yfit=lick_resolution)

	;miles_resolution = fwhm_instr*0.+2.54

	; Determine the difference in instrumental resolution
	fwhm_instr=wavelength_output_tpl*0.+2.54	; TODO: Place holder
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

;	TODO: Omitted this because mdap_spatial_binning no longer produces the
;	binning map, which is used below in mdap_spatial_radial_binning.  Need
;	to come back and fix this.

;	; Create error array using fit residuals, if required
;	if mdap_spatial_radial_binning_version gt mdap_spatial_radial_binning_version_previous $
;	   or execute_all_modules eq 1 then begin
;	   
;	    ; TODO: Why is mdap_get_erro_from_residual not used? Function obsolete?
;	    ;mdap_get_error_from_residual,residuals_ems,galaxy_minus_ems_fit_model_ems,input_errors
;
;	    sz_ems = size(log_err_ems)
;	    nl = n_elements(wavelength_output_rest_frame_log)
;	    input_errors = fltarr(sz_ems[2], nl)
;	    log_step_gal = LOG_WAV_EMS[1]-LOG_WAV_EMS[0]
;	    for i = 0, sz_ems[2]-1 do begin
;		rf_gal_lam = exp(log_wav_ems-stellar_kinematics_ems[i,0] / velscale*(log_step_gal))
;		input_errors[i,*] = interpol(log_err_ems[*,i], rf_gal_lam, $
;					     wavelength_output_rest_frame_log)
;	    endfor
;
;	    ; Perform the radial binning
;	    MDAP_SPATIAL_RADIAL_BINNING, bin_sn_ems_real, x2d_reconstructed, y2d_reconstructed, $
;					 spatial_binning_ems, xbin_ems, ybin_ems, ell, pa, $
;					 galaxy_minus_ems_fit_model_ems, input_errors, $
;					 wavelength_output_rest_frame_log, spatial_binning_rad, $
;					 r_bin, r_bin_lo, r_bin_up, r2d_bin, r2d_bin_lo, $
;					 r2d_bin_up, radially_binned_spectra, $
;					 radially_binned_errors, $
;					 output_lrange=trim_wav_range_radial_binning, $
;					 output_wav=output_wav, $
;					 n_elements_bin=nelements_within_bin_radial, $
;					 low_radial_bins_user_inputs=low_radial_bins_user_inputs, $
;					 upper_radial_bins_user_inputs=upper_radial_bins_user_inputs, $
;					 Reff_=Reff, PSFsize_=PSFsize, $
;					 add_default_bins=add_default_bins
;
;	    ; Print progress
;	    printf, 1, '[INFO] datacube '+root_name+' radial binning: ' $
;		       +MDAP_STC(n_elements(r_bin),/integer)+' bins'
;
;
;	    ; Perform the spectral fit of the radially binned data
;	    star_kin_starting_guesses_rbin = fltarr(n_elements(r_bin),4)
;	    gas_kin_starting_guesses_rbin  = fltarr(n_elements(r_bin),2)
;	    star_kin_starting_guesses_rbin[*,1] = velocity_dispersion_initial_guess[0]
;	    gas_kin_starting_guesses_rbin[*,1]=50.
;	    
;	    loglam_gal_rbin=alog(output_wav)
;	    
;	    MDAP_SPECTRAL_FITTING, radially_binned_spectra, radially_binned_errors, $
;				   loglam_gal_rbin, library_log_ems, log_wav_library_ems, $
;				   velscale, stellar_kinematics_rbin, stellar_kinematics_rbin_err, $
;				   stellar_weights_rbin, emission_line_kinematics_rbin, $
;				   emission_line_kinematics_rbin_err, $
;				   emission_line_kinematics_rbin_individual, $
;				   emission_line_kinematics_rbin_individual_err, $
;				   emission_line_intens_rbin, emission_line_intens_rbin_err, $
;				   emission_line_fluxes_rbin, emission_line_fluxes_rbin_err, $
;				   emission_line_equivW_rbin, emission_line_equivW_rbin_err, $
;				   wavelength_input=output_wav, wavelength_output_rbin, $
;				   best_fit_model_rbin, galaxy_minus_ems_fit_model_rbin, $
;				   best_template_rbin, best_template_LOSVD_conv_rbin, $
;				   reddening_rbin, reddening_rbin_err, residuals_rbin, $
;				   star_kin_starting_guesses=star_kin_starting_guesses_rbin, $
;				   gas_kin_starting_guesses=gas_kin_starting_guesses_rbin, $
;				   emission_line_file=emission_line_file_rad_binning, $
;				   extra_inputs=spectral_fit_par_rbin, $
;				   fwhm_instr_kmsec_matrix=fwhm_instr_kmsec_matrix/3., $
;				   range_v_star=[-50.,50.], range_v_gas=[-50.,50.], $
;				   mask_range=mask_range, external_library=external_library, /quiet
;
;	    ; Print progress
;	    printf,1,'[INFO] datacube '+root_name+' radial binning: spectral fitting' 
;
;	    ;LC: calculate the real S/N of binned spectra
;	    ; TODO: Is this really the "real S/N"?
;	    bin_sn_rad_real = r_bin
;	    indxx = where(wavelength_output_rbin ge w_range_for_sn_computation[0] $
;			  and wavelength_output_rbin le w_range_for_sn_computation[1] )
;
;	    for i = 0, n_elements(bin_sn_rad_real) -1 do begin
;		MDAP_CALCULATE_SPECTRUM_SN, best_fit_model_rbin[i,indxx], residuals_rbin[i,indxx], $
;					    wavelength_output_rbin[indxx], sn_per_angstorm, /rms
;		bin_sn_rad_real[i] = sn_per_angstorm[0]
;	    endfor
;
;	    ;junk = temporary(sn_per_angstorm)
;
;	    ; Measure the absorption-line indices using the radially binned spectra
;
;	    ;lick_resolution_tmp=interpol(lick_fwhm_y,resolution_x_lick,wavelength_output_rbin)
;	    ;rrr=poly_fit(wavelength_output_rbin,lick_resolution_tmp,4,yfit=lick_resolution)
;	    ;fwhm_diff_indices=sqrt(double(lick_resolution)^2.-double(fwhm_instr)^2.)*0.
;
;	    ; Initialize the instrumental resolution.  TODO: use fits extension
;	    if n_elements(fwhm_ang) ne 0 then begin			; use input file
;		fwhm_instr=interpol(fwhm_ang,ww,wavelength_output_rbin) 
;	    endif else begin						; set to a constant
;		fwhm_instr=wavelength_output_rbin*0.+2.54		; TODO: bad practice?
;	    endelse
;
;	    ;miles_resolution = fwhm_instr*0.+2.54
;	    ;fwhm_diff_indices=sqrt(miles_resolution^2-fwhm_instr^2) ;fwhm in angstrom
;
;	    ; Determine the difference in instrumental resolution
;	    lick_resolution = fwhm_instr*0.+8.4
;	    fwhm_diff_indices=sqrt(lick_resolution^2 - fwhm_instr^2)
;
;	    ; Check where the Lick instrumental resolution is larger than the data resolution
;	    indici = where(finite(fwhm_diff_indices) ne 1)
;	    if indici[0] ne -1 then fwhm_diff_indices[indici] = 0.  ; TODO: Forced to be the same
;
;	    ; Measure the absorption-line indices
;	    MDAP_MEASURE_INDICES, absorption_line_indices, wavelength_output_rbin, $
;				  galaxy_minus_ems_fit_model_rbin, best_template_rbin, $
;				  best_template_LOSVD_conv_rbin, stellar_kinematics_rbin[*,0], $
;				  residuals_rbin, fwhm_diff_indices, abs_line_indices_rbin, $
;				  abs_line_indices_errors_rbin, abs_line_indices_template_rbin, $
;				  abs_line_indices_template_losvd_rbin, dir=output_dir+'rbin_', $
;				  /noplot
;
;	    ; Print progress
;	    printf,1,'[INFO] datacube '+root_name+' radial binning: measured indices'
;
;	    ; Add the version to the header
;	    SXADDPAR, header_2d, 'BLOCK6', mdap_spatial_radial_binning_version, $
;		      'mdap_spatial_radial_binning_version'
;
;	    ; Re-execute all the remaining modules, if that wasn't already set
;	    execute_all_modules=1
;	endif
;
;	; Save progress
;	if save_intermediate_steps eq 1 then begin
;	    save, filename=root_name+mode+'_block6a.idl', /variables
;	endif

	; TODO: Here, I continue with the rest of block 6

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
;	MDAP_ADD_FITS_LAYER, output_filefits, spatial_binning_tpl, k, 'EXTNAME','Binning_map_1'
;	k=k+1

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
	p1.area_bin = area_tpl
	; p1.reddening = reddening_tpl
	p1.StoN = bin_sn_tpl_real
	p1.Nelements=nbin_tpl

	for i = 0, NINDICES-1 do begin
	    d=execute('p1.'+indices[i].name+'=ABS_LINE_INDICES[*,i]')
	    d=execute('p1.'+indices[i].name+'_err=ABS_LINE_INDICES_ERRORS[*,i]')
	endfor
	
	MWRFITS, p1, output_filefits, /silent			;?
	h1=HEADFITS(output_filefits,EXTEN = k)			;?
	SXADDPAR, h1, 'EXTNAME', 'Binning_1_data'		;Add the extension name to header
	MODFITS, output_filefits, 0, h1, exten_no=k		;?

	; Add extension with stellar kinematics
	; TODO: I removed this layer
	k=k+1
;	MDAP_ADD_FITS_LAYER, output_filefits, spatial_binning_str, k, 'EXTNAME', 'Binning_map_2'
;	k=k+1

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
	p2.area_bin = area_str
	p2.StoN = bin_sn_str_real
	p2.Nelements=nbin_str
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
	; TODO: I removed this
;	MDAP_ADD_FITS_LAYER, output_filefits, spatial_binning_ems, k, 'EXTNAME', 'Binning_map_3'
;	k=k+1

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
	p3.area_bin = area_ems
	p3.StoN = bin_sn_ems_real
	p3.Nelements=nbin_ems
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


;	;*** GET MODULES VERSION AND CHECK PREVIOUS ANALYSIS ******************************
;	check_previous_analysis = file_test(output_filefits)	; Check for output fits file
;	check_previous_session = file_test(output_idlsession)	; Check for idl session file
;
;	print, check_previous_analysis
;	print, check_previous_session
;
;	; if BOTH exist, read the fits file and determine the progress using the header
;	if check_previous_analysis eq 1 and check_previous_session eq 1 $
;	   and keyword_set(check_version) then begin
;	    tmp_header = HEADFITS(output_filefits)		; Grab the header
;	    manga_dap_version_previous = strcompress(sxpar(tmp_header,'DAPVER'),/remove_all)
;
;	    mdap_read_datacube_version=strcompress(sxpar(tmp_header,'BLOCK1'),/remove_all)
;	    mdap_read_datacube_version_previous = mdap_read_datacube_version
;
;	    mdap_spatial_binning_version=strcompress(sxpar(tmp_header,'BLOCK2'),/remove_all)
;	    mdap_spatial_binning_version_previous = mdap_spatial_binning_version
;
;	    mdap_log_rebin_version=strcompress(sxpar(tmp_header,'BLOCK3'),/remove_all)
;	    mdap_log_rebin_version_previous=mdap_log_rebin_version
;
;	    mdap_spectral_fitting_version=strcompress(sxpar(tmp_header,'BLOCK4'),/remove_all)
;	    mdap_spectral_fitting_version_previous=mdap_spectral_fitting_version
;
;	    mdap_measure_indices_version=strcompress(sxpar(tmp_header,'BLOCK5'),/remove_all)
;	    mdap_measure_indices_version_previous=mdap_measure_indices_version
;
;	    mdap_spatial_radial_binning_version=strcompress(sxpar(tmp_header,'BLOCK6'),/remove_all)
;	    mdap_spatial_radial_binning_version_previous=mdap_spatial_radial_binning_version
;	endif else begin					; Set all versions to '0'
;	    manga_dap_version_previous = '0'
;
;	    mdap_read_datacube_version='0'
;	    mdap_read_datacube_version_previous='0'
;
;	    mdap_spatial_binning_version='0'
;	    mdap_spatial_binning_version_previous='0'
;
;	    mdap_log_rebin_version='0'
;	    mdap_log_rebin_version_previous='0'
;
;	    mdap_spectral_fitting_version='0'
;	    mdap_spectral_fitting_version_previous='0'
;
;	    mdap_measure_indices_version='0'
;	    mdap_measure_indices_version_previous='0'
;
;	    mdap_spatial_radial_binning_version='0'
;	    mdap_spatial_radial_binning_version_previous='0'
;	endelse
;
;	print, manga_dap_version_previous
;
;	print, mdap_read_datacube_version
;	print, mdap_read_datacube_version_previous
;
;	print, mdap_spatial_binning_version
;	print, mdap_spatial_binning_version_previous
;
;	print, mdap_log_rebin_version
;	print, mdap_log_rebin_version_previous
;
;	print, mdap_spectral_fitting_version
;	print, mdap_spectral_fitting_version_previous
;
;	print, mdap_measure_indices_version
;	print, mdap_measure_indices_version_previous
;
;	print, mdap_spatial_radial_binning_version
;	print, mdap_spatial_radial_binning_version_previous
;
;	; update the versions of each procedure to check if these procedures should be redone
;;	MDAP_READ_DATACUBE,version=mdap_read_datacube_version	; TODO: OBSOLETE!
;	mdap_read_datacube_version=0.3
;	MDAP_READ_DRP_FITS,version=mdap_read_datacube_version	; TODO: Not all of block 1!
;	MDAP_SPATIAL_BINNING,version=mdap_spatial_binning_version
;;	MDAP_LOG_REBIN,version=mdap_log_rebin_version
;	mdap_log_rebin_version=0.2
;	MDAP_SPECTRAL_FITTING,version=mdap_spectral_fitting_version
;	MDAP_MEASURE_INDICES,version=mdap_measure_indices_version
;	MDAP_SPATIAL_RADIAL_BINNING,version=mdap_spatial_radial_binning_version
;
;	print, mdap_read_datacube_version
;	print, mdap_spatial_binning_version
;	print, mdap_log_rebin_version
;	print, mdap_spectral_fitting_version
;	print, mdap_measure_indices_version
;	print, mdap_spatial_radial_binning_version







;	; TODO: Decide how much of this is neccesary
;	;	Then put it in a procedure
;	
;	; Add the configuration data to the header
;	SXADDPAR, header, 'FITSLST', total_filelist, 'Input file with data files and parameters', $
;		  AFTER='DAPVER'
;	SXADDPAR, header, 'OROOT', output_root_dir, 'Root name of output directory'
;
;	SXADDPAR, header, 'SNRANGE', MDAP_ARRAY_STRING(w_range_for_sn_computation), $
;		          'Wavelength range for S/N computation'
;
;	SXADDPAR, header, 'TRIM1', MDAP_ARRAY_STRING(trim_wav_range_spatial_binning_1), $
;			  'Wavelength range considered in spectral fitting analysis'
;
;	SXADDPAR, header, 'TRIM2', MDAP_ARRAY_STRING(trim_wav_range_spatial_binning_2), $
;			  'Wavelength range considered in spectral fitting analysis'
;
;	SXADDPAR, header, 'TRIM3', MDAP_ARRAY_STRING(trim_wav_range_spatial_binning_3), $
;			  'Wavelength range considered in spectral fitting analysis'
;
;	SXADDPAR, header, 'TRIMR', MDAP_ARRAY_STRING(trim_wav_range_radial_binning), $
;			  'Wavelength range considered in spectral fitting analysis'
;	
;	SXADDPAR, header, 'VSCALE', MDAP_STC(velscale), 'Pixel scale in km/s'
;	
;	SXADDPAR, header, 'TPLIB1', stellar_library_spatial_binning_1, $
;			  'Stellar template library used in spectral fitting analysis'
;
;	SXADDPAR, header, 'TPLIB2', stellar_library_spatial_binning_2, $
;			  'Stellar template library used in spectral fitting analysis'
;
;	SXADDPAR, header, 'TPLIB3', stellar_library_spatial_binning_3, $
;			  'Stellar template library used in spectral fitting analysis'
;
;	SXADDPAR, header, 'RSS_SN1', MDAP_STC(sn1_rss), 'S/N level for RSS binning'
;	SXADDPAR, header, 'RSS_SNT1', MDAP_STC(sn_thr_tpl_rss), 'S/N threshold for RSS binning'
;	SXADDPAR, header, 'RSS_SN2', MDAP_STC(sn2_rss), 'S/N level for RSS binning'
;	SXADDPAR, header, 'RSS_SNT2', MDAP_STC(sn_thr_str_rss), 'S/N threshold for RSS binning'
;	SXADDPAR, header, 'RSS_SN3', MDAP_STC(sn3_rss), 'S/N level for RSS binning'
;	SXADDPAR, header, 'RSS_SNT3', MDAP_STC(sn_thr_ems_rss), 'S/N threshold for RSS binning'
;
;	SXADDPAR, header, 'CUBE_SN1', MDAP_STC(sn1_datacubes), 'S/N level for CUBE binning'
;	SXADDPAR, header, 'CUBE_SNT1', MDAP_STC(sn_thr_tpl_datacubes), $
;			  'S/N threshold for CUBE binning'
;	SXADDPAR, header, 'CUBE_SN2', MDAP_STC(sn2_datacubes), 'S/N level for CUBE binning'
;	SXADDPAR, header, 'CUBE_SNT2', MDAP_STC(sn_thr_str_datacubes), $
;			  'S/N threshold for CUBE binning'
;	SXADDPAR, header, 'CUBE_SN3', MDAP_STC(sn3_datacubes), 'S/N level for CUBE binning'
;	SXADDPAR, header, 'CUBE_SNT3', MDAP_STC(sn_thr_ems_datacubes), $
;			  'S/N threshold for CUBE binning'
;
;	SXADDPAR, header, 'SNWGT', MDAP_STC(weight_for_sn), 'S/N weighting method used'
;
;	SXADDPAR, header, 'EMLINE1', emission_line_file_spatial_binnin_1, $
;			  'File with definition of emission lines and their setup'
;	SXADDPAR, header, 'EMLINE2', emission_line_file_spatial_binnin_2, $
;			  'File with definition of emission lines and their setup'
;	SXADDPAR, header, 'EMLINE3', emission_line_file_spatial_binnin_3, $
;			  'File with definition of emission lines and their setup'
;	SXADDPAR, header, 'EMLINER', emission_line_file_radial_binning, $
;			  'File with definition of emission lines and their setup'
;	
;	SXADDPAR, header, 'ABSLINE', absorption_line_indices, $
;			  'File with definition of absorption line indices'
;
;	SXADDPAR, header, 'SSTEPS', MDAP_STC(save_intermediate_steps), $
;			  'Intermediate steps are saved'
;
;	SXADDPAR, header, 'RMNTPL', MDAP_STC(remove_null_templates), $
;			  'Remove templates with null weights during spectral fitting'
;
;;	SXADDPAR, header, 'EXTERN', external_library, 'Directory with external F90/C code'
;
;	SXADDPAR, header, 'SFPAR1', MDAP_ARRAY_STRING(spectra_fittin_parameters_patial_binning_1), $
;			  'Spectral fitting parameters'
;	SXADDPAR, header, 'SFPAR2_2m', $
;		  MDAP_ARRAY_STRING(spectra_fittin_parameters_patial_binning_2_two_moments), $
;		  'Spectral fitting parameters'
;	SXADDPAR, header, 'SFPAR2', MDAP_ARRAY_STRING(spectra_fittin_parameters_patial_binning_2), $
;			  'Spectral fitting parameters'
;	SXADDPAR, header, 'SFPAR3', MDAP_ARRAY_STRING(spectra_fittin_parameters_patial_binning_3), $
;			  'Spectral fitting parameters'
;	SXADDPAR, header, 'SFPARR', $
;		  MDAP_ARRAY_STRING(spectra_fittin_parameters_patial_binning_readial), $
;		  'Spectral fitting parameters'
;



















