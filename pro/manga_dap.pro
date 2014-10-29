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
;	MANGA_DAP, input_number, /nolog, /quiet, /plot, /dbg
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
;	/plot
;		Produce the Voronoi binning, PPXF and GANDALF plots
;
;	/dbg
;		Only attempts to fit the first spectrum of the first
;		execution plan as a test run mostly used for debugging
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
;	MDAP_CHECK_FILE_EXISTS()
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
;-------------------------------------------------------------------------------

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
	input_number, nolog=nolog, quiet=quiet, plot=plot, dbg=dbg

	manga_dap_version = '0.9'	; set the version number
	t0=systime(/seconds)/60./60.	; start time

	; BLOCK 0 --------------------------------------------------------------
	; Initialization
	;	- File I/O setup
	;	- Execution plan setup
	;	- Version control
	;-----------------------------------------------------------------------
	if ~keyword_set(quiet) then $
	    print, 'BLOCK 0 ... '
	; TODO: Put all this initialization stuff in a wrapper?

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

	; TODO: Need to pass log_file_unit to other subroutines?

	; Instantiate the log file
	if ~keyword_set(nolog) then begin
	    MDAP_LOG_INSTANTIATE, output_dir, manga_dap_version, signifier, total_filelist, $
				  input_number, datacube_name, mode, file_root, $
				  velocity_initial_guess, velocity_dispersion_initial_guess, ell, $
				  pa, Reff, log_file_unit
	endif

	MDAP_DRP_CHECK_FILETYPE, datacube_name, mode

	;-----------------------------------------------------------------------
	mdap_execute_plan_version = 0.1

	; TODO: Same as ntpl_libraries, neml_files, and nabs_files in mdap_execution_setup

	n_tpl_lib = n_elements(template_libraries)	; Number of template libraries to use
	n_ems_par = n_elements(emission_line_parameters); Number of emission line parameter files
	n_abs_par = n_elements(absorption_line_parameters);Number of absorption line parameter files

	success = MDAP_CHECK_FILE_EXISTS(template_libraries, /search)	; Aborts of fails
	if success ne 0 then $
	    message, 'Template libraries not correctly defined.'
		
	success = MDAP_CHECK_FILE_EXISTS(emission_line_parameters)
	if success ne 0 then $
	    message, 'Emission-line parameter files not correctly defined.'

	success = MDAP_CHECK_FILE_EXISTS(absorption_line_parameters)
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
	;-----------------------------------------------------------------------


	;-----------------------------------------------------------------------
	; TODO: Version control will eventually happen in
	; MDAP_GENERATE_OUTPUT_FILE_NAMES, or preferrably in a more well-named
	; procedure
	
;	manga_dap_version_previous=0

	mdap_read_drp_fits_version='0'
	mdap_calculate_sn_version='0'
	mdap_spatial_binning_version='0'

	MDAP_READ_DRP_FITS, version= mdap_read_drp_fits_version
	mdap_fiducial_bin_xy_version = '0.1'
	mdap_tpl_lib_setup_version = '0.1'
	MDAP_CALCULATE_SN, version=mdap_calculate_sn_version
	MDAP_SPATIAL_BINNING, version=mdap_spatial_binning_version

;	mdap_spatial_binning_version_previous=0
;
;	mdap_spectral_fitting_version=0
;	mdap_spectral_fitting_version_previous=0
;
;	mdap_measure_indices_version=0
;	mdap_measure_indices_version_previous=0
;
;	mdap_spatial_radial_binning_version=0
;	mdap_spatial_radial_binning_version_previous=0
;
;	execute_all_modules = 1
	;-----------------------------------------------------------------------


	if ~keyword_set(quiet) then $
	    print, 'BLOCK 0 ... DONE.'
	; END BLOCK 0 ----------------------------------------------------------



	; BLOCK 1 --------------------------------------------------------------
	;	Read the DRP fits file
	;-----------------------------------------------------------------------

	if ~keyword_set(quiet) then $
	    print, 'BLOCK 1 ... '

;	if mdap_read_datacube_version gt mdap_read_datacube_version_previous $
;	   or execute_all_modules eq 1 then begin

	; TODO: Version required for this should depend on the version of DRP
	; output being read

	; TODO: Put in the conditions for when this block is executed:
	;	- If any of the plans have the overwrite flag flipped
	;	- If any of the plans do not have a previously created output
	;	  file with the binned spectra

	; TODO: These operations are independent of the plan -------------------
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

	;-----------------------------------------------------------------------
	if ~keyword_set(quiet) then $
	    print, 'BLOCK 1 ... DONE.'
	; END BLOCK 1 ----------------------------------------------------------

	; TODO: Once the DRP fits file has been read, the required input to
	;	appropriatly manipulate ALL the template libraries is available


	; BLOCK 2 --------------------------------------------------------------
	;	- Read the template library
	;	- Match the resolution of the templates to that of the galaxy
	;	  spectra (as best as possible).
	;	- Resample the templates to match the object sampling
	;-----------------------------------------------------------------------

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
	; END BLOCK 2 ----------------------------------------------------------



	;-----------------------------------------------------------------------
	; ANALYSIS PREP --------------------------------------------------------

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

	;-----------------------------------------------------------------------
	; Now can begin cycling through execution plans
	if keyword_set(dbg) then $
	    n_plans=1				; Only run the first plan
	for i=0,n_plans-1 do begin

	    if ~keyword_set(quiet) then $
		print, 'Beginning ExecutionPlans: ', i+1

	    ; TODO: If previous file exists, read or overwrite it?

	    ; BLOCK 3 ----------------------------------------------------------
	    ; S/N CALCULATION
	    ;-------------------------------------------------------------------
	    if ~keyword_set(quiet) then $
		print, 'PLAN '+MDAP_STC(i+1,/integer)+': BLOCK 3 ...'

	    ; Determine the 'good' spectra based on a set of criteria defined by
	    ; this procedure
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
			       signal=signal, noise=noise, quiet=quiet

	    if ~keyword_set(quiet) then $
		print, 'PLAN '+MDAP_STC(i+1,/integer)+': BLOCK 3 ... DONE.'
	    ; END BLOCK 3 ------------------------------------------------------
	    ;-------------------------------------------------------------------


	    ; BLOCK 4 ----------------------------------------------------------
	    ; Spatial binning
	    ;-------------------------------------------------------------------
	    if ~keyword_set(quiet) then $
		print, 'PLAN '+MDAP_STC(i+1,/integer)+': BLOCK 4 ...'

	    ; TODO: If previous file exists, read or overwrite it?

	    MDAP_SPATIAL_BINNING, flux, ivar, mask, signal, noise, gflag, bskyx, bskyy, spaxel_dx, $
				  spaxel_dy, execution_plan[i].bin_type, $
				  execution_plan[i].bin_par, execution_plan[i].threshold_ston_bin, $
				  execution_plan[i].bin_weight_by_sn2, bin_wgts, bin_indx, $
				  bin_flux, bin_ivar, bin_mask, xbin, ybin, bin_area, bin_ston, $
				  nbin, sn_calibration=sn_calibration, plot=plot

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
			       bin_n=nbin, /read_header, quiet=quiet

	    if ~keyword_set(quiet) then $
		print, 'PLAN '+MDAP_STC(i+1,/integer)+': BLOCK 4 ... DONE.'
	    ; END BLOCK 4 ------------------------------------------------------
	    ;-------------------------------------------------------------------



	    ; ##################################################################
	    ; BEGIN computation of "model-independent" data products ###########



	    ;-------------------------------------------------------------------
	    ; Check if any analysis was requested, if not continue with the loop
	    indx = where(execution_plan[i].analysis eq 1)
	    if indx[0] eq -1 then $		; No analysis to perform
	    	continue

	    ;-------------------------------------------------------------------
	    ; Read the pre-existing template library
	    ;	The template library is required for ANY analysis
	    ; TODO: Should not be able to reach here without first creating tpl_out_fits!
	    library_key = tpl_library_keys[execution_plan[i].tpl_lib]; 'M11-MARCS'
	    tpl_out_fits = MDAP_SET_TPL_LIB_OUTPUT_FILE(output_file_root, library_key)
	   
	    MDAP_READ_RESAMPLED_TEMPLATES, tpl_out_fits, tpl_wave, tpl_flux, tpl_ivar, tpl_mask, $
					   tpl_sres, tpl_soff

	    ;-------------------------------------------------------------------
	    ; Read the emission line file; this is done at every iteration
	    ; because it's fairly quick.
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

	    ; BLOCK 5 ----------------------------------------------------------
	    ; Run a full spectral fit including stellar population, gas and
	    ; stellar kinematics
	    ;-------------------------------------------------------------------

	    if ~keyword_set(quiet) then $
		print, 'PLAN '+MDAP_STC(i+1,/integer)+': BLOCK 5 ...'

	    ; Set the starting guesses for the kinematics
	    ; TODO: get the guesses from a cross-correlation with a guess template?

	    ; TODO: convert this to a procedure that could also read in guesses
	    ; from a previous DAP-produced file
	    star_kin_starting_guesses = dblarr(n_elements(xbin),4)
	    star_kin_starting_guesses[*,0] = velocity_initial_guess		; stellar velocity
	    star_kin_starting_guesses[*,1] = velocity_dispersion_initial_guess	; stellar sigma
	    gas_kin_starting_guesses = star_kin_starting_guesses[*,0:1]		; gas velocity
	    gas_kin_starting_guesses[*,1]=50.					; gas sigma

	    ; Set ppxf_only flag; if true (1), only PPXF is run, not GANDALF
	    ppxf_only=0
	    if execution_plan[i].analysis[0] eq 1 and execution_plan[i].analysis[1] eq 0 then $
		ppxf_only=1		; Only execute the PPXF step

	    ; Perform the spectral fit
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
				   ppxf_only=ppxf_only, quiet=quiet, plot=plot, dbg=dbg

	    ; TODO: Add the spectral fitting version to the header of the output file
	    ; TODO: Add information to the log file

	    ; Write the analysis wavelength range to the header
	    MDAP_WRITE_OUTPUT, execution_plan[i].ofile, header=header, $
			       w_range_analysis=execution_plan[i].wave_range_analysis, $
			       threshold_ston_analysis=execution_plan[i].threshold_ston_analysis, $
			       /read_header, quiet=quiet
	    
	    ; Write the emission line parameters
	    if (execution_plan[i].analysis[0] eq 1 and n_elements(eml_par) ne 0 ) or $
		execution_plan[i].analysis[1] eq 1 then begin
		MDAP_WRITE_OUTPUT, execution_plan[i].ofile, eml_par=eml_par, quiet=quiet
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
			       analysis_par=execution_plan[i].analysis_par, /read_header, $
			       quiet=quiet

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
				   bestfit_gndf=bestfit_gndf, eml_model=eml_model, quiet=quiet
	    endif

	    ; Write the optimal templates and optimal template convolved with
	    ; the LOSVD; see MDAP_SPECTRAL_FITTING on how these are generated.
	    MDAP_WRITE_OUTPUT, execution_plan[i].ofile, optimal_template=best_template, $
			       losvd_optimal_template=best_template_losvd_conv, quiet=quiet

	    if ~keyword_set(quiet) then $
		print, 'PLAN '+MDAP_STC(i+1,/integer)+': BLOCK 5 ... DONE.'

	    ; END BLOCK 5 ------------------------------------------------------
	    ;-------------------------------------------------------------------


	    ; END computation of "model-independent" data products #############
	    ; ##################################################################


	endfor
	; END loop over plans --------------------------------------------------
	;-----------------------------------------------------------------------

	; close up
	if ~keyword_set(nolog) then $
	    MDAP_LOG_CLOSE, log_file_unit

END


