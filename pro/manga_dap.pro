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
;	At some point I (KBW) switched from using 'absorption-line indices' to
;	'spectral indices', primarily because the spectral-index analysis also
;	measures the strengths of spectral breaks/bandheads like D4000.  This
;	just a choice of nomenclature.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED (OUT-OF-DATE):
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
; INTERNAL SUPPORT ROUTINES (OUT-OF-DATE):
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
;	11 Nov 2014: (KBW) Inclusion of spectral index measurements, addition of
;			   infrastructure for used to check for pre-existing
;			   analysis results.  The latter is a bit messy.  Will
;			   hopefully be cleaned up when we move to python.
;-
;-------------------------------------------------------------------------------


;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
; Functions/procedures used to for the log file
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

;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------



;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
; Sets the file names

; TODO: Eventually this will generate the files names AND check them against
; existing files.  Any execution plan that matches IN DETAIL the execution plan
; of an existing file (apart from the analyses performed) will be used to input
; already completed analyses, unless the overwrite flag is flipped.   This check
; will *include* checks of the module versions used to generate the existing
; analysis.
;
;	Some of this is already built.  See below.

; TODO: Create procedures that read and write execution plans to the headers of
; fits files?

PRO MDAP_GENERATE_OUTPUT_FILE_NAMES, $
		file_root, execution_plan

	n_plans = n_elements(execution_plan)
	for i=0,n_plans-1 do begin
	    file=file_root+'BIN-'+execution_plan[i].bin_type+'-'+string(i+1,format='(I03)')+'.fits'
	    execution_plan[i].ofile = file
	endfor
END

FUNCTION MDAP_SET_TPL_LIB_OUTPUT_FILE, $
		file_root, library_key, abs_line_key=abs_line_key
	if n_elements(abs_line_key) eq 0 then $
	    return, file_root+library_key+'.fits'

	return, file_root+library_key+'_'+abs_line_key+'.fits'
END
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------


;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
; TODO: Infrastructure used to check existing data.  This needs some work...
;
; The main function is MDAP_ANALYSIS_BLOCKS_TO_PERFORM.
;
; TODO: Describe these algorithms and the decision tree regarding which blocks
; to perform given the set of desired analyses provided by the user.

FUNCTION MDAP_OUTPUT_TABLE_ROWS, $
		file, exten
	unit = fxposit(file, exten)				; Set the header unit
	hdr = headfits(unit)					; Read it
	free_lun, unit						; Free and close up the LUN
	return, fxpar(hdr, 'NAXIS2')				; Read and return number of rows
END

FUNCTION MDAP_OUTPUT_COLUMN_SIZE, $
		file, exten, col

	fxbopen, unit, file, exten
	ndim_col = fxbdimen(unit, col)
	fxbclose, unit
	free_lun, unit
	return, ndim_col
END

FUNCTION MDAP_OUTPUT_IMAGE_SIZE, $
		file, exten
	unit = fxposit(file, exten)
	hdr=headfits(unit)
	free_lun, unit
	return, [ fxpar(hdr, 'NAXIS1'), fxpar(hdr, 'NAXIS2') ]
END

FUNCTION MDAP_OUTPUT_IMAGE_SIZE_1D, $
		file, exten
	unit = fxposit(file, exten)
	hdr=headfits(unit)
	free_lun, unit
	return, fxpar(hdr, 'NAXIS1')
END

; Tests if the existing STFIT data in the file matches the expectation from the
; input analysis request
FUNCTION MDAP_CAN_USE_STFIT_DATA, $
		file, tpl_fits, analysis_par

	bin_dim = MDAP_OUTPUT_IMAGE_SIZE(file, 'FLUX')	; Dimensions of the binned spectra
	ntpl = (MDAP_OUTPUT_IMAGE_SIZE(tpl_fits, 'FLUX'))[0]

	; The STFIT table must be populated with the correct number of rows
	if MDAP_OUTPUT_TABLE_ROWS(file, 'STFIT') ne bin_dim[0] then $
	    return, 0

	; The columns of the STFIT table must have the correct size
	if MDAP_OUTPUT_COLUMN_SIZE(file, 'STFIT', 'TPLW') ne ntpl then $
	    return, 0
	if analysis_par.degree gt 0 then begin
	    if MDAP_OUTPUT_COLUMN_SIZE(file, 'STFIT', 'ADDPOLY') ne analysis_par.degree then $
		return, 0
	endif
	if analysis_par.mdegree gt 0 then begin
	    if MDAP_OUTPUT_COLUMN_SIZE(file, 'STFIT', 'MULTPOLY') ne analysis_par.mdegree then $
		return, 0
	endif
	if MDAP_OUTPUT_COLUMN_SIZE(file, 'STFIT', 'KIN') ne analysis_par.moments then $
	    return, 0
	if MDAP_OUTPUT_COLUMN_SIZE(file, 'STFIT', 'KINERR') ne analysis_par.moments then $
	    return, 0

	; The SMSK and SMOD images must be the same size as the binned spectra
	msk_dim = MDAP_OUTPUT_IMAGE_SIZE(file, 'SMSK')
	if msk_dim[0] ne bin_dim[0] or msk_dim[1] ne bin_dim[1] then $
	    return, 0
	mod_dim = MDAP_OUTPUT_IMAGE_SIZE(file, 'SMOD')
	if mod_dim[0] ne bin_dim[0] or mod_dim[1] ne bin_dim[1] then $
	    return, 0

	return, 1
END

PRO MDAP_ADD_STFIT, $
		file, tpl_fits, analysis_par, perform_block

	if MDAP_CAN_USE_STFIT_DATA(file, tpl_fits, analysis_par) eq 1 then begin
	    perform_block.spec_fit = 0	; Do not perform the spectral fitting
	endif else $
	    perform_block.spec_fit = 1	; Peform the spectral fitting block

	perform_block.ppxf_only = 1	; Only perform PPXF
END


; Tests if the existing SGFIT data in the file can be used directly
; TODO: NUMBER OF MOMENTS IS HARDWIRED!!
FUNCTION MDAP_CAN_USE_SGFIT_DATA, $
		file, tpl_fits, eml_par, analysis_par

	bin_dim = MDAP_OUTPUT_IMAGE_SIZE(file, 'FLUX')	; Dimensions of the binned spectra
	ntpl = (MDAP_OUTPUT_IMAGE_SIZE(tpl_fits, 'FLUX'))[0]	; Number of templates
	neml = n_elements(eml_par)				; Number of emission lines

	; The SGFIT table must have the correct number of rows (one per binned spectrum)
	if MDAP_OUTPUT_TABLE_ROWS(file, 'SGFIT') ne bin_dim[0] then $
	    return, 0

	; The columns of the SGFIT table must have the correct size
	if MDAP_OUTPUT_COLUMN_SIZE(file, 'SGFIT', 'TPLW') ne ntpl then $
	    return, 0
	if analysis_par.mdegree gt 0 then begin
	    if MDAP_OUTPUT_COLUMN_SIZE(file, 'SGFIT', 'MULTPOLY') ne analysis_par.mdegree then $
		return, 0
	endif
	if MDAP_OUTPUT_COLUMN_SIZE(file, 'SGFIT', 'KIN') ne 2 then $
	    return, 0
	if MDAP_OUTPUT_COLUMN_SIZE(file, 'SGFIT', 'KINERR') ne 2 then $
	    return, 0
	if analysis_par.reddening_order gt 0 then begin
	    if MDAP_OUTPUT_COLUMN_SIZE(file,'SGFIT','RED') ne analysis_par.reddening_order then $
		return, 0
	    if MDAP_OUTPUT_COLUMN_SIZE(file,'SGFIT','REDERR') ne analysis_par.reddening_order then $
		return, 0
	endif
	if MDAP_OUTPUT_COLUMN_SIZE(file, 'SGFIT', 'ELOMIT') ne neml then $
	    return, 0
	if MDAP_OUTPUT_COLUMN_SIZE(file, 'SGFIT', 'AMPL') ne neml then $
	    return, 0
	if MDAP_OUTPUT_COLUMN_SIZE(file, 'SGFIT', 'AMPLERR') ne neml then $
	    return, 0
	if MDAP_OUTPUT_COLUMN_SIZE(file, 'SGFIT', 'FLUX') ne neml then $
	    return, 0
	if MDAP_OUTPUT_COLUMN_SIZE(file, 'SGFIT', 'FLUXERR') ne neml then $
	    return, 0
	if MDAP_OUTPUT_COLUMN_SIZE(file, 'SGFIT', 'EW') ne neml then $
	    return, 0
	if MDAP_OUTPUT_COLUMN_SIZE(file, 'SGFIT', 'EWERR') ne neml then $
	    return, 0
	dim = MDAP_OUTPUT_COLUMN_SIZE(file, 'SGFIT', 'IKIN')
	if dim[0] ne neml or dim[1] ne 2 then $
	    return, 0
	dim = MDAP_OUTPUT_COLUMN_SIZE(file, 'SGFIT', 'IKINERR')
	if dim[0] ne neml or dim[1] ne 2 then $
	    return, 0

	; The SGMSK, SGMOD, and ELMOD images must be the same size as the binned spectra
	msk_dim = MDAP_OUTPUT_IMAGE_SIZE(file, 'SGMSK')
	if msk_dim[0] ne bin_dim[0] or msk_dim[1] ne bin_dim[1] then $
	    return, 0
	mod_dim = MDAP_OUTPUT_IMAGE_SIZE(file, 'SGMOD')
	if mod_dim[0] ne bin_dim[0] or mod_dim[1] ne bin_dim[1] then $
	    return, 0
	elm_dim = MDAP_OUTPUT_IMAGE_SIZE(file, 'ELMOD')
	if elm_dim[0] ne bin_dim[0] or elm_dim[1] ne bin_dim[1] then $
	    return, 0

	return, 1
END

; TODO: Order is important with respect to MDAP_ADD_STFIT.  If spec_fit is
; already turned on, this function will not turn it off even if the SGFIT data is
; valid!  This is NOT true of MDAP_ADD_STFIT.
PRO MDAP_ADD_SGFIT, $
		file, tpl_fits, eml_par, analysis_par, perform_block
	; If it's already on, don't turn it off!
	if perform_block.spec_fit eq 0 then begin
	    if MDAP_CAN_USE_SGFIT_DATA(file, tpl_fits, eml_par, analysis_par) eq 1 then begin
		perform_block.spec_fit = 0	; Do not perform the spectral fitting
	    endif else $
		perform_block.spec_fit = 1	; Must perform the spectral fitting
	endif
	perform_block.ppxf_only = 0		; Perform both PPXF and GANDALF
END


; TODO: Add a check for the size of the images (SIWAVE, SIFLUX, etc)?
FUNCTION MDAP_CAN_USE_SINDX_DATA, $
		file, abs_par

	bin_dim = MDAP_OUTPUT_IMAGE_SIZE(file, 'FLUX')	; Dimensions of the binned spectra
	nabs = n_elements(abs_par)			; Number of spectral indices

	; The SINDX table must have the correct number of rows (one per binned spectrum)
	if MDAP_OUTPUT_TABLE_ROWS(file, 'SINDX') ne bin_dim[0] then $
	    return, 0

	; The columns of the SINDX table must have the correct size
	if MDAP_OUTPUT_COLUMN_SIZE(file, 'SINDX', 'SIOMIT') ne nabs then $
	    return, 0
	if MDAP_OUTPUT_COLUMN_SIZE(file, 'SINDX', 'INDX') ne nabs then $
	    return, 0
	if MDAP_OUTPUT_COLUMN_SIZE(file, 'SINDX', 'INDXERR') ne nabs then $
	    return, 0
	if MDAP_OUTPUT_COLUMN_SIZE(file, 'SINDX', 'INDX_OTPL') ne nabs then $
	    return, 0
	if MDAP_OUTPUT_COLUMN_SIZE(file, 'SINDX', 'INDX_BOTPL') ne nabs then $
	    return, 0

	return, 1
END

FUNCTION MDAP_CAN_USE_SINDX_IMAGES, $
		file

	wave_dim = MDAP_OUTPUT_IMAGE_SIZE_1D(file, 'SIWAVE')

	if wave_dim eq 1 then $				; No wave image
	    return, 0

	si_dim = MDAP_OUTPUT_IMAGE_SIZE(file, 'FLUX')	; Dimensions of the binned spectra
	si_dim[1] = wave_dim				; Force spectrum length

	; The SIFLUX, SIIVAR, SIMASK, SIOTPL, and SIBOTPL images must have the
	; same and correct dimensions.
	flx_dim = MDAP_OUTPUT_IMAGE_SIZE(file, 'SIFLUX')
	if flx_dim[0] ne si_dim[0] or flx_dim[1] ne si_dim[1] then $
	    return, 0
	ivr_dim = MDAP_OUTPUT_IMAGE_SIZE(file, 'SIIVAR')
	if ivr_dim[0] ne si_dim[0] or ivr_dim[1] ne si_dim[1] then $
	    return, 0
	msk_dim = MDAP_OUTPUT_IMAGE_SIZE(file, 'SIMASK')
	if msk_dim[0] ne si_dim[0] or msk_dim[1] ne si_dim[1] then $
	    return, 0
	otp_dim = MDAP_OUTPUT_IMAGE_SIZE(file, 'SIOTPL')
	if otp_dim[0] ne si_dim[0] or otp_dim[1] ne si_dim[1] then $
	    return, 0
	btp_dim = MDAP_OUTPUT_IMAGE_SIZE(file, 'SIBOTPL')
	if btp_dim[0] ne si_dim[0] or btp_dim[1] ne si_dim[1] then $
	    return, 0

	return, 1
END	


; If adding the spectral index measurements, the file must have:
;	- the template weights (need the same # of templates based on the provided template file)
;	- the stellar kinematics (need the same number of binned spectra as in the file)
;		THIS TEST CURRENTLY ONLY RUNS IF THE BINNED SPECTRA IN THE FILE ARE TO BE USED!
;	- the best fitting spectra (must be same size as binned spectra)
PRO MDAP_ADD_SINDX, $
		file, abs_par, perform_block

	if MDAP_CAN_USE_SINDX_DATA(file, abs_par) eq 1 then begin
	    perform_block.spec_indx = 0		; Do not perform the spectral-index measurements
	endif else $
	    perform_block.spec_indx = 1
END


; Set which analyses to perform for a new file
PRO MDAP_NEW_FILE_BLOCKS, analysis, perform_block

	; TODO: Somehow put this in a common area with MDAP_BUILD_EXECUTION_PLANS
	perform_block.ston = 1		; Always perform the S/N block
	perform_block.bin = 1		; Always perform the spatial binning

	; Order is important for the first two if statements because of ppxf_only!
	if analysis[0] eq 1 then begin	; Wants the stellar continuum fitting
	    perform_block.spec_fit = 1		; Perform the spectral fit
	    perform_block.ppxf_only = 1		; Only ppxf
	endif

	; For a new file, the emission-line fits, REQUIRES the stellar continuum fit
	if analysis[1] eq 1 then begin	; Wants the emission-line fitting
	    perform_block.spec_fit = 1		; Perform the spectral fit
	    perform_block.ppxf_only = 0		; Not only ppxf
	endif

	; For a new file, the spectral index fits, REQUIRES the spectral fitting
	if analysis[2] eq 1 then begin	; Wants the spectral indices
	    perform_block.spec_fit = 1				; Perform the spectral fit
	    perform_block.ppxf_only = analysis[1] eq 0 ? 1 : 0	; Check if ppxf only
	    perform_block.spec_indx = 1				; Do the spectral index measurements
	endif

	; If it reaches here, only the S/N calculation and the binning are performed
END

; TODO: Should only reach here if the output file already exists!
FUNCTION MDAP_PERFORM_STON_BLOCK, $
		execution_plan, ndrp

	; Conditions to return true:
	;	A: overwrite is true
	;	B: output DRPS extension does not have the correct size

	; TODO: Should also eventually check (?):
	;	- the spaxel size
	;	- the S/N wavelength range
	;	- the on-sky positions

	; CASE A ---------------------------------------------------------------
	if execution_plan.overwrite eq 1 then $
	    return, 1

	; CASE B ---------------------------------------------------------------
	;   MDAP_CHECK_OUTPUT_FILE has already:
	;	- Checked for the existence of the DRPS extension
	;	- Checked that all the appropriate columns exist
	;   Need to:
	;	- Check the number of rows against the expected value
	if MDAP_OUTPUT_TABLE_ROWS(execution_plan.ofile, 'DRPS') ne ndrp then $
	    return, 1

	return, 0		; Do not need to redo the S/N calculation
END

; TODO: Should only reach here if the output file already exists!
FUNCTION MDAP_PERFORM_BIN_BLOCK, $
		execution_plan, ndrp

	; Conditions to return true:
	;	A: overwrite is true
	;	B: output DRPS extension does not have the correct size
	;	C: output BINS extension does not have the correct size

	; CASE A ---------------------------------------------------------------
	if execution_plan.overwrite eq 1 then $
	    return, 1

	;   MDAP_CHECK_OUTPUT_FILE has already:
	;	- Checked for the existence of the DRPS and BINS extension
	;	- Checked that all the appropriate columns exist
	;   Need to:

	; CASE B ---------------------------------------------------------------
	;	- Check the number of DRPS rows against the expected value
	if MDAP_OUTPUT_TABLE_ROWS(execution_plan.ofile, 'DRPS') ne ndrp then $
	    return, 1

	; CASE C ---------------------------------------------------------------
	;	- Check the number of BINS rows against the expected value
	nbin = (MDAP_OUTPUT_IMAGE_SIZE(execution_plan.ofile, 'FLUX'))[0]
	if MDAP_OUTPUT_TABLE_ROWS(execution_plan.ofile, 'BINS') ne nbin then $
	    return, 1

	return, 0		; Do not need to redo the spatial binning
END

; TODO: I want to be able to add extensions to a file for new blocks, without
; having to redo old ones.  Will it do this already?
FUNCTION MDAP_ANALYSIS_BLOCKS_TO_PERFORM, $
		execution_plan, ndrp, tpl_fits, eml_par, abs_par

	perform_block = { RequiredAnalysisBlock, ston:0, bin:0, spec_fit:0, ppxf_only:0, $
						 spec_indx:0 }

	; Determine if the output file already exists
	file_exists = file_test(execution_plan.ofile)

	; - If the file exists but does not look like a DAP file and the
	;   overwrite flag has been set to false, throw an error
	if file_exists eq 1 and execution_plan.overwrite eq 0 then begin
	    if MDAP_CHECK_OUTPUT_FILE(execution_plan.ofile) eq 0 then $
		message, execution_plan.ofile+' does not look like a DAP file.  To continue, ' + $
			 'set overwrite flag to true (1).'
	endif

	; - If the file does not exist, make sure all (requested) blocks are
	;   performed sequentially
	if file_exists eq 0 then begin
	    MDAP_NEW_FILE_BLOCKS, execution_plan.analysis, perform_block
	    return, perform_block
	endif

	; TODO: If overwriting, some of the existing data may not correspond to
	; the new analysis, but may not be overwritten.  How do I deal with this?
	; Via the DATEMOD keyword?

	; - The file exists, so now check the data in the different extentions

	; TODO: This test is not really necessary; everything hinges on bin
	; Peform the S/N block if:
	;	A: overwrite is true
	;	B: output DRPS extension does not have the correct size
	;	C: not determined by MDAP_PERFORM_STON_BLOCK(), but below
	perform_block.ston = MDAP_PERFORM_STON_BLOCK(execution_plan, ndrp)
	
	; Perform the BIN block if:
	;	A; overwrite is true
	;	B: output DRPS extension does not have the correct size
	;	C: output BINS extension does not have the correct size
	perform_block.bin = MDAP_PERFORM_BIN_BLOCK(execution_plan, ndrp)

	; This is redundant with the next if statement
;	if perform_block.bin eq 1 and perform_block.ston eq 0 then $
;	    perform_block.ston = 1			; Always calculate S/N if binning

	if perform_block.bin eq 1 then begin		; If rebinning, treat like a new file
	    MDAP_NEW_FILE_BLOCKS, execution_plan.analysis, perform_block
	    return, perform_block
	endif

	print, 'deciding'

	; Order is important for these statements!
	if execution_plan.analysis[0] eq 1 then begin
	    MDAP_ADD_STFIT, execution_plan.ofile, tpl_fits, execution_plan.analysis_par, $
			    perform_block
	endif
	if execution_plan.analysis[1] eq 1 then begin
	    MDAP_ADD_SGFIT, execution_plan.ofile, tpl_fits, eml_par, execution_plan.analysis_par, $
			    perform_block
	endif

	if execution_plan.analysis[2] eq 1 then begin	; Wants the spectral indices
	    ; First check that at least the STFIT results are available
	    if execution_plan.analysis[0] eq 0 and execution_plan.analysis[1] eq 0 then begin
		MDAP_ADD_STFIT, execution_plan.ofile, tpl_fits, execution_plan.analysis_par, $
				perform_block
	    endif
	    MDAP_ADD_SINDX, execution_plan.ofile, abs_par, perform_block
	endif

	return, perform_block				; All blocks should have been completed
END

FUNCTION MDAP_ALL_ANALYSIS_BLOCKS_COMPLETED, $
		perform_block

	if perform_block.ston eq 0 and $
	   perform_block.bin eq 0 and $
	   perform_block.spec_fit eq 0 and $
	   perform_block.spec_indx eq 0 then begin
	    return, 1					; All the blocks have been completed
	endif

	return, 0					; Blocks need to be done
END

FUNCTION MDAP_MORE_ANALYSIS_BLOCKS_TO_FINISH, $
		perform_block, ston=ston, bin=bin, spec_fit=spec_fit, spec_indx=spec_indx
	if keyword_set(ston) then begin
	    if perform_block.bin eq 0 and $
	       perform_block.spec_fit eq 0 and $
	       perform_block.spec_indx eq 0 then begin
		return, 0
	    endif
	endif

	if keyword_set(bin) then begin
	    if perform_block.spec_fit eq 0 and $
	       perform_block.spec_indx eq 0 then begin
		return, 0
	    endif
	endif

	if keyword_set(spec_fit) then begin
	    if perform_block.spec_indx eq 0 then begin
		return, 0
	    endif
	endif

	if keyword_set(spec_indx) then $
	    return, 0

	return, 1					; Blocks need to be done
END
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------



;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
; This is the main wrapper procedure for the DAP.
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
			      template_libraries, tpl_vacuum_wave, ems_line_keys, $
			      emission_line_parameters, ems_vacuum_wave, abs_line_keys, $
			      absorption_line_parameters, abs_vacuum_wave, signifier, bin_type, $
			      bin_par, w_range_sn, threshold_ston_bin, bin_weight_by_sn2, $
			      w_range_analysis, threshold_ston_analysis, analysis, $
			      tpl_lib_analysis, ems_par_analysis, abs_par_analysis, $
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

	; Check that the external_library exists!
	if n_elements(external_library) then begin
	    bvls_shared_lib=external_library+'bvls.so'
	    if file_test(bvls_shared_lib) eq 0 then $
		bvls_shared_lib=external_library+'bvls.dylib'
	    if file_test(bvls_shared_lib) eq 0 then $
		message, 'Shared object library does not exist!'
	endif

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
	
	mdap_read_drp_fits_version='0'
	mdap_fiducial_bin_xy_version = '0'
	mdap_tpl_lib_setup_version='0'
	mdap_calculate_sn_version='0'
	mdap_spatial_binning_version='0'
	mdap_spectral_fitting_version='0'
	mdap_spectral_index_version='0'

	MDAP_READ_DRP_FITS, version=mdap_read_drp_fits_version
	MDAP_FIDUCIAL_BIN_XY, version=mdap_fiducial_bin_xy_version
	MDAP_CREATE_DAP_TEMPLATE_LIBRARY, version=mdap_tpl_lib_setup_version
	MDAP_CALCULATE_SN, version=mdap_calculate_sn_version
	MDAP_SPATIAL_BINNING, version=mdap_spatial_binning_version
	MDAP_SPECTRAL_FITTING, version=mdap_spectral_fitting_version
	MDAP_SPECTRAL_INDEX_MEASUREMENTS, version=mdap_spectral_index_version

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
	    ndrp = (size(flux))[1]		; Number of DRP spectra

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

	; Get velocity scale of the galaxy data
	velScale=MDAP_VELOCITY_SCALE(wave, /log10)

	; Cycle through all the libraries, generating the template libraries to
	; use for the analysis, if they don't already exist!
	indx=where(use_tpl_lib eq 1)
	if indx[0] ne -1 then begin

	    for i=0,n_tpl_lib-1 do begin
		if use_tpl_lib[i] ne 1 then $
		    continue

		; TODO: Need to rethink usage of setup procedure to define library keyword?

		; TODO: What are the flux units for the templates?  Does it
		; matter because I'm renormalizing them?

		; TODO: Should the convolution conserve the integral of the template spectra?

		; Working template library used with MDAP_SPECTRAL_FITTING
		tpl_out_fits = MDAP_SET_TPL_LIB_OUTPUT_FILE(output_file_root, tpl_library_keys[i])

		; The file does not exist, so create it
		if file_test(tpl_out_fits) eq 0 then begin
		    MDAP_CREATE_DAP_TEMPLATE_LIBRARY, tpl_out_fits, tpl_library_keys[i], $
						      template_libraries[i], tpl_vacuum_wave[i], $
						      velocity_initial_guess, wave, sres, $
						      velscale, quiet=quiet

		    if ~keyword_set(nolog) then begin
			printf, log_file_unit, '[INFO] Read template library: '+tpl_library_keys[i]
			printf, log_file_unit, '[INFO] Resolution and sampling matched to: ' + $
				datacube_name
			printf, log_file_unit,'[INFO] Template library setup done using version: ' $
				+ mdap_tpl_lib_setup_version
			printf, log_file_unit, '[INFO] Results written to: ' + tpl_out_fits
		    endif

		endif

		for j=0,n_plans-1 do begin

		    ; Spectral-index analysis not performed so continue
		    if execution_plan[j].abs_par eq -1 then $
			continue

		    ; Working template library used with MDAP_MEASURE_INDICES
		    tpl_out_fits = MDAP_SET_TPL_LIB_OUTPUT_FILE(output_file_root, $
								tpl_library_keys[i], $
					abs_line_key=abs_line_keys[execution_plan[j].abs_par])

		    ; The file exists, so continue
		    if file_test(tpl_out_fits) eq 1 then $
			continue

		    ; Get the spectral resolution vector
		    abs_sres = MDAP_ABS_LINE_RESOLUTION(abs_line_keys[execution_plan[j].abs_par], $
							wave)

		    ; Create the template library with resolution matched to the
		    ; spectral-index system
		    MDAP_CREATE_DAP_TEMPLATE_LIBRARY, tpl_out_fits, tpl_library_keys[i], $
						      template_libraries[i], tpl_vacuum_wave[i], $
						      velocity_initial_guess, wave, abs_sres, $
						      velscale, quiet=quiet

		    if ~keyword_set(nolog) then begin
			printf, log_file_unit, '[INFO] Read template library: '+tpl_library_keys[i]
			printf, log_file_unit, '[INFO] Spectral index library: ' + $
					       abs_line_keys[execution_plan[j].abs_par]
			printf, log_file_unit,'[INFO] Template library setup done using version: ' $
				+ mdap_tpl_lib_setup_version
			printf, log_file_unit, '[INFO] Results written to: ' + tpl_out_fits
		    endif

		endfor	; Loop over execution plans
	    endfor	; Loop over template libraries
	endif		; If any template libraries are used

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

	    ;-------------------------------------------------------------------
	    ; Get the nominal name for the template library.  This is used to
	    ; determine which blocks to execute (via perform_block).  This file
	    ; is only used if checking to add the SINDX results.  I.e., this
	    ; needs to be here to pass to MDAP_ANALYSIS_BLOCKS_TO_PERFORM, but
	    ; it is not used unless execution_plan[i].analysis[2] eq 1.
	    if execution_plan[i].tpl_lib ne -1 then begin
		tpl_out_fits = MDAP_SET_TPL_LIB_OUTPUT_FILE(output_file_root, $
							tpl_library_keys[execution_plan[i].tpl_lib])
	    endif

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

	    ;-------------------------------------------------------------------
	    ; Read the spectral-index parameters; this is done at every
	    ; iteration because it's fairly quick.  TODO: Check that this I/O
	    ; doesn't slow things down too much!
	    if execution_plan[i].abs_par ne -1 then begin
		abs_par = MDAP_READ_ABSORPTION_LINE_PARAMETERS( absorption_line_parameters[ $
								 execution_plan[i].abs_par ])
		nabs = n_elements(abs_par)		; Number of absorption-line indices
		if abs_vacuum_wave[ execution_plan[i].abs_par ] eq 0 then begin
		    for j=0,nabs-1 do begin
			AIRTOVAC, abs_par[j].passband, lam
			abs_par[j].passband = lam
			AIRTOVAC, abs_par[j].blue_cont, lam
			abs_par[j].blue_cont = lam
			AIRTOVAC, abs_par[j].red_cont, lam
			abs_par[j].red_cont = lam
		    endfor
		endif
	    endif else $
		MDAP_ERASEVAR, abs_par	; Make sure it doesn't exist if no abs_par is defined

	    ; Setup which analysis blocks (3-6) to perform
	    perform_block = MDAP_ANALYSIS_BLOCKS_TO_PERFORM(execution_plan[i], ndrp, tpl_out_fits, $
							    eml_par, abs_par)

	    if MDAP_ALL_ANALYSIS_BLOCKS_COMPLETED(perform_block) eq 1 then begin
		if ~keyword_set(quiet) then $
		    print, 'All blocks previously completed for ExecutionPlan: ', i+1
		continue
	    endif

	    ; Alert the user which blocks will be performed
	    print, 'Blocks to perform:'
	    if perform_block.ston eq 1 then begin
		print, '    S/N calculation: YES'
	    endif else $
		print, '    S/N calculation: NO'
	    if perform_block.bin eq 1 then begin
		print, '    Spatial binning: YES'
	    endif else $
		print, '    Spatial binning: NO'
	    if perform_block.spec_fit eq 1 then begin
		print, '    Spectral fitting: YES'
		if perform_block.ppxf_only eq 1 then begin
		    print, '        PPXF ONLY: YES'
		endif else $
		    print, '        PPXF ONLY: NO'
	    endif else $
		print, '    Spectral fitting: NO'
	    if perform_block.spec_indx eq 1 then begin
		print, '    Spectral indices: YES'
	    endif else $
		print, '    Spectral indices: NO'

	    ; BLOCK 3 ----------------------------------------------------------
	    ; S/N CALCULATION
	    ;-------------------------------------------------------------------
	    if ~keyword_set(quiet) then $
		print, 'PLAN '+MDAP_STC(i+1,/integer)+': BLOCK 3 ...'

	    if perform_block.ston eq 1 then begin

		; Determine the 'good' spectra based on a set of criteria
		; defined by this procedure

		; TODO: Currently gindx is not used.  gflag is used by
		;	MDAP_CALCULATE_SN and MDAP_SPATIAL_BINNING

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
		MDAP_WRITE_OUTPUT, execution_plan[i].ofile, header=header, dx=spaxel_dx, $
				   dy=spaxel_dy, w_range_sn=execution_plan[i].wave_range_sn, $
				   xpos=bskyx, ypos=bskyy, signal=signal, noise=noise, quiet=quiet

	    endif else begin
		; TODO: spaxel_dx, spaxel_dy, bskyx, and bskyy are always
		;	calculated!  Should they be reread here?

		print, 'READING EXISTING S/N DATA'

		MDAP_DEFINE_OUTPUT, header=header, dx=spaxel_dx, dy=spaxel_dy, xpos=bskyx, $
				    ypos=bskyy, signal=signal, noise=noise
		MDAP_READ_OUTPUT, execution_plan[i].ofile, header=header, dx=spaxel_dx, $
				   dy=spaxel_dy, xpos=bskyx, ypos=bskyy, signal=signal, noise=noise
	    endelse

	    if ~keyword_set(quiet) then $
		print, 'PLAN '+MDAP_STC(i+1,/integer)+': BLOCK 3 ... DONE.'
	    ; END BLOCK 3 ------------------------------------------------------
	    ;-------------------------------------------------------------------

	    if MDAP_MORE_ANALYSIS_BLOCKS_TO_FINISH(perform_block, /ston) eq 0 then $
		continue

	    ; BLOCK 4 ----------------------------------------------------------
	    ; Spatial binning
	    ;-------------------------------------------------------------------
	    if ~keyword_set(quiet) then $
		print, 'PLAN '+MDAP_STC(i+1,/integer)+': BLOCK 4 ...'

	    if perform_block.bin eq 1 then begin
	    
		MDAP_SPATIAL_BINNING, flux, ivar, mask, signal, noise, gflag, bskyx, bskyy, $
				      spaxel_dx, spaxel_dy, execution_plan[i].bin_type, $
				      execution_plan[i].bin_par, $
				      execution_plan[i].threshold_ston_bin, $
				      execution_plan[i].bin_weight_by_sn2, bin_wgts, bin_indx, $
				      bin_flux, bin_ivar, bin_mask, xbin, ybin, bin_area, $
				      bin_ston, nbin, sn_calibration=sn_calibration, plot=plot

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
			if n_elements(execution_plan[i].bin_par) then begin
			    printf, log_file_unit, '[INFO] Bin parameter: ', $
				    execution_plan[i].bin_par
			endif
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
	    endif else begin

		print, 'READING EXISTING BIN DATA'

		; Read the binning results
		MDAP_DEFINE_OUTPUT, header=header, bin_indx=bin_indx, bin_weights=bin_wgts, $
				    bin_flux=bin_flux, bin_ivar=bin_ivar, bin_mask=bin_mask, $
				    xbin=xbin, ybin=ybin, bin_area=bin_area, bin_ston=bin_ston, $
				    bin_n=nbin
		MDAP_READ_OUTPUT, execution_plan[i].ofile, header=header, bin_indx=bin_indx, $
				  bin_weights=bin_wgts, bin_flux=bin_flux, bin_ivar=bin_ivar, $
				  bin_mask=bin_mask, xbin=xbin, ybin=ybin, bin_area=bin_area, $
				  bin_ston=bin_ston, bin_n=nbin

		; Check the size against the read in DRP data
		; TODO: Given the way the perform_block structure is created, is this necessary?
		sz = size(bin_flux)
		if sz[2] ne n_elements(wave) then $
		    message, 'Read binned spectra do not have the same length as the DRP spectra!'

	    endelse

	    if ~keyword_set(quiet) then $
		print, 'PLAN '+MDAP_STC(i+1,/integer)+': BLOCK 4 ... DONE.'
	    ; END BLOCK 4 ------------------------------------------------------
	    ;-------------------------------------------------------------------


	    if MDAP_MORE_ANALYSIS_BLOCKS_TO_FINISH(perform_block, /bin) eq 0 then $
		continue


	    ; ##################################################################
	    ; BEGIN computation of "model-independent" data products ###########


	    ; PREP FOR BLOCK 5 -------------------------------------------------
	    ;-------------------------------------------------------------------
	    ; Check if any analysis was requested, if not continue with the loop
	    indx = where(execution_plan[i].analysis eq 1)
	    if indx[0] eq -1 then $		; No analysis to perform
	    	continue

	    ;-------------------------------------------------------------------
	    ; Read the pre-existing template library
	    ;	The template library is required for ANY analysis
	    ; TODO: Should not be able to reach here without first creating tpl_out_fits!
	    tpl_out_fits = MDAP_SET_TPL_LIB_OUTPUT_FILE(output_file_root, $
							tpl_library_keys[execution_plan[i].tpl_lib])
	   
	    MDAP_READ_RESAMPLED_TEMPLATES, tpl_out_fits, tpl_wave, tpl_flux, tpl_ivar, tpl_mask, $
					   tpl_sres, tpl_soff

	    ; TODO: Need to account for tpl_soff further down the line

	    ; BLOCK 5 ----------------------------------------------------------
	    ; Run a full spectral fit including stellar population, gas and
	    ; stellar kinematics
	    ;-------------------------------------------------------------------

	    if ~keyword_set(quiet) then $
		print, 'PLAN '+MDAP_STC(i+1,/integer)+': BLOCK 5 ...'

	    if perform_block.spec_fit eq 1 then begin

		; TODO: Allow for starting guesses from a cross-correlation with
		; a guess template?

		; TODO: Allow for starting guesses from a previous DAP-produced
		; file

		; Set the starting guesses for the kinematics
		star_kin_guesses = dblarr(n_elements(xbin),4)
		star_kin_guesses[*,0] = velocity_initial_guess			; stellar velocity
		star_kin_guesses[*,1] = velocity_dispersion_initial_guess	; stellar sigma
		gas_kin_guesses = star_kin_guesses[*,0:1]			; gas velocity
		gas_kin_guesses[*,1]=50.					; gas sigma

		; Perform the spectral fit
		MDAP_SPECTRAL_FITTING, wave, bin_flux, bin_ivar, bin_mask, sres, tpl_wave, $
				       tpl_flux, tpl_ivar, tpl_mask, wavelength_output, $
				       obj_fit_mask_ppxf, weights_ppxf, add_poly_coeff_ppxf, $
				       mult_poly_coeff_ppxf, bestfit_ppxf, chi2_ppxf, $
				       obj_fit_mask_gndf, weights_gndf, mult_poly_coeff_gndf, $
				       bestfit_gndf, chi2_gndf, eml_model, stellar_kinematics, $
				       stellar_kinematics_err, emission_line_kinematics, $
				       emission_line_kinematics_err, emission_line_omitted, $
				       emission_line_kinematics_individual, $
				       emission_line_kinematics_individual_err, $
				       emission_line_intens, emission_line_intens_err, $
				       emission_line_fluxes, emission_line_fluxes_err, $
				       emission_line_EW, emission_line_EW_err, reddening_output, $
				       reddening_output_err, $
				       analysis_par=execution_plan[i].analysis_par, $
				       star_kin_starting_guesses=star_kin_guesses, $
				       gas_kin_starting_guesses=gas_kin_guesses, eml_par=eml_par, $
				       external_library=bvls_shared_lib, $
				       wave_range_analysis=execution_plan[i].wave_range_analysis, $
				       ppxf_only=perform_block.ppxf_only, quiet=quiet, plot=plot, $
				       dbg=dbg

		; TODO: Add the spectral fitting version to the header of the
		; output file, and add information to the log file

		; Write the analysis wavelength range to the header
		MDAP_WRITE_OUTPUT, execution_plan[i].ofile, header=header, $
				   w_range_analysis=execution_plan[i].wave_range_analysis, $
				threshold_ston_analysis=execution_plan[i].threshold_ston_analysis, $
				   /read_header, quiet=quiet
	    
		; Write the emission line parameters
		if n_elements(eml_par) ne 0 then $
		    MDAP_WRITE_OUTPUT, execution_plan[i].ofile, eml_par=eml_par, quiet=quiet

		; Write the stellar kinematics results; will always be done if
		; program makes it here
		MDAP_WRITE_OUTPUT, execution_plan[i].ofile, header=header, $
				   tpl_library_key=tpl_library_keys[i], $
				   obj_fit_mask_ppxf=obj_fit_mask_ppxf, weights_ppxf=weights_ppxf, $
				   add_poly_coeff_ppxf=add_poly_coeff_ppxf, $
				   mult_poly_coeff_ppxf=mult_poly_coeff_ppxf, $
				   bestfit_ppxf=bestfit_ppxf, chi2_ppxf=chi2_ppxf, $
				   stellar_kinematics_fit=stellar_kinematics, $
				   stellar_kinematics_err=stellar_kinematics_err, $
				   analysis_par=execution_plan[i].analysis_par, /read_header, $
				   quiet=quiet

		; Write the results of the emission-line fitting
		if perform_block.ppxf_only eq 0 then begin
		    MDAP_WRITE_OUTPUT, execution_plan[i].ofile, header=header, $
				       ems_line_key=ems_line_keys[i], eml_par=eml_par, $
				       obj_fit_mask_gndf=obj_fit_mask_gndf, $
				       weights_gndf=weights_gndf, $
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
				       bestfit_gndf=bestfit_gndf, eml_model=eml_model, $
				       quiet=quiet, /read_header
		endif

		; Write the optimal templates and optimal template convolved
		; with the LOSVD; see MDAP_SPECTRAL_FITTING on how these are
		; generated.
;		MDAP_WRITE_OUTPUT, execution_plan[i].ofile, optimal_template=best_template, $
;				   losvd_optimal_template=best_template_losvd_conv, quiet=quiet
	    endif else begin

		print, 'reading spec_fit data'

		; TODO: Here, I attempt to read the emission-line fitting
		; results.  The choice of the weights to use in the next block
		; critically depends on the size of the weights_gndf array.  If
		; it doesn't match the necessary size, weights_ppxf is used!
		; Add a flag to the header saying the extension is populated?
		MDAP_DEFINE_OUTPUT, obj_fit_mask_ppxf=obj_fit_mask_ppxf, $
				    weights_ppxf=weights_ppxf, $
				    add_poly_coeff_ppxf=add_poly_coeff_ppxf, $
				    mult_poly_coeff_ppxf=mult_poly_coeff_ppxf, $
				    bestfit_ppxf=bestfit_ppxf, chi2_ppxf=chi2_ppxf, $
				    stellar_kinematics_fit=stellar_kinematics, $
				    stellar_kinematics_err=stellar_kinematics_err, $
				    obj_fit_mask_gndf=obj_fit_mask_gndf, $
				    weights_gndf=weights_gndf, $
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

		MDAP_READ_OUTPUT, execution_plan[i].ofile, header=header, $
				  obj_fit_mask_ppxf=obj_fit_mask_ppxf, weights_ppxf=weights_ppxf, $
				  add_poly_coeff_ppxf=add_poly_coeff_ppxf, $
				  mult_poly_coeff_ppxf=mult_poly_coeff_ppxf, $
				  bestfit_ppxf=bestfit_ppxf, chi2_ppxf=chi2_ppxf, $
				  stellar_kinematics_fit=stellar_kinematics, $
				  stellar_kinematics_err=stellar_kinematics_err, $
				  obj_fit_mask_gndf=obj_fit_mask_gndf, $
				  weights_gndf=weights_gndf, $
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
	    endelse

	    if ~keyword_set(quiet) then $
		print, 'PLAN '+MDAP_STC(i+1,/integer)+': BLOCK 5 ... DONE.'

	    ; END BLOCK 5 ------------------------------------------------------
	    ;-------------------------------------------------------------------

	    if MDAP_MORE_ANALYSIS_BLOCKS_TO_FINISH(perform_block, /spec_fit) eq 0 then $
		continue

	    ; BLOCK 6 ----------------------------------------------------------
	    ; Perform the spectral-index measurements
	    ;-------------------------------------------------------------------

	    ; Check if spectral indices should be measured
	    if perform_block.spec_indx eq 1 then begin

		if ~keyword_set(quiet) then $
		    print, 'PLAN '+MDAP_STC(i+1,/integer)+': BLOCK 6 ...'

		;---------------------------------------------------------------
		; Get the weights to use to combine the templates, the bestfit
		; to the spectrum, the emission-line model, and the pixel mask
		; of the fit.
		; TODO: This is based on the size of the weights_gndf vector!
		if (size(weights_gndf))[2] ne (size(tpl_flux))[1] then begin
		    print, 'using ppxf results'
		    weights = weights_ppxf
		    bestfit = bestfit_ppxf
		    sz = size(bin_flux)
		    eml_model = dblarr(sz[1], sz[2])	; No emission-line model
		    fit_mask = obj_fit_mask_ppxf
		endif else begin
		    print, 'using gndf results'
		    weights = weights_gndf
		    bestfit = bestfit_gndf
		    fit_mask = obj_fit_mask_gndf
		endelse

		;---------------------------------------------------------------
		; Get the resolution-matched, emission-free spectra to use for
		; the spectral index measurements, both template and object.
		; TODO: Have a version for this function?
		; BLOCK 6a
		if MDAP_CAN_USE_SINDX_IMAGES(execution_plan[i].ofile) eq 1 then begin
		    MDAP_DEFINE_OUTPUT, si_bin_wave=si_bin_wave, si_bin_flux=si_bin_flux, $
					si_bin_ivar=si_bin_ivar, si_bin_mask=si_bin_mask, $
					si_optimal_template=si_optimal_template, $
					si_broad_optimal_template=si_broad_optimal_template
		    MDAP_READ_OUTPUT, execution_plan[i].ofile, si_bin_wave=si_bin_wave, $
				      si_bin_flux=si_bin_flux, si_bin_ivar=si_bin_ivar, $
				      si_bin_mask=si_bin_mask, $
				      si_optimal_template=si_optimal_template, $
				      si_broad_optimal_template=si_broad_optimal_template
		endif else begin
		    MDAP_SPECTRAL_INDEX_MEASUREMENTS_SPECTRA, wave, sres, bin_flux, bin_ivar, $
							      bin_mask, tpl_out_fits, weights, $
							      stellar_kinematics, $
							abs_line_keys[execution_plan[i].abs_par],$
							      si_bin_wave, si_bin_flux, $
							      si_bin_ivar, si_bin_mask, $
							      si_optimal_template, $
							      si_broad_optimal_template, $
							      bestfit=bestfit, $
							      eml_model=eml_model, $
							      fit_mask=fit_mask, $
							      remove_outliers=remove_outliers, $
						moments=execution_plan[i].analysis_par.moments
		    MDAP_WRITE_OUTPUT, execution_plan[i].ofile, si_bin_wave=si_bin_wave, $
				       si_bin_flux=si_bin_flux, si_bin_ivar=si_bin_ivar, $
				       si_bin_mask=si_bin_mask, $
				       si_optimal_template=si_optimal_template, $
				       si_broad_optimal_template=si_broad_optimal_template
		endelse

		; Perform the measurements
		MDAP_SPECTRAL_INDEX_MEASUREMENTS, abs_par, si_bin_wave, si_bin_flux, $
						  si_bin_ivar, si_bin_mask, si_optimal_template, $
						  si_broad_optimal_template, stellar_kinematics, $
						  abs_line_indx_omitted, abs_line_indx, $
						  abs_line_indx_err, abs_line_indx_otpl, $
						  abs_line_indx_botpl, dbg=dbg

		; TODO: For now using ivar for error.  Use residual instead?

		; Write the data to the output file
		MDAP_WRITE_OUTPUT, execution_plan[i].ofile, abs_par=abs_par, $
				   abs_line_key=abs_line_keys[execution_plan[i].abs_par], $
				   abs_line_indx_omitted=abs_line_indx_omitted, $
				   abs_line_indx_val=abs_line_indx, $
				   abs_line_indx_err=abs_line_indx_err, $
				   abs_line_indx_otpl=abs_line_indx_otpl, $
				   abs_line_indx_botpl=abs_line_indx_botpl, quiet=quiet, $
				   /read_header

		if ~keyword_set(quiet) then $
		    print, 'PLAN '+MDAP_STC(i+1,/integer)+': BLOCK 6 ... DONE.'

	    endif


	    ; END computation of "model-independent" data products #############
	    ; ##################################################################


	endfor
	; END loop over plans --------------------------------------------------
	;-----------------------------------------------------------------------

	; close up
	if ~keyword_set(nolog) then $
	    MDAP_LOG_CLOSE, log_file_unit

END

;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------


