;+
; NAME:
;       MDAP_OUTPUT_TABLE_ROWS()
;       MDAP_OUTPUT_COLUMN_SIZE()
;       MDAP_OUTPUT_IMAGE_SIZE()
;       MDAP_OUTPUT_IMAGE_SIZE_1D()
;
;       MDAP_CAN_USE_STFIT_DATA()
;       MDAP_CAN_USE_SGFIT_DATA()
;       MDAP_CAN_USE_ELOFIT_DATA()
;       MDAP_CAN_USE_SINDX_DATA()
;       MDAP_CAN_USE_SINDX_IMAGES()
;
;       MDAP_ADD_STFIT
;       MDAP_ADD_SGFIT
;       MDAP_ADD_ELOFIT
;       MDAP_ADD_SINDX
;
;       MDAP_PERFORM_STON_BLOCK()
;       MDAP_PERFORM_BIN_BLOCK()
;
;       MDAP_NEW_FILE_BLOCKS
;       MDAP_ANALYSIS_BLOCKS_TO_PERFORM()
;
;       MDAP_ALL_ANALYSIS_BLOCKS_COMPLETED()
;       MDAP_MORE_ANALYSIS_BLOCKS_TO_FINISH()
;
; PURPOSE:
;       A set of functions (with ()) and procedures used to execute the
;       logic in deciding which blocks to perform in MANGA_DAP.  The
;       main function is MDAP_ANALYSIS_BLOCKS_TO_PERFORM().
;
; REVISION HISTORY:
;       04 Dec 2014: (KBW) Moved from the header of manga_dap.pro.
;-
;------------------------------------------------------------------------------

@mdap_set_output_file_cols

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
        unit = fxposit(file, exten)                             ; Set the header unit
        hdr = headfits(unit)                                    ; Read it
        free_lun, unit                                          ; Free and close up the LUN
        return, fxpar(hdr, 'NAXIS2')                            ; Read and return number of rows
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

        bin_dim = MDAP_OUTPUT_IMAGE_SIZE(file, 'FLUX')  ; Dimensions of the binned spectra
        ntpl = (MDAP_OUTPUT_IMAGE_SIZE(tpl_fits, 'FLUX'))[0]

        ; The STFIT table must be populated with the correct number of rows
        if MDAP_OUTPUT_TABLE_ROWS(file, 'STFIT') ne bin_dim[0] then $
            return, 0

        cols = MDAP_SET_STFIT_COLS()

        ; The columns of the STFIT table must have the correct size
        if MDAP_OUTPUT_COLUMN_SIZE(file, 'STFIT', cols[0]) ne ntpl then $
            return, 0
        if analysis_par.degree gt 0 then begin
            if MDAP_OUTPUT_COLUMN_SIZE(file, 'STFIT', cols[1]) ne analysis_par.degree then $
                return, 0
        endif
        if analysis_par.mdegree gt 0 then begin
            if MDAP_OUTPUT_COLUMN_SIZE(file, 'STFIT', cols[2]) ne analysis_par.mdegree then $
                return, 0
        endif
        if MDAP_OUTPUT_COLUMN_SIZE(file, 'STFIT', cols[3]) ne analysis_par.moments then $
            return, 0
        if MDAP_OUTPUT_COLUMN_SIZE(file, 'STFIT', cols[4]) ne analysis_par.moments then $
            return, 0

        ; The SMSK and SMOD images must be the same size as the binned spectra
        msk_dim = MDAP_OUTPUT_IMAGE_SIZE(file, 'SMSK')
;       if msk_dim[0] ne bin_dim[0] or msk_dim[1] ne bin_dim[1] then $
        if msk_dim[0] ne bin_dim[0] || msk_dim[1] ne bin_dim[1] then $
            return, 0
        mod_dim = MDAP_OUTPUT_IMAGE_SIZE(file, 'SMOD')
;       if mod_dim[0] ne bin_dim[0] or mod_dim[1] ne bin_dim[1] then $
        if mod_dim[0] ne bin_dim[0] || mod_dim[1] ne bin_dim[1] then $
            return, 0

        return, 1
END

PRO MDAP_ADD_STFIT, $
                file, tpl_fits, analysis_par, perform_block

        if MDAP_CAN_USE_STFIT_DATA(file, tpl_fits, analysis_par) eq 1 then begin
            perform_block.spec_fit = 0  ; Do not perform the spectral fitting
        endif else $
            perform_block.spec_fit = 1  ; Peform the spectral fitting block

        perform_block.ppxf_only = 1     ; Only perform PPXF
END


; Tests if the existing SGFIT data in the file can be used directly
; TODO: NUMBER OF MOMENTS IS HARDWIRED!!
FUNCTION MDAP_CAN_USE_SGFIT_DATA, $
                file, tpl_fits, eml_par, analysis_par

        bin_dim = MDAP_OUTPUT_IMAGE_SIZE(file, 'FLUX')          ; Dimensions of the binned spectra
        ntpl = (MDAP_OUTPUT_IMAGE_SIZE(tpl_fits, 'FLUX'))[0]    ; Number of templates
        neml = n_elements(eml_par)                              ; Number of emission lines

        ; The SGFIT table must have the correct number of rows (one per binned spectrum)
        if MDAP_OUTPUT_TABLE_ROWS(file, 'SGFIT') ne bin_dim[0] then $
            return, 0

        cols = MDAP_SET_SGFIT_COLS()

        ; The columns of the SGFIT table must have the correct size
        if MDAP_OUTPUT_COLUMN_SIZE(file, 'SGFIT', cols[0]) ne ntpl then $
            return, 0
        if analysis_par.mdegree gt 0 then begin
            if MDAP_OUTPUT_COLUMN_SIZE(file, 'SGFIT', cols[1]) ne analysis_par.mdegree then $
                return, 0
        endif
        if MDAP_OUTPUT_COLUMN_SIZE(file, 'SGFIT', cols[2]) ne 2 then $
            return, 0
        if MDAP_OUTPUT_COLUMN_SIZE(file, 'SGFIT', cols[3]) ne 2 then $
            return, 0
        if analysis_par.reddening_order gt 0 then begin
            if MDAP_OUTPUT_COLUMN_SIZE(file,'SGFIT', cols[4]) ne analysis_par.reddening_order then $
                return, 0
            if MDAP_OUTPUT_COLUMN_SIZE(file,'SGFIT', cols[6]) ne analysis_par.reddening_order then $
                return, 0
        endif
        if MDAP_OUTPUT_COLUMN_SIZE(file, 'SGFIT', cols[7]) ne neml then $
            return, 0
        if MDAP_OUTPUT_COLUMN_SIZE(file, 'SGFIT', cols[8]) ne neml then $
            return, 0
        if MDAP_OUTPUT_COLUMN_SIZE(file, 'SGFIT', cols[9]) ne neml then $
            return, 0
        if MDAP_OUTPUT_COLUMN_SIZE(file, 'SGFIT', cols[12]) ne neml then $
            return, 0
        if MDAP_OUTPUT_COLUMN_SIZE(file, 'SGFIT', cols[13]) ne neml then $
            return, 0
        if MDAP_OUTPUT_COLUMN_SIZE(file, 'SGFIT', cols[14]) ne neml then $
            return, 0
        if MDAP_OUTPUT_COLUMN_SIZE(file, 'SGFIT', cols[15]) ne neml then $
            return, 0
        dim = MDAP_OUTPUT_COLUMN_SIZE(file, 'SGFIT', cols[10])
;       if dim[0] ne neml or dim[1] ne 2 then $
        if dim[0] ne neml || dim[1] ne 2 then $
            return, 0
        dim = MDAP_OUTPUT_COLUMN_SIZE(file, 'SGFIT', cols[11])
;       if dim[0] ne neml or dim[1] ne 2 then $
        if dim[0] ne neml || dim[1] ne 2 then $
            return, 0

        ; The SGMSK, SGMOD, and ELMOD images must be the same size as the binned spectra
        msk_dim = MDAP_OUTPUT_IMAGE_SIZE(file, 'SGMSK')
;       if msk_dim[0] ne bin_dim[0] or msk_dim[1] ne bin_dim[1] then $
        if msk_dim[0] ne bin_dim[0] || msk_dim[1] ne bin_dim[1] then $
            return, 0
        mod_dim = MDAP_OUTPUT_IMAGE_SIZE(file, 'SGMOD')
;       if mod_dim[0] ne bin_dim[0] or mod_dim[1] ne bin_dim[1] then $
        if mod_dim[0] ne bin_dim[0] || mod_dim[1] ne bin_dim[1] then $
            return, 0
        elm_dim = MDAP_OUTPUT_IMAGE_SIZE(file, 'ELMOD')
;       if elm_dim[0] ne bin_dim[0] or elm_dim[1] ne bin_dim[1] then $
        if elm_dim[0] ne bin_dim[0] || elm_dim[1] ne bin_dim[1] then $
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
                perform_block.spec_fit = 0      ; Do not perform the spectral fitting
            endif else $
                perform_block.spec_fit = 1      ; Must perform the spectral fitting
        endif
        perform_block.ppxf_only = 0             ; Perform both PPXF and GANDALF
END


FUNCTION MDAP_CAN_USE_ELOFIT_DATA, $
                file

        emlo_par = MDAP_DEFINE_EMISSION_LINES_ENCI_BELFIORE()

        bin_dim = MDAP_OUTPUT_IMAGE_SIZE(file, 'FLUX')  ; Dimensions of the binned spectra
        neml = n_elements(emlo_par)                     ; Number of emission-lines to fit

        ; The ELOFIT table must have the correct number of rows (one per binned spectrum)
        if MDAP_OUTPUT_TABLE_ROWS(file, 'ELOFIT') ne bin_dim[0] then $
            return, 0

        cols = MDAP_SET_ELOFIT_COLS()

        ; The columns of the ELOFIT table must have the correct size
        if MDAP_OUTPUT_COLUMN_SIZE(file, 'ELOFIT', cols[0]) ne 2 then $
            return, 0
        if MDAP_OUTPUT_COLUMN_SIZE(file, 'ELOFIT', cols[1]) ne 2 then $
            return, 0
        if MDAP_OUTPUT_COLUMN_SIZE(file, 'ELOFIT', cols[2]) ne neml then $
            return, 0
        if MDAP_OUTPUT_COLUMN_SIZE(file, 'ELOFIT', cols[3]) ne neml then $
            return, 0
        if MDAP_OUTPUT_COLUMN_SIZE(file, 'ELOFIT', cols[4]) ne neml then $
            return, 0
        dim = MDAP_OUTPUT_COLUMN_SIZE(file, 'ELOFIT', cols[5])
;       if dim[0] ne neml or dim[1] ne 2 then $
        if dim[0] ne neml || dim[1] ne 2 then $
            return, 0
        dim = MDAP_OUTPUT_COLUMN_SIZE(file, 'ELOFIT', cols[6])
;       if dim[0] ne neml or dim[1] ne 2 then $
        if dim[0] ne neml || dim[1] ne 2 then $
            return, 0
        if MDAP_OUTPUT_COLUMN_SIZE(file, 'ELOFIT', cols[7]) ne neml then $
            return, 0
        if MDAP_OUTPUT_COLUMN_SIZE(file, 'ELOFIT', cols[8]) ne neml then $
            return, 0
        if MDAP_OUTPUT_COLUMN_SIZE(file, 'ELOFIT', cols[9]) ne neml then $
            return, 0
        if MDAP_OUTPUT_COLUMN_SIZE(file, 'ELOFIT', cols[10]) ne neml then $
            return, 0
        if MDAP_OUTPUT_COLUMN_SIZE(file, 'ELOFIT', cols[11]) ne neml then $
            return, 0
        
        if MDAP_OUTPUT_COLUMN_SIZE(file, 'ELOFIT', cols[12]) ne 2 then $
            return, 0
        if MDAP_OUTPUT_COLUMN_SIZE(file, 'ELOFIT', cols[13]) ne 2 then $
            return, 0
        if MDAP_OUTPUT_COLUMN_SIZE(file, 'ELOFIT', cols[14]) ne neml then $
            return, 0
        if MDAP_OUTPUT_COLUMN_SIZE(file, 'ELOFIT', cols[15]) ne neml then $
            return, 0
        if MDAP_OUTPUT_COLUMN_SIZE(file, 'ELOFIT', cols[16]) ne neml then $
            return, 0
        dim = MDAP_OUTPUT_COLUMN_SIZE(file, 'ELOFIT', cols[17])
;       if dim[0] ne neml or dim[1] ne 2 then $
        if dim[0] ne neml || dim[1] ne 2 then $
            return, 0
        dim = MDAP_OUTPUT_COLUMN_SIZE(file, 'ELOFIT', cols[18])
;       if dim[0] ne neml or dim[1] ne 2 then $
        if dim[0] ne neml || dim[1] ne 2 then $
            return, 0
        if MDAP_OUTPUT_COLUMN_SIZE(file, 'ELOFIT', cols[19]) ne neml then $
            return, 0
        if MDAP_OUTPUT_COLUMN_SIZE(file, 'ELOFIT', cols[20]) ne neml then $
            return, 0
        if MDAP_OUTPUT_COLUMN_SIZE(file, 'ELOFIT', cols[21]) ne neml then $
            return, 0
        if MDAP_OUTPUT_COLUMN_SIZE(file, 'ELOFIT', cols[22]) ne neml then $
            return, 0
        if MDAP_OUTPUT_COLUMN_SIZE(file, 'ELOFIT', cols[23]) ne neml then $
            return, 0
        
        ; The ELOMEW and ELOMFB images must be the same size as the binned spectra
        msk_dim = MDAP_OUTPUT_IMAGE_SIZE(file, 'ELOMEW')
;       if msk_dim[0] ne bin_dim[0] or msk_dim[1] ne bin_dim[1] then $
        if msk_dim[0] ne bin_dim[0] || msk_dim[1] ne bin_dim[1] then $
            return, 0
        mod_dim = MDAP_OUTPUT_IMAGE_SIZE(file, 'ELOMFB')
;       if mod_dim[0] ne bin_dim[0] or mod_dim[1] ne bin_dim[1] then $
        if mod_dim[0] ne bin_dim[0] || mod_dim[1] ne bin_dim[1] then $
            return, 0

        return, 1
END


; If adding the spectral index measurements, the file must have:
;       - the template weights (need the same # of templates based on the provided template file)
;       - the stellar kinematics (need the same number of binned spectra as in the file)
;               THIS TEST CURRENTLY ONLY RUNS IF THE BINNED SPECTRA IN THE FILE ARE TO BE USED!
;       - the best fitting spectra (must be same size as binned spectra)
PRO MDAP_ADD_ELOFIT, $
                file, perform_block

        if MDAP_CAN_USE_ELOFIT_DATA(file) eq 1 then begin
            perform_block.eml_only = 0         ; Do not perform the emission-line-only fits
        endif else $
            perform_block.eml_only = 1
END


; TODO: Add a check for the size of the images (SIWAVE, SIFLUX, etc)?
FUNCTION MDAP_CAN_USE_SINDX_DATA, $
                file, abs_par

        bin_dim = MDAP_OUTPUT_IMAGE_SIZE(file, 'FLUX')  ; Dimensions of the binned spectra
        nabs = n_elements(abs_par)                      ; Number of spectral indices

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

        if wave_dim eq 1 then $                         ; No wave image
            return, 0

        si_dim = MDAP_OUTPUT_IMAGE_SIZE(file, 'FLUX')   ; Dimensions of the binned spectra
        si_dim[1] = wave_dim                            ; Force spectrum length

        ; The SIFLUX, SIIVAR, SIMASK, SIOTPL, and SIBOTPL images must have the
        ; same and correct dimensions.
        flx_dim = MDAP_OUTPUT_IMAGE_SIZE(file, 'SIFLUX')
;       if flx_dim[0] ne si_dim[0] or flx_dim[1] ne si_dim[1] then $
        if flx_dim[0] ne si_dim[0] || flx_dim[1] ne si_dim[1] then $
            return, 0
        ivr_dim = MDAP_OUTPUT_IMAGE_SIZE(file, 'SIIVAR')
;       if ivr_dim[0] ne si_dim[0] or ivr_dim[1] ne si_dim[1] then $
        if ivr_dim[0] ne si_dim[0] || ivr_dim[1] ne si_dim[1] then $
            return, 0
        msk_dim = MDAP_OUTPUT_IMAGE_SIZE(file, 'SIMASK')
;       if msk_dim[0] ne si_dim[0] or msk_dim[1] ne si_dim[1] then $
        if msk_dim[0] ne si_dim[0] || msk_dim[1] ne si_dim[1] then $
            return, 0
        otp_dim = MDAP_OUTPUT_IMAGE_SIZE(file, 'SIOTPL')
;       if otp_dim[0] ne si_dim[0] or otp_dim[1] ne si_dim[1] then $
        if otp_dim[0] ne si_dim[0] || otp_dim[1] ne si_dim[1] then $
            return, 0
        btp_dim = MDAP_OUTPUT_IMAGE_SIZE(file, 'SIBOTPL')
;       if btp_dim[0] ne si_dim[0] or btp_dim[1] ne si_dim[1] then $
        if btp_dim[0] ne si_dim[0] || btp_dim[1] ne si_dim[1] then $
            return, 0

        return, 1
END     


; If adding the spectral index measurements, the file must have:
;       - the template weights (need the same # of templates based on the provided template file)
;       - the stellar kinematics (need the same number of binned spectra as in the file)
;               THIS TEST CURRENTLY ONLY RUNS IF THE BINNED SPECTRA IN THE FILE ARE TO BE USED!
;       - the best fitting spectra (must be same size as binned spectra)
PRO MDAP_ADD_SINDX, $
                file, abs_par, perform_block

        if MDAP_CAN_USE_SINDX_DATA(file, abs_par) eq 1 then begin
            perform_block.spec_indx = 0         ; Do not perform the spectral-index measurements
        endif else $
            perform_block.spec_indx = 1
END


; Set which analyses to perform for a new file
PRO MDAP_NEW_FILE_BLOCKS, $
                analysis, perform_block

        ; TODO: Somehow put this in a common area with MDAP_BUILD_EXECUTION_PLANS
        perform_block.ston = 1          ; Always perform the S/N block
        perform_block.bin = 1           ; Always perform the spatial binning

        ; Order is important for the first two if statements because of ppxf_only!
        if analysis[0] eq 1 then begin  ; Wants the stellar continuum fitting
            perform_block.spec_fit = 1          ; Perform the spectral fit
            perform_block.ppxf_only = 1         ; Only ppxf
        endif

        ; For a new file, the star+gas fits, REQUIRES the stellar continuum fit
        if analysis[1] eq 1 then begin  ; Wants the star+gas fitting
            perform_block.spec_fit = 1          ; Perform the spectral fit
            perform_block.ppxf_only = 0         ; Not only ppxf
        endif

        ; The emission-line only fitting can effectively be done
        ; independent of all the others (even if this isn't the best
        ; thing to do)
        if analysis[2] eq 1 then begin  ; Wants the emission-line-only fitting
            perform_block.eml_only = 1
        endif

        ; For a new file, the spectral index fits, REQUIRES the spectral (continuum) fitting
        if analysis[3] eq 1 then begin  ; Wants the spectral indices
            perform_block.spec_fit = 1                          ; Perform the spectral fit
            perform_block.ppxf_only = analysis[1] eq 0 ? 1 : 0  ; Check if ppxf only
            perform_block.spec_indx = 1                         ; Do the spectral index measurements
        endif

        ; If it reaches here, only the S/N calculation and the binning are performed
END

; TODO: Should only reach here if the output file already exists!
FUNCTION MDAP_PERFORM_STON_BLOCK, $
                execution_plan, ndrp

        ; Conditions to return true:
        ;       A: overwrite is true
        ;       B: output DRPS extension does not have the correct size

        ; TODO: Should also eventually check (?):
        ;       - the spaxel size
        ;       - the S/N wavelength range
        ;       - the on-sky positions

        ; CASE A ---------------------------------------------------------------
        if execution_plan.overwrite eq 1 then $
            return, 1

        ; CASE B ---------------------------------------------------------------
        ;   MDAP_CHECK_OUTPUT_FILE has already:
        ;       - Checked for the existence of the DRPS extension
        ;       - Checked that all the appropriate columns exist
        ;   Need to:
        ;       - Check the number of rows against the expected value
        if MDAP_OUTPUT_TABLE_ROWS(execution_plan.ofile, 'DRPS') ne ndrp then $
            return, 1

        return, 0               ; Do not need to redo the S/N calculation
END

; TODO: Should only reach here if the output file already exists!
FUNCTION MDAP_PERFORM_BIN_BLOCK, $
                execution_plan, ndrp

        ; Conditions to return true:
        ;       A: overwrite is true
        ;       B: output DRPS extension does not have the correct size
        ;       C: output BINS extension does not have the correct size

        ; CASE A ---------------------------------------------------------------
        if execution_plan.overwrite eq 1 then $
            return, 1

        ;   MDAP_CHECK_OUTPUT_FILE has already:
        ;       - Checked for the existence of the DRPS and BINS extension
        ;       - Checked that all the appropriate columns exist
        ;   Need to:

        ; CASE B ---------------------------------------------------------------
        ;       - Check the number of DRPS rows against the expected value
        if MDAP_OUTPUT_TABLE_ROWS(execution_plan.ofile, 'DRPS') ne ndrp then $
            return, 1

        ; CASE C ---------------------------------------------------------------
        ;       - Check the number of BINS rows against the expected value
        nbin = (MDAP_OUTPUT_IMAGE_SIZE(execution_plan.ofile, 'FLUX'))[0]
        if MDAP_OUTPUT_TABLE_ROWS(execution_plan.ofile, 'BINS') ne nbin then $
            return, 1

        return, 0               ; Do not need to redo the spatial binning
END

; TODO: I want to be able to add extensions to a file for new blocks, without
; having to redo old ones.  Will it do this already?
FUNCTION MDAP_ANALYSIS_BLOCKS_TO_PERFORM, $
                execution_plan, ndrp, tpl_fits, eml_par, abs_par

        perform_block = { RequiredAnalysisBlock, ston:0, bin:0, spec_fit:0, ppxf_only:0, $
                                                 eml_only:0, spec_indx:0 }

        ; Determine if the output file already exists
        file_exists = file_test(execution_plan.ofile)

        ; - If the file exists but does not look like a DAP file and the
        ;   overwrite flag has been set to false, throw an error
;       if file_exists eq 1 and execution_plan.overwrite eq 0 then begin
        if file_exists eq 1 && execution_plan.overwrite eq 0 then begin
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
        ;       A: overwrite is true
        ;       B: output DRPS extension does not have the correct size
        ;       C: not determined by MDAP_PERFORM_STON_BLOCK(), but below
        perform_block.ston = MDAP_PERFORM_STON_BLOCK(execution_plan, ndrp)
        
        ; Perform the BIN block if:
        ;       A; overwrite is true
        ;       B: output DRPS extension does not have the correct size
        ;       C: output BINS extension does not have the correct size
        perform_block.bin = MDAP_PERFORM_BIN_BLOCK(execution_plan, ndrp)

        ; This is redundant with the next if statement
;       if perform_block.bin eq 1 and perform_block.ston eq 0 then $
;           perform_block.ston = 1                      ; Always calculate S/N if binning

        if perform_block.bin eq 1 then begin            ; If rebinning, treat like a new file
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

        if execution_plan.analysis[2] eq 1 then begin   ; Wants the emission-line-only fits
            MDAP_ADD_ELOFIT, execution_plan.ofile, perform_block
        endif

        if execution_plan.analysis[3] eq 1 then begin   ; Wants the spectral indices
            ; First check that at least the STFIT results are available
;           if execution_plan.analysis[0] eq 0 and execution_plan.analysis[1] eq 0 then begin
            if execution_plan.analysis[0] eq 0 && execution_plan.analysis[1] eq 0 then begin
                MDAP_ADD_STFIT, execution_plan.ofile, tpl_fits, execution_plan.analysis_par, $
                                perform_block
            endif
            MDAP_ADD_SINDX, execution_plan.ofile, abs_par, perform_block
        endif

        return, perform_block                           ; All blocks should have been completed
END

FUNCTION MDAP_ALL_ANALYSIS_BLOCKS_COMPLETED, $
                perform_block

;       if perform_block.ston eq 0 and $
;          perform_block.bin eq 0 and $
;          perform_block.spec_fit eq 0 and $
;          perform_block.eml_only eq 0 and $
;          perform_block.spec_indx eq 0 then begin
        if perform_block.ston eq 0 && $
           perform_block.bin eq 0 && $
           perform_block.spec_fit eq 0 && $
           perform_block.eml_only eq 0 && $
           perform_block.spec_indx eq 0 then begin
            return, 1                                   ; All the blocks have been completed
        endif

        return, 0                                       ; Blocks need to be done
END

FUNCTION MDAP_MORE_ANALYSIS_BLOCKS_TO_FINISH, $
                perform_block, ston=ston, bin=bin, spec_fit=spec_fit, eml_only=eml_only, $
                spec_indx=spec_indx
        if keyword_set(ston) then begin
;           if perform_block.bin eq 0 and $
;              perform_block.spec_fit eq 0 and $
;              perform_block.eml_only eq 0 and $
;              perform_block.spec_indx eq 0 then begin
            if perform_block.bin eq 0 && $
               perform_block.spec_fit eq 0 && $
               perform_block.eml_only eq 0 && $
               perform_block.spec_indx eq 0 then begin
                return, 0
            endif
        endif

        if keyword_set(bin) then begin
;           if perform_block.spec_fit eq 0 and $
;              perform_block.eml_only eq 0 and $
;              perform_block.spec_indx eq 0 then begin
            if perform_block.spec_fit eq 0 && $
               perform_block.eml_only eq 0 && $
               perform_block.spec_indx eq 0 then begin
                return, 0
            endif
        endif

        if keyword_set(spec_fit) then begin
;           if perform_block.eml_only eq 0 and $
;              perform_block.spec_indx eq 0 then begin
            if perform_block.eml_only eq 0 && $
               perform_block.spec_indx eq 0 then begin
                return, 0
            endif
        endif

        if keyword_set(eml_only) then begin
            if perform_block.spec_indx eq 0 then begin
                return, 0
            endif
        endif

        if keyword_set(spec_indx) then $
            return, 0

        return, 1                                       ; Blocks need to be done
END
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------



