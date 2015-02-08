;+
; NAME:
;       MDAP_ANALYSIS_SETUP_BLOCK
;
; PURPOSE:
;   MAIN BLOCK FOR MANGA_DAP:
;
;   - Reads the template library, emission-line parameters, and
;     spectral-index parameters to use
;   - Determines which of the remaining blocks to perform
;   - Interpolates guess kinematics based on the provided analysis_prior
;     in the execution_plan (if necessary)
;
; CALLING SEQUENCE:
;       MDAP_ANALYSIS_SETUP_BLOCK, output_file_root, tpl_library_keys, emission_line_parameters, $
;                                  ems_vacuum_wave, absorption_line_parameters, abs_vacuum_wave, $
;                                  ndrp, execution_plan, bskyx, bskyy, eml_par, abs_par, $
;                                  perform_block, star_kin_interp, gas_kin_interp
;
; INPUTS:
;       output_file_root string
;               Root name for all DAP output files
;
;       tpl_library_keys strarr
;               Keywords for the available template libraries
;
;       emission_line_parameters strarr
;               Available emission-line parameter files.
;
;       ems_vacuum_wave intarr
;               Integer flags that the emission-line parameters are
;               provided in vacuum.
;
;       absorption_line_parameters strarr
;               Available spectral-index parameter files.
;
;       abs_vacuum_wave intarr
;               Integer flags that the spectral-index parameters are
;               provided in vacuum.
;
;       ndrp integer
;               Number of DRP spectra that have been read
;
;       execution_plan ExecutionPlan
;               ExecutionPlan object for this set of analyses.
;
;       bskyx dblarr
;               Fiducial on-sky X position of each spectrum to analyze.
;
;       bskyy dblarr
;               Fiducial on-sky Y position of each spectrum to analyze.
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;
; OUTPUT:
;       eml_par EmissionLine[E]
;               Array of EmissionLine structures used when
;               fitting/masking emission lines.
;
;       abs_par SpectralIndex[A]
;               Array of SpectralIndex structures used to measure the
;               spectral indices.
;
;       perform_block RequiredAnalysisBlock
;               Structure defining which analysis blocks are to be
;               performed.
;
; OPTIONAL OUTPUT:
;       star_kin_interp dblarr
;               Interpolated set of stellar kinematics based on a
;               provided analysis prior.
;
;       gas_kin_interp dblarr
;               Interpolated set of gas kinematics based on a provided
;               analysis prior.
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;       MDAP_SET_TPL_LIB_OUTPUT_FILE()
;       MDAP_READ_EMISSION_LINE_PARAMETERS()
;       AIRTOVAC
;       MDAP_ERASEVAR
;       MDAP_READ_ABSORPTION_LINE_PARAMETERS()
;       MDAP_ANALYSIS_BLOCKS_TO_PERFORM()
;       MDAP_INTERPOLATE_KINEMATICS
;
; INTERNAL SUPPORT ROUTINES:
;       @mdap_analysis_block_logic
;
; REVISION HISTORY:
;       01 Feb 2015: (KBW) Pulled from manga_dap.pro
;-
;------------------------------------------------------------------------------

;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
; This is the file (mdap_analysis_block_logic.pro) that has all the
; functions that decide which analysis blocks should be performed
@mdap_analysis_block_logic
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------

PRO MDAP_ANALYSIS_SETUP_BLOCK, $
        output_file_root, tpl_library_keys, emission_line_parameters, ems_vacuum_wave, $
        absorption_line_parameters, abs_vacuum_wave, ndrp, execution_plan, bskyx, bskyy, eml_par, $
        abs_par, perform_block, star_kin_interp, gas_kin_interp

        ;---------------------------------------------------------------
        ; Get the nominal name for the template library.  This is used
        ; to determine which blocks to execute (via perform_block).
        ; This needs to be here so that tpl_out_fits can be passed to
        ; MDAP_ANALYSIS_BLOCKS_TO_PERFORM, but it may not be used.
        if execution_plan.tpl_lib ne -1 then begin
            tpl_out_fits = MDAP_SET_TPL_LIB_OUTPUT_FILE(output_file_root, $
                                                        tpl_library_keys[execution_plan.tpl_lib])
        endif

        ;---------------------------------------------------------------
        ; Read the emission line file
        if execution_plan.ems_par ne -1 then begin
            eml_par = MDAP_READ_EMISSION_LINE_PARAMETERS( $
                                                emission_line_parameters[ execution_plan.ems_par ]) 
            if ems_vacuum_wave[execution_plan.ems_par] eq 0 then begin
                neml = n_elements(eml_par)
                for j=0,neml-1 do begin
                    AIRTOVAC, eml_par[j].lambda, lam_vac
                    eml_par[j].lambda = lam_vac
                endfor
            endif
        endif else $
            MDAP_ERASEVAR, eml_par  ; Make sure it doesn't exist if no ems_par is defined

;       print, n_elements(eml_par)
;       stop

        ;---------------------------------------------------------------
        ; Read the spectral-index parameters
        if execution_plan.abs_par ne -1 then begin
            abs_par = MDAP_READ_ABSORPTION_LINE_PARAMETERS( $
                                            absorption_line_parameters[ execution_plan.abs_par ])
            if abs_vacuum_wave[ execution_plan.abs_par ] eq 0 then begin
                nabs = n_elements(abs_par)              ; Number of absorption-line indices
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
            MDAP_ERASEVAR, abs_par  ; Make sure it doesn't exist if no abs_par is defined

        ;---------------------------------------------------------------
        ; Delete the exiting file if overwrite=1
;       if execution_plan.overwrite eq 1 and FILE_TEST(execution_plan.ofile) eq 1 then begin
        if execution_plan.overwrite eq 1 && FILE_TEST(execution_plan.ofile) eq 1 then begin
            print, 'Removing existing file: '+execution_plan.ofile
            FILE_DELETE, execution_plan.ofile
        endif

        ;---------------------------------------------------------------
        ; Setup which analysis blocks (3-6) to perform
        perform_block = MDAP_ANALYSIS_BLOCKS_TO_PERFORM(execution_plan, ndrp, tpl_out_fits, $
                                                        eml_par, abs_par)

        ;---------------------------------------------------------------
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
        if perform_block.eml_only eq 1 then begin
            print, '    Emission-line-only fit: YES'
        endif else $
            print, '    Emission-line-only fit: NO'
        if perform_block.spec_indx eq 1 then begin
            print, '    Spectral indices: YES'
        endif else $
            print, '    Spectral indices: NO'

        ;---------------------------------------------------------------
        ; Interpolate data from the prior if necessary
;       if (perform_block.bin eq 1 or perform_block.spec_fit eq 1) and $
;           strlen(execution_plan.analysis_prior) ne 0 then begin
        if (perform_block.bin eq 1 || perform_block.spec_fit eq 1) $
            && strlen(execution_plan.analysis_prior) ne 0 then begin
           print, 'Interpolating kinematics based on provided prior'
           MDAP_INTERPOLATE_KINEMATICS, execution_plan.analysis_prior, bskyx, bskyy, $
                                        star_kin_interp, /velocity, /sigma, /stellar
           MDAP_INTERPOLATE_KINEMATICS, execution_plan.analysis_prior, bskyx, bskyy, $
                                        gas_kin_interp, /velocity, /sigma
        endif else begin
            ; Erase any existing interpolations
            MDAP_ERASEVAR, star_kin_interp
            MDAP_ERASEVAR, gas_kin_interp
        endelse
END

