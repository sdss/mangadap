;+
; NAME:
;       MDAP_TEMPLATE_LIBRARY_BLOCK
;
; PURPOSE:
;   MAIN BLOCK FOR MANGA_DAP:
;
;   Builds the template libraries that are required for ALL the
;   execution plans to be applied in this run of the DAP.  This includes
;   libraries with their resolution and sampling matched to the DRP
;   object spectra and to the spectral-index systems.
;   
; CALLING SEQUENCE:
;       MDAP_TEMPLATE_LIBRARY_BLOCK, manga_dap_version, output_file_root, n_tpl_lib, $
;                                    tpl_library_keys, template_libraries, tpl_vacuum_wave, $
;                                    abs_line_keys, velocity_initial_guess, wave, sres, $
;                                    execution_plan, velScale, /nolog, $
;                                    log_file_unit=log_file_unit, datacube_name=datacube_name, $
;                                    quiet=quiet
;
; INPUTS:
;       manga_dap_version MaNGADAPVersion
;               Structure used to keep track of various
;               version-controlled procedures in the DAP.
;
;       output_file_root string
;               Root name for all DAP output files from this run.
;
;       n_tpl_lib integer
;               Total number of template libraries available.
;       
;       tpl_library_keys strarr[]
;               Set of keywords used to uniquely select template
;               libraries.
;
;       template_libraries strarr[]
;               Path definitions for each of the libraries.
;
;       tpl_vacuum_wave intarr[]
;               Flag that wavelengths of each template library is (1) or
;               is not (0) in vacuum.
;
;       abs_line_keys strarr[]
;               Set of keywords used to uniquely select spectral-index
;               parameter sets.
;
;       velocity_initial_guess double
;               Initial guess for the velocity when fitting the
;               kinematics.
;
;       wave dblarr[T]
;               Wavelength of each spectral channel T.
;       
;       sres dblarr[T]
;               Spectral resolution at each wavelength channel T.
;
;       execution_plan ExecutionPlan[]
;               Array of structures providing a set of plans and
;               parameters as needed for each of the requested series of
;               analysis steps.
;
; OPTIONAL INPUTS:
;       log_file_unit LUN
;               File unit pointing to the log file
;
;       datacube_name string
;               Name of the DRP (RSS or CUBE) fits file to read.
;               Optional because this is needed if it is to be written to
;               the log file.
;
; OPTIONAL KEYWORDS:
;       /nolog
;               Suppress output to the log file.
;
;       /quiet
;               Suppress output to the screen.
;
; OUTPUT:
;       velScale double
;               Velocity scale of the logarithmically binned pixels in
;               the wave vector.  TODO: This is a little out of place as
;               output here.  C'est la vie.
;
; OPTIONAL OUTPUT:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;       MDAP_VELOCITY_SCALE()
;       MDAP_CREATE_DAP_TEMPLATE_LIBRARY
;       MDAP_SET_TPL_LIB_OUTPUT_FILE()
;       MDAP_ABS_LINE_RESOLUTION()
;
; REVISION HISTORY:
;       01 Feb 2015: (KBW) Pulled from manga_dap.pro
;-
;------------------------------------------------------------------------------

PRO MDAP_TEMPLATE_LIBRARY_BLOCK, $
            manga_dap_version, output_file_root, n_tpl_lib, tpl_library_keys, template_libraries, $
            tpl_vacuum_wave, abs_line_keys, velocity_initial_guess, wave, sres, execution_plan, $
            velScale, nolog=nolog, log_file_unit=log_file_unit, datacube_name=datacube_name, $
            quiet=quiet

        ; Number of execution plans
        n_plans = n_elements(execution_plan)

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
        indx=where(use_tpl_lib eq 1, count)
        if count ne 0 then begin

            MDAP_CREATE_DAP_TEMPLATE_LIBRARY, version=manga_dap_version.tpl_lib_setup

            for i=0,n_tpl_lib-1 do begin
                if use_tpl_lib[i] ne 1 then $
                    continue

                ; TODO: Need to rethink usage of setup procedure to define library keyword?

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
                        if n_elements(datacube_name) ne 0 then begin
                            printf, log_file_unit, '[INFO] Resolution and sampling matched to: ' + $
                                    datacube_name
                        endif
                        printf, log_file_unit,'[INFO] Template library setup done using version: ' $
                                + manga_dap_version.tpl_lib_setup
                        printf, log_file_unit, '[INFO] Results written to: ' + tpl_out_fits
                    endif

                endif

                ; Have to generate a spectral-index template set for
                ; each template library / index system combination
                ; (requested by the plans)
                for j=0,n_plans-1 do begin

                    ; Spectral-index analysis not performed so continue
                    if execution_plan[j].abs_par eq -1 then $
                        continue

                    ; This execution plan doesn't use this template library
                    if execution_plan[j].tpl_lib ne i then $
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
                                + manga_dap_version.tpl_lib_setup
                        printf, log_file_unit, '[INFO] Results written to: ' + tpl_out_fits
                    endif

                endfor  ; Loop over execution plans
            endfor      ; Loop over template libraries
        endif           ; If any template libraries are used

END

