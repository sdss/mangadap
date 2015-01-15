;+
; NAME:
;       MANGA_DAP
;
; PURPOSE:
;       Main wrapping routine for MaNGA data analysis pipeline.
;
;       For now, see documentation here:
;               https://trac.sdss.org/wiki/MANGA/Data/DAPDevelopment/DAP_v0_9_summary
;
; CALLING SEQUENCE:
;       MANGA_DAP, inptbl=inptbl, index=index, drppath=drppath, dappath=dappath, dapsrc=dapsrc, $
;                  /nolog, /quiet, /plot, /dbg
;
;           or
;
;       MANGA_DAP, par=par, drppath=drppath, dappath=dappath, dapsrc=dapsrc, /nolog, /quiet, $
;                  /plot, /dbg
;
; INPUTS:
;
; OPTIONAL INPUTS:
;       inptbl string
;           File with the input table providing the necessary input
;           parameters.  See MDAP_CREATE_INPUT_TABLE for the file
;           format.  Should not be provided along with par.  Must also
;           provide index for use.
;
;       index integer
;           Row index in inptbl with the parameters to use.
;
;       par string
;           SDSS par (yanny) file with the parameters to use.  Read
;           using YANNY_READONE(), meaning that only the first entry
;           will be read, but others can be present.
;
;       drppath string
;           Path to the DRP-produced file.  The name of the file MUST be
;           manga-[plate]-[ifudesign]-LOG[mode].fits (or .fits.gz).
;           This optional sets the path to this file.  See
;           MDAP_READ_INPUT_SETUP for the default path.
;
;       dappath string
;           Path to the root directory for the DAP output files.  The
;           output files will be placed in dappath/plate/ifudesign/ with
;           each file having the root name
;           manga-[plate]-[ifudesign]-LOG[mode]_; see MDAP_SETUP_IO.
;           See MDAP_READ_INPUT_SETUP for the default path.
;
;       dapsrc string
;           Path the the DAP source directory.  Used in
;           MDAP_EXECUTION_SETUP to set the paths to the template
;           libraries, etc.
;
; OPTIONAL KEYWORDS:
;       /nolog
;               Do not produce the log file.
;
;       /quiet
;               Limit (eliminate?) the output printed to the screen.
;
;       /plot
;               Produce the Voronoi binning, PPXF and GANDALF plots
;
;       /dbg
;               Only attempts to fit the first spectrum of the first
;               execution plan as a test run mostly used for debugging
;
; OUTPUT:
;
; OPTIONAL OUTPUT:
;
; TODO:
;       - Include MW extinction curve
;       - Somehow allow to check for overlapping objects (multiple redshifts in
;         a single IFU)
;       - If the execution_plan is supposed to overwrite old data,
;         actually start by removing the old file
;
; COMMENTS:
;       At some point I (KBW) switched from using 'absorption-line indices' to
;       'spectral indices', primarily because the spectral-index analysis also
;       measures the strengths of spectral breaks/bandheads like D4000.  This
;       just a choice of nomenclature.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED (OUT-OF-DATE):
;       MDAP_EXECUTION_SETUP
;       READCOL
;       MDAP_SETUP_IO
;       MDAP_DRP_CHECK_FILETYPE
;       MDAP_CHECK_FILE_EXISTS()
;       MDAP_BUILD_EXECUTION_PLANS
;       MDAP_PRINT_EXECUTION_PLAN
;       SXADDPAR
;       MDAP_STC()
;       MDAP_READ_DRP_FITS
;       MDAP_GET_SPAXEL_SIZE
;       MDAP_FIDUCIAL_BIN_XY
;       MDAP_SET_TPL_LIB_OUTPUT_FILE()
;       MDAP_READ_TEMPLATE_LIBRARY
;       MDAP_MATCH_RESOLUTION_TPL2OBJ
;       MDAP_TRIM_SPECTRAL_RANGE
;       MDAP_VELOCITY_SCALE()
;       MDAP_RESAMPLE_TEMPLATES
;       MDAP_NORMALIZE_TEMPLATES
;       WRITEFITS
;       MDAP_ADD_FITS_LAYER
;       MDAP_SELECT_GOOD_SPECTRA
;       MDAP_SELECT_WAVE
;       MDAP_CALCULATE_SN
;       MDAP_WRITE_OUTPUT
;       MDAP_SPATIAL_BINNING
;       MDAP_READ_EMISSION_LINE_PARAMETERS
;       MDAP_ERASEVAR
;       MDAP_SPECTRAL_FITTING
;
; INTERNAL SUPPORT ROUTINES (OUT-OF-DATE):
;       MDAP_FILE_EXTENSION_DATE()
;       MDAP_LOG_INSTANTIATE
;       MDAP_LOG_CLOSE
;       MDAP_GENERATE_OUTPUT_FILE_NAMES
;       MDAP_SET_TPL_LIB_OUTPUT_FILE
;
; REVISION HISTORY:
;       01 Sep 2014: Copied from v0_8 version by L. Coccato.
;       02 Sep 2014: (KBW) Formating and minor edits
;       03 Sep 2014: (KBW) Basic compilation errors
;       26 Oct 2014: (KBW) v0.9 (see documentation linked to above)
;       11 Nov 2014: (KBW) Inclusion of spectral index measurements, addition of
;                          infrastructure for used to check for pre-existing
;                          analysis results.  The latter is a bit messy.  Will
;                          hopefully be cleaned up when we move to python.
;       28 Nov 2014: (KBW) Allow for SDSS par file input; specify the
;                          input table using the inptbl parameter
;                          instead of in the MDAP_EXECUTION_SETUP
;                          procecure; allow the paths for the drp and
;                          dap to be changed from the default values
;                          (which are now setup in
;                          MDAP_READ_INPUT_SETUP).
;       04 Dec 2014: (KBW) Accommodate change in MDAP_FIDUCIAL_BIN_XY;
;                          Changes to bin_par to allow for the radial
;                          binning (MDAP_BUILD_EXECUTION_PLAN,
;                          MDAP_PRINT_EXECUTION_PLAN,
;                          MDAP_EXECUTION_SETUP); Logic used to
;                          determine which blocks to perform moved to
;                          mdap_analysis_block_logic.pro, include using
;                          the @ convention just before MANGA_DAP);
;                          Changes to the input of MDAP_SPATIAL_BINNING
;                          due to changes to bin_par.
;       05 Dec 2014: (KBW) Allow an analysis prior -> changes to the
;                          ExecutionPlan structure.  Moved
;                          MDAP_GENERATE_OUTPUT_FILE_NAMES to be part of
;                          building the execution plans.
;       09 Jan 2015: (KBW) Add instrumental dispersion calculation for
;                          star+gas (GANDALF) results.
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
;       printf, file_unit, 'Data products:'
;       printf, file_unit, '    Primary fits file: ', output_filefits
;       printf, file_unit, ''

PRO MDAP_LOG_CLOSE, $
                file_unit

        printf, file_unit, 'End systime(): ', systime()
        close, file_unit
        free_lun, file_unit
END

;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------


;-------------------------------------------------------------------------------
; Set the guess kinematics
PRO MDAP_INITIALIZE_GUESS_KINEMATICS, $
                nb, analysis_prior, star_kin_interp, gas_kin_interp, bin_indx, $
                velocity_initial_guess, velocity_dispersion_initial_guess, star_kin_guesses, $
                gas_kin_guesses

        star_kin_guesses = dblarr(nb, 4)
        gas_kin_guesses = dblarr(nb, 2)
        if strlen(analysis_prior) ne 0 then begin

            ; TODO: Add some limits for the guess kinematics to make
            ; sure they're not too crazy

            for j=0,nb-1 do begin
                indx = where(bin_indx eq j, count)
;               if indx[0] ne -1 then begin
                if count ne 0 then begin
                    star_kin_guesses[j,0] = mean(star_kin_interp[indx,0])
                    star_kin_guesses[j,1] = mean(star_kin_interp[indx,1])
                    gas_kin_guesses[j,0] = mean(gas_kin_interp[indx,0])
                    gas_kin_guesses[j,1] = mean(gas_kin_interp[indx,1])
                endif else begin
                    star_kin_guesses[j,0] = velocity_initial_guess
                    star_kin_guesses[j,1] = velocity_dispersion_initial_guess
                    gas_kin_guesses[j,0] = velocity_initial_guess
                    gas_kin_guesses[j,1] = 50.0
                endelse
            endfor
        endif else begin
            ; Set the starting guesses for the kinematics to a
            ; single value
            star_kin_guesses[*,0] = velocity_initial_guess              ; stellar velocity
            star_kin_guesses[*,1] = velocity_dispersion_initial_guess   ; stellar sigma
            gas_kin_guesses[*,0] = velocity_initial_guess               ; gas velocity
            gas_kin_guesses[*,1] = 50.                                  ; gas sigma
        endelse
END
;-------------------------------------------------------------------------------


;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
; Sets the template file names

FUNCTION MDAP_SET_TPL_LIB_OUTPUT_FILE, $
                file_root, library_key, abs_line_key=abs_line_key
        if n_elements(abs_line_key) eq 0 then $
            return, file_root+library_key+'.fits'

        return, file_root+library_key+'_'+abs_line_key+'.fits'
END
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------

@mdap_analysis_block_logic

;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
; This is the main wrapper procedure for the DAP.

PRO MANGA_DAP, $
        inptbl=inptbl, index=index, par=par, drppath=drppath, dappath=dappath, dapsrc=dapsrc, $
        nolog=nolog, quiet=quiet, plot=plot, dbg=dbg

        print, 'entered manga_dap'

        manga_dap_version = '0.93'      ; set the version number
        t0=systime(/seconds)/60./60.    ; start time

        ; BLOCK 0 --------------------------------------------------------------
        ; Initialization
        ;       - File I/O setup
        ;       - Execution plan setup
        ;       - Version control
        ;-----------------------------------------------------------------------
        if ~keyword_set(quiet) then $
            print, 'BLOCK 0 ... '
        ; TODO: Put all this initialization stuff in a wrapper?


        ; Execute the script the sets the user-level execution parameters
        ; total_filelist, output_root_dir, 
        MDAP_EXECUTION_SETUP, tpl_library_keys, template_libraries, tpl_vacuum_wave, $
                              ems_line_keys, emission_line_parameters, ems_vacuum_wave, $
                              abs_line_keys, absorption_line_parameters, abs_vacuum_wave, $
                              signifier, bin_par, w_range_sn, threshold_ston_bin, $
                              w_range_analysis, threshold_ston_analysis, analysis, $
                              tpl_lib_analysis, ems_par_analysis, abs_par_analysis, $
                              analysis_par, analysis_prior, overwrite_flag, dapsrc=dapsrc, $
                              save_intermediate_steps=save_intermediate_steps, $
                              remove_null_templates=remove_null_templates, $
                              external_library=external_library

;        if file_test(total_filelist) eq 0 then $
;            message, total_filelist+' does not exist!'

        ; Check the mdap_execution_setup results
        if n_elements(external_library) then begin
            bvls_shared_lib=external_library+'bvls.so'
            if file_test(bvls_shared_lib) eq 0 then $
                bvls_shared_lib=external_library+'bvls.dylib'
            if file_test(bvls_shared_lib) eq 0 then begin
                print, 'Shared object library does not exist! Continuing with internal routines.'
                ; Undefine the variables
                tempvar = size(temporary(external_library))
                tempvar = size(temporary(bvls_shared_lib))
            endif
        endif

        ; Read the input parameters
        MDAP_READ_INPUT_SETUP, inptbl=inptbl, index=index, par=par, drppath=drppath, $
                               dappath=dappath, plate, ifudesign, mode, velocity_initial_guess, $
                               velocity_dispersion_initial_guess, ell, pa, Reff, root_name, $
                               output_root_dir

        ; Set input file names and output directory
        ;     make the directory if it doesn't exist
        ; TODO: Allow for a single spectrum mode
        MDAP_SETUP_IO, root_name, output_root_dir, datacube_name, file_root, output_dir, $
                       output_file_root

        ; TODO: save_intermediate_steps and remove_null_templates are not used
        ; Check if save_intermediate_steps exits, if not set to not save intermediate steps
        if n_elements(save_intermediate_steps) eq 0 then $
            save_intermediate_steps = 0

        ; Check if remove_null_elemens exits, if not set to remove null templates
        ; TODO: Where are they removed from?
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
        ; TODO: Version control will eventually happen in
        ; MDAP_GENERATE_OUTPUT_FILE_NAMES, or preferrably in a more well-named
        ; procedure
       
        mdap_execute_plan_version = '0'
        mdap_read_drp_fits_version = '0'
        mdap_fiducial_bin_xy_version = '0'
        mdap_tpl_lib_setup_version = '0'
        mdap_calculate_sn_version = '0'
        mdap_spatial_binning_version = '0'
        mdap_spectral_fitting_version = '0'
        mdap_spectral_index_version = '0'

        MDAP_BUILD_EXECUTION_PLANS, version=mdap_execute_plan_version
        MDAP_READ_DRP_FITS, version=mdap_read_drp_fits_version
        MDAP_FIDUCIAL_BIN_XY, version=mdap_fiducial_bin_xy_version
        MDAP_CREATE_DAP_TEMPLATE_LIBRARY, version=mdap_tpl_lib_setup_version
        MDAP_CALCULATE_SN, version=mdap_calculate_sn_version
        MDAP_SPATIAL_BINNING, version=mdap_spatial_binning_version
        MDAP_SPECTRAL_FITTING, version=mdap_spectral_fitting_version
        MDAP_SPECTRAL_INDEX_MEASUREMENTS, version=mdap_spectral_index_version

        ;-----------------------------------------------------------------------


        ;-----------------------------------------------------------------------
        ; TODO: Same as ntpl_libraries, neml_files, and nabs_files in mdap_execution_setup

        n_tpl_lib = n_elements(template_libraries)      ; Number of template libraries to use
        n_ems_par = n_elements(emission_line_parameters); Number of emission line parameter files
        n_abs_par = n_elements(absorption_line_parameters);Number of absorption line parameter files

        success = MDAP_CHECK_FILE_EXISTS(template_libraries, /search)   ; Aborts of fails
        if success ne 0 then $
            message, 'Template libraries not correctly defined.'
                
        success = MDAP_CHECK_FILE_EXISTS(emission_line_parameters)
        if success ne 0 then $
            message, 'Emission-line parameter files not correctly defined.'

        success = MDAP_CHECK_FILE_EXISTS(absorption_line_parameters)
        if success ne 0 then $
            message, 'Absorption-line parameter files not correctly defined.'

        MDAP_BUILD_EXECUTION_PLANS, n_tpl_lib, n_ems_par, n_abs_par, bin_par, ell, pa, Reff, $
                                    w_range_sn, threshold_ston_bin, w_range_analysis, $
                                    threshold_ston_analysis, analysis, analysis_par, $
                                    analysis_prior, tpl_lib_analysis, ems_par_analysis, $
                                    abs_par_analysis, overwrite_flag, output_file_root, $
                                    execution_plan

        ; TODO: If any of the execution_plans are supposed to overwrite
        ; existing data, check if the files exist and remove them?

        ;MDAP_GENERATE_OUTPUT_FILE_NAMES, output_file_root, execution_plan

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

        if ~keyword_set(quiet) then $
            print, 'BLOCK 0 ... DONE.'
        ; END BLOCK 0 ----------------------------------------------------------



        ; BLOCK 1 --------------------------------------------------------------
        ;       Read the DRP fits file
        ;-----------------------------------------------------------------------

        if ~keyword_set(quiet) then $
            print, 'BLOCK 1 ... '

;       if mdap_read_datacube_version gt mdap_read_datacube_version_previous $
;          or execute_all_modules eq 1 then begin

        ; TODO: Version required for this should depend on the version of DRP
        ; output being read

        ; TODO: Put in the conditions for when this block is executed:
        ;       - If any of the plans have the overwrite flag flipped
        ;       - If any of the plans do not have a previously created output
        ;         file with the binned spectra

        ; TODO: These operations are independent of the plan -------------------
            ; TODO: Allow to read a single spectrum and skip binning step
            MDAP_READ_DRP_FITS, datacube_name, header, flux, ivar, mask, wave, sres, skyx, skyy, $
                                type=mode, unit=unit
            mask[*,*] = 0.                      ; TODO: Unmask everything for now
            ndrp = (size(flux))[1]              ; Number of DRP spectra

            ; Get the spaxel size
            MDAP_GET_SPAXEL_SIZE, header, spaxel_dx, spaxel_dy, type=mode, unit=unit
;           print, spaxel_dx, spaxel_dy

            ; Set a fiducial set of coordinates to use for each spectrum
            MDAP_FIDUCIAL_BIN_XY, wave, skyx, skyy, bskyx, bskyy, type=mode

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
            ;   1. Wavelength in angstroms
            ;   2. Resolution (lamda/delta lambda)
            ;   3. delta lamga (FWHM) in angstroms
            ;   4. delta lamga (FWHM) in km/s
            if n_elements(instrumental_fwhm_file) ne 0 then begin
                if file_test(instrumental_fwhm_file) eq 0 then begin
                    print, 'Unable to read '+instrumental_fwhm_file+'!  Using fits extention or ' $
                           + 'default.'
                endif
                print, 'Reading '+instrumental_fwhm_file
                READCOL, instrumental_fwhm_file, ww, r, fwhm_ang, fwhm_kms, /silent
            endif
            
            ; Initialize the instrumental resolution.
            if n_elements(fwhm_ang) ne 0 then begin             ; use data from input file
                if ~keyword_set(nolog) then begin
                    printf, log_file_unit, '[INFO] Spectral resolution changed using data in: '+ $
                            instrumental_fwhm_file
                endif
                sres=interpol(fwhm_ang,r,wave)          ; Interpolate to the object wavelengths
            endif else if n_elements(sres) eq 0 then begin      ; sres not available from fits file
                if ~keyword_set(nolog) then begin
                    printf, log_file_unit, '[INFO] Spectral resolution unavailable, assume R=2000'
                endif
                sres = make_array(n_elements(wave), /double, value=2000.0d)     ; Set R=2000.
            endif
            ;-------------------------------------------------------------------

;           sz=size(mask)
;           print, sz
;           print, mask[0,0:100]

        ;-----------------------------------------------------------------------
        if ~keyword_set(quiet) then $
            print, 'BLOCK 1 ... DONE.'
        ; END BLOCK 1 ----------------------------------------------------------

        ; Once the DRP fits file has been read, the required input to
        ; appropriatly manipulate ALL the template libraries is
        ; available

        ; BLOCK 2 --------------------------------------------------------------
        ;       - Read the template library
        ;       - Match the resolution of the templates to that of the galaxy
        ;         spectra (as best as possible).
        ;       - Resample the templates to match the object sampling
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
        indx=where(use_tpl_lib eq 1, count)
;       if indx[0] ne -1 then begin
        if count ne 0 then begin

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

                endfor  ; Loop over execution plans
            endfor      ; Loop over template libraries
        endif           ; If any template libraries are used

        if ~keyword_set(quiet) then $
            print, 'BLOCK 2 ... DONE.'
        ; END BLOCK 2 ----------------------------------------------------------

        ;-----------------------------------------------------------------------
        ; ANALYSIS PREP --------------------------------------------------------

        ; Initialize the extinction curve
        MW_extinction = 0.                      ; TODO: replaced with appropriate extension in the
                                                ;       input file, or total_filelist

        ; From Oliver Steele (04 Sep 2014)
        ;       Requires galRA, galDEC, equinox, lambda_obs, spectrumdata
        ;GLACTC,galRA,galDEC,equinox,gl,gb,1,/DEGREE    ; Converts RA/DEC to galactic coords 
        ;mw_ebmv = DUST_GETVAL(gl,gb,/interp)           ; Grabs E(B-V) from Schlegel dust maps
        ;FM_UNRED, lambda_obs, spectrumdata, mw_ebmv    ; De-reddens using ?? extinction curve

;       print, star_kin_starting_guesses                ; h3 and h4 initialized to 0
;       print, gas_kin_starting_guesses

        ;-----------------------------------------------------------------------
        ; Now can begin cycling through execution plans
        if keyword_set(dbg) then $
            n_plans=1                           ; Only run the first plan
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
                MDAP_ERASEVAR, eml_par  ; Make sure it doesn't exist if no ems_par is defined

;           print, n_elements(eml_par)
;           stop

            ;-------------------------------------------------------------------
            ; Read the spectral-index parameters; this is done at every
            ; iteration because it's fairly quick.  TODO: Check that this I/O
            ; doesn't slow things down too much!
            if execution_plan[i].abs_par ne -1 then begin
                abs_par = MDAP_READ_ABSORPTION_LINE_PARAMETERS( absorption_line_parameters[ $
                                                                 execution_plan[i].abs_par ])
                nabs = n_elements(abs_par)              ; Number of absorption-line indices
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
                MDAP_ERASEVAR, abs_par  ; Make sure it doesn't exist if no abs_par is defined

            ; Delete the exiting file if overwrite=1
            if execution_plan[i].overwrite eq 1 and FILE_TEST(execution_plan[i].ofile) eq 1 then $
            begin
                print, 'Removing existing file: '+execution_plan[i].ofile
                FILE_DELETE, execution_plan[i].ofile
            endif

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
            if perform_block.eml_only eq 1 then begin
                print, '    Emission-line-only fit: YES'
            endif else $
                print, '    Emission-line-only fit: NO'
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
                ;       MDAP_CALCULATE_SN and MDAP_SPATIAL_BINNING

                MDAP_SELECT_GOOD_SPECTRA, flux, ivar, mask, gflag, gindx, count=gcount

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
                ;       calculated!  Should they be reread here?

                ; Always need the S/N data

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

            ; Interpolate data from the prior if necessary
            if (perform_block.bin eq 1 or perform_block.spec_fit eq 1) and $
               strlen(execution_plan[i].analysis_prior) ne 0 then begin
               
               print, 'Interpolating kinematics for later use'
               MDAP_INTERPOLATE_KINEMATICS, execution_plan[i].analysis_prior, bskyx, bskyy, $
                                            star_kin_interp, /velocity, /sigma, /stellar
               MDAP_INTERPOLATE_KINEMATICS, execution_plan[i].analysis_prior, bskyx, bskyy, $
                                            gas_kin_interp, /velocity, /sigma
            endif


            ; BLOCK 4 ----------------------------------------------------------
            ; Spatial binning
            ;-------------------------------------------------------------------
            if ~keyword_set(quiet) then $
                print, 'PLAN '+MDAP_STC(i+1,/integer)+': BLOCK 4 ...'

            if perform_block.bin eq 1 then begin
           
                reg_flux = flux
                reg_ivar = ivar
                reg_mask = mask

                ; De-redshift the spectra before binning the spectra.
                if strlen(execution_plan[i].analysis_prior) ne 0 $
                   and execution_plan[i].bin_par.v_register eq 1 then begin

                    ; Kinematics interpolated above
                    bin_vreg = star_kin_interp[*,0]             ; Register to stellar velocities

                    ; Velocity register the spectra to the median of the input velocity
                    MDAP_VELOCITY_REGISTER, wave, reg_flux, reg_ivar, reg_mask, bin_vreg, $
                                            gflag=gflag, /register_to_median

                    ; Offset interpolations for later use
                    star_kin_interp[*,0] = star_kin_interp[*,0] - bin_vreg
                    gas_kin_interp[*,0] = gas_kin_interp[*,0] - bin_vreg

;                   for j=0,9 do $
;                       print, bin_vreg[i]
                endif else $
                    bin_vreg = dblarr(ndrp)             ; Set to 0.0

                MDAP_SPATIAL_BINNING, reg_flux, reg_ivar, reg_mask, signal, noise, gflag, bskyx, $
                                      bskyy, spaxel_dx, spaxel_dy, execution_plan[i].bin_par, $
                                      execution_plan[i].threshold_ston_bin, bin_wgts, bin_indx, $
                                      bin_flux, bin_ivar, bin_mask, xbin, ybin, bin_rad, bin_area, $
                                      bin_ston, nbin, sn_calibration=sn_calibration, plot=plot

                ; Write the version of the spatial binning to the header
                SXADDPAR, header, 'VDAPBIN', mdap_spatial_binning_version, $
                          'mdap_spatial_binning version'

                ; Add some information to the log
                if ~keyword_set(nolog) then begin
                    printf, log_file_unit, '[INFO] Type of spatial binning: ', $
                            execution_plan[i].bin_par.type
                    printf, log_file_unit, '[INFO] Spatial bining version: ' + $
                            mdap_spatial_binning_version
                    if execution_plan[i].bin_par.type ne 'NONE' then begin
                        if n_elements(execution_plan[i].bin_par) then begin
                            printf, log_file_unit, '[INFO] Bin parameter: ', $
                                    execution_plan[i].bin_par
                        endif
                        printf, log_file_unit, '[INFO] Bin threshold S/N: ', $
                                execution_plan[i].threshold_ston_bin
                        printf, log_file_unit, '[INFO] Optimal weighting by S/(N)^2: ', $
                                execution_plan[i].optimal_weighting
                    endif
                    printf, log_file_unit, '[INFO] Number of bins: ', n_elements(nbin)
                endif

                ; Write the binning results
                MDAP_WRITE_OUTPUT, execution_plan[i].ofile, header=header, $
                                   bin_par=execution_plan[i].bin_par, $
                                   threshold_ston_bin=execution_plan[i].threshold_ston_bin, $
                                   bin_vreg=bin_vreg, bin_indx=bin_indx, bin_weights=bin_wgts, $
                                   wave=wave, sres=sres, bin_flux=bin_flux, bin_ivar=bin_ivar, $
                                   bin_mask=bin_mask, xbin_rlow=xbin, ybin_rupp=ybin, $
                                   rbin=bin_rad, bin_area=bin_area, bin_ston=bin_ston, $
                                   bin_n=nbin, /read_header, quiet=quiet
            endif else begin

                ; Always need the binned data

                print, 'READING EXISTING BIN DATA'

                ; Read the binning results
                MDAP_DEFINE_OUTPUT, header=header, bin_vreg=bin_vreg, bin_indx=bin_indx, $
                                    bin_weights=bin_wgts, bin_flux=bin_flux, bin_ivar=bin_ivar, $
                                    bin_mask=bin_mask, xbin_rlow=xbin, ybin_rupp=ybin, $
                                    rbin=bin_rad, bin_area=bin_area, bin_ston=bin_ston, bin_n=nbin
                MDAP_READ_OUTPUT, execution_plan[i].ofile, header=header, bin_vreg=bin_vreg, $
                                  bin_indx=bin_indx, bin_weights=bin_wgts, bin_flux=bin_flux, $
                                  bin_ivar=bin_ivar, bin_mask=bin_mask, xbin_rlow=xbin, $
                                  ybin_rupp=ybin, rbin=bin_rad, bin_area=bin_area, $
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
            indx = where(execution_plan[i].analysis eq 1, count)
;           if indx[0] eq -1 then $             ; No analysis to perform
            if count eq 0 then $             ; No analysis to perform
                continue

            ;-------------------------------------------------------------------
            ; Read the pre-existing template library
            ;   The template library is required for ANY analysis
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

                ; TODO: Allow for starting guesses from a
                ; cross-correlation with a guess template?

                ; TODO: Stellar kinematic guess are always 4 (as
                ; required by MDAP_SPECTRAL_FITTING) as opposed to be
                ; set by execution_plan.analysis_par.moments.

                MDAP_INITIALIZE_GUESS_KINEMATICS, n_elements(nbin), $
                                                  execution_plan[i].analysis_prior, $
                                                  star_kin_interp, gas_kin_interp, bin_indx, $
                                                  velocity_initial_guess, $
                                                  velocity_dispersion_initial_guess, $
                                                  star_kin_guesses, gas_kin_guesses

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
                                   tpl_library_key=tpl_library_keys[execution_plan[i].tpl_lib], $
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
                
                    ; Get the instrumental dispersion at the fitted line centers
                    MDAP_INSTR_DISPERSION_AT_EMISSION_LINE, wave, sres, eml_par, $
                                                            emission_line_omitted, $
                                                            emission_line_kinematics_individual, $
                                                            emission_line_sinst
                
                    MDAP_WRITE_OUTPUT, execution_plan[i].ofile, header=header, $
                                       ems_line_key=ems_line_keys[execution_plan[i].ems_par], $
                                       eml_par=eml_par, obj_fit_mask_gndf=obj_fit_mask_gndf, $
                                       weights_gndf=weights_gndf, $
                                       mult_poly_coeff_gndf=mult_poly_coeff_gndf, $
                                       emission_line_kinematics_avg=emission_line_kinematics, $
                                       emission_line_kinematics_aer=emission_line_kinematics_err, $
                                       chi2_gndf=chi2_gndf, $
                                emission_line_kinematics_ind=emission_line_kinematics_individual, $
                        emission_line_kinematics_ier=emission_line_kinematics_individual_err, $
                                       emission_line_sinst=emission_line_sinst, $
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
;               MDAP_WRITE_OUTPUT, execution_plan[i].ofile, optimal_template=best_template, $
;                                  losvd_optimal_template=best_template_losvd_conv, quiet=quiet
            endif else begin

                ; TODO: This only needs to be read if the spectral index
                ; measurements will be measured

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
                                    emission_line_sinst=emission_line_sinst, $
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
                                  emission_line_sinst=emission_line_sinst, $
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

            ; END BLOCK 5 ----------------------------------------------
            ;-----------------------------------------------------------

            if MDAP_MORE_ANALYSIS_BLOCKS_TO_FINISH(perform_block, /spec_fit) eq 0 then $
                continue

            ; BLOCK 6 --------------------------------------------------
            ;   Fit the emission-line only spectrum
            ;-----------------------------------------------------------
            ; Check if the emission-line-only spectra should be fit
            if perform_block.eml_only eq 1 then begin

                if ~keyword_set(quiet) then $
                    print, 'PLAN '+MDAP_STC(i+1,/integer)+': BLOCK 6 ...'

                ;-------------------------------------------------------
                ; Get the bestfit continuum spectrum and the kinematics
                ; to use
                nb = n_elements(nbin)
                sz = size(bin_flux)
                ; Default is to set the continuum to 0
                bestfit_continuum = dblarr(sz[1], sz[2])

                ; Use pPXF continuum if available
                if (size(bestfit_ppxf))[1] eq sz[1] and (size(bestfit_ppxf))[2] eq sz[2] then $
                    bestfit_continuum = bestfit_ppxf

                ; Change to the GANDALF result if it exists
                if (size(bestfit_gndf))[1] eq sz[1] and (size(bestfit_gndf))[2] eq sz[2] then $
                    bestfit_continuum = bestfit_gndf - eml_model

                ; Get the input stellar kinematics
                if (size(stellar_kinematics))[1] eq (size(tpl_flux))[1] then begin
                    inp_kin = stellar_kinematics[*,0:1]
                endif else begin
                    MDAP_INITIALIZE_GUESS_KINEMATICS, n_elements(nbin), $
                                                      execution_plan[i].analysis_prior, $
                                                      star_kin_interp, gas_kin_interp, bin_indx, $
                                                      velocity_initial_guess, $
                                                      velocity_dispersion_initial_guess, $
                                                      star_kin_guesses, gas_kin_guesses

                    inp_kin = star_kin_guesses[*,0:1]
                endelse

                ; Perform the fit using Enci's code
                MDAP_EMISSION_LINE_ONLY_FIT, wave, bin_flux, bin_ivar, bin_mask, $
                                             bestfit_continuum, inp_kin, elo_ew_eml_model, $
                                             elo_ew_kinematics, elo_ew_kinematics_err, $
                                             elo_ew_omitted, elo_ew_kinematics_ind, $
                                             elo_ew_kinematics_ier, elo_ew_intens, elo_ew_interr, $
                                             elo_ew_fluxes, elo_ew_flxerr, elo_ew_EWidth, $
                                             elo_ew_EW_err, eml_par=emlo_par, quiet=quiet, $
                                             dbg=dbg, /enci

                ; Get the instrumental dispersion at the fitted line center
                MDAP_INSTR_DISPERSION_AT_EMISSION_LINE, wave, sres, emlo_par, elo_ew_omitted, $
                                                        elo_ew_kinematics_ind, elo_ew_sinst
                
                ; Perform the fit using Francesco's code
                MDAP_EMISSION_LINE_ONLY_FIT, wave, bin_flux, bin_ivar, bin_mask, $
                                             bestfit_continuum, inp_kin, elo_fb_eml_model, $
                                             elo_fb_kinematics, elo_fb_kinematics_err, $
                                             elo_fb_omitted, elo_fb_kinematics_ind, $
                                             elo_fb_kinematics_ier, elo_fb_intens, elo_fb_interr, $
                                             elo_fb_fluxes, elo_fb_flxerr, elo_fb_EWidth, $
                                             elo_fb_EW_err, quiet=quiet, dbg=dbg, /belfiore

                ; Get the instrumental dispersion at the fitted line center
                MDAP_INSTR_DISPERSION_AT_EMISSION_LINE, wave, sres, emlo_par, elo_fb_omitted, $
                                                        elo_fb_kinematics_ind, elo_fb_sinst
                
                ; Write the data
                MDAP_WRITE_OUTPUT, execution_plan[i].ofile, emlo_par=emlo_par, $
                                   elo_ew_kinematics_avg=elo_ew_kinematics, $
                                   elo_ew_kinematics_aer=elo_ew_kinematics_err, $
                                   elo_ew_kinematics_ind=elo_ew_kinematics_ind, $
                                   elo_ew_kinematics_ier=elo_ew_kinematics_ier, $
                                   elo_ew_sinst=elo_ew_sinst, elo_ew_omitted=elo_ew_omitted, $
                                   elo_ew_intens=elo_ew_intens, elo_ew_interr=elo_ew_interr, $
                                   elo_ew_fluxes=elo_ew_fluxes, elo_ew_flxerr=elo_ew_flxerr, $
                                   elo_ew_EWidth=elo_ew_EWidth, elo_ew_EW_err=elo_ew_EW_err, $
                                   elo_ew_eml_model=elo_ew_eml_model, $
                                   elo_fb_kinematics_avg=elo_fb_kinematics, $
                                   elo_fb_kinematics_aer=elo_fb_kinematics_err, $
                                   elo_fb_kinematics_ind=elo_fb_kinematics_ind, $
                                   elo_fb_kinematics_ier=elo_fb_kinematics_ier, $
                                   elo_fb_sinst=elo_fb_sinst, elo_fb_omitted=elo_fb_omitted, $
                                   elo_fb_intens=elo_fb_intens, elo_fb_interr=elo_fb_interr, $
                                   elo_fb_fluxes=elo_fb_fluxes, elo_fb_flxerr=elo_fb_flxerr, $
                                   elo_fb_EWidth=elo_fb_EWidth, elo_fb_EW_err=elo_fb_EW_err, $
                                   elo_fb_eml_model=elo_fb_eml_model

            endif else begin

                ; If the data is already available, the only thing that
                ; needs to be read is the fitted model.  For now, the
                ; preference is for Enci's results.

                MDAP_DEFINE_OUTPUT, elo_ew_eml_model=elo_ew_eml_model
                MDAP_READ_OUTPUT, execution_plan[i].ofile, elo_ew_eml_model=elo_ew_eml_model

            endelse

            if ~keyword_set(quiet) then $
                print, 'PLAN '+MDAP_STC(i+1,/integer)+': BLOCK 6 ... DONE.'

            ; END BLOCK 6 ------------------------------------------------------
            ;-------------------------------------------------------------------

            if MDAP_MORE_ANALYSIS_BLOCKS_TO_FINISH(perform_block, /eml_only) eq 0 then $
                continue

            ; BLOCK 7 ----------------------------------------------------------
            ; Perform the spectral-index measurements
            ;-------------------------------------------------------------------

            ; Check if spectral indices should be measured
            if perform_block.spec_indx eq 1 then begin

                if ~keyword_set(quiet) then $
                    print, 'PLAN '+MDAP_STC(i+1,/integer)+': BLOCK 7 ...'

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
                    if n_elements(elo_ew_eml_model) ne 0 then begin
                        eml_model = elo_ew_eml_model        ; Use Enci's model
                    endif else $
                        eml_model = dblarr(sz[1], sz[2])    ; No emission-line model
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
                ; BLOCK 7a; should this become it's own block?

                ; TODO: is si_bin_wave ever different from wave?

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
                    print, 'PLAN '+MDAP_STC(i+1,/integer)+': BLOCK 7 ... DONE.'

            endif


            ; END computation of "model-independent" data products #############
            ; ##################################################################


        endfor
        ; END loop over plans --------------------------------------------------
        ;-----------------------------------------------------------------------

        ; close up
        if ~keyword_set(nolog) then $
            MDAP_LOG_CLOSE, log_file_unit

        print, 'DAP FINISHED SUCCESSFULLY'
END

;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------


