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
; PROCEDURES CALLED:
;       MDAP_INITIALIZATION_BLOCK
;       MDAP_DRP_DATA_BLOCK
;       SXADDPAR
;       MDAP_TEMPLATE_LIBRARY_BLOCK
;       MDAP_ANALYSIS_SETUP_BLOCK
;       MDAP_WRITE_OUTPUT
;       MDAP_DRP_SNR_BLOCK
;       MDAP_BINNING_BLOCK
;       MDAP_SET_TPL_LIB_OUTPUT_FILE
;       MDAP_READ_RESAMPLED_TEMPLATES
;       MDAP_FULL_SPECTRAL_FIT_BLOCK
;       MDAP_EMISSION_LINE_FIT_BLOCK
;       MDAP_SPECTRAL_INDEX_BLOCK
;       MDAP_LOG_CLOSE
;
; INTERNAL SUPPORT ROUTINES:
;       @mdap_analysis_block_logic
;
; REVISION HISTORY:
;       01 Sep 2014: Copied from v0_8 version by L. Coccato.
;       02 Sep 2014: Formating and minor edits by K. Westfall (KBW)
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
;       22 Jan 2015: (KBW) Add limits for the guess kinematics
;       30 Jan 2015: (KBW) Pulled all the internal support routines, and
;                          pulled out blocks into their own routines as
;                          suggested by M. Cappellari.
;       16 Mar 2015: (KBW) Remove sn_calibration flag for
;                          MDAP_BINNING_BLOCK; noise_calib now included
;                          in bin_par
;       17 Mar 2015: (KBW) Removed remove_null_templates and
;                          save_intermediate_steps from initialization
;                          block.  Allowed for plan parameter file
;                          input.
;-
;-----------------------------------------------------------------------

;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
; This is the file (mdap_analysis_block_logic.pro) that has all the
; functions that decide which analysis blocks should be performed
@mdap_analysis_block_logic
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------

PRO MANGA_DAP, $
        inptbl=inptbl, index=index, par=par, plan=plan, drppath=drppath, dappath=dappath, $
        dapsrc=dapsrc, nolog=nolog, quiet=quiet, plot=plot, dbg=dbg

        manga_dap_version = MDAP_DEFINE_MANGA_DAP_VERSION()
        manga_dap_version.main = '0.96'     ; set the version number
        t0=systime(/seconds)/60./60.        ; start time
        thismem = memory()                  ; start memory
        maxmem = 0                          ; maximum memory

        ; BLOCK 0 ------------------------------------------------------
        ; Initialization
        ;       - File I/O setup
        ;       - Execution plan setup
        ;       - Version control
        ;---------------------------------------------------------------
        if ~keyword_set(quiet) then $
            print, 'BLOCK 0 ... '

        MDAP_INITIALIZATION_BLOCK, manga_dap_version, tpl_library_keys, template_libraries, $
                                   tpl_vacuum_wave, ems_line_keys, emission_line_parameters, $
                                   ems_vacuum_wave, abs_line_keys, absorption_line_parameters, $
                                   abs_vacuum_wave, execution_plan, plate, $
                                   ifudesign, mode, velocity_initial_guess, $
                                   velocity_dispersion_initial_guess, ell, pa, Reff, $
                                   datacube_name, output_file_root, n_tpl_lib, inptbl=inptbl, $
                                   index=index, par=par, plan=plan, drppath=drppath, $
                                   dappath=dappath, dapsrc=dapsrc, $
                                   bvls_shared_lib=bvls_shared_lib, nolog=nolog, $
                                   log_file_unit=log_file_unit, quiet=quiet
;                                  save_intermediate_steps=save_intermediate_steps, $
;                                  remove_null_templates=remove_null_templates, $
;       get_lun, unit
;       print, unit
;       free_lun, unit

        n_plans = n_elements(execution_plan)

        if ~keyword_set(quiet) then $
            print, 'BLOCK 0 ... DONE.'
        ; END BLOCK 0 --------------------------------------------------


        ; BLOCK 1 ------------------------------------------------------
        ;       Read the DRP fits file
        ;---------------------------------------------------------------

        if ~keyword_set(quiet) then $
            print, 'BLOCK 1 ... '

        MDAP_DRP_DATA_BLOCK, manga_dap_version, datacube_name, header, wave, flux, ivar, mask, $
                             sres, spaxel_dx, spaxel_dy, bskyx, bskyy, type=mode, $
                             nolog=nolog, log_file_unit=log_file_unit
;                            instrumental_fwhm_file=instrumental_fwhm_file, $
        
;       get_lun, unit
;       print, unit
;       free_lun, unit

        ndrp = (size(flux))[1]              ; Number of DRP spectra

        if ~keyword_set(quiet) then $
            print, 'BLOCK 1 ... DONE.'
        ; END BLOCK 1 --------------------------------------------------


        ;---------------------------------------------------------------
        ; Add the version to the header; can't be done in
        ; MDAP_INITIALIZATION_BLOCK because the header has to be read
        ; first.

        SXADDPAR, header, 'VDAPPLAN', manga_dap_version.execute_plan, $
                          'mdap_build_execution_plans version'


        ; BLOCK 2 ------------------------------------------------------
        ; With the DRP data read, the required input to appropriately
        ; manipulate ALL the template libraries is available.  This
        ; block:
        ;   - Read the template library
        ;   - Match the resolution of the templates to that of the
        ;     galaxy spectra (as best as possible).
        ;   - Resample the templates to match the object sampling
        ;---------------------------------------------------------------
        if ~keyword_set(quiet) then $
            print, 'BLOCK 2 ... '

        MDAP_TEMPLATE_LIBRARY_BLOCK, manga_dap_version, output_file_root, n_tpl_lib, $
                                     tpl_library_keys, template_libraries, tpl_vacuum_wave, $
                                     abs_line_keys, velocity_initial_guess, wave, sres, $
                                     execution_plan, velScale, nolog=nolog, $
                                     log_file_unit=log_file_unit, datacube_name=datacube_name, $
                                     quiet=quiet

        if ~keyword_set(quiet) then $
            print, 'BLOCK 2 ... DONE.'
        ; END BLOCK 2 --------------------------------------------------

        ;---------------------------------------------------------------
        ; TODO: Initialize the extinction curve
        MW_extinction = 0.

        ; From Oliver Steele (04 Sep 2014)
        ;       Requires galRA, galDEC, equinox, lambda_obs, spectrumdata
        ;GLACTC,galRA,galDEC,equinox,gl,gb,1,/DEGREE    ; Converts RA/DEC to galactic coords 
        ;mw_ebmv = DUST_GETVAL(gl,gb,/interp)           ; Grabs E(B-V) from Schlegel dust maps
        ;FM_UNRED, lambda_obs, spectrumdata, mw_ebmv    ; De-reddens using ?? extinction curve

        ;---------------------------------------------------------------
        ; Now can begin cycling through execution plans
        if keyword_set(dbg) then $
            n_plans=1                           ; Only run the first plan
        for i=0,n_plans-1 do begin

            if ~keyword_set(quiet) then $
                print, 'Beginning ExecutionPlans: ', i+1

            ;-----------------------------------------------------------
            ; No reading or writing occurs in this block
            ; Block:
            ;   - Reads the template library, emission-line parameters,
            ;     and spectral-index parameters to use
            ;   - Determines which of the remaining blocks to perform
            ;   - Interpolates guess kinematics based on the provided
            ;     analysis_prior in the execution_plan (if necessary)
            MDAP_ANALYSIS_SETUP_BLOCK, output_file_root, tpl_library_keys, $
                                       emission_line_parameters, ems_vacuum_wave, $
                                       absorption_line_parameters, abs_vacuum_wave, ndrp, $
                                       execution_plan[i], bskyy, bskyx, eml_par, abs_par, $
                                       perform_block, star_kin_interp, gas_kin_interp

            get_lun, unit
            print, unit
            free_lun, unit
            ;stop

            ;-----------------------------------------------------------
            ; Check there are blocks to be perform
            if MDAP_ALL_ANALYSIS_BLOCKS_COMPLETED(perform_block) eq 1 then begin
                if ~keyword_set(quiet) then $
                    print, 'All blocks previously completed for ExecutionPlan: ', i+1
                continue
            endif

            ;-----------------------------------------------------------
            ; Re-write part of the DRPS extension; no need to re-read
            ; them because they're always calculated above
            if execution_plan[i].overwrite eq 1 then begin
                ; Write basics of the DRP fits files
                MDAP_WRITE_OUTPUT, execution_plan[i].ofile, header=header, dx=spaxel_dx, $
                                   dy=spaxel_dy, xpos=bskyx, ypos=bskyy, quiet=quiet
            endif
            
            ; BLOCK 3 --------------------------------------------------
            ; S/N calculation.  If the calculation is performed, the
            ; procedure will write signal and noise to DRPS extension;
            ; if the block is not performed, the signal and noise are
            ; read from the existing file (execution_plan[i].ofile)
            ;-----------------------------------------------------------
            if ~keyword_set(quiet) then $
                print, 'PLAN '+MDAP_STC(i+1,/integer)+': BLOCK 3 ...'

            MDAP_DRP_SNR_BLOCK, manga_dap_version, execution_plan[i], perform_block, header, wave, $
                                flux, ivar, mask, gflag, signal, noise, nolog=nolog, $
                                log_file_unit=log_file_unit, quiet=quiet
            get_lun, unit
            print, unit
            free_lun, unit
            ;stop

            if ~keyword_set(quiet) then $
                print, 'PLAN '+MDAP_STC(i+1,/integer)+': BLOCK 3 ... DONE.'
            ; END BLOCK 3 ----------------------------------------------

            ; Check if there is a need to continue to the next block
            if MDAP_MORE_ANALYSIS_BLOCKS_TO_FINISH(perform_block, /ston) eq 0 then $
                continue

            ; BLOCK 4 --------------------------------------------------
            ; Spatial binning
            ;-----------------------------------------------------------
            if ~keyword_set(quiet) then $
                print, 'PLAN '+MDAP_STC(i+1,/integer)+': BLOCK 4 ...'

            MDAP_BINNING_BLOCK, manga_dap_version, execution_plan[i], perform_block, header, $
                                spaxel_dx, spaxel_dy, wave, sres, flux, ivar, mask, bskyx, bskyy, $
                                gflag, signal, noise, velocity_initial_guess, star_kin_interp, $
                                gas_kin_interp, bin_vreg, reg_velocity_initial_guess, bin_wgts, $
                                bin_indx, bin_flux, bin_ivar, bin_mask, xbin, ybin, bin_rad, $
                                bin_area, bin_ston, nbin, plot=plot, nolog=nolog, $
                                log_file_unit=log_file_unit, quiet=quiet

            get_lun, unit
            print, unit
            free_lun, unit
            ;stop

            if ~keyword_set(quiet) then $
                print, 'PLAN '+MDAP_STC(i+1,/integer)+': BLOCK 4 ... DONE.'
            ; END BLOCK 4 ----------------------------------------------
            ;-----------------------------------------------------------

;           print, 'unmasked pixels:', n_elements(where(bin_mask lt 1.))

            ; Finish if the user only wants the binned spectra
            ; Check if there is a need to continue to the next block
            if MDAP_MORE_ANALYSIS_BLOCKS_TO_FINISH(perform_block, /bin) eq 0 then $
                continue


            ;***********************************************************
            ;***********************************************************
            ;***********************************************************
            ; BEGIN computation of "model-independent" data products 


            ;-----------------------------------------------------------
            ; Read the pre-existing template library.  The template
            ; library is requiredd for all analysis except for the
            ; emission-line-only fits, if the user doesn't care about
            ; the background.  WARNING: Should not be able to reach here
            ; without first creating tpl_out_fits!

            ; TODO: What happens if only the emission-line-only analysis
            ; is selected?
            tpl_out_fits = MDAP_SET_TPL_LIB_OUTPUT_FILE(output_file_root, $
                                                        tpl_library_keys[execution_plan[i].tpl_lib])
           
            MDAP_READ_RESAMPLED_TEMPLATES, tpl_out_fits, tpl_wave, tpl_flux, tpl_ivar, tpl_mask, $
                                           tpl_sres, tpl_soff

            get_lun, unit
            print, unit
            free_lun, unit
            ;stop

            ; TODO: Need to account for tpl_soff further down the line?
            ; Or not because HARD-WIRED default is tpl_soff=0

            ; BLOCK 5 --------------------------------------------------
            ; Run a full spectral fit including stellar population, gas
            ; and stellar kinematics
            ;-----------------------------------------------------------
            if ~keyword_set(quiet) then $
                print, 'PLAN '+MDAP_STC(i+1,/integer)+': BLOCK 5 ...'

            MDAP_FULL_SPECTRAL_FIT_BLOCK, manga_dap_version, execution_plan[i], perform_block, $
                                          star_kin_interp, gas_kin_interp, bin_indx, $
                                          reg_velocity_initial_guess, $
                                          velocity_dispersion_initial_guess, tpl_library_keys, $
                                          ems_line_keys, header, wave, sres, nbin, bin_flux, $
                                          bin_ivar, bin_mask, tpl_wave, tpl_flux, tpl_ivar, $
                                          tpl_mask, obj_fit_mask_ppxf, weights_ppxf, bestfit_ppxf, $
                                          obj_fit_mask_gndf, weights_gndf, bestfit_gndf, $
                                          eml_model, stellar_kinematics, eml_par=eml_par, $
                                          bvls_shared_lib=bvls_shared_lib, quiet=quiet, plot=plot, $
                                          dbg=dbg

            get_lun, unit
            print, unit
            free_lun, unit
            ;stop

            if ~keyword_set(quiet) then $
                print, 'PLAN '+MDAP_STC(i+1,/integer)+': BLOCK 5 ... DONE.'

            ; END BLOCK 5 ----------------------------------------------
            ;-----------------------------------------------------------

            ; Check if there is a need to continue to the next block
            if MDAP_MORE_ANALYSIS_BLOCKS_TO_FINISH(perform_block, /spec_fit) eq 0 then $
                continue

            ; BLOCK 6 --------------------------------------------------
            ;   Fit the emission-line only spectrum
            ;-----------------------------------------------------------
            if ~keyword_set(quiet) then $
                print, 'PLAN '+MDAP_STC(i+1,/integer)+': BLOCK 6 ...'

            MDAP_EMISSION_LINE_FIT_BLOCK, manga_dap_version, execution_plan[i], perform_block, $
                                          star_kin_interp, gas_kin_interp, $
                                          reg_velocity_initial_guess, $
                                          velocity_dispersion_initial_guess, wave, sres, nbin, $
                                          bin_indx, bin_flux, bin_ivar, bin_mask, bestfit_ppxf, $
                                          bestfit_gndf, eml_model, stellar_kinematics, $
                                          elo_ew_eml_model, quiet=quiet, dbg=dbg

            get_lun, unit
            print, unit
            free_lun, unit
            ;stop

            if ~keyword_set(quiet) then $
                print, 'PLAN '+MDAP_STC(i+1,/integer)+': BLOCK 6 ... DONE.'

            ; END BLOCK 6 ------------------------------------------------------
            ;-------------------------------------------------------------------

            ; Check if there is a need to continue to the next block
            if MDAP_MORE_ANALYSIS_BLOCKS_TO_FINISH(perform_block, /eml_only) eq 0 then $
                continue

            ; BLOCK 7 ----------------------------------------------------------
            ;   Perform the spectral-index measurements
            ;-------------------------------------------------------------------

            if ~keyword_set(quiet) then $
                print, 'PLAN '+MDAP_STC(i+1,/integer)+': BLOCK 7 ...'

            ; Get the name of the template library file with the
            ; resolution matched to the spectral-index system.
            tpl_out_fits = MDAP_SET_TPL_LIB_OUTPUT_FILE(output_file_root, $
                                                    tpl_library_keys[execution_plan[i].tpl_lib], $
                                            abs_line_key=abs_line_keys[execution_plan[i].abs_par])

            MDAP_SPECTRAL_INDEX_BLOCK, manga_dap_version, execution_plan[i], perform_block, $
                                       (size(tpl_flux))[1], tpl_out_fits, abs_line_keys, $
                                       abs_par, obj_fit_mask_ppxf, weights_ppxf, bestfit_ppxf, $
                                       obj_fit_mask_gndf, weights_gndf, bestfit_gndf, $
                                       eml_model, stellar_kinematics, elo_ew_eml_model, wave, $
                                       sres, bin_flux, bin_ivar, bin_mask, $
                                       remove_outliers=remove_outliers, quiet=quiet, dbg=dbg

            get_lun, unit
            print, unit
            free_lun, unit
            ;stop

            if ~keyword_set(quiet) then $
                print, 'PLAN '+MDAP_STC(i+1,/integer)+': BLOCK 7 ... DONE.'

            ; END computation of "model-independent" data products #############
            ; ##################################################################


        endfor
        ; END loop over plans --------------------------------------------------
        ;-----------------------------------------------------------------------

        ; As in DRP...
        ; Track memory usage
        thismem = memory()
        maxmem = maxmem > thismem[3]
        print, 'Max memory usage = '+string(maxmem/1e6,format='(f7.1)')+' MB'
        ; ... and execution time
        print, 'Total execution time = '+string(systime(/seconds)/3600.-t0,format='(f7.3)')+' hr'

        ; close up
        if ~keyword_set(nolog) then $
            MDAP_LOG_CLOSE, log_file_unit

        print, 'DAP FINISHED SUCCESSFULLY'
END

;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------


