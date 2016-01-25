;+
; NAME:
;       MDAP_SPECTRAL_INDEX_BLOCK
;
; PURPOSE:
;   MAIN BLOCK FOR MANGA_DAP:
;   
;   Given that there is not another main execution block after this one
;   in MANGA_DAP, this procedure only executes in full if the
;   spectral-index analysis is performed.  I.e., nothing is returned to
;   the main MANGA_DAP procedure.
;
;   If the spectral-index measurements are requested:
;   - The input is used to set the weights and stellar kinematics to use
;     for the best fitting model continuum.
;   - The galaxy spectra are resolution matched to the resolution of the
;     spectral-index system, and these resolution-matched spectra are
;     written to the output file.
;   - The spectral index measurements are calculated and written to the
;     output file.
;
; CALLING SEQUENCE:
;       MDAP_SPECTRAL_INDEX_BLOCK, manga_dap_version, execution_plan, perform_block, ntpl, $
;                                  tpl_out_fits, abs_line_keys, abs_par, obj_fit_mask_ppxf, $
;                                  weights_ppxf, bestfit_ppxf, $ obj_fit_mask_gndf, weights_gndf, $
;                                  bestfit_gndf, eml_model, stellar_kinematics, elo_ew_eml_model, $
;                                  wave, sres, bin_flux, bin_ivar, bin_mask, $
;                                  remove_outliers=remove_outliers, /quiet, /dbg
;
; INPUTS:
;       manga_dap_version MaNGADAPVersion
;               Structure used to keep track of various
;               version-controlled procedures in the DAP.
;      
;       execution_plan ExecutionPlan
;               Structure providing plan and parameters needed for the
;               analyses.
;
;       perform_block RequiredAnalysisBlock
;               Structure defining which analysis blocks are to be
;               performed.
;
;       star_kin_interp dblarr[N][4]
;               Interpolated set of stellar kinematics based on a
;               provided analysis prior.
;
;       gas_kin_interp dblarr[N][2]
;               Interpolated set of gas kinematics based on a provided
;               analysis prior.
;
;       velocity_initial_guess double
;               Initial guess velocity for the galaxy.
;
;       velocity_dispersion_initial_guess double
;               Initial guess velocity dispersion for the galaxy.
;
;       ntpl integer
;               Number of template spectra.
;
;       tpl_out_fits
;               The name of the fits file containing the template spectra
;               that have been matched to the resolution of the
;               spectral-index system.
;
;       abs_line_keys strarr[]
;               Set of keywords used to uniquely select spectral-index
;               parameter sets.
;
;       abs_par SpectralIndex[A]
;               Array of SpectralIndex structures used to measure the
;               spectral indices.
;
;       weights_ppxf dblarr[B][T]
;               Template weights for each spectrum obtained by PPXF.
;
;       bestfit_ppxf dblarr[B][T]
;               Best fitting spectrum obtained by PPXF for each of the B
;               spectra.
;
;       weights_gndf dblarr[B][T]
;               Template weights for each spectrum obtained by GANDALF.
;
;       bestfit_gndf dblarr[B][T]
;               Best fitting spectrum obtained by GANDALF for each of the B
;               spectra.
;
;       eml_model dblarr[B][T]
;               Best-fitting emission-line-only model for each of the B spectra
;               obtained by GANDALF.
;
;       stellar_kinematics dblarr[B][M]
;               The best-fit stellar kinematics (M moments) for each of the B
;               fitted input galaxy spectra.
;
;       elo_ew_eml_model dblarr[B][T]
;               Best-fitting emission-line-only model determined by the
;               code provided by Enci Wang.
;
;       wave dblarr[T]
;               Wavelength of each spectral channel T.
;
;       sres dblarr[T]
;               Spectral resolution at each wavelength channel T.
;
;       bin_indx intarr[N]
;               Index of the bin that contains each DRP spectrum.
;
;       bin_flux dblarr[B][T]
;               Flux in each of the B binned spectra at each of the T
;               spectral channels.
;
;       bin_ivar dblarr[B][T]
;               Inverse variance in the binned spectra.
;
;       bin_mask dblarr[B][T]
;               Bad pixel mask for the binned spectra.
;
; OPTIONAL INPUTS:
;       remove_outliers double
;               Replace pixels in the object spectra that deviate more
;               than (remove_outliers) times the error (based on either
;               the RMS difference with respect to the best-fit model or
;               directly from the inverse variance value) with the
;               best-fitting model value before resolution matching to
;               the spectral-index system.  REQUIRES the bestfit array.
;
; OPTIONAL KEYWORDS:
;       /quiet
;               Suppress output to the screen.
;
;       /dbg
;               Only perform the analysis on the first provided bin.
;
; OUTPUT:
;
; OPTIONAL OUTPUT:
;
; COMMENTS:
;   The call to MDAP_SPECTRAL_INDEX_MEASUREMENTS_SPECTRA never uses its
;   fit_mask, replace_masked, or rms options.
;
;   In the short-term, added a mask that gives the value pixels in the
;   (broadened) optimal template spectra.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;       MDAP_SPECTRAL_INDEX_MEASUREMENTS
;       MDAP_DEFINE_OUTPUT
;       MDAP_READ_OUTPUT
;       MDAP_SPECTRAL_INDEX_MEASUREMENTS_SPECTRA
;       MDAP_WRITE_OUTPUT
;
; INTERNAL SUPPORT ROUTINES:
;       @mdap_analysis_block_logic
;
; REVISION HISTORY:
;       01 Feb 2015: Pulled from manga_dap.pro by K. Westfall (KBW)
;       12 Aug 2015: (KBW) Added arguments to allow for the
;                          MDAP_SPECTRAL_INDEX_MEASUREMENTS_SPECTRA
;                          procedure to calculate the masks for the
;                          template spectra.  Removed the masks returned
;                          by pPXF and GANDALF in favor of this template
;                          mask.  Changed arguments to
;                          MDAP_SPECTRAL_INDEX_MEASUREMENTS_SPECTRA,
;                          MDAP_READ_OUTPUT, MDAP_DEFINE_OUTPUT,
;                          MDAP_WRITE_OUTPUT to reflect new additional
;                          output.
;-
;-----------------------------------------------------------------------

;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
; This is the file (mdap_analysis_block_logic.pro) is necessary here for
; use of MDAP_CAN_USE_SINDX_IMAGES
@mdap_analysis_block_logic
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------

PRO MDAP_SPECTRAL_INDEX_BLOCK, $
                manga_dap_version, execution_plan, perform_block, star_kin_interp, gas_kin_interp, $
                velocity_initial_guess, velocity_dispersion_initial_guess, ntpl, tpl_out_fits, $
                abs_line_keys, abs_par, weights_ppxf, bestfit_ppxf, weights_gndf, bestfit_gndf, $
                eml_model, stellar_kinematics, elo_ew_eml_model, wave, sres, bin_indx, $
                bin_flux, bin_ivar, bin_mask, remove_outliers=remove_outliers, quiet=quiet, dbg=dbg

        ; Check if spectral indices should be measured
        if perform_block.spec_indx eq 1 then begin

            ;-----------------------------------------------------------
            ; Version control
            MDAP_SPECTRAL_INDEX_MEASUREMENTS, version=manga_dap_version.spectral_index

            ;-----------------------------------------------------------
            ; Get the weights to use to combine the templates, the
            ; bestfit to the full spectrum, the emission-line model, and
            ; the pixel mask of the fit.
            ; TODO: Perform a more robust check of the sizes involved
            if (size(weights_gndf))[2] ne ntpl then begin
                print, 'Spectral index measurements will use pPXF results'
                weights = weights_ppxf
                bestfit = bestfit_ppxf
                sz = size(bin_flux)
                if n_elements(elo_ew_eml_model) ne 0 then begin
                    eml_model = elo_ew_eml_model        ; Use Enci's model
                endif else $
                    eml_model = dblarr(sz[1], sz[2])    ; No emission-line model
;                fit_mask = obj_fit_mask_ppxf
            endif else begin
                print, 'Spectral index measurements will use GANDALF results'
                weights = weights_gndf
                bestfit = bestfit_gndf
;                fit_mask = obj_fit_mask_gndf
            endelse

            ;-----------------------------------------------------------
            ; Get the resolution-matched, emission-free spectra to use
            ; for the spectral index measurements, both template and
            ; object.  TODO: Make this its own block?
            if MDAP_CAN_USE_SINDX_IMAGES(execution_plan.ofile) eq 1 then begin
                MDAP_DEFINE_OUTPUT, si_bin_wave=si_bin_wave, si_bin_flux=si_bin_flux, $
                                    si_bin_ivar=si_bin_ivar, si_bin_mask=si_bin_mask, $
                                    si_otpl_flux=si_otpl, si_otpl_mask=si_otpl_mask, $
                                    si_botpl_flux=si_botpl, si_botpl_mask=si_botpl_mask
                MDAP_READ_OUTPUT, execution_plan.ofile, si_bin_wave=si_bin_wave, $
                                  si_bin_flux=si_bin_flux, si_bin_ivar=si_bin_ivar, $
                                  si_bin_mask=si_bin_mask, si_otpl_flux=si_otpl_flux, $
                                  si_otpl_mask=si_otpl_mask, si_botpl_flux=si_botpl_flux, $
                                  si_botpl_mask=si_botpl_mask
            endif else begin
                MDAP_SPECTRAL_INDEX_MEASUREMENTS_SPECTRA, bin_indx, wave, sres, bin_flux, $
                                                          bin_ivar, bin_mask, tpl_out_fits, $
                                                          execution_plan, star_kin_interp, $
                                                          gas_kin_interp, velocity_initial_guess, $
                                                          velocity_dispersion_initial_guess, $
                                                          weights, stellar_kinematics, $
                                                          abs_line_keys[execution_plan.abs_par],$
                                                          si_bin_wave, si_bin_flux, si_bin_ivar, $
                                                          si_bin_mask, si_otpl_flux, si_otpl_mask, $
                                                          si_botpl_flux, si_botpl_mask, $
                                                          bestfit=bestfit, eml_model=eml_model, $
                                                          remove_outliers=remove_outliers, $
                                                    moments=execution_plan.analysis_par.moments
;                                                          fit_mask=fit_mask, $
                MDAP_WRITE_OUTPUT, execution_plan.ofile, si_bin_wave=si_bin_wave, $
                                   si_bin_flux=si_bin_flux, si_bin_ivar=si_bin_ivar, $
                                   si_bin_mask=si_bin_mask, si_otpl_flux=si_otpl_flux, $
                                   si_otpl_mask=si_otpl_mask, si_botpl_flux=si_botpl_flux, $
                                   si_botpl_mask=si_botpl_mask, /read_header
            endelse

            ;-----------------------------------------------------------
            ; Perform the measurements
            MDAP_SPECTRAL_INDEX_MEASUREMENTS, abs_par, si_bin_wave, si_bin_flux, $
                                              si_bin_ivar, si_bin_mask, si_otpl_flux, $
                                              si_otpl_mask, si_botpl_flux, si_botpl_mask, $
                                              stellar_kinematics, spectral_index_omitted, $
                                              spectral_index_blue_cont, $
                                              spectral_index_blue_cerr, $
                                              spectral_index_red_cont, $
                                              spectral_index_red_cerr, spectral_index_raw, $
                                              spectral_index_rerr, spectral_index_rperr, $
                                              spectral_index_otpl_blue_cont, $
                                              spectral_index_otpl_red_cont, $
                                              spectral_index_otpl_index, $
                                              spectral_index_botpl_blue_cont, $
                                              spectral_index_botpl_red_cont, $
                                              spectral_index_botpl_index, $
                                              spectral_index_corr, spectral_index_cerr, $
                                              spectral_index_cperr, dbg=dbg, quiet=quiet

            ;-----------------------------------------------------------
            ; Write the data to the output file
            MDAP_WRITE_OUTPUT, execution_plan.ofile, abs_par=abs_par, $
                               abs_line_key=abs_line_keys[execution_plan.abs_par], $
                               spectral_index_omitted=spectral_index_omitted, $
                               spectral_index_blue_cont=spectral_index_blue_cont,  $
                               spectral_index_blue_cerr=spectral_index_blue_cerr, $
                               spectral_index_red_cont=spectral_index_red_cont, $
                               spectral_index_red_cerr=spectral_index_red_cerr, $
                               spectral_index_raw=spectral_index_raw, $
                               spectral_index_rerr=spectral_index_rerr, $
                               spectral_index_rperr=spectral_index_rperr, $
                               spectral_index_otpl_blue_cont=spectral_index_otpl_blue_cont, $
                               spectral_index_otpl_red_cont=spectral_index_otpl_red_cont, $
                               spectral_index_otpl_index=spectral_index_otpl_index, $
                               spectral_index_botpl_blue_cont=spectral_index_botpl_blue_cont, $
                               spectral_index_botpl_red_cont=spectral_index_botpl_red_cont, $
                               spectral_index_botpl_index=spectral_index_botpl_index, $
                               spectral_index_corr=spectral_index_corr, $
                               spectral_index_cerr=spectral_index_cerr, $
                               spectral_index_cperr=spectral_index_cperr, $
                               quiet=quiet, /read_header
        endif

        ;---------------------------------------------------------------
        ; NOTE: This is the last analysis block.  Nothing is used by
        ; later procedures in MANGA_DAP, so this procedure returns
        ; nothing.

END

