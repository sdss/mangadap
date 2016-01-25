;+
; NAME:
;       MDAP_EMISSION_LINE_FIT_BLOCK
;
; PURPOSE:
;   MAIN BLOCK FOR MANGA_DAP:
;
;   If the block is to be performed:
;       - Determine non-parametric properties of a fixed set of emission
;         lines.
;       - Determine the gas kinematics and emission properties using
;         initial guess kinematics and Francesco Belfiore's and Enci
;         Wang's codes.
;       - Writes the results to the output files
;   Else:
;       - Reads and returns only the best-fitting model from Enci Wang's
;         code.
;
; CALLING SEQUENCE:
;       MDAP_EMISSION_LINE_FIT_BLOCK, manga_dap_version, execution_plan, perform_block, $
;                                     star_kin_interp, gas_kin_interp, velocity_initial_guess, $
;                                     velocity_dispersion_initial_guess, wave, sres, nbin, $
;                                     bin_indx, bin_flux, bin_ivar, bin_mask, tpl_wave, $
;                                     bestfit_ppxf, bestfit_gndf, eml_model, stellar_kinematics, $
;                                     elo_ew_eml_model, /quiet, /dbg
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
;       wave dblarr[T]
;               Wavelength of each spectral channel T.
;       
;       sres dblarr[T]
;               Spectral resolution at each wavelength channel T.
;
;       nbin lonarr[B]
;               Number of spectra coadded in each bin.
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
;       tpl_wave dblarr[S]
;               Wavelength of all S spectral channels, which is the same for all
;               template spectra.
;
;       bestfit_ppxf dblarr[B][T]
;               Best fitting spectrum obtained by PPXF for each of the B
;               spectra.
;
;       bestfit_gndf dblarr[B][T]
;               Best fitting spectrum obtained by GANDALF for each of
;               the B spectra.
;
;       eml_model dblarr[B][T]
;               Best-fitting emission-line-only model for each of the B spectra
;               obtained by GANDALF.
;
;       stellar_kinematics dblarr[B][M]
;               The best-fit stellar kinematics (M moments) for each of the B
;               fitted input galaxy spectra.
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;       /quiet
;               Suppress output to the screen.
;
;       /dbg
;               Only perform the analysis on the first provided bin.
;
; OUTPUT:
;       elo_ew_eml_model dblarr[N][T]
;               Best-fitting emission-line-only model determined by the
;               code provided by Enci Wang.
;
; OPTIONAL OUTPUT:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; TODO:
;       - Use the non-parametric results to constrain the fits in the
;         parametric case?  And vice versa?
;       - Add zero baseline flag to the execution plan structure?
;
; PROCEDURES CALLED:
;       MDAP_EMISSION_LINE_ONLY_FIT
;       MDAP_INITIALIZE_GUESS_KINEMATICS
;       MDAP_INSTR_DISPERSION_AT_EMISSION_LINE
;       MDAP_WRITE_OUTPUT
;       MDAP_DEFINE_OUTPUT
;       MDAP_READ_OUTPUT
;       MDAP_VELOCITY_SCALE()
;       MDAP_SPECTRAL_FITTING_MASK()
;       MDAP_NONPAR_EMISSION_LINE_MEASUREMENTS
;
; REVISION HISTORY:
;       01 Feb 2015: Pulled from manga_dap.pro by K. Westfall (KBW)
;       22 Mar 2015: (KBW) Added some comments on number of kinematic
;                          moments
;       10 Jul 2015: (KBW) Corrected for many changed to
;                          MDAP_EMISSION_LINE_ONLY_FIT.  Instrumental
;                          dispersion calculation now done within that
;                          procedure.
;       07 Aug 2015: (KBW) Finished changes to input/output given
;                          changes to MDAP_EMISSION_LINE_ONLY_FIT,
;                          including the new output columns and the
;                          different weighting schemes for the average
;                          kinematics.  Added obj_fit_mask_ppxf and
;                          obj_fit_mask_gndf as parameters.  Added
;                          bestfit_mask to input of
;                          MDAP_EMISSION_LINE_ONLY_FIT, using input
;                          object masks from
;                          MDAP_FULL_SPECTRAL_FIT_BLOCK.
;       10 Aug 2015: (KBW) Previous implementation would have omitted
;                          all emission lines!  'bestfit_mask' now
;                          mimics the input into the MDAP_GANDALF_WRAP,
;                          using the function
;                          MDAP_SPECTRAL_FITTING_MASK.  This should only
;                          mask regions that would not have had the
;                          stellar continuum fit because of the fitted
;                          template.
;       18 Aug 2015: (KBW) Included call to non-parametric analysis.
;       17 Sep 2015: (KBW) Edited to allow for non-zero baseline in
;                          Gaussian fitting; also now returns the
;                          wavelength limits of the window used in the
;                          fit.
;-
;------------------------------------------------------------------------------

PRO MDAP_EMISSION_LINE_FIT_BLOCK, $
                manga_dap_version, execution_plan, perform_block, star_kin_interp, gas_kin_interp, $
                velocity_initial_guess, velocity_dispersion_initial_guess, wave, sres, nbin, $
                bin_indx, bin_flux, bin_ivar, bin_mask, tpl_wave, bestfit_ppxf, bestfit_gndf, $
                eml_model, stellar_kinematics, elo_ew_eml_model, quiet=quiet, dbg=dbg

        ; Check if the emission-line-only spectra should be fit
        if perform_block.eml_only eq 1 then begin

            MDAP_EMISSION_LINE_ONLY_FIT, version=manga_dap_version.emission_line

            ;-----------------------------------------------------------
            ; Get the bestfit continuum spectrum and kinematic guesses
            nb = n_elements(nbin)
            sz = size(bin_flux)

            ; Default is to set the continuum and mask to 0
            bestfit_continuum = dblarr(sz[1], sz[2])

            ; Use pPXF continuum if available
            if (size(bestfit_ppxf))[1] eq sz[1] && (size(bestfit_ppxf))[2] eq sz[2] then begin
                bestfit_continuum = bestfit_ppxf
            endif

            ; Change to the GANDALF result if it exists
            ; TODO: Check the size of the eml_model as well?
            if (size(bestfit_gndf))[1] eq sz[1] && (size(bestfit_gndf))[2] eq sz[2] $
               && (size(eml_model))[1] eq sz[1] && (size(eml_model))[2] eq sz[2] then begin
                bestfit_continuum = bestfit_gndf - eml_model
            endif

            ; Initialize the guess kinematics.  Always done because used
            ; by MDAP_SPECTRAL_FITTING_MASK below, not just to set
            ; inp_kin if the input stellar kinematics have the wrong
            ; size.
            MDAP_INITIALIZE_GUESS_KINEMATICS, n_elements(nbin), execution_plan.analysis_prior, $
                                              star_kin_interp, gas_kin_interp, bin_indx, $
                                              velocity_initial_guess, $
                                              velocity_dispersion_initial_guess, $
                                              star_kin_guesses, gas_kin_guesses

            ; Get the first two moments of the input stellar kinematics
            ; TODO: Check if there are at least two moments?
            if (size(stellar_kinematics))[1] eq sz[1] then begin
                inp_kin = stellar_kinematics[*,0:1]
            endif else begin
                inp_kin = star_kin_guesses[*,0:1]
            endelse

            if n_elements(tpl_wave) ne 0 then begin
                ; Get the mask that excludes regions where the fitted template is not valid
                velScale=MDAP_VELOCITY_SCALE(wave, /log10)  ; km/s/pixel of the input spectra
                fit_indx = MDAP_SPECTRAL_FITTING_MASK(wave, tpl_wave, velScale, $
                                                      execution_plan.wave_range_analysis, $
                                                      star_kin_guesses)
                bestfit_mask = make_array(sz[1], sz[2], /double, value=1.0d)
                bestfit_mask[*,fit_indx] = 0.0d
            endif else $
                bestfit_mask = make_array(sz[1], sz[2], /double, value=0.0d)

            ; Include any intrinsically masked pixels
            indx = where(bin_mask > 0, count)
            if count ne 0 then $
                bestfit_mask[indx] = 1.0d

            ;-----------------------------------------------------------
            ; Perform the non-parametric assessments of the emission lines
            MDAP_NONPAR_EMISSION_LINE_MEASUREMENTS, wave, bin_flux, bin_ivar, bin_mask, $
                                                    bestfit_continuum, bestfit_mask, inp_kin, $
                                                    nonpar_eml_omitted, nonpar_eml_flux_raw, $
                                                    nonpar_eml_flux_rerr, nonpar_eml_mom1_raw, $
                                                    nonpar_eml_mom1_rerr, nonpar_eml_mom2_raw, $
                                                    nonpar_eml_mom2_rerr, nonpar_eml_blue_flux, $
                                                    nonpar_eml_blue_ferr, nonpar_eml_red_flux, $
                                                    nonpar_eml_red_ferr, nonpar_eml_blue_fcen, $
                                                    nonpar_eml_blue_cont, nonpar_eml_blue_cerr, $
                                                    nonpar_eml_red_fcen, nonpar_eml_red_cont, $
                                                    nonpar_eml_red_cerr, nonpar_eml_flux_corr, $
                                                    nonpar_eml_flux_cerr, nonpar_eml_mom1_corr, $
                                                    nonpar_eml_mom1_cerr, nonpar_eml_mom2_corr, $
                                                    nonpar_eml_mom2_cerr, eml_bands=eml_bands, $
                                                    quiet=quiet, dbg=dbg

            MDAP_WRITE_OUTPUT, execution_plan.ofile, eml_bands=eml_bands, $
                               nonpar_eml_omitted=nonpar_eml_omitted, $
                               nonpar_eml_flux_raw=nonpar_eml_flux_raw, $
                               nonpar_eml_flux_rerr=nonpar_eml_flux_rerr, $
                               nonpar_eml_mom1_raw=nonpar_eml_mom1_raw, $
                               nonpar_eml_mom1_rerr=nonpar_eml_mom1_rerr, $
                               nonpar_eml_mom2_raw=nonpar_eml_mom2_raw, $
                               nonpar_eml_mom2_rerr=nonpar_eml_mom2_rerr, $
                               nonpar_eml_blue_flux=nonpar_eml_blue_flux, $
                               nonpar_eml_blue_ferr=nonpar_eml_blue_ferr, $
                               nonpar_eml_red_flux=nonpar_eml_red_flux, $
                               nonpar_eml_red_ferr=nonpar_eml_red_ferr, $
                               nonpar_eml_blue_fcen=nonpar_eml_blue_fcen, $
                               nonpar_eml_blue_cont=nonpar_eml_blue_cont, $
                               nonpar_eml_blue_cerr=nonpar_eml_blue_cerr, $
                               nonpar_eml_red_fcen=nonpar_eml_red_fcen, $
                               nonpar_eml_red_cont=nonpar_eml_red_cont, $
                               nonpar_eml_red_cerr=nonpar_eml_red_cerr, $
                               nonpar_eml_flux_corr=nonpar_eml_flux_corr, $
                               nonpar_eml_flux_cerr=nonpar_eml_flux_cerr, $
                               nonpar_eml_mom1_corr=nonpar_eml_mom1_corr, $
                               nonpar_eml_mom1_cerr=nonpar_eml_mom1_cerr, $
                               nonpar_eml_mom2_corr=nonpar_eml_mom2_corr, $
                               nonpar_eml_mom2_cerr=nonpar_eml_mom2_cerr, /read_header

            ;-----------------------------------------------------------
            ; Perform the fit using Enci's code
            MDAP_EMISSION_LINE_ONLY_FIT, wave, sres, bin_flux, bin_ivar, bin_mask, $
                                         bestfit_continuum, bestfit_mask, inp_kin, $
                                         elo_ew_eml_model, elo_ew_kinematics, $
                                         elo_ew_kinematics_perr, elo_ew_kinematics_serr, $
                                         elo_ew_kinematics_n, elo_ew_omitted, elo_ew_window, $
                                         elo_ew_baseline, elo_ew_base_err, $
                                         elo_ew_kinematics_ind, elo_ew_kinematics_ier, $
                                         elo_ew_sinst, elo_ew_intens, elo_ew_interr, $
                                         elo_ew_fluxes, elo_ew_flxerr, elo_ew_EWidth, $
                                         elo_ew_EW_err, eml_par=emlo_par, quiet=quiet, $
                                         dbg=dbg, /enci;, /zero_base

;           ; Get the instrumental dispersion at the fitted line center
;           MDAP_INSTR_DISPERSION_AT_EMISSION_LINE, wave, sres, emlo_par, elo_ew_omitted, $
;                                                   elo_ew_kinematics_ind, elo_ew_sinst
                
            ;-----------------------------------------------------------
            ; Perform the fit using Francesco's code
            MDAP_EMISSION_LINE_ONLY_FIT, wave, sres, bin_flux, bin_ivar, bin_mask, $
                                         bestfit_continuum, bestfit_mask, inp_kin, $
                                         elo_fb_eml_model, elo_fb_kinematics, $
                                         elo_fb_kinematics_perr, elo_fb_kinematics_serr, $
                                         elo_fb_kinematics_n, elo_fb_omitted, elo_fb_window, $
                                         elo_fb_baseline, elo_fb_base_err, $
                                         elo_fb_kinematics_ind, elo_fb_kinematics_ier, $
                                         elo_fb_sinst, elo_fb_intens, elo_fb_interr, $
                                         elo_fb_fluxes, elo_fb_flxerr, elo_fb_EWidth, $
                                         elo_fb_EW_err, quiet=quiet, dbg=dbg, /belfiore;, /zero_base

;           ; Get the instrumental dispersion at the fitted line center
;           MDAP_INSTR_DISPERSION_AT_EMISSION_LINE, wave, sres, emlo_par, elo_fb_omitted, $
;                                                   elo_fb_kinematics_ind, elo_fb_sinst

            ;-----------------------------------------------------------
            ; Write the data
            MDAP_WRITE_OUTPUT, execution_plan.ofile, emlo_par=emlo_par, $
                               elo_ew_kinematics_avg=elo_ew_kinematics, $
                               elo_ew_kinematics_aer=elo_ew_kinematics_perr, $
                               elo_ew_kinematics_ase=elo_ew_kinematics_serr, $
                               elo_ew_kinematics_n=elo_ew_kinematics_n, $
                               elo_ew_window=elo_ew_window, elo_ew_baseline=elo_ew_baseline, $
                               elo_ew_base_err=elo_ew_base_err, $
                               elo_ew_kinematics_ind=elo_ew_kinematics_ind, $
                               elo_ew_kinematics_ier=elo_ew_kinematics_ier, $
                               elo_ew_sinst=elo_ew_sinst, elo_ew_omitted=elo_ew_omitted, $
                               elo_ew_intens=elo_ew_intens, elo_ew_interr=elo_ew_interr, $
                               elo_ew_fluxes=elo_ew_fluxes, elo_ew_flxerr=elo_ew_flxerr, $
                               elo_ew_EWidth=elo_ew_EWidth, elo_ew_EW_err=elo_ew_EW_err, $
                               elo_ew_eml_model=elo_ew_eml_model, $
                               elo_fb_kinematics_avg=elo_fb_kinematics, $
                               elo_fb_kinematics_aer=elo_fb_kinematics_perr, $
                               elo_fb_kinematics_ase=elo_fb_kinematics_serr, $
                               elo_fb_kinematics_n=elo_fb_kinematics_n, $
                               elo_fb_window=elo_fb_window, elo_fb_baseline=elo_fb_baseline, $
                               elo_fb_base_err=elo_fb_base_err, $
                               elo_fb_kinematics_ind=elo_fb_kinematics_ind, $
                               elo_fb_kinematics_ier=elo_fb_kinematics_ier, $
                               elo_fb_sinst=elo_fb_sinst, elo_fb_omitted=elo_fb_omitted, $
                               elo_fb_intens=elo_fb_intens, elo_fb_interr=elo_fb_interr, $
                               elo_fb_fluxes=elo_fb_fluxes, elo_fb_flxerr=elo_fb_flxerr, $
                               elo_fb_EWidth=elo_fb_EWidth, elo_fb_EW_err=elo_fb_EW_err, $
                               elo_fb_eml_model=elo_fb_eml_model, /read_header

        endif else begin

            ;-----------------------------------------------------------
            ; If the data is already available, the only thing that
            ; needs to be read is the fitted model.  For now, the
            ; preference is for Enci's results.

            MDAP_DEFINE_OUTPUT, elo_ew_eml_model=elo_ew_eml_model
            MDAP_READ_OUTPUT, execution_plan.ofile, elo_ew_eml_model=elo_ew_eml_model

        endelse

END


