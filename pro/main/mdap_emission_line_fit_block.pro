;+
; NAME:
;       MDAP_EMISSION_LINE_FIT_BLOCK
;
; PURPOSE:
;   MAIN BLOCK FOR MANGA_DAP:
;
;   If the block is to be performed:
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
;                                     bin_indx, bin_flux, bin_ivar, bin_mask, bestfit_ppxf, $
;                                     bestfit_gndf, eml_model, stellar_kinematics, $
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
; PROCEDURES CALLED:
;       MDAP_EMISSION_LINE_ONLY_FIT
;       MDAP_INITIALIZE_GUESS_KINEMATICS
;       MDAP_INSTR_DISPERSION_AT_EMISSION_LINE
;       MDAP_WRITE_OUTPUT
;       MDAP_DEFINE_OUTPUT
;       MDAP_READ_OUTPUT
;
; REVISION HISTORY:
;       01 Feb 2015: Pulled from manga_dap.pro by K. Westfall (KBW)
;       22 Mar 2015: (KBW) Added some comments on number of kinematic
;                          moments
;-
;------------------------------------------------------------------------------

PRO MDAP_EMISSION_LINE_FIT_BLOCK, $
                manga_dap_version, execution_plan, perform_block, star_kin_interp, gas_kin_interp, $
                velocity_initial_guess, velocity_dispersion_initial_guess, wave, sres, nbin, $
                bin_indx, bin_flux, bin_ivar, bin_mask, bestfit_ppxf, bestfit_gndf, eml_model, $
                stellar_kinematics, elo_ew_eml_model, quiet=quiet, dbg=dbg
                
                
        ; Check if the emission-line-only spectra should be fit
        if perform_block.eml_only eq 1 then begin

            MDAP_EMISSION_LINE_ONLY_FIT, version=manga_dap_version.emission_line

            ;-----------------------------------------------------------
            ; Get the bestfit continuum spectrum and kinematic guesses
            nb = n_elements(nbin)
            sz = size(bin_flux)

            ; Default is to set the continuum to 0
            bestfit_continuum = dblarr(sz[1], sz[2])

            ; Use pPXF continuum if available
;           if (size(bestfit_ppxf))[1] eq sz[1] and (size(bestfit_ppxf))[2] eq sz[2] then $
            if (size(bestfit_ppxf))[1] eq sz[1] && (size(bestfit_ppxf))[2] eq sz[2] then $
                bestfit_continuum = bestfit_ppxf

            ; Change to the GANDALF result if it exists
            ; TODO: Check the size of the eml_model as well?
;           if (size(bestfit_gndf))[1] eq sz[1] and (size(bestfit_gndf))[2] eq sz[2] then $
            if (size(bestfit_gndf))[1] eq sz[1] && (size(bestfit_gndf))[2] eq sz[2] then $
                bestfit_continuum = bestfit_gndf - eml_model

            ; Get the input stellar kinematics
            ; TODO: Check if there are at least two moments?
            if (size(stellar_kinematics))[1] eq sz[1] then begin
                inp_kin = stellar_kinematics[*,0:1]
            endif else begin
                MDAP_INITIALIZE_GUESS_KINEMATICS, n_elements(nbin), execution_plan.analysis_prior, $
                                                  star_kin_interp, gas_kin_interp, bin_indx, $
                                                  velocity_initial_guess, $
                                                  velocity_dispersion_initial_guess, $
                                                  star_kin_guesses, gas_kin_guesses

                ; Only save the first two moments
                ; TODO: This is done explicitly even though
                ; MDAP_INITIALIZE_GUESS_KINEMATICS actually now only
                ; provides these two moments (as of 22 Mar 2015)
                inp_kin = star_kin_guesses[*,0:1]
            endelse

            ;-----------------------------------------------------------
            ; Perform the fit using Enci's code
            MDAP_EMISSION_LINE_ONLY_FIT, wave, bin_flux, bin_ivar, bin_mask, bestfit_continuum, $
                                         inp_kin, elo_ew_eml_model, elo_ew_kinematics, $
                                         elo_ew_kinematics_err, elo_ew_omitted, $
                                         elo_ew_kinematics_ind, elo_ew_kinematics_ier, $
                                         elo_ew_intens, elo_ew_interr, elo_ew_fluxes, $
                                         elo_ew_flxerr, elo_ew_EWidth, elo_ew_EW_err, $
                                         eml_par=emlo_par, quiet=quiet, dbg=dbg, /enci

            ; Get the instrumental dispersion at the fitted line center
            MDAP_INSTR_DISPERSION_AT_EMISSION_LINE, wave, sres, emlo_par, elo_ew_omitted, $
                                                    elo_ew_kinematics_ind, elo_ew_sinst
                
            ;-----------------------------------------------------------
            ; Perform the fit using Francesco's code
            MDAP_EMISSION_LINE_ONLY_FIT, wave, bin_flux, bin_ivar, bin_mask, bestfit_continuum, $
                                         inp_kin, elo_fb_eml_model, elo_fb_kinematics, $
                                         elo_fb_kinematics_err, elo_fb_omitted, $
                                         elo_fb_kinematics_ind, elo_fb_kinematics_ier, $
                                         elo_fb_intens, elo_fb_interr, elo_fb_fluxes, $
                                         elo_fb_flxerr, elo_fb_EWidth, elo_fb_EW_err, $
                                         quiet=quiet, dbg=dbg, /belfiore

            ; Get the instrumental dispersion at the fitted line center
            MDAP_INSTR_DISPERSION_AT_EMISSION_LINE, wave, sres, emlo_par, elo_fb_omitted, $
                                                    elo_fb_kinematics_ind, elo_fb_sinst

            ;-----------------------------------------------------------
            ; Write the data
            MDAP_WRITE_OUTPUT, execution_plan.ofile, emlo_par=emlo_par, $
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

            ;-----------------------------------------------------------
            ; If the data is already available, the only thing that
            ; needs to be read is the fitted model.  For now, the
            ; preference is for Enci's results.

            MDAP_DEFINE_OUTPUT, elo_ew_eml_model=elo_ew_eml_model
            MDAP_READ_OUTPUT, execution_plan.ofile, elo_ew_eml_model=elo_ew_eml_model

        endelse

END


