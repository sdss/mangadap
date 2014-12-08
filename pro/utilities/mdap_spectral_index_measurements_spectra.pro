;+
; NAME:
;       MDAP_SPECTRAL_INDEX_MEASUREMENTS_SPECTRA
;
; PURPOSE:
;       Prepares the binned spectra for the spectral index measurments, and
;       creates the unbroadened and broadened optimal templates.
;       
;       For the binned spectra:
;
;           1. Determine the error spectrum using either the difference with
;           respect to the input best-fitting model (/rms, bestfit=bestfit) or
;           using the input inverse variance vector.
;
;           2. If requested, replace the outlying pixels with the best-fit model
;           values (remove_outliers=remove_outliers, bestfit=bestfit).
;
;           3. If requested, remove a model of the emission lines
;           (eml_model=eml_model)
;
;           4. If requested, replace the masked region of the object spectra
;           with the value of the best-fit model.
;
;               CAREFUL!!: At least for MDAP_SPECTRAL_FITTING, some of the
;               masked pixels are masked because the MODEL is bad (i.e.
;               unobserved spectral regions of the redshifted template), not the
;               data.
;
;           5. Calculate the spectral resolution of the system over the
;           REST-FRAME wavelengths of the galaxy spectra using the MEDIAN
;           REDSHIFT obtained from the input stellar_kinematics matrix, and
;           match the resolution of the spectra to this.
;
;       For the templates:
;
;           1. Get the optimal template using the provided weights.
;
;           2. Get the LOSVD parameters using the input stellar kinematics.
;           This reverts the input velocity (cz) to the pixel-based kinematics.
;           TODO: Have a keyword that turns this off?
;           
;           3. Create the broadened optimal template using the optimal template
;           and the pixel-based LOSVD.  For use with
;           MDAP_SPECTRAL_INDEX_MEASUREMENTS, it is expected that this will
;           shift the spectrum to match the redshift of the galaxy spectrum!
;
;           4. Interpolate both the unbroadened and broadened templates to the
;           wavelength vector of the object spectra.
;
; CALLING SEQUENCE:
;       MDAP_SPECTRAL_INDEX_MEASUREMENTS_SPECTRA, wave, sres, flux, ivar, mask, tpl_lib_fits, $
;                                                 weights, stellar_kinematics, abs_line_key, $
;                                                 wave_indx, flux_indx, ivar_indx, mask_indx, $
;                                                 otpl_indx, botpl_indx, redshift, $
;                                                 bestfit=bestfit, eml_model=eml_model, $
;                                                 fit_mask=fit_mask, $
;                                                 remove_outliers=remove_outliers,
;                                                 moments=moments, /replace_masked, /oversample, $
;                                                 /rms
;
; INPUTS:
;       wave dblarr[C]
;               Wavelength at each of C spectral channels of the input spectra.
; 
;       sres dblarr[C]
;               Spectral resolution at each of the C spectral channels of the
;               input spectra.
; 
;       flux dblarr[B][C]
;               Flux in the B observed spectra with C spectral channels each.
;
;       ivar dblarr[B][C]
;               Inverse variance in the flux of the observed spectra.
;       
;       mask dblarr[B][C]
;               Bad pixel mask (0-good, 1-bad) for the object spectra.
;
;       tpl_lib_fits string
;               File name for the fits file with the resolution- and
;               sampling-matched template library.  This should have been
;               previously created using MDAP_CREATE_DAP_TEMPLATE_LIBRARY.
;
;       weights dblarr[B][T]
;               Weights that when applied to the template library provides the
;               optimal template for the fit to each of the B object spectra.
;
;       stellar_kinematics dblarr[B][M]
;               Best-fitting stellar kinematics for the B spectra, with LOSVDs
;               that have M moments.  The first element of the LOSVD MUST be the
;               redshift (cz).  If not provided, the redshift is assumed to be
;               0.
;
;       abs_line_key string
;               Keyword signifying the spectral index system, which is used to
;               set the spectral resolution.  See MDAP_ABS_LINE_RESOLUTION().
;
; OPTIONAL INPUTS:
;       bestfit dblarr[B][C]
;               Best-fitting spectrum to the object spectra.
;
;       eml_model dblarr[B][C]
;               Best-fitting emission-line model for the object spectra (can be
;               all zeros if no emission lines were fit or exist)
;       
;       fit_mask intarr[B][C]
;               Bad pixel mask (0-good, 1-bad) of the fit to the object spectra.
;               See MDAP_SPECTRAL_FITTING.  If provided, this is set to
;               mask_indx for use during the spectral index measurements.
;
;       remove_outliers double
;               Replace pixels that deviate more than (remove_outliers) times
;               the error (based on either the RMS difference with respect to
;               the best-fit model or directly from the inverse variance value)
;               with the best-fitting model value.  REQUIRES the bestfit array.

;
;       moments integer
;               Number of kinematic moments.  By default this is determined
;               based on the length of the second dimension of the
;               stellar_kinematics array.
;
; OPTIONAL KEYWORDS:
;       /replace_masked
;               Replace the masked pixels in the object spectra with the
;               best-fitting model.  Requires bestfit and fit_mask.
;               
;               CAREFUL!!: At least for MDAP_SPECTRAL_FITTING, some of the
;               masked pixels are masked because the MODEL is bad (i.e.
;               unobserved spectral regions of the redshifted template), not the
;               data.
;
;       /oversample
;               Force the convolution to oversample the template and the LOSVD
;               when creating the broadened template spectrum. See MDAP_PPXF.
;
;       /rms
;               Define the error in the spectrum (on output) by the RMS
;               difference between the data and the model.  REQUIRES the bestfit
;               array.
;
; OUTPUT:
;       wave_indx dblarr[D]
;               Wavelength vector for the resolution-matched spectra.
;
;       flux_indx dblarr[B][D]
;               Flux of the resolution-matched spectra.
;
;       ivar_indx dblarr[B][D]
;               Inverse variance of the resolution-matched spectra.
;
;       mask_indx dblarr[B][D]
;               Bad pixel mask (0-good, 1-bad) for the resolution-matched
;               spectra.
;
;       otpl_indx dblarr[B][D]
;               Optimal template spectrum, resolution-matched to the
;               spectral-index system.
;
;       botpl_indx dblarr[B][D]
;               Optimal template spectrum, resolution-matched to the
;               spectral-index system, and broadened by the best-fitting LOSVD
;               for each of the B spectra.
;
;       redshift dblarr[B]
;               The redshift of the object spectra (flux_indx) and the broadened
;               optimal template with with respect to the unbroadened optimal
;               template.
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
;       - Output the spectral resolution vector as well?
;       - Allow for the masked pixels to be replaced by the best fit data!
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;       12 Nov 2014: (KBW) Original implementation
;-
;------------------------------------------------------------------------------

PRO MDAP_SPECTRAL_INDEX_MEASUREMENTS_SPECTRA, $
                wave, sres, flux, ivar, mask, tpl_lib_fits, weights, stellar_kinematics, $
                abs_line_key, wave_indx, flux_indx, ivar_indx, mask_indx, otpl_indx, botpl_indx, $
                median_redshift, bestfit=bestfit, eml_model=eml_model, fit_mask=fit_mask, $
                remove_outliers=remove_outliers, moments=moments, replace_masked=replace_masked, $
                oversample=oversample, rms=rms

        ; Cannot get error from RMS if bestfit is not provided
        if keyword_set(rms) and n_elements(bestfit) eq 0 then $
            message, 'Must provide best-fitting spectrum to calculate error using RMS.'
        ; Cannot replace outliers if bestfit is not provided
        if n_elements(remove_outliers) ne 0 and n_elements(bestfit) eq 0 then $
            message, 'Must provide best-fitting spectrum to replace outliers.'
        ; Cannot replace the masked pixels with the model if the model is not provided
        if keyword_set(replace_masked) and n_elements(bestfit) eq 0 then $
            message, 'Must provide best-fitting spectrum to replace masked pixels.'

        ;-----------------------------------------------------------------------
        ; Manipulate the object data -------------------------------------------
        ; 0. Do NOT alter the input spectra, so copy over to some working arrays
        wave_indx = wave
        flux_indx = flux
        ivar_indx = ivar
        mask_indx = mask
        sres_indx = sres

        sz = size(flux)
        nspec = sz[1]                                   ; Number of spectra
        if keyword_set(rms) or n_elements(remove_outliers) ne 0 then begin
            for i=0,nspec-1 do begin
                ;  1. Determine the error spectrum
                mask_ = reform(mask[i,*])
                if keyword_set(rms) then begin                  ; Use the fit residual
                    ; Only use the fitted pixels if using the residual to get the error
                    if n_elements(fit_mask) ne 0 then $
                        mask_ = double(reform(fit_mask[i,*]))
                    MDAP_NOISE_FROM_RESIDUAL, reform(flux[i,*]), mask_, reform(bestfit[i,*]), sige
                endif else $                                    ; Use the inverse variance vector
                    MDAP_NOISE_FROM_IVAR, reform(ivar[i,*]), mask_, sige

                ; TODO: Don't like this switching back and forth between ivar and sige
                mask_indx[i,*] = mask_
                ivar_indx[i,*] = (sige)^(-2)

                ;  2. Replace the outliers with the model values

                ; TODO: This is different from what Lodo did.  He based the
                ; outliers on a measurement of sigma within the full spectral
                ; range of the index.

                if n_elements(remove_outliers) ne 0 then begin
                    flux_ = reform(flux_indx[i,*])
                    MDAP_REPLACE_OUTLIERS, flux_, reform(bestfit[i,*]), sige, remove_outliers
                    flux_indx[i,*] = flux_
                endif
            endfor
        endif

        ; 3. Remove the emission lines from the object spectra, if provided
        if n_elements(eml_model) ne 0 then $
            flux_indx = flux_indx - eml_model

        ; 4. Replace masked pixels with the best-fit model
        if keyword_set(replace_masked) then begin
            indx = where(mask gt 0.5)                   ; Select masked pixels
            if indx[0] ne -1 then $
                flux_indx[indx] = bestfit[indx]         ; Set them to the best-fit model
        endif

        ; 5. Match the spectral resolution of the object spectra to the
        ; absorption-line index system AT THE MEDIAN REDSHIFT OF THE SPECTRA
        c=299792.458d                                   ; Speed of light in km/s

        median_redshift = median(stellar_kinematics[*,0])/c     ; Median redshift
;       print, median_redshift

        wave_indx_rf = wave_indx/(median_redshift + 1.0d)
        abs_sres = MDAP_ABS_LINE_RESOLUTION(abs_line_key, wave_indx_rf) ; Spectral resolution

        ; TODO: This may perform some interpolation of sres_si_system onto the range
        ;       of wave_indx (?)
        MDAP_MATCH_SPECTRAL_RESOLUTION, flux_indx, ivar_indx, mask_indx, wave_indx_rf, sres_indx, $
                                        wave_indx_rf, abs_sres, soff_indx, /no_offset

;       ; 5. Mask leading and trailing pixels omitted during the spectral fit
;       if n_elements(fit_mask) ne 0 then begin
;           for i=0,nspec-1 do begin                            ; For each spectrum
;               indx = where(fit_mask eq 0)                     ; Select masked pixels
;               if indx[0] ne -1 then begin
;                   start_wave = wave[indx[0]]                  ; Beginning
;                   end_wave = wave[indx[n_elements(indx)-1]]   ; ... and end of primary wave range
;                   iindx = where(wave_indx lt start_wave and wave_indx gt end_wave)
;                   if iindx[0] ne -1 then $
;                       mask_indx[i,iindx] = 1.0                ; Mask leading and trailing pixels
;               endif
;           endfor
;       endif
        ;-----------------------------------------------------------------------
        ;-----------------------------------------------------------------------

        ;-----------------------------------------------------------------------
        ; Build the template data ----------------------------------------------

        npix = (size(flux_indx))[2]             ; Number of pixels
        otpl_indx = dblarr(nspec, npix)         ; Allocate the arrays
        botpl_indx = dblarr(nspec, npix)

        ; Read the resolution-matched templates
        MDAP_READ_RESAMPLED_TEMPLATES, tpl_lib_fits, tpl_wave, tpl_flux, tpl_ivar, tpl_mask, $
                                       tpl_sres, tpl_soff

        if n_elements(moments) eq 0 then $
            moments = (size(stellar_kinematics))[2]     ; Number of elements in stellar kinematics
        velScale = MDAP_VELOCITY_SCALE(wave, /log10)    ; Velocity scale of the object data

        for i=0,nspec-1 do begin

            ; 1. Get the optimal template using the provided weights
            optimal_template = reform(weights[i,*] # tpl_flux)  ; Create the optimal template

            ; 2. Revert the input stellar kinematics back to the pixel-based kinematics
            sol = reform(stellar_kinematics[i,0:moments-1])
            vel = sol[0]
            vel_err = 1.0
            MDAP_REVERT_PIXEL_KINEMATICS, vel, vel_err
            sol[0] = vel

            ; 3. Get the broadened optimal template.  It is expected that this will
            ; shift the spectrum to match the redshift of the galaxy spectrum!
            losvd_optimal_template = MDAP_GET_BROADENED_TEMPLATE(optimal_template, sol, velScale, $
                                                                 moments=moments, $
                                                                 oversample=oversample)

            ; 4. Interpolate both the unbroadened and broadened templates to the
            ; wavelength vector of the object spectra.
            otpl_indx[i,*] = interpol(optimal_template, tpl_wave, wave_indx)
            botpl_indx[i,*] = interpol(losvd_optimal_template, tpl_wave, wave_indx)

            ; TODO: Masks for botpl should be the same as fit_mask; however,
            ; this will not be true of otpl.  How do I account for this?

        endfor
        ;-----------------------------------------------------------------------
        ;-----------------------------------------------------------------------


END



