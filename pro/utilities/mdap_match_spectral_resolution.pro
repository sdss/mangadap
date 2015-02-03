;+
; NAME:
;       MDAP_MATCH_SPECTRAL_RESOLUTION
;
; PURPOSE:
;       Match the spectral resolution of a set of input spectra to a new
;       spectral resolution.  Both spectral resolutions (existing and target)
;       can both be wavelength dependent.
;
;       In the case where the existing resolution is LOWER than the target
;       resolution, there are three choices involved: (1) match the existing
;       resolution to the target resolution up to some constant offset that must
;       be accounted for in subsequent analyses, (2) trim the spectral range to
;       only those spectral regions where the existing resolution is better than
;       the target resolution, or (3) allow for a wavelength dependent
;       difference in the spectral resolution that must be accounted for in
;       subsequent analyses.  Option 1 is the default behavior; select option 2
;       using the keyword /no_offset; and option 3 is not currently allowed.
;
;       If using option 1, the keyword /variable_offset allows the offset to be
;       different for all input spectra.  If one expects to combine the spectra,
;       the default behavior should be used, forcing all the spectra to have a
;       constant offset.
;
;       Finally, the code masks a number of pixels at the beginning and end of
;       the spectra to remove regions affected by errors in the convolution due
;       to the censoring of the data.  The number of pixels is the FWHM of the
;       largest Gaussian applied in the convolution:
;       ceil(sig2fwhm*max(diff_sig_w)/dw).  This is currently hard-wired and
;       should be tested.
;
;       The detailed match of the spectral resolution should account for any
;       redshift between the existing and target wavelength frame.  E.g. if the
;       target frame is for a set of galaxy spectra and the existing frame is
;       for a set of templates.  This is not accounted for here, assuming that
;       the wavelength vectors of both the existing and target resolutions are
;       *in the same reference frame*!
;
;       **IMPORTANT**: Because the convolution kernal is treated as a Kronecker
;       delta function based on MDAP_MINIMUM_CNVLV_SIGMA and
;       MDAP_MIN_SIG_GAU_APPROX_DELTA, the resolution matching is not perfect
;       and will accrue error if iteratively applied!
;
;       Cannot handle masks that occur throughout the spectrum.  Only leading
;       and trailing masks of the spectral range are considered.
;
; CALLING SEQUENCE:
;       MDAP_MATCH_SPECTRAL_RESOLUTION, flux, ivar, mask, wave, sres, target_wave, target_sres, $
;                                       soff, /no_offset, /variable_offset, /log10, /quiet
;
; INPUTS:
;       flux dblarr[T][S]
;               Array containing T spectra, each with S spectral channels.
;               Replaced upon output with the spectrum with a resolution matched
;               (as best as possible) to the target value.
;
;       ivar dblarr[T][S]
;               Array containing inverse variance of the S spectral channels in
;               each of the T spectra.  Replaced upon output with the inverse
;               variances at the resolution matched spectra.
;
;       mask dblarr[T][S]
;               Bit mask for spectra.  Used only to mask spectral regions not
;               covered by all input spectra.  Value is 0 if pixel is good,
;               value is 1 if it should be masked.
;
;       wave dblarr[T][S] or dblarr[S]
;               Wavelength of each spectral channel S in angstroms for each
;               spectrum T.  Default is that the spectra are linearly sampled in
;               wavelength.  Use the /log10 keyword to set the sampling to be
;               log10-binned in wavelength.
;
;       sres dblarr[T][S] or dblarr[S]
;               Spectral resolution (R=lamda/delta lambda) for of each spectral
;               channel S for each spectrum T.  Replaced upon output with the
;               matched resolution of the spectra.
;
;       target_wave dblarr[C]
;               Wavelength coordinates of each spectral channel C over which the
;               target spectral resolution is sampled.  The vector is expected
;               to have a constant step in log10(lambda); however, the
;               coordinates are in angstroms, not log10(angstroms).
;
;       target_sres dblarr[C]
;               Target spectral resolution (R=lamda/delta lamba) as a function
;               of wavelength.
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;
;       /no_offset
;               Do not allow a baseline offset of the existing and target
;               resolution matching.  Instead, trim the input spectra to the
;               wavelength range where this offset is zero.
;
;       /variable_offset
;               Allow the offset to be different for each input spectrum, instead of
;               forcing the offset to be the same for all spectra.
;
;       /log10
;               Input spectra are log10 binned in wavelength.  Default is
;               linear.
;
;       /quiet
;               Suppress output written to the screen.
;
; OUTPUT:
;       soff dblarr[T]
;               A constant velocity dispersion artificially added to the
;               resolution of each spectrum T to ensure that the convolution
;               used to match the spectral resolution is always robust.
;
; OPTIONAL OUTPUT:
;
; COMMENTS:
;
; EXAMPLES:
;
; TODO:
;       - MDAP_CONVOL_SIGMA is very slow.  Farm this out to a C/C++ program?
;
; BUGS:
;
; PROCEDURES CALLED:
;       MDAP_MINIMUM_CNVLV_SIGMA
;       MDAP_CONVOL_SIGMA
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;       17 Sep 2014: (KBW) Original Implementation
;       23 Oct 2014: (KBW) Documentation specifies that the wavelengths should
;                          be in the same reference frame.
;       03 Nov 2014: (KBW) Copied from MDAP_MATCH_RESOLUTION_TPL2OBJ and
;                          generalized
;-
;------------------------------------------------------------------------------

PRO MDAP_MATCH_SPECTRAL_RESOLUTION_WAVE_WITH_ZERO_OFFSET, $
                positive_offset, unm_wave, waverange, zero_dispersion=zero_dispersion, quiet=quiet

        if n_elements(zero_dispersion) eq 0 then $
            zero_dispersion = 0.0               ; Base level dispersion

        sz = size(positive_offset)
        np=sz[1]                                ; Number of unmasked pixels

        indx=where(positive_offset lt zero_dispersion, count)
;       if indx[0] eq -1 then $
        if count eq 0 then $
            message, 'Template spectrum is at lower resolution at all wavelengths!'

;       if n_elements(indx) eq np then begin
        if count eq np then begin
            if ~keyword_set(quiet) then begin
                print, 'Entire spectrum is valid: ' + MDAP_STC(unm_wave[0]) + 'A - ' + $
                       MDAP_STC(unm_wave[np-1]) + 'A'
            endif
            waverange = [ unm_wave[0], unm_wave[np-1] ]
            return
        end

        if ~keyword_set(quiet) then begin
            print, 'Input valid range is ' + MDAP_STC(unm_wave[0]) + 'A - ' + $
                   MDAP_STC(unm_wave[np-1]) + 'A'
        endif

        ;-------------------------------------------------------------------------------
        ; Find the largest contiguous region of pixels that requires no offset

        ; Start counting from the first valid wavelength
        indx=where(positive_offset gt zero_dispersion, count)  ; Pixels that require an offset
        incw=[ unm_wave[0], unm_wave[np-1] ]            ; Initialize to full range
        j=0                                             ; Start at low wavelength
;       while indx[0] ne -1 and j lt np-1 do begin
        while count ne 0 and j lt np-1 do begin
            j++                                         ; Increment j
            incw[0] = unm_wave[j]                       ; Update wavelength range
            indx=where(positive_offset[j:np-1] gt zero_dispersion, count)      ; Update indices
        endwhile

        ; Start counting from the last valid wavelength
        indx=where(positive_offset gt zero_dispersion, count)    ; Pixels that require an offset
        decw=[ unm_wave[0], unm_wave[np-1] ]            ; Initialize to full range
        j=np-1                                          ; Start at high wavelength
;       while indx[0] ne -1 and j gt 0 do begin
        while count ne 0 and j gt 0 do begin
            j--                                         ; Decrement j
            decw[1] = unm_wave[j]                       ; Update wavelength range
            indx=where(positive_offset[0:j] gt zero_dispersion, count)  ; Update indices
        endwhile

        ; Select the largest region
        if decw[1]-decw[0] gt incw[1]-incw[0] then begin
            if ~keyword_set(quiet) then $
                print, 'DEC: Truncating to ' + MDAP_STC(decw[0]) + 'A - ' + MDAP_STC(decw[1]) + 'A'
            waverange=decw
        endif else begin
            if ~keyword_set(quiet) then $
                print, 'INC: Truncating to ' + MDAP_STC(incw[0]) + 'A - ' + MDAP_STC(incw[1]) + 'A'
            waverange=incw
        endelse

END

PRO MDAP_MATCH_SPECTRAL_RESOLUTION_APPLY, $
                unmasked, wave, flux, ivar, sres, interp_sres, soff, conv, conv_ivar, $
                sigma_mask, log10=log10

        sig2fwhm = 2.0d*sqrt(alog(4.0d))        ; Factor to convert sigma to FWHM
        c=299792.458d                           ; Speed of light in km/s

        ; Variance of Gaussian in angstroms needed to match input and target resolution
        diff_var_w = (wave[unmasked]/sig2fwhm)^2*(1.0d/(interp_sres)^2 $
                                            - 1.0d/(sres[unmasked])^2)

        diff_sig = (c/wave[unmasked])^2*diff_var_w + soff^2     ; Convert to km/s and add offset

        indx = where(diff_sig lt 0, count)              ; Pixels where sigma is effectively 0
        diff_sig = sqrt(diff_sig)                       ; Convert to sigma
        if ~keyword_set(log10) then $
            diff_sig = wave[unmasked]*diff_sig/c        ; Convolution sigma in angstroms
;       if indx[0] ne -1 then $
        if count ne 0 then $
            diff_sig[indx] = 0.0                        ; Set lt 0 to 0 (within range of delta func)

        ; TODO: This convolution is why there can't be intermittent masks in the
        ;       spectra.  If such masks are to be allowed, need to fix how this
        ;       is done.

        ; Perform the convolution, including the variances
        MDAP_CONVOL_SIGMA_WITH_IVAR, wave[unmasked], flux[unmasked], ivar[unmasked], $
                                     wave[unmasked], diff_sig, conv, conv_ivar

        ; Get the size of the pixel in either wavelength or km/s
        if keyword_set(log10) then begin
            dp = mdap_velocity_scale(wave, /log10)      ; Linear step in log10(wavelength)
        endif else $
            dp = wave[1] - wave[0]                      ; Linear step in wavelength

        ; TODO: Do something within mdap_convol_sigma_with_ivar?
        ; Mask the first and last few pixels due to errors in the convolution
        sigma_mask = ceil(sig2fwhm*max(diff_sig)/dp)
END

FUNCTION MDAP_MATCH_SPECTRAL_RESOLUTION_UNMASKED_PIXELS, $
                mask
            unmasked=where(mask lt 1, count);, complement=masked)        ; Select the valid pixels
;           if unmasked[0] eq -1 then $
            if count eq 0 then $
                message, 'ERROR: Entire spectrum masked!'

            ; Ignore masks littered through the full spectral range; just omit the edges
            return, (indgen(unmasked[n_elements(unmasked)-1]+1-unmasked[0])+unmasked[0])
END

; TODO: No checks are done of the input vectors to make sure they're the correct size
PRO MDAP_MATCH_SPECTRAL_RESOLUTION, $
                flux, ivar, mask, wave, sres, target_wave, target_sres, soff, no_offset=no_offset, $
                variable_offset=variable_offset, log10=log10, quiet=quiet

        ; Check the dimensionality of wave and sres
        wave_matrix = 1                         ; Default is dblarr[T][S]
        sz=size(wave)
        if sz[0] eq 1 then $
            wave_matrix = 0                     ; wave is a vector

        sres_matrix = 1                         ; Default is dblarr[T][S]
        sz=size(sres)
        if sz[0] eq 1 then $
            sres_matrix = 0                     ; sres is a vector

        if wave_matrix eq 1 and sres_matrix eq 0 then $
            message, 'Cannot handle wavelength matrix combined with resolution vector.'

        sz=size(flux)
        nt = sz[1]                              ; Number of input spectra
        soff=dblarr(nt)                         ; Initialize the baseline offset to 0.0
        sig2fwhm = 2.0d*sqrt(alog(4.0d))        ; Factor to convert sigma to FWHM
        c=299792.458d                           ; Speed of light in km/s

        sigma_mask = intarr(nt)                 ; Used to mask the first and last set of pixels

        if ~keyword_set(quiet) then $
            print, 'Number of spectra to match resolution: ', nt

        ; Allow wave and sres to be either a vector or matrix
        if wave_matrix eq 0 then $
            wave_ = wave
        if sres_matrix eq 0 then $
            sres_ = sres

        ; Determine any necessary offsets to account for an existing resolution
        ; that is lower than the target resolution
        for i=0,nt-1 do begin
            if ~keyword_set(quiet) then $
                print, '    Assessing template: ', i+1

            ; Allow wave and sres to be either a vector or matrix
            if wave_matrix eq 1 then $
                wave_ = wave[i,*]
            if sres_matrix eq 1 then $
                sres_ = sres[i,*]

            unmasked=MDAP_MATCH_SPECTRAL_RESOLUTION_UNMASKED_PIXELS(reform(mask[i,*]))
            num = n_elements(unmasked)                          ; Number of unmasked pixels
           
            ; Report some numbers to the screen
            if ~keyword_set(quiet) then begin
                print, 'Spectrum length: ', n_elements(reform(flux[i,*]))
                print, 'Number of unmasked pixels: ', num
                zero_vals = where(flux[i,unmasked] eq 0., count)  ; Number of pixels with zero flux
                print, 'Number of zero valued fluxes: ', count
;               if zero_vals[0] eq -1 then begin
;                   print, 'Number of zero valued fluxes: ', 0
;               endif else $
;                   print, 'Number of zero valued fluxes: ', n_elements(zero_vals)
            endif

            ; Object resolution interpolated to template wavelengths
            interp_sres=interpol(target_sres, target_wave, reform(wave_[unmasked]))

;           print, size(obj_sres)
;           mydevice=!D.NAME
;           set_plot, 'PS'
;           device, filename='res_plot.ps'
;           plot, tpl_wave[i,unmasked], tpl_sres[i,unmasked]
;           oplot, tpl_wave[i,unmasked], obj_sres, linestyle=2
;           device, /close
;           set_plot, mydevice
;           stop

            ; Variance (in angstroms) of Gaussian needed to match template and object resolution
            diff_var_w = (wave_[unmasked]/sig2fwhm)^2*(1.0d/(interp_sres)^2 $
                                                        - 1.0d/(sres_[unmasked])^2)
            diff_var_v = (c/wave_[unmasked])^2*diff_var_w               ; Convert to km/s

            ; Determine a constant offset (in km/s) required to ensure that the
            ; Gaussian convolution kernel is always above the minimum width
            if keyword_set(log10) then begin
                dp = mdap_velocity_scale(wave_, /log10)         ; Linear step in velocity

                positive_offset = reform((MDAP_MINIMUM_CNVLV_SIGMA(dp))^2 - diff_var_v)
                zero_dispersion = (MDAP_MIN_SIG_GAU_APPROX_DELTA(dp))^2
            endif else begin
                dp = wave_[1] - wave_[0]                ; Linear step in wavelength

                positive_offset = reform((c*MDAP_MINIMUM_CNVLV_SIGMA(dp)/wave_[unmasked])^2 $
                                         - diff_var_v)
                zero_dispersion = (c*MDAP_MIN_SIG_GAU_APPROX_DELTA(dp)/wave_[unmasked])^2
            endelse

            if ~keyword_set(no_offset) then begin               ; Allow any offset
                soff[i] = sqrt(max([ positive_offset, 0.0d ]))  ; Necessary sigma offset
                if ~keyword_set(quiet) then $
                    print, 'Offset: '+MDAP_STC(soff[i])
            endif else begin                                    ; Force offset to be 0

                soff[i] = 0.                            ; Offset is zero

                ; Find the region with zero offset needed
                MDAP_MATCH_SPECTRAL_RESOLUTION_WAVE_WITH_ZERO_OFFSET, positive_offset, $
                                                                      wave_[unmasked], waverange, $
                                                                     zero_dispersion=zero_dispersion

                ; Mask regions that require an offset
                MDAP_SELECT_WAVE, wave_, waverange, indx, complement=comp, count=count
;               if comp[0] ne -1 then $
                if count ne n_elements(wave_) then $
                    mask[i,comp] = 1.0

            endelse

        endfor

        if ~keyword_set(variable_offset) then begin
            max_offset=max(soff)                        ; Set the offset to the maximum value
            soff[*]=max_offset
            if ~keyword_set(quiet) then $
                print, 'Offset: '+MDAP_STC(max_offset)  ; Alert the user
        endif

        ;-------------------------------------------------------------------------------
        ; Convolve the templates to a common resolution wrt the object spectra
        ; with a common offset
        for i=0,nt-1 do begin
            if ~keyword_set(quiet) then begin
                print, '    Matching resolution of template: ' + MDAP_STC(i+1,/integer) + '/' + $
                       MDAP_STC(nt,/integer)
            endif

            ; List of unmasked pixels
            unmasked=MDAP_MATCH_SPECTRAL_RESOLUTION_UNMASKED_PIXELS(reform(mask[i,*]))
            num=n_elements(unmasked)    ; Number of unmasked pixels

            ; Object resolution interpolated to template wavelengths
            interp_sres=interpol(target_sres, target_wave, wave_[unmasked])

            ; Get the resolution matched spectrum and variances
            MDAP_MATCH_SPECTRAL_RESOLUTION_APPLY, unmasked, wave_, reform(flux[i,*]), $
                                                  reform(ivar[i,*]), sres_, interp_sres, soff[i], $
                                                  conv, conv_ivar, sigma_mask, log10=log10

            ; Save the data
            mask[i,unmasked[0:sigma_mask-1]] = 1.               ; Mask the edges
            mask[i,unmasked[num-1-sigma_mask:num-1]] = 1.

            ; TODO: Set the masked values to something?
            flux[i,unmasked] = conv                             ; Set the flux to the convolution
            ivar[i,unmasked] = conv_ivar                        ; Set the inverse variance

            if sres_matrix eq 1 then $
                sres[i,unmasked] = interp_sres          ; Set the new resolution

        endfor

        if sres_matrix eq 0 then $
            sres = interpol(target_sres, target_wave, wave_)    ; Set the new resolution
END


;       MDAP_SELECT_WAVE, tpl_wave[i,*], decw, indx
;       ;-------------------------------------------------------------------------------
;       ; Trim spectra (and associated vectors) if valid indices has changed
;       nind=n_elements(indx)
;       if nind ne np then begin
;           tpl_wave[i,0:nind-1] = tpl_wave[i,indx]     ; Copy the data
;           tpl_flux[i,0:nind-1] = tpl_flux[i,indx]
;           tpl_sres[i,0:nind-1] = tpl_sres[i,indx]
;           tpl_mask[i,0:nind-1] = tpl_mask[i,indx]
;           if nind lt ns then begin
;               tpl_wave[i,nind:ns-1] = 0.              ; Erase the old data
;               tpl_flux[i,nind:ns-1] = 0.
;               tpl_sres[i,nind:ns-1] = 0.
;               tpl_mask[i,nind:ns-1] = 1.
;           endif


            

