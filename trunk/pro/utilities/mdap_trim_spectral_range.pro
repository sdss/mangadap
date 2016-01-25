;+
; NAME:
;       MDAP_TRIM_SPECTRAL_RANGE
;
; PURPOSE:
;       Trim the flux, ivar, mask, sres, and wave vectors to a specified
;       wavelength range.
;
; CALLING SEQUENCE:
;       MDAP_TRIM_SPECTRAL_RANGE, flux, ivar, mask, sres, wave, waverange, $
;                                 /keep_original_wave, /mask_only
;
; INPUTS:
;       flux dblarr[]
;               A 1D or 2D array with the flux at each wavelength.  Defines the
;               dimensionality against which to check ivar and mask.  If 2D, the
;               dispersion axis is assumed to be along rows.  Replaced on output
;               with array truncated to the specified wavelength range.
;       
;       ivar dblarr[]
;               A 1D or 2D array with the inverse variances at each wavelength.
;               Dimensionality must match flux.  If 2D, the dispersion axis is
;               assumed to be along rows.  Replaced on output with array
;               truncated to the specified wavelength range.
;       
;       mask dblarr[]
;               A 1D or 2D array with the mask values at each wavelength.
;               Dimensionality must match flux.  If 2D, the dispersion axis is
;               assumed to be along rows.  Replaced on output with array
;               truncated to the specified wavelength range.
;       
;       sres dblarr[]
;               A 1D or 2D array with the spectral resolution at each
;               wavelength.  Dimensionality must be 1D with the correct number
;               of spectral channels or match the dimensionality of wave/flux.
;               If 2D, the dispersion axis is assumed to be along rows.
;               Replaced on output with array truncated to the specified
;               wavelength range.
;
;       wave dblarr[]
;               A 1D or 2D array with the wavelength coordinates.
;               Dimensionality must match flux or be 1D with the correct number
;               of spectral channels.  If wave is 2D, the dispersion axis is
;               assumed to be along rows.  If /keep_original_wave is not set,
;               the array is replaced on output by an array truncated to the
;               specified wavelength range.
;
;       waverange dblarr[2]
;               A two-element vector defining the starting and ending
;               wavelength.
;       
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;       /keep_original_wave
;               Do not replace the wave array by the trimmed version.   
;
;       /mask_only
;               Do not trim the wavelengths from the vectors/arrays, just mask
;               them
;
; OUTPUT:
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
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;       18 Sep 2014: (KBW) Original implementation
;       19 Sep 2014: (KBW) Allow to only mask the data
;       05 Jan 2015: (KBW) Allow for zero-element where statements
;-
;------------------------------------------------------------------------------

PRO MDAP_TRIM_SPECTRAL_RANGE, $
                flux, ivar, mask, sres, wave, waverange, keep_original_wave=keep_original_wave, $
                mask_only=mask_only

        sz = size(flux)
        nd = sz[0]                              ; Number of dimensions
;       if nd ne 1 and nd ne 2 then $
        if nd ne 1 && nd ne 2 then $
            message, 'Array must be either 1 or 2 dimensions.'

        nc = sz[nd]                             ; Number of spectral channels (works for 1D and 2D)

        szi=size(ivar)
;       if sz[0] ne szi[0] or sz[1] ne szi[1] or sz[2] ne szi[2] then $
        if sz[0] ne szi[0] || sz[1] ne szi[1] || sz[2] ne szi[2] then $
            message, 'Size of ivar must match the size of flux'

        szm=size(mask)
;       if sz[0] ne szm[0] or sz[1] ne szm[1] or sz[2] ne szm[2] then $
        if sz[0] ne szm[0] || sz[1] ne szm[1] || sz[2] ne szm[2] then $
            message, 'Size of mask must match the size of flux'

        szw = size(wave)
        ndw = szw[0]
        if nd lt ndw then $
            message, 'wave vector has higher dimensionality than flux array!'

;       if (ndw eq 1 and nc ne szw[1]) or $
;          (ndw eq 2 and (sz[1] ne szw[1] or nc ne szw[2])) then $
        if (ndw eq 1 && nc ne szw[1]) || $
           (ndw eq 2 && (sz[1] ne szw[1] || nc ne szw[2])) then $
            message, 'Size of wave does not match size of flux'

        szs=size(sres)
        nds = szs[0]
        if nd lt nds then $
            message, 'sres vector has higher dimensionality than flux array!'

        if nds lt ndw then $
            message, 'sres has lower dimensionality than the wave array!'

;       if (nds eq 1 and nc ne szs[1]) or $
;          (nds eq 2 and (sz[1] ne szs[1] or nc ne szs[2])) then $
        if (nds eq 1 && nc ne szs[1]) || $
           (nds eq 2 && (sz[1] ne szs[1] || nc ne szs[2])) then $
            message, 'Size of sres does not match size of flux'

        if ndw eq 1 then begin

            ; Wave vector is 1D so this only has to be done once
            MDAP_SELECT_WAVE, wave, waverange, indx, complement=comp, count=count
            ccount=n_elements(wave)-count

            ; Only applying the mask, so only need to change the mask vector
            if keyword_set(mask_only) then begin

                if ccount ne 0 then begin
                    if nd eq 1 then begin
                        mask[comp] = 1.           ; Mask pixels not within the spectral range
                    endif else $
                        mask[*,comp] = 1.
                endif

                return                                  ; Done
            endif

            if count eq 0 then $
                return

            ; Truncate the arrays
            if nd eq 1 then begin
                flux = flux[indx]
                ivar = ivar[indx]
                mask = mask[indx]
            endif else begin
                flux = flux[*,indx]
                ivar = ivar[*,indx]
                mask = mask[*,indx]
            endelse

            ; The spectral resolution array can be a different size
            if nds eq 1 then begin
                sres = sres[indx]
            endif else $
                sres = sres[*,indx]

            ; Truncate the wavelength array
;           if ~keyword_set(keep_original_wave) and ~keyword_set(mask_only) then $
            if ~keyword_set(keep_original_wave) && ~keyword_set(mask_only) then $
                wave = wave[indx]

            return                                      ; Done
        endif

        ns = sz[1]                              ; Number of spectra, if nd=2
        nu = 0                                  ; Maximum number of unmasked pixels
        for i=0,ns-1 do begin
            MDAP_SELECT_WAVE, wave[i,*], waverange, indx, complement=comp, count=count
            ccount=n_elements(wave)-count
;           if keyword_set(mask_only) and ccount ne 0 then $
            if keyword_set(mask_only) && ccount ne 0 then $
                mask[i,comp] = 1.               ; Only masking pixels

;           ng = n_elements(indx)
            ng = count
            if nu lt ng then $
                nu = ng                         ; New one is larger, so save it
        endfor

        if keyword_set(mask_only) then $
            return                              ; Mask already applied, so we're done

        wave_ = wave                            ; Save the original arrays
        flux_ = flux
        ivar_ = ivar
        mask_ = mask
        sres_ = sres

        wave = dblarr(ns, nu)                   ; Resize output arrays and initialize to 0
        flux = dblarr(ns, nu)
        ivar = dblarr(ns, nu)
        mask = dblarr(ns, nu)
        sres = dblarr(ns, nu)

        for i=0,ns-1 do begin
            ; Select the appropriate wavelengths
            MDAP_SELECT_WAVE, wave_[i,*], waverange, indx, count=ng 
;           ng=n_elements(indx)                                 ; Length of the saved data

            if ng eq 0 then $
                continue

            wave[i,0L:ng-1] = wave_[i,indx]                     ; Copy the data
            flux[i,0L:ng-1] = flux_[i,indx]
            ivar[i,0L:ng-1] = ivar_[i,indx]
            mask[i,0L:ng-1] = mask_[i,indx]
            sres[i,0L:ng-1] = sres_[i,indx]
            if ng lt nu then $
                mask[i,ng:nu-1] = 1.0                           ; Mask any pixels with no data
        endfor

        ; Replace with original wavelengths if requested
        if keyword_set(keep_original_wave) then $
            wave = wave_

END

