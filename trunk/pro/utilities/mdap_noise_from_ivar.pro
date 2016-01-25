;+
; NAME:
;       MDAP_NOISE_FROM_IVAR
;
; PURPOSE:
;       Convert an input inverse variance vector to a noise vector.  Also
;       appropriately adjust the mask to account for pixels without finite or
;       non-zero inverse variance values.
;
; CALLING SEQUENCE:
;       MDAP_NOISE_FROM_IVAR, ivar, mask, sige
;
; INPUTS:
;       ivar dblarr[C]
;               Inverse variance in each of C spectral channels
;
;       mask dblarr[C]
;               Bad pixel mask for spectrum, 0-good, 1-bad
;               ADJUSTED ON OUTPUT
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;
; OUTPUT:
;       sige dblarr[C]
;               Standard deviation of the pixel value based on the inverse
;               variance.
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
;       03 Nov 2014: Original implementation by K. Westfall
;-
;------------------------------------------------------------------------------

PRO MDAP_NOISE_FROM_IVAR, $
                ivar, mask, sige

        nc = n_elements(ivar)                                   ; Number of spectral channels
        indx = where(ivar eq 0 or finite(ivar) eq 0, count, complement=nindx)  ; Invalid variances
        if count eq nc then begin                    ; No pixels are valid
            print, 'All errors are invalid!  Ignoring errors (by setting them to unity).'
            sige = make_array(nc, /double, value=1.0d)
        endif else if count eq 0 then begin                  ; All pixels are valid
            sige = sqrt(1.0d/ivar)                              ; Errors
        endif else begin
            ; TODO: Mask doesn't need to be double
            mask[indx] = 1.0d                                   ; Mask invalid variances
            sige = dblarr(nc, /nozero)                          ; Initialize
            sige[nindx] = sqrt(1.0d/ivar[nindx])                ; Errors
            sige[indx] = 1.0d                                   ; Placeholder value (masked!)
        endelse
END

