;+
; NAME:
;       MDAP_COMBINE_SPECTRA
;
; PURPOSE:
;       Produce weighted average spectra using the supplied weights, propagating
;       the inverse variance to the result.
;
; CALLING SEQUENCE:
;       MDAP_COMBINE_SPECTRA, flux, ivar, bin_id, combined_flux, combined_ivar, combined_mask
;
; INPUTS:
;       flux dblarr[N][T]
;               A list of N spectra with T spectral channels.  All spectral
;               changes MUST have the same wavelength.  The program does not
;               account for shifts in wavelength!
;
;       ivar dblarr[N][T]
;               Inverse variance of the flux vectors.
;
;       mask dblarr[N][T]
;               Bad pixel bit mask.
;
;       bin_id intarr[N]
;               Index of the bin for each spectrum.
;
;       wgt dblarr[N]
;               Weight for each spectrum.
;
;       nbinned intarr[B]
;               Number of spectra in each of the B bins.  Length used to get the
;               number of bins.
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;
; OUTPUT:
;       combined_flux dblarr[B][T]
;               Weighted spectra.
;       
;       combined_ivar dblarr[B][T]
;               Inverse variance of weighted spectra, assuming Gaussian errors.
;
;       combined_mask dblarr[B][T]
;               Mask for combined spectra.
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
;       - Account for covariance?
;       - Add some checks for the division by ivar!
;       - allow for spectra to have different spectral resolutions?
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;       15 Sep 2014: Original implementation by K. Westfall (KBW)
;       22 Sep 2014: (KBW) Output mask (just a place-holder for now)
;       04 Dec 2014: (KBW) Account for mask of input spectra, mostly
;                          just confirmed what was in place was already
;                          doing this!
;       16 Mar 2014: (KBW) Allow for noise calibration.
;-
;------------------------------------------------------------------------------

PRO MDAP_COMBINE_SPECTRA, $
                flux, ivar, mask, bin_id, wgt, nbinned, combined_flux, combined_ivar, $
                combined_mask, noise_calib=noise_calib

        sz=size(nbinned)
        nb = sz[0]                                      ; If bins are defined, will be non-zero
        if nb eq 0 then begin                           ; No bins were produced!
            combined_flux=0
            combined_ivar=0
            combined_mask=0
            message, 'No spectra to combine!'
            return
        endif

        nb = sz[1]                                      ; Number of bins
        sz=size(flux)
        nc=sz[2]                                        ; Number of spectral channels

        combined_flux = dblarr(nb, nc)                  ; Initialize the combined spectra
        combined_ivar = dblarr(nb, nc)                  ; All pixel values initialized to 0.0
        combined_nois = dblarr(nb, nc)                  ; Noise vectors, initialized to 0s
        combined_mask = dblarr(nb, nc)                  ; All pixels unmasked (mask=0)

        for i=0,nb-1 do begin
            if nbinned[i] eq 0 then begin               ; No spectra in this bin!
                combined_flux[i,*] = 0.                 ; Set flux to zero
                combined_ivar[i,*] = 1.                 ; Set unity error
                combined_mask[i,*] = 1.                 ; Set to mask all pixels
                continue
            endif

            spec_indx=where(bin_id eq i, count)         ; Indices of spectra in this bin
            if count eq 0 then $
                message, 'No spectra with index '+string(i)
            sumwgt=dblarr(nc)                           ; Re-initialize the weights

            ; TODO: Need to optimize this
            for j=0,nbinned[i]-1 do begin               ; Combine the spectra
                ii = spec_indx[j]                       ; Spectrum to add

                ; Unmasked pixels with defined errors
                gindx = where(ivar[ii,*] gt 0. and mask[ii,*] lt 1., count)
;               if gindx[0] eq -1 then $                ; No valid pixels!
                if count eq 0 then $                    ; No valid pixels!
                    continue

                ; Add the flux
                combined_flux[i,gindx] = combined_flux[i,gindx] + wgt[ii]*flux[ii,gindx]
                ; Propagate the error
                combined_nois[i,gindx] = combined_nois[i,gindx] + wgt[ii]^2/ivar[ii,gindx]
                ; Add the weights to the sum
                sumwgt[gindx] = sumwgt[gindx] + wgt[ii]

            endfor

            gindx = where(sumwgt gt 0., gcount, complement=bindx, ncomplement=bcount)

            if gcount ne 0 then begin
                ; Get the weighted sum and its error for the pixels that were valid
                combined_flux[i,gindx] = combined_flux[i,gindx] / sumwgt[gindx]     ; Normalize
                if n_elements(noise_calib) ne 0 then begin
                    combined_nois[i,gindx] = (MDAP_CALIBRATE_NOISE(sqrt(combined_nois[i,gindx]), $
                                                                   sumwgt[gindx], $
                                                                   noise_calib))^2
                endif
                combined_ivar[i,gindx] = sumwgt[gindx]^2/combined_nois[i,gindx]
                combined_mask[i,gindx] = 0.                                     ;  ... and unmask'em
            endif

            ; Mask pixels without anything added to the sum
            if bcount ne 0 then begin
                combined_flux[i,bindx] = 0.             ; Set to zero flux
                combined_ivar[i,bindx] = 1.             ;  ... with unity error
                combined_mask[i,bindx] = 1.             ;  ... and mask'em
            endif

        endfor

END


