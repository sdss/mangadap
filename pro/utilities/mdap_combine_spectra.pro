;+
; NAME:
;	MDAP_COMBINE_SPECTRA
;
; PURPOSE:
;	Produce weighted average spectra using the supplied weights, propagating
;	the inverse variance to the result.
;
; CALLING SEQUENCE:
;	MDAP_COMBINE_SPECTRA, flux, ivar, bin_id, combined_flux, combined_ivar, combined_mask
;
; INPUTS:
;	flux dblarr[N][T]
;		A list of N spectra with T spectral channels.  All spectral
;		changes MUST have the same wavelength.  The program does not
;		account for shifts in wavelength!
;
;	ivar dblarr[N][T]
;		Inverse variance of the flux vectors.
;
;	mask dblarr[N][T]
;		Bad pixel bit mask.
;
;	bin_id intarr[N]
;		Index of the bin for each spectrum.
;
;	wgt dblarr[N]
;		Weight for each spectrum.
;
;	nbinned intarr[B]
;		Number of spectra in each of the B bins.  Length used to get the
;		number of bins.
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;
; OUTPUT:
;	combined_flux dblarr[B][T]
;		Weighted spectra.
;	
;	combined_ivar dblarr[B][T]
;		Inverse variance of weighted spectra, assuming Gaussian errors.
;
;	combined_mask dblarr[B][T]
;		Mask for combined spectra.
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
;	- Account for covariance?
;	- Add some checks for the division by ivar!
;	- include pixel mask
;	- allow for spectra to have different spectral resolutions?
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;	15 Sep 2014: (KBW) Original implementation
;	22 Sep 2014: (KBW) Output mask (just a place-holder for now)
;-
;------------------------------------------------------------------------------

PRO MDAP_COMBINE_SPECTRA, $
		flux, ivar, mask, bin_id, wgt, nbinned, combined_flux, combined_ivar, combined_mask

	sz=size(nbinned)
	nb = sz[0]					; If bins are defined, will be non-zero
	if nb eq 0 then begin				; No bins were produced!
	    combined_flux=0
	    combined_ivar=0
	    combined_mask=0
	    message, 'No spectra to combine!'
	    return
	endif

	nb = sz[1]					; Number of bins
	sz=size(flux)
	nc=sz[2]					; Number of spectral channels

	combined_flux = dblarr(nb, nc)			; Initialize the combined spectra
	combined_ivar = dblarr(nb, nc)			; All pixel values initialized to 0.0
	combined_mask = dblarr(nb, nc)			; All pixels unmasked (mask=0)

	for i=0,nb-1 do begin
	    if nbinned[i] eq 0 then begin		; No spectra in this bin!
		combined_flux[i,*] = 0.
		combined_ivar[i,*] = 1.
		combined_mask[i,*] = 1.
		continue
	    endif

	    spec_indx=where(bin_id eq i)		; Indices of spectra in this bin
	    sumwgt=dblarr(nc)				; Re-initialize the weights

	    ; TODO: Need to optimize this
	    for j=0,nbinned[i]-1 do begin		; Combine the spectra
		ii = spec_indx[j]			; Spectrum to add

		def_indx = where(ivar[ii,*] gt 0.)
		sz=size(def_indx)
		if sz[0] eq 0 then continue
;		print, def_indx

		; Add the flux
		combined_flux[i,def_indx] = combined_flux[i,def_indx] + wgt[ii]*flux[ii,def_indx]
		; Propagate the error
		combined_ivar[i,def_indx] = combined_ivar[i,def_indx] + wgt[ii]^2/ivar[ii,def_indx]
		; Add the weights to the sum
		sumwgt[def_indx] = sumwgt[def_indx] + wgt[ii]

	    endfor

	    def_indx = where(sumwgt gt 0., complement=undef_indx)

	    combined_flux[i,def_indx] = combined_flux[i,def_indx] / sumwgt[def_indx]	; Normalize
	    combined_ivar[i,def_indx] = sumwgt[def_indx]^2 / combined_ivar[i,def_indx]	; Prop. err

	    combined_flux[i,undef_indx] = 0.		; Set to zero flux
	    combined_ivar[i,undef_indx] = 1.		;  ... with unity error
	    combined_mask[i,undef_indx] = 1.		;  ... and mask'em

	endfor

END

