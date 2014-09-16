;+
; NAME:
;	MDAP_SPECTRAL_BIN_AREA
;
; PURPOSE:
;	Calculate the field-of-view area for each spatially binned spectrum.
;
; CALLING SEQUENCE:
;	MDAP_SPECTRAL_BIN_AREA, dx, dy, nbinned, binned_indx, binned_area
;
; INPUTS:
;	dx double
;		Scale arcsec/pixel in X direction
;
;	dy double
;		Scale arcsec/pixel in Y direction
;
;	nbinned intarr[B]
;		Number of spectra in each of the B bins.  Length used to get the
;		number of bins.
;
;	binned_indx intarr[N]
;		Index of the bin for each spectrum.
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;
; OUTPUT:
;	binned_area dblarr[B]
;		Cummulative area of the combined spectrum in each bin.
;
; OPTIONAL OUTPUT:
;
; COMMENTS:
;
; EXAMPLES:
;
; TODO:
;	- How is binned_area used?  Should it include the weights?
;
; BUGS:
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;	15 Sep 2014: (KBW) Original implementation
;-
;------------------------------------------------------------------------------

PRO MDAP_SPECTRAL_BIN_AREA, $
		dx, dy, nbinned, binned_indx, binned_area

	sz=size(nbinned)
	nb = sz[0]					; Number of bins
	if nb eq 0 then begin				; No bins were produced!
	    binned_area=0
	    message, 'No bins produced!'
	    return
	endif
	
	binned_area = dblarr(nb)			; Initialize bin area to 0.0

	for i=0,nb-1 do begin

	    if nbinned[i] eq 0 then $			; No spectra in this bin!
		continue

	    spec_indx=where(binned_indx eq i)			; Indices of spectra in this bin
	    binned_area[i] = n_elements(spec_indx)*dx*dy	; Total area
	endfor

END











