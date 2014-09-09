;+
; NAME:
;	MDAP_SELECT_WAVE
;
; PURPOSE:
;	Create a vector used to select out a desired wavelength range.
;
; CALLING SEQUENCE:
;	MDAP_SELECT_WAVE, wave, wstart, wend, selected
;
; INPUTS:
;	wave dblarr[]
;		List of wavelength coordinates.
;
;	wrange double[2]
;		Start (index 0) and end (index 1) of the valid wavelength range.
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;
; OUTPUT:
;	selected intarr[]
;		List of indices in 'wave' that are within the indicated
;		wavelength range.
;
; OPTIONAL OUTPUT:
;
; COMMENTS:
;
; EXAMPLES:
;
; TODO:
;	- Use splog
;
; BUGS:
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;	09 Sep 2014: (KBW) Original implementation
;-
;------------------------------------------------------------------------------

PRO MDAP_SELECT_WAVE, $
		wave, wrange, selected

	if n_elements(wrange) ne 0 then begin			; Check if the range is defined
	    if wrange[0] gt wrange[1] then begin		; Check that the order makes sense
		print, 'WARNING: Start>End in wavelength range. Swapping!'
		temp = wrange[1]				; Swap the order
		wrange[1] = wrange[0]
		wrange[0] = temp
	    endif

	    selected=where(wave ge wrange[0] and wave le wrange[1])	; Select valid pixels
	endif else $
	    selected = indgen(n_elements(wave))			; All are valid

END

