;+
; NAME:
;	MDAP_ABS_LINE_RESOLUTION
;
; PURPOSE:
;	Return a vector with the spectral resolution of the absorption-line
;	index system over a provided set of wavelengths.
;
; CALLING SEQUENCE:
;	result = MDAP_ABS_LINE_RESOLUTION(library_key, wave)
;
; INPUTS:
;	library_key string
;		A keyword for the library.  Currently the only valid value is
;		'LICK'.
;
;	wave dblarr[C]
;		Wavelength in angstroms at which to calculate the spectral
;		resolution (R=lambda/delta lambda).
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;
; OUTPUT:
;	result is a dblarr[C] with the spectral resolution at each provided
;	wavelength.
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
;	04 Nov 2014: (KBW) Original implementation
;-
;------------------------------------------------------------------------------

FUNCTION MDAP_ABS_LINE_RESOLUTION, $
		library_key, wave
	if library_key eq 'LICK' then begin
	    return, (wave/8.4d)			; Lick resolution is 8.4 angstroms
	endif else begin
	    message, 'Unknown library keyword!'
	endelse
END	

