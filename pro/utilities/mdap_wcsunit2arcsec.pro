;+
; NAME:
;	MDAP_WCSUNIT2ARCSEC
;
; PURPOSE:
;	Convert an input set of coordinates from units provided by the WCS
;	header to arcseconds.
;
; CALLING SEQUENCE:
;	MDAP_WCSUNIT2ARCSEC, unit, x, y
;
; INPUTS:
;	unit string
;		String representation of the units as determined by the CUNIT1
;		and CUNIT2 header keywords.
;
;	x double
;		X-position in units calculated from the WCS header
;
;	y double
;		Y-position in units calculated from the WCS header
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
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
; TODO:
;	- Currently only understands 'degrees'!
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;	15 Sep 2014: (KBW) Original implementation
;-
;------------------------------------------------------------------------------

PRO MDAP_WCSUNIT2ARCSEC, $
	unit, x, y

	if unit eq 'degrees' then begin			; Convert from degrees to arcseconds
	    x = x*3600.
	    y = y*3600.
	endif else begin
	    message, 'ERROR: Unkown unit type'		; Unknown unit type; fault!
	    retall
	endelse

END



