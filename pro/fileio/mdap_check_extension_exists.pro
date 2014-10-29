;+
; NAME:
;	MDAP_CHECK_EXTENSION_EXISTS
;
; PURPOSE:
;	Check if an extension exists using FXPOSIT.
;
; CALLING SEQUENCE:
;	result = MDAP_CHECK_EXTENSION_EXISTS(file, exten)
;
; INPUTS:
;	file string
;		Fits file name
;
;	exten string
;		Extension name to search for.
;
; OUTPUT:
;	Result is a flag that the extension does (1) or does not (0) exist.
;
; PROCEDURES CALLED:
;	FXPOSIT()
;
; REVISION HISTORY:
;	28 Oct 2014: (KBW) Original implementation (pulled from
;			   MDAP_WRITE_OUTPUT)
;-
;------------------------------------------------------------------------------

FUNCTION MDAP_CHECK_EXTENSION_EXISTS, $
		file, exten

	errmsg = ''				; Catch error messages
	ext=fxposit(file, exten, errmsg=errmsg)	; Attempt to find the extension
	if ext ne -1 then begin			; If found, ...
	    free_lun, ext			;	free the LUN
	    return, 1				; return TRUE
	endif
	return, 0				; return FALSE
END

