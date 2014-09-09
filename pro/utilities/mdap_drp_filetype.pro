;+
; NAME:
;	MDAP_DRP_FILETYPE
;
; PURPOSE:
;	Determine and/or test the validity of a DRP-produced file type.
;
; CALLING SEQUENCE:
;	MDAP_DRP_FILETYPE, flux, type
;
; INPUTS:
;	flux dblarr
;		Flux array read from a DRP-produced fits file.
;
;	type string
;		String representation of the file type.  If provided on input,
;		it is checked against the expected file types, either 'RSS' or
;		'CUBE'.  If not provided or if the input value is invalid, the
;		type is determined by the dimensionality of the read flux array,
;		'RSS' files have 2 dimensions and 'CUBE' files have three.
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
; TODO:
;	- Change to using splog
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

PRO MDAP_DRP_FILETYPE, $
		flux, type

	ndim=size(flux, /n_dimensions)			; Get the number of dimensions
	if ndim ne 2 and ndim ne 3 then begin
	    message, 'ERROR: Cannot determine file type because dimensionality is not 2 or 3.'
	    return
	endif

	defined_type = 1				; Assume type is defined
	if n_elements(type) eq 1 then begin		; Check the type makes sense
	    if type ne 'RSS' and type ne 'CUBE' then begin
		print, 'ERROR: Undefined type!  Must be either RSS or CUBE.'
		defined_type = 0
	    endif
	endif

	; Determine the type because it was not provided or undefined
	if n_elements(type) eq 0 or defined_type eq 0 then begin
	    if ndim eq 2 then begin
		type='RSS'
	    endif else $
		type='CUBE'
	endif

	; Finally check that the input type matches the fits file dimensionality
	if type eq 'CUBE' and ndim eq 2 then begin
	    print, 'WARNING: Input was type=CUBE, but dimensionality suggests RSS. Assuming RSS.'
	    type='RSS'
	endif
	if type eq 'RSS' and ndim eq 3 then begin
	    print, 'WARNING: Input was type=RSS, but dimensionality suggests CUBE. Assuming CUBE.'
	    type='CUBE'
	endif

END



