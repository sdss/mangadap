;+
; NAME:
;	MDAP_ARRAY_STRING
;
; PURPOSE:
;	Convert a 1D array into a comma-separated list of numbers in a string.
;
; CALLING SEQUENCE:
;	out = MDAP_ARRAY_STRING(arr, /integer)
;
; INPUTS:
;	arr dblarr
;		Array to convert
; OPTIONAL INPUT:
;	integer
;		Flag to output values as a list of integers
;
; OUTPUT:
;	out string
;		String representation
;
; TODO:
;	- What limits do MDAP_STC have in terms of precision?
;	- Somehow this also works with strings
;
; REVISION HISTORY:
;	15 Sep 2015: (KBW) Original implementation
;-
;------------------------------------------------------------------------------
FUNCTION MDAP_ARRAY_STRING, $
		arr, integer=integer

	n = n_elements(arr)
	out='['
	for i=0,n-2 do $
	    out=out+MDAP_STC(arr[i], integer=integer)+', '
	out=out+MDAP_STC(arr[n-1], integer=integer)+']'

	return, out

END
