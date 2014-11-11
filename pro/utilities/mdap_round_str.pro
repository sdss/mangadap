;+
; NAME:
;	MDAP_ROUND_STR
;
; PURPOSE:
;	Produce a string representation of a floating-point value.
;
; CALLING SEQUENCE:
;	result = MDAP_ROUND_STR(value, decimals)
;
; INPUTS:
;	value double
;		The value to convert to a string
;
;	decimals integer
;		Number of digits to keep beyond the decimal point.
;
; OUTPUT:
;	The output is a string representation of 'value'.
;
; PROCEDURES CALLED:
;	MDAP_STC()
;
; REVISION HISTORY:
;	04 Nov 2014: (KBW) Copied from L. Coccato's version and edited
;-
;------------------------------------------------------------------------------

FUNCTION MDAP_ROUND_STR, $
	value, decimals

	str = MDAP_STC(round(value*(10.^decimals))/10.^decimals)
	pos = strpos(str,'.')
	res = strmid(str, 0, pos+decimals+1)
	return, res
END


