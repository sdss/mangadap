;+
; NAME: sgn
;
; PURPOSE: Return the sign of the argument. 
;	   If the argument is zero, return zero.
;
; CALLING SEQUENCE:  out=sgn(in)
;
; PARAMETERS:
;		in	input number
;		out	-/0/+
;
; RETURN TYPE: INTEGER                               
;
; HISTORY: Drafted by A.McAllister, Feb. 1993.
;
;-
FUNCTION mdap_sgn,in

 return, fix(in gt 0) - fix(in lt 0)

end
