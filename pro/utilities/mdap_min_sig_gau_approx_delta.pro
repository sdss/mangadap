;+
; NAME:
;	MDAP_MIN_SIG_GAU_APPROX_DELTA
;
; PURPOSE:
;	Hard-wired procedure used to set the minimum allowed sigma of a Gaussian
;	before the function should be approximated as a Kronecker delta function.
;
; CALLING SEQUENCE:
;	result=MDAP_MIN_SIG_GAU_APPROX_DELTA(dx)
;
; INPUTS:
;	dx double
;		Discrete sample size
;
; OUTPUT:
;	returned double
;		Returned value is the minimum allowed sigma.
;
; REVISION HISTORY:
;	17 Sep 2014: (KBW) Original implementation based on limit in v0_8 of
;	             MDAP_CONVOL_SIGMA.
;-
;------------------------------------------------------------------------------

FUNCTION MDAP_MIN_SIG_GAU_APPROX_DELTA, $
		dx
	return, (1.0d*dx)
END

