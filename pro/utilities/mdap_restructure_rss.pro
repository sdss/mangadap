;+
; NAME:
;	MDAP_RESTRUCTURE_RSS
;
; PURPOSE:
;	Restructure the input RSS file such that the spectra are organized along
;	the second dimension.  As organized by the DRP, they are along the first
;	dimension.
;
; CALLING SEQUENCE:
;	MDAP_RESTRUCTURE_RSS, flux, ivar, mask
;
; INPUTS:
;	flux dblarr[T][N]
;		Input RSS data read from a DRP-produced 'RSS' file.  Converted
;		to new format on output.
;
;	ivar dblarr[T][N]
;		Input inverse-variance data read from a DRP-produced 'RSS' file.
;		Converted to new format on output.
;
;	mask dblarr[T][N]
;		Input pixel mask data read from a DRP-produced 'RSS' file.
;		Converted to new format on output.
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
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;	09 Sep 2014: (KBW) Original implementation
;	09 Sep 2014: (KBW) Include mask
;-
;------------------------------------------------------------------------------

PRO MDAP_RESTRUCTURE_RSS, $
		flux, ivar, mask
	flux=transpose(temporary(flux))		; Just transpose the matrices
	ivar=transpose(temporary(ivar))
	mask=transpose(temporary(mask))
END



