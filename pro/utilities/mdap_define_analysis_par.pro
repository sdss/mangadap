;+
; NAME:
;	MDAP_DEFINE_ANALYSIS_PAR
;
; PURPOSE:
;	Define the analysis parameter structure and its default values for use.
;
; CALLING SEQUENCE:
;	result = MDAP_DEFINE_ANALYSIS_PAR()
;
; REVISION HISTORY:
;	22 Oct 2014: (KBW) Original implementation
;	10 Nov 2014: (KBW) Add the default values
;-
;------------------------------------------------------------------------------
FUNCTION MDAP_DEFINE_ANALYSIS_PAR

	return, { AnalysisPar, moments:2, degree:-1, mdegree:-1, reddening:dblarr(2), $
			       reddening_order:0 }
END

