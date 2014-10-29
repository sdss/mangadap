;+
; NAME:
;	MDAP_DEFINE_ANALYSIS_PAR
;
; PURPOSE:
;	Define the analysis parameter structure for use.
;
; CALLING SEQUENCE:
;	result = MDAP_DEFINE_ANALYSIS_PAR()
;
; REVISION HISTORY:
;	22 Oct 2014: (KBW) Original implementation
;-
;------------------------------------------------------------------------------
FUNCTION MDAP_DEFINE_ANALYSIS_PAR

	return, { AnalysisPar, moments:0, degree:0, mdegree:0, reddening:dblarr(2), $
			       reddening_order:0 }
END

