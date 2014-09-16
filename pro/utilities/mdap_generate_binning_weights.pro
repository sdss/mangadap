;+
; NAME:
;	MDAP_GENERATE_BINNING_WEIGHTS
;
; PURPOSE:
;	Generate weights used when combining a set of spectra.
;
; CALLING SEQUENCE:
;	MDAP_GENERATE_BINNING_WEIGHTS, signal, noise, wgt, /weight_for_sn
;
; INPUTS:
;	signal dblarr[N]
;		Signal of each spectrum; see MDAP_CALCULATE_SPECTRUM_SN.
;
;	noise dblarr[N]
;		Noise of each spectrum; see MDAP_CALCULATE_SPECTRUM_SN.
;
; OPTIONAL INPUTS:
;	/weight_for_sn
;		Weight the spectra by (S/N)^2
;
; OPTIONAL KEYWORDS:
;
; OUTPUT:
;	wgt dblarr[N]
;		Weights to use for each spectrum
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
;	15 Sep 2014: (KBW) Original implementation
;-
;------------------------------------------------------------------------------

PRO MDAP_GENERATE_BINNING_WEIGHTS, $
		signal, noise, wgt, weight_for_sn=weight_for_sn

	if ~keyword_set(weight_for_sn) then begin
	    wgt = (signal/noise)^2			; Weight by (S/N)^2
	    return
	endif

	ns = n_elements(signal)
	wgt=dblarr(ns)
	wgt[*] = 1.0

END

