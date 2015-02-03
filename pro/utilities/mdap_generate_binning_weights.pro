;+
; NAME:
;       MDAP_GENERATE_BINNING_WEIGHTS
;
; PURPOSE:
;       Generate weights used when combining a set of spectra.
;
; CALLING SEQUENCE:
;       MDAP_GENERATE_BINNING_WEIGHTS, signal, noise, wgt, /optimal_weighting
;
; INPUTS:
;       signal dblarr[N]
;               Signal of each spectrum; see MDAP_CALCULATE_SPECTRUM_SN.
;
;       noise dblarr[N]
;               Noise of each spectrum; see MDAP_CALCULATE_SPECTRUM_SN.
;
; OPTIONAL INPUTS:
;       /optimal_weighting
;               Weight the spectra by S/N^2
;
; OPTIONAL KEYWORDS:
;
; OUTPUT:
;       wgt dblarr[N]
;               Weights to use for each spectrum
;
; OPTIONAL OUTPUT:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; TODO:
;       - Change to a function that returns wgt?
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;       15 Sep 2014: (KBW) Original implementation
;       04 Dec 2014: (KBW) Optimal weighting is by S/N^2, not (S/N)^2
;-
;------------------------------------------------------------------------------

PRO MDAP_GENERATE_BINNING_WEIGHTS, $
                signal, noise, wgt, optimal_weighting=optimal_weighting

        if keyword_set(optimal_weighting) then begin
            wgt = signal/(noise)^2                      ; Weight by S/N^2
            return
        endif

        wgt=make_array(n_elements(signal), /double, value=1.0)  ; Uniform weighting

END

