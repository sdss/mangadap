;+
; NAME:
;       MDAP_INDEX_IS_BANDHEAD
;
; PURPOSE:
;       Return a flag (0-false,1-true) if the provided absorption-line parameter
;       set defines a bandhead to be measured, not a full absorption-line index.
;
; CALLING SEQUENCE:
;       result = MDAP_INDEX_IS_BANDHEAD(abs_par)
;
; INPUTS:
;       abs_par AbsorptionIndex
;               Structure that contains the defining parameters of the
;               absorption-line index.  See
;               MDAP_READ_ABSORPTION_LINE_PARAMETERS.
;
; OUTPUT:
;       The result is a flag that the defined index is (1) or is not (0) really
;       the definition of a bandhead measurement.
;
; REVISION HISTORY:
;       05 Nov 2014: (KBW) Original implementation
;-
;------------------------------------------------------------------------------

FUNCTION MDAP_INDEX_IS_BANDHEAD, $
                abs_par

        if abs_par.passband[0] lt 1. and abs_par.passband[1] lt 1. then $
            return, 1                   ; Index is a bandhead

        return, 0                       ; Otherwise, it's not
END

