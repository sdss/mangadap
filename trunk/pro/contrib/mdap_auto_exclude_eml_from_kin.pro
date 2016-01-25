;+
; NAME:
;       MDAP_AUTO_EXCLUDE_EML_FROM_KIN
;
; PURPOSE:
;       Exclude a set of lines from the mean gas kinematics by default
;       depending on the input fitting method.  Lines excluded from Enci
;       Wang's fits are any of the OII lines because the lines are
;       unresolved and the separation of the lines are not fixed.  For
;       Francesco's code, lines that automatically have their separation
;       fixed are excluded (OII, OIII, OI, NII, SII) as is the OII line
;       when the doublet is fit as a single line.
;
;       WARNING: This function MUST match the expected output from
;                MDAP_DEFINE_EMISSION_LINES_ENCI_BELFIORE because the
;                vector indices are HARD-CODED!!!
;
; CALLING SEQUENCE:
;       result = MDAP_AUTO_EXCLUDE_EML_FROM_KIN(eml_par, /belfiore)
;
; INPUTS:
;       eml_par EmissionLine[E]
;               The parameters for each of the E emission lines.  See
;               MDAP_ASSIGN_EMISSION_LINE_PARAMETERS and
;               MDAP_DEFINE_EMISSION_LINES_ENCI_BELFIORE.
;
; OPTIONAL KEYWORDS:
;       /belfiore
;               Use Francesco Belfiore's fitting code; if not set, the
;               code assumes Enci Wang's code was used.
;
; OUTPUT:
;       The result is a vector with flags (0-no;1-yes) to automatically
;       exclude lines from the mean gas kinematics.
;
; REVISION HISTORY:
;       20 Oct 2015: Original implementation by K. Westfall (KBW)
;-
;------------------------------------------------------------------------------

FUNCTION MDAP_AUTO_EXCLUDE_EML_FROM_KIN, $
                eml_par, belfiore=belfiore

        neml = n_elements(eml_par)              ; Number of emission lines fitted
        exclude = intarr(neml)                  ; All included by default

        if keyword_set(belfiore) then begin
            exclude[1] = 1              ; Second OII line
            exclude[2] = 1              ; OII doublet fit as a single line
            exclude[4] = 1              ; First OIII line
            exclude[7] = 1              ; Second OI line
            exclude[8] = 1              ; First NII line
            exclude[10] = 1             ; Second NII line
            exclude[12] = 1             ; Second SII line
        endif else begin
            exclude[0] = 1              ; First OII line
            exclude[1] = 1              ; Second OII line
            exclude[2] = 1              ; OII doublet fit as a single line
        endelse

        return, exclude
END


