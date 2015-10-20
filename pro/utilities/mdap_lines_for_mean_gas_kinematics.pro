;+
; NAME:
;       MDAP_LINES_FOR_MEAN_GAS_KINEMATICS
;
; PURPOSE:
;       Determine which lines in the provded list to use to determine
;       the flux-weighted mean kinematics of the ionized gas.  "Good"
;       lines are those with:
;               - positive and finite fluxes and flux errors
;               - positive velocity and velocity dispersion errors
;               - positive and finite line amplitudes
;
;       If omitted is provided, omitted lines are also not considered
;       good. Duh.
;
; CALLING SEQUENCE:
;       result = MDAP_LINES_FOR_MEAN_GAS_KINEMATICS(sol_gas_A, esol_gas_A, sol_gas_F, esol_gas_F, $
;                                                   sol_gas_V, esol_gas_V, sol_gas_S, esol_gas_S,
;                                                   exclude=exclude, omitted=omitted, count=count)
;
; INPUTS:
;       sol_gas_A  dblarr[E]
;       esol_gas_A  dblarr[E]
;               Best-fitting line amplitude and its error for the E
;               emission lines.
;
;       sol_gas_F  dblarr[E]
;       esol_gas_F  dblarr[E]
;               Best-fitting line flux and its error for the E emission
;               lines.
;
;       sol_gas_V  dblarr[E]
;       esol_gas_V  dblarr[E]
;               Best-fitting line velocity and its error for the E
;               emission lines.
;
;       sol_gas_S  dblarr[E]
;       esol_gas_S  dblarr[E]
;               Best-fitting line velocity dispersion and its error for
;               the E emission lines.
;
; OPTIONAL INPUTS:
;       exclude intarr[E]
;               Flag to automatically exclude the line from the mean
;               kinematics (0-no; 1-yes).
;
;       omitted intarr[E]
;               Flag the that line fit was (1) or was not (0) omitted
;               due to a fitting issue.
;
; OPTIONAL KEYWORDS:
;
; OUTPUT:
;       Output is a list of indices with "good" line fits.
;
; OPTIONAL OUTPUT:
;       count int
;               Number of lines to include in mean gas kinematics;
;               length of output vector.
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
;       10 Dec 2014: Pulled from MDAP_GANDALF_WRAP by K. Westfall (KBW)
;                    for use with other procedures.
;       20 Oct 2015: (KBW) Add functionality to exclude specified lines;
;                          see MDAP_AUTO_EXCLUDE_EML_FROM_KIN().
;-
;------------------------------------------------------------------------------

FUNCTION MDAP_LINES_FOR_MEAN_GAS_KINEMATICS, $
                sol_gas_A, esol_gas_A, sol_gas_F, esol_gas_F, sol_gas_V, esol_gas_V, sol_gas_S, $
                esol_gas_S, exclude=exclude, omitted=omitted, count=count

        ; TODO: esol_gas_A, sol_gas_V, sol_gas_S not used!

        neml = n_elements(sol_gas_A)

        if n_elements(exclude) eq 0 then begin
            exclude_ = intarr(neml)
        endif else $
            exclude_ = exclude

        if n_elements(omitted) eq 0 then begin
            omitted_ = intarr(neml)
        endif else $
            omitted_ = omitted

        return, where( exclude_ eq 0 and omitted_ eq 0 $
                       and sol_gas_F gt 0 and esol_gas_F gt 0 and finite(sol_gas_F) eq 1 $
                       and finite(esol_gas_F) eq 1 and esol_gas_V gt 0 and esol_gas_S gt 0 $
                       and sol_gas_A gt 0 and finite(sol_gas_A) eq 1, count )
END


