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
;                                                   omitted=omitted)
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
;       omitted intarr[E]
;               Flag the that line fit was (1) or was not (0) omitted.
;
; OPTIONAL KEYWORDS:
;
; OUTPUT:
;       Output is a list of indices with "good" line fits.
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
;       10 Dec 2014: Pulled from MDAP_GANDALF_WRAP for use with other
;                    procedures by K. Westfall (KBW).
;-
;------------------------------------------------------------------------------

FUNCTION MDAP_LINES_FOR_MEAN_GAS_KINEMATICS, $
                sol_gas_A, esol_gas_A, sol_gas_F, esol_gas_F, sol_gas_V, esol_gas_V, sol_gas_S, $
                esol_gas_S, omitted=omitted, count=count

        ; TODO: esol_gas_A, sol_gas_V, sol_gas_S not used!

        if n_elements(omitted) ne 0 then begin

            return, where( omitted eq 0 and $
                           sol_gas_F gt 0 and esol_gas_F gt 0 and $
                           finite(sol_gas_F) eq 1 and finite(esol_gas_F) eq 1 and $
                           esol_gas_V gt 0 and esol_gas_S gt 0 and $
                           sol_gas_A gt 0 and finite(sol_gas_A) eq 1, count )
        endif

        return, where( sol_gas_F gt 0 and esol_gas_F gt 0 and $
                       finite(sol_gas_F) eq 1 and finite(esol_gas_F) eq 1 and $
                       esol_gas_V gt 0 and esol_gas_S gt 0 and $
                       sol_gas_A gt 0 and finite(sol_gas_A) eq 1, count )
END


