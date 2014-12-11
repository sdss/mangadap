;+
; NAME:
;       MDAP_MEAN_GAS_KINEMATICS
;
; PURPOSE:
;       Compute the weighted mean of a set of velocities and velocity
;       dispersions, and compute the propagated error.
;
; CALLING SEQUENCE:
;       MDAP_MEAN_GAS_KINEMATICS, wgts, v, ve, s, se, wv, wve, ws, wse
;
; INPUTS:
;       wgts dblarr[N]
;               Weights to use for the weighted mean.
;
;       v dblarr[N]
;       ve dblarr[N]
;               List of velocities to combine and their errors.
;
;       s dblarr[N]
;       se dblarr[N]
;               List of velocity dispersions to combine and their
;               errors.
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;
; OUTPUT:
;       wv double
;       wve double
;               Weighted-mean velocity and its propagated error.
;
;       ws double
;       wse double
;               Weighted-mean velocity dispersion and its propagated
;               error.
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
;   - Allow for error-based weighting
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;       10 Dec 2014: (KBW) Pulled from MDAP_GANDALF_WRAP for use with
;                          other procedures.
;-
;------------------------------------------------------------------------------

PRO MDAP_MEAN_GAS_KINEMATICS, $
                wgts, v, ve, s, se, wv, wve, ws, wse

        wsum = total(wgts)
        wv = total(wgts*v)/wsum
;       wve = mean(ve)
        wve = sqrt(total((wgts*ve)^2))/wsum
        ws = total(wgts*s)/wsum
;       wse = mean(se)
        wse = sqrt(total((wgts*se)^2))/wsum

END

