;+
; NAME:
;       MDAP_MEAN_GAS_KINEMATICS
;
; PURPOSE:
;       Compute the weighted mean of a set of velocities and velocity
;       dispersions, the propagated error, and the standard error using
;       the weighted standard deviation.
;
; CALLING SEQUENCE:
;       MDAP_MEAN_GAS_KINEMATICS, flux, v, ve, s, se, wv, wvpe, wvse, ws, wspe, wsse, $
;                                 /flux_weighted, /velocity_error_weighted, /quiet
;
; INPUTS:
;       flux dblarr[N]
;               Flux in the measured lines (used as weights if
;               /flux_weighted)
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
;       /flux_weighted
;               Use the provided flux values when weighting the
;               kinematic measurements.  Can be used in combination with
;               /velocity_error_weighted.
;
;       /velocity_error_weighted
;               Use the provided velocity errors when weighting the
;               kinematic measurements (both velocity and velocity
;               dispersion).  Can be used in combination with
;               /flux_weighted.
;
;       /quiet
;               Suppress output to STDOUT
;
; OUTPUT:
;       wv double
;       wvpe double
;       wvse double
;               Weighted-mean velocity, its propagated error, and its
;               standard error (weighted standardard deviation over
;               sqrt(N)).
;
;       ws double
;       wspe double
;       wsse double
;               Weighted-mean velocity dispersion, its propagated error,
;               and its standard error (weighted standardard deviation
;               over sqrt(N)).
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
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;       10 Dec 2014: Pulled from MDAP_GANDALF_WRAP for use with other
;                    procedures by K. Westfall (KBW).
;       10 Jul 2015: (KBW) Allow for error-weighting using the velocity
;                          errors, only or in combination with the flux
;                          weighting.  Allow for the returned error to
;                          be the standard error based on the
;                          error-weighted standard deviation, instead of
;                          the nominal error propagation.
;       13 Jul 2015: (KBW) Instead of allowing it as an option, always
;                          provide both the propagated error (pe) and
;                          the standard error (se)
;       12 Aug 2015: (KBW) Corrected some documentation
;-
;------------------------------------------------------------------------------

PRO MDAP_MEAN_GAS_KINEMATICS, $
                flux, v, ve, s, se, wv, wvpe, wvse, ws, wspe, wsse, flux_weighted=flux_weighted, $
                velocity_error_weighted=velocity_error_weighted, quiet=quiet

        nkin = (size(v))[2]

        ; Set the weights
        wgts = make_array(nkin, /double, value=1.0)
        if keyword_set(flux_weighted) then $
            wgts *= flux
        if keyword_set(velocity_error_weighted) then $
            wgts *= (1.0/ve/ve)
        wsum = total(wgts)

        wv = total(wgts*v)/wsum
        ws = total(wgts*s)/wsum

        if nkin gt 1 then begin
            wvse = sqrt( (total(wgts*v^2)/wsum - wv^2)/(nkin-1) )
            wsse = sqrt( (total(wgts*s^2)/wsum - ws^2)/(nkin-1) )
        endif else begin
            if ~keyword_set(quiet) then $
                print, 'WARNING: Need at least two gas kinematic measurements for standard error.'
            ; Use placeholder values
            wvse = -9999.0d
            wsse = -9999.0d
        endelse

        wvpe = sqrt(total((wgts*ve)^2))/wsum
        wspe = sqrt(total((wgts*se)^2))/wsum
END

