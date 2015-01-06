;+
; NAME:
;       MDAP_NORMALIZE_TEMPLATES
;
; PURPOSE:
;       Calculate and apply a normalization to the templates such that their
;       mean is roughly unity.
;
; CALLING SEQUENCE:
;       MDAP_NORMALIZE_TEMPLATES, tpl_flux, tpl_ivar, tpl_mask, tpl_flux_norm
;
; INPUTS:
;       tpl_flux dblarr[T][S]
;               Array containing T template spectra, each with S spectral
;               channels.  Replaced upon output with the normalized fluxes.
;
;       tpl_ivar dblarr[T][S]
;               Array containing inverse variances in the S channels of the T
;               template spectra.  Replaced upon output with the inverse
;               variances in the normalized flux.
;
;       tpl_mask dblarr[T][S]
;               Bit mask for template spectra.  Value is 0 if pixel is
;               good, value is 1 if it should be masked.
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;
; OUTPUT:
;       tpl_flux_norm double
;               The scalar normalization applied to the template fluxes.
;
; OPTIONAL OUTPUT:
;
; COMMENTS:
;
; EXAMPLES:
;
; TODO:
;       - Use an iterative mean with rejection instead of median?
;
; BUGS:
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;       19 Sep 2014: (KBW) Original Implementation
;       22 Sep 2014: (KBW) Include inverse variances
;       05 Jan 2015: (KBW) Allow for clean behavior when everything masked
;-
;------------------------------------------------------------------------------

PRO MDAP_NORMALIZE_TEMPLATES, $
                tpl_flux, tpl_ivar, tpl_mask, tpl_flux_norm

        indx = where(tpl_mask lt 1, count)
        if count eq 0 then begin
            tpl_flux_norm = 1.
        endif else begin
            tpl_flux_norm = median(tpl_flux[indx], /even)
            tpl_flux = tpl_flux / tpl_flux_norm
            tpl_ivar = (tpl_flux_norm)^2*tpl_ivar
        endelse

END









