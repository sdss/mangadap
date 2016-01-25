;+
; NAME:
;       MDAP_LOG_WAVELENGTH_SCALE
;
; PURPOSE:
;       Return the change in log (natural or base 10) wavelength per
;       pixel given a certain velocity scale.
;
; CALLING SEQUENCE:
;       result = MDAP_LOG_WAVELENGTH_SCALE(velscale, /log10)
;
; INPUTS:
;       velscale double
;               Velocity scale (dl*c/l) per pixel.
;       
; OPTIONAL KEYWORDS:
;       /LOG10
;               Input spectrum has been sample geometrically using the base 10
;               logarithm, instead of the natural logarithm.
;
; OUTPUT:
;       result is a double with the change in log (natural or base 10)
;       wavelength per pixel.
;
; REVISION HISTORY:
;       10 Feb 2015: Original implementation by K. Westfall (KBW)
;-
;------------------------------------------------------------------------------

FUNCTION MDAP_LOG_WAVELENGTH_SCALE, $
                velscale, log10=log10

        c=299792.458d                   ; Speed of light in km/s
        dl_over_l = velscale/c

        if keyword_set(log10) then $
            dl_over_l = dl_over_l / alog(10.0d)

        ; This is alog(wave[1])-alog(wave[0]) if not log10, and
        ; alog10(wave[1])-alog10(wave[0]) otherwise
        return, dl_over_l

END

