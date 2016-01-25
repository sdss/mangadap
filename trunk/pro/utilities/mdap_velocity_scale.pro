;+
; NAME:
;       MDAP_VELOCITY_SCALE
;
; PURPOSE:
;       Determine the velocity sampling of an input wavelength coordinate
;       vector.  Wavelength vector is assumed to be logarithmically sampled, but
;       its units are in angstroms.
;
; CALLING SEQUENCE:
;       returned=MDAP_VELOCITY_SCALE(wave, /log10)
;
; INPUTS:
;       wave dblarr[C]
;               Wavelength coordinates of each spectral channel C in angstroms.
;               It is expected that the spectrum has been sampled
;               geometrically.
;       
; OPTIONAL KEYWORDS:
;       /LOG10
;               Input spectrum has been sample geometrically using the base 10
;               logarithm, instead of the natural logarithm.
;
; OUTPUT:
;       returned double
;               Velocity scale of the spectrum in km/s.
;
; REVISION HISTORY:
;       18 Sep 2014: (KBW) Original implementation
;-
;------------------------------------------------------------------------------

FUNCTION MDAP_VELOCITY_SCALE, $
                wave, log10=log10

        ; TODO: Check that the wavelength values have a constant velocity
        ;       sampling

        if ~keyword_set(log10) then begin
            dl_over_l = (alog(wave[1])-alog(wave[0]))
        endif else $
            dl_over_l = (alog10(wave[1])-alog10(wave[0]))*alog(10.0d)

        ;TODO: there isn't some idlutils file that defines this?
        c=299792.458d                   ; Speed of light in km/s

        return, c*dl_over_l
END

