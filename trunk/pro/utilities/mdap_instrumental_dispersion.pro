;+
; NAME:
;       MDAP_INSTRUMENTAL_DISPERSION
;
; PURPOSE:
;       Interpolate the instrumental dispersion at a provided set of
;       redshift wavelengths and return a vector with the results.
;
; CALLING SEQUENCE:
;       result = MDAP_INSTRUMENTAL_DISPERSION(wave, sres, restwave, velocity, $
;                                             zero_instr_disp=zero_instr_disp)
;
; INPUTS:
;       wave dblarr[C]
;               Common wavelength at each of C pixels for the spectra.
;
;       sres dblarr[C]
;               Common spectral resolution (R) at each of C pixels for
;               the spectra.
;
;       restwave dblarr[L]
;               Rest wavelengths at which to sample the instrumental
;               dispersion (km/s).
;
;       velocity dblarr[L] or double
;               Velocity shift to apply to convert rest wavelength to
;               observed wavelength (km/s).  It should be either a
;               single value applied to all rest wavelengths, or a
;               wavelength-specific value for the L rest wavelengths
;               provided.
;
; OPTIONAL INPUTS:
;       zero_instr_disp integer
;               Flag to simply return zero instrumental dispersion;
;               0-false, 1 (or anything else) true.
;
; OPTIONAL KEYWORDS:
;
; OUTPUT:
;       Output is the instrumental dispersion at each of the L provided
;       redshifted wavelengths, sigma in km/s.
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
;       09 Jan 2015: (KBW) Original implementation (pulled from
;                          MDAP_SPECTRAL_FITTING)
;-
;------------------------------------------------------------------------------

FUNCTION MDAP_INSTRUMENTAL_DISPERSION, $
                wave, sres, restwave, velocity, zero_instr_disp=zero_instr_disp

        ignore = 0
        if n_elements(zero_instr_disp) ne 0 then begin
            ignore = zero_instr_disp
        endif

        if ignore ne 0 then $
            return, dblarr(n_elements(restwave))

        c=299792.458d                                   ; Speed of light in km/s
        sig2fwhm = 2.0d*sqrt(alog(4.0d))                ; Conversion from sigma to FWHM

        return, (interpol(c/sres, wave, restwave * (1.0d + velocity/c))/sig2fwhm)
END


