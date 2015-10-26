;+
; NAME:
;       MDAP_GET_PSEUDOCONTINUUM
;
; PURPOSE:
;       Measure the pseudo-continuum of the band using:
;
;                               \int_l0^l1 flux dl
;               continuum =   ---------------------
;                                 \int_l0^l1 dl
;
;       while accounting for the masked pixels (by not integrating over
;       them) and providing a formal estimate of the error.
;
; CALLING SEQUENCE:
;       MDAP_GET_PSEUDOCONTINUUM, wave, flux, ivar, mask, passband, pseudo_continuum, $
;                                 pseudo_continuum_error, wave_integral=wave_integral, err=err, $
;                                 /geometric
;
; INPUTS:
;       wave dblarr[S]
;               Wavelength of every pixel S.  The sampling does NOT need to be
;               linear in wavelength.
;
;       flux dblarr[S]
;               Flux at every pixel S.
;
;       ivar dblarr[S]
;               Inverse variance at every pixel S.
;
;       mask dblarr[S]
;               Bad pixel mask for each pixel S; 0-good, 1-bad.
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;       /geometric
;               The sampling of the wavelength vector is geometric
;               (e.g., logarithmic); required by MDAP_BIN_EDGES to set
;               the edges of the wavelength pixels.
;
; OUTPUT:
;       pseudo_continuum double
;       pseudo_continuum_error double
;               Pseudo-continuum measurement and its propagated error.
;
; OPTIONAL OUTPUT:
;       wave_integral double
;               The denominator in the pseudo-continuum measurement.
;
;       err integer
;               Error flag for calculation, returned as 0 if calculation
;               is successful, 1 otherwise.
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; TODO:
;   - If provided, have the procedure provide the denominator
;     (wave_integral) instead of calculating it.
;
; PROCEDURES CALLED:
;       MDAP_INTEGRATE_PIXELIZED_VALUE
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;       17 Aug 2015: Moved from MDAP_GET_SPECTRAL_INDEX by K. Westfall (KBW), 
;                    for use in other routines, e.g.
;                    MDAP_NONPAR_EMISSION_LINE_MEASUREMENTS
;       26 Oct 2015: (KBW) Edited to reflect changes in
;                          MDAP_INTEGRATE_PIXELIZED_VALUE.
;-
;------------------------------------------------------------------------------

PRO MDAP_GET_PSEUDOCONTINUUM, $
                wave, flux, ivar, mask, passband, pseudo_continuum, pseudo_continuum_error, $
                wave_integral=wave_integral, err=err, geometric=geometric

        n = n_elements(wave)                            ; Number of pixels
        unity = make_array(n, /double, value=1.0d)      ; Used to get denominator (=wave_integral)

        ; Ignore masked regions or regions where the inverse variance is undefined
        sige = sqrt(1.0d/ivar)
        indx = where(ivar le 0.0d, count)
        mask_ = mask
        if count ne 0 then begin
            mask_[indx] = 1.0d
            sige[indx] = 1.0d
        endif

        ; 'err' initialized to 0 in MDAP_INTEGRATE_PIXELIZED_VALUE, so
        ; not necessary to do here
        MDAP_INTEGRATE_PIXELIZED_VALUE, wave, unity, unity, mask_, passband, wave_integral, dummy, $
                                        err=err, geometric=geometric
        if err eq 1 then begin
            pseudo_continuum = -9999.0d
            pseudo_continuum_error = -9999.0d
            wave_integral = -9999.0d
            return
        endif

        MDAP_INTEGRATE_PIXELIZED_VALUE, wave, flux, sige, mask_, passband, flux_integral, $
                                        flux_integral_err, err=err, geometric=geometric
        if err eq 1 then begin
            pseudo_continuum = -9999.0d
            pseudo_continuum_error = -9999.0d
            wave_integral = -9999.0d
            return
        endif

        ; All cases when wave_integral is zero should be taken care of
        ; by having err == 1 above
        pseudo_continuum = flux_integral/wave_integral
        pseudo_continuum_error = flux_integral_err/wave_integral

;        if finite(pseudo_continuum) eq 0 || finite(pseudo_continuum_error) eq 0 then begin
;            print, wave_integral, dummy, flux_integral, flux_integral_err, err
;            stop
;        endif

END


