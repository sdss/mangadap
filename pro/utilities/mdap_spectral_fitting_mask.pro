;+
; NAME:
;       MDAP_SPECTRAL_FITTING_MASK
;
; PURPOSE:
;       Return the indices of pixels to fit at the base-level (no masks
;       for emission/sky lines) based on the original object mask,
;       template spectra, and starting guess velocities.
;
;       This should only result in leading and trailing pixels masked in
;       the spectra; i.e., no masked regions should be interspersed over
;       the full spectral range.
;
; CALLING SEQUENCE:
;       fit_indx = MDAP_SPECTRAL_FITTING_MASK(obj_wave, tpl_wave, velScale, wave_range_analysis, $
;                                             star_kin_starting_guesses)
;
; INPUTS:
;       obj_wave dblarr[C]
;               Wavelength of all C spectral channels, which is the same
;               for all N spectra.
;
;       tpl_wave dblarr[S]
;               Wavelength of all S spectral channels, which is the same
;               for all template spectra.
;
;       velScale double
;               Velocity scale of each pixel in both the template and
;               object spectra (km/s)
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;
; OUTPUT:
;       The list of pixels in the input obj_wave vector to include in
;       the spectral fitting.
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
;   - Allow maxvel and sigomit to be input parameters?
;
; PROCEDURES CALLED:
;       MDAP_SELECT_WAVE
;       MDAP_SET_INTERSECTION()
;       MDAP_STC()
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;       10 Aug 2015: Pulled from MDAP_SPECTRAL_FITTING by K. Westfall (KBW)
;-
;-----------------------------------------------------------------------

FUNCTION MDAP_SPECTRAL_FITTING_MASK, $
                obj_wave, tpl_wave, velScale, wave_range_analysis, star_kin_starting_guesses

        ; 1. Apply the wavelength range limit, if provided
        now=(size(obj_wave))[1]                             ; Number of object wavelengths
        if n_elements(wave_range_analysis) then begin
            MDAP_SELECT_WAVE, obj_wave, wave_range_analysis, fit_indx, count=count
            if count eq 0 then $
                message, 'wave_range_analysis selects no pixels!'
        endif else $
            fit_indx = indgen(now)

        ; Get the guess minimum and maxium redshift
        c=299792.458d       ; Speed of light in km/s
        maxvel = 400.       ; Maximum range in motions (km/s); including error in guess redshift
        z_min = (min(star_kin_starting_guesses[*,0])-maxvel)/c  ; Minimum redshift
        z_max = (max(star_kin_starting_guesses[*,0])+maxvel)/c  ; Maximum redshift

        ; 2. If the number of template pixels is not >= number of fitted galaxy pixels,
        ;    further limit the blue and red edges of the galaxy spectra
        now=(size(obj_wave[fit_indx]))[1]               ; Number of good object pixels
        ntw=(size(tpl_wave))[1]                         ; Number of template pixels
        if ntw lt now then begin

            ; Indices of wavelengths redward of the redshifted template
            indx=where(obj_wave gt tpl_wave[0]*(1. + z_min), count)
            if count eq 0 then $
                message, 'No overlapping wavelengths between galaxy and template!'

            ; New number of good object pixels
            now_=(size(obj_wave[indx]))[1]
            if ntw lt now_ then $
                indx=indx[0:ntw-1]                      ; Truncate the red edge as well

            ; Merge with current index
            fit_indx = MDAP_SET_INTERSECTION(indx, temporary(fit_indx), count=count)
            if count eq 0 then $
                message, 'No overlapping wavelengths between galaxy and template!'
        endif
        now=(size(obj_wave[fit_indx]))[1]               ; Update number of good object pixels

        ; 3. Limit wavelength range to avoid aliasing problems in the template convolution
        sigomit = 6.                                    ; Number of sigma to mask (+/- 3)
        nalias = fix(sigomit*maxvel/velScale)           ; Number of pixels to mask (maxvel~maxsig)
        print, 'Masking '+MDAP_STC(nalias,/integer)+' pixels at either end of the spectrum to' $
               + ' avoid convolution aliasing.'
        ; Mask to the range that should be unaffected by alias errors
        wave_range_tpl_unalias = [ tpl_wave[nalias]*(1+z_max), tpl_wave[ntw-nalias-1]/(1+z_min) ]
        MDAP_SELECT_WAVE, obj_wave, wave_range_tpl_unalias*(1. + z_min), indx
        ; Merge with current index
        fit_indx = MDAP_SET_INTERSECTION(indx, temporary(fit_indx), count=count)
        if count eq 0 then $
            message, 'No intersection between wave_range_tpl_unalias and fit_indx!'

        return, fit_indx
END

