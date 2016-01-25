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
;       23 Sep 2015: (KBW) Fixed major issues with the masked
;                          wavelengths; original approach was wrong.
;-
;-----------------------------------------------------------------------

FUNCTION MDAP_SPECTRAL_FITTING_MASK, $
                obj_wave, tpl_wave, velScale, wave_range_analysis, star_kin_starting_guesses

        ; 1. Apply the wavelength range limit, if provided
        now=(size(obj_wave))[1]                             ; Number of object wavelengths
        print, 'Original number of galaxy pixels: ', now
        if n_elements(wave_range_analysis) ne 0 then begin
;            print, wave_range_analysis
            MDAP_SELECT_WAVE, obj_wave, wave_range_analysis, fit_indx, count=count
            if count eq 0 then $
                message, 'wave_range_analysis selects no pixels!'
;            print, now, count
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
        print, 'After selecting the wavelength range to analyze: ', now
        print, 'Number of template pixels: ', ntw
        if ntw lt now then begin

            ; Indices of wavelengths redward of the redshifted template
            indx=where(obj_wave gt (tpl_wave[0]*(1. + z_min)), count)
            if count eq 0 then $
                message, 'No overlapping wavelengths between galaxy and template!'
            print, 'Pixels blueward of redshifted template: ', n_elements(obj_wave)-n_elements(indx)

            ; Merge with current index
            fit_indx = MDAP_SET_INTERSECTION(indx, temporary(fit_indx), count=count)
            if count eq 0 then $
                message, 'No overlapping wavelengths between galaxy and template!'
            now_=(size(obj_wave[fit_indx]))[1]
            print, 'Object pixels after merging with the specified wavelength range: ', now_

            ; New number of good object pixels
            if ntw lt now_ then begin
                fit_indx=fit_indx[0:ntw-1]                      ; Truncate the red edge as well
                print, 'Impose same as number of template pixels: ', n_elements(fit_indx), ntw
            endif
        endif

        ; 3. Limit wavelength range to avoid aliasing problems in the template convolution
        sigomit = 6.                                    ; Number of sigma to mask (+/- 3)
        nalias = fix(sigomit*maxvel/velScale)           ; Number of pixels to mask (maxvel~maxsig)
        ; Mask to the range that should be unaffected by alias errors
        wave_range_tpl_unalias = [ tpl_wave[nalias]*(1+z_max), tpl_wave[ntw-nalias-1]*(1+z_min) ]
        print, 'Mask to these wavelengths to avoid convolution aliasing: ', wave_range_tpl_unalias
        MDAP_SELECT_WAVE, obj_wave, wave_range_tpl_unalias, indx
        ; Merge with current index
        fit_indx = MDAP_SET_INTERSECTION(indx, temporary(fit_indx), count=count)
        if count eq 0 then $
            message, 'No intersection between wave_range_tpl_unalias and fit_indx!'

        print, 'Final wavelength range to fit: ', obj_wave[fit_indx[0]], $
               obj_wave[fit_indx[n_elements(fit_indx)-1]]
        return, fit_indx
END

