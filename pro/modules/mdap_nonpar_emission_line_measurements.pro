;+
; NAME:
;       MDAP_NONPAR_EMISSION_LINE_MEASUREMENTS
;
; PURPOSE:
;       Perform non-parametric assessments of the emission lines.  Much
;       like the emission-line-only fitting code, the fitted lines and
;       their associated bands are hard-coded for now.
;
; CALLING SEQUENCE:
;       MDAP_NONPAR_EMISSION_LINE_MEASUREMENTS, wave, flux, ivar, mask, best_fit_continuum, $
;                                               best_fit_mask, stellar_kinematics, omitted, $
;                                               flux_raw, flux_rerr, mom1_raw, mom1_rerr, $
;                                               mom2_raw, mom2_rerr, blue_flux, blue_ferr, $
;                                               red_flux, red_ferr, blue_fcen, blue_cont, $
;                                               blue_cerr, red_fcen, red_cont, red_cerr, $
;                                               flux_corr, flux_cerr, mom1_corr, mom1_cerr, $
;                                               mom2_corr, mom2_cerr, eml_bands=eml_bands, $
;                                               version=version, /geometric, /quiet, /dbg
;
; INPUTS:
;       wave dblarr[C]
;               Wavelength of all C spectral channels, which is the same
;               for all N spectra.
;
;       flux dblarr[N][C]
;               Object flux for each of N spectra with C spectral channels.
;
;       ivar dblarr[N][C]
;               Inverse variance in object flux for each of N spectra with C
;               spectral channels.
;
;       mask dblarr[N][C]
;               Pixel mask for each of N spectra with C spectral channels.
;
;       best_fit_continuum dblarr[N][C]
;               Best fitting continuum spectrum (obtained by PPXF) for
;               each of the N spectra.
;
;       best_fit_mask dblarr[N][C]
;               Bad pixel mask for pixels in the best-fitting continuum
;               data.  Pixels that are not/are masked have values of
;               0.0/1.0.
;
;       stellar_kinematics dblarr[B][K]
;               The best-fit stellar kinematics for each of the N fitted
;               (or fixed) input galaxy spectra with K fitted moments.
;               stellar_kinematics[*,0] must be  fitted velocity (cz).
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;       /geometric
;               The sampling of the wavelength vector is geometric
;               (e.g., logarithmic); required by MDAP_BIN_EDGES to set
;               the edges of the wavelength pixels.
;
;       /quiet
;           Suppress output to the terminal
;
;       /dbg
;           Run in "debug mode" by only fitting the first provided
;           spectrum.
;
; OUTPUT:
;       omitted intarr[N][E]
;               Flag setting whether or not a set of emission-line
;               measurements should be omitted due to errors
;               (0-good;1-bad).
;
;       flux_raw dblarr[N][E]
;       flux_rerr dblarr[N][E]
;               The integrated flux of the continuum-subtracted spectrum
;               over the defined passband and its error.  
;       
;       mom1_raw dblarr[N][E]
;       mom1_rerr dblarr[N][E]
;       mom2_raw dblarr[N][E]
;       mom2_rerr dblarr[N][E]
;               The first and second moments of the velocity over the
;               defined passband and their errors.
;
;       blue_flux dblarr[N][E]
;       blue_ferr dblarr[N][E]
;       red_flux dblarr[N][E]
;       red_ferr dblarr[N][E]
;               The integrated flux of the continuum-subtracted spectrum
;               over the defined blue and red sidebands, and their
;               errors.
;
;       blue_fcen dblarr[N][E]
;       blue_cont dblarr[N][E]
;       blue_cerr dblarr[N][E]
;       red_fcen dblarr[N][E]
;       red_cont dblarr[N][E]
;       red_cerr dblarr[N][E]
;               The flux-weighted centers and pseudo-continuum levels
;               and errors in the blue and red sidebands used to define
;               a linear continuum underneath the primary passband for
;               the corrected fluxes and velocity moments.
;
;       flux_corr dblarr[N][E]
;       flux_cerr dblarr[N][E]
;               The integrated flux of the continuum-subtracted spectrum
;               over the defined passband and its error, including the
;               adjusted continuum level based on the blue and red
;               pseudo-continuum measurements.  
;        
;       mom1_corr dblarr[N][E]
;       mom1_cerr dblarr[N][E]
;       mom2_corr dblarr[N][E]
;       mom2_cerr dblarr[N][E]
;               The first and second moments of the velocity over the
;               defined passband and their errors , including the
;               adjusted continuum level based on the blue and red
;               pseudo-continuum measurements.
;
; OPTIONAL OUTPUT:
;       eml_bands EmissionLineBand[E]
;               Emission-line band parameter structure used to define
;               the parameters of the non-parametric measurements.
;               These are effectively the same as the SpectralIndex
;               structure; see MDAP_ASSIGN_ABSORPTION_LINE_PARAMETERS.
;
;       version string
;               Module version.  If requested, the module is not executed and
;               only version flag is returned.
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; TODO:
;       - Convert 'omitted' to a bit mask saying why it should be omitted!
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;       18 Aug 2015: Original implementation by K. Westfall (KBW)
;       20 Aug 2015 (KBW): Include redshift factor in flux-weighted
;                          centers of sidebands.
;-
;------------------------------------------------------------------------------

PRO MDAP_EMISSION_LINES_NONPARAMETRIC_QUANTITIES, $
                wave, flux, ivar, flux_err, mask, eml_band, result, geometric=geometric, err=err, $
                quiet=quiet

        ; Initialize the result structure with the default values
        result={  flux_raw:make_array(2, /double, value=-9999.0d), $
                  mom1_raw:make_array(2, /double, value=-9999.0d), $
                  mom2_raw:make_array(2, /double, value=-9999.0d), $
                 blue_flux:make_array(2, /double, value=-9999.0d), $
                  red_flux:make_array(2, /double, value=-9999.0d), $
                 blue_fcen:-9999.0d, $
                 blue_cont:make_array(2, /double, value=-9999.0d), $
                  red_fcen:-9999.0d, $
                  red_cont:make_array(2, /double, value=-9999.0d), $
                 flux_corr:make_array(2, /double, value=-9999.0d), $
                 mom1_corr:make_array(2, /double, value=-9999.0d), $
                 mom2_corr:make_array(2, /double, value=-9999.0d) }

        ; Get the indices of the pixels in each band
        pindx = where(wave ge eml_band.bandpass[0] and wave le eml_band.bandpass[1] $
                      and mask lt 0.5, pcount)
        if pcount eq 0 then begin
            if ~keyword_set(quiet) then $
                print, 'Non-par emission line measurements: No unmasked pixels in main bandpass!'
            err = 1
            return
        endif
        bindx = where(wave ge eml_band.blueside[0] and wave le eml_band.blueside[1] $
                      and mask lt 0.5, bcount)
        rindx = where(wave ge eml_band.redside[0] and wave le eml_band.redside[1] $
                      and mask lt 0.5, rcount)

        ; Get the velocity vector (only approximate, but should be a
        ; reasonable approximation over the range relevant to a single
        ; line).
        c=299792.458d                           ; Speed of light in km/s
        v = c*(wave/eml_band.lambda - 1.0d)

        ; Integrate the continuum-subtracted flux over the main bandpass
        MDAP_INTEGRATE_PIXELIZED_VALUE, wave, flux, flux_err, mask, eml_band.bandpass, fr, fre, $
                                        err=err, geometric=geometric
        if err eq 1 then $
            return
        result.flux_raw[0] = fr
        result.flux_raw[1] = fre

        ; Get the 1st and 2nd moments
        integrand = v*flux
        integrand_err = v*flux_err
        MDAP_INTEGRATE_PIXELIZED_VALUE, wave, integrand, integrand_err, mask, eml_band.bandpass, $
                                        num, nume, err=err, geometric=geometric
        result.mom1_raw[0] = num/result.flux_raw[0]
        result.mom1_raw[1] = sqrt( (nume/result.flux_raw[0])^2 $
                                   + (result.mom1_raw[0]*result.flux_raw[1]/result.flux_raw[0])^2 )

        integrand = v*v*flux
        integrand_err = v*v*flux_err
        MDAP_INTEGRATE_PIXELIZED_VALUE, wave, integrand, integrand_err, mask, eml_band.bandpass, $
                                        num, nume, err=err, geometric=geometric
        result.mom2_raw[0] = num/result.flux_raw[0]
        result.mom2_raw[1] = sqrt( (nume/result.flux_raw[0])^2 $
                                   + (result.mom2_raw[0]*result.flux_raw[1]/result.flux_raw[0])^2 )

        result.mom2_raw[1] = sqrt( (result.mom2_raw[0])^2 + (result.mom1_raw[0])^2 )
        result.mom2_raw[0] = result.mom2_raw[0] - result.mom1_raw[0]^2
        result.mom2_raw[0] = result.mom2_raw[0]/sqrt(abs(result.mom2_raw[0]))
        result.mom2_raw[1] = abs(result.mom2_raw[1]/2.0/result.mom2_raw[0])

        ; Integrate the continuum-subtracted flux over the blue and red side bands
        MDAP_INTEGRATE_PIXELIZED_VALUE, wave, flux, flux_err, mask, eml_band.blueside, fr, fre, $
                                        err=err, geometric=geometric
        if err eq 1 then $
            return
        result.blue_flux[0] = fr
        result.blue_flux[1] = fre

        MDAP_INTEGRATE_PIXELIZED_VALUE, wave, flux, flux_err, mask, eml_band.redside, fr, fre, $
                                        err=err, geometric=geometric
        if err eq 1 then $
            return
        result.red_flux[0] = fr
        result.red_flux[1] = fre

        ;---------------------------------------------------------------
        ; Get the pseudocontinuum in the side bands

        ; Blue side band
        MDAP_GET_PSEUDOCONTINUUM, wave, flux, ivar, mask, eml_band.blueside, fr, fre, err=err, $
                                  geometric=geometric
        if err eq 1 then $
            return                              ; No wavelength integral so return with error
        result.blue_cont[0] = fr
        result.blue_cont[1] = fre

        ; Red side band
        MDAP_GET_PSEUDOCONTINUUM, wave, flux, ivar, mask, eml_band.redside, fr, fre, err=err, $
                                  geometric=geometric
        if err eq 1 then $
            return                              ; No wavelength integral so return with error
        result.red_cont[0] = fr
        result.red_cont[1] = fre

        ; TODO: Set some threshold for using this vs. using the center of the band
        ; Get the flux-weighted centers
        integrand = wave*flux
        integrand_err = wave*flux_err
        MDAP_INTEGRATE_PIXELIZED_VALUE, wave, integrand, integrand_err, mask, eml_band.blueside, $
                                        num, nume, err=err, geometric=geometric
        if abs(result.blue_flux[0]) lt 1e-6 then begin
            result.blue_fcen = total(eml_band.blueside)/2.
        endif else $
            result.blue_fcen = num/result.blue_flux[0]
        MDAP_INTEGRATE_PIXELIZED_VALUE, wave, integrand, integrand_err, mask, eml_band.redside, $
                                        num, nume, err=err, geometric=geometric
        if abs(result.red_flux[0]) lt 1e-6 then begin
            result.red_fcen = total(eml_band.redside)/2.
        endif else $
            result.red_fcen = num/result.red_flux[0]

        ;---------------------------------------------------------------
        ; Calculate the slope and intercept of the continuum line
        slope = (result.red_cont[0] - result.blue_cont[0]) / (result.red_fcen - result.blue_fcen)
        intercept = result.blue_cont[0] - result.blue_fcen*slope

        ; Calculate the continuum level at all wavelengths
        continuum = slope * wave + intercept
;        cindx = where(abs(continuum) eq 0.0, ccount)            ; Pixels with continuum == 0
;        p1 = wave/(rcen-bcen) - bcen/(rcen-bcen)
;        continuum_err = sqrt( ( blue_cont_err * (1.0d - p1))^2 + ( red_cont_err * p1 )^2 )
        integrand = flux-continuum
        ingegrand_err = flux_err

        ; Integrate the continuum-subtracted flux over the main bandpass
        MDAP_INTEGRATE_PIXELIZED_VALUE, wave, integrand, integrand_err, mask, eml_band.bandpass, $
                                        fr, fre, err=err, geometric=geometric
        if err eq 1 then $
            return
        result.flux_corr[0] = fr
        result.flux_corr[1] = fre

        ; Get the 1st and 2nd moments
        integrand = v*(flux - continuum)
        integrand_err = v*flux_err
        MDAP_INTEGRATE_PIXELIZED_VALUE, wave, integrand, integrand_err, mask, eml_band.bandpass, $
                                        num, nume, err=err, geometric=geometric
        result.mom1_corr[0] = num/result.flux_corr[0]
        result.mom1_corr[1] = sqrt( (nume/result.flux_corr[0])^2 $
                                + (result.mom1_corr[0]*result.flux_corr[1]/result.flux_corr[0])^2 )

        integrand = v*v*(flux - continuum)
        integrand_err = v*v*flux_err
        MDAP_INTEGRATE_PIXELIZED_VALUE, wave, integrand, integrand_err, mask, eml_band.bandpass, $
                                        num, nume, err=err, geometric=geometric
        result.mom2_corr[0] = num/result.flux_corr[0]
        result.mom2_corr[1] = sqrt( (nume/result.flux_corr[0])^2 $
                                + (result.mom2_corr[0]*result.flux_corr[1]/result.flux_corr[0])^2 )

        result.mom2_corr[1] = sqrt( (result.mom2_corr[0])^2 + (result.mom1_corr[0])^2 )
        result.mom2_corr[0] = result.mom2_corr[0] - result.mom1_corr[0]^2
        result.mom2_corr[0] = result.mom2_corr[0]/sqrt(abs(result.mom2_corr[0]))
        result.mom2_corr[1] = abs(result.mom2_corr[1]/2.0/result.mom2_corr[0])
END

; Omit measurements if:
;   - MDAP_EMISSION_LINES_NONPARAMETRIC_QUANTITIES returns an error
;   - The rest wavelength falls in a masked pixel, where the mask sets
;     the wavelength limits of the stellar continuum fit.
;   - TODO: ADD if some (largish) fraction of the main, blue, and red bandpasses are masked?
PRO MDAP_EMISSION_LINES_SAVE_NONPARAMETRIC_RESULT, $
                eml_band, rest_wave, mask, voff, i, j, result, err, omitted, flux_raw, flux_rerr, $
                mom1_raw, mom1_rerr, mom2_raw, mom2_rerr, blue_flux, blue_ferr, red_flux, $
                red_ferr, blue_fcen, blue_cont, blue_cerr, red_fcen, red_cont, red_cerr, $
                flux_corr, flux_cerr, mom1_corr, mom1_cerr, mom2_corr, mom2_cerr


        c=299792.458d                           ; Speed of light in km/s
        z = (1.0d + voff/c)                     ; So as to offset the flux-weighted band centers

        flux_raw[i,j] = result.flux_raw[0]
        flux_rerr[i,j] = result.flux_raw[1]
        mom1_raw[i,j] = result.mom1_raw[0] + voff
        mom1_rerr[i,j] = result.mom1_raw[1]
        mom2_raw[i,j] = result.mom2_raw[0]
        mom2_rerr[i,j] = result.mom2_raw[1]
        blue_flux[i,j] = result.blue_flux[0]
        blue_ferr[i,j] = result.blue_flux[1]
        red_flux[i,j] = result.red_flux[0]
        red_ferr[i,j] = result.red_flux[1]
        blue_fcen[i,j] = result.blue_fcen*z
        blue_cont[i,j] = result.blue_cont[0]
        blue_cerr[i,j] = result.blue_cont[1]
        red_fcen[i,j] = result.red_fcen*z
        red_cont[i,j] = result.red_cont[0]
        red_cerr[i,j] = result.red_cont[1]
        flux_corr[i,j] = result.flux_corr[0]
        flux_cerr[i,j] = result.flux_corr[1]
        mom1_corr[i,j] = result.mom1_corr[0] + voff
        mom1_cerr[i,j] = result.mom1_corr[1]
        mom2_corr[i,j] = result.mom2_corr[0]
        mom2_cerr[i,j] = result.mom2_corr[1]

        if err eq 1 || mask[(sort(abs(rest_wave-eml_band.lambda)))[0]] gt 0.5 then $
            omitted[i,j] = 1
END

PRO MDAP_NONPAR_EMISSION_LINE_MEASUREMENTS, $
                wave, flux, ivar, mask, best_fit_continuum, best_fit_mask, stellar_kinematics, $
                omitted, flux_raw, flux_rerr, mom1_raw, mom1_rerr, mom2_raw, mom2_rerr, blue_flux, $
                blue_ferr, red_flux, red_ferr, blue_fcen, blue_cont, blue_cerr, red_fcen, $
                red_cont, red_cerr, flux_corr, flux_cerr, mom1_corr, mom1_cerr, mom2_corr, $
                mom2_cerr, eml_bands=eml_bands, version=version, quiet=quiet, dbg=dbg

        version_module = '0.1'                  ; Version number
        if n_elements(version) ne 0 then begin  ; If version is defined
            version = version_module            ; ... set it to the module version
            return                              ; ... and return without doing anything
        endif

        ; Get the emission line bands
        eml_bands = MDAP_DEFINE_NONPAR_EMISSION_LINE_BANDS()

        ; Convert the inverse variance to an error, adjusting the mask for bad values
        indx = where(ivar gt 0.0d, count, complement=nindx)
        if count eq 0 then $
            message, 'All pixels are masked!'

        flux_err = ivar
        flux_err[indx] = sqrt(1.0d/ivar[indx])
        if nindx[0] ne -1 then begin
            flux_err[nindx] = 1.0d
            mask[nindx] = 1.0d
        endif

        ; Allocate the output arrays
        sz = size(flux)
        nobj = sz[1]                            ; Number of object spectra
        npix = sz[2]                            ; Number of pixels in the spectra

        neml = (size(eml_bands))[1]             ; Number of fitted lines

        omitted = intarr(nobj, neml)
        flux_raw = dblarr(nobj, neml)
        flux_rerr = dblarr(nobj, neml)
        mom1_raw = dblarr(nobj, neml)
        mom1_rerr = dblarr(nobj, neml)
        mom2_raw = dblarr(nobj, neml)
        mom2_rerr = dblarr(nobj, neml)
        blue_flux = dblarr(nobj, neml)
        blue_ferr = dblarr(nobj, neml)
        red_flux = dblarr(nobj, neml)
        red_ferr = dblarr(nobj, neml)
        blue_fcen = dblarr(nobj, neml)
        blue_cont = dblarr(nobj, neml)
        blue_cerr = dblarr(nobj, neml)
        red_fcen = dblarr(nobj, neml)
        red_cont = dblarr(nobj, neml)
        red_cerr = dblarr(nobj, neml)
        flux_corr = dblarr(nobj, neml)
        flux_cerr = dblarr(nobj, neml)
        mom1_corr = dblarr(nobj, neml)
        mom1_cerr = dblarr(nobj, neml)
        mom2_corr = dblarr(nobj, neml)
        mom2_cerr = dblarr(nobj, neml)

        ; If debugging, just analyze the first spectrum
        if keyword_set(dbg) then $
            nobj = 1

        c=299792.458d                           ; Speed of light in km/s
        for i=0,nobj-1 do begin

            ; Set the wavelength vector to the rest wavelength
            rest_wave = wave/(1+stellar_kinematics[i,0]/c)

            ; Remove the stellar continuum and copy the mask
            galaxy_eml_only = reform(flux[i,*] - best_fit_continuum[i,*])
            galaxy_mask = reform(best_fit_mask[i,*])

            ; Measure all the provided lines
            for j=0,neml-1 do begin

                MDAP_EMISSION_LINES_NONPARAMETRIC_QUANTITIES, rest_wave, galaxy_eml_only, $
                                                              reform(ivar[i,*]), $
                                                              reform(flux_err[i,*]), $
                                                              reform(mask[i,*]), eml_bands[j], $
                                                              result, err=err, quiet=quiet, $
                                                              /geometric

                MDAP_EMISSION_LINES_SAVE_NONPARAMETRIC_RESULT, eml_bands[j], rest_wave, $
                                                               galaxy_mask, $
                                                               stellar_kinematics[i,0], i, j, $
                                                               result, err, omitted, flux_raw, $
                                                               flux_rerr, mom1_raw, mom1_rerr, $
                                                               mom2_raw, mom2_rerr, blue_flux, $
                                                               blue_ferr, red_flux, red_ferr, $
                                                               blue_fcen, blue_cont, blue_cerr, $
                                                               red_fcen, red_cont, red_cerr, $
                                                               flux_corr, flux_cerr, mom1_corr, $
                                                               mom1_cerr, mom2_corr, mom2_cerr
            endfor

            ; TODO: Adjust this function and provide equivalent width
            ; measurements for these non-parametric fluxes as well?

            ;MDAP_MEASURE_EMISSION_LINE_EQUIV_WIDTH, eml_par, $
            ;                                        emission_line_kinematics_individual[i,*,0], $
            ;                                        emission_line_kinematics_individual[i,*,1], $
            ;                                        rest_wave, galaxy_no_eml, $
            ;                                        emission_line_fluxes[i,*], $
            ;                                        emission_line_fluxes_err[i,*], equiv_width, $
            ;                                        equiv_width_err

            ; TODO: Provide weighted mean of the velocity moments?

            ; TODO: Calculate the instrumental dispersion at the moment
            ; 1 location?

        endfor

END



