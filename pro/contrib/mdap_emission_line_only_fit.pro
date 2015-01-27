;+
; NAME:
;       MDAP_EMISSION_LINE_FIT_ENCI
;
; PURPOSE:
;
; CALLING SEQUENCE:
;       MDAP_EMISSION_LINE_FIT_ENCI, wave, flux, ivar, mask, best_fit_continuum, $
;                                    stellar_kinematics, eml_model, emission_line_kinematics, $
;                                    emission_line_kinematics_err, emission_line_omitted, $
;                                    emission_line_kinematics_individual, $
;                                    emission_line_kinematics_individual_err, $
;                                    emission_line_intens, emission_line_intens_err, $
;                                    emission_line_fluxes, emission_line_fluxes_err, $
;                                    emission_line_EW, emission_line_EW_err, $
;                                    eml_par=eml_par
;
; INPUTS:
;       wave dblarr[C]
;               Wavelength of all C spectral channels, which is the same for all
;               N spectra.
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
;       stellar_kinematics dblarr[B][K]
;               The best-fit stellar kinematics for each of the N fitted
;               (or fixed) input galaxy spectra with K fitted moments.
;               stellar_kinematics[*,0] must be  fitted velocity (cz).
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;
; OUTPUT:
;       eml_model dblarr[N][C]
;               Best-fitting emission-line-only model for each of the N
;               spectra.
;
;       emission_line_kinematics dblarr[N][2]
;       emission_line_kinematics_err dblarr[N][2]
;               Flux weighted means of the best-fit V and sigma for all
;               "well-fit" emission lines in the N galaxy spectra, and
;               the propagated error.
;
;       emission_line_omitted intarr[N][9]
;               Flag setting whether or not an emission-line was fit for
;               all 9 lines fit by Enci's code (0-fit;1-omitted).
;
;       emission_line_kinematics_individual dblarr[N][9][2]
;       emission_line_kinematics_individual_err dblarr[N][9][2]
;               Kinematics (V and sigma) and errors for each of the 9
;               fitted emission lines.
;
;       emission_line_intens dblarr[N][9]
;       emission_line_intens_err dblarr[N][9]
;               Best-fitting emission-line intensities and errors of the
;               9 fitted emission lines for each of the N galaxy
;               spectra.

;
;       emission_line_fluxes  dblarr[N][9]
;       emission_line_fluxes_err dblarr[N][9]
;               Gaussian line fluxes and errors for the 9 fitted
;               emission lines in the N input galaxy spectra. 
;
;       emission_line_EW dblarr[N][9]
;       emission_line_EW_err dblarr[N][9]
;               Equivalent widths and errors for the 9 fitted emission
;               lines for each of the N input galaxy spectra.  These
;               equivalent widths are computed by taking the ratio of
;               the emission_line_fluxes and the median value of the
;               emission-line-free spectrum within 5 and 10 sigma of the
;               emission line, where 'sigma' is the fitted velocity
;               dispersion
;
; OPTIONAL OUTPUT:
;       eml_par EmissionLine[9]
;               Emission-line parameter structure, of the same form used
;               by GANDALF, but with the prameters of the 9 lines fitted
;               by Enci's code.
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
;       10 Dec 2014: (KBW) Original implementation based on code
;                          provided by Enci, Lin, and Francesco.
;       22 Jan 2014: (KBW) Fixed bug in how I saved Francesco's
;                          velocities.
;-
;------------------------------------------------------------------------------

PRO MDAP_SAVE_EMISSION_LINE_FIT_ENCI_SINGLE, lambda, voff, result, kinematics, kinematics_err, $
                                             intens, intens_err, fluxes, fluxes_err, i, j

        c=299792.458d                           ; Speed of light in km/s
        kinematics[i,j,0] = (result[0]/lambda - 1.0d)*c + voff
        kinematics_err[i,j,0] = result[3]*c/lambda
        kinematics[i,j,1] = result[1]*c/result[0]
        kinematics_err[i,j,1] = kinematics[i,j,1]*sqrt( (result[3]/result[0])^2 + $
                                                        (result[4]/result[1])^2 )
        intens[i,j] = result[2]/sqrt(2.0*!pi)/result[1]
        intens_err[i,j] = intens[i,j]*sqrt( (result[5]/result[2])^2 + (result[4]/result[1])^2 )
        fluxes[i,j] = result[2]
        fluxes_err[i,j] = result[5]
END

PRO MDAP_SAVE_EMISSION_LINE_FIT_ENCI, eml_par, voff, result, omitted, kinematics, kinematics_err, $
                                      intens, intens_err, fluxes, fluxes_err, i

        j = 0

        omitted[i,j] = result.OII3727[2] gt 0.0 ? 0 : 1
        if omitted[i,j] eq 1 then begin
            MDAP_SAVE_EMISSION_LINE_SET_OMITTED, kinematics, kinematics_err, intens, intens_err, $
                                                 fluxes, fluxes_err, i, j
        endif else begin
            MDAP_SAVE_EMISSION_LINE_FIT_ENCI_SINGLE, eml_par[j].lambda, voff, result.OII3727, $
                                                     kinematics, kinematics_err, intens, $
                                                     intens_err, fluxes, fluxes_err, i, j
        endelse
        j = j+1

        omitted[i,j] = result.Hb4861[2] gt 0.0 ? 0 : 1
        if omitted[i,j] eq 1 then begin
            MDAP_SAVE_EMISSION_LINE_SET_OMITTED, kinematics, kinematics_err, intens, intens_err, $
                                                 fluxes, fluxes_err, i, j
        endif else begin
            MDAP_SAVE_EMISSION_LINE_FIT_ENCI_SINGLE, eml_par[j].lambda, voff, result.Hb4861, $
                                                     kinematics, kinematics_err, intens, $
                                                     intens_err, fluxes, fluxes_err, i, j
        endelse
        j = j+1

        omitted[i,j] = result.OIII4959[2] gt 0.0 ? 0 : 1
        if omitted[i,j] eq 1 then begin
            MDAP_SAVE_EMISSION_LINE_SET_OMITTED, kinematics, kinematics_err, intens, intens_err, $
                                                 fluxes, fluxes_err, i, j
        endif else begin
            MDAP_SAVE_EMISSION_LINE_FIT_ENCI_SINGLE, eml_par[j].lambda, voff, result.OIII4959, $
                                                     kinematics, kinematics_err, intens, $
                                                     intens_err, fluxes, fluxes_err, i, j
        endelse
        j = j+1

        omitted[i,j] = result.OIII5007[2] gt 0.0 ? 0 : 1
        if omitted[i,j] eq 1 then begin
            MDAP_SAVE_EMISSION_LINE_SET_OMITTED, kinematics, kinematics_err, intens, intens_err, $
                                                 fluxes, fluxes_err, i, j
        endif else begin
            MDAP_SAVE_EMISSION_LINE_FIT_ENCI_SINGLE, eml_par[j].lambda, voff, result.OIII5007, $
                                                     kinematics, kinematics_err, intens, $
                                                     intens_err, fluxes, fluxes_err, i, j
        endelse
        j = j+1

        omitted[i,j] = result.NII6548[2] gt 0.0 ? 0 : 1
        if omitted[i,j] eq 1 then begin
            MDAP_SAVE_EMISSION_LINE_SET_OMITTED, kinematics, kinematics_err, intens, intens_err, $
                                                 fluxes, fluxes_err, i, j
        endif else begin
            MDAP_SAVE_EMISSION_LINE_FIT_ENCI_SINGLE, eml_par[j].lambda, voff, result.NII6548, $
                                                     kinematics, kinematics_err, intens, $
                                                     intens_err, fluxes, fluxes_err, i, j
        endelse
        j = j+1

        omitted[i,j] = result.Ha6563[2] gt 0.0 ? 0 : 1
        if omitted[i,j] eq 1 then begin
            MDAP_SAVE_EMISSION_LINE_SET_OMITTED, kinematics, kinematics_err, intens, intens_err, $
                                                 fluxes, fluxes_err, i, j
        endif else begin
            MDAP_SAVE_EMISSION_LINE_FIT_ENCI_SINGLE, eml_par[j].lambda, voff, result.Ha6563, $
                                                     kinematics, kinematics_err, intens, $
                                                     intens_err, fluxes, fluxes_err, i, j
        endelse
        j = j+1

        omitted[i,j] = result.NII6583[2] gt 0.0 ? 0 : 1
        if omitted[i,j] eq 1 then begin
            MDAP_SAVE_EMISSION_LINE_SET_OMITTED, kinematics, kinematics_err, intens, intens_err, $
                                                 fluxes, fluxes_err, i, j
        endif else begin
            MDAP_SAVE_EMISSION_LINE_FIT_ENCI_SINGLE, eml_par[j].lambda, voff, result.NII6583, $
                                                     kinematics, kinematics_err, intens, $
                                                     intens_err, fluxes, fluxes_err, i, j
        endelse
        j = j+1

        omitted[i,j] = result.SII6717[2] gt 0.0 ? 0 : 1
        if omitted[i,j] eq 1 then begin
            MDAP_SAVE_EMISSION_LINE_SET_OMITTED, kinematics, kinematics_err, intens, intens_err, $
                                                 fluxes, fluxes_err, i, j
        endif else begin
            MDAP_SAVE_EMISSION_LINE_FIT_ENCI_SINGLE, eml_par[j].lambda, voff, result.SII6717, $
                                                     kinematics, kinematics_err, intens, $
                                                     intens_err, fluxes, fluxes_err, i, j
        endelse
        j = j+1

        omitted[i,j] = result.SII6731[2] gt 0.0 ? 0 : 1
        if omitted[i,j] eq 1 then begin
            MDAP_SAVE_EMISSION_LINE_SET_OMITTED, kinematics, kinematics_err, intens, intens_err, $
                                                 fluxes, fluxes_err, i, j
        endif else begin
            MDAP_SAVE_EMISSION_LINE_FIT_ENCI_SINGLE, eml_par[j].lambda, voff, result.SII6731, $
                                                     kinematics, kinematics_err, intens, $
                                                     intens_err, fluxes, fluxes_err, i, j
        endelse
        j = j+1

END

PRO MDAP_SAVE_EMISSION_LINE_SET_OMITTED, kinematics, kinematics_err, intens, intens_err, fluxes, $
                                         fluxes_err, i, j
        kinematics[i,j,0:1] = 0.0d
        kinematics_err[i,j,0:1] = 1.0d
        intens[i,j] = 0.0d
        intens_err[i,j] = 1.0d
        fluxes[i,j] = 0.0d
        fluxes_err[i,j] = 1.0d
END

; CERTAIN ORDER EXPECTED FOR BELFIORE RESULT
PRO MDAP_SAVE_EMISSION_LINE_FIT_BELFIORE, eml_par, voff, result, omitted, kinematics, $
                                          kinematics_err, intens, intens_err, fluxes, fluxes_err, i

        c=299792.458d                           ; Speed of light in km/s
        neml = n_elements(eml_par)              ; Number of emission lines fitted

        for j=0,neml-1 do begin
            omitted[i,j] = result.Ampl[j] gt 0.0 ? 0 : 1

            if omitted[i,j] eq 1 then begin
                MDAP_SAVE_EMISSION_LINE_SET_OMITTED, kinematics, kinematics_err, intens, $
                                                     intens_err, fluxes, fluxes_err, i, j
                continue
            endif

            kinematics[i,j,0] = result.Vel[j]+voff
            kinematics_err[i,j,0] = result.eVel[j]
            kinematics[i,j,1] = result.Sigma[j]
            kinematics_err[i,j,1] = result.eSigma[j]
            intens[i,j] = result.Ampl[j]
            intens_err[i,j] = result.eAmpl[j]


            fluxes[i,j] = result.Ampl[j] * sqrt(2*!pi) * result.Sigma[j] * eml_par[j].lambda $
                          * (1+(result.Vel[j]+voff)/c)/c

            fluxes_err[i,j] = sqrt( (fluxes[i,j]*result.eAmpl[j]/result.Ampl[j])^2 + $
                                    (fluxes[i,j]*result.eSigma[j]/result.Sigma[j])^2 + $
                                    (result.Ampl[j] * sqrt(2*!pi) * result.Sigma[j] * $
                                     eml_par[j].lambda * result.eVel[j]/c/c)^2 )
        endfor
END


PRO MDAP_EMISSION_LINE_ONLY_FIT, $
                wave, flux, ivar, mask, best_fit_continuum, stellar_kinematics, eml_model, $
                emission_line_kinematics, emission_line_kinematics_err, emission_line_omitted, $
                emission_line_kinematics_individual, emission_line_kinematics_individual_err, $
                emission_line_intens, emission_line_intens_err, emission_line_fluxes, $
                emission_line_fluxes_err, emission_line_EW, emission_line_EW_err, $
                eml_par=eml_par, enci=enci, belfiore=belfiore, quiet=quiet, dbg=dbg

        if keyword_set(enci) and keyword_set(belfiore) then begin
            message, 'Cannot run both Enci and Belfiore fits simultaneously.  Make separate ' $
                     + 'calls to MDAP_EMISSION_LINE_ONLY_FIT'
        endif
        if ~keyword_set(enci) and ~keyword_set(belfiore) and ~keyword_set(quiet) then $
            print, 'No mode selected for MDAP_EMISSION_LINE_ONLY_FIT.  Default is /enci.'

        ; Get the emission line parameters for the lines fit by Enci's
        ; code, in the format of the structure used by GANDALF
        eml_par = MDAP_DEFINE_EMISSION_LINES_ENCI_BELFIORE()

        ; Convert the inverse variance to an error, adjusting the mask for bad values
;       indx = where(ivar gt 0.0d, complement=nindx)
;       if indx[0] eq -1 then $
;           message, 'All pixels are masked!'
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
        c=299792.458d                           ; Speed of light in km/s
        sz = size(flux)
        nobj = sz[1]                            ; Number of object spectra
        npix = sz[2]                            ; Number of pixels in the spectra
        neml = 9                                ; Number of fitted lines

        eml_model = dblarr(nobj,npix)           ; Best-fitting emission-line model

        emission_line_kinematics = dblarr(nobj, 2)
        emission_line_kinematics_err = dblarr(nobj, 2)
        emission_line_omitted = make_array(nobj, neml, /int, value=1)       ; 0/1 -> fit/omitted
        emission_line_kinematics_individual = dblarr(nobj, neml, 2)
        emission_line_kinematics_individual_err = dblarr(nobj, neml, 2)
        emission_line_intens = dblarr(nobj, neml)
        emission_line_intens_err = dblarr(nobj, neml)
        emission_line_fluxes = dblarr(nobj, neml)
        emission_line_fluxes_err = dblarr(nobj, neml)
        emission_line_EW = dblarr(nobj, neml)
        emission_line_EW_err = dblarr(nobj, neml)

        ; TODO: Calculate a chi-square for the full fit using the
        ; eml_model
   
        ; Fit each spectrum and save the results
        if keyword_set(dbg) then $
            nobj = 1
        for i=0,nobj-1 do begin

            ; Set the wavelength vector to the rest wavelength
            ; TODO: Could also just pass velocity to MDAP_FIT_EMISSION_LINE_SPECTRUM_ENCI
            rest_wave = wave/(1+stellar_kinematics[i,0]/c)

            ; Remove the stellar continuum
            galaxy_eml_only = reform(flux[i,*] - best_fit_continuum[i,*])

            if keyword_set(belfiore) then begin
                ; Fit the emission lines
                MDAP_FIT_EMISSION_LINE_SPECTRUM_BELFIORE, rest_wave, galaxy_eml_only, $
                                                          reform(flux_err[i,*]), $
                                                          reform(mask[i,*]), $
                                                          stellar_kinematics[i,1], 0.0d, result, $
                                                          emfit, quiet=quiet

                ; Save the fitted parameters
                MDAP_SAVE_EMISSION_LINE_FIT_BELFIORE, eml_par, stellar_kinematics[i,0], result, $
                                                      emission_line_omitted, $
                                                      emission_line_kinematics_individual, $
                                                      emission_line_kinematics_individual_err, $
                                                      emission_line_intens, $
                                                      emission_line_intens_err, $
                                                      emission_line_fluxes, $
                                                      emission_line_fluxes_err, i

            endif else begin
                ; Fit the emission lines
                MDAP_FIT_EMISSION_LINE_SPECTRUM_ENCI, 0.0d, rest_wave, galaxy_eml_only, $
                                                      reform(flux_err[i,*]), reform(mask[i,*]), $
                                                      result, emfit, /OII3727, /Hb4861, /OIII5007, $
                                                      /Ha6563, /SII6717

                ; Save the fitted parameters
                MDAP_SAVE_EMISSION_LINE_FIT_ENCI, eml_par, stellar_kinematics[i,0], result, $
                                                  emission_line_omitted, $
                                                  emission_line_kinematics_individual, $
                                                  emission_line_kinematics_individual_err, $
                                                  emission_line_intens, emission_line_intens_err, $
                                                  emission_line_fluxes, emission_line_fluxes_err, i
            endelse

            ; Save the best-fitting model
            eml_model[i,*] = emfit
                
            ; Calculate the equivalent width and error
            galaxy_no_eml = reform(flux[i,*]) - emfit
            MDAP_MEASURE_EMISSION_LINE_EQUIV_WIDTH, eml_par, $
                                                    emission_line_kinematics_individual[i,*,0], $
                                                    emission_line_kinematics_individual[i,*,1], $
                                                    rest_wave, galaxy_no_eml, $
                                                    emission_line_fluxes[i,*], $
                                                    emission_line_fluxes_err[i,*], equiv_width, $
                                                    equiv_width_err
            emission_line_EW[i,*] = equiv_width
            emission_line_EW_err[i,*] = equiv_width_err

            ; Decide which lines to inlude in the mean measurement of the kinematics
            indx=MDAP_LINES_FOR_MEAN_GAS_KINEMATICS(emission_line_intens[i,*], $
                                                    emission_line_intens_err[i,*], $
                                                    emission_line_fluxes[i,*], $
                                                    emission_line_fluxes_err[i,*], $
                                                    emission_line_kinematics_individual[i,*,0], $
                                                    emission_line_kinematics_individual_err[i,*,0],$
                                                    emission_line_kinematics_individual[i,*,1], $
                                                    emission_line_kinematics_individual_err[i,*,1],$
                                                    omitted=emission_line_omitted[i,*], count=count)

            ; Get the mean kinematics; use dummy values if no good lines
;           if indx[0] eq -1 then begin
            if count eq 0 then begin
                emission_line_kinematics[i,*] = 0.0d
                emission_line_kinematics_err[i,*] = 99.0d
            endif else begin
                MDAP_MEAN_GAS_KINEMATICS, emission_line_fluxes[i,indx], $
                                          emission_line_kinematics_individual[i,indx,0], $
                                          emission_line_kinematics_individual_err[i,indx,0], $
                                          emission_line_kinematics_individual[i,indx,1], $
                                          emission_line_kinematics_individual_err[i,indx,1], $
                                          fw_vel, fw_vel_err, fw_sig, fw_sig_err
                emission_line_kinematics[i,0] = fw_vel
                emission_line_kinematics_err[i,0] = fw_vel_err
                emission_line_kinematics[i,1] = fw_sig
                emission_line_kinematics_err[i,0] = fw_sig_err
            endelse

        endfor

END



