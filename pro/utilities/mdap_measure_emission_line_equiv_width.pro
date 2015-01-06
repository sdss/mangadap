;+
; NAME:
;       MDAP_MEASURE_EMISSION_LINE_EQUIV_WIDTH
;
; PURPOSE:
;       Measure the equivalent width of a fitted emission line.  Note
;       that the definition here is DIFFERENT than the equivalent widths
;       of the spectral indices.  The two have opposite sign and the
;       details of the calculations are different.
;
;       The continuum measurement is set to be between 5-10 sigma on
;       either side of the fitted emission line.
;
;       TODO: Should check with emission-line experts to apply the
;       conventional approach.
;
; CALLING SEQUENCE:
;       MDAP_MEASURE_EMISSION_LINE_EQUIV_WIDTH, eml_par, velocity, sigma, wave, galaxy_no_eml, $
;                                               flux, flux_err, equiv_width, equiv_width_err
;
; INPUTS:
;       eml_par EmissionLine[E]
;               The parameters for the E fitted emission lines.  The
;               EmissionLine structure is defined as follows (see
;               MDAP_ASSIGN_EMISSION_LINE_PARAMETERS):
;
;               { EmissionLine, i:0L, name:'', lambda:0.0d, action:'', $
;                               kind:'', a:0.0d, v:0.0d, s:0.0d, fit:'' }
;
;       velocity dblarr[E]
;               Redshift (cz) of the E emission lines in km/s.
;
;       sigma dblarr[E]
;               Velocity dispersion of the E emission lines in km/s.
;
;       wave dblarr[C]
;               Wavelength of all C spectral channels, which is the same for all
;               N spectra.
;
;       galaxy_no_eml dblarr[C]
;               Galaxy spectrum with the emission-lines subtracted.
;
;       flux dblarr[E]
;       flux_err dblarr[E]
;               Gaussian line flux and error for the E fitted emission
;               lines in the galaxy spectrum. 
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;
; OUTPUT:
;       equiv_width dblarr[E]
;       equiv_width_err dblarr[E]
;               Equivalent width and its error for the E emission lines.
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
;       10 Dec 2014: (KBW) Pulled from MDAP_GANDALF_WRAP for use with
;                          other procedures.
;-
;------------------------------------------------------------------------------

PRO MDAP_MEASURE_EMISSION_LINE_EQUIV_WIDTH, $
                eml_par, velocity, sigma, wave, galaxy_no_eml, flux, flux_err, equiv_width, $
                equiv_width_err

        ; TODO: Perform a more robust measurement of the EW!
        ; TODO: Only fit lines within 10 sigma of the edge (mask if part of it is there)

        neml = n_elements(eml_par)
        equiv_width = dblarr(neml, /nozero)
        equiv_width_err  = dblarr(neml, /nozero)

        for i=0,neml-1 do begin
            el_range_in = MDAP_EMISSION_LINE_WINDOW(eml_par[i], velocity=velocity, nsig=5.0d, $
                                                    sigma=sigma)
            el_range_out = MDAP_EMISSION_LINE_WINDOW(eml_par[i], velocity=velocity, nsig=10.0d, $
                                                     sigma=sigma)
            el_range_lo = [ el_range_out[0], el_range_in[0] ]
            el_range_hi = [ el_range_in[1], el_range_out[1] ]

            ; TODO: Do these tests outside of the loop?  Do something
            ; else so it doesn't break the DAP 
            MDAP_SELECT_WAVE, wave, el_range_lo, indx_lo, count=count_lo
            MDAP_SELECT_WAVE, wave, el_range_hi, indx_hi, count=count_hi
;           if indx_lo[0] eq -1 or indx_hi[0] eq -1 then begin
            if count_lo eq 0 or count_hi eq 0 then begin
                ; Set dummy values
                equiv_width[i] = 0.0d
                equiv_width_err[i] = 99.0d
                continue
            endif

            indx = [ indx_lo, indx_hi ]

            cont = median(galaxy_no_eml[indx])                  ; Median continuum
            econt = robust_sigma(galaxy_no_eml[indx])           ; TODO: Overestimate?
            
            equiv_width[i] = flux[i]/cont                       ; Approximation
            equiv_width_err[i] = sqrt( (flux_err[i]/cont)^2 + (equiv_width[i]*econt/cont)^2 )
        endfor
END


