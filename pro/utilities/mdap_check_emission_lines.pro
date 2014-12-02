;+
; NAME:
;       MDAP_CHECK_EMISSION_LINES
;
; PURPOSE:
;       Run some checks on the emission line data.  Ignore any lines that are
;       outside the wavelength range, defined by the input vector.
;
; CALLING SEQUENCE:
;       MDAP_CHECK_EMISSION_LINES, eml_par, wave, vel=vel
;
; INPUTS:
;       eml_par EmissionLine[E]
;               The parameters for each of E emission lines used during the
;               fitting procedures.  The EmissionLine structure is defined as
;               follows (see MDAP_READ_EMISSION_LINE_PARAMETERS):
;
;               { EmissionLine, i:0L, name:'', lambda:0.0d, action:'', $
;                     kind:'', a:0.0d, v:0.0d, s:0.0d, fit:'' }
;
;               Once created, one selects, for example, the name of the 3rd
;               input line using: eml_par[2].name
;
;       wave dblarr[C]
;               The wavelength in angstroms for each of C spectral channels in
;               the object spectrum.
;
; OPTIONAL INPUTS:
;       vel double
;               Velocity offset of the input wavelength vector from the
;               reference frame of the defined emission lines.  'wave' is
;               assumed to be in rest wavelengths if vel is not defined.
;
; OPTIONAL KEYWORDS:
;
; OUTPUT:
;       eml_par is updated on output
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
;       30 Sep 2014: (KBW) Original implementation
;-
;------------------------------------------------------------------------------

PRO MDAP_CHECK_EMISSION_LINES, $
        eml_par, wave, velocity=velocity, quiet=quiet

        if n_elements(velocity) eq 0 then $
            velocity=0.0d                       ; Set velocity to zero if not provided

        c = 299792.458d                         ; Speed of light in km/s

        ; Find where line centroids are outside of the input wavelength range
        indx=where( eml_par.lambda lt wave[0]/(1.0d + velocity/c) or $
                    eml_par.lambda gt wave[n_elements(wave)-1]/(1.0d + velocity/c) )

        ; TODO: Report when a line is ignored?
        if indx[0] ne -1 then begin
            if ~keyword_set(quiet) then $
                print, 'These lines fall outside the spectral range:', eml_par[indx].i
            eml_par[indx].action='i'            ; Ignore the line
        endif
END

