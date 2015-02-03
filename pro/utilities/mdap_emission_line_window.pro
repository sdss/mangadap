;+
; NAME:
;       MDAP_EMISSION_LINE_WINDOW
;
; PURPOSE:
;       Use the input to define the starting and ending wavelengths of a
;       window about a redshifted emission line.
;
; CALLING SEQUENCE:
;       result = MDAP_EMISSION_LINE_WINDOW(eml_par, velocity=velocity, nsig=nsig, sigma=sigma)
;
; INPUTS:
;       eml_par EmissionLine
;               The parameters for a single emission line.  The
;               EmissionLine structure is defined as follows (see
;               MDAP_ASSIGN_EMISSION_LINE_PARAMETERS):
;
;               { EmissionLine, i:0L, name:'', lambda:0.0d, action:'', $
;                               kind:'', a:0.0d, v:0.0d, s:0.0d, fit:'' }
;
; OPTIONAL INPUTS:
;       velocity double
;               Expected redshift (cz) of the line in km/s.  Default is
;               to use eml_par.v.
;
;       nsig double
;               Number of sigma for the window. Default is 3.0.
;
;       sigma double
;               Expected velocity dispersion of the line in km/s.
;               Default is to use eml_par.s.
;
; OPTIONAL KEYWORDS:
;
; OUTPUT:
;       Output is a dblarr[2] with the start and end of the wavelength
;       window.
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

FUNCTION MDAP_EMISSION_LINE_WINDOW, $
                eml_par, velocity=velocity, nsig=nsig, sigma=sigma

        if n_elements(nsig) eq 0 then $
            nsig=3.0d                                   ; Default window size in +/- sigma

        if n_elements(sigma) eq 0 then $
            sigma=eml_par.s                             ; Set width based on input value

        ; Determine the edges of the window
        el_half_width = nsig*sigma                      ; The half-width of the mask in km/s
        sign = make_array(2, /double, value=1.0d)
        sign[0] = -1.0d                                 ; Sign for lower edge
        c = 299792.458d                                         ; Speed of light in km/s
        el_range = (1 + sign*el_half_width/c)*eml_par.lambda    ; Lower and upper edge (angstroms)

        ; Redshift the edges to match the object wavelengths, unless it's a sky line!
        if n_elements(velocity) ne 0 and eml_par.action ne 's' then $
                el_range = el_range * (1+velocity/c)            ; Redshifted wavelengths

        return, el_range                                ; Return the window range
END


