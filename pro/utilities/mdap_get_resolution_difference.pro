;+
; NAME:
;       MDAP_GET_RESOLUTION_DIFFERENCE
;
; PURPOSE:
;       Determine the difference between the resolution of the template library
;       and the galaxy data.  The procedure assumes the FWHM values are of a
;       Gaussian profile!
;
; CALLING SEQUENCE:
;       MDAP_GET_RESOLUTION_DIFFERENCE, fwhm_instr, library_key, res_var_diff
;
; INPUTS:
;       fwhm_instr dblarr[C]
;               Instrumental resolution (FHWM) in angstroms for each spectral
;               channel of the galaxy data.
;
;       library_key string
;               Keyword used to set the library being used to fit the data.
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;
; OUTPUT:
;       res_var_diff dblarr[C]
;               Quadrature difference between the galaxy and template library
;               resolution (sigma, now FWHM) in angstroms^2 for each spectral
;               channel in the galaxy data.  For a Gaussian instrumental
;               profile, this is the wavelength-dependent variance of the
;               Gaussian needed to match the template library resolution to the
;               galaxy resolution.  This is left as angstroms^2 to allow for
;               negative values (i.e. a better galaxy resolution than template
;               resolution).
;
; OPTIONAL OUTPUT:
;
; COMMENTS:
;
; EXAMPLES:
;
; TODO:
;       - Only accepts the keyword 'MARCS' so far with a hard-wired, wavelength
;         independent resolution.  This procedure will need to be expanded
;         eventually to accommodate more libraries.
;       - Assumes a GAUSSIAN profile!  Should this be further generalized?
;
; BUGS:
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;       17 Sep 2014: (KBW) Original implementation
;-
;------------------------------------------------------------------------------

PRO MDAP_GET_RESOLUTION_DIFFERENCE, $
                fwhm_instr, library_key, res_var_diff

        sz=size(fwhm_instr)
        nc=sz[1]                                        ; Number of spectral channels
        fwhm_lib=dblarr(nc, /nozero)                    ; Initialize the FWHM of the library
        if library_key eq 'MARCS' then begin
            fwhm_lib[*] = 2.73                          ; TODO: Hard-wired value.  Should somehow
                                                        ;       get this from the library itself
        endif else begin
            message, 'Unknown library keyword!'
            retall
        endelse

        fwhm2sig = 1./(2.*sqrt(alog(4.)))               ; FWHM to Gaussian sigma conversion factor

        res_var_diff = (fwhm_instr*fwhm2sig)^2 - (fwhm_lib*fwhm2sig)^2  ; Diff. in Gaussian variance

END     

