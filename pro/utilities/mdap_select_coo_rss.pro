;+
; NAME:
;       MDAP_SELECT_COO_RSS
;
; PURPOSE:
;       Determine the on-sky coordinates at a given wavelength for the
;       RSS spectra.
;
; CALLING SEQUENCE:
;       MDAP_SELECT_COO_RSS, wave, skyx, skyy, targwave, bskyx, bskyy
;
; INPUTS:
;       wave dblarr[T]
;               Wavelength vector of the on-sky coordinates.
;
;       skyx dblarr[N][T]
;               X-coordinate for each of the N spectra at each of the T
;               wavelengths.  (XPOS in the DRP LOGRSS files.)
;
;       skyy dblarr[N][T]
;               Y-coordinate for each of the N spectra at each of the T
;               wavelengths.  (YPOS in the DRP LOGRSS files.)
;
;       targwave double or dblarr[N]
;               Wavelength coordinate at which to determine the x and y
;               coordinates.  This can be a single value that is applied
;               to all spectra, or a spectrum-specific value.
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;
; OUTPUT:
;       bskyx dblarr[N]
;               X-coordinate for each of the N spectrum sampled at
;               targwave.
;
;       bskyy dblarr[N]
;               Y-coordinate for each of the N spectrum sampled at
;               targwave.
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
;       04 Dec 2014: (KBW) Original implementation
;-
;------------------------------------------------------------------------------

PRO MDAP_SELECT_COO_RSS, $
                wave, skyx, skyy, targwave, bskyx, bskyy
        sz = size(skyx)
        nx = sz[1]                                          ; Number of spectra
        if n_elements(targwave) eq 1 then begin
            wcoo = make_array(nx, /double, value=targwave)  ; Use the same targwave for all spectra
        endif else $
            wcoo = targwave 
   
        for i=0,nx-1 do begin
            bskyx[i] = interpol(skyx[i,*], wave, wcoo[i])   ; Sample the vectors
            bskyy[i] = interpol(skyy[i,*], wave, wcoo[i])
        endfor
END


