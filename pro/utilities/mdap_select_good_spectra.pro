;+
; NAME:
;       MDAP_SELECT_GOOD_SPECTRA
;
; PURPOSE:
;       Create a list of indices with spectra that satisfy the following
;       criteria for at least 20% of the length of the spectrum:
;
;               - pixel is unmasked from DRP file
;               - ivar is larger than 0
;               - ivar is not inf or NaN
;               - flux is not inf or NaN
;       
;       and that
;
;               - min(flux) ne max(flux)
;
; CALLING SEQUENCE:
;       MDAP_SELECT_GOOD_SPECTRA, flux, ivar, mask, gflag, gindx
;
; INPUTS:
;       flux dblarr[N][T]
;               Flux of N spectra read from a DRP-produced fits file.
;
;       ivar dblarr[N][T]
;               Inverse variance of spectra.
;
;       mask dblarr[N][T]
;               Bad pixel mask (0-unmasked; 1-masked)
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;
; OUTPUT:
;       gflag intarr[N]
;               Flag that spectrum is good (0=false; 1=true)
;
;       gindx intarr[]
;               Indices of good spectra.
;
; OPTIONAL OUTPUT:
;       count integer
;               Number of good spectra
;
; COMMENTS:
;
; EXAMPLES:
;
; TODO:
;       - Allow the criteria (particularly the 20% one) to change?
;
; BUGS:
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;       09 Sep 2014: Original implementation by K. Westfall (KBW)
;       05 Jan 2015: (KBW) Bug in size of indx check; return count
;       06 Feb 2015: (KBW) Incorporate pixel mask (based on DRP MASK
;                          extension, see MDAP_CONVERT_DRP_MASK)
;-
;------------------------------------------------------------------------------

PRO MDAP_SELECT_GOOD_SPECTRA, $
                flux, ivar, mask, gflag, gindx, count=count

        sz=size(flux)
        ns=sz[1]                                        ; Number of spectra
        nc=sz[1]                                        ; Number of spectral channels

        gflag=intarr(ns)                                ; Initialize all flags as false (gflag=0)
        for i=0,ns-1 do begin

            ; Select bad pixels according to the criteria listed above
            
            indx=where(mask[i,*] gt 0.5 or ivar[i,*] le 0 $
                       or finite(ivar[i,*]) eq 0 or finite(flux[i,*]) eq 0, count)

            ; Spectrum is either made up of 20% or more bad pixels, or it is constant (not real)
;           if n_elements(indici) ge 0.2*nc or min(flux[i,*]) eq max(flux[i,*]) then begin
;           if count ge 0.2*nc or min(flux[i,*]) eq max(flux[i,*]) then begin
            if count ge 0.2*nc || min(flux[i,*]) eq max(flux[i,*]) then begin
                gflag[i] = 0                            ; Flag spectrum as bad
            endif else $
                gflag[i] = 1                            ; Flag spectrum as good
        endfor

        gindx=where(gflag eq 1, count)                  ; List the indices of good spectra
        print, 'Number of good DRP spectra: '+MDAP_STC(count, /integer)+'/'+MDAP_STC(ns, /integer)
END
                

