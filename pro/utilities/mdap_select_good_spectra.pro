;+
; NAME:
;       MDAP_SELECT_GOOD_SPECTRA
;
; PURPOSE:
;       Create a list of indices with spectra that satisfy the following
;       criteria for some threshold percentage of the pixels:
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
;       MDAP_SELECT_GOOD_SPECTRA, flux, ivar, mask, gflag, gindx, count=count, goodfrac=goodfrac, $
;                                 /quiet
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
;       good_fraction_threshold floating-point
;               Threshold fraction of good pixels in a spectrum before
;               ignoring the spectrum.  Default is 0.2; i.e., any
;               spectrum with a fraction of good pixels greater than 0.2
;               will be considered a good spectrum to analyze.
;
; OPTIONAL KEYWORDS:
;       /quiet
;               Surpress output to the screen
;
; OUTPUT:
;       gflag intarr[N]
;               Flag that spectrum is good (0=false; 1=true)
;
;       gindx intarr[]
;               Indices of good spectra.
;
; OPTIONAL OUTPUT:
;       fraction_good dblarr[N]
;               Fraction of good pixels in each spectrum.
;
;       min_eq_max intarr[N]
;               Flag that the minimum flux equals the maximum flux,
;               presumably because the spectrum has no data.
;
;       count integer
;               Number of good spectra

;
; COMMENTS:
;
; EXAMPLES:
;
; TODO:
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
;       09 Feb 2015: (KBW) Debug problem with incorporation of pixel
;                          mask. Add fraction_good and min_eq_max to
;                          output.
;-
;------------------------------------------------------------------------------

PRO MDAP_SELECT_GOOD_SPECTRA, $
                flux, ivar, mask, gflag, gindx, good_fraction_threshold=good_fraction_threshold, $
                fraction_good=fraction_good, min_eq_max=min_eq_max, count=count, quiet=quiet

        if n_elements(goodfrac) eq 0 then $
            goodfrac = 0.2

        sz=size(flux)
        ns=sz[1]                                        ; Number of spectra
        nc=sz[2]                                        ; Number of spectral channels

        if ~keyword_set(quiet) then begin
            print, ['Spec', 'max=min', 'Unmasked', 'Good ivar', 'Fin. ivar', 'Fin. flux', $
                    'Total'], format='(%"%9s %9s %9s %9s %9s %9s %9s")'
        endif

        gflag=intarr(ns)                                ; Initialize all flags as false (gflag=0)
        fraction_good=dblarr(ns)                        ; Fraction of good pixels in the spectrum
        min_eq_max=intarr(ns)                           ; Flag that min(flux) eq max(flux)
        for i=0,ns-1 do begin

;           if min(flux[i,*]) eq max(flux[i,*]) then begin
;               print, [i, -1., -1., -1., -1., -1.], $
;                      format='(%"%9d %9.2f %9.2f %9.2f %9.2f %9.2f")'
;                if ~keyword_set(quiet) then begin
;                    print, [MDAP_STC(i,/integer), 'INVALID SPECTRUM (max=min)'], $
;                           format='(%"%9s %s")'
;                endif
;               continue
;           endif

            ; If reporting, determine which pixels fell into which category
            if ~keyword_set(quiet) then begin
                indx=where(mask[i,*] lt 1., count)          ; Unmasked pixels
                unmasked = double(count)*100./double(nc)
                
                indx=where(ivar[i,*] gt 0., count)          ; Positive ivar values
                posivar = double(count)*100./double(nc)
                
                indx=where(finite(ivar[i,*]) eq 1, count)   ; Finite ivar values
                finivar = double(count)*100./double(nc)

                indx=where(finite(flux[i,*]) eq 1, count)   ; Finite flux values
                finflux = double(count)*100./double(nc)
            endif

            ; Select good pixels according to the criteria listed above
            indx=where(mask[i,*] lt 1. and ivar[i,*] gt 0. and finite(ivar[i,*]) eq 1 $
                       and finite(flux[i,*]) eq 1, count)

            fraction_good[i] = double(count)/double(nc)         ; Total number of good pixels
            min_eq_max[i] = min(flux[i,*]) eq max(flux[i,*])    ; Flag the min(flux) eq max(flux)
    
            ; Spectrum is good if the total number of good pixels is
            ; greater or equal to the threshold (goodfrac)
            if fraction_good[i] ge good_fraction_threshold then $
                gflag[i] = 1                            ; Flag spectrum as good

            ; Make sure any spectrum with min(flux) eq max(flux) is also
            ; flagged as bad
            if min_eq_max[i] && gflag[i] eq 1 then begin
                print, 'WARNING: Spectrum with max(flux)=min(flux) passed.  Marking as bad!'
                gflag[i] = 0
            endif 

            ; Report on this spectrum
            if ~keyword_set(quiet) then begin
                print, [i, min_eq_max[i], unmasked, posivar, finivar, finflux, $
                        fraction_good[i]*100.], format='(%"%9d %9d %9.2f %9.2f %9.2f %9.2f %9.2f")'
            endif

        endfor

        gindx=where(gflag eq 1, count)                  ; List the indices of good spectra
        print, 'Number of good DRP spectra: ' + MDAP_STC(count, /integer) + '/' $
               + MDAP_STC(ns, /integer) + " = " + MDAP_STC(double(count)/double(ns))
END
                

