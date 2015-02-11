;+
; NAME:
;       MDAP_RESAMPLE_TEMPLATES
;
; PURPOSE:
;       Resample the template library to match the velocity scale of the object
;       spectra.
;
; CALLING SEQUENCE:
;       MDAP_RESAMPLE_TEMPLATES, tpl_flux, tpl_mask, tpl_wave, tpl_sres, velScale, /reform_sres
;
; INPUTS:
;       tpl_flux dblarr[T][S]
;               Array containing T template spectra, each with S spectral
;               channels.  Replaced upon output with the template spectrum
;               sampled with the same interval as the object spectra.
;
;       tpl_ivar dblarr[T][S]
;               Array containing inverse variances in the S spectral channels of
;               the T template spectra.  Replaced upon output with the inverse
;               variances for the template spectra sampled with the same
;               interval as the object spectra.
;
;       tpl_mask dblarr[T][S]
;               Bit mask for template spectra.  Used only to mask spectral
;               regions not covered by all templates.  Value is 0 if pixel is
;               good, value is 1 if it should be masked. Replaced with resampled
;               mask.
;
;       tpl_wave dblarr[T][S]
;               Wavelength of each spectral channel T in angstroms for each
;               template spectrum S.
;
;               ** Replaced on output by a single vector (dblarr[S]) for the
;               wavelength of the S spectral channels for ALL templates;
;               geometric mean of bin edges.
;
;       tpl_sres dblarr[T][S]
;               Spectral resolution (R=lamda/delta lambda) for of each spectral
;               channel T for each template spectrum S.  Replaced upon output
;               with the resolution resampled to the new wavelength coordinates.
;
;               If /reform_sres is set, output is a single vector restricted to
;               the common wavelength range.
;
;       velScale double
;               Velocity scale for rebinned spectra.
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;       /reform_sres
;               Assumes that the sres vector should be the same for all spectra,
;               the only difference being that the template spectra may cover
;               different wavelength ranges.
;
; OUTPUT:
;
; OPTIONAL OUTPUT:
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
;       17 Sep 2014: Original implementation by K. Westfall (KBW)
;       18 Sep 2014: (KBW) Input velocity scale instead of object wavelength
;                          vector
;       22 Sep 2014: (KBW) Include the inverse variances
;       24 Sep 2014: (KBW) Force the templates to be registered to the same
;                          wavelength solution.
;       13 Oct 2014: (KBW) Reform the tpl_sres vector to a single vector
;                          applicable to all template spectra
;       05 Jan 2015: (KBW) Allow for some templates to be fully masked
;       10 Feb 2015: (KBW) Fixed bug when resampled templates were
;                          longer than originals.  Done by using the
;                          velscale and MDAP_LOG_WAVELENGTH_SCALE() to
;                          get the number of pixels needed to cover the
;                          provided range.
;-
;------------------------------------------------------------------------------

PRO MDAP_RESAMPLE_TEMPLATES, $
                tpl_flux, tpl_ivar, tpl_mask, tpl_wave, tpl_sres, velScale, log10=log10, $
                reform_sres=reform_sres

        sz=size(tpl_flux)
        nt = sz[1]                      ; Number of template spectra

        dlogl = MDAP_LOG_WAVELENGTH_SCALE(velScale, log10=log10)

        common_wave=dblarr(2)           ; Wavelength range common to ALL template spectra

        ; Find the first spectrum that is not fully masked
        i = 0
        gtpl = where(tpl_mask[i,*] lt 1, count)
        while count eq 0 && i lt nt do begin
            i = i+1
            gtpl = where(tpl_mask[i,*] lt 1, count)
        endwhile

        ; Could not find one!
        if count eq 0 then $
            message, 'ALL templates fully masked!'

        ; Initialize the common wavelength range to the limits of the
        ; good pixels from the first good template spectrum.
        common_wave[0] = min(tpl_wave[i,gtpl])
        common_wave[1] = max(tpl_wave[i,gtpl])

        ; Determine the wavelength range common to ALL templates
        lamRange=dblarr(nt,2)           ; Wavelength range of each template spectrum
        for i=0,nt-1 do begin
            gtpl = where(tpl_mask[i,*] lt 1, ng)
            if ng eq 0 then $
                continue
            minw = min(tpl_wave[i,gtpl])
            maxw = max(tpl_wave[i,gtpl])

            lamRange[i,0] = minw                ; Initial valid wavelength
            lamRange[i,1] = maxw                ; Final valid wavelength

            if common_wave[0] lt minw then $    ; Common wavelength range
                common_wave[0] = minw
            if common_wave[1] gt maxw then $
                common_wave[1] = maxw
        endfor

;       print, lamRange
;       print, common_wave

        ; Roughly the number of pixels required to cover this spectral range logarithmically
        if keyword_set(log10) then begin
            ns = ceil((alog10(maxw)-alog10(minw))/dlogl)
        endif else $
            ns = ceil((alog(maxw)-alog(minw))/dlogl)
        print, 'Estimate ns: ', ns

        tpl_flux_ = tpl_flux            ; Save original spectra, masks, and wavelengths
        tpl_ivar_ = tpl_ivar
        tpl_mask_ = tpl_mask
        tpl_wave_ = tpl_wave

        tpl_flux=dblarr(nt, ns)         ; Re-initialize the arrays
        tpl_ivar=dblarr(nt, ns)
        logLam=dblarr(ns)               ; A SINGLE vector for the wavelength solution
        tpl_mask=dblarr(nt, ns)

        nsr = 0
        for i=0,nt-1 do begin
            print, 'Rebinning spectrum: ', i+1

            ; TODO: Must only select a single, contiguous spectral region
            ;           - What happens if the are masked regions spread through the spectral range?
            ; Select unmasked pixels
            gtpl = where(tpl_mask_[i,*] lt 1 and tpl_ivar_[i,*] gt 0., count)
            if count eq 0 then $
                continue

            ; Rebin both the spectrum and the mask
            ;    Use reform() to pass a 1D vector, as expected by MDAP_DO_LOG_REBIN

            MDAP_DO_LOG_REBIN, lamRange[i,*], reform(tpl_flux_[i,gtpl]), tpl_rebin, logLam, $
                               velscale=velScale, /log10, newrange=common_wave
;           ng=n_elements(tpl_rebin)                    ; Number of rebinned pixels
;           print, ng
            MDAP_DO_LOG_REBIN, lamRange[i,*], reform(1.0/tpl_ivar_[i,gtpl]), tpl_ivr_rebin, $
                               logLam, velscale=velScale, /log10, newrange=common_wave
;           ng=n_elements(tpl_ivr_rebin)                ; Number of rebinned pixels
;           print, ng
            MDAP_DO_LOG_REBIN, lamRange[i,*], reform(tpl_mask_[i,gtpl]), tpl_msk_rebin, logLam, $
                               velscale=velScale, /log10, newrange=common_wave
;           ng=n_elements(tpl_msk_rebin)                ; Number of rebinned pixels
;           print, n_elements(logLam)
;           print, ng

            ng=n_elements(tpl_rebin)                    ; Number of rebinned pixels
            if ng gt ns then $
                message, 'tpl_flux not big enough'
            print, 'Rebinned length: ', ng
            tpl_flux[i,0:ng-1] = tpl_rebin[0:ng-1]              ; Save the spectrum
            tpl_ivar[i,0:ng-1] = 1.0/tpl_ivr_rebin[0:ng-1]      ; Save the variances
;           tpl_wave[i,0:ng-1] = 10^(logLam[0:ng-1])            ; Save the wavelengths in angstroms

            gthalf = where(tpl_msk_rebin gt 0.5, count)        ; Find rebinned mask pixels with >0.5
;           if gthalf ne -1 then $
            if count ne 0 then $
                tpl_mask[i,gthalf]=1.0                  ;  ... and set these to masked pixels
            if ng lt ns then $
                tpl_mask[i,ng:ns-1]=1.0                 ; Mask all pixels without data

            if nsr lt ng then $                         ; Save the maximum size of the spectra
                nsr = ng                                ; Should be the same for all spectra
        endfor

        print, 'Final ns: ', nsr

        ; Remove spectral regions that are masked in ALL templates
        tpl_flux = temporary(tpl_flux[*,0:nsr-1])
        tpl_ivar = temporary(tpl_ivar[*,0:nsr-1])
        tpl_mask = temporary(tpl_mask[*,0:nsr-1])
        tpl_wave = 10^(logLam[0:nsr-1])                 ; Same for all spectra, use last instance
        
        ; Interpolate the spectral resolution to the new wavelength coordinates
        tpl_sres_ = tpl_sres
        if keyword_set(reform_sres) then begin
            tpl_sres = dblarr(nsr)
            tpl_sres[*] = interpol(tpl_sres_[0,*], tpl_wave_[0,*], tpl_wave[*])
        endif else begin
            tpl_sres = dblarr(nt,nsr)
            for i=0,nt-1 do $
                tpl_sres[i,*] = interpol(tpl_sres_[i,*], tpl_wave_[i,*], tpl_wave[*])
        endelse

END

            

