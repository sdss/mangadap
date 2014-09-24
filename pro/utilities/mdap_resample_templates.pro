;+
; NAME:
;	MDAP_RESAMPLE_TEMPLATES
;
; PURPOSE:
;	Resample the template library to match the velocity scale of the object
;	spectra.
;
; CALLING SEQUENCE:
;	MDAP_RESAMPLE_TEMPLATES, tpl_flux, tpl_mask, tpl_wave, tpl_sres, velScale
;
; INPUTS:
;	tpl_flux dblarr[T][S]
;		Array containing T template spectra, each with S spectral
;		channels.  Replaced upon output with the template spectrum
;		sampled with the same interval as the object spectra.
;
;	tpl_ivar dblarr[T][S]
;		Array containing inverse variances in the S spectral channels of
;		the T template spectra.  Replaced upon output with the inverse
;		variances for the template spectra sampled with the same
;		interval as the object spectra.
;
;	tpl_mask dblarr[T][S]
;		Bit mask for template spectra.  Used only to mask spectral
;		regions not covered by all templates.  Value is 0 if pixel is
;		good, value is 1 if it should be masked. Replaced with resampled
;		mask.
;
;	tpl_wave dblarr[T][S]
;		Wavelength of each spectral channel T in angstroms for each
;		template spectrum S.  Replaced by the resampled wavelength
;		coordinates; geometric mean of bin edges.
;
;	tpl_sres dblarr[T][S]
;		Spectral resolution (R=lamda/delta lambda) for of each spectral
;		channel T for each template spectrum S.  Replaced upon output
;		with the resolution resampled to the new wavelength coordinates.
;
;	velScale double
;		Velocity scale for rebinned spectra.
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
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
;	17 Sep 2014: (KBW) Original implementation
;	18 Sep 2014: (KBW) Input velocity scale instead of object wavelength
;			   vector
;	22 Sep 2014: (KBW) Include the inverse variances
;-
;------------------------------------------------------------------------------

PRO MDAP_RESAMPLE_TEMPLATES, $
		tpl_flux, tpl_ivar, tpl_mask, tpl_wave, tpl_sres, velScale

	sz=size(tpl_flux)
	nt = sz[1]			; Number of template spectra

	ns = 0				; New maximum size of template spectra
	lamRange=dblarr(nt,2)		; Wavelength range of template spectra
	for i=0,nt-1 do begin
	    gtpl = where(tpl_mask[i,*] lt 1)
	    ng = n_elements(gtpl)	; Number of unmasked pixels
	    if ns lt ng then $		; If larger than current maximum, save it
		ns = ng

	    lamRange[i,0] = min(tpl_wave[i,gtpl])	; Initial valid wavelength
	    lamRange[i,1] = max(tpl_wave[i,gtpl])	; Final valid wavelength
	    
	endfor

	print, 'Estimate ns: ', ns

	tpl_flux_ = tpl_flux		; Save original spectra, masks, and wavelengths
	tpl_ivar_ = tpl_flux		; Save original spectra, masks, and wavelengths
	tpl_mask_ = tpl_mask
	tpl_wave_ = tpl_wave

	tpl_flux=dblarr(nt, ns)		; Re-initialize the arrays
	tpl_ivar=dblarr(nt, ns)
	tpl_wave=dblarr(nt, ns)
	tpl_mask=dblarr(nt, ns)

	nsr = 0
	for i=0,nt-1 do begin
	    print, 'Rebinning spectrum: ', i+1

	    ; TODO: Must only select a single, contiguous spectral region
	    ;		- What happens if the are masked regions spread through the spectral range?
	    gtpl = where(tpl_mask_[i,*] lt 1 and tpl_ivar_[i,*] gt 0.)	; Select unmasked pixels

	    ; Rebin both the spectrum and the mask
	    ;    Use reform() to pass a 1D vector, as expected by MDAP_DO_LOG_REBIN
	    MDAP_DO_LOG_REBIN, lamRange[i,*], reform(tpl_flux_[i,gtpl]), tpl_rebin, logLam, $
			       velscale=velScale, /log10
	    MDAP_DO_LOG_REBIN, lamRange[i,*], reform(1.0/tpl_ivar_[i,gtpl]), tpl_ivr_rebin, $
			       logLam, velscale=velScale, /log10
	    MDAP_DO_LOG_REBIN, lamRange[i,*], reform(tpl_mask_[i,gtpl]), tpl_msk_rebin, logLam, $
			       velscale=velScale, /log10

	    ng=n_elements(tpl_rebin)			; Number of rebinned pixels
	    print, 'Rebinned length: ', ng
	    tpl_flux[i,0:ng-1] = tpl_rebin[0:ng-1]		; Save the spectrum
	    tpl_ivar[i,0:ng-1] = 1.0/tpl_ivr_rebin[0:ng-1]	; Save the variances
	    tpl_wave[i,0:ng-1] = 10^(logLam[0:ng-1])		; Save the wavelengths in angstroms

	    gthalf = where(tpl_msk_rebin gt 0.5)	; Find rebinned mask pixels with >0.5
	    if gthalf ne -1 then $
		tpl_mask[i,gthalf]=1.0			;  ... and set these to masked pixels
	    if ng lt ns then $
		tpl_mask[i,ng:ns-1]=1.0			; Mask all pixels without data

	    if nsr lt ng then $				; Save the maximum size of the spectra
		nsr = ng
	endfor

	print, 'Final ns: ', nsr

	; Remove spectral regions that are masked in ALL templates
	tpl_flux = temporary(tpl_flux[*,0:nsr-1])
	tpl_ivar = temporary(tpl_ivar[*,0:nsr-1])
	tpl_mask = temporary(tpl_mask[*,0:nsr-1])
	tpl_wave = temporary(tpl_wave[*,0:nsr-1])

	; Interpolate the spectral resolution to the new wavelength coordinates
	tpl_sres_ = tpl_sres
	tpl_sres = dblarr(nt,nsr)
	for i=0,nt-1 do $
	    tpl_sres[i,*] = (interpol(tpl_sres_[i,*], tpl_wave_[i,*], tpl_wave[i,*]))[*]

END

	    

