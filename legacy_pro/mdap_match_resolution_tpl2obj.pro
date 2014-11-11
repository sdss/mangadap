;+
; NAME:
;	MDAP_MATCH_RESOLUTION_TPL2OBJ
;
; PURPOSE:
;	Match the spectral resolution of a set of templates to the spectral
;	resolution of the object data, allowing for a wavelength-dependent
;	resolution difference.
;
;	In the case where the template resolution is LOWER than the object
;	resolution, there are three choices involved: (1) match the template
;	resolution to the object resolution up to some constant offset that must
;	be accounted for in subsequent analyses, (2) trim the template spectral
;	range to only those spectral regions where the resolution is better than
;	the galaxy resolution, or (3) allow for a wavelength dependent
;	difference in the spectral resolution that must be accounted for in
;	subsequent analyses.  Option 1 is the default behavior; select option 2
;	using the keyword /no_offset; and option 3 is not currently allowed.
;
;	If using option 1, the keyword /variable_offset allows the offset to be
;	different for all templates.  Under the expectation that the templates
;	will be combined and tested against the object spectra, the default
;	behavior is to force the offset to be the same for all templates.
;
;	Finally, the code masks a number of pixels at the beginning and end of
;	the template spectra to remove regions affected by errors in the
;	convolution due to the censoring of the data.  The number of pixels is
;	the FWHM of the largest Gaussian applied in the convolution:
;	ceil(sig2fwhm*max(diff_sig_w)/dw).  This is currently hard-wired and
;	should be tested.
;
;	The detailed match of the galaxy and template spectral resolution should
;	account for the redshift of the galaxy spectrum with respect to the
;	template.  This is not accounted for here, assuming that the wavelength
;	vectors of both the template and galaxy are *in the same reference
;	frame*!
;
; CALLING SEQUENCE:
;	MDAP_MATCH_RESOLUTION_TPL2OBJ, tpl_flux, tpl_ivar, tpl_mask, tpl_wave, $
;				       tpl_sres, wave, sres, tpl_soff, /no_offset, $
;				       /variable_offset
;
; INPUTS:
;	tpl_flux dblarr[T][S]
;		Array containing T template spectra, each with S spectral
;		channels.  Replaced upon output with the template spectrum with
;		a resolution matched (as best as possible) to the object
;		spectra.
;
;	tpl_ivar dblarr[T][S]
;		Array containing inverse variance of the S spectral channels in
;		each of the T template spectra.  Replaced upon output with the
;		inverse variances at the resolution matched (as best as
;		possible) to the object spectra.
;
;	tpl_mask dblarr[T][S]
;		Bit mask for template spectra.  Used only to mask spectral
;		regions not covered by all templates.  Value is 0 if pixel is
;		good, value is 1 if it should be masked.
;
;	tpl_wave dblarr[T][S]
;		Wavelength of each spectral channel T in angstroms for each
;		template spectrum S.  TODO: Assumed to be linear!  Should be in
;		the same reference frame as the galaxy data!
;
;	tpl_sres dblarr[T][S]
;		Spectral resolution (R=lamda/delta lambda) for of each spectral
;		channel T for each template spectrum S.  Replaced upon output
;		with the matched resolution of the template spectra.
;
;	wave dblarr[C]
;		Wavelength coordinates of each spectral channel C in the object
;		spectrum, in accordance with the DRP output.  That is, the
;		vector is expected to have a constant step in log10(lambda);
;		however, the coordinates are in angstroms, not log10(angstroms).
;		Should be in the same rerference frame as the template data!
;
;	sres dblarr[C]
;		Median spectral resolution (R=lamda/delta lamba) as a function
;		of wavelength for all fibers.
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;
;	/no_offset
;		Do not allow a baseline offset of the template resolution.
;		Instead, trim the template spectra to the wavelength range where
;		this offset is zero.
;
;	/variable_offset
;		Allow the offset to be different for each template, instead of
;		forcing the offset to be the same for all templates.
;
; OUTPUT:
;	tpl_soff dblarr[T]
;		A constant dispersion artificially added to the spectral
;		resolution of each template T to ensure that the convolution
;		used to match the spectral resolution is always robust.  This,
;		by definition, means that tpl_soff > 0 when the object spectral
;		resolution is better than than the template resolution.
;
; OPTIONAL OUTPUT:
;
; COMMENTS:
;
; EXAMPLES:
;
; TODO:
;	- MDAP_CONVOL_SIGMA is very slow.  Farm this out to a C/C++ program?
;
; BUGS:
;
; PROCEDURES CALLED:
;	MDAP_MINIMUM_CNVLV_SIGMA
;	MDAP_CONVOL_SIGMA
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;	17 Sep 2014: (KBW) Original Implementation
;	23 Oct 2014: (KBW) Documentation specifies that the wavelengths should
;			   be in the same reference frame.
;-
;------------------------------------------------------------------------------

PRO MDAP_MATCH_RESOLUTION_TPL2OBJ_WAVE_WITH_ZERO_OFFSET, $
		positive_offset, unm_wave, waverange, zero_dispersion=zero_dispersion

	if n_elements(zero_dispersion) eq 0 then $
	    zero_dispersion = 0.0		; Base level dispersion

	sz = size(positive_offset)
	np=sz[1]				; Number of unmasked pixels

	indx=where(positive_offset lt zero_dispersion)
	if indx[0] eq -1 then $
	    message, 'Template spectrum is at lower resolution at all wavelengths!'

	if n_elements(indx) eq np then begin
	    print, 'Entire spectrum is valid: ' + MDAP_STC(unm_wave[0]) + 'A - ' + $
						 MDAP_STC(unm_wave[np-1]) + 'A'
	    waverange = [ unm_wave[0], unm_wave[np-1] ]
	    return
	end

	print, 'Input valid range is ' + MDAP_STC(unm_wave[0]) + 'A - ' + $
	       MDAP_STC(unm_wave[np-1]) + 'A'

	;-------------------------------------------------------------------------------
	; Find the largest contiguous region of pixels that requires no offset

	; Start counting from the first valid wavelength
	indx=where(positive_offset gt zero_dispersion)	; Pixels that require an offset
	incw=[ unm_wave[0], unm_wave[np-1] ]		; Initialize to full range
	j=0						; Start at low wavelength
	while indx[0] ne -1 and j lt np-1 do begin
	    j++						; Increment j
	    incw[0] = unm_wave[j]			; Update wavelength range
	    indx=where(positive_offset[j:np-1] gt zero_dispersion)	; Update indices
	endwhile

	; Start counting from the last valid wavelength
	indx=where(positive_offset gt zero_dispersion)		; Pixels that require an offset
	decw=[ unm_wave[0], unm_wave[np-1] ]		; Initialize to full range
	j=np-1						; Start at high wavelength
	while indx[0] ne -1 and j gt 0 do begin
	    j--						; Decrement j
	    decw[1] = unm_wave[j]			; Update wavelength range
	    indx=where(positive_offset[0:j] gt zero_dispersion)	; Update indices
	endwhile

	; Select the largest region
	if decw[1]-decw[0] gt incw[1]-incw[0] then begin
	    print, 'DEC: Truncating to ' + MDAP_STC(decw[0]) + 'A - ' + MDAP_STC(decw[1]) + 'A'
	    waverange=decw
	endif else begin
	    print, 'INC: Truncating to ' + MDAP_STC(incw[0]) + 'A - ' + MDAP_STC(incw[1]) + 'A'
	    waverange=incw
	endelse

END


PRO MDAP_MATCH_RESOLUTION_TPL2OBJ, $
		tpl_flux, tpl_ivar, tpl_mask, tpl_wave, tpl_sres, wave, sres, tpl_soff, $
		no_offset=no_offset, variable_offset=variable_offset

	sz=size(tpl_flux)
	nt = sz[1]				; Number of template spectra
;	ns = sz[2]				; Number of spectral channels
	tpl_soff=dblarr(nt)			; Initialize the baseline offset to 0.0
	sig2fwhm = 2.0d*sqrt(alog(4.0d))	; Factor to convert sigma to FWHM
	c=299792.458d				; Speed of light in km/s

	sigma_mask = intarr(nt)

	print, 'Number of templates: ', nt

	for i=0,nt-1 do begin
	    print, '    Processing template: ', i+1

	    dw = tpl_wave[i,1] - tpl_wave[i,0]		; Linear step in wavelength
	    unmasked=where(tpl_mask[i,*] lt 1, complement=masked)	; Select the valid pixels

	    print, 'Spectrum length: ', n_elements(reform(tpl_flux[i,*]))
	    if unmasked[0] eq -1 then begin
		print, 'ERROR: Entire spectrum masked!'
		continue
	    endif

;	    print, size(unmasked)
	    print, 'Number of unmasked pixels: ', n_elements(unmasked)
	    num = n_elements(unmasked)

	    if masked[0] eq -1 then begin
	    	print, 'Number of masked pixels: ', 0
	    endif else $
	    	print, 'Number of masked pixels: ', n_elements(masked)

	    zero_vals = where(tpl_flux[i,unmasked] eq 0.)

	    if zero_vals[0] eq -1 then begin
	    	print, 'Number of zero valued fluxes: ', 0
	    endif else $
	    	print, 'Number of zero valued fluxes: ', n_elements(zero_vals)

	    ; Object resolution interpolated to template wavelengths
	    unm_wave = reform(tpl_wave[i,unmasked])
	    obj_sres=interpol(sres, wave, unm_wave)
;	    print, size(obj_sres)
;	    mydevice=!D.NAME
;	    set_plot, 'PS'
;	    device, filename='res_plot.ps'
;	    plot, tpl_wave[i,unmasked], tpl_sres[i,unmasked]
;	    oplot, tpl_wave[i,unmasked], obj_sres, linestyle=2
;	    device, /close
;	    set_plot, mydevice
;	    stop

;	    print, size(tpl_wave[i,unmasked])
	    ; Variance (in angstroms) of Gaussian needed to match template and object resolution
	    diff_var_w = (unm_wave/sig2fwhm)^2*(1.0d/(obj_sres)^2 $
			  - 1.0d/(reform(tpl_sres[i,unmasked]))^2)
;	    print, size(diff_var_w)
;	    stop
	    diff_var_v = (c/unm_wave)^2*diff_var_w		; Convert to km/s

	    ; Determine a constant offset (in km/s) required to ensure that the
	    ;     Gaussian convolution kernel is always above the minimum width
	    positive_offset = reform((c*MDAP_MINIMUM_CNVLV_SIGMA(dw)/unm_wave)^2 - diff_var_v)

	    zero_dispersion = (c*MDAP_MIN_SIG_GAU_APPROX_DELTA(dw)/unm_wave)^2

	    if ~keyword_set(no_offset) then begin		; Allow any offset
		tpl_soff[i] = max([ positive_offset, 0.0d ])
		print, 'Offset: '+MDAP_STC(tpl_soff[i])
	    endif else begin					; Force offset to be 0
		MDAP_MATCH_RESOLUTION_TPL2OBJ_WAVE_WITH_ZERO_OFFSET, positive_offset, unm_wave, $
								     waverange, $
								     zero_dispersion=zero_dispersion
		MDAP_SELECT_WAVE, reform(tpl_wave[i,*]), waverange, indx, complement=comp
		if comp[0] ne -1 then $
		    tpl_mask[i,comp] = 1.0

		indx=where(tpl_mask[i,unmasked] lt 1)	; Save the unmasked indices wrt the old mask
		if indx[0] eq -1 then begin
		    print, 'ERROR: All pixels masked!'
		    continue
		endif

		print, 'New number of unmasked pixels: ', n_elements(indx)

		unmasked = where(tpl_mask[i,*] lt 1)	; Update unmasked
		num = n_elements(unmasked)		; Update the number of unmasked pixels
		unm_wave = reform(tpl_wave[i,unmasked])	; Update the unmasked wave vector 
		obj_sres = obj_sres[indx]		; Update the object spectral resolution map
		diff_var_v = diff_var_v[indx]		; Update the variances

;		print, n_elements(unmasked), n_elements(unm_wave), n_elements(diff_var_v)

		tpl_soff[i] = 0.			; Offset is zero
	    endelse

	    indx = where(diff_var_v + tpl_soff[i] lt 0)
	    diff_sig_w = unm_wave*sqrt(diff_var_v + tpl_soff[i])/c	; Conv sigma in angs
	    diff_sig_w[indx] = 0.0			; Set lt 0 to 0 (within range of delta func)
	    tpl_soff[i] = sqrt(tpl_soff[i])		; Convert from variance to sigma

;	    print, min(diff_sig_w), max(diff_sig_w)

	    ; If a variable or no offset is allowed, go ahead and do the convolution
	    if keyword_set(variable_offset) or keyword_set(no_offset) then begin
		; Perform the convolution
		; TODO: This is slow; farm it out to a C/C++ program?
;		tpl_flux[i,unmasked] = (MDAP_CONVOL_SIGMA(unm_wave, reform(tpl_flux[i,unmasked]), $
;							  unm_wave, diff_sig_w))[*]
;		stop

		; Include the varances
		MDAP_CONVOL_SIGMA_WITH_IVAR, unm_wave, reform(tpl_flux[i,unmasked]), $
					     reform(tpl_ivar[i,unmasked]), unm_wave, $
					     diff_sig_w, conv, conv_ivar

		; TODO: Do something within mdap_convol_sigma_with_ivar?
		; Mask the first and last few pixles due to errors in the convolution
		sigma_mask[i] = ceil(sig2fwhm*max(diff_sig_w)/dw)
		tpl_mask[i,unmasked[0:sigma_mask[i]-1]] = 1.
		tpl_mask[i,unmasked[num-1-sigma_mask[i]:num-1]] = 1.
					     
		tpl_flux[i,unmasked] = conv
		tpl_ivar[i,unmasked] = conv_ivar

		; Calculate the new resolution
;		tpl_sres[i,unmasked] = 1.0d/sqrt(1.0d/(tpl_sres[i,unmasked])^2 + $
;				       (sig2fwhm*diff_sig_w/unm_wave)^2)

		tpl_sres[i,unmasked] = obj_sres
	    endif

	endfor

;	print, sigma_mask

	; If not forcing a single offset, the procedure is done
	if keyword_set(variable_offset) or keyword_set(no_offset) then $
	    return

	max_offset=max(tpl_soff)			; Set the offset to the maximum value
	tpl_soff[*]=max_offset
	print, 'Offset: '+MDAP_STC(max_offset)		; Alert the user

	;-------------------------------------------------------------------------------
	; Convolve the templates to a common resolution wrt the object spectra
	; with a common offset
	for i=0,nt-1 do begin
	    print, '    Processing template: ', i+1

	    dw = tpl_wave[i,1] - tpl_wave[i,0]		; Linear step in wavelength
	    unmasked=where(tpl_mask[i,*] lt 1)		; Select the valid pixels
	    if unmasked[0] eq -1 then $
		continue
	    num=n_elements(unmasked)

	    ; Object resolution interpolated to template wavelengths
	    unm_wave = reform(tpl_wave[i,unmasked])
	    obj_sres=interpol(sres, wave, unm_wave)

	    ; Stddev in angstroms of Gaussian needed to match template and object resolution
	    ; Same as above, just all in one step
	    diff_sig_w = unm_wave*sqrt( (tpl_soff[i]/c)^2 + 1.0/(sig2fwhm*obj_sres)^2 $
					- 1.0/(sig2fwhm*reform(tpl_sres[i,unmasked]))^2 )
	    
	    ; Perform the convolution
	    ; TODO: This is slow; farm it out to a C/C++ program?
;	    tpl_flux[i,unmasked] = (MDAP_CONVOL_SIGMA(tpl_wave[i,unmasked], tpl_flux[i,unmasked], $
;						      tpl_wave[i,unmasked], diff_sig_w))[unmasked]

	    ; Include the varances
	    MDAP_CONVOL_SIGMA_WITH_IVAR, unm_wave, reform(tpl_flux[i,unmasked]), $
					 reform(tpl_ivar[i,unmasked]), unm_wave, $
					 diff_sig_w, conv, conv_ivar

	    ; TODO: Do something within mdap_convol_sigma_with_ivar?
	    ; Mask the first and last few pixles due to errors in the convolution
	    sigma_mask[i] = ceil(sig2fwhm*max(diff_sig_w)/dw)
	    tpl_mask[i,unmasked[0:sigma_mask[i]-1]] = 1.
	    tpl_mask[i,unmasked[num-1-sigma_mask[i]:num-1]] = 1.
					     
	    tpl_flux[i,unmasked] = conv
	    tpl_ivar[i,unmasked] = conv_ivar

	    ; Calculate the new resolution
;	    tpl_sres[i,unmasked] = 1.0d/sqrt(1.0d/(tpl_sres[i,unmasked])^2 + $
;					     (sig2fwhm*diff_sig_w/unm_wave)^2)

	    tpl_sres[i,unmasked] = obj_sres
	endfor

END


;	MDAP_SELECT_WAVE, tpl_wave[i,*], decw, indx
;	;-------------------------------------------------------------------------------
;	; Trim spectra (and associated vectors) if valid indices has changed
;	nind=n_elements(indx)
;	if nind ne np then begin
;	    tpl_wave[i,0:nind-1] = tpl_wave[i,indx]	; Copy the data
;	    tpl_flux[i,0:nind-1] = tpl_flux[i,indx]
;	    tpl_sres[i,0:nind-1] = tpl_sres[i,indx]
;	    tpl_mask[i,0:nind-1] = tpl_mask[i,indx]
;	    if nind lt ns then begin
;		tpl_wave[i,nind:ns-1] = 0.		; Erase the old data
;		tpl_flux[i,nind:ns-1] = 0.
;		tpl_sres[i,nind:ns-1] = 0.
;		tpl_mask[i,nind:ns-1] = 1.
;	    endif


	    

