;+
; NAME:
;	MDAP_GET_SPECTRAL_INDEX
;
; PURPOSE:
;	TODO: Better this description.  Code will CHANGE!
;
;	Measure an absorption-line spectral index defined by a blue and red
;	continuum region and the primary passband.  Errors are calculated based
;	on Cardiel et al. 1998.
;
; CALLING SEQUENCE:
;	MDAP_GET_SPECTRAL_INDEX, wave, flux, ivar, mask, passband, blue_cont, red_cont, $
;				 equiv_width, equiv_width_err, index_mag, index_mag_err, $
;				 title=title, plbound=plbound, /plot 
;
; INPUTS:
;	wave dblarr[S]
;		Wavelength of every pixel S.  The sampling does NOT need to be
;		linear in wavelength.
;
;	flux dblarr[S]
;		Flux at every pixel S.
;
;	ivar dblarr[S]
;		Inverse variance at every pixel S.
;
;	mask dblarr[S]
;		Bad pixel mask for each pixel S; 0-good, 1-bad.  TODO: This is
;		currently NOT used.
;
;	abs_par AbsorptionIndex[]
;		Structure that contains the defining parameters of the
;		absorption-line index.
;
; OPTIONAL INPUTS:
;	title string
;		Title for the output plot.
;
;	plbound dblarr[2]
;		Y-axis boundaries of the plot, [ plbound[0]*midplot,
;		plbound[1]*midplot], where midplot is the value of the spectrum
;		at its central wavelength.  The default is plbound=[0.6,1.35].
;
; OPTIONAL KEYWORDS:
;	/plot
;		Produce a plot showing the results of the index measurement.
;
; OUTPUT:
;	equiv_width double
;	equiv_width_err double
;		Spectral index measured as an equivalent width in angstroms and
;		its error.  (formula 2 in Worthey et al. 1994)
;
;	index_mag double
;	index_mag_err double
;		Spectral index measured in magnitudes and its error.  (formula 3
;		in Worthey et al. 1994)
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
;	04 Nov 2014: (KBW) Adapted from L. Coccato's version (v0.8)
;-
;------------------------------------------------------------------------------

;-------------------------------------------------------------------------------
; Bin edge [0] is the lower edge of pixel x[0], the upper edge of pixel x[0] is
; bin edge [1], etc.  The values of x are expected to be linearly or
; geometrically spaced.
FUNCTION MDAP_BIN_EDGES, $
		x, geometric=geometric

	n = n_elements(x)
;	print, 'edges: ', n

	if n eq 1 then $
	    message, 'Cannot determine edges of a single bin!'

;	print, 'edges: ', x[0], x[n-1]
	xx = x
	if keyword_set(geometric) then $
	    xx = alog(x)
;	print, 'edges: ', xx[0], xx[n-1]
	
	edges = ([ xx[0]-(xx[1]-xx[0])/2.0d , xx[0:n-2] + (xx[1:n-1]-xx[0:n-1])/2.0, $
		   xx[n-1]+(xx[n-1]-xx[n-2])/2.0d ])

;	print, 'edges: ', edges[0], edges[n]
	if keyword_set(geometric) then $
	    edges = exp(edges)
;	print, 'edges: ', edges[0], edges[n]

	return, edges
END

;-------------------------------------------------------------------------------
; Performs the integral:
;
;	integral = \int_xrange[0]^xrange[1] y dx
;
; using a Riemann sum, where x is expected to be at the center of the pixel.
; The error is estimated by a simple propagation of the error (i.e. error in the
; sum).
;
; x is expected to be sorted and contiguous; x coordinates are assumed to be
; linearly or geometrically spaced; the provided coordinate is expected to be at
; the (geometric) center of the pixel.
PRO MDAP_INTEGRATE_PIXELIZED_VALUE, $
		x, y, ye, mask, xrange, integral, integral_err, err=err, geometric=geometric

	err=0							; Initialize error
	n = n_elements(x)					; Number of pixels
	bin_edges = MDAP_BIN_EDGES(x, geometric=geometric)	; The n+1 edges of the n pixels
	if xrange[1] lt bin_edges[0] or xrange[0] gt bin_edges[n] then begin	; No overlap
;	    print, xrange[0], xrange[1], bin_edges[0], bin_edges[n]
	    integral = 0.0d
	    integral_err = 1.0d
	    err=1						; Return an error
	    return
	endif

	; Set the values to 0.0 for masked pixels
	indx = where(mask gt 0.5)
	ym = y
	yme = ye
	if indx[0] ne -1 then begin
	    ym[indx] = 0.0d
	    yme[indx] = 0.0d
	endif

	; Interval is smaller than the nearest pixel
	indx = where(bin_edges gt xrange[0] and bin_edges lt xrange[1])
	if indx[0] eq -1 then begin
	    indx = where(bin_edges gt xrange[1])
	    integral = ym[indx[0]-1]*(xrange[1]-xrange[0])
	    integral_err = yme[indx[0]-1]*(xrange[1]-xrange[0])
	    return
	endif

	; Initialize output
	integral = 0.0d
	integral_err = 0.0d

	ni = n_elements(indx) 				; Number of edges
	; Add partial pixel at the blue end
	if indx[0] gt 0 then begin
	    integral = integral + ym[indx[0]-1]*(bin_edges[indx[0]]-xrange[0])
	    integral_err = integral_err + (yme[indx[0]-1]*(bin_edges[indx[0]]-xrange[0]))^2
	endif

	; Add full pixels
	if ni ge 2 then begin
	    integral = integral + total( ym[indx[0:ni-2]]*(bin_edges[indx[1:ni-1]] $
							  - bin_edges[indx[0:ni-2]]) )
	    integral_err = integral_err + total( (yme[indx[0:ni-2]]*(bin_edges[indx[1:ni-1]] $
								     - bin_edges[indx[0:ni-2]]))^2 )
	endif

	; Add partial pixel at the red end
	if indx[ni-1] lt n then begin
	    integral = integral + ym[indx[ni-1]]*(xrange[1]-bin_edges[indx[ni-1]])
	    integral_err = integral_err + ( yme[indx[ni-1]]*(xrange[1]-bin_edges[indx[ni-1]]) )^2
	endif

	integral_err = sqrt(integral_err)
END

;-------------------------------------------------------------------------------
; Measure the pseudo-continuum of the band using:
;
;  			\int_l0^l1 flux dl
;	continuum =   ---------------------
;			  \int_l0^l1 dl
;
; while accounting for the masked pixels (by not integrating over them) and
; providing a formal estimate of the error.
PRO MDAP_GET_PSEUDOCONTINUUM, $
		wave, flux, ivar, mask, passband, pseudo_continuum, pseudo_continuum_error, $
		wave_integral=wave_integral, err=err, geometric=geometric

	n = n_elements(wave)				; Number of pixels
	unity = make_array(n, /double, value=1.0d)

	sige = sqrt(1.0d/ivar)
	indx = where(ivar le 0.0d)
	mask_ = mask
	if indx[0] ne -1 then begin
	    mask_[indx] = 1.0d
	    sige[indx] = 1.0d
	endif

	MDAP_INTEGRATE_PIXELIZED_VALUE, wave, unity, unity, mask_, passband, wave_integral, dummy, $
					err=err, geometric=geometric
	if err eq 1 then begin		; Throw an error because there was no overlap
	    pseudo_continuum = 0.0
	    pseudo_continuum_error = 0.0
	    wave_integral = 0.0
	    return
	endif

	MDAP_INTEGRATE_PIXELIZED_VALUE, wave, flux, sige, mask_, passband, flux_integral, $
					flux_integral_err, geometric=geometric

	pseudo_continuum = flux_integral/wave_integral
	pseudo_continuum_error = flux_integral_err/wave_integral
END


;-------------------------------------------------------------------------------
PRO MDAP_GET_ABSORPTION_LINE_INDEX, $
		wave, flux, ivar, mask, passband, blueband, redband, equiv_width, equiv_width_err, $
		index_mag, index_mag_err, err=err, geometric=geometric

;	print, 'abs:', wave[0], wave[n_elements(wave)-1]
;	print, 'abs:', passband
;	print, 'abs:', blueband
;	print, 'abs:', redband
	
	; Unmasked pixels in the blue band
	bindx = where(wave ge blueband[0] and wave le blueband[1] and mask lt 0.5)
	if bindx[0] eq -1 then $
	    message, 'Blue wavelengths unavailable!'
	sn_blue = total(flux[bindx]*sqrt(ivar[bindx]))	; Approximate S/N in the blue region
	nblue = n_elements(bindx)

	bcen = (blueband[0] + blueband[1])/2.0			; Blue band center
	MDAP_GET_PSEUDOCONTINUUM, wave, flux, ivar, mask, blueband, blue_cont, blue_cont_err, $
				  wave_integral=blue_width, err=err, geometric=geometric

;	print, 'abs: ', sn_blue
;	print, 'abs: ', nblue
;	print, 'abs: ', bcen
;	print, 'abs: ', blue_cont
;	print, 'abs: ', blue_cont_err
;	print, 'err: ', err

	if err eq 1 then $
	    return				; No wavelength integral so return with error

	; Unmasked pixels in the red band
	rindx = where(wave ge  redband[0] and wave le  redband[1] and mask lt 0.5)
	if rindx[0] eq -1 then $
	    message, 'Red wavelengths unavailable!'
	sn_red = total(flux[rindx]*sqrt(ivar[rindx]))	; Approximate S/N in the red region
	nred = n_elements(rindx)

	rcen = (redband[0] + redband[1])/2.0				; Red band center
	MDAP_GET_PSEUDOCONTINUUM, wave, flux, ivar, mask, redband, red_cont, red_cont_err, $
				  wave_integral=red_width, err=err, geometric=geometric

;	print, 'abs: ', rindx
;	print, 'abs: ', sn_red
;	print, 'abs: ', nred
;	print, 'abs: ', rcen
;	print, 'abs: ', red_cont
;	print, 'abs: ', red_cont_err
;	print, 'err: ', err

	if err eq 1 then $
	    return				; No wavelength integral so return with error

	; Unmasked pixels in the main passband
	pindx = where(wave ge passband[0] and wave le passband[1] and mask lt 0.5)
	if pindx[0] eq -1 then $
	    message, 'Passband wavelengths unavailable!'
	sn_pass = total(flux[pindx]*sqrt(ivar[pindx]))	; Approximate S/N in the red region
	npass = n_elements(pindx)

	pcen = (passband[0] + passband[1])/2.0				; Passband center

	; Calculate the slope and intercept of the continuum line
	slope = (red_cont - blue_cont) / (rcen - bcen)
	intercept = blue_cont - bcen*slope

	continuum = slope * wave + intercept			; Continuum level

	; Get the equivalent width integral (index in angstrom); equation 2 Worthey et al. 1994
	integrand = 1.0d - flux / continuum			; EW integrand
	integrand_err = abs(1.0d / continuum / sqrt(ivar))
	mask_ = mask
	indx = where(ivar le 0.0d)
	if indx[0] ne -1 then begin
	    mask_[indx] = 1.0d
	    integrand_err[indx] = 1.0d
	endif

	MDAP_INTEGRATE_PIXELIZED_VALUE, wave, integrand, integrand_err, mask_, passband, $
					equiv_width, equiv_width_err, err=err, geometric=geometric
	if err eq 1 then $				; No overlap so return in error
	    return

	; Get the magnitiude-based integral; equation 3 Worthey et al. 1994
	n = n_elements(wave)				; Number of pixels
	unity = make_array(n, /double, value=1.0d)
	MDAP_INTEGRATE_PIXELIZED_VALUE, wave, unity, unity, mask_, passband, pass_width, dummy, $
					geometric=geometric

;	print, 'abs: ', pass_width

	integrand = flux / continuum				; EW integrand
	integrand_err = abs(1.0d / continuum / sqrt(ivar))
	MDAP_INTEGRATE_PIXELIZED_VALUE, wave, integrand, integrand_err, mask_, passband, $
					index_mag, index_mag_err, geometric=geometric
	
	index_mag_err = -2.5 * index_mag_err/index_mag/alog(10.0d)	; Error done first
	index_mag = -2.5 * alog10(index_mag/pass_width)			; Index in magnitudes

	; Calculate the errors (See Cardiel et al. 1998, equations 41-46)

	; This is exactly their equation assuming a constant pixel size.
	;
	;    SN_cardiel = (sn_blue+sn_red+sn_pass)/(nblue+nred+npass)/sqrt(dw)	; formula 42
	;
	; Instead, I have allowed for a variable pixel size.  Over the band, the
	; pixel size probably does not vary that much.  Try to approximate dw ~
	; Dw/N, such that equation 42 becomes
	SN_cardiel = (sn_blue+sn_red+sn_pass) * sqrt( (nblue+nred+npass) $
						      / (blue_width+red_width+pass_width) )
	
	q1=(rcen-pcen)/(rcen-bcen)
	q2=(pcen-bcen)/(rcen-bcen)
	c2= sqrt(1./pass_width + q1^2 / blue_width + q2^2 / red_width)		; formula 44

	c1=pass_width*c2							; formula 43
	equiv_width_err=(c1-c2*equiv_width)/SN_cardiel				; formula 41

	c3=2.5*c2*alog10(exp(1))						; formula 46
	index_mag_err=c3/SN_cardiel						; formula 45

	; TODO: How do these errors compare with the naive error propagation in
	;       MDAP_INTEGRATE_PIXELIZED_VALUE?

END

PRO MDAP_GET_SPECTRUM_BANDHEAD, $
		wave, flux, ivar, mask, blue_cont, red_cont, bandhead, bandhead_err, err=err

	err = 0	

	; Unmasked pixels in the red band
	rindx = where(mask lt 0.5 and wave ge red_cont[0] and wave le red_cont[1])
	if rindx[0] eq -1 then begin				; No red pixels
	    err = 1						; Set error and return
	    return
	endif
	rcont = median(flux[rindx], /even)			; Get the median continuum value
	rivar = median(ivar[rindx], /even)			; Approximate error in the median
;	print, 'bandhead: ', rcont, rivar

	; Unmasked pixels in the blue band
	bindx = where(mask lt 0.5 and wave ge blue_cont[0] and wave le blue_cont[1])
	if bindx[0] eq -1 then begin				; No blue pixels
	    err = 1						; Set error and return
	    return
	endif
	bcont = median(flux[bindx], /even)			; Get the median continuum value
	bivar = median(ivar[bindx], /even)			; Approximate error in the median
;	print, 'bandhead: ', bcont, bivar

	bandhead = rcont/bcont				; Measure the strength of the bandhead
	bandhead_err = sqrt( (bandhead/rcont)^2/rivar + (bandhead/bcont)^2/bivar ) ; ... and error

;	print, 'bandhead: ', bandhead, bandhead_err

END


PRO MDAP_GET_SPECTRAL_INDEX, $
		wave, flux, ivar, mask, abs_par, equiv_width, equiv_width_err, index_mag, $
		index_mag_err, err=err, geometric=geometric

	is_bandhead = MDAP_INDEX_IS_BANDHEAD(abs_par)
	if is_bandhead eq 1 and abs_par.units eq 'mag' then $
	    message, 'Unable to use mag units to determine index for '+abs_par.name+'!'

	if is_bandhead eq 1 then begin
;	    print, 'is bandhead'
	    MDAP_GET_SPECTRUM_BANDHEAD, wave, flux, ivar, mask, abs_par.blue_cont, $
					abs_par.red_cont, equiv_width, equiv_width_err, err=err
	endif else begin
;	    print, 'is absorption line'
	    MDAP_GET_ABSORPTION_LINE_INDEX, wave, flux, ivar, mask, abs_par.passband, $
					    abs_par.blue_cont, abs_par.red_cont, equiv_width, $
					    equiv_width_err, index_mag, index_mag_err, err=err, $
					    geometric=geometric
	endelse

END
	
;-------------------------------------------------------------------------------
; INTEGRAL MEAN OF THE BAND - the standard definition

;Inputs:
; lam [array] wavelength values of the pseudocontinua
; spc [array] spectra in the wav range of the pseudocontinua

;Outputs:
; pseudo_cont [float] value of the pseudoncontinuum
; pseudo_fit  [array] vector with as many elements as "lam" and equal
;                     to pseudo_cont

PRO MDAP_GET_PSEUDOCONTINUUM_LODO, $
		wave, flux, pseudo_cont, pseudo_fit
 
	; FORMULA (1) Worthey et al. 1994
	nwave = n_elements(wave)
	pseudo_cont = int_tabulated(wave, flux, /double) / (wave[nwave-1]-wave[0])
	pseudo_fit = make_array(nwave, /double, value=pseudo_cont)
end

;-------------------------------------------------------------------------------
PRO MDAP_GET_ABSORPTION_LINE_INDEX_LODO, $
		wave, flux, ivar, mask, passband, blue_cont, red_cont, equiv_width, $
		equiv_width_err, index_mag, index_mag_err, title=title, plbound=plbound, plot=plot 
		
	; TODO: What happens if there are intervening masks?  Interpolate across them?

	; TODO: The rebinning is meant to increase the accuracy of the integral,
	; but a Riemann sum can be done exactly.  It also changes the wavelength
	; sampling to be linear in wavelength.  Is this the primary reason this
	; resampling is done?

	; TODO: Need to mask ivar and include mask

	bindx = where(wave ge blue_cont[0] and wave le blue_cont[1])	; Pixels in the blue band
	if bindx[0] eq -1 then $
	    message, 'Blue wavelengths unavailable!'
	nblue = n_elements(bindx)					; Number of pixels
	bcen = (blue_cont[0] + blue_cont[1])/2.0			; Blue band center

	blue_width = (blue_cont[1]-blue_cont[0])			; Width of the blue band
	dw = (blue_width)/(10*nblue-1)					; Increase sampling
	wave_resample = dindgen(10*nblue) * dw + blue_cont[0]		; Set wavelength range
	flux_resample = interpol(flux, wave, wave_resample, /spline)	; Resample flux
	ivar_resample = interpol(ivar, wave, wave_resample, /spline)	; Resample the variance
	sn_blue = total(flux_resample*sqrt(ivar_resample))		; S/N in the blue region

	; Measure blue pseudo-continuum
	MDAP_GET_PSEUDOCONTINUUM, wave_resample, flux_resample, blue_pseudocont, blue_pseudofit


	rindx = where(wave ge  red_cont[0] and wave le  red_cont[1])	; Pixel in the red band
	if rindx[0] eq -1 then $
	    message, 'Red wavelengths unavailable!'
	nred = n_elements(rindx)					; Number of pixels
	rcen = (red_cont[0] + red_cont[1])/2.0				; Red band center

	red_width = (red_cont[1]-red_cont[0])				; Width of the red band
	dw = red_width/(10*nred-1)					; Increase sampling
	wave_resample = dindgen(10*nred) * dw + red_cont[0]		; Set wavelength range
	flux_resample = interpol(flux, wave, wave_resample, /spline)	; Resample flux
	ivar_resample = interpol(ivar, wave, wave_resample, /spline)	; Resample the variance
	sn_red = total(flux_resample*sqrt(ivar_resample))		; S/N in the red region

	; Measure red pseudo-continuum
	MDAP_GET_PSEUDOCONTINUUM, wave_resample, flux_resample, red_pseudocont, red_pseudofit


	pindx = where(wave ge passband[0] and wave le passband[1])	; Pixels in the band
	if pindx[0] eq -1 then $
	    message, 'Passband wavelengths unavailable!'
	npass = n_elements(pindx)					; Number of pixels
	pcen = (passband[0] + passband[1])/2.0				; Passband center

	pass_width = (passband[1]-passband[0])				; Width of the passband
	dw = pass_width/(10*npass-1)					; Increase sampling
	wave_resample = dindgen(10*npass) * dw + passband[0]		; Set wavelength range
	flux_resample = interpol(flux, wave, wave_resample, /spline)	; Resample flux
	ivar_resample = interpol(ivar, wave, wave_resample, /spline)	; Resample the variance
	sn_pass = total(flux_resample*sqrt(ivar_resample))		; S/N in the passband


	; Calculate the slope and intercept of the continuum line
	slope = (red_pseudocont - blue_pseudocont) / (rcen - bcen)
	intercept = blue_pseudocont - bcen*slope

	continuum = slope * wave_resample + intercept			; Continuum level
	integrand = 1. - flux_resample / continuum			; EW integrand

	equiv_width = int_tabulated(wave_resample, integrand)		; Equivalent width

	index_mag = int_tabulated(wave_resample, flux_resample/continuum)
	index_mag = -2.5 * alog10(index_mag/(passband[1]-passband[0]))	; Index in magnitudes


	; Calculate the errors (See Cardiel et al. 1998, equations 41-46)
	SN_cardiel = (sn_blue+sn_red+sn_pass)/(nblue+nred+npass)/sqrt(dw)	; formula 42
	
	q1=(rcen-pcen)/(rcen-bcen)
	q2=(pcen-bcen)/(rcen-bcen)
	c2= sqrt(1./pass_width + q1^2 / blue_width + q2^2 / red_width)		; formula 44

	c1=pass_width*c2							; formula 43
	equiv_width_err=(c1-c2*equiv_width)/SN_cardiel				; formula 41

	c3=2.5*c2*alog10(exp(1))						; formula 46
	index_mag_err=c3/SN_cardiel						; formula 45

	if ~keyword_set(plot) then $
	    return						; No plot, so return

	;-----------------------------------------------------------------------
	; Produce the plot -----------------------------------------------------
	BASIC_COLORS, black, white, red, green, blue, yellow, cyan, magenta, orange, mint, purple, $
		      pink, olive, lightblue, gray
	;UserSymbol,'Circle',/fill

	if n_elements(title) eq 0 then title=''
	if n_elements(plbound) eq 0 then plbound = [0.6,1.35]
	midplot=interpol(flux, wave, 0.5*(min(wave)+max(wave)))

	inddd = where( wave ge blue_cont[0]-20 and wave le red_cont[1]+20 )
	mn=min(flux(inddd))
	mx=max(flux(inddd))
	rgn=(mx-mn)*1.05

   	yrange=[mn-0.1*rgn,mx+0.1*rgn]
	plot, wave, flux, /nodata, xtitle='!6wavelength [angstrom]', yrange=yrange, ystyle=1, $
	      title=title, xrange=[ blue_cont[0]-10,red_cont(1)+10], xstyle=1
	oplot, wave, flux, thick=2
	oplot, wave[bindx], flux[bindx], thick=3
	oplot, wave[bindx], blue_pseudofit, color=blue, thick=5
	plots, bcen, blue_pseudocont, psym=8

	oplot, wave[rindx], flux[rindx], thick=3
	oplot, wave[rindx], red_pseudofit, color=red, thick=5
	plots, rcen, red_pseudocont, psym=8

	oplot, [bcen, rcen], [blue_pseudocont, red_pseudocont], linestyle=1

	oplot, wave[pindx], flux[pindx], thick=5

	oplot, [blue_cont[0], blue_cont[0]], [0, interpol(flux, wave, blue_cont[0], /spline)], $
	       linestyle=1, color=blue
	oplot, [blue_cont[1], blue_cont[1]], [0, interpol(flux, wave, blue_cont[1], /spline)], $
	       linestyle=1, color=blue
	oplot, [red_cont[0], red_cont[0]], [0, interpol(flux, wave, red_cont[0], /spline)], $
	       linestyle=1, color=red
	oplot, [red_cont[1], red_cont[1]], [0, interpol(flux, wave, red_cont[1], /spline)], $
	       linestyle=1,color=red
	oplot, [passband[0], passband[0]], [0, interpol(flux, wave, passband[0], /spline)], $
	       linestyle=1
	oplot, [passband[1], passband[1]], [0, interpol(flux, wave, passband[1], /spline)], $
	       linestyle=1

	midplot=0.5*(mx+mn)
	xyouts, min(wave[inddd])+21, mx+0.05*rgn, 'Dispersion=  ' + $
		MDAP_ROUND_STR(wave[1]-wave[0],4) + ' (ang/pxl)'
	xyouts, min(wave[inddd])+21, mx-0.00*rgn, 'Blue pscont= ' +MDAP_ROUND_STR(blue_pseudocont,4)
	xyouts, min(wave[inddd])+21, mx-0.05*rgn, 'Red pscont=  ' + MDAP_ROUND_STR(red_pseudocont,4)
	xyouts, min(wave[inddd])+21, mx-0.10*rgn, 'Eq. Width=    '+ MDAP_ROUND_STR(equiv_width,4) $
		+ ' +/- ' + MDAP_ROUND_STR(equiv_width_err,4) + ' (ang)'
	xyouts, min(wave[inddd])+21, mx-0.15*rgn, 'Index=       '+ MDAP_ROUND_STR(index_mag,5) $
		+ ' +/- ' + MDAP_ROUND_STR(index_mag_err,5) + ' (mag)'
	;-----------------------------------------------------------------------

END

