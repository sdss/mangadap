;+
; NAME:
;	MDAP_MEASURE_SPECTRAL_INDEX
;
; PURPOSE:
;	TODO: Better this description.  Code will CHANGE!
;
;	Measure an absorption-line spectral index defined by a blue and red
;	continuum region and the primary passband.  Errors are calculated based
;	on Cardiel et al. 1998.
;
; CALLING SEQUENCE:
;	MDAP_MEASURE_SPECTRAL_INDEX, wave, flux, ivar, mask, passband, blue_cont, red_cont, $
;				     equiv_width, equiv_width_err, index_mag, index_mag_err, $
;				     title=title, plbound=plbound, /plot 
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
; INTEGRAL MEAN OF THE BAND - the standard definition

;Inputs:
; lam [array] wavelength values of the pseudocontinua
; spc [array] spectra in the wav range of the pseudocontinua

;Outputs:
; pseudo_cont [float] value of the pseudoncontinuum
; pseudo_fit  [array] vector with as many elements as "lam" and equal
;                     to pseudo_cont

PRO MDAP_MEASURE_PSEUDOCONTINUUM, $
		wave, flux, pseudo_cont, pseudo_fit
 
	; FORMULA (1) Worthey et al. 1994
	nwave = n_elements(wave)
	pseudo_cont = int_tabulated(wave, flux, /double) / (wave[nwave-1]-wave[0])
	pseudo_fit = make_array(nwave, /double, value=pseudo_cont)
end

;-------------------------------------------------------------------------------
PRO MDAP_MEASURE_ABSORPTION_LINE_INDEX, $
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
	MDAP_MEASURE_PSEUDOCONTINUUM, wave_resample, flux_resample, blue_pseudocont, blue_pseudofit


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
	MDAP_MEASURE_PSEUDOCONTINUUM, wave_resample, flux_resample, red_pseudocont, red_pseudofit


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


PRO MDAP_MEASURE_SPECTRUM_BANDHEAD, $
		wave, flux, ivar, mask, blue_cont, red_cont, bandhead, bandhead_err

	rindx = where(wave ge  red_cont[0] and wave le  red_cont[1])	; Pixels in the red band
	rcont = median(flux[rindx], /even)			; Get the median continuum value
	rivar = median(ivar[rindx], /even)			; Approximate error in the median

	bindx = where(wave ge blue_cont[0] and wave le blue_cont[1])	; Pixels in the blue band
	bcont = median(flux[bindx], /even)			; Get the median continuum value
	bivar = median(ivar[bindx], /even)			; Approximate error in the median

	bandhead = rcont/bcont				; Measure the strength of the bandhead
	bandhead_err = sqrt( (bandhead/rcont)^2/rivar + (bandhead/bcont)^2/bivar ) ; ... and error

END


PRO MDAP_MEASURE_SPECTRAL_INDEX, $
		wave, flux, ivar, mask, abs_par, equiv_width, equiv_width_err, index_mag, $
		index_mag_err, title=title, plbound=plbound, plot=plot 

	is_bandhead = MDAP_INDEX_IS_BANDHEAD(abs_par)
	if is_bandhead eq 1 and abs_par.units eq 'mag' then $
	    message, 'Unable to use mag units to determine index for '+abs_par.name+'!'

	if is_bandhead eq 1 then begin
	    MDAP_MEASURE_SPECTRUM_BANDHEAD, wave, flux, ivar, mask, abs_par.blue_cont, $
					    abs_par.red_cont, equiv_width, equiv_width_err
	endif else begin
	    MDAP_MEASURE_ABSORPTION_LINE_INDEX, wave, flux, ivar, mask, abs_par.passband, $
						abs_par.blue_cont, abs_par.red_cont, equiv_width, $
						equiv_width_err, index_mag, index_mag_err, $
						title=title, plbound=plbound, plot=plot
	endelse

END
	
