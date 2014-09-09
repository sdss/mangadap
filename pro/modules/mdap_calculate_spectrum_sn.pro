;+
; NAME:
;	MDAP_CALCULATE_SPECTRUM_SN
;
; PURPOSE:
;	Calculate the median signal-to-noise (S/N) per angstrom of a spectrum.
;
; CALLING SEQUENCE:
;	mdap_calculate_spectrum_sn, flux, ivar, wave, sn_per_angstorm, signal=signal, $
;				    noise=noise, /rms, /sum
;
; INPUTS:
;	flux dblarr[]
;		Flux values of spectrum for which to compute the S/N.
;		Units are flux / angstrom .
;
;	ivar dblarr[]
;		Inverse variance of the flux.
;
;	wave dblarr[]
;		Wavelength in angstroms of each pixel in the spectrum.
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
; 	/rms
;		If set, flux is a best-fit model spectrum and ivar is the
;		residual of the data w.r.t. that model.  The S/N is then
;		calculated as:
;
;		      median(flux,/even)         max(wave)-min(wave)
;		S/N = ------------------ * sqrt( ------------------- )
;		      robust_sigma(ivar)           n_elements(wave)
;
;	/sum
;		If set, returns the S/N of the sum, not the median signal to
;		noise.
;
; OUTPUT:
;	sn_per_angstorm double
;		Mean S/N per angstrom over all wavelengths.
;
; OPTIONAL OUTPUT:
;
;	signal double
;		Mean signal, median(flux), per angstrom over all wavelengths.
;
;	noise double
;		Mean noise, median(sqrt(1.0/ivar)), per angstrom over all
;		wavelengths, if /rms is not set.  If it is, this is
;		robust_sigma(ivar).
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; TODO:
;	- use splog
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;	09 Jan 2014: Original implementation by L. Coccato
;	09 Sep 2014: (KBW) Formatting edits
;-
;------------------------------------------------------------------------------

PRO mdap_calculate_spectrum_sn, flux, ivar, wave, sn_per_angstorm, signal=signal, noise=noise, $
				rms=rms, sum=sum

	if keyword_set(rms) and keyword_set(sum) then $
	    print, 'WARNING: Cannot calculate S/N using both /rms and /sum.  Ignoring /sum.'

	if keyword_set(rms) then begin

	    signal = median(flux,/even)			; Calculate median of model signal
	    noise = robust_sigma(ivar)    		; Calculate robust stddev about model
	    lrange=(max(wavelength)-min(wavelength))	; Full wavelength range
	    sn_per_angstorm = signal / noise * sqrt(lrange/n_elements(wavelength))	; ~S/N

	endif else if keyword_set(sum) then begin

	    signal = total(flux)			; Sum of the flux
	    noise = sqrt(total(1.0/temporary(ivar)))	; Propagated error in the sum
	    sn_per_angstrom = signal/noise		; S/N

	endif else begin				; Default behavior

	    signal = median(flux,/even)		; Calculate the median signal
	    noise = median(1./sqrt(ivar),/even);	; Calculate the median noise
	    sn_per_angstorm = median(temporary(flux)*sqrt(temporary(ivar)), /even)	; Median S/N

	endelse

END


;  <S/N> = sum in quadrature of the signal to noises, divided by the  number of elements. (Matt's suggestion)
; dispersion = wavelength[1:*] - wavelength[0:n_elements(wavelength)-2]
; dispersion = [dispersion,dispersion[n_elements(dispersion)-1]]
; sn_per_angstorm = sqrt(total(spectrum/error*sqrt(dispersion))^2)
; signal = total(spectrum)*sqrt(dispersion)/n_elements(wavelength)
; noise =sqrt(total(error^2))*sqrt(dispersion)/n_elements(wavelength)

; <S/N> from Cardiel et al. (2008) formula (Equation 42)
; dispersion = wavelength[1:*] - wavelength[0:n_elements(wavelength)-2]
; dispersion = [dispersion,dispersion[n_elements(dispersion)-1]]
; sn_per_angstorm = 1./n_elements(wavelength) * int_tabulated(wavelength,spectrum/error/sqrt(dispersion))
; signal = 1./n_elements(wavelength) *int_tabulated(wavelength,spectrum/sqrt(dispersion))
; noise = 1./n_elements(wavelength) *int_tabulated(wavelength,error/sqrt(dispersion))

