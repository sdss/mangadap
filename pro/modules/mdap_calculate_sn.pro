;+
; NAME:
;	MDAP_CALCULATE_SN
;
; PURPOSE:
;	Calculate the signal-to-noise for a set of spectra using
;	MDAP_CALCUALTE_SPECTRUM_SN.
;
; CALLING SEQUENCE:
;	MDAP_CALCULATE_SN, flux, ivar, mask, wave, wsel, signal, noise, gflag=gflag, $
;			  version=version, /rms, /sum 
;
; INPUTS:
;	flux dblarr[N][T]
;		Flux values for N spectra with T spectral channels.
;
;	ivar dblarr[N][T]
;		Inverse variance in the flux.
;
;	mask dblarr[N][T]
;		Pixel mask from DRP.
;
;	wave dblarr[T]
;		Wavelength coordinates of each spectral channel.
;
;	wsel intarr[]
;		Array of indices with selected wavelengths to use in the
;		calculation.
;
; OPTIONAL INPUTS:
;	gflag intarr[N]
;		Flag (0=false; 1=true) that the spectrum is 'good' as defined by
;		MDAP_SELECT_GOOD_SPECTRA.  Spectra that are NOT good, are given
;		signal=0. and noise=1.
;
;	version string
;		Version number for calculating S/N.  If provided, the value is
;		set to the current version and the procedure is halted.
;
; OPTIONAL KEYWORDS:
; 	/rms (see MDAP_CALCULATE_SPECTRUM_SN)
;		If set, flux is a best-fit model spectrum and ivar is the
;		residual of the data w.r.t. that model.  The S/N is then
;		calculated as:
;
;		      median(flux,/even)         max(wave)-min(wave)
;		S/N = ------------------ * sqrt( ------------------- )
;		      robust_sigma(ivar)           n_elements(wave)
;
; 	/sum (see MDAP_CALCULATE_SPECTRUM_SN)
;		If set, returns the S/N of the sum, not the median signal to
;		noise.
;
; OUTPUT:
;	signal dblarr[N]
;		Calculated signal for each of N spectra.  If defined, spectra
;		with gflag=0 are omitted from the calculation and given a signal
;		of 0.
;
;	noise dblarr[N]
;		Calculated noise for each of N spectra.  If defined, spectra
;		with gflag=0 are omitted from the calculation and given a noise
;		of 1.
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
;	09 Sep 2014: (KBW) Original implementation
;-
;------------------------------------------------------------------------------

PRO MDAP_CALCULATE_SN, $ 
	flux, ivar, mask, wave, wsel, signal, noise, gflag=gflag, version=version, rms=rms, sum=sum

	version_module = '0.1'				; Version

	; If the version is requested, print it then quit
	if n_elements(version) ne 0 then begin
	    version = version_module
	    return					; Version requested so finish
	endif

	sz=size(flux)
	ns=sz[1]				; Number of spectra
	print, ns
	
	if n_elements(gflag) eq 0 then begin
	    gflag_=intarr(ns)
	    gflag_[*] = 1			; Assume all spectra are good
	endif else $
	    gflag_ = gflag

	signal=dblarr(ns)
	noise=dblarr(ns)
	for i=0,ns-1 do begin
	    if gflag_[i] eq 0 then begin
		signal[i] = 0.
		noise[i] = 1.
		continue
	    endif

	    MDAP_CALCULATE_SPECTRUM_SN, flux[i,wsel], ivar[i,wsel], mask[i,wsel], wave[wsel], $
					ston, signal=signal_, noise=noise_, rms=rms, sum=sum
	    signal[i]=signal_
	    noise[i]=noise_
	
	endfor

END


