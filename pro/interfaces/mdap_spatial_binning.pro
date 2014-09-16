;+
; NAME:
;	MDAP_SPATIAL_BINNING
;
; PURPOSE:
;	Bin an input set of spectra to a minimum S/N level.
;
; CALLING SEQUENCE:
;	MDAP_SPATIAL_BINNING, flux, ivar, signal, noise, gflag, min_sn, xcoo, ycoo, dx, dy, $
;			      spaxel_dy, binned_indx, binned_flux, binned_ivar, binned_skyx, $
;			      binned_skyy, binned_area, binned_ston, nbinnned=nbinnned, $
;			      sn_calibration=sn_calibration, weight_for_sn=weight_for_sn, /plot
;
; INPUTS:
;	flux dblarr[N][T]
;		Galaxy spectra as produced by MDAP_READ_DRP_FITS.
;
;	ivar dblarr[N][T]
;		Inverse variance of the flux
;
;	signal dblarr[N]
;		Mean galaxy signal per angstrom
;
;	noise dblarr[N]
;		Mean galaxy error per angstrom
;
;	gflag intarr[N]
;		Flag (0=false; 1=true) that the spectrum is 'good' as defined by
;		MDAP_SELECT_GOOD_SPECTRA.  Spectra that are NOT good are ignored.
;
;	min_sn double
;		Minimum S/N (per angstrom) required of the output binned spectra.  
;
;	xcoo dblarr[N]
;		Array containing the x coordinates in arcseconds (0 is the
;		center of the field of view) for each spectrum
;
;	ycoo dblarr[N]
;		Array containing the y coordinates in arcseconds (0 is the
;		center of the field of view) for each spectrum
;
;	dx double
;		Scale arcsec/pixel in X direction
;
;	dy double
;		Scale arcsec/pixel in Y direction
;
; OPTIONAL INPUTS:
;
;	sn_thr double
;		If specified, spectra with S/N lower than this value will be
;		excluded from the analysis.  
;
;	SN_CALIBRATION TODO: TYPE? or flag?
;		If provided, the estimated signal-to-noise (SN_est) is converted
;		into the real signal-to-noise using the empirical calibration
;		function defined in MDAP_CALIBRATE_SN:
;
;			tmp = SN_EST^SN_CALIBRATION[0]/sqrt(n_elements(n_elements_within_bin)
;			SN_REAL = poly(SN_EST,SN_CALIBRATION[1:*])
;
;	user_bin_map  string
;		If provided, the spatial map will be created from the fits file
;		specified by this input. The fits file must contain the CRVAL1,
;		CRVAL2, CDELT1, CDELT2, NAXIS1, NAXIS2, CRPIX1, and CRIX2 header
;		keywords (coordinate units should be in arcseconds; 0,0
;		indicates the center of the field of view).
;
; OPTIONAL KEYWORDS:
;	\plot
;		If set, some plots on X11 terminal will be shown. Not suggested
;		if the task is launched remotely. 
;
;	\weight_for_sn
;		If set, the spectra in the same spatial bin will be weighted by
;		$S/N^2$ before being added. If the voronoi binning scheme is
;		adopted, the S/N in the bin is computed via equation (3) of
;		Cappellari & Copin (2003), and the centers of the spatial bins
;		are computed by weighting spectra coordinates by $S/N^2$.  
;
; OUTPUT:
;
;	binned_indx intarr[N]
;		Indicates in which bin, i=0...B-1, each of the N spectra exist.
;
;	binned_flux dblarr[B][T]
;		The binned spectra of the spatial B bins. i-th spectrum is
;		associated to the i-th bin. 
;
;	binned_ivar dblarr[B][T]
;		Inverse variance of binned spectra.
;
;	binned_xcoo dblarr[B]
;		X-Coordinates in arcsec of the luminosity weighted centers of
;		the spatial bins. 
;
;	binned_ycoo dblarr[B]
;		Y-Coordinates in arcsec of the luminosity weighted centers of
;		the spatial bins. 
;
;	binned_area dblarr[B]
;		Area (in arcsec^2) of each spatial bin.  
;
;	binned_ston dblarr[B]
;		Mean S/N per angstrom reached in each spatial bin. 
;
; OPTIONAL OUTPUT:
;
;	nbinned intarr[B]
;		Number of spectra coadded in each bin.
;
;	version string
;		Module version. If requested, the module is not executed and only
;		the version flag is returned
;
; COMMENTS:
;
; EXAMPLES:
;
; TODO:
;	- Also include full polygon describing the bin in output? (I.e. keep
;	  binned_xvec, binned_yvec)
;	- Include something that handles the covariance
;
; BUGS:
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;	01 Sep 2014: Copied from v0_8 by L. Coccato
;	08 Sep 2014: (KBW) Formatting and comments (incomplete)
;	15 Sep 2014: (KBW) Formatting and edits due to accommodate other changes
;	16 Sep 2014: (KBW) gflag changed from optional to required parameter
;-
;------------------------------------------------------------------------------

PRO MDAP_SPATIAL_BINNING, $
		flux, ivar, mask, signal, noise, gflag, min_sn, xcoo, ycoo, dx, dy, binned_flux, $
		binned_ivar, binned_xcoo, binned_ycoo, binned_area, binned_ston, plot=plot, $
		sn_thr=sn_thr, nbinned=nbinned, sn_calibration=sn_calibration, $
		user_bin_map=user_bin_map, weight_for_sn=weight_for_sn, version=version

	version_module = '0.3'				; Version number

	if n_elements(version) ne 0 then begin		; set version and return
	    version = version_module
	    return
	endif

	sz=size(flux)
	ns=sz[1]					; Number of spectra

	; Get the S/N threshold, if provided
	sn_thr_=0.
	if n_elements(sn_thr) ne 0 then sn_thr_=sn_thr

	; Find which spectra in the 2D map are good (and bad)
	;	good = has a positive and finite noise, a finite signal,
	;	and S/N > threshold
	gflag_= where(gflag eq 1 and abs(signal/noise) ge sn_thr_, compl=bflag_)

;	ind_good = where(noise gt 0 and finite(noise) eq 1  and finite(signal) eq 1 and $
;			 abs(signal/noise) ge sn_thr_,compl=ind_bad)

	apply_voronoi_binning = 1			; Initialize by expecting to Voronoi bin

	if gflag_[0] eq -1 then begin			; No good spectra found! 
	    binned_indx=intarr(ns)
	    binned_indx[*]=-1
	    binned_xcoo=0
	    binned_ycoo=0
	    binned_ston=0
	    binned_area=0
	    apply_voronoi_binning = 0			; Just go to "combine spectra"
	endif

;	TODO: NOT DEFINED YET -----------------------------------
;	; Check the user defined binning scheme
;	endif else if if n_elements(user_bin_map) ne 0 then begin
;	    MDAP_READ_USER_DEFINED_SPATIAL_BINS, user_bin_map, header, binned_indx, success
;	    if success eq 1 then apply_voronoi_binning = 0	; Successful user definition
;	endelse
;	NOT DEFINED YET ----------------------------------------

	; Use the Voronoi binning scheme
	if apply_voronoi_binning eq 1 then begin

	    if keyword_set(plot) then begin			; setup plot
		r = GET_SCREEN_SIZE()
		window, xsize=r[0]*0.4, ysize=r[1]*0.8, retain=2
		loadct, 32
	    endif

	    MDAP_VORONOI_2D_BINNING, xcoo[gflag_], ycoo[gflag_], signal[gflag_], noise[gflag_], $
				     min_sn, binned_indx, binned_xvec, binned_yvec, binned_xcoo, $
				     binned_ycoo, binned_ston, nbinned, scale, plot=plot, $
				     sn_calibration=sn_calibration, /quiet, $
				     weight_for_sn=weight_for_sn

	    ; Set binned_indx to same length as flux
	    MDAP_INSERT_FLAGGED, gflag_, binned_indx, ns
	endif

	print, 'Number of spatial bins: ', mdap_stc(n_elements(nbinned),/integer)

	; Generate the weights to use in combining the spectra

	; TODO: Should this be output from the voronoi routine.  If so, what
	;	should be done in the case of the user-supplied binning scheme.
	;	Does the user supply the weights as well?
	MDAP_GENERATE_BINNING_WEIGHTS, signal, noise, wgt, weight_for_sn=weight_for_sn

	; Combine the spectra
	MDAP_COMBINE_SPECTRA, flux, ivar, mask, binned_indx, wgt, nbinned, binned_flux, binned_ivar

	; Determine the effective on-sky area of each combined spectrum
	MDAP_SPECTRAL_BIN_AREA, dx, dy, nbinned, binned_indx, binned_area

	; TODO: How is binned_area used?  Should it include the weights?

END


