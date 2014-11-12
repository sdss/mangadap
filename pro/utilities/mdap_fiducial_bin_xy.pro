;+
; NAME:
;	MDAP_FIDUCIAL_BIN_XY
;
; PURPOSE:
;	Compute a single set of coordinates for each spectrum.  Currently set xy
;	coordinate at central pixel.
;
; CALLING SEQUENCE:
;	MDAP_FIDUCIAL_BIN_XY, skyx, skyy, bskyx, bskyy, version=version
;
; INPUTS:
;	skyx dblarr[N][T]
;		Spatial X coordinate in arcseconds for each spectral channel,
;		with 0 at the center of the field of view.
;
;	skyy dblarr[N][T]
;		Spatial Y coordinate in arcseconds for each spectral channel,
;		with 0 at the center of the field of view.
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;
; OUTPUT:
;	bskyx dblarr[N]
;		Spatial X coordinate in arcseconds to use globally for each
;		spectrum.
;
;	bskyy dblarr[N]
;		Spatial Y coordinate in arcseconds to use globally for each
;		spectrum.
;
; OPTIONAL OUTPUT:
;	version string
;		Module version.  If requested, the module is not executed and
;		only version flag is returned.
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
;	15 Sep 2014: (KBW) Original implementation
;	12 Nov 2014: (KBW) Added version
;-
;------------------------------------------------------------------------------

PRO MDAP_FIDUCIAL_BIN_XY, $
		skyx, skyy, bskyx, bskyy, version=version

	version_module = '0.1'			; Version number
	if n_elements(version) ne 0 then begin	; If version is defined
	    version = version_module		; ... set it to the module version
	    return				; ... and return without doing anything
	endif

	sz=size(skyx)				; Size of the skyx array
	ns=sz[1]				; Number of spectra
	nc=sz[2]				; Number of spectral channels

	bskyx=dblarr(ns)			; Initialize the arrays
	bskyy=dblarr(ns)

	for i=0,ns-1 do begin			; Set fiducial coo to coo at central channel
	    bskyx[i] = skyx[i,nc/2]
	    bskyy[i] = skyy[i,nc/2]
	endfor

END

