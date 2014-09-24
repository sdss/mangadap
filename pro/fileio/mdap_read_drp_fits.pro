;+
; NAME:
;   MDAP_READ_DRP_FITS
;
; PURPOSE:
;	Read in DRP produced fits files.  Primarily a wrapper function that
;	calls separate routines for reading RSS or CUBE formats into a single
;	output format.
;
; CALLING SEQUENCE:
;	MDAP_READ_DRP_FITS,$
;		file, header, flux, ivar, wave, skyx, skyy, cdelt1, cdelt2, type=type, $
;		version=version
;
; INPUTS:
;	file string
;		Name of the fits file containing the input in RSS or CUBE data.
;
;		Fits extensions for both formats are explained in the TRM:
;			https://trac.sdss.org/wiki/MANGA/TRM/TRM_v1.0/datamodel
;
;		Flux units are expected to be 1e-17 ergs/s/cm^{-2}/angstrom/pixel
;
;		Wavelength solution is expected to be log10-linear.
;
; OPTIONAL INPUTS:
;	type string
;		File type.  Must be either 'RSS' or 'CUBE'.  If not provided,
;		the type is determined by the dimensionality of the input fits
;		image: RSS=2, CUBE=3.
;
; OPTIONAL KEYWORDS:
;
; OUTPUT:
;	header
;		The fits header, which will be altered by DAP procedures.
;
;	flux dblarr[N][T]
;		Array of N galaxy spectra, each with spectral channels (pixels).
;		As explained in mdap_read_drp_rss.pro, RSS spectra are ordered
;		as in the input fits file.  As explained in
;		mdap_read_drp_cube.pro, CUBE files have N=NX*NY where NX and NY
;		are the number of spaxels in the X and Y on-sky dimension,
;		respectively.  The array is ordered such that the (i,j)th
;		spectrum is in index [i*NY+j,*] of the array.
;
;	ivar dblarr[N][T]
;		Inverse variances of the spectra.
;		
;	mask dblarr[N][T]
;		Bad pixel mask.
;		
;	wave dblarr[T]
;		Wavelength coordinates of each pixel, in accordance with the DRP
;		output.  That is, the vector is expected to have a constant step
;		in log10(lambda); however, the coordinates are in angstroms, not
;		log10(angstroms).
;
;	sres dblarr[T]
;		Median spectral resolution (R=lamda/delta lamba) as a function
;		of wavelength for all fibers.
;
;	skyx dblarr[N][T]
;		Spatial X coordinate in arcseconds for each spectral channel,
;		with 0 at the center of the field of view.
;
;	skyy dblarr[N][T]
;		Spatial Y coordinate in arcseconds for each spectral channel,
;		with 0 at the center of the field of view.
;
; OPTIONAL OUTPUT:
;	unit string
;		String representation of the units as determined by the CUNIT1
;		and CUNIT2 header keywords.
;
;	version string
;		Module version. If requested, code simply returns the version
;		flag.
;
; TODO:
;	- accommodate covariance matrices
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;	HEADFITS()
;	READFITS()
;	SXPAR()
;	MDAP_DRP_FILETYPE
;	MDAP_CUBE_ONSKY_XY
;	MDAP_RSS_ONSKY_XY
;	MDAP_RESTRUCTURE_CUBE
;	MDAP_RESTRUCTURE_RSS
;	MDAP_ONSKY_XY_ARCS2CENTER
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;	09 Sep 2014: (KBW) Original implementation. Adapted from mdap_read_datacube
;	09 Sep 2014: (KBW) Include mask and instrumental resolution
;	12 Sep 2014: (KBW) New version of MDAP_DRP_FILETYPE
;	12 Sep 2014: (KBW) Remove cdelt? as parameters.  Pixel scale is now
;		           determined using MDAP_GET_SPAXEL_SIZE
;	15 Sep 2014: (KBW) Change to useage of MDAP_OFFSET_XY_CENTER,
;			   MDAP_WCS_UNITS, and MDAP_WCSUNIT2ARCSEC
;	15 Sep 2014: (KBW) Added unit as an optional output
;-
;------------------------------------------------------------------------------

PRO MDAP_READ_DRP_FITS,$
		file, header, flux, ivar, mask, wave, sres, skyx, skyy, type=type, unit=unit, $
		version=version

	version_module = '0.1'				; Version

	; If the version is requested, print it then quit
	if n_elements(version) ne 0 then begin
	    version = version_module
	    return
	endif

	header=HEADFITS(file, exten=1)			; Read the header (EXTENSION IS IMPORTANT!)
	flux=double(READFITS(file, exten_no=1))		; Flux is in extension 1
	ivar=double(READFITS(file, exten_no=2)) 	; Inverse variance is in extension 2
	mask=double(READFITS(file, exten_no=3)) 	; Pixel mask is in extension 3
	wave=double(READFITS(file, exten_no=4))		; Wavelength solution is in extension 4
	sres=double(READFITS(file, exten_no=5))		; Spectral resolution is in extension 5

	MDAP_DRP_FILETYPE, header, type			; Determine/Confirm the file type

	if type eq 'CUBE' then begin
	    MDAP_CUBE_ONSKY_XY, header, skyx, skyy	; Calculate on-sky spaxel coordinates
	    MDAP_RESTRUCTURE_CUBE, flux, ivar, mask	; Restructure the CUBE

	    MDAP_OFFSET_XY_CENTER, header, skyx, skyy	; Offset xy center to target RA/DEC
	    MDAP_WCS_UNITS, header, wcsunit		; Get WCS units
	    MDAP_WCSUNIT2ARCSEC, wcsunit, skyx, skyy	; Convert distance to arcseconds

	    if n_elements(unit) ne 0 then $		; Set unit value, if requested
		unit=wcsunit

	endif else begin
	    MDAP_RSS_ONSKY_XY, file, skyx, skyy		; Get the on-sky spaxel coordinates
	    MDAP_RESTRUCTURE_RSS, flux, ivar, mask	; Restructure the RSS

	    ; XY coordinates of RSS data are already in arcsec offset from IFU center

	endelse

	; These operations are not necessary for the RSS cubes!
;	MDAP_OFFSET_XY_CENTER, header, skyx, skyy	; Offset xy center to target RA/DEC
;	MDAP_WCS_UNITS, header, wcsunit			; Get WCS units
;	MDAP_WCSUNIT2ARCSEC, wcsunit, skyx, skyy	; Convert distance to arcseconds
;
;	if n_elements(unit) ne 0 then $
;	    unit=wcsunit

	; Plot the spectra
;	sz=size(flux)
;	print, sz
;	for i=0,sz[1]-1 do begin
;	    if stddev(flux[i,*]) gt 0 then $
;		plot, wave, flux[i,*]
;	endfor

;	plot, wave, sres
;	print, sres

end



