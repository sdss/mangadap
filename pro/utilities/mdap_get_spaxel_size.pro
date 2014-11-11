;+
; NAME:
;	MDAP_GET_SPAXEL_SIZE
;
; PURPOSE:
;	Determine the spaxel size in arcsec using the WCS information or set it
;	to a fixed value for RSS data.  For the latter, a square pixel is
;	adopted with equivalent collecting area of a circular fiber with a
;	diameter of 2 arcseconds and a nominal PSF such that the effective area
;	is increased to 2.5 arcsec.
;
; CALLING SEQUENCE:
;	MDAP_GET_SPAXEL_SIZE, header, dx, dy, type=type
;
; INPUTS:
;	header hdu
;		Header read from the DRP-produced fits file (CUBE or RSS).
;
; OPTIONAL INPUTS:
;	type string
;		The type of fits file that was read.  If not provided, this is
;		determined from the fits header using MDAP_DRP_FILETYPE.
;
; OPTIONAL KEYWORDS:
;
; OUTPUT:
;	dx double
;		Size of the spaxel in x.
;
;	dy double
;		Size of the spaxel in y.
;
; OPTIONAL OUTPUT:
;
; COMMENTS:
;
; EXAMPLES:
;
; TODO:
;	- Does not check that the eigenvectors point in the right direction!
;	- READ THE FULL WCS SYSTEM AND GET THE PIXEL SCALE THAT WAY!
;
; BUGS:
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;	12 Sep 2014: (KBW) Original implementation
;-
;------------------------------------------------------------------------------

PRO MDAP_GET_SPAXEL_SIZE, $
	header, dx, dy, type=type, unit=unit

	if n_elements(type) eq 0 then $				; If type is not provided,
	    MDAP_DRP_FILETYPE, header, type			; get it from the header

	if type eq 'CUBE' then begin
	
	    cd_matrix = dblarr(2,2)				; Read the CD matrix

	    cd_matrix[0,0] = SXPAR(header, "CD1_1")		; Diagonal values
	    cd_matrix[1,1] = SXPAR(header, "CD2_2")

	    cd_matrix[0,1] = SXPAR(header, "CD1_2")		; Non-zero values indicate rotation
	    cd_matrix[1,0] = SXPAR(header, "CD2_1")

	    eval = hqr(elmhes(cd_matrix), /double)		; Get its eigenvalues

	    dx = abs(real_part(eval[0]))			; Pixel size is abs val. of real pt
	    dy = abs(real_part(eval[1]))

	    if n_elements(unit) eq 0 then $
		MDAP_WCS_UNITS, header, unit			; Get WCS units

	    MDAP_WCSUNIT2ARCSEC, unit, dx, dy			; Convert from WCS unit to arcsec

	endif else begin
	    ; TODO: Find a better solution to this
	    dx = sqrt(!PI)*(2.5/2.0)				; area = pi * (2.5/2.0)^2
	    dy = dx
	endelse

END

