;+
; NAME:
;	MDAP_ONSKY_XY_ARCS2CENTER
;
; PURPOSE:
;	Compute the distance in arcseconds of each spaxel from the center of the
;	targetted object using the keywords 'OBJRA' and 'OBJDEC' from the
;	supplied header.  As this is a calculation of distance, the fomula
;	includes the approximate correction for the cosine of the declination.
;
; CALLING SEQUENCE:
;	MDAP_ONSKY_XY_ARCS2CENTER, header, skyx, skyy
;
; INPUTS:
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
;	- Check that other procedures expect input that has the distance from
;	  the center!
;
; BUGS:
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;	09 Sep 2014: (KBW) Original Implementation
;-
;------------------------------------------------------------------------------

PRO MDAP_ONSKY_XY_ARCS2CENTER, $
		header, skyx, skyy

	; TODO: Does not check that keywords exist!
	targx = SXPAR(header, 'OBJRA')			; Read the target RA ...
	targy = SXPAR(header, 'OBJDEC')			; ... and declination

	skyx = (skyx-targx)*3600.*cos(targy*!DtoR)	; Convert to arcsecond distance wrt object
	skyy = (skyy-targy)*3600.

END

