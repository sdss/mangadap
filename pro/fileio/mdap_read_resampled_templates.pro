;+
; NAME:
;	MDAP_READ_RESAMPLED_TEMPLATES
;
; PURPOSE:
;
; CALLING SEQUENCE:
;	MDAP_READ_RESAMPLED_TEMPLATES, tpl_out_fits, tpl_wave, tpl_flux, tpl_ivar, tpl_mask, $
;				       tpl_sres, tpl_soff
;
; INPUTS:
;	tpl_out_fits string
;		File name for output.
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;
; OUTPUT:
;	tpl_wave dblarr[S]
;		The vector with the wavelength coordinate of every S spectral
;		channel.  The spectra are binned to a constant step in
;		log10(wavelength), to match the DRP-produced spectra.
;
;	tpl_flux dblarr[T][S]
;		The normalized flux of the template library.
;
;	tpl_ivar dblarr[T][S]
;		The inverse variance of the normalized flux.
;
;	tpl_mask dblarr[T][S]
;		The bad pixel mask of the template library (0/1=good/bad)
;
;	tpl_sres dblarr[S]
;		The spectral resolution (R= delta lambda/lambda) of the template
;		library matched to the object spectra assuming a redshift of
;		z_guess.
;
;	tpl_soff dblarr[T]
;		The list of offsets in the spectral resolution applied to the
;		template needed to ensure that the form of the spectral
;		resolution was consistent with the object spectra (in the case
;		where the object spectrum resolution is better than the template
;		library; see MDAP_MATCH_RESOLUTION_TPL2OBJ).
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
;	READFITS()
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;	23 Oct 2014: (KBW) Original implementation
;-
;------------------------------------------------------------------------------

PRO MDAP_READ_RESAMPLED_TEMPLATES, tpl_out_fits, tpl_wave, tpl_flux, tpl_ivar, tpl_mask, tpl_sres, $
				   tpl_soff

	; TODO: No check are performed!  Should check that all the extensions exist
	    
	tpl_wave=READFITS(tpl_out_fits, exten_no=1)	; Read the wavelength
	tpl_flux=READFITS(tpl_out_fits, exten_no=2)	; Read the flux
	tpl_ivar=READFITS(tpl_out_fits, exten_no=3)	; Read the inverse variance
	tpl_mask=READFITS(tpl_out_fits, exten_no=4)	; Read the mask
	tpl_sres=READFITS(tpl_out_fits, exten_no=5)	; Read the spectral resolution
	tpl_soff=READFITS(tpl_out_fits, exten_no=6)	; Read the resolution offset
END


