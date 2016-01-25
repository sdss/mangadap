;+
; NAME:
;       MDAP_WRITE_RESAMPLED_TEMPLATES
;
; PURPOSE:
;       Write the resolution-matched and resampled template library to a single
;       fits file.
;
; CALLING SEQUENCE:
;       MDAP_WRITE_RESAMPLED_TEMPLATES, tpl_out_fits, library_key, z_guess, $
;                                       mdap_tpl_lib_setup_version, tpl_wave, tpl_flux, tpl_ivar, $
;                                       tpl_mask, tpl_sres, tpl_soff
;
; INPUTS:
;       tpl_out_fits string
;               File name for output.
;
;       library_key string
;               Template library identifier to write to the fits file header.
;
;       z_guess double
;               Guess value of the redshift used to match the spectral
;               resolution of the object spectrum to the rest wavelength
;               spectral resolution of the template library.
;
;       tpl_flux_norm double
;               The normalization applied to all template spectra.
;
;       mdap_tpl_lib_setup_version string
;               Version of the procedures to perform the resolution matching and
;               resampling.
;
;       tpl_wave dblarr[S]
;               The vector with the wavelength coordinate of every S spectral
;               channel.  The spectra are binned to a constant step in
;               log10(wavelength), to match the DRP-produced spectra.
;
;       tpl_flux dblarr[T][S]
;               The normalized flux of the template library.
;
;       tpl_ivar dblarr[T][S]
;               The inverse variance of the normalized flux.
;
;       tpl_mask dblarr[T][S]
;               The bad pixel mask of the template library (0/1=good/bad)
;
;       tpl_sres dblarr[S]
;               The spectral resolution (R= delta lambda/lambda) of the template
;               library matched to the object spectra assuming a redshift of
;               z_guess.
;
;       tpl_soff dblarr[T]
;               The list of offsets in the spectral resolution applied to the
;               template needed to ensure that the form of the spectral
;               resolution was consistent with the object spectra (in the case
;               where the object spectrum resolution is better than the template
;               library; see MDAP_MATCH_RESOLUTION_TPL2OBJ).
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
; BUGS:
;
; PROCEDURES CALLED:
;       MKHDR
;       SXADDPAR
;       WRITEFITS
;       MDAP_FITS_ZERO_NAXES
;       MDAP_FITS_HEADER_ADD_DATE
;
; INTERNAL SUPPORT ROUTINES:
;       MDAP_ADD_TPL_EXTENSION
;
; REVISION HISTORY:
;       23 Oct 2014: (KBW) Original implementation
;-
;------------------------------------------------------------------------------

PRO MDAP_ADD_TPL_EXTENSION, file, data, extname, comment
        MKHDR, hdr, data, /image                        ; Make a minimal image extension header
        SXADDPAR, hdr, 'EXTNAME', extname, comment      ; Add the extension name
        WRITEFITS, file, data, hdr, /append             ; Write the data to an appended extension
END

PRO MDAP_WRITE_RESAMPLED_TEMPLATES, $
                tpl_out_fits, library_key, z_guess, tpl_flux_norm, mdap_tpl_lib_setup_version, $
                tpl_wave, tpl_flux, tpl_ivar, tpl_mask, tpl_sres, tpl_soff

        ; This routine should NOT be called if tpl_out_fits already exists

        ; Identical to MDAP_INITIALIZE_FITS_FILE, but with additional header key words
        MKHDR, hdr, dblarr(1), /extend          ; Create a bare-bones header with extensions
        MDAP_FITS_ZERO_NAXES, hdr               ; Get rid of the axes
        MDAP_FITS_HEADER_ADD_DATE, hdr          ; Add the DAP date

        ; Add some useful keywords
        ; TODO: Any more?
        SXADDPAR, hdr, 'LIBKEY', library_key, 'Library identifier'
        SXADDPAR, hdr, 'ZGUESS', z_guess, 'Guess redshift used to match spectral resolution'
        SXADDPAR, hdr, 'FLXNRM', tpl_flux_norm, 'Flux normalization (units from input library)'
        SXADDPAR, hdr, 'VDAPTPL', mdap_tpl_lib_setup_version, 'mdap_tpl_lib_setup version'

        WRITEFITS, tpl_out_fits, 0, hdr                 ; Write the primary header

        ; Wavelength extension
        MDAP_ADD_TPL_EXTENSION, tpl_out_fits, tpl_wave, 'WAVE', 'Wavelength coordinates (ang)'
        MDAP_ADD_TPL_EXTENSION, tpl_out_fits, tpl_flux, 'FLUX', 'Normalized flux'
        MDAP_ADD_TPL_EXTENSION, tpl_out_fits, tpl_ivar, 'IVAR', 'Inverse variance of flux'
        MDAP_ADD_TPL_EXTENSION, tpl_out_fits, tpl_mask, 'MASK', 'Bad pixel mask (0/1-good/bad)'
        MDAP_ADD_TPL_EXTENSION, tpl_out_fits, tpl_sres, 'SRES', 'Spectral resolution (R=dl/l)'
        MDAP_ADD_TPL_EXTENSION, tpl_out_fits, tpl_soff, 'SOFF', 'Spectral resolution offset'

END


