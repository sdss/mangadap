;+
; NAME:
;       MDAP_READ_TEMPLATE_LIBRARY
;
; PURPOSE:
;       Read the full template library into a single array with spectra selected
;       by the array row.
;
; CALLING SEQUENCE:
;       MDAP_READ_TEMPLATE_LIBRARY, library_key, library_files, tpl_flux, $
;                                   tpl_ivar, tpl_mask, tpl_wave, tpl_sres
;
; INPUTS:
;       library_key string
;               Keyword used to set the library being used to fit the data.
;
;               TODO: For now this is separate from the library_files, but they
;               should eventually be the same.  This keyword is used to set the
;               spectral resolution, until this is determined from the library
;               files themselves.
;
;       library_files string
;               String used to search for matching files.  Each file is expected
;               to be a 1D spectrum with the wavelength solution (CRVAL1,
;               CRPIX1, CDELT1) in the header. 
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;
; OUTPUT:
;       tpl_flux dblarr[T][S]
;               Array containing T template spectra, each with S spectral channels.
;
;       tpl_ivar dblarr[T][S]
;               Array with inverse variances in the S spectral channels of the T
;               template spectra.
;
;       tpl_mask dblarr[T][S]
;               Bit mask for template spectra.  Used only to mask spectral
;               regions not covered by all templates.  Value is 0 if pixel is
;               good, value is 1 if it should be masked.
;
;       tpl_wave dblarr[T][S]
;               Wavelength of each spectral channel T for each template spectrum S.
;
;       tpl_sres dblarr[T][S]
;               Spectral resolution (R=lamda/delta lambda) for of each spectral
;               channel T for each template spectrum S.
;
; OPTIONAL OUTPUT:
;
; COMMENTS:
;
; EXAMPLES:
;
; TODO:
;       - Include inverse variances, if available
;
; BUGS:
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;       17 Sep 2014: (KBW) Original implementation
;       10 Oct 2014: (KBW) Changed keyword from 'MARCS' TO 'M11-MARCS'
;-
;------------------------------------------------------------------------------

;--------------------------------------------------------------------------------
; Determine the size of the output array
;       Procedure checks each file in list exists, reads is header, checks that
;       the number of dimensions of each file is 1, and finds the maximum number
;       of spectral channels
PRO MDAP_OUTPUT_TPL_NC, $
                list, nc

        sz=size(list)
        ns = sz[1]                                      ; Number of files to read

        nc = 0

        for i=0,ns-1 do begin
            if file_test(list[i]) eq 0 then begin               ; Check that the file exists
                message, 'ERROR: Cannot read '+list[i]
                retall
            endif

            header = HEADFITS(list[i])                          ; Read the full header
            naxis = SXPAR(header, 'NAXIS', count=nnaxis)        ; Get the number of axes
            if nnaxis eq 0 then begin
                message, 'ERROR: Cannot read header of '+list[i]
                retall
            endif
            if naxis ne 1 then begin
                message, 'ERROR: '+list[i]+' is not one dimension!'
                retall
            endif

            naxis1 = SXPAR(header, 'NAXIS1', count=nnaxis)      ; Get the length of the first axis
            if nnaxis eq 0 then begin
                message, 'ERROR: Cannot number of pixels in '+list[i]
                retall
            endif

            if nc lt naxis1 then $                              ; Get the maximum length
                nc = naxis1

        endfor
END
;--------------------------------------------------------------------------------


;------------------------------------------------------------------------------
;       Return the spectral resolution of the library
;
;       library_key is a keyword
;
;       returned value is the constant spectral resolution in angstroms
;
;       TODO: This will need vast improvements once this operation is more
;       complicated.

FUNCTION MDAP_GET_TEMPLATE_RESOLUTION, $
                library_key
        if library_key eq 'M11-MARCS' then begin
            return, 2.73
        endif else begin
            message, 'Unknown library keyword!'
        endelse
END     
;------------------------------------------------------------------------------

PRO MDAP_READ_TEMPLATE_LIBRARY, $
                library_key, library_files, tpl_flux, tpl_ivar, tpl_mask, tpl_wave, tpl_sres

        ; Get the resolution of the template library (FWHM)
        fwhm_tpl = MDAP_GET_TEMPLATE_RESOLUTION(library_key)

        ; Find the list of files matching the search string
        list = file_search(library_files, count=ntpl)

        MDAP_OUTPUT_TPL_NC, list, nc            ; Get the number of spectral channels for output

        tpl_flux = dblarr(ntpl, nc)             ; Initialize arrays
        tpl_ivar = make_array(ntpl, nc, /double, value=1.0d)    ; TODO: Read values if available
        tpl_mask = dblarr(ntpl, nc)
        tpl_wave = dblarr(ntpl, nc)
        tpl_sres = dblarr(ntpl, nc)

        for i=0,ntpl-1 do begin

            MDAP_READ_1DSPEC_FITS, list[i], flux, wave  ; Read the spectral data

            sz=size(flux)

            tpl_flux[i,0L:sz[1]-1] = flux[*]            ; Flux
            tpl_wave[i,0L:sz[1]-1] = wave[*]            ; Wavelength
            tpl_sres[i,0L:sz[1]-1] = wave[*]/fwhm_tpl   ; Spectral resolution = lambda/delta lambda

            if sz[1] lt nc then begin
                dw=tpl_wave[i,1]-tpl_wave[i,0]  ; Pixel size
                add_wave = (dindgen(nc-sz[1])+1.0d)*dw+tplwave[i,sz[1]-1]
                tpl_wave[i,sz[1]:nc-1] = add_wave[*]    ; Add the wavelengths for the other pixels

                tpl_mask[i,sz[1]:nc-1] = 1.0            ; Mask pixels that have no data
            endif

        endfor

END

