;+
; NAME:
;       MDAP_READ_1DSPEC_FITS
;
; PURPOSE:
;       Read a one-dimensional spectrum from a fits file, including the
;       wavelength coordinates.
;
; CALLING SEQUENCE:
;       MDAP_READ_1DSPEC_FITS, file, flux, wave
;
; INPUTS:
;       file string
;               Fits file to read.
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;
; OUTPUT:
;       flux dblarr[C]
;               Flux (or flux density) values for each of C spectral channels.
;
;       wave dblarr[C]
;               Wavelength of each pixel.  If header does not have the correct
;               keywords, the returned vector is the just the pixel number,
;               1...C.
;
; OPTIONAL OUTPUT:
;
; COMMENTS:
;
; EXAMPLES:
;
; TODO:
;       - Determine if the spectrum is linearly or logarithmically binned?
;       - Allow for WCS keyword names to be input
;       - Handle spectra that are more than one dimension
;
; BUGS:
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;       17 Sep 2014: (KBW) Original implementation
;-
;------------------------------------------------------------------------------

PRO MDAP_READ_1DSPEC_FITS, $
                file, flux, wave

        if file_test(file) eq 0 then begin                      ; Test that the file exists
            message, 'ERROR: Cannot read '+file
            retall                                              ; FAULT
        endif

        flux=READFITS(file)                                     ; Read the data

        sz=size(flux)                                           ; Get the number of dimensions
        if sz[0] ne 1 then begin                                ; Is it one-dimensional?
            message, file+' has '+MDAP_STC(sz[0],/integer)+' dimensions!'       ; Nope
            retall                                              ; FAULT
        endif

        nc = sz[1]                                              ; Number of spectral channels
        wave = (findgen(nc+1))[1:*]                             ; Initialize wave to the pixel #

        header=HEADFITS(file)                                   ; Read the header
        crval=SXPAR(header, 'CRVAL1', count=ncrval)             ; Try to read the CRVAL1 ...
        crpix=SXPAR(header, 'CRPIX1', count=ncrpix)             ; ... CRPIX1
        cdelt=SXPAR(header, 'CDELT1', count=ncdelt)             ; ... and CDELT1

        if ncrval eq 0 or ncrpix eq 0 or ncdelt eq 0 then begin ; Is there a known wavelength sol.?
            print, 'WARNING: Could not determine wavelength solution.  Returning pixel coordinates.'
            return
        endif

        wave = (wave - crpix) * cdelt + crval                   ; Set the wavelengh of each pixel
                                                                ; TODO: Assumes linear sampling
END


