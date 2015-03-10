;+
; NAME:
;
; PURPOSE:
;
; CALLING SEQUENCE:
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
; BUGS:
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;
;-
;------------------------------------------------------------------------------

PRO REBIN_STELIB

        lib_files = 'M11_STELIB_*z_?.??.fits'
        list = file_search(lib_files, count=ntpl)

        print, ntpl

        ; Bin the spectra by 2 pixels so that the resolution is ~3.1
        ; pixels instead of ~6.2 pixels

;       ntpl=1
        for i=0,ntpl-1 do begin

            ; Get the output file name
            ofile = strmid(list[i], 0, strpos(list[i], '.fits'))+'_s.fits'

            ; Read the data and header
            data = readfits(list[i], hdr)

            ; Read the WCS coords
            crpix = fxpar(hdr, 'CRPIX1')
            crval = fxpar(hdr, 'CRVAL1')
            cdelt = fxpar(hdr, 'CDELT1')

            npix = n_elements(data)

            wave = crval + (dindgen(npix)+1-crpix)*cdelt

            ; Assumes linear sampling!
            l0 = crval + (1-crpix)*cdelt            ; Wavelength of first pixel
            l1 = l0+cdelt
            
            new_crpix = 1
            new_crval = (l0+l1)/2.
            new_cdelt = 2.0*cdelt
            new_npix = npix/2

            new_wave = new_crval + (dindgen(new_npix)+1-new_crpix)*new_cdelt

            new_data = interpol(data, wave, new_wave, /spline)

            FXADDPAR, hdr, 'CRPIX1', new_crpix
            FXADDPAR, hdr, 'CRVAL1', new_crval
            FXADDPAR, hdr, 'CDELT1', new_cdelt

            metallicity = fxpar(hdr, 'Z')
;           print, metallicity
;           print, n_elements(new_wave)

;           if metallicity lt 0.019 or metallicity gt 0.021 then begin
;               print, metallicity
                indx=where(new_wave lt 7900)
;               print, n_elements(new_wave[indx])
                new_data = new_data[indx]
;           endif

            WRITEFITS, ofile, new_data, hdr
        endfor
end

