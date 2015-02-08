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

PRO MDAP_REWRITE_DRP_HEADERS, $
                ifile, ofile
        print, ifile
        print, ofile

        if FILE_TEST(ofile) eq 1 then begin
            print, 'Removing existing output file: ', ofile
            FILE_DELETE, ofile
        endif

        gzi = STRPOS(ifile, '.gz')
;       if gzi ge 0 and STRPOS(ofile, '.gz') lt 0 then $
        if gzi ge 0 && STRPOS(ofile, '.gz') lt 0 then $
            ofile += '.gz'

        FILE_COPY, ifile, ofile

        gzi = STRPOS(ofile, '.gz')
        if gzi ge 0 then begin
            oofile = STRMID(ofile, 0, gzi)
            cmnd = 'gunzip '+ofile
            print, cmnd
            SPAWN, cmnd
        endif else $
            oofile = ofile

        print, oofile
        FITS_INFO, oofile, N_ext=next, extname=extname, /silent

        print, 'Number of extensions:', next
        print, 'Removing SIMPLE keyword in all non-Primary headers'

        for i=0,next-1 do begin
            extname[i] = strcompress(extname[i], /remove_all)
            print, 'Fixing extention '+extname[i]

            ext=i+1
            data=readfits(oofile, hdr, exten=ext)

            print, fxpar(hdr, 'SIMPLE')

            SXDELPAR, hdr, 'SIMPLE'
            MODFITS, oofile, hata, hdr, exten_no=ext
        endfor

        print, 'DONE with fixes'

        if gzi gt 0 then begin
            cmnd = 'gzip '+oofile
            print, cmnd
            SPAWN, cmnd
        endif
END


