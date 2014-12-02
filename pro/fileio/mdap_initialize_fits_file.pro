;+
; NAME:
;       MDAP_INITIALIZE_FITS_FILE
;
; PURPOSE:
;       If the file does not exist, create a bare-bones primary header that
;       accepts extensions and write the file.
;
; CALLING SEQUENCE:
;       MDAP_INITIALIZE_FITS_FILE, file
;
; INPUTS:
;       file string
;               Name of the file
;
; PROCEDURES CALLED:
;       MKHDR
;       MDAP_FITS_ZERO_NAXES
;       MDAP_FITS_HEADER_ADD_DATE
;       WRITEFITS
;
; REVISION HISTORY:
;       23 Oct 2014: (KBW) Original Implementation
;-
;-------------------------------------------------------------------------------

PRO MDAP_INITIALIZE_FITS_FILE, $
                file
        if file_test(file) eq 1 then $          ; File already exists!
            return

        MKHDR, header, dblarr(1), /extend       ; Create a bare-bones header with extensions
        MDAP_FITS_ZERO_NAXES, header            ; Get rid of the axes
        MDAP_FITS_HEADER_ADD_DATE, header       ; Add the DAP date
        WRITEFITS, file, 0, header              ; Write it
END

