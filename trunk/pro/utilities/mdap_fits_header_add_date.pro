;+
; NAME:
;       MDAP_FITS_HEADER_ADD_DATE
;
; PURPOSE:
;       Add the date (in fits form) to the header.
;
; CALLING SEQUENCE:
;       MDAP_FITS_HEADER_ADD_DATE, header, /modified
;
; INPUTS:
;       header HDU
;               Header to modify
;
; OPTIONAL KEYWORDS:
;       /modified
;               If provided, then add/change the keyword DATEMOD, otherwise
;               add/change DATEDAP.
;
; PROCEDURES CALLED:
;       SXADDPAR
;       DATE_CONV()
;
; REVISION HISTORY:
;       23 Oct 2014: (KBW) Original implementation
;-
;-------------------------------------------------------------------------------
PRO MDAP_FITS_HEADER_ADD_DATE, $
                header, modified=modified
        if keyword_set(modified) then begin
            SXADDPAR, header, 'DATEMOD', DATE_CONV( systime(/julian,/utc), 'F'), $
                      'UTC date last modified by DAP'
        endif else begin
            SXADDPAR, header, 'DATEDAP', DATE_CONV( systime(/julian,/utc), 'F'), $
                      'UTC date DAP created file'
        endelse
END

