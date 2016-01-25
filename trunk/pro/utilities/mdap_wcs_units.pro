;+
; NAME:
;       MDAP_WCS_UNITS
;
; PURPOSE:
;       Determine the units for the WCS coordinates in the provided header.
;
; CALLING SEQUENCE:
;       MDAP_WCS_UNITS, header, unit
;
; INPUTS:
;       header hdu
;               Header read from the DRP-produced fits file.
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;
; OUTPUT:
;       unit string
;               String representation of the units as determined by the CUNIT1
;               and CUNIT2 header keywords.
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
;       15 Sep 2014: (KBW) Original Implementation
;-
;------------------------------------------------------------------------------

PRO MDAP_WCS_UNITS, $
        header, unit

        unit=sxpar(header, 'CUNIT1')                    ; Pull the units from the header
        if !ERR eq -1 then begin                        ; Keyword not found!
            print, 'CUNIT1 unavailable in header.  Assuming units are degrees.'
            unit='degrees'                              ; Assume its degrees
        endif else begin
            unit=strcompress(unit, /remove_all)         ; Make sure there aren't any spaces
            if strcompress(sxpar(header, 'CUNIT2'), /remove_all) ne unit then $
                print, 'CUNIT1 != CUNIT2!  Assuming units are the same for both axes.'
        endelse

END

