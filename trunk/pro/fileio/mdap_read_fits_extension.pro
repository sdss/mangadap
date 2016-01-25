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

PRO MDAP_READ_FITS_EXTENSION, $
                file, exten, data, header=header, to_double=to_double

        errmsg = ''                     ; Define the error message so it will be collected

        ; Open file at 'FLUX' extension
        fits_unit = FXPOSIT(file, exten, errmsg=errmsg)
        if fits_unit eq -1 then $
            message, 'Could not read extension '+exten+': '+errmsg

        if n_elements(header) ne 0 then begin
            data = READFITS(fits_unit, header)      ; Read the header if requested
        endif else $
            data = READFITS(fits_unit)              ; Read the data only

        ; Convert it to a double if requested
        if keyword_set(to_double) then $
            data = double(data)

        print, 'Successfully read '+exten+' extension'

        ; Free the LUN
        free_lun, fits_unit
END

