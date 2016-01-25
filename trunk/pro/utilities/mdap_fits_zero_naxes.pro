;+
; NAME:
;       MDAP_FITS_ZERO_NAXES
;
; PURPOSE:
;       Remove NAXIS and up to 10 NAXISn keywords, and then add NAXIS=0 to the
;       header
;
; CALLING SEQUENCE:
;       MDAP_FITS_ZERO_NAXES, header
;
; INPUTS:
;       header HDU
;               Header to alter
;
; PROCEDURES CALLED:
;       SXDELPAR
;
; INTERNAL SUPPORT ROUTINES:
;       MDAP_FITS_DEL_NAXES
;
; REVISION HISTORY:
;
;-
;------------------------------------------------------------------------------

;-------------------------------------------------------------------------------
PRO MDAP_FITS_DEL_NAXES, $
                header
        ; TODO: Do this, or actually read NAXIS and delete that number of NAXISn keywords?
        SXDELPAR, header, 'NAXIS'
        for i=1,10 do $
            SXDELPAR, header, 'NAXIS'+strcompress(string(i),/remove_all)
END

;-------------------------------------------------------------------------------
; Delete NAXIS and up to 10 NAXISn keywords then add NAXIS = 0 to the header
PRO MDAP_FITS_ZERO_NAXES, $
                header
        MDAP_FITS_DEL_NAXES, header
        SXADDPAR, header, 'NAXIS', 0, 'Primary HDU', after='BITPIX'
END

