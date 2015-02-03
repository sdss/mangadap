;+
; NAME:
;       MDAP_SET_TPL_LIB_OUTPUT_FILE
;
; PURPOSE:
;
; CALLING SEQUENCE:
;       MDAP_SET_TPL_LIB_OUTPUT_FILE, file_root, library_key, abs_line_key=abs_line_key
;
; INPUTS:
;       file_root string
;               Root name for all DAP output.
;
;       library_key string
;               String keyword representation of the template library.
;
; OPTIONAL INPUTS:
;       abs_line_key string
;               String keyword representation of the spectral-index
;               library.
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
;       30 Jan 2015: (KBW) Pulled from manga_dap.pro
;-
;------------------------------------------------------------------------------

FUNCTION MDAP_SET_TPL_LIB_OUTPUT_FILE, $
                file_root, library_key, abs_line_key=abs_line_key
        if n_elements(abs_line_key) eq 0 then $
            return, file_root+library_key+'.fits'

        return, file_root+library_key+'_'+abs_line_key+'.fits'
END

