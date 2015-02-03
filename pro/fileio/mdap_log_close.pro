;+
; NAME:
;       MDAP_LOG_CLOSE
;
; PURPOSE:
;       Close the log file unit, after adding some closing information.
;
; CALLING SEQUENCE:
;       MDAP_LOG_CLOSE, file_unit
;
; INPUTS:
;       file_unit unit
;               Unit to close.
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
;       30 Jan 2015: (KBW) Pulled from manga_dap.pro
;-
;------------------------------------------------------------------------------

PRO MDAP_LOG_CLOSE, $
                file_unit

        printf, file_unit, 'End systime(): ', systime()
        close, file_unit
        free_lun, file_unit
END


