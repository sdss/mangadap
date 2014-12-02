;+
; NAME:
;       MDAP_DRP_CHECK_TYPE
;
; PURPOSE:
;       Check that the input file type is the actual file type; if not, write a
;       warning message and fix the file type.
;
; CALLING SEQUENCE:
;       MDAP_DRP_CHECK_TYPE, file, type
;
; INPUTS:
;       file string
;               File name of DRP-produced fits file to check.
;
;       type string
;               Fits file type.  Must be either 'RSS' or 'CUBE'.  If not
;               correct, replaced with correct value on output.
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
;       MDAP_DRP_FILETYPE
;       HEADFITS()
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;       09 Oct 2014: (KBW) Original implementation
;-
;------------------------------------------------------------------------------

PRO MDAP_DRP_CHECK_FILETYPE, $
                file, type

        header = headfits(file, exten=1)                ; Read the header
        MDAP_DRP_FILETYPE, header, type                 ; Get/check the file type

END

