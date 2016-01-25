;+
; NAME:
;       MDAP_DRP_FILETYPE
;
; PURPOSE:
;       Determine and/or test the validity of a DRP-produced file type.
;
; CALLING SEQUENCE:
;       MDAP_DRP_FILETYPE, header, type
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
;       type string
;               String representation of the file type.  If provided on input,
;               it is checked against the expected file types, either 'RSS' or
;               'CUBE'.  If not provided or if the input value is invalid, the
;               type is determined by the keyword 'NAXIS' in the input header:
;               for 'RSS' files, NAXIS=2; for 'CUBE' files, NAXIS=3.
;
; OPTIONAL OUTPUT:
;
; COMMENTS:
;
; EXAMPLES:
;
; TODO:
;       - Change to using splog
;
; BUGS:
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;       09 Sep 2014: (KBW) Original implementation
;       12 Sep 2014: (KBW) Convert input to be the fits header, not the flux
;                          data
;-
;------------------------------------------------------------------------------

PRO MDAP_DRP_FILETYPE, $
                header, type

        ndim = SXPAR(header,'NAXIS')                    ; Get the number of dimensions
;       ndim=size(flux, /n_dimensions)

;       if ndim ne 2 and ndim ne 3 then begin
        if ndim ne 2 && ndim ne 3 then begin
            message, 'ERROR: Cannot determine file type because dimensionality is not 2 or 3.'
            retall
        endif

        defined_type = 1                                ; Assume type is defined
        if n_elements(type) ne 0 then begin             ; Check the type makes sense
;           if type ne 'RSS' and type ne 'CUBE' then begin
            if type ne 'RSS' && type ne 'CUBE' then begin
                print, 'ERROR: Undefined type!  Must be either RSS or CUBE.  Obtained from header.'
                defined_type = 0
            endif
        endif

        ; Determine the type because it was not provided or undefined
;       if n_elements(type) eq 0 or defined_type eq 0 then begin
        if n_elements(type) eq 0 || defined_type eq 0 then begin
            if ndim eq 2 then begin
                type='RSS'
            endif else $
                type='CUBE'
        endif

        ; Finally check that the input type matches the fits file dimensionality
;       if type eq 'CUBE' and ndim eq 2 then begin
        if type eq 'CUBE' && ndim eq 2 then begin
            print, 'WARNING: Input was type=CUBE, but dimensionality suggests RSS. Assuming RSS.'
            type='RSS'
        endif
;       if type eq 'RSS' and ndim eq 3 then begin
        if type eq 'RSS' && ndim eq 3 then begin
            print, 'WARNING: Input was type=RSS, but dimensionality suggests CUBE. Assuming CUBE.'
            type='CUBE'
        endif

END



