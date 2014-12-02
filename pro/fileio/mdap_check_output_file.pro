;+
; NAME:
;       MDAP_CHECK_OUTPUT_FILE
;
; PURPOSE:
;       Check that the provided file has the correct extension names and
;       the correct columns in each binary table.
;
; CALLING SEQUENCE:
;       MDAP_CHECK_OUTPUT_FILE, file
;
; INPUTS:
;       file string
;               Fits file name
;
; OUTPUT:
;       Returned value is 1 for true (file checks out as valid), 0 for
;       false.
;
; PROCEDURES CALLED:
;       @mdap_set_output_file_cols
;       FITS_INFO
;       FXBCOLNUM
;
; INTERNAL SUPPORT ROUTINES:
;       MDAP_FIND_EXTENSION_NAME()
;       MDAP_CHECK_BINTABLE_COLUMNS()
;
; REVISION HISTORY:
;       01 Dec 2014: (KBW) Original Implementation (comments added)
;-
;------------------------------------------------------------------------------

;-------------------------------------------------------------------------------
; Find the extension name in the list of available extensions
FUNCTION MDAP_FIND_EXTENSION_NAME, $
                extname, ext
        next = n_elements(extname)              ; Number of extensions
        for i=0,next-1 do begin                 ; Search through all extensions
            if extname[i] eq ext then $
                return, 1                       ; Extension found so return true
        endfor
        print, ext+' extension not found!'      ; Warn the user
        return, 0                               ; Extension name not found
END

@mdap_set_output_file_cols

;-------------------------------------------------------------------------------
; Check that the table has all the necessary columns
FUNCTION MDAP_CHECK_BINTABLE_COLUMNS, $
                file, exten, col_list

        fxbopen, unit, file, exten

        ncol = n_elements(col_list)             ; Number of column names
        errmsg = ''                             ; Suppress error messages
        for i=0,ncol-1 do begin
            if FXBCOLNUM(unit, col_list[i]) eq 0 then begin
                print, col_list[i]+' unavailable in extension '+exten+'!'
                fxbclose, unit                  ; Close up ...
                free_lun, unit                  ; ... and free the LUN
                return, 0                       ; All columns not found
            endif
        endfor

        fxbclose, unit                          ; Close up ...
        free_lun, unit                          ; ... and free the LUN
        return, 1                               ; All columns found
END

;-------------------------------------------------------------------------------
; Check that the file at least has the necessary properties of a DAP-produced
; file
FUNCTION MDAP_CHECK_OUTPUT_FILE, $
                file

        ; Get the number of extensions and the extension names
        FITS_INFO, file, N_ext=next, extname=extname, /silent
        print, next
        extname = extname[1:next]
        for i=0,next-1 do $
            extname[i] = strcompress(extname[i], /remove_all)

        ; Ensure that the file has the correct number of extensions
        ; TODO: Not necessary if extensions are read based on their name, not their index
;       if next ne 17 then begin
;           print, file + ' has '+MDAP_STC(next)+' extensions!'
;           return, 0
;       endif

        ; Check that the extensions have the correct names
        if MDAP_FIND_EXTENSION_NAME(extname, 'DRPS') eq 0 then  return, 0
        if MDAP_FIND_EXTENSION_NAME(extname, 'BINS') eq 0 then  return, 0
        if MDAP_FIND_EXTENSION_NAME(extname, 'WAVE') eq 0 then  return, 0
        if MDAP_FIND_EXTENSION_NAME(extname, 'SRES') eq 0 then  return, 0
        if MDAP_FIND_EXTENSION_NAME(extname, 'FLUX') eq 0 then  return, 0
        if MDAP_FIND_EXTENSION_NAME(extname, 'IVAR') eq 0 then  return, 0
        if MDAP_FIND_EXTENSION_NAME(extname, 'MASK') eq 0 then  return, 0
        if MDAP_FIND_EXTENSION_NAME(extname, 'ELPAR') eq 0 then return, 0
        if MDAP_FIND_EXTENSION_NAME(extname, 'STFIT') eq 0 then return, 0
        if MDAP_FIND_EXTENSION_NAME(extname, 'SMSK') eq 0 then  return, 0
        if MDAP_FIND_EXTENSION_NAME(extname, 'SMOD') eq 0 then  return, 0
        if MDAP_FIND_EXTENSION_NAME(extname, 'SGFIT') eq 0 then return, 0
        if MDAP_FIND_EXTENSION_NAME(extname, 'SGMSK') eq 0 then return, 0
        if MDAP_FIND_EXTENSION_NAME(extname, 'SGMOD') eq 0 then return, 0
        if MDAP_FIND_EXTENSION_NAME(extname, 'ELMOD') eq 0 then return, 0
        if MDAP_FIND_EXTENSION_NAME(extname, 'SIPAR') eq 0 then return, 0
        if MDAP_FIND_EXTENSION_NAME(extname, 'SINDX') eq 0 then return, 0
        if MDAP_FIND_EXTENSION_NAME(extname, 'SIWAVE') eq 0 then return, 0
        if MDAP_FIND_EXTENSION_NAME(extname, 'SIFLUX') eq 0 then return, 0
        if MDAP_FIND_EXTENSION_NAME(extname, 'SIIVAR') eq 0 then return, 0
        if MDAP_FIND_EXTENSION_NAME(extname, 'SIMASK') eq 0 then return, 0
        if MDAP_FIND_EXTENSION_NAME(extname, 'SIOTPL') eq 0 then return, 0
        if MDAP_FIND_EXTENSION_NAME(extname, 'SIBOTPL') eq 0 then return, 0
;       if MDAP_FIND_EXTENSION_NAME(extname, 'OTPL') eq 0 then  return, 0
;       if MDAP_FIND_EXTENSION_NAME(extname, 'BOTPL') eq 0 then return, 0

        ; Check that the binary tables have the correct columns
        cols = MDAP_SET_DRPS_COLS()
        if MDAP_CHECK_BINTABLE_COLUMNS(file, 'DRPS', cols) eq 0 then return, 0

        cols = MDAP_SET_BINS_COLS()
        if MDAP_CHECK_BINTABLE_COLUMNS(file, 'BINS', cols) eq 0 then return, 0

        cols = MDAP_SET_ELPAR_COLS()
        if MDAP_CHECK_BINTABLE_COLUMNS(file, 'ELPAR', cols) eq 0 then return, 0

        cols = MDAP_SET_STFIT_COLS()
        if MDAP_CHECK_BINTABLE_COLUMNS(file, 'STFIT', cols) eq 0 then return, 0

        cols = MDAP_SET_SGFIT_COLS()
        if MDAP_CHECK_BINTABLE_COLUMNS(file, 'SGFIT', cols) eq 0 then return, 0

        cols = MDAP_SET_SIPAR_COLS()
        if MDAP_CHECK_BINTABLE_COLUMNS(file, 'SIPAR', cols) eq 0 then return, 0

        cols = MDAP_SET_SINDX_COLS()
        if MDAP_CHECK_BINTABLE_COLUMNS(file, 'SINDX', cols) eq 0 then return, 0

        return, 1                               ; File satisfied all checks
END


