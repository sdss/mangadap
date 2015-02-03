;+
; NAME:
;       MDAP_CHECK_FILE
;
; PURPOSE:
;       Check that a file (or list of files) exist.  If the keyword /search is
;       provided, the input list is used to search for files instead of being a
;       file name itself.
;
; CALLING SEQUENCE:
;       result=MDAP_CHECK_FILE(list, /search)
;
; INPUTS:
;       list string[]
;               An array of files to check for using FILE_TEST() or search for
;               using  FILE_SEARCH().
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;       /search
;               Strings in list are not files but search strings for used to get
;               a list of files.
;
; OUTPUT:
;       result integer
;               Success of the search: 0 for success, 1 for failure.  Success is
;               returned if list has no entries or all strings in list result in
;               at least one existing file; failure returned otherwise.
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
;       09 Oct 2014: (KBW) Original implementation
;-
;------------------------------------------------------------------------------

FUNCTION MDAP_CHECK_FILE_EXISTS, $
                list, search=search

        n_list = n_elements(list)
        if n_list eq 0 then $
            return, 0                           ; No searches performed, return success

        for i=0,n_list-1 do begin
            if keyword_set(search) then begin
                found_list = file_search(list[i], count=nfiles)
            endif else $
                nfiles = file_test(list[i])

            if nfiles eq 0 then $
                return, 1                       ; No files resulted from the search
        endfor
                
        return, 0                               ; All searches resulted in a list of files
END

