;+
; NAME:
;       MDAP_SELECTED_KEYWORD_INDEX
;
; PURPOSE:
;       Given a set of keywords and a set of selected keywords, return a
;       vector with the selected keyword into the index of the keyword
;       in the list.
;
; CALLING SEQUENCE:
;       result = MDAP_SELECTED_KEYWORD_INDEX(keys, selected_key)
;
; INPUTS:
;       keys strarr[]
;               List of valid keywords.
;
;       selected_key strarr[]
;               List of selected keywords.
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;
; OUTPUT:
;       Result is a intarr[] with the list of indices from the provided
;       list of keywords.  The default of 0 is returned if the
;       selected_key has no length; if a key with non-zero length is not
;       found, the function throws a message.
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
;       10 Feb 2015: Original implementation by K. Westfall (KBW)
;-
;------------------------------------------------------------------------------

FUNCTION MDAP_SELECTED_KEYWORD_INDEX, $
                keys, selected_key

        nkeys = n_elements(keys)                    ; Number of available keys
        nselected = n_elements(selected_key)        ; NUmber to select

        ; Make sure there's something to do!
        if nkeys eq 0 || nselected eq 0 then $
            message, 'Incorrect usage of MDAP_SELECTED_KEYWORD_INDEX'

        ; Default output is that the index is zero
        index_for_analysis = intarr(nselected)

        ; Search through all the selected keyword for the appropriate index
        for i=0,nselected-1 do begin

            ; If the keyword is not set (has no length) just use the
            ; default index (0)
            if strlen(selected_key) eq 0 then $
                continue

            ; Try to find the key 
            found = 0
            for j=0,nkeys-1 do begin
                if selected_key[i] eq keys[j] then begin
                    index_for_analysis[i] = j
                    print, 'Selected '+keys[j]+'; index='+MDAP_STC(j)
                    found = 1
                    break
                endif
            endfor

            ; Key was not found, force a fault
            if found eq 0 then begin
                print, 'OPTIONS ARE:'
                for j=0,nkeys-1 do $
                    print, '    '+keys[j]
                message, 'Could not find '+selected_key[i]+' among the valid keywords!'
            endif
        endfor

        ; Return the appropriate vector of indices
        return, index_for_analysis

END

