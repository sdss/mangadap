;+
; NAME:
;       MDAP_INSERT_FLAGGED
;
; PURPOSE:
;       Insert -1 values for indices not returned by a where() statement into a
;       vector produced using an operation based on the where() vector.
;
; CALLING SEQUENCE:
;       MDAP_INSERT_FLAGGED, where_result, where_produced, original_vector_size
;
; INPUTS:
;       where_result intarr[]
;               Index array resulting from a call to where.
;
;       where_produced intarr[]
;               Array resulting from some operation based on the where result.
;               Replaced by vector with null results set to -1.
;
;       original_vector_size int
;               Original length of the vector used to produce where_result
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
;       15 Sep 2014: (KBW) Original implementation
;-
;------------------------------------------------------------------------------

PRO MDAP_INSERT_FLAGGED, $
                where_result, where_produced, original_vector_size

        where_produced_=intarr(original_vector_size)
        n=n_elements(where_result)

        j=0
        for i=0,original_vector_size-1 do begin
            if j ge n then begin
                where_produced_[i] = -1
                continue
            endif

            if where_result[j] eq i then begin
                where_produced_[i] = where_produced[j]
                ++j
            endif else $
                where_produced_[i] = -1
        endfor

        where_produced = where_produced_
END

