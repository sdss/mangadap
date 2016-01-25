;+
; NAME:
;       MDAP_COORDINATE_GRID
;
; PURPOSE:
;       Provided with an input list of coordinates, create a grid that
;       covers the full coordinate range at the provided sampling rate,
;       allowing for consideration of an auxilliary set of coordinates
;       and absolute and fractional buffers.
;
;       If both an absolute and fractional buffer are requested, the
;       fractional buffer is applied before the absolute one such that
;
;           Dx = (Dx_range * fbuffer)+abuffer
;
; CALLING SEQUENCE:
;       MDAP_COORDINATE_GRID, xcoo, dx, xs, nx, xaux=xaux, abuffer=abuffer, fbuffer=fbuffer
;
; INPUTS:
;       xcoo dblarr[N]
;               List of coordinates to cover with the grid.
;
;       dx double
;               Sample width of the grid.
;
; OPTIONAL INPUTS:
;       xaux dblarr[A]
;               List of auxilliary coordinates that must also be covered
;               by the grid.
;
;       abuffer double
;               An absolute buffer to add to the total width of the
;               grid.  If fbuffer is also provided, this absolute buffer
;               is applied after the fractional buffer.  Can be
;               negative.
;
;       fbuffer double
;               A fractional buffer used to change the size of the grid.
;               If abuffer is also provided, this fractional buffer is
;               applied before the absolute buffer.  Can be less than
;               one.
;
; OPTIONAL KEYWORDS:
;
; OUTPUT:
;       xs double
;               The starting edge of the grid.
;
;       nx long
;               Number of cells in the grid.
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
;       04 Dec 2014: (KBW) Original implementation
;-
;------------------------------------------------------------------------------

PRO MDAP_COORDINATE_GRID, $
                coo, d, s, n, aux=aux, abuffer=abuffer, fbuffer=fbuffer

        ; Minimum coordinate
        s = min(coo)
        if n_elements(aux) ne 0 then begin
            s_ = min(aux)
            if s gt s_ then $
                s = s_
        endif

        ; Maximum coordinate
        e = max(coo)
        if n_elements(aux) ne 0 then begin
            e_ = max(aux)
            if e lt e_ then $
                e = e_
        endif

        Delta = (e-s)                       ; Full range

        if n_elements(fbuffer) ne 0 then $  ; Add a fractional buffer
            Delta = Delta * fbuffer

        if n_elements(abuffer) ne 0 then $  ; Add an absolute buffer
            Delta = Delta + abuffer

        n = ceil(Delta/d)                   ; Number pixels for the grid
        Delta = n*d                         ; Account for round off
        s = (e+s-Delta)/2.                  ; Offset initial value to account for new width
END



