;+
; NAME:
;       MDAP_ERASEVAR
;
; PURPOSE:
;       Free the memory of a variable and force it to become undefined.
;
; CALLING SEQUENCE:
;       MDAP_ERASEVAR, var
;
; INPUTS:
;       var
;               The variable to erase
;
; REVISION HISTORY:
;       30 Sep 2014: (KBW) Original implementation
;       09 Oct 2014: (KBW) First checks that var exists
;-
;------------------------------------------------------------------------------

PRO MDAP_ERASEVAR, $
        var

        if n_elements(var) then $
            tmpvar=size(temporary(var))         ; Erase using temporary()
END

