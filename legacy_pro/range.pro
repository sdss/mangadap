;#############################################################################
;+
; NAME:
;   RANGE
;
; PURPOSE:
;   RANGE(x1,x2) = x1,x1+1,...,x2. In this case x1, x2 should be integers.
;   RANGE(x1,x2,n) = x1,x1+dx,...,x2 with N integer the result has length N.
;   The result will have the type of x1, but if three parameters are used
;   the result will be at least of type float.
;
; CALLING SEQUENCE:
;   a = RANGE(x1,x2,n,/OPEN,/LOG,/ASINH,SOFT=soft)
; 
; KEYWORDS:
;   OPEN: Set this keyword to exclude the extremes of the interval.
;   LOG: Set this keyword to have the values logarithmically spaced
;   ASINH: Set this keyword to have the values spaced according to an
;     inverse hyperbolic sine function (Lupton et al. 1999, AJ, 118, 1406).
;   SOFT: When using the /ASINH keyword, the optional keyword 
;     0 < SOFT < 1 (default SOFT=0.3) can be assigned to the fractional 
;     value which defines the transition between a nearly linear spacing 
;     (x < soft) and a nearly logarithmic spacing (x > soft).
;
; EXAMPLES:
;   IDL> print, range(0,1,5)
;       0.000000  0.250000  0.500000  0.750000  1.00000
;   IDL> print, range(0,1,5,/OPEN)
;       0.100000  0.300000  0.500000  0.700000  0.90000
;   IDL> print, range(1,1e4,5,/LOG)
;       1.00000   10.0000   100.000   1000.00   10000.0
;   IDL> print, range(1,1e4,5,/LOG)
;
; MODIFICATION HISTORY:
;   V1.0: Michele Cappellari, Leiden, 16 October 2001
;   V1.1: added /OPEN keyword, MC, Leiden, 9 March 2003
;   V1.2: added /LOG keyword, MC, Oxford, 28 November 2008
;   V1.3: added /ASINH and SOFT keywords. Updeted documentation.
;     MC, Oxford, 28 September 2009
;-
;----------------------------------------------------------------------
FUNCTION range, xx1, xx2, n, OPEN=open, LOG=log, ASINH=asinh, SOFT=soft
compile_opt idl2
on_error, 2

if keyword_set(log) then begin
    x1 = alog(xx1)
    x2 = alog(xx2)
endif else begin    
    x1 = xx1
    x2 = xx2
endelse

t = size(x1,/TYPE)
if keyword_set(asinh) then begin 
    if n_elements(soft) gt 0 then s = soft else s = 0.3d
    v = indgen(n, TYPE=t)/(n-1.0)
    sh = sinh(v/s)
    v = x1 + (x2-x1)/sh[n-1]*sh
endif else begin
    if keyword_set(open) then $ ; Open interval: excluding the extremes
        v = x1 + (x2-x1)*(0.5+indgen(n, TYPE=t))/n $
    else case n_params() of
        2: if x1 lt x2 then $
                v = x1 + indgen(x2-x1+1, TYPE=t) $
            else $
                v = x1 - indgen(x1-x2+1, TYPE=t)
        3: v = x1 + (x2-x1)/(n-1.0)*indgen(n, TYPE=t)
        else: message, '2 or 3 parameters are needed'
    endcase
    if keyword_set(log) then v = exp(v)
endelse    

return, v
END
;----------------------------------------------------------------------
