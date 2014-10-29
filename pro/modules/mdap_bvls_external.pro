;------------------------------------------------------------------------------------
PRO mdap_bvls_external, A, b, bnd, x, external_library, IERR=ierr
;
; IDL Interface to the compiled Fortran 90 BVLS routine
; By Michele Cappellari, Leiden, 17 April 2002
;
;print, 'external library used'
s = size(A)
IF s[0] ne 2 then message, 'A must have two dimensions'
if s[3] ne 5 then message, 'A must be DOUBLE'
m = long(s[1])
n = long(s[2])
IF (N LT 2 OR M LT 2 OR N_ELEMENTS(B) NE M $
 OR (SIZE(BND))[1] NE 2 OR (SIZE(BND))[2] NE N) THEN $
   MESSAGE, 'Wrong input arrays size'

a = double(a)
b = double(b)
x = dblarr(n)
rnorm = 0d0
nsetp = long(0)
w = dblarr(n)
index = lonarr(n) ; This corresponds to Fortran INTEGER*4
ierr = long(0)

; WARNING: Use the full path when calling the external library
;so_file=file_search('./','bvls.so',/FULLY_QUALIFY_PATH)

tmp = CALL_EXTERNAL(external_library+'bvls.dylib', 'bvls', $
    	A, M, N, B, BND, X, RNORM, NSETP, W, INDEX, IERR)

END
