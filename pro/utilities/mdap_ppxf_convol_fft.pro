;----------------------------------------------------------------------------
FUNCTION mdap_ppxf_convol_fft, f, k
compile_opt idl2, hidden

nf = n_elements(f)
nk = n_elements(k)
n = 2L^ceil(alog(nf+nk/2)/alog(2))
f1 = dblarr(n)
k1 = f1
f1[0] = f
k1[0] = rotate(k,2)
k1 = shift(k1,-(nk-1)/2)
con = n*double(fft(fft(f1,-1)*fft(k1,-1),1))

return, con[0:nf-1]
END
;----------------------------------------------------------------------------
