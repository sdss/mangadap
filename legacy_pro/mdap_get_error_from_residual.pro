pro mdap_get_error_from_residual,residuals,galaxy,errors


   errors=residuals*0.
   sz=size(galaxy)
   nbins=sz[2]/20
   fake_x = findgen(sz[2])
   rms_per_bin = fltarr(nbins)
   bin_center = fltarr(nbins)
   bins = mdap_range(0,sz[2],nbins+1)
   for i = 0, sz[1]-1 do begin
      res=residuals[i,*]
      for kk = 0, nbins-1 do begin
         indici = where(fake_x ge bins[kk] and fake_x lt bins[kk+1])
         rms_per_bin[kk] = stddev(res[indici])  
         bin_center[kk] = (bins[kk+1]+ bins[kk])*0.5
      endfor
      rms = interpol(rms_per_bin,bin_center,findgen(sz[2]))
      errors[i,*] = temporary(rms)
   endfor

end
