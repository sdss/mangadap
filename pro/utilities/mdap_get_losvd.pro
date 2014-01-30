pro mdap_get_losvd,sol,velscale,losvd
   dx = ceil(5d*sol[1]/velscale) 
   n = 2*dx + 1
   x = range(dx,-dx,n)          ; Evaluate the Gaussian using steps of 1/factor pixel
   losvd = dblarr(n)
   w = (x )/(sol[1]/velscale)
   w2 = w*w
   losvd = exp(-0.5d*w2)/(sqrt(2d*!dpi)*sol[1]/velscale) ; Normalized total(Gaussian)=1
   poly = 1d + sol[2]/Sqrt(3d)*(w*(2d*w2-3d)) $          ; H3
          + sol[3]/Sqrt(24d)*(w2*(4d*w2-12d)+3d)         ; H4
   losvd = losvd*poly
end
