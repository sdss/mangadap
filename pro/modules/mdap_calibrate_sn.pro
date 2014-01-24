function mdap_calibrate_sn,sn_est,Nspaxl,coeff


   tmp = sn_est^coeff[0]/sqrt(Nspaxl)
   fun = poly(tmp,coeff[1:*])

return,fun
end
