;+
; NAME:
;
; PURPOSE:
;
; CALLING SEQUENCE:
;
; INPUTS:
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
;
;-
;------------------------------------------------------------------------------


; Pulled from old version of spatial_binning
   test = file_test(user_bin_map)
   if test eq 0 then begin
      print,user_bin_map,' not found. User input ignored'
      goto,user_map_failure
   endif

   user_bin_map = readfits(user_bin_map,header_user,/silent)
   x_user = sxpar(header_user,'CRVAL1',COUNT=cn_CRVAL1) + sxpar(header_user,'CDELT1',COUNT=cn_CDELT1)*(findgen(sxpar(header_user,'NAXIS1',COUNT=cn_NAXIS1))- sxpar(header_user,'CRPIX1',COUNT=cn_CRPIX1)+1)
   y_user = sxpar(header_user,'CRVAL2',COUNT=cn_CRVAL2) + sxpar(header_user,'CDELT2',COUNT=cn_CDELT2)*(findgen(sxpar(header_user,'NAXIS2',COUNT=cn_NAXIS2))- sxpar(header_user,'CRPIX2',COUNT=cn_CRPIX2)+1)
   x2d_user=x_user # (y_user *0 + 1)
   y2d_user=(x_user*0+1)#y_user

   if cn_CRVAL1 eq 0 or cn_CDELT1 eq 0 or cn_NAXIS1 eq 0 or cn_CRVAL2 eq 0 or cn_CDELT2 eq 0 or cn_NAXIS2 eq 0 or cn_CRPIX1 eq 0 or cn_CRPIX2 eq 0 then begin
      print, user_bin_map,' does not have proper header keyword. CRVAL1, CRVAL2, CDELT1, CDELT2, NAXIS1, NAXIS2, CRPIX1, and CRIX2 are required. User input ignored'
      goto,user_map_failure
   endif

   ;x2d_user
   ;y2d_user
   binNum_=fltarr(n_elements(x2d[ind_good]))*0.-1.
   user_input_bins = user_bin_map(uniq(user_bin_map,sort(user_bin_map)))
   for i = 0, n_elements(user_input_bins)-1 do begin
      for j = 0, n_elements(binNum_) -1 do begin
         inside_bin_ith = where(user_bin_map eq user_input_bins[i],complement=outise_bin_ith)
         if inside_bin_ith[0] ne -1 then xB=x2d_user[inside_bin_ith]   ;X-coordinates of the points inside the i-th user bin
         if inside_bin_ith[0] ne -1 then yB=y2d_user[inside_bin_ith]   ;Y-coordinates of the points inside the i-th user bin
         dist_B = sqrt( (x2d[ind_good[j]]-xB)^2 + (y2d[ind_good[j]]-yB)^2) ;distances from the j-th datapoint to the points inside the i-th user bin
         if outise_bin_ith[0] ne -1 then xO=x2d_user[outise_bin_ith]   ;X-coordinates of the points outside the i-th user bin
         if outise_bin_ith[0] ne -1 then yO=y2d_user[outise_bin_ith]   ;Y-coordinates of the points outside the i-th user bin
         if outise_bin_ith[0] ne -1 then dist_O = sqrt( (x2d[ind_good[j]]-xO)^2 + (y2d[ind_good[j]]-yO)^2) else  dist_O = 10.^10.;distances from the j-th datapoint to the points outside the i-th user bin
         
         if min(dist_B) lt min(dist_O) then begin  ; the j-th datapoint is located inside the i-th user bin
            binNum_[J] = user_input_bins[I]
            ;kk=1
         endif
         if min(dist_B) eq min(dist_O) then begin  ; In the case it was not already assigned before
            if binNum_[J] eq -1 then  binNum_[J] = user_input_bins[I]
            ;kk=1
         endif
          ;user_input_bins[I] eq 0 then stop
      endfor
      ;if kk eq 1 then print, user_input_bins[I]
      ;if kk eq 0 then print, -10.*user_input_bins[I]
      ;kk=0
   endfor
;stop
   ;computation of bin_sn and coordinate once I know which is the bin each data point belongs to
   valid_user_bins = user_input_bins(where(user_input_bins ge 0))
   n_valid_user_bins =n_elements(valid_user_bins)
   bin_sn = fltarr(n_valid_user_bins)-99
   xNode  = fltarr(n_valid_user_bins)
   yNode = fltarr(n_valid_user_bins)
   binNum__=binNum_ 
   for i = 0, n_valid_user_bins-1 do begin
       sel_bin = valid_user_bins[i]
       if sel_bin eq -1 then continue
       indici = where(binNum_ eq sel_bin[0])
       if indici[0] eq -1 then continue
       bin_sn[i] = total(signal[ind_good[indici]])/sqrt(total(noise[ind_good[indici]]^2))
       ;bin_sn[i] = total(signal[indici])/sqrt(total(noise[indici]^2))
       if n_elements(SN_CALIBRATION) ne 0 then begin   ;SN empirical calibration applied if needed
          tmp = bin_sn[i]
          bin_sn[i] = mdap_calibrate_sn(tmp,n_elements(indici),sn_calibration)
       endif
       xNode[i] = total(x2d[ind_good[indici]]*signal[ind_good[indici]])/total(signal[ind_good[indici]])
       yNode[i] = total(y2d[ind_good[indici]]*signal[ind_good[indici]])/total(signal[ind_good[indici]])
       ;xNode[i] = total(x2d[indici]*signal[indici])/total(signal[indici])
       ;yNode[i] = total(y2d[indici]*signal[indici])/total(signal[indici])
       
   endfor

   xnode=xnode(where(bin_sn gt -99))
   ynode=ynode(where(bin_sn gt -99))
   bin_sn = bin_sn(where(bin_sn gt -99))
   indici_bins = where(binNum_(uniq(binNum_,sort(binNum_))) ne -1)
   nbins = n_elements(indici_bins)


   list_bins_within_binNum_=binNum_(uniq(binNum_,sort(binNum_)))
   if list_bins_within_binNum_[0] eq -1 then list_bins_within_binNum_=list_bins_within_binNum_[1:*]
   ii=0
   for i = 0,nbins-1 do begin
      indici = where(binNum_ eq list_bins_within_binNum_[i])
      if indici[0] ne -1 then begin
         binNum_[indici] = ii
         ii = ii+1
      endif
   endfor
;stop



