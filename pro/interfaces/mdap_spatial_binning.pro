pro mdap_spatial_binning,data,error,signal,noise,min_sn,x2d,y2d,stepx,stepy,$
    binning_map,spectra,errors,xNode,yNode,area_bins,bin_sn,plot=plot,sn_thr=sn_thr,$
    x2d_reconstructed=x2d_reconstructed,y2d_reconstructed=y2d_reconstructed, $
    nelements_within_bin=nelements_within_bin,SN_CALIBRATION=SN_CALIBRATION,user_bin_map=user_bin_map,weight_for_sn=weight_for_sn,$
    version=version


;*** INPUTS ***
;
;data       Galaxy spectra as produced by mdap_read_datacube.pro.
;           [NxMxT dbl array] if DATACUBE format or [NxT dbl array]if RSS format.
;
;
;error      Errors associated to data, produced by mdap_read_datacube.pro.
;           [NxMxT dbl array] if DATACUBE format or [NxT dbl array]if RSS format.
;
;signal     Mean galaxy signal per \AA, produced by mdap_read_datacube.pro.
;           [NxM dbl array] if DATACUBE format or [N dbl array]if RSS format.
;
;noise      Mean galaxy error per \AA,  produced by mdap_read_datacube.pro. 
;           [NxM dbl array] if DATACUBE format or [N dbl array]if RSS format.
;
;min_sn     [float] Minimum S/N (per \AA) required for the output binned spectra.  
;
;x2d        Array containing the x coordinates in arcseconds (0 is the center of the field of view), produced by  mdap_read_datacube.pro.    
;           [NxM dbl array] if DATACUBE format or [N dbl array]if RSS format.
;
;y2d        Array containing the y coordinates in arcseconds (0 is the center of the field of view), produced by  mdap_read_datacube.pro.    
;           [NxM dbl array] if DATACUBE format or [N dbl array]if RSS format.
;
;stepx      [float] Scale arcsec/pixel along X direction, computed by  mdap_read_datacube.pro. 
;
;stepy      [float] Scale arcsec/pixel along Y direction, computed by  mdap_read_datacube.pro. 
;
;
;*** OPTIONAL INPUTS ***
;
; sn_thr            [float] If specified, spectra with S/N lower than this value will be excluded from the analysis.  
;
; x2d_reconstructed [N'xM' array] Two-dimesional map of X coordinates where the output spatial_binning should be created. Required and used only  
;                    if the input data are in RSS format. 
;
; y2d_reconstructed [N'xM' array] Two-dimesional map of Y coordinates where the output spatial_binning should be created. Required and used only 
;                    if the input data are in RSS format. 
;
;SN_CALIBRATION    If provided, the estimated signal-to-noise (SN_est) is converted into the real signal-to-noise using the empirical
;                  calibration function defined in mdap_calibrate_sn.pro:
;
;                    tmp = SN_EST^SN_CALIBRATION[0]/sqrt(n_elements(n_elements_within_bin)
;                    SN_REAL = poly(SN_EST,SN_CALIBRATION[1:*])
;
;user_bin_map  string   If provided, the spatial map will be created
;                      from the fits fiel specified by this input. The
;                      fits file must contain the CRVAL1, CRVAL2,
;                      CDELT1, CDELT2, NAXIS1, NAXIS2, CRPIX1, and
;                      CRIX2 header keywords (coordinate units in
;                      arcseconds, 0,0 indicates the center of the
;                      field of view.
;
;
;*** INPUT KEYWORDS ***
;
;\plot               If set, some plots on X11 terminal will be shown. Not suggested if the task is launched remotely. 
;
;weight_for_sn      If set, the spectra in the same spatial bin will be
;                   weighted by $S/N^2$ before being added. If The voronoi binning scheme
;                   is adopted, the S/N in the bin is computed via equation (3) of
;                   Cappellari & Copin (2003), and the centers of the
;                   spatial bins are computed by weighting spectra coordinates by $S/N^2$.
;
;*** OUTPUTS ***
;
;binning_map  Two dimensional map showing the binning sheme. Pixels beloning to the i-th bin have value i (i=0, 1, ..., Nbins). 
;             Pixels associated to no spatial bin have value -1.  
;             [NxM dbl array] if inputs are in DATACUBE format or [N'xM' dbl array] if inputs are in RSS format (interpolated 
;             over x2d_reconstructed and y2d_reconstructed).
;
;spectra     [Nbins x T dbl array]    The binned spectra of the spatial Nbins. i-th spectrum is associated to the i-th bin. 
;
;errors      [Nbins x T dbl array]    Errors vectors associate do the binned spectra. 
;
;xNode       [Nbins elements array]   X-Coordinates in arcsec of the luminosity weighted centers of the spatial bins. 
;
;yNode       [Nbins elements array]   Y-Coordinates in arcsec of the luminosity weighted centers of the spatial bins. 
;
;area_bins   [Nbins elements array]   Area (in arcsec$^2$) of each spatial bin.  
;
;bin_sn      [Nbins elements array]   Mean S/N per \AA reached in each spatial bin. 
;

;*** OPTIONAL OUTPUTS ***
;
; nelements_within_bin [Nbins elements array] number of spaxels (in the case of DATACUBE format) or number of fibres (in the case of RSS format) coadded in each spatial bin.
;
; version  [string]            Module version. If requested, the module is not execute and only version flag is returned
;

version_module = '0.2'
if n_elements(version) ne 0 then begin
 version = version_module
 goto, end_module
endif

if keyword_set(plot) then begin
  r = GET_SCREEN_SIZE()
  window, xsize=r[0]*0.4, ysize=r[1]*0.8, retain=2
  loadct, 32
endif

;min_non_Null = min(Noise(where(noise gt 0)))
;noise_=noise
;indici = where(noise eq 0)
;if indici[0] ne -1 then noise_[indici]=min_Non_null
;signal_=signal
;indici = where(signal lt 0)
;if indici[0] ne -1 then signal_[indici] = 0.


;-- Undersand if I am dealing with RSS spectra or a datacube
sz=size(data)  ; sz[0] = 2 --> RSS; sz[0] = 3 --> DATACUBE

sn_thr_=0.
if n_elements(sn_thr) ne 0 then sn_thr_=sn_thr

;ind_good = where(noise gt 0 and finite(noise) eq 1 and signal gt 0 and finite(signal) eq 1 and signal/noise ge sn_thr_,compl=ind_bad)
ind_good = where(noise gt 0 and finite(noise) eq 1  and finite(signal) eq 1 and abs(signal/noise) ge sn_thr_,compl=ind_bad)

;if indici[0] ne -2 then noise[indici] = 10.^10.
if ind_good[0] eq -1 then begin

   if sz[0] eq 3 then begin   ;datacube
      xnode=0
      ynode=0
      bin_sn = 0
      binning_map = signal*0.
      spectra_ = total(data,1)
      spectra = total(temporary(spectra_),1)
      errors = spectra*0.+1.
      area_bins = n_elements(binning_map)*stepx[0]*stepy[0]
      goto, end_module
   endif

   if sz[0] eq 2 then begin  ;RSS spectra
      xnode=0
      ynode=0
      bin_sn = 0
      binning_map = signal*0.
      spectra = total(data,2)
      errors = spectra*0.+1.
      area_bins = spectra*stepx[0]*stepy[0]
      goto, end_module
   endif

endif


if n_elements(user_bin_map) ne 0 then begin    ;USER INPUT SPATIAL BINNING MAP
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

endif else begin  ; VORONOI SPATIAL BINNING SCHEME
   user_map_failure:
  ; if keyword_set(weight_for_sn) then WVT = 1
   mdap_voronoi_2d_binning,x2d[ind_good],y2d[ind_good],signal[ind_good],noise[ind_good],min_sn,binNum_, $
      xNode, yNode, xBar, yBar, bin_sn, nPixels, scale, plot=plot,SN_CALIBRATION=SN_CALIBRATION,/QUIET,weight_for_sn=weight_for_sn;,/no_cvt
   nbins=n_elements(xnode)
  
endelse


print, 'Number of spatial bins: ',mdap_stc(nbins,/integer)

;   window,0,retain=2
;   plot, x2d,y2d,psym=4,/iso
;   oplot, x2d[ind_good],y2d[ind_good],psym=4,color=200
;  for i = 0,max(binnum_) do begin
;     indici = where(binNum_ eq i)
;     xyouts,x2d[ind_good[indici]],y2d[ind_good[indici]],mdap_stc(i),charsize=1.2;,color=220
;  endfor

;stop
nelements_within_bin=fltarr(nbins)

;sz=size(data)
;binNum_ = binNum_-min(binNum_)+1.
;binning_map = reform(binnum_,sz[1],sz[2])
;junk = congrid(reform(binning_map,sz[1],sz[2],1),sz[1],sz[2],sz[3])

if sz[0] eq 3 then begin  ;DATACUBE
   binning_map=x2d*0.-1
   binning_map[ind_good]=binNum_
   x1d_p=findgen(sz[1])
   y1d_p=findgen(sz[2])
   x2d_p=(y1d_p*0.+1.)#x1d_p
   y2d_p=y1d_p#(x1d_p*0.+1.)
   spectra = dblarr(sz[3],nbins)
   errors =  dblarr(sz[3],nbins)
   data_=reform(data,sz[1]*sz[2],sz[3])
   error_=reform(error,sz[1]*sz[2],sz[3])
   binning_map_=reform(binning_map,sz[1]*sz[2])
   area_bins=dblarr(nbins)
;junkx=x2d[uniq(x2d,sort(x2d))]
;junky=y2d[uniq(y2d,sort(y2d))]
;stepx=junkx[1]-junkx[0]
   for k = 0, nbins-1 do begin
      ind = where(binning_map_ eq k)
      if keyword_set(weight_for_sn) then begin
         w = signal[ind_good[ind]]/noise[ind_good[ind]]^2
         tw = total(w)/n_elements(ind)
         ;g = 1./error_[ind,*]^2
         ;GG=errors[*,k]
         e2=error_[ind,*]^2
         for dd = 0, n_elements(ind)-1 do begin
            spectra[*,k] =  spectra[*,k] + data_[ind[dd],*]*w[dd] / tw 
            ;GG = GG + (1./g[dd,*] * (w[dd]/tw)^2)^(-1)
            errors[*,k] =errors[*,k] + e2[dd,*]* (w[dd]/tw)^2
         endfor
          errors[*,k]=sqrt(errors[*,k])
       ;  stop
       ; errors[*,k] = 1/sqrt(GG)*n_elements(ind)   ;If I put the factor *n_elements(ind) I get it right
       
      endif else begin
         spectra[*,k] = total(data_[ind,*],1)
         errors[*,k] = sqrt(total(error_[ind,*]^2.,1))
      endelse
      indici = where( finite(errors[*,k]) ne 1, compl=indici_ok)
    if indici[0] ne -1 then errors[indici,k] = median(errors[indici_ok,k])      
      nelements_within_bin[k] = n_elements(ind)
   endfor
   area_bins=nelements_within_bin*stepx[0]*stepy[0] ;area in arcsec2
endif 
if sz[0] eq 2 then begin  ;RSS
   binning_map=x2d_reconstructed*0.-1
;stop
   spectra = dblarr(sz[2],nbins)
   errors =  dblarr(sz[2],nbins)
   ;data_=reform(data,sz[1]*sz[2],sz[3])
   ;error_=reform(error,sz[1]*sz[2],sz[3])
 
   for k = 0, nbins-1 do begin
      ind = where(binNum_ eq k)
      if ind[0] eq -1 then continue
      if keyword_set(weight_for_sn) then begin
         w = signal[ind_good[ind]]/noise[ind_good[ind]]^2
         tw = total(w)/n_elements(ind)
         e2=error[ind_good[ind],*]^2;  e=1.g^2;;;;g = 1./error[ind_good[ind],*]^2
         ;GG=errors[*,k]
         for dd = 0, n_elements(ind)-1 do begin
            spectra[*,k] =  spectra[*,k] + data[ind_good[ind[dd]],*]*w[dd] / tw 
            errors[*,k] = errors[*,k] + e2[dd,*]*(w[dd] / tw )^2
         endfor
         ;errors[*,k] = 1/sqrt(GG)*n_elements(ind)   ;If I put the factor *n_elements(ind) I get it right
          errors[*,k] = sqrt(errors[*,k])
     endif else begin
         spectra[*,k] = total(data[ind_good[ind],*],1)
         errors[*,k] = sqrt(total(error[ind_good[ind],*]^2.,1))
      endelse
    ; stop
     indici = where( finite(errors[*,k]) ne 1, compl=indici_ok)
     if indici[0] ne -1 then errors[indici,k] = median(errors[indici_ok,k])
     nelements_within_bin[k] = n_elements(ind)
  endfor

;stop
   reconstructed_R2d = sqrt(x2d_reconstructed^2+y2d_reconstructed^2)
   sz = size(reconstructed_R2d)

 ;plot,x2d_reconstructed,y2d_reconstructed,/nodata,/iso
 ;oplot,xnode,ynode,psym=7,color=200
   for j = 0, sz[2]-1 do begin
      for i = 0, sz[1]-1 do begin
         dist = sqrt((x2d_reconstructed[i,j] - xnode)^2+(y2d_reconstructed[i,j]- ynode)^2)
         ;dist_rss = sqrt((x2d_reconstructed[i,j] - x2d[])^2+(y2d_reconstructed[i,j]- y2d)^2)
         dist_rss = sqrt((x2d_reconstructed[i,j] - x2d[ind_good])^2+(y2d_reconstructed[i,j]- y2d[ind_good])^2)
         junk = where(dist eq min(dist))
         binning_map[i,j] = junk[0]   ;IF MORE BINS ARE ENCLOSED IN THE RESOLUTION ELEMENT, ONLY THE CLOSEST TO THE CENTER RULES. IN THIS CASE THE MAP IS USELESS, BECAUSE IT FAILS IN ALLOCATING ALL THE BINS IN THE RECONSTRUCTED F.O.V.
         if ind_bad[0] ne -1 then dist_bad = sqrt((x2d_reconstructed[i,j] - x2d[ind_bad])^2+(y2d_reconstructed[i,j]- y2d[ind_bad])^2) else dist_bad=10.^10.
         ;if min(dist) ge min(dist_bad) then binning_map[i,j] = -1
    
         if min(dist_rss) ge min(dist_bad) or min(dist_rss) ge 1.5 then binning_map[i,j] = -1   ;If I am closer to a bad pixel than a good pixel, or if I am more distant than 3" to a fibre location
         ;if binning_map[i,j] ne -1 then plots,x2d_reconstructed[i,j] ,y2d_reconstructed[i,j] ,psym=4
         ;if binning_map[i,j] eq -1 then plots,x2d_reconstructed[i,j] ,y2d_reconstructed[i,j] ,psym=5
      endfor
   endfor
   area_bins = fltarr(nbins)
   for i = 0,nbins-1 do  area_bins[i] = nelements_within_bin[i]*!PI*1.25^2 ;area in arcsec2, CONSIDERING A FIBRE OF 2"5 DIAMETER
   
endif 


end_module:

end
