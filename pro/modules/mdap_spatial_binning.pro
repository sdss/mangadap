pro mdap_spatial_binning,data,error,signal,noise,min_sn,x2d,y2d,stepx,stepy,$
    binning_map,spectra,errors,xNode,yNode,area_bins,bin_sn,plot=plot,sn_thr=sn_thr,$
    x2d_reconstructed=x2d_reconstructed,y2d_reconstructed=y2d_reconstructed, $
    nelements_within_bin=nelements_within_bin,SN_CALIBRATION=SN_CALIBRATION,$
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
;
;*** INPUT KEYWORDS ***
;
;\plot               If set, some plots on X11 terminal will be shown. Not suggested if the task is launched remotely. 
;
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

version_module = '0.1'
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

ind_good = where(noise gt 0 and finite(noise) eq 1 and signal gt 0 and finite(signal) eq 1 and signal/noise ge sn_thr_,compl=ind_bad)

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


mdap_voronoi_2d_binning,x2d[ind_good],y2d[ind_good],signal[ind_good],noise[ind_good],min_sn,binNum_, $
    xNode, yNode, xBar, yBar, bin_sn, nPixels, scale, plot=plot,SN_CALIBRATION=SN_CALIBRATION, /QUIET;,/no_cvt
print, 'Number of spatial bins: ',mdap_stc(n_elements(xNode),/integer)

;   window,0,retain=2
;   plot, x2d,y2d,psym=4,/iso
;   oplot, x2d[ind_good],y2d[ind_good],psym=4,color=200
;  for i = 0,max(binnum_) do begin
;     indici = where(binNum_ eq i)
;     xyouts,x2d[ind_good[indici]],y2d[ind_good[indici]],mdap_stc(i),charsize=1.2;,color=220
;  endfor

;stop
nbins=n_elements(xnode)
nelements_within_bin=fltarr(n_elements(xnode))

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
      spectra[*,k] = total(data_[ind,*],1)
      errors[*,k] = sqrt(total(error_[ind,*]^2.,1))
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
      spectra[*,k] = total(data[ind_good[ind],*],1)
      errors[*,k] = sqrt(total(error[ind_good[ind],*]^2.,1))
      nelements_within_bin[k] = n_elements(ind)
   endfor
;stop
   reconstructed_R2d = sqrt(x2d_reconstructed^2+y2d_reconstructed^2)
   sz = size(reconstructed_R2d)
   for j = 0, sz[2]-1 do begin
      for i = 0, sz[1]-1 do begin
         dist = sqrt((x2d_reconstructed[i,j] - xnode)^2+(y2d_reconstructed[i,j]- ynode)^2)
         junk = where(dist eq min(dist))
         binning_map[i,j] = junk[0]   ;IF MORE BINS ARE ENCLOSED IN THE RESOLUTION ELEMENT, ONLY THE CLOSEST TO THE CENTER RULES. IN THIS CASE THE MAP IS USELESS, BECAUSE IT FAILS IN ALLOCATING ALL THE BINS IN THE RECONSTRUCTED F.O.V.
         dist_bad = sqrt((x2d_reconstructed[i,j] - x2d[ind_bad])^2+(y2d_reconstructed[i,j]- y2d[ind_bad])^2)
         if min(dist) ge min(dist_bad) then binning_map[i,j] = -1
      endfor
   endfor
   area_bins = fltarr(nbins)
   for i = 0,nbins-1 do  area_bins[i] = nelements_within_bin[i]*!PI*1.25^2 ;area in arcsec2, CONSIDERING A FIBRE OF 2"5 DIAMETER
   
endif 


end_module:

end
