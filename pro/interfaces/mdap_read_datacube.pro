pro mdap_read_datacube,datacube_name,data,error,wavelength,$
      x2d,y2d,signal,noise,cdelt1,cdelt2,header2d,lrange=lrange,keep_original_step=keep_original_step,$
     x2d_reconstructed=x2d_reconstructed,y2d_reconstructed=y2d_reconstructed,$
     signal2d_reconstructed=signal2d_reconstructed,noise2d_reconstructed=noise2d_reconstructed,version=version,$
      number_of_fibres=number_of_fibres,use_total=use_total
      
;*** INPUTS ***
;
;datacube_name  [string] Name of the fits file containing the input in datacube or RSS formats. 
;                         For inputs in datacube format, the file must be a fits file with the following
;                         extensions:
;
;                            1. flux in 1e-17 ergs/s/cm$^{-2}$/\AA. This extension must be a 3D
;                               array, with the wavelength direction along the 3rd axis.
; 
;                            2. Inverse variance associated to the first extension
;
;                            3. Wavelength solution (1D dbl array). Constant linear step
;                              (preferred)? Constant Log-step? Constant ln-step?
;
;                            4. RSS format (row-stacked spectra) (NOT USED)
;                  
;                            5. coverage map; (NOT USED)
;
;                         For inputs in RSS format, the file must be a fits file with the following
;                         extensions:
;
;                            1. flux in 1e-17 ergs/s/cm$^{-2}$/\AA. This extension must be a 3D
;                               array, with the wavelength direction along the 3rd axis.
; 
;                            2. Inverse variance associated to the first extension
;
;                            3. Wavelength solution (1D dbl array). Constant linear step
;                              (preferred)? Constant Log-step? Constant ln-step?
;                            
;                            4. X position table.  Since position is a function of wavelength this is an 
;                               array of the same size as the flux array.  Coordinates should be in arcseconds 
;                               relative to the center of the data cube (i.e., roughly the center of the galaxy).
;                            5. Y position table.

;
;*** OPTIONAL INPUTS ***
;
;lrange   [2 elem vector]  Wavelentgh range (in angstrom) where to extract the information for Signal 
;                         and Noise. Default: use the entire spectra range.
;
;*** OPTIONAL KEYWORDS ***
;
;/keep_original_step      If set, the wavelength output vector will be the same as the one
;                         define from the input fits file. The default
;                         is to re-construct (and  therefore re-inrpolate the galaxy and error
;                         spectra) the output wavelength vector with constant ang/pixel step
;                         using the minimum ang/pixel step that is stored in the wavelength
;                         solution. For MANGA, it is suggested to set this keyword.
;
; /use_total              If set, the signal is the total of the counts in the selected wavelength range, the noise is the sum in
;                         quadrature of the noises in the selected range. Useful for emission lines.
;
;*** OUTPUTS ***
;
;data       [NxMxT dbl array]  Datacube of galaxy spectra. Spectra are resampled over the vector wavelength.
;        or [NxT dbl array]    If input is in RSS format, data is a [NxT dbl array], where N is the number of RSS spectra.
;
;error      [NxMxT dbl array]  Errors associated to data.
;        or [NxT dbl array]    If input is in RSS format, error is a [NxT dbl array], where N is the number of RSS spectra.
;
;wavelength [T dbl array]      Wavalenght value where data and errors are computed. The vector is constructed 
;                              with constant linear step in ang/pixel (unless the \keep_original_step keyword is selected). 
;                              If input spectra have a logarithmic sampling, the minimum available step is used 
;                              (e.g. log_lam[1]-log_lam[0]).
;
;x2d        [NxM dbl array]    Array containing the x coordinates in arcseconds (0 is the center of the field of view).
;        or [N dbl array]      If inputs are in RSS format, it is a N elements float vector
;
;y2d        [NxM dbl array]    Array containing the y coordinates in arcseconds (0 is the center of the field of view).
;        or [N dbl array]      If inputs are in RSS format, it is a N elements float vector
;
;signal     [NxM dbl array]    Mean galaxy signal per \AA, obtained considering all the wavelength range 
;        or [N dbl array]     (or only the range specified by lrange). Calculation is done on original spectra, not those 
;                              resampled over the vector wavelenght.
;                              The signal for each i,j-th spectrum is calculated as median(signal), unless the keyword /use\_total is set.
;                              If inputs are in RSS format, it is a N elements float vector
;
;noise      [NxM dbl array]    Mean galaxy error per \AA, obtained considering all the wavelength range (or only the range 
;        or [N dbl array]      specified by lrange). Calculation is done on original spectra, not those resampled over the 
;                              vector wavelenght.
;                              The noise for each i,j-th spectrum is calculated as median(noise), unless the keyword /use\_total is set.
;                              If inputs are in RSS format, it is a N elements float vector
;
;cdelt1     [double]           spatial sampling along x direction (arcsec/pixel). 
;                              If inputs are in RSS format, it is set to 0.5 
;
;cdelt2     [double]           spatial sampling along y direction (arcsec/pixel).    
;                              If inputs are in RSS format, it is set to 0.5 
;
;header2d   [str array]        The header for two-dimensional fields
;                              If inputs are in RSS format, TO BE DONE
;
;
;*** OPTIONAL OUTPUTS ***
; x2d_reconstructed     [NxM array]   If input is in DATACUBE format, x2d_reconstructed = x2d.
;                   or  [N'xM' array] If inputs are in RSS format, the x2d coordinates are resampled over a 2D map with fixed 0"5/pixel sampling
;                                     and define the  x2d_reconstructed map.
;
; y2d_reconstructed     [NxM array]   If input is in DATACUBE format, y2d_reconstructed = y2d.
;                   or  [N'xM' array] If inputs are in RSS format, the y2d coordinates are resampled over a 2D map with fixed 0"5/pixel sampling
;                                     and define the  y2d_reconstructed map.
;
; signal2d_reconstructed    [NxM array]   If input is in DATACUBE format, signal2d_reconstructed = signal.
;                        or [N'xM' array] If inputs are in RSS format, the signal is resampled over the 2D map defined by
;                                         x2d_reconstructed  and y2d_reconstructed. 
;
; noise2d_reconstructed     [NxM array]   If input is in DATACUBE format, noise2d_reconstructed = noise.
;                        or [N'xM' array] If inputs are in RSS format, the noise is resampled over the 2D map defined by
;                                         x2d_reconstructed  and y2d_reconstructed. 
;
; version  [string]      Module version. If requested, the module is not execute and only version flag is returned
;

version_module = '0.2'
if n_elements(version) ne 0 then begin
 version = version_module
 goto, end_module
endif

data=double(readfits(datacube_name,header,exten_no=0)) ; data)
inv=double(readfits(datacube_name,exten_no=1))  
error=sqrt(1./temporary(inv))                      ;
wavelength=double(readfits(datacube_name,exten_no=2)) ;lambdas
lam_sel=indgen(n_elements(wavelength))
if n_elements(lrange) ne 0 then lam_sel = where(wavelength ge lrange[0] and wavelength le lrange[1])
sz=size(data)
if sz[0] eq 3 then RSS = 0   ; INPUT IS IN DATACUBE FORMAT
if sz[0] eq 2 then RSS = 1   ; INPUT IS IN RSS FORMAT

if RSS eq 1 then begin
   RESOLUTION_ELEMENT = 0.5
   data_ = transpose(temporary(data))
   error_ = transpose(temporary(error))
   x2d_=readfits(datacube_name,exten_no=3)
   y2d_=readfits(datacube_name,exten_no=4)

   x2d=median(temporary(x2d_[lam_sel,*]),/even,dimension=1)
   y2d=median(temporary(y2d_[lam_sel,*]),/even,dimension=1)

   ;--
   ;   KEEPING FIBRES, DITHERS, AND EXPOSURES SEPARATED
   sz=size(data_)
   data  = temporary(data_)
   error = temporary(error_)
   ;  
  



;   ;--
;   ;   ;   ; PER FIBRE/DITHER VERSION
;   sz=size(data_)
;   NUMBER_OF_DITHERS = 3.
;   exposures = sz[1]/number_of_fibres/NUMBER_OF_DITHERS
;   data = fltarr(number_of_fibres*NUMBER_OF_DITHERS,sz[2])
;   error2 = fltarr(number_of_fibres*NUMBER_OF_DITHERS,sz[2])
;   x2d=fltarr(number_of_fibres*NUMBER_OF_DITHERS)
;   y2d=fltarr(number_of_fibres*NUMBER_OF_DITHERS)
;   for i = 0, number_of_fibres*NUMBER_OF_DITHERS-1 do begin
;      tmpx=0
;      tmpy=0
;      
;       for j = 0, exposures-1 do begin
;           data[i,*]=data[i,*]+data_[i+j*number_of_fibres*NUMBER_OF_DITHERS,*]
;        
;           ;print, i, i+j*number_of_fibres*NUMBER_OF_DITHERS
;           error2[i,*]=error2[i,*]+error_[i+j*number_of_fibres*NUMBER_OF_DITHERS,*]^2
;           ;x2d[i] = x2d[i] + x2d_[i+j*number_of_fibres]
;           ;y2d[i] = y2d[i] + y2d_[i+j*number_of_fibres]
;           tmpx=[tmpx,x2d_[i+j*number_of_fibres*NUMBER_OF_DITHERS]]
;           tmpy=[tmpy,y2d_[i+j*number_of_fibres*NUMBER_OF_DITHERS]]
;      endfor     
;      ; plot, data[i,*]
;       x2d[i] = median(tmpx[1:*],/even);x2d[i]/exposures
;       y2d[i] = median(tmpy[1:*],/even);y2d[i]/exposures
;      ; stop
;   endfor
;   error=sqrt(temporary(error2))
;;window,3,retain=2
;;plot, x2d_, y2d_,psym=4,/iso;,xrange=[-3,3],yrange=[-3,3]
;;oplot, x2d_[0:126], y2d_[0:126],psym=6,color=200
;;oplot, x2d_[127:127*2-1], y2d_[127:127*2-1],psym=6,color=200
;;oplot, x2d_[127*2:127*3-1], y2d_[127*2:127*3-1],psym=6,color=200
;;
;;if exposures gt 1 then begin
;;   oplot, x2d_[127*3:127*4-1], y2d_[127*3:127*4-1],psym=6,color=200
;;   oplot, x2d_[127*4:127*5-1], y2d_[127*4:127*5-1],psym=6,color=200
;;   oplot, x2d_[127*5:127*6-1], y2d_[127*5:127*6-1],psym=6,color=200
; ;  if exposures gt 2 then begin
; ;     oplot, x2d_[127*6:127*7-1], y2d_[127*6:127*7-1],psym=6,color=200
;;      oplot, x2d_[127*7:127*8-1], y2d_[127*7:127*8-1],psym=6,color=200
;;      oplot, x2d_[127*8:127*9-1], y2d_[127*8:127*9-1],psym=6,color=200
;;   en;dif
;;endif;
;
;;for i = 0, n_elements(x2d)-1 do xyouts,x2d[i],y2d[i],mdap_stc(i+1),color=220,charsize=1.2;
;
;   if number_of_fibres eq 127 then begin  ;fix coordinate bug in 127 fibre bundle
;   ;   plot, x2d, y2d,psym=4,/iso
;      indices_to_fix=0
;      corner_a=0
;      corner_b=0
;      corner_c=0
;      corner_d=0
;      corner_e=0
;      corner_f=0
;      for kdit = 0, NUMBER_OF_DITHERS-1 do begin
;         indices_to_fix = [indices_to_fix,68+number_of_fibres*kdit-1]
;         corner_a =  [corner_a,67+number_of_fibres*kdit-1]
;         corner_b =  [corner_b,48+number_of_fibres*kdit-1]
;         corner_c =  [corner_c,47+number_of_fibres*kdit-1]
;         corner_d =  [corner_d,69+number_of_fibres*kdit-1]
;         corner_e =  [corner_e,72+number_of_fibres*kdit-1]
;         corner_f =  [corner_f,73+number_of_fibres*kdit-1]
;      endfor
;      indices_to_fix =indices_to_fix[1:*]
;      corner_a  = corner_a[1:*]
;      corner_b  = corner_b[1:*]
;      corner_c  = corner_c[1:*]
;      corner_d  = corner_d[1:*]
;      corner_e  = corner_e[1:*]
;      corner_f  = corner_f[1:*]
;
;      ;for kdit = 0, number_of_fibres*NUMBER_OF_DITHERS -1 do begin
;      for counter = 0, n_elements(indices_to_fix)-1 do begin
;          for jexp = 0, exposures-1 do begin
;            a = corner_a[counter]+jexp*number_of_fibres*NUMBER_OF_DITHERS
;            b = corner_b[counter]+jexp*number_of_fibres*NUMBER_OF_DITHERS
;            c = corner_c[counter]+jexp*number_of_fibres*NUMBER_OF_DITHERS
;            d = corner_d[counter]+jexp*number_of_fibres*NUMBER_OF_DITHERS
;            e = corner_e[counter]+jexp*number_of_fibres*NUMBER_OF_DITHERS
;            f = corner_f[counter]+jexp*number_of_fibres*NUMBER_OF_DITHERS
;            tmpx = [x2d_[a],x2d_[b],x2d_[c],x2d_[d],x2d_[e],x2d_[f]]
;            tmpy = [y2d_[a],y2d_[b],y2d_[c],y2d_[d],y2d_[e],y2d_[f]]
;         ;   print, a,b,c,d,e,f
;         ;     plot, x2d, y2d,psym=4,/iso
;         ;     plots,[x2d_[a],x2d_[b],x2d_[c],x2d_[d],x2d_[e],x2d_[f]],[y2d_[a],y2d_[b],y2d_[c],y2d_[d],y2d_[e],y2d_[f]],psym=4+jexp,color=200
;          ;top
;          endfor
;          xf = median(tmpx,/even)
;          yf = median(tmpy,/even)
;          ;top
;         if abs(x2d[indices_to_fix[counter]] - xf) gt 4 then x2d[indices_to_fix[counter]] = xf
;         if abs(y2d[indices_to_fix[counter]] - yf) gt 4 then y2d[indices_to_fix[counter]] = yf
;      endfor
;
;   endif



;   ;
;   ; PER FIBRE VERSION
;    sz=size(data_)
; 
;    exposures = sz[1]/number_of_fibres
;    data = fltarr(number_of_fibres,sz[2])
;    error2 = fltarr(number_of_fibres,sz[2])
;    x2d=fltarr(number_of_fibres)
;    y2d=fltarr(number_of_fibres)

;   for i = 0, number_of_fibres-1 do begin
;      tmpx=0
;      tmpy=0
;       for j = 0, exposures-1 do begin
;           data[i,*]=data[i,*]+data_[i+j*number_of_fibres,*]
;           error2[i,*]=error2[i,*]+error_[i+j*number_of_fibres,*]^2
;           ;x2d[i] = x2d[i] + x2d_[i+j*number_of_fibres]
;           ;y2d[i] = y2d[i] + y2d_[i+j*number_of_fibres]
;           tmpx=[tmpx,x2d_[i+j*number_of_fibres]]
;           tmpy=[tmpy,y2d_[i+j*number_of_fibres]]
;      endfor     
;      ; plot, data[i,*]
;       x2d[i] = median(tmpx[1:*],/even);x2d[i]/exposures
;       y2d[i] = median(tmpy[1:*],/even);y2d[i]/exposures
;      ; stop
;   endfor
;   error=sqrt(temporary(error2))
;   if number_of_fibres eq 127 then begin  ;fix coordinate bug in 127 fibre bundle
;       xf=median(x2d[[67,48,47,69,72,73]-1],/even)
;       yf=median(y2d[[67,48,47,69,72,73]-1],/even)
;       if abs(x2d[67] - xf) gt 1.5 then x2d[67] = xf
;       if abs(y2d[67] - yf) gt 1.5 then y2d[67] = yf
;   endif
     ;plot, x2d,y2d,psym=4,/iso
     ;for jj = 0, n_elements(x2d)-1 do xyouts,x2d[jj],y2d[jj],mdap_stc(jj+1),charsize=1.1,color=210
     ;   stop
   ;--


endif
sz=size(data)

; removing spectra if bad defined for more than 20% of wavelength extension
if RSS eq 0 then begin      ;datacube
for j = 0, sz[2]-1 do begin
   for i = 0 , sz[1] -1 do begin
      indici = where(error[i,j,*] le 0 or finite(error[i,j,*]) ne 1 or finite(data[i,j,*]) ne 1)
      ;if n_elements(indici) ge 0.2 * n_elements(data[i,j,*]) then begin
      if n_elements(indici) ge 0.2 * n_elements(data[i,j,*]) or min(data[i,j,*]) eq max(data[i,j,*]) then begin
         data[i,j,*] = data[i,j,*]*0./0.
         error[i,j,*] = error[i,j,*]*0.+1.
      endif
   endfor
endfor
endif
if RSS eq 1 then begin      ;datacube
   for i = 0 , sz[1] -1 do begin
      indici = where(error[i,*] le 0 or finite(error[i,*]) ne 1 or finite(data[i,*]) ne 1)
   ;   if n_elements(indici) ge 0.2 * n_elements(data[i,*]) then begin
      if n_elements(indici) ge 0.2 * n_elements(data[i,*]) or min(data[i,*]) eq max(data[i,*]) then begin
         data[i,*] = data[i,*]*0./0.
         error[i,*] = error[i,*]*0.+1.
      endif
   endfor
endif

; defyning S/N from data, before any re-sampling
if RSS eq 0 then begin      ;datacube
   signal = fltarr(sz[1],sz[2])
   noise = fltarr(sz[1],sz[2])
   for j = 0, sz[2]-1 do begin
      for i = 0, sz[1]-1 do begin
         if ~keyword_set(use_total) then begin
            mdap_calculate_spectrum_sn,data[i,j,lam_Sel],error[i,j,lam_Sel],wavelength[lam_Sel],junk,signal=sss,noise=nnn
            signal[i,j]=sss[0]
            noise[i,j] = nnn[0]
         endif else begin
            signal[i,j]=total(data[i,j,lam_Sel])
            noise[i,j] = sqrt(total(error[i,j,lam_Sel]^2))
         endelse
      endfor
   endfor
endif else begin        ; rss
  ; stop
   
   if ~keyword_set(use_total) then begin
      signal = median(data[*,lam_Sel],dimension=2,/even)
      noise  = median(error[*,lam_Sel],dimension=2,/even)
   endif else begin
      nrss = n_elements(data[*,0])
      signal=fltarr(nrss)
      noise=fltarr(rss)
      for j = 0, nrss-1 do begin
          signal[j] = total(data[j,lam_Sel])
          noise[j] = sqrt(total(error[j,lam_Sel]^2))
      endfor
   endelse

   indici = where(finite(noise) eq 0 or finite(signal) eq 0  or noise le 0,compl=pos)
   if indici[0] ne -1 then noise(indici) = max(noise[pos])
   if indici[0] ne -1 then signal(indici) = 0./0.
endelse

;-- resampling over a new wavelenght range with a constant ang/pxl step, if required.
if ~keyword_set(keep_original_step) then begin
 
   min_step_lin = min( wavelength[1:n_elements(wavelength)-2] - wavelength[0:n_elements(wavelength)-1])
   nw=round((max(wavelength)-min(wavelength))/min_step_lin)
   wavelength_linear_constant_step = dindgen(nw)*min_step_lin+min(wavelength)

   if RSS eq 0 then begin  ;datacube
      data_interpolated = dblarr(sz[1],sz[2],n_elements(wavelength_linear_constant_step))
      error_interpolated = dblarr(sz[1],sz[2],n_elements(wavelength_linear_constant_step))
      for i = 0, sz[1]-1 do begin
         for j = 0, sz[2]-1 do begin
            data_interpolated[i,j,*] = interpol(data[i,j,*],wavelength,wavelength_linear_constant_step)
            error_interpolated[i,j,*] = interpol(error[i,j,*],wavelength,wavelength_linear_constant_step)
         endfor
      endfor
      wavelength = temporary(wavelength_linear_constant_step)
      data = temporary(data_interpolated)
      error = temporary(error_interpolated)

   endif else begin   ;RSS

      data_interpolated = dblarr(n_elements(wavelength_linear_constant_step),sz[2])
      error_interpolated = dblarr(n_elements(wavelength_linear_constant_step),sz[2])
      for i = 0, sz[1]-1 do begin
          data_interpolated[i,*] = interpol(data[i,*],wavelength,wavelength_linear_constant_step)
          error_interpolated[i,*] = interpol(error[i,*],wavelength,wavelength_linear_constant_step)
      endfor

      wavelength = temporary(wavelength_linear_constant_step)
      data = temporary(data_interpolated)
      error = temporary(error_interpolated)
   endelse

endif
;--


indx=where(error le 0 or finite(error) ne 1,compl=pos)
;if indx[0] ne -1 then error[indx] = max(error[pos])
if indx[0] ne -1 then error[indx] = median(error[pos])


sz=size(data)

; spatial coordinates
if RSS eq 0 then begin; IF I HAVE DATACUBES
   cdelt1=sxpar(header,'CDELT1')
   cdelt2=sxpar(header,'CDELT2')
   x=(findgen(sz[1])-(sz[1]-1)/2.)*cdelt1
   y=(findgen(sz[2])-(sz[2]-1)/2.)*cdelt2
   x2d=x#(y*0+1.)
   y2d=(x*0.+1.)#y
   x2d_reconstructed = x2d
   y2d_reconstructed = y2d
   signal2d_reconstructed = signal
   noise2d_reconstructed = noise

   header2d=["SIMPLE  =                    T / Written by IDL:",$
             "BITPIX  =                  -32 / Number of bits per data pixel",$
             "NAXIS   =                    3 / Number of data axes",$
             "EXTEND  =                    T / FITS data may contain extensions",$
             "CRPIX1  =              "+mdap_stc((sz[1]+1)/2.)+" /",$
             "CTYPE1  =                'arcsec' /",$
             "CRVAL1  =                  0.0    /",$
             "CDELT1  =          "+string(cdelt1)+" /",$
             "CRPIX2  =              "+mdap_stc((sz[2]+1)/2.)+" /",$
             "CTYPE2  =                'arcsec' /",$
             "CRVAL2  =                  0.0    /",$
             "CDELT2  =          "+string(cdelt2)+" /",$
             "OBJECT  =   "+datacube_name+"  /"]
endif else begin  ; IF I HAVE RSS SPECTRA, I NEED TO INTERPOLATE SOME INFO ON A 2D ARRAY
   cdelt1 = RESOLUTION_ELEMENT
   cdelt2 = RESOLUTION_ELEMENT
   max_d = round(max([abs(x2d),abs(y2d)]))+2.
   nbinsx = max_d*2/cdelt1
   nbinsx = fix(nbinsx/2)*2+1
   nbinsy = max_d*2/cdelt2
   nbinsy = fix(nbinsy/2)*2+1

   x=cdelt1*(findgen(nbinsx)-(nbinsx-1)/2)
   y=cdelt2*(findgen(nbinsy)-(nbinsy-1)/2)

   x2d_reconstructed=x#(y*0.+1.)
   y2d_reconstructed=(x*0.+1.)#y
   
  
  ; mdap_interpolate_2dmaps,signal,x2d,y2d,x2d_reconstructed,y2d_reconstructed,x2d_reconstructed,y2d_reconstructed,signal2d_reconstructed
  ; signal2d_reconstructed = reform(signal2d_reconstructed,nbinsx,nbinsy)
  ; mdap_interpolate_2dmaps,noise,x2d,y2d,x2d_reconstructed,y2d_reconstructed,x2d_reconstructed,y2d_reconstructed,noise2d_reconstructed
  ; noise2d_reconstructed = reform(noise2d_reconstructed,nbinsx,nbinsy)
  signal2d_reconstructed = x2d_reconstructed/0.
  noise2d_reconstructed =  x2d_reconstructed/0.
  for j = 0, nbinsy-1 do begin
     for i = 0, nbinsx-1 do begin
        dist = sqrt((x2d_reconstructed[i,j]-x2d)^2+(y2d_reconstructed[i,j]-y2d)^2)
        min_dist = min(dist)
        indici = where(dist le min_dist)
        if min_dist le 1.5 then signal2d_reconstructed[i,j] = mean(signal[indici])
        if min_dist le 1.5 then noise2d_reconstructed[i,j] = mean(noise[indici])
     endfor
  endfor
 
   header2d=["SIMPLE  =                    T / Written by IDL:",$
             "BITPIX  =                  -32 / Number of bits per data pixel",$
             "NAXIS   =                    3 / Number of data axes",$
             "EXTEND  =                    T / FITS data may contain extensions",$
             "CRPIX1  =              "+mdap_stc((nbinsx+1)/2.)+" /",$
             "CTYPE1  =                'arcsec' /",$
             "CRVAL1  =                  0.0    /",$
             "CDELT1  =          "+string(cdelt1)+" /",$
             "CRPIX2  =              "+mdap_stc((nbinsy+1)/2.)+" /",$
             "CTYPE2  =                'arcsec' /",$
             "CRVAL2  =                  0.0    /",$
             "CDELT2  =          "+string(cdelt2)+" /",$
             "OBJECT  =   "+datacube_name+"  /"]

endelse




;-- PATCH added from version 0.7 onwards, to center the coordinates
;   onto the galaxy center. Header informations are updated accordingly.
; get image center from signal2d_reconstructed
sz_=size(signal2d_reconstructed)
gcntrd,signal2d_reconstructed,sz_[1]/2.,sz_[2]/2.,xcen,ycen,3.
xcen_pix=xcen
ycen_pix=ycen
x0=bilinear(x2d_reconstructed,xcen_pix,ycen_pix) 
xcen_pix=xcen
ycen_pix=ycen
y0=bilinear(y2d_reconstructed,xcen_pix,ycen_pix)
if abs(x0) ge 4 then x0 = 0.
if abs(y0) ge 4 then y0 = 0.

x2d=x2d-x0
y2d=y2d-y0
x2d_reconstructed=x2d_reconstructed-x0
y2d_reconstructed=y2d_reconstructed-y0
header2d[6] = "CRVAL1  =                 "+mdap_stc(-x0)+"    /"
header2d[10] = "CRVAL2  =                "+mdap_stc(-y0)+"    /"

;-------------------------------------------------





end_module:
end
