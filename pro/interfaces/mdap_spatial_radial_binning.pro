; pro test
; restore,'ma008_cgcg122-022_mdap_session.idl'
; 
; bin_sn_ems_real = bin_sn_ems
; indxx = where(wavelength_output_rest_frame_log ge 4500 and wavelength_output_rest_frame_log le 7000)
; for i = 0, n_elements(bin_sn_ems) -1 do begin
;    mdap_calculate_spectrum_sn,best_fit_model_ems[i,indxx],residuals_ems[i,indxx],$
;     wavelength_output_rest_frame_log[indxx],sn_per_angstorm,/rms
;     bin_sn_ems_real[i] = sn_per_angstorm[0]
; endfor
; 
; mdap_get_error_from_residual,residuals_ems,galaxy_minus_ems_fit_model_ems,input_errors
; 
; 
; mdap_spatial_radial_binning,bin_sn_ems_real,x2d,y2d,spatial_binning_ems,xbin_ems,ybin_ems,ell,pa,$
;       galaxy_minus_ems_fit_model_ems,input_errors,wavelength_output_rest_frame_log,$
;       spatial_binning_scheme,r_bin,r_bin_lo,r_bin_up,r2d_bin,r2d_bin_lo,r2d_bin_up,binned_spectra,binned_errors,$
;       output_lrange=[3600,8500],output_wav=output_wav,n_elements_bin=n_elements_bin
; 
; writefits,'test.fits', spatial_binning_scheme
; stop
; end

pro mdap_spatial_radial_binning,signal_to_noise,x2d,y2d,reference_2dmap,xbin,ybin,ell,pa,$
      input_spectra,input_errors,wav,$
      spatial_binning_scheme,r_bin,r_bin_lo,r_bin_up,r2d_bin,r2d_bin_lo,r2d_bin_up,binned_spectra,binned_errors,$
      output_lrange=output_lrange,output_wav=output_wav,n_elements_bin=n_elements_bin,$
      low_radial_bins_user_inputs=low_radial_bins_user_inputs,upper_radial_bins_user_inputs=upper_radial_bins_user_inputs,$
      Reff_=Reff_,PSFsize_=PSFsize_,add_default_bins=add_default_bins,$
      version=version


;INPUTS
;
; signal_to_noise  S/N measured of the input spectra (measures on best
;                  fit model and residuals) 
; x2d              output of mdap_read_datacube
; y2d              output of mdap_read_datacube
; reference_2dmap  output of mdap_spatial_binning (ems run) (pixels
;                  with no signal must have value -1)
; xbin             output of mdap_spatial_binning (ems run)
; ybin             output of mdap_spatial_binning (ems run) 
; ell              Galaxy mean ellipticity (el = 1-b/a)
; PA               Galaxy mean position angle (0=north ; 90=east)
; input_spectra    output of mdap_spectral_fitting (ems run): galaxy_minus_ems_fit_model_ems.
;                  galaxy spectra with no emission lines, log rebinned, and rest frame.
;                  Dimension [NxM]  N=number of bins, M=number of wavelength points per spectrum. 
;                  N is the same number of xbin and ybin
; input_errors     Errors associated to the input_spectra 
; wav              Wavelenght (rest frame) where the input_spectra and
;                  input_errors are sampled.
;
; OPTIONAL INPUT 
; output_lrange    range to trim the output spectra. Default: no trimming
;
; OUTPUTS
;spatial_binning_scheme.  same dimensions as signal. Pixels belonging
;                         to the same radial bin have the same value. Pixels with no signal (no
;                         bin assigned) have value -1, pixels of the same bin have values that
;                         range from 0 to Nbins -1
;
; r_bin 
; r_bin_lo
; r_bin_up
; r2d_bin 
; r2d_bin_lo
; r2d_bin_up
; binned_spectra  TxM array T is the number of radial bins
; binned_errors  TxM array T is the number of radial bins
;
; radial_bins_user_inputs=radial_bins_user_inputs,Reff=Reff,PSFsize=PSFsize, TO BE DESCRIBED...
;
; OPTIONAL OUTPUTS
; output_wav
;
version_module = '0.7' ;16/04/2014

if n_elements(version) ne 0 then begin
 version = version_module
 goto, end_module
endif


;- Check the number of input spectra, and decide the number of spectra
;  per bin accordingly.

sz = size(input_spectra)
if sz[1] eq 1 then begin   ; only 1 input spectrum.... I do not bin.
   spatial_binning_scheme = 0./0.
   spatial_binning_scheme[where(finite(signal) eq 1)] = 0
   r_bin = [0]
   r_bin_lo = [0]
   r_bin_up = max(sqrt(x2d^2+y2d^2))
   binned_spectra  = reform(input_spectra)
   binned_errors  = reform(input_errors)
   binned_spectra  = double(binned_spectra[output_lrange])
   binned_errors  = double(binned_errors[output_lrange])
   n_elements_per_r_bin  = n_elements(where(finite(signal) eq 1))
   r2d_bin = spatial_binning_scheme*0.+r_bin
   r2d_bin_lo = spatial_binning_scheme*0.+r_bin_lo
   r2d_bin_up = spatial_binning_scheme*0.+r_bin_up
   nbins =1
   goto, end_module
endif


spatial_binning_scheme = reference_2dmap*0
r2d_bin  = reference_2dmap*0/0.
r2d_bin_lo = reference_2dmap*0/0.
r2d_bin_up = reference_2dmap*0/0.

indici_null = where(reference_2dmap eq -1,complement=ind_not_null)
if indici_null[0] ne -1 then begin
   spatial_binning_scheme[indici_null] = -1
endif

;-- rotate coordianates to align galaxy photometric major axis to
;   vertical direction
x2d_r =  x2d * cos(pa*!pi/180.) + y2d*sin(pa*!pi/180.)
y2d_r = -x2d * sin(pa*!pi/180.) + y2d*cos(pa*!pi/180.)
xbin_r =  xbin * cos(pa*!pi/180.) + ybin*sin(pa*!pi/180.)  ;clockwise rotation of PA
ybin_r = -xbin * sin(pa*!pi/180.) + ybin*cos(pa*!pi/180.)
;--

;-- computation of semimajor axis (aligned vertically, with ellipticity
;   as input value
abin = sqrt(xbin_r^2/(1.-ell)^2 + ybin_r^2)
a2d=sqrt(x2d_r^2/(1.-ell)^2 + y2d_r^2)
a2d_uniq = a2d[uniq(a2d,sort(a2d))]
;--

;-- definition of radial bins
IF n_elements(low_radial_bins_user_inputs) eq 0 or n_elements(upper_radial_bins_user_inputs) eq 0 or n_elements(low_radial_bins_user_inputs) ne n_elements(upper_radial_bins_user_inputs)  then begin
   nbins = 6
   if sz[1] lt 50 then nbins =5
   if sz[1] lt 15 then nbins =3
   if sz[1] lt 10 then nbins =2  
   
   junk = mdap_range(1,max(abin)*1.1,nbins+1,/log)
   junk[0] = 0
   r_bin_up=junk[1:*]
   r_bin_lo=junk[0:nbins-1]
endif else begin   ; USER INPUT

   if n_elements(Reff_) eq 0 then Reff = max(abin)/2. else Reff=Reff_ ; If Reff is not defined, I set it to be halph of the maximum semimajor axis of the outermost datapoint.
   if n_elements(PSFsize_) eq 0 then PSFsize = 2. else PSFsize=PSFsize_ ; If PSFsize is not defined, I set it to 2"

   junk_sz=size(radial_bins_user_inputs)
   nbins = n_elements(low_radial_bins_user_inputs)
   low_radial_bins_user_inputs_=fltarr(nbins)
   upper_radial_bins_user_inputs_=fltarr(nbins)
   for i = 0, nbins-1 do begin
      d = execute("low_radial_bins_user_inputs_[i] = "+string(low_radial_bins_user_inputs[i]))
      d = execute("upper_radial_bins_user_inputs_[i] = "+ string(upper_radial_bins_user_inputs[i]))
   endfor
  ; stop
 
   r_bin_lo = low_radial_bins_user_inputs_
   r_bin_up = upper_radial_bins_user_inputs_


   if keyword_set(add_default_bins) then begin

      added_nbins = 6
      if sz[1] lt 50 then added_nbins =5
      if sz[1] lt 15 then added_nbins =3
      if sz[1] lt 10 then added_nbins =2  
   
      junk = mdap_range(1,max(abin)*1.1,added_nbins+1,/log)
      junk[0] = 0
      added_r_bin_up=junk[1:*]
      added_r_bin_lo=junk[0:added_nbins-1]
      r_bin_lo = [r_bin_lo,added_r_bin_lo]
      r_bin_up = [r_bin_up,added_r_bin_up]
      nbins = nbins+added_nbins
   endif

endelse
r_bin=(r_bin_up+r_bin_lo)/2.

;--


;-- start elliptical binning on the r_bins
output_wav = wav
indici_wave_out=indgen(n_elements(wav))
if n_elements(output_lrange) ne 0 then begin
    indici_wave_out = where(wav ge output_lrange[0] and wav le output_lrange[1]) 
    output_wav = wav[indici_wave_out]
endif
n_pixel_wrange_final = n_elements(output_wav)
binned_spectra = dblarr(n_pixel_wrange_final,nbins)
binned_errors = dblarr(n_pixel_wrange_final,nbins)

;basic_colors, black, white, red, green, blue, yellow, cyan, magenta, orange, mint, purple, pink, olive, lightblue, gray   
;window,0,retain=2
;color=[red,green,yellow,blue,red,green,yellow,blue,red,green,yellow,blue,red,green,yellow,blue]
;plot,xbin,ybin,/iso,psym=4
n_elements_bin=fltarr(nbins)
k=0
for i = 0, nbins-1 do begin
   indici = where(abin ge r_bin_lo[i] and abin lt r_bin_up[i])
   if indici[0] eq -1 then continue
   r_bin[i] = median(abin[indici],/even)   ;setting the bin center as the median of semimajor axis
                                           ;of the data within that bin
 ;  oplot,xbin[indici],ybin[indici],color=color[i],psym=7
   
   for kk = 0, n_elements(indici)-1 do begin      ;setting the output binning scheme 2d map
      where_i = where(reference_2dmap eq indici[kk])
      if where_i[0] ne -1 then begin 
         spatial_binning_scheme[where_i] = k
         r2d_bin[where_i]  = r_bin[i]
         r2d_bin_lo[where_i] = r_bin_lo[i]
         r2d_bin_up[where_i] = r_bin_up[i]
      endif
   endfor
   n_elements_bin[i] = n_elements(INDICI)
   if indici[0] eq -1 then n_elements_bin[i] = 0

   ;Adding up spectra within the bin
   spc=total(input_spectra(indici,*),1)
   err=sqrt(total(input_errors(indici,*)^2.,1))
   binned_spectra[*,i] = spc[indici_wave_out]
   binned_errors[*,i] = err[indici_wave_out]
   k=k+1

endfor
indici = where(n_elements_bin ne 0)
binned_spectra = binned_spectra[*,indici]
binned_errors = binned_errors[*,indici]
n_elements_bin = n_elements_bin[indici]
r_bin=r_bin[indici]
r_bin_up=r_bin_up[indici]
r_bin_lo=r_bin_lo[indici]
end_module:

end
 
