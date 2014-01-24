;######################################################################
;
; Copyright (C) 2001-2003, Michele Cappellari
; E-mail: cappellari@strw.leidenuniv.nl
;
; This software is provided as is without any warranty whatsoever.
; Permission to use, for non-commercial purposes is granted.
; Permission to modify for personal or internal use is granted,
; provided this copyright and disclaimer are included unchanged
; at the beginning of the file. All other rights are reserved.
;
;######################################################################
;
;
;----------------------------------------------------------------------
pro mdap_do_log10_rebin, lamRange, spec, specNew, logLam, $
    OVERSAMPLE=oversample, VELSCALE=velScale
compile_opt idl2
on_error, 2

if n_elements(lamRange) ne 2 then message, 'lamRange must contain two elements'
s = size(spec)
if s[0] ne 1 then message, 'input spectrum must be a vector'
n = s[1]
if n_elements(oversample) ne 0 then m = n*oversample else m = n

dLam = (lamRange[1]-lamRange[0])/(n-1d)          ; Assume constant dLam
lim = lamRange/dLam + [-0.5d,0.5d]               ; All in units of dLam
borders = range(lim[0],lim[1],n+1)               ; Linearly
logLim = alog10(lim)

c = 299792.458d                                  ; Speed of light in km/s
if n_elements(velScale) gt 0 then begin          ; Velocity scale is set by user
    logScale = velScale/c/alog(10.0d)
    m = floor((logLim[1]-logLim[0])/logScale)    ; Number of output pixels
    logLim[1] = logLim[0] + m*logScale
endif else $
    velScale = (logLim[1]-logLim[0])/m*c*alog(10.0d) ; Only for output

newBorders = 10^(range(logLim[0],logLim[1],m+1)) ; Logarithmically
k = floor(newBorders-lim[0]) < (n-1) > 0
if n_params() ge 4 then $  ; Optional output log(wavelength): log of geometric mean!
    logLam = alog10(sqrt((newBorders*shift(newBorders,1))[1:*])*dLam)

specNew = dblarr(m,/NOZERO)
for j=0,m-1 do begin                             ; Do analytic integration
    a = newBorders[j] - borders[k[j]]
    b = borders[k[j+1]+1] - newBorders[j+1]
    specNew[j] = total(spec[k[j]:k[j+1]]) - a*spec[k[j]] - b*spec[k[j+1]]
endfor

end


pro mdap_log_rebin,spectra,errors,wavelength,library,fwhm_diff,$
    log_spc,log_err,log_wav,library_log,log_wav_library,$
    input_velscale=input_velscale,wave_range=wave_range,flux=flux,$
    version=version,gal_wavelength_log_step=gal_wavelength_log_step,quiet=quiet
  ;  templ_wavelength_log_step=templ_wavelength_log_step,flux=flux

; ** INPUTS **
;spectra      NxT dbl array   The N galaxy spectra to resample.
;
;errors       NxT dbl array   Errors associated to the N spectra.
;
;wavelength   dbl array       Array with T elements, that specifies the wavelengths where galaxy spectra and errors are defined;
;
;library      string          Directory containing the template spectra. Template spectra must be in fits format, defined over a 
;                             linear wavelength with constant ang/pix step, and must contain all the header information needed to  
;                             recover the wavelenght domain (i.e. CRPIX1, CDELT1, and CRVAL1).
;
;fwhm_diff    dbl array       Vector containing the FWHM(lambda) that the stellar spectra must be convolved for. At the moment, 
;                             the value median(fwhm_diff/wavelength*c) is used for broadening. Implementation to use the full 
;                             LSF(lambda) information are foreseen.
;
; ** OPTIONAL INPUTS **
;input_velscale  flt          Constant km/sec/pixel to be used when rebinning the input spectra. If not provided, the value will be
;                             automatically set by the procedure.
;
;wave_range   2 elem array    If specified, the galaxy spectra will be trimmed to this wavelength range (units: angstrom). Default: 
;                             use the entire input wavelength range. Stellar spectra will be trimmed by 
;                             wave_range[0] - 250 ang and wave_range[1] + 250 ang.
;
; ** OPTIONAL KEYWORDS **
;
; /flux                        If set, flux conservation is applied to the log resampling. **Do not use** for template fitting.
; /gal_wavelength_log_step &  Set this keyword if the input galaxy spectra are logarithmically sampled (i.e. wavelength has a 
;                             logarithmic progression).
; /quiet                       If set, message prompt is suppressed.
;
; ** OUTPUTS **
;
;log_spc    N x TT dbl array  The logarithmically resampled (ln-lambda) N galaxy spectra, over the wavelength range ln(wave_range).
;
;log_err    N x TT dbl array  The errors associated to the log_spc. Errors are rebinned using the following formulas:
;                                 lamrange=minmax(wavelength)
;                                 mdap_do_log_rebin,lamrange,errors^2,log_err2,loglam,velscale=velscale
;                                 log_err = sqrt(log_err2)
;                             where mdap_do_log_rebin.pro is the original procedure by M. Cappellari.
;
;log_wav    TT dlb array      The values of the ln-lambda over which log_spc and log_err are defined.
;
;library_log  W x M dbl array The W stellar template spectra, log resampled.
;
;log_wav_library M dbl array  The values of log-lambda over which the stellar templates are defined.

;*** OPTIONAL OUTPUTS ***
;
; version  [string]            Module version. If requested, the module is not execute and only version flag is returned
;

version_module = '0.1'
if n_elements(version) ne 0 then begin
 version = version_module
 goto, end_module
endif




c = 299792.458d                                  ; Speed of light in km/s

fwhm_sigma_conversion = double(2.*sqrt(alog(4.)))
sz = size(spectra)
if sz[0] eq 1 then nspec = 1
if sz[0] eq 2 then nspec = sz[2]


step_w_gal = wavelength[1]-wavelength[0]
;-- log rebin of the input galaxy spectra and errors
if n_elements(input_velscale) ne 0 then velscale=input_velscale
lamrange=[min(wavelength),max(wavelength)]
;print, 'velscale=',velscale
;loglam = alog(wavelength)  ;  alog10(wavelength)  
FOR i = 0, nspec-1 DO BEGIN

   if not keyword_set(gal_wavelength_log_step) then begin ; galaxy spectra have constant ang/pxl step
     ; mdap_do_log10_rebin,lamrange,spectra[*,i],spcnew,loglam,velscale=velscale,flux=flux
      mdap_do_log_rebin,lamrange,spectra[*,i]*step_w_gal,spcnew,loglam,velscale=velscale,flux=flux
      if i eq 0 then begin 
         log_spc = dblarr(n_elements(loglam),nspec)
         log_err = dblarr(n_elements(loglam),nspec)
      endif
      log_spc[*,i] = temporary(spcnew)
      
      ;mdap_do_log10_rebin,lamrange,errors[*,i]^2.,err2new,loglam,velscale=velscale,flux=flux ;
      mdap_do_log_rebin,lamrange,errors[*,i]^2.,err2new,loglam,velscale=velscale,flux=flux ;
      ind_null = where(finite(err2new) eq 0,complement=ok)
      if ind_null[0] ne -1 then err2new[ind_null] = max(err2new[ok])
      log_err[*,i] = sqrt(temporary(err2new))   
   endif else begin; galaxy spectra have logarithimc ang/pxl step
      loglam = alog(wavelength) ;  alog10(wavelength)  
      tmp_spc = spectra[*,i]
      tmp_err = errors[*,i]
      now_velscale_is_defined:
      if n_elements(velscale) ne 0 then begin
         ;if velscale ne (loglam[1]-loglam[0])*c*alog(10.d) then begin 
         if velscale ne (loglam[1]-loglam[0])*c then begin 
            if ~keyword_set(quiet) and i eq 0 then print, 'changing velscale from '+mdap_stc((loglam[1]-loglam[0])*c)+' km/sec/pxl to '+mdap_stc(velscale)+ ' km/sec/pxl';
            ;log_scale = double(velScale/c/alog(10.0d))
            log_scale = double(velScale/c)
            npix = round((max(loglam)-min(loglam))/log_scale)+1
            loglam_new = dindgen(npix) *log_scale + min(loglam)
            tmp_spc = interpol(tmp_spc,loglam,loglam_new)
            tmp_err = interpol(tmp_err,loglam,loglam_new)
            loglam=loglam_new
            velscale = (loglam[1]-loglam[0])*c
;            print, i
         endif

      ;endif else velscale = (loglam[1]-loglam[0])*c*alog(10.d)  
      endif else begin 
         ;velscale = (loglam[1]-loglam[0])*c*alog(10.d)  
         velscale = (loglam[1]-loglam[0])*c
         if ~keyword_set(quiet) then print, 'velscale not defined.... doing it now'
         goto, now_velscale_is_defined
      endelse

     ; stop
      if i eq 0 then begin 
         log_spc = dblarr(n_elements(loglam),nspec)
         log_err = dblarr(n_elements(loglam),nspec)
      endif
      log_spc[*,i] = tmp_spc
      log_err[*,i] = tmp_err
   endelse

   
ENDFOR
;--

;-- log rebin of the template spectra
;step_w_gal = (max(wavelength)-min(wavelength))/n_elements(wavelength)
fwhm_diff_kms = fwhm_diff/wavelength*c
res = file_search(library,count=ntemplates)
stddev_diff_pxl = median(fwhm_diff_kms / velscale / fwhm_sigma_conversion,/even) ;ROOM FOR IMPROVEMENT: convolution as function of lambda
if stddev_diff_pxl ge 0.0001  then lsf = psf_Gaussian(NPIXEL=8*stddev_diff_pxl, ST_DEV=stddev_diff_pxl, /NORM, NDIM=1)
if stddev_diff_pxl lt 0.0001  then lsf = 0
for i = 0, ntemplates-1 do begin
    mdap_readfits_lambda,res[i],spc,lam
    lamrange=[min(lam),max(lam)]
    ;mdap_do_log10_rebin,lamrange,spc,spcnew,loglam_templ,velscale=velscale;,flux=flux
    mdap_do_log_rebin,lamrange,spc,spcnew,loglam_templ,velscale=velscale,flux=flux;,/flux
    if i eq 0 then library_log = dblarr(n_elements(loglam_templ),ntemplates)
    ;spcnew_conv = spcnew;convol(spcnew,lsf,/edge_truncate) ;ROOM FOR IMPROVEMENT: convolution as function of lambda
    if stddev_diff_pxl ne 0 then spcnew_conv = convol(spcnew,lsf,/edge_truncate) ;ROOM FOR IMPROVEMENT: convolution as function of lambda
    if stddev_diff_pxl eq 0 then spcnew_conv = spcnew
    library_log[*,i]=temporary(spcnew_conv)
endfor
;--

;-- trim galaxy
wave_range_=[min(wavelength),max(wavelength)]
if n_elements(wave_range) ne 0 then wave_range_=wave_range

;ind = where(loglam ge alog10(wave_range_[0]) and loglam le alog10(wave_range_[1]))
ind = where(loglam ge alog(wave_range_[0]) and loglam le alog(wave_range_[1]))
log_spc = log_spc[ind,*]
log_err = log_err[ind,*]
log_wav = loglam[ind]
;--
;-- trim stars
;ind = where(loglam_templ ge alog10(wave_range_[0]-250.) and loglam_templ le alog10(wave_range_[1]+250.))
ind = where(loglam_templ ge alog(wave_range_[0]-250.) and loglam_templ le alog(wave_range_[1]+250.))
library_log = library_log[ind,*]
log_wav_library = loglam_templ[ind]
;--
;stop
end_module:
end
