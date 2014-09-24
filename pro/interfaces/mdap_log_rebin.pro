;+
; NAME:
;	MDAP_LOG_REBIN
;
; PURPOSE:
;	TODO
;
; CALLING SEQUENCE:
;	MDAP_LOG_REBIN, flux, ivar, wave, tpl_lib, res_var_diff, flux_log, ivar_log, wave_log, $
;			tpl_lib_log, tpl_lib_wave_log, input_velscale=input_velscale, $
;			wave_range=wave_range, flux=flux, $
;			gal_wavelength_log_step=gal_wavelength_log_step, quiet=quiet, $
;			version=version,
;
; INPUTS:

	flux dblarr[N][C]
		N galaxy spectra, each with C spectral channels

	ivar dblarr[N][C]
		Inverse variance of the N galaxy spectra, each with C spectral channels

	wave dblarr[C]
		The wavelength coordinate of each of the C spectral channels

	tpl_lib string
		Directory containing the template spectra.

library      string          Directory containing the template spectra. Template spectra must be in fits format, defined over a 
                             linear wavelength with constant ang/pix step, and must contain all the header information needed to  
                             recover the wavelenght domain (i.e. CRPIX1, CDELT1, and CRVAL1).

fwhm_diff    dbl array       Vector containing the FWHM(lambda) that the stellar spectra must be convolved for. At the moment, 
                             the value median(fwhm_diff/wavelength*c) is used for broadening. Implementation to use the full 
                             LSF(lambda) information are foreseen.



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





pro mdap_log_rebin,spectra,errors,wavelength,library,fwhm_diff,$
    log_spc,log_err,log_wav,library_log,log_wav_library,$
    input_velscale=input_velscale,wave_range=wave_range,flux=flux,$
    version=version,gal_wavelength_log_step=gal_wavelength_log_step,quiet=quiet
  ;  templ_wavelength_log_step=templ_wavelength_log_step,flux=flux

; ** INPUTS **
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


; Set module
; set wavelength range; using input or existing range




version_module = '0.2'
if n_elements(version) ne 0 then begin
 version = version_module
 goto, end_module
endif


wave_range_=[min(wavelength),max(wavelength)]
if n_elements(wave_range) ne 0 then wave_range_=wave_range


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

;  ; ;-- log rebin of the template spectra
;  ;step_w_gal = (max(wavelength)-min(wavelength))/n_elements(wavelength)
;  fwhm_diff_kms = fwhm_diff/wavelength*c
;  res = file_search(library,count=ntemplates)
;  stddev_diff_pxl = median(fwhm_diff_kms / velscale / fwhm_sigma_conversion,/even) ;ROOM FOR IMPROVEMENT: convolution as function of lambda
;  if stddev_diff_pxl ge 0.0001  then lsf = psf_Gaussian(NPIXEL=8*stddev_diff_pxl, ST_DEV=stddev_diff_pxl, /NORM, NDIM=1)
;  if stddev_diff_pxl lt 0.0001  then lsf = 0
;  for i = 0, ntemplates-1 do begin
;      mdap_readfits_lambda,res[i],spc,lam
;      lamrange=[min(lam),max(lam)]
;      ;mdap_do_log10_rebin,lamrange,spc,spcnew,loglam_templ,velscale=velscale;,flux=flux
;      mdap_do_log_rebin,lamrange,spc,spcnew,loglam_templ,velscale=velscale,flux=flux;,/flux
;      if i eq 0 then library_log = dblarr(n_elements(loglam_templ),ntemplates)
;      ;spcnew_conv = spcnew;convol(spcnew,lsf,/edge_truncate) ;ROOM FOR IMPROVEMENT: convolution as function of lambda
;      if stddev_diff_pxl ne 0 then spcnew_conv = convol(spcnew,lsf,/edge_truncate) ;ROOM FOR IMPROVEMENT: convolution as function of lambda
;      if stddev_diff_pxl eq 0 then spcnew_conv = spcnew
;      library_log[*,i]=temporary(spcnew_conv)
;  endfor
;--
;-- log rebin of the template spectra
;step_w_gal = (max(wavelength)-min(wavelength))/n_elements(wavelength)
;fwhm_diff_kms = fwhm_diff/wavelength*c
res = file_search(library,count=ntemplates)
;stddev_diff_pxl = fwhm_diff_kms/ velscale
for i = 0, ntemplates-1 do begin
    mdap_readfits_lambda,res[i],spc,lam
    lamrange=[min(lam),max(lam)]

    fwhm_diff_ = interpol(fwhm_diff,wavelength,lam);/lam*c
    mdap_do_log_rebin,lamrange,spc,spcnew,loglam_templ,velscale=velscale,flux=flux;,/flux
    

    ;-- trim stars
    ind = where(loglam_templ ge alog(wave_range_[0]-200.) and loglam_templ le alog(wave_range_[1]+15.))
    if ind[0] ne -1 then begin
       loglam_templ = loglam_templ[ind]
       spcnew = spcnew[ind]
    endif
    ;-- 

    if i eq 0 then library_log = dblarr(n_elements(loglam_templ),ntemplates)
    ;-- convolution as function of lambda
    stddev_diff_ang = interpol(fwhm_diff_,alog(lam),loglam_templ)/2.3548
    t0=systime(/seconds)
    spcnew_conv = mdap_convol_sigma(float(exp(loglam_templ)),float(spcnew),float(exp(loglam_templ)),float(stddev_diff_ang));spcnew
    ;spcnew_conv = mdap_convol_sigma(exp(loglam_templ),spcnew,exp(loglam_templ),stddev_diff_ang);spcnew
    ;print, 'star ',mdap_stc(i+1,/integer),' of ',mdap_stc(ntemplates,/integer), '   (N=',mdap_stc(n_elements(spcnew_conv),/integer),' pixels), secs=',mdap_stc(systime(/seconds)-t0)
    library_log[*,i]=temporary(spcnew_conv)
    
endfor
log_wav_library=temporary(loglam_templ)
;--

;-- trim galaxy

;ind = where(loglam ge alog10(wave_range_[0]) and loglam le alog10(wave_range_[1]))
ind = where(loglam ge alog(wave_range_[0]) and loglam le alog(wave_range_[1]))
log_spc = log_spc[ind,*]
log_err = log_err[ind,*]
log_wav = loglam[ind]
;--
;stop
end_module:
end
