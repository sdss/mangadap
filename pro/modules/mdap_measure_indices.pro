pro mdap_measure_indices,wavelength,spectra, best_template,best_template_LOSVD,$
         stellar_velocity,residuals,fwhm_diff_indices_,abs_line_indices,$
         abs_line_indices_errors,abs_line_indices_template,abs_line_indices_template_losvd,$
         dir=dir,version=version,noplot=noplot,remove_outliers=remove_outliers


; INPUTS
; wavelength          [N dblarray].     Vector containing the wavelenghts of the input spectra. The dispersion can be also not constant.
; 
; spectra             [T x N dblarray]. Vector containing the T input galaxy spextra, with emission line removed. Spectra are defined over the vector wavelength.
;
; best_template       [T x N dblarray]. Vector containing the T best fitting stellar templates obtained when fitting the kinematics of the input spectra. 
;                                       Spectra are defined over the vector wavelength.
;
; best_template_LOSVD [T x N dblarray]. Vector containing the T best fitting stellar templates, convolved with the best-fitting stellar LOSVD, obtained when 
;                                       fitting the kinematics of the input spectra. Spectra are defined over the vector wavelength. 
;
; stellar_velocity    [T flt array].    Vector containing the best fitting stellar velocity for the T input spectra in km/sec. 
;                                       P.S. Set it to zero if the input spectra are in rest-frame.

; residuals           [T x N dblarray]. Vector containing the residuals from the best fit model to the input galaxy spectra. 
;                                       Residuals are defined over the vector wavelength.  
;
; fwhm_diff_indices_  [N dblarray].     Vector that specifies the FWHM (in angstrom) that should be used to broaden the spectra, 
;                                       best_template, and best_templat\_LOSVD to match the spectral resolution of the spectral indices system.\\
;
; OPTIONAL INPUTS
; dir=dir             [string]          Directory where to store the ps files showing the measurements 
;
; OPTIONAL INPUTS     integer   
;      
; remove_outliers     integer          If defined and positive, pixels that deviate more than remove_outliers[*] * rms will be
;                                      substituted by the best fit
;                                      model. Default: 0
;
; OUTPUTS
; abs_line_indices    [Tx37 dbl array]  Absorption line indices of the T inut spectra, corrected for intrinsic broadening. The measured 37 
;                                       indices are defined in Table \ref{dap_tab:indices}.\\
;
; abs_line_indices_errors          [T x 37 dbl array]. Errors associated to abs\_line\_ indices. 
;
; abs_line_indices_template        [T x 37 dbl array]Absorption line indices measured on best\_template.
;
; abs_line_indices\_template_losvd [T x 37 dbl array]Absorption line indices measured on best\_template\_ LOSVD.
;
;  OPTIONAL OUTPUTS
; version              [ string] It specifies the module version. If requested, the module is not execute and only version flag is returned.\\


version_module = '0.1'
if n_elements(version) ne 0 then begin
 version = version_module
 goto, end_module
endif

output_dir=''
if n_elements(dir) ne 0 then output_dir = dir

lambda_step = wavelength[1] - wavelength[0]
fwhm_diff_indices=fwhm_diff_indices_/lambda_step  ;from angstrom to pixels
fwhm_sigma_conversion = double(2.*sqrt(alog(4.)));2.35482
sigma = fwhm_diff_indices/fwhm_sigma_conversion
c=299792.48d
sz = size(spectra)
mdap_read_indices_definitions,indices=indices

;definition of outputs   spc        final value         
abs_line_indices=fltarr(sz[1],n_elements(indices.name))/0.;-999.
abs_line_indices_errors=fltarr(sz[1],n_elements(indices.name))/0.;-999.
abs_line_indices_template=fltarr(sz[1],n_elements(indices.name))/0.;-999.
abs_line_indices_template_losvd=fltarr(sz[1],n_elements(indices.name))/0.;-999.

if ~keyword_set(noplot) then set_plot,'ps'
FOR i = 0, sz[1]-1 DO BEGIN
   
   spc = spectra[i,*]
   res = residuals[i,*]
   template = best_template[i,*]
   template_LOSVD = best_template_LOSVD[i,*]
   v = stellar_velocity[i]

   
    ;--- noise computation from rms vector
    
  
    nbins = fix(n_elements(spc)/20)
   ; stop
    bins = range(0, n_elements(spc),nbins+1)
    fake_x = findgen(n_elements(spc))
    rms_per_bin = fltarr(nbins)
    bin_center = fltarr(nbins)
     
    for kk = 0, nbins-1 do begin
       indici = where(fake_x ge bins[kk] and fake_x lt bins[kk+1])
       ;rms_per_bin[kk] = sqrt(stddev(residuals[indici])^2 + median(residuals[indici])^2)  ;convert systematics into rms...not sure if it is the correct things to do
       rms_per_bin[kk] = stddev(res[indici])  
       bin_center[kk] = (bins[kk+1]+ bins[kk])*0.5
    endfor
    noise = interpol(rms_per_bin,bin_center,findgen(n_elements(spc)))
    ;---
   
   ;  ;--identification of 3-sigma outliers and replacement with stellar_best-fit_model. 
    ksigma = -1
    if n_elements(remove_outliers) ne 0 then ksigma = remove_outliers[0]

    if ksigma gt 0 then begin
       for k = 0, n_elements(indices.name)-1 do begin
          blue_cont = (1.d +v/c) * indices[k].blue_cont
          red_cont = (1.d +v/c) * indices[k].red_cont
          ind_to_clean = where(wavelength ge blue_cont[0] and wavelength le red_cont[1])
          if ind_to_clean[0] ne -1 then begin
             star_best = spc[ind_to_clean] - res[ind_to_clean] 
             ind_outliers = where(abs(res[ind_to_clean]) ge ksigma*robust_sigma(res[ind_to_clean] ),comple=no_outliers )
             if ind_outliers[0] ne -1 and n_elements(no_outliers) ge n_elements(ind_outliers)  then spc[ind_to_clean[ind_outliers]] = star_best[ind_outliers]
          endif
     endfor
    endif
   ;--

     ; these lines triggers the sigma clipping procedure, which is done index by index (to have a better handle of the RMS around the index bandpass).
   ;--

   ;--convolution to match Lick system. NEEDS TO BE OPTIMIZED.
   if max(sigma) le 0 then begin
      x_pixels=findgen(n_elements(wavelength))
      spc_conv=spc;mdap_convol_sigma(x_pixels,spc,x_pixels,sigma) 
      template_conv = template;mdap_convol_sigma(x_pixels,template,x_pixels,sigma) 
      template_LOSVD_conv = template_LOSVD;mdap_convol_sigma(x_pixels,template_LOSVD,x_pixels,sigma) 
   endif else begin
      x_pixels=findgen(n_elements(wavelength))
      spc_conv=mdap_convol_sigma(x_pixels,spc,x_pixels,sigma) 
      template_conv = mdap_convol_sigma(x_pixels,template,x_pixels,sigma) 
      template_LOSVD_conv = mdap_convol_sigma(x_pixels,template_LOSVD,x_pixels,sigma) 
   endelse
   ;--
   if ~keyword_set(noplot) then device, filename=output_dir+'indices_spectrum_'+mdap_stc(i+1,/integer)+'.ps',/color,xsize=15,ysize=12


  

   FOR k = 0, n_elements(indices.name)-1 DO BEGIN



   ;bandpass velocity shift 
   ;WARNING: relativistic correction needed.
   passband = (1.d +v/c) * indices[k].passband
   blue_cont = (1.d +v/c) * indices[k].blue_cont
   red_cont = (1.d +v/c) * indices[k].red_cont
   
   ;measurement of k-th index on i-th spectrum
   ;check if index bandpasses are within the wav-range  
   if indices[k].blue_cont[0] ge min(wavelength) and red_cont[1] le max(wavelength) then begin
    ;  print, 'measuring ',indices[k].name

  
      ;.. actual computation module
  ;stop
        mdap_do_measure_index,spc_conv,wavelength,passband ,blue_cont,$
                red_cont,equiv_w1,index_mag1,noise=noise,error=error,title='galaxy, spc '+mdap_stc(i+1,/integer)+' index: '+indices[k].name,noplot=noplot
       
        mdap_do_measure_index,template_conv,wavelength,indices[k].passband ,indices[k].blue_cont,$
                indices[k].red_cont,equiv_w2,index_mag2,title='template, spc '+mdap_stc(i+1,/integer)+' index: '+indices[k].name,noplot=noplot
       
        mdap_do_measure_index,template_LOSVD_conv,wavelength,indices[k].passband ,indices[k].blue_cont,$
                indices[k].red_cont,equiv_w3,index_mag3,title='template_LOSVD, spc '+mdap_stc(i+1,/integer)+' index: '+indices[k].name,noplot=noplot
       
      ;..  
        
        if indices[k].units eq 'mag' then correction = index_mag2[0] - index_mag3[0]
        if indices[k].units eq 'ang' then correction = equiv_w2[0] / equiv_w3[0]

        if indices[k].units eq 'mag' then begin
           abs_line_indices[i,k] = index_mag1 + correction[0]
           abs_line_indices_template[i,k] = index_mag2[0]
           abs_line_indices_template_losvd[i,k] = index_mag3[0]
           abs_line_indices_errors[i,k]=error[1] 
        endif

        if indices[k].units eq 'ang' then begin
           ind = equiv_w1 / correction[0]
        
           abs_line_indices[i,k] = equiv_w1 / correction[0]
           abs_line_indices_template[i,k] = equiv_w2[0]
           abs_line_indices_template_losvd[i,k] = equiv_w3[0]
           abs_line_indices_errors[i,k]=error[0] 
        endif


     endif else begin
      ;  print, 'skipping ',indices[k].name
     endelse

  

   ENDFOR
   if ~keyword_set(noplot) then device, /close
 ;  print, 'Abs indices measured in spectrum '+output_dir+mdap_stc(i+1,/integer)+'/'+mdap_stc(sz[1],/integer)
ENDFOR
if ~keyword_set(noplot) then set_plot, 'x'
end_module:
end
