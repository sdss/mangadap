; "NEW - V1.3" comments denotes V1.3 modifications. Notably these are:
;
; 1) Keyword FOR_ERRORS have been added and is passed onto GANDALF
; 2) We now call MASK_EMISSION_LINES while specifying the desired
;    rest-frame wavelength to be fitted
; 3) MASK_EMISSION_LINES itself now checks for lines outside this
;    wavelength range, excluding them
; 4) We now keep excluding the Na D region, which can be affected by
;    interstellar absorption
; 5) Set the starting guess velocity for the gas not just to the
;    receding velocity of the stars but also allow for velocity
;    offset, as specfied in the emission-line setup. This makes it
;    easier to fit red or blue wings
; 6) We now append to the solution output the goodpixels array (1
;    good, 0, masked), which can be used later in making plots
;    with SHOW_FIT.PRO

pro remove_indices_from_goodpixels,goodpixels,indici_to_remove_from_good_pixels

for i = 0, n_elements(indici_to_remove_from_good_pixels)-1 do begin
   rm = where(indici_to_remove_from_good_pixels[i] eq goodpixels)
  if rm[0] ne -1 then remove,rm,goodpixels
endfor



end



function mask_emission_lines,npix,Vsys,emission_setup,velscale,l0_gal,lstep_gal,$
                             sigma=sigma,l_rf_range=l_rf_range,log10=log10

; Return a list of goodpixels to fit that excludes regions potentially
; affected by gas emission and by sky lines. Unless the log10 keyword
; is specified, wavelength values are assumed to be ln-rebinned, and
; are defined by the l0_gal, lstep_gal, npix parameters. The position of
; gas and sky emission lines is set by the input emission_setup
; structure and the width of the mask by the sigma parameter. If a
; sigma value is not passed than the width of each line is taken from
; the emission_setup structure.
;
; The rest-frame fitting wavelength range can be manually restricted
; using the l_rf_range keyword to pass min and max observed
; wavelength. Typically used to exclude regions at either side of
; spectra.


; speed of light
c = 299792.458d
; define good pixels array
goodpixels = range(0,npix-1) 
; if set, exclude regions at either ends of the spectra using the keyword l_rf_range
if keyword_set(l_rf_range) then begin
    pix0     = ceil((alog(l_rf_range[0])-l0_gal)/lstep_gal+Vsys/velscale)
    pix1     = ceil((alog(l_rf_range[1])-l0_gal)/lstep_gal+Vsys/velscale)
    if keyword_set(log10) then begin
        pix0     = ceil((alog10(l_rf_range[0])-l0_gal)/lstep_gal+Vsys/velscale)
        pix1     = ceil((alog10(l_rf_range[1])-l0_gal)/lstep_gal+Vsys/velscale)
    endif
    goodpixels = range(max([pix0,0]),min([pix1,npix-1])) ; NEW - V1.3 
endif

tmppixels  = goodpixels
; looping over the listed emission-lines and mask those tagged with an
; 'm' for mask. Mask sky lines at rest-frame wavelength
for i = 0,n_elements(emission_setup.i)-1 do begin
    if (emission_setup.action[i] eq 'm') then begin
        ;print,'--> masking ' + emission_setup.name[i]
        if (emission_setup.name[i] ne 'sky') then $
          meml_cpix = ceil((alog(emission_setup.lambda[i])-l0_gal)/lstep_gal+Vsys/velscale)
        if (emission_setup.name[i] ne 'sky' and keyword_set(log10)) then $
          meml_cpix = ceil((alog10(emission_setup.lambda[i])-l0_gal)/lstep_gal+Vsys/velscale)
        ; sky lines are at rest-frame
        if (emission_setup.name[i] eq 'sky') then $
          meml_cpix = ceil((alog(emission_setup.lambda[i])-l0_gal)/lstep_gal) 
        if (emission_setup.name[i] eq 'sky' and keyword_set(log10)) then $
          meml_cpix = ceil((alog10(emission_setup.lambda[i])-l0_gal)/lstep_gal) 
        ; set the width of the mask in pixels using either
        ; 3 times the sigma of each line in the emission-line setup 
        ; or the provided sigma value 
        if keyword_set(sigma) then msigma = 3*sigma/velscale
        if not keyword_set(sigma) then msigma = 3*emission_setup.s[i]/velscale
        meml_bpix = meml_cpix - msigma
        meml_rpix = meml_cpix + msigma
        w = where(goodpixels ge meml_bpix and goodpixels le meml_rpix) 
        if (w[0] ne -1) then begin
            tmppixels[w] = -1 
        endif else begin 
          ;  print,'this line is outside your wavelength range. We shall ignore it' ; NEW - V1.3
            emission_setup.action[i] = 'i'                                         ; NEW - V1.3 
        endelse
    endif
endfor
w = where(tmppixels ne -1)
goodpixels = goodpixels[w]

return,goodpixels
end

function remouve_detected_emission,galaxy,bestfit,emission_templates,sol_gas_A,AoN_thresholds,$
                                   AoN=AoN,goodpixel=goodpixel
; Given the galaxy spectrum, the best fit, the emission-line
; amplitudes, a vector with the A/N threshold for each line, and the
; array containing the spectra of each best-fitting emission-line
; templates, this function simply compute the residual-noise lvl,
; compares it the amplitude of each lines, and remouves from the
; galaxy spectrum only the best-matching emission-line templates for
; which the correspoding A/N exceed the input threshold.  This is a
; necessary step prior to measurements of the strength of the stellar
; absorption line features
;
; A list of goodpixels may be optionally input, for instance if any
; pixel was excluded by sigma-clipping during the continuum and
; emission-line fitting, or by excluding pixels on either ends of the
; spectra
;
; Also optionally outputs the computed A/N ratios.

; Get the Residual Noise lvl.
resid = galaxy - bestfit
if keyword_set(goodpixel) then resid = resid[goodpixels]
resid_noise = robust_sigma(resid, /ZERO)
; A/N of each line
AoN = sol_gas_A/resid_noise
; Create neat spectrum, that is, a spectrum where only detected
; emission has been remouved
neat_galaxy = galaxy
for i = 0,n_elements(sol_gas_A)-1 do begin
    if (AoN[i] ge AoN_thresholds[i]) then neat_galaxy = neat_galaxy - emission_templates[*,i]
endfor

return,neat_galaxy
end


PRO mdap_gandalf_wrap,templates,loglam_templates,galaxy, loglam_gal, noise,velscale,start_, sol,$
       EMISSION_SETUP_FILE=EMISSION_SETUP_FILE, $
       gas_intens,gas_fluxes,gas_ew,gas_intens_err,gas_fluxes_err,gas_ew_err,$
       BESTFIT=bestFit, BIAS=bias,  MDEGREE=MDEGREE,DEGREE=degree, ERROR=error,$
       MOMENTS=moments, reddening=reddening,err_reddening=err_reddening,$
       VSYST=VSYST, WEIGHTS=weights, BF_COMP2 = bf_comp2,$
       quiet=quiet,FOR_ERRORS=FOR_ERRORS,$
       fix_star_kin=fix_star_kin,fix_gas_kin=fix_gas_kin,$
       range_v_star=range_v_star,range_s_star=range_s_star,range_v_gas=range_v_gas,range_s_gas=range_s_gas,$
       mask_range=mask_range,fitted_pixels=fitted_pixels,external_library=external_library



; Example of an IDL wrapper calling PPXF and GANDALF in order to
; derive the stellar and gaseous kinematics in the case of MANGA
; SPECTRA. The stellar continuum is matched with a combination of
; stellar population models, whereas emission-lines are represented by
; Gaussian templates, with interdepencies regulated by the input
; emission-setup file


c = 299792.4580d ; Speed of light in km/s

; Preliminary Checks !
IF NOT KEYWORD_SET(EMISSION_SETUP_FILE) THEN message,'Please enter the filename with the emission-line setup'
IF NOT KEYWORD_SET(MDEGREE)             THEN mdegree=0


w = where(noise eq 0 or finite(noise) eq 0) & if (w[0] ne -1) then noise[w] = max(noise)*1000.

if n_elements(mask_range) ne 0 then begin 
   indici_to_remove_from_good_pixels = [0];where(
   for i = 0, n_elements(mask_range)-2,2 do begin
      indici_to_remove_from_good_pixels = [indici_to_remove_from_good_pixels,where(loglam_gal ge alog(mask_range[i]) and loglam_gal le alog(mask_range[i+1]))]
   endfor
   indici_to_remove_from_good_pixels = indici_to_remove_from_good_pixels[1:*]
endif


int_disp = 0.
;-


l0_gal   = loglam_gal[0];sxpar(hdr_gal,'CRVAL1')
lstep_gal = loglam_gal[1]-loglam_gal[0];sxpar(hdr_gal,'CD1_1')
l0_templ = loglam_templates[0]
lstep_templ = loglam_templates[1] - loglam_templates[0]




offset = 0.
if n_elements(VSYST) ne 0 then offset = vsyst

; ; F) Preamble to PPXF
; print,'--> and initial guesses for pPXF using the SDSS redshifts'
; ; this is the velocity scale spanned by each log_lambda interval
; ; You need to include a alog(10) factor as well.
; lstep_gal = sxpar(hdr_gal,'CD1_1')
; velscale  = c*lstep_gal*alog(10.0d)
; ; Initial V and sigma guesses, in km/s. For V we use the receiding
; ; velocity derived by the SDSS (in the low-Z regime !!!), plus the
; ; previously derived velocity offset. For sigma we stick at least to
; ; 150 km/s
; SDSS_z = sxpar(hdr_gal,'Z')
; SDSS_s = sxpar(hdr_gal,'VEL_DIS')

; if (SDSS_s gt velscale/2.0d) then $
;   start  = [alog(SDSS_z+1)*c+offset,SDSS_s] else $
;   start  = [alog(SDSS_z+1)*c+offset,150.0d]
SDSS_z = start_[0]/c
start=[alog(SDSS_z+1)*c+offset,start_[1]]
; G) masking emission-lines
; read in emission-line setup files
;print,'--> reading in the emission-line setup file, and creating the corresponding structure'
eml_file = emission_setup_file
readcol,eml_file,eml_i,eml_name,eml_lambda,eml_action,eml_kind,eml_a,eml_v,eml_s,eml_fit,$
  f='(i,a,f,a,a,f,f,f,a)',skipline=2,comment='#',/silent
; creating emission setup structure, to be used for masking the
; emission-line contaminated regions and fitting the emission lines in
; GANDALF
emission_setup = create_struct('i',eml_i,'name',eml_name,'lambda',eml_lambda,'action',eml_action,$
                               'kind',eml_kind,'a',eml_a,'v',eml_v,'s',eml_s,'fit',eml_fit)





; call the function mask_emission_lines, which loops over the listed
; emission-lines and mask those tagged with an 'm' for mask. Use
; emission-line masks with a sigma of 200 km/s, and use only the
; wavelength region corresponding to the rest-frame extent of the
; templates (this relies on the fairly good SDSS guess for the
; redshift z).

;;  ;l_rf_range = 10^[l0_templ,l0_templ+n_pix_templ*sxpar(hdr_templ,'CD1_1')]                       ; NEW - V1.3
;;  l_rf_range = exp(l0_templ,l0_templ+n_pix_templ*sxpar(hdr_templ,'CD1_1'))                       ; NEW - V1.3
n_pix_templ = n_elements(loglam_templates)
lstep_templ = loglam_templates[1]-loglam_templates[0]
;l_rf_range = 10^[l0_templ,l0_templ+n_pix_templ*lstep_templ]                       ; NEW - V1.3
l_rf_range = exp([l0_templ,l0_templ+n_pix_templ*lstep_templ])                       ; NEW - V1.3
;goodpixels = mask_emission_lines(n_elements(galaxy),alog(SDSS_z+1)*c,emission_setup,velscale,$ ; NEW - V1.3
;                                 l0_gal,lstep_gal,sigma=200.0,/log10,l_rf_range=l_rf_range)    ; NEW - V1.3
goodpixels = mask_emission_lines(n_elements(galaxy),alog(SDSS_z+1)*c,emission_setup,velscale,$ ; NEW - V1.3
                                 l0_gal,lstep_gal,sigma=250.0,l_rf_range=l_rf_range)    ; NEW - V1.3
; H) PPXF fit! Fit only V and sigma
; A constant additive polynomial (degree=0) is used in together with
; the multiplicative polynomials (always recommended).

if n_elements(indici_to_remove_from_good_pixels) ne 0 then remove_indices_from_goodpixels,goodpixels,indici_to_remove_from_good_pixels



;stop
; If the galaxy spectrum have zero values in some regions, I increase the noise in
; those regions.
indici_neg = where(galaxy le 0)
if indici_neg[0] ne -1 then noise(indici_neg) = max(noise)
if ~keyword_set(fix_star_kin) then begin
;stop
mdap_ppxf, templates, galaxy, noise, velscale, start, ppxfsol,$
    goodpixels=goodpixels,bias=bias, moments=moments, degree=degree, mdegree=mdegree,$
       range_v_star=range_v_star,range_s_star=range_s_star,ERROR=ERROR_stars,/quiet,bestfit=junk,$
       external_library=external_library
endif else begin
   ppxfsol = [start_[0]+offset,start_[1],start_[2],start_[3],0.,0.,start_[4],start_[5]]
   ERROR_stars = fltarr(6)+99.
endelse
;,VSYST=vsyst


;print,'--> Fit to the stellar continuum masking regions potentially affected by gas emission!'
;if keyword_set(DEBUG) then pause

; J) Preamble to GANDALF
; Switch the tag 'action' from 'm' (for mask) to 'f' (for fit) for all
; lines we wish to fit, i.e, do not touch the ones with an
; 'i'. Likewise keep masking the regions affected by interestella NaD
; absorption and which were originally affected by sky lines.
emission_setup_orig = emission_setup
i_lines = where(emission_setup.action eq 'm' and emission_setup.name ne 'sky' ) ; NEW - V1.3
;i_lines = where(emission_setup.action eq 'm' and emission_setup.name ne 'sky' and emission_setup.name ne 'NaI') ; NEW - V1.3
emission_setup.action[i_lines] = 'f'

; Re-assign the goodpixels array, masking what is still tagged with an 'm'
;goodpixels = mask_emission_lines(n_elements(galaxy),alog(SDSS_z+1)*c,emission_setup,velscale,$
;                                 l0_gal,lstep_gal,sigma=200.0,/log10,l_rf_range=l_rf_range) ; NEW - V1.3
goodpixels = mask_emission_lines(n_elements(galaxy),alog(SDSS_z+1)*c,emission_setup,velscale,$
                                 l0_gal,lstep_gal,sigma=200.0,l_rf_range=l_rf_range) ; NEW - V1.3
;print,'--> lift the emission-line mask...'

; Prepare emission_setup structure for GANDALF, which should only
; deal with the lines we fit
i_f = where(emission_setup.action eq 'f') 
dummy = emission_setup
emission_setup = create_struct('i',dummy.i[i_f],'name',dummy.name[i_f],$
                               'lambda',dummy.lambda[i_f],'action',dummy.action[i_f],$
                               'kind',dummy.kind[i_f],'a',dummy.a[i_f],$
                               'v',dummy.v[i_f],'s',dummy.s[i_f],$
                               'fit',dummy.fit[i_f])
; Assign the stellar systemic velocity as initial guess for the gas
; kinematics, adding also the velocity offset specified in the
; emission-line setup. This makes it easier to add blue or red wings to
; the line profile.
emission_setup.v = ppxfsol[0]-offset + emission_setup.v ; NEW - V1.3
; save the stellar kinematics
sol_star = [ppxfsol[0],ppxfsol[1],ppxfsol[2],ppxfsol[3],0,0,-1]

; H) Call Gandalf, giving it only the stellar kinematics as input
; sol. Now include reddening 

mdegree_=mdegree
if n_elements(reddening) ne 0 then junk = temporary(mdegree_)

;stop
;warning, if the noise vector is not defined, do not run gandalf, but set
;the gandalf outputs to dummy values to have the workflow continue
;without crashing.
if n_elements(indici_to_remove_from_good_pixels) ne 0 then remove_indices_from_goodpixels,goodpixels,indici_to_remove_from_good_pixels
mdap_gandalf, templates, galaxy, noise, velscale, ppxfsol, emission_setup, $
  l0_gal, lstep_gal, GOODPIXELS=goodpixels, INT_DISP=50., $
  BESTFIT=bestfit, EMISSION_TEMPLATES=emission_templates, WEIGHTS=weights, $
  L0_TEMPL=l0_templ,DEGREE=-1, MDEGREE=mdegree_, $
  /FOR_ERRORS, ERROR=esol,$
  REDDENING=REDDENING,$
  fix_gas_kin=fix_gas_kin,$
  range_v_gas=range_v_gas,range_s_gas=range_s_gas,quiet=quiet,$
       external_library=external_library
sol=temporary(ppxfsol)
fitted_pixels=goodpixels
  ;REDDENING=[0.05]
;stop
;*** SONO ARRIVATO QUI, devo aggiustare gli outputs in modo che vengano passati giusti a mdap_spectral_fitting.pro ***
;*** DEVO ANCHE INSERIRE I FIX_STAR_KIN, RANGE_V_STAR, ECC ALL'INTERNO DI MDAP_PPXF E MDAP_GANDALF.PRO

;if keyword_set(DEBUG) then pause
;print,'--> ... and fit simultaneously the stellar continuum and the ionised-gas emission lines!'

; K) Make the unconvolved optimal template. This can be useful for the
; line-strength analysis in order to work out the necessary correction
; due measured indices due to kinematic broadening.
n_templ=n_elements(templates[0,*])
 nweights = weights[0:n_templ-1]/total(weights[0:n_templ-1])
if n_elements(GOODPIXELS) ne 0 then begin
   chi2=total((galaxy[GOODPIXELS]-bestfit[GOODPIXELS])^2/noise[GOODPIXELS]^2)
endif else begin
   chi2=total((galaxy-bestfit)^2/noise^2)
endelse
; otemplate = dblarr(n_elements(templates[*,0]))
; for j=0,n_templ-1 do otemplate = otemplate + templates[*,j]*nweights[j]
;print,'--> Creating the unconvolved optimal template spectrum...'

; L) Call the routine that remouves only the detected emission,
; assuming a constant A/N cut of 4. The routine will also output the
; A/N of each line.

; Best fitting amplitudes
i_l = where(emission_setup.kind eq 'l') 
sol_gas_A = sol[dindgen(n_elements(i_l))*4+1]
; constant A/N=4 threshold
AoN_thresholds = dblarr(n_elements(i_l)) + 4.0
spec_neat = remouve_detected_emission(galaxy,bestfit,emission_templates,sol_gas_A,AoN_thresholds,AoN=sol_gas_AoN)
bf_comp2 = galaxy-spec_neat

;print,'--> ... and cleaning the galaxy spectrum from any detected gas emission line'

; M) Add to the emission setup structure the flux, amplitude and
; kinematics of each line, and call the result the fit_results
; structure. Add to it also the A/N values, the stellar kinematics,
; and the normalised template weights.  This is to save not only the
; emission-line fitting results but also the conditions under which
; the fit was performed.
i_l = where(emission_setup.kind eq 'l',compl=i_d) 
sol_gas_F = sol[dindgen(n_elements(i_l))*4+0]
sol_gas_A = sol[dindgen(n_elements(i_l))*4+1]
sol_gas_V = sol[dindgen(n_elements(i_l))*4+2]
sol_gas_S = sol[dindgen(n_elements(i_l))*4+3]
;

if n_elements(reddening) ne 0 then sol_EBmV  = sol[n_elements(i_l)*4:*] else sol_EBmV =-999
if n_elements(reddening) ne 0 then reddening=sol_EBmV 
dummy = emission_setup
if keyword_set(FOR_ERRORS) then begin
    esol_gas_F = esol[dindgen(n_elements(i_l))*4+0]
    esol_gas_A = esol[dindgen(n_elements(i_l))*4+1]
    esol_gas_V = esol[dindgen(n_elements(i_l))*4+2]
    esol_gas_S = esol[dindgen(n_elements(i_l))*4+3]
    ;
    if n_elements(reddening) ne 0 then esol_EBmV  = esol[n_elements(i_l)*4:*] else esol_EBmV  = -999
    

    fit_results = create_struct('i',dummy.i[i_l],'name',dummy.name[i_l],'lambda',dummy.lambda[i_l],$
                                'action',dummy.action[i_l],'kind',dummy.kind[i_l],'a',dummy.a[i_l],$
                                'v',dummy.v[i_l],'s',dummy.s[i_l],'fit',dummy.fit[i_l],$
                                'Flux',sol_gas_F,'Ampl',sol_gas_A,'Vel',sol_gas_V,'Sigma',sol_gas_S,$
                                'eFlux',esol_gas_F,'eAmpl',esol_gas_A,'eVel',esol_gas_V,'eSigma',esol_gas_S,$
                                'AoN',sol_gas_AoN,'EBmV',sol_EBmV,'eEBmV',esol_EBmV,'Vel_stars',sol_star[0]-offset,$
                                'Sigma_stars',sol_star[1],'Norm_Weights',nweights)
endif else begin
    fit_results = create_struct('i',dummy.i[i_l],'name',dummy.name[i_l],'lambda',dummy.lambda[i_l],$
                                'action',dummy.action[i_l],'kind',dummy.kind[i_l],'a',dummy.a[i_l],$
                                'v',dummy.v[i_l],'s',dummy.s[i_l],'fit',dummy.fit[i_l],$
                                'Flux',sol_gas_F,'Ampl',sol_gas_A,'Vel',sol_gas_V,'Sigma',sol_gas_S,$
                                'AoN',sol_gas_AoN,'EBmV',sol_EBmV,'Vel_stars',sol_star[0]-offset,$
                                'Sigma_stars',sol_star[1],'Norm_Weights',nweights)
endelse
;i_l = where(emission_setup.kind eq 'l')

;sol = [sol_star[0]-offset,sol_star[1:3],chi2,sol_gas_V,sol_gas_S]
;error = [ERROR_stars[0:3],esol_gas_V,esol_gas_S]
gas_intens_ = sol_gas_A
gas_intens_err_ = esol_gas_A 
;gas_fluxes = gas_intens * sol[8]/velscale*sqrt(2.*!pi)
gas_fluxes_ = sol_gas_F
gas_fluxes_err_ = esol_gas_F

where_gas_is_not_null = where(sol_gas_F gt 0 and esol_gas_F gt 0 and finite(sol_gas_F) eq 1 and finite(esol_gas_F) eq 1)

if moments eq 6 then sol = [sol_star[0]-offset,sol_star[1:moments-1],chi2,total(sol_gas_V[where_gas_is_not_null]*gas_fluxes_[where_gas_is_not_null])/total(gas_fluxes_[where_gas_is_not_null]),total(sol_gas_S[where_gas_is_not_null]*gas_fluxes_[where_gas_is_not_null])/total(gas_fluxes_[where_gas_is_not_null])]
if moments eq 6 then error = [ERROR_stars[0:moments-1],mean(esol_gas_V[where_gas_is_not_null]),mean(esol_gas_S[where_gas_is_not_null])]

if moments eq 4 then sol = [sol_star[0]-offset,sol_star[1:moments-1],0,0,chi2,total(sol_gas_V[where_gas_is_not_null]*gas_fluxes_[where_gas_is_not_null])/total(gas_fluxes_[where_gas_is_not_null]),total(sol_gas_S[where_gas_is_not_null]*gas_fluxes_[where_gas_is_not_null])/total(gas_fluxes_[where_gas_is_not_null])]
if moments eq 4 then error = [ERROR_stars[0:moments-1],0,0,mean(esol_gas_V[where_gas_is_not_null]),mean(esol_gas_S[where_gas_is_not_null])]

if moments eq 2 then sol = [sol_star[0]-offset,sol_star[1:moments-1],0,0,0,0,chi2,total(sol_gas_V[where_gas_is_not_null]*gas_fluxes_[where_gas_is_not_null])/total(gas_fluxes_[where_gas_is_not_null]),total(sol_gas_S[where_gas_is_not_null]*gas_fluxes_[where_gas_is_not_null])/total(gas_fluxes_[where_gas_is_not_null])]
if moments eq 2 then error = [ERROR_stars[0:moments-1],0,0,0,0,mean(esol_gas_V[where_gas_is_not_null]),mean(esol_gas_S[where_gas_is_not_null])]
;-- compute gas EW (as in gandalf)
;stop

;l0_gal   = loglam_gal[0];sxpar(hdr_gal,'CRVAL1')
;lstep_gal = loglam_gal[1]-loglam_gal[0];sxpar(hdr_gal,'CD1_1')
;l0_templ = loglam_templates[0]
;lstep_templ = loglam_templates[1] - loglam_templates[0]

rf_l  = exp(loglam_gal -sol[0]/velscale*lstep_gal)
gas_ew_=gas_fluxes_*0.
gas_ew_err_=gas_fluxes_*0.
for j = 0, n_elements(dummy.lambda[i_l])-1 do begin
    rf_l_line  = exp(dummy.lambda[i_l])
    S_line     = sol[6] 
    F_obs_line = gas_fluxes_[j]
    j_buffer   = where(abs(rf_l-rf_l_line) lt 10*(S_line/c)*rf_l_line and abs(rf_l-rf_l_line) gt  5*(S_line/c)*rf_l_line)
    if j_buffer[0] ne -1 then begin
       C_line     = median(spec_neat[j_buffer])
       gas_ew_[j]  = gas_fluxes[j]/C_line
       gas_ew_err_[j]  = gas_fluxes_err[j]/C_line+abs(gas_fluxes[j])/C_line^2.*robust_sigma(spec_neat[j_buffer])
    endif
endfor

gas_intens = [0]
gas_intens_err =[0]
gas_ew = [0]
gas_ew_err =[0]
gas_fluxes = [0]
gas_fluxes_err =[0]
k=0
kk=0
for i = 0,n_elements(eml_lambda)-1 do begin
   
   was_it_in_range = where(abs(emission_setup.lambda - eml_lambda[i]) le 0.001)
   if was_it_in_range[0] eq -1 then begin
      gas_intens = [gas_intens,0]
      gas_intens_err = [gas_intens_err,0]
      gas_fluxes = [gas_fluxes,0]
      gas_fluxes_err = [gas_fluxes_err,0]
      gas_ew = [gas_ew,0]
      gas_ew_err = [gas_ew_err,0]
      continue
   endif

   if emission_setup.KIND[k] eq 'l' then begin
      gas_intens = [gas_intens,gas_intens_[kk]]
      gas_intens_err = [gas_intens_err,gas_intens_err_[kk]]
      gas_fluxes = [gas_fluxes,gas_fluxes_[kk]]
      gas_fluxes_err = [gas_fluxes_err,gas_fluxes_err_[kk]]
      gas_ew = [gas_ew,gas_ew_[kk]]
      gas_ew_err = [gas_ew_err,gas_ew_err_[kk]]
      kk = kk+1
   endif else begin

    ; If the line was a doublet, I assigne
    ; the emission equal to the emission
    ; of the main line in the doubled
    ; times the line ratio
      doublet_of = strmid(emission_setup.KIND[k],1)
      indici = where(emission_setup.i[i_l] eq doublet_of)
      value = emission_setup.A[i]
      gas_intens = [gas_intens,gas_intens_[indici[0]]*value]
      gas_intens_err = [gas_intens_err,gas_intens_err_[indici[0]]*value]
      gas_fluxes = [gas_fluxes,gas_fluxes_[indici[0]]*value]
      gas_fluxes_err = [gas_fluxes_err,gas_fluxes_err_[indici[0]]*value]
      gas_ew = [gas_ew,gas_ew_[indici[0]]*value]
      gas_ew_err = [gas_ew_err,gas_ew_err_[indici[0]]*value]
    endelse
    k = k+1

 
   
endfor
gas_intens = gas_intens[1:*]
gas_intens_err =gas_intens_err[1:*] 
gas_ew = gas_ew[1:*] 
gas_ew_err =gas_ew_err[1:*]
gas_fluxes =gas_fluxes[1:*]
gas_fluxes_err =gas_fluxes_err[1:*]

if n_elements(reddening) ne 0 then reddening = sol_EBmV;,esol_EBmV]
if n_elements(reddening) ne 0 then err_reddening = esol_EBmV


END

