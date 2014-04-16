pro mdap_spectral_fitting,galaxy,noise,loglam_gal,templates,loglam_templates,velscale,$
     stellar_kinematics,stellar_kinematics_err,stellar_weights,emission_line_kinematics,emission_line_kinematics_err,$
      emission_line_kinematics_individual,emission_line_kinematics_individual_err,$
     emission_line_intens,emission_line_intens_err,emission_line_fluxes,emission_line_fluxes_err,emission_line_equivW,emission_line_equivW_err,wavelength_input=wavelength_input,$
     wavelength_output,best_fit_model,galaxy_minus_ems_fit_model,best_template,best_template_LOSVD_conv,reddening_output,reddening_output_err,residuals,$
     star_kin_starting_guesses=star_kin_starting_guesses,gas_kin_starting_guesses=gas_kin_starting_guesses,$
     emission_line_file=emission_line_file,quiet=quiet,$
     extra_inputs=extra_inputs,use_previos_guesses=use_previos_guesses,$
     fix_star_kin=fix_star_kin,fix_gas_kin=fix_gas_kin,$
     range_v_star=range_v_star,range_s_star=range_s_star,range_v_gas=range_v_gas,range_s_gas=range_s_gas,rest_frame_log=rest_frame_log,$
     version=version,mask_range=mask_range,$
     external_library=external_library,fwhm_instr_kmsec_matrix=fwhm_instr_kmsec_matrix
;,$
    ; emission_line_intens=emission_line_intens,emission_line_intens_rms=emission_line_intens_rms;,$
    ; best_fit_model_log=best_fit_model_log,emission_model_log=emission_model_log,$
    ; optimal_template_log=optimal_template_log


; *** INPUTS ***
; galaxy     [MxN  array]  It contains the N galaxy spectra to fit, logarithmically 
;            sampled (natural log). Units: 1e-17 erg/s/cm^2/Angstrom. Spectra are defined of the
;            M elements array loglam_gal (see below).
;
; noise      [MxN array]  It contains the N error vectors for the
;            N galaxy spectra. Same units as galaxy. Error vectors are defined of the
;             M elements array loglam_gal (see below).
;
; loglam_gal [M dblarray] It contains the log wavelength values where
;            galaxy and noise spectra are sampled.The vector must have a constant 
;            log(angstrom) sampling. 
;
; templates  [MM x NN array]. It contains the NN stellar template
;            spectra, logarithmically sampled at the same km/sec/pixel as the
;            galaxy spectra. Same units as galaxy, except an arbitrary
;            multiplicative factor.
;
; LOGLAM_TEMPLATES [MM dblarray] It contains the log wavelength values where
;            templates are sampled. It must have a constant log(angstrom) sampling.
;
; velscale   sampling of the input spectra, in km/sec/pixel.
;
;
; *** OPTIONAL INPUTS ***
   

; extra_inputs  [string array]  A string array containing other inputs
;               that might be used in the fitting procedure, such as
;               the number of polinomyal degree. Variable will be
;               initialized with the IDL execute command.
;                    for i = 0, n_elements(extra_inputs)-1 do d = execute(extra_inputs[i])
;               EXAMPLE:  extra_inputs=['MOMENTS=2','DEGREE=-1','BIAS=0','reddening=0','LAMBDA=exp(loglam_gal)']       
;
;                Warning: The reddening (stars and/or gas) fit is
;                 optional, and it is performed by the Gandalf module. If
;                 the reddening fit is required, MDEGREE and DEGREE are
;                 used as the input value for the pPXF run, but
;                 automatically set to 0 and -1 respectively in the
;                 Gandalf execution.
; 
; star_kin_starting_guesses [ N x 4 fltarray] The stellar kinematics starting guesses for V, 
;                           sigma, H3, and H4 for the N galaxy spectra to fit.
;                           Default values are 0. for V, H3, and H4, and 50 km/sec for sigma. 
;                           Starting guesses values are overrided by
;                           the \use_previos_guesses keyword, if set.
;
; gas_kin_starting_guesses [ N x 2 fltarray] The emission line kinematics starting guesses for V, 
;                           sigma, for the N galaxy spectra to fit.
;                           Default values are 0 km/sec for V, and 50 km/sec for sigma. 
;                           Starting guesses values are overrided by
;                           the \use_previos_guesses keyword, if set.
;
; emission_line_file [string] string containing the name of the file with the information of the
;                    emission lines to be fitted.
;                    The format must be compatible with the module that performs the fit of the emission lines. 
;                    If Gandalf or S-Gandalf are used, the input file must be an ascii file with the 
;                    following 9 columns as in the following example (comments starts with ''\#''):
; 
;                     #  ID     CODE    wavelength   action  line/    Intensity	V_g/i	sig_g/i	fit-kind
;                     #                 [angstrom]   f/i/m   doublet                           f/tN
;                     #   0	HeII   3203.15       m	      l         1.000	   0	10	t25
;                     #   1     NeV    3345.81       m        l         1.000      0    10      t25
;                         2     NeV    3425.81       m        l         1.000      0    10      t25
;                         3	OII    3726.03       m	      l 	1.000	   0	10	t25
;               
;                 Col 1: ID. Unique integer identifyer of the emission line.
;
;                 Col 2: CODE. String. Name of the element. Special characters as ``[``, ``.'' are not permitted.;
;
;                 Col 3: wav. Float. Rest frame wavelength of the emission line to fit. Warning, the
;                             name of the emission line in the final DAP output will be defined by the string: CODE  +''_''+ROUND(wav).
;
;                 Col 4: action. String. Possible values are ``f'' (fit), ``m'' (mask), and ``i''
;                          (ignore). All lines are masked in the pPXF run, and then fit in the Gandalf run,
;                          unless they are marked with ``i'', in this latter case they will be ignored.;
;
;                 Col 5: line. String Possible values are ``l'' (line) or ``dN'' (doublet). In this
;                              case, N indicates the line ID to which the  doublet is linked, linked to the line with
;                              ID=N. For example, if emission line with ID=4 has line=d3, then the 
;                              emission line with ID=3 must have line = l.
;
;                 Col 6: int. Float. Relative intensity of the gas emission (positive) or absorption 
;                             (negative) lines. It is set to 1 for lines (line = ``l''). For doublets, (line = ``dN'') 
;                             it indicates the ration elements of the doublets.
;
;                 Col 7: vel. Float. Guess of the velocity offset with respect the galaxy systemic velocity (ignored in the DAP)
;
;                 Col 8: sigma. Float. Guess of the Velocity dispersion (ignored in the DAP)
;
;                 Col 9: mode. String. Possible values are ``f'' (fit) or ``tN'' (tied). If a line has mode = tN, its kinematics is tied 
;                                       to the line with ID=N. In this case, the line with 00=N must have mode = f.
;
;
; range_v_star  [2 elements array]. It specifies the boundaries for the stellar best fit velocity (in km/sec). Default: starting_guess +/- 2000 km/sec.
;
; range_s_star  [2 elements array]. It specifies the boundaries for the stellar best fit velocity dispersion (in km/sec). Default: 21 < sigma < 499 km/sec.
;
; range_v_gas   [2 elements array]. It specifies the boundaries for the emission line best fit velocity (in km/sec). Default: starting_guess +/- 2000 km/sec.
;
; range_s_gas   [2 elements array]. It specifies the boundaries for the emission line best fit velocity dispersion (in km/sec). Default: starting_guess +/- 2000 km/sec.
;
; wavelength_input [QQ elements array]. If specified, it will be used to create wavelength_output, i.e. the wavelength
;                  vector (constant ang/pixel step, in linear units) to interpolate the final results on. 
;                  if keyword /rest_frame_log is set, the vector is set to exp(loglam_templates), and user inpiut will be overwritten
;                  In this case QQ = MM. The default is to use the wavelenth vector specified by the template stars.
;
; external_library [string] String that specifies the path to the external FORTRAN library, which contains the fortran versions of mdap_bvls.pro. 
;                  If not specified, or if the path is invalid, the default internal IDL mdap_bvls code is used. 
;
; mask_range   If defined, it specifies the wavelength ranges to avoid
;            in the fit. It must contain an even number of entries, in
;           angstrom. E.g. l0,l1,l2,l3,....l(2n-1),l(2n) will mask all the pixels where the 
;             l0 < exp(loglam_gal) <l1, or l2 < exp(loglam_gal) <l3, or l(2n-1) < exp(loglam_gal) <l(2n)
;
; fwhm_instr_ kmsec_matrix  2xW elements vector. It defines the instrumental FWHM as function of wavelength. 
;                               fwhm_instr_kmsec_matrix[0,*] specifies the wavelength (in angstrom, at which the instrumental FWHM is measured. 
;                               fwhm_instr_kmsec_matrix[1,*] specifies values of the instrumental FWHM (in km/sec) measured at fwhm_instr_kmsec_matrix[0,*]. 
;                               If undefined, a constant instrumental FWHM=122.45 km/sec is adopted. \\
       
; *** OPTIONAL KEYWORDS ***
;
; /use_previos_guesses  If set, the starting guesses for spectrum i-th
;                       will be the best fit values from spectrum
;                      (i-1)-th (i>0). Input starting guesses will be ignored.
;
;
; /fix_star_kin         If set, the stellar kinematics are not
;                       fitted. The return value is that of the starting guesses 
;
; /fix_gas_kin          If set, the emission-lines kinematics are not
;                       fitted. The return value is that of the starting guesses 
;
; /quiet                If set, some information are not printed on screen
;
;
; /rest_frame_log  If set, the output spectra (galaxy_minus_ems_fit_model, best_fit_model, residuals, 
;                  best_template, and best_template_LOSVD_conv) are produced at rest-frame wavelength.
;
; *** OUTPUTS ***
;
;  stellar_kinematics     [N x 5 flt array]  It contains the best fit values of V, sigma, h3, h4, and chi2/DOF 
;                         for each of the N fitted input galaxy spectra. If \fix_star_kin is set, the array is not defined.
;
;  stellar_kinematics_err [N x 5 flt arrary]  It contains the errors to the best fit values of V, sigma, h3, h4, and chi2/DOF 
;                         for each of the N fitted input galaxy spectra
;
;  stellar_weights        [N x NN dbl array]. It contains the weights of the NN templates for each of the N input galaxy spectra
;
;
;
;  emission_line_kinematics [N x 2 flt array]  It contains the best fit values of V, sigma (emission lines) for each of the N fitted
;                            input galaxy spectra. If \fix_gas_kin is set, the array is not defined.
;
;
; emission_line_kinematics_err [N x 2 flt array]  It contains the errors to the best fit values of V, sigma (emission lines) for each of the N fitted
;                            input galaxy spectra. If \fix_gas_kin is set, the array is not defined.
;
;emission_line_kinematics_individual,emission_line_kinematics_individual_err
; kinematics and errors of individual emission lines.
; emission_line_fluxes  [N x T flt array]  It contains the fluxes of the T fitted emission lines for each of the N input galaxy spectra. 
;                       Values are corrected for reddening
;
; emission_line_fluxes_err  [N x T flt array]  Errors associated to emission_line_fluxes
;
;
; emission_line_equivW      [N x T flt array]  It contains the Equivalent widths of the T fitted emission lines for each of the N input galaxy spectra.
;                           Equivalent widths are computed by the ratio of emission_line_fluxes and the median value of the stellar spectrum within 5 and
;                           10 sigma from the emission line.Sigma is the emission line velocity dispersion.
; 
; emission_line_equivW_err  [N x T flt array] Errors associated to emission_line_equivW  
;

; wavelength_output     [QQ elements flt array] It will contain the linear wavelength values over which the output spectra are sampled. 
;                       Default: it is set to wavelength_input (if defined), or automatically computed with the smallest lambda/pixel 
;                       step obtained  from exp(loglam_gal). It is set to exp(LOGLAM\_TEMPLATES) if the keyword /rest\_ frame\_log is set.
;
; best_fit_model        [N x QQ flt array] It will contain the best fit models for each of the input galaxy spectra (dereddended if 
;                       MW_extinction is not zero), sampled over
;                       wavelength_output. The spectrum is in rest frame if  /rest_frame_log is set.
;
; galaxy_minus_ems_fit_model [N x QQ flt array] It will contain the input galaxy spectra minus the emission lines best fit models 
;                            (dereddended if MW_extinction is not zero),  sampled over wavelength_output, for each of the N
;                            input spectra.  The spectrum is in rest frame if  /rest_frame_log is set.
;
; best_template [N x QQ flt array] It will contain the best fitting template for each of the N input galaxy spectra
;                sampled over wavelength_output (rest frame wavelength). 
;
; best_template_LOSVD_conv [N x QQ flt array] It will contain the best fitting template for each of the N input galaxy spectra
;                convolved by best fitting LOSVD and sampled over wavelength_output (rest frame wavelength).
;
; reddening_output [N x 2 elements array]. Best fit values for the reddening of stars (reddening_output[*,0]) and gas (reddening_output[*,1]) 
;                                          for all the N galaxy spectra. 
;                                          If the reddening fit is not requred, the output value is set to 0. (reddening_output[*,0:1] = [0,0]). 
;                                          If only the reddening of the stars is fitted, the reddening of the gas is set to 0 (reddening_output[*,1] = 0).
;
; reddening_output_err &[N x 2 elements array]. Errors associated to reddening\_ output. If not fitted, the error is automatically set to 99.
; 
; residuals [N x QQ flt array] It contains the difference between the observed galaxy spectra (dereddened if the MW_reddening is defined) 
;            and thebest_fit_model, sampled over wavelength_output. The spectrum is in rest frame if  /rest_frame_log is set.
;
;
;

;*** OPTIONAL OUTPUTS ***
;
; version  [string]            Module version. If requested, the module is not execute and only version flag is returned
;
version_module = '0.7'
if n_elements(version) ne 0 then begin
 version = version_module
 goto, end_module
endif


c=299792.458d
;-- set_starting guesses and initializing extra variables
sz=size(galaxy)
if sz[0] eq 1 then galaxy = reform(galaxy,sz[1],1)
sz=size(galaxy)

if n_elements(star_kin_starting_guesses) eq 0 then begin
   star_kin_starting_guesses = fltarr(sz[2],2) ;start guess for velocity,sigma = 0
   star_kin_starting_guesses[*,1]=50. ;set start guess for sigma
endif
if n_elements(gas_kin_starting_guesses) eq 0 then star_kin_starting_guesses = star_kin_starting_guesses

if n_elements(extra_inputs) ne 0 then begin
   for k = 0, n_elements(extra_inputs)-1 do begin
       d = execute(extra_inputs[k])
       if not keyword_set(quiet) then  print, 'setting variable '+extra_inputs[k]
   endfor
endif

log_0_gal=min(loglam_gal)   ;log-wavelength start (for galaxy spectra)
log_step_gal=loglam_gal[1]-loglam_gal[0] ;log-wavelength step (for galaxy spectra)



;--

;-- reading emission line file WARNING: this section needs adjustments
;   for different sofware requirements (ppxf, gandalf, s-gandalf, etc)
readcol,emission_line_file,junk,em_name,wav,junk,junk,emss,format='I,A,F,A,A,F',comment='#',/silent
emss = mdap_sgn(emss)


INT_DISP = fltarr(n_elements(emss))
if n_elements(fwhm_instr_kmsec_matrix) eq 0 then INT_DISP =INT_DISP+52.
;--

;-- initialization of output variables

sztempl=size(templates)
if sztempl[0] eq 1 then templates = reform(templates,sztempl[1],1)
sztempl=size(templates)

;if not keyword_set(fix_star_kin) then begin
   stellar_kinematics=fltarr(sz[2],5)
   stellar_kinematics_err=fltarr(sz[2],4)
   stellar_weights=dblarr(sz[2],sztempl[2])
;endif
emission_line_kinematics=fltarr(sz[2],2)
emission_line_kinematics_err=fltarr(sz[2],2)
emission_line_intens=fltarr(sz[2],n_elements(wav))
emission_line_intens_err=fltarr(sz[2],n_elements(wav))
emission_line_fluxes=fltarr(sz[2],n_elements(wav))
emission_line_fluxes_err=fltarr(sz[2],n_elements(wav))
emission_line_equivW=fltarr(sz[2],n_elements(wav))
emission_line_equivW_err=fltarr(sz[2],n_elements(wav))
emission_line_kinematics_individual=fltarr(sz[2],n_elements(wav),2)
emission_line_kinematics_individual_err=fltarr(sz[2],n_elements(wav),2)

if keyword_set(rest_frame_log) then begin
      wavelength_output=exp(loglam_templates);exp(loglam_gal-sol[0]/velscale*(log_step_gal))
     ; if n_elements(wavelength_input) ne 0 then wavelength_output=wavelength_input
endif else begin

   if n_elements(wavelength_input) ne 0 then begin
      wavelength_output=wavelength_input
   endif else begin
                                ;tmp = double(10.^(loglam_gal))
      tmp = double(exp(loglam_gal))
      ;stp = double(tmp[1:n_elements(tmp)-1]-tmp[0:n_elements(tmp)-2])
      stp = tmp[1]-tmp[0];double(tmp[1:n_elements(tmp)-1]-tmp[0:n_elements(tmp)-2])
      n= (max(tmp) - min(tmp) ) / stp
      wavelength_output = fltarr(round(n))*stp + min(tmp)
   endelse
endelse
;w_step = (wavelength_output[1]-wavelength_output[0])

nw = n_elements(wavelength_output)
best_fit_model=fltarr(sz[2],nw, /NOZERO) 
best_template=fltarr(sz[2],nw, /NOZERO) 
best_template_LOSVD_conv=fltarr(sz[2],nw, /NOZERO) 
galaxy_minus_ems_fit_model=fltarr(sz[2],nw, /NOZERO) 
residuals=fltarr(sz[2],nw, /NOZERO) 
;best_fit_model=dblarr(sz[2],nw) 
;best_template=dblarr(sz[2],nw) 
;best_template_LOSVD_conv=dblarr(sz[2],nw) 
;galaxy_minus_ems_fit_model=dblarr(sz[2],nw) 
;residuals=dblarr(sz[2],nw) 


; do I really need it? best_fit_model_log=dblarr(sz[2],n_elements(loglam_gal))
; do I really need it? emission_model_log=dblarr(sz[2],n_elements(loglam_gal))
; do I really need it? optimal_template_log=dblarr(n_elements(loglam_templates),sz[2])
;-- 

;-- fitting loop

reddening_output = fltarr(sz[2],2)
reddening_output_err = fltarr(sz[2],2)
;dv = (loglam_templates[0]-loglam_gal[0])*c*alog(10.)
dv = (loglam_templates[0]-loglam_gal[0])*c

if not keyword_set(quiet) then print, 'Fitting '+mdap_stc(sz[2],/integer)+' spectra, please, wait...'
print, 'Fitting '+mdap_stc(sz[2],/integer)+' spectra, please, wait...'
;  'Spectrum ',mdap_stc(i+1,/integer),'/',mdap_stc(sz[2],/integer),' fitted'
;window,0,retain=2

if n_elements(mdegree) ne 0 then MDEGREE_=MDEGREE
if n_elements(reddening) ne 0 then junk = temporary(MDEGREE_)
;emission_line_file_ori = emission_line_file

FOR i = 0, sz[2]-1 DO BEGIN  ;loop over all the spectra
;FOR i = 45, 45 DO BEGIN  ;loop over all the spectra

   if n_elements(reddening) ne 0 then ebv=reddening
   start = [star_kin_starting_guesses[i,0],star_kin_starting_guesses[i,1],star_kin_starting_guesses[i,2],star_kin_starting_guesses[i,3],gas_kin_starting_guesses[i,0],gas_kin_starting_guesses[i,1]]

   if i gt 0 and keyword_set(use_previos_guesses) then begin
      start = [sol[0],sol[1],sol[2],sol[3],sol[7],sol[8]]
   endif

   galaxy_=galaxy[*,i]
   noise_=noise[*,i]


 ;  ;-- correct for galactic reddening, if provided
 ;  IF KEYWORD_SET(MW_extinction) THEN BEGIN
 ;     ;dereddening_attenuation = DUST_CALZETTI(l0_gal,lstep_gal,n_elements(galaxy_),-MW_extinction,0.0d,/log10)
 ;     dereddening_attenuation = mdap_dust_calzetti(log_0_gal,log_step_gal,n_elements(galaxy_),-MW_extinction,0.0d)
 ;     galaxy_ = galaxy_*temporary(dereddening_attenuation)
 ;  ENDIF

  ;-- fitting with mdap_sgandalf
; Print, ' Fitting spectrum ',mdap_stc(i+1,/integer),'/',mdap_stc(sz[2],/integer)   ;TEST LINE


;GOODPIXELS=where(galaxy_/noise_ ge .5 )
   if ~keyword_set(quiet) then print, 'fitting spectrum',i
;print, 'fitting spectrum',i
;MDEGREE_=MDEGREE
;if n_elements(reddening) ne 0 then junk = temporary(MDEGREE_)
;   mdap_sgandalf, templates,[[alog(wav)],[emss]], loglam_gal, galaxy_, noise_, velScale, start, sol, $
;       gas_intens,gas_fluxes,gas_ew,gas_intens_err,gas_fluxes_err,gas_ew_err,$
;       BESTFIT=bestFit, BIAS=bias,  MDEGREE=MDEGREE_,DEGREE=degree, ERROR=error, $
;       MOMENTS=moments, LAMBDA=LAMBDA,reddening=ebv,$
;       POLYWEIGHTS=polyweights, VSYST=dv, WEIGHTS=weights, BF_COMP2 = bf_comp2, BF_COMP1 = bf_comp1,$
;       OPT_TEMPL=ot1,mpoly=mpoly,ADDITIVE_POL=additive_pol,/quiet,$
;       fix_star_kin=fix_star_kin,fix_gas_kin=fix_gas_kin,$
;      range_v_star=range_v_star,range_s_star=range_s_star,range_v_gas=range_v_gas,range_s_gas=range_s_gas
;

    
   INT_DISP=interpol(fwhm_instr_kmsec_matrix[1,*],alog(fwhm_instr_kmsec_matrix[0,*]),alog(wav*(1.+start[2]/c)))/2.3548
status = 1
   mdap_gandalf_wrap,templates,loglam_templates,galaxy_,loglam_gal,noise_,velScale, start, sol,$
       sol_gas_V,esol_gas_V,sol_gas_S,esol_gas_S, $
       EMISSION_SETUP_FILE=emission_line_file, $
       gas_intens,gas_fluxes,gas_ew,gas_intens_err,gas_fluxes_err,gas_ew_err,$
       BESTFIT=bestFit, BIAS=bias,  MDEGREE=MDEGREE,DEGREE=degree, ERROR=error, $
       MOMENTS=moments, reddening=ebv,err_reddening=err_reddening,$
       VSYST=dv, WEIGHTS=weights, BF_COMP2 = bf_comp2,$
       FOR_ERRORS=1,$    ; ERROR COMPUTATION MUST BE TRIGGERED!!!
       fix_star_kin=fix_star_kin,fix_gas_kin=fix_gas_kin,$
       range_v_star=range_v_star,range_s_star=range_s_star,range_v_gas=range_v_gas,range_s_gas=range_s_gas,quiet=quiet,$
       mask_range=mask_range,fitted_pixels=fitted_pixels,external_library=external_library,INT_DISP=INT_DISP,status = status
   
   if status eq 0 then begin  ; if the fit failed, I set dummy outputs to not to stop the workflow.
     sol = fltarr(9)/0.
     error = fltarr(8)+99.
     gas_intens=fltarr(n_elements(wav))
     gas_fluxes=fltarr(n_elements(wav))
     gas_ew=fltarr(n_elements(wav))
     gas_intens_err=fltarr(n_elements(wav))+99.
     gas_fluxes_err=fltarr(n_elements(wav))+99.
     gas_ew_err=fltarr(n_elements(wav))+99.
     bestFit=loglam_gal*0.
     if n_elements(reddening) ne 0 then ebv= ebv*0.
     if n_elements(reddening) ne 0 then err_reddening= ebv*0.+99

     weights = dblarr(sz[2],sztempl[2])
     bf_comp2=loglam_gal*0.
   endif


;stop

; If the error vector is flat (i.e. errors are not reliable), I rescale the formal errors for sqrt(chi2/dof), as instructed by mpfit and ppxf.
if min(noise_) eq max(noise_) then error[0:5] = error[0:5] * sqrt(sol[6]) ; If the error vector is flat (i.e. errors are not reliable), I rescale the formal errors for sqrt(chi2/dof), as instructed by mpfit and ppxf.  ma005_142.790030+22.746507datacubes_block5.idl

   ; plots for checks... remove these lines when running on remote server
    ;plot,exp(loglam_gal), galaxy_,title='GANDALF + '+string(i),xrange=[3600,10350],xstyle=1
    ;oplot,exp(loglam_gal),bestfit,color=200
    ;if n_elements(fitted_pixels) ne 0 then OPLOT,exp(loglam_gal[fitted_pixels]),bestfit[fitted_pixels],color=200
    
    ;print,'start', start
    ;print,'sol',  sol
    ;print,'error', error
   ;if n_elements(ebv) ne 0 then  print, 'reddening = ',ebv
   ;if n_elements(ebv) eq 0 then  print, 'reddening = Not FITTED'
   ; print, ''
   ; print, minmax(exp(loglam_gal)).reset

   ; do I really need it? best_fit_model_log[i,*]=bestfit
   ; do I really need it? emission_model_log[i,*]=bf_comp2
  ; print, 'spectrum '+string(i)
;top

   ;-- storing outputs
   ;reddening_output[i] = (n_elements(ebv) eq 0) ?  0 : ebv[0]; reddening_output[i] = 0 : reddening_output[i] = ebv[0]
   if n_elements(ebv) eq 2 then reddening_output[i,*]= ebv
   if n_elements(ebv) eq 1 then reddening_output[i,*]= [ebv[0],0.]
   if n_elements(ebv) eq 2 then reddening_output_err[i,*]= err_reddening
   if n_elements(ebv) eq 1 then reddening_output_err[i,*]= [err_reddening[0],999.]
   if n_elements(ebv) eq 0 then reddening_output[i,*]= [0,0]
   if n_elements(ebv) eq 0 then reddening_output_err[i,*]= [99,99]

   bf_template = (templates # weights[0:sztempl[2]-1])
  ; do I really need it?  optimal_template_log[*,i]=bf_template
   mdap_get_losvd,sol,velscale,losvd
   bf_template_LOSVD = convol(bf_template,losvd)

   ;if not keyword_set(fix_star_kin) then begin
      stellar_kinematics[i,*]=[sol[0:3],sol[6]]  ; store V,sigma,h3, h4, (stars) and chi2/DOF; not storing h5 and h6.
      stellar_kinematics_err[i,*]=error[0:3]  ; store errors for V,sigma,h3, and h4 (stars); not storing h5 and h6.
      stellar_weights[i,*]= weights[0:sztempl[2]-1]
   ;endif
   emission_line_kinematics[i,*]=sol[7:8]
   emission_line_kinematics_err[i,*]=error[6:7]
;   stop
   emission_line_intens[i,*]= gas_intens; corrected for reddening, if fitted
   emission_line_intens_err[i,*]= gas_intens_err
   emission_line_equivW[i,*]= gas_ew
   emission_line_equivW_err[i,*]= gas_ew_err
   emission_line_fluxes[i,*]= gas_fluxes     ; corrected for reddening, if fitted
   emission_line_fluxes_err[i,*] = gas_fluxes_err

   emission_line_kinematics_individual[i,*,0]=sol_gas_V
   emission_line_kinematics_individual_err[i,*,0]=esol_gas_V
   emission_line_kinematics_individual[i,*,1]=sol_gas_S
   emission_line_kinematics_individual_err[i,*,1]=esol_gas_S
   ;-- check how many emission lines have intensity > 4*intensity_error
   ;      (i.e. 3 sigma detections)
   ;test = abs(emission_line_intens[i,*]) - 4.*abs(emission_line_intens_err[i,*])
   test = abs(gas_intens) - abs(gas_intens_err)
   indici = where(test gt 0 and gas_intens gt 0) 
;stop
  ; stop
;    window,0,retain=2,xsize=1900
;      residuals_ = galaxy_-bestfit
;     plot,exp(loglam_gal), bf_comp2,title=string(i),xrange=[3500,7500],yrange=[0.0,max(bf_comp2)],xstyle=1,ystyle=1
;    oplot,exp(loglam_gal),residuals_,color=200,psym=3
;
;    for kk = 0, n_elements(indici)-1 do oplot,[wav[indici[kk]]*(1+sol[0]/c),wav[indici[kk]]*(1+sol[0]/c)],[0,1],color=230,thick=2,linestyle=1;

; xyouts,wav[indici[kk]]*(1+sol[0]/c),gas_intens[indici[kk]],stc(abs(gas_intens[indici[kk]]/gas_intens_err[indici[kk]])),orientatio=90

   ; oplot, exp(loglam_gal),bf_comp2,color=200
   ; for kk = 0, n_elements(gas_intens)-1 do print,gas_intens[kk],gas_intens_err[kk] , abs(gas_intens[kk]/gas_intens_err[kk]),wav[kk]*(1+sol[0]/c)
   ; for kk = 0, n_elements(gas_intens)-1 do xyouts,wav[kk]*(1+sol[0]/c),gas_intens[kk],stc(abs(gas_intens[kk]/gas_intens_err[kk])),orientatio=90




;stop

    if n_elements(indici) lt 3 then begin  ;I want at least 3 emission lines (with positive intensity) to be detected at 1 sigma level....
       ;emission_line_intens[i,*] = 0.
       ;emission_line_intens_err[i,*] = 0.
       ;emission_line_fluxes[i,*] = 0./0.
       ;emission_line_fluxes_err[i,*] = 0./0.
       ;emission_line_equivW[i,*] = 0./0.
       ;emission_line_equivW_err[i,*] = 0./0.
       emission_line_kinematics[i,*] = 0./0.
       emission_line_kinematics_err[i,*] = 0./0.
    endif             
   ;--

;stop
   if keyword_set(rest_frame_log) then begin
      rf_gal_lam = exp(loglam_gal-sol[0]/velscale*(log_step_gal))
      bestFit_interp = interpol(bestFit,rf_gal_lam,wavelength_output)
      best_fit_model[i,*] = temporary(bestFit_interp)

      galaxy_minus_ems_fit_model_interp = interpol(galaxy_-bf_comp2,rf_gal_lam,wavelength_output);
      galaxy_minus_ems_fit_model[i,*] = temporary(galaxy_minus_ems_fit_model_interp)

      residuals_interp = interpol(galaxy_-bestfit,rf_gal_lam,wavelength_output)
      residuals[i,*] = temporary(residuals_interp)

      ;bf_template = interpol(bf_template,loglam_templates,wavelength_output)
      best_template[i,*] = temporary(bf_template)

      ;bf_template_LOSVD = interpol(bf_template_LOSVD,loglam_templates,wavelength_output)
      best_template_LOSVD_conv[i,*] =temporary(bf_template_LOSVD)

   endif else begin
      ;bestFit_linear = interpol(bestFit,10.^(loglam_gal),wavelength_output);/(wavelength_output[1]-wavelength_output[0])
      bestFit_linear = float(interpol(bestFit,exp(loglam_gal),wavelength_output)) ;/w_step
      best_fit_model[i,*] = temporary(bestFit_linear)
      galaxy_minus_ems_fit_model_log=galaxy_-bf_comp2
     ;galaxy_minus_ems_fit_model_linear = interpol(galaxy_minus_ems_fit_model_log,10.^(loglam_gal),wavelength_output);/(wavelength_output[1]-wavelength_output[0])
      galaxy_minus_ems_fit_model_linear = float(interpol(galaxy_minus_ems_fit_model_log,exp(loglam_gal),wavelength_output)) ;/w_step
      galaxy_minus_ems_fit_model[i,*] = temporary(galaxy_minus_ems_fit_model_linear)

      residuals_ = galaxy_-bestfit
      ;residuals_linear = interpol(residuals_,10.^(loglam_gal),wavelength_output);/(wavelength_output[1]-wavelength_output[0])
      residuals_linear = float(interpol(residuals_,exp(loglam_gal),wavelength_output)) ;/w_step
      residuals[i,*] = temporary(residuals_linear)

      ;bf_template_linear = interpol(bf_template,10.^(loglam_templates),wavelength_output);/(wavelength_output[1]-wavelength_output[0])
      bf_template_linear = float(interpol(bf_template,exp(loglam_templates),wavelength_output)) ;/w_step
      best_template[i,*] = temporary(bf_template_linear)

      ;bf_template_LOSVD_linear = interpol(bf_template_LOSVD,10.^(loglam_templates),wavelength_output);/(wavelength_output[1]-wavelength_output[0])
      bf_template_LOSVD_linear = float(interpol(bf_template_LOSVD,exp(loglam_templates),wavelength_output)) ;/w_step
      best_template_LOSVD_conv[i,*] =temporary(bf_template_LOSVD_linear)
   endelse

;  endelse

;stop

ENDFOR
;--
;junk = temporary(reddening)

;stop
end_module:
end




  ;-- 
;   ;-- storing outputs
;   ;reddening_output[i] = (n_elements(ebv) eq 0) ?  0 : ebv[0]; reddening_output[i] = 0 : reddening_output[i] = ebv[0]
;   if n_elements(ebv) ne 0 then reddening_output[i] = ebv[0]
;   if n_elements(MW_extinction) ne 0 then reddening_output[i] = reddeningn_output[i] - MW_extinction[0] 
;   bf_template = (templates # weights[0:sztempl[2]-1])
;   optimal_template_log[*,i]=bf_template
;   mdap_get_losvd,sol,velscale,losvd
;   bf_template_LOSVD = convol(bf_template,losvd)
;   stellar_kinematics[i,*]=sol[0:6]
;   stellar_weights[i,*]= weights[0:sztempl[2]-1]
;   emission_line_kinematics[i,*]=sol[7:8]
;   emission_line_fluxes[i,*]= weights[sztempl[2]:*]*emss;
;
;   ;bestFit_linear = mdap_rebin_spectrum(bestFit,10.^(loglam_gal),wavelength_output)/(wavelength_output[1]-wavelength_output[0])
;   bestFit_linear = mdap_rebin_spectrum(bestFit,exp(loglam_gal),wavelength_output)/w_step
;   best_fit_model[i,*] = temporary(bestFit_linear)
;
;   ;bf_template_linear = mdap_rebin_spectrum(bf_template,10.^(loglam_templates),wavelength_output)/(wavelength_output[1]-wavelength_output[0])
;   bf_template_linear = mdap_rebin_spectrum(bf_template,exp(loglam_templates),wavelength_output)/w_step
;   best_template[i,*] = temporary(bf_template_linear)
;
;   ;bf_template_LOSVD_linear = mdap_rebin_spectrum(bf_template_LOSVD,10.^(loglam_templates),wavelength_output)/(wavelength_output[1]-wavelength_output[0])
;   bf_template_LOSVD_linear = mdap_rebin_spectrum(bf_template_LOSVD,exp(loglam_templates),wavelength_output)/w_step
;   best_template_LOSVD_conv[i,*] = temporary(bf_template_LOSVD_linear)
;
;   galaxy_minus_ems_fit_model_log=galaxy_-bf_comp2
;   ;galaxy_minus_ems_fit_model_linear = mdap_rebin_spectrum(galaxy_minus_ems_fit_model_log,10.^(loglam_gal),wavelength_output)/(wavelength_output[1]-wavelength_output[0])
;   galaxy_minus_ems_fit_model_linear = mdap_rebin_spectrum(galaxy_minus_ems_fit_model_log,exp(loglam_gal),wavelength_output)/w_step
;   galaxy_minus_ems_fit_model[i,*] = temporary(galaxy_minus_ems_fit_model_linear)
;
;   residuals_ = galaxy_-bestfit
;   ;residuals_linear = mdap_rebin_spectrum(residuals_,10.^(loglam_gal),wavelength_output)/(wavelength_output[1]-wavelength_output[0])
;   residuals_linear = mdap_rebin_spectrum(residuals_,exp(loglam_gal),wavelength_output)/w_step
;   residuals[i,*] = temporary(residuals_linear)

  ;-- 



;  if keyword_set(flux) then begin
;     ;bestFit_linear = mdap_rebin_spectrum(bestFit,10.^(loglam_gal),wavelength_output)
;     bestFit_linear = mdap_rebin_spectrum(bestFit,exp(loglam_gal),wavelength_output)
;  endif else begin
;     ;bestFit_linear = interpol(bestFit,10.^(loglam_gal),wavelength_output)
;     bestFit_linear = interpol(bestFit,exp(loglam_gal),wavelength_output)
;  endelse
;  best_fit_model[i,*] = temporary(bestFit_linear)
;;
;
;  if keyword_set(flux) then begin
;     ;bf_template_linear = mdap_rebin_spectrum(bf_template,10.^(loglam_templates),wavelength_output)/(wavelength_output[1]-wavelength_output[0])
;     bf_template_linear = mdap_rebin_spectrum(bf_template,exp(loglam_templates),wavelength_output)/(wavelength_output[1]-wavelength_output[0])
;  endif else begin
;     ;bf_template_linear = interpol(bf_template,10.^(loglam_templates),wavelength_output)
;     bf_template_linear = interpol(bf_template,exp(loglam_templates),wavelength_output)
;  endelse
;  best_template[i,*] = temporary(bf_template_linear)
;
;  if keyword_set(flux) then begin
;     ;bf_template_LOSVD_linear = mdap_rebin_spectrum(bf_template_LOSVD,10.^(loglam_templates),wavelength_output)/(wavelength_output[1]-wavelength_output[0])
;     bf_template_LOSVD_linear = mdap_rebin_spectrum(bf_template_LOSVD,exp(loglam_templates),wavelength_output)/(wavelength_output[1]-wavelength_output[0])
;  endif else begin
;     ;bf_template_LOSVD_linear =interpol(bf_template_LOSVD,10.^(loglam_templates),wavelength_output)
;     bf_template_LOSVD_linear =interpol(bf_template_LOSVD,exp(loglam_templates),wavelength_output)
;  endelse
;  best_template_LOSVD_conv[i,*] = temporary(bf_template_LOSVD_linear)
;
;  galaxy_minus_ems_fit_model_log=galaxy_-bf_comp2
;  if keyword_set(flux) then begin
;     ;galaxy_minus_ems_fit_model_linear = mdap_rebin_spectrum(galaxy_minus_ems_fit_model_log,10.^(loglam_gal),wavelength_output)/(wavelength_output[1]-wavelength_output[0])
;     galaxy_minus_ems_fit_model_linear = mdap_rebin_spectrum(galaxy_minus_ems_fit_model_log,exp(loglam_gal),wavelength_output)/(wavelength_output[1]-wavelength_output[0])
;  endif else begin
;     ;galaxy_minus_ems_fit_model_linear = interpol(galaxy_minus_ems_fit_model_log,10.^(loglam_gal),wavelength_output)
;     galaxy_minus_ems_fit_model_linear = interpol(galaxy_minus_ems_fit_model_log,exp(loglam_gal),wavelength_output)
;  endelse
;  galaxy_minus_ems_fit_model[i,*] = temporary(galaxy_minus_ems_fit_model_linear)
;
;  residuals_ = galaxy_-bestfit
;  if keyword_set(flux) then begin
;     ;residuals_linear = mdap_rebin_spectrum(residuals_,10.^(loglam_gal),wavelength_output)/(wavelength_output[1]-wavelength_output[0])
;     residuals_linear = mdap_rebin_spectrum(residuals_,exp(loglam_gal),wavelength_output)/(wavelength_output[1]-wavelength_output[0])/(wavelength_output[1]-wavelength_output[0])
;  endif else begin
;     ;residuals_linear = interpol(residuals_,10.^(loglam_gal),wavelength_output)
;     residuals_linear = interpol(residuals_,exp(loglam_gal),wavelength_output)
;  endelse
;  residuals[i,*] = temporary(residuals_linear)
;
;   ;-- 


;stop
