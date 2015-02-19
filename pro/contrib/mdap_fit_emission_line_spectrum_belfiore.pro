;+
; NAME:
;       MDAP_FIT_EMISSION_LINE_SPECTRUM_BELFIORE
;
; PURPOSE:
;
; CALLING SEQUENCE:
;       MDAP_FIT_EMISSION_LINE_SPECTRUM_BELFIORE, wave, galaxy_eml_only, flux_err, mask, $
;                                                 star_sigma, redshift, El, eml_model, /quiet
;
; INPUTS:
;       wave dblarr[C]
;           Wavelength of each of the C spectral channels.
;
;       galaxy_eml_only dblarr[C]
;           Flux of the galaxy spectrum with the stellar continuum
;           subtracted in each of the C spectral channels.
;
;       flux_err dblarr[C]
;           Error in the continuum-subtracted galaxy spectrum.
;
;       mask dblarr[C]
;           Bad pixel mask for (1-bad;0-good) the continuum-subtrated
;           galaxy spectrum.
;
;       star_sigma double
;           Velocity dispersion of the stars found for this fiber.  Used
;           as a guess value when fitting the gas dispersion.
;
;       redshift double
;           Systemic redshift (z).  All velocities (cz) are measured
;           with respect to this redshift.
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;       /quiet
;           Suppress output written to stdout.
;
; OUTPUT:
;       El struct
;           Structure contiaining the results of the line fits.  Defined as:
;            
;           El = { name:strarr(nlines2), lambda:fltarr(nlines2), Ampl:fltarr(nlines2), $
;                  Vel:fltarr(nlines2), Sigma:fltarr(nlines2), eAmpl:fltarr(nlines2), $
;                  eVel:fltarr(nlines2), eSigma:fltarr(nlines2) }
;
;           where nlines2 is the number of fitted lines.  For now this
;           is hard-coded to be the 9 lines defined by
;           MDAP_DEFINE_EMISSION_LINES_ENCI_BELFIORE.
;
;       eml_model dblarr[C]
;           The best-fitting model emission-line spectrum.
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
;   --From Francesco Belfiore's fb_mdap_spectral_fitting_v3.pro
;       -- Mar 2013: (FB) Original implementation. Adapted from fb_mdap
;                         used for p_manga data
;       21 Oct 2014: (FB) Added documentation. Slightly modified the
;                         output structures to allow saving of the
;                         coefficients of additive and multiplicative
;                         polynomials
;       19 Nov 2014: (FB) fixed a bug with the peakfitting of tied lines
;       27 Nov 2014: (FB) better estimate of the starting guess for the
;                         velocity of the lines
;       09 Dec 2014: (FB) new first guess for sigma of emission lines is
;                         80 km s-1. Uncommented the tied lines option
;   --
;   11 Dec 2014: (KBW) Edited down to basic emission-line-only fitting
;                      aspects for inclusion in the DAP.
;   19 Feb 2015: (KBW) Fixed an error in where counting (line 184)
;-
;------------------------------------------------------------------------------


;----------------------------------------------------------------------------
; Define the lines to fit
PRO MDAP_DEFINE_EMISSION_LINES_BELFIORE, $
                n_l, w_l, fit_l, con_l

        n_l = [ 'OII_3727', 'Hb', 'OIII_4959', 'OIII_5007', 'NII_6548', 'Ha', 'NII_6583', $
                'SII_6716', 'SII_6731' ]
        ;w_l = [ 3727.38, 4861.33, 4958.92, 5006.84, 6548.03, 6562.80, 6583.41, 6716.47, 6730.85 ]
        fit_l = [ 1, 1, 1, 1, 1, 1, 1, 1, 1 ]
        con_l = [ 0, 0, 2, 2, 1, 0, 1, 0, 0 ]

        ; Force the rest wavelengths to be the same between Enci and Belfiore fits
        eml_par = MDAP_DEFINE_EMISSION_LINES_ENCI_BELFIORE()
        w_l = eml_par.lambda

END
;----------------------------------------------------------------------------


;----------------------------------------------------------------------------
;function defining one or more Gaussians, input for MPFIT
function gausfit, p, x=x, y=y, err=err, nlines2=nlines2, wl=wl, redshift=redshift

  c=299792.458d  ; speed of light in KM/s
  g_tot=0.0
  for i=0, nlines2-1 do begin
  g_p= p[3*i+2]*exp( -( (  x - wl[i]*(1+p[3*i]/c+redshift) )^2 )/ (2.*p[3*i+1]^2)  )
  g_tot=g_tot+g_p
  endfor

  return, ((y-g_tot)/err)
end
;----------------------------------------------------------------------------


;----------------------------------------------------------------------------
; procedure performing the actual line fitting using MPFIT
PRO peak_fitting, resid, noise_, wave, star_sigma, redshift, w_l, n_l, nlines2, l, eflg=eflg, $
                  quiet=quiet
  c=299792.458d
  
  ; **** Set Parameters for MPFIT
  err=noise_
  ;x = exp(loglam_gal)
  x = wave
  y = resid
  
  ; *** Declare Emission Lines
  
  L=REPLICATE({n:' ', wl:0.0, s:0.0, a:0.0, v:0.0, es:0.0, ea:0.0, ev:0.0}, nlines2)
  l.n=N_L
  l.wl=W_L
  
  ; *** Set emission lines starting parameters
  ;
  ;emission lines sigma and velocity are set equal to those of the stars
  sigma2=fltarr(nlines2)
  v2=fltarr(nlines2)
  ;NOTE sometimes the stars are badly fitted, especially when the continuum is very low. In this case the velocity/sigma for the stars can  be nonsense.
  for i=0, nlines2-1 do begin
;   if star_sigma le 200 and star_sigma ge 70 then begin
    if star_sigma le 200 && star_sigma ge 70 then begin
      sigma2[i]=l[i].wl*star_sigma/c
    endif else begin
      sigma2[i]=l[i].wl*80/c
    endelse
  endfor
  
  mask=600.0                                                                          ;draw a mask of 600 km s-1 around each line
  dw=mask*l.wl/c                                                                      ;width of the mask in angstrom
 
                                                      ;guess for the line velocity using the max of the flux within the mask
  
  for i=0, nlines2-1 do begin
      ;subscripts of the wav range within the mask
      wrange=where(x gt l[i].wl*(1+redshift)-dw[i] and x lt l[i].wl*(1+redshift)+dw[i], count)
      if count eq 0 then begin
        amp = 0.0
        vguess = redshift*c
        continue
      endif
      amp=max(resid[wrange], w_max)                                                       ;max of the flux within the mask, and its subscripts (w_max)
      x_max=x[wrange[w_max]]                                                                      ;wavelength of max subscript
      vguess=(x_max - l[i].wl*(1+redshift))/(l[i].wl*(1+redshift))*c
;      if vguess gt -400 and vguess lt 400 then begin
       if vguess gt -400 && vguess lt 400 then begin
          v2[i]=vguess
       endif else begin
          v2[i]=0
       endelse

  endfor

  p_start=fltarr(3*nlines2)
  for i=0, nlines2-1 do begin
    ;subscripts of the wav range withing the mask
    wrange=where(x gt l[i].wl*(1+redshift)-dw[i] and x lt l[i].wl*(1+redshift)+dw[i], count)
    if count eq 0 then begin
      amp = 0.0
      vguess = redshift*c
      continue
    endif
    amp=max(resid[wrange], w_max)                                                       ;max of the flux within the mask, and its subscripts (w_max)
    p_start[3*i]   =  v2[i]
    p_start[1+3*i] =  sigma2[i]

    if amp gt 0 then begin
      p_start[2+3*i] =  amp
    endif else begin
      p_start[2+3*i] = 0
     endelse 
    endfor
    

    ;print, p_start
    ; *** Constrain parameters
    
    n_parameters = 3*nlines2
    parinfo=replicate({fixed:0,limited:[0,0],limits:[0.,0.], tied:''},n_parameters)
    
    ;all amplitudes need to be positive!
    for i=0, nlines2-1 do begin
      parinfo[3*i+2].limited=[1,0]
      parinfo[3*i+2].limits=[0]
    endfor
    
    
    ;Tie the relative amplitude and sigmas and velocities of multiplet lines (here only [OIII] and [NII])
    wOIII_4959=where( l.n eq 'OIII_4959', count_4959)
    wOIII_5007=where(l.n eq 'OIII_5007', count_5007)
    
    
    WNII_6548=where(l.n eq 'NII_6548', count_6548)
    wNII_6583=where(l.n eq 'NII_6583', count_6583)
    ;print, 'OIII 4959', WOIII_4959, 'OIII 5007', WOIII_5007

;   if wOIII_4959[0] ne -1 and wOIII_5007[0] ne -1 then begin
;   if count_4959 ne 0 and count_5007 ne 0 then begin
    if count_4959 ne 0 && count_5007 ne 0 then begin
      ;ratio of amplitudes fixed to theoretical
      parinfo[3*wOIII_4959+2].tied='0.330*p['+string(3*wOIII_5007+2)+']'
      ;fix to have same velocities
      parinfo[3*wOIII_4959].tied='p['+string(3*wOIII_5007)+']'
      ;fix to have same sigma
      parinfo[3*wOIII_4959+1].tied='p['+string(3*wOIII_5007+1)+']'
    endif
    
;   if WNII_6548[0] ne -1 and wNII_6583[0] ne -1 then begin
;   if count_6548 ne 0 and count_6583 ne 0 then begin
    if count_6548 ne 0 && count_6583 ne 0 then begin
      ;ratio of amplitudes fixed to theoretical
      parinfo[3*WNII_6548+2].tied='0.340*p['+string(3*wNII_6583+2)+']'
      ;fix to have same velocities
      parinfo[3*WNII_6548].tied='p['+string(3*wNII_6583)+']'
      ;fix to have same sigma
      parinfo[3*WNII_6548+1].tied='p['+string(3*wNII_6583+1)+']'
    endif
    
    ; TODO: KBW add the SII lines to this list?

    for i=0, nlines2-1 do begin              ;all sigmas between 70 and 200 km sec-1
      parinfo[3*i+1].limited=[1,1]
      parinfo[3*i+1].limits=[70, 200]*l[i].wl/c
    endfor
    
    for i=0, nlines2-1 do begin              ;all velocities between -400 and 400 km sec-1
      parinfo[3*i].limited=[1,1]
      parinfo[3*i].limits=[-400, 400]
    endfor
    
    if nlines2 gt 1 then begin                  ;if there are more than one lines
      for i=0, nlines2-2 do begin              ;all velocities are tied
        parinfo[3*(i+1)].tied=strcompress('p[0]', /remove_all)
      endfor
    endif
    
    
    ; *** Fit function
    p = mpfit('gausfit',p_start, functargs={x:x,y:y,err:err, nlines2:nlines2, wl:l.wl, redshift:redshift}, bestnorm=chi2, parinfo=parinfo, perror=perror, status=status, errmsg=errmsg, quiet=1)

    ;g_tot=0.0
    eflg = 0
    if status LE 0 then begin
        print, errmsg
        eflg = 1
        return
    endif
    ;;print, p
    ;print, 'fitting lines', l.n
    ;for kk=0, nlines2-1 do begin
    ;  g_p= p[3*kk+2]*exp( -( (  x - l[kk].wl*(1+p[3*kk]/c+redshift) )^2 )/ (2.*p[3*kk+1]^2)  )
    ;  g_tot=g_tot+g_p
    ;endfor
    
    
    ; *** Calculate the reduced chi^2
    n_data = n_elements(lambda)
    n_pars = n_elements(p)
    chi2_reduced = chi2/(n_data-n_pars-1)
    
    
    ; *** Save Output
    
    for i=0, nlines2-1 do begin
      l[i].s=p[1+3*i]/l[i].wl*c
      l[i].es = perror[3*i+1]/l[i].wl*c
      l[i].v=  p[3*i]
      l[i].ev = perror[3*i]
      l[i].a=p[2+3*i]
      l[i].ea=perror[2+3*i]
      
    endfor

    if ~keyword_set(quiet) then $
        print, 'sigma, vel, ampl', l.s, l.v, l.a
  END
;----------------------------------------------------------------------------


;----------------------------------------------------------------------------
; Pulled out section of the code that fits a single emission-line-only
; spectrum
PRO MDAP_FIT_EMISSION_LINE_SPECTRUM_BELFIORE, $
                wave, galaxy_eml_only, flux_err, mask, star_sigma, redshift, El, eml_model, $
                quiet=quiet

    ; Define the lines to fit; !! HARD-CODED !!
    MDAP_DEFINE_EMISSION_LINES_BELFIORE, n_l, w_l, fit_l, con_l
    nlines2 = n_elements(n_l)                       ; Number of lines
  
    ; Output structure
    El = { name:strarr(nlines2), lambda:fltarr(nlines2), Ampl:fltarr(nlines2), $
           Vel:fltarr(nlines2), Sigma:fltarr(nlines2), eAmpl:fltarr(nlines2), $
           eVel:fltarr(nlines2), eSigma:fltarr(nlines2) }

    indx = where(mask lt 1., count)
;   if indx[0] eq -1 then $
    if count eq 0 then $
        message, 'Entire spectrum masked!'

    ;---------------------------------------------------------------
    ; Fit emission lines

    ; TODO: Unecessary with current HARD-CODED lines
    wnofit=where(fit_l gt 0, count)
    if count eq 0 then $
        message, 'No lines to fit!'
    w_l=w_l[wnofit]
    n_l=n_l[wnofit]
    con_l=con_l[wnofit]
    fit_l=fit_l[wnofit]
    
    ; Select untied lines
    untied=where(con_l eq 0, count)
    n_tied=max(con_l)
   
    ; *** Call the emission line fitting prcedure ***
    
    ;*** Fit the untied lines
    if count ne 0 then begin
      
    for kl=0, n_elements(untied)-1 do begin
      s_l=untied[kl]
      wl_m=w_l[untied[kl]]
      nl_m=n_l[untied[kl]]
      
      ; *** if the line is strong (fit_l = 1) then always fit it ***
    
        peak_fitting, galaxy_eml_only[indx], flux_err[indx], wave[indx], star_sigma, redshift, $
                      wl_m, nl_m, 1, l, eflg=eflg, quiet=quiet
        if eflg eq 1 then begin
            El.name[s_l]=l[0].n
            El.lambda[s_l]=l[0].wl
            El.ampl[s_l]=0.
            El.vel[s_l]=0.
            El.sigma[s_l]=0.
            El.eampl[s_l]=1.
            El.evel[s_l]=1.
            El.esigma[s_l]=1.
        endif else begin
            El.name[s_l]=l[0].n
            El.lambda[s_l]=l[0].wl
            El.ampl[s_l]=l[0].a
            El.vel[s_l]=l[0].v
            El.sigma[s_l]=l[0].s
            El.eampl[s_l]=l[0].ea
            El.evel[s_l]=l[0].ev
            El.esigma[s_l]=l[0].es
            if ~keyword_set(quiet) then $
                print, 'Succesfully fitted ', nl_m, ' untied lines'
        endelse
    endfor
    
    endif
    ;now I've fitted all the untied lines

    
    ;*** Fitting of the kinematically tied lines (only the velocity is tied!)***
    print, n_tied
    for kl=1, n_tied do begin
      wgroup=where(con_l eq kl, count)
;     while wgroup[0] eq -1 do begin
      while count eq 0 do begin
        if ~keyword_set(quiet) then $
            print, 'Error in your input emission line file!'
        kl=kl+1
        wgroup=where(con_l eq kl, count)
      endwhile
      
      wl_m=w_l[wgroup]
      nl_m=n_l[wgroup]
      n_m=N_elements(wgroup)
      
        peak_fitting, galaxy_eml_only[indx], flux_err[indx], wave[indx], star_sigma, redshift, $
                      wl_m, nl_m, n_m, l, eflg=eflg, quiet=quiet
        
        if eflg eq 1 then begin
            El.name[wgroup]=l.n
            El.lambda[wgroup]=l.wl
            El.ampl[wgroup]=0.
            El.vel[wgroup]=0.
            El.sigma[wgroup]=0.
            El.eampl[wgroup]=1.
            El.evel[wgroup]=1.
            El.esigma[wgroup]=1.
        endif else begin
            El.name[wgroup]=l.n
            El.lambda[wgroup]=l.wl
            El.ampl[wgroup]=l.a
            El.vel[wgroup]=l.v
            El.sigma[wgroup]=l.s
            El.eampl[wgroup]=l.ea
            El.evel[wgroup]=l.ev
            El.esigma[wgroup]=l.es
            if ~keyword_set(quiet) then $
                print, 'I fitted the ', nl_m, ' lines'
        endelse
    endfor
   
   ; Get the full model of the emission-line spectrum (ignoring masks)
   c=299792.458d  ; speed of light in KM/s
    eml_model=wave*0.0
   for kk=0, nlines2-1 do begin
     a1=El.ampl[kk]
     a1=a1[0]
     l1=El.lambda[kk]
     l1=l1[0]
     v1=El.vel[kk]
     v1=v1[0]
     s1=El.sigma[kk]
     s1=s1[0]*l1/c

     g_p= a1*exp( -( (  wave - l1*(1+v1/c+redshift) )^2 )/ (2.*s1^2)  )

     eml_model=eml_model+g_p
   endfor

END
;----------------------------------------------------------------------------





