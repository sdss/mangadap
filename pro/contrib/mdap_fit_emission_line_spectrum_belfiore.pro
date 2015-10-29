;+
; NAME:
;       MDAP_FIT_EMISSION_LINE_SPECTRUM_BELFIORE
;
; PURPOSE:
;
; CALLING SEQUENCE:
;       MDAP_FIT_EMISSION_LINE_SPECTRUM_BELFIORE, wave, galaxy_eml_only, flux_err, mask, $
;                                                 star_sigma, redshift, El, eml_model, /zero_base, $
;                                                 /quiet
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
;       /zero_base
;           Force the baseline to be zero
;        
;       /quiet
;           Suppress output written to stdout.
;
; OUTPUT:
;       El struct
;           Structure contiaining the results of the line fits.  Defined as:
;            
;           El = { name:strarr(nlines2), lambda:dblarr(nlines2), Ampl:dblarr(nlines2), $
;                  Vel:dblarr(nlines2), Sigma:dblarr(nlines2), eAmpl:dblarr(nlines2), $
;                  eVel:dblarr(nlines2), eSigma:dblarr(nlines2) }
;
;           where nlines2 is the number of fitted lines.  For now this
;           is hard-coded to be the lines defined by
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
; TODO:
;   - Handle when MPFITEXPR() status = 5 in some other way than saving
;     the result?
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
;       11 Dec 2014: Edited down to basic emission-line-only fitting
;                    aspects for inclusion in the DAP by K. Westfall
;                    (KBW).
;       19 Feb 2015: (KBW) Fixed an error in where counting (line 184)
;       13 Aug 2015: (KBW) Fit OII as a doublet, and add the fit to the
;                          OI doublet.  The velocities of the OI, OII,
;                          OIII, NII, and SII lines are now tied; the
;                          SII lines were NOT tied in the previous
;                          version.  As before, only the flux ratios of
;                          the OIII and NII lines are fixed.
;                          Documentation changes.  Added warning if
;                          MPFIT() reaches the maximum number of
;                          iterations.
;       16 Sep 2015: (KBW) Allow for a *local* non-zero baseline.
;                          Implementing this required moving away from
;                          Francesco's original implementation a bit.
;                          PEAK_FITTING used to fit the full spectrum;
;                          now it only fits in a window around the (set
;                          of) emission line(s) to be fit.  The window
;                          is hard-wired (as in Enci's code) to be +/-
;                          25 angstroms about the range in rest
;                          wavelength of the fitted lines.  Add a
;                          keyword to fix the baseline to zero.
;       20 Oct 2015: (KBW) Add a second fit of the OII doublet that
;                          treats it as a single line.  The output model
;                          (eml_model) provides only the single Gaussian
;                          fit to the OII doublet!  Formatting edits.
;                          Added gausmod function.  Added functionality
;                          to MDAP_DEFINE_EMISSION_LINES_BELFIORE to
;                          indicate whether or not the dispersions of
;                          two "connected" lines should be tied during
;                          the fitting.  Implemented an iterative fit to
;                          avoid status=2 errors from MPFIT, which I've
;                          taken to mean that the input and output
;                          parameters are the same.  Change to double
;                          precision. Correction to output velocity
;                          dispersion measurement.
;-
;------------------------------------------------------------------------------


;----------------------------------------------------------------------------
; Define the lines to fit
PRO MDAP_DEFINE_EMISSION_LINES_BELFIORE, $
                n_l, w_l, fit_l, con_l, sig_l

        ; line names
        n_l = [ 'OII_3727', 'OII_3729', 'OIId_3727', 'Hb', 'OIII_4959', 'OIII_5007', 'OI_6300', $
                'OI_6363', 'NII_6548', 'Ha', 'NII_6583', 'SII_6716', 'SII_6731' ]
        ; flag to fit line (0-no;1-yes)
        fit_l = [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 ]
        ; define sets of lines to be fit simultaneously; velocities will be tied; running index
        con_l = [ 3, 3, 0, 0, 2, 2, 4, 4, 1, 1, 1, 5, 5 ]
        ; flag to tie dispersions (as well as velocities) for simultaneous line fits (0-no;1-yes)
        sig_l = [ 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ]

        ; Force the rest wavelengths to be the same between Enci and Belfiore fits
        eml_par = MDAP_DEFINE_EMISSION_LINES_ENCI_BELFIORE()
        w_l = eml_par.lambda

END
;----------------------------------------------------------------------------


;----------------------------------------------------------------------------
;function used to generate model
function gausmod, p, x=x, nlines2=nlines2, wl=wl, redshift=redshift

    c=299792.458d  ; speed of light in KM/s
    g_tot=make_array(n_elements(x), /double, value=p[3*nlines2])      ; initialize to background
    for i=0, nlines2-1 do begin
        g_p= p[3*i+2]*exp( -( (  x - wl[i]*(1.0+p[3*i]/c+redshift) )^2 )/ (2.*p[3*i+1]^2)  )
        g_tot=g_tot+g_p
    endfor

    return, g_tot
end
;----------------------------------------------------------------------------
;----------------------------------------------------------------------------
;function defining one or more Gaussians, input for MPFIT
function gausfit, p, x=x, y=y, err=err, nlines2=nlines2, wl=wl, redshift=redshift
    g_tot = gausmod(p, x=x, nlines2=nlines2, wl=wl, redshift=redshift)
    return, ((y-g_tot)/err)
end
;----------------------------------------------------------------------------


;----------------------------------------------------------------------------
; procedure performing the actual line fitting using MPFIT
PRO peak_fitting, resid, noise_, wave, star_sigma, redshift, w_l, n_l, nlines2, l, $
                  tie_sig=tie_sig, zero_base=zero_base, eflg=eflg, quiet=quiet

    c=299792.458d                   ; Speed of light in km/s
  
    ; Build the output structure
    L=REPLICATE({n:' ', wl:0.0, ws:0.0, we:0.0,  b:0.0,  s:0.0,  a:0.0,  v:0.0, $
                                                eb:0.0, es:0.0, ea:0.0, ev:0.0}, nlines2)
    ; The spectral range to fit:
    bin = 25.                       ; Half-width of the fitting window
    wave_min = min(w_l)-bin         ; Wavelength range to fit
    wave_max = max(w_l)+bin
;   print, wave_min, wave_max
    fit_indx = where( wave gt wave_min and wave lt wave_max, count) ; Pixels to include in the fit
    eflg = 0
    if count eq 0 then begin
        print, 'No relevant wavelengths in spectrum!'
        eflg = 1
        return
    endif

    ; Define some hard limits
    mask_width = 600.               ; mask used for guess parameters in km/s
    maxvel = 500.                   ; maximum allowed velocity
    velrange = [-maxvel, maxvel]    ; allowed velocity range
    sigrange = [ 50., 450. ]        ; allowed velocity dispersion range
    default_sigma = 80.             ; default velocity dispersion
    default_amplitude = 1.          ; default amplitude
;   xtol = 1d-6
  
    ; ------------------------------------------------------------------
    ; **** Set Parameters for MPFIT
    err=noise_[fit_indx]
    x = wave[fit_indx]
    y = resid[fit_indx]

;   plot, wave, resid
;   oplot, x, y, color=200
;   stop
  
    ; *** Initialize the output structure values
    l.n=N_L
    l.wl=W_L
    l[*].ws = wave_min
    l[*].we = wave_max

    ; ------------------------------------------------------------------
    ; *** Set emission lines starting parameters
    ;
    ; emission lines sigma and velocity are set equal to those of the stars
    ;
    ; NOTE sometimes the stars are badly fitted, especially when the
    ; continuum is very low. In this case the velocity/sigma for the stars
    ; can  be nonsense.
    dw=mask_width*l.wl/c                    ;width of the mask in angstrom

    np_total = 3*nlines2+1
    p_start=dblarr(np_total)           ; Last element is the baseline
    
    parinfo=replicate({fixed:0,limited:[0,0],limits:[0.,0.], tied:''}, np_total)

    for i=0, nlines2-1 do begin
        if star_sigma ge sigrange[0] && star_sigma le sigrange[1] then begin
            sigma2=l[i].wl*star_sigma/c
        endif else $
            sigma2=l[i].wl*default_sigma/c

        ;subscripts of the wav range within the mask
        wrange=where(x gt l[i].wl*(1+redshift)-dw[i] and x lt l[i].wl*(1+redshift)+dw[i], count)
        ; TODO: this error case needs to be handled better
        if count eq 0 then begin
            a2 = 0.0
            v2 = redshift*c
            continue
        endif
        ;max of the flux within the mask, and its subscript (w_max)
        a2=max(y[wrange], w_max)

;       print, l[i].wl*(1+redshift)-dw[i], l[i].wl*(1+redshift)+dw[i]
;       print, y[wrange[w_max]], a2[i]

;       plot, x, y
;       oplot, x[wrange], y[wrange], color='00FF00'x
;       oplot, [x[wrange[0]]], [y[wrange[0]]], psym=2, color='FF0000'x
;       oplot, [x[wrange[count-1]]], [y[wrange[count-1]]], psym=2, color='FF0000'x
;       oplot, [x[wrange[w_max]]], [a2[i]], psym=2, color='0000FF'x
;       stop

        if a2 le 0.0 then $
            a2 = default_amplitude

        ;wavelength of max subscript
        x_max = x[wrange[w_max]]
        ; velocity at max flux
        v2 = (x_max - l[i].wl*(1+redshift))/(l[i].wl*(1+redshift))*c
        if abs(v2) gt maxvel then $
            v2 = 0.

        ; Impose velocity range
        ; TODO: expand this range?
        p_start[3*i]   =  v2
        parinfo[3*i].limited=[1,1]
        parinfo[3*i].limits=velrange

        ; Impose velocity dispersion range
        ; TODO: expand this range?
        p_start[1+3*i] =  sigma2
        parinfo[1+3*i].limited=[1,1]
        parinfo[1+3*i].limits=l[i].wl*sigrange/c

        ; All amplitudes need to be positive!
        p_start[2+3*i] =  a2
        parinfo[2+3*i].limited=[1,0]
        parinfo[2+3*i].limits=[0.,0]

    endfor

;   print, p_start
;   gstrt = gausmod(p_start, x=x, nlines2=nlines2, wl=l.wl, redshift=redshift)
;   plot, x, y
;   oplot, x, gstrt, color='FF0000'x

    ; ------------------------------------------------------------------
    ; *** Apply additional constraints

;    ; Tie the relative amplitude and sigmas and velocities of multiplet
;    ; lines (here only [OIII] and [NII])
;    wOIII_4959 = where(l.n eq 'OIII_4959', count_4959)
;    wOIII_5007 = where(l.n eq 'OIII_5007', count_5007)
;    
;    ; KBW: To tie fluxes, must tie both amplitude and velocity dispersion!
;    if count_4959 ne 0 && count_5007 ne 0 then begin
;        ;ratio of amplitudes fixed to theoretical
;        parinfo[3*wOIII_4959+2].tied=strcompress('0.330*p['+string(3*wOIII_5007+2)+']',/remove_all)
;        ;fix to have same velocities
;        parinfo[3*wOIII_4959].tied=strcompress('p['+string(3*wOIII_5007)+']',/remove_all)
;        ;fix to have same sigma
;        parinfo[3*wOIII_4959+1].tied=strcompress('p['+string(3*wOIII_5007+1)+']',/remove_all)
;    endif
;    
;    WNII_6548 = where(l.n eq 'NII_6548', count_6548)
;    wNII_6583 = where(l.n eq 'NII_6583', count_6583)
;
;    if count_6548 ne 0 && count_6583 ne 0 then begin
;        ;ratio of amplitudes fixed to theoretical
;        parinfo[3*WNII_6548+2].tied=strcompress('0.340*p['+string(3*wNII_6583+2)+']',/remove_all)
;        ;fix to have same velocities
;        parinfo[3*WNII_6548].tied=strcompress('p['+string(3*wNII_6583)+']',/remove_all)
;        ;fix to have same sigma
;        parinfo[3*WNII_6548+1].tied=strcompress('p['+string(3*wNII_6583+1)+']',/remove_all)
;    endif
    
    ; if there is more than one line, tie all the velocities, and tie the dispersions if requested
    if nlines2 gt 1 then begin
        for i=0, nlines2-2 do begin
            parinfo[3*(i+1)].tied=strcompress('p[0]', /remove_all)
            p_start[3*(i+1)] = p_start[0]
            if n_elements(tie_sig) gt 0 && tie_sig gt 0 then begin
                parinfo[1+3*(i+1)].tied=strcompress('p[1]', /remove_all)
                p_start[1+3*(i+1)] = p_start[1]
            endif
        endfor
    endif
    
;   gstrt = gausmod(p_start, x=x, nlines2=nlines2, wl=l.wl, redshift=redshift)
;   oplot, x, gstrt, color='0000FF'x, linestyle=2

    if keyword_set(zero_base) then begin
        ; Force the baseline to be zero
        p_start[np_total-1] = 0.
        parinfo[np_total-1].fixed=1
    endif else begin
        ; Otherwise set it to the sigma-clipped mean within the spectral range to fit
        MEANCLIP, y, meanbg, clipsig=3.
        p_start[np_total-1] = meanbg
    endelse

;   print, parinfo.tied
   
    ; ------------------------------------------------------------------
    ; *** Fit function

    ; To avoid getting the same thing that's put in, I run the fit in a
    ; while...do loop for a maximum of 5 iterations
    status = 2
    ntry = 0
    p_start_orig = p_start
    while (ntry lt 10 && status eq 2) do begin
        p_start = (1.0+0.1*randomu(seed, np_total))*p_start_orig
        for i=0,np_total-1 do begin
            if parinfo[i].limited[0] eq 1 && p_start[i] lt parinfo[i].limits[0] then $
                p_start[i] = parinfo[i].limits[0]+0.01*abs(parinfo[i].limits[0])
            if parinfo[i].limited[1] eq 1 && p_start[i] gt parinfo[i].limits[1] then $
                p_start[i] = parinfo[i].limits[1]-0.01*abs(parinfo[i].limits[1])
        endfor

;       print, p_start_orig
;       print, p_start
        status=0
        p = mpfit('gausfit',p_start, functargs={x:x,y:y,err:err, nlines2:nlines2, wl:l.wl, $
                  redshift:redshift}, bestnorm=chi2, parinfo=parinfo, perror=perror, $
                  status=status, errmsg=errmsg, quiet=1);, xtol=xtol)
        ntry++
    endwhile
    if ntry eq 10 then $
        print, 'WARNING: peak_fitting(): Hit maximum number of fit restarts!'

    ; TODO: Only print this if !quiet
    if status eq 5 then $
        print, 'WARNING: Maximum number of iterations reached in MPFIT()!'

;   print, p

;   gmod = gausmod(p, x=x, nlines2=nlines2, wl=l.wl, redshift=redshift)
;   oplot, x, gmod, color='00FF00'x
;   stop
    for i=0,np_total-1 do begin
        all_same = 1
        if abs(p[i]-p_start[i]) gt 1e-4 then begin
            all_same = 0
            break
        endif
    endfor
    if all_same eq 1 then begin
        print, 'WARNING: Input and output parameters are all within 1e-4!'
;       print, status
;       print, errmsg
        print, l.n
    endif

    ; Expects that the input and output parameters are the same if
    ; status=2, warn if not true
    if status eq 2 && all_same eq 0 then $
        print, 'WARNING: status is 2 but the input/output parameters are different'

    eflg = 0
    if status LE 0 then begin
        print, errmsg
;       print, parinfo[*].limited[0]
;       print, parinfo[*].limits[0]
;       print, p_start
;       print, parinfo[*].limits[1]
;       print, parinfo[*].limited[1]
        eflg = 1
        return
    endif

    ; ------------------------------------------------------------------
    ; *** Calculate the reduced chi^2
    ; TODO: Return this?
    n_data = n_elements(x);lambda)          ; Number of fitted data points
    np_free = np_total
    if keyword_set(zero_base) then $
        np_free -= 1                        ; Baseline was fixed during fit
    chi2_reduced = chi2/(n_data-np_free-1)  ; !! NEVER USED !!
    
    ; *** Save Output
    for i=0, nlines2-1 do begin
;        l[i].s=p[1+3*i]*c/l[i].wl
;        l[i].es = perror[3*i+1]*c/l[i].wl
        lobs = l[i].wl*(1.+p[3*i]/c+redshift)   ; fitted wavelength
        l[i].s=p[1+3*i]*c/lobs                  ; Velocity dispersion
        l[i].es = perror[1+3*i]*c/lobs
        l[i].v=  p[3*i]                         ; Velocity relative to the input redshift
        l[i].ev = perror[3*i]
        l[i].a=p[2+3*i]                         ; Amplitude
        l[i].ea=perror[2+3*i]

        ; Same baseline saved for all lines
        if ~keyword_set(zero_base) then begin
            l[i].b=p[np_total-1]
            l[i].eb=perror[np_total-1]
        endif
      
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
                zero_base=zero_base, quiet=quiet

    ; Define the lines to fit; !! HARD-CODED !!
    MDAP_DEFINE_EMISSION_LINES_BELFIORE, n_l, w_l, fit_l, con_l, sig_l
    nlines2 = n_elements(n_l)                       ; Number of lines
  
    ; Output structure
    El = {  name:strarr(nlines2), lambda:dblarr(nlines2),   Lmin:dblarr(nlines2), $
            Lmax:dblarr(nlines2),   Base:dblarr(nlines2),   Ampl:dblarr(nlines2), $
             Vel:dblarr(nlines2),  Sigma:dblarr(nlines2),  eBase:dblarr(nlines2), $
           eAmpl:dblarr(nlines2),   eVel:dblarr(nlines2), eSigma:dblarr(nlines2) }

    indx = where(mask lt 1., count)
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
   
    ; *** Call the emission line fitting procedure ***
    
    ;*** Fit the untied lines
    if count ne 0 then begin
      
    for kl=0, n_elements(untied)-1 do begin
      s_l=untied[kl]
      wl_m=w_l[untied[kl]]
      nl_m=n_l[untied[kl]]
      
      ; *** if the line is strong (fit_l = 1) then always fit it ***
    
        peak_fitting, galaxy_eml_only[indx], flux_err[indx], wave[indx], star_sigma, redshift, $
                      wl_m, nl_m, 1, l, zero_base=zero_base, eflg=eflg, quiet=quiet
        if eflg eq 1 then begin
            El.name[s_l]=l[0].n
            El.lambda[s_l]=l[0].wl
            El.Lmin[s_l]=l[0].ws
            El.Lmax[s_l]=l[0].we
            El.base[s_l]=0.
            El.ampl[s_l]=0.
            El.vel[s_l]=0.
            El.sigma[s_l]=0.
            El.ebase[s_l]=0.
            El.eampl[s_l]=1.
            El.evel[s_l]=1.
            El.esigma[s_l]=1.
        endif else begin
            El.name[s_l]=l[0].n
            El.lambda[s_l]=l[0].wl
            El.Lmin[s_l]=l[0].ws
            El.Lmax[s_l]=l[0].we
            El.base[s_l]=l[0].b
            El.ampl[s_l]=l[0].a
            El.vel[s_l]=l[0].v
            El.sigma[s_l]=l[0].s
            El.ebase[s_l]=l[0].eb
            El.eampl[s_l]=l[0].ea
            El.evel[s_l]=l[0].ev
            El.esigma[s_l]=l[0].es
            if ~keyword_set(quiet) then $
                print, 'Successfully fitted ', nl_m, ' untied lines'
        endelse
    endfor
    
    endif
    ;now I've fitted all the untied lines

    
    ;*** Fitting of the kinematically tied lines (only the velocity is tied!)***
    if ~keyword_set(quiet) then $
        print, 'Number of tied lines:', n_tied
    for kl=1, n_tied do begin
      wgroup=where(con_l eq kl, count)
      while count eq 0 do begin
        if ~keyword_set(quiet) then $
            print, 'Error in your input emission line file!'
        kl=kl+1
        wgroup=where(con_l eq kl, count)
      endwhile
      
      wl_m=w_l[wgroup]
      nl_m=n_l[wgroup]
      n_m=N_elements(wgroup)

      if total(sig_l[wgroup]) gt 0 then begin
        tie_sig = 1
      endif else $
        tie_sig = 0
      
        peak_fitting, galaxy_eml_only[indx], flux_err[indx], wave[indx], star_sigma, redshift, $
                      wl_m, nl_m, n_m, l, tie_sig=tie_sig, zero_base=zero_base, eflg=eflg, $
                      quiet=quiet
        
        if eflg eq 1 then begin
            El.name[wgroup]=l.n
            El.lambda[wgroup]=l.wl
            El.Lmin[wgroup]=l.ws
            El.Lmax[wgroup]=l.we
            El.base[wgroup]=0.
            El.ampl[wgroup]=0.
            El.vel[wgroup]=0.
            El.sigma[wgroup]=0.
            El.ebase[wgroup]=1.
            El.eampl[wgroup]=1.
            El.evel[wgroup]=1.
            El.esigma[wgroup]=1.
        endif else begin
            El.name[wgroup]=l.n
            El.lambda[wgroup]=l.wl
            El.Lmin[wgroup]=l.ws
            El.Lmax[wgroup]=l.we
            El.base[wgroup]=l.b
            El.ampl[wgroup]=l.a
            El.vel[wgroup]=l.v
            El.sigma[wgroup]=l.s
            El.ebase[wgroup]=l.eb
            El.eampl[wgroup]=l.ea
            El.evel[wgroup]=l.ev
            El.esigma[wgroup]=l.es
            if ~keyword_set(quiet) then $
                print, 'I fitted the ', nl_m, ' lines'
        endelse
    endfor
   
   ; Get the full model of the emission-line spectrum (ignoring masks and the baseline)
   c=299792.458d  ; speed of light in KM/s
    eml_model=wave*0.0
    ; ignore the OII line fit as two separate lines
    startk=2
   for kk=startk, nlines2-1 do begin
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





