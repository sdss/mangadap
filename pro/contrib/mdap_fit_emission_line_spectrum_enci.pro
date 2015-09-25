;+
; NAME:
;       MDAP_FIT_EMISSION_LINE_SPECTRUM_ENCI
;
; PURPOSE:
;       Enci's code to fit an emission-line only spectrum.
;
; CALLING SEQUENCE:
;       MDAP_FIT_EMISSION_LINE_SPECTRUM_ENCI, shift0, lambda, pure_emi, err, mask, fitwin, result, $
;                                             emfit, OII3727=OII3727, Hb4861=Hb4861, $
;                                             OIII5007=OIII5007, OI6300=OI6300, Ha6563=Ha6563, $
;                                             SII6717=SII6717, zero_base=zero_base, plotsp=plotsp
;
; INPUTS:
;       shift0 double
;               Input guess redshift of the spectrum.
;
;       lambda dblarr[C]
;               Wavelength at each of the C pixels.
;
;       pure_emi dblarr[C] 
;               Flux in each of C pixels.
;
;       err dblarr[C] 
;               Flux error in each of the C pixels.
;
;       mask dblarr[C]
;               Bad pixel mask for spectrum (1-bad; 0-good)
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;       /OII3727
;               Flag to fit the OII lines (simultaneously fit both)
;
;       /Hb4861
;               Flag to fit the H-beta line
;
;       /OIII5007
;               Flag to fit the OIII lines (simultaneously fit both)
;
;       /OI6300
;               Flag to fit the OI lines (simultaneously fit both)
;
;       /Ha6563
;               Flag to fit the H-alpha and NII lines (simultaneously
;               fit all three).
;
;       /SII6717
;               Flag to fit the SII lines (simultaneously fit both)
;
;       /zero_base
;               Force the baseline to be 0
;
;       /plotsp
;               Flag to show plots of the fit to each line.
;
; OUTPUT:
;       fitwin structure
;               Structure containing the limits of the window used
;               during the fit of each (set of) line(s).
;
;       result structure
;               Structure containing the fitting results.  Each line has
;               a 6-element dblarr in the structure.  The elements are
;               the three parameters of the fitted Gaussian using the
;               IDL task GAUSS1(): (line centroid, sigma, and area=flux)
;               and the formal (MPFIT) errors in those parameters.
;
;       emfit dblarr[C]
;               The best-fitting emission-line spectrum.
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
;       10 Dec 2014: Incorporation of code written by Enci into the DAP
;                    by K. Westfall (KBW).  Allow for a mask.
;       13 Aug 2015: (KBW) Fit OII as a doublet, and add the fit to the
;                          OI doublet.  Added 'status' and 'errmsg' to
;                          calls to MPFITEXPR().  If MPFITEXPR() ends in
;                          error, the error is printed and the fitted
;                          flux is set to 0.  Also checks if status = 5,
;                          indicating that MPFITEXPR() reached the
;                          maximum number of iterations.
;       16 Sep 2015: (KBW) Allow for non-zero baseline; returned with
;                          fit
;       17 Sep 2015: (KBW) Return wavelength limits of the fitting
;                          window used for each (set of) line(s)
;-
;------------------------------------------------------------------------------

PRO MDAP_FIT_EMISSION_LINE_SPECTRUM_ENCI, $
                shift0, lambda, pure_emi, err, mask, fitwin, result, emfit, OII3727=OII3727, $
                Hb4861=Hb4861, OIII5007=OIII5007, OI6300=OI6300, Ha6563=Ha6563, SII6717=SII6717, $
                zero_base=zero_base, plotsp=plotsp

;pro emission_fit1,shift0,lambda,pure_emi,err,result,emfit,$
;        OII3727=OII3727,Hb4861=Hb4861,OIII5007=OIII5007,Ha6563=Ha6563,$
;		SII6717=SII6717,plotsp=plotsp
;
; We fit all these lines with x0+guassian function.
; In the case the continue component is not subtracted well, our fittings 
; still make sence. 

; sigma are forced between 50km/s to 450km/s.
; input sigma 100 km/s

; velocity are forced between -400km/s to 400km/s
; input velocity = 0 km/s

; All fittings are free, except for the bound of OIII and NII doublets. 
; 

; Save the wavelength limits of the fitting window(s)
fitwin={  OII3727:fltarr(2), OII3729:fltarr(2),  Hb4861:fltarr(2), OIII4959:fltarr(2), $
         OIII5007:fltarr(2),  OI6300:fltarr(2),  OI6363:fltarr(2),  NII6548:fltarr(2), $
           Ha6563:fltarr(2), NII6583:fltarr(2), SII6717:fltarr(2),  SII6731:fltarr(2) }

; Number of parameters to save: 2*(3 gauss parameters + baseline) (value and error)
npar = 8
result={  OII3727:fltarr(npar), OII3729:fltarr(npar),  Hb4861:fltarr(npar), OIII4959:fltarr(npar), $
         OIII5007:fltarr(npar),  OI6300:fltarr(npar),  OI6363:fltarr(npar),  NII6548:fltarr(npar), $
           Ha6563:fltarr(npar), NII6583:fltarr(npar), SII6717:fltarr(npar),  SII6731:fltarr(npar) }

; Model spectrum (without baseline!)
emfit=pure_emi*0.0

shift=shift0+1
; Speed of light in km/s
c0=299792.458d

bin=25

if keyword_set(plotsp) then begin
        !p.multi=[0,5,1,0,0]
endif

;;----------------------------------------OII3727
;if keyword_set(OII3727) then begin
;   indx=where(lambda gt 3727-bin and lambda lt 3727+bin and mask lt 0.5)
;   if n_elements(indx) gt 5 then begin
;   lam1=lambda[indx]
;   emi1=pure_emi[indx]
;   err1=err[indx]
;        meanclip,emi1,meansp,sig,clipsig=3
;        meansp=0.0
;
;   expr1='p[0]+gauss1(x,p[1:3])'
;   par=replicate({value:0.D, fixed:0,limited:[0,0], limits:[0.D,0.D]},4)
;   par(*).value = [meansp,3727.*shift,3727.0*100.0/c0,600.]
;   par(0).fixed=1
;   par(1).limited(0) = 1
;   par(1).limits(0) = par(1).value-3727.0*400.0/c0
;   par(1).limited(1) = 1
;   par(1).limits(1) = par(1).value+3727.0*400.0/c0
;   par(2).limited(0) = 1
;   par(2).limits(0) = 3727.0*50.0/c0
;   par(2).limited(1) = 1
;   par(2).limits(1) = 3727.0*450.0/c0
;   par(3).limited(0)=1
;   par(3).limits(0)=0
;   result1=mpfitexpr(expr1,lam1,emi1,err1,parinfo=par,perror=perror,bestnorm=bestnorm,yfit=yfit,/quiet)
;        result.OII3727=[result1[1:3],perror[1:3]]
;
;   ;emfit[indx]=yfit
;   ; KBW: allow for mask, but then provide best fitting Gaussian over full window
;   indx=where(lambda gt 3727-bin and lambda lt 3727+bin, count)
;   if count ne 0 then $
;       emfit[indx] = emfit[indx] + result1[0]+GAUSS1(lambda[indx],result1[1:3])
;if keyword_set(plotsp) then begin
;
;   plot,lam1,emi1,title='OII3727'
;   ;oplot,lam1,result1[0]+gauss1(lam1,result1[1:3]),color=djs_icolor('red')
;        oplot,lam1,yfit,color=djs_icolor('red'),linestyle=2
;endif
;
;endif ;hb
;endif

;---------------------------------OII_3727+3729----------------------
if keyword_set(OII3727) then begin

    ; Set the fitting window; same for both lines
    fitwin.OII3727 = [ 3727-bin, 3730+bin ]
    fitwin.OII3729 = [ 3727-bin, 3730+bin ]

   indx=where(lambda gt fitwin.OII3727[0] and lambda lt fitwin.OII3727[1] and mask lt 0.5)
   if n_elements(indx) gt 5 then begin
   lam1=lambda[indx]
   emi1=pure_emi[indx]
   err1=err[indx]
        meanclip,emi1,meansp,sig,clipsig=3
        ; Initial guess background is 0.0
        meansp=0.0

   expr1='p[0]+gauss1(x,p[1:3])+gauss1(x,p[4:6])'
   par=replicate({value:0.D, fixed:0, tied:'', limited:[0,0], limits:[0.D,0.D]},7)
   par(*).value = [meansp,3727.*shift,3727.*100.0/c0,300.,3730.*shift,3730.*100.0/c0,900.]
   ; Force the baseline to be zero
   if keyword_set(zero_base) then $
        par(0).fixed=1
;   par(3).tied = 'p(6)/3.'
   par(1).limited(0) = 1
   par(1).limits(0) = par(1).value-3727.*400.0/c0
   par(1).limited(1) = 1
   par(1).limits(1) = par(1).value+3727.*400.0/c0
   par(4).limited(0) = 1
   par(4).limits(0) = par(4).value-3730.*400.0/c0
   par(4).limited(1) = 1
   par(4).limits(1) = par(4).value+3730.*400.0/c0
   par(2).limited(0)=1
   par(2).limits(0)=3727.*50.0/c0
   par(2).limited(1) = 1
   par(2).limits(1) = 3727.0*450.0/c0
   par(5).limited(0)=1
   par(5).limits(0)=3730.0*50.0/c0
   par(5).limited(1) = 1
   par(5).limits(1) = 3730.0*450.0/c0
   par(3).limited(0) = 1
   par(3).limits(0) = 0
   par(6).limited(0) = 1
   par(6).limits(0) = 0
   result1=mpfitexpr(expr1, lam1, emi1, err1, parinfo=par, perror=perror, bestnorm=bestnorm, $
                     yfit=yfit, status=status, errmsg=errmsg, /quiet)

   if status gt 0 then begin
;        result.OII3727=[result1[1:3],perror[1:3]]
;        result.OII3729=[result1[4:6],perror[4:6]]
        result.OII3727=[result1[0:3],perror[0:3]]
        result.OII3729=[result1[0],result1[4:6],perror[0],perror[4:6]]
   endif else $
        print, errmsg

    if status eq 5 then $
        print, 'WARNING: Maximum number of iterations reached in MPFITEXPR()!'

   ;emfit[indx]=yfit
   ; KBW: allow for mask, but then provide best fitting Gaussian over full window
   indx=where(lambda gt fitwin.OII3727[0] and lambda lt fitwin.OII3727[1], count)
   if count ne 0 then $
    emfit[indx] = emfit[indx] + GAUSS1(lambda[indx],result1[1:3]) $
                              + GAUSS1(lambda[indx],result1[4:6])
;    emfit[indx] = emfit[indx] + result1[0] + GAUSS1(lambda[indx],result1[1:3]) $
;                              + GAUSS1(lambda[indx],result1[4:6])
if keyword_set(plotsp) then begin
   plot,lam1,emi1,title='OII_3727+3729';,$;ytitle='flux(unit of 10^-15)'
;               position=[0.25,0.35,0.40,0.55]
   ;oplot,lam1,result1[0]+gauss1(lam1,result1[1:3]),color=djs_icolor('red')
   ;oplot,lam1,result1[0]+gauss1(lam1,result1[4:6]),color=djs_icolor('red')
        oplot,lam1,yfit,color=djs_icolor('red'),linestyle=2
endif
endif
endif




;----------------------------------------Hbeta4861
if keyword_set(Hb4861) then begin

    ; Set the fitting window; same for both lines
    fitwin.Hb4861 = [ 4861-bin, 4861+bin ]

   indx=where(lambda gt fitwin.Hb4861[0] and lambda lt fitwin.Hb4861[1] and mask lt 0.5)
   if n_elements(indx) gt 5 then begin
   lam1=lambda[indx]
   emi1=pure_emi[indx]
   err1=err[indx]
        meanclip,emi1,meansp,sig,clipsig=3
        ; Initial guess background is 0.0
        meansp=0.0

   expr1='p[0]+gauss1(x,p[1:3])'
   par=replicate({value:0.D, fixed:0,limited:[0,0], limits:[0.D,0.D]},4)
   par(*).value = [meansp,4861.*shift,4861.0*100.0/c0,600.]
   ; Force the baseline to be zero
   if keyword_set(zero_base) then $
        par(0).fixed=1
   par(1).limited(0) = 1
   par(1).limits(0) = par(1).value-4861.0*400.0/c0
   par(1).limited(1) = 1
   par(1).limits(1) = par(1).value+4861.0*400.0/c0
   par(2).limited(0) = 1
   par(2).limits(0) = 4861.0*50.0/c0
   par(2).limited(1) = 1
   par(2).limits(1) = 4861.0*450.0/c0
   par(3).limited(0)=1
   par(3).limits(0)=0
   result1=mpfitexpr(expr1, lam1, emi1, err1, parinfo=par, perror=perror, bestnorm=bestnorm, $
                     yfit=yfit, status=status, errmsg=errmsg, /quiet)

    if status gt 0 then begin
        result.Hb4861=[result1[0:3],perror[0:3]]
;       result.Hb4861=[result1[1:3],perror[1:3]]
    endif else $
        print, errmsg

    if status eq 5 then $
        print, 'WARNING: Maximum number of iterations reached in MPFITEXPR()!'

   ;emfit[indx]=yfit
   ; KBW: allow for mask, but then provide best fitting Gaussian over full window
   indx=where(lambda gt fitwin.Hb4861[0] and lambda lt fitwin.Hb4861[1], count)
   if count ne 0 then $
        emfit[indx] = emfit[indx] + GAUSS1(lambda[indx],result1[1:3])
;       emfit[indx] = emfit[indx] + result1[0] + GAUSS1(lambda[indx],result1[1:3])
if keyword_set(plotsp) then begin

   plot,lam1,emi1,title='H_beta'
   ;oplot,lam1,result1[0]+gauss1(lam1,result1[1:3]),color=djs_icolor('red')
        oplot,lam1,yfit,color=djs_icolor('red'),linestyle=2
endif

endif ;hb
endif


;---------------------------------OIII5007----------------------
;if (keyword_set(OIII5007) or keyword_set(OIII4959)) then begin
if (keyword_set(OIII5007) || keyword_set(OIII4959)) then begin

    ; Set the fitting window; same for both lines
    fitwin.OIII4959 = [ 4959-bin, 5007+bin ]
    fitwin.OIII5007 = [ 4959-bin, 5007+bin ]

   indx=where(lambda gt fitwin.OIII4959[0] and lambda lt fitwin.OIII4959[1] and mask lt 0.5)
   if n_elements(indx) gt 5 then begin
   lam1=lambda[indx]
   emi1=pure_emi[indx]
   err1=err[indx]
        meanclip,emi1,meansp,sig,clipsig=3
        ; Initial guess background is 0.0
        meansp=0.0

   expr1='p[0]+gauss1(x,p[1:3])+gauss1(x,p[4:6])'
   par=replicate({value:0.D, fixed:0, tied:'', limited:[0,0], limits:[0.D,0.D]},7)
   par(*).value = [meansp,4959.*shift,4959.*100.0/c0,300.,5007.*shift,5007.*100.0/c0,900.]
   ; Force the baseline to be zero
   if keyword_set(zero_base) then $
        par(0).fixed=1
   par(3).tied = 'p[6]/3.'
   par(1).limited(0) = 1
   par(1).limits(0) = par(1).value-4959.*400.0/c0
   par(1).limited(1) = 1
   par(1).limits(1) = par(1).value+4959.*400.0/c0
   par(4).limited(0) = 1
   par(4).limits(0) = par(4).value-5007.*400.0/c0
   par(4).limited(1) = 1
   par(4).limits(1) = par(4).value+5007.*400.0/c0
        par(2).limited(0)=1
        par(2).limits(0)=4959.*50.0/c0
   par(2).limited(1) = 1
   par(2).limits(1) = 4959.0*450.0/c0
        par(5).limited(0)=1
        par(5).limits(0)=5007.0*50.0/c0
   par(5).limited(1) = 1
   par(5).limits(1) = 5007.0*450.0/c0
        par(3).limited(0) = 1
        par(3).limits(0) = 0
        par(6).limited(0) = 1
        par(6).limits(0) = 0
   result1=mpfitexpr(expr1, lam1, emi1, err1, parinfo=par, perror=perror, bestnorm=bestnorm, $
                     yfit=yfit, status=status, errmsg=errmsg, /quiet)

    if status gt 0 then begin
        result.OIII4959=[result1[0:3],perror[0:3]]
        result.OIII5007=[result1[0],result1[4:6],perror[0],perror[4:6]]
;       result.OIII4959=[result1[1:3],perror[1:3]]
;       result.OIII5007=[result1[4:6],perror[4:6]]
    endif else $
        print, errmsg

    if status eq 5 then $
        print, 'WARNING: Maximum number of iterations reached in MPFITEXPR()!'
    
   ;emfit[indx]=yfit
   ; KBW: allow for mask, but then provide best fitting Gaussian over full window
   indx=where(lambda gt fitwin.OIII4959[0] and lambda lt fitwin.OIII4959[1], count)
   if count ne 0 then $
    emfit[indx] = emfit[indx] + GAUSS1(lambda[indx],result1[1:3]) $
                              + GAUSS1(lambda[indx],result1[4:6])
;   emfit[indx] = emfit[indx] + result1[0] + GAUSS1(lambda[indx],result1[1:3]) $
;                             + GAUSS1(lambda[indx],result1[4:6])
if keyword_set(plotsp) then begin
   plot,lam1,emi1,title='OIII_4959+5007';,$;ytitle='flux(unit of 10^-15)'
;               position=[0.25,0.35,0.40,0.55]
   ;oplot,lam1,result1[0]+gauss1(lam1,result1[1:3]),color=djs_icolor('red')
   ;oplot,lam1,result1[0]+gauss1(lam1,result1[4:6]),color=djs_icolor('red')
        oplot,lam1,yfit,color=djs_icolor('red'),linestyle=2
endif
endif
endif



;---------------------------------OI_6302+6365----------------------
if keyword_set(OI6300) then begin

    ; Set the fitting window; same for both lines
    fitwin.OI6300 = [ 6302-bin, 6365+bin ]
    fitwin.OI6363 = [ 6302-bin, 6365+bin ]

   indx=where(lambda gt fitwin.OI6300[0] and lambda lt fitwin.OI6300[1] and mask lt 0.5)
   if n_elements(indx) gt 5 then begin
   lam1=lambda[indx]
   emi1=pure_emi[indx]
   err1=err[indx]
        meanclip,emi1,meansp,sig,clipsig=3
        ; Initial guess background is 0.0
        meansp=0.0

   expr1='p[0]+gauss1(x,p[1:3])+gauss1(x,p[4:6])'
   par=replicate({value:0.D, fixed:0, tied:'', limited:[0,0], limits:[0.D,0.D]},7)
   par(*).value = [meansp,6302.*shift,6302.*100.0/c0,300.,6365.*shift,6365.*100.0/c0,900.]
   ; Force the baseline to be zero
   if keyword_set(zero_base) then $
        par(0).fixed=1
;   par(3).tied = 'p(6)/3.'
   par(1).limited(0) = 1
   par(1).limits(0) = par(1).value-6302.*400.0/c0
   par(1).limited(1) = 1
   par(1).limits(1) = par(1).value+6302.*400.0/c0
   par(4).limited(0) = 1
   par(4).limits(0) = par(4).value-6365.*400.0/c0
   par(4).limited(1) = 1
   par(4).limits(1) = par(4).value+6365.*400.0/c0
   par(2).limited(0)=1
   par(2).limits(0)=6302.*50.0/c0
   par(2).limited(1) = 1
   par(2).limits(1) = 6302.0*450.0/c0
   par(5).limited(0)=1
   par(5).limits(0)=6365.0*50.0/c0
   par(5).limited(1) = 1
   par(5).limits(1) = 6365.0*450.0/c0
   par(3).limited(0) = 1
   par(3).limits(0) = 0
   par(6).limited(0) = 1
   par(6).limits(0) = 0
   result1=mpfitexpr(expr1, lam1, emi1, err1, parinfo=par, perror=perror, bestnorm=bestnorm, $
                     yfit=yfit, status=status, errmsg=errmsg, /quiet)

    if status gt 0 then begin
        result.OI6300=[result1[0:3],perror[0:3]]
        result.OI6363=[result1[0],result1[4:6],perror[0],perror[4:6]]
;       result.OI6300=[result1[1:3],perror[1:3]]
;       result.OI6363=[result1[4:6],perror[4:6]]
    endif else $
        print, errmsg

    if status eq 5 then $
        print, 'WARNING: Maximum number of iterations reached in MPFITEXPR()!'

   ;emfit[indx]=yfit
   ; KBW: allow for mask, but then provide best fitting Gaussian over full window
   indx=where(lambda gt fitwin.OI6300[0] and lambda lt fitwin.OI6300[1], count)
   if count ne 0 then $
    emfit[indx] = emfit[indx] + GAUSS1(lambda[indx],result1[1:3]) $
                              + GAUSS1(lambda[indx],result1[4:6])
;   emfit[indx] = emfit[indx] + result1[0] + GAUSS1(lambda[indx],result1[1:3]) $
;                             + GAUSS1(lambda[indx],result1[4:6])
if keyword_set(plotsp) then begin
   plot,lam1,emi1,title='OII_6300+6363';,$;ytitle='flux(unit of 10^-15)'
;               position=[0.25,0.35,0.40,0.55]
   ;oplot,lam1,result1[0]+gauss1(lam1,result1[1:3]),color=djs_icolor('red')
   ;oplot,lam1,result1[0]+gauss1(lam1,result1[4:6]),color=djs_icolor('red')
        oplot,lam1,yfit,color=djs_icolor('red'),linestyle=2
endif
endif
endif


;----------------------------------------------------------Halpha + NII
;if keyword_set(Ha6563) or keyword_set(NII6583) then begin
if keyword_set(Ha6563) || keyword_set(NII6583) then begin

    ; Set the fitting window; same for both lines
    fitwin.NII6548 = [ 6548-bin, 6583+bin ]
    fitwin.Ha6563 =  [ 6548-bin, 6583+bin ]
    fitwin.NII6583 = [ 6548-bin, 6583+bin ]

   indx=where(lambda gt fitwin.NII6548[0] and lambda lt fitwin.NII6548[1] and mask lt 0.5)
   if n_elements(indx) gt 5 then begin
   lam1=lambda[indx]
   emi1=pure_emi[indx]
   err1=err[indx]
        meanclip,emi1,meansp,sig,clipsig=3
        ; Initial guess background is 0.0
        meansp=0.0
;   print,'Halpha, meansp:',meansp
   expr1='p[0]+gauss1(x,p[1:3])+gauss1(x,p[4:6])+gauss1(x,p[7:9])'
   par=replicate({value:0.D, fixed:0, tied:' ', limited:[0,0], limits:[0.D,0.D]},10)
   par(*).value = [meansp,6548.*shift,6548.*100.0/c0,200.,6563*shift,6563.*100.0/c0,2000.,6583.*shift,6583.0*100.0/c0,500.]
   ; Force the baseline to be zero
   if keyword_set(zero_base) then $
        par(0).fixed=1
   par(3).tied = '0.348116 * p[9]'
   par(1).limited(0) = 1
   par(1).limits(0) = par(1).value-6548.*400/c0
   par(1).limited(1) = 1
   par(1).limits(1) = par(1).value+6548.*400/c0
        par(2).limited(0)=1
        par(2).limits(0)=6548.0*50/c0
   par(2).limited(1) = 1
   par(2).limits(1) = 6548.0*450/c0
        par(3).limited(0)=1
        par(3).limits(0)=0
        par(4).limited(0) = 1
        par(4).limits(0) = par(4).value-6563.*400/c0
        par(4).limited(1) = 1
        par(4).limits(1) = par(4).value+6563.*400/c0
        par(5).limited(0)=1
        par(5).limits(0)=6563.0*50/c0
        par(5).limited(1) = 1
   par(5).limits(1) = 6563.0*450/c0
        par(6).limited(0)=1
        par(6).limits(0)=0
   par(7).limited(0) = 1
   par(7).limits(0) = par(7).value-6583.*400/c0
   par(7).limited(1) = 1
   par(7).limits(1) = par(7).value+6583.*400/c0
        par(8).limited(0) = 1
   par(8).limits(0) = 6583.0*50/c0
   par(8).limited(1) = 1
   par(8).limits(1) = 6583.0*450/c0
        par(9).limited(0)=1
        par(9).limits(0)=0

   result1=mpfitexpr(expr1, lam1, emi1, err1, parinfo=par, perror=perror, bestnorm=bestnorm, $
                     yfit=yfit, status=status, errmsg=errmsg, /quiet)

;   TODO: Sort Gaussians by wavelength to make sure fits ordered properly

    if status gt 0 then begin
        result.NII6548=[result1[0:3],perror[0:3]]
        result.Ha6563=[result1[0],result1[4:6],perror[0],perror[4:6]]
        result.NII6583=[result1[0],result1[7:9],perror[0],perror[7:9]]
;       result.NII6548=[result1[1:3],perror[1:3]]
;       result.Ha6563=[result1[4:6],perror[4:6]]
;       result.NII6583=[result1[7:9],perror[7:9]]
    endif else $
        print, errmsg

    if status eq 5 then $
        print, 'WARNING: Maximum number of iterations reached in MPFITEXPR()!'

   ;emfit[indx]=yfit
   ; KBW: allow for mask, but then provide best fitting Gaussian over full window
   indx=where(lambda gt fitwin.NII6548[0] and lambda lt fitwin.NII6548[1], count)
   if count ne 0 then $
    emfit[indx] = emfit[indx] + GAUSS1(lambda[indx],result1[1:3]) $
                              + GAUSS1(lambda[indx],result1[4:6]) $
                              + GAUSS1(lambda[indx],result1[7:9])
;   emfit[indx] = emfit[indx] + result1[0] + GAUSS1(lambda[indx],result1[1:3]) $
;                 + GAUSS1(lambda[indx],result1[4:6]) + GAUSS1(lambda[indx],result1[7:9])
if keyword_set(plotsp) then begin
   plot,lam1,emi1,title='Ha+NII'
   ;oplot,lam1,result1[0]+gauss1(lam1,result1[1:3]),color=djs_icolor('red')
   ;oplot,lam1,result1[0]+gauss1(lam1,result1[4:6]),color=djs_icolor('red')
   ;oplot,lam1,result1[0]+gauss1(lam1,result1[7:9]),color=djs_icolor('red')
   oplot,lam1,yfit,color=djs_icolor('red'),linestyle=2
endif
endif
endif

;---------------------------------SII_6717+6731----------------------
;if (keyword_set(SII6717) or keyword_set(SII6731)) then begin
if (keyword_set(SII6717) || keyword_set(SII6731)) then begin

    ; Set the fitting window; same for both lines
    fitwin.SII6717 = [ 6717-bin, 6731+bin ]
    fitwin.SII6731 = [ 6717-bin, 6731+bin ]

   indx=where(lambda gt fitwin.SII6717[0] and lambda lt fitwin.SII6717[1] and mask lt 0.5)
   if n_elements(indx) gt 5 then begin
   lam1=lambda[indx]
   emi1=pure_emi[indx]
   err1=err[indx]
        meanclip,emi1,meansp,sig,clipsig=3
        ; Initial guess background is 0.0
        meansp=0.0

   expr1='p[0]+gauss1(x,p[1:3])+gauss1(x,p[4:6])'
   par=replicate({value:0.D, fixed:0, tied:'', limited:[0,0], limits:[0.D,0.D]},7)
   par(*).value = [meansp,6717.*shift,6717.*100.0/c0,300.,6731.*shift,6731.*100.0/c0,900.]
   ; Force the baseline to be zero
   if keyword_set(zero_base) then $
        par(0).fixed=1
;   par(3).tied = 'p(6)/3.'
   par(1).limited(0) = 1
   par(1).limits(0) = par(1).value-6717.*400.0/c0
   par(1).limited(1) = 1
   par(1).limits(1) = par(1).value+6717.*400.0/c0
   par(4).limited(0) = 1
   par(4).limits(0) = par(4).value-6731.*400.0/c0
   par(4).limited(1) = 1
   par(4).limits(1) = par(4).value+6731.*400.0/c0
   par(2).limited(0)=1
   par(2).limits(0)=6717.*50.0/c0
   par(2).limited(1) = 1
   par(2).limits(1) = 6717.0*450.0/c0
   par(5).limited(0)=1
   par(5).limits(0)=6731.0*50.0/c0
   par(5).limited(1) = 1
   par(5).limits(1) = 6731.0*450.0/c0
   par(3).limited(0) = 1
   par(3).limits(0) = 0
   par(6).limited(0) = 1
   par(6).limits(0) = 0
   result1=mpfitexpr(expr1, lam1, emi1, err1, parinfo=par, perror=perror, bestnorm=bestnorm, $
                     yfit=yfit, status=status, errmsg=errmsg, /quiet)

    if status gt 0 then begin
        result.SII6717=[result1[0:3],perror[0:3]]
        result.SII6731=[result1[0],result1[4:6],perror[0],perror[4:6]]
;       result.SII6717=[result1[1:3],perror[1:3]]
;       result.SII6731=[result1[4:6],perror[4:6]]
    endif else $
        print, errmsg

    if status eq 5 then $
        print, 'WARNING: Maximum number of iterations reached in MPFITEXPR()!'
        
   ;emfit[indx]=yfit
   ; KBW: allow for mask, but then provide best fitting Gaussian over full window
   indx=where(lambda gt fitwin.SII6717[0] and lambda lt fitwin.SII6717[1], count)
   if count ne 0 then $
    emfit[indx] = emfit[indx] + GAUSS1(lambda[indx],result1[1:3]) $
                              + GAUSS1(lambda[indx],result1[4:6])
;   emfit[indx] = emfit[indx] + result1[0] + GAUSS1(lambda[indx],result1[1:3]) $
;                             + GAUSS1(lambda[indx],result1[4:6])
if keyword_set(plotsp) then begin
   plot,lam1,emi1,title='SII_6717+6731';,$;ytitle='flux(unit of 10^-15)'
;               position=[0.25,0.35,0.40,0.55]
   ;oplot,lam1,result1[0]+gauss1(lam1,result1[1:3]),color=djs_icolor('red')
   ;oplot,lam1,result1[0]+gauss1(lam1,result1[4:6]),color=djs_icolor('red')
        oplot,lam1,yfit,color=djs_icolor('red'),linestyle=2
endif
endif
endif


end

