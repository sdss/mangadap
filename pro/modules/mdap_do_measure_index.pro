;NAME 
;mdap_do_measure_index
; INPUTS
;spc_      [array]    spectra to calculate the line strenght index
;lambda_    [array]   lambda vector (same number of elements of
;                     spc_). In angstrom
;passband  [array]    2 elements array defining the index passband boundaries
;blue_cont [array]    2 elements array defining the blue
;                     pseudocontinua boundaries
;red_cont  [array]    2 elements array defining the red
;                     pseudocontinua boundaries
;
; OPTIONAL INPUTS:
;norm      [float]    If set, spectra is normalized to its value at
;                     norm angstrom
;title     [string]   Title to write into the plot  
;
;plbound   [array]    Two elements array specifying the boundaries of
;                     the plot (Y axis), which will be set to [plbound[0]*midplot,plbound[1]*midplot]
;                     where midplot is the value of the spectrum at
;                     the wavelenght middle range (=0.5*[min(lambda)+max(lambda)])
;                     Default plbound=[0.6,1.35]; midplot=1.
;
;rebin     [float]    If set, the input spectrum will be rebinned
;                     according to a new step (ang/pixel, defined by
;                     rebin). The starting point of lambda will
;                     remains  unchaged. The input "lambda_" and "spc_"
;                     parameters are not overwritten.
;                     Default: no rebinning.
;noise     [float]    It is useful only when errors need to be
;                     retrieved. Default: noise = sqrt(input spectrum) 
; OPTIONAL KEYWORS
;\noplot              If set, all the plotting commands (plot, oplot,
;                     plots and xyouts) in the routine are not executed
; OUTPUTS:
;ew        [float]    Line equivalent width in angstrom (formula 2 in Worthey et al. 1994)
;index_mag [float]    Line strenght index value in magnitudes (formula 3 in Worthey et al. 1994)
;
;
; OPTIONAL OUTPUTS
;
;errors    [float]    This variable will contain the errors on the
;                     indices computed using Cardiel et al. 1998, A&AS, 127, 597
;                     Equations 41 -46

;HISTORY:
; Ver 0.0 Sept 2013


pro mdap_measure_pseudocontinua,lam,center,spc,pseudo_cont,pseudo_fit
; INTEGRAL MEAN OF THE BAND - the standard definition

;Inputs:
; lam [array] wavelength values of the pseudocontinua
; spc [array] spectra in the wav range of the pseudocontinua

;Outputs:
; pseudo_cont [float] value of the pseudoncontinuum
; pseudo_fit  [array] vector with as many elements as "lam" and equal
;                     to pseudo_cont
 
pseudo_cont=1./(max(lam)-min(lam))*int_tabulated(lam,spc)    ; FORMULA (1) Worthey et al. 1994
pseudo_fit=fltarr(n_elements(lam))*0.+pseudo_cont[0]
end

pro mdap_do_measure_index,spc_,lambda_,passband,blue_cont,red_cont,ew,index_mag,cz=cz,norm=norm,title=title,$
plbound=plbound,noplot=noplot,rebin=rebin,errors=errors,noise=noise,silent=silent

;loadct, 33,/silent
if ~keyword_set(noplot) then begin
   basic_colors, black, white, red, green, blue, yellow, cyan, magenta, orange, mint, purple, pink, olive, lightblue, gray   
   UserSymbol,'Circle',/fill
endif
spc=spc_
lambda=lambda_

;If I have to resample the input spectrum
;print,'Input spectrum has an ang/pixel step of: '+strcompress(string(lambda[1]-lambda[0]),/remove_all)
if n_elements(rebin) ne 0 then begin
   if rebin ne 0 then begin
      NTOT=long((max(lambda)-min(lambda))/rebin)+1
      lll = findgen(NTOT)*rebin+min(lambda)
      sss=interpol(spc,lambda,lll,/spline)
      lambda=lll
      spc=sss
   endif
;stop
;   print,'Input spectrum rebinned. New ang/pixel step is: '+strcompress(string(lambda[1]-lambda[0]),/remove_all)
endif



;- IF I WANT TO NORMALIZE THE SPECTRUM
normalization=1
if n_elements(norm) ne 0 then begin
 normalization=interpol(spc,lambda,norm,/spline)
 spc=spc/normalization[0]
endif
;-

;compute pseudocontinua
ind_blue=where(lambda ge  blue_cont[0] and lambda le blue_cont[1])
ind_red=where(lambda ge  red_cont[0] and lambda le red_cont[1])

blu_lam=mean(blue_cont)
red_lam=mean(red_cont)

lambda_new=findgen(10.*n_elements(ind_blue))*(blue_cont[1]- blue_cont[0])/(10.*n_elements(ind_blue)-1)+blue_cont[0]
spc_new_rebin=interpol(spc,lambda,lambda_new,/spline)
mdap_measure_pseudocontinua,lambda_new,blu_lam,spc_new_rebin,blu_pseudocont,blu_pseudofit ;blue pseudocontinuum value)

lambda_new=findgen(10.*n_elements(ind_red))*(red_cont[1]- red_cont[0])/(10.*n_elements(ind_red)-1)+red_cont[0]
spc_new_rebin=interpol(spc,lambda,lambda_new,/spline)
mdap_measure_pseudocontinua,lambda_new,red_lam,spc_new_rebin,red_pseudocont,red_pseudofit ;red pseudocontinuum value)



if not keyword_set(noplot) then begin
   if n_elements(title) eq 0 then title=''
   if n_elements(plbound) eq 0 then plbound = [0.6,1.35]
   midplot=interpol(spc,lambda,0.5*(min(lambda)+max(lambda)))
   if n_elements(norm) ne 0 then midplot=1.
   
 ;  plot, lambda,spc,/nodata,color=black,xtitle='!6wavelength [angstrom]',yrange=[plbound[0]*midplot,plbound[1]*midplot],ystyle=1,title=title,$
 ;        xrange=[min(lambda)-30,max(lambda)+30],xstyle=1

   inddd=where(lambda ge blue_cont[0]-20 and lambda le red_cont(1)+20)
   mn=min(spc(inddd))
   mx=max(spc(inddd))
   rgn=(mx-mn)*1.05

   yrange=[mn-0.1*rgn,mx+0.1*rgn]
   plot, lambda,spc,/nodata,xtitle='!6wavelength [angstrom]',yrange=yrange,ystyle=1,title=title,$
         xrange=[ blue_cont[0]-10,red_cont(1)+10],xstyle=1;,color=black
   oplot, lambda,spc,thick=2;,color=black
   oplot,lambda[ind_blue] , spc[ind_blue],thick=3;,color=black
   oplot,lambda[ind_blue],blu_pseudofit,color=blue,thick=5
   plots,blu_lam,blu_pseudocont,psym=8;,color=black
   oplot,lambda[ind_red] , spc[ind_red],thick=3;,color=black
   oplot,lambda[ind_red],red_pseudofit,color=red,thick=5
   plots,red_lam,red_pseudocont,psym=8;,color=black
   oplot,[blu_lam,red_lam],[blu_pseudocont,red_pseudocont],linestyle=1;,color=black
endif

indici_band=where(lambda ge passband[0] and lambda le passband[1])

if not keyword_set(noplot) then begin
   oplot,lambda[indici_band],spc[indici_band],thick=5;,color=black
   oplot,[blue_cont[0],blue_cont[0]],[0,interpol(spc,lambda,blue_cont[0],/spline)],linestyle=1,color=blue
   oplot,[blue_cont[1],blue_cont[1]],[0,interpol(spc,lambda,blue_cont[1],/spline)],linestyle=1,color=blue
   oplot,[red_cont[0],red_cont[0]],[0,interpol(spc,lambda,red_cont[0],/spline)],linestyle=1,color=red
   oplot,[red_cont[1],red_cont[1]],[0,interpol(spc,lambda,red_cont[1],/spline)],linestyle=1,color=red
   oplot, [passband[0],passband[0]],[0,interpol(spc,lambda,passband[0],/spline)],linestyle=1;,color=black
   oplot, [passband[1],passband[1]],[0,interpol(spc,lambda,passband[1],/spline)],linestyle=1;,color=black
endif


;angular coefficient and intercepta of the line connection the 2
;pseudocontinua values
m=(red_pseudocont-blu_pseudocont)/(red_lam-blu_lam)
q=blu_pseudocont-blu_lam*m[0]

; REBIN OF THE SPECTRUM IN THE REGION OF INTEREST IN ORDER TO KEEP THE
; CORRECT STARTING AND ENDING POINTS
lambda_new=findgen(10.*n_elements(indici_band))*(passband[1]-passband[0])/(10.*n_elements(indici_band)-1)+passband[0]

cont=m[0]*lambda_new+q[0]
spc_new_rebin=interpol(spc,lambda,lambda_new,/spline)
integranda=1.-spc_new_rebin/cont

;nel caso in cui ci siano punti dello spettro nulli non li
;considero nell'integrale.... un po' artificioso, non so se al voglio
;tenere come feature.
indici_not_null = where(spc_new_rebin gt 0.0,complement=indici_null__)

if n_elements(indici_not_null) le 10 then begin
   ew = 0./0.
   index_mag = 0./0.
   errors=[ew ,index_mag]
   goto, end_mdap_do_measure_index
endif

ew=int_tabulated(lambda_new[indici_not_null],integranda[indici_not_null])    ;equivalenth width
index_mag_=int_tabulated(lambda_new[indici_not_null],spc_new_rebin[indici_not_null]/cont[indici_not_null])
;index_mag=-2.5*alog10(1./(passband[1]-passband[0])*index_mag_)    ;line index in magnitude
if indici_null__[0] ne -1 then index_mag=-2.5*alog10(1./(max(lambda_new[indici_not_null])-min(lambda_new[indici_not_null]))*index_mag_)    ;line index in magnitude
if indici_null__[0] eq -1 then index_mag=-2.5*alog10(1./(passband[1]-passband[0])*index_mag_)    ;line index in magnitude




;;; ERROR COMPUTATION (See Cardiel et al. 1998, equations 41-46)
;
; spectrum and err need to be resampled in uniform ang/pxl for error
; computation. Input vectors WILL NOT be overwritten
lambda__ = mdap_range(min(lambda_),max(lambda_),n_elements(lambda)*2)
spc__ = interpol(spc_,lambda_,lambda__)
if n_elements(noise) ne 0 then noise__ = interpol(noise,lambda_,lambda__)

step_cardiel=lambda__[1]-lambda__[0]    ;step in angstrom/pixel 
Nb_indici = where(lambda__ ge blue_cont[0] and lambda__ le blue_cont[1])
Ni_indici = where(lambda__ ge passband[0] and lambda__ le passband[1])
Nr_indici = where(lambda__ ge red_cont[0] and lambda__ le red_cont[1])
N_indici=[Nb_indici,Ni_indici,Nr_indici]
lambda_cardiel=lambda__(N_indici)
Ntot_cardiel=n_elements(N_indici)

spectrum_cardiel=spc__(N_indici);*normalization[0]
noise_cardiel=sqrt(spectrum_cardiel)
if n_elements(noise) ne 0 then noise_cardiel=noise__(N_indici)

SN_cardiel=1./sqrt(step_cardiel)*total(spectrum_cardiel/noise_cardiel)/Ntot_cardiel   ;formula 42

;print, sn_cardiel
;stop
dl_b=blue_cont[1]-blue_cont[0]
dl_c=passband[1]-passband[0]
dl_r=red_cont[1]-red_cont[0];

l_b=mean(blue_cont)
l_c=mean(passband)
l_r=mean(red_cont)
q1=(l_r-l_c)/(l_r-l_b)
q2=(l_c-l_b)/(l_r-l_b)
c2= sqrt(1./dl_c+q1^2/dl_b+q2^2/dl_r); formula 44

c1=dl_c*c2 ;formula 43
error_angstrom=(c1-c2*ew)/SN_cardiel   ;formula 41
c3=2.5*c2*alog10(exp(1))  ;formula 46
error_mag=c3/SN_cardiel   ;formula 45
errors=[error_angstrom,error_mag]
;; 


if not keyword_set(noplot) then begin
   inddd=where(lambda ge blue_cont[0]-20 and lambda le red_cont(1)+20)
   mn=min(spc(inddd))
   mx=max(spc(inddd))
   rgn=(mx-mn)*1.05
   midplot=0.5*(mx+mn)
  ; xyouts,min(lambda),1.28*midplot,'Dispersion=  '+strcompress(string(lambda(1)-lambda(0)),/remove_all)+' (ang/pxl)',color=black
  ; xyouts,min(lambda),1.25*midplot,'Blue pscont= '+strcompress(string(blu_pseudocont),/remove_all),color=black
  ; xyouts,min(lambda),1.22*midplot,'Red pscont=  '+strcompress(string(red_pseudocont),/remove_all),color=black
  ; xyouts,min(lambda),1.19*midplot,'Eq. Width=       '+strcompress(string(ew),/remove_all)+' (ang)',color=black
  ; xyouts,min(lambda),1.16*midplot,'Index=       '+strcompress(string(index_mag),/remove_all)+' (mag)',color=black
;   xyouts,min(lambda(inddd))+1,mx+0.05*rgn,'Dispersion=  '+strcompress(string(lambda(1)-lambda(0)),/remove_all)+' (ang/pxl)',color=black
;   xyouts,min(lambda(inddd))+1,mx-0.00*rgn,'Blue pscont= '+strcompress(string(blu_pseudocont),/remove_all),color=black
;   xyouts,min(lambda(inddd))+1,mx-0.05*rgn,'Red pscont=  '+strcompress(string(red_pseudocont),/remove_all),color=black
;   xyouts,min(lambda(inddd))+1,mx-0.10*rgn,'Eq. Width=       '+strcompress(string(ew),/remove_all)+' (ang)',color=black
;   xyouts,min(lambda(inddd))+1,mx-0.15*rgn,'Index=       '+strcompress(string(index_mag),/remove_all)+' (mag)',color=black
   xyouts,min(lambda(inddd))+21,mx+0.05*rgn,'Dispersion=  '+mdap_round_str(lambda(1)-lambda(0),4)+' (ang/pxl)';,color=black
   xyouts,min(lambda(inddd))+21,mx-0.00*rgn,'Blue pscont= '+mdap_round_str(blu_pseudocont,4);,color=black
   xyouts,min(lambda(inddd))+21,mx-0.05*rgn,'Red pscont=  '+mdap_round_str(red_pseudocont,4);,color=black
   xyouts,min(lambda(inddd))+21,mx-0.10*rgn,'Eq. Width=       '+mdap_round_str(ew,4)+' +/- '+mdap_round_str(errors[0],4)+' (ang)';,color=black
   xyouts,min(lambda(inddd))+21,mx-0.15*rgn,'Index=       '+mdap_round_str(index_mag,5)+' +/- '+mdap_round_str(errors[1],5)+' (mag)';,color=black
;   xyouts,min(lambda(inddd))+21,mx-0.20*rgn,'S/N=       '+mdap_round_str(SN_cardiel,3),color=black
;   xyouts,min(lambda(inddd))+21,mx-0.25*rgn,'mean noise=       '+mdap_round_str(mean(sqrt(noise_cardiel)),3),color=black
endif
;device, /close
;set_plot ,'x'

;print, 'EW is: ',index
;print, 'Index is ',index_mag, +' (mag)'
;print, '...'
end_mdap_do_measure_index:
end
