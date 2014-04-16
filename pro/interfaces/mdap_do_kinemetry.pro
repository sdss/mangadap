pro mdap_do_kinemetry,image,x2d,y2d,x,y,velocity,velocity_err,$ 
                      PA_kin,PA_kin_std, q_kin, q_kin_std, Vsyst, Vsyst_std, Rad_kin, Vrot, Vrot_err,Vexp, Vexp_err,$
                      ;gal_center_x,gal_center_y
                      version=version

; This module implements the kinemetry.pro by D. Kraijnovic to extract the rotation curve and outflow/inflows. 
; The following steps are executed: 
;
; 1) The galaxy image is used to determine the galaxy center. If
;    the center is outside a $3''\times 3''$ box aroud the center of
;    the field of view, the galaxy center is automatically set to
;    (0,0), i.e. the center of the field of view. The galaxy kinematic
;    center is set to be the galaxy photometric center in the
;    kinemetry.pro runs.
;
; 2) A first kinemetry run is executed, to get the The galaxy
;    systemic velocity is kept constant at all radii (a first
;    kinemetry.pro run is executed to measure the systemic velocities
;    Vs_i at each i-th radius. The galaxy systemic velocity (Vsyst) and
;    its 1sigma error (Vsyst_std) are computed by the median and the
;    standard deviation of the systemic velocities Vs_i).
;
; 3) A second kinemetry run is executed, fixing the kinematic center
;    and the systemic velocity. In this run, the kinematic position
;    angles (PA_i) and flattening (q_i) are measured for each
;    radius. The galaxy mean kinematic position angle (PA_kin) and its
;    error (PA_kin_std) are measured as the median and standard
;    deviation of all the PA_i previously measured.  The galaxy mean
;    kinematic axial ratio (q_kin) and its error (q_kin_std) are
;    measured as the median and standard deviation of all the q_i
;    previously measured.
;
; 4) A third and final kinemetry run is executed, fixing the
;    kinematic center, position angle, and axial ratio determined in the
;    previous runs. This last run determines the rotation velocity, and
;    expansion velocity (i.e. inflows/outflows) for several radii.
;
;
;INPUTS
; image         [NxM array] Galaxy image. It is used to determine the location of the center
;
; x2d           [NxM array] X coordinates corresponding to image (0,0) is the center of the field of view
;
; y2d           [NxM array] Y coordinates corresponding to image (0,0) is the center of the field of view.
;
; x             [T elements array] X coordinates at which the velocites are measured.
;
; y             [T elements array] X coordinates at which the velocites are measured.
;
; velocity      [T elements array] The measured velocities (in km/sec).
;
; velocity_err  [T elements array] The measured velocity errors (in km/sec).

;OUTPUTS
; PA_kin        [Float]   Median kinematic position angle used to determine Vrot and Vexp (third kinemetry.pro run). 
;
; PA_kin_std    [Float]   Standard deviation of the kinematic position angles measured in the second kinemetry run.
;
; q_kin         [Float]   Median kinematic axial ratio used to determine Vrot and Vexp (third kinemetry.pro run). 
;
; q_kin_std     [Float]   Standard deviation of the kinematic axial ratio measured in the second kinemetry run. 
;
; Vsyst         [Float]   Systemic velocity used to determine Vrot and Vexp (third kinemetry.pro run) 
;
; Vsyst_std     [Float]   Standard deviation of the systemic velocities determined in the first kinemetry run.
;
; Rad_kin       [W elements array]  Semi major axis of the ellipses where Vrot and Vexp are measured. 
;
; Vrot          [W elements array]  Rotational velocity measured at Rad_kin.  
;
; Vrot_err      [W elements array] Errors on Vrot. 
;
; Vexp          [W elements array] Outflow/Inflow velocity measured at Rad_kin.  
;
; Vexp_err      [W elements array] Errors on Vexp. 
;
; gal_center_x  [Float] X coordinate of the galaxy center (0,0 is the
; center of the field of view). OBSOLETE
;
; gal_center_y  [Float] Y coordinate of the galaxy center (0,0 is the
; center of the field of view). OBSOLETE
;
;*** OPTIONAL OUTPUTS ***
;
; version  [string]            Module version. If requested, the module is not execute and only version flag is returned

; Version 0.1: 21 Feb L. Coccato


version_module='0.7'
if n_elements(version) ne 0 then begin
 version = version_module
 goto, end_kinemetry
endif


;-------------- removed in version 0.7 -----------
; get image center from signal2d_reconstructed
; sz=size(image)
; gcntrd,image,sz[1]/2.,sz[2]/2.,xcen,ycen,3.
; xcen_pix=xcen
; ycen_pix=ycen
; x0=bilinear(x2d,xcen_pix,ycen_pix) 
; xcen_pix=xcen
; ycen_pix=ycen
; y0=bilinear(y2d,xcen_pix,ycen_pix)
; if abs(x0) ge 4 then x0 = 0.
; if abs(y0) ge 4 then y0 = 0.
;-------------------------------------------------
x0 = 0.
y0 = 0.
gal_center_x = x0
gal_center_y = y0

bad_data=where(finite(velocity) ne 1 or finite(velocity_err) ne 1 or velocity_err le 0 or velocity_err ge 10000., complement=good)
;first round to get Vsystemic
if good[0] eq -1 then begin ;no good points to run kinemetry
   goto, failure_mode
endif
;stop
d = execute('MDAP_KINEMETRY,x[good]-x0,y[good]-y0,velocity[good],rad,pa,q,cf,ERROR=velocity_err[good],NPA=9,NQ=9,NTRM=2,scale=1')
if d eq 0 then goto, failure_mode
   vsyst = median(cf[*,0],/even)
   vsyst_std = robust_sigma(cf[*,0])

;second round performing the analysis keeping the systemic velocity at 0
d = execute('MDAP_KINEMETRY,x[good]-x0,y[good]-y0,velocity[good]-vsyst,rad,pa,q,cf,ERROR=velocity_err[good], NPA=15,NQ=15,NTRM=2, ER_PA=ER_PA,ER_Q=ER_Q,ER_CF=ER_CF,COVER=0.2,/VSYS,scale=1')
if d eq 0 then goto, failure_mode



; fixing the geometry and final run.
indici = where(pa ge 180)
if indici[0] ne -1 then begin
   pa[indici] = pa[indici]-180.
   cf[indici,2] = -cf[indici,2]
   cf[indici,1] = -cf[indici,1]
endif

vm = mean(cf[*,2])
if vm le 0 then cf[*,2] = -cf[*,2]


  
pa_var=temporary(pa)
q_var=temporary(q)
PA_kin = median(pa_var,/EVEN)
q_kin = median(q_var,/even)
PA_kin_std = robust_sigma(pa_var)
q_kin_std = robust_sigma(q_var)

d = execute('MDAP_KINEMETRY,x[good]-x0,y[good]-y0,velocity[good]-vsyst,rad_f,pa_f,q_f,cf_f,ERROR=velocity_err[good],NTRM=2, ER_PA=ER_PA_f,ER_Q=ER_Q_f,ER_CF=ER_CF_f,COVER=0.2,/VSYS,PAQ=[PA_kin-90,q_kin],scale=1')
  
if d eq 0 then begin
   dd=execute('MDAP_KINEMETRY,x[good]-x0,y[good]-y0,velocity[good]-vsyst,rad_f,pa_f,q_f,cf_f,ERROR=velocity_err[good],NTRM=2, ER_PA=ER_PA_f,ER_Q=ER_Q_f,ER_CF=ER_CF_f,COVER=0.3,/VSYS,PAQ=[PA_kin-90,q_kin],scale=1')
   if dd eq 0 then goto, failure_mode
endif

   indici = where(pa_f ge 180) 
   if indici[0] ne -1 then begin
      pa_f[indici] = pa_f[indici]-180.
      cf_f[indici,2] = -cf_f[indici,2]
   endif

vm = mean(cf_f[*,2])
if vm le 0 then cf_f[*,2] = -cf_f[*,2]
if vm le 0 then PA_kin=PA_kin-180


Rad_kin=rad_f
Vrot = cf_f[*,2]
Vrot_err = er_cf_f[*,2]
Vexp = cf_f[*,1]
Vexp_err = er_cf_f[*,1]

goto, end_kinemetry
failure_mode:
PA_kin = [0]
PA_kin_std =[99]
q_kin = [0]
q_kin_std =[99]
Vsyst = [0]
Vsyst_std =[99]
Rad_kin = [0]
Vrot = [0]
Vrot_err =[99]
Vexp = [0]
Vexp_err =[99]
gal_center_x  = [0]
gal_center_y  = [0]



end_kinemetry:
end
