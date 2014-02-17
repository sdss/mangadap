pro mdap_do_kinemetry,image,x2d,y2d,x,y,velocity,velocity_err,$ 
                      PA_kin,PA_kin_std, q_kin, q_kin_std, Vsyst, Vsyst_std, Rad_kin, Vrot, Vrot_err,Vexp, Vexp_err,$
                      gal_center_x,gal_center_y

  
; This module implements the kinemetry.pro by D. Kraijnovic to extract the rotation curve and outflow/inflows. 
; The following criteria are applied: 
;
;  1) Center of rotation is fixed to be consistent with the location of the peack of luminosity (computed through centroid on an input image).    
;   
;  2) The galaxy systemic velocity is kept constant at all radii (a first kinemetry.pro run is executed to compute the average systemic velocity).
;
;  3) The kinematic position angle and kinematic axial ratio are kept constant: a first kinemetry run is done leaving PA and q free to vary, 
;     then a second kinemetry.pro run is performed fixing PA and q to the median of the values computed in the previous run.
;
; Only 2 terms of the kinematic expansion are fitted, which correspond
; to the rotation and outflows/inflows terms of the velocity field.


;INPUTS
; image   NxN array. Galaxy image. It is used to determine the location of the center
; x2d     NxN array. X coordinates corresponding to image (0,0 is the center of the field of view
; y2d
; x
; y
; velocity
; velocity_err

;OUTPUTS
; PA_kin
; PA_kin_std
; q_kin
; q_kin_std
; Vsyst
; Vsyst_std
; Rad_kin
; Vrot
; Vrot_err
; Vexp
; Vexp_err
; gal_center_x 
; gal_center_y 




; get image center from signal2d_reconstructed
sz=size(image)
gcntrd,image,sz[1]/2.,sz[2]/2.,xcen,ycen,3.
xcen_pix=xcen
ycen_pix=ycen
x0=bilinear(x2d,xcen_pix,ycen_pix) 
xcen_pix=xcen
ycen_pix=ycen
y0=bilinear(y2d,xcen_pix,ycen_pix)

gal_center_x = x0
gal_center_y = y0

bad_data=where(finite(velocity) ne 1 or finite(velocity_err) ne 1 or velocity_err le 0, complement=good)
;first round to get Vsystemic
if good[0] eq -1 then begin ;no good points to run kinemetry
   stop
endif

MDAP_KINEMETRY,x[good]-x0,y[good]-y0,velocity[good],rad,pa,q,cf,ERROR=velocity_err[good],NPA=9,NQ=9,NTRM=2
   vsyst = median(cf[*,0],/even)
   vsyst_std = robust_sigma(cf[*,0])

;second round performing the analysis keeping the systemic velocity at 0
   MDAP_KINEMETRY,x[good]-x0,y[good]-y0,velocity[good]-vsyst,rad,pa,q,cf,ERROR=velocity_err[good],$
       NPA=15,NQ=15,NTRM=2, ER_PA=ER_PA,ER_Q=ER_Q,ER_CF=ER_CF,COVER=0.2,/VSYS



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

d = execute('MDAP_KINEMETRY,x[good]-x0,y[good]-y0,velocity[good]-vsyst,rad_f,pa_f,q_f,cf_f,ERROR=velocity_err[good],NTRM=2, ER_PA=ER_PA_f,ER_Q=ER_Q_f,ER_CF=ER_CF_f,COVER=0.2,/VSYS,PAQ=[PA_kin-90,q_kin]')
  
if d eq 0 then MDAP_KINEMETRY,x[good]-x0,y[good]-y0,velocity[good]-vsyst,rad_f,pa_f,q_f,cf_f,ERROR=velocity_err[good],NTRM=2, ER_PA=ER_PA_f,ER_Q=ER_Q_f,ER_CF=ER_CF_f,COVER=0.3,/VSYS,PAQ=[PA_kin,q_kin]

   indici = where(pa_f ge 180) 
   if indici[0] ne -1 then begin
      pa_f[indici] = pa_f[indici]-180.
      cf_f[indici,2] = -cf_f[indici,2]
   endif

vm = mean(cf_f[*,2])
if vm le 0 then cf_f[*,2] = -cf_f[*,2]


Rad_kin=rad_f
Vrot = cf_f[*,2]
Vrot_err = er_cf_f[*,2]
Vexp = cf_f[*,1]
Vexp_err = er_cf_f[*,1]

end
