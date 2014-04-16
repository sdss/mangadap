pro mdap_do_k_rprofiles,radii,vel_,vel_err,sigma,sigma_err,xbin,ybin,flux,ell,pa, $
lambda_profile,lambda_profile_err,vsigma_profile,vsigma_profile_err,sigma_profile,sigma_profile_err,add_last=add_last
 

; INPUTS
; vel_      [N elm. array] Velocities of the N input spectra.
;
; vel_err   [N elm. array] Velocity errors of the N input spectra. Systemic velocity needs to be already subtracted.
;
; sigma     [N elm. array] Velocity dispersions of the N input spectra.
;
; sigma_err [N elm. array] Velocity dispersion errors of the N input spectra.
;
; xbin      [N elm. array] X coordinate of the N spectra.
;
; ybin      [N elm. array] Y coordinate of the N spectra.
;
; flux      [N elm. array] Flux of the N spectra.
;
; ell       [float]        Mean galaxy ellipticity.
; 
; pa        [float]        Mean galaxy position angle (0=north, 90=east)
;
;
; radii     [M elem. array] upper boundaries at which radial profiles need
;           to be computed. If keyword add_last is set, the outermost
;           bin is appended to radii, and it in output it contains M+1 elements.
;
; KEYWORDS
;
; \add_last  If set, a last bin is added to the profiles, which
;            includes all the points outside the last radii value.

; OUTPUTS 
;
; lambda_profile  [M elem. array] radial (cumulative) profile of the mean lambda parameter. If \add_last is set, it contains M+1 elements.
;
; lambda_profile_err  [M elem. array] errors on lambda_profile. If \add_last is set, it contains M elements.
;
;
; vsigma_profile  [M elem. array] radial (cumulative) profile of the mean V/Sigma parameter. If \add_last is set, it contains M+1 elements.
;
;
; vsigma_profile_err   [M elem. array] errors on vsigma_profile . If \add_last is set, it contains M+1 elements.
;
;
; sigma_profile   [M elem. array] radial (cumulative) profile of the mean Sigma parameter. If \add_last is set, it contains M+1 elements.
;
;
; sigma_profile_err   [M elem. array] errors on sigma_profile. If \add_last is set, it contains M+1 elements.
;



;-- rotate coordianates to align galaxy photometric major axis to
;   vertical direction
xbin_r =  xbin * cos(pa*!pi/180.) + ybin*sin(pa*!pi/180.)  ;clockwise rotation of PA
ybin_r = -xbin * sin(pa*!pi/180.) + ybin*cos(pa*!pi/180.)
;--

;-- computation of semimajor axis (aligned vertically, with ellipticity
;   as input value
abin = sqrt(xbin_r^2/(1.-ell)^2 + ybin_r^2)
;--

if keyword_set(add_last) then begin
  if max(radii) lt max(abin) then radii=[radii,max(abin)]
endif
nrad=n_elements(radii)

lambda_profile_err = fltarr(nrad)-999
vsigma_profile_err = fltarr(nrad)-999
sigma_profile_err = fltarr(nrad)-999
lambda_profile = fltarr(nrad)-999
vsigma_profile = fltarr(nrad)-999
sigma_profile = fltarr(nrad)-999


vel=abs(vel_)
FOR i = 1 ,nrad-1 DO BEGIN

   sel = where(abin lt radii[i])
   if sel[0] eq -1 then continue

   F = flux[sel]
   R = abin[sel]
   V=vel[sel]
   dV = vel_err[sel]
   S=sigma[sel]
   dS=sigma_err[sel]

   lambda_profile[i] = total(F * R * V) / total(F * R * sqrt(V^2+S^2))
   vsigma_profile[i] = total(F * V/S) / total(F)
   sigma_profile[i] = total(F * S) / total(F)

   W = total(R*F*sqrt(V^2+S^2))
   lambda_profile_err[i] = sqrt(total((F*R/W*dV)^2))
 
   T = total(F)
   vsigma_profile_err[i] = sqrt(total((F*dV/S/T)^2 + (F*V*dS/S^2/T)^2))

   sigma_profile_err[i] = sqrt(total((F*dS/T)^2))

ENDFOR

end
