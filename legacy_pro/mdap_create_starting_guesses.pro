pro mdap_create_starting_guesses,input_map,input_xbin,input_ybin,x2d,y2d,$
 output_xbin,output_ybin,output_start_guess,velocity_initial_guess,velocity_dispersion_initial_guess,H3_initial_guess,H4_initial_guess,h3h4=h3h4

; This interface is responsible to transform the kinematics measured on
; a spatial binning scheme into starting guesses for the fit of another
; spatial binning scheme. The main module that performs the computation
; is mdap\_interpolate\_2dmaps.pro
;
; INPUT
; input_map  Nx4 (or Nx2)  elements array containing the measurements  of velocity, velocity dispersion 
;            of the N spectra to be used to get the starting guesses. If the keyword /h3h4 is set, input_map
;            must be Nx4 elements array, and previous measurements of h3 and h4 must be given.
;              input_map[*,0]: previous velocity measurements
;              input_map[*,1]: previous velocity dispersion measurements
;              input_map[*,2]: (optional) previous H3 measurements
;              input_map[*,3]: (optional) previous H4 measurements
;
; input_xbin  N elements array. X coordinates of the center of the spatial bins where the measurements in input_map refer to.
;
; input_ybin  N elements array. Y coordinates of the center of the spatial bins where the measurements in input_map refer to.
;
; x2d         N'xM' elements array. two dimensional map of X coordinates of a regular grid mapping the full field of view. 
;
; y2d         N'xM' elements array. two dimensional map of Y coordinates of a regular grid mapping the full field of view.
; 
; output_xbin M elements array. X coordinates of the center of the spatial bins where the starting guesses are to be computed
;
; output_ybin M elements array. Y coordinates of the center of the spatial bins where the starting guesses are to be computed
;
; velocity_initial_guess  Float. Guess for the velocity to be used in the case of previous measurements are NaN, or not defined.
;
; velocity_dispersion_initial_guess  Float. Guess for the velocity to be used in the case of previous measurements are NaN, or not defined, 
;                                    or outside a fiducial interval 21<sigma<499 km/sec.
;
; H3_initial_guess  Float. Guess for the velocity to be used in the case of previous measurements are NaN, or not defined, or outside a fiducial
;                          interval -0.399 < h3 < 0.399.
;
; H4_initial_guess  Float. Guess for the velocity to be used in the case of previous measurements are NaN, or not defined, or outside a fiducial
;                          interval -0.399 < h4 < 0.399.. 
;
;
; OPTIONAL KEYWORDS
;
;/h3h4    If set, the guesses of Gauss Hermite moments (h3, h4) will be computed  (otherwise set to 
;        to H3\_initial\_guess and H4\_initial\_guess). The input input_map must contain h3 and h4 measurements.
;
;
; OUTPUTS
; output_start_guess  Mx4 elements array containing the interpolated velues of of velocity, velocity dispersion of the N spectra 
;             to be used as starting guesses. If the keyword /h3h4 is
;             set, it will be Mx4 elements array.
;               output_start_guess[*,0]: velocity interpolated on the  output_xbin, output_ybin, grid to be used as starting guesses.
;               output_start_guess[*,1]: velocity dispersion interpolated on the  output_xbin, output_ybin, grid to be used as starting guesses.
;               output_start_guess[*,2]: (optional) H3 interpolated on the  output_xbin, output_ybin, grid to be used as starting guesses.
;               output_start_guess[*,3]: (optional) H4 interpolated on the  output_xbin, output_ybin, grid to be used as starting guesses.
;
;


if n_elements(input_xbin) gt 2 then begin
   mdap_interpolate_2dmaps,input_map[*,0],input_xbin,input_ybin,x2d,y2d,output_xbin,output_ybin,vel_starting_guesses
   mdap_interpolate_2dmaps,input_map[*,1],input_xbin,input_ybin,x2d,y2d,output_xbin,output_ybin,disp_starting_guesses
   if keyword_set(h3h4) then mdap_interpolate_2dmaps,input_map[*,2],input_xbin,input_ybin,x2d,y2d,output_xbin,output_ybin,H3_starting_guesses
   if keyword_set(h3h4) then mdap_interpolate_2dmaps,input_map[*,3],input_xbin,input_ybin,x2d,y2d,output_xbin,output_ybin,H4_starting_guesses
endif else begin
   vel_starting_guesses = velocity_initial_guess[0]
   disp_starting_guesses = velocity_dispersion_initial_guess[0]
   if keyword_set(h3h4) then H3_starting_guesses = H3_initial_guess[0]
   if keyword_set(h3h4) then H4_starting_guesses = H4_initial_guess[0]
endelse

if keyword_set(h3h4) then output_start_guess = fltarr(n_elements(output_xbin),4)
if not keyword_set(h3h4) then output_start_guess = fltarr(n_elements(output_xbin),2)

output_start_guess[*,0]=vel_starting_guesses
output_start_guess[*,1]=disp_starting_guesses
if keyword_set(h3h4) then output_start_guess[*,2]=H3_starting_guesses
if keyword_set(h3h4) then output_start_guess[*,3]=H4_starting_guesses

; I make sure that there are no NULL values
indici = where(finite(output_start_guess[*,0]) eq 0)
if indici[0] ne -1 then output_start_guess[indici,0] = velocity_initial_guess[0]
indici = where(finite(output_start_guess[*,1]) eq 0)
if indici[0] ne -1 then output_start_guess[indici,1] = velocity_dispersion_initial_guess[0]


; I constrain 21. < sigma < 499.  km/sec
indici = where(output_start_guess[*,1] le 21.)
if indici[0] ne -1 then output_start_guess[indici,1] = 21.
indici = where(output_start_guess[*,1] ge 499.)
if indici[0] ne -1 then output_start_guess[indici,1] = 499.

if keyword_set(h3h4) then begin
   ; I constrain -0.399 < H3 < 0.399
   indici = where(output_start_guess[*,2] le -.399 or finite(output_start_guess[*,2]) eq 0)
   if indici[0] ne -1 then output_start_guess[indici,2] = -.399
   indici = where(output_start_guess[*,2] ge .399 or finite(output_start_guess[*,2]) eq 0)
   if indici[0] ne -1 then output_start_guess[indici,2] = .399

   ; I constrain -0.399 < H4 < 0.399
   indici = where(output_start_guess[*,3] le -.399 or finite(output_start_guess[*,3]) eq 0)
   if indici[0] ne -1 then output_start_guess[indici,3] = -.399
   indici = where(output_start_guess[*,3] ge .399 or finite(output_start_guess[*,3]) eq 0)
   if indici[0] ne -1 then output_start_guess[indici,3] = .399
endif


end
