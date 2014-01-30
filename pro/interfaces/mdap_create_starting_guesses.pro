pro mdap_create_starting_guesses,input_map,input_xbin,input_ybin,x2d,y2d,$
 output_xbin,output_ybin,output_start_guess,velocity_initial_guess,velocity_dispersion_initial_guess,H3_initial_guess,H4_initial_guess,h3h4=h3h4


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
