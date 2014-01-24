pro mdap_display_2dmap_ps,map,x2d,y2d,filename=filename,extra_key=extra_key,invert_pl_range=invert_pl_range

if n_elements(extra_key) gt 0 then begin
   for i= 0, n_elements(extra_key)-1 do d=execute(extra_key[i])
endif


map_=map
indici =where(finite(map_) eq 0)
;map_[indici] = 1e10
set_plot,'ps'
basic_colors_rev, black, white, red, green, blue, yellow, cyan, magenta, orange, mint, purple, pink, olive, lightblue, gray   
   device,filename=filename,/color,bits_per_pixel=100,xsize=12,ysize=12
;window,1,retain=2,xsize=500,ysize=500,colors=255
;device, decomposed=0
;executed=0
;re_execute_for_ps:

if n_elements(value_range) eq 0 then begin
   value_range_ = minmax(map[where(finite(map) gt 0)])
   rang = value_range_[1]-value_range_[0]
   value_range = [min(map[where(finite(map) gt 0)])-0.05*rang , max(map[where(finite(map) gt 0)])+0.05*rang]
endif
;if keyword_set(invert_pl_range) then  value_range =  [value_range[1], value_range[0]] 
;stop
;map_=map
;indici =where(finite(map_) eq 0)
basic_colors_rev, black, white, red, green, blue, yellow, cyan, magenta, orange, mint, purple, pink, olive, lightblue, gray   

if keyword_set(invert_pl_range) then begin
   if indici[0] ne -1 then map_[indici] = value_range[1] +10
mdap_place_image,-map_,x2d,y2d,itable='mdap_color.tab',min=-value_range[1],max=-value_range[0],$
       ptable=1,title=title,xtitle=xtitle,ytitle=ytitle,position=[0.1,0.2,0.98,-.9]
endif else begin
   if indici[0] ne -1 then map_[indici] = value_range[0] -10

  mdap_place_image,map_,x2d,y2d,itable='mdap_color.tab',min=value_range[0],max=value_range[1],$
        ptable=1,title=title,xtitle=xtitle,ytitle=ytitle,position=[0.1,0.2,0.98,-.9]
endelse



ypos = min(y2d)-.15*(max(y2d)-min(y2d))  
xpos=min(x2d)
xlen=max(x2d)-min(x2d)
ylen=0.05*(max(y2d)-min(y2d))
if keyword_set(invert_pl_range) then value_range =  [value_range[1], value_range[0]] 
mdap_place_cbar,xpos,ypos,xlen,ylen,value_range,ctable='mdap_color.tab',label=label

   device, /close
   set_plot, 'x'

end
