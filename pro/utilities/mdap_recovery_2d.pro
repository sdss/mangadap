pro mdap_recovery_2d,input,binning,field,map2d,voronoi_binning=voronoi_binning,list=list,fits_file=fits_file

; input     name of the output manga_dap fits file
;
; binning   a flag that specifies the binning scheme we want to
;           extract the info from. Depending on the binning scheme,
;           only some quantities are available (e.g. the stellar
;           kinematics is available only on the str binning scheme)
;                1 = tpl binning scheme
;                2 = str binning scheme
;                3 = ems binning scheme
;                4 = radial binning scheme
;
; field     the quantity we want to map. e.g. VEL, VEL_ERR, CA4227.
;           The list of quantities is available with the keyword list.
;
;map2d the reconstructed 2d map.
;
;fits_file   the name of the output fits file which will contain the 2D
;           reconstructed map
;
;voronoi_binning   if set, the geometry of the voronoi binningh scheme
;                  is mantained. In the case of radial binning, this is automatically
;                  set,and the radial binning geometry is mantained.
;                  The default is to use for each (x,y) the value of the
;                  neighbours. The closest points are 
;                    1- the closest point to x,y. This defines the smallest distance
;                    2- all the points closer to smallest_distance+sqrt(stepx^2+stepy^2)/2. 
;                       (this is to smooth the borders, and to merge
;                       RSS spectra which are closer to the spatial resolution.
;
; list    if set, it prompts out the list of available fields for the selected binning


;input='test_datacube.fits'
;field='hb'
mdap_readfits_2d,input,image,header=hrd,x2d=x2d,y2d=y2d
junk = readfits(input,hdr,exten_no=0)
if binning eq 1 then begin
   binning_map1 = readfits(input,exten_no=2)
   data1=mrdfits(input,3)
 ;  d=execute('quantity = data1.'+field)
endif

if binning eq 2 then begin
   binning_map1 = readfits(input,exten_no=4)
   data1=mrdfits(input,5)
 ;  d=execute('quantity = data1.'+field)
endif

if binning eq 3 then begin
   binning_map1 = readfits(input,exten_no=6)
   data1=mrdfits(input,7)
  ; d=execute('quantity = data1.'+field)
endif

if binning eq 4 then begin
   binning_map1 = readfits(input,exten_no=8)
   data1=mrdfits(input,9)
;   d=execute('quantity = data1.'+field)
endif
   
uniq_x = x2d(uniq(x2d,sort(x2d)))
stepx = uniq_x[1]-uniq_x[0]
uniq_y = y2d(uniq(y2d,sort(y2d)))
stepy = uniq_y[1]-uniq_y[0]
min_dist=sqrt(stepx^2+stepy^2)/2.

if keyword_set(list) then begin
   help, data1,/str
   print, 'field to plot?'
   read,field
end
d=execute('quantity = data1.'+field)
if binning eq 4 then voronoi_binning = 1
sz=size(image)
map2d = fltarr(sz[1],sz[2])/0.
if ~keyword_set(voronoi_binning) then begin ; the field at position (x,y) is equal
   for j = 0, sz[2] -1 do begin             ; to the mean of the N closest measurements.
      for i = 0, sz[1]-1 do begin           ; They are the closest point + everything within +/-min_dist
         dist = sqrt( (data1.x-x2d[i,j])^2 + (data1.y-y2d[i,j])^2 )
         indici = where(dist le min(dist)+min_dist)
         if indici[0] ne -1 then map2d[i,j] =mean(quantity[indici])
      endfor
   endfor
   indici = where(binning_map1 eq -1)
   if indici[0] ne -1 then map2d[indici] = 0./.0
endif else begin                            ; I mantain the voronoi binning geometry. Use it only for datacubes, as it is not precise for RSS spectra.
   for i = 0, max(binning_map1) do begin
      indici = where(binning_map1 eq i)
     if indici[0] ne -1 then  map2d[indici] = quantity[i]
   endfor
endelse

if n_elements(fits_file) ne 0 then writefits,fits_file,map2d,hdr
if n_elements(voronoi_binning) ne 0 then junk = temporary(voronoi_binning)
;stop
end
