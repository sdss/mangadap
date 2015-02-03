;pro esempio;
;
;set_plot, 'ps'
;device, filename='idl.ps',/color,bits_per_pixel=100
;plot,[0,2],[0,1],position=[0.1,0.1,0.8,0.7],title='!3test'
;place_cbar,0.1,0.9,1,.1,[-100.,100.],ctable='mycolor.tab',label='Values [colors]',csize=.5;,marks=[-100,-50,-10,0,10,50,100]
;device, /close
;set_plot, 'x'
;end
;
;pro esempio2
; 
;plot,[0,2],[0,1],position=[0.1,0.1,0.8,0.7],title='!3test',xrange=[0,2],yrange=[0,1]
;place_cbar,0.1,1.2,1,.1,[-210.,210.],ctable=18,label='Values [units]',marks=[-200.,-100.,0.,100.,200.];
;
;end


pro mdap_place_cbar,x,y,xlen,ylen,zrange,ctable=ctable,label=label,marks=marks,csize=csize,invert_ctable=invert_ctable

basic_colors, black, white, red, green, blue, yellow, cyan, magenta, $
orange, mint, purple, pink, olive, lightblue, gray   


;x,y         starting X,Y positions for the color bar (in DATA units)
;xlen,ylen   lenght of the color bar (in DATA units)
;zrange      [vmin,vma] minimun and maximum value to plot.
;

; OPTIONAL INPUTS
; ctable=ctable     if [string]: name of the colortable to use (e.g. mycolor.tab)
;                 if [float] or [integer]: number of the color table to load
;                with loadct . Default: ctable = 0         
;   
;label=label       label to associate to the colorbar. Default = no label
;
; marks=marks   marks to put in the colorbar (default=automatic)
;
;csize=csize     size for the labels and numbers of the colorbar (default = .9)
; invert_ctable  if set, the color scheme is inverted
; Version: 2.0  bug on negative coordinates solved
; Version: 2.1  xstyle=1 missing

if n_elements(marks) eq 0 then marks = 0
if n_elements(csize) eq 0 then csize=.9
;get sizes from the current device
dimY=!D.Y_size 
dimX=!D.X_size 

ctable_=0
if n_elements(ctable) ne 0 then ctable_=ctable
sz=size(ctable_)
if sz[1] eq 7 then loadct,41,file=ctable_,/silent
if sz[1] eq 2 or sz[1] eq 4 then loadct,ctable_,/silent

label_=' '
if n_elements(label) ne 0 then label_=label

npt=1000.
colorbar=fltarr(npt,10)
vmax=zrange[1]
vmin=zrange[0]
value_vector=findgen(npt)*abs((vmax-vmin))/npt+vmin
if keyword_set(invert_ctable) then value_vector=reverse(value_vector);findgen(npt)*abs((-vmax+vmin))/npt-vmax
for i = 0, 9 do begin
    colorbar(*,i)=value_vector
endfor

if !d.name eq 'X' then begin
   dimY=!D.Y_size 
   dimX=!D.X_size 
   frx=xlen/(!x.crange[1]-!x.crange[0])
   fry=ylen/(!y.crange[1]-!y.crange[0])
   xsize=dimX*(!x.window[1]-!x.window[0])*frx
   ysize=dimY*(!y.window[1]-!y.window[0])*fry
   img=congrid(colorbar,xsize+1,ysize+1)
   tv,bytscl(img,/nan),x,y,/DATA 

  ; x0=!x.window[0]+(!x.window[1]-!x.window[0])*abs(x)/(!x.crange[1]-!x.crange[0])
  ; y0=!y.window[0]+(!y.window[1]-!y.window[0])*abs(y)/(!y.crange[1]-!y.crange[0])
   x0=!x.window[0]+(!x.window[1]-!x.window[0])*(x-!x.crange[0])/(!x.crange[1]-!x.crange[0])
   y0=!y.window[0]+(!y.window[1]-!y.window[0])*(y-!y.crange[0])/(!y.crange[1]-!y.crange[0])

   x1=x0+(!x.window[1]-!x.window[0])*(xlen)/(!x.crange[1]-!x.crange[0])
   y1=y0+(!y.window[1]-!y.window[0])*(ylen)/(!y.crange[1]-!y.crange[0])
   loadct,0,/silent
;stop
   plot, [0,1],[0,1],/noerase,xrange=[vmin,vmax],yrange=[0,1],xtitle=label_,/nodata,$
     position=[x0,y0,x1,y1],/norm,ytickname=replicate(' ',20),yticks=1,ytickv=[0,1],$
     xtickv=marks,xticks=(n_elements(marks)-1),xticklen=1.05,color=255,charsize=csize,xstyle=1


endif
if !d.name eq 'PS' then begin
   tv,bytscl(colorbar),x,y,xsize=xlen,ysize=ylen,/data
   
   x0=!x.window[0]+(!x.window[1]-!x.window[0])*(x-!x.crange[0])/(!x.crange[1]-!x.crange[0])
   y0=!y.window[0]+(!y.window[1]-!y.window[0])*(y-!y.crange[0])/(!y.crange[1]-!y.crange[0])

  ; x1=!x.window[0]+(!x.window[1]-!x.window[0])*(xlen+x)/(!x.crange[1]-!x.crange[0])
  ; y1=!y.window[0]+(!y.window[1]-!y.window[0])*(ylen+y)/(!y.crange[1]-!y.crange[0])
   x1=x0+(!x.window[1]-!x.window[0])*(xlen)/(!x.crange[1]-!x.crange[0])
   y1=y0+(!y.window[1]-!y.window[0])*(ylen)/(!y.crange[1]-!y.crange[0])
   loadct,0,/silent
   plot, [0,1],[0,1],/noerase,xrange=[vmin,vmax],yrange=[0,1],xtitle=label_,/nodata,$
     position=[x0,y0,x1,y1],/norm,ytickname=replicate(' ',20),yticks=1,ytickv=[0,1],$
     xtickv=marks,xticks=(n_elements(marks)-1),xticklen=1.05,charsize=csize,xstyle=1


endif

;print, xsize
;stop
end
