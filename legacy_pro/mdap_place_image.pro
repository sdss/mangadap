pro mdap_place_image,ima2d,x2d,y2d,fits_file=fits_file,itable=itable,ptable=ptable,$
                extra_inputs=extra_inputs,$        
  title=title,xtitle=xtitle,ytitle=ytitle,min=min,max=max,position=position

; INPUTS
; ima2d 2D array       2D array to display
;
; x2d                  2D array containing the X coordinates to
;                      display. same number of elements of ima2d
;
; y2d                  2D array containing the Y coordinates to
;                      display. same number of elements of ima2d
;

; OPTIONAL INPUTS
; fits_file  [string]  if provided, the ima2d, x2d and y2d are read
;                      from fits_file
;
;
; extra_inputs [string] vector of optional inputs to execute. Only
;                      those matching with plot/bytscl are used.
;
; itable=itable     Color table for the image to display 
;                   if [string]: name of the colortable to use (e.g. mdap_color.tab)
;                   if [float] or [integer]: number of the color table to load
;                   with loadct . Default: ctable = 0         
;
; ptable=ptable     Color table for the plot around the image. Only
;                   IDL standard tables are supported. Default = 0
;
; position          same of IDL plot. If position[2] or position[3]
;                   are negative, they will be recomputed in order to preserve the
;                   object scale
;
; History
; Version 1.0 20 Dec 2011  L. Coccato
; Version 2.0 27 Jun 2012  L. Coccato. Aspect ratio is maintained if
;                           position[2] or position[3] are negative.
; Version 2.1 27 Jun 2012  L. Coccato. Graphic keywords added, no need
;                          to include them in extra_inputs (which is 
;                          manteined for compatibility)
;                  
if n_elements(fits_file) ne 0 then readfits_2d,fits_file, ima2d,x2d=x2d,y2d=y2d

if n_elements(extra_inputs) ne 0 then begin
   for i = 0, n_elements(extra_inputs)-1 do begin
      d=execute(extra_inputs[i])
   endfor
endif

if n_elements(position) ne 0 then begin
   if position[2] le 0 or position[3] le 0 then begin ;KEEP STARTING POSITION, BUT RECOMPUTING END POSITION TO KEEP ASPECT RATIO
      
      extr =.99-max(position[0:1])
      if position[2] gt 0 then extr = position[2]-max(position[0:1])
      if position[3] gt 0 then extr = position[3]-max(position[0:1])

      maj1 = max([max(x2d)-min(x2d),max(y2d)-min(y2d)])
      position[2] = position[0]+ (max(x2d)-min(x2d)) / maj1 * extr
      position[3] = position[1]+ (max(y2d)-min(y2d)) / maj1 * extr
   endif

endif


if !d.name eq 'PS' then begin


   ptable_=0
   if n_elements(ptable) ne 0 then ptable_=ptable
   loadct,ptable_,/silent
   plot, [0,1],[0,1],xrange=[min(x2d),max(x2d)],yrange=[min(y2d),max(y2d)],/norm,$
         position=position,xstyle=1,ystyle=1,noerase=noerase,$
         ytickname=replicate(' ',30),xtickname=replicate(' ',30)


   xlen=max(x2d)-min(x2d)
   ylen=max(y2d)-min(y2d)
   x=min(x2d)
   y=min(y2d)


   itable_=0
   if n_elements(itable) ne 0 then itable_=itable
   sz=size(itable_)
   if sz[1] eq 7 then loadct,41,file=itable_,/silent
   if sz[1] eq 2 or sz[1] eq 4 then loadct,itable_,/silent

   tv,bytscl(ima2d,/nan,min=min,max=max),x,y,xsize=xlen,ysize=ylen,/data

   x0=(!x.window[1]-!x.window[0])*(x-!x.crange[0])/(!x.crange[1]-!x.crange[0])
   y0=(!y.window[1]-!y.window[0])*(y-!y.crange[0])/(!y.crange[1]-!y.crange[0])

   x1=x0+(!x.window[1]-!x.window[0])*(xlen)/(!x.crange[1]-!x.crange[0])
   y1=y0+(!y.window[1]-!y.window[0])*(ylen)/(!y.crange[1]-!y.crange[0])

   loadct,ptable_,/silent

   plot, [0,1],[0,1],/noerase,xrange=[min(x2d),max(x2d)],yrange=[min(y2d),max(y2d)],$
         position=position,/nodata,xstyle=1,ystyle=1,/norm,title=title,xtitle=xtitle,ytitle=ytitle,$
         ytickname=ytickname,xtickname=xtickname,charsize=charsize
endif

if !d.name eq 'X' then begin

   ptable_=0
   if n_elements(ptable) ne 0 then ptable_=ptable
   loadct,ptable_,/silent
   plot, [0,1],[0,1],xrange=[min(x2d),max(x2d)],yrange=[min(y2d),max(y2d)],/norm,$
         position=position,xstyle=1,ystyle=1,color=0,$
         ytickname=replicate(' ',30),xtickname=replicate(' ',30),noerase=noerase,/nodata
   


   ptable_=0
   if n_elements(ptable) ne 0 then ptable_=ptable
   loadct,ptable_,/silent

   xlen=max(x2d)-min(x2d)
   ylen=max(y2d)-min(y2d)
   x=min(x2d)
   y=min(y2d)

   dimY=!D.Y_size 
   dimX=!D.X_size 
   frx=xlen/(!x.crange[1]-!x.crange[0])
   fry=ylen/(!y.crange[1]-!y.crange[0])
   xsize=dimX*(!x.window[1]-!x.window[0])*frx
   ysize=dimY*(!y.window[1]-!y.window[0])*fry
   ;stop
   img=congrid(ima2d,xsize+1,ysize+1)

   itable_=0
   if n_elements(itable) ne 0 then itable_=itable
   sz=size(itable_)
   if sz[1] eq 7 then loadct,41,file=itable_,/silent
   if sz[1] eq 2 or sz[1] eq 4 then loadct,itable_,/silent

   tv,bytscl(img,/nan,min=min,max=max),x,y,/DATA 


    x0=!x.window[0]+(!x.window[1]-!x.window[0])*(x-!x.crange[0])/(!x.crange[1]-!x.crange[0])
   y0=!y.window[0]+(!y.window[1]-!y.window[0])*(y-!y.crange[0])/(!y.crange[1]-!y.crange[0])

   x1=x0+(!x.window[1]-!x.window[0])*(xlen)/(!x.crange[1]-!x.crange[0])
   y1=y0+(!y.window[1]-!y.window[0])*(ylen)/(!y.crange[1]-!y.crange[0])

   loadct,ptable_,/silent

   plot, [0,1],[0,1],/noerase,xrange=[min(x2d),max(x2d)],yrange=[min(y2d),max(y2d)],$
         position=position,/nodata,xstyle=1,ystyle=1,/norm,title=title,xtitle=xtitle,ytitle=ytitle,$
         ytickname=ytickname,xtickname=xtickname,charsize=charsize
endif

end
