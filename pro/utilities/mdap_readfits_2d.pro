pro extr_str,file,head,record_value,head_string=head_string;,silent=silent

; By Lodovico Coccato
; version: 29 Mar 2005

; Description:
; Extract the 'head' value from the header of the image 'file' and
; stored it in 'record_value'. The heaedr can be directly supplied by
; the string input head_string


; Example:
; extr_str,'galassia.fits','CCDNAME',results
;
; The procedure will look at the header of the file galassia.fits and 
; search for the CCDNAME.
; If it find that record, the results will be the value of that 
; field (ex:  'EEV13             ' ). If the record is not found the 
; results will be the string 'No record found in header'.


; --- Start program ---



; -- Inputs parameters --
; file    :name of the input file (string)
; head    : keywords in the header to look for (string)

record_name_length=8   ; length of the record name in the image header (spaces included)
record_value_length=22 ; lenght of the record value name in the image header (spaces included)
if n_elements(head_string) eq 0 then begin
   junk=READFITS(file,header,/noscale,/silent)
endif else begin
   header = head_string
endelse
head = STRCOMPRESS(head,/remove_all)

record_value='No record found in header'

FOR i = 0, n_elements(header)-1 DO BEGIN
   key = STRMID(header(i),0,record_name_length)
   key = STRCOMPRESS(key,/remove_all)
   IF (key EQ head) THEN BEGIN
      record_value = STRMID(header(i),record_name_length+1,record_value_length)
      goto, jump 
   ENDIF
ENDFOR
jump:

end
; NAME:
; mdap_readfits_2d
;
; PURPOSE: 
; Read a fits image using readfits, and it returns also the matrices
; of X and Y coordinates (world coordinate system) by reading, pixel
; scale, reference pixel, and reference pixel value from the image
; hader.

; USAGE:
; readfits_2d,filename,image,[header=header],[x2d=x2d],[y2d=y2d],[/noscale],
;              [x_pixel2d=x_pixel2d],[y_pixel2d=y_pixel2d]   
; INPUTS:
; 
; filename     [string]    name of the fits file to read
;
; OUTPUTS:
; image        [2D array]  Two dimensional array read from the fits
;                          file with readfits 
;
; OPTIONAL OUTPUTS:
;
; header   [string array]
;
; x2d          [2D array]  Matrix of x world coordinate values,
;                          computed from image header
;
; y2d          [2D array]  Matrix of y world coordinate values,
;                          computed from image header
;
; x_pixel2d    [2D array]  Matrix of x in pixel coordinate
;
; y_pixel2d    [2D array]  Matrix of y in pixel coordinate
;

; HISTORY:
; Version 1.0  29 Oct 2008 by L. Coccato
; Version 1.1  30 Oct 2008 by L. Coccato. x_pixel2d and y_pixel2d added



; OPTIONAL KEYWORDS:
; /noscale    If set, the keyword /noscale will be applied to readfits.



; START MAIN CODE

pro mdap_readfits_2d,filename,image,header=header,x2d=x2d,y2d=y2d,noscale=noscale,$
x_pixel2d=x_pixel2d,y_pixel2d=y_pixel2d,silent=silent

if keyword_set(noscale) then begin
   image=readfits(filename,header,/silent,/noscale)
endif else begin
   image=readfits(filename,header,/silent)
endelse

; CALCULATION OF X and Y MATRICES FROM HEADER

; Extract information on dimensions of lenslet and reference
; coordinates.
ref_pix1=0./0.
ref_pix2=0./0.
step_pix1=0./0.
step_pix2=0./0.
ref_val1=0./0.
ref_val2=0./0.
if keyword_set(silent) then begin 
   extr_str,filename,'CRPIX1',a
   if a ne 'No record found in header' then d=execute('ref_pix1='+strcompress(a,/remove_all)) 
   extr_str,filename,'CRPIX2',a
   if a ne 'No record found in header' then d=execute('ref_pix2='+strcompress(a,/remove_all)) 
   extr_str,filename,'CDELT1',a
   if a ne 'No record found in header' then d=execute('step_pix1='+strcompress(a,/remove_all))
   extr_str,filename,'CDELT2',a
   if a ne 'No record found in header' then d=execute('step_pix2='+strcompress(a,/remove_all))
   extr_str,filename,'CRVAL1',a
   if a ne 'No record found in header' then d=execute('ref_val1='+strcompress(a,/remove_all))
   extr_str,filename,'CRVAL2',a
   if a ne 'No record found in header' then d=execute('ref_val2='+strcompress(a,/remove_all))
   
endif else begin
   extr_str,filename,'CRPIX1',a
   d=execute('ref_pix1='+strcompress(a,/remove_all)) 
   extr_str,filename,'CRPIX2',a
   d=execute('ref_pix2='+strcompress(a,/remove_all)) 
   extr_str,filename,'CDELT1',a
   d=execute('step_pix1='+strcompress(a,/remove_all))
   extr_str,filename,'CDELT2',a
   d=execute('step_pix2='+strcompress(a,/remove_all))
   extr_str,filename,'CRVAL1',a
   d=execute('ref_val1='+strcompress(a,/remove_all))
   extr_str,filename,'CRVAL2',a
   d=execute('ref_val2='+strcompress(a,/remove_all))
endelse


sz=size(image)
if n_elements(ref_pix1) eq 0 or finite(ref_pix1) eq 0 then goto, mmm

startx=-(ref_pix1-1.)*step_pix1+ref_val1 
x=indgen(sz[1])*step_pix1 +startx

starty=-(ref_pix2-1.)*step_pix2+ref_val2 
y=indgen(sz[2])*step_pix2 +starty


x2d=x#(y*0.0+1) 
y2d=(x*0.0+1)#y 

mmm:
x=indgen(sz[1]) +1
y=indgen(sz[2]) +1

x_pixel2d=x#(y*0.0+1) 
y_pixel2d=(x*0.0+1)#y 


end
