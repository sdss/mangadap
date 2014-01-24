pro mdap_readfits_lambda,filename,output_fits,lambda_vector,header,lambda_direction=lambda_direction,OCRPIX1=OCRPIX1,OCRPIX2=OCRPIX2,OCDELT1=OCDELT1,OCDELT2=OCDELT2,OCRVAL1=OCRVAL1,OCRVAL2=OCRVAL2,exten_no=exten_no,help=help



; NEEDS TO BE RE-ADAPTED: read_header must be replaced. useless manga
; keywords must be removed.

; PURPOSE:

; It reads a fits file using readfits function, but it
; returns also the header, and the vector with lambda values (wavelenghts).
; Useful for long slit spectra or 2D spectra arranged as long slit.


; INPUTS:
;
; filename      [string]     Name of the fits file to read 
;
;
; OUTPUTS:
;
; output_fits   [2D array]   The image, read by readfits function (/noscale and /silent
;                            keywords will be used). Its dimensions are
;                            sz1 and sz2
;
; lambda_vector [1D array]   The vector of lambda values. Its dimension
;                            is sz1 id wavelenght direction is along
;                            x(default) or sz2 if wavelenght direction
;                            is along y axis.
; header    [string array]   A string array containing image
;                            header. Read by readfits.
;


; OPTIONAL KEYWORDS:
;
; lambda_direction [string]  If set to 'x' (default), then wavelength
;                            direction is assumed to be x axis. If set
;                            to 'y', then wavelength direction is
;                            assumed to be y axis. All other values
;                            are converted to 'x'

; OCRPIX1          [string]  alternative value for standard header 'CRPIX1': reference pixel lambda coordinate               
; OCRPIX2          [string]  alternative value for standard header 'CRPIX2': reference pixel slit coordinate             
; OCDELT1          [string]  alternative value for standard header 'CDELT1': pixel step along lambda                         
; OCDELT2          [string]  alternative value for standard header 'CDELT2': pixel step along slit                       
; OCRVAL1          [string]  alternative value for standard header 'CRVAL1': pixel value along lambda at reference point   
; OCRVAL2          [string]  alternative value for standard header 'CRVAL2': pixel value along slit at reference point     
; exten_no         int       extension number of the fits file to read (default = 0)
;
;
; DEPENDENCIES:
;
; It requires read_header.pro and readfits.pro package.

; HISTORY:

; version 1.0   24 october 2007 by Lodovico Coccato
; version 1.1   13 January 2008 by Lodovico Coccato. Some errors in the comments fixed
; version 1.2   02 May 2008 by Lodovico Coccato. Help keyword added
; version 1.3   22 Oct 2008 by L. Coccato. Added /double and removed
;               /noscale 

; version 1.4   23 Jul 2008 by L. Coccato. Added exten_no optional input




; START CODE
;if keyword_set(help) then goto, print_help
if n_elements(lambda_direction) eq 0 then lambda_direction='x'
if lambda_direction ne 'y' then lambda_direction='x'

output_fits=readfits(filename,header,/silent,exten_no=exten_no)


sz=size(output_fits)



;if w direction is along x axis

if lambda_direction eq 'x' then begin
   CRPIX1='CRPIX1'  ;reference pixel lambda coordinate               
   CRPIX2='CRPIX2'  ;reference pixel slit coordinate             
   CDELT1='CDELT1'  ;pixel step along lambda                       
   CDELT2='CDELT2'  ;pixel step along slit                     
   CRVAL1='CRVAL1'  ;pixel value along lambda at reference point   
   CRVAL2='CRVAL2'  ;pixel value along slit at reference point     
   sz1=sz(1)
endif

;if w direction is along y axis
if lambda_direction eq 'y' then begin
   CRPIX2='CRPIX1'  ;reference pixel slit coordinate
   CRPIX1='CRPIX2'  ;reference pixel lambda coordinate
   CDELT2='CDELT1'  ;pixel step along slit
   CDELT1='CDELT2'  ;pixel step along lambda
   CRVAL2='CRVAL1'  ;pixel value along slit at reference point
   CRVAL1='CRVAL2'  ;pixel value along lambda at reference point
   sz1=sz(2)
endif

;If I provide different strings for the (default) header entries
if n_elements(OCRPIX1) ne 0 then CRPIX1 = OCRPIX1
if n_elements(OCRPIX2) ne 0 then CRPIX2 = OCRPIX2
if n_elements(OCDELT1) ne 0 then CDELT1 = OCDELT1
if n_elements(OCDELT2) ne 0 then CDELT2 = OCDELT2
if n_elements(OCRVAL1) ne 0 then CRVAL1 = OCRVAL1
if n_elements(OCRVAL2) ne 0 then CRVAL2 = OCRVAL2

read_header,filename,CRPIX1,a
d=execute('ref_pix1='+strcompress(a,/remove_all))
;read_header,filename,CRPIX2,a
;d=execute('ref_pix2='+strcompress(a,/remove_all))
read_header,filename,CDELT1,a
d=execute('step_pix1='+strcompress(a,/remove_all))
;read_header,filename,CDELT2,a
;d=execute('step_pix2='+strcompress(a,/remove_all))
read_header,filename,CRVAL1,a
d=execute('ref_val1='+strcompress(a,/remove_all))
;read_header,filename,CRVAL2,a
;d=execute('ref_val2='+strcompress(a,/remove_all))
sz=size(galaxy)
ref_pix1 = double(ref_pix1)
step_pix1=double(step_pix1)
ref_val1=double(ref_val1)

starty=-(ref_pix1-1.d)*step_pix1+ref_val1
lambda_vector=double(findgen(sz1)*step_pix1 + starty)

goto, fine

;print_help:
;
;print, '; PURPOSE:'
;print, ''
;print, '; It reads a fits file using readfits function, but it'
;print, '; returns also the header, and the vector with lambda values (wavelenghts).'
;print, '; Useful for long slit spectra or 2D spectra arranged as long slit.'
;print, ';'
;print, ';'
;print, '; USAGE:'
;print, ';   readfits_lambda,filename,output_fits,lambda_vector,header,[lambda_direction=lambda_direction],[OCRPIX1=OCRPIX1],'
;print, ';              [OCRPIX2=OCRPIX2],[OCDELT1=OCDELT1],[OCDELT2=OCDELT2],[OCRVAL1=OCRVAL1],[OCRVAL2=OCRVAL2],[help=help]'
;print, ';'
;print, ';'
;print, '; INPUTS:'
;print, ';'
;print, '; filename      [string]     Name of the fits file to read '
;print, ';'
;print, ';'
;print, '; OUTPUTS:'
;print, ';'
;print, '; output_fits   [2D array]   The image, read by readfits function (/noscale and /silent'
;print, ';                            keywords will be used). Its dimensions are'
;print, ';                            sz1 and sz2'
;print, ';'
;print, '; lambda_vector [1D array]   The vector of lambda values. Its dimension'
;print, ';                            is sz1 id wavelenght direction is along'
;print, ';                            x(default) or sz2 if wavelenght direction'
;print, ';                            is along y axis.'
;print, '; header    [string array]   A string array containing image'
;print, ';                            header. Read by readfits.'
;print, ';'
;print, ';'
;print, ';'
;print, '; OPTIONAL KEYWORDS:'
;print, ';'
;print, '; lambda_direction [string]  If set to "x" (default), then wavelength'
;print, ';                            direction is assumed to be x axis. If set'
;print, ';                            to "y", then wavelength direction is'
;print, ';                            assumed to be y axis. All other values'
;print, ';                            are converted to "x"'
;print, ';'
;print, '; OCRPIX1          [string]  alternative value for standard header "CRPIX1": reference pixel lambda coordinate '              
;print, '; OCRPIX2          [string]  alternative value for standard header "CRPIX2": reference pixel slit coordinate  '           
;print, '; OCDELT1          [string]  alternative value for standard header "CDELT1": pixel step along lambda         '                
;print, '; OCDELT2          [string]  alternative value for standard header "CDELT2": pixel step along slit          '             
;print, '; OCRVAL1          [string]  alternative value for standard header "CRVAL1": pixel value along lambda at reference point '  
;print, '; OCRVAL2          [string]  alternative value for standard header "CRVAL2": pixel value along slit at reference point '    

fine:
end
