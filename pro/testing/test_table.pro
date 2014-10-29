
pro test_table

;file='test.fits'
;mkhdr,hdu,4,[0],/extend
;sxdelpar, hdu, 'NAXIS1'
;sxaddpar, hdu, 'NAXIS', 0
;writefits, file, 0, hdu
;
;fxbhmake, hdu, 10, 'TEST',  'test table'
;junk_double = 1.0d
;junk_int = fix(0.1)
;fxbaddcol, 1, hdu, junk_double, 'XPOS', 'test 1' 
;fxbaddcol, 2, hdu, junk_double, 'YPOS', 'test 2' 
;fxbaddcol, 3, hdu, junk_double, 'ZPOS', 'test 3' 
;
;fxbcreate, tbl, file, hdu
;fxbfinish, tbl
;
;fxbopen, tbl, file, 'TEST', access='RW'
;
;fxbgrow, tbl, hdu, 20
;
;xpos=dindgen(20)
;ypos=dindgen(20)+5
;fxbwritm, tbl, ['XPOS', 'YPOS'], xpos, ypos
;
;fxbfinish, tbl


;      ; (1) First, create sample values
;      x = findgen(43) & y = findgen(43)+1 & psf = randomn(seed,23,21,43)
;
;      ; (2) Create primary header, write it to disk, and make extension header
;      fxhmake,header,/initialize,/extend,/date
;      fxwrite,'sample.fits',header
;      fxbhmake,header,43,'TESTEXT','Test binary table extension'
;
;      ; (3) Fill extension header with desired column names
;      fxbaddcol,1,header,x[0],'X'             ;Use first element in each array
;      fxbaddcol,2,header,y[0],'Y'             ;to determine column properties
;      fxbaddcol,3,header,psf[*,*,0],'PSF'
;
;      ; (4) Write extension header to FITS file
;      fxbcreate,unit,'sample.fits',header
;
;      for i=0,42 do begin
;	fxbwrite,unit, x[i], 1, i+1
;	fxbwrite,unit, y[i], 2, i+1
;	fxbwrite,unit, psf[*,*,i], 3, i+1
;      endfor
;
;      ; (5) Use FXBWRITM to write all data to the extension in a single call
;;      fxbwritm,unit,['X','Y','PSF'], x, y, psf
;;      fxbfinish,unit                 ;Close the file
;
;	file='sample.fits'
;	fits_info, file
;

; INITIALIZE THE HEADER
      ; Create primary header, write it to disk, and make extension header
      file = 'sample.fits'
      fxhmake,header,/initialize,/extend,/date
      fxwrite,file,header
      fxbhmake,header,1,'TESTEXT','Test binary table extension'

      ; (3) Fill extension header with desired column names
      fxbaddcol,1,header,1.0d,'X'             ;Use first element in each array
      fxbaddcol,2,header,1.0d,'Y'             ;to determine column properties
      fxbaddcol,3,header,1.0d,'PSF'

      ; (4) Write extension header to FITS file
      fxbcreate,unit,'sample.fits',header
      fxbfinish,unit

; APPEND AN IMAGE EXTENSION
    dummy=dblarr(1)
    MKHDR, header, dummy, /image
    SXADDPAR, header, 'EXTNAME', 'TESTIMG', 'Test image extention'
    WRITEFITS, file, dummy, header, /append


fits_info, file
stop

	FXBOPEN, unit, file, 'TESTEXT', header, access='RW'
	print, header
	print, fxbdimen(unit, 'X')
	print, fxbdimen(unit, 'Y')
	print, fxbdimen(unit, 'PSF')
	fxbfinish, unit

; CREATE THE DATA ARRAYS
      x = dindgen(43)
      y = dindgen(43)+1
      psf = randomn(seed,23,21,43)

; ADJUST THE TABLE HEADER
      fxbhmake,header,1,'TESTEXT','Test binary table extension'
      fxbaddcol,1,header,x[0],'X'             ;Use first element in each array
      fxbaddcol,2,header,y[0],'Y'             ;to determine column properties
      fxbaddcol,3,header,psf[*,*,0],'PSF'

      MODFITS, file, 0, header, extname='TESTEXT'

	FXBOPEN, unit, file, 'TESTEXT', header, access='RW'
	print, header
	print, fxbdimen(unit, 'X')
	print, fxbdimen(unit, 'Y')
	print, fxbdimen(unit, 'PSF')
	fxbfinish, unit

fits_info, file
stop

; ADJUST THE TABLE SIZE
	FXBOPEN, unit, file, 'TESTEXT', header, access='RW'
	FXBGROW, unit, header, 43
	fxbfinish, unit

fits_info, file
stop

; WRITE THE DATA
	FXBOPEN, unit, file, 'TESTEXT', header, access='RW'
      fxbwritm,unit,['X','Y'], x, y
      FXBFINISH, unit

stop

	FXBOPEN, unit, file, 'TESTEXT', header, access='RW'
      fxbwritm,unit,['PSF'], psf
      FXBFINISH, unit

      x = dindgen(43)+3
      y = dindgen(43)-4

	FXBOPEN, unit, file, 'TESTEXT', header, access='RW'
      fxbwritm,unit,['X','Y'], x, y
      FXBFINISH, unit

fits_info, file
stop

; ADD DATA TO THE IMAGE EXTENSION
	img=randomn(see,100,100)
	mkhdr, header, img, /image
	sxaddpar, header, 'EXTNAME', 'TESTIMG', 'Test image'
	modfits, file, img, header, extname='TESTIMG'

;      for i=0,42 do begin
;	fxbwrite,unit, x[i], 1, i+1
;	fxbwrite,unit, y[i], 2, i+1
;	fxbwrite,unit, psf[*,*,i], 3, i+1
;      endfor

fits_info, file

;
;
;
;errmsg = ''
;unit=fxposit(file, 'TEST', errmsg=errmsg)
;print, unit
;
;fxbopen, unit, file, 'TEST' 
;fxbread, unit, xpos, 'XPOS'

END



