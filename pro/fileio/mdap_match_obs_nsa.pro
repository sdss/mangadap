;+
; NAME:
;	MDAP_MATCH_OBS_NSA
;
; PURPOSE:
;	Provided a list of DRP fits files, crossmatch the OBJRA and OBJDEC with
;	the NSA catalog to grab the NSAID numbers.
;
; CALLING SEQUENCE:
;	MDAP_MATCH_OBS_NSA, fitslst, nsa_cat, sep, ofile
;
; INPUTS:
;	fitslst string
;		File where the first column is a list of DRP produced fits
;		files.
;
;	nsa_cat string
;		File with the NSA catalog to search.  CANNOT handle gzipped
;		files.
;
;	sep double
;		Maximum separation to use in matching the coordinates in arcsec.
;		Value is typecast as a double, even if it is not entered as
;		such.  Code DOES NOT handle multple matches.  The first target
;		within this separation is considered to be the matching target.
;
;	ofile string
;		File for procedure output, listing the fits file name, MaNGA ID,
;		OBJRA, OBJDEC, NSA ID, NSA RA, and NSA DEC.
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;
; OUTPUT:
;
; OPTIONAL OUTPUT:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; TODO:
;	- Figure out how to handle gzipped binary tables
;	- CONSTANTS isn't in my path; other default definitions of pi?
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;	03 Sep 2014: (KBW) Original impelementation
;	05 Sep 2014: (KBW) Use N_ELEMENTS instead of ARRAY_LENGTH
;-
;------------------------------------------------------------------------------

PRO MDAP_MATCH_OBS_NSA, $
		fitslst, nsa_cat, sep, ofile

	READCOL, fitslst, fitsfile, comment='#', format='A', /silent	; Read the fits files
	nfiles=n_elements(fitsfile)					; Number of files

	; Prepare to read the data
	mangaid = MAKE_ARRAY(nfiles, /string)
	objra = MAKE_ARRAY(nfiles, /double)
	objdec = MAKE_ARRAY(nfiles, /double)
	for i=0,nfiles-1 do begin
	    header=HEADFITS(fitsfile[i])				; Read the header
	    mangaid[i]=STRCOMPRESS(FXPAR(header,'MANGAID'), /remove_all); ... and the MaNGAID
	    objra[i]=FXPAR(header,'OBJRA')				; ... and the RA
	    objdec[i]=FXPAR(header,'OBJDEC')				; ... and the DEC
;	    print, fitsfile[i], mangaid[i], objra[i], objdec[i]
	endfor

	; Read the NSA catalog data
	get_lun, catdb
	FXBOPEN, catdb, nsa_cat, 1, header
	FXBREAD, catdb, nsaid, 'NSAID'
	FXBREAD, catdb, nsara, 'RA'
	FXBREAD, catdb, nsadec, 'DEC'
	FXBCLOSE, catdb
	free_lun, catdb

	n_nsa=n_elements(nsaid)						; Catalog length

	; Match the objra and objdec to nsara and nsadec
	nsa_index = MAKE_ARRAY(nfiles, /long)				; Index in nsa for each obs

	;CONSTANTS, pi=pi			; Get the value of pi
	pi=3.14159d
	d2r=pi/180				; Conversion from deg to radians
	d=(double(sep)/3600.)^2			; sep is in arcsec, convert to deg^2
	
	; Brute force it
	for i=0,nfiles-1 do begin
	    cosd=cos(objdec[i]*d2r)
	    for j=0L,n_nsa-1 do begin
		dd=((nsara[j]-objra[i])*cosd)^2+(nsadec[j]-objdec[i])^2
;		print, dd, d
		if dd lt d then begin
		    nsa_index[i] = j
		    break
		endif
	    endfor
	    if j eq n_nsa then begin
		nsa_index[i] = -1
		print, 'Could not find entry in NSA for file '+fitsfile[i]
	    endif
	endfor

	; Write the file
	OPENW, unit, ofile, /get_lun
	PRINTF, unit, '# '+SYSTIME()
	PRINTF, unit, '#', 'FITSFILE', 'MANGAID', 'OBJRA', 'OBJDEC', 'NSAID', 'NSARA', 'NSADEC', $
		   format='( A1, A59, A15, A11, A11, A15, A11, A11 )'
	for i=0,nfiles-1 do $
	    if (nsa_index[i] ge 0) then begin
		PRINTF, unit, fitsfile[i], mangaid[i], objra[i], objdec[i], nsaid[nsa_index[i]], $
			   nsara[nsa_index[i]], nsadec[nsa_index[i]], $
			   format='( A60, A15, F11.5, F11.6, I15, F11.5, F11.6 )'
	    endif else begin
		PRINTF, unit, fitsfile[i], mangaid[i], objra[i], objdec[i], -1, -1.0, -1.0, $
			   format='( A60, A15, F11.5, F11.6, I15, F11.5, F11.6 )'
	    endelse
	CLOSE, unit
	free_lun, unit

END

