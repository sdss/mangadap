;+
; NAME:
;	MDAP_CREATE_INPUT_TABLE
;
; PURPOSE:
;	Create the input file used by MANGA_DAP.  It contains the list of galaxy
;	data to process with additional input properties used by the wrapped
;	procedures.  Columns are:
;
;	    1. [string] Name of the fits file, WITHOUT the .fits extension
;	    2. [float]  Galaxy mean velocity (km/sec)
;	    3. [float]  Galaxy mean velocity dispersion (km/sec)
;	    4. [float]  Galaxy mean ellipticity
;	    5. [float]  Galaxy mean Position Angle
;	    6. [float]  MANGA Bundle size
;	    7. [float]  Galaxy effective radius
;	    8. [string] Dataformat. Valid entries are: CUBE or RSS
;
;	The procedure uses a list of NSA IDs, matched to a list of DRP fits
;	files (see MATCH_OBS_NSA), to pull the relevant data from the NSA
;	database and create the file.
;
; CALLING SEQUENCE:
;	MDAP_CREATE_INPUT_TABLE, ifile, nsa_cat, ofile, def_vdisp=def_vdisp
;
; INPUTS:
;
;	ifile string
;		File with columns (1) fits file, (2-4) not used, (5) NSA ID.
;		This file can be generated using MATCH_OBS_NSA, or generated
;		manually.  Any NSA ID that is less than zero will be ignored.
;
;	nsa_cat string
;		NSA fits database.  E.g. nsa_v0_1_2.fits
;
;	ofile string
;		Output file name following the format given above
;
; OPTIONAL INPUTS:
;	def_vdisp double
;		Default value for the velocity dispersion. If not set, default
;		is set to -1.
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
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;	03 Sep 2014: (KBW) Original implementation
;	05 Sep 2014: (KBW) Use N_ELEMENTS instead of ARRAY_LENGTH
;	09 Sep 2014: (KBW) Change mode to either CUBE or RSS
;-
;------------------------------------------------------------------------------

PRO MDAP_CREATE_INPUT_TABLE, $
		ifile, nsa_cat, ofile, def_vdisp=def_vdisp

	; Read the input file
	READCOL, ifile, fitsfile, mangaid, objra, objdec, nsaid, comment='#', format='A,A,D,D,L'
;	n_manga = ARRAY_LENGTH(nsaid)						; Number of files
	n_manga = n_elements(nsaid)						; Number of files

	; read IFUDSGN from each fits file to get the IFU size
	ifun = MAKE_ARRAY(n_manga, /integer)				; Array to contain IFU size
	for i=0,n_manga-1 do begin
	    header=HEADFITS(fitsfile[i])
	    ifun[i]=ROUND(float(FXPAR(header,'IFUDSGN'))/100.)
	endfor

	; Read the NSA catalog to get the relevant quantities
	FXBOPEN, catdb, nsa_cat, 1, header		; Open the file

	FXBREAD, catdb, nsaid_, 'NSAID'			; NSA ID
;	n_nsa = ARRAY_LENGTH(nsaid_)			; Size of catalog
	n_nsa = n_elements(nsaid_)			; Size of catalog

	FXBREAD, catdb, redshift_, 'Z'			; Redshift

	; VDISP is not alway available
	if ~keyword_set(def_vdisp) then $
	    def_vdisp=-1.0d

	fxberr=''
 	if FXBCOLNUM(catdb, 'VDISP', errmsg=fxberr) eq 0 then begin
	    print, 'WARNING: VDISP column not available!'
	    vdisp_ = MAKE_ARRAY(n_nsa, /double, value=def_vdisp); Placeholder
	endif else begin
	    FXBREAD, catdb, vdisp_, 'VDISP'			; Velocity dispersion
	endelse

	FXBREAD, catdb, boa_, 'SERSIC_BA'		; B/A from Sersic fit
	FXBREAD, catdb, pa_, 'SERSIC_PHI'		; Position angle from Sersic fit
	FXBREAD, catdb, reff_, 'SERSIC_TH50'		; Half-light radius fro Sersic fit
	FXBCLOSE, catdb					; Close up

	ii_ = SORT(nsaid_)				; Sorted indices for full catalog
	ii = SORT(nsaid)				; Sorted indices for fits files
	jj=ii						; Indices in DRP output in NSA catalog

	; Find the NSA ID of the target galaxies and save the relevant data
	j = 0L
	for i=0,n_manga-1 do begin
	    if nsaid[ii[i]] lt 0 then begin		; No entry in NSA
		jj[ii[i]] = -1
		continue
	    endif

;	    print, i, ii[i], nsaid[ii[i]]
	    while (j lt n_nsa) and (nsaid_[ii_[j]] lt nsaid[ii[i]]) do ++j
	    
	    if (nsaid[ii[i]] ne nsaid_[ii_[j]]) then begin
		PRINT, 'Could not find NSAID: '+nsaid[ii[i]]
		RETURN
	    endif

	    jj[ii[i]] = ii_[j]
;	    print, format='(%"Found %d at index %d")', nsaid[ii[i]], jj[i]
	endfor

;	print, jj

	c=299792.458d			; TODO: there isn't some idlutils file that defines this?

	; Write the file
	OPENW, 1, ofile
	PRINTF, 1, '# '+SYSTIME()
	PRINTF, 1, '#', 'FITSFILE', 'V', 'VDISP', 'ELL', 'PA', 'IFU', 'Reff', 'FORMAT', $
		   format='( A1, A59, A14, A10, A10, A10, A4, A10, A10 )'
	for i=0,n_manga-1 do begin

	    ; Use the file name to get the dataformat	
	    if strpos(fitsfile[i], 'RSS') ne -1 then begin
		frmt='RSS'
	    endif else begin
		frmt='CUBE'
	    endelse

	    ; Output is different if the NSAID is or is not present
	    ;	redshift is converted to km/s on output
	    ;   b/a is converted to ellipticity on output
	    root = STRMID(fitsfile[i], 0, STRPOS(fitsfile[i], '.fits', /reverse_search))
	    if jj[i] ge 0 then begin 
		PRINTF, 1, root, c*redshift_[jj[i]], vdisp_[jj[i]], 1.0-boa_[jj[i]], $
			   pa_[jj[i]], ifun[i], reff_[jj[i]], frmt, $
			   format='( A60, E14.6, E10.2, E10.2, E10.2, I4, E10.2, A10 )'
	    endif else begin
		PRINTF, 1, root, -1, -1, -1, -1, ifun[i], -1, frmt, $
			   format='( A60, I14, I10, I10, I10, I4, I10, A10 )'
	    endelse
	endfor
	CLOSE, 1

END		    
	    
