;+
; NAME:
;	MDAP_CHECK_BINTABLE_COLUMNS
;
; PURPOSE:
;	Check that the specified binary table has the listed columns.
;
; CALLING SEQUENCE:
;	result = MDAP_CHECK_BINTABLE_COLUMNS(file, exten_no, col_list)
;
; INPUTS:
;	file string
;		Name of the fits file to check.
;
;	exten_no integer
;		Extension number.
;
;	col_list strarr[]
;		Array of column names to check for.
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;
; OUTPUT:
;	Result is 0 if one or more columns is missing, 1 if all exist.
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
;	FITS_OPEN
;	FITS_CLOSE
;	FITS_HELP
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;	13 Nov 2014: (KBW) Original implementation
;-
;------------------------------------------------------------------------------

FUNCTION MDAP_CHECK_BINTABLE_COLUMNS, $
		file, exten_no, col_list

	compile_opt idl2

	if file_test(file) ne 1 then begin
	    print, 'File does not exist!'
	    return, 0
	endif

	FITS_OPEN, file, unit
	if unit.nextend lt exten_no then begin
	    FITS_CLOSE, unit
	    print, 'File does not have exten_no='+MDAP_STC(exten_no, /integer)
	    return, 0
	endif

	if unit.xtension[exten_no] ne 'BINTABLE' then begin
	    FITS_CLOSE, unit
	    print, 'XTENSION is not a BINTABLE'
	    return, 0
	endif
	    
	FITS_READ, unit, dummy, hdr, /header_only, /no_pdu, exten_no=exten_no
	FITS_CLOSE, unit

	TBINFO, hdr, tb_str

	ncol_tab = n_elements(tb_str.ttype)
	ncol_inp = n_elements(col_list)

	if ncol_tab lt ncol_inp then begin
	    print, 'Table has insufficient number of columns.'
	    return, 0
	endif

	all_exist=1
	for i=0,ncol_inp-1 do begin
	    found=0
	    for j=0,ncol_tab-1 do begin
		if tb_str.ttype[j] eq col_list[i] then begin
		    found=1
		    break
		endif
	    endfor
	    if found eq 0 then begin
		print, 'Column '+col_list[i]+' does not exist!'
		all_exist=0
		break
	    endif
	endfor
	return, all_exist
END

