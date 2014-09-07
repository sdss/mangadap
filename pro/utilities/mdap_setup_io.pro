;+
; NAME:
;	MDAP_SETUP_IO
;
; PURPOSE:
;	Produce file and directory I/O strings.
;
; CALLING SEQUENCE:
;	MDAP_SETUP_IO, $
;		root_name, output_root_dir, datacube_name, file_root, output_dir, $
;		output_file_root, output_filefits, output_idlsession
;
; INPUTS:
;	root_name string
;		Root name of the path to the input fits file. Read from
;		total_filelist in the MANGA_DAP configuration file.
;		E.g: ./redux/v1_0_0/7443/stack/manga-7443-12701-LOGCUBE
;
;	output_root_dir string
;		Path of output directory for this version of the DAP. (Read from
;		MANGA_DAP configuration file; defined by environmental variable
;		$MANGA_SPECTRO_ANALYSIS.)  E.g: ./analysis/v0_9_0
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;
; OUTPUT:
;	datacube_name string
;		Name of fits file to process
;
;	file_root string
;		Root name for output files.  E.g: manga-7443-12701-LOGCUBE
;
;	output_dir string
;		Path to output directory.
;		E.g:./analysis/v0_9_0/manga-7443-12701-LOGCUBE
;
;	output_file_root
;		Root name for all output files.
;		E.g: ./analysis/v0_9_0/manga-7443-12701-LOGCUBE/manga-7443-12701-LOGCUBE_
;
;	output_filefits
;		Name of fits output file.  E.g: 
;		./analysis/v0_9_0/manga-7443-12701-LOGCUBE/manga-7443-12701-LOGCUBE_high_level.fits
;
;	output_idlsession
;		Name of idl session file output. E.g.:
;		./analysis/v0_9_0/manga-7443-12701-LOGCUBE/manga-7443-12701-LOGCUBE_mdap_session.idl
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
;	05 Sep 2014: (KBW) Original implementation
;-
;------------------------------------------------------------------------------

PRO MDAP_SETUP_IO, $
	root_name, output_root_dir, datacube_name, file_root, output_dir, output_file_root, $
	output_filefits, output_idlsession

	datacube_name=root_name+'.fits'  		; Set the input file
	if FILE_TEST(datacube_name) eq 0 then $
	    datacube_name=root_name+'.fits.gz'		; Try a compressed file
	if FILE_TEST(datacube_name) eq 0 then begin
	    print, 'Fits file not found!'
	    return
	endif
	print, datacube_name+' found'

	file_root = STRMID(root_name, STRPOS(root_name, '/', /reverse_search)+1)
	print, file_root

	; Initialize the output directory
	output_dir=output_root_dir+'/'+file_root
	; Initialize the output file root name; this is placed in front of every output file
	output_file_root=output_dir+'/'+file_root+'_'

	output_filefits=output_file_root+'high_level.fits'	; Output file name
	output_idlsession=output_file_root+'mdap_session.idl'	; Output idl session

	print, output_dir
	print, output_file_root
	print, output_filefits
	print, output_idlsession

;	output_dir=output_root_dir+'results_'+mode+'/'+root_name+'/'+root_name+'_'

	; If the output directory does not exist create it
	; TODO: Check for previous output
	if FILE_TEST(output_dir) eq 0 then $
	    SPAWN, 'mkdir -p '+output_dir

END

