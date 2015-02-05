;+
; NAME:
;       MDAP_SETUP_IO
;
; PURPOSE:
;       Produce file and directory I/O strings.
;
; CALLING SEQUENCE:
;       MDAP_SETUP_IO, $
;               root_name, output_root_dir, datacube_name, file_root, output_dir, $
;               output_file_root
;
; INPUTS:
;       root_name string
;               Root name of the path to the input fits file. Read from
;               total_filelist in the MANGA_DAP configuration file.
;               E.g: ./redux/v1_0_0/7443/stack/manga-7443-12701-LOGCUBE
;
;       output_root_dir string
;               Path of output directory for this version of the DAP. (Read from
;               MANGA_DAP configuration file; defined by environmental variable
;               $MANGA_SPECTRO_ANALYSIS.)  E.g: ./analysis/v0_9_0
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;
; OUTPUT:
;       datacube_name string
;               Name of fits file to process
;
;       file_root string
;               Root name for output files.  E.g: manga-7443-12701-LOGCUBE
;
;       output_dir string
;               Path to output directory.  This is now exactly the same
;               as output_root_dir.
;
;               OBSOLETE: output_root_dir/file_root
;               CURRENT: output_root_dir
;
;       output_file_root
;               Root name for all output files.
;
;               output_dir/file_root_
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
;       05 Sep 2014: (KBW) Original implementation
;       09 Sep 2014: (KBW) No longer outputs the output fits file or idl session
;       28 Nov 2014: (KBW) Slight adjustment to output_dir and
;                          output_file_root for new directory structure
;-
;------------------------------------------------------------------------------

PRO MDAP_SETUP_IO, $
        root_name, output_root_dir, datacube_name, file_root, output_dir, output_file_root

        datacube_name=root_name+'.fits'                 ; Set the input file
        if FILE_TEST(datacube_name) eq 0 then $
            datacube_name=root_name+'.fits.gz'          ; Try a compressed file
        if FILE_TEST(datacube_name) eq 0 then begin
            print, 'Tried: '+root_name+'.fits'
            print, 'Tried: '+root_name+'.fits.gz'
            message, 'Fits file not found!'
        endif
        print, datacube_name+' found'

        ; Root name of the output file (removes all the path
        ; information)
        file_root = STRMID(root_name, STRPOS(root_name, '/', /reverse_search)+1)
        print, 'File root: '+file_root

        ; Initialize the output directory
        ;output_dir=output_root_dir+'/'+file_root
        output_dir=output_root_dir
        ; Initialize the output file root name; this is placed in front of every output file
        output_file_root=output_dir+'/'+file_root+'_'

        print, 'Output dir: '+output_dir
        print, 'Output file root: '+output_file_root

        ; If the output directory does not exist create it
        if FILE_TEST(output_dir, /directory) eq 0 then $
            SPAWN, 'mkdir -p '+output_dir

END

