;+
; NAME:
;       MDAP_DRPALL_INFO
;
; PURPOSE:
;       Pull information from the DRPall file required as input for the DAP.
;       FOR NOW, this actually pulls information from the provided NSA file!
;
; CALLING SEQUENCE:
;       MDAP_DRA_ALL_INFO, list, nsa_cat, match_file, input_table
;
; INPUTS:
;       list string
;               File with the list of DRP files to process.
;       
;       match_file string
;               File for the match of the DRP files to the NSA ID.
;
;       input_table string
;               File for the DAP input file data.
;
; OPTIONAL INPUTS:
;       nsa_cat string
;               NSA catalog file from which to draw information.  Default is
;               nsa_cat=getenv('MANGAWORK_DIR')+'/manga/target/temp/12-nsa_v1b_0_0_v2.fits.gz'
;
; OUTPUT:
;       Output is in the form of the data written to the match_file and the
;       input_table.
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
;       13 Nov 2014: (KBW) Original implementation (changed to force hard-wired file names)
;-
;------------------------------------------------------------------------------

PRO MDAP_DRPALL_INFO
;               list, match_file, input_table, nsa_cat=nsa_cat

        ; TODO: These are hardwired files !!
        output_root=getenv('MANGA_SPECTRO_ANALYSIS')+'/'+getenv('MANGADAP_VER')

        drp_list=output_root+'/drp_completed.lst'
        nsa_match=output_root+'/manga_drp_nsa_match.inp'
        input_table=output_root+'/manga_dap_table.inp'

        nsa_cat=getenv('MANGAWORK_DIR')+'/manga/target/temp/12-nsa_v1b_0_0_v2.fits.gz'

        ; Test the required files exist
        if file_test(drp_list) eq 0 then $
            message, 'Completed DRP files does not exist: '+drp_list
        if file_test(nsa_cat) eq 0 then $
            message, 'NSA catalog does not exist: '+nsa_cat

        ; Match the observed galaxies to their NSA counterparts
        ; TODO: What happens when a galaxy is not found?
        MDAP_MATCH_OBS_NSA, drp_list, nsa_cat, 10, nsa_match
        ; Grab information from the NSA catalog to create the DAP input table
        MDAP_CREATE_INPUT_TABLE, nsa_match, nsa_cat, input_table, def_vdisp=100.0

END

