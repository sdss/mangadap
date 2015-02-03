;+
; NAME:
;       MDAP_LOG_INSTANTIATE
;
; PURPOSE:
;       Instantiate the main log file for the DAP.
;
; CALLING SEQUENCE:
;       MDAP_LOG_INSTANTIATE, output_dir, manga_dap_version, signifier, datacube_name, mode, $
;                             file_root, guess_vel, guess_sig, ell, pa, Reff, file_unit, $
;                             inptbl=inptbl, index=index, par=par
;
; INPUTS:
;       output_dir string
;               Output directory for current run of the DAP.
;
;       manga_dap_version string
;               DAP version number.
;
;       signifier string
;               Unique signifier of the mdap_execution_setup procedure
;               used.
;
;       datacube_name string
;               Name of the datacube to analyze.
;
;       mode string
;               DRP file type: CUBE or RSS.
;
;       file_root string
;               Root name for the DAP output files.
;
;       guess_vel double
;               Guess velocity to use.
;
;       guess_sig double
;               Guess velocity dispersion to use.
;
;       ell double
;               Ellipticity of the galaxy isophotes; used to get the
;               in-plane radius.
;
;       pa double
;               Position angle of the galaxy isophotes; used to get the
;               in-plane radius.
;
;       Reff double
;               Effective radius of the galaxy; used to renormalize the
;               galaxy radii.
;
; OPTIONAL INPUTS:
;       inptbl string
;           File with the input table providing the necessary input
;           parameters.  See MDAP_CREATE_INPUT_TABLE for the file
;           format.  Should not be provided along with par.  Must also
;           provide index for use.
;
;       index integer
;           Row index in inptbl with the parameters to use.
;
;       par string
;           SDSS par (yanny) file with the parameters to use.  Read
;           using YANNY_READONE(), meaning that only the first entry
;           will be read, but others can be present.
;
; OPTIONAL KEYWORDS:
;
; OUTPUT:
;       file_unit unit
;               File unit to write to.  Unit is returned open.
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
;           30 Jan 2015: (KBW) Pulled from manga_dap.pro
;-
;------------------------------------------------------------------------------

;-------------------------------------------------------------------------------
; Get a date string to use as the extension of the log file name
FUNCTION MDAP_FILE_EXTENSION_DATE

        vec = fix(date_conv( systime(/julian,/utc), 'V'))
        ext = MDAP_STC(vec[0], /integer)+'_'+$ 
              MDAP_STC(vec[1], /integer)+'_'+$
              MDAP_STC(vec[2], /integer)+'_'+$
              MDAP_STC(vec[3], /integer)+'_'+$
              MDAP_STC(vec[4], /integer)

        return, ext
END

PRO MDAP_LOG_INSTANTIATE, $
                output_dir, manga_dap_version, signifier, datacube_name, mode, file_root, $
                guess_vel, guess_sig, ell, pa, Reff, file_unit, inptbl=inptbl, index=index, par=par

        ofile = output_dir+'/mdap_'+MDAP_FILE_EXTENSION_DATE()+'.log'
        openw, file_unit, ofile, /get_lun
        
        cd, current=dir
        printf, file_unit, 'DAP called from: ', dir
        printf, file_unit, 'Start systime(): ', systime()
        printf, file_unit, ''
        printf, file_unit, 'DAP Version: ', manga_dap_version
        printf, file_unit, ''
        printf, file_unit, 'MDAP_EXECUTION_SETUP signifier: ', signifier
        if n_elements(inptbl) ne 0 and n_elements(index) ne 0 then begin
            printf, file_unit, 'Input table: ', inptbl
            printf, file_unit, 'Line to process: ', (index+1)
        endif else if n_elements(par) ne 0 then $
            printf, file_unit, 'Parameter file: ', par
        printf, file_unit, ''
        printf, file_unit, 'Fits file to process: ', datacube_name
        printf, file_unit, 'File type: ', mode
        printf, file_unit, 'Output directory: ', output_dir
        printf, file_unit, 'Root name for output files: ', file_root
        printf, file_unit, ''
        printf, file_unit, 'Input parameters (or guesses):'
        printf, file_unit, '            Velocity: ', guess_vel
        printf, file_unit, '     Vel. dispersion: ', guess_sig
        printf, file_unit, '         Ellipticity: ', ell
        printf, file_unit, '      Position Angle: ', pa
        printf, file_unit, '    Effective radius: ', Reff
        printf, file_unit, ''

END

;output_filefits, 
;       printf, file_unit, 'Data products:'
;       printf, file_unit, '    Primary fits file: ', output_filefits
;       printf, file_unit, ''

