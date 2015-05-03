;+
; NAME:
;       MDAP_READ_INPUT_SETUP
;
; PURPOSE:
;       Setup the information used by the DAP to define the input DRP
;       file, its path, the output path, and some analysis parameters.
;
; CALLING SEQUENCE:
;       MDAP_READ_INPUT_SETUP, inptbl=inptbl, index=index, par=par, drppath=drppath, $
;                              dappath=dappath, plate, ifudesign, mode, velocity_initial_guess, $
;                              velocity_dispersion_initial_guess, ell, pa, Reff, root_name, $
;                              output_root_dir
;
; INPUTS:
;
; OPTIONAL INPUTS:
;       inptbl string
;               Input ASCII file with the necessary information
;               tabulated in a file.  Returned data is set by the index
;               parameter.  Both must be provided to use this output
;               capability.  See MDAP_CREATE_INPUT_TABLE.
;
;       index integer
;               Index of the data row (starting from 0) from inptbl to
;               output.
;
;       par string
;               SDSS par (Yanny) file with the required information in
;               the DAPPAR struct.
;
;       drppath string
;               Path to the DRP file to use instead of the default path.
;               The default path is
;
;               def_drppath = getenv('MANGA_SPECTRO_REDUX') + '/' + $
;                             getenv('MANGADRP_VER') + '/' + $
;                             MDAP_STC(plate, /integer) + '/stack'
;
;       dappath string
;               Path for the DAP output to use instead of the default
;               path.  The default path is
;
;               def_dappath = getenv('MANGA_SPECTRO_ANALYSIS') + '/' $
;                             + getenv('MANGADRP_VER') + '/' $
;                             + getenv('MANGADAP_VER')
;
; OPTIONAL KEYWORDS:
;
; OUTPUT:
;       plate integer
;               Plate designation of the observed data.
;
;       ifudesign integer
;               IFU design of the observation, which is the size of the
;               fiber bundle appended to the ID of the bundle.  E.g.,
;               12703.
;
;       mode string
;               Mode of the DRP file.  Valid values are either 'CUBE' or
;               'RSS'.
;
;       velocity_initial_guess double
;               Initial guess for the velocity (cz) of the system in
;               km/s.
;
;       velocity_dispersion_initial_guess double
;               Initial guess for the velocity dispersion of the system
;               in km/s.
;
;       ell double
;               Ellipticity (1-b/a) of the object used to define
;               in-plane radius (semi-major radius).
;
;       pa double
;               On-sky position angle (angle from N through E) of the
;               elliptical isophotal contours used to define the
;               in-plane (semi-major) radius.
;
;       Reff double
;               Effective radius of the galaxy in arcsec.
;
;       root_name string
;               Root name of the DRP file.  This is
;
;               drppath + '/manga-' + MDAP_STC(plate, /integer) + $
;               '-' + MDAP_STC(ifudesign, /integer) + '-LOG' + mode
;
;       output_root_dir string
;               Root name of the output directory.  This is
;
;               dappath + '/' + MDAP_STC(plate, /integer) + '/' + $
;               MDAP_STC(ifudesign, /integer) + '/'
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
;       28 Nov 2014: Original implementation by K. Westfall (KBW)
;       29 Apr 2015: (KBW) Change output DAP path to include DRP version
;-
;-----------------------------------------------------------------------

PRO MDAP_READ_INPUT_SETUP, $
        inptbl=inptbl, index=index, par=par, drppath=drppath, dappath=dappath, plate, ifudesign, $
        mode, velocity_initial_guess, velocity_dispersion_initial_guess, ell, pa, Reff, $
        root_name, output_root_dir

        ; Check the argument input
;       if n_elements(inptbl) eq 0 and n_elements(par) eq 0 then $
        if n_elements(inptbl) eq 0 && n_elements(par) eq 0 then $
            message, 'Must provide either inptbl or par!'

;       if n_elements(inptbl) ne 0 and n_elements(index) eq 0 then $
        if n_elements(inptbl) ne 0 && n_elements(index) eq 0 then $
            message, 'Must provide index of line number to use in '+inptbl+'!'

;       if n_elements(inptbl) ne 0 and n_elements(par) ne 0 then $
        if n_elements(inptbl) ne 0 && n_elements(par) ne 0 then $
            message, 'Cannot provide both an input table and a parameter file.  Pick one!'

        ; Get the parameter data
        if n_elements(inptbl) ne 0 then begin
            if file_test(inptbl) eq 0 then $
                message, 'Cannot open '+inptbl+'!'
            READCOL, inptbl, plate_, ifudesign_, mode_, vel_, vdisp_, ell_, pa_, reff_, /silent, $
                     format='L,L,A,F,F,F,F,F', comment='#'
            nt = n_elements(plate_)
            if index ge nt then $
                message, inptbl+' does not have a row index '+MDAP_STC(index, /integer)

            plate = plate_[index]
            ifudesign = ifudesign_[index]
            mode = mode_[index]
            velocity_initial_guess = vel_[index]
            velocity_dispersion_initial_guess = vdisp_[index]
            ell = ell_[index]
            pa = pa_[index]
            Reff = reff_[index]

        endif
        
        if n_elements(par) ne 0 then begin
            if file_test(par) eq 0 then $
                message, 'Cannot open '+par+'!'
            
            params = YANNY_READONE(par)
            plate = params.plate
            ifudesign = params.ifudesign
            mode = params.mode
            velocity_initial_guess = params.vel
            velocity_dispersion_initial_guess = params.vdisp
            ell = params.ell
            pa = params.pa
            Reff = params.reff
        endif

        ; Define the DRP input path
        if n_elements(drppath) eq 0 then $
            drppath = getenv('MANGA_SPECTRO_REDUX') + '/' + getenv('MANGADRP_VER') + '/' $
                      + MDAP_STC(plate, /integer) + '/stack'

        ; Define the DAP output path
        if n_elements(dappath) eq 0 then $
            dappath = getenv('MANGA_SPECTRO_ANALYSIS') + '/' + getenv('MANGADRP_VER') + '/' $
                      + getenv('MANGADAP_VER')

        
        root_name = drppath + '/manga-' + MDAP_STC(plate, /integer) + '-' + $
                    MDAP_STC(ifudesign, /integer) + '-LOG' + mode

        output_root_dir = dappath + '/' + MDAP_STC(plate, /integer) + '/' + $
                          MDAP_STC(ifudesign, /integer)

END



            
