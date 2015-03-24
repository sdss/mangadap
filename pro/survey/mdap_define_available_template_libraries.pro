;===============================================================================
; STELLAR TEMPLATE LIBRARY FILES
;===============================================================================
;       The stellar template library files should be a list of 1D fits files,
;       and be associated with one of the following library keys:
;
;                                Spectral
;                 KEY    resolution (ang)       
;       -------------    ----------------
;           M11-MARCS                2.73
;          M11-STELIB                3.40
;          M11-ELODIE                0.55
;           M11-MILES                2.54
;               MILES                2.50
;              STELIB                3.40
;             MIUSCAT                2.51
;
;       See dapsrc/pro/fileio/mdap_read_template_library.pro.  When
;       adding additional libraries, you must add the keyword to the
;       MDAP_GET_TEMPLATE_RESOLUTION() function with the appropriate
;       resolution.
;      
;       TODO: Include the resolution as a keyword in the fits file
;       headers of the library spectra?
;
;       The fits files must have CRVAL1, CRPIX1, and CDELT1 keywords used to
;       define the wavelength coordinates of each pixel:
;
;               pix = (findgen(npixels+1))[1:*]
;               wave = (pix - CRPIX1) * CDELT1 + CRVAL1
;
;       The reference frame of the template wavelengths must also be defined as
;       either vacuum or air via the tpl_vacuum_wave vector.  Set the value to
;       be 1(true) if the wavelengths are in vacuum, 0(false) otherwise.  It is
;       expected that the DRP spectra are in vacuum wavelengths.  The DAP will
;       therefore use the IDL routine AIRTOVAC to convert the template
;       wavelengths to vacuum if tpl_vacuum_wave = 0.
;
;===============================================================================
;
; REVISION HISTORY:
;       ?? ??? 2015: Original implementation by K. Westfall (KBW)
;       24 Mar 2015: (KBW) Added MIUSCAT library
;
; dapsrc is an optional input to define the DAP source path instead of
; using environmental varaibles.

PRO MDAP_DEFINE_AVAILABLE_TEMPLATE_LIBRARIES, $
        tpl_library_keys, template_libraries, tpl_vacuum_wave, dapsrc=dapsrc

        ; Define the DAP source path
        if n_elements(dapsrc) eq 0 then $
            dapsrc = getenv('MANGADAP_DIR')

        ;-----------------------------------------------------------------------
        ; Define the set of template libraries.  The format expected for these
        ; files is described above.
        ntpl_libraries = 7
        tpl_library_keys = strarr(ntpl_libraries)
        template_libraries = strarr(ntpl_libraries)
        tpl_vacuum_wave = intarr(ntpl_libraries)

        tpl_library_keys[0] = 'M11-MARCS'
        template_libraries[0] = dapsrc+'/external/templates/m11_marcs/*_s.fits'
        tpl_vacuum_wave[0] = 0

        tpl_library_keys[1] = 'M11-STELIB'
        template_libraries[1] = dapsrc+'/external/templates/m11_stelib/*_s.fits'
        tpl_vacuum_wave[1] = 0

        tpl_library_keys[2] = 'M11-ELODIE'
        template_libraries[2] = dapsrc+'/external/templates/m11_elodie/*.fits'
        tpl_vacuum_wave[2] = 0

        tpl_library_keys[3] = 'M11-MILES'
        template_libraries[3] = dapsrc+'/external/templates/m11_miles/*.fits'
        tpl_vacuum_wave[3] = 0

        tpl_library_keys[4] = 'MILES'
        template_libraries[4] = dapsrc+'/external/templates/miles/*.fits'
        ; TODO: Unknown if this library is in vacuum or in air
        tpl_vacuum_wave[4] = 0

        tpl_library_keys[5] = 'STELIB'
        template_libraries[5] = dapsrc+'/external/templates/stelib/*.fits'
        ; TODO: Unknown if this library is in vacuum or in air
        tpl_vacuum_wave[5] = 0

        tpl_library_keys[6] = 'MIUSCAT'
        template_libraries[6] = dapsrc+'/external/templates/miuscat/*.fits'
        ; TODO: Unknown if this library is in vacuum or in air
        tpl_vacuum_wave[6] = 0

END

