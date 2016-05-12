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
;     M11-STELIB-ZSOL                3.40
;          M11-ELODIE                0.55
;           M11-MILES                2.54
;               MILES                2.50
;           MILES-AVG                2.50
;          MILES-THIN                2.50
;              STELIB                3.40
;             MIUSCAT                2.51
;        MIUSCAT-THIN                2.51
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
;       07 Oct 2015: (KBW) Added M11-STELIB-ZSOL library
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
        ntpl_libraries = 11
        tpl_library_keys = strarr(ntpl_libraries)
        template_libraries = strarr(ntpl_libraries)
        tpl_vacuum_wave = intarr(ntpl_libraries)

        i = 0
        tpl_library_keys[i] = 'M11-MARCS'
        template_libraries[i] = dapsrc+'/data/stellar_templates/m11_marcs/*_s.fits'
        tpl_vacuum_wave[i] = 0
        i++

        tpl_library_keys[i] = 'M11-STELIB'
        template_libraries[i] = dapsrc+'/data/stellar_templates/m11_stelib/*_s.fits'
        tpl_vacuum_wave[i] = 0
        i++

        tpl_library_keys[i] = 'M11-STELIB-ZSOL'
        template_libraries[i] = dapsrc+'/data/stellar_templates/m11_stelib_zsol/*_s.fits'
        tpl_vacuum_wave[i] = 0
        i++

        tpl_library_keys[i] = 'M11-ELODIE'
        template_libraries[i] = dapsrc+'/data/stellar_templates/m11_elodie/*.fits'
        tpl_vacuum_wave[i] = 0
        i++

        tpl_library_keys[i] = 'M11-MILES'
        template_libraries[i] = dapsrc+'/data/stellar_templates/m11_miles/*.fits'
        tpl_vacuum_wave[i] = 0
        i++

        tpl_library_keys[i] = 'MILES'
        template_libraries[i] = dapsrc+'/data/stellar_templates/miles/*.fits'
        ; TODO: Unknown if this library is in vacuum or in air
        tpl_vacuum_wave[i] = 0
        i++

        tpl_library_keys[i] = 'MILES-AVG'
        template_libraries[i] = dapsrc+'/data/stellar_templates/miles_avg/*.fits'
        ; TODO: Unknown if this library is in vacuum or in air
        tpl_vacuum_wave[i] = 0
        i++

        tpl_library_keys[i] = 'MILES-THIN'
        template_libraries[i] = dapsrc+'/data/stellar_templates/miles_thin/*.fits'
        ; TODO: Unknown if this library is in vacuum or in air
        tpl_vacuum_wave[i] = 0
        i++

        tpl_library_keys[i] = 'STELIB'
        template_libraries[i] = dapsrc+'/data/stellar_templates/stelib/*.fits'
        ; TODO: Unknown if this library is in vacuum or in air
        tpl_vacuum_wave[i] = 0
        i++

        tpl_library_keys[i] = 'MIUSCAT'
        template_libraries[i] = dapsrc+'/data/stellar_templates/miuscat/*.fits'
        ; TODO: Unknown if this library is in vacuum or in air
        tpl_vacuum_wave[i] = 0
        i++

        tpl_library_keys[i] = 'MIUSCAT-THIN'
        template_libraries[i] = dapsrc+'/data/stellar_templates/miuscat_thin/*.fits'
        ; TODO: Unknown if this library is in vacuum or in air
        tpl_vacuum_wave[i] = 0
        i++

END

