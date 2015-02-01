;+
; NAME:
;       MDAP_DRP_DATA_BLOCK
;
; PURPOSE:
;   MAIN BLOCK FOR MANGA_DAP:
;       - Read the DRP-based data arrays
;
; CALLING SEQUENCE:
;       MDAP_DRP_DATA_BLOCK, manga_dap_version, datacube_name, header, wave, flux, ivar, mask, $
;                            sres, spaxel_dx, spaxel_dy, bskyx, bskyy, type=type, $
;                            instrumental_fwhm_file=instrumental_fwhm_file, nolog=nolog, 
;                            log_file_unit=log_file_unit
;
; INPUTS:

;       manga_dap_version MaNGADAPVersion
;               Structure used to keep track of various
;               version-controlled procedures in the DAP.
;
;       datacube_name string
;               Name of the DRP (RSS or CUBE) fits file to read.
;      
;
; OPTIONAL INPUTS:
;       type string
;               Type of DRP fits file to read, either 'RSS' or 'CUBE'.
;               If provided, it's checked.  If not provided, it's
;               determined.
;
;       instrumental_fwhm_file string
;
;               File with the instrumental resolution contained in four
;               colums: wavelength, resolution, FWHM in angstroms, FWHM
;               in km/s.  This OVERRIDES the existing spectral resolution
;               information provided in the DRP fits file.
;
;       log_file_unit LUN
;               File unit pointing to the log file
;
; OPTIONAL KEYWORDS:
;       /nolog
;               Suppress output to the log file.
;
; OUTPUT:
;       header strarr[]
;               Fits header
;            
;       wave dblarr[T]
;               Wavelength of each spectral channel T.
;       
;       flux dblarr[N][T]
;               Galaxy spectra as produced by MDAP_READ_DRP_FITS.
;
;       ivar dblarr[N][T]
;               Inverse variance of the flux
;
;       mask dblarr[N][T]
;               Bad pixel mask.
;
;       sres dblarr[T]
;               Spectral resolution at each wavelength channel T.
;
;       spaxel_dx double
;               Spaxel size in on-sky X (arcsec).
;
;       spaxel_dy double
;               Spaxel size in on-sky Y (arcsec).
;
;       bskyx dblarr[N]
;               Array containing the fiducial x coordinates in
;               arcseconds (0 is the center of the field of view) for
;               each spectrum.  See MDAP_FIDUCIAL_BIN_XY.
;
;       bskyy dblarr[N]
;               Array containing the fiducial y coordinates in
;               arcseconds (0 is the center of the field of view) for
;               each spectrum. See MDAP_FIDUCIAL_BIN_XY.
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
;       MDAP_READ_DRP_FITS
;       MDAP_FIDUCIAL_BIN_XY
;       MDAP_GET_SPAXEL_SIZE
;       SXADDPAR
;       MDAP_STC()
;       READCOL
;
; REVISION HISTORY:
;       01 Feb 2015: (KBW) Pulled from manga_dap.pro
;-
;------------------------------------------------------------------------------

PRO MDAP_DRP_DATA_BLOCK, $
                manga_dap_version, datacube_name, header, wave, flux, ivar, mask, sres, spaxel_dx, $
                spaxel_dy, bskyx, bskyy, type=type, instrumental_fwhm_file=instrumental_fwhm_file, $
                nolog=nolog, log_file_unit=log_file_unit

        ;---------------------------------------------------------------
        ; Version control
        MDAP_READ_DRP_FITS, version=manga_dap_version.read_drp_fits
        MDAP_FIDUCIAL_BIN_XY, version=manga_dap_version.fiducial_bin_xy

        ; TODO: Version required for this should depend on the version of DRP
        ; output being read

        ; TODO: Put in the conditions for when this block is executed:
        ;       - If any of the plans have the overwrite flag flipped
        ;       - If any of the plans do not have a previously created output
        ;         file with the binned spectra

        ; TODO: Allow to read a single spectrum and skip binning step
        MDAP_READ_DRP_FITS, datacube_name, header, flux, ivar, mask, wave, sres, skyx, skyy, $
                            type=type, unit=unit

        mask[*,*] = 0.                      ; TODO: Unmask everything for now

        ; Get the spaxel size
        MDAP_GET_SPAXEL_SIZE, header, spaxel_dx, spaxel_dy, type=type, unit=unit

        ; Set a fiducial set of coordinates to use for each spectrum
        MDAP_FIDUCIAL_BIN_XY, wave, skyx, skyy, bskyx, bskyy, type=type

        ; Add/update the version number for DAP read version
        SXADDPAR, header, 'VDAPREAD', manga_dap_version.read_drp_fits, 'mdap_read_drp_fits version'
        SXADDPAR, header, 'VDAPFXY', manga_dap_version.fiducial_bin_xy, $
                  'mdap_fiducial_bin_xy version'

        ; Print information to the log file
        if ~keyword_set(nolog) then begin
            printf, log_file_unit, '[INFO] Read DRP data using version: ' + $
                    manga_dap_version.read_drp_fits
            printf, log_file_unit, '[INFO] Input data type: '+type
            sz = size(flux)
            printf, log_file_unit, '[INFO] Total number of spectra: '+MDAP_STC(sz[1])
            printf, log_file_unit, '[INFO] Number of spectral channels: '+MDAP_STC(sz[2])
            printf, log_file_unit, '[INFO] Fiducial bin centers calcualted using version: ' + $
                    manga_dap_version.fiducial_bin_xy
        endif

        ;###########################################################
        ; TODO: This block of code is essentially obsolete.  We
        ; expect to always read the spectral resolution from the
        ; header of the DRP-produced fits file.

        ; Change the spectral resolution according to an input file
        ; Columns must be:
        ;   1. Wavelength in angstroms
        ;   2. Resolution (lamda/delta lambda)
        ;   3. delta lamga (FWHM) in angstroms
        ;   4. delta lamga (FWHM) in km/s
        if n_elements(instrumental_fwhm_file) ne 0 then begin
            if file_test(instrumental_fwhm_file) eq 0 then begin
                print, 'Unable to read '+instrumental_fwhm_file+'!  Using fits extention or ' $
                       + 'default.'
            endif
            print, 'Reading '+instrumental_fwhm_file
            READCOL, instrumental_fwhm_file, ww, r, fwhm_ang, fwhm_kms, /silent
        endif
            
        ; Initialize the instrumental resolution.
        if n_elements(fwhm_ang) ne 0 then begin             ; use data from input file
            if ~keyword_set(nolog) then begin
                printf, log_file_unit, '[INFO] Spectral resolution changed using data in: '+ $
                        instrumental_fwhm_file
            endif
            sres=interpol(fwhm_ang,r,wave)      ; Interpolate to the object wavelengths
        endif 
        
        if n_elements(sres) eq 0 then begin     ; sres not read or  available from fits file
            if ~keyword_set(nolog) then begin
                printf, log_file_unit, '[INFO] Spectral resolution unavailable, assume R=2000'
            endif
            sres = make_array(n_elements(wave), /double, value=2000.0d)     ; Set R=2000.
        endif
        ;###########################################################

END

