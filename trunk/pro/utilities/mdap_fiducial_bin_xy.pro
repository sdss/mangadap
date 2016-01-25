;+
; NAME:
;       MDAP_FIDUCIAL_BIN_XY
;
; PURPOSE:
;       Compute a single set of coordinates for each spectrum. For CUBE
;       files, coordinates are expected to be independent of wavelength
;       so that the coordinates can just be set to the first value.  For
;       the RSS files, fidicial coorinates are set at the central
;       wavelength using MDAP_SELECT_COO_RSS.  Previous usage was to set
;       the fiducial coordinates to the wavelength of the central pixel;
;       these are not necessarily the same thing given that spectra are
;       typically logarithmically binned in wavelength.
;
; CALLING SEQUENCE:
;       MDAP_FIDUCIAL_BIN_XY, wave, skyx, skyy, bskyx, bskyy, version=version
;
; INPUTS:
;       wave dblarr[T]
;               Wavelength coordinate of each spectral channel for ALL
;               spectra.
;
;       skyx dblarr[N][T]
;               Spatial X coordinate in arcseconds for each spectral channel,
;               with 0 at the center of the field of view.
;
;       skyy dblarr[N][T]
;               Spatial Y coordinate in arcseconds for each spectral channel,
;               with 0 at the center of the field of view.
;
; OPTIONAL INPUTS:
;       type string
;               Type of the DRP file.  Should be CUBE or RSS.  In
;               practice, if type is not CUBE (either by a test of
;               equality or because type is not provided), then the type
;               is assumed to be RSS.
;
; OPTIONAL KEYWORDS:
;
; OUTPUT:
;       bskyx dblarr[N]
;               Spatial X coordinate in arcseconds to use globally for each
;               spectrum.
;
;       bskyy dblarr[N]
;               Spatial Y coordinate in arcseconds to use globally for each
;               spectrum.
;
; OPTIONAL OUTPUT:
;       version string
;               Module version.  If requested, the module is not executed and
;               only version flag is returned.
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
;       15 Sep 2014: (KBW) Original implementation
;       12 Nov 2014: (KBW) Added version
;       04 Dec 2014: (KBW) Add wave, type.  Wave needed for RSS.  For
;                          CUBE files, coordinates are expected to be
;                          independent of wavelength so that the
;                          coordinates can just be set to the first
;                          value.  For the RSS files, fidicial
;                          coorinates now set at the central wavelength
;                          using MDAP_SELECT_COO_RSS.
;-
;------------------------------------------------------------------------------

PRO MDAP_FIDUCIAL_BIN_XY, $
                wave, skyx, skyy, bskyx, bskyy, type=type, version=version

        version_module = '0.2'                  ; Version number
        if n_elements(version) ne 0 then begin  ; If version is defined
            version = version_module            ; ... set it to the module version
            return                              ; ... and return without doing anything
        endif

        if type eq 'CUBE' then begin
            bskyx = skyx[*,0]
            bskyy = skyy[*,0]
            return
        endif else begin
            targwave = (min(wave)+max(wave))/2.0
            MDAP_SELECT_COO_RSS, wave, skyx, skyy, targwave, bskyx, bskyy
        endelse

END

;PRO MDAP_FIDUCIAL_BIN_XY, $
;                skyx, skyy, bskyx, bskyy, version=version
;
;        version_module = '0.1'                  ; Version number
;        if n_elements(version) ne 0 then begin  ; If version is defined
;            version = version_module            ; ... set it to the module version
;            return                              ; ... and return without doing anything
;        endif
;
;        sz=size(skyx)                           ; Size of the skyx array
;        ns=sz[1]                                ; Number of spectra
;        nc=sz[2]                                ; Number of spectral channels
;
;        bskyx=dblarr(ns)                        ; Initialize the arrays
;        bskyy=dblarr(ns)
;
;        for i=0,ns-1 do begin                   ; Set fiducial coo to coo at central channel
;            bskyx[i] = skyx[i,nc/2]
;            bskyy[i] = skyy[i,nc/2]
;        endfor
;
;END


