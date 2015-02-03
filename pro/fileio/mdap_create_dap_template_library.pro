;+
; NAME:
;       MDAP_CREATE_DAP_TEMPLATE_LIBRARY
;
; PURPOSE:
;       Read in a template library and manipulate it such that it can be
;       directly used with the analysis procedures of the DAP.  This is
;       primarily a wrapper for a number of other DAP routines.
;
; CALLING SEQUENCE:
;       MDAP_CREATE_DAP_TEMPLATE_LIBRARY, ofile, library_key, template_libraries, tpl_vacuum_wave, $
;                                         velocity_offset, wave, sres, velscale, quiet=quiet
;
; INPUTS:
;       ofile string
;               Output file for the manipulated template spectra.
;       
;       library_key string
;               Keyword used to identify the template library.  See
;               MDAP_READ_TEMPLATE_LIBRARY.
;       
;       template_libraries string
;               String identifier for the template library files.
;       
;       tpl_vacuum_wave integer
;               Flag that the template library is (1) or is not (0) at vacuum
;               wavelengths.
;
;       velocity_offset double
;               Offset to apply when matching the spectral resolution of the
;               template library to the provided spectral resolution.
;
;       wave dblarr[C]
;               Wavelength vector used to trim the spectral range of the library
;               and to match the spectral resolution to the provided input.
;
;       sres dblarr[C]
;               Spectral resolution (R=lamda/delta lamba) to which to match the
;               template library spectral resolution.
;
;       velscale double
;               Velocity scale of the pixel used when resampling the resolution
;               matched template library.
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;       /quiet
;               Suppress output to the screen
;
; OUTPUT:
;       The DAP-ready template library is written to disk with the name 'ofile'.
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
;       04 Nov 2014: (KBW) Original implementation (pulled from MANGA_DAP)
;-
;------------------------------------------------------------------------------

PRO MDAP_CREATE_DAP_TEMPLATE_LIBRARY, $
                ofile, library_key, template_libraries, tpl_vacuum_wave, velocity_offset, wave, $
                sres, velscale, quiet=quiet, version=version

        version_module = '0.9'                  ; Version number
        if n_elements(version) ne 0 then begin  ; If version is defined
            version = version_module            ; ... set it to the module version
            return                              ; ... and return without doing anything
        endif

        ; Read the template library
        MDAP_READ_TEMPLATE_LIBRARY, library_key, template_libraries, tpl_flux, tpl_ivar, tpl_mask, $
                                    tpl_wave, tpl_sres

        ; Convert the wavelengths to vacuum if necessary
        ; TODO: Does not account for the effect on the spectral resolution
        if tpl_vacuum_wave eq 0 then $
            AIRTOVAC, tpl_wave

        ; Need to limit template wavelength range before running the resolution
        ; matching procedure to avoid extrapolation errors.  As with the
        ; matching, this trimming accounts for the redshift of the galaxy.

        c=299792.458d                                   ; Speed of light in km/s
        redshift = velocity_offset / c                  ; Redshift
        if ~keyword_set(quiet) then begin
            print, 'Trimming templates to ' + MDAP_STC(min(wave/(1.+redshift))) + $
                    'A to ' + MDAP_STC(max(wave/(1.+redshift))) + 'A'
        endif
        MDAP_TRIM_SPECTRAL_RANGE, tpl_flux, tpl_ivar, tpl_mask, tpl_sres, tpl_wave, $
                                  ([min(wave/(1.+redshift)), max(wave/(1.+redshift))])

        ; Match the resolution of the templates to the galaxy data accounting
        ; for the redshift of the galaxy.
        MDAP_MATCH_SPECTRAL_RESOLUTION, tpl_flux, tpl_ivar, tpl_mask, tpl_wave, tpl_sres, $
                                        wave/(1.+redshift), sres, tpl_soff, /no_offset

        ; Resample the templates to match the object spectra: On input, tpl_wave
        ; expected to have the same size as tpl_flux (wavelength defined
        ; separately for each spectrum); on output, tpl_wave and tpl_sres are
        ; single vectors common to ALL template spectra.  Sampling is forced to
        ; match galaxy sampling, via velscale
        MDAP_RESAMPLE_TEMPLATES, tpl_flux, tpl_ivar, tpl_mask, tpl_wave, tpl_sres, velScale, $
                                 /reform_sres

        ; Normalize the templates and save the normalization
        MDAP_NORMALIZE_TEMPLATES, tpl_flux, tpl_ivar, tpl_mask, tpl_flux_norm

        ; Write the template data to a fits file
        MDAP_WRITE_RESAMPLED_TEMPLATES, ofile, library_key, redshift, tpl_flux_norm, $
                                        version_module, tpl_wave, tpl_flux, tpl_ivar, tpl_mask, $
                                        tpl_sres, tpl_soff
END


