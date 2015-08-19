;+
; NAME:
;       MDAP_SPECTRAL_INDEX_MEASUREMENTS
;
; PURPOSE:
;       Measure absorption-line indices using MDAP_GET_SPECTRAL_INDEX.
;       The indices are measured for the data, as well as the broadened
;       and unbroadened version of the optimal template determined by
;       the stellar continuum fit.  When the spectral range of the
;       templates are sufficient, these template-based measurements are
;       used to correct the index measurement from the galaxy spectrum
;       for the velocity dispersion.  The corrections are different if
;       the index is measured in angstroms or magnitudes.  For
;       magnitudes:
;
;           losvd_correction = otpl_index - botpl_index
;           spectral_index = spectral_index_raw + losvd_correction
;
;       For angstroms:
;
;           losvd_correction = otpl_index / botpl_index
;           spectral_index = spectral_index_raw * losvd_correction
;
;
; CALLING SEQUENCE:
;       MDAP_SPECTRAL_INDEX_MEASUREMENTS, abs_par, wave, flux, ivar, mask, otpl_flux, otpl_mask, $
;                                         botpl_flux, botpl_mask, stellar_kinematics, $
;                                         spectral_index_omitted, spectral_index_blue_cont, $
;                                         spectral_index_blue_cerr, spectral_index_red_cont, $
;                                         spectral_index_red_cerr, spectral_index_raw, $
;                                         spectral_index_rerr, spectral_index_rperr, $
;                                         spectral_index_otpl_blue_cont, $
;                                         spectral_index_otpl_red_cont, spectral_index_otpl_index, $
;                                         spectral_index_botpl_blue_cont, $
;                                         spectral_index_botpl_red_cont, $
;                                         spectral_index_botpl_index, spectral_index_corr, $
;                                         spectral_index_cerr, spectral_index_cperr, $
;                                         version=version, /dbg, /quiet
; 
; INPUTS:
;       abs_par SpectralIndex[I]
;               Array of SpectralIndex structures.  See
;               MDAP_READ_ABSORPTION_LINE_PARAMETERS.
;
;       wave dblarr[C]
;               Wavelength at each of C spectral channels of the input spectra.
; 
;       flux dblarr[B][C]
;               Flux in the B observed spectra with C spectral channels each.
;               These are the spectra to use for the spectral index meaurements.
;               The spectra must be prepared for the measurements BEFORE being
;               passed to this function.  See
;               MDAP_SPECTRAL_INDEX_MEASUREMENTS_SPECTRA.
;
;       ivar dblarr[B][C]
;               Inverse variance in the flux of the observed spectra.
;       
;       mask dblarr[B][C]
;               Bad pixel mask (0-good, 1-bad) for the object spectra.
;
;       otpl_flux dblarr[B][C]
;       otpl_mask dblarr[B][C]
;               Optimal template spectrum, resolution-matched to the
;               spectral-index system, and mask.
;
;       botpl_flux dblarr[B][C]
;       botpl_mask dblarr[B][C]
;               Optimal template spectrum, resolution-matched to the
;               spectral-index system, and broadened by the best-fitting
;               LOSVD for each of the B spectra, and its bad pixel mask.
;
;       stellar_kinematics dblarr[B][M]
;               Best-fitting stellar kinematics with M moments for each of the B
;               object spectra.
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;       /dbg
;               Only fit the set of indices to the first spectrum as a debugging
;               test.
;
;       /quiet
;               Suppress output to STDOUT
;
; OUTPUT:
;       spectral_index_omitted intarr[B][I]
;               Flag that the index was (1) or was not (0) omitted because it
;               fell outside the spectral range of the observations.
;
;       spectral_index_blue_cont dblarr[B][I]
;       spectral_index_blue_cerr dblarr[B][I]
;       spectral_index_red_cont dblarr[B][I]
;       spectral_index_red_cerr dblarr[B][I]
;               The blue and red pseudo-continuum measurements and their
;               propagated error for each of the B binned spectra and
;               each of the I spectral indices.
;
;       spectral_index_raw dblarr[B][I]
;       spectral_index_rerr dblarr[B][I]
;       spectral_index_rperr dblarr[B][I]
;               Spectral index measured directly from the galaxy spectra
;               following Worthey et al. (1994), the nominal error
;               following Cardiel et al. (1998) (err), and an estimate
;               of error from nominal error propagation (perr).  The
;               index is defined either in angstroms or magnitudes as
;               specified by abs_par.
;
;       spectral_index_otpl_blue_cont dblarr[B][I]
;       spectral_index_otpl_red_cont dblarr[B][I]
;       spectral_index_otpl_index dblarr[B][I]
;               The blue and red pseudo-continuum measurements and the
;               index value for the (unbroadened) optimal template for
;               each of the B binned spectra and each of the I spectral
;               indices.
;
;       spectral_index_botpl_blue_cont dblarr[B][I]
;       spectral_index_botpl_red_cont dblarr[B][I]
;       spectral_index_botpl_index dblarr[B][I]
;               The blue and red pseudo-continuum measurements and the
;               index value for the broadened optimal template for each
;               of the B binned spectra and each of the I spectral
;               indices.
;
;       spectral_index_corr dblarr[B][I]
;       spectral_index_cerr dblarr[B][I]
;       spectral_index_cperr dblarr[B][I]
;               Spectral indices for the galaxy spectra after correcting
;               for Doppler broadening using the correction derived by
;               comparing the index measurements from the broadened and
;               unbroadened versions of the optimal template.
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
; TODO:
;       - Allow for a tag that will ignore the MDAP_REVERT_PIXEL_KINEMATICS
;         call(s).
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;       31 Oct 2014: Copied from L. Coccato's version (PDAP v0_8) by K.
;                    Westfall (KBW)
;       12 Aug 2015: (KBW) Change default values in output arrays to
;                          -9999.0.  Added a LARGE number of additional
;                          outputs: pseudo-continua, raw indices,
;                          propagated errors, and correction details.
;-
;----------------------------------------------------------------------

PRO MDAP_SPECTRAL_INDEX_MEASUREMENTS, $
                abs_par, wave, flux, ivar, mask, otpl_flux, otpl_mask, botpl_flux, botpl_mask, $
                stellar_kinematics, spectral_index_omitted, spectral_index_blue_cont, $
                spectral_index_blue_cerr, spectral_index_red_cont, spectral_index_red_cerr, $
                spectral_index_raw, spectral_index_rerr, spectral_index_rperr, $
                spectral_index_otpl_blue_cont, spectral_index_otpl_red_cont, $
                spectral_index_otpl_index, spectral_index_botpl_blue_cont, $
                spectral_index_botpl_red_cont, spectral_index_botpl_index, spectral_index_corr, $
                spectral_index_cerr, spectral_index_cperr, version=version, dbg=dbg, quiet=quiet

        ; Just return the version if requested
        version_module = '1.1'                  ; Version number
        if n_elements(version) ne 0 then begin  ; If version is defined
            version = version_module            ; ... set it to the module version
            return                              ; ... and return without doing anything
        endif

        ;---------------------------------------------------------------
        ; SETUP --------------------------------------------------------

        ; On input, the mask expected to be directly applicable to flux
        ; and botpl.  HOWEVER, it is expected to be applicable to the
        ; DEREDSHIFTED wavelengths of otpl, the optimal template.

        ; TODO: Do some checks of the input arrays?

        c=299792.458d                           ; Speed of light in km/s
        redshift = stellar_kinematics[*,0]/c    ; Redshift of each spectrum

        ; Allocate the output arrays, setting default values to -9999
        nspec = (size(flux))[1]                         ; Number of spectra
        npix = (size(flux))[2]                          ; Number of wavelength channels
        nabs = n_elements(abs_par)                      ; Number of spectral indices

        ; Dummy errors for use with otpl and botpl
        tpl_ivar = make_array((size(flux))[2], /double, value=1.0d)

        spectral_index_omitted = intarr(nspec, nabs)    ; Flag that index was omitted

        spectral_index_blue_cont = make_array(nspec, nabs, /double, value=-9999.0d)  ; Blue contin
        spectral_index_blue_cerr = make_array(nspec, nabs, /double, value=-9999.0d)  ; and error
        spectral_index_red_cont = make_array(nspec, nabs, /double, value=-9999.0d)   ; Red continuum
        spectral_index_red_cerr = make_array(nspec, nabs, /double, value=-9999.0d)   ; and error

        spectral_index_raw = make_array(nspec, nabs, /double, value=-9999.0d)    ; Raw measurement
        spectral_index_rerr = make_array(nspec, nabs, /double, value=-9999.0d)    ; error formula
        spectral_index_rperr = make_array(nspec, nabs, /double, value=-9999.0d)   ; propagated

        spectral_index_otpl_blue_cont = make_array(nspec, nabs, /double, value=-9999.0d)
        spectral_index_otpl_red_cont = make_array(nspec, nabs, /double, value=-9999.0d)
        spectral_index_otpl_index = make_array(nspec, nabs, /double, value=-9999.0d)

        spectral_index_botpl_blue_cont = make_array(nspec, nabs, /double, value=-9999.0d)
        spectral_index_botpl_red_cont = make_array(nspec, nabs, /double, value=-9999.0d)
        spectral_index_botpl_index = make_array(nspec, nabs, /double, value=-9999.0d)

        spectral_index_corr = make_array(nspec, nabs, /double, value=-9999.0d)  ; Corrected index
        spectral_index_cerr = make_array(nspec, nabs, /double, value=-9999.0d)  ; error formula
        spectral_index_cperr = make_array(nspec, nabs, /double, value=-9999.0d) ; propaged error

        ;---------------------------------------------------------------

        ;---------------------------------------------------------------
        ; Measure all spectral indices for all spectra -----------------
        if keyword_set(dbg) then $              ; Only analyze the first spectrum
            nspec = 1
        for i=0,nspec-1 do begin

;            ; Set the mask for otpl to be the DEREDSHIFTED version of the input
;            ; mask
;            otpl_mask = interpol(reform(mask[i,*]), wave, wave/(redshift[i]+1.0d))
;            indx = where(otpl_mask gt 0.5, count, complement=nindx, ncomplement=ncount)
;            if count ne 0 then $
;                otpl_mask[indx] = 1.0d
;            if ncount ne 0 then $
;                otpl_mask[nindx] = 0.0d

            ; Get the REDSHIFTED passbands for getting the index measurements of
            ; the object data
            abs_par_z = abs_par
            abs_par_z.passband = abs_par_z.passband*(redshift[i]+1.0d)
            abs_par_z.blue_cont = abs_par_z.blue_cont*(redshift[i]+1.0d)
            abs_par_z.red_cont = abs_par_z.red_cont*(redshift[i]+1.0d)

            ; Measure each index
            for j=0,nabs-1 do begin
;               print, 'index:', j+1
;
;               print, 'def:'
;               print, abs_par_z[j].passband
;               print, abs_par_z[j].blue_cont
;               print, abs_par_z[j].red_cont

                ; Measure the index based on the object data (use redshifted windows)
                err = 0
                MDAP_GET_SPECTRAL_INDEX, wave, reform(flux[i,*]), reform(ivar[i,*]), $
                                         reform(mask[i,*]), abs_par_z[j], blue_cont, $
                                         blue_cont_err, red_cont, red_cont_err, $
                                         equiv_width, equiv_width_err, equiv_width_perr, $
                                         index_mag, index_mag_err, index_mag_perr, err=err, $
                                         /geometric, /quiet; quiet=quiet

                if err eq 1 then begin                  ; Measurement returned error!
                    spectral_index_omitted[i,j] = 1
                    continue                            ; Continue without saving anything
                endif

                ; Save the raw results
                spectral_index_blue_cont[i,j] = blue_cont
                spectral_index_blue_cerr[i,j] = blue_cont_err
                spectral_index_red_cont[i,j] = red_cont
                spectral_index_red_cerr[i,j] = red_cont_err
                if abs_par[j].units eq 'ang' then begin
                    spectral_index_raw[i,j] = equiv_width
                    spectral_index_rerr[i,j] = equiv_width_err
                    spectral_index_rperr[i,j] = equiv_width_perr
                endif
                if abs_par[j].units eq 'mag' then begin
                    spectral_index_raw[i,j] = index_mag
                    spectral_index_rerr[i,j] = index_mag_err
                    spectral_index_rperr[i,j] = index_mag_perr
                endif

                ; Measure the index based on the optimal template (use rest-frame windows)
                MDAP_GET_SPECTRAL_INDEX, wave, reform(otpl_flux[i,*]), tpl_ivar, $
                                         reform(otpl_mask[i,*]), abs_par[j], blue_cont, $
                                         blue_cont_err, red_cont, red_cont_err, $
                                         equiv_width, equiv_width_err, equiv_width_perr, $
                                         index_mag, index_mag_err, index_mag_perr, err=err, $
                                         /geometric, /quiet; quiet=quiet

                if err eq 1 then begin                  ; Measurement returned error!
                    spectral_index_omitted[i,j] = 1
                    continue                            ; Continue without saving template results
                endif

                ; Save the raw results
                spectral_index_otpl_blue_cont[i,j] = blue_cont
                spectral_index_otpl_red_cont[i,j] = red_cont
                if abs_par[j].units eq 'ang' then $
                    spectral_index_otpl_index[i,j] = equiv_width
                if abs_par[j].units eq 'mag' then $
                    spectral_index_otpl_index[i,j] = index_mag

                ; Measure the index based on the broadened optimal template (use redshifted windows)
                MDAP_GET_SPECTRAL_INDEX, wave, reform(botpl_flux[i,*]), tpl_ivar, $
                                         reform(botpl_mask[i,*]), abs_par_z[j], blue_cont, $
                                         blue_cont_err, red_cont, red_cont_err, $
                                         equiv_width, equiv_width_err, equiv_width_perr, $
                                         index_mag, index_mag_err, index_mag_perr, err=err, $
                                         /geometric, /quiet; quiet=quiet

                if err eq 1 then begin                  ; Measurement returned error!
                    spectral_index_omitted[i,j] = 1
                    continue                            ; Continue without saving template results
                endif

                ; Save the raw results
                spectral_index_botpl_blue_cont[i,j] = blue_cont
                spectral_index_botpl_red_cont[i,j] = red_cont
                if abs_par[j].units eq 'ang' then $
                    spectral_index_botpl_index[i,j] = equiv_width
                if abs_par[j].units eq 'mag' then $
                    spectral_index_botpl_index[i,j] = index_mag

                ; Correct the index measurements using the difference between the
                ; measurements performed on the optimal template and its
                ; broadened version.  The correction depends on the desired
                ; units.
                if abs_par[j].units eq 'ang' then begin
                    losvd_correction = spectral_index_otpl_index[i,j] $
                                       / spectral_index_botpl_index[i,j]
                    spectral_index_corr[i,j] = spectral_index_raw[i,j] * losvd_correction
                    spectral_index_cerr[i,j] = spectral_index_rerr[i,j] * abs(losvd_correction)
                    spectral_index_cperr[i,j] = spectral_index_rperr[i,j] * abs(losvd_correction)
                endif
                if abs_par[j].units eq 'mag' then begin
                    losvd_correction = spectral_index_otpl_index[i,j] $
                                       - spectral_index_botpl_index[i,j]
                    spectral_index_corr[i,j] = spectral_index_raw[i,j] + losvd_correction
                    spectral_index_cerr[i,j] = spectral_index_rerr[i,j]
                    spectral_index_cperr[i,j] = spectral_index_rperr[i,j]
                endif

;               print, abs_line_indx[i,j], abs_line_indx_err[i,j], abs_line_indx_otpl[i,j], $
;                      abs_line_indx_botpl[i,j]
;               stop

            endfor ; End loop over indices

        endfor ; End loop over spectra
        ;---------------------------------------------------------------

END

