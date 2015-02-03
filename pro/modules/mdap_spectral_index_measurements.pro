;+
; NAME:
;       MDAP_SPECTRAL_INDEX_MEASUREMENTS
;
; PURPOSE:
;       Measure absorption-line indices.
;
;       TODO: Improve this description.
;
; CALLING SEQUENCE:
;       MDAP_SPECTRAL_INDEX_MEASUREMENTS, tpl_lib_fits, weights, stellar_kinematics, abs_par, $
;                                         obj_wave, obj_flux, obj_ivar, obj_mask, $
;                                         abs_line_indx_omitted, abs_line_indx, abs_line_indx_err, $
;                                         abs_line_indx_otpl, abs_line_indx_botpl, $
;                                         moments=moments, version=version, oversample=oversample, $
;                                         dbg=dbg
;
; INPUTS:
;       stellar_kinematics dblarr[B][M]
;               Best-fitting stellar kinematics with M moments for each of the B
;               object spectra.
;
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
; OPTIONAL INPUTS:
;
;       (REMOVED!) output_file_root string
;               Root name for PS file output.
;
; OPTIONAL KEYWORDS:
;       (REMOVED!) /plot
;               Create the output PS plots.
;
;       /dbg
;               Only fit the set of indices to the first spectrum as a debugging
;               test.
;
; OUTPUT:
;       abs_line_omitted intarr[B][I]
;               Flag that the index was (1) or was not (0) omitted because it
;               fell outside the spectral range of the observations.
;
;       abs_line_indx dblarr[B][I]
;               Absorption-line-index measurements for the I absorption lines in
;               the B *object* spectra, corrected for intrinsic broadening.
;
;       abs_line_indx_err dblarr[B][I]
;               Error in each of the index measurements.
;
;       abs_line_indx_otpl dblarr[B][I]
;               Absorption-line-index measurements based on the B
;               optimal_template spectra.
;               
;       abs_line_indx_botpl dblarr[B][I]
;               Absorption-line-index measurements based on the B
;               losvd_optimal_template spectra.
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
;       31 Oct 2014: (KBW) Copied from L. Coccato's version (PDAP v0_8)
;-
;------------------------------------------------------------------------------

PRO MDAP_SPECTRAL_INDEX_MEASUREMENTS, $
                abs_par, wave, flux, ivar, mask, otpl, botpl, stellar_kinematics, $
                abs_line_indx_omitted, abs_line_indx, abs_line_indx_err, abs_line_indx_otpl, $
                abs_line_indx_botpl, version=version, dbg=dbg

        ; On input, the mask expected to be directly applicable to flux and
        ; botpl.  HOWEVER, it is expected to be applicable to the DEREDSHIFTED
        ; wavelengths of otpl, the optimal template.

        ; TODO: Plotting option has been removed.  There has to be a better way to create QA plots.
        ;output_file_root=output_file_root, plot=plot, 
        ; TODO: Do some checks of the input arrays?

        version_module = '0.9'                  ; Version number
        if n_elements(version) ne 0 then begin  ; If version is defined
            version = version_module            ; ... set it to the module version
            return                              ; ... and return without doing anything
        endif

        ;-----------------------------------------------------------------------
        ; SETUP ----------------------------------------------------------------

        ; Set the inverse variance vectors to dummy values for use with otpl and
        ; botpl
        tpl_ivar = make_array((size(flux))[2], /double, value=1.0d)

        c=299792.458d                           ; Speed of light in km/s
        redshift = stellar_kinematics[*,0]/c    ; Redshift of each spectrum

        ; Determine the minimum and maximum wavelength available to the index
        ; measurements; force the index to be common to all spectra (flux, otpl,
        ; botpl)

;       min_wave = min(wave)            ; Grab the minimum
;       if min_wave lt min(tpl_wave) then $     ; Adjust it if necessary
;           min_wave = min(tpl_wave)
;       max_wave = max(wave)            ; Grab the maximum
;       if max_wave gt max(tpl_wave) then $     ; Adjust it if necessary
;           max_wave = max(tpl_wave)

        ; Allocate the output arrays
        nspec = (size(flux))[1]                         ; Number of spectra
        nabs = n_elements(abs_par)                      ; Number of spectral indices
        abs_line_indx_omitted = intarr(nspec, nabs)     ; Flag that index was omitted
        abs_line_indx = dblarr(nspec, nabs)             ; Index measurements (object spectra)
        abs_line_indx_err = dblarr(nspec, nabs)         ; ... error
        abs_line_indx_otpl = dblarr(nspec, nabs)        ; Index measurements (optimal template)
        abs_line_indx_botpl = dblarr(nspec, nabs)       ; Index measurements (broadened op template)

;       if keyword_set(plot) and n_elements(output_file_root) eq 0 then begin
;           print, 'No output root name provide for PS plots.  Using ./mdap_abs_line'
;           output_file_root = './mdap_abs_line'
;           set_plot,'ps'
;       endif
        ;-----------------------------------------------------------------------


        ;-----------------------------------------------------------------------
        ; Measure all spectral indices for all spectra -------------------------
        if keyword_set(dbg) then $              ; Only analyze the first spectrum
            nspec = 1
        for i=0,nspec-1 do begin

;           ; Initialize the plot file
;           if keyword_set(plot) then begin
;               device, filename=output_file_root+'indx_'+mdap_stc(i+1,/integer)+'.ps', /color, $
;                       xsize=15, ysize=12
;           endif

            ; Set the mask for otpl to be the DEREDSHIFTED version of the input
            ; mask
            otpl_mask = interpol(reform(mask[i,*]), wave, wave/(redshift[i]+1.0d))
            indx = where(otpl_mask gt 0.5, count, complement=nindx, ncomplement=ncount)
;           if indx[0] ne -1 then $
            if count ne 0 then $
                otpl_mask[indx] = 1.0d
            if ncount ne 0 then $
                otpl_mask[nindx] = 0.0d

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

;               ; Check if the index can be measured
;               if abs_par[j].red_cont[1] lt min_wave or $
;                  abs_par[j].blue_cont[0] gt max_wave or $
;                  abs_par_z[j].red_cont[1] lt min_wave or $
;                  abs_par_z[j].blue_cont[0] gt max_wave then begin
;                   abs_line_indx_omitted[i,j] = 1      ; Flag the index as omitted
;                   print, 'OMITTED!'
;                   continue                            ; Index not in wavelength range, go to next
;               endif

                ; Measure the index based on the object data (use redshifted windows)
                err = 0
                MDAP_GET_SPECTRAL_INDEX, wave, reform(flux[i,*]), reform(ivar[i,*]), $
                                         reform(mask[i,*]), abs_par_z[j], obj_equiv_width, $
                                         obj_equiv_width_err, obj_index_mag, obj_index_mag_err, $
                                         err=err, /geometric
                if err eq 1 then begin                          ; Measurement returned error!
                    abs_line_indx_omitted[i,j] = 1
                    continue
                endif
                ; Measure the index based on the optimal template (use rest-frame windows)
                MDAP_GET_SPECTRAL_INDEX, wave, reform(otpl[i,*]), tpl_ivar, otpl_mask, $
                                         abs_par[j], otpl_equiv_width, otpl_equiv_width_err, $
                                         otpl_index_mag, otpl_index_mag_err, err=err, /geometric
                if err eq 1 then begin                          ; Measurement returned error!
                    abs_line_indx_omitted[i,j] = 1
                    continue
                endif
                ; Measure the index based on the broadened template (use redshifted windows)
                MDAP_GET_SPECTRAL_INDEX, wave, reform(botpl[i,*]), tpl_ivar, reform(mask[i,*]), $
                                         abs_par_z[j], botpl_equiv_width, botpl_equiv_width_err, $
                                         botpl_index_mag, botpl_index_mag_err, err=err, /geometric
                if err eq 1 then begin                          ; Measurement returned error!
                    abs_line_indx_omitted[i,j] = 1
                    continue
                endif

                ; TODO: Instead of providing corrected values, provide direct
                ; measurements and the corrections?

                ; Correct the index measured using the difference between the
                ; measurements performed on the optimal template and its
                ; broadened version.  The correction depends on the desired
                ; units.  This also saves the results.
                if abs_par[j].units eq 'mag' then begin
                    losvd_correction = otpl_index_mag - botpl_index_mag
                    abs_line_indx[i,j] = obj_index_mag + losvd_correction
                    abs_line_indx_otpl[i,j] = otpl_index_mag
                    abs_line_indx_botpl[i,j] = botpl_index_mag
                    abs_line_indx_err[i,j] = obj_index_mag_err
                endif

                if abs_par[j].units eq 'ang' then begin
                    losvd_correction = otpl_equiv_width / botpl_equiv_width
                    abs_line_indx[i,j] = obj_equiv_width / losvd_correction
                    abs_line_indx_otpl[i,j] = otpl_equiv_width
                    abs_line_indx_botpl[i,j] = botpl_equiv_width
                    abs_line_indx_err[i,j] = obj_equiv_width_err / abs(losvd_correction)
                endif

;               print, abs_line_indx[i,j], abs_line_indx_err[i,j], abs_line_indx_otpl[i,j], $
;                      abs_line_indx_botpl[i,j]
;               stop

            endfor ; End loop over indices

;           if keyword_set(plot) then $
;               device, /close

        endfor ; End loop over spectra
        ;-----------------------------------------------------------------------

;       if keyword_set(plot) then $
;           set_plot, 'x'

END

