;+
; NAME:
;       MDAP_BINNING_BLOCK
;
; PURPOSE:
;   MAIN BLOCK FOR MANGA_DAP:
;
;   If the block is to be performed:
;       - Velocity register the DRP spectra, if requested
;       - Stack spectra according to the predefined binning parameters
;       - Write the binned data to the output file
;   Else:
;       - Read the binned data from the existing output file
;
; CALLING SEQUENCE:
;       MDAP_BINNING_BLOCK, manga_dap_version, execution_plan, perform_block, header, spaxel_dx, $
;                           spaxel_dy, wave, sres, flux, ivar, mask, bskyx, bskyy, gflag, signal, $
;                           noise, velocity_initial_guess, star_kin_interp, gas_kin_interp, $
;                           bin_vreg, reg_velocity_initial_guess, bin_wgts, bin_indx, bin_flux, $
;                           bin_ivar, bin_mask, xbin, ybin, bin_rad, bin_area, bin_ston, nbin, $
;                           /plot, /nolog, log_file_unit=log_file_unit, /quiet
;
; INPUTS:
;       manga_dap_version MaNGADAPVersion
;               Structure used to keep track of various
;               version-controlled procedures in the DAP.
;      
;       execution_plan ExecutionPlan
;               Structure providing plan and parameters needed for the
;               analyses.
;
;       perform_block RequiredAnalysisBlock
;               Structure defining which analysis blocks are to be
;               performed.
;
;       header strarr[]
;               Fits header
;            
;       spaxel_dx double
;               Spaxel size in on-sky X (arcsec).
;
;       spaxel_dy double
;               Spaxel size in on-sky Y (arcsec).
;
;       wave dblarr[T]
;               Wavelength of each spectral channel T.
;       
;       sres dblarr[T]
;               Spectral resolution at each wavelength channel T.
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
;       bskyx dblarr[N]
;               Array containing the x coordinates in arcseconds (0 is the
;               center of the field of view) for each spectrum
;
;       bskyy dblarr[N]
;               Array containing the y coordinates in arcseconds (0 is the
;               center of the field of view) for each spectrum
;
;       gflag intarr[N]
;               Flag (0=false; 1=true) that the spectrum is 'good' as defined by
;               MDAP_SELECT_GOOD_SPECTRA.  Spectra that are NOT good are ignored.
;
;       signal dblarr[N]
;               Mean galaxy signal per angstrom
;
;       noise dblarr[N]
;               Mean galaxy error per angstrom
;
;       velocity_initial_guess double
;               Initial guess velocity for the galaxy.
;
;       star_kin_interp dblarr[N][4]
;               Interpolated set of stellar kinematics based on a
;               provided analysis prior.
;
;       gas_kin_interp dblarr[N][2]
;               Interpolated set of gas kinematics based on a provided
;               analysis prior.
;
; OPTIONAL INPUTS:
;       log_file_unit LUN
;               File unit pointing to the log file
;
; OPTIONAL KEYWORDS:
;       /to_surface_brightness
;               Use the aperture area (spaxel_dx*spaxel_dy) to convert
;               the binned spectra to the average surface brightness in
;               the bin in on-sky units (e.g., arcsec^2).
;
;       /plot
;               Show the plot produced by the spatial binning procedure.
;
;       /nolog
;               Suppress output to the log file.
;
;       /quiet
;               Suppress output to the screen.
;
; OUTPUT:
;       bin_vreg dblarr[N]
;               Velocity used to register each of the DRP spectra.
;
;       reg_velocity_initial_guess double
;               Initial velocity guess that takes into account the change
;               made by registering the velocities of the DRP spectra.
;
;       bin_wgts dblarr[N]
;               Weights used when stacking each DRP spectrum.
;       
;       bin_indx intarr[N]
;               Index of the bin that contains each DRP spectrum.
;
;       bin_flux dblarr[B][T]
;               Flux in each of the B binned spectra at each of the T
;               spectral channels.
;
;       bin_ivar dblarr[B][T]
;               Inverse variance in the binned spectra.
;
;       bin_mask dblarr[B][T]
;               Bad pixel mask for the binned spectra.
;
;       xbin dblarr[B]
;               If STON binning, x-coordinates in arcsec of the
;               luminosity-weighted centers of the spatial bins.
;
;               If RADIAL binning, lower edge of the radial bin.
; 
;       ybin dblarr[B]
;               If STON binning, y-coordinates in arcsec of the
;               luminosity-weighted centers of the spatial bins.
;
;               If RADIAL binning, upper edge of the radial bin.
; 
;       bin_rad dblarr[B]
;               If STON binning, luminosity-weighted sky-plane radius of
;               the spectra in the bin.
;
;               IF RADIAL binning, luminosity-weighted in-plane radius
;               of the spectra in the radial bin.
;
;       bin_area dblarr[B]
;               Area (in arcsec^2) of each spatial bin.  This is just
;               the sum of the area of each spectrum in the bin.  This
;               does NOT properly account for the overlapping area for
;               the RSS spectra.
;
;       bin_ston dblarr[B]
;               Mean S/N per angstrom reached in each spatial bin. 
;
;       nbin lonarr[B]
;               Number of spectra coadded in each bin.
;
; OPTIONAL OUTPUT:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; TODO:
;
; PROCEDURES CALLED:
;       MDAP_SPATIAL_BINNING
;       MDAP_VELOCITY_REGISTER
;       SXADDPAR
;       MDAP_WRITE_OUTPUT
;       MDAP_DEFINE_OUTPUT
;       MDAP_READ_OUTPUT
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;       01 Feb 2015: Pulled from manga_dap.pro by K. Westfall (KBW)
;       16 Mar 2015: (KBW) Include adjustments for bin_par.noise_calib
;       17 Jun 2015: (KBW) Fixed a bug in logging; may be rampant!
;       14 Aug 2015: (KBW) Allowed binned flux to be converted from mean
;                          flux per aperture (e.g., fiber) to per unit
;                          area (e.g., arcsec^2) using the provided
;                          spaxel dx and dy.  This is assumed to be
;                          appropriate for both the CUBE and RSS
;                          spectra.
;       22 Sep 2015: (KBW) Fixed bug in propagation of error for flux
;                          per unit area (calculated error, not inverse
;                          variance!)
;-
;------------------------------------------------------------------------------

PRO MDAP_BINNING_BLOCK, $
                manga_dap_version, execution_plan, perform_block, header, spaxel_dx, spaxel_dy, $
                wave, sres, flux, ivar, mask, bskyx, bskyy, gflag, signal, noise, $
                velocity_initial_guess, star_kin_interp, gas_kin_interp, bin_vreg, $
                reg_velocity_initial_guess, bin_wgts, bin_indx, bin_flux, bin_ivar, bin_mask, $
                xbin, ybin, bin_rad, bin_area, bin_ston, nbin, $
                to_surface_brightness=to_surface_brightness, plot=plot, nolog=nolog, $
                log_file_unit=log_file_unit, quiet=quiet

        if perform_block.bin eq 1 then begin

            ;-----------------------------------------------------------
            ; Version control
            MDAP_SPATIAL_BINNING, version=manga_dap_version.spatial_binning

            ;-----------------------------------------------------------
            ; Do not overwrite original arrays
            reg_flux = flux
            reg_ivar = ivar
            reg_mask = mask

            reg_velocity_initial_guess = velocity_initial_guess

            ;-----------------------------------------------------------
            ; De-redshift the spectra before binning the spectra.
;           if strlen(execution_plan.analysis_prior) ne 0 $
;              and execution_plan.bin_par.v_register eq 1 $
;              and n_elements(star_kin_interp) ne 0 then begin
            if strlen(execution_plan.analysis_prior) ne 0 $
               && execution_plan.bin_par.v_register eq 1 $
               && n_elements(star_kin_interp) ne 0 then begin

                ; Velocity register the spectra to the median of the
                ; input stellar velocity.  Code is HARD-WIRED to
                ; register to the median.  Only de-redshift the "good"
                ; spectra defined by gflag (0-bad; 1-good; see
                ; MDAP_SELECT_GOOD_SPECTRA).
                bin_vreg = star_kin_interp[*,0]
                MDAP_VELOCITY_REGISTER, wave, reg_flux, reg_ivar, reg_mask, bin_vreg, gflag=gflag, $
                                        /register_to_median

                ; Offset interpolations for later use; these are okay to
                ; overwrite because they are re-interpolated (or erased)
                ; for each plan
                star_kin_interp[*,0] = star_kin_interp[*,0] - bin_vreg
                if n_elements(gas_kin_interp) ne 0 then $
                    gas_kin_interp[*,0] = gas_kin_interp[*,0] - bin_vreg

                ; Reset the velocity guess to the median of star_kin_interp
                reg_velocity_initial_guess = median(star_kin_interp[*,0])

;               for j=0,9 do $
;                   print, bin_vreg[i]
            endif else $
                bin_vreg = dblarr((size(flux))[1])             ; Set to 0.0

;           print, 'unmasked REG pixels:', n_elements(where(reg_mask lt 1.))
;           print, 'non-zero REG pixels:', n_elements(where(reg_flux gt 0.))
;           ns = (size(reg_flux))[1]
;           for i=0,ns-1 do begin
;               if gflag[i] eq 1 then begin
;                   plot, wave, reg_flux[i,*]
;                   stop
;               endif
;           endfor
;           stop

            MDAP_SPATIAL_BINNING, reg_flux, reg_ivar, reg_mask, signal, noise, gflag, bskyx, $
                                  bskyy, spaxel_dx, spaxel_dy, execution_plan.bin_par, $
                                  execution_plan.threshold_ston_bin, bin_wgts, bin_indx, bin_flux, $
                                  bin_ivar, bin_mask, xbin, ybin, bin_rad, bin_area, bin_ston, $
                                  nbin, plot=plot, quiet=quiet

            if keyword_set(to_surface_brightness) then begin
                aperture_area = spaxel_dx*spaxel_dy
                bin_flux = bin_flux / aperture_area
                indx = where(bin_ivar gt 0.0)
                bin_ivar[indx] = bin_ivar[indx] * aperture_area * aperture_area

                ; Change the units in the header
                SXADDPAR, header, 'BUNIT', '1E-17 erg/s/cm^2/Ang/arcsec^2', $
                      'Flux units are per arcsec^2'
                SXADDPAR, header, 'APAREA', aperture_area, 'Area of each spectrum area in arcsec^2'
            endif

            ; Write the version of the spatial binning to the header
            SXADDPAR, header, 'VDAPBIN', manga_dap_version.spatial_binning, $
                      'mdap_spatial_binning version'

            ; Add some information to the log
            if ~keyword_set(nolog) then begin
                printf, log_file_unit, '[INFO] Spatial bining version: ' + $
                        manga_dap_version.spatial_binning
                printf, log_file_unit, '[INFO] Type of spatial binning: ', $
                        execution_plan.bin_par.type
                if keyword_set(to_surface_brightness) then begin
                    printf, log_file_unit, '[INFO] Flux converted to on-sky area.'
                    printf, log_file_unit, '[INFO]    Aperture area: ', spaxel_dx*spaxel_dy
                endif

                if execution_plan.bin_par.type ne 'NONE' then begin
                    printf, log_file_unit, '    Velocity register the spectra before binning: ', $
                            execution_plan.bin_par.v_register
                    printf, log_file_unit, '    Use S/(N)^2 weighting when combining spectra: ', $
                            execution_plan.bin_par.optimal_weighting
                    printf, log_file_unit, '    Change nominal noise based on S/N calibration: ', $
                            execution_plan.bin_par.noise_calib
                endif
                if execution_plan.bin_par.type eq 'STON' then $
                    printf, log_file_unit, '    Minimum S/N per bin: ', execution_plan.bin_par.ston
                if execution_plan.bin_par.type eq 'RADIAL' then begin
                    printf, log_file_unit, '    X center: ', execution_plan.bin_par.cx
                    printf, log_file_unit, '    Y center: ', execution_plan.bin_par.cy
                    printf, log_file_unit, '    Position angle: ', execution_plan.bin_par.pa
                    printf, log_file_unit, '    Ellipticity: ', execution_plan.bin_par.ell
                    printf, log_file_unit, '    Starting radius: ', execution_plan.bin_par.rs
                    printf, log_file_unit, '    Ending radius (-1 for maximum): ', $
                            execution_plan.bin_par.re
                    printf, log_file_unit, '    Number of radial bins: ', execution_plan.bin_par.nr
                    printf, log_file_unit, '    Logarithmic bins: ', execution_plan.bin_par.rlog
                    printf, log_file_unit, '    Radius scale parameter: ', $
                            execution_plan.bin_par.rscale
                endif

                printf, log_file_unit, '[INFO] Number of bins: ', n_elements(nbin)
            endif

            ; Write the binning results
            MDAP_WRITE_OUTPUT, execution_plan.ofile, header=header, $
                               bin_par=execution_plan.bin_par, $
                               threshold_ston_bin=execution_plan.threshold_ston_bin, $
                               bin_vreg=bin_vreg, bin_indx=bin_indx, bin_weights=bin_wgts, $
                               wave=wave, sres=sres, bin_flux=bin_flux, bin_ivar=bin_ivar, $
                               bin_mask=bin_mask, xbin_rlow=xbin, ybin_rupp=ybin, rbin=bin_rad, $
                               bin_area=bin_area, bin_ston=bin_ston, bin_n=nbin, $; /read_header, $
                               quiet=quiet
        endif else begin

            ; Always need the binned data
            print, 'READING EXISTING BIN DATA'

            ; Read the binning results
            MDAP_DEFINE_OUTPUT, header=header, bin_vreg=bin_vreg, bin_indx=bin_indx, $
                                bin_weights=bin_wgts, bin_flux=bin_flux, bin_ivar=bin_ivar, $
                                bin_mask=bin_mask, xbin_rlow=xbin, ybin_rupp=ybin, rbin=bin_rad, $
                                bin_area=bin_area, bin_ston=bin_ston, bin_n=nbin

            MDAP_READ_OUTPUT, execution_plan.ofile, header=header, bin_vreg=bin_vreg, $
                              bin_indx=bin_indx, bin_weights=bin_wgts, bin_flux=bin_flux, $
                              bin_ivar=bin_ivar, bin_mask=bin_mask, xbin_rlow=xbin, $
                              ybin_rupp=ybin, rbin=bin_rad, bin_area=bin_area, bin_ston=bin_ston, $
                              bin_n=nbin

            ; As a sanity check, compare the size of the read flux against the DRP data
            sz = size(bin_flux)
            if sz[2] ne n_elements(wave) then $
                message, 'Read binned spectra do not have the same length as the DRP spectra!'

        endelse

END


