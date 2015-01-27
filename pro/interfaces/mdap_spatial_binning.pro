;+
; NAME:
;       MDAP_SPATIAL_BINNING
;
; PURPOSE:
;       Bin an input set of spectra using the binning parameters
;       provided.  See MDAP_DEFINE_BIN_PAR() for a description of the
;       available options.
;
;       No shifting of the spectrum is done!  The output wavelength
;       solution is the same as the input wavelength solution and
;       assumed to be the same for all spectra.  Shifts to the spectra
;       can be applied using MDAP_REDSHIFT_REGISTER_SPECTRA.
;
; CALLING SEQUENCE:
;       MDAP_SPATIAL_BINNING, flux, ivar, mask, signal, noise, gflag, xcoo, ycoo, dx, dy, $
;                             bin_par, threshold_ston_bin, bin_weights, binned_indx, binned_flux, $
;                             binned_ivar, binned_mask, binned_x_rl, binned_y_ru, binned_wrad, $
;                             binned_area, binned_ston, nbinned, sn_calibration=sn_calibration, $
;                             version=version, plot=plot
;
; INPUTS:
;       flux dblarr[N][T]
;               Galaxy spectra as produced by MDAP_READ_DRP_FITS.
;
;       ivar dblarr[N][T]
;               Inverse variance of the flux
;
;       mask dblarr[N][T]
;               Bad pixel mask.
;
;       signal dblarr[N]
;               Mean galaxy signal per angstrom
;
;       noise dblarr[N]
;               Mean galaxy error per angstrom
;
;       gflag intarr[N]
;               Flag (0=false; 1=true) that the spectrum is 'good' as defined by
;               MDAP_SELECT_GOOD_SPECTRA.  Spectra that are NOT good are ignored.
;
;       xcoo dblarr[N]
;               Array containing the x coordinates in arcseconds (0 is the
;               center of the field of view) for each spectrum
;
;       ycoo dblarr[N]
;               Array containing the y coordinates in arcseconds (0 is the
;               center of the field of view) for each spectrum
;
;       dx double
;               Scale arcsec/pixel in X direction
;
;       dy double
;               Scale arcsec/pixel in Y direction
;
;       bin_par BinPar[P]
;               Structure containing the parameters used to perform the
;               spatial binning.  See MDAP_DEFINE_BIN_PAR() for a
;               description of the structure.
;
;       threshold_ston_bin double
;               The S/N threshold for the inclusion of any DRP spectrum in the
;               binning process.
;               
; OPTIONAL INPUTS:
;
;       SN_CALIBRATION TODO: TYPE? or flag?
;               If provided, the estimated signal-to-noise (SN_est) is converted
;               into the real signal-to-noise using the empirical calibration
;               function defined in MDAP_CALIBRATE_SN:
;
;                       tmp = SN_EST^SN_CALIBRATION[0]/sqrt(n_elements(n_elements_within_bin)
;                       SN_REAL = poly(SN_EST,SN_CALIBRATION[1:*])
;
; OPTIONAL KEYWORDS:
;       \plot
;               If set, some plots on X11 terminal will be shown. Not suggested
;               if the task is launched remotely. 
;
; OUTPUT:
;       bin_weights dblarr[N]
;               Weights used for each spectrum for binning.             
;
;       binned_indx lonarr[N]
;               Indicates in which bin, i=0...B-1, each of the N spectra were
;               placed.
;
;       binned_flux dblarr[B][T]
;               The binned spectra of the spatial B bins. i-th spectrum is
;               associated to the i-th bin. 
;
;       binned_ivar dblarr[B][T]
;               Inverse variance of binned spectra.
;
;       binned_mask dblarr[B][T]
;               Pixel mask.
;
;       binned_x_rl dblarr[B]
;               If STON binning, x-coordinates in arcsec of the
;               luminosity-weighted centers of the spatial bins.
;
;               If RADIAL binning, lower edge of the radial bin.
; 
;       binned_y_ru dblarr[B]
;               If STON binning, y-coordinates in arcsec of the
;               luminosity-weighted centers of the spatial bins.
;
;               If RADIAL binning, upper edge of the radial bin.
; 
;       binned_wrad dblarr[B]
;               If STON binning, luminosity-weighted sky-plane radius of
;               the spectra in the bin.
;
;               IF RADIAL binning, luminosity-weighted in-plane radius
;               of the spectra in the radial bin.
;
;       binned_area dblarr[B]
;               Area (in arcsec^2) of each spatial bin.  This is just
;               the sum of the area of each spectrum in the bin.  This
;               does NOT properly account for the overlapping area for
;               the RSS spectra.
;
;       binned_ston dblarr[B]
;               Mean S/N per angstrom reached in each spatial bin. 
;
;       nbinned lonarr[B]
;               Number of spectra coadded in each bin.
;
; OPTIONAL OUTPUT:
;       version string
;               Module version. If requested, the module is not executed and only
;               the version flag is returned
;
; COMMENTS:
;
; EXAMPLES:
;
; TODO:
;       - Also include full polygon describing the bin in output? (I.e. keep
;         binned_xvec, binned_yvec)
;       - Include something that handles the covariance
;
;       Allow for user-defined binning
;       user_bin_map  string
;               If provided, the spatial map will be created from the fits file
;               specified by this input. The fits file must contain the CRVAL1,
;               CRVAL2, CDELT1, CDELT2, NAXIS1, NAXIS2, CRPIX1, and CRIX2 header
;               keywords (coordinate units should be in arcseconds; 0,0
;               indicates the center of the field of view).
;
; BUGS:
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;       01 Sep 2014: Copied from v0_8 by L. Coccato
;       08 Sep 2014: (KBW) Formatting and comments (incomplete)
;       15 Sep 2014: (KBW) Formatting and edits due to accommodate other changes
;       16 Sep 2014: (KBW) gflag changed from optional to required parameter
;       22 Sep 2014: (KBW) Output mask for combined spectra (TODO: This is just a place-holder)
;       13 Oct 2014: (KBW) Changed input/output format
;       04 Dec 2014: (KBW) Allow for radial binning; binned_indx changed
;                          to long; change to weight_by_sn2
;-
;------------------------------------------------------------------------------

PRO MDAP_SPATIAL_BINNING, $
                flux, ivar, mask, signal, noise, gflag, xcoo, ycoo, dx, dy, bin_par, $
                threshold_ston_bin, bin_weights, binned_indx, binned_flux, binned_ivar, $
                binned_mask, binned_x_rl, binned_y_ru, binned_wrad, binned_area, binned_ston, $
                nbinned, sn_calibration=sn_calibration, version=version, plot=plot

        version_module = '0.4'                          ; Version number

        if n_elements(version) ne 0 then begin          ; set version and return
            version = version_module
            return
        endif

        ; Find which spectra in the 2D map are good (and bad)
        ;       good = has a positive and finite noise, a finite signal,
        ;       and S/N > threshold
        gindx= where(gflag eq 1 and abs(signal/noise) ge threshold_ston_bin, count, compl=bindx)

;        if gindx[0] eq -1 then $                        ; No good spectra so return and fail
        if count eq 0 then $                        ; No good spectra so return and fail
            message, 'No good spectra!'

        print, ' unmasked pixels: ', n_elements(where(mask lt 1.))
        print, ' non-zero pixels: ', n_elements(where(flux gt 0.))

        ngood = n_elements(gindx)
        sz=size(flux)
        ns=sz[1]                                        ; Number of spectra

        if bin_par.type eq 'NONE' then begin                ; No binning
            bin_weights = dblarr(ns)                    ; Initialize weights to 0
            bin_weights[gindx] = 1.0d                   ; Set weights to unity
            binned_indx = make_array(ns, /long, value=-1)
            binned_indx[gindx] = indgen(ngood)
            binned_flux = flux[gindx,*]
            binned_ivar = ivar[gindx,*]
            binned_mask = mask[gindx,*]
            binned_x_rl = xcoo[gindx,*]
            binned_y_ru = ycoo[gindx,*]
            binned_wrad = sqrt( binned_x_rl^2 + binned_y_ru^2 )
            binned_area = make_array(ngood, /double, value=dx*dy)
            binned_ston = signal[gindx]/noise[gindx]
            nbinned = make_array(ngood, /long, value=1)
            return
        endif

        if bin_par.optimal_weighting eq 1 then $        ; Flag to use S/(N)^2 weighting
            optimal_weighting = 1

        ; Will only perform one of the following:
        if bin_par.type eq 'ALL' then begin

            MDAP_ALL_SPECTRA_BINNING, xcoo[gindx], ycoo[gindx], signal[gindx], noise[gindx], $
                                      binned_indx, binned_x_rl, binned_y_ru, binned_ston, $
                                      nbinned, sn_calibration=sn_calibration, $
                                      optimal_weighting=optimal_weighting

            ; Approximate the luminosity-weighted radius
            binned_wrad = sqrt( binned_x_rl^2 + binned_y_ru^2 )
        endif 
        
        if bin_par.type eq 'STON' then begin             ; Use the Voronoi binning scheme
            if keyword_set(plot) then begin                     ; setup plot
;               mydevice=!D.NAME
;               set_plot, 'PS'
;               device, filename='bin_plot.ps'
                screenr = GET_SCREEN_SIZE()
                window, xsize=screenr[0]*0.4, ysize=screenr[1]*0.8, retain=2
                MDAP_BASIC_COLORS, black, white, red, green, blue, yellow, cyan, magenta, orange, $
                                   mint, purple, pink, olive, lightblue, gray   
                
;               loadct, 32
            endif

            ; TODO: **** For now do NOT use weight_by_sn2.  Something
            ; seems inconsistent when using this option in the voronoi
            ; binning and the S/N calculation in
            ; MDAP_CALCULATE_BIN_SN().  This function does not seem to
            ; account for the change in the terms in signal and noise.
            ; Need to fix this!!!

            MDAP_VORONOI_2D_BINNING, xcoo[gindx], ycoo[gindx], signal[gindx], noise[gindx], $
                                     bin_par.ston, binned_indx, binned_xvec, binned_yvec, $
                                     binned_x_rl, binned_y_ru, binned_ston, v_area, v_scale, $
                                     sn_calibration=sn_calibration, $
                                     optimal_weighting=optimal_weighting, plot=plot, /quiet

            ; Approximate the luminosity-weighted radius
            binned_wrad = sqrt( binned_x_rl^2 + binned_y_ru^2 )

            nbins = n_elements(binned_x_rl)
            nbinned = lonarr(nbins)
            for i=0,nbins-1 do begin
                indx = where(binned_indx eq i, count)
;               if indx[0] eq -1 then $
                if count eq 0 then $
                    continue
                nbinned[i] = n_elements(indx)
            endfor

;           if keyword_set(plot) then begin                     ; close plot
;               device, /close
;               set_plot, mydevice
;           endif
        endif
        
        if bin_par.type eq 'RADIAL' then begin
            MDAP_RADIAL_BINNING, xcoo[gindx], ycoo[gindx], signal[gindx], noise[gindx], bin_par, $
                                 binned_indx, binned_x_rl, binned_y_ru, binned_wrad, binned_ston, $
                                 nbinned, sn_calibration=sn_calibration, $
                                 optimal_weighting=optimal_weighting 
        endif

        print, 'N per bin:'
        for i=0,n_elements(nbinned)-1 do $
            print, i+1, ':', nbinned[i]

        ; Remove any bins with zero spectra
        indx = where(nbinned eq 0, count, complement=nindx)
        if count ne 0 then begin
            print, 'Some bins had zero spectra!'

            ; Update the binned arrays
            nbinned = nbinned[nindx]
            binned_x_rl = binned_x_rl[nindx]
            binned_y_ru = binned_y_ru[nindx]
            binned_wrad = binned_wrad[nindx]
            binned_ston = binned_ston[nindx]

            ; Update the index designations for each of the DRP spectra
            for i=0,n_elements(nindx)-1 do begin
                bnindx = where(binned_indx eq nindx[i])
                binned_indx[bnindx] = i
            endfor
        endif

;       TODO: NOT DEFINED YET -----------------------------------
;       ; Check the user defined binning scheme, fault if check fails!
;       if bin_par.type eq 'USER'then begin
;           MDAP_READ_USER_DEFINED_SPATIAL_BINS, user_bin_map, header, binned_indx, success
;       endif
;       NOT DEFINED YET ----------------------------------------

        ; Set binned_indx to same length as flux
        MDAP_INSERT_FLAGGED, gindx, binned_indx, ns
        print, 'Number of spatial bins: ', MDAP_STC(n_elements(nbinned),/integer)

        ; Generate the weights to use in combining the spectra
        ; TODO: Should this be input/output from the voronoi routine in the 'STON' case?
        MDAP_GENERATE_BINNING_WEIGHTS, signal, noise, bin_weights, $
                                       optimal_weighting=optimal_weighting
        indx = where(binned_indx lt 0, count)
        if count ne 0 then $
            bin_weights[indx] = 0.0d

        ; Combine the spectra
        MDAP_COMBINE_SPECTRA, flux, ivar, mask, binned_indx, bin_weights, nbinned, binned_flux, $
                              binned_ivar, binned_mask

        ; Determine the effective on-sky area of each combined spectrum
        ; TODO: does not account for overlapping regions!!
        MDAP_SPATIAL_BIN_AREA, dx, dy, nbinned, binned_indx, binned_area

        ; TODO: How is binned_area used?  Should it include the weights?

END


