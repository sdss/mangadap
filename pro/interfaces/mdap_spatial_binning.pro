;+
; NAME:
;       MDAP_SPATIAL_BINNING
;
; PURPOSE:
;       Bin an input set of spectra to a minimum S/N level.
;
; CALLING SEQUENCE:
;       MDAP_SPATIAL_BINNING, flux, ivar, mask, signal, noise, gflag, xcoo, ycoo, dx, dy, $
;                             bin_type, bin_par, threshold_ston_bin, bin_weight_by_sn2, $
;                             bin_weights, binned_indx, binned_flux, binned_ivar, binned_mask, $
;                             binned_xcoo, binned_ycoo, binned_area, binned_ston, nbinned, $
;                             sn_calibration=sn_calibration, version=version, plot=plot
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
;       bin_type string
;               Type of binning to apply.
;               Valid values are:
;                       'NONE' - No binning is performed.
;                       'ALL' - All spectra are combined into a single bin.
;                       'STON' - Spectra are binned, using the Voronoi binning
;                                scheme, to a minimum S/N.  This type requires a
;                                single parameter (provided via bin_par), which
;                                is the minimum S/N level.
;                       'RAD' - Spectra are binned in radius.  This type
;                               requires five parameters, the x and y center,
;                               ellipticity, position angle, and width of the
;                               radial bins.  
;                       
;                               TODO: When binning, the spectra are deredshifted
;                               before creating the bin, meaning that a set of
;                               velocities for each spectrum must be available!
;                       
;                               TODO: Currently only the bin width is input!!
;                       
;       bin_par double
;               Parameter(s) required for the binning.
;               TODO: Change this to a pointer!
;
;       threshold_ston_bin double
;               The S/N threshold for the inclusion of any DRP spectrum in the
;               binning process.
;
;       bin_weight_by_sn2 integer
;               Flag to weight each spectrum by S/N^2 when producing the binned
;               spectra.  If true (eq 1), weights are S/N^2; if false (eq 0),
;               the spectra are added with uniform weights.
;               See MDAP_GENERATE_BINNING_WEIGHTS.
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
;       binned_indx intarr[N]
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
;       binned_xcoo dblarr[B]
;               X-Coordinates in arcsec of the luminosity-weighted centers of
;               the spatial bins. 
;
;       binned_ycoo dblarr[B]
;               Y-Coordinates in arcsec of the luminosity-weighted centers of
;               the spatial bins. 
;
;       binned_area dblarr[B]
;               Area (in arcsec^2) of each spatial bin.  
;
;       binned_ston dblarr[B]
;               Mean S/N per angstrom reached in each spatial bin. 
;
;       nbinned intarr[B]
;               Number of spectra coadded in each bin.
;
; OPTIONAL OUTPUT:
;
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
;-
;------------------------------------------------------------------------------

PRO MDAP_SPATIAL_BINNING, $
                flux, ivar, mask, signal, noise, gflag, xcoo, ycoo, dx, dy, bin_type, bin_par, $
                threshold_ston_bin, bin_weight_by_sn2, bin_weights, binned_indx, binned_flux, $
                binned_ivar, binned_mask, binned_xcoo, binned_ycoo, binned_area, binned_ston, $
                nbinned, sn_calibration=sn_calibration, version=version, plot=plot

        version_module = '0.3'                          ; Version number

        if n_elements(version) ne 0 then begin          ; set version and return
            version = version_module
            return
        endif

        ; Find which spectra in the 2D map are good (and bad)
        ;       good = has a positive and finite noise, a finite signal,
        ;       and S/N > threshold
        gindx= where(gflag eq 1 and abs(signal/noise) ge threshold_ston_bin, compl=bindx)

        if gindx[0] eq -1 then $                        ; No good spectra so return and fail
            message, 'No good spectra!'

        ngood = n_elements(gindx)
        sz=size(flux)
        ns=sz[1]                                        ; Number of spectra

        if bin_type eq 'NONE' then begin                ; No binning
            bin_weights = dblarr(ns)                    ; Initialize weights to 0
            bin_weights[gindx] = 1.0d                   ; Set weights to unity
            binned_indx = make_array(ns, /integer, value=-1)
            binned_indx[gindx] = indgen(ngood)
            binned_flux = flux[gindx,*]
            binned_ivar = ivar[gindx,*]
            binned_mask = mask[gindx,*]
            binned_xcoo = xcoo[gindx,*]
            binned_ycoo = ycoo[gindx,*]
            binned_area = make_array(ngood, /double, value=dx*dy)
            binned_ston = signal[gindx]/noise[gindx]
            nbinned = make_array(ngood, /integer, value=1)
            return
        endif

;       TODO: NOT DEFINED YET -----------------------------------
;       ; Check the user defined binning scheme
;       endif else if if n_elements(user_bin_map) ne 0 then begin
;           MDAP_READ_USER_DEFINED_SPATIAL_BINS, user_bin_map, header, binned_indx, success
;           if success eq 1 then apply_voronoi_binning = 0      ; Successful user definition
;       endelse
;       NOT DEFINED YET ----------------------------------------

        if bin_weight_by_sn2 eq 1 then $                        ; Flag to use S/N^2 weighting
            weight_for_sn = 1

        if bin_type eq 'ALL' then begin
            binned_indx = intarr(ngood)                         ; All spectra are in bin i=0
            nbinned = ngood
        endif else if bin_type eq 'STON' then begin             ; Use the Voronoi binning scheme
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

            MDAP_VORONOI_2D_BINNING, xcoo[gindx], ycoo[gindx], signal[gindx], noise[gindx], $
                                     bin_par, binned_indx, binned_xvec, binned_yvec, binned_xcoo, $
                                     binned_ycoo, binned_ston, nbinned, scale, $
                                     sn_calibration=sn_calibration, weight_for_sn=weight_for_sn, $
                                     plot=plot, /quiet

;           if keyword_set(plot) then begin                     ; close plot
;               device, /close
;               set_plot, mydevice
;           endif
        endif else if bin_type eq 'RAD' then $
            message, 'Cannot use radial binning scheme yet!'

        ; Set binned_indx to same length as flux
        MDAP_INSERT_FLAGGED, gindx, binned_indx, ns

        print, 'Number of spatial bins: ', MDAP_STC(n_elements(nbinned),/integer)

        ; Generate the weights to use in combining the spectra
        ; TODO: Should this be input/output from the voronoi routine in the 'STON' case?
        MDAP_GENERATE_BINNING_WEIGHTS, signal, noise, bin_weights, weight_for_sn=weight_for_sn
        bin_weights[where(binned_indx lt 0)] = 0.0d

        ; Combine the spectra
        MDAP_COMBINE_SPECTRA, flux, ivar, mask, binned_indx, bin_weights, nbinned, binned_flux, $
                              binned_ivar, binned_mask

        ; Determine the effective on-sky area of each combined spectrum
        ; TODO: does not account for overlapping regions!!
        MDAP_SPATIAL_BIN_AREA, dx, dy, nbinned, binned_indx, binned_area

        ; TODO: How is binned_area used?  Should it include the weights?

END


