;+
; NAME:
;       MDAP_RADIAL_BINNING
;
; PURPOSE:
;       Determine the radial bins for a set of spectra with the provided
;       set of on-sky x and y coordinates.  The coordinates should be in
;       arcseconds from some reference position.  This reference
;       position is at (bin_par.cx, bin_par.cy) with respect to the
;       galaxy center.  The x coordinates must increase toward the east
;       and decrease toward the west.  The position angle must be
;       defined as the angle from North through east.
;
; CALLING SEQUENCE:
;       MDAP_RADIAL_BINNING, fibx, fiby, signal, noise, bin_par, binned_indx, binned_rlow, $
;                            binned_rupp, binned_wrad, binned_ston, nbinned
;
; INPUTS:
;       fibx dblarr[N]
;               X coordinate of the fiber in arcsec in the reference
;               frame of the IFU for each of the N spectra
;
;       fiby dblarr[N]
;               Y coordinate of the fiber in arcsec in the reference
;               frame of the IFU for each of the N spectra
;
;       signal dblarr[N]
;               Mean galaxy signal per angstrom
;
;       noise dblarr[N]
;               Mean galaxy error per angstrom
;
;       bin_par BinPar structure
;               Structure containing the needed parameters for the
;               radial binning.  See MDAP_DEFINE_BIN_PAR().
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;
; OUTPUT:
;       binned_indx intarr[N]
;               Indicates in which bin, i=0...B-1, each of the N spectra were
;               placed.
;
;       binned_rlow dblarr[B]
;               Lower edge of the B radial bins.
;
;       binned_rupp dblarr[B]
;               Upper edge of the B radial bins.
;
;       binned_wrad dblarr[B]
;               Signal-weighted radius of the spectra in each radial bin B.
;
;       binned_ston dblarr[B]
;               S/N of the combined spectra in the bin.
;
;       nbinned intarr[B]
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
;       - Add versioning?
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;       03 Dec 2014: Adapted from some C++ code by K. Westfall (KBW)
;       16 Mar 2015: (KBW) Change to calibrated S/N calculation
;       14 Aug 2015: (KBW) Edited the documentation to clarify the x
;                          coordinate convention (+x is toward east, -x
;                          toward west).
;-
;-----------------------------------------------------------------------

;-----------------------------------------------------------------------
; Calculate the radius assuming a planar geometry, projected into the
; plane of the sky.  A circle in this plane will project into an ellipse
; with a certain position angle (angle from North, at 0 deg, through
; East, at 90 deg) and ellipticity (1-b/a).  The provided coordinates
; (x,y) are in the focal plane, which is pointed at (xc,yc) with on-sky
; North at an angle on_sky_rotation relative to the y-axis of the focal
; plane.
;
; The solution of the circle cartesian coordinates are computed using
; the following set of linear equations:
;
;     |   1.0  0.0   0.0  0.0  0.0  0.0 |       |  x |        |   x |
;     |   0.0  1.0   0.0  0.0  0.0  0.0 |       |  y |        |   y |
; A = |  coso sino  -1.0  0.0  0.0  0.0 | * x = | xs | =  b = | -xc |
;     | -sino coso   0.0 -1.0  0.0  0.0 |       | ys |        | -yc |
;     |   0.0  0.0  sinp cosp -1.0  0.0 |       | xd |        | 0.0 |
;     |   0.0  0.0 -cosp sinp  0.0 -b/a |       | yd |        | 0.0 |
;
;
; where (xs,ys) are the sky-plane coordinates of the fiber at (x,y), and
; (xd,yd) are the in-plane coordinates.  The solution to the linear set
; of equations is determined by first determining the LU decomposition
; of A, then solving the system (and improving it) using the
; fiber-by-fiber coordinate solution (b).
; 
; Once x is solved, the in-plane radius and azimuth are then:
;
;       Rd^2 = xd^2 + yd^2
;       thetad = tan^-1(-yd/xd)
;
; The latter are returned.  All angles are input and returned in
; degrees.

PRO MDAP_MAJOR_AXIS_POLAR_COO, $
                x, y, xc, yc, on_sky_rotation, position_angle, ellipticity, R, theta

        pi = 3.141592653589793d                                     ; Define pi

        boa = 1.0d - ellipticity                                    ; For convenience
        sinor = sin(on_sky_rotation*pi/180.0d)
        cosor = cos(on_sky_rotation*pi/180.0d)
        sinpa = sin(position_angle*pi/180.0d)
        cospa = cos(position_angle*pi/180.0d)

        A = [ [   1.0d,  0.0d,   0.0d,  0.0d,  0.0d, 0.0d ], $      ; Define A
              [   0.0d,  1.0d,   0.0d,  0.0d,  0.0d, 0.0d ], $
              [  cosor, sinor,  -1.0d,  0.0d,  0.0d, 0.0d ], $
              [ -sinor, cosor,   0.0d, -1.0d,  0.0d, 0.0d ], $
              [   0.0d,  0.0d,  sinpa, cospa, -1.0d, 0.0d ], $
              [   0.0d,  0.0d, -cospa, sinpa,  0.0d, -boa ] ]

        ALU = A
        LA_LUDC, ALU, indx                          ; Get its LU decomposition

        b =   [   x[0],  y[0],    -xc,   -yc,  0.0d, 0.0d ]         ; Define b

        nx = n_elements(x)                          ; Allocate space for output data
        R = dblarr(nx,/nozero)
        theta = dblarr(nx,/nozero)
        for i=0,nx-1 do begin

            b[0] = x[i]                             ; Modify b for this spectrum
            b[1] = y[i]

            x0 = LA_LUSOL(ALU, indx, b)             ; Solve
            xi = LA_LUMPROVE(A, ALU, indx, b, x0)   ; ... and improve the solution to x

            R[i] = sqrt(xi[4]*xi[4] + xi[5]*xi[5])  ; Get radius
            theta[i] = atan(-xi[5],xi[4])*180.0d/pi ; ... and azimuth (in degrees)
        endfor
END

PRO MDAP_RADIAL_BINNING, $
                fibx, fiby, signal, noise, bin_par, binned_indx, binned_rlow, binned_rupp, $
                binned_wrad, binned_ston, nbinned

        ; Compute the "major-axis" coordinates of each fiber
        MDAP_MAJOR_AXIS_POLAR_COO, fibx, fiby, bin_par.cx, bin_par.cy, 0.0d, bin_par.pa, $
                               bin_par.ell, R, theta

        ; Scale R if requested
;        print, 'RSCALE: ', bin_par.rscale
        if bin_par.rscale gt 0.0 then $
            R = R / bin_par.rscale

;        print, bin_par.cx, bin_par.cy, 0.0d, bin_par.pa, bin_par.ell
;
;        for i=0,9 do $
;            print, fibx[i], fiby[i], R[i], theta[i]
;
;        stop

;        print, 'input num: ', bin_par.nr
;        print, 'input end: ', bin_par.re

        ; Get the maximum radius to use in the binning, if not provided
        if bin_par.re lt 0.0 then $
            bin_par.re = max(R)+0.1d

;        print, 'max, re: ', max(R), bin_par.re

;       print, 'RS: ', bin_par.rs
;       print, 'minR', R[r_sort_indx[0]], R[r_sort_indx[1]], R[r_sort_indx[2]], R[r_sort_indx[3]] 
;       print, 'RE: ', bin_par.re
;       print, 'maxR', max(R)

;        print, 'input start: ', bin_par.rs

        ; Get the lower and upper edges of the radial bins
        if bin_par.rlog eq 1 then begin
            ; Minimum r must be positive for logarithmic binning
            minr = 0.1/bin_par.rscale           ; 0.1 arcsec in units of the scale radius
            if bin_par.rs lt minr then begin
                indx = where(R gt minr)
;                print, 'min > 0.1/rscale: ', min(R[indx])
                bin_par.rs = min(R[indx])/1.1
            endif
            bin_edges = MDAP_RANGE(bin_par.rs, bin_par.re, bin_par.nr+1, /log)
        endif else $
            bin_edges = MDAP_RANGE(bin_par.rs, bin_par.re, bin_par.nr+1)

;        print, 'min, rs: ', min(R), bin_par.rs

        binned_rlow = bin_edges[0:bin_par.nr-1]
        binned_rupp = bin_edges[1:bin_par.nr]

;        print, bin_par.rs, bin_par.re, bin_par.rlog
        print, 'Radial bin limits:'
        for i=0,bin_par.nr-1 do $
            print, binned_rlow[i], binned_rupp[i]

;       stop

        binned_wrad = dblarr(bin_par.nr,/nozero)        ; Allocate array with the weighted radius
        binned_ston = dblarr(bin_par.nr)                ; ... the S/N per bin
        nbinned = lonarr(bin_par.nr)                    ; ... and the number in bin

        n = n_elements(fibx)                            ; Number of spectra
        binned_indx = make_array(n, /long, value=-1)    ; Allocate binned index array

        if bin_par.optimal_weighting eq 1 then $        ; Flag to use S/(N)^2 weighting
            optimal_weighting = 1

        ; Get the binned indices
        print, 'Signal-weighted radius, S/N, number in bin:'
        for i=0,bin_par.nr-1 do begin
            indx = where(R gt binned_rlow[i] and R lt binned_rupp[i], count)   ; Spectra in the bin
;           if indx[0] eq -1 then begin                 ; No spectra in the bin
            if count eq 0 then begin                 ; No spectra in the bin
                binned_wrad[i] = -1.0d                  ; Set a dummy value
                continue                                ; nbinned and ston already 0
            endif

            binned_indx[indx] = i                       ; Set the bin the spectra are in

            wsum = total(signal[indx])                  ; Calculate the weighted radius
            binned_wrad[i] = total(signal[indx]*R[indx])/wsum

            ; Calculate the S/N
            binned_ston[i] = MDAP_CALCULATE_BIN_SN(signal[indx], noise[indx], $
                                                   noise_calib=bin_par.noise_calib, $
                                                   optimal_weighting=optimal_weighting)

            nbinned[i] = n_elements(indx)               ; Add the number of spectra in the bin

            print, binned_wrad[i], binned_ston[i], nbinned[i]

        endfor
;       stop
END



