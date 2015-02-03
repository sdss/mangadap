;+
; NAME:
;
; PURPOSE:
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;
; OUTPUT:
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
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;
;-
;------------------------------------------------------------------------------

PRO MDAP_VELOCITY_FIELD, $
                ifile, ofile, map=map

        MDAP_DEFINE_OUTPUT, xpos=drpx, dx=dx, ypos=drpy, dy=dy, bin_indx=drp_bin, xbin=x_s, $
                            ybin=y_s, stellar_kinematics_fit=kin_s
        MDAP_READ_OUTPUT, ifile, xpos=drpx, dx=dx, ypos=drpy, dy=dy, bin_indx=drp_bin, xbin=x_s, $
                          ybin=y_s, stellar_kinematics_fit=kin_s

        ; Force the grid to be square
        if abs(dy-dx) gt 1e-6  then begin
            if dy gt dx then begin
                dy = dx
            endif else $
                dx = dy
        endif

        MDAP_COORDINATE_GRID, x_s, dx, xs, nx, aux=drpx, fbuffer=1.1
        print, xs, dx, nx

        MDAP_COORDINATE_GRID, y_s, dy, ys, ny, aux=drpy, fbuffer=1.1
        print, ys, dy, ny

        redshift=kin_s[*,0]


        for i=0,9 do $
            print, x_s[i], y_s[i], redshift[i]

        MDAP_PROPERTY_MAP_RENDER, x_s, y_s, redshift, xs, dx, nx, ys, dy, ny, map, $
                                  default=mean(redshift), /tps

;        MDAP_PROPERTY_MAP_RENDER, x_s, y_s, redshift, xs, dx, nx, ys, dy, ny, map, $
;                                  default=0.0

        MDAP_PROPERTY_MAP_INTERPOLATE, map, xs, dx, ys, dy, x_s, y_s, interp

        for i=0,9 do $
            print, x_s[i], y_s[i], redshift[i], interp[i]

        print, mean(redshift-interp)

        writefits, ofile, map

END



