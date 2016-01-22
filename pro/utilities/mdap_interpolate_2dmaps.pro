;+
; NAME:
;       MDAP_INTERPOLATE_2DMAPS
;
; PURPOSE:
;       This routine:
;
;       - uses the GRID_TPS IDL function to interpolate a set of N input values,
;         defined over an irregular grid, on to a regular grid based on the
;         input coordinates and sampling values.
;
;       - interpolates the values over this regular grid and re-interpolates
;         them over an irregular grid using the IDL BILINEAR function.
;
; CALLING SEQUENCE:
;       MDAP_INTERPOLATE_2DMAPS, x_in, y_in, z_in, dx, dy, x_out, y_out, z_out
;
; INPUTS:
;       x_in dblarr[N]
;               X-coordinates for input values.
;
;       y_in dblarr[N]
;               Y-coordinates for input values.
;
;       z_in dblarr[N]
;               'Z' values to interpolate.
;
;       dx double
;               Grid size in X for interpolated grid.
;
;       dy double
;               Grid size in Y for interpolated grid.
;
;       x_out dblarr[M]
;               X-coordinates for output values.
;
;       y_out dblarr[M]
;               Y-coordinates for output values.
;
; OPTIONAL INPUTS:
;       zlim dblarr[2]
;               Low, zlim[0], and high, zlim[1], limits for the output
;               interpolated values.  
;
;       default double
;               Default value in the case the interpolation cannot be
;               done or is outside the provided limits.  If not
;               specified, the default value is 0.
;
; OPTIONAL KEYWORDS:
;       /quiet
;               Don't print any messages.
;
; OUTPUT:
;       z_out dblarr[M]
;               Interpolated 'Z' values.
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
;       23 Sep 2014: Copied from v0_8 by L. Coccato
;       24 Sep 2014: (KBW) Formatting and some edits; full_grid option
;                          removed (may bring it back later);
;                          intermediate regular grid is set to encompass
;                          all of x_in, y_in, x_out, and y_out and
;                          sampled using the input pixel size, done to
;                          allow for better integration of RSS spectra.
;       03 Dec 2014: (KBW) Allow to impose limits on the output 'Z'
;                          values
;-
;------------------------------------------------------------------------------

; Get the intermediate grid to use in one coordinate
PRO MDAP_INTERPOLATE_2DMAPS_GRID, $
                x, dx, buffer, xs, nx
        xs = min(x)                         ; Minimum x coordinate
        xe = max(x)                         ; Maximum x coordinate
        Dx = buffer*(xe-xs)                 ; Full span in x including buffer
        nx = ceil(Dx/dx)                    ; Number x pixels for the intermediate grid
        xs = (xe+xs-Dx)/2.                  ; Offset initial value to account for buffer
END

PRO MDAP_INTERPOLATE_2DMAPS, $
                x_in, y_in, z_in, dx, dy, x_out, y_out, z_out, zlim=zlim, default=default, $
                quiet=quiet

        ; TODO: Check the dimensionality of the input arrays

        ; Initialize the output
        nout = n_elements(x_out)                ; Number of output values
        z_out = dblarr(nout, /nozero)           ; Allocate the output array, but do not initialize

        ; Do not consider input data that is not finite
        indx = where( finite(z_in) eq 1, count )
;       if indx[0] eq -1 or n_elements(indx) lt 7 then $        ; Limit of 7 is for use of GRID_TPS
        if count lt 7 then $        ; Limit of 7 is for use of GRID_TPS
            if ~keyword_set(quiet) then $
                print, 'ERROR: Insufficient number of valid measurements!'
            if n_elements(default) ne 0 then begin
                z_out[*] = default
            endif else $
                z_out[*] = 0.0
            return
        endif
                
        ; Get the size of the intermediate grid
        buffer = 1.1                            ; Use a 10% buffer
                                                ; TODO: Make this an optional input?

        MDAP_INTERPOLATE_2DMAPS_GRID, [x_in,x_out], dx, buffer, xs, nx  ; Get the x grid
        MDAP_INTERPOLATE_2DMAPS_GRID, [y_in,y_out], dy, buffer, ys, ny  ; Get the y grid

        ; Interpolate the input values to the intermediate grid
        z_in_grid = GRID_TPS(x_in[indx], y_in[indx], z_in[indx], ngrid=[nx,ny], start=[xs,ys], $
                             delta=[dx,dy])

        ; Interpolate the output values based on the regular grid
        for i=0,nout-1 do $
            z_out[i] = bilinear(z_in_grid, (x_out[i]-xs)/dx, (y_out[i]-ys)/dy)

END

