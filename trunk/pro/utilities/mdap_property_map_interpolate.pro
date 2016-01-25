;+
; NAME:
;       MDAP_PROPERTY_MAP_INTERPOLATE
;
; PURPOSE:
;       Use bilinear interpolation to sample an input property map at
;       the provided coordinates.
;
; CALLING SEQUENCE:
;       MDAP_PROPERTY_MAP_INTERPOLATE, z_map, map_xs, map_dx, map_ys, map_dy, x_out, y_out, $
;                                      z_out, llim=llim, ulim=ulim, default=default
;
; INPUTS:
;       z_map dblarr[NX][NY]
;               Mapped data in a grid of size NX by NY.
;
;       map_xs double
;               X-coordinate of the first pixel in z_map.
;
;       map_dx double
;               Pixel width in X in z_map.
;
;       map_ys double
;               Y-coordinate of the first pixel in z_map.
;
;       map_dy double
;               Pixel width in Y in z_map.
;
;       x_out dblarr[N]
;               Array with X-coordinate at which to sample z_map.
;
;       y_out dblarr[N]
;               Array with Y-coordinate at which to sample z_map.
;
;
; OPTIONAL INPUTS:
;       llim double
;               Lower limit for the sampled quantities.
;
;       ulim double
;               Upper limit for the sampled quantities.
;
;       default double
;               Default value for coordinates falling outside of the
;               grid, or outside of the provided limits (llim and/or
;               ulim).  Set to 0.0 if not provided.
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
;       04 Dec 2014: (KBW) Original implementation, based on L.
;                          Coccato's 2D interpolation code.
;-
;------------------------------------------------------------------------------

PRO MDAP_PROPERTY_MAP_INTERPOLATE, $
                z_map, map_xs, map_dx, map_ys, map_dy, x_out, y_out, z_out, llim=llim, ulim=ulim, $
                default=default

        if n_elements(default) eq 0 then $  ; Define default if not provided
            default = 0.0d

        n = n_elements(x_out)               ; Number of data points to interpolate
        z_out=dblarr(n)                     ; Allocate the array

        ; Interpolate the output values based on the regular grid
        for i=0,n-1 do begin
            z_out[i] = bilinear(z_map, (x_out[i]-map_xs)/map_dx, (y_out[i]-map_ys)/map_dy, $
                                missing=default)
            if i lt 10 then $
                print, (x_out[i]-map_xs)/map_dx, (y_out[i]-map_ys)/map_dy, z_out[i]
            if n_elements(llim) ne 0 then begin
                if z_out[i] lt llim then $
                    z_out[i] = default
            endif
            if n_elements(ulim) ne 0 then begin
                if z_out[i] gt ulim then $
                    z_out[i] = default
            endif
        endfor
END

