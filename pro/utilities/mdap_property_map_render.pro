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

PRO MDAP_PROPERTY_MAP_RENDER_PIXELS, $
                x_samp, y_samp, z_samp, map_xs, map_dx, map_nx, map_ys, map_dy, map_ny, z_map, $
                default=default

        if n_elements(default) eq 0 then $                  ; Define default if not provided
            default = 0.0d

        z_map = make_array(map_nx, map_ny, /double, value=default)  ; Allocate space for map

        i = floor((x_samp-map_xs)/map_dx)
        j = floor((y_samp-map_ys)/map_dy)

        k = where(x_samp gt map_xs and i lt map_nx and y_samp gt map_ys and j lt map_ny, count)
        
;       if k[0] ne -1 then $
        if count ne 0 then $
            z_map[i[k],j[k]] = z_samp[k]
END


PRO MDAP_PROPERTY_MAP_ADD_APERTURE, $
                ap_x, ap_y, ap_z, ap_rx, ap_ry, z_map, n_map

        ; Get the size of the image
        sz = size(z_map)
        nx = sz[1]
        ny = sz[2]

        ; Get the search box
        sx = floor(ap_x-ap_rx-1)
        ex = ceil(ap_x+ap_rx+1)
        sy = floor(ap_y-ap_ry-1)
        ey = ceil(ap_y+ap_ry+1)

        for i=sx,ex do begin
            if sx lt 0 or ex ge nx then $
                continue
            ir = ((i - ap_x)/ap_rx)^2
            for j=sy,ey do begin
                if sy lt 0 or ey ge ny then $
                    continue

                if ir + ((j-ap_y)/ap_ry)^2 le 1.0d then begin
                    z_map[i,j] = z_map[i,j] + ap_z
                    n_map[i,j] = n_map[i,j] + 1
                endif
            endfor
        endfor
END


PRO MDAP_PROPERTY_MAP_RENDER_APERTURE, $
                x_samp, y_samp, z_samp, map_xs, map_dx, map_nx, map_ys, map_dy, map_ny, z_map, $
                aperture_radius, default=default

        z_map = dblarr(map_nx, map_ny)              ; Allocate space for map, force = 0.0
        nap = intarr(map_nx, map_ny)

        na = n_elements(x_samp)                     ; Number of apertures
        for i=0,na-1 do begin
            MDAP_PROPERTY_MAP_ADD_APERTURE, (x_samp[i]-map_xs)/map_dx, (y_samp[i]-map_ys)/map_dy, $
                                            z_samp[i], aperture_radius/map_dx, $
                                            aperture_radius/map_dy, z_map, nap
        endfor

        indx = where(nap gt 0, count, complement=nindx, ncomplement=ncount)
;       if indx[0] ne -1 then $
        if count ne 0 then $
            z_map[indx] = z_map[indx]/nap[indx]     ; Average the values

;       if n_elements(default) ne 0 and nindx[0] ne -1 then $
        if n_elements(default) ne 0 and ncount ne 0 then $
            z_map[nindx] = default                  ; Set default values for pixels with no data

END


PRO MDAP_PROPERTY_MAP_RENDER_TPS, $
                x_samp, y_samp, z_samp, map_xs, map_dx, map_nx, map_ys, map_dy, map_ny, z_map, $
                default=default

        z_map = GRID_TPS(x_samp, y_samp, z_samp, ngrid=[map_nx,map_ny], start=[map_xs,map_ys], $
                             delta=[map_dx,map_dy])

END


PRO MDAP_PROPERTY_MAP_RENDER, $
                x_samp, y_samp, z_samp, map_xs, map_dx, map_nx, map_ys, map_dy, map_ny, z_map, $
                default=default, aperture_radius=aperture_radius, tps=tps

        if keyword_set(tps) then begin
            MDAP_PROPERTY_MAP_RENDER_TPS, x_samp, y_samp, z_samp, map_xs, map_dx, map_nx, map_ys, $
                                          map_dy, map_ny, z_map, default=default
            return
        endif

        if n_elements(aperture_radius) ne 0 then begin
            MDAP_PROPERTY_MAP_RENDER_APERTURE, x_samp, y_samp, z_samp, map_xs, map_dx, map_nx, $
                                               map_ys, map_dy, map_ny, z_map, aperture_radius, $
                                               default=default
            return
        endif

        MDAP_PROPERTY_MAP_RENDER_PIXELS, x_samp, y_samp, z_samp, map_xs, map_dx, map_nx, map_ys, $
                                         map_dy, map_ny, z_map, default=default

END


