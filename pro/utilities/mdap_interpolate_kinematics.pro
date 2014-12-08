;+
; NAME:
;       MDAP_INTERPOLATE_KINEMATICS
;
; PURPOSE:
;       Given a DAP output file, read the desired kinematic product and
;       interpolate the predicted value at the provided on-sky
;       coordinates.
;
;       TODO: Currently only interpolates the stellar or (average)
;       emission-line velocities.
;
; CALLING SEQUENCE:
;       MDAP_INTERPOLATE_KINEMATICS, dap_file, xcoo, ycoo, interpolated_kin, /velocity, /stellar
;
; INPUTS:
;       dap_file string
;               DAP output file with the kinematics to interpolate.
;
;       xcoo dblarr[N]
;               A set of on-sky X coordinates at which to sample the
;               kinematics in the provided DAP file.
;
;       ycoo dblarr[N]
;               A set of on-sky Y coordinates at which to sample the
;               kinematics in the provided DAP file.
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;       /velocity
;               Sample the velocity component.
;
;       /stellar
;               Use the stellar kinematics instead of the average
;               emission-line kinematics.
;
; OUTPUT:
;       interpolated_kin dblarr[N]
;               Kinematics interpolated to the provided on-sky
;               coordinates.
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
;       - Add other kinematic moments and check that the moments exist
;         in the input file.
;       - Allow for different interpolation schemes?
;       - Allow to output the rendered map?
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;       06 Dec 2014: (KBW) Original implementation
;-
;------------------------------------------------------------------------------

PRO MDAP_INTERPOLATE_KINEMATICS, $
                dap_file, xcoo, ycoo, interpolated_kin, velocity=velocity, sigma=sigma, $
                stellar=stellar

        ; Check that there's something to return!
        if ~keyword_set(velocity) and ~keyword_set(sigma) then $
            message, 'No property to interplate!'

        ; Read the bin positions and the requested kinematics
        if keyword_set(stellar) then begin
            MDAP_DEFINE_OUTPUT, xbin_rlow=xbin, dx=dx, ybin_rupp=ybin, dy=dy, $
                                stellar_kinematics_fit=kin
            MDAP_READ_OUTPUT, dap_file, xbin_rlow=xbin, dx=dx, ybin_rupp=ybin, dy=dy, $
                              stellar_kinematics_fit=kin
        endif else begin
            MDAP_DEFINE_OUTPUT, xbin_rlow=xbin, dx=dx, ybin_rupp=ybin, dy=dy, $
                                emission_line_kinematics_avg=kin
            MDAP_READ_OUTPUT, dap_file, xbin_rlow=xbin, dx=dx, ybin_rupp=ybin, dy=dy, $
                              emission_line_kinematics_avg=kin
        endelse

;       for i=0,9 do $
;           print, xbin[i], ybin[i], kin[i,0], kin[i,1]

        ; Create the coordinate grid for the interpolated map, ensuring
        ; that the grid covers both the input coordinate and the read
        ; bin coordinates.
        MDAP_COORDINATE_GRID, xbin, dx, xs, nx, aux=xcoo, fbuffer=1.05
        MDAP_COORDINATE_GRID, ybin, dy, ys, ny, aux=ycoo, fbuffer=1.05

;       print, xs, dx, nx, xs+dx*(nx-1)
;       print, ys, dy, ny, ys+dy*(ny-1)

        ; Get the interpolated values
        if keyword_set(velocity) then begin
            ; Render a gridded map of the velocity; default is to use
            ; the smooth GRID_TPS rendering
            MDAP_PROPERTY_MAP_RENDER, xbin, ybin, reform(kin[*,0]), xs, dx, nx, ys, dy, ny, map, $
                                      default=mean(kin[*,0]), /tps

;           print, size(map)
;           for i=0,9 do $
;               print, map[i,i]

            ; Interpolate the velocity at the input coordinates
            MDAP_PROPERTY_MAP_INTERPOLATE, map, xs, dx, ys, dy, xcoo, ycoo, interpolated_velocity

;           for i=0,9 do $
;               print, xcoo[i], ycoo[i], interpolated_velocity[i]

        endif

        if keyword_set(sigma) then begin
            ; Render a gridded map of the velocity dispersion; default is to use
            ; the smooth GRID_TPS rendering
            MDAP_PROPERTY_MAP_RENDER, xbin, ybin, reform(kin[*,1]), xs, dx, nx, ys, dy, ny, map, $
                                      default=mean(kin[*,1]), /tps
            ; Interpolate the velocity at the input coordinates
            MDAP_PROPERTY_MAP_INTERPOLATE, map, xs, dx, ys, dy, xcoo, ycoo, interpolated_sigma
        endif
       
        ; Copy them to the returned array
        if keyword_set(velocity) and ~keyword_set(sigma) then $
            interpolated_kin = interpolated_velocity
        if ~keyword_set(velocity) and keyword_set(sigma) then $
            interpolated_kin = interpolated_sigma
        if keyword_set(velocity) and keyword_set(sigma) then $
            interpolated_kin = [ [interpolated_velocity], [interpolated_sigma] ]


END


