;+
; NAME:
;       MDAP_INTERPOLATE_KINEMATICS
;
; PURPOSE:
;       Given a DAP output file, read the desired kinematic product and
;       interpolate the predicted value at the provided on-sky
;       coordinates.
;
;       Nominal behavior is to base the interpolation on a regridded
;       interpolation of the binned data.  However, if there are too few
;       bins for usage of the GRID_TPS function in
;       MDAP_PROPERTY_MAP_RENDER, then the binned results are first
;       redistributed into the pre-binned coordinates using
;       MDAP_DISTRIBUTE_BINNED_DATA.  The map is then rendered and
;       interpolated from to get the data at the input xcoo and ycoo.
;
;       TODO: Need to test the approach when the number of bins is
;       small.
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
;       min_bins integer
;               Minimum number of bins before forcing the interpolation
;               to be based on a "redistribution" of the data.  This
;               MUST be 7 or larger.  Default value is 7.
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
;       - Function will not work if the provided prior is a RADIAL file.
;         Check for this and use "redistributed interpolation"?
;       - Impose limits on the interpolated values?
; 
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;       06 Dec 2014: Original implementation by K. Westfall
;       19 Feb 2015: (KBW) Fix behavior when small number (or one) bin
;                          in prior file.
;-
;------------------------------------------------------------------------------

PRO MDAP_INTERPOLATE_KINEMATICS, $
                dap_file, xcoo, ycoo, interpolated_kin, velocity=velocity, sigma=sigma, $
                stellar=stellar, min_bins=min_bins

        ; Check that there's something to return!
        if ~keyword_set(velocity) && ~keyword_set(sigma) then $
            message, 'No property to interplate!'

        ; Get/Check the minimum number of bins to use for the rendering
        if n_elements(min_bins) eq 0 then $
            min_bins = 7
        if min_bins lt 7 then begin
            print, 'WARNING: Minimum number of bins must be 7 or more.  Changing to 7.'
            min_bins = 7
        endif

        ; Get the kinematics to interpolate
        if keyword_set(stellar) then begin
            MDAP_DEFINE_OUTPUT, stellar_kinematics_fit=kin
            MDAP_READ_OUTPUT, dap_file, stellar_kinematics_fit=kin

            ; TODO: Probably need a more robust test here:
            ; Stellar kinematics aren't available!
            if (size(kin))[0] ne 2 then begin
                print, 'Stellar kinematics not available in prior DAP file!  Ignoring.'
                tempvar = temporary(kin)
            endif

        endif else begin

            ; Try to use the GANDALF results
            MDAP_DEFINE_OUTPUT, emission_line_kinematics_avg=kin
            MDAP_READ_OUTPUT, dap_file, emission_line_kinematics_avg=kin

            ; TODO: Probably need a more robust test here:
            ; If these aren't available, try to get Enci's
            if (size(kin))[0] ne 2 then begin
                MDAP_DEFINE_OUTPUT, elo_ew_kinematics_avg=kin
                MDAP_READ_OUTPUT, dap_file, elo_ew_kinematics_avg=kin
            endif

            ; TODO: Probably need a more robust test here:
            ; Gas kinematics aren't available!
            if (size(kin))[0] ne 2 then begin
                print, 'Gas kinematics not available in prior DAP file!  Ignoring.'
                tempvar = temporary(kin)
            endif

        endelse

        ; If not available, return without setting the interpolated kinematics
        if n_elements(kin) eq 0 then $
            return
;           message, 'No kinematics available to interpolate!'

        ; Select the appropriate kinematic property
        if keyword_set(velocity) then $
            cz = reform(kin[*,0])
        if keyword_set(sigma) then $
            sig = reform(kin[*,1])

        ; Get the bin positions
        if (size(kin))[1] lt min_bins then begin           ; Insufficient number of bins
            ; Get the original DRP x and y positions, and the bin index for each DRP spectrum
            MDAP_DEFINE_OUTPUT, xpos=xbin, ypos=ybin, bin_indx=bin_indx
            MDAP_READ_OUTPUT, dap_file, xpos=xbin, ypos=ybin, bin_indx=bin_indx

            if n_elements(xbin) lt min_bins then $
                message, 'Insufficient number of spectra in DRP file for interpolation!'

            if keyword_set(velocity) then begin
                MDAP_DISTRIBUTE_BINNED_DATA, cz, bin_indx, cz_drp
                cz = cz_drp
            endif
            if keyword_set(sigma) then begin
                MDAP_DISTRIBUTE_BINNED_DATA, sig, bin_indx, sig_drp
                sig = sig_drp
            endif
        endif else begin
            ; Use the provided bin positions directly
            MDAP_DEFINE_OUTPUT, xbin_rlow=xbin, ybin_rupp=ybin
            MDAP_READ_OUTPUT, dap_file, xbin_rlow=xbin, ybin_rupp=ybin
        endelse

        ; Get the spaxel size
        MDAP_DEFINE_OUTPUT, dx=dx, dy=dy
        MDAP_READ_OUTPUT, dap_file, dx=dx, dy=dy

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
            MDAP_PROPERTY_MAP_RENDER, xbin, ybin, cz, xs, dx, nx, ys, dy, ny, map, $
                                      default=mean(cz), /tps

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
            MDAP_PROPERTY_MAP_RENDER, xbin, ybin, sig, xs, dx, nx, ys, dy, ny, map, $
                                      default=mean(sig), /tps
            ; Interpolate the velocity at the input coordinates
            MDAP_PROPERTY_MAP_INTERPOLATE, map, xs, dx, ys, dy, xcoo, ycoo, interpolated_sigma
        endif
       
        ; Copy them to the returned array
        if keyword_set(velocity) && ~keyword_set(sigma) then $
            interpolated_kin = interpolated_velocity
        if ~keyword_set(velocity) && keyword_set(sigma) then $
            interpolated_kin = interpolated_sigma
        if keyword_set(velocity) && keyword_set(sigma) then $
            interpolated_kin = [ [interpolated_velocity], [interpolated_sigma] ]

END


