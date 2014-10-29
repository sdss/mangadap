;+
; NAME:
;       MDAP_UPDATE_STARTING_GUESSES
;
; PURPOSE:
;       Transform the kinematic results for one spatial binning inot starting
;       guesses for another.  The main module that performs the computation is
;       MDAP_INTERPOLATE_2DMAPS.
;
; CALLING SEQUENCE:
;       MDAP_UPDATE_STARTING_GUESSES, input_map, input_xbin, input_ybin, x2d, y2d, output_xbin, $
;                                     output_ybin, output_start_guess, velocity_initial_guess, $
;                                     velocity_dispersion_initial_guess, H3_initial_guess, $
;                                     H4_initial_guess, h3h4=h3h4
;
; INPUTS:
;       input_map dblarr[B][2 or 4]
;               Array with kinematic measurements for each of the N spatial
;               bins.  If the keyword /h3h4 is set, input_map must be Nx4
;               elements array, and previous measurements of h3 and h4 must be
;               given.  The input map columns are:
;
;                       input_map[*,0]: previous velocity measurements
;                       input_map[*,1]: previous velocity dispersion measurements
;                       input_map[*,2]: (optional) previous H3 measurements
;                       input_map[*,3]: (optional) previous H4 measurements
;
;       input_xbin dblarr[B]
;               Array with the X coordinates of the center of the N spatial
;               bins.
;
;       input_ybin dblarr[B]
;               Array with the Y coordinates of the center of the N spatial
;               bins.
;
;       dx double
;               Sampling size in X by interpolation algorithm,
;               MDAP_INTERPOLATE_2DMAPS.
;
;       dy double
;               Sampling size in Y by interpolation algorithm,
;               MDAP_INTERPOLATE_2DMAPS.
;
;       output_xbin dblarr[b]
;               X coordinates at which to interpolate new starting guesses.
;
;       output_ybin dblarr[b]
;               Y coordinates at which to interpolate new starting guesses.
;
;       def_velocity  double
;               Default velocity to use if previous measurements are NaN or are
;               not defined.
;
;       def_veldisp  double
;               Default velocity dispersion to use if previous measurements are
;               NaN, not defined, or <21 km/ or > 499 km/s.
;
;       dev_h3 double
;               Default h3 value to use if previous measurements are
;               NaN, not defined, or <-0.399 or >0.399
;
;       dev_h4 double
;               Default h4 value to use if previous measurements are
;               NaN, not defined, or <-0.399 or >0.399
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;
;       /h3h4
;               If set, the guesses of Gauss Hermite moments (h3, h4) will be
;               computed.  Otherwise, they are set to the default values.  If
;               calculation requested, the input_map must contain h3 and h4
;               measurements.
;
; OUTPUT:
;       output_map dblarr[b][2 or 4]
;               Interpolated velocity, velocity dispersion, h3 and h4 of the b
;               spectra to be used as starting guesses.  The h3 and h4 columns
;               are only included if /h3h4 is set.
;
; OPTIONAL OUTPUT:
;
; COMMENTS:
;
; EXAMPLES:
;
; TODO:
;       - Input a 4X2 array with the limits for the kinematic measurements.
;
; BUGS:
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;       23 Sep 2014: Copied from version 0_8 by L. Coccato
;       24 Sep 2014: (KBW) Formatting and edits; input spaxel size instead of
;                    original map size
;-
;------------------------------------------------------------------------------

PRO MDAP_IMPOSE_LIMITS, $
        vec, range

        ; Range must be two elements
        if n_elements(range) ne 2 then $
            message, 'Input range must have two elements.'

        ; Impose lower limit
        indx = where(vec lt range[0])
        if indx[0] ne -1 then $
            range[indx] = range[0]
       
        ; Impose upper limit
        indx = where(vec gt range[1])
        if indx[0] ne -1 then $
            range[indx] = range[1]
END
       

PRO MDAP_UPDATE_STARTING_GUESSES, $
                input_map, input_xbin, input_ybin, spaxel_dx, spaxel_dy, output_xbin, output_ybin, $
                def_velocity, dev_veldisp, output_map, def_h3=def_h3, def_h4=def_h4, h3h4=h3h4

        sz=size(input_map)
        if sz[0] ne 2 then $
            message, 'Input value array must have two dimensions!'
        if sz[2] ne 2 and sz[2] ne 4 then $
            message, 'Input value array must have 2 or 4 columns!'
        if keyword_set(h3h4) and sz[2] ne 4 then $
            message, 'Input value array must have 4 columns in order to create h3/h4 guesses.'

        if keyword_set(h3h4) and (n_elements(def_h3) eq 0 or n_elements(def_h4) eq 0) then $
            message, 'Must provided default h3 and h4 values if interpolating these properties.'

        ns = sz[1]                                      ; Number of spectra
        nc = sz[2]                                      ; Number of map columns

        ; Interpolate the data, use the default values if the interpolation
        ; cannot be done and impose hard-wired limits

        ; Velocity
        MDAP_INTERPOLATE_2DMAPS, input_xbin, input_ybin, input_map[*,0], spaxel_dx, spaxel_dy, $
                                 output_xbin, output_ybin, interp_velocity, default=def_velocity

        ; Velocity dispersion
        MDAP_INTERPOLATE_2DMAPS, input_map[*,1], input_xbin, input_ybin, skyx, skyy, $
                                 output_xbin, output_ybin, interp_veldisp, default=def_veldisp
        MDAP_IMPOSE_LIMITS, interp_veldisp, ([21., 499.])

        if keyword_set(h3h4) then begin
            ; H3
            MDAP_INTERPOLATE_2DMAPS, input_map[*,2], input_xbin, input_ybin, skyx, skyy, $
                                     output_xbin, output_ybin, interp_h3, default=def_h3
            MDAP_IMPOSE_LIMITS, interp_h3, ([-0.399., 0.399])

            ; H4
            MDAP_INTERPOLATE_2DMAPS, input_map[*,3], input_xbin, input_ybin, skyx, skyy, $
                                     output_xbin, output_ybin, interp_h4, default=def_h4
            MDAP_IMPOSE_LIMITS, interp_h4, ([-0.399., 0.399])
        endif

        ; Save the results
        output_map = dblarr(ns, nc)             ; Initialize the array
        output_map[*,0]=interp_velocity[*]      ; Velocity
        output_map[*,1]=interp_veldisp[*]       ; Velocity dispersion
        if keyword_set(h3h4) then begin
            output_map[*,2]=interp_h3[*]        ; H3
            output_map[*,3]=interp_h4[*]        ; H4
        endif

END

