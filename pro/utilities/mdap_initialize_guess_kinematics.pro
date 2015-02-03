;+
; NAME:
;       MDAP_INITIALIZE_GUESS_KINEMATICS
;
; PURPOSE:
;
; CALLING SEQUENCE:
;       MDAP_INITIALIZE_GUESS_KINEMATICS, nb, analysis_prior, star_kin_interp, gas_kin_interp, $
;                                         bin_indx, velocity_initial_guess, $
;                                         velocity_dispersion_initial_guess, star_kin_guesses, $
;                                         gas_kin_guesses, vrange=vrange, srange=srange
;
; INPUTS:
;      nb integer
;               Number of binned spectra.
;
;      analysis_prior string
;               String representation of the analysis prior to use.
;               This is only checked to see if it has non-zero length.
;               If it does, the procedure assumes that the
;               MDAP_INTERPOLATE_KINEMATICS has been used to produce an
;               interpolated set of kinematics to use for the guess
;               kinematics.
;
;      star_kin_interp dblarr[nb][2]
;               Array with stellar kinematics (velocity and velocity
;               dispersion) interpolated from a previously available
;               dataset to use as a prior.
;
;      gas_kin_interp dblarr[nb][2]
;               Array with gas kinematics (velocity and velocity
;               dispersion) interpolated from a previously available
;               dataset to use as a prior.
;
;      bin_indx dblarr[]
;               Array that indicates which DRP spectra went into each
;               bin.
;
;      velocity_initial_guess double
;               Default value to use for the guess velocity.
;
;      velocity_dispersion_initial_guess double
;               Default value to use for the guess velocity dispersion
;               for the *stars*.  The guess velocity dispersion of the
;               gas is always taken to be 50 km/s.
;
; OPTIONAL INPUTS:
;      vrange dblarr[2]
;               The minimum (vrange[0]) and maximum (vrange[1]) to allow
;               for the guess velocities withe respect to the nominal
;               guess velocity (velocity_initial_guess).  So the limits
;               are actually = vrange+velocity_initial_guess!! Default:
;               vrange = [-500.,500.].
;
;      srange dblarr[2]
;               The minimum (srange[0]) and maximum (srange[1]) allowed
;               guess velocity dispersion.  Default: srange=[0.d,500d].
;
; OPTIONAL KEYWORDS:
;
; OUTPUT:
;      star_kin_guesses dblarr[nb][4]
;               Guess values for the stellar kinematics (v, sigma, h3,
;               h4) for each binned spectrum.
;
;      gas_kin_guesses dblarr[nb][2]
;               Guess values for the gas kinematics (v, sigma) for each
;               binned spectrum.
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
;       30 Jan 2015: (KBW) Pulled from manga_dap.pro
;-
;------------------------------------------------------------------------------

PRO MDAP_INITIALIZE_GUESS_KINEMATICS, $
                nb, analysis_prior, star_kin_interp, gas_kin_interp, bin_indx, $
                velocity_initial_guess, velocity_dispersion_initial_guess, star_kin_guesses, $
                gas_kin_guesses, vrange=vrange, srange=srange

        ; Set automatic vrange to +/- 500 km/s around the input guess velocity
        if n_elements(vrange) ne 2 then $
            vrange = [-500.0d, 500.0d]

        ; Set automatic srange to 0-500 km/s
        if n_elements(srange) ne 2 then $
            srange = [0.0d, 500.0d]

        ; Set the guess values for the stellar kinematics
        star_kin_guesses = dblarr(nb, 4)
        if strlen(analysis_prior) ne 0 and n_elements(star_kin_interp) ne 0 then begin

            for j=0,nb-1 do begin
                indx = where(bin_indx eq j, count)
                if count ne 0 then begin
                    star_kin_guesses[j,0] = mean(star_kin_interp[indx,0])   ; interpolated kin
                    star_kin_guesses[j,1] = mean(star_kin_interp[indx,1])
                endif else begin
                    star_kin_guesses[j,0] = velocity_initial_guess
                    star_kin_guesses[j,1] = velocity_dispersion_initial_guess
                endelse
            endfor

            indx = where(star_kin_guesses[*,0] lt velocity_initial_guess+vrange[0] or $
                         star_kin_guesses[*,0] gt velocity_initial_guess+vrange[1], count)
            if count ne 0 then $
                star_kin_guesses[indx,0] = velocity_initial_guess

            indx = where(star_kin_guesses[*,1] lt srange[0] or $
                         star_kin_guesses[*,1] gt srange[1], count)
            if count ne 0 then $
                star_kin_guesses[indx,1] = velocity_dispersion_initial_guess

        endif else begin
            ; Set the starting guesses for the kinematics to a
            ; single value
            star_kin_guesses[*,0] = velocity_initial_guess              ; stellar velocity
            star_kin_guesses[*,1] = velocity_dispersion_initial_guess   ; stellar sigma
        endelse


        ; Set the guess values for the gas kinematics
        ; TODO: Allow for gas dispersion limits to be different than stellar dispersion limits?
        gas_kin_guesses = dblarr(nb, 2)
        if strlen(analysis_prior) ne 0 and n_elements(gas_kin_interp) ne 0 then begin

            for j=0,nb-1 do begin
                indx = where(bin_indx eq j, count)
                if count ne 0 then begin
                    gas_kin_guesses[j,0] = mean(gas_kin_interp[indx,0])   ; interpolated kin
                    gas_kin_guesses[j,1] = mean(gas_kin_interp[indx,1])
                endif else begin
                    gas_kin_guesses[j,0] = velocity_initial_guess  ; gas velocity
                    gas_kin_guesses[j,1] = 50.                     ; gas sigma
                endelse
            endfor

            indx = where(gas_kin_guesses[*,0] lt velocity_initial_guess+vrange[0] or $
                         gas_kin_guesses[*,0] gt velocity_initial_guess+vrange[1], count)
            if count ne 0 then $
                gas_kin_guesses[indx,0] = velocity_initial_guess

            indx = where(gas_kin_guesses[*,1] lt srange[0] or $
                         gas_kin_guesses[*,1] gt srange[1], count)
            if count ne 0 then $
                gas_kin_guesses[indx,1] = 50.

        endif else begin
            ; Set the starting guesses for the kinematics to a
            ; single value
            gas_kin_guesses[*,0] = velocity_initial_guess      ; gas velocity
            gas_kin_guesses[*,1] = 50.                         ; gas sigma
        endelse

END


