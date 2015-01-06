;+
; NAME:
;       MDAP_VELOCITY_REGISTER
;
; PURPOSE:
;       Shift the input spectra by the provided velocities.  The input
;       wavelength vector is not changed.  Regions of the input spectra
;       shifted out of the input wavelength range are lost.  Unobserved
;       regions are set to have zero flux, unity inverse variance, and
;       masked.  If requested, the velocities are reset to differences
;       with respect to the median velocity.
;
;       TODO: Include a set of limits, such that velocities outside of
;       those limits are ignored and set to the median velocity offset?
;
; CALLING SEQUENCE:
;       MDAP_VELOCITY_REGISTER, wave, flux, ivar, mask, velocity, gflag=gflag, $
;                               /register_to_median, /bad_ivar
;
; INPUTS:
;       wave dblarr[C]
;               Wavelength of each of the C spectral channels.  This
;               vector does NOT change when registering the velocities
;               of the input spectra.
;
;       flux dblarr[N][C]
;               Flux in each of C spectral channels for each of N
;               spectra.  On output, the spectra will have been shifted
;               based on the input velocity vector.
;
;       ivar dblarr[N][C]
;               Inverse variance in each of C spectral channels for each
;               of N spectra.  On output, these variances will have been
;               shifted based on the input velocity vector.  The inverse
;               variances are based on a simple interpolation and do not
;               incorporate any error in that interpolation.
;
;       mask dblarr[N][C]
;               Bad pixel mask associated with the C spectral channels
;               in each of the N spectra.  These are shifted along with
;               flux and ivar.  Pixels that are interpolated to <0.5 are
;               considered to be unmasked.  Other pixels are masked.
;               Any pixel that was not in the original wavelength range
;               are masked.
;
;       velocity dblarr[N]
;               Estimated velocity (cz in km/s) of each spectrum.  The
;               input spectra are de-redshifted by the input value.  If
;               register_to_median is set, this vector contains the
;               offset velocities relative to the median.
;
; OPTIONAL INPUTS:
;       gflag intarr[N]
;               Flag that a spectrum is "good".  Any spectrum with
;               gflag=0 will be ignored.
;
; OPTIONAL KEYWORDS:
;       /register_to_median
;               Register the spectra to the median of the input
;               velocities.  If set, the velocity vector values are
;               offset by the median of the vector.
;
;       /bad_ivar
;               Mask bad (negative valued) inverse variance
;               interpolations.
;
; OUTPUT:
;       On output, flux, ivar, mask, and possibly velocity will have
;       changed based on the velocity offsets applied.
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
;       06 Dec 2014: (KBW) Original implementation
;-
;------------------------------------------------------------------------------

PRO MDAP_VELOCITY_REGISTER, $
                wave, flux, ivar, mask, velocity, gflag=gflag, $
                register_to_median=register_to_median, bad_ivar=bad_ivar

        ; Offset the velocities to the median value
        if keyword_set(register_to_median) then begin
            median_velocity = median(velocity)
            print, 'Median velocity offset: ', median_velocity
            velocity = velocity-median_velocity
        endif

        ; Convert velocity to redshift
        c=299792.458d                           ; Speed of light in km/s
        redshift = velocity/c

        ; Get the dimensions of the input data
        ; TODO: Provide some checks?
        sz = size(flux)
        ns = sz[1]                              ; Number of spectra
        nw = sz[2]                              ; Number of wavelength channels

        print, ns, nw

        ; Copy the input spectra for interpolation
        flux_ = flux
        ivar_ = ivar
        mask_ = mask

        reg = make_array(ns, /integer, value=1)
        if n_elements(gflag) eq ns then $
            reg = gflag

        for i=0,ns-1 do begin
            if reg[i] eq 0 then $               ; Flag to not register spectrum
                continue
            print, 'Velocity registering spectrum: ' + MDAP_STC(i+1, /integer) + '/' $
                   + MDAP_STC(ns, /integer) + ', velocity=' + MDAP_STC(velocity[i])

            ; Interpolate the spectrum at its deredshifted wavelengths
            ; to the input wavelength vector such that the wavelength
            ; vector does not change
            flux[i,*] = interpol(flux_[i,*], wave/(1.0d + redshift[i]), wave)
            ivar[i,*] = interpol(ivar_[i,*], wave/(1.0d + redshift[i]), wave)
            mask[i,*] = interpol(mask_[i,*], wave/(1.0d + redshift[i]), wave)

            ; Fix the mask to be either 1 (masked) or 0 (unmasked)
            indx = where(mask lt 0.5, count, complement=nindx, ncomplement=ncount)
;           if indx[0] ne -1 then $
            if count ne 0 then $
                mask[i,indx] = 0.0
;           if nindx[0] ne -1 then $
            if ncount ne 0 then $
                mask[i,nindx] = 1.0

            ; Mask bad ivar interpolations
            if keyword_set(bad_ivar) then begin
                indx = where(ivar le 0., count)
;               if indx[0] ne -1 then $
                if count ne 0 then $
                    mask[i,indx] = 1.0
            endif

            ; Find the extrapolated pixels and mask them
            indx = where(wave lt wave[0]/(1.0d + redshift[i]) $
                         or wave gt wave[nw-1]/(1.0d + redshift[i]), count)
;           if indx[0] ne -1 then begin
            if count ne 0 then begin
                flux[i,indx] = 0.0d
                ivar[i,indx] = 1.0d
                mask[i,indx] = 1.0d
            endif
        endfor

END

