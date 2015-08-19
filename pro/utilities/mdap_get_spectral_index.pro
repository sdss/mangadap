;+
; NAME:
;       MDAP_GET_SPECTRAL_INDEX
;
; PURPOSE:
;       Measure a spectral index, either defined as an absorption-line
;       index (using a blue and red continuum region and a primary
;       passband; see equations 2 and 3 from Worthey et al. 1994) or as
;       a spectral bandhead (only defined using a blue and red
;       passband).   For the former, errors are calculated based on both
;       a naive propagation of the error and the equations provided by
;       Cardiel et al. 1998.
;
; CALLING SEQUENCE:
;       MDAP_GET_SPECTRAL_INDEX, wave, flux, ivar, mask, abs_par, blue_cont, blue_cont_err, $
;                                red_cont, red_cont_err, equiv_width, equiv_width_err, $
;                                equiv_width_perr, index_mag, index_mag_err, index_mag_perr, $
;                                err=err, /geometric, /quiet
;
; INPUTS:
;       wave dblarr[S]
;               Wavelength of every pixel S.  The sampling does NOT need to be
;               linear in wavelength.
;
;       flux dblarr[S]
;               Flux at every pixel S.
;
;       ivar dblarr[S]
;               Inverse variance at every pixel S.
;
;       mask dblarr[S]
;               Bad pixel mask for each pixel S; 0-good, 1-bad.
;
;       abs_par AbsorptionIndex[]
;               Structure that contains the defining parameters of the
;               absorption-line index.
;
; OPTIONAL INPUTS:
;       err integer
;               Error flag for calculation, returned as 0 if calculation
;               is successful, 1 otherwise.
;
;
; OPTIONAL KEYWORDS:
;       /geometric
;               The sampling of the wavelength vector is geometric
;               (e.g., logarithmic); required by MDAP_BIN_EDGES to set
;               the edges of the wavelength pixels.
;
;       /quiet
;               Suppress output to STDOUT
;
; OUTPUT:
;       blue_cont double
;       blue_cont_err double
;       red_cont double
;       red_cont_err double
;               The blue and red pseudo-continuum measurements and their
;               propagated error.
;
;       equiv_width double
;       equiv_width_err double
;       equiv_width_perr double
;               Spectral index measured as an equivalent width in
;               angstroms (formula 2 in Worthey et al. 1994) and two
;               estimates of its error.  The first (_err) follows the
;               calculation provided by Cardiel et al. 1998); the second
;               (_perr) follows from a nominal propagation of error
;               through the calculation of the index.
;
;       index_mag double
;       index_mag_err double
;       index_mag_perr double
;               Spectral index measured in magnitudes (formula 3 in
;               Worthey et al. 1994) and two estimates of its error.
;               The first (_err) follows the calculation provided by
;               Cardiel et al. 1998); the second (_perr) follows from a
;               nominal propagation of error through the calculation of
;               the index.
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
;       MDAP_INDEX_IS_BANDHEAD()
;
; INTERNAL SUPPORT ROUTINES:
;       MDAP_GET_SPECTRUM_BANDHEAD
;       MDAP_GET_ABSORPTION_LINE_INDEX
;       MDAP_GET_PSEUDOCONTINUUM
;       MDAP_INTEGRATE_PIXELIZED_VALUE
;       MDAP_BIN_EDGES()
;
; REVISION HISTORY:
;       04 Nov 2014: Adapted from L. Coccato's version (v0.8) by K.
;                    Westfall (KBW)
;       11 Aug 2014: (KBW) Fixed a minor (but inconsequential) bug in
;                          MDAP_BIN_EDGES.  Restructured
;                          MDAP_GET_ABSORPTION_LINE_INDEX and allowed it
;                          to return the naive propagated error, instead
;                          of just the Cardiel-based error.  Adjusted
;                          MDAP_GET_SPECTRAL_INDEX to also return the
;                          propagated error.  MDAP_GET_SPECTRUM_BANDHEAD
;                          now uses MDAP_GET_PSEUDOCONTINUUM to get the
;                          flux weighted integral over the band, just as
;                          done for calculating the pseudocontinua on
;                          either side of the absorption-line index
;                          passband.  Added magnitude calculation to
;                          MDAP_GET_SPECTRUM_BANDHEAD, defined as the
;                          color between the two passbands.
;                          MDAP_GET_SPECTRAL_INDEX,
;                          MDAP_GET_ABSORPTION_LINE_INDEX, and
;                          MDAP_GET_SPECTRUM_BANDHEAD all now return the
;                          pseudocontinua as well and their (propagated)
;                          errors.  Added /quiet keywords.  Corrected
;                          some of the above documentation.
;       17 Aug 2014: (KBW) Moved MDAP_INTEGRATE_PIXELIZED_VALUE and
;                          MDAP_GET_PSEUDOCONTINUUM to their own files:
;                          utilities/mdap_integrate_pixelized_value.pro
;                          and utilities/mdap_get_pseudocontinuum.pro.
;                          Also moved MDAP_BIN_EDGES to be a support
;                          function in
;                          utilities/mdap_integrate_pixelized_value.pro.
;-
;------------------------------------------------------------------------------

;-------------------------------------------------------------------------------
; Bin edge [0] is the lower edge of pixel x[0], the upper edge of pixel x[0] is
; bin edge [1], etc.  The values of x are expected to be linearly or
; geometrically spaced.
;FUNCTION MDAP_BIN_EDGES, $
;                x, geometric=geometric
;
;        ; Get the number of bin centers
;        n = n_elements(x)
;        if n eq 1 then $
;            message, 'Cannot determine edges of a single bin!'
;
;        ; Convert the bin centers to their logarithm if the bins are
;        ; logarithmically spaced
;        xx = x
;        if keyword_set(geometric) then $
;            xx = alog(x)
;
;        ; Calculate the n+1 edges of the bins
;        edges = ([ xx[0]-(xx[1]-xx[0])/2.0d , xx[0:n-2] + (xx[1:n-1]-xx[0:n-2])/2.0, $
;                   xx[n-1]+(xx[n-1]-xx[n-2])/2.0d ])
;
;        ; Convert them back to linear space, if necessary
;        if keyword_set(geometric) then $
;            edges = exp(edges)
;
;        ; Return the bin edges
;        return, edges
;END

;-------------------------------------------------------------------------------
; Performs the integral:
;
;       integral = \int_xrange[0]^xrange[1] y dx
;
; using a Riemann sum, where x is expected to be at the center of the pixel.
; The error is estimated by a simple propagation of the error (i.e. error in the
; sum).
;
; x is expected to be sorted and contiguous; x coordinates are assumed to be
; linearly or geometrically spaced; the provided coordinate is expected to be at
; the (geometric) center of the pixel.
;PRO MDAP_INTEGRATE_PIXELIZED_VALUE, $
;                x, y, ye, mask, xrange, integral, integral_err, err=err, geometric=geometric
;
;        err=0                                                   ; Initialize error
;        n = n_elements(x)                                       ; Number of pixels
;        bin_edges = MDAP_BIN_EDGES(x, geometric=geometric)      ; The n+1 edges of the n pixels
;        if xrange[1] lt bin_edges[0] || xrange[0] gt bin_edges[n] then begin    ; No overlap
;;           print, xrange[0], xrange[1], bin_edges[0], bin_edges[n]
;            integral = 0.0d
;            integral_err = 1.0d
;            err=1                                               ; Return an error
;            return
;        endif
;
;        ; Set the values to 0.0 for masked pixels
;        indx = where(mask gt 0.5, count)
;        ym = y
;        yme = ye
;        if count ne 0 then begin
;            ym[indx] = 0.0d
;            yme[indx] = 0.0d
;        endif
;
;        ; Interval is smaller than the nearest pixel
;        indx = where(bin_edges gt xrange[0] and bin_edges lt xrange[1], ni)
;        if ni eq 0 then begin
;            indx = where(bin_edges gt xrange[1])
;            integral = ym[indx[0]-1]*(xrange[1]-xrange[0])
;            integral_err = yme[indx[0]-1]*(xrange[1]-xrange[0])
;            return
;        endif
;
;        ; Initialize output
;        integral = 0.0d
;        integral_err = 0.0d
;
;        ; Add partial pixel at the blue end
;        if indx[0] gt 0 then begin
;            integral = integral + ym[indx[0]-1]*(bin_edges[indx[0]]-xrange[0])
;            integral_err = integral_err + (yme[indx[0]-1]*(bin_edges[indx[0]]-xrange[0]))^2
;        endif
;
;        ; Add full pixels
;        if ni ge 2 then begin
;            integral = integral + total( ym[indx[0:ni-2]]*(bin_edges[indx[1:ni-1]] $
;                                                          - bin_edges[indx[0:ni-2]]) )
;            integral_err = integral_err + total( (yme[indx[0:ni-2]]*(bin_edges[indx[1:ni-1]] $
;                                                                     - bin_edges[indx[0:ni-2]]))^2 )
;        endif
;
;        ; Add partial pixel at the red end
;        if indx[ni-1] lt n then begin
;            integral = integral + ym[indx[ni-1]]*(xrange[1]-bin_edges[indx[ni-1]])
;            integral_err = integral_err + ( yme[indx[ni-1]]*(xrange[1]-bin_edges[indx[ni-1]]) )^2
;        endif
;
;        integral_err = sqrt(integral_err)
;END

;-------------------------------------------------------------------------------
; Measure the pseudo-continuum of the band using:
;
;                       \int_l0^l1 flux dl
;       continuum =   ---------------------
;                         \int_l0^l1 dl
;
; while accounting for the masked pixels (by not integrating over them) and
; providing a formal estimate of the error.
;PRO MDAP_GET_PSEUDOCONTINUUM, $
;                wave, flux, ivar, mask, passband, pseudo_continuum, pseudo_continuum_error, $
;                wave_integral=wave_integral, err=err, geometric=geometric
;
;        n = n_elements(wave)                            ; Number of pixels
;        unity = make_array(n, /double, value=1.0d)      ; Used to get denominator (=wave_integral)
;
;        ; Ignore masked regions or regions where the inverse variance is undefined
;        sige = sqrt(1.0d/ivar)
;        indx = where(ivar le 0.0d, count)
;        mask_ = mask
;        if count ne 0 then begin
;            mask_[indx] = 1.0d
;            sige[indx] = 1.0d
;        endif
;
;        MDAP_INTEGRATE_PIXELIZED_VALUE, wave, unity, unity, mask_, passband, wave_integral, dummy, $
;                                        err=err, geometric=geometric
;        ; MDAP_INTEGRATE_PIXELIZED_VALUE returns err=1 if the passband
;        ; did not overlap with the provided wavelength vector.
;        if err eq 1 then begin
;            pseudo_continuum = -9999.0d
;            pseudo_continuum_error = -9999.0d
;            wave_integral = -9999.0d
;            return
;        endif
;
;        MDAP_INTEGRATE_PIXELIZED_VALUE, wave, flux, sige, mask_, passband, flux_integral, $
;                                        flux_integral_err, geometric=geometric
;        ; If err != 1 above, then it won't equal zero here
;        pseudo_continuum = flux_integral/wave_integral
;        pseudo_continuum_error = flux_integral_err/wave_integral
;END


;-------------------------------------------------------------------------------
; _err is the error from Cardiel et al. 1998
; _perr is the naive propagation of the error using the provided ivar vector
PRO MDAP_GET_ABSORPTION_LINE_INDEX, $
                wave, flux, ivar, mask, passband, blueband, redband, blue_cont, blue_cont_err, $
                red_cont, red_cont_err, equiv_width, equiv_width_err, equiv_width_perr, index_mag, $
                index_mag_err, index_mag_perr, err=err, geometric=geometric, quiet=quiet

;       print, 'abs:', wave[0], wave[n_elements(wave)-1]
;       print, 'abs:', passband
;       print, 'abs:', blueband
;       print, 'abs:', redband
      
        ;---------------------------------------------------------------
        ; Check that there are unmasked pixels in all three bands

        ; Ignore masked regions or regions where the inverse variance is undefined
        mask_ = mask
        mindx = where(ivar le 0.0d, mcount)
        if mcount ne 0 then $
            mask_[mindx] = 1.0d

        ; Blue band
        bindx = where(wave ge blueband[0] and wave le blueband[1] and mask_ lt 0.5, count)
        if count eq 0 then begin
            if ~keyword_set(quiet) then $
                print, 'Blue wavelengths unavailable!'
            err = 1
            return
        endif

        ; Red band
        rindx = where(wave ge  redband[0] and wave le  redband[1] and mask_ lt 0.5, count)
        if count eq 0 then begin
            if ~keyword_set(quiet) then $
                print, 'Red wavelengths unavailable!'
            err = 1
            return
        endif

        ; Main passband
        pindx = where(wave ge passband[0] and wave le passband[1] and mask_ lt 0.5, count)
        if count eq 0 then begin
            if ~keyword_set(quiet) then $
                print, 'Passband wavelengths unavailable!'
            err = 1
            return
        endif
        ;---------------------------------------------------------------

;       print, 'abs: ', sn_blue
;       print, 'abs: ', nblue
;       print, 'abs: ', bcen
;       print, 'abs: ', blue_cont
;       print, 'abs: ', blue_cont_err
;       print, 'err: ', err

;       print, 'abs: ', rindx
;       print, 'abs: ', sn_red
;       print, 'abs: ', nred
;       print, 'abs: ', rcen
;       print, 'abs: ', red_cont
;       print, 'abs: ', red_cont_err
;       print, 'err: ', err

        ;---------------------------------------------------------------
        ; Get the pseudocontinuum in the off bands
        ; Blue band
        MDAP_GET_PSEUDOCONTINUUM, wave, flux, ivar, mask, blueband, blue_cont, blue_cont_err, $
                                  wave_integral=blue_width, err=err, geometric=geometric
        if err eq 1 then $
            return                              ; No wavelength integral so return with error

        ; Red band
        MDAP_GET_PSEUDOCONTINUUM, wave, flux, ivar, mask, redband, red_cont, red_cont_err, $
                                  wave_integral=red_width, err=err, geometric=geometric

        if err eq 1 then $
            return                              ; No wavelength integral so return with error
        ;---------------------------------------------------------------


        ; Calculate the slope and intercept of the continuum line
        bcen = (blueband[0] + blueband[1])/2.0                  ; Blue band center
        rcen = (redband[0] + redband[1])/2.0                    ; Red band center
        pcen = (passband[0] + passband[1])/2.0                  ; Passband center
        slope = (red_cont - blue_cont) / (rcen - bcen)
        intercept = blue_cont - bcen*slope

        ; Calculate the continuum level at all wavelengths
        continuum = slope * wave + intercept
        cindx = where(abs(continuum) eq 0.0, ccount)            ; Pixels with continuum == 0
        p1 = wave/(rcen-bcen) - bcen/(rcen-bcen)
        continuum_err = sqrt( ( blue_cont_err * (1.0d - p1))^2 + ( red_cont_err * p1 )^2 )

        ; Get the equivalent width integral (index in angstrom); equation 2 Worthey et al. 1994
        integrand = 1.0d - flux / continuum                     ; EW integrand
;       integrand_err = abs(1.0d / continuum / sqrt(ivar))      ; And its error
        integrand_err = sqrt( 1.0d / continuum^2 / ivar + (continuum_err * flux / continuum^2)^2 )

        ; Mask regions where the continuum is zero
        if ccount gt 0 then begin
            integrand[cindx] = 0.0d
            integrand_err[cindx] = 1.0d
            mask_[cindx] = 1.0
        endif

        ; Ignore masked regions or regions where the inverse variance is undefined
        if mcount ne 0 then $
            integrand_err[mindx] = 1.0d

        ; Determine the equivalent width of the line
        MDAP_INTEGRATE_PIXELIZED_VALUE, wave, integrand, integrand_err, mask_, passband, $
                                        equiv_width, equiv_width_perr, err=err, geometric=geometric
        if err eq 1 then $                              ; No overlap so return in error
            return

        ; Get the magnitude-based integral; equation 3 Worthey et al. 1994
        integrand = flux / continuum                            ; EW integrand
;       integrand_err = abs(1.0d / continuum / sqrt(ivar))
        if ccount gt 0 then $
            integrand[cindx] = 0.0d
        MDAP_INTEGRATE_PIXELIZED_VALUE, wave, integrand, integrand_err, mask_, passband, $
                                        index_mag, index_mag_perr, geometric=geometric
        
        ; Convert to a magnitude
        n = n_elements(wave)                            ; Number of pixels
        unity = make_array(n, /double, value=1.0d)
        MDAP_INTEGRATE_PIXELIZED_VALUE, wave, unity, unity, mask_, passband, pass_width, dummy, $
                                        geometric=geometric
        ; !! pass_width cancels out in the calculation of the propagated error !!
        index_mag_perr = abs(2.5 * index_mag_perr/index_mag/alog(10.0d))    ; Error done first
        index_mag = -2.5 * alog10(index_mag/pass_width)                     ; Index in magnitudes


        ;---------------------------------------------------------------
        ; Calculate the errors following Cardiel et al. 1998, eqns 41-46

        ; Get the S/N and number of elements in the different bands
        sn_blue = total(flux[bindx]*sqrt(ivar[bindx]))  ; Approximate S/N in the blue region
        nblue = n_elements(bindx)
        sn_red = total(flux[rindx]*sqrt(ivar[rindx]))   ; Approximate S/N in the red region
        nred = n_elements(rindx)
        sn_pass = total(flux[pindx]*sqrt(ivar[pindx]))  ; Approximate S/N in the main passband
        npass = n_elements(pindx)

        ; Their equation (eqn 42) assuming a constant pixel size is:
        ;
        ;  SN_cardiel = (sn_blue+sn_red+sn_pass)/(nblue+nred+npass)/sqrt(dw)
        ;
        ; Instead, I have allowed for a variable pixel size.  Over the band, the
        ; pixel size probably does not vary that much such that one can approximate dw ~
        ; Dw/N, such that equation 42 becomes

;        SN_cardiel = (sn_blue+sn_red+sn_pass) * sqrt( (nblue+nred+npass) $
;                                                      / (blue_width+red_width+pass_width) )
        SN_cardiel = sn_blue/sqrt(nblue*blue_width) + sn_red/sqrt(nred*red_width) $
                     + sn_pass/sqrt(npass*pass_width)
        
        q1=(rcen-pcen)/(rcen-bcen)
        q2=(pcen-bcen)/(rcen-bcen)
        c2= sqrt(1./pass_width + q1^2 / blue_width + q2^2 / red_width)          ; formula 44

        c1=pass_width*c2                                                        ; formula 43
        equiv_width_err=(c1-c2*equiv_width)/SN_cardiel                          ; formula 41

        c3=2.5*c2*alog10(exp(1))                                                ; formula 46
        index_mag_err=c3/SN_cardiel                                             ; formula 45
        ;---------------------------------------------------------------

END

; Calculate the spectral index as a bandhead
PRO MDAP_GET_SPECTRUM_BANDHEAD, $
                wave, flux, ivar, mask, blueband, redband, blue_cont, blue_cont_err, red_cont, $
                red_cont_err, bandhead, bandhead_err, bandhead_mag, bandhead_mag_err, err=err, $
                geometric=geometric, quiet=quiet

        bandhead = -9999.0d
        bandhead_err = -9999.0d
        bandhead_mag = -9999.0d
        bandhead_mag_err = -9999.0d

        err = 0 
        MDAP_GET_PSEUDOCONTINUUM, wave, flux, ivar, mask, blueband, blue_cont, blue_cont_err, $
                                  err=err, geometric=geometric
        if err eq 1 then $
            return                              ; No wavelength integral so return with error
        MDAP_GET_PSEUDOCONTINUUM, wave, flux, ivar, mask, redband, red_cont, red_cont_err, $
                                  err=err, geometric=geometric
        if err eq 1 then $
            return                              ; No wavelength integral so return with error

        ; Measure the strength of the bandhead and its error
        bandhead = red_cont/blue_cont
        bandhead_err = sqrt( (red_cont_err/blue_cont)^2 + (bandhead*blue_cont_err/blue_cont)^2 )

        if bandhead lt 0.0 then begin
            if ~keyword_set(quiet) then $
                print, 'Spectrum bandhead less than zero.  Cannot convert to magnitude!'
            err = 1
            return
        endif

        ; Convert the bandhead to a color in magnitudes
        bandhead_mag = -2.5 * alog10(bandhead)                          ; Index in magnitudes
        bandhead_mag_err = abs(2.5 * bandhead_err/bandhead/alog(10.0d)) ; ... and error

;       print, 'bandhead: ', bandhead, bandhead_err

END

; Primary wrapper routine
PRO MDAP_GET_SPECTRAL_INDEX, $
                wave, flux, ivar, mask, abs_par, blue_cont, blue_cont_err, red_cont, red_cont_err, $
                equiv_width, equiv_width_err, equiv_width_perr, index_mag, index_mag_err, $
                index_mag_perr, err=err, geometric=geometric, quiet=quiet

        is_bandhead = MDAP_INDEX_IS_BANDHEAD(abs_par)
        if is_bandhead eq 1 then begin
;            print, 'IS BANDHEAD'
            MDAP_GET_SPECTRUM_BANDHEAD, wave, flux, ivar, mask, abs_par.blue_cont, $
                                        abs_par.red_cont, blue_cont, blue_cont_err, red_cont, $
                                        red_cont_err, equiv_width, equiv_width_err, index_mag, $
                                        index_mag_err, err=err, geometric=geometric, quiet=quiet
            equiv_width_perr = equiv_width_err
            index_mag_perr = index_mag_err
        endif else begin
;            print, 'IS ABSORPTION'
            MDAP_GET_ABSORPTION_LINE_INDEX, wave, flux, ivar, mask, abs_par.passband, $
                                            abs_par.blue_cont, abs_par.red_cont, blue_cont, $
                                            blue_cont_err, red_cont, red_cont_err, equiv_width, $
                                            equiv_width_err, equiv_width_perr, index_mag, $
                                            index_mag_err, index_mag_perr, err=err, $
                                            geometric=geometric, quiet=quiet
        endelse

END
        
; Calculate the spectral index as a bandhead
;PRO MDAP_GET_SPECTRUM_BANDHEAD, $
;                wave, flux, ivar, mask, blue_cont, red_cont, bandhead, bandhead_err, err=err
;
;        err = 0 
;
;        ; Unmasked pixels in the red band
;        rindx = where(mask lt 0.5 and wave ge red_cont[0] and wave le red_cont[1], count)
;        if count eq 0 then begin                            ; No red pixels
;            err = 1                                             ; Set error and return
;            return
;        endif
;        rcont = median(flux[rindx], /even)                      ; Get the median continuum value
;        rivar = median(ivar[rindx], /even)                      ; Approximate error in the median
;;       print, 'bandhead: ', rcont, rivar
;
;        ; Unmasked pixels in the blue band
;        bindx = where(mask lt 0.5 and wave ge blue_cont[0] and wave le blue_cont[1], count)
;        if count eq 0 then begin                            ; No blue pixels
;            err = 1                                             ; Set error and return
;            return
;        endif
;        bcont = median(flux[bindx], /even)                      ; Get the median continuum value
;        bivar = median(ivar[bindx], /even)                      ; Approximate error in the median
;;       print, 'bandhead: ', bcont, bivar
;
;        bandhead = rcont/bcont                          ; Measure the strength of the bandhead
;        bandhead_err = sqrt( (bandhead/rcont)^2/rivar + (bandhead/bcont)^2/bivar ) ; ... and error
;
;;       print, 'bandhead: ', bandhead, bandhead_err
;
;END


;-------------------------------------------------------------------------------
; INTEGRAL MEAN OF THE BAND - the standard definition

;Inputs:
; lam [array] wavelength values of the pseudocontinua
; spc [array] spectra in the wav range of the pseudocontinua

;Outputs:
; pseudo_cont [float] value of the pseudoncontinuum
; pseudo_fit  [array] vector with as many elements as "lam" and equal
;                     to pseudo_cont

;PRO MDAP_GET_PSEUDOCONTINUUM_LODO, $
;                wave, flux, pseudo_cont, pseudo_fit
; 
;        ; FORMULA (1) Worthey et al. 1994
;        nwave = n_elements(wave)
;        pseudo_cont = int_tabulated(wave, flux, /double) / (wave[nwave-1]-wave[0])
;        pseudo_fit = make_array(nwave, /double, value=pseudo_cont)
;end
;
;;-------------------------------------------------------------------------------
;PRO MDAP_GET_ABSORPTION_LINE_INDEX_LODO, $
;                wave, flux, ivar, mask, passband, blue_cont, red_cont, equiv_width, $
;                equiv_width_err, index_mag, index_mag_err, title=title, plbound=plbound, plot=plot 
;                
;        ; TODO: What happens if there are intervening masks?  Interpolate across them?
;
;        ; TODO: The rebinning is meant to increase the accuracy of the integral,
;        ; but a Riemann sum can be done exactly.  It also changes the wavelength
;        ; sampling to be linear in wavelength.  Is this the primary reason this
;        ; resampling is done?
;
;        ; TODO: Need to mask ivar and include mask
;
;        bindx = where(wave ge blue_cont[0] and wave le blue_cont[1], nblue); Pixels in the blue band
;;       if bindx[0] eq -1 then $
;        if nblue eq 0 then $
;            message, 'Blue wavelengths unavailable!'
;;       nblue = n_elements(bindx)                                       ; Number of pixels
;        bcen = (blue_cont[0] + blue_cont[1])/2.0                        ; Blue band center
;
;        blue_width = (blue_cont[1]-blue_cont[0])                        ; Width of the blue band
;        dw = (blue_width)/(10*nblue-1)                                  ; Increase sampling
;        wave_resample = dindgen(10*nblue) * dw + blue_cont[0]           ; Set wavelength range
;        flux_resample = interpol(flux, wave, wave_resample, /spline)    ; Resample flux
;        ivar_resample = interpol(ivar, wave, wave_resample, /spline)    ; Resample the variance
;        sn_blue = total(flux_resample*sqrt(ivar_resample))              ; S/N in the blue region
;
;        ; Measure blue pseudo-continuum
;        MDAP_GET_PSEUDOCONTINUUM, wave_resample, flux_resample, blue_pseudocont, blue_pseudofit
;
;
;        rindx = where(wave ge  red_cont[0] and wave le  red_cont[1], nred)  ; Pixel in the red band
;;       if rindx[0] eq -1 then $
;        if nred eq 0 then $
;            message, 'Red wavelengths unavailable!'
;;       nred = n_elements(rindx)                                        ; Number of pixels
;        rcen = (red_cont[0] + red_cont[1])/2.0                          ; Red band center
;
;        red_width = (red_cont[1]-red_cont[0])                           ; Width of the red band
;        dw = red_width/(10*nred-1)                                      ; Increase sampling
;        wave_resample = dindgen(10*nred) * dw + red_cont[0]             ; Set wavelength range
;        flux_resample = interpol(flux, wave, wave_resample, /spline)    ; Resample flux
;        ivar_resample = interpol(ivar, wave, wave_resample, /spline)    ; Resample the variance
;        sn_red = total(flux_resample*sqrt(ivar_resample))               ; S/N in the red region
;
;        ; Measure red pseudo-continuum
;        MDAP_GET_PSEUDOCONTINUUM, wave_resample, flux_resample, red_pseudocont, red_pseudofit
;
;
;        pindx = where(wave ge passband[0] and wave le passband[1], npass)  ; Pixels in the band
;;       if pindx[0] eq -1 then $
;        if npass eq 0 then $
;            message, 'Passband wavelengths unavailable!'
;;       npass = n_elements(pindx)                                       ; Number of pixels
;        pcen = (passband[0] + passband[1])/2.0                          ; Passband center
;
;        pass_width = (passband[1]-passband[0])                          ; Width of the passband
;        dw = pass_width/(10*npass-1)                                    ; Increase sampling
;        wave_resample = dindgen(10*npass) * dw + passband[0]            ; Set wavelength range
;        flux_resample = interpol(flux, wave, wave_resample, /spline)    ; Resample flux
;        ivar_resample = interpol(ivar, wave, wave_resample, /spline)    ; Resample the variance
;        sn_pass = total(flux_resample*sqrt(ivar_resample))              ; S/N in the passband
;
;
;        ; Calculate the slope and intercept of the continuum line
;        slope = (red_pseudocont - blue_pseudocont) / (rcen - bcen)
;        intercept = blue_pseudocont - bcen*slope
;
;        continuum = slope * wave_resample + intercept                   ; Continuum level
;        integrand = 1. - flux_resample / continuum                      ; EW integrand
;
;        equiv_width = int_tabulated(wave_resample, integrand)           ; Equivalent width
;
;        index_mag = int_tabulated(wave_resample, flux_resample/continuum)
;        index_mag = -2.5 * alog10(index_mag/(passband[1]-passband[0]))  ; Index in magnitudes
;
;
;        ; Calculate the errors (See Cardiel et al. 1998, equations 41-46)
;        SN_cardiel = (sn_blue+sn_red+sn_pass)/(nblue+nred+npass)/sqrt(dw)       ; formula 42
;        
;        q1=(rcen-pcen)/(rcen-bcen)
;        q2=(pcen-bcen)/(rcen-bcen)
;        c2= sqrt(1./pass_width + q1^2 / blue_width + q2^2 / red_width)          ; formula 44
;
;        c1=pass_width*c2                                                        ; formula 43
;        equiv_width_err=(c1-c2*equiv_width)/SN_cardiel                          ; formula 41
;
;        c3=2.5*c2*alog10(exp(1))                                                ; formula 46
;        index_mag_err=c3/SN_cardiel                                             ; formula 45
;
;        if ~keyword_set(plot) then $
;            return                                              ; No plot, so return
;
;        ;-----------------------------------------------------------------------
;        ; Produce the plot -----------------------------------------------------
;        BASIC_COLORS, black, white, red, green, blue, yellow, cyan, magenta, orange, mint, purple, $
;                      pink, olive, lightblue, gray
;        ;UserSymbol,'Circle',/fill
;
;        if n_elements(title) eq 0 then title=''
;        if n_elements(plbound) eq 0 then plbound = [0.6,1.35]
;        midplot=interpol(flux, wave, 0.5*(min(wave)+max(wave)))
;
;        ; TODO: Where should always provide indicies given checks above?
;        inddd = where( wave ge blue_cont[0]-20 and wave le red_cont[1]+20 )
;        mn=min(flux(inddd))
;        mx=max(flux(inddd))
;        rgn=(mx-mn)*1.05
;
;        yrange=[mn-0.1*rgn,mx+0.1*rgn]
;        plot, wave, flux, /nodata, xtitle='!6wavelength [angstrom]', yrange=yrange, ystyle=1, $
;              title=title, xrange=[ blue_cont[0]-10,red_cont(1)+10], xstyle=1
;        oplot, wave, flux, thick=2
;        oplot, wave[bindx], flux[bindx], thick=3
;        oplot, wave[bindx], blue_pseudofit, color=blue, thick=5
;        plots, bcen, blue_pseudocont, psym=8
;
;        oplot, wave[rindx], flux[rindx], thick=3
;        oplot, wave[rindx], red_pseudofit, color=red, thick=5
;        plots, rcen, red_pseudocont, psym=8
;
;        oplot, [bcen, rcen], [blue_pseudocont, red_pseudocont], linestyle=1
;
;        oplot, wave[pindx], flux[pindx], thick=5
;
;        oplot, [blue_cont[0], blue_cont[0]], [0, interpol(flux, wave, blue_cont[0], /spline)], $
;               linestyle=1, color=blue
;        oplot, [blue_cont[1], blue_cont[1]], [0, interpol(flux, wave, blue_cont[1], /spline)], $
;               linestyle=1, color=blue
;        oplot, [red_cont[0], red_cont[0]], [0, interpol(flux, wave, red_cont[0], /spline)], $
;               linestyle=1, color=red
;        oplot, [red_cont[1], red_cont[1]], [0, interpol(flux, wave, red_cont[1], /spline)], $
;               linestyle=1,color=red
;        oplot, [passband[0], passband[0]], [0, interpol(flux, wave, passband[0], /spline)], $
;               linestyle=1
;        oplot, [passband[1], passband[1]], [0, interpol(flux, wave, passband[1], /spline)], $
;               linestyle=1
;
;        midplot=0.5*(mx+mn)
;        xyouts, min(wave[inddd])+21, mx+0.05*rgn, 'Dispersion=  ' + $
;                MDAP_ROUND_STR(wave[1]-wave[0],4) + ' (ang/pxl)'
;        xyouts, min(wave[inddd])+21, mx-0.00*rgn, 'Blue pscont= ' +MDAP_ROUND_STR(blue_pseudocont,4)
;        xyouts, min(wave[inddd])+21, mx-0.05*rgn, 'Red pscont=  ' + MDAP_ROUND_STR(red_pseudocont,4)
;        xyouts, min(wave[inddd])+21, mx-0.10*rgn, 'Eq. Width=    '+ MDAP_ROUND_STR(equiv_width,4) $
;                + ' +/- ' + MDAP_ROUND_STR(equiv_width_err,4) + ' (ang)'
;        xyouts, min(wave[inddd])+21, mx-0.15*rgn, 'Index=       '+ MDAP_ROUND_STR(index_mag,5) $
;                + ' +/- ' + MDAP_ROUND_STR(index_mag_err,5) + ' (mag)'
;        ;-----------------------------------------------------------------------
;
;END
;
