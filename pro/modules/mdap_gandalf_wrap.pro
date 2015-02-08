;+
; NAME:
;       MDAP_GANDALF_WRAP
;
; PURPOSE:
;       Wrapper that calls PPXF and GANDALF in order to derive the stellar and
;       gaseous kinematics for an input set of spectra.  TODO: Why is this
;       specific to MaNGA?  The stellar continuum is matched with a combination
;       of provided templates, whereas emission-lines are represented by
;       Gaussian functions, with interdepencies regulated by the input
;       emission-setup file.
;
; CALLING SEQUENCE:
;       MDAP_GANDALF_WRAP, tpl_wave, tpl_flux, obj_wave, obj_flux, obj_ivar, obj_mask, velScale, $
;                          start, fitted_pixels_ppxf, weights_ppxf, add_poly_coeff_ppxf, $
;                          mult_poly_coeff_ppxf, bestfit_ppxf, chi2_ppxf, fitted_pixels_gndf, $
;                          weights_gndf, mult_poly_coeff_gndf, bestfit_gndf, chi2_gndf, eml_model, $
;                          sol, err, gas_intens, gas_intens_err, gas_vel, gas_vel_err, gas_sig, $
;                          gas_sig_err, gas_flux, gas_flux_err, gas_ew, gas_ew_err, $
;                          eml_par=eml_par, vsyst=vsyst, bias=bias, mdegree=mdegree, $
;                          degree=degree, moments=moments, reddening=ebv, $
;                          range_v_star=range_v_star, range_s_star=range_s_star, $
;                          range_v_gas=range_v_gas, range_s_gas=range_s_gas, $
;                          region_mask=region_mask, external_library=external_library, $
;                          int_disp=instr_disp, err_reddening=err_reddening, $
;                          ppxf_status=ppxf_status, gandalf_status=gandalf_status, $
;                          ppxf_only=ppxf_only, /oversample, /for_errors, /fix_star_kin, $
;                          /fix_gas_kin, /quiet, /plot
;
; INPUTS:
;       tpl_wave dblarr[S]
;               The wavelength in angstrom for each of S spectral channels,
;               which must be the same for all T template spectra.
;       
;       tpl_flux dblarr[T][S]
;               The fluxes for a library of T template spectra, each with S
;               spectral channels.  The sampling of the spectra is expected to
;               be logarithmic using base 10.  TODO: Always base 10?  All pixels
;               are expected to be valid! TODO: Allow for mask input?  The
;               velocity sampling of both the template and object spectra MUST
;               be the same.
;
;       obj_wave dblarr[C]
;               The wavelength in angstroms for each of C spectral channels in
;               the object spectrum.
;       
;       obj_flux dblarr[C]
;               A single object (galaxy) spectrum with C spectral channels.  The
;               sampling of the spectra is expected to be logarithmic using base
;               10.  The velocity sampling of both the template and object
;               spectra MUST be the same.
;
;       obj_ivar dblarr[C]
;               Inverse variance in the flux of the object spectrum.
;
;       obj_mask dblarr[C]
;               Bad pixel mask for the object spectrum.  Good/bad pixels have
;               values of 0/1.
;               
;       velScale double
;               Velocity sampling of the spectra.  Must be the same for both the
;               object and template spectra.
;
;               TODO: Better to just calculate this using MDAP_VELOCITY_SCALE?
;
;       start dblarr[6]
;               Initial guesses for the:
;                   start[0] stellar veocity (km/sec)
;                   start[1] stellar velocity dispersion (km/sec)
;                   start[2] stellar h3 Gauss Hermite moment
;                   start[3] stellar h4 Gauss Hermite moment stellar velocity dispersion
;                   start[4] gas velocity (km/sec)
;                   start[5] gas velocity dispersion (km/sec).
;
; OPTIONAL INPUTS:
;       eml_par EmissionLine[E]
;               The parameters for each of E emission lines used during the
;               fitting procedures.  The EmissionLine structure is defined as
;               follows (see MDAP_READ_EMISSION_LINE_PARAMETERS):
;
;               { EmissionLine, i:0L, name:'', lambda:0.0d, action:'', $
;                     kind:'', a:0.0d, v:0.0d, s:0.0d, fit:'' }
;
;               Once created, one selects, for example, the name of the 3rd
;               input line using: eml_par[2].name
;
;               If not provided, the code runs pPXF and then returns without
;               determining any gas paramters.
;
;       vsyst double
;               On input, pPXF expects the template and galaxy to have the same
;               wavelength coordinate system.  This accounts for the fact that
;               this may not be true, calculated (MDAP_SPECTRAL_FITTING) as the
;               velocity offset of the template with respect to the initial
;               wavelength of the galaxy.  That is, we determine the
;               pseudo-velocity shift required to match the observed wavelengths
;               of the template library and the galaxy spectrum.
;
;       bias float
;               (COPIED FROM pPXF:) This parameter biases the (h3,h4,...)
;               measurements towards zero (Gaussian LOSVD) unless their
;               inclusion significantly decreses the error in the fit.  Set this
;               to BIAS=0.0 not to bias the fit: the solution (including
;               [V,sigma]) will be noisier in that case.  The default BIAS
;               should provide acceptable results in most cases, but it would be
;               safe to test it with Monte Carlo simulations. This keyword
;               precisely corresponds to the parameter \lambda in the Cappellari
;               & Emsellem (2004) paper. Note that the penalty depends on the
;               *relative* change of the fit residuals, so it is insensitive to
;               proper scaling of the NOISE vector. A nonzero BIAS can be safely
;               used even without a reliable NOISE spectrum, or with equal
;               weighting for all pixels.
;
;       mdegree integer
;               Degree of *multiplicative* polynomial to be used by both  pPXF
;               and GANDALF.  GANDALF will only use this polynomial if the
;               reddening is NOT fit.  Default: mdegree=0 (no multiplicative
;               polynomials are used).
;
;       degree integer
;               Degree of *additive* (TODO: Check this!) polynomial to be used
;               in pPXF only.  Default: degree=-1 (no additive polynomials are
;               used).
;
;       moments integer
;               (Copied from pPXF:) Order of the Gauss-Hermite moments to fit.
;               Set this keyword to 4 to fit [h3, h4] and to 6 to fit [h3, h4,
;               h5, h6].  Note that in all cases the G-H moments are fitted
;               (nonlinearly) *together* with [V, sigma].  Default behavior is
;               to set moments=2 such that only [V, sigma] are fitted and the
;               other parameters are returned as zero.  If moments=0 then only
;               the templates and the continuum additive polynomials are fitted
;               and the WEIGHTS are returned in output.
;
;       reddening dblarr[1 or 2]
;               If it exists upon input, the stellar reddening (if it is
;               dblarr[1]) and the gas reddening (balmer decrement, if it is
;               dblarr[2]) are fit.  On output, it is replaced by the best-fit
;               reddening.
;
;       range_v_star dblarr[2]
;               Lower [0] and upper [1] limits on the best-fitting stellar
;               velocity in km/s.  Default values are the starting_guess +/-
;               2000 km/s.
;
;       range_s_star dblarr[2]
;               Lower [0] and upper [1] limits on the best-fitting stellar
;               velocity dispersion in km/s.  Default values are the
;               21<sigma<499 km/s.
;
;       range_v_gas dblarr[2]
;               Lower [0] and upper [1] limits on the best-fitting gas velocity
;               in km/s.  Default values are the starting_guess +/- 2000 km/s.
;
;       range_s_gas dblarr[2]
;               Lower [0] and upper [1] limits on the best-fitting gas velocity
;               dispersion in km/s.  Default values are the starting_guess +/-
;               2000 km/s. TODO: Are these limits right?
;
;       region_mask dblarr[2*N]
;               If defined, it provides specifies a set of wavelength ranges
;               (lower and upper limits in angstroms) to omit from the fit.  It
;               must contain an even number of entries such that the ranges are:
;                   lower_wave_limit[i] = region_mask[2*i]
;                   upper_wave_limit[i] = region_mask[2*i+1]
;               for i=0..N-1.
;
;               TODO: Omit this and use the input pixel mask?
;
;       external_library string
;               Path to the external FORTRAN library, which contains the fortran
;               versions of mdap_bvls.pro.  If not specified, or if the path is
;               invalid, the default internal IDL mdap_bvls code is used. 
;
;       int_disp dblarr[E]
;               Instrumental velocity dispersion (in km/s) for all E emission
;               lines measured at the observed wavelength of the object spectra.
;               For MDAP_SPECTRAL_FITTING, this is genereated using the rest
;               wavelengths from an input file and the gas starting velocity
;               guess.
;
;       ppxf_only integer
;               Flag to run the GANDALF wrapper only using PPXF to fit the stellar
;               kinematics and optimal template.  Emission lines will be masked
;               if eml_par is provided; however, the emissions will not be fit
;               using GANDALF.  0/1 = false/true.
;
; OPTIONAL KEYWORDS:
;       /oversample 
;               (Copied from pPXF:) Set this keyword to oversample the template
;               by a factor 30x before convolving it with a well sampled LOSVD.
;               This can be useful to extract proper velocities, even when sigma
;               < 0.7*velScale and the dispersion information becomes totally
;               unreliable due to undersampling.  IMPORTANT: One should sample
;               the spectrum more finely is possible, before resorting to the
;               use of this keyword! 
;
;       /for_errors
;               Compute errors in the emission-line data using GANDALF.  TODO:
;               This is now "Mandatory for the DAP workflow."  Should it be?
;
;       /fix_star_kin
;               Fix the stellar kinematics to the starting guesses during the
;               fit.
;
;       /fix_gas_kin
;               Fix the gas kinematics to the starting guesses during the fit.
;
;       /quiet
;               Execute GANDALF in quiet mode.  TODO: pPXF if always quiet.
;               Allow p_quiet and g_quiet?  Or quiet/not quiet for both?
;
;       /plot
;               Produce the plots generated by PPXF and GANDALF
;
; OUTPUT:
;       fitted_pixels_ppxf intarr[]
;               Indices of the pixels used by PPXF.
;
;       weights_ppxf dblarr[T]
;               The weights applied to each template spectrum to create the
;               optimal composite template as determined by PPXF.
;
;       add_poly_coeff_ppxf dblarr[P]
;               Legendre polynomial weights of the additive polynomial if used in PPXF.
;
;       mult_poly_coeff_ppxf dblarr[P]
;               Legendre polynomial weights of the multiplicative polynomial if used by PPXF.
;
;       bestfit_ppxf dblarr[C]
;               The best fit model determine by PPXF, defined over the
;               wavelength range of the object data (with C spectral channels).
;
;       chi2_ppxf double
;               Chi^2/DOF for the PPXF fit.
;
;       fitted_pixels_gndf intarr[]
;               Indices of the pixels used by GANDALF.
;
;       weights_gndf dblarr[T]
;               The weights applied to each template spectrum to create the
;               optimal composite template as determined by GANDALF.
;
;       mult_poly_coeff_gndf dblarr[P]
;               Legendre polynomial weights of the multiplicative polynomial if used by GANDALF.
;
;       bestfit_gndf dblarr[C]
;               The best fit model determined by GANDALF, defined over the
;               wavelength range of the object data (with C spectral channels).
;
;       chi2_gndf double
;               Chi^2/DOF for the GANDALF fit.
;
;       eml_model dblarr[C]
;               Best-fitting model of the emission lines.
;
;       sol dblarr[9]
;               The best-fit kinematic parameters:
;                       sol[0]: stellar velocity (km/s).
;                       sol[1]: stellar velocity dispersion  (km/s).
;                       sol[2]: stellar h3 Gauss-Hermite moment (km/s).
;                       sol[3]: stellar h4 Gauss-Hermite moment (km/s).
;                       sol[4]: stellar h5 Gauss-Hermite moment.
;                       sol[5]: stellar h6 Gauss-Hermite moment.
;                       sol[6]: chi^2 (NOT USED)
;                       sol[7]: mean flux weighted velocity of the emission lines (km/s).
;                       sol[8]: mean flux weighted velocity dispersion (km/s).  
;
;       err dblarr[8]
;               Error in the best-fit kinematic solution.  Elements are:
;                       err[0]: error on the stellar velocity (km/s).
;                       err[1]: error on the stellar velocity dispersion  (km/s).
;                       err[2]: error on the stellar h3 Gauss-Hermite moment.
;                       err[3]: error on the stellar h4 Gauss-Hermite moment.
;                       err[4]: error on the stellar h5 Gauss-Hermite moment (not used).
;                       err[5]: error on the stellar h6 Gauss-Hermite moment (not used).
;                       err[6]: error on the gas mean velocity (km/s).
;                       err[7]: error on the gas mean velocity dispersion (km/s). 
;
;               TODO: Use 9-element vector and ignore the chi-square (i=6)
;
;       gas_intens     dblarr[E]
;       gas_intens_err dblarr[E]
;               Intensity (corrected for reddening) of the E emission lines and
;               their errors.  The intensities of lines in a multiplet are
;               constrained by the flux ratio defined in the eml_par structure.
;
;       gas_vel     dblarr[E]
;       gas_vel_err dblarr[E]
;               Velocity of the E emission lines and their errors.  TODO: Are
;               the tied velocities in this?
;
;       gas_sig     dblarr[E]
;       gas_sig_err dblarr[E]
;               Velocity dispersions of the E emission lines and their errors.
;               TODO: Are the tied velocity dispersions in this?
;
;       gas_flux      dblarr[E]
;       gas_flux_err  dblarr[E]
;               Integrated fluxes (corrected for reddening) of the E emission
;               lines and their errors.  The intensities of lines in a multiplet
;               are constrained by the flux ratio defined in the eml_par
;               structure.
;
;       gas_EW     dblarr[E]
;       gas_EW_err dblarr[E]
;               Equivalent widths (corrected for reddening) of the E emission
;               lines and their errors.  The intensities of lines in a multiplet
;               are constrained by the flux ratio defined in the eml_par
;               structure.  Equivalent widths are computed by comparing the
;               emission-line flux with the median flux of the stellar continuum
;               in the spectral region
;                       l_i-10*FWHM_i < l < l_i-5*FWHM_i
;               and
;                       l_i+5*FWHM_i < l < l_i+10*FWHM_i,
;               where l_i is the central wavelength of the i-th emission line,
;               and FWHM_i is its measured FWHM (intrinsic plus instrumental).
;
; OPTIONAL OUTPUT:
;       reddening dblarr[1 or 2]
;               If it exists upon input, the stellar reddening (if it is
;               dblarr[1]) and the gas reddening (balmer decrement, if it is
;               dblarr[2]) are fit.  On output, it is replaced by the best-fit
;               reddening.
;
;       err_reddening dblarr[1 or 2]
;               Error in the reddening if fitted.
;
;       ppxf_status integer
;               If 0, pPXF failed; if 1, pPXF was successful.  If pPXF fails,
;               GANDALF is still executed, but the stellar kinematics ere fixed
;               to the starting guesses.
;
;       gandalf_status integer
;               If 0, GANDALF failed; if 1, GANDALF was successful.
;
; COMMENTS:
;
; EXAMPLES:
;
; TODO:
;       - Change to using structures!
;       - Add a log10 keyword
;
; BUGS:
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;       29 Sep 2014: Copied from v0_8 by L. Coccato (last edited 16 Apr 2014)
;       29 Sep 2014: (KBW) Formatting and edits
;-
;------------------------------------------------------------------------------

; "NEW - V1.3" comments denotes V1.3 modifications. Notably these are:
;
; 1) Keyword FOR_ERRORS have been added and is passed onto GANDALF
; 2) We now call MASK_EMISSION_LINES while specifying the desired
;    rest-frame wavelength to be fitted
; 3) MASK_EMISSION_LINES itself now checks for lines outside this
;    wavelength range, excluding them
; 4) We now keep excluding the Na D region, which can be affected by
;    interstellar absorption
; 5) Set the starting guess velocity for the gas not just to the
;    receding velocity of the stars but also allow for velocity
;    offset, as specfied in the emission-line setup. This makes it
;    easier to fit red or blue wings
; 6) We now append to the solution output the goodpixels array (1
;    good, 0, masked), which can be used later in making plots
;    with SHOW_FIT.PRO


;---------------------------------------------------------------------------------------------------
;       Description:
;
;PRO MDAP_GANDALFW_REMOVE_INDICES_FROM_GOODPIXELS, $
;               goodpixels, indx
;
;       ; TODO: Revisit what's going on here
;       n=n_elements(indx)
;       for i=0,n-1 do begin
;           rm = where(indx[i] eq goodpixels)
;           if rm[0] ne -1 then $
;               REMOVE, rm, goodpixels
;       endfor
;
;       kk = goodpixels[uniq(goodpixels,sort(goodpixels))]
;       goodpixels=kk
;end
;---------------------------------------------------------------------------------------------------




;---------------------------------------------------------------------------------------------------
;       Return a list of goodpixels to fit that excludes regions potentially
;       affected by gas emission and by sky lines.  Unless the log10 keyword is
;       specified, wavelength values are assumed to be ln-rebinned, and are
;       defined by the l0_gal, lstep_gal, npix parameters.  The position of gas
;       and sky emission lines is set by the input eml_par structure and the
;       width of the mask by the sigma parameter.  If a sigma value is not
;       passed than the width of each line is taken from the eml_par structure.
;
;       The rest-frame wavelength range can be manually restricted using the
;       l_rf_range keyword to pass min and max observed wavelengths.  Typically
;       used to exclude regions at either end of spectra.

;FUNCTION MDAP_GANDALFW_MASK_EMISSION_LINES, $
;               npix, Vsys, eml_par, velscale, l0_gal, lstep_gal, sigma=sigma, $
;               l_rf_range=l_rf_range, log10=log10
;
;       c = 299792.458d                         ; Speed of light in km/s
;
;       goodpixels = mdap_range(0,npix-1)               ; Initialize good-pixel array
;
;       ; if set, exclude regions at either ends of the spectra using the keyword l_rf_range
;       ; TODO: THIS WON'T WORK ANYMORE
;       if keyword_set(l_rf_range) then begin
;           pix0     = ceil((alog10(l_rf_range[0])-l0_gal)/lstep_gal+Vsys/velscale)
;           pix1     = ceil((alog10(l_rf_range[1])-l0_gal)/lstep_gal+Vsys/velscale)
;
;           goodpixels = mdap_range(max([pix0,0]),min([pix1,npix-1])) ; NEW - V1.3 
;       endif
;
;       tmppixels  = goodpixels
;
;       ; looping over the listed emission-lines and mask those tagged with an
;       ; 'm' for mask. Mask sky lines at rest-frame wavelength
;       for i = 0,n_elements(eml_par)-1 do begin
;           if (eml_par[i].action eq 'm') then begin
;               if (eml_par[i].name ne 'sky') then begin
;                   meml_cpix = ceil((alog10(eml_par[i].lambda)-l0_gal)/lstep_gal+Vsys/velscale)
;               endif else $
;                   meml_cpix = ceil((alog10(eml_par[i].lambda)-l0_gal)/lstep_gal) 
;       
;               ; set the width of the mask in pixels using either 3 times the
;               ; sigma of each line in the emission-line setup or the provided
;               ; sigma value 
;               if keyword_set(sigma) then begin
;                   msigma = 3*sigma/velscale
;               endif else $
;                   msigma = 3*eml_par[i].s/velscale
;
;               meml_bpix = meml_cpix - msigma
;               meml_rpix = meml_cpix + msigma
;               w = where(goodpixels ge meml_bpix and goodpixels le meml_rpix) 
;               if (w[0] ne -1) then begin
;                   tmppixels[w] = -1 
;               endif else $
;                   eml_par[i].action = 'i'             ; Ignore lines outside wavelength range
;           endif
;       endfor
;
;       return, goodpixels[where(tmppixels ne -1)]
;END

; Mimic the output from pPXF
; TODO: does not check that sol_in has the expected size!
PRO MDAP_MIMIC_PPXF_KIN, $
                sol_in, sol_out, chi2, moments=moments, mdegree=mdegree
        if n_elements(moments) ne 0 then begin
            sol_out = [sol_in[0:moments-1], chi2]
        endif else begin
            moments=2
            sol_out = [sol_in[0:1], chi2]
        endelse

        if n_elements(mdegree) eq 0 then $
            return

        sol_out=[sol_out[*], dblarr(mdegree)]
END

PRO MDAP_SAVE_MULT_LEGENDRE, $
                sol, mult

        max_moments=6
        nfit = n_elements(sol)-1
        mult = nfit gt max_moments ? sol[max_moments+1:n_elements(sol)-1] : -1
        
END

; Save the output kinematics to an array of a fixed size
; TODO: Does not check that the input from pPXF has the correct size
PRO MDAP_SAVE_STAR_KIN, $
                sol_in, sol_out, voff, moments=moments, chi=chi
        max_moments=6
        if n_elements(moments) eq 0 then $
            moments=2

        sol_out = [sol_in[0]-voff, sol_in[1:moments-1], dblarr(max_moments-moments)]
        if keyword_set(chi) then $
            sol_out = [sol_out, sol_in[max_moments]]

END

FUNCTION MDAP_GANDALFW_SOL, $
                sol_star, gas_vel, gas_sig

        return, ([ sol_star, gas_vel, gas_sig ])
END


;-------------------------------------------------------------------------------
; Create default or failed output for GANDALF
PRO MDAP_GANDALFW_DEFAULT, $
                sol_star, err_star, nc, ntpl, fitted_pixels_gndf, weights_gndf, $
                mult_poly_coeff_gndf, bestfit_gndf, chi2_gndf, eml_model, sol, err, $
                neml_f, gas_intens, gas_intens_err, gas_vel, gas_vel_err, gas_sig, gas_sig_err, $
                gas_flux, gas_flux_err, gas_ew, gas_ew_err, mdegree=mdegree, $
                reddening=reddening, err_reddening=err_reddening

        sol = MDAP_GANDALFW_SOL(sol_star, !VALUES.D_NAN, !VALUES.D_NAN)
        err = MDAP_GANDALFW_SOL(err_star, 99.0d, 99.0d) ; TODO: Set errors ton NaN?

        fitted_pixels_gndf=-1                           ; No pixels fitted by GANDALF
        weights_gndf = dblarr(ntpl)                             ; No template weights
        if n_elements(mdegree) ne 0 then begin
            mult_poly_coeff_gndf = dblarr(mdegree)              ; No polynomial weights
        endif else $
            MDAP_ERASEVAR, mult_poly_coeff_gndf         ; Remove variable
        bestfit_gndf = dblarr(nc)                               ; No GANDALF fit
        chi2_gndf = -1                                  ; No chi^2
        eml_model = dblarr(nc)                          ; No emission-line model

        gas_intens = neml_f gt 0 ? dblarr(neml_f) : -1
        gas_intens_err = neml_f gt 0 ? make_array(neml_f, /double, value=99.0d) : -1

        gas_vel = neml_f gt 0 ? make_array(neml_f, /double, value=!VALUES.D_NAN) : -1
        gas_vel_err = neml_f gt 0 ? make_array(neml_f, /double, value=99.0d) : -1

        gas_sig = neml_f gt 0 ? make_array(neml_f, /double, value=!VALUES.D_NAN) : -1
        gas_sig_err = neml_f gt 0 ? make_array(neml_f, /double, value=99.0d) : -1

        gas_flux = neml_f gt 0 ? dblarr(neml_f) : -1
        gas_flux_err = neml_f gt 0 ? make_array(neml_f, /double, value=99.0d) : -1
            
        gas_ew = neml_f gt 0 ? dblarr(neml_f) : -1
        gas_ew_err = neml_f gt 0 ? make_array(neml_f, /double, value=99.0d) : -1

        if n_elements(reddening) ne 0 then begin
            nred = n_elements(reddening)
            reddening = dblarr(nred)                    ; Set reddening to 0 with def error
            err_reddening = make_array(nred, /double, value=99.0d)
        endif
END
;---------------------------------------------------------------------------------------------------
;       Adjust the mask vector of the object spectrum to mask the emission lines
;       provided in the eml_par structure.  This masking is performed based on
;       the ACTION parameter in the structure.  The ACTION value can be either:
;
;               'i': ignore the line, as if the line were commented out.
;
;               'f': fit the line and mask the line when fitting the stellar
;               continuum.
;
;               'm': mask the line when fitting the stellar continuum but do NOT
;               fit the line itself
;
;               's': defines a sky line that should be masked.  When masked, the
;               wavelength of the line is NOT adjusted for the redshift of the
;               object spectrum.
;
;       If provided, masked lines will use the velocity value to redshift the
;       emission line to match the input frame.
;
;       TODO: The value of eml_par.v is NOT used to do this.
;
;       If provided, the width of the masking region is defined as 3 times sigma
;       (in km/s) on either side of the line center.  If not provided, the mask
;       size uses the eml_par.s value.
;       
;       TODO: Set a minimum mask size?
;
;       If mask_fit is defined, action='f' lines will also be masked.
;

PRO MDAP_GANDALFW_MASK_EMISSION_LINES, $
                eml_par, obj_wave, obj_mask, velocity=velocity, sigma=sigma, mask_fit=mask_fit

        if n_elements(eml_par) lt 1 then $      ; No emission lines; TODO: Is this check necessary?
            return

;       print, size(eml_par)
;       print, size(eml_par.action)
;       print, where(eml_par.action ne 'i')
;       print, n_elements(where(eml_par.action ne 'i'))
;       print, eml_par.action

        if keyword_set(mask_fit) then begin
            indx = where(eml_par.action ne 'i', count) ; Mask 'f', 'm', and 's' lines
        endif else $
            indx = where(eml_par.action ne 'i' and eml_par.action ne 'f', count) ; Mask 'm' and 's'

;       print, indx
;       stop
            
;       if indx[0] eq -1 then $
        if count eq 0 then $
            return                              ; No lines to mask

        nm = n_elements(indx)                   ; Number of lines to mask
        print, 'Number of lines to mask: ', nm
        
        ; TODO: Do NOT include sky lines in list of 'emission lines.' Create a
        ; separate structure that contains the observed wavelengths of sky lines
        ; that should be masked.

        for i=0,nm-1 do begin

            ; TODO: Number of sigma currently hard-wired to be 3.  Allow this to change?
            if eml_par[indx[i]].action eq 's' then begin
                el_range = MDAP_EMISSION_LINE_WINDOW(eml_par[indx[i]], sigma=sigma)
            endif else $
                el_range = MDAP_EMISSION_LINE_WINDOW(eml_par[indx[i]], velocity=velocity, $
                                                     sigma=sigma)

            MDAP_SELECT_WAVE, obj_wave, el_range, wave_indx, count=count    ; Get the indices
;           if wave_indx[0] ne -1 then $
            if count ne 0 then $
                obj_mask[wave_indx] = 1.0d                             ; Mask them
        endfor

END
;---------------------------------------------------------------------------------------------------


;---------------------------------------------------------------------------------------------------
;
;       Given the object spectrum, the best fit, the emission-line amplitudes, a
;       vector with the A/N threshold for each line, and the array containing
;       the spectra of each best-fitting emission-line template, this function
;       computes the residual-noise lvl, compares it the amplitude of each
;       lines, and removes from the galaxy spectrum only the best-matching
;       emission-line templates for which the correspoding A/N exceeds the input
;       threshold.  This is a necessary step prior to measurements of the
;       strength of the stellar absorption line features
;
;       A list of goodpixels may be provided.  This is useful if any pixel was
;       excluded by sigma-clipping during the continuum and emission-line
;       fitting or on either ends of the spectra.
;
;       Also optionally outputs the computed A/N ratios.

; TODO: Return to this!

FUNCTION MDAP_GANDALFW_EMISSION_LINE_NOISE, $
                galaxy, wave, noise, bestfit, eml_par, goodpixels=goodpixels, velocity=velocity, $
                sigma=sigma, resid=resid, noise_mask=noise_mask
;       print, 'found it!'

        i_f = where(eml_par.action eq 'f', count)                       ; Fitted emission lines
        if count eq 0 then $
            message, 'No emission lines to fit!'
        i_l = where(eml_par.action eq 'f' and eml_par.kind eq 'l', neml_l) ; Independent lines
        if neml_l eq 0 then $
            message, 'No fitted lines, all doublets!'
;       neml_l = n_elements(i_l)                                ; Number of fitted emission lines

        fit_noise = dblarr(neml_l, /nozero)                     ; Noise near the fitted line(s)

        if keyword_set(resid) then $
            resid_gal = galaxy - bestfit                        ; Get the residual of the fit

        if keyword_set(noise_mask) then begin
            for i=0,neml_l-1 do begin
                el_range = MDAP_EMISSION_LINE_WINDOW(eml_par[i_l[i]], velocity=velocity, $
                                                     sigma=sigma)
                MDAP_SELECT_WAVE, wave, el_range, indx, count=count     ; Get the indices
;               if indx[0] eq -1 then begin
                if count eq 0 then begin
                    fit_noise[i] = 99.0d
                    continue
                endif

                dindx = where(fix(strmid(eml_par[i_f].kind,1)) eq eml_par[i_l[i]].i, count)
;               if dindx[0] ne -1 then begin
                if count ne 0 then begin
                    ntied = n_elements(dindx)
                    for j=0,ntied-1 do begin
                        el_range = MDAP_EMISSION_LINE_WINDOW(eml_par[i_f[dindx[j]]], $
                                                             velocity=velocity, sigma=sigma)
                        MDAP_SELECT_WAVE, wave, el_range, aindx, count=acount    ; Get the indices
;                       if aindx[0] ne -1 then $
                        if acount ne 0 then $
                            indx = [indx, aindx]
                    endfor
                endif

                if n_elements(goodpixels) ne 0 then begin
                    indx = MDAP_SET_INTERSECTION(temporary(indx), goodpixels, count=count)
                    if count eq 0 then begin
                        fit_noise[i] = 99.0d
                        continue
                    endif 
                endif 
                fit_noise[i_f[i]] = (keyword_set(resid) ? robust_sigma(resid_gal[indx], /zero) : $
                                                          median(noise[indx]))
            endfor
        endif else $
            fit_noise[*] = (keyword_set(resid) ? robust_sigma(resid_gal[goodpixels], /zero) : $
                                                 median(noise[goodpixels]))

        return, fit_noise

END

FUNCTION MDAP_GANDALFW_REMOVE_DETECTED_EMISSION, $
                galaxy, bestfit, emission_templates, sol_gas_A, sol_gas_N, AoN_thresholds, $
                AoN=AoN, positive=positive

        neml = n_elements(sol_gas_A)            ; Number of emission lines
;       print, neml, n_elements(sol_gas_N)
        AoN = sol_gas_A/sol_gas_N               ; A/N of each line
;       print, n_elements(AoN)
;       print, n_elements(AoN_thresholds)

        indx = where(sol_gas_N lt 0, count)     ; Remove values with negative errors
;       if indx[0] ne -1 then $
        if count ne 0 then $
            AoN[indx] = 0.0d

        if keyword_set(positive) then begin     ; Remove emission lines with negative amplitudes
            indx = where(sol_gas_A lt 0, count)
;           if indx[0] ne -1 then $
            if count ne 0 then $
                AoN[indx] = 0.0d
        endif

        ; Remove the model emission lines from the galaxy spectrum
        neat_galaxy = galaxy
        for i = 0,neml-1 do begin
;           print, i, neml, AoN[i], AoN_thresholds[i]
            if (AoN[i] gt AoN_thresholds[i]) then $
                neat_galaxy = neat_galaxy - emission_templates[*,i]
                ; TODO: what is the order of the emission line templates
        endfor

        return, neat_galaxy
END
;---------------------------------------------------------------------------------------------------

; The number of output emission line fits is the same size as eml_par_fit.
; Doublets in eml_par_fit are given the same results as their linked lines with
; adjustments to their fluxes, etc, based on the provided line ratios.

; TODO: How are doublets treated in GANDALF?

; TODO: There was a bug in the v0_8 version of the code, that may or may not be
;       addressed by the codde below.  See line 739-740 of Brett's version, and
;       his e-mail from 5 Nov 2014.

PRO MDAP_GANDALFW_SAVE_RESULTS, $
                eml_par, sol_gas_A, esol_gas_A, sol_gas_V, esol_gas_V, sol_gas_S, esol_gas_S, $
                sol_gas_F, esol_gas_F, sol_gas_EW, esol_gas_EW, gas_intens, gas_intens_err, $
                gas_vel, gas_vel_err, gas_sig, gas_sig_err, gas_flux, gas_flux_err, gas_ew, $
                gas_ew_err

        i_f = where(eml_par.action eq 'f', neml_f)              ; Fitted emission lines
;       neml_f = n_elements(i_f)                                ; Number of fitted emission lines
        if neml_f eq 0 then $
            message, 'No emission lines to fit!'
        i_l = where(eml_par[i_f].kind eq 'l', neml_l)           ; Independent lines
        if neml_l eq 0 then $
            message, 'No fitted lines, all doublets!'

        ; Initialize the vectors
        gas_intens = dblarr(neml_f,/nozero)
        gas_intens_err = dblarr(neml_f,/nozero)
        gas_vel = dblarr(neml_f,/nozero)
        gas_vel_err = dblarr(neml_f,/nozero)
        gas_sig = dblarr(neml_f,/nozero)
        gas_sig_err = dblarr(neml_f,/nozero)
        gas_flux = dblarr(neml_f,/nozero)
        gas_flux_err = dblarr(neml_f,/nozero)
        gas_ew = dblarr(neml_f,/nozero)
        gas_ew_err = dblarr(neml_f,/nozero)

;       neml_l = n_elements(i_l)                                ; Number of independent lines
        for i=0,neml_l-1 do begin

            ; Copy the fitted lines
            gas_intens[i_l[i]] = sol_gas_A[i]
            gas_intens_err[i_l[i]] = esol_gas_A[i]
            gas_vel[i_l[i]] = sol_gas_V[i]
            gas_vel_err[i_l[i]] = esol_gas_V[i]
            gas_sig[i_l[i]] = sol_gas_S[i]
            gas_sig_err[i_l[i]] = esol_gas_S[i]
            gas_flux[i_l[i]] = sol_gas_F[i]
            gas_flux_err[i_l[i]] = esol_gas_F[i]
            gas_ew[i_l[i]] = sol_gas_EW[i]
            gas_ew_err[i_l[i]] = esol_gas_EW[i]

            ; Find any lines that are doublets of this line, and set the data
            ; after adjusting for the line ratio

            ; TODO: Set the error differently since these lines are not independently fit

            ; TODO: What about lines that are kinematically tied?

            ; TODO: Report the *measured* values, not the reddening corrected
            ;       values.  Need to figure out how to do this in the context of
            ;       two lines that are tied.
            indx = where(fix(strmid(eml_par[i_f].kind,1)) eq eml_par[i_f[i_l[i]]].i, count)
            if count ne 0 then begin
                gas_intens[indx] = eml_par[i_f[indx]].a*gas_intens[i_l[i]]
                gas_intens_err[indx] = eml_par[i_f[indx]].a*gas_intens_err[i_l[i]]
                gas_vel[indx] = gas_vel[i_l[i]]
                gas_vel_err[indx] = gas_vel_err[i_l[i]]
                gas_sig[indx] = gas_sig[i_l[i]]
                gas_sig_err[indx] = gas_sig_err[i_l[i]]
                gas_flux[indx] = eml_par[i_f[indx]].a*gas_flux[i_l[i]]
                gas_flux_err[indx] = eml_par[i_f[indx]].a*gas_flux_err[i_l[i]]
                gas_ew[indx] = eml_par[i_f[indx]].a*gas_ew[i_l[i]]
                gas_ew_err[indx] = eml_par[i_f[indx]].a*gas_ew_err[i_l[i]]
            endif
        endfor
END

PRO MDAP_GANDALF_WRAP, $
                tpl_wave, tpl_flux, obj_wave, obj_flux, obj_ivar, obj_mask, velScale, start, $
                fitted_pixels_ppxf, weights_ppxf, add_poly_coeff_ppxf, mult_poly_coeff_ppxf, $
                bestfit_ppxf, chi2_ppxf, fitted_pixels_gndf, weights_gndf, mult_poly_coeff_gndf, $
                bestfit_gndf, chi2_gndf, eml_model, sol, err, gas_intens, gas_intens_err, gas_vel, $
                gas_vel_err, gas_sig, gas_sig_err, gas_flux, gas_flux_err, gas_ew, gas_ew_err, $
                eml_par=eml_par, vsyst=vsyst, bias=bias, mdegree=mdegree, degree=degree, $
                moments=moments, reddening=ebv, range_v_star=range_v_star, $
                range_s_star=range_s_star, range_v_gas=range_v_gas, range_s_gas=range_s_gas, $
                region_mask=region_mask, external_library=external_library, int_disp=int_disp, $
                err_reddening=err_reddening, ppxf_status=ppxf_status, $
                gandalf_status=gandalf_status, ppxf_only=ppxf_only, oversample=oversample, $
                for_errors=for_errors, fix_star_kin=fix_star_kin, fix_gas_kin=fix_gas_kin, $
                quiet=quiet, plot=plot
                
        ; TODO: How degenerate are the reddening, err_reddening, and for_error checks?

        ; Preliminary Checks !
        if ~keyword_set(mdegree) then $
            mdegree=0

        ; TODO: Does this need to be done here.  Probably not because pPXF does the same thing.
        if n_elements(moments) eq 0 then $
            moments=2                                   ; Defaults to fit just V,sigma

        if n_elements(ppxf_only) eq 0 then $
            ppxf_only = 0

;       print, 'mdegree: ', mdegree
;       print, 'moments: ', moments

        sz = size(tpl_flux)
        ntpl = sz[1]                                    ; Number of templates

        ; Initialize the error vector and the mask
        ; TODO: Could I just continue using ivar?
        nc = n_elements(obj_flux)                               ; Number of spectral channels
        obj_mask_ = obj_mask                                    ; Copy the input mask; TODO: needed?

        MDAP_NOISE_FROM_IVAR, obj_ivar, obj_mask_, obj_sige

;       indx = where(obj_ivar eq 0 or finite(obj_ivar) eq 0, complement=nindx)  ; Invalid variances
;       if n_elements(indx) eq nc then begin                    ; No pixels are valid
;           print, 'All errors are invalid!  Ignoring errors (by setting them to unity).'
;           obj_sige = make_array(nc, /double, value=1.0d)
;       endif else if indx[0] eq -1 then begin                  ; All pixels are valid
;           obj_sige = sqrt(1.0d/obj_ivar)                      ; Errors
;       endif else begin
;           obj_mask_[indx] = 1.0d                              ; Mask invalid variances
;           ; TODO: Mask doesn't need to be double
;           obj_sige = dblarr(nc, /nozero)                      ; Initialize
;           obj_sige[nindx] = sqrt(1.0d/obj_ivar[nindx])        ; Errors
;           obj_sige[indx] = 1.0d                               ; Placeholder value (masked!)
;       endelse

        ; Mask the input wavelength ranges
        ; TODO: Is this still necessary?
        if n_elements(region_mask) ne 0 then begin
            ; TODO: MDAP_SPECTRAL_FITTING checks this also!
            if n_elements(region_mask) MOD 2 ne 0 then $
                message, 'Input region_mask does not have an even number of entries!'
            nm = n_elements(region_mask)/2
            for i=0,nm-1 do begin
                MDAP_SELECT_WAVE, obj_wave, region_mask[2*i:2*i+1], indx, count=count
                if count ne 0 then $
                    obj_mask_[indx] = 1.0d                      ; Add these regions to the mask
            endfor
        endif

        ; Set the starting guesses for pPXF
        start_ppxf = start[0:1]                                 ; Copy the first two values
        voff = (n_elements(vsyst) ne 0) ? vsyst : 0.0d          ; Offset velocity
        start_ppxf[0] = start[0]+voff                           ; Offset the guess velocity

;       print, "Start pPXF: ", start_ppxf
        
        ; Mask the emission lines:
        ;   Uses MDAP_GANDALFW_MASK_EMISSION_LINES to loop over the emission
        ;   lines and mask those with 'm', including sky lines.  The wavelengths
        ;   masked include the velocity shift, if they're not sky lines.  The
        ;   accuracy of the mask depends on a relatively good guess for the
        ;   velocity.  The width of the mask is set to +/-3.0*250 km/s.
        ; TODO: Allow this to change?

        ; TODO: Do NOT include sky lines in list of 'emission lines.' Create a
        ; separate structure that contains the observed wavelengths of sky
        ; lines that should be masked.

        ; Copy to the continuum mask (used in fitting the continuum, not the emission lines
        cnt_mask_ = obj_mask_
        neml = n_elements(eml_par)                      ; Number of defined emission lines
        if neml ne 0 then begin
            MDAP_GANDALFW_MASK_EMISSION_LINES, eml_par, obj_wave, cnt_mask_, velocity=start[0], $
                                               sigma=250.0d, /mask_fit

            ; Remove ignored ems lines from the int_disp vector
            ; TODO: int_disp_ no longer has the same number of elements as eml_par!
            indx = where(eml_par.action ne 'i' and eml_par.action ne 'sky', count)
            int_disp_ = (count ne 0) ? int_disp[indx] : int_disp
        endif else $
            int_disp_ = int_disp

;       print, 'Number of unmasked pixels, continuum fit: ', n_elements(where(cnt_mask_ lt 1.0))

        ; If the object flux is <= 0, then increase the noise in those regions.
        ; TODO: get rid of this
        indx = where(obj_flux le 0, count, complement=nindx, ncomplement=ncount)
;       if indx[0] ne -1 then begin
        if count ne 0 then begin
;           if nindx[0] ne -1 then begin
            if ncount ne 0 then begin
                maxerr = max(obj_sige[nindx])
            endif else $
                maxerr = max(obj_sige)
            obj_sige[indx] = maxerr
        endif

        ; Run pPXF, fitting only V and sigma.  A constant additive polynomial
        ; (degree=0) is used in together with the multiplicative polynomials
        ; (always recommended).
        ppxf_status = 1
        if ~keyword_set(fix_star_kin) then begin        ; TODO: How much of the above is necessary
                                                        ;       if fix_star_kin has been set?

            print, 'Starting pPXF'
            fitted_pixels_ppxf = where(cnt_mask_ lt 1., count)      ; Select the unmasked pixels
;           if fitted_pixels_ppxf[0] eq -1 then $
            if count eq 0 then $
                message, 'No pixels to fit by pPXF!'

;           plot, obj_wave, obj_flux, color=200
;           oplot, obj_wave[goodpixels], obj_flux[goodpixels]
;           stop

;           plot, obj_wave, obj_sige
;           stop

;           for i=0,ntpl-1 do begin
;               print, i+1, ntpl
;               plot, tpl_wave, tpl_flux[i,*]
;               stop
;           endfor

            ; TODO: Do checks of ranges AFTER running pPXF?

;           print, size(tpl_flux)
;           print, size(obj_flux)
;           print, size(obj_sige)
;           moments=0
;           mdegree=0
            ; TODO: Keep the additive and multiplicative continua using polyweights and sol
            ; TODO: Don't like passing transpose!
            MDAP_PPXF, transpose(tpl_flux), obj_flux, obj_sige, velScale, start_ppxf, sol_ppxf, $
                       bestfit=bestfit_ppxf, bias=bias, degree=degree, error=err_ppxf, $
                       goodpixels=fitted_pixels_ppxf, mdegree=mdegree, moments=moments, $
                       oversample=oversample, weights=weights_ppxf, range_v_star=range_v_star, $
                       range_s_star=range_s_star, external_library=external_library, $
                       polyweights=add_poly_coeff_ppxf, plot=plot
            ; TODO: Check the success of pPXF!
            ppxf_status = 0
        endif

        print, 'Done pPXF'

        ; If PPXF ended in error, or PPXF was not run then create default output
        if ppxf_status ne 0 then begin

            ; Mimic the pPXF using the starting guesses
            MDAP_MIMIC_PPXF_KIN, start, sol_ppxf, 99.0d, moments=moments, mdegree=mdegree
            err_ppxf = make_array(moments, /double, value=99.0d)        ; Set place-holder error

            weights_ppxf = dblarr(ntpl)                                 ; No template weights
            bestfit_ppxf = dblarr(nc)                                   ; No fit

            fitted_pixels_ppxf=-1                                       ; No pixels fit
            add_poly_coeff_ppxf=-1                                      ; No additive continuum

            ; TODO: If ppxf fails, shouldn't this just return?

        endif

        MDAP_SAVE_STAR_KIN, sol_ppxf, sol_star, voff, moments=moments, /chi
        MDAP_SAVE_STAR_KIN, err_ppxf, err_star, 0.0d, moments=moments

        ; Save the multiplicative polynomial terms
        MDAP_SAVE_MULT_LEGENDRE, sol_ppxf, mult_poly_coeff_ppxf

        chi2_ppxf = sol_star[n_elements(sol_star)-1]

;       print, sol_ppxf
;       print, err_ppxf
;
;       print, sol_star
;       print, err_star
;
;       print, weights_ppxf
;       stop

        ; No emission lines, or do not want to fit them, so return with just the result of pPXF
;       if n_elements(eml_par) eq 0 or ppxf_only eq 1 then begin
        if n_elements(eml_par) eq 0 || ppxf_only eq 1 then begin
            if n_elements(eml_par) ne 0 then begin
                i_f = where(eml_par.action eq 'f', neml_f)
;               neml_f = i_f[0] eq -1 ? 0 : n_elements(i_f)
            endif else $
                neml_f = 0

            MDAP_GANDALFW_DEFAULT, sol_star, err_star, nc, ntpl, fitted_pixels_gndf, weights_gndf, $
                                   mult_poly_coeff_gndf, bestfit_gndf, chi2_gndf, eml_model, sol, $
                                   err, neml_f, gas_intens, gas_intens_err, gas_vel, gas_vel_err, $
                                   gas_sig, gas_sig_err, gas_flux, gas_flux_err, gas_ew, $
                                   gas_ew_err, mdegree=mdegree, reddening=reddening, $
                                   err_reddening=err_reddening
            return
        endif

        ; TODO: Adjust the continuum mask?

        ; Now create a mask that includes the emission lines to fit
;       print, 'Number of unmasked pixels: ', n_elements(where(obj_mask_ lt 1.0))
        ems_mask_ = obj_mask_                           ; Emission line mask
        MDAP_GANDALFW_MASK_EMISSION_LINES, eml_par, obj_wave, ems_mask_, velocity=start[0], $
                                           sigma=250.0d

;       print, 'Number of unmasked pixels, emission-line fit: ', n_elements(where(ems_mask_ lt 1.0))

        ; Prepare emission_setup structure for GANDALF, which should only deal
        ;   with the lines we fit

        ; TODO: This is done in GANDALF as well.  Does it need to be done here?
        i_f = where(eml_par.action eq 'f', nf)
;       if i_f[0] eq -1 then $
        if nf eq 0 then $
            message, 'No lines to fit!'

;       nf = n_elements(i_f)
;       print, 'Number of lines to fit: ', nf

;       eml_par_fit = MDAP_ASSIGN_EMISSION_LINE_PARAMETERS(eml_par[i_f].i, eml_par[i_f].name, $
;                                                          eml_par[i_f].lambda, $
;                                                          eml_par[i_f].action, $
;                                                          eml_par[i_f].kind, eml_par[i_f].a, $
;                                                          eml_par[i_f].v, eml_par[i_f].s, $
;                                                          eml_par[i_f].fit)

        eml_par_fit = create_struct('i',eml_par[i_f].i,'name',eml_par[i_f].name,$
                               'lambda',eml_par[i_f].lambda,'action',eml_par[i_f].action,$
                               'kind',eml_par[i_f].kind,'a',eml_par[i_f].a,$
                               'v',eml_par[i_f].v,'s',eml_par[i_f].s,$
                               'fit',eml_par[i_f].fit)


        ; Assign input initial guess for the gas kinematics, adding also the
        ; velocity offset specified in the emission-line setup.
;       eml_par_fit[*].v = start[4]                     ; Gas velocity
;       eml_par_fit[*].s = start[5]                     ; Gas velocity dispersion
        eml_par_fit.v[*] = start[4]                     ; Gas velocity
        eml_par_fit.s[*] = start[5]                     ; Gas velocity dispersion

; mdegree_=mdegree
; if n_elements(reddening) ne 0 then junk = temporary(mdegree_)

        ; Include a multiplicative polynomial ONLY if the reddening is not calculated
        ; TODO: Throw an error instead
        if n_elements(reddening) eq 0 then $
            mdegree_=mdegree
        ; Otherwise, mdegree_ is undefined

        ; Call GANDALF using only the stellar kinematics as input sol
        fitted_pixels_gndf = where(ems_mask_ lt 1., count)      ; Select the unmasked pixels
;       if fitted_pixels_gndf[0] eq -1 then $
        if count eq 0 then $
            message, 'No pixels to fit by GANDALF!'

        ; TODO: Need to check that my eml_par works the same as their emission_setup
        l0_gal = alog10(obj_wave[0])
        lstep_gal = alog10(obj_wave[1])-alog10(obj_wave[0])
        l0_templ = alog10(tpl_wave[0])

;       print, l0_gal, lstep_gal, velScale, velScale/c/alog(10.0d), l0_templ

        print, 'Starting GANDALF'

;       stop
        sol_gas = sol_ppxf
        MDAP_GANDALF, transpose(tpl_flux), obj_flux, obj_sige, velScale, sol_gas, eml_par_fit, $
                      l0_gal, lstep_gal, goodpixels=fitted_pixels_gndf, int_disp=int_disp_, $
                      bestfit=bestfit_gndf, emission_templates=emission_templates, $
                      weights=weights_gndf, l0_templ=l0_templ, degree=-1, mdegree=mdegree_, $
                      /for_errors, error=esol_gas, reddening=reddening, fix_gas_kin=fix_gas_kin, $
                      range_v_gas=range_v_gas, range_s_gas=range_s_gas, $
                      external_library=external_library, plot=plot, /log10, $
                      mult_poly_coeff=mult_poly_coeff_gndf;, quiet=quiet

        print, 'DONE GANDALF'

        ; TODO: Incude status
        gandalf_status = 0

        ; TODO: Output three structures:
        ;       PolynomialFit : type[add, mult], procedure[ppxf, gandalf], order, coeffs
        ;       StellarContinuumFit : template weights, sol[v,sig,h3,h4,...], chi2
        ;       EmissionLineFit: intensity, sol[v,sig,h3,h4,...], flux, EW, chi2
        ;               Save individual fits and mean values

        ; TODO: Save continuum used in EW measurement? (Can get it from EW and flux)

        if gandalf_status ne 0 then begin               ; GANDALF failed!
            neml_f = n_elements(i_f)                    ; Number of *fitted* emission lines
            MDAP_GANDALFW_DEFAULT, sol_star, err_star, nc, ntpl, fitted_pixels_gndf, weights_gndf, $
                                   mult_poly_coeff_gndf, bestfit_gndf, chi2_gndf, eml_model, sol, $
                                   err, neml_f, gas_intens, gas_intens_err, gas_vel, gas_vel_err, $
                                   gas_sig, gas_sig_err, gas_flux, gas_flux_err, gas_ew, $
                                   gas_ew_err, mdegree=mdegree, reddening=reddening, $
                                   err_reddening=err_reddening
            return                                      ; Done
        endif

;       print, size(sol_ppxf)
;       print, size(sol_star)
;       print, size(sol_gas)

        ; Select the independently fitted lines
        i_l = where(eml_par.action eq 'f' and eml_par.kind eq 'l', count);, complement=i_d) 
        ; TODO: Do this check earlier!
;       if i_l[0] eq -1 then $
        if count eq 0 then $
            message, 'All lines are doublets!'
        neml_l = n_elements(i_l)

        ; Copy over the data to new variables
        sol_gas_F = sol_gas[indgen(neml_l)*4+0] ; Fluxes
        sol_gas_A = sol_gas[indgen(neml_l)*4+1] ; Amplitudes
        sol_gas_V = sol_gas[indgen(neml_l)*4+2] ; Velocities
        sol_gas_S = sol_gas[indgen(neml_l)*4+3] ; Sigmas
        ; TODO: Higher moments in the gas?

;       print, sol_gas_F
;       stop

        ; Save the errors as well
        if keyword_set(for_errors) then begin
            esol_gas_F = esol_gas[indgen(neml_l)*4+0]
            esol_gas_A = esol_gas[indgen(neml_l)*4+1]
            esol_gas_V = esol_gas[indgen(neml_l)*4+2]
            esol_gas_S = esol_gas[indgen(neml_l)*4+3]
        endif

        ; TODO: Just save directly to reddening?
        if n_elements(reddening) ne 0 then begin        ; Reddening
            sol_EBmV = sol_gas[neml_l*4:*]
        endif else $
            sol_EBmV =-999.0d                           ; Default value; TODO: Why not -99?
        if keyword_set(for_errors) then begin
            if n_elements(reddening) ne 0 then begin
                esol_EBmV  = esol_gas[n_elements(i_l)*4:*]
            endif else $
                esol_EBmV  = -999.0d                    ; Default value; TODO: Why not -99?
        endif

        ; Get the noise in the amplitude either using the residual or the error
        ; TODO: Currently hard-wired to use the residual
        sol_gas_N = MDAP_GANDALFW_EMISSION_LINE_NOISE(obj_flux, obj_wave, obj_sige, bestfit_gndf, $
                                                      eml_par, goodpixels=fitted_pixels_gndf, $
                                                      /resid)
        AoN_thresholds = dblarr(neml_l)                 ; Set threshold=0, i.e. subtract all lines
;       print, neml_l, n_elements(sol_gas_F), n_elements(sol_gas_N)
;       print, size(emission_templates)
;       nemlt = (size(emission_templates))[2]
;       for i=0,nemlt-1 do begin
;           plot, obj_wave, emission_templates[*,i]
;           stop
;       endfor

        ; Remove the model emission lines from the galaxy spectrum
        spec_neat = MDAP_GANDALFW_REMOVE_DETECTED_EMISSION(obj_flux, bestfit_gndf, $
                                                           emission_templates, sol_gas_A, $
                                                           sol_gas_N, AoN_thresholds, $
                                                           AoN=sol_gas_AoN, /positive)

;       plot, obj_wave, spec_neat
;       stop

        ; TODO: Difference between this and the sum of all the (valid) emission-line templates?
        eml_model = obj_flux-spec_neat          ; model of the emission lines

;       plot, obj_wave, eml_model
;       stop

        ; TODO: Set criteria for whether or not a line is included in the mean kinematics
        indx = MDAP_LINES_FOR_MEAN_GAS_KINEMATICS(sol_gas_A, esol_gas_A, sol_gas_F, esol_gas_F, $
                                                  sol_gas_V, esol_gas_V, sol_gas_S, esol_gas_S, $
                                                  count=count)

        ; Set the mean kinematics over the valid lines
;       if indx[0] eq -1 then begin
        if count eq 0 then begin
            print, 'No valid lines for use in mean kinematics'
;           fw_vel = !VALUES.D_NAN
;           fw_sig = !VALUES.D_NAN
            fw_vel = start[4]
            fw_sig = start[5]

            fw_vel_err = 99.0d
            fw_sig_err = 99.0d

        endif else begin
            MDAP_MEAN_GAS_KINEMATICS, sol_gas_F[indx], sol_gas_V[indx], esol_gas_V[indx], $
                                      sol_gas_S[indx], esol_gas_S[indx], fw_vel, fw_vel_err, $
                                      fw_sig, fw_sig_err
        endelse

;       print, fw_vel, fw_sig
;       stop

        ; Weights are saved via input/output to GANDALF

        ; Save the kinematics
        sol = MDAP_GANDALFW_SOL(sol_star, fw_vel, fw_sig)               ; Final sol vector
        err = MDAP_GANDALFW_SOL(err_star, fw_vel_err, fw_sig_err)       ; ... and error

        ; Recalculate the chi-square including the emission-line fits
        ; TODO: Have GANDALF return this!
        ; TODO: goodpixels should always be ne 0, otherwise nothing is fit!
        ; TODO: Keep the pPXF-only chi^2
;       chi2=total(((obj_flux[goodpixels]-bestfit[goodpixels])/obj_sige[goodpixels])^2)
        resid = (obj_flux[fitted_pixels_gndf]-bestfit_gndf[fitted_pixels_gndf]) / $
                obj_sige[fitted_pixels_gndf]
        chi2_gndf = robust_sigma(resid, /zero)^2

;       print, sol
;       print, err

        ; Measure the equivalent widths
        MDAP_MEASURE_EMISSION_LINE_EQUIV_WIDTH, eml_par[i_l], fw_vel, fw_sig, obj_wave, spec_neat, $
                                                sol_gas_F, esol_gas_F, sol_gas_EW, esol_gas_EW

        ; Save the full set of results
        MDAP_GANDALFW_SAVE_RESULTS, eml_par, sol_gas_A, esol_gas_A, sol_gas_V, esol_gas_V, $
                                    sol_gas_S, esol_gas_S, sol_gas_F, esol_gas_F, sol_gas_EW, $
                                    esol_gas_EW, gas_intens, gas_intens_err, gas_vel, gas_vel_err, $
                                    gas_sig, gas_sig_err, gas_flux, gas_flux_err, gas_ew, gas_ew_err

        ; Save the reddening data
        if n_elements(reddening) ne 0 then begin
            reddening = sol_EBmV
            err_reddening = esol_EBmV
        endif

END

