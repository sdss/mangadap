;+
; NAME:
;       MDAP_SPECTRAL_FITTING
;
; PURPOSE:
;       TODO: Fill this in
;
;               Warning: The reddening (stars and/or gas) fit is optional, and
;               it is performed by the Gandalf module. If the reddening fit is
;               required, MDEGREE and DEGREE are used as the input value for the
;               pPXF run, but automatically set to 0 and -1 respectively in the
;               Gandalf execution.
;
;               WARNING: the name of the emission line in the final DAP output
;               will be defined by the string: CODE +''_''+ROUND(wav).
; 
; CALLING SEQUENCE:
;       MDAP_SPECTRAL_FITTING, obj_wave, obj_flux, obj_ivar, obj_mask, obj_sres, tpl_wave, $
;                              tpl_flux, tpl_ivar, tpl_mask, wavelength_output, obj_fit_mask_ppxf, $
;                              weights_ppxf, add_poly_coeff_ppxf, mult_poly_coeff_ppxf, $
;                              bestfit_ppxf, chi2_ppxf, obj_fit_mask_gndf, weights_gndf, $
;                              mult_poly_coeff_gndf, bestfit_gndf, chi2_gndf, eml_model, $
;                              stellar_kinematics, stellar_kinematics_err, $
;                              emission_line_kinematics, emission_line_kinematics_err, $
;                              emission_line_omitted, emission_line_kinematics_individual, $
;                              emission_line_kinematics_individual_err, emission_line_intens, $
;                              emission_line_intens_err, emission_line_fluxes, $
;                              emission_line_fluxes_err, emission_line_EW, emission_line_EW_err, $
;                              reddening_output, reddening_output_err, analysis_par=analysis_par, $
;                              default_velocity=default_velocity, $
;                              default_velocity_dispersion=default_velocity_dispersion, $
;                              star_kin_starting_guesses=star_kin_starting_guesses, $
;                              gas_kin_starting_guesses=gas_kin_starting_guesses, eml_par=eml_par, $
;                              range_v_star=range_v_star, range_s_star=range_s_star, $
;                              range_v_gas=range_v_gas, range_s_gas=range_s_gas, $
;                              wavelength_input=wavelength_input, $
;                              external_library=external_library, region_mask=region_mask, $
;                              wave_range_analysis=wave_range_analysis, ppxf_only=ppxf_only, $
;                              version=version, /use_previous_guesses, /fix_star_kin, $
;                              /fix_gas_kin, /quiet, /rest_frame_log, /oversample, /plot, /dbg
;
; INPUTS:
;       obj_wave dblarr[C]
;               Wavelength of all C spectral channels, which is the same for all
;               N spectra.
;
;       obj_flux dblarr[N][C]
;               Object flux for each of N spectra with C spectral channels.
;
;       obj_ivar dblarr[N][C]
;               Inverse variance in object flux for each of N spectra with C
;               spectral channels.
;
;       obj_mask dblarr[N][C]
;               Pixel mask for each of N spectra with C spectral channels.
;
;       obj_sres dblarr[C]
;               Spectral resolution at each of C spectral channels, which is the
;               same for all N spectra.
;
;       tpl_wave dblarr[S]
;               Wavelength of all S spectral channels, which is the same for all
;               template spectra.
;
;       tpl_flux dblarr[T][S]
;               Template flux for the library of T spectra with S spectral
;               channels.
;
;       tpl_ivar dblarr[T][S]
;               Inverse variance in the template flux for each of the T spectra
;               with S spectral channels.
;
;       tpl_mask dblarr[T][S]
;               Pixel mask for each of T spectra with S spectral channels.
;
;       star_kin_starting_guesses dblarr[N][2]
;               The stellar kinematics starting guesses for V, sigma, H3, and H4
;               for the N galaxy spectra to fit.  If not provided, the default
;               guess is V=H3=H4=0 and sigma=50 km/s.  Starting guess values
;               are overridden by the \use_previous_guesses keyword, if set.
;
;       gas_kin_starting_guesses dblarr[N][2]
;               The emission-line kinematics starting guesses for V and sigma
;               for the N galaxy spectra to fit.  If not provided, the default
;               guess is V=0 km/s and sigma = 50 km/s.  Starting guess values
;               are overridden by the \use_previous_guesses keyword, if set.
;
; OPTIONAL INPUTS:
;       analysis_par AnalysisPar structure
;               A structure that defines parameters used by PPXF and GANDALF in
;               the fitting procedure.  See its definition in
;               MDAP_DEFINE_ANALYSIS_PAR.pro
;
;       default_velocity double
;               Default velocity to use if *_kin_starting_guesses are
;               not supplied or they have the wrong size.
; 
;       default_velocity_dispersion double
;               Default velocity dispersion to use if
;               *_kin_starting_guesses are not supplied or they have the
;               wrong size.
;
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
;       range_v_star dblarr[2]
;               Allowed range for the stellar velocity in km/s.  The default
;               is set to the starting guess +/- 2000 km/s.
;
;       range_s_star dblarr[2]
;               Allowed range for the stellar velocity dispersion in km/s.  The
;               default is set 21 < sigma < 499 km/s.
;
;       range_v_gas dblarr[2]
;               Allowed range for the gas velocity in km/s.  The default is set
;               to the starting guess +/- 2000 km/sec.
;
;       range_s_gas dblarr[2]
;               Allowed range for the gas velocity dispersion in km/s.  The
;               default is set to the starting guess +/- 2000 km/s.
;
;               TODO: Is the above correct?  What about H3 and H4 for the stars?
;
;       wavelength_input dblarr[QQ]
;               If specified, it will be used to create wavelength_output, i.e.
;               the wavelength vector (constant ang/pixel step, in linear units)
;               to interpolate the final results on.  If keyword /rest_frame_log
;               is set, the vector is set to exp(loglam_templates), and user
;               input will be overwritten. In this case QQ = MM. The default is
;               to use the wavelength vector specified by the template stars.
;
;       external_library string
;               Path to the external FORTRAN library, which contains the fortran
;               versions of mdap_bvls.pro.  If not specified, or if the path is
;               invalid, the default internal IDL mdap_bvls code is used. 
;
;       region_mask dblarr[2*N]
;               If defined, it provides specifies a set of wavelength ranges
;               (lower and upper limits in angstroms) to omit from the fit.  It
;               must contain an even number of entries such that the ranges are:
;                   lower_wave_limit[i] = region_mask[2*i]
;                   upper_wave_limit[i] = region_mask[2*i+1]
;               for i=0..N-1.
;
;       wave_range_analysis dblarr[2]
;               Works the same as region_mask, but globally to all
;               spectra.  This is used here, region_mask is used by
;               MDAP_GANDALF_WRAP.
;
;       ppxf_only integer
;               Flag to run the GANDALF wrapper only using PPXF to fit the stellar
;               kinematics and optimal template.  Emission lines will be masked
;               if eml_par is provided; however, the emissions will not be fit
;               using GANDALF.  0/1 = false/true.
;
; OPTIONAL KEYWORDS:
;       /use_previous_guesses
;               Set the starting guesses for the i-th spectrum based on the
;               result obtained for the (i-1)-th spectrum.  **ANY input starting
;               guesses will be ignored (TODO: except for the fit to the first
;               spectrum?).**
;
;       /fix_star_kin
;               Fix the stellar kinematics during the fit (TODO: To the 'guess'
;               values?).  The returned values are the same as the starting
;               guess.
;
;       /fix_gas_kin
;               Fix the gas kinematics are during the fit (TODO: To the 'guess'
;               values?).  The returned values are the same as the starting
;               guess.
;
;       /quiet
;               Suppress information printed to the screen.
;
;       /rest_frame_log
;               Shift all output spectra (galaxy_minus_ems_fit_model,
;               best_fit_model, residuals, best_template, and
;               best_template_losvd_conv) to the rest frame.
;
;       /oversample
;               Turn on oversampling in pPXF.
;
;       /plot
;               Produce the plots generated by PPXF and GANDALF
;
;       /dbg
;               Only fit the first spectrum then return, used as a test
;               during debugging.
;
; OUTPUT:
;       wavelength_output dblarr[QQ]
;               The linear wavelength coordinates over which the output spectra
;               are sampled.  Default behavior is to use wavelength_input (if
;               defined) or automatically compute the vector using the smallest
;               lambda/pixel step obtained from exp(loglam_gal).  If
;               /rest_frame_log is set, then wavelength_output is set to
;               exp(loglam_templates).  TODO: Why is this needed?
;
;       obj_fit_mask_ppxf dblarr[N][QQ]
;               Bad pixel mask for pixels fitted by PPXF.  Pixels included/not
;               included in the fit are given values of 0.0/1.0.
;
;       weights_ppxf dblarr[N][T]
;               Template weights for each spectrum obtained by PPXF.
;
;       add_poly_coeff_ppxf dblarr[N][AD]
;               Additive order AD legendre polynomial coefficients obtained for
;               each of N spectra by PPXF.
;
;       mult_poly_coeff_ppxf dblarr[N][MD]
;               Multiplicative order MD legendre polynomial coefficients obtained
;               for each of N spectra by PPXF.
;
;       bestfit_ppxf dblarr[N][QQ]
;               Best fitting spectrum obtained by PPXF for each of the N
;               spectra.
;
;       chi2_ppxf dblarr[N]
;               Chi-square per degree of freedom obtained by the PPXF fit to
;               each of the N spectra.
;
;       obj_fit_mask_gndf dblarr[N][QQ]
;               Bad pixel mask for pixels fitted by GANDALF.  Pixels
;               included/not included in the fit are given values of 0.0/1.0.
;
;       weights_gndf dblarr[N][T]
;               Template weights for each spectrum obtained by GANDALF.
;
;       mult_poly_coeff_gndf dblarr[N][MD]
;               Multiplicative order MD legendre polynomial coefficients obtained
;               for each of N spectra by GANDALF.
;
;       bestfit_gndf dblarr[N][QQ]
;               Best fitting spectrum obtained by GANDALF for each of the N
;               spectra.
;
;       chi2_gndf dblarr[N]
;               Chi-square per degree of freedom obtained by the GANDALF fit to
;               each of the N spectra.
;
;       eml_model dblarr[N][QQ]
;               Best-fitting emission-line-only model for each of the N spectra
;               obtained by GANDALF.
;
;       stellar_kinematics dblarr[N][M]
;               The best-fit M moments for each of the N fitted input
;               galaxy spectra.  If \fix_star_kin is set,
;               MDAP_GANDALF_WRAP does NOT run the pPXF step.  This
;               array is then just a copy of the kinematic guesses.
;
;       stellar_kinematics_err dblarr[N][M]
;               Formal estimates of the errors in the best-fit M moments
;               for the stellar kinematics for each of the N fitted
;               input galaxy spectra.
;
;       emission_line_kinematics dblarr[N][2]
;               The best-fit V and sigma for the emission lines in the N
;               galaxy spectra. If \fix_gas_kin is set, output values
;               are the same as the input
;
;       emission_line_kinematics_err dblarr[N][2]
;               Formal estimates of the errors in the best-fit V and
;               sigma for the gas kinematics for each of the N fitted
;               input galaxy spectra.
;
;       emission_line_omitted intarr[N][E]
;               Flag setting whether or not an emission-line was fit for all E
;               emission lines in the eml_par structure.  0 means the
;               emission-line was fit, 1 means it was not.  TODO: Alter this
;               based on the success of the fit?
;
;       emission_line_kinematics_individual dblarr[N][E][2]
;       emission_line_kinematics_individual_err dblarr[N][E][2]
;               Kinematics and errors for each fitted emission line. TODO:
;               Better description?
;
;       emission_line_intens dblarr[N][E]
;       emission_line_intens_err dblarr[N][E]
;               Best-fitting emission-line intensities and errors of all fitted
;               emission lines for each of the  the N galaxy spectra.
;
;       emission_line_fluxes  dblarr[N][E]
;               Reddening-corrected integrated fluxes (TODO: based on a Gaussian
;               line?) of all fitted emission lines for each of the N input
;               galaxy spectra. 
;
;       emission_line_fluxes_err dblarr[N][E]
;               Errors in all fitted emission line fluxes.  TODO: How are these calculated?
;
;       emission_line_EW dblarr[N][E]
;               Equivalent widths of all fitted emission lines for each of the N
;               input galaxy spectra.  These equivalent widths are computed by
;               taking the ratio of the emission_line_fluxes and the median
;               value of the stellar spectrum within 5 and 10 sigma of the
;               emission line, where 'sigma' is the velocity dispersion (TODO:
;               Accounting for instrumental broadening? Or the observed value?)
;               emission line velocity dispersion.
; 
;       emission_line_EW_err dblarr[N][E]
;               Calculated errors in the equivalent width measurements.  TODO:
;               How is this done?
;
;       reddening_output dblarr[N][2]
;               Best-fit values for stellar reddening (reddening_output[*,0])
;               and gas reddening (reddening_output[*,1]) for all N galaxy
;               spectra.  If the reddening fit is not performed, the output
;               value is set to 0. (reddening_output[*,0:1] = [0,0]).  If only
;               the reddening of the stars is fitted, the reddening of the gas
;               is set to 0 (reddening_output[*,1] = 0).
;
;
;       reddening_output_err dblarr[N][2]
;               Errors associated to reddening_ output.  If not fitted, the
;               error is automatically set to 99. TODO: How are the errors
;               determined?
; 
; OPTIONAL OUTPUT:
;       version string
;               Module version.  If requested, the module is not executed and
;               only version flag is returned.
;
; COMMENTS:
;
; EXAMPLES:
;
; TODO:
;       - Put more of the output into structures!
;       - Include a log10 keyword
;
; BUGS:
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;       24 Sep 2014: Copied from v0_8 by L. Coccato
;       24 Sep 2014: Formatting + many edits by K. Westfall (KBW)
;       16 Oct 2014: (KBW) Rearranged input and added more output
;       11 Nov 2014: (KBW) Do not produce or output the optimal templates.
;                          These were originally used as input for the spectral
;                          index measurements; however, their resolutions had to
;                          be matched to the index system.  To make things more
;                          efficient, I now create a resolution matched set of
;                          templates at the start.  The spectral index
;                          measurements are then based on the optimal templates
;                          created by combining these templates with the weights
;                          determined here.
;       22 Jan 2015: (KBW) Better handling if default guesses are used.
;                          Added default_velocity keyword.
;       22 Mar 2015: (KBW) Ensure stellar kinematic moments allow for
;                          something other than moments=4.  Error in
;                          reporting of additive polynomial.  degree=0
;                          is valid as a constant offset.  degree=-1 is
;                          for when no additive polynomials are allowed.
;                          Other minor edits.
;       10 Jul 2015: (KBW) Converted placeholder values to -9999.0d
;       10 Aug 2015: (KBW) Moved selection of pixels to fit to
;                          pro/utilities/mdap_spectral_fitting_mask.pro
;                          to allow for that function to be used
;                          elsewhere.
;       15 Nov 2015: (KBW) Make sure that the starting velocity
;                          dispersion of the stars is larger than the
;                          minimum allowed by pPXF
;       08 Mar 2015: (KBW) Add oversample as input parameter
;-
;------------------------------------------------------------------------------

PRO MDAP_SPECTRAL_FITTING_CHECK_DIMENSIONS, $
                flux, ivar, mask, wave, sres

        ; The only things that have to be provided are the flux and the wavelength
;       if n_elements(flux) eq 0 or n_elements(wave) eq 0 then $
        if n_elements(flux) eq 0 || n_elements(wave) eq 0 then $
            message, 'Flux and/or wavelength arrays not defined!'

        sz=size(flux)                           ; All array must match the size of flux

        sz_ = size(wave)
        if sz_[0] ne 1 then $
            message, 'wave must be one dimensional!'
;       if (sz_[0] eq 1 and sz[2] ne sz_[1]) or (sz_[0] eq 2 and sz[2] ne sz_[2]) then $
        if (sz_[0] eq 1 && sz[2] ne sz_[1]) || (sz_[0] eq 2 && sz[2] ne sz_[2]) then $
            message, 'flux and wave must have the same number of spectral channels.'

        if n_elements(sres) ne 0 then begin 
            sz_ = size(sres)
            if sz_[0] ne 1 then $
                message, 'sres must be one dimensional!'
;           if (sz_[0] eq 1 and sz[2] ne sz_[1]) or (sz_[0] eq 2 and sz[2] ne sz_[2]) then $
            if (sz_[0] eq 1 && sz[2] ne sz_[1]) || (sz_[0] eq 2 && sz[2] ne sz_[2]) then $
                message, 'flux and sres must have the same number of spectral channels.'
        endif

        if n_elements(ivar) eq 0 then begin
            sz_ = size(ivar)
;           if sz[0] ne sz_[0] or sz[1] ne sz_[1] or (sz_[0] eq 2 and sz[2] ne sz_[2]) then $
            if sz[0] ne sz_[0] || sz[1] ne sz_[1] || (sz_[0] eq 2 && sz[2] ne sz_[2]) then $
                message, 'flux and ivar sizes must match!'
        endif

        if n_elements(mask) eq 0 then begin
            sz_ = size(mask)
;           if sz[0] ne sz_[0] or sz[1] ne sz_[1] or (sz_[0] eq 2 and sz[2] ne sz_[2]) then $
            if sz[0] ne sz_[0] || sz[1] ne sz_[1] || (sz_[0] eq 2 && sz[2] ne sz_[2]) then $
                message, 'flux and mask sizes must match!'
        endif

        ; Ensure that the spectra, and associated arrays, are two-dimensional
        if sz[0] eq 1 then begin
            flux=reform(flux, 1, sz[1])                 ; Arrays with one row
            ivar=reform(ivar, 1, sz[1])
            mask=reform(mask, 1, sz[1])
        endif

;       print, 'Dimensions are fine.'
END

; WARNING: Default velocity is 0!
; Will flag if use_defaults is used (use_defaults=1)
PRO MDAP_SPECTRAL_FITTING_INIT_GUESSES, $
                guesses, ns, default_velocity, default_velocity_dispersion

        use_defaults = 1
        if n_elements(guesses) ne 0 then begin
            sz=size(guesses)
            if sz[0] ne 2 then begin
                print, 'Kinematic guess arrays must be two-dimensional!  Using default.'
            endif else $
                use_defaults=0
        endif

;       print, use_defaults
;       print, size(guesses)
;       stop

        ; Use the provided default values
        if use_defaults eq 1 then begin
            if n_elements(default_velocity) eq 0 then $
                message, 'Must provided default velocity or guess array!'
            if n_elements(default_velocity_dispersion) eq 0 then $
                message, 'Must provided default velocity disperison or guess array!'
            guesses = dblarr(ns, 2)
            guesses[*,0] = default_velocity
            guesses[*,1] = default_velocity_dispersion
        endif
END

; TODO: Remove this.  Just always output on input wavelength grid
PRO MDAP_SPECTRAL_FITTING_SET_OUTPUT_WAVELENGTHS, $
                obj_wave, tpl_wave, wavelength_output, same_wave_as_input, $
                wavelength_input=wavelength_input, rest_frame_log=rest_frame_log

        if keyword_set(rest_frame_log) then begin
            wavelength_output=tpl_wave                  ; Use template wavelengths
            same_wave_as_input = 0
            print, 'Output wavelengths set to template wavelengths.'
            return
        endif

        if n_elements(wavelength_input) ne 0 then begin
            wavelength_output=wavelength_input          ; Use input
            same_wave_as_input = 0
            print, 'Using wavelength_input for output wavelengths'
            return
        endif

        wavelength_output=obj_wave                      ; Use object wavelengths
        same_wave_as_input = 1
;       print, 'Input and output wavelengths are the same.'

;       nw = n_elements(obj_wave)
;       dw = obj_wave[1]-obj_wave[0]
;       n= ceil(obj_wave[nw-1] - obj_wave[0]) / dw)
;       wavelength_output = dindgen(n)*dw + obj_wave[0]
END

FUNCTION MDAP_SPECTRAL_FITTING_START_KIN, $
                star_kin_starting_guesses, gas_kin_starting_guesses
        return, [ star_kin_starting_guesses, gas_kin_starting_guesses ]
END

PRO MDAP_SPECTRAL_FITTING, $
                obj_wave, obj_flux, obj_ivar, obj_mask, obj_sres, tpl_wave, tpl_flux, tpl_ivar, $
                tpl_mask, star_kin_starting_guesses, gas_kin_starting_guesses, wavelength_output, $
                obj_fit_mask_ppxf, weights_ppxf, add_poly_coeff_ppxf, $
                mult_poly_coeff_ppxf, bestfit_ppxf, chi2_ppxf, obj_fit_mask_gndf, weights_gndf, $
                mult_poly_coeff_gndf, bestfit_gndf, chi2_gndf, eml_model, stellar_kinematics, $
                stellar_kinematics_err, emission_line_kinematics, emission_line_kinematics_err, $
                emission_line_omitted, emission_line_kinematics_individual, $
                emission_line_kinematics_individual_err, emission_line_intens, $
                emission_line_intens_err, emission_line_fluxes, emission_line_fluxes_err, $
                emission_line_EW, emission_line_EW_err, reddening_output, reddening_output_err, $
                analysis_par=analysis_par, default_velocity=default_velocity, $
                default_dispersion=default_dispersion, eml_par=eml_par, $
                range_v_star=range_v_star, range_s_star=range_s_star, range_v_gas=range_v_gas, $
                range_s_gas=range_s_gas, wavelength_input=wavelength_input, $
                external_library=external_library, region_mask=region_mask, $
                wave_range_analysis=wave_range_analysis, ppxf_only=ppxf_only, version=version, $
                use_previous_guesses=use_previous_guesses, fix_star_kin=fix_star_kin, $
                fix_gas_kin=fix_gas_kin, quiet=quiet, rest_frame_log=rest_frame_log, $
                oversample=oversample, plot=plot, dbg=dbg

        version_module = '0.9'                  ; Version number
        if n_elements(version) ne 0 then begin  ; If version is defined
            version = version_module            ; ... set it to the module version
            return                              ; ... and return without doing anything
        endif

        ; Check that the number of elements in region_mask is even!
        if n_elements(region_mask) MOD 2 ne 0 then $
            message, 'Input region_mask does not have an even number of entries!'

        velScale=MDAP_VELOCITY_SCALE(obj_wave, /log10)  ; km/s/pixel of the input spectra

        ; Set the number of moments to fit for the stellar kinematics
        ; TODO: If running through DAP should never actual enter the if
        ; statement
        moments=analysis_par.moments
        if n_elements(moments) eq 0 then $
            moments=2                                   ; Defaults to fit just V,sigma

        ; Check the dimesions of the input arrays:
        ;       - ivar and mask do not need to exist
        ;       - if they do, they must have the same number of spectra as flux
        ;       - All arrays must have the same number of spectral channels
        ;       - wave and sres (for object spectra) must be 1D
        ;       - flux, ivar, and mask arrays must be 2D, even if there is only one spectrum (row)
        MDAP_SPECTRAL_FITTING_CHECK_DIMENSIONS, obj_flux, obj_ivar, obj_mask, obj_wave, obj_sres
        MDAP_SPECTRAL_FITTING_CHECK_DIMENSIONS, tpl_flux, tpl_ivar, tpl_mask, tpl_wave

        ; TODO: Get rid of this option.  Force everything to be on the
        ; same wavelength system as the input object wavelength vector.
        MDAP_SPECTRAL_FITTING_SET_OUTPUT_WAVELENGTHS, obj_wave, tpl_wave, wavelength_output, $
                                                      same_wave_as_input, $
                                                      wavelength_input=wavelength_input, $
                                                      rest_frame_log=rest_frame_log

        ; Intialize some useful dimensions
        sz=size(obj_flux)
        nobj = sz[1]                            ; Number of object spectra

        sz=size(tpl_flux)
        ntpl = sz[1]                            ; Number of template spectra
        ntpl_w = sz[2]                          ; Number of spectral channels in the template

        neml = n_elements(eml_par)              ; Number of emission lines: 0 if eml_par not defined

        nw = n_elements(wavelength_output)      ; Number of spectral channels FOR OUTPUT

        print, 'Number of object spectra: ', nobj
        print, 'Number of template spectra: ', ntpl
        print, 'Number of emission lines: ', neml
        print, 'Number of spectral channels for output: ', nw

        ; Initialize the output arrays
        obj_fit_mask_ppxf = make_array(nobj, nw, /double, value=1.0d)
        weights_ppxf = dblarr(nobj, ntpl)
        degree=analysis_par.degree
        if n_elements(degree) ne 0 then begin
            if degree ge 0 then $
                add_poly_coeff_ppxf = dblarr(nobj, degree+1)
        endif
        mdegree=analysis_par.mdegree
        if n_elements(mdegree) ne 0 then begin
            if mdegree gt 0 then $
                mult_poly_coeff_ppxf = dblarr(nobj, mdegree)
        endif
        bestfit_ppxf = dblarr(nobj, nw)
        chi2_ppxf = dblarr(nobj)

        if ppxf_only eq 0 then begin
            obj_fit_mask_gndf = make_array(nobj, nw, /double, value=1.0d)
            weights_gndf = dblarr(nobj, ntpl)
            if analysis_par.reddening_order eq 0 && n_elements(mdegree) ne 0 then begin
                if mdegree gt 0 then $
                    mult_poly_coeff_gndf = dblarr(nobj, mdegree)
            endif
            bestfit_gndf = dblarr(nobj, nw)
            chi2_gndf = dblarr(nobj)
            eml_model = dblarr(nobj, nw)
        endif

        stellar_kinematics = dblarr(nobj, moments)
        stellar_kinematics_err= dblarr(nobj, moments)

        ; TODO: Allow for more than 2 kinematic moments for the gas?
        if ppxf_only eq 0 then begin
            emission_line_kinematics = dblarr(nobj, 2)
            emission_line_kinematics_err = dblarr(nobj, 2)
            emission_line_omitted = make_array(nobj, neml, /int, value=1)       ; 0/1 -> fit/omitted
            emission_line_kinematics_individual = dblarr(nobj, neml, 2)
            emission_line_kinematics_individual_err = dblarr( nobj, neml, 2)
            emission_line_intens = dblarr(nobj, neml)
            emission_line_intens_err = dblarr(nobj, neml)
            emission_line_fluxes = dblarr(nobj, neml)
            emission_line_fluxes_err = dblarr(nobj, neml)
            emission_line_EW = dblarr(nobj, neml)
            emission_line_EW_err = dblarr(nobj, neml)
            
            reddening_output = dblarr(nobj, 2)
            reddening_output_err = dblarr(nobj, 2)
        endif

        ; Initialize the starting guesses
;       MDAP_SPECTRAL_FITTING_INIT_GUESSES, star_kin_starting_guesses, nobj, /h3h4, $
;                                           use_defaults=use_defaults
;       if use_defaults eq 1 && n_elements(default_velocity) ne 0 then $
;           star_kin_starting_guesses[*,0] = default_velocity
;       MDAP_SPECTRAL_FITTING_INIT_GUESSES, gas_kin_starting_guesses, nobj, $
;                                           use_defaults=use_defaults
;       if use_defaults eq 1 && n_elements(default_velocity) ne 0 then $
;           gas_kin_starting_guesses[*,0] = default_velocity

        ; Starting guesses are both of size [nobj,2]
        MDAP_SPECTRAL_FITTING_INIT_GUESSES, star_kin_starting_guesses, nobj, default_velocity, $
                                            default_dispersion
        MDAP_SPECTRAL_FITTING_INIT_GUESSES, gas_kin_starting_guesses, nobj, default_velocity, $
                                            default_dispersion

        ; Limit the fitted pixels
        fit_indx = MDAP_SPECTRAL_FITTING_MASK(obj_wave, tpl_wave, velScale, wave_range_analysis, $
                                              star_kin_starting_guesses)

;-----------------------------------------------------------------------
; Moved to pro/utilities/mdap_spectral_fitting_mask.pro
;       ; 1. Apply the wavelength range limit, if provided
;       now=(size(obj_wave))[1]
;       if n_elements(wave_range_analysis) then begin
;           MDAP_SELECT_WAVE, obj_wave, wave_range_analysis, fit_indx, count=count
;           if count eq 0 then $
;               message, 'wave_range_analysis selects no pixels!'
;       endif else $
;           fit_indx = indgen(now)
;
;       c=299792.458d                                   ; Speed of light in km/s
;       maxvel = 400.                                   ; Maximum expected motions (km/s)
;       z_min = (min(star_kin_starting_guesses[*,0])-maxvel)/c  ; Minimum redshift
;       z_max = (max(star_kin_starting_guesses[*,0])+maxvel)/c  ; Maximum redshift
;
;       ; 2. If the number of template pixels is not >= number of fitted galaxy pixels,
;       ;    further limit the blue and red edges of the galaxy spectra
;       now=(size(obj_wave[fit_indx]))[1]               ; Number of object pixels
;       ntw=(size(tpl_wave))[1]                         ; Number of template pixels
;       if ntw lt now then begin
;
;           ; Expected minimum redshift of all spectra
;           z_min = (min(star_kin_starting_guesses[*,0])-maxvel)/c
;
;           ; Indices of wavelengths redward of the redshifted template
;           indx=where(obj_wave gt tpl_wave[0]*(1. + z_min), count)
;           if count eq 0 then $
;               message, 'No overlapping wavelengths between galaxy and template!'
;
;           now_=(size(obj_wave[indx]))[1]
;           if ntw lt now_ then $
;               indx=indx[0:ntw-1]                      ; Truncate the red edge as well
;           ; Merge with current index
;           fit_indx = MDAP_SET_INTERSECTION(indx, temporary(fit_indx), count=count)
;           if fit_indx[0] eq -1 then $
;           if count eq 0 then $
;               message, 'No overlapping wavelengths between galaxy and template!'
;       endif
;       now=(size(obj_wave[fit_indx]))[1]               ; Update number of object pixels
;
;       ; 3. Limit wavelength range to avoid aliasing problems in the template convolution
;       ; TODO: Allow these to be input parameters?
;       sigomit = 6.                                    ; Number of sigma to mask (+/- 3)
;       nalias = fix(sigomit*maxvel/velScale)           ; Number of pixels to mask (maxvel~maxsig)
;       print, 'Masking '+MDAP_STC(nalias,/integer)+' pixels at either end of the spectrum to' $
;              + ' avoid convolution aliasing.'
;       ; Mask to the range that should be unaffected by alias errors
;       wave_range_tpl_unalias = [ tpl_wave[nalias]*(1+z_max), tpl_wave[ntw-nalias-1]/(1+z_min) ]
;       MDAP_SELECT_WAVE, obj_wave, wave_range_tpl_unalias*(1. + z_min), indx
;       ; Merge with current index
;       fit_indx = MDAP_SET_INTERSECTION(indx, temporary(fit_indx), count=count)
;       if count eq 0 then $
;           message, 'No intersection between wave_range_tpl_unalias and fit_indx!'
;-----------------------------------------------------------------------

        ; Truncate the vectors to the valid pixels
        obj_wave_lim = obj_wave[fit_indx]
        obj_sres_lim = obj_sres[fit_indx]
        obj_flux_lim = obj_flux[*,fit_indx]
        obj_ivar_lim = obj_ivar[*,fit_indx]
        obj_mask_lim = obj_mask[*,fit_indx]

;        print, n_elements(tpl_wave), n_elements(obj_wave_lim)
;        plot, obj_wave, obj_flux[0,*]
;        oplot, obj_wave_lim, obj_flux_lim[0,*], color=200
;        stop

        ; Set to ignore emission lines that are not within the fitted wavelength range
        eml_par_lim = eml_par
        if n_elements(eml_par) ne 0 && ppxf_only eq 0 then begin
            MDAP_CHECK_EMISSION_LINES, eml_par_lim, obj_wave_lim, $
                                       velocity=gas_kin_starting_guesses[0,0]
        endif

        print, 'Number of emission lines to fit: ', n_elements(eml_par_lim)

        ; On input, pPXF expects the template and galaxy to have the same
        ;   wavelength coordinate system.  To account for the fact that this may not
        ;   be true, we calculate the velocity offset of the template with respect
        ;   to the initial wavelength of the galaxy.  That is, we determine the
        ;   pseudo-velocity shift required to match the observed wavelengths of the
        ;   template library and the galaxy spectrum.
        voff = (alog10(tpl_wave[0]) - alog10(obj_wave_lim[0]))*velScale / $
               (alog10(obj_wave_lim[1])-alog10(obj_wave_lim[0]))

        ; Initialize the starting guesses vector used by MDAP_GANDALF_WRAP,
        ;   which DOES NOT change these values during its execution meaning this
        ;   only needs to be done once.
        start = MDAP_SPECTRAL_FITTING_START_KIN(reform(star_kin_starting_guesses[0,*]), $
                                                reform(gas_kin_starting_guesses[0,*]) )
        ; the returned vector has *four* elements!!

        ; Make sure that the velocity dispersion of the stars is larger
        ; than the minimum allowed by pPXF:
        if start[1] lt velScale/10. then $
            start[1] = velScale/9.

        if ~keyword_set(quiet) then begin
            print, "Voff: ", voff
            print, 'Number of spectra to fit: ', nobj
            print, "Start: ", start
        endif

        ; Initialize the instrumental dispersion by interpolating the value from
        ;   the spectral resolution at the guess redshifts of the emission lines
        instr_disp = MDAP_INSTRUMENTAL_DISPERSION(obj_wave_lim, obj_sres_lim, eml_par_lim.lambda, $
                                                  start[2], $
                                                  zero_instr_disp=analysis_par.zero_instr_disp)

        ; Begin loop over all the object spectra
        c=299792.458d                                   ; Speed of light in km/s (needed below)
        if keyword_set(dbg) then $
            nobj=1                                      ; Only fit the first spectrum
        for i=0,nobj-1 do begin

            if ~keyword_set(quiet) then $
                print, 'Fitting spectrum '+mdap_stc(i+1,/integer)+' of '+mdap_stc(nobj,/integer)
          
            ; If not fitting the reddening, (re)set the values
            if analysis_par.reddening_order gt 0 then $
                ebv=reddening

            ; Update the starting guesses and instrumental dispersion with the
            ;   previous fit values, if requested
            if i gt 0 && keyword_set(use_previous_guesses) then begin
                ; Returned solution fro GANDALF_WRAP has 9 elements:
                ; 0:5 = stellar kin
                ; 6 = chi-square
                ; 7:8 = gas kin
                ; start vector should only have:
                ; 0:1 stellar v and sigma
                ; 2:3 gas v and sigma
                start = MDAP_SPECTRAL_FITTING_START_KIN(sol[0:1], sol[7:8])

                instr_disp = MDAP_INSTRUMENTAL_DISPERSION(obj_wave_lim, obj_sres_lim, $
                                                          eml_par_lim.lambda, start[2], $
                                                    zero_instr_disp=analysis_par.zero_instr_disp)
            endif

            ; TODO: correct for galactic reddening, if provided
;           galaxy = reform(obj_flux_lim[i,*])
;           if keyword_set(MW_extinction) then begin
;               dereddening_attenuation = DUST_CALZETTI(l0_gal, lstep_gal, n_elements(galaxy_), $
;                                                       -MW_extinction, 0.0d, /log10)
;               dereddening_attenuation = mdap_dust_calzetti(log_0_gal, log_step_gal, $
;                                                            n_elements(galaxy_), -MW_extinction, $
;                                                            0.0d)
;               galaxy = galaxy*temporary(dereddening_attenuation)
;           endif

            ; TODO: Need to make these status flags actually do something!
            ppxf_status = 1                     ; Initialize the status: 0=success, 1=failure
            gandalf_status = 1

            ; TODO: Create a default behavior for things like mdegree, degree, moments, etc...
            ; TODO: Need to reform() because of expected orientation for pPXF and GANDALF
            ; TODO: Include a log10 keyword
            MDAP_GANDALF_WRAP, tpl_wave, tpl_flux, obj_wave_lim, reform(obj_flux_lim[i,*]), $
                               reform(obj_ivar_lim[i,*]), reform(obj_mask_lim[i,*]), velScale, $
                               start, fitted_pixels_ppxf, weights_ppxf_i, add_poly_coeff_ppxf_i, $
                               mult_poly_coeff_ppxf_i, bestfit_ppxf_i, chi2_ppxf_i, $
                               fitted_pixels_gndf, weights_gndf_i, mult_poly_coeff_gndf_i, $
                               bestfit_gndf_i, chi2_gndf_i, eml_model_i, sol, err, gas_intens, $
                               gas_intens_err, gas_vel, gas_vel_err, gas_sig, gas_sig_err, $
                               gas_flux, gas_flux_err, gas_ew, gas_ew_err, eml_par=eml_par_lim, $
                               vsyst=voff, bias=bias, mdegree=mdegree, degree=degree, $
                               moments=moments, reddening=ebv, range_v_star=range_v_star, $
                               range_s_star=range_s_star, range_v_gas=range_v_gas, $
                               range_s_gas=range_s_gas, region_mask=region_mask, $
                               external_library=external_library, int_disp=instr_disp, $
                               err_reddening=err_reddening, ppxf_status=ppxf_status, $
                               gandalf_status=gandalf_status, oversample=oversample, $
                               for_errors=1, fix_star_kin=fix_star_kin, fix_gas_kin=fix_gas_kin, $
                               ppxf_only=ppxf_only, quiet=quiet, plot=plot

            ; On output the best fit and emline model are at obj_wave_lim

;           indx = where(weights_gndf_i[0:ntpl-1] gt 0)
;           print, 'Templates with non-zero weights: ', indx

;           print, 'Solution: ', sol
;           print, 'Error:    ', err

;           plot, obj_wave_lim, bestfit
;           oplot, obj_wave_lim, eml_model

            ; Impose some limits on the results
            ; TODO: Put this into a procedure

            ; If the error vector is flat (i.e. errors are not
            ; reliable), rescale the formal errors by sqrt(chi2/dof), as
            ; instructed by MPFIT and pPXF.  Output sol should ALWAYS
            ; have chi2/dof at index=6
            if min(obj_ivar_lim[i,*]) eq max(obj_ivar_lim[i,*]) then begin
                err = err * sqrt(sol[6])
                gas_vel_err = gas_vel_err * sqrt(sol[6]) 
                gas_sig_err = gas_sig_err * sqrt(sol[6]) 
            endif

            ; If there are fewer than 3 fitted lines with intensity S/N of > 1,
            ; delete *all* kinematics measurements for the gas.

            ; TODO: Let this be user-defined?  Seems somewhat extreme.

            line_ston = abs(gas_intens/gas_intens_err)
            indx = where(gas_intens gt 0 and line_ston gt 1.0d)
            if n_elements(indx) lt 3 then begin

;               print, 'Bad lines!'

                ; TODO: This is inconsistent with the default behavior in GANDALF_WRAP
                gas_vel[*] = -9999.0d       ; Individual gas velocities
                gas_vel_err[*] = -9999.0d   ; Individual gas velocity errors
                gas_sig[*] = -9999.0d       ; Individual gas velocity dispersions
                gas_sig_err[*] = -9999.0d   ; Individual gas velocity disperison errors

                sol[7] = -9999.0d           ; Mean Gas velocity
                err[6] = -9999.0d           ; Mean Gas velocity error

                sol[8] = -9999.0d           ; Mean Gas velocity dispersion
                err[7] = -9999.0d           ; Mean Gas velocity dispersion error
            endif             

            ; plots for checks... remove these lines when running on remote server
;           plot, obj_wave, obj_flux[i,*], title='GANDALF + '+string(i), xrange=[3600,10350], $
;                 xstyle=1
;           oplot, obj_wave, bestfit, color=200
;           if n_elements(fitted_pixels) ne 0 then $
;               oplot, obj_wave[fitted_pixels], bestfit[fitted_pixels], color=200
    
            ; Store output stellar kinematics
            if keyword_set(fix_star_kin) then begin
                stellar_kinematics[i,0:1] = start[0:1]
                stellar_kinematics_err[i,*] = -9999.0d
            endif else begin
                stellar_kinematics[i,*]=sol[0:moments-1]
                stellar_kinematics_err[i,*]=err[0:moments-1]
            endelse

            ; Store output gas kinematics
            if ppxf_only eq 0 then begin

                emission_line_kinematics[i,*]=sol[7:8]          ; V,sigma for emission lines
                emission_line_kinematics_err[i,*]=err[6:7]      ; V,sigma errors for emission lines

                i_f = where(eml_par_lim.action eq 'f', count)
;               if i_f[0] ne -1 then begin
                if count ne 0 then begin
                    emission_line_omitted[i,i_f] = 0            ; Flag fitted lines
                    emission_line_intens[i,i_f]= gas_intens     ; Intensity values, dereddened
                    emission_line_intens_err[i,i_f]= gas_intens_err     ; Error in intensity
                    emission_line_EW[i,i_f]= gas_ew             ; Emission equivalent width (EW)
                    emission_line_EW_err[i,i_f]= gas_ew_err     ; error in EW
                    emission_line_fluxes[i,i_f]= gas_flux       ; Gas fluxes, dereddened
                    emission_line_fluxes_err[i,i_f] = gas_flux_err      ; errors in gas fluxes
                    emission_line_kinematics_individual[i,i_f,0]=gas_vel        ; vel for lines
                    emission_line_kinematics_individual_err[i,i_f,0]=gas_vel_err; vel errors
                    emission_line_kinematics_individual[i,i_f,1]=gas_sig        ; sigma for lines
                    emission_line_kinematics_individual_err[i,i_f,1]=gas_sig_err; sigma error
                endif

                if n_elements(ebv) eq 2 then begin
                    reddening_output[i,*]= ebv
                    reddening_output_err[i,*]= err_reddening
                endif else if n_elements(ebv) eq 1 then begin
                    reddening_output[i,*]= [ebv[0],-9999.0d]
                    reddening_output_err[i,*]= [err_reddening[0],-9999.0d]
                endif else begin
                    reddening_output[i,*]= [-9999.0d,-9999.0d]
                    reddening_output_err[i,*]= [-9999.0d,-9999.0d]
                endelse

            endif

            ; Convert the "velocities" (velScale*pixel_shift) to redshifts
            ; TODO: This means that one can no longer create the correct LOSVD
            ;       using these velocities!!!
            v = stellar_kinematics[i,0]
            ve = stellar_kinematics_err[i,0]
            MDAP_CONVERT_PIXEL_KINEMATICS, v, ve
            stellar_kinematics[i,0] = v
            stellar_kinematics_err[i,0] = ve

            if ppxf_only eq 0 then begin
                v = emission_line_kinematics[i,0]
                ve = emission_line_kinematics_err[i,0]
                MDAP_CONVERT_PIXEL_KINEMATICS, v, ve
                emission_line_kinematics[i,0] = v
                emission_line_kinematics_err[i,0] = ve
                
                if i_f[0] ne -1 then begin
                    v = emission_line_kinematics_individual[i,i_f,0]
                    ve = emission_line_kinematics_individual_err[i,i_f,0]
                    MDAP_CONVERT_PIXEL_KINEMATICS, v, ve
                    emission_line_kinematics_individual[i,i_f,0] = v
                    emission_line_kinematics_individual_err[i,i_f,0] = ve
                endif
            endif

;           print, 'fit mask'
            ; Unmask the fitted pixels
            if fitted_pixels_ppxf[0] ne -1 then $
                obj_fit_mask_ppxf[i,fit_indx[fitted_pixels_ppxf]] = 0.0d
            if fitted_pixels_gndf[0] ne -1 then $
                obj_fit_mask_gndf[i,fit_indx[fitted_pixels_gndf]] = 0.0d

            weights_ppxf[i,*] = weights_ppxf_i                          ; Save pPXF weights
            chi2_ppxf[i] = chi2_ppxf_i                                  ; Save pPXF chi^2
            if n_elements(degree) ne 0 then begin
                if degree gt 0 then $
                    add_poly_coeff_ppxf[i,*] = add_poly_coeff_ppxf_i    ; Save additive poly
            endif
            if n_elements(mdegree) ne 0 then begin
                if mdegree gt 0 then $
                    mult_poly_coeff_ppxf[i,*] = mult_poly_coeff_ppxf_i  ; Save multiplicative poly
            endif

            if ppxf_only eq 0 then begin
                weights_gndf[i,*] = weights_gndf_i[0:ntpl-1]            ; Save GANDALF weights
                chi2_gndf[i] = chi2_gndf_i                              ; Save GANDALF chi^2
                if analysis_par.reddening_order eq 0 && n_elements(mdegree) ne 0 then begin
                    if mdegree gt 0 then $
                        mult_poly_coeff_gndf[i,*] = mult_poly_coeff_gndf_i  ; Save multiplicative
                endif
            endif

            ;--------------------
            ; TODO: This whole procedure, generating the best-fit template, the
            ; LOSVD, and the convolution should use routines that are common with
            ; the fitting algorithm!

            ; bf_template = (templates # weights[0:sztempl[2]-1])
            ; TODO: Does the size need to be specified, or is it already correct (=ntpl)?

            ; V offset between template and galaxy initial wave NOT included in SOL
            ; so bf_template and bf_template_losvd have wavelengths corresponding to tpl_wave

            ; TODO: how is bf_template used?  Used in block 5 for index measurements
;           print, size(reform(weights[0:ntpl-1]))
;           print, size(tpl_flux)
            ; TODO: Check for errors
;           best_template_i = ppxf_only eq 1 ? reform(weights_ppxf[i,*] # tpl_flux) : $
;                                              reform(weights_gndf[i,*] # tpl_flux)
;           print, size(bf_template)

            ; TODO: So far MDAP_GET_BROADENED_TEMPLATE is ONLY used here!
;           best_template_losvd_i = MDAP_GET_BROADENED_TEMPLATE(best_template_i, sol, velscale, $
;                                                               moments=moments, $
;                                                               oversample=oversample)
            ;--------------------
            ; TODO: Put this all in a procedure!

            ; Output/input wavelengths are the same; copy the spectra and go to next spectrum
            if same_wave_as_input eq 1 then begin

                ; Best-fit PPXF model
                bestfit_ppxf[i,fit_indx] = bestfit_ppxf_i

                if ppxf_only eq 0 then begin
                    bestfit_gndf[i,fit_indx] = bestfit_gndf_i    ; Best-fit GANDALF model
                    eml_model[i,fit_indx] = eml_model_i ; Best-fit emission-line-only model
                endif

                continue                                        ; Go to next spectrum
            endif
            
            ; Generate the output spectra
            ;rf_gal_lam = exp(loglam_gal-sol[0]/velscale*(log_step_gal))
            if keyword_set(rest_frame_log) then begin
                obj_wave_lim_ = obj_wave_lim/(1.0d + sol[0]/c)          ; Shift wave to V=0
            endif

            ; Best-fit PPXF model mapped to the output wavelength vector
            bestfit_ppxf[i,fit_indx] = interpol(bestfit_ppxf_i, obj_wave_lim_, wavelength_output)

            ; Best-fit GANDALF model and emission-line only model
            if ppxf_only eq 0 then begin
                bestfit_gndf[i,fit_indx] = interpol(bestfit_gndf_i,obj_wave_lim_,wavelength_output)
                eml_model[i,fit_indx] = interpol(eml_model_i, obj_wave_lim_, wavelength_output)
            endif

        endfor
END


