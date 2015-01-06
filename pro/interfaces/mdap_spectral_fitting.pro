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
;                              best_template, best_template_losvd_conv, stellar_kinematics, $
;                              stellar_kinematics_err, emission_line_kinematics, $
;                              emission_line_kinematics_err, emission_line_omitted, $
;                              emission_line_kinematics_individual, $
;                              emission_line_kinematics_individual_err, emission_line_intens, $
;                              emission_line_intens_err, emission_line_fluxes, $
;                              emission_line_fluxes_err, emission_line_EW, emission_line_EW_err, $
;                              reddening_output, reddening_output_err, analysis_par=analysis_par, $
;                              star_kin_starting_guesses=star_kin_starting_guesses, $
;                              gas_kin_starting_guesses=gas_kin_starting_guesses, eml_par=eml_par, $
;                              range_v_star=range_v_star, range_s_star=range_s_star, $
;                              range_v_gas=range_v_gas, range_s_gas=range_s_gas, $
;                              wavelength_input=wavelength_input, $
;                              external_library=external_library, region_mask=region_mask, $
;                              wave_range_analysis=wave_range_analysis, $
;                              ppxf_only=ppxf_only, /use_previous_guesses, /fix_star_kin, $
;                              /fix_gas_kin, /quiet, /rest_frame_log, /plot, /dbg, version=version
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
; OPTIONAL INPUTS:
;       =====================================OLD========================================
;       extra_inputs strarr[]
;
;       TODO: Change this to a structure that sets all possible extra parameters!
;
;               A string array containing other inputs to be used in the fitting
;               procedure, such as the polynomial order. The list is initialized
;               using the execute command:
;
;                   for i = 0, n_elements(extra_inputs)-1 do $
;                       d = execute(extra_inputs[i])
;
;               where, for example,
;
;                    extra_inputs=['MOMENTS=2', 'DEGREE=-1', 'BIAS=0', 'reddening=0', $
;                                  'LAMBDA=exp(loglam_gal)']       
;
;               Warning: The reddening (stars and/or gas) fit is optional, and
;               it is performed by the Gandalf module. If the reddening fit is
;               required, MDEGREE and DEGREE are used as the input value for the
;               pPXF run, but automatically set to 0 and -1 respectively in the
;               Gandalf execution.
;
;               TODO: Are these basically all the optional parameters/keywords
;               for pPXF?
;
;               MDEGREE
;               BIAS
;               DEGREE
;               reddening
;               lambda
;       =====================================OLD========================================
;
;       analysis_par AnalysisPar structure
;               A structure that defines parameters used by PPXF and GANDALF in
;               the fitting procedure.  See its definition in
;               MDAP_DEFINE_ANALYSIS_PAR.pro
; 
;       star_kin_starting_guesses dblarr[N][4]
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
;               TODO: Don't understand this.  Is it necessary?
;
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
;       best_template dblarr[N][QQ]
;               The best-fitting template (sum of the weighted template in the
;               library) for each of the N galaxy spectra, sampled over
;               wavelength_output (rest frame wavelength).  The weights used
;               are:
;                       weights_ppxf if ppxf_only=1 or gandalf_status=1
;                       weights_gndf otherwise
;
;       best_template_losvd_conv dblarr[N][QQ]
;               The best-fitting template (sum of the weighted templates in the
;               library) for each of the N galaxy spectra, convolved with the
;               best-fitting LOSVD and sampled over wavelength_output (rest
;               frame wavelength).  These are the best_template spectra
;               convolved with the LOSVDs in stellar_kinematics.
;
;       stellar_kinematics dblarr[N][5]
;               The best-fit V, sigma, h3, h4, and chi2/DOF for each of the N
;               fitted input galaxy spectra.  If \fix_star_kin is set, the array
;               is not defined.  TODO: Which is it?  Is it set to the input
;               values (as described above) or undefined?
;
;       stellar_kinematics_err dblarr[N][4]
;               Estimates of the errors in the best-fit V, sigma, h3, h4 for the
;               stellar kinematics for each of the N fitted input galaxy
;               spectra.  TODO: How are these error generated?
;
;       emission_line_kinematics dblarr[N][2]
;               The best-fit V and sigma for the emission lines in the N galaxy
;               spectra. If \fix_gas_kin is set, the array is not defined.
;               TODO: Which is it?  Is it set to the input values (as described
;               above) or undefined?  Is chi^2 included?
;
;       emission_line_kinematics_err dblarr[N][2]
;               Estimates of the errors in the best-fit V and sigma for the gas
;               kinematics for each of the N fitted input galaxy spectra.  TODO:
;               How are these error generated?
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
;       24 Sep 2014: (KBW) Formatting + many edits
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
;-
;------------------------------------------------------------------------------

PRO MDAP_SPECTRAL_FITTING_CHECK_DIMENSIONS, $
                flux, ivar, mask, wave, sres

        ; The only things that have to be provided are the flux and the wavelength
        if n_elements(flux) eq 0 or n_elements(wave) eq 0 then $
            message, 'Flux and/or wavelength arrays not defined!'

        sz=size(flux)                           ; All array must match the size of flux

        sz_ = size(wave)
        if sz_[0] ne 1 then $
            message, 'wave must be one dimensional!'
        if (sz_[0] eq 1 and sz[2] ne sz_[1]) or (sz_[0] eq 2 and sz[2] ne sz_[2]) then $
            message, 'flux and wave must have the same number of spectral channels.'

        if n_elements(sres) ne 0 then begin 
            sz_ = size(sres)
            if sz_[0] ne 1 then $
                message, 'sres must be one dimensional!'
            if (sz_[0] eq 1 and sz[2] ne sz_[1]) or (sz_[0] eq 2 and sz[2] ne sz_[2]) then $
                message, 'flux and sres must have the same number of spectral channels.'
        endif

        if n_elements(ivar) eq 0 then begin
            sz_ = size(ivar)
            if sz[0] ne sz_[0] or sz[1] ne sz_[1] or (sz_[0] eq 2 and sz[2] ne sz_[2]) then $
                message, 'flux and ivar sizes must match!'
        endif

        if n_elements(mask) eq 0 then begin
            sz_ = size(mask)
            if sz[0] ne sz_[0] or sz[1] ne sz_[1] or (sz_[0] eq 2 and sz[2] ne sz_[2]) then $
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

PRO MDAP_SPECTRAL_FITTING_INIT_GUESSES, $
                guesses, ns, h3h4=h3h4

        use_defaults = 1
        if n_elements(guesses) ne 0 then begin
            sz=size(guesses)
            if sz[0] ne 2 then begin
                print, 'Kinematic guess arrays must be two-dimensional!  Using default.'
            endif else if (keyword_set(h3h4) and sz[2] ne 4) $
                           or (~keyword_set(h3h4) and sz[2] ne 2) then begin
                print, 'Incorrect number of kinematic parameters!  Using default.'
            endif else $
                use_defaults=0
        endif

        if use_defaults eq 1 then begin
            if keyword_set(h3h4) then begin
                ndim = 4
            endif else $
                ndim = 2
            guesses = dblarr(ns, ndim)          ; set the array size and initialize v=h3=h4=0.0d
            guesses[*,1] = 50.0d                ; ... but initalize sigma=50 km/s
        endif
END

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
                tpl_mask, wavelength_output, obj_fit_mask_ppxf, weights_ppxf, add_poly_coeff_ppxf, $
                mult_poly_coeff_ppxf, bestfit_ppxf, chi2_ppxf, obj_fit_mask_gndf, weights_gndf, $
                mult_poly_coeff_gndf, bestfit_gndf, chi2_gndf, eml_model, stellar_kinematics, $
                stellar_kinematics_err, emission_line_kinematics, emission_line_kinematics_err, $
                emission_line_omitted, emission_line_kinematics_individual, $
                emission_line_kinematics_individual_err, emission_line_intens, $
                emission_line_intens_err, emission_line_fluxes, emission_line_fluxes_err, $
                emission_line_EW, emission_line_EW_err, reddening_output, reddening_output_err, $
                analysis_par=analysis_par, star_kin_starting_guesses=star_kin_starting_guesses, $
                gas_kin_starting_guesses=gas_kin_starting_guesses, eml_par=eml_par, $
                range_v_star=range_v_star, range_s_star=range_s_star, range_v_gas=range_v_gas, $
                range_s_gas=range_s_gas, wavelength_input=wavelength_input, $
                external_library=external_library, region_mask=region_mask, $
                wave_range_analysis=wave_range_analysis, $
                use_previous_guesses=use_previous_guesses, fix_star_kin=fix_star_kin, $
                fix_gas_kin=fix_gas_kin, quiet=quiet, rest_frame_log=rest_frame_log, $
                plot=plot, ppxf_only=ppxf_only, dbg=dbg, version=version

        version_module = '0.9'                  ; Version number
        if n_elements(version) ne 0 then begin  ; If version is defined
            version = version_module            ; ... set it to the module version
            return                              ; ... and return without doing anything
        endif

        ; Check that the number of elements in region_mask is even!
        if n_elements(region_mask) MOD 2 ne 0 then $
            message, 'Input region_mask does not have an even number of entries!'

        velScale=MDAP_VELOCITY_SCALE(obj_wave, /log10)  ; km/s/pixel of the input spectra

;       print, 'velscale: ', velScale

        ; TODO: Need to know the moments to set size of the stellar_kinematics vector
        ; TODO: Does this need to be done here.  Probably not because pPXF does the same thing.
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

        MDAP_SPECTRAL_FITTING_SET_OUTPUT_WAVELENGTHS, obj_wave, tpl_wave, wavelength_output, $
                                                      same_wave_as_input, $
                                                      wavelength_input=wavelength_input, $
                                                      rest_frame_log=rest_frame_log

;       ; Declare the extra input parameters
;       ; TODO: What is the full list of possible values
;       if n_elements(extra_inputs) ne 0 then begin
;           for k = 0, n_elements(extra_inputs)-1 do begin
;               if strlen(extra_inputs[k]) eq 0 then $
;                   continue
;
;               d = execute(extra_inputs[k])
;               if ~keyword_set(quiet) then $
;                   print, 'setting variable '+extra_inputs[k]
;
;           endfor
;       endif

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

;       stop

        ; Initialize the output arrays
        ; TODO: Convert this to a structure?

        obj_fit_mask_ppxf = make_array(nobj, nw, /double, value=1.0d)
        weights_ppxf = dblarr(nobj, ntpl)
        degree=analysis_par.degree
        if n_elements(degree) ne 0 then begin
            if degree gt 0 then $
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
            if analysis_par.reddening_order eq 0 and n_elements(mdegree) ne 0 then begin
                if mdegree gt 0 then $
                    mult_poly_coeff_gndf = dblarr(nobj, mdegree)
            endif
            bestfit_gndf = dblarr(nobj, nw)
            chi2_gndf = dblarr(nobj)
            eml_model = dblarr(nobj, nw)
        endif

;       best_template = dblarr(nobj, nw)
;       best_template_losvd_conv = dblarr(nobj, nw)

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
        MDAP_SPECTRAL_FITTING_INIT_GUESSES, star_kin_starting_guesses, nobj, /h3h4
        MDAP_SPECTRAL_FITTING_INIT_GUESSES, gas_kin_starting_guesses, nobj

;       print, star_kin_starting_guesses
;       print, gas_kin_starting_guesses

        ;-----------------------------------------------------------------------
        ; Limited the fitted pixels
        ;-----------------------------------------------------------------------
        ; 1. Apply the wavelength range limit, if provided
        now=(size(obj_wave))[1]
        if n_elements(wave_range_analysis) then begin
            MDAP_SELECT_WAVE, obj_wave, wave_range_analysis, fit_indx, count=count
;           if fit_indx[0] eq -1 then $
            if count eq 0 then $
                message, 'wave_range_analysis selects no pixels!'
        endif else $
            fit_indx = indgen(now)

        c=299792.458d                                   ; Speed of light in km/s
        maxvel = 400.                                   ; Maximum expected motions (km/s)
        z_min = (min(star_kin_starting_guesses[*,0])-maxvel)/c  ; Minimum redshift
        z_max = (min(star_kin_starting_guesses[*,0])+maxvel)/c  ; Maximum redshift

        ; 2. If the number of template pixels is not >= number of fitted galaxy pixels,
        ;    further limit the blue and red edges of the galaxy spectra
        now=(size(obj_wave[fit_indx]))[1]               ; Number of object pixels
        ntw=(size(tpl_wave))[1]                         ; Number of template pixels
        if ntw lt now then begin

            ; Expected minimum redshift of all spectra
            z_min = (min(star_kin_starting_guesses[*,0])-maxvel)/c

            ; Indices of wavelengths redward of the redshifted template
            indx=where(obj_wave gt tpl_wave[0]*(1. + z_min), count)
            if count eq 0 then $
                message, 'No overlapping wavelengths between galaxy and template!'

            now_=(size(obj_wave[indx]))[1]
            if ntw lt now_ then $
                indx=indx[0:ntw-1]                      ; Truncate the red edge as well
            ; Merge with current index
            fit_indx = MDAP_SET_INTERSECTION(indx, temporary(fit_indx), count=count)
;           if fit_indx[0] eq -1 then $
            if count eq 0 then $
                message, 'No overlapping wavelengths between galaxy and template!'
        endif
        now=(size(obj_wave[fit_indx]))[1]               ; Update number of object pixels

        ; 3. Limit wavelength range to avoid aliasing problems in the template convolution
        ; TODO: Allow these to be input parameters?
        sigomit = 6.                                    ; Number of sigma to mask (+/- 3)
        nalias = fix(sigomit*maxvel/velScale)           ; Number of pixels to mask (maxvel~maxsig)
        ; Mask to the range that should be unaffected by alias errors
        wave_range_tpl_unalias = [ tpl_wave[nalias]*(1+z_max), tpl_wave[ntw-nalias-1]/(1+z_min) ]
        MDAP_SELECT_WAVE, obj_wave, wave_range_tpl_unalias*(1. + z_min), indx
        ; Merge with current index
        fit_indx = MDAP_SET_INTERSECTION(indx, temporary(fit_indx), count=count)
;       if fit_indx[0] eq -1 then $
        if count eq 0 then $
            message, 'No intersection between wave_range_tpl_unalias and fit_indx!'

        ; Truncate the vectors to the valid pixels
        obj_wave_lim = obj_wave[fit_indx]
        obj_sres_lim = obj_sres[fit_indx]
        obj_flux_lim = obj_flux[*,fit_indx]
        obj_ivar_lim = obj_ivar[*,fit_indx]
        obj_mask_lim = obj_mask[*,fit_indx]

        ; Set to ignore emission lines that are not within the fitted wavelength range
        eml_par_lim = eml_par
        if n_elements(eml_par) ne 0 and ppxf_only eq 0 then begin
            MDAP_CHECK_EMISSION_LINES, eml_par_lim, obj_wave_lim, $
                                       velocity=gas_kin_starting_guesses[0]
        endif

        print, n_elements(eml_par_lim)

;       stop

        ; On input, pPXF expects the template and galaxy to have the same
        ;   wavelength coordinate system.  To account for the fact that this may not
        ;   be true, we calculate the velocity offset of the template with respect
        ;   to the initial wavelength of the galaxy.  That is, we determine the
        ;   pseudo-velocity shift required to match the observed wavelengths of the
        ;   template library and the galaxy spectrum.

        voff = (alog10(tpl_wave[0]) - alog10(obj_wave_lim[0]))*velScale / $
               (alog10(obj_wave_lim[1])-alog10(obj_wave_lim[0]))

        c=299792.458d                                   ; Speed of light in km/s (needed below)
;       voff = (1.0d - obj_wave_lim[0]/tpl_wave[0])*c   ; Offset velocity in km/s

        ; Initialize the starting guesses vector used by MDAP_GANDALF_WRAP,
        ;   which DOES NOT change these values during its execution meaning this
        ;   only needs to be done once.
        ; TODO: Make this a structure?
        start = MDAP_SPECTRAL_FITTING_START_KIN(reform(star_kin_starting_guesses[0,*]), $
                                                reform(gas_kin_starting_guesses[0,*]) )

        if ~keyword_set(quiet) then begin
            print, "Voff: ", voff
            print, 'Number of spectra to fit: ', nobj
            print, "Start: ", start
        endif

        ; Initialize the instrumental dispersion by interpolating the value from
        ;   the spectral resolution at the guess redshifts of the emission lines
        ; TODO: Make this a function?
        sig2fwhm = 2.0d*sqrt(alog(4.0d))                ; Conversion from sigma to FWHM
        instr_disp = interpol(c/obj_sres_lim, obj_wave_lim, $
                              eml_par_lim.lambda * (1.0d + start[4]/c))/sig2fwhm

        ; Begin loop over all the object spectra
        if keyword_set(dbg) then $
            nobj=1                                      ; Only fit the first spectrum
        for i=0,nobj-1 do begin
;       for i=0,0 do begin
;           i = 0

            if ~keyword_set(quiet) then $
                print, 'Fitting spectrum '+mdap_stc(i+1,/integer)+' of '+mdap_stc(nobj,/integer)
          
            ; TODO: Necessary at every iteration?
            ; If not fitting the reddening, (re)set the values
            if analysis_par.reddening_order gt 0 then $
                ebv=reddening

            ; Update the starting guesses and instrumental dispersion with the
            ;   previous fit values, if requested
            if i gt 0 and keyword_set(use_previous_guesses) then begin
                start = MDAP_SPECTRAL_FITTING_START_KIN(sol[0:3], sol[7:8])
                instr_disp = interpol(c/obj_sres_lim, obj_wave_lim, $
                                      eml_par_lim.lambda * (1.0d + start[4]/c))/sig2fwhm
            endif

;           indx = where(eml_par_lim.action eq 'f')
;           if indx[0] eq -1 then begin
;               neml_f = 0
;           endif else $
;               neml_f = n_elements(indx)

;   galaxy_=galaxy[*,i]
;   noise_=noise[*,i]

 ;  ;-- correct for galactic reddening, if provided
 ;  TODO: RETURN TO THIS, does this need to be done for every spectrum?

 ;  IF KEYWORD_SET(MW_extinction) THEN BEGIN
 ;     ;dereddening_attenuation = DUST_CALZETTI(l0_gal,lstep_gal,n_elements(galaxy_),-MW_extinction,0.0d,/log10)
 ;     dereddening_attenuation = mdap_dust_calzetti(log_0_gal,log_step_gal,n_elements(galaxy_),-MW_extinction,0.0d)
 ;     galaxy_ = galaxy_*temporary(dereddening_attenuation)
 ;  ENDIF

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

            ; If the error vector is flat (i.e. errors are not reliable),
            ;   rescale the formal errors for sqrt(chi2/dof), as instructed
            ;   by mpfit and ppxf.
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
                gas_vel[*] = !VALUES.D_NAN      ; Individual gas velocities
                gas_vel_err[*] = 99.0d          ; Individual gas velocity errors
                gas_sig[*] = !VALUES.D_NAN      ; Individual gas velocity dispersions
                gas_sig_err[*] = 99.0d          ; Individual gas velocity disperison errors

                sol[7] = !VALUES.D_NAN          ; Mean Gas velocity
                err[6] = 99.0d          ; Mean Gas velocity error

                sol[8] = !VALUES.D_NAN          ; Mean Gas velocity dispersion
                err[7] = 99.0d          ; Mean Gas velocity dispersion error
            endif             

            ; plots for checks... remove these lines when running on remote server
;           plot, obj_wave, obj_flux[i,*], title='GANDALF + '+string(i), xrange=[3600,10350], $
;                 xstyle=1
;           oplot, obj_wave, bestfit, color=200
;           if n_elements(fitted_pixels) ne 0 then $
;               oplot, obj_wave[fitted_pixels], bestfit[fitted_pixels], color=200
    
            ; TODO: Just found out that the indexing of arrays in IDL is
            ; super-fucked.  Will need to check that this hasn't messed me
            ; up somewhere...

            ; Store output
            ; TODO: Put input and output into structures?

            stellar_kinematics[i,*]=sol[0:moments-1]            ; V,sigma,h3,h4 for stars
            stellar_kinematics_err[i,*]=err[0:moments-1]        ; V,sigma,h3,h4 errors

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
                    reddening_output[i,*]= [ebv[0],0.0d]
                    reddening_output_err[i,*]= [err_reddening[0],99.0d]
                endif else begin
                    reddening_output[i,*]= [0.0d,0.0d]
                    reddening_output_err[i,*]= [99.0d,99.0d]
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

            weights_ppxf[i,*] = weights_ppxf_i                          ; Set the weights
            chi2_ppxf[i] = chi2_ppxf_i
            if n_elements(degree) ne 0 then begin
                if degree gt 0 then $
                    add_poly_coeff_ppxf[i,*] = add_poly_coeff_ppxf_i
            endif
            if n_elements(mdegree) ne 0 then begin
                if mdegree gt 0 then $
                    mult_poly_coeff_ppxf[i,*] = mult_poly_coeff_ppxf_i
            endif

            if ppxf_only eq 0 then begin
                weights_gndf[i,*] = weights_gndf_i[0:ntpl-1]
                chi2_gndf[i] = chi2_gndf_i
                if analysis_par.reddening_order eq 0 and n_elements(mdegree) ne 0 then begin
                    if mdegree gt 0 then $
                        mult_poly_coeff_gndf[i,*] = mult_poly_coeff_gndf_i
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
;           print, 'save models'
            if same_wave_as_input eq 1 then begin

                ; Best-fit PPXF model
                bestfit_ppxf[i,fit_indx] = bestfit_ppxf_i

                if ppxf_only eq 0 then begin
                    bestfit_gndf[i,fit_indx] = bestfit_gndf_i    ; Best-fit GANDALF model
                    eml_model[i,fit_indx] = eml_model_i ; Best-fit emission-line-only model
                endif

;               ; Best fitting template (does not include LOSVD effects or reddening)
;               best_template[i,fit_indx] = interpol(best_template_i, tpl_wave, obj_wave_lim)
;
;               ; Best fitting galaxy continuum (does not include reddening)
;               best_template_losvd_conv[i,fit_indx] = interpol(best_template_losvd_i, tpl_wave, $
;                                                               obj_wave_lim)

                continue                                        ; Go to next spectrum
            endif
;           endif else begin
            
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

;           ; Best fitting template (does not include LOSVD effects or reddening)
;           best_template[i,fit_indx] = interpol(best_template_i, tpl_wave, wavelength_output)
;           
;           ; Best fitting galaxy continuum (does not include reddening)
;           best_template_losvd_conv[i,fit_indx] = interpol(best_template_losvd_i, tpl_wave, $
;                                                           wavelength_output)

;           endelse
        endfor
END


