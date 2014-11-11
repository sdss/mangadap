;+
; NAME:
;	MDAP_ABSORPTION_LINE_INDEX_MEASUREMENTS
;
; PURPOSE:
;	Measure absorption-line indices.
;
;	TODO: Better this description.
;
; CALLING SEQUENCE:
;	MDAP_ABSORPTION_LINE_INDEX_MEASUREMENTS, tpl_lib_fits, weights, stellar_kinematics, $
;						 abs_par_file, abs_vacuum_wave, obj_wave, obj_sres,$
;						 obj_flux, obj_ivar, obj_mask, bestfit, eml_model, $
;						 abs_line_indx_omitted, abs_line_indx, $
;						 abs_line_indx_err, abs_line_indx_otpl, $
;						 abs_line_indx_botpl, moments=moments, $
;						 remove_outliers=remove_outliers, version=version, $
;						 output_file_root=output_file_root, /oversample, $
;						 /rms
;
; INPUTS:
;	tpl_lib_fits string
;		File name for the fits file with the resolution- and
;		sampling-matched template library.  This should have been
;		previously created using MDAP_CREATE_DAP_TEMPLATE_LIBRARY.
;
;	weights dblarr[B][T]
;		Weights that when applied to the template library provides the
;		optimal template for the fit to each of the B object spectra.
;
;	stellar_kinematics dblarr[B][M]
;		Best-fitting stellar kinematics with M moments for each of the B
;		object spectra.
;
;	abs_par_file string
;		Name of the file containing the indices to measure.  File has 9
;		columns:
;
;	# ID Name        pass_0    pass_1   blcnt_0   blcnt_1   rdcnt_0   rdcnt_1 units
;	#                   ang       ang       ang       ang       ang       ang
;  	   0 D4000      0000.00   0000.00   3750.00   3950.00   4050.00   4250.00   ang 
;	   1 CaII0p39   3899.50   4003.50   3806.50   3833.80   4020.70   4052.40   ang
;	   2 HDeltaA    4083.50   4122.25   4041.60   4079.75   4128.50   4161.00   ang
;	   3 HDeltaF    4091.00   4112.25   4057.25   4088.50   4114.75   4137.25   ang
;
;		1 - ID: ID number for the index (integer)
;
;		2 - Name: Name for the index (string)
;
;		3,4 - pass_0,pass_1: Starting and ending wavelength over which
;		to get the equivalent width (TODO: CHECK). (float)
;		
;		5,6 - blcnt_0,blcnt_1: Starting and ending wavelength over which
;		to measure the continuum level blueward of the line index. (float)
;		
;		7,8 - rdcnt_0,rdcnt_1: Starting and ending wavelength over which
;		to measure the continuum level redward of the line index. (float)
;		
;		9 - units: Units to return for the line index
;
;	obj_wave dblarr[C]
;		Wavelength at each of C spectral channels of the input spectra.
; 
;	obj_sres dblarr[C]
;		Spectral resolution at each of the C spectral channels of the
;		input spectra.
; 
; 	obj_flux dblarr[B][C]
;		Flux in the B observed spectra with C spectral channels each.
;
;	obj_ivar dblarr[B][C]
;		Inverse variance in the flux of the observed spectra.
;	
;	obj_mask dblarr[B][C]
;		Bad pixel mask (0-good, 1-bad) for the object spectra.
;
;	bestfit dblarr[B][C]
;		Best-fitting spectrum to the object spectra.
;
;	eml_model dblarr[B][C]
;		Best-fitting emission-line model for the object spectra (can be
;		all zeros if no emission lines were fit or exist)
;	
; OPTIONAL INPUTS:
;	moments integer
;		Number of kinematic moments.  By default this is determined
;		based on the length of the second dimension of
;		stellar_kinematics.
;
;	output_file_root string
;		Root name for PS file output.
;
;	remove_outliers double
;		Replace pixels that deviate more than remove_outliers*rms with
;		the best-fitting spectrum.
;
; OPTIONAL KEYWORDS:
;	/plot
;		Create the output PS plots.  CURRENTLY REMOVED!
;
;	/oversample
;		Force the convolution to oversample the template and the LOSVD
;		when creating the broadened template spectrum. See MDAP_PPXF.
;
;	/rms
;		Use the RMS difference between the best-fit model and the 
;
;	/dbg
;		Only fit the set of indices to the first spectrum as a debugging
;		test.
;
; OUTPUT:
;	abs_line_omitted intarr[B][I]
;		Flag that the index was (1) or was not (0) omitted because it
;		fell outside the spectral range of the observations.
;
;	abs_line_indx dblarr[B][I]
;		Absorption-line-index measurements for the I absorption lines in
;		the B *object* spectra, corrected for intrinsic broadening.
;
;	abs_line_indx_err dblarr[B][I]
;		Error in each of the index measurements.
;
;	abs_line_indx_otpl dblarr[B][I]
;		Absorption-line-index measurements based on the B
;		optimal_template spectra.
;		
;	abs_line_indx_botpl dblarr[B][I]
;		Absorption-line-index measurements based on the B
;		losvd_optimal_template spectra.
;
; OPTIONAL OUTPUT:
;	version string
;		Module version.  If requested, the module is not executed and
;		only version flag is returned.
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; TODO:
;	- Allow for a tag that will ignore the MDAP_REVERT_PIXEL_KINEMATICS
;	  call(s).
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;	31 Oct 2014: (KBW) Copied from L. Coccato's version (PDAP v0_8)
;-
;------------------------------------------------------------------------------

PRO MDAP_ABSORPTION_LINE_INDEX_MEASUREMENTS, $
		tpl_lib_fits, weights, stellar_kinematics, abs_par, obj_wave, obj_sres, obj_flux, $
		obj_ivar, obj_mask, bestfit, eml_model, abs_line_indx_omitted, abs_line_indx, $
		abs_line_indx_err, abs_line_indx_otpl, abs_line_indx_botpl, moments=moments, $
		output_file_root=output_file_root, remove_outliers=remove_outliers, $
		version=version, oversample=oversample, rms=rms, dbg=dbg

	; TODO: Plotting option has been removed.  There has to be a better way to create QA plots.
	;, plot=plot

	; TODO: Do some checks of the input arrays?

	version_module = '0.9'			; Version number
	if n_elements(version) ne 0 then begin	; If version is defined
	    version = version_module		; ... set it to the module version
	    return				; ... and return without doing anything
	endif

	;-----------------------------------------------------------------------
	; SETUP ----------------------------------------------------------------
	; Get the template spectra that have previously been matched to the
	; spectral resolution of the absorption-line index system. See MANGA_DAP.
	MDAP_READ_RESAMPLED_TEMPLATES, tpl_lib_fits, tpl_wave, tpl_flux, tpl_ivar, tpl_mask, $
				       tpl_sres, tpl_soff

	; Manipulate the object data:
	; 0. Do NOT alter the input spectra, so copy over to some working arrays
	obj_wave_indx = obj_wave
	obj_flux_indx = obj_flux
	obj_ivar_indx = obj_ivar
	obj_mask_indx = obj_mask
	obj_sres_indx = obj_sres

	sz = size(obj_flux)
	nspec = sz[1]					; Number of spectra
	for i=0,nspec-1 do begin

	    ;  1. Determine the error spectrum
	    obj_mask_ = reform(obj_mask[i,*])
	    if keyword_set(rms) then begin			; Use the fit residual
		MDAP_NOISE_FROM_RESIDUAL, reform(obj_flux[i,*]), obj_mask_, reform(bestfit[i,*]), $
					  obj_sige
	    endif else $					; Use the inverse variance vector
		MDAP_NOISE_FROM_IVAR, reform(obj_ivar[i,*]), obj_mask_, obj_sige

	    ; TODO: Don't like this switching back and forth between ivar and sige

	    obj_mask_indx[i,*] = obj_mask_
	    obj_ivar_indx[i,*] = (obj_sige)^(-2)

	    ;  2. Replace the outliers with the model

	    ; TODO: This is different from what Lodo did.  He based the outliers
	    ;       on a measurement of sigma within the full spectral range of
	    ;       the index.

	    obj_flux_ = reform(obj_flux[i,*])
	    if n_elements(remove_outliers) ne 0 then $
		MDAP_REPLACE_OUTLIERS, obj_flux_, reform(bestfit[i,*]), obj_sige, remove_outliers

	    obj_flux_indx[i,*] = obj_flux_
	endfor

	; 3. Remove the emission lines from the object spectra
	obj_flux_indx = obj_flux_indx - eml_model

	; 4. Match the spectral resolution of the object spectra to the
	; absorption-line index system.

	; TODO: This may perform some interpolation of tpl_sres onto the range
	;       of obj_wave (?)
	MDAP_MATCH_SPECTRAL_RESOLUTION, obj_flux_indx, obj_ivar_indx, obj_mask_indx, $
					obj_wave_indx, obj_sres_indx, tpl_wave, tpl_sres, $
					obj_soff_indx, /no_offset

	; Get the number of moments
	if n_elements(moments) eq 0 then $
	    moments = (size(stellar_kinematics))[2]	; Number of elements in stellar kinematics

	; Velocity scale of the object data
	velScale = MDAP_VELOCITY_SCALE(obj_wave, /log10)

;	if keyword_set(plot) and n_elements(output_file_root) eq 0 then begin
;	    print, 'No output root name provide for PS plots.  Using ./mdap_abs_line'
;	    output_file_root = './mdap_abs_line'
;	endif

	c=299792.48d				; Speed of light in km/s

	; Determine the minimum and maximum wavelength available to the index
	; measurements; force the index to be common to both the object and
	; template wavelengths
	min_wave = min(obj_wave_indx)		; Grab the minimum
	if min_wave lt min(tpl_wave) then $	; Adjust it if necessary
	    min_wave = min(tpl_wave)
	max_wave = max(obj_wave_indx)		; Grab the maximum
	if max_wave gt max(tpl_wave) then $	; Adjust it if necessary
	    max_wave = max(tpl_wave)

	; Allocate the output arrays
	nabs = n_elements(abs_par)
	abs_line_indx_omitted = intarr(nspec, nabs)	; Flag that index was omitted
	abs_line_indx = dblarr(nspec, nabs)		; Index measurements (object spectra)
	abs_line_indx_err = dblarr(nspec, nabs)		; ... error
	abs_line_indx_otpl = dblarr(nspec, nabs)	; Index measurements (optimal template)
	abs_line_indx_botpl = dblarr(nspec, nabs) 	; Index measurements (broadened op template)

	; TODO: These are dummy vectors!
	otpl_ivar = make_array((size(tpl_flux))[2], /double, value=1.0)
	otpl_mask = dblarr((size(tpl_flux))[2])

	;-----------------------------------------------------------------------

;	if keyword_set(plot) then $
;	    set_plot,'ps'

	;-----------------------------------------------------------------------
	; Measure all spectral indices for all spectra -------------------------
	if keyword_set(dbg) then $				; Only analyze the first spectrum
	    nspec = 1
	for i=0,nspec-1 do begin

	    ; Get the redshifted absorption-line windows
	    redshift = stellar_kinematics[i,0]/c
	    print, 'redshift: ', redshift
	    abs_par_z = abs_par
	    for i=0,nabs-1 do begin
		abs_par_z.passband = abs_par[i].passband*(1.+redshift)
		abs_par_z.blue_cont = abs_par[i].blue_cont*(1.+redshift)
		abs_par_z.red_cont = abs_par[i].red_cont*(1.+redshift)
	    endfor

	    optimal_template = reform(weights[i,*] # tpl_flux)	; Create the optimal template

	    ; Get the pixel-based kinematics
	    sol = reform(stellar_kinematics[i,0:moments-1])
	    vel = sol[0]
	    vel_err = 1.0
	    MDAP_REVERT_PIXEL_KINEMATICS, vel, vel_err
	    sol[0] = vel

	    ; Get the broadened optimal template
	    losvd_optimal_template = MDAP_GET_BROADENED_TEMPLATE(optimal_template, sol, velScale, $
								 moments=moments, $
								 oversample=oversample)

;	    ; Initialize the plot file
;	    if keyword_set(plot) then begin
;		device, filename=output_file_root+'indx_'+mdap_stc(i+1,/integer)+'.ps', /color, $
;			xsize=15, ysize=12
;	    endif

	    ; Measure each index
	    for j=0,nabs-1 do begin
		print, 'index:', j+1

		; Get the REDSHIFTED passbands for getting the index
		; measurements of the object data
		abs_par_z = abs_par[j]
		abs_par_z.passband = abs_par_z.passband*(redshift+1.0d)
		abs_par_z.blue_cont = abs_par_z.blue_cont*(redshift+1.0d)
		abs_par_z.red_cont = abs_par_z.red_cont*(redshift+1.0d)

		print, 'def:'
		print, abs_par_z.passband
		print, abs_par_z.blue_cont
		print, abs_par_z.red_cont

		; Check if the index can be measured
		if abs_par[j].red_cont[1] lt min_wave or $
		   abs_par[j].blue_cont[0] gt max_wave or $
		   abs_par_z.red_cont[1] lt min_wave or $
		   abs_par_z.blue_cont[0] gt max_wave then begin
		    abs_line_indx_omitted[i,j] = 1	; Flag the index as omitted
		    print, 'OMITTED!'
		    continue				; Index not in wavelength range, go to next
		endif

		; Measure the index based on the object data
		MDAP_MEASURE_SPECTRAL_INDEX, obj_wave_indx, reform(obj_flux_indx[i,*]), $
					     reform(obj_ivar_indx[i,*]), $
					     reform(obj_mask_indx[i,*]), $
					     abs_par_z, obj_equiv_width, $
					     obj_equiv_width_err, obj_index_mag, obj_index_mag_err
		; Measure the index based on the optimal template
		MDAP_MEASURE_SPECTRAL_INDEX, tpl_wave, optimal_template, otpl_ivar, otpl_mask, $
					     abs_par[j], otpl_equiv_width, otpl_equiv_width_err, $
					     otpl_index_mag, otpl_index_mag_err
		; Measure the index based on the broadened template
		MDAP_MEASURE_SPECTRAL_INDEX, tpl_wave, losvd_optimal_template, otpl_ivar, $
					     otpl_mask, abs_par[j], botpl_equiv_width, $
					     botpl_equiv_width_err, botpl_index_mag, $
					     botpl_index_mag_err

		; Correct the index measured using the difference between the
		; measurements performed on the optimal template and its
		; broadened version.  The correction depends on the desired
		; units.  This also saves the results.
		if abs_par[j].units eq 'mag' then begin
		    losvd_correction = otpl_index_mag - botpl_index_mag
		    abs_line_indx[i,j] = obj_index_mag + losvd_correction
		    abs_line_indx_otpl[i,j] = otpl_index_mag
		    abs_line_indx_botpl[i,j] = botpl_index_mag
		    abs_line_indx_err[i,j] = obj_index_mag_err
		endif

		if abs_par[j].units eq 'ang' then begin
		    losvd_correction = otpl_equiv_width / botpl_equiv_width
		    abs_line_indx[i,j] = obj_equiv_width / losvd_correction
		    abs_line_indx_otpl[i,j] = otpl_equiv_width
		    abs_line_indx_botpl[i,j] = botpl_equiv_width
		    abs_line_indx_err[i,j] = obj_equiv_width_err / abs(losvd_correction)
		endif

	    endfor ; End loop over indices

;	    if keyword_set(plot) then $
;		device, /close

	endfor ; End loop over spectra
	;-----------------------------------------------------------------------

;	if keyword_set(plot) then $
;	    set_plot, 'x'

END

