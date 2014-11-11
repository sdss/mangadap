;+
; NAME:
;	MDAP_READ_OUTPUT
;
; PURPOSE:
;	Complement to MDAP_WRITE_OUTPUT, which reads DAP-produced output into
;	the provided data arrays.  See that procedure for the data model.
;
; CALLING SEQUENCE:
;
;	MDAP_READ_OUTPUT, file, header=header, dx=dx, dy=dy, w_range_sn=w_range_sn, xpos=xpos, $
;			  ypos=ypos, signal=signal, noise=noise, bin_type=bin_type, $
;			  bin_par=bin_par, threshold_ston_bin=threshold_ston_bin, $
;			  bin_indx=bin_indx, bin_weights=bin_weights, wave=wave, sres=sres, $
;			  bin_flux=bin_flux, bin_ivar=bin_ivar, bin_mask=bin_mask, xbin=xbin, $
;			  ybin=ybin, bin_area=bin_area, bin_ston=bin_ston, bin_n=bin_n, $
;			  bin_flag=bin_flag, w_range_analysis=w_range_analysis, $
;			  threshold_ston_analysis=threshold_ston_analysis, $
;			  analysis_par=analysis_par, weights_ppxf=weights_ppxf, $
;			  add_poly_coeff_ppxf=add_pol_coeff_ppxf, $
;			  mult_poly_coeff_ppxf=mult_pol_coeff_ppxf, $
;			  stellar_kinematics_fit=stellar_kinematics_fit, $
;			  stellar_kinematics_err=stellar_kinematics_err, chi2_ppxf=chi2_ppxf, $
;			  obj_fit_mask_ppxf=obj_fit_mask_ppxf, bestfit_ppxf=bestfit_ppxf, $
;			  weights_gndf=weights_gndf, mult_poly_coeff_gndf=mult_poly_coeff_gndf, $
;			  emission_line_kinematics_avg=emission_line_kinematics_avg, $
;			  emission_line_kinematics_aer=emission_line_kinematics_aer, $
;			  chi2_gndf=chi2_gndf, $
;			  emission_line_kinematics_ind=emission_line_kinematics_ind, $
;			  emission_line_kinematics_ier=emission_line_kinematics_ier, $
;			  emission_line_omitted=emission_line_omitted, $
;			  emission_line_intens=emission_line_intens, $
;			  emission_line_interr=emission_line_interr, $
;			  emission_line_fluxes=emission_line_fluxes, $
;			  emission_line_flxerr=emission_line_flxerr, $
;			  emission_line_EWidth=emission_line_EWidth, $
;			  emission_line_EW_err=emission_line_EW_err, reddening_val=reddening_val, $
;			  reddening_err=reddening_err, obj_fit_mask_gndf=obj_fit_mask_gndf, $
;			  bestfit_gndf=bestfit_gndf, eml_model=eml_model, $
;			  optimal_template=optimal_template, $
;			  losvd_optimal_template=losvd_optimal_template
;
; INPUTS:
;	file string
;		Name of the fits file of a DAP-produced fits file.
;
; OPTIONAL INPUTS:
;
;	IF DEFINED, the following inputs will be replaced with data read from
;	the DAP-produced fits file.
;
;	header fits HDU 
;		Header keywords to include with primary header of file.
;	
;	dx double
;		Spaxel size in X, written to header.  Ignored if header not provided.
;
;	dy double
;		Spaxel size in Y, written to header.  Ignored if header not provided.
;
;	w_range_sn dblarr[2]
;		Wavelength range used to calculate the signal and noise per DRP
;		spectrum, written to header.  Ignored if header not provided.
;	
;	xpos dblarr[ndrp]
;		Fiducial X position of every DRP spectrum.  Written to 'DRPS'
;		extension.
;
;	ypos dblarr[ndrp]
;		Fiducial Y position of every DRP spectrum.  Written to 'DRPS'
;		extension.
;
;	signal dblarr[ndrp]
;		Mean signal per pixel in every DRP spectrum.  Written to 'DRPS'
;		extension.
;
;	noise dblarr[ndrp]
;		Mean noise per pixel in every DRP spectrum.  Written to 'DRPS'
;		extension.
;
;	bin_type string
;		Type of binning algorithm used.  Written to header.
;
;	bin_par double
;		Single binning parameter.  Written to header.
;
;	threshold_ston_bin double
;		Threshold for inclusion of a DRP spectrum in any bin. Written to
;		header.
;
;	bin_indx intarr[ndrp]
;		Index (0...B-1) of the bin which contains each DRP spectrum.
;		Written to 'DRPS' extension.
;
;	bin_weights dblarr[ndrp]
;		Weight of each DRP spectrum in its respective bin.  Written to
;		'DRPS' extension.
;
;	wave dblarr[C]
;		Wavelength coordinate of each pixel in the binned spectra.
;		Written as a 1D image to the 'WAVE' extension.
;		
;	sres dblarr[C]
;		Spectral resolution of each pixel in the binned spectra.
;		Written as a 1D image to the 'WAVE' extension.
;
;	bin_flux dblarr[B][C]
;
;		Flux in each channel C of each binned spectrum B.  Written as a
;		2D image to the 'FLUX' extension.
;
;	bin_ivar dblarr[B][C]
;		Inverse variances of each channel C of each binned spectrum B.
;		Written as a 2D image to the 'IVAR' extension.
;
;	bin_mask dblarr[B][C]
;		Pixel mask of each channel C of each binned spectrum B.  Written
;		as a 2D image to the 'MASK' extension.
;
;	xbin dblarr[B]
;		Luminosity-weighted X position of the bin.  Written as a column
;		in 'BINS' extension.
;
;	ybin dblarr[B]
;		Luminosity-weighted Y position of the bin.  Written as a column
;		in 'BINS' extension.
;
;	bin_area dblarr[B]
;		Area of each bin B.  Written as a column in 'BINS' extension.
;
;	bin_ston dblarr[B]
;		S/N of each bin B.  Written as a column in 'BINS' extension.
;
;	bin_n intarr[B]
;		Number of spectra in each bin B.  Written as a column in 'BINS'
;		extension.
;
;	bin_flag intarr[B]
;		Analysis flag for each of the B binned spectra.  Written as a
;		colum in 'BINS' extension.
;
;	w_range_analysis dblarr[2]
;		The nominal wavelength range used in the analysis.  This may be
;		further limited to accommodate some requirements of PPXF and
;		GANDALF.  The true fitted pixels are available from
;		obj_fit_mask_ppxf and obj_fit_mask_gndf.
;
;	obj_fit_mask_ppxf dblarr[B][C]
;		Bad pixel mask used in PPXF fit to all B spectra.  0 for a
;		fitted pixel, 1 for a masked one.
;
;	weights_ppxf dblarr[B][T]
;		Weights of the template spectra in the PPXF fit for each of the
;		B spectra.
;
;	add_poly_coeff_ppxf dblarr[B][]
;		Coefficients of additive polynomials used in the PPXF fit for
;		all B spectra.
;
;	mult_poly_coeff_ppxf dblarr[B][]
;		Coefficients of multiplicative polynomials used in the PPXF fit
;		for all B spectra.
;
;	bestfit_ppxf dblarr[B][C]
;		Best-fit results from PPXF for all B spectra.
;
;	chi2_ppxf dblarr[B]
;		Chi-square per degree of freedom from the PPXF fits to all B
;		spectra.
;
;	stellar_kinematics_fit dblarr[B][]
;	stellar_kinematics_err dblarr[B][]
;		Stellar kinematics (V, sig, [h3, h4, h5, h6]) and their errors
;		from the PPXF fit to all B spectra.
;
;	extra_inputs strarr[]
;		Extra parameters set for the spectral fitting.  TODO: Change
;		this to a structure!
;
; OPTIONAL KEYWORDS:
;
; OUTPUT:
;
; OPTIONAL OUTPUT:
;
; COMMENTS:
;
; EXAMPLES:
;
; TODO:
;	- MAKE SURE THINGS ARE ABSTRACTED SUCH THAT I'M READING EXTENSIONS BY
;	  THEIR NAME, NOT THEIR NUMBER!
;
;	- Include:?
;	eml_par	EmissionLine[]
;		An array of EmissionLine structures used to define the emission
;		lines to be fit by GANDALF.
;
;
; BUGS:
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;	28 Oct 2014: (KBW) Original Implementation
;-
;------------------------------------------------------------------------------

; The file that declares the functions used to define the names of the columns
; for the binary tables.
@mdap_set_output_file_cols

;-------------------------------------------------------------------------------
; Determine if the header needs to be read
FUNCTION MDAP_READ_NEED_HEADER, $
		header, dx, dy, w_range_sn, bin_type, bin_par, threshold_ston_bin, $
		w_range_analysis, threshold_ston_analysis, analysis_par, tpl_library_key, $
		ems_line_key, abs_line_key

	if n_elements(header)     ne 0 then              return, 1
	if n_elements(dx)         ne 0 then              return, 1
	if n_elements(dy)         ne 0 then              return, 1
	if n_elements(w_range_sn) ne 0 then              return, 1
	if n_elements(bin_type)   ne 0 then              return, 1
	if n_elements(bin_par)    ne 0 then              return, 1
	if n_elements(threshold_ston_bin) ne 0 then      return, 1
	if n_elements(w_range_analysis) ne 0 then        return, 1
	if n_elements(threshold_ston_analysis) ne 0 then return, 1
	if n_elements(analysis_par) ne 0 then            return, 1
	if n_elements(tpl_library_key) ne 0 then         return, 1
	if n_elements(ems_line_key) ne 0 then            return, 1
	if n_elements(abs_line_key) ne 0 then            return, 1

	return, 0
END

;-------------------------------------------------------------------------------
; Read the header data
PRO MDAP_READ_SET_HEADER_DATA, $
		header, dx, dy, w_range_sn, bin_type, bin_par, threshold_ston_bin, $
		w_range_analysis, threshold_ston_analysis, analysis_par, tpl_library_key, $
		ems_line_key, abs_line_key

	dx = SXPAR(header, 'SPAXDX', /silent)
	dy = SXPAR(header, 'SPAXDY', /silent)
	w_range_sn = [ SXPAR(header, 'SNWAVE1', /silent),SXPAR(header, 'SNWAVE2', /silent) ]
	bin_type = SXPAR(header, 'BINTYPE', /silent)
	bin_par = SXPAR(header, 'BINPAR', /silent)
	threshold_ston_bin = SXPAR(header, 'MINSNBIN', /silent)
	w_range_analysis = [ SXPAR(header, 'FITWAVE1', /silent),SXPAR(header, 'FITWAVE2', /silent) ]
	threshold_ston_analysis = SXPAR(header, 'MINSNFIT', /silent)

	analysis_par = MDAP_DEFINE_ANALYSIS_PAR()
	analysis_par.moments = SXPAR(header, 'MOMENTS', /silent)
	analysis_par.degree = SXPAR(header, 'DEGREE', /silent)
	analysis_par.mdegree = SXPAR(header, 'MDEGREE', /silent)
	analysis_par.reddening_order = SXPAR(header, 'REDORD', /silent)
	if analysis_par.reddening_order ne 0 then begin
	    reddening = SXPAR(header, 'REDORD', /silent)
	    analysis_par.reddening = dblarr(2)
	    for i=0,analysis_par.reddening_order-1 do $
		analysis_par.reddening[i] = reddening[i]
	endif

	tpl_library_key = SXPAR(header, 'TPLKEY', /silent)
	ems_line_key = SXPAR(header, 'EMLKEY', /silent)
	abs_line_key = SXPAR(header, 'SPIKEY', /silent)

END

;-------------------------------------------------------------------------------
; Check if data from the DRPS extension is desired
FUNCTION MDAP_READ_WANT_DRPS, $
		xpos, ypos, signal, noise, bin_indx, bin_weights

	if n_elements(xpos)        ne 0 then return, 1
	if n_elements(ypos)        ne 0 then return, 1
	if n_elements(signal)      ne 0 then return, 1
	if n_elements(noise)       ne 0 then return, 1
	if n_elements(bin_indx)    ne 0 then return, 1
	if n_elements(bin_weights) ne 0 then return, 1

	return, 0
END

; Read the full DRPS binary table (more efficient than reading one column at a
; time?)
PRO MDAP_READ_DRPS, $
		file, xpos, ypos, signal, noise, bin_indx, bin_weights

	cols = MDAP_SET_DRPS_COLS()

	fxbopen, unit, file, 'DRPS'
	fxbreadm, unit, cols, xpos, ypos, signal, noise, bin_indx, bin_weights
	fxbclose, unit
	free_lun, unit
END

;-------------------------------------------------------------------------------
; Check if data from the BINS extension is desired
FUNCTION MDAP_READ_WANT_BINS, $
		xbin, ybin, bin_area, bin_ston, bin_n, bin_flag

	if n_elements(xpos)        ne 0 then return, 1
	if n_elements(ypos)        ne 0 then return, 1
	if n_elements(signal)      ne 0 then return, 1
	if n_elements(noise)       ne 0 then return, 1
	if n_elements(bin_indx)    ne 0 then return, 1
	if n_elements(bin_weights) ne 0 then return, 1

	return, 0
END

; Read the full BINS binary table (more efficient than reading one column at a
; time?)
PRO MDAP_READ_BINS, $
		file, xbin, ybin, bin_area, bin_ston, bin_n, bin_flag

	cols = MDAP_SET_BINS_COLS()

	fxbopen, unit, file, 'BINS'
	fxbreadm, unit, cols, xbin, ybin, bin_area, bin_ston, bin_n, bin_flag
	fxbclose, unit
	free_lun, unit
END

; Read an image extension
FUNCTION MDAP_READ_IMG, $
		file, exten

	errmsg = ''					; Suppress errormessages
	unit=fxposit(file, exten, errmsg=errmsg)	; Find extension and open
	if unit eq -1 then $				; Check there was no error
	    message, errmsg
	
	data=READFITS(unit)				; Read the data
	free_lun, unit					; Close the file and free the LUN
	return, data					; Return the data
END

;-------------------------------------------------------------------------------
; Check if data from the STFIT extension is desired
FUNCTION MDAP_READ_WANT_STFIT, $
		weights_ppxf, add_poly_coeff_ppxf, mult_poly_coeff_ppxf, stellar_kinematics_fit, $
		stellar_kinematics_err, chi2_ppxf

	if n_elements(weights_ppxf)           ne 0 then return, 1
	if n_elements(add_poly_coeff_ppxf)    ne 0 then return, 1
	if n_elements(mult_poly_coeff_ppxf)   ne 0 then return, 1
	if n_elements(stellar_kinematics_fit) ne 0 then return, 1
	if n_elements(stellar_kinematics_err) ne 0 then return, 1
	if n_elements(chi2_ppxf)              ne 0 then return, 1

	return, 0
END

; Read the full STFIT binary table (more efficient than reading one column at a
; time?)
PRO MDAP_READ_STFIT, $
		file, weights_ppxf, add_poly_coeff_ppxf, mult_poly_coeff_ppxf, $
		stellar_kinematics_fit, stellar_kinematics_err, chi2_ppxf

	cols = MDAP_SET_STFIT_COLS()

	fxbopen, unit, file, 'STFIT'
	fxbreadm, unit, cols, weights_ppxf, add_poly_coeff_ppxf, mult_poly_coeff_ppxf, $
			      stellar_kinematics_fit, stellar_kinematics_err, chi2_ppxf
	fxbclose, unit
	free_lun, unit

	; TODO: 2D data needs to be transposed!

	weights_ppxf = transpose(temporary(weights_ppxf))
	add_poly_coeff_ppxf = transpose(temporary(add_poly_coeff_ppxf))
	mult_poly_coeff_ppxf = transpose(temporary(mult_poly_coeff_ppxf))
	stellar_kinematics_fit = transpose(temporary(stellar_kinematics_fit))
	stellar_kinematics_err = transpose(temporary(stellar_kinematics_err))

END

;-------------------------------------------------------------------------------
; Check if data from the SGFIT extension is desired
FUNCTION MDAP_READ_WANT_SGFIT, $ 
		weights_gndf, mult_poly_coeff_gndf, emission_line_kinematics_avg, $
		emission_line_kinematics_aer, chi2_gndf, emission_line_kinematics_ind, $
		emission_line_kinematics_ier, emission_line_omitted, emission_line_intens, $
		emission_line_interr, emission_line_fluxes, emission_line_flxerr, $
		emission_line_EWidth, emission_line_EW_err, reddening_val, reddening_err

	if n_elements(weights_gndf)                 ne 0 then return, 1
	if n_elements(mult_poly_coeff_gndf)         ne 0 then return, 1
	if n_elements(emission_line_kinematics_avg) ne 0 then return, 1
	if n_elements(emission_line_kinematics_aer) ne 0 then return, 1
	if n_elements(chi2_gndf)                    ne 0 then return, 1
	if n_elements(emission_line_kinematics_ind) ne 0 then return, 1
	if n_elements(emission_line_kinematics_ier) ne 0 then return, 1
	if n_elements(emission_line_omitted)        ne 0 then return, 1
	if n_elements(emission_line_intens)    ne 0 then return, 1
	if n_elements(emission_line_interr)    ne 0 then return, 1
	if n_elements(emission_line_fluxes)    ne 0 then return, 1
	if n_elements(emission_line_flxerr)    ne 0 then return, 1
	if n_elements(emission_line_EWidth)    ne 0 then return, 1
	if n_elements(emission_line_EW_err)    ne 0 then return, 1
	if n_elements(reddening_val)                ne 0 then return, 1
	if n_elements(reddening_err)                ne 0 then return, 1

	return, 0
END

; Read the full SGFIT binary table (more efficient than reading one column at a
; time?)
PRO MDAP_READ_SGFIT, $
		file, weights_gndf, mult_poly_coeff_gndf, emission_line_kinematics_avg, $
		emission_line_kinematics_aer, chi2_gndf, emission_line_kinematics_ind, $
		emission_line_kinematics_ier, emission_line_omitted, emission_line_intens, $
		emission_line_interr, emission_line_fluxes, emission_line_flxerr, $
		emission_line_EWidth, emission_line_EW_err, reddening_val, reddening_err

	cols = MDAP_SET_SGFIT_COLS()

	fxbopen, unit, file, 'SGFIT'
	fxbreadm, unit, cols, weights_gndf, mult_poly_coeff_gndf, emission_line_kinematics_avg, $
			      emission_line_kinematics_aer, chi2_gndf, reddening_val, $
			      reddening_err, emission_line_omitted, emission_line_intens, $
			      emission_line_interr, emission_line_kinematics_ind, $
			      emission_line_kinematics_ier, emission_line_fluxes, $
			      emission_line_flxerr, emission_line_EWidth, emission_line_EW_err

	fxbclose, unit
	free_lun, unit

	; TODO: 2D data needs to be transposed!

	weights_gndf = transpose(temporary(weights_gndf))
	mult_poly_coeff_gndf = transpose(temporary(mult_poly_coeff_gndf))
	emission_line_kinematics_avg = transpose(temporary(emission_line_kinematics_avg))
	emission_line_kinematics_aer = transpose(temporary(emission_line_kinematics_aer))
	reddening_val = transpose(temporary(reddening_val))
	reddening_err = transpose(temporary(reddening_err))
	emission_line_intens = transpose(temporary(emission_line_intens))
	emission_line_interr = transpose(temporary(emission_line_interr))
	emission_line_kinematics_ind = transpose(temporary(emission_line_kinematics_ind), [2,0,1])
	emission_line_kinematics_ier = transpose(temporary(emission_line_kinematics_ier), [2,0,1])
	emission_line_fluxes = transpose(temporary(emission_line_fluxes))
	emission_line_flxerr = transpose(temporary(emission_line_flxerr))
	emission_line_EWidth = transpose(temporary(emission_line_EWidth))
	emission_line_EW_err = transpose(temporary(emission_line_EW_err))

END

;-------------------------------------------------------------------------------
; Check if data from the SIPAR extension is desired
FUNCTION MDAP_READ_WANT_SIPAR, $
		abs_par

	if n_elements(abs_par) ne 0 then return, 1

	return, 0
END

; Read the full SIPAR binary table (more efficient than reading one column at a
; time?)
PRO MDAP_READ_SIPAR, $
		file, abs_par

	cols = MDAP_SET_SIPAR_COLS()

	fxbopen, unit, file, 'SIPAR'
	fxbreadm, unit, cols, name, passband, blue_cont, red_cont, unit
	fxbclose, unit
	free_lun, unit

	nindx = n_elements(name)
	id = indgen(nindx)

	abs_par = MDAP_ASSIGN_ABSORPTION_LINE_PARAMETERS(id, name, reform(passband[0,*]), $
							 reform(passband[1,*]), $
							 reform(blue_cont[0,*]), $
							 reform(blue_cont[1,*]), $
							 reform(red_cont[0,*]), $
							 reform(red_cont[1,*]), unit)

END

;-------------------------------------------------------------------------------
; Check if data from the SINDX extension is desired
FUNCTION MDAP_READ_WANT_SINDX, $ 
		abs_line_indx_omitted, abs_line_indx_val, abs_line_indx_err, abs_line_indx_otpl, $
		abs_line_indx_botpl

	if n_elements(abs_line_indx_omitted) ne 0 then return, 1
	if n_elements(abs_line_indx_val)     ne 0 then return, 1
	if n_elements(abs_line_indx_err)     ne 0 then return, 1
	if n_elements(abs_line_indx_otpl)    ne 0 then return, 1
	if n_elements(abs_line_indx_botpl)   ne 0 then return, 1

	return, 0
END

; Read the full SINDX binary table (more efficient than reading one column at a
; time?)
PRO MDAP_READ_SINDX, $
		file, abs_line_indx_omitted, abs_line_indx_val, abs_line_indx_err, $
		abs_line_indx_otpl, abs_line_indx_botpl

	cols = MDAP_SET_SINDX_COLS()

	fxbopen, unit, file, 'SINDX'
	fxbreadm, unit, cols, abs_line_indx_omitted, abs_line_indx_val, abs_line_indx_err, $
			      abs_line_indx_otpl, abs_line_indx_botpl
	fxbclose, unit
	free_lun, unit

	; TODO: 2D data needs to be transposed!

	abs_line_indx_omitted = transpose(temporary(abs_line_indx_omitted))
	abs_line_indx_val = transpose(temporary(abs_line_indx_val))
	abs_line_indx_err = transpose(temporary(abs_line_indx_err))
	abs_line_indx_otpl = transpose(temporary(abs_line_indx_otpl))
	abs_line_indx_botpl = transpose(temporary(abs_line_indx_botpl))

END



;-------------------------------------------------------------------------------
; Order should not really be important here
; TODO: Add information such that eml_par can be read?
PRO MDAP_READ_OUTPUT, $
		file, header=header, dx=dx, dy=dy, w_range_sn=w_range_sn, xpos=xpos, ypos=ypos, $
		signal=signal, noise=noise, bin_type=bin_type, bin_par=bin_par, $
		threshold_ston_bin=threshold_ston_bin, bin_indx=bin_indx, bin_weights=bin_weights, $
		wave=wave, sres=sres, bin_flux=bin_flux, bin_ivar=bin_ivar, bin_mask=bin_mask, $
		xbin=xbin, ybin=ybin, bin_area=bin_area, bin_ston=bin_ston, bin_n=bin_n, $
		bin_flag=bin_flag, w_range_analysis=w_range_analysis, $
		threshold_ston_analysis=threshold_ston_analysis, tpl_library_key=tpl_library_key, $
		ems_line_key=ems_line_key, analysis_par=analysis_par, weights_ppxf=weights_ppxf, $
		add_poly_coeff_ppxf=add_pol_coeff_ppxf, mult_poly_coeff_ppxf=mult_pol_coeff_ppxf, $
		stellar_kinematics_fit=stellar_kinematics_fit, $
		stellar_kinematics_err=stellar_kinematics_err, chi2_ppxf=chi2_ppxf, $
		obj_fit_mask_ppxf=obj_fit_mask_ppxf, bestfit_ppxf=bestfit_ppxf, $
		weights_gndf=weights_gndf, mult_poly_coeff_gndf=mult_poly_coeff_gndf, $
		emission_line_kinematics_avg=emission_line_kinematics_avg, $
		emission_line_kinematics_aer=emission_line_kinematics_aer, $
		chi2_gndf=chi2_gndf, emission_line_kinematics_ind=emission_line_kinematics_ind, $
		emission_line_kinematics_ier=emission_line_kinematics_ier, $
		emission_line_omitted=emission_line_omitted, $
		emission_line_intens=emission_line_intens, $
		emission_line_interr=emission_line_interr, $
		emission_line_fluxes=emission_line_fluxes, $
		emission_line_flxerr=emission_line_flxerr, $
		emission_line_EWidth=emission_line_EWidth, $
		emission_line_EW_err=emission_line_EW_err, $
		reddening_val=reddening_val, reddening_err=reddening_err, $
		obj_fit_mask_gndf=obj_fit_mask_gndf, bestfit_gndf=bestfit_gndf, $
		eml_model=eml_model, optimal_template=optimal_template, $
		losvd_optimal_template=losvd_optimal_template, abs_par=abs_par, $
		abs_line_key=abs_line_key, abs_line_indx_omitted=abs_line_indx_omitted, $
		abs_line_indx_val=abs_line_indx_val, abs_line_indx_err=abs_line_indx_err, $
		abs_line_indx_otpl=abs_line_indx_otpl, abs_line_indx_botpl=abs_line_indx_botpl

	if file_test(file) eq 0 then $				; Make sure file exists
	    message, 'File does not exist!'

	; Check that the file looks like a DAP file
	if MDAP_CHECK_OUTPUT_FILE(file) eq 0 then $
	    message, 'File does not have DAP-like (hard-coded) properties!'

	; If the header is needed/requested, read it
	if MDAP_READ_NEED_HEADER(header, dx, dy, w_range_sn, bin_type, bin_par, $
				 threshold_ston_bin, w_range_analysis, threshold_ston_analysis, $
				 analysis_par, tpl_library_key, ems_line_key, $
				 abs_line_key) ne 0 then begin

	    header=headfits(file, exten=0)
	    MDAP_READ_SET_HEADER_DATA, header, dx, dy, w_range_sn, bin_type, bin_par, $
				       threshold_ston_bin, w_range_analysis, $
				       threshold_ston_analysis, analysis_par, tpl_library_key, $
				       ems_line_key, abs_line_key
	endif

	; Read the DRPS extension if any of its vectors are requested
	if MDAP_READ_WANT_DRPS(xpos, ypos, signal, noise, bin_indx, bin_weights) eq 1 then $
	    MDAP_READ_DRPS, file, xpos, ypos, signal, noise, bin_indx, bin_weights

	; Read the BINS extension if any of its vectors are requested
	if MDAP_READ_WANT_BINS(xbin, ybin, bin_area, bin_ston, bin_n, bin_flag) eq 1 then $
	    MDAP_READ_BINS, file, xbin, ybin, bin_area, bin_ston, bin_n, bin_flag

	; Read the wavelength vector if requested
	if n_elements(wave) ne 0 then $
	    wave=MDAP_READ_IMG(file, 'WAVE')

	; Read the spectral resolution vector if requested
	if n_elements(sres) ne 0 then $
	    sres=MDAP_READ_IMG(file, 'SRES')

	; Read the binned spectra, if requested
	if n_elements(bin_flux) ne 0 then $
	    bin_flux=MDAP_READ_IMG(file, 'FLUX')

	; Read the inverse variances, if requested
	if n_elements(bin_ivar) ne 0 then $
	    bin_ivar=MDAP_READ_IMG(file, 'IVAR')

	; Read the bad pixel mask, if requested
	if n_elements(bin_mask) ne 0 then $
	    bin_mask=MDAP_READ_IMG(file, 'MASK')

	; Read the emission-line parameters
	; TODO: to allow this need to write more information to the DAP output...
;	if n_elements(eml_par) ne 0 then $
;	    eml_par = MDAP_READ_EMISSION_LINE_PARAMETERS(file)

	; Read the STFIT extension if any of its vectors are requested
	if MDAP_READ_WANT_STFIT(weights_ppxf, add_poly_coeff_ppxf, mult_poly_coeff_ppxf, $
				stellar_kinematics_fit, stellar_kinematics_err, $
				chi2_ppxf) eq 1 then begin

	    MDAP_READ_STFIT, file, weights_ppxf, add_poly_coeff_ppxf, mult_poly_coeff_ppxf, $
			     stellar_kinematics_fit, stellar_kinematics_err, chi2_ppxf
	endif

	; Read the fit_pixel mask resulting from PPXF only, if requested
	if n_elements(obj_fit_mask_ppxf) ne 0 then $
	    obj_fit_mask_ppxf = MDAP_READ_IMG(file, 'SMSK')

	; Read the best-fitting PPXF-only model, if requested
	if n_elements(bestfit_ppxf) ne 0 then $
	    bestfit_ppxf = MDAP_READ_IMG(file, 'SMOD')

	; Read the SGFIT extension if any of its vectors are requested
	if MDAP_READ_WANT_SGFIT(weights_gndf, mult_poly_coeff_gndf, emission_line_kinematics_avg, $
				emission_line_kinematics_aer, chi2_gndf, $
				emission_line_kinematics_ind, emission_line_kinematics_ier, $
				emission_line_omitted, emission_line_intens, emission_line_interr, $
				emission_line_fluxes, emission_line_flxerr, emission_line_EWidth, $
				emission_line_EW_err, reddening_val, reddening_err) eq 1 then begin

	    MDAP_READ_SGFIT, file, weights_gndf, mult_poly_coeff_gndf, $
			     emission_line_kinematics_avg, emission_line_kinematics_aer, $
			     chi2_gndf, emission_line_kinematics_ind, $
			     emission_line_kinematics_ier, emission_line_omitted, $
			     emission_line_intens, emission_line_interr, emission_line_fluxes, $
			     emission_line_flxerr, emission_line_EWidth, emission_line_EW_err, $
			     reddening_val, reddening_err

	endif

	; Read the fit_pixel mask resulting from GANDALF, if requested
	if n_elements(obj_fit_mask_gndf) ne 0 then $
	    obj_fit_mask_gndf = MDAP_READ_IMG(file, 'SGMSK')

	; Read the best-fitting GANDALF model, if requested
	if n_elements(bestfit_gndf) ne 0 then $
	    bestfit_gndf = MDAP_READ_IMG(file, 'SGMOD')

	; Read the best-fitting emission-line-only model, if requested
	if n_elements(eml_model) ne 0 then $
	    eml_model = MDAP_READ_IMG(file, 'ELMOD')

	; Read the optimal template, if requested
	if n_elements(optimal_template) ne 0 then $
	    optimal_tempalte = MDAP_READ_IMG(file, 'OTPL')

	; Read the losvd-broadened optimal template, if requested
	if n_elements(losvd_optimal_template) ne 0 then $
	    losvd_optimal_tempalte = MDAP_READ_IMG(file, 'BOTPL')

	; Read the SIPAR extension if abs_par is requested
	if MDAP_READ_WANT_SIPAR(abs_par) eq 1 then $
	    MDAP_READ_SIPAR, file, abs_par

	; Read the SINDX extension if any of its vectors are requested
	if MDAP_READ_WANT_SINDX(abs_line_indx_omitted, abs_line_indx_val, abs_line_indx_err, $
				abs_line_indx_otpl, abs_line_indx_botpl) eq 1 then begin

	    MDAP_READ_SINDX, file, abs_line_indx_omitted, abs_line_indx_val, abs_line_indx_err, $
			     abs_line_indx_otpl, abs_line_indx_botpl
	endif

END



