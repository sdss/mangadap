;+
; NAME:
;       MDAP_DEFINE_OUTPUT
;
; PURPOSE:
;       Define the output variables to read.  Meant to be called before
;       MDAP_READ_OUTPUT to ensure the variables that are to be read have been
;       previously defined, such that n_element() tests return non-zero values.
;       Input variables are set to '0', IF THEY DON'T ALREADY EXIST!
;
; CALLING SEQUENCE:
;       MDAP_DEFINE_OUTPUT, header=header, dx=dx, dy=dy, w_range_sn=w_range_sn, xpos=xpos, $
;                           ypos=ypos, signal=signal, noise=noise, bin_type=bin_type, $
;                           bin_par=bin_par, threshold_ston_bin=threshold_ston_bin, $
;                           bin_indx=bin_indx, bin_weights=bin_weights, wave=wave, sres=sres, $
;                           bin_flux=bin_flux, bin_ivar=bin_ivar, bin_mask=bin_mask, xbin=xbin, $
;                           ybin=ybin, bin_area=bin_area, bin_ston=bin_ston, bin_n=bin_n, $
;                           bin_flag=bin_flag, w_range_analysis=w_range_analysis, $
;                           threshold_ston_analysis=threshold_ston_analysis, $
;                           analysis_par=analysis_par, weights_ppxf=weights_ppxf, $
;                           add_poly_coeff_ppxf=add_pol_coeff_ppxf, $
;                           mult_poly_coeff_ppxf=mult_pol_coeff_ppxf, $
;                           stellar_kinematics_fit=stellar_kinematics_fit, $
;                           stellar_kinematics_err=stellar_kinematics_err, $
;                           chi2_ppxf=chi2_ppxf, obj_fit_mask_ppxf=obj_fit_mask_ppxf, $
;                           bestfit_ppxf=bestfit_ppxf, weights_gndf=weights_gndf, $
;                           mult_poly_coeff_gndf=mult_poly_coeff_gndf, $
;                           emission_line_kinematics_avg=emission_line_kinematics_avg, $
;                           emission_line_kinematics_aer=emission_line_kinematics_aer, $
;                           chi2_gndf=chi2_gndf, $
;                           emission_line_kinematics_ind=emission_line_kinematics_ind, $
;                           emission_line_kinematics_ier=emission_line_kinematics_ier, $
;                           emission_line_omitted=emission_line_omitted, $
;                           emission_line_intens=emission_line_intens, $
;                           emission_line_interr=emission_line_interr, $
;                           emission_line_fluxes=emission_line_fluxes, $
;                           emission_line_flxerr=emission_line_flxerr, $
;                           emission_line_EWidth=emission_line_EWidth, $
;                           emission_line_EW_err=emission_line_EW_err, $
;                           reddening_val=reddening_val, reddening_err=reddening_err, $
;                           obj_fit_mask_gndf=obj_fit_mask_gndf, bestfit_gndf=bestfit_gndf, $
;                           eml_model=eml_model, optimal_template=optimal_template, $
;                           losvd_optimal_template=losvd_optimal_template
;
; OPTIONAL INPUTS:
;       All optional outputs are set to string values ('0') simply to make sure
;       they are defined.  Using MDAP_READ_OUTPUT, these vectors will be:
;
;       header fits HDU 
;               Header keywords to include with primary header of file.
;       
;       dx double
;               Spaxel size in X, written to header.  Ignored if header not provided.
;
;       dy double
;               Spaxel size in Y, written to header.  Ignored if header not provided.
;
;       w_range_sn dblarr[2]
;               Wavelength range used to calculate the signal and noise per DRP
;               spectrum, written to header.  Ignored if header not provided.
;       
;       xpos dblarr[ndrp]
;               Fiducial X position of every DRP spectrum.  Written to 'DRPS'
;               extension.
;
;       ypos dblarr[ndrp]
;               Fiducial Y position of every DRP spectrum.  Written to 'DRPS'
;               extension.
;
;       signal dblarr[ndrp]
;               Mean signal per pixel in every DRP spectrum.  Written to 'DRPS'
;               extension.
;
;       noise dblarr[ndrp]
;               Mean noise per pixel in every DRP spectrum.  Written to 'DRPS'
;               extension.
;
;       bin_type string
;               Type of binning algorithm used.  Written to header.
;
;       bin_par double
;               Single binning parameter.  Written to header.
;
;       threshold_ston_bin double
;               Threshold for inclusion of a DRP spectrum in any bin. Written to
;               header.
;
;       bin_indx intarr[ndrp]
;               Index (0...B-1) of the bin which contains each DRP spectrum.
;               Written to 'DRPS' extension.
;
;       bin_weights dblarr[ndrp]
;               Weight of each DRP spectrum in its respective bin.  Written to
;               'DRPS' extension.
;
;       wave dblarr[C]
;               Wavelength coordinate of each pixel in the binned spectra.
;               Written as a 1D image to the 'WAVE' extension.
;               
;       sres dblarr[C]
;               Spectral resolution of each pixel in the binned spectra.
;               Written as a 1D image to the 'WAVE' extension.
;
;       bin_flux dblarr[B][C]
;
;               Flux in each channel C of each binned spectrum B.  Written as a
;               2D image to the 'FLUX' extension.
;
;       bin_ivar dblarr[B][C]
;               Inverse variances of each channel C of each binned spectrum B.
;               Written as a 2D image to the 'IVAR' extension.
;
;       bin_mask dblarr[B][C]
;               Pixel mask of each channel C of each binned spectrum B.  Written
;               as a 2D image to the 'MASK' extension.
;
;       xbin dblarr[B]
;               Luminosity-weighted X position of the bin.  Written as a column
;               in 'BINS' extension.
;
;       ybin dblarr[B]
;               Luminosity-weighted Y position of the bin.  Written as a column
;               in 'BINS' extension.
;
;       bin_area dblarr[B]
;               Area of each bin B.  Written as a column in 'BINS' extension.
;
;       bin_ston dblarr[B]
;               S/N of each bin B.  Written as a column in 'BINS' extension.
;
;       bin_n intarr[B]
;               Number of spectra in each bin B.  Written as a column in 'BINS'
;               extension.
;
;       bin_flag intarr[B]
;               Analysis flag for each of the B binned spectra.  Written as a
;               colum in 'BINS' extension.
;
;       w_range_analysis dblarr[2]
;               The nominal wavelength range used in the analysis.  This may be
;               further limited to accommodate some requirements of PPXF and
;               GANDALF.  The true fitted pixels are available from
;               obj_fit_mask_ppxf and obj_fit_mask_gndf.
;
;       obj_fit_mask_ppxf dblarr[B][C]
;               Bad pixel mask used in PPXF fit to all B spectra.  0 for a
;               fitted pixel, 1 for a masked one.
;
;       weights_ppxf dblarr[B][T]
;               Weights of the template spectra in the PPXF fit for each of the
;               B spectra.
;
;       add_poly_coeff_ppxf dblarr[B][]
;               Coefficients of additive polynomials used in the PPXF fit for
;               all B spectra.
;
;       mult_poly_coeff_ppxf dblarr[B][]
;               Coefficients of multiplicative polynomials used in the PPXF fit
;               for all B spectra.
;
;       bestfit_ppxf dblarr[B][C]
;               Best-fit results from PPXF for all B spectra.
;
;       chi2_ppxf dblarr[B]
;               Chi-square per degree of freedom from the PPXF fits to all B
;               spectra.
;
;       stellar_kinematics_fit dblarr[B][]
;       stellar_kinematics_err dblarr[B][]
;               Stellar kinematics (V, sig, [h3, h4, h5, h6]) and their errors
;               from the PPXF fit to all B spectra.
;
;       extra_inputs strarr[]
;               Extra parameters set for the spectral fitting.  TODO: Change
;               this to a structure!
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
; BUGS:
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;       28 Oct 2014: (KBW) Original Implementation
;-
;------------------------------------------------------------------------------

;-------------------------------------------------------------------------------
; Define all the supplied variables
PRO MDAP_DEFINE_OUTPUT, $
                header=header, dx=dx, dy=dy, w_range_sn=w_range_sn, xpos=xpos, ypos=ypos, $
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
                eml_model=eml_model, abs_par=abs_par, abs_line_key=abs_line_key, $
                abs_line_indx_omitted=abs_line_indx_omitted, abs_line_indx_val=abs_line_indx_val, $
                abs_line_indx_err=abs_line_indx_err, abs_line_indx_otpl=abs_line_indx_otpl, $
                abs_line_indx_botpl=abs_line_indx_botpl, si_bin_wave=si_bin_wave, $
                si_bin_flux=si_bin_flux, si_bin_ivar=si_bin_ivar, si_bin_mask=si_bin_mask, $
                si_optimal_template=si_optimal_template, $
                si_broad_optimal_template=si_broad_optimal_template

        if n_elements(header) eq 0 then header = '0'
        if n_elements(dx) eq 0 then dx = '0'
        if n_elements(dy) eq 0 then dy = '0'
        if n_elements(w_range_sn) eq 0 then w_range_sn = '0'
        if n_elements(xpos) eq 0 then xpos = '0'
        if n_elements(ypos) eq 0 then ypos = '0'
        if n_elements(signal) eq 0 then signal = '0'
        if n_elements(noise) eq 0 then noise = '0'
        if n_elements(bin_type) eq 0 then bin_type = '0'
        if n_elements(bin_par) eq 0 then bin_par = '0'
        if n_elements(threshold_ston_bin) eq 0 then threshold_ston_bin = '0'
        if n_elements(bin_indx) eq 0 then bin_indx = '0'
        if n_elements(bin_weights) eq 0 then bin_weights = '0'
        if n_elements(wave) eq 0 then wave = '0'
        if n_elements(sres) eq 0 then sres = '0'
        if n_elements(bin_flux) eq 0 then bin_flux = '0'
        if n_elements(bin_ivar) eq 0 then bin_ivar = '0'
        if n_elements(bin_mask) eq 0 then bin_mask = '0'
        if n_elements(xbin) eq 0 then xbin = '0'
        if n_elements(ybin) eq 0 then ybin = '0'
        if n_elements(bin_area) eq 0 then bin_area = '0'
        if n_elements(bin_ston) eq 0 then bin_ston = '0'
        if n_elements(bin_n) eq 0 then bin_n = '0'
        if n_elements(bin_flag) eq 0 then bin_flag = '0'
        if n_elements(w_range_analysis) eq 0 then w_range_analysis = '0'
        if n_elements(threshold_ston_analysis) eq 0 then threshold_ston_analysis = '0'
        if n_elements(tpl_library_key) eq 0 then tpl_library_key = '0'
        if n_elements(ems_line_key) eq 0 then ems_line_key = '0'
        if n_elements(analysis_par) eq 0 then analysis_par = '0'
        if n_elements(weights_ppxf) eq 0 then weights_ppxf = '0'
        if n_elements(add_poly_coeff_ppxf) eq 0 then add_pol_coeff_ppxf = '0'
        if n_elements(mult_poly_coeff_ppxf) eq 0 then mult_pol_coeff_ppxf = '0'
        if n_elements(stellar_kinematics_fit) eq 0 then stellar_kinematics_fit = '0'
        if n_elements(stellar_kinematics_err) eq 0 then stellar_kinematics_err = '0'
        if n_elements(chi2_ppxf) eq 0 then chi2_ppxf = '0'
        if n_elements(obj_fit_mask_ppxf) eq 0 then obj_fit_mask_ppxf = '0'
        if n_elements(bestfit_ppxf) eq 0 then bestfit_ppxf = '0'
        if n_elements(weights_gndf) eq 0 then weights_gndf = '0'
        if n_elements(mult_poly_coeff_gndf) eq 0 then mult_poly_coeff_gndf = '0'
        if n_elements(emission_line_kinematics_avg) eq 0 then emission_line_kinematics_avg = '0'
        if n_elements(emission_line_kinematics_aer) eq 0 then emission_line_kinematics_aer = '0'
        if n_elements(chi2_gndf) eq 0 then chi2_gndf = '0'
        if n_elements(emission_line_kinematics_ind) eq 0 then emission_line_kinematics_ind = '0'
        if n_elements(emission_line_kinematics_ier) eq 0 then emission_line_kinematics_ier = '0'
        if n_elements(emission_line_omitted) eq 0 then emission_line_omitted = '0'
        if n_elements(emission_line_intens) eq 0 then emission_line_intens = '0'
        if n_elements(emission_line_interr) eq 0 then emission_line_interr = '0'
        if n_elements(emission_line_fluxes) eq 0 then emission_line_fluxes = '0'
        if n_elements(emission_line_flxerr) eq 0 then emission_line_flxerr = '0'
        if n_elements(emission_line_EWidth) eq 0 then emission_line_EWidth = '0'
        if n_elements(emission_line_EW_err) eq 0 then emission_line_EW_err = '0'
        if n_elements(reddening_val) eq 0 then reddening_val = '0'
        if n_elements(reddening_err) eq 0 then reddening_err = '0'
        if n_elements(obj_fit_mask_gndf) eq 0 then obj_fit_mask_gndf = '0'
        if n_elements(bestfit_gndf) eq 0 then bestfit_gndf = '0'
        if n_elements(eml_model) eq 0 then eml_model = '0'
        if n_elements(abs_par) eq 0 then abs_par = '0'
        if n_elements(abs_line_key) eq 0 then abs_line_key = '0'
        if n_elements(abs_line_indx_omitted) eq 0 then abs_line_indx_omitted = '0'
        if n_elements(abs_line_indx_val) eq 0 then abs_line_indx_val = '0'
        if n_elements(abs_line_indx_err) eq 0 then abs_line_indx_err = '0'
        if n_elements(abs_line_indx_otpl) eq 0 then abs_line_indx_otpl = '0'
        if n_elements(abs_line_indx_botpl) eq 0 then abs_line_indx_botpl = '0'
        if n_elements(si_bin_wave) eq 0 then si_bin_wave = '0'
        if n_elements(si_bin_flux) eq 0 then si_bin_flux = '0'
        if n_elements(si_bin_ivar) eq 0 then si_bin_ivar = '0'
        if n_elements(si_bin_mask) eq 0 then si_bin_mask = '0'
        if n_elements(si_optimal_template) eq 0 then si_optimal_template = '0'
        if n_elements(si_broad_optimal_template) eq 0 then si_broad_optimal_template  = '0'
END


