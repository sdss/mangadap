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
;                           ypos=ypos, fraction_good=fraction_good, min_eq_max=min_eq_max, $
;                           signal=signal, noise=noise, bin_par=bin_par, $
;                           threshold_ston_bin=threshold_ston_bin, bin_vreg=bin_vreg, $
;                           bin_indx=bin_indx, bin_weights=bin_weights, wave=wave, sres=sres, $
;                           bin_flux=bin_flux, bin_ivar=bin_ivar, bin_mask=bin_mask, $
;                           xbin_rlow=xbin_rlow, ybin_rupp=ybin_rupp, rbin=rbin, $
;                           bin_area=bin_area, bin_ston=bin_ston, bin_n=bin_n, bin_flag=bin_flag, $
;                           w_range_analysis=w_range_analysis, $ 
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
;                           emission_line_sinst=emission_line_sinst, $
;                           emission_line_omitted=emission_line_omitted, $
;                           emission_line_intens=emission_line_intens, $
;                           emission_line_interr=emission_line_interr, $
;                           emission_line_fluxes=emission_line_fluxes, $
;                           emission_line_flxerr=emission_line_flxerr, $
;                           emission_line_EWidth=emission_line_EWidth, $
;                           emission_line_EW_err=emission_line_EW_err, $
;                           reddening_val=reddening_val, reddening_err=reddening_err, $
;                           obj_fit_mask_gndf=obj_fit_mask_gndf, bestfit_gndf=bestfit_gndf, $
;                           eml_model=eml_model, nonpar_eml_omitted=nonpar_eml_omitted, $
;                           nonpar_eml_flux_raw=nonpar_eml_flux_raw, $
;                           nonpar_eml_flux_rerr=nonpar_eml_flux_rerr, $
;                           nonpar_eml_mom1_raw=nonpar_eml_mom1_raw, $
;                           nonpar_eml_mom1_rerr=nonpar_eml_mom1_rerr, $
;                           nonpar_eml_mom2_raw=nonpar_eml_mom2_raw, $
;                           nonpar_eml_mom2_rerr=nonpar_eml_mom2_rerr, $
;                           nonpar_eml_blue_flux=nonpar_eml_blue_flux, $
;                           nonpar_eml_blue_ferr=nonpar_eml_blue_ferr, $
;                           nonpar_eml_red_flux=nonpar_eml_red_flux, $
;                           nonpar_eml_red_ferr=nonpar_eml_red_ferr, $
;                           nonpar_eml_blue_fcen=nonpar_eml_blue_fcen, $
;                           nonpar_eml_blue_cont=nonpar_eml_blue_cont, $
;                           nonpar_eml_blue_cerr=nonpar_eml_blue_cerr, $
;                           nonpar_eml_red_fcen=nonpar_eml_red_fcen, $
;                           nonpar_eml_red_cont=nonpar_eml_red_cont, $
;                           nonpar_eml_red_cerr=nonpar_eml_red_cerr, $
;                           nonpar_eml_flux_corr=nonpar_eml_flux_corr, $
;                           nonpar_eml_flux_cerr=nonpar_eml_flux_cerr, $
;                           nonpar_eml_mom1_corr=nonpar_eml_mom1_corr, $
;                           nonpar_eml_mom1_cerr=nonpar_eml_mom1_cerr, $
;                           nonpar_eml_mom2_corr=nonpar_eml_mom2_corr, $
;                           nonpar_eml_mom2_cerr=nonpar_eml_mom2_cerr, $
;                           elo_ew_kinematics_avg=elo_ew_kinematics_avg, $
;                           elo_ew_kinematics_aer=elo_ew_kinematics_aer, $
;                           elo_ew_kinematics_ase=elo_ew_kinematics_ase, $
;                           elo_ew_kinematics_n=elo_ew_kinematics_n, $
;                           elo_ew_window=elo_ew_window, elo_ew_baseline=elo_ew_baseline, $
;                           elo_ew_base_err=elo_ew_base_err, $
;                           elo_ew_kinematics_ind=elo_ew_kinematics_ind, $
;                           elo_ew_kinematics_ier=elo_ew_kinematics_ier, $
;                           elo_ew_sinst=elo_ew_sinst, elo_ew_omitted=elo_ew_omitted, $
;                           elo_ew_intens=elo_ew_intens, elo_ew_interr=elo_ew_interr, $
;                           elo_ew_fluxes=elo_ew_fluxes, elo_ew_flxerr=elo_ew_flxerr, $
;                           elo_ew_EWidth=elo_ew_EWidth, elo_ew_EW_err=elo_ew_EW_err, $
;                           elo_ew_eml_model=elo_ew_eml_model, $
;                           elo_fb_kinematics_avg=elo_fb_kinematics_avg, $
;                           elo_fb_kinematics_aer=elo_fb_kinematics_aer, $
;                           elo_fb_kinematics_ase=elo_fb_kinematics_ase, $
;                           elo_fb_kinematics_n=elo_fb_kinematics_n, $
;                           elo_fb_window=elo_fb_window, elo_fb_baseline=elo_fb_baseline, $
;                           elo_fb_base_err=elo_fb_base_err, $
;                           elo_fb_kinematics_ind=elo_fb_kinematics_ind, $
;                           elo_fb_kinematics_ier=elo_fb_kinematics_ier, $
;                           elo_fb_sinst=elo_fb_sinst, elo_fb_omitted=elo_fb_omitted, $
;                           elo_fb_intens=elo_fb_intens, elo_fb_interr=elo_fb_interr, $
;                           elo_fb_fluxes=elo_fb_fluxes, elo_fb_flxerr=elo_fb_flxerr, $
;                           elo_fb_EWidth=elo_fb_EWidth, elo_fb_EW_err=elo_fb_EW_err, $
;                           elo_fb_eml_model=elo_fb_eml_model, abs_par=abs_par, $
;                           abs_line_key=abs_line_key, $
;                           spectral_index_omitted=spectral_index_omitted, $
;                           spectral_index_blue_cont=spectral_index_blue_cont, $
;                           spectral_index_blue_cerr=spectral_index_blue_cerr, $
;                           spectral_index_red_cont=spectral_index_red_cont, $
;                           spectral_index_red_cerr=spectral_index_red_cerr, $
;                           spectral_index_raw=spectral_index_raw, $
;                           spectral_index_rerr=spectral_index_rerr, $
;                           spectral_index_rperr=spectral_index_rperr, $
;                           spectral_index_otpl_blue_cont=spectral_index_otpl_blue_cont, $
;                           spectral_index_otpl_red_cont=spectral_index_otpl_red_cont, $
;                           spectral_index_otpl_index=spectral_index_otpl_index, $
;                           spectral_index_botpl_blue_cont=spectral_index_botpl_blue_cont, $
;                           spectral_index_botpl_red_cont=spectral_index_botpl_red_cont, $
;                           spectral_index_botpl_index=spectral_index_botpl_index, $
;                           spectral_index_corr=spectral_index_corr, $
;                           spectral_index_cerr=spectral_index_cerr, $
;                           spectral_index_cperr=spectral_index_cperr, si_bin_wave=si_bin_wave, $
;                           si_bin_flux=si_bin_flux, si_bin_ivar=si_bin_ivar, $
;                           si_bin_mask=si_bin_mask, si_otpl_flux=si_otpl_flux, $
;                           si_otpl_mask=si_otpl_mask, si_botpl_flux=si_botpl_flux, $
;                           si_botpl_mask=si_botpl_mask
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
;       fraction_good dblarr[ndrp]
;               Fraction of good pixels in each of the DRP spectra.
;               Written to 'DRPS' extension.
;
;       min_eq_max intarr[ndrp]
;               Flag (0-false; 1-true) that the minimum and maximum flux
;               values are the same, presumably meaning that the
;               spectrum has no data.  Written to 'DRPS' extension.
;
;       signal dblarr[ndrp]
;               Mean signal per pixel in every DRP spectrum.  Written to 'DRPS'
;               extension.
;
;       noise dblarr[ndrp]
;               Mean noise per pixel in every DRP spectrum.  Written to 'DRPS'
;               extension.
;
;       bin_par BinPar
;               Structure containing all the binning parameters.
;
;       threshold_ston_bin double
;               Threshold for inclusion of a DRP spectrum in any bin. Written to
;               header.
;
;       bin_vreg dblarr[ndrp]
;               Velocity of each of the DRP spectra used to de-redshift
;               the spectra before binning.  Written to 'DRPS'.
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
;       xbin_rlow dblarr[B]
;               For the result of a RADIAL binning, this is the lower
;               limit of the radial bin; otherwise, this is the
;               luminosity-weighted X position of the bin.  Written as a
;               column in 'BINS' extension.  
;
;       ybin_rupp dblarr[B]
;               For the result of a RADIAL binning, this is the upper
;               limit of the radial bin; otherwise, this is the
;               luminosity-weighted Y position of the bin.  Written as a
;               column in 'BINS' extension.  
;
;       rbin dblarr[B]
;               For the result of a RADIAL binning, this is the
;               luminosity weighted radius of the spectra in the bin;
;               otherwise, this is the radius of the luminosity weighted
;               on-sky coordinates.  Written as a column in teh 'BINS'
;               extension.
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
;       obj_fit_mask_ppxf dblarr[B][C]
;               Bad pixel mask used in PPXF fit to all B spectra.  0 for a
;               fitted pixel, 1 for a masked one.
;
;       bestfit_ppxf dblarr[B][C]
;               Best fitting spectrum obtained by PPXF for each of the B
;               spectra.
;
;       weights_gndf dblarr[B][T]
;               Template weights for each spectrum obtained by GANDALF.
;
;       mult_poly_coeff_gndf dblarr[B][MD]
;               Multiplicative order MD legendre polynomial coefficients obtained
;               for each of N spectra by GANDALF.
;
;       emission_line_kinematics_avg dblarr[B][2]
;               The best-fit V and sigma for the emission lines in the B galaxy
;               spectra.
;
;       emission_line_kinematics_aer dblarr[B][2]
;               Estimates of the errors in the best-fit V and sigma for the gas
;               kinematics for each of the N fitted input galaxy spectra.
;
;       chi2_gndf dblarr[B]
;               Chi-square per degree of freedom obtained by the GANDALF fit to
;               each of the B spectra.
;
;       emission_line_kinematics_ind dblarr[B][E][2]
;       emission_line_kinematics_ier dblarr[B][E][2]
;               Kinematics and errors for each fitted emission line.
;
;       emission_line_sinst dblarr[B][E]
;               Instrumental dispersion at the fitted line centers of
;               the emission lines.
;
;       emission_line_omitted intarr[B][E]
;               Flag setting whether or not an emission-line was fit for all E
;               emission lines in the eml_par structure.  0 means the
;               emission-line was fit, 1 means it was not.
;
;       emission_line_intens dblarr[B][E]
;       emission_line_interr dblarr[B][E]
;               Best-fitting emission-line intensities and errors of all fitted
;               emission lines for each of the B galaxy spectra.
;
;       emission_line_fluxes dblarr[B][E]
;       emission_line_flxerr dblarr[B][E]
;               Reddening-corrected integrated fluxes and errors for all fitted
;               emission lines for each of the N input galaxy spectra. 
;
;       emission_line_EWidth dblarr[B][E]
;       emission_line_EW_err dblarr[B][E]
;               Equivalent widths and errors of all fitted emission lines for
;               each of the B input galaxy spectra.  These equivalent widths are
;               computed by taking the ratio of the emission_line_fluxes and the
;               median value of the stellar spectrum within 5 and 10 sigma of
;               the emission line, where 'sigma' is the velocity dispersion.
;
;       reddening_val dblarr[B][2]
;       reddening_err dblarr[B][2]
;               Best-fit values and errors for stellar reddening
;               (reddening_val[*,0]) and gas reddening (reddening_val[*,1]) for
;               all B galaxy spectra.  If the reddening fit is not performed,
;               the output value is set to 0. (reddening_output[*,0:1] = [0,0]).
;               If only the reddening of the stars is fitted, the reddening of
;               the gas is set to 0 (reddening_output[*,1] = 0).
;
;       obj_fit_mask_gndf dblarr[B][C]
;               Bad pixel mask used in the GANDALF fit to all B spectra.  0 for
;               a fitted pixel, 1 for a masked one.
;
;       bestfit_gndf dblarr[B][C]
;               Best fitting spectrum obtained by GANDALF for each of the N
;               spectra.
;
;       eml_model dblarr[B][C]
;               Best-fitting emission-line-only model for each of the N spectra
;               obtained by GANDALF.
;
;
;       nonpar_eml_omitted intarr[N][E]
;               Flag setting whether or not a set of emission-line
;               measurements should be omitted due to errors
;               (0-good;1-bad).
;
;       nonpar_eml_flux_raw dblarr[N][E]
;       nonpar_eml_flux_rerr dblarr[N][E]
;               The integrated flux of the continuum-subtracted spectrum
;               over the defined passband and its error.  
;       
;       nonpar_eml_mom1_raw dblarr[N][E]
;       nonpar_eml_mom1_rerr dblarr[N][E]
;       nonpar_eml_mom2_raw dblarr[N][E]
;       nonpar_eml_mom2_rerr dblarr[N][E]
;               The first and second moments of the velocity over the
;               defined passband and their errors.
;
;       nonpar_eml_blue_flux dblarr[N][E]
;       nonpar_eml_blue_ferr dblarr[N][E]
;       nonpar_eml_red_flux dblarr[N][E]
;       nonpar_eml_red_ferr dblarr[N][E]
;               The integrated flux of the continuum-subtracted spectrum
;               over the defined blue and red sidebands, and their
;               errors.
;
;       nonpar_eml_blue_fcen dblarr[N][E]
;       nonpar_eml_blue_cont dblarr[N][E]
;       nonpar_eml_blue_cerr dblarr[N][E]
;       nonpar_eml_red_fcen dblarr[N][E]
;       nonpar_eml_red_cont dblarr[N][E]
;       nonpar_eml_red_cerr dblarr[N][E]
;               The flux-weighted centers and pseudo-continuum levels
;               and errors in the blue and red sidebands used to define
;               a linear continuum underneath the primary passband for
;               the corrected fluxes and velocity moments.
;
;       nonpar_eml_flux_corr dblarr[N][E]
;       nonpar_eml_flux_cerr dblarr[N][E]
;               The integrated flux of the continuum-subtracted spectrum
;               over the defined passband and its error, including the
;               adjusted continuum level based on the blue and red
;               pseudo-continuum measurements.  
;        
;       nonpar_eml_mom1_corr dblarr[N][E]
;       nonpar_eml_mom1_cerr dblarr[N][E]
;       nonpar_eml_mom2_corr dblarr[N][E]
;       nonpar_eml_mom2_cerr dblarr[N][E]
;               The first and second moments of the velocity over the
;               defined passband and their errors , including the
;               adjusted continuum level based on the blue and red
;               pseudo-continuum measurements.
;
;       elo_ew_kinematics_avg dblarr[B][8]
;       elo_ew_kinematics_aer dblarr[B][8]
;       elo_ew_kinematics_ase dblarr[B][8]
;       elo_ew_kinematics_n intarr[B]
;               The weighted-mean kinematics (avg), propagated error
;               (aer), standard error (ase; based on the weighted
;               standard deviation), and number of used emission lines,
;               for the set of fitted emission lines determined using
;               Enci Wang's code.  There are 8 elements, 2 kinematic
;               measurements (v, sigma) for each of the 4 weighting
;               types (flux, velocity error, both, unweighted); see
;               MDAP_EMISSION_LINE_ONLY_FIT.
;
;       elo_ew_window intarr[N][E][2]
;               Wavelength limits used in the fit of each (set of)
;               emission line(s) for Enci Wang's code.
;
;       elo_ew_baseline intarr[N][E]
;       elo_ew_base_err intarr[N][E]
;               Constant baseline fit beneath each (group of) emission
;               line(s) and its error determined using Enci Wang's code.
;
;       elo_ew_kinematics_ind dblarr[B][O][2]
;       elo_ew_kinematics_ierr dblarr[B][O][2]
;               Kinematics of the individual lines (and errors)
;               determined using Enci Wang's code.
;
;       elo_ew_sinst dblarr[B][O]
;               The instrumental dispersion at the measured centroid of
;               the emission lines returned by Enci Wang's code.
;
;       elo_ew_omitted intarr[B][O]
;               Flag that a line has been omitted from the emission-line
;               analysis using Enci Wang's code.  All fits are actually
;               analyzed; this flag is tripped (set to 1) if fitted flux
;               is not greater than 0 (see
;               MDAP_SAVE_EMISSION_LINE_FIT_ENCI).
;
;       elo_ew_intens dblarr[B][O]
;       elo_ew_interr dblarr[B][O]
;               Measured intensity (and error) of each of the lines fit
;               by Enci Wang's code.  Note this is a post calculation;
;               Enci's code actually uses the Gaussian area (flux) as
;               the fitting parameter.
;
;       elo_ew_fluxes dblarr[B][O]
;       elo_ew_flxerr dblarr[B][O]
;               Measured flux (and error) of each of the lines fit
;               by Enci Wang's code.
;
;       elo_ew_EWidth dblarr[B][O]
;       elo_ew_EW_err dblarr[B][O]
;               Measured equivalent width (and error) of each of the
;               lines fit by Enci Wang's code.  Post-calculated.
;
;       elo_ew_eml_model dblarr[B][C]
;               The best-fitting emission-line-only spectrum determined
;               by Enci Wang's code.
;
;       elo_fb_kinematics_avg dblarr[B][8]
;       elo_fb_kinematics_aer dblarr[B][8]
;       elo_fb_kinematics_ase dblarr[B][8]
;       elo_fb_kinematics_n intarr[B]
;               The weighted-mean kinematics (avg), propagated error
;               (aer), standard error (ase; based on the weighted
;               standard deviation), and number of used emission lines,
;               for the set of fitted emission lines determined using
;               Francesco Belfiore's code.  There are 8 elements, 2
;               kinematic measurements (v, sigma) for each of the 4
;               weighting types (flux, velocity error, both,
;               unweighted); see MDAP_EMISSION_LINE_ONLY_FIT.
;
;       elo_fb_window intarr[N][E][2]
;               Wavelength limits used in the fit of each (set of)
;               emission line(s) for Francesco Belfiore's code.
;
;       elo_fb_baseline intarr[N][E]
;       elo_fb_base_err intarr[N][E]
;               Constant baseline fit beneath each (group of) emission
;               line(s) and its error determined using Francesco
;               Belfiore's code.
;
;       elo_fb_kinematics_ind dblarr[B][O][2]
;       elo_fb_kinematics_ierr dblarr[B][O][2]
;               Kinematics of the individual lines (and errors)
;               determined using Francesco Belfiore's code.
;
;       elo_fb_sinst dblarr[B][O]
;               The instrumental dispersion at the measured centroid of
;               the emission lines returned by Francesco Belfiore's
;               code.
;
;       elo_fb_omitted intarr[B][O]
;               Flag that a line has been omitted from the emission-line
;               analysis using Francesco Belfiore's code.  All fits are
;               actually analyzed; this flag is tripped (set to 1) if
;               fitted amplitude is not greater than 0 (see
;               MDAP_SAVE_EMISSION_LINE_FIT_BELFIORE).
;
;       elo_fb_intens dblarr[B][O]
;       elo_fb_interr dblarr[B][O]
;               Measured intensity (and error) of each of the lines fit
;               by Francesco Belfiore's code.
;
;       elo_fb_fluxes dblarr[B][O]
;       elo_fb_flxerr dblarr[B][O]
;               Measured flux (and error) of each of the lines fit by
;               Francesco Belfiore's code.  Note this is a post
;               calculation; Francesco's code actually uses the
;               intensity as the fitting parameter.
;
;       elo_fb_EWidth dblarr[B][O]
;       elo_fb_EW_err dblarr[B][O]
;               Measured equivalent width (and error) of each of the
;               lines fit by Francesco Belfiore's code.
;               Post-calculated.
;
;       elo_fb_eml_model dblarr[B][C]
;               The best-fitting emission-line-only spectrum determined
;               by Francesco Belfiore's code.
;
;       abs_par SpectralIndex[I]
;               Array of structures with the spectral index parameters.  See
;               MDAP_READ_ABSORPTION_LINE_PARAMETERS.
;
;       abs_line_key string
;               Keyword signifying the spectral index file.
;
;       spectral_index_omitted intarr[B][I]
;               Flag that the index was (1) or was not (0) omitted because it
;               fell outside the spectral range of the observations.
;
;       spectral_index_blue_cont dblarr[B][I]
;       spectral_index_blue_cerr dblarr[B][I]
;       spectral_index_red_cont dblarr[B][I]
;       spectral_index_red_cerr dblarr[B][I]
;               The blue and red pseudo-continuum measurements and their
;               propagated error for each of the B binned spectra and
;               each of the I spectral indices.
;
;       spectral_index_raw dblarr[B][I]
;       spectral_index_rerr dblarr[B][I]
;       spectral_index_rperr dblarr[B][I]
;               Spectral index measured directly from the galaxy spectra
;               following Worthey et al. (1994), the nominal error
;               following Cardiel et al. (1998) (err), and an estimate
;               of error from nominal error propagation (perr).  The
;               index is defined either in angstroms or magnitudes as
;               specified by abs_par.
;
;       spectral_index_otpl_blue_cont dblarr[B][I]
;       spectral_index_otpl_red_cont dblarr[B][I]
;       spectral_index_otpl_index dblarr[B][I]
;               The blue and red pseudo-continuum measurements and the
;               index value for the (unbroadened) optimal template for
;               each of the B binned spectra and each of the I spectral
;               indices.
;
;       spectral_index_botpl_blue_cont dblarr[B][I]
;       spectral_index_botpl_red_cont dblarr[B][I]
;       spectral_index_botpl_index dblarr[B][I]
;               The blue and red pseudo-continuum measurements and the
;               index value for the broadened optimal template for each
;               of the B binned spectra and each of the I spectral
;               indices.
;
;       spectral_index_corr dblarr[B][I]
;       spectral_index_cerr dblarr[B][I]
;       spectral_index_cperr dblarr[B][I]
;               Spectral indices for the galaxy spectra after correcting
;               for Doppler broadening using the correction derived by
;               comparing the index measurements from the broadened and
;               unbroadened versions of the optimal template.
;
;       si_bin_flux dblarr[B][D]
;               Flux of the binned spectra that have had their resolution
;               matched tot the spectral index system.
;
;       si_bin_ivar dblarr[B][D]
;               Inverse variance of the binned spectra that have had their
;               resolution matched tot the spectral index system.
;
;       si_bin_mask dblarr[B][D]
;               Bad pixel mask (0-good, 1-bad) of the binned spectra that have
;               had their resolution matched tot the spectral index system.
;
;       si_otpl_flux dblarr[B][D]
;       si_otpl_mask dblarr[B][D]
;               Optimal template spectrum, resolution-matched to the
;               spectral-index system, and its bad pixel mask.  The bad
;               pixel mask is just a de-redshifted version of the
;               broadened optimal template mask.
;
;       si_botpl_flux dblarr[B][D]
;       si_botpl_mask dblarr[B][D]
;               Optimal template spectrum, resolution-matched to the
;               spectral-index system, and broadened by the best-fitting
;               LOSVD for each of the B spectra, and its bad pixel mask.
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
;       28 Oct 2014: Original implementation by K. Westfall (KBW)
;       08 Dec 2014: (KBW) Accommodate radial binning and velocity
;                          registration
;       12 Dec 2014: (KBW) New format incorporating emission-line only
;                          results
;       09 Jan 2015: (KBW) Include instrumental dispersion for GANDALF fit
;       09 Feb 2015: (KBW) Include fraction of good pixels and min(flux)
;                          == max(flux) flag in DRPS extension.
;       12 Aug 2015: (KBW) Adjusted for new parameters.
;       18 Aug 2015: (KBW) Added non-parametric emission-line data
;       17 Sep 2015: (KBW) Added emission-line window and baseline data
;-
;------------------------------------------------------------------------------


;-------------------------------------------------------------------------------
; Define all the supplied variables
PRO MDAP_DEFINE_OUTPUT, $
                header=header, dx=dx, dy=dy, w_range_sn=w_range_sn, xpos=xpos, ypos=ypos, $
                fraction_good=fraction_good, min_eq_max=min_eq_max, signal=signal, noise=noise, $
                bin_par=bin_par, threshold_ston_bin=threshold_ston_bin, bin_vreg=bin_vreg, $
                bin_indx=bin_indx, bin_weights=bin_weights, wave=wave, sres=sres, $
                bin_flux=bin_flux, bin_ivar=bin_ivar, bin_mask=bin_mask, xbin_rlow=xbin_rlow, $
                ybin_rupp=ybin_rupp, rbin=rbin, bin_area=bin_area, bin_ston=bin_ston, bin_n=bin_n, $
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
                emission_line_sinst=emission_line_sinst, $
                emission_line_omitted=emission_line_omitted, $
                emission_line_intens=emission_line_intens, $
                emission_line_interr=emission_line_interr, $
                emission_line_fluxes=emission_line_fluxes, $
                emission_line_flxerr=emission_line_flxerr, $
                emission_line_EWidth=emission_line_EWidth, $
                emission_line_EW_err=emission_line_EW_err, $
                reddening_val=reddening_val, reddening_err=reddening_err, $
                obj_fit_mask_gndf=obj_fit_mask_gndf, bestfit_gndf=bestfit_gndf, $
                eml_model=eml_model, elo_ew_kinematics_avg=elo_ew_kinematics_avg, $
                elo_ew_kinematics_aer=elo_ew_kinematics_aer, $
                elo_ew_kinematics_ase=elo_ew_kinematics_ase, $
                elo_ew_kinematics_n=elo_ew_kinematics_n, elo_ew_window=elo_ew_window, $
                elo_ew_baseline=elo_ew_baseline, elo_ew_base_err=elo_ew_base_err, $
                elo_ew_kinematics_ind=elo_ew_kinematics_ind, $
                elo_ew_kinematics_ier=elo_ew_kinematics_ier, elo_ew_sinst=elo_ew_sinst, $
                elo_ew_omitted=elo_ew_omitted, elo_ew_intens=elo_ew_intens, $
                elo_ew_interr=elo_ew_interr, elo_ew_fluxes=elo_ew_fluxes, $
                elo_ew_flxerr=elo_ew_flxerr, elo_ew_EWidth=elo_ew_EWidth, $
                elo_ew_EW_err=elo_ew_EW_err, elo_ew_eml_model=elo_ew_eml_model, $
                elo_fb_kinematics_avg=elo_fb_kinematics_avg, $
                elo_fb_kinematics_aer=elo_fb_kinematics_aer, $
                elo_fb_kinematics_ase=elo_fb_kinematics_ase, $
                elo_fb_kinematics_n=elo_fb_kinematics_n, elo_fb_window=elo_fb_window, $
                elo_fb_baseline=elo_fb_baseline, elo_fb_base_err=elo_fb_base_err, $
                elo_fb_kinematics_ind=elo_fb_kinematics_ind, $
                elo_fb_kinematics_ier=elo_fb_kinematics_ier, elo_fb_sinst=elo_fb_sinst, $
                elo_fb_omitted=elo_fb_omitted, elo_fb_intens=elo_fb_intens, $
                elo_fb_interr=elo_fb_interr, elo_fb_fluxes=elo_fb_fluxes, $
                elo_fb_flxerr=elo_fb_flxerr, elo_fb_EWidth=elo_fb_EWidth, $
                elo_fb_EW_err=elo_fb_EW_err, elo_fb_eml_model=elo_fb_eml_model, abs_par=abs_par, $
                abs_line_key=abs_line_key, spectral_index_omitted=spectral_index_omitted, $
                spectral_index_blue_cont=spectral_index_blue_cont, $
                spectral_index_blue_cerr=spectral_index_blue_cerr, $
                spectral_index_red_cont=spectral_index_red_cont, $
                spectral_index_red_cerr=spectral_index_red_cerr, $
                spectral_index_raw=spectral_index_raw, spectral_index_rerr=spectral_index_rerr, $
                spectral_index_rperr=spectral_index_rperr, $
                spectral_index_otpl_blue_cont=spectral_index_otpl_blue_cont, $
                spectral_index_otpl_red_cont=spectral_index_otpl_red_cont, $
                spectral_index_otpl_index=spectral_index_otpl_index, $
                spectral_index_botpl_blue_cont=spectral_index_botpl_blue_cont, $
                spectral_index_botpl_red_cont=spectral_index_botpl_red_cont, $
                spectral_index_botpl_index=spectral_index_botpl_index, $
                spectral_index_corr=spectral_index_corr, spectral_index_cerr=spectral_index_cerr, $
                spectral_index_cperr=spectral_index_cperr, si_bin_wave=si_bin_wave, $
                si_bin_flux=si_bin_flux, si_bin_ivar=si_bin_ivar, si_bin_mask=si_bin_mask, $
                si_otpl_flux=si_otpl_flux, si_otpl_mask=si_otpl_mask, si_botpl_flux=si_botpl_flux, $
                si_botpl_mask=si_botpl_mask

        if n_elements(header) eq 0 then header = '0'
        if n_elements(dx) eq 0 then dx = '0'
        if n_elements(dy) eq 0 then dy = '0'
        if n_elements(w_range_sn) eq 0 then w_range_sn = '0'
        if n_elements(xpos) eq 0 then xpos = '0'
        if n_elements(ypos) eq 0 then ypos = '0'
        if n_elements(fraction_good) eq 0 then fraction_good = '0'
        if n_elements(min_eq_max) eq 0 then min_eq_max = '0'
        if n_elements(signal) eq 0 then signal = '0'
        if n_elements(noise) eq 0 then noise = '0'
        if n_elements(bin_par) eq 0 then bin_par = '0'
        if n_elements(threshold_ston_bin) eq 0 then threshold_ston_bin = '0'
        if n_elements(bin_vreg) eq 0 then bin_vreg = '0'
        if n_elements(bin_indx) eq 0 then bin_indx = '0'
        if n_elements(bin_weights) eq 0 then bin_weights = '0'
        if n_elements(wave) eq 0 then wave = '0'
        if n_elements(sres) eq 0 then sres = '0'
        if n_elements(bin_flux) eq 0 then bin_flux = '0'
        if n_elements(bin_ivar) eq 0 then bin_ivar = '0'
        if n_elements(bin_mask) eq 0 then bin_mask = '0'
        if n_elements(xbin_rlow) eq 0 then xbin_rlow = '0'
        if n_elements(ybin_rupp) eq 0 then ybin_rupp = '0'
        if n_elements(rbin) eq 0 then rbin = '0'
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
        if n_elements(emission_line_sinst) eq 0 then emission_line_sinst = '0'
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
        if n_elements(nonpar_eml_omitted) eq 0 then nonpar_eml_omitted = '0'
        if n_elements(nonpar_eml_flux_raw) eq 0 then nonpar_eml_flux_raw = '0'
        if n_elements(nonpar_eml_flux_rerr) eq 0 then nonpar_eml_flux_rerr = '0'
        if n_elements(nonpar_eml_mom1_raw) eq 0 then nonpar_eml_mom1_raw = '0'
        if n_elements(nonpar_eml_mom1_rerr) eq 0 then nonpar_eml_mom1_rerr = '0'
        if n_elements(nonpar_eml_mom2_raw) eq 0 then nonpar_eml_mom2_raw = '0'
        if n_elements(nonpar_eml_mom2_rerr) eq 0 then nonpar_eml_mom2_rerr = '0'
        if n_elements(nonpar_eml_blue_flux) eq 0 then nonpar_eml_blue_flux = '0'
        if n_elements(nonpar_eml_blue_ferr) eq 0 then nonpar_eml_blue_ferr = '0'
        if n_elements(nonpar_eml_red_flux) eq 0 then nonpar_eml_red_flux = '0'
        if n_elements(nonpar_eml_red_ferr) eq 0 then nonpar_eml_red_ferr = '0'
        if n_elements(nonpar_eml_blue_fcen) eq 0 then nonpar_eml_blue_fcen = '0'
        if n_elements(nonpar_eml_blue_cont) eq 0 then nonpar_eml_blue_cont = '0'
        if n_elements(nonpar_eml_blue_cerr) eq 0 then nonpar_eml_blue_cerr = '0'
        if n_elements(nonpar_eml_red_fcen) eq 0 then nonpar_eml_red_fcen = '0'
        if n_elements(nonpar_eml_red_cont) eq 0 then nonpar_eml_red_cont = '0'
        if n_elements(nonpar_eml_red_cerr) eq 0 then nonpar_eml_red_cerr = '0'
        if n_elements(nonpar_eml_flux_corr) eq 0 then nonpar_eml_flux_corr = '0'
        if n_elements(nonpar_eml_flux_cerr) eq 0 then nonpar_eml_flux_cerr = '0'
        if n_elements(nonpar_eml_mom1_corr) eq 0 then nonpar_eml_mom1_corr = '0'
        if n_elements(nonpar_eml_mom1_cerr) eq 0 then nonpar_eml_mom1_cerr = '0'
        if n_elements(nonpar_eml_mom2_corr) eq 0 then nonpar_eml_mom2_corr = '0'
        if n_elements(nonpar_eml_mom2_cerr) eq 0 then nonpar_eml_mom2_cerr = '0'
        if n_elements(elo_ew_kinematics_avg) eq 0 then elo_ew_kinematics_avg = '0'
        if n_elements(elo_ew_kinematics_aer) eq 0 then elo_ew_kinematics_aer = '0'
        if n_elements(elo_ew_kinematics_ase) eq 0 then elo_ew_kinematics_ase = '0'
        if n_elements(elo_ew_kinematics_n) eq 0 then elo_ew_kinematics_n = '0'
        if n_elements(elo_ew_window) eq 0 then elo_ew_window = '0'
        if n_elements(elo_ew_baseline) eq 0 then elo_ew_baseline = '0'
        if n_elements(elo_ew_base_err) eq 0 then elo_ew_base_err = '0'
        if n_elements(elo_ew_kinematics_ind) eq 0 then elo_ew_kinematics_ind = '0'
        if n_elements(elo_ew_kinematics_ier) eq 0 then elo_ew_kinematics_ier = '0'
        if n_elements(elo_ew_sinst) eq 0 then elo_ew_sinst = '0'
        if n_elements(elo_ew_omitted) eq 0 then elo_ew_omitted = '0'
        if n_elements(elo_ew_intens) eq 0 then elo_ew_intens = '0'
        if n_elements(elo_ew_interr) eq 0 then elo_ew_interr = '0'
        if n_elements(elo_ew_fluxes) eq 0 then elo_ew_fluxes = '0'
        if n_elements(elo_ew_flxerr) eq 0 then elo_ew_flxerr = '0'
        if n_elements(elo_ew_EWidth) eq 0 then elo_ew_EWidth = '0'
        if n_elements(elo_ew_EW_err) eq 0 then elo_ew_EW_err = '0'
        if n_elements(elo_ew_eml_model) eq 0 then elo_ew_eml_model = '0'
        if n_elements(elo_fb_kinematics_avg) eq 0 then elo_fb_kinematics_avg = '0'
        if n_elements(elo_fb_kinematics_aer) eq 0 then elo_fb_kinematics_aer = '0'
        if n_elements(elo_fb_kinematics_ase) eq 0 then elo_fb_kinematics_ase = '0'
        if n_elements(elo_fb_kinematics_n) eq 0 then elo_fb_kinematics_n = '0'
        if n_elements(elo_fb_window) eq 0 then elo_fb_window = '0'
        if n_elements(elo_fb_baseline) eq 0 then elo_fb_baseline = '0'
        if n_elements(elo_fb_base_err) eq 0 then elo_fb_base_err = '0'
        if n_elements(elo_fb_kinematics_ind) eq 0 then elo_fb_kinematics_ind = '0'
        if n_elements(elo_fb_kinematics_ier) eq 0 then elo_fb_kinematics_ier = '0'
        if n_elements(elo_fb_sinst) eq 0 then elo_fb_sinst = '0'
        if n_elements(elo_fb_omitted) eq 0 then elo_fb_omitted = '0'
        if n_elements(elo_fb_intens) eq 0 then elo_fb_intens = '0'
        if n_elements(elo_fb_interr) eq 0 then elo_fb_interr = '0'
        if n_elements(elo_fb_fluxes) eq 0 then elo_fb_fluxes = '0'
        if n_elements(elo_fb_flxerr) eq 0 then elo_fb_flxerr = '0'
        if n_elements(elo_fb_EWidth) eq 0 then elo_fb_EWidth = '0'
        if n_elements(elo_fb_EW_err) eq 0 then elo_fb_EW_err = '0'
        if n_elements(elo_fb_eml_model) eq 0 then elo_fb_eml_model = '0'
        if n_elements(abs_par) eq 0 then abs_par = '0'
        if n_elements(abs_line_key) eq 0 then abs_line_key = '0'
        if n_elements(spectral_index_omitted) eq 0 then spectral_index_omitted = '0'
        if n_elements(spectral_index_blue_cont) eq 0 then spectral_index_blue_cont = '0'
        if n_elements(spectral_index_blue_cerr) eq 0 then spectral_index_blue_cerr = '0'
        if n_elements(spectral_index_red_cont) eq 0 then spectral_index_red_cont = '0'
        if n_elements(spectral_index_red_cerr) eq 0 then spectral_index_red_cerr = '0'
        if n_elements(spectral_index_raw) eq 0 then spectral_index_raw = '0'
        if n_elements(spectral_index_rerr) eq 0 then spectral_index_rerr = '0'
        if n_elements(spectral_index_rperr) eq 0 then spectral_index_rperr = '0'
        if n_elements(spectral_index_otpl_blue_cont) eq 0 then spectral_index_otpl_blue_cont = '0'
        if n_elements(spectral_index_otpl_red_cont) eq 0 then spectral_index_otpl_red_cont = '0'
        if n_elements(spectral_index_otpl_index) eq 0 then spectral_index_otpl_index = '0'
        if n_elements(spectral_index_botpl_blue_cont) eq 0 then spectral_index_botpl_blue_cont = '0'
        if n_elements(spectral_index_botpl_red_cont) eq 0 then spectral_index_botpl_red_cont = '0'
        if n_elements(spectral_index_botpl_index) eq 0 then spectral_index_botpl_index = '0'
        if n_elements(spectral_index_corr) eq 0 then spectral_index_corr = '0'
        if n_elements(spectral_index_cerr) eq 0 then spectral_index_cerr = '0'
        if n_elements(spectral_index_cperr) eq 0 then spectral_index_cperr = '0'
        if n_elements(si_bin_wave) eq 0 then si_bin_wave = '0'
        if n_elements(si_bin_flux) eq 0 then si_bin_flux = '0'
        if n_elements(si_bin_ivar) eq 0 then si_bin_ivar = '0'
        if n_elements(si_bin_mask) eq 0 then si_bin_mask = '0'
        if n_elements(si_otpl_flux) eq 0 then si_otpl_flux = '0'
        if n_elements(si_otpl_mask) eq 0 then si_otpl_mask = '0'
        if n_elements(si_botpl_flux) eq 0 then si_botpl_flux = '0'
        if n_elements(si_botpl_mask) eq 0 then si_botpl_mask = '0'
END


