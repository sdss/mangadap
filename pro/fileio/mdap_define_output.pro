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
;                           ypos=ypos, signal=signal, noise=noise, bin_par=bin_par, $
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
;                           eml_model=eml_model, elo_ew_kinematics_fit=elo_ew_kinematics_fit, $
;                           elo_ew_kinematics_err=elo_ew_kinematics_err, $
;                           elo_ew_kinematics_ind=elo_ew_kinematics_ind, $
;                           elo_ew_kinematics_ier=elo_ew_kinematics_ier, $
;                           elo_ew_sinst=elo_ew_sinst, elo_ew_omitted=elo_ew_omitted, $
;                           elo_ew_intens=elo_ew_intens, elo_ew_interr=elo_ew_interr, $
;                           elo_ew_fluxes=elo_ew_fluxes, elo_ew_flxerr=elo_ew_flxerr, $
;                           elo_ew_EWidth=elo_ew_EWidth, elo_ew_EW_err=elo_ew_EW_err, $
;                           elo_ew_eml_model=elo_ew_eml_model, $
;                           elo_fb_kinematics_avg=elo_fb_kinematics_avg, $
;                           elo_fb_kinematics_aer=elo_fb_kinematics_aer, $
;                           elo_fb_kinematics_ind=elo_fb_kinematics_ind, $
;                           elo_fb_kinematics_ier=elo_fb_kinematics_ier, $
;                           elo_fb_sinst=elo_fb_sinst, elo_fb_omitted=elo_fb_omitted, $
;                           elo_fb_intens=elo_fb_intens, elo_fb_interr=elo_fb_interr, $
;                           elo_fb_fluxes=elo_fb_fluxes, elo_fb_flxerr=elo_fb_flxerr, $
;                           elo_fb_EWidth=elo_fb_EWidth, elo_fb_EW_err=elo_fb_EW_err, $
;                           elo_fb_eml_model=elo_fb_eml_model, abs_par=abs_par, $
;                           abs_line_key=abs_line_key, $
;                           abs_line_indx_omitted=abs_line_indx_omitted, $
;                           abs_line_indx_val=abs_line_indx_val, $
;                           abs_line_indx_err=abs_line_indx_err, $
;                           abs_line_indx_otpl=abs_line_indx_otpl, $
;                           abs_line_indx_botpl=abs_line_indx_botpl, si_bin_wave=si_bin_wave, $
;                           si_bin_flux=si_bin_flux, si_bin_ivar=si_bin_ivar, $
;                           si_bin_mask=si_bin_mask, si_optimal_template=si_optimal_template, $
;                           si_broad_optimal_template=si_broad_optimal_template
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
;       elo_ew_kinematics_avg dblarr[B][2]
;       elo_ew_kinematics_aer dblarr[B][2]
;               The mean kinematics (and errors) of the emission lines
;               determined using Enci Wang's code.
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
;       elo_fb_kinematics_avg dblarr[B][2]
;       elo_fb_kinematics_aer dblarr[B][2]
;               The mean kinematics (and errors) of the emission lines
;               determined using Francesco Belfiore's code.
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
;       abs_line_indx_omitted intarr[B][I]
;               Flag that the spectral index I has (1) or has not (0) been
;               omitted from the measurements of binned spectrum B.  Spectral
;               indices are omitted if any part of their definition occurs
;               outside of the spectral range of the spectrum.
;
;       abs_line_indx_val dblarr[B][I]
;               Value of spectral index I for spectrum B.  These values have
;               been corrected for the velocity dispersion using the factor
;               derived from the measurements based on the broadened and
;               unbroadened optimal template.
;
;       abs_line_indx_err dblarr[B][I]
;               Error in spectral index I for spectrum B.
;
;       abs_line_indx_otpl dblarr[B][I]
;               Spectral index I measured using the optimal template for
;               spectrum B.
;
;       abs_line_indx_botpl dblarr[B][I]
;               Spectral index I measured using the broadened optimal template
;               for spectrum B.
;
;       si_bin_wave dblarr[D]
;               Wavelengths of the pixels in the binned spectra that have had
;               their resolution matched tot the spectral index system.  NOTE: D
;               can be (and is likely) different from C because of how the
;               resolution matching process censors the data.  See
;               MDAP_MATCH_SPECTRAL_RESOLUTION.
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
;       si_optimal_template dblarr[B][F]
;               The best-fitting template (sum of the weighted template in the
;               library) for each of the B galaxy spectra with the resolution
;               matched to that of the spectral index system.  NOTE: F can be
;               (and is likely) different from both D and C because of how the
;               resolution matching process censors the data.  See
;               MDAP_MATCH_SPECTRAL_RESOLUTION.
;
;       si_broad_optimal_template dblarr[B][F]
;               The best-fitting template (sum of the weighted template in the
;               library) for each of the B galaxy spectra with the resolution
;               matched to that of the spectral index system, and broadened by
;               the best-fitting line-of-sight velocity distribution.  TODO:
;               Does this include the velocity shift?
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
;       08 Dec 2014: (KBW) Accommodate radial binning and velocity
;                          registration
;       12 Dec 2014: (KBW) New format incorporating emission-line only
;                          results
;       09 Jan 2015: (KBW) Include instrumental dispersion for GANDALF fit
;-
;------------------------------------------------------------------------------


;-------------------------------------------------------------------------------
; Define all the supplied variables
PRO MDAP_DEFINE_OUTPUT, $
                header=header, dx=dx, dy=dy, w_range_sn=w_range_sn, xpos=xpos, ypos=ypos, $
                signal=signal, noise=noise, bin_par=bin_par, $
                threshold_ston_bin=threshold_ston_bin, bin_vreg=bin_vreg, bin_indx=bin_indx, $
                bin_weights=bin_weights, wave=wave, sres=sres, bin_flux=bin_flux, $
                bin_ivar=bin_ivar, bin_mask=bin_mask, xbin_rlow=xbin_rlow, ybin_rupp=ybin_rupp, $
                rbin=rbin, bin_area=bin_area, bin_ston=bin_ston, bin_n=bin_n, bin_flag=bin_flag, $
                w_range_analysis=w_range_analysis, $
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
                elo_ew_kinematics_ind=elo_ew_kinematics_ind, $
                elo_ew_kinematics_ier=elo_ew_kinematics_ier, elo_ew_sinst=elo_ew_sinst, $
                elo_ew_omitted=elo_ew_omitted, elo_ew_intens=elo_ew_intens, $
                elo_ew_interr=elo_ew_interr, elo_ew_fluxes=elo_ew_fluxes, $
                elo_ew_flxerr=elo_ew_flxerr, elo_ew_EWidth=elo_ew_EWidth, $
                elo_ew_EW_err=elo_ew_EW_err, elo_ew_eml_model=elo_ew_eml_model, $
                elo_fb_kinematics_avg=elo_fb_kinematics_avg, $
                elo_fb_kinematics_aer=elo_fb_kinematics_aer, $
                elo_fb_kinematics_ind=elo_fb_kinematics_ind, $
                elo_fb_kinematics_ier=elo_fb_kinematics_ier, elo_fb_sinst=elo_fb_sinst, $
                elo_fb_omitted=elo_fb_omitted, elo_fb_intens=elo_fb_intens, $
                elo_fb_interr=elo_fb_interr, elo_fb_fluxes=elo_fb_fluxes, $
                elo_fb_flxerr=elo_fb_flxerr, elo_fb_EWidth=elo_fb_EWidth, $
                elo_fb_EW_err=elo_fb_EW_err, elo_fb_eml_model=elo_fb_eml_model, abs_par=abs_par, $
                abs_line_key=abs_line_key, abs_line_indx_omitted=abs_line_indx_omitted, $
                abs_line_indx_val=abs_line_indx_val, abs_line_indx_err=abs_line_indx_err, $
                abs_line_indx_otpl=abs_line_indx_otpl, abs_line_indx_botpl=abs_line_indx_botpl, $
                si_bin_wave=si_bin_wave, si_bin_flux=si_bin_flux, si_bin_ivar=si_bin_ivar, $
                si_bin_mask=si_bin_mask, si_optimal_template=si_optimal_template, $
                si_broad_optimal_template=si_broad_optimal_template

        if n_elements(header) eq 0 then header = '0'
        if n_elements(dx) eq 0 then dx = '0'
        if n_elements(dy) eq 0 then dy = '0'
        if n_elements(w_range_sn) eq 0 then w_range_sn = '0'
        if n_elements(xpos) eq 0 then xpos = '0'
        if n_elements(ypos) eq 0 then ypos = '0'
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
        if n_elements(elo_ew_kinematics_avg) eq 0 then elo_ew_kinematics_avg = '0'
        if n_elements(elo_ew_kinematics_aer) eq 0 then elo_ew_kinematics_aer = '0'
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


