;+
; NAME:
;       MDAP_WRITE_OUTPUT
;
; PURPOSE:
;
;       This is an attempt at a flexible writing routine that will create the
;       primary DAP output files.  The main draw back at this point is that it
;       cannot shrink the binary tables.  Hopefully this can be done in
;       python...
;
;       EXTENSIONS ARE HARD-CODED:
;               0  - Empty
;               1  - DRPS    : BINTABLE
;               2  - BINS    : BINTABLE
;               3  - WAVE    : IMAGE
;               4  - SRES    : IMAGE
;               5  - FLUX    : IMAGE
;               6  - IVAR    : IMAGE
;               7  - MASK    : IMAGE
;               8  - ELPAR   : BINTABLE
;               9  - STFIT   : BINTABLE
;               10 - SMSK    : IMAGE
;               11 - SMOD    : IMAGE
;               12 - SGFIT   : BINTABLE
;               13 - SGMSK   : IMAGE
;               14 - SGMOD   : IMAGE
;               15 - ELMOD   : IMAGE
;               16 - ELOPAR  : BINTABLE
;               17 - ELOFIT  : BINTABLE
;               18 - ELOMEW  : IMAGE
;               19 - ELOMFB  : IMAGE
;               20 - SIWAVE  : IMAGE
;               21 - SIFLUX  : IMAGE
;               22 - SIIVAR  : IMAGE
;               23 - SIMASK  : IMAGE
;               24 - SIOTPL  : IMAGE
;               25 - SIBOTPL : IMAGE
;               26 - SIPAR   : BINTABLE
;               27 - SINDX   : BINTABLE
;
;       ANALYSIS FLAGS:
;                0 - spectrum not analyzed
;                1 - stellar continuum fit
;               91 - stellar continuum fit failed
;                2 - star+gas fit
;               92 - star+gas fit failed
;                3 - emission-line only fit
;               93 - emission-line only fit failed
;                4 - spectral index measurements
;               94 - spectral index measurements failed
;
;       Output data format:
;
;           Header keywords:
;               - Number of DRP spectra
;               - INPUT spaxel size in x and y (for each spaxel?)
;               - Wavelength range used for S/N calculation
;
;               - Number of binned spectra
;               - Binning type and parameters
;               - S/N threshold for bin inclusion
;
;               - Wavelength range of analysis
;               - S/N threshold for analysis
;
;               - Template library file and/or keyword
;               - Emission line file and/or keyword
;               - Absorption-line file and/or keyword
;
;           Binary tables:
;               - For each DRP spectrum: DRPS
;                 9 columns, ndrp rows
;                   - Fiducial X, Y
;                   - Fraction of good pixels in the spectrum
;                   - min(flux) == max(flux) flag
;                   - signal and noise
;                   - Registered velocity (taken out)
;                   - Index of bin to which DRP spectrum is assigned
;                   - Weight in binned spectrum
;
;               - For each binned spectrum: BINS
;                 7 columns, nbin rows
;                   - Luminosity weighted X, Y, R, area, S/N, and number of spectra
;                   - analysis flag(s)
;
;               - Gas emission-line structure (save this?): ELPAR
;                 4 columns, neml rows
;                   - gas line name: NNNNNN-LLLL, e.g. Ha-----6562, ArIII--7136
;                   - rest wavelength (in vacuum)
;                   - tied to?, tied type(t,v,s,''), doublet with?
;
;               - For stellar analysis of each binned spectrum: STFIT
;                 1 column, 1 [ntpl] column, 1 [nadd] column, 1 [nmult] column,
;                    2 [moments] columns, nbin rows
;                   - star-only chi^2/DOF
;                       TODO: change to chi2 and DOF separately
;                   - star-only template weights (vector)
;                   - star-only additive and multiplicative polynomials (vectors)
;                   - stellar kinematics and errors (vectors)
;
;               - For star+gas analysis of each binned spectrum: SGFIT
;                 1 column, 1 [ntpl] column, 1 [nmult] column, 2 [moments] columns,
;                 2 [nred] columns, 8 [neml] columns, 2 [neml][moments] columns, nbin rows
;                   - star+gas chi^2/DOF
;                       TODO: change to chi2 and DOF separately
;                   - star+gas template weights (vector)
;                   - star+gas multiplicative polynomials (vectors)
;                   - mean gas kinematics and errors (vectors)
;                   - reddening and errors (vectors)
;
;                   - For each emission line (only good fits):
;                       - gas line omission flag (vector)
;                       - gas kinematics and errors (vectors)
;                       - instrumental dispersion at the fitted wavelengths (vector)
;                       - gas intensities and errors (vectors)
;                       - gas fluxes and errors (vectors)
;                       - gas equivalent widths and errors (vectors)
;                       ;TODO: add gas reddening correction
;
;               - Emission-line-only line fit details(save this?): ELOPAR
;                 2 columns, nemlo rows
;                   - gas line name: NNNNNN-LLLL, e.g. Ha-----6562, ArIII--7136
;                   - rest wavelength (in vacuum)
;
;               - For the emission-line only fit of each binned spectrum: ELOFIT
;                 4 [moments] columns, 16 [nemlo] columns, 4 [nemlo][moments] columns, nbin rows
;                   - mean gas kinematics and errors, EW and FB (vectors)
;                   - For each emission line (only good fits):
;                       - gas line omission flag, EW and FB (vector)
;                       - gas kinematics and errors, EW and FB (vector)
;                       - instrumental dispersion at the fitted wavelength, EW and FB (vector)
;                       - gas intensities and errors, EW and FB (vector)
;                       - gas fluxes and errors, EW and FB (vector)
;                       - gas equivalent widths and errors, EW and FB (vector)
;
;               - Spectral index parameters: SIPAR
;                 2 columns, 3 [2] columns
;                   - Name
;                   - Passband (vector)
;                   - Blue continuum (vector)
;                   - Red continuum (vector)
;                   - Units
; 
;               - For spectral index meaurements of each binned spectrum: SINDX
;                 5 [nabs] columns
;                   - For each spectral index
;                       - omission flag (vector)
;                       - index value (vector)
;                       - index error (vector)
;                       - optimal template index value (vector)
;                       - broadened optimal template index value (vector)
;
;           Image Extensions:
;               - For binned spectra:
;                       - wavelength vector (1D): WAVE
;                       - spectral-resolution vector (1D): SRES
;                       - flux (2D): FLUX
;                       - inverse variance (2D): IVAR
;                       - mask (2D): MASK
;
;                       - star-only fit mask (2D): SMSK
;                       - best-fitting star-only model (2D): SMOD
;
;                       - star+gas fit mask (2D): SGMSK
;                       - best-fitting star+gas model (2D): SGMOD
;                       - best-fitting emission-line only model (2D): ELMOD
;                       - Best-fitting model from Enci Wang's
;                         emission-line only fit (2D): ELOEW
;                       - Best-fitting model from Francesco Belfiore's
;                         emission-line only fit (2D): ELOFB
;
;                       - spectral-index resolution matched wavelength (1D): SIWAVE
;                       -                                   flux (2D): SIFLUX
;                       -                                   ivar (2D): SIIVAR
;                       -                                   mask (2D): SIMASK
;                       -                                   optimal template (2D): SIOTPL
;                       -                                   optimal broadened template (2D): SIBOTPL
;
; CALLING SEQUENCE:
;       MDAP_WRITE_OUTPUT, file, header=header, dx=dx, dy=dy, w_range_sn=w_range_sn, xpos=xpos, $
;                          ypos=ypos, fraction_good=fraction_good, min_eq_max=min_eq_max, $
;                          signal=signal, noise=noise, bin_par=bin_par, $
;                          threshold_ston_bin=threshold_ston_bin, bin_vreg=bin_vreg, $
;                          bin_indx=bin_indx, bin_weights=bin_weights, wave=wave, sres=sres, $
;                          bin_flux=bin_flux, bin_ivar=bin_ivar, bin_mask=bin_mask, $
;                          xbin_rlow=xbin_rlow, ybin_rupp=ybin_rupp, rbin=rbin, bin_area=bin_area, $
;                          bin_ston=bin_ston, bin_n=bin_n, bin_flag=bin_flag, $
;                          w_range_analysis=w_range_analysis, $
;                          threshold_ston_analysis=threshold_ston_analysis, $
;                          tpl_library_key=tpl_library_key, eml_line_key=eml_line_key, $
;                          eml_par=eml_par, analysis_par=analysis_par, weights_ppxf=weights_ppxf, $
;                          add_poly_coeff_ppxf=add_poly_coeff_ppxf, $
;                          mult_poly_coeff_ppxf=mult_poly_coeff_ppxf, $
;                          stellar_kinematics_fit=stellar_kinematics_fit, $
;                          stellar_kinematics_err=stellar_kinematics_err, chi2_ppxf=chi2_ppxf, $
;                          obj_fit_mask_ppxf=obj_fit_mask_ppxf, bestfit_ppxf=bestfit_ppxf, $
;                          weights_gndf=weights_gndf, mult_poly_coeff_gndf=mult_poly_coeff_gndf, $
;                          emission_line_kinematics_avg=emission_line_kinematics_avg, $
;                          emission_line_kinematics_aer=emission_line_kinematics_aer, $
;                          chi2_gndf=chi2_gndf, $
;                          emission_line_kinematics_ind=emission_line_kinematics_ind, $
;                          emission_line_kinematics_ier=emission_line_kinematics_ier, $
;                          emission_line_sinst=emission_line_sinst, $
;                          emission_line_omitted=emission_line_omitted, $
;                          emission_line_intens=emission_line_intens, $
;                          emission_line_interr=emission_line_interr, $
;                          emission_line_fluxes=emission_line_fluxes, $
;                          emission_line_flxerr=emission_line_flxerr, $
;                          emission_line_EWidth=emission_line_EWidth, $
;                          emission_line_EW_err=emission_line_EW_err, $
;                          reddening_val=reddening_val, reddening_err=reddening_err, $
;                          obj_fit_mask_gndf=obj_fit_mask_gndf, bestfit_gndf=bestfit_gndf, $
;                          eml_model=eml_model, emlo_par=emlo_par, $
;                          elo_ew_kinematics_fit=elo_ew_kinematics_fit, $
;                          elo_ew_kinematics_err=elo_ew_kinematics_err, $
;                          elo_ew_kinematics_ind=elo_ew_kinematics_ind, $
;                          elo_ew_kinematics_ier=elo_ew_kinematics_ier, $
;                          elo_ew_sinst=elo_ew_sinst, elo_ew_omitted=elo_ew_omitted, $
;                          elo_ew_intens=elo_ew_intens, elo_ew_interr=elo_ew_interr, $
;                          elo_ew_fluxes=elo_ew_fluxes, elo_ew_flxerr=elo_ew_flxerr, $
;                          elo_ew_EWidth=elo_ew_EWidth, elo_ew_EW_err=elo_ew_EW_err, $
;                          elo_ew_eml_model=elo_ew_eml_model, $
;                          elo_fb_kinematics_avg=elo_fb_kinematics_avg, $
;                          elo_fb_kinematics_aer=elo_fb_kinematics_aer, $
;                          elo_fb_kinematics_ind=elo_fb_kinematics_ind, $
;                          elo_fb_kinematics_ier=elo_fb_kinematics_ier, $
;                          elo_fb_sinst=elo_fb_sinst, elo_fb_omitted=elo_fb_omitted, $
;                          elo_fb_intens=elo_fb_intens, elo_fb_interr=elo_fb_interr, $
;                          elo_fb_fluxes=elo_fb_fluxes, elo_fb_flxerr=elo_fb_flxerr, $
;                          elo_fb_EWidth=elo_fb_EWidth, elo_fb_EW_err=elo_fb_EW_err, $
;                          elo_fb_eml_model=elo_fb_eml_model, abs_par=abs_par, $
;                          abs_line_key=abs_line_key, abs_line_indx_omitted=abs_line_indx_omitted, $
;                          abs_line_indx_val=abs_line_indx_val, $
;                          abs_line_indx_err=abs_line_indx_err, $
;                          abs_line_indx_otpl=abs_line_indx_otpl, $
;                          abs_line_indx_botpl=abs_line_indx_botpl, si_bin_wave=si_bin_wave, $
;                          si_bin_flux=si_bin_flux, si_bin_ivar=si_bin_ivar, $
;                          si_bin_mask=si_bin_mask, si_optimal_template=si_optimal_template, $
;                          si_broad_optimal_template= si_broad_optimal_template, $
;                          /read_header, /quiet
;
; INPUTS:
;       file string
;               Name of the fits file for the output.  Will create/modify
;               existing data and/or create new extensions.
;
; OPTIONAL INPUTS:
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
;               Structure containing the binning parameters. Written to the header.
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

;   TODO: Radial binning parameters

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
;       threshold_ston_analysis double
;               S/N threshold for analysis of a binned spectrum.
;
;       tpl_library_key string
;               Keyword used to signify the template library.
;
;       ems_line_key string
;               Keyword used to signify the emission-line parameter file.
;
;       eml_par EmissionLine[]
;               An array of EmissionLine structures used to define the emission
;               lines to be fit by GANDALF.
;
;       analysis_par AnalysisPar structure
;               A structure that defines parameters used by PPXF and GANDALF in
;               the fitting procedure.  See its definition in
;               MDAP_DEFINE_ANALYSIS_PAR.pro
;
;       weights_ppxf dblarr[B][T]
;               Template weights for each spectrum obtained by PPXF.
;
;       add_poly_coeff_ppxf dblarr[B][AD]
;               Additive order AD legendre polynomial coefficients obtained for
;               each of N spectra by PPXF.
;
;       mult_poly_coeff_ppxf dblarr[B][MD]
;               Multiplicative order MD legendre polynomial coefficients obtained
;               for each of N spectra by PPXF.
;
;       stellar_kinematics_fit dblarr[B][K]
;               The best-fit stellar kinematics for each of the N fitted (or
;               fixed) input galaxy spectra.
;
;       stellar_kinematics_err dblarr[B][K]
;               Estimates of the errors in the best-fit stellar kinematics for
;               each of the B fitted input galaxy spectra.
;
;       chi2_ppxf dblarr[B]
;               Chi-square per degree of freedom obtained by the PPXF fit to
;               each of the B spectra.
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
;       emlo_par EmissionLine[O]
;               Structure with the emission-line parameters used during
;               the emission-line only fit.  The only things written to
;               the file are the name (the same as done for the GANDALF
;               parameter structure) and the rest wavelength of the line
;               (in vacuum).  The lines fit during the emission-line
;               only analysis are hard-coded (see
;               MDAP_DEFINE_EMISSION_LINES_ENCI_BELFIORE).
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
;       /read_header
;               Flag to read the primary header (exten=0) from the existing
;               ofile and edit it.  If a header is also provided to the
;               procedure, this tag overwrites it!
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
;       - If string table elements change, the overwrite command in mdap_setup
;         will not work.  FXBWRITM will likely throw an error.
;
; BUGS:
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;       15 Oct 2014: Original Implementation K. Westfall (KBW)
;       08 Dec 2014: (KBW) Write data associated with the radial binning
;                          and velocity registration.
;       12 Dec 2014: (KBW) New format incorporating emission-line only
;                          results
;       09 Jan 2015: (KBW) Include instrumental dispersion for GANDALF fit
;       09 Feb 2015: (KBW) Include fraction of good pixels and min(flux)
;                          == max(flux) flag in DRPS extension.
;-
;------------------------------------------------------------------------------

; Read the binary table column names
@mdap_set_output_file_cols

;------------------------------------------------------------------------------
; Determine the number of rows needed for the DRPS extention by getting the
; maximum length of the input vectors
FUNCTION MDAP_WRITE_OUTPUT_NUMBER_OF_DRPS, $
                xpos=xpos, ypos=ypos, fraction_good=fraction_good, min_eq_max=min_eq_max, $
                signal=signal, noise=noise, bin_vreg=bin_vreg, bin_indx=bin_indx, $
                bin_weights=bin_weights

        length = intarr(9)
        if n_elements(xpos) ne 0 then           length[0] = n_elements(xpos)
        if n_elements(ypos) ne 0 then           length[1] = n_elements(ypos)
        if n_elements(fraction_good) ne 0 then  length[2] = n_elements(fraction_good)
        if n_elements(min_eq_max) ne 0 then     length[3] = n_elements(min_eq_max)
        if n_elements(signal) ne 0 then         length[4] = n_elements(signal)
        if n_elements(noise) ne 0 then          length[5] = n_elements(noise)
        if n_elements(bin_vreg) ne 0 then       length[6] = n_elements(bin_vreg)
        if n_elements(bin_indx) ne 0 then       length[7] = n_elements(bin_indx)
        if n_elements(bin_weights) ne 0 then    length[8] = n_elements(bin_weights)
        return, max(length)
END
        
;------------------------------------------------------------------------------
; Determine the number of rows needed for the BINS extention by getting the
; maximum length of the input vectors
FUNCTION MDAP_WRITE_OUTPUT_NUMBER_OF_BINS, $
                bin_flux=bin_flux, bin_ivar=bin_ivar, bin_mask=bin_mask, xbin_rlow=xbin_rlow, $
                ybin_rupp=ybin_rupp, rbin=rbin, bin_area=bin_area, bin_ston=bin_ston, bin_n=bin_n, $
                bin_flag=bin_flag

        length = intarr(10)
        if n_elements(bin_flux) ne 0 then  length[0] = (size(bin_flux))[1]
        if n_elements(bin_ivar) ne 0 then  length[1] = (size(bin_ivar))[1]
        if n_elements(bin_mask) ne 0 then  length[2] = (size(bin_mask))[1]
        if n_elements(xbin_rlow) ne 0 then length[3] = n_elements(xbin_rlow)
        if n_elements(ybin_rupp) ne 0 then length[4] = n_elements(ybin_rupp)
        if n_elements(rbin) ne 0 then      length[5] = n_elements(rbin)
        if n_elements(bin_area) ne 0 then  length[6] = n_elements(bin_area)
        if n_elements(bin_ston) ne 0 then  length[7] = n_elements(bin_ston)
        if n_elements(bin_n) ne 0 then     length[8] = n_elements(bin_n)
        if n_elements(bin_flag) ne 0 then  length[9] = n_elements(bin_flag)
        return, max(length)

END

;-------------------------------------------------------------------------------
; Add the bin type and its parameter(s) to the header
PRO MDAP_WRITE_OUTPUT_UPDATE_HEADER_BIN, $
                header, bin_par
        SXADDPAR, header, 'BINTYPE', bin_par.type, 'Type of binning performed (keyword)'

        if bin_par.type ne 'RADIAL' then begin
            SXADDPAR, header, 'BINSCL', 1.0d, 'Coordinate scaling for binning (arcsec)'
        endif else $
            SXADDPAR, header, 'BINSCL', bin_par.rscale, 'Radial binning: Radial scaling (arcsec)'

        if bin_par.type eq 'NONE' then $            ; Nothing left to print
            return

        if bin_par.v_register eq 1 then begin
            SXADDPAR, header, 'BINVRG', 'True', 'Velocity register spectra before binning'
        endif else $
            SXADDPAR, header, 'BINVRG', 'False', 'Velocity register spectra before binning'

        if bin_par.optimal_weighting eq 1 then begin
            SXADDPAR, header, 'BINWGT', 'Optimal', 'Weighting when combining spectra'
        endif else $
            SXADDPAR, header, 'BINWGT', 'Uniform', 'Weighting when combining spectra'

        if bin_par.type eq 'ALL' then $             ; Nothing left to print
            return

        if bin_par.type eq 'STON' then begin
            SXADDPAR, header, 'BINSN', bin_par.ston, 'Voronoi binning: Minimum S/N'
            return                                  ; Nothing left to print
        endif

        ; Binning type is 'RADIAL'
        SXADDPAR, header, 'BINCX', bin_par.cx, 'Radial binning: Ellipse X-center'
        SXADDPAR, header, 'BINCY', bin_par.cy, 'Radial binning: Ellipse Y-center'
        SXADDPAR, header, 'BINPA', bin_par.pa, 'Radial binning: Ellipse position angle'
        SXADDPAR, header, 'BINELL', bin_par.ell, 'Radial binning: Ellipticity'
END

;-------------------------------------------------------------------------------
; Correct the naxis keywords in the header based on the number of elements in sz
; (number of axes) and the value in each element of sz (e.g., NAXIS1 = sz[0])
PRO MDAP_WRITE_OUTPUT_FIX_NAXES, $
                header, sz=sz

        MDAP_FITS_ZERO_NAXES, header            ; Remove all NAXIS parameter and set NAXIS=0

        if n_elements(sz) eq 0 then $
            return                              ; No axes!

        ndim = n_elements(sz)                   ; sz is a vector with the size in each dimension
        SXADDPAR, header, 'NAXIS', ndim, 'Number of axes'       ; Number of dimensions
        for i=0,ndim-1 do begin
            a=strcompress(string(i+1), /remove_all)
            SXADDPAR, header, 'NAXIS'+a, sz[0], 'Length of axis '+a     ; Size of ath dimension
        endfor
END

;-------------------------------------------------------------------------------
; Update the primary header and write it.
;       file and header are the fits file name and its current header
;   - Parameters provided are added
;   - Any NAXISn parameters are removed and NAXIS is set to 0
;   - Any XTENSION parameter is removed
;   - If the file exists, the modification date is added, and the header is
;     modified
;   - If the file does not exist, it is written with the header in extension 0
;
PRO MDAP_WRITE_OUTPUT_UPDATE_HEADER, $
                file, header, ndrp=ndrp, dx=dx, dy=dy, w_range_sn=w_range_sn, nbin=nbin, $
                bin_par=bin_par, threshold_ston_bin=threshold_ston_bin, $
                w_range_analysis=w_range_analysis, $
                threshold_ston_analysis=threshold_ston_analysis, analysis_par=analysis_par, $
                tpl_library_key=tpl_library_key, ems_line_key=ems_line_key, $
                abs_line_key=abs_line_key

        ; Add the number of DRP spectra
        if n_elements(ndrp) ne 0 then $
            SXADDPAR, header, 'NS_DRP', ndrp, 'Number of DRP produced spectra'

        ; Add the spaxel sizes
        if n_elements(dx) ne 0 then $
            SXADDPAR, header, 'SPAXDX', dx, 'Spaxel size in X (arcsec)'
        if n_elements(dy) ne 0 then $
            SXADDPAR, header, 'SPAXDY', dy, 'Spaxel size in Y (arcsec)'

        ; Add the wavelength range used to calculate the S/N per pixel
        if n_elements(w_range_sn) ne 0 then begin
            if n_elements(w_range_sn) ne 2 then $
                message, 'Wave range must have two elements!'
            SXADDPAR, header, 'SNWAVE1', w_range_sn[0], $
                      'Starting wavelength for S/N calculation (ang)'
            SXADDPAR, header, 'SNWAVE2', w_range_sn[1], $
                      'Ending wavelength for S/N calculation (ang)'
        endif

        ; Add the number of binned spectra
        if n_elements(nbin) ne 0 then $
            SXADDPAR, header, 'NS_BIN', nbin, 'Number of binned spectra'

        ; Add the bin type used
        if n_elements(bin_par) ne 0 then $
            MDAP_WRITE_OUTPUT_UPDATE_HEADER_BIN, header, bin_par

        ; Add the threshold S/N for inclusion of spectra in any bin
        if n_elements(threshold_ston_bin) ne 0 then begin
            SXADDPAR, header, 'MINSNBIN', threshold_ston_bin, $
                      'Minimum S/N of spectrum to include in any bin'
        endif

        ; Add the wavelength range used for the analysis
        if n_elements(w_range_analysis) ne 0 then begin
            if n_elements(w_range_analysis) ne 2 then $
                message, 'Wave range must have two elements!'
            SXADDPAR, header, 'FITWAVE1', w_range_analysis[0], $
                      'Starting wavelength for analysis (ang)'
            SXADDPAR, header, 'FITWAVE2', w_range_analysis[1], $
                      'Ending wavelength for analysis (ang)'
        endif

        ; Add the threshold S/N for analysis of the spectra
        if n_elements(threshold_ston_analysis) ne 0 then begin
            SXADDPAR, header, 'MINSNFIT', threshold_ston_analysis, $
                      'Minimum S/N of spectrum to analyze'
        endif

        ; Print the extra parameters
        ; TODO: Move to an update_analysis_par_header procedure
        if n_elements(analysis_par) ne 0 then begin
            SXADDPAR, header, 'MOMENTS', analysis_par.moments, 'Number of stellar kinematic moments'
            SXADDPAR, header, 'DEGREE', analysis_par.degree, 'Additive polynomial degree (PPXF)'
            SXADDPAR, header, 'MDEGREE', analysis_par.mdegree, 'Multiplicative polynomial degree'
            SXADDPAR, header, 'REDORD', analysis_par.reddening_order, 'Reddening order (0,1,2)'
            if analysis_par.reddening_order gt 0 then begin $
                str='[' +MDAP_STC(analysis_par.reddening[0])
                if analysis_par.reddening_order eq 2 then $
                    str=str+','+MDAP_STC(analysis_par.reddening[1])
                str=str+']'
                SXADDPAR, header, 'REDINP', str, 'Input reddening parameters'
            endif
            SXADDPAR, header, 'ZSINST', analysis_par.zero_instr_disp, $
                              'Force zero sigma_instr in gas fit (0-false)'
        endif

        ; Add the template, emission-line, and absorption-line keywords
        if n_elements(tpl_library_key) ne 0 then $
            SXADDPAR, header, 'TPLKEY', tpl_library_key, 'Template library identifier'
        if n_elements(ems_line_key) ne 0 then $
            SXADDPAR, header, 'EMLKEY', ems_line_key, 'Emission-line parameter identifier'
        if n_elements(abs_line_key) ne 0 then $
            SXADDPAR, header, 'SPIKEY', abs_line_key, 'Spectral-index parameter identifier'

        MDAP_WRITE_OUTPUT_FIX_NAXES, header             ; Remove all NAXIS elements, setting NAXIS=0
        ; This is the PRIMARY HDU so...
        SXDELPAR, header, 'XTENSION'                    ; ... remove extension keyword
        SXDELPAR, header, 'EXTNAME'                     ; ... remove extension name
        SXADDPAR, header, 'EXTEND', 'T', 'FITS data may contain extensions'; ... allow extensions

        ; Update the file header
        if file_test(file) eq 1 then begin
            MDAP_FITS_HEADER_ADD_DATE, header, /modified        ; Add/change last date modified
            MODFITS, file, 0, header, exten_no=0                ; Modify the header only
        endif else begin
            MDAP_FITS_HEADER_ADD_DATE, header                   ; Add parameter with original date
            WRITEFITS, file, 0, header                          ; Write the header (no data)
        endelse
END

;-------------------------------------------------------------------------------
; Add a binary table column that is single valued (ie. not a vector or matrix).
; The optional keywords are used to set the type.  Currently supported are
; double (dbl), integer(int), long integer (lng), string (str with a size set by
; dummystr).
PRO MDAP_FXBADDCOL_VALUE, $
                indx, hdr, colname, comment, dbl=dbl, int=int, lng=lng, str=str, dummystr=dummystr
        if keyword_set(dbl) then begin
            FXBADDCOL, indx, hdr, 1.0d, colname, comment                ; double
        endif else if keyword_set(int) then begin
            FXBADDCOL, indx, hdr, 0, colname, comment                   ; integer
        endif else if keyword_set(lng) then begin
            FXBADDCOL, indx, hdr, 0L, colname, comment                  ; long
        endif else if keyword_set(str) then begin
            if n_elements(dummystr) eq 0 then $
                dummystr = ' '
            FXBADDCOL, indx, hdr, dummystr, colname, comment            ; string of length dummystr
        endif else $
            message, 'Unknown type for table column!'
END

;-------------------------------------------------------------------------------
; Add a binary table column that is a vector of length nn for each row.  If the
; size is 0 or less, the column is added but with only a single element.  The
; optional keywords are used to set the type.  Currently supported are double
; (dbl), integer(int), long integer (lng), string (str with a size set by
; dummystr).
PRO MDAP_FXBADDCOL_VECTOR, $
                indx, hdr, nn, colname, comment, dbl=dbl, int=int, str=str, dummystr=dummystr
        
        if nn gt 0 then begin
            if keyword_set(dbl) then begin
                FXBADDCOL, indx, hdr, dblarr(nn), colname, comment
            endif else if keyword_set(int) then begin
                FXBADDCOL, indx, hdr, intarr(nn), colname, comment
            endif else if keyword_set(str) then begin
                if n_elements(dummystr) eq 0 then $
                    dummystr = ' '
                FXBADDCOL, indx, hdr, make_array(nn, /string, value=dummystr), colname, comment
            endif else $
                message, 'Unknown type for table column!'
        endif else begin
            MDAP_FXBADDCOL_VALUE, indx, hdr, colname, comment, dbl=dbl, int=int, str=str, $
                                  dummystr=dummystr
        endelse

END

;-------------------------------------------------------------------------------
; Add a binary table column that is a matrix of size mm (column) X nn (rows) for
; each row.  If nn and mm are not both greater than 0, the column is added but
; with only a 1X1 matrix.  The optional keywords are used to set the type.
; Currently supported are double (dbl), integer(int), long integer (lng), string
; (str with a size set by dummystr).
PRO MDAP_FXBADDCOL_MATRIX, $
                indx, hdr, nn, mm, colname, comment, dbl=dbl, int=int, str=str, dummystr=dummystr
        
;       if nn gt 0 and mm gt 0 then begin
        if nn gt 0 && mm gt 0 then begin
            if keyword_set(dbl) then begin
                FXBADDCOL, indx, hdr, dblarr(nn, mm), colname, comment
            endif else if keyword_set(int) then begin
                FXBADDCOL, indx, hdr, intarr(nn, mm), colname, comment
            endif else if keyword_set(str) then begin
                if n_elements(dummystr) eq 0 then $
                    dummystr = ' '
                FXBADDCOL, indx, hdr, make_array(nn, mm, /string, value=dummystr), colname, comment
            endif else $
                message, 'Unknown type for table column!'
        endif else begin
            MDAP_FXBADDCOL_MATRIX, indx, hdr, 1, 1, colname, comment, dbl=dbl, int=int, str=str, $
                                   dummystr=dummystr
        endelse

END

;-------------------------------------------------------------------------------
; Initialize the DRPS extension
;       - If the file and extension exist, nothing is done
;       - If the file does not exist, MDAP_INITIALIZE_FITS_FILE will create it
;       - If the extension does not exist, the table is instantiated with all
;         the columns and with a single row
PRO MDAP_WRITE_OUTPUT_DRPS_INITIALIZE, $
                file

        ; Initialize the file (does nothing if the file already exists)
        MDAP_INITIALIZE_FITS_FILE, file

        if MDAP_CHECK_EXTENSION_EXISTS(file, 'DRPS') eq 0 then begin
            ; Create a base level header; only one row for now!
            FXBHMAKE, bth, 1, 'DRPS', 'Binary table with properties of DRP spectra'
            MDAP_FITS_HEADER_ADD_DATE, bth, /modified           ; Add/change last date modified

            ; FROM mdap_set_output_file_cols.pro (included above)
            ; TODO: Make this a loop, with a function that also provides the comments, column types
            cols = MDAP_SET_DRPS_COLS()
            
            ; Create the table columns using placeholders to define the column data type
            MDAP_FXBADDCOL_VALUE, 1, bth, cols[0], ' Fiducial X position of spaxel (arcsec)', /dbl
            MDAP_FXBADDCOL_VALUE, 2, bth, cols[1], ' Fiducial Y position of spaxel (arcsec)', /dbl
            MDAP_FXBADDCOL_VALUE, 3, bth, cols[2], ' Fraction of good pixels in spectrum', /dbl
            MDAP_FXBADDCOL_VALUE, 4, bth, cols[3], ' min(flux)=max(flux) in spectrum', /int
            MDAP_FXBADDCOL_VALUE, 5, bth, cols[4], ' Mean flux/pixel from SNWAVE1->SNWAVE2', /dbl
            MDAP_FXBADDCOL_VALUE, 6, bth, cols[5], ' Mean noise/pixel from SNWAVE1->SNWAVE2', /dbl
            MDAP_FXBADDCOL_VALUE, 7, bth, cols[6], ' Velocity used for deredshifting', /dbl
            MDAP_FXBADDCOL_VALUE, 8, bth, cols[7], ' Index of bin for spectrum (-1 for no bin)',/int
            MDAP_FXBADDCOL_VALUE, 9, bth, cols[8], ' Weight of spectrum in bin', /dbl

            FXBCREATE, tbl, file, bth                   ; Create the binary table extension
            FXBFINISH, tbl                              ; Close up
;           print, tbl
            free_lun, tbl
        endif
END

;-------------------------------------------------------------------------------
; Initialize the BINS extension
;       - If the file and extension exist, nothing is done
;       - If the file does not exist, MDAP_INITIALIZE_FITS_FILE will create it
;       - If the extension does not exist, the table is instantiated with all
;         the columns and with a single row
PRO MDAP_WRITE_OUTPUT_BINS_INITIALIZE, $
                file

        ; Initialize the file (does nothing if the file already exists)
        MDAP_INITIALIZE_FITS_FILE, file

        if MDAP_CHECK_EXTENSION_EXISTS(file, 'BINS') eq 0 then begin

            ; Create a base level header; only one row for now!
            FXBHMAKE, bth, 1, 'BINS', 'Binary table with properties of binned spectra'
            MDAP_FITS_HEADER_ADD_DATE, bth, /modified           ; Add/change last date modified
            
            ; FROM mdap_set_output_file_cols.pro (included above)
            ; TODO: Make this a loop, with a function that also provides the comments, column types
            cols = MDAP_SET_BINS_COLS()
            
            ; Create the table columns using placeholders to define the column data type
            MDAP_FXBADDCOL_VALUE, 1, bth, cols[0], $
                                 ' L-weighted X position (arcsec) or bin lower radius (BINSCL)', $
                                  /dbl
            MDAP_FXBADDCOL_VALUE, 2, bth, cols[1], $
                                 ' L-weighted Y position (arcsec) or bin upper radius (BINSCL)', $
                                 /dbl
            MDAP_FXBADDCOL_VALUE, 3, bth, cols[2], ' Luminosity weighted radius (BINSCL)', /dbl
            MDAP_FXBADDCOL_VALUE, 4, bth, cols[3], ' Area of bin (arcsec^2)', /dbl
            MDAP_FXBADDCOL_VALUE, 5, bth, cols[4], ' S/N of bin', /dbl
            MDAP_FXBADDCOL_VALUE, 6, bth, cols[5], ' Number of spectra in bin', /lng
            MDAP_FXBADDCOL_VALUE, 7, bth, cols[6], ' Analysis flag', /int

            FXBCREATE, tbl, file, bth                   ; Create the binary table extension
            FXBFINISH, tbl                              ; Close up
;           print, tbl
            free_lun, tbl
        endif
END

;-------------------------------------------------------------------------------
; Initialize the ELPAR extension
;       - If the file and extension exist, nothing is done
;       - If the file does not exist, MDAP_INITIALIZE_FITS_FILE will create it
;       - If the extension does not exist, the table is instantiated with all
;         the columns and with a single row
PRO MDAP_WRITE_OUTPUT_ELPAR_INITIALIZE, $
                file

        ; Initialize the file (does nothing if the file already exists)
        MDAP_INITIALIZE_FITS_FILE, file

        if MDAP_CHECK_EXTENSION_EXISTS(file, 'ELPAR') eq 0 then begin
            ; Create a base level header; only one row for now!
            FXBHMAKE, bth, 1, 'ELPAR', 'Binary table with the emission-line parameters'
            MDAP_FITS_HEADER_ADD_DATE, bth, /modified           ; Add/change last date modified
            
            ; FROM mdap_set_output_file_cols.pro (included above)
            ; TODO: Make this a loop, with a function that also provides the comments, column types
            cols = MDAP_SET_ELPAR_COLS()
            
            ; Create the table columns using placeholders to define the column data type
            MDAP_SET_EMISSION_LINE_NAME_DEFAULT, dstr
            MDAP_FXBADDCOL_VALUE, 1, bth, cols[0], ' Emission-line identifier', /str, $
                                  dummystr=dstr
            MDAP_FXBADDCOL_VALUE, 2, bth, cols[1], ' Rest wavelength (ang)', /dbl
            MDAP_FXBADDCOL_VALUE, 3, bth, cols[2], ' Kinematics are tied to line in this row',/int
            MDAP_SET_EMISSION_LINE_TIED_DEFAULT, dstr
            MDAP_FXBADDCOL_VALUE, 4, bth, cols[3], ' Tied type: v-vel, s-sig, t-both', /str, $
                                  dummystr=dstr
            MDAP_FXBADDCOL_VALUE, 5, bth, cols[4], ' Line is a doublet of line in this row', /int

            FXBCREATE, tbl, file, bth                   ; Create the binary table extension
            FXBFINISH, tbl                              ; Close up
;           print, tbl
            free_lun, tbl
        endif
END

;-------------------------------------------------------------------------------
; Compare the current size of a table to the requested size input to an
; initialization.  If the new size is defined (!= -1) and different from the
; existing size, select to modify the existing table to the new size.  cur_size
; is replaced with the new size for the table (if different from the current
; size).
PRO MDAP_WRITE_OUTPUT_COMPARE_TABLE_SIZE, $
                cur_size, inp_size, modify

        modify = 0                              ; Initialize to NOT modify
        nn = n_elements(cur_size)               ; Get the number of columns to check
        for i=0,nn-1 do begin
            ; If the new size is defined (!= -1) and different from the existing
            ; size, select to modify the existing table to the new size
;           if (inp_size[i] ne -1 and cur_size[i] ne inp_size[i]) then begin
            if (inp_size[i] ne -1 && cur_size[i] ne inp_size[i]) then begin
                cur_size[i] = inp_size[i]
                modify = 1
            endif
        endfor
END

;-------------------------------------------------------------------------------
; Initialize the STFIT extension
;       - If the file does not exist, MDAP_INITIALIZE_FITS_FILE will create it
;       - If the extension exists, check that the sizes match the input
;       - If the size is the same as the existing extension, finish
;       - Otherwise, modify the fits extension, either by creating it or
;         modifying its column prperties
PRO MDAP_WRITE_OUTPUT_STELLAR_KIN_INITIALIZE, $
                file, weights_ppxf=weights_ppxf, add_poly_coeff_ppxf=add_poly_coeff_ppxf, $
                mult_poly_coeff_ppxf=mult_poly_coeff_ppxf, $
                stellar_kinematics_fit=stellar_kinematics_fit, $
                stellar_kinematics_err=stellar_kinematics_err

        ; Initialize the file (does nothing if the file already exists)
        MDAP_INITIALIZE_FITS_FILE, file

        ; Check if the extension exists
        extension_exists = MDAP_CHECK_EXTENSION_EXISTS(file, 'STFIT')

        ; Determine the dimensions of the input vectors
        inp_size = make_array(4, /int, value=-1)

        if n_elements(weights_ppxf) ne 0 then $
            inp_size[0] = (size(weights_ppxf))[2]               ; Number of templates
        if n_elements(add_poly_coeff_ppxf) ne 0 then $
            inp_size[1] = (size(add_poly_coeff_ppxf))[2]        ; Order of additive polynomial
        if n_elements(mult_poly_coeff_ppxf) ne 0 then $
            inp_size[2] = (size(mult_poly_coeff_ppxf))[2]       ; Order of multiplicative polynomial

        ; Number of kinematic moments
        if n_elements(stellar_kinematics_fit) ne 0 then begin
            inp_size[3] = (size(stellar_kinematics_fit))[2]
        endif else if n_elements(stellar_kinematics_err) ne 0 then $
            inp_size[3] = (size(stellar_kinematics_err))[2]

        ; FROM mdap_set_output_file_cols.pro (included above)
        ; TODO: Make this a loop, with a function that also provides the comments, column types
        cols = MDAP_SET_STFIT_COLS()
            
        modify = 0
        nrows = 1                       ; TODO: Is nrows needed?
        if extension_exists eq 1 then begin             ; If the extension exists...

            ; Get the dimensions of the existing table
            fxbopen, unit, file, 'STFIT', bth           ; Open the extension
            nrows = fxpar(bth, 'NAXIS2')                ; Get the existing number of rows
            cur_size = intarr(4)
            cur_size[0] = fxbdimen(unit, cols[0])       ; Number of TPLW elements
            cur_size[1] = fxbdimen(unit, cols[1])       ; Number of ADDPOLY elements
            cur_size[2] = fxbdimen(unit, cols[2])       ; Number of MULTPOLY elements
            cur_size[3] = fxbdimen(unit, cols[3])       ; Number of KIN elements, same as KINERR
            fxbclose, unit                              ; Close the file
;           print, unit
            free_lun, unit                              ; Free the LUN

            ; Compare the current size to the existing size and decide if the
            ; size needs to be modified
            MDAP_WRITE_OUTPUT_COMPARE_TABLE_SIZE, cur_size, inp_size, modify
        endif else begin                                ; If it does not exist ...
            cur_size = inp_size                         ; Set the modified size to the input size
            modify = 1
        endelse

        ; Extension exists and all the columns have the correct size, return
        if modify eq 0 then $
            return

        ; Create the header
        FXBHMAKE, bth, nrows, 'STFIT', 'Binary table with stellar kinematics fit results'
        MDAP_FITS_HEADER_ADD_DATE, bth, /modified               ; Add/change last date modified
            
        ; Create the table columns using placeholders to define the column data type and size
        MDAP_FXBADDCOL_VECTOR, 1, bth, cur_size[0], cols[0], ' Template weights (PPXF)', /dbl
        MDAP_FXBADDCOL_VECTOR, 2, bth, cur_size[1], cols[1], $
                               ' Coefficients of additive polynomial (PPXF)', /dbl
        MDAP_FXBADDCOL_VECTOR, 3, bth, cur_size[2], cols[2], $
                               ' Coefficients of multiplicative polynomial (PPXF)', /dbl
        MDAP_FXBADDCOL_VECTOR, 4, bth, cur_size[3], cols[3], ' Stellar kinematics', /dbl
        MDAP_FXBADDCOL_VECTOR, 5, bth, cur_size[3], cols[4], ' Stellar kinematic errors', /dbl
        MDAP_FXBADDCOL_VALUE, 6, bth, cols[5], ' Reduced chi-square (chi-square/DOF; PPXF)', /dbl

        if extension_exists eq 0 then begin             ; If the extension does not exist...
            FXBCREATE, tbl, file, bth                   ; Create the binary table extension
                                                        ; TODO: Write fake data to it?
            FXBFINISH, tbl                              ; Close up
;           print, tbl
            free_lun, tbl
        endif else begin                                ; If the extension does exist...
            bytdb = bytarr(fxpar(bth, 'NAXIS1'), fxpar(bth, 'NAXIS2'))  ; Allocate data
            MODFITS, file, bytdb, bth, extname='STFIT'  ; Modify the header and allocate the data
        endelse

END

;-------------------------------------------------------------------------------
; Initialize the SGFIT extension
;       - If the file does not exist, MDAP_INITIALIZE_FITS_FILE will create it
;       - If the extension exists, check that the sizes match the input
;       - If the size is the same as the existing extension, finish
;       - Otherwise, modify the fits extension, either by creating it or
;         modifying its column prperties
PRO MDAP_WRITE_OUTPUT_STAR_AND_GAS_FIT_INITIALIZE, $
                file, eml_par=eml_par, weights_gndf=weights_gndf, $
                mult_poly_coeff_gndf=mult_poly_coeff_gndf, $
                emission_line_kinematics_avg=emission_line_kinematics_avg, $
                emission_line_kinematics_aer=emission_line_kinematics_aer, $
                emission_line_kinematics_ind=emission_line_kinematics_ind, $
                emission_line_kinematics_ier=emission_line_kinematics_ier, $
                emission_line_sinst=emission_line_sinst, $
                emission_line_omitted=emission_line_omitted, $
                emission_line_intens=emission_line_intens, $
                emission_line_interr=emission_line_interr, $
                emission_line_fluxes=emission_line_fluxes, $
                emission_line_flxerr=emission_line_flxerr, $
                emission_line_EWidth=emission_line_EWidth, $
                emission_line_EW_err=emission_line_EW_err, reddening_val=reddening_val, $
                reddening_err=reddening_err

        ; Initialize the file (does nothing if the file already exists)
        MDAP_INITIALIZE_FITS_FILE, file

        ; Check if the extension exists
        extension_exists = MDAP_CHECK_EXTENSION_EXISTS(file, 'SGFIT')

        ; Determine the dimensions of the input vectors
        inp_size = make_array(5, /int, value=-1)
        if n_elements(weights_gndf) ne 0 then $
            inp_size[0] = (size(weights_gndf))[2]               ; Number of templates
        if n_elements(mult_poly_coeff_gndf) ne 0 then $
            inp_size[1] = (size(mult_poly_coeff_gndf))[2]       ; Order of multiplicative polynomial

        ; Number of kinematic moments
        if n_elements(emission_line_kinematics_avg) ne 0 then begin
            inp_size[2]=(size(emission_line_kinematics_avg))[2]
        endif else if n_elements(emission_line_kinematics_aer) ne 0 then begin
            inp_size[2]=(size(emission_line_kinematics_aer))[2]
        endif else if n_elements(emission_line_kinematics_ind) ne 0 then begin
            inp_size[2]=(size(emission_line_kinematics_ind))[3]
        endif else if n_elements(emission_line_kinematics_ier) ne 0 then $
            inp_size[2]=(size(emission_line_kinematics_ier))[3]

        ; Number of reddening coefficients
        if n_elements(reddening_val) ne 0 then begin
            inp_size[3]=(size(reddening_val))[2]
        endif else if n_elements(reddening_err) ne 0 then $
            inp_size[3]=(size(reddening_err))[2]

        ; Number of emission-lines
        if n_elements(emission_line_kinematics_ind) ne 0 then begin
            inp_size[4]=(size(emission_line_kinematics_ind))[2]
        endif else if n_elements(emission_line_kinematics_ier) ne 0 then begin
            inp_size[4]=(size(emission_line_kinematics_ier))[2]
        endif else if n_elements(emission_line_sinst) ne 0 then begin
            inp_size[4]=(size(emission_line_sinst))[2]
        endif else if n_elements(emission_line_omitted) ne 0 then begin
            inp_size[4]=(size(emission_line_omitted))[2]
        endif else if n_elements(emission_line_intens) ne 0 then begin
            inp_size[4]=(size(emission_line_intens))[2]
        endif else if n_elements(emission_line_interr) ne 0 then begin
            inp_size[4]=(size(emission_line_interr))[2]
        endif else if n_elements(emission_line_fluxes) ne 0 then begin
            inp_size[4]=(size(emission_line_fluxes))[2]
        endif else if n_elements(emission_line_flxerr) ne 0 then begin
            inp_size[4]=(size(emission_line_flxerr))[2]
        endif else if n_elements(emission_line_EWidth) ne 0 then begin
            inp_size[4]=(size(emission_line_EWidth))[2]
        endif else if n_elements(emission_line_EW_err) ne 0 then $
            inp_size[4]=(size(emission_line_EW_err))[2]

        ; TODO: Check the number of emission lines against the ELPAR extension?
;       if n_elements(eml_par) ne 0 and inp_size[4] gt 0 then $
        if n_elements(eml_par) ne 0 && inp_size[4] gt 0 then $
            if n_elements(eml_par) ne inp_size[4] then $
                message, 'Number of line structures does not match the number of line results!'

        ; FROM mdap_set_output_file_cols.pro (included above)
        ; TODO: Make this a loop, with a function that also provides the comments, column types
        cols = MDAP_SET_SGFIT_COLS()
            
        modify = 0
        nrows = 1                       ; TODO: Is nrows needed?
        if extension_exists eq 1 then begin             ; If the extension exists...

            ; Get the dimensions of the existing table
            fxbopen, unit, file, 'SGFIT', bth           ; Open the file
            nrows = fxpar(bth, 'NAXIS2')                ; Number of rows
            cur_size = intarr(5)
            cur_size[0] = fxbdimen(unit, cols[0])       ; Number of TPLW elements
            cur_size[1] = fxbdimen(unit, cols[1])       ; Number of MULTPOLY elements
            cur_size[2] = fxbdimen(unit, cols[2])       ; Number of KIN elements, same as KINERR
            cur_size[3] = fxbdimen(unit, cols[5])       ; Number of reddening elements
            cur_size[4] = fxbdimen(unit, cols[7])       ; Number of EML elements
            fxbclose, unit                              ; Close the file
;           print, unit
            free_lun, unit                              ; Free the LUN

            ; Compare the current size to the existing size and decide if the
            ; size needs to be modified
            MDAP_WRITE_OUTPUT_COMPARE_TABLE_SIZE, cur_size, inp_size, modify
        endif else begin                                ; If it does not exist ...
            cur_size = inp_size                         ; Set the modified size to the input size
            modify = 1
        endelse

        ; Extension exists and all the columns have the correct size, return
        if modify eq 0 then $
            return

        ; Create the header
        FXBHMAKE, bth, nrows, 'SGFIT', 'Binary table with star+gas fit results'
        MDAP_FITS_HEADER_ADD_DATE, bth, /modified               ; Add/change last date modified
            
        ; Create the table columns using placeholders to define the column data type and size
        MDAP_FXBADDCOL_VECTOR, 1, bth, cur_size[0], cols[0], ' Template weights (GANDALF)', /dbl
        MDAP_FXBADDCOL_VECTOR, 2, bth, cur_size[1], cols[1], $
                               ' Coefficients of multiplicative polynomial (GANDALF)', /dbl
        MDAP_FXBADDCOL_VECTOR, 3, bth, cur_size[2], cols[2], ' Average emission-line kinematics', $
                               /dbl
        MDAP_FXBADDCOL_VECTOR, 4, bth, cur_size[2], cols[3], $
                               ' Average emission-line kinematic errors', /dbl
        MDAP_FXBADDCOL_VALUE, 5, bth, cols[4], ' Reduced chi-square (chi-square/DOF; GANDALF)', /dbl

        MDAP_FXBADDCOL_VECTOR, 6, bth, cur_size[3], cols[5], ' Reddening (GANDALF)', /dbl
        MDAP_FXBADDCOL_VECTOR, 7, bth, cur_size[3], cols[6], ' Reddening error (GANDALF)', /dbl
        MDAP_FXBADDCOL_VECTOR, 8, bth, cur_size[4], cols[7], ' Emission-line omission flag', /int
        MDAP_FXBADDCOL_VECTOR, 9, bth, cur_size[4], cols[8], ' Emission-line amplitude', /dbl
        MDAP_FXBADDCOL_VECTOR, 10, bth, cur_size[4], cols[9], ' Emission-line amplitude error', /dbl
        MDAP_FXBADDCOL_MATRIX, 11, bth, cur_size[4], cur_size[2], cols[10], $
                               ' Individual emission-line kinematics', /dbl
        MDAP_FXBADDCOL_MATRIX, 12, bth, cur_size[4], cur_size[2], cols[11], $
                               ' Individual emission-line kinematic errors', /dbl
        MDAP_FXBADDCOL_VECTOR, 13, bth, cur_size[4], cols[12], $
                               ' Instrumental dispersion (sigma) at emission-line centroids', /dbl
        MDAP_FXBADDCOL_VECTOR, 13, bth, cur_size[4], cols[13], ' Emission-line flux', /dbl
        MDAP_FXBADDCOL_VECTOR, 14, bth, cur_size[4], cols[14], ' Emission-line flux error', /dbl
        MDAP_FXBADDCOL_VECTOR, 15, bth, cur_size[4], cols[15], ' Emission-line equivalent width', $
                               /dbl
        MDAP_FXBADDCOL_VECTOR, 16, bth, cur_size[4], cols[16], $
                               ' Emission-line equivalent width error', /dbl

        if extension_exists eq 0 then begin             ; If the extension does not exist...
            FXBCREATE, tbl, file, bth                   ; Create the binary table extension
                                                        ; TODO: Write fake data to it?
            FXBFINISH, tbl                              ; Close up
;           print, tbl
            free_lun, tbl
        endif else begin                                ; If the extension does exist...
            bytdb = bytarr(fxpar(bth, 'NAXIS1'), fxpar(bth, 'NAXIS2'))
            MODFITS, file, bytdb, bth, extname='SGFIT'  ; Modify the header and allocate the data
        endelse
END

;-------------------------------------------------------------------------------
; Initialize the ELOPAR extension
;       - If the file and extension exist, nothing is done
;       - If the file does not exist, MDAP_INITIALIZE_FITS_FILE will create it
;       - If the extension does not exist, the table is instantiated with all
;         the columns and with a single row
PRO MDAP_WRITE_OUTPUT_ELOPAR_INITIALIZE, $
                file

        ; Initialize the file (does nothing if the file already exists)
        MDAP_INITIALIZE_FITS_FILE, file

        if MDAP_CHECK_EXTENSION_EXISTS(file, 'ELOPAR') eq 0 then begin
            ; Create a base level header; only one row for now!
            FXBHMAKE, bth, 1, 'ELOPAR', 'Binary table with emission-lines fit by EM-only spectra'
            MDAP_FITS_HEADER_ADD_DATE, bth, /modified           ; Add/change last date modified
            
            ; FROM mdap_set_output_file_cols.pro (included above)
            ; TODO: Make this a loop, with a function that also provides the comments, column types
            cols = MDAP_SET_ELOPAR_COLS()
            
            ; Create the table columns using placeholders to define the column data type
            MDAP_SET_EMISSION_LINE_NAME_DEFAULT, dstr
            MDAP_FXBADDCOL_VALUE, 1, bth, cols[0], ' Emission-line identifier', /str, $
                                  dummystr=dstr
            MDAP_FXBADDCOL_VALUE, 2, bth, cols[1], ' Rest wavelength (ang)', /dbl

            FXBCREATE, tbl, file, bth                   ; Create the binary table extension
            FXBFINISH, tbl                              ; Close up
;           print, tbl
            free_lun, tbl
        endif
END

;-------------------------------------------------------------------------------
; Initialize the SGFIT extension
;       - If the file does not exist, MDAP_INITIALIZE_FITS_FILE will create it
;       - If the extension exists, check that the sizes match the input
;       - If the size is the same as the existing extension, finish
;       - Otherwise, modify the fits extension, either by creating it or
;         modifying its column prperties
PRO MDAP_WRITE_OUTPUT_EMISSION_LINE_FIT_INITIALIZE, $
                file, emlo_par=emlo_par, elo_ew_kinematics_avg=elo_ew_kinematics_avg, $
                elo_ew_kinematics_aer=elo_ew_kinematics_aer, $
                elo_ew_kinematics_ind=elo_ew_kinematics_ind, $
                elo_ew_kinematics_ier=elo_ew_kinematics_ier, elo_ew_sinst=elo_ew_sinst, $
                elo_ew_omitted=elo_ew_omitted, elo_ew_intens=elo_ew_intens, $
                elo_ew_interr=elo_ew_interr, elo_ew_fluxes=elo_ew_fluxes, $
                elo_ew_flxerr=elo_ew_flxerr, elo_ew_EWidth=elo_ew_EWidth, $
                elo_ew_EW_err=elo_ew_EW_err, elo_fb_kinematics_avg=elo_fb_kinematics_avg, $
                elo_fb_kinematics_aer=elo_fb_kinematics_aer, $
                elo_fb_kinematics_ind=elo_fb_kinematics_ind, $
                elo_fb_kinematics_ier=elo_fb_kinematics_ier, elo_fb_sinst=elo_fb_sinst, $
                elo_fb_omitted=elo_fb_omitted, elo_fb_intens=elo_fb_intens, $
                elo_fb_interr=elo_fb_interr, elo_fb_fluxes=elo_fb_fluxes, $
                elo_fb_flxerr=elo_fb_flxerr, elo_fb_EWidth=elo_fb_EWidth, $
                elo_fb_EW_err=elo_fb_EW_err

        ; Initialize the file (does nothing if the file already exists)
        MDAP_INITIALIZE_FITS_FILE, file

        ; Check if the extension exists
        extension_exists = MDAP_CHECK_EXTENSION_EXISTS(file, 'ELOFIT')

        ; Determine the dimensions of the input vectors
        inp_size = make_array(2, /int, value=-1)

        ; Number of kinematic moments
        if n_elements(elo_ew_kinematics_avg) ne 0 then begin
            inp_size[0]=(size(elo_ew_kinematics_avg))[2]
        endif else if n_elements(elo_ew_kinematics_aer) ne 0 then begin
            inp_size[0]=(size(elo_ew_kinematics_aer))[2]
        endif else if n_elements(elo_ew_kinematics_ind) ne 0 then begin
            inp_size[0]=(size(elo_ew_kinematics_ind))[3]
        endif else if n_elements(elo_ew_kinematics_ier) ne 0 then begin
            inp_size[0]=(size(elo_ew_kinematics_ier))[3]
        endif else if n_elements(elo_fb_kinematics_avg) ne 0 then begin
            inp_size[0]=(size(elo_fb_kinematics_avg))[2]
        endif else if n_elements(elo_fb_kinematics_aer) ne 0 then begin
            inp_size[0]=(size(elo_fb_kinematics_aer))[2]
        endif else if n_elements(elo_fb_kinematics_ind) ne 0 then begin
            inp_size[0]=(size(elo_fb_kinematics_ind))[3]
        endif else if n_elements(elo_fb_kinematics_ier) ne 0 then $
            inp_size[0]=(size(elo_fb_kinematics_ier))[3]

        ; Number of emission-lines
        if n_elements(elo_ew_kinematics_ind) ne 0 then begin
            inp_size[1]=(size(elo_ew_kinematics_ind))[2]
        endif else if n_elements(elo_ew_kinematics_ier) ne 0 then begin
            inp_size[1]=(size(elo_ew_kinematics_ier))[2]
        endif else if n_elements(elo_ew_sinst) ne 0 then begin
            inp_size[1]=(size(elo_ew_sinst))[2]
        endif else if n_elements(elo_ew_omitted) ne 0 then begin
            inp_size[1]=(size(elo_ew_omitted))[2]
        endif else if n_elements(elo_ew_intens) ne 0 then begin
            inp_size[1]=(size(elo_ew_intens))[2]
        endif else if n_elements(elo_ew_interr) ne 0 then begin
            inp_size[1]=(size(elo_ew_interr))[2]
        endif else if n_elements(elo_ew_fluxes) ne 0 then begin
            inp_size[1]=(size(elo_ew_fluxes))[2]
        endif else if n_elements(elo_ew_flxerr) ne 0 then begin
            inp_size[1]=(size(elo_ew_flxerr))[2]
        endif else if n_elements(elo_ew_EWidth) ne 0 then begin
            inp_size[1]=(size(elo_ew_EWidth))[2]
        endif else if n_elements(elo_ew_EW_err) ne 0 then begin
            inp_size[1]=(size(elo_ew_EW_err))[2]
        endif else if n_elements(elo_fb_kinematics_ind) ne 0 then begin
            inp_size[1]=(size(elo_fb_kinematics_ind))[2]
        endif else if n_elements(elo_fb_kinematics_ier) ne 0 then begin
            inp_size[1]=(size(elo_fb_kinematics_ier))[2]
        endif else if n_elements(elo_fb_sinst) ne 0 then begin
            inp_size[1]=(size(elo_fb_sinst))[2]
        endif else if n_elements(elo_fb_omitted) ne 0 then begin
            inp_size[1]=(size(elo_fb_omitted))[2]
        endif else if n_elements(elo_fb_intens) ne 0 then begin
            inp_size[1]=(size(elo_fb_intens))[2]
        endif else if n_elements(elo_fb_interr) ne 0 then begin
            inp_size[1]=(size(elo_fb_interr))[2]
        endif else if n_elements(elo_fb_fluxes) ne 0 then begin
            inp_size[1]=(size(elo_fb_fluxes))[2]
        endif else if n_elements(elo_fb_flxerr) ne 0 then begin
            inp_size[1]=(size(elo_fb_flxerr))[2]
        endif else if n_elements(elo_fb_EWidth) ne 0 then begin
            inp_size[1]=(size(elo_fb_EWidth))[2]
        endif else if n_elements(elo_fb_EW_err) ne 0 then $
            inp_size[1]=(size(elo_fb_EW_err))[2]

        ; TODO: Check the number of emission lines against the ELOPAR extension?
;       if n_elements(emlo_par) ne 0 and inp_size[1] gt 0 then $
        if n_elements(emlo_par) ne 0 && inp_size[1] gt 0 then $
            if n_elements(emlo_par) ne inp_size[1] then $
                message, 'Number of line structures does not match the number of line results!'

        ; FROM mdap_set_output_file_cols.pro (included above)
        ; TODO: Make this a loop, with a function that also provides the comments, column types
        cols = MDAP_SET_ELOFIT_COLS()
            
        modify = 0
        nrows = 1                       ; TODO: Is nrows needed?
        if extension_exists eq 1 then begin             ; If the extension exists...

            ; Get the dimensions of the existing table
            fxbopen, unit, file, 'ELOFIT', bth           ; Open the file
            nrows = fxpar(bth, 'NAXIS2')                ; Number of rows
            cur_size = intarr(2)
            cur_size[0] = fxbdimen(unit, cols[0])       ; Number of EW_KIN elements (num of moments)
            cur_size[1] = fxbdimen(unit, cols[2])       ; Number of emission lines
            fxbclose, unit                              ; Close the file
;           print, unit
            free_lun, unit                              ; Free the LUN

            ; Compare the current size to the existing size and decide if the
            ; size needs to be modified
            MDAP_WRITE_OUTPUT_COMPARE_TABLE_SIZE, cur_size, inp_size, modify
        endif else begin                                ; If it does not exist ...
            cur_size = inp_size                         ; Set the modified size to the input size
            modify = 1
        endelse

        ; Extension exists and all the columns have the correct size, return
        if modify eq 0 then $
            return

        ; Create the header
        FXBHMAKE, bth, nrows, 'ELOFIT', 'Binary table with emission-line-only fit results'
        MDAP_FITS_HEADER_ADD_DATE, bth, /modified               ; Add/change last date modified
            
        ; Create the table columns using placeholders to define the column data type and size
        MDAP_FXBADDCOL_VECTOR, 1, bth, cur_size[0], cols[0], $
                               ' Enci: Average emission-line kinematics', /dbl
        MDAP_FXBADDCOL_VECTOR, 2, bth, cur_size[0], cols[1], $
                               ' Enci: Average emission-line kinematic errors', /dbl
        MDAP_FXBADDCOL_VECTOR, 3, bth, cur_size[1], cols[2], $
                               ' Enci: Emission-line omission flag', /int
        MDAP_FXBADDCOL_VECTOR, 4, bth, cur_size[1], cols[3], $
                               ' Enci: Emission-line amplitude', /dbl
        MDAP_FXBADDCOL_VECTOR, 5, bth, cur_size[1], cols[4], $
                               ' Enci: Emission-line amplitude error', /dbl
        MDAP_FXBADDCOL_MATRIX, 6, bth, cur_size[1], cur_size[0], cols[5], $
                               ' Enci: Individual emission-line kinematics', /dbl
        MDAP_FXBADDCOL_MATRIX, 7, bth, cur_size[1], cur_size[0], cols[6], $
                               ' Enci: Individual emission-line kinematic errors', /dbl
        MDAP_FXBADDCOL_VECTOR, 8, bth, cur_size[1], cols[7], $
                               ' Enci: Instrumental dispersion at fitted centroid', /dbl
        MDAP_FXBADDCOL_VECTOR, 9, bth, cur_size[1], cols[8], $
                               ' Enci: Emission-line flux', /dbl
        MDAP_FXBADDCOL_VECTOR, 10, bth, cur_size[1], cols[9], $
                               ' Enci: Emission-line flux error', /dbl
        MDAP_FXBADDCOL_VECTOR, 11, bth, cur_size[1], cols[10], $
                               ' Enci: Emission-line equivalent width', /dbl
        MDAP_FXBADDCOL_VECTOR, 12, bth, cur_size[1], cols[11], $
                               ' Enci: Emission-line equivalent width error', /dbl
        MDAP_FXBADDCOL_VECTOR, 13, bth, cur_size[0], cols[12], $
                               ' Belfiore: Average emission-line kinematics', /dbl
        MDAP_FXBADDCOL_VECTOR, 14, bth, cur_size[0], cols[13], $
                               ' Belfiore: Average emission-line kinematic errors', /dbl
        MDAP_FXBADDCOL_VECTOR, 15, bth, cur_size[1], cols[14], $
                               ' Belfiore: Emission-line omission flag', /int
        MDAP_FXBADDCOL_VECTOR, 16, bth, cur_size[1], cols[15], $
                               ' Belfiore: Emission-line amplitude', /dbl
        MDAP_FXBADDCOL_VECTOR, 17, bth, cur_size[1], cols[16], $
                               ' Belfiore: Emission-line amplitude error', /dbl
        MDAP_FXBADDCOL_MATRIX, 18, bth, cur_size[1], cur_size[0], cols[17], $
                               ' Belfiore: Individual emission-line kinematics', /dbl
        MDAP_FXBADDCOL_MATRIX, 19, bth, cur_size[1], cur_size[0], cols[18], $
                               ' Belfiore: Individual emission-line kinematic errors', /dbl
        MDAP_FXBADDCOL_VECTOR, 20, bth, cur_size[1], cols[19], $
                               ' Belfiore: Instrumental dispersion at fitted centroid', /dbl
        MDAP_FXBADDCOL_VECTOR, 21, bth, cur_size[1], cols[20], $
                               ' Belfiore: Emission-line flux', /dbl
        MDAP_FXBADDCOL_VECTOR, 22, bth, cur_size[1], cols[21], $
                               ' Belfiore: Emission-line flux error', /dbl
        MDAP_FXBADDCOL_VECTOR, 23, bth, cur_size[1], cols[22], $
                               ' Belfiore: Emission-line equivalent width', /dbl
        MDAP_FXBADDCOL_VECTOR, 24, bth, cur_size[1], cols[23], $
                               ' Belfiore: Emission-line equivalent width error', /dbl

        if extension_exists eq 0 then begin             ; If the extension does not exist...
            FXBCREATE, tbl, file, bth                   ; Create the binary table extension
                                                        ; TODO: Write fake data to it?
            FXBFINISH, tbl                              ; Close up
;           print, tbl
            free_lun, tbl
        endif else begin                                ; If the extension does exist...
            bytdb = bytarr(fxpar(bth, 'NAXIS1'), fxpar(bth, 'NAXIS2'))
            MODFITS, file, bytdb, bth, extname='ELOFIT' ; Modify the header and allocate the data
        endelse
END

;-------------------------------------------------------------------------------
; Set the default size for the spectral index names, which is just at most 10
; spaces for the index name.
PRO MDAP_SET_SPECTRAL_INDEX_NAME_DEFAULT, $
                str
        str='NNNNNNNNNN'
END

;-------------------------------------------------------------------------------
; Setup the default size for the spectral index unit type.
PRO MDAP_SET_SPECTRAL_INDEX_UNIT_DEFAULT, $
                str
        str='NNN'
END

;-------------------------------------------------------------------------------
; Initialize the SIPAR extension
;       - If the file and extension exist, nothing is done
;       - If the file does not exist, MDAP_INITIALIZE_FITS_FILE will create it
;       - If the extension does not exist, the table is instantiated with all
;         the columns and with a single row
PRO MDAP_WRITE_OUTPUT_SIPAR_INITIALIZE, $
                file

        ; Initialize the file (does nothing if the file already exists)
        MDAP_INITIALIZE_FITS_FILE, file

        if MDAP_CHECK_EXTENSION_EXISTS(file, 'SIPAR') eq 0 then begin
            ; Create a base level header; only one row for now!
            FXBHMAKE, bth, 1, 'SIPAR', 'Binary table with the spectral index parameters'
            MDAP_FITS_HEADER_ADD_DATE, bth, /modified           ; Add/change last date modified
            
            ; FROM mdap_set_output_file_cols.pro (included above)
            ; TODO: Make this a loop, with a function that also provides the comments, column types
            cols = MDAP_SET_SIPAR_COLS()
            
            ; Create the table columns using placeholders to define the column data type
            MDAP_SET_SPECTRAL_INDEX_NAME_DEFAULT, dstr
            MDAP_FXBADDCOL_VALUE, 1, bth, cols[0], ' Spectral index identifier', /str, $
                                  dummystr=dstr
            MDAP_FXBADDCOL_VECTOR, 2, bth, 2, cols[1], ' Primary passband of the index (ang)', $
                                   /dbl
            MDAP_FXBADDCOL_VECTOR, 3, bth, 2, cols[2], $
                                   ' Continuum band blueward of the index (ang)', /dbl
            MDAP_FXBADDCOL_VECTOR, 4, bth, 2, cols[3], $
                                   ' Continuum band redward of the index (ang)', /dbl
            MDAP_SET_SPECTRAL_INDEX_UNIT_DEFAULT, dstr
            MDAP_FXBADDCOL_VALUE, 5, bth, cols[4], ' Spectral index unit (mag or ang)', /str, $
                                  dummystr=dstr

            FXBCREATE, tbl, file, bth                   ; Create the binary table extension
            FXBFINISH, tbl                              ; Close up
;           print, tbl
            free_lun, tbl
        endif
END

;-------------------------------------------------------------------------------
; Initialize the SGFIT extension
;       - If the file does not exist, MDAP_INITIALIZE_FITS_FILE will create it
;       - If the extension exists, check that the sizes match the input
;       - If the size is the same as the existing extension, finish
;       - Otherwise, modify the fits extension, either by creating it or
;         modifying its column prperties
PRO MDAP_WRITE_OUTPUT_SPECTRAL_INDICES_INITIALIZE, $
                file, abs_par=abs_par, abs_line_indx_omitted=abs_line_indx_omitted, $
                abs_line_indx_val=abs_line_indx_val, abs_line_indx_err=abs_line_indx_err, $
                abs_line_indx_otpl=abs_line_indx_otpl, abs_line_indx_botpl=abs_line_indx_botpl

        ; Initialize the file (does nothing if the file already exists)
        MDAP_INITIALIZE_FITS_FILE, file

        ; Check if the extension exists
        extension_exists = MDAP_CHECK_EXTENSION_EXISTS(file, 'SINDX')

        ; Determine the dimensions of the input vectors
        inp_size = make_array(1, /int, value=-1)

        ; Number of spectral indices
        if n_elements(abs_line_indx_omitted) ne 0 then begin
            inp_size[0]=(size(abs_line_indx_omitted))[2]
        endif else if n_elements(abs_line_indx_val) ne 0 then begin
            inp_size[0]=(size(abs_line_indx_val))[2]
        endif else if n_elements(abs_line_indx_err) ne 0 then begin
            inp_size[0]=(size(abs_line_indx_err))[2]
        endif else if n_elements(abs_line_indx_otpl) ne 0 then begin
            inp_size[0]=(size(abs_line_indx_otpl))[2]
        endif else if n_elements(abs_line_indx_botpl) ne 0 then $
            inp_size[0]=(size(abs_line_indx_botpl))[2]

        ; TODO: Check the number of spectral indices against the SIPAR extension?
;       if n_elements(abs_par) ne 0 and inp_size[0] gt 0 then $
        if n_elements(abs_par) ne 0 && inp_size[0] gt 0 then $
            if n_elements(abs_par) ne inp_size[0] then $
                message, 'Number of index structures does not match the number of index results!'

        ; FROM mdap_set_output_file_cols.pro (included above)
        ; TODO: Make this a loop, with a function that also provides the comments, column types
        cols = MDAP_SET_SINDX_COLS()
            
        modify = 0
        nrows = 1                       ; TODO: Is nrows needed?
        if extension_exists eq 1 then begin             ; If the extension exists...

            ; Get the dimensions of the existing table
            fxbopen, unit, file, 'SINDX', bth           ; Open the file
            nrows = fxpar(bth, 'NAXIS2')                ; Number of rows
            cur_size = intarr(1)
            cur_size[0] = fxbdimen(unit, cols[0])       ; Number of spectral index elements
            fxbclose, unit                              ; Close the file
;           print, unit
            free_lun, unit                              ; Free the LUN

            ; Compare the current size to the existing size and decide if the
            ; size needs to be modified
            MDAP_WRITE_OUTPUT_COMPARE_TABLE_SIZE, cur_size, inp_size, modify
        endif else begin                                ; If it does not exist ...
            cur_size = inp_size                         ; Set the modified size to the input size
            modify = 1
        endelse

        ; Extension exists and all the columns have the correct size, return
        if modify eq 0 then $
            return

        ; Create the header
        FXBHMAKE, bth, nrows, 'SINDX', 'Binary table with spectral index measurements'
        MDAP_FITS_HEADER_ADD_DATE, bth, /modified               ; Add/change last date modified
            
        ; Create the table columns using placeholders to define the column data type and size
        MDAP_FXBADDCOL_VECTOR, 1, bth, cur_size[0], cols[0], ' Spectral index omission flag', /int
        MDAP_FXBADDCOL_VECTOR, 2, bth, cur_size[0], cols[1], ' Spectral index value', /dbl
        MDAP_FXBADDCOL_VECTOR, 3, bth, cur_size[0], cols[2], ' Spectral index error', /dbl
        MDAP_FXBADDCOL_VECTOR, 4, bth, cur_size[0], cols[3], $
                               ' Spectral index based on the optimal template', /dbl
        MDAP_FXBADDCOL_VECTOR, 5, bth, cur_size[0], cols[4], $
                               ' Spectral index based on the broadened optimal template', /dbl

        if extension_exists eq 0 then begin             ; If the extension does not exist...
            FXBCREATE, tbl, file, bth                   ; Create the binary table extension
                                                        ; TODO: Write fake data to it?
            FXBFINISH, tbl                              ; Close up
;           print, tbl
            free_lun, tbl
        endif else begin                                ; If the extension does exist...
            bytdb = bytarr(fxpar(bth, 'NAXIS1'), fxpar(bth, 'NAXIS2'))
            MODFITS, file, bytdb, bth, extname='SINDX'  ; Modify the header and allocate the data
        endelse
END

;-------------------------------------------------------------------------------
; Attempt to match the size of the input file extention to the input number of
; rows.  Currently this only works if the length is equal to or greater than the
; existing number of rows.
PRO MDAP_WRITE_OUTPUT_TBL_MATCH_SIZE, $
                file, length, extension

        FXBOPEN, unit, file, extension, hdr, access='RW'        ; Open the file

        nrows = fxpar(hdr, 'NAXIS2')                            ; Get the number of rows
        if nrows eq 0 then $
            message, 'NAXIS2 not available in header of '+file

        if length eq nrows then begin                           ; Length and nrows are the same
            FXBFINISH, unit                                         ; Close up
            free_lun, unit
            return
        endif

        if length lt nrows then $                               ; Currently cannot shrink tables
            message, 'Cannot shrink table.  Must overwrite table.'

        if length gt nrows then $
            FXBGROW, unit, hdr, length                          ; Grow the table

        FXBFINISH, unit                                         ; Close up
;       print, unit
        free_lun, unit
END

;-------------------------------------------------------------------------------
; Setup emission-line names, six spaces for the name and 5 spaces for
; the truncated wavelength: NNNNNN-LLLLL (e.g., Ha-----6562)
FUNCTION MDAP_SET_EMISSION_LINE_NAME, $
                eml_par

        neml = n_elements(eml_par)                      ; Number of lines
        elname = strarr(neml)                           ; Array to hold the lines
        for i=0,neml-1 do begin
            sln=strlen(eml_par[i].name)                 ; Length of the name string
            wl=strcompress(string(fix(eml_par[i].lambda)), /remove_all) ; Wavelength string
            slw=strlen(wl)                              ; Length of the wavelength string
            dash=''                                     ; Add dashes between name and wavelength
            for j=0,11-sln-slw-1 do $
                dash = dash+'-'
            elname[i]= eml_par[i].name+dash+wl
        endfor
        return, elname                                  ; Return the string
END

;-------------------------------------------------------------------------------
; Set the default size for the emission-line names, six spaces for the name and
; 5 spaces for the truncated wavelength: NNNNNN-LLLLL (e.g., Ha-----6562)
PRO MDAP_SET_EMISSION_LINE_NAME_DEFAULT, $
                str
        str='NNNNNNNNNNNN'
END

;-------------------------------------------------------------------------------
; Setup the emission-line tied type value.  The type is only a single character,
; but because of IDL weirdness (or more likely my use of it), the type has to be
; 3 characters.
FUNCTION MDAP_SET_EMISSION_LINE_TIED, $
                type
        if strlen(type) ne 1 then $
            message, 'Type is not the correct length!'
        return, ' '+type+' '
END

;-------------------------------------------------------------------------------
; Setup the default size for the emission-line tied type.
PRO MDAP_SET_EMISSION_LINE_TIED_DEFAULT, $
                str
        str='NNN'
END

;-------------------------------------------------------------------------------
; Setup the spectral index names, which is just at most 10 spaces for the index
; name.
FUNCTION MDAP_SET_SPECTRAL_INDEX_NAME, $
                abs_par

        nabs = n_elements(abs_par)                      ; Number of indices
        siname = strarr(nabs)                           ; Array to hold the line names

        ; Get the length of the default string
        ;dummystr='----------' ;MDAP_SET_SPECTRAL_INDEX_NAME_DEFAULT
        MDAP_SET_SPECTRAL_INDEX_NAME_DEFAULT, dstr
        deflen = strlen(dstr)

        ; Get the names
        for i=0,nabs-1 do begin
            if strlen(abs_par[i].name) gt deflen then $ ; Make sure it's not too long
                message, 'Name of spectral index is too long!'

            siname[i] = abs_par[i].name                 ; Set the name
        endfor

        return, siname                                  ; Return the string
END

;-------------------------------------------------------------------------------
; Setup the spectral index unit name. The unit should only be only 'mag' or
; 'ang'.
FUNCTION MDAP_SET_SPECTRAL_INDEX_UNIT, $
                unit
        if strlen(unit) ne 3 then $
            message, 'Unit is not the correct length!'
        return, unit
END

;-------------------------------------------------------------------------------
; Check the input vectors to write to the DRPS extension.  Returns a flag that
; there is something to write, and the size of the vectors to write.  Will throw
; an error if the input vectors are not the same size
PRO MDAP_WRITE_OUTPUT_DRPS_CHECK_INPUTS, $
                something_to_write, ninp, xpos=xpos, ypos=ypos, fraction_good=fraction_good, $
                min_eq_max=min_eq_max, signal=signal, noise=noise, bin_vreg=bin_vreg, $
                bin_indx=bin_indx, bin_weights=bin_weights

        ; Check that ndrp matches the size of one of the existing inputs
        ; TODO: Assumes all input vectors have the same length!
        nel = intarr(9)
        nel[0] = n_elements(xpos)
        nel[1] = n_elements(ypos)
        nel[2] = n_elements(fraction_good)
        nel[3] = n_elements(min_eq_max)
        nel[4] = n_elements(signal)
        nel[5] = n_elements(noise)
        nel[6] = n_elements(bin_vreg)
        nel[7] = n_elements(bin_indx)
        nel[8] = n_elements(bin_weights)
        ninp = max(nel)

        for i=0,8 do begin
;           if nel[i] gt 0 and nel[i] ne ninp then $
            if nel[i] gt 0 && nel[i] ne ninp then $
                message, 'Input vectors for DRPS have different lengths!'
        endfor

        something_to_write = ninp ne 0
END

;-------------------------------------------------------------------------------
; Check the input vectors to write to the BINS extension.  Returns a flag that
; there is something to write, and the size of the vectors to write.  Will throw
; an error if the input vectors are not the same size
PRO MDAP_WRITE_OUTPUT_BINS_CHECK_INPUTS, $
                something_to_write, ninp, xbin_rlow=xbin_rlow, ybin_rupp=ybin_rupp, rbin=rbin, $
                bin_area=bin_area, bin_ston=bin_ston, bin_n=bin_n, bin_flag=bin_flag

        ; Check that ndrp matches the size of one of the existing inputs
        ; TODO: Assumes all input vectors have the same length!
        nel = intarr(7)
        nel[0] = n_elements(xbin_rlow)
        nel[1] = n_elements(ybin_rupp)
        nel[2] = n_elements(rbin)
        nel[3] = n_elements(bin_area)
        nel[4] = n_elements(bin_ston)
        nel[5] = n_elements(bin_n)
        nel[6] = n_elements(bin_flag)
        ninp = max(nel)

        for i=0,5 do begin
;           if nel[i] gt 0 and nel[i] ne ninp then $
            if nel[i] gt 0 && nel[i] ne ninp then $
                message, 'Input vectors for BINS have different lengths!'
        endfor

        something_to_write = ninp ne 0
END

;-------------------------------------------------------------------------------
; Check the input vectors to write to the STFIT extension.  Returns a flag that
; there is something to write, and the size of the vectors to write.  Will throw
; an error if the input vectors are not the same size
PRO MDAP_WRITE_OUTPUT_STELLAR_KIN_CHECK_INPUTS, $
                something_to_write, ninp, weights_ppxf=weights_ppxf, $
                add_poly_coeff_ppxf=add_poly_coeff_ppxf, $
                mult_poly_coeff_ppxf=mult_poly_coeff_ppxf, $
                stellar_kinematics_fit=stellar_kinematics_fit, $
                stellar_kinematics_err=stellar_kinematics_err, chi2_ppxf=chi2_ppxf

        ; Check that ndrp matches the size of one of the existing inputs
        ; TODO: Assumes all input vectors have the same length!
        nel = intarr(6)
        nel[0] = n_elements(weights_ppxf) eq 0 ? 0 : (size(weights_ppxf))[1]
        nel[1] = n_elements(add_poly_coeff_ppxf) eq 0 ? 0 : (size(add_poly_coeff_ppxf))[1]
        nel[2] = n_elements(mult_poly_coeff_ppxf) eq 0 ? 0 : (size(mult_poly_coeff_ppxf))[1]
        nel[3] = n_elements(stellar_kinematics_fit) eq 0 ? 0 : (size(stellar_kinematics_fit))[1]
        nel[4] = n_elements(stellar_kinematics_err) eq 0 ? 0 : (size(stellar_kinematics_err))[1]
        nel[5] = n_elements(chi2_ppxf)
        ninp = max(nel)

        for i=0,5 do begin
;           if nel[i] gt 0 and nel[i] ne ninp then $
            if nel[i] gt 0 && nel[i] ne ninp then $
                message, 'Input vectors or STFIT have different lengths!'
        endfor

        something_to_write = ninp ne 0
END

;-------------------------------------------------------------------------------
; Check the input vectors to write to the STFIT extension.  Returns a flag that
; there is something to write, and the size of the vectors to write.  Will throw
; an error if the input vectors are not the same size
PRO MDAP_WRITE_OUTPUT_STAR_AND_GAS_FIT_CHECK_INPUTS, $
                something_to_write, ninp, weights_gndf=weights_gndf, $
                mult_poly_coeff_gndf=mult_poly_coeff_gndf, $
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
                reddening_val=reddening_val, reddening_err=reddening_err

        ; Check that ndrp matches the size of one of the existing inputs
        ; TODO: Assumes all input vectors have the same length!
        nel = intarr(17)
        nel[0] = n_elements(weights_gndf) eq 0 ? 0 : (size(weights_gndf))[1]
        nel[1] = n_elements(mult_poly_coeff_gndf) eq 0 ? 0 : (size(mult_poly_coeff_gndf))[1]
        nel[2] = n_elements(emission_line_kinematics_avg) eq 0 ? 0 : $
                 (size(emission_line_kinematics_avg))[1]
        nel[3] = n_elements(emission_line_kinematics_aer) eq 0 ? 0 : $
                 (size(emission_line_kinematics_aer))[1]
        nel[4] = n_elements(chi2_gndf) eq 0 ? 0 : (size(chi2_gndf))[1]
        nel[5] = n_elements(emission_line_kinematics_ind) eq 0 ? 0 : $
                 (size(emission_line_kinematics_ind))[1]
        nel[6] = n_elements(emission_line_kinematics_ier) eq 0 ? 0 : $
                 (size(emission_line_kinematics_ier))[1]
        nel[7] = n_elements(emission_line_sinst) eq 0 ? 0 : (size(emission_line_sinst))[1]
        nel[8] = n_elements(emission_line_omitted) eq 0 ? 0 : (size(emission_line_omitted))[1]
        nel[9] = n_elements(emission_line_intens) eq 0 ? 0 : (size(emission_line_intens))[1]
        nel[10] = n_elements(emission_line_interr) eq 0 ? 0 : (size(emission_line_interr))[1]
        nel[11] = n_elements(emission_line_fluxes) eq 0 ? 0 : (size(emission_line_fluxes))[1]
        nel[12] = n_elements(emission_line_flxerr) eq 0 ? 0 : (size(emission_line_flxerr))[1]
        nel[13] = n_elements(emission_line_EWidth) eq 0 ? 0 : (size(emission_line_EWidth))[1]
        nel[14] = n_elements(emission_line_EW_err) eq 0 ? 0 : (size(emission_line_EW_err))[1]
        nel[15] = n_elements(reddening_val) eq 0 ? 0 : (size(reddening_val))[1]
        nel[16] = n_elements(reddening_err) eq 0 ? 0 : (size(reddening_err))[1]
        ninp = max(nel)

        for i=0,16 do begin
;           if nel[i] gt 0 and nel[i] ne ninp then $
            if nel[i] gt 0 && nel[i] ne ninp then $
                message, 'Input vectors for SGFIT have different lengths!'
        endfor

        something_to_write = ninp ne 0
END

;-------------------------------------------------------------------------------
; Check the input vectors to write to the STFIT extension.  Returns a flag that
; there is something to write, and the size of the vectors to write.  Will throw
; an error if the input vectors are not the same size
PRO MDAP_WRITE_OUTPUT_EMISSION_LINE_FIT_CHECK_INPUTS, $
                something_to_write, ninp, elo_ew_kinematics_avg=elo_ew_kinematics_avg, $
                elo_ew_kinematics_aer=elo_ew_kinematics_aer, $
                elo_ew_kinematics_ind=elo_ew_kinematics_ind, $
                elo_ew_kinematics_ier=elo_ew_kinematics_ier, elo_ew_sinst=elo_ew_sinst, $
                elo_ew_omitted=elo_ew_omitted, elo_ew_intens=elo_ew_intens, $
                elo_ew_interr=elo_ew_interr, elo_ew_fluxes=elo_ew_fluxes, $
                elo_ew_flxerr=elo_ew_flxerr, elo_ew_EWidth=elo_ew_EWidth, $
                elo_ew_EW_err=elo_ew_EW_err, elo_fb_kinematics_avg=elo_fb_kinematics_avg, $
                elo_fb_kinematics_aer=elo_fb_kinematics_aer, $
                elo_fb_kinematics_ind=elo_fb_kinematics_ind, $
                elo_fb_kinematics_ier=elo_fb_kinematics_ier, elo_fb_sinst=elo_fb_sinst, $
                elo_fb_omitted=elo_fb_omitted, elo_fb_intens=elo_fb_intens, $
                elo_fb_interr=elo_fb_interr, elo_fb_fluxes=elo_fb_fluxes, $
                elo_fb_flxerr=elo_fb_flxerr, elo_fb_EWidth=elo_fb_EWidth, $
                elo_fb_EW_err=elo_fb_EW_err

        ; Check that ndrp matches the size of one of the existing inputs
        ; TODO: Assumes all input vectors have the same length!
        nel = intarr(24)
        nel[0] = n_elements(elo_ew_kinematics_avg) eq 0 ? 0 : (size(elo_ew_kinematics_avg))[1]
        nel[1] = n_elements(elo_ew_kinematics_aer) eq 0 ? 0 : (size(elo_ew_kinematics_aer))[1]
        nel[2] = n_elements(elo_ew_kinematics_ind) eq 0 ? 0 : (size(elo_ew_kinematics_ind))[1]
        nel[3] = n_elements(elo_ew_kinematics_ier) eq 0 ? 0 : (size(elo_ew_kinematics_ier))[1]
        nel[4] = n_elements(elo_ew_sinst) eq 0 ? 0 : (size(elo_ew_sinst))[1]
        nel[5] = n_elements(elo_ew_omitted) eq 0 ? 0 : (size(elo_ew_omitted))[1]
        nel[6] = n_elements(elo_ew_intens) eq 0 ? 0 : (size(elo_ew_intens))[1]
        nel[7] = n_elements(elo_ew_interr) eq 0 ? 0 : (size(elo_ew_interr))[1]
        nel[8] = n_elements(elo_ew_fluxes) eq 0 ? 0 : (size(elo_ew_fluxes))[1]
        nel[9] = n_elements(elo_ew_flxerr) eq 0 ? 0 : (size(elo_ew_flxerr))[1]
        nel[10] = n_elements(elo_ew_EWidth) eq 0 ? 0 : (size(elo_ew_EWidth))[1]
        nel[11] = n_elements(elo_ew_EW_err) eq 0 ? 0 : (size(elo_ew_EW_err))[1]
        nel[12] = n_elements(elo_fb_kinematics_avg) eq 0 ? 0 : (size(elo_fb_kinematics_avg))[1]
        nel[13] = n_elements(elo_fb_kinematics_aer) eq 0 ? 0 : (size(elo_fb_kinematics_aer))[1]
        nel[14] = n_elements(elo_fb_kinematics_ind) eq 0 ? 0 : (size(elo_fb_kinematics_ind))[1]
        nel[15] = n_elements(elo_fb_kinematics_ier) eq 0 ? 0 : (size(elo_fb_kinematics_ier))[1]
        nel[16] = n_elements(elo_fb_sinst) eq 0 ? 0 : (size(elo_fb_sinst))[1]
        nel[17] = n_elements(elo_fb_omitted) eq 0 ? 0 : (size(elo_fb_omitted))[1]
        nel[18] = n_elements(elo_fb_intens) eq 0 ? 0 : (size(elo_fb_intens))[1]
        nel[19] = n_elements(elo_fb_interr) eq 0 ? 0 : (size(elo_fb_interr))[1]
        nel[20] = n_elements(elo_fb_fluxes) eq 0 ? 0 : (size(elo_fb_fluxes))[1]
        nel[21] = n_elements(elo_fb_flxerr) eq 0 ? 0 : (size(elo_fb_flxerr))[1]
        nel[22] = n_elements(elo_fb_EWidth) eq 0 ? 0 : (size(elo_fb_EWidth))[1]
        nel[23] = n_elements(elo_fb_EW_err) eq 0 ? 0 : (size(elo_fb_EW_err))[1]
        ninp = max(nel)

        for i=0,23 do begin
;           if nel[i] gt 0 and nel[i] ne ninp then $
            if nel[i] gt 0 && nel[i] ne ninp then $
                message, 'Input vectors for ELOFIT have different lengths!'
        endfor

        something_to_write = ninp ne 0
END

;-------------------------------------------------------------------------------
; Check the input vectors to write to the STFIT extension.  Returns a flag that
; there is something to write, and the size of the vectors to write.  Will throw
; an error if the input vectors are not the same size
PRO MDAP_WRITE_OUTPUT_SPECTRAL_INDICES_CHECK_INPUTS, $
                something_to_write, ninp, abs_line_indx_omitted=abs_line_indx_omitted, $
                abs_line_indx_val=abs_line_indx_val, abs_line_indx_err=abs_line_indx_err, $
                abs_line_indx_otpl=abs_line_indx_otpl, abs_line_indx_botpl=abs_line_indx_botpl

        ; Check that ndrp matches the size of one of the existing inputs
        ; TODO: Assumes all input vectors have the same length!
        nel = intarr(5)
        nel[0] = n_elements(abs_line_indx_omitted) eq 0 ? 0 : (size(abs_line_indx_omitted))[1]
        nel[1] = n_elements(abs_line_indx_val) eq 0 ? 0 : (size(abs_line_indx_val))[1]
        nel[2] = n_elements(abs_line_indx_err) eq 0 ? 0 : (size(abs_line_indx_err))[1]
        nel[3] = n_elements(abs_line_indx_otpl) eq 0 ? 0 : (size(abs_line_indx_otpl))[1]
        nel[4] = n_elements(abs_line_indx_botpl) eq 0 ? 0 : (size(abs_line_indx_botpl))[1]
        ninp = max(nel)

        for i=0,4 do begin
;           if nel[i] gt 0 and nel[i] ne ninp then $
            if nel[i] gt 0 && nel[i] ne ninp then $
                message, 'Input vectors for SINDX have different lengths!'
        endfor

        something_to_write = ninp ne 0
END

;-------------------------------------------------------------------------------
; Write/update the DRPS extension
;       - Initialize the extension
;       - Check that there are inputs to write
;       - Check that the size is as expected, if provided (ndrp)
;       - Resize (number of rows) the extension, if necessary
;       - Write/update the data
PRO MDAP_WRITE_OUTPUT_UPDATE_DRPS, $
                file, ndrp=ndrp, xpos=xpos, ypos=ypos, fraction_good=fraction_good, $
                min_eq_max=min_eq_max, signal=signal, noise=noise, bin_vreg=bin_vreg, $
                bin_indx=bin_indx, bin_weights=bin_weights, quiet=quiet

;           get_lun, unit
;           print, unit
;           if unit ne 100 then stop
;           free_lun, unit
                

        MDAP_WRITE_OUTPUT_DRPS_INITIALIZE, file         ; Initialize the extension

        ; Check that one of the vectors are input
        MDAP_WRITE_OUTPUT_DRPS_CHECK_INPUTS, something_to_write, ninp, xpos=xpos, ypos=ypos, $
                                             fraction_good=fraction_good, min_eq_max=min_eq_max, $
                                             signal=signal, noise=noise, bin_vreg=bin_vreg, $
                                             bin_indx=bin_indx, bin_weights=bin_weights

        ; If there is nothing to write, return
        if something_to_write eq 0 then begin
            if ~keyword_set(quiet) then $
                print, 'MDAP_WRITE_OUTPUT_UPDATE_DRPS: Nothing to update'
            return
        endif

        ; Check that the length of the vectors matched the expected value, if input
        if n_elements(ndrp) ne 0 then begin
            if ndrp ne ninp then begin
                message, 'Input vectors do not have the expected size! ' + MDAP_STC(ndrp,/integer) $
                         + ' ' + MDAP_STC(ninp,/integer) + ')'
            endif
        endif else $
            ndrp = ninp

        ; Ensure the table has the correct number of rows
        ; TODO: This will fail if the table is longer than ndrp!
        MDAP_WRITE_OUTPUT_TBL_MATCH_SIZE, file, ndrp, 'DRPS'

        ; FROM mdap_set_output_file_cols.pro (included above)
        cols = MDAP_SET_DRPS_COLS()
            
        FXBOPEN, tbl, file, 'DRPS', access='RW'         ; Open the file

        ; Write columns if they were provided
        ; TODO: Is this too inefficient?
        if n_elements(xpos) ne 0 then           FXBWRITM, tbl, [cols[0]], xpos
        if n_elements(ypos) ne 0 then           FXBWRITM, tbl, [cols[1]], ypos
        if n_elements(fraction_good) ne 0 then  FXBWRITM, tbl, [cols[2]], fraction_good
        if n_elements(min_eq_max) ne 0 then     FXBWRITM, tbl, [cols[3]], min_eq_max
        if n_elements(signal) ne 0 then         FXBWRITM, tbl, [cols[4]], signal
        if n_elements(noise) ne 0 then          FXBWRITM, tbl, [cols[5]], noise
        if n_elements(bin_vreg) ne 0 then       FXBWRITM, tbl, [cols[6]], bin_vreg
        if n_elements(bin_indx) ne 0 then       FXBWRITM, tbl, [cols[7]], bin_indx
        if n_elements(bin_weights) ne 0 then    FXBWRITM, tbl, [cols[8]], bin_weights

        FXBFINISH, tbl                                  ; Close the file
;       print, tbl
        free_lun, tbl
END

;-------------------------------------------------------------------------------
; Write/update the BINS extension
;       - Initialize the extension
;       - Check that there are inputs to write
;       - Check that the size is as expected, if provided (ndrp)
;       - Resize (number of rows) the extension, if necessary
;       - Write/update the data
PRO MDAP_WRITE_OUTPUT_UPDATE_BINS, $
                file, nbin=nbin, xbin_rlow=xbin_rlow, ybin_rupp=ybin_rupp, rbin=rbin, $
                bin_area=bin_area, bin_ston=bin_ston, bin_n=bin_n, bin_flag=bin_flag, quiet=quiet

        MDAP_WRITE_OUTPUT_BINS_INITIALIZE, file         ; Initialize the extension

        ; Check that one of the vectors are input
        MDAP_WRITE_OUTPUT_BINS_CHECK_INPUTS, something_to_write, ninp, xbin_rlow=xbin_rlow, $
                                             ybin_rupp=ybin_rupp, rbin=rbin, bin_area=bin_area, $
                                             bin_ston=bin_ston, bin_n=bin_n, bin_flag=bin_flag

        ; If there is nothing to write, return
        if something_to_write eq 0 then begin
            if ~keyword_set(quiet) then $
                print, 'MDAP_WRITE_OUTPUT_UPDATE_BINS: Nothing to update'
            return
        endif

        ; Check that the length of the vectors matched the expected value, if input
        if n_elements(nbin) ne 0 then begin
            if nbin ne ninp then begin
                message, 'Input vectors do not have the expected size! ' + MDAP_STC(nbin,/integer) $
                         + ' ' + MDAP_STC(ninp,/integer) + ')'
            endif
        endif else $
            nbin = ninp

        ; Ensure the table has the correct number of rows
        ; TODO: This will fail if the table is longer than nbin!
        MDAP_WRITE_OUTPUT_TBL_MATCH_SIZE, file, nbin, 'BINS'

        ; FROM mdap_set_output_file_cols.pro (included above)
        cols = MDAP_SET_BINS_COLS()
            
        FXBOPEN, tbl, file, 'BINS', access='RW'         ; Open the file

        ; Write columns if they were provided
        ; TODO: Is this too inefficient?
        if n_elements(xbin_rlow) ne 0 then FXBWRITM, tbl, [cols[0]], xbin_rlow
        if n_elements(ybin_rupp) ne 0 then FXBWRITM, tbl, [cols[1]], ybin_rupp
        if n_elements(rbin) ne 0 then      FXBWRITM, tbl, [cols[2]], rbin
        if n_elements(bin_area) ne 0 then  FXBWRITM, tbl, [cols[3]], bin_area
        if n_elements(bin_ston) ne 0 then  FXBWRITM, tbl, [cols[4]], bin_ston
        if n_elements(bin_n) ne 0 then     FXBWRITM, tbl, [cols[5]], bin_n
        if n_elements(bin_flag) ne 0 then  FXBWRITM, tbl, [cols[6]], bin_flag

        FXBFINISH, tbl                                  ; Close the file
;       print, tbl
        free_lun, tbl
END

;-------------------------------------------------------------------------------
; Write/update the ELPAR extension
;       - Initialize the extension
;       - Check that there are inputs to write
;       - Check that the size is as expected, if provided (ndrp)
;       - Resize (number of rows) the extension, if necessary
;       - Write/update the data
PRO MDAP_WRITE_OUTPUT_UPDATE_ELPAR, $
                file, eml_par=eml_par, quiet=quiet

        MDAP_WRITE_OUTPUT_ELPAR_INITIALIZE, file        ; Initialize the extension

        ; Check that there is something to write
        neml = n_elements(eml_par)                      ; Number of elements to output
        if neml eq 0 then begin
            if ~keyword_set(quiet) then $
                print, 'MDAP_WRITE_OUTPUT_UPDATE_ELPAR: Nothing to update'
            return
        endif

        ; Ensure the table has the correct number of rows
        ; TODO: This will fail if the table is longer than ndrp!
        MDAP_WRITE_OUTPUT_TBL_MATCH_SIZE, file, neml, 'ELPAR'

        FXBOPEN, tbl, file, 'ELPAR', hdr, access='RW'   ; Open the extension

        elname = MDAP_SET_EMISSION_LINE_NAME(eml_par)   ; Set the emission-line identifiers
        elwave = eml_par[*].lambda                      ; Copy the wavelengths

        ; TODO: FXBWRITM does not work with TiedType if it is only one character.  Why?
        eltied = intarr(neml, /nozero)                  ; Index of tied line (-1 if not tied)
        elttyp = strarr(neml)                           ; Type of tied (-,t,v,s)
        eldoub = intarr(neml, /nozero)                  ; Index of doublet line (-1 if not doublet)
        for i=0,neml-1 do begin

            ; Check if any of the kinematics are tied
            tt=strmid(eml_par[i].fit,0,1)
;           if tt eq 't' or tt eq 'v' or tt eq 's' then begin
            if tt eq 't' || tt eq 'v' || tt eq 's' then begin
;               j_refline = where(eml_par[*].i eq fix(strmid(eml_par[i].fit,1)))
;               eltied[i] = j_refline[0]                ; Set to -1 if not found!
                j_refline = where(eml_par[*].i eq fix(strmid(eml_par[i].fit,1)), count)
                eltied[i] = count eq 0 ? -1 : j_refline[0]                ; Set to -1 if not found!
                ; Set the type of kinematic tying
                elttyp[i] = MDAP_SET_EMISSION_LINE_TIED(tt)
            endif else begin
                eltied[i] = -1                          ; Set untied values
                elttyp[i] = MDAP_SET_EMISSION_LINE_TIED('-')
            endelse

            ; Check if any of the lines are doubles
            if strmid(eml_par[i].kind,0,1) eq 'd' then begin
;               j_refline = where(eml_par[*].i eq fix(strmid(eml_par[i].kind,1)))
;               eldoub[i] = j_refline[0]                ; Set to -1 if not found!
                j_refline = where(eml_par[*].i eq fix(strmid(eml_par[i].kind,1)), count)
                eldoub[i] = count eq 0 ? -1 : j_refline[0]                ; Set to -1 if not found!
            endif else $
                eldoub[i] = -1
        endfor
        
;       for i=0,neml-1 do $
;           print, '+'+elttyp[i]+'+', n_elements(elttyp[i])
;       print, hdr

        ; FROM mdap_set_output_file_cols.pro (included above)
        cols = MDAP_SET_ELPAR_COLS()
            
        ; Write ALL the columns simultaneously
        FXBWRITM, tbl, [ cols[0], cols[1], cols[2], cols[3], cols[4] ], elname, elwave, eltied, $
                  elttyp, eldoub

        FXBFINISH, tbl                                  ; Close the file
;       print, tbl
        free_lun, tbl
END

;-------------------------------------------------------------------------------
; Write/update the STFIT extension
PRO MDAP_WRITE_OUTPUT_UPDATE_STELLAR_KIN, $
                file, nbin=nbin, weights_ppxf=weights_ppxf, $
                add_poly_coeff_ppxf=add_poly_coeff_ppxf, $
                mult_poly_coeff_ppxf=mult_poly_coeff_ppxf, $
                stellar_kinematics_fit=stellar_kinematics_fit, $
                stellar_kinematics_err=stellar_kinematics_err, chi2_ppxf=chi2_ppxf, $
                quiet=quiet


        ; Initialize the extension
        MDAP_WRITE_OUTPUT_STELLAR_KIN_INITIALIZE, file, weights_ppxf=weights_ppxf, $
                                                  add_poly_coeff_ppxf=add_poly_coeff_ppxf, $
                                                  mult_poly_coeff_ppxf=mult_poly_coeff_ppxf, $
                                                  stellar_kinematics_fit=stellar_kinematics_fit, $
                                                  stellar_kinematics_err=stellar_kinematics_err


        ; Check that one of the vectors are input and that they have the same length
        MDAP_WRITE_OUTPUT_STELLAR_KIN_CHECK_INPUTS, something_to_write, ninp, $
                                                    weights_ppxf=weights_ppxf, $
                                                    add_poly_coeff_ppxf=add_poly_coeff_ppxf, $
                                                    mult_poly_coeff_ppxf=mult_poly_coeff_ppxf, $
                                                    stellar_kinematics_fit=stellar_kinematics_fit, $
                                                    stellar_kinematics_err=stellar_kinematics_err, $
                                                    chi2_ppxf=chi2_ppxf

        ; If there is nothing to write, return
        if something_to_write eq 0 then begin
            if ~keyword_set(quiet) then $
                print, 'MDAP_WRITE_OUTPUT_UPDATE_STELLAR_KIN: Nothing to update'
            return
        endif

        ; Check that the length of the vectors matched the expected value, if input
        if n_elements(nbin) ne 0 then begin
            if nbin ne ninp then $
                message, 'Input vectors do not have the expected size!'
        endif else $
            nbin = ninp

        ; Ensure the table has the correct number of rows
        ; TODO: This will fail if the table is longer than ndrp!
        MDAP_WRITE_OUTPUT_TBL_MATCH_SIZE, file, nbin, 'STFIT'

        ; FROM mdap_set_output_file_cols.pro (included above)
        cols = MDAP_SET_STFIT_COLS()
            
        FXBOPEN, tbl, file, 'STFIT', access='RW'        ; Open the file

        ; Write columns if they were provided; column/row ordering has to be flipped
        ; TODO: Is this too inefficient?
        ; TODO: Need to redo column/row reordering
        if n_elements(weights_ppxf) ne 0 then $
            FXBWRITM, tbl, [cols[0]], transpose(weights_ppxf)
        if n_elements(add_poly_coeff_ppxf) ne 0 then $
            FXBWRITM, tbl, [cols[1]], transpose(add_poly_coeff_ppxf)
        if n_elements(mult_poly_coeff_ppxf) ne 0 then $
            FXBWRITM, tbl, [cols[2]], transpose(mult_poly_coeff_ppxf)
        if n_elements(stellar_kinematics_fit) ne 0 then $
            FXBWRITM, tbl, [cols[3]], transpose(stellar_kinematics_fit)
        if n_elements(stellar_kinematics_err) ne 0 then $
            FXBWRITM, tbl, [cols[4]], transpose(stellar_kinematics_err)
        if n_elements(chi2_ppxf) ne 0 then $
            FXBWRITM, tbl, [cols[5]], chi2_ppxf

        FXBFINISH, tbl                                  ; Close the file
;       print, tbl
        free_lun, tbl

;       FXBOPEN, tbl, file, 'STFIT'
;       FXBREAD, tbl, testw, 'TPLW'
;       print, size(testw)
;       print, testw[*,0]
;       FXBCLOSE, tbl
;       free_lun, tbl
;       stop
END

;-------------------------------------------------------------------------------
; Write/update the SGFIT extension
; TODO: Expects reddening and reddening error to have the same size!
PRO MDAP_WRITE_OUTPUT_UPDATE_STAR_AND_GAS_FIT, $
                file, nbin=nbin, eml_par=eml_par, weights_gndf=weights_gndf, $
                mult_poly_coeff_gndf=mult_poly_coeff_gndf, $
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
                reddening_val=reddening_val, reddening_err=reddening_err, quiet=quiet

        ; Initialize the extension
        MDAP_WRITE_OUTPUT_STAR_AND_GAS_FIT_INITIALIZE, file, eml_par=eml_par, $
                                        weights_gndf=weights_gndf, $
                                        mult_poly_coeff_gndf=mult_poly_coeff_gndf, $
                                        emission_line_kinematics_avg=emission_line_kinematics_avg, $
                                        emission_line_kinematics_aer=emission_line_kinematics_aer, $
                                        emission_line_kinematics_ind=emission_line_kinematics_ind, $
                                        emission_line_kinematics_ier=emission_line_kinematics_ier, $
                                        emission_line_sinst=emission_line_sinst, $
                                        emission_line_omitted=emission_line_omitted, $
                                        emission_line_intens=emission_line_intens, $
                                        emission_line_interr=emission_line_interr, $
                                        emission_line_fluxes=emission_line_fluxes, $
                                        emission_line_flxerr=emission_line_flxerr, $
                                        emission_line_EWidth=emission_line_EWidth, $
                                        emission_line_EW_err=emission_line_EW_err, $
                                        reddening_val=reddening_val, $
                                        reddening_err=reddening_err

        ; Check that one of the vectors are input and that they have the same length
        MDAP_WRITE_OUTPUT_STAR_AND_GAS_FIT_CHECK_INPUTS, something_to_write, ninp, $
                                        weights_gndf=weights_gndf, $
                                        mult_poly_coeff_gndf=mult_poly_coeff_gndf, $
                                        emission_line_kinematics_avg=emission_line_kinematics_avg, $
                                        emission_line_kinematics_aer=emission_line_kinematics_aer, $
                                        chi2_gndf=chi2_gndf, $
                                        emission_line_kinematics_ind=emission_line_kinematics_ind, $
                                        emission_line_kinematics_ier=emission_line_kinematics_ier, $
                                        emission_line_sinst=emission_line_sinst, $
                                        emission_line_omitted=emission_line_omitted, $
                                        emission_line_intens=emission_line_intens, $
                                        emission_line_interr=emission_line_interr, $
                                        emission_line_fluxes=emission_line_fluxes, $
                                        emission_line_flxerr=emission_line_flxerr, $
                                        emission_line_EWidth=emission_line_EWidth, $
                                        emission_line_EW_err=emission_line_EW_err, $
                                        reddening_val=reddening_val, reddening_err=reddening_err

        ; If there is nothing to write, return
        if something_to_write eq 0 then begin
            if ~keyword_set(quiet) then $
                print, 'MDAP_WRITE_OUTPUT_UPDATE_STAR_AND_GAS_FIT: Nothing to update'
            return
        endif

        ; Check that the length of the vectors matched the expected value, if input
        if n_elements(nbin) ne 0 then begin
            if nbin ne ninp then $
                message, 'Input vectors do not have the expected size!'
        endif else $
            nbin = ninp

;       fits_info,file
;       stop

        ; Ensure the table has the correct number of rows
        ; TODO: This will fail if the table is longer than nbin!
        MDAP_WRITE_OUTPUT_TBL_MATCH_SIZE, file, nbin, 'SGFIT'

        ; FROM mdap_set_output_file_cols.pro (included above)
        cols = MDAP_SET_SGFIT_COLS()
            
;       fits_info, file
;       stop

        FXBOPEN, tbl, file, 'SGFIT', bth, access='RW'   ; Open the file

        ; Write columns if they were provided; column/row ordering has to be adjusted
        ; TODO: Is this too inefficient?
        ; TODO: Need to redo column/row reordering
        if n_elements(weights_gndf) ne 0 then $
            FXBWRITM, tbl, [cols[0]], transpose(weights_gndf)
        if n_elements(mult_poly_coeff_gndf) ne 0 then $
            FXBWRITM, tbl, [cols[1]], transpose(mult_poly_coeff_gndf)
        if n_elements(emission_line_kinematics_avg) ne 0 then $
            FXBWRITM, tbl, [cols[2]], transpose(emission_line_kinematics_avg)
        if n_elements(emission_line_kinematics_aer) ne 0 then $
            FXBWRITM, tbl, [cols[3]], transpose(emission_line_kinematics_aer)
        if n_elements(chi2_gndf) ne 0 then $
            FXBWRITM, tbl, [cols[4]], chi2_gndf
        if n_elements(reddening_val) ne 0 then $
            FXBWRITM, tbl, [cols[5]], transpose(reddening_val)
        if n_elements(reddening_err) ne 0 then $
            FXBWRITM, tbl, [cols[6]], transpose(reddening_err)
        if n_elements(emission_line_omitted) ne 0 then $
            FXBWRITM, tbl, [cols[7]], emission_line_omitted
        if n_elements(emission_line_intens) ne 0 then $
            FXBWRITM, tbl, [cols[8]], transpose(emission_line_intens)
        if n_elements(emission_line_interr) ne 0 then $
            FXBWRITM, tbl, [cols[9]], transpose(emission_line_interr)
        if n_elements(emission_line_kinematics_ind) ne 0 then $
            FXBWRITM, tbl, [cols[10]], transpose(emission_line_kinematics_ind, [1, 2, 0])
        if n_elements(emission_line_kinematics_ier) ne 0 then $
            FXBWRITM, tbl, [cols[11]], transpose(emission_line_kinematics_ier, [1, 2, 0])
        if n_elements(emission_line_sinst) ne 0 then $
            FXBWRITM, tbl, [cols[12]], transpose(emission_line_sinst)
        if n_elements(emission_line_fluxes) ne 0 then $
            FXBWRITM, tbl, [cols[13]], transpose(emission_line_fluxes)
        if n_elements(emission_line_flxerr) ne 0 then $
            FXBWRITM, tbl, [cols[14]], transpose(emission_line_flxerr)
        if n_elements(emission_line_EWidth) ne 0 then $
            FXBWRITM, tbl, [cols[15]], transpose(emission_line_EWidth)
        if n_elements(emission_line_EW_err) ne 0 then $
            FXBWRITM, tbl, [cols[16]], transpose(emission_line_EW_err)

        FXBFINISH, tbl                                  ; Close the file
;       print, tbl
        free_lun, tbl

;       FXBOPEN, tbl, file, 'SGFIT'
;       FXBREAD, tbl, testw, 'IKIN'
;       print, size(testw)
;       print, 'read:', testw[*,*,0]
;       print, 'data:', emission_line_kinematics_ind[0,*,*]
;       FXBCLOSE, tbl
;       free_lun, tbl
;       stop
END

;-------------------------------------------------------------------------------
; Write/update the ELOPAR extension
;       - Initialize the extension
;       - Check that there are inputs to write
;       - Write/update the data
PRO MDAP_WRITE_OUTPUT_UPDATE_ELOPAR, $
                file, emlo_par=emlo_par, quiet=quiet

        MDAP_WRITE_OUTPUT_ELOPAR_INITIALIZE, file       ; Initialize the extension

        ; Check that there is something to write
        neml = n_elements(emlo_par)                      ; Number of elements to output
        if neml eq 0 then begin
            if ~keyword_set(quiet) then $
                print, 'MDAP_WRITE_OUTPUT_UPDATE_ELOPAR: Nothing to update'
            return
        endif

        ; Ensure the table has the correct number of rows
        ; TODO: This will fail if the table is longer than ndrp!
        MDAP_WRITE_OUTPUT_TBL_MATCH_SIZE, file, neml, 'ELOPAR'

        FXBOPEN, tbl, file, 'ELOPAR', hdr, access='RW'   ; Open the extension

        elname = MDAP_SET_EMISSION_LINE_NAME(emlo_par)   ; Set the emission-line identifiers
        elwave = emlo_par[*].lambda                      ; Copy the wavelengths

        ; FROM mdap_set_output_file_cols.pro (included above)
        cols = MDAP_SET_ELOPAR_COLS()
            
        ; Write ALL the columns simultaneously
        FXBWRITM, tbl, [ cols[0], cols[1] ], elname, elwave

        FXBFINISH, tbl                                  ; Close the file
;       print, tbl
        free_lun, tbl
END

;-------------------------------------------------------------------------------
; Write/update the ELOFIT extension
; TODO: Expects reddening and reddening error to have the same size!
PRO MDAP_WRITE_OUTPUT_UPDATE_EMISSION_LINE_FIT, $
                file, nbin=nbin, emlo_par=emlo_par, elo_ew_kinematics_avg=elo_ew_kinematics_avg, $
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
                elo_fb_EW_err=elo_fb_EW_err

        ; Initialize the extension
        MDAP_WRITE_OUTPUT_EMISSION_LINE_FIT_INITIALIZE, file, emlo_par=emlo_par, $
                        elo_ew_kinematics_avg=elo_ew_kinematics_avg, $
                        elo_ew_kinematics_aer=elo_ew_kinematics_aer, $
                        elo_ew_kinematics_ind=elo_ew_kinematics_ind, $
                        elo_ew_kinematics_ier=elo_ew_kinematics_ier, elo_ew_sinst=elo_ew_sinst, $
                        elo_ew_omitted=elo_ew_omitted, elo_ew_intens=elo_ew_intens, $
                        elo_ew_interr=elo_ew_interr, elo_ew_fluxes=elo_ew_fluxes, $
                        elo_ew_flxerr=elo_ew_flxerr, elo_ew_EWidth=elo_ew_EWidth, $
                        elo_ew_EW_err=elo_ew_EW_err, elo_fb_kinematics_avg=elo_fb_kinematics_avg, $
                        elo_fb_kinematics_aer=elo_fb_kinematics_aer, $
                        elo_fb_kinematics_ind=elo_fb_kinematics_ind, $
                        elo_fb_kinematics_ier=elo_fb_kinematics_ier, elo_fb_sinst=elo_fb_sinst, $
                        elo_fb_omitted=elo_fb_omitted, elo_fb_intens=elo_fb_intens, $
                        elo_fb_interr=elo_fb_interr, elo_fb_fluxes=elo_fb_fluxes, $
                        elo_fb_flxerr=elo_fb_flxerr, elo_fb_EWidth=elo_fb_EWidth, $
                        elo_fb_EW_err=elo_fb_EW_err

        ; Check that one of the vectors are input and that they have the same length
        MDAP_WRITE_OUTPUT_EMISSION_LINE_FIT_CHECK_INPUTS, something_to_write, ninp, $
                        elo_ew_kinematics_avg=elo_ew_kinematics_avg, $
                        elo_ew_kinematics_aer=elo_ew_kinematics_aer, $
                        elo_ew_kinematics_ind=elo_ew_kinematics_ind, $
                        elo_ew_kinematics_ier=elo_ew_kinematics_ier, elo_ew_sinst=elo_ew_sinst, $
                        elo_ew_omitted=elo_ew_omitted, elo_ew_intens=elo_ew_intens, $
                        elo_ew_interr=elo_ew_interr, elo_ew_fluxes=elo_ew_fluxes, $
                        elo_ew_flxerr=elo_ew_flxerr, elo_ew_EWidth=elo_ew_EWidth, $
                        elo_ew_EW_err=elo_ew_EW_err, elo_fb_kinematics_avg=elo_fb_kinematics_avg, $
                        elo_fb_kinematics_aer=elo_fb_kinematics_aer, $
                        elo_fb_kinematics_ind=elo_fb_kinematics_ind, $
                        elo_fb_kinematics_ier=elo_fb_kinematics_ier, elo_fb_sinst=elo_fb_sinst, $
                        elo_fb_omitted=elo_fb_omitted, elo_fb_intens=elo_fb_intens, $
                        elo_fb_interr=elo_fb_interr, elo_fb_fluxes=elo_fb_fluxes, $
                        elo_fb_flxerr=elo_fb_flxerr, elo_fb_EWidth=elo_fb_EWidth, $
                        elo_fb_EW_err=elo_fb_EW_err

        ; If there is nothing to write, return
        if something_to_write eq 0 then begin
            if ~keyword_set(quiet) then $
                print, 'MDAP_WRITE_OUTPUT_UPDATE_EMISSION_LINE_FIT: Nothing to update'
            return
        endif

        ; Check that the length of the vectors matched the expected value, if input
        if n_elements(nbin) ne 0 then begin
            if nbin ne ninp then $
                message, 'Input vectors do not have the expected size!'
        endif else $
            nbin = ninp

;       fits_info,file
;       stop

        ; Ensure the table has the correct number of rows
        ; TODO: This will fail if the table is longer than nbin!
        MDAP_WRITE_OUTPUT_TBL_MATCH_SIZE, file, nbin, 'ELOFIT'

        ; FROM mdap_set_output_file_cols.pro (included above)
        cols = MDAP_SET_ELOFIT_COLS()
            
;       fits_info, file
;       stop

        FXBOPEN, tbl, file, 'ELOFIT', bth, access='RW'   ; Open the file

        ; Write columns if they were provided; column/row ordering has to be adjusted
        ; TODO: Is this too inefficient?
        ; TODO: Need to redo column/row reordering
        if n_elements(elo_ew_kinematics_avg) ne 0 then $
            FXBWRITM, tbl, [cols[0]], transpose(elo_ew_kinematics_avg)
        if n_elements(elo_ew_kinematics_aer) ne 0 then $
            FXBWRITM, tbl, [cols[1]], transpose(elo_ew_kinematics_aer)
        if n_elements(elo_ew_omitted) ne 0 then $
            FXBWRITM, tbl, [cols[2]], elo_ew_omitted
        if n_elements(elo_ew_intens) ne 0 then $
            FXBWRITM, tbl, [cols[3]], transpose(elo_ew_intens)
        if n_elements(elo_ew_interr) ne 0 then $
            FXBWRITM, tbl, [cols[4]], transpose(elo_ew_interr)
        if n_elements(elo_ew_kinematics_ind) ne 0 then $
            FXBWRITM, tbl, [cols[5]], transpose(elo_ew_kinematics_ind, [1, 2, 0])
        if n_elements(elo_ew_kinematics_ier) ne 0 then $
            FXBWRITM, tbl, [cols[6]], transpose(elo_ew_kinematics_ier, [1, 2, 0])
        if n_elements(elo_ew_sinst) ne 0 then $
            FXBWRITM, tbl, [cols[7]], transpose(elo_ew_sinst)
        if n_elements(elo_ew_fluxes) ne 0 then $
            FXBWRITM, tbl, [cols[8]], transpose(elo_ew_fluxes)
        if n_elements(elo_ew_flxerr) ne 0 then $
            FXBWRITM, tbl, [cols[9]], transpose(elo_ew_flxerr)
        if n_elements(elo_ew_EWidth) ne 0 then $
            FXBWRITM, tbl, [cols[10]], transpose(elo_ew_EWidth)
        if n_elements(elo_ew_EW_err) ne 0 then $
            FXBWRITM, tbl, [cols[11]], transpose(elo_ew_EW_err)
        if n_elements(elo_fb_kinematics_avg) ne 0 then $
            FXBWRITM, tbl, [cols[12]], transpose(elo_fb_kinematics_avg)
        if n_elements(elo_fb_kinematics_aer) ne 0 then $
            FXBWRITM, tbl, [cols[13]], transpose(elo_fb_kinematics_aer)
        if n_elements(elo_fb_omitted) ne 0 then $
            FXBWRITM, tbl, [cols[14]], elo_fb_omitted
        if n_elements(elo_fb_intens) ne 0 then $
            FXBWRITM, tbl, [cols[15]], transpose(elo_fb_intens)
        if n_elements(elo_fb_interr) ne 0 then $
            FXBWRITM, tbl, [cols[16]], transpose(elo_fb_interr)
        if n_elements(elo_fb_kinematics_ind) ne 0 then $
            FXBWRITM, tbl, [cols[17]], transpose(elo_fb_kinematics_ind, [1, 2, 0])
        if n_elements(elo_fb_kinematics_ier) ne 0 then $
            FXBWRITM, tbl, [cols[18]], transpose(elo_fb_kinematics_ier, [1, 2, 0])
        if n_elements(elo_fb_sinst) ne 0 then $
            FXBWRITM, tbl, [cols[19]], transpose(elo_fb_sinst)
        if n_elements(elo_fb_fluxes) ne 0 then $
            FXBWRITM, tbl, [cols[20]], transpose(elo_fb_fluxes)
        if n_elements(elo_fb_flxerr) ne 0 then $
            FXBWRITM, tbl, [cols[21]], transpose(elo_fb_flxerr)
        if n_elements(elo_fb_EWidth) ne 0 then $
            FXBWRITM, tbl, [cols[22]], transpose(elo_fb_EWidth)
        if n_elements(elo_fb_EW_err) ne 0 then $
            FXBWRITM, tbl, [cols[23]], transpose(elo_fb_EW_err)

        FXBFINISH, tbl                                  ; Close the file
;       print, tbl
        free_lun, tbl

END

;-------------------------------------------------------------------------------
; Write/update the SIPAR extension
;       - Initialize the extension
;       - Check that there are inputs to write
;       - Check that the size is as expected, if provided (ndrp)
;       - Resize (number of rows) the extension, if necessary
;       - Write/update the data
PRO MDAP_WRITE_OUTPUT_UPDATE_SIPAR, $
                file, abs_par=abs_par, quiet=quiet

        MDAP_WRITE_OUTPUT_SIPAR_INITIALIZE, file        ; Initialize the extension

        ; Check that there is something to write
        nabs = n_elements(abs_par)                      ; Number of elements to output
        if nabs eq 0 then begin
            if ~keyword_set(quiet) then $
                print, 'MDAP_WRITE_OUTPUT_UPDATE_SIPAR: Nothing to update'
            return
        endif

        ; Ensure the table has the correct number of rows
        ; TODO: This will fail if the table is longer than ndrp!
        MDAP_WRITE_OUTPUT_TBL_MATCH_SIZE, file, nabs, 'SIPAR'

        ; FROM mdap_set_output_file_cols.pro (included above)
        cols = MDAP_SET_SIPAR_COLS()

        FXBOPEN, tbl, file, 'SIPAR', hdr, access='RW'   ; Open the extension

        siname = MDAP_SET_SPECTRAL_INDEX_NAME(abs_par)  ; Set the emission-line identifiers
        siunit = abs_par[*].units                       ; Copy the unit

        ; Copy the band definitions
        passband = dblarr(2,nabs)                       ; Initialize the passband array
        blueband = dblarr(2,nabs)                       ; Initialize the blue continuum band array
        redband = dblarr(2,nabs)                        ; Initialize the red continuum band array
        for i=0,nabs-1 do begin
            passband[*,i] = abs_par[i].passband         ; Copy the passband
            blueband[*,i] = abs_par[i].blue_cont        ; Copy the blue continuum band
            redband[*,i] = abs_par[i].red_cont          ; Copy the red continuum band
        endfor

        ; Write ALL the columns simultaneously
        FXBWRITM, tbl, [ cols[0], cols[1], cols[2], cols[3], cols[4] ], siname, passband, $
                  blueband, redband, siunit

        FXBFINISH, tbl                                  ; Close the file
;       print, tbl
        free_lun, tbl
END

;-------------------------------------------------------------------------------
; Write/update the SINDX extension
PRO MDAP_WRITE_OUTPUT_UPDATE_SPECTRAL_INDICES, $
                file, nbin=nbin, abs_par=abs_par, abs_line_indx_omitted=abs_line_indx_omitted, $
                abs_line_indx_val=abs_line_indx_val, abs_line_indx_err=abs_line_indx_err, $
                abs_line_indx_otpl=abs_line_indx_otpl, abs_line_indx_botpl=abs_line_indx_botpl

        ; Initialize the extension
        MDAP_WRITE_OUTPUT_SPECTRAL_INDICES_INITIALIZE, file, abs_par=abs_par, $
                                                      abs_line_indx_omitted=abs_line_indx_omitted, $
                                                      abs_line_indx_val=abs_line_indx_val, $
                                                      abs_line_indx_err=abs_line_indx_err, $
                                                      abs_line_indx_otpl=abs_line_indx_otpl, $
                                                      abs_line_indx_botpl=abs_line_indx_botpl

        ; Check that one of the vectors are input and that they have the same length
        MDAP_WRITE_OUTPUT_SPECTRAL_INDICES_CHECK_INPUTS, something_to_write, ninp, $
                                                      abs_line_indx_omitted=abs_line_indx_omitted, $
                                                      abs_line_indx_val=abs_line_indx_val, $
                                                      abs_line_indx_err=abs_line_indx_err, $
                                                      abs_line_indx_otpl=abs_line_indx_otpl, $
                                                      abs_line_indx_botpl=abs_line_indx_botpl

        ; If there is nothing to write, return
        if something_to_write eq 0 then begin
            if ~keyword_set(quiet) then $
                print, 'MDAP_WRITE_OUTPUT_UPDATE_SPECTRAL_INDICES: Nothing to update'
            return
        endif

        ; Check that the length of the vectors matched the expected value, if input
        if n_elements(nbin) ne 0 then begin
            if nbin ne ninp then $
                message, 'Input vectors do not have the expected size!'
        endif else $
            nbin = ninp

;       fits_info,file
;       stop

        ; Ensure the table has the correct number of rows
        ; TODO: This will fail if the table is longer than nbin!
        MDAP_WRITE_OUTPUT_TBL_MATCH_SIZE, file, nbin, 'SINDX'

        ; FROM mdap_set_output_file_cols.pro (included above)
        cols = MDAP_SET_SINDX_COLS()

;       fits_info, file
;       stop

        FXBOPEN, tbl, file, 'SINDX', bth, access='RW'   ; Open the file

        ; Write columns if they were provided; column/row ordering has to be adjusted
        ; TODO: Is this too inefficient?
        ; TODO: Need to redo column/row reordering
        if n_elements(abs_line_indx_omitted) ne 0 then $
            FXBWRITM, tbl, [cols[0]], transpose(abs_line_indx_omitted)
        if n_elements(abs_line_indx_val) ne 0 then $
            FXBWRITM, tbl, [cols[1]], transpose(abs_line_indx_val)
        if n_elements(abs_line_indx_err) ne 0 then $
            FXBWRITM, tbl, [cols[2]], transpose(abs_line_indx_err)
        if n_elements(abs_line_indx_otpl) ne 0 then $
            FXBWRITM, tbl, [cols[3]], transpose(abs_line_indx_otpl)
        if n_elements(abs_line_indx_botpl) ne 0 then $
            FXBWRITM, tbl, [cols[4]], transpose(abs_line_indx_botpl)

        FXBFINISH, tbl                                  ; Close the file
;       print, tbl
        free_lun, tbl

END


;-------------------------------------------------------------------------------
PRO MDAP_WRITE_OUTPUT_UPDATE_IMAGE, $
                file, exten, data, comment

        if MDAP_CHECK_EXTENSION_EXISTS(file, exten) eq 0 then begin
            MKHDR, hdr, data, /image                    ; Make a minimal image extension header
            SXADDPAR, hdr, 'EXTNAME', exten, comment    ; Set the extension name
            WRITEFITS, file, data, hdr, /append         ; Write the data to an appended extension
        endif else begin
            MKHDR, hdr, data, /image                    ; Make a minimal image extension header
            SXADDPAR, hdr, 'EXTNAME', exten, comment    ; Set the extension name
            MODFITS, file, data, hdr, extname=exten     ; Extension already exists, modify it
        endelse
END

;-------------------------------------------------------------------------------
PRO MDAP_WRITE_OUTPUT_BLANK_IMAGE_EXTENSION, $
                file, exten, comment

        if MDAP_CHECK_EXTENSION_EXISTS(file, exten) eq 0 then begin
            dummy = dblarr(1)                           ; Set a dummy array
            MKHDR, hdr, dummy, /image                   ; Make a minimal image extension header
            SXADDPAR, hdr, 'EXTNAME', exten, comment    ; Set the extension name
            WRITEFITS, file, dummy, hdr, /append        ; Write the data to an appended extension
        endif 
END

;-------------------------------------------------------------------------------
; ORDER OF THE OPERATIONS IS IMPORTANT
PRO MDAP_WRITE_OUTPUT, $
                file, header=header, dx=dx, dy=dy, w_range_sn=w_range_sn, xpos=xpos, ypos=ypos, $
                fraction_good=fraction_good, min_eq_max=min_eq_max, signal=signal, noise=noise, $
                bin_par=bin_par, threshold_ston_bin=threshold_ston_bin, bin_vreg=bin_vreg, $
                bin_indx=bin_indx, bin_weights=bin_weights, wave=wave, sres=sres, $
                bin_flux=bin_flux, bin_ivar=bin_ivar, bin_mask=bin_mask, xbin_rlow=xbin_rlow, $
                ybin_rupp=ybin_rupp, rbin=rbin, bin_area=bin_area, bin_ston=bin_ston, bin_n=bin_n, $
                bin_flag=bin_flag, w_range_analysis=w_range_analysis, $
                threshold_ston_analysis=threshold_ston_analysis, tpl_library_key=tpl_library_key, $
                ems_line_key=ems_line_key, eml_par=eml_par, analysis_par=analysis_par, $
                weights_ppxf=weights_ppxf, add_poly_coeff_ppxf=add_poly_coeff_ppxf, $
                mult_poly_coeff_ppxf=mult_poly_coeff_ppxf, $
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
                eml_model=eml_model, emlo_par=emlo_par, $
                elo_ew_kinematics_avg=elo_ew_kinematics_avg, $
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
                si_broad_optimal_template=si_broad_optimal_template, read_header=read_header, $
                quiet=quiet

;           uu = 1
;           print, uu
;           get_lun, unit
;           print, unit
;           if unit ne 100 then stop
;           free_lun, unit

;       if file_test(file) eq 1 then $
;           fits_info, file

        ; Check the number of DRP spectra
        ndrp = MDAP_WRITE_OUTPUT_NUMBER_OF_DRPS(xpos=xpos, ypos=ypos, fraction_good=fraction_good, $
                                                min_eq_max=min_eq_max, signal=signal, noise=noise, $
                                                bin_vreg=bin_vreg, bin_indx=bin_indx, $
                                                bin_weights=bin_weights)
        if ndrp eq 0 then begin
            tmp=temporary(ndrp)
        endif else begin
            if ~keyword_set(quiet) then $
                print, 'Number of DRP spectra: ', ndrp
        endelse

        ; Check the number of bins
        nbin = MDAP_WRITE_OUTPUT_NUMBER_OF_BINS(bin_flux=bin_flux, bin_ivar=bin_ivar, $
                                                bin_mask=bin_mask, xbin_rlow=xbin_rlow, $
                                                ybin_rupp=ybin_rupp, rbin=rbin, bin_area=bin_area, $
                                                bin_ston=bin_ston, bin_n=bin_n, bin_flag=bin_flag)
        if nbin eq 0 then begin
            tmp=temporary(nbin)
        endif else begin
            if ~keyword_set(quiet) then $
                print, 'Number of BIN spectra: ', nbin
        endelse

        ; Read the existing header
        ; TODO: This will overwrite any input header!
        if keyword_set(read_header) then begin
            if file_test(file) eq 0 then begin
                print, 'Output file does not exist; cannot read existing header. Ignoring keyword.'
            endif else $
                header = headfits(file, exten=0)
        endif

        ; Update header
        ; If file does not exist, it will be created
        if n_elements(header) ne 0 then begin
            MDAP_WRITE_OUTPUT_UPDATE_HEADER, file, header, ndrp=ndrp, dx=dx, dy=dy, $
                                             w_range_sn=w_range_sn, nbin=nbin, bin_par=bin_par, $
                                             threshold_ston_bin=threshold_ston_bin, $
                                             w_range_analysis=w_range_analysis, $
                                             threshold_ston_analysis=threshold_ston_analysis, $
                                             analysis_par=analysis_par, $
                                             tpl_library_key=tpl_library_key, $
                                             ems_line_key=ems_line_key, abs_line_key=abs_line_key
        endif

;           uu += 1
;           print, uu
;           get_lun, unit
;           print, unit
;           if unit ne 100 then stop
;           free_lun, unit

        ; Update properties of DRP spectra
        MDAP_WRITE_OUTPUT_UPDATE_DRPS, file, ndrp=ndrp, xpos=xpos, ypos=ypos, $
                                       fraction_good=fraction_good, min_eq_max=min_eq_max, $
                                       signal=signal, noise=noise, bin_vreg=bin_vreg, $
                                       bin_indx=bin_indx, bin_weights=bin_weights, quiet=quiet

;           uu += 1
;           print, uu
;           get_lun, unit
;           print, unit
;           if unit ne 100 then stop
;           free_lun, unit

        ; Update properties of binned spectra
        MDAP_WRITE_OUTPUT_UPDATE_BINS, file, nbin=nbin, xbin_rlow=xbin_rlow, ybin_rupp=ybin_rupp, $
                                       rbin=rbin, bin_area=bin_area, bin_ston=bin_ston, $
                                       bin_n=bin_n, bin_flag=bin_flag, quiet=quiet

;           uu += 1
;           print, uu
;           get_lun, unit
;           print, unit
;           if unit ne 100 then stop
;           free_lun, unit

        ; Update the images for the binned spectra
        if n_elements(wave) ne 0 then begin
            MDAP_WRITE_OUTPUT_UPDATE_IMAGE, file, 'WAVE', wave, 'Wavelength (angstroms)'
        endif else $
            MDAP_WRITE_OUTPUT_BLANK_IMAGE_EXTENSION, file, 'WAVE', 'Wavelength (angstroms)'

;           uu += 1
;           print, uu
;           get_lun, unit
;           print, unit
;           if unit ne 100 then stop
;           free_lun, unit

        ; Update the spectral resolution vector
        if n_elements(sres) ne 0 then begin
            MDAP_WRITE_OUTPUT_UPDATE_IMAGE, file, 'SRES', sres, 'Spectral resolution (R, FWHM)'
        endif else $
            MDAP_WRITE_OUTPUT_BLANK_IMAGE_EXTENSION, file, 'SRES', 'Spectral resolution (R, FWHM)'

;           uu += 1
;           print, uu
;           get_lun, unit
;           print, unit
;           if unit ne 100 then stop
;           free_lun, unit

;       fits_info, file
;       stop

        ; Update the binned flux
        if n_elements(bin_flux) ne 0 then begin
            MDAP_WRITE_OUTPUT_UPDATE_IMAGE, file, 'FLUX', bin_flux, 'Binned flux (DRP units)'
        endif else $
            MDAP_WRITE_OUTPUT_BLANK_IMAGE_EXTENSION, file, 'FLUX', 'Binned flux (DRP units)'

;           uu += 1
;           print, uu
;           get_lun, unit
;           print, unit
;           if unit ne 100 then stop
;           free_lun, unit

;       fits_info, file
;       stop

        ; Update the inverse variances
        if n_elements(bin_ivar) ne 0 then begin
            MDAP_WRITE_OUTPUT_UPDATE_IMAGE, file, 'IVAR', bin_ivar, 'Inverse variance'
        endif else $
            MDAP_WRITE_OUTPUT_BLANK_IMAGE_EXTENSION, file, 'IVAR', 'Inverse variance'

;           uu += 1
;           print, uu
;           get_lun, unit
;           print, unit
;           if unit ne 100 then stop
;           free_lun, unit

;       fits_info, file
;       stop

        ; Update the bad pixel mask of the binned spectra
        if n_elements(bin_mask) ne 0 then begin
            MDAP_WRITE_OUTPUT_UPDATE_IMAGE, file, 'MASK', bin_mask, 'Bad pixel mask (0/1=good/bad)'
        endif else $
            MDAP_WRITE_OUTPUT_BLANK_IMAGE_EXTENSION, file, 'MASK', 'Bad pixel mask (0/1=good/bad)'

;           uu += 1
;           print, uu
;           get_lun, unit
;           print, unit
;           if unit ne 100 then stop
;           free_lun, unit

;       fits_info, file
;       stop

        ; Update the emission-line parameters
        MDAP_WRITE_OUTPUT_UPDATE_ELPAR, file, eml_par=eml_par, quiet=quiet

;           uu += 1
;           print, uu
;           get_lun, unit
;           print, unit
;           if unit ne 100 then stop
;           free_lun, unit

        ; Update the stellar kinematics
        MDAP_WRITE_OUTPUT_UPDATE_STELLAR_KIN, file, nbin=nbin, weights_ppxf=weights_ppxf, $
                                              add_poly_coeff_ppxf=add_poly_coeff_ppxf, $
                                              mult_poly_coeff_ppxf=mult_poly_coeff_ppxf, $
                                              stellar_kinematics_fit=stellar_kinematics_fit, $
                                              stellar_kinematics_err=stellar_kinematics_err, $
                                              chi2_ppxf=chi2_ppxf, quiet=quiet

;           uu += 1
;           print, uu
;           get_lun, unit
;           print, unit
;           if unit ne 100 then stop
;           free_lun, unit

        ; PPXF-only mask
        if n_elements(obj_fit_mask_ppxf) ne 0 then begin
            MDAP_WRITE_OUTPUT_UPDATE_IMAGE, file, 'SMSK', obj_fit_mask_ppxf, $
                                            'Pixel mask for PPXF fit (0/1=fit/omitted)'
        endif else begin
            MDAP_WRITE_OUTPUT_BLANK_IMAGE_EXTENSION, file, 'SMSK', $
                                               'Pixel mask for PPXF fit (0/1=fit/omitted)'
        endelse

;           uu += 1
;           print, uu
;           get_lun, unit
;           print, unit
;           if unit ne 100 then stop
;           free_lun, unit

        ; PPXF-only best fit
        if n_elements(bestfit_ppxf) ne 0 then begin
            MDAP_WRITE_OUTPUT_UPDATE_IMAGE, file, 'SMOD', bestfit_ppxf, $
                                            'Best-fit model for the stellar continuum'
        endif else begin
            MDAP_WRITE_OUTPUT_BLANK_IMAGE_EXTENSION, file, 'SMOD', $
                                               'Best-fit model for the stellar continuum'
        endelse

;           uu += 1
;           print, uu
;           get_lun, unit
;           print, unit
;           if unit ne 100 then stop
;           free_lun, unit

        ; Update the emission-line fitting results
        MDAP_WRITE_OUTPUT_UPDATE_STAR_AND_GAS_FIT, file, nbin=nbin, eml_par=eml_par, $
                                                    weights_gndf=weights_gndf, $
                                                    mult_poly_coeff_gndf=mult_poly_coeff_gndf, $
                                    emission_line_kinematics_avg=emission_line_kinematics_avg, $
                                    emission_line_kinematics_aer=emission_line_kinematics_aer, $
                                                    chi2_gndf=chi2_gndf, $
                                    emission_line_kinematics_ind=emission_line_kinematics_ind, $
                                    emission_line_kinematics_ier=emission_line_kinematics_ier, $
                                                    emission_line_sinst=emission_line_sinst, $
                                                    emission_line_omitted=emission_line_omitted, $
                                                    emission_line_intens=emission_line_intens, $
                                                    emission_line_interr=emission_line_interr, $
                                                    emission_line_fluxes=emission_line_fluxes, $
                                                    emission_line_flxerr=emission_line_flxerr, $
                                                    emission_line_EWidth=emission_line_EWidth, $
                                                    emission_line_EW_err=emission_line_EW_err, $
                                                    reddening_val=reddening_val, $
                                                    reddening_err=reddening_err, quiet=quiet

;           uu += 1
;           print, uu
;           get_lun, unit
;           print, unit
;           if unit ne 100 then stop
;           free_lun, unit

        ; GANDALF mask
        if n_elements(obj_fit_mask_gndf) ne 0 then begin
            MDAP_WRITE_OUTPUT_UPDATE_IMAGE, file, 'SGMSK', obj_fit_mask_gndf, $
                                            'Pixel mask for GANDALF fit (0/1=fit/omitted)'
        endif else begin
            MDAP_WRITE_OUTPUT_BLANK_IMAGE_EXTENSION, file, 'SGMSK', $
                                               'Pixel mask for GANDALF fit (0/1=fit/omitted)'
        endelse

;           uu += 1
;           print, uu
;           get_lun, unit
;           print, unit
;           if unit ne 100 then stop
;           free_lun, unit

        ; GANDALF best-fit model
        if n_elements(bestfit_gndf) ne 0 then begin
            MDAP_WRITE_OUTPUT_UPDATE_IMAGE, file, 'SGMOD', bestfit_gndf, $
                                            'Best-fit model of stellar+gas spectrum'
        endif else begin
            MDAP_WRITE_OUTPUT_BLANK_IMAGE_EXTENSION, file, 'SGMOD', $
                                                     'Best-fit model of stellar+gas spectrum'
        endelse

;           uu += 1
;           print, uu
;           get_lun, unit
;           print, unit
;           if unit ne 100 then stop
;           free_lun, unit

        ; GANDALF best-fit emission-only model
        if n_elements(eml_model) ne 0 then begin
            MDAP_WRITE_OUTPUT_UPDATE_IMAGE, file, 'ELMOD', eml_model, $
                                            'Best-fit emission-line-only model'
        endif else begin
            MDAP_WRITE_OUTPUT_BLANK_IMAGE_EXTENSION, file, 'ELMOD', $
                                                     'Best-fit emission-line-only model'
        endelse

;           uu += 1
;           print, uu
;           get_lun, unit
;           print, unit
;           if unit ne 100 then stop
;           free_lun, unit

        ; Update the emission-line parameters
        MDAP_WRITE_OUTPUT_UPDATE_ELOPAR, file, emlo_par=emlo_par, quiet=quiet

;           uu += 1
;           print, uu
;           get_lun, unit
;           print, unit
;           if unit ne 100 then stop
;           free_lun, unit

        ; Update the emission-line fitting results
        MDAP_WRITE_OUTPUT_UPDATE_EMISSION_LINE_FIT, file, nbin=nbin, emlo_par=emlo_par, $
                        elo_ew_kinematics_avg=elo_ew_kinematics_avg, $
                        elo_ew_kinematics_aer=elo_ew_kinematics_aer, $
                        elo_ew_kinematics_ind=elo_ew_kinematics_ind, $
                        elo_ew_kinematics_ier=elo_ew_kinematics_ier, elo_ew_sinst=elo_ew_sinst, $
                        elo_ew_omitted=elo_ew_omitted, elo_ew_intens=elo_ew_intens, $
                        elo_ew_interr=elo_ew_interr, elo_ew_fluxes=elo_ew_fluxes, $
                        elo_ew_flxerr=elo_ew_flxerr, elo_ew_EWidth=elo_ew_EWidth, $
                        elo_ew_EW_err=elo_ew_EW_err, elo_fb_kinematics_avg=elo_fb_kinematics_avg, $
                        elo_fb_kinematics_aer=elo_fb_kinematics_aer, $
                        elo_fb_kinematics_ind=elo_fb_kinematics_ind, $
                        elo_fb_kinematics_ier=elo_fb_kinematics_ier, elo_fb_sinst=elo_fb_sinst, $
                        elo_fb_omitted=elo_fb_omitted, elo_fb_intens=elo_fb_intens, $
                        elo_fb_interr=elo_fb_interr, elo_fb_fluxes=elo_fb_fluxes, $
                        elo_fb_flxerr=elo_fb_flxerr, elo_fb_EWidth=elo_fb_EWidth, $
                        elo_fb_EW_err=elo_fb_EW_err

;           uu += 1
;           print, uu
;           get_lun, unit
;           print, unit
;           if unit ne 100 then stop
;           free_lun, unit

        ; The best-fitting emission-line model from the
        ; emission-line-only analysis provided by Enci Wang
        if n_elements(elo_ew_eml_model) ne 0 then begin
            MDAP_WRITE_OUTPUT_UPDATE_IMAGE, file, 'ELOMEW', elo_ew_eml_model, $
                                                'Emission-line-only model (Enci Wang)'
        endif else begin
            MDAP_WRITE_OUTPUT_BLANK_IMAGE_EXTENSION, file, 'ELOMEW', $
                                                'Emission-line-only model (Enci Wang)'
        endelse

;           uu += 1
;           print, uu
;           get_lun, unit
;           print, unit
;           if unit ne 100 then stop
;           free_lun, unit

        ; The best-fitting emission-line model from the
        ; emission-line-only analysis provided by Francesco Belfiore
        if n_elements(elo_fb_eml_model) ne 0 then begin
            MDAP_WRITE_OUTPUT_UPDATE_IMAGE, file, 'ELOMFB', elo_fb_eml_model, $
                                                'Emission-line-only model (Francesco Belfiore)'
        endif else begin
            MDAP_WRITE_OUTPUT_BLANK_IMAGE_EXTENSION, file, 'ELOMFB', $
                                                'Emission-line-only model (Francesco Belfiore)'
        endelse

;           uu += 1
;           print, uu
;           get_lun, unit
;           print, unit
;           if unit ne 100 then stop
;           free_lun, unit

        ; Wavelength vector of the binned spectra that have been resolution
        ; matched to the spectral index system.
        if n_elements(si_bin_wave) ne 0 then begin
            MDAP_WRITE_OUTPUT_UPDATE_IMAGE, file, 'SIWAVE', si_bin_wave, $
                                                'Spectral-index resolution matched: wavelength'
        endif else begin
            MDAP_WRITE_OUTPUT_BLANK_IMAGE_EXTENSION, file, 'SIWAVE', $
                                                'Spectral-index resolution matched: wavelength'
        endelse

;           uu += 1
;           print, uu
;           get_lun, unit
;           print, unit
;           if unit ne 100 then stop
;           free_lun, unit

        ; Binned spectra that have been resolution matched to the spectral index
        ; system.
        if n_elements(si_bin_flux) ne 0 then begin
            MDAP_WRITE_OUTPUT_UPDATE_IMAGE, file, 'SIFLUX', si_bin_flux, $
                                                'Spectral-index resolution matched: binned flux'
        endif else begin
            MDAP_WRITE_OUTPUT_BLANK_IMAGE_EXTENSION, file, 'SIFLUX', $
                                                'Spectral-index resolution matched: binned flux'
        endelse

;           uu += 1
;           print, uu
;           get_lun, unit
;           print, unit
;           if unit ne 100 then stop
;           free_lun, unit

        ; Inverse variance of the binned spectra that have been resolution
        ; matched to the spectral index system.
        if n_elements(si_bin_ivar) ne 0 then begin
            MDAP_WRITE_OUTPUT_UPDATE_IMAGE, file, 'SIIVAR', si_bin_ivar, $
                                        'Spectral-index resolution matched: inverse variance'
        endif else begin
            MDAP_WRITE_OUTPUT_BLANK_IMAGE_EXTENSION, file, 'SIIVAR', $
                                        'Spectral-index resolution matched: inverse variance'
        endelse

;           uu += 1
;           print, uu
;           get_lun, unit
;           print, unit
;           if unit ne 100 then stop
;           free_lun, unit

        ; Bad pixel mask of the binned spectra that have been resolution matched
        ; to the spectral index system.
        if n_elements(si_bin_mask) ne 0 then begin
            MDAP_WRITE_OUTPUT_UPDATE_IMAGE, file, 'SIMASK', si_bin_mask, $
                                            'Spectral-index resolution matched: bad pixel mask'
        endif else begin
            MDAP_WRITE_OUTPUT_BLANK_IMAGE_EXTENSION, file, 'SIMASK', $
                                            'Spectral-index resolution matched: bad pixel mask'
        endelse

;           uu += 1
;           print, uu
;           get_lun, unit
;           print, unit
;           if unit ne 100 then stop
;           free_lun, unit

        ; Optimal template fit to each binned spectrum and resolution matched to
        ; the spectral index system.
        if n_elements(si_optimal_template) ne 0 then begin
            MDAP_WRITE_OUTPUT_UPDATE_IMAGE, file, 'SIOTPL', si_optimal_template, $
                                        'Spectral-index resolution matched: optimal template'
        endif else begin
            MDAP_WRITE_OUTPUT_BLANK_IMAGE_EXTENSION, file, 'SIOTPL', $
                                        'Spectral-index resolution matched: optimal template'
        endelse

;           uu += 1
;           print, uu
;           get_lun, unit
;           print, unit
;           if unit ne 100 then stop
;           free_lun, unit

        ; Broadened optimal template fit to each binned spectrum and resolution
        ; matched to the spectral index system.
        if n_elements(si_broad_optimal_template) ne 0 then begin
            MDAP_WRITE_OUTPUT_UPDATE_IMAGE, file, 'SIBOTPL', si_broad_optimal_template, $
                                'Spectral-index resolution matched: broadened optimal template'
        endif else begin
            MDAP_WRITE_OUTPUT_BLANK_IMAGE_EXTENSION, file, 'SIBOTPL', $
                                'Spectral-index resolution matched: broadened optimal template'
        endelse

;           uu += 1
;           print, uu
;           get_lun, unit
;           print, unit
;           if unit ne 100 then stop
;           free_lun, unit

        ; Update the spectral index parameters
        MDAP_WRITE_OUTPUT_UPDATE_SIPAR, file, abs_par=abs_par, quiet=quiet

;           uu += 1
;           print, uu
;           get_lun, unit
;           print, unit
;           if unit ne 100 then stop
;           free_lun, unit

        ; Update the spectral index measurements
        MDAP_WRITE_OUTPUT_UPDATE_SPECTRAL_INDICES, file, nbin=nbin, abs_par=abs_par, $
                                                   abs_line_indx_omitted=abs_line_indx_omitted, $
                                                   abs_line_indx_val=abs_line_indx_val, $
                                                   abs_line_indx_err=abs_line_indx_err, $
                                                   abs_line_indx_otpl=abs_line_indx_otpl, $
                                                   abs_line_indx_botpl=abs_line_indx_botpl

;           uu += 1
;           print, uu
            get_lun, unit
            print, unit
;           if unit ne 100 then stop
            free_lun, unit

;       fits_info, file
;       stop

END



