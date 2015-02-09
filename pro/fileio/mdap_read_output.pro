;+
; NAME:
;       MDAP_READ_OUTPUT
;
; PURPOSE:
;       Complement to MDAP_WRITE_OUTPUT, which reads DAP-produced output into
;       the provided data arrays.  See that procedure for the data model.
;
; CALLING SEQUENCE:
;       MDAP_READ_OUTPUT, file, header=header, dx=dx, dy=dy, w_range_sn=w_range_sn, xpos=xpos, $
;                         ypos=ypos, fraction_good=fraction_good, min_eq_max=min_eq_max, $
;                         signal=signal, noise=noise, bin_par=bin_par, $
;                         threshold_ston_bin=threshold_ston_bin, bin_vreg=bin_vreg, $
;                         bin_indx=bin_indx, bin_weights=bin_weights, wave=wave, sres=sres, $
;                         bin_flux=bin_flux, bin_ivar=bin_ivar, bin_mask=bin_mask, $
;                         xbin_rlow=xbin_rlow, ybin_rupp=ybin_rupp, rbin=rbin, bin_area=bin_area, $
;                         bin_ston=bin_ston, bin_n=bin_n, bin_flag=bin_flag, $
;                         w_range_analysis=w_range_analysis, $ 
;                         threshold_ston_analysis=threshold_ston_analysis, $
;                         tpl_library_key=tpl_library_key, ems_line_key=ems_line_key, $
;                         analysis_par=analysis_par, weights_ppxf=weights_ppxf, $
;                         add_poly_coeff_ppxf=add_pol_coeff_ppxf, $
;                         mult_poly_coeff_ppxf=mult_pol_coeff_ppxf, $
;                         stellar_kinematics_fit=stellar_kinematics_fit, $
;                         stellar_kinematics_err=stellar_kinematics_err, chi2_ppxf=chi2_ppxf, $
;                         obj_fit_mask_ppxf=obj_fit_mask_ppxf, bestfit_ppxf=bestfit_ppxf, $
;                         weights_gndf=weights_gndf, mult_poly_coeff_gndf=mult_poly_coeff_gndf, $
;                         emission_line_kinematics_avg=emission_line_kinematics_avg, $
;                         emission_line_kinematics_aer=emission_line_kinematics_aer, $
;                         chi2_gndf=chi2_gndf, $
;                         emission_line_kinematics_ind=emission_line_kinematics_ind, $
;                         emission_line_kinematics_ier=emission_line_kinematics_ier, $
;                         emission_line_sinst=emission_line_sinst, $
;                         emission_line_omitted=emission_line_omitted, $
;                         emission_line_intens=emission_line_intens, $
;                         emission_line_interr=emission_line_interr, $
;                         emission_line_fluxes=emission_line_fluxes, $
;                         emission_line_flxerr=emission_line_flxerr, $
;                         emission_line_EWidth=emission_line_EWidth, $
;                         emission_line_EW_err=emission_line_EW_err, $
;                         reddening_val=reddening_val, reddening_err=reddening_err, $
;                         obj_fit_mask_gndf=obj_fit_mask_gndf, bestfit_gndf=bestfit_gndf, $
;                         eml_model=eml_model, elo_ew_kinematics_fit=elo_ew_kinematics_fit, $
;                         elo_ew_kinematics_err=elo_ew_kinematics_err, $
;                         elo_ew_kinematics_ind=elo_ew_kinematics_ind, $
;                         elo_ew_kinematics_ier=elo_ew_kinematics_ier, $
;                         elo_ew_sinst=elo_ew_sinst, elo_ew_omitted=elo_ew_omitted, $
;                         elo_ew_intens=elo_ew_intens, elo_ew_interr=elo_ew_interr, $
;                         elo_ew_fluxes=elo_ew_fluxes, elo_ew_flxerr=elo_ew_flxerr, $
;                         elo_ew_EWidth=elo_ew_EWidth, elo_ew_EW_err=elo_ew_EW_err, $
;                         elo_ew_eml_model=elo_ew_eml_model, $
;                         elo_fb_kinematics_avg=elo_fb_kinematics_avg, $
;                         elo_fb_kinematics_aer=elo_fb_kinematics_aer, $
;                         elo_fb_kinematics_ind=elo_fb_kinematics_ind, $
;                         elo_fb_kinematics_ier=elo_fb_kinematics_ier, $
;                         elo_fb_sinst=elo_fb_sinst, elo_fb_omitted=elo_fb_omitted, $
;                         elo_fb_intens=elo_fb_intens, elo_fb_interr=elo_fb_interr, $
;                         elo_fb_fluxes=elo_fb_fluxes, elo_fb_flxerr=elo_fb_flxerr, $
;                         elo_fb_EWidth=elo_fb_EWidth, elo_fb_EW_err=elo_fb_EW_err, $
;                         elo_fb_eml_model=elo_fb_eml_model, abs_par=abs_par, $
;                         abs_line_key=abs_line_key, abs_line_indx_omitted=abs_line_indx_omitted, $
;                         abs_line_indx_val=abs_line_indx_val, $
;                         abs_line_indx_err=abs_line_indx_err, $
;                         abs_line_indx_otpl=abs_line_indx_otpl, $
;                         abs_line_indx_botpl=abs_line_indx_botpl, si_bin_wave=si_bin_wave, $
;                         si_bin_flux=si_bin_flux, si_bin_ivar=si_bin_ivar, $
;                         si_bin_mask=si_bin_mask, si_optimal_template=si_optimal_template, $
;                         si_broad_optimal_template= si_broad_optimal_template
;
; INPUTS:
;       file string
;               Name of the fits file of a DAP-produced fits file.
;
; OPTIONAL INPUTS:
;
;       IF DEFINED, the following inputs will be replaced with data read from
;       the DAP-produced fits file.
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
;               Structure containing the binning parameters.  TODO:
;               Currently cannot read any of the specifics related to
;               the definition of the bin edges.
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
;       - MAKE SURE THINGS ARE ABSTRACTED SUCH THAT I'M READING EXTENSIONS BY
;         THEIR NAME, NOT THEIR NUMBER!
;
;       - Include:?
;       eml_par EmissionLine[]
;               An array of EmissionLine structures used to define the emission
;               lines to be fit by GANDALF.
;
;
; BUGS:
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;       28 Oct 2014: Original Implementation by K. Westfall (KBW)
;       06 Dec 2014: (KBW) Accommodate changes for radial binning and
;                          velocity registration
;       12 Dec 2014: (KBW) New format incorporating emission-line only
;                          results
;       09 Jan 2015: (KBW) Include instrumental dispersion for GANDALF fit
;       09 Feb 2015: (KBW) Include fraction of good pixels and min(flux)
;                          == max(flux) flag in DRPS extension.
;-
;------------------------------------------------------------------------------

; The file that declares the functions used to define the names of the columns
; for the binary tables.
@mdap_set_output_file_cols

;-------------------------------------------------------------------------------
; Determine if the header needs to be read
FUNCTION MDAP_READ_NEED_HEADER, $
                header, dx, dy, w_range_sn, bin_par, threshold_ston_bin, w_range_analysis, $
                threshold_ston_analysis, analysis_par, tpl_library_key, ems_line_key, abs_line_key

        if n_elements(header)     ne 0 then              return, 1
        if n_elements(dx)         ne 0 then              return, 1
        if n_elements(dy)         ne 0 then              return, 1
        if n_elements(w_range_sn) ne 0 then              return, 1
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
                header, dx, dy, w_range_sn, bin_par, threshold_ston_bin, w_range_analysis, $
                threshold_ston_analysis, analysis_par, tpl_library_key, ems_line_key, abs_line_key

        dx = SXPAR(header, 'SPAXDX', /silent)
        dy = SXPAR(header, 'SPAXDY', /silent)
        w_range_sn = [ SXPAR(header, 'SNWAVE1', /silent),SXPAR(header, 'SNWAVE2', /silent) ]

        ; TODO: Move to a read_bin_par_header procedure?
        bin_par = MDAP_DEFINE_BIN_PAR()
        bin_par.type = SXPAR(header, 'BINTYPE', /silent)
        if bin_par.type ne 'NONE' then begin
            if SXPAR(header, 'BINVRG', /silent) eq 'True' then $
                bin_par.v_register = 1
            if SXPAR(header, 'BINWGT', /silent) eq 'Optimal' then $
                bin_par.optimal_weighting = 1
            if bin_par.type eq 'STON' then $
                bin_par.ston = SXPAR(header, 'BINSN', /silent)
            if bin_par.type eq 'RADIAL' then begin
                bin_par.cx = SXPAR(header, 'BINCX', /silent)
                bin_par.cy = SXPAR(header, 'BINCY', /silent)
                bin_par.pa = SXPAR(header, 'BINPA', /silent)
                bin_par.ell = SXPAR(header, 'BINELL', /silent)
                bin_par.rscale = SXPAR(header, 'BINSCL', /silent)
            endif
        endif

        threshold_ston_bin = SXPAR(header, 'MINSNBIN', /silent)
        w_range_analysis = [ SXPAR(header, 'FITWAVE1', /silent),SXPAR(header, 'FITWAVE2', /silent) ]
        threshold_ston_analysis = SXPAR(header, 'MINSNFIT', /silent)

        ; TODO: Move to a read_analysis_par_header procedure?
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
                xpos, ypos, fraction_good, min_eq_max, signal, noise, bin_vreg, bin_indx, $
                bin_weights

        if n_elements(xpos) ne 0 then           return, 1
        if n_elements(ypos) ne 0 then           return, 1
        if n_elements(fraction_good)ne 0 then   return, 1
        if n_elements(min_eq_max) ne 0 then     return, 1
        if n_elements(signal) ne 0 then         return, 1
        if n_elements(noise) ne 0 then          return, 1
        if n_elements(bin_vreg) ne 0 then       return, 1
        if n_elements(bin_indx) ne 0 then       return, 1
        if n_elements(bin_weights) ne 0 then    return, 1

        return, 0
END

; Read the full DRPS binary table (more efficient than reading one column at a
; time?)
PRO MDAP_READ_DRPS, $
                file, xpos, ypos, fraction_good, min_eq_max, signal, noise, bin_vreg, bin_indx, $
                bin_weights

        cols = MDAP_SET_DRPS_COLS()

        ; !! ORDER IS IMPORTANT !!
        fxbopen, unit, file, 'DRPS'
        fxbreadm, unit, cols, xpos, ypos, fraction_good, min_eq_max, signal, noise, bin_vreg, $
                  bin_indx, bin_weights
        fxbclose, unit
        free_lun, unit
END

;-------------------------------------------------------------------------------
; Check if data from the BINS extension is desired
FUNCTION MDAP_READ_WANT_BINS, $
                xbin_rlow, ybin_rupp, rbin, bin_area, bin_ston, bin_n, bin_flag

        if n_elements(xbin_rlow) ne 0 then return, 1
        if n_elements(ybin_rupp) ne 0 then return, 1
        if n_elements(rbin)      ne 0 then return, 1
        if n_elements(bin_area)  ne 0 then return, 1
        if n_elements(bin_ston)  ne 0 then return, 1
        if n_elements(bin_n)     ne 0 then return, 1
        if n_elements(bin_flag)  ne 0 then return, 1

        return, 0
END

; Read the full BINS binary table (more efficient than reading one column at a
; time?)
PRO MDAP_READ_BINS, $
                file, xbin_rlow, ybin_rupp, rbin, bin_area, bin_ston, bin_n, bin_flag

        cols = MDAP_SET_BINS_COLS()

        fxbopen, unit, file, 'BINS'
        fxbreadm, unit, cols, xbin_rlow, ybin_rupp, rbin, bin_area, bin_ston, bin_n, bin_flag
        fxbclose, unit
        free_lun, unit
END

; Read an image extension
FUNCTION MDAP_READ_IMG, $
                file, exten

        errmsg = ''                                     ; Suppress errormessages
        unit=fxposit(file, exten, errmsg=errmsg)        ; Find extension and open
        if unit eq -1 then $                            ; Check there was no error
            message, errmsg
        
        data=READFITS(unit)                             ; Read the data
        free_lun, unit                                  ; Close the file and free the LUN
        return, data                                    ; Return the data
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

        if (size(weights_ppxf))[0] eq 2 then $
            weights_ppxf = transpose(temporary(weights_ppxf))
        if (size(add_poly_coeff_ppxf))[0] eq 2 then $
            add_poly_coeff_ppxf = transpose(temporary(add_poly_coeff_ppxf))
        if (size(mult_poly_coeff_ppxf))[0] eq 2 then $
            mult_poly_coeff_ppxf = transpose(temporary(mult_poly_coeff_ppxf))
        if (size(stellar_kinematics_fit))[0] eq 2 then $
            stellar_kinematics_fit = transpose(temporary(stellar_kinematics_fit))
        if (size(stellar_kinematics_err))[0] eq 2 then $
            stellar_kinematics_err = transpose(temporary(stellar_kinematics_err))
END

;-------------------------------------------------------------------------------
; Check if data from the SGFIT extension is desired
FUNCTION MDAP_READ_WANT_SGFIT, $ 
                weights_gndf, mult_poly_coeff_gndf, emission_line_kinematics_avg, $
                emission_line_kinematics_aer, chi2_gndf, emission_line_kinematics_ind, $
                emission_line_kinematics_ier, emission_line_sinst, emission_line_omitted, $
                emission_line_intens, emission_line_interr, emission_line_fluxes, $
                emission_line_flxerr, emission_line_EWidth, emission_line_EW_err, reddening_val, $
                reddening_err

        if n_elements(weights_gndf)                 ne 0 then return, 1
        if n_elements(mult_poly_coeff_gndf)         ne 0 then return, 1
        if n_elements(emission_line_kinematics_avg) ne 0 then return, 1
        if n_elements(emission_line_kinematics_aer) ne 0 then return, 1
        if n_elements(chi2_gndf)                    ne 0 then return, 1
        if n_elements(emission_line_kinematics_ind) ne 0 then return, 1
        if n_elements(emission_line_kinematics_ier) ne 0 then return, 1
        if n_elements(emission_line_sinst)          ne 0 then return, 1
        if n_elements(emission_line_omitted)        ne 0 then return, 1
        if n_elements(emission_line_intens)         ne 0 then return, 1
        if n_elements(emission_line_interr)         ne 0 then return, 1
        if n_elements(emission_line_fluxes)         ne 0 then return, 1
        if n_elements(emission_line_flxerr)         ne 0 then return, 1
        if n_elements(emission_line_EWidth)         ne 0 then return, 1
        if n_elements(emission_line_EW_err)         ne 0 then return, 1
        if n_elements(reddening_val)                ne 0 then return, 1
        if n_elements(reddening_err)                ne 0 then return, 1

        return, 0
END

;-------------------------------------------------------------------------------
; Read the full SGFIT binary table (more efficient than reading one column at a
; time?)
PRO MDAP_READ_SGFIT, $
                file, weights_gndf, mult_poly_coeff_gndf, emission_line_kinematics_avg, $
                emission_line_kinematics_aer, chi2_gndf, emission_line_kinematics_ind, $
                emission_line_kinematics_ier, emission_line_sinst, emission_line_omitted, $
                emission_line_intens, emission_line_interr, emission_line_fluxes, $
                emission_line_flxerr, emission_line_EWidth, emission_line_EW_err, reddening_val, $
                reddening_err

        cols = MDAP_SET_SGFIT_COLS()

        ; !! ORDER IS IMPORTANT !!
        fxbopen, unit, file, 'SGFIT'
        fxbreadm, unit, cols, weights_gndf, mult_poly_coeff_gndf, emission_line_kinematics_avg, $
                              emission_line_kinematics_aer, chi2_gndf, reddening_val, $
                              reddening_err, emission_line_omitted, emission_line_intens, $
                              emission_line_interr, emission_line_kinematics_ind, $
                              emission_line_kinematics_ier, emission_line_sinst, $
                              emission_line_fluxes, emission_line_flxerr, emission_line_EWidth, $
                              emission_line_EW_err

        fxbclose, unit
        free_lun, unit

        ; TODO: 2D data needs to be transposed!

        if (size(weights_gndf))[0] eq 2 then $
            weights_gndf = transpose(temporary(weights_gndf))
        if (size(mult_poly_coeff_gndf))[0] eq 2 then $
            mult_poly_coeff_gndf = transpose(temporary(mult_poly_coeff_gndf))
        if (size(emission_line_kinematics_avg))[0] eq 2 then $
            emission_line_kinematics_avg = transpose(temporary(emission_line_kinematics_avg))
        if (size(emission_line_kinematics_aer))[0] eq 2 then $
            emission_line_kinematics_aer = transpose(temporary(emission_line_kinematics_aer))
        if (size(reddening_val))[0] eq 2 then $
            reddening_val = transpose(temporary(reddening_val))
        if (size(reddening_err))[0] eq 2 then $
            reddening_err = transpose(temporary(reddening_err))
        if (size(emission_line_intens))[0] eq 2 then $
            emission_line_intens = transpose(temporary(emission_line_intens))
        if (size(emission_line_interr))[0] eq 2 then $
            emission_line_interr = transpose(temporary(emission_line_interr))
        if (size(emission_line_kinematics_ind))[0] eq 3 then $
            emission_line_kinematics_ind=transpose(temporary(emission_line_kinematics_ind), [2,0,1])
        if (size(emission_line_kinematics_ier))[0] eq 3 then $
            emission_line_kinematics_ier=transpose(temporary(emission_line_kinematics_ier), [2,0,1])
        if (size(emission_line_sinst))[0] eq 2 then $
            emission_line_sinst = transpose(temporary(emission_line_sinst))
        if (size(emission_line_fluxes))[0] eq 2 then $
            emission_line_fluxes = transpose(temporary(emission_line_fluxes))
        if (size(emission_line_flxerr))[0] eq 2 then $
            emission_line_flxerr = transpose(temporary(emission_line_flxerr))
        if (size(emission_line_EWidth))[0] eq 2 then $
            emission_line_EWidth = transpose(temporary(emission_line_EWidth))
        if (size(emission_line_EW_err))[0] eq 2 then $
            emission_line_EW_err = transpose(temporary(emission_line_EW_err))

END

;-------------------------------------------------------------------------------
; Check if data from the ELOFIT extension is desired
FUNCTION MDAP_READ_WANT_ELOFIT, $
                elo_ew_kinematics_avg, elo_ew_kinematics_aer, elo_ew_kinematics_ind, $
                elo_ew_kinematics_ier, elo_ew_sinst, elo_ew_omitted, elo_ew_intens, $
                elo_ew_interr, elo_ew_fluxes, elo_ew_flxerr, elo_ew_EWidth, elo_ew_EW_err, $
                elo_fb_kinematics_avg, elo_fb_kinematics_aer, elo_fb_kinematics_ind, $
                elo_fb_kinematics_ier, elo_fb_sinst, elo_fb_omitted, elo_fb_intens, elo_fb_interr, $
                elo_fb_fluxes, elo_fb_flxerr, elo_fb_EWidth, elo_fb_EW_err

        if n_elements(elo_ew_kinematics_avg) eq 0 then return, 1
        if n_elements(elo_ew_kinematics_aer) eq 0 then return, 1
        if n_elements(elo_ew_kinematics_ind) eq 0 then return, 1
        if n_elements(elo_ew_kinematics_ier) eq 0 then return, 1
        if n_elements(elo_ew_sinst) eq 0 then return, 1
        if n_elements(elo_ew_omitted) eq 0 then return, 1
        if n_elements(elo_ew_intens) eq 0 then return, 1
        if n_elements(elo_ew_interr) eq 0 then return, 1
        if n_elements(elo_ew_fluxes) eq 0 then return, 1
        if n_elements(elo_ew_flxerr) eq 0 then return, 1
        if n_elements(elo_ew_EWidth) eq 0 then return, 1
        if n_elements(elo_ew_EW_err) eq 0 then return, 1
        if n_elements(elo_fb_kinematics_avg) eq 0 then return, 1
        if n_elements(elo_fb_kinematics_aer) eq 0 then return, 1
        if n_elements(elo_fb_kinematics_ind) eq 0 then return, 1
        if n_elements(elo_fb_kinematics_ier) eq 0 then return, 1
        if n_elements(elo_fb_sinst) eq 0 then return, 1
        if n_elements(elo_fb_omitted) eq 0 then return, 1
        if n_elements(elo_fb_intens) eq 0 then return, 1
        if n_elements(elo_fb_interr) eq 0 then return, 1
        if n_elements(elo_fb_fluxes) eq 0 then return, 1
        if n_elements(elo_fb_flxerr) eq 0 then return, 1
        if n_elements(elo_fb_EWidth) eq 0 then return, 1
        if n_elements(elo_fb_EW_err) eq 0 then return, 1

        return, 0
END

;-------------------------------------------------------------------------------
; Read the full ELOFIT binary table (more efficient than reading one column at a
; time?)
PRO MDAP_READ_ELOFIT, $
                file, elo_ew_kinematics_avg, elo_ew_kinematics_aer, elo_ew_kinematics_ind, $
                elo_ew_kinematics_ier, elo_ew_sinst, elo_ew_omitted, elo_ew_intens, $
                elo_ew_interr, elo_ew_fluxes, elo_ew_flxerr, elo_ew_EWidth, elo_ew_EW_err, $
                elo_fb_kinematics_avg, elo_fb_kinematics_aer, elo_fb_kinematics_ind, $
                elo_fb_kinematics_ier, elo_fb_sinst, elo_fb_omitted, elo_fb_intens, elo_fb_interr, $
                elo_fb_fluxes, elo_fb_flxerr, elo_fb_EWidth, elo_fb_EW_err

        cols = MDAP_SET_ELOFIT_COLS()

        fxbopen, unit, file, 'ELOFIT'
        fxbreadm, unit, cols, elo_ew_kinematics_avg, elo_ew_kinematics_aer, elo_ew_omitted, $
                              elo_ew_intens, elo_ew_interr, elo_ew_kinematics_ind, $
                              elo_ew_kinematics_ier, elo_ew_sinst, elo_ew_fluxes, elo_ew_flxerr, $
                              elo_ew_EWidth, elo_ew_EW_err, elo_fb_kinematics_avg, $
                              elo_fb_kinematics_aer, elo_fb_omitted, elo_fb_intens, elo_fb_interr, $
                              elo_fb_kinematics_ind, elo_fb_kinematics_ier, elo_fb_sinst, $
                              elo_fb_fluxes, elo_fb_flxerr, elo_fb_EWidth, elo_fb_EW_err
        fxbclose, unit
        free_lun, unit

        ; TODO: 2D data needs to be transposed!

        if (size(elo_ew_kinematics_avg))[0] eq 2 then $
            elo_ew_kinematics_avg = transpose(temporary(elo_ew_kinematics_avg))
        if (size(elo_ew_kinematics_aer))[0] eq 2 then $
            elo_ew_kinematics_aer = transpose(temporary(elo_ew_kinematics_aer))
        if (size(elo_ew_intens))[0] eq 2 then $
            elo_ew_intens = transpose(temporary(elo_ew_intens))
        if (size(elo_ew_interr))[0] eq 2 then $
            elo_ew_interr = transpose(temporary(elo_ew_interr))
        if (size(elo_ew_kinematics_ind))[0] eq 3 then $
            elo_ew_kinematics_ind = transpose(temporary(elo_ew_kinematics_ind), [2,0,1])
        if (size(elo_ew_kinematics_ier))[0] eq 3 then $
            elo_ew_kinematics_ier = transpose(temporary(elo_ew_kinematics_ier), [2,0,1])
        if (size(elo_ew_sinst))[0] eq 2 then $
            elo_ew_sinst = transpose(temporary(elo_ew_sinst))
        if (size(elo_ew_fluxes))[0] eq 2 then $
            elo_ew_fluxes = transpose(temporary(elo_ew_fluxes))
        if (size(elo_ew_flxerr))[0] eq 2 then $
            elo_ew_flxerr = transpose(temporary(elo_ew_flxerr))
        if (size(elo_ew_EWidth))[0] eq 2 then $
            elo_ew_EWidth = transpose(temporary(elo_ew_EWidth))
        if (size(elo_ew_EW_err))[0] eq 2 then $
            elo_ew_EW_err = transpose(temporary(elo_ew_EW_err))
        if (size(elo_fb_kinematics_avg))[0] eq 2 then $
            elo_fb_kinematics_avg = transpose(temporary(elo_fb_kinematics_avg))
        if (size(elo_fb_kinematics_aer))[0] eq 2 then $
            elo_fb_kinematics_aer = transpose(temporary(elo_fb_kinematics_aer))
        if (size(elo_fb_intens))[0] eq 2 then $
            elo_fb_intens = transpose(temporary(elo_fb_intens))
        if (size(elo_fb_interr))[0] eq 2 then $
            elo_fb_interr = transpose(temporary(elo_fb_interr))
        if (size(elo_fb_kinematics_ind))[0] eq 3 then $
            elo_fb_kinematics_ind = transpose(temporary(elo_fb_kinematics_ind), [2,0,1])
        if (size(elo_fb_kinematics_ier))[0] eq 3 then $
            elo_fb_kinematics_ier = transpose(temporary(elo_fb_kinematics_ier), [2,0,1])
        if (size(elo_fb_sinst))[0] eq 2 then $
            elo_fb_sinst = transpose(temporary(elo_fb_sinst))
        if (size(elo_fb_fluxes))[0] eq 2 then $
            elo_fb_fluxes = transpose(temporary(elo_fb_fluxes))
        if (size(elo_fb_flxerr))[0] eq 2 then $
            elo_fb_flxerr = transpose(temporary(elo_fb_flxerr))
        if (size(elo_fb_EWidth))[0] eq 2 then $
            elo_fb_EWidth = transpose(temporary(elo_fb_EWidth))
        if (size(elo_fb_EW_err))[0] eq 2 then $
            elo_fb_EW_err = transpose(temporary(elo_fb_EW_err))

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

        if (size(abs_line_indx_omitted))[0] eq 2 then $
            abs_line_indx_omitted = transpose(temporary(abs_line_indx_omitted))
        if (size(abs_line_indx_val))[0] eq 2 then $
            abs_line_indx_val = transpose(temporary(abs_line_indx_val))
        if (size(abs_line_indx_err))[0] eq 2 then $
            abs_line_indx_err = transpose(temporary(abs_line_indx_err))
        if (size(abs_line_indx_otpl))[0] eq 2 then $
            abs_line_indx_otpl = transpose(temporary(abs_line_indx_otpl))
        if (size(abs_line_indx_botpl))[0] eq 2 then $
            abs_line_indx_botpl = transpose(temporary(abs_line_indx_botpl))

END



;-------------------------------------------------------------------------------
; Order should not really be important here
; TODO: Add information such that eml_par can be read?
PRO MDAP_READ_OUTPUT, $
                file, header=header, dx=dx, dy=dy, w_range_sn=w_range_sn, xpos=xpos, ypos=ypos, $
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
                si_broad_optimal_template=si_broad_optimal_template

        if file_test(file) eq 0 then $                          ; Make sure file exists
            message, 'File does not exist!'

        ; Check that the file looks like a DAP file
        if MDAP_CHECK_OUTPUT_FILE(file) eq 0 then $
            message, 'File does not have DAP-like (hard-coded) properties!'

        ; If the header is needed/requested, read it
        if MDAP_READ_NEED_HEADER(header, dx, dy, w_range_sn, bin_par, threshold_ston_bin, $
                                 w_range_analysis, threshold_ston_analysis, analysis_par, $
                                 tpl_library_key, ems_line_key, abs_line_key) ne 0 then begin

            header=headfits(file, exten=0)
            MDAP_READ_SET_HEADER_DATA, header, dx, dy, w_range_sn, bin_par, threshold_ston_bin, $
                                       w_range_analysis, threshold_ston_analysis, analysis_par, $
                                       tpl_library_key, ems_line_key, abs_line_key
        endif

        ; Read the DRPS extension if any of its vectors are requested
        if MDAP_READ_WANT_DRPS(xpos, ypos, fraction_good, min_eq_max, signal, noise, bin_vreg, $
                               bin_indx, bin_weights) eq 1 then begin
            MDAP_READ_DRPS, file, xpos, ypos, fraction_good, min_eq_max, signal, noise, bin_vreg, $
                            bin_indx, bin_weights
        endif

        ; Read the BINS extension if any of its vectors are requested
        if MDAP_READ_WANT_BINS(xbin_rlow, ybin_rupp, rbin, bin_area, bin_ston, bin_n, $
                               bin_flag) eq 1 then begin
            MDAP_READ_BINS, file, xbin_rlow, ybin_rupp, bin_area, bin_ston, bin_n, bin_flag
        endif

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
;       if n_elements(eml_par) ne 0 then $
;           eml_par = MDAP_READ_EMISSION_LINE_PARAMETERS(file)

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
                                emission_line_sinst, emission_line_omitted, emission_line_intens, $
                                emission_line_interr, emission_line_fluxes, emission_line_flxerr, $
                                emission_line_EWidth, emission_line_EW_err, reddening_val, $
                                reddening_err) eq 1 then begin

            MDAP_READ_SGFIT, file, weights_gndf, mult_poly_coeff_gndf, $
                             emission_line_kinematics_avg, emission_line_kinematics_aer, $
                             chi2_gndf, emission_line_kinematics_ind, $
                             emission_line_kinematics_ier, emission_line_sinst, $
                             emission_line_omitted, emission_line_intens, emission_line_interr, $
                             emission_line_fluxes, emission_line_flxerr, emission_line_EWidth, $
                             emission_line_EW_err, reddening_val, reddening_err

        endif

        ; Read the fit_pixel mask resulting from GANDALF, if requested
        if n_elements(obj_fit_mask_gndf) ne 0 then $
            obj_fit_mask_gndf = MDAP_READ_IMG(file, 'SGMSK')

        ; Read the best-fitting GANDALF model, if requested
        if n_elements(bestfit_gndf) ne 0 then $
            bestfit_gndf = MDAP_READ_IMG(file, 'SGMOD')

        ; Read the best-fitting emission-line model from GANDALF, if requested
        if n_elements(eml_model) ne 0 then $
            eml_model = MDAP_READ_IMG(file, 'ELMOD')

        ; Read the ELOFIT extension if any of its vectors are requested
        if MDAP_READ_WANT_ELOFIT(elo_ew_kinematics_avg, elo_ew_kinematics_aer, $
                                 elo_ew_kinematics_ind, elo_ew_kinematics_ier, elo_ew_sinst, $
                                 elo_ew_omitted, elo_ew_intens, elo_ew_interr, elo_ew_fluxes, $
                                 elo_ew_flxerr, elo_ew_EWidth, elo_ew_EW_err, $
                                 elo_fb_kinematics_avg, elo_fb_kinematics_aer, $
                                 elo_fb_kinematics_ind, elo_fb_kinematics_ier, elo_fb_sinst, $
                                 elo_fb_omitted, elo_fb_intens, elo_fb_interr, elo_fb_fluxes, $
                                 elo_fb_flxerr, elo_fb_EWidth, elo_fb_EW_err) eq 1 then begin

            MDAP_READ_ELOFIT, file, elo_ew_kinematics_avg, elo_ew_kinematics_aer, $
                              elo_ew_kinematics_ind, elo_ew_kinematics_ier, elo_ew_sinst, $
                              elo_ew_omitted, elo_ew_intens, elo_ew_interr, elo_ew_fluxes, $
                              elo_ew_flxerr, elo_ew_EWidth, elo_ew_EW_err, elo_fb_kinematics_avg, $
                              elo_fb_kinematics_aer, elo_fb_kinematics_ind, elo_fb_kinematics_ier, $
                              elo_fb_sinst, elo_fb_omitted, elo_fb_intens, elo_fb_interr, $
                              elo_fb_fluxes, elo_fb_flxerr, elo_fb_EWidth, elo_fb_EW_err
        endif

        ; Read the best-fitting emission-line-only model from Enci Wang's code, if requested
        if n_elements(elo_ew_eml_model) ne 0 then $
            elo_ew_eml_model = MDAP_READ_IMG(file, 'ELOMEW')

        ; Read the best-fitting emission-line-only model from Francesco
        ; Belfiore's code, if requested
        if n_elements(elo_fb_eml_model) ne 0 then $
            elo_fb_eml_model = MDAP_READ_IMG(file, 'ELOMFB')

        ; Read the wavelength vector for the spectral-index-system
        ; resolution-matched spectra, if requested
        if n_elements(si_bin_wave) ne 0 then $
            si_bin_wave=MDAP_READ_IMG(file, 'SIWAVE')

        ; Read the flux for the spectral-index-system resolution-matched
        ; spectra, if requested
        if n_elements(si_bin_flux) ne 0 then $
            si_bin_flux=MDAP_READ_IMG(file, 'SIFLUX')

        ; Read the inverse varance for the spectral-index-system
        ; resolution-matched spectra, if requested
        if n_elements(si_bin_ivar) ne 0 then $
            si_bin_ivar=MDAP_READ_IMG(file, 'SIIVAR')

        ; Read the bad pixel mask for the spectral-index-system
        ; resolution-matched spectra, if requested
        if n_elements(si_bin_mask) ne 0 then $
            si_bin_mask=MDAP_READ_IMG(file, 'SIMASK')

        ; Read the optimal template matched to the resolution of the
        ; spectral-index system, if requested
        if n_elements(si_optimal_template) ne 0 then $
            si_optimal_template=MDAP_READ_IMG(file, 'SIOTPL')

        ; Read the broadened optimal template matched to the resolution of the
        ; spectral-index system, if requested
        if n_elements(si_broad_optimal_template) ne 0 then $
            si_broad_optimal_template=MDAP_READ_IMG(file, 'SIBOTPL')

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



