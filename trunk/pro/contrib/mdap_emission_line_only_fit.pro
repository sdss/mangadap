;+
; NAME:
;       MDAP_EMISSION_LINE_ONLY_FIT
;
; PURPOSE:
;       Fit emission lines using both Enci Wang's code and Francesco
;       Belfiore's codes that were provided in Dec 2014.
;
;       The returned velocity dispersion *are* corrected for the
;       instrumental resolution provided in the sres vector.  The
;       corrections are applied for all fits, regardless of whether or
;       not the fit is flagged for omission.  The corrected dispersions
;       are allowed to be "negative"; in this case, the instrumental
;       dispersion is *larger* than the measured dispersion.
;
;       To recover the measured dispersion for each line, compute:
;
;           sigma^2_obs = sigma_kin*abs(sigma_kin) + sigma^2_inst
;
;       where sigma_kin is the value returned in
;       emission_line_kinematics_individual[i,j,1] and sigma_inst is
;       returned in emission_line_sinst[i,j]
;
; CALLING SEQUENCE:
;       MDAP_EMISSION_LINE_ONLY_FIT, wave, sres, flux, ivar, mask, best_fit_continuum, $
;                                    best_fit_mask, stellar_kinematics, eml_model, $
;                                    emission_line_kinematics, emission_line_kinematics_perr, $
;                                    emission_line_kinematics_serr, emission_line_kineamtics_n, $
;                                    emission_line_omitted, emission_line_window, $
;                                    emission_line_kinematics_individual, $
;                                    emission_line_kinematics_individual_err, emission_line_sinst, $
;                                    emission_line_intens, emission_line_intens_err, $
;                                    emission_line_fluxes, emission_line_fluxes_err, $
;                                    emission_line_EW, emission_line_EW_err, $
;                                    eml_par=eml_par, version=version, /enci, /zero_base, /quiet, $
;                                    /dbg 
;
; INPUTS:
;       wave dblarr[C]
;               Wavelength of all C spectral channels, which is the same for all
;               N spectra.
;
;       sres dblarr[C]
;               Spectral resolution (R=lambda/delta lambda) of the
;               spectrum at the wavelengths provided by the wavelength
;               vector.
;
;       flux dblarr[N][C]
;               Object flux for each of N spectra with C spectral channels.
;
;       ivar dblarr[N][C]
;               Inverse variance in object flux for each of N spectra with C
;               spectral channels.
;
;       mask dblarr[N][C]
;               Pixel mask for each of N spectra with C spectral channels.
;
;       best_fit_continuum dblarr[N][C]
;               Best fitting continuum spectrum (obtained by PPXF) for
;               each of the N spectra.
;
;       best_fit_mask dblarr[N][C]
;               Bad pixel mask for pixels in the best-fitting continuum
;               data.  Pixels that are not/are masked have values of
;               0.0/1.0.
;
;       stellar_kinematics dblarr[B][K]
;               The best-fit stellar kinematics for each of the N fitted
;               (or fixed) input galaxy spectra with K fitted moments.
;               stellar_kinematics[*,0] must be  fitted velocity (cz).
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;       /enci
;           Use Enci Wang's fitting code
;
;       /belfiore
;           Use Francesco Belfiore's fitting code
;
;       /zero_base
;           Impose a zero baseline during the fit
;
;       /quiet
;           Suppress output to the terminal
;
;       /dbg
;           Run in "debug mode" by only fitting the first provided
;           spectrum.
;
; OUTPUT:
;       eml_model dblarr[N][C]
;               Best-fitting emission-line-only model for each of the N
;               spectra.
;
;       emission_line_kinematics dblarr[N][8]
;       emission_line_kinematics_perr dblarr[N][8]
;       emission_line_kinematics_serr dblarr[N][8]
;       emission_line_kinematics_n intarr[N]
;               Weighted means of the best-fit V and sigma for all
;               "well-fit" emission lines in the N galaxy spectra.
;               There are eight columns because there are four weighting
;               methods:
;                   - [0:1] - flux weighted
;                   - [2:3] - velocity-error weighted
;                   - [4:5] - both flux and velocity-error weighted
;                   - [6:7] - uniformly weighted
;               The errors are calculated via nominal propagation
;               (_perr) and by using the weighted standard deviation to
;               calculate the standard error (_serr); the number of
;               measurements used in the mean is given in _n.
;
;       emission_line_omitted intarr[N][E]
;               Flag setting whether or not an emission-line was fit for
;               all E lines fit by Enci's code (0-fit;1-omitted).
;
;       emission_line_window intarr[N][E][2]
;               Wavelength limits used in the fit of each (set of)
;               emission line(s).
;
;       emission_line_baseline intarr[N][E]
;       emission_line_base_err intarr[N][E]
;               Constant baseline fit beneath each (group of) emission
;               line(s) and its error.
;
;       emission_line_kinematics_individual dblarr[N][E][2]
;       emission_line_kinematics_individual_err dblarr[N][E][2]
;               Kinematics (V and sigma) and errors for each of the E
;               fitted emission lines.
;
;       emission_line_sinst dblarr[N][E]
;               The instrumental dispersion at the locations of the
;               fitted emission lines.
;
;       emission_line_intens dblarr[N][E]
;       emission_line_intens_err dblarr[N][E]
;               Best-fitting emission-line intensities and errors of the
;               E fitted emission lines for each of the N galaxy
;               spectra.
;
;       emission_line_fluxes  dblarr[N][E]
;       emission_line_fluxes_err dblarr[N][E]
;               Gaussian line fluxes and errors for the E fitted
;               emission lines in the N input galaxy spectra. 
;
;       emission_line_EW dblarr[N][E]
;       emission_line_EW_err dblarr[N][E]
;               Equivalent widths and errors for the E fitted emission
;               lines for each of the N input galaxy spectra.  These
;               equivalent widths are computed by taking the ratio of
;               the emission_line_fluxes and the median value of the
;               emission-line-free spectrum within 5 and 10 sigma of the
;               emission line, where 'sigma' is the fitted velocity
;               dispersion
;
; OPTIONAL OUTPUT:
;       eml_par EmissionLine[E]
;               Emission-line parameter structure, of the same form used
;               by GANDALF, but with the prameters of the E lines fitted
;               by Enci's code.
;
;       version string
;               Module version.  If requested, the module is not executed and
;               only version flag is returned.
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; TODO:
;       - Allow the emission lines to be fit with a non-zero continuum.
;       - THE DOUBLET LINES FROM FRANCESCO'S CODE SHOULD NOT BE TREATED
;         AS IF THEY PROVIDE INDEPENDENT MEASUREMENTS OF THE VELOCITY!!
;
; PROCEDURES CALLED:
;       MDAP_DEFINE_EMISSION_LINES_ENCI_BELFIORE()
;       MDAP_FIT_EMISSION_LINE_SPECTRUM_BELFIORE
;       MDAP_FIT_EMISSION_LINE_SPECTRUM_ENCI
;       MDAP_MEASURE_EMISSION_LINE_EQUIV_WIDTH
;       MDAP_INSTRUMENTAL_DISPERSION()
;       MDAP_AUTO_EXCLUDE_EML_FROM_KIN()
;       MDAP_LINES_FOR_MEAN_GAS_KINEMATICS()
;       MDAP_MEAN_GAS_KINEMATICS
;
; INTERNAL SUPPORT ROUTINES:
;       MDAP_SAVE_EMISSION_LINE_FIT_CHECK_OMIT()
;       MDAP_SAVE_EMISSION_LINE_FIT_ENCI_SINGLE
;       MDAP_SAVE_EMISSION_LINE_FIT_ENCI
;       MDAP_SAVE_EMISSION_LINE_FIT_BELFIORE
;
; REVISION HISTORY:
;       10 Dec 2014: Original implementation based on code provided by
;                    Enci, Lin, and Francesco; by K. Westfall (KBW).
;       22 Jan 2015: (KBW) Fixed bug in how I saved Francesco's
;                          velocities.
;       13 Jul 2015: (KBW) Added output for different weighting and
;                          standard error based on weighted standard
;                          deviation.  Moved calculation of instrumental
;                          dispersion here so that this could be taken
;                          into account when determining the weighted
;                          mean values.
;       07 Aug 2015: (KBW) MDAP_SAVE_EMISSION_LINE_SET_OMITTED is not
;                          used anymore, commented it out.  Added
;                          best_fit_mask as an input parameter; used to
;                          define regions where the emission-line fits
;                          should be omitted because the continuum
;                          subtraction will be incorrect.
;       13 Aug 2015: (KBW) Edited Enci's and Francesco's code to fit the
;                          OII lines as a doublet (instead of a single
;                          ine, and include fits of the OI doublet at
;                          6300 A.  Edited this code accordingly.
;       16 Sep 2015: (KBW) Edited both contributed codes to allow for a
;                          non-zero background; edited this code
;                          accordingly.
;       17 Sep 2015: (KBW) Return the window limits used for the fit.
;       20 Oct 2015: (KBW) Include a fit to the OII doublet as a single
;                          line.  Improved flux guess in Enci's code.
;                          Automatically exclude some lines from the
;                          mean gas kinematics; see
;                          MDAP_AUTO_EXCLUDE_EML_FROM_KIN().
;-
;------------------------------------------------------------------------------

; Decide if the line should be omitted.
; Omit line if:
;   - amplitude is not greater than 0
;   - nearest pixel to line center is masked (mask is greater than 0)
FUNCTION MDAP_SAVE_EMISSION_LINE_FIT_CHECK_OMIT, $
                amplitude, line_center, rest_wave, galaxy_mask
        if amplitude gt 0.0 && galaxy_mask[(sort(abs(rest_wave-line_center)))[0]] lt 0.5 then $
            return, 0
        return, 1
END

FUNCTION MDAP_SAVE_EMISSION_LINE_FIT_ENCI_RESC
        np = 4
        cols = { base:0, base_e:0+np, $
                  cen:1,  cen_e:1+np, $
                  sig:2,  sig_e:2+np, $
                 area:3, area_e:3+np }
        return, cols
END

PRO MDAP_SAVE_EMISSION_LINE_FIT_ENCI_SINGLE, lambda, voff, fitwin, result, cols, fitwindow, $
                                             baseline, baseline_err, kinematics, kinematics_err, $
                                             intens, intens_err, fluxes, fluxes_err, i, j

        c=299792.458d                           ; Speed of light in km/s
;        kinematics[i,j,0] = (result[0]/lambda - 1.0d)*c + voff
;        kinematics_err[i,j,0] = result[3]*c/lambda
;        kinematics[i,j,1] = result[1]*c/result[0]
;        kinematics_err[i,j,1] = kinematics[i,j,1]*sqrt( (result[3]/result[0])^2 + $
;                                                        (result[4]/result[1])^2 )
;        intens[i,j] = result[2]/sqrt(2.0*!pi)/result[1]
;        intens_err[i,j] = intens[i,j]*sqrt( (result[5]/result[2])^2 + (result[4]/result[1])^2 )
;        fluxes[i,j] = result[2]
;        fluxes_err[i,j] = result[5]

        fitwindow[i,j,0] = fitwin[0]*(1.0d + voff/c)
        fitwindow[i,j,1] = fitwin[1]*(1.0d + voff/c)

        baseline[i,j] = result[cols.base]
        baseline_err[i,j] = result[cols.base_e]

        kinematics[i,j,0] = (result[cols.cen]/lambda - 1.0d)*c + voff
        kinematics_err[i,j,0] = result[cols.cen_e]*c/lambda
        kinematics[i,j,1] = result[cols.sig]*c/result[cols.cen]
        kinematics_err[i,j,1] = abs(kinematics[i,j,1]) $
                                        * sqrt( (result[cols.cen_e]/result[cols.cen])^2 $
                                              + (result[cols.sig_e]/result[cols.sig])^2 )

        intens[i,j] = result[cols.area]/sqrt(2.0*!pi)/result[cols.sig]
        intens_err[i,j] = abs(intens[i,j]) $
                                        * sqrt( (result[cols.area_e]/result[cols.area])^2 $
                                              + (result[cols.sig_e]/result[cols.sig])^2 )
        fluxes[i,j] = result[cols.area]
        fluxes_err[i,j] = result[cols.area_e]
END

PRO MDAP_SAVE_EMISSION_LINE_FIT_ENCI, eml_par, rest_wave, galaxy_mask, voff, fitwin, result, $
                                      omitted, fitwindow, baseline, baseline_err, kinematics, $
                                      kinematics_err, intens, intens_err, fluxes, fluxes_err, i

        cols = MDAP_SAVE_EMISSION_LINE_FIT_ENCI_RESC()

        j = 0

        ;---------------------------------------------------------------
        ; OII 3727
;        omitted[i,j] = result.OII3727[2] gt 0.0 ? 0 : 1
        omitted[i,j] = MDAP_SAVE_EMISSION_LINE_FIT_CHECK_OMIT(result.OII3727[cols.area], $
                                                              result.OII3727[cols.cen], rest_wave, $
                                                              galaxy_mask)
;       if omitted[i,j] eq 1 then begin
;           MDAP_SAVE_EMISSION_LINE_SET_OMITTED, kinematics, kinematics_err, intens, intens_err, $
;                                                fluxes, fluxes_err, i, j
;       endif else begin
            MDAP_SAVE_EMISSION_LINE_FIT_ENCI_SINGLE, eml_par[j].lambda, voff, fitwin.OII3727, $
                                                     result.OII3727, cols, fitwindow, baseline, $
                                                     baseline_err, kinematics, kinematics_err, $
                                                     intens, intens_err, fluxes, fluxes_err, i, j
;       endelse
        j = j+1
        ;---------------------------------------------------------------

        ;---------------------------------------------------------------
        ; OII 3729
;        omitted[i,j] = result.OII3727[2] gt 0.0 ? 0 : 1
        omitted[i,j] = MDAP_SAVE_EMISSION_LINE_FIT_CHECK_OMIT(result.OII3729[cols.area], $
                                                              result.OII3729[cols.cen], rest_wave, $
                                                              galaxy_mask)
;       if omitted[i,j] eq 1 then begin
;           MDAP_SAVE_EMISSION_LINE_SET_OMITTED, kinematics, kinematics_err, intens, intens_err, $
;                                                fluxes, fluxes_err, i, j
;       endif else begin
            MDAP_SAVE_EMISSION_LINE_FIT_ENCI_SINGLE, eml_par[j].lambda, voff, fitwin.OII3729, $
                                                     result.OII3729, cols, fitwindow, baseline, $
                                                     baseline_err, kinematics, kinematics_err, $
                                                     intens, intens_err, fluxes, fluxes_err, i, j
;       endelse
        j = j+1
        ;---------------------------------------------------------------

        ;---------------------------------------------------------------
        ; OII 3727 double fit as a single line
;        omitted[i,j] = result.OII3727[2] gt 0.0 ? 0 : 1
        omitted[i,j] = MDAP_SAVE_EMISSION_LINE_FIT_CHECK_OMIT(result.OIId3727[cols.area], $
                                                              result.OIId3727[cols.cen], $
                                                              rest_wave, galaxy_mask)
;       if omitted[i,j] eq 1 then begin
;           MDAP_SAVE_EMISSION_LINE_SET_OMITTED, kinematics, kinematics_err, intens, intens_err, $
;                                                fluxes, fluxes_err, i, j
;       endif else begin
            MDAP_SAVE_EMISSION_LINE_FIT_ENCI_SINGLE, eml_par[j].lambda, voff, fitwin.OIId3727, $
                                                     result.OIId3727, cols, fitwindow, baseline, $
                                                     baseline_err, kinematics, kinematics_err, $
                                                     intens, intens_err, fluxes, fluxes_err, i, j
;       endelse
        j = j+1
        ;---------------------------------------------------------------

        ;---------------------------------------------------------------
        ; H-beta
;        omitted[i,j] = result.Hb4861[2] gt 0.0 ? 0 : 1
        omitted[i,j] = MDAP_SAVE_EMISSION_LINE_FIT_CHECK_OMIT(result.Hb4861[cols.area], $
                                                              result.Hb4861[cols.cen], rest_wave, $
                                                              galaxy_mask)
;       if omitted[i,j] eq 1 then begin
;           MDAP_SAVE_EMISSION_LINE_SET_OMITTED, kinematics, kinematics_err, intens, intens_err, $
;                                                fluxes, fluxes_err, i, j
;       endif else begin
            MDAP_SAVE_EMISSION_LINE_FIT_ENCI_SINGLE, eml_par[j].lambda, voff, fitwin.Hb4861, $
                                                     result.Hb4861, cols, fitwindow, baseline, $
                                                     baseline_err, kinematics, kinematics_err, $
                                                     intens, intens_err, fluxes, fluxes_err, i, j
;       endelse
        j = j+1
        ;---------------------------------------------------------------

        ;---------------------------------------------------------------
        ; OIII 4959
;       omitted[i,j] = result.OIII4959[2] gt 0.0 ? 0 : 1
        omitted[i,j] = MDAP_SAVE_EMISSION_LINE_FIT_CHECK_OMIT(result.OIII4959[cols.area], $
                                                              result.OIII4959[cols.cen], $
                                                              rest_wave, galaxy_mask)
;       if omitted[i,j] eq 1 then begin
;           MDAP_SAVE_EMISSION_LINE_SET_OMITTED, kinematics, kinematics_err, intens, intens_err, $
;                                                fluxes, fluxes_err, i, j
;       endif else begin
            MDAP_SAVE_EMISSION_LINE_FIT_ENCI_SINGLE, eml_par[j].lambda, voff, fitwin.OIII4959, $
                                                     result.OIII4959, cols, fitwindow, baseline, $
                                                     baseline_err, kinematics, kinematics_err, $
                                                     intens, intens_err, fluxes, fluxes_err, i, j
;       endelse
        j = j+1
        ;---------------------------------------------------------------

        ;---------------------------------------------------------------
        ; OIII 5007
;        omitted[i,j] = result.OIII5007[2] gt 0.0 ? 0 : 1
        omitted[i,j] = MDAP_SAVE_EMISSION_LINE_FIT_CHECK_OMIT(result.OIII5007[cols.area], $
                                                              result.OIII5007[cols.cen], $
                                                              rest_wave, galaxy_mask)
;       if omitted[i,j] eq 1 then begin
;           MDAP_SAVE_EMISSION_LINE_SET_OMITTED, kinematics, kinematics_err, intens, intens_err, $
;                                                fluxes, fluxes_err, i, j
;       endif else begin
            MDAP_SAVE_EMISSION_LINE_FIT_ENCI_SINGLE, eml_par[j].lambda, voff, fitwin.OIII5007, $
                                                     result.OIII5007, cols, fitwindow, baseline, $
                                                     baseline_err, kinematics, kinematics_err, $
                                                     intens, intens_err, fluxes, fluxes_err, i, j
;       endelse
        j = j+1
        ;---------------------------------------------------------------

        ;---------------------------------------------------------------
        ; OI 6300
        omitted[i,j] = MDAP_SAVE_EMISSION_LINE_FIT_CHECK_OMIT(result.OI6300[cols.area], $
                                                              result.OI6300[cols.cen], $
                                                              rest_wave, galaxy_mask)
        MDAP_SAVE_EMISSION_LINE_FIT_ENCI_SINGLE, eml_par[j].lambda, voff, fitwin.OI6300, $
                                                 result.OI6300, cols, fitwindow, baseline, $
                                                 baseline_err, kinematics, kinematics_err, intens, $
                                                 intens_err, fluxes, fluxes_err, i, j
        j = j+1
        ;---------------------------------------------------------------
        
        ;---------------------------------------------------------------
        ; OI 6363
        omitted[i,j] = MDAP_SAVE_EMISSION_LINE_FIT_CHECK_OMIT(result.OI6363[cols.area], $
                                                              result.OI6363[cols.cen], rest_wave, $
                                                              galaxy_mask)
        MDAP_SAVE_EMISSION_LINE_FIT_ENCI_SINGLE, eml_par[j].lambda, voff, fitwin.OI6363, $
                                                 result.OI6363, cols, fitwindow, baseline, $
                                                 baseline_err, kinematics, kinematics_err, intens, $
                                                 intens_err, fluxes, fluxes_err, i, j
        j = j+1
        ;---------------------------------------------------------------
        

        ;---------------------------------------------------------------
        ; NII 6548
;        omitted[i,j] = result.NII6548[2] gt 0.0 ? 0 : 1
        omitted[i,j] = MDAP_SAVE_EMISSION_LINE_FIT_CHECK_OMIT(result.NII6548[cols.area], $
                                                              result.NII6548[cols.cen], rest_wave, $
                                                              galaxy_mask)
;       if omitted[i,j] eq 1 then begin
;           MDAP_SAVE_EMISSION_LINE_SET_OMITTED, kinematics, kinematics_err, intens, intens_err, $
;                                                fluxes, fluxes_err, i, j
;       endif else begin
            MDAP_SAVE_EMISSION_LINE_FIT_ENCI_SINGLE, eml_par[j].lambda, voff, fitwin.NII6548, $
                                                     result.NII6548, cols, fitwindow, baseline, $
                                                     baseline_err, kinematics, kinematics_err, $
                                                     intens, intens_err, fluxes, fluxes_err, i, j
;       endelse
        j = j+1
        ;---------------------------------------------------------------

        ;---------------------------------------------------------------
        ; H-alpha
;        omitted[i,j] = result.Ha6563[2] gt 0.0 ? 0 : 1
        omitted[i,j] = MDAP_SAVE_EMISSION_LINE_FIT_CHECK_OMIT(result.Ha6563[cols.area], $
                                                              result.Ha6563[cols.cen], rest_wave, $
                                                              galaxy_mask)
;       if omitted[i,j] eq 1 then begin
;           MDAP_SAVE_EMISSION_LINE_SET_OMITTED, kinematics, kinematics_err, intens, intens_err, $
;                                                fluxes, fluxes_err, i, j
;       endif else begin
            MDAP_SAVE_EMISSION_LINE_FIT_ENCI_SINGLE, eml_par[j].lambda, voff, fitwin.Ha6563, $
                                                     result.Ha6563, cols, fitwindow, baseline, $
                                                     baseline_err, kinematics, kinematics_err, $
                                                     intens, intens_err, fluxes, fluxes_err, i, j
;       endelse
        j = j+1
        ;---------------------------------------------------------------

        ;---------------------------------------------------------------
        ; NII 6583
;        omitted[i,j] = result.NII6583[2] gt 0.0 ? 0 : 1
        omitted[i,j] = MDAP_SAVE_EMISSION_LINE_FIT_CHECK_OMIT(result.NII6583[cols.area], $
                                                              result.NII6583[cols.cen], rest_wave, $
                                                              galaxy_mask)
;       if omitted[i,j] eq 1 then begin
;           MDAP_SAVE_EMISSION_LINE_SET_OMITTED, kinematics, kinematics_err, intens, intens_err, $
;                                                fluxes, fluxes_err, i, j
;       endif else begin
            MDAP_SAVE_EMISSION_LINE_FIT_ENCI_SINGLE, eml_par[j].lambda, voff, fitwin.NII6583, $
                                                     result.NII6583, cols, fitwindow, baseline, $
                                                     baseline_err, kinematics, kinematics_err, $
                                                     intens, intens_err, fluxes, fluxes_err, i, j
;       endelse
        j = j+1
        ;---------------------------------------------------------------

        ;---------------------------------------------------------------
        ; SII 6717
;        omitted[i,j] = result.SII6717[2] gt 0.0 ? 0 : 1
        omitted[i,j] = MDAP_SAVE_EMISSION_LINE_FIT_CHECK_OMIT(result.SII6717[cols.area], $
                                                              result.SII6717[cols.cen], rest_wave, $
                                                              galaxy_mask)
;       if omitted[i,j] eq 1 then begin
;           MDAP_SAVE_EMISSION_LINE_SET_OMITTED, kinematics, kinematics_err, intens, intens_err, $
;                                                fluxes, fluxes_err, i, j
;       endif else begin
            MDAP_SAVE_EMISSION_LINE_FIT_ENCI_SINGLE, eml_par[j].lambda, voff, fitwin.SII6717, $
                                                     result.SII6717, cols, fitwindow, baseline, $
                                                     baseline_err, kinematics, kinematics_err, $
                                                     intens, intens_err, fluxes, fluxes_err, i, j
;       endelse
        j = j+1
        ;---------------------------------------------------------------

        ;---------------------------------------------------------------
        ; SII 6731
;        omitted[i,j] = result.SII6731[2] gt 0.0 ? 0 : 1
        omitted[i,j] = MDAP_SAVE_EMISSION_LINE_FIT_CHECK_OMIT(result.SII6731[cols.area], $
                                                              result.SII6731[cols.cen], rest_wave, $
                                                              galaxy_mask)
;       if omitted[i,j] eq 1 then begin
;           MDAP_SAVE_EMISSION_LINE_SET_OMITTED, kinematics, kinematics_err, intens, intens_err, $
;                                                fluxes, fluxes_err, i, j
;       endif else begin
            MDAP_SAVE_EMISSION_LINE_FIT_ENCI_SINGLE, eml_par[j].lambda, voff, fitwin.SII6731, $
                                                     result.SII6731, cols, fitwindow, baseline, $
                                                     baseline_err, kinematics, kinematics_err, $
                                                     intens, intens_err, fluxes, fluxes_err, i, j
;       endelse
        j = j+1
        ;---------------------------------------------------------------

END

;PRO MDAP_SAVE_EMISSION_LINE_SET_OMITTED, kinematics, kinematics_err, intens, intens_err, fluxes, $
;                                         fluxes_err, i, j
;        kinematics[i,j,0:1] = 0.0d
;        kinematics_err[i,j,0:1] = 1.0d
;        intens[i,j] = 0.0d
;        intens_err[i,j] = 1.0d
;        fluxes[i,j] = 0.0d
;        fluxes_err[i,j] = 1.0d
;END

; CERTAIN ORDER EXPECTED FOR BELFIORE RESULT
PRO MDAP_SAVE_EMISSION_LINE_FIT_BELFIORE, eml_par, voff, result, rest_wave, galaxy_mask, omitted, $
                                          fitwindow, baseline, baseline_err, kinematics, $
                                          kinematics_err, intens, intens_err, fluxes, fluxes_err, i

        c=299792.458d                           ; Speed of light in km/s
        neml = n_elements(eml_par)              ; Number of emission lines fitted

        for j=0,neml-1 do begin

;           omitted[i,j] = result.Ampl[j] gt 0.0 ? 0 : 1

            line_center = eml_par[j].lambda * (1+result.Vel[j]/c)
            omitted[i,j] = MDAP_SAVE_EMISSION_LINE_FIT_CHECK_OMIT(result.Ampl[j], line_center, $
                                                                  rest_wave, galaxy_mask)

;           if omitted[i,j] eq 1 then begin
;               MDAP_SAVE_EMISSION_LINE_SET_OMITTED, kinematics, kinematics_err, intens, $
;                                                    intens_err, fluxes, fluxes_err, i, j
;               continue
;           endif

;           fitwindow[i,j,0] = result.Lmin[j] * (1.0d + (result.Vel[j]+voff)/c)
;           fitwindow[i,j,1] = result.Lmax[j] * (1.0d + (result.Vel[j]+voff)/c)

            ; Does not include fitted velocity, just redshift offset
            fitwindow[i,j,0] = result.Lmin[j] * (1.0d + voff/c)
            fitwindow[i,j,1] = result.Lmax[j] * (1.0d + voff/c)

            baseline[i,j] = result.Base[j]
            baseline_err[i,j] = result.eBase[j]

            kinematics[i,j,0] = result.Vel[j]+voff
            kinematics_err[i,j,0] = result.eVel[j]
            kinematics[i,j,1] = result.Sigma[j]
            kinematics_err[i,j,1] = result.eSigma[j]

            intens[i,j] = result.Ampl[j]
            intens_err[i,j] = result.eAmpl[j]

;           fluxes[i,j] = result.Ampl[j] * sqrt(2*!pi) * result.Sigma[j] * eml_par[j].lambda $
;                         * (1.0+(result.Vel[j]+voff)/c)/c

            ; Does not include redshift offset, just fitted velocity
            fluxes[i,j] = result.Ampl[j] * sqrt(2*!pi) * result.Sigma[j] * eml_par[j].lambda $
                          * (1.0+result.Vel[j]/c)/c

            fluxes_err[i,j] = sqrt( (fluxes[i,j]*result.eAmpl[j]/result.Ampl[j])^2 + $
                                    (fluxes[i,j]*result.eSigma[j]/result.Sigma[j])^2 + $
                                    (result.Ampl[j] * sqrt(2*!pi) * result.Sigma[j] * $
                                     eml_par[j].lambda * result.eVel[j]/c/c)^2 )
        endfor
END


PRO MDAP_EMISSION_LINE_ONLY_FIT, $
                wave, sres, flux, ivar, mask, best_fit_continuum, best_fit_mask, $
                stellar_kinematics, eml_model, emission_line_kinematics, $
                emission_line_kinematics_perr, emission_line_kinematics_serr, $
                emission_line_kinematics_n, emission_line_omitted, emission_line_window, $
                emission_line_baseline, emission_line_base_err, $
                emission_line_kinematics_individual, emission_line_kinematics_individual_err, $
                emission_line_sinst, emission_line_intens, emission_line_intens_err, $
                emission_line_fluxes, emission_line_fluxes_err, emission_line_EW, $
                emission_line_EW_err, eml_par=eml_par, version=version, enci=enci, $
                belfiore=belfiore, zero_base=zero_base, quiet=quiet, dbg=dbg

        version_module = '1.2'                  ; Version number
        if n_elements(version) ne 0 then begin  ; If version is defined
            version = version_module            ; ... set it to the module version
            return                              ; ... and return without doing anything
        endif

        if keyword_set(enci) && keyword_set(belfiore) then begin
            message, 'Cannot run both Enci and Belfiore fits simultaneously.  Make separate ' $
                     + 'calls to MDAP_EMISSION_LINE_ONLY_FIT'
        endif
        if ~keyword_set(enci) && ~keyword_set(belfiore) && ~keyword_set(quiet) then $
            print, 'No mode selected for MDAP_EMISSION_LINE_ONLY_FIT.  Default is /enci.'

        ; Get the emission line parameters for the lines fit by Enci's
        ; code, in the format of the structure used by GANDALF
        eml_par = MDAP_DEFINE_EMISSION_LINES_ENCI_BELFIORE()

        ; Convert the inverse variance to an error, adjusting the mask for bad values
        indx = where(ivar gt 0.0d, count, complement=nindx, ncomplement=ncnt)
;       print, count, ncnt
;       stop
        if count eq 0 then $
            message, 'All pixels are masked!'

        flux_err = ivar
        flux_err[indx] = sqrt(1.0d/ivar[indx])
        if ncnt ne 0 then begin
            flux_err[nindx] = 1.0d
            mask[nindx] = 1.0d
        endif

        ; Allocate the output arrays
        c=299792.458d                           ; Speed of light in km/s
        sz = size(flux)
        nobj = sz[1]                            ; Number of object spectra
        npix = sz[2]                            ; Number of pixels in the spectra

        neml = (size(eml_par))[1]               ; Number of fitted lines

        eml_model = dblarr(nobj,npix)           ; Best-fitting emission-line model

        nstat = 4                           ; Number of ways the weighted kinematics are determined

        emission_line_kinematics = dblarr(nobj, 2*nstat)
        emission_line_kinematics_perr = dblarr(nobj, 2*nstat)
        emission_line_kinematics_serr = dblarr(nobj, 2*nstat)
        emission_line_kinematics_n = intarr(nobj)
        emission_line_omitted = make_array(nobj, neml, /int, value=1)       ; 0/1 -> fit/omitted
        emission_line_window = dblarr(nobj, neml, 2)
        emission_line_baseline = dblarr(nobj, neml)
        emission_line_base_err = dblarr(nobj, neml)
        emission_line_kinematics_individual = dblarr(nobj, neml, 2)
        emission_line_kinematics_individual_err = dblarr(nobj, neml, 2)
        emission_line_sinst = dblarr(nobj, neml)
        emission_line_intens = dblarr(nobj, neml)
        emission_line_intens_err = dblarr(nobj, neml)
        emission_line_fluxes = dblarr(nobj, neml)
        emission_line_fluxes_err = dblarr(nobj, neml)
        emission_line_EW = dblarr(nobj, neml)
        emission_line_EW_err = dblarr(nobj, neml)

        mean_kin_auto_exclude = MDAP_AUTO_EXCLUDE_EML_FROM_KIN(eml_par, belfiore=belfiore)

        ; TODO: Calculate a chi-square for the full fit using the
        ; eml_model
   
        ; Fit each spectrum and save the results
        if keyword_set(dbg) then $
            nobj = 1
        for i=0,nobj-1 do begin

            ; Set the wavelength vector to the rest wavelength
            ; TODO: Could also just pass velocity to MDAP_FIT_EMISSION_LINE_SPECTRUM_ENCI
            rest_wave = wave/(1+stellar_kinematics[i,0]/c)

            ; Remove the stellar continuum
            galaxy_eml_only = reform(flux[i,*] - best_fit_continuum[i,*])
            galaxy_mask = reform(best_fit_mask[i,*])

            if keyword_set(belfiore) then begin
                ; Fit the emission lines
                MDAP_FIT_EMISSION_LINE_SPECTRUM_BELFIORE, rest_wave, galaxy_eml_only, $
                                                          reform(flux_err[i,*]), $
                                                          reform(mask[i,*]), $
                                                          stellar_kinematics[i,1], 0.0d, result, $
                                                          emfit, zero_base=zero_base, /quiet
                                                          ; quiet=quiet

                ; Save the fitted parameters
                MDAP_SAVE_EMISSION_LINE_FIT_BELFIORE, eml_par, stellar_kinematics[i,0], result, $
                                                      rest_wave, galaxy_mask, $
                                                      emission_line_omitted, emission_line_window, $
                                                      emission_line_baseline, $
                                                      emission_line_base_err, $
                                                      emission_line_kinematics_individual, $
                                                      emission_line_kinematics_individual_err, $
                                                      emission_line_intens, $
                                                      emission_line_intens_err, $
                                                      emission_line_fluxes, $
                                                      emission_line_fluxes_err, i
            endif else begin
                ; Fit the emission lines
                MDAP_FIT_EMISSION_LINE_SPECTRUM_ENCI, 0.0d, rest_wave, galaxy_eml_only, $
                                                      reform(flux_err[i,*]), reform(mask[i,*]), $
                                                      fitwin, result, emfit, /OII3727, /OIId3727, $
                                                      /Hb4861, /OIII5007, /OI6300, /Ha6563, $
                                                      /SII6717, zero_base=zero_base

                ; Save the fitted parameters
                MDAP_SAVE_EMISSION_LINE_FIT_ENCI, eml_par, rest_wave, galaxy_mask, $
                                                  stellar_kinematics[i,0], fitwin, result, $
                                                  emission_line_omitted, emission_line_window, $
                                                  emission_line_baseline, $
                                                  emission_line_base_err, $
                                                  emission_line_kinematics_individual, $
                                                  emission_line_kinematics_individual_err, $
                                                  emission_line_intens, emission_line_intens_err, $
                                                  emission_line_fluxes, emission_line_fluxes_err, i
            endelse

            ; Save the best-fitting model
            eml_model[i,*] = emfit
                
            ; Calculate the equivalent width and error
            galaxy_no_eml = reform(flux[i,*]) - emfit
            MDAP_MEASURE_EMISSION_LINE_EQUIV_WIDTH, eml_par, $
                                                    emission_line_kinematics_individual[i,*,0], $
                                                    emission_line_kinematics_individual[i,*,1], $
                                                    rest_wave, galaxy_no_eml, $
                                                    emission_line_fluxes[i,*], $
                                                    emission_line_fluxes_err[i,*], equiv_width, $
                                                    equiv_width_err
            emission_line_EW[i,*] = equiv_width
            emission_line_EW_err[i,*] = equiv_width_err

            ; Calculate the instrumental dispersion
            emission_line_sinst[i,*] = MDAP_INSTRUMENTAL_DISPERSION(wave, sres, eml_par[*].lambda, $
                                                        emission_line_kinematics_individual[i,*,0])

            ; Decide which lines to inlude in the mean measurement of the kinematics
            indx=MDAP_LINES_FOR_MEAN_GAS_KINEMATICS(emission_line_intens[i,*], $
                                                    emission_line_intens_err[i,*], $
                                                    emission_line_fluxes[i,*], $
                                                    emission_line_fluxes_err[i,*], $
                                                    emission_line_kinematics_individual[i,*,0], $
                                                    emission_line_kinematics_individual_err[i,*,0],$
                                                    emission_line_kinematics_individual[i,*,1], $
                                                    emission_line_kinematics_individual_err[i,*,1],$
                                                    exclude=mean_kin_auto_exclude, $
                                                    omitted=emission_line_omitted[i,*], count=count)

            ; Convert the measured velocity dispersions to the square
            ; and correct that value for the instrumental dispersion.
            ; This is done after the check above in case lines to
            ; include in the average depend on criteria based on the
            ; observe dispersion, not the corrected one.
            emission_line_kinematics_individual_err[i,*,1] = $
                                            2 * emission_line_kinematics_individual[i,*,1] $
                                              * emission_line_kinematics_individual_err[i,*,1]
            emission_line_kinematics_individual[i,*,1] = $
                                            (emission_line_kinematics_individual[i,*,1])^2 $
                                             - (emission_line_sinst[i,*])^2

            ; Get the mean kinematics
            emission_line_kinematics_n[i] = count
            if count eq 0 then begin

                ; Set the mean kinematics to placeholder values
                emission_line_kinematics[i,*] = -9999.0d
                emission_line_kinematics_perr[i,*] = -9999.0d
                emission_line_kinematics_serr[i,*] = -9999.0d

            endif else begin

                ; TODO:
                ; THE DOUBLET LINES FROM FRANCESCO'S CODE SHOULD NOT BE
                ; TREATED AS IF THEY PROVIDE INDEPENDENT MEASUREMENTS OF
                ; THE VELOCITY!!

                ;-------------------------------------------------------
                ; flux-weighted
                MDAP_MEAN_GAS_KINEMATICS, emission_line_fluxes[i,indx], $
                                          emission_line_kinematics_individual[i,indx,0], $
                                          emission_line_kinematics_individual_err[i,indx,0], $
                                          emission_line_kinematics_individual[i,indx,1], $
                                          emission_line_kinematics_individual_err[i,indx,1], $
                                          fw_vel, fw_vel_err, fw_vel_serr, fw_sig, fw_sig_err, $
                                          fw_sig_serr, /flux_weighted

                ; Convert sigma^2 quantities to sigma, allowing for
                ; negative values
                fw_sig = fw_sig / sqrt(abs(fw_sig))
                if (abs(fw_sig) > 0.0) then begin
                    fw_sig_err = fw_sig_err/2.0/abs(fw_sig)
                    fw_sig_serr = fw_sig_serr/2.0/abs(fw_sig)
                endif else begin
                    fw_sig_err = -9999.0d
                    fw_sig_serr = -9999.0d
                endelse

                emission_line_kinematics[i,0] = fw_vel              ; Weighted-mean velocity
                emission_line_kinematics_perr[i,0] = fw_vel_err     ; propagated error
                emission_line_kinematics_serr[i,0] = fw_vel_serr    ; standard error
                emission_line_kinematics[i,1] = fw_sig              ; Weighted-mean velocity disp
                emission_line_kinematics_perr[i,1] = fw_sig_err     ; propagated error
                emission_line_kinematics_serr[i,1] = fw_sig_serr    ; standard error

                ;-------------------------------------------------------
                ; velocity-error-weighted
                MDAP_MEAN_GAS_KINEMATICS, emission_line_fluxes[i,indx], $
                                          emission_line_kinematics_individual[i,indx,0], $
                                          emission_line_kinematics_individual_err[i,indx,0], $
                                          emission_line_kinematics_individual[i,indx,1], $
                                          emission_line_kinematics_individual_err[i,indx,1], $
                                          fw_vel, fw_vel_err, fw_vel_serr, fw_sig, fw_sig_err, $
                                          fw_sig_serr, /velocity_error_weighted

                ; Convert sigma^2 quantities to sigma, allowing for
                ; negative values
                fw_sig = fw_sig / sqrt(abs(fw_sig))
                if (abs(fw_sig) > 0.0) then begin
                    fw_sig_err = fw_sig_err/2.0/abs(fw_sig)
                    fw_sig_serr = fw_sig_serr/2.0/abs(fw_sig)
                endif else begin
                    fw_sig_err = -9999.0d
                    fw_sig_serr = -9999.0d
                endelse

                emission_line_kinematics[i,2] = fw_vel              ; Weighted-mean velocity
                emission_line_kinematics_perr[i,2] = fw_vel_err     ; propagated error
                emission_line_kinematics_serr[i,2] = fw_vel_serr    ; standard error
                emission_line_kinematics[i,3] = fw_sig              ; Weighted-mean velocity disp
                emission_line_kinematics_perr[i,3] = fw_sig_err     ; propagated error
                emission_line_kinematics_serr[i,3] = fw_sig_serr    ; standard error

                ;-------------------------------------------------------
                ; flux- and velocity-error-weighted
                MDAP_MEAN_GAS_KINEMATICS, emission_line_fluxes[i,indx], $
                                          emission_line_kinematics_individual[i,indx,0], $
                                          emission_line_kinematics_individual_err[i,indx,0], $
                                          emission_line_kinematics_individual[i,indx,1], $
                                          emission_line_kinematics_individual_err[i,indx,1], $
                                          fw_vel, fw_vel_err, fw_vel_serr, fw_sig, fw_sig_err, $
                                          fw_sig_serr, /flux_weighted, /velocity_error_weighted

                ; Convert sigma^2 quantities to sigma, allowing for
                ; negative values
                fw_sig = fw_sig / sqrt(abs(fw_sig))
                if (abs(fw_sig) > 0.0) then begin
                    fw_sig_err = fw_sig_err/2.0/abs(fw_sig)
                    fw_sig_serr = fw_sig_serr/2.0/abs(fw_sig)
                endif else begin
                    fw_sig_err = -9999.0d
                    fw_sig_serr = -9999.0d
                endelse

                emission_line_kinematics[i,4] = fw_vel              ; Weighted-mean velocity
                emission_line_kinematics_perr[i,4] = fw_vel_err     ; propagated error
                emission_line_kinematics_serr[i,4] = fw_vel_serr    ; standard error
                emission_line_kinematics[i,5] = fw_sig              ; Weighted-mean velocity disp
                emission_line_kinematics_perr[i,5] = fw_sig_err     ; propagated error
                emission_line_kinematics_serr[i,5] = fw_sig_serr    ; standard error

                ;-------------------------------------------------------
                ; uniform weighting
                MDAP_MEAN_GAS_KINEMATICS, emission_line_fluxes[i,indx], $
                                          emission_line_kinematics_individual[i,indx,0], $
                                          emission_line_kinematics_individual_err[i,indx,0], $
                                          emission_line_kinematics_individual[i,indx,1], $
                                          emission_line_kinematics_individual_err[i,indx,1], $
                                          fw_vel, fw_vel_err, fw_vel_serr, fw_sig, fw_sig_err, $
                                          fw_sig_serr

                ; Convert sigma^2 quantities to sigma, allowing for
                ; negative values
                fw_sig = fw_sig / sqrt(abs(fw_sig))
                if (abs(fw_sig) > 0.0) then begin
                    fw_sig_err = fw_sig_err/2.0/abs(fw_sig)
                    fw_sig_serr = fw_sig_serr/2.0/abs(fw_sig)
                endif else begin
                    fw_sig_err = -9999.0d
                    fw_sig_serr = -9999.0d
                endelse

                emission_line_kinematics[i,6] = fw_vel              ; Weighted-mean velocity
                emission_line_kinematics_perr[i,6] = fw_vel_err     ; propagated error
                emission_line_kinematics_serr[i,6] = fw_vel_serr    ; standard error
                emission_line_kinematics[i,7] = fw_sig              ; Weighted-mean velocity disp
                emission_line_kinematics_perr[i,7] = fw_sig_err     ; propagated error
                emission_line_kinematics_serr[i,7] = fw_sig_serr    ; standard error

            endelse

            ; Convert sigma^2 quantities to sigma, allowing for negative
            ; values
            indx = where(abs(emission_line_kinematics_individual[i,*,1]) gt 0.0, count, $
                         complement=nindx, ncomplement=ncount)
            if count gt 0 then begin
                emission_line_kinematics_individual[i,indx,1] = $
                                        emission_line_kinematics_individual[i,indx,1] $
                                        / sqrt(abs(emission_line_kinematics_individual[i,indx,1]))
                emission_line_kinematics_individual_err[i,indx,1] = $
                                        emission_line_kinematics_individual_err[i,indx,1] $
                                        / 2.0 / abs(emission_line_kinematics_individual[i,indx,1])
            endif
            if ncount gt 0 then begin
                emission_line_kinematics_individual[i,nindx,1] = 0.0d
                emission_line_kinematics_individual_err[i,nindx,1] = -9999.0d
            endif

        endfor

END



