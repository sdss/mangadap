;+
; NAME:
;       MDAP_CALCULATE_BIN_SN
;
; PURPOSE:
;       Calculate the signal-to-noise of spectra in a bin.  All signal
;       and noise values are used for the calculation.  If a set of
;       calibration coefficients are provided, they are used to alter
;       the S/N from the nominal value.
;
;       For uniform weighting, the nominal signal-to-noise is:
;
;           S/N = Sum(s)/sqrt(Sum(n^2))
;
;       For optimal weighting, the nominal signal-to-noise is:
;
;           S/N = sqrt(Sum((s/n)^2))
;
;
; CALLING SEQUENCE:
;       result = MDAP_CALCULATE_BIN_SN(signal[indx], noise[indx], sn_calibration=sn_calibration)
;
; INPUTS:
;       signal dblarr[N]
;               Signal in each of N spectra.
;
;       noise dblarr[N]
;               Noise in each of N spectra.
;
; OPTIONAL INPUTS:
;       sn_calibration dblarr[C]
;               Set of C coefficents used for the calibration of the S/N
;               measurement.  See MDAP_CALIBRATE_SN().
;
; OUTPUT:
;       Returns the nominal or calibrated average S/N per pixel.
;
; PROCEDURES CALLED:
;       MDAP_CALIBRATE_SN
;
; REVISION HISTORY:
;       04 Dec 2014: (KBW) Pulled from MDAP_VORONOI_2D_BINNING, include
;                          tag for optimal weighting
;-
;------------------------------------------------------------------------------

FUNCTION MDAP_CALCULATE_BIN_SN, $
                signal, noise, sn_calibration=sn_calibration, optimal_weighting=optimal_weighting

        if keyword_set(optimal_weighting) then begin
            sn_est = sqrt(total((signal/noise)^2))          ; eqn (3) Cappellari & Copin (2003)
        endif else $
            sn_est = total(signal)/sqrt(total(noise^2))     ; eqn (2) Cappellari & Copin (2003)

        if n_elements(sn_calibration) eq 0 then $
            return, sn_est

        return, mdap_calibrate_sn(sn_est, n_elements(signal), sn_calibration)
END

