;+
; NAME:
;       MDAP_CALCULATE_BIN_SN
;
; PURPOSE:
;       Calculate the signal-to-noise of spectra in a bin.  All signal
;       and noise values are used for the calculation.
;
;       For uniform weighting, the nominal signal-to-noise is:
;
;           S/N = Sum(s)/sqrt(Sum(n^2))
;
;       For optimal weighting, the nominal signal-to-noise is:
;
;           S/N = sqrt(Sum((s/n)^2))
;
;       When noise_calib is set the code calls MDAP_CALIBRATE_NOISE to
;       apply (or not apply depending on the value of noise_calib) a
;       calibration of the noise vector to account for covariance among
;       the spaxels.
;
; CALLING SEQUENCE:
;       result = MDAP_CALCULATE_BIN_SN(signal[indx], noise[indx], noise_calib=noise_calib, $
;                                      /optimal_weighting)
;
; INPUTS:
;       signal dblarr[N]
;               Signal in each of N spectra.
;
;       noise dblarr[N]
;               Noise in each of N spectra.
;
; OPTIONAL INPUTS:
;       noise_calib int
;               Apply a calibration of the signal-to-noise that accounts
;               approximately accounts for the covariance between
;               values.  See the equation above determined for MaNGA.
;
; OPTIONAL KEYWORDS:
;       /optimal_weighting
;               Calculate the S/N assuming the values are weighted
;               optimally.
; 
; OUTPUT:
;       Returns the nominal or calibrated average S/N per pixel.
;
; PROCEDURES CALLED:
;       MDAP_CALIBRATE_SN
;
; REVISION HISTORY:
;       04 Dec 2014: Pulled from MDAP_VORONOI_2D_BINNING by K. Westfall
;                    (KBW), include tag for optimal weighting
;       16 Mar 2015: (KBW) Changed to allow for noise_calib instead of
;                          PDAP-based S/N calibration.
;-
;------------------------------------------------------------------------------

FUNCTION MDAP_CALCULATE_BIN_SN, $
                signal, noise, noise_calib=noise_calib, optimal_weighting=optimal_weighting

        if keyword_set(optimal_weighting) then begin
            if n_elements(noise_calib) ne 0 && noise_calib ne 0 then $
                message, 'Cannot optimally weight S/N and calculate noise calibration.'
            return, sqrt(total((signal/noise)^2))       ; eqn (3) Cappellari & Copin (2003)
        endif

        sum_noise = sqrt(total(noise^2))
        if n_elements(noise_calib) ne 0 then $
            sum_noise = mdap_calibrate_noise(sum_noise, n_elements(noise), noise_calib)

        return, total(signal)/sum_noise                 ; eqn (2) Cappellari & Copin (2003)

END

