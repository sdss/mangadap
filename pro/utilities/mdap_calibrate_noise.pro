;+
; NAME:
;       MDAP_CALIBRATE_NOISE
;
; PURPOSE:
;       Calibrate a noise value/vector based on a calibration of the
;       noise to account for covariance.  The coefficient is HARD-WIRED
;       here based on a calibration done for MaNGA data cubes.  The
;       results are explained here:
;
;       https://trac.sdss.org/wiki/MANGA/Projects/Covariance_cube#SNinSpatiallyBinnedSpectra
;
;       There, we show:
;
;       S/N_calib = S/N_{no covar} / (1 + 1.62*log10(N_bin))
;
;       where N_bin is the number of binned spaxels.
;
; CALLING SEQUENCE:
;       result = MDAP_CALIBRATE_NOISE(noise_nocovar, nbin, noise_calib)
;
; INPUTS:
;       noise_nocovar double or dblarr
;               Nominal calculation of the noise that does not include
;               covariance.
;
;       n_bin int
;               Number of bins that have been combined.
;
;       noise_calib int
;               Type of calibration to apply:
;                   0 - No calibration
;                   1 - Calibration given above
;
; OUTPUT:
;       Result is the calibrated noise value/vector.
;
; REVISION HISTORY:
;       16 Mar 2015: Original implementation by K. Westfall (KBW)
;-
;------------------------------------------------------------------------------

FUNCTION MDAP_CALIBRATE_NOISE, $
                noise_nocovar, nbin, noise_calib
        if noise_calib eq 0 then $
            return, noise_nocovar

        return, noise_nocovar * (1.0 + 1.62 * alog10(nbin))
END


