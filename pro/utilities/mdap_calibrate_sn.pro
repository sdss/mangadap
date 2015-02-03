;+
; NAME:
;       MDAP_CALIBRATE_SN
;
; PURPOSE:
;       Convert the nominal calculation of the S/N to a value that
;       should account of effects that were not included in the formal
;       calculation of the errors.  The calibration is done according to:
;
;               S/N = c_i*(S/N_nom^c_0/sqrt(npix))^(i-1)
;
;       where i=1,C-1.
;
; CALLING SEQUENCE:
;       result = MDAP_CALIBRATE_SN(sn_nom, npix, coeff)
;
; INPUTS:
;       sn_nom double
;               The nominal calculation of the signal-to-noise:
;
;                   sn_nom = total(signal)/sqrt(total(noise^2))
;
;       npix integer
;               Number of pixel used in the signal-to-noise calculation.
;
;       coeff dblarr[C]
;           A list of C coefficients needed for the calculation of the
;           calibrated S/N; this is the c vector described in the
;           PURPOSE summary.
;
; OUTPUT:
;       Result is the calibrated S/N.
;
; REVISION HISTORY:
;       04 Dec 2014: (KBW) Copied from the L. Coccato's version.
;-
;------------------------------------------------------------------------------

FUNCTION MDAP_CALIBRATE_SN, $
                sn_nom, npix, coeff
        x = sn_nom^coeff[0]/sqrt(npix)
        return, poly(x, coeff[1:*])
END


