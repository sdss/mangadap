;+
; NAME:
;       MDAP_REVERT_PIXEL_KINEMATICS
;
; PURPOSE:
;       Revert kinematics from redshifts back to the pixel shifts from
;       PPXF/GANDALF.  PPXF/GANDALF determines the velocity offset by making the
;       approximation that every pixel (logarithmically binned in wavelength) is
;       a constant change in velocity.  Over large velocity shifts, this
;       approximation can become poor.  An e-mail summary from Michele
;       Cappellari:
;
;
;           "The velocity scale per pixel is input as
;
;               velScale = c*logLamScale = c*(logLam[1] - logLam[0])
;       
;           so the output velocity becomes
;
;               vel = c*ln(lam1/lam0)
;
;           which implies that the relation between PPXF output velocity and
;           redshift is
;
;               1 + z = exp(vel/c)
;
;           which reduces in the low-redshift case to z~vel/c."
;
;
;       Put another way by Renbin Yan:
;
;           "We should be writing 1+z = exp(pixelshift*pixelwidth/c)."
;
;       This function reverts the cz velocities to the 'vel' values expected by
;       PPXF/GANDALF.
;
;       **IMPORTANT**: After the the kinematics have been converted to cz using
;       MDAP_CONVERT_PIXEL_KINEMATICS, you must *first* use this procedure to
;       get the "pixel-based velocities" *before* using
;       MDAP_GET_BROADENED_TEMPLATE to generate the best-fitting stellar
;       continuum using the the correct LOSVD kernel and the optimal template.
;
; CALLING SEQUENCE:
;       MDAP_REVERT_PIXEL_KINEMATICS
;
; INPUTS:
;       cz double
;               The redshift obtained after converting the 'vel' returned by
;               PPXF/GANDALF to cz.
;
;       cz_err double
;               The error in 'vel'.
;
; OUTPUT:
;       The input cz and cz_err are converted to the pixel shifts vel and
;       vel_err upon output.  Note that the error is NOT a propagation of the
;       error in the calculation used to calculate vel.  It is the inverse of
;       the operation used to convert the 'vel' error to the error in cz applied
;       by MDAP_CONVERT_PIXEL_KINEMATICS.
;
; REVISION HISTORY:
;       05 Nov 2014: (KBW) Original implementation
;-
;------------------------------------------------------------------------------

PRO MDAP_REVERT_PIXEL_KINEMATICS, $
        cz, cz_err

        c=299792.458d                           ; Speed of light in km/s
;       cz_err = abs( cz_err / (cz/c+1.0d) )    ; error in vel (using error propagation)
        cz = c*alog(cz/c+1.0d)                  ; vel
        cz_err = cz_err/exp(cz/c)               ; NOT an error propagation!  (inverse of convert)
END


