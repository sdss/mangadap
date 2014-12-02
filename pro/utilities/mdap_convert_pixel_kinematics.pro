;+
; NAME:
;       MDAP_CONVERT_PIXEL_KINEMATICS
;
; PURPOSE:
;       Convert kinematics from PPXF/GANDALF from pixel shifts to redshifts.
;       PPXF/GANDALF determines the velocity offset by making the approximation
;       that every pixel (logarithmically binned in wavelength) is a constant
;       change in velocity.  Over large velocity shifts, this approximation can
;       become poor.  An e-mail summary from Michele Cappellari:
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
;       This function converts the 'vel' values provided by PPXF/GANDALF to cz
;       velocities.
;
;       **IMPORTANT**: After this conversion, you must *first* use
;       MDAP_REVERT_PIXEL_KINEMATICS to get the "pixel-based velocities"
;       *before* using MDAP_GET_BROADENED_TEMPLATE to generate the best-fitting
;       stellar continuum using the the correct LOSVD kernel and the optimal
;       template.
;
; CALLING SEQUENCE:
;       MDAP_CONVERT_PIXEL_KINEMATICS, v, v_err
;
; INPUTS:
;       v double
;               The 'vel' returned by PPXF/GANDALF, which is the pixel shift
;               times the pixel width.
;
;       v_err double
;               The error in 'vel'
;
; OUTPUT:
;       The input vel and vel_err are converted to redshift (cz) and error upon
;       output.
;
; REVISION HISTORY:
;       05 Nov 2014: (KBW) Original implementation, moved from
;                          MDAP_SPECTRAL_FITTING
;-
;------------------------------------------------------------------------------

PRO MDAP_CONVERT_PIXEL_KINEMATICS, $
        v, v_err

        c=299792.458d                   ; Speed of light in km/s
        v_err = abs(exp(v/c)*v_err)     ; error in cz (done here because v is replaced next)
        v = (exp(v/c)-1.0d)*c           ; cz
END


