;+
; NAME:
;       MDAP_MINIMUM_CNVLV_SIGMA
;
; PURPOSE:
;       Hard-wired procedure used to set the minimum allowed sigma to use in
;       application of a convolution.  This is different from
;       MDAP_MIN_SIG_GAU_APPROX_DELTA!  This function should be used to set an
;       absolute minimum sigma that is sent to MDAP_CONVOL_SIGMA, whereas
;       MDAP_MIN_SIG_GAU_APPROX_DELTA tells MDAP_CONVOL_SIGMA when it should
;       approximate the Gaussian kernel by a Kronecker delta function.  For
;       example, if you never want the kernel to be approximated, this function
;       should return the same value as MDAP_MIN_SIG_GAU_APPROX_DELTA or higher.
;       If you are okay with the approximation being used, this function should
;       return a value that is lower than that returned by
;       MDAP_MIN_SIG_GAU_APPROX_DELTA.
;
; CALLING SEQUENCE:
;       result=MDAP_MINIMUM_CNVLV_SIGMA(dx)
;
; INPUTS:
;       dx double
;               Discrete sample size
;
; OUTPUT:
;       returned double
;               Returned value is the minimum allowed sigma.
;
; REVISION HISTORY:
;       17 Sep 2014: (KBW) Original implementation
;-
;------------------------------------------------------------------------------

FUNCTION MDAP_MINIMUM_CNVLV_SIGMA, $
                dx
;       return, (MDAP_MIN_SIG_GAU_APPROX_DELTA(dx))
        return, 0.0d
END

