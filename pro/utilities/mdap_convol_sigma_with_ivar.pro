;+
; NAME:
;       MDAP_CONVOL_SIGMA_WITH_IVAR
;
; PURPOSE:
;       Convolution of a vector with a Gaussian function with variable sigma.
;
; CALLING SEQUENCE:
;       MDAP_CONVOL_SIGMA_WITH_IVAR, x, vector, ivar, x_sigma, sigma
;
; INPUTS:
;       x dblarr[]
;               Values for which "vector" is defined. Units: pixels (or
;               angstrom, but with constant ang/step sampling).
;
;       vector dblarr[]
;               Value to be convolved at each x.
;
;       ivar dblarr[]
;               Inverse variances in vector
;
;       x_sigma dblarr[]
;               Values for which "sigma" is defined
;
;       sigma dblarr[]
;               Values of the dispersion of the convolving Gaussian function at
;               each x_sigma.  Units: same as X
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;
; OUTPUT:
;       Convolved vector with the same number of elements as vector.
;
;       conv dblarr[]
;               Result of convolution.
;
;       conv_ivar dblarr[]
;               Approximation of the inverse variance of conv.
;
; OPTIONAL OUTPUT:
;
; COMMENTS:
;
; EXAMPLES:
;
; TODO:
;       - Go back and make sure that I don't divide by zero!
;
; BUGS:
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
; v1.0   30/07/2009 by L. Coccato original code
; v2.0   19/08/2009 by L. Coccato bug on kernel area fixed. Boundaries
;                         kernel area where overestimated with total(kernel). 
;                         Analitic expression used instead. Double
;                         precision used.
; v2.1   13/11/2009 by L. Coccato.  minor changes in the text.                    
; v3.0   17 Sep 2014: (KBW) Formatting, editing, comments
; v3.1   22 Sep 2014: (KBW) Include rough (probably incorrect) approximation of
;                           errors
;-
;------------------------------------------------------------------------------

;-------------------------------------------------------------------------------
; Return the integral over a Kronecker delta function
FUNCTION MDAP_KRONEKER_PRODUCT, $
                array, element

        fun=fltarr(n_elements(array))                   ; Set all values to zero
        fun[element]=array[element]                     ; Except at location of delta function

        return, fun
END
;-------------------------------------------------------------------------------


PRO MDAP_CONVOL_SIGMA_WITH_IVAR, $
                x, vector, ivar, x_sigma, sigma, conv, conv_ivar

        sigma_new = interpol(sigma, x_sigma, x)         ; Interpolate sigma to the x coo

        nx = n_elements(x)                              ; Number of vector elements

        conv = dblarr(nx)                               ; Initialize the arrays to zero
        conv_ivar = dblarr(nx)

        gau_denom = 2.*sigma_new^2                      ; Denominator in exp() of Gaussian
        gau_norm = sigma_new*sqrt(2.*!pi)               ; Normalization of Gaussian

        dx = (x[nx-1]-x[0])/double(nx-1)                ; Sampling step
        min_sig = MDAP_MIN_SIG_GAU_APPROX_DELTA(dx)     ; Minimum sigma allowed before using
                                                        ;     Kronecker delta approximation

        divbyreal=reform(where(ivar gt 0, ndr, complement=divbyzero, ncomplement=ndz))
;       ndr = n_elements(divbyreal)

        ; Convolution is:
        ;       f*g(x) = integral_-inf^inf f(X) g(x-X) dX
        ;  below is the Riemann sum approximation
;       plot, x, vector
;       stop
        for i=0L,ndr-1 do begin

            if sigma_new[divbyreal[i]] gt min_sig then begin

                ; Convolution integrand
                integrand = double(exp( -(x-x[divbyreal[i]])^2 / gau_denom[divbyreal[i]]) $
                                   / gau_norm[divbyreal[i]] * vector[divbyreal[i]] * dx)
                if check_math(/noclear) eq 32 then $
                    junk=check_math()                   ; Check for underflow error and clear them

                ; Simple error propagation
                integrand_e = double((exp( -(x-x[divbyreal[i]])^2 / gau_denom[divbyreal[i]]) $
                                      / gau_norm[divbyreal[i]] * dx)^2 / ivar[divbyreal[i]])
                if check_math() eq 32 then $
                    junk=check_math()                   ; Check for underflow error and clear them

            endif else begin
                integrand = MDAP_KRONEKER_PRODUCT(vector,divbyreal[i])
                integrand_e = dblarr(nx)
                integrand_e[divbyreal] = MDAP_KRONEKER_PRODUCT(1.0/ivar[divbyreal],i)
;               print, nx, n_elements(integrand_e)
            endelse

            conv=conv+temporary(integrand)
            conv_ivar=conv_ivar+temporary(integrand_e)

        endfor
        ; Force regions where the convolution sigma is too close to the minimum
        ; sigma to have exactly the input value
        ind = where(sigma_new le 1.05*min_sig, count)
;       if ind[0] ne -1 then begin
        if count ne 0 then begin
            conv[ind]=vector[ind]
            conv_ivar[ind]=1.0/ivar[ind]
        endif

;       oplot, x, conv, color=200
;       stop
;       plot, x, conv-vector
;       stop

        conv_ivar[divbyreal] = 1.0/conv_ivar[divbyreal]
;       if divbyzero[0] ne -1 then $
        if ndz ne 0 then $
            conv_ivar[divbyzero] = 0.
;       stop

END


