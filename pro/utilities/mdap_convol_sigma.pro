;+
; NAME:
;       MDAP_CONVOL_SIGMA
;
; PURPOSE:
;       Convolution of a vector with a Gaussian function with variable sigma.
;
; CALLING SEQUENCE:
;       Return = MDAP_CONVOL_SIGMA(x, vector, x_sigma, sigma)
;
; INPUTS:
;       x dblarr[]
;               Values for which "vector" is defined. Units: pixels (or
;               angstrom, but with constant ang/step sampling).
;
;       vector dblarr[]
;               Value to be convolved at each x.
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
; OPTIONAL OUTPUT:
;
; COMMENTS:
;
; EXAMPLES:
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
;-
;------------------------------------------------------------------------------

;-------------------------------------------------------------------------------
; Return the integral over a Kronecker delta function
FUNCTION MDAP_KRONEKER_PRODUCT, $
                array, element

        fun=dblarr(n_elements(array))                   ; Set all values to zero
        fun[element]=array[element]                     ; Except at location of delta function

        return, fun
END
;-------------------------------------------------------------------------------


FUNCTION MDAP_CONVOL_SIGMA, $
                x, vector, x_sigma, sigma

        sigma_new = interpol(sigma, x_sigma, x)         ; Interpolate sigma to the x coo

        nx = n_elements(x)                              ; Number of vector elements

        conv = dblarr(nx)                               ; Initialize the array to zero

        gau_denom = 2.*sigma_new^2                      ; Denominator in exp() of Gaussian
        gau_norm = sigma_new*sqrt(2.*!pi)               ; Normalization of Gaussian

        dx = (x[nx-1]-x[0])/double(nx-1)                ; Sampling step
;       print, dx, x[1]-x[0]
        min_sig = MDAP_MIN_SIG_GAU_APPROX_DELTA(dx)     ; Minimum sigma allowed before using
                                                        ;     Kronecker delta approximation
;       print, 'min_sig: ', min_sig

        ; Convolution is:
        ;       f*g(x) = integral_-inf^inf f(X) g(x-X) dX
        ;  below is the Riemann sum approximation

;       plot, x, vector
;       stop
        for i=0L,nx-1 do begin

            if sigma_new[i] gt min_sig then begin
                integrand = exp( -(x-x[i])^2 / gau_denom[i])/gau_norm[i] * vector[i] * dx

                ; This operation nearly *always* produces underflow.  At the end
                ; of the loop, clear the math error.

;               ; Check for underflow
;               if check_math() eq 32 then $
;                   integrand[*] = 0.0d
;               if check_math() eq 32 then begin
;                   print, 'UNDERFLOW!'
;                   integrand[*] = 0.0d
;               endif
            endif else $
                integrand = MDAP_KRONEKER_PRODUCT(vector,i)

            conv=conv+integrand

        endfor

        if check_math(/noclear) eq 32 then $
            junk=check_math()                           ; Clear underflow errors
        
        ; Force regions where the convolution sigma is too close to the minimum
        ; sigma to have exactly the input value
        ind = where(sigma_new le 1.05*min_sig, count)
;       if ind[0] ne -1 then conv[ind]=vector[ind]
        if count ne 0 then conv[ind]=vector[ind]

;       oplot, x, vector, color=200
;       stop

;       plot, x, conv-vector
;       stop

        RETURN, conv
END


