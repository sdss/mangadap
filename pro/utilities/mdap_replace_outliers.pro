;+
; NAME:
;       MDAP_REPLACE_OUTLIERS
;
; PURPOSE:
;       Replace outlying data, as determined by the model data, the error
;       vector, and the number of sigma, with the values of a model (or just
;       another vector).
;
; CALLING SEQUENCE:
;       MDAP_REPLACE_OUTLIERS, data, model, sigma, nsigma
;
; INPUTS:
;       data dblarr[D]
;               Data vector.  Possibly changed on output with some pixels
;               replaced with model values.
;
;       model dblarr[D]
;               Values to determine the residuals of the data and used to
;               replace any outlying data.
;
;       sigma dblarr[D]
;               Standard deviation of each datum about its value in data.
;
;       nsigma double
;               Number of sigma for allowed deviations
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;
; OUTPUT:
;
; OPTIONAL OUTPUT:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; TODO:
;       - Allow the model vector and values used to replace data to be different
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;
;-
;------------------------------------------------------------------------------

PRO MDAP_REPLACE_OUTLIERS, $
                data, model, sigma, nsigma

        chi = abs((data-model)/sigma)           ; Absolute deviation in units of sigma
        indx = where(chi gt nsigma)             ; Indices of outlying data
        if indx[0] = -1 then $                  ; No outliers, so just return
            return
        data[indx] = model[indx]                ; Replace data with model values for outliers
END


