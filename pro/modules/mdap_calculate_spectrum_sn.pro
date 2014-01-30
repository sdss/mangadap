pro mdap_calculate_spectrum_sn,spectrum,error,wavelength,sn_per_angstorm,$
                               signal=signal,noise=noise,double=double,rms=rms

; This procedure computes the mean S/N per angstrom of a given
; spectrum. 

; INPUTS
;spectrum   [N elements array] the spectrum to compute the SN. Units = Flux / angstrom
;error      [N elements array] the vector array associated to spectrum
;wavelength [N elements array] the vector over which spectrum and
;                              error are defined. Units = angstrom 
; KEYWORDS
; /double   If set, the computation is done in double precision.
; /rms      If set, the error input spectrum is assumed to be the
;           residual from a best fit model, the input spectrum is
;           assumed to be the best fit model, and the SN is computed 
;           as:
;                           median(spectrum,/even)            max(wavelength)- min(wavelength)
;                 S/N = --------------------------- * sqrt( ---------------------------------- )
;                           robust_sigma(error)                    n_elements(wavelength)
;
;
;
; OUTPUTS
;
; sn_per_angstorm [double/float]] Mean S/N per angstrom (computed over
;                                 the entire wavelength range).

; OPTIONAL OUTPUTS
;
; signal [double/float]  Mean signal per angstrom (computed over
;                        the entire wavelength range).
;
;                            signal = median(spectrum)
;
;                        unless the keyword /rms is set.
;
;
; noise  [double/float]  Mean noise per angstrom (computed over
;                        the entire wavelength range).
;
;                            noise = median(error)
;
;                        unless the keyword /rms is set.
;
;
;   <S/N> = median(spectrum/error)
; 
; 9 Jan 2014 by L. Coccato

lrange=  (max(wavelength)-min(wavelength))
if ~keyword_set(rms) then begin

  ; signal = int_tabulated(wavelength,spectrum,double=double) / lrange
  ; noise  =  sqrt(int_tabulated(wavelength,error^2,double=double) / lrange)

   signal = median(spectrum,/even)
   noise = median(error,/even)
   sn_per_angstorm = median(spectrum/error,/even)   ;equivalent of median( flux * sqrt(ivar) ,/even)
endif else begin
   signal = median(spectrum,/even)
   noise = robust_sigma(error)    
   sn_per_angstorm = signal / noise * sqrt(lrange/n_elements(wavelength))
endelse

;  <S/N> = sum in quadrature of the signal to noises, divided by the  number of elements. (Matt's suggestion)
; dispersion = wavelength[1:*] - wavelength[0:n_elements(wavelength)-2]
; dispersion = [dispersion,dispersion[n_elements(dispersion)-1]]
; sn_per_angstorm = sqrt(total(spectrum/error*sqrt(dispersion))^2)
; signal = total(spectrum)*sqrt(dispersion)/n_elements(wavelength)
; noise =sqrt(total(error^2))*sqrt(dispersion)/n_elements(wavelength)

; <S/N> from Cardiel et al. (2008) formula (Equation 42)
; dispersion = wavelength[1:*] - wavelength[0:n_elements(wavelength)-2]
; dispersion = [dispersion,dispersion[n_elements(dispersion)-1]]
; sn_per_angstorm = 1./n_elements(wavelength) * int_tabulated(wavelength,spectrum/error/sqrt(dispersion))
; signal = 1./n_elements(wavelength) *int_tabulated(wavelength,spectrum/sqrt(dispersion))
; noise = 1./n_elements(wavelength) *int_tabulated(wavelength,error/sqrt(dispersion))
end
