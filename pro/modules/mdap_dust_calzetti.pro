function mdap_dust_calzetti,l0_gal,lstep_gal,npix,ebv,vstar,LOG10=log10
; This procedure uses the dust model of Calzetti et al. (2000, ApJ,
; 533, 682), and for a given E(B-V) value returns the flux attenuation
; array, which can be used to get reddened templates. Here the spectra
; are assumed to be binned on a ln-rebinned wavelentgh grid as defined
; by input l0_gal,lstep_gal,npix parameters. The input receiding
; velocity vstar, is used to derive the dust reddening in the galaxy
; rest-frame.
; 
; Can be used also to de-reddened the galaxy spectra by the Milky-Way
; dust extinction, using as E(B-V) the opposite of the Schlegel et
; al. values found in NED and vstar = 0.
;
; Initial version kindly provided by S. Kaviray, Oxford, 2006.

; reconstruct the wavelength array in Anstroms, and compute rest-frame
; values
lambda = exp(dindgen(npix)*lstep_gal + l0_gal)
if keyword_set(log10) then lambda = 10^(dindgen(npix)*lstep_gal + l0_gal)
lambda = lambda/exp(vstar/299792.458d)

; array to hold k(lambda) values
k = fltarr(n_elements(lambda))           

for i=0,n_elements(lambda)-1 do begin
     ; convert wavelength units from angstroms to micrometres
     l = lambda(i)/1e4                   
     ; assign k values
     if (l ge 0.63 and l le 2.2) then k(i) = 2.659*(-1.857+1.040/l)+4.05
     if (l lt 0.63)              then k(i) = 2.659*(-2.156+1.509/l-0.198/l^2+0.011/l^3)+4.05
     if (l gt 2.2)               then k(i) = 0.0
endfor

return,(10^(-0.4*ebv*k))
; this should be then multiplied by the spectrum flux array
end
