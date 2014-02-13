; NAME:
;
; mdap_convol_sigma
;
;
; AUTHOR: 
;
; L. Coccato
;
; PURPOSE:
;
; Convolution of a vector with a Gaussian function with variable sigma
;
;
; USAGE:
;
; Return = mdap_convol_sigma(x,vector,x_sigma,sigma)
;
;
; RETURN VALUE
;
; Convolved vector, same elements of "vector"
;
; INPUT:
;
; x         [1D array]  Values for which "vector" is defined 
;                       Units: pixels (or angstrom, but with constant
;                       ang/step sampling).
;
; vector    [1D array]  vector to be convolved. Same elements of "x"
; 
; x_sigma   [1D array]  Values for which "sigma" is defined
;
; sigma     [1D array]  Values of the dispersion of the convolving
;                       Gaussian function. Same elements of "x_sigma"
;                       Units: same as X
;
;
;
; HISTORY:
;
; v1.0   30/07/2009 by L. Coccato original code
; v2.0   19/08/2009 by L. Coccato bug on kernel area fixed. Boundaries
;                         kernel area where overestimated with total(kernel). 
;                         Analitic expression used instead. Double
;                         precision used.
; v2.1   13/11/2009 by L. Coccato.  minor changes in the text.                    
; INSERT NEW VERSION NUMBER
;
;%%%

function mdap_kroneker_product,array,element

  fun=fltarr(n_elements(array));*0.
  fun[element]=array[element]
return,fun
end
function mdap_convol_sigma,x,vector,x_sigma,sigma



sigma_new=interpol(sigma,x_sigma,x)
;print, 'N=',n_elements(vector)
conv=fltarr(n_elements(x))

;sigma_new2=2.d*sigma_new^2
sigma_new2=2.*sigma_new^2
;sigma_new_sqrt2pi = sigma_new*sqrt(2.d*double(!pi))
sigma_new_sqrt2pi = sigma_new*sqrt(2.*!pi)
med_step = (max(x)-min(x))/n_elements(x)
min_step = 1.*med_step
for i = 0, n_elements(vector)-2 do begin
   area = (x[i+1]-x[i])*vector[i]
   if sigma_new(i) gt min_step then kernel = exp(-(x-x(i))^2/sigma_new2(i))*area/sigma_new_sqrt2pi[i] else kernel = mdap_kroneker_product(vector,i)
   conv=conv+temporary(kernel)
endfor
area = (x[i]-x[i-1])*vector[i]
if sigma_new(i) gt min_step then kernel = exp(-(x-x(i))^2/sigma_new2(i))*area/sigma_new_sqrt2pi[i]  else kernel = mdap_kroneker_product(vector,i)
conv=conv+temporary(kernel)


ind = where(sigma_new le 1.05*min_step)
if ind[0] ne -1 then conv[ind]=vector[ind]


return,conv
end



