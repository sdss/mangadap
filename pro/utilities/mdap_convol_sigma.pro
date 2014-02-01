; NAME:
;
; convol_sigma
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
; Return = convol_sigma(x,vector,x_sigma,sigma,[/plot_example],[/help])
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
; OPTIONAL KEYWORDS:
; 
;
;                 on "vector" are done. Return value = -1
;
; /help           If set, the help screen is printed. Return value = -1
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
;
;%%%


function mdap_convol_sigma,x,vector,x_sigma,sigma;,help=help;,slow=slow

;if keyword_set(help) then begin
;   Result = FILE_WHICH('convol_sigma.pro')
;   readcol,Result,h,format='A,A',delimiter='@',/silent
;   end_of_help = where(h eq ';%%%' )
;   for i = 0, end_of_help(0) do print, h(i)
 ;  conv=-1
 ;  goto, fine
;endif

sigma_new=interpol(sigma,x_sigma,x)

conv=0.d

sigma_new2=2.d*sigma_new^2
sigma_new_sqrt2pi = sigma_new*sqrt(2.d*double(!pi))
for i = 0, n_elements(vector)-2 do begin
   area = (x[i+1]-x[i])*vector[i]
   kernel=exp(-(x-x(i))^2/sigma_new2(i))*area/sigma_new_sqrt2pi[i] 
   conv=conv+kernel
endfor
area = (x[i]-x[i-1])*vector[i]
kernel=exp(-(x-x(i))^2/sigma_new2(i))*area/sigma_new_sqrt2pi[i] 
conv+=kernel
;stop
;goto,fine

;fine:
return,conv
end



