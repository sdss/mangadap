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


function mdap_convol_sigma,x,vector,x_sigma,sigma,help=help;,slow=slow

if keyword_set(help) then begin
   Result = FILE_WHICH('convol_sigma.pro')
   readcol,Result,h,format='A,A',delimiter='@',/silent
   end_of_help = where(h eq ';%%%' )
   for i = 0, end_of_help(0) do print, h(i)
   conv=-1
   goto, fine
endif

if keyword_set(plot_example) then begin
   exec_example;,slow=slow
   conv=-1
   goto, fine
end
sigma_new=interpol(sigma,x_sigma,x)

step=double(x(1)-x(0))
conv=0.d

;if keyword_set(slow) then goto,do_slow

area = step*vector
for i = 0, n_elements(vector)-1 do begin
   kernel= exp(-(x-x(i))^2./(2.*sigma_new(i)^2.))
  ; kernel=kernel/total(kernel)*(area(i))      ;v1.0  wrong!!!
   kernel=kernel/(sigma_new(i)*sqrt(2.d*double(!pi)))*(area(i))  ;v2.0 correct!!!
   conv=conv+kernel
endfor
;stop
goto,fine
;
;
;this should be faster. NO: IT IS NOT!!!!
;do_slow:
; area = step*vector
; area_2d = (1.+area*0.) ## area
; x_2d = (1+x*0.) ## x  
; xi = x ## (1+x*0.) 
; sigma_2d= (1.+sigma_new*0.) ## sigma_new
; kernel_2d = area_2d/(sigma_2d*sqrt(2.*!pi)) * exp(-(x_2d-xi)^2./(2.*sigma_2d^2.))
; conv = total(kernel_2d,1)
;

fine:
return,conv
end

