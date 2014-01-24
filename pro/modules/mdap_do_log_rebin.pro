;######################################################################
;
; Copyright (C) 2001-2011, Michele Cappellari
; E-mail: cappellari_at_astro.ox.ac.uk
;
; This software is provided as is without any warranty whatsoever.
; Permission to use, for non-commercial purposes is granted.
; Permission to modify for personal or internal use is granted,
; provided this copyright and disclaimer are included unchanged
; at the beginning of the file. All other rights are reserved.
;
;######################################################################

; NAME:
;   MDAP_DO_LOG_REBIN
;
; PURPOSE:
;   Logarithmically rebin a spectrum, while rigorously conserving the flux. 
;   Basically the photons in the spectrum are simply ridistributed according 
;   to a new grid of pixels, with non-uniform size in the spectral direction.
;
;   This routine makes the `standard' zero-order assumption that the spectrum
;   is *constant* within each pixels. It is possible to perform log-rebinning
;   by assuming the spectrum is represented by a piece-wise polynomial of
;   higer degree, while still obtaining a uniquely defined linear problem,
;   but this reduces to a deconvolution and amplifies noise.
;
;   This same routine can be used to compute approximate errors
;   of the log-rebinned spectrum. To do this type the command
;
;       MDAP_DO_LOG_REBIN, lamRange, err^2, err2New
;
;   and the desired errors will be given by SQRT(err2New).
;   NB: This rebinning of the error-spectrum is very *approximate* as 
;   it does not consider the correlation introduced by the rebinning!
;
; CALLING SEQUENCE:
;   MDAP_DO_LOG_REBIN, lamRange, spec, specNew, logLam, $
;       OVERSAMPLE=oversample, VELSCALE=velScale, /FLUX
;
; INPUTS:
;   LAMRANGE: two elements vector containing the central wavelength
;       of the first and last pixels in the spectrum, which is assumed
;       to have constant wavelength scale! E.g. from the values in the
;       standard FITS keywords: LAMRANGE = CRVAL1 + [0,CDELT1*(NAXIS1-1)].
;       It must be LAMRANGE[0] < LAMRANGE[1].
;   SPEC: input spectrum.
;
; OUTPUTS:
;   SPECNEW: logarithmically rebinned spectrum.
;   LOGLAM: log(lambda) (*natural* logarithm: ALOG) of the central 
;       wavelength of each pixel. This is the log of the geometric 
;       mean of the borders of each pixel.
;
; KEYWORDS:
;   FLUX: Set this keyword to preserve total flux. In this case the 
;       log rebinning changes the pixels flux in proportion to their 
;       dLam so the following command will show large differences 
;       beween the spectral shape before and after LOG_REBIN:
;       
;           plot, exp(logLam), specNew  ; Plot log-rebinned spectrum
;           oplot, range(lamRange[0],lamRange[1],n_elements(spec)), spec
;           
;       By defaul, when this keyword is *not* set, the above two lines 
;       produce two spectra that almost perfectly overlap each other.
;   OVERSAMPLE: Oversampling can be done, not to loose spectral resolution, 
;       especally for extended wavelength ranges and to avoid aliasing.
;       Default: OVERSAMPLE=1 ==> Same number of output pixels as input.
;   VELSCALE: velocity scale in km/s per pixels. If this variable is
;       not defined, then it will contain in output the velocity scale.
;       If this variable is defined by the user it will be used
;       to set the output number of pixels and wavelength scale.
;
; MODIFICATION HISTORY:
;   V1.0: Using interpolation. Michele Cappellari, Leiden, 22 October 2001
;   V2.0: Analytic flux conservation. MC, Potsdam, 15 June 2003
;   V2.1: Allow a velocity scale to be specified by the user.
;       MC, Leiden, 2 August 2003
;   V2.2: Output the optional logarithmically spaced wavelength at the
;       geometric mean of the wavelength at the border of each pixel.
;       Thanks to Jesus Falcon-Barroso. MC, Leiden, 5 November 2003
;   V2.21: Verify that lamRange[0] < lamRange[1]. 
;       MC, Vicenza, 29 December 2004
;   V2.22: Modified the documentation after feedback from James Price.
;       MC, Oxford, 21 October 2010
;   V2.3: By default now preserve the shape of the spectrum, not the 
;       total flux. This seems what most users expect from the procedure.
;       Set the keyword /FLUX to preserve flux like in previous version.
;       MC, Oxford, 30 November 2011
;----------------------------------------------------------------------


;----------------------------------------------------------------------
pro mdap_do_log_rebin, lamRange, spec, specNew, logLam, $
    OVERSAMPLE=oversample, VELSCALE=velScale, FLUX=flux


compile_opt idl2
on_error, 2

if n_elements(lamRange) ne 2 then message, 'lamRange must contain two elements'
if lamRange[0] ge lamRange[1] then message, 'It must be lamRange[0] < lamRange[1]'
s = size(spec)
if s[0] ne 1 then message, 'input spectrum must be a vector'
n = s[1]
if n_elements(oversample) ne 0 then m = n*oversample else m = n

dLam = (lamRange[1]-lamRange[0])/(n-1d)          ; Assume constant dLam
lim = lamRange/dLam + [-0.5d,0.5d]               ; All in units of dLam
borders = range(lim[0],lim[1],n+1)               ; Linearly
logLim = alog(lim)

c = 299792.458d                                  ; Speed of light in km/s
if n_elements(velScale) gt 0 then begin          ; Velocity scale is set by user
    logScale = velScale/c
    m = floor((logLim[1]-logLim[0])/logScale)    ; Number of output pixels
    logLim[1] = logLim[0] + m*logScale
endif else $
    velScale = (logLim[1]-logLim[0])/m*c         ; Only for output

newBorders = exp(range(logLim[0],logLim[1],m+1)) ; Logarithmically
k = floor(newBorders-lim[0]) < (n-1) > 0
if n_params() ge 4 then $  ; Optional output log(wavelength): log of geometric mean!
    logLam = alog(sqrt(newBorders[1:*]*newBorders)*dLam)
 
specNew = dblarr(m,/NOZERO)
for j=0,m-1 do begin                             ; Do analytic integration
    a = newBorders[j] - borders[k[j]]
    b = borders[k[j+1]+1] - newBorders[j+1]
    specNew[j] = total(spec[k[j]:k[j+1]]) - a*spec[k[j]] - b*spec[k[j+1]]
endfor

if ~keyword_set(flux) then specNew /= newBorders[1:*] - newBorders

end
;----------------------------------------------------------------------

