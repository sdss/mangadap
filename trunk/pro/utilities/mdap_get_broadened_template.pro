;+
; NAME:
;       MDAP_GET_BROADENED_TEMPLATE
;
; PURPOSE:
;       Broaden (and shift) by the line-of-sight velocity distribution (LOSVD)
;       using the parameters in the input sol vector.  This code is pulled
;       directly from pPXF.
;
; CALLING SEQUENCE:
;       result = MDAP_GET_BROADENED_TEMPLATE(template, sol, velscale, moments=moments, $
;                                            vsyst=vsyst, /oversample)
;
; INPUTS:
;       template dblarr[]
;               Template spectrum to be convolved.
;
;       sol dblarr[]
;               Vector with the parameters of the LOSVD.  There must be at least
;               4 entries:
;                       0: Velocity (km/s)
;                       1: Sigma (km/s)
;                       2: h3
;                       3: h4
;               If present, h5 and h6 are used as well (elements 4 and 5)
;
;       velscale double
;               Pixel scale in km/s of input template spectrum.
;
; OPTIONAL INPUTS:
;       moments integer
;               Sets the number of moments to use in creating the LOSVD.
;               Default is 4.
;
;       vsyst double
;               Reference velocity for the measurements in sol in km/s.  Default
;               is 0.0d.
;
; OPTIONAL KEYWORDS:
;       /OVERSAMPLE
;               Set this keyword to oversample the template by a factor 30x
;               before convolving it with a well sampled LOSVD. This can be
;               useful to extract proper velocities, even when sigma <
;               0.7*velScale and the dispersion information becomes totally
;               unreliable due to undersampling.  IMPORTANT: One should sample
;               the spectrum more finely is possible, before resorting to the
;               use of this keyword! 
;
; OUTPUT:
;       result dblarr[]
;               The template broadened by the LOSVD defined by the parameters in
;               sol.
;
; OPTIONAL OUTPUT:
;
; COMMENTS:
;
; EXAMPLES:
;
; TODO:
;       - May be different from best fit returned by pPXF because this
;         does not include the multiplicative (or additive) polynomials.
;
; BUGS:
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;        8 Oct 2014: Original Implementation by K. Westfall (KBW)
;       14 Aug 2015: (KBW) Convolution performed over pixels with given
;                          velocity width, so needed to change the
;                          normalization accordingly.  Now return
;                          template * velscale, not just template.
;-
;------------------------------------------------------------------------------

FUNCTION MDAP_GET_BROADENED_TEMPLATE, $
                template, sol, velscale, moments=moments, vsyst=vsyst, $
                oversample=oversample

        npars = n_elements(moments) ne 0 ? moments : 4          ; Number of LOSVD parameters
        sz=size(sol)
        if sz[1] lt npars then $
            message, 'Incorrect number of kinematic parameters!'

        sz=size(template)
        if sz[0] ne 1 then $
            message, 'Only supports 1D templates'

        if keyword_set(oversample) then factor = 30 else factor = 1

        voff = n_elements(vsyst) ne 0 ? vsyst : 0.0d            ; Offset velocity in km/s

        ; sol = [vel,sigma,h3,h4,h5,h6,...]                     ; Velocities are in km/s
        dx = ceil((abs(voff)+abs(sol[0])+5d*sol[1])/velscale)   ; Sample LOSVD to vsyst+vel+5*sigma
        n = 2*dx*factor + 1                                     ; Number of pixels required
        v = mdap_range(dx,-dx,n)*velscale                       ; Velocity range for LOSVD samples

        losvd = dblarr(n,/NOZERO)                               ; LOSVD array
        
        w = (v - (voff+sol[0]))/sol[1]
        w2 = w^2
        losvd = exp(-0.5d*w2)/(sqrt(2d*!dpi)*sol[1])            ; Normalized total(Gaussian)=1

        ; Hermite polynomials normalized as in Appendix A of van der Marel & Franx (1993).
        ; These coefficients are given e.g. in Appendix C of Cappellari et al. (2002)
        ;
        if npars gt 2 then begin
            poly = 1d + sol[2]/sqrt(3d)*(w*(2d*w2-3d)) $                        ; H3
                      + sol[3]/sqrt(24d)*(w2*(4d*w2-12d)+3d)                    ; H4

            if npars eq 6 then $
                poly = poly + sol[4]/sqrt(60d)*(w*(w2*(4d*w2-20d)+15d)) $       ; H5
                            + sol[5]/sqrt(720d)*(w2*(w2*(8d*w2-60d)+90d)-15d)   ; H6

            losvd = losvd*poly
        endif

        ; Convolve the spectrum, oversampling if requested
        if factor gt 1 then pix = mdap_range(0d,sz[1]-1d,sz[1]*factor) ; Oversampled pixels range
        broadened_template = dblarr(sz[1],/NOZERO)
        if factor eq 1 then begin                       ; No oversampling of the template spectrum
            broadened_template = mdap_ppxf_convol_fft(template,losvd)
        endif else begin                                ; Oversample the template then convolve
            tpl = interpolate(template,pix,CUBIC=-0.5)  ; Sinc-like interpolation
            ; Convolve oversampled template then rebin to orig size
            broadened_template = rebin(mdap_ppxf_convol_fft(tpl,losvd),sz[1])
        endelse

        return, broadened_template*velscale
end

