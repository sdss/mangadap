;#############################################################################
;
; Copyright (C) 2001-2012, Michele Cappellari
; E-mail: cappellari_at_astro.ox.ac.uk
;
; Updated versions of the software are available from my web page
; http://purl.org/cappellari/idl
;
; If you have found this software useful for your research,
; I would appreciate an acknowledgment to the use of the
; "Penalized Pixel-Fitting method by Cappellari & Emsellem (2004)".
;
; This software is provided as is without any warranty whatsoever.
; Permission to use, for non-commercial purposes is granted.
; Permission to modify for personal or internal use is granted,
; provided this copyright and disclaimer are included unchanged
; at the beginning of the file. All other rights are reserved.
;
;#############################################################################
;+
; NAME:
;   MDAP_SGANDALF
;
; SYNTAX:
; mdap_sgandalf, templates1,templates2, loglam,galaxy, noise, velScale,
; start, sol,gas_intens,gas_fluxes,gas_ew, 
;     [BESTFIT=bestFit], [BIAS=bias], [\CLEAN], [DEGREE=degree], [ERROR=error], 
;     [GOODPIXELS=goodPixels], [MDEGREE=mdegree], [MOMENTS=moments],
;     [POLYWEIGHTS=polyweights], [\QUIET, [VSYST=vsyst], [WEIGHTS=weights],
;     [REGUL=regul], [LAMBDA=lambda], [REDDENING=reddening],[ADDITIVE_POL=additive_pol],
;     [BF_COMP1 = bf_comp1], [BF_COMP2 = bf_comp2], [T1=t1],[T2=t2],[mpoly=mpoly]
;
; PURPOSE:
;   Extract galaxy stellar and ionized gas kinematics (V, sigma, h3,
;   h4, h5, h6, Vgas, Sgas)
;   or the stellar population by fitting a template to an observed
;   spectrum in pixel space, using an implementation of the Penalized
;   Pixel-Fitting.
;   The original Penalized Pixel-Fitting method is described in: 
;   Cappellari M., & Emsellem E., 2004, PASP, 116, 138.
;
;   The following key optional features are also available:
;   1) An optimal template, positive linear combination of different
;      input templates, can be fitted together with the kinematics.
;   2) One can enforce smoothness on the template weights during the fit.
;      This is useful to attach a physical meaning to the weights
;      e.g. in terms of the star formation history of a galaxy.
;   3) Additive and/or multiplicative polynomials can be included to adjust
;      the continuum shape of the template to the observed spectrum.
;   4) Iterative sigma clipping can be used to clean the spectrum.
;   6) One can derive an estimate of the reddening in the spectrum.
;
; CALLING SEQUENCE:
; SGANDALF, templates1,templates2, loglam,galaxy, noise, velScale, start, sol,gas_intens,gas_fluxes,gas_ew, 
;     [BESTFIT=bestFit], [BIAS=bias], [\CLEAN], [DEGREE=degree], [ERROR=error], 
;     [GOODPIXELS=goodPixels], [MDEGREE=mdegree], [MOMENTS=moments],
;     [POLYWEIGHTS=polyweights], [\PLOT], [\QUIET, [VSYST=vsyst], [WEIGHTS=weights],
;     [REGUL=regul], [LAMBDA=lambda], [REDDENING=reddening],[ADDITIVE_POL=additive_pol],
;     [BF_COMP1 = bf_comp1], [BF_COMP2 = bf_comp2], [OPT_TEMPL = opt_templ],[mpoly=mpoly]
;
;
; INPUT PARAMETERS:
;   TEMPLATES1: vector containing the spectrum of a single template star or more
;       commonly an array of dimensions TEMPLATES[nPixels,nTemplates] containing
;       different templates to be optimized during the fit of the kinematics.
;       nPixels has to be >= the number of galaxy pixels.
;     - To apply linear regularization to the WEIGHTS via the keyword REGUL,
;       TEMPLATES should be an array of two TEMPLATES[nPixels,nAge], three
;       TEMPLATES[nPixels,nAge,nMetal] or four TEMPLATES[nPixels,nAge,nMetal,nAlpha]
;       dimensions, depending on the number of population variables one wants to study.
;       This can be useful to try to attach a physical meaning to the output WEIGHTS, in
;       term of the galaxy star formation history and chmemical composition distribution.
;       In that case the templates may represent single stellar population SSP models
;       and should be arranged in sequence of increasing age, metallicity or alpha along 
;       the second, third or fourth dimension of the array respectively.
;     - TEMPLATES and GALAXY do not need to span the same wavelength range. However
;       an error will be returned by SGANDALF, if the velocity shift in pixels,
;       required to match the galaxy with the templates, becomes larger than
;       nPixels. In that case one has to truncate either the galaxy or the
;       templates to make the two rest-frame spectral ranges more similar.
;   TEMPLATES2: [N elements X 2 elements] vector containing 
;       TEMPLATES2[*,0]  the values of the wavelenghts of the N emission lines to
;                        fit, in logarithmic units (ln or Log10 as in the input galaxy
;                        spectrum).
;       TEMPLATES2[*,1]  The N signs of the gas lines to fit, +1 for
;                        emission lines, -1 for absorption lines (e.g. NaI.) 
;   GALAXY: vector containing the spectrum of the galaxy to be measured. The
;       star and the galaxy spectra have to be logarithmically rebinned but the
;       continuum does *not* have to be subtracted. The rebinning may be
;       performed with the LOG_REBIN routine that is distributed with SGANDALF.
;     - For high redshift galaxies, one should bring the spectra close to the
;       restframe wavelength, before doing the SGANDALF fit, to prevent too large
;       velocity shifts of the templates. This can be done by dividing the
;       observed wavelenghts by (1+z), where z is a rough estimate of the
;       galaxy redshift, before the logarithmic rebinning.
;     - GALAXY can also be an array of dimensions GALAXY[nGalPixels,2] containing
;       two spectra to be fitted, at the same time, with a reflection-symmetric
;       LOSVD. This is useful for spectra taken at point-symmetric spatial
;       positions with respect to the center of an equilibrium stellar system.
;       For a discussion of the usefulness of this two-sided fitting
;       see e.g. Section 3.6 of Rix & White (1992, MNRAS, 254, 389).
;     - IMPORTANT: 1) For the two-sided fitting the VSYST keyword has to be used.
;       2) Make sure the spectra are rescaled to be not too many order of
;       magnitude different from unity, to avoid over or underflow problems
;       in the calculation. E.g. units of erg/(s cm^2 A) may cause problems!
;   NOISE: vector containing the 1*sigma error (per pixel) in the galaxy spectrum.
;       If GALAXY is a Nx2 array, NOISE has to be an array with the same dimensions.
;     - IMPORTANT: the penalty term of the sgandalf method is based on the *relative*
;       change of the fit residuals. For this reason the penalty will work as
;       expected even if no reliable estimate of the NOISE is available
;       (see Cappellari & Emsellem [2004] for details).
;       If no reliable noise is available this keyword can just be set to:
;           NOISE = galaxy*0+1 ; Same weight for all pixels
;   VELSCALE: velocity scale of the spectra in km/s per pixel. It has to be the
;       same for both the galaxy and the template spectra.
;   START: 4 elements vector: [velStart_stars,sigmaStart_stars,velStart_gas,sigmaStart_gas] 
;       with the initial estimate for the velocity and the velocity dispersion in km/s.
;     - Unless a good initial guess is available, it is recommended to set the starting
;       sigma >= 3*velScale in km/s (i.e. 3 pixels). In fact when the LOSVD is severely
;       undersampled, and far from the true solution, the chi^2 of the fit becomes weakly
;       sensitive to small variations in sigma (see sgandalf paper). In some instances the
;       near-constancy of chi^2 may cause premature convergence of the optimization.
;     - In the case of two-sided fitting a good starting value for the
;       velocity is velStart=0.0 (in this case VSYST will generally be nonzero).
;       Alternatively on should keep in mind that velStart refers to the first
;       input galaxy spectrum, while the second will have velocity -velStart.
;
;
;
; OPTIONAL INPUTS AND KEYWORDS:
;   BIAS: This parameter biases the (h3,h4,...) measurements towards zero
;       (Gaussian LOSVD) unless their inclusion significantly decreses the
;       error in the fit. Set this to BIAS=0.0 not to bias the fit: the
;       solution (including [V,sigma]) will be noisier in that case. The
;       default BIAS should provide acceptable results in most cases, but it
;       would be safe to test it with Monte Carlo simulations. This keyword
;       precisely corresponds to the parameter \lambda in the Cappellari &
;       Emsellem (2004) paper. Note that the penalty depends on the *relative*
;       change of the fit residuals, so it is insensitive to proper scaling
;       of the NOISE vector. A nonzero BIAS can be safely used even without a
;       reliable NOISE spectrum, or with equal weighting for all pixels.
;   /CLEAN: set this keyword to use the iterative sigma clipping method
;       described in Section 2.1 of Cappellari et al. (2002, ApJ, 578, 787).
;       This is useful to remove from the fit unmasked bad pixels, residual
;       gas emissions or cosmic rays.
;     - IMPORTANT: This is recommended *only* if a reliable estimate of the
;       NOISE spectrum is available. See also note below for SOL.
;   DEGREE: degree of the *additive* Legendre polynomial used to correct
;       the template continuum shape during the fit.
;       Default: DEGREE = -1, i.e. no additive polynomial are fitted.
;   GOODPIXELS: integer vector containing the indices of the good pixels in the
;       GALAXY spectrum (in increasing order). Only these pixels are included in
;       the fit. If the /CLEAN keyword is set, in output this vector will be updated
;       to contain the indices of the pixels that were actually used in the fit.
;     - IMPORTANT: in all likely situations this keyword *has* to be specified.
;   LAMBDA: When the keyword REDDENING is used, the user has to pass in this
;       keyword a vector with the same dimensions of GALAXY, giving the restframe
;       wavelength in Angstrom of every pixel in the input galaxy spectrum.
;       If one uses my LOG_REBIN routine to rebin the spectrum before the SGANDALF fit:
;           LOG_REBIN, lamRange, galaxy, galaxyNew, logLam
;       the wavelength can be obtained as lambda = EXP(logLam).
;   MDEGREE: degree of the *multiplicative* Legendre polynomial (with mean of 1)
;       used to correct the continuum shape during the fit (default: 0). The
;       zero degree multiplicative polynomial is always included in the fit as
;       it corresponds to the weights assigned to the templates.
;       Note that the computation time is longer with multiplicative
;       polynomials than with the same number of additive polynomials.
;     - IMPORTANT: Multiplicative polynomials cannot be used when
;       the REDDENING keyword is set.
;   MOMENTS: Order of the Gauss-Hermite moments to fit. Set this keyword to 4
;       to fit [h3,h4] and to 6 to fit [h3,h4,h5,h6]. Note that in all cases
;       the G-H moments are fitted (nonlinearly) *together* with [V,sigma].
;     - If MOMENTS=2 or MOMENTS is not set then only [V,sigma] are
;       fitted and the other parameters are returned as zero.
;     - If MOMENTS=0 then only the templates and the continuum additive
;       polynomials are fitted and the WEIGHTS are returned in output.
;   POLYWEIGHTS: vector with the weights of the additive Legendre polynomials.
;       The best fitting additive polynomial can be explicitly evaluated as
;           x = mdap_range(-1d,1d,n_elements(galaxy))
;           apoly = 0d ; Additive polynomial
;           for j=0,DEGREE do apoly += legendre(x,j)*polyWeights[j]
;     - When doing a two-sided fitting (see help for GALAXY parameter), the additive
;       polynomials are allowed to be different for the left and right spectrum.
;       In that case the output weights of the additive polynomials alternate between
;       the first (left) spectrum and the second (right) spectrum.
;   /QUIET: set this keyword to suppress verbose output of the best fitting
;       parameters at the end of the fit.
;   REDDENING: Set this keyword to an initail estimate of the reddening E(B-V)>=0
;       to fit a positive reddening together with the kinematics and the templates.
;       After the fit the input estimate is replaced with the best fitting E(B-V) value.
;       The fit assumes the exctinction curve of Calzetti et al. (2000, ApJ, 533, 682)
;       but any other prescriptions could be trivially implemented by modifying the
;       function SGANDALF_REDDENING_CURVE below.
;     - IMPORTANT: The MDEGREE keyword cannot be used when REDDENING is set.
;   REGUL: If this keyword is nonzero, the program applies second-degree
;       linear regularization to the WEIGHTS during the SGANDALF fit.
;       Regularization is done in one, two or three dimensions depending on whether
;       the array of TEMPLATES has two, three or four dimensions respectively.
;       Large REGUL values correspond to smoother WEIGHTS output. The WEIGHTS tend
;       to a linear trend for large REGUL. When this keyword is nonzero the solution
;       will be a trade-off between smoothness of WEIGHTS and goodness of fit.
;     - The effect of the regularization scheme is to enforce the numerical second 
;       derivatives between neighbouring weights (in every dimension) to be equal 
;       to -w[j-1]+2*w[j]-w[j+1]=0 with an error Delta=1/REGUL. It may be helpful 
;       to define REGUL=1/Delta and view Delta as the regularization error.
;     - IMPORTANT: Delta needs to be of the same order of magnitude as the typical 
;       WEIGHTS to play an effect on the regularization. One way to achieve this is: 
;       (i) divide the TEMPLATES array by a scalar in such a way that the typical 
;       template has a median of one (e.g. TEMPLATES/=median(TEMPLATES)); 
;       (ii) do the same for the input GALAXY spectrum (e.g. GALAXY/=median(GALAXY)). 
;       In this situation Delta and REGUL should be *roughly* of order unity.  
;     - Here is a possible recipe for chosing the regularization parameter REGUL:
;          (i) Perform an un-regularized fit (REGUL=0) and then rescale the input
;              NOISE spectrum so that Chi^2/DOF = Chi^2/N_ELEMENTS(goodPixels) = 1.
;              This is achieved by rescaling the input NOISE spectrum as
;              NOISE = NOISE*sqrt(Chi^2/DOF) = NOISE*sqrt(SOL[6]);
;         (ii) Increase REGUL and iteratively redo the sgandalf fit until the Chi^2
;              increases from the unregularized Chi^2 = N_ELEMENTS(goodPixels)
;              value by DeltaChi^2 = sqrt(2*n_elements(goodPixels)).
;       The derived regularization corresponds to the maximum one still consistent
;       with the observations and the derived star formation history will be the
;       smoothest (minimum curvature) that is still consistent with the observations.
;     - For a detailed explanation see Section 18.5 of Press et al. (1992,
;       Numerical Recipes 2nd ed.) available here http://www.nrbook.com/a/bookfpdf.php.
;       The adopted implementation corresponds to their equation (18.5.10).
;   VSYST: galaxy systemic velocity (zero by default). The input initial guess
;       and the output velocities are measured with respect to this velocity.
;       The value assigned to this keyword is *crucial* for the two-sided fitting.
;       In this case VSYST can be determined from a previous normal one-sided
;       fit to the galaxy velocity profile. After that initial fit, VSYST
;       can be defined as the measured velocity at the galaxy center.
;       More accurately VSYST is the value which has to be subtracted to obtain
;       a nearly anti-symmetric velocity profile at the two opposite sides of
;       the galaxy nucleus.
;     - IMPORTANT: this value is generally *different* from the systemic
;       velocity one can get from the literature. Do not try to use that!
;   WEIGHTS: a named variable to receive the value of the weights by which each stellar
;       template and ionized gas template was multiplied to best fit
;       the galaxy spectrum.
;        Stellar weights are WEIGHTS[0:Ntemplates1-1]
;        Gas emission lines intensities are  WEIGHTS[Ntemplates1 : *]
;        N.B. Gas intensities are de-reddened (i.e. the intensity in
;        the input spectrum is lower than WEIGHTS[Ntemplates1 : *],
;        because it account for the reddening at that wavelength.
;        If reddening is not fitted, then intensities and weights
;        are the same.
;

; OUTPUT PARAMETER: 
;   SOL: 9+MDEGREE elements vector containing in output the values of
;       [Vel_star,Sigma_Star,h3,h4,h5,h6,Chi^2/DOF,Vel_gas,Sigma_gas] of the best fitting solution, where DOF
;       is the number of Degrees of Freedom (number of fitted spectral pixels).
;       Vel is the velocity, Sigma is the velocity dispersion, h3-h6 are the
;       Gauss-Hermite coefficients. The model parameter are fitted simultaneously.
;     - I hardcoded the following safety limits on the fitting parameters:
;         a) Vel is constrained to be +/-2000 km/s from the first input guess
;         b) velScale/10 < Sigma < 1000 km/s
;         c) -0.3 < [h3,h4,...] < 0.3 (limits are extreme value for real galaxies)
;     - IMPORTANT: if Chi^2/DOF is not ~1 it means that the errors are not
;       properly estimated, or that the template is bad and it is *not* safe
;       to set the /CLEAN keyword.
;     - When MDEGREE > 1 then SOL contains in output the 9+MDEGREE elements
;       [Vel_star,Sigma_Star,h3,h4,h5,h6,Chi^2/DOF,Vel_gas,Sigma_gas,cx1,cx2,...,cxn], where cx1,cx2,...,cxn
;       are the coefficients of the multiplicative Legendre polynomials
;       of order 1,2,...,n. The polynomial can be explicitly evaluated as:
;           x = mdap_range(-1d,1d,n_elements(galaxy))
;           mpoly = 1d ; Multiplicative polynomial
;           for j=1,MDEGREE do mpoly += legendre(x,j)*sol[6+j]
;     - When the reddening correction is used, SOL contains 10
;       elements [Vel_star,Sigma_Star,h3,h4,h5,h6,Chi^2/DOF,Vel_gas,Sigma_gas,ebv], where ebv 
;       is the best fitting reddening value.
;
; gas_intens [dbl array]  It contains the intensities of the emission lines (corrected for reddening if the reddening is fitted).
;
; gas_fluxes  [dbl array]  It contains the fluxes of the emission lines (corrected for reddening if the reddening is fitted).
;
; gas_ew [dbl array]  It contains the equivalent widths of the emission lines (corrected for reddening if the reddening is fitted).
;
; gas_intens_err [dbl array] It contains the errors associated to gas_intens.
;
; gas_fluxes_err [dbl array] It contains the errors associated to gas_fluxes.
;
; gas_ew_err [dbl array] It contains the errors associated to gas_ew.
;
; OPTIONAL OUTPUTS
; BESTFIT: a named variable to receive a vector with the best fitting
;       model: this is a linear combination of the templates,
;       multiplied by multiplicative pols (if any) or reddening
;       corrected (if required), convolved with the best fitting
;       LOSVD, with added polynomial continuum terms
;       and the best fitting gas emission lines.
; BF_COMP1 a named variable to receive a vector with the best fitting
;       model for the stellar component: this is a linear combination
;       of the stellar templates,
;       multiplied by multiplicative pols (if any) or reddening
;       corrected (if required), convolved with the best fitting
;       LOSVD. It does not contain additive pols.
; BF_COMP2  a named variable to receive a vector with the best fitting
;        ionized gas kinematics, convolved with the gas LOSVD, and
;        mutiplied by the reddening curve (if applicable).
; OPT_TEMPL  a named variable to reveive a vector containing the optimal
;       template. This is the linear combination of the templates,
;       multiplied by multiplicative pols (if any) or reddening corrected
;       (if required), NOT CONVOLVED with the stellar LOSVD. It does
;       not contain additive pols.
; MPOLY a named variable to receive a vector with the multiplicative
;       pol that has been multipliedto BF_COMP1 and BESTFIT.
; ADDITIVE_POL a named variable to receive a vector with the additive
;      pol that has been added to BESTFIT 
; ERROR: a named variable that will contain a vector of *formal* errors
;       (1*sigma) for the fitted parameters in the output vector SOL.
;

;-
;----------------------------------------------------------------------------
FUNCTION mdap_sgandalf_reddening_curve, lambda, ebv
compile_opt idl2, hidden
;
; Reddening curve of Calzetti et al. (2000, ApJ, 533, 682; here C+00).
; This is reliable between 0.12 and 2.2 micrometres.
; - LAMBDA is the restframe wavelength in Angstrom of each pixel in the
; input galaxy spectrum (1 Angstrom = 1d-4 micrometres)
; - EBV is the assumed E(B-V) colour excess to redden the spectrum.
; In output the vector FRAC gives the fraction by which the flux at each
; wavelength has to be multiplied, to model the dust reddening effect.

k1 = lambda*0
lam = 1e4/lambda ; Convert Angstrom to micrometres and take 1/lambda
rv = 4.05d ; C+00 equation (5)

w1 = where(lambda ge 6300d, m1, COMPLEMENT=w2, NCOMPLEMENT=m2)
; C+00 equation (3) but extrapolate for lam>2.2
if m1 gt 0 then k1[w1] = rv + 2.659d*(1.040d*lam[w1] - 1.857d)
; C+00 equation (4) but extrapolate for lam<0.12
if m2 gt 0 then k1[w2] = rv + $
    2.659d*(1.509d*lam[w2] - 0.198d*lam[w2]^2 + 0.011d*lam[w2]^3 - 2.156d)
fact = 10d^(-0.4d*ebv*(k1>0))  ; Calzetti+00 equation (2) with opposite sign

return, fact ; The model spectrum has to be multiplied by this vector
END
;----------------------------------------------------------------------------
FUNCTION mdap_sgandalf_convol_fft, f, k
compile_opt idl2, hidden

nf = n_elements(f)
nk = n_elements(k)
n = 2L^ceil(alog(nf+nk/2)/alog(2))
f1 = dblarr(n)
k1 = f1
f1[0] = f
k1[0] = rotate(k,2)
k1 = shift(k1,-(nk-1)/2)
con = n*double(fft(fft(f1,-1)*fft(k1,-1),1))

return, con[0:nf-1]
END
;----------------------------------------------------------------------------
FUNCTION mdap_sgandalf_BVLS_Solve, A, b, npoly
compile_opt idl2, hidden

; No need to enforce positivity constraints if fitting one single template:
; use faster linear least-squares solution instead of BVLS.
;
s = size(a)
if s[0] eq 1 then $ ; A is a vector, not an array
    soluz = total(A*b)/total(A*A) $
else if s[2] eq npoly+1 then $ ; Fitting a single template
    soluz = la_least_squares(transpose(A),b) $
else begin               ; Fitting multiple templates
    bnd = dblarr(2,s[2],/NOZERO)
    mx = (machar()).xmax
    if npoly gt 0 then bnd[0,0:npoly-1] = -mx ; No bounds on Legendre polynomials
    bnd[0,npoly:*] = 0d  ; Positivity constraints on the templates (and sky spectra)
    bnd[1,*] = mx
    mdap_bvls, A, B, bnd, soluz, ITMAX=15*s[2], IERR=ierr
    if ierr ne 0 then message, 'BVLS Error n. ' + strtrim(ierr,2)
endelse

return, soluz
END
;----------------------------------------------------------------------------
FUNCTION mdap_sgandalf_fitfunc_optimal_template, pars, loglam=loglam,$
    BESTFIT=bestFit, BIAS=bias, CLEAN=clean, DEGREE=degree, $
    GALAXY=galaxy, GOODPIXELS=goodPixels, MDEGREE=mdegree, $
    NOISE=noise, QUIET=quiet, SKY=sky, STAR1=star1,GAS_EMS_LIST=gas_ems_list, VSYST=vsyst, WEIGHTS=weights, $
    REGUL=regul, REG_DIM=reg_dim, LAMBDA=lambda,$
    BF_COMP1 = bf_comp1, BF_COMP2 = bf_comp2,MOMENTS=moments,ADDITIVE_POL=additive_pol;,MPOLY=mpoly
compile_opt idl2, hidden


;s = size(galaxy)
;nspec = s[0]
;stop
npix = n_elements(galaxy);s[1]
;npars = n_elements(pars) - mdegree;*nspec  ; Parameters of the LOSVD only
;if n_elements(lambda) gt 1 then npars = npars - 1 ; Fitting reddening
pars_comp1 = pars[0:moments-1]
pars_comp2 = pars[moments:moments+1]

if n_elements(pars) gt moments+2  then  pars_pol = pars[moments+2:*]  ;parameters for multiplicative pols or reddening.

if n_elements(lambda) gt 1 and n_elements(pars_pol) gt 1 then stop  ; If I fit the reddening, only 1 parameter should be present (no mult. pols)
npars1=n_elements(pars_comp1)

;npars2=n_elements(pars_comp2) ==2!!
; pars = [vel,sigma,h3,h4,...,m1,m2,...]    ; Velocities are in pixels
;
dx = ceil(abs(vsyst)+abs(pars_comp1[0])+5d*pars_comp1[1]) ; Sample the Gaussian and GH at least to vsyst+vel+5*sigma
n = 2*dx + 1
x = mdap_range(dx,-dx,n)   ; Evaluate the Gaussian using steps of 1/factor pixel
losvd_comp1 = dblarr(n,/NOZERO)
;losvd_comp2 = dblarr(n,/NOZERO)

;for k=0,nspec-1 do begin    ; nspec=2 for two-sided fitting, otherwise nspec=1
s = 1d;(k eq 0) ? 1d : -1d ; s=+1 for left spectrum, s=-1 for right one
;LOSVD FOR STELLAR COMPONENT
vel = vsyst + pars_comp1[0]
w = (x - vel)/pars_comp1[1]
w2 = w*w
losvd_comp1 = exp(-0.5d*w2)/(sqrt(2d*!dpi)*pars_comp1[1]) ; Normalized total(Gaussian)=1

if npars1 gt 2 then begin
    poly = 1d + s*pars_comp1[2]/Sqrt(3d)*(w*(2d*w2-3d)) $     ; H3
              + pars_comp1[3]/Sqrt(24d)*(w2*(4d*w2-12d)+3d)   ; H4
    if npars1 eq 6 then $
        poly = poly + s*pars_comp1[4]/Sqrt(60d)*(w*(w2*(4d*w2-20d)+15d)) $    ; H5
                    + pars_comp1[5]/Sqrt(720d)*(w2*(w2*(8d*w2-60d)+90d)-15d)  ; H6
    losvd_comp1 = temporary(losvd_comp1)*poly
endif
;
x = mdap_range(-1d,1d,npix) ; X needs to be within [-1,1] for Legendre Polynomials
mpoly = 1d  ; The loop below can be null if mdegree < 1

for j=1,mdegree do  mpoly = mpoly + legendre(x,j) * pars_pol[j-1]

; Multiplicative polynomials do not make sense when fitting reddening.
; In that case one has to assume the spectrum is well calibrated.
;
;if n_elements(lambda) gt 1 then mpoly = mdap_sgandalf_reddening_curve(lambda, pars[npars])
if n_elements(lambda) gt 1 then  mpoly = mdap_sgandalf_reddening_curve(lambda, pars_pol[0])

; Fill the columns of the design matrix of the least-squares problem
;
s1 = size(star1)
;s2 = size(gas_ems_list)
;ss = size(sky)
;nsky = (ss[0] lt 2) ? ss[0] : ss[2] ; Number of sky spectra
ntemp1 = (s1[0] eq 2) ? s1[2] : 1     ; Number of template spectra for stars
; ntemp2 = (s2[0] eq 2) ? s2[2] : 1     ; Number of template spectra for emission lines
ntemp2 = n_elements(gas_ems_list[*,0])
;nrows = (degree + 1 + nsky)*nspec + ntemp
nrows1 = degree + 1 + ntemp1
nrows2 = ntemp2
ncols = npix
if regul gt 0 then begin
    nr = n_elements(reg_dim)
    reg2 = reg_dim - 2
    case nr of
        1: nreg = reg2
        2: nreg = 2*product(reg2) + 2*total(reg2) ; Rectangle sides have one finite difference
        3: nreg = 3*product(reg2) + 4*total(reg2) $ ; Hyper-rectangle edges have one finite difference
                + 4*( product(reg2[[0,1]]) + product(reg2[[0,2]]) + product(reg2[[1,2]]) )
    endcase
    ncols = ncols + nreg
endif
c_comp1 = dblarr(ncols,nrows1,/nozero)  ; This array is used for estimating predictions
c_comp2 = dblarr(ncols,nrows2,/nozero) 

;for j=0,degree do $ ; Fill first columns of the Design Matrix
;    if nspec eq 2 then begin
;        leg = legendre(x,j)
;        c[0,2*j] = [leg,leg*0d]   ; Additive polynomials for left spectrum
;        c[0,2*j+1] = [leg*0d,leg] ; Additive polynomials for right spectrum
;    endif else c[0,j] = legendre(x,j)

; Fill first columns of the Design Matrix
for j=0,degree do  c_comp1[0,j] = legendre(x,j)  ;only for stellar component,

;pix = mdap_range(0d,s[1]-1d,s[1])
;if factor gt 1 then pix = mdap_range(0d,s[1]-1d,s[1]*factor) ; Oversampled pixels range
;tmp = dblarr(s[1],nspec,/NOZERO)
tmp_comp1 = dblarr(s1[1],/NOZERO)
;for j=0,ntemp1-1 do begin
;    if factor eq 1 then $ ; No oversampling of the template spectrum
;        for k=0,nspec-1 do tmp[*,k] = mdap_sgandalf_convol_fft(star[*,j],losvd[*,k]) $
;    else begin             ; Oversample the template spectrum before convolution
;        st = interpolate(star[*,j],pix,CUBIC=-0.5)   ; Sinc-like interpolation
;        for k=0,nspec-1 do tmp[*,k] = rebin(mdap_sgandalf_convol_fft(st,losvd[*,k]),s[1])
;    endelse
;    c[0,(degree+1)*nspec+j] = (mpoly*tmp[0:npix-1,*])[*] ; reform into a vector
;endfor

for j=0,ntemp1-1 do begin
                ; Oversample the template spectrum before convolution
    ;st = interpolate(star1[*,j],pix,CUBIC=-0.5)   ; Sinc-like interpolation
    tmp_comp1 = mdap_sgandalf_convol_fft(star1[*,j],losvd_comp1)
    c_comp1[0,(degree+1)+j] = (mpoly*tmp_comp1[0:npix-1,*])[*] ; reform into a vector
endfor

; tmp_comp2 = dblarr(s2[1],/NOZERO)
; for j=0,ntemp2-1 do begin
;     ;st = interpolate(gas_ems_list[*,j],pix,CUBIC=-0.5)   ; Sinc-like interpolation
;     tmp_comp2 = mdap_sgandalf_convol_fft(gas_ems_list[*,j],losvd_comp2)
;     ;stop
;     c_comp2[0,j] = tmp_comp2[0:npix-1];,*]);[*] ; reform into a vector
;  ;   c_comp2[0,j] = 0.
;    ; if max( c_comp2[*,j]) ne 0 then c_comp2[*,j] = c_comp2[*,j] / max(c_comp2[*,j]);tmp_comp2[0:npix-1])
;     
; endfor
;stop
;tmp_comp2 = dblarr(s1[1],/NOZERO)
x = dindgen(npix) 
;stop

pos = gas_ems_list[*,0]+pars_comp2[0]

w2=2.*pars_comp2[1]*pars_comp2[1]


for j=0,ntemp2-1 do begin
   if pos[j] ge 0 and pos[j] le npix then $
         
          c_comp2[0,j] = gas_ems_list[j,1] * exp(-(x-pos[j])^2/w2) else c_comp2[0,j] = dblarr(npix);begin
  
   if n_elements(lambda) gt 1 then c_comp2[0,j] = c_comp2[*,j] * mpoly
endfor


;for j=0,ntemp2-1 do c_comp2[*,j] = exp(-(x-(gas_ems_list[j]+pars_comp2[0]))^2./(2.*pars_comp2[1]^2.))
 ;stop
; Add second-degree 1D, 2D or 3D linear regularization
; Press W.H., et al., 1992, Numerical Recipes, 2nd ed. equation (18.5.10)
;
if regul gt 0 then begin
    i = indgen(reg_dim) + (degree+1);*nspec
    p = npix;*nspec
    diff = [-1d,2d,-1d]*regul
    dim = size(reg_dim,/DIM) > 1
    case dim of
        1: for j=1,reg_dim-2 do c[p++,i[j+[-1,0,1]]] = diff
        2: for k=0,reg_dim[1]-1 do $
               for j=0,reg_dim[0]-1 do begin
                   if j ne 0 && j ne reg_dim[0]-1 then c_comp1[p++,i[j+[-1,0,1],k]] = diff
                   if k ne 0 && k ne reg_dim[1]-1 then c_comp1[p++,i[j,k+[-1,0,1]]] = diff
               endfor
        3: for l=0,reg_dim[2]-1 do $
               for k=0,reg_dim[1]-1 do $
                   for j=0,reg_dim[0]-1 do begin
                       if j ne 0 && j ne reg_dim[0]-1 then c_comp1[p++,i[j+[-1,0,1],k,l]] = diff
                       if k ne 0 && k ne reg_dim[1]-1 then c_comp1[p++,i[j,k+[-1,0,1],l]] = diff
                       if l ne 0 && l ne reg_dim[2]-1 then c_comp1[p++,i[j,k,l+[-1,0,1]]] = diff
                   endfor
    endcase
endif

;for j=0,nsky-1 do begin
;    skyj = sky[*,j]
;    k = (degree+1)*nspec + ntemp
;    if nspec eq 2 then begin
;        c[0,k+2*j] = [skyj,skyj*0d]   ; Sky for left spectrum
;        c[0,k+2*j+1] = [skyj*0d,skyj] ; Sky for right spectrum
;    endif else c[0,k+j] = skyj
;endfor



a_comp1 = dblarr(ncols,nrows1,/nozero)  ; This array is used for estimating predictions
a_comp2 = dblarr(ncols,nrows2,/nozero) 
;for j=0,nrows1-1 do a_comp1[0,j] = c_comp1[0:npix-1,j]/noise ; Weight all columns with errors
;for j=0,nrows2-1 do a_comp2[0,j] = c_comp2[0:npix-1,j]/noise ; Weight all columns with errors
for j=0,nrows1-1 do a_comp1[0,j] = c_comp1[*,j]/noise ; Weight all columns with errors
for j=0,nrows2-1 do a_comp2[0,j] = c_comp2[*,j]/noise ; Weight all columns with errors

;if regul gt 0 then begins
;    aa = a[[goodPixels,mdap_range(npix*nspec,ncols-1)],*]
;    bb = [galaxy[goodPixels]/noise[goodPixels],replicate(0d,nreg)]
;endif else begin
;    aa = a[goodPixels,*]
;    bb = galaxy[goodPixels]/noise[goodPixels]
;endelse
;stop
if regul gt 0 then begin
   ; aa_comp1 = a_comp1[[goodPixels,mdap_range(npix,ncols-1)],*]
   ; aa_comp2 = a_comp2[[goodPixels,mdap_range(npix,ncols-1)],*]
   ; aa1 = [[aa_comp1],[aa_comp2*0.]]
   ; aa2 = [[aa_comp1*0.],[aa_comp2]]
   ; aa = temporary(aa1)+temporary(aa2)
    aa = [[a_comp1[[goodPixels,mdap_range(npix,ncols-1)],*]],[a_comp2[[goodPixels,mdap_range(npix,ncols-1)],*]]]
    bb = [galaxy[goodPixels]/noise[goodPixels],replicate(0d,nreg)]
endif else begin
   ; aa1 = [[a_comp1[goodPixels,*]],[a_comp2[goodPixels,*]*0.]]
   ; aa2 = [[a_comp1[goodPixels,*]*0.],[a_comp2[goodPixels,*]]]
   ; stop
    ;aa1 = [[a_comp1[goodPixels,*]],[dblarr(n_elements(goodPixels),ntemp2)]]
    ;aa2 = [[dblarr(n_elements(goodPixels),ntemp1+degree+1)],[a_comp2[goodPixels,*]]]
    ;aa=temporary(aa1)+temporary(aa2)
    aa = [[a_comp1[goodPixels,*]],[a_comp2[goodPixels,*]]];[dblarr(n_elements(goodPixels),ntemp2)]]
    bb = galaxy[goodPixels]/noise[goodPixels]
   ; stop
endelse


; Select the spectral region to fit and solve the overconditioned system
; using SVD/BVLS. Use unweighted array for estimating bestfit predictions.
; Iterate to exclude pixels deviating more than 3*sigma if /CLEAN keyword is set.
;
npoly = (degree+1);*nspec ; Number of additive polynomials in the fit

repeat begin
   
    weights = mdap_sgandalf_BVLS_Solve(aa,bb,npoly)
    
   ; weights[ntemp1+5]=2.
   ; if npoly 

    bf_comp2 = c_comp2 #  weights[npoly+ntemp1 : *]
    if npoly eq 0 then begin
       bf_comp1 = c_comp1 #  weights[0 : npoly+ntemp1-1] 
       additive_pol = fltarr(npix);bf_comp1*0.
       bestfit = bf_comp1 + bf_comp2 ;
    endif else begin
       bf_comp1 = c_comp1[*,npoly : *] #  weights[npoly : npoly+ntemp1-1] 
       additive_pol = c_comp1[*,0:npoly-1] # weights[0:npoly-1]
       bestfit = bf_comp1 + bf_comp2 + additive_pol
    endelse
     
   ; bf_comp2 = c_comp2 #  weights[npoly+ntemp1 : *]
   ; bestfit = bf_comp1 + bf_comp2 + ;

;    plot,galaxy
;    oplot, bf_Comp1,color=200
;    oplot, bf_comp2,color=200
;stop
   ; weights1 = mdap_sgandalf_BVLS_Solve(a_comp1[goodPixels,*],bb,npoly)
   ; bf_comp1 = c_comp1 #  weights1
   ; weights2 = mdap_sgandalf_BVLS_Solve(a_comp2[goodPixels,*],(galaxy[goodPixels]-bf_comp1)/noise[goodPixels],npoly)
    
   ; bf_comp1 = c_comp1 #  weights[0:ntemp1-1]
  ;  bf_comp2 = c_comp2 #  weights2;[ntemp1:*]
  ;  bestfit = bf_comp1 + bf_comp2

  ;  weights=[weights1,weights2]
;stop
   ; if regul gt 0 then begin
   ;    bb = [(galaxy[goodPixels]-bf_comp1)/noise[goodPixels],replicate(0d,nreg)]
   ; endif else begin
   ;    bb=(galaxy[goodPixels]-bf_comp1)/noise[goodPixels]
   ; endelse
   ; weights2 = mdap_sgandalf_BVLS_Solve(a_comp2[goodPixels,*],bb,npoly)
   ; weights[ntemp1]=weights2
   ; bf_comp2 = c_comp2 #  weights[ntemp1:*]
   
 ;   bestfit = bf_comp1 + bf_comp2
;    bestfit = c[0:npix*nspec-1,*] # weights
    err = (galaxy[goodPixels]-bestfit[goodPixels])/noise[goodPixels]
    if keyword_set(clean) then begin
        tmp = where(abs(err) gt 3, m, COMPLEM=w) ; select errors larger than 3*sigma
        if (m ne 0) then begin
            if ~keyword_set(quiet) then print, 'Outliers:', m
            goodPixels = goodPixels[w]
        endif
    endif else break
endrep until (m eq 0)

; Penalize the solution towards (h3,h4,...)=0 if the inclusion of
; these additional terms does not significantly decrease the error.
; Only for the stellar component
if npars1 gt 2 && bias ne 0 then $
    err = err + bias*robust_sigma(err, /ZERO)*sqrt(total(pars_comp1[2:*]^2))

return, err
END
;----------------------------------------------------------------------------
PRO mdap_sgandalf, templates1,templates2, loglam, galaxy, noise, velScale, start, sol, $
    gas_intens,gas_fluxes,gas_ew,gas_intens_err,gas_fluxes_err,gas_ew_err,$
    BESTFIT=bestFit, BIAS=bias, CLEAN=clean, DEGREE=degree, ERROR=error, $
    GOODPIXELS=goodPixels, MDEGREE=mdegree, MOMENTS=moments, $
    POLYWEIGHTS=polyweights, QUIET=quiet, VSYST=vsyst, WEIGHTS=weights, $
    REGUL=regul, LAMBDA=lambda, REDDENING=reddening,ADDITIVE_POL=additive_pol,$
    BF_COMP1 = bf_comp1, BF_COMP2 = bf_comp2, OPT_templ=opt_templ,mpoly=mpoly,$;, PLOT=plot
    fix_star_kin=fix_star_kin,fix_gas_kin=fix_gas_kin,$
    range_v_star=range_v_star,range_s_star=range_s_star,range_v_gas=range_v_gas,range_s_gas=range_s_gas,$
    emission_lines_fluxes=emission_lines_fluxes

compile_opt idl2
on_error, 2

; Do extensive checking of possible input errors
;
; templates1 : stellar spectra
; templates2 : gaussian emission lines for ionized gas fitting.
;
c=299792.458d

s1 = size(templates1)
if s1[0] ge 3 then begin
    reg_dim = s1[2:s1[0]]
    star1 = reform(templates1,s1[1],product(reg_dim))
    s1 = size(star1)
endif else begin
    reg_dim = s1[2]
    star1 = templates1
endelse

gas_ems_list = templates2
x = dindgen(n_elements(galaxy))
gas_ems_list_pix= [[templates2[*,0]*0.],[templates2[*,1]]];dindgen(n_elements(galaxy),2)
gas_ems_list_pix[*,0] = (gas_ems_list[*,0]-min(loglam))*max(x) / (max(loglam)-min(loglam))

if n_elements(regul) eq 0 then regul = 0
s2 = size(galaxy)
s3 = size(noise)
;s4 = size(sky)
if (s1[0] gt 2 || s2[0] gt 2 || s3[0] gt 2) then message, 'Wrong input dimensions'
if ~array_equal(s2,s3) then message, 'GALAXY and NOISE must have the same size/type'
if (s1[1] lt s2[1]) then message, 'STAR length cannot be smaller than GALAXY'
if n_elements(reddening) gt 0 then begin
    if ~array_equal(size(lambda),s2) then $
        message, 'LAMBDA and GALAXY must have the same size/type'
    if n_elements(mdegree) gt 0 then $
        message, 'MDEGREE cannot be used with REDDENING keyword'
endif else lambda = 0
if ~array_equal(noise gt 0, 1) then message, 'NOISE must be a positive vector'
;if s4[0] gt 0 && s4[1] ne s2[1] then message, 'SKY must have the same size as GALAXY'
if n_elements(degree) le 0 then degree = -1 else degree = degree > (-1)
if n_elements(mdegree) le 0 then mdegree = 0 else mdegree = mdegree > 0
;if keyword_set(oversample) then factor = 30 else factor = 1
nGood = n_elements(goodPixels)
if nGood le 0 then begin
    nGood = s2[1]
    goodPixels = indgen(nGood)
endif else begin
    if ~array_equal((goodPixels[1:*] - goodPixels) gt 0, 1) then $
        message, 'GOODPIXELS is not monotonic or contains duplicated values'
    if goodPixels[0] lt 0 || goodPixels[nGood-1] gt s2[1]-1 then $
        message, 'GOODPIXELS are outside the data range'
endelse
if n_elements(bias) le 0 then bias = 0.7d*sqrt(500d/n_elements(goodPixels)) ; pPXF paper pg.144 left
if n_elements(moments) eq 0 then moments = 2 else $
    if total(moments eq [0,2,4,6]) eq 0 then message, 'MOMENTS should be 0, 2, 4 or 6'
;if moments ge 2 && n_elements(start) ne 4 then message, 'START must have 4 elements [V_star,sigma_star,V_gas,sigma_gas]'
;if s2[0] eq 2 then goodPixels = [goodPixels,s2[1]-1+goodPixels]  ; two-sided fitting of LOSVD
if n_elements(vsyst) eq 0 then begin
  ;  if s2[0] eq 2 then message, 'VSYST must be defined for two-sided fitting'
    vsyst = 0d
endif

gaspars = 2
if moments gt 0 then begin

    ; Explicitly specify the step for the numerical derivatives (pixel/100)
    ; in MPFIT routine and force safety limits on the fitting parameters.
    ;
    if mdegree gt 0 then begin
        start1 = dblarr(moments+gaspars+mdegree)  ; Set [h3,h4,...] to zero as initial guess
        parinfo = replicate({step:1d-2,limits:[0d,0d],limited:[1,1],fixed:0}, moments+gaspars+mdegree)
        parinfo[moments+gaspars:*].limits = [-1d,1d] ; force <100% corrections
        parinfo[moments+gaspars:*].step = 1d-3
    endif else if n_elements(reddening) gt 0 then begin
        start1 = dblarr(moments+gaspars+1)  ; Set [h3,h4,...] to zero as initial guess
        start1[moments+gaspars] = reddening   ; No multiplicative pols are fitted.
        parinfo = replicate({step:1d-2,limits:[0d,0d],limited:[1,1],fixed:0}, moments+gaspars+1)
        parinfo[moments+gaspars].limits = [0d,10d] ; force positive E(B-V) < 10 mag
        parinfo[moments+gaspars].step = 1d-3
    endif else begin
        start1 = dblarr(moments+gaspars)  ; Set [h3,h4,...] to zero as initial guess
        parinfo = replicate({step:1d-2,limits:[0d,0d],limited:[1,1],fixed:0}, moments+gaspars)
    endelse
    start1[0] = start[0:1]/velScale  ; Convert velocity scale to pixels
    start1[moments] =  start[4:5]/velScale   ; start[moments:moments+1]/velScale  
    parinfo[0].limits = start1[0] + [-2d3,2d3]/velScale ; hard-coded +/-2000 km/s from first guess
    parinfo[1].limits = [0.2d,5d2/velScale] ; hard-coded velScale/20<sigma<500 km/s

    if n_elements(range_v_star) gt 0 then parinfo[0].limits = start1[0] + [range_v_star[0],range_v_star[1]]/velScale 
    if n_elements(range_s_star) gt 0 then parinfo[1].limits = [max([0.2d,start1[1] + range_s_star[0]/velscale]), (start1[1] + range_s_star[1]/velScale)] 

    parinfo[moments].limits = start1[moments] + [-2d3,2d3]/velScale ; hard-coded +/-2000 km/s from first guess
    parinfo[moments+1].limits = [0.1d,5d2/velScale] ; hard-coded velScale/20<sigma<500 km/s

    if n_elements(range_v_gas) gt 0 then parinfo[moments].limits = start1[moments] + [range_v_gas[0],range_v_gas[1]]/velScale 
    if n_elements(range_s_gas) gt 0 then parinfo[moments+1].limits = [max([0.1d,start1[moments+1] + range_s_gas[0]/velscale]), (start1[moments+1] + range_s_gas[1]/velScale)] 


    if moments gt 2 then begin
        parinfo[2:moments-1].limits = [-0.4d,0.4d] ; -0.4<[h3,h4,...]<0.4
        parinfo[2:moments-1].step = 1d-3
       
        start1[2:moments-1] = start[2:moments-1]
     endif
    ;if keyword_set(fix_star_kin) then stop
    if keyword_set(fix_star_kin) then parinfo[0:moments-1].fixed = 1
    if keyword_set(fix_gas_kin) then parinfo[moments:moments+1].fixed = 1
    if s1[1] le 2*(abs(vsyst/velScale)+abs(start1[0])+5d*start1[1]) then $
        message, 'Velocity shift too big: Adjust wavelength ranges of spectrum and templates'


    ;start1 =[ Vstar, Sstar, <h3,h4,h5,h6>,Vgas,Sgas,<M1,... Mmdegree, or reddening>]

    ; Here the actual calculation starts.
    ; If required, once the minimum is found, clean the pixels deviating
    ; more than 3*sigma from the best fit and repeat the minimization
    ; until the set of cleaned pixels does not change any more.
    ;
    good = goodPixels
    for j=0,4 do begin ; Do at most five cleaning iterations
        ;if n_elements(sky) eq 0 then $
            functArgs = {BIAS:bias, DEGREE:degree,  GALAXY:double(galaxy), $
                GOODPIXELS:goodPixels, MDEGREE:mdegree, NOISE:double(noise), $
                STAR1:double(star1),GAS_EMS_LIST:double(gas_ems_list_pix), VSYST:vsyst/velScale, $
                REGUL:regul, REG_DIM:reg_dim, LAMBDA:lambda,MOMENTS:moments,LOGLAM:loglam} 
       ; else $
       ;     functArgs = {BIAS:bias, DEGREE:degree, FACTOR:factor, GALAXY:double(galaxy), $
       ;         GOODPIXELS:goodPixels, MDEGREE:mdegree, NOISE:double(noise), $
       ;         SKY:double(sky), STAR:double(star), VSYST:vsyst/velScale, $
       ;         REGUL:regul, REG_DIM:reg_dim, LAMBDA:lambda}
        res = mpfit('mdap_sgandalf_fitfunc_optimal_template', start1, ERRMSG=errmsg, $
            FTOL=1d-4, FUNCTARGS=functArgs, NFEV=ncalls, PARINFO=parinfo, $
            PERROR=error, /QUIET)
        if errmsg ne '' then message, errmsg
        if ~keyword_set(clean) then break
        goodOld = goodPixels
        goodPixels = good
        tmp = mdap_sgandalf_fitfunc_optimal_template(res, BIAS=bias, /CLEAN, $
            DEGREE=degree, GALAXY=galaxy, $
            GOODPIXELS=goodPixels, MDEGREE=mdegree, NOISE=noise, $
            QUIET=quiet, STAR1=double(star1),GAS_EMS_LIST=double(gas_ems_list_pix), VSYST=vsyst/velScale, $
            REGUL=regul, REG_DIM=reg_dim, LAMBDA=lambda,MOMENTS=moments,LOGLAM=loglam)
        if array_equal(goodOld,goodPixels) then break
    endfor

endif else begin

    if mdegree gt 0 then message, 'MDEGREE>0 not implemented with MOMENTS=0'
    res = start
    res[0] = res[0:1]/velScale  ; Convert velocity scale to pixels
    res[moments] = res[moments:moments+1]/velScale  ; Convert velocity scale to pixels
    error = start*0d              ; No output error on parameters
    ncalls = 1

endelse

; Evaluate scatter at the bestfit (with BIAS=0)
; and also get the output BESTFIT and WEIGHTS.
;
err = mdap_sgandalf_fitfunc_optimal_template(res, BESTFIT=bestFit, BIAS=0, DEGREE=degree, $
     GALAXY=galaxy, GOODPIXELS=goodPixels, MDEGREE=mdegree, $
    NOISE=noise,  STAR1=double(star1),GAS_EMS_LIST=double(gas_ems_list_pix), VSYST=vsyst/velScale, WEIGHTS=weights, $
    REGUL=regul, REG_DIM=reg_dim, LAMBDA=lambda,BF_COMP1 = bf_comp1, BF_COMP2 = bf_comp2,MOMENTS=moments,LOGLAM=loglam,ADDITIVE_POL=additive_pol)
chi2 = robust_sigma(err, /ZERO)^2       ; Robust computation of Chi^2/DOF.
sol = dblarr(7+gaspars+mdegree)


;stop
;sol[0] = res[0:(n_elements(start)>moments)-1]
;sol[0] = res[0:1];kmoments+gaspar-1]
if moments eq 0 then begin
   sol[0:5] = [0,0,0,0,0,0]
;   sol[6] = chi2
   sol[7:8] = res[0:1]*velScale
  if mdegree gt 0 then sol[9:*] = res[2:*]
endif

if moments eq 2 then begin
   sol[0:1] = res[0:1]*velScale
   sol[2:5] = [0,0,0,0]
;   sol[6] = chi2
   sol[7:8] = res[2:3]*velScale
  if mdegree gt 0  then sol[9:*] = res[4:*]
endif

if moments eq 4 then begin
   sol[0:1] = res[0:1]*velScale
   sol[2:3] = res[2:3]
   sol[4:5] = [0,0]
;   sol[6]=chi2
   sol[7:8] = res[4:5]*velscale
   if mdegree gt 0 then sol[9:*] = res[6:*]
endif

if moments eq 6 then begin
   sol[0:1] = res[0:1]*velScale
   sol[2:5] = res[2:5]
;   sol[6]=chi2
   sol[7:8] = res[6:7]*velscale
   if mdegree gt 0 then sol[9:*] = res[8:*]
endif


;sol[0] = sol[0:1]*velScale ; Bring velocity scale back to km/s
;sol[moments] = sol[moments:moments+1]*velScale ; Bring velocity scale back to km/s

error[0:1] = error[0:1]*velScale ; Convert errors to km/s
error[moments:moments+1] = error[moments:moments+1]*velScale ;
sol[6] = chi2
;stop
;if mdegree ge 1 then sol[7] = res[moments:*]
if n_elements(reddening) gt 0 then reddening = res[moments+gaspars] ; Replace input with best fit
if degree ge 0 then polyweights = weights[0:(degree+1)-1] ; output weights for the additive polynomials
weights = weights[(degree+1):*] ; output weights for the templates (and emission lines) only


;weights[sztempl[2]:*]*emss

; Print final results on the screen.
;
sz=size(stars1)
if ~keyword_set(quiet) then begin
    print, 'feval', 'V', 'sigma', 'h3', 'h4', 'h5', 'h6','Chi2/DOF', 'Vgas','Sgas', FORMAT='(10A10)'
    print, ncalls, sol[0:8], FORMAT='(i10,2f10.1,5f10.3,2f10.1)'
    nw = n_elements(weights)
    if n_elements(reddening) gt 0 then print, 'Reddening E(B-V): ', reddening, FORMAT='(a,g0.3)'
    ;print, 'Stellar weights: ', total(weights gt 0), ' / ', nw, FORMAT='(a,g0,a,g0)'
    ;if n_elements(weights) le 20 then begin
    if degree ge 0 then begin
       print,'Additive pol. weights:'
       print, polyweights
    endif
    if mdegree gt 0 then begin
       print,'Multiplicative pol. coeff:'
       print, sol[9:*]
    endif

        print, 'Templates weights:'
        print, weights[0:s1[2]-1];, FORMAT='(10g11.3)'
    ;endif
        print, ''
        print,'GAS Intensities:'
        print,''
        print, weights[s1[2]:*]*templates2[*,1];, FORMAT='(10g11.3)'
endif

;Construct optimal stellar template, not LOSVD convolved.
x = mdap_range(-1d,1d,s2[1])
mpoly = 1d                   ; Multiplicative polynomial
if lambda[0] eq 0 then for j=1,MDEGREE do mpoly += legendre(x,j)*sol[7+gaspars-1+j]
if n_elements(lambda) gt 1 then mpoly =  mdap_sgandalf_reddening_curve(lambda, reddening[0])
templates1_=templates1
for i = 0, sz[2]-1 do templates1_[0,i] = templates1[*,i]*mpoly ;multipl pols only on template1, because templ2 is gas.
opt_templ = templates1_ #  weights[0:s1[2]-1]



gas_intens=weights[s1[2]:*]*templates2[*,1]
gas_intens_err = gas_intens*.0

step_ln = loglam[1]-loglam[0]
line_observed_position = (alog(exp(templates2[*,0])+sol[7]/c*exp(templates2[*,0])));  ;log-wave range were to compute statistics (depends on the line FWHM)
interval = 2.3548*sol[8]/velscale*step_ln  ;FWHM in log_angstrom

for j = 0, n_elements(gas_intens)-1 do begin
   ;stop
    indici = where(loglam ge line_observed_position[j]-3.*interval and loglam le line_observed_position[j]+4.*interval)   ;  3*FWHM < line < 4*FWHM
    if indici[0] ne -1 then gas_intens_err[j] = stddev(galaxy[indici]-bestFit[indici]) 
    ;stop
endfor


;---------------------------------------------------------------------------------------- 
; NO NEED TO CORRECT THE EMISSION LINES FOR REDDENING 
; IT IS ALREADY DONE WITHIN mdap_sgandalf_fitfunc_optimal_template
; see line: 
;       if n_elements(lambda) gt 1 then c_comp2[0,j] = c_comp2[*,j] * mpoly
;
; if n_elements(reddening) gt 0 then begin
;    l0_gal=min(loglam)
;    lstep_gal=loglam[1]-loglam[0]
;    Vstar = sol[0]
;    reddening_attenuation = $
;       mdap_dust_calzetti(l0_gal,lstep_gal,n_elements(galaxy),reddening,Vstar)        
;    ;reddening_attenuation = $
;    ;   mdap_dust_calzetti(l0_gal,lstep_gal,n_elements(galaxy),reddening,Vstar,\log10)        
;    ob_lambda = exp(dindgen(n_elements(galaxy))*lstep_gal + l0_gal)
;    rf_lambda = ob_lambda/exp(Vstar/c)
;    reddening_attenuation_emission = $
;       interpol(reddening_attenuation,rf_lambda,exp(templates2[*,0]))
;    gas_intens = gas_intens * reddening_attenuation_emission
; endif
;---------------------------------------------------------------------------------------- 

;gas_fluxes = gas_intens * sol[8]/velscale*sqrt(2.*!pi)
gas_fluxes = gas_intens * sol[8] * sqrt(2*!pi) * exp(templates2[*,0]) * exp(sol[7]/c)/c   ;same formula in gandalf
gas_fluxes_err = gas_intens_err * sol[8] * sqrt(2*!pi) * exp(templates2[*,0]) * exp(sol[7]/c)/c   ;not sure about the formula, check.
;-- compute gas EW (as in gandalf)
rf_l  = exp(loglam -sol[0]/velscale*step_ln)
gas_ew=gas_fluxes*0.
gas_ew_err=gas_fluxes*0.
for j = 0, n_elements(templates2[*,0])-1 do begin
    rf_l_line  = exp(templates2[j,0])
    S_line     = sol[8] 
    F_obs_line = gas_fluxes[j]
    j_buffer   = where(abs(rf_l-rf_l_line) lt 10*(S_line/c)*rf_l_line and abs(rf_l-rf_l_line) gt  5*(S_line/c)*rf_l_line)
    if j_buffer[0] ne -1 then begin
       C_line     = median(bf_comp1[j_buffer])
       gas_ew[j]  = gas_fluxes[j]/C_line
       gas_ew_err[j]  = gas_fluxes_err[j]/C_line+abs(gas_fluxes[j])/C_line^2.*robust_sigma(bf_comp1[j_buffer])
    endif
endfor
;--
if n_elements(reddening) and n_elements(lambda) gt 0 then junk = temporary(mdegree)
END
;----------------------------------------------------------------------------
