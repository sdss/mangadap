;
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
;   PPXF
;
; PURPOSE:
;   Extract galaxy stellar kinematics (V, sigma, h3, h4, h5, h6)
;   or the stellar population by fitting a template to an observed
;   spectrum in pixel space, using the Penalized Pixel-Fitting (pPXF)
;   method described in: Cappellari M., & Emsellem E., 2004, PASP, 116, 138.
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
;   5) It is possible to fit a mirror-symmetric LOSVD to two spectra at
;      the same time. This is useful for spectra taken at point-symmetric
;      spatial positions with respect to the center of an equilibrium
;      stellar system.
;   6) One can include sky spectra in the fit, to deal with cases
;      where the sky dominates the observed spectrum and an accurate
;      sky subtraction is critical.
;   7) One can derive an estimate of the reddening in the spectrum.
;
; CALLING SEQUENCE:
;   PPXF, templates, galaxy, noise, velScale, start, sol, $
;       BESTFIT=bestFit, BIAS=bias, /CLEAN, DEGREE=degree, ERROR=error, $
;       GOODPIXELS=goodPixels, LAMBDA=lambda, MDEGREE=mdegree, MOMENTS=moments, $
;       /OVERSAMPLE, /PLOT, POLYWEIGHTS=polyWeights, /QUIET, REDDENING=reddening, $
;       REGUL=regul, SKY=sky, VSYST=vsyst, WEIGHTS=weights
;
; INPUT PARAMETERS:
;   TEMPLATES: vector containing the spectrum of a single template star or more
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
;       an error will be returned by PPXF, if the velocity shift in pixels,
;       required to match the galaxy with the templates, becomes larger than
;       nPixels. In that case one has to truncate either the galaxy or the
;       templates to make the two rest-frame spectral ranges more similar.
;   GALAXY: vector containing the spectrum of the galaxy to be measured. The
;       star and the galaxy spectra have to be logarithmically rebinned but the
;       continuum does *not* have to be subtracted. The rebinning may be
;       performed with the LOG_REBIN routine that is distributed with PPXF.
;     - For high redshift galaxies, one should bring the spectra close to the
;       restframe wavelength, before doing the PPXF fit, to prevent too large
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
;     - IMPORTANT: the penalty term of the pPXF method is based on the *relative*
;       change of the fit residuals. For this reason the penalty will work as
;       expected even if no reliable estimate of the NOISE is available
;       (see Cappellari & Emsellem [2004] for details).
;       If no reliable noise is available this keyword can just be set to:
;           NOISE = galaxy*0+1 ; Same weight for all pixels
;   VELSCALE: velocity scale of the spectra in km/s per pixel. It has to be the
;       same for both the galaxy and the template spectra.
;   START: two elements vector [velStart,sigmaStart] with the initial estimate
;       for the velocity and the velocity dispersion in km/s.
;     - Unless a good initial guess is available, it is recommended to set the starting
;       sigma >= 3*velScale in km/s (i.e. 3 pixels). In fact when the LOSVD is severely
;       undersampled, and far from the true solution, the chi^2 of the fit becomes weakly
;       sensitive to small variations in sigma (see pPXF paper). In some instances the
;       near-constancy of chi^2 may cause premature convergence of the optimization.
;     - In the case of two-sided fitting a good starting value for the
;       velocity is velStart=0.0 (in this case VSYST will generally be nonzero).
;       Alternatively on should keep in mind that velStart refers to the first
;       input galaxy spectrum, while the second will have velocity -velStart.
;
; KEYWORDS:
;   BESTFIT: a named variable to receive a vector with the best fitting
;       template: this is a linear combination of the templates, convolved with
;       the best fitting LOSVD, with added polynomial continuum terms.
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
;       the template continuum shape during the fit (default: 4).
;       Set DEGREE = -1 not to include any additive polynomial.
;   ERROR: a named variable that will contain a vector of *formal* errors
;       (1*sigma) for the fitted parameters in the output vector SOL. This 
;       option can be used when speed is essential, to obtain an order of 
;       magnitude estimate of the uncertainties, but we *strongly* recommend to 
;       run Monte Carlo simulations to obtain more reliable errors. In fact these 
;       errors can be severely underestimated in the region where the penalty 
;       effect is most important (sigma < 2*velScale).
;     - These errors are meaningless unless Chi^2/DOF~1 (see parameter SOL below).
;       However if one *assume* that the fit is good, a corrected estimate of the
;       errors is: errorCorr = error*sqrt(chi^2/DOF) = error*sqrt(sol[6]).
;     - IMPORTANT: when running Monte Carlo simulations to determine the error,
;       the penalty (BIAS) should be set to zero, or better to a very small value.
;       See Section 3.4 of Cappellari & Emsellem (2004) for an explanation.
;   GOODPIXELS: integer vector containing the indices of the good pixels in the
;       GALAXY spectrum (in increasing order). Only these pixels are included in
;       the fit. If the /CLEAN keyword is set, in output this vector will be updated
;       to contain the indices of the pixels that were actually used in the fit.
;     - IMPORTANT: in all likely situations this keyword *has* to be specified.
;   LAMBDA: When the keyword REDDENING is used, the user has to pass in this
;       keyword a vector with the same dimensions of GALAXY, giving the restframe
;       wavelength in Angstrom of every pixel in the input galaxy spectrum.
;       If one uses my LOG_REBIN routine to rebin the spectrum before the PPXF fit:
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
;   /OVERSAMPLE: Set this keyword to oversample the template by a factor 30x
;       before convolving it with a well sampled LOSVD. This can be useful to
;       extract proper velocities, even when sigma < 0.7*velScale and the 
;       dispersion information becomes totally unreliable due to undersampling. 
;       IMPORTANT: One should sample the spectrum more finely is possible, 
;       before resorting to the use of this keyword! 
;   /PLOT: set this keyword to plot the best fitting solution and the residuals
;       at the end of the fit.
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
;       function PPXF_REDDENING_CURVE below.
;     - IMPORTANT: The MDEGREE keyword cannot be used when REDDENING is set.
;   REGUL: If this keyword is nonzero, the program applies second-degree
;       linear regularization to the WEIGHTS during the PPXF fit.
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
;         (ii) Increase REGUL and iteratively redo the pPXF fit until the Chi^2
;              increases from the unregularized Chi^2 = N_ELEMENTS(goodPixels)
;              value by DeltaChi^2 = sqrt(2*n_elements(goodPixels)).
;       The derived regularization corresponds to the maximum one still consistent
;       with the observations and the derived star formation history will be the
;       smoothest (minimum curvature) that is still consistent with the observations.
;     - For a detailed explanation see Section 18.5 of Press et al. (1992,
;       Numerical Recipes 2nd ed.) available here http://www.nrbook.com/a/bookfpdf.php.
;       The adopted implementation corresponds to their equation (18.5.10).
;   SKY: vector containing the spectrum of the sky to be included in the fit, or array
;       of dimensions SKY[nPixels,nSky] containing different sky spectra to add to
;       the model of the observed GALAXY spectrum. The SKY has to be log-rebinned as
;       the GALAXY spectrum and needs to have the same number of pixels.
;     - The sky is generally subtracted from the data before the PPXF fit. However,
;       for oservations very heavily dominated by the sky spectrum, where a very
;       accurate sky subtraction is critical, it may be useful *not* to subtract
;       the sky from the spectrum, but to include it in the fit using this keyword.
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
;   WEIGHTS: a named variable to receive the value of the weights by which each
;       template was multiplied to best fit the galaxy spectrum. The optimal
;       template can be computed with an array-vector multiplication:
;           TEMP = TEMPLATES # WEIGHTS (in IDL syntax)
;     - When the SKY keyword is used WEIGHTS[0:nTemplates-1] contains the weights
;       for the templates, while WEIGHTS[nTemplates:*] gives the ones for the sky.
;       In that case the best fitting galaxy template and sky are given by:
;           TEMP = TEMPLATES # WEIGHTS[0:nTemplates-1]
;           BESTSKY = SKY # WEIGHTS[nTemplates:*]
;     - When doing a two-sided fitting (see help for GALAXY parameter) *together*
;       with the SKY keyword, the sky weights are allowed to be different for the
;       left and right spectrum. In that case the output sky weights alternate
;       between the first (left) spectrum and the second (right) spectrum.
;
; OUTPUT PARAMETER:
;   SOL: seven elements vector containing in output the values of
;       [Vel,Sigma,h3,h4,h5,h6,Chi^2/DOF] of the best fitting solution, where DOF
;       is the number of Degrees of Freedom (number of fitted spectral pixels).
;       Vel is the velocity, Sigma is the velocity dispersion, h3-h6 are the
;       Gauss-Hermite coefficients. The model parameter are fitted simultaneously.
;     - I hardcoded the following safety limits on the fitting parameters:
;         a) Vel is constrained to be +/-2000 km/s from the first input guess
;         b) velScale/10 < Sigma < 1000 km/s
;         c) -0.3 < [h3,h4,...] < 0.3 (limits are extreme value for real galaxies)
;     - In the case of two-sided LOSVD fitting the output values refer
;       to the first input galaxy spectrum, while the second spectrum will
;       have by construction kinematics parameters [-Vel,Sigma,-h3,h4,-h5,h6].
;       If VSYST is nonzero (as required for two-sided fitting), then the
;       output velocity is measured with respect to VSIST.
;     - IMPORTANT: if Chi^2/DOF is not ~1 it means that the errors are not
;       properly estimated, or that the template is bad and it is *not* safe
;       to set the /CLEAN keyword.
;     - When MDEGREE > 1 then SOL contains in output the 7+MDEGREE elements
;       [Vel,Sigma,h3,h4,h5,h6,Chi^2/DOF,cx1,cx2,...,cxn], where cx1,cx2,...,cxn
;       are the coefficients of the multiplicative Legendre polynomials
;       of order 1,2,...,n. The polynomial can be explicitly evaluated as:
;           x = mdap_range(-1d,1d,n_elements(galaxy))
;           mpoly = 1d ; Multiplicative polynomial
;           for j=1,MDEGREE do mpoly += legendre(x,j)*sol[6+j]
;
;--------------------------------
; IMPORTANT: Proper usage of pPXF
;--------------------------------
;
; The PPXF routine can give sensible quick results with the default BIAS
; parameter, however, like in any penalized/filtered/regularized method, the
; optimal amount of penalization generally depends on the problem under study.
;
; The general rule here is that the penalty should leave the line-of-sight
; velocity-distribution (LOSVD) virtually unaffected, when it is well
; sampled and the signal-to-noise ratio (S/N) is sufficiently high.
;
; EXAMPLE: If you expect an LOSVD with up to a high h4 ~ 0.2 and your
; adopted penalty (BIAS) biases the solution towards a much lower h4 ~ 0.1,
; even when the measured sigma > 3*velScale and the S/N is high, then you
; are *misusing* the pPXF method!
;
; THE RECIPE: The following is a simple practical recipe for a sensible
; determination of the penalty in pPXF:
;
; 1. Choose a minimum (S/N)_min level for your kinematics extraction and
;    spatially bin your data so that there are no spectra below (S/N)_min;
;
; 2. Perform a fit of your kinematics *without* penalty (PPXF keyword BIAS=0).
;    The solution will be noisy and may be affected by spurious solutions,
;    however this step will allow you to check the expected mean ranges in
;    the Gauss-Hermite parameters [h3,h4] for the galaxy under study;
;
; 3. Perform a Monte Carlo simulation of your spectra, following e.g. the
;    included ppxf_simulation_example.pro routine. Adopt as S/N in the simulation 
;    the chosen value (S/N)_min and as input [h3,h4] the maximum representative
;    values measured in the non-penalized pPXF fit of the previous step;
;
; 4. Choose as penalty (BIAS) the *largest* value such that, for sigma > 3*velScale,
;    the mean difference delta between the output [h3,h4] and the input [h3,h4]
;    is well within (e.g. delta~rms/3) the rms scatter of the simulated values
;    (see e.g. Fig.2 of Emsellem et al. 2004, MNRAS, 352, 721).
;
;--------------------------------
;
; REQUIRED ROUTINES:
;       RANGE: by M. Cappellari (included in the PPXF distribution)
;       BVLS: by M. Cappellari from http://purl.org/cappellari/idl
;       MPFIT: by C.B. Markwardt http://purl.com/net/mpfit
;       ROBUST_SIGMA: by H. Freudenreich, which is included in
;           The IDL Astronomy User's Library http://idlastro.gfsc.nasa.gov/
;       CGPLOT: by David Fanning http://www.idlcoyote.com/, included in
;           The IDL Astronomy User's Library http://idlastro.gfsc.nasa.gov/
;
; MODIFICATION HISTORY:
;   V1.0 -- Created by Michele Cappellari, Leiden, 10 October 2001.
;   V3.47 -- First released version. MC, Leiden, 8 December 2003
;   V3.5 -- Included /OVERSAMPLE option. MC, Leiden, 11 December 2003
;   V3.6 -- Added MDEGREE option for multiplicative polynomials.
;           Linear implementation: fast, works well in most cases, but
;           can fail in certain cases. MC, Leiden, 19 March 2004
;   V3.7 -- Revised implementation of MDEGREE option. Nonlinear implementation:
;           straightforward, robust, but slower. MC, Leiden, 23 March 2004
;   V3.71 -- Updated documentation. MC, Leiden, 31 March 2004
;   V3.72 -- Corrected program stop after fit when MOMENTS=2.
;           Bug was introduced in V3.7. MC, Leiden, 28 April 2004
;   V3.73 -- Corrected bug: keyword ERROR was returned in pixels
;           instead of km/s. Decreased lower limit on fitted dispersion.
;           Thanks to Igor V. Chilingarian. MC, Leiden, 7 August 2004
;   V4.0 -- Introduced optional two-sided fitting assuming a reflection-symmetric
;           LOSVD for two input spectra. MC, Vicenza, 16 August 2004
;   V4.1 -- Corrected implementation of two-sided fitting of the LOSVD.
;           Thanks to Stefan van Dongen for reporting problems.
;           MC, Leiden, 3 September 2004
;   V4.11 -- Increased maximum number of iterations ITMAX in BVLS.
;           Thanks to Jesus Falcon-Barroso for reporting problems.
;           Introduced error message when velocity shift is too big.
;           Corrected output when MOMENTS=0. MC, Leiden, 21 September 2004
;   V4.12 -- Handle special case where a single template without additive
;           polynomials is fitted to the galaxy. MC, Leiden, 11 November 2004
;   V4.13 -- Updated documentation. MC, Vicenza, 30 December 2004
;   V4.14 -- Make sure input NOISE is a positive vector. MC, Leiden, 12 January 2005
;   V4.15 -- Verify that GOODPIXELS is monotonic and does not contain duplicated values.
;           After feedback from Richard McDermid. MC, Leiden, 10 February 2005
;   V4.16 -- Print number of nonzero templates. Do not print outliers in /QUIET mode.
;           MC, Leiden, 20 January 2006
;   V4.17 -- Updated documentation with important note on penalty determination.
;           MC, Oxford, 6 October 2007
;   V4.2 -- Introduced optional fitting of SKY spectrum. Many thanks to
;           Anne-Marie Weijmans for testing. MC, Oxford, 15 March 2008
;   V4.21 -- Use LA_LEAST_SQUARES (IDL 5.6) instead of SVDC when fitting
;           a single template. Please let me know if you need to use PPXF
;           with an older IDL version. MC, Oxford, 17 May 2008
;   V4.22 -- Added keyword POLYWEIGHTS. MC, Windhoek, 3 July 2008
;   V4.23 -- Corrected error message for too big velocity shift.
;           MC, Oxford, 27 November 2008
;   V4.3 -- Introduced REGUL keyword to perform linear regularization of WEIGHTS
;           in one or two dimensions. MC, Oxford, 4 Mach 2009
;   V4.4 -- Introduced Calzetti et al. (2000) PPXF_REDDENING_CURVE function to
;           estimate the reddening from the fit. MC, Oxford, 18 September 2009
;   V4.5 -- Dramatic speed up in the convolution of long spectra.
;           MC, Oxford, 13 April 2010
;   V4.6 -- Important fix to /CLEAN procedure: bad pixels are now properly
;           updated during the 3sigma iterations. MC, Oxford, 12 April 2011
;   V4.61 -- Use Coyote Graphics (http://www.idlcoyote.com/) by David W. Fanning.
;           The required routines are now included in NASA IDL Astronomy Library.
;           MC, Oxford, 29 July 2011
;   V4.62 -- Included option for 3D regularization and updated documentation of
;           REGUL keyword. MC, Oxford, 17 October 2011
;   V4.63 -- Do not change TEMPLATES array in output when REGUL is nonzero.
;           From feedback of Richard McDermid. MC, Oxford 25 October 2011
;   V4.64 -- Increased oversampling factor to 30x, when the /OVERSAMPLE keyword
;           is used. Updated corresponding documentation. Thanks to Nora 
;           Lu"tzgendorf for test cases illustrating errors in the recovered 
;           velocity when the sigma is severely undersampled.
;           MC, Oxford, 9 December 2011
;   V4.65 -- Expanded documentation of REGUL keyword. MC, Oxford, 15 November 2012
;-
;----------------------------------------------------------------------------
FUNCTION ppxf_reddening_curve, lambda, ebv
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
FUNCTION ppxf_BVLS_Solve, A, b, npoly,$
       external_library=external_library
compile_opt idl2, hidden

; No need to enforce positivity constraints if fitting one single template:
; use faster linear least-squares solution instead of BVLS.
;
s = size(a)
if s[0] eq 1 then $ ; A is a vector, not an array
    soluz = total(A*b)/total(A^2) $
else if s[2] eq npoly+1 then $ ; Fitting a single template
    soluz = la_least_squares(transpose(A),b) $
else begin               ; Fitting multiple templates
    bnd = dblarr(2,s[2],/NOZERO)
    mx = (machar()).xmax
    if npoly gt 0 then bnd[0,0:npoly-1] = -mx ; No bounds on Legendre polynomials
    bnd[0,npoly:*] = 0d  ; Positivity constraints on the templates (and sky spectra)
    bnd[1,*] = mx
    if external_library[0] eq 'none' then begin
        mdap_BVLS, A, B, bnd, soluz, ITMAX=15*s[2], IERR=ierr 
    endif else $
        mdap_bvls_external, A, B, bnd, soluz, external_library, IERR=ierr

    if ierr ne 0 then message, 'BVLS Error n. ' + strtrim(ierr,2)
endelse

return, soluz
END
;----------------------------------------------------------------------------
FUNCTION ppxf_fitfunc_optimal_template, pars, $
    BESTFIT=bestFit, BIAS=bias, CLEAN=clean, DEGREE=degree, $
    FACTOR=factor, GALAXY=galaxy, GOODPIXELS=goodPixels, MDEGREE=mdegree, $
    NOISE=noise, QUIET=quiet, SKY=sky, STAR=star, VSYST=vsyst, WEIGHTS=weights, $
    REGUL=regul, REG_DIM=reg_dim, LAMBDA=lambda,$
       external_library=external_library
compile_opt idl2, hidden

s = size(galaxy)
nspec = s[0]
npix = s[1]
npars = n_elements(pars) - mdegree*nspec  ; Parameters of the LOSVD only
if n_elements(lambda) gt 1 then npars = npars - 1 ; Fitting reddening

; pars = [vel,sigma,h3,h4,...,m1,m2,...]    ; Velocities are in pixels
;
dx = ceil(abs(vsyst)+abs(pars[0])+5d*pars[1]) ; Sample the Gaussian and GH at least to vsyst+vel+5*sigma
n = 2*dx*factor + 1
x = mdap_range(dx,-dx,n)   ; Evaluate the Gaussian using steps of 1/factor pixel
losvd = dblarr(n,nspec,/NOZERO)
for k=0,nspec-1 do begin    ; nspec=2 for two-sided fitting, otherwise nspec=1
    s = (k eq 0) ? 1d : -1d ; s=+1 for left spectrum, s=-1 for right one
    vel = vsyst + s*pars[0]
    w = (x - vel)/pars[1]
    w2 = w^2
    losvd[*,k] = exp(-0.5d*w2)/(sqrt(2d*!dpi)*pars[1]) ; Normalized total(Gaussian)=1

    ; Hermite polynomials normalized as in Appendix A of van der Marel & Franx (1993).
    ; These coefficients are given e.g. in Appendix C of Cappellari et al. (2002)
    ;
    if npars gt 2 then begin
        poly = 1d + s*pars[2]/Sqrt(3d)*(w*(2d*w2-3d)) $     ; H3
                  + pars[3]/Sqrt(24d)*(w2*(4d*w2-12d)+3d)   ; H4
        if npars eq 6 then $
            poly = poly + s*pars[4]/Sqrt(60d)*(w*(w2*(4d*w2-20d)+15d)) $    ; H5
                        + pars[5]/Sqrt(720d)*(w2*(w2*(8d*w2-60d)+90d)-15d)  ; H6
        losvd[*,k] = losvd[*,k]*poly
    endif
endfor

; The zeroth order multiplicative term is already included in the
; linear fit of the templates. The polynomial below has mean of 1.
;
x = mdap_range(-1d,1d,npix) ; X needs to be within [-1,1] for Legendre Polynomials
mpoly = 1d  ; The loop below can be null if mdegree < 1
for j=1,mdegree do $
    if nspec eq 2 then $ ; Different multiplicative poly for left and right spectra
        mpoly = mpoly + legendre(x,j) # pars[npars+2*j-[2,1]] $
    else mpoly = mpoly + legendre(x,j) * pars[npars+j-1]

; Multiplicative polynomials do not make sense when fitting reddening.
; In that case one has to assume the spectrum is well calibrated.
;
if n_elements(lambda) gt 1 then mpoly = ppxf_reddening_curve(lambda, pars[npars])

; Fill the columns of the design matrix of the least-squares problem
;
s = size(star)
ss = size(sky)
nsky = (ss[0] lt 2) ? ss[0] : ss[2] ; Number of sky spectra
ntemp = (s[0] eq 2) ? s[2] : 1     ; Number of template spectra
nrows = (degree + 1 + nsky)*nspec + ntemp
ncols = npix*nspec
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
c = dblarr(ncols,nrows)  ; This array is used for estimating predictions

for j=0,degree do $ ; Fill first columns of the Design Matrix
    if nspec eq 2 then begin
        leg = legendre(x,j)
        c[0,2*j] = [leg,leg*0d]   ; Additive polynomials for left spectrum
        c[0,2*j+1] = [leg*0d,leg] ; Additive polynomials for right spectrum
    endif else c[0,j] = legendre(x,j)

if factor gt 1 then pix = mdap_range(0d,s[1]-1d,s[1]*factor) ; Oversampled pixels range
tmp = dblarr(s[1],nspec,/NOZERO)
for j=0,ntemp-1 do begin
    if factor eq 1 then $ ; No oversampling of the template spectrum
        for k=0,nspec-1 do tmp[*,k] = mdap_ppxf_convol_fft(star[*,j],losvd[*,k]) $
    else begin             ; Oversample the template spectrum before convolution
        st = interpolate(star[*,j],pix,CUBIC=-0.5)   ; Sinc-like interpolation
        for k=0,nspec-1 do tmp[*,k] = rebin(mdap_ppxf_convol_fft(st,losvd[*,k]),s[1])
    endelse
    c[0,(degree+1)*nspec+j] = (mpoly*tmp[0:npix-1,*])[*] ; reform into a vector
endfor

; Add second-degree 1D, 2D or 3D linear regularization
; Press W.H., et al., 1992, Numerical Recipes, 2nd ed. equation (18.5.10)
;
if regul gt 0 then begin
    i = indgen(reg_dim) + (degree+1)*nspec
    p = npix*nspec
    diff = [-1d,2d,-1d]*regul
    dim = size(reg_dim,/DIM) > 1
    case dim of
        1: for j=1,reg_dim-2 do c[p++,i[j+[-1,0,1]]] = diff
        2: for k=0,reg_dim[1]-1 do $
               for j=0,reg_dim[0]-1 do begin
                   if j ne 0 && j ne reg_dim[0]-1 then c[p++,i[j+[-1,0,1],k]] = diff
                   if k ne 0 && k ne reg_dim[1]-1 then c[p++,i[j,k+[-1,0,1]]] = diff
               endfor
        3: for l=0,reg_dim[2]-1 do $
               for k=0,reg_dim[1]-1 do $
                   for j=0,reg_dim[0]-1 do begin
                       if j ne 0 && j ne reg_dim[0]-1 then c[p++,i[j+[-1,0,1],k,l]] = diff
                       if k ne 0 && k ne reg_dim[1]-1 then c[p++,i[j,k+[-1,0,1],l]] = diff
                       if l ne 0 && l ne reg_dim[2]-1 then c[p++,i[j,k,l+[-1,0,1]]] = diff
                   endfor
    endcase
endif

for j=0,nsky-1 do begin
    skyj = sky[*,j]
    k = (degree+1)*nspec + ntemp
    if nspec eq 2 then begin
        c[0,k+2*j] = [skyj,skyj*0d]   ; Sky for left spectrum
        c[0,k+2*j+1] = [skyj*0d,skyj] ; Sky for right spectrum
    endif else c[0,k+j] = skyj
endfor

a = c                     ; This array is used for the actual solution of the system
for j=0,nrows-1 do a[0,j] = c[0:npix*nspec-1,j]/noise ; Weight all columns with errors

if regul gt 0 then begin
    aa = a[[goodPixels,mdap_range(npix*nspec,ncols-1)],*]
    bb = [galaxy[goodPixels]/noise[goodPixels],replicate(0d,nreg)]
endif else begin
    aa = a[goodPixels,*]
    bb = galaxy[goodPixels]/noise[goodPixels]
endelse

; Select the spectral region to fit and solve the overconditioned system
; using SVD/BVLS. Use unweighted array for estimating bestfit predictions.
; Iterate to exclude pixels deviating more than 3*sigma if /CLEAN keyword is set.
;
npoly = (degree+1)*nspec ; Number of additive polynomials in the fit

repeat begin
    weights = ppxf_BVLS_Solve(aa,bb,npoly,$
       external_library=external_library)
    bestfit = c[0:npix*nspec-1,*] # weights
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
;
if npars gt 2 && bias ne 0 then $
    err = err + bias*robust_sigma(err, /ZERO)*sqrt(total(pars[2:*]^2))

return, err
END
;----------------------------------------------------------------------------
PRO mdap_ppxf, templates, galaxy, noise, velScale, start, sol, $
    BESTFIT=bestFit, BIAS=bias, CLEAN=clean, DEGREE=degree, ERROR=error, $
    GOODPIXELS=goodPixels, MDEGREE=mdegree, MOMENTS=moments, OVERSAMPLE=oversample, $
    POLYWEIGHTS=polyweights, PLOT=plot, QUIET=quiet, SKY=sky, VSYST=vsyst, WEIGHTS=weights, $
    REGUL=regul, LAMBDA=lambda, REDDENING=reddening,$
       fix_star_kin=fix_star_kin,range_v_star=range_v_star,range_s_star=range_s_star,$
       external_library=external_library
compile_opt idl2
on_error, 2

; Do extensive checking of possible input errors
;
s1 = size(templates)
if s1[0] ge 3 then begin
    reg_dim = s1[2:s1[0]]
    star = reform(templates,s1[1],product(reg_dim))
    s1 = size(star)
endif else begin
    reg_dim = s1[2]
    star = templates
endelse
if n_elements(regul) eq 0 then regul = 0
s2 = size(galaxy)
s3 = size(noise)
s4 = size(sky)
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
if s4[0] gt 0 && s4[1] ne s2[1] then message, 'SKY must have the same size as GALAXY'
if n_elements(degree) le 0 then degree = 4 else degree = degree > (-1)
if n_elements(mdegree) le 0 then mdegree = 0 else mdegree = mdegree > 0
if keyword_set(oversample) then factor = 30 else factor = 1
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
if moments ge 2 && n_elements(start) ne 2 then message, 'START must have two elements [V,sigma]'
if s2[0] eq 2 then goodPixels = [goodPixels,s2[1]-1+goodPixels]  ; two-sided fitting of LOSVD
if n_elements(vsyst) eq 0 then begin
    if s2[0] eq 2 then message, 'VSYST must be defined for two-sided fitting'
    vsyst = 0d
endif


if n_elements(external_library) eq 0 then external_library = 'none'

if moments gt 0 then begin

    ; Explicitly specify the step for the numerical derivatives (pixel/100)
    ; in MPFIT routine and force safety limits on the fitting parameters.
    ;
    if mdegree gt 0 then begin
        start1 = dblarr(moments+mdegree*s2[0])  ; Set [h3,h4,...] to zero as initial guess
        parinfo = replicate({step:1d-2,limits:[0d,0d],limited:[1,1],fixed:0}, moments+mdegree*s2[0])
        parinfo[moments:*].limits = [-1d,1d] ; force <100% corrections
        parinfo[moments:*].step = 1d-3
    endif else if n_elements(reddening) gt 0 then begin
        start1 = dblarr(moments+1)  ; Set [h3,h4,...] to zero as initial guess
        start1[moments] = reddening
        parinfo = replicate({step:1d-2,limits:[0d,0d],limited:[1,1],fixed:0}, moments+1)
        parinfo[moments].limits = [0d,10d] ; force positive E(B-V) < 10 mag
        parinfo[moments].step = 1d-3
    endif else begin
        start1 = dblarr(moments)  ; Set [h3,h4,...] to zero as initial guess
        parinfo = replicate({step:1d-2,limits:[0d,0d],limited:[1,1],fixed:0}, moments)
    endelse
    start1[0] = start/velScale  ; Convert velocity scale to pixels
    parinfo[0].limits = start1[0] + [-2d3,2d3]/velScale ; +/-2000 km/s from first guess
    parinfo[1].limits = [0.1d,1d3/velScale] ; hard-coded velScale/10<sigma<1000 km/s
    if n_elements(range_v_star) gt 0 then parinfo[0].limits = start1[0] + [range_v_star[0],range_v_star[1]]/velScale 
    if n_elements(range_s_star) gt 0 then parinfo[1].limits = [max([0.1d,start1[1] + range_s_star[0]/velscale]), start1[1] + range_s_star[1]/velScale] 

    if keyword_set(fix_star_kin) then parinfo[0:1].fixed = 1

    if moments gt 2 then begin
        parinfo[2:moments-1].limits = [-0.4d,0.4d] ; -0.4<[h3,h4,...]<0.4
        parinfo[2:moments-1].step = 1d-3
        if keyword_set(fix_star_kin) then parinfo[2:moments-1].fixed = 1
    endif
    if s1[1] le 2*(abs(vsyst/velScale)+abs(start1[0])+5d*start1[1]) then $
        message, 'Velocity shift too big: Adjust wavelength ranges of spectrum and templates'

;    print, start1
;    print, parinfo.limits
;    stop

    ; Here the actual calculation starts.
    ; If required, once the minimum is found, clean the pixels deviating
    ; more than 3*sigma from the best fit and repeat the minimization
    ; until the set of cleaned pixels does not change any more.
    ;
    good = goodPixels
    for j=0,4 do begin ; Do at most five cleaning iterations
        if n_elements(sky) eq 0 then $
            functArgs = {BIAS:bias, DEGREE:degree, FACTOR:factor, GALAXY:double(galaxy), $
                GOODPIXELS:goodPixels, MDEGREE:mdegree, NOISE:double(noise), $
                STAR:double(star), VSYST:vsyst/velScale, $
                REGUL:regul, REG_DIM:reg_dim, LAMBDA:lambda,$
                external_library:external_library} $
        else $
            functArgs = {BIAS:bias, DEGREE:degree, FACTOR:factor, GALAXY:double(galaxy), $
                GOODPIXELS:goodPixels, MDEGREE:mdegree, NOISE:double(noise), $
                SKY:double(sky), STAR:double(star), VSYST:vsyst/velScale, $
                REGUL:regul, REG_DIM:reg_dim, LAMBDA:lambda,$
                external_library:external_library}
        res = mpfit('ppxf_fitfunc_optimal_template', start1, ERRMSG=errmsg, $
            FTOL=1d-4, FUNCTARGS=functArgs, NFEV=ncalls, PARINFO=parinfo, $
            PERROR=error, /QUIET)
        if errmsg ne '' then message, errmsg
        if ~keyword_set(clean) then break
        goodOld = goodPixels
        goodPixels = good
        tmp = ppxf_fitfunc_optimal_template(res, BIAS=bias, /CLEAN, $
            DEGREE=degree, FACTOR=factor, GALAXY=galaxy, $
            GOODPIXELS=goodPixels, MDEGREE=mdegree, NOISE=noise, $
            QUIET=quiet, SKY=sky, STAR=star, VSYST=vsyst/velScale, $
            REGUL=regul, REG_DIM=reg_dim, LAMBDA=lambda,$
            external_library=external_library)
        if array_equal(goodOld,goodPixels) then break
    endfor

endif else begin

    if mdegree gt 0 then message, 'MDEGREE>0 not implemented with MOMENTS=0'
    res = start
    res[0:1] = res[0:1]/velScale  ; Convert velocity scale to pixels
    error = start*0d              ; No output error on parameters
    ncalls = 1

endelse

; Evaluate scatter at the bestfit (with BIAS=0)
; and also get the output BESTFIT and WEIGHTS.
;
err = ppxf_fitfunc_optimal_template(res, BESTFIT=bestFit, BIAS=0, DEGREE=degree, $
    FACTOR=factor, GALAXY=galaxy, GOODPIXELS=goodPixels, MDEGREE=mdegree, $
    NOISE=noise, SKY=sky, STAR=star, VSYST=vsyst/velScale, WEIGHTS=weights, $
    REGUL=regul, REG_DIM=reg_dim, LAMBDA=lambda,$
       external_library=external_library)
chi2 = robust_sigma(err, /ZERO)^2       ; Robust computation of Chi^2/DOF.
sol = dblarr(7+mdegree*s2[0])
sol[0] = res[0:(n_elements(start)>moments)-1]
sol[0] = sol[0:1]*velScale ; Bring velocity scale back to km/s
error[0:1] = error[0:1]*velScale ; Convert errors to km/s
sol[6] = chi2
if mdegree ge 1 then sol[7] = res[moments:*]
if n_elements(reddening) gt 0 then reddening = res[moments] ; Replace input with best fit
if degree ge 0 then polyweights = weights[0:(degree+1)*s2[0]-1] ; output weights for the additive polynomials
weights = weights[(degree+1)*s2[0]:*] ; output weights for the templates (or sky) only

; Print final results on the screen.
;
if ~keyword_set(quiet) then begin
    print, 'feval', 'V', 'sigma', 'h3', 'h4', 'h5', 'h6', 'Chi2/DOF', FORMAT='(8A10)'
    print, ncalls, sol[0:6], FORMAT='(i10,2f10.1,5f10.3)'
    nw = n_elements(weights)
    if n_elements(reddening) gt 0 then print, 'Reddening E(B-V): ', reddening, FORMAT='(a,g0.3)'
    print, 'Nonzero Templates: ', total(weights gt 0), ' / ', nw, FORMAT='(a,g0,a,g0)'
    if n_elements(weights) le 20 then begin
        print, 'Templates weights:'
        print, weights, FORMAT='(10g11.3)'
    endif
endif

; Plot final data-model comparison if required.
;
if keyword_set(plot) then begin
    mdap_basic_colors, black, white, red, green, blue, yellow, cyan, magenta, orange, mint, $
                       purple, pink, olive, lightblue, gray   
    screenr = GET_SCREEN_SIZE()
    window, xsize=screenr[0]*0.5, ysize=screenr[1]*0.5, retain=2
                       

    mn = min(bestfit[goodPixels], MAX=mx)
    resid = mn + galaxy - bestfit
    plot, galaxy, XTITLE='pixels', YTITLE='counts', /XSTYLE, /YNOZERO, $
        YRANGE=[min(resid[goodPixels]),mx], YSTYLE=2, XRANGE=[-0.02,1.02]*s2[1]*s2[0]
    oplot, bestfit, COLOR=red, THICK=2
    n = n_elements(goodPixels)
    oplot, goodPixels, replicate(mn,n*s2[0]), PSYM=3, COLOR=green;, /OVERPLOT
    oplot, goodPixels, resid[goodPixels], PSYM=4, COLOR=green, SYMSIZE=0.3;, /OVERPLOT
    w = where((goodPixels[1:*] - goodPixels) gt 1, m)
    for j=0,m-1 do begin
        x = mdap_range(goodPixels[w[j]],goodPixels[w[j]+1])
        oplot, x, resid[x], COLOR=blue;, /OVERPLOT
    endfor
    w = (m gt 0) ? [0,w,w+1,n-1] : [0,n-1]  ; Add first and last point
    for j=0,n_elements(w)-1 do $
        oplot, goodPixels[w[[j,j]]], [mn,bestfit[goodPixels[w[j]]]], COLOR=green;, /OVERPLOT
    wait, 0.01 ; Ensure screen refresh under windows
endif
;if keyword_set(plot) then begin
;    mn = min(bestfit[goodPixels], MAX=mx)
;    resid = mn + galaxy - bestfit
;    cgplot, galaxy, XTITLE='pixels', YTITLE='counts', /XSTYLE, /YNOZERO, $
;        YRANGE=[min(resid[goodPixels]),mx], YSTYLE=2, XRANGE=[-0.02,1.02]*s2[1]*s2[0]
;    cgplot, bestfit, COLOR='red', THICK=2, /OVERPLOT
;    n = n_elements(goodPixels)
;    cgplot, goodPixels, replicate(mn,n*s2[0]), PSYM=3, COLOR='forest green', /OVERPLOT
;    cgplot, goodPixels, resid[goodPixels], PSYM=4, COLOR='forest green', SYMSIZE=0.3, /OVERPLOT
;    w = where((goodPixels[1:*] - goodPixels) gt 1, m)
;    for j=0,m-1 do begin
;        x = mdap_range(goodPixels[w[j]],goodPixels[w[j]+1])
;        cgplot, x, resid[x], COLOR='blue', /OVERPLOT
;    endfor
;    w = (m gt 0) ? [0,w,w+1,n-1] : [0,n-1]  ; Add first and last point
;    for j=0,n_elements(w)-1 do $
;        cgplot, goodPixels[w[[j,j]]], [mn,bestfit[goodPixels[w[j]]]], COLOR='forest green', /OVERPLOT
;    wait, 0.01 ; Ensure screen refresh under windows
;endif

END
;----------------------------------------------------------------------------
