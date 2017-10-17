################################################################################
#
# Copyright (C) 2001-2017, Michele Cappellari
# E-mail: michele.cappellari_at_physics.ox.ac.uk
#
# Updated versions of the software are available from my web page
# http://purl.org/cappellari/software
#
# If you have found this software useful for your research,
# I would appreciate an acknowledgment to the use of the
# "Penalized Pixel-Fitting method by Cappellari & Emsellem (2004)
#  as upgraded in Cappellari (2017)".
#
# This software is provided as is without any warranty whatsoever.
# Permission to use, for non-commercial purposes is granted.
# Permission to modify for personal or internal use is granted,
# provided this copyright and disclaimer are included unchanged
# at the beginning of the file. All other rights are reserved.
#
################################################################################
#+
# NAME:
#   ppxf()
#
# PURPOSE:
#       Extract galaxy stellar kinematics (V, sigma, h3, h4, h5, h6,...)
#       or the stellar population and gas emission by fitting a template
#       to an observed spectrum in pixel space, using the
#       Penalized Pixel-Fitting (pPXF) method originally described in
#       Cappellari M., & Emsellem E., 2004, PASP, 116, 138
#           http://adsabs.harvard.edu/abs/2004PASP..116..138C
#       and upgraded in Cappellari M., 2017, MNRAS, 466, 798
#           http://adsabs.harvard.edu/abs/2017MNRAS.466..798C
#
#   The following key optional features are also available:
#   1)  An optimal template, positive linear combination of different input
#       templates, can be fitted together with the kinematics.
#   2)  One can enforce smoothness on the template weights during the fit. This
#       is useful to attach a physical meaning to the weights e.g. in terms of
#       the star formation history of a galaxy.
#   3)  One can fit multiple kinematic components for both the stars and the gas
#       emission lines. Both the stellar and gas LOSVD can be penalized and can
#       be described by a general Gauss-Hermite series.
#   4)  Any parameter of the LOSVD (e.g. sigma) for any kinematic component can
#       either be fitted, or held fixed to a given value, while other parameters
#       are fitted. Alternatively, parameters can be constrained to lie within
#       given limits or even tied by simple relations to other parameters.
#   5)  Additive and/or multiplicative polynomials can be included to adjust the
#       continuum shape of the template to the observed spectrum.
#   6)  Iterative sigma clipping can be used to clean the spectrum.
#   7)  It is possible to fit a mirror-symmetric LOSVD to two spectra at the
#       same time. This is useful for spectra taken at point-symmetric spatial
#       positions with respect to the center of an equilibrium stellar system.
#   8)  One can include sky spectra in the fit, to deal with cases where the sky
#       dominates the observed spectrum and an accurate sky subtraction is
#       critical.
#   9)  One can derive an estimate of the reddening in the spectrum.
#  10)  The covariance matrix can be input instead of the error spectrum, to
#       account for correlated errors in the spectral pixels.
#  11)  One can specify the weights fraction between two kinematics components,
#       e.g. to model bulge and disk contributions.
#  12)  One can use templates with higher resolution than the galaxy, to
#       improve the accuracy of the LOSVD extraction at low dispersion.
#
# CALLING SEQUENCE:
#
#   from ppxf import ppxf
#
#   pp = ppxf(templates, galaxy, noise, velscale, start,
#             bias=None, bounds=None, clean=False, component=0, degree=4,
#             fixed=None, fraction=None, gas_component=None, gas_names=None,
#             goodpixels=None, lam=None, linear=False, mask=None,
#             method='capfit', mdegree=0, moments=2, plot=False, quiet=False,
#             reddening=None, reg_ord=2, reg_dim=None, regul=0, sigma_diff=0,
#             sky=None, templates_rfft=None, tied=None, trig=False,
#             velscale_ratio=None, vsyst=0)
#
#   print(pp.sol)  # print best-fitting kinematics (V, sigma, h3, h4)
#   pp.plot()      # Plot best fit and gas lines
#
# INPUT PARAMETERS:
#   TEMPLATES: vector containing the spectrum of a single template star or more
#       commonly an array of dimensions TEMPLATES[nPixels, nTemplates]
#       containing different templates to be optimized during the fit of the
#       kinematics. nPixels has to be >= the number of galaxy pixels.
#     - To apply linear regularization to the WEIGHTS via the keyword REGUL,
#       TEMPLATES should be an array of two TEMPLATES[nPixels, nAge], three
#       TEMPLATES[nPixels, nAge, nMetal] or four
#       TEMPLATES[nPixels, nAge, nMetal, nAlpha] dimensions, depending on the
#       number of population variables one wants to study.
#       This can be useful to try to attach a physical meaning to the output
#       WEIGHTS, in term of the galaxy star formation history and chemical
#       composition distribution.
#       In that case the templates may represent single stellar population SSP
#       models and should be arranged in sequence of increasing age, metallicity
#       or alpha along the second, third or fourth dimension of the array
#       respectively.
#   GALAXY: vector containing the spectrum of the galaxy to be measured. The
#       star and the galaxy spectra have to be logarithmically rebinned but the
#       continuum should *not* be subtracted. The rebinning may be performed
#       with the LOG_REBIN routine that is distributed with PPXF.
#     - For high redshift galaxies, one should bring the spectra close to the
#       restframe wavelength, before doing the PPXF fit. This can be done by
#       dividing the observed wavelength by (1 + z), where z is a rough estimate
#       of the galaxy redshift, before the logarithmic rebinning.
#       See Section 2.4 of Cappellari (2017) for details.
#     - GALAXY can also be an array of dimensions GALAXY[nGalPixels, 2]
#       containing two spectra to be fitted, at the same time, with a
#       reflection-symmetric LOSVD. This is useful for spectra taken at
#       point-symmetric spatial positions with respect to the center of an
#       equilibrium stellar system.
#       For a discussion of the usefulness of this two-sided fitting see e.g.
#       Section 3.6 of Rix & White (1992, MNRAS, 254, 389).
#     - IMPORTANT: (1) For the two-sided fitting the VSYST keyword has to be
#       used. (2) Make sure the spectra are rescaled to be not too many order of
#       magnitude different from unity, to avoid over or underflow problems in
#       the calculation. E.g. units of erg/(s cm^2 A) may cause problems!
#   NOISE: vector containing the 1*sigma error (per pixel) in the galaxy
#       spectrum, or covariance matrix describing the correlated errors in the
#       galaxy spectrum. Of course this vector/matrix must have the same units
#       as the galaxy spectrum.
#     - If GALAXY is a Nx2 array, NOISE has to be an array with the same
#       dimensions.
#     - When NOISE has dimensions NxN it is assumed to contain the covariance
#       matrix with elements sigma(i, j). When the errors in the spectrum are
#       uncorrelated it is mathematically equivalent to input in PPXF an error
#       vector NOISE=errvec or a NxN diagonal matrix NOISE=np.diag(errvec**2)
#       (note squared!).
#     - IMPORTANT: the penalty term of the pPXF method is based on the
#       *relative* change of the fit residuals. For this reason the penalty will
#       work as expected even if no reliable estimate of the NOISE is available
#       (see Cappellari & Emsellem [2004] for details).
#       If no reliable noise is available this keyword can just be set to:
#           NOISE = np.ones_like(galaxy)  # Same weight for all pixels
#   VELSCALE: velocity scale of the spectra in km/s per pixel. It has to be the
#       same for both the galaxy and the template spectra.
#       An exception is when the VELSCALE_RATIO keyword is used, in which case
#       one can input TEMPLATES with smaller VELSCALE than GALAXY.
#     - VELSCALE is *defined* in pPXF by VELSCALE = c*Delta[np.log(lambda)],
#       which is approximately VELSCALE ~ c*Delta(lambda)/lambda.
#       See Section 2.3 of Cappellari (2017) for details.
#   START: Vector, or list/array of vectors [start1, start2, ...], with the 
#       initial estimate for the LOSVD parameters.
#     - When LOSVD parameters are not held fixed, each vector only needs to
#       contain START = [velStart, sigmaStart] the initial guess for the
#       velocity and the velocity dispersion in km/s. The starting values for
#       h3-h6 (if they are fitted) are all set to zero by default.
#       In other words, when MOMENTS=4,
#           START = [velStart, sigmaStart]
#       is interpreted as
#           START = [velStart, sigmaStart, 0, 0]
#     - When the LOSVD for some kinematic components is held fixed (see FIXED
#       keyword), all values for [Vel, Sigma, h3, h4,...] can be provided.
#     - Unless a good initial guess is available, it is recommended to set the
#       starting sigma >= 3*velscale in km/s (i.e. 3 pixels). In fact when the
#       LOSVD is severely undersampled, and far from the true solution, the
#       chi^2 of the fit becomes weakly sensitive to small variations in sigma
#       (see pPXF paper). In some instances the near-constancy of chi^2 may
#       cause premature convergence of the optimization.
#     - In the case of two-sided fitting a good starting value for the velocity
#       is velStart=0.0 (in this case VSYST will generally be nonzero).
#       Alternatively on should keep in mind that velStart refers to the first
#       input galaxy spectrum, while the second will have velocity -velStart.
#     - With multiple kinematic components START must be a list of starting
#       values, one for each different component.
#     - EXAMPLE: We want to fit two kinematic components. We fit 4 moments for
#       the first component and 2 moments for the second one as follows
#           component = [0, 0, ... 0, 1, 1, ... 1]
#           moments = [4, 2]
#           start = [[V1, sigma1], [V2, sigma2]]
#
# KEYWORDS:
#   BIAS: This parameter biases the (h3, h4, ...) measurements towards zero
#       (Gaussian LOSVD) unless their inclusion significantly decreases the
#       error in the fit. Set this to BIAS=0.0 not to bias the fit: the solution
#       (including [V, sigma]) will be noisier in that case. The default BIAS
#       should provide acceptable results in most cases, but it would be safe to
#       test it with Monte Carlo simulations. This keyword precisely corresponds
#       to the parameter \lambda in the Cappellari & Emsellem (2004) paper. Note
#       that the penalty depends on the *relative* change of the fit residuals,
#       so it is insensitive to proper scaling of the NOISE vector. A nonzero
#       BIAS can be safely used even without a reliable NOISE spectrum, or with
#       equal weighting for all pixels.
#   BOUNDS: Lower and upper bounds for every kinematic parameter.
#       This is an array, or list of arrays, with the same dimensions as START,
#       except for the last one, which is two. In practice, for every elements
#       of START one needs to specify a pair of values [lower, upper].
#     - EXAMPLE: We want to fit two kinematic components, with 4 moments for the
#       first component and 2 for the second (e.g. stars and gas). In this case
#           moments = [4, 2]
#           start_stars = [V1, sigma1, 0, 0]
#           start_gas = [V2, sigma2]
#           start = [start_stars, start_gas]
#       then we can specify boundaries for each kinematic parameter as
#           bounds_stars = [[V1_lo, V1_up], [sigma1_lo, sigma1_up], [-0.3, 0.3], [-0.3, 0.3]]
#           bounds_gas = [[V2_lo, V2_up], [sigma2_lo, sigma2_up]]
#           bounds = [bounds_stars, bounds_gas]
#   COMPONENT: When fitting more than one kinematic component, this keyword
#       should contain the component number of each input template. In principle
#       every template can belong to a different kinematic component.
#     - EXAMPLE: We want to fit the first 50 templates to component 0 and the
#       last 10 templates to component 1. In this case
#           component = [0]*50 + [1]*10
#       which, in Python syntax, is equivalent to
#           component = [0, 0, ... 0, 1, 1, ... 1]
#     - This keyword is especially useful when fitting both emission (gas) and
#       absorption (stars) templates simultaneously (see example for MOMENTS
#       keyword).
#   CLEAN: set this keyword to use the iterative sigma clipping method described
#       in Section 2.1 of Cappellari et al. (2002, ApJ, 578, 787).
#       This is useful to remove from the fit unmasked bad pixels, residual gas
#       emissions or cosmic rays.
#     - IMPORTANT: This is recommended *only* if a reliable estimate of the
#       NOISE spectrum is available. See also note below for .CHI2.
#   DEGREE: degree of the *additive* Legendre polynomial used to correct the
#       template continuum shape during the fit (default: 4).
#       Set DEGREE = -1 not to include any additive polynomial.
#   FIXED: Boolean specifying whether a given kinematic parameter has to be held
#       fixed with the value given in START.
#       This is an array, or list, with the same dimensions as START.
#     - EXAMPLE: We want to fit two kinematic components, with 4 moments for the
#       first component and 2 for the second. In this case
#           moments = [4, 2]
#           start = [[V1, sigma1, h3, h4], [V2, sigma2]]
#       then we can held fixed e.g. the sigma (only) of both components using
#           fixed = [[0, 1, 0, 0], [0, 1]]
#     - NOTE: Setting a negative MOMENTS for a kinematic component is entirely
#       equivalent to setting "fixed = 1" for all parameters of the given
#       kinematic component. In other words
#           moments = [-4, 2]
#       is equivalent to
#           moments = [4, 2]
#           fixed = [[1, 1, 1, 1], [0, 0]]
#   FRACTION: This keyword allows one to fix the ratio between the first two
#       kinematic components. This is a scalar defined as follows
#           FRACTION = np.sum(WEIGHTS[COMPONENT == 0]) \
#                    / np.sum(WEIGHTS[COMPONENT < 2])
#       This is useful e.g. to try to kinematically decompose bulge and disk.
#     - IMPORTANT: The TEMPLATES and GALAXY spectra should be normalized with
#       mean ~ 1 (as order of magnitude) for the FRACTION keyword to work as
#       expected. A warning is printed if this is not the case and the resulting
#       output FRACTION is inaccurate.
#     - The remaining kinematic components (COMPONENT > 1) are left free, and
#       this allows, for example, to still include gas emission line components.
#   GAS_COMPONENT: boolean vector, of the same size as COMPONENT, set to True
#       where the given COMPONENT describes a gas emission line. If given, pPXF
#       provides the pp.gas_flux and pp.gas_flux_error in output.
#     - EXAMPLE: If the first nstar components are stellar templates and the
#       rest are gas emission lines, one will set
#           gas_component = component >= nstar
#     - This keyword is also used to plot the gas lines with a different color.
#   GAS_NAMES: string array specifying the names of the emission lines (e.g.
#       gas_names=["Hbeta", "[OIII]",...], one per gas line. The length of
#       this vector must match the number of nonzero elements in GAS_COMPONENT.
#       This vector is only used to print the line names on the console.
#   GOODPIXELS: integer vector containing the indices of the good pixels in the
#       GALAXY spectrum (in increasing order). Only these pixels are included in
#       the fit.
#     - IMPORTANT: in all likely situations this keyword *has* to be specified.
#   LAM: When the keyword REDDENING is used, the user has to pass in this
#       keyword a vector with the same dimensions of GALAXY, giving the
#       restframe wavelength in Angstrom of every pixel in the input galaxy
#       spectrum. If one uses my LOG_REBIN routine to rebin the spectrum before
#       the PPXF fit:
#           from ppxf_util import log_rebin
#           specNew, logLam, velscale = log_rebin(lamRange, galaxy)
#       the wavelength can be obtained as lam = np.exp(logLam).
#     - When LAM is given, the wavelength is shown in the plot, instead of the
#       pixels.
#   LINEAR: set to True to keep *all* nonlinear parameters fixed and *only*
#       perform a linear fit for the templates and additive polynomials weights.
#       The output solution is a copy of the input one and the errors are zero.
#   MASK: Boolean vector of length GALAXY.size specifying with True the pixels
#       that should be included in the fit. This keyword is just an alternative
#       way of specifying the GOODPIXELS.
#   METHOD: {'capfit', 'trf', 'dogbox', 'lm'}, optional.
#       Algorithm to perform the non-linear minimization step. Default 'capfit'.
#       See documentation of scipy.optimize.least_squares for method != 'capfit'.
#   MDEGREE: degree of the *multiplicative* Legendre polynomial (with mean of 1)
#       used to correct the continuum shape during the fit (default: 0). The
#       zero degree multiplicative polynomial is always included in the fit as
#       it corresponds to the weights assigned to the templates.
#       Note that the computation time is longer with multiplicative polynomials
#       than with the same number of additive polynomials.
#     - IMPORTANT: Multiplicative polynomials cannot be used when the REDDENING
#       keyword is set.
#   MOMENTS: Order of the Gauss-Hermite moments to fit. Set this keyword to 4 to
#       fit [h3, h4] and to 6 to fit [h3, h4, h5, h6]. Note that in all cases
#       the G-H moments are fitted (non-linearly) *together* with [V, sigma].
#     - If MOMENTS=2 or MOMENTS is not set then only [V, sigma] are fitted and
#       the other parameters are returned as zero.
#     - If MOMENTS is negative then the kinematics of the given COMPONENT are
#       kept fixed to the input values.
#       NOTE: Setting a negative MOMENTS for a kinematic component is entirely
#       equivalent to setting `fixed = 1` for all parameters of the given
#       kinematic component.
#     - EXAMPLE: We want to keep fixed component 0, which has an LOSVD described
#       by [V, sigma, h3, h4] and is modeled with 100 spectral templates;
#       At the same time we fit [V, sigma] for COMPONENT=1, which is described
#       by 5 templates (this situation may arise when fitting stellar templates
#       with pre-determined stellar kinematics, while fitting the gas emission).
#       We should give in input to ppxf() the following parameters:
#           component = [0]*100 + [1]*5   # --> [0, 0, ... 0, 1, 1, 1, 1, 1]
#           moments = [-4, 2]
#           start = [[V, sigma, h3, h4], [V, sigma]]
#   VELSCALE_RATIO: Integer. Gives the integer ratio > 1 between the VELSCALE of
#       the GALAXY and the TEMPLATES. When this keyword is used, the templates
#       are convolved by the LOSVD at their native resolution, and only
#       subsequently are integrated over the pixels and fitted to GALAXY.
#       This is useful for accurate recovery of the LOSVD below VELSCALE when
#       templates with higher resolution than the galaxy spectrum are available.
#     - Note that in realistic situations the uncertainty in the knowledge and
#       variations of the intrinsic line-spread function become the limiting
#       factor in recovering the LOSVD well below VELSCALE.
#   PLOT: set this keyword to plot the best fitting solution and the residuals
#       at the end of the fit.
#     - One can also call the class function pp.plot() after the call to ppxf.
#   QUIET: set this keyword to suppress verbose output of the best fitting
#       parameters at the end of the fit.
#   REDDENING: Set this keyword to an initial estimate of the reddening
#       E(B-V) >= 0 to fit a positive reddening together with the kinematics and
#       the templates. The fit assumes the extinction curve of Calzetti et al.
#       (2000, ApJ, 533, 682) but any other prescriptions could be trivially
#       implemented by modifying the function REDDENING_CURVE below.
#     - IMPORTANT: The MDEGREE keyword cannot be used when REDDENING is set.
#   REGUL: If this keyword is nonzero, the program applies first or second order
#       linear regularization to the WEIGHTS during the PPXF fit.
#       Regularization is done in one, two or three dimensions depending on
#       whether the array of TEMPLATES has two, three or four dimensions
#       respectively.
#       Large REGUL values correspond to smoother WEIGHTS output. When this
#       keyword is nonzero the solution will be a trade-off between smoothness
#       of WEIGHTS and goodness of fit.
#     - When fitting multiple kinematic COMPONENT the regularization is applied
#       only to the first COMPONENT=0, while additional components are not
#       regularized. This is useful when fitting stellar population together
#       with gas emission lines. In that case the SSP spectral templates must be
#       given first and the gas emission templates are given last. In this
#       situation one has to use the REG_DIM keyword (below), to give PPXF the
#       dimensions of the population parameters (e.g. nAge, nMetal, nAlpha).
#       An usage example is given in ppxf_example_population_gas_sdss.py.
#     - The effect of the regularization scheme is the following:
#       With REG_ORD=1 it enforces the numerical first derivatives between
#       neighbouring weights (in the 1-dim case) to be equal to
#       `w[j] - w[j+1] = 0` with an error Delta=1/REGUL.
#       With REG_ORD=2 it enforces the numerical second derivatives between
#       neighboring weights (in the 1-dim case) to be equal to
#       `w[j-1] - 2*w[j] + w[j+1] = 0` with an error Delta=1/REGUL.
#       It may be helpful to define REGUL=1/Delta and view Delta as the
#       regularization error.
#     - IMPORTANT: Delta needs to be smaller but of the same order of magnitude
#       of the typical WEIGHTS to play an effect on the regularization. One
#       quick way to achieve this is:
#           (i) Divide the full TEMPLATES array by a scalar in such a way that
#               the typical template has a median of one:
#               TEMPLATES /= np.median(TEMPLATES);
#          (ii) Do the same for the input GALAXY spectrum:
#               GALAXY /= np.median(GALAXY).
#               In this situation a sensible guess for Delta will be a few
#               percent (e.g. 0.01 --> REGUL=100).
#     - Alternatively, for a more rigorous definition of the parameter REGUL:
#           (a) Perform an un-regularized fit (REGUL=0) and then rescale the
#               input NOISE spectrum so that
#               Chi^2/DOF = Chi^2/N_ELEMENTS(goodPixels) = 1.
#               This is achieved by rescaling the input NOISE spectrum as
#               NOISE = NOISE*sqrt(Chi**2/DOF) = NOISE*sqrt(pp.chi2);
#          (b) Increase REGUL and iteratively redo the pPXF fit until the Chi^2
#              increases from the unregularized value Chi^2 = len(goodPixels)
#              value by DeltaChi^2 = sqrt(2*len(goodPixels)).
#       The derived regularization corresponds to the maximum one still
#       consistent with the observations and the derived star formation history
#       will be the smoothest (minimum curvature or minimum variation) that is
#       still consistent with the observations.
#     - For a detailed explanation see Section 3.5 of Cappellari (2017).
#   REG_ORD: Order of the derivative that is minimized by the regularization.
#       The following two rotationally symmetric estimators are supported:
#       REG_ORD=1: minimizes integral of squared gradient: Grad[w] @ Grad[w].
#       REG_ORD=2: minimizes integral of squared curvature: Laplacian[w]**2.
#   REG_DIM: When using regularization with more than one kinematic component
#       (using the COMPONENT keyword), the regularization is only applied to the
#       first one (COMPONENT=0). This is useful to fit the stellar population
#       and gas emission together.
#       In this situation one has to use the REG_DIM keyword, to give PPXF the
#       dimensions of the population parameters (e.g. nAge, nMetal, nAlpha).
#       One should creates the initial array of population templates like
#       e.g. TEMPLATES[nPixels, nAge, nMetal, nAlpha] and define
#           reg_dim = TEMPLATES.shape[1:] = [nAge, nMetal, nAlpha]
#       The array of stellar templates is then reshaped into a 2-dim array as
#           TEMPLATES = TEMPLATES.reshape(TEMPLATES.shape[0], -1)
#       and the gas emission templates are appended as extra columns at the end.
#       An usage example is given in ppxf_example_population_gas_sdss.py.
#     - When using regularization with a single component (the COMPONENT keyword
#       is not used or contains identical values), the number of population
#       templates along different dimensions (e.g. nAge, nMetal, nAlpha) is
#       inferred from the dimensions of the TEMPLATES array and this keyword is
#       not necessary.
#   SIGMA_DIFF: Quadratic difference in km/s defined as
#           sigma_diff**2 = sigma_inst**2 - sigma_temp**2
#       between the instrumental dispersion of the galaxy spectrum and the
#       instrumental dispersion of the template spectra.
#       This keyword is useful when the templates have higher resolution than the
#       galaxy and they were not convolved to match the instrumental dispersion
#       of the galaxy spectrum. In this situation the convolution is done by
#       pPXF with increased accuracy, using an analytic Fourier Transform.
#   SKY: vector containing the spectrum of the sky to be included in the fit, or
#       array of dimensions SKY[nPixels, nSky] containing different sky spectra
#       to add to the model of the observed GALAXY spectrum. The SKY has to be
#       log-rebinned as the GALAXY spectrum and needs to have the same number of
#       pixels.
#     - The sky is generally subtracted from the data before the PPXF fit.
#       However, for observations very heavily dominated by the sky spectrum,
#       where a very accurate sky subtraction is critical, it may be useful
#       *not* to subtract the sky from the spectrum, but to include it in the
#       fit using this keyword.
#   TEMPLATES_RFFT: When calling pPXF many times with identical set of templates,
#       one can use this keyword to pass the real FFT of the templates, computed
#       in a previous pPXF call, stored in the pp.templates_rfft attribute.
#       This keyword mainly exists to demonstrate that there is no need for it...
#       IMPORTANT: Use this keyword only if you understand what you are doing!
#   TIED: A list of string expressions. Each expression "ties" the parameter to
#       other free or fixed parameters.  Any expression involving constants and
#       the parameter array p[j] are permitted. Since they are totally
#       constrained, tied parameters are considered to be fixed; no errors are
#       computed for them.
#       This is an array, or list of arrays, with the same dimensions as START.
#       In practice, for every elements of START one needs to specify either
#       an empty string '' implying that the parameter is free, or a string
#       expression involving some of the variables p[j], where the index j
#       represents the index of the flattened list of kinematic parameters.
#     - EXAMPLE: We want to fit three kinematic components, with 4 moments for
#       the first component and 2 moments for the second and third (e.g. stars
#       and two gas components). In this case
#           moments = [4, 2, 2]
#           start = [[V1, sigma1, 0, 0], [V2, sigma2], [V3, sigma3]]
#       then we can force the equality constraint V2 = V3 as follows
#           tied = [['', '', '', ''], ['', ''], ['p[4]', '']]
#       or we can force the equality constraint sigma2 = sigma3 as follows
#           tied = [['', '', '', ''], ['', ''], ['', 'p[5]']]
#     - NOTE: One could in principle use the `tied` keyword to completely tie
#       the LOSVD of two kinematic components. However this same effect is more
#       efficient achieved by assigning them to the same kinematic component
#       using the `component` keyword.
#   TRIG: Set `trig=True` to use trigonometric series as alternative to Legendre
#       polynomials, for both the additive and multiplicative polynomials.
#       When `trig=True` the fitted series below has N=degree/2 or N=mdegree/2
#           A_0 + \sum_{n=1}^{N} [A_n*cos(n*th) + B_n*sin(n*th)]
#     - IMPORTANT: The trigonometric series has periodic boundary conditions.
#       This is sometimes a desirable property, but this expansion is not as
#       flexible as the Legendre polynomials.
#   VSYST: galaxy systemic velocity (zero by default). The input initial guess
#       and the output velocities are measured with respect to this velocity.
#       The value assigned to this keyword is *crucial* for the two-sided
#       fitting. In this case VSYST can be determined from a previous normal
#       one-sided fit to the galaxy velocity profile. After that initial fit,
#       VSYST can be defined as the measured velocity at the galaxy center.
#       More accurately VSYST is the value which has to be subtracted to obtain
#       a nearly anti-symmetric velocity profile at the two opposite sides of
#       the galaxy nucleus.
#     - IMPORTANT: this value is generally *different* from the systemic
#       velocity one can get from the literature. Do not try to use that!
#
# OUTPUT PARAMETERS (stored as attributes of the PPXF class):
#   .BESTFIT: a named variable to receive a vector with the best fitting
#       template: this is a linear combination of the templates, convolved with
#       the best fitting LOSVD, multiplied by the multiplicative polynomials and
#       with subsequently added polynomial continuum terms.
#     - A version of this vector, *without* LOSVD convolution, is given by
#           BESTFIT = (TEMPLATES @ WEIGHTS)*mpoly + apoly,
#       where the expressions to evaluate mpoly and apoly are given in the
#       documentation of MPOLYWEIGHTS and POLYWEIGHTS respectively.
#   .GAS_BESTFIT: If `gas_component` is not None, this attribute returns the
#       best-fitting gas spectrum alone after the fit. The stellar spectrum
#       alone can be computed as stellar_spectrum = pp.bestfit - pp.gas_bestfit
#   .CHI2: The reduced chi^2 (=chi^2/DOF) of the fit.
#     - IMPORTANT: if Chi^2/DOF is not ~1 it means that the errors are not
#       properly estimated, or that the template is bad and it is *not* safe to
#       set the /CLEAN keyword.
#   .GAS_FLUX: Vector with the integrated flux of all lines set as True in the
#       input GAS_COMPONENT keyword. If a line is composed of a doublet, the
#       flux is that of both lines.
#     - IMPORTANT: pPXF makes no assumptions about the input flux units:
#       The returned .gas_flux has the same units as the value one would obtain
#       by just summing the values of the pixels of the gas emission.
#       This implies that, if the spectrum is in units of erg/(cm**2 s A), the
#       pPXF value should be multiplied by the pixel size in Anstrom at the line
#       wavelength to obtain the integrated line flux in units of erg/(cm**2 s).
#   .GAS_FLUX_ERROR: *formal* uncertainty (1*sigma) for the quantity pp.gas_flux.
#     - This error is approximate as it ignores the covariance between the gas
#       flux and any non-linear parameter. Bootstrapping can be used for more
#       accurate errors.
#     - These errors are meaningless unless Chi^2/DOF~1. However if one
#       *assumes* that the fit is good, a corrected estimate of the errors is:
#           gas_flux_error_corr = gas_flux_error*sqrt(chi^2/DOF)
#                               = pp.gas_flux_error*sqrt(pp.chi2).
#   .GOODPIXELS: integer vector containing the indices of the good pixels in the
#       fit. This vector is the same as the input GOODPIXELS if the CLEAN
#       keyword is *not* set, otherwise it will be updated by removing the
#       detected outliers.
#   .ERROR: this variable contain a vector of *formal* uncertainty (1*sigma) for
#       the fitted parameters in the output vector SOL. This option can be used
#       when speed is essential, to obtain an order of magnitude estimate of the
#       uncertainties, but we *strongly* recommend to run bootstrapping
#       simulations to obtain more reliable errors. In fact these errors can be
#       severely underestimated in the region where the penalty effect is most
#       important (sigma < 2*velscale).
#     - These errors are meaningless unless Chi^2/DOF~1. However if one
#       *assumes* that the fit is good, a corrected estimate of the errors is:
#           error_corr = error*sqrt(chi^2/DOF) = pp.error*sqrt(pp.chi2).
#     - IMPORTANT: when running Monte Carlo simulations to determine the error,
#       the penalty (BIAS) should be set to zero, or better to a very small
#       value. See Section 3.4 of Cappellari & Emsellem (2004) for an
#       explanation.
#   .POLYWEIGHTS: When DEGREE >= 0 contains the weights of the additive Legendre
#       polynomials of order 0, 1, ... DEGREE. The best fitting additive
#       polynomial can be explicitly evaluated as
#           from numpy.polynomial import legendre
#           x = np.linspace(-1, 1, len(galaxy))
#           apoly = legendre.legval(x, pp.polyweights)
#     - When `trig=True` the polynomial is evaluated as
#           apoly = pp.trigval(x, pp.polyweights)
#     - When doing a two-sided fitting (see help for GALAXY parameter), the
#       additive polynomials are allowed to be different for the left and right
#       spectrum. In that case the output weights of the additive polynomials
#       alternate between the first (left) spectrum and the second (right)
#       spectrum.
#   .MATRIX: prediction matrix[nPixels, DEGREE+nTemplates] of the linear system.
#     - pp.matrix[nPixels, :DEGREE] contains the additive polynomials if
#       DEGREE >= 0.
#     - pp.matrix[nPixels, DEGREE:] contains the templates convolved by the
#       LOSVD and multiplied by the multiplicative polynomials if MDEGREE > 0.
#     - pp.matrix[nPixels, -nGas:] contains the nGas emission line templates if
#       given. In the latter case the best fitting gas emission line spectrum is
#           lines = pp.matrix[:, -nGas:] @ pp.weights[-nGas:]
#     - The *formal* errors on pp.weights can be obtained as follows
#           import numpy as np
#           design_matrix = pp.matrix/pp.noise[:, None]
#           covariance_matrix = np.linalg.inv(design_matrix.T @ design_matrix)
#           weights_err = np.sqrt(np.diag(covariance_matrix)[DEGREE+1:])
#       Note however that the analytic errors computed in this way are only
#       meaningful when all pp.weights are nonzero in the pPXF solution.
#       This is generally not the case, but it may happen in special situations
#       (e.g. when fitting gas emission lines using a single stellar template).
#       For reliable errors Bootstrapping is generally recommended over this.
#   .MPOLYWEIGHTS: When MDEGREE > 0 this contains in output the coefficients of
#       the multiplicative Legendre polynomials of order 1, 2, ... MDEGREE.
#       The polynomial can be explicitly evaluated as:
#           from numpy.polynomial import legendre
#           x = np.linspace(-1, 1, len(galaxy))
#           mpoly = legendre.legval(x, np.append(1, pp.mpolyweights))
#     - When `trig=True` the polynomial is evaluated as
#           mpoly = pp.trigval(x, np.append(1, pp.mpolyweights))
#   .REDDENING: Best fitting E(B-V) value if the REDDENING keyword is set.
#   .SOL: Vector containing in output the parameters of the kinematics.
#       If MOMENTS=2 this contains [Vel, Sigma]
#       If MOMENTS=4 this contains [Vel, Sigma, h3, h4]
#       If MOMENTS=N this contains [Vel, Sigma, h3,... hN]
#     - When fitting multiple kinematic COMPONENT, pp.sol contains a list  with
#       the solution for all different components, one after the other, sorted
#       by COMPONENT: pp.sol = [sol1, sol2,...].
#     - Vel is the velocity, Sigma is the velocity dispersion, h3-h6 are the
#       Gauss-Hermite coefficients. The model parameters are fitted
#       simultaneously.
#     - IMPORTANT: The precise relation between the output pPXF velocity and
#       redshift is Vel = c*np.log(1 + z).
#       See Section 2.3 of Cappellari (2017) for a detailed explanation.
#     - These are the default safety limits on the fitting parameters:
#           a) Vel is constrained to be +/-2000 km/s from the first input guess
#           b) velscale/100 < Sigma < 1000 km/s
#           c) -0.3 < [h3, h4, ...] < 0.3 (limits are extreme value for real
#               galaxies)
#       They can be changed using the BOUNDS keyword.
#     - In the case of two-sided LOSVD fitting the output values refer to the
#       first input galaxy spectrum, while the second spectrum will have by
#       construction kinematics parameters [-Vel, Sigma, -h3, h4, -h5, h6].
#       If VSYST is nonzero (as required for two-sided fitting), then the output
#       velocity is measured with respect to VSIST.
#   .STATUS: Contains the output status of the optimization. Positive values
#       generally represent success (see scipy.optimize.least_squares
#       documentation).
#   .WEIGHTS: receives the value of the weights by which each template was
#       multiplied to best fit the galaxy spectrum. The optimal template can be
#       computed with an array-vector multiplication:
#           BESTEMP = TEMPLATES @ WEIGHTS
#     - See .MATRIX, for a discussion about the formal errors on .WEIGHTS.
#     - These weights do not include the weights of the additive polynomials
#       which are separately stored in pp.polyweights.
#     - When the SKY keyword is used WEIGHTS[:nTemplates] contains the weights
#       for the templates, while WEIGHTS[nTemplates:] gives the ones for the
#       sky. In that case the best fitting galaxy template and sky are given by:
#           BESTEMP = TEMPLATES @ WEIGHTS[:nTemplates]
#           BESTSKY = SKY @ WEIGHTS[nTemplates:]
#     - When doing a two-sided fitting (see help for GALAXY parameter)
#       *together* with the SKY keyword, the sky weights are allowed to be
#       different for the left and right spectrum. In that case the output sky
#       weights alternate between the first (left) spectrum and the second
#       (right) spectrum.
#
#--------------------------------
# IMPORTANT: Proper usage of pPXF
#--------------------------------
#
# The PPXF routine can give sensible quick results with the default BIAS
# parameter, however, like in any penalized/filtered/regularized method, the
# optimal amount of penalization generally depends on the problem under study.
#
# The general rule here is that the penalty should leave the line-of-sight
# velocity-distribution (LOSVD) virtually unaffected, when it is well sampled
# and the signal-to-noise ratio (S/N) is sufficiently high.
#
# EXAMPLE: If you expect an LOSVD with up to a high h4 ~ 0.2 and your adopted
# penalty (BIAS) biases the solution towards a much lower h4 ~ 0.1, even when
# the measured sigma > 3*velscale and the S/N is high, then you
# are *misusing* the pPXF method!
#
# THE RECIPE: The following is a simple practical recipe for a sensible
# determination of the penalty in pPXF:
#
# 1. Choose a minimum (S/N)_min level for your kinematics extraction and
#   spatially bin your data so that there are no spectra below (S/N)_min;
#
# 2. Perform a fit of your kinematics *without* penalty (PPXF keyword BIAS=0).
#   The solution will be noisy and may be affected by spurious solutions,
#   however this step will allow you to check the expected mean ranges in the
#   Gauss-Hermite parameters [h3, h4] for the galaxy under study;
#
# 3. Perform a Monte Carlo simulation of your spectra, following e.g. the
#   included ppxf_simulation_example.py routine. Adopt as S/N in the simulation
#   the chosen value (S/N)_min and as input [h3, h4] the maximum representative
#   values measured in the non-penalized pPXF fit of the previous step;
#
# 4. Choose as penalty (BIAS) the *largest* value such that, for
#   sigma > 3*velscale, the mean difference delta between the output [h3, h4]
#   and the input [h3, h4] is well within (e.g. delta~rms/3) the rms scatter of
#   the simulated values (see an example in Fig.2 of Emsellem et al. 2004,
#   MNRAS, 352, 721).
#
#--------------------------------
#
# REQUIRED ROUTINES:
#       CAPFIT: file capfit.py included in the distribution
#
# MODIFICATION HISTORY:
#   V1.0.0: Created by Michele Cappellari, Leiden, 10 October 2001.
#   V3.4.7: First released version. MC, Leiden, 8 December 2003
#   V3.5.0: Included /OVERSAMPLE option. MC, Leiden, 11 December 2003
#   V3.6.0: Added MDEGREE option for multiplicative polynomials.
#           Linear implementation: fast, works well in most cases, but can fail
#           in certain cases. MC, Leiden, 19 March 2004
#   V3.7.0: Revised implementation of MDEGREE option. Nonlinear
#           implementation: straightforward, robust, but slower.
#           MC, Leiden, 23 March 2004
#   V3.7.1: Updated documentation. MC, Leiden, 31 March 2004
#   V3.7.2: Corrected program stop after fit when MOMENTS=2. Bug was
#           introduced in V3.7.0. MC, Leiden, 28 April 2004
#   V3.7.3: Corrected bug: keyword ERROR was returned in pixels instead of
#           km/s. Decreased lower limit on fitted dispersion. Thanks to Igor
#           V. Chilingarian. MC, Leiden, 7 August 2004
#   V4.0.0: Introduced optional two-sided fitting assuming a reflection
#           symmetric LOSVD for two input spectra. MC, Vicenza, 16 August 2004
#   V4.1.0: Corrected implementation of two-sided fitting of the LOSVD.
#           Thanks to Stefan van Dongen for reporting problems.
#           MC, Leiden, 3 September 2004
#   V4.1.1: Increased maximum number of iterations ITMAX in BVLS.
#           Thanks to Jesus Falcon-Barroso for reporting problems.
#           Introduced error message when velocity shift is too big.
#           Corrected output when MOMENTS=0. MC, Leiden, 21 September 2004
#   V4.1.2: Handle special case where a single template without additive
#           polynomials is fitted to the galaxy. MC, Leiden, 11 November 2004
#   V4.1.3: Updated documentation. MC, Vicenza, 30 December 2004
#   V4.1.4: Make sure input NOISE is a positive vector.
#           MC, Leiden, 12 January 2005
#   V4.1.5: Verify that GOODPIXELS is monotonic and does not contain
#           duplicated values. After feedback from Richard McDermid.
#           MC, Leiden, 10 February 2005
#   V4.1.6: Print number of nonzero templates. Do not print outliers in
#           /QUIET mode. MC, Leiden, 20 January 2006
#   V4.1.7: Updated documentation with important note on penalty
#           determination. MC, Oxford, 6 October 2007
#   V4.2.0: Introduced optional fitting of SKY spectrum. Many thanks to
#           Anne-Marie Weijmans for testing. MC, Oxford, 15 March 2008
#   V4.2.1: Use LA_LEAST_SQUARES (IDL 5.6) instead of SVDC when fitting
#           a single template. Please let me know if you need to use PPXF
#           with an older IDL version. MC, Oxford, 17 May 2008
#   V4.2.2: Added keyword POLYWEIGHTS. MC, Windhoek, 3 July 2008
#   V4.2.3: Corrected error message for too big velocity shift.
#           MC, Oxford, 27 November 2008
#   V4.3.0: Introduced REGUL keyword to perform linear regularization of
#           WEIGHTS in one or two dimensions. MC, Oxford, 4 Mach 2009
#   V4.4.0: Introduced Calzetti et al. (2000) PPXF_REDDENING_CURVE function to
#           estimate the reddening from the fit. MC, Oxford, 18 September 2009
#   V4.5.0: Dramatic speed up in the convolution of long spectra.
#           MC, Oxford, 13 April 2010
#   V4.6.0: Important fix to /CLEAN procedure: bad pixels are now properly
#           updated during the 3sigma iterations. MC, Oxford, 12 April 2011
#   V4.6.1: Use Coyote Graphics (http://www.idlcoyote.com/) by David W.
#           Fanning. The required routines are now included in NASA IDL
#           Astronomy Library. MC, Oxford, 29 July 2011
#   V4.6.2: Included option for 3D regularization and updated documentation of
#           REGUL keyword. MC, Oxford, 17 October 2011
#   V4.6.3: Do not change TEMPLATES array in output when REGUL is nonzero.
#           From feedback of Richard McDermid. MC, Oxford 25 October 2011
#   V4.6.4: Increased oversampling factor to 30x, when the /OVERSAMPLE keyword
#           is used. Updated corresponding documentation. Thanks to Nora
#           Lu"tzgendorf for test cases illustrating errors in the recovered
#           velocity when the sigma is severely undersampled.
#           MC, Oxford, 9 December 2011
#   V4.6.5: Expanded documentation of REGUL keyword.
#           MC, Oxford, 15 November 2012
#   V4.6.6: Uses CAP_RANGE to avoid potential naming conflicts.
#           MC, Paranal, 8 November 2013
#   V5.0.0: Translated from IDL into Python and tested against the original
#           version. MC, Oxford, 6 December 2013
#   V5.0.1: Minor cleaning and corrections. MC, Oxford, 12 December 2013
#   V5.1.0: Allow for a different LOSVD for each template. Templates can be
#           stellar or can be gas emission lines. A PPXF version adapted for
#           multiple kinematic components existed for years. It was updated in
#           JAN/2012 for the paper by Johnston et al. (2013, MNRAS). This
#           version merges those changes with the public PPXF version, making
#           sure that all previous PPXF options are still supported.
#           MC, Oxford, 9 January 2014
#   V5.1.1: Fixed typo in the documentation of nnls_flags.
#           MC, Dallas Airport, 9 February 2014
#   V5.1.2: Replaced REBIN with INTERPOLATE with /OVERSAMPLE keyword. This is
#           to account for the fact that the Line Spread Function of the
#           observed galaxy spectrum already includes pixel convolution. Thanks
#           to Mike Blanton for the suggestion. MC, Oxford, 6 May 2014
#   V5.1.3: Allow for an input covariance matrix instead of an error spectrum.
#           MC, Oxford, 7 May 2014
#   V5.1.4: Support both Python 2.6/2.7 and Python 3.x. MC, Oxford, 25 May 2014
#   V5.1.5: Fixed deprecation warning. MC, Oxford, 21 June 2014
#   V5.1.6: Catch an additional input error. Updated documentation for Python.
#           Included templates `matrix` in output. Modified plotting colours.
#           MC, Oxford, 6 August 2014
#   V5.1.7: Relaxed requirement on input maximum velocity shift.
#           Minor reorganization of the code structure.
#           MC, Oxford, 3 September 2014
#   V5.1.8: Fixed program stop with `reddening` keyword. Thanks to Masatao
#           Onodera for reporting the problem. MC, Utah, 10 September 2014
#   V5.1.9: Pre-compute FFT and oversampling of templates. This speeds up the
#           calculation for very long or highly-oversampled spectra. Thanks to
#           Remco van den Bosch for reporting situations where this optimization
#           may be useful. MC, Las Vegas Airport, 13 September 2014
#   V5.1.10: Fixed bug in saving output introduced in previous version.
#           MC, Oxford, 14 October 2014
#   V5.1.11: Reverted change introduced in V5.1.2. Thanks to Nora Lu"tzgendorf
#           for reporting problems with oversample. MC, Sydney, 5 February 2015
#   V5.1.12: Use color= instead of c= to avoid new Matplotlib 1.4 bug.
#           MC, Oxford, 25 February 2015
#   V5.1.13: Updated documentation. MC, Oxford, 24 April 2015
#   V5.1.14: Fixed deprecation warning in Numpy 1.10.
#           MC, Oxford, 19 October 2015
#   V5.1.15: Updated documentation. Thanks to Peter Weilbacher for
#           corrections. MC, Oxford, 22 October 2015
#   V5.1.16: Fixed potentially misleading typo in documentation of MOMENTS.
#           MC, Oxford, 9 November 2015
#   V5.1.17: Expanded explanation of the relation between output velocity and
#           redshift. MC, Oxford, 21 January 2016
#   V5.1.18: Fixed deprecation warning in Numpy 1.11. Changed order from 1 to 3
#           during oversampling. Warn if sigma is under-sampled.
#           MC, Oxford, 20 April 2016
#   V5.2.0: Included `bounds`, `fixed` and `fraction` keywords.
#           MC, Baltimore, 26 April 2016
#   V5.3.0: Included `velscale_ratio` keyword to pass a set of templates with
#           higher resolution than the galaxy spectrum.
#           Changed `oversample` keyword to require integers not Booleans.
#           MC, Oxford, 9 May 2016
#   V5.3.1: Use wavelength in plot when available. Make plot() a class function.
#           Changes suggested and provided by Johann Cohen-Tanugi (LUPM).
#           MC, Oxford, 18 May 2016
#   V5.3.2: Backward compatibility change: allow `start` to be smaller than
#           `moments`. After feedback by Masato Onodera (NAOJ).
#           Updated documentation of `bounds` and `fixed`.
#           MC, Oxford, 22 May 2016
#   V5.3.3: Fixed Python 2 compatibility. Thanks to Masato Onodera (NAOJ).
#           MC, Oxford 24 May 2016
#   V6.0.0: Compute the Fourier Transform of the LOSVD analytically:
#           Major improvement in velocity accuracy when sigma < velscale.
#           Removed OVERSAMPLE keyword, which is now unnecessary.
#           Removed limit on velocity shift of templates. 
#           Simplified FFT zero padding. Updated documentation.
#           MC, Oxford, 28 July 2016
#   V6.0.1: Allow MOMENTS to be an arbitrary integer.
#           Allow for scalar MOMENTS with multiple kinematic components.
#           MC, Oxford, 10 August 2016
#   V6.0.2: Improved formatting of printed output. MC, Oxford, 15 August 2016
#   V6.0.3: Return usual Chi**2/DOF instead of Biweight estimate.
#           MC, Oxford, 1 December 2016
#   V6.0.4: Re-introduced `linear` keyword to only perform a linear fit and
#           skip the non-linear optimization. MC, Oxford, 30 January 2017
#   V6.0.5: Consistently use new _format_output() function both with/without
#           the `linear` keyword. Added .status attribute.
#           Changes suggested by Kyle B. Westfall (Santa Cruz).
#           MC, Oxford, 21 February 2017
#   V6.0.6: Added _linear_fit() and _nonlinear_fit() functions to better
#           clarify the code structure. Included `templates_rfft` keyword.
#           Updated documentation. Some code simplifications. 
#           MC, Oxford, 23 February 2017
#   V6.0.7: Use next_fast_len() for optimal rfft() zero padding.
#         - Included keyword `gas_component` in the .plot() method, to
#           distinguish gas emission lines in best-fitting plots.
#         - Improved plot of residuals for noisy spectra.
#         - Simplified regularization implementation.
#           MC, Oxford, 13 March 2017
#   V6.1.0: Introduced `trig` keyword to use a trigonometric series as
#           alternative to Legendre polynomials. MC, Oxford, 15 March 2017
#   V6.2.0: Improved curvature criterion for regularization when dim > 1.
#           MC, Oxford, 27 March 2017
#   V6.3.0: Included `reg_ord` keyword to allow for both first and second order
#           regularization. MC, Oxford, 30 March 2017
#   V6.3.1: Fixed program stop when fitting two galaxy spectra with
#           reflection-symmetric LOSVD. MC, Oxford, 13 April 2017
#   V6.3.2: Fixed possible program stop introduced in V6.0.7 and consequently 
#           removed unnecessary function _templates_rfft(). Many thanks to 
#           Jesus Falcon-Barroso for a very clear and useful bug report!
#           MC, Oxford, 4 May 2017
#   V6.4.0: Introduced `tied` keyword to tie parameters during fitting.
#           Included discussion of formal errors of .weights.
#           MC, Oxford, 12 May 2017
#   V6.4.1: _linear_fit() does not return unused status any more, for
#           consistency with the correspinding change to cap_mpfit.
#           MC, Oxford, 25 May 2017
#   V6.4.2: Fixed removal of bounds in solution, introduced in V6.4.1.
#           Thanks to Kyle B. Westfall (Santa Cruz) for reporting this.
#         - Included `method` keyword to use Scipy's least_squares()
#           as alternative to MPFIT. 
#         - Force float division in pixel conversion of `start` and `bounds`.
#           MC, Oxford, 2 June 2017
#   V6.5.0: Replaced MPFIT with CAPFIT for a Levenberg-Marquardt method with
#           fixed or tied variables, which rigorously accounts for box
#           constraints. MC, Oxford, 23 June 2017
#   V6.6.0: Print and return gas fluxes and errors, if requested, with the new
#           `gas_component` and `gas_names` keywords. MC, Oxford, 27 June 2017
#   V6.6.1: Included note on .gas_flux output units. Thanks to Xihan Ji
#           (Tsinghua University) for the feedback. MC, Oxford, 4 August 2017
#   V6.6.2: Fixed program stop with a 2-dim templates array and regularization.
#           Thanks to Adriano Poci (Macquarie University) for the clear report
#           and the fix. MC, Oxford, 15 September 2017
#   V6.6.3: Reduced bounds on multiplicative polynomials and clipped to positive
#           values. Thanks to Xihan Ji (Tsinghua University) for providing an
#           example of slightly negative gas emission lines, when the spectrum
#           contains essentially just noise.
#         - Improved visualization of masked pixels.
#           MC, Oxford, 25 September 2017
#   V6.6.4: Check for NaN in `galaxy` and check all `bounds` have two elements.
#           Allow `start` to be either a list or an array or vectors.
#           MC, Oxford, 5 October 2017
#   V6.6.5: Raise an error if any template is identically zero in fitted range.
#           This can happen if a gas line entirely falls within a masked region.
#         - Corrected `gas_flux_error` normalization, when input not normalized.
#         - Return .gas_bestfit attribute when gas_component is not None.
#         - Fixed program stop with `linear` keyword.
#           MC, Oxford, 16 October 2017
#
################################################################################

from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial import legendre, hermite
from scipy import optimize, linalg, misc, fftpack

#import capfit
from . import capfit

################################################################################

def trigvander(x, deg):
    """
    Analogue to legendre.legvander(), but for a trigonometric
    series rather than Legendre polynomials:

    `deg` must be an even integer

    """
    u = np.pi*x[:, None]   # [-pi, pi] interval
    j = np.arange(1, deg//2 + 1)
    mat = np.ones((x.size, deg + 1))
    mat[:, 1:] = np.hstack([np.cos(j*u), np.sin(j*u)])

    return mat

################################################################################

def trigval(x, c):
    """
    Analogue to legendre.legval(), but for a trigonometric
    series rather than Legendre polynomials:

    Evaluate a trigonometric series with coefficients `c` at points `x`.

    """
    return trigvander(x, c.size - 1).dot(c)

################################################################################

def nnls_flags(A, b, npoly):
    """
    Solves min||A*x - b|| with
    x[j] >= 0 for j >= npoly
    x[j] free for j < npoly
    where A[m, n], b[m], x[n], flag[n]

    """
    m, n = A.shape
    AA = np.hstack([A, -A[:, :npoly]])
    x = optimize.nnls(AA, b)[0]
    x[:npoly] -= x[n:]

    return x[:n]

################################################################################

def rebin(x, factor):
    """
    Rebin a one-dimensional vector by averaging
    in groups of "factor" adjacent values

    """
    return np.mean(x.reshape(-1, factor), axis=1)

################################################################################

def robust_sigma(y, zero=False):
    """
    Biweight estimate of the scale (standard deviation).
    Implements the approach described in
    "Understanding Robust and Exploratory Data Analysis"
    Hoaglin, Mosteller, Tukey ed., 1983, Chapter 12B, pg. 417

    """
    y = np.ravel(y)
    d = y if zero else y - np.median(y)

    mad = np.median(np.abs(d))
    u2 = (d/(9.0*mad))**2  # c = 9
    good = u2 < 1.0
    u1 = 1.0 - u2[good]
    num = y.size * ((d[good]*u1**2)**2).sum()
    den = (u1*(1.0 - 5.0*u2[good])).sum()
    sigma = np.sqrt(num/(den*(den - 1.0)))  # see note in above reference

    return sigma

################################################################################

def reddening_curve(lam1, ebv):
    """
    Reddening curve of Calzetti et al. (2000, ApJ, 533, 682; here C+00).
    This is reliable between 0.12 and 2.2 micrometres.
    - LAMBDA is the restframe wavelength in Angstrom of each pixel in the
      input galaxy spectrum (1 Angstrom = 1e-4 micrometres)
    - EBV is the assumed E(B-V) colour excess to redden the spectrum.
      In output the vector FRAC gives the fraction by which the flux at each
      wavelength has to be multiplied, to model the dust reddening effect.

    """
    lam = 1e4/lam1  # Convert Angstrom to micrometres and take 1/lambda
    rv = 4.05  # C+00 equation (5)

    # C+00 equation (3) but extrapolate for lam > 2.2
    # C+00 equation (4) but extrapolate for lam < 0.12
    k1 = np.where(lam >= 6300,
                  rv + 2.659*(1.040*lam - 1.857),
                  rv + 2.659*(1.509*lam - 0.198*lam**2 + 0.011*lam**3 - 2.156))
    fact = 10**(-0.4*ebv*k1.clip(0))  # Calzetti+00 equation (2) with opposite sign

    return fact # The model spectrum has to be multiplied by this vector

################################################################################

def _bvls_solve(A, b, npoly):

    # No need to enforce positivity constraints if fitting one single template:
    # use faster linear least-squares solution instead of NNLS.
    m, n = A.shape
    if m == 1:                  # A is a vector, not an array
        soluz = A.dot(b)/A.dot(A)
    elif n == npoly + 1:        # Fitting a single template
        soluz = linalg.lstsq(A, b)[0]
    else:                       # Fitting multiple templates
        soluz = nnls_flags(A, b, npoly)

    return soluz

################################################################################

def _losvd_rfft(pars, nspec, moments, nl, ncomp, vsyst, factor, sigma_diff):
    """
    Analytic Fourier Transform (of real input) of the Gauss-Hermite LOSVD.
    Equation (38) of Cappellari M., 2017, MNRAS, 466, 798
    http://adsabs.harvard.edu/abs/2017MNRAS.466..798C

    """
    losvd_rfft = np.empty((nl, ncomp, nspec), dtype=complex)
    p = 0
    for j, mom in enumerate(moments):  # loop over kinematic components
        for k in range(nspec):  # nspec=2 for two-sided fitting, otherwise nspec=1
            s = 1 if k == 0 else -1  # s=+1 for left spectrum, s=-1 for right one
            vel, sig = vsyst + s*pars[0 + p], pars[1 + p]
            a, b = [vel, sigma_diff]/sig
            w = np.linspace(0, np.pi*factor*sig, nl)
            losvd_rfft[:, j, k] = np.exp(1j*a*w - 0.5*(1 + b**2)*w**2)

            if mom > 2:
                n = np.arange(3, mom + 1)
                nrm = np.sqrt(misc.factorial(n)*2**n)   # vdMF93 Normalization
                coeff = np.append([1, 0, 0], (s*1j)**n * pars[p - 1 + n]/nrm)
                poly = hermite.hermval(w, coeff)
                losvd_rfft[:, j, k] *= poly
        p += mom

    return np.conj(losvd_rfft)

################################################################################

def _regularization(a, npoly, npix, nspec, reg_dim, reg_ord, regul):
    """
    Add first or second order 1D, 2D or 3D linear regularization.
    Equation (25) of Cappellari M., 2017, MNRAS, 466, 798
    http://adsabs.harvard.edu/abs/2017MNRAS.466..798C

    """
    b = a[:, npoly : npoly + np.prod(reg_dim)].reshape(-1, *reg_dim)
    p = npix*nspec

    if reg_ord == 1:   # Minimize integral of (Grad[w] @ Grad[w])
        diff = np.array([1, -1])*regul
        if reg_dim.size == 1:
            for j in range(reg_dim[0] - 1):
                b[p, j : j + 2] = diff
                p += 1
        elif reg_dim.size == 2:
            for k in range(reg_dim[1]):
                for j in range(reg_dim[0]):
                    if j < reg_dim[0] - 1:
                        b[p, j : j + 2, k] = diff
                        p += 1
                    if k < reg_dim[1] - 1:
                        b[p, j, k : k + 2] = diff
                        p += 1
        elif reg_dim.size == 3:
            for q in range(reg_dim[2]):
                for k in range(reg_dim[1]):
                    for j in range(reg_dim[0]):
                        if j < reg_dim[0] - 1:
                            b[p, j : j + 2, k, q] = diff
                            p += 1
                        if k < reg_dim[1] - 1:
                            b[p, j, k : k + 2, q] = diff
                            p += 1
                        if q < reg_dim[2] - 1:
                            b[p, j, k, q : q + 2] = diff
                            p += 1
    elif reg_ord == 2:   # Minimize integral of Laplacian[w]**2
        diff = np.array([1, -2, 1])*regul
        if reg_dim.size == 1:
            for j in range(1, reg_dim[0] - 1):
                b[p, j - 1 : j + 2] = diff
                p += 1
        elif reg_dim.size == 2:
            for k in range(reg_dim[1]):
                for j in range(reg_dim[0]):
                    if 0 < j < reg_dim[0] - 1:
                        b[p, j - 1 : j + 2, k] = diff
                    if 0 < k < reg_dim[1] - 1:
                        b[p, j, k - 1 : k + 2] += diff
                    p += 1
        elif reg_dim.size == 3:
            for q in range(reg_dim[2]):
                for k in range(reg_dim[1]):
                    for j in range(reg_dim[0]):
                        if 0 < j < reg_dim[0] - 1:
                            b[p, j - 1 : j + 2, k, q] = diff
                        if 0 < k < reg_dim[1] - 1:
                            b[p, j, k - 1 : k + 2, q] += diff
                        if 0 < q < reg_dim[2] - 1:
                            b[p, j, k, q - 1 : q + 2] += diff
                        p += 1

################################################################################

class ppxf(object):

    def __init__(self, templates, galaxy, noise, velscale, start,
                 bias=None, bounds=None, clean=False, component=0, degree=4,
                 fixed=None, fraction=None, gas_component=None, gas_names=None,
                 goodpixels=None, lam=None, linear=False, mask=None,
                 method='capfit', mdegree=0, moments=2, plot=False, quiet=False,
                 reddening=None, reg_ord=2, reg_dim=None, regul=0, sigma_diff=0,
                 sky=None, templates_rfft=None, tied=None, trig=False,
                 velscale_ratio=None, vsyst=0):

        # Do extensive checking of possible input errors
        #
        self.galaxy = galaxy
        self.nspec = galaxy.ndim     # nspec=2 for reflection-symmetric LOSVD
        self.npix = galaxy.shape[0]  # total pixels in the galaxy spectrum
        self.noise = noise
        self.clean = clean
        self.fraction = fraction
        self.degree = max(degree, -1)
        self.mdegree = max(mdegree, 0)
        self.method = method
        self.quiet = quiet
        self.sky = sky
        self.vsyst = vsyst/velscale
        self.regul = regul
        self.lam = lam
        self.nfev = 0
        self.reddening = reddening
        self.reg_dim = np.asarray(reg_dim)
        self.reg_ord = reg_ord
        self.star = templates.reshape(templates.shape[0], -1)
        self.npix_temp, self.ntemp = self.star.shape
        self.factor = 1   # default value
        self.sigma_diff = sigma_diff/velscale
        self.status = 0   # Initialize status as failed
        self.velscale = velscale

        if method != 'capfit':
            assert method in ['trf', 'dogbox', 'lm'], \
                "`method` must be 'capfit', 'trf', 'dogbox' or 'lm'"
            assert tied is None, "Parameters can only be tied with method='capfit'"
            assert fixed is None, "Parameters can only be fixed with method='capfit'"
            if method == 'lm':
                assert bounds is None, "Bounds not supported with method='lm'"

        if trig:
            assert degree < 0 or  degree % 2 == 0, \
                "`degree` must be even with trig=True"
            assert mdegree < 0 or mdegree % 2 == 0, \
                "`mdegree` must be even with trig=True"
            self.polyval = trigval
            self.polyvander = trigvander
        else:
            self.polyval = legendre.legval
            self.polyvander = legendre.legvander

        if velscale_ratio is not None:
            assert isinstance(velscale_ratio, int), \
                "VELSCALE_RATIO must be an integer"
            self.npix_temp -= self.npix_temp % velscale_ratio
            # Make size multiple of velscale_ratio
            self.star = self.star[:self.npix_temp, :]
            # This is the size after rebin()
            self.npix_temp //= velscale_ratio
            self.factor = velscale_ratio

        component = np.atleast_1d(component)
        assert component.dtype == int, "COMPONENT must be integers"

        if component.size == 1 and self.ntemp > 1:  # component is a scalar
            # all templates have the same LOSVD
            self.component = np.zeros(self.ntemp, dtype=int)
        else:
            assert component.size == self.ntemp, \
                "There must be one kinematic COMPONENT per template"
            self.component = component

        if gas_component is None:
            self.gas_component = None
        else:
            self.gas_component = np.asarray(gas_component)
            assert self.gas_component.dtype == bool, \
                "`gas_component` must be boolean"
            assert self.gas_component.size == component.size, \
                "`gas_component` and `component` must have the same size"
            if gas_names is None:
                self.gas_names = np.full(np.sum(gas_component), 'Unknown')
            else:
                assert gas_component.sum() == len(gas_names), \
                    "There must be one name per gas emission line template"
                self.gas_names = gas_names

        tmp = np.unique(component)
        self.ncomp = tmp.size
        assert np.array_equal(tmp, np.arange(self.ncomp)), \
            "COMPONENT must range from 0 to NCOMP-1"

        if fraction is not None:
            assert 0 < fraction < 1, "Must be `0 < fraction < 1`"
            assert self.ncomp >= 2, \
                "At least 2 COMPONENTs are needed with FRACTION keyword"

        if regul > 0 and reg_dim is None:
            assert self.ncomp == 1, \
                "REG_DIM must be specified with more than one component"
            self.reg_dim = np.asarray(templates.shape[1:])

        moments = np.atleast_1d(moments)
        if moments.size == 1:
            # moments is scalar: all LOSVDs have same number of G-H moments
            moments = np.full(self.ncomp, moments, dtype=int)

        self.fixall = moments < 0  # negative moments --> keep entire LOSVD fixed
        self.moments = np.abs(moments)

        assert tmp.size == self.moments.size, \
            "MOMENTS must be an array of length NCOMP"

        if regul is None:
            self.regul = 0

        if sky is not None:
            assert sky.shape[0] == galaxy.shape[0], \
                "GALAXY and SKY must have the same size"

        assert galaxy.ndim < 3 and noise.ndim < 3, \
            "Wrong GALAXY or NOISE input dimensions"

        if noise.ndim == 2 and noise.shape[0] == noise.shape[1]:
            # NOISE is a 2-dim covariance matrix
            assert noise.shape[0] == galaxy.shape[0], \
                "Covariance Matrix must have size xpix*npix"
            # Cholesky factor of symmetric, positive-definite covariance matrix
            noise = linalg.cholesky(noise, lower=1)
            # Invert Cholesky factor
            self.noise = linalg.solve_triangular(noise, np.identity(noise.shape[0]), lower=1)
        else:   # NOISE is an error spectrum
            assert galaxy.shape == noise.shape, \
                "GALAXY and NOISE must have the same size"
            assert np.all((noise > 0) & np.isfinite(noise)), \
                "NOISE must be a positive vector"
            if self.nspec == 2:   # reflection-symmetric LOSVD
                self.noise = self.noise.T.reshape(-1)
                self.galaxy = self.galaxy.T.reshape(-1)

        assert np.all(np.isfinite(galaxy)), 'GALAXY must be finite'

        assert self.npix_temp >= galaxy.shape[0], \
            "TEMPLATES length cannot be smaller than GALAXY"

        if reddening is not None:
            assert lam is not None, "LAM must be given with REDDENING keyword"
            assert mdegree < 1, "MDEGREE cannot be used with REDDENING keyword"

        if lam is not None:
            assert lam.shape == galaxy.shape, \
                "GALAXY and LAM must have the same size"

        if mask is not None:
            assert mask.dtype == bool, "MASK must be a boolean vector"
            assert mask.shape == galaxy.shape, \
                "GALAXY and MASK must have the same size"
            assert goodpixels is None, \
                "GOODPIXELS and MASK cannot be used together"
            goodpixels = np.flatnonzero(mask)

        if goodpixels is None:
            self.goodpixels = np.arange(galaxy.shape[0])
        else:
            assert np.all(np.diff(goodpixels) > 0), \
                "GOODPIXELS is not monotonic or contains duplicated values"
            assert goodpixels[0] >= 0 and goodpixels[-1] < galaxy.shape[0], \
                "GOODPIXELS are outside the data range"
            self.goodpixels = goodpixels

        m1 = np.max(np.abs(self.star), 0)
        m2 = np.max(np.abs(self.star[self.goodpixels, :]), 0)
        assert np.all(m2 > m1/1e3), \
            "TEMPLATES cannot be identically zero in fitted range"

        if bias is None:
            # Cappellari & Emsellem (2004) pg.144 left
            self.bias = 0.7*np.sqrt(500./self.goodpixels.size)
        else:
            self.bias = bias

        if self.ncomp == 1:
            start1 = [start]
        else:
            assert hasattr(start, "__len__"), \
                "START must be a list/array of vectors [start1, start2,...]"
            assert len(start) == self.ncomp, \
                "There must be one START per COMPONENT"
            start1 = list(start)  # Make a copy in both Python 2 and 3

        # Pad with zeros when `start[j]` has fewer elements than `moments[j]`
        for j, (st, mo) in enumerate(zip(start1, self.moments)):
            st = np.asarray(st, dtype=float)   # Make sure starting guess is float
            start1[j] = np.pad(st, (0, mo - len(st)), 'constant')

        if bounds is not None:
            if self.ncomp == 1:
                bounds = [bounds]
            assert list(map(len, bounds)) == list(map(len, start1)), \
                "BOUNDS and START must have the same shape"
            assert np.all([hasattr(c, "__len__") and len(c) == 2
                           for a in bounds for c in a]), \
                "All BOUNDS must have two elements [lb, ub]"

        if fixed is not None:
            if self.ncomp == 1:
                fixed = [fixed]
            assert list(map(len, fixed)) == list(map(len, start1)), \
                "FIXED and START must have the same shape"

        if tied is not None:
            if self.ncomp == 1:
                tied = [tied]
            assert list(map(len, tied)) == list(map(len, start1)), \
                "TIED and START must have the same shape"

        if galaxy.ndim == 2:
            # two-sided fitting of LOSVD
            assert vsyst != 0, "VSYST must be defined for two-sided fitting"
            self.goodpixels = np.append(self.goodpixels, galaxy.shape[0] + self.goodpixels)

        self.npad = fftpack.next_fast_len(self.star.shape[0])
        if templates_rfft is None:
            # Pre-compute FFT of real input of all templates
            self.templates_rfft = np.fft.rfft(self.star, self.npad, axis=0)
        else:
            self.templates_rfft = templates_rfft

        # Convert velocity from km/s to pixels
        for j, s in enumerate(start1):
            start1[j][:2] = s[:2]/velscale

        if linear:
            assert mdegree <= 0, "Must be `mdegree` <= 0 with `linear`=True"
            if reddening is None:
                params = np.concatenate(start1)   # Flatten list
            else:
                params = np.append(start1, reddening)
            perror = np.zeros_like(params)
            self.method = 'linear'
            self.status = 1   # Status irrelevant for linear fit
            self.njev = 0     # Jacobian is not evaluated
        else:
            params, perror = self._nonlinear_fit(start1, bounds, fixed, tied, clean)

        self.bias = 0   # Evaluate residuals without bias
        err = self._linear_fit(params)
        self.chi2 = np.sum(err**2)/(err.size - params.size)   # Chi**2/DOF
        self._format_output(params, perror)
        if plot:   # Plot final data-model comparison if required.
            self.plot()

################################################################################

    def _nonlinear_fit(self, start0, bounds0, fixed0, tied0, clean):
        """
        This function implements the procedure described in
        Section 3.4 of Cappellari M., 2017, MNRAS, 466, 798
        http://adsabs.harvard.edu/abs/2017MNRAS.466..798C

        """
        ngh = self.moments.sum()
        npars = ngh + self.mdegree * self.nspec

        if self.reddening is not None:
            npars += 1

        # Explicitly specify the step for the numerical derivatives
        # and force safety limits on the fitting parameters.
        #
        # Set [h3, h4, ...] and mult. polynomials to zero as initial guess
        # and constrain -0.3 < [h3, h4, ...] < 0.3
        #
        start = np.zeros(npars)
        fixed = np.full(npars, False)
        tied = np.full(npars, '', dtype=object)
        step = np.full(npars, 0.001)
        bounds = np.tile([-0.3, 0.3], (npars, 1))

        p = 0
        for j, st in enumerate(start0):
            if bounds0 is None:
                bn = [st[0] + np.array([-2e3, 2e3])/self.velscale,  # V bounds
                      [0.01, 1e3/self.velscale]]                # sigma bounds
            else:
                bn = np.array(bounds0[j][:2], dtype=float)/self.velscale
            for k in range(self.moments[j]):
                if self.fixall[j]:  # Negative moment --> keep entire LOSVD fixed
                    fixed[k + p] = True
                elif fixed0 is not None:  # Keep individual LOSVD parameters fixed
                    fixed[k + p] = fixed0[j][k]
                if tied0 is not None:
                    tied[k + p] = tied0[j][k]
                if k < 2:
                    start[k + p] = st[k].clip(*bn[k])
                    bounds[k + p] = bn[k]
                    step[k + p] = 0.01
                else:
                    start[k + p] = st[k]
                    if bounds0 is not None:
                        bounds[k + p] = bounds0[j][k]
            p += self.moments[j]

        if self.mdegree > 0:
            for j in range(ngh, npars):
                bounds[j] = [-0.5, 0.5]  # Force <50% corrections
        elif self.reddening is not None:
            start[ngh] = self.reddening
            bounds[ngh] = [0., 10.]  # Force positive E(B-V) < 10 mag

        if self.method == 'lm':
            step = 0.01  # only a scalar is supported
            bounds = np.tile([-np.inf, np.inf], (npars, 1))   # No bounds

        # Here the actual calculation starts.
        # If required, once the minimum is found, clean the pixels deviating
        # more than 3*sigma from the best fit and repeat the minimization
        # until the set of cleaned pixels does not change any more.
        #
        ftol = 1e-4
        good = self.goodpixels.copy()
        for j in range(5):  # Do at most five cleaning iterations
            self.clean = False  # No cleaning during chi2 optimization
            if self.method == 'capfit':
                res = capfit.capfit(self._linear_fit, start, ftol=ftol,
                                    bounds=bounds.T, abs_step=step,
                                    x_scale='jac', tied=tied, fixed=fixed)
                params = res.x
                perror = res.x_err
            else:
                res = optimize.least_squares(self._linear_fit, start, ftol=ftol,
                                             bounds=bounds.T, diff_step=step,
                                             x_scale='jac', method=self.method)
                params = res.x
                perror = capfit.cov_err(res.jac)[1]
            if not clean:
                break
            good_old = self.goodpixels.copy()
            self.goodpixels = good.copy()  # Reset goodpixels
            self.clean = True  # Do cleaning during linear fit
            self._linear_fit(params)
            if np.array_equal(good_old, self.goodpixels):
                break

        self.status = res.status
        self.njev = res.njev

        return params, perror

################################################################################

    def _linear_fit(self, pars, fjac=None):
        """
        This function implements the procedure described in
        Sec.3.3 of Cappellari M., 2017, MNRAS, 466, 798
        http://adsabs.harvard.edu/abs/2017MNRAS.466..798C

        """
        # pars = [vel_1, sigma_1, h3_1, h4_1, ... # Velocities are in pixels.
        #         ...                             # For all kinematic components
        #         vel_n, sigma_n, h3_n, h4_n, ...
        #         m1, m2, ...]                    # Multiplicative polynomials

        nspec = self.nspec
        npix = self.npix
        ngh = pars.size - self.mdegree*nspec  # Parameters of the LOSVD only
        if self.reddening is not None:
            ngh -= 1  # Fitting reddening

        nl = self.templates_rfft.shape[0]
        losvd_rfft = _losvd_rfft(pars, nspec, self.moments, nl, self.ncomp,
                                 self.vsyst, self.factor, self.sigma_diff)

        # The zeroth order multiplicative term is already included in the
        # linear fit of the templates. The polynomial below has mean of 1.
        # X needs to be within [-1, 1] for Legendre Polynomials
        x = np.linspace(-1, 1, npix)
        if self.mdegree > 0:
            if nspec == 2:  # Different multiplicative poly for left/right spectra
                mpoly1 = self.polyval(x, np.append(1.0, pars[ngh :: 2]))
                mpoly2 = self.polyval(x, np.append(1.0, pars[ngh + 1 :: 2]))
                mpoly = np.append(mpoly1, mpoly2).clip(0.1)
            else:
                mpoly = self.polyval(x, np.append(1.0, pars[ngh:])).clip(0.1)
        else:
            mpoly = 1.0

        # Multiplicative polynomials do not make sense when fitting reddening.
        # In that case one has to assume the spectrum is well calibrated.
        #
        if self.reddening is not None:
            mpoly = reddening_curve(self.lam, pars[ngh])

        skydim = len(np.shape(self.sky))  # This can be zero
        if skydim == 0:
            nsky = 0
        elif skydim == 1:
            nsky = 1  # Number of sky spectra
        else:
            nsky = np.shape(self.sky)[1]

        npoly = (self.degree + 1)*nspec  # Number of additive polynomials in fit
        ncols = npoly + nsky*nspec + self.ntemp
        nrows = npix*nspec
        if self.regul > 0:
            if self.reg_ord == 1:
                nr = self.reg_dim.size
                nreg = nr*np.prod(self.reg_dim)
            elif self.reg_ord == 2:
                nreg = np.prod(self.reg_dim)
            nrows += nreg

        if self.fraction is not None:
            nrows += 1

        # This array is used for estimating predictions
        c = np.zeros((npix*nspec, ncols))

        if self.degree >= 0:  # Fill first columns of the Design Matrix
            vand = self.polyvander(x, self.degree)
            c[: npix, : npoly//nspec] = vand
            if nspec == 2:
                c[npix :, npoly//nspec : npoly] = vand  # poly for right spectrum

        tmp = np.empty((nspec, self.npix_temp))
        for j, templates_rfft in enumerate(self.templates_rfft.T):  # columns loop
            for k in range(nspec):
                pr = templates_rfft * losvd_rfft[:, self.component[j], k]
                tt = np.fft.irfft(pr, self.npad)
                if self.factor == 1:  # Template has same resolution as galaxy
                    tmp[k, :] = tt[:self.npix_temp]
                else:                 # Template has higher resolution than galaxy
                    tmp[k, :] = rebin(tt[:self.npix_temp*self.factor], self.factor)
            c[:, npoly + j] = mpoly*tmp[:, :npix].ravel()  # reform into a vector

        if nsky > 0:
            k = npoly + self.ntemp
            c[: npix, k : k + nsky] = self.sky
            if nspec == 2:
                c[npix :, k + nsky : k + 2*nsky] = self.sky  # Sky for right spectrum

        # This array is used for the system solution
        a = np.zeros((nrows, ncols))

        if self.noise.ndim > 1 and self.noise.shape[0] == self.noise.shape[1]:
            # input NOISE is a npix*npix covariance matrix
            a[:npix*nspec, :] = self.noise.dot(c)
            b = self.noise.dot(self.galaxy)
        else:
            # input NOISE is a 1sigma error vector
            a[:npix*nspec, :] = c/self.noise[:, None] # Weight columns with errors
            b = self.galaxy/self.noise

        if self.regul > 0:
            _regularization(a, npoly, npix, nspec, self.reg_dim, self.reg_ord, self.regul)

        # Equation (30) of Cappellari (2017)
        if self.fraction is not None:
            ff = a[-1, -self.ntemp:]
            ff[self.component == 0] = self.fraction - 1
            ff[self.component == 1] = self.fraction
            ff *= 1e9

        # Select the spectral region to fit and solve the over-conditioned system
        # using SVD/BVLS. Use unweighted array for estimating bestfit predictions.
        # Iterate to exclude pixels deviating >3*sigma if /CLEAN keyword is set.

        m = 1
        while m > 0:
            if self.regul > 0 or self.fraction is not None:
                if self.regul == 0:
                    nreg = 1
                aa = a[np.append(self.goodpixels, np.arange(npix*nspec, nrows)), :]
                bb = np.append(b[self.goodpixels], np.zeros(nreg))
            else:
                aa = a[self.goodpixels, :]
                bb = b[self.goodpixels]
            self.weights = _bvls_solve(aa, bb, npoly)
            self.bestfit = c.dot(self.weights)
            if self.noise.ndim > 1 and self.noise.shape[0] == self.noise.shape[1]:
                # input NOISE is a npix*npix covariance matrix
                err = (self.noise.dot(self.galaxy - self.bestfit))[self.goodpixels]
            else:
                # input NOISE is a 1sigma error vector
                err = ((self.galaxy - self.bestfit)/self.noise)[self.goodpixels]
            if self.clean:
                w = np.abs(err) < 3  # select residuals smaller than 3*sigma
                m = err.size - w.sum()
                if m > 0:
                    self.goodpixels = self.goodpixels[w]
                    if not self.quiet:
                        print('Outliers:', m)
            else:
                break

        self.matrix = c          # Return LOSVD-convolved templates matrix

        # Penalize the solution towards (h3, h4, ...) = 0 if the inclusion of
        # these additional terms does not significantly decrease the error.
        # The lines below implement eq.(8)-(9) in Cappellari & Emsellem (2004)
        #
        if np.any(self.moments > 2) and self.bias > 0:
            D2 = p = 0
            for mom in self.moments:  # loop over kinematic components
                if mom > 2:
                    D2 += np.sum(pars[2 + p : mom + p]**2)  # eq.(8) CE04
                p += mom
            err += self.bias*robust_sigma(err, zero=True)*np.sqrt(D2)  # eq.(9) CE04

        self.nfev += 1

        return err

################################################################################

    def plot(self):
        """
        Produces a plot of the pPXF best fit.
        One can call pp.plot() after pPXF terminates.

        """
        if self.lam is None:
            plt.xlabel("Pixels")
            x = np.arange(self.galaxy.size)
        else:
            plt.xlabel(r"Wavelength [$\AA$]")
            x = self.lam

        ll, rr = np.min(x), np.max(x)
        resid = self.galaxy - self.bestfit
        mn = np.min(self.bestfit[self.goodpixels])
        mn -= np.percentile(np.abs(resid[self.goodpixels]), 99)
        mx = np.max(self.bestfit[self.goodpixels])
        resid += mn   # Offset residuals to avoid overlap
        mn1 = np.min(resid[self.goodpixels])
        plt.ylabel("Relative Flux")
        plt.xlim([ll, rr] + np.array([-0.02, 0.02])*(rr - ll))
        plt.ylim([mn1, mx] + np.array([-0.05, 0.05])*(mx - mn1))
        plt.plot(x, self.galaxy, 'k')
        plt.plot(x[self.goodpixels], resid[self.goodpixels], 'd',
                 color='LimeGreen', mec='LimeGreen', ms=4)

        w = np.flatnonzero(np.diff(self.goodpixels) > 1)
        ww = np.hstack([0, w, w + 1, -1]) if w.size > 0 else [0, -1]
        for gj in self.goodpixels[ww]:
            plt.plot(x[[gj, gj]], [mn, self.bestfit[gj]], 'LimeGreen')
        for wj in w:
            j = slice(self.goodpixels[wj], self.goodpixels[wj+1] + 1)
            plt.plot(x[j], resid[j], 'b')

        if self.gas_component is None:
            plt.plot(x, self.bestfit, 'r', linewidth=2)
            plt.plot(x[self.goodpixels], self.goodpixels*0 + mn, '.k', ms=1)
        else:
            stars_spectrum = self.bestfit - self.gas_bestfit
            plt.plot(x, self.gas_bestfit + mn, c='magenta', linewidth=2)
            plt.plot(x, self.bestfit, c='orange', linewidth=2)
            plt.plot(x, stars_spectrum, 'r', linewidth=2)

################################################################################

    def _format_output(self, params, perror):
        """
        Store the best fitting parameters in the output solution
        and print the results on the console if quiet=False

        """
        p = 0
        self.sol = []
        self.error = []
        for mom in self.moments:
            params[p : p + 2] *= self.velscale  # Bring velocity scale back to km/s
            self.sol.append(params[p : p + mom])
            perror[p : p + 2] *= self.velscale  # Bring velocity scale back to km/s
            self.error.append(perror[p : p + mom])
            p += mom
        if self.mdegree > 0:
            self.mpolyweights = params[p:]
        if self.reddening is not None:
            self.reddening = params[-1]  # Replace input with best fit
        if self.degree >= 0:
            # output weights for the additive polynomials
            self.polyweights = self.weights[: (self.degree + 1)*self.galaxy.ndim]
        # output weights for the templates (or sky) only
        self.weights = self.weights[(self.degree + 1)*self.galaxy.ndim :]

        if not self.quiet:
            nmom = np.max(self.moments)
            txt = ["Vel", "sigma"] + ["h" + str(j) for j in range(3, nmom+1)]
            print(("Best Fit:" + "{:>10}"*nmom).format(*txt))
            for j, (sol, mom) in enumerate(zip(self.sol, self.moments)):
                print((" comp. {}:" + "{:10.0f}"*2 + "{:10.3f}"*(mom-2)).format(j, *sol))
            print("chi2/DOF: {:.4g}".format(self.chi2))
            print("method = " + self.method, "; Jac calls:", self.njev,
                  '; Func calls:', self.nfev, '; Status:', self.status)
            nw = self.weights.size
            if self.reddening is not None:
                print("Reddening E(B-V): {:.3g}".format(self.reddening))
            print('Nonzero Templates: ', np.sum(self.weights > 0), ' / ', nw)
            if self.weights.size <= 20:
                print('Templates weights:')
                print(("{:10.3g}"*self.weights.size).format(*self.weights))

        if self.fraction is not None:
            fracFit = np.sum(self.weights[self.component == 0])\
                      / np.sum(self.weights[self.component < 2])
            if not self.quiet:
                print("Weights Fraction w[0]/w[0+1]: {:.3g}".format(fracFit))
            if abs(fracFit - self.fraction) > 0.01:
                print("Warning: FRACTION is inaccurate. TEMPLATES and GALAXY "
                      "should have mean ~ 1 when using the FRACTION keyword")

        if self.gas_component is not None:
            gas = self.gas_component
            spectra = self.matrix[:, self.degree + 1:]    # Remove polynomials
            integ = np.sum(spectra[:, gas], 0)
            self.gas_flux = integ*self.weights[gas]
            design_matrix = spectra[:, gas]/self.noise[:, None]
            self.gas_flux_error = integ*capfit.cov_err(design_matrix)[1]
            self.gas_bestfit = spectra[:, gas].dot(self.weights[gas])

            if not self.quiet:
                print('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
                print('gas_component   name       flux       err      V     sig')
                print('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
                for j, comp in enumerate(self.component[gas]):
                    print('Comp: %d  %12s %10.4g  %8.2g  %6.0f  %4.0f' %
                          (comp, self.gas_names[j], self.gas_flux[j], self.gas_flux_error[j],
                           self.sol[comp][0], self.sol[comp][1]))
                print('---------------------------------------------------------')

        if self.ncomp ==1:
            self.sol = self.sol[0]
            self.error = self.error[0]

################################################################################
