###############################################################################
#
# Copyright (C) 2001-2018, Michele Cappellari
# E-mail: michele.cappellari_at_physics.ox.ac.uk
#
# This software is provided as is without any warranty whatsoever.
# Permission to use, for non-commercial purposes is granted.
# Permission to modify for personal or internal use is granted,
# provided this copyright and disclaimer are included unchanged
# at the beginning of the file. All other rights are reserved.
#
###############################################################################

# This file contains the following independent programs:
#
# 1) log_rebin() to rebin a spectrum logarithmically
# 2) determine_goodpixels() to mask gas emission lines for pPXF
# 3) vac_to_air() to convert vacuum to air wavelengths
# 4) air_to_vac() to convert air to vacuum wavelengths
# 5) emission_lines() to create gas emission line templates for pPXF
# 6) gaussian_filter1d() to convolve a spectrum with a variable sigma
# 7) plot_weights_2d() to plot an image of the 2-dim weights
# 8) convolve_gauss_hermite() to accurately convolve a spectrum with a LOSVD

from __future__ import print_function

import numpy as np
from scipy import special, fftpack
import matplotlib.pyplot as plt

from . import ppxf

###############################################################################
#
# NAME:
#   LOG_REBIN
#
# MODIFICATION HISTORY:
#   V1.0.0: Using interpolation. Michele Cappellari, Leiden, 22 October 2001
#   V2.0.0: Analytic flux conservation. MC, Potsdam, 15 June 2003
#   V2.1.0: Allow a velocity scale to be specified by the user.
#       MC, Leiden, 2 August 2003
#   V2.2.0: Output the optional logarithmically spaced wavelength at the
#       geometric mean of the wavelength at the border of each pixel.
#       Thanks to Jesus Falcon-Barroso. MC, Leiden, 5 November 2003
#   V2.2.1: Verify that lamRange[0] < lamRange[1].
#       MC, Vicenza, 29 December 2004
#   V2.2.2: Modified the documentation after feedback from James Price.
#       MC, Oxford, 21 October 2010
#   V2.3.0: By default now preserve the shape of the spectrum, not the
#       total flux. This seems what most users expect from the procedure.
#       Set the keyword /FLUX to preserve flux like in previous version.
#       MC, Oxford, 30 November 2011
#   V3.0.0: Translated from IDL into Python. MC, Santiago, 23 November 2013
#   V3.1.0: Fully vectorized log_rebin. Typical speed up by two orders of magnitude.
#       MC, Oxford, 4 March 2014
#   V3.1.1: Updated documentation. MC, Oxford, 16 August 2016

def log_rebin(lamRange, spec, oversample=1, velscale=None, flux=False):
    """
    Logarithmically rebin a spectrum, while rigorously conserving the flux.
    Basically the photons in the spectrum are simply redistributed according
    to a new grid of pixels, with non-uniform size in the spectral direction.
    
    When the flux keyword is set, this program performs an exact integration 
    of the original spectrum, assumed to be a step function within the 
    linearly-spaced pixels, onto the new logarithmically-spaced pixels. 
    The output was tested to agree with the analytic solution.

    :param lamRange: two elements vector containing the central wavelength
        of the first and last pixels in the spectrum, which is assumed
        to have constant wavelength scale! E.g. from the values in the
        standard FITS keywords: LAMRANGE = CRVAL1 + [0, CDELT1*(NAXIS1 - 1)].
        It must be LAMRANGE[0] < LAMRANGE[1].
    :param spec: input spectrum.
    :param oversample: can be used, not to loose spectral resolution,
        especally for extended wavelength ranges and to avoid aliasing.
        Default: OVERSAMPLE=1 ==> Same number of output pixels as input.
    :param velscale: velocity scale in km/s per pixels. If this variable is
        not defined, then it will contain in output the velocity scale.
        If this variable is defined by the user it will be used
        to set the output number of pixels and wavelength scale.
    :param flux: (boolean) True to preserve total flux. In this case the
        log rebinning changes the pixels flux in proportion to their
        dLam so the following command will show large differences
        beween the spectral shape before and after LOG_REBIN:
           plt.plot(exp(logLam), specNew)  # Plot log-rebinned spectrum
           plt.plot(np.linspace(lamRange[0], lamRange[1], spec.size), spec)
        By defaul, when this is False, the above two lines produce
        two spectra that almost perfectly overlap each other.
    :return: [specNew, logLam, velscale] where logLam is the natural
        logarithm of the wavelength and velscale is in km/s.

    """
    lamRange = np.asarray(lamRange)
    assert len(lamRange) == 2, 'lamRange must contain two elements'
    assert lamRange[0] < lamRange[1], 'It must be lamRange[0] < lamRange[1]'
    assert spec.ndim == 1, 'input spectrum must be a vector'
    n = spec.shape[0]
    m = int(n*oversample)

    dLam = np.diff(lamRange)/(n - 1.)        # Assume constant dLam
    lim = lamRange/dLam + [-0.5, 0.5]        # All in units of dLam
    borders = np.linspace(*lim, num=n+1)     # Linearly
    logLim = np.log(lim)

    c = 299792.458                           # Speed of light in km/s
    if velscale is None:                     # Velocity scale is set by user
        velscale = np.diff(logLim)/m*c       # Only for output
    else:
        logScale = velscale/c
        m = int(np.diff(logLim)/logScale)    # Number of output pixels
        logLim[1] = logLim[0] + m*logScale

    newBorders = np.exp(np.linspace(*logLim, num=m+1)) # Logarithmically
    k = (newBorders - lim[0]).clip(0, n-1).astype(int)

    specNew = np.add.reduceat(spec, k)[:-1]  # Do analytic integral
    specNew *= np.diff(k) > 0    # fix for design flaw of reduceat()
    specNew += np.diff((newBorders - borders[k])*spec[k])

    if not flux:
        specNew /= np.diff(newBorders)

    # Output log(wavelength): log of geometric mean
    logLam = np.log(np.sqrt(newBorders[1:]*newBorders[:-1])*dLam)

    return specNew, logLam, velscale

###############################################################################
#
# NAME:
#   DETERMINE_GOODPIXELS
#
# MODIFICATION HISTORY:
#   V1.0.0: Michele Cappellari, Leiden, 9 September 2005
#   V1.0.1: Made a separate routine and included additional common emission lines.
#       MC, Oxford 12 January 2012
#   V2.0.0: Translated from IDL into Python. MC, Oxford, 10 December 2013
#   V2.0.1: Updated line list. MC, Oxford, 8 January 2014
#   V2.0.2: Use redshift instead of velocity as input for higher accuracy at large z.
#       MC, Lexington, 31 March 2015

def determine_goodpixels(logLam, lamRangeTemp, z):
    """
    Generates a list of goodpixels to mask a given set of gas emission
    lines. This is meant to be used as input for PPXF.

    :param logLam: Natural logarithm np.log(wave) of the wavelength in
        Angstrom of each pixel of the log rebinned *galaxy* spectrum.
    :param lamRangeTemp: Two elements vectors [lamMin2, lamMax2] with the minimum
        and maximum wavelength in Angstrom in the stellar *template* used in PPXF.
    :param z: Estimate of the galaxy redshift.
    :return: vector of goodPixels to be used as input for pPXF

    """
#                     -----[OII]-----    Hdelta   Hgamma   Hbeta   -----[OIII]-----   [OI]    -----[NII]-----   Halpha   -----[SII]-----
    lines = np.array([3726.03, 3728.82, 4101.76, 4340.47, 4861.33, 4958.92, 5006.84, 6300.30, 6548.03, 6583.41, 6562.80, 6716.47, 6730.85])
    dv = np.full_like(lines, 800)  # width/2 of masked gas emission region in km/s
    c = 299792.458 # speed of light in km/s

    flag = False
    for line, dvj in zip(lines, dv):
        flag |= (np.exp(logLam) > line*(1 + z)*(1 - dvj/c)) \
              & (np.exp(logLam) < line*(1 + z)*(1 + dvj/c))

    flag |= np.exp(logLam) > lamRangeTemp[1]*(1 + z)*(1 - 900/c)   # Mask edges of
    flag |= np.exp(logLam) < lamRangeTemp[0]*(1 + z)*(1 + 900/c)   # stellar library

    return np.flatnonzero(~flag)

###############################################################################

def _wave_convert(lam):
    """
    Convert between vacuum and air wavelengths using
    equation (1) of Ciddor 1996, Applied Optics 35, 1566
        http://dx.doi.org/10.1364/AO.35.001566

    :param lam - Wavelength in Angstroms
    :return: conversion factor

    """
    lam = np.asarray(lam)
    sigma2 = (1e4/lam)**2
    fact = 1 + 5.792105e-2/(238.0185 - sigma2) + 1.67917e-3/(57.362 - sigma2)

    return fact

###############################################################################

def vac_to_air(lam_vac):
    """
    Convert vacuum to air wavelengths

    :param lam_vac - Wavelength in Angstroms
    :return: lam_air - Wavelength in Angstroms

    """
    return lam_vac/_wave_convert(lam_vac)

###############################################################################

def air_to_vac(lam_air):
    """
    Convert air to vacuum wavelengths

    :param lam_air - Wavelength in Angstroms
    :return: lam_vac - Wavelength in Angstroms

    """
    return lam_air*_wave_convert(lam_air)

###############################################################################
# NAME:
#   GAUSSIAN
#
# MODIFICATION HISTORY:
#   V1.0.0: Written using analytic pixel integration.
#       Michele Cappellari, Oxford, 10 August 2016
#   V2.0.0: Define lines in frequency domain for a rigorous
#       convolution within pPXF at any sigma, including sigma=0.
#       Introduced `pixel` keyword for optional pixel convolution.
#       MC, Oxford, 26 May 2017

def gaussian(logLam_temp, line_wave, FWHM_gal, pixel=True):
    """
    Instrumental Gaussian line spread function (LSF), optionally integrated
    within the pixels. The function is normalized in such a way that
    
            line.sum() = 1
    
    When the LSF is not severey undersampled, and when pixel=False, the output
    of this function is nearly indistinguishable from a normalized Gaussian:
    
      x = (logLam_temp[:, None] - np.log(line_wave))/dx
      gauss = np.exp(-0.5*(x/xsig)**2)
      gauss /= np.sqrt(2*np.pi)*xsig

    However, to deal rigorously with the possibility of severe undersampling,
    this Gaussian is defined analytically in frequency domain and transformed
    numerically to time domain. This makes the convolution within pPXF exact
    to machine precision regardless of sigma (including sigma=0).
    
    :param logLam_temp: np.log(wavelength) in Angstrom
    :param line_wave: Vector of lines wavelength in Angstrom
    :param FWHM_gal: FWHM in Angstrom. This can be a scalar or the name of
        a function wich returns the instrumental FWHM for given wavelength.
        In this case the sigma returned by pPXF will be the intrinsic one,
        namely the one corrected for instrumental dispersion, in the same
        way as the stellar kinematics is returned.
      - To measure the *observed* dispersion, ignoring the instrumental
        dispersison, one can set FWHM_gal=0. In this case the Gaussian
        line templates reduce to Dirac delta functions. The sigma returned
        by pPXF will be the same one would measure by fitting a Gaussian
        to the observed spectrum (exept for the fact that this function
        accurately deals with pixel integration).
    :param pixel: set to True to perform integration over the pixels.
    :return: LSF computed for every logLam_temp

    """
    line_wave = np.asarray(line_wave)

    if callable(FWHM_gal):
        FWHM_gal = FWHM_gal(line_wave)

    n = logLam_temp.size
    npad = fftpack.next_fast_len(n)
    nl = npad//2 + 1  # Expected length of rfft

    dx = (logLam_temp[-1] - logLam_temp[0])/(n - 1)
    x0 = (np.log(line_wave) - logLam_temp[0])/dx
    xsig = FWHM_gal/2.355/line_wave/dx    # sigma in pixels units
    w = np.linspace(0, np.pi, nl)[:, None]

    # Gaussian with sigma=xsig and center=x0,
    # optionally convolved with an unitary pixel UnitBox[]
    # analytically defined in frequency domain
    # and numerically transformed to time domain
    rfft = np.exp(-0.5*(w*xsig)**2 - 1j*w*x0)
    if pixel:
        rfft *= np.sinc(w/(2*np.pi))
    line = np.fft.irfft(rfft, n=npad, axis=0)

    return line[:n, :]

###############################################################################
# NAME:
#   EMISSION_LINES
#
# MODIFICATION HISTORY:
#   V1.0.0: Michele Cappellari, Oxford, 7 January 2014
#   V1.1.0: Fixes [OIII] and [NII] doublets to the theoretical flux ratio.
#       Returns line names together with emission lines templates.
#       MC, Oxford, 3 August 2014
#   V1.1.1: Only returns lines included within the estimated fitted wavelength range.
#       This avoids identically zero gas templates being included in the PPXF fit
#       which can cause numerical instabilities in the solution of the system.
#       MC, Oxford, 3 September 2014
#   V1.2.0: Perform integration over the pixels of the Gaussian line spread function
#       using the new function gaussian(). Thanks to Eric Emsellem for the suggestion.
#       MC, Oxford, 10 August 2016
#   V1.2.1: Allow FWHM_gal to be a function of wavelength. MC, Oxford, 16 August 2016
#   V1.2.2: Introduced `pixel` keyword for optional pixel convolution.
#       MC, Oxford, 3 August 2017
#   V1.3.0: New `tie_balmer` keyword to assume intrinsic Balmer decrement.
#       New `limit_doublets` keyword to limit ratios of [OII] & [SII] doublets.
#       New `vacuum` keyword to return wavelengths in vacuum.
#       MC, Oxford, 31 October 2017

def emission_lines(logLam_temp, lamRange_gal, FWHM_gal, pixel=True,
                   tie_balmer=False, limit_doublets=False, vacuum=False):
    """
    Generates an array of Gaussian emission lines to be used as gas templates in PPXF.

    Generally, these templates represent the instrumental line spread function
    (LSF) at the set of wavelengths of each emission line. In this case, pPXF
    will return the intrinsic (i.e. true) dispersion of the gas lines.

    Alternatively, one can input FWHM_gal=0, in which case pPXF will return a
    dispersion which includes both the intrumental and the intrinsic disperson.

    Additional lines can be easily added by editing the code of this procedure,
    which is meant as a template to be modified by the users where needed.

    For accuracy the Gaussians are integrated over the pixels boundaries.
    This can be changed by setting `pixel`=False.

    The [OI], [OIII] and [NII] doublets are fixed at theoretical flux ratio~3.

    The [OII] and [SII] doublets can be restricted to physical range of ratios.

    The Balmet Series can be fixed to the theoretically predicted decrement.

    :param logLam_temp: is the natural log of the wavelength of the templates in
        Angstrom. logLam_temp should be the same as that of the stellar templates.
    :param lamRange_gal: is the estimated rest-frame fitted wavelength range
        Typically lamRange_gal = np.array([np.min(wave), np.max(wave)])/(1 + z),
        where wave is the observed wavelength of the fitted galaxy pixels and
        z is an initial rough estimate of the galaxy redshift.
    :param FWHM_gal: is the instrumantal FWHM of the galaxy spectrum under study
        in Angstrom. One can pass either a scalar or the name "func" of a function
        func(wave) which returns the FWHM for a given vector of input wavelengths.
    :param pixel: Set this to False to ignore pixels integration (default True).
    :param tie_balmer: Set this to True to tie the Balmer lines according to
        a theoretical decrement (case B recombination T=1e4 K, n=100 cm^-3).
    :param limit_doublets: Set this to True to limit the rato of the
        [OII] and [SII] doublets to the ranges allowed by atomic physics.
      - IMPORTANT: when using this keyword, the two output fluxes (flux_1 and
        flux_2) provided by pPXF for the two lines of the doublet, do *not*
        represent the actual fluxes of the two lines, but the fluxes of the two
        input *doublets* of which the fit is a linear combination.
        If the two doublets templates have line ratios rat_1 and rat_2, and
        pPXF prints fluxes flux_1 and flux_2, the actual ratio and flux of the
        fitted doublet will be
            flux_total = flux_1 + flux_1
            ratio_total = (rat_1*flux_1 + rat_2*flux_2)/flux_total
      - EXAMPLE: For the [SII] doublet, the adopted ratios for the templates are
            ratio_d1 = flux([SII]6716/6731) = 0.44
            ratio_d2 = flux([SII]6716/6731) = 1.43.
        If pPXF prints (and returns in pp.gas_flux)
            flux([SII]6731_d1) = flux_1
            flux([SII]6731_d2) = flux_2
        the total flux and true lines ratio of the [SII] doublet are
            flux_total = flux_1 + flux_2
            ratio_total([SII]6716/6731) = (0.44*flux_1 + 1.43*flux_2)/flux_total
      - Similarly, for [OII], the adopted ratios for the templates are
            ratio_d1 = flux([OII]3729/3726) = 0.28
            ratio_d2 = flux([OII]3729/3726) = 1.47.
        If pPXF prints (and returns in pp.gas_flux)
            flux([OII]3726_d1) = flux_1
            flux([OII]3726_d2) = flux_2
        the total flux and true lines ratio of the [OII] doublet are
            flux_total = flux_1 + flux_2
            ratio_total([OII]3729/3726) = (0.28*flux_1 + 1.47*flux_2)/flux_total
    :param vacuum: set to True to assume wavelengths are given in vacuum.
        By default the wavelengths are assumed to be measured in air.
    :return: emission_lines, line_names, line_wave

    """
    if tie_balmer:

        # Balmer decrement for Case B recombination (T=1e4 K, ne=100 cm^-3)
        # Table 4.4 of Dopita & Sutherland 2003 https://www.amazon.com/dp/3540433627
        # Balmer: Htheta  Heta     Hzeta    Heps    Hdelta   Hgamma    Hbeta   Halpha
        wave = [3797.90, 3835.39, 3889.05, 3970.07, 4101.76, 4340.47, 4861.33, 6562.80]  # air wavelengths
        if vacuum:
            wave = air_to_vac(wave)
        gauss = gaussian(logLam_temp, wave, FWHM_gal, pixel)
        emission_lines = gauss.dot([0.0530, 0.0731, 0.105, 0.159, 0.259, 0.468, 1, 2.86])
        line_names = ['Balmer']
        line_wave = np.mean(wave)

    else:

        # Use fewer lines here, as the weak ones are difficult to measure
        # Balmer:    Hdelta   Hgamma    Hbeta   Halpha
        line_wave = [4101.76, 4340.47, 4861.33, 6562.80]  # air wavelengths
        if vacuum:
            line_wave = air_to_vac(line_wave)
        line_names = ['Hdelta', 'Hgamma', 'Hbeta', 'Halpha']
        emission_lines = gaussian(logLam_temp, line_wave, FWHM_gal, pixel)


    if limit_doublets:

        # The line ratio of this doublet lam3729/lam3726 is constrained by
        # atomic physics to lie in the range 0.28--1.47 (e.g. fig.5.8 of
        # Osterbrock & Ferland 2005 https://www.amazon.co.uk/dp/1891389343/).
        # We model this doublet as a linear combination of two doublets with the
        # maximum and minimum ratios, to limit the ratio to the desired range.
        #       -----[OII]-----
        wave = [3726.03, 3728.82]    # air wavelengths
        if vacuum:
            wave = air_to_vac(wave)
        names = ['[OII]3726_d1', '[OII]3726_d2']
        gauss = gaussian(logLam_temp, wave, FWHM_gal, pixel)
        doublets = gauss.dot([[1, 1], [0.28, 1.47]])  # produces *two* doublets
        emission_lines = np.column_stack([emission_lines, doublets])
        line_names = np.append(line_names, names)
        line_wave = np.append(line_wave, wave)

        # The line ratio of this doublet lam6716/lam6731 is constrained by
        # atomic physics to lie in the range 0.44--1.43 (e.g. fig.5.8 of
        # Osterbrock & Ferland 2005 https://www.amazon.co.uk/dp/1891389343/).
        # We model this doublet as a linear combination of two doublets with the
        # maximum and minimum ratios, to limit the ratio to the desired range.
        #       -----[SII]-----
        wave = [6716.47, 6730.85]    # air wavelengths
        if vacuum:
            wave = air_to_vac(wave)
        names = ['[SII]6731_d1', '[SII]6731_d2']
        gauss = gaussian(logLam_temp, wave, FWHM_gal, pixel)
        doublets = gauss.dot([[0.44, 1.43], [1, 1]])  # produces *two* doublets
        emission_lines = np.column_stack([emission_lines, doublets])
        line_names = np.append(line_names, names)
        line_wave = np.append(line_wave, wave)

    else:

        # Here the doublets are free to have any ratio
        #       -----[OII]-----    -----[SII]-----
        wave = [3726.03, 3728.82, 6716.47, 6730.85]  # air wavelengths
        if vacuum:
            wave = air_to_vac(wave)
        names = ['[OII]3726', '[OII]3729', '[SII]6716', '[SII]6731']
        gauss = gaussian(logLam_temp, wave, FWHM_gal, pixel)
        emission_lines = np.column_stack([emission_lines, gauss])
        line_names = np.append(line_names, names)
        line_wave = np.append(line_wave, wave)


    # To keep the flux ratio of a doublet fixed, we place the two lines in a single template
    #       -----[OIII]-----
    wave = [4958.92, 5006.84]    # air wavelengths
    if vacuum:
        wave = air_to_vac(wave)
    doublet = gaussian(logLam_temp, wave, FWHM_gal, pixel).dot([0.33, 1])
    emission_lines = np.column_stack([emission_lines, doublet])
    line_names = np.append(line_names, '[OIII]5007_d')  # single template for this doublet
    line_wave = np.append(line_wave, wave[1])

    # To keep the flux ratio of a doublet fixed, we place the two lines in a single template
    #        -----[OI]-----
    wave = [6300.30, 6363.67]    # air wavelengths
    if vacuum:
        wave = air_to_vac(wave)
    doublet = gaussian(logLam_temp, wave, FWHM_gal, pixel).dot([1, 0.33])
    emission_lines = np.column_stack([emission_lines, doublet])
    line_names = np.append(line_names, '[OI]6300_d')  # single template for this doublet
    line_wave = np.append(line_wave, wave[0])

    # To keep the flux ratio of a doublet fixed, we place the two lines in a single template
    #       -----[NII]-----
    wave = [6548.03, 6583.41]    # air wavelengths
    if vacuum:
        wave = air_to_vac(wave)
    doublet = gaussian(logLam_temp, wave, FWHM_gal, pixel).dot([0.33, 1])
    emission_lines = np.column_stack([emission_lines, doublet])
    line_names = np.append(line_names, '[NII]6583_d')  # single template for this doublet
    line_wave = np.append(line_wave, wave[1])

    # Only include lines falling within the estimated fitted wavelength range.
    # This is important to avoid instabilities in the pPXF system solution
    #
    w = (line_wave > lamRange_gal[0]) & (line_wave < lamRange_gal[1])
    emission_lines = emission_lines[:, w]
    line_names = line_names[w]
    line_wave = line_wave[w]

    print('Emission lines included in gas templates:')
    print(line_names)

    return emission_lines, line_names, line_wave

###############################################################################
# NAME:
#   GAUSSIAN_FILTER1D
#
# MODIFICATION HISTORY:
#   V1.0.0: Written as a replacement for the Scipy routine with the same name,
#       to be used with variable sigma per pixel. MC, Oxford, 10 October 2015

def gaussian_filter1d(spec, sig):
    """
    Convolve a spectrum by a Gaussian with different sigma for every pixel.
    If all sigma are the same this routine produces the same output as
    scipy.ndimage.gaussian_filter1d, except for the border treatment.
    Here the first/last p pixels are filled with zeros.
    When creating a template library for SDSS data, this implementation
    is 60x faster than a naive for loop over pixels.

    :param spec: vector with the spectrum to convolve
    :param sig: vector of sigma values (in pixels) for every pixel
    :return: spec convolved with a Gaussian with dispersion sig

    """
    sig = sig.clip(0.01)  # forces zero sigmas to have 0.01 pixels
    p = int(np.ceil(np.max(3*sig)))
    m = 2*p + 1  # kernel size
    x2 = np.linspace(-p, p, m)**2

    n = spec.size
    a = np.zeros((m, n))
    for j in range(m):   # Loop over the small size of the kernel
        a[j, p:-p] = spec[j:n-m+j+1]

    gau = np.exp(-x2[:, None]/(2*sig**2))
    gau /= np.sum(gau, 0)[None, :]  # Normalize kernel

    conv_spectrum = np.sum(a*gau, 0)

    return conv_spectrum

###############################################################################
# MODIFICATION HISTORY:
#   V1.0.0: Written. Michele Cappellari, Oxford, 25 November 2016
#   V1.0.1: Set `edgecolors` keyword in pcolormesh.
#       MC, Oxford, 14 March 2017

def plot_weights_2d(xgrid, ygrid, weights, xlabel="log Age (yr)",
                    ylabel="[M/H]", title="Mass Fraction", nodots=False,
                    colorbar=True, **kwargs):
    """
    Plot an image of the 2-dim weights, as a function of xgrid and ygrid.
    This function allows for non-uniform spacing in x or y.

    """
    assert weights.ndim == 2, "`weights` must be 2-dim"
    assert xgrid.shape == ygrid.shape == weights.shape, \
        'Input arrays (xgrid, ygrid, weights) must have the same shape'

    x = xgrid[:, 0]  # Grid centers
    y = ygrid[0, :]
    xb = (x[1:] + x[:-1])/2  # internal grid borders
    yb = (y[1:] + y[:-1])/2
    xb = np.hstack([1.5*x[0] - x[1]/2, xb, 1.5*x[-1] - x[-2]/2])  # 1st/last border
    yb = np.hstack([1.5*y[0] - y[1]/2, yb, 1.5*y[-1] - y[-2]/2])

    # pcolormesh() is used below to allow for irregular spacing in the
    # sampling of the stellar population parameters (e.g. metallicity)

    ax = plt.gca()
    pc = plt.pcolormesh(xb, yb, weights.T, edgecolors='face', **kwargs)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    if not nodots:
        plt.plot(xgrid, ygrid, 'w,')
    if colorbar:
        plt.colorbar(pc)
        plt.sca(ax)  # Activate main plot before returning

    return pc

###############################################################################
# MODIFICATION HISTORY:
#   V1.0.0: Written. Michele Cappellari, Oxford, 8 February 2018

def convolve_gauss_hermite(templates, start, velscale, npix,
                           velscale_ratio=1, sigma_diff=0, vsyst=0):
    """
    Convolve a spectrum, or a set of spectra, arranged into columns of an array,
    with a LOSVD parametrized by the Gauss-Hermite series.

    This is intended to reproduce what pPXF does for the convolution and it
    uses the analytic Fourier Transform of the LOSVD introduced in

        Cappellari (2017) http://adsabs.harvard.edu/abs/2017MNRAS.466..798C

    EXAMPLE:
        ...
        pp = ppxf(templates, galaxy, noise, velscale, start,
                  degree=4, mdegree=4, vsyst=dv, velscale_ratio=velscale_ratio)

        spec = convolve_gauss_hermite(templates, pp.sol, velscale, galaxy.size,
                                      velscale_ratio=velscale_ratio, vsyst=dv)

        # The spectrum below is equal to pp.bestfit to machine precision
        spectrum = (spec @ pp.weights)*pp.mpoly + pp.apoly

    :param spectra: log rebinned spectra
    :param start: parameters of the LOSVD [vel, sig, h3, h4,...]
    :param velscale: velocity scale c*dLogLam in km/s
    :param npix: number of output pixels
    :return: vector or array with convolved spectra

    """
    npix_temp = templates.shape[0]
    templates = templates.reshape(npix_temp, -1)
    start = np.array(start)  # make copy
    start[:2] /= velscale
    vsyst /= velscale
    npad = fftpack.next_fast_len(npix_temp)

    if velscale_ratio != 1:
        assert isinstance(velscale_ratio, int), "VELSCALE_RATIO must be an integer"
        npix_temp -= npix_temp % velscale_ratio
        templates = templates[:npix_temp, :]  # Make size multiple of velscale_ratio

    templates_rfft = np.fft.rfft(templates, npad, axis=0)
    lvd_rfft = ppxf.losvd_rfft(start, 1, start.shape, templates_rfft.shape[0],
                               1, vsyst, velscale_ratio, sigma_diff)

    conv_temp = np.fft.irfft(templates_rfft*lvd_rfft[:, 0], npad, axis=0)
    conv_temp = ppxf.rebin(conv_temp[:npix*velscale_ratio, :], velscale_ratio)

    return conv_temp

################################################################################
