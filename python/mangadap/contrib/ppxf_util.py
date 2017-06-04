###############################################################################
#
# Copyright (C) 2001-2017, Michele Cappellari
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

from __future__ import print_function

import numpy as np
from scipy import special, fftpack
import matplotlib.pyplot as plt

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

def log_rebin(lamRange, spec, oversample=False, velscale=None, flux=False):
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
        standard FITS keywords: LAMRANGE = CRVAL1 + [0,CDELT1*(NAXIS1-1)].
        It must be LAMRANGE[0] < LAMRANGE[1].
    :param spec: input spectrum.
    :param oversample: Oversampling can be done, not to loose spectral resolution,
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
    :return: [specNew, logLam, velscale]

    """
    lamRange = np.asarray(lamRange)
    assert len(lamRange) == 2, 'lamRange must contain two elements'
    assert lamRange[0] < lamRange[1], 'It must be lamRange[0] < lamRange[1]'
    s = spec.shape
    assert len(s) == 1, 'input spectrum must be a vector'
    n = s[0]
    if oversample:
        m = int(n*oversample)
    else:
        m = int(n)

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
    dv = lines*0 + 800 # width/2 of masked gas emission region in km/s
    c = 299792.458 # speed of light in km/s

    flag = np.zeros_like(logLam, dtype=bool)
    for line, dvj in zip(lines, dv):
        flag |= (np.exp(logLam) > line*(1 + z)*(1 - dvj/c)) \
              & (np.exp(logLam) < line*(1 + z)*(1 + dvj/c))

    flag |= np.exp(logLam) > lamRangeTemp[1]*(1 + z)*(1 - 900/c)   # Mask edges of
    flag |= np.exp(logLam) < lamRangeTemp[0]*(1 + z)*(1 + 900/c)   # stellar library

    return np.where(flag == 0)[0]

###############################################################################

def vac_to_air(lam_vac):
    """
    Convert vacuum to air wavelengths using
    equation (1) of Ciddor 1996, Applied Optics 35, 1566
        http://dx.doi.org/10.1364/AO.35.001566

    :param lam_vac - Wavelength in Angstroms
    :return: lam_air - Wavelength in Angstroms

    """
    sigma2 = (1e4/lam_vac)**2
    fact = 1 + 5.792105e-2/(238.0185 - sigma2) + 1.67917e-3/(57.362 - sigma2)

    return lam_vac/fact

###############################################################################

def air_to_vac(lam_air):
    """
    Convert air to vacuum wavelengths using
    equation (1) of Ciddor 1996, Applied Optics 35, 1566
        http://dx.doi.org/10.1364/AO.35.001566
    :param lam_air - Wavelength in Angstroms
    :return: lam_vac - Wavelength in Angstroms

    """
    sigma2 = (1e4/lam_air)**2
    fact = 1 + 5.792105e-2/(238.0185 - sigma2) + 1.67917e-3/(57.362 - sigma2)

    return lam_air*fact

###############################################################################
# NAME:
#   EMLINE
#
# MODIFICATION HISTORY:
#   V1.0.0: Written using analytic pixel integration.
#       Michele Cappellari, Oxford, 10 August 2016
#   V2.0.0: Define lines in frequency domain for a rigorous
#       convolution within pPXF at any sigma. 
#       Introduced `pixel` keyword for optional pixel convolution.
#       MC, Oxford, 26 May 2017

def emline(logLam_temp, line_wave, FWHM_gal, pixel=True):
    """
    Instrumental Gaussian line spread function (LSF), 
    optionally integrated within the pixels. The function 
    is normalized in such a way that
    
            line.sum() = 1
    
    When the LSF is not severey undersampled, and when 
    pixel=False, the output of this function is nearly 
    indistinguishable from a normalized Gaussian:
    
      x = (logLam_temp - np.log(line_wave))/dx
      gauss = np.exp(-0.5*(x/xsig)**2)
      gauss /= np.sqrt(2*np.pi)*xsig

    However, to deal rigorously with the possibility of severe 
    undersampling, this Gaussian is defined analytically in 
    frequency domain and transformed numerically to time domain. 
    This makes the convolution exact within pPXF regardless of sigma.
    
    :param logLam_temp: np.log(wavelength) in Angstrom
    :param line_wave: Vector of lines wavelength in Angstrom
    :param FWHM_gal: FWHM in Angstrom. This can be a scalar or the
        name of a function wich returns the FWHM for given wavelength.
    :param pixel: set to True to perform integration over the pixels.
    :return: LSF computed for every logLam_temp

    """
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
#       using the new function emline(). Thanks to Eric Emsellem for the suggestion.
#       MC, Oxford, 10 August 2016
#   V1.2.1: Allow FWHM_gal to be a function of wavelength. MC, Oxford, 16 August 2016

def emission_lines(logLam_temp, lamRange_gal, FWHM_gal):
    """
    Generates an array of Gaussian emission lines to be used as gas templates in PPXF.
    These templates represent the instrumental line spread function (LSF) at the
    set of wavelengths of each emission line.

    Additional lines can be easily added by editing the code of this procedure,
    which is meant as a template to be modified by the users where needed.

    For accuracy the Gaussians are integrated over the pixels boundaries.
    This integration is only useful for quite unresolved Gaussians but one should
    keep in mind that, if the LSF is not well resolved, the input spectrum is not
    properly sampled and one is wasting useful information from the spectrograph!

    The [OI], [OIII] and [NII] doublets are fixed at theoretical flux ratio~3.

    :param logLam_temp: is the natural log of the wavelength of the templates in
        Angstrom. logLam_temp should be the same as that of the stellar templates.
    :param lamRange_gal: is the estimated rest-frame fitted wavelength range
        Typically lamRange_gal = np.array([np.min(wave), np.max(wave)])/(1 + z),
        where wave is the observed wavelength of the fitted galaxy pixels and
        z is an initial rough estimate of the galaxy redshift.
    :param FWHM_gal: is the instrumantal FWHM of the galaxy spectrum under study
        in Angstrom. One can pass either a scalar or the name "func" of a function
        func(wave) which returns the FWHM for a given vector of input wavelengths.
    :return: emission_lines, line_names, line_wave

    """

    # Balmer Series:      Hdelta   Hgamma    Hbeta   Halpha
    line_wave = np.array([4101.76, 4340.47, 4861.33, 6562.80])  # air wavelengths
    line_names = np.array(['Hdelta', 'Hgamma', 'Hbeta', 'Halpha'])
    emission_lines = emline(logLam_temp, line_wave, FWHM_gal)

    #                 -----[OII]-----    -----[SII]-----
    lines = np.array([3726.03, 3728.82, 6716.47, 6730.85])  # air wavelengths
    names = np.array(['[OII]3726', '[OII]3729', '[SII]6716', '[SII]6731'])
    gauss = emline(logLam_temp, lines, FWHM_gal)
    emission_lines = np.append(emission_lines, gauss, 1)
    line_names = np.append(line_names, names)
    line_wave = np.append(line_wave, lines)

    # To keep the flux ratio of a doublet fixed, we place the two lines in a single template
    #                 -----[OIII]-----
    lines = np.array([4958.92, 5006.84])    # air wavelengths
    doublet = 0.33*emline(logLam_temp, lines[0], FWHM_gal) + emline(logLam_temp, lines[1], FWHM_gal)
    emission_lines = np.append(emission_lines, doublet, 1)
    line_names = np.append(line_names, '[OIII]5007d') # single template for this doublet
    line_wave = np.append(line_wave, lines[1])

    # To keep the flux ratio of a doublet fixed, we place the two lines in a single template
    #                  -----[OI]-----
    lines = np.array([6300.30, 6363.67])    # air wavelengths
    doublet = emline(logLam_temp, lines[0], FWHM_gal) + 0.33*emline(logLam_temp, lines[1], FWHM_gal)
    emission_lines = np.append(emission_lines, doublet, 1)
    line_names = np.append(line_names, '[OI]6300d') # single template for this doublet
    line_wave = np.append(line_wave, lines[0])

    # To keep the flux ratio of a doublet fixed, we place the two lines in a single template
    #                 -----[NII]-----
    lines = np.array([6548.03, 6583.41])    # air wavelengths
    doublet = 0.33*emline(logLam_temp, lines[0], FWHM_gal) + emline(logLam_temp, lines[1], FWHM_gal)
    emission_lines = np.append(emission_lines, doublet, 1)
    line_names = np.append(line_names, '[NII]6583d') # single template for this doublet
    line_wave = np.append(line_wave, lines[1])

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
    assert xgrid.shape == ygrid.shape == weights.shape, \
        'Input arrays (xgrid, ygrid, weights) must have the same size'
    assert xgrid.ndim == 2, '(xgrid, ygrid, weights) must be 2-dim arrays'

    x = xgrid[:, 0]
    y = ygrid[0, :]
    xb = (x[1:] + x[:-1])/2  # grid borders
    yb = (y[1:] + y[:-1])/2
    xb = np.hstack([1.5*x[0] - x[1]/2, xb, 1.5*x[-1] - x[-2]/2])
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
