# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Implements an emission-line fitting function using pPXF.

*License*:
    Copyright (c) 2018, SDSS-IV/MaNGA Pipeline Group
        Licensed under BSD 3-clause license - see LICENSE.rst

*Source location*:
    $MANGADAP_DIR/python/mangadap/contrib/xjmc.py

*Imports and python version compliance*:
    ::

        from __future__ import division
        from __future__ import print_function
        from __future__ import absolute_import
        from __future__ import unicode_literals

        import sys
        if sys.version > '3':
            long = int

        import warnings
        import numpy as np

        from scipy import fftpack
        from scipy.ndimage import rank_filter

        import astropy.constants

        from ppxf import ppxf, capfit

*Class usage examples*:
        Add examples

*Revision history*:
    | **31 Aug 2017**: First commit by Xihan Ji (XJ)
    | **06 Feb 2018**: K. Westfall (KBW) added documentation
    | **09 Feb 2018**: (KBW) Return the bin matching vector
    | **20 Mar 2018**: (KBW) Correct error in carrying around pixels
        rejected during fit
    | **22 May 2018**: (KBW) Change import to ppxf package.
    | **29 May 2018**: (KBW) Remove original function from Xihan and
        rename *_edit function.

"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import sys
if sys.version > '3':
    long = int

import warnings
import numpy as np

from scipy import fftpack
from scipy.ndimage import rank_filter

import astropy.constants

from ppxf import ppxf, capfit

# For debugging
from matplotlib import pyplot as plt
# from mangadap.proc.ppxffit import PPXFModel

def calculate_noise(residuals, width=101):
    """
    Robust determination of the error spectrum as 1/2 of the interval
    enclosing 68% of the values in a given running window.  Created by
    Michele Cappellari (23 September 2017)

    Args:

        residuals (numpy.ndarray): Vector of residuals between the model
            and data
        width (int): (**Optional**) Size of the window for the
            statistics.  Should be an odd number and larger and 10.
            Default is 101.
    
    Returns:
        numpy.ndarray: Vector with the same size as the input residuals
        with half of the 68% confidence interval of the residuals within
        a box of size `width` centered at each position.

    Raises:
        ValueError: Raised if the width is not an odd number or if it is
            less than 11 elements.
    """
    if width%2 != 1:
        raise ValueError('Must provided odd number width.')
    if width < 10:
        raise ValueError('Width must be 11 or higher.')
    
    p = (1 - 0.6827)/2  # 1sigma
    lo, up = round(width*p), round(width*(1 - p))

    upper = rank_filter(residuals, rank=up, size=width)
    lower = rank_filter(residuals, rank=lo, size=width)
    noise = (upper - lower)/2

    return noise


def _ppxf_component_setup(component, gas_template, start, single_gas_component=False,
                          gas_start=None):
    """
    Setup the component, moment, and starting point arrays for a
    ppxf-fitting interation.

    This primarily just sets all the gas components to a single
    component if requested.  It also restructures the moment and
    starting point arrays as necessary to reflect this change.

    ..warning::

        This function requires that all stellar components have
        component numbers that are **less than** any gas components.
        However, this **is not checked** by the function.

    Args:
        component (numpy.ndarray): The full list of components for all
            templates.  The total number of components is NCOMP.
        gas_template (numpy.ndarray): A boolean array identifying the
            gas templates.
        start (list, numpy.ndarray): The starting kinematics to use for
            each object spectrum.  Shape is (NOBJ,); each element has
            shape (NCOMP,); each component element has shape (NMOM,),
            such that the number of moments can be different for each
            component.  NCOMP and NMOM should be the same for all object
            spectra.
        single_gas_component (bool): (**Optional**) Flag to force all
            gas components to have the same kinematics.
        gas_start (list, numpy.ndarray): (**Optional**) The starting
            kinematics to use for the gas components.  Shape must be
            (NOBJ,2).

    Returns:
        numpy.ndarray: Three arrays: (1) The list of downselected and
        renumbered, if necessary, kinematic components, (2) the
        downselected and reordered, if necessary, number of moments to
        fit, and (3) the downselected and reordered starting guesses for
        each component.

    Raises:
        ValueError: Raised if the provided gas starting kinematics is
            not correct.
    """
    # Number of object spectra to fit
    nobj = len(start)

    if gas_start is not None and np.asarray(gas_start).shape != (nobj,2):
        raise ValueError('Provided gas starting kinematics has incorrect shape.')

    # Get the gas component numbers
    gas_components = np.unique(component[gas_template])

    # Component assignments for each template
    _component = component.copy()
    n_gas_components = 1 if single_gas_component else len(gas_components)
    if n_gas_components == 1:
        _component[gas_template] = np.amax(component[~gas_template])+1

    # Shape of start is (nobj,); each element of start has shape (ncomp,);
    # each component has shape (nmom,), which can be different for each
    # component
    stellar_components = np.arange(np.amin(gas_components))
    _moments = np.array([ len(s) for s in start[0] ])
    _moments[stellar_components] *= -1
    if n_gas_components == 1:
        # Force a single gas component
        _moments = np.append( _moments[_moments < 0], [ 2 ])

    # Reset the starting values for the kinematics
    _start = np.empty(len(start), dtype=object)
    for i,s in enumerate(start):
        _start[i] = s[stellar_components].tolist() if len(stellar_components) > 0 else []
        if gas_start is None:
            _gas_start = [np.mean(np.asarray(s[gas_components]), axis=0).tolist()] \
                                if n_gas_components == 1 else [s[gas_components]]
        else:
            _gas_start = [gas_start[i].tolist()]*n_gas_components
        _start[i] += _gas_start
        
    return _component, _moments, _start

# NOT USED
#def _ppxf_fit_matrix(start, npix_obj, templates_rfft, npix_tpl, npad, component, vsyst, moments,
#                     velscale, velscale_ratio):
#    """
#    Mimic the first bit of ppxf._linear_fit()
#
#    Limited as follows:
#        - no additive or multiplicative polynomial
#        - no reddening
#        - only single galaxy spectrum (ppxf.nspec = 1)
#        - no sigma_diff
#    """
#    # Template and component shape
#    ntpl, nl = templates_rfft.shape
#    if len(component) != ntpl:
#        raise ValueError('Component length must match the number of templates.')
#    ncomp = np.max(component)+1
#
#    # Convert from km/s to pixel coordinates
#    _start = np.array(start) #, dtype=object)
#    for j in range(ncomp):
#        _start[j][:2] /= velscale
#    _vsyst = vsyst/velscale
#
#    # Parameter vector:
#    # pars = [vel_1, sigma_1, h3_1, h4_1, ... # Velocities are in pixels.
#    #         ...                             # For all kinematic components
#    #         vel_n, sigma_n, h3_n, h4_n, ...
#    pars = _start.ravel()
#
#    # FFT of the LOSVD
#    losvd_rfft = ppxf._losvd_rfft(pars, 1, np.abs(moments), nl, ncomp, _vsyst, velscale_ratio, 0.)
#
#    # Construct and return the convolved templates
#    c = np.zeros((ntpl, npix_obj), dtype=float)
#    for i, tpl_rfft in enumerate(templates_rfft):
#        pr = tpl_rfft * losvd_rfft[:, component[i], 0]
#        tt = np.fft.irfft(pr, npad)
#        c[i,:] = tt[:npix_obj] if velscale_ratio == 1 else \
#                    ppxf.rebin(tt[:npix_tpl*velscale_ratio], velscale_ratio)[:npix_obj]
#    return c


def _validate_templates_components(templates, gas_template, component, moments, mask, start, tied,
                                   tpl_to_use, velscale, velscale_ratio=None, vsyst=0):
    """
    Check that templates are valid in the sense that they have non-zero
    components in valid pixel regions.  If all templates are valid, the
    returned data is identical to the input.  If not, the returned data
    are restructured as necessary for running pPXF.
    
    The start vector is used to set the guess offset (including vsyst)
    between the provided mask and the template.  The tpl_to_use vector
    is used to set a tempate as invalid.

    ..warning::
        Unclear what happens if an entire component is lost!

    Args:
        templates (numpy.ndarray): Templates library to use for fitting.
            Shape is (NTPLPIX,NTPL).
        gas_template (numpy.ndarray): Boolean vector that selects the
            gas templates.  Shape is (NTPL,).
        component (numpy.ndarray): Integer vector identifying the
            kinematic component for each template.  Shape is (NTPL,).
        moments (numpy.ndarray): Integer vector with the number of
            kinematic moments for each component.  Shape is (NCOMP,).
        mask (numpy.ndarray): Boolean vector that selects the pixels in
            the object spectrum to fit (i.e., mask=True for pixels to
            fit and mask=False for pixels to ignore).  As in pPXF, the
            length is expected to be less than or equal to NTPLPIX in
            the template spectra.
        start (list): The starting kinematics for each kinematic
            component.  Length is NCOMP.
        tied (list): The pPXF tying list establishing which parameters
            should be tied during the fit.  Length is NCOMP, and each
            component has length NMOM.
        tpl_to_use (numpy.ndarray):  In addition to removing
            unconstrained templates, this boolean vector selects that
            should even be considered in the fit.  Shape is (NTPL,).
        velscale (float): The pixel scale of the object spectrum to fit
            in km/s.
        velscale_ratio (int): (**Optional**) The ratio between the
            object and template pixel scale.
        vsyst (float): (**Optional**) The pseudo velocity shift between
            the template and object spectra just due to the difference
            in the starting wavelength.

    Returns:
        numpy.ndarray: Eight arrays are returned: (1) Boolean vector
        with which templates were valid (length is NTPL), (2) the valid
        set of templates to pass to pPXF [shape is (NTPLPIX,_NTPL)], (3)
        the boolean vector selecting gas templates (length is _NTPL),
        (4) an integer vector mapping the input component number to the
        output component number (i.e., component_map[0] is the original
        component number for the downselected 0th component; length is
        _NCOMP), (5) the new component number for the downselected
        templates (length is _NTPL), (6) the number of moments for the
        new components (length is _NCOMP), (7) the starting kinematics
        for each new component (length is _NCOMP), and (8) the parameter
        tying object reordered as necessary for the new component list
        (**this needs to be checked**).

    Raises:
        ValueError: Raised if no templates are valid.

    """
#    # Get the FFT of the templates
#    npix_tpl = templates.shape[1]
#    npad = fftpack.next_fast_len(npix_tpl)
#    templates_rfft = np.fft.rfft(templates, npad, axis=1)
#
#    # Get the design matrix without any polynomial adjustment
#    c = _ppxf_fit_matrix(start, len(mask), templates_rfft, npix_tpl//velscale_ratio, npad,
#                         component, vsyst, moments, velscale, velscale_ratio)
#
#    # Determine which templates have constraining data
#    valid = np.max(np.absolute(c * mask.astype(float)[None,:]), axis=1) > 1e3*np.finfo(float).eps
#    valid &= tpl_to_use
#    ncomp = np.max(component)+1
#
#    for i in range(c.shape[0]):
#        if not gas_template[i]:
#            continue
#        plt.plot(c[i,:], zorder=1, linestyle='-' if valid[i] else ':')
#    plt.fill_between(np.arange(c.shape[1]),0, mask.astype(float)/4,
#                     color='k', alpha=0.3, zorder=0, lw=0)
#    plt.show()

    # Pulled from the pPXF test for consistency
    npix_tpl = templates.shape[1]
    vmed = np.median([a[0] for a in start])
    dx = int((vsyst + vmed)/velscale)  # Approximate velocity shift
    c = templates if velscale_ratio is None else \
            np.mean(templates.reshape(-1, npix_tpl//velscale_ratio, velscale_ratio), axis=2)
    c = c[:,:len(mask)]
    m1 = np.max(np.abs(c), axis=1)
    c = np.roll(c, dx, axis=1)
    m2 = np.max(np.abs(c * mask.astype(float)[None,:]), axis=1)
    valid = m2 > m1/1e3
    valid &= tpl_to_use
    ncomp = np.max(component)+1

#    for i in range(c.shape[0]):
#        if not gas_template[i]:
#            continue
#        plt.plot(c[i,:], zorder=1, linestyle='-' if valid[i] else ':')
#    plt.fill_between(np.arange(c.shape[1]),0, mask.astype(float)/4,
#                     color='k', alpha=0.3, zorder=0, lw=0)
#    plt.show()

    # All templates are valid, just return the input
    if np.all(valid):
        return valid, templates, gas_template, np.arange(ncomp), component, moments, start, tied

    # None (!) of the templates are valid
    if np.sum(valid) == 0:
        raise ValueError('No valid templates in fit!')

    # Templates to fit
    _templates = templates[valid,:]
    # Set which are gas templates
    _gas_template = gas_template[valid]
    # Sequential template numbers
    _component = component[valid]
    # Map between old and new component numbers
    component_map = np.unique(_component)

    if not np.array_equal(component_map, np.arange(ncomp)):
        # Subset of components are not in sequence so they need to be remapped
        remapped_component = _component.copy()
        # Starting kinematics for each component
        _start = []
        _tied = None if tied is None else []
        for i, cm in enumerate(component_map):
            indx = _component == cm
            remapped_component[indx] = i
            _start += [start[cm]]
            if tied is not None:
                _tied += [tied[cm]]
        _component = remapped_component
        _moments = moments[component_map]
    else:
        _start = start
        _tied = tied
        _moments = moments

    return valid, _templates, _gas_template, component_map, _component, _moments, _start, _tied


def _reorder_solution(ppsol, pperr, component_map, moments, start=None, fill_value=-999.):
    """
    Provided the best-fitting parameters from a ppxf instance,
    restructure the parameters to account for components that were
    removed during the template validation.

    Args:
        ppsol (list): The pPXF solution parameters
        pperr (list): The formal errors on the pPXF solution parameters
        component_map (numpy.ndarray): An integer vector mapping the
            input component number to the output component number; i.e.,
            component_map[0] is the original component number for the
            downselected 0th component; length is _NCOMP.
        moments (list): The number of moments for the *original* set of
            components; length is NCOMP.
        start (list): (**Optional**) The starting kinematics for the
            *original* set of components.  If provided, the starting
            values are included in the output parameter object for
            components that were not included in the pPPXF fit;
            otherwise, the unfit components are given placeholder
            parameter values.
        fill_value (float): (**Optional**) The placeholder value to give
            parameters and errors for components not included in the
            pPXF fit.

    Returns:
        list: Two lists with the rearranged best-fitting parameters and
        errors.
    """
    # If all the components were fit, just return the input
    ncomp = len(moments)
    if np.array_equal(component_map, np.arange(ncomp)):
        return ppsol, pperr
    _ncomp = len(component_map)

    # Rearrange the output to match the input components
    new_start = []
    new_error = []
    for i in range(ncomp):
        indx = [component_map == i]
        if np.sum(indx) == 0:
            new_start += ([[fill_value]*abs(moments[i])] if start is None else [start[i]])
            new_error += [[fill_value]*abs(moments[i])]
            continue
        indx = np.arange(_ncomp)[indx][0]
        new_start += [ppsol[indx]]
        new_error += [pperr[indx]]
    return new_start, new_error


def _fit_iteration(templates, wave, flux, noise, velscale, start, moments, component, gas_template,
                   tpl_to_use=None, reject_boxcar=101, velscale_ratio=None, degree=-1, mdegree=0,
                   reddening=None, tied=None, mask=None, vsyst=0, plot=False, quiet=True, sigma_rej=3.):
    """
    Run a single fit+rejection iteration of the pPXF fit for all input
    spectra with the provided set of constraints/options.

    Args: 
        templates (numpy.ndarray): Full template library.  Shape is
            (NTPL,NTPLPIX).
        wave (numpy.ndarray): Wavelength vector.  Shape is (NPIX,).
        flux (numpy.ndarray): Object spectra to fit.  Shape is
            (NSPEC,NPIX).
        noise (numpy.ndarray): Error in the object spectra to fit.
            Shape is (NSPEC,NPIX).
        velscale (float): The pixel scale of the object spectra in km/s.
        start (list): The starting kinematics for each kinematic
            component.  Length is NCOMP.
        moments (numpy.ndarray): Integer vector with the number of
            kinematic moments for each component.  Shape is (NCOMP,).
        component (numpy.ndarray): Integer vector identifying the
            kinematic component for each template.  Shape is (NTPL,).
        gas_template (numpy.ndarray): Boolean vector that selects the
            gas templates.  Shape is (NTPL,).
        tpl_to_use (numpy.ndarray):  (**Optional**) Boolean vector
            selecting templates to consider during the fit.  Shape is
            (NTPL,).  If None, all templates are used in the fit.
        reject_boxcar (int): (**Optional**) Size of the window for the
            rejection statistics.  Should be an odd number and larger
            and 10.  Default is 101.  If None, no rejection iteration is
            performed.
        velscale_ratio (int): (**Optional**) The ratio between the
            object and template pixel scale.
        degree (int): (**Optional**) Order of the additive polynomial to
            include in the fit.  Not included by default.
        mdegree (int): (**Optional**) Order of the multiplicative
            polynomial to fit.  Not included by default.
        reddening (float): (**Optional**) Initial E(B-V) guess for
            reddening (uses ppxf-default Calzetti 2000 model).  No
            attentuation fit by default.
        tied (list): (**Optional**) List of parameters to tie during the
            fit.
        mask (numpy.ndarray): (**Optional**) Boolean vector that selects
            the pixels in the object spectra to fit (i.e., mask=True for
            pixels to fit and mask=False for pixels to ignore).  Shape
            is (NSPEC,NPIX).  All pixels are fit by default.
        vsyst (float): (**Optional**) The pseudo velocity shift between
            the template and object spectra just due to the difference
            in the starting wavelength.
        plot (bool): (**Optional**) Show the pPXF fit plot at each
            iteration.  Default is to skip the plot.
        quiet (bool): (**Optional**) Suppress output to the terminal
            (default).
        sigma_rej (float): (**Optional**) Sigma values used for the rejection.
            Defalut is 3.
            

    Returns:
        numpy.ndarray: Eleven arrays are returned: (1) The best-fitting
        model for each spectrum [shape is (NSPEC,NPIX)]; (2) the
        best-fitting emission-line-only model for each spectrum [shape
        is (NSPEC,NPIX)]; (3) a boolean array that is True for all
        spectral pixels included in the fit [shape is (NSPEC,NPIX)]; (4)
        the best-fitting weight for each template in each spectrum
        [shape is (NSPEC,NTPL)]; (5) the error in the best-fitting
        template weights [shape is (NSPEC,NTPL)]; (6) the coefficients
        of the additive polynomial for each spectrum [shape is
        (NSPEC,DEGREE+1)]; (7) the coefficients of the multiplicative
        polynomial for each spectrum [shape is (NSPEC,MDEGREE)]; (8) the
        best-fitting reddening values for each spectrum [shape is
        (NSPEC,)]; (9) the input kinematics for each fit [shape is
        (NSPEC,sum(MOMENTS)]; (10) the best-fit kinematics for each fit
        [shape is (NSPEC,sum(MOMENTS)]; (11) the formal error in the
        best-fit kinematics for each fit [shape is (NSPEC,sum(MOMENTS)].

    """
    # Some useful shape numbers
    nspec = flux.shape[0]
    ntpl = templates.shape[0]
    nkin = np.sum(np.absolute(moments))

    # Establish which templates should be used
    _tpl_to_use = np.ones((nspec,ntpl), dtype=bool) if tpl_to_use is None else tpl_to_use.copy()

    # Initialize the output
    model = np.zeros(flux.shape, dtype=float)
    eml_model = np.zeros(flux.shape, dtype=float)
    model_mask = np.ones(flux.shape, dtype=bool) if mask is None else mask.copy()
    tpl_wgts = np.zeros((nspec,ntpl), dtype=float)
    tpl_wgts_err = np.zeros((nspec,ntpl), dtype=float)
    addcoef = None if degree < 0 else np.zeros((nspec,degree+1), dtype=float)
    multcoef = None if mdegree < 1 else np.zeros((nspec,mdegree), dtype=float)
    ebv = None if reddening is None else np.zeros(nspec, dtype=float)
    kininp = np.zeros((nspec,nkin), dtype=float)
    kin = np.zeros((nspec,nkin), dtype=float)
    kin_err = np.zeros((nspec,nkin), dtype=float)

    # For debugging
    linear=False
#    linear=True
#    reject_boxcar=None
#    mdegree=0

    # Fit each spectrum individually
    for i in range(nspec):

        # Report progress
        print('Fitting spectrum: {0}/{1}'.format(i+1,nspec), end='\r')

        # Run the first fit
        if plot:
            plt.clf()
        valid_templates, _templates, _gas_template, component_map, _component, _moments, _start, \
                _tied = _validate_templates_components(templates, gas_template, component, moments,
                                                       model_mask[i,:], start[i], tied,
                                                       _tpl_to_use[i,:], velscale,
                                                       velscale_ratio=velscale_ratio, vsyst=vsyst)
        pp = ppxf.ppxf(_templates.T, flux[i,:], noise[i,:], velscale, _start,
                       velscale_ratio=velscale_ratio, plot=plot, moments=_moments, degree=degree,
                       mdegree=mdegree, lam=wave, reddening=reddening, tied=_tied,
                       mask=model_mask[i,:], vsyst=vsyst, component=_component,
                       gas_component=_gas_template, quiet=quiet, linear=linear)
        if plot:
            plt.show()

        # Reject 3-sigma outliers and refit, if requested by a provided
        # boxcar width
        if reject_boxcar is not None:
            # - Calculate residuals
            resid = flux[i,:] - pp.bestfit
            # - Select pixels included in the fit and not fit by
            # emission lines
            reject_pixels = list(set(pp.goodpixels)
                                    & set(np.arange(len(resid))[pp.gas_bestfit < 1e-6]))
            # - Calculate the 1-sigma confidence interval
            NOISE = calculate_noise(resid[reject_pixels], width=reject_boxcar)
            # - Reject pixels with > 3-sigma residuals
            model_mask[i,reject_pixels] &= (abs(resid[reject_pixels]) < (sigma_rej*NOISE))
        
            # Reorder the output; sets any omitted components to have
            # the starting values from the original input
            sol, err = _reorder_solution(pp.sol, pp.error, component_map, moments, start=start[i])

            # Refit using best-fit kinematics from previous fit
            if plot:
                plt.clf()
            valid_templates, _templates, _gas_template, component_map, _component, _moments, \
                    _start, _tied = _validate_templates_components(templates, gas_template,
                                                                   component, moments,
                                                                   model_mask[i,:], sol, tied,
                                                                   _tpl_to_use[i,:], velscale,
                                                                   velscale_ratio=velscale_ratio,
                                                                   vsyst=vsyst)
            pp = ppxf.ppxf(_templates.T, flux[i,:], noise[i,:], velscale, _start,
                           velscale_ratio=velscale_ratio, plot=plot, moments=_moments,
                           degree=degree, mdegree=mdegree, lam=wave, reddening=reddening,
                           tied=_tied, mask=model_mask[i,:], vsyst=vsyst, component=_component,
                           gas_component=_gas_template, quiet=quiet, linear=linear)
            if plot:
                plt.show()

#        # npoly assumed to be zero
#        print('GAS TEMPLATE SUM:')
#        print(pp.matrix.shape)
#        print(pp.lam.shape)
#        s = np.sum(pp.matrix[:,pp.gas_component], axis=0)
#        ws = np.sum(pp.matrix[:,pp.gas_component]*pp.weights[pp.gas_component], axis=0)
#        wi = np.sum((pp.matrix[:,pp.gas_component]*pp.weights[None,pp.gas_component])[:-1,:]
#                            * np.diff(pp.lam)[:,None], axis=0)
#
#        mw = np.ma.divide(np.sum((pp.lam[:,None] * pp.matrix[:,pp.gas_component] 
#                                    * pp.weights[None,pp.gas_component])[:-1,:]
#                                        * np.diff(pp.lam)[:,None], axis=0),wi)
#        fmw = np.mean(pp.lam)
#        gc = np.unique(pp.component[pp.gas_component]).astype(int)
#        print(gc)
#        print([np.exp(pp.sol[i][0]/astropy.constants.c.to('km/s').value) for i in gc])
#        print(s)
#        print(ws)
#        print(wi)
#        print(np.ma.divide(wi,ws))
#        print(mw)
#        print(mw/fmw)
#        print(np.ma.divide(wi,ws)*(fmw/mw))
#
#        print('GAS FLUXES:')
#        print(pp.gas_flux)

        # Reorder the output; sets any omitted components to a default
        # value of -999.
        sol, err = _reorder_solution(pp.sol, pp.error, component_map, moments)

        # Save the results
        model[i,:] = pp.bestfit
        eml_model[i,:] = np.dot(pp.matrix[:,degree+1:][:,_gas_template], pp.weights[_gas_template])

        tpl_wgts[i,valid_templates] = pp.weights.copy()
        design_matrix = pp.matrix/pp.noise[:, None]
        _, tpl_wgts_err[i,valid_templates] = capfit.cov_err(design_matrix)

        if degree > -1:
            addcoef[i,:] = pp.polyweights.copy()
        if mdegree > 0:
            multcoef[i,:] = pp.mpolyweights.copy()
        if reddening is not None:
            ebv[i] = pp.reddening

        kininp[i,:] = np.concatenate(tuple(start[i]))
        kin[i,:] = np.concatenate(tuple(sol))
        kin_err[i,:] = np.concatenate(tuple(err))
        
#        plt.plot(flux[i,:], color='b', lw=2)
#        plt.plot(np.ma.MaskedArray(flux[i,:], mask=np.invert(model_mask[i,:])), color='k', lw=2)
#        plt.plot(model[i,:], color='C3', lw=1)
#        plt.plot(model[i,:]-eml_model[i,:], color='C2', lw=1)
#        plt.show()

    # Done
    print('Fitting spectrum: {0}/{0}'.format(nspec))

    return model, eml_model, model_mask, tpl_wgts, tpl_wgts_err, addcoef, multcoef, ebv, \
                    kininp, kin, kin_err


def _combine_stellar_templates(templates, gas_template, wgts, component):
    """
    Construct the optimal stellar template for each fitted spectrum.

    Args:
        templates (numpy.ndarray): Templates library to use for fitting.
            Shape is (NTPLPIX,NTPL).
        gas_template (numpy.ndarray): Boolean vector that selects the
            gas templates.  Shape is (NTPL,).
        wgts (numpy.ndarray):  The best-fitting template weights for
            each template in each spectrum.  Shape is (NSPEC,NTPL).
        component (numpy.ndarray): Integer vector identifying the
            kinematic component for each template.  Shape is (NTPL,).

    Returns:
        numpy.ndarray: Five arrays are returned: (1) the weights of only
        the stellar templates used in constructing each optimal template
        [shape is (NSPEC,NTPL)]; (2) the optimal stellar templates for
        each spectrum with the gas templates appended [shape is
        (NSPEC+NGASTPL,NPIXTPL)]; (3) a boolean array that selects the
        gas templates [shape is (NSPEC+NGASTPL,); (4) boolean array
        selecting which templates to use with each spectrum [shape is
        (NSPEC,NSPEC+NGASTPL)]; (5) integer array setting the component
        associated with each template [shape is (NSPEC+NGASTPL)].

    Raises:
        ValueError: Raised if there is more than one stellar component.

    """
    # Get the best-fit templates for each bin
    stellar_wgts = wgts.copy()
    stellar_wgts[:,gas_template] = 0.0
    optimal_template = np.dot(stellar_wgts, templates)

    # Update the list of templates
    _templates = np.append(np.atleast_2d(optimal_template), templates[gas_template,:], axis=0)
    _gas_template = np.ones(_templates.shape[0], dtype=bool)
    _gas_template[:optimal_template.shape[0]] = False
   
    # Set which templates to use with each spectrum
    nspec = wgts.shape[0]
    _tpl_to_use = np.zeros((nspec,_templates.shape[0]), dtype=bool)
    _tpl_to_use[:,nspec:] = True
    _tpl_to_use[np.arange(nspec),np.arange(nspec)] = True

    # Update the components
    _component = np.append(np.zeros(nspec).astype(int), component[gas_template])

    # Check that the components make sense
    if np.any(np.unique(_component) != np.arange(np.amax(_component)+1)):
        raise ValueError('Problem constructing new component vector.  Likely more than one '
                         'stellar component.')

    # Return new template list
    return stellar_wgts, _templates, _gas_template, _tpl_to_use, _component


def emline_fitter_with_ppxf(templates, wave, flux, noise, mask, velscale, velscale_ratio,
                            inp_component, gas_template, inp_moments, inp_start, tied=None,
                            degree=-1, mdegree=0, reddening=None, reject_boxcar=101, vsyst=0,
                            tpl_to_use=None, flux_binned=None, noise_binned=None, mask_binned=None,
                            x_binned=None, y_binned=None, x=None, y=None, plot=False, quiet=False,
                            debug=False, sigma_rej=3.):

    """
    Main calling function for fitting stellar-continuum and nebular
    emission lines in many spectra using pPXF.

    This is a generalization of :func:`emline_fitter_with_ppxf` provided
    by Xihan Ji and Michele Cappellari.
    
    The function *does not* fit for the stellar kinematics; these are
    fixed during the fit (as provided in `inp_start`) and should have
    resulted from a previous fit to the stellar-continuum only.  The
    number of stellar kinematic moments must have been the same for all
    spectra, and there can only be one stellar component.

    The templates are expected to be an integer number of velscale_ratio
    in length; see
    :func:`mangadap.proc.ppxffit.PPXFFit.check_templates`.

    The fitting procedure can be performed for a single set of spectra,
    or a set of spectra that are remapped from an input set of binned
    spectra.  The latter operation is selected by providing the binned
    flux, error, and mask arrays and the binned and unbinned on-sky
    coordinates (see argument description below).

    When the binned spectra are provided, the fitting procedure is as
    follows:
        - The binned spectra are fit with the stellar components fixed
          to the provided kinematics in the starting value array and
          with all the gas templates part of a single kinematic
          component.  This first fit iteration includes any rejection
          iteration according to the provided boxcar width.
        - The nearest binned spectrum is then found for each provided
          spectrum based on the provided coordinates.  The fit to the
          binned spectrum is then used to set (1) the single optimal
          stellar template and (2) the initial guess for the gas
          kinematics to use for the subsequent fits to each individual
          spectrum
        - Each spectrum is fit for the first time with the stellar
          kinematics fixed to the result for the associated binned
          spectrum and, again, with all the gas templates free but part
          of the same kinematic component.  This fit iteration includes
          any rejection iteration according to the provided boxcar
          width.
        - Finally the spectra are fit without any rejection iteration
          and allowing the gas templates to be associated with multiple
          kinematic components as requested by the user.

    When no binned spectra are provided, the procedure is virtually the
    same except the initial fit to the binned spectra is skipped.  An
    initial fit to the spectra is performed to construct the optimal
    template, instead of basing the optimal template on the initial fit
    to the binned spectrum.

    Args:
        templates (numpy.ndarray): Templates library to use for fitting.
            Shape is (NTPLPIX,NTPL).
        wave (numpy.ndarray): Wavelength vector.  Shape is (NPIX,).
        flux (numpy.ndarray): Object spectra to fit.  Shape is
            (NSPEC,NPIX).
        noise (numpy.ndarray): Error in the object spectra to fit.
            Shape is (NSPEC,NPIX).
        mask (numpy.ndarray): Boolean vector that selects the pixels in
            the object spectra to fit (i.e., mask=True for pixels to fit
            and mask=False for pixels to ignore).  Shape is
            (NSPEC,NPIX).  All pixels are fit by default.
        velscale (float): The pixel scale of the object spectra in km/s.
        velscale_ratio (int): The ratio between the object and template
            pixel scale; must be an integer.
        inp_component (numpy.ndarray): Integer vector identifying the
            kinematic component for each template.  Shape is (NTPL,).
        gas_template (numpy.ndarray): Boolean vector that selects the
            gas templates.  Shape is (NTPL,).
        inp_moments (numpy.ndarray): Integer vector with the number of
            kinematic moments for each component.  Shape is (NCOMP,).
        inp_start (list, numpy.ndarray): The starting kinematics to use
            for each spectrum.  Shape is (NSPEC,); each element has
            shape (NCOMP,); each component element has shape (NMOM,),
            such that the number of moments can be different for each
            component.  NCOMP and NMOM should be the same for all object
            spectra.
        tied (list): (**Optional**) List of parameters to tie during the
            fit.  Shape is (NCOMP,).
        degree (int): (**Optional**) Order of the additive polynomial to
            include in the fit.  Not included by default.
        mdegree (int): (**Optional**) Order of the multiplicative
            polynomial to fit.  Not included by default.
        reddening (float): (**Optional**) Initial E(B-V) guess for
            reddening (uses ppxf-default Calzetti 2000 model).  No
            attentuation fit by default.
        reject_boxcar (int): (**Optional**) Size of the window for the
            rejection statistics.  Should be an odd number and larger
            and 10.  Default is 101.  If None, no rejection iterations
            are performed.
        vsyst (float): (**Optional**) The pseudo velocity shift between
            the template and object spectra just due to the difference
            in the starting wavelength. Default is 0 km/s.
        tpl_to_use (numpy.ndarray):  (**Optional**) Boolean vector
            selecting templates to consider during the fit.  If None,
            all templates are used in the fit.  If provided, the shape
            must be (NBIN,NTPL) when providing the binned data and
            (NSPEC,NTPL) when no binned data are provided.
        flux_binned (numpy.ndarray): (**Optional**) Binned spectra with
            previous fits to the stellar kinematics.  See purpose above.
        noise_binned (numpy.ndarray): (**Optional**) Error in the binned
            spectra.
        mask_binned (numpy.ndarray): (**Optional**) Mask for the binned
            spectra.
        x_binned (numpy.ndarray): (**Optional**) On-sky bin x
            coordinate; shape is (NBIN,).
        y_binned (numpy.ndarray): (**Optional**) On-sky bin y
            coordinate; shape is (NBIN,).
        x (numpy.ndarray): (**Optional**) On-sky spectrum x coordinates;
            shape is (NSPEC,)
        y (numpy.ndarray): (**Optional**) On-sky spectrum y coordinates;
            shape is (NSPEC,).
        plot (bool): (**Optional**) Show the pPXF fit plot at each
            iteration.  Default is to skip the plot.
        quiet (bool): (**Optional**) Suppress output to the terminal
            (default).
        debug (bool): (**Optional**) Run in debugging mode.  Currently,
            all this does is perform the initial setup and then return
            empty vectors of the correct shape.  No fits are performed.
   
    Returns:
        numpy.ndarray: Eleven arrays are returned: (1) The best-fitting
        model for each spectrum [shape is (NSPEC,NPIX)]; (2) the
        best-fitting emission-line-only model for each spectrum [shape
        is (NSPEC,NPIX)]; (3) a boolean array that is True for all
        spectral pixels included in the fit [shape is (NSPEC,NPIX)]; (4)
        the best-fitting weight for each template in each spectrum
        [shape is (NSPEC,NTPL)]; (5) the error in the best-fitting
        template weights [shape is (NSPEC,NTPL)]; (6) the coefficients
        of the additive polynomial for each spectrum [shape is
        (NSPEC,DEGREE+1)], None if not fit; (7) the coefficients of the
        multiplicative polynomial for each spectrum [shape is
        (NSPEC,MDEGREE)]; (8) the best-fitting reddening values for each
        spectrum [shape is (NSPEC,)], None if not fit; (9) the input
        kinematics for each fit [shape is (NSPEC,sum(MOMENTS)], None if
        not fit; (10) the best-fit kinematics for each fit [shape is
        (NSPEC,sum(MOMENTS)]; (11) the formal error in the best-fit
        kinematics for each fit [shape is (NSPEC,sum(MOMENTS)].

    Raises:
        NotImplementedError: Raised if the number of stellar components
            is larger than 1.
        ValueError: Raised if: (1) some of the necessary bin-remapping
            data has not been provided (see function description); (2)
            the template flag object does not have the correct shape;
            (3) the wavelength and flux vectors do not match; and (4)
            the input list of kinematics does not have the correct
            length.
    """

    # Check that there is either one or zero stellar components
    if np.sum(gas_template) != templates.shape[0] and np.amax(inp_component[~gas_template]) != 0:
        raise NotImplementedError('Can only fit one stellar component!')

    # Confirm which binned datasets were provided (mask_binned can be
    # None)
    binned_data_provided = np.array([ x is not None for x in [flux_binned, noise_binned,
                                                              x_binned, y_binned, x, y] ])
    if np.any(binned_data_provided) and not np.all(binned_data_provided):
        raise ValueError('To use bin-remapping mode, must provide all of the following: '
                         'flux_binned, noise_binned, x_binned, y_binned, x, y')

    # Set the fitting mode
    mode = 'fitBins' if np.all(binned_data_provided) else 'noBins'

    # Check the shape of provided template flags
    if tpl_to_use is not None:
        if mode == 'fitBins' and tpl_to_use.shape[0] != flux_binned.shape[0]:
            raise ValueError('Template flag rows does not match number of binned spectra.')
        if mode == 'noBins' and tpl_to_use.shape[0] != flux.shape[0]:
            raise ValueError('Template flag rows does not match number of spectra.')
        if tpl_to_use.shape[1] != templates.shape[0]:
            raise ValueError('Template flags do not match the number of templates.')

    # Check the shape of the input flux, wavelength, and start objects
    nspec, nwave = flux.shape
    if mode == 'noBins' and len(inp_start) != nspec:
        raise ValueError('Input starting kinematic arrays do not match the input spectra.')
    if len(wave) != nwave:
        raise ValueError('Mismatch of wavelength vector with provided spectra.')

    # Get the number of templates and the number of kinematic parameters
    ntpl, npix_tpl = templates.shape
    nkin = np.sum(np.absolute(inp_moments))

    # Instantiate the output
    model_flux = np.zeros(flux.shape, dtype=float)
    model_eml_flux = np.zeros(flux.shape, dtype=float)
    model_mask = np.zeros(flux.shape, dtype=bool)

    tpl_wgt = np.zeros((nspec,ntpl), dtype=float)
    tpl_wgt_err = np.zeros((nspec,ntpl), dtype=float)

    addcoef = None if degree < 0 else np.zeros((nspec,degree+1), dtype=float)
    multcoef = None if mdegree < 1 else np.zeros((nspec,mdegree), dtype=float)
    ebv = None if reddening is None else np.zeros(nspec, dtype=float)

    kininp = np.zeros((nspec,nkin), dtype=float)
    kin = np.zeros((nspec,nkin), dtype=float)
    kin_err = np.zeros((nspec,nkin), dtype=float)

    # If debugging, just return the initialized output
    if debug:
        warnings.warn('JUST DEBUGGING.  NO EMISSION-LINE FITS PERFORMED!!')
        if mode == 'fitBins':
            # - Get the index of the nearest bin for every spaxel
            nearest_bin = np.argmin(np.square(x[:,None] - x_binned)
                                        + np.square(y[:,None] - y_binned), axis=1)
        else:
            nearest_bin = np.arange(nspec)
        kininp = np.array([np.concatenate(tuple(inp_start[0]))]*nspec)
        kin = kininp
        kinerr = kin/10
        return model_flux, model_eml_flux, model_mask, tpl_wgt, tpl_wgt_err, addcoef, multcoef, \
                    ebv, kininp, kin, kin_err, nearest_bin

    # First fit the binned data
    if mode == 'fitBins':

        # If binned data provided, the input starting kinematics arrays
        # must be sized appropriately for the binned data 
        nbin = flux_binned.shape[0]
        if len(inp_start) != nbin:
            raise ValueError('Input starting kinematic arrays do not match the binned spectra.')
        if len(wave) != flux_binned.shape[1]:
            raise ValueError('Mismatch of wavelength vector with provided binned spectra.')

        # First fit the binned data:
        # - Stellar components with fixed kinematics
        # - All gas templates in a single component
        component, moments, start = _ppxf_component_setup(inp_component, gas_template, inp_start,
                                                          single_gas_component=True)
        _, _, _, binned_tpl_wgts, _, _, _, _, _, binned_kin, _ \
                    = _fit_iteration(templates, wave, flux_binned, noise_binned, velscale, start,
                                     moments, component, gas_template, tpl_to_use=tpl_to_use,
                                     reject_boxcar=reject_boxcar, velscale_ratio=velscale_ratio,
                                     degree=degree, mdegree=mdegree, reddening=reddening,
                                     mask=mask_binned, vsyst=vsyst, plot=plot, quiet=quiet,
                                     sigma_rej=sigma_rej)

        valid_bin_str_fit = np.sum(binned_tpl_wgts[:,np.invert(gas_template)], axis=1) > 0
        valid_bin_fit = np.sum(binned_tpl_wgts, axis=1) > 0
        nvalid_bin_str_fit = np.sum(valid_bin_str_fit)

        # - Create new template set with the optimal stellar template
        #   for each spectrum, only using those with valid stellar
        #   weights
        stellar_wgts, _templates, _gas_template, _tpl_to_use, _component \
                    = _combine_stellar_templates(templates, gas_template,
                                                 binned_tpl_wgts[valid_bin_str_fit,:],
                                                 inp_component)

        # - Get the index of the nearest bin for every spaxel
        nearest_bin = np.argmin(np.square(x[:, None] - x_binned[valid_bin_str_fit]) 
                                    + np.square(y[:, None] - y_binned[valid_bin_str_fit]), axis=1)

        # - Use the binned data as the starting guess for the nearest
        #   spaxel
        stellar_wgts = stellar_wgts[nearest_bin,:]
        _templates = np.append(_templates[:nvalid_bin_str_fit,:][nearest_bin,:],
                                    _templates[nvalid_bin_str_fit:,:], axis=0)
        _gas_template = np.zeros(_templates.shape[0], dtype=bool)
        _gas_template[nspec:] = True

        _tpl_to_use = np.zeros((nspec,_templates.shape[0]), dtype=bool)
        _tpl_to_use[:,_gas_template] = True
        _tpl_to_use[np.arange(nspec),np.arange(nspec)] = True

        _component = np.append(np.zeros(nspec, dtype=int), _component[nvalid_bin_str_fit:])

        # - Use the starting positions from the fit to the bins for the
        #   fit to the individual spaxels
        n_gas_comp = len(np.unique(_component[_gas_template]))
        _start = inp_start[nearest_bin]
        gas_start = binned_kin[valid_bin_str_fit,np.absolute(moments[0]):][nearest_bin,:]
        for i in range(nspec):
            _start[i,1:] = np.array([ [gas_start[i].tolist()]*n_gas_comp ])

    else:
        nearest_bin = np.arange(nspec)
        stellar_wgts = None
        _templates = templates
        _gas_template = gas_template
        _tpl_to_use = tpl_to_use
        _component = inp_component
        _start = inp_start


    # Fit the data iteratively:
    #  - Fit with all the gas templates as part of one component
    component, moments, start = _ppxf_component_setup(_component, _gas_template, _start,
                                                      single_gas_component=True)
    model_mask = mask.copy()
    _, _, model_mask, tpl_wgts, _, _, _, _, _, kin, _ \
            = _fit_iteration(_templates, wave, flux, noise, velscale, start, moments, component,
                             _gas_template, tpl_to_use=_tpl_to_use, reject_boxcar=reject_boxcar,
                             velscale_ratio=velscale_ratio, degree=degree, mdegree=mdegree,
                             reddening=reddening, mask=model_mask, vsyst=vsyst, plot=plot,
                             quiet=quiet, sigma_rej=sigma_rej)

    if mode == 'noBins':
        # - If no previous fit to binned spectra, create new template
        #   set with the optimal stellar template for each spectrum
        # TODO: Template modes:
        #   - use all
        #   - use anything that was non-zero in the first fit
        #   - use single, optimal template
        stellar_wgts, _templates, _gas_template, _tpl_to_use, _component \
                    = _combine_stellar_templates(_templates, _gas_template, tpl_wgts, _component)

    # - Use this result to reset the starting estimates for the gas kinematics
    all_stellar_moments = np.sum(np.absolute(moments[:-1]))
    gas_start = kin[:,all_stellar_moments:]
    component, moments, start = _ppxf_component_setup(_component, _gas_template, _start,
                                                      gas_start=gas_start)

    # - Refit without rejection but with the tying in place
    model_flux, model_eml_flux, model_mask, tpl_wgts, tpl_wgts_err, addcoef, multcoef, ebv, \
            kininp, kin, kin_err = _fit_iteration(_templates, wave, flux, noise, velscale, start,
                                                  moments, component, _gas_template,
                                                  tpl_to_use=_tpl_to_use, reject_boxcar=None,
                                                  velscale_ratio=velscale_ratio, degree=degree,
                                                  mdegree=mdegree, reddening=reddening, tied=tied,
                                                  mask=model_mask, vsyst=vsyst, plot=plot,
                                                  quiet=quiet, sigma_rej=sigma_rej)

    # - Use the single output weight to renormalize the individual
    #   stellar template weights (only one of the weights for the
    #   non-gas templates should be non-zero); stellar-weight errors are
    #   always returned as 0.
    _tpl_wgts = stellar_wgts * np.sum(tpl_wgts[:,np.invert(_gas_template)], axis=1)[:,None]
    _tpl_wgts[:,gas_template] = tpl_wgts[:,_gas_template]
    _tpl_wgts_err = np.zeros((nspec,ntpl), dtype=float)
    _tpl_wgts_err[:,gas_template] = tpl_wgts_err[:,_gas_template]

    return model_flux, model_eml_flux, model_mask, _tpl_wgts, _tpl_wgts_err, addcoef, multcoef, \
                    ebv, kininp, kin, kin_err, nearest_bin

