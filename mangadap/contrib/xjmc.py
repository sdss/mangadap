# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Implements an emission-line fitting function using pPXF.

Revision history
----------------

    | **31 Aug 2017**: First commit by Xihan Ji (XJ)
    | **06 Feb 2018**: K. Westfall (KBW) added documentation
    | **09 Feb 2018**: (KBW) Return the bin matching vector
    | **20 Mar 2018**: (KBW) Correct error in carrying around pixels
        rejected during fit
    | **22 May 2018**: (KBW) Change import to ppxf package.
    | **29 May 2018**: (KBW) Remove original function from Xihan and
        rename ``*_edit`` function.

----

.. include license and copyright
.. include:: ../copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../links.rst
"""

import warnings
import numpy as np

from scipy import fftpack
from scipy.ndimage import rank_filter

import astropy.constants

from ppxf import ppxf, capfit

# For debugging
from matplotlib import pyplot as plt
# from mangadap.proc.ppxffit import PPXFModel

def ppxf_tied_parameters(component, vgrp, sgrp, moments):
    r"""
    Construct the object used to tie kinematic parameters in pPXF.

    .. note::
        - The components and kinematics groups can potentially have
          reduntant information.  E.g., If all sigma groups also have
          tied velocities, they'll be part of a component and will not
          required a tied parameter.

    .. warning::
        - high-order moments are not currently tied, but this is
          straight-forward to add.

    Args:
        component (array-like):
            The kinematic component assigned to each template.  Shape is
            :math:`(N_{\rm tpl},)`.
        vgrp (array-like):
            The velocity group assigned to each template.  Shape is
            :math:`(N_{\rm tpl},)`.  Can be None for no groups.
        sgrp (array-like):
            The velocity-dispersion (sigma) group assigned to each
            template.  Shape is :math:`(N_{\rm tpl},)`.  Can be None for
            no groups.
        moments (array-like):
            The number of moments assigned to each component.  Can be
            negative for fixed parameters.  Shape is :math:`(N_{\rm
            comp},)`.

    Returns:
        list: The tied object to pass to ppxf.  Will be None if both
        vgrp or sgrp are None or when the consolidation of the groups
        and components result in no tying being necessary.

    Raises:
        ValueError:
            Raised if the components or the groups are not an
            uninterupted sequence from 0..N-1, if the moments are
            provided for each component, or the length of the `vgrp` or
            `sgrp` arrays do not match the input component array.
    """
    if vgrp is None and sgrp is None:
        # Nothing to tie!
        return None

    # Check input
    ncomp = np.amax(component)+1
    if not np.array_equal(np.arange(ncomp), np.unique(component)):
        raise ValueError('Components must range from 0..N-1.')
    if len(moments) != ncomp:
        raise ValueError('Must provide the number of moments for each component.')
    ntpl = len(component)
    if vgrp is not None and not np.array_equal(np.arange(np.amax(vgrp)+1), np.unique(vgrp)):
        raise ValueError('Velocity groups must range from 0..N-1.')
    if vgrp is not None and len(vgrp) != ntpl:
        raise ValueError('Must provided a velocity group for each template.')
    if sgrp is not None and not np.array_equal(np.arange(np.amax(sgrp)+1), np.unique(sgrp)):
        raise ValueError('Sigma groups must range from 0..N-1.')
    if sgrp is not None and len(sgrp) != ntpl:
        raise ValueError('Must provided a sigma group for each template.')

    # Build the tied parameters as a vector
    tied = np.full(np.sum(np.absolute(moments)), '', dtype=object)
    tpli = np.arange(ntpl)
    for i in range(ncomp):
        # Do not allow tying to fixed components?
        if moments[i] < 0:
            continue

        # Velocity group of this component
        if vgrp is not None:
            indx = np.unique(component[tpli[vgrp == i]])
            if len(indx) > 1:
                parn = [ 0 + np.sum(np.absolute(moments[:j])) for j in indx ]
                tied[parn[1:]] = 'p[{0}]'.format(parn[0])
            
        # Sigma group of this component
        if sgrp is not None:
            indx = np.unique(component[tpli[sgrp == i]])
            if len(indx) > 1:
                parn = [ 1 + np.sum(np.absolute(moments[:j])) for j in indx ]
                tied[parn[1:]] = 'p[{0}]'.format(parn[0])

    # Check if anything is actually tied
    if np.array_equal(np.unique(tied), ['']):
        return None

    # Return after restructuring the vector into a list
    s = [np.sum(np.absolute(moments[:i])) for i in range(ncomp)]
    return [tied[s[i]:s[i]+np.absolute(moments[i])].tolist() for i in range(ncomp)]


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
        component (:obj:`numpy.ndarray`):
            The full list of components for all templates.  The total
            number of components is NCOMP.
        gas_template (:obj:`numpy.ndarray`):
            A boolean array identifying the gas templates.
        start (:obj:`list`, :obj:`numpy.ndarray`):
            The starting kinematics to use for each object spectrum.
            Shape is (NOBJ,); each element has shape (NCOMP,); each
            component element has shape (NMOM,), such that the number of
            moments can be different for each component.  NCOMP and NMOM
            should be the same for all object spectra.
        single_gas_component (:obj:`bool`, optional):
            Flag to force all gas components to have the same
            kinematics.
        gas_start (:obj:`list`, :obj:`numpy.ndarray`, optional):
            The starting kinematics to use for the gas components.
            Shape must be (NOBJ,2).

    Returns:
        Three numpy arrays are returned:
            - (1) The list of downselected and renumbered, if necessary,
              kinematic components, 
            - (2) the downselected and reordered, if necessary, number
              of moments to fit, and 
            - (3) the downselected and reordered starting guesses for
              each component.

    Raises:
        ValueError:
            Raised if the provided gas starting kinematics is not
            correct.
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
        _moments = np.append(_moments[_moments < 0], [2])

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

def _reset_components(c, valid):
    # Select valid components
    _c = c[valid]
    # Map between old and new component numbers
    c_map, inv = np.unique(_c, return_inverse=True)
    # Return the new component numbers and the old-to-new map
    return np.arange(len(c_map))[inv], c_map


def _validate_templates_components(templates, gas_template, component, vgrp, sgrp, moments, mask,
                                   start, tpl_to_use, velscale, velscale_ratio=None, vsyst=0):
    r"""
    Check that templates are valid in the sense that they have non-zero
    components in valid pixel regions.
    
    If all templates are valid, the returned data is identical to the
    input.  If not, the returned data are restructured as necessary for
    running pPXF.
    
    The start vector is used to set the guess offset (including vsyst)
    between the provided mask and the template.  The tpl_to_use vector
    is used to set a tempate as invalid.

    Args:
        templates (numpy.ndarray):
            Templates library to use for fitting.  Shape is
            (NTPLPIX,NTPL).
        gas_template (numpy.ndarray):
            Boolean vector that selects the gas templates.  Shape is
            (NTPL,).
        component (numpy.ndarray):
            Integer vector identifying the kinematic component for each
            template.  Shape is (NTPL,).
        vgrp (array-like):
            The integer velocity group associated with each template.
            Shape is :math:`(N_{\rm tpl},)` , but can be None.
        sgrp (array-like):
            The integer sigma group associated with each template.
            Shape is :math:`(N_{\rm tpl},)` , but can be None. 
        moments (numpy.ndarray):
            Integer vector with the number of kinematic moments for each
            component.  Shape is (NCOMP,).
        mask (numpy.ndarray):
            Boolean vector that selects the pixels in the object
            spectrum to fit (i.e., mask=True for pixels to fit and
            mask=False for pixels to ignore).  As in pPXF, the length is
            expected to be less than or equal to NTPLPIX in the template
            spectra.
        start (list):
            The starting kinematics for each kinematic component.
            Length is NCOMP.
        tpl_to_use (numpy.ndarray):
            In addition to removing unconstrained templates, this
            boolean vector selects that should even be considered in the
            fit.  Shape is (NTPL,).
        velscale (float):
            The pixel scale of the object spectrum to fit in km/s.
        velscale_ratio (:obj:`int`, optional):
            The ratio between the object and template pixel scale.
        vsyst (:obj:`float`, optional): 
            The pseudo velocity shift between the template and object
            spectra just due to the difference in the starting
            wavelength.

    Returns:
        tuple: Ten arrays are returned: 
            - (1) Boolean vector with which templates were valid (length
              is NTPL), 
            - (2) the valid set of templates to pass to pPXF [shape is
              (NTPLPIX,_NTPL)], 
            - (3) the boolean vector selecting gas templates (length is
              _NTPL),
            - (4) an integer vector mapping the input component number
              to the output component number (i.e., component_map[0] is
              the original component number for the downselected 0th
              component; length is _NCOMP), 
            - (5) the new component number for the downselected
              templates (length is _NTPL), 
            - (6) the new velocity group number for the downselected
              templates (length is _NTPL), 
            - (7) the new sigma group number for the downselected
              templates (length is _NTPL), 
            - (8) the number of moments for the new components (length
              is _NCOMP), 
            - (9) the starting kinematics for each new component (length
              is _NCOMP), and 
            - (10) the parameter tying object reordered as necessary for
              the new component list (**this needs to be checked**).

    Raises:
        ValueError:
            Raised if no templates are valid.

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
        return valid, templates, gas_template, np.arange(ncomp), component, \
                    vgrp, sgrp, moments, start

    # None (!) of the templates are valid
    if not np.any(valid):
        raise ValueError('No valid templates in fit!')

    # Templates to fit
    _templates = templates[valid,:]
    # Set which are gas templates
    _gas_template = gas_template[valid]

    # Reset components, moments, and starting kinematics
    _component, component_map = _reset_components(component, valid)
    _moments = moments[component_map]
    _start = [start[cm] for cm in component_map]
    # Reset velocity groups
    _vgrp = None if vgrp is None else _reset_components(vgrp, valid)[0]
    # Reset velocity dispersion groups
    _sgrp = None if sgrp is None else _reset_components(sgrp, valid)[0]

    return valid, _templates, _gas_template, component_map, _component, \
                _vgrp, _sgrp, _moments, _start


def _reorder_solution(ppsol, pperr, component_map, moments, start=None, fill_value=-999.):
    """
    Provided the best-fitting parameters from a ppxf instance,
    restructure the parameters to account for components that were
    removed during the template validation.

    Args:
        ppsol (list):
            The pPXF solution parameters
        pperr (list):
            The formal errors on the pPXF solution parameters
        component_map (numpy.ndarray):
            An integer vector mapping the input component number to the
            output component number; i.e., component_map[0] is the
            original component number for the downselected 0th
            component; length is _NCOMP.
        moments (list):
            The number of moments for the *original* set of components;
            length is NCOMP.
        start (:obj:`list`, optional): 
            The starting kinematics for the *original* set of
            components.  If provided, the starting values are included
            in the output parameter object for components that were not
            included in the pPPXF fit; otherwise, the unfit components
            are given placeholder parameter values.
        fill_value (:obj:`float`, optional):
            The placeholder value to give parameters and errors for
            components not included in the pPXF fit.

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
        indx = component_map == i
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
                   reddening=None, vgrp=None, sgrp=None, mask=None, vsyst=0, plot=False, quiet=True,
                   sigma_rej=3.):
    r"""
    Run a single fit+rejection iteration of the pPXF fit for all input
    spectra with the provided set of constraints/options.

    Args: 
        templates (numpy.ndarray):
            Full template library.  Shape is (NTPL,NTPLPIX).
        wave (numpy.ndarray):
            Wavelength vector.  Shape is (NPIX,).
        flux (numpy.ndarray):
            Object spectra to fit.  Shape is (NSPEC,NPIX).
        noise (numpy.ndarray):
            Error in the object spectra to fit.  Shape is (NSPEC,NPIX).
        velscale (float):
            The pixel scale of the object spectra in km/s.
        start (list):
            The starting kinematics for each kinematic component.
            Length is NCOMP.
        moments (numpy.ndarray):
            Integer vector with the number of kinematic moments for each
            component.  Shape is (NCOMP,).
        component (numpy.ndarray):
            Integer vector identifying the kinematic component for each
            template.  Shape is (NTPL,).
        gas_template (numpy.ndarray):
            Boolean vector that selects the gas templates.  Shape is
            (NTPL,).
        tpl_to_use (numpy.ndarray, optional):
            Boolean vector selecting templates to consider during the
            fit.  Shape is (NTPL,).  If None, all templates are used in
            the fit.
        reject_boxcar (:obj:`int`, optional):
            Size of the window for the rejection statistics.  Should be
            an odd number and larger and 10.  Default is 101.  If None,
            no rejection iteration is performed.
        velscale_ratio (:obj:`int`, optional):
            The ratio between the object and template pixel scale.
        degree (:obj:`int`, optional):
            Order of the additive polynomial to include in the fit.  Not
            included by default.
        mdegree (:obj:`int`, optional):
            Order of the multiplicative polynomial to fit.  Not included
            by default.
        reddening (:obj:`float`, optional):
            Initial E(B-V) guess for reddening (uses ppxf-default
            Calzetti 2000 model).  No attentuation fit by default.
        vgrp (array-like, optional):
            The integer velocity group associated with each template.
            Shape is :math:`(N_{\rm tpl},)`. 
        sgrp (array-like, optional):
            The integer sigma group associated with each template.
            Shape is :math:`(N_{\rm tpl},)`. 
        mask (numpy.ndarray, optional):
            Boolean vector that selects the pixels in the object spectra
            to fit (i.e., mask=True for pixels to fit and mask=False for
            pixels to ignore).  Shape is (NSPEC,NPIX).  All pixels are
            fit by default.
        vsyst (:obj:`float`, optional):
            The pseudo velocity shift between the template and object
            spectra just due to the difference in the starting
            wavelength.
        plot (:obj:`bool`, optional):
            Show the pPXF fit plot at each iteration.  Default is to
            skip the plot.
        quiet (:obj:`bool`, optional):
            Suppress output to the terminal (default).
        sigma_rej (:obj:`float`, optional):
            Sigma values used for the rejection.  Default is 3.
            
    Returns:
        tuple: Eleven arrays are returned:

            - (1) The best-fitting model for each spectrum [shape is
              (NSPEC,NPIX)];
            - (2) the best-fitting emission-line-only model for each
              spectrum [shape is (NSPEC,NPIX)]; 
            - (3) a boolean array that is True for all spectral pixels
              included in the fit [shape is (NSPEC,NPIX)]; 
            - (4) the best-fitting weight for each template in each
              spectrum [shape is (NSPEC,NTPL)]; 
            - (5) the error in the best-fitting template weights [shape
              is (NSPEC,NTPL)]; 
            - (6) the coefficients of the additive polynomial for each
              spectrum [shape is (NSPEC,DEGREE+1)]; -(7) the
              coefficients of the multiplicative polynomial for each
              spectrum [shape is (NSPEC,MDEGREE)];
            - (8) the best-fitting reddening values for each spectrum
              [shape is (NSPEC,)]; 
            - (9) the input kinematics for each fit [shape is
              (NSPEC,sum(MOMENTS)]; 
            - (10) the best-fit kinematics for each fit [shape is
              (NSPEC,sum(MOMENTS)]; 
            - (11) the formal error in the best-fit kinematics for each
              fit [shape is (NSPEC,sum(MOMENTS)].

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

        # Confirm that all templates and components are valid and
        # rearrange component arrays, if necessary
        valid_templates, _templates, _gas_template, component_map, _component, _vgrp, _sgrp, \
                _moments, _start \
                        = _validate_templates_components(templates, gas_template, component, vgrp,
                                                         sgrp, moments, model_mask[i,:], start[i], 
                                                         _tpl_to_use[i,:], velscale,
                                                         velscale_ratio=velscale_ratio,
                                                         vsyst=vsyst)

        # Construct the parameter tying structure
        tied = ppxf_tied_parameters(_component, _vgrp, _sgrp, _moments)

        # Run the first fit
        if plot:
            plt.clf()
        pp = ppxf.ppxf(_templates.T, flux[i,:], noise[i,:], velscale, _start,
                       velscale_ratio=velscale_ratio, plot=plot, moments=_moments, degree=degree,
                       mdegree=mdegree, lam=wave, reddening=reddening, tied=tied,
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

            # Confirm that all templates and components are valid and
            # rearrange component arrays, if necessary
            valid_templates, _templates, _gas_template, component_map, _component, _vgrp, _sgrp, \
                    _moments, _start \
                            = _validate_templates_components(templates, gas_template, component,
                                                             vgrp, sgrp, moments, model_mask[i,:],
                                                             sol, _tpl_to_use[i,:], velscale,
                                                             velscale_ratio=velscale_ratio,
                                                             vsyst=vsyst)

            # Construct the parameter tying structure
            tied = ppxf_tied_parameters(_component, _vgrp, _sgrp, _moments)

            # Refit using best-fit kinematics from previous fit as
            # initial guesses
            if plot:
                plt.clf()
            pp = ppxf.ppxf(_templates.T, flux[i,:], noise[i,:], velscale, _start,
                           velscale_ratio=velscale_ratio, plot=plot, moments=_moments,
                           degree=degree, mdegree=mdegree, lam=wave, reddening=reddening,
                           tied=tied, mask=model_mask[i,:], vsyst=vsyst, component=_component,
                           gas_component=_gas_template, quiet=quiet, linear=linear)
            if plot:
                plt.show()

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

#        print('nspec:{0}, chi2:{1}'.format(i+1, pp.chi2*pp.ndof))
        
    # Done
    print('Fitting spectrum: {0}/{0}'.format(nspec))

    return model, eml_model, model_mask, tpl_wgts, tpl_wgts_err, addcoef, multcoef, ebv, \
                    kininp, kin, kin_err


def _set_to_using_optimal_templates(nspec, ntpl, wgts, component, vgrp=None, sgrp=None):
    # Set the gas templates
    _gas_template = np.ones(ntpl, dtype=bool)
    _gas_template[:nspec] = False
    ngas = np.sum(_gas_template)
   
    # Set which templates to use with each spectrum; do not use optimal
    # templates where the stellar continuum was given no weight
    _tpl_to_use = np.zeros((nspec,ntpl), dtype=bool)
    _tpl_to_use[:,nspec:] = True
    _tpl_to_use[np.arange(nspec),np.arange(nspec)] = np.sum(wgts[:nspec,:], axis=1) > 0

    # Update the components
    _component = np.append(np.zeros(nspec, dtype=int), component[-ngas:])

    # Check that the components make sense
    if np.any(np.unique(_component) != np.arange(np.amax(_component)+1)):
        raise ValueError('Problem constructing new component vector.  Likely more than one '
                         'stellar component.')

    # Update the kinematic groups
    _vgrp = None if vgrp is None else np.append(np.zeros(nspec, dtype=int), vgrp[-ngas:])
    _sgrp = None if sgrp is None else np.append(np.zeros(nspec, dtype=int), sgrp[-ngas:])

    return _gas_template, _tpl_to_use, _component, _vgrp, _sgrp


def _combine_stellar_templates(templates, gas_template, wgts, component, vgrp=None, sgrp=None):
    r"""
    Reconstruct the set of templates used to fit each spectrum.

    Based on an initial fit to the spectra, construct a single optimal
    stellar template to use for each spectrum, and combine those with
    the existing gas templates.  Each optimal stellar-continuum template
    is set to the zeroth component.

    .. warning::

        - There can only be one stellar component per fit.
        - It is possible for some optimal stellar-continuum templates to
          be 0 everywhere.  In this case, the template is excluded from
          the array selecting which templates to include in the fit;
          however, that template *is* assigned a component.  These
          details are handled by :func:`_ppxf_component_setup`.
  
    .. todo::

        Allow for different template construct modes?  E.g., use all,
        use anything that was non-zero in the first fit, or use single
        optimal template (as done now).

    Args:
        templates (numpy.ndarray):
            Templates library to use for fitting.  Shape is
            (NTPLPIX,NTPL).
        gas_template (numpy.ndarray):
            Boolean vector that selects the gas templates.  Shape is
            (NTPL,).
        wgts (numpy.ndarray):
            The best-fitting template weights for each template in each
            spectrum.  Shape is (NSPEC,NTPL).
        component (numpy.ndarray):
            Integer vector identifying the kinematic component for each
            template.  Shape is (NTPL,).
        vgrp (array-like, optional):
            The integer velocity group associated with each template.
            Shape is :math:`(N_{\rm tpl},)` , but can be None.
        sgrp (array-like, optional):
            The integer sigma group associated with each template.
            Shape is :math:`(N_{\rm tpl},)` , but can be None. 

    Returns:
        tuple: Seven arrays are returned:

            - (1) the weights of only the stellar templates used in
              constructing each optimal template [shape is
              (NSPEC,NTPL)]; 
            - (2) the optimal stellar templates for each spectrum with
              the gas templates appended [shape is
              (NSPEC+NGASTPL,NPIXTPL)]; 
            - (3) a boolean array that selects the gas templates [shape
              is (NSPEC+NGASTPL,); 
            - (4) boolean array selecting which templates to use with
              each spectrum [shape is (NSPEC,NSPEC+NGASTPL)]
            - (5) integer array setting the component associated with
              each template [shape is (NSPEC+NGASTPL)].
            - (6) integer array setting the velocity group associated
              with each template [shape is (NSPEC+NGASTPL)].
            - (7) integer array setting the sigma group associated with
              each template [shape is (NSPEC+NGASTPL)].

    Raises:
        ValueError:
            Raised if there is more than one stellar component.

    """
    # Get the best-fit templates for each bin
    stellar_wgts = wgts.copy()
    stellar_wgts[:,gas_template] = 0.0
    optimal_template = np.dot(stellar_wgts, templates)

    # Update the list of templates
    _templates = np.append(np.atleast_2d(optimal_template), templates[gas_template,:], axis=0)

    # Update the relevant vectors
    _gas_template, _tpl_to_use, _component, _vgrp, _sgrp \
            = _set_to_using_optimal_templates(stellar_wgts.shape[0], _templates.shape[0],
                                              stellar_wgts, component, vgrp=vgrp, sgrp=sgrp)

    # Return new template list
    return stellar_wgts, _templates, _gas_template, _tpl_to_use, _component, _vgrp, _sgrp


def emline_fitter_with_ppxf(templates, wave, flux, noise, mask, velscale, velscale_ratio,
                            inp_component, gas_template, inp_moments, inp_start, vgrp=None,
                            sgrp=None, degree=-1, mdegree=0, reddening=None, reject_boxcar=101,
                            vsyst=0, tpl_to_use=None, binid=None, flux_binned=None,
                            noise_binned=None, mask_binned=None, x_binned=None, y_binned=None,
                            x=None, y=None, plot=False, quiet=False, debug=False, sigma_rej=3.):
    r"""
    Main calling function for fitting stellar-continuum and nebular
    emission lines in many spectra using pPXF.

    This is a generalization of the original function provided by Xihan
    Ji and Michele Cappellari.
    
    The function *does not* fit for the stellar kinematics; these are
    fixed during the fit (as provided in `inp_start`) and should have
    resulted from a previous fit to the stellar-continuum only; see
    :class:`mangadap.proc.ppxffit.PPXFFit`.  The number of stellar
    kinematic moments must have been the same for all spectra, and there
    can only be one stellar component.

    The templates are expected to be an integer number of
    `velscale_ratio` in length; see
    :func:`mangadap.proc.ppxffit.PPXFFit.check_templates`.

    The fitting procedure can be performed for a single set of spectra,
    or a set of spectra that are remapped from an input set of binned
    spectra.  The latter operation is selected by providing the binned
    flux, error, and mask arrays.  The individual spectra can be mapped
    to a binned spectrum using binned and unbinned on-sky coordinates or
    by providing the bin indices directly (see argument description
    below).

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

    To tie the kinematics assigned to each template, you have to assign
    them to velocity or velocity dispersion (sigma) groups using `vgrp`
    and `sgrp`, respectively.  These arrays are used to construct the
    ppxf `tied` argument, appropriate for each spectrum fit.

    Each template is identified as a gas template using the provided
    `gas_template` argument.

    .. todo::

        - Allow mask(s) to be optional
        - Update the docs

    Args:
        templates (numpy.ndarray):
            Templates library to use for fitting.  Shape is
            (NTPLPIX,NTPL).
        wave (numpy.ndarray):
            Wavelength vector.  Shape is (NPIX,).
        flux (numpy.ndarray):
            Object spectra to fit.  Shape is (NSPEC,NPIX).
        noise (numpy.ndarray):
            Error in the object spectra to fit.  Shape is (NSPEC,NPIX).
        mask (numpy.ndarray):
            Boolean vector that selects the pixels in the object spectra
            to fit (i.e., mask=True for pixels to fit and mask=False for
            pixels to ignore).  Shape is (NSPEC,NPIX).  All pixels are
            fit by default.
        velscale (float):
            The pixel scale of the object spectra in km/s.
        velscale_ratio (int):
            The ratio between the object and template pixel scale; must
            be an integer.
        inp_component (numpy.ndarray):
            Integer vector identifying the kinematic component for each
            template.  Shape is (NTPL,).
        gas_template (numpy.ndarray):
            Boolean vector that selects the gas templates.  Shape is
            (NTPL,).
        inp_moments (numpy.ndarray):
            Integer vector with the number of kinematic moments for each
            component.  Shape is (NCOMP,).
        inp_start (list, numpy.ndarray):
            The starting kinematics to use for each spectrum.  Shape is
            (NSPEC,); each element has shape (NCOMP,); each component
            element has shape (NMOM,), such that the number of moments
            can be different for each component.  NCOMP and NMOM should
            be the same for all object spectra.
        vgrp (array-like, optional):
            The integer velocity group associated with each template.
            Shape is :math:`(N_{\rm tpl},)`. 
        sgrp (array-like, optional):
            The integer sigma group associated with each template.
            Shape is :math:`(N_{\rm tpl},)`. 
        degree (:obj:`int`, optional):
            Order of the additive polynomial to include in the fit.  Not
            included by default.
        mdegree (:obj:`int`, optional):
            Order of the multiplicative polynomial to fit.  Not included
            by default.
        reddening (:obj:`float`, optional):
            Initial E(B-V) guess for reddening (uses ppxf-default
            Calzetti 2000 model).  No attentuation fit by default.
        reject_boxcar (:obj:`int`, optional):
            Size of the window for the rejection statistics.  Should be
            an odd number and larger and 10.  Default is 101.  If None,
            no rejection iterations are performed.
        vsyst (:obj:`float`, optional):
            The pseudo velocity shift between the template and object
            spectra just due to the difference in the starting
            wavelength. Default is 0 km/s.
        tpl_to_use (numpy.ndarray, optional):
            Boolean vector selecting templates to consider during the
            fit.  If None, all templates are used in the fit.  If
            provided, the shape must be (NBIN,NTPL) when providing the
            binned data and (NSPEC,NTPL) when no binned data are
            provided.
        binid (numpy.ndarray, optional):
            Bin index associated with each spectrum.  Ignored if binned
            spectra are not provided.  If binned spectra are provided,
            and this is None, coordinates must be provided and they are
            used to associate each bin with a binned spectrum just based
            by proximity.  If provided, this associates each spectrum to
            the binned spectrum to use for the stellar kinematics.
            Shape is (NSPEC,).
        flux_binned (numpy.ndarray, optional):
            Binned spectra with previous fits to the stellar kinematics.
            See purpose above.
        noise_binned (numpy.ndarray, optional):
            Error in the binned spectra.
        mask_binned (numpy.ndarray, optional):
            Mask for the binned spectra.
        x_binned (numpy.ndarray, optional):
            On-sky bin x coordinate; shape is (NBIN,).
        y_binned (numpy.ndarray, optional):
            On-sky bin y coordinate; shape is (NBIN,).
        x (numpy.ndarray, optional):
            On-sky spectrum x coordinates; shape is (NSPEC,)
        y (numpy.ndarray, optional):
            On-sky spectrum y coordinates; shape is (NSPEC,).
        plot (:obj:`bool`, optional):
            Show the pPXF fit plot at each iteration.  Default is to
            skip the plot.
        quiet (:obj:`bool`, optional):
            Suppress output to the terminal (default).
        debug (:obj:`bool`, optional):
            Run in debugging mode.  Currently, all this does is perform
            the initial setup and then return empty vectors of the
            correct shape.  No fits are performed.
   
    Returns:
        Eleven arrays are returned:
            - (1) The best-fitting model for each spectrum [shape is
              (NSPEC,NPIX)]; 
            - (2) the best-fitting emission-line-only model for each
              spectrum [shape is (NSPEC,NPIX)]; 
            - (3) a boolean array that is True for all spectral pixels
              included in the fit [shape is (NSPEC,NPIX)]; 
            - (4) the best-fitting weight for each template in each
              spectrum [shape is (NSPEC,NTPL)]; 
            - (5) the error in the best-fitting template weights [shape
              is (NSPEC,NTPL)]; 
            - (6) the coefficients of the additive polynomial for each
              spectrum [shape is (NSPEC,DEGREE+1)], None if not fit; 
            - (7) the coefficients of the multiplicative polynomial for
              each spectrum [shape is (NSPEC,MDEGREE)]; 
            - (8) the best-fitting reddening values for each spectrum
              [shape is (NSPEC,)], None if not fit; 
            - (9) the input kinematics for each fit [shape is
              (NSPEC,sum(MOMENTS)], None if not fit; 
            - (10) the best-fit kinematics for each fit [shape is
              (NSPEC,sum(MOMENTS)]; 
            - (11) the formal error in the best-fit kinematics for each
              fit [shape is (NSPEC,sum(MOMENTS)].

    Raises:
        NotImplementedError:
            Raised if the number of stellar components is larger than 1.
        ValueError:
            Raised if: (1) only some of the necessary bin-remapping data
            has not been provided (see function description); (2) the
            template flag object does not have the correct shape; (3)
            the wavelength and flux vectors do not match; (4) any of the
            spectra are fully masked (i.e., there's nothing to fit); or
            (5) the input list of kinematics does not have the correct
            length.
    """

    # Check that there is either one or zero stellar components
    if np.sum(gas_template) != templates.shape[0] and np.amax(inp_component[~gas_template]) != 0:
        raise NotImplementedError('Can only fit one stellar component!')

    # There must be data for all spectra to proceed
    # TODO: Mask cannot be None
    if np.any(np.invert(np.any(mask, axis=1))):
        raise ValueError('All spectra must have at least some pixels to fit!')

    # Confirm necessary binned data provided, if any is provided
    binned_spectra_provided = flux_binned is not None and noise_binned is not None
    can_match_binned_spectra = binid is not None \
                                or np.all([a is not None for a in [x_binned, y_binned, x, y]])
    if binned_spectra_provided and not can_match_binned_spectra:
        raise ValueError('Cannot match binned and unbinned spectra; provide the ID matches '
                         'directly or provide coordinates for a proximity match.')
    if not binned_spectra_provided and can_match_binned_spectra:
        warnings.warn('Binned spectra not (or incorrectly) provided; ignoring some arguments.')

    # If binned data is provided, get the binning match
    mode = 'noBins'         # Fitting mode
    _binid = None
    if binned_spectra_provided and can_match_binned_spectra:
        # Set fitting mode
        mode = 'fitBins'
        # Associate individual and binned spectra
        _binid = np.argmin(np.square(x[:,None]-x_binned) + np.square(y[:,None]-y_binned), axis=1) \
                            if binid is None else binid
        # Determine which individual spectra are components of binned
        # spectra actually made up of more than one spectrum.
        # TODO: Only do this if binid is provided directly?
        uniq, inv, cnt = np.unique(_binid, return_inverse=True, return_counts=True)
        component_of_bin = cnt[inv] > 1 

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
    if mode == 'fitBins' and len(_binid) != nspec:
        raise ValueError('Length if array matching spectra to bins is incorrect.')
    if len(wave) != nwave:
        raise ValueError('Mismatch of wavelength vector with provided spectra.')

    # Get the number of templates and the number of kinematic parameters
    ntpl, npix_tpl = templates.shape
    nkin = np.sum(np.absolute(inp_moments))

    # Check that the number of components is correct
    if len(inp_component) != ntpl:
        raise ValueError('Must assign each template to a kinematic component.')
    if vgrp is not None and len(vgrp) != ntpl:
        raise ValueError('If provided, must assign each template to a velocity group.')
    if sgrp is not None and len(sgrp) != ntpl:
        raise ValueError('If provided, must assign each template to a sigma group.')

    # Instantiate the output
    model_flux = np.zeros(flux.shape, dtype=float)
    model_eml_flux = np.zeros(flux.shape, dtype=float)
#    model_mask = np.zeros(flux.shape, dtype=bool)
    model_mask = mask.copy()

    tpl_wgts = np.zeros((nspec,ntpl), dtype=float)
    tpl_wgts_err = np.zeros((nspec,ntpl), dtype=float)

    addcoef = None if degree < 0 else np.zeros((nspec,degree+1), dtype=float)
    multcoef = None if mdegree < 1 else np.zeros((nspec,mdegree), dtype=float)
    ebv = None if reddening is None else np.zeros(nspec, dtype=float)

    kininp = np.zeros((nspec,nkin), dtype=float)
    kin = np.zeros((nspec,nkin), dtype=float)
    kin_err = np.zeros((nspec,nkin), dtype=float)

    # If debugging, just return the initialized output
    if debug:
        warnings.warn('JUST DEBUGGING.  NO EMISSION-LINE FITS PERFORMED!!')
        kininp = np.array([np.concatenate(tuple(inp_start[0]))]*nspec)
        kin = kininp
        kinerr = kin/10
        return model_flux, model_eml_flux, model_mask, tpl_wgt, tpl_wgt_err, addcoef, multcoef, \
                    ebv, kininp, kin, kin_err, _binid #nearest_bin

    # Fit the binned data
    if mode == 'fitBins':

        # For the binned data, the input starting kinematics arrays
        # must have the appropriate shape
        nbin = flux_binned.shape[0]
        if len(inp_start) != nbin:
            raise ValueError('Input starting kinematic arrays do not match the binned spectra.')
        if len(wave) != flux_binned.shape[1]:
            raise ValueError('Mismatch of wavelength vector with provided binned spectra.')

        # First fit iteration:
        # - Stellar components with fixed kinematics
        # - All gas templates in a single component (no parameter tying)
        component, moments, start = _ppxf_component_setup(inp_component, gas_template, inp_start,
                                                          single_gas_component=True)
        _, _, _mask, binned_tpl_wgts, _, _, _, _, _, binned_kin, _ \
                    = _fit_iteration(templates, wave, flux_binned, noise_binned, velscale, start,
                                     moments, component, gas_template, tpl_to_use=tpl_to_use,
                                     reject_boxcar=reject_boxcar, velscale_ratio=velscale_ratio,
                                     degree=degree, mdegree=mdegree, reddening=reddening,
                                     mask=mask_binned, vsyst=vsyst, plot=plot, quiet=quiet,
                                     sigma_rej=sigma_rej)

#        # - Determine which bins have applied non-zero weights to any
#        #   stellar template
#        valid_bin_str_fit = np.sum(binned_tpl_wgts[:,np.invert(gas_template)], axis=1) > 0
#        nvalid_bin_str_fit = np.sum(valid_bin_str_fit)

        # - Create a new template set that includes the optimal stellar
        #   template for each *binned* spectrum based on the
        #   best-fitting weights from above.  Propagate the changes to
        #   the components and groups.
        stellar_wgts, _templates, _gas_template, _tpl_to_use, _component, _vgrp, _sgrp \
                    = _combine_stellar_templates(templates, gas_template,
#                                                 binned_tpl_wgts[valid_bin_str_fit,:],
                                                 binned_tpl_wgts, inp_component, vgrp, sgrp)

#        # - Get the index of the nearest bin for every spaxel
#        nearest_bin = np.argmin(np.square(x[:, None] - x_binned[valid_bin_str_fit]) 
#                                    + np.square(y[:, None] - y_binned[valid_bin_str_fit]), axis=1)

        # - Use the binned data as the starting guess for the spaxels
#        stellar_wgts = stellar_wgts[nearest_bin,:]
#        _templates = np.append(_templates[:nvalid_bin_str_fit,:][nearest_bin,:],
#                                    _templates[nvalid_bin_str_fit:,:], axis=0)

        # - Restructure the relevant arrays to prepare for the fits to
        #   the *individual* spectra.  Propagate the changes to the
        #   components and groups.
        stellar_wgts = stellar_wgts[_binid,:]
        _templates = np.append(_templates[:nbin,:][_binid,:], _templates[nbin:,:], axis=0)
        _gas_template, _tpl_to_use, _component, _vgrp, _sgrp \
                = _set_to_using_optimal_templates(nspec, _templates.shape[0], stellar_wgts,
                                                  _component, vgrp=_vgrp, sgrp=_sgrp)
        tpl_wgts = np.ones((nspec, nspec+np.sum(_gas_template)), dtype=float)

        # - For bins that are also individual spaxels, update the mask
        # to include the rejected pixels
        if not np.all(component_of_bin):
            single_spaxel_bin = np.invert(component_of_bin)
            model_mask[single_spaxel_bin,:] &= _mask[_binid[single_spaxel_bin],:]

        # - Use the best-fit parameters from the binned spectra as the
        #   starting guesses for the fits to the the individual spaxels
        n_gas_comp = len(np.unique(_component[_gas_template]))
        _start = inp_start[_binid]
        gas_start = binned_kin[:,np.absolute(moments[0]):][_binid,:]
        for i in range(nspec):
            _start[i,1:] = np.array([ [gas_start[i].tolist()]*n_gas_comp ])
    else:
        # Binned spectra were not provided, prepare to fit the
        # individual spectra by just pointing to the input
        _binid = np.arange(nspec)
        component_of_bin = np.ones(nspec, dtype=bool)
        stellar_wgts = None
        _templates = templates
        _gas_template = gas_template
        _tpl_to_use = tpl_to_use
        _component = inp_component
        _vgrp = vgrp
        _sgrp = sgrp
        _start = inp_start

    # First fit to the individual spectra
    #  - Fit with all the gas templates as part of one component (no
    #    parameter tying).  When binned spectra have been fit
    #    previously, this only fits spectra that are components of a bin
    #    with more than one spectrum!

    component, moments, start = _ppxf_component_setup(_component, _gas_template, _start,
                                                      single_gas_component=True)

    # Restructure the kinematics array so that it can accept only
    # fitting the individual spectra that themselves consisted of an
    # entire bin
    kin = start.copy()
    kin = np.array([ k for k in kin ]).reshape(nspec,-1)

    _, _, model_mask[component_of_bin,:], tpl_wgts[component_of_bin,:], _, _, _, _, _, \
            kin[component_of_bin,:], _ \
                    = _fit_iteration(_templates, wave, flux[component_of_bin,:],
                                     noise[component_of_bin,:], velscale, start[component_of_bin],
                                     moments, component, _gas_template,
                                     tpl_to_use=_tpl_to_use[component_of_bin,:],
                                     reject_boxcar=reject_boxcar, velscale_ratio=velscale_ratio,
                                     degree=degree, mdegree=mdegree, reddening=reddening,
                                     mask=model_mask[component_of_bin,:], vsyst=vsyst, plot=plot,
                                     quiet=quiet, sigma_rej=sigma_rej)

    if mode == 'noBins':
        # If binned spectra were not provided, the fit above is the
        # first fit to the individual spectra using all
        # stellar-continuum templates.  Use the result of that fit to
        # construct the optimal template for each spectrum, and
        # restructure the relevant component arrays.
        stellar_wgts, _templates, _gas_template, _tpl_to_use, _component, _vgrp, _sgrp \
                    = _combine_stellar_templates(_templates, _gas_template, tpl_wgts, _component,
                                                 _vgrp, _sgrp)

    # Second fit to the individual spectra
    # - Use first fit to to reset the starting estimates for the gas kinematics
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
                                                  mdegree=mdegree, reddening=reddening,
                                                  vgrp=_vgrp, sgrp=_sgrp, mask=model_mask,
                                                  vsyst=vsyst, plot=plot, quiet=quiet,
                                                  sigma_rej=sigma_rej)

    # - Use the single output weight to renormalize the individual
    #   stellar template weights (only one of the weights for the
    #   non-gas templates should be non-zero); stellar-weight errors are
    #   always returned as 0.
    _tpl_wgts = stellar_wgts * np.sum(tpl_wgts[:,np.invert(_gas_template)], axis=1)[:,None]
    _tpl_wgts[:,gas_template] = tpl_wgts[:,_gas_template]
    _tpl_wgts_err = np.zeros((nspec,ntpl), dtype=float)
    _tpl_wgts_err[:,gas_template] = tpl_wgts_err[:,_gas_template]

    return model_flux, model_eml_flux, model_mask, _tpl_wgts, _tpl_wgts_err, addcoef, multcoef, \
                    ebv, kininp, kin, kin_err, _binid

