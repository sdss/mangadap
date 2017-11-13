#!/usr/bin/env python

import os
import warnings
import numpy as np
import matplotlib.pyplot as plt
import astropy.constants
import copy

from scipy import fftpack
from scipy.ndimage import rank_filter

from captools.ppxf import ppxf, _losvd_rfft, rebin
from captools import capfit

## Function used to calculate 1 sigma errors of the residuals
## Currently the pixels in the same window are assigned with
## the same error
#def calculate_noise(residuals, width=140):
#
#    # Number of windows
#    n_win = (residuals.shape[0]//width)
#
#    # Assign each window with equal number of pixels
#    cut = n_win*width
#    resid_temp = residuals[:cut].reshape(n_win, -1)
#
#    # Calculate noise (standard deviation) for each window
#    stdd = np.std(resid_temp, axis=1, ddof=1)
#    stdd = stdd.reshape(n_win, -1)
#
#    # Create the noise array with the same size as galaxy array
#    # There are some pixels uncovered by any window, which are 
#    # assigned with the same sigma as the nearest window
#    noise = np.repeat(stdd, width, axis=1).ravel()
#    uncover = np.repeat(stdd[-1], residuals.shape[0]-cut)
#    noise = np.append(noise, uncover)
#    
#    return noise

# Created by Michele Cappellari (23 September 2017)
def calculate_noise(residuals, width=101):

    """
    Robust determination of the error spectrum as 1/2 of the
    interval enclosing 68% of the values in a given running window

    """
    if width%2 != 1:
        raise ValueError('Must provided odd number width.')
    if width < 10:
        raise ValueError('Width must be 11 or higher.')
    
#    assert width%2 == 1, "width must be odd"
#    assert width > 10, "width %d too small"%width

    p = (1 - 0.6827)/2  # 1sigma
    lo, up = round(width*p), round(width*(1 - p))

    upper = rank_filter(residuals, rank=up, size=width)
    lower = rank_filter(residuals, rank=lo, size=width)
    noise = (upper - lower)/2

    return noise


# Stellar components must all have numbers less than any gas component;
# generalized to n stellar components (for what it's worth)
def _ppxf_component_setup(component, gas_template, start, single_gas_component=False,
                          gas_start=None):

    nobj = len(start)
#    print('nobj: ', nobj)

    if gas_start is not None and np.asarray(gas_start).shape != (nobj,2):
        raise ValueError('Provided gas starting kinematics has incorrect shape.')

#    print('input component: ', component)

    # Get the gas component numbers
    gas_components = np.unique(component[gas_template])

#    print('gas components (from components): ', gas_components)

    # Component assignments for each template
    _component = component.copy()
    n_gas_components = 1 if single_gas_component else len(gas_components)
    if n_gas_components == 1:
        _component[gas_template] = np.amax(component[~gas_template])+1

#    print('n gas components for output', n_gas_components)
#    print('output component: ', _component)

    # start shape is (nspec,); each element of start has shape (ncomp,);
    # each component has shape (nmom,), which can be different for each
    # component
    stellar_components = np.arange(np.amin(gas_components))
#    print('stellar components: ', stellar_components)
    _moments = np.array([ len(s) for s in start[0] ])
    _moments[stellar_components] *= -1
    if n_gas_components == 1:
        # Force a single gas component
        _moments = np.append( _moments[_moments < 0], [ 2 ])
#    print('output moments: ', _moments)

#    print('type gas_start: ', type(gas_start))

    _start = np.empty(len(start), dtype=object)
    for i,s in enumerate(start):
#        print(s)
        _start[i] = s[stellar_components].tolist() if len(stellar_components) > 0 else []
        if gas_start is None:
            _gas_start = [np.mean(np.asarray(s[gas_components]), axis=0).tolist()] \
                                if n_gas_components == 1 else [s[gas_components]]
        else:
            _gas_start = [gas_start[i].tolist()]*n_gas_components
        _start[i] += _gas_start
#        print(_start[i])
#    print(start[0:5])
#    print(_start[0:5])
        
    return _component, _moments, _start


def _ppxf_fit_matrix(start, npix_obj, templates_rfft, npix_tpl, npad, component, vsyst, moments,
                     velscale, velscale_ratio):
    """
    Mimic the first bit of ppxf._linear_fit()

    Limited as follows:
        - no additive or multiplicative polynomial
        - no reddening
        - only single galaxy spectrum (ppxf.nspec = 1)
        - no sigma_diff
    """
    # Template and component shape
    ntpl, nl = templates_rfft.shape
    if len(component) != ntpl:
        raise ValueError('Component length must match the number of templates.')
    ncomp = np.max(component)+1

    # Convert from km/s to pixel coordinates
    _start = np.array(start) #, dtype=object)
    for j in range(ncomp):
        _start[j][:2] /= velscale
    _vsyst = vsyst/velscale

    # Parameter vector:
    # pars = [vel_1, sigma_1, h3_1, h4_1, ... # Velocities are in pixels.
    #         ...                             # For all kinematic components
    #         vel_n, sigma_n, h3_n, h4_n, ...
    pars = _start.ravel()

    # FFT of the LOSVD
    losvd_rfft = _losvd_rfft(pars, 1, np.abs(moments), nl, ncomp, _vsyst, velscale_ratio, 0.)

    # Construct and return the convolved templates
    c = np.zeros((ntpl, npix_obj), dtype=float)
    for i, tpl_rfft in enumerate(templates_rfft):
        pr = tpl_rfft * losvd_rfft[:, component[i], 0]
        tt = np.fft.irfft(pr, npad)
        c[i,:] = tt[:npix_obj] if velscale_ratio == 1 else \
                    rebin(tt[:npix_tpl*velscale_ratio], velscale_ratio)[:npix_obj]
    return c


def _validate_templates_components(templates, gas_template, component, moments, mask, start, tied,
                                   tpl_to_use, velscale, velscale_ratio=None, vsyst=0):
    """
    Return the validated templates and reordered components.
    
    Check that templates are valid in the sense that they have non-zero
    components in valid pixel regions.  The start vector is used to set
    the guess offset (including vsyst) between the provided mask and the
    template.  The tpl_to_use vector is used to set a tempate as
    invalid.
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
    print(component[valid])
    valid &= tpl_to_use
    print(component[valid])
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
#        return valid, None, None, None, None, None, None, None

    # Templates to fit
    _templates = templates[valid,:]
    # Set which are gas templates
    _gas_template = gas_template[valid]
    # Sequential template numbers
    _component = component[valid]
    # Map between old and new component numbers
    component_map = np.unique(_component)

    print(component)
    print(_component)
    print(component_map)

    if not np.array_equal(component_map, np.arange(ncomp)):
        # Subset of components are not in sequence so they need to be remapped
        remapped_component = _component.copy()
        # Starting kinematics for each component
        _start = []
        _tied = None if tied is None else []
        for i, cm in enumerate(component_map):
            indx = _component == cm
            print(i)
            print(np.sum(indx))
            print(cm)
            print(start[cm])
            remapped_component[indx] = i
            _start += [start[cm]]
            if tied is not None:
                _tied += [tied[cm]]
        _component = remapped_component
        _moments = moments[component_map]
        print(_component)
        print(_moments)
    else:
        _start = start
        _tied = tied
        _moments = moments

    return valid, _templates, _gas_template, component_map, _component, _moments, _start, _tied


def _reorder_solution(ppsol, pperr, component_map, moments, start=None, fill_value=-999.):
    ncomp = len(moments)
    if np.array_equal(component_map, np.arange(ncomp)):
        return ppsol, pperr

    _ncomp = len(component_map)

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


# Run a single fit+rejection iteration
def _fit_iteration(templates, wave, flux, noise, velscale, start, moments, component, gas_template,
                   tpl_to_use=None, reject_boxcar=101, velscale_ratio=None, degree=-1, mdegree=0,
                   reddening=None, tied=None, mask=None, vsyst=0, plot=False, quiet=True):
    """
    mask should be True for the pixels to fit
    """

#        if templates_rfft is None:
#            # Pre-compute FFT of real input of all templates
#            self.templates_rfft = np.fft.rfft(self.star, self.npad, axis=0)
#        else:
#            self.templates_rfft = templates_rfft
#

#    med_flux = np.median(flux, axis=1)
#    _flux = flux / med_flux[:,None]
#    _noise = noise / med_flux[:,None]

    nspec = flux.shape[0]
    ntpl = templates.shape[0]
    nkin = np.sum(np.absolute(moments))

    _tpl_to_use = np.ones((nspec,ntpl), dtype=bool) if tpl_to_use is None else tpl_to_use.copy()
    print(_tpl_to_use.shape)
    print(np.sum(_tpl_to_use, axis=1))
    print(_tpl_to_use[0,:])

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

#    _start = start.copy()

    # For debugging
    linear=False
#    linear=True
#    reject_boxcar=None
#    mdegree=0

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
        pp = ppxf(_templates.T, flux[i,:], noise[i,:], velscale, _start,
                  velscale_ratio=velscale_ratio, plot=plot, moments=_moments, degree=degree,
                  mdegree=mdegree, lam=wave, reddening=reddening, tied=_tied, mask=model_mask[i,:],
                  vsyst=vsyst, component=_component, gas_component=_gas_template, quiet=quiet,
                  linear=linear)
        if plot:
            plt.show()

        # Reject 3-sigma outliers and refit, if requested by a provided
        # boxcar width
        if reject_boxcar is not None:
            # - Calculate residuals
            resid = flux[i,:] - pp.bestfit
            # - Select pixels includes in the fit and not fit by
            # emission lines
            reject_pixels = list(set(pp.goodpixels)
                                    & set(np.arange(len(resid))[pp.gas_bestfit < 1e-6]))
            # - Calculate the 1-sigma confidence interval
            NOISE = calculate_noise(resid[reject_pixels], width=reject_boxcar)
            # - Reject pixels with > 3-sigma residuals
            model_mask[i,reject_pixels] &= (abs(resid[reject_pixels]) < (3*NOISE))
        
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
            pp = ppxf(_templates.T, flux[i,:], noise[i,:], velscale, _start,
                      velscale_ratio=velscale_ratio, plot=plot, moments=_moments, degree=degree,
                      mdegree=mdegree, lam=wave, reddening=reddening, tied=_tied,
                      mask=model_mask[i,:], vsyst=vsyst, component=_component,
                      gas_component=_gas_template, quiet=quiet, linear=linear)
            if plot:
                plt.show()

        # Reorder the output; sets any omitted components to a default
        # value of -999.
        sol, err = _reorder_solution(pp.sol, pp.error, component_map, moments)

        # Save the results
        model[i,:] = pp.bestfit
        eml_model[i,:] = np.dot(pp.matrix[:,degree+1:][:,_gas_template], pp.weights[_gas_template])

#        plt.plot(flux[i,:], color='b', lw=2)
#        plt.plot(np.ma.MaskedArray(flux[i,:], mask=np.invert(model_mask[i,:])), color='k', lw=2)
#        plt.plot(model[i,:], color='C3', lw=1)
#        plt.plot(model[i,:]-eml_model[i,:], color='C2', lw=1)
#        plt.show()

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
        
    print('Fitting spectrum: {0}/{0}'.format(nspec))

    return model, eml_model, model_mask, tpl_wgts, tpl_wgts_err, addcoef, multcoef, ebv, \
                    kininp, kin, kin_err


def _combine_stellar_templates(templates, gas_template, wgts, component):

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
        raise ValueError('Problem constructing new component vector.')

    # Return new template list
    return stellar_wgts, _templates, _gas_template, _tpl_to_use, _component


# Requirements:
#   - All spectra are expected to have been fit with the same number of
#     stellar moments
#   - Only one component for the non-gas templates
#   - If tpl_to_use is provided, the shape must be (Nbin,Ntpl) when
#     providing the binned data and (Nspec,Ntpl) when no binned data are
#     provided.
#   - The templates are expected to be an integer number of
#     velscale_ratio in length; see PPXFFit.check_templates()

def emline_fitter_with_ppxf_edit(templates, wave, flux, noise, mask, velscale, velscale_ratio,
                                 inp_component, gas_template, inp_moments, inp_start,
                                 tied=None, degree=-1, mdegree=0, reddening=None,
                                 reject_boxcar=101, vsyst=0, tpl_to_use=None, flux_binned=None,
                                 noise_binned=None, mask_binned=None, x_binned=None, y_binned=None,
                                 x=None, y=None, plot=False, quiet=False, debug=False):

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

    # Instantiate the output
    nspec, nwave = flux.shape
    if len(wave) != nwave:
        raise ValueError('Mismatch of wavelength vector with provided spectra.')
#    masked_fraction = np.sum(mask, axis=1)/nwave
#    indx = masked_fraction > 1/3
#    model_binid = np.full(nspec, -1, dtype=int)
#    model_binid[indx] = np.arange(np.sum(indx))

    model_flux = np.zeros(flux.shape, dtype=float)
    model_eml_flux = np.zeros(flux.shape, dtype=float)
    model_mask = np.zeros(flux.shape, dtype=bool)

    ntpl, npix_tpl = templates.shape
    nkin = np.sum(np.absolute(inp_moments))

    tpl_wgt = np.zeros((nspec,ntpl), dtype=float)
    tpl_wgt_err = np.zeros((nspec,ntpl), dtype=float)

    addcoef = None if degree < 0 else np.zeros((nspec,degree+1), dtype=float)
    multcoef = None if mdegree < 1 else np.zeros((nspec,mdegree), dtype=float)
    ebv = None if reddening is None else np.zeros(nspec, dtype=float)

    kin = np.zeros((nspec,nkin), dtype=float)
    kin_err = np.zeros((nspec,nkin), dtype=float)

    if debug:
        warnings.warn('JUST DEBUGGING.  NO EMISSION-LINE FITS PERFORMED!!')
        return model_flux, model_eml_flux, model_mask, tpl_wgt, tpl_wgt_err, addcoef, multcoef, \
                    ebv, kin, kin_err

    # First fit the binned data
    if mode == 'fitBins':

        # If provided the binned data, the input starting kinematics
        # arrays must be sized appropriately for the binned data 
        nbin = flux_binned.shape[0]
        if len(inp_start) != nbin:
            raise ValueError('Input starting kinematic arrays do not match the binned spectra.')
        if len(wave) != flux_binned.shape[1]:
            raise ValueError('Mismatch of wavelength vector with provided binned spectra.')


        # First fit the binned data:
        # - Stellar components with fixed kinematics
        # - All gas templates in a single component
#        # DEBUG #####
#        binkin = np.zeros((nbin,nkin), dtype=float)
#        for i in range(nbin):
#            binkin[i,:] = np.concatenate(tuple(inp_start[i]))
#        print(binkin[0,:])
#        all_stellar_moments = np.sum(np.absolute(inp_moments[:-1]))
#        gas_start = binkin[:,all_stellar_moments:]
#        component, moments, start = _ppxf_component_setup(inp_component, gas_template, inp_start,
#                                                          gas_start=gas_start)
#        print(component)
#        print(moments)
#        # END DEBUG #
        component, moments, start = _ppxf_component_setup(inp_component, gas_template, inp_start,
                                                          single_gas_component=True)
        _, _, _, binned_tpl_wgts, _, _, _, _, _, binned_kin, _ \
                    = _fit_iteration(templates, wave, flux_binned, noise_binned, velscale, start,
                                     moments, component, gas_template, tpl_to_use=tpl_to_use,
                                     reject_boxcar=reject_boxcar, velscale_ratio=velscale_ratio,
                                     degree=degree, mdegree=mdegree, reddening=reddening,
                                     mask=mask_binned, vsyst=vsyst, plot=plot, quiet=quiet)

        # - Create new template set with the optimal stellar template
        #   for each spectrum
        stellar_wgts, _templates, _gas_template, _tpl_to_use, _component \
                    = _combine_stellar_templates(templates, gas_template, binned_tpl_wgts,
                                                 inp_component)

        # - Get the index of the nearest bin for every spaxel
        nearest_bin = np.argmin((x[:, None] - x_binned)**2 + (y[:, None] - y_binned)**2, axis=1)

        # - Use the binned data as the starting guess for the nearest
        #   spaxel
        stellar_wgts = stellar_wgts[nearest_bin,:]
        _templates = np.append(_templates[:nbin,:][nearest_bin,:], _templates[nbin:,:], axis=0)
        _gas_template = np.zeros(_templates.shape[0], dtype=bool)
        _gas_template[nspec:] = True

        _tpl_to_use = np.zeros((nspec,_templates.shape[0]), dtype=bool)
        _tpl_to_use[:,_gas_template] = True
        _tpl_to_use[np.arange(nspec),np.arange(nspec)] = True

        print('Templates to use for each spectrum')
        print(np.sum(_tpl_to_use, axis=1))
        print(_tpl_to_use.shape)
        print(_tpl_to_use[0,:])

        _component = np.append(np.zeros(nspec, dtype=int), _component[nbin:])
        print(_component)
        print(len(_component))

        # - Use the starting positions from the fit to the bins for the
        #   fit to the individual spaxels
        n_gas_comp = len(np.unique(_component[_gas_template]))
        _start = inp_start[nearest_bin]
        gas_start = binned_kin[:,np.absolute(moments[0]):][nearest_bin,:]
        for i in range(nspec):
            _start[i,1:] = np.array([ [gas_start[i].tolist()]*n_gas_comp ])
    else:
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

    _, _, _, tpl_wgts, _, _, _, _, _, kin, _ \
            = _fit_iteration(_templates, wave, flux, noise, velscale, start, moments, component,
                             _gas_template, tpl_to_use=_tpl_to_use, reject_boxcar=reject_boxcar,
                             velscale_ratio=velscale_ratio, degree=degree, mdegree=mdegree,
                             reddening=reddening, mask=mask, vsyst=vsyst, plot=plot, quiet=quiet)

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
                                                  mask=mask, vsyst=vsyst, plot=plot, quiet=quiet)

    # - Use the single output weight to renormalize the individual
    #   stellar template weights (only one of the weights for the
    #   non-gas templates should be non-zero); stellar-weight errors are
    #   always returned as 0.
    _tpl_wgts = stellar_wgts * np.sum(tpl_wgts[:,np.invert(_gas_template)], axis=1)[:,None]
    _tpl_wgts[:,gas_template] = tpl_wgts[:,_gas_template]
    _tpl_wgts_err = np.zeros((nspec,ntpl), dtype=float)
    _tpl_wgts_err[:,gas_template] = tpl_wgts_err[:,_gas_template]

    return model_flux, model_eml_flux, model_mask, _tpl_wgts, _tpl_wgts_err, addcoef, multcoef, \
                    ebv, kininp, kin, kin_err

    


# Function: emline_fitter_with_ppxf
# Input is:
#   - wave0: wavelength vector (in vacuum); shape is (Nwave,)
#   - flux: observed, unbinned flux; masked array with shape
#     (Nspaxel,Nwave)
#   - noise: error in observed, unbinned flux; masked array with
#     shape (Nspaxel,Nwave)
#   - sres: spectral resolution (R=lambda/delta lambda) as a
#     function of wavelength for each unbinned spectrum; shape
#     is (Nspaxel,Nwave)
#   - flux_binned: binned flux; masked array with shape
#     (Nbin,Nwave)
#   - noise_binned: noise in binned flux; masked array with
#     shape (Nbin,Nwave)
#   - sres_binned: spectral resolution (R=lambda/delta lambda)
#     as a function of wavelength for each binned spectrum;
#     shape is (Nbin,Nwave)
#   - velscale: Velocity step per pixel
#   - velscale_ratio: Ratio of velocity step per pixel in the
#     observed data versus in the template data
#   - dv: Velocity offset between the galaxy and template data
#     due to the difference in the initial wavelength of the
#     spectra
#   - stars_vel: Velocity of the stellar component; shape is
#     (Nbins,)
#   - stars_sig: Velocity dispersion of the stellar component;
#     shape is (Nbins,)
#   - templates: Template flux; shape is (Ntpl,Ntplwave)
#   - gas_tpl: Boolean array for extracting gas templates; shape
#     is (Ntpl,)
#   - guess_vel: Initial guess velocity for the gas components;
#     shape is (Nbins,)
#   - guess_sig: Initial guess velocity dispersion for the gas
#     components; shape is (Nbins,)
#   - gas_names: Name of the gas templats; shape is (Ngastpl,)
#   - eml_wave: Rest wavelength of gas templates (in vacuum);
#     shape is (Ngastpl,)
#   - eml_component: Component vector for ppxf (each kinimatic
#     group has a unique index, and an index below zero will
#     automatically fix the moments); shape is (Nkin_group,)
#   - eml_tied: Tied vector for ppxf (to decide the tying 
#     relations among gas components); shape is (Nkin_group,
#     Nmoment)
#   - template_sres: spectral resolution (R=lambda/delta
#     lambda) as a function of wavelength for all the templates
#     templates; shape is (Ntplwave,)
#   - degree: Additive polynomial order
#   - mdegree: Multiplicative polynomial order
#   - x: On-sky spaxel x coordinates; shape is (Nspaxel,)
#   - y: On-sky spaxel y coordinates; shape is (Nspaxel,)
#   - x_binned: On-sky bin x coordinate; shape is (Nbin,)
#   - y_binned: On-sky bin y coordinate; shape is (Nbin,)
#   - mask_binned: Masks for binned spectra; shape is (Nbin,
#     Nwave)
#   - mask_drp: Masks for unbinned spectra; shape is (Nbin,
#     Nwave)
#   - nsa_z: NSA (or guess) redshift of the galaxy
#
# Output is:
#   - model_flux: stellar-continuum + emission-line model; shape
#     is (Nmod, Nwave); first axis is ordered by model ID number
#   - model_eml_flux: model emission-line flux only; shape is
#     (Nmod, Nwave); first axis is ordered by model ID number
#   - model_mask: boolean or bit mask for fitted models; shape
#     is (Nmod, Nwave); first axis is ordered by model ID number
#   - model_binid: ID numbers assigned to each spaxel with a
#     fitted model; any spaxel without a model should have
#     model_binid = -1; the number of >-1 IDs must be Nmod;
#     shape is (Nspaxel,)
#   - weights: weights of the templates (the best fit template 
#     can be obtained by setting: bestemp = templates @ weights);
#     shape is (Nmod, Ntpl)
#   - weights_err: 1 sigma formal error of weights; shape is
#     (Nmod, Ntpl)
#   - mweights: coefficients of multiplicative Legendre polynomials
#     of order > 0; shape is (Nmod, mdegree)
#   - eml_kin: Kinematics (velocity and velocity dispersion) of
#     each emission line; shape is (Nmod,Neml,Nkin)
#   - eml_kinerr: Error in the kinematics of each emission line
#   - eml_sigmacorr: Quadrature corrections required to obtain
#     the astrophysical velocity dispersion; shape is
#     (Nmod,Neml); corrections are expected to be applied as
#     follows:
#       sigma = numpy.ma.sqrt( numpy.square(eml_kin[:,:,1])
#                               - eml_sigmacorr)
#     NOTE THAT: eml_sigmacorr here are actually in unit of 
#     (km/s)^2 and could have negative values
# --------------------------------------------------------------

def emline_fitter_with_ppxf(wave0, flux, noise, sres, flux_binned, noise_bin,
                            velscale, velscale_ratio, dv, vel, sig, templates, gas_tpl, 
                            vel_gas, sig_gas, gas_names, eml_wave, eml_component, 
                            eml_tied, templates_sres, degree, mdegree, x, y, xbin, 
                            ybin, mask_binned, mask_drp, nsa_z, debug=False,
                            mode='bin_to_spaxel'):

    if debug:
        warnings.warn('JUST DEBUGGING.  NO EMISSION-LINE FITS PERFORMED!!')
        n_spaxels, nwave = flux.shape
        n_eml = np.sum(gas_tpl)
        masked_fraction = np.sum(mask_drp, axis=1)/nwave
        indx = masked_fraction > 1/3
        model_binid = np.full(n_spaxels, -1, dtype=int)
        model_binid[indx] = np.arange(np.sum(indx))

        model_flux = np.zeros(flux.shape, dtype=float)
        model_eml_flux = np.zeros(flux.shape, dtype=float)
        model_mask = np.zeros(flux.shape, dtype=float)
        
        eml_flux = np.zeros([n_spaxels,n_eml], dtype=float)
        eml_fluxerr = np.zeros([n_spaxels,n_eml], dtype=float)
        eml_kin = np.zeros([n_spaxels,n_eml,2], dtype=float)
        eml_kinerr = np.zeros([n_spaxels,n_eml,2], dtype=float)
    
        eml_sigmacorr = np.zeros([n_spaxels,n_eml], dtype=float)
        return model_flux, model_eml_flux, model_mask, model_binid, eml_flux, eml_fluxerr, \
                    eml_kin, eml_kinerr, eml_sigmacorr

    # templates should have the shape (Nwave, Ntpl) for ppxf
    templates = templates.T
    
    # Total number of bins
    n_binnum = flux_binned.shape[0]
    
    # Total number of spaxels
    n_spaxels = flux.shape[0]
    
    # Total number of emission lines
    n_eml = np.sum(gas_tpl)
    
    # Total number of stellar templates
    n_temps = np.sum(~gas_tpl)

    gas_vel = np.zeros(n_binnum)
    gas_sig = np.zeros(n_binnum)
    gas_flux_binned = np.zeros([n_binnum,n_eml])
    star_weights = np.zeros([n_binnum, n_temps])
    star_weights_err = np.zeros([n_binnum, n_temps])
    bestfit_template = np.zeros([n_binnum, templates.shape[0]])
    
    # Testing use
    #beg = n_binnum if os.path.isfile(npz_file) else 0
    
    # 'bin' mode contains only one iteration, fitting the binned data
    # 'bin_to_spaxel' mode contains two iterations, fitting both the binned and
    # unbinned data
    if isinstance(mode, str) and mode == 'bin':
        beg1 = 0
        beg2 = n_spaxels
        n = n_binnum
        model_shape = flux_binned.shape
        model_binid = None
    elif isinstance(mode, str) and mode == 'bin_to_spaxel':
        beg1 = 0
        beg2 = 0
        n = n_spaxels
        model_shape = flux.shape
        model_binid = np.zeros(n_spaxels)
        id_num = 0
        # Get the index of the nearest bin for every spaxel
        nearest_bin = np.argmin((x[:, None] - xbin)**2 + (y[:, None] - ybin)**2, axis=1)
    
    # Initialize output arrays
    eml_flux = np.zeros([n,n_eml])
    eml_fluxerr = np.zeros([n,n_eml])
    eml_kin = np.zeros([n,n_eml,2])
    eml_kinerr = np.zeros([n,n_eml,2])
    weights = np.zeros([n,templates.shape[1]])
    weights_err = np.zeros([n,templates.shape[1]])
    mweights = np.zeros([n,mdegree])

    eml_sigmacorr = np.zeros([n,n_eml])

    model_flux = np.zeros(model_shape)
    model_eml_flux = np.zeros(model_shape)
    model_mask = np.zeros(model_shape)
    
    # First iteration
    for i in range(beg1, n_binnum):
        print('Fitting bin: {0}/{1}'.format(i+1,n_binnum))
        
        # Normalize spectrum to avoid numerical issues
        galaxy = flux_binned[i,:]/np.median(flux_binned[i,:])
        
        noise_binned = noise_bin[i,:]/np.median(flux_binned[i,:])
        
        # Use the first moment of Halpha line as the initial velocity
        # guess, and use the stellar velocity dispersion as initial
        # sigma guess for gas components
        #kin = [vel_gas[i], 100.]
        kin = [vel_gas[i], sig[i]]
        
        # Assign all emission lines with single component to tie their
        # kinematics (and the same for stellar templates)
        component = np.asarray([0]*n_temps + [1]*n_eml) if mode == 'bin_to_spaxel'\
        else np.append(np.zeros(n_temps,dtype=int), eml_component+1)
        
        # The negative moments automatically fix the stellar kinematics
        moments = [-2, 2] if mode == 'bin_to_spaxel'\
        else [-2] + [2]*np.max(component+1)
        
        # Tie emission lines if one selects the 'bin' mode
        tied = None if mode == 'bin_to_spaxel' else eml_tied
        
        # Initial guess for components (for stellar components, these
        # are fixed values)
        start = [[vel[i], sig[i]]]
        start += [kin for j in range(len(moments)-1)]
        
        # Tie kinematics of all emission lines while fixing the stellar kinematics in the first
        # fit
        pp = ppxf(templates, galaxy, noise_binned, velscale, start, velscale_ratio=velscale_ratio,
                  plot=False, moments=moments, degree=degree, mdegree=mdegree, tied=tied,
                  mask=mask_binned[i,:], vsyst=dv, lam=wave0, component=component, quiet=True,
                  gas_component=component>0, gas_names=gas_names)
        
        # Calculate the residuals
        resid = galaxy - pp.bestfit
        
        # Calculate the noises
        NOISE = calculate_noise(resid, width=reject_boxcar)
        
        # Mask out pixels with residuals > 3 sigmas
        mask = (abs(resid) < (3*NOISE)) & mask_binned[i,:]
        
        # Add a three-sigma cut in the second fit while tying the kinematics as the first one
        pp = ppxf(templates, galaxy, noise_binned, velscale, start, velscale_ratio=velscale_ratio,
                  mask=mask, plot=False, moments=moments, degree=degree, mdegree=mdegree,
                  vsyst=dv, lam=wave0, component=component, quiet=True, tied=tied,
                  gas_component=component>0, gas_names=gas_names)
        
        gas_vel[i] = pp.sol[1][0]
        gas_sig[i] = pp.sol[1][1]
        # Scale back the fluxes; note that the fluxes here are caculated by direct summing rather
        # than integration! If one uses ppxf_util to construct input templates, the ouput fluxes
        # should be corrected using pixel step near the lines
        gas_flux_binned[i] = pp.gas_flux*np.median(flux_binned[i,:])
        star_weights[i] = pp.weights[0:n_temps]
        bestfit_template[i,:] = templates[:,~gas_tpl] @ star_weights[i]
        
        design_matrix = pp.matrix/pp.noise[:, None]
        star_weights_err[i] = capfit.cov_err(design_matrix)[1][0:n_temps]
        
        # The ouput of the fitting are parsed in the following lines
        # (for 'bin' mode)
        if mode == 'bin':
            eml_kin[i][:,0] = np.full(eml_kin.shape[1], pp.sol[1][0])
            eml_kin[i][:,1] = pp.sol[:,1][eml_component+1]

            # The errors are calculated according to docstring in ppxf.py
            eml_kinerr[i][:,0] = np.full(eml_kin.shape[1], pp.error[1][0])*np.sqrt(pp.chi2)
            eml_kinerr[i][:,1] = pp.error[:,1][eml_component+1]*np.sqrt(pp.chi2)

            # Full models of flux and models of emission lines
            model_flux[i,:] = pp.bestfit*np.median(flux_binned[i,:])

            spectra = pp.matrix[:,degree+1:]
            gas_spectrum = spectra[:,pp.gas_component].dot(pp.weights[pp.gas_component])
            model_eml_flux[i,:] = gas_spectrum*np.median(flux_binned[i,:])

            # Masks for models
            model_mask[i,:] = mask

            # Calculate the sigma corrections
            FWHM = wave0/np.max(sres, axis=0)
            eml_wave_new = eml_wave*np.exp(pp.sol[1][0]/astropy.constants.c.to('km/s').value)
            index_new = np.argmin(abs(wave0-eml_wave_new[:,None]),axis=1)
            FWHM_diff2 = FWHM[index_new]**2 - FWHM[index]**2
            eml_sigmacorr[i] = FWHM_diff2/((2.355*step)**2)*(velscale**2)

            # Weights and weight errors of gas tempaltes
            # Note that the star weights errors should never be used
            weights[i] = np.append(star_weights[i], pp.weights[1:])
            weights_err[i] = np.append(star_weights_err[i], pp.gas_flux_error)

            # Coefficients of multiplicative Legendre polynomials of order > 0
            mweights[i] = pp.mpolyweights
        
        #plt.clf()
        #pp.plot()
        #plt.show()
    
#    # Save the results (only for testing)
#    np.savez_compressed('/Users/mac/Documents/New_test_results/first_'+plate+'_'+ifu+'.npz', 
#                        gas_vel = gas_vel, gas_sig = gas_sig,
#                        gas_flux_binned = gas_flux_binned,
#                        bestfit_template = bestfit_template)

    # Second iteration
    # Load data from the first iteration (only for testing)
#    first_iter = np.load('/Users/mac/Documents/New_test_results/first_'+plate+'_'+ifu+'.npz')
#    gas_vel = first_iter['gas_vel']
#    gas_sig = first_iter['gas_sig']
#    bestfit_template = first_iter['bestfit_template']
    
     # Testing use
#    if beg == 0:
#        print('writing {0}'.format(npz_file))
#        np.savez_compressed(npz_file, gas_vel = gas_vel, gas_sig = gas_sig,
#                            gas_flux_binned = gas_flux_binned, bestfit_template = bestfit_template)
#    else:
#        print('reading {0}'.format(npz_file))
#        inp = np.load(npz_file)
#        gas_vel = inp['gas_vel']
#        gas_sig = inp['gas_sig']
#        bestfit_template = inp['bestfit_template']
    
#    for i in range(1084,n_spaxels):
    for i in range(beg2, n_spaxels):
        print('Fitting spaxel: {0}/{1}'.format(i+1,n_spaxels))
        
        # Here start the first half of the second iteration
        
        mask_DRP = mask_drp[i,:]
        
        # Continue the loop if the current spaxel has no spectrum
        # TODO: decide if 1/3 is appropriate to indicate insufficient
        # pixels?
        if sum(~mask_DRP) > flux[i,:].shape[0]/3:
            model_binid[i] = -1
            continue
        else:
            model_binid[i] = id_num
            id_num += 1
        
        # Normalize spectrum to avoid numerical issues
        galaxy = flux[i,:]/np.median(flux[i,:])
                
        noise_spaxel = noise[i,:]/np.median(flux[i,:])
        noise_spaxel[(noise_spaxel<=0.)|(~np.isfinite(noise_spaxel))] = 1.
        
        # This is the index of the nearest bin to the spaxel being used
        k = nearest_bin[i]
        
        gas_templates = templates[gas_tpl]
        templates_sec = np.column_stack([bestfit_template[k], gas_templates])
        
        # Only one stellar template is used in the second iteration
        n_temps = 1
        
        # This is the gas kinematics [V, sigma] of the nearest bin which
        # will be used as the initial guess for emission-line fitting
        kin = [gas_vel[k], gas_sig[k]]
        
        # Combine the input gas components and stellar ones
        component = np.asarray([0]*n_temps + [1]*n_eml)
        
        moments = [-2, 2]
        
        # Adopt the same value in the first iteration for star component
        start = [[vel[k], sig[k]]]
        # Adopt the same starting value for all gas component
        start += [kin for j in range(len(moments)-1)]
        
        # In the third fit, the tying of kinematics is the same as the first two while the initial
        # guess for gas components, input noise and mask are different
        pp = ppxf(templates_sec, galaxy, noise_spaxel, velscale, start, velscale_ratio=velscale_ratio,
                  mask=mask_DRP, plot=False, moments=moments, degree=degree, mdegree=mdegree,
                  vsyst=dv, lam=wave0, component=component, quiet=True,
                  gas_component=component>0, gas_names=gas_names)
        
        # Calculate the residuals
        resid = galaxy - pp.bestfit
        
        # Calculate the noises
        NOISE = calculate_noise(resid, width=reject_boxcar)
        
        # Mask out pixels with residuals > 3 sigmas
        mask = (abs(resid) < (3*NOISE))&mask_DRP
        #print('Outliers: '+str(3109 - sum(mask)))
        
        # The forth fit adds a three-sigma cut after the third one
        pp = ppxf(templates_sec, galaxy, noise_spaxel, velscale, start, velscale_ratio=velscale_ratio,
                  mask=mask, plot=False, moments=moments, degree=degree, mdegree=mdegree,
                  vsyst=dv, lam=wave0, component=component, quiet=True,
                  gas_component=component>0, gas_names=gas_names)
        
        # Here start the second half of the second iteration
        
        # This is the gas kinematics [V, sigma] of the current spaxel from the last 
        # call of ppxf which will be used as the initial guess for emission-line fitting
        kin = [pp.sol[1][0], pp.sol[1][1]]
        
        # Combine the input gas components and stellar ones
        # Now each kinematic group has a unique component index (>0)
        component = np.append(np.zeros(n_temps,dtype=int), eml_component+1)
        
        moments = [-2] + [2]*np.max(component+1)
        
        # Adopt the same value in the first iteration for star component
        start = [[vel[k], sig[k]]]
        # Adopt the same starting value for all gas component
        start += [kin for j in range(len(moments)-1)]
        
        # In the fifth fit, the velocities and sigmas of emission lines are tied as desired
        pp = ppxf(templates_sec, galaxy, noise_spaxel, velscale, start, velscale_ratio=velscale_ratio,
                  mask=mask_DRP, quiet=True, plot=False, moments=moments, degree=degree,
                  mdegree=mdegree, vsyst=dv, lam=wave0, component=component, tied=eml_tied,
                  gas_component=component>0, gas_names=gas_names)
        
        # The ouput of the fitting are parsed in the following lines
        
        eml_kin[i][:,0] = np.full(eml_kin.shape[1], pp.sol[1][0])
        eml_kin[i][:,1] = pp.sol[:,1][eml_component+1]
        
        # The errors are calculated according to docstring in ppxf.py
        eml_kinerr[i][:,0] = np.full(eml_kin.shape[1], pp.error[1][0])*np.sqrt(pp.chi2)
        eml_kinerr[i][:,1] = pp.error[:,1][eml_component+1]*np.sqrt(pp.chi2)
        
        # Integral fluxes and their errors of all emission lines (not in the final output)
        # Reset pixel-step for each emline (since ppxf returns flux caculated by direct summing)
        eml_wave0 = eml_wave*(1+nsa_z)
        index = np.argmin(abs(wave0-eml_wave0[:,None]),axis=1)
        step = wave0[index+1] - wave0[index]
        # Remain to be examined if /(1+nsa_z) correction is needed...
        eml_flux[i] = pp.gas_flux*np.median(flux[i,:])*step
        eml_fluxerr[i] = pp.gas_flux_error*np.sqrt(pp.chi2)*np.median(flux[i,:])*step
        
        # Full models of flux and models of emission lines
        model_flux[i,:] = pp.bestfit*np.median(flux[i,:])
        
        spectra = pp.matrix[:,degree+1:]
        gas_spectrum = spectra[:,pp.gas_component].dot(pp.weights[pp.gas_component])
        model_eml_flux[i,:] = gas_spectrum*np.median(flux[i,:])
        
        # Masks for models
        model_mask[i,:] = mask
        
        # Calculate the sigma corrections
        FWHM = wave0/np.max(sres, axis=0)
        eml_wave_new = eml_wave*np.exp(pp.sol[1][0]/astropy.constants.c.to('km/s').value)
        index_new = np.argmin(abs(wave0-eml_wave_new[:,None]),axis=1)
        FWHM_diff2 = FWHM[index_new]**2 - FWHM[index]**2
        eml_sigmacorr[i] = FWHM_diff2/((2.355*step)**2)*(velscale**2)
        
        # Weights and weight errors of gas tempaltes
        # Note that the star weights errors should never be used
        weights[i] = np.append(star_weights[k]*pp.weights[0], pp.weights[1:])
        weights_err[i] = np.append(star_weights_err[k]*pp.weights[0], pp.gas_flux_error)
        
        # Coefficients of multiplicative Legendre polynomials of order > 0
        mweights[i] = pp.mpolyweights
        
        #plt.clf()
        #pp.plot()
        #plt.show()
    
#    # Save the results (only for testing)
#    np.savez_compressed('/Users/mac/Documents/New_test_results/second_new'+plate+'_'+ifu+'.npz', 
#                        model_flux = model_flux, model_eml_flux = model_eml_flux, model_mask = model_mask,
#                        model_binid = model_binid, eml_flux = eml_flux, eml_fluxerr = eml_fluxerr,
#                        eml_kin = eml_kin, eml_kinerr = eml_kinerr, eml_sigmacorr = eml_sigmacorr)
#    
#    second_iter = np.load('/Users/mac/Documents/New_test_results/second_new'+plate+'_'+ifu+'.npz')
#    model_flux = second_iter['model_flux']
#    model_eml_flux = second_iter['model_eml_flux']
#    model_mask = second_iter['model_mask']
#    model_binid = second_iter['model_binid']
#    eml_flux = second_iter['eml_flux']
#    eml_fluxerr = second_iter['eml_fluxerr']
#    eml_kin = second_iter['eml_kin']
#    eml_kinerr = second_iter['eml_kinerr']
#    eml_sigmacorr = second_iter['eml_sigmacorr']
    
    return model_flux, model_eml_flux, model_mask, model_binid, weights, weights_err, mweights,\
           eml_kin, eml_kinerr, eml_sigmacorr


