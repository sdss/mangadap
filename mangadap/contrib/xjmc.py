# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Implements an emission-line fitting function using pPXF.

----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

import time
import warnings

from IPython import embed

import numpy as np

from scipy import fftpack
from scipy.ndimage import rank_filter

import astropy.constants

from ppxf import ppxf, capfit

# For debugging
from matplotlib import pyplot as plt
# from mangadap.proc.ppxffit import PPXFModel


# NOTE: Pulled from mangadap.proc.ppxffit.PPXFFit.ppxf_tpl_obj_voff
def ppxf_vsyst(tpl_wave, obj_wave, velscale, velscale_ratio=None):
    """
    Determine the pseudo velocity offset between the template and
    object spectra due to the difference in the starting wavelengths.

    Calculation is independent of the base of the logarithm used the
    sampling of the spectra, but the wavelengths *must* be
    geometrically binned.
    
    .. todo::

        Compute ``velscale_ratio`` directly from the input?

    Args:
        tpl_wave (`numpy.ndarray`_):
            Wavelength vector for the template library to fit to the
            object spectrum.
        obj_wave (`numpy.ndarray`_):
            Wavelength vector for the object spectrum to be fit.
        velscale (:obj:`float`):
            Velocity step per pixel in km/s for the
            **object** spectrum.
        velscale_ratio (:obj:`int`, optional): 
            The **integer** ratio between the velocity scale of the
            pixel in the galaxy data to that of the template data.
            This is used only when constructing the template library.
            If None, assumes that the velocity scales are identical;
            i.e., ``velscale_ratio=None`` is identical to
            ``velscale_ratio=1``.
    
    Returns:
        :obj:`float`: Velocity offset in km/s between the initial
        wavelengths of the template and object spectra.
    """
    dlogl = np.log(obj_wave[0])-np.log(tpl_wave[0]) if velscale_ratio is None \
                else np.log(obj_wave[0])-np.mean(np.log(tpl_wave[0:velscale_ratio]))
    return dlogl*velscale / np.diff(np.log(obj_wave[0:2]))[0]


def ppxf_tied_parameters(component, vgrp, sgrp, moments):
    r"""
    Construct the object used to tie kinematic parameters in pPXF.

    .. note::

        The components and kinematics groups can potentially have
        redundant information. E.g., If all sigma groups also have
        tied velocities, they'll be part of a component and will not
        required a tied parameter.

    .. warning::

        High-order moments are not currently tied, but this is
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
        :obj:`list`: The tied object to pass to ``ppxf``. Will be
        None if both vgrp or sgrp are None or when the consolidation
        of the groups and components result in no tying being
        necessary.

    Raises:
        ValueError:
            Raised if the components or the groups are not an
            uninterupted sequence from 0..N-1, if the moments are
            provided for each component, or the length of the
            ``vgrp`` or ``sgrp`` arrays do not match the input
            component array.
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
        residuals (`numpy.ndarray`_):
            Vector of residuals between the model and data
        width (:obj:`int`, optional):
            Size of the window for the statistics. Should be an odd
            number and larger and 10.
    
    Returns:
        `numpy.ndarray`_: Vector with the same size as the input
        residuals with half of the 68% confidence interval of the
        residuals within a box of size `width` centered at each
        position.

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
        component (`numpy.ndarray`_):
            The full list of components for all templates.  The total
            number of components is NCOMP.
        gas_template (`numpy.ndarray`_):
            A boolean array identifying the gas templates.
        start (:obj:`list`, `numpy.ndarray`_):
            The starting kinematics to use for each object spectrum.
            Shape is (NOBJ,); each element has shape (NCOMP,); each
            component element has shape (NMOM,), such that the number of
            moments can be different for each component.  NCOMP and NMOM
            should be the same for all object spectra.
        single_gas_component (:obj:`bool`, optional):
            Flag to force all gas components to have the same
            kinematics.
        gas_start (:obj:`list`, `numpy.ndarray`_, optional):
            The starting kinematics to use for the gas components.
            Shape must be (NOBJ,2).

    Returns:
        :obj:`tuple`: Three numpy arrays are returned:

            #. The list of downselected and renumbered, if necessary,
               kinematic components,
            #. the downselected and reordered, if necessary, number
               of moments to fit, and
            #. the downselected and reordered starting guesses for
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


def _reset_components(c, valid):
    r"""
    Down-select and appropriately reindex a set of ppxf kinematic components.
    
    Args:
        c (`numpy.ndarray`_):
            Integer array with the kinematic components assigned to
            each template. Shape is :math:`(N_{\rm tpl},)`.
        valid (`numpy.ndarray`_):
            Boolean array selecting the valid templates for component
            reassignment.

    Returns:
        :obj:`tuple`: Returns the down-selected and re-ordered set of
        components and an integer array with the mapping between the
        old and new components; i.e., see ``component_map`` in
        :func:`_reorder_solution`.
    """
    # Select valid components
    _c = c[valid]
    # Map between old and new component numbers
    c_map, inv = np.unique(_c, return_inverse=True)
    # Return the new component numbers and the old-to-new map
    return np.arange(len(c_map))[inv], c_map


def _reset_kinematic_constraints(comp, valid, A, b):
    r"""
    Down-select the kinematic constraints set for each component base on the
    boolean flagging of valid templates.

    Args:
        comp (`numpy.ndarray`_):
            Integer array with the kinematic components assigned to
            each template. Shape is :math:`(N_{\rm tpl},)`.
        valid (`numpy.ndarray`_):
            Boolean array selecting the valid templates for component
            reassignment.
        A (`numpy.ndarray`_):
            Linear constraint matrix used by ppxf to impose constraints on
            the fitted kinematics; see the ppxf keyword argument
            ``constr_kinematics``. Shape is :math:`(N_{\rm constr},2 N_{\rm
            comp})`; i.e., each row of A constitutes an imposed constraint
            and each column pair provides the coefficient for each kinematic
            moment. This method currently only allows for 2 moments in the
            fit.  Can be None.
        b (`numpy.ndarray`_):
            Linear constraint vector used by ppxf to impose constraints one
            the fitted kinematics; see the ppxf keyword argument
            ``constr_kinematics``. Shape is :math:`(N_{\rm constr},)`. Can be
            None.

    Returns:
        :obj:`tuple`: Return the adjusted ``A`` and ``b`` objects after
        accounting for any components removed because their constituent
        templates were all invalid. If ``A`` **or** ``b`` are None, both
        objects are returned as None; otherwise, both are adjusted arrays.
        Note that if the arrays do not require adjustment (i.e., no
        components have been removed), the input ``A`` and ``b`` are
        returned, *not* a copy of these arrays.
    """
    if A is None or b is None:
        # Nothing to do
        return None, None

    # Find the components removed
    removed_components = list(set(np.unique(comp)) - set(np.unique(comp[valid])))

    if len(removed_components) == 0:
        # No components removed so just return the input
        return A, b

    # Number of components
    ncomp = A.shape[1]//2

    # Find the good columns in A ...
    good_columns = np.logical_not(np.isin(np.repeat(np.arange(ncomp), 2), removed_components))
    # ... and the good rows
    good_rows = np.all(A[:,np.logical_not(good_columns)] == 0., axis=1)

    if not np.any(good_rows) or not np.any(good_columns):
        # No valid constraints remain
        return None, None

    # Return the remaining constraints
    return A[np.ix_(good_rows, good_columns)], b[good_rows]


def _good_templates(templates, gas_template, mask, start, velscale, velscale_ratio, vsyst):
    """
    Determine which of the templates to use during the fit.

    Good template are *any* stellar-continuum template and any gas
    template that is non-zero over the expected fitting range.

    The details of this code should be pulled directly from ppxf for
    consistency.

    Args:
        templates (`numpy.ndarray`_):
            Templates library to use for fitting.  Shape is
            (NTPLPIX,NTPL).
        gas_template (`numpy.ndarray`_):
            Boolean vector that selects the gas templates.  Shape is
            (NTPL,).
        mask (`numpy.ndarray`_):
            Boolean vector that selects the pixels in the object
            spectrum to fit (i.e., mask=True for pixels to fit and
            mask=False for pixels to ignore).  As in pPXF, the length is
            expected to be less than or equal to NTPLPIX in the template
            spectra.
        start (list):
            The starting kinematics for each kinematic component.
            Length is NCOMP.
        velscale (float):
            The pixel scale of the object spectrum to fit in km/s.
        velscale_ratio (:obj:`int`, optional):
            The ratio between the object and template pixel scale.
        vsyst (:obj:`float`, optional): 
            The pseudo velocity shift between the template and object
            spectra just due to the difference in the starting
            wavelength.

    Returns:
        `numpy.ndarray`_: Boolean array flagging good templates,
        which means anything that is not a gas template and any gas
        template that has non-zero values in the fitting region.
    """
    valid = np.ones(templates.shape[0], dtype=bool)
    if not np.any(gas_template):
        # No gas templates, which should never happen in this module!
        return valid

    # Starts at line 1545 of ppxf version 7.4.0
    vmed = np.median([a[0] for a in start])
    dx = int(np.round((vsyst + vmed)/velscale))  # Approximate velocity shift
    gtpl = templates[gas_template,:].T
    tmp = ppxf.rebin(gtpl, velscale_ratio)
    gas_peak = np.max(np.abs(tmp), axis=0)
    tmp = np.roll(tmp, dx, axis=0)
    good_peak = np.max(np.abs(tmp[np.where(mask)[0],:]), axis=0)
    valid[gas_template] = np.logical_not(good_peak <= gas_peak/1e3)
    return valid


def _validate_templates_components(templates, gas_template, component, vgrp, sgrp, moments, mask,
                                   start, tpl_to_use, velscale, velscale_ratio=None, vsyst=0,
                                   constr_kinem=None):
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
        templates (`numpy.ndarray`_):
            Templates library to use for fitting.  Shape is
            (NTPLPIX,NTPL).
        gas_template (`numpy.ndarray`_):
            Boolean vector that selects the gas templates.  Shape is
            (NTPL,).
        component (`numpy.ndarray`_):
            Integer vector identifying the kinematic component for each
            template.  Shape is (NTPL,).
        vgrp (array-like):
            The integer velocity group associated with each template.
            Shape is :math:`(N_{\rm tpl},)` , but can be None.
        sgrp (array-like):
            The integer sigma group associated with each template.
            Shape is :math:`(N_{\rm tpl},)` , but can be None. 
        moments (`numpy.ndarray`_):
            Integer vector with the number of kinematic moments for each
            component.  Shape is (NCOMP,).
        mask (`numpy.ndarray`_):
            Boolean vector that selects the pixels in the object
            spectrum to fit (i.e., mask=True for pixels to fit and
            mask=False for pixels to ignore).  As in pPXF, the length is
            expected to be less than or equal to NTPLPIX in the template
            spectra.
        start (list):
            The starting kinematics for each kinematic component.
            Length is NCOMP.
        tpl_to_use (`numpy.ndarray`_):
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
        constr_kinem (:obj:`dict`, optional):
            A dictionary with the constraints to apply to the kinematics of
            each component; see ``constr_kinem`` in ppxf.

    Returns:

        :obj:`tuple`: Eleven objects are returned. The number of new valid
        templates is listed as _NTPL:

            #. Boolean vector with which templates were valid (length is
               NTPL),

            #. the valid set of templates to pass to pPXF [shape is
               (NTPLPIX,_NTPL)],

            #. the boolean vector selecting gas templates (length is _NTPL),

            #. an integer vector mapping the input component number to the
               output component number (i.e., component_map[0] is the
               original component number for the downselected 0th component;
               length is _NCOMP),

            #. the new component number for the downselected templates (length
               is _NTPL),

            #. the new velocity group number for the downselected templates
               (length is _NTPL),

            #. the new sigma group number for the downselected templates
               (length is _NTPL),

            #. the number of moments for the new components (length
               is _NCOMP),

            #. the starting kinematics for each new component (length
               is _NCOMP),

            #. the parameter tying object reordered as necessary for the new
               component list, and

            #. the kinematics constraint dictionary for the down-selected set
               of templates.


    Raises:
        ValueError:
            Raised if no templates are valid.

    """
    # Find the valid templates
    valid = _good_templates(templates, gas_template, mask, start, velscale, velscale_ratio, vsyst)
    valid &= tpl_to_use

    # None (!) of the templates are valid
    if not np.any(valid):
        raise ValueError('No valid templates in fit!')

    # Number of kinematic components
    ncomp = np.max(component)+1

    # All templates are valid, just return the input
    if np.all(valid):
        return valid, templates, gas_template, np.arange(ncomp), component, \
                    vgrp, sgrp, moments, start, constr_kinem

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

    _constr_kinem = None
    if constr_kinem is not None:
        _constr_kinem = {}
        if 'A_ineq' in constr_kinem.keys():
            A_ineq, b_ineq = _reset_kinematic_constraints(component, valid, constr_kinem['A_ineq'],
                                                          constr_kinem['b_ineq'])
            if A_ineq is not None:
                _constr_kinem['A_ineq'], _constr_kinem['b_ineq'] = A_ineq, b_ineq
        if 'A_eq' in constr_kinem.keys():
            A_eq, b_eq = _reset_kinematic_constraints(component, valid, constr_kinem['A_eq'],
                                                      constr_kinem['b_eq'])
            if A_eq is not None:
                _constr_kinem['A_eq'], _constr_kinem['b_eq'] = A_eq, b_eq
        if len(_constr_kinem.keys()) == 0:
            _constr_kinem = None

    return valid, _templates, _gas_template, component_map, _component, \
                _vgrp, _sgrp, _moments, _start, _constr_kinem


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
        component_map (`numpy.ndarray`_):
            An integer vector mapping the input component number to the
            output component number; i.e., component_map[0] is the
            original component number for the downselected 0th
            component; length is NCOMP.
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


# TODO: Change 'mask' to 'gpm'
def _fit_iteration(tpl_wave, templates, wave, flux, noise, velscale, start, moments, component,
                   gas_template, tpl_to_use=None, reject_boxcar=101, velscale_ratio=None,
                   degree=-1, mdegree=0, reddening=None, vgrp=None, sgrp=None, constr_kinem=None,
                   mask=None, plot=False, quiet=True, sigma_rej=3., starting_spectrum=None,
                   ppxf_faults='flag'):
    r"""
    Run a single fit+rejection iteration of pPXF for all input
    spectra with the provided set of constraints/options.

    Args: 
        tpl_wave (`numpy.ndarray`_):
            Wavelength vector for the template library to fit to the
            object spectrum.
        templates (`numpy.ndarray`_):
            Full template library.  Shape is (NTPL,NTPLPIX).
        wave (`numpy.ndarray`_):
            Wavelength vector.  Shape is (NPIX,).
        flux (`numpy.ndarray`_):
            Object spectra to fit.  Shape is (NSPEC,NPIX).
        noise (`numpy.ndarray`_):
            Error in the object spectra to fit.  Shape is (NSPEC,NPIX).
        velscale (float):
            The pixel scale of the object spectra in km/s.
        start (list):
            The starting kinematics for each kinematic component.
            Length is NCOMP.
        moments (`numpy.ndarray`_):
            Integer vector with the number of kinematic moments for each
            component.  Shape is (NCOMP,).
        component (`numpy.ndarray`_):
            Integer vector identifying the kinematic component for each
            template.  Shape is (NTPL,).
        gas_template (`numpy.ndarray`_):
            Boolean vector that selects the gas templates.  Shape is
            (NTPL,).
        tpl_to_use (`numpy.ndarray`_, optional):
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
        constr_kinem (:obj:`dict`, optional):
            A dictionary with the constraints to apply to the kinematics of
            each component; see ``constr_kinem`` in ppxf.
        mask (`numpy.ndarray`_, optional):
            Boolean vector that selects the pixels in the object spectra
            to fit (i.e., mask=True for pixels to fit and mask=False for
            pixels to ignore).  Shape is (NSPEC,NPIX).  All pixels are
            fit by default.
        vsyst (:obj:`float`, optional):
            The pseudo velocity shift between the template and object
            spectra just due to the difference in the starting
            wavelength.
        plot (:obj:`bool`, optional):
            Show the pPXF fit plot at each iteration.
        quiet (:obj:`bool`, optional):
            Suppress output to the terminal.
        sigma_rej (:obj:`float`, optional):
            Sigma values used for the rejection.
        starting_spectrum (:obj:`int`, optional):
            Select a spectrum index to start the iteration. Only used
            for debugging! If None, starts with the first spectrum.
        ppxf_faults (:obj:`str`, optional):
            Dictates how exceptions raised by ppxf are treated.
            Allowed values are:

                - ``'flag'``: Log the fault and continue. Any
                  spectrum that faults is flagged as such in the last
                  object returned by the method (see below).
                - ``'raise'``: Re-raise the exception returned by
                  ppxf.
            
    Returns:
        :obj:`tuple`: Twelve arrays are returned:

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
            - (12) a boolean array flagging spectra that caused ppxf
              to fault; True means the fit faulted, False means it
              was successful.

    """

    if ppxf_faults not in ['flag', 'raise']:
        raise ValueError('Keyword ppxf_faults must be \'flag\' or \'raise\'.')

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
    fault = np.zeros(nspec, dtype=bool)

    # Use the mask to set the starting and ending pixels and the psuedo
    # velocity offset
    ps = np.zeros(nspec, dtype=int)
    pe = np.full(nspec, flux.shape[1]-1, dtype=int)
    indx = np.any(model_mask, axis=1) & np.logical_not(np.all(model_mask, axis=1))
    if np.any(indx):
        ps[indx], pe[indx] = np.atleast_2d(np.array([np.where(_m)[0][[0,-1]]
                                                        for _m in model_mask[indx]])).T
    pe += 1
    vsyst = np.array([-ppxf_vsyst(tpl_wave, wave[s:e], velscale, velscale_ratio=velscale_ratio)
                        for s, e in zip(ps, pe)])

    #-------------------------------------------------------------------
    # For debugging
    linear=False
#    linear=True
#    reject_boxcar=None
#    mdegree=0
    if starting_spectrum is None:
        starting_spectrum = 0
    else:
        fault[:starting_spectrum] = True
    #-------------------------------------------------------------------

    #-------------------------------------------------------------------
    # Fit each spectrum individually
    for i in range(starting_spectrum, nspec):
        if not np.any(model_mask[i,ps[i]:pe[i]]):
            # Nothing to fit!
            continue

        # Report progress
        print('Fitting spectrum: {0}/{1}'.format(i+1,nspec), end='\r')

        # Confirm that all templates and components are valid and
        # rearrange component arrays, if necessary
        valid_templates, _templates, _gas_template, component_map, _component, _vgrp, _sgrp, \
                _moments, _start, _constr_kinem \
                        = _validate_templates_components(templates, gas_template, component, vgrp,
                                                         sgrp, moments, model_mask[i,ps[i]:pe[i]],
                                                         start[i], _tpl_to_use[i,:], velscale,
                                                         velscale_ratio=velscale_ratio,
                                                         vsyst=vsyst[i], constr_kinem=constr_kinem)

        # Construct the parameter tying structure
        tied = ppxf_tied_parameters(_component, _vgrp, _sgrp, _moments)

        # Run the first fit.
        # NOTE: lsq_box is the default in ppxf 7.4.0, so no need to
        # define it here.
        if plot:
            plt.clf()
        try:
            # Note that use of constr_kinem requires method='capfit'; if
            # constr_templ is ever used, this would require a switch to
            # linear_method='lsq_lin'
            pp = ppxf.ppxf(_templates.T, flux[i,ps[i]:pe[i]], noise[i,ps[i]:pe[i]], velscale,
                           _start, velscale_ratio=velscale_ratio, plot=plot, moments=_moments,
                           degree=degree, mdegree=mdegree, lam=wave[ps[i]:pe[i]],
                           reddening=reddening, tied=tied,
                           constr_kinem={} if _constr_kinem is None else _constr_kinem,
                           mask=model_mask[i,ps[i]:pe[i]], vsyst=vsyst[i], component=_component,
                           gas_component=_gas_template, quiet=quiet, method='capfit',
                           linear=linear, linear_method='lsq_box')
        except Exception as e:
            if ppxf_faults == 'raise':
                raise e
            warnings.warn('pPXF fault: "{0}".  Logging fault and continuing.'.format(str(e)))
            fault[i] = True
            continue
        if plot:
            plt.show()

        # Reject 3-sigma outliers and refit, if requested by a provided
        # boxcar width
        if reject_boxcar is not None:
            # - Calculate residuals
            resid = flux[i,ps[i]:pe[i]] - pp.bestfit
            # - Select pixels included in the fit and not fit by
            # emission lines
            reject_pixels = list(set(pp.goodpixels)
                                    & set(np.arange(len(resid))[pp.gas_bestfit < 1e-6]))
            # - Calculate the 1-sigma confidence interval
            rms = calculate_noise(resid[reject_pixels], width=reject_boxcar)
            # - Reject pixels with > 3-sigma residuals
            model_mask[i,reject_pixels+ps[i]] &= (np.absolute(resid[reject_pixels]) < sigma_rej*rms)
        
            # Reorder the output; sets any omitted components to have
            # the starting values from the original input
            sol, err = _reorder_solution(pp.sol, pp.error, component_map, moments, start=start[i])

            # Confirm that all templates and components are valid and
            # rearrange component arrays, if necessary
            valid_templates, _templates, _gas_template, component_map, _component, _vgrp, _sgrp, \
                    _moments, _start, _constr_kinem \
                            = _validate_templates_components(templates, gas_template, component,
                                                             vgrp, sgrp, moments,
                                                             model_mask[i,ps[i]:pe[i]],
                                                             sol, _tpl_to_use[i,:], velscale,
                                                             velscale_ratio=velscale_ratio,
                                                             vsyst=vsyst[i],
                                                             constr_kinem=constr_kinem)

            # Construct the parameter tying structure
            tied = ppxf_tied_parameters(_component, _vgrp, _sgrp, _moments)

            # Refit using best-fit kinematics from previous fit as
            # initial guesses
            if plot:
                plt.clf()
            try:
                # Note that use of constr_kinem requires method='capfit'; if
                # constr_templ is ever used, this would require a switch to
                # linear_method='lsq_lin'
                pp = ppxf.ppxf(_templates.T, flux[i,ps[i]:pe[i]], noise[i,ps[i]:pe[i]], velscale,
                               _start, velscale_ratio=velscale_ratio, plot=plot, moments=_moments,
                               degree=degree, mdegree=mdegree, lam=wave[ps[i]:pe[i]],
                               reddening=reddening, tied=tied,
                               constr_kinem={} if _constr_kinem is None else _constr_kinem,
                               mask=model_mask[i,ps[i]:pe[i]], vsyst=vsyst[i],
                               component=_component, gas_component=_gas_template, quiet=quiet,
                               method='capfit', linear=linear, linear_method='lsq_box') 
            except Exception as e:
                if ppxf_faults == 'raise':
                    raise e
                warnings.warn('pPXF fault: "{0}".  Logging fault and continuing.'.format(str(e)))
                fault[i] = True
                continue
            if plot:
                plt.show()

        # Reorder the output; sets any omitted components to a default
        # value of -999.
        sol, err = _reorder_solution(pp.sol, pp.error, component_map, moments)

        # Save the results
        model[i,ps[i]:pe[i]] = pp.bestfit
        eml_model[i,ps[i]:pe[i]] = np.dot(pp.matrix[:,degree+1:][:,_gas_template],
                                          pp.weights[_gas_template])

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

    # Done
    print('Fitting spectrum: {0}/{0}'.format(nspec))

    return model, eml_model, model_mask, tpl_wgts, tpl_wgts_err, addcoef, multcoef, ebv, \
                    kininp, kin, kin_err, fault


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
        templates (`numpy.ndarray`_):
            Templates library to use for fitting.  Shape is
            (NTPLPIX,NTPL).
        gas_template (`numpy.ndarray`_):
            Boolean vector that selects the gas templates.  Shape is
            (NTPL,).
        wgts (`numpy.ndarray`_):
            The best-fitting template weights for each template in each
            spectrum.  Shape is (NSPEC,NTPL).
        component (`numpy.ndarray`_):
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


def emline_fitter_with_ppxf(tpl_wave, templates, wave, flux, noise, mask, velscale, velscale_ratio,
                            inp_component, gas_template, inp_moments, inp_start, vgrp=None,
                            sgrp=None, constr_kinem=None, degree=-1, mdegree=0, reddening=None,
                            reject_boxcar=101, tpl_to_use=None, binid=None, flux_binned=None,
                            noise_binned=None, mask_binned=None, x_binned=None, y_binned=None,
                            x=None, y=None, plot=False, quiet=False, debug=False, sigma_rej=3.,
                            ppxf_faults='flag'):
    r"""
    Main calling function for simultaneously fitting stellar-continuum and
    nebular emission lines in many spectra using pPXF.

    This is a generalization of a function originally provided by Xihan Ji
    and Michele Cappellari.

    The function *does not* fit the stellar kinematics; these are fixed
    during the fit (as provided in ``inp_start``) and should have resulted
    from a previous fit to the stellar-continuum alone; see e.g.,
    :class:`~mangadap.proc.ppxffit.PPXFFit`. The number of stellar kinematic
    moments must have been the same for all spectra, and there can only be
    one stellar component.

    The templates are expected to be an integer number of ``velscale_ratio``
    in length; see :func:`~mangadap.proc.ppxffit.PPXFFit.check_templates`.

    The fitting procedure can be performed for a single set of spectra, or a
    set of spectra that are remapped from an input set of binned spectra. The
    latter operation is selected by providing the binned flux, error, and
    mask arrays. The individual spectra can be mapped to a binned spectrum
    using binned and unbinned on-sky coordinates or by providing the bin
    indices directly (see argument description below).

    When the binned spectra are provided, the fitting procedure is as
    follows:

        - If a direct association of binned spectrum to remapped
          spectrum is not provided (via ``binid``), the Cartesian
          coordinates are used to identify the closest association
          between binned and remapped spectra.
        
        - The binned spectra are fit with the stellar components
          fixed to the provided kinematics in the starting value
          array and with all the gas templates part of a single
          kinematic component. This first fit iteration includes any
          rejection iterations according to the provided boxcar
          width.

        - The fit to the binned spectrum is then used to set (1) the
          single optimal stellar template and (2) the initial guess
          for the gas kinematics to use for the subsequent fits to
          each remapped spectrum.

        - Each spectrum is fit for the first time with the stellar
          kinematics fixed to the result for the associated binned
          spectrum and, again, with all the gas templates free but
          part of the same kinematic component. This fit iteration
          includes any rejection iterations according to the provided
          boxcar width.

        - Finally the spectra are fit without any rejection iteration and
          allowing the gas templates to be associated with multiple kinematic
          components as requested by the user (using the ``vgrp``, ``sgrp``,
          and ``constr_kinem`` arguments).

    When no binned spectra are provided, the procedure is virtually
    the same except the initial fit to the binned spectra is skipped.
    An initial fit to the spectra is performed to construct the
    optimal template, instead of basing the optimal template on the
    initial fit to the binned spectrum.

    To tie the kinematics assigned to each template, you have to assign them
    to velocity and/or velocity dispersion (sigma) groups using ``vgrp`` and
    ``sgrp``, respectively. These arrays are used to construct the ppxf
    ``tied`` argument, appropriate for each spectrum fit; see
    :func:`ppxf_tied_parameters`. To impose constraints on the component
    kinematics, use ``constr_kinem``; see the ppxf documentation.

    Each template is identified as a gas template using the provided
    ``gas_template`` argument.

    .. todo::

        - Allow mask(s) to be optional
        - Skip the last step if the gas are already part of the same
          kinematic component as dictated by the input tying data.

    Args:
        tpl_wave (`numpy.ndarray`_):
            Wavelength array for the spectral templates.  Shape is
            :math:`(N_{\rm tpl, pix},)`.
        templates (`numpy.ndarray`_):
            Templates library to use for fitting. Shape is
            :math:`(N_{\rm tpl},N_{\rm tpl, pix})`.
        wave (`numpy.ndarray`_):
            Wavelength vector.  Shape is :math:`(N_{\rm pix},)`.
        flux (`numpy.ndarray`_):
            Object spectra to fit. 
            Shape is :math:`(N_{\rm spec},N_{\rm pix})`.
        noise (`numpy.ndarray`_):
            Error in the object spectra to fit. Shape must match
            ``flux``.
        mask (`numpy.ndarray`_):
            Boolean vector that selects the pixels in the object
            spectra to fit (i.e., mask=True for pixels to fit and
            mask=False for pixels to ignore). Shape must match
            ``flux``. All pixels are fit by default.
        velscale (:obj:`float`):
            The pixel scale of the object spectra in km/s.
        velscale_ratio (:obj:`int`):
            The ratio between the object and template pixel scale; must
            be an integer.
        inp_component (`numpy.ndarray`_):
            Integer vector identifying the kinematic component for
            each template. Shape is :math:`(N_{\rm tpl},)`. The
            unique components *must* be 0..:math:`N_{\rm comp}-1`,
            without skipping any numbers.
        gas_template (`numpy.ndarray`_):
            Boolean vector that selects the gas templates. Shape is
            :math:`(N_{\rm tpl},)`.
        inp_moments (`numpy.ndarray`_):
            Integer vector with the number of kinematic moments for each
            component. Shape is :math:`(N_{\rm comp},)`. WARNING: These must
            currently all be 2.
        inp_start (:obj:`list`, `numpy.ndarray`_):
            The starting kinematics to use for each spectrum. Shape
            is :math:`(N_{\rm spec},)`; each element has shape
            :math:`(N_{\rm comp},)`, and each component element has
            shape :math:`(N_{\rm mom},)`. The object types allow each
            component to have its own number of moments, but each
            spectrum should have the same :math:`(N_{\rm comp},)` and
            :math:`(N_{\rm mom},)` for each component.
        vgrp (array-like, optional):
            The integer velocity group associated with each template.
            Shape is :math:`(N_{\rm tpl},)`. 
        sgrp (array-like, optional):
            The integer sigma group associated with each template.
            Shape is :math:`(N_{\rm tpl},)`. 
        constr_kinem (:obj:`dict`, optional):
            A dictionary with the constraints to apply to the kinematics of
            each component; see ``constr_kinem`` in ppxf.
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
            Size of the window for the rejection statistics. Should
            be an odd number and larger than 10. Default is 101. If
            None, no rejection iterations are performed.
        tpl_to_use (`numpy.ndarray`_, optional):
            Boolean vector selecting templates to consider during the
            fit. If None, all templates are used in the fit. If
            provided, the shape must be :math:`(N_{\rm bin},N_{\rm
            tpl})` when providing the binned data and :math:`(N_{\rm
            spec},N_{\rm tpl})` when no binned data are provided.
        binid (`numpy.ndarray`_, optional):
            Bin index associated with each spectrum.  Ignored if binned
            spectra are not provided.  If binned spectra are provided,
            and this is None, coordinates must be provided and they are
            used to associate each bin with a binned spectrum just based
            by proximity.  If provided, this associates each spectrum to
            the binned spectrum to use for the stellar kinematics.
            Shape is :math:`(N_{\rm spec},)`.
        flux_binned (`numpy.ndarray`_, optional):
            Binned spectra with previous fits to the stellar
            kinematics. See purpose above. Shape must be
            :math:`(N_{\rm bin},N_{\rm pix})`.
        noise_binned (`numpy.ndarray`_, optional):
            Error in the binned spectra. Shape must match
            ``flux_binned``.
        mask_binned (`numpy.ndarray`_, optional):
            Good-pixel mask for the binned spectra. Shape must match
            ``flux_binned``.
        x_binned (`numpy.ndarray`_, optional):
            On-sky bin x coordinate; shape is :math:`(N_{\rm bin},)`.
        y_binned (`numpy.ndarray`_, optional):
            On-sky bin y coordinate; shape is :math:`(N_{\rm bin},)`.
        x (`numpy.ndarray`_, optional):
            On-sky spectrum x coordinates; shape is :math:`(N_{\rm
            spec},)`.
        y (`numpy.ndarray`_, optional):
            On-sky spectrum y coordinates; shape is :math:`(N_{\rm
            spec},)`.
        plot (:obj:`bool`, optional):
            Show the pPXF fit plot at each iteration.
        quiet (:obj:`bool`, optional):
            Suppress output to the terminal.
        debug (:obj:`bool`, optional):
            Run in debugging mode.  Currently, all this does is perform
            the initial setup and then return empty vectors of the
            correct shape.  No fits are performed.
        ppxf_faults (:obj:`str`, optional):
            Dictates how exceptions raised by ppxf are treated.
            Allowed values are:

                - ``'flag'``: Log the fault and continue. Any
                  spectrum that faults is flagged as such in the last
                  object returned by the method (see below).
                - ``'raise'``: Re-raise the exception returned by
                  ppxf.
   
    Returns:
        :obj:`tuple`:  The following arrays are returned:

            #. The best-fitting model for each spectrum.
               Shape is :math:`(N_{\rm spec},N_{\rm pix})`.

            #. The best-fitting emission-line-only model for each
               spectrum. Shape is :math:`(N_{\rm spec},N_{\rm pix})`.
               
            #. Boolean array that is True for all spectral pixels
               included in the fit spectrum. Shape is :math:`(N_{\rm
               spec},N_{\rm pix})`.

            #. Best-fitting weight for each template in each
               spectrum. Shape is :math:`(N_{\rm spec}, N_{\rm
               tpl})`.

            #. Error in the best-fitting template weights.
               Shape is :math:`(N_{\rm spec}, N_{\rm tpl})`.

            #. Coefficients of the additive polynomial for each
               spectrum; None if not fit. Shape is :math:`(N_{\rm
               spec}, d+1`, where :math:`d` is ``degree``.
            
            #. Coefficients of the multiplicative polynomial for
               each spectrum; None if not fit. Shape is
               :math:`(N_{\rm spec}, m`, where :math:`m` is
               ``mdegree``.

            #. Best-fitting reddening values for each spectrum; None
               if not fit. Shape is :math:`(N_{\rm spec},)`.

            #. Input kinematics for each fit. Shape is :math:`(N_{\rm
               spec},\sum_i k_i)`, where :math:`k_i` is the number of
               kinematic moments for component :math:`i`.

            #. Best-fit kinematics for each fit. Shape is
               :math:`(N_{\rm spec},\sum_i k_i)`, where :math:`k_i`
               is the number of kinematic moments for component
               :math:`i`.

            #. Formal error in the best-fit kinematics for each
               spectrum. Shape is :math:`(N_{\rm spec},\sum_i k_i)`,
               where :math:`k_i` is the number of kinematic moments
               for component :math:`i`.

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
    ngas = np.sum(gas_template)
    if ngas != templates.shape[0] and np.amax(inp_component[~gas_template]) != 0:
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
    # NOTE: Reminder that `and` takes precedence over `or`, so this is the same
    # has having parenthesis around the first two and last two conditionals in
    # the if statement.
    if ntpl == ngas and not np.all(inp_moments == 2) \
            or ntpl != ngas and not np.all(inp_moments[1:] == 2):
        raise ValueError('Currently, moments for all gas components must be set to 2.')
    nkin = np.sum(np.absolute(inp_moments))

    # Check that the number of components is correct
    if len(inp_component) != ntpl:
        raise ValueError('Must assign each template to a kinematic component.')
    if vgrp is not None and len(vgrp) != ntpl:
        raise ValueError('If provided, must assign each template to a velocity group.')
    if sgrp is not None and len(sgrp) != ntpl:
        raise ValueError('If provided, must assign each template to a sigma group.')

    # Check that the input kinematic constraints are correct
    ncomp = np.amax(inp_component)+1
    if constr_kinem is not None:
        ckeys = list(constr_kinem.keys())
        for k in ['eq', 'ineq']:
            a = 'A_{0}'.format(k)
            b = 'b_{0}'.format(k)
            if np.sum(np.isin([a, b], ckeys)) not in [0, 2]:
                raise ValueError('Must provide both {0} and {1} in constr_kinem.'.format(a,b))
            if a in ckeys:
                if constr_kinem[a].shape[1] != nkin:
                    raise ValueError('{0} has incorrect number of columns.'.format(a))
                if constr_kinem[a].shape[0] != constr_kinem[b].size:
                    raise ValueError('Number of rows in {0} does not match {1}.'.format(a,b))

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
        # - All gas templates in a single component (no parameter tying necessary)
        component, moments, start = _ppxf_component_setup(inp_component, gas_template, inp_start,
                                                          single_gas_component=True)
        _, _, _mask, binned_tpl_wgts, _, _, _, _, _, binned_kin, _, fault \
                    = _fit_iteration(tpl_wave, templates, wave, flux_binned, noise_binned,
                                     velscale, start, moments, component, gas_template,
                                     tpl_to_use=tpl_to_use, reject_boxcar=reject_boxcar,
                                     velscale_ratio=velscale_ratio, degree=degree, mdegree=mdegree,
                                     reddening=reddening, mask=mask_binned, plot=plot, quiet=quiet,
                                     sigma_rej=sigma_rej, ppxf_faults=ppxf_faults) 

        # - Create a new template set that includes the optimal stellar
        #   template for each *binned* spectrum based on the
        #   best-fitting weights from above.  Propagate the changes to
        #   the components and groups.
        stellar_wgts, _templates, _gas_template, _tpl_to_use, _component, _vgrp, _sgrp \
                    = _combine_stellar_templates(templates, gas_template,
                                                 binned_tpl_wgts, inp_component, vgrp, sgrp)

        # - Restructure the relevant arrays to prepare for the fits to
        #   the *individual* spectra.  Propagate the changes to the
        #   components and groups.
        fault = fault[_binid]
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
        fault = np.zeros(nspec, dtype=bool)
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
    #    parameter tying necessary). When binned spectra have been fit
    #    previously, this only fits spectra that are components of a bin with
    #    more than one spectrum!
    component, moments, start = _ppxf_component_setup(_component, _gas_template, _start,
                                                      single_gas_component=True)

    # Restructure the kinematics array so that it can accept only
    # fitting the individual spectra that themselves consisted of an
    # entire bin
    kin = start.copy()
    kin = np.array([ k for k in kin ]).reshape(nspec,-1)

    indx = component_of_bin & np.logical_not(fault)
    _, _, model_mask[indx,:], tpl_wgts[indx,:], _, _, _, _, _, kin[indx,:], _, fault[indx] \
            = _fit_iteration(tpl_wave, _templates, wave, flux[indx,:], noise[indx,:], velscale,
                             start[indx], moments, component, _gas_template,
                             tpl_to_use=_tpl_to_use[indx,:], reject_boxcar=reject_boxcar,
                             velscale_ratio=velscale_ratio, degree=degree, mdegree=mdegree,
                             reddening=reddening, mask=model_mask[indx,:], plot=plot, quiet=quiet,
                             sigma_rej=sigma_rej, ppxf_faults=ppxf_faults) #, starting_spectrum=1133)

    if mode == 'noBins':
        # If binned spectra were not provided, the fit above is the
        # first fit to the individual spectra using all
        # stellar-continuum templates.  Use the result of that fit to
        # construct the optimal template for each spectrum, and
        # restructure the relevant component arrays.
        stellar_wgts, _templates, _gas_template, _tpl_to_use, _component, _vgrp, _sgrp \
                    = _combine_stellar_templates(_templates, _gas_template, tpl_wgts, _component,
                                                 _vgrp, _sgrp)
        tpl_wgts = np.ones((nspec, nspec+np.sum(_gas_template)), dtype=float)

    # Second fit to the individual spectra
    # - Use first fit to to reset the starting estimates for the gas kinematics
    all_stellar_moments = np.sum(np.absolute(moments[:-1]))
    gas_start = kin[:,all_stellar_moments:]
    component, moments, start = _ppxf_component_setup(_component, _gas_template, _start,
                                                      gas_start=gas_start)

    # - Refit without rejection but with the tying and constraints in place
    kin = np.zeros((nspec,nkin), dtype=float)
    indx = np.logical_not(fault)
    model_flux[indx,:], model_eml_flux[indx,:], model_mask[indx,:], tpl_wgts[indx,:], \
        tpl_wgts_err, _addcoef, _multcoef, _ebv, kininp[indx,:], kin[indx,:], \
        kin_err[indx,:], fault[indx] \
                = _fit_iteration(tpl_wave, _templates, wave, flux[indx,:], noise[indx,:], velscale,
                                 start[indx], moments, component, _gas_template,
                                 tpl_to_use=_tpl_to_use[indx,:], reject_boxcar=None,
                                 velscale_ratio=velscale_ratio, degree=degree, mdegree=mdegree,
                                 reddening=reddening, vgrp=_vgrp, sgrp=_sgrp,
                                 constr_kinem=constr_kinem, mask=model_mask[indx,:], plot=plot,
                                 quiet=quiet, sigma_rej=sigma_rej, ppxf_faults=ppxf_faults)

    # Save the low-order continuum coefficients
    if degree > -1:
        addcoef[indx,:] = _addcoef
    if mdegree > 0:
        multcoef[indx,:] = _multcoef
    if reddening is not None:
        ebv[indx] = _ebv

    # - Use the single output weight to renormalize the individual
    #   stellar template weights (only one of the weights for the
    #   non-gas templates should be non-zero). Stellar-template
    #   weight errors are always returned as 0; all weights for
    #   failed fits are set to 0.

    _tpl_wgts = stellar_wgts * np.sum(tpl_wgts[:,np.invert(_gas_template)], axis=1)[:,None]
    _tpl_wgts[np.logical_not(indx),:] = 0.
    _indx = np.ix_(indx,gas_template)
    _tpl_wgts[_indx] = tpl_wgts[np.ix_(indx,_gas_template)]
    _tpl_wgts_err = np.zeros((nspec,ntpl), dtype=float)
    _tpl_wgts_err[_indx] = tpl_wgts_err[:,_gas_template]

    return model_flux, model_eml_flux, model_mask, _tpl_wgts, _tpl_wgts_err, addcoef, multcoef, \
                    ebv, kininp, kin, kin_err, _binid, fault

