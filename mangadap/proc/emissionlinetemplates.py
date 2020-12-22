# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
r"""

Class that constructs a set of emission-line templates, primarily for
use with :class:`~mangadap.proc.sasuke.Sasuke`.

Class usage examples
--------------------

To construct emission-line templates, you need a wavelength vector,
the instrumental resolution at which to construct the templates, and
an emission line database (see
:class:`mangadap.par.emissionlinedb.EmissionLineDB`). A simple
construction would be::
        
    # Imports
    import numpy
    from mangadap.par.emissionlinedb import EmissionLineDB
    from mangadap.proc.emissionlinetemplates import EmissionLineTemplates

    wave = numpy.logspace(*(numpy.log10([3600,10300]), 4563)
    sigma_inst = 30                     # Instrumental resolution in km/s

    tpl = EmissionLineTemplates(wave, sigma_inst, emldb=EmissionLineDB.from_key('ELPMILES'))

----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

import os
import time
import logging

from IPython import embed

import numpy
from scipy import interpolate

import astropy.constants

from ..util.sampling import spectrum_velocity_scale, spectral_coordinate_step
from ..util.log import log_output
from ..util import lineprofiles
from .spectralfitting import EmissionLineFit

# For debugging
from matplotlib import pyplot


class EmissionLineTemplates:
    r"""
    Construct a set of emission-line templates based on an emission-line
    database, primarily for use in the :class:`~mangadap.proc.sasuke.Sasuke`
    emission-line fitter.

    The templates are constructed based on the constraints provided by
    the emission-line database.  See
    :class:`~mangadap.par.emissionlinedb.EmissionLinePar` for the
    structure of each row in the database and an explanation for each of
    its columns.
    
    Only lines with ``action=f`` are included in any template.  The array
    :attr:`tpli` provides the index of the template that contains each
    line in the emission-line database.  Lines that are not assigned to
    any template --- either because they do not have ``action=f`` or their
    center lies outside the wavelength range in :attr:`wave` --- are
    given an index of -1.

    Only lines with ``mode=a`` (i.e., flux, velocity, and velocity
    dispersion are all tied) are included in the same template.

    Lines with tied velocities are assigned the same velocity component
    (:attr:`vgrp`) and lines with the tied velocity dispersions are
    assigned the same sigma component (:attr:`sgrp`).

    .. warning::

        The construction of templates for use with
        :class:`~mangadap.proc.sasuke.Sasuke` does
        *not* allow one to tie fluxes while leaving the velocities
        and/or velocity dispersions as independent.

    Args:
        wave (array-like):
            A single wavelength vector with the wavelengths for the
            template spectra.  The wavelengths are expected to be sample
            either linearly or geometrically (see :attr:`log`).
        sigma_inst (:obj:`float`, array-like):
            The single value or value as a function of wavelength for
            the instrumental dispersion (km/s) to use for the template
            construction.
        log (:obj:`bool`, optional):
            Flag that the wavelengths have been sampled geometrically.
        base (:obj:`float`, optional):
            Base for the geometric sampling.
        emldb (:class:`mangadap.par.emissionlinedb.EmissionLineDB`, optional):
            Emission-line database that is parsed to setup which lines
            to include in the same template because they are modeled as
            having the same velocity, velocity dispersion and flux
            ratio.  If not provided, no templates are constructed in
            instantiation; to build the templates using an existing
            instantiation, use :func:`build_templates`.
        flux_density (:obj:`bool`, optional):
            Return spectrum in units of flux density (flux per
            angstrom).  Default is to return the spectrum in units of
            flux per pixel.
        loggers (:obj:`list`, optional):
            List of `logging.Logger`_ objects to log progress; ignored
            if quiet=True.  Logging is done using
            :func:`mangadap.util.log.log_output`.  Default is no
            logging.  This can be reset in some methods.
        quiet (:obj:`bool`, optional):
            Suppress all terminal and logging output.

    Attributes:
        wave (`numpy.ndarray`_):
            Array with the wavelength (angstroms) of each pixel for all
            the constructed templates.  Shape is :math:`(N_{\rm pix},)`.
        sigma_inst (`scipy.interpolate.interp1d`_):
            The object used to interpolate the instrumental dispersion
            (km/s) at the rest wavelength of each spectral line.
        log (:obj:`bool`):
            Flag that the spectrum is sampled geometrically.
        base (:obj:`float`):
            Base for the geometric sampling.
        emldb (:class:`mangadap.par.emissionlinedb.EmissionLineDB`):
            Database with the emission-line parameters.  Shape is
            :math:`(N_{\rm pix},)`.
        ntpl (:obj:`int`):
            Total number of templates.
        flux (`numpy.ndarray`_):
            The template spectra.  Shape is :math:`(N_{\rm tpl},N_{\rm
            pix})`.
        tpli (`numpy.ndarray`_):
            The index of the template containing each emission line.
            Any emission-lines with ``tpli==-1`` means that the emission
            line was not included in any template, which should only
            occur for lines with ``action==i`` in the database.  Shape
            is :math:`(N_{\rm eml},)`.
        comp (`numpy.ndarray`_):
            The kinematic component assigned to each template.  Templates
            with the same component number are forced to have the same
            velocity and velocity dispersion by pPXF.  Shape is
            :math:`(N_{\rm tpl},)`.
        vgrp (`numpy.ndarray`_):
            The velocity group assigned to each template.  Templates in
            the same velocity group have their velocity parameters tied
            in pPXF, but the velocity dispersion parameters are
            independent.  Shape is :math:`(N_{\rm tpl},)`.
        sgrp (`numpy.ndarray`_):
            The velocity disperison (sigma) group assigned to each
            template.  Templates in the same sigma group have their
            velocity dispersion parameters tied in pPXF, but the
            velocity parameters are independent.  Shape is
            :math:`(N_{\rm tpl},)`.
        A_ineq (`numpy.ndarray`_):
            A matrix used to constrain the kinematic parameters between
            kinematic components. See the ``constr_kinem`` for pPXF (version
            7.0 and later). Shape is :math:`(N_{\rm constr}, N_{\rm comp})`;
            i.e., the number of constraints applied by the number of
            kinematic components.
        b_ineq (`numpy.ndarray`_):
            A vector used to constrain the kinematic parameters between
            kinematic components. See the ``constr_kinem`` for pPXF (version
            7.0 and later). Shape is :math:`(N_{\rm constr},)`; i.e., the
            number of constraints applied.
        eml_sigma_inst (`numpy.ndarray`_):
            The instrumental dispersion (km/s) at the rest wavelength
            of each emission line. This is mostly used to aid the
            velocity dispersion corrections determined by
            :class:`~mangadap.proc.sasuke.Sasuke`. Shape is
            :math:`(N_{\rm eml},)`.
        loggers (:obj:`list`):
            List of `logging.Logger`_ objects to log progress; ignored
            if quiet=True.  Logging is done using
            :func:`mangadap.util.log.log_output`.
        quiet (:obj:`bool`):
            Suppress all terminal and logging output.
    """
    def __init__(self, wave, sigma_inst, log=True, base=10, emldb=None, flux_density=True,
                 loggers=None, quiet=False):

        self.loggers=None
        self.quiet=None

        self.wave = numpy.asarray(wave)
        if len(self.wave.shape) != 1:
            raise ValueError('Provided wavelengths must be a single vector.')

        _sinst = numpy.full(self.wave.size, sigma_inst, dtype=float) \
                    if isinstance(sigma_inst, (int, float)) else numpy.asarray(sigma_inst)
        if _sinst.shape != self.wave.shape:
            raise ValueError('Provided sigma_inst must be a single number or a vector with the'
                             'same length as the wavelength vector.')
        self.sigma_inst = interpolate.interp1d(self.wave, _sinst, assume_sorted=True,
                                               bounds_error=False, fill_value=0)
        self.log = log
        self.base = base

        self.emldb = None           # Original database
        self.line_flux = None       # Line flux in the relevant template
        self.flux_density = None    # Templates are in flux per Angstrom, not flux per pixel
        self.ntpl = None            # Number of templates
        self.flux = None            # Template fluxes
        self.tpli = None            # Template associated with each emission line
        self.comp = None            # Kinematic component associated with each template
        self.vgrp = None            # Velocity group associated with each template
        self.sgrp = None            # Sigma group associated with each template
        self.A_ineq = None          # Linear inequality constraint matrix for the kinematics
        self.b_ineq = None          # Linear inequality constraint vector for the kinematics
        self.eml_sigma_inst = None  # Instrumental dispersion at the center of each line

        # If emission-line database is provided, use it to build the templates
        if emldb is not None:
            self.build_templates(emldb, flux_density=flux_density, loggers=loggers, quiet=quiet)

    def _parse_emission_line_database(self):
        r"""
        Parse the input emission-line database; see
        :class:`mangadap.par.emissionlinedb.EmissionLinePar`.

        Only lines with `action=f` are included in any template.  The
        array :attr:`tpli` provides the index of the template that
        contains each line in the emission-line database.  Lines that
        are not assigned to any template --- either because they do not
        have `action=f` or their center lies outside the wavelength
        range in :attr:`wave` --- are given an index of -1.

        Only lines with `mode=a` (i.e., tie flux, velocity, and velocity
        dispersion) are included in the same template.

        Lines with tied velocities are assigned the same velocity
        component (:attr:`vgrp`) and lines with the tied velocity
        dispersions are assigned the same sigma component
        (:attr:`sgrp`).

        .. warning::

            The construction of templates for use with
            :class:`~mangadap.proc.sasuke.Sasuke` does *not* allow
            one to tie fluxes while leaving the velocities and/or
            velocity dispersions as independent.

        This is an entirely internal procedure, taking no arguments and
        only assigning results to self.

        Raises:
            ValueError: Raised if the parsing of the database leaves
                some templates without an assigned component, velocity
                group, and/or sigma group.  This implies an error in the
                construction of the emission-line database, not the code
                itself.

        """
        # Get the list of lines to ignore
        ignore_line = (self.emldb['action'] != 'f') | (self.emldb['restwave'] < self.wave[0]) \
                        | (self.emldb['restwave'] > self.wave[-1])
        if numpy.any(ignore_line) and not self.quiet:
            log_output(self.loggers, 1, logging.INFO,
                       'Lines ignored because of action or wavelength limits: {0}'.format(
                       ', '.join(self.emldb['name'][ignore_line])))

        # Determine which parameters are tied by equality:
        tie_eq = numpy.array([p is not None and '=' in p 
                              for p in self.emldb['tie_par'].ravel()]).reshape(
                                  self.emldb['tie_par'].shape)
        tie_eq[ignore_line,:] = False

        # Determine which parameters re tied by inequality:
        tie_ineq = numpy.array([p is not None and '=' not in p
                                for p in self.emldb['tie_par'].ravel()]).reshape(
                                  self.emldb['tie_par'].shape) & numpy.logical_not(tie_eq)
        tie_ineq[ignore_line,:] = False

        # Ensure the problem is solvable.
        # TODO: I'm sure there's a more robust way to do this...
        if numpy.any(numpy.all(tie_eq | tie_ineq, axis=0)):
            raise ValueError('There must be at least one line that has an independent velocity '
                             'and another one that has an independent velocity dispsersion; '
                             'however, these do not have to be satisfied by a single line.')

        # Check that any tied line includes the index of the line it's tied to.
        untied = self.emldb['tie_index'] < 0
        indx = numpy.any(tie_eq | tie_ineq, axis=1) & untied
        if numpy.any(indx):
            raise ValueError('Index of reference line unspecified for {0}.'.format(
                             ', '.join(self.emldb['name'][indx])))

        # TODO: Need a check that the user hasn't input, e.g., '=0.3' for
        # anything other than the flux.
        if numpy.any([f is not None and '=' in f and len(f.strip('=')) > 0 
                        for f in self.emldb['tie_par'][:,1:].ravel()]):
            raise NotImplementedError('DAP cannot yet impose equality constraints for kinematics. '
                                      'Kinematics can only be forced to be identical (using \'=\')'
                                      ' or constrained to be within a fractional tolerance.')

        # Set the normalized line flux.  All default to unity.
        self.line_flux = numpy.ones(self.emldb.neml, dtype=float)
        # Specify those that are tied by equality
        self.line_flux[tie_eq[:,0]] = [float(f.strip('=')) 
                                        for f in self.emldb['tie_par'][tie_eq[:,0]][:,0]]

        # Currently cannot deal with lines that have tied fluxes but
        # independent kinematics
        if numpy.any(tie_eq[:,0] & numpy.logical_not(numpy.all(tie_eq[:,1:], axis=1))):
            raise NotImplementedError('DAP cannot currently accept lines with tied fluxes but '
                                      'independent kinematics.')

        # The total number of templates to construct is the number of lines in
        # the database minus the number of lines to ignore and the number of
        # lines where all parameters are tied by equality
        tied_all = numpy.all(tie_eq, axis=1)
        self.ntpl = self.emldb.size - numpy.sum(ignore_line) - numpy.sum(tied_all)

        # Initialize the components and tied kinematic groups
        self.comp = numpy.full(self.ntpl, -1, dtype=int)
        self.vgrp = numpy.full(self.ntpl, -1, dtype=int)
        self.sgrp = numpy.full(self.ntpl, -1, dtype=int)
        self.tpli = numpy.full(self.emldb.size, -1, dtype=int)

        # Reference lines are those that are:
        #   - *not* ignored (action != 'i'),
        #   - completely unconstrained (emldb['tie_index'] is indefined;
        #     i.e., less than 0), and/or
        #   - lines that are only tied by inequalities,
        # All reference lines are in separate templates, kinematic components,
        # velocity groups, and sigma groups
        ref_line = (untied | (numpy.logical_not(untied)
                              & numpy.logical_not(numpy.any(tie_eq[:,1:], axis=1)))) \
                    & numpy.logical_not(ignore_line)
        nref = numpy.sum(ref_line)

        self.tpli[ref_line] = numpy.arange(nref)
        self.comp[:nref] = numpy.arange(nref)
        self.vgrp[:nref] = numpy.arange(nref)
        self.sgrp[:nref] = numpy.arange(nref)

        finished = ref_line | ignore_line
        while numpy.sum(finished) != self.emldb.size:
            # Find the indices of lines that are tied to finished lines
            start_sum = numpy.sum(finished)
            for i in range(self.emldb.size):
                if finished[i]:
                    continue
                indx = numpy.where(self.emldb['index'] == self.emldb['tie_index'][i])[0]
                if not finished[indx]:
                    continue

                finished[i] = True

                # All line properties are tied such that the line will be part
                # of an existing template
                if numpy.all(tie_eq[i]):
                    self.tpli[i] = self.tpli[indx]
                # Both kinematic parameters are tied such that the line is part
                # of a different template but part of an existing kinematic
                # component
                elif numpy.all(tie_eq[i,1:]):
                    self.tpli[i] = numpy.amax(self.tpli)+1
                    self.comp[self.tpli[i]] = self.comp[self.tpli[indx]]
                    self.vgrp[self.tpli[i]] = self.vgrp[self.tpli[indx]]
                    self.sgrp[self.tpli[i]] = self.sgrp[self.tpli[indx]]
                # Line is part of a different template and kinematic component
                # with an untied sigma, but tied to an existing velocity group
                elif tie_eq[i,1]:
                    self.tpli[i] = numpy.amax(self.tpli)+1
                    self.comp[self.tpli[i]] = numpy.amax(self.comp)+1
                    self.sgrp[self.tpli[i]] = numpy.amax(self.sgrp)+1
                    self.vgrp[self.tpli[i]] = self.vgrp[self.tpli[indx]]
                # Line is part of a different template and kinematic component
                # with an untied velocity, but tied to an existing sigma group
                elif tie_eq[i,2]:
                    self.tpli[i] = numpy.amax(self.tpli)+1
                    self.comp[self.tpli[i]] = numpy.amax(self.comp)+1
                    self.vgrp[self.tpli[i]] = numpy.amax(self.vgrp)+1
                    self.sgrp[self.tpli[i]] = self.sgrp[self.tpli[indx]]

            # If the loop ends up with the same number of parsed lines
            # that it started with, there must be an error in the
            # construction of the input database.
            if start_sum == numpy.sum(finished):
                raise ValueError('Unable to parse the input database.  Check tying parameters.')

        # All templates must have an assigned component, velocity group and
        # dispersion group. Otherwise, something has gone wrong or the input
        # database hasn't been constructed correctly.
        if numpy.any(self.comp < 0) or numpy.any(self.vgrp < 0) or numpy.any(self.sgrp < 0):
            raise ValueError('Templates without an assigned component.  Check the input database.')

        # If there are no inequalities set in the database, we're done
        if not numpy.any(tie_ineq):
            self.A_ineq = None
            self.b_ineq = None
            return

        # Gut check
        assert not numpy.any(tie_eq & tie_ineq), \
                'CODING ERROR: Components are tied with both equality and inequality constraints'

        # The DAP currently doesn't deal with any template weight (flux)
        # inequalities:
        if numpy.any(tie_ineq[:,0]):
            raise NotImplementedError('DAP currently does not allow inequality constraints on '
                                      'line fluxes.')

        # Check that lines in the same component do not have different
        # constraints.
        # TODO: I'm not actually sure that it's possible for this to fault, but
        # it's here just in case.
        ncomp = numpy.amax(self.comp)+1
        # Component for each line
        compi = numpy.full(self.emldb.size, -1, dtype=int)
        indx = numpy.logical_not(ignore_line)
        compi[indx] = self.comp[self.tpli[indx]]
        for i in range(ncomp):
            indx = (compi == i) & tie_ineq[:,1]
            if numpy.any(indx) and len(numpy.unique(self.emldb['tie_par'][indx,1])) != 1:
                raise ValueError('Inequality constraints for velocities of lines '
                                 + '{0} ({1}) '.format(self.emldb['index'][indx],
                                                       self.emldb['name'][indx])
                                 + 'must be identical because they are in the same '
                                 + 'kinematic component.')
            indx = (compi == i) & tie_ineq[:,2]
            if numpy.any(indx) and len(numpy.unique(self.emldb['tie_par'][indx,2])) != 1:
                raise ValueError('Inequality constraints for dispersions of lines '
                                 + '{0} ({1}) '.format(self.emldb['index'][indx],
                                                       self.emldb['name'][indx])
                                 + 'must be identical because they are in the same '
                                 + 'kinematic component.')

        # The number of constraints is twice the number of specified
        # inequalities in the database
        indx = compi != -1
        nconstr = 2 * numpy.sum(tie_ineq[indx,1:])

        # Get the fractional constraint on each parameter
        tie_comp_par = numpy.zeros((ncomp, 2), dtype=float)
        # NOTE: This loop is necessary because of recursive tying of lines. For
        # example, in the OII doublet, one of the lines has its velocity tied
        # to the H-alpha line and the other has its velocity and velocity
        # dispersion tied to the first line. To make sure that the velocity
        # dispersion of both OII lines is within some inequality constraint of
        # the H-alpha dispersion, you need to pick the correct constraint. The
        # checks above should catch if this recursion leads to, e.g., different
        # constraints on the dispersion, so below I pick the correct constraint
        # by picking the maximum value. This works because the `=` constraint
        # is interpreted first as having a 0.0 fractional constraint.
        for i in range(ncomp):
            indx = compi == i
            _par = numpy.array([0.0 if f is None or len(f.strip('=')) == 0 else float(f.strip('='))
                                for f in self.emldb['tie_par'][indx,1:].ravel()]).reshape(
                                        numpy.sum(indx),2)
            tie_comp_par[i,:] = numpy.amax(_par, axis=0)

        # The definition of the fractional constraint is such that the value
        # MUST be larger than one. Determine if any of the values don't meet
        # this constraint and raise an exception if so.
        # TODO: Do this when reading the EmissionLineDB...
        if numpy.any(numpy.any((tie_comp_par < 0) 
                               | ((tie_comp_par > 0) & numpy.logical_not(tie_comp_par > 1)),
                               axis=1)):
            raise ValueError('Any defined inequality constraints must be >1!')

        # NOTE: For reference, this didn't work (compared to the explicit loop
        # above) because it was just assigning the most recent parameter with
        # the same component. Depending on the order of the line list, this
        # meant that some fractional constraints were set to 0, leading to a
        # singluar A_ineq array.
#        tie_comp_par[compi[indx],:] = numpy.array([0.0 if f is None or len(f.strip('=')) == 0 
#                                                 else float(f.strip('='))
#                                            for f in self.emldb['tie_par'][indx,1:].ravel()
#                                            ]).reshape(numpy.sum(indx),2)

        # Set a vector indicating the ineq connections between components
        tie_comp = numpy.full(ncomp, -1, dtype=int)
        indx_match = self.emldb.tie_index_match()
        indx = (indx_match >= 0) & numpy.any(tie_ineq[:,1:], axis=1) & (compi != -1)
        tie_comp[compi[indx]] = compi[indx_match[indx]]

        # Build the inequality for each component
        # NOTE: Inequality constraint is *always* defined with b_ineq=0 for
        # what the DAP can currently accommodate.
        self.b_ineq = numpy.zeros(nconstr, dtype=float)
        self.A_ineq = numpy.zeros((nconstr,2*ncomp), dtype=float)
        constr = 0
        for i in range(ncomp):
            if tie_comp[i] < 0:
                continue
            for j in range(2):
                if tie_comp_par[i,j] > 1:
                    self.A_ineq[constr,2*tie_comp[i]+j] = 1/tie_comp_par[i,j]
                    self.A_ineq[constr,2*i+j] = -1.
                    constr += 1
                    self.A_ineq[constr,2*tie_comp[i]+j] = -tie_comp_par[i,j]
                    self.A_ineq[constr,2*i+j] = 1.
                    constr += 1

    def build_templates(self, emldb, flux_density=True, profile='FFTGaussianLSF', loggers=None,
                        quiet=False):
        r"""
        Build the set of templates for a given emission-line database.
        The function uses the current values in :attr:`wave` and
        :attr:`sigma_inst`.  Any existing templates from a previous call
        to :func:`build_templates` or from the object instantiation will
        be overwritten using the provided emission-line database.

        See :func:`check_database` for the requirements of the provided
        emission-line database, and see
        :func:`_parse_emission_line_database` for how the database is
        interpretted when constructing the templates.

        The function constructs and returns the following attributes:
        :attr:`flux`, :attr:`comp`, :attr:`vgrp`, and :attr:`sgrp` 

        .. todo::

            - Warn the user if any line is undersampled; i.e., the FWHM
              of the line is less than 2.1 or :math:`\sigma < 0.9`.

            - Warn the user if any line grouped in the same template
              falls outside the spectral range.

        Args:
            emldb (:class:`mangadap.par.emissionlinedb.EmissionLineDB`):
                Emission-line database.
            flux_density (:obj:`bool`, optional):
                Return spectrum in units of flux density (flux per
                angstrom).  Default is to return the spectrum in units
                of flux per pixel.
            loggers (:obj:`list`, optional):
                List of `logging.Logger`_ objects to log progress;
                ignored if quiet=True.  Logging is done using
                :func:`mangadap.util.log.log_output`.
            quiet (:obj:`bool`, optional):
                Suppress all terminal and logging output.

        Returns:
            :obj:`tuple`: Returns 6 `numpy.ndarray`_ objects: (1) the set of
            templates with shape :math:`N_{\rm tpl}\times N_{\rm wave}`, (2)
            the kinematic component assignement for each template, (3) the
            velocity group associated with each template, (4) the sigma group
            associated with each template, (5) the A matrix used to impose
            constraints on the kinematics, and (6) the b vector used to
            impose constraints on the kinematics; for the latter two objects,
            see the ``constr_kinem`` keyword argument for pPXF.
        """
        #---------------------------------------------------------------
        # Initialize the reporting
        if loggers is not None:
            self.loggers = loggers
        self.quiet = quiet
        self.flux_density = flux_density

        #---------------------------------------------------------------
        # Check the database can be used with this class
        EmissionLineFit.check_emission_line_database(emldb, wave=self.wave)
        # Save a pointer to the database
        self.emldb = emldb
        # Parse the database for construction of the templates
        self._parse_emission_line_database()

        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO,
                       'Number of emission lines to fit: {0}'.format(numpy.sum(self.tpli>-1)))
            log_output(self.loggers, 1, logging.INFO,
                       'Number of emission-line templates: {0}'.format(len(self.comp)))
            log_output(self.loggers, 1, logging.INFO,
                       'Number of emission-line kinematic components: {0}'.format(
                                                                    numpy.amax(self.comp)+1))
            log_output(self.loggers, 1, logging.INFO,
                       'Number of emission-line kinematic constraints: {0}'.format(
                                                0 if self.b_ineq is None else self.b_ineq.size))
            log_output(self.loggers, 1, logging.INFO,
                       'Number of emission-line velocity groups: {0}'.format(
                                                                    numpy.amax(self.vgrp)+1))
            log_output(self.loggers, 1, logging.INFO,
                       'Number of emission-line sigma groups: {0}'.format(
                                                                    numpy.amax(self.sgrp)+1))

        # Get the instrumental dispersion at the center of each line
        self.eml_sigma_inst = self.sigma_inst(self.emldb['restwave'])

        # Convert from wavelengths to pixel coordinates
        _wave = numpy.log(self.wave)/numpy.log(self.base) if self.log else self.wave
        _dw = spectral_coordinate_step(self.wave, log=self.log, base=self.base)
        _restwave = numpy.log(self.emldb['restwave'])/numpy.log(self.base) if self.log \
                            else self.emldb['restwave']

        # Rest wavelength in pixel units
        _restwave = (_restwave - _wave[0])/_dw
        # Flux to pixel units; less accurate when spectrum is
        # logarithmically binned
        dl = self.emldb['restwave']*(numpy.power(self.base,_dw/2)-numpy.power(self.base,-_dw/2)) \
                        if self.log else _dw
        _flux = self.line_flux / dl if self.flux_density else self.line_flux
        # Dispersion in pixel units
        _sigma = self.eml_sigma_inst * self.emldb['restwave'] \
                        / astropy.constants.c.to('km/s').value / dl

        # Construct the templates
        _profile = eval('lineprofiles.'+profile)()  # Declare an instance of the desired profile
        pix = numpy.arange(self.wave.size)
        self.flux = numpy.zeros((self.ntpl,self.wave.size), dtype=float)
        for i in range(self.ntpl):
            # Find all the lines associated with this template:
            index = numpy.arange(self.emldb.size)[self.tpli == i]
            # Add each line to the template
            for j in index:
                # Use the first three moments of the line to set the
                # parameters
                p = _profile.parameters_from_moments(_flux[j], _restwave[j], _sigma[j])
                # Add the line to the flux in this template
                self.flux[i,:] += _profile(pix, p)

        return self.flux, self.comp, self.vgrp, self.sgrp, self.A_ineq, self.b_ineq




