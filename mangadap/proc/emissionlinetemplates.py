# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
r"""

Class that constructs a set of emission-line templates, primarily for
use with :class:`mangadap.proc.sasuke.Sasuke`.

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

Revision history
----------------

    | **08 Sep 2017**: Originally pulled from
        :mod:`mangadap.proc.sasuke` by K. Westfall (KBW)

----

.. include license and copyright
.. include:: ../copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../links.rst
"""

import os
import time
import logging

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
    database, primarily for use in :class:`mangadap.proc.sasuke.Sasuke`.

    The templates are constructed based on the constraints provided by
    the emission-line database.  See
    :class:`mangadap.par.emissionlinedb.EmissionLinePar` for the
    structure of each row in the database and an explanation for each of
    its columns.  The selected profile type for each line **must** have
    a ``parameters_from_moments`` method that returns the parameters of
    the line provided the first three moments (moments 0, 1, and 2).

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

        The construction of templates for use with :class:`Sasuke` does
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
        wave (numpy.ndarray):
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
        flux (numpy.ndarray):
            The template spectra.  Shape is :math:`(N_{\rm tpl},N_{\rm
            pix})`.
        tpli (numpy.ndarray):
            The index of the template containing each emission line.
            Any emission-lines with ``tpli==-1`` means that the emission
            line was not included in any template, which should only
            occur for lines with ``action==i`` in the database.  Shape
            is :math:`(N_{\rm eml},)`.
        comp (numpy.ndarray):
            The component number assigned to each template.  Templates
            with the same component number are forced to have the same
            velocity and velocity dispersion by pPXF.  Shape is
            :math:`(N_{\rm tpl},)`.
        vgrp (numpy.ndarray):
            The velocity group assigned to each template.  Templates in
            the same velocity group have their velocity parameters tied
            in pPXF, but the velocity dispersion parameters are
            independent.  Shape is :math:`(N_{\rm tpl},)`.
        sgrp (numpy.ndarray):
            The velocity disperison (sigma) group assigned to each
            template.  Templates in the same sigma group have their
            velocity dispersion parameters tied in pPXF, but the
            velocity parameters are independent.  Shape is
            :math:`(N_{\rm tpl},)`.
        eml_sigma_inst (numpy.ndarray):
            The instrumental dispersion (km/s) at the rest wavelength of
            each emission line.  This is mostly used to aid the velocity
            dispersion corrections determined by :class:`Sasuke`.  Shape
            is :math:`(N_{\rm eml},)`.
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
                                               bounds_error=False, fill_value=-1)
        self.log = log
        self.base = base

#        self.dv = numpy.full(self.wave.size, spectrum_velocity_scale(wave), dtype=float) if log \
#                    else astropy.constants.c.to('km/s').value*spectral_coordinate_step(wave)/wave

        self.emldb = None           # Original database
        self.flux_density = None    # Templates are in flux per Angstrom, not flux per pixel
        self.ntpl = None            # Number of templates
        self.flux = None            # Template fluxes
        self.tpli = None            # Template associated with each emission line
        self.comp = None            # Kinematic component associated with each template
        self.vgrp = None            # Velocity group associated with each template
        self.sgrp = None            # Sigma group associated with each template
        self.eml_sigma_inst = None  # Instrumental dispersion at the center of each line

        # If emission-line database is provided, use it to build the templates
        if emldb is not None:
            self.build_templates(emldb, flux_density=flux_density, loggers=loggers, quiet=quiet)


    def _tied_index(self, i, primary=False):
        r"""
        Return the index of the line to which this one is tied.
        
        The result may be that this line is tied to one that is also
        tied to a second line.  If that's the case, the ``primary``
        keyword can be use to trace the parameter tying back to the
        independent line.

        Args:
            i (:obj:`int`):
                The index of the line in the database.
            primary (:obj:`bool`, optional):
                Trace the line tying all the way back to the independent
                (primary) line.

        Returns:
            :obj:`int`: The index of the line to which this one is tied.

        Raises:
            ValueError:
                Raised if the primary option is selected and the line
                does not trace back to a primary line.  This represents
                a poorly constructed emission-line database and should
                be avoided!

        """
        db_rows = numpy.arange(self.emldb.size)
        indx = db_rows[self.emldb['index'] == int(self.emldb['mode'][i][1:])][0]
        if not primary:
            return indx
        max_iter = 100
        j = 0
        while self.emldb['mode'][indx] != 'f' and j < max_iter:
            indx = db_rows[self.emldb['index'] == int(self.emldb['mode'][indx][1:])][0]
            j+=1
        if j == max_iter:
            raise ValueError('Line {0} (index={1}) does not trace back to a primary line!'.format(
                                i, self.emldb['index'][i]))
        return indx


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

            The construction of templates for use with :class:`Sasuke`
            does *not* allow one to tie fluxes while leaving the
            velocities and/or velocity dispersions as independent.

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
        ignore_line = self.emldb['action'] != 'f'

        # The total number of templates to construct is the number of
        # lines in the database minus the number of lines with mode=aN
        tied_all = numpy.array([m[0] == 'a' for m in self.emldb['mode']])
        self.ntpl = self.emldb.size - numpy.sum(ignore_line) - numpy.sum(tied_all)

        # Initialize the components
        self.comp = numpy.zeros(self.ntpl, dtype=int)-1
        self.vgrp = numpy.zeros(self.ntpl, dtype=int)-1
        self.sgrp = numpy.zeros(self.ntpl, dtype=int)-1

        # All the primary lines go into individual templates, kinematic
        # components, velocity groups, and sigma groups
        self.tpli = numpy.zeros(self.emldb.size, dtype=int)-1
        primary_line = (self.emldb['mode'] == 'f') & numpy.invert(ignore_line)
        nprimary = numpy.sum(primary_line)
        self.tpli[primary_line] = numpy.arange(nprimary)
        self.comp[:nprimary] = numpy.arange(nprimary)
        self.vgrp[:nprimary] = numpy.arange(nprimary)
        self.sgrp[:nprimary] = numpy.arange(nprimary)

        finished = primary_line | ignore_line
        while numpy.sum(finished) != self.emldb.size:
            # Find the indices of lines that are tied to finished lines
            start_sum = numpy.sum(finished)
            for i in range(self.emldb.size):
                if finished[i]:
                    continue
                indx = self._tied_index(i)
                if not finished[indx]:
                    continue

                finished[i] = True

                # Mode=a: Line is part of an existing template
                if self.emldb['mode'][i][0] == 'a':
                    self.tpli[i] = self.tpli[indx]
                # Mode=k: Line is part of a different template but an
                # existing kinematic component
                if self.emldb['mode'][i][0] == 'k':
                    self.tpli[i] = numpy.amax(self.tpli)+1
                    self.comp[self.tpli[i]] = self.comp[self.tpli[indx]]
                    self.vgrp[self.tpli[i]] = self.vgrp[self.tpli[indx]]
                    self.sgrp[self.tpli[i]] = self.sgrp[self.tpli[indx]]
                # Mode=v: Line is part of a different template and
                # kinematic component with an untied sigma, but tied to
                # an existing velocity group
                if self.emldb['mode'][i][0] == 'v':
                    self.tpli[i] = numpy.amax(self.tpli)+1
                    self.comp[self.tpli[i]] = numpy.amax(self.comp)+1
                    self.sgrp[self.tpli[i]] = numpy.amax(self.sgrp)+1
                    self.vgrp[self.tpli[i]] = self.vgrp[self.tpli[indx]]
                # Mode=s: Line is part of a different template and
                # kinematic component with an untied velocity, but tied
                # to an existing sigma group
                if self.emldb['mode'][i][0] == 's':
                    self.tpli[i] = numpy.amax(self.tpli)+1
                    self.comp[self.tpli[i]] = numpy.amax(self.comp)+1
                    self.vgrp[self.tpli[i]] = numpy.amax(self.vgrp)+1
                    self.sgrp[self.tpli[i]] = self.sgrp[self.tpli[indx]]

            # If the loop ends up with the same number of parsed lines
            # that it started with, there must be an error in the
            # construction of the input database.
            if start_sum == numpy.sum(finished):
                raise ValueError('Unable to parse the input database.  Check tying parameters.')

        # Debug:
        if numpy.any(self.comp < 0) or numpy.any(self.vgrp < 0) or numpy.any(self.sgrp < 0):
            raise ValueError('Templates without an assigned component.  Check the input database.')


    def check_database(self, emldb):
        r"""

        Check that the provided emission-line database can be used with
        the :class:`EmissionLineTemplates` class.  Most checks are
        performed by
        :func:`mangadap.proc.spectralfitting.EmissionLineFit.check_emission_line_database`.

        Additional checks specific to :class:`EmissionLineTemplates`
        are:

            - Any lines with mode `w` is treated as `f` and a warning is
              provided.
            - The :class:`EmissionLineTemplates` object *cannot* be used
              with mode `x`; any lines with this mode will cause a
              ValueError to be raised..

        This function does *not* check if the initial parameters
        provided by the database are consistent with other elements in
        the database because they are not used to construct the
        templates.

        Args:
            emldb (:class:`mangadap.par.emissionlinedb.EmissionLineDB`):
                Emission-line database.

        Raises:
            TypeError: Raised if the provided object is not an instance
                of :class:`mangadap.par.emissionlinedb.EmissionLineDB`.
            ValueError: Raised if any line has a mode of `x` or if the
                database does not provide a valid definition for any
                templates.
            NameError: Raised if a defined profile type is not known.
        """
        EmissionLineFit.check_emission_line_database(emldb, wave=self.wave, check_par=False)

        # Check that no lines only tie the fluxes
        if numpy.any([m[0] == 'x' for m in emldb['mode']]):
            raise ValueError('Cannot tie only fluxes in an EmissionLineTemplates object.')

        # Warn user of any lines with mode=w
        if numpy.any([m[0] == 'w' for m in emldb['mode']]):
            warnings.warn('Any line with mode=w treated the same as mode=f.')


    def build_templates(self, emldb, flux_density=True, loggers=None, quiet=False):
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
            numpy.ndarray: Returns 4 arrays: (1) the set of templates
            with shape :math:`N_{\rm tpl}\times N_{\rm wave}`, (2) the
            kinematic component assignement for each template, (3) the
            velocity group associated with each template, and (4) the
            sigma group assocated with each template.

        """
        #---------------------------------------------------------------
        # Initialize the reporting
        if loggers is not None:
            self.loggers = loggers
        self.quiet = quiet
        self.flux_density = flux_density

        #---------------------------------------------------------------
        # Check the database can be used with this class
        self.check_database(emldb)
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
                       'Number of emission-line velocity groups: {0}'.format(
                                                                    numpy.amax(self.vgrp)+1))
            log_output(self.loggers, 1, logging.INFO,
                       'Number of emission-line sigma groups: {0}'.format(
                                                                    numpy.amax(self.sgrp)+1))

        # Get the instrumental dispersion at the center of each line
        self.eml_sigma_inst = self.sigma_inst(self.emldb['restwave'])
#        pyplot.plot(self.sigma_inst.x, self.sigma_inst.y, lw=2)
#        pyplot.plot(self.sigma_inst.x, 2.5*astropy.constants.c.to('km/s').value/self.sigma_inst.x/2.355)
#        pyplot.scatter(self.emldb['restwave'], self.eml_sigma_inst, marker='.', s=100, lw=0,
#                       color='C2')
#        pyplot.show()

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
        _flux = self.emldb['flux'] / dl if self.flux_density else self.emldb['flux']
        # Dispersion in pixel units
        _sigma = self.eml_sigma_inst * self.emldb['restwave'] \
                        / astropy.constants.c.to('km/s').value / dl

        # Construct the templates
        pix = numpy.arange(self.wave.size)
        self.flux = numpy.zeros((self.ntpl,self.wave.size), dtype=float)
        for i in range(self.ntpl):
            # Find all the lines associated with this template:
            index = numpy.arange(self.emldb.size)[self.tpli == i]
            # Add each line to the template
            for j in index:
                # Declare an instance of the desired profile
                profile = eval('lineprofiles.'+self.emldb['profile'][j])()
                # Use the first three moments of the line to set the
                # parameters
                p = profile.parameters_from_moments(_flux[j], _restwave[j], _sigma[j])
#                print(p)
#                f = profile(pix, p)
#                d = numpy.power(self.base, _wave + _dw/2) - numpy.power(self.base, _wave - _dw/2)
#                print(numpy.sum(d*f))
#                print(numpy.sum(f))
#                print(self.eml_sigma_inst[j])
#                print(numpy.amax(f))
#                print(numpy.amax(f)/2)
#                pyplot.plot(self.wave, f)
#                pyplot.show()
#                # Warn the user if the line is undersampled.
#                srt = numpy.argsort(numpy.absolute(v))
#                if self.eml_sigma_inst[j]/self.dv[srt[0]] < 0.9:
#                    warnings.warn('{0} line is undersampled!'.format(self.emldb['name'][j]))
                # Add the line to the flux in this template
                self.flux[i,:] += profile(pix, p)

                #line_spec = profile(pix,p)
                #line_spec *= (_flux[j]/numpy.sum(line_spec))
                #print(numpy.sum(line_spec), _flux[j], numpy.sum(line_spec)/_flux[j])
                #self.flux[i,:] += line_spec

        return self.flux, self.comp, self.vgrp, self.sgrp




