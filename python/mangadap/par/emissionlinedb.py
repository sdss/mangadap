# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Container class for a database of emission-line parameters, as well as
support classes and functions.

*License*:
    Copyright (c) 2016, SDSS-IV/MaNGA Pipeline Group
        Licensed under BSD 3-clause license - see LICENSE.rst

*Source location*:
    $MANGADAP_DIR/python/mangadap/par/emissionlinedb.py

*Imports and python version compliance*:
    ::

        from __future__ import division
        from __future__ import print_function
        from __future__ import absolute_import
        from __future__ import unicode_literals

        import sys
        import warnings
        if sys.version > '3':
            long = int
        
        import os.path
        from os import environ
        import glob
        import numpy

        from pydl.goddard.astro import airtovac
        from pydl.pydlutils.yanny import yanny
        from .parset import ParSet, ParDatabase
        from .spectralfeaturedb import available_spectral_feature_databases, SpectralFeatureDBDef
        from ..proc.util import select_proc_method
        from ..config.defaults import default_dap_source

*Class usage examples*:
    To define an emission line::

        from mangadap.par.emissionlinedb import EmissionLinePar
        p = EmissionLinePar(44, 'Ha', 6564.632, action='f', line='l',
                            flux=1.0, vel=0.0, sig=10., mode='f')

    More often, however, you will want to define an emission-line
    database using an SDSS parameter file.  To do so, you can use one of
    the default set of available emission-line databases (see
    :func:`available_emission_line_databases`)::
    
        from mangadap.par.emissionlinedb import EmissionLineDB
        p = EmissionLineDB('STRONG')

    The above call requires that the ``$MANGADAP_DIR`` environmental
    variable is set.  If it is not, you can define it's location, as
    in::

        from mangadap.par.emissionlinedb import EmissionLineDB
        p = EmissionLineDB('STRONG', dapsrc='/path/to/dap/source')

    Finally, you can create your own `SDSS-style parameter file`_ with
    your own emission lines to fit.  Example files are provided in
    ``$MANGADAP_DIR/data/emission_lines`` with a companion
    ``README`` file.  With your own file, you have to point to the file
    using :class:`EmissionLineDBDef`, which you can then pass to
    :class:`EmissionLineDB`::

        from mangadap.par.emissionlinedb import EmissionLineDBDef
        from mangadap.par.emissionlinedb import EmissionLineDB
        d = EmissionLineDBDef(key='USER',
                              file_path='/path/to/parameter/file',
                              in_vacuum=True)
        p = EmissionLineDB('USER', emldb_list=d)

    The reference frame of the emission-line wavelengths must be defined
    as either vacuum or air, using 'in_vacuum'.  It is expected that the
    object spectra to be fit are calibrated to vacuum wavelengths.  If
    'in_vacuum' is false, this class will use
    :func:`mangadap.util.idlutils.airtovac` to convert the emission-line
    wavelengths to vacuum.

*Revision history*:
    | **17 Mar 2016**: Original implementation by K. Westfall (KBW)
    | **11 May 2016**: (KBW) Switch to using `pydl.pydlutils.yanny`_ and
        `pydl.goddard.astro.airtovac`_ instead of internal functions
    | **13 Jul 2016**: (KBW) Include log_bounded, blueside, and redside
        in database.

.. _pydl.pydlutils.yanny: http://pydl.readthedocs.io/en/stable/api/pydl.pydlutils.yanny.yanny.html
.. _pydl.goddard.astro.airtovac: http://pydl.readthedocs.io/en/stable/api/pydl.goddard.astro.airtovac.html#pydl.goddard.astro.airtovac
.. _SDSS-style parameter file: http://www.sdss.org/dr12/software/par/

"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import sys
import warnings
if sys.version > '3':
    long = int

import os.path
from os import environ
import glob
import numpy

from pydl.goddard.astro import airtovac
from pydl.pydlutils.yanny import yanny
from .parset import ParSet, ParDatabase
from .spectralfeaturedb import available_spectral_feature_databases, SpectralFeatureDBDef
from ..proc.util import select_proc_method
from ..config.defaults import default_dap_source

# Add strict versioning
# from distutils.version import StrictVersion


class EmissionLinePar(ParSet):
    """
    Parameter object that defines a set of emission-line parameters used
    by various algorithms in the DAP.

    .. todo::
        - Specify these algorithms
        - provide some basic printing functions for user-level
          interaction

    See :class:`mangadap.par.parset.ParSet` for attributes and raised
    exceptions.

    Args:
        index (int) : An index used to refer to the line in the *line*
            and *mode* attributes.
        name (str) : A name for the line.
        restwave (float) : The rest wavelength of the line in angstroms
            *in vacuum*.
        action (str) : (**Optional**) Describes how the line should be
            treated.  Possible values are:

            - ``'i'``: ignore the line, as if the line were commented
              out.
            - ``'f'``': fit the line and/or mask the line when fitting
              the stellar continuum.
            - ``'m'``': mask the line when fitting the stellar continuum
              but do NOT fit the line itself
            - ``'s'``: defines a sky line that should be masked.  When
              masked, the wavelength of the line is NOT adjusted for the
              redshift of the object spectrum.

            Default is ``'f'``.
        flux (float) : (**Optional**) **Relative** flux of the emission
            (positive) or absorption (negative) lines.  This should most
            often be unity if ``line='l'`` and indicates the ratio of
            line flux if ``line='dN'``.  Default is 1.0.
        mode (float) : (**Optional**) Fitting mode for the line.
            Possible values are:
                - ``'f'``: Fit the line independently of all others in
                  its own window.
                - ``'wN'``: Fit the line with untied parameters, but use
                  a window that includes both this line and the line
                  with index N.
                - ``'xN'``: Fit the line with its flux tied to the line
                  with index N.
                - ``'vN'``: Fit the line with the velocity tied to the
                  line with index N.
                - ``'sN'``: Fit the line with the velocity dispersion
                  tied to the line with index N.
                - ``'kN'``: Fit the line with the velocity and velocity
                  dispersion tied to the line with index N.
                - ``'aN'``: Fit the line with the flux, velocity, and
                  velocity dispersion tied to the line with index N.
            Default is ``'f'``.
        profile (str): (**Optional**)  The class definition of the
            profile shape.  This is kept as a string until it is used.
            Once it is used, it is converted to the class name using::
        
                eval(profile)
            
            This line will fail if the profile type has not been
            defined.  Default is ``'GaussianLineProfile'``.
        ncomp (int): (**Optional**) The number of components (number of
            separate line profiles) to use when fitting the line.
            Default is 1.
        output_model (bool): (**Optional**) Flag to include the
            best-fitting model of the line in the emission-line model
            spectrum. Default is True.
        par (numpy.ndarray) : (**Optional**) A list of the initial guess
            for the line profile parameters.  The number of parameters
            must match the struct declaration at the top of the file.
            The initial parameters are automatically adjusted to provide
            any designated flux ratios, and the center is automatically
            adjusted to the provided redshift for the spectrum.  For
            example, for a GaussianLineProfile, this is typically set to
            ``[1.0, 0.0, 100.0]``.
        fix (numpy.ndarray): (**Optional**) A list of flags for fixing
            the input guess parameters during the fit.  Use 0 for a free
            parameter, 1 for a fixed parameter.  The parameter value is
            only fixed **after** adjusted in the flux and or center
            based on the redshift and the implied tied parameters.  For
            a free set of parameters using a GaussianLineProfile, this
            is set to ``[ 0, 0, 0 ]``.
        lobnd (numpy.ndarray): (**Optional**) A list of lower bounds for
            the parameters.  For each parameter, use None to indicate no
            lower bound.  For a GaussianLineProfile with positive flux
            and standard deviation, this is set to ``[ 0.0, None, 0.0
            ]``.
        hibnd (numpy.ndarray): (**Optional**) A list of upper bounds for
            the parameters.  For each parameter, use None to indicate no
            upper bound.  For a GaussianLineProfile with maximum
            standard deviation of 500 km/s, this is set to ``[ None,
            None, 500.0 ]``.
        log_bnd (numpy.ndarray): (**Optional**) A list of flags used
            when determining if a fit parameter is near the imposed
            boundary.  If true, the fraction of the boundary range used
            is done in logarithmic, not linear, separation.
        blueside (numpy.ndarray): (**Optional**) A two-element vector
            with the starting and ending wavelength for a bandpass
            blueward of the emission line, which is used to set the
            continuum level near the emission line when calculating the
            equivalent width.
        redside (numpy.ndarray): (**Optional**) A two-element vector
            with the starting and ending wavelength for a bandpass
            redward of the emission line, which is used to set the
            continuum level near the emission line when calculating the
            equivalent width.
    """
    def __init__(self, index, name, restwave, action=None, flux=None, mode=None, profile=None,
                 ncomp=None, output_model=None, par=None, fix=None, lobnd=None, hibnd=None,
                 log_bnd=None, blueside=None, redside=None):
        
        in_fl = [ int, float ]
        action_options = [ 'i', 'f', 'm', 's']
        arr_like = [ numpy.ndarray, list ]
        
        l = [ name, action, mode ]
        for i in range(len(l)):
            if l[i] is None:
                continue
            l[i] = l[i].strip()

        pars =     [ 'index', 'name', 'restwave', 'action', 'flux', 'mode', 'profile', 'ncomp',
                     'output_model', 'par', 'fix', 'lobnd', 'hibnd', 'log_bnd', 'blueside',
                     'redside' ]
        values =   [   index,   l[0],   restwave,     l[1],   flux,   l[2],   profile,   ncomp,
                       output_model, par,   fix,   lobnd,   hibnd, log_bnd, blueside, redside ]
        defaults = [ None, None, None, 'f', 1.0, 'f', 'GaussianLineProfile', 1, True, None, None,
                     None, None, None, None, None ]
        options =  [ None, None, None, action_options,  None,  None, None, None, None, None, None,
                     None, None, None, None, None ]
        dtypes =   [ int, str, in_fl, str, in_fl, str, str, int, bool, arr_like, arr_like,
                     arr_like, arr_like, arr_like, arr_like, arr_like ]

        ParSet.__init__(self, pars, values=values, defaults=defaults, options=options,
                        dtypes=dtypes)
        self._check()


    def _check(self):
        """
        Check the parameter list:

            - Amplitude has to be larger than zero.
            - *mode* must be either ``'f'``, ``'wN'``, ``'xN'``,
              ``'vN'``, ``'sN'``, ``'kN'``, ``'aN'``

        .. todo::
            - Add check to __setitem__()?
            - Add check of profile type

        Raises:
            ValueError: Raised if one of the conditions above are not
                met.
        """
        if self.data['action'] != 'i' and not self.data['flux'] > 0:
            warnings.warn('Emission-line fluxes must be larger than 0.  Ignoring line with ' \
                          'index {0}'.format(self.data['index']))
            self.data['action'] = 'i'
        if self.data['mode'][0] not in ['f', 'w', 'x', 'v', 's', 'k', 'a']:
            raise ValueError('Mode must be either independent (f), fit independently but in the ' \
                             'same window (w), with a tied flux ratio (x), with a tied velocity ' \
                             '(v), with a tied velocity dispersion (s), with both kinematics ' \
                             'tied (k), or with the fluxes and kinematics tied (a).')
        # Check the lengths of all the arrays
        npar = len(self.data['par'])
        if numpy.any(numpy.array([len(self.data['fix']), len(self.data['lobnd']),
                                  len(self.data['hibnd']), len(self.data['log_bnd'])]) != npar):
            raise ValueError('Number of parameters must be the same for par, fix, lobnd, hibnd, '
                             'and log_bnd.')
        if numpy.any(numpy.array([ len(self.data['blueside']), len(self.data['redside'])]) != 2):
            raise ValueError('Bandpasses must be two-element vectors!')


def available_emission_line_databases(dapsrc=None):
    """
    Return the list of database keys and file names for the available
    emission-line databases.  The currently available libraries are:
    
    +-------------+-----+----------------------------------+
    |         KEY |   N | Description                      |
    +=============+=====+==================================+
    |    STANDARD |  62 |  Original line list with nearly  |
    |             |     |  all strong and weak lines       |
    +-------------+-----+----------------------------------+
    |      STRONG |  13 | A subset of mostly strong lines. |
    +-------------+-----+----------------------------------+
    |    EXTENDED |  21 |  Include HeI/II, SIII, Hgam-Heps |
    +-------------+-----+----------------------------------+

    This is a simple wrapper for
    :func:`mangadap.par.spectralfeaturedb.available_spectral_feature_databases`.

    Args:
        dapsrc (str): (**Optional**) Root path to the DAP source
            directory.  If not provided, the default is defined by
            :func:`mangadap.config.defaults.default_dap_source`.

    Returns:
        list: An list of :func:`SpectralFeatureDBDef` objects, each of
        which defines a unique emission-line database.

    .. todo::
        - Add backup function for Python 2.
        - Somehow add a python call that reads the databases and
          constructs the table for presentation in sphinx so that the
          text above doesn't have to be edited with changes in the
          available databases.
        
    """
    return available_spectral_feature_databases('emission_lines', dapsrc=dapsrc)


class EmissionLineDB(ParDatabase):
    """
    Basic container class for the database of emission-line parameters.
    See :class:`mangadap.parset.ParDatabase` for additional attributes.

    Args:
        database_key (str): Keyword selecting the database to use.
        emldb_list (list): (**Optional**) List of
            :class:`SpectralFeatureDBDef` objects that defines the
            unique key for the database and the path to the source SDSS
            parameter file.
        dapsrc (str): (**Optional**) Root path to the DAP source
            directory.  If not provided, the default is defined by
            :func:`mangadap.config.defaults.default_dap_source`.

    Attributes:
        database (str): Keyword of the selected database to use.
        neml (int): Number of emission lines in the database

    """
    def __init__(self, database_key, emldb_list=None, dapsrc=None):

        # TODO: The approach here (read using yanny, set to par
        # individually, then covert back to record array using
        # ParDatabase) is stupid...

        # Get the details of the selected database
        self.database = select_proc_method(database_key, SpectralFeatureDBDef,
                                           method_list=emldb_list,
                                           available_func=available_emission_line_databases,
                                           dapsrc=dapsrc)
        
        # Check that the database exists
        if not os.path.isfile(self.database['file_path']):
            raise FileNotFoundError('Database file {0} does not exist!'.format(
                                                                    self.database['file_path']))

        # Read the yanny file
#        par = yanny(self.database['file_path'])
        par = yanny(filename=self.database['file_path'], raw=True)
        if len(par['DAPEML']['index']) == 0:
            raise ValueError('Could not find DAPEML entries in {0}!'.format(
                                                                    self.database['file_path']))

        # Setup the array of emission line database parameters
        self.neml = len(par['DAPEML']['index'])
        parlist = []
        for i in range(self.neml):
            invac = par['DAPEML']['waveref'][i] == 'vac'
            # Convert None's to +/- inf
            lobnd = [ -numpy.inf if p == 'None' else float(p) \
                            for p in par['DAPEML']['lower_bound'][i]]
            hibnd = [ numpy.inf if p == 'None' else float(p) \
                            for p in par['DAPEML']['upper_bound'][i] ]
            parlist += [ EmissionLinePar(par['DAPEML']['index'][i], par['DAPEML']['name'][i],
                    par['DAPEML']['lambda'][i] if invac else airtovac(par['DAPEML']['lambda'][i]),
                                         action=par['DAPEML']['action'][i],
                                         flux=par['DAPEML']['relative_flux'][i],
                                         mode=par['DAPEML']['mode'][i],
                                         profile=par['DAPEML']['profile'][i],
                                         ncomp=par['DAPEML']['ncomp'][i],
                                         output_model=bool(par['DAPEML']['output_model'][i]),
                                         par=par['DAPEML']['par'][i],
                                         fix=par['DAPEML']['fix'][i],
                                         lobnd=lobnd, hibnd=hibnd,
                                         log_bnd=par['DAPEML']['log_bounded'][i],
                                blueside=par['DAPEML']['blueside'][i] if invac \
                                        else airtovac(numpy.array(par['DAPEML']['blueside'][i])),
                                redside=par['DAPEML']['redside'][i] if invac \
                                        else airtovac(numpy.array(par['DAPEML']['redside'][i])) ) ]

        ParDatabase.__init__(self, parlist)

        # Ensure that all indices are unique
        if len(numpy.unique(self.data['index'])) != self.neml:
            raise ValueError('Indices in {0} database are not all unique!'.format(
                                                                            self.database['key']))
        


    


