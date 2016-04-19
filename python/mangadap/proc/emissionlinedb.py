# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Container class for a database of emission-line parameters, as well as
support classes and functions.

*License*:
    Copyright (c) 2016, SDSS-IV/MaNGA Pipeline Group
        Licensed under BSD 3-clause license - see LICENSE.rst

*Source location*:
    $MANGADAP_DIR/python/mangadap/proc/emissionlinedb.py

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
            try:
                from configparser import ConfigParser
            except ImportError:
                warnings.warn('Unable to import configparser!  Beware!')
            try:
                from configparser import ExtendedInterpolation
            except ImportError:
                warnings.warn('Unable to import ExtendedInterpolation!  Some configurations will fail!')
        else:
            try:
                from ConfigParser import ConfigParser
            except ImportError:
                warnings.warn('Unable to import ConfigParser!  Beware!')
            try:
                from ConfigParser import ExtendedInterpolation
            except ImportError:
                warnings.warn('Unable to import ExtendedInterpolation!  Some configurations will fail!')
        
        import os.path
        from os import environ
        import glob
        import numpy
        
        from ..config.defaults import default_dap_source
        from ..util.idlutils import airtovac
        from ..util.yanny import yanny
        from ..par.parset import ParSet
        from .util import _select_proc_method

.. warning::

    Because of the use of the ``ExtendedInterpolation`` in
    `configparser.ConfigParser`_,
    :func:`available_emission_line_databases` is not python 2
    compiliant.

*Class usage examples*:
    To define an emission line::

        from mangadap.proc.emissionlinedb import EmissionLinePar
        p = EmissionLinePar(44, 'Ha', 6564.632, action='f', line='l',
                            amp=1.0, vel=0.0, sig=10., mode='f')

    More often, however, you will want to define an emission-line
    database using an SDSS parameter file.  To do so, you can use one of
    the default set of available emission-line databases (see
    :func:`available_emission_line_databases`)::
    
        from mangadap.proc.emissionlinedb import EmissionLineDB
        p = EmissionLineDB('STRONG')

    The above call requires that the ``$MANGADAP_DIR`` environmental
    variable is set.  If it is not, you can define it's location, as
    in::

        from mangadap.proc.emissionlinedb import EmissionLineDB
        p = EmissionLineDB('STRONG', dapsrc='/path/to/dap/source')

    Finally, you can create your own :class:`mangadap.util.yanny`
    parameter file with your own emission lines to fit.  Example files
    are provided in ``$MANGADAP_DIR/external/emission_lines`` with a
    companion ``README`` file.  With your own file, you have to point to
    the file using :class:`EmissionLineDBDef`, which you can then pass
    to :class:`EmissionLineDB`::

        from mangadap.proc.emissionlinedb import EmissionLineDBDef
        from mangadap.proc.emissionlinedb import EmissionLineDB
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

.. _configparser.ConfigParser: https://docs.python.org/3/library/configparser.html#configparser.ConfigParser

"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import sys
import warnings
if sys.version > '3':
    long = int
    try:
        from configparser import ConfigParser
    except ImportError:
        warnings.warn('Unable to import configparser!  Beware!', ImportWarning)
    try:
        from configparser import ExtendedInterpolation
    except ImportError:
        warnings.warn('Unable to import ExtendedInterpolation!  Some configurations will fail!',
                      ImportWarning)
else:
    try:
        from ConfigParser import ConfigParser
    except ImportError:
        warnings.warn('Unable to import ConfigParser!  Beware!', ImportWarning)
    try:
        from ConfigParser import ExtendedInterpolation
    except ImportError:
        warnings.warn('Unable to import ExtendedInterpolation!  Some configurations will fail!',
                      ImportWarning)

import os.path
from os import environ
import glob
import numpy

from ..config.defaults import default_dap_source
from ..util.idlutils import airtovac
from ..util.yanny import yanny
from ..par.parset import ParSet, ParDatabase
from .util import _select_proc_method
from .spectralfeaturedb import available_spectral_feature_databases, SpectralFeatureDBDef

__author__ = 'Kyle B. Westfall'
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
        line (str) : (**Optional**) Type of line, which can be either
            ``l`` for a line or ``dN`` for a doublet.  For doublets, the
            ``N`` indicates the *index* of the other line in the
            doublet.  The line to which the doublet is tied should have
            ``line='l'``; for example, if emission line with ``index=4``
            has ``line='d3'``, then the emission line with ``index=3``
            must have ``line='l'``.  Default is ``'l'``.
        amp (float) : (**Optional**) Relative intensity of the gas
            emission (positive) or absorption (negative) lines with
            respect to the doublet.  Therefore, this should most often
            be unity if ``line='l'`` and indicate the ratio of line
            **intensity** if ``line='dN'``.   Default is 1.0.
        vel (float) : (**Optional**) Guess for the velocity offset with
            respect to the galaxy systemic velocity (km/s).  Default is
            0.0.
        sig (float) : (**Optional**) Guess for the velocity dispersion
            (km/s). Default is 10.0.
        mode (float) : (**Optional**) Fitting mode for the line, which
            can be either ``'f'`` to fit the line independently or
            ``'tN'`` to set both the velocity and dispersion to be tied
            to a fitted line (``mode='f'``) with ``index=N``.  One can
            also force only the velocities to be tied using ``'vN'`` or
            only the velocity dispersions to be tied using ``'sN'``.
            Default is ``'f'``.
    """
    def __init__(self, index, name, restwave, action=None, line=None, amp=None, vel=None, sig=None,
                 mode=None):
        
        in_fl = [ int, float ]
        action_options = [ 'i', 'f', 'm', 's']
        
        l = [ name, action, line, mode ]
        for i in range(len(l)):
            if l[i] is None:
                continue
            l[i] = l[i].strip()

        pars =     [ 'index', 'name', 'restwave', 'action', 'line', 'amp', 'vel', 'sig', 'mode' ]
        values =   [   index,   l[0],   restwave,     l[1],   l[2],   amp,   vel,   sig,   l[3] ]
        defaults = [    None,   None,       None,      'f',    'l',   1.0,   0.0,  10.0,    'f' ]
        options =  [    None,   None,       None, action_options, None, None, None, None,  None ]
        dtypes =   [     int,    str,      in_fl,      str,    str, in_fl, in_fl, in_fl,    str ]

        ParSet.__init__(self, pars, values=values, defaults=defaults, options=options,
                        dtypes=dtypes)
        self._check()


    def _check(self):
        """
        Check the parameter list:

            - *line* must be either ``'l'`` or ``'dN'``.
            - Amplitude has to be larger than zero.
            - Velocity dispersion has to be greater than 0.
            - *mode* must be either ``'f'``, ``'tN'``, ``'vN'``, ``'sN'``

        .. todo::
            - Add check to __setitem__()?

        Raises:
            ValueError: Raised if one of the conditions above are not
                met.
        """
        if self.data['line'][0] not in ['l', 'd']:
            raise ValueError('Line must be either independent (l) or a double (d).')
        if self.data['action'] != 'i' and not self.data['amp'] > 0:
            warnings.warn('Emission-line amplitudes must be larger than 0.  Ignoring line with ' \
                          'index {0}'.format(self.data['index']))
            self.data['action'] = 'i'
        if not self.data['sig'] > 0:
            raise ValueError('The line velocity dispersion must be larger than 0.')
        if self.data['mode'][0] not in ['f', 't', 'v', 's']:
            raise ValueError('Mode must be either independent (f), fit with the velocity and '
                             'velocity dispersion tied (t), only the velocity tied (v), or only '
                             'the velocity dispersion tied (s).')


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
    :func:`mangadap.proc.available_spectral_feature_databases`.

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
        version (str): Version number
        database (str): Keyword of the selected database to use.
        neml (int): Number of emission lines in the database

    """
    def __init__(self, database_key, emldb_list=None, dapsrc=None):

        # TODO: The approach here (read using yanny, set to par
        # individually, then covert back to record array using
        # ParDatabase) is stupid...

        self.version = '1.0'

        # Get the details of the selected database
        self.database = _select_proc_method(database_key, SpectralFeatureDBDef,
                                            method_list=emldb_list,
                                            available_func=available_emission_line_databases,
                                            dapsrc=dapsrc)
        
        # Check that the database exists
        if not os.path.isfile(self.database['file_path']):
            raise FileNotFoundError('Database file {0} does not exist!'.format(
                                                                    self.database['file_path']))

        # Read the yanny file
        par = yanny(self.database['file_path'])
        if len(par['DAPEML']['index']) == 0:
            raise ValueError('Could not find DAPEML entries in {0}!'.self.database['file_path'])

        # Setup the array of emission line database parameters
        self.neml = len(par['DAPEML']['index'])
        parlist = []
        for i in range(self.neml):
            invac = par['DAPEML']['waveref'][i] == 'vac'
            parlist += [ EmissionLinePar(par['DAPEML']['index'][i], par['DAPEML']['name'][i],
                    par['DAPEML']['lambda'][i] if invac else airtovac(par['DAPEML']['lambda'][i]),
                                         action=par['DAPEML']['action'][i],
                                         line=par['DAPEML']['line'][i],
                                         amp=par['DAPEML']['ampl'][i],
                                         vel=par['DAPEML']['vel'][i],
                                         sig=par['DAPEML']['sig'][i],
                                         mode=par['DAPEML']['mode'][i]) ]

        ParDatabase.__init__(self, parlist)

        # Ensure that all indices are unique
        if len(numpy.unique(self.data['index'])) != self.neml:
            raise ValueError('Indices in {0} database are not all unique!'.format(
                                                                            self.database['key']))
        


    


