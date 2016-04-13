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
        from ..config.util import validate_emission_line_config
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
from ..config.util import validate_emission_line_config
from ..util.idlutils import airtovac
from ..util.yanny import yanny
from ..par.parset import ParSet
from .util import _select_proc_method

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
        action (str) : (Optional) Describes how the line should be
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
        line (str) : (Optional) Type of line, which can be either ``l``
            for a line or ``dN`` for a doublet.  For doublets, the ``N``
            indicates the *index* of the other line in the doublet.  The
            line to which the doublet is tied should have ``line='l'``;
            for example, if emission line with ``index=4`` has
            ``line='d3'``, then the emission line with ``index=3`` must
            have ``line='l'``.  Default is ``'l'``.
        amp (float) : (Optional) Relative intensity of the gas emission
            (positive) or absorption (negative) lines with respect to
            the doublet.  Therefore, this should most often be unity if
            ``line='l'`` and indicate the ratio of line **intensity** if
            ``line='dN'``.   Default is 1.0.
        vel (float) : (Optional) Guess for the velocity offset with
            respect to the galaxy systemic velocity (km/s).  Default is
            0.0.
        sig (float) : (Optional) Guess for the velocity dispersion
            (km/s). Default is 10.0.
        mode (float) : (Optional) Fitting mode for the line, which can
            be either ``'f'`` to fit the line independently or ``'tN'``
            to set both the velocity and dispersion to be tied to a
            fitted line (``mode='f'``) with ``index=N``.  One can also
            force only the velocities to be tied using ``'vN'`` or only
            the velocity dispersions to be tied using ``'sN'``.  Default
            is ``'f'``.

    """
    def __init__(self, index, name, restwave, action=None, line=None, amp=None, vel=None, sig=None,
                 mode=None):
        
        in_fl = [ int, float ]
        action_options = [ 'i', 'f', 'm', 's']

        n, a, l, m = name.strip(), action.strip(), line.strip(), mode.strip()

        pars =     [ 'index', 'name', 'restwave', 'action', 'line', 'amp', 'vel', 'sig', 'mode' ]
        values =   [   index,      n,   restwave,        a,      l,   amp,   vel,   sig,      m ]
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
        if not self.data['amp'] > 0:
            raise ValueError('The line amplitude must be larger than 0.')
        if not self.data['sig'] > 0:
            raise ValueError('The line velocity dispersion must be larger than 0.')
        if self.data['mode'][0] not in ['f', 't', 'v', 's']:
            raise ValueError('Mode must be either independent (f), fit with the velocity and '
                             'velocity dispersion tied (t), only the velocity tied (v), or only '
                             'the velocity dispersion tied (s).')


class EmissionLineDBDef(ParSet):
    """
    Class with parameters used to define an emission-line database.
    Options and defaults in ParSet base class are set to None.
    """
    def __init__(self, key, file_path, in_vacuum=None ):

        pars =     [ 'key', 'file_path', 'in_vacuum' ]
        values =   [   key,   file_path,   in_vacuum ]
        defaults = [  None,        None,        True ]
        dtypes =   [   str,         str,        bool ]

        ParSet.__init__(self, pars, values=values, defaults=defaults, dtypes=dtypes)


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

    .. warning::

        Function is currently only valid for Python 3.2 or greater!

    Args:
        dapsrc (str): (Optional) Root path to the DAP source
            directory.  If not provided, the default is defined by
            :func:`mangadap.config.defaults.default_dap_source`.

    Returns:
        list: An list of :func:`EmissionLineDBDef`.  objects, each of
        which defines a unique emission-line database.

    Raises:
        NotADirectoryError: Raised if the provided or default
            *dapsrc* is not a directory.
        FileNotFoundError: Raised if the database parameter file is not
            found.
        KeyError: Raised if the emission-line keywords are not all
            unique.
        NameError: Raised if either ConfigParser or
            ExtendedInterpolation are not correctly imported.  The
            latter is a *Python 3 only module*!

    .. todo::
        - Add backup function for Python 2.
        - Somehow add a python call that reads the databases and
          constructs the table for presentation in sphinx so that the
          text above doesn't have to be edited with changes in the
          available databases.
        
    """
    # Check the source directory exists
    dapsrc = default_dap_source() if dapsrc is None else str(dapsrc)
    if not os.path.isdir(dapsrc):
        raise NotADirectoryError('{0} does not exist!'.format(dapsrc))

    # Check the configuration files exist
    ini_files = glob.glob(dapsrc+'/python/mangadap/config/emission_lines/*.ini')
    if len(ini_files) == 0:
        raise IOError('Could not find any configuration files in {0} !'.format(
                      dapsrc+'/python/mangadap/config/emission_lines'))

    # Build the list of library definitions
    eml_databases = []
    for f in ini_files:
        # Read the config file
        cnfg = ConfigParser(environ, allow_no_value=True, interpolation=ExtendedInterpolation())
        cnfg.read(f)
        # Ensure it has the necessary elements to define the
        # emission-line database
        validate_emission_line_config(cnfg)
        # Append the definition of the emission-line database
        eml_databases +=  [ EmissionLineDBDef(key=cnfg['default']['key'],
                                              file_path=cnfg['default']['file_path'],
                                              in_vacuum=cnfg['default'].getboolean('in_vacuum'))]

    # Check the keywords of the libraries are all unique
    if len(numpy.unique( numpy.array([eml['key'] for eml in eml_databases]) )) \
            != len(eml_databases):
        raise KeyError('Emission-line database keywords are not all unique!')

    # Return the default list of databases
    return eml_databases


class EmissionLineDB:
    """
    Basic container class for the database of emission-line parameters.

    .. todo::
        - Need to figure out is it's better to have an array of
          EmissionLinePar objects, or if I should convert self.data to a
          numpy record array.

    Args:
        database_key (str): Keyword selecting the database to use.
        emldb_list (list): (Optional) List of :class:`EmissionLineDBDef`
            objects that define the parameters required to setup the
            emission-line database.  The *database_key* must select one
            of the objects in this list.
        dapsrc (str): (Optional) Root path to the DAP source directory.
            If not provided, the default is defined by
            :func:`mangadap.config.defaults.default_dap_source`.

    Attributes:
        version (str): Version number
        database (str): Keyword of the selected database to use.
        dbparset (:class:`EmissionLineDBDef`): Parameter set required to
            setup the emission-line database.
        neml (int): Number of emission lines in the database
        data (:class:`numpy.array`) : Array of :class:`EmissionLinePar`
            instances, one per emission line in the database.

    """
    def __init__(self, database_key, emldb_list=None, dapsrc=None):

        self.version = '1.0'

        # Define the library properties
        self.database = None
        self.dbparset = None
        self.neml = None
        self.data = None

        # Read and process the database
        self.read_database(database_key, emldb_list=emldb_list, dapsrc=dapsrc)
        

    def __getitem__(self, index):
        """Select a specified parameter set.

        Args:
            index (int): Index to return.  Must be in the appropriate
                range.

        Returns:
            :class:`EmissionLinePar`: Parameter set of the selected
            emission line.

        Raises:
            IndexError: Raised if the provided index is out of bounds.
        """
        if index < 0 or index >= self.neml:
            raise IndexError('Index out of range.')
        return self.data[index]


    def _check(self):
        """
        Check that the database is correctly defined:

            - All the indices must be unique.

        .. todo::
            - Other checks needed?

        Raises:
            ValueError: Raised if the indices are not all unique.
        """
        if len(numpy.unique( numpy.array([d['index'] for d in self.data]))) != self.neml:
            raise ValueError('Indices in {0} database are not all unique!'.format(self.database))

    
    def read_database(self, database_key, emldb_list=None, dapsrc=None):
        """
        Select and read the database from the provided list.  Used to
        set :attr:`database` and :attr:`dbparset`; see
        :func:`mangadap.proc.util._select_proc_method`.

        Args:
            database_key (str): Keyword selecting the database to use.
            emldb_list (list): (Optional) List of
                :class:`EmissionLineDBDef` objects that define the
                parameters required to setup the emission-line database.
                The *database_key* must select one of the objects in
                this list.
            dapsrc (str): (Optional) Root path to the DAP source
                directory.  If not provided, the default is defined by
                :func:`mangadap.config.defaults.default_dap_source`.
        """
        # Get the details of the selected database
        self.dbparset = _select_proc_method(database_key, EmissionLineDBDef, emldb_list=emldb_list,
                                            available_func=available_emission_line_databases,
                                            dapsrc=dapsrc)
        self.database = database_key
        
        # Check that the database exists
        if not os.path.isfile(self.dbparset['file_path']):
            raise FileNotFoundError('Database file {0} does not exist!'.format(
                                                                    self.dbparset['file_path']))

        # Read the yanny file
        par = yanny(self.dbparset['file_path'])
        if len(par['DAPEML']['index']) == 0:
            raise ValueError('Could not find DAPEML entries in {0}!'.self.dbparset['file_path'])

        # Setup the array of emission line database parameters
        self.neml = len(par['DAPEML']['index'])
        self.data = numpy.empty(self.neml, dtype=object)
        for i in range(self.neml):
            self.data[i] = EmissionLinePar(par['DAPEML']['index'][i],
                                           par['DAPEML']['name'][i],
                                    par['DAPEML']['lambda'][i] if self.dbparset['in_vacuum'] \
                                        else airtovac(par['DAPEML']['lambda'][i]),
                                           action=par['DAPEML']['action'][i],
                                           line=par['DAPEML']['line'][i],
                                           amp=par['DAPEML']['ampl'][i],
                                           vel=par['DAPEML']['vel'][i],
                                           sig=par['DAPEML']['sig'][i],
                                           mode=par['DAPEML']['mode'][i])

        # Check that the database is setup correctly
        self._check()


