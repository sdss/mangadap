# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Container class for the database of absorption-line indices to measure.

*License*:
    Copyright (c) 2015, SDSS-IV/MaNGA Pipeline Group
        Licensed under BSD 3-clause license - see LICENSE.rst

*Source location*:
    $MANGADAP_DIR/python/mangadap/proc/absorptionindexdb.py

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
        from ..config.util import validate_absorption_index_config
        from ..util.idlutils import airtovac
        from ..util.yanny import yanny
        from ..par.parset import ParSet
        from ..par.bandpassfilter import BandPassFilterPar
        from .util import _select_proc_method

.. warning::

    Because of the use of the ``ExtendedInterpolation`` in
    `configparser.ConfigParser`_,
    :func:`available_absorption_index_databases` is not python 2
    compiliant.

*Class usage examples*:
    Absorption-line index databases are defined using SDSS parameter
    files.  To define a database, you can use one of the default set of
    available absorption-line index databases (see
    :func:`available_absorption_index_databases`)::
    
        from mangadap.proc.absorptionindexdb import AbsorptionIndexDB
        p = AbsorptionIndexDB('LICK')

    The above call requires that the ``$MANGADAP_DIR`` environmental
    variable is set.  If it is not, you can define it's location, as
    in::

        from mangadap.proc.absorptionindexdb import AbsorptionIndexDB
        p = AbsorptionIndexDB('LICK', dapsrc='/path/to/dap/source')

    Finally, you can create your own SDSS parameter file (see
    :class:`mangadap.util.yanny`) with your own absorption-line indices
    to use.  Example files are provided in
    ``$MANGADAP_DIR/external/absorption_indices`` with a companion
    ``README`` file.  With your own file, you have to point to the file
    using :class:`AbsorptionIndexDBDef`, which you can then pass to
    :class:`AbsorptionIndexDB`::

        from mangadap.proc.absorptionindexdb import AbsorptionIndexDBDef
        from mangadap.proc.absorptionindexdb import AbsorptionIndexDB
        d = AbsorptionIndexDBDef(key='USER',
                                 file_path='/path/to/parameter/file')
        p = AbsorptionIndexDB('USER', indxdb_list=d)

    The reference frame of the absorption-line index can be different
    *for each index*; a column in the SDSS parameter file is used to
    specify either air or vacuum wavelengths.  This is different from,
    e.g., the definition of an emission-line database (see
    :class:`mangadap.proc.emissionlinedb.EmissionLineDB`), which only
    allows you to set the reference frame for the entire database.
    Conversely, :class:`AbsorptionIndexDB` will convert the input index
    definition to vacuum on a case-by-case basis based on the SDSS
    parameter file entries.  This means that the
    :class:`AbsorptionIndexDB` object only provides vacuum wavelengths.
    
*Revision history*:
    | **18 Mar 2016**: Original implementation by K. Westfall (KBW)

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
from ..config.util import validate_absorption_index_config
from ..util.idlutils import airtovac
from ..util.yanny import yanny
from ..par.parset import ParSet
from ..par.bandpassfilter import BandPassFilterPar
from .util import _select_proc_method

__author__ = 'Kyle B. Westfall'
# Add strict versioning
# from distutils.version import StrictVersion

class AbsorptionIndexDBDef(ParSet):
    """
    Class with parameters used to define an absorption-line index
    database.  Options and defaults in ParSet base class are set to
    None.
    """
    def __init__(self, key, file_path):

        pars =     [ 'key', 'file_path' ]
        values =   [   key,   file_path ]
        defaults = [  None,        None ]
        dtypes =   [   str,         str ]

        ParSet.__init__(self, pars, values=values, defaults=defaults, dtypes=dtypes)


def available_absorption_index_databases(dapsrc=None):
    """
    Return the list of database keys and file names for the available
    absorption-line index databases.  The currently available databases
    are:
    
    +-------------+-----+----------------------------------+
    |         KEY |   N | Description                      |
    +=============+=====+==================================+
    |        LICK |  21 | The standard Lick indices from   |
    |             |     | Trager et al. (1998).            |
    +-------------+-----+----------------------------------+
    |    STANDARD |  38 | Includes additional Balmer and   |
    |             |     | IMF-sensitive indices            |
    +-------------+-----+----------------------------------+

    .. warning::
        Function is currently only valid for Python 3.2 or greater!

    Args:
        dapsrc (str): (Optional) Root path to the DAP source
            directory.  If not provided, the default is defined by
            :func:`mangadap.config.defaults.default_dap_source`.

    Returns:
        list: An list of :class:`AbsorptionIndexDBDef` objects, each of
        which defines a unique set of absorption-line indices (see
        :class:`mangadap.par.bandpassfilter.BandPassFilterPar`).

    Raises:
        NotADirectoryError: Raised if the provided or default
            *dapsrc* is not a directory.
        FileNotFoundError: Raised if the database parameter file is not
            found.
        KeyError: Raised if the keywords are not all unique.
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
    ini_files = glob.glob(dapsrc+'/python/mangadap/config/absorption_indices/*.ini')
    if len(ini_files) == 0:
        raise IOError('Could not find any configuration files in {0} !'.format(
                      dapsrc+'/python/mangadap/config/absorption_indices'))

    # Build the list of library definitions
    indx_databases = []
    for f in ini_files:
        # Read the config file
        cnfg = ConfigParser(environ, allow_no_value=True, interpolation=ExtendedInterpolation())
        cnfg.read(f)
        # Ensure it has the necessary elements to define the
        # absorption-line index database
        validate_absorption_index_config(cnfg)
        # Append the definition of the absorption-line index database
        indx_databases += [ AbsorptionIndexDBDef(key=cnfg['default']['key'],
                                                 file_path=cnfg['default']['file_path']) ]

    # Check the keywords of the libraries are all unique
    if len(numpy.unique( numpy.array([indx['key'] for indx in indx_databases]) )) \
            != len(indx_databases):
        raise KeyError('Absorption-line index database keywords are not all unique!')

    # Return the default list of databases
    return indx_databases


class AbsorptionIndexDB:
    """
    Basic container class for the database of absorption-line index
    parameters.

    .. todo::
        - Need to figure out is it's better to have an array of
          BandPassFilterPar objects, or if I should convert self.data to
          a numpy record array.

    Args:
        database_key (str): Keyword selecting the database to use.
        indxdb_list (list): (Optional) List of
            :class:`AbsorptionIndexDBDef` objects that define the
            parameters required to setup the absorption-line index
            database.  The *database_key* must select one of the objects
            in this list.
        dapsrc (str): (Optional) Root path to the DAP source directory.
            If not provided, the default is defined by
            :func:`mangadap.config.defaults.default_dap_source`.

    Attributes:
        version (str): Version number
        database (str): Keyword of the selected database to use.
        dbparset (:class:`AbsorptionIndexDBDef`): Parameter set required
            to setup the absorption-line index database.
        nindx (int): Number of absorption-line indices in the database.
        data (:class:`numpy.array`) : Array of
            :class:`mangadap.par.bandpassfilter.BandPassFilterPar`
            instances, one per absorption-line index.

    """
    def __init__(self, database_key, indxdb_list=None, dapsrc=None):

        self.version = '1.0'

        # Define the library properties
        self.database = None
        self.dbparset = None
        self.nindx = None
        self.data = None

        # Read and process the database
        self.read_database(database_key, indxdb_list=indxdb_list, dapsrc=dapsrc)
        

    def __getitem__(self, index):
        """Select and return a single parameter set.

        Args:
            index (int): Index to return.  Must be in the appropriate
                range.

        Returns:
            :class:`mangadap.par.bandpassfilter.BandPassFilterPar`:
                Parameter set of the selected absorption-line index

        Raises:
            IndexError: Raised if the provided index is out of bounds.
        """
        if index < 0 or index >= self.nindx:
            raise IndexError('Index out of range.')
        return self.data[index]


    def _check(self):
        """
        Check that the database is correctly defined:

            - All the indices must be unique.

        Raises:
            ValueError: Raised if the indices are not all unique.

        .. todo::
            - Other checks needed?

        """
        if len(numpy.unique( numpy.array([d['index'] for d in self.data]))) != self.nindx:
            raise ValueError('Indices in {0} database are not all unique!'.format(self.database))


    def read_database(self, database_key, indxdb_list=None, dapsrc=None):
        """
        Select and read the database from the provided list.  Used to
        set :attr:`database` and :attr:`dbparset`; see
        :func:`mangadap.proc.util._select_proc_method`.

        Args:
            database_key (str): Keyword selecting the database to use.
            indxdb_list (list): (Optional) List of
                :class:`AbsorptionIndexDBDef` objects that define the
                parameters required to setup the absorption-line index
                database.  The *database_key* must select one of the
                objects in this list.
            dapsrc (str): (Optional) Root path to the DAP source
                directory.  If not provided, the default is defined by
                :func:`mangadap.config.defaults.default_dap_source`.
        """
        # Get the details of the selected database
        self.dbparset = _select_proc_method(database_key, AbsorptionIndexDBDef,
                                            indxdb_list=indxdb_list,
                                            available_func=available_absorption_index_databases,
                                            dapsrc=dapsrc)
        self.database = database_key
        
        # Check that the database exists
        if not os.path.isfile(self.dbparset['file_path']):
            raise FileNotFoundError('Database file {0} does not exist!'.format(
                                                                    self.dbparset['file_path']))

        # Read the yanny file
        par = yanny(self.dbparset['file_path'])
        if len(par['DAPABI']['index']) == 0:
            raise ValueError('Could not find DAPABI entries in {0}!'.self.dbparset['file_path'])

        # Setup the array of absorption-line index database parameters
        self.nindx = len(par['DAPABI']['index'])
        self.data = numpy.empty(self.nindx, dtype=object)
        for i in range(self.nindx):
            invac = par['DAPABI']['waveref'][i] == 'vac'
            comp = par['DAPABI']['component'][i] != 0
            self.data[i] = BandPassFilterPar(par['DAPABI']['index'][i],
                                             par['DAPABI']['name'][i],
                                par['DAPABI']['blueside'][i] if invac \
                                        else airtovac(numpy.array(par['DAPABI']['blueside'][i])),
                                par['DAPABI']['redside'][i] if invac \
                                        else airtovac(numpy.array(par['DAPABI']['redside'][i])),
                        primary=par['DAPABI']['primary'][i] if invac \
                                        else airtovac(numpy.array(par['DAPABI']['primary'][i])),
                                             units=par['DAPABI']['units'][i],
                                             component=comp)

        # Check that the database is correctly defined
        self._check()

