# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Container class for the database of bandhead indices to measure.

*License*:
    Copyright (c) 2015, SDSS-IV/MaNGA Pipeline Group
        Licensed under BSD 3-clause license - see LICENSE.rst

*Source location*:
    $MANGADAP_DIR/python/mangadap/proc/bandheadindexdb.py

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
        from ..config.util import validate_bandhead_index_config
        from ..util.idlutils import airtovac
        from ..util.yanny import yanny
        from ..par.parset import ParSet
        from ..par.bandpassfilter import BandPassFilterPar
        from .util import _select_proc_method

.. warning::

    Because of the use of the ``ExtendedInterpolation`` in
    `configparser.ConfigParser`_,
    :func:`available_bandhead_index_databases` is not python 2
    compiliant.

*Class usage examples*:
    Bandhead index databases are defined using SDSS parameter files.  To
    define a database, you can use one of the default set of available
    bandhead index databases (see
    :func:`available_bandhead_index_databases`)::
    
        from mangadap.proc.bandheadindexdb import BandheadIndexDB
        p = BandheadIndexDB('STANDARD')

    The above call requires that the ``$MANGADAP_DIR`` environmental
    variable is set.  If it is not, you can define it's location, as
    in::

        from mangadap.proc.bandheadindexdb import BandheadIndexDB
        p = BandheadIndexDB('STANDARD', dapsrc='/path/to/dap/source')

    Finally, you can create your own SDSS parameter file (see
    :class:`mangadap.util.yanny`) with your own bandhead indices
    to use.  Example files are provided in
    ``$MANGADAP_DIR/external/bandhead_indices`` with a companion
    ``README`` file.  With your own file, you have to point to the file
    using :class:`BandheadIndexDBDef`, which you can then pass to
    :class:`BandheadIndexDB`::

        from mangadap.proc.bandheadindexdb import BandheadIndexDBDef
        from mangadap.proc.bandheadindexdb import BandheadIndexDB
        d = BandheadIndexDBDef(key='USER',
                               file_path='/path/to/parameter/file')
        p = BandheadIndexDB('USER', indxdb_list=d)

    The reference frame of the bandhead index can be different *for each
    index*; a column in the SDSS parameter file is used to specify
    either air or vacuum wavelengths.  This is different from, e.g., the
    definition of an emission-line database (see
    :class:`mangadap.proc.emissionlinedb.EmissionLineDB`), which only
    allows you to set the reference frame for the entire database.
    Conversely, :class:`BandheadIndexDB` will convert the input index
    definition to vacuum on a case-by-case basis based on the SDSS
    parameter file entries.  This means that the
    :class:`BandheadIndexDB` object only provides vacuum wavelengths.
    
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
from ..config.util import validate_bandhead_index_config
from ..util.idlutils import airtovac
from ..util.yanny import yanny
from ..par.parset import ParSet
from ..par.bandpassfilter import BandPassFilterPar
from .util import _select_proc_method

__author__ = 'Kyle B. Westfall'
# Add strict versioning
# from distutils.version import StrictVersion

class BandheadIndexDBDef(ParSet):
    """
    Class with parameters used to define an bandhead index database.
    Options and defaults in ParSet base class are set to None.
    """
    def __init__(self, key, file_path):

        pars =     [ 'key', 'file_path' ]
        values =   [   key,   file_path ]
        defaults = [  None,        None ]
        dtypes =   [   str,         str ]

        ParSet.__init__(self, pars, values=values, defaults=defaults, dtypes=dtypes)


def available_bandhead_index_databases(dapsrc=None):
    """
    Return the list of database keys and file names for the available
    bandhead-line index databases.  The currently available databases
    are:
    
    +-------------+-----+----------------------------------+
    |         KEY |   N | Description                      |
    +=============+=====+==================================+
    |    STANDARD |   3 | D4000, Dn4000, and Ti0 (8860 A). |
    +-------------+-----+----------------------------------+

    .. warning::
        Function is currently only valid for Python 3.2 or greater!

    Args:
        dapsrc (str): (Optional) Root path to the DAP source
            directory.  If not provided, the default is defined by
            :func:`mangadap.config.defaults.default_dap_source`.

    Returns:
        list: An list of :class:`BandheadIndexDBDef` objects, each of
        which defines a unique set of bandhead indices (see
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
    ini_files = glob.glob(dapsrc+'/python/mangadap/config/bandhead_indices/*.ini')
    if len(ini_files) == 0:
        raise IOError('Could not find any configuration files in {0} !'.format(
                      dapsrc+'/python/mangadap/config/bandhead_indices'))

    # Build the list of library definitions
    indx_databases = []
    for f in ini_files:
        # Read the config file
        cnfg = ConfigParser(environ, allow_no_value=True, interpolation=ExtendedInterpolation())
        cnfg.read(f)

        # Ensure it has the necessary elements to define the bandhead
        # index database
        validate_bandhead_index_config(cnfg)
        # Append the definition of the bandhead index database
        indx_databases += [ BandheadIndexDBDef(key=cnfg['default']['key'],
                                               file_path=cnfg['default']['file_path']) ]

    # Check the keywords of the libraries are all unique
    if len(numpy.unique( numpy.array([indx['key'] for indx in indx_databases]) )) \
            != len(indx_databases):
        raise KeyError('Bandhead index database keywords are not all unique!')

    # Return the default list of databases
    return indx_databases


class BandheadIndexDB:
    """
    Basic container class for the database of bandhead index parameters.

    .. todo::
        - Need to figure out is it's better to have an array of
          BandPassFilterPar objects, or if I should convert self.data to
          a numpy record array.

    Args:
        database_key (str): Keyword selecting the database to use.
        indxdb_list (list): (Optional) List of
            :class:`BandheadIndexDBDef` objects that define the
            parameters required to setup the bandhead index database.
            The *database_key* must select one of the objects in this
            list.
        dapsrc (str): (Optional) Root path to the DAP source directory.
            If not provided, the default is defined by
            :func:`mangadap.config.defaults.default_dap_source`.

    Attributes:
        version (str): Version number
        database (str): Keyword of the selected database to use.
        dbparset (:class:`BandheadIndexDBDef`): Parameter set required
            to setup the bandhead index database.
        nindx (int): Number of bandhead indices in the database.
        data (:class:`numpy.array`) : Array of
            :class:`mangadap.par.bandpassfilter.BandPassFilterPar`
            instances, one per bandhead index.

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
            Parameter set of the selected bandhead index

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

        .. todo::
            - Other checks needed?

        Raises:
            ValueError: Raised if the indices are not all unique.

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
                :class:`BandheadIndexDBDef` objects that define the
                parameters required to setup the bandhead index
                database.  The *database_key* must select one of the
                objects in this list.
            dapsrc (str): (Optional) Root path to the DAP source
                directory.  If not provided, the default is defined by
                :func:`mangadap.config.defaults.default_dap_source`.

        Raises:
            KeyError: Raised if the selected keyword is not among the
                provided list or if the provided list has more than one
                identical keyword.
            TypeError: Raised if the input *indxdb_list* object is not a
                list or a :class:`BandheadIndexDBDef`.

        """
        # Get the details of the selected database
        self.dbparset = _select_proc_method(database_key, BandheadIndexDBDef,
                                            indxdb_list=indxdb_list,
                                            available_func=available_bandhead_index_databases,
                                            dapsrc=dapsrc)
        self.database = database_key
        
        # Check that the database exists
        if not os.path.isfile(self.dbparset['file_path']):
            raise FileNotFoundError('Database file {0} does not exist!'.format(
                                                                    self.dbparset['file_path']))

        # Read the yanny file
        par = yanny(self.dbparset['file_path'])
        if len(par['DAPBHI']['index']) == 0:
            raise ValueError('Could not find DAPBHI entries in {0}!'.self.dbparset['file_path'])

        # Setup the array of bandhead index database parameters
        self.nindx = len(par['DAPBHI']['index'])
        self.data = numpy.empty(self.nindx, dtype=object)
        for i in range(self.nindx):
            invac = par['DAPBHI']['waveref'][i] == 'vac'
            self.data[i] = BandPassFilterPar(par['DAPBHI']['index'][i],
                                             par['DAPBHI']['name'][i],
                                par['DAPBHI']['blueside'][i] if invac \
                                        else airtovac(numpy.array(par['DAPBHI']['blueside'][i])),
                                par['DAPBHI']['redside'][i] if invac \
                                        else airtovac(numpy.array(par['DAPBHI']['redside'][i])),
                                             integrand=par['DAPBHI']['integrand'][i],
                                             order=par['DAPBHI']['order'][i])

        # Check that the database is correctly defined
        self._check()

