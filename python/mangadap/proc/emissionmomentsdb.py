# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Container class for the database of emission-line bandpass filters used
for non-parameteric measurements of emission-line flux and velocity
moments.

*License*:
    Copyright (c) 2015, SDSS-IV/MaNGA Pipeline Group
        Licensed under BSD 3-clause license - see LICENSE.rst

*Source location*:
    $MANGADAP_DIR/python/mangadap/proc/emissionmomentsdb.py

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
        from ..config.util import validate_emission_bandpass_filter_config
        from ..util.idlutils import airtovac
        from ..util.yanny import yanny
        from ..par.parset import ParSet
        from ..par.bandpassfilter import BandPassFilterPar
        from .util import _select_proc_method

.. warning::

    Because of the use of the ``ExtendedInterpolation`` in
    `configparser.ConfigParser`_,
    :func:`available_emission_bandpass_filter_databases` is not python 2
    compiliant.

*Class usage examples*:
    Emission-line moment databases are defined using SDSS parameter
    files.  To define a database, you can use one of the default set of
    available emission-line moment databases (see
    :func:`available_emission_bandpass_filter_databases`)::
    
        from mangadap.proc.emissionmomentsdb import EmissionMomentsDB
        p = EmissionMomentsDB('STRONG')

    The above call requires that the ``$MANGADAP_DIR`` environmental
    variable is set.  If it is not, you can define it's location, as
    in::

        from mangadap.proc.emissionmomentsdb import EmissionMomentsDB
        p = EmissionMomentsDB('STRONG', dapsrc='/path/to/dap/source')

    Finally, you can create your own SDSS parameter file (see
    :class:`mangadap.util.yanny`) with your own emission line passbands
    to use.  Example files are provided in
    ``$MANGADAP_DIR/external/emission_bandpass_filters`` with a
    companion ``README`` file.  With your own file, you have to point to
    the file using :class:`EmissionMomentsDBDef`, which you can then
    pass to :class:`EmissionMomentsDB`::

        from mangadap.proc.emissionmomentsdb import EmissionMomentsDBDef
        from mangadap.proc.emissionmomentsdb import EmissionMomentsDB
        d = EmissionMomentsDBDef(key='USER',
                                 file_path='/path/to/parameter/file',
                                 in_vacuum=True)
        p = EmissionMomentsDB('USER', emldb_list=d)

    The reference frame of the emission-line passband wavelengths must
    be defined as either vacuum or air, using 'in_vacuum'.  It is
    expected that the object spectra to be analyzed are calibrated to
    vacuum wavelengths.  If 'in_vacuum' is false, this class will use
    :func:`mangadap.util.idlutils.airtovac` to convert the emission-line
    bandpass-filter wavelengths to vacuum.

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
from ..config.util import validate_emission_bandpass_filter_config
from ..util.idlutils import airtovac
from ..util.yanny import yanny
from ..par.parset import ParSet
from ..par.bandpassfilter import BandPassFilterPar
from .util import _select_proc_method

__author__ = 'Kyle B. Westfall'
# Add strict versioning
# from distutils.version import StrictVersion

class EmissionMomentsDBDef(ParSet):
    """
    Class with parameters used to define an emission-line moment
    database.  Options and defaults in ParSet base class are set to
    None.
    """
    def __init__(self, key, file_path, in_vacuum=None ):

        pars =     [ 'key', 'file_path', 'in_vacuum' ]
        values =   [   key,   file_path,   in_vacuum ]
        defaults = [  None,        None,        True ]
        dtypes =   [   str,         str,        bool ]

        ParSet.__init__(self, pars, values=values, defaults=defaults, dtypes=dtypes)


def available_emission_bandpass_filter_databases(dapsrc=None):
    """
    Return the list of database keys and file names for the available
    emission-line moment databases.  The currently available databases
    are:
    
    +-------------+-----+----------------------------------+
    |         KEY |   N | Description                      |
    +=============+=====+==================================+
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
        list: An list of :class:`EmissionMomentsDBDef` objects, each of
        which defines a unique set of emission-line bandpass filters
        (see :class:`mangadap.par.bandpassfilter.BandPassFilterPar`),
        which are used to construct the flux and velocity moments for a
        set of emission lines.

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
    ini_files = glob.glob(dapsrc+'/python/mangadap/config/emission_bandpass_filters/*.ini')
    if len(ini_files) == 0:
        raise IOError('Could not find any configuration files in {0} !'.format(
                      dapsrc+'/python/mangadap/config/emission_bandpass_filters'))

    # Build the list of library definitions
    eml_databases = []
    for f in ini_files:
        # Read the config file
        cnfg = ConfigParser(environ, allow_no_value=True, interpolation=ExtendedInterpolation())
        cnfg.read(f)
        # Ensure it has the necessary elements to define the bandpass
        # filter database
        validate_emission_bandpass_filter_config(cnfg)
        # Append the definition of the emission-line bandpass filter database
        eml_databases += [ EmissionMomentsDBDef(key=cnfg['default']['key'],
                                                file_path=cnfg['default']['file_path'],
                                                in_vacuum=cnfg['default'].getboolean('in_vacuum'))]

    # Check the keywords of the libraries are all unique
    if len(numpy.unique( numpy.array([eml['key'] for eml in eml_databases]) )) \
            != len(eml_databases):
        raise KeyError('Emission-line bandpass filter database keywords are not all unique!')

    # Return the default list of databases
    return eml_databases


class EmissionMomentsDB:
    """
    Basic container class for the database of emission-line bandpass
    filter parameters.

    .. todo::
        - Need to figure out is it's better to have an array of
          BandPassFilterPar objects, or if I should convert self.data to
          a numpy record array.

    Args:
        database_key (str): Keyword selecting the database to use.
        emldb_list (list): (Optional) List of
            :class:`EmissionMomentsDBDef` objects that define the
            parameters required to setup the emission-line bandpass
            filter database.  The *database_key* must select one of the
            objects in this list.
        dapsrc (str): (Optional) Root path to the DAP source directory.
            If not provided, the default is defined by
            :func:`mangadap.config.defaults.default_dap_source`.

    Attributes:
        version (str): Version number
        database (str): Keyword of the selected database to use.
        dbparset (:class:`EmissionMomentsDBDef`): Parameter set required
            to setup the emission-line bandpass filter database.
        neml (int): Number of emission-line bandpass filters in the
            database
        data (:class:`numpy.array`) : Array of
            :class:`mangadap.par.bandpassfilter.BandPassFilterPar`
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
            :class:`mangadap.par.bandpassfilter.BandPassFilterPar`:
            Parameter set of the selected emission-line bandpass filter.

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
                :class:`EmissionMomentsDBDef` objects that define the
                parameters required to setup the emission-line bandpass
                filter database.  The *database_key* must select one of
                the objects in this list.
            dapsrc (str): (Optional) Root path to the DAP source
                directory.  If not provided, the default is defined by
                :func:`mangadap.config.defaults.default_dap_source`.
        """
        # Get the details of the selected database
        self.dbparset = _select_proc_method(database_key, EmissionMomentsDBDef,
                                            emldb_list=emldb_list,
                                        available_func=available_emission_bandpass_filter_databases,
                                            dapsrc=dapsrc)
        self.database = database_key
        
        # Check that the database exists
        if not os.path.isfile(self.dbparset['file_path']):
            raise FileNotFoundError('Database file {0} does not exist!'.format(
                                                                    self.dbparset['file_path']))

        # Read the yanny file
        par = yanny(self.dbparset['file_path'])
        if len(par['DAPELB']['index']) == 0:
            raise ValueError('Could not find DAPELB entries in {0}!'.self.dbparset['file_path'])

        # Setup the array of emission-line bandpass-filter database parameters
        self.neml = len(par['DAPELB']['index'])
        self.data = numpy.empty(self.neml, dtype=object)
        for i in range(self.neml):
            self.data[i] = BandPassFilterPar(par['DAPELB']['index'][i],
                                             par['DAPELB']['name'][i],
                                par['DAPELB']['blueside'][i] if self.dbparset['in_vacuum'] \
                                        else airtovac(numpy.array(par['DAPELB']['blueside'][i])),
                                par['DAPELB']['redside'][i] if self.dbparset['in_vacuum'] \
                                        else airtovac(numpy.array(par['DAPELB']['redside'][i])),
                        restwave=par['DAPELB']['lambda'][i] if self.dbparset['in_vacuum'] \
                                        else airtovac(par['DAPELB']['lambda'][i]),
                        primary=par['DAPELB']['primary'][i] if self.dbparset['in_vacuum'] \
                                        else airtovac(numpy.array(par['DAPELB']['primary'][i])))

        # Check that the database is correctly defined
        self._check()

