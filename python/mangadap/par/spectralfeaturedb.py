# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""

Container class for databases of spectral features.  This is the base
class used by :class:`mangadap.par.artifactdb.ArtifactDB`,
:class:`mangadap.par.emissionlinedb.EmissionLineDB`,
:class:`mangadap.par.emissionmomentsdb.EmissionMomentsDB`,
:class:`mangadap.par.absorptionindexdb.AbsorptionIndexDB`, and
:class:`mangadap.par.bandheadindexdb.BandheadIndexDB`.

*License*:
    Copyright (c) 2015, SDSS-IV/MaNGA Pipeline Group
        Licensed under BSD 3-clause license - see LICENSE.rst

*Source location*:
    $MANGADAP_DIR/python/mangadap/par/spectralfeaturedb.py

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

        from .parset import ParSet
        from ..config.defaults import default_dap_source

.. warning::

    Because of the use of the ``ExtendedInterpolation`` in
    `configparser.ConfigParser`_,
    :func:`available_absorption_index_databases` is not python 2
    compiliant.

*Class usage examples*:
    Add usage examples
    
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

from .parset import ParSet
from ..config.defaults import default_dap_source

__author__ = 'Kyle B. Westfall'
# Add strict versioning
# from distutils.version import StrictVersion

class SpectralFeatureDBDef(ParSet):
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


def validate_spectral_feature_config(cnfg):
    """ 
    Validate the `configparser.ConfigParser`_ object that is meant to
    define a spectral-feature database.

    Args:
        cnfg (`configparser.ConfigParser`_): Object meant to contain
            defining parameters of a spectral-feature database.

    Raises:
        KeyError: Raised if required keyword does not exist.
        ValueError: Raised if key has unacceptable value.

    """
    # Check for required keywords
    if 'key' not in cnfg.options('default'):
        raise KeyError('Keyword \'key\' must be provided.')
    if 'file_path' not in cnfg.options('default'):
        raise KeyError('Keyword \'file_path\' must be provided.')


def available_spectral_feature_databases(sub_directory, dapsrc=None):
    """
    Generic routine for finding available spectral-feature database
    definitions in the config directory path.

    .. warning::
        Function is currently only valid for Python 3.2 or greater!

    Args:
        dapsrc (str): (Optional) Root path to the DAP source
            directory.  If not provided, the default is defined by
            :func:`mangadap.config.defaults.default_dap_source`.

    Returns:
        list: An list of :class:`SpectralFeatureDBDef` objects.

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
    """
    # Check the source directory exists
    dapsrc = default_dap_source() if dapsrc is None else str(dapsrc)
    if not os.path.isdir(dapsrc):
        raise NotADirectoryError('{0} does not exist!'.format(dapsrc))

    # Check the configuration files exist
    ini_path = os.path.join(dapsrc, 'python/mangadap/config', sub_directory, '*.ini')
    ini_files = glob.glob(ini_path)
    if len(ini_files) == 0:
        raise IOError('Could not find any configuration files in {0} !'.format(ini_path))

    # Build the list of library definitions
    databases = []
    for f in ini_files:
        # Read the config file
        cnfg = ConfigParser(environ, allow_no_value=True, interpolation=ExtendedInterpolation())
        cnfg.read(f)
        # Ensure it has the necessary elements to define the
        # absorption-line index database
        validate_spectral_feature_config(cnfg)
        # Append the definition of the absorption-line index database
        databases += [ SpectralFeatureDBDef(key=cnfg['default']['key'],
                                            file_path=cnfg['default']['file_path']) ]

    # Check the keywords of the libraries are all unique
    if len(numpy.unique( numpy.array([db['key'] for db in databases]) )) != len(databases):
        raise KeyError('Spectral-feature database keywords are not all unique!')

    # Return the default list of databases
    return databases


