# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""

Container class for databases of spectral features.  This is the base
class used by :class:`mangadap.par.artifactdb.ArtifactDB`,
:class:`mangadap.par.emissionlinedb.EmissionLineDB`,
:class:`mangadap.par.emissionmomentsdb.EmissionMomentsDB`,
:class:`mangadap.par.absorptionindexdb.AbsorptionIndexDB`, and
:class:`mangadap.par.bandheadindexdb.BandheadIndexDB`.

Revision history
----------------

    | **18 Mar 2016**: Original implementation by K. Westfall (KBW)
    | **25 Feb 2017**: (KBW) Change to using
        :class:`mangadap.util.parser.DefaultConfig`

----

.. include license and copyright
.. include:: ../copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../links.rst
"""

import os
import glob
import numpy

from .parset import ParSet
from ..config.defaults import dap_source_dir
from ..util.parser import DefaultConfig

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
    Validate the :class:`mangadap.util.parser.DefaultConfig` object with
    the spectral-feature database parameters.

    Args:
        cnfg (:class:`mangadap.util.parser.DefaultConfig`): Object with
            spectral-feature database parameters.

    Raises:
        KeyError: Raised if required keyword does not exist.

    """
    # Check for required keywords
    for k in [ 'key', 'file_path']:
        if k not in cnfg:
            raise KeyError('No keyword \'{0}\' in spectral feature config!'.format(k))


def available_spectral_feature_databases(sub_directory, dapsrc=None):
    """
    Generic routine for finding available spectral-feature database
    definitions in the config directory path.

    .. warning::
        Function is currently only valid for Python 3.2 or greater!

    Args:
        dapsrc (str): (Optional) Root path to the DAP source
            directory.  If not provided, the default is defined by
            :func:`mangadap.config.defaults.dap_source_dir`.

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
    dapsrc = dap_source_dir() if dapsrc is None else str(dapsrc)
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
        # Read and validate the config file
        cnfg = DefaultConfig(f=f, interpolate=True)
        validate_spectral_feature_config(cnfg)
        # Append the definition of the absorption-line index database
        databases += [ SpectralFeatureDBDef(key=cnfg['key'], file_path=cnfg['file_path']) ]

    # Check the keywords of the libraries are all unique
    if len(numpy.unique( numpy.array([db['key'] for db in databases]) )) != len(databases):
        raise KeyError('Spectral-feature database keywords are not all unique!')

    # Return the default list of databases
    return databases


