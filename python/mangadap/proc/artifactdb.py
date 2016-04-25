# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""

Container class for a database with a list of spectral regions with
known artifacts that should be ignored during analysis of the data.
These can be anything, but is only currently used to define spectral
regions with poorly subtracted sky lines.  They are also currently
independent of spatial position and expected to be applied for all
spectra in an RSS or CUBE file.

*License*:
    Copyright (c) 2016, SDSS-IV/MaNGA Pipeline Group
        Licensed under BSD 3-clause license - see LICENSE.rst

*Source location*:
    $MANGADAP_DIR/python/mangadap/proc/artifactdb.py

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
        import numpy

        from ..util.idlutils import airtovac
        from ..util.yanny import yanny
        from ..par.parset import ParSet, ParDatabase
        from .spectralfeaturedb import available_spectral_feature_databases, SpectralFeatureDBDef
        from .util import _select_proc_method

*Class usage examples*:
    Add example usage!

*Revision history*:
    | **16 Apr 2016**: Original implementation by K. Westfall (KBW)
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
import numpy

from ..util.idlutils import airtovac
from ..util.yanny import yanny
from ..par.parset import ParSet, ParDatabase
from .spectralfeaturedb import available_spectral_feature_databases, SpectralFeatureDBDef
from .util import _select_proc_method

__author__ = 'Kyle B. Westfall'
# Add strict versioning
# from distutils.version import StrictVersion


class ArtifactPar(ParSet):
    """
    Parameter object that defines a set of artifacts to be ignored
    during analysis.
    
    See :class:`mangadap.par.parset.ParSet` for attributes and raised
    exceptions.

    Args:
        index (int) : An index used to refer to the line in the *line*
            and *mode* attributes.
        name (str) : A name for the line.
        waverange (numpy.ndarray, list) : A two-element vector with the
            starting and ending wavelength (angstroms in VACUUM) where
            the artifact affects the data.
    """
    def __init__(self, index, name, waverange):
        
        arr_like = [ numpy.ndarray, list ]

        _name = name.strip()

        pars =     [ 'index', 'name', 'waverange' ]
        values =   [   index,  _name,   waverange ]
        dtypes =   [     int,    str,    arr_like ]

        ParSet.__init__(self, pars, values=values, dtypes=dtypes)


    def _check(self):
        """
        Check the parameter list:
            
            - Make sure the waverange only has two elements.

        Raises:
            ValueError: Raised if one of the conditions above are not
                met.
        """
        if len(self.data['waverange']) != 2:
            raise ValueError('Wavelength range must have two and only two elements.')


def available_artifact_databases(dapsrc=None):
    """
    Return the list of database keys and file names for the available
    artifact databases.  The currently available databases are:
    
    +-------------+-----+----------------------------------+
    |         KEY |   N | Description                      |
    +=============+=====+==================================+
    |         SKY |   1 | Poorly subtracted sky lines      |
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
    return available_spectral_feature_databases('artifacts', dapsrc=dapsrc)


class ArtifactDB(ParDatabase):
    """
    Basic container class for the database of artifact parameters.  See
    :class:`mangadap.parset.ParDatabase` for additional attributes.

    Args:
        database_key (str): Keyword selecting the database to use.
        artdb_list (list): (**Optional**) List of
            :class:`SpectralFeatureDBDef` objects that defines the
            unique key for the database and the path to the source SDSS
            parameter file.
        dapsrc (str): (**Optional**) Root path to the DAP source
            directory.  If not provided, the default is defined by
            :func:`mangadap.config.defaults.default_dap_source`.

    Attributes:
        version (str): Version number
        database (:class:`mangadap.par.ParSet`): Database parameters.
        nart (int): Number of artifacts in the database

    """
    def __init__(self, database_key, artdb_list=None, dapsrc=None):

        # TODO: The approach here (read using yanny, set to par
        # individually, then covert back to record array using
        # ParDatabase) is stupid...

        self.version = '1.0'

        # Get the details of the selected database
        self.database = _select_proc_method(database_key, SpectralFeatureDBDef,
                                            method_list=artdb_list,
                                            available_func=available_artifact_databases,
                                            dapsrc=dapsrc)
        
        # Check that the database exists
        if not os.path.isfile(self.database['file_path']):
            raise FileNotFoundError('Database file {0} does not exist!'.format(
                                                                    self.database['file_path']))

        # Read the yanny file
        par = yanny(self.database['file_path'])
        if len(par['DAPART']['index']) == 0:
            raise ValueError('Could not find DAPART entries in {0}!'.self.database['file_path'])

        # Setup the array of emission line database parameters
        self.nart = len(par['DAPART']['index'])
        parlist = []
        for i in range(self.nart):
            invac = par['DAPART']['waveref'][i] == 'vac'
            parlist += [ ArtifactPar(par['DAPART']['index'][i], par['DAPART']['name'][i],
                                     numpy.asarray(par['DAPART']['waverange'][i]) \
                                      if invac else airtovac(par['DAPEML']['waverange'][i]) )]

        ParDatabase.__init__(self, parlist)

        # Ensure that all indices are unique
        if len(numpy.unique(self.data['index'])) != self.nart:
            raise ValueError('Indices in {0} database are not all unique!'.format(
                                                                            self.database['key']))
        


    


