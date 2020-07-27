# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Container class for a database with a list of spectral regions with
known artifacts that should be ignored during analysis of the data.
These can be anything, but is only currently used to define spectral
regions with poorly subtracted sky lines.  They are also currently
independent of spatial position and expected to be applied for all
spectra in an RSS or CUBE file.

----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

import os
import numpy

from pydl.goddard.astro import airtovac
from pydl.pydlutils.yanny import yanny
from .parset import KeywordParSet, ParDatabase
from .spectralfeaturedb import SpectralFeatureDB

# Add strict versioning
# from distutils.version import StrictVersion

class ArtifactPar(KeywordParSet):
    """
    Parameter object that defines a set of artifacts to be ignored
    during analysis.

    The defined parameters are:

    .. include:: ../tables/artifactpar.rst
    
    See :class:`mangadap.par.parset.ParSet` for attributes and raised
    exceptions.
    """
    def __init__(self, index=None, name=None, waverange=None):
        
        arr_like = [ numpy.ndarray, list ]

        _name = None if name is None else name.strip()

        pars =     [ 'index', 'name', 'waverange' ]
        values =   [   index,  _name,   waverange ]
        dtypes =   [     int,    str,    arr_like ]

        descr = ['An index used to refer to the artifact',
                 'A name for the artifact',
                 'A two-element vector with the starting and ending wavelength (angstroms ' \
                 'in **vacuum**) where the artifact affects the data.']

        super(ArtifactPar, self).__init__(pars, values=values, dtypes=dtypes, descr=descr)

    def _check(self):
        """
        Check the parameter list.

        Only check performed is to make sure the wavelength range
        only has two elements.

        Raises:
            ValueError:
                Raised if one of the conditions above are not met.
        """
        if len(self.data['waverange']) != 2:
            raise ValueError('Wavelength range must have two and only two elements.')


class ArtifactDB(SpectralFeatureDB):
    """
    Basic container class for the database of artifacts.

    Each row of the database is parsed using :class:`ArtifactPar`.

    The primary instantiation requires the SDSS parameter file with
    the artifact data. To instantiate using a keyword (and
    optionally a directory that holds the parameter files), use the
    :func:`mangadap.par.spectralfeaturedb.SpectralFeatureDB.from_key`
    class method.  See the base class for additional attributes.

    Args:
        parfile (:obj:`str`):
            The SDSS parameter file with the artifact database.

    Attributes:
        key (:obj:`str`):
            Database signifying keyword
        file (:obj:`str`):
            File with the artifact data
        size (:obj:`int`):
            Number of artifacts in the database. 
    """
    default_data_dir = 'artifacts'
    def _parse_yanny(self):
        """
        Parse the yanny file (provided by :attr:`file`) for the artifact
        database.

        Returns:
            :obj:`list`: The list of
            :class:`mangadap.par.parset.ParSet` instances for each
            line of the database.
        """
        # Read the yanny file
        par = yanny(filename=self.file, raw=True)
        if len(par['DAPART']['index']) == 0:
            raise ValueError('Could not find DAPART entries in {0}!'.format(self.file))

        # Setup the array of emission line database parameters
        self.size = len(par['DAPART']['index'])
        parlist = []
        for i in range(self.size):
            invac = par['DAPART']['waveref'][i] == 'vac'
            parlist += [ ArtifactPar(index=par['DAPART']['index'][i], name=par['DAPART']['name'][i],
                                     waverange=numpy.asarray(par['DAPART']['waverange'][i]) \
                                      if invac else airtovac(par['DAPEML']['waverange'][i]) )]
        return parlist
