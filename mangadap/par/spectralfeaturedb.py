# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""

Container class for databases of spectral features. This implements
the base class used by :class:`mangadap.par.artifactdb.ArtifactDB`,
:class:`mangadap.par.emissionlinedb.EmissionLineDB`,
:class:`mangadap.par.emissionmomentsdb.EmissionMomentsDB`,
:class:`mangadap.par.absorptionindexdb.AbsorptionIndexDB`, and
:class:`mangadap.par.bandheadindexdb.BandheadIndexDB`.

----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
from pathlib import Path

from IPython import embed

import numpy

from pydl.pydlutils.yanny import yanny

from ..config import defaults
from .parset import ParDatabase
from ..util.parser import DefaultConfig
from ..proc import util

# TODO: Enable searching of local directory

class SpectralFeatureDB(ParDatabase):
    """
    Basic container class for the parameters databases of spectral
    features.  This is the base class for all of the following:

        - class used by :class:`~mangadap.par.artifactdb.ArtifactDB`
        - :class:`~mangadap.par.emissionlinedb.EmissionLineDB`
        - :class:`~mangadap.par.emissionmomentsdb.EmissionMomentsDB`
        - :class:`~mangadap.par.absorptionindexdb.AbsorptionIndexDB`
        - :class:`~mangadap.par.bandheadindexdb.BandheadIndexDB`

    See :class:`~mangadap.par.parset.ParDatabase` for additional
    attributes.

    Each derived class must define its own default directory that
    contains the relevant databases, the class that defines the base
    :class:`~mangadap.par.parset.ParSet` for each
    :class:`~mangadap.par.parset.ParDatabase` entry, and the method
    that parses the parameter file into the parameter list.

    The primary instantiation requires the SDSS parameter file. To
    instantiate using a keyword (and optionally a directory that
    holds the parameter files), use :func:`from_key`.

    Args:
        parfile (:obj:`str`, `Path`_):
            The SDSS parameter file with the emission-line database.

    Attributes:
        key (:obj:`str`):
            Database signifying keyword
        file (:obj:`str`):
            File with the emission-line data
        size (:obj:`int`):
            Number of elements in the database. 
    """
    default_data_dir = None
    def __init__(self, parfile):
        # TODO: The approach here (read using yanny, set to par
        # individually, then covert back to record array using
        # ParDatabase) is stupid...
        _parfile = Path(parfile).resolve()
        if not _parfile.exists():
            raise FileNotFoundError(f'{_parfile} does not exist!')

        self.key = util.get_database_key(str(_parfile))
        self.file = str(_parfile)
        self.size = None

        # _parse_yanny() has to set `size` for each subclass.
        ParDatabase.__init__(self, self._parse_yanny())
        self._validate()

    def _validate(self):
        # Ensure that all indices are unique
        if len(numpy.unique(self.data['index'])) != self.size:
            raise ValueError(f'Database indices for {self.key} are not all unique!')

    def _parse_yanny(self):
        raise NotImplementedError(f'_parse_yanny not defined for {self.__class__.__name__}.')

    @classmethod
    def default_path(cls):
        """
        Return the default path with the emission-line databases.
        """
        if cls.default_data_dir is None:
            raise ValueError(f'Default data directory undefined for {cls.__class__.__name__}.')
        return defaults.dap_data_root() / cls.default_data_dir

    @classmethod
    def available_databases(cls, directory_path=None):
        """
        Return the list of available database files.

        Args:
            directory_path (:obj:`str`, optional):
                 Root path with the database files. If None, uses the
                 default directory defined by :func:`default_path`.

        Returns:
            :obj:`dict`: A dictionary with the database files and
            associated keyword.

        Raises:
            NotADirectoryError:
                Raised if the provided or default directory does not
                exist.
            ValueError:
                Raised if the keywords found for all the ``*.par``
                files are not unique.
        """
        if directory_path is None:
            directory_path = cls.default_path()
        if not directory_path.exists():
            raise NotADirectoryError(f'{directory_path} not found!')
        files = list(directory_path.glob('*.par'))
        keys = [util.get_database_key(str(f)) for f in files]
        if len(keys) != len(numpy.unique(keys)):
            raise ValueError(f'Keys read for par files in {directory_path} are not unique!  Names '
                             'of par files must be case-insensitive and unique.')
        return {k: f for k,f in zip(keys, files)}

    @classmethod
    def from_key(cls, key, directory_path=None):
        r"""
        Instantiate the object using a keyword.

        Args:
            key (:obj:`str`):
                Keyword selecting the database to use.
            directory_path (:obj:`str`, optional):
                Root path with the database parameter files. If None,
                uses the default set by :func:`available_databases`.
                Note that the file search includes *any* file with a
                ``.par`` extension. The root of the file should be
                case-insensitive.
        """
        databases = cls.available_databases(directory_path=directory_path)
        available_keys = list(databases.keys())
        if key not in available_keys:
            raise KeyError(f'No database found with key {key}.  Keywords found '
                           f'are: {available_keys}')
        return cls(databases[key])


