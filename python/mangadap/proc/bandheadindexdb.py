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

        import os.path
        import numpy

        from pydl.goddard.astro import airtovac
        from pydl.pydlutils.yanny import yanny
        from ..par.parset import ParDatabase
        from ..proc.bandpassfilter import BandPassFilterPar
        from .spectralfeaturedb import available_spectral_feature_databases, SpectralFeatureDBDef
        from .util import _select_proc_method


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

    Finally, you can create your own `SDSS-style parameter file`_ with
    your own bandhead indices to use.  Example files are provided in
    ``$MANGADAP_DIR/data/bandhead_indices`` with a companion
    ``README`` file.  With your own file, you have to point to the file
    using :class:`SpectralFeatureDBDef`, which you can then pass to
    :class:`BandheadIndexDB`::

        from mangadap.proc.spectralfeaturedb import SpectralFeatureDBDef
        from mangadap.proc.bandheadindexdb import BandheadIndexDB
        d = SpectralFeatureDBDef(key='USER',
                                 file_path='/path/to/parameter/file')
        p = BandheadIndexDB('USER', indxdb_list=d)

    The reference frame of the bandhead index can be different *for each
    index*; a column in the SDSS parameter file is used to specify
    either air or vacuum wavelengths.  When reading the input parameter
    file, :class:`BandheadIndexDB` will convert the input index
    definition to vacuum on a case-by-case basis based on the SDSS
    parameter file entries, meaning :class:`BandheadIndexDB` only
    provides vacuum wavelengths.
    
*Revision history*:
    | **18 Mar 2016**: Original implementation by K. Westfall (KBW)
    | **11 May 2016**: (KBW) Switch to using `pydl.pydlutils.yanny`_ and
        `pydl.goddard.astro.airtovac`_ instead of internal functions

.. _pydl.pydlutils.yanny: http://pydl.readthedocs.io/en/stable/api/pydl.pydlutils.yanny.yanny.html
.. _pydl.goddard.astro.airtovac: http://pydl.readthedocs.io/en/stable/api/pydl.goddard.astro.airtovac.html#pydl.goddard.astro.airtovac
.. _SDSS-style parameter file: http://www.sdss.org/dr12/software/par/

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

#from ..util.idlutils import airtovac
#from ..util.yanny import yanny
from pydl.goddard.astro import airtovac
from pydl.pydlutils.yanny import yanny
from ..par.parset import ParDatabase
from ..proc.bandpassfilter import BandPassFilterPar
from .spectralfeaturedb import available_spectral_feature_databases, SpectralFeatureDBDef
from .util import _select_proc_method

__author__ = 'Kyle B. Westfall'
# Add strict versioning
# from distutils.version import StrictVersion

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

    This is a simple wrapper for
    :func:`mangadap.proc.available_spectral_feature_databases`.

    Args:
        dapsrc (str): (**Optional**) Root path to the DAP source
            directory.  If not provided, the default is defined by
            :func:`mangadap.config.defaults.default_dap_source`.

    Returns:
        list: An list of :class:`SpectralFeatureDBDef` objects, each of
        which defines a unique set of bandhead indices (see
        :class:`mangadap.proc.bandpassfilter.BandPassFilterPar`).

    .. todo::
        - Add backup function for Python 2.
        - Somehow add a python call that reads the databases and
          constructs the table for presentation in sphinx so that the
          text above doesn't have to be edited with changes in the
          available databases.

    """
    return available_spectral_feature_databases('bandhead_indices', dapsrc=dapsrc)
    

class BandheadIndexDB(ParDatabase):
    """
    Basic container class for the database of bandhead-index parameters.
    See :class:`mangadap.parset.ParDatabase` for additional attributes.

    Args:
        database_key (str): Keyword selecting the database to use.
        indxdb_list (list): (**Optional**) List of
            :class:`SpectralFeatureDBDef` objects that define the
            parameters required to setup the bandhead-index database.
            The *database_key* must select one of the objects in this
            list.
        dapsrc (str): (**Optional**) Root path to the DAP source
            directory.  If not provided, the default is defined by
            :func:`mangadap.config.defaults.default_dap_source`.

    Attributes:
        version (str): Version number
        database (:class:`mangadap.par.ParSet`): Database parameters.
        nindx (int): Number of artifacts in the database

    """
    def __init__(self, database_key, indxdb_list=None, dapsrc=None):


        # TODO: The approach here (read using yanny, set to par
        # individually, then covert back to record array using
        # ParDatabase) is stupid...

        self.version = '1.0'

        # Get the details of the selected database
        self.database = _select_proc_method(database_key, SpectralFeatureDBDef,
                                            method_list=indxdb_list,
                                            available_func=available_bandhead_index_databases,
                                            dapsrc=dapsrc)
        
        # Check that the database exists
        if not os.path.isfile(self.database['file_path']):
            raise FileNotFoundError('Database file {0} does not exist!'.format(
                                                                    self.database['file_path']))

        # Read the yanny file
#        par = yanny(self.database['file_path'])
        par = yanny(filename=self.database['file_path'], raw=True)
        if len(par['DAPBHI']['index']) == 0:
            raise ValueError('Could not find DAPBHI entries in {0}!'.format(
                                                                    self.database['file_path']))

        # Check if any of the bands are dummy bands and warn the user
        self.dummy = numpy.any(numpy.array(par['DAPBHI']['blueside']) < 0, axis=1)
        self.dummy |= numpy.any(numpy.array(par['DAPBHI']['redside']) < 0, axis=1)
        if numpy.sum(self.dummy) > 0:
            warnings.warn('Bands with negative wavelengths are used to insert dummy values.'
                          '  Ignoring input bands with indices: {0}'.format(
                                                numpy.array(par['DAPBHI']['index'])[self.dummy]))

        # Setup the array of bandhead index database parameters
        self.nindx = len(par['DAPBHI']['index'])
        parlist = []
        for i in range(self.nindx):
            invac = par['DAPBHI']['waveref'][i] == 'vac'
            parlist += [ BandPassFilterPar(par['DAPBHI']['index'][i],
                                           par['DAPBHI']['name'][i],
                                par['DAPBHI']['blueside'][i] if invac \
                                        else airtovac(numpy.array(par['DAPBHI']['blueside'][i])),
                                par['DAPBHI']['redside'][i] if invac \
                                        else airtovac(numpy.array(par['DAPBHI']['redside'][i])),
                                           integrand=par['DAPBHI']['integrand'][i],
                                           order=par['DAPBHI']['order'][i]) ]

        ParDatabase.__init__(self, parlist)

        # Ensure that all indices are unique
        if len(numpy.unique(self.data['index'])) != self.nindx:
            raise ValueError('Indices in {0} database are not all unique!'.format(
                                                                            self.database['key']))

