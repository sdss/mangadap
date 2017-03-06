# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Container class for the database of emission-line bandpass filters used
for non-parameteric measurements of emission-line flux and velocity
moments.

.. todo::
    Combine this with the main emissionlinemoments.py file.

*License*:
    Copyright (c) 2015, SDSS-IV/MaNGA Pipeline Group
        Licensed under BSD 3-clause license - see LICENSE.rst

*Source location*:
    $MANGADAP_DIR/python/mangadap/par/emissionmomentsdb.py

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
        from .parset import ParDatabase
        from .spectralfeaturedb import available_spectral_feature_databases, SpectralFeatureDBDef
        from ..proc.bandpassfilter import BandPassFilterPar
        from ..proc.util import _select_proc_method

*Class usage examples*:
    Emission-line moment databases are defined using SDSS parameter
    files.  To define a database, you can use one of the default set of
    available emission-line moment databases (see
    :func:`available_emission_bandpass_filter_databases`)::
    
        from mangadap.par.emissionmomentsdb import EmissionMomentsDB
        p = EmissionMomentsDB('STRONG')

    The above call requires that the ``$MANGADAP_DIR`` environmental
    variable is set.  If it is not, you can define it's location, as
    in::

        from mangadap.par.emissionmomentsdb import EmissionMomentsDB
        p = EmissionMomentsDB('STRONG', dapsrc='/path/to/dap/source')

    Finally, you can create your own `SDSS-style parameter file`_ with
    your own emission line passbands to use.  Example files are provided
    in ``$MANGADAP_DIR/data/emission_bandpass_filters`` with a
    companion ``README`` file.  With your own file, you have to point to
    the file using :class:`SpectralFeatureDBDef`, which you can then
    pass to :class:`EmissionMomentsDB`::

        from mangadap.par.spectralfeaturedb import SpectralFeatureDBDef
        from mangadap.par.emissionmomentsdb import EmissionMomentsDB
        d = SpectralFeatureDBDef(key='USER',
                                 file_path='/path/to/parameter/file')
        p = EmissionMomentsDB('USER', emldb_list=d)

    The reference frame of the emission-line passband wavelengths must
    be defined as either vacuum or air, using 'in_vacuum'.  It is
    expected that the object spectra to be analyzed are calibrated to
    vacuum wavelengths.  If 'in_vacuum' is false, this class will use
    :func:`mangadap.util.idlutils.airtovac` to convert the emission-line
    bandpass-filter wavelengths to vacuum.

*Revision history*:
    | **17 Mar 2016**: Original implementation by K. Westfall (KBW)
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

from pydl.goddard.astro import airtovac
from pydl.pydlutils.yanny import yanny
from .parset import ParDatabase
from .spectralfeaturedb import available_spectral_feature_databases, SpectralFeatureDBDef
from ..proc.bandpassfilter import BandPassFilterPar
from ..proc.util import _select_proc_method

__author__ = 'Kyle B. Westfall'
# Add strict versioning
# from distutils.version import StrictVersion

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

    This is a simple wrapper for
    :func:`mangadap.par.spectralfeaturedb.available_spectral_feature_databases`.

    Args:
        dapsrc (str): (**Optional**) Root path to the DAP source
            directory.  If not provided, the default is defined by
            :func:`mangadap.config.defaults.default_dap_source`.

    Returns:
        list: An list of :class:`SpectralFeatureDBDef` objects, each of
        which defines a unique set of absorption-line indices (see
        :class:`mangadap.proc.bandpassfilter.BandPassFilterPar`).

    .. todo::
        - Add backup function for Python 2.
        - Somehow add a python call that reads the databases and
          constructs the table for presentation in sphinx so that the
          text above doesn't have to be edited with changes in the
          available databases.

    """
    return available_spectral_feature_databases('emission_bandpass_filters', dapsrc=dapsrc)
    

class EmissionMomentsDB(ParDatabase):
    """
    Basic container class for the database of emission-line bandpass
    filter parameters.  See :class:`mangadap.parset.ParDatabase` for
    additional attributes.

    Args:
        database_key (str): Keyword selecting the database to use.
        emldb_list (list): (**Optional**) List of
            :class:`SpectralFeatureDBDef` objects that define the
            parameters required to setup the emission-line bandpass
            filter database.  The *database_key* must select one of the
            objects in this list.
        dapsrc (str): (**Optional**) Root path to the DAP source
            directory.  If not provided, the default is defined by
            :func:`mangadap.config.defaults.default_dap_source`.

    Attributes:
        version (str): Version number
        database (str): Keyword of the selected database to use.
        neml (int): Number of emission-line bandpass filters in the
            database
    """
    def __init__(self, database_key, emldb_list=None, dapsrc=None):

        self.version = '1.0'

        # Get the details of the selected database
        self.database = _select_proc_method(database_key, SpectralFeatureDBDef,
                                            method_list=emldb_list,
                                        available_func=available_emission_bandpass_filter_databases,
                                            dapsrc=dapsrc)
        
        # Check that the database exists
        if not os.path.isfile(self.database['file_path']):
            raise FileNotFoundError('Database file {0} does not exist!'.format(
                                                                    self.database['file_path']))

        # Read the yanny file
#        par = yanny(self.database['file_path'])
        par = yanny(filename=self.database['file_path'], raw=True)
        if len(par['DAPELB']['index']) == 0:
            raise ValueError('Could not find DAPELB entries in {0}!'.format(
                                                                    self.database['file_path']))

        # Check if any of the bands are dummy bands and warn the user
        self.dummy = numpy.any(numpy.array(par['DAPELB']['blueside']) < 0, axis=1)
        self.dummy |= numpy.any(numpy.array(par['DAPELB']['redside']) < 0, axis=1)
        self.dummy |= numpy.any(numpy.array(par['DAPELB']['primary']) < 0, axis=1)
        if numpy.sum(self.dummy) > 0:
            warnings.warn('Bands with negative wavelengths are used to insert dummy values.'
                          '  Ignoring input bands with indices: {0}'.format(
                                                numpy.array(par['DAPELB']['index'])[self.dummy]))

        # Setup the array of absorption-line index database parameters
        self.neml = len(par['DAPELB']['index'])
        parlist = []
        for i in range(self.neml):
            invac = par['DAPELB']['waveref'][i] == 'vac'
            parlist += [ BandPassFilterPar(par['DAPELB']['index'][i],
                                           par['DAPELB']['name'][i],
                                par['DAPELB']['blueside'][i] if invac \
                                        else airtovac(numpy.array(par['DAPELB']['blueside'][i])),
                                par['DAPELB']['redside'][i] if invac \
                                        else airtovac(numpy.array(par['DAPELB']['redside'][i])),
                                restwave=par['DAPELB']['lambda'][i] if invac \
                                        else airtovac(numpy.array(par['DAPELB']['lambda'][i])),
                                primary=par['DAPELB']['primary'][i] if invac \
                                        else airtovac(numpy.array(par['DAPELB']['primary'][i]))) ]

        ParDatabase.__init__(self, parlist)

        # Ensure that all indices are unique
        if len(numpy.unique(self.data['index'])) != self.neml:
            raise ValueError('Indices in {0} database are not all unique!'.format(
                                                                            self.database['key']))

