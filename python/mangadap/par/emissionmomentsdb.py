# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Container class for the database of emission-line bandpass filters used
for non-parameteric measurements of emission-line flux and velocity
moments.

.. todo::
    Combine this with the main emissionlinemoments.py file.

Class usage examples
--------------------

Emission-line moment databases are defined using SDSS parameter
files. To define a database, you can use one of the default set of
available emission-line moment databases::

    from mangadap.par.emissionmomentsdb import EmissionMomentsDB
    print(EmissionMomentsDB.available_databases())
    elmom = EmissionMomentsDB.from_key('ELBMPL9')

The above call uses the :func:`EmissionMoments.from_key` method to
define the database using its keyword and the database provided with
the MaNGA DAP source distribution. You can also define the database
directly for an SDSS-style parameter file::

    from mangadap.par.emissionmomentsdb import EmissionMomentsDB
    elmom = EmissionMomentsDB('/path/to/emission/moments/database/myelm.par')

The above will read the file and set the database keyword to 'MYELM'
(i.e., the capitalized root name of the ``*.par`` file). See
:ref:`emissionlines` for the format of the parameter file.

Revision history
----------------

    | **17 Mar 2016**: Original implementation by K. Westfall (KBW)
    | **11 May 2016**: (KBW) Switch to using `pydl.pydlutils.yanny`_ and
        `pydl.goddard.astro.airtovac`_ instead of internal functions
    | **06 Oct 2017**: (KBW) Add function to return channel names.
    | **02 Dec 2019**: (KBW) Significantly reworked to use the new
        base class.

----

.. include license and copyright
.. include:: ../copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../links.rst
"""
import os
import warnings

import numpy

from pydl.goddard.astro import airtovac
from pydl.pydlutils.yanny import yanny

from .spectralfeaturedb import SpectralFeatureDB
from ..proc.bandpassfilter import BandPassFilterPar

class EmissionMomentsDB(SpectralFeatureDB):
    r"""
    Basic container class for the database of emission-line moments
    to calculate.

    Each row of the database is parsed using
    :class:`mangadap.proc.bandpassfilter.BandPassFilterPar`.  For the
    format of the input file, see :ref:`emissionlines-moments`.

    The primary instantiation requires the SDSS parameter file with
    the bandpass data. To instantiate using a keyword (and optionally
    a directory that holds the parameter files), use the
    :func:`mangadap.par.spectralfeaturedb.SpectralFeatureDB.from_key`
    class method.  See the base class for additional attributes.

    Args:
        parfile (:obj:`str`):
            The SDSS parameter file with the database.

    Attributes:
        key (:obj:`str`):
            Database signifying keyword
        file (:obj:`str`):
            File with the data
        size (:obj:`int`):
            Number of features in the database. 
        dummy (`numpy.ndarray`_):
            Boolean array flagging bandpasses as dummy placeholders.
    """
    default_data_dir = 'emission_bandpass_filters'
    def _parse_yanny(self):
        """
        Parse the yanny file (provided by :attr:`file`) for the
        emission-line moment database.

        Returns:
            :obj:`list`: The list of
            :class:`mangadap.par.parset.ParSet` instances for each
            line of the database.
        """
        # Read the yanny file
        par = yanny(filename=self.file, raw=True)
        if len(par['DAPELB']['index']) == 0:
            raise ValueError('Could not find DAPELB entries in {0}!'.format(self.file))

        # Check if any of the bands are dummy bands and warn the user
        self.dummy = numpy.any(numpy.array(par['DAPELB']['blueside']) < 0, axis=1)
        self.dummy |= numpy.any(numpy.array(par['DAPELB']['redside']) < 0, axis=1)
        self.dummy |= numpy.any(numpy.array(par['DAPELB']['primary']) < 0, axis=1)
        if numpy.sum(self.dummy) > 0:
            warnings.warn('Bands with negative wavelengths are used to insert dummy values.'
                          '  Ignoring input bands with indices: {0}'.format(
                                                numpy.array(par['DAPELB']['index'])[self.dummy]))

        # Setup the array of absorption-line index database parameters
        self.size = len(par['DAPELB']['index'])
        parlist = []
        for i in range(self.size):
            invac = par['DAPELB']['waveref'][i] == 'vac'
            parlist += [ BandPassFilterPar(index=par['DAPELB']['index'][i],
                                           name=par['DAPELB']['name'][i],
                                blueside=par['DAPELB']['blueside'][i] if invac \
                                        else airtovac(numpy.array(par['DAPELB']['blueside'][i])),
                                redside=par['DAPELB']['redside'][i] if invac \
                                        else airtovac(numpy.array(par['DAPELB']['redside'][i])),
                                restwave=par['DAPELB']['lambda'][i] if invac \
                                        else airtovac(numpy.array(par['DAPELB']['lambda'][i])),
                                primary=par['DAPELB']['primary'][i] if invac \
                                        else airtovac(numpy.array(par['DAPELB']['primary'][i]))) ]
        return parlist

    def channel_names(self, dicttype=True):
        """
        Return a dictionary with the channel names as the dictionary
        key and the channel number as the dictionary value. If
        ``dicttype`` is False, a list is returned with just the
        string keys.
        """
        channels = [ '{0}-{1}'.format(self.data['name'][i], int(self.data['restwave'][i])) 
                            for i in range(self.size) ]
        return { n:i for i,n in enumerate(channels) } if dicttype else channels
