# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Container class for the database of absorption-line indices to measure.

Class usage examples
--------------------

Absorption-line index databases are defined using SDSS parameter files. To
define a database, you can use one of the default set of available
absorption-line index databases::

    from mangadap.par.absorptionindexdb import AbsorptionIndexDB
    print(AbsorptionIndexDB.available_databases())
    absdb = AbsorptionIndexDB.from_key('LICKINDX')

The above call uses the
:func:`~mangadap.par.spectralfeaturedb.SpectralFeatureDB.from_key`
method to define the database using its keyword and the database
provided with the MaNGA DAP source distribution. You can also define
the database directly for an SDSS-style parameter file::

    from mangadap.par.absorptionindexdb import AbsorptionIndexDB
    absdb = AbsorptionIndexDB('/path/to/absorption/index/database/myabs.par')

The above will read the file and set the database keyword to
'MYABS' (i.e., the capitalized root name of the ``*.par`` file).
See :ref:`spectralindices` for the format of the parameter file.
    
----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import os
import warnings
import numpy

from pydl.goddard.astro import airtovac
from pydl.pydlutils.yanny import yanny

from .spectralfeaturedb import SpectralFeatureDB
from ..proc.bandpassfilter import BandPassFilterPar

class AbsorptionIndexDB(SpectralFeatureDB):
    r"""
    Basic container class for the database of absorption-line
    indices.

    Each row of the database is parsed using
    :class:`mangadap.proc.bandpassfilter.BandPassFilterPar`.  For the
    format of the input file, see :ref:`spectralindices-absorption`.

    The primary instantiation requires the SDSS parameter file with
    the bandpass data. To instantiate using a keyword (and optionally
    a directory that holds the parameter files), use the
    :func:`~mangadap.par.spectralfeaturedb.SpectralFeatureDB.from_key`
    class method. See the base class for additional attributes.

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
    default_data_dir = 'absorption_indices'
    def _parse_yanny(self):
        """
        Parse the yanny file (provided by :attr:`file`) for the
        bandhead database.

        Returns:
            :obj:`list`: The list of
            :class:`mangadap.par.parset.ParSet` instances for each
            line of the database.
        """
        # Read the yanny file
        par = yanny(filename=self.file, raw=True)
        if len(par['DAPABI']['index']) == 0:
            raise ValueError('Could not find DAPABI entries in {0}!'.format(self.file))

        # Check if any of the bands are dummy bands and warn the user
        self.dummy = numpy.any(numpy.array(par['DAPABI']['blueside']) < 0, axis=1)
        self.dummy |= numpy.any(numpy.array(par['DAPABI']['redside']) < 0, axis=1)
        self.dummy |= numpy.any(numpy.array(par['DAPABI']['primary']) < 0, axis=1)
        if numpy.sum(self.dummy) > 0:
            warnings.warn('Bands with negative wavelengths are used to insert dummy values.'
                          '  Ignoring input bands with indices: {0}'.format(
                                                numpy.array(par['DAPABI']['index'])[self.dummy]))

        # Setup the array of absorption-line index database parameters
        self.size = len(par['DAPABI']['index'])
        parlist = []
        for i in range(self.size):
            invac = par['DAPABI']['waveref'][i] == 'vac'
            comp = par['DAPABI']['component'][i] != 0
            parlist += [ BandPassFilterPar(index=par['DAPABI']['index'][i],
                                           name=par['DAPABI']['name'][i],
                                blueside=par['DAPABI']['blueside'][i] if invac \
                                        else airtovac(numpy.array(par['DAPABI']['blueside'][i])),
                                redside=par['DAPABI']['redside'][i] if invac \
                                        else airtovac(numpy.array(par['DAPABI']['redside'][i])),
                                primary=par['DAPABI']['primary'][i] if invac \
                                        else airtovac(numpy.array(par['DAPABI']['primary'][i])),
                                           units=par['DAPABI']['units'][i],
                                           integrand='flambda',
                                           component=comp) ]
        return parlist

    def channel_names(self, offset=0):
        """
        Return a dictionary with the channel names as the dictionary
        key and the channel number as the dictionary value. An
        ``offset`` can be added to the channel number; i.e., if the
        offset is 2, the channel numbers will be a running number
        starting with 2.
        """
        return {self.data['name'][i] : i + offset for i in range(self.size)}
