# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Container class for a database of emission-line parameters, as well as
support classes and functions.

Class usage examples
--------------------

To define an emission line::

    from mangadap.par.emissionlinedb import EmissionLinePar
    p = EmissionLinePar(index=44, name='Ha', restwave=6564.632, action='f', line='l',
                        flux=1.0, vel=0.0, sig=10., mode='f')

More often, however, you will want to define an emission-line
database using an SDSS parameter file. To do so, you can use one of
the default set of available emission-line databases::

    from mangadap.par.emissionlinedb import EmissionLineDB
    print(EmissionLineDB.available_databases())
    emldb = EmissionLineDB.from_key('ELPMPL9')

The above call uses the
:func:`~mangadap.par.spectralfeaturedb.SpectralFeatureDB.from_key`
method to define the database using its keyword and the database
provided with the MaNGA DAP source distribution. You can also define
the database directly for an SDSS-style parameter file::

    from mangadap.par.emissionlinedb import EmissionLineDB
    emldb = EmissionLineDB('/path/to/emission/line/database/myeml.par')

The above will read the file and set the database keyword to
'MYEML' (i.e., the capitalized root name of the ``*.par`` file).
See :ref:`emissionlines` for the format of the parameter file.

----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

import os
import glob

from IPython import embed

import numpy

from pydl.goddard.astro import airtovac
from pydl.pydlutils.yanny import yanny

from .parset import KeywordParSet
from .spectralfeaturedb import SpectralFeatureDB
from ..proc import util
from ..util.datatable import DataTable


class EmissionLinePar(KeywordParSet):
    r"""
    Metadata used to define an emission line for the DAP to fit
    parametrically.

    These parameters are only for a single emission line; see
    :class:`~mangadap.par.emissionlinedb.EmissionLineDB` for setting a full
    line list of parameters.

    .. todo::
        - Point to/add a description the paramters.

    """
    def __init__(self, index=None, name=None, restwave=None, action=None, tie_index=None,
                 tie_par=None, blueside=None, redside=None):
        
        in_fl = [ int, float ]
        action_options = [ 'i', 'f', 'm', 's']
        arr_like = [ numpy.ndarray, list ]

        _name = None if name is None else name.strip()
        _action = None if action is None else action.strip()
        
        pars =     ['index', 'name', 'restwave', 'action', 'tie_index', 'tie_par', 'blueside',
                    'redside']
        values =   [index, _name, restwave, _action, tie_index, tie_par, blueside, redside]
        defaults = [None, None, None, 'f', None, None, None, None]
        options =  [None, None, None, action_options, None, None, None, None]
        dtypes =   [int, str, in_fl, str, int, arr_like, arr_like, arr_like]

        descr = ['An index used to refer to the line for tying parameters; must be >0.',
                 'A name for the line.',
                 'The rest wavelength of the line in angstroms *in vacuum*.',
                 'Describes how the line should be treated.  See ' \
                    ':ref:`emission-line-modeling-action`. Default is ``f``.',
                 'Index of the line to which parameters for this line are tied.  ' \
                    'See :ref:`emission-line-modeling-tying`.',
                 'Details of how each model parameter for this line is tied to the reference ' \
                    'line.  See :ref:`emission-line-modeling-tying`.',
                 'A two-element vector with the starting and ending wavelength for a bandpass ' \
                    'blueward of the emission line, which is used to set the continuum level ' \
                    'near the emission line when calculating the equivalent width.',
                 'A two-element vector with the starting and ending wavelength for a bandpass ' \
                    'redward of the emission line, which is used to set the continuum level ' \
                    'near the emission line when calculating the equivalent width.']

        super().__init__(pars, values=values, defaults=defaults, options=options, dtypes=dtypes,
                         descr=descr)

        self._check()

    def _check(self):
        if self.data['index'] is not None and self.data['index'] <= 0:
            raise ValueError('Index must be larger than 0.')


class EmissionLineDefinitionTable(DataTable):
    """
    A :class:`~mangadap.util.datatable.DataTable` with the data from an
    :class:`~mangadap.par.emissionlinedb.EmissionLineDB` object.

    Table includes:

    .. include:: ../tables/emissionlinedefinitiontable.rst
    
    Args:
        name_len (:obj:`int`, optional):
            The maximum length of any of the emission line names.
        tie_len (:obj:`int`, optional):
            The maximum length of a tying string. Not sure if this is needed.
        shape (:obj:`int`, :obj:`tuple`, optional):
            The shape of the initial array. If None, the data array
            will not be instantiated; use
            :func:`~mangadap.util.datatable.DataTable.init` to
            initialize the data array after instantiation.
    """
    def __init__(self, name_len=1, tie_len=1, shape=None):
        # NOTE: This requires python 3.7 to make sure that this is an "ordered"
        # dictionary.
        datamodel = dict(ID=dict(typ=int, shape=None, descr='Emission line ID number'),
                         NAME=dict(typ='<U{0:d}'.format(name_len), shape=None,
                                   descr='Name of the emission line'),
                         RESTWAVE=dict(typ=float, shape=None,
                                       descr='Rest wavelength of the emission line'),
                         ACTION=dict(typ='<U1', shape=None,
                                     descr='Action to take for this emission line; see '
                                           ':ref:`emission-line-modeling-action`'),
                         TIE_ID=dict(typ=int, shape=None,
                                     descr='ID of the line to which this one is tied.'),
                         TIE_FLUX=dict(typ='<U{0:d}'.format(tie_len), shape=None,
                                      descr='Tying parameter for the flux of each line.'),
                         TIE_VEL=dict(typ='<U{0:d}'.format(tie_len), shape=None,
                                      descr='Tying parameter for the velocity of each line.'),
                         TIE_SIG=dict(typ='<U{0:d}'.format(tie_len), shape=None,
                                      descr='Tying parameter for the dispersion of each line.'))

        keys = list(datamodel.keys())
        super().__init__(keys, [datamodel[k]['typ'] for k in keys],
                         element_shapes=[datamodel[k]['shape'] for k in keys],
                         descr=[datamodel[k]['descr'] for k in keys], shape=shape)


class EmissionLineDB(SpectralFeatureDB):
    r"""
    Basic container class for the database of emission-line parameters.

    Each row of the database is parsed using
    :class:`~mangadap.par.emissionlinedb.EmissionLinePar`. For the format of
    the input file, see :ref:`emissionlines-modeling`.

    The primary instantiation requires the SDSS parameter file with the
    emission-line data. To instantiate using a keyword, use the
    :func:`~mangadap.par.spectralfeaturedb.SpectralFeatureDB.from_key` method
    of the base class; see the base class for instantiation and additional
    attributes.

    Attributes:
        key (:obj:`str`):
            Database signifying keyword
        file (:obj:`str`):
            File with the emission-line data
        size (:obj:`int`):
            Number of emission lines in the database. 
    """

    default_data_dir = 'emission_lines'
    """
    Name of the directory in the mangadap data path with the distributed
    emission-line databases.
    """

    def _parse_yanny(self):
        """
        Parse the yanny file (provided by :attr:`file`) for the emission-line
        database.

        Returns:
            :obj:`list`: The list of :class:`~mangadap.par.parset.ParSet`
            instances for each line of the database.
        """
        # Read the yanny file
        par = yanny(filename=self.file, raw=True)
        if len(par['DAPEML']['index']) == 0:
            raise ValueError('Could not find DAPEML entries in {0}!'.format(self.file))

        # Setup the array of emission line database parameters
        self.size = len(par['DAPEML']['index'])
        parlist = []
        for i in range(self.size):
            invac = par['DAPEML']['waveref'][i] == 'vac'

            tie_index = -1 if par['DAPEML']['tie'][i][0] == 'None' \
                            else int(par['DAPEML']['tie'][i][0])
            tie_par = [None if t == 'None' else t for t in par['DAPEML']['tie'][i][1:]]

            parlist += [EmissionLinePar(index=par['DAPEML']['index'][i],
                                       name=par['DAPEML']['name'][i],
                                       restwave=par['DAPEML']['restwave'][i] if invac 
                                                else airtovac(par['DAPEML']['restwave'][i]),
                                       action=par['DAPEML']['action'][i], tie_index=tie_index,
                                       tie_par=tie_par,
                                       blueside=par['DAPEML']['blueside'][i] if invac else \
                                            airtovac(numpy.array(par['DAPEML']['blueside'][i])),
                                       redside=par['DAPEML']['redside'][i] if invac else \
                                            airtovac(numpy.array(par['DAPEML']['redside'][i])))]
        return parlist

    def channel_names(self, dicttype=True):
        """
        Return a dictionary with the channel names as the dictionary
        key and the channel number as the dictionary value. If
        ``dicttype`` is False, a list is returned with just the
        string keys.
        """
        channels = ['{0}-{1}'.format(self.data['name'][i], int(self.data['restwave'][i])) 
                        for i in range(self.size)]
        return {n:i for i,n in enumerate(channels)} if dicttype else channels

    def to_datatable(self, quiet=False):
        """
        Constuct an
        :class:`~mangadap.par.emissionlinedb.EmissionLineDefinitionTable`
        instance with the line database.

        Args:
            quiet (:obj:`bool`, optional):
                Suppress terminal output

        Returns:
            :class:`~mangadap.par.emissionlinedb.EmissionLineDefinitionTable`:
            Table with the emission-line metadata.
        """
        name_len = numpy.amax([len(n) for n in self.data['name']])
        tie_len = numpy.amax([0 if t is None else len(t) 
                                for t in numpy.ravel(self.data['tie_par'])])

        # Instatiate the table data that will be saved defining the set
        # of emission-line moments measured
        db = EmissionLineDefinitionTable(name_len=name_len, tie_len=tie_len, shape=self.size)

        hk = ['ID', 'NAME', 'RESTWAVE', 'ACTION', 'TIE_ID']
        mk = ['index', 'name', 'restwave', 'action', 'tie_index']
        for _hk, _mk in zip(hk,mk):
            db[_hk] = self.data[_mk]

        db['TIE_FLUX'] = numpy.array([str(d) for d in self.data['tie_par'][:,0]])
        db['TIE_VEL'] = numpy.array([str(d) for d in self.data['tie_par'][:,1]])
        db['TIE_SIG'] = numpy.array([str(d) for d in self.data['tie_par'][:,2]])

        return db

    def tie_index_match(self):
        """
        Return a vector with the row indices for the tied parameters, which
        is more useful than the ID number of the tied line.
        """
        tie_indx = numpy.full(self.size, -1, dtype=int)
        indx = self.data['tie_index'] > 0
        tie_indx[indx] = [numpy.where(self.data['index'] == i)[0][0]
                            for i in self.data['tie_index'][indx]]
        return tie_indx

    @property
    def neml(self):
        """
        Number of emission lines in the database.
        """
        return self.nsets


