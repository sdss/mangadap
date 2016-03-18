# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Container class that defines a bandpass filter.

*License*:
    Copyright (c) 2015, Kyle B. Westfall
    Licensed under BSD 3-clause license - see LICENSE.rst

*Source location*:
    $MANGADAP_DIR/python/mangadap/par/bandpassfilter.py

*Python2/3 compliance*::

    from __future__ import division
    from __future__ import print_function
    from __future__ import absolute_import
    from __future__ import unicode_literals
    
    import sys
    if sys.version > '3':
        long = int

*Imports*::

    import numpy
    from .parset import ParSet

*Class usage examples*:

    To define an bandpass filter::

        from mangadap.par.bandpassfilter import BandPassFilterPar
        p = BandPassFilterPar(44, 'Ha', [6483.0,6513.0],
                             [6623.0,6653.0], restwave=6564.632,
                             primary=[6557.6,6571.6]) 

    However, this class is mostly to provide a base class used by
    :class:`mangadap.proc.emissionmomentsdb.EmissionMomentsDB`,
    :class:`mangadap.proc.absorptionindexdb.AbsorptionIndexDB`, and
    :class:`mangadap.proc.bandheadindexdb.BandheadIndexDB`; it is not
    really meant to be used as given above.

..todo::
    - Add a function that prints the filter.

*Revision history*:
    | **18 Mar 2016**: Original implementation by K. Westfall (KBW)

"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import sys
if sys.version > '3':
    long = int

import numpy
from .parset import ParSet

__author__ = 'Kyle B. Westfall'
# Add strict versioning
# from distutils.version import StrictVersion

class BandPassFilterPar(ParSet):
    r"""

    Parameter object that defines a set of band-pass filters.  This is a
    base class for similar objects used calculating fluxes and flux
    moments of emission lines
    (:class:`mangadap.proc.emissionmomentsdb.EmissionMomentsDB`),
    absorption line indices
    (:class:`mangadap.proc.absorptionindexdb.AbsorptionIndexDB`), and
    spectral continuum bandheads
    (:class:`mangadap.proc.bandheadindexdb.BandheadIndexDB`) in the DAP.

    All wavelengths are expected to be IN VACUUM, to match the expected
    application of this filter to SDSS MaNGA spectra.

    See :class:`mangadap.par.parset.ParSet` for attributes and raised
    exceptions.

    Args:
        index (int) : An index used to refer to the bandpass filter.
        name (str) : A name for the bandpass filter.
        blueside (numpy.ndarray, list) : A two-element vector with the
            starting and ending wavelength (angstroms in VACUUM) for a
            passband to the blue of the primary spectral feature.
        redside (numpy.ndarray, list): A two-element vector with the
            starting and ending wavelength (angstroms in VACUUM) for a
            passband to the red of the primary spectral feature.
        restwave (float) : (Optional) The rest wavelength of the line in
            the primary passband in angstroms *in vacuum*.  This is used
            to convert the first and second moments to velocity (km/s)
            for emission lines.
        primary (numpy.ndarray, list) : (Optional) A two-element vector
            with the starting and ending wavelength (angstroms in
            VACUUM) for the primary passband surrounding the spectral
            feature of interest.  This is used by
            :class:`mangadap.proc.emissionmomentsdb.EmissionMomentsDB`
            and
            :class:`mangadap.proc.absorptionindexdb.AbsorptionIndexDB`,
            but it is irrelevant for
            :class:`mangadap.proc.bandheadindexdb.BandheadIndexDB`.
        units (str) : (Optional) Define the unit for the spectral index
            as either angstroms ('ang') or magnitudes ('mag').
            Currently only used by
            :class:`mangadap.proc.absorptionindexdb.AbsorptionIndexDB`.
        component (bool) : (Optional) Flag that the bandpass definition
            is a component of a larger set of bandpass filters.  This is
            currently only used by
            :class:`mangadap.proc.absorptionindexdb.AbsorptionIndexDB`
            to combine index measurements into a single index.  If True,
            all components with the same *name* are summed to form the
            composite index.
        integrand (str) : (Optional) Currently only used by
            :class:`mangadap.proc.bandheadindexdb.BandheadIndexDB`.
            Define the integrand over the passband used to construct and
            index as either flux per unit frequency (``'fnu'``),
            :math:`F_\nu`, or flux per unit wavelength (``'flambda'``),
            :math:`F_\lambda`.
        order (str): (Optional) Currently only used by
            :class:`mangadap.proc.bandheadindexdb.BandheadIndexDB`.
            Define the order to use when constructing the index.  The
            options are either a ratio of red-to-blue or blue-to-red,
            which are respectively selected using ``'r_b'`` or
            ``'b_r'``.

    """
    def __init__(self, index, name, blueside, redside, restwave=None, primary=None, units=None,
                 component=None, integrand=None, order=None):
        
        in_fl = [ int, float ]
        ar_like = [ numpy.ndarray, list ]
        n = name.strip()
        wave_opts = [ 'vac', 'air' ]
        unit_opts = [ 'ang', 'mag' ]
        integ_opts = [ 'flambda', 'fnu' ]
        order_opts = [ 'b_r', 'r_b' ]
        
        pars  =    [ 'index', 'name', 'restwave', 'primary', 'blueside', 'redside',   'units',
                     'component', 'integrand',    'order' ]
        values  =  [   index,      n,   restwave,   primary,   blueside,   redside,     units,
                       component,   integrand,      order ]
        options  = [    None,   None,       None,      None,       None,      None, unit_opts,
                            None,  integ_opts, order_opts ]
        dtypes  =  [     int,    str,      in_fl,   ar_like,    ar_like,   ar_like,       str,
                            bool,         str,        str ]

        ParSet.__init__(self, pars, values=values, options=options, dtypes=dtypes)
        self._check()


    def __repr__(self):
        """Return a string representation of the bandpass filter."""
        summary = 'index     : {0}'.format(self['index'])
        summary += '\nname      : {0}'.format(self['name'])
        if self['restwave'] is not None:
            summary += '\nrestwave  : {0}'.format(self['restwave'])
        if self['primary'] is not None:
            summary += '\nprimary   : {0}'.format(self['primary'])
        summary += '\nblueside  : {0}'.format(self['blueside'])
        summary += '\nredside   : {0}'.format(self['redside'])
        if self['units'] is not None:
            summary += '\nunits     : {0}'.format(self['units'])
        if self['component'] is not None:
            summary += '\ncomponent : {0}'.format(self['component'])
        if self['integrand'] is not None:
            summary += '\nintegrand : {0}'.format(self['integrand'])
        if self['order'] is not None:
            summary += '\norder     : {0}'.format(self['order'])
        return summary


    def _check(self):
        """
        Check the parameter list.
            - Make sure arrays only have two elements.

        .. todo::
            - Add check to __setitem__()?

        Raises:
            ValueError : Raised if one of the conditions above are not
                met.
        """
        if self.data['primary'] is not None and len(self.data['primary']) != 2:
            raise ValueError('Primary passband must have two and only two elements.')
        if len(self.data['blueside']) != 2:
            raise ValueError('Blue sideband must have two and only two elements.')
        if len(self.data['redside']) != 2:
            raise ValueError('Red sideband must have two and only two elements.')



