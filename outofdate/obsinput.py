# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
r"""
Define a parameter instance that holds the input information needed to
run the DAP for a specific MaNGA observation.

Class usage examples
--------------------

To define a set of observation parameteters::

    from mangadap.par.obsinput import ObsInputPar
    p = ObsInputPar(plate=7443, ifudesign=12701, mode='CUBE', vel=6139,
                    vdisp=100, ell=0.3416, pa=150.24, reff=5.70)

Or declare the parameter object by reading an `SDSS-style parameter
file`_:

    from mangadap.par.obsinput import ObsInputPar
    p = ObsInputPar.from_par_file('mangadap-7443-12701-LOGCUBE.par')

Revision history
----------------

    | **15 Mar 2016**: Original implementation by K. Westfall (KBW)
    | **11 May 2016**: (KBW) Switch to using `pydl.pydlutils.yanny`_
        instead of internal yanny reader.
    | **24 Aug 2017**: (KBW) Provide new 'valid' attributes to signify
        if the input data was changed to the default value.
    | **11 Nov 2017**: (KBW) Change velocity minimum to -500 km/s to
        match :class:`mangadap.survey.drpcomplete.DRPComplete`

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
from pydl.pydlutils.yanny import yanny

from .parset import KeywordParSet
from ..drpfits import DRPFits

class ObsInputPar(KeywordParSet):
    r"""
    Parameter object that defines the input observation and guess
    parameters to be used by the DAP.

    The defined parameters are:

    .. include:: ../tables/obsinputpar.rst

    See :class:`mangadap.par.parset.ParSet` for attributes and raised
    exceptions.

    .. todo::
        Set default velocity to 0?

    Attributes:
        valid_vdisp (:obj:`bool`):
            Input velocity dispersion was >0.
        valid_ell (:obj:`bool`):
            Input ellipticiticy was within [0,1].
        valid_pa (:obj:`bool`):
            Input position angle was within [0,360).
        valid_reff (:obj:`bool`):
            Input effective radius was >0.
    """
    def __init__(self, plate=None, ifudesign=None, mode=None, vel=None, vdisp=None, ell=None,
                 pa=None, reff=None):

        in_fl = [ int, float ]
        mode_keys = DRPFits.mode_options()

        pars =     [ 'plate', 'ifudesign',    'mode', 'vel', 'vdisp', 'ell',  'pa', 'reff' ]
        values =   [   plate,   ifudesign,      mode,   vel,   vdisp,   ell,    pa,   reff ]
        defaults = [    None,        None,    'CUBE',  None,   100.0,   0.0,   0.0,    1.0 ]
        options =  [    None,        None, mode_keys,  None,    None,  None,  None,   None ]
        dtypes =   [     int,         int,       str, in_fl,   in_fl, in_fl, in_fl,  in_fl ]
        descr = [ 'Plate number', 'IFU designation',
                  'DRP 3D mode; see :func:`mangadap.drpfits.DRPFits.mode_options`.',
                  'Systemic velocity (km/s)', 'Guess velocity dispersion (km/s)',
                  r'Isophotal ellipticity (:math:`\varepsilon = 1-b/a`)',
                  'Position angle (degrees) of the isophotal major axis from N through E',
                  'Effective radius (arcsec)']

        super(ObsInputPar, self).__init__(pars, values=values, defaults=defaults, options=options,
                                          dtypes=dtypes, descr=descr)

        self.valid_vdisp = True
        self.valid_ell = True
        self.valid_pa = True
        self.valid_reff = True

        self._check()

    def _check(self):
        r"""
        Check that the values in the parameter list adhere to some
        hard-coded limits.

            - Velocity (:math:`cz`) and velocity dispersion have to be
              greater than -500.
            - Ellipticity has to be :math:`0 \leq \varepsilon < 1`
            - Position angle has to be :math:`0 \leq \phi < 360`
            - Effective radius has to be greater than zero.

        Raises:
            ValueError: Raised if velocity (:math:`cz`) is less than
                zero.
        """
        if self.data['vel'] is not None and not self.data['vel'] > -500:
            raise ValueError('Velocity must be > -500!')
        if self.data['vdisp'] < 0:
            warnings.warn('Velocity dispersion less than 0; using default of 100 km/s.')
            self.data['vdisp'] = 100.0
            self.valid_vdisp = False
        if self.data['ell'] < 0 or self.data['ell'] > 1:
            warnings.warn('Ellipticity must be 0 <= ell <= 1; setting to 0.')
            self.data['ell'] = 0.
            self.valid_ell = False
        if self.data['pa'] < 0 or self.data['pa'] >= 360:
            warnings.warn('Imposing 0-360 range on position angle.')
            while self.data['pa'] < 0:
                self.data['pa'] += 360
            while self.data['pa'] >= 360:
                self.data['pa'] -= 360
            self.valid_pa = False
        if not self.data['reff'] > 0:
            warnings.warn('Effective radius must be larger than 0; setting to unity.')
            self.data['reff'] = 1.0
            self.valid_reff = False

    # Add a "to_par_file" method?

    @classmethod
    def from_par_file(cls, f):
        """
        Read the observation parameters from the provided yanny file.

        Args:
            f (str) : Name of the file to read

        Returns:
            :class:`ObsInputPar`: Derived instance of
            :class:`mangadap.par.ParSet` with the input observational
            parameters of the DRP data product to analyze with the DAP.

        Raises:
            FileNotFoundError: Raised if the provided file does not
                exist.
            ValueError: Raised if the input yanny file has more than
                one entry of DAPPAR.
            KeyError: Raised if selected keys are not in provided file.

        """
        if not os.path.isfile(f):
            raise FileNotFoundError('Could not open {0}!'.format(f))
   
        # Read the file
        par = yanny(filename=f, raw=True)

        # Check the number of entries
        if len(par['DAPPAR']['plate']) > 1:
            raise ValueError('File must contain only instance of DAPPAR!')

        # Return the ObsInputPar instance
        return cls(par['DAPPAR']['plate'][0],
                   par['DAPPAR']['ifudesign'][0],
                   mode=par['DAPPAR']['mode'][0],
                   vel=par['DAPPAR']['vel'][0],
                   vdisp=par['DAPPAR']['vdisp'][0],
                   ell=par['DAPPAR']['ell'][0],
                   pa=par['DAPPAR']['pa'][0],
                   reff=par['DAPPAR']['reff'][0])

