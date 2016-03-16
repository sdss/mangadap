# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
r"""
Define a parameter instance that holds the input information needed to
run the DAP for a specific MaNGA observation.

*Source location*:
    $MANGADAP_DIR/python/mangadap/par/obsinput.py

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
    import os.path
    from mangadap.par.parset import ParSet
    from mangadap.util.yanny import yanny
    from mangadap.proc.options import drp_3dmode_options

*Class usage examples*:

    To define a set of observation parameteters::

        from mangadap.par.obsinput import ObsInputPar
        p = ObsInputPar(plate=7443, ifudesign=12701, mode='CUBE', vel=6139,
                        vdisp=100, ell=0.3416, pa=150.24, reff=5.70)

    Or declare the parameter object by reading a yanny parameter file
    using the helper function :func:`read_obs_input`::

        from mangadap.par.obsinput import read_obs_input
        p = read_obs_input('mangadap-7443-12701-LOGCUBE.par')

*Revision history*:
    | **15 Mar 2016**: Original implementation by K. Westfall (KBW)

"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import sys
if sys.version > '3':
    long = int

import numpy
import os.path
from mangadap.par.parset import ParSet
from mangadap.util.yanny import yanny
from mangadap.proc.options import drp_3dmode_options
import warnings

__author__ = 'Kyle B. Westfall'


class ObsInputPar(ParSet):
    """
    Parameter object that defines the input observation and guess
    parameters to be used by the DAP.

    Although technically optional in this version, use of the class by
    the DAP requires that *vel* must be defined.

    See :class:`mangadap.par.parset.ParSet` for attributes and raised
    exceptions.

    Args:
        plate (int) : Plate number; default is None.
        ifudesign (int) : IFU designation; default is None.
        mode (str) : DRP 3D mode; see
            :func:`mangadap.proc.options.drp_3dmode_options`; default is
            CUBE.
        vel (float) : Systemic velocity (km/s); default is None.
        vdisp (float) : Guess velocity dispersion (km/s); default is
            100.
        ell (float) : Isophotal ellipticity (:math:`\varepsilon =
            1-b/a`); default is 0 (circular).
        pa (float) : Position angle (degrees) of the isophotal major
            axis; default is 0.
        reff (float) : Effective radius (arcsec); default is 1.0.

    """
    def __init__(self, plate, ifudesign, mode=None, vel=None, vdisp=None, ell=None, pa=None,
                 reff=None):

        in_fl = [ int, float ]
        mode_keys = drp_3dmode_options()

        pars =     [ 'plate', 'ifudesign',    'mode', 'vel', 'vdisp', 'ell',  'pa', 'reff' ]
        values =   [   plate,   ifudesign,      mode,   vel,   vdisp,   ell,    pa,   reff ]
        defaults = [    None,        None,    'CUBE',  None,   100.0,   0.0,   0.0,    1.0 ]
        options =  [    None,        None, mode_keys,  None,    None,  None,  None,   None ]
        dtypes =   [     int,         int,       str, in_fl,   in_fl, in_fl, in_fl,  in_fl ]

        ParSet.__init__(self, pars, values=values, defaults=defaults, options=options,
                        dtypes=dtypes)
        self._check()


    def _check(self):
        """
        Check that the values in the parameter list adhere to some
        hard-coded limits.

            - Velocity (:math:`cz`) and velocity dispersion have to be
              greater than 0.
            - Ellipticity has to be :math:`0 \leq \varepsilon < 1'
            - Position angle has to be :math:`0 \leq \phi < 360`
            - Effective radius has to be greater than zero.

        Raises:
            ValueError : Raised if velocity (:math:`cz`) is less than
                zero.
        """
        if self.data['vel'] < 0:
            raise ValueError('Velocity must be > 0!')
        if self.data['vdisp'] < 0:
            warnings.warn('Velocity dispersion less than 0; using default of 100 km/s.')
            self.data['vdisp'] = 100.0
        if self.data['ell'] < 0:
            warnings.warn('Ellipticity less than 0; setting to 0.')
            self.data['ell'] = 0.
        if not self.data['ell'] < 1:
            warnings.warn('Ellipticity greater or equal to 1; setting to 0.')
            self.data['ell'] = 0.
        if self.data['pa'] < 0 or self.data['pa'] >= 360:
            warnings.warn('Imposing 0-360 range on position angle.')
            while self.data['pa'] < 0:
                self.data['pa'] += 360
            while self.data['pa'] >= 360:
                self.data['pa'] -= 360
        if not self.data['reff'] > 0:
            warnings.warn('Effective radius must be larger than 0; setting to unity.')
            self.data['reff'] = 1.0


def read_obs_input(f):
    """ Read the observation parameters from the provided yanny file.

    Args:
        f (str) : Name of the file to read

    Returns:
        :class:`ObsInputPar` : Derived instance of
            :class:`mangadap.par.ParSet` with the input observational
            parameters of the DRP data product to analyze with the DAP.

    Raises:
        FileNotFoundError : Raised if the provided file does not exist.

        ValueError : Raised if the input yanny file has more than one
            entry of DAPPAR.

        KeyError : Raised if selected keys are not in provided file.

    """
    if not os.path.isfile(f):
        raise FileNotFoundError('Could not open {0}!'.format(f))
   
    # Read the file
    par = yanny(f)

    # Check the number of entries
    if len(par['DAPPAR']['plate']) > 1:
        raise ValueError('File must contain only instance of DAPPAR!')

    # Return the ObsInputPar instance
    return ObsInputPar(
                        par['DAPPAR']['plate'][0],
                        par['DAPPAR']['ifudesign'][0],
                        mode=par['DAPPAR']['mode'][0],
                        vel=par['DAPPAR']['vel'][0],
                        vdisp=par['DAPPAR']['vdisp'][0],
                        ell=par['DAPPAR']['ell'][0],
                        pa=par['DAPPAR']['pa'][0],
                        reff=par['DAPPAR']['reff'][0]
                      )
    


