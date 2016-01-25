"""

Define some utility classes used to hold parameters.

*Source location*:
    $MANGADAP_DIR/python/mangadap/util/par.py

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
    from mangadap.util.options import binning_options, bin_weighting_options
    from mangadap.util.options import spectral_analysis_options

*Class usage examples*:

    .. todo::

        Add some usage comments here!

*Revision history*:
    | **16 Jun 2015**: Original implementation by K. Westfall (KBW)

"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import sys
if sys.version > '3':
    long = int

import numpy
from mangadap.par.ParSet import ParSet
from mangadap.proc.options import binning_options, bin_weighting_options

__author__ = 'Kyle B. Westfall'


class SpatialBinningPar(ParSet):
    """
    Base class with parameters used in the spatial binning procedures.
    """
    def __init__(self, snr_threshold=None, v_register=None, weighting=None, noise_calib=None):

        in_fl = [ int, float ]

        # Get the binning options
        bin_keys = binning_options()
        # Get the weighting options
        wgt_keys = bin_weighting_options()

        pars =     [ 'bin_type', 'snr_threshold', 'v_register', 'weighting', 'noise_calib' ]
        values =   [       None,   snr_threshold,   v_register,   weighting,   noise_calib ]
        defaults = [       None,             0.0,        False,   'Uniform',          True ]
        options =  [   bin_keys,            None,         None,    wgt_keys,          None ]
        dtypes =   [        str,           in_fl,         bool,         str,          bool ]

        ParSet.__init__(self, pars, values=values, defaults=defaults, options=options,
                        dtypes=dtypes)


class NoBinningPar(SpatialBinningPar):
    """
    Derived class with the parameters needed to perform the NONE binning.
    """
    def __init__(self, snr_threshold=None):
        SpatialBinningPar.__init__(self, snr_threshold=snr_threshold)
        self['bin_type'] = 'NONE'

        # Delete the unneeded keys
        del self['v_register']
        del self['weighting']
        del self['noise_calib']


class AllBinningPar(SpatialBinningPar):
    """
    Derived class with the parameters needed to perform the NONE binning.
    """
    def __init__(self, snr_threshold=None, v_register=None, weighting=None, noise_calib=None):
        SpatialBinningPar.__init__(self, snr_threshold=snr_threshold, v_register=v_register,
                                         weighting=weighting, noise_calib=noise_calib)
        self['bin_type'] = 'ALL'


class SNRBinningPar(SpatialBinningPar):
    """
    Derived class with the parameters needed to perform the SNR (Voronoi) binning
    """
    def __init__(self, snr_threshold=None, v_register=None, weighting=None, noise_calib=None,
                 target_snr=5.):
        SpatialBinningPar.__init__(self, snr_threshold=snr_threshold, v_register=v_register,
                                         weighting=weighting, noise_calib=noise_calib)
        self['bin_type'] = 'SNR'
        self.add('target_snr', target_snr, dtype=[ int, float ])    # Target S/N


class RadialBinningPar(SpatialBinningPar):
    """
    Derived class with the parameters needed to perform the RADIAL binning.
    """
    def __init__(self, snr_threshold=None, v_register=None, weighting=None, noise_calib=None,
                       cx=0.0, cy=0.0, pa=0.0, ell=0.0, rs=0.0, re=-1.0, nr=5, rlog=False,
                       rscale=1.0):
        in_fl = [ int, float ]
        SpatialBinningPar.__init__(self, snr_threshold=snr_threshold, v_register=v_register,
                                         weighting=weighting, noise_calib=noise_calib)
        self['bin_type'] = 'RADIAL'
        self.add('cx', cx, dtype=[ int, float ])            # X-center
        self.add('cy', cy, dtype=[ int, float ])            # Y-center
        self.add('pa', pa, dtype=[ int, float ])            # Position angle
        self.add('ell', ell, dtype=[ int, float ])          # Ellipticity (1-b/a)
        self.add('rs', rs, dtype=[ int, float ])            # Starting radius
        self.add('re', re, dtype=[ int, float ])            # Ending radius
        self.add('nr', nr, dtype=int)                       # Number of radial bins
        self.add('rlog', rlog, dtype=bool)                  # Geometrically space the bins
        self.add('rscale', rscale, dtype=[ int, float ])    # Scale factor for radius


