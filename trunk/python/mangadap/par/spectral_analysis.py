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
from mangadap.proc.options import spectral_analysis_options

__author__ = 'Kyle B. Westfall'


class SpectralAnalysisPar(ParSet):
    """
    Base class with parameters used in the spectral analysis procedures.
    """
    def __init__(self, lamrange=None, snr_threshold=None, prior=None):

        # Get the analysis options
        analysis_keys = spectral_analysis_options()
        in_fl = [ int, float ]

        pars =   [ 'analysis_type',    'lamrange', 'snr_threshold',  'prior' ]
        values = [            None,      lamrange,   snr_threshold,    prior ]
        defaults = [          None,          None,             0.0,     None ]
        options = [  analysis_keys,          None,            None,     None ]
        dtypes = [             str, numpy.ndarray,           in_fl,      str ]
    
        ParSet.__init__(self, pars, values=values, defaults=defaults, options=options,
                        dtypes=dtypes)


class StellarContinuumAnalysisPar(SpectralAnalysisPar):
    """
    Derived class with parameters used in the stellar continuum fitting.
    Parameters currently mirror those in ppxf.

    pPXF parameters not included:
        Added in analysis call
            goodpixels, lam
        Ignored for now:
            plot, quiet, sky, vsyst, component, reg_dim

        reg_dim is always None because component is not yet allowed.

        template library options are left free!  will have to be able to
        find template library among defined set

    """
    def __init__(self, lamrange=None, snr_threshold=None, prior=None, tpl_lib='MILES',
                 eml_par='STANDARD', bias=None, clean=False, degree=-1, mdegree=0, moments=2,
                 oversample=False, regul=0, reddening=None):

        in_fl = [ int, float ]
        SpectralAnalysisPar.__init__(self, lamrange=lamrange, snr_threshold=snr_threshold,
                                        prior=prior)
        self['analysis_type'] = 'STRCNT'
        self.add('tpl_lib', tpl_lib, dtype='str')       # Template library keyword
        self.add('eml_par', eml_par, dtype='str')       # Emission-line database keyword
        self.add('bias', bias, dtype=in_fl)             # Bias fit
        self.add('clean', clean, dtype=bool)            # Clean using sigma clip
        self.add('degree', degree, dtype=int)           # Degree of additive polynomial
        self.add('mdegree', mdegree, dtype=int)         # Degree of multiplicative polynomial
        self.add('moments', moments, options=[2, 4, 6], dtype=int)    # Number of kinematic moments
        self.add('oversample', oversample, dtype=bool)  # Oversample the templates 30x
        self.add('regul', regul, dtype=in_fl)           # Regularization factor
        self.add('reddening', reddening, dtype=in_fl)   # Initial guess for reddening to fit


class EmissionLineAnalysisPar(SpectralAnalysisPar):
    """
    Derived class with parameters used in the emission-line fitting.

    Should the template library key be included?  Design thing... Should
    all the SpectralAnalysisPar objects be independent?

    allow for moments, LSF, etc

    """
    def __init__(lamrange=None, snr_threshold=None, prior=None, eml_par='STANDARD'):
        in_fl = [ int, float ]
        SpectralAnalysisPar.__init__(self, lamrange=lamrange, snr_threshold=snr_threshold,
                                        prior=prior)
        self['analysis_type'] = 'EMLINE'
        self.add('eml_par', eml_par, dtype='str')       # Emission-line database keyword


class SpectralIndexAnalysisPar(SpectralAnalysisPar):
    """
    Derived class with parameters used in the spectral-index measurements.
    """
    def __init__(self, lamrange=None, snr_threshold=None, prior=None, tpl_lib='MILES',
                 eml_par='STANDARD', sindx_par = 'LICK'):

        in_fl = [ int, float ]
        SpectralAnalysisPar.__init__(self, lamrange=lamrange, snr_threshold=snr_threshold,
                                        prior=prior)
        self['analysis_type'] = 'SPINDX'
        self.add('tpl_lib', tpl_lib, dtype='str')       # Template library keyword
        self.add('eml_par', eml_par, dtype='str')       # Emission-line database keyword
        self.add('sindx_par', sindx_par, dtype='str')   # Spectral-index database keyword


