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
from mangadap.util.options import binning_options, bin_weighting_options
from mangadap.util.options import spectral_analysis_options

__author__ = 'Kyle B. Westfall'

class ParSet:
    """
    Generic base class to handle and manipulate a list of operational
    parameters.
    """
    def __init__(self, pars, values=None, defaults=None, options=None, dtypes=None):
        # Check that the list of input parameters is a list of strings
        if type(pars) != list:
            raise Exception('Input parameter keys must be provided as a list.')
        for key in pars:
            if type(key) != str:
                raise Exception('Input parameter keys must be strings.')
        
        # Get the length of the parameter list and make sure the list
        # has unique values
        self.npar = len(pars)
        if len(numpy.unique(numpy.array(pars))) != self.npar:
            raise Exception('All input parameter keys must be unique.')

        # Check that the other lists, if provided, have the correct type
        # and length
        if values is not None and (type(values) != list or len(values) != self.npar):
            raise Exception('Values must be a list with the same length as the number of keys.')
        if defaults is not None and (type(defaults) != list or len(defaults) != self.npar):
            raise Exception('Defaults must be a list with the same length as the number of keys.')
        if options is not None and (type(options) != list or len(options) != self.npar):
            raise Exception('Options must be a list with the same length as the number of keys.')
        if dtypes is not None and (type(dtypes) != list or len(dtypes) != self.npar):
            raise Exception('Data types must be a list with the same length as the number of keys.')

        # Set up dummy lists for no input
        if values is None:
            values = [None]*self.npar
        if defaults is None:
            defaults = [None]*self.npar
        if options is None:
            options = [None]*self.npar
        if dtypes is None:
            dtypes = [None]*self.npar

        # Set the valid options
        self.options = dict([ (p, [o]) if o is not None and type(o) != list else (p, o) \
                                       for p, o in zip(pars, options) ])
        # Set the valid types
        self.dtype = dict([ (p, [t]) if t is not None and type(t) != list else (p, t) \
                                     for p, t in zip(pars, dtypes) ])

        # Set the data dictionary using the internal functions
        self.data = {}
        for p, d, v in zip(pars, defaults, values):
            if v is None:
                self.__setitem__(p, d)
                continue
            self.__setitem__(p, v)


    def __getitem__(self, key):
        return self.data[key]


    def __setitem__(self, key, value):
        if value is None:
            self.data[key] = value
            return

        if self.options[key] is not None and value not in self.options[key]:
            raise ValueError('Input value invalid: {0}.\nOptions are: {1}'.format(value,
                                                                                self.options[key]))
        if self.dtype[key] is not None and type(value) not in self.dtype[key]:
            raise ValueError('Input value incorrect type: {0}.\nOptions are: {1}'.format(value,
                                                                                self.dtype[key]))
        self.data[key] = value
        

    def __iter__(self):
        return iter(self.data.values())

    
    def add(self, key, value, options=None, dtype=None):
        self.options[key] = [options] if options is not None and type(options) != list else options
        self.dtype[key] = [dtype] if dtype is not None and type(dtype) != list else dtype
        self.__setitem__(key, value)


class TemplateLibraryParSet(ParSet):
    """
    Class with parameters used to define the template library.
    Options and defaults in ParSet base class are set to None.
    """
    def __init__(self, key, file_search, fwhm, in_vacuum, wave_limit, lower_flux_limit): 
        # Perform some checks of the input
        in_fl = [ int, float ]
        
        pars =   [ 'key', 'file_search', 'fwhm', 'in_vacuum',  'wave_limit', 'lower_flux_limit' ]
        values = [   key,   file_search,   fwhm,   in_vacuum,    wave_limit,   lower_flux_limit ]
        dtypes = [   str,           str,  in_fl,        bool, numpy.ndarray,              in_fl ]

        ParSet.__init__(self, pars, values=values, dtypes=dtypes)


class EmissionLineParSet(ParSet):
    """
    Class with parameters used to define an emission-line database.
    Options and defaults in ParSet base class are set to None.
    """
    def __init__(self, key, file_path):

        pars =   [ 'key', 'file_path' ]
        values = [   key,   file_path ]
        dtypes = [   str,         str ]

        ParSet.__init__(self, pars, values=values, dtypes=dtypes)


class SpectralIndexParSet(ParSet):
    """
    Class with parameters used to define a spectral-index database.
    Options and defaults in ParSet base class are set to None.
    """
    def __init__(self, key, file_path):

        pars =   [ 'key', 'file_path' ]
        values = [   key,   file_path ]
        dtypes = [   str,         str ]

        ParSet.__init__(self, pars, values=values, dtypes=dtypes)


class SpatialBinningParSet(ParSet):
    """
    Base class with parameters used in the spatial binning procedures.
    """
    def __init__(self, snr_lamrange=None, snr_threshold=None, v_register=None, weighting=None,
                 noise_calib=None):

        in_fl = [ int, float ]

        # Get the binning options
        bin_keys = binning_options()
        # Get the weighting options
        wgt_keys = bin_weighting_options()

        pars =     [ 'bin_type',                    'snr_lamrange', 'snr_threshold', 'v_register',
                     'weighting', 'noise_calib' ]
        values =   [       None,                      snr_lamrange,   snr_threshold,   v_register,
                       weighting,   noise_calib ]
        defaults = [       None, numpy.ndarray([5560.00, 6942.00]),            0.0,         False,
                       'Uniform',          True ]
        options =  [   bin_keys,                              None,            None,         None,
                        wgt_keys,          None ]
        dtypes =   [        str,                     numpy.ndarray,           in_fl,         bool,
                             str,          bool ]

        ParSet.__init__(self, pars, values=values, defaults=defaults, options=options,
                        dtypes=dtypes)


class NoBinningParSet(SpatialBinningParSet):
    """
    Derived class with the parameters needed to perform the NONE binning.
    """
    def __init__(self, snr_lamrange=None, snr_threshold=None):
        SpatialBinningParSet.__init__(self, snr_lamrange=snr_lamrange, snr_threshold=snr_threshold)
        self['bin_type'] = 'NONE'


class AllBinningParSet(SpatialBinningParSet):
    """
    Derived class with the parameters needed to perform the NONE binning.
    """
    def __init__(self, snr_lamrange=None, snr_threshold=None, v_register=None, weighting=None,
                 noise_calib=None):
        SpatialBinningParSet.__init__(self, snr_lamrange=snr_lamrange, snr_threshold=snr_threshold,
                                      v_register=v_register, weighting=weighting,
                                      noise_calib=noise_calib)
        self['bin_type'] = 'ALL'


class SNRBinningParSet(SpatialBinningParSet):
    """
    Derived class with the parameters needed to perform the SNR (Voronoi) binning
    """
    def __init__(self, snr_lamrange=None, snr_threshold=None, v_register=None, weighting=None,
                 noise_calib=None, target_snr=5.):
        SpatialBinningParSet.__init__(self, snr_lamrange=snr_lamrange, snr_threshold=snr_threshold,
                                      v_register=v_register, weighting=weighting,
                                      noise_calib=noise_calib)
        self['bin_type'] = 'SNR'
        self.add('target_snr', target_snr, dtype=[ int, float ])    # Target S/N


class RadialBinningParSet(SpatialBinningParSet):
    """
    Derived class with the parameters needed to perform the RADIAL binning.
    """
    def __init__(self, snr_lamrange=None, snr_threshold=None, v_register=None, weighting=None,
                 noise_calib=None, cx=0.0, cy=0.0, pa=0.0, ell=0.0, rs=0.0, re=-1.0, nr=5,
                 rlog=False, rscale=1.0):
        in_fl = [ int, float ]
        SpatialBinningParSet.__init__(self, snr_lamrange=snr_lamrange, snr_threshold=snr_threshold,
                                      v_register=v_register, weighting=weighting,
                                      noise_calib=noise_calib)
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

class SpectralAnalysisParSet(ParSet):
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


class StellarContinuumAnalysisParSet(SpectralAnalysisParSet):
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
        SpectralAnalysisParSet.__init__(self, lamrange=lamrange, snr_threshold=snr_threshold,
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


class EmissionLineAnalysisParSet(SpectralAnalysisParSet):
    """
    Derived class with parameters used in the emission-line fitting.

    Should the template library key be included?  Design thing... Should
    all the SpectralAnalysisParSet objects be independent?

    allow for moments, LSF, etc

    """
    def __init__(lamrange=None, snr_threshold=None, prior=None, eml_par='STANDARD'):
        in_fl = [ int, float ]
        SpectralAnalysisParSet.__init__(self, lamrange=lamrange, snr_threshold=snr_threshold,
                                        prior=prior)
        self['analysis_type'] = 'EMLINE'
        self.add('eml_par', eml_par, dtype='str')       # Emission-line database keyword


class SpectralIndexAnalysisParSet(SpectralAnalysisParSet):
    """
    Derived class with parameters used in the spectral-index measurements.
    """
    def __init__(self, lamrange=None, snr_threshold=None, prior=None, tpl_lib='MILES',
                 eml_par='STANDARD', sindx_par = 'LICK'):

        in_fl = [ int, float ]
        SpectralAnalysisParSet.__init__(self, lamrange=lamrange, snr_threshold=snr_threshold,
                                        prior=prior)
        self['analysis_type'] = 'SPINDX'
        self.add('tpl_lib', tpl_lib, dtype='str')       # Template library keyword
        self.add('eml_par', eml_par, dtype='str')       # Emission-line database keyword
        self.add('sindx_par', sindx_par, dtype='str')   # Spectral-index database keyword


