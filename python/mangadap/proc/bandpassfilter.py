# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Container class that defines a bandpass filter.

*License*:
    Copyright (c) 2015, SDSS-IV/MaNGA Pipeline Group
        Licensed under BSD 3-clause license - see LICENSE.rst

*Source location*:
    $MANGADAP_DIR/python/mangadap/par/bandpassfilter.py

*Imports and python version compliance*:
    ::

        from __future__ import division
        from __future__ import print_function
        from __future__ import absolute_import
        from __future__ import unicode_literals
    
        import sys
        if sys.version > '3':
            long = int

        import warnings
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

*Revision history*:
    | **18 Mar 2016**: Original implementation by K. Westfall (KBW)
    | **20 Apr 2016**: (KBW) Include measurements

"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import sys
if sys.version > '3':
    long = int

import warnings
import numpy
from ..par.parset import ParSet
from ..util.instrument import _pixel_borders

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
        Check the parameter list:
            
            - Make sure arrays only have two elements.

        .. todo::
            - Add check to __setitem__()?

        Raises:
            ValueError: Raised if one of the conditions above are not
                met.
        """
        if self.data['primary'] is not None and len(self.data['primary']) != 2:
            raise ValueError('Primary passband must have two and only two elements.')
        if len(self.data['blueside']) != 2:
            raise ValueError('Blue sideband must have two and only two elements.')
        if len(self.data['redside']) != 2:
            raise ValueError('Red sideband must have two and only two elements.')


def passband_median(x, y, passband=None):
    """
    Determine the median of the y values within the passband
    """

    # Ensure that y is a numpy array and compress it if it's a masked
    # array
    _x = x if isinstance(x, numpy.ma.MaskedArray) else numpy.ma.MaskedArray(x)
    _y = y if isinstance(y, numpy.ma.MaskedArray) else numpy.ma.MaskedArray(y)
    if len(_x.shape) != 1 or len(_y.shape) != 1:
        raise ValueError('Can only provide one-dimensional vectors.')
    mask = numpy.ma.getmaskarray(_x) | numpy.ma.getmaskarray(_y)
    _y[mask] = numpy.ma.masked
    _y = _y.compressed()
    if _y.size == 0:
        warnings.warn('No elements remain in the passband.')
        return 0.0

    if passband is None:
        return numpy.median(_y)

    _x[mask] = numpy.ma.masked
    _x = _x.compressed()

    indx = numpy.array([ numpy.arange(_x.size)[numpy.logical_and(_x > p[0], _x < p[1])]
                                for p in passband ])
    nonzero = numpy.array([ len(ii) > 0 for ii in indx ])
    if numpy.any(nonzero):
        warnings.warn('Returning empty passbands with median values of 0!')
    return numpy.array([ 0.0 if len(ii) == 0 else numpy.median(_y[ii]) for ii in indx ])


def pixel_fraction_in_passband(x, passband, dx):
    """
    Get the width of each x interval to include in the passband
        integral
    """
    return ((x[1:]-passband[0])/dx).clip(0,1) + ((passband[1]-x[:-1])/dx).clip(0,1) - 1.0


def passband_integral(x, y, passband=None, borders=False, log=False, base=10.0):
    """
    Determine the integral of a discrete set of data over a passband
    accounting for fractional pixels.

    If borders = True, the x values are the borders of the interval over
    which y is measured.  These can be non-uniform.  In this case, the
    number of x values must be N+1 for N y values.
    
    Otherwise, the x values are expected to be uniformly sampled either
    linearly or geometrically.  In this case, x can be either a two
    element vector giving the (geometric) centers of the first and last
    sample interval or a vector with the N samples.  In either case,
    :func:`mangadap.util.instrument._pixel_borders` is used to determine
    the borders of the sample intervals.

    If the passband is None, the integral is determined across all data.
    """
    # Ensure that y is a numpy array, fill any masked values with 0s,
    # and check its shape
    _y = y.filled(0.0) if isinstance(y, numpy.ma.MaskedArray) else numpy.asarray(y) 
    _x = numpy.asarray(x)
    if len(_x.shape) != 1 or len(_y.shape) != 1:
        raise ValueError('Can only provide one-dimensional vectors.')
    # Check the input and calculate the pixel borders
    if borders and len(_x) != len(_y)+1:
        raise ValueError('Must provide N+1 borders for N samples.')
    if not borders:
        _x, dx = _pixel_borders( numpy.array([x[0],x[-1]]), _y.size, log=log, base=base)
    if passband is None:
        return numpy.dot(_y, numpy.diff(_x))
    _passband = numpy.asarray(passband)
    if len(_passband.shape) > 2:
        raise ValueError('Can only handle passband arrays of 1 or 2 dimensions.')
    if _passband.shape[-1] != 2:
        raise ValueError('Passband(s) must have shape N_bands x 2.')

    # Interval for each sample
    dx = numpy.diff(_x)

    # Get the integral over a single passband
    if len(_passband.shape) == 1:
        inband_dx = pixel_fraction_in_passband(_x, _passband, dx)*dx
    else:
        # Get the integral over all passbands
        inband_dx = numpy.empty((_passband.shape[0],y.size), dtype=numpy.float)
        for i,pb in enumerate(_passband):
            inband_dx[i,:] = pixel_fraction_in_passband(_x, pb, dx)*dx
    # Sum of y*dx
    return numpy.dot(_y, inband_dx.T)


def passband_integrated_width(x, y, passband=None, log=False, base=10.0):
    """
    Determine the integrated width of the passband, accounting for masked pixels.
    """
    unity = y.copy()
    unity[numpy.invert(numpy.ma.getmaskarray(y))] = 1.0
    return passband_integral(x, unity, passband=passband, log=log, base=base)


def passband_integrated_mean(x, y, passband=None, err=None, log=False, base=10.0):
    """
    Determine the integrated mean over a (set of) passband(s).  Nominal
    errors are returned if err is provided.
    """
    # Get the passband interval, accounting for masked pixels
    interval = passband_integrated_width(x, y, passband=passband, log=log, base=base)
    interval = numpy.ma.MaskedArray(interval, mask=numpy.invert(numpy.absolute(interval) > 0))

    # Get the integral of y over the passband
    integral = numpy.ma.MaskedArray(passband_integral(x, y, passband=passband, log=log, base=base))
    
    # Return the integrated mean, and the error if possible
    if err is None:
        return integral/interval, None
    err_integral = numpy.ma.MaskedArray(passband_integral(x, numpy.square(err), passband=passband,
                                        log=log, base=base))
    return integral/interval, numpy.ma.sqrt(err_integral)/interval


def passband_weighted_mean(x, y, z, passband=None, yerr=None, zerr=None, log=False, base=10.0):
    """
    Determine the y-weighted integral of z over the passband.  Nominal
    errors are returned if err is provided.
    """
    # Get the passband interval, accounting for masked pixels
    weighted_integral = passband_integral(x, z*y, passband=passband, log=log, base=base)
    integral = passband_integral(x, y, passband=passband, log=log, base=base)
    integral = numpy.ma.MaskedArray(integral, mask=numpy.invert(numpy.absolute(integral) > 0))

    if yerr is None and zerr is None:
        return weighted_integral/integral, None

    weighted_integral = numpy.ma.MaskedArray(weighted_integral,
                                        mask=numpy.invert(numpy.absolute(weighted_integral) > 0))

    werr = numpy.zeros(z.shape)
    if yerr is not None:
        werr += z*numpy.square(yerr)
    if zerr is not None:
        werr += y*numpy.square(zerr)
    weighted_integral_err = passband_integral(x, werr, passband=passband, log=log, base=base)

    integral_err = 0 if yerr is None else passband_integral(x, numpy.square(yerr),
                                                            passband=passband, log=log, base=base)
    weighted_mean = weighted_integral/integral
    weighted_mean_err = numpy.ma.sqrt(
                            weighted_integral_err*numpy.square(weighted_mean/weighted_integral) \
                                    + integral_err*numpy.square(weighted_mean/integral))
    return weighted_mean, weighted_mean_err

