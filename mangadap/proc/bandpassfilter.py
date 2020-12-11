# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Container class that defines a bandpass filter.

Class usage examples
--------------------

To define an bandpass filter::

    from mangadap.par.bandpassfilter import BandPassFilterPar
    p = BandPassFilterPar(index=44, name='Ha', blueside=[6483.0,6513.0],
                          redside=[6623.0,6653.0], restwave=6564.632,
                          primary=[6557.6,6571.6]) 

However, this class is mostly to provide a base class used by
:class:`mangadap.par.emissionmomentsdb.EmissionMomentsDB`,
:class:`mangadap.par.absorptionindexdb.AbsorptionIndexDB`, and
:class:`mangadap.par.bandheadindexdb.BandheadIndexDB`; it is not
really meant to be used as given above.

----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

import warnings

from IPython import embed

import numpy

from ..par.parset import KeywordParSet
from ..util.sampling import grid_borders

from matplotlib import pyplot

# Add strict versioning
# from distutils.version import StrictVersion

class BandPassFilterPar(KeywordParSet):
    r"""
    Parameter object that defines a set of band-pass filters.  This is a
    base class for similar objects used calculating fluxes and flux
    moments of emission lines
    (:class:`mangadap.par.emissionmomentsdb.EmissionMomentsDB`),
    absorption line indices
    (:class:`mangadap.par.absorptionindexdb.AbsorptionIndexDB`), and
    spectral continuum bandheads
    (:class:`mangadap.par.bandheadindexdb.BandheadIndexDB`) in the DAP.

    All wavelengths are expected to be **in vacuum**, to match the
    expected application of this filter to SDSS MaNGA spectra.

    See :class:`mangadap.par.parset.ParSet` for attributes and raised
    exceptions.

    The defined parameters are:

    .. include:: ../tables/bandpassfilterpar.rst
    """
    def __init__(self, index=None, name=None, blueside=None, redside=None, restwave=None,
                 primary=None, units=None, component=None, integrand=None, order=None):
        
        in_fl = [ int, float ]
        ar_like = [ numpy.ndarray, list ]
        n = None if name is None else name.strip()
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

        descr = ['An index used to refer to the bandpass filter',
                 'A name for the bandpass filter',
                 'A two-element vector with the starting and ending wavelength (angstroms in ' \
                    '**vacuum**) for a passband to the blue of the primary spectral feature.',
                 'A two-element vector with the starting and ending wavelength (angstroms in ' \
                    '**vacuum**) for a passband to the red of the primary spectral feature.',
                 'The rest wavelength of the line in the primary passband in angstroms **in ' \
                    'vacuum**.  This is used to convert the first and second moments to ' \
                    'velocity (km/s) for emission lines.',
                 'A two-element vector with the starting and ending wavelength (angstroms in ' \
                    '**vacuum**) for the primary passband surrounding the spectral feature of ' \
                    'interest.  This is used by ' \
                    ':class:`mangadap.par.emissionmomentsdb.EmissionMomentsDB` and ' \
                    ':class:`mangadap.par.absorptionindexdb.AbsorptionIndexDB`, but it is ' \
                    'irrelevant for :class:`mangadap.par.bandheadindexdb.BandheadIndexDB`.',
                 'Define the unit for the spectral index as either angstroms (``ang``) or ' \
                    'magnitudes (``mag``).  Currently only used by ' \
                    ':class:`mangadap.par.absorptionindexdb.AbsorptionIndexDB`.',
                 'Flag that the bandpass definition is a component of a larger set of bandpass ' \
                    'filters.  This is currently only used by ' \
                    ':class:`mangadap.par.absorptionindexdb.AbsorptionIndexDB` to combine index ' \
                    'measurements into a single index.  If True, all components with the same ' \
                    '*name* are summed to form the composite index.',
                 r'Currently only used by :class:`mangadap.par.bandheadindexdb.BandheadIndexDB`.' \
                    r'  Define the integrand over the passband used to construct and index as ' \
                    r'either flux per unit frequency (``fnu``), :math:`F_\nu`, or flux per unit ' \
                    r'wavelength (``flambda``), :math:`F_\lambda`.',
                 'Currently only used by :class:`mangadap.par.bandheadindexdb.BandheadIndexDB`.' \
                    '  Define the order to use when constructing the index.  The options are ' \
                    'either a ratio of red-to-blue or blue-to-red, which are respectively ' \
                    'selected using ``r_b`` or ``b_r``.']

        super(BandPassFilterPar, self).__init__(pars, values=values, options=options,
                                                dtypes=dtypes, descr=descr)
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
        if self.data['blueside'] is not None and len(self.data['blueside']) != 2:
            raise ValueError('Blue sideband must have two and only two elements.')
        if self.data['redside'] is not None and len(self.data['redside']) != 2:
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
        return numpy.ma.median(_y)

    _x[mask] = numpy.ma.masked
    _x = _x.compressed()

    indx = numpy.array([numpy.where(numpy.logical_and(_x > p[0], _x < p[1]))[0]
                                for p in passband], dtype=object)
    nonzero = numpy.array([ len(ii) > 0 for ii in indx ])
    if not numpy.all(nonzero):
        warnings.warn('Returning empty passbands with median values of 0!')
    return numpy.array([ 0.0 if len(ii) == 0 else numpy.ma.median(_y[ii]) for ii in indx ])


def pixel_fraction_in_passband(x, passband, dx=None):
    r"""
    Determine the width of each x interval to include in one or more
    passband integrals.

    Args:
        x (array-like):
            1D array with pixel borders.
        passband (array-like):
            1D or 2D array with the passbands to integrate over. If
            more than one passband is provided, the shape of the
            input should be :math:`(N_{\rm passband},2)`.
        dx (array-like, optional):
            1D array with the size of each pixel. If not provided,
            calculated from `x`. Can be provided to speed
            calculations for multiple calls.

    Returns:
        `numpy.ndarray`_: The array with the fraction of each pixel
        within each passband. The output shape is :math:`(N_{\rm
        pixel},N_{\rm passband})`, where :math:`N_{\rm pixel}` is the
        1 less than the number of pixel borders.
    
    """
    _x = numpy.atleast_1d(x)
    if _x.ndim != 1:
        raise ValueError('Pixel coordinates must be a 1D vector.')
    _dx = numpy.diff(_x) if dx is None else numpy.atleast_1d(dx)
    if _dx.ndim != 1:
        raise ValueError('Pixel sizes must be a 1D vector.')
    if _x.size != _dx.size+1:
        raise ValueError('Input coordinates must be the pixel borders.')
    _passband = numpy.atleast_1d(passband)
    if _passband.ndim == 1:
        return ((_x[1:]-_passband[0])/_dx).clip(0,1)+((_passband[1]-_x[:-1])/_dx).clip(0,1)-1.0
    
    return ((x[1:,None]-passband[None,:,0])/dx[:,None]).clip(0,1) \
                + ((passband[None,:,1]-x[:-1,None])/dx[:,None]).clip(0,1) - 1.0


def passband_integral(x, y, passband=None, borders=False, log=False, base=10.0, quad=False):
    r"""
    Determine the integral of a discrete set of data over a passband
    accounting for fractional pixels.

    The integral is done by a simple sum:
    
    .. math::
        S(y) = \sum_i \delta_i y_i {\rm d}x_i
        
    where :math:`\delta_i` is the fraction of pixel :math:`i` within the
    integration limits as determined by
    :func:`pixel_fraction_in_passband`.  Given the errors in :math:`y`
    (:math:`\epsilon_i`), the propagated error in :math:`S` is

    .. math::

        [\epsilon S(y)]^2 = \sum_i (\delta_i \epsilon_i {\rm d}x_i)^2, 

    which can be calculated using the ``quad`` keyword. For example,
    given values to integrate, ``y``, and errors, ``e``, calculate
    the passband integral and its error as follows::

        integral = passband_integral(x, y)
        integral_error = passband_integral(x, e, quad=True)

    .. warning::

        Although the equation for the formal error propagation is
        correct, a comparison of the formal errors to errors from a
        Monte Carlo simulation show that the former are pretty poor
        estimates of the latter.

    Args:
        x (array-like):
            The 1D array of x coordinates for the integration. If
            borders=True, this array provides the borders of the
            interval over which y is measured. These can be
            non-uniform. In this case, the number of x values must be
            :math:`N+1` for :math:`N` 1D ``y`` values (see below).
            Otherwise, the `x` values are expected to be uniformly
            sampled either linearly or geometrically. In this case,
            ``x`` can be either a two element vector giving the
            (geometric) centers of the first and last sample interval
            or a vector with the :math:`N` samples. In either case,
            :func:`mangadap.util.sampling.grid_borders` is used to
            determine the borders of the sample intervals.
        y (array-like):
            The 1D or 2D array of measured values to be integrated.
            The first axis must always provide the values to be
            integrated. If those values are specific to the passband,
            ``y`` should be 2D with shape :math:`(N_y,N_{\rm
            passband})`.
        passband (array-like, optional):
            The 1D or 2D array with the set of passbands (in units of
            x) to integrate. The length of the last dimension must
            always be 2, providing the starting and ending x
            coordinate. If 2D, the array must have shape
            :math:`(N_{\rm passband},2)`. If not provided, the
            integral is just performed over the full :math:`y`
            vector. If ``y`` is 2D, the last axis of ``y`` and the
            first axis here must match.
        borders (:obj:`bool`, optional):
            The :math:`x` values provided are the border coordinates
            over which the :math:`y` values have been determined.
        log (:obj:`bool`, optional):
            The coordinates are logarithmically sampled. Ignored if
            ``borders`` is True.
        base (:obj:`float`, optional):
            The base of the logarithm used in the geometric sampling
            if ``log`` is True.
        quad (:obj:`bool`, optional):
            Perform the quadrature sum instead of the direct sum.  This
            is used for error calculations.
        
    Returns:
        :obj:`float`, `numpy.ndarray`_: The passband integral for
        each passband.

    Raises:
        ValueError:
            Raised if ``x`` is not a 1D array, if the shapes of ``x``
            and ``y`` are not appropriately matched, if ``passbands``
            has more than 2 dimensions, if the last axis of
            ``passbands`` is not 2 elements long, or if the shapes of
            ``passbands`` and ``y`` are not appropriately matched.
    """
    # Ensure that y is a numpy array, fill any masked values with 0s,
    # and check its shape
    _y = y.filled(0.0) if isinstance(y, numpy.ma.MaskedArray) else numpy.atleast_1d(y) 

    if passband is None and _y.ndim > 1:
        raise ValueError('If no passbands provided, must provide single vector to integrate.')

    ny = _y.shape[0]

    # Check the input coordinates and calculate the borders
    _x = numpy.atleast_1d(x)
    if len(_x.shape) != 1:
        raise ValueError('Coordinates must be a 1D vector.')
    if borders and len(_x) != ny+1:
        raise ValueError('Must provide N+1 borders for N samples.')
    if not borders:
        _x, dx = grid_borders(numpy.array([x[0],x[-1]]), ny, log=log, base=base)

    if passband is None:
        # No passband, so just get the full integral; y is checked to
        # be 1D above
        return numpy.sqrt(numpy.dot(numpy.square(_y), numpy.square(numpy.diff(_x)))) if quad \
                    else numpy.dot(_y, numpy.diff(_x))

    # Check the input passbands
    _passband = numpy.atleast_2d(passband)
    if _passband.ndim > 2:
        raise ValueError('Can only handle passband arrays of 1 or 2 dimensions.')
    if _passband.shape[-1] != 2:
        raise ValueError('Passband(s) must have shape N_bands x 2.')
    
    # Reshape y if necessary
    _y = numpy.atleast_2d(_y)
    if _y.shape[0] == 1:
        _y = numpy.tile(_y[0], (_passband.shape[0],1)).T

    # Check its shape
    if _y.shape[1] != _passband.shape[0]:
        raise ValueError('For 2D array input, second axis must match the number of passbands.')
        
    # Get the (fractional) interval of each pixel in each passband
    dx = numpy.diff(_x)
    inband_dx = pixel_fraction_in_passband(_x, _passband, dx)*dx[:,None]

    # Sum of y*dx, or sqrt(sum((y*dx)^2))
    integral = numpy.sqrt(numpy.sum(numpy.square(_y*inband_dx), axis=0)) if quad \
                    else numpy.sum(_y*inband_dx, axis=0)

    # Reshape the returned array if only one passband provided
    return integral if numpy.atleast_1d(passband).ndim == 2 else integral[0]


def passband_integrated_width(x, y, passband=None, borders=False, log=False, base=10.0):
    """
    Determine the integrated width of the passband, accounting for masked pixels.
    """
    unity = y.copy()
    unity[numpy.invert(numpy.ma.getmaskarray(y))] = 1.0
    return passband_integral(x, unity, passband=passband, borders=borders, log=log, base=base)


def passband_integrated_mean(x, y, err=None, passband=None, borders=False, log=False, base=10.0):
    """
    Determine the integrated mean over a (set of) passband(s).  Nominal
    errors are returned if err is provided.
    """
    # Get the passband interval, accounting for masked pixels
    interval = passband_integrated_width(x, y, passband=passband, borders=borders, log=log,
                                         base=base)
    interval = numpy.ma.MaskedArray(interval, mask=numpy.invert(numpy.absolute(interval) > 0))

    # Get the integral of y over the passband
    integral = numpy.ma.MaskedArray(passband_integral(x, y, passband=passband, borders=borders,
                                                      log=log, base=base))
    if err is None:
        # No error
        return numpy.ma.divide(integral, interval), None

    if len(y) != len(err):
        raise ValueError('Value and error have different sizes.')

    # Return the integrated mean and error
    integral_err = numpy.ma.MaskedArray(passband_integral(x, err, passband=passband,
                                                          borders=borders, log=log, base=base,
                                                          quad=True))
    return numpy.ma.divide(integral, interval), numpy.ma.divide(integral_err, interval)


def passband_weighted_mean(x, y, z, passband=None, borders=False, yerr=None, zerr=None, log=False,
                           base=10.0):
    r"""
    Determine the y-weighted mean of z over the passband; i.e.,

    .. math::

        \mu = \frac{\int_{x_1}^{x_2} y z dx}{\int_{x_1}_{x_2} y dx},
        
    where the integrals are determined by the discrete sum; i.e.,
    :math:`\mu = S(y z)/S(y)` where :math:`S(y)` is provided by
    :func:`passband_integral`.

    Because the error in the numerator and denominator of :math:`\mu`
    are correlated, they are calculated here by explicit error
    propagation:

    .. math::

        \epsilon_\mu^2 = \frac{1}{S(y)^2} (S((z-\mu)^2 \epsilon_y^2) +
        S(y^2\epsilon_z^2)).

    .. todo::
        -doc the args

    Returns:
        `numpy.ma.MaskedArray`_: Two masked arrays with the
        passband-weighted mean and its error. The error is returned
        as `None` if no errors are provided.
    """
    # Get the passband interval, accounting for masked pixels
    weighted_integral = passband_integral(x, z*y, passband=passband, borders=borders, log=log,
                                          base=base)
    integral = passband_integral(x, y, passband=passband, borders=borders, log=log, base=base)
    integral = numpy.ma.MaskedArray(integral, mask=numpy.invert(numpy.absolute(integral) > 0))

    if yerr is None and zerr is None:
        return numpy.ma.divide(weighted_integral,integral), None

    # Construct error
    weighted_integral_err = numpy.zeros_like(weighted_integral, dtype=float)

    # Get the integrand for the y errors
    if yerr is not None:
        if yerr.shape != y.shape:
            raise ValueError('Incorrect shape for y errors.')
        ye_integrand = (z[:,None]-numpy.atleast_1d(integral)[None,:])*yerr[:,None]
        weighted_integral_err += numpy.square(passband_integral(x, ye_integrand, passband=passband,
                                                                borders=borders, log=log,
                                                                base=base, quad=True))

    # Get the integrand for the z errors
    if zerr is not None:
        if zerr.shape != z.shape:
            raise ValueError('Incorrect shape for z errors.')
        ze_integrand = _zerr*y
        weighted_integral_err += numpy.square(passband_integral(x, ze_integrand, passband=passband,
                                                                borders=borders, log=log,
                                                                base=base, quad=True))
    
    # Finalize the error
    weighted_integral_err = numpy.sqrt(weighted_integral_err)

    # Return the passband weighted mean and its error
    return numpy.ma.divide(weighted_integral, integral), \
                numpy.ma.divide(weighted_integral_err, numpy.absolute(integral))


def passband_weighted_sdev(x, y, z, passband=None, borders=False, yerr=None, zerr=None, log=False,
                           base=10.0):
    r"""
    Determine the y-weighted standard deviation of z over the passband;
    i.e.,

    .. math::

        \sigma^2 = \frac{\int_{x_1}^{x_2} y z^2 dx}{\int_{x_1}_{x_2} y
        dx} - \mu^2,
        
    where :math:`\mu` is the y-weighted mean of z over the passband (see
    :func:`passband_weighted_mean`) and the integrals are determined by
    the discrete sum; i.e., :math:`\sigma^2 = S(y z^2)/S(y) - \mu^2`
    where :math:`S(y)` is provided by :func:`passband_integral`.

    Because the error in the various components of the :math:`\sigma`
    calculation are correlated, they are calculated here by explicit
    error propagation:

    .. math::

        \epsilon_{\sigma^2}^2 = \frac{1}{S(y)^2} (S([(z-\mu)^2 + \mu^2 -
        \mu_2]^2 \epsilon_y^2) + S(4 (z-\mu)^2 y^2\epsilon_z^2)),
        
    where :math:`\mu_2` is the first term in the :math:`\sigma^2`
    calculation above and
    
    .. math::

        \epsilon_\sigma = \frac{\epsilon_{\sigma^2}}{2 \sigma}.

    Returns:
        `numpy.ma.MaskedArray`_: Two masked arrays with the
        passband-weighted standard deviation and its error. The error
        is returned as `None` if no errors are provided.
    """
    mu = passband_weighted_mean(x, y, z, passband=passband, borders=borders, log=log,
                                base=base)[0]
    mu2 = passband_weighted_mean(x, y, z*z, passband=passband, borders=borders, log=log,
                                 base=base)[0]
    sigma = numpy.ma.sqrt(mu2 - numpy.square(mu))

    if yerr is None and zerr is None:
        return sigma, None

    # Construct error
    sigma_error = numpy.zeros_like(sigma, dtype=float) if isinstance(sigma, numpy.ndarray) else 0.

    # Get the integrand for the y errors
    if yerr is not None:
        if yerr.shape != y.shape:
            raise ValueError('Incorrect shape for y errors.')
        ye_integrand = (numpy.square(z[:,None]-numpy.atleast_1d(mu)[None,:])
                            + numpy.atleast_1d(numpy.square(mu)-mu2)[None,:])*yerr[:,None]
        sigma_error += numpy.square(passband_integral(x, ye_integrand, passband=passband,
                                                      borders=borders, log=log, base=base,
                                                      quad=True))

    # Get the integrand for the z errors
    if zerr is not None:
        if zerr.shape != z.shape:
            raise ValueError('Incorrect shape for z errors.')
        ze_integrand = 2*(z[:,None]-numpy.atleast_1d(mu)[None,:])*(y*zerr)[:,None]
        sigma_error += numpy.square(passband_integral(x, ze_integrand, passband=passband,
                                                      borders=borders, log=log, base=base,
                                                      quad=True))

    # Finalize the error in sigma
    sigma_error = numpy.sqrt(sigma_error) \
                        / numpy.absolute(passband_integral(x, y, passband=passband,
                                         borders=borders, log=log, base=base) * 2*sigma)

    return sigma, sigma_error


def pseudocontinuum(x, y, passband=None, err=None, log=True, weighted_center=False):
    """
    Get the pseudocontinua in a set of passbands for a single vector (y)

    Returns:
        :obj:`float`, `numpy.ndarray`_: Return five arrays or floats:
        (1) The center of each passband, (2) the mean continuum
        level, (3) the propagated error in the continuum level (will
        be None if no errors are provided), (4) flag that part of the
        passband was masked, (5) flag that the passband was fully
        masked or empty.
    """
    # Calculate the pseudo-continua in the sidebands
    continuum, continuum_error = passband_integrated_mean(x, y, err=err, passband=passband,
                                                          log=log)
    if weighted_center:
        band_center = passband_weighted_mean(x, y, x, passband=passband, log=log)[0]
    else:
        band_center = numpy.mean(x) if passband is None else numpy.mean(passband, axis=1)

    # Calculate the fraction of the band that is covered by unmasked
    # pixels
    interval_frac = passband_integrated_width(x, y, passband=passband, log=log) \
                            / numpy.diff(passband, axis=1).ravel()

    return band_center, continuum, continuum_error, interval_frac < 1.0, \
                numpy.invert(interval_frac > 0.0)


def emission_line_equivalent_width(wave, flux, bluebands, redbands, line_centroid, line_flux,
                                   ivar=None, mask=None, log=True, redshift=None,
                                   line_flux_err=None, include_band=None):
    """
    Compute the equivalent width for emission lines provided the
    spectra and the previously measured line flux.

    Args:
        wave (`numpy.ndarray`_):
            Vector with the observed wavelength of each spectrum in
            angstroms.
        flux (`numpy.ndarray`_):
            Array (1 or 2) with the observed flux density (units in per
            angstrom) with size Nspec x Nwave.
        blueside (`numpy.ndarray`_):
            Wavelength limits for the blue sidebands in angstroms, with
            size Nbands x 2.
        redside (`numpy.ndarray`_):
            Wavelength limits for the red sidebands in angstroms, with
            size Nbands x 2.
        line_centroid (`numpy.ndarray`_):
            Wavelengths at which to sample the continuum for the
            equivalent width measurement.  Can be anything, but should
            typically be the *observed* wavelength of line center in
            angstroms, with size Nspec x Nband.
        line_flux (`numpy.ndarray`_):
            Integrated flux of the emission feature with size Nspec x
            Nband.
        ivar (`numpy.ndarray`_, optional): 
            Inverse variance in the observed flux, with shape that
            matches flux.  Default ignores error calculation.
            **Currently the equivalent width errors do not account for
            errors in the continuum characterization beneath the line.
            So this ivar array is ignored!**
        mask (`numpy.ndarray`_, optional): 
            Boolean bad-pixel mask: True values are considered to be bad
            pixels.  Default assumes all pixels are good.
        log (:obj:`bool`, optional): 
            Boolean that the spectra are logarithmically sampled in
            wavelength.
        redshift (`numpy.ndarray`_, optional):
            Redshift of each spectrum used to appropriately shift the
            definition of the bandpasses.  If a single vector is
            provided, the length must match the number of provided
            spectra; if an array is provided, the shape must be
            (Nspec,Nband).  If None, all measurements are done assuming
            the redshift is 0 (i.e., that the observed and rest frames
            are identical).
        line_flux_err (`numpy.ndarray`_, optional):
            Errors in the line flux, with size Nspec x Nband.  Default
            is to ignore the error propagation.
        include_band (`numpy.ndarray`_, optional): 
            Boolean array with size Nspec x Nband used to select which
            bands to use in the calculation for each spectrum.  Default
            is to include all bands for all spectra. 

    Returns:
        :obj:`tuple`: Returns the following arrays, all with shape
        (Nspec, Nbands):
            
            #. the passband median in the blue and red sidebands
            #. a boolean array if the sideband medians were measured
               and have positive values
            #. The continuum value used to calculate the equivalent
               width
            #. The equivalent with measurements and their errors; the
               errors are 0 if no line flux errors were provided

    Raises:
        ValueError: Raised if the shapes of the input arrays are not
            correct.

    """
    # Errors in the flux are currently not
    if ivar is not None:
        warnings.warn('Equivalent width does not propagate error from continuum measurement.')
    noise = None

    # Convert the flux to a masked array, if necessary
    if mask is not None and mask.shape != flux.shape:
        raise ValueError('Input mask must have the same shape as the flux array.')
    _flux = numpy.ma.atleast_2d(flux if isinstance(flux, numpy.ma.MaskedArray) \
                                        else numpy.ma.MaskedArray(flux, mask=mask))
    if len(wave) != _flux.shape[1]:
        raise ValueError('Wavelength vector does not match shape of the flux array.')

#    # Convert the ivar to 1-sigma error, if available
#    if ivar is None:
#        noise = None
#    else:
#        if ivar.shape != _flux.shape:
#            raise ValueError('Input ivar array must be the same shape as the flux array.')
#        _ivar = ivar if isinstance(ivar, numpy.ma.MaskedArray) else \
#                    numpy.ma.MaskedArray(ivar, mask=numpy.ma.getmaskarray(_flux))
#        noise = numpy.ma.sqrt(1.0 /_ivar)

    # Check band definitions
    if bluebands.shape[1] != 2:
        raise ValueError('Input bands must have shape Nband x 2.')
    if bluebands.shape != redbands.shape:
        raise ValueError('Input side bands do not have the same shape.')

    # Check the shapes
    nspec = flux.shape[0]
    nbands = bluebands.shape[0]

    expected_shape = (nspec, nbands)
    if line_centroid.shape != expected_shape:
        raise ValueError('Line centroid array must have shape: {0}'.format(expected_shape))
    if line_flux.shape != expected_shape:
        raise ValueError('Line flux array must have shape: {0}'.format(expected_shape))
    if line_flux_err is not None and line_flux_err.shape != expected_shape:
        raise ValueError('Line flux error array must have shape: {0}'.format(expected_shape))

    # Check the input band flags
    _include_band = numpy.ones(expected_shape, dtype=numpy.bool) \
                        if include_band is None else include_band
    if _include_band.shape != expected_shape:
        raise ValueError('Bands flags array must have shape: {0}'.format(expected_shape))

    # Check the input redshifts
    _redshift = numpy.zeros(expected_shape, dtype=numpy.float) if redshift is None else redshift
    if _redshift.ndim == 1:
        if len(_redshift) != nspec:
            raise ValueError('Must provide at least one redshift per input spectrum.')
        _redshift = numpy.array([_redshift]*nbands).T
    if _redshift.ndim == 2 and _redshift.shape != expected_shape:
        raise ValueError('Provided redshift array does not match the number of spectra and bands.')

    # Initialize the output data
    bmed = numpy.zeros(expected_shape, dtype=float)
    rmed = numpy.zeros(expected_shape, dtype=float)
    pos = _include_band.copy()
    ewcont = numpy.zeros(expected_shape, dtype=float)
    ew = numpy.zeros(expected_shape, dtype=float)
    ewerr = numpy.zeros(expected_shape, dtype=float)

    # Measure the pseudo-continuum in the sidebands
    for i in range(nspec):

        print('Measuring emission-line equivalent widths in spectrum: {0}/{1}'.format(i+1,nspec),
              end='\r')

        # Spectra to use
        spec = _flux[i,:]
        _noise = None if noise is None else noise[i,:]

        # Bands to use
        indx = _include_band[i,:]

        # Shift the bands to the appropriate redshift
        _bluebands = bluebands[indx,:]*(1.0+_redshift[i,indx,None])
        _redbands = redbands[indx,:]*(1.0+_redshift[i,indx,None])

        # Center of each band (for all bands!)
        bcen = numpy.mean(bluebands*(1.0+_redshift[i,:,None]), axis=1)
        rcen = numpy.mean(redbands*(1.0+_redshift[i,:,None]), axis=1)

        # Median of each band
        bmed[i,indx] = passband_median(wave, _flux[i,:], passband=_bluebands)
        rmed[i,indx] = passband_median(wave, _flux[i,:], passband=_redbands)

        # Construct the continuum level at the center of the main
        # band
        pos[i,:] = indx & (bmed[i,:] > 0) & (rmed[i,:] > 0)

        # Continuum at line centroid
        m = (rmed[i,pos[i,:]] - bmed[i,pos[i,:]])/(rcen[pos[i,:]] - bcen[pos[i,:]])
        ewcont[i,pos[i,:]] = m*(line_centroid[i,pos[i,:]] - bcen[pos[i,:]]) + bmed[i,pos[i,:]]

        # Compute the equivalent width
        ew[i,pos[i,:]] = line_flux[i,pos[i,:]] / ewcont[i,pos[i,:]]
        ewerr[i,pos[i,:]] = 0.0 if line_flux_err is None \
                                else line_flux_err[i,pos[i,:]] / ewcont[i,pos[i,:]]

        # Correct for the redshift
        ew[i,pos[i,:]] /= (1.0+_redshift[i,pos[i,:]])
        ewerr[i,pos[i,:]] /= (1.0+_redshift[i,pos[i,:]])

#        pyplot.step(wave, _flux[i,:], where='mid', color='k', lw=0.5, zorder=1)
#        pyplot.scatter(bcen[indx], bmed[i,indx], marker='.', s=50, color='b', lw=0, zorder=3)
#        pyplot.scatter(rcen[indx], rmed[i,indx], marker='.', s=50, color='r', lw=0, zorder=3)
#        pyplot.scatter(line_centroid[i,indx], ewcont[i,indx], marker='.', s=50, color='g',
#                       lw=0, zorder=3)
#        pyplot.show()

    print('Measuring emission-line equivalent widths in spectrum: {0}/{0}'.format(nspec))
    return bmed, rmed, pos, ewcont, ew, ewerr


