# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""

Binning!

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
from scipy import sparse

from astropy.io import fits

from vorbin.voronoi_2d_binning import voronoi_2d_binning

from ..par.parset import KeywordParSet
from ..util.geometry import SemiMajorAxisCoo
from ..util.covariance import Covariance


# BASE CLASS -----------------------------------------------------------

class SpatialBinning():
    """
    Base class for spatially binning data.
    """
    def __init__(self, bintype, par=None):
        self.bintype = bintype
        self.par = par

    @staticmethod
    def bin_index(x, y, gpm=None, par=None):
        """
        Undefined in base class. Just provides expected calling
        sequence.
        """
        pass
# ----------------------------------------------------------------------


# GLOBAL BINNING -------------------------------------------------------
class GlobalBinning(SpatialBinning):
    """
    Class that performs the global binning.
    """
    def __init__(self):
        SpatialBinning.__init__(self, 'global')

    @staticmethod
    def bin_index(x, y, gpm=None, par=None):
        """
        Bin the data and return the indices of the bins.

        The calling sequence is just needed to match the expectation
        for a :class:`SpatialBinning` object.

        Args:
            x (`numpy.ndarray`_):
                Fiducial on-sky X position of each spectrum. X
                increases with RA. Only used to set the size of the
                returned array.
            y (`numpy.ndarray`_):
                **Not used.**
            gpm (`numpy.ndarray`_, optional):
                Good value (pixel) mask used to select data that should be
                included in the bin.  If None, all values are considered valid.
                If provided, masked values (gpm=False) are given a bin index of
                -1.
            par (object, optional):
                **Not used.**

        Returns:
            `numpy.ndarray`_: An integer bin index for each spectrum.
        """
        binid = numpy.zeros(x.shape, dtype=int)
        if gpm is not None:
            if gpm.size != x.size:
                raise ValueError(f'Provided GPM has incorrect shape; found {gpm.size}, ' \
                                 f'expected {x.size}.')
            binid[numpy.logical_not(gpm)] = -1
        return binid
# ----------------------------------------------------------------------

# RADIAL BINNING -------------------------------------------------------
class RadialBinningPar(KeywordParSet):
    """
    Class with parameters used by the radial binning algorithm.  See
    :class:`mangadap.par.parset.ParSet` for attributes.

    The defined parameters are:

    .. include:: ../tables/radialbinningpar.rst
    """
    def __init__(self, center=[0., 0.], pa=-1, ell=-1, radius_scale=-1, radii=[-1,-1,10],
                 log_step=False):
        in_fl = [ int, float ]
        ar_like = [ numpy.ndarray, list ]
        
        pars =   [ 'center',  'pa', 'ell', 'radius_scale', 'radii', 'log_step' ]
        values = [   center,    pa,   ell,   radius_scale,   radii,   log_step ]
        dtypes = [  ar_like, in_fl, in_fl,          in_fl, ar_like,       bool ]

        descr = ['A two-element array defining the center to use in the definition of the ' \
                    'elliptical bins.  This is defined as a sky-right offset in arcseconds from ' \
                    'the nominal center of the object.',
                 'Sets the position angle, defined from N through E of the major axis of the ' \
                    'isophotal ellipse used to define the elliptical bins.  Set to <0 to ' \
                    'indicate that the value should be filled by the metadata in the reduction ' \
                    'assessment class; see :func:`fill`.',
                 'Sets the ellipticity (1-b/a) of the isophotal ellipse use to define the ' \
                    'elliptical bins.  Set to <0 to indicate that the value should be filled by ' \
                    'the metadata in the reduction assessment class; see :func:`fill`.',
                 'Defines a scale factor to use when defining the radial bins.  For example, ' \
                    'you might want to scale to the a certain number of effective radii or ' \
                    'physical scale in kpc.  For no scale, use 1.0.  Set to <0 to indicate that ' \
                    'the value should be filled by the metadata in the reduction assessment ' \
                    'class; see :func:`fill`.',
                 'A three-element array defining the starting and ending radius for the bin ' \
                    'edges and the number of bins to create.  If the starting radius is -1, ' \
                    'the inner-most radius is set to 0 when not using log bins or 0.1 arcsec ' \
                    'when using logarithmically spaced bins.  If the ending radius is -1, the ' \
                    'outer-most radius is set by the spaxel at the largest radius.',
                 'A flag that the radial bins should be a geometric series.']

        super().__init__(pars, values=values, dtypes=dtypes, descr=descr)
    
    def toheader(self, hdr):
        """
        Copy some of the parameters to a header.
        Args:
            hdr (`astropy.io.fits.Header`_):
                Header object to write to.

        Returns:
            `astropy.io.fits.Header`_: Edited header object

        Raises:
            TypeError:
                Raised if input is not an `astropy.io.fits.Header`_
                object.
        """
        if not isinstance(hdr, fits.Header):
            raise TypeError('Input is not a astropy.io.fits.Header object!')

        hdr['BINCX'] = (self['center'][0], 'Radial binning on-sky X center')
        hdr['BINCY'] = (self['center'][1], 'Radial binning on-sky Y center')
        hdr['BINPA'] = (self['pa'], 'Radial binning position angle')
        hdr['BINELL'] = (self['ell'], 'Radial binning ellipticity (1-b/a)')
        hdr['BINSCL'] = (self['radius_scale'], 'Radial binning radius scale (arcsec)')
        hdr['BINRAD'] = (str(list(self['radii'])), 'Radial binning radius sampling')
        hdr['BINLGR'] = (str(self['log_step']), 'Geometric step used by radial binning')
        return hdr

    def fromheader(self, hdr):
        """
        Copy the information from the header.

        Args:
            hdr (`astropy.io.fits.Header`_):
                Header object to read from.
        """
        if not isinstance(hdr, fits.Header):
            raise TypeError('Input is not a astropy.io.fits.Header object!')

        self['center'] = [ hdr['BINCX'], hdr['BINCY'] ]
        self['pa'] = hdr['BINPA']
        self['ell'] = hdr['BINELL']
        self['radius_scale'] = hdr['BINSCL']
        self['radii'] = eval(hdr['BINRAD'])
        self['log_step'] = bool(hdr['BINLGR'])

    # TODO: This should actually use the cube metadata, not the rdxqa
    def fill(self, rdxqa):
        """
        Use the metadata in the reduction assessment results to fill out the
        binning parameters.
        """
        if self['pa'] < 0:
            self['pa'] = 0. if rdxqa.pa is None else rdxqa.pa
        if self['ell'] < 0:
            self['ell'] = 0. if rdxqa.ell is None else rdxqa.ell
        if self['radius_scale'] < 0:
            self['radius_scale'] = 1.0 if rdxqa.reff is None else rdxqa.reff


class RadialBinning(SpatialBinning):
    """
    Class that performs the radial binning.
    """
    def __init__(self, par=None):
        SpatialBinning.__init__(self, 'radial', par=par)
        self.sma_coo = None
        self.to_center = False
        self.rs = None
        self.dr = None

    def _r_start_step(self, r):
        """
        Bin starting radii are::
            
            rs = self.par['radii'][0]
            re = self.par['radii'][1]
            nr = int(self.par['radii'][2])
            starting_radii = numpy.logspace(numpy.log10(rs), numpy.log10(re), num=nr,
                                            endpoint=False) if self.par['log_step'] else \
                             numpy.linspace(rs, re, num=nr, endpoint=False)

        For ending radii, just swap ``rs`` and ``re`` above.

        Args:
            r (`numpy.ndarray`_):
                Fiducial radius of each spectrum.
        """
        # Set minimum radius, if necessary
        if self.par['radii'][0] < 0:
            ndec = int(numpy.floor(numpy.absolute(numpy.log10(0.5/self.par['radius_scale']))))
            ndec = 2 if ndec < 2 else ndec
            self.par['radii'][0] = numpy.around(0.5/self.par['radius_scale'], decimals=ndec) \
                                        if self.par['log_step'] else 0.0
            self.to_center = True

        # Set maximum radius, if necessary
        if self.par['radii'][1] < 0:
            self.par['radii'][1] = numpy.around(numpy.amax(r)+0.1, decimals=2)

        # Get the logarithmically stepped values
        if self.par['log_step']:
            self.rs = numpy.log10(self.par['radii'][0])
            self.dr = (numpy.log10(self.par['radii'][1]) - self.rs)/self.par['radii'][2]
            return

        # Get the linearly stepped values
        self.rs = self.par['radii'][0]
        self.dr = (self.par['radii'][1] - self.rs)/self.par['radii'][2]
        
    def bin_index(self, x, y, gpm=None, par=None):
        """
        Bin the data and return the indices of the bins.

        Args:
            x (`numpy.ndarray`_):
                Fiducial on-sky X position of each spectrum. X
                increases with RA.
            y (`numpy.ndarray`_):
                Fiducial on-sky Y position of each spectrum. Y
                increases with DEC.
            gpm (`numpy.ndarray`_, optional):
                Good value (pixel) mask used to select data that should be
                included in the bin.  If None, all values are considered valid.
                If provided, masked values (gpm=False) are given a bin index of
                -1.
            par (:class:`RadialBinningPar`, optional):
                Binning parameters.  Cannot be None.

        Returns:
            `numpy.ndarray`_: An integer bin index for each spectrum.

        Raises:
            ValueError:
                Raised if the sizes of ``x`` and ``y`` do not match
                or if ``par`` is None.
        """
        # Check the input
        if x.size != y.size:
            raise ValueError('Dimensionality of x and y coordinates do not match!')
        if gpm is not None and gpm.size != x.size:
            raise ValueError(f'Provided GPM has incorrect shape; found {gpm.size}, ' \
                             f'expected {x.size}.')

        # Perform any necessary initialization
        if par is not None:
            self.par = par
            if self.sma_coo is not None:
                del self.sma_coo
                self.sma_coo = None
        if self.par is None:
            raise ValueError('Required parameters for RadialBinning have not been defined.')
        if self.sma_coo is None:
            self.sma_coo = SemiMajorAxisCoo(xc=self.par['center'][0], yc=self.par['center'][1],
                                            rot=0.0, pa=self.par['pa'], ell=self.par['ell'])

        # Get the radius and azimuth
        r, theta = self.sma_coo.polar(x,y)

        # Get the starting radius and step per bin
        r /= self.par['radius_scale']
        indx = r >= 0.0
        if self.par['log_step']:
            indx = r > 0
        self._r_start_step(r[indx])

        # Determine the bin indices
        if self.par['log_step']:
            r[indx] = numpy.log10(r[indx])

        binid = numpy.floor((r - self.rs)/self.dr).astype(int)
        if self.to_center:
            binid[(binid<0) & indx] = 0          # Go all the way to R=0
        else:
            binid[r<self.rs] = -1            # or not
        binid[binid >= int(self.par['radii'][2])] = -1  # Remove any points outside the last bin

        if gpm is not None:
            # Remove any masked pixels
            binid[gpm] = -1

        return binid

    def bin_area(self):
        """
        Return the nominal area for each elliptical bin.

        Returns:
            `numpy.ndarray`_: Array with the analytic on-sky area of
            each elliptical bin in square arcsec given its radial
            limits.

        Raises:
            ValueError:
                Raised if :attr:`par` is None.
        """
        if self.par is None:
            raise ValueError('Required parameters not defined.')

        rs = self.par['radii'][0]
        re = self.par['radii'][1]
        nr = int(self.par['radii'][2])
        bin_radii = numpy.append(
                        numpy.logspace(numpy.log10(rs), numpy.log10(re), num=nr, endpoint=False) \
                        if self.par['log_step'] \
                        else numpy.linspace(rs, re, num=nr, endpoint=False), re) \
                                * self.par['radius_scale']

        return numpy.pi*(1.0-self.par['ell'])*numpy.diff(numpy.square(bin_radii))
# ----------------------------------------------------------------------


# VORONOI BINNING ------------------------------------------------------
class VoronoiBinningPar(KeywordParSet):
    r"""
    Class with parameters used by the Voronoi binning algorithm.

    See :class:`mangadap.par.parset.ParSet` for attributes.  See
    `vorbin`_ for the main algorithm.

    The defined parameters are:

    .. include:: ../tables/voronoibinningpar.rst
    """
    def __init__(self, target_snr=10., signal=None, noise=None, covar=None):
        in_fl = [ int, float ]
        covar_type = [ float, numpy.ndarray, Covariance, sparse.spmatrix ]
        ar_like = [ numpy.ndarray, list ]
        
        pars =   [ 'target_snr', 'signal', 'noise',    'covar' ]
        values = [   target_snr,   signal,   noise,      covar ]
        dtypes = [        in_fl,  ar_like, ar_like, covar_type ]

        descr = ['The target S/N for each bin.',
                 'The array of signal measurements for each on-sky position to bin.  See ' \
                    ':func:`fill` to fill this based on the data in the reduction assessments ' \
                    'object.',
                 'The array of noise measurements for each on-sky position to bin.  If not ' \
                    'provided, ``covar`` must be provided and be a full covariance matrix.  See ' \
                    ':func:`fill` to fill this based on the data in the reduction assessments ' \
                    'object.',
                 r'Covariance matrix or calibration normalization.  For the latter, the value ' \
                    r'is used to renormalize using :math:`n_{\rm calib} = n_{\rm nominal} ' \
                    r'(1 + \alpha\ \log\ N_{\rm bin})`, where :math:`N_{\rm bin}` is the number ' \
                    r'of binned spaxels and :math:`\alpha` is the value provided. See ' \
                    r':func:`fill` to fill this based on the data in the reduction assessments ' \
                    r'object.']

        super().__init__(pars, values=values, dtypes=dtypes, descr=descr)

    def toheader(self, hdr):
        """
        Copy some of the parameters to a header.

        Args:
            hdr (`astropy.io.fits.Header`_):
                Header object to write to.

        Returns:
            `astropy.io.fits.Header`_: Edited header object

        Raises:
            TypeError:
                Raised if input is not an `astropy.io.fits.Header`_
                object.
        """
        if not isinstance(hdr, fits.Header):
            raise TypeError('Input is not a astropy.io.fits.Header object!')

        hdr['BINSNR'] = (self['target_snr'], 'Voronoi binning minimum S/N')
        if isinstance(self['covar'], float):
            hdr['BINCOV'] = ('calib', 'Voronoi binning S/N covariance type')
            hdr['NCALIB'] = (self['covar'], 'Noise calibration value')
        elif isinstance(self['covar'], (numpy.ndarray, Covariance, sparse.spmatrix)):
            hdr['BINCOV'] = ('full', 'Voronoi binning S/N covariance type')
        else:
            hdr['BINCOV'] = ('none', 'Voronoi binning S/N covariance type')
        return hdr

    def fromheader(self, hdr):
        """
        Copy the information from the header.

        Args:
            hdr (`astropy.io.fits.Header`_):
                Header object to read from.
        """
        if not isinstance(hdr, fits.Header):
            raise TypeError('Input is not a astropy.io.fits.Header object!')

        self['target_snr'] = hdr['BINSNR']
        if hdr['BINCOV'] == 'calib':
            self['covar'] = hdr['NCALIB']

    def fill(self, rdxqa):
        """
        Use the data in the reduction assessment results to fill out the binning
        parameters.
        """
        self['signal'] = rdxqa['SPECTRUM'].data['SIGNAL']
        if rdxqa.correlation is None or isinstance(self['covar'], (float, numpy.floating)):
            self['noise'] = numpy.sqrt(rdxqa['SPECTRUM'].data['VARIANCE'])
            return
        # Overwrite any existing calibration coefficient
        rdxqa.correlation.revert_correlation()
        covar = rdxqa.correlation.toarray()
        rdxqa.correlation.to_correlation()
        i, j = numpy.meshgrid(numpy.arange(covar.shape[0]), numpy.arange(covar.shape[1]))
        self['covar'] = Covariance(inp=sparse.coo_matrix((covar[covar > 0].ravel(),
                                   (i[covar > 0].ravel(), j[covar > 0].ravel())),
                                   shape=covar.shape).tocsr())


class VoronoiBinning(SpatialBinning):
    """
    Class that wraps the contributed voronoi binning code.
    """
    def __init__(self, par=None):
        SpatialBinning.__init__(self, 'voronoi', par=par)
        self.covar = None

    def sn_calculation_no_covariance(self, index, signal, noise): 
        """
        S/N calculation built to match the requirements of
        :func:`mangadap.contrib.voronoi_2d_binning`.

        This calculation is identical to the one provided directly by
        the code.

        Args:
            index (`numpy.ndarray`_):
                Indices of the measurements in a single bin.
            signal (`numpy.ndarray`_):
                The signal measurements.
            noise (`numpy.ndarray`_):
                The noise measurements.

        Returns:
            :obj:`float`: The nominal signal-to-noise reached by
            summing the measurements selected by ``index``.
        """
        return  numpy.sum(signal[index]) / numpy.sqrt(numpy.sum(numpy.square(noise[index])))

    def sn_calculation_covariance_matrix(self, index, signal, noise):
        """
        Calculate the S/N using a full covariance matrix.

        The method uses the internal :attr:`covar`.

        Args:
            index (`numpy.ndarray`_):
                Indices of the measurements in a single bin.
            signal (`numpy.ndarray`_):
                The signal measurements.
            noise (`numpy.ndarray`_):
                The noise measurements.

        Returns:
            :obj:`float`: The nominal signal-to-noise reached by
            summing the measurements selected by ``index``, including
            any covariance.
        """
        _index = numpy.atleast_1d(index)
        return numpy.sum(signal[index])/numpy.sqrt(numpy.sum(self.covar[_index,:][:,_index]))

    def sn_calculation_calibrate_noise(self, index, signal, noise):
        r"""
        Calculate the S/N using a calibration of the S/N following:

        .. math::

            N_{\rm calib} = N_{\rm nominal} (1 + \alpha\ \log N_{\rm bin})

        where :math:`N_{\rm bin}` is the number of binned spaxels and
        :math:`\alpha` is an empirically derived constant used to adjust
        the noise for the affects of covariance.

        The calibration constant, :math:`\alpha`, is kept internally
        using :attr:`covar`.

        Args:
            index (:obj:`int`, `numpy.ndarray`_):
                Indices of the measurements in a single bin.
            signal (`numpy.ndarray`_):
                The signal measurements.
            noise (`numpy.ndarray`_):
                The noise measurements.

        Returns:
            :obj:`float`: The nominal signal-to-noise reached by
            summing the measurements selected by ``index``, after
            applying the covariance calibration.
        """
        nbin = 1 if isinstance(index, (int, numpy.integer)) else len(index)
        return numpy.sum(signal[index]) / (numpy.sqrt(numpy.sum(numpy.square(noise[index]))) \
                        * (1.0 + self.covar*numpy.log10(nbin)))

    def bin_index(self, x, y, gpm=None, par=None):
        """
        Bin the data and return the indices of the bins.

        Args:
            x (`numpy.ndarray`_):
                Fiducial on-sky X position of each spectrum. X
                increases with RA.
            y (`numpy.ndarray`_):
                Fiducial on-sky Y position of each spectrum. Y
                increases with DEC.
            gpm (`numpy.ndarray`_, optional):
                Good value (pixel) mask used to select data that should be
                included in the bin.  If None, all values are considered valid.
                If provided, masked values (gpm=False) are given a bin index of
                -1.
            par (:class:`RadialBinningPar`, optional):
                Binning parameters.  Cannot be None.

        Returns:
            `numpy.ndarray`_: An integer bin index for each spectrum.

        Raises:
            ValueError:
                Raised if the sizes of ``x`` and ``y`` do not match,
                if ``par`` is None, if various checks of the signal,
                noise, or covariance elements are incorrectly
                matched.
        """
        # Check the position input
        if x.size != y.size:
            raise ValueError('Dimensionality of x and y coordinates do not match!')
        if gpm is not None and gpm.size != x.size:
            raise ValueError(f'Provided GPM has incorrect shape; found {gpm.size}, ' \
                             f'expected {x.size}.')

        # Initialize the parameter object
        if par is not None:
            self.par = par
        if self.par is None:
            raise ValueError('Required parameters for Voronoi binning have not been defined.')

        # Check the signal input
        if self.par['signal'] is None:
            raise ValueError('Signal measurements not provided for Voronoi binning!')
        if self.par['signal'].size != x.size:
            raise ValueError('Dimensionality of signal does not match on-sky coordinates.')

        # Construct the covariance input if provided
        self.covar = None
        sn_func = self.sn_calculation_no_covariance
        if self.par['covar'] is not None:
            if isinstance(self.par['covar'], Covariance):
                # Check dimensionality
                if self.par['covar'].dim == 3:
                    raise ValueError('Input Covariance object must be two-dimensional')
                # Make sure the object is a full covariance matrix
                if self.par['covar'].is_correlation:
                    self.par['covar'].revert_correlation()
            # Fill the full array if necessary
            if not isinstance(self.par['covar'], (numpy.ndarray,float)):
                self.covar = self.par['covar'].toarray()
                sn_func = self.sn_calculation_covariance_matrix
            else:
                self.covar = self.par['covar']
                sn_func = self.sn_calculation_calibrate_noise

            if isinstance(self.covar, numpy.ndarray):
                # Check that the size matches the input signal array
                if self.covar.shape[0] != self.par['signal'].size:
                    raise ValueError('Dimensionality of covariance does not match signal.')
                _noise = numpy.sqrt(numpy.diag(self.covar))

        # Make sure the noise vector is available
        if self.covar is None or isinstance(self.covar, float):
            _noise = self.par['noise']
        if _noise is None:
            raise ValueError('Could not construct noise measurements for Voronoi binning.')
        if _noise.size != self.par['signal'].size:
            raise ValueError('Dimensionality of noise does not match signal.')

        # Down-select to the valid spaxels
        revert_covar = False
        if gpm is None:
            _x = x
            _y = y
            _signal = self.par['signal']
        else:
            _x = x[gpm]
            _y = y[gpm]
            _signal = self.par['signal'][gpm]
            _noise = _noise[gpm]
            if self.covar is not None and isinstance(self.covar, numpy.ndarray):
                revert_covar = True
                _covar = self.covar.copy()
                self.covar = self.covar[gpm,:][:,gpm]

        # All spaxels have S/N greater than threshold, so return each
        # spaxel in its own "bin"
        if numpy.min(_signal/_noise) > self.par['target_snr']:
            warnings.warn('All pixels have enough S/N. Binning is not needed')
            if gpm is None:
                return numpy.arange(x.size, dtype=int)
            binid = numpy.full(x.size, -1, dtype=int)
            binid[gpm] = numpy.arange(numpy.sum(gpm))
            return binid

        # Cannot reach the S/N using all spaxels, so return all spaxels
        # in a single bin
        sn_total = sn_func(numpy.arange(_signal.size), _signal, _noise)
        if sn_total < self.par['target_snr']:
            warnings.warn('Cannot reach target S/N using all data.')
            if gpm is None:
                return numpy.zeros(x.size, dtype=int)
            binid = numpy.full(x.size, -1, dtype=int)
            binid[gpm] = numpy.zeros(numpy.sum(gpm))
            return binid

        # Call the contributed code and return the bin index
#        _binid, xNode, yNode, xBar, yBar, sn, area, scale = \
#            voronoi_2d_binning(_x, _y, _signal, _noise, self.par['target_snr'],
#                               sn_func=sn_func, plot=False) #True, quiet=False)
        try:
            _binid, xNode, yNode, xBar, yBar, sn, area, scale = \
                voronoi_2d_binning(_x, _y, _signal, _noise, self.par['target_snr'],
                                   sn_func=sn_func, plot=False) #True, quiet=False)
        # TODO: Specify an exception type?
        except:
            warnings.warn('Binning algorithm has raised an exception.  Assume this is because '
                          'all the spaxels should be in the same bin.')
            _binid = numpy.zeros(_signal.size)
        if gpm is None:
            return _binid
        binid = numpy.full(x.size, -1, dtype=int)
        binid[gpm] = _binid
        return binid
# ----------------------------------------------------------------------


# SQUARE BINNING -------------------------------------------------------
# Written by Kate Rubin 9/27/18
# Based on KBW's test_new_binning_scheme.py and the RadialBinning class
class SquareBinningPar(KeywordParSet):
    """
    Class with parameters used by the square binning algorithm.  See
    :class:`mangadap.par.parset.ParSet` for attributes.

    The defined parameters are:

    .. include:: ../tables/squarebinningpar.rst

    """
    def __init__(self, binsz=2.0):
        pars = ['binsz']
        values = [binsz]
        dtypes = [float]
        descr = ['Desired bin size in arcsec']

        super().__init__(pars, values=values, dtypes=dtypes, descr=descr)

    def toheader(self, hdr):
        """
        Copy some of the parameters to a header.

        Args:
            hdr (`astropy.io.fits.Header`_):
                Header object to write to.

        Returns:
            `astropy.io.fits.Header`_: Edited header object

        Raises:
            TypeError:
                Raised if input is not an `astropy.io.fits.Header`_
                object.
        """
        if not isinstance(hdr, fits.Header):
            raise TypeError('Input is not a astropy.io.fits.Header object!')
        hdr['BINSIZE'] = (self['binsz'], 'Square bin size (arcsec)')
        return hdr

    def fromheader(self, hdr):
        """
        Copy the information from the header.

        Args:
            hdr (`astropy.io.fits.Header`_):
                Header object to write to.
        """
        if not isinstance(hdr, fits.Header):
            raise TypeError('Input is not a astropy.io.fits.Header object!')

        self['binsz'] = hdr['BINSIZE']


class SquareBinning(SpatialBinning):
    """
    Class to perform binning of full cube in square apertures
    Length of aperture side is given in arcsec with binsz
    """
    def __init__(self, par=None):
        SpatialBinning.__init__(self, 'square', par=par)
        self.binsz = None

    def bin_index(self, x, y, gpm=None, par=None):
        """
        Bin the data and return the indices of the bins.

        Args:
            x (`numpy.ndarray`_):
                Fiducial on-sky X position of each spectrum. X
                increases with RA.
            y (`numpy.ndarray`_):
                Fiducial on-sky Y position of each spectrum. Y
                increases with DEC.
            gpm (`numpy.ndarray`_, optional):
                Good value (pixel) mask used to select data that should be
                included in the bin.  If None, all values are considered valid.
                If provided, masked values (gpm=False) are given a bin index of
                -1.
            par (:class:`SquareBinningPar`, optional):
                Binning parameters.  Cannot be None.

        Returns:
            `numpy.ndarray`_: An integer bin index for each spectrum.

        Raises:
            ValueError:
                Raised if the sizes of ``x`` and ``y`` do not match or
                if ``par`` is None.
        """
        _x = numpy.asarray(x)
        if len(_x.shape) != 1:
            raise ValueError('On-sky coordinates must be one-dimensional.')
        if len(y) != _x.size:
            raise ValueError('Input coordinates are of different lengths.')
        _y = numpy.asarray(y)

        binid = numpy.full(_x.size, -1, dtype=int)
        if gpm is not None:
            if gpm.size != _x.size:
                raise ValueError(f'Provided GPM has incorrect shape; found {gpm.size}, ' \
                                 f'expected {x.size}.')
            _x = _x[gpm]
            _y = _y[gpm]

        # Perform any necessary initialization
        if par is not None:
            self.par = par
            if self.binsz is not None:
                del self.binsz
                self.binsz = None
        if self.par is None:
            raise ValueError('Required parameters for SquareBinning have not been defined.')
        if self.binsz is None:
            self.binsz = self.par['binsz']



        _binid = numpy.full(_x.size, -1, dtype=int)
        binsz = self.binsz
        minx = numpy.min(_x)-0.25
        maxx = numpy.max(_x)+0.25
        miny = numpy.min(_y)-0.25
        maxy = numpy.max(_y)+0.25

        ctbin = 0

        # Start from the center, work out to edge of each quadrant

        # Quadrant 1
        nslicex = numpy.ceil((maxx) / binsz)
        nslicey = numpy.ceil((maxy) / binsz)

        x_lo = numpy.linspace(0.0,(nslicex*binsz),nslicex+1.0)
        y_lo = numpy.linspace(0.0,(nslicey*binsz),nslicey+1.0)


        # Find which spaxels land in each aperture
        for xi in x_lo:
            for yj in y_lo:
                indx = (_x >= xi) & (_x < xi+binsz) & (_y >= yj) & (_y < yj+binsz)
                if any(indx):
                    _binid[indx] = ctbin
                    ctbin = ctbin + 1


        # Quadrant 2
        nslicex = numpy.ceil((maxx) / binsz)
        nslicey = numpy.floor((miny) / binsz)

        x_lo = numpy.linspace(0.0, (nslicex * binsz), nslicex + 1.0)
        y_lo = numpy.linspace(0.0, (nslicey * binsz), numpy.abs(nslicey) + 1.0)

        # Find which spaxels land in each aperture
        for xi in x_lo:
            for yj in y_lo:
                indx = (_x >= xi) & (_x < xi + binsz) & (_y <= yj) & (_y > yj - binsz)
                if any(indx):
                    _binid[indx] = ctbin
                    ctbin = ctbin + 1

        # Quadrant 3
        nslicex = numpy.floor((minx) / binsz)
        nslicey = numpy.ceil((maxy) / binsz)

        x_lo = numpy.linspace(0.0, (nslicex * binsz), numpy.abs(nslicex) + 1.0)
        y_lo = numpy.linspace(0.0, (nslicey * binsz), numpy.abs(nslicey) + 1.0)

        # Find which spaxels land in each aperture
        for xi in x_lo:
            for yj in y_lo:
                indx = (_x <= xi) & (_x > xi - binsz) & (_y >= yj) & (_y < yj + binsz)
                if any(indx):
                    _binid[indx] = ctbin
                    ctbin = ctbin + 1

        # Quadrant 4
        nslicex = numpy.floor((minx) / binsz)
        nslicey = numpy.floor((miny) / binsz)

        x_lo = numpy.linspace(0.0, (nslicex * binsz), numpy.abs(nslicex) + 1.0)
        y_lo = numpy.linspace(0.0, (nslicey * binsz), numpy.abs(nslicey) + 1.0)

        # Find which spaxels land in each aperture
        for xi in x_lo:
            for yj in y_lo:
                indx = (_x <= xi) & (_x > xi - binsz) & (_y <= yj) & (_y > yj - binsz)
                if any(indx):
                    _binid[indx] = ctbin
                    ctbin = ctbin + 1

        if gpm is None:
            return _binid
        binid[gpm] = _binid
        return binid


