# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""

Binning!

Revision history
----------------

    | **01 Apr 2016**: Implementation begun by K. Westfall (KBW)
    | **22 May 2018**: (KBW) Import vorbin package version of
        voronoi_2d_binning

----

.. include license and copyright
.. include:: ../copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../links.rst
"""
import warnings
import numpy
from scipy import sparse
from astropy.io import fits

from vorbin.voronoi_2d_binning import voronoi_2d_binning

from ..par.parset import KeywordParSet
from ..util.geometry import SemiMajorAxisCoo
from ..util.covariance import Covariance

from matplotlib import pyplot

# Add strict versioning
# from distutils.version import StrictVersion


# BASE CLASS -----------------------------------------------------------

class SpatialBinning():
    """
    Base class for spatially binning data.
    """
    def __init__(self, bintype, par=None):
        self.bintype = bintype
        self.par = par


# ----------------------------------------------------------------------


# GLOBAL BINNING -------------------------------------------------------

class GlobalBinning(SpatialBinning):
    """
    Class that performs the global binning.
    """
    def __init__(self):
        SpatialBinning.__init__(self, 'global')

    @staticmethod
    def bin_index(x, y, par=None):
        return numpy.full(x.shape, 0, dtype=int)


# ----------------------------------------------------------------------


# RADIAL BINNING -------------------------------------------------------

class RadialBinningPar(KeywordParSet):
    """
    Class with parameters used by the radial binning algorithm.  See
    :class:`mangadap.par.parset.ParSet` for attributes.

    The defined parameters are:

    .. include:: ../tables/radialbinningpar.rst
    """
    def __init__(self, center=None, pa=None, ell=None, radius_scale=None, radii=None,
                 log_step=None):
        in_fl = [ int, float ]
        ar_like = [ numpy.ndarray, list ]
        
        pars =   [ 'center',  'pa', 'ell', 'radius_scale', 'radii', 'log_step' ]
        values = [   center,    pa,   ell,   radius_scale,   radii,   log_step ]
        dtypes = [  ar_like, in_fl, in_fl,          in_fl, ar_like,       bool ]

        descr = ['A two-element array defining the center to use in the definition of the ' \
                    'elliptical bins.  This is defined as a sky-right offset in arcseconds from ' \
                    'the nominal center of the object.',
                 'Sets the position angle, defined from N through E of the major axis of the ' \
                    'isophotal ellipse used to define the elliptical bins.',
                 'Sets the ellipticity (1-b/a) of the isophotal ellipse use to define the ' \
                    'elliptical bins.',
                 'Defines a scale factor to use when defining the radial bins.  For example, ' \
                    'you might want to scale to the a certain number of effective radii or ' \
                    'physical scale in kpc.  For no scale, use 1.0.',
                 'A three-element array defining the starting and ending radius for the bin ' \
                    'edges and the number of bins to create.  If the starting radius is -1, ' \
                    'the inner-most radius is set to 0 when not using log bins or 0.1 arcsec ' \
                    'when using logarithmically spaced bins.  If the ending radius is -1, the ' \
                    'outer-most radius is set by the spaxel at the largest radius.',
                 'A flag that the radial bins should be a geometric series.']

        super(RadialBinningPar, self).__init__(pars, values=values, dtypes=dtypes, descr=descr)
    
    def toheader(self, hdr):
        """
        Copy some of the parameters to a header.

        Args:
            hdr (`astropy.io.fits.Header`_): Header object to write to.

        Returns:
            `astropy.io.fits.Header`_: Edited header object

        Raises:
            TypeError: Raised if input is not an
                `astropy.io.fits.Header`_ object.
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

        hdr (`astropy.io.fits.Header`_): Header object to write to.
        """
        if not isinstance(hdr, fits.Header):
            raise TypeError('Input is not a astropy.io.fits.Header object!')

        self['center'] = [ hdr['BINCX'], hdr['BINCY'] ]
        self['pa'] = hdr['BINPA']
        self['ell'] = hdr['BINELL']
        self['radius_scale'] = hdr['BINSCL']
        self['radii'] = eval(hdr['BINRAD'])
        self['log_step'] = bool(hdr['BINLGR'])


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

        For ending radii, just swap rs and re above
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
        

    def bin_index(self, x, y, par=None):
        """
        Bin the data and return the indices of the bins.
        """
        # Check the input
        if x.size != y.size:
            raise ValueError('Dimensionality of x and y coordinates do not match!')

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

#        print(binid)
#        print(numpy.sum(binid < 0))
#
#        print(self.rs, self.dr, self.to_center, self.par['log_step'])
#        pyplot.scatter(r, theta, marker='.', color='k', s=40, lw=0, c=binid, cmap='viridis')
#        for i in range(self.par['radii'][2]+1):
#            rr = self.rs + i*self.dr
#            pyplot.plot( [rr, rr], [0,360], color='r')
#        pyplot.show()

        return binid


    def bin_area(self):
        """
        Return the nominal area for each elliptical bin.
        """
        if self.par is None:
            raise ValueError('Required parameters not defined.')

        rs = self.par['radii'][0]
        re = self.par['radii'][1]
        nr = int(self.par['radii'][2])
        bin_radii = numpy.append(
                        numpy.logspace(numpy.log10(rs), numpy.log10(re), num=nr, endpoint=False) \
                        if self.par['log_step'] else numpy.linspace(rs, re, num=nr, endpoint=False),
                        re)*self.par['radius_scale']
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
    def __init__(self, target_snr=None, signal=None, noise=None, covar=None):
        in_fl = [ int, float ]
        covar_type = [ float, numpy.ndarray, Covariance, sparse.spmatrix ]
        ar_like = [ numpy.ndarray, list ]
        
        pars =   [ 'target_snr', 'signal', 'noise',    'covar' ]
        values = [   target_snr,   signal,   noise,      covar ]
        dtypes = [        in_fl,  ar_like, ar_like, covar_type ]

        descr = ['The target S/N for each bin.',
                 'The array of signal measurements for each on-sky position to bin.',
                 'The array of noise measurements for each on-sky position to bin.  If not ' \
                    'provided, ``covar`` must be provided and be a full covariance matrix.',
                 r'Covariance matrix or calibration normalization.  For the latter, the value ' \
                    r'is used to renormalize using :math:`n_{\rm calib} = n_{\rm nominal} ' \
                    r'(1 + \alpha\ \log\ N_{\rm bin})`, where :math:`N_{\rm bin}` is the number ' \
                    r'of binned spaxels and :math:`\alpha` is the value provided.']

        super(VoronoiBinningPar, self).__init__(pars, values=values, dtypes=dtypes, descr=descr)

    def toheader(self, hdr):
        """
        Copy some of the parameters to a header.

        Args:
            hdr (`astropy.io.fits.Header`_): Header object to write to.

        Returns:
            `astropy.io.fits.Header`_: Edited header object

        Raises:
            TypeError: Raised if input is not an
                `astropy.io.fits.Header`_ object.
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

        hdr (`astropy.io.fits.Header`_): Header object to write to.
        """
        if not isinstance(hdr, fits.Header):
            raise TypeError('Input is not a astropy.io.fits.Header object!')

        self['target_snr'] = hdr['BINSNR']
        if hdr['BINCOV'] == 'calib':
            self['covar'] = hdr['NCALIB']


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
        """
        return  numpy.sum(signal[index]) / numpy.sqrt(numpy.sum(numpy.square(noise[index])))


    def sn_calculation_covariance_matrix(self, index, signal, noise):
        """
        Calculate the S/N using a full covariance matrix.
        """
#        _index = numpy.arange(signal.size)[index]
#        i, j = map( lambda x: x.ravel(), numpy.meshgrid(_index, _index) )
#        return numpy.sum(signal[index])/numpy.sqrt(numpy.sum(self.covar[i,j]))
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
        """
        return numpy.sum(signal[index]) \
                / (numpy.sqrt(numpy.sum(numpy.square(noise[index]))) \
                    * (1.0 + self.covar*numpy.log10(len(signal[index]))))


    def bin_index(self, x, y, par=None):
        # Check the position input
        if x.size != y.size:
            raise ValueError('Dimensionality of x and y coordinates do not match!')

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
#                self.par['covar'].show()
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

        # All spaxels have S/N greater than threshold, so return each
        # spaxel in its own "bin"
        if numpy.min(self.par['signal']/_noise) > self.par['target_snr']:
            warnings.warn('All pixels have enough S/N. Binning is not needed')
            return numpy.arange(self.par['signal'].size)

        # Cannot reach the S/N using all spaxels, so return all spaxels
        # in a single bin
        sn_total = sn_func(numpy.arange(self.par['signal'].size), self.par['signal'], _noise)
        if sn_total < self.par['target_snr']:
            warnings.warn('Cannot reach target S/N using all data.')
            return numpy.zeros(self.par['signal'].size)

        # Call the contributed code and return the bin index
#        pyplot.scatter(x, _noise, s=30, marker='.', color='k')
#        pyplot.show()
        try:
            binid, xNode, yNode, xBar, yBar, sn, area, scale = \
                voronoi_2d_binning(x, y, self.par['signal'], _noise, self.par['target_snr'],
                                   sn_func=sn_func, plot=False) #True, quiet=False)
#            pyplot.show()
#            pyplot.scatter(numpy.sqrt(numpy.square(xBar)+numpy.square(yBar)), sn, marker='.',
#                           color='k', s=40, lw=0)
#            pyplot.show()
        except:
            warnings.warn('Binning algorithm has raised an exception.  Assume this is because '
                          'all the spaxels should be in the same bin.')
            binid = numpy.zeros(self.par['signal'].size)

        return binid


# ----------------------------------------------------------------------


# SQUARE BINNING -------------------------------------------------------

#Written by Kate Rubin 9/27/18
#Based on KBW's test_new_binning_scheme.py and the RadialBinning class

class SquareBinningPar(KeywordParSet):
    """
    Class with parameters used by the square binning algorithm.  See
    :class:`mangadap.par.parset.ParSet` for attributes.

    The defined parameters are:

    .. include:: ../tables/squarebinningpar.rst

    """

    def __init__(self, binsz=None):
        in_fl = [int, float]

        pars = ['binsz']
        values = [binsz]
        dtypes = [float]
        descr = ['Desired bin size in arcsec']

        super(SquareBinningPar, self).__init__(pars, values=values, dtypes=dtypes, descr=descr)

    def toheader(self, hdr):
        """
        Copy some of the parameters to a header.

        Args:
            hdr (`astropy.io.fits.Header`_): Header object to write to.

        Returns:
            `astropy.io.fits.Header`_: Edited header object

        Raises:
            TypeError: Raised if input is not an
                `astropy.io.fits.Header`_ object.
        """
        if not isinstance(hdr, fits.Header):
            raise TypeError('Input is not a astropy.io.fits.Header object!')

        hdr['BINSIZE'] = (self['binsz'], 'Square bin size (arcsec)')
        return hdr

    def fromheader(self, hdr):
        """
        Copy the information from the header.

        hdr (`astropy.io.fits.Header`_): Header object to write to.
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

    def bin_index(self, x, y, par=None):
        """
            Bin the data and return the indices of the bins.
        """
        _x = numpy.asarray(x)
        if len(_x.shape) != 1:
            raise ValueError('On-sky coordinates must be one-dimensional.')
        nspaxels = _x.size
        if len(y) != nspaxels:
            raise ValueError('Input coordinates are of different lengths.')
        _y = numpy.asarray(y)


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



        binid = numpy.full(nspaxels, -1, dtype=int)
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
                    binid[indx] = ctbin
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
                    binid[indx] = ctbin
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
                    binid[indx] = ctbin
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
                    binid[indx] = ctbin
                    ctbin = ctbin + 1


        return binid


