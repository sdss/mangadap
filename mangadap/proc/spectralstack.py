# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""

Stack some spectra!

Revision history
----------------

    | **01 Apr 2016**: Implementation begun by K. Westfall (KBW)
    | **10 Nov 2016**: (KBW) Include dispersion cube in DRP file to
        construct spectral resolution of binned spectra.
    | **30 Nov 2016**: (KBW) Use, e.g., [None,:] slicing instead of,
        e.g., numpy.array([v]*n) when necessary for array arithmetic
    | **01 Dec 2016**: (KBW) Allow stacking the spectral resolution to
        be turned off via :class:`SpectralStackPar`.
    | **06 Dec 2016**: (KBW) In
        :func:`mangadap.proc.spectralstack.SpectralStack._set_rebin_transfer_matrix`,
        the number of bins is set by the number of unique indices;
        before based on the maximum unique index, meaning that "missing"
        bins were included.  They're now excluded, forcing classes like
        :class:`mangadap.proc.spatiallybinnedspectra.SpatiallyBinnedSpectra`
        to keep track of missing bins.
    | **30 Aug 2017**: (KBW) Switch from
        :func:`mangadap.util.instrument.resample_vector` to
        :func:`mangadap.util.instrument.resample1d` in untested function
        :func:`SpectralStack.register`.
    | **30 Aug 2018**: (KBW) Switch from resample1d to
        :func:`mangadap.util.sampling.Resample`

----

.. include license and copyright
.. include:: ../copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../links.rst
"""

import numpy
from scipy import sparse
from astropy.io import fits
import astropy.constants

from ..par.parset import KeywordParSet
from ..util.covariance import Covariance
from ..util.filter import interpolate_masked_vector
from ..util.sampling import Resample, spectral_coordinate_step

from matplotlib import pyplot, rc

# Add strict versioning
# from distutils.version import StrictVersion

class SpectralStackPar(KeywordParSet):
    r"""
    Class with parameters used to set how to stack a set of spectra.
    See :class:`mangadap.par.parset.ParSet` for attributes.

    .. todo::
        Allow for sigma rejection.

    The defined parameters are:

    .. include:: ../tables/spectralstackpar.rst
    """
    def __init__(self, operation=None, vel_register=None, vel_offsets=None, covar_mode=None,
                 covar_par=None):
        in_fl = [ int, float ]
        ar_like = [ numpy.ndarray, list ]
        op_options = SpectralStack.operation_options()
        covar_options = SpectralStack.covariance_mode_options()
        
        pars =     ['operation', 'vel_register', 'vel_offsets',  'covar_mode',   'covar_par']
        values =   [  operation,   vel_register,   vel_offsets,    covar_mode,     covar_par]
        defaults = [     'mean',          False,          None,        'none',          None]
        options =  [ op_options,           None,          None, covar_options,          None]
        dtypes =   [        str,           bool,       ar_like,           str, in_fl+ar_like]
        descr = ['Operation to perform for the stacked spectrum.  See ' \
                    ':func:`SpectralStack.operation_options` for the available operation options.',
                 'Flag to velocity register the spectra before adding them based on a provided ' \
                    'prior measurement of the velocities.',
                 'List of velocity offsets to apply to the spectra to stack.',
                 'Describes how to incorporate covariance into the spectral stacking.  ' \
                    'See :func:`SpectralStack.covariance_mode_options` for the available options.',
                 'The parameter(s) needed to perform a given method of handling the ' \
                    'covariance.  See :func:`SpectralStack.covariance_mode_options` for the ' \
                    'available options.']
        super(SpectralStackPar, self).__init__(pars, values=values, defaults=defaults,
                                               options=options, dtypes=dtypes, descr=descr)

    def toheader(self, hdr):
        """
        Copy the information to a header.

        Args:
            hdr (`astropy.io.fits.Header`_):
                Header object to write to.

        Returns:
            `astropy.io.fits.Header`_ : Edited header object.
        """
        hdr['STCKOP'] = (self['operation'], 'Stacking operation')
        hdr['STCKVREG'] = (str(self['vel_register']), 'Spectra shifted in velocity before stacked')
        hdr['STCKCRMD'] = (str(self['covar_mode']), 'Stacking treatment of covariance')
        hdr['STCKCRPR'] = (str(self['covar_par']), 'Covariance parameter(s)')
        return hdr

    def fromheader(self, hdr):
        """
        Copy the information from the header

        Args:
            hdr (`astropy.io.fits.Header`_):
                Header object to read from.
        """
        self['operation'] = hdr['STCKOP']
        self['vel_register'] = eval(hdr['STCKVREG'])
        self['covar_mode'] = hdr['STCKCRMD']
        self['covar_par'] = eval(hdr['STCKCRPR'])


class SpectralStack:
    r"""
    Class whose primary function is to stack a set of spectra while
    treating covariance between spectra.

    See :func:`covariance_mode_options` for available methods of
    accounting for covariance.

    Class also approximates the resulting spectral resolution of the
    stacked data.

    .. todo::
        List attributes

    """
    def __init__(self):
        # Keep the matrix used to bin the spectra. This will be a
        # scipy.sparse.csr_matrix.
        self.rebin_T = None

        # Internal arrays for callback
        self.wave = None
        self.flux = None
        self.fluxsqr = None
        self.fluxmean = None
        self.fluxsdev = None
        self.npix = None
        self.ivar = None
        self.sres = None
        self.covar = None

    @staticmethod
    def _check_covariance_type(covariance_mode, covar, ivar):
        """
        Check that the covariance variable has the correct type for
        the given mode.

        Args:
            covariance_mode (:obj:`str`):
                Covariance handling mode; see
                :func:`covariance_mode_options`.
            covar (None, :obj:`float`, :class:`mangadap.util.covariance.Covariance`):
                The object to check the type against the covariance
                handling mode.
            ivar (None, object):
                Inverse variance. Only check performed is whether or
                not this is None.
            
        Returns:
            :obj:`bool`: Flag that type is correct.
        """
        if covariance_mode == 'none':
            return True
        if covariance_mode in ['calibrate', 'channels', 'wavelengths'] and ivar is None:
            return False
        if covariance_mode == 'calibrate' and not isinstance(covar, float):
            return False
        if covariance_mode == 'calibrate':
            return True
        if not isinstance(covar, Covariance):
            return False
        return covar.dim == 3

    @staticmethod
    def _check_covariance_shape(covariance_mode, covar, nwave, nspec):
        """
        Check that the input covariance object has the correct shape
        for the given mode.

        Args:
            covariance_mode (:obj:`str`):
                Covariance handling mode; see
                :func:`covariance_mode_options`.
            covar (None, :obj:`float`, :class:`mangadap.util.covariance.Covariance`):
                The object to check the type against the covariance
                handling mode.
            nwave (:obj:`int`):
                Number of wavelength channels.
            nspec (:obj:`int`):
                Number of spectra.

        Returns:
            :obj:`bool`: Flag that the covariance data has the
            correct shape.
        """
        if covariance_mode in [ 'none', 'calibrate']:
            return True
        if covariance_mode in [ 'full', 'approx_correlation' ] \
                and (covar.dim != 3 or covar.shape[-1] != nwave):
            return False
        if covar.shape[0] != nspec:
            return False
        return True

    @staticmethod
    def _get_input_mask(flux, ivar=None, mask=None, dtype=bool):
        """
        Construct the baseline mask using the input arrays.

        Args:
            flux (`numpy.ndarray`_, `numpy.ma.MaskedArray`_):
                Input flux array.
            ivar (`numpy.ndarray`_):
                Inverse variance array. Mask excludes any inverse
                variance values that are not greater than 0.
            mask (`numpy.ndarray`_):
                Baseline boolean mask array.
            dtype (data type):
                Output data type for array.
        
        Returns:
            `numpy.ndarray`_: Boolean mask array converted to
            ``dtype``.
        """
        inp_mask = numpy.zeros(flux.shape, dtype=bool) if mask is None else mask
        if isinstance(flux, numpy.ma.MaskedArray):
            inp_mask |= numpy.ma.getmaskarray(flux)
        if ivar is not None:
            inp_mask |= numpy.invert(ivar>0)
        return inp_mask.astype(dtype)

    @staticmethod
    def _check_input_sres(sres, nspec):
        r"""
        Check the shape and type of the input spectral resolution.

        If the input is a masked array, interpolate over the mask.

        Args:
            sres (`numpy.ndarray`_, `numpy.ma.MaskedArray`_):
                Input spectral resolution data. If None, None is
                returned.
            nspec (:obj:`int`):
                Expected number of spectra.

        Returns:
            `numpy.ndarray`_: The interpolated spectral resolution
            vectors. If the input ``sres`` is a single vector, the
            spectral resolution is repeated for each spectrum. Shape
            is :math:`(N_{\rm spec}, N_{\rm wave})`. If the input
            ``sres`` is None, instead None is returned.
        """
        if sres is None:
            return None
        if sres.ndim > 2:
            raise ValueError('Spectral resolution must be 2D or less.')
        if sres.ndim == 2 and nspec > 0 and sres.shape[0] != nspec:
            raise ValueError('Spectral resolution array is not the correct shape.')
        _sres = sres
        if isinstance(sres, numpy.ma.MaskedArray):
            _sres = numpy.apply_along_axis(interpolate_masked_vector, 1,
                                           sres.reshape(1,-1) if sres.ndim == 1 else sres)
            if sres.ndim == 1:
                _sres = _sres[0,:]
        if _sres.ndim == 2 or nspec == 0:
            return _sres
        return numpy.array([_sres]*nspec)

    def _set_rebin_transfer_matrix(self, binid, binwgt=None):
        r"""
        Construct the transfer matrix that rebins the spectra.

        The shape of the transfer matrix is :math:`(N_{\rm bin}
        \times N_{\rm spec})` with the expectation that the spectrum
        flux array has shape :math:`(N_{\rm spec} \times N_{\rm
        wave})`.

        The binned spectra are calculated by matrix multiplication,
        :math:`\mathbf{B} = \mathbf{T} \times \mathbf{F}` such that
        the covariance matrix can be calculated as :math:`\mathbf{C}
        = \mathbf{T} \times \mathbf{\Sigma} \times \mathbf{T}^{\rm
        T}`, where :math:`\mathbf{\Sigma}` is the covariance matrix
        in the flux array, :math:`\mathbf{F}`.

        If weighting, the sum of the weights is normalized to the
        number of points included in the bin.

        Args:
            binid (`numpy.ndarray`_):
                Index, one per spectrum in the flux array, for the
                binned spectrum. Indices of less than one are
                ignored.
            binwgt (`numpy.ndarray`_, optional):
                List of weights for the spectra. If not provided, the
                weights are uniform.
        """
        nspec = binid.size
        valid = binid > -1
        unique_bins = numpy.unique(binid[valid])
        nbin = len(unique_bins)
        self.rebin_T = numpy.zeros((nbin,nspec), dtype=numpy.float)
        for j in range(nbin):
            indx = binid == unique_bins[j]
            self.rebin_T[j,indx] = 1.0 if binwgt is None else \
                                    binwgt[indx]*numpy.sum(indx)/numpy.sum(binwgt[indx])
        self.rebin_T = sparse.csr_matrix(self.rebin_T)

    def _stack_without_covariance(self, flux, ivar=None, sres=None):
        """
        Stack the spectra, ignoring any covariance that may or may not
        exist.

        Sets :attr:`flux`, :attr:`fluxsqr`, :attr:`npix`, :attr:`ivar`,
        :attr:`sres`.

        The stored data is always based on the **sum** of the spectra
        in the stack.
        
        Args:
            flux (`numpy.ma.MaskedArray`_):
                Flux array.
            ivar (:obj:`numpy.ma.MaskedArray`, optional):
                The inverse variance array.
            sres (`numpy.ma.MaskedArray`_, optional):
                1D or 2D spectral resolution as a function of
                wavelength for all or each input spectrum. Default is
                to ignore any spectral resolution data.
        """
        # Calculate the sum of the flux, flux^2, and determine the
        # number of pixels in the sum
        # NOTE: In these numpy.ma expressions, masked values are
        # ignored, and output values are masked only if all values in
        # the arithmetic expression are masked
        rt = self.rebin_T.toarray()
        self.flux = numpy.ma.dot(rt, flux)
        self.fluxsqr = numpy.ma.dot(rt, numpy.square(flux))
        self.npix = numpy.ma.dot(rt, numpy.invert(numpy.ma.getmaskarray(flux))).astype(int)

        # Calculate the spectral resolution
        self.sres = None
        if sres is not None:
            Tc = numpy.sum(rt, axis=1)
            Tc[numpy.invert(Tc>0)] = 1.0
            self.sres = numpy.ma.power(numpy.ma.dot(rt, numpy.ma.power(sres, -2))/Tc[:,None], -0.5)

        # No inverse variance so we're done
        if ivar is None:
            self.ivar = None
            return

        # Calculate the propagated inverse variance (assumes no
        # covariance)
        self.ivar = numpy.ma.power(numpy.ma.dot(numpy.power(rt, 2.0), numpy.ma.power(ivar, -1.)),
                                   -1.)

    def _stack_with_covariance(self, flux, covariance_mode, covar, ivar=None, sres=None):
        """
        Stack the spectra and incorporate covariance.
        
        Args:
            flux (`numpy.ma.MaskedArray`_):
                Flux array.
            covariance_mode (:obj:`str`):
                Covariance handling mode; see
                :func:`covariance_mode_options`.
            covar (None, :obj:`float`, :class:`mangadap.util.covariance.Covariance`):
                The relevant covariance object that must match the
                needs of the covariance handling mode.
            ivar (:obj:`numpy.ma.MaskedArray`, optional):
                The inverse variance array. Must not be None if
                ``covariance_mode`` is 'channels' or 'wavelengths'.
            sres (`numpy.ma.MaskedArray`_, optional):
                1D or 2D spectral resolution as a function of
                wavelength for all or each input spectrum. Default is
                to ignore any spectral resolution data.
        """
        # First stack without covariance. This sets self.flux,
        # self.fluxsqr, self.npix, self.ivar, self.sres.
        self._stack_without_covariance(flux, ivar=ivar, sres=sres)

        # If calibrating the noise based on the equation,
        #   Noise_corrected = Noise_nocovar * (1 + f * log10(N))
        # apply it and return
        if covariance_mode == 'calibrate':
            if self.ivar is not None:
                self.ivar /= numpy.square(1.0 + covar*numpy.ma.log10(self.npix))
            return

        # Check that the code knows what to do otherwise
        # TODO: Isn't this check done elsewhere?
        if covariance_mode not in ['channels', 'wavelengths', 'approx_correlation', 'full']:
            raise ValueError('Unknown covariance mode: {0}'.format(covariance_mode))

        # Recalibrate the error based on a selection of covariance
        # channels
        recalibrate_ivar = covariance_mode in ['channels', 'wavelengths']
        if recalibrate_ivar and self.ivar is None:
            raise ValueError('Must provide ivar to recalibrate based on covar.')

        # Setup for output
        nbin = self.flux.shape[0]
        nchan = covar.shape[-1]
        self.covar = numpy.empty(nchan, dtype=sparse.csr.csr_matrix)
        variance_ratio = numpy.ma.zeros( (nbin,nchan), dtype=numpy.float) \
                            if recalibrate_ivar else None

        # Calculate the covariance in the stack
        for i in range(nchan):
            j = covar.input_indx[i]
            cov_array = covar.with_lower_triangle(channel=j)
            self.covar[i] = sparse.triu(self.rebin_T.dot(cov_array.dot(self.rebin_T.T))).tocsr()

            # Get the variance ratio
            if recalibrate_ivar:
                variance_ratio[:,i] = self.covar[i].diagonal() * self.ivar[:,j]
                variance_ratio[numpy.ma.getmaskarray(self.ivar)[:,j],i] = numpy.ma.masked

        self.covar = Covariance(self.covar, input_indx=covar.input_indx)

        # Set ivar by recalibrating the existing data
        if recalibrate_ivar:
            ratio = numpy.ma.median( variance_ratio, axis=1 )
            # TODO: Add a debug mode for this?
#            self._recalibrate_ivar_figure(ratio, ofile='ivar_calibration.pdf')
            self.ivar = numpy.ma.power(ratio, -1.0)[:,None] * self.ivar
            return

        # Get the inverse variance from the full covariance matrix
        self.ivar = numpy.ma.power(self.covar.variance(), -1.).filled(0.0)


    def _recalibrate_ivar_figure(self, ratio, ofile=None):
        """
        Make a plot showing the computed inverse variance
        recalibration based on the covariance matrices.

        Args:
            ratio (`numpy.ndarray`_):
                Ratio between the error with and without covariance.
            ofile (:obj:`str`, optional):
                File for the plot.  If None, output to the screen.
        """
        font = { 'size' : 20 }
        rc('font', **font)

        w,h = pyplot.figaspect(1)
        fig = pyplot.figure(figsize=(1.5*w,1.5*h))

        ax = fig.add_axes([0.2,0.3,0.7,0.4])
        ax.minorticks_on()
        ax.tick_params(which='major', length=10, direction='in', top=True, right=True)
        ax.tick_params(which='minor', length=5, direction='in', top=True, right=True)
        ax.grid(True, which='major', color='0.8', zorder=0, linestyle='-', lw=0.5)

        nbin = numpy.sum(self.rebin_T.toarray(), axis=1)
        mod_nbin = numpy.arange(0, numpy.amax(nbin), 0.1)+1
        ax.scatter(nbin, numpy.sqrt(ratio), marker='.', s=50, lw=0, color='k', zorder=3)
        ax.plot(mod_nbin, 1+1.62*numpy.log10(mod_nbin), color='C3', zorder=2)

        ax.text(0.5, -0.18, r'$N_{\rm bin}$', horizontalalignment='center',
                verticalalignment='center', transform=ax.transAxes)
        ax.text(-0.14, 0.5, r'$f_{\rm covar}$',
                horizontalalignment='center', verticalalignment='center',
                transform=ax.transAxes, rotation='vertical')
        ax.text(0.96, 0.1, '8249-12705', horizontalalignment='right',
                verticalalignment='center', transform=ax.transAxes)

        if ofile is None:
            pyplot.show()
        else:
            fig.canvas.print_figure(ofile, bbox_inches='tight')
        fig.clear()
        pyplot.close(fig)

    def _covar_in_mean(self):
        """
        Compute the covariance in the mean spectrum by propagating the
        division by the number of pixels through the covariance matrix.

        Returns:
            :class:`mangadap.util.covariance.Covariance`: Covariance
            in the mean spectrum. Returns None if :attr:`covar` is
            None.
        """
        if self.covar is None:
            return None
        nchan = self.covar.shape[-1]
        nbin = self.flux.shape[0]
        inpix = numpy.ma.power(self.npix, -1.)
        covar = numpy.empty(nchan, dtype=sparse.csr.csr_matrix)
        for i in range(nchan):
            j = self.covar.input_indx[i]
            _inpix = inpix[:,j,None]*inpix[None,:,j]
            covar[i] = sparse.triu(self.covar.toarray(channel=j) * _inpix).tocsr()
        return Covariance(covar, input_indx=self.covar.input_indx)
            
    def _get_stack_mean(self):
        """
        Convert the summed stack to the mean stacked spectra using
        the internal data.

        Returns:
            :obj:`tuple`: See the return statement for :func:`stack`.
        """
        return self.wave, self.fluxmean, self.fluxsdev, self.npix, \
                        self.ivar * numpy.square(self.npix), self.sres, self._covar_in_mean()

    @staticmethod
    def operation_options():
        """
        Return the allowed stacking operations.
        
        Current operations are:
            - ``mean``: Construct the mean of the spectra
            - ``sum``: Construct the spectrum sum.

        Returns:
            :obj:`list`: List of available operations.
        """
        return ['mean', 'sum']

    @staticmethod
    def covariance_mode_options(par_needed=False):
        r"""
        Return the list of allowed covariance options.

        The two parameters ``covar_mode`` and ``covar_par`` (see
        :class:`SpectralStackPar`) set how covariance is accounted
        for in the stacking procedure. The valid options are:

            - ``none``: The noise in the stacked spectrum is a
              nominal propagation of the error assuming no
              covariance. No parameters needed.

            - ``calibrate``: Where :math:`N_{\rm bin}` is the number
              of binned spaxels and :math:`\alpha` as a provided
              parameter, the spectral noise is calibrated following:

              .. math::

                    n_{\rm calib} = n_{\rm nominal} (1 + \alpha \log\
                    N_{\rm bin})

            - ``channels``: The noise vector of each stacked spectrum
              is adjusted based on the mean ratio of the nominal and
              formally correct calculations of the noise measurements
              over a number of spectral channels. The channels are
              drawn from across the full spectral range. The number
              of channels to use is a defined parameter. The
              covariance matrix must be provided to :func:`stack`.

            - ``wavelengths``: Functionally equivalent to
              ``channels``; however, the channels to use are set by a
              list of provided wavelengths. The covariance matrix
              must be provided to :func:`stack`.

            - ``approx_correlation``: Approximate the covariance
              matrix using a Gaussian description of the correlation
              between pixels. See
              :func:`mangadap.datacube.datacube.DataCube.approximate_correlation_matrix`.
              The value of :math:`\sigma` provides for the Gaussian
              desciption of :math:`\rho_{ij}` in the correlation
              matrix. The covariance matrix must be provided to
              :func:`stack`.

            - ``full``: The full covariance cube is calculated and
              the noise vectors are constructed using the formally
              correct calculation. No parameters needed. The
              covariance matrix must be provided to :func:`stack`.

        Returns:
            :obj:`list`: List of the allowed modes.
        """
        modes = ['calibrate', 'approx_correlation', 'channels', 'wavelengths']
        if par_needed:
            return modes
        return modes + ['none', 'full']

    @staticmethod
    def parse_covariance_parameters(mode, par):
        """
        Parse the parameters needed for the treatment of the
        covariance when stacking spectra.
    
        Args:
            mode (:obj:`str`):
                Mode to use. Must be an allowed mode; see
                :func:`covariance_mode_options`.
            par (:obj:`str`):
                String representation of the parameters for the
                specified mode.

        Returns:
            :obj:`float`, :obj:`list`: Parameters parsed from the
            input string for the designated covariance mode.

        Raises:
            TypeError:
                Raised if the input parameter could not be converted
                to a float as needed by the specified mode.
            ValueError:
                Raised if the mode is not recognized.
        """
        mode_options = SpectralStack.covariance_mode_options()
        if mode not in mode_options:
            raise ValueError('Mode not among valid options: {0}.\nOptions are: {1}'.format(mode,
                                                                                    mode_options))
        if mode in ['none', 'full']:
            return None
        if mode in ['calibrate', 'approx_correlation']:
            try:
                return float(par)
            except:
                raise TypeError('Could not convert to float: {0}'.format(par))
        if mode == 'channels':
            return int(par) #[ int(e.strip()) for e in par.split(',') ]
        if mode == 'wavelengths':
            return [ float(e.strip()) for e in par.split(',') ]

    @staticmethod
    def min_max_wave(wave, voff):
        r"""
        Determine the minimum and maximum of all shifted wavelength
        ranges.

        Args:
            wave (array-like):
                Original wavelengths.  Should be 1D.
            voff (:obj:`float`, array-like):
                The velocity shift in km/s to apply to the wavelength
                vector (i.e., the velocity shift should be
                :math:`v_{\rm off} = -cz`). Each element is applied
                to the wavelength vector to determine the maximum
                wavelength range allowed for all spectra. Should be
                1D.

        Returns:
            :obj:`tuple`: Two floats with the minimum and maximum
            redshifted wavelengths.

        Raises:
            ValueError:
                Raised if either ``wave`` or ``voff`` have the wrong
                dimensionality.
        """
        _wave = numpy.atleast_1d(wave)
        if _wave.ndim != 1:
            raise ValueError('Wavelength vector should be 1D!')
        _voff = numpy.atleast_1d(voff)
        if _voff.ndim != 1:
            raise ValueError('Velocity should be a float or 1D vector.')
        _wave = numpy.array([numpy.amin(_wave), numpy.amax(_wave)])
        doppler_shifted = _wave[None,:]*(1.+_voff[:,None]/astropy.constants.c.to('km/s').value)
        return numpy.amin(doppler_shifted), numpy.amax(doppler_shifted)

    @staticmethod
    def register(wave, voff, flux, ivar=None, mask=None, sres=None, log=False, base=10.0,
                 keep_range=False, flim=0.5):
        r"""
        Register a set of spectra to the same wavelength range given a
        set of velocity offsets.

        .. warning::
            **THIS FUNCTION IS UNTESTED!**

        .. todo::
            - Allow for correction for deredshifting flux.

        Args:
            wave (`numpy.ndarray`_):
                Single wavelength vector for all input spectra. Must
                be 1D with shape :math:`(N_{\rm wave},)`.
            voff (:obj:`float`, array-like):
                The velocity shift in km/s to apply to the wavelength
                vector (i.e., the velocity shift should be
                :math:`v_{\rm off} = -cz`). Should be 1D at most.
            flux (`numpy.ndarray`_, `numpy.ma.MaskedArray`_):
                Flux array to register. Must be 2D with shape
                :math:`(N_{\rm spec},N_{\rm wave})`.
            ivar (`numpy.ndarray`_, `numpy.ma.MaskedArray`_, optional):
                Inverse variance in the spectrum fluxes. Shape must
                match ``flux``.
            mask (`numpy.ndarray`_, optional):
                Boolean array with values in ``flux`` to mask.
                Default assumes nothing is masked. If ``flux`` and/or
                ``ivar`` are masked array, this is included in a
                union with those masks.  Shape must match ``flux``.
            sres (`numpy.ndarray`_, `numpy.ma.MaskedArray`_, optional):
                1D or 2D spectral resolution as a function of
                wavelength for all or each input spectrum. Default is
                to ignore any spectral resolution data. If a masked
                array, the masked pixels are interpolated.
            log (:obj:`bool`, optional):
                Flag that the wavelength vector is sampled
                geometrically in wavelength.
            base (:obj:`float`, optional):
                If the wavelength vector is geometrically sampled,
                this is the base of the logarithmic step.
            keep_range (:obj:`bool`, optional):
                When registering the wavelengths of the shifted
                spectra, keep the identical spectral range as input.
            flim (:obj:`float`, optional):
                Mask any pixels in the output flux arrays that are
                covered by less than this fraction by the input
                masked flux array.

        Returns:
            :obj:`tuple`:Returns four arrays: (1) the new wavelength
            vector common to all spectra; (2) the new masked flux
            array; (3) the new masked inverse variance array, which
            will be None if no inverse variances are provided to the
            function; and (4) the new spectral resolution vectors,
            which will be None if no spectral resolution vectors are
            provided.

        Raises:
            ValueError:
                Raised if the wavelength or velocity offset vectors are
                not one-dimensional, if the flux array is not
                two-dimensional, if the inverse variance or mask arrays
                do not have the same shape as the flux array, or if the
                number of wavelengths does not match the second axis of
                the flux array.
        """
        # Check the input
        if len(wave.shape) != 1:
            raise ValueError('Input wavelength array must be one-dimensional.')
        if len(flux.shape) != 2:
            raise ValueError('Input flux array must be two-dimensional.  To register a single ' \
                             'flux vector, use mangadap.util.sampling.Resample.')
        if flux.shape[1] != wave.size:
            raise ValueError('Flux array shape does not match the wavelength data.')
        if ivar is not None and ivar.shape != flux.shape:
            raise ValueError('Input inverse variance array must have the same shape as flux array.')
        if mask is not None and mask.shape != flux.shape:
            raise ValueError('Input mask array must have the same shape as flux array.')
        if sres is not None and sres.shape != wave.shape and sres.shape != flux.shape:
            raise ValueError('Input spectral resolution data has incorrect shape.')

        # Get the mask
        inp_mask = SpectralStack._get_input_mask(flux, ivar=ivar, mask=mask, dtype=float)

        # Get the spectral resolution (always a masked array)
        inp_sres = SpectralStack._check_input_sres(sres, flux.shape[0])

        # Output spectral range
        outRange = [wave[0], wave[-1]] if keep_range \
                        else list(SpectralStack.min_max_wave(wave, voff))
        # Sampling (logarithmic or linear)
        dw = spectral_coordinate_step(wave, log=log, base=base)

        # Resample the flux and error
        ferr = None if ivar is None else numpy.ma.power(ivar, -0.5)
        resamp = Resample(flux, e=ferr, mask=inp_mask, x=wave, inLog=log, newRange=outRange,
                          newdx=dw, base=base)

        # Mask pixels covered below the designate fraction
        _mask = resamp.outf < flim

        # Get the inverse variance
        if ivar is None:
            _ivar = None
        else:
            _ivar = numpy.ma.power(resamp.oute, -2)
            _ivar[_mask] = numpy.ma.masked

        # Interpolate the spectral resolution at the location of the new
        # wavelength vector
        if inp_sres is None:
            _sres = None
        else:
            interp = interpolate.interp1d(wave, inp_sres, axis=1, assume_sorted=True,
                                          fill_value='extrapolate')
            _sres = interp(resamp.outx)

        # Return the result
        return resamp.outx, numpy.ma.MaskedArray(resamp.outy, mask=_mask), _ivar, _sres

    @staticmethod
    def build_covariance_data(cube, covariance_mode, covariance_par):
        """
        Construct the covariance data needed by the specified
        covariance mode.

        Args:
            cube (:class:`mangadap.datacube.datacube.DataCube`):
                Datacube object to bin.
            covariance_mode (:obj:`str`):
                Mode to use. Must be an allowed mode; see
                :func:`covariance_mode_options`.
            covariance_par (:obj:`int`, :obj:`float`, :obj:`list`):
                Parameters required by the specified mode.

        Returns:
            None, :obj:`float`,
            :class:`mangadap.util.covariance.Covariance`: Covariance
            data relevant to the selected mode.

        Raises:
            ValueError:
                Raised if the covariance mode is not valid or if a
                covariance matrix cannot be computed for the provided
                datacube.
            TypeError:
                Raised if the covariance parameters do not have the
                correct type.
        """
        # Check that the covariance data is correct
        if covariance_mode is not None \
                and covariance_mode not in SpectralStack.covariance_mode_options():
            raise ValueError('Unrecognized method for covariance: {0}'.format(covariance_mode))
        if covariance_mode is None:
            covariance_mode = 'none'
        if covariance_mode == 'none':
            return None
        
        if covariance_mode == 'calibrate':
            return covariance_par

        # TODO: I don't think this catches every use case...
        if not cube.can_compute_covariance:
            raise ValueError('Cannot compute covariance for the provided datacube.')

        if covariance_mode in [ 'channels', 'wavelengths' ]:
            if covariance_mode == 'channels':
                if not isinstance(covariance_par, int):
                    raise TypeError('Unable to handle \'channels\' covariance parameter with ' \
                                    'type: {0}'.format(type(covariance_par)))
                _covariance_par = numpy.linspace(0,cube.nwave-1,num=covariance_par).astype(int)
            else:
                _covariance_par = covariance_par
                if isinstance(_covariance_par, float):
                    _covariance_par = numpy.array([covariance_par])
                if isinstance(_covariance_par, list):
                    _covariance_par = numpy.array(covariance_par)
                if not isinstance(_covariance_par, numpy.ndarray):
                    raise TypeError('Unable to handle covariance parameter of type: {0}'.format(
                                    type(covariance_par)))

                # Convert the wavelengths to channel numbers
                _covariance_par = numpy.unique(numpy.argsort(numpy.absolute(
                                                cube.wave[:,None]-_covariance_par[None,:]),
                                                             axis=1)[0,:])
            return cube.covariance_cube(channels=_covariance_par)

        if covariance_mode == 'approx_correlation':
            if not isinstance(covariance_par, float):
                raise TypeError('For approximate correlation matrix, provide sigma as a float.')
            return cube.covariance_cube(sigma_rho=covariance_par)

        if covariance_mode == 'full':
            warnings.warn('Sit back.  This is going to take a while...')
            return cube.covariance_cube()

        raise ValueError('Unrecognized covariance method: {0}'.format(covariance_mode))


    def stack_datacube(self, cube, binid, par=None):
        r"""
        Wrapper function for :func:`stack` that accepts a datacube
        object.

        .. todo::
            - List the required elements of ``par``.

        Args:
            cube (:class:`mangadap.datacube.datacube.DataCube`):
                Datacube to stack.
            binid (:obj:`numpy.ndarray`):
                Indices of the bin in which to place each spectrum.
                Shape must be the flattened spatial shape of the
                datacube; i.e., :math:`(N_{\rm spec},)`.
            par (:class:`SpectralStackPar`, optional):
                Set of parameters used to define how to stack the
                spectra. See :func:`stack`. Does not need to be
                provided on initialization, but a parameter set is
                required for all the stacking routines.

        Returns:
            :obj:`tuple`: See the return statement for :func:`stack`.
        """
        flux = cube.copy_to_masked_array(flag=cube.do_not_stack_flags())
        ivar = cube.copy_to_masked_array(attr='ivar', flag=cube.do_not_stack_flags())
        sres = cube.copy_to_array(attr='sres')
        covar = SpectralStack.build_covariance_data(cube, par['covar_mode'], par['covar_par'])

        # Make sure all the inverse variance values are valid
        indx = numpy.invert(ivar > 0)
        if numpy.any(indx):
            flux.mask |= indx
            ivar.mask |= indx

        return self.stack(cube.wave, flux, binid=binid, ivar=ivar, log=True, keep_range=True) \
                    if par is None else \
                    self.stack(cube.wave, flux, operation=par['operation'], binid=binid, ivar=ivar,
                               sres=sres, voff=par['vel_offsets'], log=True,
                               covariance_mode=par['covar_mode'], covar=covar, keep_range=True)

    def stack(self, wave, flux, operation='mean', binid=None, binwgt=None, ivar=None, mask=None,
              sres=None, voff=None, log=False, base=10.0, covariance_mode=None, covar=None,
              keep_range=False):
        r"""
        Stack a set of spectra.

        The method also sets all the returned variables to the
        internal attributes; however, they are always kept as the
        sum, not the mean, of the spectra, regardless of the value of
        ``operation``.

        .. todo:
            - Allow for renormalization of spectra before stacking.
              Where and how, TBD.
            - If only one spectrum returned, return flux array as a
              single vector?

        Args:
            wave (`numpy.ndarray`_):
                Single wavelength vector for all input spectra.
            flux (`numpy.ndarray`_, `numpy.ma.MaskedArray`_):
                Spectrum flux values. Shape must be :math:`(N_{\rm
                spec}, N_{\rm wave})`.
            operation (:obj:`str`, optional):
                Stacking operation to perform. See
                :func:`operation_options`.
            binid (`numpy.ndarray`_, optional):
                Indices of the bin in which to place each spectrum.
                Shape is :math:`(N_{\rm spec},)`. If not provided,
                all the spectra will be combined.
            binwgt (`numpy.ndarray`_, optional):
                Weights for each of the spectra. Shape is
                :math:`(N_{\rm spec},)`. If None, weights are
                uniform.
            ivar (`numpy.ndarray`_, `numpy.ma.MaskedArray`_, optional):
                Inverse variance in the spectrum fluxes. Shape must
                match ``flux``.
            mask (`numpy.ndarray`_, optional):
                Binary mask values for the spectrum fluxes
                (False=unmasked; True=masked). Shape must match
                ``flux``. Default assumes no pixel mask.
            sres (`numpy.ndarray`_, `numpy.ma.MaskedArray`_, optional):
                1D or 2D spectral resolution as a function of
                wavelength for all or each input spectrum. Shape must
                be :math:`(N_{\rm wave},)` or :math:`(N_{\rm
                spec},N_{\rm wave})`. If provided as a masked array,
                masked pixels are replaced with the interpolated
                resolution at that location.
            voff (`numpy.ndarray`_, optional):
                Vector with velocity offsets to apply to each
                spectrum before stacking. See :func:`register`.
            log (:obj:`bool`, optional):
                Flag that the wavelength vector is geometrically stepped
                in wavelength.
            base (:obj:`float`, optional):
                If the wavelength vector is geometrically stepped, this
                is the base of the logarithmic step.
            covariance_mode (:obj:`str`, optional):
                Keyword for method to use for dealing with covariance;
                see :func:`covariance_mode_options`.
            covar (:obj:`float`, :class:`mangadap.util.covariance.Covariance`, optional):
                Covariance object to use, which must match the
                expectation from the covariance mode.  See
                :func:`covariance_mode_options` and
                :func:`_check_covariance_type`.
            keep_range (:obj:`float`, optional):
                When registering the wavelengths of the shifted spectra,
                keep the identical spectral range as input.

        Returns:
            :obj:`tuple`: Returns seven objects: the wavelength
            vector, the stacked flux, the standard deviation in the
            stacked flux, the number of spectra in each stacked
            pixel, the stacked inverse variance, the stacked spectral
            resolution, and the stacked covariance.
            
        Raises:
            ValueError:
                Raised if the wavelength vector is not one-dimensional,
                if the flux array is not two-dimensional, if the inverse
                variance or mask arrays do not have the same shape as
                the flux array, or if the number of wavelengths does not
                match the second axis of the flux array.  Also raised if
                the covariance mode is not recognized; see
                :func:`covariance_mode_options`.
            TypeError:
                Raised if the object type of ``covar`` does not match
                what is expected given ``covariance_mode``.
            NotImplementedError:
                Raised if the requesting to both velocity register
                the spectra and correct the error vectors for spatial
                covariance.
        """
        # Check the input shapes
        if len(wave.shape) != 1:
            raise ValueError('Input wavelength vector must be one-dimensional!')
        nwave = wave.size
        if len(flux.shape) != 2:
            raise ValueError('Can only stack two-dimensional matrices.')
        if flux.shape[1] != nwave:
            raise ValueError('Flux array shape does not match the wavelength data.')
        if ivar is not None and flux.shape != ivar.shape:
            raise ValueError('Shape of the inverse-variance array must match the flux array.')
        if mask is not None and flux.shape != mask.shape:
            raise ValueError('Shape of the mask array must match the flux array.')
        if sres is not None and sres.shape != wave.shape and sres.shape != flux.shape:
            raise ValueError('Input spectral resolution data has incorrect shape.')

        nspec = flux.shape[0]
        if binid is not None and binid.size != nspec:
            raise ValueError('Length of binid must match the number of input spectra.')
        if binwgt is not None and binwgt.size != nspec:
            raise ValueError('Length of binwgt must match the number of input spectra.')

        # Check that the covariance data is correct
        if covariance_mode is not None \
                and covariance_mode not in SpectralStack.covariance_mode_options():
            raise ValueError('Unrecognized method for covariance: {0}'.format(covariance_mode))
        if covariance_mode is None:
            covariance_mode = 'none'
        if not SpectralStack._check_covariance_type(covariance_mode, covar, ivar):
            raise TypeError('Incorrect covariance and/or inverse variance object type for input ' \
                            'mode:\n mode: {0}\n input covar type: {1}\n input ivar' \
                            ' type: {2}'.format(covariance_mode, type(covar), type(ivar)))
        if not SpectralStack._check_covariance_shape(covariance_mode, covar, nwave, nspec):
            raise ValueError('Covariance object has incorrect shape for use with specified mode.')
        if isinstance(covar, Covariance) and voff is not None:
            raise NotImplementedError('Currently cannot both velocity register and apply ' \
                                      'covariance matrix calculation!')
        
        # Get the masked, velocity registered flux and inverse variance
        # arrays
        if voff is None:
            _mask = SpectralStack._get_input_mask(flux, ivar=ivar, mask=mask)
            _flux = numpy.ma.MaskedArray(flux, mask=_mask)
            _ivar = None if ivar is None else numpy.ma.MaskedArray(ivar, mask=_mask)
            _sres = SpectralStack._check_input_sres(sres, flux.shape[0])
            self.wave = wave
        else:
            self.wave, _flux, _ivar, _sres = register(wave, voff, flux, ivar=ivar, mask=mask,
                                                      sres=sres, log=log, base=base,
                                                      keep_range=keep_range)

        # Calculate the transfer matrix
        self._set_rebin_transfer_matrix(numpy.zeros(nspec, dtype=numpy.int) 
                                            if binid is None else binid, binwgt=binwgt)

        # Stack the spectra with or without covariance
        if covariance_mode == 'none':
            self._stack_without_covariance(_flux, ivar=_ivar, sres=_sres)
        else:
            self._stack_with_covariance(_flux, covariance_mode, covar, ivar=_ivar, sres=_sres)

        # Calculate the standard deviation in the flux, even if the flux
        # operation is to sum the data
        self.fluxmean = self.flux/self.npix
        self.fluxsdev = numpy.ma.sqrt((self.fluxsqr/self.npix - numpy.square(self.fluxmean))
                                            * self.npix*numpy.ma.power((self.npix-1), -1.))

        # Interpolate across any masked pixels in the spectral
        # resolution array, if provided and requested.
        if self.sres is not None:
            if numpy.sum(self.sres.mask) == 0:
                self.sres = self.sres.data
            else:
                outshape = self.sres.shape
                self.sres = numpy.apply_along_axis(interpolate_masked_vector, 1,
                                                   self.sres.reshape(1,-1) if self.sres.ndim == 1
                                                        else self.sres.reshape(outshape[0], -1)
                                                   ).reshape(outshape)

        # Return the stacked data
        if operation == 'sum':
            return self.wave, self.flux, self.fluxsdev, self.npix, self.ivar, self.sres, self.covar
        return self._get_stack_mean()


        
