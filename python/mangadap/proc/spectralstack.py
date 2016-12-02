# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""

Stack some spectra!

*License*:
    Copyright (c) 2015, SDSS-IV/MaNGA Pipeline Group
        Licensed under BSD 3-clause license - see LICENSE.rst

*Source location*:
    $MANGADAP_DIR/python/mangadap/proc/spatialbins.py

*Imports and python version compliance*:
    ::

        from __future__ import division
        from __future__ import print_function
        from __future__ import absolute_import
        from __future__ import unicode_literals

        import sys
        if sys.version > '3':
            long = int
        
        import numpy
        from scipy import sparse
        from astropy.io import fits
        import astropy.constants

        from ..par.parset import ParSet
        from ..util.covariance import Covariance
        from ..util.misc import inverse_with_zeros

*Class usage examples*:

    .. todo::
        Add examples

*Revision history*:
    | **01 Apr 2016**: Implementation begun by K. Westfall (KBW)
    | **10 Nov 2016**: (KBW) Include dispersion cube in DRP file to
        construct spectral resolution of binned spectra.
    | **30 Nov 2016**: (KBW) Use, e.g., [None,:] slicing instead of,
        e.g., numpy.array([v]*n) when necessary for array arithmetic
    | **01 Dec 2016**: (KBW) Allow stacking the spectral resolution to
        be turned off via :class:`SpectralStackPar`.

.. _astropy.io.fits.hdu.hdulist.HDUList: http://docs.astropy.org/en/v1.0.2/io/fits/api/hdulists.html
.. _astropy.io.fits.Header: http://docs.astropy.org/en/stable/io/fits/api/headers.html#header
.. _numpy.ma.MaskedArray: http://docs.scipy.org/doc/numpy-1.10.1/reference/maskedarray.baseclass.html#numpy.ma.MaskedArray

"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import sys
if sys.version > '3':
    long = int

import numpy
from scipy import sparse
from astropy.io import fits
import astropy.constants

from ..par.parset import ParSet
from ..util.covariance import Covariance
from ..util.misc import inverse_with_zeros
from ..util.constants import constants

from matplotlib import pyplot

__author__ = 'Kyle B. Westfall'
# Add strict versioning
# from distutils.version import StrictVersion

class SpectralStackPar(ParSet):
    r"""
    Class with parameters used to set how to stack a set of spectra.
    See :class:`mangadap.par.parset.ParSet` for attributes.

    .. todo::
        Allow for sigma rejection.

    Args:
        operation (str): Operation to perform for the stacked spectrum.
            See :func:`SpectralStack.operation_options` for the
            available operation options.
        vel_register (bool): Flag to velocity register the spectra
            before adding them based on a provided prior measurement of
            the velocities.
        vel_offsets (list, numpy.ndarray): List of velocity offsets to
            apply to the spectra to stack.
        covar_mode (str): Describes how to incorporate covariance into
            the spectral stacking.  See
            :func:`SpectralStack.covariance_mode_options` for the
            available options.
        covar_par (int, float, numpy.ndarray, list): The parameter(s)
            needed to perform a given method of handling the covariance.
            See :func:`SpectralStack.covariance_mode_options` for the
            available options.

    """
    def __init__(self, operation, vel_register, vel_offsets, covar_mode, covar_par, stack_sres):
        in_fl = [ int, float ]
        ar_like = [ numpy.ndarray, list ]
        op_options = SpectralStack.operation_options()
        covar_options = SpectralStack.covariance_mode_options()
        
        pars =     [ 'operation', 'vel_register', 'vel_offsets',  'covar_mode',   'covar_par',
                        'stack_sres' ]
        values =   [   operation,   vel_register,   vel_offsets,    covar_mode,     covar_par,
                          stack_sres ]
        defaults = [      'mean',          False,          None,        'none',          None,
                                True ]
        options =  [  op_options,           None,          None, covar_options,          None,
                                None ]
        dtypes =   [         str,           bool,       ar_like,           str, in_fl+ar_like,
                                bool ]

        ParSet.__init__(self, pars, values=values, defaults=defaults, options=options,
                        dtypes=dtypes)


    def toheader(self, hdr):
        """
        Copy the information to a header.

        Args:
            hdr (`astropy.io.fits.Header`_): Header object to write to.
        Returns:
            `astropy.io.fits.Header`_ : Edited header object.
        """
        hdr['STCKOP'] = (self['operation'], 'Stacking operation')
        hdr['STCKVREG'] = (str(self['vel_register']), 'Spectra shifted in velocity before stacked')
        hdr['STCKCRMD'] = (str(self['covar_mode']), 'Stacking treatment of covariance')
        hdr['STCKCRPR'] = (str(self['covar_par']), 'Covariance parameter(s)')
        hdr['STCKRES'] = (str(self['stack_sres']), 'Spectral resolution stacked')
        return hdr


    def fromheader(self, hdr):
        """
        Copy the information from the header

        hdr (`astropy.io.fits.Header`_): Header object to write to.
        """
        self['operation'] = hdr['STCKOP']
        self['vel_register'] = bool(hdr['STCKVREG'])
        self['covar_mode'] = hdr['STCKCRMD']
        self['covar_par'] = eval(hdr['STCKCRPR'])
        self['stack_sres'] = eval(hdr['STCKRES'])


class SpectralStack():
    r"""
    Class whose primary function is to stack a set of spectra while
    treating covariance between spectra.

    See :func:`covariance_mode_options` for available methods of
    accounting for covariance.

    Class also approximates the resulting spectral resolution of the
    stacked data.

    """
    def __init__(self):
        # Keep the matrix used to bin the spectra
        self.rebin_T = None

        # Internal arrays for callback
        self.wave = None
        self.flux = None
        self.fluxsqr = None
        self.flux_sdev = None
        self.npix = None
        self.ivar = None
        self.sres = None
        self.covar = None


    @staticmethod
    def _check_covariance_type(covariance_mode, covar, ivar):
        """
        Check that the covariance variable has the correct type for the given mode.

        Args:
            covariance_mode (str): Covariance handling mode; see
                :func:`covariance_mode_options`.

            covar (None, float,
                :class:`mangadap.util.covariance.Covariance`): The
                object to check the type against the covariance handling
                mode.

        Returns:
            bool: Flag that type is correct.
        """
        if covariance_mode == 'none':
            return True
        if covariance_mode in [ 'calibrate', 'channels', 'wavelengths' ] and ivar is None:
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

        Check that the input covariance object has the correct shape for
        the given mode.

        Args:

        """
        if covariance_mode in [ 'none', 'calibrate']:
            return True
        if covariance_mode in [ 'full', 'approx_correlation' ] \
                and (covar.dim != 3 or covar.shape[0] != nwave):
            return False
        if covar.shape[-1] != nspec:
            return False
        return True


    @staticmethod
    def _get_input_mask(flux, ivar=None, mask=None, dtype=bool):
        inp_mask = numpy.full(flux.shape, False, dtype=bool) if mask is None else mask
        if isinstance(flux, numpy.ma.MaskedArray):
            inp_mask |= numpy.ma.getmaskarray(flux)
        if ivar is not None:
            inp_mask |= ~(ivar>0)
        return inp_mask.astype(dtype)

     
    @staticmethod
    def _get_input_sres(sres, nspec):
        if sres is None:
            return None
        if len(sres.shape) == 2 and nspec > 0 and sres.shape[0] != nspec:
            raise ValueError('spectral resolution array is not the correct shape.')
        if (len(sres.shape) == 2 and sres.shape[0] == nspec) or nspec == 0:
            return sres
        return numpy.array([sres]*nspec)


    def _set_rebin_transfer_matrix(self, binid, binwgt=None):
        r"""
        Construct the transfer matrix that rebins the spectra.  The
        output shape is :math:`(N_{\rm bin} \times N_{\rm spec})` with
        the expectation that the spectrum flux array has shape
        :math:`(N_{\rm spec} \times N_{\rm wave})`.  The binned spectra
        are calculated by matrix multiplication, :math:`\mathbf{B} =
        \mathbf{T} \times \mathbf{F}` such that the covariance matrix
        can be calculated as :math:`\mathbf{C} = \mathbf{T} \times
        \mathbf{\Sigma} \times \mathbf{T}^{\rm T}`, where
        :math:`\mathbf{\Sigma}` is the covariance matrix in the flux
        array, :math:`\mathbf{F}`.

        If weighting, the sum of the weights is normalized to the number
        of points included in the bin.

        Args:
            binid (numpy.ndarray): List if indices, one per spectrum in
                the flux array, for the binned spectrum.  Indices of
                less than one are ignored.
            binwgt (numpy.ndarray): (**Optional**) List of weights for
                the spectra.  If not provided, the weights are uniform.
        """
        nspec = binid.size
        valid = binid > -1
        unique_bins = numpy.unique(binid[valid])
        nbin = numpy.amax(unique_bins)+1        # Allow for missing bin numbers

        self.rebin_T = numpy.zeros((nbin,nspec), dtype=numpy.float)
        for j in range(nbin):
            indx = binid == j
            self.rebin_T[j,indx] = 1.0 if binwgt is None else \
                                    binwgt[indx]*numpy.sum(indx)/numpy.sum(binwgt[indx])
        self.rebin_T = sparse.csr_matrix(self.rebin_T)


    def _stack_without_covariance(self, flux, ivar=None, sres=None):
        """

        Stack the spectra, ignoring any covariance that may or may not
        exist.

        Sets :attr:`flux`, :attr:`fluxsqr`, :attr:`npix`, :attr:`ivar`, :attr:`sres`.

        The stored data is always based on the SUM of the spectra in the
        stack.
        
        Args:
            flux (numpy.ma.MaskedArray): Flux array.
            ivar (numpy.ma.MaskedArray): (**Optional**) The inverse
                variance array.
            sres (numpy.ndarray): (**Optional**) 1D or 2D spectral
                resolution as a function of wavelength for all or each
                input spectrum.  Default is to ignore any spectral
                resolution data.
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

        if sres is not None:
            Tc = rt.sum(axis=1).ravel()
            Tc[~(Tc>0)] = 1.0
            self.sres = numpy.ma.power(numpy.ma.dot(rt, numpy.ma.power(sres, -2))/Tc[:,None], -0.5)

        if ivar is None:
            return

        # No covariance so:
        self.ivar = numpy.ma.power( numpy.ma.dot(numpy.power(rt, 2.0), numpy.ma.power(ivar, -1.)),
                                   -1.)
        
#        pyplot.plot(self.wave, ivar[20*44+20,:])
#        pyplot.plot(self.wave, self.ivar[0,:])
#        pyplot.show()


    def _stack_with_covariance(self, flux, covariance_mode, covar, ivar=None, sres=None):
        """
        Stack the spectra and incorporate covariance.
        
        Return the
        sum of the flux, the sum of the flux squared, the number of
        points in each sum, and the covariance.

        ivar must not be none if covariance_mode is channels or wavelengths

        Size has to match self.rebin_T
        """
        # If calibrating, first stack without covariance, and the
        # renormalize the error.  Covar is expected to be a single float
        self._stack_without_covariance(flux, ivar=ivar, sres=sres)

        if covariance_mode == 'calibrate':
            if ivar is not None:
                self.ivar /= numpy.square(1.0 + covar*numpy.ma.log10(self.npix))
            return

        nwave = flux.shape[1]
        if covariance_mode in [ 'channels', 'wavelengths' ]:
            nchan = covar.shape[0]
            nbin = self.flux.shape[0]
            self.covar = numpy.empty(nchan, dtype=sparse.csr.csr_matrix)
            variance_ratio = numpy.ma.zeros( (nbin,nchan), dtype=numpy.float)
            for i in range(nchan):
                j = covar.input_indx[i]
                cov_array = covar._with_lower_triangle(plane=j)
                self.covar[i] = sparse.triu(self.rebin_T.dot(
                                             cov_array.dot(self.rebin_T.T))).tocsr()
#                pyplot.scatter(numpy.sqrt(self.covar[i].diagonal()), numpy.sqrt(1./self.ivar[:,j]),
#                               marker='.', s=30, lw=0, color='k')
                variance_ratio[:,i] = self.covar[i].diagonal() * self.ivar[:,j]
#                variance_ratio.mask[:,i] = self.ivar.mask[:,j]
                variance_ratio[numpy.ma.getmaskarray(self.ivar)[:,j],i] = numpy.ma.masked
#            pyplot.show()

#            print(numpy.sum(variance_ratio.mask))
#            for i in range(nbin):
##                pyplot.scatter(numpy.arange(nchan), variance_ratio[i,:], marker='.', s=30, lw=0)
#                pyplot.plot(numpy.arange(nchan), variance_ratio[i,:], lw=0.5)
#            pyplot.xlabel(r'Covariance Channel')
#            pyplot.ylabel(r'$C_{ii}\ I_{ii}$')
#            pyplot.show()
#            exit()

            # Below commented because there's a better way to do this
            # with, e.g., [:,None] indexing
#            ratio = numpy.array([numpy.ma.mean( variance_ratio, axis=1 )]*nwave).T
#            ratio = numpy.array([numpy.ma.median( variance_ratio, axis=1 )]*nwave).T
#            self.ivar = (numpy.ma.power(_ratio, -1.0).ravel()*self.ivar.ravel()).reshape(-1,nwave)
            ratio = numpy.ma.median( variance_ratio, axis=1 )
            self.ivar = numpy.ma.power(ratio, -1.0)[:,None] * self.ivar
            self.covar = Covariance(inp=self.covar, input_indx=covar.input_indx)
#            j = self.covar.input_indx[0]
#            self.covar.show(plane=j)
            return

        if covariance_mode in [ 'approx_correlation', 'full' ]:
            self.covar = numpy.empty(nwave, dtype=sparse.csr.csr_matrix)
            for i in range(nwave):
                self.covar[i] = sparse.triu(self.rebin_T.dot(
                                             covar[i].cov.dot(self.rebin_T.T))).tocsr()
            self.covar = Covariance(inp=self.covar)
            return

        raise ValueError('Unknown covariance mode: {0}'.format(covariance_mode))


    def _covar_in_mean(self):
        """

        Compute the covariance in the mean spectrum by propagating the
        division by the number of pixels through the covariance matrix.
       
        """
        if self.covar is None:
            return None
        
        nchan = self.covar.shape[0]
        nbin = self.flux.shape[0]
        inpix = numpy.ma.power(self.npix, -1.)
        covar = numpy.empty(nchan, dtype=sparse.csr.csr_matrix)
        for i in range(nchan):
            j = self.covar.input_indx[i]
            _inpix = numpy.ma.MaskedArray( [inpix[:,j]]*nbin ).T
            _inpix = (_inpix.ravel() * _inpix.T.ravel()).reshape(nbin,nbin)
#            pyplot.imshow(_inpix, origin='lower', interpolation='nearest')
#            pyplot.colorbar()
#            pyplot.show()
#            self.covar.show(plane=j)
#            pyplot.imshow((self.covar.toarray(plane=j).ravel()*_inpix.ravel()).reshape(nbin,nbin),
#                          origin='lower', interpolation='nearest')
#            pyplot.colorbar()
#            pyplot.show()
            covar[i] = sparse.triu((self.covar.toarray(plane=j).ravel()
                                            * _inpix.ravel()).reshape(nbin,nbin)).tocsr()
#            pyplot.imshow(covar[i].toarray(), origin='lower', interpolation='nearest')
#            pyplot.colorbar()
#            pyplot.show()
        return Covariance(inp=covar, input_indx=self.covar.input_indx)
            

    @staticmethod
    def operation_options():
        """
        Return the allowed stacking operations.  Current operations are:
        
            ``mean``: Construct the mean of the spectra

            ``sum``: Construct the spectrum sum.

        Returns:
            list: List of available operations.
        """
        return ['mean', 'sum']


    @staticmethod
    def covariance_mode_options(par_needed=False):
        r"""
        Accounting for covariance:  The two parameters covar_mode and
        covar_par set how covariance is accounted for in the stacking
        procedure.  The valid options are:

            ``none``: The noise in the stacked spectrum is a nominal
            propagation of the error assuming no covariance.  No
            parameters needed.

            ``calibrate``: The spectral noise is calibrated following:

            .. math::

                n_{\rm calib} = n_{\rm nominal} (1 + \alpha \log\
                N_{\rm bin})

            where :math:`N_{\rm bin}` is the number of binned spaxels.
            The value of :math:`\alpha` must be provided as a parameter.
     
            ``channels``: The noise vector of each stacked spectrum is
            adjusted based on the mean ratio of the nominal and formally
            correct calculations of the noise measurements over a number
            of spectral channels.  The channels are drawn from across
            the full spectral range.  The number of channels to use is a
            defined parameter.  The covariance matrix must be provided
            to :func:`stack`.

            ``wavelengths``: Functionally equivalent to ``channels``;
            however, the channels to use is set by a list of provided
            wavelengths.  The covariance matrix must be provided to
            :func:`stack`.

            ``approx_correlation``: Approximate the covariance matrix
            using a Gaussian description of the correlation between
            pixels.  See
            :func:`mangadap.drpfits.DRPFits.covariance_matrix`.  The
            value of :math:`\sigma` provides for the Gaussian desciption
            of :math:`\rho_{ij}` in the correlation matrix.  The
            covariance matrix must be provided to :func:`stack`.

            ``full``: The full covariance cube is calculated and the
            noise vectors are constructed using the formally correct
            calculation.  No parameters needed.  The covariance matrix
            must be provided to :func:`stack`.

        Returns:
            list: List of the allowed modes.
        """
        modes = [ 'calibrate', 'approx_correlation', 'channels', 'wavelengths' ]
        if par_needed:
            return modes
        return modes + [ 'none', 'full' ]


    @staticmethod
    def parse_covariance_parameters(mode, par):
        """
        Parse the parameters needed for the treatment of the covariance when
        stacking spectra.
    
        Args:
            mode (str): Mode to use.  Must be an allowed mode; see
                :func:`covariance_mode_options`.
            par (str): String representation of the parameters for the
                specified mode.

        Returns:
            float or list: Parameters parsed from the input string for
            the designated covariance mode.

        Raises:
            TypeError: Raised if the input parameter could not be
                converted to a float as needed by the specified mode.
            ValueError: Raised if the mode is not recognized.
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
        """
        Determine the minimum and maximum of all shifted wavelength
        ranges.

        Args:
            wave (numpy.ndarray): Original wavelengths.

            voff (numpy.ndarray): Velocity offsets in km/s.  Each
                element is applied to the wavelength vector to determine
                the maximum wavelength range.  Does not need to be
                one-dimensional.

        Returns:
            float: Two floats with the minimum and maximum redshifted
            wavelengths.
        """
        _wave = numpy.array([numpy.amin(wave), numpy.amax(wave)])
        redshifted = _wave[None,:]*(1.+voff[:,None]/astropy.constants.c.to('km/s').value)
        return numpy.amin(redshifted), numpy.amax(redshifted)


    #TODO: Untested!
    @staticmethod
    def register(wave, voff, flux, ivar=None, mask=None, sres=None, log=False, base=10.0,
                 keep_range=False):
        """

        Register a set of spectra to the same wavelength range given a
        set of velocity offsets.

        .. warning::
            **THIS FUNCTION IS UNTESTED!**

        Args:
            wave (numpy.ndarray): Single wavelength vector for all input
                spectra.
            voff (numpy.ndarray): Vector with velocity offsets to apply
                to each spectrum.
            flux (numpy.ndarray): Spectrum flux values.  Can be a
                masked array.
            ivar (numpy.ndarray): (**Optional**) Inverse variance in the
                spectrum fluxes.  Can be a masked array.
            mask (numpy.ndarray): (**Optional**) Binary mask values for
                the spectrum fluxes; 0 (False) is unmasked, anything
                else is masked.  Default assumes no pixel mask.
            sres (numpy.ndarray): (**Optional**) 1D or 2D spectral
                resolution as a function of wavelength for all or each
                input spectrum.  Default is to ignore any spectral
                resolution data.
            log (bool): (**Optional**) Flag that the wavelength vector
                is geometrically stepped in wavelength.
            base (float): (**Optional**) If the wavelength vector is
                geometrically stepped, this is the base of the
                logarithmic step.  Default is 10.0.
            keep_range (float): (**Optional**)  When registering the
                wavelengths of the shifted spectra, keep the identical
                spectral range as input.

        Returns:
            numpy.ndarray, `numpy.ma.MaskedArray`_: Returns four arrays:
            (1) the new wavelength vector, common to all spectra; (2)
            the new flux array; (3) the new inverse variance array,
            which will be None if no inverse variances are provided to
            the function; and (4) the new spectral resolution vectors,
            which will also be None if no spectral resolution vectors
            are provided.

        Raises:
            ValueError: Raised if the wavelength or velocity offset
                vectors are not one-dimensional, if the flux array is
                not two-dimensional, if the inverse variance or mask
                arrays do not have the same shape as the flux array, or
                if the number of wavelengths does not match the second
                axis of the flux array.
        """
        # Check the input
        if len(wave.shape) != 1:
            raise ValueError('Input wavelength array must be one-dimensional.')
        if len(flux.shape) != 2:
            raise ValueError('Input flux array must be two-dimensional.  To register a single ' \
                             'flux vector, use mangadap.util.instrument.resample_vector.')
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

        # Get the spectral resolution
        inp_sres = SpectralStack._get_input_sres(sres, flux.shape[0])

        # Input and output spectral range
        inRange = [wave[0], wave[-1]]
        outRange = inRange if keep_range else list(min_max_wave(wave, voff))
        # Sampling (logarithmic or linear)
        dw = (numpy.log(wave[1]) - numpy.log(wave[0]))/numpy.log(base) \
                    if log else (wave[1] - wave[0])
        # Output number of pixels
        if keep_range:
            nwave = wave.size
        else:
            nwave, outRange = resample_vector_npix(outRange=outRange, dx=dw, log=log, base=base)
        
        #Initialize the output arrays
        _flux = numpy.empty((nspec,nwave), dtype=numpy.float)
        _ivar = numpy.zeros((nspec,nwave), dtype=numpy.float)
        _mask = numpy.zeros((nspec,nwave), dtype=numpy.float)
        _sres = None if sres is None else numpy.zeros((nspec,nwave), dtype=numpy.float)
        var = None if ivar is None else inverse_with_zeros(ivar, absolute=False)

        # Resample each spectrum
        for i in range(nspec):
            _wave, _flux[i,:] = resample_vector(flux[i,:].ravel(), xRange=inRange, inLog=log,
                                                newRange=outRange, newpix=nwave, base=base)
            _wave, _mask[i,:] = resample_vector(inp_mask[i,:].ravel(), xRange=inRange, inLog=log,
                                                newRange=outRange, newpix=nwave, base=base,
                                                ext_value=1)
            if var is not None:
                _wave, _ivar[i,:] = resample_vector(var[i,:].ravel(), xRange=inRange, inLog=log,
                                                    newRange=outRange, newpix=nwave, base=base)
            if inp_sres is not None:
                _wave, _sres[i,:] = resample_vector(inp_sres[i,:].ravel(), xRange=inRange,
                                                    inLog=log, newRange=outRange, newpix=nwave,
                                                    base=base)
        indx = _mask > 0.5
        _mask[indx] = 1.0
        _mask[numpy.invert(indx)] = 0.0
        _mask = _mask.astype(bool)

        return _wave, numpy.ma.MaskedArray(_flux, mask=_mask), \
               numpy.ma.MaskedArray(inverse_with_zeros(_ivar, absolute=False), mask=_mask), \
               _sres


    @staticmethod
    def build_covariance_data_DRPFits(drpf, covariance_mode, covariance_par):

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

        if covariance_mode in [ 'channels', 'wavelengths' ]:
            if covariance_mode == 'channels':
                if not isinstance(covariance_par, int):
                    raise TypeError('Unable to handle \'channels\' covariance parameter with ' \
                                    'type: {0}'.format(type(covariance_par)))
                _covariance_par = numpy.linspace(0,drpf.nwave-1,num=covariance_par).astype(int)
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
                                                drpf['WAVE'].data[:,None]-_covariance_par[None,:]),
                                                             axis=1)[0,:])
#                _wave = numpy.array([drpf['WAVE'].data]*len(_covariance_par)).T
#                _chan = numpy.array([_covariance_par]*drpf['WAVE'].data.size)
#                _covariance_par = numpy.unique(numpy.argsort(
#                                               numpy.absolute(_wave-_chan), axis=1)[0,:])

            return drpf.covariance_cube(channels=_covariance_par)

        if covariance_mode == 'approx_correlation':
            if not isinstance(covariance_par, float):
                raise TypeError('For approximate correlation matrix, provide sigma as a float.')
            return drpf.covariance_cube(sigma_rho=covariance_par)

        if covariance_mode == 'full':
            warnings.warn('Sit back.  This is going to take a while...')
            return drpf.covariance_cube()

        raise ValueError('Unrecognized covariance method: {0}'.format(covariance_mode))


    def stack_DRPFits(self, drpf, binid, par=None):
        """
        Wrapper function for :func:`stack` that uses a DRPFits file.

        Args:
            par (ParSet or dict): (**Optional**) Set of parameters used
                to define how the stack the spectra.  See :func:`stack`.
                Does not need to be provided on initialization, but a
                parameter set is required for all the stacking routines.

        Returns:
            numpy.ndarray, :class:`mangadap.util.covariance.Covariance`:
            Returns six elements.  See :func:`stack`.

        """
        wave = drpf['WAVE'].data
        flux = drpf.copy_to_masked_array(flag=drpf.do_not_stack_flags())
        ivar = drpf.copy_to_masked_array(ext='IVAR', flag=drpf.do_not_stack_flags())
        sres = None if par is None or not par['stack_sres'] \
                else drpf.spectral_resolution(toarray=True)
        covar = None if par is None else \
                    self.build_covariance_data_DRPFits(drpf, par['covar_mode'], par['covar_par'])

        return self.stack(wave, flux, binid=binid, ivar=ivar, sres=sres, log=True,
                          keep_range=True) if par is None else \
                    self.stack(wave, flux, operation=par['operation'], binid=binid, ivar=ivar,
                               sres=sres, voff=par['vel_offsets'], log=True,
                               covariance_mode=par['covar_mode'], covar=covar, keep_range=True)


    def stack(self, wave, flux, operation='mean', binid=None, binwgt=None, ivar=None, mask=None,
              sres=None, voff=None, log=False, base=10.0, covariance_mode=None, covar=None,
              keep_range=False):
        """
        Stack a set of spectra.  If binid is None, all the spectra in
        the array are stacked into a single output spectrum.

        Register a set of spectra to the same wavelength range given a
        set of velocity offsets.

        Args:
            wave (numpy.ndarray): Single wavelength vector for all input
                spectra.
            flux (numpy.ndarray): Spectrum flux values.  Can be a
                masked array.
            operation (str): (**Optional**) Stacking operation to
                perform.  See :func:`operation_options`.  Default is
                ``mean``.
            binid (numpy.ndarray): (**Optional**) Indices of the bin in
                which to place each spectrum.  If not provided, all the
                spectra will be combined.
            binwgt (numpy.ndarray): (**Optional**) Weights for each of
                the spectra.  If not provided, all weights are uniform.
            ivar (numpy.ndarray): (**Optional**) Inverse variance in the
                spectrum fluxes.  Can be a masked array.
            mask (numpy.ndarray): (**Optional**) Binary mask values for
                the spectrum fluxes; 0 (False) is unmasked, anything
                else is masked.  Default assumes no pixel mask.
            sres (numpy.ndarray): (**Optional**) 1D or 2D spectral
                resolution as a function of wavelength for all or each
                input spectrum.  Default is to ignore any spectral
                resolution data.
            voff (numpy.ndarray): (**Optional**) Vector with velocity
                offsets to apply to each spectrum.  Default is no
                velocity offsets
            log (bool): (**Optional**) Flag that the wavelength vector
                is geometrically stepped in wavelength.
            base (float): (**Optional**) If the wavelength vector is
                geometrically stepped, this is the base of the
                logarithmic step.  Default is 10.0.
            covariance_mode (str): (**Optional**) Keyword for method to
                use for dealing with covariance; see
                :func:`covariance_mode_options`.  Default is to ignore
                covariance.
            covar (None, float, Covariance): (**Optional**) Covariance
                object to use, which must match the expectation from the
                covariance mode.  See :func:`covariance_mode_options`
                and :func:`_check_covariance_type`.
            keep_range (float): (**Optional**)  When registering the
                wavelengths of the shifted spectra, keep the identical
                spectral range as input.

        Returns:
            numpy.ndarray, `numpy.ma.MaskedArray`_:
            

        return self.wave, self.flux, sdev, self.npix, self.ivar, self.covar

        Raises:
            ValueError: Raised if the wavelength vector is not
                one-dimensional, if the flux array is not
                two-dimensional, if the inverse variance or mask arrays
                do not have the same shape as the flux array, or if the
                number of wavelengths does not match the second axis of
                the flux array.  Also raised if the covariance mode is
                not recognized; see :func:`covariance_mode_options`.
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
            _sres = SpectralStack._get_input_sres(sres, flux.shape[0])
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

#        pyplot.plot(wave, ivar[20*44+20,:])
#        pyplot.plot(self.wave, self.ivar[0,:])
#        pyplot.show()

        # Calculate the standard deviation in the flux, even if the flux
        # operation is to sum the data
        mean = self.flux/self.npix
        sdev = numpy.ma.sqrt((self.fluxsqr/self.npix - numpy.square(mean))
                                * self.npix*numpy.ma.power((self.npix-1), -1.))

        # If summing, then the stacking procedure is done
        if operation == 'sum':
            return self.wave, self.flux, sdev, self.npix, self.ivar, self.sres, self.covar

        # If stacking to the mean, calculate the inverse variance of the
        # mean
        covar = self._covar_in_mean()
#        print(covar.input_indx)
#        covar.show(plane=0)

        return self.wave, mean, sdev, self.npix, self.ivar * numpy.square(self.npix), \
                    self.sres, covar


        
