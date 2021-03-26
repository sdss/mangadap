"""
Implement the derived class for MaNGA datacubes.

----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import time
import os
import warnings

from IPython import embed

import numpy
import scipy
import astropy
from matplotlib import pyplot


from astropy.io import fits
from astropy.wcs import WCS

from pydl.goddard.astro import airtovac
from ..util import sampling

from ..config import defaults
from ..util.drpfits import DRPFits, DRPFitsBitMask
from ..util.parser import DefaultConfig
from ..util.constants import DAPConstants
from ..util.covariance import Covariance
from ..spectra import MaNGARSS
from .datacube import DataCube

class MUSEDataCube(DataCube):
    r"""
    Container class for a MUSE datacube.  Copied by KHRR from MaNGADataCube

    For additional description and attributes, see the two base
    classes.

    :func:`~mangadap.datacube.datacube.DataCube.from_config` to
    instantiate from a configuration file.

    Note that ``z``, ``vdisp``, ``ell``, ``pa``, and ``reff`` are
    saved to the :attr:`~mangadap.datacube.datacube.DataCube.meta`
    dictionary.

    Args:
        ifile (:obj:`str`):
            The file with the MUSE datacube. The name is expected to
            follow the nominal naming convention, which is used to
            parse the plate, ifudesign, and whether or not the
            datacube has been binned logarithmically in wavelength.
        z (:obj:`float`, optional):
            Estimated bulk redshift. If None, some of the DAP
            analysis modules will fault.
        vdisp (:obj:`float`, optional):
            Estimated velocity dispersion. If None, some of the DAP
            analysis modules will assume an initial guess of 100
            km/s.
        ell (:obj:`float`, optional):
            Characteristic isophotal ellipticity (1-b/a). If None,
            some of the DAP modules will issue a warning and continue
            by assuming ``ell=0``.
        pa (:obj:`float`, optional):
            Characteristic isophotal position angle (through E from
            N). If None, some of the DAP modules will issue a warning
            and continue by assuming ``pa=0``.
        reff (:obj:`float`, optional):
            Effective (half-light) radius in arcsec. If None, some of
            the DAP modules will issue a warning and continue by
            assuming ``reff=1``.
        sres_ext (:obj:`str`, optional):
            The extension to use when constructing the spectral
            resolution vectors. See
            :func:`~mangadap.util.drpfits.DRPFits.spectral_resolution`.
        sres_fill (:obj:`bool`, optional):
            Fill masked values by interpolation. Default is to leave
            masked pixels in returned array.
        covar_ext (:obj:`str`, optional):
            Extension to use as the single spatial correlation matrix
            for all wavelength channels, read from the DRP file. This does not
            work with MUSE yet.  For
            generating the covariance matrix directly for an
            arbitrary wavelength channel using the RSS file, see
            :func:`mangadap.datacube.datacube.DataCube.covariance_matrix`.
    """
    def __init__(self, ifile, objra=None, objdec=None, z=None, vdisp=None, ell=None, \
                 pa=None, reff=None, sres_ifile=None,
                 sres_fill=True, covar_ext=None, plate=None, ifudesign=None, ebvgal=None):


        if not os.path.isfile(ifile):
            raise FileNotFoundError('File does not exist: {0}'.format(ifile))

        self.filename = os.path.abspath(ifile)
        self.plate = plate
        self.ifudesign = ifudesign
        self.drpver = '0v0'

        # Collect the metadata into a dictionary
        meta = {}
        meta['z'] = z
        meta['vdisp'] = vdisp
        meta['ell'] = ell
        meta['pa'] = pa
        meta['reff'] = reff
        meta['OBJRA'] = objra
        meta['OBJDEC'] = objdec


        with fits.open(ifile) as hdu:
            # Read covariance first
            #covar = None if covar_ext is None \
            #            else Covariance.from_fits(hdu, ivar_ext=None, covar_ext=covar_ext,
            #                                      impose_triu=True, correlation=True)

            covar = None

            print('Reading MUSE datacube data ...', end='\r')
            prihdr = hdu[0].header
            datahdr = hdu['DATA'].header
            flux = 10.0**(-3.0) * hdu['DATA'].data   ##  this is in units of 10**(-20)*erg/s/cm**2/Angstrom
                                                     ##  correct to units of 10**(-17)*erg/s/cm**2/Angstrom

            variance = 10.0**(-6.0) * hdu['STAT'].data  ## correct to units of 10**(-17)*erg/s/cm**2/Angstrom
            ivar = 1.0/variance
            error = variance**0.5

            ## Add objra and objdec to prihdr
            #prihdr['OBJRA'] = objra
            #prihdr['OBJDEC'] = objdec
            prihdr['EBVGAL'] = ebvgal


            datadim = numpy.shape(flux)
            spatial_shape = datadim[1:]
            nbins = numpy.prod(spatial_shape)

            # Define mask
            mask = numpy.full(datadim, False, dtype=bool)
            mask[numpy.logical_not(numpy.isfinite(ivar)) | numpy.logical_not(numpy.isfinite(flux))] = True

            # Read in linear air wavelength array and convert to vacuum
            wave = datahdr['CRVAL3']+(numpy.arange(datadim[0]))*datahdr['CD3_3']
            wave = airtovac(wave)

            # Read in spectral resolution
            # for specres, we want 1sigma LSF in Angstroms
            # MUSE sres file gives FWHM LSF in Angstroms
            lsf = numpy.genfromtxt(sres_ifile, comments='#')
            lsf_wave_air = lsf[:, 0]
            lsf_wave = airtovac(lsf_wave_air)
            lsf_sig = lsf[:, 1] / (2.0 * numpy.sqrt(2.0 * numpy.log(2.0)))

            # Interpolate onto wavelength array?
            # Convert lsf wavelength array into vacuum (done above), then interpolate using resampled cube
            sres_func = scipy.interpolate.interp1d(lsf_wave, lsf_sig, assume_sorted=True,
                                                   fill_value='extrapolate')

            # Number of angstroms per pixel
            ang_per_pix = numpy.array([sampling.angstroms_per_pixel(wave, log=False, regular=False)])
            fullRange = numpy.array([numpy.amin(wave), numpy.amax(wave)]).astype(float)

            # Compute desired spectral step
            l7400 = numpy.argmin(numpy.abs(wave-7400.0))
            dlogl = numpy.log10(wave[l7400]) - numpy.log10(wave[l7400-1])

            # Get the number of pixels needed
            npix, newRange = sampling.grid_npix(rng=fullRange, dx=dlogl, log=True)

            resampled_flux = numpy.zeros((npix, spatial_shape[0], spatial_shape[1]), dtype=numpy.float64)
            resampled_error = numpy.zeros((npix, spatial_shape[0], spatial_shape[1]), dtype=numpy.float64)
            resampled_wave = numpy.zeros((npix), dtype=numpy.float64)
            resampled_mask = numpy.zeros((npix, spatial_shape[0], spatial_shape[1]), dtype=bool)

            spatial_index = numpy.array([(i, j) for i, j in zip(*numpy.unravel_index(numpy.arange(nbins), spatial_shape))])

            # Loop through each bin
            for i in range(nbins):
                # Rebin the observed wavelength range
                # Do I need to add xBorders?  I guess not...
                # What do we do with the covar this produces?  -- Kyle says nothing

                r = sampling.Resample(flux[:,spatial_index[i][0],spatial_index[i][1]],
                                        e=error[:,spatial_index[i][0],spatial_index[i][1]],
                                        mask=mask[:,spatial_index[i][0],spatial_index[i][1]],
                                        x=wave, inLog=False,
                                        newRange=newRange, newLog=True,
                                        newdx=dlogl, step=True, covar=False)

                # Save the result
                if(i==0):
                    resampled_wave = r.outx

                resampled_flux[:,spatial_index[i][0],spatial_index[i][1]] = r.outy
                resampled_error[:,spatial_index[i][0],spatial_index[i][1]] = r.oute
                indx = r.outf < 0.9

                resampled_flux[:,spatial_index[i][0],spatial_index[i][1]][indx] = 0.0
                resampled_error[:,spatial_index[i][0],spatial_index[i][1]][indx] = 1.0
                resampled_mask[:,spatial_index[i][0],spatial_index[i][1]][indx] = True


            # Recreate ivar
            resampled_ivar = 1.0 / resampled_error**2

            # Resample the spectral resolution by simple interpolation.
            sres = numpy.ma.MaskedArray(sres_func(resampled_wave))

            # Convert to lambda / dlambda (as in DRPFits, with dlambda the FWHM of the resolution)
            sres = numpy.ma.divide(resampled_wave, sres) / DAPConstants.sig2fwhm

            # Set one corner to be masked to force first bin index to end up as -1
            # resampled_mask[:,spatial_index[0][0],spatial_index[0][1]] = True


            # NOTE: Need to keep a log of what spectral resolution
            # vectors were used so that the same vectors are read from
            # the RSS file, if/when it is loaded.
            #self.sres_ext, sres = DRPFits.spectral_resolution(hdu, ext=sres_ext, fill=sres_fill)
            #self.sres_fill = sres_fill
            #sres = sres.filled(0.0)


            # NOTE: Transposes are done here because of how the data is
            # read from the fits file. The covariance is NOT transposed
            # accordingly because of how the correlation matrices are
            # stored by the DRP. Some care needs to be taken here
            # though because the transpose changes the array memory
            # storage from C to Fortran contiguous. However, everything
            # is expected to be okay because numpy array-flattening
            # always performs a C-like flattening, even if the memory
            # storage is Fortran contiguous.

            # Is it a problem that the header now has the wrong wavelength info????

            DataCube.__init__(self, resampled_flux.T, wave=resampled_wave,
                              ivar=resampled_ivar.T, mask=resampled_mask.T, bitmask=None,
                              sres=sres, covar=None, wcs=WCS(header=datahdr, fix=True),
                              pixelscale=0.2, log=True, meta=meta, prihdr=prihdr,
                              fluxhdr=datahdr)
        print('Reading MUSE datacube data ... DONE')


    def file_path(self):

        return self.filename


    def mean_sky_coordinates(self, center_coo=None, offset='OBJ'):
        """
        Calculate the sky coordinates for each spaxel.

        This is a simple wrapper for
        :func:`mangadap.datacube.datacube.DataCube.mean_sky_coordinates`
        that passes the relevant set of center coordinates; see that
        function for the default behavior. Any provided coordinates
        (see ``center_coo``) take precedence.

        Args:
            center_coo (:obj:`tuple`, optional):
                A two-tuple with the coordinates in right-ascension
                and declination for the coordinate-frame origin. If
                None, no offset is performed.
            offset (:obj:`str`, optional):
                The offset coordinates to use from the primary
                header. Must be either 'obj' or 'ifu'. If None, no
                header coordinates are used.

        Returns:
            :obj:`tuple`: Two `numpy.ndarray`_ objects with the RA
            and declination of each pixel in degrees, or its offset
            from the center in arcseconds. In both cases the shape of
            the returned arrays matches the spatial shape of the
            datacube.
        """
        if center_coo is not None and offset is not None:
            warnings.warn('Center coordinates provided directly.  Ignore other offset requests.')
        if offset is not None and offset.upper() not in ['OBJ', 'IFU']:
            raise ValueError('Offset must be None, obj or ifu.')
        _center_coo = (self.meta['{0}RA'.format(offset.upper())],
                       self.meta['{0}DEC'.format(offset.upper())]) \
                           if center_coo is None and offset is not None else center_coo
        return super().mean_sky_coordinates(center_coo=_center_coo)


    @classmethod
    def from_config(cls, cfgfile, directory_path=None):
        """
        Construct a :class:`mangadap.datacube.muse.MUSEDataCube` object from a
        configuration file.

        The format of the configuration file is:

        .. todo::

            Fill this in.

        Args:
            cfgfile (:obj:`str`):
                Configuration file
            directory_path (:obj:`str`, optional):
                The exact path to the reduced MUSE data file. Overrides any value in the
                configuration file.
        """
        # Read the configuration file
        cfg = DefaultConfig(cfgfile, interpolate=True)

        cubefil = cfg.get('cubefil')
        sresfil = cfg.get('sresfil')

        # Set the attributes, forcing a known type
        plate = cfg.getint('plate')
        ifu = cfg.getint('ifu')
        if plate is None or ifu is None:
            raise ValueError('Configuration file must define the plate and IFU.')
        log = cfg.getbool('log', default=False)

        # Overwrite what's in the file with the method keyword arguments
        _directory_path = cfg.get('directory_path') if directory_path is None else directory_path

        ifile = _directory_path + cubefil

        # Get the other possible keywords
        # TODO: Come up with a better way to do this
        kwargs = {}
        kwargs['plate'] = plate
        kwargs['ifudesign'] = ifu
        kwargs['z'] = cfg.getfloat('z')
        kwargs['objra'] = cfg.getfloat('objra')
        kwargs['objdec'] = cfg.getfloat('objdec')
        kwargs['sres_ifile'] = sresfil
        kwargs['sres_fill'] = cfg.getbool('sres_fill', default=True)
        kwargs['covar_ext'] = cfg.get('covar_ext', default=None)
        kwargs['vdisp'] = cfg.getfloat('vdisp')
        kwargs['ell'] = cfg.getfloat('ell')
        kwargs['pa'] = cfg.getfloat('pa')
        kwargs['reff'] = cfg.getfloat('reff')
        kwargs['ebvgal'] = cfg.getfloat('ebvgal')

        return cls(ifile, **kwargs)