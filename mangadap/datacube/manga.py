"""
Implement the derived class for MaNGA datacubes.
"""


import os
import warnings

from IPython import embed

import numpy

from astropy.io import fits
from astropy.wcs import WCS

from ..config import defaults
from ..drpfits import DRPFitsBitMask
from ..util.constants import DAPConstants
from ..util.covariance import Covariance
from ..util.filter import interpolate_masked_vector
from .datacube import DataCube

class MaNGADataCube(DataCube):
    r"""
    Container class for a MaNGA datacube.

    For additional description and attributes, see
    :class:`mangadap.datacube.datacube.DataCube`.

    See :func:`from_plateifu` to instantiate using the plate and ifu
    numbers.

    Args:
        ifile (:obj:`str`):
            The file with the MaNGA datacube. The name is expected to
            follow the nominal naming convention, which is used to
            parse the plate, ifudesign, and whether or not the
            datacube has been binned logarithmically in wavelength.
        sres_ext (:obj:`str`, optional):
            The base extension name to use when constructing the
            spectral resolution vectors. See
            :func:`spectral_resolution`. Should be either 'SPECRES'
            or 'DISP'.
        sres_pre (:obj:`bool`, optional):
            Read the pre-pixelized version of the spectral
            resolution, instead of the post-pixelized version.
        sres_fill (:obj:`bool`, optional):
            Fill masked values by interpolation. Default is to leave
            masked pixels in returned array.
        covar_ext (:obj:`str`, optional):
            Extension to use as the single spatial correlation matrix
            for all wavelength channels, read from the DRP file. For
            generating the covariance matrix directly for an
            arbitrary wavelength channel using the RSS file, see
            :func:`covariance_matrix`.
    """
    def __init__(self, ifile, sres_ext='DISP', sres_pre=True, sres_fill=True, covar_ext=None):

        if not os.path.isfile(ifile):
            raise FileNotFoundError('File does not exist: {0}'.format(ifile))

        # Parse the relevant information from the filename
        self.directory_path, _ifile = os.path.split(os.path.abspath(ifile)) 
        self.plate, self.ifudesign, log = _ifile.split('-')[1:4]
        self.plate = int(self.plate)
        self.ifudesign = int(self.ifudesign)
        log = 'LOG' in log

        # Try to define the BitMask object
        try:
            bitmask = DRPFitsBitMask(mode='CUBE')
        except:
            warnings.warn('Unable to define bit mask for MaNGA datacube.  Can only distinguish '
                          ' between masked (values greater than 0) and unmasked (values of 0).')
            bitmask = None

        # Open the file and initialize the base class
        with fits.open(ifile) as hdu:
            print('Reading file data ...', end='\r')
            self.prihdr = hdu[0].header
            self.fluxhdr = hdu['FLUX'].header
            sres = MaNGADataCube.spectral_resolution(hdu, ext=sres_ext, pre=sres_pre,
                                                     fill=sres_fill).filled(0.0)
            covar = None if covar_ext is None \
                        else Covariance.from_fits(hdu, ivar_ext=None, covar_ext=covar_ext,
                                                  impose_triu=True, correlation=True)

            # NOTE: Transposes are done here because of how the data is
            # read from the fits file. The covariance is NOT transposed
            # accordingly because of how the correlation matrices are
            # stored by the DRP. Some care needs to be taken here
            # though because the transpose changes the array memory
            # storage from C to Fortran contiguous. However, everything
            # is expected to be okay because numpy array flattening
            # always performs a C-like flattening, even if the memory
            # storage Fortran contiguous.
            super(MaNGADataCube, self).__init__(hdu['FLUX'].data.T, wave=hdu['WAVE'].data,
                                                ivar=hdu['IVAR'].data.T, mask=hdu['MASK'].data.T,
                                                bitmask=bitmask, sres=sres.T, covar=covar,
                                                wcs=WCS(header=self.fluxhdr, fix=True),
                                                pixelscale=0.5, log=log)
        print('Reading file data ... DONE')

        # Try to use the header to set the DRP version
        self.drpver = self.prihdr['VERSDRP3'] if 'VERSDRP3' in self.prihdr \
                        else defaults.drp_version()
        # Reduction path is always set to the default. A warning is
        # thrown if the default reduction path is not the same as the
        # path expected for the file.
        self.redux_path = defaults.drp_redux_path(drpver=self.drpver)
        if self.directory_path != defaults.drp_directory_path(self.plate, drpver=self.drpver,
                                                              redux_path=self.redux_path):
            warnings.warn('Default reduction path does not match file path.  May not be able to '
                          'find paired RSS file if requested.')

    @staticmethod
    def spectral_resolution_extension(hdu, ext=None, pre=False):
        """
        Determine the spectral resolution channel to use.

        Precedence is to use the DISP extension if it exists, and the
        SPECRES extension otherwise.

        Args:
            hdu (`astropy.io.fits.HDUList`):
                The opened MaNGA DRP file.
            ext (:obj:`str`, optional):
                Specify the extension with the spectral estimate to use.
                Should be in [ None, 'DISP', 'SPECRES'].  The default is
                None, which means it will return, in order of
                precedence, the data in 'DISP', 'SPECRES', or a None
                value if neither are present.
            pre (:obj:`bool`, optional):
                Read the pre-pixelized version of the spectral
                resolution, instead of the post-pixelized version.  This
                prepends 'PRE' to the extension name.

        Returns:
            :obj:`str`: The name of the preferred extension to use.
        """
        sres_ext = [h.name for h in hdu if h.name in ['PREDISP', 'DISP', 'PRESPECRES', 'SPECRES']]
        _ext = ext
        if ext is None:
            _ext = 'PREDISP' if pre else 'DISP'
            if _ext not in sres_ext:
                _ext = 'PRESPECRES' if pre else 'SPECRES'
        return None if _ext not in sres_ext else _ext

    @staticmethod
    def spectral_resolution(hdu, ext=None, pre=False, fill=False, median=False):
        """
        Return the spectral resolution at each spatial and spectral
        position.

        See :func:`spectral_resolution_extension` for a description
        of the precedence used when ``ext`` is None.

        Args:
            hdu (`astropy.io.fits.HDUList`):
                The opened MaNGA DRP file.
            ext (:obj:`str`, optional):
                Specify the extension with the spectral estimate to use.
                Should be in [ None, 'DISP', 'SPECRES'].  The default is
                None, which means it will return, in order of
                precedence, the data in 'DISP', 'SPECRES', or a None
                value if neither are present.
            pre (:obj:`bool`, optional):
                Read the pre-pixelized version of the spectral
                resolution, instead of the post-pixelized version.  This
                prepends 'PRE' to the extension name.
            fill (:obj:`bool`, optional):
                Fill masked values by interpolation.  Default is to
                leave masked pixels in returned array.
            median (:obj:`bool`, optional):
                Return a single vector with the median spectral
                resolution instead of a per spectrum array.  When using
                the `SPECRES` extension, this just returns the vector
                provided by the DRP file; when using the `DISP`
                extension, this performs a masked median across the
                array and then interpolates any wavelengths that were
                masked in all vectors.

        Returns:
            `numpy.ma.MaskedArray`_ : Even if interpolated such that
            there should be no masked values, the function returns a
            masked array.  Array contains the spectral resolution
            (:math:`R = \lambda/\Delta\lambda`) pulled from the DRP
            file.
        """
        # Determine which spectral resolution element to use
        _ext = MaNGADataCube.spectral_resolution_extension(hdu, ext=ext, pre=pre)
        # If no valid extension, raise an exception
        if ext is None and _ext is None:
            raise ValueError('No valid spectral resolution extension.')
        if ext is not None and _ext is None:
            raise ValueError('No extension: {0}'.format(ext))
            
        # Build the spectral resolution vectors
        if 'SPECRES' in _ext:
            sres = numpy.ma.MaskedArray(hdu[_ext].data, mask=numpy.invert(hdu[_ext].data > 0))
            if fill:
                sres = numpy.ma.MaskedArray(interpolate_masked_vector(sres))
            if not median:
                spatial_shape = hdu['FLUX'].data.shape[1:][::-1]
                sres = numpy.ma.tile(sres, spatial_shape + (1,)).T
            return sres

        # Otherwise dealing with the DISP cube
        sres = numpy.ma.MaskedArray(hdu[_ext].data)
        # Mask any non-positive value
        sres[numpy.invert(sres > 0)] = numpy.ma.masked
        # Convert from sigma in angstroms to spectral resolution
        # (based on FWHM)
        sres = numpy.ma.divide(hdu['WAVE'].data[:,None,None], sres) / DAPConstants.sig2fwhm
        # Interpolate over any masked values
        if fill:
            outshape = sres.shape
            sres = numpy.ma.MaskedArray(
                        numpy.ma.apply_along_axis(interpolate_masked_vector, 0,
                                                  sres.reshape(outshape[0], -1))).reshape(outshape)
        if median:
            sres = numpy.ma.median(sres.reshape(outshape[0],-1), axis=1)
        return sres

    @classmethod
    def from_plateifu(cls, plate, ifudesign, log=True, drpver=None, redux_path=None,
                      directory_path=None, **kwargs):
        """
        Construct MaNGA datacube object based on its plate-ifu
        designation.

        The provided plate, ifudesign, and boolean for the log
        binning are used to construct the name of the file. The DRP
        version and reduction and directory paths are used to
        construct the directory with the file.

        Args:
            plate (:obj:`int`):
                Plate number
            ifudesign (:obj:`int`):
                IFU design
            log (:obj:`bool`, optional):
                Use the datacube that is logarithmically binned in
                wavelength.
            drpver (:obj:`str`, optional):
                DRP version, which is used to define the default DRP
                redux path. Default is defined by
                :func:`mangadap.config.defaults.drp_version`
            redux_path (:obj:`str`, optional):
                The path to the top level directory containing the
                DRP output files for a given DRP version. Default is
                defined by
                :func:`mangadap.config.defaults.drp_redux_path`.
            directory_path (:obj:`str`, optional):
                The exact path to the DRP file. Default is defined by
                :func:`mangadap.config.defaults.drp_directory_path`.
                Providing this ignores anything provided for
                ``drpver`` or ``redux_path``.
            **kwargs:
                Keyword arguments passed directly to the primary
                instantiation method; see :class:`MaNGADataCube`.
        """
        # Set the attributes, forcing a known type
        _plate = int(plate)
        _ifudesign = int(ifudesign)

        # Setup the directory path.
        if directory_path is None:
            directory_path = defaults.drp_directory_path(_plate, drpver=drpver,
                                                         redux_path=redux_path)
        return cls(os.path.join(directory_path,
                                MaNGADataCube.file_name(_plate, _ifudesign, log=log)), **kwargs)

    # TODO: Include a class method that instantiates from (or wraps a Marvin Cube)

    @staticmethod
    def file_name(plate, ifudesign, log=True):
        """
        Return the name of the DRP datacube file.

        Args:
            plate (:obj:`int`):
                Plate number
            ifudesign (:obj:`int`):
                IFU design
            log (:obj:`bool`, optional):
                Use the datacube that is logarithmically binned in
                wavelength.

        Returns:
            :obj:`str`: The relevant file name.
        """
        mode = '{0}CUBE'.format('LOG' if log else 'LIN')
        return '{0}.fits.gz'.format(defaults.manga_fits_root(plate, ifudesign, mode))

    def file_path(self):
        """Return the full path to the DRP datacube file."""
        return os.path.join(self.directory_path,
                            self.file_name(self.plate, self.ifudesign, log=self.log))

    @staticmethod
    def do_not_fit_flags():
        """Return the maskbit names that should not be fit."""
        return ['DONOTUSE', 'FORESTAR']

    @staticmethod
    def do_not_stack_flags():
        """Return the maskbit names that should not be stacked."""
        return ['DONOTUSE', 'FORESTAR']

    @staticmethod
    def default_paths(plate, ifudesign, log=True, drpver=None, redux_path=None,
                      directory_path=None, output_file=None):
        """
        Return the primary directory and file name with the DRP fits
        LOG-binned file.

        Args:
            plate (:obj:`int`):
                Plate number
            ifudesign (:obj:`int`):
                IFU design
            log (:obj:`bool`, optional):
                Use the datacube that is logarithmically binned in
                wavelength.
            drpver (:obj:`str`, optional):
                DRP version, which is used to define the default DRP
                redux path. Default is defined by
                :func:`mangadap.config.defaults.drp_version`
            redux_path (:obj:`str`, optional):
                The path to the top level directory containing the
                DRP output files for a given DRP version. Default is
                defined by
                :func:`mangadap.config.defaults.drp_redux_path`.
            directory_path (:obj:`str`, optional):
                The exact path to the DRP file. Default is defined by
                :func:`mangadap.config.defaults.drp_directory_path`.
                Providing this ignores anything provided for
                ``drpver`` or ``redux_path``.
            output_file (:obj:`str`, optional):
                The name of the file with the DRP data. Default set
                by :func:`MaNGADataCube.file_name`.

        Returns:
            :obj:`str`: Two strings with the path to and name of the DRP
            data file.
        """
        _directory_path = defaults.drp_directory_path(plate, drpver=drpver,
                                                      redux_path=redux_path) \
                                if directory_path is None else directory_path
        _output_file = MaNGADataCube.file_name(plate, ifudesign, log=log) \
                            if output_file is None else output_file
        return _directory_path, _output_file

    def mean_sky_coordinates(self, center_coo=None, offset=None):
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
        _center_coo = (self.prihdr['{0}RA'.format(offset.upper())],
                       self.prihdr['{0}DEC'.format(offset.upper())]) \
                           if center_coo is None else center_coo
        return super(MaNGADataCube, self).mean_sky_coordinates(center_coo=_center_coo)
