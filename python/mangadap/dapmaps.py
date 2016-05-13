# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""

Defines the class that constructs and interfaces with the DAP maps file
output.

There has to be a better way to construct the file(s)!

*License*:
    Copyright (c) 2015, SDSS-IV/MaNGA Pipeline Group
        Licensed under BSD 3-clause license - see LICENSE.rst

*Source location*:
    $MANGADAP_DIR/python/mangadap/dapmaps.py

*Imports and python version compliance*:
    ::

*Class usage examples*:
    Add some usage comments here!

*Revision history*:
    | **09 May 2016**: Original Implementation by K. Westfall (KBW)
"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import sys
if sys.version > '3':
    long = int

import numpy
import logging
import resource
import time
import os

from astropy.wcs import WCS
from astropy.io import fits

from .util.log import init_DAP_logging, module_logging, log_output
from .util.fileio import write_hdu
from .drpfits import DRPFits
from .config.defaults import default_dap_method_path, default_dap_file_name
from .proc.reductionassessments import ReductionAssessment
from .proc.spatiallybinnedspectra import SpatiallyBinnedSpectra
from .proc.stellarcontinuummodel import StellarContinuumModel
from .proc.emissionlinemoments import EmissionLineMoments
from .proc.emissionlinemodel import EmissionLineModel
from .proc.spectralindices import SpectralIndices

from matplotlib import pyplot

__author__ = 'Kyle Westfall'


class construct_maps_file:
    """
    Construct a DAP maps file based on the input data.

    Set as a class for coherency reasons, but should not be used as an
    object!

    """
    def __init__(self, drpf, rdxqa=None, binned_spectra=None, stellar_continuum=None,
                 emission_line_moments=None, emission_line_model=None, spectral_indices=None,
                 dapsrc=None, dapver=None, analysis_path=None, directory_path=None,
                 output_file=None, clobber=True):

        # The output method directory is, for now, the combination of
        # the binned_spectrum and stellar_continuum method keys
        if directory_path is None and (binned_spectra is None or stellar_continuum is None):
            raise ValueError('Could not define output directory path.')

        self.method = '{0}-{1}'.format(binned_spectra.method['key'],
                                       stellar_continuum.method['key']) \
                                if directory_path is None else None

        # Type checking
        if not isinstance(drpf, DRPFits):
            raise TypeError('Input must have type DRPFits.')
        if rdxqa is not None and not isinstance(rdxqa, ReductionAssessment):
            raise TypeError('Input must have type ReductionAssessment.')
        if binned_spectra is not None and not isinstance(binned_spectra, SpatiallyBinnedSpectra):
            raise TypeError('Input must have type SpatiallyBinnedSpectra.')
        if stellar_continuum is not None and not isinstance(stellar_continuum,
                                                            StellarContinuumModel):
            raise TypeError('Input must have type StellarContinuumModel.')
        if emission_line_moments is not None and not isinstance(emission_line_moments,
                                                                EmissionLineMoments):
            raise TypeError('Input must have type EmissionLineMoments.')
        if emission_line_model is not None and not isinstance(emission_line_model,
                                                              EmissionLineModel):
            raise TypeError('Input must have type EmissionLineModel.')
        if spectral_indices is not None and not isinstance(spectral_indices, SpectralIndices):
            raise TypeError('Input must have type SpectralIndices.')

        self.drpf = drpf

        # Set the output directory path
        # TODO: Get DAP version from __version__ string
        self.directory_path = default_dap_method_path(self.method, self.drpf.plate,
                                                      drpver=self.drpf.drpver, dapver=dapver,
                                                      analysis_path=analysis_path) \
                                                if directory_path is None else str(directory_path)

        # Set the output file
        self.output_file = default_dap_file_name(self.drpf.plate, self.drpf.ifudesign,
                                                 self.drpf.mode, self.method) \
                                    if output_file is None else str(output_file)

        ofile = os.path.join(self.directory_path, self.output_file)
        if os.path.isfile(ofile) and not clobber:
            # TODO: Perform some checks to make sure the existing file
            # has the correct content?
            warnings.warn('File exists: {0}!  Set clobber=True to overwrite.'.format(ofile))
            return

        # Report
        print('Output path: ', self.directory_path)
        print('Output file: ', self.output_file)

        # Initialize the primary header
        prihdr = self._primary_header()

        # Get the base map header
        self.multichannel_maphdr = self._map_header(nchannels=2)
        self.singlechannel_maphdr = self._map_header(nchannels=1)

        # Construct the hdu list for each input object.
        # Reduction assessments:
        self.rdxqalist = self.reduction_assessment_maps(prihdr, rdxqa)
        # Binned spectra:
        self.bspeclist = self.binned_spectra_maps(prihdr, binned_spectra)
        # Stellar-continuum fits:
        self.contlist = self.stellar_continuum_maps(prihdr, stellar_continuum)
        # Emission-line moments:
        self.emlmomlist = self.emission_line_moment_maps(prihdr, emission_line_moments)
        # Emission-line models:
        self.elmodlist = self.emission_line_model_maps(prihdr, emission_line_model)
        # Spectral indices:
        self.sindxlist = self.spectral_index_maps(prihdr, spectral_indices)

        # Save the data to the hdu attribute
        self.hdu = fits.HDUList([ fits.PrimaryHDU(header=prihdr),
                                  *self.rdxqalist,
                                  *self.bspeclist,
                                  *self.contlist,
                                  *self.emlmomlist,
                                  *self.elmodlist,
                                  *self.sindxlist ])

        # Write the file
        if not os.path.isdir(self.directory_path):
            os.makedirs(self.directory_path)
        write_hdu(self.hdu, ofile, clobber=clobber, checksum=True)



    @staticmethod
    def _clean_primary_header(hdr):
        # Remove some keys that are incorrect for DAP data
        hdr.remove('BSCALE')
        hdr.remove('BZERO')
        hdr.remove('BUNIT')
        hdr.remove('MASKNAME')


    def _primary_header(self):
        hdr = self.drpf.hdu['PRIMARY'].header.copy()
        self._clean_primary_header(hdr)
        # Add Authors
        hdr['AUTHOR'] = 'K Westfall & B Andrews <kyle.westfall@port.co.uk, andrewsb@pitt.edu>'
        # Add versioning
        # Add DAP quality
        # Other meta data?
        return hdr


    def _clean_map_header(self, hdr, nchannels=1):

        # Change header keywords to the default values for the third axis
        if nchannels > 1:
            hdr['NAXIS'] = 3
            hdr.remove('CTYPE3')
            hdr.remove('CUNIT3')
            hdr['CTYPE3'] = ' '
            hdr['CUNIT3'] = ' '
            hdr['CRPIX3'] = 1
            hdr['CRVAL3'] = 1.
            hdr['CD3_3']  = 1.
        else:
            hdr['NAXIS'] = 2
            #hdr.remove('NAXIS3')
            hdr.remove('CTYPE3')
            hdr.remove('CUNIT3')
            hdr.remove('CRPIX3')
            hdr.remove('CRVAL3')
            hdr.remove('CD3_3')

        # Remove everything but the WCS information
        w = WCS(header=hdr)
        hdr = w.to_header().copy()

        # Fix the DATE-OBS keyword:
        hdr.comments['DATE-OBS'] = 'Date of median exposure'
        hdr.comments['MJD-OBS'] = '[d] MJD for DATE-OBS'

        # Add back in the BSCALE and BZERO values
        hdr['BSCALE'] = 1.
        hdr['BZERO'] = 0.

        # Add back in the default CTYPE3, CUNIT3
        if nchannels > 1:
            hdr['CTYPE3'] = (' ', 'Undefined type')
            hdr['CUNIT3'] = (' ', 'Undefined units')

        return hdr


    def _map_header(self, nchannels=1):
        hdr = self.drpf.hdu['FLUX'].header.copy()
        hdr = self._clean_map_header(hdr, nchannels=nchannels)
        # Add Authors
        hdr['AUTHOR'] = 'K Westfall & B Andrews <kyle.westfall@port.co.uk, andrewsb@pitt.edu>'
        # Add versioning
        # Add DAP quality
        # Other meta data?
        return hdr


    @staticmethod
    def _mask_data_type(bit_type):
        if bit_type == numpy.uint64:
            return ('FLAG64BIT', '64-bit mask')
        if bit_type == numpy.uint32:
            return ('FLAG32BIT', '32-bit mask')
        if bit_type == numpy.uint16:
            return ('FLAG16BIT', '16-bit mask')
        if bit_type == numpy.uint8:
            return ('FLAG8BIT', '8-bit mask')
        if bit_type == numpy.bool:
            return ('MASKZERO', 'Binary mask; zero values are good/unmasked')
        raise ValueError('Invalid bit_type: {0}!'.format(str(bit_type)))

    
    @staticmethod
    def _add_channel_names(hdr, names, units=None):
        _hdr = hdr.copy()
        nchannels = len(names)
        ndig = int(numpy.log10(nchannels))+1
        if units is not None and len(units) != nchannels:
            raise ValueError('Length of column names and units are not the same.')

        for i in range(nchannels):
            _hdr['C'+'{0}'.format(i+1).zfill(ndig)] = (names[i], 'Data in channel {0}'.format(i+1))
            if units is not None:
                _hdr['U'+'{0}'.format(i+1).zfill(ndig)] \
                            = (units[i], 'Units of data in channel {0}'.format(i+1))
        return _hdr


    @staticmethod
    def _finalize_map_header(hdr, ext, bunit=None, hduclas2='DATA', err=False, qual=False,
                             nchannels=None, bit_type=None):

        _hdr = hdr.copy()

        if bunit is not None:
            _hdr['BUNIT'] = (bunit, 'Unit of pixel value')

        # Add the common HDUCLASS keys
        _hdr['HDUCLASS'] = ('SDSS', 'SDSS format class')
        _hdr['HDUCLAS1'] = ('IMAGE' if nchannels is None or nchannels == 1 else 'CUBE', \
                                        'Data format')
        if hduclas2 == 'DATA':
            _hdr['HDUCLAS2'] = 'DATA'
            if err:
                _hdr['ERRDATA'] = (ext+'_IVAR', 'Associated inv. variance extension')
            if qual:
                _hdr['QUALDATA'] = (ext+'_MASK', 'Associated quality extension')
            return _hdr

        if hduclas2 == 'ERROR':
            _hdr['HDUCLAS2'] = 'ERROR'
            _hdr['HDUCLAS3'] = ('INVMSE', 'Value is inverse mean-square error')
            _hdr['SCIDATA'] = (ext, 'Associated data extension')
            if qual:
                _hdr['QUALDATA'] = (ext+'_MASK', 'Associated quality extension')
            return _hdr

        if hduclas2 == 'QUALITY':
            _hdr['HDUCLAS2'] = 'QUALITY'
            _hdr['HDUCLAS3'] = construct_maps_file._mask_data_type(bit_type)
            _hdr['SCIDATA'] = (ext, 'Associated data extension')
            if err:
                _hdr['ERRDATA'] = (ext+'_IVAR', 'Associated inv. variance extension')
            return _hdr
            
        raise ValueError('HDUCLAS2 must be DATA, ERROR, or QUALITY.')


#    def _single_reconstructed_hdu(self, binned_spectra, inext, incol, unique_bins, reconstruct,
#                                  indx, outext, bunit=None, hduclas2='DATA', err=False,
#                                  qual=False, bit_type=None):
#        data = numpy.zeros(binned_spectra.spatial_shape, dtype=numpy.float)
#        data.ravel()[indx] = binned_spectra[inext].data[incol][unique_bins[reconstruct[indx]]]
#        return fits.ImageHDU(data=data.T, name=outext,
#                             header=self._finalize_map_header(self.singlechannel_maphdr, outext,
#                                                              bunit=bunit, hduclas2=hduclas2,
#                                                              err=err, qual=qual,
#                                                              bit_type=bit_type))


    def reduction_assessment_maps(self, prihdr, rdxqa):
        """
        Transposes are necessary to keep x,y orientation consistent with
        DRP output.
        """

        # Add data to the primary header
        if rdxqa is not None:
            prihdr['RDXQAKEY'] = (rdxqa.method['key'])
            prihdr['PA'] = (rdxqa.pa, 'Isophotal position angle')
            prihdr['ELL'] = (rdxqa.ell, 'Isophotal ellipticity (1-b/a)')

        # On-sky coordinates
        hdr = self._add_channel_names(self.multichannel_maphdr, ['On-sky X', 'On-sky Y'])
        data = None if rdxqa is None else \
                    rdxqa['SPECTRUM'].data['SKY_COO'].reshape(*self.drpf.spatial_shape, -1).T
        hdus = [ fits.ImageHDU(data=data,
                               header=self._finalize_map_header(hdr, 'SPXL_SKYCOO', bunit='arcsec',
                                                                nchannels=2),
                               name='SPXL_SKYCOO') ]

        # Elliptical coordinates
        hdr = self._add_channel_names(self.multichannel_maphdr,
                                      ['Elliptical radius', 'Elliptical azimuth'],
                                      ['arcsec', 'degrees'])
        if rdxqa is not None:
            hdr['PA'] = (rdxqa.pa, 'Isophotal position angle')
            hdr['ELL'] = (rdxqa.ell, 'Isophotal ellipticity (1-b/a)')
        data = None if rdxqa is None else \
                    rdxqa['SPECTRUM'].data['ELL_COO'].reshape(*self.drpf.spatial_shape, -1).T
        hdus += [ fits.ImageHDU(data=data,
                                header=self._finalize_map_header(hdr, 'SPXL_ELLCOO', bunit='arcsec',
                                                                 nchannels=2),
                                name='SPXL_ELLCOO') ]

        # Spectral coverage and SNR assessments of the DRP data
        data = None if rdxqa is None else \
                    rdxqa['SPECTRUM'].data['FGOODPIX'].reshape(self.drpf.spatial_shape).T
        hdus += [ fits.ImageHDU(data=data,
                                header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                 'SPXL_SPCOV'),
                                name='SPXL_SPCOV') ]
        data = None if rdxqa is None else \
                    rdxqa['SPECTRUM'].data['SIGNAL'].reshape(self.drpf.spatial_shape).T
        hdus += [ fits.ImageHDU(data=data,
                                header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                 'SPXL_MFLUX',
                                                                bunit='1E-17 erg/s/cm^2/ang/spaxel',
                                                                 err=True),
                                name='SPXL_MFLUX') ]

        if rdxqa is None:
            ivar = None
        else:
            ivar = rdxqa['SPECTRUM'].data['VARIANCE'].reshape(self.drpf.spatial_shape).T
            indx = ivar > 0
            ivar[indx] = 1.0/ivar[indx]
            ivar[numpy.invert(indx)] = 0.0
        hdus += [ fits.ImageHDU(data=data,
                                header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                 'SPXL_MFLUX', hduclas2='ERROR',
                                                        bunit='(1E-17 erg/s/cm^2/ang/spaxel)^{-2}'),
                                name='SPXL_MFLUX_IVAR') ]
        data = None if rdxqa is None else \
                    rdxqa['SPECTRUM'].data['SNR'].reshape(self.drpf.spatial_shape).T
        hdus += [ fits.ImageHDU(data=data,
                                header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                 'SPXL_SNR'),
                                name='SPXL_SNR') ]

        return hdus

#        pyplot.imshow(radius.T, origin='lower', interpolation='nearest')
#        pyplot.show()
#        pyplot.imshow(self.drpf['FLUX'].data[:,:,1000].T, origin='lower', interpolation='nearest')
#        pyplot.show()
#        exit()


    def binned_spectra_maps(self, prihdr, binned_spectra):
        """
        Transposes are necessary to keep x,y orientation consistent with
        DRP output.
        """

        # Construct and return the empty hdus
        if binned_spectra is None:
            hdus = [ fits.ImageHDU(data=None,
                                    header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                     'BINID'),
                                    name='BINID') ]

            hdus += [ fits.ImageHDU(data=None,
                                    header=self._finalize_map_header(hdr, 'BIN_LWSKYCOO',
                                                                     bunit='arcsec', nchannels=2),
                                    name='BIN_LWSKYCOO') ]

            hdus += [ fits.ImageHDU(data=None,
                                    header=self._finalize_map_header(hdr, 'BIN_LWELLCOO',
                                                                     bunit='arcsec', nchannels=2),
                                    name='BIN_LWELLCOO') ]

            hdus += [ fits.ImageHDU(data=None,
                                    header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                     'BIN_AREA', bunit='arcsec^2'),
                                    name='BIN_AREA') ]

            hdus += [ fits.ImageHDU(data=None,
                                    header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                     'BIN_FAREA'),
                                    name='BIN_FAREA') ]

            hdus += [ fits.ImageHDU(data=None,
                                    header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                     'BIN_MFLUX',
                                                                bunit='1E-17 erg/s/cm^2/ang/spaxel',
                                                                     err=True),
                                    name='BIN_MFLUX') ]
           
            hdus += [ fits.ImageHDU(data=None,
                                    header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                     'BIN_MFLUX', hduclas2='ERROR',
                                                        bunit='(1E-17 erg/s/cm^2/ang/spaxel)^{-2}'),
                                    name='BIN_MFLUX_IVAR') ]

            hdus += [ fits.ImageHDU(data=None,
                                    header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                     'BIN_SNR'),
                                    name='BIN_SNR') ]

            return hdus

        # Add data to the primary header
        if binned_spectra.reff is not None:
            prihdr['REFF'] = (binned_spectra.reff, 'Effective radius (arcsec)')
        prihdr['BINKEY'] = (binned_spectra.method['key'], 'Spectal binning method keyword')
        prihdr['BINMINSN'] = (binned_spectra.method['minimum_snr'],
                                    'Minimum S/N of spectrum to include')
        prihdr['FSPCOV'] = (binned_spectra['PRIMARY'].header['FSPCOV'],
                                    'Minimum allowed valid spectral coverage fraction')
        prihdr['NBINS'] = (binned_spectra.nbins, 'Number of unique spatial bins')
        if len(binned_spectra.missing_bins) > 0:
            prihdr['EMPTYBIN'] = (str(binned_spectra.missing_bins), 'List of bins with no data')

        # Bin index
        bin_indx = binned_spectra['BINID'].data.copy().ravel()
        hdus = [ fits.ImageHDU(data=bin_indx.reshape(binned_spectra.spatial_shape).T,
                               header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                'BINID'),
                               name='BINID') ]
#        pyplot.imshow(bin_indx.reshape(binned_spectra.spatial_shape).T, origin='lower',
#                      interpolation='nearest')
#        pyplot.show()

        # Determine how to reconstruct the image data
        unique_bins, reconstruct = numpy.unique(bin_indx, return_inverse=True)
        indx = bin_indx > -1

        # Get the luminosity weighted Cartesian coordinates
        data = numpy.zeros((*binned_spectra.spatial_shape, 2), dtype=numpy.float)
        data.reshape(-1,2)[indx,0] = binned_spectra['BINS'].data['LW_SKY_COO'][
                                                            unique_bins[reconstruct[indx]],0]
        data.reshape(-1,2)[indx,1] = binned_spectra['BINS'].data['LW_SKY_COO'][
                                                            unique_bins[reconstruct[indx]],1]
        hdr = self._add_channel_names(self.multichannel_maphdr,
                                      ['Lum. weighted on-sky X', 'Lum. weighted on-sky Y'])
        hdus += [ fits.ImageHDU(data=data.T,
                                header=self._finalize_map_header(hdr, 'BIN_LWSKYCOO',
                                                                 bunit='arcsec', nchannels=2),
                                name='BIN_LWSKYCOO') ]

#        pyplot.imshow(data.T[1,:,:], origin='lower', interpolation='nearest')
#        pyplot.show()

        # Get the luminosity-weighted polar coordinates
        data = numpy.zeros((*binned_spectra.spatial_shape, 2), dtype=numpy.float)
        data.reshape(-1,2)[indx,0] = binned_spectra['BINS'].data['LW_ELL_COO'][
                                                            unique_bins[reconstruct[indx]],0]
        data.reshape(-1,2)[indx,1] = binned_spectra['BINS'].data['LW_ELL_COO'][
                                                            unique_bins[reconstruct[indx]],1]
        hdr = self._add_channel_names(self.multichannel_maphdr,
                                      ['Lum. weighted elliptical radius',
                                       'Lum. weighted elliptical azimuth'],
                                      ['arcsec', 'degrees'])
        if binned_spectra.rdxqa is not None:
            hdr['PA'] = (binned_spectra.rdxqa.pa, 'Isophotal position angle')
            hdr['ELL'] = (binned_spectra.rdxqa.ell, 'Isophotal ellipticity (1-b/a)')
        hdus += [ fits.ImageHDU(data=data.T,
                                header=self._finalize_map_header(hdr, 'BIN_LWELLCOO',
                                                                 bunit='arcsec', nchannels=2),
                                name='BIN_LWELLCOO') ]

#        pyplot.imshow(lwcoo.T[0,:,:], origin='lower', interpolation='nearest')
#        pyplot.show()

        # Bin area
        data = numpy.zeros(binned_spectra.spatial_shape, dtype=numpy.float)
        data.ravel()[indx] = binned_spectra['BINS'].data['AREA'][unique_bins[reconstruct[indx]]]
        hdus += [ fits.ImageHDU(data=data.T,
                                header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                 'BIN_AREA', bunit='arcsec^2'),
                                name='BIN_AREA') ]
        # Fractional bin area
        data = numpy.zeros(binned_spectra.spatial_shape, dtype=numpy.float)
        data.ravel()[indx] = binned_spectra['BINS'].data['AREA_FRAC'][
                                                                unique_bins[reconstruct[indx]]]
        hdus += [ fits.ImageHDU(data=data.T,
                                header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                 'BIN_FAREA'),
                                name='BIN_FAREA') ]


        # Mean flux in the bin
        data = numpy.zeros(binned_spectra.spatial_shape, dtype=numpy.float)
        data.ravel()[indx] = binned_spectra['BINS'].data['SIGNAL'][unique_bins[reconstruct[indx]]]
        hdus += [ fits.ImageHDU(data=data.T,
                                header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                 'BIN_MFLUX',
                                                             bunit='1E-17 erg/s/cm^2/ang/spaxel',
                                                                 err=True),
                                name='BIN_MFLUX') ]
        # Inverse variance in the mean flux
        data = numpy.zeros(binned_spectra.spatial_shape, dtype=numpy.float)
        data.ravel()[indx] = binned_spectra['BINS'].data['VARIANCE'][unique_bins[reconstruct[indx]]]
        pos = data > 0
        data[pos] = 1.0/data[pos]
        data[numpy.invert(pos)] = 0.0
        hdus += [ fits.ImageHDU(data=data.T,
                                header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                 'BIN_MFLUX', hduclas2='ERROR',
                                                        bunit='(1E-17 erg/s/cm^2/ang/spaxel)^{-2}'),
                                name='BIN_MFLUX_IVAR') ]
        # Signal-to-noise in the bin
        data = numpy.zeros(binned_spectra.spatial_shape, dtype=numpy.float)
        data.ravel()[indx] = binned_spectra['BINS'].data['SNR'][unique_bins[reconstruct[indx]]]
        hdus += [ fits.ImageHDU(data=data.T,
                                header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                 'BIN_SNR'),
                                name='BIN_SNR') ]

        return hdus


    def stellar_continuum_maps(self, prihdr, stellar_continuum):
        """
        Transposes are necessary to keep x,y orientation consistent with
        DRP output.
        """

        # Construct and return the empty hdus
        if stellar_continuum is None:
            hdus = [ fits.ImageHDU(data=None,
                                    header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                     'STELLAR_VEL', err=True,
                                                                     qual=True, bunit='km/s'),
                                    name='STELLAR_VEL') ]
            hdus += [ fits.ImageHDU(data=None,
                                    header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                     'STELLAR_VEL',
                                                                     bunit='(km/s)^{-2}',
                                                                     hduclas2='ERROR', qual=True),
                                    name='STELLAR_VEL_IVAR') ]
            hdus += [ fits.ImageHDU(data=None,
                                    header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                     'STELLAR_VEL',
                                                                     hduclas2='QUALITY', err=True,
                                                                     bit_type=numpy.bool),
                                    name='STELLAR_VEL_MASK') ]
            hdus = [ fits.ImageHDU(data=None,
                                    header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                     'STELLAR_SIGMA', err=True,
                                                                     qual=True, bunit='km/s'),
                                    name='STELLAR_SIGMA') ]
            hdus += [ fits.ImageHDU(data=None,
                                    header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                     'STELLAR_SIGMA',
                                                                     bunit='(km/s)^{-2}',
                                                                     hduclas2='ERROR', qual=True),
                                    name='STELLAR_SIGMA_IVAR') ]
            hdus += [ fits.ImageHDU(data=None,
                                    header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                     'STELLAR_SIGMA',
                                                                     hduclas2='QUALITY', err=True,
                                                                     bit_type=numpy.bool),
                                    name='STELLAR_SIGMA_MASK') ]
            return hdus

        # Add data to the primary header
        prihdr['VSTEP'] = (stellar_continuum['PRIMARY'].header['VSTEP'],
                                'Velocity step per spectral channel.')
        prihdr['SCKEY'] = (stellar_continuum.method['key'],
                                'Stellar-continuum modeling method keyword')
        prihdr['SCMINSN'] = (stellar_continuum.method['minimum_snr'],
                                'Minimum S/N of spectrum to include')
        # Add the guess kinematics, if provided
        try:
            prihdr['INPVEL'] = (stellar_continuum['PRIMARY'].header['INPVEL'],
                                    'Initial guess velocity')
        except:
            pass

        try:
            prihdr['INPSIG'] = (stellar_continuum['PRIMARY'].header['INPSIG'],
                                    'Initial guess velocity dispersion')
        except:
            pass
#        prihdr['NMOD'] = (stellar_continuum.nmodels, 'Number of unique stellar-continuum models')
#        if len(self.missing_models) > 0:
#            prihdr['EMPTYMOD'] = (str(self.missing_models), 'List of models with no data')

        # Bin index
        bin_indx = stellar_continuum['BINID'].data.copy().ravel()

        # Determine how to reconstruct the image data
        unique_bins, reconstruct = numpy.unique(bin_indx, return_inverse=True)
        indx = bin_indx > -1

        # Stellar velocity
        data = numpy.zeros(stellar_continuum.spatial_shape, dtype=numpy.float)
        data.ravel()[indx] = stellar_continuum['PAR'].data['KIN'][unique_bins[reconstruct[indx]],0]
        hdus = [ fits.ImageHDU(data=data.T,
                               header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                'STELLAR_VEL', bunit='km/s',
                                                                err=True, qual=True),
                               name='STELLAR_VEL') ]
        # Inverse variance
        data = numpy.zeros(stellar_continuum.spatial_shape, dtype=numpy.float)
        data.ravel()[indx] = stellar_continuum['PAR'].data['KINERR'][
                                                                unique_bins[reconstruct[indx]],0]
        pos = data > 0
        data[pos] = numpy.square(1.0/data[pos])
        data[numpy.invert(pos)] = 0.0
        hdus += [ fits.ImageHDU(data=data.T,
                                header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                 'STELLAR_VEL', hduclas2='ERROR',
                                                                 bunit='(km/s)^{-2}', qual=True),
                                name='STELLAR_VEL_IVAR') ]
        # Bitmask
        hdus += [ fits.ImageHDU(data=numpy.zeros(data.shape, dtype=numpy.uint8),
                                header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                 'STELLAR_VEL',
                                                                 hduclas2='QUALITY', err=True,
                                                                 bit_type=numpy.bool),
                                name='STELLAR_VEL_MASK') ]

        # Stellar velocity dispersion
        data = numpy.zeros(stellar_continuum.spatial_shape, dtype=numpy.float)
        data.ravel()[indx] = stellar_continuum['PAR'].data['KIN'][unique_bins[reconstruct[indx]],1]
        hdus += [ fits.ImageHDU(data=data.T,
                                header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                 'STELLAR_SIGMA', bunit='km/s',
                                                                 err=True, qual=True),
                                name='STELLAR_SIGMA') ]
        # Inverse variance
        data = numpy.zeros(stellar_continuum.spatial_shape, dtype=numpy.float)
        data.ravel()[indx] = stellar_continuum['PAR'].data['KINERR'][
                                                                unique_bins[reconstruct[indx]],1]
        pos = data > 0
        data[pos] = numpy.square(1.0/data[pos])
        data[numpy.invert(pos)] = 0.0
        hdus += [ fits.ImageHDU(data=data.T,
                                header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                 'STELLAR_SIGMA', hduclas2='ERROR',
                                                                 bunit='(km/s)^{-2}', qual=True),
                                name='STELLAR_SIGMA_IVAR') ]
        # Bitmask
        hdus += [ fits.ImageHDU(data==numpy.zeros(data.shape, dtype=numpy.uint8),
                                header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                 'STELLAR_SIGMA',
                                                                 hduclas2='QUALITY', err=True,
                                                                 bit_type=numpy.bool),
                                name='STELLAR_SIGMA_MASK') ]

        return hdus


    def emission_line_moment_maps(self, prihdr, emission_line_moments):
        """
        Transposes are necessary to keep x,y orientation consistent with
        DRP output.
        """

        # Construct and return the empty hdus
        if emission_line_moments is None:
            hdus = [ fits.ImageHDU(data=None,
                                    header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                     'EMLINE_SFLUX', err=True,
                                                                     qual=True,
                                                                 bunit='1E-17 erg/s/cm^2/spaxel'),
                                    name='EMLINE_SFLUX') ]
            hdus += [ fits.ImageHDU(data=None,
                                    header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                     'EMLINE_SFLUX',
                                                            bunit='(1E-17 erg/s/cm^2/spaxel)^{-2}',
                                                                     hduclas2='ERROR', qual=True),
                                    name='EMLINE_SFLUX_IVAR') ]
            hdus += [ fits.ImageHDU(data=None,
                                    header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                     'EMLINE_SFLUX',
                                                                     hduclas2='QUALITY', err=True,
                                                                     bit_type=numpy.bool),
                                    name='EMLINE_SFLUX_MASK') ]
            return hdus

        # Add data to the primary header
        prihdr['ELMKEY'] = (emission_line_moments.database['key'],
                                'Emission-line moments database keyword')
        prihdr['ELMMINSN'] = (emission_line_moments.database['minimum_snr'],
                                'Minimum S/N of spectrum to include')
        prihdr['ARTDB'] = (emission_line_moments.database['artifacts'], 'Artifact database keyword')
        prihdr['MOMDB'] = (emission_line_moments.database['passbands'],
                                'Emission-line moments database keyword')

        # Bin index
        bin_indx = emission_line_moments.binned_spectra['BINID'].data.copy().ravel()

        # Determine how to reconstruct the image data
        unique_bins, reconstruct = numpy.unique(bin_indx, return_inverse=True)
        indx = bin_indx > -1

        # Integrated flux
        data = numpy.zeros((*emission_line_moments.spatial_shape, emission_line_moments.nmom),
                           dtype=numpy.float)
        data.reshape(-1, emission_line_moments.nmom)[indx,:] \
                = emission_line_moments['ELMMNTS'].data['FLUX'][unique_bins[reconstruct[indx]],:]
        # TODO: Convert wavelengths to air for column names?
        cols = [ '{0}-{1:.0f}'.format(n,w) \
                        for n,w in zip(emission_line_moments['ELMBAND'].data['NAME'],
                                       emission_line_moments['ELMBAND'].data['RESTWAVE']) ]
        hdr = self._add_channel_names(self.multichannel_maphdr, cols)
        hdus = [ fits.ImageHDU(data=data.T,
                               header=self._finalize_map_header(hdr, 'EMLINE_SFLUX', err=True,
                                                                qual=True,
                                                                bunit='1E-17 erg/s/cm^2/spaxel'),
                                    name='EMLINE_SFLUX') ]
        # Inverse variance
        data = numpy.zeros((*emission_line_moments.spatial_shape, emission_line_moments.nmom),
                           dtype=numpy.float)
        data.reshape(-1, emission_line_moments.nmom)[indx,:] \
                = emission_line_moments['ELMMNTS'].data['FLUXERR'][unique_bins[reconstruct[indx]],:]
        pos = data > 0
        data[pos] = numpy.square(1.0/data[pos])
        data[numpy.invert(pos)] = 0.0
        hdus += [ fits.ImageHDU(data=data.T,
                                header=self._finalize_map_header(hdr, 'EMLINE_SFLUX',
                                                                 hduclas2='ERROR', qual=True,
                                                            bunit='(1E-17 erg/s/cm^2/spaxel)^{-2}'),
                                name='EMLINE_SFLUX_IVAR') ]
        # Bitmask
        bit_type = emission_line_moments.bitmask.minimum_uint_dtype()
        data = numpy.zeros((*emission_line_moments.spatial_shape, emission_line_moments.nmom),
                           dtype=bit_type)
        data.reshape(-1, emission_line_moments.nmom)[indx,:] \
                = emission_line_moments['ELMMNTS'].data['MASK'][unique_bins[reconstruct[indx]],:]
        hdus += [ fits.ImageHDU(data=data.T,
                                header=self._finalize_map_header(hdr, 'EMLINE_SFLUX',
                                                                 hduclas2='QUALITY', err=True,
                                                                 bit_type=bit_type),
                                name='EMLINE_SFLUX_MASK') ]

        return hdus


    def emission_line_model_maps(self, prihdr, emission_line_model):
        """
        Transposes are necessary to keep x,y orientation consistent with
        DRP output.
        """

        # Construct and return the empty hdus
        if emission_line_model is None:
            hdus = [ fits.ImageHDU(data=None,
                                   header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                    'EMLINE_GFLUX', err=True,
                                                                    qual=True,
                                                                 bunit='1E-17 erg/s/cm^2/spaxel'),
                                   name='EMLINE_GFLUX') ]
            hdus += [ fits.ImageHDU(data=None,
                                    header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                     'EMLINE_GFLUX',
                                                            bunit='(1E-17 erg/s/cm^2/spaxel)^{-2}',
                                                                     hduclas2='ERROR', qual=True),
                                    name='EMLINE_GFLUX_IVAR') ]
            hdus += [ fits.ImageHDU(data=None,
                                    header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                     'EMLINE_GFLUX',
                                                                     hduclas2='QUALITY', err=True,
                                                                     bit_type=numpy.bool),
                                    name='EMLINE_GFLUX_MASK') ]
            hdus = [ fits.ImageHDU(data=None,
                                   header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                    'EMLINE_GVEL', err=True,
                                                                    qual=True,
                                                                    bunit='km/s'),
                                   name='EMLINE_GVEL') ]
            hdus += [ fits.ImageHDU(data=None,
                                    header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                     'EMLINE_GVEL',
                                                                     bunit='(km/s)^{-2}',
                                                                     hduclas2='ERROR', qual=True),
                                    name='EMLINE_GVEL_IVAR') ]
            hdus += [ fits.ImageHDU(data=None,
                                    header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                     'EMLINE_GVEL',
                                                                     hduclas2='QUALITY', err=True,
                                                                     bit_type=numpy.bool),
                                    name='EMLINE_GVEL_MASK') ]
            hdus = [ fits.ImageHDU(data=None,
                                   header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                    'EMLINE_GSIGMA', err=True,
                                                                    qual=True,
                                                                    bunit='km/s'),
                                   name='EMLINE_GSIGMA') ]
            hdus += [ fits.ImageHDU(data=None,
                                    header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                     'EMLINE_GSIGMA',
                                                                     bunit='(km/s)^{-2}',
                                                                     hduclas2='ERROR', qual=True),
                                    name='EMLINE_GSIGMA_IVAR') ]
            hdus += [ fits.ImageHDU(data=None,
                                    header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                     'EMLINE_GSIGMA',
                                                                     hduclas2='QUALITY', err=True,
                                                                     bit_type=numpy.bool),
                                    name='EMLINE_GSIGMA_MASK') ]
            return hdus

        # Add data to the primary header
        prihdr['ELFKEY'] = (emission_line_model.method['key'],
                                'Emission-line modeling method keyword')
        prihdr['ELFMINSN'] = (emission_line_model.method['minimum_snr'],
                                'Minimum S/N of spectrum to include')
        prihdr['ARTDB'] = (emission_line_model.method['artifacts'],
                                'Artifact database keyword')
        prihdr['EMLDB'] = (emission_line_model.method['emission_lines'],
                                'Emission-line database keyword')

        # Bin index
        bin_indx = emission_line_model.binned_spectra['BINID'].data.copy().ravel()

        # Determine how to reconstruct the image data
        unique_bins, reconstruct = numpy.unique(bin_indx, return_inverse=True)
        indx = bin_indx > -1

        # Channel identifiers (units are all the same)
        # TODO: Convert wavelengths to air for column names?
        cols = [ '{0}-{1:.0f}'.format(n,w) \
                        for n,w in zip(emission_line_model['PAR'].data['NAME'],
                                       emission_line_model['PAR'].data['RESTWAVE']) ]
        hdr = self._add_channel_names(self.multichannel_maphdr, cols)

        # Integrated flux
        data = numpy.zeros((*emission_line_model.spatial_shape, emission_line_model.neml),
                           dtype=numpy.float)
        data.reshape(-1, emission_line_model.neml)[indx,:] \
                = emission_line_model['EMLDATA'].data['FLUX'][unique_bins[reconstruct[indx]],:]
        hdus = [ fits.ImageHDU(data=data.T,
                               header=self._finalize_map_header(hdr, 'EMLINE_GFLUX', err=True,
                                                                qual=True,
                                                                bunit='1E-17 erg/s/cm^2/spaxel'),
                               name='EMLINE_GFLUX') ]
        # Inverse variance
        data = numpy.zeros((*emission_line_model.spatial_shape, emission_line_model.neml),
                           dtype=numpy.float)
        data.reshape(-1, emission_line_model.neml)[indx,:] \
                = emission_line_model['EMLDATA'].data['FLUXERR'][unique_bins[reconstruct[indx]],:]
        pos = data > 0
        data[pos] = numpy.square(1.0/data[pos])
        data[numpy.invert(pos)] = 0.0
        hdus += [ fits.ImageHDU(data=data.T,
                                header=self._finalize_map_header(hdr, 'EMLINE_GFLUX',
                                                                 hduclas2='ERROR', qual=True,
                                                            bunit='(1E-17 erg/s/cm^2/spaxel)^{-2}'),
                                name='EMLINE_GFLUX_IVAR') ]
        # Bitmask
        bit_type = emission_line_model.bitmask.minimum_uint_dtype()
        data = numpy.zeros((*emission_line_model.spatial_shape, emission_line_model.neml),
                           dtype=bit_type)
        data.reshape(-1, emission_line_model.neml)[indx,:] \
                = emission_line_model['EMLDATA'].data['MASK'][unique_bins[reconstruct[indx]],:]
        hdus += [ fits.ImageHDU(data=data.T,
                                header=self._finalize_map_header(hdr, 'EMLINE_GFLUX',
                                                                 hduclas2='QUALITY', err=True,
                                                                 bit_type=bit_type),
                                name='EMLINE_GFLUX_MASK') ]

        # Velocity
        data = numpy.zeros((*emission_line_model.spatial_shape, emission_line_model.neml),
                           dtype=numpy.float)
        data.reshape(-1, emission_line_model.neml)[indx,:] \
                = emission_line_model['EMLDATA'].data['KIN'][unique_bins[reconstruct[indx]],:,0]
        hdus += [ fits.ImageHDU(data=data.T,
                                header=self._finalize_map_header(hdr, 'EMLINE_GVEL', err=True,
                                                                 qual=True,
                                                                 bunit='km/s'),
                                name='EMLINE_GVEL') ]
        # Inverse variance
        data = numpy.zeros((*emission_line_model.spatial_shape, emission_line_model.neml),
                           dtype=numpy.float)
        data.reshape(-1, emission_line_model.neml)[indx,:] \
                = emission_line_model['EMLDATA'].data['KINERR'][unique_bins[reconstruct[indx]],:,0]
        pos = data > 0
        data[pos] = numpy.square(1.0/data[pos])
        data[numpy.invert(pos)] = 0.0
        hdus += [ fits.ImageHDU(data=data.T,
                                header=self._finalize_map_header(hdr, 'EMLINE_GVEL',
                                                                 hduclas2='ERROR', qual=True,
                                                                 bunit='(km/s)^{-2}'),
                                name='EMLINE_GVEL_IVAR') ]
        # Bitmask
        bit_type = emission_line_model.bitmask.minimum_uint_dtype()
        data = numpy.zeros((*emission_line_model.spatial_shape, emission_line_model.neml),
                           dtype=bit_type)
        data.reshape(-1, emission_line_model.neml)[indx,:] \
                = emission_line_model['EMLDATA'].data['MASK'][unique_bins[reconstruct[indx]],:]
        hdus += [ fits.ImageHDU(data=data.T,
                                header=self._finalize_map_header(hdr, 'EMLINE_GVEL',
                                                                 hduclas2='QUALITY', err=True,
                                                                 bit_type=bit_type),
                                name='EMLINE_GVEL_MASK') ]

        # Velocity dispersion
        data = numpy.zeros((*emission_line_model.spatial_shape, emission_line_model.neml),
                           dtype=numpy.float)
        data.reshape(-1, emission_line_model.neml)[indx,:] \
                = emission_line_model['EMLDATA'].data['KIN'][unique_bins[reconstruct[indx]],:,1]
        hdus += [ fits.ImageHDU(data=data.T,
                                header=self._finalize_map_header(hdr, 'EMLINE_GSIGMA', err=True,
                                                                 qual=True,
                                                                 bunit='km/s'),
                                name='EMLINE_GSIGMA') ]
        # Inverse variance
        data = numpy.zeros((*emission_line_model.spatial_shape, emission_line_model.neml),
                           dtype=numpy.float)
        data.reshape(-1, emission_line_model.neml)[indx,:] \
                = emission_line_model['EMLDATA'].data['KINERR'][unique_bins[reconstruct[indx]],:,1]
        pos = data > 0
        data[pos] = numpy.square(1.0/data[pos])
        data[numpy.invert(pos)] = 0.0
        hdus += [ fits.ImageHDU(data=data.T,
                                header=self._finalize_map_header(hdr, 'EMLINE_GSIGMA',
                                                                 hduclas2='ERROR', qual=True,
                                                                 bunit='(km/s)^{-2}'),
                                name='EMLINE_GSIGMA_IVAR') ]
        # Bitmask
        bit_type = emission_line_model.bitmask.minimum_uint_dtype()
        data = numpy.zeros((*emission_line_model.spatial_shape, emission_line_model.neml),
                           dtype=bit_type)
        data.reshape(-1, emission_line_model.neml)[indx,:] \
                = emission_line_model['EMLDATA'].data['MASK'][unique_bins[reconstruct[indx]],:]
        hdus += [ fits.ImageHDU(data=data.T,
                                header=self._finalize_map_header(hdr, 'EMLINE_GSIGMA',
                                                                 hduclas2='QUALITY', err=True,
                                                                 bit_type=bit_type),
                                name='EMLINE_GSIGMA_MASK') ]
        return hdus


    def spectral_index_maps(self, prihdr, spectral_indices):
        """
        Transposes are necessary to keep x,y orientation consistent with
        DRP output.
        """

        # Construct and return the empty hdus
        if spectral_indices is None:
            hdus = [ fits.ImageHDU(data=None,
                                   header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                    'SPECINDEX', err=True,
                                                                    qual=True),
                                   name='SPECINDEX') ]
            hdus += [ fits.ImageHDU(data=None,
                                    header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                     'SPECINDEX', hduclas2='ERROR',
                                                                     qual=True),
                                    name='SPECINDEX_IVAR') ]
            hdus += [ fits.ImageHDU(data=None,
                                    header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                     'SPECINDEX',
                                                                     hduclas2='QUALITY', err=True,
                                                                     bit_type=numpy.bool),
                                    name='SPECINDEX_MASK') ]
            return hdus

        # Add data to the primary header
        prihdr['SIKEY'] = (spectral_indices.database['key'],
                                'Spectral-index database keyword')
        prihdr['SIMINSN'] = (spectral_indices.database['minimum_snr'],
                                'Minimum S/N of spectrum to include')
        prihdr['FWHM'] = (spectral_indices.database['fwhm'],
                                'FWHM of index system resolution (ang)')
        prihdr['ARTDB'] = (spectral_indices.database['artifacts'],
                                'Artifact database keyword')
        prihdr['ABSDB'] = (spectral_indices.database['absindex'],
                                'Absorption-index database keyword')
        prihdr['BHDDB'] = (spectral_indices.database['bandhead'],
                                'Bandhead-index database keyword')

        # Bin index
        bin_indx = spectral_indices.binned_spectra['BINID'].data.copy().ravel()

        # Determine how to reconstruct the image data
        unique_bins, reconstruct = numpy.unique(bin_indx, return_inverse=True)
        indx = bin_indx > -1

        # Channel identifiers and units
        hdr = self._add_channel_names(self.multichannel_maphdr,
                                      spectral_indices['SIPAR'].data['NAME'],
                                      units=spectral_indices['SIPAR'].data['UNIT'])
        errunits = [ '({0})'.format(u)+'^{-2}' if len(u) > 0 else ''
                            for u in spectral_indices['SIPAR'].data['UNIT'] ]
        errhdr = self._add_channel_names(self.multichannel_maphdr,
                                         spectral_indices['SIPAR'].data['NAME'], units=errunits)
        maskhdr = self._add_channel_names(self.multichannel_maphdr,
                                          spectral_indices['SIPAR'].data['NAME'])

        # Index values
        data = numpy.zeros((*spectral_indices.spatial_shape, spectral_indices.nindx),
                           dtype=numpy.float)
        data.reshape(-1, spectral_indices.nindx)[indx,:] \
                = spectral_indices['SINDX'].data['INDX'][unique_bins[reconstruct[indx]],:]
        hdus = [ fits.ImageHDU(data=data.T,
                               header=self._finalize_map_header(hdr, 'SPECINDEX', err=True,
                                                                qual=True),
                               name='SPECINDEX') ]
        # Inverse variance
        data = numpy.zeros((*spectral_indices.spatial_shape, spectral_indices.nindx),
                           dtype=numpy.float)
        data.reshape(-1, spectral_indices.nindx)[indx,:] \
                = spectral_indices['SINDX'].data['INDXERR'][unique_bins[reconstruct[indx]],:]
        pos = data > 0
        data[pos] = numpy.square(1.0/data[pos])
        data[numpy.invert(pos)] = 0.0
        hdus += [ fits.ImageHDU(data=data.T,
                                header=self._finalize_map_header(errhdr, 'SPECINDEX',
                                                                 hduclas2='ERROR', qual=True),
                                name='SPECINDEX_IVAR') ]
        # Bitmask
        bit_type = spectral_indices.bitmask.minimum_uint_dtype()
        data = numpy.zeros((*spectral_indices.spatial_shape, spectral_indices.nindx),
                           dtype=bit_type)
        data.reshape(-1, spectral_indices.nindx)[indx,:] \
                = spectral_indices['SINDX'].data['MASK'][unique_bins[reconstruct[indx]],:]
        hdus += [ fits.ImageHDU(data=data.T,
                                header=self._finalize_map_header(hdr, 'SPECINDEX',
                                                                 hduclas2='QUALITY', err=True,
                                                                 bit_type=bit_type),
                                name='SPECINDEX_MASK') ]

        return hdus




