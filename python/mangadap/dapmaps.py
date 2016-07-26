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
    | **26 Jul 2016**: (KBW) Flag the Gaussian-fitted flux as unreliable
        if the summed flux is not within a factor of two

.. todo::
    Do something different than add an empty extension when the data are
    not available?
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
import astropy.constants

from .util.bitmask import BitMask
from .util.log import log_output
from .util.fileio import write_hdu
from .drpfits import DRPFits, DRPQuality3DBitMask
from .dapqual import DAPQualityBitMask
from .config.defaults import default_dap_method, default_dap_method_path, default_dap_file_name
from .config.defaults import default_dap_source
from .proc.reductionassessments import ReductionAssessment
from .proc.spatiallybinnedspectra import SpatiallyBinnedSpectra
from .proc.stellarcontinuummodel import StellarContinuumModel
from .proc.emissionlinemoments import EmissionLineMoments
from .proc.emissionlinemodel import EmissionLineModel
from .proc.spectralindices import SpectralIndices

from matplotlib import pyplot

__author__ = 'Kyle Westfall'

class DAPMapsBitMask(BitMask):
    """
    .. todo::
        - Force read IDLUTILS version as opposed to internal one?
    """
    def __init__(self, dapsrc=None):
        dapsrc = default_dap_source() if dapsrc is None else str(dapsrc)
        BitMask.__init__(self, ini_file=os.path.join(dapsrc, 'python', 'mangadap', 'config',
                                                     'bitmasks', 'dap_maps_bits.ini'))


class construct_maps_file:
    """
    Construct a DAP maps file based on the input data.

    Set as a class for coherency reasons, but should not be used as an
    object!

    Should force all intermediate objects to be provided.

    """
    def __init__(self, drpf, rdxqa=None, binned_spectra=None, stellar_continuum=None,
                 emission_line_moments=None, emission_line_model=None, spectral_indices=None,
                 nsa_redshift=None, dapsrc=None, dapver=None, analysis_path=None,
                 directory_path=None, output_file=None, clobber=True, loggers=None, quiet=False):

        # The output method directory is, for now, the combination of
        # the binned_spectrum and stellar_continuum method keys
        if directory_path is None and (binned_spectra is None or stellar_continuum is None):
            raise ValueError('Could not define output directory path.')

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

        # Initialize the reporting
        self.loggers = None if loggers is None else loggers
        self.quiet = quiet

        self.drpf = drpf
        self.spatial_shape = self.drpf.spatial_shape

        self.nsa_redshift = nsa_redshift

        # Set the output directory path
        # TODO: Get DAP version from __version__ string
        self.method = default_dap_method(binned_spectra=binned_spectra,
                                         stellar_continuum=stellar_continuum) \
                                if directory_path is None else None
        self.directory_path = default_dap_method_path(self.method, plate=self.drpf.plate,
                                                      ifudesign=self.drpf.ifudesign,
                                                      drpver=self.drpf.drpver, dapver=dapver,
                                                      analysis_path=analysis_path) \
                                                if directory_path is None else str(directory_path)

        # Set the output file
        self.output_file = default_dap_file_name(self.drpf.plate, self.drpf.ifudesign,
                                                 self.method, mode='MAPS') \
                                    if output_file is None else str(output_file)

        # Report
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, '-'*50)
            log_output(self.loggers, 1, logging.INFO, 'CONSTRUCTING OUTPUT MAPS:')
            log_output(self.loggers, 1, logging.INFO, '-'*50)
            log_output(self.loggers, 1, logging.INFO, 'Output path: {0}'.format(
                                                                            self.directory_path))
            log_output(self.loggers, 1, logging.INFO, 'Output file: {0}'.format(
                                                                            self.output_file))
            log_output(self.loggers, 1, logging.INFO, 'Output maps have shape {0}'.format(
                                                                            self.spatial_shape))
            log_output(self.loggers, 1, logging.INFO, 'NSA redshift: {0}'.format(self.nsa_redshift))

        ofile = os.path.join(self.directory_path, self.output_file)
        if os.path.isfile(ofile) and not clobber:
            # TODO: Perform some checks to make sure the existing file
            # has the correct content?
            warnings.warn('Output file exists!  Set clobber=True to overwrite.')
            return

        # Initialize the primary header
        prihdr = self._initialize_primary_header()

        # Get the base map headers
        self.multichannel_maphdr = self._map_header(nchannels=2)
        self.singlechannel_maphdr = self._map_header(nchannels=1)

        # Initialize the pixel mask
        self.bitmask = DAPMapsBitMask(dapsrc=dapsrc)
        self.common_mask = self._initialize_mask(binned_spectra)

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
        prihdr = self._finalize_primary_header(prihdr, binned_spectra, dapsrc=dapsrc)
        self.hdu = fits.HDUList([ fits.PrimaryHDU(header=prihdr),
                                  *self.rdxqalist,
                                  *self.bspeclist,
                                  *self.contlist,
                                  *self.emlmomlist,
                                  *self.elmodlist,
                                  *self.sindxlist
                                ])

        #---------------------------------------------------------------
        # TEMPORARY FLAGS:
        # Flag the Gaussian-fitted flux as unreliable if the summed flux
        # is not within a factor of two.
        factor = 5.0
        indx = (self.hdu['BINID'].data > -1) \
                & (self.hdu['EMLINE_GFLUX_MASK'].data == 0) \
                & (self.hdu['EMLINE_SFLUX_MASK'].data == 0) \
                & ( (self.hdu['EMLINE_GFLUX'].data < self.hdu['EMLINE_SFLUX'].data/factor)
                    | (self.hdu['EMLINE_GFLUX'].data > self.hdu['EMLINE_SFLUX'].data*factor) )
        print('unreliable Gaussian flux compared to summed flux: ', numpy.sum(indx))
        if numpy.sum(indx):
            self.hdu['EMLINE_GFLUX_MASK'].data[indx] \
                    = self.bitmask.turn_on(self.hdu['EMLINE_GFLUX_MASK'].data[indx], 'UNRELIABLE')
            self.hdu['EMLINE_GVEL_MASK'].data[indx] \
                    = self.bitmask.turn_on(self.hdu['EMLINE_GVEL_MASK'].data[indx], 'UNRELIABLE')
            self.hdu['EMLINE_GSIGMA_MASK'].data[indx] \
                    = self.bitmask.turn_on(self.hdu['EMLINE_GSIGMA_MASK'].data[indx], 'UNRELIABLE')

        # Flag the stellar velocity dispersions measured in spectra with
        # S/N<10 as unreliable
        indx = (self.hdu['BINID'].data > -1) & (self.hdu['BIN_SNR'].data < 10) \
                    & (self.hdu['STELLAR_SIGMA_MASK'].data == 0)
        print('unreliable sigma because of low S/N: ', numpy.sum(indx))
        if numpy.sum(indx):
            self.hdu['STELLAR_SIGMA_MASK'].data[indx] \
                    = self.bitmask.turn_on(self.hdu['STELLAR_SIGMA_MASK'].data[indx], 'UNRELIABLE')
        #---------------------------------------------------------------

        # Write the file
        if not os.path.isdir(self.directory_path):
            os.makedirs(self.directory_path)
        write_hdu(self.hdu, ofile, clobber=clobber, checksum=True, loggers=self.loggers,
                  quiet=self.quiet)
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, '-'*50)


    @staticmethod
    def _convert_to_Newtonian_velocity(v, redshift, ivar=False):
        if ivar:
            return v*numpy.square(1.0+redshift)
        return (v-astropy.constants.c.to('km/s').value*redshift)/(1.0+redshift)


    @staticmethod
    def _clean_primary_header(hdr):
        # Remove some keys that are incorrect for DAP data
        hdr.remove('BSCALE')
        hdr.remove('BZERO')
        hdr.remove('BUNIT')
#        hdr.remove('MASKNAME')


    def _initialize_primary_header(self):
        hdr = self.drpf.hdu['PRIMARY'].header.copy()

        # Clean out header for info not pertinent to the DAP
        self._clean_primary_header(hdr)

        # Change MASKNAME
        hdr['MASKNAME'] = 'MANGA_DAPPIXMASK'

        # Add versioning
        return hdr


    def _finalize_primary_header(self, prihdr, binned_spectra, dapsrc=None):

        # Initialize the DAP quality flag
        self.dapqualbm = DAPQualityBitMask(dapsrc=dapsrc)
        self.drp3qualbm = DRPQuality3DBitMask()
        self.dapqual = self.dapqualbm.minimum_dtype()(0)    # Type casting original flag to 0
        if self.drp3qualbm.flagged(self.drpf['PRIMARY'].header['DRP3QUAL'], flag='CRITICAL'):
            if not self.quiet:
                log_output(self.loggers, 1, logging.INFO, 'DRP File is flagged as CRITICAL!')
            self.dapqual = self.dapqualbm.turn_on(self.dapqual, 'CRITICAL')
            self.dapqual = self.dapqualbm.turn_on(self.dapqual, 'DRPCRIT')

        #Excluded from 2.0 tag to allow for dependency of DAP to be on
        #IDLUTILS v5_5_25
        if binned_spectra is not None and binned_spectra.method['binclass'] is not None \
                and binned_spectra.method['binclass'].bintype == 'voronoi' \
                and binned_spectra.nbins == 1:
            self.dapqual = self.dapqualbm.turn_on(self.dapqual, 'SINGLEBIN')

        # Determine if there's a foreground star
        if numpy.sum(self.drpf.bitmask.flagged(self.drpf['MASK'].data, flag='FORESTAR')) > 0:
            self.dapqual = self.dapqualbm.turn_on(self.dapqual, 'FORESTAR')

        # Commit the quality flag to the header
        prihdr['DAPQUAL'] = (self.dapqual, 'DAP quality bitmask')

        # Finalize authors
        prihdr['AUTHOR'] = 'K Westfall, B Andrews <kyle.westfall@port.co.uk, andrewsb@pitt.edu>'

        return prihdr



    def _initialize_mask(self, binned_spectra=None):
        """
        For each flag that is copied from the DRP3PIXMASK, find which
        pixels are masked in the cube, and flag the same bits using
        DAPMapsBitMask all the pixels in a given spectrum are masked.
        """
        copied_drp_pix_flags = [ 'NOCOV', 'LOWCOV', 'DEADFIBER', 'FORESTAR', 'DONOTUSE' ]
        mask = numpy.zeros(self.spatial_shape, dtype=self.bitmask.minimum_dtype())
        for flag in copied_drp_pix_flags:
            drpmask = numpy.all(self.drpf.bitmask.flagged(self.drpf['MASK'].data, flag=flag),
                                axis=2)
            if numpy.sum(drpmask) == 0:
                continue
            mask[drpmask] = self.bitmask.turn_on(mask[drpmask], flag)

        if binned_spectra is None:
            return mask

#        mask = self._binning_mask(mask, binned_spectra)
#        pyplot.imshow(mask, origin='lower', interpolation='nearest')
#        pyplot.show()
        return self._binning_mask(mask, binned_spectra)


    def _binning_mask(self, _mask, binned_spectra):
        """
        Copy relevant pixel masks from the binned spectra to the map
        pixels.

        LOW_SPECCOV, LOW_SNR, NONE_IN_STACK are all consolidated into NOVALUE

        """
        binned_spectra_flags = [ 'LOW_SPECCOV', 'LOW_SNR', 'NONE_IN_STACK' ]
        mask = _mask.copy()
        for flag in binned_spectra_flags:
            bsmask = numpy.all(binned_spectra.bitmask.flagged(binned_spectra['MASK'].data,
                                                              flag=flag), axis=2)
            if numpy.sum(bsmask) == 0:
                continue
            mask[bsmask] = self.bitmask.turn_on(mask[bsmask], 'NOVALUE')

#        pyplot.imshow(mask, origin='lower', interpolation='nearest')
#        pyplot.show()
        return self._consolidate_donotuse(mask)


    def _consolidate_donotuse(self, mask):
        flgd = self.bitmask.flagged(mask, flag=[ 'LOWCOV', 'FORESTAR', 'NOVALUE', 'FITFAILED',
                                                 'NEARBOUND', 'MATHERROR' ])
        mask[flgd] = self.bitmask.turn_on(mask[flgd], 'DONOTUSE')
        return mask


    @staticmethod
    def _clean_map_header(hdr, nchannels=1):

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
        # Change the pixel mask name
        hdr['MASKNAME'] = 'MANGA_DAPPIXMASK'
        # Add versioning
        # Add DAP quality
        # Other meta data?
        return hdr


    @staticmethod
    def _mask_data_type(bit_type):
        if bit_type in [numpy.uint64, numpy.int64]:
            return ('FLAG64BIT', '64-bit mask')
        if bit_type in [numpy.uint32, numpy.int32]:
            return ('FLAG32BIT', '32-bit mask')
        if bit_type in [numpy.uint16, numpy.int16]:
            return ('FLAG16BIT', '16-bit mask')
        if bit_type in [numpy.uint8, numpy.int8]:
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

        .. todo::
            Add the spaxel correlation data.

        """

        # Add data to the primary header
        if rdxqa is not None:
            prihdr = rdxqa._initialize_header(prihdr)
#            prihdr['RDXQAKEY'] = (rdxqa.method['key'])
#            prihdr['PA'] = (rdxqa.pa, 'Isophotal position angle')
#            prihdr['ELL'] = (rdxqa.ell, 'Isophotal ellipticity (1-b/a)')

        # On-sky coordinates
        hdr = self._add_channel_names(self.multichannel_maphdr, ['On-sky X', 'On-sky Y'])
        data = None if rdxqa is None else \
                    rdxqa['SPECTRUM'].data['SKY_COO'].reshape(*self.drpf.spatial_shape, -1).T
        hdus = [ fits.ImageHDU(data=data,
                               header=self._finalize_map_header(hdr, 'SPX_SKYCOO', bunit='arcsec',
                                                                nchannels=2),
                               name='SPX_SKYCOO') ]

        # Elliptical coordinates
        hdr = self._add_channel_names(self.multichannel_maphdr,
                                      ['Elliptical radius', 'Elliptical azimuth'],
                                      units=['arcsec', 'degrees'])
        data = None if rdxqa is None else \
                    rdxqa['SPECTRUM'].data['ELL_COO'].reshape(*self.drpf.spatial_shape, -1).T
        hdus += [ fits.ImageHDU(data=data,
                                header=self._finalize_map_header(hdr, 'SPX_ELLCOO', bunit='arcsec',
                                                                 nchannels=2),
                                name='SPX_ELLCOO') ]

#        # Spectral coverage
#        data = None if rdxqa is None else \
#                    rdxqa['SPECTRUM'].data['FGOODPIX'].reshape(self.drpf.spatial_shape).T
#        hdus += [ fits.ImageHDU(data=data,
#                                header=self._finalize_map_header(self.singlechannel_maphdr,
#                                                                 'SPX_SPCOV'),
#                                name='SPX_SPCOV') ]

        # SNR assessments of the DRP data
        data = None if rdxqa is None else \
                    rdxqa['SPECTRUM'].data['SIGNAL'].reshape(self.drpf.spatial_shape).T
        hdus += [ fits.ImageHDU(data=data,
                                header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                 'SPX_MFLUX',
                                                                bunit='1E-17 erg/s/cm^2/ang/spaxel',
                                                                 err=True),
                                name='SPX_MFLUX') ]

        if rdxqa is None:
            ivar = None
        else:
            ivar = rdxqa['SPECTRUM'].data['VARIANCE'].reshape(self.drpf.spatial_shape).T
            indx = ivar > 0
            ivar[indx] = 1.0/ivar[indx]
            ivar[numpy.invert(indx)] = 0.0
        hdus += [ fits.ImageHDU(data=data,
                                header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                 'SPX_MFLUX', hduclas2='ERROR',
                                                        bunit='(1E-17 erg/s/cm^2/ang/spaxel)^{-2}'),
                                name='SPX_MFLUX_IVAR') ]
        data = None if rdxqa is None else \
                    rdxqa['SPECTRUM'].data['SNR'].reshape(self.drpf.spatial_shape).T
        hdus += [ fits.ImageHDU(data=data,
                                header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                 'SPX_SNR'),
                                name='SPX_SNR') ]

        return hdus

#        pyplot.imshow(radius.T, origin='lower', interpolation='nearest')
#        pyplot.show()
#        pyplot.imshow(self.drpf['FLUX'].data[:,:,1000].T, origin='lower', interpolation='nearest')
#        pyplot.show()

#    def binning_mask(self, binned_spectra):
#        """
#        Return the basic binning masks:
#            


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
                                                                     err=True, qual=True),
                                    name='BIN_MFLUX') ]
            hdus += [ fits.ImageHDU(data=None,
                                    header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                     'BIN_MFLUX', hduclas2='ERROR',
                                                        bunit='(1E-17 erg/s/cm^2/ang/spaxel)^{-2}',
                                                                     qual=True),
                                    name='BIN_MFLUX_IVAR') ]
            hdus += [ fits.ImageHDU(data=None,
                                    header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                     'BIN_MFLUX',
                                                                     hduclas2='QUALITY', err=True,
                                                                     bit_type=numpy.bool),
                                    name='BIN_MFLUX_MASK') ]
            hdus += [ fits.ImageHDU(data=None,
                                    header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                     'BIN_SNR'),
                                    name='BIN_SNR') ]

            return hdus

        # Add data to the primary header
        # TODO: Apply this to each extension instead of the primary
        # header to allow for extension specific binning?
        prihdr = binned_spectra._initialize_header(prihdr)
        if binned_spectra.is_unbinned:
            prihdr['BINTYPE'] = ('None', 'Binning method')
        else:
            prihdr = binned_spectra._add_method_header(prihdr)
        prihdr = binned_spectra._get_reddening_hdr(prihdr)

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
                                      units=['arcsec', 'degrees'])
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
                                                                 err=True, qual=True),
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
                                                                 qual=True,
                                                        bunit='(1E-17 erg/s/cm^2/ang/spaxel)^{-2}'),
                                name='BIN_MFLUX_IVAR') ]

        # Bitmask
        mask = self.common_mask.copy()
        hdus += [ fits.ImageHDU(data=mask.T,
                                header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                 'BIN_MFLUX',
                                                                 hduclas2='QUALITY', err=True,
                                                                 bit_type=mask.dtype.type),
                                name='BIN_MFLUX_MASK') ]

        # Signal-to-noise in the bin
        data = numpy.zeros(binned_spectra.spatial_shape, dtype=numpy.float)
        data.ravel()[indx] = binned_spectra['BINS'].data['SNR'][unique_bins[reconstruct[indx]]]
        hdus += [ fits.ImageHDU(data=data.T,
                                header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                 'BIN_SNR'),
                                name='BIN_SNR') ]

        return hdus


    def _stellar_continuum_mask_to_map_mask(self, stellar_continuum, unique_bins, reconstruct,
                                            indx, for_dispersion=False):
        """
        Propagate and consolidate the stellar-continuum masks to the map
        pixel mask.

        DIDNOTUSE, FORESTAR propagated from already existing mask
            (self.common_mask)

        LOW_SNR, NO_FIT, INSUFFICIENT_DATA propagated to NOVALUE

        FIT_FAILED, NEAR_BOUND propagated to FITFAILED and NEARBOUND

        NEGATIVE_WEIGHTS consolidated into UNRELIABLE; if
            dispersion=True, include BAD_SIGMA in this
        """

        # Construct the existing mask data
        bit_type = stellar_continuum.bitmask.minimum_dtype()
        output_shape = stellar_continuum.spatial_shape
        sc_mask = numpy.zeros(output_shape, dtype=bit_type)
        sc_mask.ravel()[indx] \
                = stellar_continuum['PAR'].data['MASK'][unique_bins[reconstruct[indx]]]

        # Copy the common mask
        mask = self.common_mask.copy()

        # Find the low S/N bins
        snmask = numpy.all(stellar_continuum.bitmask.flagged(stellar_continuum['MASK'].data,
                                                             flag='LOW_SNR'), axis=2)
        if numpy.sum(snmask) > 0:
            mask[snmask] = self.bitmask.turn_on(mask[snmask], 'NOVALUE')

        # Consolidate to NOVALUE
        flgd = stellar_continuum.bitmask.flagged(sc_mask, flag=['NO_FIT', 'INSUFFICIENT_DATA' ])
        mask[flgd] = self.bitmask.turn_on(mask[flgd], 'NOVALUE')

        # Copy FIT_FAILED
        flgd = stellar_continuum.bitmask.flagged(sc_mask, flag='FIT_FAILED')
        mask[flgd] = self.bitmask.turn_on(mask[flgd], 'FITFAILED')

        # Copy NEAR_BOUND
        flgd = stellar_continuum.bitmask.flagged(sc_mask, flag='NEAR_BOUND')
        mask[flgd] = self.bitmask.turn_on(mask[flgd], 'NEARBOUND')

        # Consolidate to UNRELIABLE
        flgd = stellar_continuum.bitmask.flagged(sc_mask, flag='NEGATIVE_WEIGHTS')
        mask[flgd] = self.bitmask.turn_on(mask[flgd], 'UNRELIABLE')
        if for_dispersion:
            flgd = stellar_continuum.bitmask.flagged(sc_mask, flag='BAD_SIGMA')
            mask[flgd] = self.bitmask.turn_on(mask[flgd], 'UNRELIABLE')

        return self._consolidate_donotuse(mask)


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
            hdus += [ fits.ImageHDU(data=None,
                                    header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                     'STELLAR_SIGMACORR'),
                                    name='STELLAR_SIGMACORR') ]
            hdus += [ fits.ImageHDU(data=None,
                                    header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                     'STELLAR_CONT_FRESID'),
                                    name='STELLAR_CONT_FRESID') ]
            hdus += [ fits.ImageHDU(data=None,
                                    header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                     'STELLAR_CONT_RCHI2'),
                                    name='STELLAR_CONT_RCHI2') ]
            return hdus

        # Add data to the primary header
        prihdr = stellar_continuum._initialize_header(prihdr)
        prihdr = stellar_continuum._add_method_header(prihdr)

        # Bin index
        bin_indx = stellar_continuum['BINID'].data.copy().ravel()

        # Determine how to reconstruct the image data
        unique_bins, reconstruct = numpy.unique(bin_indx, return_inverse=True)
        indx = bin_indx > -1

        # Get the velocity mask
        vel_mask = self._stellar_continuum_mask_to_map_mask(stellar_continuum, unique_bins,
                                                            reconstruct, indx)
        # Add the bad sigma measurements
        sig_mask = self._stellar_continuum_mask_to_map_mask(stellar_continuum, unique_bins,
                                                            reconstruct, indx, for_dispersion=True)

        # Stellar velocity
        data = numpy.zeros(stellar_continuum.spatial_shape, dtype=numpy.float)
        data.ravel()[indx] = stellar_continuum['PAR'].data['KIN'][unique_bins[reconstruct[indx]],0]
        # Convert redshift to Newtonian velocity
        data = self._convert_to_Newtonian_velocity(data, self.nsa_redshift)
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
        data = self._convert_to_Newtonian_velocity(data, self.nsa_redshift, ivar=True)
        hdus += [ fits.ImageHDU(data=data.T,
                                header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                 'STELLAR_VEL', hduclas2='ERROR',
                                                                 bunit='(km/s)^{-2}', qual=True),
                                name='STELLAR_VEL_IVAR') ]
        # Bitmask
        hdus += [ fits.ImageHDU(data=vel_mask.T,
                                header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                 'STELLAR_VEL',
                                                                 hduclas2='QUALITY', err=True,
                                                                 bit_type=vel_mask.dtype.type),
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
        hdus += [ fits.ImageHDU(data=sig_mask.T,
                                header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                 'STELLAR_SIGMA',
                                                                 hduclas2='QUALITY', err=True,
                                                                 bit_type=sig_mask.dtype.type),
                                name='STELLAR_SIGMA_MASK') ]

        # Stellar velocity dispersion correction
        data = numpy.zeros(stellar_continuum.spatial_shape, dtype=numpy.float)
        data.ravel()[indx] \
                = stellar_continuum['PAR'].data['SIGMACORR'][unique_bins[reconstruct[indx]]]
        hdus += [ fits.ImageHDU(data=data.T,
                               header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                'STELLAR_SIGMACORR', bunit='km/s'),
                               name='STELLAR_SIGMACORR') ]

        # Continuum fit statistics
        s = numpy.array([False, True, False, True, False])
        data = numpy.zeros(stellar_continuum.spatial_shape+(2,), dtype=numpy.float)
        data.reshape(-1,2)[indx,:] \
                = stellar_continuum['PAR'].data['FABSRESID'][unique_bins[reconstruct[indx]],:][:,s]
        hdr = self._add_channel_names(self.multichannel_maphdr, ['68th percentile',
                                                                 '99th percentile'])
        hdus += [ fits.ImageHDU(data=data.T,
                                header=self._finalize_map_header(hdr, 'STELLAR_CONT_FRESID'),
                                name='STELLAR_CONT_FRESID') ]

        data = numpy.zeros(stellar_continuum.spatial_shape, dtype=numpy.float)
        data.ravel()[indx] = stellar_continuum['PAR'].data['RCHI2'][unique_bins[reconstruct[indx]]]
        hdus += [ fits.ImageHDU(data=data.T,
                                header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                 'STELLAR_CONT_RCHI2'),
                                name='STELLAR_CONT_RCHI2') ]

        return hdus


    def _emission_line_moment_mask_to_map_mask(self, emission_line_moments, unique_bins,
                                               reconstruct, indx):
        """
        Propagate and consolidate the emission-line moment masks to the map
        pixel mask.

        DIDNOTUSE, FORESTAR propagated from already existing mask (self.common_mask)

        MAIN_EMPTY, BLUE_EMPTY, RED_EMPTY, UNDEFINED_BANDS propagated to NOVALUE

        MAIN_JUMP, BLUE_JUMP, RED_JUMP, JUMP_BTWN_SIDEBANDS propagated to FITFAILED

        NO_ABSORPTION_CORRECTION propaged to NOCORRECTION

        DIVBYZERO propagated to MATHERROR

        Second moments not provided so no need to include UNDEFINED_MOM2

        Need to further assess *_INCOMP to see if these lead to bad values (BADVALUE).
        """

        # Construct the existing mask data
        bit_type = emission_line_moments.bitmask.minimum_dtype()
        output_shape = (*emission_line_moments.spatial_shape, emission_line_moments.nmom)
        elm_mask = numpy.zeros(output_shape, dtype=bit_type)
        elm_mask.reshape(-1, emission_line_moments.nmom)[indx,:] \
                = emission_line_moments['ELMMNTS'].data['MASK'][unique_bins[reconstruct[indx]],:]

        # Copy the common mask
        mask = numpy.array([self.common_mask.copy()]*emission_line_moments.nmom).transpose(1,2,0)

        # Consolidate to NOVALUE
        flgd = emission_line_moments.bitmask.flagged(elm_mask, flag=['MAIN_EMPTY', 'BLUE_EMPTY',
                                                                'RED_EMPTY', 'UNDEFINED_BANDS' ])
        mask[flgd] = self.bitmask.turn_on(mask[flgd], 'NOVALUE')

        # Consolidate to FITFAILED
        flgd = emission_line_moments.bitmask.flagged(elm_mask, flag=['MAIN_JUMP', 'BLUE_JUMP',
                                                                 'RED_JUMP', 'JUMP_BTWN_SIDEBANDS'])
        mask[flgd] = self.bitmask.turn_on(mask[flgd], 'FITFAILED')

        # Consolidate to NOCORRECTION
        flgd = emission_line_moments.bitmask.flagged(elm_mask, flag=['NO_ABSORPTION_CORRECTION'])
        mask[flgd] = self.bitmask.turn_on(mask[flgd], 'NOCORRECTION')

        # Consolidate to MATHERROR
        flgd = emission_line_moments.bitmask.flagged(elm_mask, flag=['DIVBYZERO'])
        mask[flgd] = self.bitmask.turn_on(mask[flgd], 'MATHERROR')

        return self._consolidate_donotuse(mask)


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
            hdus += [ fits.ImageHDU(data=None,
                                    header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                     'EMLINE_SEW', err=True,
                                                                     qual=True, bunit='ang'),
                                    name='EMLINE_SEW') ]
            hdus += [ fits.ImageHDU(data=None,
                                    header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                     'EMLINE_SEW',
                                                                     bunit='(ang)^{-2}',
                                                                     hduclas2='ERROR', qual=True),
                                    name='EMLINE_SEW_IVAR') ]
            hdus += [ fits.ImageHDU(data=None,
                                    header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                     'EMLINE_SEW',
                                                                     hduclas2='QUALITY', err=True,
                                                                     bit_type=numpy.bool),
                                    name='EMLINE_SEW_MASK') ]
            return hdus

        # Add data to the primary header
        prihdr = emission_line_moments._initialize_header(prihdr)

        # Bin index
        bin_indx = emission_line_moments.binned_spectra['BINID'].data.copy().ravel()

        # Determine how to reconstruct the image data
        unique_bins, reconstruct = numpy.unique(bin_indx, return_inverse=True)
        indx = bin_indx > -1

        # The mask is common to all extensions
        mask = self._emission_line_moment_mask_to_map_mask(emission_line_moments, unique_bins,
                                                           reconstruct, indx)

        # Integrated flux
        data = numpy.zeros((*emission_line_moments.spatial_shape, emission_line_moments.nmom),
                           dtype=numpy.float)
        data.reshape(-1, emission_line_moments.nmom)[indx,:] \
                = emission_line_moments['ELMMNTS'].data['FLUX'][unique_bins[reconstruct[indx]],:]
        # TODO: Convert wavelengths to air for column names?
        cols = [ '{0}-{1}'.format(n,int(w)) \
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
        hdus += [ fits.ImageHDU(data=mask.T,
                                header=self._finalize_map_header(hdr, 'EMLINE_SFLUX',
                                                                 hduclas2='QUALITY', err=True,
                                                                 bit_type=mask.dtype.type),
                                name='EMLINE_SFLUX_MASK') ]

        # Equivalent width
        data = numpy.zeros((*emission_line_moments.spatial_shape, emission_line_moments.nmom),
                           dtype=numpy.float)
        data.reshape(-1, emission_line_moments.nmom)[indx,:] \
                = emission_line_moments['ELMMNTS'].data['EW'][unique_bins[reconstruct[indx]],:]
        # TODO: Convert wavelengths to air for column names?
        cols = [ '{0}-{1}'.format(n,int(w)) \
                        for n,w in zip(emission_line_moments['ELMBAND'].data['NAME'],
                                       emission_line_moments['ELMBAND'].data['RESTWAVE']) ]
        hdr = self._add_channel_names(self.multichannel_maphdr, cols)
        hdus += [ fits.ImageHDU(data=data.T,
                                header=self._finalize_map_header(hdr, 'EMLINE_SEW', err=True,
                                                                 qual=True, bunit='ang'),
                                name='EMLINE_SEW') ]
        # Inverse variance
        data = numpy.zeros((*emission_line_moments.spatial_shape, emission_line_moments.nmom),
                           dtype=numpy.float)
        data.reshape(-1, emission_line_moments.nmom)[indx,:] \
                = emission_line_moments['ELMMNTS'].data['EWERR'][unique_bins[reconstruct[indx]],:]
        pos = data > 0
        data[pos] = numpy.square(1.0/data[pos])
        data[numpy.invert(pos)] = 0.0
        hdus += [ fits.ImageHDU(data=data.T,
                                header=self._finalize_map_header(hdr, 'EMLINE_SEW',
                                                                 hduclas2='ERROR', qual=True,
                                                                 bunit='(ang)^{-2}'),
                                name='EMLINE_SEW_IVAR') ]
        # Bitmask
        hdus += [ fits.ImageHDU(data=mask.T,
                                header=self._finalize_map_header(hdr, 'EMLINE_SEW',
                                                                 hduclas2='QUALITY', err=True,
                                                                 bit_type=mask.dtype.type),
                                name='EMLINE_SEW_MASK') ]
        return hdus


    def _emission_line_model_mask_to_map_mask(self, emission_line_model, unique_bins, reconstruct,
                                              indx, dispersion=False):
        """
        Propagate and consolidate the emission-line moment masks to the map
        pixel mask.

        INSUFFICIENT_DATA propagated to NOVALUE

        FIT_FAILED, UNDEFINED_COVAR propagated to FITFAILED

        NEAR_BOUND copied to NEARBOUND

        UNDEFINED_SIGMA propagated to UNRELIABLE
        """

        # Construct the existing mask data
        bit_type = emission_line_model.bitmask.minimum_dtype()
        elf_mask = numpy.zeros((*emission_line_model.spatial_shape, emission_line_model.neml),
                               dtype=bit_type)
        elf_mask.reshape(-1, emission_line_model.neml)[indx,:] \
                = emission_line_model['EMLDATA'].data['MASK'][unique_bins[reconstruct[indx]],:]

        # Copy the common mask
        mask = numpy.array([self.common_mask.copy()]*emission_line_model.neml).transpose(1,2,0)

        # Consolidate to NOVALUE
        flgd = emission_line_model.bitmask.flagged(elf_mask, flag=['INSUFFICIENT_DATA'])
        mask[flgd] = self.bitmask.turn_on(mask[flgd], 'NOVALUE')

        # Consolidate to FITFAILED
        flgd = emission_line_model.bitmask.flagged(elf_mask, flag=['FIT_FAILED', 'UNDEFINED_COVAR'])
        mask[flgd] = self.bitmask.turn_on(mask[flgd], 'FITFAILED')

        # Copy NEAR_BOUND
        flgd = emission_line_model.bitmask.flagged(elf_mask, flag='NEAR_BOUND')
        mask[flgd] = self.bitmask.turn_on(mask[flgd], 'NEARBOUND')

        # Consolidate to UNRELIABLE; currently only done for the dispersion
        if dispersion:
            flgd = emission_line_model.bitmask.flagged(elf_mask, flag=['UNDEFINED_SIGMA'])
            mask[flgd] = self.bitmask.turn_on(mask[flgd], 'UNRELIABLE')

        return self._consolidate_donotuse(mask)


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
            hdus += [ fits.ImageHDU(data=None,
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
            hdus += [ fits.ImageHDU(data=None,
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
            hdus += [ fits.ImageHDU(data=None,
                                   header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                    'EMLINE_GSIGMACORR',
                                                                    bunit='km/s'),
                                   name='EMLINE_GSIGMACORR') ]
            return hdus

        # Add data to the primary header
        prihdr = emission_line_model._initialize_header(prihdr)
        prihdr = emission_line_model._add_method_header(prihdr)

        # Bin index
        bin_indx = emission_line_model.binned_spectra['BINID'].data.copy().ravel()

        # Determine how to reconstruct the image data
        unique_bins, reconstruct = numpy.unique(bin_indx, return_inverse=True)
        indx = bin_indx > -1

        # The extensions share a common mask, except for the dispersion
        mask = self._emission_line_model_mask_to_map_mask(emission_line_model, unique_bins,
                                                          reconstruct, indx)

        # Channel identifiers (units are all the same)
        # TODO: Convert wavelengths to air for column names?
        cols = [ '{0}-{1}'.format(n,int(w)) \
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
        # TODO: Need to propagate bad inverse variances to mask
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
        hdus += [ fits.ImageHDU(data=mask.T,
                                header=self._finalize_map_header(hdr, 'EMLINE_GFLUX',
                                                                 hduclas2='QUALITY', err=True,
                                                                 bit_type=mask.dtype.type),
                                name='EMLINE_GFLUX_MASK') ]

        # Velocity
        data = numpy.zeros((*emission_line_model.spatial_shape, emission_line_model.neml),
                           dtype=numpy.float)
        data.reshape(-1, emission_line_model.neml)[indx,:] \
                = emission_line_model['EMLDATA'].data['KIN'][unique_bins[reconstruct[indx]],:,0]
        data = self._convert_to_Newtonian_velocity(data, self.nsa_redshift)
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
        data = self._convert_to_Newtonian_velocity(data, self.nsa_redshift, ivar=True)
        hdus += [ fits.ImageHDU(data=data.T,
                                header=self._finalize_map_header(hdr, 'EMLINE_GVEL',
                                                                 hduclas2='ERROR', qual=True,
                                                                 bunit='(km/s)^{-2}'),
                                name='EMLINE_GVEL_IVAR') ]
        # Bitmask
        hdus += [ fits.ImageHDU(data=mask.T,
                                header=self._finalize_map_header(hdr, 'EMLINE_GVEL',
                                                                 hduclas2='QUALITY', err=True,
                                                                 bit_type=mask.dtype.type),
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
        mask = self._emission_line_model_mask_to_map_mask(emission_line_model, unique_bins,
                                                          reconstruct, indx, dispersion=True)
        hdus += [ fits.ImageHDU(data=mask.T,
                                header=self._finalize_map_header(hdr, 'EMLINE_GSIGMA',
                                                                 hduclas2='QUALITY', err=True,
                                                                 bit_type=mask.dtype.type),
                                name='EMLINE_GSIGMA_MASK') ]

        # Instrumental dispersion
        data = numpy.zeros((*emission_line_model.spatial_shape, emission_line_model.neml),
                           dtype=numpy.float)
        data.reshape(-1, emission_line_model.neml)[indx,:] \
                = emission_line_model['EMLDATA'].data['SINST'][unique_bins[reconstruct[indx]],:]
        hdus += [ fits.ImageHDU(data=data.T,
                                header=self._finalize_map_header(hdr, 'EMLINE_INSTSIGMA',
                                                                 bunit='km/s'),
                                name='EMLINE_INSTSIGMA') ]

        return hdus


    def _spectral_index_mask_to_map_mask(self, spectral_indices, unique_bins, reconstruct, indx):
        """
        Propagate and consolidate the spectral index masks to the map
        pixel mask.

        DIDNOTUSE, FORESTAR propagated from already existing mask (self.common_mask)

        MAIN_EMPTY, BLUE_EMPTY, RED_EMPTY propagated to NOVALUE

        NO_DISPERSION_CORRECTION propagated to NOCORRECTION

        DIVBYZERO propagated to MATHERROR

        Need to further assess *_INCOMP to see if these lead to bad values (BADVALUE).
        """

        # Construct the existing mask data
        bit_type = spectral_indices.bitmask.minimum_dtype()
        si_mask = numpy.zeros((*spectral_indices.spatial_shape, spectral_indices.nindx),
                              dtype=bit_type)
        si_mask.reshape(-1, spectral_indices.nindx)[indx,:] \
                = spectral_indices['SINDX'].data['MASK'][unique_bins[reconstruct[indx]],:]

        # Copy the common mask
        mask = numpy.array([self.common_mask.copy()]*spectral_indices.nindx).transpose(1,2,0)

        # Consolidate to NOVALUE
        flgd = spectral_indices.bitmask.flagged(si_mask, flag=['MAIN_EMPTY', 'BLUE_EMPTY',
                                                                                    'RED_EMPTY' ])
        mask[flgd] = self.bitmask.turn_on(mask[flgd], 'NOVALUE')

        # Consolidate to NOCORRECTION
        flgd = spectral_indices.bitmask.flagged(si_mask, flag=['NO_DISPERSION_CORRECTION'])
        mask[flgd] = self.bitmask.turn_on(mask[flgd], 'NOCORRECTION')

        # Consolidate to MATHERROR
        flgd = spectral_indices.bitmask.flagged(si_mask, flag=['DIVBYZERO'])
        mask[flgd] = self.bitmask.turn_on(mask[flgd], 'MATHERROR')

        return self._consolidate_donotuse(mask)


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
            hdus += [ fits.ImageHDU(data=None,
                                    header=self._finalize_map_header(self.singlechannel_maphdr,
                                                                     'SPECINDEX_CORR'),
                                    name='SPECINDEX_CORR') ]
            return hdus

        if spectral_indices['PRIMARY'].header['SICORR']:
            raise ValueError('Cannot construct maps file with dispersion corrections applied to '
                             'indices!')

        # Add data to the primary header
        prihdr = spectral_indices._initialize_header(prihdr)

        # Bin index
        bin_indx = spectral_indices.binned_spectra['BINID'].data.copy().ravel()

        # Determine how to reconstruct the image data
        unique_bins, reconstruct = numpy.unique(bin_indx, return_inverse=True)
        indx = bin_indx > -1

        # The extensions share a common mask
        mask = self._spectral_index_mask_to_map_mask(spectral_indices, unique_bins, reconstruct,
                                                     indx)

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
        hdus += [ fits.ImageHDU(data=mask.T,
                                header=self._finalize_map_header(hdr, 'SPECINDEX',
                                                                 hduclas2='QUALITY', err=True,
                                                                 bit_type=mask.dtype.type),
                                name='SPECINDEX_MASK') ]

        # Index dispersion corrections
        data = numpy.zeros((*spectral_indices.spatial_shape, spectral_indices.nindx),
                           dtype=numpy.float)
        data.reshape(-1, spectral_indices.nindx)[indx,:] \
                = spectral_indices['SINDX'].data['INDX_DISPCORR'][unique_bins[reconstruct[indx]],:]
        hdus += [ fits.ImageHDU(data=data.T,
                                header=self._finalize_map_header(hdr, 'SPECINDEX_CORR'),
                                name='SPECINDEX_CORR') ]

        return hdus




