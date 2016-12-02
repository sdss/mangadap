# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""

Defines the class that constructs the DAP cube file output.

*License*:
    Copyright (c) 2016, SDSS-IV/MaNGA Pipeline Group
        Licensed under BSD 3-clause license - see LICENSE.rst

*Source location*:
    $MANGADAP_DIR/python/mangadap/dapcube.py

*Imports and python version compliance*:
    ::

*Class usage examples*:
    Add some usage comments here!

*Revision history*:
    | **13 Jul 2016**: Original Implementation by K. Westfall (KBW)

.. todo::
    - Do something different than add an empty extension when the data
      are not available?
    - Develop DAPFits to be a common interface between the maps and
      model cube files.
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

# For versioning
import scipy
import astropy
import pydl

from astropy.wcs import WCS
from astropy.io import fits

from .util.bitmask import BitMask
from .util.log import log_output
from .util.fileio import write_hdu
from .util.extinction import apply_reddening
from .mangafits import MaNGAFits
from .drpfits import DRPFits, DRPQuality3DBitMask
from .dapqual import DAPQualityBitMask
from .config.defaults import default_dap_method, default_dap_method_path, default_dap_file_name
from .config.defaults import default_dap_source
from .proc.spatiallybinnedspectra import SpatiallyBinnedSpectra
from .proc.stellarcontinuummodel import StellarContinuumModel
from .proc.emissionlinemodel import EmissionLineModel
from .dapmaps import DAPQualityBitMask, construct_maps_file

from matplotlib import pyplot

__author__ = 'Kyle Westfall'


class DAPCubeBitMask(BitMask):
    """
    .. todo::
        - Force read IDLUTILS version as opposed to internal one?
    """
    def __init__(self, dapsrc=None):
        dapsrc = default_dap_source() if dapsrc is None else str(dapsrc)
        BitMask.__init__(self, ini_file=os.path.join(dapsrc, 'python', 'mangadap', 'config',
                                                     'bitmasks', 'dap_cube_bits.ini'))


class construct_cube_file:
    """
    Construct a DAP cube file based on the input data.

    Set as a class for coherency reasons, but should not be used as an
    object!

    Should force all intermediate objects to be provided.

    """
    def __init__(self, drpf, binned_spectra=None, stellar_continuum=None, emission_line_model=None,
                 dapsrc=None, dapver=None, analysis_path=None, directory_path=None,
                 output_file=None, clobber=True, loggers=None, quiet=False):

        # The output method directory is, for now, the combination of
        # the binned_spectrum and stellar_continuum method keys
        if directory_path is None and (binned_spectra is None or stellar_continuum is None):
            raise ValueError('Could not define output directory path.')

        # Type checking
        if not isinstance(drpf, DRPFits):
            raise TypeError('Input must have type DRPFits.')
        if binned_spectra is not None and not isinstance(binned_spectra, SpatiallyBinnedSpectra):
            raise TypeError('Input must have type SpatiallyBinnedSpectra.')
        if stellar_continuum is not None and not isinstance(stellar_continuum,
                                                            StellarContinuumModel):
            raise TypeError('Input must have type StellarContinuumModel.')
        if emission_line_model is not None and not isinstance(emission_line_model,
                                                              EmissionLineModel):
            raise TypeError('Input must have type EmissionLineModel.')

        # Initialize the reporting
        self.loggers = None if loggers is None else loggers
        self.quiet = quiet

        self.drpf = drpf
        self.spatial_shape = self.drpf.spatial_shape

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
                                                 self.method, mode='LOGCUBE') \
                                    if output_file is None else str(output_file)

        # Report
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, '-'*50)
            log_output(self.loggers, 1, logging.INFO, 'CONSTRUCTING OUTPUT MODEL CUBE:')
            log_output(self.loggers, 1, logging.INFO, '-'*50)
            log_output(self.loggers, 1, logging.INFO, 'Output path: {0}'.format(
                                                                            self.directory_path))
            log_output(self.loggers, 1, logging.INFO, 'Output file: {0}'.format(
                                                                            self.output_file))

        ofile = os.path.join(self.directory_path, self.output_file)
        if os.path.isfile(ofile) and not clobber:
            # TODO: Perform some checks to make sure the existing file
            # has the correct content?
            warnings.warn('Output file exists!  Set clobber=True to overwrite.')
            return

        # Initialize the primary header
        prihdr = self._initialize_primary_header()

        # Get the base map header
        self.base_cubehdr = self._cube_header()

        # Initialize the pixel mask
        self.bitmask = DAPCubeBitMask(dapsrc=dapsrc)
        self.common_mask = self._initialize_mask(binned_spectra)

        # Construct the hdu list for each input object.
        # Binned and model spectra:
        prihdr, speclist, spec_specarr = self.model_and_data_cube(prihdr, binned_spectra,
                                                                  stellar_continuum,
                                                                  emission_line_model)
        # Emission-line only models:
        prihdr, elmodlist, elmod_specarr \
                = self.emission_line_model_cube(prihdr, binned_spectra, emission_line_model)

        # Save the data to the hdu attribute
        prihdr = self._finalize_primary_header(prihdr, binned_spectra, dapsrc=dapsrc)
#        lst = [ fits.PrimaryHDU(header=prihdr) ]
#        lst += speclist
#        lst += elmodlist
#        self.hdu = fits.HDUList(lst)
        self.hdu = fits.HDUList([ fits.PrimaryHDU(header=prihdr),
                                  *speclist,
                                  *elmodlist
                                ])

        # Restructure the cubes to match the DRP data
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, 'Restructuring data to match DRP.')
        self.spectral_arrays = numpy.append(spec_specarr, elmod_specarr)
        MaNGAFits.restructure_cube(self.hdu, ext=self.spectral_arrays, inverse=True)

        # Write the file
        if not os.path.isdir(self.directory_path):
            os.makedirs(self.directory_path)
        write_hdu(self.hdu, ofile, clobber=clobber, checksum=True, loggers=self.loggers,
                  quiet=self.quiet)
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, '-'*50)


    # Make common functionality with DAP Maps files
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
        hdr['MASKNAME'] = 'MANGA_DAPSPECMASK'

        # Add versioning
        hdr['VERSPY'] = ('.'.join([ str(v) for v in sys.version_info[:3]]), 'Python version')
        hdr['VERSNP'] = (numpy.__version__, 'Numpy version')
        hdr['VERSSCI'] = (scipy.__version__, 'Scipy version')
        hdr['VERSAST'] = (astropy.__version__, 'Astropy version')
        hdr['VERSPYDL'] = (pydl.__version__, 'pydl version')
        # TODO: Hard-coded! Should read DAP version programmatically
        hdr['VERSDAP'] = ('2.0.2', 'MaNGA DAP version')
        
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
        prihdr['AUTHOR'] = 'K Westfall, B Andrews <westfall@ucolick.org, andrewsb@pitt.edu>'

        return prihdr


    @staticmethod
    def _clean_cube_header(hdr):

        # Remove everything but the WCS information
        w = WCS(header=hdr)
        hdr = w.to_header().copy()

        # Fix the DATE-OBS keyword:
        hdr.comments['DATE-OBS'] = 'Date of median exposure'
        hdr.comments['MJD-OBS'] = '[d] MJD for DATE-OBS'

        # Add back in the BSCALE and BZERO values; BUNIT added during
        # "finalize"
        hdr['BSCALE'] = 1.
        hdr['BZERO'] = 0.

        return hdr


    def _cube_header(self):
        """
        .. todo::
            Add versioning!
        """
        hdr = self.drpf.hdu['FLUX'].header.copy()
        hdr = self._clean_cube_header(hdr)
        # Add Authors
        hdr['AUTHOR'] = 'K Westfall & B Andrews <kyle.westfall@ucolick.org, andrewsb@pitt.edu>'
        # Change the pixel mask name
        hdr['MASKNAME'] = 'MANGA_DAPSPECMASK'
        # Add versioning
        # Add DAP quality
        # Other meta data?
        return hdr


#    @staticmethod
#    def _mask_data_type(bit_type):
#        if bit_type in [numpy.uint64, numpy.int64]:
#            return ('FLAG64BIT', '64-bit mask')
#        if bit_type in [numpy.uint32, numpy.int32]:
#            return ('FLAG32BIT', '32-bit mask')
#        if bit_type in [numpy.uint16, numpy.int16]:
#            return ('FLAG16BIT', '16-bit mask')
#        if bit_type in [numpy.uint8, numpy.int8]:
#            return ('FLAG8BIT', '8-bit mask')
#        if bit_type == numpy.bool:
#            return ('MASKZERO', 'Binary mask; zero values are good/unmasked')
#        raise ValueError('Invalid bit_type: {0}!'.format(str(bit_type)))

    
    @staticmethod
    def _finalize_cube_header(hdr, ext, bunit=None, hduclas2='DATA', err=False, qual=False,
                              bit_type=None, prepend=True):

        _hdr = hdr.copy()

        if bunit is not None:
            _hdr['BUNIT'] = (bunit, 'Unit of pixel value')

        # Add the common HDUCLASS keys
        _hdr['HDUCLASS'] = ('SDSS', 'SDSS format class')
        _hdr['HDUCLAS1'] = ('CUBE', 'Data format')
        if hduclas2 == 'DATA':
            _hdr['HDUCLAS2'] = 'DATA'
            if err:
                _hdr['ERRDATA'] = (ext+'_IVAR' if prepend else 'IVAR',
                                    'Associated inv. variance extension')
            if qual:
                _hdr['QUALDATA'] = (ext+'_MASK' if prepend else 'MASK',
                                    'Associated quality extension')
            return _hdr

        if hduclas2 == 'ERROR':
            _hdr['HDUCLAS2'] = 'ERROR'
            _hdr['HDUCLAS3'] = ('INVMSE', 'Value is inverse mean-square error')
            _hdr['SCIDATA'] = (ext, 'Associated data extension')
            if qual:
                _hdr['QUALDATA'] = (ext+'_MASK' if prepend else 'MASK',
                                    'Associated quality extension')
            return _hdr

        if hduclas2 == 'QUALITY':
            _hdr['HDUCLAS2'] = 'QUALITY'
            _hdr['HDUCLAS3'] = construct_maps_file._mask_data_type(bit_type)
            _hdr['SCIDATA'] = (ext, 'Associated data extension')
            if err:
                _hdr['ERRDATA'] = (ext+'_IVAR' if prepend else 'IVAR',
                                    'Associated inv. variance extension')
            return _hdr
            
        raise ValueError('HDUCLAS2 must be DATA, ERROR, or QUALITY.')


    def _initialize_mask(self, binned_spectra):
        """
        Initialize the mask based on the DRP cube mask.

        IVAR_INVALID:  Test of invalid inverse variance based on the DRP
        file.  THIS SHOULD BE A TEST THAT IS INCLUDED IN THE
        BINNED_SPECTRA OBJECT.  This should match the flags from the
        stellar continuum model object.

        Copy the FORESTAR flags from the DRP file to the this one
        """
        mask = numpy.zeros(binned_spectra.shape, dtype=self.bitmask.minimum_dtype())
        indx = ~(binned_spectra['IVAR'].data > 0) | ~(numpy.isfinite(binned_spectra['IVAR'].data))
        mask[indx] = self.bitmask.turn_on(mask[indx], 'IVARINVALID')

        indx = binned_spectra.bitmask.flagged(binned_spectra['MASK'].data, flag='FORESTAR')
        mask[indx] = self.bitmask.turn_on(mask[indx], 'FORESTAR')

        return mask


#    def _include_stellar_continuum_mask(self, mask, stellar_continuum):
#        """
#        From the stellar continuum model object:
#            - copy ARTIFACT
#            - copy INVALID_ERROR into IVAR_INVALID
#            - consolidate DIDNOTUSE, LOW_SNR, OUTSIDE_RANGE, EML_REGION,
#              TPL_PIXELS, TRUNCATED, PPXF_REJECT, INVALID_ERROR into SC_IGNORED
#            - copy FIT_FAILED, NEAR_BOUND into SC_FAILED
#        """
#        indx = self.bitmask.flagged(self.drpf['MASK'].data, flag='ARTIFACT')
#        mask[indx] = self.bitmask.turn_on(mask[indx], 'ARTIFACT')
#
#        indx = self.bitmask.flagged(self.drpf['MASK'].data, flag='INVALID_ERROR')
#        mask[indx] = self.bitmask.turn_on(mask[indx], 'IVAR_INVALID')
#
#        flags = [ 'DIDNOTUSE', 'LOW_SNR', 'OUTSIDE_RANGE', 'EML_REGION', 'TPL_PIXELS', 'TRUNCATED',
#                  'PPXF_REJECT', 'INVALID_ERROR' ]
#        indx = self.bitmask.flagged(self.drpf['MASK'].data, flag=flags)
#        mask[indx] = self.bitmask.turn_on(mask[indx], 'SC_IGNORED')
#
#        indx = self.bitmask.flagged(self.drpf['MASK'].data, flag='FIT_FAILED')
#        mask[indx] = self.bitmask.turn_on(mask[indx], 'SC_FAILED')
#
#        return mask


    def _include_model_and_data_mask(self, mask, binned_spectra, stellar_continuum,
                                     emission_line_model):
        """
        For the binned spectra:
            - consolidate DIDNOTUSE, LOW_SPECCOV, LOW_SNR from
              binned_spectra into IGNORED
            - consolidate NONE_IN_STACK from binned_spectra into
              FLUXINVALID

        For the model spectra:
            - copy INVALID_ERROR from stellar_continuum into
              IVARINVALID

            - copy ARTIFACTs from both stellar_continuum and
              emission_line_model

            - from stellar_continuum, consolidate DIDNOTUSE, LOW_SNR,
              OUTSIDE_RANGE, EML_REGION, TPL_PIXELS, TRUNCATED,
              PPXF_REJECT, INVALID_ERROR into a list of pixels ignored
              by the stellar-continuum fit
            - from emission_line_model, consolidate DIDNOTUSE, LOW_SNR,
              OUTSIDE_RANGE into a list of pixels ignored by the
              emission-line fit
            - flag pixels as FITIGNORED if the pixel is ignored by
              **both** the stellar-continuum and emission-line fits

            - from stellar_continuum, consolidate FIT_FAILED and
              NEAR_BOUND into a list of failed pixels from the
              stellar-continuum fit
            - do the same for the emission-line fitting mask
            - flag pixels as FIT_FAILED if the pixel failed in
              **either** the stellar-continuum or emission-line fits
        """
        # Draw from binned_spectra
        flags = [ 'DIDNOTUSE', 'LOW_SPECCOV', 'LOW_SNR' ]
        indx = binned_spectra.bitmask.flagged(binned_spectra['MASK'].data, flag=flags)
        mask[indx] = self.bitmask.turn_on(mask[indx], 'IGNORED')

        indx = binned_spectra.bitmask.flagged(binned_spectra['MASK'].data, flag='NONE_IN_STACK')
        mask[indx] = self.bitmask.turn_on(mask[indx], 'FLUXINVALID')

        # Draw from stellar_continuum
        indx = stellar_continuum.bitmask.flagged(stellar_continuum['MASK'].data,
                                                 flag='INVALID_ERROR')
        mask[indx] = self.bitmask.turn_on(mask[indx], 'IVARINVALID')

        # Draw from stellar_continuum and emission_line_model
        indx = stellar_continuum.bitmask.flagged(stellar_continuum['MASK'].data, flag='ARTIFACT') \
                    | emission_line_model.bitmask.flagged(emission_line_model['MASK'].data,
                                                          flag='ARTIFACT')
        mask[indx] = self.bitmask.turn_on(mask[indx], 'ARTIFACT')

        flags = [ 'DIDNOTUSE', 'LOW_SNR', 'OUTSIDE_RANGE', 'EML_REGION', 'TPL_PIXELS',
                  'TRUNCATED', 'PPXF_REJECT', 'INVALID_ERROR' ]
        sc_indx = stellar_continuum.bitmask.flagged(stellar_continuum['MASK'].data, flag=flags)
        flags = [ 'DIDNOTUSE', 'LOW_SNR', 'OUTSIDE_RANGE' ]
        el_indx = emission_line_model.bitmask.flagged(emission_line_model['MASK'].data, flag=flags)
        mask[sc_indx & el_indx] = self.bitmask.turn_on(mask[sc_indx & el_indx], 'FITIGNORED')

#        t_sc = numpy.ma.MaskedArray(stellar_continuum['FLUX'].data, mask=~sc_indx)
#        t_el = numpy.ma.MaskedArray(stellar_continuum['FLUX'].data, mask=~el_indx)
#        t_bt = numpy.ma.MaskedArray(stellar_continuum['FLUX'].data, mask=(el_indx & sc_indx))
#        pyplot.step(stellar_continuum['WAVE'].data, t_sc.data[22,22,:], where='mid', lw=4.0, color='0.5', zorder=1)
#        pyplot.step(stellar_continuum['WAVE'].data, t_sc[22,22,:], where='mid', lw=2.0, color='r', zorder=2)
#        pyplot.step(stellar_continuum['WAVE'].data, t_el[22,22,:], where='mid', lw=1.0, color='b', zorder=3)
#        pyplot.step(stellar_continuum['WAVE'].data, t_bt.data[22,22,:], where='mid', lw=4.0, color='0.5', zorder=1)
#        pyplot.step(stellar_continuum['WAVE'].data, t_bt[22,22,:], where='mid', lw=2.0, color='r', zorder=2)
#        pyplot.show()

        # TODO: What do I do with BAD_SIGMA?
        sc_indx = stellar_continuum.bitmask.flagged(stellar_continuum['MASK'].data,
                                                    flag=['FIT_FAILED', 'NEAR_BOUND'])
        el_indx = emission_line_model.bitmask.flagged(emission_line_model['MASK'].data,
                                                      flag=['FIT_FAILED', 'NEAR_BOUND'])
        mask[sc_indx | el_indx] = self.bitmask.turn_on(mask[sc_indx | el_indx], 'FITFAILED')

        return mask


    def model_and_data_cube(self, prihdr, binned_spectra, stellar_continuum, emission_line_model):
        # Construct and return the empty hdus
        if binned_spectra is None or stellar_continuum is None or emission_line_model is None:
            hdus = [ fits.ImageHDU(data=None,
                                   header=self._finalize_cube_header(self.base_cubehdr, 'FLUX',
                                                                bunit='1E-17 erg/s/cm^2/ang/spaxel',
                                                                     err=True, qual=True),
                                   name='FLUX') ]
            hdus += [ fits.ImageHDU(data=None,
                                    header=self._finalize_cube_header(self.base_cubehdr, 'FLUX',
                                                                      hduclas2='ERROR',
                                                        bunit='(1E-17 erg/s/cm^2/ang/spaxel)^{-2}',
                                                                      qual=True, prepend=False),
                                    name='IVAR') ]
            hdus += [ fits.ImageHDU(data=None,
                                    header=self._finalize_cube_header(self.base_cubehdr, 'FLUX',
                                                                      hduclas2='QUALITY', err=True,
                                                                      bit_type=numpy.bool,
                                                                      prepend=False),
                                    name='MASK') ]
            hdus += [ fits.ImageHDU(data=None, name='WAVE') ]
            hdus += [ fits.ImageHDU(data=None, name='REDCORR') ]
            hdus += [ fits.ImageHDU(data=None,
                                    header=self._finalize_cube_header(self.base_cubehdr, 'MODEL',
                                                               bunit='1E-17 erg/s/cm^2/ang/spaxel',
                                                                      qual=True, prepend=False),
                                    name='MODEL') ]
            return prihdr, hdus, []

        # Add data to the primary header
        prihdr = binned_spectra.rdxqa._initialize_header(prihdr)
        prihdr = binned_spectra._initialize_header(prihdr)
        prihdr = binned_spectra._add_method_header(prihdr)
        prihdr = binned_spectra._get_reddening_hdr(prihdr)
        prihdr = stellar_continuum._initialize_header(prihdr)
        prihdr = stellar_continuum._add_method_header(prihdr)
        prihdr = emission_line_model._initialize_header(prihdr)
        prihdr = emission_line_model._add_method_header(prihdr)

        # Get the reddened data
        data, ivar = apply_reddening(binned_spectra['FLUX'].data, binned_spectra['REDCORR'].data,
                                     deredden=False, ivar=binned_spectra['IVAR'].data)

        # Get the model
        model = stellar_continuum['FLUX'].data + emission_line_model['FLUX'].data \
                    + emission_line_model['BASE'].data
        # and redden it
        model = apply_reddening(model, binned_spectra['REDCORR'].data, deredden=False)

        # Get the mask data
        mask = self.common_mask.copy()
        mask = self._include_model_and_data_mask(mask, binned_spectra, stellar_continuum,
                                                 emission_line_model)

        # Construct the list of HDUs
        hdus = [ fits.ImageHDU(data=data,
                               header=self._finalize_cube_header(self.base_cubehdr, 'FLUX',
                                                            bunit='1E-17 erg/s/cm^2/ang/spaxel',
                                                                 err=True, qual=True),
                               name='FLUX') ]
        hdus += [ fits.ImageHDU(data=ivar,
                                header=self._finalize_cube_header(self.base_cubehdr, 'FLUX',
                                                                  hduclas2='ERROR',
                                                        bunit='(1E-17 erg/s/cm^2/ang/spaxel)^{-2}',
                                                                  qual=True, prepend=False),
                                name='IVAR') ]
        hdus += [ fits.ImageHDU(data=mask,
                                header=self._finalize_cube_header(self.base_cubehdr, 'FLUX',
                                                                  hduclas2='QUALITY', err=True,
                                                                  bit_type=mask.dtype.type,
                                                                  prepend=False),
                                name='MASK') ]
        hdus += [ fits.ImageHDU(data=binned_spectra['WAVE'].data, name='WAVE') ]
        hdus += [ fits.ImageHDU(data=binned_spectra['REDCORR'].data, name='REDCORR') ]
        hdus += [ fits.ImageHDU(data=model,
                                header=self._finalize_cube_header(self.base_cubehdr, 'MODEL',
                                                             bunit='1E-17 erg/s/cm^2/ang/spaxel',
                                                                  qual=True, prepend=False),
                                name='MODEL') ]
        return prihdr, hdus, ['FLUX', 'IVAR', 'MASK', 'MODEL']


    def _include_emission_line_model_mask(self, mask, emission_line_model):
        """
        From the emission-line model object:
            - copy ARTIFACT
            - consolidate DIDNOTUSE, LOW_SNR, OUTSIDE_RANGE into IGNORED_EL
            - copy FIT_FAILED, NEAR_BOUND into ELFAILED
        """
        indx = emission_line_model.bitmask.flagged(emission_line_model['MASK'].data,
                                                   flag='ARTIFACT')
        mask[indx] = self.bitmask.turn_on(mask[indx], 'ARTIFACT')

        indx = emission_line_model.bitmask.flagged(emission_line_model['MASK'].data,
                                                   flag=['DIDNOTUSE', 'LOW_SNR', 'OUTSIDE_RANGE'])
        mask[indx] = self.bitmask.turn_on(mask[indx], 'ELIGNORED')

        indx = emission_line_model.bitmask.flagged(emission_line_model['MASK'].data,
                                                   flag=['FIT_FAILED', 'NEAR_BOUND'])
        mask[indx] = self.bitmask.turn_on(mask[indx], 'ELFAILED')

        return mask


    def emission_line_model_cube(self, prihdr, binned_spectra, emission_line_model):
        # Construct and return the empty hdus
        if emission_line_model is None:
            hdus = [ fits.ImageHDU(data=None,
                                   header=self._finalize_cube_header(self.base_cubehdr, 'EMLINE',
                                                               bunit='1E-17 erg/s/cm^2/ang/spaxel',
                                                                     qual=True),
                                   name='EMLINE') ]
            hdus += [ fits.ImageHDU(data=None,
                                    header=self._finalize_cube_header(self.base_cubehdr, 'EMLINE',
                                                            bunit='1E-17 erg/s/cm^2/ang/spaxel',
                                                                     qual=True),
                                    name='EMLINE_BASE') ]
            hdus += [ fits.ImageHDU(data=None,
                                    header=self._finalize_cube_header(self.base_cubehdr, 'EMLINE',
                                                                      hduclas2='QUALITY',
                                                                      bit_type=numpy.bool),
                                    name='EMLINE_MASK') ]
            return prihdr, hdus, []

        # Add data to the primary header
        prihdr = emission_line_model._initialize_header(prihdr)
        prihdr = emission_line_model._add_method_header(prihdr)

        # Get the reddened line model and baselines
        model = apply_reddening(emission_line_model['FLUX'].data, binned_spectra['REDCORR'].data,
                                deredden=False)
        base = apply_reddening(emission_line_model['BASE'].data, binned_spectra['REDCORR'].data,
                               deredden=False)

        # Get the mask data
        mask = self.common_mask.copy()
        mask = self._include_emission_line_model_mask(mask, emission_line_model)

        # Construct the list of HDUs
        hdus = [ fits.ImageHDU(data=model,
                               header=self._finalize_cube_header(self.base_cubehdr, 'EMLINE',
                                                               bunit='1E-17 erg/s/cm^2/ang/spaxel',
                                                                     qual=True),
                               name='EMLINE') ]
        hdus += [ fits.ImageHDU(data=base,
                                header=self._finalize_cube_header(self.base_cubehdr, 'EMLINE',
                                                            bunit='1E-17 erg/s/cm^2/ang/spaxel',
                                                                  qual=True),
                                name='EMLINE_BASE') ]
        hdus += [ fits.ImageHDU(data=mask,
                                header=self._finalize_cube_header(self.base_cubehdr, 'EMLINE',
                                                                  hduclas2='QUALITY',
                                                                  bit_type=mask.dtype.type),
                                name='EMLINE_MASK') ]
        return prihdr, hdus, ['EMLINE', 'EMLINE_BASE', 'EMLINE_MASK']
        
