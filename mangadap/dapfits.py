# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Defines a class used to interface with the final maps files produced by
the MaNGA Data Analysis Pipeline (DAP).

.. todo::
    - Allow DAPFits to read/access both the MAPS and the LOGCUBE files,
      not just the former.

----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
from pathlib import Path
import logging
import warnings

from IPython import embed

import numpy

from astropy.io import fits
import astropy.units
from astropy.cosmology import FlatLambdaCDM

from .util.constants import DAPConstants
from .util.fitsutil import DAPFitsUtil
from .util.dapbitmask import DAPBitMask
from .util.log import log_output
from .util.covariance import Covariance
from .config import defaults

# NOTE: These are only required for type checking.
from .datacube import DataCube
from .proc.reductionassessments import ReductionAssessment
from .proc.spatiallybinnedspectra import SpatiallyBinnedSpectra
from .proc.stellarcontinuummodel import StellarContinuumModel
from .proc.emissionlinemoments import EmissionLineMoments
from .proc.emissionlinemodel import EmissionLineModel
from .proc.spectralindices import SpectralIndices


#-----------------------------------------------------------------------
class DAPQualityBitMask(DAPBitMask):
    """
    Class used to flag/interpret global quality bit.

    The defined bits are listed at :ref:`metadatamodel-dapqual`.

    .. todo::
        - Force read IDLUTILS version as opposed to internal one?
    """
    cfg_root = 'dap_quality_bits'


class DAPMapsBitMask(DAPBitMask):
    """
    Class used to flag/interpret bits in DAP MAPS files.

    The defined bits are listed at :ref:`metadatamodel-dappixmask`.

    .. todo::
        - Force read IDLUTILS version as opposed to internal one?
    """
    cfg_root = 'dap_maps_bits'


class DAPCubeBitMask(DAPBitMask):
    """
    Class used to flag/interpret bits in DAP model LOGCUBE files.

    The defined bits are listed at :ref:`metadatamodel-dapspecmask`.

    .. todo::
        - Force read IDLUTILS version as opposed to internal one?
    """
    cfg_root = 'dap_cube_bits'


#-----------------------------------------------------------------------
class construct_maps_file:
    """
    Construct a DAP MAPS file.

    Set as a class to enforce a namespace. Should not be used as an
    object!

    Should force all intermediate objects to be provided.

    """
    def __init__(self, cube, method=None, rdxqa=None, binned_spectra=None, stellar_continuum=None,
                 emission_line_moments=None, emission_line_model=None, spectral_indices=None,
                 redshift=None, output_path=None, output_file=None, overwrite=True, loggers=None,
                 quiet=False, single_precision=False):

        #---------------------------------------------------------------
        # Initialize the reporting
        self.loggers = None if loggers is None else loggers
        self.quiet = quiet

        self.float_dtype = 'float32' if single_precision else 'float'

        #---------------------------------------------------------------
        # Check input types
        # TODO: Do we need this?
        confirm_dap_types(cube, rdxqa, binned_spectra, stellar_continuum, emission_line_moments,
                          emission_line_model, spectral_indices)

        #---------------------------------------------------------------
        # Set the output paths
        self.cube = cube
        self.meta = cube.meta
        self.method, self.directory_path, self.output_file \
                = default_paths(cube, file_type='MAPS',
                        rdxqa_method=None if rdxqa is None else rdxqa.method['key'],
                        binning_method=None if binned_spectra is None 
                                       else binned_spectra.method['key'],
                        stelcont_method=None if stellar_continuum is None
                                        else stellar_continuum.method['key'],
                        elmom_method=None if emission_line_moments is None
                                     else emission_line_moments.database['key'],
                        elmodel_method=None if emission_line_model is None
                                       else emission_line_model.method['key'],
                        sindex_method=None if spectral_indices is None
                                      else spectral_indices.database['key'],
                        method=method, output_path=output_path, output_file=output_file)

#        self._set_paths(directory_path, dapver, analysis_path, output_file, binned_spectra,
#                        stellar_continuum, emission_line_model)

        # Save input for reference
        self.spatial_shape = self.cube.spatial_shape
        self.redshift = redshift

        # Report
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, '-'*50)
            log_output(loggers, 1, logging.INFO, f'{"CONSTRUCTING OUTPUT MAPS":^50}')
            log_output(self.loggers, 1, logging.INFO, '-'*50)
            log_output(self.loggers, 1, logging.INFO, f'Output path: {self.directory_path}')
            log_output(self.loggers, 1, logging.INFO, f'Output file: {self.output_file}')
            log_output(self.loggers, 1, logging.INFO, f'Output map shape: {self.spatial_shape}')
            log_output(self.loggers, 1, logging.INFO, f'Redshift: {self.redshift:.4f}')

        #---------------------------------------------------------------
        # Check if the file already exists
        ofile = self.directory_path / self.output_file
        if ofile.exists() and not overwrite:
            # TODO: Perform some checks to make sure the existing file
            # has the correct content?
            warnings.warn('Output file exists!  Set overwrite=True to overwrite.')
            return

        #---------------------------------------------------------------
        # Initialize the primary header
        prihdr = DAPFitsUtil.initialize_dap_primary_header(self.cube, maskname='MANGA_DAPPIXMASK')

        # Add the DAP method
        prihdr['DAPTYPE'] = (self.method, 'DAP analysis method')
        # Add the format of this file
        prihdr['DAPFRMT'] = ('MAPS', 'DAP data file format')

        # Get the base map headers
        self.multichannel_maphdr \
                = DAPFitsUtil.build_map_header(self.cube.fluxhdr,
                                'K Westfall <westfall@ucolick.org> & SDSS-IV Data Group',
                                               multichannel=True, maskname='MANGA_DAPPIXMASK')

        self.singlechannel_maphdr \
                = DAPFitsUtil.build_map_header(self.cube.fluxhdr,
                                'K Westfall <westfall@ucolick.org> & SDSS-IV Data Group',
                                               maskname='MANGA_DAPPIXMASK')

        #---------------------------------------------------------------
        # Initialize the pixel mask
        self.bitmask = DAPMapsBitMask()

        # Get the mask for the binned spectra
        self.bin_mask = self._build_binning_mask(binned_spectra)

        #---------------------------------------------------------------
        # Construct the hdu list for each input object.
        # Reduction assessments:
        rdxqalist = self.reduction_assessment_maps(prihdr, rdxqa)
        # Construct the BINID extension
        binidlist = combine_binid_extensions(self.cube, binned_spectra, stellar_continuum,
                                             emission_line_moments, emission_line_model,
                                             spectral_indices, dtype='int32')
        # Binned spectra:
        bspeclist = self.binned_spectra_maps(prihdr, binned_spectra)
        # Stellar-continuum fits:
        contlist = self.stellar_continuum_maps(prihdr, stellar_continuum)
        # Emission-line moments:
        emlmomlist = self.emission_line_moment_maps(prihdr, emission_line_moments)
        # Emission-line models:
        elmodlist = self.emission_line_model_maps(prihdr, emission_line_model)

        # Spectral indices:
        sindxlist = self.spectral_index_maps(prihdr, spectral_indices)

        # Save the data to the hdu attribute
        prihdr = add_snr_metrics_to_header(prihdr, self.cube, rdxqalist[1].data[:,:,1].ravel())
        
        prihdr = finalize_dap_primary_header(prihdr, self.cube, self.meta, binned_spectra,
                                             stellar_continuum, loggers=self.loggers,
                                             quiet=self.quiet)
        self.hdu = fits.HDUList([ fits.PrimaryHDU(header=prihdr),
                                  *rdxqalist,
                                  *binidlist,
                                  *bspeclist,
                                  *contlist,
                                  *emlmomlist,
                                  *elmodlist,
                                  *sindxlist
                                ])

        extensions = [ h.name for h in self.hdu ]

        #---------------------------------------------------------------
        # TEMPORARY FLAGS:
        # Flag the Gaussian-fitted flux as unreliable if the summed flux
        # is not within the factor below
#        if emission_line_moments is not None and emission_line_model is not None \
#                    and self.hdu['EMLINE_GFLUX'].data.shape != self.hdu['EMLINE_SFLUX'].data.shape:
#            warnings.warn('Cannot compare emission-line moment and model fluxes!')
#        elif emission_line_moments is not None and emission_line_model is not None \
#                    and self.hdu['EMLINE_GFLUX'].data.shape == self.hdu['EMLINE_SFLUX'].data.shape:
#            factor = 5.0
#            indx = (self.hdu['EMLINE_GFLUX_MASK'].data == 0) \
#                    & (self.hdu['EMLINE_SFLUX_MASK'].data == 0) \
#                    & ( (self.hdu['EMLINE_GFLUX'].data < self.hdu['EMLINE_SFLUX'].data/factor)
#                        | (self.hdu['EMLINE_GFLUX'].data > self.hdu['EMLINE_SFLUX'].data*factor) ) \
#                    & ( (emission_line_moments['BINID'].data > -1)
#                        & (emission_line_model['BINID'].data > -1) )[:,:,None]
#            print('unreliable Gaussian flux compared to summed flux: ', numpy.sum(indx))
#            if numpy.sum(indx) > 0:
#                self.hdu['EMLINE_GFLUX_MASK'].data[indx] \
#                    = self.bitmask.turn_on(self.hdu['EMLINE_GFLUX_MASK'].data[indx], 'UNRELIABLE')
#                self.hdu['EMLINE_GVEL_MASK'].data[indx] \
#                    = self.bitmask.turn_on(self.hdu['EMLINE_GVEL_MASK'].data[indx], 'UNRELIABLE')
#                self.hdu['EMLINE_GSIGMA_MASK'].data[indx] \
#                    = self.bitmask.turn_on(self.hdu['EMLINE_GSIGMA_MASK'].data[indx], 'UNRELIABLE')

            # TODO: Add if EMLINE_GSIGMA < EMLINE_INSTSIGMA !!

        # Flag the stellar velocity dispersions measured in spectra with
        # S/N<10 as unreliable
#        indx = (binned_spectra['BINID'].data > -1) & (self.hdu['BIN_SNR'].data < 10) \
#                    & (stellar_continuum['BINID'].data > -1) \
#                    & (self.hdu['STELLAR_SIGMA_MASK'].data == 0)
#        print('unreliable sigma because of low S/N: ', numpy.sum(indx))
#        if numpy.sum(indx):
#            self.hdu['STELLAR_SIGMA_MASK'].data[indx] \
#                    = self.bitmask.turn_on(self.hdu['STELLAR_SIGMA_MASK'].data[indx], 'UNRELIABLE')

        # TODO: FIX THIS !!!!!!!!!!!!!
        # Flag any inverse variances that are not positive as DONOTUSE
        # and MATHERROR
        ext = ['BIN_MFLUX', 'STELLAR_VEL', 'STELLAR_SIGMA', 'EMLINE_SFLUX', 'EMLINE_SEW',
               'EMLINE_GFLUX', 'EMLINE_GEW', 'EMLINE_GVEL', 'EMLINE_GSIGMA', 'SPECINDEX']
        for e in ext:
            if f'{e}_MASK' not in extensions \
                    or f'{e}_IVAR' not in extensions \
                    or self.hdu[f'{e}_MASK'].data is None \
                    or self.hdu[f'{e}_IVAR'].data is None:
                continue
            indx = numpy.logical_not(self.bitmask.flagged(self.hdu['{0}_MASK'.format(e)].data)) \
                            & numpy.logical_not(self.hdu['{0}_IVAR'.format(e)].data > 0)
            if numpy.sum(indx) > 0:
                self.hdu['{0}_MASK'.format(e)].data[indx] \
                    = self.bitmask.turn_on(self.hdu['{0}_MASK'.format(e)].data[indx], 'MATHERROR')
                self.hdu['{0}_MASK'.format(e)].data[indx] \
                    = self.bitmask.turn_on(self.hdu['{0}_MASK'.format(e)].data[indx], 'DONOTUSE')

        #---------------------------------------------------------------
        # Check that the path exists
        if not self.directory_path.exists():
            self.directory_path.mkdir(parents=True)
        # Write the maps file
        DAPFitsUtil.write(self.hdu, str(ofile), overwrite=overwrite, checksum=True,
                          loggers=self.loggers, quiet=self.quiet)
        # End
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, '-'*50)

    def _consolidate_donotuse(self, mask):
        return self.bitmask.consolidate(mask, ['NOCOV', 'LOWCOV', 'DEADFIBER', 'FORESTAR',
                                               'NOVALUE', 'MATHERROR', 'FITFAILED', 'NEARBOUND'],
                                        'DONOTUSE')

    def _build_binning_mask(self, binned_spectra):
        """
        Consolidate the map mask values into NOVALUE.
        """
        # Consolidate everything in the map mask into NOVALUE
        indx = binned_spectra.bitmask.flagged(binned_spectra['MAPMASK'].data)
        mask = numpy.zeros(indx.shape, dtype=self.bitmask.minimum_dtype())
        mask[indx] = self.bitmask.turn_on(mask[indx], 'NOVALUE')

        # Isolate FORESTAR
        indx = binned_spectra.bitmask.flagged(binned_spectra['MAPMASK'].data, flag='FORESTAR')
        mask[indx] = self.bitmask.turn_on(mask[indx], 'FORESTAR')

        # Marginalize the binned spectra mask just for NONE_IN_STACK and
        # convert it to NOVALUE
        bin_mask = DAPFitsUtil.marginalize_mask(binned_spectra['MASK'].data,
                                                inp_flags='NONE_IN_STACK',
                                                inp_bitmask=binned_spectra.bitmask, 
                                                out_flag='NOVALUE', out_bitmask=self.bitmask,
                                                dispaxis=1)

        # Reconstruct the full mask with just this flag
        bin_mask = DAPFitsUtil.reconstruct_map(self.spatial_shape,
                                               binned_spectra['BINID'].data.ravel(), bin_mask)

        # Add these bits to the full mask
        indx = self.bitmask.flagged(bin_mask, 'NOVALUE')
        mask[indx] = self.bitmask.turn_on(mask[indx], 'NOVALUE')
        
        # Return the mask after consolidating flags into DONOTUSE
        return self._consolidate_donotuse(mask)


    def _stellar_continuum_mask_to_map_mask(self, stellar_continuum, sc_mask, for_dispersion=False):
        """

        sc_mask constains the reconstructed map of the masks in each
        stellar continuum fit (using the stellar continuum bitmask).

        Propagate and consolidate the stellar-continuum masks to the map
        pixel mask.

        DIDNOTUSE, FORESTAR propagated from already existing mask
            (self.bin_mask)

        LOW_SNR, NO_FIT, INSUFFICIENT_DATA propagated to NOVALUE

        FIT_FAILED, NEAR_BOUND propagated to FITFAILED and NEARBOUND

        NEGATIVE_WEIGHTS consolidated into UNRELIABLE; if
            dispersion=True, include BAD_SIGMA in this
        """
        # Consolidate everything in the map mask into NOVALUE
        indx = stellar_continuum.bitmask.flagged(stellar_continuum['MAPMASK'].data)
        mask = numpy.zeros(indx.shape, dtype=self.bitmask.minimum_dtype())
        mask[indx] = self.bitmask.turn_on(mask[indx], 'NOVALUE')

        # Isolate FORESTAR
        indx = stellar_continuum.bitmask.flagged(stellar_continuum['MAPMASK'].data,
                                                 flag='FORESTAR')
        mask[indx] = self.bitmask.turn_on(mask[indx], 'FORESTAR')
        
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

        # Convert MIN_SIGMA into a NEARBOUND for the dispersion
        if for_dispersion:
            flgd = stellar_continuum.bitmask.flagged(sc_mask, flag='BAD_SIGMA')
            mask[flgd] = self.bitmask.turn_on(mask[flgd], 'UNRELIABLE')

            flgd = stellar_continuum.bitmask.flagged(sc_mask, flag='MIN_SIGMA')
            mask[flgd] = self.bitmask.turn_on(mask[flgd], 'NEARBOUND')

        return self._consolidate_donotuse(mask)


    def _emission_line_moment_mask_to_map_mask(self, emission_line_moments, elm_mask):
        """
        Propagate and consolidate the emission-line moment masks to the map
        pixel mask.

        DIDNOTUSE, FORESTAR propagated from already existing mask (self.bin_mask)

        MAIN_EMPTY, BLUE_EMPTY, RED_EMPTY, UNDEFINED_BANDS propagated to NOVALUE

        MAIN_JUMP, BLUE_JUMP, RED_JUMP, JUMP_BTWN_SIDEBANDS propagated to FITFAILED

        MAIN_INCOMP, BLUE_INCOMP, RED_INCOMP propagated to UNRELIABLE

        NO_ABSORPTION_CORRECTION propagated to NOCORRECTION

        DIVBYZERO propagated to MATHERROR

        Second moments not provided so no need to include UNDEFINED_MOM2
        """
        # Consolidate everything in the map mask into NOVALUE
        indx = emission_line_moments.bitmask.flagged(emission_line_moments['MAPMASK'].data)
        mask = numpy.zeros(indx.shape, dtype=self.bitmask.minimum_dtype())
        mask[indx] = self.bitmask.turn_on(mask[indx], 'NOVALUE')

        # Isolate FORESTAR
        indx = emission_line_moments.bitmask.flagged(emission_line_moments['MAPMASK'].data,
                                                     flag='FORESTAR')
        mask[indx] = self.bitmask.turn_on(mask[indx], 'FORESTAR')

        # Reshape the mask to include all emission lines
        mask = numpy.array([mask]*emission_line_moments.nmom).transpose(1,2,0)

        # Consolidate to NOVALUE
        flgd = emission_line_moments.bitmask.flagged(elm_mask,
                                                     flag=['MAIN_EMPTY', 'BLUE_EMPTY', 'RED_EMPTY',
                                                           'UNDEFINED_BANDS' ])
        mask[flgd] = self.bitmask.turn_on(mask[flgd], 'NOVALUE')

        # Consolidate to FITFAILED
        flgd = emission_line_moments.bitmask.flagged(elm_mask, flag=['MAIN_JUMP', 'BLUE_JUMP',
                                                                 'RED_JUMP', 'JUMP_BTWN_SIDEBANDS'])
        mask[flgd] = self.bitmask.turn_on(mask[flgd], 'FITFAILED')

        # Consolidate to UNRELIABLE
        flgd = emission_line_moments.bitmask.flagged(elm_mask, flag=['MAIN_INCOMP', 'BLUE_INCOMP',
                                                                     'RED_INCOMP' ])
        mask[flgd] = self.bitmask.turn_on(mask[flgd], 'UNRELIABLE')

        # Consolidate to NOCORRECTION
        flgd = emission_line_moments.bitmask.flagged(elm_mask, flag=['NO_ABSORPTION_CORRECTION'])
        mask[flgd] = self.bitmask.turn_on(mask[flgd], 'NOCORRECTION')

        # Consolidate to MATHERROR
        flgd = emission_line_moments.bitmask.flagged(elm_mask, flag=['DIVBYZERO'])
        mask[flgd] = self.bitmask.turn_on(mask[flgd], 'MATHERROR')

        return self._consolidate_donotuse(mask)


    def _emission_line_model_mask_to_map_mask(self, emission_line_model, elf_mask,
                                              for_dispersion=False):
        """
        Propagate and consolidate the emission-line moment masks to the map
        pixel mask.

        INSUFFICIENT_DATA propagated to NOVALUE

        FIT_FAILED, UNDEFINED_COVAR propagated to FITFAILED

        NEAR_BOUND copied to NEARBOUND

        UNDEFINED_SIGMA propagated to UNRELIABLE
        """
        # Consolidate everything in the map mask into NOVALUE
        indx = emission_line_model.bitmask.flagged(emission_line_model['MAPMASK'].data)
        mask = numpy.zeros(indx.shape, dtype=self.bitmask.minimum_dtype())
        mask[indx] = self.bitmask.turn_on(mask[indx], 'NOVALUE')

        # Isolate FORESTAR
        indx = emission_line_model.bitmask.flagged(emission_line_model['MAPMASK'].data,
                                                   flag='FORESTAR')
        mask[indx] = self.bitmask.turn_on(mask[indx], 'FORESTAR')

        # Reshape the mask to include all emission lines
        mask = numpy.array([mask]*emission_line_model.neml).transpose(1,2,0)

        # Consolidate to NOVALUE
        flgd = emission_line_model.bitmask.flagged(elf_mask, flag=['INSUFFICIENT_DATA'])
        mask[flgd] = self.bitmask.turn_on(mask[flgd], 'NOVALUE')

        # Consolidate to FITFAILED
        flgd = emission_line_model.bitmask.flagged(elf_mask, flag=['FIT_FAILED', 'UNDEFINED_COVAR'])
        mask[flgd] = self.bitmask.turn_on(mask[flgd], 'FITFAILED')

        # Copy NEAR_BOUND
        flgd = emission_line_model.bitmask.flagged(elf_mask, flag='NEAR_BOUND')
        mask[flgd] = self.bitmask.turn_on(mask[flgd], 'NEARBOUND')

        # Consolidate dispersion flags to NEARBOUND and UNRELIABLE
        if for_dispersion:
            flgd = emission_line_model.bitmask.flagged(elf_mask, flag='MIN_SIGMA')
            mask[flgd] = self.bitmask.turn_on(mask[flgd], 'NEARBOUND')

            flgd = emission_line_model.bitmask.flagged(elf_mask, flag=['UNDEFINED_SIGMA',
                                                                       'BAD_SIGMA'])
            mask[flgd] = self.bitmask.turn_on(mask[flgd], 'UNRELIABLE')

        return self._consolidate_donotuse(mask)


    def _spectral_index_mask_to_map_mask(self, spectral_indices, si_mask):
        """
        Propagate and consolidate the spectral index masks to the map
        pixel mask.

        DIDNOTUSE, FORESTAR propagated from already existing mask (self.bin_mask)

        MAIN_EMPTY, BLUE_EMPTY, RED_EMPTY, UNDEFINED_BANDS propagated to NOVALUE

        MAIN_INCOMP, BLUE_INCOMP, RED_INCOMP propagated to UNRELIABLE

        NO_DISPERSION_CORRECTION propagated to NOCORRECTION

        DIVBYZERO propagated to MATHERROR

        Need to further assess all the INCOMP flags to see if these lead to bad
        values (BADVALUE).
        """
        # Consolidate everything in the map mask into NOVALUE
        indx = spectral_indices.bitmask.flagged(spectral_indices['MAPMASK'].data)
        mask = numpy.zeros(indx.shape, dtype=self.bitmask.minimum_dtype())
        mask[indx] = self.bitmask.turn_on(mask[indx], 'NOVALUE')

        # Isolate FORESTAR
        indx = spectral_indices.bitmask.flagged(spectral_indices['MAPMASK'].data, flag='FORESTAR')
        mask[indx] = self.bitmask.turn_on(mask[indx], 'FORESTAR')

        # Reshape the mask to include all the spectral indices
        mask = numpy.array([mask]*spectral_indices.nindx).transpose(1,2,0)

        # Consolidate to NOVALUE
        flgd = spectral_indices.bitmask.flagged(si_mask, flag=['MAIN_EMPTY', 'BLUE_EMPTY',
                                                               'RED_EMPTY', 'UNDEFINED_BANDS' ])
        mask[flgd] = self.bitmask.turn_on(mask[flgd], 'NOVALUE')

        # Consolidate to UNRELIABLE
        flgd = spectral_indices.bitmask.flagged(si_mask, flag=['MAIN_INCOMP', 'BLUE_INCOMP',
                                                               'RED_INCOMP' ])
        mask[flgd] = self.bitmask.turn_on(mask[flgd], 'UNRELIABLE')

        # Consolidate to NOCORRECTION
        # TODO: What do I do about this!
        flgd = spectral_indices.bitmask.flagged(si_mask, flag=['NO_DISPERSION_CORRECTION'])
        mask[flgd] = self.bitmask.turn_on(mask[flgd], 'NOCORRECTION')

        # Consolidate to MATHERROR
        flgd = spectral_indices.bitmask.flagged(si_mask, flag=['DIVBYZERO'])
        mask[flgd] = self.bitmask.turn_on(mask[flgd], 'MATHERROR')

#        if abs_only:
#            # Get the number of absorption-line indices. Indices *must*
#            # be ordered as [0:nabs], [nabs:nabs+nbhd]
#            nabs, nbhd = SpectralIndices.count_indices(spectral_indices.absdb,
#                                                       spectral_indices.bhddb)
#            mask[:,:,nabs:] = self.bitmask.turn_on(mask[:,:,nabs:], 'NOVALUE')

        return self._consolidate_donotuse(mask)

    @staticmethod
    def _get_kpc_per_arcsec(z):
        if z <= 0:
            warnings.warn('Systemic velocity <=0; radii in h^-1 kpc not provided.')
            return -1.

        H0 = 100 * astropy.units.km / astropy.units.s / astropy.units.Mpc
        cosmo = FlatLambdaCDM(H0=H0, Om0=0.3)
        return numpy.radians(1/3600) * 1e3 * cosmo.angular_diameter_distance(z).value

    def reduction_assessment_maps(self, prihdr, rdxqa):
        """
        Constructs the 'SPX_SKYCOO', 'SPX_ELLCOO', 'SPX_MFLUX',
        'SPX_MFLUX_IVAR', and 'SPX_SNR' map extensions.

        .. todo::
            Add the spaxel correlation data.
        """
        #---------------------------------------------------------------
        # Extensions filled from the reduction assessments
        ext = ['SPX_SKYCOO', 'SPX_ELLCOO', 'SPX_MFLUX', 'SPX_MFLUX_IVAR', 'SPX_SNR']

        if rdxqa is None:
            # Construct and return the empty hdus
            return DAPFitsUtil.empty_hdus(ext)

        #---------------------------------------------------------------
        # Add data to the primary header
        prihdr = rdxqa._initialize_primary_header(hdr=prihdr)

        #---------------------------------------------------------------
        # Get the extension headers
        hdr = [ DAPFitsUtil.finalize_dap_header(self.multichannel_maphdr, 'SPX_SKYCOO',
                                                bunit='arcsec', multichannel=True,
                                                channel_names=['On-sky X', 'On-sky Y']),
                DAPFitsUtil.finalize_dap_header(self.multichannel_maphdr, 'SPX_ELLCOO',
                                                multichannel=True,
                                                channel_names=['Elliptical radius', 'R/Re',
                                                               'R h/kpc', 'Elliptical azimuth'],
                                                channel_units=['arcsec', '', 'kpc/h',
                                                               'degrees']),
                DAPFitsUtil.finalize_dap_header(self.singlechannel_maphdr, 'SPX_MFLUX',
                                                bunit='1E-17 erg/s/cm^2/ang/spaxel', err=True),
                DAPFitsUtil.finalize_dap_header(self.singlechannel_maphdr, 'SPX_MFLUX',
                                                hduclas2='ERROR',
                                                bunit='(1E-17 erg/s/cm^2/ang/spaxel)^{-2}'),
                DAPFitsUtil.finalize_dap_header(self.singlechannel_maphdr, 'SPX_SNR')
              ]

        #---------------------------------------------------------------
        # Get the data arrays
        # On-sky coordinates
        minimum_value = numpy.finfo(self.float_dtype).eps
        spx_skycoo = rdxqa['SPECTRUM'].data['SKY_COO'].copy().reshape(*self.spatial_shape, -1)
        spx_skycoo[numpy.absolute(spx_skycoo) < minimum_value] = 0.0

        # Elliptical coordinates
        spx_ellcoo = rdxqa['SPECTRUM'].data['ELL_COO'].copy().reshape(*self.spatial_shape, -1)
        spx_ellcoo = numpy.repeat(spx_ellcoo, (3,1), axis=2)

        # Calculate the radius normalized by the effective radius
        if 'reff' in self.meta.keys() and self.meta['reff'] is not None:
            spx_ellcoo[:,:,1] /= self.meta['reff']

        # Calculate the radius in units of h^-1 kpc
        hkpc_per_arcsec = self.__class__._get_kpc_per_arcsec(self.meta['z'])
        if hkpc_per_arcsec > 0:
            spx_ellcoo[:,:,2] *= hkpc_per_arcsec
        else:
            spx_ellcoo[:,:,2] = -1
        spx_ellcoo[numpy.absolute(spx_ellcoo) < minimum_value] = 0.0

        # Bin signal
        signal = rdxqa['SPECTRUM'].data['SIGNAL'].copy().reshape(self.spatial_shape)
        signal[numpy.absolute(signal) < minimum_value] = 0.0
        # Bin inverse variance
        ivar = rdxqa['SPECTRUM'].data['VARIANCE'].copy().reshape(self.spatial_shape)
        ivar = numpy.ma.power(ivar, -1).filled(0.0)
        ivar[numpy.absolute(ivar) < minimum_value] = 0.0
        # Bin S/N
        snr = rdxqa['SPECTRUM'].data['SNR'].copy().reshape(self.spatial_shape)
        snr[numpy.absolute(snr) < minimum_value] = 0.0

        # Organize the extension data
        data = [ spx_skycoo.astype(self.float_dtype), spx_ellcoo.astype(self.float_dtype),
                 signal.astype(self.float_dtype), ivar.astype(self.float_dtype),
                 snr.astype(self.float_dtype) ]

        #---------------------------------------------------------------
        # Return the map hdus
        return DAPFitsUtil.list_of_image_hdus(data, hdr, ext)

    
    def binned_spectra_maps(self, prihdr, binned_spectra):
        """
        Constructs the 'BIN_LWSKYCOO', 'BIN_LWELLCOO', 'BIN_AREA',
        'BIN_FAREA', 'BIN_MFLUX', 'BIN_MFLUX_IVAR', 'BIN_MFLUX_MASK',
        and 'BIN_SNR' map extensions.
        """
        #---------------------------------------------------------------
        ext = ['BIN_LWSKYCOO', 'BIN_LWELLCOO', 'BIN_AREA', 'BIN_FAREA', 'BIN_MFLUX',
               'BIN_MFLUX_IVAR', 'BIN_MFLUX_MASK', 'BIN_SNR']

        if binned_spectra is None:
            # Construct and return the empty hdus
            return DAPFitsUtil.empty_hdus(ext)

        #---------------------------------------------------------------
        # Add data to the primary header
        # TODO: Apply this to each extension instead of the primary
        # header to allow for extension specific binning?
        prihdr = binned_spectra._initialize_primary_header(hdr=prihdr)
        if binned_spectra.is_unbinned:
            prihdr['BINTYPE'] = ('None', 'Binning method')
        else:
            prihdr = binned_spectra._add_method_header(prihdr)
        prihdr = binned_spectra._add_reddening_header(prihdr)

        #---------------------------------------------------------------
        # Get the extension headers
        hdr = [ DAPFitsUtil.finalize_dap_header(self.multichannel_maphdr, 'BIN_LWSKYCOO',
                                                bunit='arcsec', multichannel=True,
                                                channel_names=['Lum. weighted on-sky X',
                                                               'Lum. weighted on-sky Y']),
                DAPFitsUtil.finalize_dap_header(self.multichannel_maphdr, 'BIN_LWELLCOO',
                                                multichannel=True,
                                                channel_names= ['Lum. weighted elliptical radius',
                                                                'R/Re', 'R h/kpc',
                                                                'Lum. weighted elliptical azimuth'],
                                                channel_units=['arcsec', '', 'kpc/h', 'degrees']),
                DAPFitsUtil.finalize_dap_header(self.singlechannel_maphdr, 'BIN_AREA',
                                                bunit='arcsec^2'),
                DAPFitsUtil.finalize_dap_header(self.singlechannel_maphdr, 'BIN_FAREA'),
                DAPFitsUtil.finalize_dap_header(self.singlechannel_maphdr, 'BIN_MFLUX',
                                                bunit='1E-17 erg/s/cm^2/ang/spaxel', err=True,
                                                qual=True),
                DAPFitsUtil.finalize_dap_header(self.singlechannel_maphdr, 'BIN_MFLUX',
                                                hduclas2='ERROR', qual=True,
                                                bunit='(1E-17 erg/s/cm^2/ang/spaxel)^{-2}'),
                DAPFitsUtil.finalize_dap_header(self.singlechannel_maphdr, 'BIN_MFLUX',
                                                hduclas2='QUALITY', err=True,
                                                bit_type=self.bitmask.minimum_dtype()),
                DAPFitsUtil.finalize_dap_header(self.singlechannel_maphdr, 'BIN_SNR')
               ]

        #---------------------------------------------------------------
        # Get the data arrays
        arr = [ binned_spectra['BINS'].data['LW_SKY_COO'][:,0],
                binned_spectra['BINS'].data['LW_SKY_COO'][:,1],
                binned_spectra['BINS'].data['LW_ELL_COO'][:,0],
                binned_spectra['BINS'].data['LW_ELL_COO'][:,0],
                binned_spectra['BINS'].data['LW_ELL_COO'][:,0],
                binned_spectra['BINS'].data['LW_ELL_COO'][:,1],
                binned_spectra['BINS'].data['AREA'],
                binned_spectra['BINS'].data['AREA_FRAC'],
                binned_spectra['BINS'].data['SIGNAL'],
                numpy.ma.power(binned_spectra['BINS'].data['VARIANCE'].copy(), -1).filled(0.0),
                binned_spectra['BINS'].data['SNR']
              ]

        dtypes = [self.float_dtype]*len(arr)

        # Bin index
        bin_indx = binned_spectra['BINID'].data.copy().ravel()

        # Remap the data to the DRP spatial shape
        arr = list(DAPFitsUtil.reconstruct_map(self.spatial_shape, bin_indx, arr, dtype=dtypes))

        # Calculate the radius normalized by the effective radius
        if 'reff' in self.meta.keys() and self.meta['reff'] is not None:
            arr[3] = (arr[3]/self.meta['reff']).astype(self.float_dtype)
        # Calculate the radius in units of h^-1 kpc
        hkpc_per_arcsec = self.__class__._get_kpc_per_arcsec(self.meta['z'])
        arr[4] = (arr[4]*hkpc_per_arcsec).astype(self.float_dtype) if hkpc_per_arcsec > 0 \
                        else (0*arr[4]-1).astype(self.float_dtype)

        # Organize the extension data
        data = [ numpy.array(arr[0:2]).transpose(1,2,0), numpy.array(arr[2:6]).transpose(1,2,0) ] \
                    + arr[6:-1] + [ self.bin_mask.copy(), arr[-1] ]

        #---------------------------------------------------------------
        # Return the map hdus
        return DAPFitsUtil.list_of_image_hdus(data, hdr, ext)


    @staticmethod
    def _join_maps(a):
        return numpy.array(a).transpose(1,2,0)

    def stellar_continuum_maps(self, prihdr, stellar_continuum):
        """
        Construct the 'STELLAR_VEL', 'STELLAR_VEL_IVAR',
        'STELLAR_VEL_MASK', 'STELLAR_SIGMA', 'STELLAR_SIGMA_IVAR',
        'STELLAR_SIGMA_MASK', 'STELLAR_SIGMACORR',
        'STELLAR_FOM' maps extensions.
        """
        #---------------------------------------------------------------
        ext = [ 'STELLAR_VEL', 'STELLAR_VEL_IVAR', 'STELLAR_VEL_MASK', 'STELLAR_SIGMA',
                'STELLAR_SIGMA_IVAR', 'STELLAR_SIGMA_MASK', 'STELLAR_SIGMACORR', 'STELLAR_FOM' ]

        if stellar_continuum is None:
            # Construct and return the empty hdus
            return DAPFitsUtil.empty_hdus(ext)

        #---------------------------------------------------------------
        # Add data to the primary header
        prihdr = stellar_continuum._initialize_primary_header(hdr=prihdr)
        prihdr = stellar_continuum._add_method_header(prihdr)

        # Figure-of-merit units
        fomunits = ['']*9
        fomunits[0] = '1E-17 erg/s/cm^2/ang/spaxel'

        #---------------------------------------------------------------
        # Get the extension headers
        hdr = [ DAPFitsUtil.finalize_dap_header(self.singlechannel_maphdr, 'STELLAR_VEL',
                                                bunit='km/s', err=True, qual=True),
                DAPFitsUtil.finalize_dap_header(self.singlechannel_maphdr, 'STELLAR_VEL',
                                                hduclas2='ERROR', bunit='(km/s)^{-2}', qual=True),
                DAPFitsUtil.finalize_dap_header(self.singlechannel_maphdr, 'STELLAR_VEL',
                                                hduclas2='QUALITY', err=True,
                                                bit_type=self.bitmask.minimum_dtype()),
                DAPFitsUtil.finalize_dap_header(self.singlechannel_maphdr, 'STELLAR_SIGMA',
                                                bunit='km/s', err=True, qual=True),
                DAPFitsUtil.finalize_dap_header(self.singlechannel_maphdr, 'STELLAR_SIGMA',
                                                hduclas2='ERROR', bunit='(km/s)^{-2}', qual=True),
                DAPFitsUtil.finalize_dap_header(self.singlechannel_maphdr, 'STELLAR_SIGMA',
                                                hduclas2='QUALITY', err=True,
                                                bit_type=self.bitmask.minimum_dtype()),
#                DAPFitsUtil.finalize_dap_header(self.singlechannel_maphdr, 'STELLAR_SIGMACORR',
#                                                bunit='km/s'),
                DAPFitsUtil.finalize_dap_header(self.multichannel_maphdr, 'STELLAR_SIGMACORR',
                                                multichannel=True, bunit='km/s',
                                                channel_names=['resolution difference', 'fit']),
                DAPFitsUtil.finalize_dap_header(self.multichannel_maphdr, 'STELLAR_FOM',
                                                multichannel=True,
                                                channel_names=['rms', 'frms', 'rchi2',
                                                               '68th perc frac resid',
                                                               '99th perc frac resid',
                                                               'max frac resid',
                                                               '68th perc per pix chi',
                                                               '99th perc per pix chi',
                                                               'max per pix chi'],
                                                channel_units=fomunits)
              ]


        #---------------------------------------------------------------
        # Get the data arrays
        arr = [ DAPFitsUtil.redshift_to_Newtonian_velocity(
                                    stellar_continuum['PAR'].data['KIN'][:,0], self.redshift),
                DAPFitsUtil.redshift_to_Newtonian_velocity(numpy.ma.power(
                                                    stellar_continuum['PAR'].data['KINERR'][:,0],
                                                    -2).filled(0.0), self.redshift, ivar=True),
                stellar_continuum['PAR'].data['KIN'][:,1],
                numpy.ma.power(stellar_continuum['PAR'].data['KINERR'][:,1], -2).filled(0.0),
                stellar_continuum['PAR'].data['SIGMACORR_EMP'],
                stellar_continuum['PAR'].data['SIGMACORR_SRES'],
                stellar_continuum['PAR'].data['RMS'],
                stellar_continuum['PAR'].data['FRMS'],
                stellar_continuum['PAR'].data['RCHI2'],
                stellar_continuum['PAR'].data['FRMSGRW'][:,1],
                stellar_continuum['PAR'].data['FRMSGRW'][:,3],
                stellar_continuum['PAR'].data['FRMSGRW'][:,4],
                stellar_continuum['PAR'].data['CHIGRW'][:,1],
                stellar_continuum['PAR'].data['CHIGRW'][:,3],
                stellar_continuum['PAR'].data['CHIGRW'][:,4],
                stellar_continuum['PAR'].data['MASK']
              ]

        # Data types
        dtypes = [self.float_dtype]*(len(arr)-1) + [arr[-1].dtype.name]

        # Bin index
        bin_indx = stellar_continuum['BINID'].data.copy().ravel()

        # Remap the data to the DRP spatial shape
        arr = list(DAPFitsUtil.reconstruct_map(self.spatial_shape, bin_indx, arr, dtype=dtypes))
        
        # Get the masks
        vel_mask = self._stellar_continuum_mask_to_map_mask(stellar_continuum, arr[-1].copy())
        sig_mask = self._stellar_continuum_mask_to_map_mask(stellar_continuum, arr[-1].copy(),
                                                            for_dispersion=True)

        # Organize the extension data
        data = arr[0:2] + [ vel_mask ] + arr[2:4] + [ sig_mask ] \
                    + [ self._join_maps(arr[4:6]), self._join_maps(arr[6:-1]) ]

        #---------------------------------------------------------------
        # Return the map hdus
        return DAPFitsUtil.list_of_image_hdus(data, hdr, ext)


    def emission_line_moment_maps(self, prihdr, emission_line_moments):
        """
        Construct the 'EMLINE_SFLUX', 'EMLINE_SFLUX_IVAR',
        'EMLINE_SFLUX_MASK', 'EMLINE_SEW', 'EMLINE_SEW_IVAR', and
        'EMLINE_SEW_MASK' maps extensions.
        """
        #---------------------------------------------------------------
        ext = [ 'EMLINE_SFLUX', 'EMLINE_SFLUX_IVAR', 'EMLINE_SFLUX_MASK', 'EMLINE_SEW',
                'EMLINE_SEW_CNT', 'EMLINE_SEW_IVAR', 'EMLINE_SEW_MASK']

        if emission_line_moments is None:
            # Construct and return the empty hdus
            return DAPFitsUtil.empty_hdus(ext)

        #---------------------------------------------------------------
        # Add data to the primary header
        prihdr = emission_line_moments._initialize_primary_header(hdr=prihdr)

        #---------------------------------------------------------------
        # Get the extension headers

        # Need multichannel map header if more than one moment
        multichannel = emission_line_moments.nmom > 1
        # Build the channel names
#        names = [ '{0}-{1}'.format(n,int(w)) \
#                        for n,w in zip(emission_line_moments['ELMBAND'].data['NAME'],
#                                       emission_line_moments['ELMBAND'].data['RESTWAVE']) ]
        names = emission_line_moments.channel_names()
        # Create the basic header for all extensions
        base_hdr = DAPFitsUtil.add_channel_names(self.multichannel_maphdr if multichannel
                                                 else self.singlechannel_maphdr, names)

        hdr = [ DAPFitsUtil.finalize_dap_header(base_hdr, 'EMLINE_SFLUX', err=True, qual=True,
                                                bunit='1E-17 erg/s/cm^2/spaxel',
                                                multichannel=multichannel),
                DAPFitsUtil.finalize_dap_header(base_hdr, 'EMLINE_SFLUX', hduclas2='ERROR',
                                                qual=True, bunit='(1E-17 erg/s/cm^2/spaxel)^{-2}',
                                                multichannel=multichannel),
                DAPFitsUtil.finalize_dap_header(base_hdr, 'EMLINE_SFLUX', hduclas2='QUALITY',
                                                err=True, bit_type=self.bitmask.minimum_dtype(),
                                                multichannel=multichannel),
                DAPFitsUtil.finalize_dap_header(base_hdr, 'EMLINE_SEW', err=True, qual=True,
                                                bunit='ang', multichannel=multichannel),
                DAPFitsUtil.finalize_dap_header(base_hdr, 'EMLINE_SEW_CNT',
                                                bunit='1E-17 erg/s/cm^2/ang/spaxel',
                                                multichannel=multichannel),
                DAPFitsUtil.finalize_dap_header(base_hdr, 'EMLINE_SEW', hduclas2='ERROR', qual=True,
                                                bunit='(ang)^{-2}', multichannel=multichannel),
                DAPFitsUtil.finalize_dap_header(base_hdr, 'EMLINE_SEW', hduclas2='QUALITY',
                                                err=True, bit_type=self.bitmask.minimum_dtype(),
                                                multichannel=multichannel)
              ]

        #---------------------------------------------------------------
        # Get the data arrays
        arr = [ emission_line_moments['ELMMNTS'].data['FLUX'][:,m]
                    for m in range(emission_line_moments.nmom) ]
        arr += [ numpy.ma.power(emission_line_moments['ELMMNTS'].data['FLUXERR'][:,m],
                                -2).filled(0.0) for m in range(emission_line_moments.nmom) ]
        arr += [ emission_line_moments['ELMMNTS'].data['EW'][:,m]
                    for m in range(emission_line_moments.nmom) ]
        arr += [ emission_line_moments['ELMMNTS'].data['EWCONT'][:,m]
                    for m in range(emission_line_moments.nmom) ]
        arr += [ numpy.ma.power(emission_line_moments['ELMMNTS'].data['EWERR'][:,m],
                                -2).filled(0.0) for m in range(emission_line_moments.nmom) ]
        arr += [ emission_line_moments['ELMMNTS'].data['MASK'][:,m]
                    for m in range(emission_line_moments.nmom) ]

        dtypes = [self.float_dtype]*(len(arr)-emission_line_moments.nmom) \
                        + [a.dtype.name for a in arr[-emission_line_moments.nmom:]]

        # Bin index
        bin_indx = emission_line_moments['BINID'].data.copy().ravel()

        # Remap the data to the DRP spatial shape
        arr = list(DAPFitsUtil.reconstruct_map(self.spatial_shape, bin_indx, arr, dtype=dtypes))

        data = [ numpy.array(arr[emission_line_moments.nmom*i:
                                 emission_line_moments.nmom*(i+1)]).transpose(1,2,0) 
                        for i in range(6) ]
        
        # Get the mask
        elm_mask = self._emission_line_moment_mask_to_map_mask(emission_line_moments,
                                                               data[-1].copy())

        # Organize the extension data
        data = data[:2] + [elm_mask] + data[2:-1] + [elm_mask]

        #---------------------------------------------------------------
        # Return the map hdus
        return DAPFitsUtil.list_of_image_hdus(data, hdr, ext)


    def emission_line_model_maps(self, prihdr, emission_line_model):
        """
        Construct the 'EMLINE_GFLUX', 'EMLINE_GFLUX_IVAR',
        'EMLINE_GFLUX_MASK', 'EMLINE_GEW', 'EMLINE_GEW_CNT',
        'EMLINE_GEW_IVAR', 'EMLINE_GEW_MASK', 'EMLINE_GVEL',
        'EMLINE_GVEL_IVAR', 'EMLINE_GVEL_MASK', 'EMLINE_GSIGMA',
        'EMLINE_GSIGMA_IVAR', 'EMLINE_GSIGMA_MASK',
        'EMLINE_INSTSIGMA', 'EMLINE_TPLSIGMA', 'EMLINE_GA',
        'EMLINE_GANR', 'EMLINE_FOM', 'EMLINE_LFOM' map extensions.
        """
        #---------------------------------------------------------------
        ext = [ 'EMLINE_GFLUX', 'EMLINE_GFLUX_IVAR', 'EMLINE_GFLUX_MASK', 'EMLINE_GEW',
                'EMLINE_GEW_CNT', 'EMLINE_GEW_IVAR', 'EMLINE_GEW_MASK', 'EMLINE_GVEL',
                'EMLINE_GVEL_IVAR', 'EMLINE_GVEL_MASK', 'EMLINE_GSIGMA', 'EMLINE_GSIGMA_IVAR',
                'EMLINE_GSIGMA_MASK', 'EMLINE_INSTSIGMA', 'EMLINE_TPLSIGMA', 'EMLINE_GA',
                'EMLINE_GANR', 'EMLINE_FOM', 'EMLINE_LFOM' ]

        if emission_line_model is None:
            # Construct and return the empty hdus
            return DAPFitsUtil.empty_hdus(ext)

        #---------------------------------------------------------------
        # Add data to the primary header
        prihdr = emission_line_model._initialize_primary_header(hdr=prihdr)
        prihdr = emission_line_model._add_method_header(prihdr)

        #---------------------------------------------------------------
        # Get the extension headers

        # Need multichannel map header if more than one moment
        multichannel = emission_line_model.neml > 1

        # Create the basic header for all extensions
        if 'data' in emission_line_model['PAR'].__dict__:
            names = [ '{0}-{1}'.format(n,int(w)) \
                        for n,w in zip(emission_line_model['PAR'].data['NAME'],
                                       emission_line_model['PAR'].data['RESTWAVE']) ]
            base_hdr = DAPFitsUtil.add_channel_names(self.multichannel_maphdr if multichannel
                                                     else self.singlechannel_maphdr, names)
        else:
            # TODO: This should throw a ValueError, not just a warning
            warnings.warn('Emission-line model does not include channel names!')
            names = None
            base_hdr = self.multichannel_maphdr if multichannel else self.singlechannel_maphdr

        # Get the figure of merit data to output
        fom_names, fom_units, fom_data = emission_line_model.fit_figures_of_merit()
       
        hdr = [ DAPFitsUtil.finalize_dap_header(base_hdr, 'EMLINE_GFLUX', err=True, qual=True,
                                                bunit='1E-17 erg/s/cm^2/spaxel',
                                                multichannel=multichannel),
                DAPFitsUtil.finalize_dap_header(base_hdr, 'EMLINE_GFLUX', hduclas2='ERROR',
                                                qual=True, bunit='(1E-17 erg/s/cm^2/spaxel)^{-2}',
                                                multichannel=multichannel),
                DAPFitsUtil.finalize_dap_header(base_hdr, 'EMLINE_GFLUX', hduclas2='QUALITY',
                                                err=True, bit_type=self.bitmask.minimum_dtype(),
                                                multichannel=multichannel),
                DAPFitsUtil.finalize_dap_header(base_hdr, 'EMLINE_GEW', err=True, qual=True,
                                                bunit='ang', multichannel=multichannel),
                DAPFitsUtil.finalize_dap_header(base_hdr, 'EMLINE_GEW_CNT',
                                                bunit='1E-17 erg/s/cm^2/ang/spaxel',
                                                multichannel=multichannel),
                DAPFitsUtil.finalize_dap_header(base_hdr, 'EMLINE_GEW', hduclas2='ERROR', qual=True,
                                                bunit='(ang)^{-2}', multichannel=multichannel),
                DAPFitsUtil.finalize_dap_header(base_hdr, 'EMLINE_GEW', hduclas2='QUALITY',
                                                err=True, bit_type=self.bitmask.minimum_dtype(),
                                                multichannel=multichannel),
                DAPFitsUtil.finalize_dap_header(base_hdr, 'EMLINE_GVEL', err=True, qual=True,
                                                bunit='km/s', multichannel=multichannel),
                DAPFitsUtil.finalize_dap_header(base_hdr, 'EMLINE_GVEL', hduclas2='ERROR',
                                                qual=True, bunit='(km/s)^{-2}',
                                                multichannel=multichannel),
                DAPFitsUtil.finalize_dap_header(base_hdr, 'EMLINE_GVEL', hduclas2='QUALITY',
                                                err=True, bit_type=self.bitmask.minimum_dtype(),
                                                multichannel=multichannel),
                DAPFitsUtil.finalize_dap_header(base_hdr, 'EMLINE_GSIGMA', err=True, qual=True,
                                                bunit='km/s', multichannel=multichannel),
                DAPFitsUtil.finalize_dap_header(base_hdr, 'EMLINE_GSIGMA', hduclas2='ERROR',
                                                qual=True, bunit='(km/s)^{-2}',
                                                multichannel=multichannel),
                DAPFitsUtil.finalize_dap_header(base_hdr, 'EMLINE_GSIGMA', hduclas2='QUALITY',
                                                err=True, bit_type=self.bitmask.minimum_dtype(),
                                                multichannel=multichannel),
                DAPFitsUtil.finalize_dap_header(base_hdr, 'EMLINE_INSTSIGMA', bunit='km/s',
                                                multichannel=multichannel),
                DAPFitsUtil.finalize_dap_header(base_hdr, 'EMLINE_TPLSIGMA', bunit='km/s',
                                                multichannel=multichannel),
                DAPFitsUtil.finalize_dap_header(base_hdr, 'EMLINE_GA',
                                                bunit='1E-17 erg/s/cm^2/ang/spaxel',
                                                multichannel=multichannel),
                DAPFitsUtil.finalize_dap_header(base_hdr, 'EMLINE_GANR',
                                                multichannel=multichannel),
                DAPFitsUtil.finalize_dap_header(self.multichannel_maphdr, 'EMLINE_FOM',
                                                multichannel=True, channel_names=fom_names,
                                                channel_units=fom_units),
                DAPFitsUtil.finalize_dap_header(base_hdr, 'EMLINE_LFOM', multichannel=multichannel)
              ]

        #---------------------------------------------------------------
        # Get the data arrays
        n_arr_with_eml_channels = 0
        # Fluxes
        arr = [ emission_line_model['EMLDATA'].data['FLUX'][:,m]
                    for m in range(emission_line_model.neml) ]
        n_arr_with_eml_channels += 1
        # Flux errors
        arr += [ numpy.ma.power(emission_line_model['EMLDATA'].data['FLUXERR'][:,m],
                                -2).filled(0.0) for m in range(emission_line_model.neml) ]
        n_arr_with_eml_channels += 1
        # Equivalent width
        arr += [ emission_line_model['EMLDATA'].data['EW'][:,m]
                    for m in range(emission_line_model.neml) ]
        n_arr_with_eml_channels += 1
        # Equivalent width continuum
        arr += [ emission_line_model['EMLDATA'].data['EWCONT'][:,m]
                    for m in range(emission_line_model.neml) ]
        n_arr_with_eml_channels += 1
        # Equivalent width errors
        arr += [ numpy.ma.power(emission_line_model['EMLDATA'].data['EWERR'][:,m],
                                -2).filled(0.0) for m in range(emission_line_model.neml) ]
        n_arr_with_eml_channels += 1
        # Velocities
        arr += [ DAPFitsUtil.redshift_to_Newtonian_velocity(
                        emission_line_model['EMLDATA'].data['KIN'][:,m,0], self.redshift)
                    for m in range(emission_line_model.neml) ]
        n_arr_with_eml_channels += 1
        # Velocity errors
        arr += [ DAPFitsUtil.redshift_to_Newtonian_velocity(
                        numpy.ma.power(emission_line_model['EMLDATA'].data['KINERR'][:,m,0],
                                       -2).filled(0.0), self.redshift, ivar=True) \
                    for m in range(emission_line_model.neml) ]
        n_arr_with_eml_channels += 1
        # Velocity dispersions
        arr += [ emission_line_model['EMLDATA'].data['KIN'][:,m,1]
                    for m in range(emission_line_model.neml) ]
        n_arr_with_eml_channels += 1
        # Velocity dispersion errors
        arr += [ numpy.ma.power(emission_line_model['EMLDATA'].data['KINERR'][:,m,1],
                                -2).filled(0.0) for m in range(emission_line_model.neml) ]
        n_arr_with_eml_channels += 1
        # Instrumental dispersions
        arr += [ emission_line_model['EMLDATA'].data['SIGMAINST'][:,m]
                    for m in range(emission_line_model.neml) ]
        n_arr_with_eml_channels += 1
        # Template velocity dispersions
        arr += [ emission_line_model['EMLDATA'].data['SIGMATPL'][:,m]
                    for m in range(emission_line_model.neml) ]
        n_arr_with_eml_channels += 1
        # Amplitude
        arr += [ emission_line_model['EMLDATA'].data['AMP'][:,m]
                    for m in range(emission_line_model.neml) ]
        n_arr_with_eml_channels += 1
        # A/N
        arr += [ emission_line_model['EMLDATA'].data['ANR'][:,m]
                    for m in range(emission_line_model.neml) ]
        n_arr_with_eml_channels += 1
        # Line reduced chi-square
        mp = 3  # Number of model parameters, TODO: get this from emission_line_model...
        reduced_chi2 = numpy.ma.divide(emission_line_model['EMLDATA'].data['LINE_CHI2'],
                                       emission_line_model['EMLDATA'].data['LINE_NSTAT'] 
                                            - mp).filled(-999.)
        arr += [ reduced_chi2[:,m] for m in range(emission_line_model.neml) ]
        n_arr_with_eml_channels += 1
        # Full spectrum fit-quality metrics (if available)
        nfom = len(fom_names)
        arr += [ fom_data[:,m] for m in range(nfom) ]
        # Mask data
        arr += [ emission_line_model['EMLDATA'].data['MASK'][:,m]
                    for m in range(emission_line_model.neml) ]

        # Data types
        dtypes = [self.float_dtype]*(n_arr_with_eml_channels*emission_line_model.neml + nfom) \
                        + [a.dtype.name for a in arr[-emission_line_model.neml:]]

        # Bin index
        bin_indx = emission_line_model['BINID'].data.copy().ravel()

        # Remap the data to the DRP spatial shape
        arr = list(DAPFitsUtil.reconstruct_map(self.spatial_shape, bin_indx, arr, dtype=dtypes))

        # Join the maps for each emission line in the set of quantities 
        data = [ self._join_maps(arr[emission_line_model.neml*i:emission_line_model.neml*(i+1)])
                        for i in range(n_arr_with_eml_channels) ]
        data += [ self._join_maps(arr[-emission_line_model.neml-nfom:-emission_line_model.neml]) ]
        data += [ self._join_maps(arr[-emission_line_model.neml:]) ]

        # Get the masks
        base_mask = self._emission_line_model_mask_to_map_mask(emission_line_model, data[-1].copy())
        sig_mask = self._emission_line_model_mask_to_map_mask(emission_line_model, data[-1].copy(),
                                                              for_dispersion=True)

        # Organize the extension data
        data = data[:2] + [ base_mask ] + data[2:5] + [ base_mask ] + data[5:7] + [ base_mask ] \
                + data[7:9] + [ sig_mask ] + data[9:-3] + [data[-2], data[-3]]

        #---------------------------------------------------------------
        # Return the map hdus
        return DAPFitsUtil.list_of_image_hdus(data, hdr, ext)


    def spectral_index_maps(self, prihdr, spectral_indices):
        """
        Construct the 'SPECINDEX', 'SPECINDEX_IVAR', 'SPECINDEX_MASK',
        and 'SPECINDEX_CORR'.
        """
        #---------------------------------------------------------------
#        ext = [ 'SPECINDEX', 'SPECINDEX_IVAR', 'SPECINDEX_MASK', 'SPECINDEX_CORR',
#                'SPECINDEX_BCEN', 'SPECINDEX_BCNT', 'SPECINDEX_RCEN', 'SPECINDEX_RCNT',
#                'SPECINDEX_MODEL' ]
        ext = [ 'SPECINDEX', 'SPECINDEX_IVAR', 'SPECINDEX_MASK', 'SPECINDEX_CORR',
                'SPECINDEX_MODEL', 'SPECINDEX_BF', 'SPECINDEX_BF_IVAR', 'SPECINDEX_BF_MASK',
                'SPECINDEX_BF_CORR', 'SPECINDEX_BF_MODEL', 'SPECINDEX_WGT', 'SPECINDEX_WGT_IVAR',
                'SPECINDEX_WGT_MASK', 'SPECINDEX_WGT_CORR', 'SPECINDEX_WGT_MODEL']

        if spectral_indices is None:
            # Construct and return the empty hdus
            return DAPFitsUtil.empty_hdus(ext)

        #---------------------------------------------------------------
        # Add data to the primary header
        prihdr = spectral_indices._initialize_primary_header(hdr=prihdr)

        #---------------------------------------------------------------
        # Get the extension headers

        # Need multichannel map header if more than one index
        multichannel = spectral_indices.nindx > 1
        # Create the basic header for all extensions
        base_hdr = self.multichannel_maphdr if multichannel else self.singlechannel_maphdr

        errunits = [ '({0})'.format(u)+'^{-2}' if len(u) > 0 else ''
                            for u in spectral_indices['SIPAR'].data['UNIT'] ]

        hdr = [# Worthey/Trager indices
               DAPFitsUtil.finalize_dap_header(base_hdr, 'SPECINDEX', err=True, qual=True,
                                               multichannel=multichannel,
                                              channel_names=spectral_indices['SIPAR'].data['NAME'],
                                              channel_units=spectral_indices['SIPAR'].data['UNIT']),
               DAPFitsUtil.finalize_dap_header(base_hdr, 'SPECINDEX', hduclas2='ERROR', qual=True,
                                               multichannel=multichannel,
                                              channel_names=spectral_indices['SIPAR'].data['NAME'],
                                               channel_units=errunits),
               DAPFitsUtil.finalize_dap_header(base_hdr, 'SPECINDEX', hduclas2='QUALITY', err=True,
                                               bit_type=self.bitmask.minimum_dtype(),
                                               multichannel=multichannel,
                                              channel_names=spectral_indices['SIPAR'].data['NAME']),
               DAPFitsUtil.finalize_dap_header(base_hdr, 'SPECINDEX_CORR',
                                               multichannel=multichannel,
                                              channel_names=spectral_indices['SIPAR'].data['NAME']),
               DAPFitsUtil.finalize_dap_header(base_hdr, 'SPECINDEX_MODEL',
                                               multichannel=multichannel,
                                              channel_names=spectral_indices['SIPAR'].data['NAME']),
               # Burstein/Faber indices
               DAPFitsUtil.finalize_dap_header(base_hdr, 'SPECINDEX_BF', err=True, qual=True,
                                               multichannel=multichannel,
                                              channel_names=spectral_indices['SIPAR'].data['NAME'],
                                              channel_units=spectral_indices['SIPAR'].data['UNIT']),
               DAPFitsUtil.finalize_dap_header(base_hdr, 'SPECINDEX_BF', hduclas2='ERROR', qual=True,
                                               multichannel=multichannel,
                                              channel_names=spectral_indices['SIPAR'].data['NAME'],
                                               channel_units=errunits),
               DAPFitsUtil.finalize_dap_header(base_hdr, 'SPECINDEX_BF', hduclas2='QUALITY', err=True,
                                               bit_type=self.bitmask.minimum_dtype(),
                                               multichannel=multichannel,
                                              channel_names=spectral_indices['SIPAR'].data['NAME']),
               DAPFitsUtil.finalize_dap_header(base_hdr, 'SPECINDEX_BF_CORR',
                                               multichannel=multichannel,
                                              channel_names=spectral_indices['SIPAR'].data['NAME']),
               DAPFitsUtil.finalize_dap_header(base_hdr, 'SPECINDEX_BF_MODEL',
                                               multichannel=multichannel,
                                              channel_names=spectral_indices['SIPAR'].data['NAME']),
               # Continuum in the main passband
               DAPFitsUtil.finalize_dap_header(base_hdr, 'SPECINDEX_WGT', err=True, qual=True,
                                               multichannel=multichannel,
                                              channel_names=spectral_indices['SIPAR'].data['NAME'],
                                              channel_units=spectral_indices['SIPAR'].data['UNIT']),
               DAPFitsUtil.finalize_dap_header(base_hdr, 'SPECINDEX_WGT', hduclas2='ERROR', qual=True,
                                               multichannel=multichannel,
                                              channel_names=spectral_indices['SIPAR'].data['NAME'],
                                               channel_units=errunits),
               DAPFitsUtil.finalize_dap_header(base_hdr, 'SPECINDEX_WGT', hduclas2='QUALITY', err=True,
                                               bit_type=self.bitmask.minimum_dtype(),
                                               multichannel=multichannel,
                                              channel_names=spectral_indices['SIPAR'].data['NAME']),
               DAPFitsUtil.finalize_dap_header(base_hdr, 'SPECINDEX_WGT_CORR',
                                               multichannel=multichannel,
                                              channel_names=spectral_indices['SIPAR'].data['NAME']),
               DAPFitsUtil.finalize_dap_header(base_hdr, 'SPECINDEX_WGT_MODEL',
                                               multichannel=multichannel,
                                              channel_names=spectral_indices['SIPAR'].data['NAME'])
              ]

        #---------------------------------------------------------------
        # Get the data arrays
        arr = [ spectral_indices['SINDX'].data['INDX'][:,m]
                    for m in range(spectral_indices.nindx) ]
        arr += [ numpy.ma.power(spectral_indices['SINDX'].data['INDX_ERR'][:,m], -2.).filled(0.0)
                    for m in range(spectral_indices.nindx) ]
        arr += [ spectral_indices['SINDX'].data['INDX_CORR'][:,m]
                    for m in range(spectral_indices.nindx) ]
        arr += [ spectral_indices['SINDX'].data['INDX_MOD'][:,m]
                    for m in range(spectral_indices.nindx) ]

        arr += [ spectral_indices['SINDX'].data['INDX_BF'][:,m]
                    for m in range(spectral_indices.nindx) ]
        arr += [ numpy.ma.power(spectral_indices['SINDX'].data['INDX_BF_ERR'][:,m], -2.).filled(0.0)
                    for m in range(spectral_indices.nindx) ]
        arr += [ spectral_indices['SINDX'].data['INDX_BF_CORR'][:,m]
                    for m in range(spectral_indices.nindx) ]
        arr += [ spectral_indices['SINDX'].data['INDX_BF_MOD'][:,m]
                    for m in range(spectral_indices.nindx) ]

        arr += [ spectral_indices['SINDX'].data['AWGT'][:,m]
                    for m in range(spectral_indices.nindx) ]
        arr += [ numpy.ma.power(spectral_indices['SINDX'].data['AWGT_ERR'][:,m], -2.).filled(0.0)
                    for m in range(spectral_indices.nindx) ]
        arr += [ spectral_indices['SINDX'].data['AWGT_CORR'][:,m]
                    for m in range(spectral_indices.nindx) ]
        arr += [ spectral_indices['SINDX'].data['AWGT_MOD'][:,m]
                    for m in range(spectral_indices.nindx) ]

        arr += [ spectral_indices['SINDX'].data['MASK'][:,m]
                    for m in range(spectral_indices.nindx) ]

        dtypes = [self.float_dtype]*(len(arr)-spectral_indices.nindx) \
                        + [a.dtype.name for a in arr[-spectral_indices.nindx:]]

        # Bin index
        bin_indx = spectral_indices['BINID'].data.copy().ravel()

        # Remap the data to the DRP spatial shape
        arr = list(DAPFitsUtil.reconstruct_map(self.spatial_shape, bin_indx, arr, dtype=dtypes))

        data = [ numpy.array(arr[spectral_indices.nindx*i:
                                 spectral_indices.nindx*(i+1)]).transpose(1,2,0) \
                        for i in range(13) ]

        # Get the mask
        si_mask = self._spectral_index_mask_to_map_mask(spectral_indices, data[-1].copy())

        # Organize the extension data
        data = data[:2] + [si_mask] + data[2:6] + [si_mask] + data[6:10] + [si_mask] + data[10:-1]

        #---------------------------------------------------------------
        # Return the map hdus
        return DAPFitsUtil.list_of_image_hdus(data, hdr, ext)


#-----------------------------------------------------------------------
class construct_cube_file:
    """
    Construct a DAP cube file based on the input data.

    Set as a class for coherency reasons, but should not be used as an
    object!

    Should force all intermediate objects to be provided.

    """
    def __init__(self, cube, method=None, rdxqa=None, binned_spectra=None, stellar_continuum=None,
                 emission_line_moments=None, emission_line_model=None, spectral_indices=None,
                 output_path=None, output_file=None, overwrite=True, loggers=None, quiet=False,
                 single_precision=False):

        #---------------------------------------------------------------
        # Initialize the reporting
        self.loggers = None if loggers is None else loggers
        self.quiet = quiet

        self.float_dtype = 'float32' if single_precision else 'float'

        #---------------------------------------------------------------
        # Check input types
        confirm_dap_types(cube, None, binned_spectra, stellar_continuum, None,
                          emission_line_model, None)

        #---------------------------------------------------------------
        # Set the output paths
        self.cube = cube
        self.meta = cube.meta
        self.method, self.directory_path, self.output_file \
                = default_paths(cube, file_type='LOGCUBE',
                        rdxqa_method=None if rdxqa is None else rdxqa.method['key'],
                        binning_method=None if binned_spectra is None 
                                       else binned_spectra.method['key'],
                        stelcont_method=None if stellar_continuum is None
                                        else stellar_continuum.method['key'],
                        elmom_method=None if emission_line_moments is None
                                     else emission_line_moments.database['key'],
                        elmodel_method=None if emission_line_model is None
                                       else emission_line_model.method['key'],
                        sindex_method=None if spectral_indices is None
                                      else spectral_indices.database['key'],
                        method=method, output_path=output_path, output_file=output_file)

        # Save input for reference
        self.shape = self.cube.shape
        self.spatial_shape = self.cube.spatial_shape

        # Report
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, '-'*50)
            log_output(loggers, 1, logging.INFO, '{0:^50}'.format('CONSTRUCTING OUTPUT MODEL CUBE'))
            log_output(self.loggers, 1, logging.INFO, '-'*50)
            log_output(self.loggers, 1, logging.INFO,
                       'Output path: {0}'.format(self.directory_path))
            log_output(self.loggers, 1, logging.INFO,
                       'Output file: {0}'.format(self.output_file))
            log_output(self.loggers, 1, logging.INFO,
                       'Output cubes have shape {0}'.format(self.shape))

        #---------------------------------------------------------------
        # Check if the file already exists
        ofile = self.directory_path / self.output_file
        if ofile.exists() and not overwrite:
            # TODO: Perform some checks to make sure the existing file
            # has the correct content?
            warnings.warn('Output file exists!  Set overwrite=True to overwrite.')
            return

        #---------------------------------------------------------------
        # Construct the 3D data arrays
        binned_spectra_3d_hdu = None if binned_spectra is None \
                                        else binned_spectra.construct_3d_hdu()
        stellar_continuum_3d_hdu = None if stellar_continuum is None \
                                        else stellar_continuum.construct_3d_hdu()
        emission_line_model_3d_hdu = None if emission_line_model is None \
                                        else emission_line_model.construct_3d_hdu()

        #---------------------------------------------------------------
        # Initialize the primary header
        prihdr = DAPFitsUtil.initialize_dap_primary_header(self.cube, maskname='MANGA_DAPSPECMASK')
        # Add the DAP method
        prihdr['DAPTYPE'] = (self.method, 'DAP analysis method')
        # Add the format of this file
        prihdr['DAPFRMT'] = ('LOGCUBE', 'DAP data file format')

        # Get the base map header
        self.base_cubehdr = DAPFitsUtil.build_cube_header(self.cube,
                                'K Westfall <westfall@ucolick.org> & SDSS-IV Data Group',
                                                        maskname='MANGA_DAPSPECMASK')

        #---------------------------------------------------------------
        # Initialize the pixel mask
        self.bitmask = DAPCubeBitMask()

        #---------------------------------------------------------------
        # Construct the hdu list for each input object.
        # This adds the binning-specific flags to the mask
        # Binned spectra: ['FLUX', 'IVAR', 'MASK', 'WAVE', 'REDCORR']
        prihdr, binlist = self.binned_data_cube(prihdr, binned_spectra, binned_spectra_3d_hdu)

        # Model spectra: [ 'MODEL', 'MODEL_MASK', 'EMLINE', 'STELLAR', 'STELLAR_MASK' ]
        prihdr, modlist = self.model_cubes(prihdr, binned_spectra, stellar_continuum,
                                           stellar_continuum_3d_hdu, emission_line_model,
                                           emission_line_model_3d_hdu)

        # Get the BINIDs
        binidlist = combine_binid_extensions(self.cube, binned_spectra, stellar_continuum, None,
                                             emission_line_model, None, dtype='int32')

        #---------------------------------------------------------------
        # Save the data to the hdu attribute
        prihdr = finalize_dap_primary_header(prihdr, self.cube, self.meta, binned_spectra,
                                             stellar_continuum, loggers=self.loggers,
                                             quiet=self.quiet)

        # Ensure extensions are in the correct order
        self.hdu = fits.HDUList([ fits.PrimaryHDU(header=prihdr),
                                  *binlist,   # FLUX, IVAR, MASK, LSF, WAVE, REDCORR
                                  *modlist,   # MODEL, MODEL_MASK, EMLINE, STELLAR, STELLAR_MASK
                                  *binidlist  # BINID
                                ])
        #---------------------------------------------------------------

        # Check that the path exists
        if not self.directory_path.exists():
            self.directory_path.mkdir(parents=True)

        # Write the model cube file
        DAPFitsUtil.write(self.hdu, str(ofile), overwrite=overwrite, checksum=True,
                          loggers=self.loggers, quiet=self.quiet)
        # End
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, '-'*50)

    def _get_data_mask(self, binned_spectra, binned_spectra_3d_hdu):
        """
        For the binned spectra:
            - consolidate ANY flag except NO_STDDEV from binned_spectra
              into IGNORED; done use SpatiallyBinnedSpectra.do_not_fit_flags
            - propagate NONE_IN_STACK from binned_spectra into
              FLUXINVALID
            - propagate IVARINVALID and FORESTAR
        """
        # Initialize the mask
        mask = numpy.zeros(binned_spectra_3d_hdu['MASK'].data.shape,
                           dtype=self.bitmask.minimum_dtype())

        # Consolidate DIDNOTUSE, FORESTAR, LOW_SPECCOV, LOW_SNR from
        # binned_spectra into IGNORED
        flags = binned_spectra.do_not_fit_flags()
        indx = binned_spectra.bitmask.flagged(binned_spectra_3d_hdu['MASK'].data, flag=flags)
        mask[indx] = self.bitmask.turn_on(mask[indx], 'IGNORED')

        # Propagate NONE_IN_STACK from binned_spectra into FLUXINVALID
        indx = binned_spectra.bitmask.flagged(binned_spectra_3d_hdu['MASK'].data,
                                              flag='NONE_IN_STACK')
        mask[indx] = self.bitmask.turn_on(mask[indx], 'FLUXINVALID')

        # Propagate IVARINVALID
        indx = binned_spectra.bitmask.flagged(binned_spectra_3d_hdu['MASK'].data,
                                              flag='IVARINVALID')
        mask[indx] = self.bitmask.turn_on(mask[indx], 'IVARINVALID')

        # Propagate FORESTAR
        indx = binned_spectra.bitmask.flagged(binned_spectra_3d_hdu['MASK'].data, flag='FORESTAR')
        mask[indx] = self.bitmask.turn_on(mask[indx], 'FORESTAR')

        return mask

    def _get_stellar_continuum_mask(self, stellar_continuum, stellar_continuum_3d_hdu):
        """
        For the stellar continuum models:
            - propagate FORESTAR
            - propagate INVALID_ERROR into IVARINVALID
            - propagate ARTIFACTs
            - consolidate DIDNOTUSE, FORESTAR, LOW_SNR, ARTIFACT,
              OUTSIDE_RANGE, EML_REGION, TPL_PIXELS, TRUNCATED,
              PPXF_REJECT, INVALID_ERROR, NO_FIT into FITIGNORED
            - consolidate FIT_FAILED into FITFAILED
        """
        # Initialize
        mask = numpy.zeros(stellar_continuum_3d_hdu['MASK'].data.shape,
                           dtype=self.bitmask.minimum_dtype())

        # Propagate FORESTAR
        indx = stellar_continuum.bitmask.flagged(stellar_continuum_3d_hdu['MASK'].data,
                                                 flag='FORESTAR')
        mask[indx] = self.bitmask.turn_on(mask[indx], 'FORESTAR')

        # Propagate invalid errors
        indx = stellar_continuum.bitmask.flagged(stellar_continuum_3d_hdu['MASK'].data,
                                                 flag='INVALID_ERROR')
        mask[indx] = self.bitmask.turn_on(mask[indx], 'IVARINVALID')

        # Propagate artifacts
        indx = stellar_continuum.bitmask.flagged(stellar_continuum_3d_hdu['MASK'].data,
                                                 flag='ARTIFACT')
        mask[indx] = self.bitmask.turn_on(mask[indx], 'ARTIFACT')

        # Set which pixel were ignored in the fit
        flags = ['DIDNOTUSE', 'FORESTAR', 'LOW_SNR', 'ARTIFACT', 'OUTSIDE_RANGE', 'EML_REGION',
                 'TPL_PIXELS', 'TRUNCATED', 'PPXF_REJECT', 'INVALID_ERROR', 'NO_FIT' ]
        indx = stellar_continuum.bitmask.flagged(stellar_continuum_3d_hdu['MASK'].data, flag=flags)
        mask[indx] = self.bitmask.turn_on(mask[indx], 'FITIGNORED')

        # Set which pixels failed during the fit
        indx = stellar_continuum.bitmask.flagged(stellar_continuum_3d_hdu['MASK'].data,
                                                 flag='FIT_FAILED')
        mask[indx] = self.bitmask.turn_on(mask[indx], 'FITFAILED')

        return self._include_no_model_mask(mask)

    def _get_emission_line_model_mask(self, emission_line_model, emission_line_model_3d_hdu):
        """
        For the emission-line models:
            - propagate FORESTAR
            - propagate ARTIFACT
            - consolidate DIDNOTUSE, FORESTAR, LOW_SNR, ARTIFACT,
              OUTSIDE_RANGE, TPL_PIXELS, TRUNCATED, PPXF_REJECT,
              INVALID_ERROR into FITIGNORED
        """
        if emission_line_model is None:
            return None

        # Initialize
        mask = numpy.zeros(emission_line_model_3d_hdu['MASK'].data.shape,
                           dtype=self.bitmask.minimum_dtype())

        # Flag FORESTAR
        indx = emission_line_model.bitmask.flagged(emission_line_model_3d_hdu['MASK'].data,
                                                   flag='FORESTAR')
        mask[indx] = self.bitmask.turn_on(mask[indx], 'FORESTAR')

        # Flag artifacts
        indx = emission_line_model.bitmask.flagged(emission_line_model_3d_hdu['MASK'].data,
                                                   flag='ARTIFACT')
        mask[indx] = self.bitmask.turn_on(mask[indx], 'ARTIFACT')

        # Set which pixel were ignored in the fit
        flags = ['DIDNOTUSE', 'FORESTAR', 'LOW_SNR', 'ARTIFACT', 'OUTSIDE_RANGE', 'TPL_PIXELS',
                 'TRUNCATED', 'PPXF_REJECT', 'INVALID_ERROR' ]
        indx = emission_line_model.bitmask.flagged(emission_line_model_3d_hdu['MASK'].data,
                                                   flag=flags)
        mask[indx] = self.bitmask.turn_on(mask[indx], 'FITIGNORED')

        return self._include_no_model_mask(mask)

    def binned_data_cube(self, prihdr, binned_spectra, binned_spectra_3d_hdu):
        """
        Constructs the 'FLUX', 'IVAR', 'MASK', 'LSF', 'WAVE', and
        'REDCORR' model cube extensions, and begins the construction
        of the 'MASK'.

        Returns the primary header, a list of ImageHDUs, and the mask
        array.
        """
        #---------------------------------------------------------------
        ext = ['FLUX', 'IVAR', 'MASK', 'LSF', 'WAVE', 'REDCORR']

        if binned_spectra is None:
            # Construct and return the empty hdus
            return prihdr, DAPFitsUtil.empty_hdus(ext)

        #---------------------------------------------------------------
        # Add data to the primary header
        # TODO: Just copy header from maps file?
        prihdr = binned_spectra.rdxqa._initialize_primary_header(hdr=prihdr)
        prihdr = binned_spectra._initialize_primary_header(hdr=prihdr)
        prihdr = binned_spectra._add_method_header(prihdr)
        prihdr = binned_spectra._add_reddening_header(prihdr)

        # Save the original extension from the DRP file used to
        # construct the datacube spectral resolution data.
        lsfhdr = self.base_cubehdr.copy()
        if hasattr(binned_spectra.cube, 'sres_ext') and binned_spectra.cube.sres_ext is not None:
            lsfhdr['DRPEXT'] = (binned_spectra.cube.sres_ext, 'Source ext from DRP file')

        #---------------------------------------------------------------
        # Get the extension headers; The WAVE and REDCORR extensions
        # have no header.
        hdr = [DAPFitsUtil.finalize_dap_header(self.base_cubehdr, 'FLUX',
                                               bunit='1E-17 erg/s/cm^2/ang/spaxel', err=True,
                                               qual=True),
               DAPFitsUtil.finalize_dap_header(self.base_cubehdr, 'FLUX', hduclas2='ERROR',
                                               bunit='(1E-17 erg/s/cm^2/ang/spaxel)^{-2}',
                                               qual=True, prepend=False),
               DAPFitsUtil.finalize_dap_header(self.base_cubehdr, 'FLUX', hduclas2='QUALITY',
                                               err=True, bit_type=self.bitmask.minimum_dtype(),
                                               prepend=False),
               DAPFitsUtil.finalize_dap_header(lsfhdr, 'LSF', bunit='angstrom'),
               None, None]

        #---------------------------------------------------------------
        # Get the data arrays
        # Reddened flux data
        flux, ivar = binned_spectra.galext.apply(binned_spectra_3d_hdu['FLUX'].data,
                                                 deredden=False,
                                                 ivar=binned_spectra_3d_hdu['IVAR'].data)

        # Convert spectral resolution to dispersion of line-spread
        # function in angstroms
        lsf = numpy.ma.divide(binned_spectra['WAVE'].data[None,None,:],
                              binned_spectra_3d_hdu['SPECRES'].data).filled(0.0) \
                                / DAPConstants.sig2fwhm

        # Return the primary header, the list of ImageHDUs, and the mask
        data = [flux.astype(self.float_dtype), ivar.astype(self.float_dtype),
                self._get_data_mask(binned_spectra, binned_spectra_3d_hdu),
                lsf.astype(self.float_dtype),
                binned_spectra['WAVE'].data.copy().astype(self.float_dtype),
                binned_spectra['REDCORR'].data.copy().astype(self.float_dtype)]
        return prihdr, DAPFitsUtil.list_of_image_hdus(data, hdr, ext)

    def model_cubes(self, prihdr, binned_spectra, stellar_continuum, stellar_continuum_3d_hdu,
                    emission_line_model, emission_line_model_3d_hdu):
        """
        Constructs the 'MODEL', 'MODEL_MASK', 'EMLINE', 'STELLAR', and
        'STELLAR_MASK' model cube extensions, adds the model information
        to the header, and appends data to the mask.
        """
        #---------------------------------------------------------------
        ext = ['MODEL', 'MODEL_MASK', 'EMLINE', 'STELLAR', 'STELLAR_MASK']

        if binned_spectra is None and stellar_continuum is None and emission_line_model is None:
            # Construct and return empty hdus
            return prihdr, DAPFitsUtil.empty_hdus(ext)

        #---------------------------------------------------------------
        # Add the model data to the primary header
        if stellar_continuum is not None:
            prihdr = stellar_continuum._initialize_primary_header(hdr=prihdr)
            prihdr = stellar_continuum._add_method_header(prihdr)
        if emission_line_model is not None:
            prihdr = emission_line_model._initialize_primary_header(hdr=prihdr)
            prihdr = emission_line_model._add_method_header(prihdr)

        #---------------------------------------------------------------
        # Get the extension headers
        hdr = [ DAPFitsUtil.finalize_dap_header(self.base_cubehdr, 'MODEL',
                                                bunit='1E-17 erg/s/cm^2/ang/spaxel', qual=True),
                DAPFitsUtil.finalize_dap_header(self.base_cubehdr, 'MODEL', hduclas2='QUALITY',
                                                bit_type=self.bitmask.minimum_dtype()),
                DAPFitsUtil.finalize_dap_header(self.base_cubehdr, 'EMLINE',
                                                bunit='1E-17 erg/s/cm^2/ang/spaxel'),
                DAPFitsUtil.finalize_dap_header(self.base_cubehdr, 'STELLAR',
                                                bunit='1E-17 erg/s/cm^2/ang/spaxel', qual=True),
                DAPFitsUtil.finalize_dap_header(self.base_cubehdr, 'STELLAR', hduclas2='QUALITY',
                                                bit_type=self.bitmask.minimum_dtype())
              ]

        #---------------------------------------------------------------
        # Get the data arrays
        # Best-fitting composite model
        model = stellar_continuum_3d_hdu['FLUX'].data.copy() if emission_line_model is None \
                    else emission_line_model_3d_hdu['FLUX'].data

        # with reddening
        model = binned_spectra.galext.apply(model, deredden=False)
        if stellar_continuum is None:
            scnt = None
            scnt_mask = None
        else:
            scnt = binned_spectra.galext.apply(stellar_continuum_3d_hdu['FLUX'].data,
                                               deredden=False)
            scnt_mask = self._get_stellar_continuum_mask(stellar_continuum,
                                                         stellar_continuum_3d_hdu)
        if emission_line_model is None:
            line = None
            line_mask = None
        else:
            line = binned_spectra.galext.apply(emission_line_model_3d_hdu['EMLINE'].data,
                                               deredden=False)
            line_mask = self._get_emission_line_model_mask(emission_line_model,
                                                           emission_line_model_3d_hdu)

        # Compile the data arrays
        data = [model.astype(self.float_dtype),
                line_mask, None if line is None else line.astype(self.float_dtype),
                None if scnt is None else scnt.astype(self.float_dtype), scnt_mask]

        # Return the primary header and the list of HDUs
        return prihdr, DAPFitsUtil.list_of_image_hdus(data, hdr, ext)

    def _set_no_model(self, mask):
        indx = numpy.arange(len(mask))[numpy.logical_not(mask > 0)]
        if len(indx) == 0:
            return mask
        mask[:indx[0]] = self.bitmask.turn_on(mask[:indx[0]], 'NOMODEL')
        mask[indx[-1]+1:] = self.bitmask.turn_on(mask[indx[-1]+1:], 'NOMODEL')
        return mask

    def _include_no_model_mask(self, mask):
        """
        Assign the NOMODEL bit to spectral regions outside of the fitted
        spectral range.

        Modeled off of StellarContiuumModel.reset_continuum_mask_window.
        """
        return mask if numpy.sum(mask > 0) == 0 else \
                    numpy.apply_along_axis(self._set_no_model, 2, mask)


#-----------------------------------------------------------------------
def default_paths(cube, file_type='MAPS', rdxqa_method=None, binning_method=None,
                  stelcont_method=None, elmom_method=None, elmodel_method=None, sindex_method=None,
                  method=None, output_path=None, output_file=None):
    """
    Set the default directory and file name for the output file.

    Args:
        cube (:class:`mangadap.datacube.datacube.DataCube`):
            Datacube to analyze.
        file_type (:obj:`str`, optional):
            The DAP output file type.  Should be 'MAPS' or 'LOGCUBE'.
        rdxqa_method (:obj:`str`, optional):
            The method key for the basic assessments of the datacube.
        binning_method (:obj:`str`, optional):
            The method key for the spatial binning.
        stelcont_method (:obj:`str`, optional):
            The method key for the stellar-continuum fitting method.  If
            None, not included in output file name.
        elmom_method (:obj:`str`, optional):
            The method key for the emission-line moments method.  If None,
            not included in the output file name.
        elmodel_method (:obj:`str`, optional):
            The method key for the emission-line modeling method.  If None,
            not included in the output file name.
        sindex_method (:obj:`str`, optional):
            The method key for the spectral index measurements.  If None,
            not included in the output file name.
        method (:obj:`str`, optional):
            The method to use instead of the default.  The default method is
            to simply concatenate the list of methods used by eadch module.
            The order of the keys is the order of operations (rdxqa,
            binning, stellar continuum, emission-line model, spectral
            indices).
        output_path (:obj:`str`, `Path`_, optional):
            The path for the output file.  If None, the current working
            directory is used.
        output_file (:obj:`str`, optional):
            The name of the output file. The full path of the output file
            will be :attr:`directory_path`/:attr:`output_file`.  If None,
            the default is to combine ``cube.output_root`` and the method
            (see ``method``).

    Returns:
        :obj:`tuple`: Returns a `Path`_ with the output directory and a
        :obj:`str` with the output file name.
    """
    directory_path = Path('.').resolve() if output_path is None \
                            else Path(output_path).resolve()
    if method is not None:
        _output_file = f'{cube.output_root}-{file_type}-{method}.fits.gz' \
                            if output_file is None else output_file
        return method, directory_path, _output_file

    # TODO: This is ridiculous...
    method = None
    if rdxqa_method is not None:
        method = f'{rdxqa_method}'
    if binning_method is not None:
        method = f'{binning_method}' if method is None else f'{method}-{binning_method}'
    if stelcont_method is not None:
        method = f'{stelcont_method}' if method is None else f'{method}-{stelcont_method}'
    if elmom_method is not None and elmodel_method is not None:
        method = f'{elmodel_method}-{elmom_method}' if method is None \
                    else f'{method}-{elmodel_method}-{elmom_method}'
    elif elmom_method is None and elmodel_method is not None:
        method = f'{elmodel_method}' if method is None else f'{method}-{elmodel_method}'
    elif elmom_method is not None and elmodel_method is None:
        method = f'{elmom_method}' if method is None else f'{method}-{elmom_method}'
    if sindex_method is not None:
        method = f'{sindex_method}' if method is None else f'{method}-{sindex_method}'
    _output_file = f'{cube.output_root}-{file_type}-{method}.fits.gz' \
                        if output_file is None else output_file
    return method, directory_path, _output_file


def combine_binid_extensions(cube, binned_spectra, stellar_continuum, emission_line_moments,
                             emission_line_model, spectral_indices, dtype=None):

    """
    Combine the bin IDs from the different analysis steps into a single
    multichannel extension.
    """
    if cube is None:
        raise ValueError('DataCube must be provided.')

    hdr = DAPFitsUtil.finalize_dap_header( DAPFitsUtil.build_map_header(cube.fluxhdr,
                                        'K Westfall <westfall@ucolick.org> & SDSS-IV Data Group',
                                                    multichannel=True, maskname='MANGA_DAPPIXMASK'),
                                          'BINID', multichannel=True,
                                          channel_names=[ 'Binned spectra', 'Stellar continua',
                                                          'Em. line moments', 'Em. line models',
                                                          'Spectral indices' ])

    # Build the ID data for each analysis product
    binid = numpy.zeros(cube.spatial_shape+(5,), dtype=int)
    if binned_spectra is not None:
        binid[:,:,0] = binned_spectra['BINID'].data
    if stellar_continuum is not None:
        binid[:,:,1] = stellar_continuum['BINID'].data
    if emission_line_moments is not None:
        binid[:,:,2] = emission_line_moments['BINID'].data
    if emission_line_model is not None:
        binid[:,:,3] = emission_line_model['BINID'].data
    if spectral_indices is not None:
        binid[:,:,4] = spectral_indices['BINID'].data

    _dtype = 'int' if dtype is None else dtype

    return DAPFitsUtil.list_of_image_hdus([ binid.astype(_dtype) ], [ hdr ], [ 'BINID' ])


def confirm_dap_types(cube, rdxqa, binned_spectra, stellar_continuum, emission_line_moments,
                      emission_line_model, spectral_indices):

    if not isinstance(cube, DataCube):
        raise TypeError('Input must have type DataCube.')
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


# TODO: Need to abstract this for non-DRP datacubes.
def finalize_dap_primary_header(prihdr, cube, metadata, binned_spectra, stellar_continuum,
                                loggers=None, quiet=False):

    # Initialize the DAP quality flag
    dapqualbm = DAPQualityBitMask()
    dapqual = dapqualbm.minimum_dtype()(0)          # Type casting original flag to 0
    # TODO: Move this to the datacube base class?
    if cube.meta['drpcrit']:
        if not quiet:
            log_output(loggers, 1, logging.INFO, 'DRP file is flagged CRITICAL!')
        dapqual = dapqualbm.turn_on(dapqual, ['CRITICAL', 'DRPCRIT'])

    # Flag the file as CRITICAL if the stellar continuum fits are bad
    # for all spectra
    if stellar_continuum is not None:
        mask = stellar_continuum.bitmask.flagged(stellar_continuum['PAR'].data['MASK'],
                                                 ['NO_FIT', 'FIT_FAILED', 'INSUFFICIENT_DATA',
                                                  'NEAR_BOUND'])
        if numpy.all(mask):
            dapqual = dapqualbm.turn_on(dapqual, ['CRITICAL', 'DAPCRIT'])

    # TODO: Do similar checks for the other major modules

    # Signify that the Voronoi binning resulted in a single bin
    if binned_spectra is not None and binned_spectra.method['binclass'] is not None \
            and binned_spectra.method['binclass'].bintype == 'voronoi' \
            and binned_spectra.nbins == 1:
        dapqual = dapqualbm.turn_on(dapqual, 'SINGLEBIN')

    # Input photometric geometry and scale were invalid
    if 'ell' not in metadata.keys() or metadata['ell'] is None or \
            'pa' not in metadata.keys() or metadata['pa'] is None or \
            'reff' not in metadata.keys() or metadata['reff'] is None:
        dapqual = dapqualbm.turn_on(dapqual, 'BADGEOM')

    # Propagate some flags
    # TODO: This is a bit of a hack.  I should be more careful about how
    # "propagated" flags are defined.

    # Determine if there's a foreground star
    if cube.propagate_flags() is not None:
        flags = cube.propagate_flags()
        # TODO: This will barf if flags is a numpy array
        flags = flags if isinstance(flags, list) else [flags]
        for flag in flags:
            if numpy.any(cube.bitmask.flagged(cube.mask, flag=flag)):
                if flag in dapqualbm.keys():
                    dapqual = dapqualbm.turn_on(dapqual, flag)
                warnings.warn(f'{flag} is not a flag in {dapqualbm.__class__.__name__} and '
                              'cannot be propagated!')

    # Commit the quality flag to the header
    if dapqualbm.flagged(dapqual, 'CRITICAL'):
        warnings.warn('DAP files are flagged CRITICAL.')
    prihdr['DAPQUAL'] = (dapqual, 'DAP quality bitmask')

    # Finalize authors
    prihdr['AUTHOR'] = 'K Westfall <westfall@ucolick.org>'
    return prihdr


def add_snr_metrics_to_header(hdr, cube, r_re):
    """
    For all valid spaxels within 1 Re < R < 1.5 Re calculate the median
    S/N and combined S/N (including covariance) in the griz bands.

    Adds the following keywords (8) to the header:
        - SNR?MED: Median S/N within 1-1.5 Re in each of the griz bands.
        - SNR?RING: The combined S/N of all spaxels within 1-1.5 Re in
          the griz bands, incuding covariance

    The inclusion of covariance is based on the exact calculation of the
    correlation matrix at the flux-weighted center of each band, and
    assuming no change in the correlation between the response-weighted
    correlation within the broad band.

    Args:
        hdr (`astropy.io.fits.Header`_):
            Header object to edit
        cube (:class:`mangadap.datacube.datacube.DataCube`):
            Datacube
        r_re (`numpy.ndarray`_):
            *Flattened* array with the semi-major axis
            radius in units of the effective radius for all spaxels
            in the datacube. Shape should be (Nx*Ny,), where the
            shape of the datacube is (Nx,Ny,Nwave).

    Returns:
        `astropy.io.fits.Header`_: The edited header object.

    Raises:
        FileNotFoundError: Raised if any of the response function files
            cannot be found.
    """
    root = defaults.dap_data_root()
    filter_response_file = [root / 'filter_response' / f 
                            for f in ['gunn_2001_g_response.db', 'gunn_2001_r_response.db',
                                      'gunn_2001_i_response.db', 'gunn_2001_z_response.db']]
    for f in filter_response_file:
        if not f.exists():
            raise FileNotFoundError('{0} does not exist!'.format(f))

    # Set the header keywords
    key_med = ['SNRGMED', 'SNRRMED', 'SNRIMED', 'SNRZMED' ]
    com_med = [ 'Median g-band SNR from 1-1.5 Re', 'Median r-band SNR from 1-1.5 Re',
                'Median i-band SNR from 1-1.5 Re', 'Median z-band SNR from 1-1.5 Re' ]
                    
    key_ring = ['SNRGRING', 'SNRRRING', 'SNRIRING', 'SNRZRING' ]
    com_ring = [ 'g-band SNR in 1-1.5 Re bin', 'r-band SNR in 1-1.5 Re bin',
                 'i-band SNR in 1-1.5 Re bin', 'z-band SNR in 1-1.5 Re bin' ]

    # Run the S/N calculation; this is the same calculation as done in
    # mangadap.proc.spatialbinning.VoronoiBinning.sn_calculation_covariance_matrix
    nfilter = len(filter_response_file)
    flags = cube.do_not_use_flags()
#    flags = ['DONOTUSE', 'FORESTAR']
    for i in range(nfilter):
        response_func = numpy.genfromtxt(str(filter_response_file[i]))[:,:2]
        signal, variance, snr = [a.ravel() for a in 
                                 cube.flux_stats(response_func=response_func, flag=flags)]
        covar = None
        if cube.can_compute_covariance:
            covar_wave = cube.central_wavelength(response_func=response_func, flag=flags)
            covar_channel = numpy.argsort(numpy.absolute(cube.wave-covar_wave))[0]
            covar = cube.covariance_matrix(covar_channel)

        # Get the spaxels within the radius limits
        indx = numpy.arange(r_re.size)[(r_re > 1) & (r_re < 1.5) & numpy.logical_not(signal.mask)]
        if len(indx) == 0:
            hdr[key_med[i]] = (0., com_med[i])
            hdr[key_ring[i]] = (0., com_ring[i])
            continue

        if covar is None:
            warnings.warn('Covariance not available!  Continuing without it.')
            covar = Covariance.from_variance(variance, correlation=True)
            
        # Get the appropriate covariance pixels to select
        ci, cj = map( lambda x: x.ravel(), numpy.meshgrid(indx, indx) )

        # Use the covariance matrix from the single wavelength channel
        # calculation, but renormalize it to the mean variance over the
        # response function
        covar = numpy.identity(len(variance), dtype=float) if covar is None \
                        else covar.apply_new_variance(variance).toarray()
        # Set the median S/N ...
        hdr[key_med[i]] = (numpy.ma.median(snr[indx]), com_med[i])
        # ... and the combined S/N
        hdr[key_ring[i]] = (numpy.ma.sum(signal[indx])/numpy.sqrt(numpy.sum(covar[ci,cj])),
                            com_ring[i])

    return hdr






