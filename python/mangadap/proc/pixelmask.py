# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
A class heirarchy for pixel masks.

*License*:
    Copyright (c) 2015, SDSS-IV/MaNGA Pipeline Group
        Licensed under BSD 3-clause license - see LICENSE.rst

*Source location*:
    $MANGADAP_DIR/python/mangadap/proc/pixelmask.py

*Imports and python version compliance*:
    ::

        from __future__ import division
        from __future__ import print_function
        from __future__ import absolute_import
        from __future__ import unicode_literals

        import sys
        import warnings
        if sys.version > '3':
            long = int

        import numpy
        import astropy.constants

        from ..util.bitmask import BitMask
        from .artifactdb import ArtifactDB
        from .emissionlinedb import EmissionLineDB

*Class usage examples*:
        Add examples

*Revision history*:
    | **18 Apr 2016**: Original implementation K. Westfall (KBW)

"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import sys
import warnings
if sys.version > '3':
    long = int

import numpy
import astropy.constants

from ..util.bitmask import BitMask
from .artifactdb import ArtifactDB
from .emissionlinedb import EmissionLineDB

from matplotlib import pyplot

__author__ = 'Kyle B. Westfall'
# Add strict versioning
# from distutils.version import StrictVersion


class PixelMask:
    def __init__(self):
        self.shape = None


    def _empty_mask(self, wave, nspec=None):
        self.shape = (len(wave),) if nspec is None else (nspec,len(wave))
        return numpy.full(self.shape, False, dtype=numpy.bool)


    def _mask_wavelength_ranges(self, wave, waverange, nspec=None):
        """
        Mask a set of wavelength bands
        """
        self.shape = (len(wave),) if nspec is None else (nspec,len(wave))
        mask = numpy.any( numpy.array([ numpy.logical_and(wave>l, wave<u) \
                                        for l,u in zip(waverange[:,0], waverange[:,1])]), axis=0 )
        return mask if nspec is None else numpy.array([mask]*nspec)
        

class SpectralPixelMask(PixelMask):
    """

    Container that produces a mask for the stellar continuum based on a
    set of emission lines and artifacts.

    """
    def __init__(self, artdb=None, emldb=None, waverange=None):
        if artdb is not None and not isinstance(artdb, ArtifactDB):
            raise TypeError('Must provide EmissionLineDB for emission-lines to mask.')
        self.artdb = artdb

        if emldb is not None and not isinstance(emldb, EmissionLineDB):
            raise TypeError('Must provide EmissionLineDB for emission-lines to mask.')
        self.emldb = emldb

        if waverange is not None and len(waverange) != 2:
            raise ValueError('Provided wavelength range must have two and only two elements.')
        self.waverange = waverange


    def _waverange_mask(self, wave, nspec=None):
        if self.waverange is None:
            return self._empty_mask(wave, nspec=nspec)
        return self._mask_wavelength_ranges(wave, numpy.array([waverange]), nspec=nspec)


    def _artifact_mask(self, wave, nspec=None):
        if self.artdb is None:
            return self._empty_mask(wave, nspec=nspec)
        return self._mask_wavelength_ranges(wave, self.artdb['waverange'], nspec=nspec)


    def _check_eml_kin_argument(self, kin):
        if kin is None:
            return kin

        if isinstance(kin, float):
            return numpy.full(self.emldb.neml, kin, dtype=numpy.float)
        if isinstance(kin, (list, numpy.ndarray)):
            if len(kin) != self.emldb.neml:
                raise ValueError('Provided vector does not have a matching length.')
            return numpy.asarray(kin).astype(numpy.float)
       

    def _get_emission_line_bands(self, velocity_offsets=0.0, sigma=250.0, nsigma=3.0):

        # Mask everything but the lines to ignore
        indx = numpy.invert(self.emldb['action'] == 'i')
        nbands = numpy.sum(indx)
        # No lines to mask
        if nbands == 0:
            return None

        # Get the number of standard deviations to cover with the mask
        _nsigma = self._check_eml_kin_argument(nsigma)
        if _nsigma is None:
            raise ValueError('Must provide the number of sigma to cover with the mask.')

        # Get the mask centers
        _velocity_offsets = self._check_eml_kin_argument(velocity_offsets)
        _velocity_offsets[ self.emldb['action'] == 's' ] = 0.0
        center = self.emldb['restwave'][indx] if _velocity_offsets is None else \
                    self.emldb['restwave'][indx] \
                        * (1.0 + _velocity_offsets[indx]/astropy.constants.c.to('km/s').value)
        
        # Get the mask widths
        _sigma = self._check_eml_kin_argument(sigma)
        halfwidth = _nsigma[indx] * (self.emldb['sig'][indx] if _sigma is None else _sigma[indx])

        return numpy.array([ center*(1.0-halfwidth/astropy.constants.c.to('km/s').value),
                             center*(1.0+halfwidth/astropy.constants.c.to('km/s').value) ]).T


    def _emission_line_mask(self, wave, nspec=None, velocity_offsets=0.0, sigma=250.0, nsigma=3.0):
        """
        Velocities are spectrum-dependent.  Sigma and nsigma are line-dependent.
        """
        if self.emldb is None:
            return self._empty_mask(wave, nspec=nspec)


        if (isinstance(velocity_offsets, (list, numpy.ndarray)) and len(velocity_offsets) > 1): 
            _velocity_offsets = numpy.asarray(velocity_offsets)
            if len(_velocity_offsets) != nspec:
                raise ValueError('Velocity offsets do not match the number of spectra.')
             
            mask = self._empty_mask(wave, nspec=nspec)
            for i in range(len(_velocity_offsets)):
                waverange = self._get_emission_line_bands(velocity_offsets=_velocity_offsets[i],
                                                          sigma=sigma, nsigma=nsigma)
                mask[i,:] = self._mask_wavelength_ranges(wave, waverange)
            return mask
        
        waverange = self._get_emission_line_bands(velocity_offsets=velocity_offsets, sigma=sigma,
                                                  nsigma=nsigma)
        return self._mask_wavelength_ranges(wave, waverange, nspec=nspec)
        


    def boolean(self, wave, nspec=None, velocity_offsets=0.0, sigma=250.0, nsigma=3.0):
        if nspec is not None and not nspec > 0:
            raise ValueError('Number of spectra must be larger than 0!')
        return self._waverange_mask(wave, nspec=nspec) | self._artifact_mask(wave, nspec=nspec) \
                | self._emission_line_mask(wave, nspec=nspec, velocity_offsets=velocity_offsets,
                                           sigma=sigma, nsigma=nsigma)
        

    def bits(self, bitmask, wave, nspec=None, mask=None, velocity_offsets=0.0, sigma=250.0,
             nsigma=3.0, waverange_flag='OUTSIDE_RANGE', art_flag='ARTIFACT',
             eml_flag='EML_REGION'):

        # Check the bitmask type
        if not isinstance(bitmask, BitMask):
            raise TypeError('Must provide object of type BitMask.')
        if nspec is not None and not nspec > 0:
            raise ValueError('Number of spectra must be larger than 0!')

        # Get the wavelength range mask
        wavemask = self._waverange_mask(wave, nspec=nspec)
        # Check that the input mask has the same size
        if mask is not None and wavemask.shape != mask.shape:
            raise ValueError('Input mask does not have the correct shape.')

        # Get the artifact mask
        artmask = self._artifact_mask(wave, nspec=nspec)

        # Get the emission-line mask
        emlmask = self._emission_line_mask(wave, nspec=nspec, velocity_offsets=velocity_offsets,
                                           sigma=sigma, nsigma=nsigma)

        # Construct and return the mask
        _mask = numpy.zeros(wavemask.shape, dtype=bitmask.minimum_dtype()) \
                    if mask is None else mask

        _mask[wavemask] = bitmask.turn_on(_mask[wavemask], flag=waverange_flag)
        _mask[artmask] = bitmask.turn_on(_mask[artmask], flag=art_flag)
        _mask[emlmask] = bitmask.turn_on(_mask[emlmask], flag=eml_flag)
        return _mask


