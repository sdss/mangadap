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

        from .bitmask import BitMask
        from ..par.artifactdb import ArtifactDB
        from ..par.emissionlinedb import EmissionLineDB

*Class usage examples*:
        Add examples

*Revision history*:
    | **18 Apr 2016**: Original implementation K. Westfall (KBW)
    | **30 Nov 2016**: (KBW) Generalization of :class:`PixelMask` and
        minor edits

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

from .bitmask import BitMask
from ..par.artifactdb import ArtifactDB
from ..par.emissionlinedb import EmissionLineDB

from matplotlib import pyplot

__author__ = 'Kyle B. Westfall'
# Add strict versioning
# from distutils.version import StrictVersion


class PixelMask:
    """
    Base class for a general 1D or 2D pixel mask.

    Attributes:
        shape (tuple): Shape of the array for the mask.
    """
    def __init__(self):
        self.shape = None


    def _set_shape(self, x, ny=None):
        """
        Set :attr:`shape` of the object

        Args:
            x (numpy.ndarray): Vector with the x-coordinates
            ny (int): (**Optional**) Size of the second dimension.
                Default is that there is no second dimension.
        """
        self.shape = (len(x),) if ny is None else (ny,len(x))


    def _empty_mask(self, x, ny=None):
        """
        Return an empty mask with the correct shape.

        Args:
            x (numpy.ndarray): Coordinate vector
            ny (int): (**Optional**) Size of the second dimension.
                Default is that there is only one dimension.

        Returns:
            numpy.ndarray : Boolean mask of the correct shape.
        """
        self._set_shape(x, ny=ny)
        return numpy.full(self.shape, False, dtype=numpy.bool)


    def _mask_coordinate_ranges(self, x, rng, ny=None):
        """
        Flag any x coordinates between a set of range limits as True.
        The mask is repeated in the second dimension, if requested.

        Args:
            x (numpy.ndarray): Coordinate vector
            rng (list, numpy.ndarray): (List of) Coordinate ranges that
                should be masked.
            ny (int): (**Optional**) Size of the second dimension.
                Default is that there is only one dimension.

        Returns:
            numpy.ndarray : Boolean mask of the correct shape.
        """
        self._set_shape(x, ny=ny)
        _rng = numpy.atleast_2d(rng)
        mask = numpy.any( numpy.array([ numpy.logical_and(x>l, x<u) \
                                        for l,u in zip(_rng[:,0], _rng[:,1])]), axis=0)
        return mask if len(self.shape) == 1 else numpy.array([mask]*self.shape[0])
        

class SpectralPixelMask(PixelMask):
    """
    Container that produces a mask for the stellar continuum based on a
    set of emission lines and artifacts.

    Args:
        artdb (:class:`mangadap.proc.artifactdb.ArtifactDB`):
            (**Optional**) Database with the list of artifacts to mask.
        emldb (:class:`mangadap.proc.emissionlinedb.EmissionLineDB`):
            (**Optional**) Database with the list of emission lines to
            mask.
        waverange (numpy.ndarray): (**Optional**) Any pixels **outside**
            this wavelength range are masked.

    Attributes:
        artdb (:class:`mangadap.proc.artifactdb.ArtifactDB`):
            Database with the list of artifacts to mask.
        emldb (:class:`mangadap.proc.emissionlinedb.EmissionLineDB`):
            Database with the list of emission lines to
            mask.
        waverange (numpy.ndarray): Any pixels **outside**
            this wavelength range are masked.
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
        """
        Mask the pixels **not** within the selected wavelength range.
        
        Args:
            wave (numpy.ndarray): Wavelength coordinate vector
            nspec (int): (**Optional**) Number of spectra to mask.
                Default is just one.

        Returns:
            numpy.ndarray : Boolean mask of the correct shape.
        """
        if self.waverange is None:
            return self._empty_mask(wave, ny=nspec)
        return numpy.invert(self._mask_coordinate_ranges(wave, self.waverange, ny=nspec))


    def _artifact_mask(self, wave, nspec=None):
        """
        Mask the pixels in the wavelength range(s) defined by the
        artifact database.
        
        Args:
            wave (numpy.ndarray): Wavelength coordinate vector
            nspec (int): (**Optional**) Number of spectra to mask.
                Default is just one.

        Returns:
            numpy.ndarray : Boolean mask of the correct shape.
        """
        if self.artdb is None:
            return self._empty_mask(wave, ny=nspec)
        return self._mask_coordinate_ranges(wave, self.artdb['waverange'], ny=nspec)


    def _check_eml_kin_argument(self, kin):
        """
        Check and return the correct kinematics vectors.  If a single
        number is provided, the function returns the number repeated for
        each emission line.

        Args:
            kin (float, list, numpy.ndarray): An input set of kinematics
                to used by the mask.

        Returns:
            numpy.ndarray : A 1D float array of the correct shape with
                the kinematics 

        """
        if kin is None:
            return None

        if isinstance(kin, (int,float)):
            return numpy.full(self.emldb.neml, kin, dtype=numpy.float)
        if isinstance(kin, (list, numpy.ndarray)):
            if len(kin) != self.emldb.neml:
                raise ValueError('Provided vector does not have a matching length.')
            return numpy.atleast_1d(kin).astype(numpy.float)
       

    def _get_emission_line_bands(self, velocity_offsets=0.0, sigma=250.0, nsigma=3.0):
        r"""

        Set the emission-line masks based on the redshift and dispersion
        of the line.

        Args:
            velocity_offsets (float, numpy.ndarray): (**Optional**) One
                or more velocity offsets to apply to the emission-line
                bands. Default is to apply no velocity offset.
            sigma (float, numpy.ndarray): (**Optional**) One or more
                velocity dispersions to use for setting the width of the
                emission-line band.  Default is a width based on a
                dispersion of 250 km/s.
            nsigma (float, numpy.ndarray): (**Optional**) One or more
                numbers that sets the width of the band in units of the
                provided velocity dipsersions.

        Returns:
            numpy.ndarray : A :math:`N_l \times 2` array with the
                wavelength limits of the emission-line bands.

        """

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
        Mask the pixels in the wavelength range(s) defined by the
        emission-line database that has been adjusted by a set of
        velocities and dispersions.
        
        Currently, the velocity offsets are applied to all lines for
        each spectrum, whereas sigma and nsigma are applied to all
        spectra for each line.

        Args:
            wave (numpy.ndarray): Wavelength coordinate vector
            nspec (int): (**Optional**) Number of spectra to mask.
                Default is just one.
            velocity_offsets (float, numpy.ndarray): (**Optional**) One
                or more velocity offsets to apply to the emission-line
                bands on a spectrum-by-spectrum basis. Default is to
                apply no velocity offset.
            sigma (float, numpy.ndarray): (**Optional**) One or more
                velocity dispersions to use for setting the width of the
                emission-line band on a line-by-line basis.  Default is
                a width based on a dispersion of 250 km/s.
            nsigma (float, numpy.ndarray): (**Optional**) One or more
                numbers that sets the width of the band in units of the
                provided velocity dipsersions band on a line-by-line
                basis.

        Returns:
            numpy.ndarray : Boolean mask of the correct shape.
        """
        if self.emldb is None:
            return self._empty_mask(wave, ny=nspec)


        if isinstance(velocity_offsets, (list, numpy.ndarray)) and len(velocity_offsets) > 1: 
            _velocity_offsets = numpy.asarray(velocity_offsets)
            if len(_velocity_offsets) != nspec:
                raise ValueError('Velocity offsets do not match the number of spectra.')
             
            mask = self._empty_mask(wave, ny=nspec)
            for i in range(len(_velocity_offsets)):
                waverange = self._get_emission_line_bands(velocity_offsets=_velocity_offsets[i],
                                                          sigma=sigma, nsigma=nsigma)
                mask[i,:] = self._mask_coordinate_ranges(wave, waverange)
            return mask

        _velocity_offsets = velocity_offsets[0] \
                if isinstance(velocity_offsets, (list, numpy.ndarray)) \
                        and len(velocity_offsets) == 1 else velocity_offsets
        
        waverange = self._get_emission_line_bands(velocity_offsets=_velocity_offsets, sigma=sigma,
                                                  nsigma=nsigma)
        return self._mask_coordinate_ranges(wave, waverange, ny=nspec)
        


    def boolean(self, wave, nspec=None, velocity_offsets=0.0, sigma=250.0, nsigma=3.0):
        """
        Construct the full boolean mask that includes the desired
        wavelength range, omitting artifacts, and omitting emission
        lines.

        Args:
            wave (numpy.ndarray): Wavelength coordinate vector
            nspec (int): (**Optional**) Number of spectra to mask.
                Default is just one.
            velocity_offsets (float, numpy.ndarray): (**Optional**) One
                or more velocity offsets to apply to the emission-line
                bands on a spectrum-by-spectrum basis. Default is to
                apply no velocity offset.
            sigma (float, numpy.ndarray): (**Optional**) One or more
                velocity dispersions to use for setting the width of the
                emission-line band on a line-by-line basis.  Default is
                a width based on a dispersion of 250 km/s.
            nsigma (float, numpy.ndarray): (**Optional**) One or more
                numbers that sets the width of the band in units of the
                provided velocity dipsersions band on a line-by-line
                basis.

        Returns:
            numpy.ndarray : Boolean mask of the correct shape.
        """
        if nspec is not None and not nspec > 0:
            raise ValueError('Number of spectra must be larger than 0!')
        return self._waverange_mask(wave, nspec=nspec) | self._artifact_mask(wave, nspec=nspec) \
                | self._emission_line_mask(wave, nspec=nspec, velocity_offsets=velocity_offsets,
                                           sigma=sigma, nsigma=nsigma)
        

    def bits(self, bitmask, wave, nspec=None, mask=None, velocity_offsets=0.0, sigma=250.0,
             nsigma=3.0, waverange_flag='OUTSIDE_RANGE', art_flag='ARTIFACT',
             eml_flag='EML_REGION'):
        """
        Construct a bit mask that signifies pixels as outside the
        desired wavelength range, as being affected by an artifact, and
        being designated as an emission-line region.

        Args:
            wave (numpy.ndarray): Wavelength coordinate vector
            nspec (int): (**Optional**) Number of spectra to mask.
                Default is just one.
            velocity_offsets (float, numpy.ndarray): (**Optional**) One
                or more velocity offsets to apply to the emission-line
                bands on a spectrum-by-spectrum basis. Default is to
                apply no velocity offset.
            sigma (float, numpy.ndarray): (**Optional**) One or more
                velocity dispersions to use for setting the width of the
                emission-line band on a line-by-line basis.  Default is
                a width based on a dispersion of 250 km/s.
            nsigma (float, numpy.ndarray): (**Optional**) One or more
                numbers that sets the width of the band in units of the
                provided velocity dipsersions band on a line-by-line
                basis.

        Returns:
            numpy.ndarray : Bit mask of the correct shape.
        """

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


