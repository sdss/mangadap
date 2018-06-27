# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
A class hierarchy that performs the spectral-index measurements.

*License*:
    Copyright (c) 2015, SDSS-IV/MaNGA Pipeline Group
        Licensed under BSD 3-clause license - see LICENSE.rst

*Source location*:
    $MANGADAP_DIR/python/mangadap/proc/spectralindices.py

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

        import glob
        import os.path
        import numpy
        from astropy.io import fits

        from ..par.parset import ParSet
        from ..par.artifactdb import ArtifactDB
        from ..par.absorptionindexdb import AbsorptionIndexDB
        from ..par.bandheadindexdb import BandheadIndexDB
        from ..config.defaults import dap_source_dir, default_dap_file_name
        from ..config.defaults import default_dap_method, default_dap_method_path
        from ..config.defaults import default_dap_common_path
        from ..util.instrument import SpectralResolution, match_spectral_resolution
        from ..util.instrument import spectral_coordinate_step, spectrum_velocity_scale
        from ..util.fitsutil import DAPFitsUtil
        from ..util.fileio import init_record_array, rec_to_fits_type
        from ..util.log import log_output
        from ..util.bitmask import BitMask
        from ..util.pixelmask import SpectralPixelMask
        from ..util.parser import DefaultConfig
        from .spatiallybinnedspectra import SpatiallyBinnedSpectra
        from .templatelibrary import TemplateLibrary
        from .stellarcontinuummodel import StellarContinuumModel
        from .emissionlinemodel import EmissionLineModel
        from .bandpassfilter import passband_integral, passband_integrated_width
        from .bandpassfilter import passband_integrated_mean
        from .bandpassfilter import passband_weighted_mean
        from .util import select_proc_method, flux_to_fnu


*Class usage examples*:
    Add examples!


*Notes*:
    
    If neither stellar-continuum nor emission-line models are provided:
        - Indices are measure on the binned spectra
        - No velocity-dispersion corrections are calculated

    If a stellar-continuum model is provided without an emission-line
    model:
        - Indices are measured on the binned spectra
        - Velocity-dispersion corrections are computed for any binned
          spectrum with a stellar-continuum fit based on the optimal
          template

    If an emission-line model is provided without a stellar-continuum
    model:
        - Indices are measured on the relevant (binned or unbinned)
          spectra; spectra with emission-line fits have the model
          emission lines subtracted from them before these measurements.
        - If the emission-line model includes data regarding the
          stellar-continuum fit (template spectra and template weights),
          corrections are calculated for spectra with emission-line
          models based on the continuum fits; otherwise, no corrections
          are calculated.

    If both stellar-continuum and emission-line models are provided, and
    if the stellar-continuum and emission-line fits are performed on the
    same spectra:
        - Indices are measured on the relevant (binned or unbinned)
          spectra; spectra with emission-line fits have the model
          emission lines subtracted from them before these measurements.
        - Velocity-dispersion corrections are based on the
          stellar-continuum templates and weights

    If both stellar-continuum and emission-line models are provided, and
    if the stellar-continuum and emission-line fits are performed on
    different spectra:
        - The behavior is exactly as if the stellar-continuum model was
          not provided.

*Revision history*:
    | **20 Apr 2016**: Implementation begun by K. Westfall (KBW)
    | **09 May 2016**: (KBW) Add subtraction of emission-line models
    | **11 Jul 2016**: (KBW) Allow to not apply dispersion corrections
        for index measurements
    | **28 Jul 2016**: (KBW) Fixed error in initialization of guess
        redshift when stellar continuum is provided.
    | **23 Feb 2017**: (KBW) Use DAPFitsUtil read and write functions.
    | **27 Feb 2017**: (KBW) Use DefaultConfig
    | **02 Feb 2018**: (KBW) Allow for stellar-continuum and
        emission-line models to be performed on different spectra (i.e.,
        allow for the hybrid binning scheme).  Adjust for change to
        :func:`mangadap.proc.stellarcontinuummodel.StellarContinuumModel.fill_to_match`.
    | **15 Mar 2018**: (KBW) Correct the indices measured in angstroms
        for redshift.  Keep the indices as measured by the best-fitting
        model.

.. _astropy.io.fits.hdu.hdulist.HDUList: http://docs.astropy.org/en/v1.0.2/io/fits/api/hdulists.html
.. _glob.glob: https://docs.python.org/3.4/library/glob.html
.. _numpy.recarray: https://docs.scipy.org/doc/numpy/reference/generated/numpy.recarray.html


"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import sys
import warnings
if sys.version > '3':
    long = int

import glob
import os
import logging

import numpy
from astropy.io import fits

from ..par.parset import ParSet
from ..par.artifactdb import ArtifactDB
from ..par.absorptionindexdb import AbsorptionIndexDB
from ..par.bandheadindexdb import BandheadIndexDB
from ..config.defaults import dap_source_dir, default_dap_file_name
from ..config.defaults import default_dap_method, default_dap_method_path
from ..config.defaults import default_dap_common_path
from ..util.instrument import SpectralResolution, match_spectral_resolution
from ..util.instrument import spectral_coordinate_step, spectrum_velocity_scale
from ..util.fitsutil import DAPFitsUtil
from ..util.fileio import init_record_array, rec_to_fits_type
from ..util.log import log_output
from ..util.bitmask import BitMask
from ..util.pixelmask import SpectralPixelMask
from ..util.parser import DefaultConfig
from .spatiallybinnedspectra import SpatiallyBinnedSpectra
from .templatelibrary import TemplateLibrary
from .stellarcontinuummodel import StellarContinuumModel
from .emissionlinemodel import EmissionLineModel
from .bandpassfilter import passband_integral, passband_integrated_width
from .bandpassfilter import passband_integrated_mean
from .bandpassfilter import passband_weighted_mean, pseudocontinuum
from .util import select_proc_method, flux_to_fnu

import astropy.constants
from matplotlib import pyplot
#from memory_profiler import profile

class SpectralIndicesDef(ParSet):
    """
    A class that holds the parameters necessary to perform the
    spectral-index measurements.

    Args:
        key (str): Keyword used to distinguish between different
            spectral-index databases.
        minimum_snr (bool): Minimum S/N of spectrum to fit
        fwhm (int, float): Resolution FWHM in angstroms at which to make
            the measurements.
        compute_corrections (bool): Flag to compute dispersion corrections
            to indices.  Dispersion corrections are always calculated!
        artifacts (str): String identifying the artifact database to use
        absindex (str): String identifying the absorption-index database
            to use
        bandhead (str): String identifying the bandhead-index database
            to use
    """
    def __init__(self, key, minimum_snr, fwhm, compute_corrections, artifacts, absindex, bandhead):
        in_fl = [ int, float ]

        pars =     [ 'key', 'minimum_snr', 'fwhm', 'compute_corrections', 'artifacts', 'absindex',
                        'bandhead' ]
        values =   [   key,   minimum_snr,   fwhm,   compute_corrections,   artifacts,   absindex,
                          bandhead ]
        dtypes =   [   str,         in_fl,  in_fl,                  bool,         str,        str,
                               str ]

        ParSet.__init__(self, pars, values=values, dtypes=dtypes)


def validate_spectral_indices_config(cnfg):
    """ 
    Validate :class:`mangadap.util.parser.DefaultConfig` object with
    spectral-index measurement parameters.

    Args:
        cnfg (:class:`mangadap.util.parser.DefaultConfig`): Object with
            parameters to validate.

    Raises:
        KeyError: Raised if any required keywords do not exist.
        ValueError: Raised if keys have unacceptable values.
        FileNotFoundError: Raised if a file is specified but could not
            be found.
    """
    # Check for required keywords
    if not cnfg.keyword_specified('key'):
        raise KeyError('Keyword \'key\' must be provided.')

    if not cnfg.keyword_specified('absorption_indices') \
                and not cnfg.keyword_specified('bandhead_indices'):
        raise ValueError('Must provide either an absorption-index database or a bandhead-index ' \
                         'database, or both.')


def available_spectral_index_databases(dapsrc=None):
    """
    Return the list of available spectral index databases

    Available database combinations:

    .. todo::
        Fill in

    Args:
        dapsrc (str): (**Optional**) Root path to the DAP source
            directory.  If not provided, the default is defined by
            :func:`mangadap.config.defaults.dap_source_dir`.

    Returns:
        list: A list of :func:`SpectralIndicesDef` objects, each
        defining a spectral-index database to measure.

    Raises:
        NotADirectoryError: Raised if the provided or default
            *dapsrc* is not a directory.
        OSError/IOError: Raised if no spectral-index configuration files
            could be found.
        KeyError: Raised if the spectral-index database keywords are not
            all unique.

    .. todo::
        - Somehow add a python call that reads the databases and
          constructs the table for presentation in sphinx so that the
          text above doesn't have to be edited with changes in the
          available databases.
        
    """
    # Check the source directory exists
    dapsrc = dap_source_dir() if dapsrc is None else str(dapsrc)
    if not os.path.isdir(dapsrc):
        raise NotADirectoryError('{0} does not exist!'.format(dapsrc))

    # Check the configuration files exist
    search_dir = os.path.join(dapsrc, 'python/mangadap/config/spectral_indices')
    ini_files = glob.glob(os.path.join(search_dir, '*.ini'))
    if len(ini_files) == 0:
        raise IOError('Could not find any configuration files in {0} !'.format(search_dir))

    # Build the list of library definitions
    index_set_list = []
    for f in ini_files:
        # Read the config file
        cnfg = DefaultConfig(f=f)
        # Ensure it has the necessary elements to define the template
        # library
        validate_spectral_indices_config(cnfg)

        index_set_list += [ SpectralIndicesDef(cnfg['key'],
                                               cnfg.getfloat('minimum_snr', default=0.), 
                                               cnfg.getfloat('resolution_fwhm', default=-1),
                                               cnfg.getbool('compute_sigma_correction',
                                                            default=False),
                                               cnfg['artifact_mask'], cnfg['absorption_indices'],
                                               cnfg['bandhead_indices']) ]

    # Check the keywords of the libraries are all unique
    if len(numpy.unique(numpy.array([index['key'] for index in index_set_list]))) \
                != len(index_set_list):
        raise KeyError('Spectral-index set keywords are not all unique!')

    # Return the default list of assessment methods
    return index_set_list


class SpectralIndicesBitMask(BitMask):
    r"""

    Derived class that specifies the mask bits for the spectral-index
    measurements.  See :class:`mangadap.util.bitmask.BitMask` for
    attributes.

    A list of the bits and meanings are provided by the base class
    function :func:`mangadap.util.bitmask.BitMask.info`; i.e.,::

        from mangadap.proc.spectralindices import SpectralIndicesBitMask
        bm = SpectralIndicesBitMask()
        bm.info()

    """
    def __init__(self, dapsrc=None):
        dapsrc = dap_source_dir() if dapsrc is None else str(dapsrc)
        BitMask.__init__(self, ini_file=os.path.join(dapsrc, 'python', 'mangadap', 'config',
                                                     'bitmasks', 'spectral_indices_bits.ini'))


# TODO: These two should have the same base class
class AbsorptionLineIndices:
    """
    Measure a set of spectral indices in a single spectrum.

    By default, the center of the two side-bands is the flux-weighted
    center; set weighted_center=False to get use the unweighted center
    of the band.
    """
    def __init__(self, wave, flux, bluebands, redbands, mainbands, err=None, log=True, units=None,
                 weighted_center=True):

        # Check the shape of the input spectrum
        if len(flux.shape) != 1:
            raise ValueError('Input flux must be a single vector!')
        if wave.shape != flux.shape:
            raise ValueError('Input flux and wavelength vectors must have the same shape.')
        
        if len(bluebands.shape) != 2:
            raise ValueError('Band definitions must be two-dimensional: Nindx x 2.')
        if bluebands.shape != redbands.shape or bluebands.shape != mainbands.shape:
            raise ValueError('Input bands must have identical shape.')

        self.nindx = bluebands.shape[0]
        # Make sure the units exist
        if units is None:
            warnings.warn('Input units not provided.  Assuming angstroms.')
            self.units = numpy.full(self.nindx, 'ang', dtype=object)
        else:
            self.units = units
        self.order = None

        # Get the two pseudocontinua and the flux-weighted band centers
        # pseudocontinuum returns center, continuum, and continuum_err
        # as MaskedArrays; by default center is the flux-weighted center
        self.blue_center, self.blue_continuum, self.blue_continuum_err, self.blue_incomplete, \
            self.blue_empty = pseudocontinuum(wave, flux, passband=bluebands, err=err, log=log,
                                              weighted_center=weighted_center)

        self.red_center, self.red_continuum, self.red_continuum_err, self.red_incomplete, \
            self.red_empty = pseudocontinuum(wave, flux, passband=redbands, err=err, log=log,
                                             weighted_center=weighted_center)

        # Get the parameters for the linear continuum across the
        # primary passband
        self.continuum_m = (self.red_continuum - self.blue_continuum) \
                                / (self.red_center - self.blue_center)
        self.continuum_b = self.blue_continuum - self.blue_center * self.continuum_m

#        print('Number of non-finite continuum m,b: {0} {1}'.format(
#                                numpy.sum(numpy.invert(numpy.isfinite(self.continuum_m))),
#                                numpy.sum(numpy.invert(numpy.isfinite(self.continuum_b)))))

        # Compute the continuum normalized indices.  This has to be done
        # in a for loop because the continuum is index-dependent
        self.index = numpy.zeros(self.nindx, dtype=numpy.float)
        self.index_err = numpy.zeros(self.nindx, dtype=numpy.float)
        self.main_incomplete = numpy.zeros(self.nindx, dtype=numpy.bool)
        self.main_empty = numpy.zeros(self.nindx, dtype=numpy.bool)
        self.divbyzero = numpy.zeros(self.nindx, dtype=numpy.bool)

        lg10 = numpy.log(10.0)

        for i,m in enumerate(mainbands):
            if self.blue_empty[i] or self.red_empty[i]:
                continue
                
            # From Worthey et al. 1994, eqns. 2 and 3
            cont = self.continuum_b[i] + self.continuum_m[i]*wave
            integrand = 1.0 - flux/cont if self.units[i] == 'ang' else flux/cont

            # Calculate the integral over the passband
            self.index[i] = passband_integral(wave, integrand, passband=m, log=log)
            if err is not None:
                self.index_err[i] = numpy.sqrt(passband_integral(wave, numpy.square(err/cont),
                                               passband=m, log=log))

            # Get the fraction of the band covered by the spectrum and
            # flag bands that are only partially covered or empty
            interval = passband_integrated_width(wave, flux, passband=m, log=log)
            interval_frac = interval / numpy.diff(m)[0]
            self.main_incomplete[i] = interval_frac < 1.0
            self.main_empty[i] = numpy.invert(interval_frac > 0.0)

#            # TEST: Change calculation
#            cont = self.continuum_b[i] + self.continuum_m[i]*wave
#            fint = passband_integral(wave, flux, passband=m, log=log)
#            cint = passband_integral(wave, cont, passband=m, log=log)
#
##            print(interval*(1-fint/cint) - self.index[i])
#            self.index[i] = interval*(1.0-fint/cint if self.units[i] == 'ang' else fint/cint)
#            if err is not None:
#                eint = numpy.sqrt(passband_integral(wave, numpy.square(err), passband=m, log=log))
#                self.index_err[i] = abs(interval*eint/cint)

            # Convert the index to magnitudes using Worthey et al. 1994,
            # eqn. 3: The passband interval cancels out of the error
            # propagation.  The error calculation is done first so as to
            # not replace the linear calculation of the index.
            if self.units[i] == 'mag':
                self.divbyzero[i] = not ((numpy.absolute(interval) > 0) & (self.index[i] > 0))
                if not self.divbyzero[i]:
                    # If err is None, self.index_err=0
                    self.index_err[i] = numpy.absolute(2.5*self.index_err[i]/self.index[i]/lg10)
                    self.index[i] = -2.5 * numpy.log10(self.index[i]/interval)


class BandheadIndices:
    """
    Measure a set of bandhead indices in a single spectrum.
    """
    def __init__(self, wave, flux, bluebands, redbands, err=None, log=True, order=None):

        # Check the shape of the input spectrum
        if len(flux.shape) != 1:
            raise ValueError('Input flux must be a single vector!')
        if wave.shape != flux.shape:
            raise ValueError('Input flux and wavelength vectors must have the same shape.')
        
        if len(bluebands.shape) != 2:
            raise ValueError('Band definitions must be two-dimensional: Nindx x 2.')
        if bluebands.shape != redbands.shape:
            raise ValueError('Input bands must have identical shape.')

        self.nindx = bluebands.shape[0]
        self.units = None
        # Check input order
        if order is None:
            warnings.warn('Input order not specified.  Assuming r_b.')
            self.order = numpy.full(self.nindx, 'r_b', dtype=object)
        else:
            self.order = order

        # Get the two pseudocontinua and the flux-weighted band centers
        # pseudocontinuum returns center, continuum, and continuum_err
        # as MaskedArrays
        self.blue_center, self.blue_continuum, self.blue_continuum_err, self.blue_incomplete, \
            self.blue_empty = pseudocontinuum(wave, flux, passband=bluebands, err=err, log=log)

        self.red_center, self.red_continuum, self.red_continuum_err, self.red_incomplete, \
            self.red_empty = pseudocontinuum(wave, flux, passband=redbands, err=err, log=log)

        # Determine which indices have both a valid index and index
        # error calculation
        self.main_incomplete = None
        self.main_empty = None

#        self.divbyzero = numpy.invert(numpy.absolute(self.red_continuum.filled(0.0)) > 0) \
#                            | numpy.invert(numpy.absolute(self.blue_continuum.filled(0.0)) > 0)
        self.divbyzero = numpy.zeros(self.nindx, dtype=bool)

        # Calculate the index in the correct order
        blue_n = order == 'b_r'
        self.index = numpy.ma.zeros(self.nindx, dtype=numpy.float)
        self.index[blue_n] = numpy.ma.divide(self.blue_continuum[blue_n],
                                             self.red_continuum[blue_n]).filled(0.0)
        self.divbyzero[blue_n] = numpy.invert(numpy.absolute(self.red_continuum[blue_n])>0.0)
        self.index[numpy.invert(blue_n)] = numpy.ma.divide(self.red_continuum[numpy.invert(blue_n)],
                                              self.blue_continuum[numpy.invert(blue_n)]).filled(0.0)
        self.divbyzero[numpy.invert(blue_n)] \
                = numpy.invert(numpy.absolute(self.blue_continuum[numpy.invert(blue_n)])>0.0)

        # Calculate the index errors
        berr = numpy.ma.zeros(self.nindx) if self.blue_continuum_err is None \
                                            else self.blue_continuum_err
        rerr = numpy.ma.zeros(self.nindx) if self.red_continuum_err is None \
                                            else self.red_continuum_err
        # Error is independent of ratio order when written in this way
        self.index_err = numpy.ma.sqrt(
                            numpy.square(numpy.ma.divide(berr*self.index,self.blue_continuum)) +
                            numpy.square(numpy.ma.divide(rerr*self.index,self.red_continuum))
                                      ).filled(0.0)



#        n = numpy.ma.zeros(self.nindx, dtype=numpy.float)
#        d = numpy.ma.zeros(self.nindx, dtype=numpy.float)
#        n[blue_n] = self.blue_continuum[blue_n]
#        d[blue_n] = self.red_continuum[blue_n]
#        n[~blue_n] = self.red_continuum[~blue_n]
#        d[~blue_n] = self.blue_continuum[~blue_n]
#
#        nerr = numpy.ma.zeros(self.nindx)
#        derr = numpy.ma.zeros(self.nindx)
#        if self.blue_continuum_err is not None:
#            nerr[blue_n] = self.blue_continuum_err[blue_n]
#            derr[~blue_n] = self.blue_continuum_err[~blue_n]
#        if self.red_continuum_err is not None:
#            nerr[~blue_n] = self.red_continuum_err[~blue_n]
#            derr[blue_n] = self.red_continuum_err[blue_n]
#
#        # Determine which indices have both a valid index and index
#        # error calculation
#        self.main_incomplete = None
#        self.main_empty = None
#        self.divbyzero = ~((numpy.absolute(d) > 0) & (numpy.absolute(n) > 0))
#
#        # Calculate the indices and their nominal errors
#        self.index = numpy.zeros(self.nindx, dtype=numpy.float)
#        self.index[~self.divbyzero] = n[~self.divbyzero]/d[~self.divbyzero]
#        self.index_err = numpy.zeros(self.nindx, dtype=numpy.float)
#        self.index_err[~self.divbyzero] = numpy.sqrt(
#                            numpy.square(nerr[~self.divbyzero]*self.index[~self.divbyzero]
#                                            / n[~self.divbyzero])
#                          + numpy.square(derr[~self.divbyzero]*self.index[~self.divbyzero]
#                                            / d[~self.divbyzero]) )


class SpectralIndices:
    r"""
    Class that computes and interfaces with the spectral-index
    measurements.

    If neither stellar-continuum nor emission-line models are provided:
        - Indices are measure on the binned spectra
        - No velocity-dispersion corrections are calculated

    If a stellar-continuum model is provided without an emission-line
    model:
        - Indices are measured on the binned spectra
        - Velocity-dispersion corrections are computed for any binned
          spectrum with a stellar-continuum fit based on the optimal
          template

    If an emission-line model is provided without a stellar-continuum
    model:
        - Indices are measured on the relevant (binned or unbinned)
          spectra; spectra with emission-line fits have the model
          emission lines subtracted from them before these measurements.
        - If the emission-line model includes data regarding the
          stellar-continuum fit (template spectra and template weights),
          corrections are calculated for spectra with emission-line
          models based on the continuum fits; otherwise, no corrections
          are calculated.

    If both stellar-continuum and emission-line models are provided, and
    if the stellar-continuum and emission-line fits are performed on the
    same spectra:
        - Indices are measured on the relevant (binned or unbinned)
          spectra; spectra with emission-line fits have the model
          emission lines subtracted from them before these measurements.
        - Velocity-dispersion corrections are based on the
          stellar-continuum templates and weights

    If both stellar-continuum and emission-line models are provided, and
    if the stellar-continuum and emission-line fits are performed on
    different spectra:
        - The behavior is exactly as if the stellar-continuum model was
          not provided.

    **Detail what should be provided in terms of the redshift**

    Args:
        database_key (str): Keyword used to select the specfic list of
            indices to measure and how they should be measured;  see
            :class:`SpectralIndicesDef`.
        binned_spectra
            (:class:`mangadap.proc.spatiallybinnedspectra.SpatiallyBinnedSpectra`):
            The binned spectra for the measurements.
        redshift (float, numpy.ndarray): (**Optional**) A single or
            spectrum-dependent redshift, :math:`z`, to use for shifting
            the index bands.  Default is to measure the indices at their
            provided wavelengths (i.e., :math:`z=0`).
        stellar_continuum
            (:class:`mangadap.proc.stellarcontinuummodel.StellarContinuumModel`):
            (**Optional**) The stellar-continuum model as applied to the
            binned spectra.
        emission_line_model
            (:class:`mangadap.proc.emissionlinemodel.EmissionLineModel`):
            (**Optional**) The emission-line model as applied to either
            the binned spectra or the unbinned spaxels.
        database_list (list): (**Optional**) List of
            :class:`SpectralIndicesDef` objects that define one or more
            methods to use for the spectral-index measurements.  The
            default list is provided by the config files in the DAP
            source directory and compiled into this list using
            :func:`available_spectral_index_databases`.
        artifact_list (list): (**Optional**) List of
            :class:`mangadap.par.spectralfeaturedb.SpectralFeatureDBDef`
            objects that define the unique key for the artifact database
            (see :mod:`mangadap.par.artifactdb`).
        absorption_index_list (list): (**Optional**) List of
            :class:`mangadap.par.spectralfeaturedb.SpectralFeatureDBDef`
            objects that define the unique key for the absorption-index
            database (see :mod:`mangadap.par.absorptionindexdb`).
        bandhead_index_list (list): (**Optional**) List of
            :class:`mangadap.par.spectralfeaturedb.SpectralFeatureDBDef`
            objects that define the unique key for the bandhead-index
            database (see :mod:`mangadap.par.bandheadindexdb`).
        dapsrc (str): (**Optional**) Root path to the DAP source
            directory.  If not provided, the default is defined by
            :func:`mangadap.config.defaults.dap_source_dir`.
        dapver (str): (**Optional**) The DAP version to use for the
            analysis, used to override the default defined by
            :func:`mangadap.config.defaults.default_dap_version`.
        analysis_path (str): (**Optional**) The top-level path for the
            DAP output files, used to override the default defined by
            :func:`mangadap.config.defaults.default_analysis_path`.
        directory_path (str): The exact path to the directory with DAP
            output that is common to number DAP "methods".  See
            :attr:`directory_path`.
        output_file (str): (**Optional**) Exact name for the output
            file.  The default is to use
            :func:`mangadap.config.defaults.default_dap_file_name`.
        hardcopy (bool): (**Optional**) Flag to write the HDUList
            attribute to disk.  Default is True; if False, the HDUList
            is only kept in memory and would have to be reconstructed.
        tpl_symlink_dir (str): (**Optional**) Create a symbolic link to
            the created template library file in the supplied directory.
            Default is to produce no symbolic link.
        clobber (bool): (**Optional**) Overwrite any existing files.
            Default is to use any existing file instead of redoing the
            analysis and overwriting the existing output.
        checksum (bool): (**Optional**) Use the checksum in the fits
            header to confirm that the data has not been corrupted.  The
            checksum is **always** written to the fits header when the
            file is created; this argument does not toggle that
            functionality.
        loggers (list): (**Optional**) List of `logging.Logger`_ objects
            to log progress; ignored if quiet=True.  Logging is done
            using :func:`mangadap.util.log.log_output`.  Default is no
            logging.
        quiet (bool): (**Optional**) Suppress all terminal and logging
            output.  Default is False.

    """
    def __init__(self, database_key, binned_spectra, redshift=None, stellar_continuum=None,
                 emission_line_model=None, database_list=None, artifact_list=None,
                 absorption_index_list=None, bandhead_index_list=None, dapsrc=None, dapver=None,
                 analysis_path=None, directory_path=None, output_file=None, hardcopy=True,
                 tpl_symlink_dir=None, clobber=False, checksum=False, loggers=None, quiet=False):

        self.loggers = None
        self.quiet = False

        # Define the database properties
        self.database = None
        self.artdb = None
        self.pixelmask = None
        self.absdb = None
        self.bhddb = None
        self._define_databases(database_key, database_list=database_list,
                               artifact_list=artifact_list,
                               absorption_index_list=absorption_index_list,
                               bandhead_index_list=bandhead_index_list)

        self.nabs, self.nbhd = self.count_indices(self.absdb, self.bhddb)
        self.nindx = self.nabs + self.nbhd
        self.compute_corrections = False
        self.correct_indices = False            # !! HARDCODED !!

        self.binned_spectra = None
        self.redshift = None
        self.stellar_continuum = None
        self.emission_line_model = None

        # Define the output directory and file
        self.directory_path = None              # Set in _set_paths
        self.output_file = None
        self.hardcopy = None

        # Initialize the objects used in the assessments
        self.bitmask = SpectralIndicesBitMask(dapsrc=dapsrc)

        self.hdu = None
        self.checksum = checksum
        self.spatial_shape = None
        self.nspec = None
        self.spatial_index = None
        self.image_arrays = None
        self._assign_image_arrays()

        self.nbins = None
        self.missing_bins = None

        # Run the assessments of the DRP file
        self.measure(binned_spectra, redshift=redshift, stellar_continuum=stellar_continuum,
                     emission_line_model=emission_line_model, dapsrc=dapsrc, dapver=dapver,
                     analysis_path=analysis_path, directory_path=directory_path,
                     output_file=output_file, hardcopy=hardcopy, tpl_symlink_dir=tpl_symlink_dir,
                     clobber=clobber, loggers=loggers, quiet=quiet)


    def __getitem__(self, key):
        return self.hdu[key]


    def _define_databases(self, database_key, database_list=None, artifact_list=None,
                          absorption_index_list=None, bandhead_index_list=None, dapsrc=None):
        r"""

        Select the database of indices

        """
        # Grab the specific database
        self.database = select_proc_method(database_key, SpectralIndicesDef,
                                           method_list=database_list,
                                           available_func=available_spectral_index_databases,
                                           dapsrc=dapsrc)

        # Instantiate the artifact, absorption-index, and bandhead-index
        # databases
        self.artdb = None if self.database['artifacts'] is None else \
                ArtifactDB(self.database['artifacts'], artdb_list=artifact_list, dapsrc=dapsrc)
        # TODO: Generalize the name of this object
        self.pixelmask = SpectralPixelMask(artdb=self.artdb)

        self.absdb = None if self.database['absindex'] is None else \
                AbsorptionIndexDB(self.database['absindex'], indxdb_list=absorption_index_list,
                                  dapsrc=dapsrc)

        self.bhddb = None if self.database['bandhead'] is None else \
                BandheadIndexDB(self.database['bandhead'], indxdb_list=bandhead_index_list,
                                dapsrc=dapsrc)


    def _set_paths(self, directory_path, dapver, analysis_path, output_file):
        """
        Set the :attr:`directory_path` and :attr:`output_file`.  If not
        provided, the defaults are set using, respectively,
        :func:`mangadap.config.defaults.default_dap_common_path` and
        :func:`mangadap.config.defaults.default_dap_file_name`.

        Args:
            directory_path (str): The exact path to the DAP
                spectral-index file.  See :attr:`directory_path`.
            dapver (str): DAP version.
            analysis_path (str): The path to the top-level directory
                containing the DAP output files for a given DRP and DAP
                version.
            output_file (str): The name of the file with spectral-index
                moment measurements.  See :func:`measure`.
        """
        # Set the output directory path
        method = default_dap_method(binned_spectra=self.binned_spectra,
                                    stellar_continuum=self.stellar_continuum)
        self.analysis_path = default_analysis_path if analysis_path is None else str(analysis_path)
        self.directory_path = default_dap_method_path(method, plate=self.binned_spectra.drpf.plate,
                                                      ifudesign=self.binned_spectra.drpf.ifudesign,
                                                      ref=True,
                                                      drpver=self.binned_spectra.drpf.drpver,
                                                      dapver=dapver,
                                                      analysis_path=self.analysis_path) \
                                        if directory_path is None else str(directory_path)

        # Set the output file
        ref_method = '{0}-{1}'.format(self.binned_spectra.rdxqa.method['key'],
                                      self.binned_spectra.method['key'])
        if self.stellar_continuum is not None:
            ref_method = '{0}-{1}'.format(ref_method, self.stellar_continuum.method['key'])
        if self.emission_line_model is not None:
            ref_method = '{0}-{1}'.format(ref_method, self.emission_line_model.method['key'])
        ref_method = '{0}-{1}'.format(ref_method, self.database['key'])

        self.output_file = default_dap_file_name(self.binned_spectra.drpf.plate,
                                                 self.binned_spectra.drpf.ifudesign, ref_method) \
                                        if output_file is None else str(output_file)


    def _initialize_primary_header(self, hdr=None, measurements_binid=None):
        """
        Initialize the header of :attr:`hdu`.

        Returns:
            astropy.io.fits.Header : Edited header object.

        """
        # Copy the from the DRP and clean it
        if hdr is None:
            hdr = self.binned_spectra.drpf.hdu['PRIMARY'].header.copy()
            hdr = DAPFitsUtil.clean_dap_primary_header(hdr)
        
        hdr['AUTHOR'] = 'Kyle B. Westfall <westfall@ucolick.org>'
        hdr['SIKEY'] = (self.database['key'], 'Spectral-index database keyword')
        hdr['SIMINSN'] = (self.database['minimum_snr'], 'Minimum S/N of spectrum to include')
        hdr['SIFWHM'] = (self.database['fwhm'], 'FWHM of index system resolution (ang)')
        hdr['ARTDB'] = (self.database['artifacts'], 'Artifact database keyword')
        hdr['ABSDB'] = (self.database['absindex'], 'Absorption-index database keyword')
        hdr['BHDDB'] = (self.database['bandhead'], 'Bandhead-index database keyword')
        if self.stellar_continuum is not None:
            hdr['SCKEY'] = (self.stellar_continuum.method['key'], 'Stellar-continuum model keyword')
        if self.emission_line_model is not None:
            hdr['ELFKEY'] = (self.emission_line_model.method['key'],
                                'Emission-line modeling method keyword')
        hdr['NBINS'] = (self.nbins, 'Number of unique spatial bins')
        hdr['SICORR'] = (self.compute_corrections, 'Velocity dispersion corrections computed')
        hdr['SIREBIN'] = (measurements_binid is not None, 'Bin IDs disconnected from SC binning')
        return hdr


#    def _initialize_mask(self, good_snr):
#        """
#
#        Initialize the mask be setting the DIDNOTUSE, FORESTAR, and LOW_SNR masks
#
#        """
#        # Initialize to all zeros
#        mask = numpy.zeros(self.shape, dtype=self.bitmask.minimum_dtype())
#
#        # Turn on the flag stating that the pixel wasn't used
##        print('Mask size:', self.binned_spectra['MASK'].data.size)
#        indx = self.binned_spectra.bitmask.flagged(self.binned_spectra['MASK'].data,
#                                                   flag=self.binned_spectra.do_not_fit_flags())
##        print('Masked as DIDNOTUSE:', numpy.sum(indx))
#        mask[indx] = self.bitmask.turn_on(mask[indx], 'DIDNOTUSE')
#
#        # Turn on the flag stating that the pixel has a foreground star
#        indx = self.binned_spectra.bitmask.flagged(self.binned_spectra['MASK'].data,
#                                                   flag='FORESTAR')
##        print('Masked as FORESTAR: ', numpy.sum(indx))
#        mask[indx] = self.bitmask.turn_on(mask[indx], 'FORESTAR')
#
#        # Turn on the flag stating that the S/N in the spectrum was
#        # below the requested limit
#        indx = numpy.array([numpy.invert(good_snr).reshape(self.spatial_shape).T]*self.nwave).T
#        mask[indx] = self.bitmask.turn_on(mask[indx], flag='LOW_SNR')
#
#        return mask


    def _index_database_dtype(self, name_len):
        r"""
        Construct the record array data type for the output fits
        extension.
        """
        return [ ('TYPE','<U10'),
                 ('ID',numpy.int),
                 ('NAME','<U{0:d}'.format(name_len)),
                 ('PASSBAND', numpy.float, (2,)),
                 ('BLUEBAND', numpy.float, (2,)),
                 ('REDBAND', numpy.float, (2,)),
                 ('UNIT', '<U3'),
                 ('COMPONENT', numpy.uint8),
                 ('INTEGRAND', '<U7'),
                 ('ORDER', '<U3')
               ]

    
    def _compile_database(self):
        """
        Compile the database with the specifications of each index.
        """
        name_len = 0
        if self.absdb is not None:
            for n in self.absdb['name']:
                if name_len < len(n):
                    name_len = len(n)
        if self.bhddb is not None:
            for n in self.bhddb['name']:
                if name_len < len(n):
                    name_len = len(n)

        # Instatiate the table data that will be saved defining the set
        # of indices measured
        passband_database = init_record_array(self.nindx, self._index_database_dtype(name_len))

        t = 0 if self.absdb is None else self.absdb.nsets
        if self.absdb is not None:
            passband_database['TYPE'][:t] = 'absorption'
            hk = [ 'ID', 'NAME', 'PASSBAND', 'BLUEBAND', 'REDBAND', 'UNIT', 'COMPONENT' ]
            ak = [ 'index', 'name', 'primary', 'blueside', 'redside', 'units', 'component' ]
            for _hk, _ak in zip(hk,ak):
                passband_database[_hk][:t] = self.absdb[_ak]
        if self.bhddb is not None:
            passband_database['TYPE'][t:] = 'bandhead'
            hk = [ 'ID', 'NAME', 'BLUEBAND', 'REDBAND', 'INTEGRAND', 'ORDER' ]
            ak = [ 'index', 'name', 'blueside', 'redside', 'integrand', 'order' ]
            for _hk, _ak in zip(hk,ak):
                passband_database[_hk][t:] = self.bhddb[_ak]
        return passband_database


    def _assign_image_arrays(self):
        """
        Set :attr:`image_arrays`, which contains the list of extensions
        in :attr:`hdu` that are on-sky image data.
        """
        self.image_arrays = [ 'BINID' ]


    def _get_missing_bins(self, unique_bins=None):
        if unique_bins is None:
            good_snr = self.binned_spectra.above_snr_limit(self.database['minimum_snr'])
            return numpy.sort(self.binned_spectra.missing_bins + 
                        self.binned_spectra['BINS'].data['BINID'][numpy.invert(good_snr)].tolist())
        return SpatiallyBinnedSpectra._get_missing_bins(unique_bins)


    def _assign_redshifts(self, redshift, measure_on_unbinned_spaxels, good_snr,
                          default_redshift=None):
        """
        Set the redshift to use for each spectrum for the spectral index
        measurements.
        
        In terms of precedence, directly provided redshifts override
        those in any available StellarContinuumModel.

        If self.stellar_continuum and redshift are None, the default
        redshift is used (or 0.0 if this is also None).

        To get the stellar kinematics, the function calls
        :func:`mangadap.proc.stellarcontinuummodel.StellarContinuumModel.matched_kinematics`.
        It is expected that the stellar kinematics were fixed to these
        values during any emission-line modeling that may have altered
        the continuum fit itself (e.g., :class:`mangadap.proc.Sasuke`).

        In this function, the provided redshift must be a single value
        or None; therefore, the means of any vectors should be provided
        instead of the full vector.

        The function is borrows heavily from
        :func:`mangadap.proc.emissionlinemodel.EmissionLineModel._assign_input_kinematics`.

        Args:

            redshift (float, numpy.ndarray): Redshifts (:math:`z`) to
                use for each spectrum.  If None, the default

            default_redshift (float): (**Optional**) Only used if there
                are stellar kinematics available.  Provides the default
                redshift to use for spectra without stellar
                measurements; see :arg:`redshift` in
                :func:`mangadap.proc.stellarcontinuummodel.StellarContinuumModel.matched_kinematics`.
                If None (default), the median of the unmasked stellar
                velocities will be used.
        """

        # Construct the binid matrix if measuring on unbinned spaxels
        if measure_on_unbinned_spaxels:
            binid = numpy.full(self.binned_spectra.spatial_shape, -1, dtype=int)
            binid.ravel()[good_snr] = numpy.arange(self.nbins)
            missing = []
            nspec = self.binned_spectra.drpf.nspec
        else:
            binid = self.binned_spectra['BINID'].data
            missing = self.binned_spectra.missing_bins
            nspec = self.binned_spectra.nbins

        #---------------------------------------------------------------
        # Get the redshift measured for the stars and use them if no
        # default value is provided
        if self.stellar_continuum is not None and redshift is None:
            self.redshift, _ = self.stellar_continuum.matched_kinematics(
                                                binid, redshift=default_redshift,
                                                nearest=True, missing=missing)
            if measure_on_unbinned_spaxels:
                tmp = self.redshift.copy()
                self.redshift = numpy.zeros(nspec, dtype=float)
                self.redshift[good_snr] = tmp
            return

        #---------------------------------------------------------------
        # Use the default value(s)
        _redshift = numpy.atleast_1d(redshift)
        if len(_redshift) not in [ 1, nspec ]:
            raise ValueError('Provided redshift must be either a single value or match the '
                             'number of binned spectra or the number of unbinned spaxels.')
        self.redshift = numpy.full(nspec, redshift, dtype=float) \
                                if len(_redshift) == 1 else _redshift.copy()


    def _flag_good_spectra(self, measure_on_unbinned_spaxels):
        if measure_on_unbinned_spaxels:
            fgoodpix = self.binned_spectra.check_fgoodpix()
            good_snr = self.binned_spectra.rdxqa['SPECTRUM'].data['SNR'] \
                                > self.database['minimum_snr']
            return fgoodpix & good_snr
        return self.binned_spectra.above_snr_limit(self.database['minimum_snr'])


    @staticmethod
    def spectra_for_index_measurements(binned_spectra, measure_on_unbinned_spaxels=False,
                                       pixelmask=None, select=None, resolution_fwhm=None,
                                       emission_line_model=None):
        """
        Compile the set of spectra for the spectral-index measurements.

        If the input fwhm is > 0, this function will match the spectral
        resolution of the data to the spectral-index system, based on
        the provided FWHM.  It also subtracts the emission-line model if
        provided.

        .. todo::
            Allow resolution_fwhm to be wavelength dependent, provided
            via a vector.

        Args:
            binned_spectra
                (:class:`mangadap.proc.spatiallybinnedspectra.SpatiallyBinnedSpectra`):
                The binned spectra object.  The returned spectra are
                either the binned spectra or the DRP spectra internal to
                this object.
            measure_on_unbinned_spaxels (bool): (**Optional**) Flag to
                return the unbinned spaxels as opposed to the binned
                spectra.  Default is to use the binned spectra.
            pixelmask
                (:class:`mangadap.util.pixelmask.SpectralPixelMask`):
                (**Optional**) Defines the pixels that should
                automatically be masked during the measurements.  By
                default, nothing is masked in addition to that specified
                by the data mask.
            select (numpy.ndarray): (**Optional**) Boolean vector
                selecting which spectra to return.  The length must
                match the number of binned spectra or the total number
                of DRP spectra, depending on the provided
                :arg:`measure_on_unbinned_spaxels`.  Default is to
                return all spectra.

            resolution_fwhm (float): (**Optional**)
                Wavelength-independent FWHM of the resolution element at
                which to measure the indices.  If
                > 0, the spectra are resolution matched from the input
                resolution to the provided resolution; otherwise, the
                resolution is not altered.

            emission_line_model
                (:class:`mangadap.proc.emissionlinemodel.EmissionLineModel`):
                (**Optional**) Object providing the emission-line model.
                The emission-line model must match the selection of the
                spectra (binned or unbinned) to fit, as given by the
                fitting method (`deconstruct_bins`).

        Returns:
            numpy.ndarray, numpy.ma.MaskedArray: Three arrays are
            returned: (1) the common wavelength vector, (2) the masked
            flux array, and (3) the masked inverse variance array.

        Raises:
            ValueError: Raised if the emission-line model and spectra
                selection (binned vs. unbinned) do not match.
        """

        # Check that the spectrum selection and emission-line model are
        # consistent
        if measure_on_unbinned_spaxels and emission_line_model is not None \
                and not emission_line_model.method['deconstruct_bins']:
            raise ValueError('Cannot use this emission-line model with unbinned spaxels.')

        # Get the main data arrays
        if measure_on_unbinned_spaxels:
            # TODO: Should probably make this a function within
            # SpatiallyBinnedSpectra, particularly because of the
            # dereddening
            flags = binned_spectra.drpf.do_not_fit_flags()
            binid = numpy.arange(binned_spectra.drpf.nspec).reshape(
                                                                binned_spectra.drpf.spatial_shape)
            wave = binned_spectra.drpf['WAVE'].data
            flux = binned_spectra.drpf.copy_to_masked_array(flag=flags)
            ivar = binned_spectra.drpf.copy_to_masked_array(ext='IVAR', flag=flags)
            flux, ivar = binned_spectra.galext.apply(flux, ivar=ivar, deredden=True)
            missing = None
            sres = binned_spectra.drpf.spectral_resolution(
                        ext= 'SPECRES' if binned_spectra.method['spec_res'] == 'cube' else None,
                        toarray=True, pre=binned_spectra.method['prepixel_sres'], fill=True).data
        else:
            flags = binned_spectra.do_not_fit_flags()
            binid = binned_spectra['BINID'].data
            wave = binned_spectra['WAVE'].data
            flux = binned_spectra.copy_to_masked_array(flag=flags)
            ivar = binned_spectra.copy_to_masked_array(ext='IVAR', flag=flags)
            missing = binned_spectra.missing_bins
            sres = binned_spectra.copy_to_array(ext='SPECRES')

        # Total number of spectra
        nspec = flux.shape[0]

        # Mask any pixels in the pixel mask
        if pixelmask is not None:
            indx = pixelmask.boolean(wave, nspec=nspec)
            flux[indx] = numpy.ma.masked
            ivar[indx] = numpy.ma.masked

        # Remove the emission lines if provided        
#        warnings.warn('DEBUG')
        if emission_line_model is not None:
#            pyplot.imshow(flux, origin='lower', interpolation='nearest', aspect='auto')
#            pyplot.show()
            eml_model = emission_line_model.fill_to_match(binid, missing=missing)
#            pyplot.imshow(eml_model, origin='lower', interpolation='nearest', aspect='auto')
#            pyplot.show()
            no_eml = numpy.invert(numpy.ma.getmaskarray(flux)) & numpy.ma.getmaskarray(eml_model)
            flux -= eml_model
            flux.mask[no_eml] = False
#            pyplot.imshow(flux, origin='lower', interpolation='nearest', aspect='auto')
#            pyplot.show()
#            exit()

        # Make sure ivar mask is identical to flux mask
        ivar.mask = flux.mask.copy()

        # Set the selected spectra
        _select = numpy.ones(nspec, dtype=bool) if select is None else select

        # Adjust the spectral resolution
        _resolution_fwhm = -1 if resolution_fwhm is None else resolution_fwhm
        if _resolution_fwhm > 0:
            flux, ivar = adjust_spectral_resolution(wave, flux, ivar, sres, _resolution_fwhm)

        return wave, flux[select,:], ivar[select,:]


    @staticmethod
    def adjust_spectral_resolution(wave, flux, ivar, sres, resolution_fwhm):
        """
        flux and ivar are expected to be masked arrays
        """
        # Revert flux and ivar to unmasked arrays
        mask = numpy.ma.getmaskarray(flux).copy()
        flux = numpy.asarray(flux)
        ivar = numpy.asarray(ivar)
        nspec = flux.shape[0]

        # For masked pixels, interpolate across them using simple linear
        # interpolation.
        # TODO: This approach needs to be vetted.
        for i in range(nspec):
            mp = mask[i,:]              # Masked pixels
            up = numpy.invert(mp)       # Unmasked pixels
            if numpy.sum(mp) == 0:      # No masked pixels so continue
                continue
            if numpy.sum(up) == 0:      # No unmasked pixels!
                flux[i,:] = 0.0
                continue
            # Linearly interpolate both the flux and inverse variance
            flux[i,mp] = numpy.interp(wave[mp], wave[up], flux[i,up])
            ivar[i,mp] = numpy.interp(wave[mp], wave[up], ivar[i,up])

#        pyplot.step(wave, flux[0,:], where='mid', linestyle='-', color='r', lw=0.5)
#        pyplot.show()

        # Use :func:`mangadap.instrument.match_spectral_resolution` to
        # match the spectral resolution of the binned spectra to the
        # spectral-index system
        new_sres = wave/resolution_fwhm
        
        _flux, _sres, sigoff, _mask, _ivar \
                = match_spectral_resolution(wave, flux, sres, wave, new_sres, ivar=ivar,
                                            log10=True, new_log10=True)

        # FOR DEBUGGING
#        warnings.warn('NOT MATCHING SPECTRAL RESOLUTION!')
#        new_flux = flux.copy()
#        new_mask = mask.copy()

#        pyplot.step(wave, flux[0,:], where='mid', linestyle='-', color='k', lw=0.5)
#        pyplot.step(wave, new_flux[0,:], where='mid', linestyle='-', color='r', lw=2.5)
#        pyplot.show()
#        pyplot.step(wave, ivar[0,:]/numpy.mean(ivar[0,:]), where='mid', linestyle='-', color='k',
#                    lw=0.5)
#        pyplot.step(wave, new_ivar[0,:]/numpy.mean(new_ivar[0,:]), where='mid', linestyle='-', color='r', lw=2.5)
#        pyplot.show()

        # From instrument.py:  "Any pixel that had a resolution that was
        # lower than the target resolution (up to some tolerance defined
        # by *min_sig_pix*) is returned as masked."  This is exactly
        # what should be masked with the SPECRES_LOW bit.
        if numpy.sum(_mask) > 0:
            mask[_mask > 0] = True
        return numpy.ma.MaskedArray(flux,mask=mask), numpy.ma.MaskedArray(ivar,mask=mask)


    @staticmethod
    def unit_selection(absdb, bhddb):
        """
        Return boolean arrays selection which indices are unitless, in
        angstrom units, or in magnitude units.
        """
        nabs, nbhd = SpectralIndices.count_indices(absdb, bhddb)
        
        # Flag the indices as either having magnitude or angstrom units
        angu = numpy.zeros(nabs + nbhd, dtype=bool)
        magu = numpy.zeros(nabs + nbhd, dtype=bool)
        ules = numpy.zeros(nabs + nbhd, dtype=bool)

        # Flag absorption-line indices as either angstrom or magnitude
        if nabs > 0:
            angu[:nabs] = absdb['units'] == 'ang'
            magu[:nabs] = absdb['units'] == 'mag'

        # Bandhead indices are unitless
        if nbhd > 0:
            ules[nabs:] = True

        return ules, angu, magu


    def _resolution_matched_templates(self, dapsrc=None, dapver=None, analysis_path=None,
                                      tpl_symlink_dir=None):
        """
        Get a version of the template library that has had its
        resolution matched to that of the spectral-index database.
        """
        velocity_offset=self.stellar_continuum.method['fitpar']['template_library'].velocity_offset
        return self.stellar_continuum.get_template_library(dapsrc=dapsrc, dapver=dapver,
                                                           analysis_path=analysis_path,
                                                           tpl_symlink_dir=tpl_symlink_dir,
                                                           velocity_offset=velocity_offset,
                                                           resolution_fwhm=self.database['fwhm'])
    

    @staticmethod
    def calculate_dispersion_corrections(absdb, bhddb, wave, flux, continuum, continuum_dcnvlv,
                                         redshift=None, redshift_dcnvlv=None, bitmask=None):
        """
        Calculate the dispersion corrections using the best-fitting
        template models.

        Allow the "deconvolved" continuum spectra to be at a different
        redshift than the best-fitting continuum.
        """
#        pyplot.step(wave, flux[0,:], where='mid', linestyle='-', color='k', lw=0.5, zorder=1)
#        pyplot.plot(wave, continuum[0,:], linestyle='-', color='g', lw=1.0, zorder=3, alpha=0.5)
#        pyplot.plot(wave, continuum_dcnvlv[0,:], linestyle='-', color='b', lw=1.0, zorder=3,
#                    alpha=0.5)
#        pyplot.show()

        # Make sure continuum includes the flux masks
        nspec = flux.shape[0]
        _flux = numpy.ma.MaskedArray(flux)
        _continuum = numpy.ma.MaskedArray(continuum)
        _continuum_dcnvlv = numpy.ma.MaskedArray(continuum_dcnvlv)
        _continuum[_flux.mask] = numpy.ma.masked
        _continuum_dcnvlv[_flux.mask] = numpy.ma.masked

        # Measure the indices for both sets of spectra
        indx = SpectralIndices.measure_indices(absdb, bhddb, wave, _continuum, redshift=redshift,
                                               bitmask=bitmask)
        dcnvlv_indx = SpectralIndices.measure_indices(absdb, bhddb, wave, _continuum_dcnvlv,
                                                      redshift=redshift_dcnvlv, bitmask=bitmask)

#        print('Flagged D4000:')
#        for bit in ['MAIN_EMPTY', 'BLUE_EMPTY', 'RED_EMPTY', 'DIVBYZERO' ]:
#            print('indx {0}:        {1}/{2}'.format(bit,
#                                            numpy.sum(bitmask.flagged(indx['MASK'][:,43],
#                                            flag=bit)), indx['MASK'].shape[0]))
#            print('dcnvlv_indx {0}: {1}/{2}'.format(bit,
#                                            numpy.sum(bitmask.flagged(dcnvlv_indx['MASK'][:,43],
#                                            flag=bit)), dcnvlv_indx['MASK'].shape[0]))

        # Do not apply the correction if any of the bands were empty or
        # would have had to divide by zero
        bad_indx = bitmask.flagged(indx['MASK'], flag=[ 'MAIN_EMPTY', 'BLUE_EMPTY', 'RED_EMPTY',
                                                        'DIVBYZERO' ])
        bad_indx |= bitmask.flagged(dcnvlv_indx['MASK'], flag=[ 'MAIN_EMPTY', 'BLUE_EMPTY',
                                                                'RED_EMPTY', 'DIVBYZERO' ])

        print('Number of indices: ', bad_indx.size)
        print('Bad indices: ', numpy.sum(bad_indx))

        # Determine which indices have good measurements
        ules, angu, magu = SpectralIndices.unit_selection(absdb, bhddb)
        good_les = numpy.invert(bad_indx) & (numpy.array([ules]*nspec)) \
                        & (numpy.absolute(indx['INDX']) > 0)
        good_ang = numpy.invert(bad_indx) & (numpy.array([angu]*nspec)) \
                        & (numpy.absolute(indx['INDX']) > 0)
        good_mag = numpy.invert(bad_indx) & (numpy.array([magu]*nspec))

        print('Good unitless indices: ', numpy.sum(good_les))
        print('Good angstrom indices: ', numpy.sum(good_ang))
        print('Good magnitude indices: ', numpy.sum(good_mag))

        # Save the *good* indices for the best-fitting models
        model_indices = numpy.zeros(indx['INDX'].shape, dtype=numpy.float)
        model_indices[good_les | good_ang | good_mag] \
                = indx['INDX'][good_les | good_ang | good_mag]

        # Determine the dispersion corrections
        corrections = numpy.zeros(indx['INDX'].shape, dtype=numpy.float)
        corrections[good_les] = dcnvlv_indx['INDX'][good_les] / indx['INDX'][good_les]
        corrections[good_ang] = dcnvlv_indx['INDX'][good_ang] / indx['INDX'][good_ang]
        corrections[good_mag] = dcnvlv_indx['INDX'][good_mag] - indx['INDX'][good_mag]

        # Return the results
        return model_indices, corrections, good_les, good_ang, good_mag


    @staticmethod
    def apply_dispersion_corrections(indx, indxcorr, err=None, unit=None):
        """
        Apply a set of dispersion corrections.  Errors in the dispersion
        corrections are assumed to be negligible.

        Args:
            indx (array-like): Indices to correct.
            indxcorr (array-like): Index corrections.
            err (array-like): (**Optional**) Error in the indices.
            unit (str): (**Optional**) Unit of the index; must be either
                magnitudes (mag) or angstroms (ang) or None.  Default is
                None.  Unitless corrections and angstrom corrections are
                treated identically.
    
        Returns:
            numpy.ndarray: The corrected indices and errors are
            returned.  If no errors are returned, the second returned
            object is None.

        Raises:
            ValueError: Raised if the unit is not ang or mag.

        """
        if unit not in [ None, 'ang', 'mag' ]:
            raise ValueError('Unit must be None, ang, or mag.')

        if unit in [ None, 'ang' ]:
            _indx = indx * indxcorr
            _err = None if err is None else err * numpy.absolute(indxcorr)
        else:
            _indx = indx + indxcorr
            _err = None if err is None else err.copy()
        return _indx, _err


    @staticmethod
    def count_indices(absdb, bhddb):
        r"""
        Count the total number (absorption-line and bandhead) indices.
        """
        return 0 if absdb is None else absdb.nsets, 0 if bhddb is None else bhddb.nsets


    @staticmethod
    def output_dtype(nindx, bitmask=None):
        r"""
        Construct the record array data type for the output fits
        extension.
        """
        return [ ('BINID',numpy.int),
                 ('BINID_INDEX',numpy.int),
                 ('REDSHIFT', numpy.float),
                 ('MASK', numpy.bool if bitmask is None else bitmask.minimum_dtype(), (nindx,)),
                 ('BCEN', numpy.float, (nindx,)), 
                 ('BCONT', numpy.float, (nindx,)), 
                 ('BCONTERR', numpy.float, (nindx,)),
                 ('RCEN', numpy.float, (nindx,)), 
                 ('RCONT', numpy.float, (nindx,)), 
                 ('RCONTERR', numpy.float, (nindx,)), 
                 ('INDX', numpy.float, (nindx,)), 
                 ('INDXERR', numpy.float, (nindx,)),
                 ('MODEL_INDX', numpy.float, (nindx,)), 
                 ('INDX_DISPCORR', numpy.float, (nindx,))
               ]


    @staticmethod
    def check_and_prep_input(wave, flux, ivar=None, mask=None, redshift=None, bitmask=None):

        # Check the bitmask if provided
        if bitmask is not None and not isinstance(bitmask, BitMask):
            raise TypeError('Input bitmask must have type BitMask.')

        # Check the input wavelength and flux shapes
        if len(wave.shape) != 1:
            raise ValueError('Input wavelengths must be a single vector.')
        if len(wave) != flux.shape[1]:
            raise ValueError('Wavelength vector does not match shape of the flux array.')

        # Check the mask shape
        if mask is not None and mask.shape != flux.shape:
            raise ValueError('Input mask must have the same shape as the flux array.')

        # Check the input redshifts
        nspec = flux.shape[0]
        _redshift = numpy.zeros(nspec, dtype=numpy.float) if redshift is None else redshift
        if len(_redshift) != nspec:
            raise ValueError('Must provide one redshift per input spectrum (flux.shape[0]).')

        # Convert the input arrays to masked arrays if they aren't
        # already, and compare the array shapes
        _flux = flux if isinstance(flux, numpy.ma.MaskedArray) else \
                    numpy.ma.MaskedArray(flux, mask=(None if mask is None else mask > 0))

        if ivar is None:
            noise = None
        else:
            if ivar.shape != flux.shape:
                raise ValueError('Input ivar array must be the same shape as the flux array.')
            _ivar = ivar if isinstance(ivar, numpy.ma.MaskedArray) else \
                        numpy.ma.MaskedArray(ivar, mask=(None if mask is None else mask > 0))
            noise = numpy.ma.sqrt(1.0 /_ivar)

        return _flux, noise, _redshift


    @staticmethod
    def set_masks(measurements, blue_incomplete, blue_empty, red_incomplete, red_empty, divbyzero,
                  main_incomplete=None, main_empty=None, bitmask=None):

#        print('blue incomplete: {0}/{1}'.format(numpy.sum(blue_incomplete), blue_incomplete.size))
#        print('blue empty: {0}/{1}'.format(numpy.sum(blue_empty), blue_empty.size))
#        print('red incomplete: {0}/{1}'.format(numpy.sum(red_incomplete), red_incomplete.size))
#        print('red empty: {0}/{1}'.format(numpy.sum(red_empty), red_empty.size))
#        print('divbyzero: {0}/{1}'.format(numpy.sum(divbyzero), divbyzero.size))
#        if main_incomplete is not None:
#            print('main incomplete: {0}/{1}'.format(numpy.sum(main_incomplete),
#                                                    main_incomplete.size))
#        if main_empty is not None:
#            print('main empty: {0}/{1}'.format(numpy.sum(main_empty), main_empty.size))

        measurements['MASK'][blue_incomplete] = True if bitmask is None else \
                bitmask.turn_on(measurements['MASK'][blue_incomplete], 'BLUE_INCOMP')
        measurements['MASK'][blue_empty] = True if bitmask is None else \
                bitmask.turn_on(measurements['MASK'][blue_empty], 'BLUE_EMPTY')
        measurements['MASK'][red_incomplete] = True if bitmask is None else \
                bitmask.turn_on(measurements['MASK'][red_incomplete], 'RED_INCOMP')
        measurements['MASK'][red_empty] = True if bitmask is None else \
                bitmask.turn_on(measurements['MASK'][red_empty], 'RED_EMPTY')
        measurements['MASK'][divbyzero] = True if bitmask is None else \
                bitmask.turn_on(measurements['MASK'][divbyzero], 'DIVBYZERO')
        if main_incomplete is not None:
            measurements['MASK'][main_incomplete] = True if bitmask is None else \
                    bitmask.turn_on(measurements['MASK'][main_incomplete], 'MAIN_INCOMP')
        if main_empty is not None:
            measurements['MASK'][main_empty] = True if bitmask is None else \
                    bitmask.turn_on(measurements['MASK'][main_empty], 'MAIN_EMPTY')
        return measurements


    @staticmethod
    def save_results(results, measurements, good, err=False, bitmask=None):
        if not isinstance(results, (AbsorptionLineIndices, BandheadIndices)):
            raise TypeError('Input must be of type AbsorptionLineIndices or BandheadIndices')

        # Save the data
        measurements['BCEN'][good] = results.blue_center
        measurements['BCONT'][good] = results.blue_continuum
        if err:
            measurements['BCONTERR'][good] = results.blue_continuum_err
        measurements['RCEN'][good] = results.red_center
        measurements['RCONT'][good] = results.red_continuum
        if err:
            measurements['RCONTERR'][good] = results.red_continuum_err
        measurements['INDX'][good] = results.index
        if err:
            measurements['INDXERR'][good] = results.index_err

        # Reshape the flags
        blue_incomplete = numpy.zeros(len(good), dtype=numpy.bool)
        blue_incomplete[good] = results.blue_incomplete
        blue_empty = numpy.zeros(len(good), dtype=numpy.bool)
        blue_empty[good] = results.blue_empty
        red_incomplete = numpy.zeros(len(good), dtype=numpy.bool)
        red_incomplete[good] = results.red_incomplete
        red_empty = numpy.zeros(len(good), dtype=numpy.bool)
        red_empty[good] = results.red_empty
        main_incomplete = None
        if results.main_incomplete is not None:
            main_incomplete = numpy.zeros(len(good), dtype=numpy.bool)
            main_incomplete[good] = results.main_incomplete
        main_empty = None
        if results.main_empty is not None:
            main_empty = numpy.zeros(len(good),dtype=numpy.bool)
            main_empty[good] = results.main_empty
        divbyzero = numpy.zeros(len(good), dtype=numpy.bool)
        divbyzero[good] = results.divbyzero

        # Set the masks
        return SpectralIndices.set_masks(measurements, blue_incomplete, blue_empty, red_incomplete,
                                         red_empty, divbyzero, main_incomplete=main_incomplete,
                                         main_empty=main_empty, bitmask=bitmask)


    @staticmethod
    def measure_indices(absdb, bhddb, wave, flux, ivar=None, mask=None, redshift=None,
                        bitmask=None):
        r"""
        Measure the spectral indices in a set of spectra.

        Args:
            absdb
                (:class:`mangadap.par.aborptionlinedb.AbsorptionIndexDB`):
                Database with the absorption-line index definitions.
                Can be None.
            bhddb
                (:class:`mangadap.par.bandheadindexdb.BandheadIndexDB`):
                Database with the bandhead index definitions.
            wave (array-like): 1D vector with the wavelength of each
                pixel.  *Assumed to be logarithmically binned in
                radius.*
            flux (array-like): 2D array with the flux, ordered as
                :math:`N_{\rm spec}\times N_{\rm wave}`.  Can be a
                numpy.ma.MaskedArray; masked pixels will be ignored in
                the measurement.
            ivar (array-like): (**Optional**) Inverse variance in the
                flux.  Must match flux array shape.  Used to calculate
                propagated errors in the index.  Default is that errors
                are ignored.
            mask (array-like): (**Optional**) Boolean array flagging to
                ignore (mask=True) or include (mask=False) each flux
                measurement in the index calculation.
            redshift (array-like): (**Optional**) Redshift to use for
                each spectrum when determining the index.  Must have the
                correct length compared to the flux array.  Default is
                to assume the spectra are at rest wavelength.
            bitmask (:class:`mangadap.util.bitmask.BitMask`):
                (**Optional**)  If an index is flagged for some reason
                (see :func:`set_masks`), this object is used to set the
                mask value; this should typically be
                :class:`SpectralIndicesBitMask` object.  If not
                provided, the masked values are all set to True (False
                otherwise).

        Returns:
            `numpy.recarray`_: A record array with the following
            columns, each with one element per spectrum:

                 0. ``BINID``: Bin identifier
                 1. ``BINID_INDEX``: Index of the bin identifier
                 2. ``REDSHIFT``: Redshift used for measurement
                 3. ``MASK``: Boolean or maskbit value for index
                 4. ``BCEN``: Blue passband center
                 5. ``BCONT``: Blue passband pseudo-continuum
                 6. ``BCONTERR``: Error in the above
                 7. ``RCEN``: Red passband center
                 8. ``RCONT``: Red passband pseudo-continuum
                 9. ``RCONTERR``: Error in the above
                10. ``INDX``: Index value
                11. ``INDXERR``: Error in the above
                12. ``MODEL_INDX``: Index measured on the best-fitting
                                    stellar-continuum model
                13. ``INDX_DISPCORR``: Index dispersion correction

            This function does not add the ``BINID``, ``BINID_INDEX``,
            or ``INDX_DISPCORR`` values.  Each element in columns 3-12
            are vectors with a length of :math:`N_{\rm index}` --- the
            total number of indices calcualte (see
            :func:`count_indices`).

        """
        # Check the input databases
        if absdb is not None and not isinstance(absdb, AbsorptionIndexDB):
            raise TypeError('Input database must have type AbsorptionIndexDB.')
        if bhddb is not None and not isinstance(bhddb, BandheadIndexDB):
            raise TypeError('Input database must have type BandheadIndexDB.')

        # Get the number of indices
        nabs, nbhd = SpectralIndices.count_indices(absdb, bhddb)
        nindx = nabs+nbhd
        if nindx == 0:
            raise ValueError('No indices to measure!')

        nspec = flux.shape[0]
        # Check the input and initialize the output
        measurements = init_record_array(nspec, SpectralIndices.output_dtype(nindx,bitmask=bitmask))
        _flux, noise, measurements['REDSHIFT'] \
                    = SpectralIndices.check_and_prep_input(wave, flux, ivar=ivar, mask=mask,
                                                           redshift=redshift, bitmask=bitmask)
#        print(_flux.shape)
#        print(None if noise is None else noise.shape)
#        print(wave.shape)
#        print((numpy.array([wave]*nspec)).shape)

        # Create the f_nu spectra
        # The conversion is multiplicative, meaning the calculation of
        # the error calculation can use exactly the same function 
        flux_fnu = flux_to_fnu(numpy.array([wave]*nspec), _flux)
        noise_fnu = None if noise is None else flux_to_fnu(numpy.array([wave]*nspec), noise)

#        nu = 1e-13*astropy.constants.c.to('nm/s').value/wave    # In THz
#        pyplot.plot(wave, _flux[0,:]*wave, linestyle='-', color='r', lw=1, zorder=2)
#        pyplot.plot(wave, flux_fnu[0,:]*nu, linestyle='-', color='0.5', lw=3, zorder=1)
#        pyplot.show()

        # Get the list of good indices of each type
        abs_fnu = numpy.zeros(nabs, dtype=numpy.bool) if absdb is None \
                        else numpy.invert(absdb.dummy) & (absdb['integrand'] == 'fnu')
        good_abs_fnu = numpy.zeros(nindx, dtype=numpy.bool)
        if nabs > 0:
            good_abs_fnu[:nabs][abs_fnu] = True
        abs_flambda = numpy.zeros(nabs, dtype=numpy.bool) if absdb is None \
                        else numpy.invert(absdb.dummy) & (absdb['integrand'] == 'flambda')
        good_abs_flambda = numpy.zeros(nindx, dtype=numpy.bool)
        if nabs > 0:
            good_abs_flambda[:nabs][abs_flambda] = True

        bhd_fnu = numpy.zeros(nbhd, dtype=numpy.bool) if bhddb is None \
                        else numpy.invert(bhddb.dummy) & (bhddb['integrand'] == 'fnu')
        good_bhd_fnu = numpy.zeros(nindx, dtype=numpy.bool)
        if nbhd > 0:
            good_bhd_fnu[nabs:][bhd_fnu] = True
        bhd_flambda = numpy.zeros(nbhd, dtype=numpy.bool) if bhddb is None \
                        else numpy.invert(bhddb.dummy) & (bhddb['integrand'] == 'flambda')
        good_bhd_flambda = numpy.zeros(nindx, dtype=numpy.bool)
        if nbhd > 0:
            good_bhd_flambda[nabs:][bhd_flambda] = True

        # Mask any dummy indices
        dummy = numpy.zeros(nindx, dtype=numpy.bool)
        dummy[:nabs] = absdb.dummy
        dummy[nabs:] = bhddb.dummy
        if numpy.any(dummy):
            measurements['MASK'][:,dummy] = True if bitmask is None else \
                    bitmask.turn_on(measurements['MASK'][:,dummy], 'UNDEFINED_BANDS')

        # No valid indices
        if numpy.all(dummy):
            return measurements

        # Perform the measurements on each spectrum
        for i in range(nspec):

            print('Measuring spectral indices in spectrum: {0}/{1}'.format(i+1,nspec), end='\r')

            # -----------------------------------
            # Measure the absorption-line indices
            if nabs > 0:

                # Shift the bands
                _bluebands = absdb['blueside']*(1.0+measurements['REDSHIFT'][i])
                _redbands = absdb['redside']*(1.0+measurements['REDSHIFT'][i])
                _mainbands = absdb['primary']*(1.0+measurements['REDSHIFT'][i])

                # Integrate over F_nu
                if numpy.sum(good_abs_fnu) > 0:
                    # Make the measurements ...
                    _noise = None if noise_fnu is None else noise_fnu[i,:]
                    results = AbsorptionLineIndices(wave, flux_fnu[i,:], _bluebands[abs_fnu],
                                                    _redbands[abs_fnu], _mainbands[abs_fnu],
                                                    err=_noise, units=absdb['units'][abs_fnu])
                    # ... and save them
                    measurements[i] = SpectralIndices.save_results(results, measurements[i],
                                                                   good_abs_fnu,
                                                                   err=_noise is not None,
                                                                   bitmask=bitmask)

                # Integrate over F_lambda
                if numpy.sum(good_abs_flambda) > 0:
                    # Make the measurements ...
                    _noise = None if noise is None else noise[i,:]
                    results = AbsorptionLineIndices(wave, _flux[i,:], _bluebands[abs_flambda],
                                                    _redbands[abs_flambda], _mainbands[abs_flambda],
                                                    err=_noise, units=absdb['units'][abs_flambda])
                    # ... and save them
                    measurements[i] = SpectralIndices.save_results(results, measurements[i],
                                                                   good_abs_flambda,
                                                                   err=_noise is not None,
                                                                   bitmask=bitmask)

            # -----------------------------------
            # Measure the bandhead indices
            if nbhd > 0:

                # Shift the bands
                _bluebands = bhddb['blueside']*(1.0+measurements['REDSHIFT'][i])
                _redbands = bhddb['redside']*(1.0+measurements['REDSHIFT'][i])

                # Integrate over F_nu
                if numpy.sum(good_bhd_fnu) > 0:
                    # Make the measurements ...
                    _noise = None if noise_fnu is None else noise_fnu[i,:]
                    results = BandheadIndices(wave, flux_fnu[i,:], _bluebands[bhd_fnu],
                                              _redbands[bhd_fnu], err=_noise,
                                              order=bhddb['order'][bhd_fnu])
                    # ... and save them
                    measurements[i] = SpectralIndices.save_results(results, measurements[i],
                                                                   good_bhd_fnu,
                                                                   err=_noise is not None,
                                                                   bitmask=bitmask)

                # Integrate over F_lambda
                if numpy.sum(good_bhd_flambda) > 0:
                    # Make the measurements ...
                    _noise = None if noise is None else noise[i,:]
                    results = BandheadIndices(wave, _flux[i,:], _bluebands[bhd_flambda],
                                              _redbands[bhd_flambda], err=_noise,
                                              order=bhddb['order'][bhd_flambda])
                    # ... and save them
                    measurements[i] = SpectralIndices.save_results(results, measurements[i],
                                                                   good_bhd_flambda,
                                                                   err=_noise is not None,
                                                                   bitmask=bitmask)

#            pyplot.scatter(_interval, interval/_interval, marker='.', color='k', s=50)
#            pyplot.show()

#            center = numpy.mean(_sidebands, axis=1).reshape(-1,self.nindx)
#            pyplot.step(wave, flux[i,:], where='mid', color='k', lw=0.5, linestyle='-')
#            pyplot.scatter(center[0,:], measurements['BCONT'][i,:], marker='.',
#                           s=100, color='b', lw=0)
#            pyplot.scatter(center[1,:], measurements['RCONT'][i,:], marker='.',
#                           s=100, color='r', lw=0)
#            pyplot.show()

        print('Measuring spectral indices in spectrum: {0}/{0}'.format(nspec))

        # Correct the indices with angstrom units to rest-frame
        _, angu, _ = SpectralIndices.unit_selection(absdb, bhddb)
        measurements['INDX'][:,angu] /= (1+measurements['REDSHIFT'][:,None])

#        x = numpy.append(center, numpy.array([[-100]*self.absdb.nsets]).T,axis=1).ravel()
#        _x = numpy.ma.MaskedArray(x, mask=x==-100)
#        print(_x)
#        y = numpy.append(pseudocontinuum, numpy.array([[-100]*self.absdb.nsets]).T,axis=1).ravel()
#        _y = numpy.ma.MaskedArray(y, mask=y==-100)
#        print(type(_y))
#        print(numpy.sum(_y.mask))
#        pyplot.plot(_x, _y, color='g', lw=2, linestyle='-')
#        pyplot.show()

#        print('Total masked: {0}/{1}'.format(numpy.sum(measurements['MASK'] > 0),
#                                             measurements['MASK'].size))

        # Return the data
        return measurements


    def file_name(self):
        """Return the name of the output file."""
        return self.output_file


    def file_path(self):
        """Return the full path to the output file."""
        if self.directory_path is None or self.output_file is None:
            return None
        return os.path.join(self.directory_path, self.output_file)


    def measure(self, binned_spectra, redshift=None, stellar_continuum=None,
                emission_line_model=None, dapsrc=None, dapver=None, analysis_path=None,
                directory_path=None, output_file=None, hardcopy=True, tpl_symlink_dir=None,
                clobber=False, loggers=None, quiet=False):
        """
        Measure the spectral indices using the binned spectra and the
        internal spectral index database, and construct the internal
        data structure.

        If neither stellar-continuum nor emission-line models are
        provided:
            - Indices are measure on the binned spectra
            - No velocity-dispersion corrections are calculated

        If a stellar-continuum model is provided without an
        emission-line model:
            - Indices are measured on the binned spectra
            - Velocity-dispersion corrections are computed for any
              binned spectrum with a stellar-continuum fit based on the
              optimal template

        If an emission-line model is provided without a
        stellar-continuum model:
            - Indices are measured on the relevant (binned or unbinned)
              spectra; spectra with emission-line fits have the model
              emission lines subtracted from them before these
              measurements.
            - If the emission-line model includes data regarding the
              stellar-continuum fit (template spectra and template
              weights), corrections are calculated for spectra with
              emission-line models based on the continuum fits;
              otherwise, no corrections are calculated.

        If both stellar-continuum and emission-line models are provided,
        and if the stellar-continuum and emission-line fits are
        performed on the same spectra:
            - Indices are measured on the relevant (binned or unbinned)
              spectra; spectra with emission-line fits have the model
              emission lines subtracted from them before these
              measurements.
            - Velocity-dispersion corrections are based on the
              stellar-continuum templates and weights

        If both stellar-continuum and emission-line models are provided,
        and if the stellar-continuum and emission-line fits are
        performed on different spectra:
            - The behavior is exactly as if the stellar-continuum model
              was not provided.

        Args:
            binned_spectra
                (:class:`mangadap.proc.spatiallybinnedspectra.SpatiallyBinnedSpectra`):
                The binned spectra for the measurements.
            redshift (float, numpy.ndarray): (**Optional**) A single or
                spectrum-dependent redshift, :math:`z`, to use for
                shifting the index bands.  Default is to measure the
                indices at their provided wavelengths (i.e.,

                :math:`z=0`).  If providing spectrum-dependent values,
                the number of values must be the same as the number of
                stpectrum bins (i.e., binned_spectra.nbins) if either
                the emission-line model is not provided or it was not
                determined by deconstructing the bins; the number of
                values must be the same as the number of DRP spectra if
                the opposite is true (an emission-line model is provided
                that deconstructed the bins for its fit).


            stellar_continuum
                (:class:`mangadap.proc.stellarcontinuummodel.StellarContinuumModel`):
                (**Optional**) The stellar-continuum model as applied to
                the binned spectra.
            emission_line_model
                (:class:`mangadap.proc.emissionlinemodel.EmissionLineModel`):
                (**Optional**) The emission-line model as applied to
                either the binned spectra or the unbinned spaxels.
            dapsrc (str): (**Optional**) Root path to the DAP source
                directory.  If not provided, the default is defined by
                :func:`mangadap.config.defaults.dap_source_dir`.
            dapver (str): (**Optional**) The DAP version to use for the
                analysis, used to override the default defined by
                :func:`mangadap.config.defaults.default_dap_version`.
            analysis_path (str): (**Optional**) The top-level path for
                the DAP output files, used to override the default
                defined by
                :func:`mangadap.config.defaults.default_analysis_path`.
            directory_path (str): The exact path to the directory with
                DAP output that is common to number DAP "methods".  See
                :attr:`directory_path`.
            output_file (str): (**Optional**) Exact name for the output
                file.  The default is to use
                :func:`mangadap.config.defaults.default_dap_file_name`.
            hardcopy (bool): (**Optional**) Flag to write the HDUList
                attribute to disk.  Default is True; if False, the
                HDUList is only kept in memory and would have to be
                reconstructed.
            tpl_symlink_dir (str): (**Optional**) Create a symbolic link
                to the created template library file in the supplied
                directory.  Default is to produce no symbolic link.
            clobber (bool): (**Optional**) Overwrite any existing files.
                Default is to use any existing file instead of redoing
                the analysis and overwriting the existing output.
            loggers (list): (**Optional**) List of `logging.Logger`_
                objects to log progress; ignored if quiet=True.  Logging
                is done using :func:`mangadap.util.log.log_output`.
                Default is no logging.
            quiet (bool): (**Optional**) Suppress all terminal and
                logging output.  Default is False.
        """
        # Initialize the reporting
        if loggers is not None:
            self.loggers = loggers
        self.quiet = quiet

        # SpatiallyBinnedSpectra object always needed
        if binned_spectra is None:
            raise ValueError('Must provide spectra object for fitting.')
        if not isinstance(binned_spectra, SpatiallyBinnedSpectra):
            raise TypeError('Must provide a valid SpatiallyBinnedSpectra object!')
        if binned_spectra.hdu is None:
            raise ValueError('Provided SpatiallyBinnedSpectra object is undefined!')
        self.binned_spectra = binned_spectra

        # Check stellar-continuum model object, if provided
        self.stellar_continuum = None
        if stellar_continuum is not None:
            if not isinstance(stellar_continuum, StellarContinuumModel):
                raise TypeError('Provided stellar continuum must have StellarContinuumModel type!')
            if stellar_continuum.hdu is None:
                raise ValueError('Provided StellarContinuumModel is undefined!')
            self.stellar_continuum = stellar_continuum

        # Check emission-line model object, if provided
        self.emission_line_model = None
        if emission_line_model is not None:
            if not isinstance(emission_line_model, EmissionLineModel):
                raise TypeError('Provided emission line models must be of type EmissionLineModel.')
            if emission_line_model.hdu is None:
                raise ValueError('Provided EmissionLineModel is undefined!')
            self.emission_line_model = emission_line_model

        # What stellar continuum is available?
        #  - Assume the stellar continuum can always be extracted if a
        #    StellarContinuumModel object is provided
        #  - Check if the EmissionLineModel fitter has the appropriate
        #    function; True for Sasuke, not currently true for Elric
        eml_stellar_continuum_available = self.emission_line_model is not None \
                and callable(self.emission_line_model.method['fitclass'].construct_continuum_models)

        # Determine if the velocity-dispersion corrections can be
        # determined
        self.compute_corrections = self.database['compute_corrections']
        if self.compute_corrections and self.stellar_continuum is None \
                and not eml_stellar_continuum_available:
            warnings.warn('Cannot compute dispersion corrections; no continuum model available.')
            self.compute_corrections = False

        # Can only correct the indices if the corrections are provided
        if not self.compute_corrections and self.correct_indices:
            warnings.warn('Cannot apply corrections because they are not being computed.')
            self.correct_indices = False

        # What spectra to use?
        #  - Assume StellarContinuumModel always fits the binned spectra
        #  - The EmissionLineModel fits the binned spectra or unbinned
        #    spaxels as specified by its deconstruct_bins flag
        measure_on_unbinned_spaxels = self.emission_line_model is not None \
                and self.emission_line_model.method['deconstruct_bins']

        self.spatial_shape =self.binned_spectra.spatial_shape
        self.nspec = self.binned_spectra.drpf.nspec if measure_on_unbinned_spaxels \
                            else self.binned_spectra.nbins
        self.spatial_index = self.binned_spectra.spatial_index.copy()
        
        #---------------------------------------------------------------
        # Get the good spectra
        good_snr = self._flag_good_spectra(measure_on_unbinned_spaxels)

        # Set the number of bins measured and missing bins
        self.nbins = numpy.sum(good_snr)
        self.missing_bins = [] if measure_on_unbinned_spaxels else self._get_missing_bins()
        
        # Get the redshifts to apply
        self._assign_redshifts(redshift, measure_on_unbinned_spaxels, good_snr)

        # Report
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, '-'*50)
            log_output(loggers, 1, logging.INFO, '{0:^50}'.format('SPECTRAL-INDEX MEASUREMENTS'))
            log_output(self.loggers, 1, logging.INFO, '-'*50)
            log_output(self.loggers, 1, logging.INFO, 'Measurements for {0}'.format(
                        'unbinned spaxels' if measure_on_unbinned_spaxels else 'binned spectra'))
            log_output(self.loggers, 1, logging.INFO, 'Number of spectra: {0}'.format(self.nspec))
            if not measure_on_unbinned_spaxels and len(self.binned_spectra.missing_bins) > 0:
                log_output(self.loggers, 1, logging.INFO, 'Missing bins: {0}'.format(
                                                            len(self.binned_spectra.missing_bins)))
            log_output(self.loggers, 1, logging.INFO, 'With good S/N and to measure: {0}'.format(
                                                            numpy.sum(good_snr)))
            
        # Make sure there are good spectra
        if numpy.sum(good_snr) == 0:
            raise ValueError('No good spectra for measurements!')

        #---------------------------------------------------------------
        # (Re)Set the output paths
        self._set_paths(directory_path, dapver, analysis_path, output_file)

        #---------------------------------------------------------------
        # Check that the file path is defined
        ofile = self.file_path()
        if ofile is None:
            raise ValueError('File path for output file is undefined!')

        # Report
        if not self.quiet:
            log_output(self.loggers,1,logging.INFO,'Output path: {0}'.format(self.directory_path))
            log_output(self.loggers,1,logging.INFO,'Output file: {0}'.format(self.output_file))
        
        #---------------------------------------------------------------
        # If the file already exists, and not clobbering, just read the
        # file
        if os.path.isfile(ofile) and not clobber:
            self.hardcopy = True
            if not self.quiet:
                log_output(self.loggers, 1, logging.INFO, 'Using existing file')
            self.read(checksum=self.checksum)
            if not self.quiet:
                log_output(self.loggers, 1, logging.INFO, '-'*50)
            return

#        pyplot.scatter(numpy.arange(self.nbins), self.redshift, marker='.', s=50, color='k', lw=0)
#        pyplot.show()

        #---------------------------------------------------------------
        # Get the spectra to use for the measurements
        wave, flux, ivar = SpectralIndices.spectra_for_index_measurements(self.binned_spectra,
                                        measure_on_unbinned_spaxels=measure_on_unbinned_spaxels,
                                                    pixelmask=self.pixelmask, select=good_snr,
                                                    resolution_fwhm=self.database['fwhm'],
                                                    emission_line_model=self.emission_line_model)

        #---------------------------------------------------------------
        # Perform the measurements on the galaxy spectra
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO,
                       'Measuring spectral indices in observed spectra...')
        measurements = self.measure_indices(self.absdb, self.bhddb, wave, flux, ivar=ivar,
                                            redshift=self.redshift[good_snr], bitmask=self.bitmask)

#        print('Flagged D4000:')
#        for bit in ['MAIN_EMPTY', 'BLUE_EMPTY', 'RED_EMPTY', 'DIVBYZERO' ]:
#            print('measurements {0}: {1}/{2}'.format(bit,
#                        numpy.sum(self.bitmask.flagged(measurements['MASK'][:,43], flag=bit)),
#                        measurements['MASK'].shape[0]))

        #---------------------------------------------------------------
        # Determine the velocity dispersion corrections, if requested.
        if self.compute_corrections:

            # Get the template spectra to use
            replacement_templates = None if self.database['fwhm'] < 0 \
                    else self._resolution_matched_templates(dapsrc=dapsrc, dapver=dapver,
                                                            analysis_path=analysis_path,
                                                            tpl_symlink_dir=tpl_symlink_dir)
            # Have to use the corrected velocity dispersion if templates
            # have been broadened to a new resolution; otherwise, the
            # corrections are to a zero dispersion model **at the
            # resolution of the templates**
            corrected_dispersion = self.database['fwhm'] > 0

            # Set the bin IDs to match the stellar continuum to:
            binid = numpy.arange(binned_spectra.drpf.nspec).reshape(
                            binned_spectra.drpf.spatial_shape) if measure_on_unbinned_spaxels \
                                                               else binned_spectra['BINID'].data

            # Get two versions of the best-fitting continuum:
            #   - exactly the best fitting, continuum-only model
            #   - the same without convolving with the velocity
            #     dispersion
            # TODO: To also correct for difference in definition of
            # velocity, may also want to deredshift the model and
            # perform the measurements at rest wavelengths; not
            # currently possible though
            fill_to_match_f = self.emission_line_model.fill_continuum_to_match \
                                if eml_stellar_continuum_available else \
                                self.stellar_continuum.fill_to_match
            if not self.quiet:
                log_output(self.loggers, 1, logging.INFO, 'Constructing models with LOSVD')
            continuum = fill_to_match_f(binid, replacement_templates=replacement_templates,
                                        corrected_dispersion=corrected_dispersion) 
            if not self.quiet:
                log_output(self.loggers, 1, logging.INFO, 'Constructing models without LOSVD')
            continuum_dcnvlv = fill_to_match_f(binid, replacement_templates=replacement_templates,
                                               redshift_only=True) #, deredshift=True)

#            pyplot.imshow(flux, origin='lower', interpolation='nearest', aspect='auto')
#            pyplot.colorbar()
#            pyplot.show()
#
#            pyplot.imshow(continuum[good_snr,:], origin='lower', interpolation='nearest',
#                          aspect='auto')
#            pyplot.colorbar()
#            pyplot.show()
#
#            pyplot.imshow(continuum_dcnvlv[good_snr,:], origin='lower', interpolation='nearest',
#                          aspect='auto')
#            pyplot.colorbar()
#            pyplot.show()
#
#            indx = numpy.argmax(numpy.ma.mean(flux, axis=1))
#            print(indx)
#            pyplot.plot(wave, flux[indx,:])
#            pyplot.plot(wave, continuum[good_snr,:][indx,:])
#            pyplot.plot(wave, continuum_dcnvlv[good_snr,:][indx,:])
#            pyplot.show()

            # Get the corrections by performing the measurements on the
            # best-fitting continuum models, with and without the
            # velocity dispersion broadening
            if not self.quiet:
                log_output(self.loggers, 1, logging.INFO,
                           'Calculating dispersion corrections using stellar continuum model...')
            measurements['MODEL_INDX'], measurements['INDX_DISPCORR'], \
                    good_les, good_ang, good_mag \
                        = SpectralIndices.calculate_dispersion_corrections(self.absdb, self.bhddb,
                                        wave, flux, continuum[good_snr,:],
                                        continuum_dcnvlv[good_snr,:],
                                        redshift=self.redshift[good_snr],
                                        redshift_dcnvlv=self.redshift[good_snr],
                                        bitmask=self.bitmask)

            # Flag bad corrections
            bad_correction = numpy.invert(good_les) & numpy.invert(good_ang) \
                                & numpy.invert(good_mag)
            measurements['MASK'][bad_correction] = self.bitmask.turn_on(
                                                            measurements['MASK'][bad_correction],
                                                            'NO_DISPERSION_CORRECTION')
            # Apply the corrections
            if self.correct_indices:
                # Measured indices
                measurements['INDX'][good_les], measurements['INDXERR'][good_les] \
                        = SpectralIndices.apply_dispersion_corrections(
                                                    measurements['INDX'][good_les],
                                                    measurements['INDX_DISPCORR'][good_les],
                                                    err=measurements['INDXERR'][good_les])
                measurements['INDX'][good_ang], measurements['INDXERR'][good_ang] \
                        = SpectralIndices.apply_dispersion_corrections(
                                                    measurements['INDX'][good_ang],
                                                    measurements['INDX_DISPCORR'][good_ang],
                                                    err=measurements['INDXERR'][good_ang],
                                                    unit='ang')
                measurements['INDX'][good_mag], measurements['INDXERR'][good_mag] \
                        = SpectralIndices.apply_dispersion_corrections(
                                                    measurements['INDX'][good_mag],
                                                    measurements['INDX_DISPCORR'][good_mag],
                                                    err=measurements['INDXERR'][good_mag],
                                                    unit='mag')
                # Model indices
                measurements['MODEL_INDX'][good_les], _ \
                        = SpectralIndices.apply_dispersion_corrections(
                                                    measurements['MODEL_INDX'][good_les],
                                                    measurements['INDX_DISPCORR'][good_les])
                measurements['MODEL_INDX'][good_ang], _ \
                        = SpectralIndices.apply_dispersion_corrections(
                                                    measurements['MODEL_INDX'][good_ang],
                                                    measurements['INDX_DISPCORR'][good_ang],
                                                    unit='ang')
                measurements['MODEL_INDX'][good_mag], _ \
                        = SpectralIndices.apply_dispersion_corrections(
                                                    measurements['MODEL_INDX'][good_mag],
                                                    measurements['INDX_DISPCORR'][good_mag],
                                                    unit='mag')

        #---------------------------------------------------------------
        # Set the number of bins measured, missing bins, and bin IDs
        measurements['BINID'] = numpy.arange(self.nbins) if measure_on_unbinned_spaxels \
                                    else self.binned_spectra['BINS'].data['BINID'][good_snr]
        measurements['BINID_INDEX'] = numpy.arange(self.nbins) if measure_on_unbinned_spaxels \
                                        else numpy.arange(self.binned_spectra.nbins)[good_snr]
        measurements_binid = None
        if measure_on_unbinned_spaxels:
            measurements_binid = numpy.full(self.binned_spectra.spatial_shape, -1, dtype=int)
            measurements_binid.ravel()[good_snr] = numpy.arange(self.nbins)

        #---------------------------------------------------------------
        # Initialize the header keywords
        self.hardcopy = hardcopy
        pri_hdr = self._initialize_primary_header(measurements_binid=measurements_binid)
        map_hdr = DAPFitsUtil.build_map_header(self.binned_spectra.drpf,
                                               'K Westfall <westfall@ucolick.org>')
        # Get the spatial map mask
        map_mask = numpy.zeros(self.spatial_shape, dtype=self.bitmask.minimum_dtype())

        # Account for measurements on individual spaxels
        if measurements_binid is None:
            # Add any spaxel not used because it was flagged by the
            # binning step
            indx = self.binned_spectra['MAPMASK'].data > 0
            map_mask[indx] = self.bitmask.turn_on(map_mask[indx], 'DIDNOTUSE')
            # Isolate any spaxels with foreground stars
            indx = self.binned_spectra.bitmask.flagged(self.binned_spectra['MAPMASK'].data,
                                                       'FORESTAR')
            map_mask[indx] = self.bitmask.turn_on(map_mask[indx], 'FORESTAR')
            # Get the bins that were below the S/N limit
            indx = numpy.invert(DAPFitsUtil.reconstruct_map(self.spatial_shape,
                                                        self.binned_spectra['BINID'].data.ravel(),
                                                        good_snr, dtype='bool')) & (map_mask == 0)
            map_mask[indx] = self.bitmask.turn_on(map_mask[indx], 'LOW_SNR')

            # Get the bin ids with measured indices
            bin_indx = DAPFitsUtil.downselect_bins(self.binned_spectra['BINID'].data.ravel(),
                                                measurements['BINID']).reshape(self.spatial_shape)
        else:
            # Assume any model with a binid less than zero is from a
            # spaxel that was not used
            indx = measurements_binid < 0
            map_mask[indx] = self.bitmask.turn_on(map_mask[indx], 'DIDNOTUSE')

            # The number of valid bins MUST match the number of
            # measurements
            nvalid = numpy.sum(numpy.invert(indx))
            if nvalid != len(measurements):
                raise ValueError('Provided id does not match the number of measurements.')

            # Get the bin ids with fitted models
            bin_indx = measurements_binid

        # Compile the information on the suite of measured indices
        passband_database = self._compile_database()

        # Save the data to the hdu attribute
        self.hdu = fits.HDUList([ fits.PrimaryHDU(header=pri_hdr),
                                  fits.ImageHDU(data=bin_indx, header=map_hdr, name='BINID'),
                                  fits.ImageHDU(data=map_mask, header=map_hdr, name='MAPMASK'),
                                  fits.BinTableHDU.from_columns( [ fits.Column(name=n,
                                                    format=rec_to_fits_type(passband_database[n]),
                            array=passband_database[n]) for n in passband_database.dtype.names ],
                                                               name='SIPAR'),
                                  fits.BinTableHDU.from_columns( [ fits.Column(name=n,
                                                    format=rec_to_fits_type(measurements[n]),
                                array=measurements[n]) for n in measurements.dtype.names ],
                                                               name='SINDX')
                                ])

        # Write the data, if requested
        if self.hardcopy:
            if not os.path.isdir(self.directory_path):
                os.makedirs(self.directory_path)
            self.write(clobber=clobber)
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, '-'*50)


    def write(self, clobber=False):
        """
        Write the hdu object to the file.
        """
        DAPFitsUtil.write(self.hdu, self.file_path(), clobber=clobber, checksum=True,
                          loggers=self.loggers, quiet=self.quiet)


    def read(self, ifile=None, strict=True, checksum=False):
        """
        Read an existing file with a previously binned set of spectra.
        """
        if ifile is None:
            ifile = self.file_path()
        if not os.path.isfile(ifile):
            raise FileNotFoundError('File does not exist!: {0}'.format(ifile))

        if self.hdu is not None:
            self.hdu.close()

#        self.hdu = fits.open(ifile, checksum=checksum)
        self.hdu = DAPFitsUtil.read(ifile, checksum=checksum)

        # Confirm that the internal method is the same as the method
        # that was used in writing the file
        if self.hdu['PRIMARY'].header['SIKEY'] != self.database['key']:
            if strict:
                raise ValueError('Keywords in header does not match specified method keyword!')
            else:
                warnings.warn('Keywords in header does not match specified method keyword!')
        # TODO: "strict" should also check other aspects of the file to
        # make sure that the details of the method are also the same,
        # not just the keyword

#        if not self.quiet:
#            log_output(self.loggers, 1, logging.INFO, 'Reverting to python-native structure.')
#        DAPFitsUtil.restructure_map(self.hdu, ext=self.image_arrays)

        self.nbins = self.hdu['PRIMARY'].header['NBINS']
        unique_bins = numpy.unique(self.hdu['BINID'].data.ravel()) \
                            if self.hdu['PRIMARY'].header['SIREBIN'] else None
        self.missing_bins = self._get_missing_bins(unique_bins=unique_bins)


