# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
A class hierarchy that performs the spectral-index measurements.

Notes
-----

If neither stellar-continuum nor emission-line models are provided:

    - Indices are measure on the binned spectra
    - No velocity-dispersion corrections are calculated

If a stellar-continuum model is provided without an emission-line
model:

    - Indices are measured on the binned spectra
    - Velocity-dispersion corrections are computed for any binned
      spectrum with a stellar-continuum fit based on the optimal
      template

If an emission-line model is provided without a stellar-continuum model:

    - Indices are measured on the relevant (binned or unbinned) spectra;
      spectra with emission-line fits have the model emission lines
      subtracted from them before these measurements.
    - If the emission-line model includes data regarding the
      stellar-continuum fit (template spectra and template weights),
      corrections are calculated for spectra with emission-line models
      based on the continuum fits; otherwise, no corrections are
      calculated.

If both stellar-continuum and emission-line models are provided, and if
the stellar-continuum and emission-line fits are performed on the same
spectra:

    - Indices are measured on the relevant (binned or unbinned) spectra;
      spectra with emission-line fits have the model emission lines
      subtracted from them before these measurements.
    - Velocity-dispersion corrections are based on the stellar-continuum
      templates and weights

If both stellar-continuum and emission-line models are provided, and if
the stellar-continuum and emission-line fits are performed on different
spectra:

    - The behavior is exactly as if the stellar-continuum model was not
      provided.

----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

import glob
import os
import logging
import warnings

from IPython import embed

import numpy
from astropy.io import fits
import astropy.constants

from ..par.parset import KeywordParSet
from ..par.artifactdb import ArtifactDB
from ..par.absorptionindexdb import AbsorptionIndexDB
from ..par.bandheadindexdb import BandheadIndexDB
from ..config import defaults
from ..util.datatable import DataTable
from ..util.resolution import SpectralResolution, match_spectral_resolution
from ..util.sampling import spectral_coordinate_step, spectrum_velocity_scale
from ..util.fitsutil import DAPFitsUtil
from ..util.fileio import init_record_array, rec_to_fits_type
from ..util.log import log_output
from ..util.bitmask import BitMask
from ..util.dapbitmask import DAPBitMask
from ..util.pixelmask import SpectralPixelMask
from ..util.parser import DefaultConfig
from .spatiallybinnedspectra import SpatiallyBinnedSpectra
from .templatelibrary import TemplateLibrary
from .stellarcontinuummodel import StellarContinuumModel
from .emissionlinemodel import EmissionLineModel
from .bandpassfilter import passband_integral, passband_integrated_width, pseudocontinuum
from . import util

from matplotlib import pyplot


class SpectralIndicesDef(KeywordParSet):
    """
    A class that holds the parameters necessary to perform the
    spectral-index measurements.

    The defined parameters are:

    .. include:: ../tables/spectralindicesdef.rst
    """
    def __init__(self, key=None, minimum_snr=None, fwhm=None, compute_corrections=None,
                 artifacts=None, absindex=None, bandhead=None):
        in_fl = [ int, float ]

        pars =     [ 'key', 'minimum_snr', 'fwhm', 'compute_corrections', 'artifacts', 'absindex',
                        'bandhead' ]
        values =   [   key,   minimum_snr,   fwhm,   compute_corrections,   artifacts,   absindex,
                          bandhead ]
        dtypes =   [   str,         in_fl,  in_fl,                  bool,         str,        str,
                               str ]
        descr = ['Keyword used to distinguish between different spectral-index databases.',
                 'Minimum S/N of spectrum to fit',
                 'Resolution FWHM in angstroms at which to make the measurements.',
                 'Flag to compute dispersion corrections to indices.  Dispersion corrections ' \
                    'are always calculated!',
                 'String identifying the artifact database to use',
                 'String identifying the absorption-index database to use',
                 'String identifying the bandhead-index database to use']

        super(SpectralIndicesDef, self).__init__(pars, values=values, dtypes=dtypes, descr=descr)


def validate_spectral_indices_config(cnfg):
    """ 
    Validate :class:`~mangadap.util.parser.DefaultConfig` object with
    spectral-index measurement parameters.

    Args:
        cnfg (:class:`~mangadap.util.parser.DefaultConfig`):
            Object with parameters to validate.

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


def available_spectral_index_databases():
    """
    Return the list of available spectral index databases

    Available database combinations:

    .. todo::
        Fill in

    Returns:
        list: A list of :func:`SpectralIndicesDef` objects, each
        defining a spectral-index database to measure.

    Raises:
        IOError:
            Raised if no spectral-index configuration files could be
            found.
        KeyError:
            Raised if the spectral-index database keywords are not
            all unique.

    .. todo::
        - Somehow add a python call that reads the databases and
          constructs the table for presentation in sphinx so that the
          text above doesn't have to be edited with changes in the
          available databases.
        
    """
    # Check the configuration files exist
    search_dir = os.path.join(defaults.dap_config_root(), 'spectral_indices')
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

        index_set_list += [ SpectralIndicesDef(key=cnfg['key'],
                                               minimum_snr=cnfg.getfloat('minimum_snr', default=0.), 
                                               fwhm=cnfg.getfloat('resolution_fwhm', default=-1),
                                    compute_corrections=cnfg.getbool('compute_sigma_correction',
                                                                     default=False),
                                               artifacts=cnfg['artifact_mask'],
                                               absindex=cnfg['absorption_indices'],
                                               bandhead=cnfg['bandhead_indices']) ]

    # Check the keywords of the libraries are all unique
    if len(numpy.unique(numpy.array([index['key'] for index in index_set_list]))) \
                != len(index_set_list):
        raise KeyError('Spectral-index set keywords are not all unique!')

    # Return the default list of assessment methods
    return index_set_list


class SpectralIndicesBitMask(DAPBitMask):
    r"""
    Derived class that specifies the mask bits for the spectral-index
    measurements. The maskbits defined are:

    .. include:: ../tables/spectralindicesbitmask.rst
    """
    cfg_root = 'spectral_indices_bits'


class SpectralIndicesDefinitionTable(DataTable):
    """
    Table with the definitions of the spectral indices.

    Table includes:

    .. include:: ../tables/spectralindicesdefinitiontable.rst

    Args:
        name_len (:obj:`int`):
            The maximum length of any of the spectral index names.
        shape (:obj:`int`, :obj:`tuple`, optional):
            The shape of the initial array. If None, the data array
            will not be instantiated; use :func:`init` to initialize
            the data array after instantiation.
    """
    def __init__(self, name_len=1, shape=None):
        # NOTE: This should require python 3.7 to make sure that this
        # is an "ordered" dictionary.
        datamodel = dict(TYPE=dict(typ='<U10', shape=None,
                                   descr='Type of spectral index, either absorption or bandhead.'),
                         ID=dict(typ=int, shape=None, descr='ID number for the spectral index'),
                         NAME=dict(typ='<U{0:d}'.format(name_len), shape=None,
                                   descr='Unique name for the spectral index'),
                         PASSBAND=dict(typ=float, shape=(2,),
                                       descr='Lower and upper rest wavelength of the main index '
                                             'passband (absorption-line indices only)'),
                         BLUEBAND=dict(typ=float, shape=(2,),
                                       descr='Lower and upper rest wavelength of the blue '
                                             'sideband used to define the linear continuum'),
                         REDBAND=dict(typ=float, shape=(2,),
                                       descr='Lower and upper rest wavelength of the red '
                                             'sideband used to define the linear continuum'),
                         UNIT=dict(typ='<U3', shape=None, descr='Index units'),
                         COMPONENT=dict(typ=bool, shape=None,
                                        descr='Flag if the index is a component of a '
                                              'multi-component index (ignored)'),
                         INTEGRAND=dict(typ='<U7', shape=None,
                                        descr='The flux integrand for the index.'),
                         ORDER=dict(typ='<U3', shape=None,
                                    descr='The numerator-denominator order of the index '
                                          '(bandhead/color indices only)'))

        keys = list(datamodel.keys())
        super(SpectralIndicesDefinitionTable,
                self).__init__(keys, [datamodel[k]['typ'] for k in keys],
                               element_shapes=[datamodel[k]['shape'] for k in keys],
                               descr=[datamodel[k]['descr'] for k in keys],
                               shape=shape) 


class SpectralIndicesDataTable(DataTable):
    """
    Primary data table with the results of the spectral-index
    measurements.

    Table includes:

    .. include:: ../tables/spectralindicesdatatable.rst

    Args:
        nindx (:obj:`int`):
            Number of spectral indices to measure.
        bitmask (:class:`~mangadap.util.bitmask.BitMask`, optional):
            Object used to flag mask bits. If None, flags are simply
            boolean.
        shape (:obj:`int`, :obj:`tuple`, optional):
            The shape of the initial array. If None, the data array
            will not be instantiated; use :func:`init` to initialize
            the data array after instantiation.
    """
    def __init__(self, nindx=1, bitmask=None, shape=None):
        # NOTE: This should require python 3.7 to make sure that this
        # is an "ordered" dictionary.
        datamodel = dict(BINID=dict(typ=int, shape=None, descr='Spectrum/Bin ID number'),
                         BINID_INDEX=dict(typ=int, shape=None,
                                          descr='Index of the spectrum in the list of '
                                                'provided spectra.'),
                         REDSHIFT=dict(typ=float, shape=None,
                                       descr='Redshift used for shifting the passbands'),
                         MASK=dict(typ=bool if bitmask is None else bitmask.minimum_dtype(),
                                   shape=(nindx,),
                                   descr='Bad-value boolean or bit mask value for the moments'),
                         BCEN=dict(typ=float, shape=(nindx,), descr='Center of the blue sideband'),
                         BCONT=dict(typ=float, shape=(nindx,),
                                    descr='Pseudo-continuum in the blue sideband'),
                         BCONT_ERR=dict(typ=float, shape=(nindx,),
                                        descr='Error in the blue-sideband pseudo-continuum'),
                         BCONT_MOD=dict(typ=float, shape=(nindx,),
                                        descr='Pseudo-continuum in the blue sideband of the '
                                              'best-fitting model spectrum'),
                         BCONT_CORR=dict(typ=float, shape=(nindx,),
                                         descr='Multiplicative correction to apply to the '
                                               'pseudo-continuum in the blue sideband to match '
                                               'the measurement for a spectrum with no Doppler '
                                               'broadening'),
                         RCEN=dict(typ=float, shape=(nindx,), descr='Center of the red sideband'),
                         RCONT=dict(typ=float, shape=(nindx,),
                                    descr='Pseudo-continuum in the red sideband'),
                         RCONT_ERR=dict(typ=float, shape=(nindx,),
                                        descr='Error in the red-sideband pseudo-continuum'),
                         RCONT_MOD=dict(typ=float, shape=(nindx,),
                                        descr='Pseudo-continuum in the red sideband of the '
                                              'best-fitting model spectrum'),
                         RCONT_CORR=dict(typ=float, shape=(nindx,),
                                         descr='Multiplicative correction to apply to the '
                                               'pseudo-continuum in the red sideband to match '
                                               'the measurement for a spectrum with no Doppler '
                                               'broadening'),
                         MCONT=dict(typ=float, shape=(nindx,),
                                    descr='Interpolated continuum at the center of the main index '
                                          'passband (absorption-line indices only)'),
                         MCONT_ERR=dict(typ=float, shape=(nindx,),
                                        descr='Error in the continuum of the main passband'),
                         MCONT_MOD=dict(typ=float, shape=(nindx,),
                                        descr='Continuum in the main passband of the best-fitting '
                                              'model spectrum'),
                         MCONT_CORR=dict(typ=float, shape=(nindx,),
                                         descr='Multiplicative correction to apply to the '
                                               'continuum in the main passband to match the '
                                               'measurement for a spectrum with no Doppler '
                                               'broadening'),
                         AWGT=dict(typ=float, shape=(nindx,),
                                   descr='Weight needed to construct an index measured for the '
                                         'sum (or mean) of a set of spectra.  This identical to '
                                         'data pulled from either BCONT, RCONT, or MCONT, as '
                                         'appropriate to the calculation.'),
                         AWGT_ERR=dict(typ=float, shape=(nindx,),
                                        descr='Error in the aggregation weight'),
                         AWGT_MOD=dict(typ=float, shape=(nindx,),
                                        descr='Weight as measured using the best-fitting '
                                              'model spectrum'),
                         AWGT_CORR=dict(typ=float, shape=(nindx,),
                                         descr='Multiplicative correction to apply to the '
                                               'weight needed to match the measurement for a '
                                               'spectrum with no Doppler broadening'),
                         INDX=dict(typ=float, shape=(nindx,),
                                   descr='Index measurement for the observed spectrum.  '
                                         'Absorption lines are measured using the definition in '
                                         'Worthey et al. (1994)'),
                         INDX_ERR=dict(typ=float, shape=(nindx,),
                                        descr='Error in the spectral index'),
                         INDX_MOD=dict(typ=float, shape=(nindx,),
                                        descr='Spectral index as measured using the best-fitting '
                                              'model spectrum'),
                         INDX_CORR=dict(typ=float, shape=(nindx,),
                                         descr='Index correction needed to match the measurement '
                                               'for a spectrum with no Doppler broadening'),
                         INDX_BF=dict(typ=float, shape=(nindx,),
                                      descr='Index measurement for the observed spectrum.  '
                                            'Absorption lines are measured using the definition '
                                            'used by Burstein et al. (1984) and Faber et al. '
                                            '(1985); see also Faber et al. (1977)'),
                         INDX_BF_ERR=dict(typ=float, shape=(nindx,),
                                        descr='Error in the spectral index'),
                         INDX_BF_MOD=dict(typ=float, shape=(nindx,),
                                        descr='Spectral index as measured using the best-fitting '
                                              'model spectrum'),
                         INDX_BF_CORR=dict(typ=float, shape=(nindx,),
                                         descr='Index correction needed to match the measurement '
                                               'for a spectrum with no Doppler broadening'))

        keys = list(datamodel.keys())
        super(SpectralIndicesDataTable,
                self).__init__(keys, [datamodel[k]['typ'] for k in keys],
                               element_shapes=[datamodel[k]['shape'] for k in keys],
                               descr=[datamodel[k]['descr'] for k in keys],
                               shape=shape)


# TODO: SpecralIndices.save_results requires that AbsorptionLineIndices
# and BandheadIndices have the following attributes in common:
#        blue_center
#        blue_continuum
#        blue_continuum_err
#        blue_incomplete
#        blue_empty
#
#        red_center
#        red_continuum
#        red_continuum_err
#        red_incomplete
#        red_empty
#
#        index
#        index_err
#
#        main_incomplete
#        main_empty
#        divbyzero
#
# Put these in a common base class?


class AbsorptionLineIndices:
    r"""
    Measure a set of absorption-line indices and metrics in a single
    spectrum.

    .. include:: ../include/absindices.rst

    .. note::

        The calculations performed by this class are agnostic as to
        the redshift of the spectrum provided. It is expected that
        you will have either de-redshifted the spectrum or redshifted
        the passband definitions appropriately to measure the index
        in a redshifted spectrum. However, note that, for the latter
        case, the calculated indices with angstrom units should be
        divided by :math:`1+z` (where :math:`z` is the redshift) to
        calculate the *rest-frame* spectral index; no factor is
        needed for the indices in magnitude units.
    
    Args:
        wave (`numpy.ndarray`_):
            Wavelength vector.
        flux (array-like):
            Flux vector. Units must be appropriate to the definition
            of the provided index passbands. Can be a
            `numpy.ma.MaskedArray`_; masked pixels will be ignored in
            the measurement.
        bluebands (`numpy.ndarray`_):
            Array with the definition of the blue sidebands. Shape
            must be :math:`(N_{\rm index},2)`, where the second axis
            provides the starting and ending wavelength of each band.
        redbands (`numpy.ndarray`_):
            Array with the definition of the red sidebands. Shape
            must be :math:`(N_{\rm index},2)`, where the second axis
            provides the starting and ending wavelength of each band.
        mainbands (`numpy.ndarray`_):
            Array with the definition of the primary passbands. Shape
            must be :math:`(N_{\rm index},2)`, where the second axis
            provides the starting and ending wavelength of each band.
        err (`numpy.ndarray`_, optional):
            Error in the flux, used for a nominal propagation of the
            index error.
        log (`numpy.ndarray`_, optional):
            Flag that the provided flux vector is sampled
            logarithmically in wavelength.
        units (`numpy.ndarray`_, optional):
            Array with the units for each index to compute. Elements
            of the array must be either ``'ang'`` or ``'mag'``. If
            ``units`` is None, ``'ang'`` is assumed for all indices.

    Attributes:
        nindx (:obj:`int`):
            Number of indices
        units (`numpy.ndarray`_):
            String array with the units of each index (``'ang'`` or
            ``'mag'``).
        blue_center (`numpy.ndarray`_):
            Center of the blue sideband.
        blue_continuum (`numpy.ndarray`_):
            Flux in the blue sideband.
        blue_continuum_err (`numpy.ndarray`_):
            Propagated error in the blue sideband flux.
        blue_incomplete (`numpy.ndarray`_):
            Boolean array flagging if the blue sideband contained
            *any* masked pixels.
        blue_empty (`numpy.ndarray`_):
            Boolean array flagging if the blue sideband was
            completely empty.
        red_center (`numpy.ndarray`_):
            Center of the red sideband.
        red_continuum (`numpy.ndarray`_):
            Flux in the red sideband.
        red_continuum_err (`numpy.ndarray`_):
            Propagated error in the red sideband flux.
        red_incomplete (`numpy.ndarray`_):
            Boolean array flagging if the red sideband contained
            *any* masked pixels.
        red_empty (`numpy.ndarray`_):
            Boolean array flagging if the red sideband was completely
            empty.
        continuum_m (`numpy.ndarray`_):
            Slope of the linear continuum defined by the two
            sidebands.
        continuum_b (`numpy.ndarray`_):
            Intercept of the linear continuum defined by the two
            sidebands.
        main_continuum (`numpy.ndarray`_):
            The continuum interpolated to the center of the main
            passband based on the linear continuum constructed using
            the blue and red sidebands.
        main_continuum_err (`numpy.ndarray`_):
            The propagated error in the continuum interpolated to the
            center of the main passband.
        main_flux (`numpy.ndarray`_):
            The integral of the flux over the main passband.
        main_flux_err (`numpy.ndarray`_):
            Propagated error in the integral of the flux over the
            main passband.
        main_incomplete (`numpy.ndarray`_):
            Boolean array flagging if the main passband contained
            *any* masked pixels.
        main_empty (`numpy.ndarray`_):
            Boolean array flagging if the main passband was
            completely empty.
        index (`numpy.ndarray`_):
            Computed index following the Worthey et al. definition.
        index_err (`numpy.ndarray`_):
            Error in the WT index from nominal error propagation.
            Does *not* include error in the linear continuum.
        index_bf (`numpy.ndarray`_):
            Computed index following the Burstein et al. definition.
        index_bf_err (`numpy.ndarray`_):
            Error in the BF index from nominal error propagation.
            Does *not* include error in the linear continuum.
        divbyzero (`numpy.ndarray`_):
            Boolean array flagging if the index computation
            encountered a division by 0 error; only flagged for
            indices with magnitude units.
    """
    def __init__(self, wave, flux, bluebands, redbands, mainbands, err=None, log=True, units=None):

#        weighted_center (:obj:`bool`, optional):
#            For the construction of the linear continuum, compute the
#            flux-weighted center of the two side bands. If False, the
#            unweighted center is used (i.e., :math:`(\lambda_1 +
#            \lambda_2)/2`, where :math:`\lambda_1, \lambda_2` define
#            the edges of the sideband.)
#                 weighted_center=True):

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

        # Get the two pseudocontinua and the bandpass centers.
        # `pseudocontinuum` returns center, continuum, and
        # continuum_err as MaskedArrays. The optional keyword for the
        # center calculation defaults to False, meaning the band
        # centers should just be the unweighted center.
        self.blue_center, self.blue_continuum, self.blue_continuum_err, self.blue_incomplete, \
                self.blue_empty = pseudocontinuum(wave, flux, passband=bluebands, err=err, log=log)

        self.red_center, self.red_continuum, self.red_continuum_err, self.red_incomplete, \
                self.red_empty = pseudocontinuum(wave, flux, passband=redbands, err=err, log=log)

        # Get the parameters for the linear continuum across the
        # primary passband
        self.continuum_m = (self.red_continuum - self.blue_continuum) \
                                / (self.red_center - self.blue_center)
        self.continuum_b = self.blue_continuum - self.blue_center * self.continuum_m

        # Calculate the continuum at the center of the main band
        self.main_center = numpy.mean(mainbands, axis=1)
        dx = (self.main_center - self.blue_center) / (self.red_center - self.blue_center)
        self.main_continuum = dx * self.red_continuum + (1-dx) * self.blue_continuum
        # Above is the same as:
        #   self.main_continuum = self.continuum_m * self.main_center + self.continuum_b
        self.main_continuum_err = None if err is None else \
                                    numpy.sqrt(numpy.square(dx * self.red_continuum_err) 
                                               + numpy.square((1-dx) * self.blue_continuum_err))

#        print('Number of non-finite continuum m,b: {0} {1}'.format(
#                                numpy.sum(numpy.invert(numpy.isfinite(self.continuum_m))),
#                                numpy.sum(numpy.invert(numpy.isfinite(self.continuum_b)))))

        # Compute the continuum normalized indices.  This has to be done
        # in a for loop because the continuum is index-dependent
        self.main_flux = numpy.zeros(self.nindx, dtype=float)
        self.main_flux_err = numpy.zeros(self.nindx, dtype=float)
        self.main_incomplete = numpy.zeros(self.nindx, dtype=bool)
        self.main_empty = numpy.zeros(self.nindx, dtype=bool)
        self.index = numpy.zeros(self.nindx, dtype=float)
        self.index_err = numpy.zeros(self.nindx, dtype=float)
        self.index_bf = numpy.zeros(self.nindx, dtype=float)
        self.index_bf_err = numpy.zeros(self.nindx, dtype=float)
        # TODO: Need two of these? One for each of the index definitions?
        self.divbyzero = numpy.zeros(self.nindx, dtype=bool)

        # For convenience:
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
                self.index_err[i] = passband_integral(wave, err/cont, passband=m, log=log,
                                                      quad=True)

            # Calculate the integrated flux over passband
            self.main_flux[i] = passband_integral(wave, flux, passband=m, log=log)
            if err is not None:
                self.main_flux_err[i] = passband_integral(wave, err, passband=m, log=log,
                                                          quad=True)

            # Get the fraction of the band covered by the spectrum and
            # flag bands that are only partially covered or empty
            # TODO: Can interval ever be negative?
            interval = passband_integrated_width(wave, flux, passband=m, log=log)
            interval_frac = interval / numpy.diff(m)[0]
            self.main_incomplete[i] = interval_frac < 1.0
            self.main_empty[i] = numpy.invert(interval_frac > 0.0)

            # Common to both calculations of the BF indices
            # TODO: Need to catch divbyzero here, as well!
            bf = self.main_flux[i] * util.inverse(self.main_continuum[i])
            bf_err = 0.0 if err is None else \
                        numpy.sqrt(numpy.square(self.main_flux_err[i] 
                                                * util.inverse(self.main_flux[i]))
                                      + numpy.square(self.main_continuum_err[i]
                                                * util.inverse(self.main_continuum[i])))

            if self.units[i] == 'mag':
                # Calculation of the index in mag units requires
                # division by the band interval. If the full interval
                # is masked, this leads to a division by 0. This
                # catches that issue.
                # TODO: This is more general than a division by zero.
                # It should catch any computation that would result in
                # a NaN or Inf. ``index`` and ``bf`` need to be larger
                # than 0 for the log computation to work.
                self.divbyzero[i] = not (numpy.absolute(interval) > 0 and self.index[i] > 0 
                                            and bf > 0)
                if not self.divbyzero[i]:

                    # Convert the index to magnitudes using Worthey et
                    # al. 1994, eqn. 3: The passband interval cancels
                    # out of the error propagation. The error
                    # calculation is done first so as to not replace
                    # the linear calculation of the index. NOTE: If err
                    # is None, self.index_err=0, so there's no need to
                    # check index_err here.
                    self.index_err[i] = numpy.absolute(2.5*self.index_err[i]/self.index[i]/lg10)
                    self.index[i] = -2.5 * numpy.log10(self.index[i]/interval)

                    # Similarly for the BF definition:
                    self.index_bf[i] = -2.5 * numpy.log10(bf / interval)
                    self.index_bf_err[i] = 2.5 * bf_err / lg10
            else:
                # Complete the BF definition:
                self.index_bf[i] = interval - bf
                self.index_bf_err[i] = numpy.absolute(bf*bf_err)


class BandheadIndices:
    r"""
    Measure a set of bandhead, or "color", indices in a single
    spectrum.

    .. include:: ../include/bhdindices.rst

    .. note::

        The calculations performed by this class are agnostic as to
        the redshift of the spectrum provided. It is expected that
        you will have either de-redshifted the spectrum or redshifted
        the passband definitions appropriately to measure the index
        in a redshifted spectrum. Because the indices are flux
        ratios, indices calculated for a redshifted spectrum (with
        appropriately defined passbands) are identical to those
        calculated shifted to the rest frame.

        Also, the integrand is always the provided flux vector,
        meaning that its units should be appropriate to *all* of the
        indices to be calculated. For example, the D4000 index uses
        continua measured by the integral :math:`\int F_\nu {\rm
        d}\lambda`, meaning the spectrum provided to this class
        should be in per-frequency units.

    Args:
        wave (`numpy.ndarray`_):
            Wavelength vector.
        flux (array-like):
            Flux vector. Units must be appropriate to the definition
            of the provided index passbands. Can be a
            `numpy.ma.MaskedArray`_; masked pixels will be ignored in
            the measurement.
        bluebands (`numpy.ndarray`_):
            Array with the definition of the blue sidebands. Shape
            must be :math:`(N_{\rm index},2)`, where the second axis
            provides the starting and ending wavelength of each band.
        redbands (`numpy.ndarray`_):
            Array with the definition of the red sidebands. Shape
            must be :math:`(N_{\rm index},2)`, where the second axis
            provides the starting and ending wavelength of each band.
        err (`numpy.ndarray`_, optional):
            Error in the flux, used for a nominal propagation of the
            index error.
        log (`numpy.ndarray`_, optional):
            Flag that the provided flux vector is sampled
            logarithmically in wavelength.
        order (`numpy.ndarray`_, optional):
            String array signifying which band is in the numerator
            and which in the denominator. The string must be
            ``'r_b'`` --- signifying the index is calculated as the
            red continuum divided by the blue continuum --- or
            ``'b_r'``. If None, ``'r_b'`` is assumed for all
            measurements.
    
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
        if len(self.order) != self.nindx:
            raise ValueError('Must provide the order for each index.')

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

        self.divbyzero = numpy.zeros(self.nindx, dtype=bool)

        # Calculate the index in the correct order
        blue_n = order == 'b_r'
        self.index = numpy.ma.zeros(self.nindx, dtype=float)
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

    .. todo::

        **Detail what should be provided in terms of the redshift.**

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
        artifact_path (:obj:`str`, optional):
            Path to the directory with artifact parameter files.  If
            None, defined by
            :attr:`mangadap.par.artifactdb.ArtifactDB.default_data_dir`.
        absorption_index_path (:obj:`str`, optional):
            Path to the directory with absorption-line index parameter
            files.  If None, defined by
            :attr:`mangadap.par.absorptionindexdb.AbsorptionIndexDB.default_data_dir`.
        bandhead_index_path (:obj:`str`, optional):
            Path to the directory with bandhead index parameter files.
            If None, defined by
            :attr:`mangadap.par.bandheadindexdb.BandheadIndexDB.default_data_dir`.
        dapsrc (str): (**Optional**) Root path to the DAP source
            directory.  If not provided, the default is defined by
            :func:`mangadap.config.defaults.dap_source_dir`.
        dapver (str): (**Optional**) The DAP version to use for the
            analysis, used to override the default defined by
            :func:`mangadap.config.defaults.dap_version`.
        analysis_path (str): (**Optional**) The top-level path for the
            DAP output files, used to override the default defined by
            :func:`mangadap.config.defaults.dap_analysis_path`.
        directory_path (str): The exact path to the directory with DAP
            output that is common to number DAP "methods".  See
            :attr:`directory_path`.
        output_file (str): (**Optional**) Exact name for the output
            file.  The default is to use
            :func:`mangadap.config.defaults.dap_file_name`.
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
                 emission_line_model=None, database_list=None, artifact_path=None,
                 absorption_index_path=None, bandhead_index_path=None, dapsrc=None, dapver=None,
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
                               artifact_path=artifact_path,
                               absorption_index_path=absorption_index_path,
                               bandhead_index_path=bandhead_index_path)

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
        self.bitmask = SpectralIndicesBitMask()

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
                     emission_line_model=emission_line_model, dapver=dapver,
                     analysis_path=analysis_path, directory_path=directory_path,
                     output_file=output_file, hardcopy=hardcopy, tpl_symlink_dir=tpl_symlink_dir,
                     clobber=clobber, loggers=loggers, quiet=quiet)

    def __getitem__(self, key):
        return self.hdu[key]

    def _define_databases(self, database_key, database_list=None, artifact_path=None,
                          absorption_index_path=None, bandhead_index_path=None):
        r"""

        Select the database of indices

        """
        # Grab the specific database
        self.database = util.select_proc_method(database_key, SpectralIndicesDef,
                                                method_list=database_list,
                                                available_func=available_spectral_index_databases)

        # Instantiate the artifact, absorption-index, and bandhead-index
        # databases
        self.artdb = None if self.database['artifacts'] is None else \
                ArtifactDB.from_key(self.database['artifacts'], directory_path=artifact_path)
        # TODO: Generalize the name of this object
        self.pixelmask = SpectralPixelMask(artdb=self.artdb)

        self.absdb = None if self.database['absindex'] is None else \
                AbsorptionIndexDB.from_key(self.database['absindex'],
                                           directory_path=absorption_index_path)

        self.bhddb = None if self.database['bandhead'] is None else \
                BandheadIndexDB.from_key(self.database['bandhead'],
                                         directory_path=bandhead_index_path)

    def _set_paths(self, directory_path, dapver, analysis_path, output_file):
        """
        Set the :attr:`directory_path` and :attr:`output_file`.  If not
        provided, the defaults are set using, respectively,
        :func:`mangadap.config.defaults.dap_common_path` and
        :func:`mangadap.config.defaults.dap_file_name`.

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
        continuum_templates = 'None' if self.stellar_continuum is None \
                            else self.stellar_continuum.method['fitpar']['template_library_key']
        eml_templates = 'None' if self.emission_line_model is None \
                            else self.emission_line_model.method['continuum_tpl_key']
        method = defaults.dap_method(self.binned_spectra.method['key'], continuum_templates,
                                     eml_templates)
        self.analysis_path = defaults.dap_analysis_path(drpver=self.binned_spectra.cube.drpver,
                                                        dapver=dapver) \
                                    if analysis_path is None else str(analysis_path)
        self.directory_path \
                = defaults.dap_method_path(method, plate=self.binned_spectra.cube.plate,
                                           ifudesign=self.binned_spectra.cube.ifudesign, ref=True,
                                           drpver=self.binned_spectra.cube.drpver, dapver=dapver,
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

        self.output_file = defaults.dap_file_name(self.binned_spectra.cube.plate,
                                                  self.binned_spectra.cube.ifudesign, ref_method) \
                                        if output_file is None else str(output_file)

    def _initialize_primary_header(self, hdr=None, measurements_binid=None):
        """
        Construct the primary header for the reference file.

        Args:
            hdr (`astropy.io.fits.Header`_, optional):
                Input base header for added keywords. If None, uses
                the :attr:`cube` header (if there is one) and then
                cleans the header using
                :func:`mangadap.util.fitsutil.DAPFitsUtil.clean_dap_primary_header`.
            measurements_binid (`numpy.ndarray`_, optional):
                Bin IDs for spectral index measurements. Only use is
                to check if this is None, and the boolean result is
                saved to the header to indicate if the spectral
                indices are disconnected from the stellar-continuum
                measurements.

        Returns:
            `astropy.io.fits.Header`_: Initialized header object.
        """
        # Copy the from the DRP and clean it
        if hdr is None:
            hdr = self.binned_spectra.cube.prihdr.copy()
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

        t = 0 if self.absdb is None else self.absdb.size
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
        r"""
        Set the redshift to use for each spectrum for the spectral index
        measurements.
        
        In terms of precedence, directly provided redshifts override
        those in any available StellarContinuumModel.

        If :attr:`stellar_continuum` and redshift are None, the default
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

            redshift (:obj:`float`, numpy.ndarray):
                Redshifts (:math:`z`) to use for each spectrum.  If
                None, the default
            measure_on_unbinned_spaxels (:obj:`bool`):
                Flag that method expects to measure moments on unbinned
                spaxels.
            good_snr (numpy.ndarray):
                Boolean array setting which spectra have sufficient S/N
                for the measurements.
            default_redshift (:obj:`float`, optional):
                Only used if there are stellar kinematics available.
                Provides the default redshift to use for spectra
                without stellar measurements; see ``redshift`` in
                :func:`mangadap.proc.stellarcontinuummodel.StellarContinuumModel.matched_kinematics`.
                If None (default), the median of the unmasked stellar
                velocities will be used.

        """

        # Construct the binid matrix if measuring on unbinned spaxels
        if measure_on_unbinned_spaxels:
            binid = numpy.full(self.binned_spectra.spatial_shape, -1, dtype=int)
            binid.ravel()[good_snr] = numpy.arange(self.nbins)
            missing = []
            nspec = self.binned_spectra.cube.nspec
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
            return self.binned_spectra.check_fgoodpix()
#            fgoodpix = self.binned_spectra.check_fgoodpix()
#            good_snr = self.binned_spectra.rdxqa['SPECTRUM'].data['SNR'] \
#                                > self.database['minimum_snr']
#            return fgoodpix & good_snr
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
                ``measure_on_unbinned_spaxels``.  Default is to
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
                and emission_line_model.method['deconstruct_bins'] == 'ignore':
            raise ValueError('Cannot use this emission-line model with unbinned spaxels.')

        # Get the main data arrays
        if measure_on_unbinned_spaxels:
            # TODO: Should probably make this a function within
            # SpatiallyBinnedSpectra, particularly because of the
            # dereddening
            flags = binned_spectra.cube.do_not_fit_flags()
            binid = numpy.arange(binned_spectra.cube.nspec).reshape(binned_spectra.spatial_shape)
            wave = binned_spectra['WAVE'].data
            flux = binned_spectra.cube.copy_to_masked_array(flag=flags)
            ivar = binned_spectra.cube.copy_to_masked_array(attr='ivar', flag=flags)
            flux, ivar = binned_spectra.galext.apply(flux, ivar=ivar, deredden=True)
            missing = None
            sres = binned_spectra.cube.copy_to_array(attr='sres')
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

        # Use :func:`mangadap.util.resolution.match_spectral_resolution`
        # to match the spectral resolution of the binned spectra to the
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

        # From resolution.py:  "Any pixel that had a resolution that was
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

    def _resolution_matched_templates(self, dapver=None, analysis_path=None, tpl_symlink_dir=None):
        """
        Get a version of the template library that has had its
        resolution matched to that of the spectral-index database.
        """
        velocity_offset=self.stellar_continuum.method['fitpar']['template_library'].velocity_offset
        return self.stellar_continuum.get_template_library(dapver=dapver,
                                                           analysis_path=analysis_path,
                                                           tpl_symlink_dir=tpl_symlink_dir,
                                                           velocity_offset=velocity_offset,
                                                           resolution_fwhm=self.database['fwhm'])

    @staticmethod
    def calculate_dispersion_corrections(absdb, bhddb, wave, flux, continuum, continuum_dcnvlv,
                                         redshift=None, redshift_dcnvlv=None, bitmask=None):
        r"""
        Calculate the spectral-index velocity-dispersion corrections.

        The indices are measured twice, once on the best-fitting
        model (``continuum``; produces :math:`{\mathcal I}_{\rm
        mod}`), and once on the "deconvolved" best-fitting model
        (``continuum_dcnvlv``); i.e., the best-fitting model but with
        0 velocity dispersion, :math:`{\mathcal I}_{\sigma=0}`. The
        corrections, :math:`\delta {\rm mathcal I}` are then
        constructed as:

        .. math::

            \delta{\mathcal I} = \left\{
                \begin{array}{ll}
                    {\mathcal I}_{\sigma=0}\ /\ {\mathcal I}_{\rm mod}, & \mbox{for unitless indices} \\[3pt]
                    {\mathcal I}_{\sigma=0}\ /\ {\mathcal I}_{\rm mod}, & \mbox{for angstrom units} \\[3pt]
                    {\mathcal I}_{\sigma=0} - {\mathcal I}_{\rm mod}, & \mbox{for magnitude units}
                \end{array}\right.,

        where "bandhead" (color) indices are treated as unitless. The
        construction of the two continuum spectra can be such that
        they have different redshifts (see ``redshift`` and
        ``redshift_dcnvlv``).  The corrected index values are then:

        .. math::

            {\mathcal I}^c = \left\{
                \begin{array}{ll}
                    {\mathcal I}\ \delta{\mathcal I}, & \mbox{for unitless indices} \\[3pt]
                    {\mathcal I}\ \delta{\mathcal I}, & \mbox{for angstrom units} \\[3pt]
                    {\mathcal I} + \delta{\mathcal I}, & \mbox{for magnitude units} \\[3pt]
                \end{array}\right..

        For absorption-line indices, the index corrections are
        calculated for both the Worthey/Trager definition and the
        Burstein/Faber definition.

        Similar corrections are calculated for the continuum in the
        main passband, as well as the two sidebands:

        .. math::

            \delta C = C_{\sigma=0} / C_{\rm mod},

        where :math:`C` is used here as a generic continuum
        measurement (cf. :class:`AbsorptionLineIndices` documentation
        nomenclature).

        Args:
            absdb (:class:`~mangadap.par.aborptionlinedb.AbsorptionIndexDB`):
                Database with the absorption-line index definitions.
            bhddb (:class:`~mangadap.par.bandheadindexdb.BandheadIndexDB`):
                Database with the bandhead index definitions.
            wave (`numpy.ndarray`_):
                1D vector with the wavelength of each pixel. *Assumed
                to be logarithmically binned in radius.*
            flux (`numpy.ndarray`_):
                2D array with the flux, ordered as :math:`N_{\rm
                spec}\times N_{\rm wave}`. Can be a
                numpy.ma.MaskedArray; masked pixels will be ignored
                in the measurement. This is only used to get the
                masked pixels to ignore during the index computation
                in the continuum model.s
            continuum (`numpy.ndarray`_):
                Best-fitting continuum-only models to the spectra.
            continuum_dcnvlv (`numpy.ndarray`_):
                Best-fitting continuum-only models to the spectra,
                but with a velocity dispersion of 0 km/s.
            redshift (array-like, optional):
                Redshift to use for each spectrum when determining
                the index. Must have the correct length compared to
                the flux array. Default is to assume the spectra are
                at rest wavelength.
            bitmask (:class:`mangadap.util.bitmask.BitMask`, optional):
                If an index is flagged for some reason (see
                :func:`set_masks`), this object is used to set the
                mask value; this should typically be
                :class:`SpectralIndicesBitMask` object. If not
                provided, the masked values are all set to True
                (False otherwise).

        Returns:
            :obj:`tuple`: Returns 16 (!) `numpy.ndarray`_ objects:
            The model blue pseudo-continua and corrections, the model
            red pseudo-continua and corrections, the model main
            passband continua and corrections, the model weights and
            corrections, the model Worthey/Trager indices and
            corrections, the model Burstein/Faber indices and
            corrections, and boolean arrays selecting the good
            unitless indices, angstrom-unit indices, and
            magnitude-unit indices, and absorption-line indices.
        """
        # Make sure continuum includes the flux masks
        nspec = flux.shape[0]
        _flux = numpy.ma.MaskedArray(flux)
        _continuum = numpy.ma.MaskedArray(continuum)
        _continuum_dcnvlv = numpy.ma.MaskedArray(continuum_dcnvlv)
        _continuum[numpy.ma.getmaskarray(_flux)] = numpy.ma.masked
        _continuum_dcnvlv[numpy.ma.getmaskarray(_flux)] = numpy.ma.masked

        # Measure the indices for both sets of spectra
        indx = SpectralIndices.measure_indices(absdb, bhddb, wave, _continuum, redshift=redshift,
                                               bitmask=bitmask)
        dcnvlv_indx = SpectralIndices.measure_indices(absdb, bhddb, wave, _continuum_dcnvlv,
                                                      redshift=redshift_dcnvlv, bitmask=bitmask)

        # Do not compute corrections if any of the bands were empty or
        # would have had to divide by zero
        bad_flags = ['MAIN_EMPTY', 'BLUE_EMPTY', 'RED_EMPTY', 'DIVBYZERO']
        bad_indx = bitmask.flagged(indx['MASK'], flag=bad_flags)
        bad_indx |= bitmask.flagged(dcnvlv_indx['MASK'], flag=bad_flags)

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

        #---------------------------------------------------------------
        # Sideband model values and corrections (same computation
        # regardless of the units).

        good = good_les | good_ang | good_mag
        bcont_mod = numpy.zeros(indx['BCONT'].shape, dtype=float)
        bcont_mod[good] = indx['BCONT'][good]
        bcont_corr = numpy.zeros(indx['BCONT'].shape, dtype=float)
        bcont_corr[good] = dcnvlv_indx['BCONT'][good] / indx['BCONT'][good]

        rcont_mod = numpy.zeros(indx['RCONT'].shape, dtype=float)
        rcont_mod[good] = indx['RCONT'][good]
        rcont_corr = numpy.zeros(indx['RCONT'].shape, dtype=float)
        rcont_corr[good] = dcnvlv_indx['RCONT'][good] / indx['RCONT'][good]

        # The weights are just continuum values drawn from different
        # sources, so the model and correction compuations are done in
        # the same way.
        awgt_mod = numpy.zeros(indx['AWGT'].shape, dtype=float)
        awgt_mod[good] = indx['AWGT'][good]
        awgt_corr = numpy.zeros(indx['AWGT'].shape, dtype=float)
        awgt_corr[good] = dcnvlv_indx['AWGT'][good] / indx['AWGT'][good]

        #---------------------------------------------------------------
        # Continuum in the main passband (same computation regardless
        # of the units).

        # Only valid for the absorption-line indices
        nabs, nbhd = SpectralIndices.count_indices(absdb, bhddb)
        is_abs = numpy.ones(nabs+nbhd, dtype=bool)
        is_abs[nabs:] = False
        _good = good & is_abs[None,:]

        # Save the *good* indices for the best-fitting models
        mcont_mod = numpy.zeros(indx['MCONT'].shape, dtype=float)
        mcont_mod[_good] = indx['MCONT'][_good]

        # Determine the dispersion corrections
        mcont_corr = numpy.zeros(indx['MCONT'].shape, dtype=float)
        mcont_corr[_good] = dcnvlv_indx['MCONT'][_good] / indx['MCONT'][_good]

        #---------------------------------------------------------------
        # Worthey/Trager indices and bandhead indices
        # Save the *good* indices for the best-fitting models
        indx_mod = numpy.zeros(indx['INDX'].shape, dtype=float)
        indx_mod[good] = indx['INDX'][good]

        # Determine the dispersion corrections
        indx_corr = numpy.zeros(indx['INDX'].shape, dtype=float)
        indx_corr[good_les] = dcnvlv_indx['INDX'][good_les] / indx['INDX'][good_les]
        indx_corr[good_ang] = dcnvlv_indx['INDX'][good_ang] / indx['INDX'][good_ang]
        indx_corr[good_mag] = dcnvlv_indx['INDX'][good_mag] - indx['INDX'][good_mag]

        #---------------------------------------------------------------
        # Burstein/Faber indices. This will include the bandhead/color
        # indices; however, data in INDX_BF is currently just a copy of
        # what's in INDX.

        # Save the *good* indices for the best-fitting models
        indx_bf_mod = numpy.zeros(indx['INDX_BF'].shape, dtype=float)
        indx_bf_mod[good_les | good_ang | good_mag] \
                = indx['INDX_BF'][good_les | good_ang | good_mag]

        # Determine the dispersion corrections
        indx_bf_corr = numpy.zeros(indx['INDX_BF'].shape, dtype=float)
        indx_bf_corr[good_les] = dcnvlv_indx['INDX_BF'][good_les] / indx['INDX_BF'][good_les]
        indx_bf_corr[good_ang] = dcnvlv_indx['INDX_BF'][good_ang] / indx['INDX_BF'][good_ang]
        indx_bf_corr[good_mag] = dcnvlv_indx['INDX_BF'][good_mag] - indx['INDX_BF'][good_mag]

        # Return the results
        return bcont_mod, bcont_corr, rcont_mod, rcont_corr, mcont_mod, mcont_corr, \
               awgt_mod, awgt_corr, indx_mod, indx_corr, indx_bf_mod, indx_bf_corr, \
               good_les, good_ang, good_mag, is_abs

    @staticmethod
    def apply_dispersion_corrections(indx, indxcorr, err=None, unit=None):
        """
        Apply a set of dispersion corrections.  Errors in the dispersion
        corrections are assumed to be negligible.

        Args:
            indx (array-like):
                Indices to correct.
            indxcorr (array-like):
                Index corrections.
            err (array-like, optional):
                Error in the indices.
            unit (:obj:`str`, optional):
                Unit of the index; must be either magnitudes
                (``'mag'``) or angstroms (``'ang'``) or None. Default
                is None. Unitless corrections and angstrom
                corrections are treated identically.
    
        Returns:
            :obj:`tuple`: Returns two objects with the corrected
            indices and errors. If no errors are returned, the second
            returned object is None.

        Raises:
            ValueError:
                Raised if the unit is not ang or mag.
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
        return 0 if absdb is None else absdb.size, 0 if bhddb is None else bhddb.size

#    @staticmethod
#    def output_dtype(nindx, bitmask=None):
#        r"""
#        Construct the record array data type for the output fits
#        extension.
#        """
#        return [ ('BINID', numpy.int),
#                 ('BINID_INDEX',numpy.int),
#                 ('REDSHIFT', numpy.float),
#                 ('MASK', numpy.bool if bitmask is None else bitmask.minimum_dtype(), (nindx,)),
#                 ('BCEN', numpy.float, (nindx,)), 
#                 ('BCONT', numpy.float, (nindx,)), 
#                 ('BCONT_ERR', numpy.float, (nindx,)),
#                 ('BCONT_MOD', numpy.float, (nindx,)),
#                 ('BCONT_CORR', numpy.float, (nindx,)),
#                 ('RCEN', numpy.float, (nindx,)), 
#                 ('RCONT', numpy.float, (nindx,)), 
#                 ('RCONT_ERR', numpy.float, (nindx,)),
#                 ('RCONT_MOD', numpy.float, (nindx,)),
#                 ('RCONT_CORR', numpy.float, (nindx,)),
#                 ('MCONT', numpy.float, (nindx,)), 
#                 ('MCONT_ERR', numpy.float, (nindx,)),
#                 ('MCONT_MOD', numpy.float, (nindx,)), 
#                 ('MCONT_CORR', numpy.float, (nindx,)),
#                 ('AWGT', numpy.float, (nindx,)), 
#                 ('AWGT_ERR', numpy.float, (nindx,)),
#                 ('AWGT_MOD', numpy.float, (nindx,)), 
#                 ('AWGT_CORR', numpy.float, (nindx,)),
#                 ('INDX', numpy.float, (nindx,)), 
#                 ('INDX_ERR', numpy.float, (nindx,)),
#                 ('INDX_MOD', numpy.float, (nindx,)), 
#                 ('INDX_CORR', numpy.float, (nindx,)),
#                 ('INDX_BF', numpy.float, (nindx,)), 
#                 ('INDX_BF_ERR', numpy.float, (nindx,)),
#                 ('INDX_BF_MOD', numpy.float, (nindx,)), 
#                 ('INDX_BF_CORR', numpy.float, (nindx,))
#               ]

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
        _redshift = numpy.zeros(nspec, dtype=float) if redshift is None else redshift
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
        """
        Save the index measuremements to the measurement database.
        """
        is_abs = isinstance(results, AbsorptionLineIndices)
        is_bhd = isinstance(results, BandheadIndices)
        if not is_abs and not is_bhd:
            raise TypeError('Input must be of type AbsorptionLineIndices or BandheadIndices')

        # Save the data, general to both absorption-line and bandhead
        # indices
        measurements['BCEN'][good] = results.blue_center
        measurements['BCONT'][good] = results.blue_continuum
        if err:
            measurements['BCONT_ERR'][good] = results.blue_continuum_err
            
        measurements['RCEN'][good] = results.red_center
        measurements['RCONT'][good] = results.red_continuum
        if err:
            measurements['RCONT_ERR'][good] = results.red_continuum_err

        measurements['INDX'][good] = results.index
        if err:
            measurements['INDX_ERR'][good] = results.index_err

        # Type specific
        if is_abs:
            # Absorption-line indices also save the 2nd index
            # definition and the continuum in the main band
            measurements['INDX_BF'][good] = results.index_bf
            measurements['MCONT'][good] = results.main_continuum
            if err:
                measurements['INDX_BF_ERR'][good] = results.index_bf_err
                measurements['MCONT_ERR'][good] = results.main_continuum_err
        else:
            # Bandhead indices have 1 definition, but they're saved in
            # both places.
            measurements['INDX_BF'][good] = results.index
            if err:
                measurements['INDX_BF_ERR'][good] = results.index_err

        # Reshape the flags
        blue_incomplete = numpy.zeros(len(good), dtype=bool)
        blue_incomplete[good] = results.blue_incomplete

        blue_empty = numpy.zeros(len(good), dtype=bool)
        blue_empty[good] = results.blue_empty
        
        red_incomplete = numpy.zeros(len(good), dtype=bool)
        red_incomplete[good] = results.red_incomplete
        
        red_empty = numpy.zeros(len(good), dtype=bool)
        red_empty[good] = results.red_empty
        
        main_incomplete = None
        if results.main_incomplete is not None:
            main_incomplete = numpy.zeros(len(good), dtype=bool)
            main_incomplete[good] = results.main_incomplete

        main_empty = None
        if results.main_empty is not None:
            main_empty = numpy.zeros(len(good), dtype=bool)
            main_empty[good] = results.main_empty
            
        divbyzero = numpy.zeros(len(good), dtype=bool)
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

        The returned object is :class:`SpectralIndicesDataTable`; see
        the documentation of this object for the list of returned
        quantities.

        See the documentation of :class:`AbsorptionLineIndices` and
        :class:`BandheadIndices` for the definitions of the index
        measurements. For the purpose of their storage here, ``INDX``
        refers to :math:`{\mathcal I}_{\rm WT}` for the
        absorption-line indices. There is only one definition for the
        bandhead/color indices, and they are included in *both*
        ``INDX`` and ``INDX_BF`` (see :func:`save_results`).

        The weights to use when combining indices are simply copies
        of the data in other columns, collected for convenience (and
        to ease the construction of the main DAP output file). For
        all absorption-line indices, ``AWGT`` is identically
        ``MCONT``; for bandhead indices, ``AWGT`` is ``BCONT`` and
        ``RCONT`` for indices with order ``'r_b'`` and ``'b_r'``
        order, respectively.

        This function does not add the ``BINID``, ``BINID_INDEX``,
        ``*_MOD`` or ``*_CORR`` values. All except the first 3
        columns vectors with a length of :math:`N_{\rm index}` ---
        the total number of indices calculated (see
        :func:`count_indices`).

        Args:
            absdb (:class:`mangadap.par.aborptionlinedb.AbsorptionIndexDB`):
                Database with the absorption-line index definitions.
                Can be None.
            bhddb (:class:`mangadap.par.bandheadindexdb.BandheadIndexDB`):
                Database with the bandhead index definitions.
            wave (array-like):
                1D vector with the wavelength of each pixel. *Assumed
                to be logarithmically binned in radius.*
            flux (array-like):
                2D array with the flux, ordered as :math:`N_{\rm
                spec}\times N_{\rm wave}`. Can be a
                numpy.ma.MaskedArray; masked pixels will be ignored
                in the measurement.
            ivar (array-like, optional):
                Inverse variance in the flux. Must match flux array
                shape. Used to calculate propagated errors in the
                index. Default is that errors are ignored.
            mask (array-like, optional):
                Boolean array flagging to ignore (mask=True) or
                include (mask=False) each flux measurement in the
                index calculation.
            redshift (array-like, optional):
                Redshift to use for each spectrum when determining
                the index. Must have the correct length compared to
                the flux array. Default is to assume the spectra are
                at rest wavelength.
            bitmask (:class:`mangadap.util.bitmask.BitMask`, optional):
                If an index is flagged for some reason (see
                :func:`set_masks`), this object is used to set the
                mask value; this should typically be
                :class:`SpectralIndicesBitMask` object. If not
                provided, the masked values are all set to True
                (False otherwise).

        Returns:
            :class:`SpectralIndicesDataTable`: The object with the
            index-measurement data.
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
        measurements = SpectralIndicesDataTable(nindx=nindx, bitmask=bitmask, shape=nspec)
#        measurements = init_record_array(nspec, SpectralIndices.output_dtype(nindx,bitmask=bitmask))
        _flux, noise, measurements['REDSHIFT'] \
                    = SpectralIndices.check_and_prep_input(wave, flux, ivar=ivar, mask=mask,
                                                           redshift=redshift, bitmask=bitmask)

        # Create the f_nu spectra. NOTE: The conversion is
        # multiplicative, meaning the calculation of the errors can use
        # exactly the same function
        flux_fnu = util.flux_to_fnu(numpy.array([wave]*nspec), _flux)
        noise_fnu = None if noise is None else util.flux_to_fnu(numpy.array([wave]*nspec), noise)

        # Get the list of good indices of each type
        abs_fnu = numpy.zeros(nabs, dtype=bool) if absdb is None \
                        else numpy.invert(absdb.dummy) & (absdb['integrand'] == 'fnu')
        good_abs_fnu = numpy.zeros(nindx, dtype=bool)
        if nabs > 0:
            good_abs_fnu[:nabs][abs_fnu] = True
        abs_flambda = numpy.zeros(nabs, dtype=bool) if absdb is None \
                        else numpy.invert(absdb.dummy) & (absdb['integrand'] == 'flambda')
        good_abs_flambda = numpy.zeros(nindx, dtype=bool)
        if nabs > 0:
            good_abs_flambda[:nabs][abs_flambda] = True

        bhd_fnu = numpy.zeros(nbhd, dtype=bool) if bhddb is None \
                        else numpy.invert(bhddb.dummy) & (bhddb['integrand'] == 'fnu')
        good_bhd_fnu = numpy.zeros(nindx, dtype=bool)
        if nbhd > 0:
            good_bhd_fnu[nabs:][bhd_fnu] = True
        bhd_flambda = numpy.zeros(nbhd, dtype=bool) if bhddb is None \
                        else numpy.invert(bhddb.dummy) & (bhddb['integrand'] == 'flambda')
        good_bhd_flambda = numpy.zeros(nindx, dtype=bool)
        if nbhd > 0:
            good_bhd_flambda[nabs:][bhd_flambda] = True

        # Mask any dummy indices
        dummy = numpy.zeros(nindx, dtype=bool)
        if absdb is not None:
            dummy[:nabs] = absdb.dummy
        if bhddb is not None:
            dummy[nabs:] = bhddb.dummy
        if numpy.any(dummy):
            measurements['MASK'][:,dummy] = True if bitmask is None else \
                    bitmask.turn_on(measurements['MASK'][:,dummy], 'UNDEFINED_BANDS')

        # No valid indices
        if numpy.all(dummy):
            return measurements

#        warnings.simplefilter("error", RuntimeWarning)

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
                    # NOTE: This basically should never happen so I throw a warning.
                    warnings.warn('It is odd to construct an absorption-line index by '
                                  'integrating F_nu over lambda; just sayin...')
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

                    # Observed wavelength converted to rest wavelength
                    measurements['BCEN'][i,good_abs_fnu] /= (1+measurements['REDSHIFT'][i])
                    measurements['RCEN'][i,good_abs_fnu] /= (1+measurements['REDSHIFT'][i])
                    # Observed F_nu continuum should be divided by
                    # (1+z) to convert to the rest-frame values
                    measurements['BCONT'][i,good_abs_fnu] /= (1+measurements['REDSHIFT'][i])
                    measurements['BCONT_ERR'][i,good_abs_fnu] /= (1+measurements['REDSHIFT'][i])
                    measurements['RCONT'][i,good_abs_fnu] /= (1+measurements['REDSHIFT'][i])
                    measurements['RCONT_ERR'][i,good_abs_fnu] /= (1+measurements['REDSHIFT'][i])
                    measurements['MCONT'][i,good_abs_fnu] /= (1+measurements['REDSHIFT'][i])
                    measurements['MCONT_ERR'][i,good_abs_fnu] /= (1+measurements['REDSHIFT'][i])

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

                    # Observed wavelength converted to rest wavelength
                    measurements['BCEN'][i,good_abs_flambda] /= (1+measurements['REDSHIFT'][i])
                    measurements['RCEN'][i,good_abs_flambda] /= (1+measurements['REDSHIFT'][i])
                    # Observed F_lambda continuum should be multiplied
                    # by (1+z) to convert to the rest-frame values
                    measurements['BCONT'][i,good_abs_flambda] *= (1+measurements['REDSHIFT'][i])
                    measurements['BCONT_ERR'][i,good_abs_flambda] *= (1+measurements['REDSHIFT'][i])
                    measurements['RCONT'][i,good_abs_flambda] *= (1+measurements['REDSHIFT'][i])
                    measurements['RCONT_ERR'][i,good_abs_flambda] *= (1+measurements['REDSHIFT'][i])
                    measurements['MCONT'][i,good_abs_flambda] *= (1+measurements['REDSHIFT'][i])
                    measurements['MCONT_ERR'][i,good_abs_flambda] *= (1+measurements['REDSHIFT'][i])

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

                    # Observed wavelength converted to rest wavelength
                    measurements['BCEN'][i,good_bhd_fnu] /= (1+measurements['REDSHIFT'][i])
                    measurements['RCEN'][i,good_bhd_fnu] /= (1+measurements['REDSHIFT'][i])
                    # Observed F_nu continuum should be divided by
                    # (1+z) to convert to the rest-frame values
                    measurements['BCONT'][i,good_bhd_fnu] /= (1+measurements['REDSHIFT'][i])
                    measurements['BCONT_ERR'][i,good_bhd_fnu] /= (1+measurements['REDSHIFT'][i])
                    measurements['RCONT'][i,good_bhd_fnu] /= (1+measurements['REDSHIFT'][i])
                    measurements['RCONT_ERR'][i,good_bhd_fnu] /= (1+measurements['REDSHIFT'][i])

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

                    # Observed wavelength converted to rest wavelength
                    measurements['BCEN'][i,good_bhd_flambda] /= (1+measurements['REDSHIFT'][i])
                    measurements['RCEN'][i,good_bhd_flambda] /= (1+measurements['REDSHIFT'][i])
                    # Observed F_lambda continuum should be multiplied
                    # by (1+z) to convert to the rest-frame values
                    measurements['BCONT'][i,good_bhd_flambda] *= (1+measurements['REDSHIFT'][i])
                    measurements['BCONT_ERR'][i,good_bhd_flambda] *= (1+measurements['REDSHIFT'][i])
                    measurements['RCONT'][i,good_bhd_flambda] *= (1+measurements['REDSHIFT'][i])
                    measurements['RCONT_ERR'][i,good_bhd_flambda] *= (1+measurements['REDSHIFT'][i])

        print('Measuring spectral indices in spectrum: {0}/{0}'.format(nspec))

        # Correct the indices (and their errors) with angstrom units to
        # rest-frame
        _, angu, _ = SpectralIndices.unit_selection(absdb, bhddb)
        measurements['INDX'][:,angu] /= (1+measurements['REDSHIFT'][:,None])
        measurements['INDX_ERR'][:,angu] /= (1+measurements['REDSHIFT'][:,None])
        measurements['INDX_BF'][:,angu] /= (1+measurements['REDSHIFT'][:,None])
        measurements['INDX_BF_ERR'][:,angu] /= (1+measurements['REDSHIFT'][:,None])

        # Collect the weights. The weight is the same for both
        # definitions of absorption-line indices, but should really
        # only be applied to the 2nd definition (INDX_BF).
        is_abs = numpy.ones(nindx, dtype=bool)
        is_abs[nabs:] = False
        if nabs > 0:
            measurements['AWGT'][:,is_abs] = measurements['MCONT'][:,is_abs]
            measurements['AWGT_ERR'][:,is_abs] = measurements['MCONT_ERR'][:,is_abs]
        if nbhd > 0:
            use = numpy.logical_not(is_abs)
            use[nabs:] &= (bhddb['order'] == 'r_b')
            measurements['AWGT'][:,use] = measurements['BCONT'][:,use]
            measurements['AWGT_ERR'][:,use] = measurements['BCONT_ERR'][:,use]
            use = numpy.logical_not(is_abs)
            use[nabs:] &= (bhddb['order'] == 'b_r')
            measurements['AWGT'][:,use] = measurements['RCONT'][:,use]
            measurements['AWGT_ERR'][:,use] = measurements['RCONT_ERR'][:,use]

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
                emission_line_model=None, dapver=None, analysis_path=None, directory_path=None,
                output_file=None, hardcopy=True, tpl_symlink_dir=None, clobber=False, loggers=None,
                quiet=False):
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
            dapver (str): (**Optional**) The DAP version to use for the
                analysis, used to override the default defined by
                :func:`mangadap.config.defaults.dap_version`.
            analysis_path (str): (**Optional**) The top-level path for
                the DAP output files, used to override the default
                defined by
                :func:`mangadap.config.defaults.dap_analysis_path`.
            directory_path (str): The exact path to the directory with
                DAP output that is common to number DAP "methods".  See
                :attr:`directory_path`.
            output_file (str): (**Optional**) Exact name for the output
                file.  The default is to use
                :func:`mangadap.config.defaults.dap_file_name`.
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
                and self.emission_line_model.method['deconstruct_bins'] != 'ignore'

        self.spatial_shape =self.binned_spectra.spatial_shape
        self.nspec = self.binned_spectra.cube.nspec if measure_on_unbinned_spaxels \
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

        #---------------------------------------------------------------
        # Determine the velocity dispersion corrections, if requested.
        if self.compute_corrections:

            # Get the template spectra to use
            replacement_templates = None if self.database['fwhm'] < 0 \
                    else self._resolution_matched_templates(dapver=dapver,
                                                            analysis_path=self.analysis_path,
                                                            tpl_symlink_dir=tpl_symlink_dir)
            # Have to use the corrected velocity dispersion if templates
            # have been broadened to a new resolution; otherwise, the
            # corrections are to a zero dispersion model **at the
            # resolution of the templates**
            corrected_dispersion = self.database['fwhm'] > 0

            # Set the bin IDs to match the stellar continuum to:
            binid = numpy.arange(binned_spectra.cube.nspec).reshape(binned_spectra.spatial_shape) \
                        if measure_on_unbinned_spaxels else binned_spectra['BINID'].data

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

            measurements['BCONT_MOD'], measurements['BCONT_CORR'], measurements['RCONT_MOD'], \
                measurements['RCONT_CORR'], measurements['MCONT_MOD'], measurements['MCONT_CORR'], \
                measurements['AWGT_MOD'], measurements['AWGT_CORR'], \
                measurements['INDX_MOD'], measurements['INDX_CORR'], \
                measurements['INDX_BF_MOD'], measurements['INDX_BF_CORR'], \
                good_les, good_ang, good_mag, is_abs \
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
            # TODO: Get rid of this? I.e., correct_indices is *always*
            # False as currently coded...
            if self.correct_indices:
                #-------------------------------------------------------
                # Continuum in the two sidebands and the weights
                good = good_les | good_ang | good_mag
                # Measured continuum
                measurements['BCONT'][good], measurements['BCONT_ERR'][good], \
                        = SpectralIndices.apply_dispersion_corrections(
                                                    measurements['BCONT'][good],
                                                    measurements['BCONT_CORR'][good],
                                                    err=measurements['BCONT_ERR'][good])
                measurements['RCONT'][good], measurements['RCONT_ERR'][good], \
                        = SpectralIndices.apply_dispersion_corrections(
                                                    measurements['RCONT'][good],
                                                    measurements['RCONT_CORR'][good],
                                                    err=measurements['RCONT_ERR'][good])
                measurements['AWGT'][good], measurements['AWGT_ERR'][good], \
                        = SpectralIndices.apply_dispersion_corrections(
                                                    measurements['AWGT'][good],
                                                    measurements['AWGT_CORR'][good],
                                                    err=measurements['AWGT_ERR'][good])
                # Model continuum
                measurements['BCONT_MOD'][good], _ \
                        = SpectralIndices.apply_dispersion_corrections(
                                                    measurements['BCONT_MOD'][good],
                                                    measurements['BCONT_CORR'][good])
                measurements['RCONT_MOD'][good], _ \
                        = SpectralIndices.apply_dispersion_corrections(
                                                    measurements['RCONT_MOD'][good],
                                                    measurements['RCONT_CORR'][good])
                measurements['AWGT_MOD'][good], _ \
                        = SpectralIndices.apply_dispersion_corrections(
                                                    measurements['AWGT_MOD'][good],
                                                    measurements['AWGT_CORR'][good])

                #-------------------------------------------------------
                # Worthey/Trager indices
                # Model indices
                measurements['INDX'][good_les], measurements['INDX_ERR'][good_les] \
                        = SpectralIndices.apply_dispersion_corrections(
                                                    measurements['INDX'][good_les],
                                                    measurements['INDX_CORR'][good_les],
                                                    err=measurements['INDX_ERR'][good_les])
                measurements['INDX'][good_ang], measurements['INDX_ERR'][good_ang] \
                        = SpectralIndices.apply_dispersion_corrections(
                                                    measurements['INDX'][good_ang],
                                                    measurements['INDX_CORR'][good_ang],
                                                    err=measurements['INDX_ERR'][good_ang],
                                                    unit='ang')
                measurements['INDX'][good_mag], measurements['INDX_ERR'][good_mag] \
                        = SpectralIndices.apply_dispersion_corrections(
                                                    measurements['INDX'][good_mag],
                                                    measurements['INDX_CORR'][good_mag],
                                                    err=measurements['INDX_ERR'][good_mag],
                                                    unit='mag')
                # Model indices
                measurements['INDX_MOD'][good_les], _ \
                        = SpectralIndices.apply_dispersion_corrections(
                                                    measurements['INDX_MOD'][good_les],
                                                    measurements['INDX_CORR'][good_les])
                measurements['INDX_MOD'][good_ang], _ \
                        = SpectralIndices.apply_dispersion_corrections(
                                                    measurements['INDX_MOD'][good_ang],
                                                    measurements['INDX_CORR'][good_ang],
                                                    unit='ang')
                measurements['INDX_MOD'][good_mag], _ \
                        = SpectralIndices.apply_dispersion_corrections(
                                                    measurements['INDX_MOD'][good_mag],
                                                    measurements['INDX_CORR'][good_mag],
                                                    unit='mag')

                #-------------------------------------------------------
                # Burstein/Faber indices
                # Model indices
                measurements['INDX_BF'][good_les], measurements['INDX_BF_ERR'][good_les] \
                        = SpectralIndices.apply_dispersion_corrections(
                                                    measurements['INDX_BF'][good_les],
                                                    measurements['INDX_BF_CORR'][good_les],
                                                    err=measurements['INDX_BF_ERR'][good_les])
                measurements['INDX_BF'][good_ang], measurements['INDX_BF_ERR'][good_ang] \
                        = SpectralIndices.apply_dispersion_corrections(
                                                    measurements['INDX_BF'][good_ang],
                                                    measurements['INDX_BF_CORR'][good_ang],
                                                    err=measurements['INDX_BF_ERR'][good_ang],
                                                    unit='ang')
                measurements['INDX_BF'][good_mag], measurements['INDX_BF_ERR'][good_mag] \
                        = SpectralIndices.apply_dispersion_corrections(
                                                    measurements['INDX_BF'][good_mag],
                                                    measurements['INDX_BF_CORR'][good_mag],
                                                    err=measurements['INDX_BF_ERR'][good_mag],
                                                    unit='mag')
                # Model indices
                measurements['INDX_BF_MOD'][good_les], _ \
                        = SpectralIndices.apply_dispersion_corrections(
                                                    measurements['INDX_BF_MOD'][good_les],
                                                    measurements['INDX_BF_CORR'][good_les])
                measurements['INDX_BF_MOD'][good_ang], _ \
                        = SpectralIndices.apply_dispersion_corrections(
                                                    measurements['INDX_BF_MOD'][good_ang],
                                                    measurements['INDX_BF_CORR'][good_ang],
                                                    unit='ang')
                measurements['INDX_BF_MOD'][good_mag], _ \
                        = SpectralIndices.apply_dispersion_corrections(
                                                    measurements['INDX_BF_MOD'][good_mag],
                                                    measurements['INDX_BF_CORR'][good_mag],
                                                    unit='mag')

                #-------------------------------------------------------
                # Continuum in the main passband; specific to the
                # absorption-line indices
                _good_les = good_les & is_abs[None,:]
                _good_ang = good_ang & is_abs[None,:]
                _good_mag = good_mag & is_abs[None,:]
                good = _good_les | _good_ang | _good_mag
                # Measured continuum
                measurements['MCONT'][good], measurements['MCONT_ERR'][good], \
                        = SpectralIndices.apply_dispersion_corrections(
                                                    measurements['MCONT'][good],
                                                    measurements['MCONT_CORR'][good],
                                                    err=measurements['MCONT_ERR'][good])
                # Model continuum
                measurements['MCONT_MOD'][good], _ \
                        = SpectralIndices.apply_dispersion_corrections(
                                                    measurements['MCONT_MOD'][good],
                                                    measurements['MCONT_CORR'][good])


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
        map_hdr = DAPFitsUtil.build_map_header(self.binned_spectra.cube.fluxhdr,
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
        self.hdu = fits.HDUList([fits.PrimaryHDU(header=pri_hdr),
                                 fits.ImageHDU(data=bin_indx, header=map_hdr, name='BINID'),
                                 fits.ImageHDU(data=map_mask, header=map_hdr, name='MAPMASK'),
                                 fits.BinTableHDU.from_columns([fits.Column(name=n,
                                                    format=rec_to_fits_type(passband_database[n]),
                            array=passband_database[n]) for n in passband_database.dtype.names ],
                                                               name='SIPAR'),
                                  measurements.to_hdu(name='SINDX')])

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

        self.nbins = self.hdu['PRIMARY'].header['NBINS']
        unique_bins = numpy.unique(self.hdu['BINID'].data.ravel()) \
                            if self.hdu['PRIMARY'].header['SIREBIN'] else None
        self.missing_bins = self._get_missing_bins(unique_bins=unique_bins)

