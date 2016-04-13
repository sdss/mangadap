# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""

Contains utility functions needed for configuration.

*License*:
    Copyright (c) 2015, SDSS-IV/MaNGA Pipeline Group
        Licensed under BSD 3-clause license - see LICENSE.rst

*Source location*:
    $MANGADAP_DIR/python/mangadap/config/util.py

*Imports and python version compliance*:
    ::

        from __future__ import division
        from __future__ import print_function
        from __future__ import absolute_import
        from __future__ import unicode_literals
        
        import sys
        if sys.version > '3':
            long = int
            try:
                from configparser import ConfigParser
            except ImportError:
                print('WARNING: Unable to import configparser!  Beware!')
        else:
            try:
                from ConfigParser import ConfigParser
            except ImportError:
                print('WARNING: Unable to import ConfigParser!  Beware!')
    
        import numpy
        import os.path
    
        from .defaults import default_dap_source

*Revision history*:
    | **29 Jan 2016**: Original implementation by K. Westfall (KBW)
    | **17 Feb 2016**: (KBW) Added try/except block for ConfigParser
    | **16 Mar 2016**: (KBW) Moved :func:`_read_dap_mask_bits` here;
        used to be in :mod:`mangadap.util.bitmask`

.. _configparser.ConfigParser: https://docs.python.org/3/library/configparser.html#configparser.ConfigParser

"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import sys
if sys.version > '3':
    long = int
    try:
        from configparser import ConfigParser
    except ImportError:
        print('WARNING: Unable to import configparser!  Beware!')
else:
    try:
        from ConfigParser import ConfigParser
    except ImportError:
        print('WARNING: Unable to import ConfigParser!  Beware!')

import numpy
import os.path

from .defaults import default_dap_source
from ..proc.spectralstack import SpectralStack

__author__ = 'Kyle B. Westfall'


def validate_reduction_assessment_config(cnfg):
    """ 
    Validate the `configparser.ConfigParser`_ object that is meant to
    define a reduction assessment method.

    Args:
        cnfg (`configparser.ConfigParser`_): Object meant to contain
            defining parameters of the reduction assessment method as
            needed by
            :class:`mangadap.proc.reductionassessments.ReductionAssessmentsDef`.

    Returns:
        bool: Booleans that specify how the reduction assessment should
        be constructed.  The flags specify to use (1) the wavelength
        range, (2) a bandpass filter parameter file, or (3) a file with
        a filter response function.

    Raises:
        KeyError: Raised if required keyword does not exist.
        ValueError: Raised if keys have unacceptable values.
        FileNotFoundError: Raised if a file is specified but could not
            be found.
    """
    # Check for required keywords
    if 'key' not in cnfg.options('default'):
        raise KeyError('Keyword \'key\' must be provided.')
    if 'wave_limits' not in cnfg.options('default') \
        and 'par_file' not in cnfg.options('default') \
        and 'response_function_file' not in cnfg.options('default'):
        raise KeyError('Method undefined.  Must provide \'wave_limits\' or \'par_file\' '
                       'or \'response_function_file\'.')

    def_range = ('wave_limits' in cnfg.options('default') \
                    and cnfg['default']['wave_limits'] is not None)

    def_par = ('par_file' in cnfg.options('default') and cnfg['default']['par_file'] is not None)

    def_response = ('response_function_file' in cnfg.options('default') \
                        and cnfg['default']['response_function_file'] is not None)

    if numpy.sum([ def_range, def_par, def_response] ) != 1:
        raise ValueError('Method undefined.  Must provide one and only one of \'wave1\' and '
                         '\'wave2\' or \'par_file\' or \'response_function_file\'.')

    if def_par and not os.path.isfile(cnfg['default']['par_file']):
        raise FileNotFoundError('par_file does not exist: {0}'.format(cnfg['default']['par_file']))
    if def_response and not os.path.isfile(cnfg['default']['response_function_file']):
        raise FileNotFoundError('response_function_file does not exist: {0}'.format(
                                cnfg['default']['response_function_file']))

    return def_range, def_par, def_response


def validate_spatial_binning_scheme_config(cnfg):
    """ 
    Validate the `configparser.ConfigParser`_ object that is meant to
    define a spatial binning scheme.

    Args:
        cnfg (`configparser.ConfigParser`_): Object meant to contain
            defining parameters of the binning method as needed by
            :class:`mangadap.proc.spatiallybinnedspectra.SpatiallyBinnedSpectraDef`.

    Returns:
        str: Name of the binning method to apply.

    Raises:
        KeyError: Raised if any required keywords do not exist.
        ValueError: Raised if keys have unacceptable values.
        FileNotFoundError: Raised if a file is specified but could not
            be found.
    """
    # Check for required keywords
    if 'key' not in cnfg.options('default'):
        raise KeyError('Keyword \'key\' must be provided.')
    if 'method' not in cnfg.options('default'):
        raise KeyError('Keyword \'method\' must be provided.')

    if 'minimum_snr' not in cnfg.options('default') or cnfg['default']['minimum_snr'] is None:
        cnfg['default']['minimum_snr']= '0.0'

    if 'operation' not in cnfg.options('default') or cnfg['default']['operation'] is None:
        cnfg['default']['operation'] = 'mean'

    if 'velocity_register' not in cnfg.options('default') \
            or cnfg['default']['velocity_register'] is None:
        cnfg['default']['velocity_register'] = 'False'

    if 'stack_covariance_mode' not in cnfg.options('default') \
            or cnfg['default']['stack_covariance_mode'] is None:
        cnfg['default']['stack_covariance_mode'] = 'none'

    covar_par_needed_modes = SpectralStack.covariance_mode_options(par_needed=True)
    if cnfg['default']['stack_covariance_mode'] in covar_par_needed_modes \
            and ('stack_covariance_par' not in cnfg.options('default') \
                    or cnfg['default']['stack_covariance_par'] is None):
        raise ValueError('For covariance mode = {0}, must provide a parameter!'.format(
                         cnfg['default']['stack_covariance_mode']))
        
    if 'stack_covariance_par' not in cnfg.options('default') \
            or cnfg['default']['stack_covariance_par'] is None:
        cnfg['default']['stack_covariance_par'] = 'none'

    if cnfg['default']['method'] in [ 'none', 'global' ]:
        return

    if cnfg['default']['method'] == 'voronoi':
        if 'target_snr' not in cnfg.options('default'):
            raise KeyError('Keyword \'target_snr\' must be provided for Voronoi binning.')
        return

    if cnfg['default']['method'] == 'radial':
        if 'center' not in cnfg.options('default'):
            raise KeyError('Keyword \'center\' must be provided for radial binning.')
        if 'pa' not in cnfg.options('default'):
            raise KeyError('Keyword \'pa\' must be provided for radial binning.')
        if 'ell' not in cnfg.options('default'):
            raise KeyError('Keyword \'ell\' must be provided for radial binning.')
        if 'radii' not in cnfg.options('default'):
            raise KeyError('Keyword \'radii\' must be provided for radial binning.')

        if 'radius_scale' not in cnfg.options('default') or cnfg['default']['radius_scale'] is None:
            cnfg['default']['radius_scale'] = '1.0'
        if 'log_step' not in cnfg.options('default') or cnfg['default']['log_step'] is None:
            cnfg['default']['log_step'] = 'False'
        return

    raise ValueError('{0} is not a recognized binning method.'.format(cnfg['default']['method']))


def validate_spectral_template_config(cnfg):
    """ 
    Validate the `configparser.ConfigParser`_ object that is meant to
    define a template library.

    Args:
        cnfg (`configparser.ConfigParser`_): Object meant to contain
            defining parameters of the template library needed by
            :class:`mangadap.proc.templatelibrary.TemplateLibraryDef`.

    Raises:
        KeyError: Raised if required keyword does not exist.
        ValueError: Raised if key has unacceptable value.

    """
    # Check for required keywords
    if 'key' not in cnfg.options('default'):
        raise KeyError('Keyword \'key\' must be provided.')
    if 'file_search' not in cnfg.options('default'):
        raise KeyError('Keyword \'file_search\' must be provided.')
    if 'fwhm' not in cnfg.options('default') and 'sres_ext' not in cnfg.options('default'):
        raise KeyError('Must provided keyword \'fwhm\' or \'sres_ext\'.')

    # Put in default values
    if 'in_vacuum' not in cnfg.options('default') or cnfg['default']['in_vacuum'] is None:
        cnfg['default']['in_vacuum'] = 'False'
    if 'wave_limit' not in cnfg.options('default') or cnfg['default']['wave_limit'] is None:
        cnfg['default']['wave_limit'] = 'None, None'
    if 'lower_flux_limit' not in cnfg.options('default') \
      or cnfg['default']['lower_flux_limit'] is None:
        cnfg['default']['lower_flux_limit'] = 'None'
    if 'log10' not in cnfg.options('default') or cnfg['default']['log10'] is None:
        cnfg['default']['log10'] = 'False'


def validate_emission_line_config(cnfg):
    """ 

    Validate the `configparser.ConfigParser`_ object that is meant to
    define an emission-line database.  `'in_vacuum'` can be present in
    the config file; however, it's not necessary.

    Args:
        cnfg (`configparser.ConfigParser`_): Object meant to contain
            defining parameters of the emission-line database needed by
            :class:`mangadap.proc.emissionlinedb.EmissionLineDBDef`.

    Raises:
        KeyError: Raised if required keyword does not exist.
        ValueError: Raised if key has unacceptable value.

    """
    # Check for required keywords
    if 'key' not in cnfg.options('default'):
        raise KeyError('Keyword \'key\' must be provided.')
    if 'file_path' not in cnfg.options('default'):
        raise KeyError('Keyword \'file_path\' must be provided.')


def validate_emission_bandpass_filter_config(cnfg):
    """ 

    Validate the `configparser.ConfigParser`_ object that is meant to
    define an emission-line passband database.  This is currently
    identical to :func:`validate_emission_line_config`.  `'in_vacuum'`
    can be present in the config file; however, it's not necessary.

    Args:
        cnfg (`configparser.ConfigParser`_): Object meant to contain
            defining parameters of the emission-line passband database
            needed by
            :class:`mangadap.proc.emissionlinedb.EmissionMomentsDBDef`.

    Raises:
        KeyError: Raised if required keyword does not exist.
        ValueError: Raised if key has unacceptable value.
    """
    # Check for required keywords
    if 'key' not in cnfg.options('default'):
        raise KeyError('Keyword \'key\' must be provided.')
    if 'file_path' not in cnfg.options('default'):
        raise KeyError('Keyword \'file_path\' must be provided.')


def validate_absorption_index_config(cnfg):
    """ 
    Validate the `configparser.ConfigParser`_ object that is meant to
    define an absorption-line index database.  This is currently
    identical to :func:`validate_emission_line_config`; however,
    `'in_vacuum'` should **not** be a keyword.

    Args:
        cnfg (`configparser.ConfigParser`_): Object meant to contain
            defining parameters of the absorption-line index database
            needed by
            :class:`mangadap.proc.absorptionindexdb.AbsorptionIndexDBDef`.

    Raises:
        KeyError: Raised if required keyword does not exist.
        ValueError: Raised if key has unacceptable value.
    """
    # Check for required keywords
    if 'key' not in cnfg.options('default'):
        raise KeyError('Keyword \'key\' must be provided.')
    if 'file_path' not in cnfg.options('default'):
        raise KeyError('Keyword \'file_path\' must be provided.')


def validate_bandhead_index_config(cnfg):
    """ 
    Validate the `configparser.ConfigParser`_ object that is meant to
    define a bandhead-index database.  This is currently identical to
    :func:`validate_absorption_index_config`.

    Args:
        cnfg (`configparser.ConfigParser`_): Object meant to contain
            defining parameters of the bandhead index database needed by
            :class:`mangadap.proc.bandheadindexdb.BandheadIndexDBDef`.

    Raises:
        KeyError: Raised if required keyword does not exist.
        ValueError: Raised if key has unacceptable value.
    """
    # Check for required keywords
    if 'key' not in cnfg.options('default'):
        raise KeyError('Keyword \'key\' must be provided.')
    if 'file_path' not in cnfg.options('default'):
        raise KeyError('Keyword \'file_path\' must be provided.')


"""
def read_dap_mask_bits(f, dapsrc=None):
    Read a config file that defines a set of mask bits to be used with
    :class:`mangadap.util.bitmask.BitMask`.

    Args:
        f (str): File name.  MUST be found in
            ${dapsrc}/python/mangadap/config/bitmasks/

        dapsrc (str): (Optional) Root of the DAP source directory tree.
            Default is to use
            :func:`mangadap.config.defaults.default_dap_source`

    Returns:
        keys (list): List of bitmask keywords.  They must be ordered by
            their bit value.
        descr (list): A list of descriptions for each keyword.

    Raise:
        FileNotFoundError: Raised if the file does not exist.
        ValueError: Raised if the list of bits do not have a sequential
            list of values or if the minimum mask value is not 0.

    dapsrc = default_dap_source() if dapsrc is None else str(dapsrc)

    _f = os.path.join(dapsrc, 'python', 'mangadap', 'config', 'bitmasks', f)
    if not os.path.isfile(_f):
        raise FileNotFoundError('File not found: {0}'.format(_f))

    cnfg = ConfigParser()
    cnfg.read(_f)

    keys = numpy.array(cnfg.sections())
    vals = numpy.zeros( keys.size, dtype=numpy.int )
    descr = numpy.zeros( keys.size, dtype=object )
    for i,k in enumerate(keys):
        vals[i] = cnfg[k]['value']
        descr[i] = cnfg[k]['descr']

    if numpy.amin(vals) != 0:
        raise ValueError('Bit values must start from 0.')
    if numpy.amax(vals) != keys.size-1:
        raise ValueError('Bit values must be sequential starting from 0.')
    srt = numpy.argsort(vals)

    return keys[srt], descr[srt]
"""
        
    

