# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""

Contains utility functions needed for configuration.

*License*:
    Copyright (c) 2015, Kyle B. Westfall, Brett H. Andrews
    Licensed under BSD 3-clause license - see LICENSE.rst

*Source location*:
    $MANGADAP_DIR/python/mangadap/config/util.py

*Python2/3 compliance*::

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

*Imports*::

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

__author__ = 'Kyle B. Westfall'


def validate_spectral_template_config(cnfg):
    """ 
    Validate the `config.ConfigParser`_ object that is meant to define a
    template library.

    Args:
        cnfg (`config.ConfigParser`_): Object meant to contain defining
            parameters of the template library needed by
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
    Validate the `config.ConfigParser`_ object that is meant to define
    an emission-line database.

    ``'in_vacuum'`` can be provided; however, it's not necessary.

    Args:
        cnfg (`config.ConfigParser`_): Object meant to contain defining
            parameters of the emission-line database needed by
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
    Validate the `config.ConfigParser`_ object that is meant to define
    an emission-line passband database.  This is currently identical to
    :func:`validate_emission_line_config`.

    ``'in_vacuum'`` can be provided; however, it's not necessary.

    Args:
        cnfg (`config.ConfigParser`_): Object meant to contain defining
            parameters of the emission-line passband database needed by
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
    Validate the `config.ConfigParser`_ object that is meant to define
    an absorption-line index database.  This is currently identical to
    :func:`validate_emission_line_config`; however, ``'in_vacuum'``
    should **not** be a keyword.

    Args:
        cnfg (`config.ConfigParser`_): Object meant to contain defining
            parameters of the absorption-line index database needed by
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
    Validate the `config.ConfigParser`_ object that is meant to define a
    bandhead-index database.  This is currently identical to
    :func:`validate_absorption_index_config`.

    Args:
        cnfg (`config.ConfigParser`_): Object meant to contain defining
            parameters of the bandhead index database needed by
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


def read_dap_mask_bits(f, dapsrc=None):
    """
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
    """

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
        
    

