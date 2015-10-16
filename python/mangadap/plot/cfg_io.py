# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""Methods for reading and manipulating config files."""

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import copy

import pandas as pd
from mangadap.plot import plotdap

try:
    from ConfigParser import RawConfigParser
except ImportError:
    from configparser import RawConfigParser

def read_config(filename):
    """Read in config file.

    Parses comma-separated strings or multi-line strings as lists.

    Args:
        filename (str): Full path to file.

    Returns:
        dict
    """
    config = RawConfigParser()
    config.read(filename)
    d = {}
    for section in config.sections():
        d[section] = {}
        for k in config.options(section):
            tmp = config.get(section, k)
            if ',' in tmp:
                tmp = [it.strip() for it in tmp.split(',')]
            elif '\n' in tmp:
                tmp = [it.strip() for it in tmp.split('\n')]
            d[section][k] = tmp
    return d

def read_plottypes_config(filename):
    """Read in plottypes config file.

    Args:
        filename (str): Full path to file.

    Returns:
        list
    """
    ptypes = read_config(filename)
    plottypes = ptypes['general']['plottypes']
    if isinstance(plottypes, str):
        plottypes = [plottypes]
    return plottypes


# FIX: is there a package to parse dtypes? (e.g., json.loads)

def string_to_float(d):
    """Convert values in a dictionary from strings to floats."""
    for k, v in d.items():
        d[k] = float(v)
    return d

def convert_to_number(s):
    try:
        s = int(s)
    except (TypeError, ValueError):
        try:
            s = float(s)
        except (TypeError, ValueError):
            raise TypeError
    return s

def convert_to_number_list(item):
    try:
        item = [int(it) for it in item]

    except (TypeError, ValueError):
        try:
            item = [float(it) for it in item]
        except (TypeError, ValueError):
            pass
    finally:
        return item

def convert_dtype(item):
    """Convert value from string to boolean or None."""
    try:
        if item in ('True', 'true', 'False', 'false'):
            item = item.lower() in ('True', 'true')
        elif item in ('None', 'none'):
            item = None
        item = convert_to_number(item)
    except (TypeError, ValueError):
        try:
            item = convert_to_number_list(item)
        except (TypeError, ValueError):
            pass
        finally:
            return item
    finally:
        return item

def specind_units(dapdata, si):
    """Get units for spectral index measurements.

    Args:
        dapdata: dap.DAP object.
        si (list): Spectral indices.
    Returns:
        Series
    """
    si_tex = [plotdap.pretty_specind_units(dapdata.sipar.unit[s]) for s in si]
    return pd.Series(si_tex, index=si)

def merge_dicts(d1, d2, sections=None):
    """Append keys and values from one dictionary into another.

    Args:
        d1 (dict): All keys and values from this dictionary will be in output.
        d2 (dict): Config parser output dictionary (two-level dictionary)
        sections (list): Top-level keys from d2 to merge into d1. Defaults to
            None.

    Returns:
       dict
    """
    out = copy.deepcopy(d1)
    if sections is None:
        sections = d2.keys()
    for section in sections:
        for k, v in d2[section].items():
            out[k] = v
    return out

def cblabels_to_series(d):
    """Convert colorbar labels to Series.

    If cblabels is a string, first convert it to a list of length
    len(columns). Then convert the list to a Series.

    Args:
        d (dict): Config parser output dictionary.

    Returns:
        Series
    """
    if isinstance(d['colorbar']['cblabels'], str):
        labels = [d['colorbar']['cblabels'] for _ in d['data']['columns']]
    else:
        labels = d['colorbar']['cblabels']
    return pd.Series(labels, index=d['data']['columns'])

def convert_config_dtypes(d, plottype=None, dapdata=None):
    """Convert dtypes of config parser dictionary to useful forms.

    Args:
        d (dict): Config parser output dictionary.
        plottype (str): Type of plot. Defaults to None.
        dapdata: dap.DAP object. Defaults to None.

    Returns:
        dict
    """
    if 'columns' in d['data']:
        columns = d['data']['columns']

    if 'titles' in d['data']:
        d['data']['titles'] = pd.Series(d['data']['titles'], index=columns)
    
    if 'colorbar' in d:
        if plottype is 'specind':
            d['colorbar']['cblabels'] = specind_units(dapdata, columns)
        else:
            d['colorbar']['cblabels'] = cblabels_to_series(d)
    
    if 'snr' in d['data']:
        d['snr'] = string_to_float(d['snr'])

    dout = {}
    for section in d:
        for k, v in d[section].items():
            dout[k] = convert_dtype(v)
    return dout
