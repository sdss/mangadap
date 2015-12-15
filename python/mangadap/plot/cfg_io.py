# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""Methods for reading and manipulating config files."""

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)
import ast

import pandas as pd
from mangadap.plot import plotdap

try:
    from ConfigParser import ConfigParser
except ImportError:
    from configparser import ConfigParser

def read_plottypes_config(filename):
    """Read in plottypes config file.

    Args:
        filename (str): Full path to file.

    Returns:
        list
    """
    ptypes = read_config(filename)
    plottypes = ptypes['plottypes']
    if isinstance(plottypes, str):
        plottypes = [plottypes]
    return plottypes

def tobool(inp):
    if inp in ('True', 'true', 'False', 'false'):
        return inp.lower() in ('True', 'true')

def tolist(inp):
    if ',' in inp:
        return [it.strip() for it in inp.split(',') if it != '']
    elif '\n' in inp:
        return [it.strip() for it in inp.split('\n') if it != '']
    else:
        return [inp]

def tolist_of_lists(inp):
    if '[' in inp:
        return list(ast.literal_eval(inp))

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

def read_config(filename):
    """Read in config file and assign dtypes.

    Args:
        filename (str): Full path to file.

    Returns:
        dict
    """
    config = ConfigParser()
    config.read(filename)
    d = {}
    for section in config.sections():
        for k in config.options(section):
            string = config.get(section, k)
            if section == 'str':
                d[k] = string
            elif section == 'list_str':
                d[k] = tolist(string)
            elif section == 'series_str':
                d[k] = pd.Series(tolist(string), index=d['columns'])
            elif section == 'int':
                d[k] = config.getint(section, k)
            elif section == 'list_int':
                d[k] = [int(it) for it in tolist(string)]
            elif section == 'series_int':
                d[k] = pd.Series([int(it) for it in tolist(string)],
                                 index=d['columns'])
            elif section == 'float':
                d[k] = config.getfloat(section, k)
            elif section == 'list_float':
                d[k] = [float(it) for it in tolist(string)]
            elif section == 'series_float':
                d[k] = pd.Series([float(it) for it in tolist(string)],
                                 index=d['columns'])
            elif section == 'bool':
                d[k] = config.getboolean(section, k)
            elif section == 'list_bool':
                d[k] = [tobool(it) for it in tolist(string)]
            elif section == 'series_bool':
                d[k] = pd.Series([tobool(it) for it in tolist(string)],
                                 index=d['columns'])
            elif section == 'series_list_float':
                d[k] = pd.Series(tolist_of_lists(string), index=d['columns'])
            elif section == 'list_list_str':
                d[k] = tolist_of_lists(string)
    return d

def make_kws(d):
    """Construct keyword arg dictionaries.

    Args:
        d (dict): plot_kws as read in from config file.

    Returns:
        dict
    """
    cb_kws_given = False
    for k in d:
        if 'cb-' in k:
            cb_kws_given = True

    if cb_kws_given:
        d['cb_kws'] = {}
        for k in list(d):
            if 'cb-' in k:
                k_cb = k.split('cb-')[-1]
                d['cb_kws'][k_cb] = d.pop(k)
    return d

