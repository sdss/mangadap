from __future__ import division, print_function, absolute_import

import copy

import pandas as pd
import plotdap

try:
    from ConfigParser import RawConfigParser
except ImportError:
    from configparser import RawConfigParser

def read_config(filename):
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

def string_to_float(d):
    """Convert values in a dictionary from strings to floats.""" 
    for k, v in d.items():
        d[k]= float(v)
    return d

def specind_units(dapdata, si):
    si_tex = [plotdap.pretty_specind_units(dapdata.sipar.unit[s]) for s in si]
    return pd.Series(si_tex, index=si)

def merge_dicts(d1, d2, sections=None):
    """Append keys and values ."""
    out = copy.deepcopy(d1)
    if sections is None:
        sections = d2.keys()
    for section in sections:
        for k, v in d2[section].items():
            out[k] = v
    return out

def cblabels_to_series(d):
    """
    If cblabels is a string, first convert it to a list of length
    len(columns). Then convert the list to a Series.
    """
    if type(d['colorbar']['cblabels']) is str:
        labels = [d['colorbar']['cblabels'] for _ in d['data']['columns']]
    else:
        labels = d['colorbar']['cblabels']
    return pd.Series(labels, index=d['data']['columns'])


def convert_config_dtypes(mp_kws, d, plottype=None, dapdata=None):
    columns = d['data']['columns']
    d['snr'] = string_to_float(d['snr'])
    d['data']['titles'] = pd.Series(d['data']['titles'], index=columns)
    if plottype is 'specind':
        d['colorbar']['cblabels'] = specind_units(dapdata, columns)
    else:
        d['colorbar']['cblabels'] = cblabels_to_series(d)
    return merge_dicts(mp_kws, d)
