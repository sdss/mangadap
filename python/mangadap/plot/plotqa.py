# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""Plot DAP QA."""

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import os
from os.path import join
import sys
import copy

from imp import reload

import matplotlib as mpl
mpl.use('Agg')

from mangadap.plot import dap
from mangadap.plot import plotdap
from mangadap.plot import cfg_io
from mangadap.plot import util

import __main__ as main

if hasattr(main, '__file__'):
    # run as script
    try:
        file_list = sys.argv[1]
        plottypes_list = sys.argv[2]
    except IndexError:
        raise IndexError('Usage: python plotqa.py file_list plottypes_list')
else:
    # interactive session
    # DRPQA file
    # file_list = join(os.getenv('MANGA_SPECTRO_ANALYSIS'),
    #                  os.getenv('MANGADRP_VER'), os.getenv('MANGADAP_VER'),
    #                  '7495', '3704', 'CUBE_files_to_plot.txt')
    file_list = join(os.getenv('MANGA_MPL3'),
                     '7443', '1901', 'CUBE_files_to_plot.txt')
    plottypes_list = 'dapqa_plottypes.ini'

file_kws_all = util.read_file_list(file_list)
cfg_dir = util.make_config_path(plottypes_list)
paths_cfg = join(cfg_dir, 'sdss_paths.ini')

for file_kws in file_kws_all:
    path_data = util.make_data_path(paths_cfg, file_kws)
    gal = dap.DAP(path_data, paths_cfg, file_kws, verbose=False)
    gal.get_all_ext()
    mg_kws = util.make_mg_kws(gal, file_kws)

    plottypes = cfg_io.read_plottypes_config(join(cfg_dir, plottypes_list))
    for plottype in plottypes:
        cfg = cfg_io.read_config(join(cfg_dir, plottype + '.ini'))
        plot_kws = cfg_io.convert_config_dtypes(cfg, plottype, dapdata=gal)
        plotdap.make_plots(plottype=plottype, dapdata=gal,
                           mg_kws=mg_kws, plot_kws=plot_kws)


cfg = cfg_io.read_config(join(cfg_dir, plottype + '.ini'))
plot_kws = cfg_io.convert_config_dtypes(cfg, plottype, dapdata=gal)


# There should be a limited number of options.
# I would like to allow for None options or expanding a single option to
# multiple similar plots.

# assemble cb_kws after reading it in

# MOVE THIS FUNCTION TO cfg_io.py AND TEST IT

"""Read in config file and assign dtypes."""
config = ConfigParser()
config.read(join(cfg_dir, 'kinematics_maps_dtypes.ini'))
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
            d[k] =  pd.Series([tobool(it) for it in tolist(string)],
                              index=d['columns'])



def tobool(inp):
    if inp in ('True', 'true', 'False', 'false'):
        return inp.lower() in ('True', 'true')


def tolist(inp):
    if ',' in inp:
        return [it.strip() for it in inp.split(',') if inp != '']
    elif '\n' in inp:
        return [it.strip() for it in inp.split('\n') if inp != '']


# TO DO
# add option to pass in all cb_kws
# cfg_io can't read in dictionaries
# velocity: symmetric
# must test (esp, flux_band_maps) with MPL4 FITS files
# interactive mode that doesn't kill plots


# DRPQA file
# python3 plotqa.py $MANGA_MPL4/$MANGADRP_VER/$MANGADAP_VER/7443/6104/CUBE_files_to_plot.txt drpqa_plottypes.ini

# Normal DAP file
# python3 plotqa.py $MANGA_MPL4/$MANGADRP_VER/$MANGADAP_VER/7443/1901/CUBE_files_to_plot.txt drpqa_plottypes.ini
# python3 plotqa.py $MANGA_MPL3/7443/1901/CUBE_files_to_plot.txt dapqa_plottypes.ini


