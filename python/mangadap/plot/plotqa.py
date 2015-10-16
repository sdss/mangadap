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
#mpl.use('Agg')

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
    #file_list = join(os.getenv('MANGA_MPL4'), os.getenv('MANGADRP_VER'),
    #                 os.getenv('MANGADAP_VER'),
    #                 '7443', '6104', 'CUBE_files_to_plot.txt')
    file_list = join(os.getenv('MANGA_MPL3'),
                     '7443', '1901', 'CUBE_files_to_plot.txt')
    plottypes_list = 'drpqa_plottypes.ini'

file_kws_all = util.read_file_list(file_list)
cfg_dir = util.make_config_path(plottypes_list)
paths_cfg = join(cfg_dir, 'sdss_paths.ini')

for file_kws in file_kws_all:
    path_data = util.make_data_path(paths_cfg, file_kws)
    # Read DAP file
    gal = dap.DAP(path_data, paths_cfg, file_kws, verbose=False)
    gal.get_all_ext()
    mg_kws = util.make_mg_kws(gal, file_kws)


reload(cfg_io)
plottype = 'emline'
cfg = cfg_io.read_config(join(cfg_dir, plottype + '.ini'))
plot_kws = cfg_io.convert_config_dtypes(cfg, plottype, dapdata=gal)

plotdap.plot_emlines(dapdata=gal, bin=0, mg_kws=mg_kws)

plotdap.plot_emline_multi(dapdata=gal, bin=0, mg_kws=mg_kws)

reload(plotdap)
plotdap.plot_emline(dapdata=gal, bin=0, mg_kws=mg_kws, nii=False, win_cen=6565.)


    plottypes = cfg_io.read_plottypes_config(join(cfg_dir, plottypes_list))
    for plottype in plottypes:
        cfg = cfg_io.read_config(join(cfg_dir, plottype + '.ini'))
        plot_kws = cfg_io.convert_config_dtypes(cfg, plottype, dapdata=gal)
        plotdap.make_plots(dapdata=gal, mg_kws=mg_kws, **plot_kws)


# DRPQA file
# python plotqa.py $MANGA_MPL4/$MANGADRP_VER/$MANGADAP_VER/7443/6104/CUBE_files_to_plot.txt drpqa_plottypes.ini

# Normal DAP file
# python plotqa.py $MANGA_MPL4/$MANGADRP_VER/$MANGADAP_VER/7443/1901/CUBE_files_to_plot.txt drpqa_plottypes.ini


