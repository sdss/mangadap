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

from mangadap import dap_access
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
    # Utah MPL4 file
    file_list = join(os.getenv('MANGA_SPECTRO_ANALYSIS'),
                     os.getenv('MANGADRP_VER'), os.getenv('MANGADAP_VER'),
                     'full', '8082', '3702', 'CUBE_files_to_plot.txt')
    # Portsmouth MPL4 file
    # file_list = join(os.getenv('MANGA_SANDBOX_MPL4'),
    #                  os.getenv('MANGADRP_VER'), os.getenv('MANGADAP_VER'),
    #                  os.getenv('MANGADAP_TPL'), '7443', '1901',
    #                  'CUBE_files_to_plot.txt')
    plottypes_list = 'dapqa_plottypes.ini'

file_kws_all = util.read_file_list(file_list)
cfg_dir = util.make_config_path(plottypes_list)
paths_cfg = join(cfg_dir, 'sdss_paths.ini')

for file_kws in file_kws_all:
    path_data = util.make_data_path(paths_cfg, file_kws)
    gal = dap_access.DAPAccess(path_data, file_kws, paths_cfg, verbose=True)
    gal.get_all_ext()
    mg_kws = util.make_mg_kws(gal, file_kws)
    
    ptypes_list = util.check_h3_h4(gal.bintype, gal.binsn, plottypes_list)
    plottypes = cfg_io.read_plottypes_config(join(cfg_dir, ptypes_list))
    for plottype in plottypes:
        cfg = cfg_io.read_config(join(cfg_dir, plottype + '.ini'))
        plot_kws = cfg_io.make_kws(cfg)
        plot_kws['main'] = hasattr(main, '__file__')
        plotdap.make_plots(plottype=plottype, dapdata=gal, mg_kws=mg_kws,
                           plot_kws=plot_kws)


from mangadap.plot import plotdap
reload(plotdap)
cfg = cfg_io.read_config(join(cfg_dir, plottype + '.ini'))
plot_kws = cfg_io.make_kws(cfg)
plot_kws['main'] = hasattr(main, '__file__')
plot_kws['columns'] = ['Ha6564']
plot_kws['make_multi'] = False
plotdap.make_plots(plottype=plottype, dapdata=gal, mg_kws=mg_kws, plot_kws=plot_kws)


x = 10.**(MaxNLocator(7).tick_values(np.log10(0.1), np.log10(12.5)))

LogLocator(base=10, subs=[1.0], numdecs=1, numticks=7).tick_values(3.1, 12.5)

sub = [1.0, 2.0, 3.0, 6.0]
ax.xaxis.set_ticks(np.arange(0.15, 12.5, 0.712123))
ticker.LogLocator(subs=subs).tick_values(0.15, 12.5)

"""
The problem is that LogLocator is taking vmin and vmax to be 1e-2 and 1e3 and
won't let me specify them.
"""

# MPL-4 file
# python3 plotqa.py $MANGA_SPECTRO_ANALYSIS/$MANGADRP_VER/$MANGADAP_VER/full/7443/3702/CUBE_files_to_plot.txt dapqa_plottypes.ini

# MPL-4 Sandbox
# python3 plotqa.py $MANGA_SANDBOX_MPL4/$MANGADRP_VER/$MANGADAP_VER/mpl4_m11stelibzsol/7443/1901/CUBE_files_to_plot.txt dapqa_plottypes.ini
