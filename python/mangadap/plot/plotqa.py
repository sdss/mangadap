from __future__ import division, print_function, absolute_import

import os
from os.path import join
import sys
import copy

from imp import reload

import dap
import plotdap
import cfg_io
import util

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
    file_list = join(os.getenv('MANGA_SPECTRO_ANALYSIS'), '7443', '1901',
                     'CUBE_files_to_plot.txt')
    plottypes_list = 'drpqa_plottypes.ini'


file_kws_all = util.read_file_list(file_list)
cfg_dir = util.make_config_path(plottypes_list)

#cfg_dir = join(os.getenv('MANGADAP_DIR'), 'python', 'mangadap', 'plot',
#                   'config')
paths_cfg = join(cfg_dir, 'sdss_paths.ini')

reload(dap)
reload(plotdap)
for file_kws in file_kws_all:
    path_data = util.make_data_path(paths_cfg, file_kws)

    # Read DAP file
    gal = dap.DAP(path_data, paths_cfg, file_kws)
    gal.get_all_ext()
    mg_kws = copy.deepcopy(file_kws)
    mg_kws['mangaid'] = gal.mangaid

    # plottypes = ['snr', 'flux_ew', 'specind']
    # plottypes = ['drpqa']
    plottypes = cfg_io.read_plottypes_config(join(cfg_dir, plottypes_list))

    for plottype in plottypes:
        cfg = cfg_io.read_config(join(cfg_dir, plottype + '.ini'))
        plot_kws = cfg_io.convert_config_dtypes(cfg, plottype, dapdata=gal)
        print(plot_kws)
        plotdap.make_plots(dapdata=gal, mg_kws=mg_kws, **plot_kws)





# TO DO
# kinematics (barfing on val_no_measure=0)
# specind
# binnum

# gradients (emflux, specind)
# spectra
# emline zoomins

# python plotqa.py $MANGA_SPECTRO_ANALYSIS/7443/1901/CUBE_files_to_plot.txt dapqa_plottypes.ini

"""
Traceback (most recent call last):
  File "plotqa.py", line 57, in <module>
    plotdap.make_plots(dapdata=gal, mg_kws=mg_kws, **plot_kws)
  File "/Users/andrews/manga/mangadap/trunk/python/mangadap/plot/plotdap.py", line 665, in make_plots
    snr_thresh=snr_thresh)
  File "/Users/andrews/manga/mangadap/trunk/python/mangadap/plot/plotdap.py", line 80, in make_image
    no_data = make_mask_no_data(im, no_measure)
  File "/Users/andrews/manga/mangadap/trunk/python/mangadap/plot/plotdap.py", line 123, in make_mask_no_data
    no_data = np.isnan(data)
TypeError: ufunc 'isnan' not supported for the input types, and the inputs could not be safely coerced to any supported types according to the casting rule ''safe''
"""

