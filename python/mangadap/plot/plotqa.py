from __future__ import division, print_function, absolute_import

import os
from os.path import join
import sys
import copy

from imp import reload

from sdss.files import base_path
import dap
import plotdap
import cfg_io
import util


if hasattr('__main__', '__file__'):
    # run as script
    try:
        file_list = sys.argv[1]
        plottypes_list = sys.argv[2]
    except IndexError:
        raise IndexError('Usage: python plotqa.py file_list plottypes_list')
else:
    # interactive session
    file_list = join(os.getenv('$MANGA_SPECTRO_ANALYSIS'), 'trunk_mpl3',
                     '7443', '12701', 'CUBE_files_to_plot.txt')
    plottypes_list = 'drpqa_plottypes.ini'


file_kws_all = util.read_file_list(file_list)
cfg_dir = util.make_config_path(plottypes_list)

reload(dap)
reload(plotdap)
for file_kws in file_kws_all:
    path_data = util.make_data_path(file_kws)

    # Read DAP file
    gal = dap.DAP(path_data, file_kws)
    gal.get_all_ext()
    mg_kws = copy.deepcopy(file_kws)
    mg_kws['mangaid'] = gal.mangaid
    
    # plottypes = ['snr', 'flux_ew', 'specind']
    # plottypes = ['drpqa']
    plottypes = cfg_io.read_plottypes_config(join(cfg_dir, plottypes_list))

    for plottype in plottypes:
        cfg = cfg_io.read_config(join(cfg_dir, plottype + '.ini'))
        plot_kws = cfg_io.convert_config_dtypes(cfg, plottype, dapdata=gal)
        plotdap.make_plots(dapdata=gal, mg_kws=mg_kws, **plot_kws)





# TO DO

# python plotqa.py $MANGA_SPECTRO_ANALYSIS/trunk_mpl3/7443/12701/CUBE_files_to_plot.txt drpqa_plottypes.ini
