from __future__ import (division, print_function, absolute_import,
                        unicode_literals)
import os
from os.path import join
import copy
import numpy as np
import matplotlib.pyplot as plt

import pandas as pd
import seaborn as sns

from imp import reload

from mangadap import dap_access
from mangadap.plot import util
from mangadap.plot import cfg_io
from mangadap.stack import selections
from mangadap.stack import stack


# Set paths
path_mangadap = join(os.getenv('MANGADAP_DIR'), 'python', 'mangadap')
path_config = join(path_mangadap, 'stack', 'config')

# Read in meta-data sources
paths_cfg = join(path_mangadap, 'plot', 'config', 'sdss_paths.ini')
drpall = util.read_drpall(paths_cfg)
metadata_refs = dict(drpall=drpall)

# DO THIS VIA sdss_paths.ini
path_data = join(os.getenv('MANGA_SPECTRO_ANALYSIS'), os.getenv('MANGADRP_VER'),
                 os.getenv('MANGADAP_VER'), 'full')

# Read config file
cfg_name = 'test.ini'
cfg = cfg_io.read_config(join(path_config, cfg_name))
sample_conditions = [v for k, v in cfg.items() if 'sample_condition' in k]
bin_conditions = [v for k, v in cfg.items() if 'bin_condition' in k]
stack_values = [v for k, v in cfg.items() if 'stack_value' in k]

# READ THESE IN FROM CONFIG FILE
path_out = join(os.getenv('HOME'), 'tmp', 'dap_stack_out',
                os.getenv('MANGADRP_VER'), os.getenv('MANGADAP_VER'))
gal_ids = ['mangaid', 'RA', 'DEC']
value_names = ['Ha', 'D4000']
gal_columns = gal_ids + value_names


"""SAMPLE SELECTION"""
pifu_bool = selections.do_selection(sample_conditions, metadata_refs)
plateifus = drpall['plateifu'][pifu_bool]


"""LOOP OVER GALAXIES"""
objs = []
vals = []
# SPECIFY PLATES BY NUMBER VIA CONFIG FILE
plateifus = ['7443-3702', '7443-6102']
for plateifu in plateifus:
    """Get Bin Data"""
    filename = 'manga-{}-LOGCUBE_BIN-RADIAL-004.fits'.format(plateifu)
    file_kws = util.parse_fits_filename(filename)
    path_gal = join(path_data, file_kws['plate'], file_kws['ifudesign'])
    gal = dap_access.DAPAccess(path_gal, file_kws)
    gal.get_all_ext()
    galdata_refs = dict(dapdata=gal.__dict__)
    
    """Bin Selection"""
    bins_selected = selections.do_selection(bin_conditions, galdata_refs)
    bins_notnan = [selections.cfg_to_notnan(sv, galdata_refs)
                   for sv in stack_values]
    bins = selections.join_logical_and([bins_selected] + bins_notnan)
    
    """Save data"""
    objs.append([gal.mangaid, gal.header['OBJRA'], gal.header['OBJDEC']])
    vals.append([selections.cfg_to_data(sv, galdata_refs) for sv in stack_values])


"""Galaxy-internal combine"""
# from mangadap.stack import stack
# reload(stack)
vals_combined = [[stack.mean(col) for col in val] for val in vals]
indiv_gals = [obj + val for obj, val in zip(objs, vals_combined)]
galaxies = pd.DataFrame(indiv_gals, columns=gal_columns, index=plateifus)
galaxies.index.name = 'plateifu'

"""Cross-sample combine"""
vals_T = [list(it) for it in zip(*vals)] # "transpose" list
vals_concat = [pd.concat(it) for it in vals_T] # concat data by type
df_bin = pd.concat(vals_concat, axis=1, keys=value_names)
cross_sample = pd.Series([stack.mean(val) for val in vals_concat],
                         index=value_names)


"""Return Results"""
fout_stem = join(path_out, cfg_name.split('.ini')[0])
csvkwargs = dict(sep='\t', float_format='%10.5f')
fout_galaxies = fout_stem + '_gal.txt'
fout_cross_sample = fout_stem + '_sample.txt'

print('')
galaxies.to_csv(fout_galaxies, **csvkwargs)
print('Wrote {}'.format(fout_galaxies))
cross_sample.to_csv(fout_cross_sample, **csvkwargs)
print('Wrote {}'.format(fout_cross_sample))







