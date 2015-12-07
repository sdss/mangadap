"""
Generate data to populate Individual Galaxies table of
stkin_<Plate Century Group>.html pages.

Run as package from within /Users/andrewsb/manga/mangadap/trunk/python/mangadap/

python -m www.table_stkin_mpl4
"""
from __future__ import division, print_function, absolute_import

import os
from os.path import join
import re
import json

import numpy as np

from mangadap.plot import util

def list_dirs(path):
    out = []
    for it in os.listdir(path):
        if os.path.isdir(join(path, it)) & it.isdigit():
            out.append(it)
    return out

# Read in DRPall file
path_cfg = join(os.getenv('MANGADAP_DIR'), 'python', 'mangadap', 'plot',
                'config', 'sdss_paths.ini')
drpall = util.read_drpall(path_cfg)
plateifu = drpall['plateifu']

path_data = join(os.getenv('MANGA_SPECTRO_ANALYSIS'), os.getenv('MANGADRP_VER'),
                 os.getenv('MANGADAP_VER'), 'full')

plates = list_dirs(path_data)
#plates = ['7495']

table_json = {}
for plate in plates:
    plategroup = plate[:2] + '00'
    if plategroup not in table_json:
        table_json[plategroup] = []
    ifus_str = np.array(list_dirs(join(path_data, plate)))
    #ifus_str = np.array(['12701', '12702', '12703', '12704', '12705', '1901', '1902',
    #        '3701', '3702', '3703', '3704', '6101', '6102', '6103', '6104',
    #        '9101', '9102'])
    ifus_int = [int(ifu) for ifu in ifus_str]
    ifus_argsort = np.argsort(ifus_int)
    ifus_sort = ifus_str[ifus_argsort]
    for ifu in ifus_sort:
        pifu = '-'.join((plate, ifu))
        table_json[plategroup].append({'plateifu':pifu})

path_out = join(os.getenv('MANGA_SPECTRO_ANALYSIS'), 'test', 'andrews',
                os.getenv('MANGADRP_VER'), os.getenv('MANGADAP_VER'), 'full')
with open(os.path.join(path_out, 'table_stkin.json'), 'w') as outfile:
    json.dump(table_json, outfile, indent=4)
