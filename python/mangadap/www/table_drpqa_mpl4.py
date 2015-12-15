"""
Generate data to populate Individual Galaxies table of drpqa_mpl4.html.

Run as package from within /Users/andrewsb/manga/mangadap/trunk/python/mangadap/

python -m www.table_drpqa_mpl4
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
drp3qual = drpall['drp3qual']

path_data = join(os.getenv('MANGA_SPECTRO_ANALYSIS'), os.getenv('MANGADRP_VER'),
                 os.getenv('MANGADAP_VER'))
path_out = join(os.getenv('MANGA_SPECTRO_ANALYSIS'), 'test', 'andrews',
                os.getenv('MANGADRP_VER'), os.getenv('MANGADAP_VER'))

plates = list_dirs(path_data)
#plates = ['7495']

# table_json = []
# for plate in plates:
#     ifus = list_dirs(join(path_data, plate))
#     # ifus = ['12701', '12702', '12703', '12704', '12705', '1901', '1902',
#     #         '3701', '3702', '3703', '3704', '6101', '6102', '6103', '6104',
#     #         '9101', '9102']
#     for ifu in ifus:
#         pifu = '-'.join((plate, ifu))
#         qual = np.asscalar(drp3qual[plateifu == pifu][0])
#         table_json.append({'plateifu':pifu, 'drp3qual':qual})
# 
# path_out = join(os.getenv('MANGA_SPECTRO_ANALYSIS'), 'test', 'andrews',
#                 os.getenv('MANGADRP_VER'), os.getenv('MANGADAP_VER'))
# with open(os.path.join(path_out, 'table_drpqa_mpl4.json'), 'w') as outfile:
#     json.dump(table_json, outfile, indent=4)


table_json = []
page = 0
for i, plate in enumerate(plates):
    ifus = list_dirs(join(path_data, plate))
    # ifus = ['12701', '12702', '12703', '12704', '12705', '1901', '1902',
    #         '3701', '3702', '3703', '3704', '6101', '6102', '6103', '6104',
    #         '9101', '9102']
    for ifu in ifus:
        pifu = '-'.join((plate, ifu))
        qual = np.asscalar(drp3qual[plateifu == pifu][0])
        table_json.append({'plateifu':pifu, 'drp3qual':qual})
    if (i % 6 == 5) or (i == len(plates) - 1):
        fout = os.path.join(path_out, 'table_drpqa_mpl4_{}.json'.format(page))
        with open(fout, 'w') as outfile:
            json.dump(table_json, outfile, indent=4)
        table_json = []
        page += 1


