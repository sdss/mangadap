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

from plot import util

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

path_data = join(os.getenv('MANGA_SPECTRO_ANALYSIS'),
                 os.getenv('MANGADAP_VER'))
plates = list_dirs(path_data)

table_json = []
for plate in plates:
    ifus = list_dirs(join(path_data, plate))
    for ifu in ifus:
        pifu = '-'.join((plate, ifu))
        qual = drp3qual[plateifu == pifu][0]
        table_json.append({'plateifu':pifu, 'drp3qual':qual})


with open(os.path.join(path_data, 'table_drpqa_mpl4.json'), 'w') as outfile:
    json.dump(table_json, outfile, indent=4)
