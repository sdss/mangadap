"""
Generate data to populate Individual Galaxies table of dapqa.html.
"""
from __future__ import print_function

import os
import json

path_data = os.path.join(os.getenv('MANGA_SPECTRO_ANALYSIS'), 'trunk_mpl3')


def list_dirs(path):
    out = []
    for it in os.listdir(path):
        if os.path.isdir(os.path.join(path, it)):
            out.append(it)
    return out


plates = list_dirs(path_data)

dap_modes = ['CUBE', 'RSS']
binnings = {'CUBE': ['STON-001', 'NONE-002', 'RADIAL-003', 'RADIAL-004', 'ALL-005', 'ALL-006', 'ALL-007'], 
            'RSS': ['RADIAL-001', 'RADIAL-002', 'ALL-003', 'ALL-004', 'ALL-005']}

table = []

for plate in plates:
    ifus = list_dirs(os.path.join(path_data, plate))
    for ifu in ifus:
        for dap_mode in dap_modes:
            for binning in binnings[dap_mode]:
                table.append([plate, ifu, dap_mode, binning])


table_json = []
for item in table:
    if 'RADIAL' in item[3]:
        plotname = 'gradients'
    elif 'ALL' in item[3]:
        plotname = ''
    else:
        plotname = 'maps'
    table_json.append({'plate':item[0], 'ifudsgn':item[1], 'mode':item[2], 'binning':item[3], 'spec':'0000', 'plotname':plotname})


with open(os.path.join(path_data, 'table_data.json'), 'w') as outfile:
    json.dump(table_json, outfile, indent=4)
