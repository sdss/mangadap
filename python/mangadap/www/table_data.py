"""
Generate data to populate dapqa.html.
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

dap_modes = ['CUBE']
binnings = ['STON-001', 'NONE-002'] #, 'RADIAL-003', 'RADIAL-004', 'ALL-005', 'ALL-006', 'ALL-007']

table = []

for plate in plates:
    ifus = list_dirs(os.path.join(path_data, plate))
    for ifu in ifus:
       for dap_mode in dap_modes:
           for binning in binnings:
            table.append([plate, ifu, dap_mode, binning])

for i in range(10):
    print(table[i])
                 

f = open(os.path.join(path_data, 'table_data.txt'), 'w')
f.write('[')
for i, item in enumerate(table):
    f.write("{'plate':%s, 'ifudsgn':'%s', 'mode':'%s', 'binning':'%s', 'spec':'0000', 'plotname':'maps'}" % tuple(item))
    if i != len(table)-1:
        f.write(',\n')

f.write(']')
f.close()


table_json = []
for item in table:
    table_json.append({'plate':item[0], 'ifudsgn':item[1], 'mode':item[2], 'binning':item[3], 'spec':'0000', 'plotname':'maps'})


with open(os.path.join(path_data, 'table_data.json'), 'w') as outfile:
    json.dump(table_json, outfile, indent=4)
