#!/usr/bin/env python -W ignore::DeprecationWarning

"""
FILE
    make_qa_file_list.py

DESCRIPTION

    Make a list of files that plot_qa_wrap.py will process.

    Usage: python make_qa_file_list [-overwrite]
   
"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import os
import sys
import errno
import numpy as np



try:
    file_list = sys.argv[1]
    if '-overwrite' in sys.argv:
        overwrite = True
    else:
        overwrite = False
    if '-test_dir' in sys.argv:
        test_dir = True
    else:
        test_dir = False
except:
    file_list = 'qa_file_list.txt'
    overwrite = False
    test_dir = False


home = os.path.expanduser('~')
manga_dap_ver = os.getenv('MANGADAP_VER')
path_analysis = os.getenv('MANGA_SPECTRO_ANALYSIS')
path_analysis_plots = os.getenv('MANGA_SPECTRO_ANALYSIS_PLOTS')

if test_dir:
    path_dap_ver = '/'.join([path_analysis, manga_dap_ver, ''])
    path_dap_ver_plots = '/'.join([path_analysis_plots, manga_dap_ver, ''])
    path_file_list = path_dap_ver_plots
else:
    path_dap_ver = '/'.join([path_analysis, manga_drp_ver, manga_dap_ver, ''])
    path_file_list = path_dap_ver


files_out = []
manga_ids = []

pids_in = os.listdir(path_dap_ver)
pids = [it for it in pids_in if os.path.isdir(path_dap_ver + it)]

# get the IFUs from each plate
for pid in pids:
    ifudesigns_in = os.listdir('/'.join([path_dap_ver, pid]))
    ifudesigns = [it for it in ifudesigns_in
        if os.path.isdir('/'.join([path_dap_ver, pid, it]))]
    for ifudesign in ifudesigns:
        manga_ids.append('-'.join([pid, ifudesign]))

# get the DAP output file names from each galaxy
for manga_id in manga_ids:
    file_stem = '-'.join(['manga', manga_id, 'LOGCUBE_BIN'])
    path_tmp = '/'.join([path_dap_ver] + manga_id.split('-') + [''])
    files_in_tmp = os.listdir(path_tmp)
    files_in = [it for it in files_in_tmp if os.path.isfile(path_tmp + it)]
    files_tmp = np.array([it for it in files_in 
                         if it[:len(file_stem)] == file_stem])
    exec_num = np.array([it.split('.fits')[0][-3:] for it in files_tmp])
    ind_sort = np.argsort(exec_num)
    files = files_tmp[ind_sort]
    for it in files:
        files_out.append(it)


try:
    os.makedirs(path_file_list)
    print('\nCreated directory: %s\n' % path_file_list)
except OSError as e:
    if e.errno != errno.EEXIST:
        raise e
    pass


if not overwrite:
    if not os.path.isfile(path_file_list + file_list):
        write = True
    else:
        write = False
else:
    write = True
    if os.path.isfile(path_file_list + file_list):
        print()
        print('Overwriting file list...')
        print()

if write:
    np.savetxt(path_file_list + file_list, np.array(files_out), fmt='%s')
    print('Wrote: %s' % path_file_list + file_list)
else:
    print('%s already exists. Use overwrite keyword to remake it.' % file_list)
