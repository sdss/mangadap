#!/usr/bin/env python -W ignore::DeprecationWarning

"""
FILE
    make_qa_file_list.py

DESCRIPTION

    Make a list of files that plot_qa_wrap.py will process.

    Usage: python make_qa_file_list.py path [output file name] [-overwrite]
   
"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import os
import sys
import errno
import numpy as np

print()

#----- Read in command line arguments -----
try:    
    path_out = sys.argv[1]
except:
    print('Usage: python make_qa_file_list.py path [output file name] [-overwrite]')

    
try:
    if sys.argv[2][0] is not '-':
        file_list = sys.argv[2]
    else:
        file_list = 'qa_file_list.txt'
except IndexError:
    file_list = 'qa_file_list.txt'


if '-overwrite' in sys.argv:
    overwrite = True
else:
    overwrite = False

#----------------------------------------


#----- Set Path ------
if path_out[-1] is not '/':
    path_out += '/'

print('Path: %s\n' % path_out)
#---------------------



#----- Read the file names -----
pid, ifudesign = path_out.split('/')[-3:-1]
file_stem = '-'.join(['manga', pid, ifudesign, 'LOG'])

files_tmp0 = os.listdir(path_out)
files_tmp1 = np.array([it for it in files_tmp0 if file_stem in it])


exec_num = np.array([it.split('.fits')[0][-3:] for it in files_tmp1])
ind_sort1 = np.argsort(exec_num)
files_tmp2 = files_tmp1[ind_sort1]
mode_id = np.array([it.split('LOG')[1].split('_')[0] for it in files_tmp2])
ind_sort2 = np.argsort(mode_id)
files_out = files_tmp2[ind_sort2]
#------------------------------


# #----- Read the file names -----
# files_out = []
# manga_pids = []
# 
# pids_in = os.listdir(path_out)
# pids = [it for it in pids_in if os.path.isdir(path_out + it)]
# 
# # get the IFUs from each plate
# for pid in pids:
#     ifudesigns_in = os.listdir('/'.join([path_out, pid]))
#     ifudesigns = [it for it in ifudesigns_in
#         if os.path.isdir('/'.join([path_out, pid, it]))]
#     for ifudesign in ifudesigns:
#         manga_pids.append('-'.join([pid, ifudesign]))
# 
# # get the DAP output file names from each galaxy
# for manga_pid in manga_pids:
#     file_stem = '-'.join(['manga', manga_pid, 'LOG'])
#     path_tmp = path_out + '/'.join(manga_pid.split('-') + [''])
#     files_in_tmp = os.listdir(path_tmp)
#     files_in = [it for it in files_in_tmp if os.path.isfile(path_tmp + it)]
#     files_tmp0 = np.array([it for it in files_in 
#                           if it[:len(file_stem)] == file_stem])
#     exec_num = np.array([it.split('.fits')[0][-3:] for it in files_tmp0])
#     ind_sort1 = np.argsort(exec_num)
#     files_tmp1 = files_tmp0[ind_sort1]
#     mode_id = np.array([it.split('LOG')[1].split('_')[0] for it in files_tmp1])
#     ind_sort2 = np.argsort(mode_id)
#     files = files_tmp1[ind_sort2]
#     for it in files:
#         files_out.append(it)
# 
# #------------------------------


#----- Create directories if necessary -----
try:
    os.makedirs(path_out)
    print('\nCreated directory: %s\n' % path_out)
except OSError as e:
    if e.errno != errno.EEXIST:
        raise e
    pass
#-------------------------------------------


#----- Write file list ------
if not overwrite:
    if not os.path.isfile(path_out + file_list):
        write = True
    else:
        write = False
else:
    write = True
    if os.path.isfile(path_out + file_list):
        print()
        print('Overwriting file list...')
        print()

if write:
    np.savetxt(path_out + file_list, files_out, fmt='%s')
    print('Wrote: %s' % path_out + file_list)
else:
    print('%s already exists. Use overwrite keyword to remake it.' % file_list)

#----------------------------