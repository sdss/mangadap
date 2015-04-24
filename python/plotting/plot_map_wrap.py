#!/usr/bin/env python -W ignore::DeprecationWarning

"""
FILE
    plot_map_wrap.py

DESCRIPTION

    Generate a map of a given quantity.
   
    Usage: python plot_map_wrap.py qa_file_list.txt

"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import os
import sys
import errno
import copy
import numpy as np

from imp import reload

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.backends.backend_pdf import PdfPages

import plot_qa
from plot_qa import PlotQA
from plot_qa import FitError


dict_tmp = {}
try:
    dict_tmp.iteritems
except AttributeError:
    def iterkeys(d):
        return iter(d.keys())
    def itervalues(d):
        return iter(d.values())
    def iteritems(d):
        return iter(d.items())
else:
    def iterkeys(d):
        return d.iterkeys()
    def itervalues(d):
        return d.itervalues()
    def iteritems(d):
        return d.iteritems()

try:
    file_list = sys.argv[1]
    parname = sys.argv[2]
    if '-diverging' in sys.argv:
        diverging_cmap = True
    else:
        diverging_cmap = False
except:
    print('Usage: python plot_map_wrap.py qa_file_list.txt parname [-diverging]')
    sys.exit(2)



#----- Set Path and File Names -----
home = os.path.expanduser('~')
manga_dap_ver = os.getenv('MANGADAP_VER')
path_analysis = os.getenv('MANGA_SPECTRO_ANALYSIS')
path_analysis_plots = os.getenv('MANGA_SPECTRO_ANALYSIS_PLOTS')
path_dap_ver = '/'.join([path_analysis, manga_drp_ver, manga_dap_ver, ''])
path_file_list = path_dap_ver

if test_dir:
    path_dap_ver = '/'.join([path_analysis, manga_dap_ver, ''])
    path_dap_ver_plots = '/'.join([path_analysis_plots, manga_dap_ver, ''])
    path_file_list = path_dap_ver_plots


print('Input file list directory: %s' % (path_file_list))
print()
print('Input file list: %s' % (file_list))
print()
files = np.genfromtxt(path_file_list + file_list, dtype='str')
files = np.atleast_1d(files)
#-----------------------------------



#----- Make plots for each DAP output file ------

for dap_file in files:
    
    #----- Read in galaxy ID and DAP parameters -----
        
    stem_file = dap_file.strip('.fits')
    done_file = stem_file + '.done'
    
    ig1, pid, ifudesign, mode_in, binning_type, exec_num = stem_file.split('-')
    
    manga_id = '-'.join([pid, ifudesign])
    mode = mode_in.split('_')[0]
    
    path_galaxy = path_dap_ver + ('/').join([pid, ifudesign, ''])
    if test_dir:
        path_galaxy_plots = path_dap_ver_plots + ('/').join([pid, ifudesign, ''])
    else:
        path_galaxy_plots = path_galaxy + 'plots/'
    
    # create plot directories if necessary
    try:
        os.makedirs(path_galaxy_plots)
        print('\nCreated directory: %s\n' % path_galaxy_plots)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise e
        pass
    
    #------------------------------------------------- 
    
    
    #----- Plot -----
    # reload(plot_qa)
    # from plot_qa import PlotQA
    
    qa = PlotQA(path_galaxy + dap_file)
    print(manga_id, binning_type)
    qa.select_wave_range()
    qa.set_axis_lims()
    #qa.ch1, qa.ch1_r = qa.define_custom_cubehelix(rot=-2, start=1, gamma=1)
    

    # change mapval to quantity of your choosing
    mapval = qa.__dict__[parname]
    
    
    if diverging_cmap:
        cmap = cm.RdBu_r
    else:
        cmap = qa.linearL
    
    
    kwargs = dict(cblabel=parname,
                  cmap=cmap,
                  title_text=parname,
                  nodots=True)
    
    
    qa.plot_map(mapval, **kwargs)
    fout =  ('_').join([stem_file, parname, 'map']) + '.png'
    plt.savefig(path_galaxy_plots + fout)
    print('Wrote: %s' % fout)
    
    #---------------

