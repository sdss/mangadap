#!/usr/bin/env python -W ignore::DeprecationWarning

"""
FILE
    plot_qa_wrap.py

DESCRIPTION

    Generate QA plots.
   
"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import os
import sys
import copy
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.backends.backend_pdf import PdfPages

import plot_qa
from plot_qa import PlotQA

try:
    manga_id = sys.argv[1]
except:
    print('Usage: %s manga_id' % sys.argv[0]); sys.exit(1)


#manga_id = '8131-6102'
pid, ifudesign = manga_id.split('-')

#----- Set Path and File Names -----
manga_dap_ver = os.getenv('MANGADAP_VER')
path_analysis = os.getenv('MANGA_SPECTRO_ANALYSIS')
path_dap_ver = '/'.join([path_analysis, manga_dap_ver, ''])
path_galaxy = path_dap_ver + ('/').join([pid, ifudesign, ''])
path_plots = path_galaxy + 'plots/'

if not os.path.isdir(path_plots):
    os.mkdir(path_plots)
    print('\nCreated directory: %s\n' % path_plots)

#-----------------------------------

#----- Analysis Parameters -----

# Make this more flexible.
# Input name of fits file?
# Fit all DAP output fits files?

mode = 'LOGCUBE'
binning_type = 'STON'
binning_num = '001'
#-------------------------------

#----- DAP Output File Name -----
galaxy_id = ('-').join(['manga', manga_id])
root_name = ('-').join([galaxy_id, mode])
binning_id = ('-').join(['BIN', binning_type, binning_num])
dap_file =  ('_').join([root_name, binning_id]) + '.fits'
#--------------------------------

#----- Calculate Fit Metrics -----
# Remove
reload(plot_qa)
from plot_qa import PlotQA

qa = PlotQA(path_galaxy + dap_file)
print('\nTemplate Library:', qa.tpl_lib)
qa.select_wave_range()
qa.set_axis_lims()
qa.calc_chisq()
qa.calc_resid()
qa.calc_resid_data()
#---------------------------------

#----- Plot Parameters -----

#--- Maps ---
mapname = [
    'chisq',
    'stvel',
    'stvdisp',
    'resid',
    'sth3',
    'sth4',]

chisq_args = dict(val=qa.chisq_bin,
                  kwargs=dict(cblabel=r'$\chi_{\rm red}^2$',
                              cmap=qa.cubehelix111,
                              title_text=r'$\chi_{\rm red}^2$'))
resid_args = dict(val=qa.resid_data_bin_percent99,
                  kwargs=dict(cblabel=r'99th percentile |resid| / galaxy',
                              cmap=qa.cubehelix111,
                              title_text='99th percentile |resid| / galaxy'))
stvel_args = dict(val=qa.stvel,
                  kwargs=dict(cblabel=r'v$_\star$ [km/s]',
                              cmap=cm.RdBu_r,
                              title_text=r'v$_\star$'))
stvdisp_args = dict(val=qa.stvdisp,
                    kwargs=dict(cblabel=r'$\sigma_\star$ [km/s]',
                                cmap=qa.cubehelix111, 
                                title_text=r'$\sigma_\star$'))
sth3_args = dict(val=qa.sth3,
                 kwargs=dict(cblabel='h3',
                             cbrange_symmetric=True,
                             cmap=cm.RdBu_r,
                             title_text='h3'))
sth4_args = dict(val=qa.sth4,
                 kwargs=dict(cblabel='h4',
                             cbrange_symmetric=True,
                             cmap=cm.RdBu_r,
                             title_text='h4'))

map_kwargs = dict(chisq=chisq_args,
                  resid=resid_args,
                  stvel=stvel_args,
                  stvdisp=stvdisp_args,
                  sth3=sth3_args,
                  sth4=sth4_args)

map_kwargs_interp = copy.deepcopy(map_kwargs)
for v in map_kwargs_interp.itervalues():
    v['kwargs']['interpolated'] = True

bin_num_map_kwargs = dict(cblabel=r'$\chi_{\rm red}^2$',
                          cmap=qa.cubehelix111,
                          title_text=r'$\chi_{\rm red}^2$',
                          nodots=True,
                          spaxel_num=True,
                          figsize=(15, 12))

#-----------


#--- Spectra ----

lam_lim = [3600, 7000]

# set the same min and max flux for each spectrum
p32 = np.array([np.percentile(qa.galaxy[i], 32) for i in xrange(qa.n_bins)])
p50 = np.array([np.percentile(qa.galaxy[i], 50) for i in xrange(qa.n_bins)])
p68 = np.array([np.percentile(qa.galaxy[i], 68) for i in xrange(qa.n_bins)])

sig1 = (np.sort(p68) - np.sort(p32))[-1] / 2.

fluxmin = 0.
fluxmax = np.max(p50) + 4. * sig1

residmin = -0.25
residmax = 0.25

resid_kwargs = dict(lw=1,
                    leg_lab=['fit'],
                    xlim=lam_lim,
                    ylim=[[fluxmin, fluxmax], [residmin, residmax]],
                    masks=True)

all_spec_obs_kwargs = dict(rest_frame=False,
                           xlim=lam_lim,
                           ylim=[fluxmin, 2.])
all_spec_rest_kwargs = dict(rest_frame=True,
                            xlim=lam_lim,
                            ylim=[fluxmin, 2.])

all_resid_obs_kwargs = dict(resid=True,
                            rest_frame=False,
                            xlim=lam_lim,
                            ylim=[residmin, residmax])

all_resid_rest_kwargs = dict(resid=True,
                             rest_frame=True,
                             xlim=lam_lim,
                             ylim=[residmin, residmax])
#----------------


#-----------------------------


#----- Plots -----
print('')


# reload(plot_qa)
# from plot_qa import PlotQA
# 
# qa = PlotQA(path_galaxy + dap_file)
# print('\nTemplate Library:', qa.tpl_lib)
# qa.select_wave_range()
# qa.set_axis_lims()
# qa.calc_chisq()
# qa.calc_resid()
# qa.calc_resid_data()


# Plot bin numbers on top of chisq map
qa.plot_map(qa.chisq_bin, **bin_num_map_kwargs)
fout =  ('_').join([root_name, binning_id, 'bin', 'num']) + '.pdf'
plt.savefig(path_plots + fout)
print('Wrote: %s' % fout)


# Plot binned maps of chisq, resid/galaxy, stellar kinematics
qa.plot_multi_map(mapname, map_kwargs)
fout =  ('_').join([root_name, binning_id, 'maps']) + '.png'
plt.savefig(path_plots + fout)
print('Wrote: %s' % fout)

# Plot interpolated maps of chisq, resid/galaxy, stellar kinematics
qa.plot_multi_map(mapname, map_kwargs_interp)
fout =  ('_').join([root_name, binning_id, 'maps', 'interp']) + '.png'
plt.savefig(path_plots + fout)
print('Wrote: %s' % fout)


# Overplot all residuals (observed frame)
qa.plot_all_spec(**all_resid_obs_kwargs)
fout =  ('_').join([root_name, binning_id, 'all', 'resid', 'obs']) + '.png'
plt.savefig(path_plots + fout, dpi=200)
print('Wrote: %s' % fout)


# Overplot all residuals (rest frame)
qa.plot_all_spec(**all_resid_rest_kwargs)
fout =  ('_').join([root_name, binning_id, 'all', 'resid', 'rest']) + '.png'
plt.savefig(path_plots + fout, dpi=200)
print('Wrote: %s' % fout)


# Overplot all spectra (observed frame)
qa.plot_all_spec(**all_spec_obs_kwargs)
fout =  ('_').join([root_name, binning_id, 'all', 'spec', 'obs']) + '.png'
plt.savefig(path_plots + fout, dpi=200)
print('Wrote: %s' % fout)


# Overplot all spectra (rest frame)
qa.plot_all_spec(**all_spec_rest_kwargs)
fout =  ('_').join([root_name, binning_id, 'all', 'spec', 'rest']) + '.png'
plt.savefig(path_plots + fout, dpi=200)
print('Wrote: %s' % fout)


# Plot all spectra individually (one file)
fout =  ('_').join([root_name, binning_id, 'resid', 'all', 'bins']) + '.pdf'
print('Writing: %s ...' % fout)
with PdfPages(path_plots + fout) as pdf:
    for bin in xrange(qa.n_bins):
        qa.plot_resid(bin=bin, **resid_kwargs)
        pdf.savefig()
        plt.close()
        plt.clf()
        if bin % 25 == 0:
            print('    bin %s...done' % bin)

plt.close()
print('Wrote: %s' % fout)



print('')
#--------------


