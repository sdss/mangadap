#!/usr/bin/env python -W ignore::DeprecationWarning

"""
FILE
    plot_qa_wrap.py

DESCRIPTION

    Generate QA plots.
   
    usage: python plot_qa_wrap.py [file_list] [-overwrite] 

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
    if '-overwrite' in sys.argv:
        overwrite = True
except:
    file_list = 'qa_file_list.txt'
    overwrite = False

print()
print('Overwrite plots:', overwrite)
print()

#----- Set Path and File Names -----
home = os.path.expanduser('~')
manga_dap_ver = os.getenv('MANGADAP_VER')
path_analysis = os.getenv('MANGA_SPECTRO_ANALYSIS')
path_analysis_plots = os.getenv('MANGA_SPECTRO_ANALYSIS_PLOTS')
path_dap_ver = '/'.join([path_analysis, manga_dap_ver, ''])
path_dap_ver_plots = '/'.join([path_analysis_plots, manga_dap_ver, ''])

print('Input file list directory: %s' % (path_dap_ver_plots))
print()
print('Input file list: %s' % (file_list))
print()
files = np.genfromtxt(path_dap_ver_plots + file_list, dtype='str')
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
    path_galaxy_plots = path_dap_ver_plots + ('/').join([pid, ifudesign, ''])
    
    # create plot directories if necessary
    try:
        os.makedirs(path_galaxy_plots)
        print('\nCreated directory: %s\n' % path_galaxy_plots)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise e
        pass
    
    # check for done file
    if not overwrite:
        if os.path.isfile(path_galaxy_plots + done_file):
            print('%-12s %-6s: Plots already exist.' % (manga_id, binning_type) + 
                  '  Use overwrite keyword to remake them')
            continue
    
    #------------------------------------------------- 
    
    
    #----- Plots for each binning type  ----- 
    if binning_type == 'NONE':
        plot_map = True
        overplot_all = True
        plot_all_spec_as_pdf = False 
    elif binning_type == 'STON':
        plot_map = True
        overplot_all = True
        plot_all_spec_as_pdf = True
    elif binning_type == 'RADIAL':
        plot_map = False
        overplot_all = True
        plot_all_spec_as_pdf = True
    elif binning_type == 'ALL':
        plot_map = False
        overplot_all = False
        plot_all_spec_as_pdf = True
    #----------------------------------------
    
    
    #----- Calculate Fit Metrics -----
    # Remove
    #reload(plot_qa)
    #from plot_qa import PlotQA
    try:
        qa = PlotQA(path_galaxy + dap_file)
        print(manga_id, binning_type)
        # print('Template Library:', qa.tpl_lib)
        qa.select_wave_range()
        qa.set_axis_lims()
        qa.calc_chisq()
        qa.calc_resid()
        qa.calc_resid_data()
        qa.ch1, qa.ch1_r = qa.define_custom_cubehelix(rot=-2, start=1, gamma=1)
    except FitError:
        print('\n%-12s %-6s: *NOT* fit by DAP\n' % (manga_id, binning_type))
        continue
    
    #---------------------------------
    
    
    
    #----- Plot Parameters -----
     
    #--- Emission Line Kinematics Map Parameters ---
    emkin_mapname = [
      'signal',
      'noise',
      'snr',
      'emvel',
      'emvdisp',
      'emvdisp_halpha',]
    
    signal_args = dict(val=qa.signal,
                      kwargs=dict(cblabel=r'signal',
                                  cmap=qa.ch1,
                                  title_text=r'DRPS signal',
                                  nodots=True))
    noise_args = dict(val=qa.noise,
                      kwargs=dict(cblabel=r'noise',
                                  cmap=qa.ch1,
                                  title_text='DRPS noise',
                                  nodots=True))
    snr_args = dict(val=qa.snr,
                        kwargs=dict(cblabel=r'S/N',
                                    cmap=qa.ch1, 
                                    title_text=r'DRPS S/N',
                                    nodots=True))
    emvel_args = dict(val=qa.emvel_ew,
                      kwargs=dict(cblabel=r'v$_{\rm gas}$ [km/s]',
                                  cmap=cm.coolwarm,
                                  title_text=r'v_gas (Enci)',
                                  nodots=True))
    emvdisp_args = dict(val=qa.emvdisp_ew,
                     kwargs=dict(cblabel=r'$\sigma_{\rm gas}$ [km/s]',
                                 cbrange=[35, 150],
                                 cmap=qa.ch1,
                                 title_text=r'$\sigma$_gas (Enci)',
                                 nodots=True))
    emvdisp_halpha_args = dict(val=qa.emvdisp_halpha_ew,
                     kwargs=dict(cblabel=r'$\sigma_{\rm inst}$ H$\alpha$ [km/s]',
                                 cmap=qa.ch1,
                                 title_text=r'$\sigma_{\rm inst}$ H$\alpha$ (Enci)',
                                 nodots=True))
    
    emkin_map_kwargs = dict(signal=signal_args,
                            noise=noise_args,
                            emvel=emvel_args,
                            snr=snr_args,
                            emvdisp=emvdisp_args,
                            emvdisp_halpha=emvdisp_halpha_args)
    #---------------------------------
    
    #--- Stellar Kinematics Map Parameters ---
    # stkin_mapname = [
    #   'chisq',
    #   'stvel',
    #   'stvdisp',
    #   'resid',
    #   'sth3',
    #   'sth4',]
    stkin_mapname = [
      'chisq',
      'stvel',
      'stvdisp',
      'resid',
      'emvel',
      'emvdisp',]
    
    chisq_args = dict(val=qa.chisq_bin,
                      kwargs=dict(cblabel=r'$\chi_{\rm red}^2$',
                                  cmap=qa.ch1,
                                  title_text=r'$\chi_{\rm red}^2$',
                                  nodots=True))
    resid_args = dict(val=qa.resid_data_bin_percent99,
                      kwargs=dict(cblabel='99th percentile |resid| / galaxy',
                                  cmap=qa.ch1,
                                  title_text='99th percentile |resid| / galaxy',
                                  nodots=True))
    stvel_args = dict(val=qa.stvel,
                      kwargs=dict(cblabel=r'$v_\star$ [km/s]',
                                  cmap=cm.coolwarm,
                                  title_text=r'v_star',
                                  nodots=True))
    stvdisp_args = dict(val=qa.stvdisp,
                        kwargs=dict(cblabel=r'$\sigma_\star$ [km/s]',
                                    cbrange=[35, 150],
                                    cmap=qa.ch1, 
                                    title_text=r'$\sigma$_star',
                                  nodots=True))
    # sth3_args = dict(val=qa.sth3,
    #                  kwargs=dict(cblabel='h3',
    #                              cbrange_symmetric=True,
    #                              cmap=cm.coolwarm,
    #                              title_text='h3',
    #                               nodots=True))
    # sth4_args = dict(val=qa.sth4,
    #                  kwargs=dict(cblabel='h4',
    #                              cbrange_symmetric=True,
    #                              cmap=cm.coolwarm,
    #                              title_text='h4',
    #                               nodots=True))
    
    # stkin_map_kwargs = dict(chisq=chisq_args,
    #                         resid=resid_args,
    #                         stvel=stvel_args,
    #                         stvdisp=stvdisp_args,
    #                         sth3=sth3_args,
    #                         sth4=sth4_args)
    
    stkin_map_kwargs = dict(chisq=chisq_args,
                            resid=resid_args,
                            stvel=stvel_args,
                            stvdisp=stvdisp_args,
                            emvel=emvel_args,
                            emvdisp=emvdisp_args)
    
    stkin_map_kwargs_interp = copy.deepcopy(stkin_map_kwargs)
    for v in itervalues(stkin_map_kwargs_interp):
        v['kwargs']['interpolated'] = True
        #v['kwargs']['cbrange_clip'] = False
    
    bin_num_map_kwargs = dict(cblabel=r'$\chi_{\rm red}^2$',
                              cmap=qa.ch1,
                              title_text=r'$\chi_{\rm red}^2$',
                              nodots=True,
                              spaxel_num=True,
                              figsize=(15, 12))
    
    #----------------------------------------
    
    
    #--- Emission Line Fluxes Map Parameters ---
    emflux_mapname = [
      'oii3727',
      'hbeta',
      'nii6583',
      'oiii5007',
      'halpha',
      'sii6717',]
    
    oii3727_args = dict(val=qa.oii3727_ew,
                        kwargs=dict(cblabel=r'Flux [10$^{-17}$ erg/s/cm$^2$]',
                                    cmap=qa.ch1,
                                    title_text=r'[OII] $\lambda$3727 (Enci)',
                                    nodots=True))
    hbeta_args = dict(val=qa.hbeta_ew,
                      kwargs=dict(cblabel=r'Flux [10$^{-17}$ erg/s/cm$^2$]',
                                  cmap=qa.ch1,
                                  title_text=r'H$\beta$ (Enci)',
                                  nodots=True))
    nii6583_args = dict(val=qa.nii6583_ew,
                        kwargs=dict(cblabel=r'Flux [10$^{-17}$ erg/s/cm$^2$]',
                                    cmap=qa.ch1, 
                                    title_text=r'[NII] $\lambda$6583 (Enci)',
                                    nodots=True))
    oiii5007_args = dict(val=qa.oiii5007_ew,
                         kwargs=dict(cblabel=r'Flux [10$^{-17}$ erg/s/cm$^2$]',
                                     cmap=qa.ch1,
                                     title_text=r'[OIII] $\lambda$5007 (Enci)',
                                     nodots=True))
    halpha_args = dict(val=qa.halpha_ew,
                       kwargs=dict(cblabel=r'Flux [10$^{-17}$ erg/s/cm$^2$]',
                                   cmap=qa.ch1,
                                   title_text=r'H$\beta$ (Enci)',
                                   nodots=True))
    sii6717_args = dict(val=qa.sii6717_ew,
                     kwargs=dict(cblabel=r'Flux [10$^{-17}$ erg/s/cm$^2$]',
                                 cmap=qa.ch1,
                                 title_text=r'[SII] $\lambda$6717 (Enci)',
                                 nodots=True))
    
    
    emflux_map_kwargs = dict(oii3727=oii3727_args,
                             hbeta=hbeta_args,
                             oiii5007=oiii5007_args,
                             nii6583=nii6583_args,
                             halpha=halpha_args,
                             sii6717=sii6717_args)
    #---------------------------------
    
    
    
    #--- Spectra Parameters ----
    
    lam_lim = [3600, 9300]
    
    # set the same min and max flux for each spectrum
    p32 = np.array([np.percentile(qa.galaxy[i], 32) for i in range(qa.n_bins)])
    p50 = np.array([np.percentile(qa.galaxy[i], 50) for i in range(qa.n_bins)])
    p68 = np.array([np.percentile(qa.galaxy[i], 68) for i in range(qa.n_bins)])
    
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
    # qa.ch1, qa.ch1_r = qa.define_custom_cubehelix(rot=-2, start=1, gamma=1)
    
    
    
    #plot_map = False
    if plot_map:
    
        # Plot bin numbers on top of chisq map
        qa.plot_map(qa.chisq_bin, **bin_num_map_kwargs)
        fout =  ('_').join([stem_file, 'bin', 'num']) + '.png'
        plt.savefig(path_galaxy_plots + fout)
        print('Wrote: %s' % fout)
        
        # Plot binned maps of chisq, resid/galaxy, stellar kinematics
        qa.plot_multi_map(stkin_mapname, stkin_map_kwargs)
        fout =  ('_').join([stem_file, 'stkin', 'maps']) + '.png'
        plt.savefig(path_galaxy_plots + fout)
        print('Wrote: %s' % fout)
        
        # Plot interpolated maps of chisq, resid/galaxy, stellar kinematics
        qa.plot_multi_map(stkin_mapname, stkin_map_kwargs_interp)
        fout =  ('_').join([stem_file, 'stkin', 'maps', 'interp']) + '.png'
        plt.savefig(path_galaxy_plots + fout)
        print('Wrote: %s' % fout)
        
        # Plot binned maps of emission line fluxes
        qa.plot_multi_map(emflux_mapname, emflux_map_kwargs)
        fout =  ('_').join([stem_file, 'emflux', 'maps']) + '.png'
        plt.savefig(path_galaxy_plots + fout)
        print('Wrote: %s' % fout)
        
        # Plot binned maps of signal to noise and emission line kinematics
        qa.plot_multi_map(emkin_mapname, emkin_map_kwargs)
        fout =  ('_').join([stem_file, 'emkin', 'maps']) + '.png'
        plt.savefig(path_galaxy_plots + fout)
        print('Wrote: %s' % fout)
    
    
    #overplot_all = False
    if overplot_all:
    
        # Overplot all residuals (observed frame)
        qa.plot_all_spec(**all_resid_obs_kwargs)
        fout =  ('_').join([stem_file, 'all', 'resid', 'obs']) + '.png'
        plt.savefig(path_galaxy_plots + fout, dpi=200)
        print('Wrote: %s' % fout)
        
        # Overplot all residuals (rest frame)
        qa.plot_all_spec(**all_resid_rest_kwargs)
        fout =  ('_').join([stem_file, 'all', 'resid', 'rest']) + '.png'
        plt.savefig(path_galaxy_plots + fout, dpi=200)
        print('Wrote: %s' % fout)
        
        # Overplot all spectra (observed frame)
        qa.plot_all_spec(**all_spec_obs_kwargs)
        fout =  ('_').join([stem_file, 'all', 'spec', 'obs']) + '.png'
        plt.savefig(path_galaxy_plots + fout, dpi=200)
        print('Wrote: %s' % fout)
        
        # Overplot all spectra (rest frame)
        qa.plot_all_spec(**all_spec_rest_kwargs)
        fout =  ('_').join([stem_file, 'all', 'spec', 'rest']) + '.png'
        plt.savefig(path_galaxy_plots + fout, dpi=200)
        print('Wrote: %s' % fout)
        
    
    plot_spec_pdf = True #False
    if plot_spec_pdf:
    
        # Plot spectra individually (one file)
        if plot_all_spec_as_pdf:
            spec_to_plot = np.arange(qa.n_bins)
        else:
            spec_to_plot = np.arange(0, qa.n_bins, 10)
        
        fout =  ('_').join([stem_file, 'resid', 'all', 'bins']) + '.pdf'
        print('Writing: %s ...' % fout)
        with PdfPages(path_galaxy_plots + fout) as pdf:
            for bin in spec_to_plot:
                fig = qa.plot_resid(bin=bin, **resid_kwargs)
                try:
                    pdf.savefig(fig)
                except TypeError as e:
                    print('TypeError')
                plt.close()
                plt.clf()
                if bin % 25 == 0:
                    print('    bin %s...done' % bin)
        
        plt.close()
        print('Wrote: %s' % fout)
    
    open(path_galaxy_plots + done_file, 'a').close()
    
    print('')
    #--------------

#---------------------------------------------



