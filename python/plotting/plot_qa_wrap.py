#!/usr/bin/env python -W ignore::DeprecationWarning

"""
FILE
    plot_qa_wrap.py

DESCRIPTION

    Generate QA plots.
   
    Usage: python plot_qa_wrap.py path [file_list] [-overwrite]
        [-plot_spec_pdf] [-no_plot_spec_png] [-stkin_interp]

    ARGUMENTS:
        path
            Path with DAP output files.

        [file list] 
            File with the list of DAP files; default is qa_file_list.txt.

        [-overwrite]
            Overwrite any existing output file; default is False.

        [-plot_spec_pdf]
            Plot multiple spectra as separate pages in a single pdf file;
            default is False.

        [-no_plot_spec_png]
            Do NOT plot spectra as separate png files; default is False.

        [-stkin_interp]
            Plot an interpolated version of the stellar kinematics maps;
            default is False.

REVISION HISTORY:
    24 Feb 2015 - Original implementation by B. Andrews
    28 Apr 2015: (KBW) Some documentation edits and use of glob; added mode


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

from astropy.stats import sigma_clip

import matplotlib as mpl
mpl.use('Agg')
print('backend:', mpl.get_backend())
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.backends.backend_pdf import PdfPages

# for Utah cluster
# from plotting.plot_qa import PlotQA
# from plotting.plot_qa import FitError

# for Brett
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


#----- Read in command line arguments -----
try:    
    path_gal = sys.argv[1]
except:
    raise Exception('Usage: python plot_qa_wrap.py path [file list] [-overwrite]' + 
                    '[-no_plot_spec_png] [-plot_spec_pdf] [-no_stkin_interp]')

 
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
if '-plot_spec_pdf' in sys.argv:
    plot_spec_pdf = True
else:
    plot_spec_pdf = False
if '-no_plot_spec_png' in sys.argv:
    plot_spec_png = False
else:
    plot_spec_png = True
if '-stkin_interp' in sys.argv:
    stkin_interp = True
else:
    stkin_interp = False


print()
print('Overwrite plots:', overwrite)
print()
#----------------------------------------


#----- Set Path ------
if path_gal[-1] is not '/':
    path_gal += '/'

print('Path: %s\n' % path_gal)
#---------------------


# #----- Set Path and File Names -----
# manga_drp_ver = os.getenv('MANGADRP_VER')
# manga_dap_ver = os.getenv('MANGADAP_VER')
# path_analysis = os.getenv('MANGA_SPECTRO_ANALYSIS')
# 
# if test_dir:
#     if path_dap is None:
#         path_dap = '/'.join([path_analysis, manga_dap_ver, ''])
#     path_analysis_plots = os.getenv('MANGA_SPECTRO_ANALYSIS_PLOTS')
#     path_dap_plots = '/'.join([path_analysis_plots, manga_dap_ver, ''])
#     if path_file_list is None:
#         path_file_list = path_dap_plots
# else:
#     if path_dap is None:
#         path_dap = '/'.join([path_analysis, manga_drp_ver, manga_dap_ver, ''])
#     if path_file_list is None:
#         path_file_list = path_dap
# 
# 
# print('Input file list directory: %s' % (path_file_list))
# print()
# print('Input file list: %s' % (file_list))
# print()
# #-------------------------------


#----- Read the file names -----
files = np.genfromtxt(path_gal + file_list, dtype='str')
files = np.atleast_1d(files)
#-------------------------------



#----- Make plots for each DAP output file ------

for dap_file in files:
    
    #----- Read in galaxy ID and DAP parameters -----
    
    stem_file = dap_file.strip('.fits')
    done_file = stem_file + '.done'
    
    ig1, pid, ifudesign, mode_in, binning_type, exec_num = stem_file.split('-')
    
    manga_pid = '-'.join([pid, ifudesign])
    mode = mode_in.split('_')[0]
    
    path_gal_plots = path_gal + 'plots/'
    path_gal_plots_spec = path_gal_plots + 'spectra/'
    
    # create plot directories if necessary
    try:
        os.makedirs(path_gal_plots_spec)
        print('\nCreated directory: %s\n' % path_gal_plots_spec)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise e
        pass
    
    # check for done file
    if not overwrite:
        if os.path.isfile(path_gal_plots + done_file):
            print('%-12s %-6s: Plots already exist.' % (manga_pid, binning_type) + 
                  '  Use overwrite keyword to remake them')
            continue
    
    #------------------------------------------------- 
    
    
    #----- Plots for each binning type  ----- 
    if mode == 'LOGCUBE':
        if binning_type == 'NONE':
            plot_map = True
            plot_gradients = False
            overplot_all = True
            plot_all_spec_as_pdf = False
            plot_all_spec_as_png = False
            plot_h3_h4 = False
        elif binning_type == 'STON':
            plot_map = True
            plot_gradients = False
            overplot_all = True
            plot_all_spec_as_pdf = True
            plot_all_spec_as_png = True
            plot_h3_h4 = True
        elif binning_type == 'RADIAL':
            plot_map = False
            plot_gradients = True
            overplot_all = True
            plot_all_spec_as_pdf = True
            plot_all_spec_as_png = True
            plot_h3_h4 = False
        elif binning_type == 'ALL':
            plot_map = False
            plot_gradients = False
            overplot_all = False
            plot_all_spec_as_pdf = True
            plot_all_spec_as_png = True
            plot_h3_h4 = False
    elif mode == 'LOGRSS':
        # add this line
        # if binning_type == 'RADIAL':
        plot_map = False
        plot_gradients = True
        overplot_all = True
        plot_all_spec_as_pdf = True
        plot_all_spec_as_png = True
        plot_h3_h4 = False
    #----------------------------------------
    
    
    #----- Calculate Fit Metrics -----
    # Remove
    #reload(plot_qa)
    #from plot_qa import PlotQA
    try:
        qa = PlotQA(path_gal + dap_file)
        print('pid-ifu: %s, binning: %s' % (manga_pid, binning_type))
        print('Template Library:', qa.tpl_lib)
        qa.select_wave_range()
        qa.set_axis_lims()
        qa.calc_chisq()
        qa.calc_resid()
        qa.calc_resid_data()
        #qa.ch1, qa.ch1_r = qa.define_custom_cubehelix(rot=-2, start=1, gamma=1)
        #qa.ch4, qa.ch4_r = qa.define_custom_cubehelix(rot=1.5, start=0.5, gamma=1)
    except FitError:
        print('\n%-12s %-6s: *NOT* fit by DAP\n' % (manga_pid, binning_type))
        continue
    
    #---------------------------------
    
    
    
    #----- Plot Parameters -----
     
    
    #--- Kinematics Map Parameters ---
    
    if plot_h3_h4:
        kin_mapname = [
          'stvel',
          'stvdisp',
          'sth3',
          'emvel',
          'emvdisp',
          'sth4',]
    else:
        kin_mapname = [
          'stvel',
          'stvdisp',
          'chisq',
          'emvel',
          'emvdisp',
          'resid',]
    
    # plot sigma_star and sigma_gas on same scale
    ind_stvdisp = np.where((qa.stvdisp / qa.stvdisperr >= 3.) & (qa.stvdisperr > 0.))
    # qa.emvdisp is 0 for all spaxels
    # ind_emvdisp = np.where(qa.emvdisp_ew / qa.emvdisperr_ew >= 3. & (qa.stvdisperr > 0.))
    vdisp_range = sigma_clip(np.concatenate((qa.stvdisp[ind_stvdisp],
                             qa.emvdisp_ew)), sig=3)
    cbrange_vdisp = [vdisp_range.min(), vdisp_range.max()]
    
    stvel_args = dict(val=qa.stvel_rest, val_err=qa.stvelerr_rest,
                      kwargs=dict(cblabel=r'$v_\star$ [km/s]',
                                  cmap=cm.coolwarm,
                                  cbrange_symmetric=True,
                                  title_text=r'v_star',
                                  nodots=True))
    stvdisp_args = dict(val=qa.stvdisp, val_err=qa.stvdisperr,
                        kwargs=dict(cblabel=r'$\sigma_\star$ [km/s]',
                                    cbrange=cbrange_vdisp,
                                    cbrange_clip=False,
                                    cmap=qa.linearL, 
                                    title_text=r'$\sigma$_star',
                                  nodots=True))
    sth3_args = dict(val=qa.sth3, val_err=qa.sth3err,
                     kwargs=dict(cblabel='h3',
                                 cbrange_symmetric=True,
                                 cmap=cm.coolwarm,
                                 title_text='h3',
                                  nodots=True))
    sth4_args = dict(val=qa.sth4, val_err=qa.sth4err,
                     kwargs=dict(cblabel='h4',
                                 cbrange_symmetric=False,
                                 cmap=qa.linearL,
                                 title_text='h4',
                                  nodots=True))
    
    emvel_args = dict(val=qa.emvel_rest_ew, val_err=qa.emvelerr_rest_ew,
                      kwargs=dict(cblabel=r'v$_{\rm gas}$ [km/s]',
                                  cmap=cm.coolwarm,
                                  cbrange_symmetric=True,
                                  title_text=r'v_gas (Enci)',
                                  nodots=True))
    emvdisp_args = dict(val=qa.emvdisp_ew, val_err=qa.emvdisperr_ew,
                     kwargs=dict(cblabel=r'$\sigma_{\rm gas}$ [km/s]',
                                 cbrange=cbrange_vdisp,
                                 cbrange_clip=False,
                                 cmap=qa.linearL,
                                 title_text=r'$\sigma$_gas (Enci)',
                                 nodots=True))
    
    chisq_args = dict(val=qa.chisq_bin,
                      kwargs=dict(cblabel=r'$\chi_{\rm red}^2$',
                                  cmap=qa.linearL,
                                  title_text=r'$\chi_{\rm red}^2$',
                                  nodots=True))
    resid_args = dict(val=qa.resid_data_bin_percent99,
                      kwargs=dict(cblabel='99th percentile |resid| / galaxy',
                                  cmap=qa.linearL,
                                  title_text='99th percentile |resid| / galaxy',
                                  nodots=True))
    
    
    kin_map_kwargs = dict(stvel=stvel_args,
                          stvdisp=stvdisp_args,
                          sth3=sth3_args,
                          sth4=sth4_args,
                          emvel=emvel_args,
                          emvdisp=emvdisp_args,
                          chisq=chisq_args,
                          resid=resid_args,)
    
    kin_map_kwargs_interp = copy.deepcopy(kin_map_kwargs)
    for v in itervalues(kin_map_kwargs_interp):
        v['kwargs']['interpolated'] = True
    
    bin_num_map_kwargs = dict(cblabel=r'$\chi_{\rm red}^2$',
                              cmap=qa.linearL,
                              title_text=r'$\chi_{\rm red}^2$',
                              nodots=True,
                              spaxel_num=True,
                              figsize=(15, 12))
    
    #----------------------------------------
    
    
    #--- Emission Line Fluxes Map Parameters ---
    emflux_mapname = [
      'oii3727',
      'hbeta',
      'oiii5007',
      'halpha',
      'nii6583',
      'sii6717',]
    oii3727_args = dict(val=qa.oii3727_ew, val2=qa.oii3727_fb,
                        val_err=qa.oii3727err_ew, val2_err=qa.oii3727err_fb,
                        kwargs=dict(title_text=r'[OII] $\lambda$3727 (Enci)'))
    hbeta_args = dict(val=qa.hbeta_ew, val2=qa.hbeta_fb,
                      val_err=qa.hbetaerr_ew, val2_err=qa.hbetaerr_fb,
                      kwargs=dict(title_text=r'H$\beta$ (Enci)'))
    oiii5007_args = dict(val=qa.oiii5007_ew, val2=qa.oiii5007_fb,
                         val_err=qa.oiii5007err_ew, val2_err=qa.oiii5007err_fb,
                         kwargs=dict(title_text=r'[OIII] $\lambda$5007 (Enci)'))
    halpha_args = dict(val=qa.halpha_ew, val2=qa.halpha_fb,
                       val_err=qa.halphaerr_ew, val2_err=qa.halphaerr_fb,
                       kwargs=dict(title_text=r'H$\alpha$ (Enci)'))
    nii6583_args = dict(val=qa.nii6583_ew, val2=qa.nii6583_fb,
                        val_err=qa.nii6583err_ew, val2_err=qa.nii6583err_fb,
                        kwargs=dict(title_text=r'[NII] $\lambda$6583 (Enci)'))
    sii6717_args = dict(val=qa.sii6717_ew, val2=qa.sii6717_fb,
                        val_err=qa.sii6717err_ew, val2_err=qa.sii6717err_fb,
                        kwargs=dict(title_text=r'[SII] $\lambda$6717 (Enci)'))

    
    emflux_map_kwargs = dict(oii3727=oii3727_args,
                             hbeta=hbeta_args,
                             oiii5007=oiii5007_args,
                             halpha=halpha_args,
                             nii6583=nii6583_args,
                             sii6717=sii6717_args)
    
    for v in itervalues(emflux_map_kwargs):
        v['kwargs']['cblabel'] = r'Flux [10$^{-17}$ erg/s/cm$^2$]'
        v['kwargs']['cmap'] = qa.linearL
        v['kwargs']['nodots'] = True
        v['kwargs']['val_no_measure'] = 0.
        v['kwargs']['cbrange_clip'] = False
    
    #---------------------------------
    
    
    #--- SNR Map Parameters ---
    snr_mapname = [
      'signal',
      'noise',
      'snr',
      'halpha',
      'resid',
      'chisq',]
    
    
    signal_args = dict(val=qa.signal,
                      kwargs=dict(cblabel=r'signal',
                                  cmap=qa.linearL,
                                  title_text=r'DRPS signal',
                                  nodots=True))
    noise_args = dict(val=qa.noise,
                      kwargs=dict(cblabel=r'noise',
                                  cmap=qa.linearL,
                                  title_text='DRPS noise',
                                  nodots=True))
    snr_args = dict(val=qa.snr,
                        kwargs=dict(cblabel=r'S/N',
                                    cmap=qa.linearL, 
                                    title_text=r'DRPS S/N',
                                    nodots=True))
    
    
    snr_map_kwargs = dict(signal=signal_args,
                          noise=noise_args,
                          snr=snr_args,
                          chisq=chisq_args,
                          resid=resid_args,
                          halpha=halpha_args)
    #---------------------------------
    
    
    
    #--- Spectra Parameters ----
    
    lam_lim = [3600, 9300]
    
    # set the same min and max flux for each spectrum
    # p32 = np.array([np.percentile(qa.galaxy[i], 32) for i in range(qa.n_bins)])
    # p50 = np.array([np.percentile(qa.galaxy[i], 50) for i in range(qa.n_bins)])
    # p68 = np.array([np.percentile(qa.galaxy[i], 68) for i in range(qa.n_bins)])
    # sig1 = (np.sort(p68) - np.sort(p32))[-1] / 2.
    # fluxmax = np.max(p50) + 4. * sig1
    
    fluxmin = 0.
    fluxmax = None
    
    residmin = -0.25
    residmax = 0.25
    
    resid_kwargs = dict(xlim=lam_lim,
                        ylim=[[fluxmin, fluxmax], None], 
                        #ylim=[[fluxmin, fluxmax], [residmin, residmax]],
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
    
    #reload(plot_qa)
    #from plot_qa import PlotQA
     
    #qa = PlotQA(path_gal + dap_file)
    #print('\nTemplate Library:', qa.tpl_lib)
    #qa.select_wave_range()
    #qa.set_axis_lims()
    #qa.calc_chisq()
    #qa.calc_resid()
    #qa.calc_resid_data()
    #qa.plot_multi_map(emflux_mapname, emflux_map_kwargs)
    
    
    #plot_map = False
    if plot_map:
    
        # Plot bin numbers on top of chisq map
        qa.plot_map(qa.chisq_bin, **bin_num_map_kwargs)
        fout =  ('_').join([stem_file, 'bin', 'num']) + '.png'
        plt.savefig(path_gal_plots + fout)
        print('Wrote: %s' % fout)
        
        # Plot binned maps of chisq, resid/galaxy, stellar kinematics
        qa.plot_multi_map(kin_mapname, kin_map_kwargs)
        fout =  ('_').join([stem_file, 'kin', 'maps']) + '.png'
        plt.savefig(path_gal_plots + fout)
        print('Wrote: %s' % fout)
        
        if stkin_interp:
            # Plot interpolated maps of chisq, resid/galaxy, stellar kinematics
            qa.plot_multi_map(kin_mapname, kin_map_kwargs_interp)
            fout =  ('_').join([stem_file, 'kin', 'maps', 'interp']) + '.png'
            plt.savefig(path_gal_plots + fout)
            print('Wrote: %s' % fout)
        
        # Plot binned maps of emission line fluxes
        qa.plot_multi_map(emflux_mapname, emflux_map_kwargs)
        fout =  ('_').join([stem_file, 'emflux', 'maps']) + '.png'
        plt.savefig(path_gal_plots + fout)
        print('Wrote: %s' % fout)
        
        # Plot binned maps of signal to noise and emission line kinematics
        qa.plot_multi_map(snr_mapname, snr_map_kwargs)
        fout =  ('_').join([stem_file, 'snr', 'maps']) + '.png'
        plt.savefig(path_gal_plots + fout)
        print('Wrote: %s' % fout)
    
    
    plot_gradients = False
    if plot_gradients:
    
        # Plot emission line flux gradients
        qa.plot_radial_gradients(emflux_mapname, emflux_map_kwargs)
        fout =  ('_').join([stem_file, 'emflux', 'gradients']) + '.png'
        plt.savefig(path_gal_plots + fout)
        print('Wrote: %s' % fout)
    
    
    overplot_all = False
    if overplot_all:
    
        # Overplot all residuals (observed frame)
        qa.plot_all_spec(**all_resid_obs_kwargs)
        fout =  ('_').join([stem_file, 'all', 'resid', 'obs']) + '.png'
        plt.savefig(path_gal_plots + fout, dpi=200)
        print('Wrote: %s' % fout)
        
        # Overplot all residuals (rest frame)
        qa.plot_all_spec(**all_resid_rest_kwargs)
        fout =  ('_').join([stem_file, 'all', 'resid', 'rest']) + '.png'
        plt.savefig(path_gal_plots + fout, dpi=200)
        print('Wrote: %s' % fout)
        
        # Overplot all spectra (observed frame)
        qa.plot_all_spec(**all_spec_obs_kwargs)
        fout =  ('_').join([stem_file, 'all', 'spec', 'obs']) + '.png'
        plt.savefig(path_gal_plots + fout, dpi=200)
        print('Wrote: %s' % fout)
        
        # Overplot all spectra (rest frame)
        qa.plot_all_spec(**all_spec_rest_kwargs)
        fout =  ('_').join([stem_file, 'all', 'spec', 'rest']) + '.png'
        plt.savefig(path_gal_plots + fout, dpi=200)
        print('Wrote: %s' % fout)
        
    
    #plot_spec_pdf = False
    if plot_spec_pdf:
    
        # Plot spectra individually (one file)
        if plot_all_spec_as_pdf:
            spec_to_plot = np.arange(qa.n_bins)
        else:
            spec_to_plot = np.arange(0, qa.n_bins, 50)
                
        fout =  ('_').join([stem_file, 'resid', 'all', 'bins']) + '.pdf'
        print('Writing: %s ...' % fout)
        with PdfPages(path_gal_plots + fout) as pdf:
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
    
    
    #plot_spec_png = False
    if plot_spec_png:
    
        # Plot spectra in individually files
        if plot_all_spec_as_png:
            spec_to_plot = np.arange(qa.n_bins)
        else:
            #spec_to_plot = np.arange(0, qa.n_bins, 50)
            if qa.n_bins < 100.:
                spec_to_plot = np.arange(qa.n_bins)
            else:
                spec_to_plot = np.linspace(0, qa.n_bins-1, 100).astype(int)
        
        fout =  ('_').join([stem_file, 'resid', 'all', 'bins']) + '.png'
        for bin in spec_to_plot:
            fig = qa.plot_resid(bin=bin, **resid_kwargs)
            fout =  ('_').join([stem_file, 'spec']) + '-%0.4i.png' % bin
            plt.savefig(path_gal_plots_spec + fout, dpi=200)
            plt.close()
            print('Wrote: %s' % fout)
    
    
    open(path_gal_plots + done_file, 'a').close()
    
    print('')
    #--------------

#---------------------------------------------





