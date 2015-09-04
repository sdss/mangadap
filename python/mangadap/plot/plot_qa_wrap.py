#!/usr/bin/env python
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

from __future__ import division, print_function, absolute_import

import os
from os.path import join
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

from plot_qa import PlotQA
from plot_qa import FitError


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
        file_list = 'CUBE_files_to_plot.txt'
except IndexError:
    file_list = 'CUBE_files_to_plot.txt'


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
path_gal = join(path_gal, '')
print('Path: {}\n'.format(path_gal))
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
files = np.genfromtxt(join(path_gal, file_list), dtype='str')
files = np.atleast_1d(files)
#-------------------------------

# central wavelength for the emission line zoom in plots
win_cen = np.array([3727., 4861., 4985., 6565., 6565., 6723.])
win_names = ['oii', 'hbeta', 'oiii', 'nii', 'halpha', 'sii']
#-----

#----- Make plots for each DAP output file ------

for dap_file in files:
    
    #----- Read in galaxy ID and DAP parameters -----

    stem_file = dap_file.strip('.fits')
    done_file = ''.join((stem_file, '.done'))
    
    ig1, plate, ifudesign, mode_in, binning_type, niter = stem_file.split('-')
    manga_pid = '-'.join([plate, ifudesign])
    mode = mode_in.split('_')[0].strip('LOG')

    fin_kws = dict(plate=plate, ifudesign=ifudesign, mode=mode,
                   bintype=binning_type, niter=niter)
    
    path_gal_plots = join(path_gal, 'plots', '')
    path_gal_plots_spec = join(path_gal_plots, 'spectra')
    path_gal_plots_maps = join(path_gal_plots, 'maps')
    path_gal_plots_gradients = join(path_gal_plots, 'gradients')
    
    # create plot directories if necessary
    try:
        os.makedirs(path_gal_plots_spec)
        print('\nCreated directory: {}\n'.format(path_gal_plots_spec))
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise e
        pass
    
    try:
        os.makedirs(path_gal_plots_maps)
        print('\nCreated directory: {}\n'.format(path_gal_plots_maps))
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise e
        pass
    
    try:
        os.makedirs(path_gal_plots_gradients)
        print('\nCreated directory: {}\n'.format(path_gal_plots_gradients))
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise e
        pass
    
    # check for done file
    if not overwrite:
        if os.path.isfile(join(path_gal_plots, done_file)):
            print('{:12} {:6}: '.format(manga_pid, binning_type) + 
                  'Plots already exist. Use overwrite keyword to remake them')
            continue
    
    #------------------------------------------------- 
    
    
    #----- Plots for each binning type  ----- 
    if mode == 'CUBE':
        if binning_type == 'NONE':
            plot_map = True
            plot_indiv_map = True
            plot_gradients = False
            plot_indiv_gradients = False
            overplot_all = True
            plot_all_spec_as_pdf = False
            plot_all_spec_as_png = False
            plot_h3_h4 = False
        elif binning_type == 'STON':
            plot_map = True
            plot_indiv_map = True
            plot_gradients = False
            plot_indiv_gradients = False
            overplot_all = True
            plot_all_spec_as_pdf = True
            plot_all_spec_as_png = True
            plot_h3_h4 = True
        elif binning_type == 'RADIAL':
            plot_map = False
            plot_indiv_map = False
            plot_gradients = True
            plot_indiv_gradients = True
            overplot_all = True
            plot_all_spec_as_pdf = True
            plot_all_spec_as_png = True
            plot_h3_h4 = False
        elif binning_type == 'ALL':
            plot_map = False
            plot_indiv_map = False
            plot_gradients = False
            plot_indiv_gradients = False
            overplot_all = False
            plot_all_spec_as_pdf = True
            plot_all_spec_as_png = True
            plot_h3_h4 = False
    elif mode == 'RSS':
        if binning_type == 'RADIAL':
            plot_map = False
            plot_indiv_map = False
            plot_gradients = True
            plot_indiv_gradients = True
            overplot_all = True
            plot_all_spec_as_pdf = True
            plot_all_spec_as_png = True
            plot_h3_h4 = False
        if binning_type == 'ALL':
            plot_map = False
            plot_indiv_map = False
            plot_gradients = False
            plot_indiv_gradients = False
            overplot_all = False
            plot_all_spec_as_pdf = True
            plot_all_spec_as_png = True
            plot_h3_h4 = False
    #----------------------------------------
    
    
    #----- Calculate Fit Metrics -----
    # Remove
    #reload(plot_qa)
    #from plot_qa import PlotQA
    try:
        qa = PlotQA(join(path_gal, dap_file))
        
        # UNCOMMENT
        # qa = PlotQA(fin_kws)
        print('pid-ifu: {}, binning: {}'.format(manga_pid, binning_type))
        print('Template Library:', qa.tpl_lib)
        qa.select_wave_range()
        qa.set_axis_lims()
        qa.calc_chisq()
        qa.calc_resid()
        qa.calc_resid_data()
        #qa.ch1, qa.ch1_r = qa.define_custom_cubehelix(rot=-2, start=1, gamma=1)
        #qa.ch4, qa.ch4_r = qa.define_custom_cubehelix(rot=1.5, start=0.5, gamma=1)
    except FitError:
        print('\n{:12} {:6}: *NOT* fit by DAP\n'.format(manga_pid,
              binning_type))
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
    # qa.emvdisp is 0 or 99 for all spaxels
    # ind_emvdisp = np.where(qa.emvdisp_ew / qa.emvdisperr_ew >= 3. & (qa.stvdisperr > 0.))
    vdisp_range = sigma_clip(np.concatenate((qa.stvdisp[ind_stvdisp],
                             qa.emvdisp_ew)), sig=3)
    cbrange_vdisp = [vdisp_range.min(), vdisp_range.max()]
    
    stvel_args = dict(val=qa.stvel_rest, val_err=qa.stvelerr_rest,
                      kwargs=dict(cblabel=r'v$_\star$ [km/s]',
                                  cmap=cm.coolwarm,
                                  cbrange_symmetric=True,
                                  title_text=r'v_star',
                                  nodots=True))
    stvdisp_args = dict(val=qa.stvdisp, val_err=qa.stvdisperr,
                        kwargs=dict(cblabel=r'$\sigma_\star$ [km/s]',
                                    #cbrange=cbrange_vdisp,
                                    #cbrange_clip=False,
                                    cbrange_clip=True,
                                    snr_thresh=3.,
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
                                  title_text=r'v_gas (Wang)',
                                  nodots=True))
    emvdisp_args = dict(val=qa.emvdisp_ew,
                     kwargs=dict(cblabel=r'$\sigma_{\rm gas}$ [km/s]',
                                 #cbrange=cbrange_vdisp,
                                 #cbrange_clip=False,
                                 cbrange_clip=True,
                                 cmap=qa.linearL,
                                 title_text=r'$\sigma$_gas (Wang)',
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
    for v in kin_map_kwargs_interp.values():
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
                        kwargs=dict(title_text=r'[OII] $\lambda$3727'))
    hbeta_args = dict(val=qa.hbeta_ew, val2=qa.hbeta_fb,
                      val_err=qa.hbetaerr_ew, val2_err=qa.hbetaerr_fb,
                      kwargs=dict(title_text=r'H$\beta$'))
    oiii5007_args = dict(val=qa.oiii5007_ew, val2=qa.oiii5007_fb,
                         val_err=qa.oiii5007err_ew, val2_err=qa.oiii5007err_fb,
                         kwargs=dict(title_text=r'[OIII] $\lambda$5007'))
    halpha_args = dict(val=qa.halpha_ew, val2=qa.halpha_fb,
                       val_err=qa.halphaerr_ew, val2_err=qa.halphaerr_fb,
                       kwargs=dict(title_text=r'H$\alpha$'))
    nii6583_args = dict(val=qa.nii6583_ew, val2=qa.nii6583_fb,
                        val_err=qa.nii6583err_ew, val2_err=qa.nii6583err_fb,
                        kwargs=dict(title_text=r'[NII] $\lambda$6583'))
    sii6717_args = dict(val=qa.sii6717_ew, val2=qa.sii6717_fb,
                        val_err=qa.sii6717err_ew, val2_err=qa.sii6717err_fb,
                        kwargs=dict(title_text=r'[SII] $\lambda$6717'))
    
    
    emflux_map_kwargs = dict(oii3727=oii3727_args,
                             hbeta=hbeta_args,
                             oiii5007=oiii5007_args,
                             halpha=halpha_args,
                             nii6583=nii6583_args,
                             sii6717=sii6717_args)
    
    for v in emflux_map_kwargs.values():
        v['kwargs']['cblabel'] = r'Flux [10$^{-17}$ erg/s/cm$^2$]'
        v['kwargs']['cmap'] = qa.linearL
        v['kwargs']['nodots'] = True
        v['kwargs']['val_no_measure'] = 0.
        v['kwargs']['cbrange_clip'] = False
        v['kwargs']['show_flux_contours'] = False
    
    emflux_ew_map_kwargs = copy.deepcopy(emflux_map_kwargs)
    emflux_fb_map_kwargs = copy.deepcopy(emflux_map_kwargs)
    
    for v in emflux_fb_map_kwargs.values():
        v['val'] = v['val2']
        v['val_err'] = v['val2_err']

    for dd, nm in zip([emflux_ew_map_kwargs, emflux_fb_map_kwargs],
                      ['Wang', 'Belfiore']):
        for vv in dd.values():
            vv['kwargs']['title_text'] += ' (%s)' % nm
            del vv['val2']
            del vv['val2_err']
    
    #---------------------------------
    
    #--- Spectral Index Map Parameters ---



    spec_ind_mapname = [
      'D4000',
      'HDeltaA',
      'Hb',
      'Mgb'
      'Fe5270_5335'
      'CaII0p86']
    # need name of spectal index (not a combination of indices) to get units
    spec_ind_name = [
      'D4000',
      'HDeltaA',
      'Hb',
      'Mgb',
      'Fe5270',
      'CaII0p86A']
    D4000_args = dict(val=qa.D4000, val_err=qa.D4000err,
                      kwargs=dict(title_text=r'D4000'))
    HDeltaA_args = dict(val=qa.hbeta_ew, val_err=qa.hbetaerr_ew,
                        kwargs=dict(title_text=r'H$\delta$A'))
    Hb_args = dict(val=qa.Hb, val_err=qa.Hberr,
                   kwargs=dict(title_text=r'H$\beta$'))
    Mgb_args = dict(val=qa.Mgb, val_err=qa.Mgberr,
                    kwargs=dict(title_text=r'Mg $b$'))
    Fe5270_5335_args = dict(val=qa.Fe5270_5335, val_err=qa.Fe5270_5335err,
                            kwargs=dict(title_text=r'0.72$\cdot$Fe5270 + ' +
                                        r'0.28$\cdot$Fe5335'))
    CaII0p86_args = dict(val=qa.CaII0p86, val_err=qa.CaII0p86,
                         kwargs=dict(title_text=r'(CaII0p86A + CaII0p86B + ' + 
                                     'CaII0p86C)/3.'))
    
    spec_ind_map_kwargs = dict(D4000=D4000_args,
                               HDeltaA=HDeltaA_args,
                               Hb=Hb_args,
                               Mgb=Mgb_args,
                               Fe5270_5335=Fe5270_5335_args,
                               CaII0p86=CaII0p86_args)
    
    for name, v in zip(spec_ind_name, spec_ind_map_kwargs.values()):
        spec_ind_units = qa.get_spec_ind_units(name, qa.sipar.SINAME, qa.sipar.UNIT)
        if spec_ind_units == 'ang':
            si_cblabel = r'$\AA$'
        elif spec_ind_units == 'mag':
            si_cblabel = 'Mag'
        else:
            raise('Unknown spectral index units.')
            si_cblabel = None
        v['kwargs']['cblabel'] = si_cblabel
        v['kwargs']['cmap'] = qa.linearL
        v['kwargs']['cbrange_symmetric'] = True
        v['kwargs']['nodots'] = True
        v['kwargs']['val_no_measure'] = 0.
        v['kwargs']['cbrange_clip'] = True
        v['kwargs']['show_flux_contours'] = False
    
        
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
    #
    #qa = PlotQA(path_gal + dap_file)
    #print('\nTemplate Library:', qa.tpl_lib)
    #qa.select_wave_range()
    #qa.set_axis_lims()
    #qa.calc_chisq()
    #qa.calc_resid()
    #qa.calc_resid_data()
    #qa.plot_multi_map(emflux_mapname, emflux_map_kwargs)
    
    # plot_map = False
    if plot_map:
    
        # Plot bin numbers on top of chisq map
        qa.plot_map(qa.chisq_bin, **bin_num_map_kwargs)
        fout =  ('_').join([stem_file, 'bin', 'num']) + '.png'
        plt.savefig(path_gal_plots + fout)
        print('Wrote: %s' % fout)
        fout =  ('_').join([stem_file, 'bin', 'num']) + '.pdf'
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
        
        # Plot binned maps of emission line fluxes (Enci Wang's code)
        qa.plot_multi_map(emflux_mapname, emflux_ew_map_kwargs)
        fout =  ('_').join([stem_file, 'emflux', 'ew','maps']) + '.png'
        plt.savefig(path_gal_plots + fout)
        print('Wrote: %s' % fout)
        
        # Plot binned maps of emission line fluxes (Francesco Belfiore's code)
        qa.plot_multi_map(emflux_mapname, emflux_fb_map_kwargs)
        fout =  ('_').join([stem_file, 'emflux', 'fb', 'maps']) + '.png'
        plt.savefig(path_gal_plots + fout)
        print('Wrote: %s' % fout)
        
        # Plot binned maps of signal to noise and emission line kinematics
        qa.plot_multi_map(snr_mapname, snr_map_kwargs)
        fout =  ('_').join([stem_file, 'snr', 'maps']) + '.png'
        plt.savefig(path_gal_plots + fout)
        print('Wrote: %s' % fout)
    
    
    # plot_indiv_map = False
    if plot_indiv_map:
            
        for k, v in kin_map_kwargs.items():
            val_err = None
            if 'val_err' in v:
                val_err = v['val_err']
            qa.plot_map(v['val'], z_err=val_err, **v['kwargs'])
            fout =  ('_').join([stem_file, k, 'map']) + '.png'
            plt.savefig(path_gal_plots_maps + fout)
        print('Wrote: individual kin maps')
        
        if stkin_interp:
            for k, v in kin_map_kwargs_interp.items():
                val_err = None
                if 'val_err' in v:
                    val_err = v['val_err']
                qa.plot_map(v['val'], z_err=val_err, **v['kwargs'])
                fout =  ('_').join([stem_file, k, 'map', 'interp']) + '.png'
                plt.savefig(path_gal_plots_maps + fout)
            print('Wrote: individual kin interp maps')
        
        for k, v in emflux_ew_map_kwargs.items():
            val_err = None
            if 'val_err' in v:
                val_err = v['val_err']
            qa.plot_map(v['val'], z_err=val_err, **v['kwargs'])
            fout =  ('_').join([stem_file, k, 'ew', 'map']) + '.png'
            plt.savefig(path_gal_plots_maps + fout)
        print('Wrote: individual emflux_ew maps')
        
        for k, v in emflux_fb_map_kwargs.items():
            val_err = None
            if 'val_err' in v:
                val_err = v['val_err']
            qa.plot_map(v['val'], z_err=val_err, **v['kwargs'])
            fout =  ('_').join([stem_file, k, 'fb', 'map']) + '.png'
            plt.savefig(path_gal_plots_maps + fout)
        print('Wrote: individual emflux_fb maps')
        
        #--------
        #
        # Spectral Index Plots
        #
        #--------
        for k, v in spec_ind_map_kwargs.items():
            val_err = None
            if 'val_err' in v:
                val_err = v['val_err']
            qa.plot_map(v['val'], z_err=val_err, **v['kwargs'])
            fout =  ('_').join([stem_file, k, 'map']) + '.png'
            plt.savefig(join(path_gal_plots_maps, fout))
        print('Wrote: individual spectral index maps')

        for k, v in snr_map_kwargs.items():
            if k is 'halpha':
                pass
            else:
                val_err = None
                if 'val_err' in v:
                    val_err = v['val_err']
                qa.plot_map(v['val'], z_err=val_err, **v['kwargs'])
                fout =  ('_').join([stem_file, k, 'map']) + '.png'
                plt.savefig(path_gal_plots_maps + fout)
        print('Wrote: individual snr maps')
    
    
    plot_gradients = False
    if plot_gradients:
    
        # Plot emission line flux gradients
        qa.plot_multi_radial_gradients(emflux_mapname, emflux_map_kwargs)
        fout =  ('_').join([stem_file, 'emflux', 'gradients']) + '.png'
        plt.savefig(path_gal_plots + fout)
        print('Wrote: %s' % fout)
    
    if plot_indiv_gradients:
    
        for k, v in emflux_map_kwargs.items():
            qa.plot_radial_gradients(k, emflux_map_kwargs, c_ind=[2, 0],
                                     leglabels=['F. Belfiore', 'E. Wang'])
            fout =  ('_').join([stem_file, k, 'gradient']) + '.png'
            plt.savefig(path_gal_plots_gradients + fout)
        print('Wrote: individual emflux gradients')
    
    
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
        
    
    plot_spec_pdf = False
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
    
    
    plot_spec_png = False
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


    plot_emline_windows = False
    if plot_emline_windows:

        # Plot spectra in individually files
        if plot_all_spec_as_png:
            spec_to_plot = np.arange(qa.n_bins)
        else:
            if qa.n_bins < 100.:
                spec_to_plot = np.arange(qa.n_bins)
            else:
                spec_to_plot = np.linspace(0, qa.n_bins-1, 100).astype(int)

        for bin in spec_to_plot:
            # multi-panel plot
            fout =  ('_').join([stem_file, 'spec']) + '-%0.4i_%s.png' % (bin, 'emlines')
            if not os.path.isfile(fout):
                qa.plot_emline_multi(bin=bin)
                plt.savefig(path_gal_plots_spec + fout, dpi=200)
                plt.close()
                print('Wrote: %s' % fout)

            # individual plots
            for ii, (wc, wn) in enumerate(zip(win_cen, win_names)):
                nii = False
                if ii == 3:
                    nii = True
                fout =  ('_').join([stem_file, 'spec']) + '-%0.4i_%s.png' % (bin, wn)
                if not os.path.isfile(fout):
                    qa.plot_emline(bin=bin, nii=nii, win_cen=wc)
                    plt.savefig(path_gal_plots_spec + fout, dpi=200)
                    plt.close()
                    #print('Wrote: %s' % fout)
    
    open(path_gal_plots + done_file, 'a').close()
    
    print('')
    #--------------

#---------------------------------------------



# def airtovac(wave_air):
#     """Convert from air to vacuum wavelengths.
# 
#     Based on Ciddor (1996) (implemented in IDL AIRTOVAC)
# 
#     Args:
#         wave_vac (array): wavelengths in air
# 
#     Returns:
#         wavelengths in vacuum
#     """
#     sigma2 = (1e4 / wave_air)**2. # Convert to wavenumber squared
#     # Compute conversion factor
#     factor = (1. +  5.792105e-2/(238.0185 - sigma2) + 
#               (1.67917e-3 / ( 57.362 - sigma2)))
#     return wave_air * factor
# 
# 
# def set_emline_param(self):
#     keys = ['oii3727', 'hbeta', 'oiii4959', 'oiii5007', 'nii6548',
#         'halpha', 'nii6583', 'sii6717', ' sii6731']
#     names = ['[OII]3727', r'H$\beta$', '[OIII]4959', '[OIII]5007',
#         '[NII]6548', r'H$\alpha$', '[NII]6583', '[SII]6716', '[SII]6731']
#     wave_air = np.array([3727., 4861.32, 4958.83, 5006.77, 6547.96, 6562.80, 6583.34,
#         6716.31, 6730.68])
#     wave_vac = airtovac(wave_air)
#     
#     a = [names, wave_vac]
#     b = [list(x) for x in zip(*a)]
#     self.emline_par = pd.DataFrame(b, index=keys, columns=['name', 'wave'])
# 
# set_emline_param(qa)




#for w in win_cen:
#    plot_emline(qa, win_cen=w)
#qa.plot_emline(qa, bin=0, nii=True, win_cen=win_cen[3])
#plot_emline_multi(qa, bin=0)
