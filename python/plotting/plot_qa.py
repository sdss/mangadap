"""
FILE
    plot_qa.py

DESCRIPTION

    Generate quality assessment plots.
   
"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import os
import numpy as np
import copy
import fitsio

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.patches as patches
from matplotlib.ticker import MaxNLocator
import cubehelix

import warnings
try:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        import seaborn as sns
except ImportError:
    seaborn_installed = False


class PlotQA:
    def __init__(self, filename):
        self.dap_file = filename
        self.manga_id = ('-').join(self.dap_file.split('-')[1:3])
        self.setup()

    def setup(self):
        self.read_fits()
        self.deredshift()
        self.plot_setup()

    def fields(self):
        ret = []
        for nm in dir(self):
           if not nm.startswith('__') and not callable(getattr(self, nm)):
              ret.append(nm)
        return ret

    def read_fits(self):
        fin = fitsio.FITS(self.dap_file)
        self.drps = fin['DRPS'].read()
        self.bin_prop = fin['BINS'].read()
        self.wave_obs = fin['WAVE'].read()
        self.galaxy = fin['FLUX'].read().T
        self.ivar = fin['IVAR'].read().T
        self.model = fin['SMOD'].read().T
        self.stfit = fin['STFIT'].read()
        self.smsk = fin['SMSK'].read().T
        self.stfit = fin['STFIT'].read().T
        self.n_bins, self.n_pix = self.galaxy.shape

        # bin properties
        self.xpos, self.ypos, ig1, ig2, ig3, ig4, self.binvr, self.binid, \
            self.binw = [self.drps[name] for name in self.drps.dtype.names]
        self.binxrl, self.binyru, self.binr, self.bina, self.binsn, \
            self.nbin, self.binf = \
            [self.bin_prop[name] for name in self.bin_prop.dtype.names]

        # stellar kinematics
        self.tplw, self.addpoly, self.multpoly, self.kin, self.kinerr, \
            self.rchi2 = [self.stfit[name] for name in self.stfit.dtype.names]
        self.stvel, self.stvdisp, self.sth3, self.sth4 = self.kin.T

        h = fitsio.read_header(self.dap_file, 0)
        self.tpl_lib = h['TPLKEY']


    def deredshift(self):
        c = 299792.458
        self.v_star = np.array([item[3][0] for item in self.stfit]) 
        self.z = self.v_star / c
        self.z_median = np.median(self.z)

        # de-redshift spectra
        wave_obs_grid = np.ones((self.n_bins, self.n_pix)) * self.wave_obs
        self.wave_rest = (wave_obs_grid.T / (1. + self.z)).T
        self.wave_rest_med = self.wave_obs / (1. + self.z_median)
        self.dlam = ((self.wave_rest_med[1] - self.wave_rest_med[0]) / 
                     self.wave_rest_med[0])

        # conserve flux by multiplying by (1+z)
        self.galaxy_rest = (self.galaxy.T * (1. + self.z)).T
        self.ivar_rest = (self.ivar.T * (1. + self.z)).T
        self.model_rest = (self.model.T * (1. + self.z)).T


    def select_wave_range(self, lam_good=None):
        '''Default wavelength range:

        includes:
        - blue-ward of [SII]6731,
        - CaT region
        
        excludes:
        - sky lines at 4230, 5575
        '''
        if lam_good is None:
            self.lam_good = np.array([[3650., 4200.], # ends at 4230 sky line
                                      [4250., 5400.], # ends at 5575 sky line
                                      [5450., 6732.]]) # ends at [SII]
                                      #[8470., 8770.]]) # CaT
        else:
            self.lam_good = lam_good

        n_regions = len(self.lam_good)
        ind = []
        for i in xrange(n_regions):
            ind.append(np.where(
                       (self.wave_rest > self.lam_good[i, 0]) & 
                       (self.wave_rest < self.lam_good[i, 1])))
        
        ind0 = np.concatenate([ind[i][0] for i in xrange(n_regions)])
        ind1 = np.concatenate([ind[i][1] for i in xrange(n_regions)])
        ind_lexsort = np.lexsort((ind1, ind0))
        self.ind_lam = (ind0[ind_lexsort], ind1[ind_lexsort])

    #-------------------------------------------------------------------------



    #----- Derive Quantities -------------------------------------------------

    def calc_metric_per_bin(self, metric_pix):
        '''
        calculate average metric (e.g., chisq or residuals) for each binned
        spectrum while masking bad pixels
        '''
        metric_bin = np.ones(self.n_bins) * np.nan
        for i in xrange(self.n_bins):
            indt = np.where(self.ind_lam[0] == i)[0]
            metric_bin[i] = np.nanmean(metric_pix[i][self.ind_lam[1][indt]])
        return metric_bin


    def calc_chisq(self):
        '''
        chisq_pix: chisq at every pixel
        chisq_bin: chisq for each binned spectrum
        '''
        self.chisq_pix = (self.galaxy - self.model)**2. * self.ivar
        self.chisq_pix[np.where(self.smsk==1)] = np.nan
        self.chisq_bin = self.calc_metric_per_bin(self.chisq_pix)


    def calc_resid(self):
        '''
        resid_pix: residual at every pixel
        resid_bin: residual for each binned spectrum
        '''
        self.resid_pix = self.galaxy - self.model
        self.resid_pix[np.where(self.smsk==1)] = np.nan
        self.resid_bin = self.calc_metric_per_bin(self.resid_pix)
        

    def calc_resid_mod(self):
        '''
        resid_mod_pix: residual / model at every pixel
        resid_mod_bin: residual / model for each binned spectrum
        '''
        self.resid_mod_pix = (self.galaxy - self.model) / self.model
        self.resid_mod_pix[np.where(self.smsk==1)] = np.nan
        self.resid_mod_bin = self.calc_metric_per_bin(self.resid_mod_pix)
        

    def calc_resid_err(self):
        '''
        resid_err_pix: residual / error at every pixel
        resid_err_bin: residual / error for each binned spectrum
        '''
        self.resid_err_pix = (self.galaxy - self.model) * np.sqrt(self.ivar)
        self.resid_err_pix[np.where(self.smsk==1)] = np.nan
        self.resid_err_bin = self.calc_metric_per_bin(self.resid_err_pix)

    def calc_resid_data(self):
        '''
        resid_data_pix: residual / data at every pixel
        resid_data_bin: residual / data for each binned spectrum
        resid_data_bin_percent99: 99th percentile of resid_data_bin
        resid_data_bin_percent68: 68th percentile of resid_data_bin
        '''
        self.resid_data_pix = np.abs(self.galaxy - self.model) / self.galaxy
        self.resid_data_pix[np.where(self.smsk==1)] = np.nan
        self.resid_data_bin = self.calc_metric_per_bin(self.resid_data_pix)
        self.resid_data_bin_percent99 = np.array(
            [np.percentile(self.resid_data_pix[i][np.where(self.smsk[i]==0)], 99.)
            for i in xrange(self.n_bins)])
        self.resid_data_bin_percent68 = np.array(
            [np.percentile(self.resid_data_pix[i][np.where(self.smsk[i]==0)], 68.)
            for i in xrange(self.n_bins)])

    #-------------------------------------------------------------------------
    


    #----- Plot Setup --------------------------------------------------------

    def define_custom_cubehelix(self, rot=1, start=1, gamma=1):
        self.cubehelix_cust = cubehelix.cmap(reverse=True, rot=rot,
                                             start=start, gamma=gamma)


    def plot_setup(self):
        self.cubehelix111 = cubehelix.cmap(reverse=True, rot=1, start=1, gamma=1)

    #-------------------------------------------------------------------------



    #----- Plots -------------------------------------------------------------
    

    def find_map_lims(x):
        return np.array([np.floor(np.min(x)), np.ceil(np.max(x))])


    def set_axis_lims(self):
        ind_b = np.where(self.binid >= 0)
        self.xy_lim = np.max((self.xpos[ind_b], self.ypos[ind_b])) * 1.2
        
        # Don't need to set these unless comparing different analyses of the
        # same galaxy

        #self.chisq_lim = self.find_map_lims(self.chisq_bin)
        #self.stvel_lim = 
        #self.stvdisp_lim = 
        #self.sth3_lim = 
        #self.sth4_lim = 
        

    def count_models(x):
        if x is not None:
            try:
                if len(x.shape) == 2:
                    x = [x]
                    n = 1
            except AttributeError:
                n = len(x)
        else:
            n = 0
        return x, n
    
    
    def generate_colors(colors=None, n_models=0):
        if colors is None:
            n_colors = 1
            if n_models > 0:
                n_colors = n_models
            ctmp = sns.color_palette('bright', n_colors)
            colors = ctmp
            colors[1] = ctmp[2]
            colors[2] = ctmp[1]
        return colors
    

    def plot_hist(self, x, label, axlim=None):
        '''
        Plot distribution of metric (e.g., chisq).
        '''
        sns.set_context('talk')
        fig, ax = plt.subplots(figsize=(8, 6))
        sns.axlabel(label, '$N$')
        ax.set_title(label)
        c = sns.color_palette()[0]
        sns.distplot(x, kde=False, color=c)
        sns.distplot(x, kde=False, color=c,
                     hist_kws=dict(histtype='step', alpha=1, lw=2))
        if axlim is not None:
            ax.axis(axlim)


    def plot_multi_map(self,
                       map_order,
                       args,
                       n_ax=6,
                       figsize=(20, 12)):
        '''
        Plot multiple maps at once.
        '''

        fig = plt.figure(figsize=figsize)
        sns.set_context('poster', rc={'lines.linewidth': 2})
        bigAxes = fig.add_axes([0.04, 0.05, 0.9, 0.88], frameon=False)
        bigAxes.set_xticks([])
        bigAxes.set_yticks([])
        bigAxes.set_xlabel('arcsec')
        bigAxes.set_ylabel('arcsec')
        bigAxes.set_title(self.manga_id)

        for i in xrange(n_ax):
            dx = 0.31 * i
            dy = 0.45
            if i >= (n_ax / 2):
                dx = 0.31 * (i - n_ax / 2)
                dy = 0
            left, bottom = (0.08+dx, 0.1+dy)
            sns.set_context('talk', rc={'lines.linewidth': 2})
            sns.set_style(rc={'axes.facecolor': '#A8A8A8'})
            ax = fig.add_axes([left, bottom, 0.2, 0.3333])
            d = args[map_order[i]]
            k = d['kwargs']
            self.plot_map(d['val'], fig=fig, ax=ax, axloc=[left, bottom],
                          **d['kwargs'])


    def plot_map(self,
                 z,
                 interpolated=False,
                 flux=None,
                 cblabel=None,
                 cbrange=None,
                 cbrange_symmetric=False,
                 n_ticks=7,
                 cmap=cm.RdBu_r,
                 xy_max=None,
                 title_text=None,
                 spaxel_num=False,
                 nodots=False,
                 fontsize=7,
                 n_colors=64,
                 fig=None,
                 ax=None,
                 axloc=None,
                 figsize=(10, 8)):
        '''
        Plot map
        '''
        if ax is None:
            sns.set_context('poster', rc={'lines.linewidth': 2})
        else:
            sns.set_context('talk', rc={'lines.linewidth': 2})
        
        sns.set_style(rc={'axes.facecolor': '#A8A8A8'})

        if ax is None:
            fig = plt.figure(figsize=figsize)
            ax = fig.add_axes([0.12, 0.1, 2/3., 5/6.])
            #ax = fig.add_axes([0.15, 0.15, 0.6, 0.75])

            sns.axlabel('arcscec', 'arcsec')

        if title_text is not None:
            ax.set_title(title_text)

        ax.grid(False, which='both', axis='both')

        ind_tmp = np.where(self.binid >= 0)[0]
        ind_b = self.binid[ind_tmp]
        bin_num = np.arange(len(z), dtype=int).astype(str)

        spaxel_size = 0.5 # arcsec

        if not nodots:
            ax.plot(self.binxrl, self.binyru, '.k', markersize=3, zorder=10)

        # colorbar range
        cbrange_tmp = [z.min(), z.max()]
        if cbrange is None:
            cbrange = cbrange_tmp
        else:
            for i in xrange(len(cbrange)):
                if cbrange[i] is None:
                    cbrange[i] = cbrange_tmp[i]
        if cbrange_symmetric:
            cb_max = np.max(np.abs(cbrange))
            cbrange = [-cb_max, cb_max]
        cbdelta = cbrange[1] - cbrange[0]
    
        # plot spaxels
        if interpolated:
            levels = np.linspace(cbrange[0], cbrange[1], n_colors)
            p = ax.tricontourf(self.binxrl, self.binyru, z, levels=levels,
                               cmap=cmap)
        else:
            for i in xrange(len(ind_tmp)):
                ctmp = cmap((z[ind_b][i] - cbrange[0]) / cbdelta)
                delta = spaxel_size / 2.
                rect = patches.Rectangle((self.xpos[ind_tmp][i] - delta,
                                         self.ypos[ind_tmp][i] - delta),
                                         width=spaxel_size, height=spaxel_size,
                                         color=ctmp)
                ax.add_patch(rect)
                # dummy plot for colorbar
                p = ax.scatter(self.xpos[ind_tmp][0], self.ypos[ind_tmp][0],
                               cmap=cmap, c=z[ind_b][0], s=0, vmin=cbrange[0],
                               vmax=cbrange[1])
        if spaxel_num:
            for i in xrange(self.n_bins):
                if (i >= 100) and (self.nbin[i] <= 2):
                    fontsize_tmp = fontsize - 1
                elif self.nbin[i] == 1:
                    fontsize_tmp = fontsize
                else:
                    fontsize_tmp = fontsize + 2
                ctmp = cmap((z[i] - cbrange[0]) / cbdelta)
                ax.text(self.binxrl[i] - (spaxel_size / 4),
                        self.binyru[i] - (spaxel_size / 4),
                        str(i), fontsize=fontsize_tmp,
                        color=(1.-ctmp[0], 1.-ctmp[1], 1.-ctmp[2], ctmp[3]))

        # overplot flux contours (which flux? in what band?)
        if flux is not None:
            ax.tricontour(self.binxrl, self.binyru,
                          -2.5*np.log10(flux/np.max(flux).ravel()),
                          levels=np.arange(20), colors='k') # 1 mag contours

        # colorbar
        if axloc is None:
            cax = fig.add_axes([0.85, 0.1, 0.02, 5/6.])
        else:
            cax = fig.add_axes([axloc[0]+0.21, axloc[1], 0.01, 0.3333])
        try:
            ticks = MaxNLocator(n_ticks).tick_values(cbrange[0], cbrange[1])
        except AttributeError:
            print('AttributeError: MaxNLocator instance has no attribute' +
                  ' "tick_values" ')
            cb = fig.colorbar(p, cax)
        else:
            cb = fig.colorbar(p, cax, ticks=ticks)
        if cblabel is not None:
            cb.set_label(cblabel)

        # axis limits
        #xlim = ax.get_xlim()
        #ylim = ax.get_ylim()
        if xy_max is None:
            #xy_max = np.max(np.abs([xlim, ylim]))
            xy_max = self.xy_lim
        ax.axis([-xy_max, xy_max, -xy_max, xy_max])

        sns.set_style(rc={'axes.facecolor': '#EAEAF2'})



    def plot_spec(self,
                  rest_frame=True,
                  show_ivar=True,
                  colors=None,
                  kwargs={'alpha':0.75}, 
                  bin=0,
                  figsize=(20, 12)):
        '''
        Plot spectrum, model fits, and inverse variance.
        '''
        # Uncomment when you implement reading in multiple DAP analysis files
        # models, n_models = self.count_models(models)
        # c = generate_colors(colors, n_models)
        # labels = ['M11-STELIB', 'MILES']
        n_models = 1
        c = sns.color_palette('bright', n_models)
        labels = [self.tpl_lib]

        if rest_frame:
            wave = self.wave_rest[bin]
            gal = self.galaxy_rest
            ivar = self.ivar_rest
            models = [self.model_rest] 
        else:
            wave = self.wave_obs
            gal = self.galaxy
            ivar = self.ivar
            models = [self.model]

        sns.set_context('poster', rc={"lines.linewidth": 2})
        fig = plt.figure(figsize=figsize)
        ax = fig.add_axes([0.1, 0.1, 0.85, 0.8])
        ax.set_xlabel(r'$\lambda [\AA]$')
        ax.set_ylabel('Flux [10$^{-17}$ erg/s/cm$^2$]')
        ax.set_title('%s bin %i' % (self.manga_id, bin))
        pgal = ax.plot(wave, gal[bin], color='gray')[0]
        pmod = []
        for i in xrange(n_models):
            pmod.append(ax.plot(wave, models[i][bin], color=c[i], **kwargs)[0])
        xylim = ax.axis()
        ax.axis(xylim)
        if show_ivar:
            pivar = ax.plot(wave, 1./np.sqrt(ivar[bin]) * 10., color='r')[0]
        leg = plt.legend([pgal] + pmod + [pivar],
                         ['galaxy'] + labels + [r'error $\times$ 10'], loc=2)



    def plot_all_spec(self,
                      resid=False,
                      rest_frame=True,
                      xlim=None,
                      ylim=None,
                      figsize=(20, 12)):
        '''
        Plot all spectra (normalized to their median value).
        '''
        if rest_frame:
            gal = self.galaxy_rest
            lframe = 'rest frame'
        else:
            wave = self.wave_obs
            gal = self.galaxy
            lframe = 'observed frame'
            xlim = np.array(xlim) * (1. + self.z_median)

        if resid:
            ylabel = r'$\Delta$'
        else:
            ylabel = 'Normalized Flux [10$^{-17}$ erg/s/cm$^2$]'

        sns.set_context('poster', rc={"lines.linewidth": 2})
        fig = plt.figure(figsize=figsize)
        ax = fig.add_axes([0.1, 0.1, 0.85, 0.8])
        ax.set_title('%s' % (self.manga_id))
        ax.set_xlabel('%s ' % lframe + r'$\lambda [\AA]$')
        ax.set_ylabel(ylabel)       

        for bin in xrange(self.n_bins):
            if rest_frame:
                wave = self.wave_rest[bin]
            if resid:
                y = self.resid_pix[bin]
            else:
                y = gal[bin] / np.median(gal[bin])
            ax.plot(wave, y, color='k', lw=0.2)

        if xlim is not None:
            ax.set_xlim(xlim)
        if ylim is not None:
            ax.set_ylim(ylim)



    def plot_resid(self,
                   bin=0,
                   resid_model=False,
                   resid_err=False,
                   rest_frame=True,
                   kwargs={'alpha':0.75},
                   lw=2,
                   leg_lab=None,
                   xlim=None,
                   ylim=None,
                   masks=True,
                   figsize=(20, 12)):
        # Uncomment when you implement reading in multiple DAP analysis files
        # models, n_models = self.count_models(models)
        # c = generate_colors(colors, n_models)
        # labels = ['M11-STELIB', 'MILES']
        n_models = 1
        c = sns.color_palette('bright', n_models)
        labels = [self.tpl_lib]
        
        n_ax = 2
        if resid_model:
            n_ax += 1
        if resid_err:
            n_ax += 1

        if ylim is not None:
            ylim = np.array(ylim)
            if ylim.shape == (2,):
                ylim = np.ones((4, 2)) * ylim

        if rest_frame:
            wave = self.wave_rest[bin]
            gal = self.galaxy_rest
            ivar = self.ivar_rest
            models = [self.model_rest]
        else:
            wave = self.wave_obs
            gal = self.galaxy
            ivar = self.ivar
            models = [self.model]
        
        residuals = [gal - models[j] for j in xrange(n_models)]

        ind = np.where((wave >= xlim[0]) & (wave <= xlim[1]))[0]

        if masks:
            ind_split = np.where(np.diff(self.smsk[bin]) != 0)[0]
            smsk_vals = np.array_split(self.smsk[bin], ind_split + 1)
            smsk_val = np.array([item[0] for item in smsk_vals], dtype=int)
            if 0 not in ind_split:
                ind_split = np.append([0], ind_split)
            if self.smsk.shape[1] not in ind_split:
                ind_split = np.append(ind_split, [self.smsk.shape[1] - 1])


        sns.set_context('poster', rc={"lines.linewidth": lw})
        fig = plt.figure(figsize=figsize)

        ind_mod = [np.where(models[k][bin] > 0.)[0] for k in xrange(n_models)]
        for i in xrange(n_ax):
            if i == 0:
                bottom = 0.32
                height = 0.63
            else:
                bottom = 0.1
                height = 0.2
            ax = fig.add_axes([0.1, bottom, 0.85, height])

            if xlim is not None:
                ax.set_xlim(xlim)
            if ylim is not None:
                ax.set_ylim(ylim[i])

            if i in np.arange(n_ax-1):
                plt.setp(ax.get_xticklabels(), visible=False)

            if i == n_ax-1:
                ax.set_xlabel(r'$\lambda_{\rm rest} [\AA]$')
                ax.set_ylabel(r'$\Delta / \sigma$')

            if i != 0:
                ax.plot([xlim[0], xlim[1]], [0, 0], color='gray', alpha=0.75, 
                        lw=2, zorder=1)

            p = []

            # spectrum
            if i == 0:
                ax.set_title('%s bin %i' % (self.manga_id, bin))
                ax.set_ylabel('Flux [10$^{-17}$ erg/s/cm$^2$]')
                p.append(ax.plot(wave[ind], gal[bin][ind], color='gray')[0])
                y = models

            # residuals
            elif i == 1:
                ax.set_ylabel(r'$\Delta$')
                y = residuals

            # residuals / model
            elif i == 2:
                ax.set_ylabel(r'$\Delta$ / model')                    
                y = [residuals[j] / models[j] for j in xrange(n_models)]

            # residuals / error
            elif i == 3:
                y = residuals * np.sqrt(ivar[bin])
            
            for j in xrange(n_models):
                p.append(ax.plot(wave[ind_mod[j]], y[j][bin][ind_mod[j]],
                         color=c[j], **kwargs)[0])
                if i == 0:
                    plt.legend(p, ['galaxy'] + leg_lab, loc=2)

            # plot masks
            if masks:
                mskargs = dict(facecolor='gray', alpha=0.25)
                lglam = np.log10(wave)
                dlglam = (lglam[1] - lglam[0]) / 2.
                ylower, yupper = -100, 10000
                for m in xrange(len(smsk_val)):
                    x = np.array([10**(lglam[ind_split[m]+1] - dlglam),
                                 10**(lglam[ind_split[m+1]] + dlglam)])
                    if smsk_val[m]:
                        ax.fill_between(x, ylower, yupper, **mskargs)
