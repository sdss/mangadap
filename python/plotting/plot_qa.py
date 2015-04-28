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

from astropy.io import fits
from astropy.stats import sigma_clip

import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import ListedColormap
import matplotlib.patches as patches
from matplotlib.ticker import MaxNLocator

import warnings
try:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        import seaborn as sns
    seaborn_installed = True
except ImportError:
    seaborn_installed = False


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


class FitError(Exception):
    pass


class PlotQA:
    def __init__(self, filename):
        self.dap_file = filename
        self.manga_pid = ('-').join(self.dap_file.split('-')[1:3])
        self.dap_mode = self.dap_file.split('LOG')[1].split('_')[0]
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
        fin = fits.open(self.dap_file)
        self.drps = fin[fin.index_of('DRPS')].data
        self.bin_prop = fin[fin.index_of('BINS')].data
        self.wave_obs = fin[fin.index_of('WAVE')].data
        self.galaxy = fin[fin.index_of('FLUX')].data.T
        self.ivar = fin[fin.index_of('IVAR')].data.T
        self.smod = fin[fin.index_of('SMOD')].data.T
        self.stfit = fin[fin.index_of('STFIT')].data
        self.smsk = fin[fin.index_of('SMSK')].data.T
        self.elopar = fin[fin.index_of('ELOPAR')].data.T
        self.elofit = fin[fin.index_of('ELOFIT')].data.T
        self.elomew = fin[fin.index_of('ELOMEW')].data.T
        self.elomfb = fin[fin.index_of('ELOMFB')].data.T

        # spaxel, pixel and bin counts
        self.n_bins, self.n_pix = self.galaxy.shape
        self.n_spaxels = len(self.drps)
        if self.dap_mode == 'CUBE':
            self.sqrt_n_spaxels = np.sqrt(self.n_spaxels)
            if np.mod(self.sqrt_n_spaxels, 1) != 0:
                print('Sqrt of number of bins in cube is not an integer.')
            else:
                self.sqrt_n_spaxels = int(self.sqrt_n_spaxels)

        # DRP measurements
        (self.xpos, self.ypos, ig1, ig2, self.signal, self.noise, self.binvr,
            self.binid, self.binw) = \
            [self.drps[name] for name in self.drps.dtype.names]
        self.signal[np.where(self.signal == 0.)] = np.nan
        self.noise[np.where(self.noise == 1.)] = np.nan
        self.snr = self.signal / self.noise

        # bin properties
        (self.binxrl, self.binyru, self.binr, self.bina, self.binsn,
            self.nbin, self.binf) = \
            [self.bin_prop[name] for name in self.bin_prop.dtype.names]

        if (self.smod == 0.).all():
            raise FitError

        # stellar kinematics
        self.tplw, self.addpoly, self.multpoly, self.kin, self.kinerr, \
            self.rchi2 = [self.stfit[name] for name in self.stfit.dtype.names]
        if self.kin.shape[1] == 2:
            self.stvel, self.stvdisp = self.kin.T
            self.sth3 = self.stvel * np.nan
            self.sth4 = self.stvel * np.nan
        elif self.kin.shape[1] == 4:
            self.stvel, self.stvdisp, self.sth3, self.sth4 = self.kin.T
        
        # emission line kinematics
        self.emvel_ew, self.emvdisp_ew = self.elofit['KIN_EW'].T
        self.emvdisp_halpha_ew = self.elofit['sinst_ew'][:, 5]

        # emission line fluxes
        (self.oii3727_ew, self.hbeta_ew, self.oiii4959_ew, self.oiii5007_ew,
            self.nii6548_ew, self.halpha_ew, self.nii6583_ew, self.sii6717_ew,
            self.sii6731_ew) = self.elofit['FLUX_EW'].T

        # emission line + stellar continuum fits
        self.fullfitew = self.smod + self.elomew
        self.fullfitfb = self.smod + self.elomfb

        # Need to account for bins where certain emission lines weren't fit
        # self.elofit['ELOMIT_EW'] # = 1/0 where 1 is not fit and 0 is fit

        h = fin[0].header
        self.tpl_lib = h['TPLKEY']
        try:
            self.manga_id = h['MANGAID']
        except KeyError:
            self.manga_id = 'n/a'

        fin.close()


    # def read_val(self, extname, colname=None):
    #     fin = fits.open(self.dap_file)
    #     extval = fin[fin.index_of(extname)].data
    #     fin.close()
    # 
    #     if colname is not None:
    #         value = extval[colname]
    #     else:
    #         value = extval
    # 
    #     return value
 

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
        self.smod_rest = (self.smod.T * (1. + self.z)).T
        self.fullfitew_rest = (self.fullfitew.T * (1. + self.z)).T
        self.fullfitfb_rest = (self.fullfitfb.T * (1. + self.z)).T

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
        for i in range(n_regions):
            ind.append(np.where(
                       (self.wave_rest > self.lam_good[i, 0]) & 
                       (self.wave_rest < self.lam_good[i, 1])))
        
        ind0 = np.concatenate([ind[i][0] for i in range(n_regions)])
        ind1 = np.concatenate([ind[i][1] for i in range(n_regions)])
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
        for i in range(self.n_bins):
            indt = np.where(self.ind_lam[0] == i)[0]
            metric_bin[i] = np.nanmean(metric_pix[i][self.ind_lam[1][indt]])
        return metric_bin


    def calc_chisq(self):
        '''
        chisq_pix: chisq at every pixel
        chisq_bin: chisq for each binned spectrum
        '''
        self.chisq_pix = (self.galaxy - self.smod)**2. * self.ivar
        self.chisq_pix[np.where(self.smsk==1)] = np.nan
        self.chisq_bin = self.calc_metric_per_bin(self.chisq_pix)


    def calc_resid(self):
        '''
        resid_pix: residual at every pixel
        resid_bin: residual for each binned spectrum
        '''
        self.resid_pix = self.galaxy - self.smod
        self.resid_pix[np.where(self.smsk==1)] = np.nan
        self.resid_bin = self.calc_metric_per_bin(self.resid_pix)
        

    def calc_resid_mod(self):
        '''
        resid_mod_pix: residual / model at every pixel
        resid_mod_bin: residual / model for each binned spectrum
        '''
        self.resid_mod_pix = (self.galaxy - self.smod) / self.smod
        self.resid_mod_pix[np.where(self.smsk==1)] = np.nan
        self.resid_mod_bin = self.calc_metric_per_bin(self.resid_mod_pix)
        

    def calc_resid_err(self):
        '''
        resid_err_pix: residual / error at every pixel
        resid_err_bin: residual / error for each binned spectrum
        '''
        self.resid_err_pix = (self.galaxy - self.smod) * np.sqrt(self.ivar)
        self.resid_err_pix[np.where(self.smsk==1)] = np.nan
        self.resid_err_bin = self.calc_metric_per_bin(self.resid_err_pix)

    def calc_resid_data(self):
        '''
        resid_data_pix: residual / data at every pixel
        resid_data_bin: residual / data for each binned spectrum
        resid_data_bin_percent99: 99th percentile of resid_data_bin
        resid_data_bin_percent68: 68th percentile of resid_data_bin
        '''
        with np.errstate(divide='ignore', invalid='ignore'):
            self.resid_data_pix = np.abs(self.galaxy - self.smod) / self.galaxy
        self.resid_data_pix[np.where(self.smsk==1)] = np.nan
        self.resid_data_bin = self.calc_metric_per_bin(self.resid_data_pix)
        self.resid_data_bin_percent99 = np.array(
            [np.percentile(self.resid_data_pix[i][np.where(self.smsk[i]==0)], 99.)
            for i in range(self.n_bins)])
        self.resid_data_bin_percent68 = np.array(
            [np.percentile(self.resid_data_pix[i][np.where(self.smsk[i]==0)], 68.)
            for i in range(self.n_bins)])

    #-------------------------------------------------------------------------
    

    #----- FITS file ---------------------------------------------------------
    #def fits_out(self):


    #-------------------------------------------------------------------------



    #----- Plot Setup --------------------------------------------------------

    def plot_setup(self):
        self.cubehelix111, self.cubehelix111_r = \
            self.define_custom_cubehelix(rot=1, start=1, gamma=1)
        self.linearL, self.linearL_r = self.linear_Lab()
        self.set_axis_lims()
        self.cm_gray = mpl.colors.ListedColormap(['#303030'])

    def reverse_cmap(self, x):
        def out(y):
            return x(1. - y)
        return out


    def define_custom_cubehelix(self, rot=-1.5, start=0.5, gamma=1):
        cdict = mpl._cm.cubehelix(gamma=gamma, s=start, r=rot)
        c_r_dict = {}
        for k in iterkeys(cdict):
            c_r_dict[k] = self.reverse_cmap(cdict[k])
        cmap = LinearSegmentedColormap('cubehelix_cust', cdict)
        cmap_r = LinearSegmentedColormap('cubehelix_cust_r', c_r_dict)
        return (cmap, cmap_r)


    def linear_Lab(self):
        try:
            LinL_file = os.path.join(os.environ['MANGADAP_DIR'], 'python', 'plotting',
                                     'Linear_L_0-1.csv')
            LinL = np.loadtxt(LinL_file, delimiter=',')
        except:
            print('Could not find or open Linear_L_0-1.csv')
            return self.define_custom_cubehelix()

        b3 = LinL[:, 2] # value of blue at sample n
        b2 = LinL[:, 2] # value of blue at sample n
        b1 = np.linspace(0, 1, len(b2)) # position of sample n - ranges from 0 to 1
        
        # setting up columns for list
        g3 = LinL[:, 1]
        g2 = LinL[:, 1]
        g1 = np.linspace(0, 1, len(g2))
        
        r3 = LinL[:, 0]
        r2 = LinL[:, 0]
        r1 = np.linspace(0, 1, len(r2))
        
        # creating list
        R = zip(r1, r2, r3)
        G = zip(g1, g2, g3)
        B = zip(b1, b2, b3)
        
        # transposing list
        RGB = zip(R, G, B)
        rgb = zip(*RGB)
        
        # creating dictionary
        k = ['red', 'green', 'blue']
        LinearL = dict(zip(k, rgb)) # makes a dictionary from 2 lists

        LinearL_r = {}
        for k in iterkeys(LinearL):
            LinearL_r[k] = self.reverse_cmap(LinearL[k])

        cmap = LinearSegmentedColormap('linearL', LinearL)
        cmap_r = LinearSegmentedColormap('linearL_r', LinearL_r)

        return (cmap, cmap_r)

    #-------------------------------------------------------------------------



    #----- Plots -------------------------------------------------------------
    

    def find_map_lims(x):
        return np.array([np.floor(np.min(x)), np.ceil(np.max(x))])


    def set_axis_lims(self):
        ind_b = np.where(self.binid >= 0)
        self.xy_lim = np.max((np.abs(self.xpos[ind_b]), 
                             np.abs(self.ypos[ind_b]))) * 1.1
        
        # Don't need to set these unless comparing different analyses of the
        # same galaxy

        #self.chisq_lim = self.find_map_lims(self.chisq_bin)
        #self.stvel_lim = 
        #self.stvdisp_lim = 
        #self.sth3_lim = 
        #self.sth4_lim = 
        

    # def count_models(x):
    #     if x is not None:
    #         try:
    #             if len(x.shape) == 2:
    #                 x = [x]
    #                 n = 1
    #         except AttributeError:
    #             n = len(x)
    #     else:
    #         n = 0
    #     return x, n
    
    
    # def generate_colors(colors=None, n_models=0):
    #     if colors is None:
    #         n_colors = 1
    #         if n_models > 0:
    #             n_colors = n_models
    #         ctmp = sns.color_palette('bright', n_colors)
    #         colors = ctmp
    #         colors[1] = ctmp[2]
    #         colors[2] = ctmp[1]
    #     return colors
    

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
        if seaborn_installed:
            sns.set_context('poster', rc={'lines.linewidth': 2})
        bigAxes = fig.add_axes([0.04, 0.05, 0.9, 0.88], frameon=False)
        bigAxes.set_xticks([])
        bigAxes.set_yticks([])
        bigAxes.set_xlabel('arcsec', fontsize=20)
        bigAxes.set_ylabel('arcsec', fontsize=20)
        bigAxes.set_title(
            'pid-ifu %s     manga-id %s' % (self.manga_pid, self.manga_id),
            fontsize=20)

        for i in range(n_ax):
            dx = 0.31 * i
            dy = 0.45
            if i >= (n_ax / 2):
                dx = 0.31 * (i - n_ax / 2)
                dy = 0
            left, bottom = (0.08+dx, 0.1+dy)
            if seaborn_installed:
                sns.set_context('talk', rc={'lines.linewidth': 2})
                sns.set_style(rc={'axes.facecolor': '#A8A8A8'})

            ax = fig.add_axes([left, bottom, 0.2, 0.3333])

            if not seaborn_installed:
                ax.set_axis_bgcolor('#A8A8A8')
                ax.grid(False, which='both', axis='both')

            d = args[map_order[i]]
            k = d['kwargs']
            if not np.isnan(d['val']).all():
                self.plot_map(d['val'], fig=fig, ax=ax, axloc=[left, bottom],
                              **d['kwargs'])
                # self.plot_map(d['val'], fig=fig, ax=ax, axloc=[left, bottom],
                #               **d['kwargs'])





    def plot_map(self,
                 z,
                 val_no_measure=None,
                 interpolated=False,
                 flux=None,
                 cblabel=None,
                 cbrange=None,
                 cbrange_clip=True,
                 cbrange_symmetric=False,
                 n_ticks=7,
                 cmap=cm.coolwarm,
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
        if seaborn_installed:
            if ax is None:
                sns.set_context('poster', rc={'lines.linewidth': 2})
            else:
                sns.set_context('talk', rc={'lines.linewidth': 2})
            sns.set_style(rc={'axes.facecolor': '#A8A8A8'})

        if ax is None:
            fig = plt.figure(figsize=figsize)
            ax = fig.add_axes([0.12, 0.1, 2/3., 5/6.])
            #ax = fig.add_axes([0.15, 0.15, 0.6, 0.75])

            ax.set_xlabel('arcscec')
            ax.set_ylabel('arcsec')

        if title_text is not None:
            ax.set_title(title_text, fontsize=20)

        if not seaborn_installed:
            ax.set_axis_bgcolor('#A8A8A8')

        ax.grid(False, which='both', axis='both')

        kwargs = {}

        ind_tmp = np.where(self.binid >= 0)[0]
        ind_b = self.binid[ind_tmp]
        bin_num = np.arange(len(z), dtype=int).astype(str)

        spaxel_size = 0.5 # arcsec
        delta = spaxel_size / 2.

        if not nodots:
            ax.plot(-self.binxrl, self.binyru, '.k', markersize=3, zorder=10)
        
        if val_no_measure is None:
            val_no_measure = 0.
        
        # colorbar range
        if len(z) == self.n_spaxels:
            ind_cb = ind_tmp
        elif len(z) == self.n_bins:
            ind_cb = ind_b

        cbrange_tmp = [z[ind_cb].min(), z[ind_cb].max()]

        if cbrange_clip:
            zclip = sigma_clip(z[ind_cb], sig=3)
            cbrange_tmp = [zclip.min(), zclip.max()]

        if cbrange is None:
            cbrange = cbrange_tmp
        else:
            for i in range(len(cbrange)):
                if cbrange[i] is None:
                    cbrange[i] = cbrange_tmp[i]
                else:
                    if i == 0:
                        kwargs['vmin'] = cbrange[i]
                    elif i == 1:
                        kwargs['vmax'] = cbrange[i]

        if cbrange_symmetric:
            cb_max = np.max(np.abs(cbrange))
            cbrange = [-cb_max, cb_max]

        cbdelta = cbrange[1] - cbrange[0]


        #----- Re-orient Data -----

        if len(z) == self.n_spaxels:
            im1 = np.reshape(z, (self.sqrt_n_spaxels, self.sqrt_n_spaxels))
            im = im1.T[::-1]

        elif len(z) == self.n_bins:
            # xpos and ypos are oriented to start in the lower left and go up
            # then right
            # Re-orient data so that it starts in the upper left and goes
            # right then down

            xpos1 = np.reshape(self.xpos, (self.sqrt_n_spaxels,
                               self.sqrt_n_spaxels))
            xpos2 = xpos1.T[::-1]

            ypos1 = np.reshape(self.ypos, (self.sqrt_n_spaxels,
                               self.sqrt_n_spaxels))
            ypos2 = ypos1.T[::-1]

            binid1 = np.reshape(self.binid, (self.sqrt_n_spaxels,
                                self.sqrt_n_spaxels))
            binid2 = binid1.T[::-1]
            binid3 = np.ma.array(binid2, mask=(binid2 == -1))

            # create a masked array of the data
            im = np.empty((self.sqrt_n_spaxels, self.sqrt_n_spaxels)) * np.nan
            for i in range(self.sqrt_n_spaxels):
                for j in range(self.sqrt_n_spaxels):
                    if not binid3.mask[i, j]:
                        im[i, j] = z[binid3.data[i, j]]

        ind_no_measure = np.where(im == val_no_measure)
        mask_no_data = np.isnan(im)
        mask_no_data[im == val_no_measure] = True
        im_mask_no_data = np.ma.array(im, mask=mask_no_data)


        # plot map
        if interpolated:
            levels = np.linspace(cbrange[0], cbrange[1], n_colors)
            #ind_cb = np.arange(300)
            p = ax.tricontourf(-self.binxrl[ind_cb], self.binyru[ind_cb],
                               z[ind_cb], levels=levels, cmap=cmap)
        else:
            extent = np.array([-self.xpos.max() - delta,
                              -self.xpos.min() + delta,
                              self.ypos.min() - delta,
                              self.ypos.max() + delta])

            if len(ind_no_measure[0]) > 0:
                for i0, i1 in zip(ind_no_measure[0], ind_no_measure[1]):
                    ax.add_patch(mpl.patches.Rectangle(
                        (-(xpos2[i0, i1] + delta), ypos2[i0, i1] - delta),
                        spaxel_size, spaxel_size, hatch='///', linewidth=5,
                        fill=True, fc='#D6D6E5', zorder=10))

            if cbrange_clip:
                if 'vmin' not in kwargs.keys():
                    kwargs['vmin'] = zclip.min()
                if 'vmax' not in kwargs.keys():
                    kwargs['vmax'] = zclip.max()

            p = ax.imshow(im_mask_no_data, interpolation='none',
                          extent=extent, cmap=cmap, **kwargs)


        if spaxel_num:
            for i in range(self.n_bins):
                if (i >= 100) and (i < 1000) and (self.nbin[i] <= 2):
                    fontsize_tmp = fontsize - 1
                elif (i >= 1000) and (self.nbin[i] <= 2):
                    fontsize_tmp = fontsize - 2
                elif self.nbin[i] == 1:
                    fontsize_tmp = fontsize
                else:
                    fontsize_tmp = fontsize + 2
                ctmp = cmap((z[i] - cbrange[0]) / cbdelta)
                ax.text(-(self.binxrl[i] + (spaxel_size / 4)),
                        self.binyru[i] - (spaxel_size / 4),
                        str(i), fontsize=fontsize_tmp,
                        color=(1.-ctmp[0], 1.-ctmp[1], 1.-ctmp[2], ctmp[3]))

        # overplot flux contours (which flux? in what band?)
        if flux is not None:
            ax.tricontour(-self.binxrl, self.binyru,
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

        if seaborn_installed:
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

        if seaborn_installed:
            c = sns.color_palette('bright', n_models)
        else:
            c = ['b', 'lime', 'r', 'DarkOrchid', 'gold']

        labels = [self.tpl_lib]

        if rest_frame:
            wave = self.wave_rest[bin]
            gal = self.galaxy_rest
            ivar = self.ivar_rest
            models = [self.smod_rest] 
        else:
            wave = self.wave_obs
            gal = self.galaxy
            ivar = self.ivar
            models = [self.smod]

        if seaborn_installed:
            sns.set_context('poster', rc={"lines.linewidth": 2})

        fig = plt.figure(figsize=figsize)
        ax = fig.add_axes([0.1, 0.1, 0.85, 0.8])
        ax.set_xlabel(r'$\lambda [\AA]$')
        ax.set_ylabel('Flux [10$^{-17}$ erg/s/cm$^2$]')
        ax.set_title('pid-ifu %s     manga-id %s     bin %i' % (self.manga_pid,
                     self.manga_id, bin))
        pgal = ax.plot(wave, gal[bin], color='gray')[0]
        pmod = []
        for i in range(n_models):
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

        if seaborn_installed:
            sns.set_context('poster', rc={"lines.linewidth": 2})
            sns.set_style(rc={'axes.facecolor': '#EAEAF2'})

        fig = plt.figure(figsize=figsize)
        ax = fig.add_axes([0.1, 0.1, 0.85, 0.8])
        ax.set_title('pid-ifu %s     manga-id %s' % (self.manga_pid, self.manga_id))
        ax.set_xlabel('%s ' % lframe + r'$\lambda [\AA]$')
        ax.set_ylabel(ylabel)       

        for bin in range(self.n_bins):
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
                   rest_frame=True,
                   kwargs={'alpha':0.75},
                   lw=1,
                   xlim=None,
                   ylim=None,
                   masks=True,
                   figsize=(20, 12)):

        if seaborn_installed:
            c = sns.color_palette('bright', 5)
            sns.set_context('poster', rc={"lines.linewidth": lw})
        else:
            c = ['b', 'lime', 'r', 'DarkOrchid', 'gold']

        n_ax = 2

        if rest_frame:
            wave = self.wave_rest[bin]
            gal = self.galaxy_rest[bin]
            ivar = self.ivar_rest[bin]
            stmodel = self.smod_rest[bin]
            fullfitew = self.fullfitew_rest[bin]
            fullfitfb = self.fullfitfb_rest[bin]
        else:
            wave = self.wave_obs
            gal = self.galaxy[bin]
            ivar = self.ivar[bin]
            stmodel = self.smod[bin]
            fullfitfb = self.fullfitfb[bin]
        
        residew = gal - fullfitew
        residfb = gal - fullfitfb

        ind_stmodel = np.where(stmodel > 0.)[0]
        ind = np.where((wave >= xlim[0]) & (wave <= xlim[1]))[0]

        ylim_tmp = copy.deepcopy(ylim)
        if ylim_tmp[0] is not None:
            if ylim_tmp[0][1] is None:
                p50 = np.percentile(gal[ind], 50)
                ylim_tmp[0][1] = p50 * 3.
                if ylim_tmp[1] is None:
                    ymax_tmp = 0.1 + p50 * 0.2
                    ylim_tmp[1] = [-ymax_tmp, ymax_tmp]

        if masks:
            ind_split = np.where(np.diff(self.smsk[bin]) != 0)[0]
            smsk_vals = np.array_split(self.smsk[bin], ind_split + 1)
            smsk_val = np.array([item[0] for item in smsk_vals], dtype=int)
            if 0 not in ind_split:
                ind_split = np.append([0], ind_split)
            if self.smsk.shape[1] not in ind_split:
                ind_split = np.append(ind_split, [self.smsk.shape[1] - 1])


        fig = plt.figure(figsize=figsize)

        for i in range(n_ax):
            if i == 0:
                bottom = 0.32
                height = 0.63
            else:
                bottom = 0.1
                height = 0.2
            ax = fig.add_axes([0.1, bottom, 0.85, height])

            if xlim is not None:
                ax.set_xlim(xlim)
            if ylim_tmp is not None:
                ax.set_ylim(ylim_tmp[i])

            if i in np.arange(n_ax-1):
                plt.setp(ax.get_xticklabels(), visible=False)

            if i == n_ax-1:
                ax.set_xlabel(r'$\lambda_{\rm rest} [\AA]$')
                ax.set_ylabel(r'$\Delta / \sigma$')

            if i != 0:
                ax.plot([xlim[0], xlim[1]], [0, 0], color='gray', alpha=0.75, 
                        lw=2, zorder=1)

            p = []
            labels = []

            # spectrum
            if i == 0:
                ax.set_title('pid-ifu %s     manga-id %s     bin %i' % (
                             self.manga_pid, self.manga_id, bin))
                flux_units = 'Flux'
                if self.dap_mode == 'CUBE':
                    flux_units += ' / spaxel'
                elif self.dap_mode == 'RSS':
                    flux_units += ' / fiber'
                ax.set_ylabel('%s [10$^{-17}$ erg/s/cm$^2$]' % flux_units)
                p.append(ax.plot(wave[ind], gal[ind], color='#808080')[0])
                labels.append('galaxy')

                # stellar continuum fit
                p.append(ax.plot(wave[ind_stmodel], stmodel[ind_stmodel],
                         color=c[1])[0])
                labels.append('stellar continuum fit')

                y = [fullfitfb, fullfitew]
                fit_labels = ['F. Belfiore emline + cont fit', 
                              'E. Wang emline + cont fit']

            #  residuals
            elif i == 1:
                ax.set_ylabel(r'$\Delta$')
                y = [residfb, residew]

            # emission line + stellar continuum fits
            for j in range(len(y)):
                p.append(ax.plot(wave[ind_stmodel], y[j][ind_stmodel],
                         color=c[2-2*j])[0])
                if i == 0:
                    labels.append(fit_labels[j])

            if i == 0:
                plt.legend(p, labels, loc=2)

            # plot masks
            if masks:
                mskargs = dict(facecolor='gray', alpha=0.25)
                lglam = np.log10(wave)
                dlglam = (lglam[1] - lglam[0]) / 2.
                ylower, yupper = ylim_tmp[i]
                for m in range(len(smsk_val)):
                    x = np.array([10**(lglam[ind_split[m]+1] - dlglam),
                                 10**(lglam[ind_split[m+1]] + dlglam)])
                    if smsk_val[m]:
                        ax.fill_between(x, ylower, yupper, **mskargs)

        return fig
