# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""Functions for plotting DAP output."""

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import sys
import copy

import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.ticker import MaxNLocator

from astropy.stats import sigma_clip

from mangadap.plot import util

import warnings
try:
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        import seaborn as sns
except ImportError:
    print('Seaborn could not be imported. Continuing...')

def set_extent(xpos, ypos, delta):
    """Set extent of map."""
    return np.array([-xpos.max() - delta, -xpos.min() + delta,
                     ypos.min() - delta, ypos.max() + delta])

def reorient(x):
    """Reorient XPOS and YPOS.

    XPOS and YPOS are oriented to start in the lower left and go up then
    right. Re-orient data so that it starts in the upper left and goes right
    then down.

    Args:
        x (array): Positional values.
        sqrtn (int): Square root of the number of spaxels.

    Returns:
        array: Square 2-D array of size (sqrtn, sqrtn).
    """
    sqrtn = int(np.sqrt(len(x)))
    x_sq = np.reshape(x, (sqrtn, sqrtn))
    return x_sq.T[::-1]

def make_image(val, err, xpos, ypos, binid, delta, val_no_measure,
               snr_thresh):
    """Make masked array of image.

    Args:
        val (array): Values.
        err (array): Errors.
        xpos (array): x-coordinates of bins.
        ypos (array): y-coordinates of bins.
        binid (array): Bin ID numbers.
        delta (float): Half of the spaxel size in arcsec.
        val_no_measure (float): Value that corresponds to no measurement.
        snr_thresh (float): Signal-to-noise theshold below which is not
           considered a measurement.

    Returns:
        tuple: (masked array of image,
                tuple of (x, y) coordinates of bins with no measurement)
    """
    # create a masked array of the data
    im = val[binid]
    im[binid == -1] = np.nan

    im_err = err[binid]
    im_err[binid == -1] = np.nan

    no_measure = make_mask_no_measurement(im, im_err, val_no_measure,
                                          snr_thresh)
    no_data = make_mask_no_data(im, no_measure)
    image = np.ma.array(im, mask=no_data)

    # spaxels with data but no measurement
    xpos_re = reorient(xpos)
    ypos_re = reorient(ypos)
    xy_nomeasure = (-(xpos_re[no_measure]+delta), ypos_re[no_measure]-delta)
    return image, xy_nomeasure

def make_mask_no_measurement(data, err, val_no_measure, snr_thresh):
    """Mask invalid measurements within a data array.

    Args:
        data (array): Data.
        err (array): Error.
        val_no_measure (float): Value in data array that corresponds to no
            measurement.
        snr_thresh (float): Signal-to-noise threshold for keeping a valid
            measurement.

    Returns:
        array: Boolean array for mask (i.e., True corresponds to value to be
            masked out).
    """
    no_measure = (data == val_no_measure)
    if all([[i is None for i in x] for x in err]):
        err = None
    if err is not None:
        no_measure[err == 0.] = True
        no_measure[(np.abs(data / err) < snr_thresh)] = True
    return no_measure

def make_mask_no_data(data, mask_no_measurement):
    """Mask entries with no data or invalid measurements.

    Args:
        data (array): Data.
        mask_no_measure (array): Mask for entries with no measurement.

    Returns:
        array: Boolean array for mask (i.e., True corresponds to value to be
            masked out).
    """
    data = data.astype(float)
    no_data = np.isnan(data)
    no_data[mask_no_measurement] = True
    return no_data

def set_vmin_vmax(d, cbrange):
    """Set minimum and maximum values of the color map."""
    if 'vmin' not in d.keys():
        d['vmin'] = cbrange[0]
    if 'vmax' not in d.keys():
        d['vmax'] = cbrange[1]
    return d

def cbrange_sigclip(image, sigma):
    """Sigma clip colorbar range.

    Args:
        image (masked array): Image.
        sigma (float): Sigma to clip.

    Returns:
        list: Colorbar range.
    """
    imclip = sigma_clip(image.data[~image.mask], sig=sigma)
    try:
        cbrange = [imclip.min(), imclip.max()]
    except ValueError:
        cbrange = [image.min(), image.max()]
    return cbrange

def cbrange_user_defined(cbrange, cbrange_user):
    """Set user-specified colorbar range.

    Args:
        cbrange (list): Input colorbar range.
        cbrange_user (list): User-specified colorbar range. If a value is
            None, then the colorbar uses the previous value.

    Returns:
        list: Colorbar range.
    """
    for i in range(2):
        if cbrange_user[i] is not None:
            cbrange[i] = cbrange_user[i]
    return cbrange

def set_cbrange(image, cbrange=None, sigclip=None, symmetric=False):
    """Set colorbar range.

    Args:
        image (masked array): Image.
        cbrange (list): User-specified colorbar range. Default is None.
        sigclip (float): Sigma value for sigma clipping. If None, then do not
            clip. Default is None.
        symmetric (boolean): If True, make colorbar symmetric around zero.
            Default is False.

    Returns:
        list: Colorbar range.
    """

    if sigclip is not None:
        cbr = cbrange_sigclip(image, sigclip)
    else:
        cbr = [image.min(), image.max()]

    if cbrange is not None:
        cbr = cbrange_user_defined(cbr, cbrange)

    if symmetric:
        cb_max = np.max(np.abs(cbr))
        cbr = [-cb_max, cb_max]

    return cbr

def make_draw_colorbar_kws(image, cb_kws):
    """Make keyword args dictionary to pass to draw_colorbar.

    Args:
        image (masked array): Image to display.
        cb_kws (dict): Keyword args to set and draw colorbar.

    Returns:
        dict: draw_colorbar keyword args
    """
    keys = ('cbrange', 'sigclip', 'symmetric')
    cbrange_kws = {k: cb_kws.pop(k, None) for k in keys}
    cb_kws['cbrange'] = set_cbrange(image, **cbrange_kws)
    return cb_kws

def draw_colorbar(fig, mappable, axloc=None, cbrange=None, n_ticks=7,
                  label_kws=None, tick_params_kws=None):
    """Make colorbar.

    Args:
        fig: plt.figure object.
        mappable: Plotting element to map to colorbar.
        axloc (list): Specify (left, bottom, width, height) of colorbar axis.
            Defaults to None.
        cbrange (list): Colorbar min and max.
        n_ticks (int): Number of ticks on colorbar.
        label_kws (dict): Keyword args to set colorbar label. Default is None.
        tick_params_kws (dict): Keyword args to set colorbar tick parameters.
            Default is None.

    Returns:
        tuple: (plt.figure object, plt.figure axis object)
    """
    label_kws = util.none_to_empty_dict(label_kws)
    tick_params_kws = util.none_to_empty_dict(tick_params_kws)

    if axloc is not None:
        cax = fig.add_axes(axloc)
    else:
        cax = None

    try:
        ticks = MaxNLocator(n_ticks).tick_values(*cbrange)
    except AttributeError:
        print('AttributeError: MaxNLocator instance has no attribute'
              ' "tick_values" ')
        cb = fig.colorbar(mappable, cax)
    else:
        cb = fig.colorbar(mappable, cax, ticks=ticks)

    if label_kws['label'] is not None:
        cb.set_label(**label_kws)

    if tick_params_kws is not None:
        cb.ax.tick_params(**tick_params_kws)

    return fig, cb

def pretty_specind_units(units):
    """Convert units of spectral index to colorbar label.

    Args:
        units (str): 'ang' or 'mag'

    Returns:
        str
    """
    if units == 'ang':
        cblabel = r'$\AA$'
    elif units == 'mag':
        cblabel = 'Mag'
    else:
        raise ValueError('Unknown spectral index units.')
    return cblabel

def set_cmaps(cmaps, n_plots):
    if cmaps is None:
        try:
            cmap, cmap_r = util.linear_Lab()
        except IOError:
            cmap = cm.Blues_r
        cmaps = [cmap for _ in range(n_plots)]

    if len(cmaps) == 1:
        cmaps = [cmaps[0] for _ in range(n_plots)]

    return cmaps

def set_map_background_color(spaxel_size, color='#A8A8A8'):
    """Set default parameters for a single panel plot.

    Args:
        spaxel_size (float): Size of spaxels in arcsec.
        color (str): Background color. Default is '#A8A8A8' (gray).

    Returns:
        tuple: axes keyword args, patch keyword args
    """
    ax_kws = dict(facecolor=color)
    patch_kws = dict(width=spaxel_size, height=spaxel_size, hatch='xx',
                     linewidth=0, fill=True, fc=color, ec='w', zorder=10)
    return ax_kws, patch_kws

def set_map_par(cmap, title, cblabel, titlefontsize=28, cbfontsize=20):
    """Set default parameters for a single panel plot.

    Args:
        cmap (list): Colormap.
        title (str): Plot title.
        cblabel (str): Color bar label.
        titlefontsize (int): Title font size. Default is 28.
        cbfontsize (int): Color bar label and color bar tick label font size.
            Default is 20.

    Returns:
        tuple: title keyword args, imshow keyword args, colorbar keyword args
    """
    imshow_kws = dict(cmap=cmap)
    title_kws = dict(fontsize=titlefontsize, label=title)
    cb_kws = dict(axloc=[0.82, 0.1, 0.02, 5/6.],
                  cbrange=None, sigclip=3, symmetric=True,
                  label_kws=dict(label=cblabel, size=cbfontsize),
                  tick_params_kws=dict(labelsize=cbfontsize))
    return title_kws, imshow_kws, cb_kws

def map_ax_setup(fig=None, ax=None, fig_kws=None, facecolor='#EAEAF2'):
    """Basic axes setup for maps.

    Args:
        fig: Matplotlib plt.figure object. Use if creating subplot of a
            multi-panel plot. Default is None.
        ax: Matplotlib plt.figure axis object. Use if creating subplot of a
            multi-panel plot. Default is None.
        fig_kws (dict): Keyword args to pass to plt.figure. Default is None.
        facecolor (str): Axis facecolor. Default is '#EAEAF2'.

    Returns:
        tuple: (plt.figure object, plt.figure axis object)
    """
    fig_kws = util.none_to_empty_dict(fig_kws)

    if 'seaborn' in sys.modules:
        if ax is None:
            sns.set_context('poster', rc={'lines.linewidth': 2})
        else:
            sns.set_context('talk', rc={'lines.linewidth': 2})
        sns.set_style(rc={'axes.facecolor': facecolor})

    if ax is None:
        fig = plt.figure(**fig_kws)
        ax = fig.add_axes([0.12, 0.1, 2/3., 5/6.])
        ax.set_xlabel('arcscec')
        ax.set_ylabel('arcsec')

    if 'seaborn' not in sys.modules:
        ax.set_axis_bgcolor(facecolor)

    ax.grid(False, which='both', axis='both')
    return fig, ax

def make_map_title(analysis_id):
    """Make a map title to identify galaxy and analysis run.

    Args:
        analysis_id (dict): Identifying info about galaxy and analysis run.

    Returns:
        str: Map title.
    """
    return '     '.join(('pid-ifu {plate}-{ifudesign}', 'manga-id {mangaid}',
                         '{bintype}-{niter}')).format(**analysis_id)

def make_big_axes(fig, axloc=(0.04, 0.05, 0.9, 0.88), xlabel=None, ylabel=None,
                  labelsize=20, title_kws=None, mg_kws=None):
    """Make a big axes for x- and y-labels and a large plot title.

    Args:
        fig: Matplotlib plt.figure object.
        axloc (list): Specify [left, bottom, width, height] of axis. Defaults
            to [0.04, 0.05, 0.9, 0.88].
        xlabel: x-axis label. Defaults to None.
        ylabel: y-axis label. Defaults to None.
        labelsize (int): x- and y-axis labelsize. Defaults to 20.
        title_kws (dict): Keyword args to pass to ax.set_title. Default is None
        mg_kws (dict): MaNGA analysis ID keyword args to pass to
            make_map_title. Default is None.

    Returns:
        plt.figure axis object
    """
    title_kws = util.none_to_empty_dict(title_kws)
    mg_kws = util.none_to_empty_dict(mg_kws)

    # make axis without a frame or x and y ticks
    bigAxes = fig.add_axes(axloc, frameon=False)
    bigAxes.set_xticks([])
    bigAxes.set_yticks([])

    if xlabel is not None:
        bigAxes.set_xlabel(xlabel, fontsize=labelsize)
    if ylabel is not None:
        bigAxes.set_ylabel(ylabel, fontsize=labelsize)

    # set title
    if ('label' not in title_kws) and (mg_kws is not None):
        title_kws['label'] = make_map_title(mg_kws)
    if 'label' in title_kws:
        bigAxes.set_title(**title_kws)
    return bigAxes

def show_bin_num(dapdata, val, ax, imshow_kws, fontsize=6):
    """Display bin number on map.

    Args:
        dapdata:
        val (array):
        ax:
        image:
        fontsize (int): Nominal font size. Defaults to 6.

    Returns:
        plt.figure axis object
    """
    binxrl = dapdata.bins.binxrl.values
    binyru = dapdata.bins.binyru.values
    nbin = dapdata.bins.nbin.values
    for i, (x, y, nb, v) in enumerate(zip(binxrl, binyru, nbin, val)):
        fontsize_tmp = set_bin_num_fontsize(fontsize, i, nb)
        color = set_bin_num_color(v, imshow_kws)
        ax.text(-x, y, str(i), fontsize=fontsize_tmp, color=color,
                horizontalalignment='center', verticalalignment='center',
                zorder=10)
    return ax

def set_bin_num_color(value, imshow_kws):
    """Set color of bin numbers to be inverse of bin color.

    Args:
        value (float): Map value for a bin.
        imshow_kws (dict): Keyword args passed to ax.imshow that create the
            colormap.

    Returns:
        tuple: Text color for bin number.
    """
    cmap, vmin, vmax = (imshow_kws[k] for k in ['cmap', 'vmin', 'vmax'])
    ctmp = cmap((value - vmin) / (vmax - vmin))
    color = (1.-ctmp[0], 1.-ctmp[1], 1.-ctmp[2], ctmp[3])  # invert
    return color

def set_bin_num_fontsize(fontsize, number, nbin):
    """Set font size of bin numbers.

    Args:
        fontsize (int): Nominal font size. Defaults to 7.
        number (int): Bin number.
        nbin (int): Number of bins.

    Returns:
        int: fontsize
    """
    if (number >= 10) and (number < 100) and (nbin <= 2):
        fontsize_out = fontsize - 2
    elif (number >= 100) and (number < 1000) and (nbin <= 2):
        fontsize_out = fontsize - 3
    elif (number >= 1000) and (nbin <= 2):
        fontsize_out = fontsize - 4
    elif nbin == 1:
        fontsize_out = fontsize
    else:
        fontsize_out = fontsize + 2
    return fontsize_out


def plot_bindot(binxrl, binyru):
    """SNIPPETS ONLY
    """
    bindot_args = ()
    if show_bindot:
        bindot_args = (-binxrl, binyru)


def plot_map(image, extent, xy_nomeasure=None, fig=None, ax=None,
             dapdata=None, fig_kws=None, ax_kws=None, title_kws=None,
             patch_kws=None, imshow_kws=None, cb_kws=None, binnum_kws=None,
             bindot_args=()):
    """Make map.

    Args:
        image (masked array): Image to display.
        extent (array): Minimum and maximum x- and y-values.
        xy_nomeasure (tuple): x- and y-coordinates of spaxels without
            measurements.
        fig: Matplotlib plt.figure object. Use if creating subplot of a
            multi-panel plot. Default is None.
        ax: Matplotlib plt.figure axis object. Use if creating subplot of a
            multi-panel plot. Default is None.
        fig_kws (dict): Keyword args to pass to plt.figure. Default is None.
        ax_kws (dict): Keyword args to draw axis. Default is None.
        title_kws (dict): Keyword args to pass to ax.set_title. Default is
            None.
        patch_kws (dict): Keyword args to pass to ax.add_patch. Default is
            None.
        imshow_kws (dict): Keyword args to pass to ax.imshow. Default is None.
        cb_kws (dict): Keyword args to set and draw colorbar. Default is None.
        binnum_kws (dict): Keyword args to pass to show_bin_num. Default is
            None.
        bindot_args (tuple): x- and y-coordinates of bins.

    Returns:
        tuple: (plt.figure object, plt.figure axis object)
    """
    fig_kws = util.none_to_empty_dict(fig_kws)
    ax_kws = util.none_to_empty_dict(ax_kws)
    title_kws = util.none_to_empty_dict(title_kws)
    patch_kws = util.none_to_empty_dict(patch_kws)
    imshow_kws = util.none_to_empty_dict(imshow_kws)
    cb_kws = util.none_to_empty_dict(cb_kws)
    binnum_kws = util.none_to_empty_dict(binnum_kws)

    fig, ax = map_ax_setup(fig, ax, fig_kws=fig_kws, **ax_kws)

    if title_kws['label'] is not None:
        ax.set_title(**title_kws)

    drawcb_kws = make_draw_colorbar_kws(image, cb_kws)
    imshow_kws = set_vmin_vmax(imshow_kws, drawcb_kws['cbrange'])

    # Plot regions with no measurement as hatched
    if xy_nomeasure is not None:
        for xh, yh in zip(*xy_nomeasure):
            ax.add_patch(mpl.patches.Rectangle((xh, yh), **patch_kws))

    p = ax.imshow(image, interpolation='none', extent=extent, **imshow_kws)

    fig, cb = draw_colorbar(fig, p, **drawcb_kws)

    if binnum_kws:
        ax = show_bin_num(dapdata=dapdata, ax=ax, imshow_kws=imshow_kws,
                          **binnum_kws)

    if bindot_args:
        ax.plot(*bindot_args, color='k', marker='.', markersize=3, ls='None',
                zorder=10)

    return fig, ax

def plot_multi_map(all_panel_kws, fig_kws=None, patch_kws=None, mg_kws=None):
    """Make multi-panel map plot.

    Args:
        all_panel_kws (dict): Collection of keyword args for each panel
        fig_kws (dict): Keyword args to pass to plt.figure. Default is None.
        patch_kws (dict): Keyword args to pass to ax.add_patch. Default is
            None.
        mg_kws (dict): MaNGA analysis ID keyword args to pass to
            make_map_title. Default is None.

    Returns:
        tuple: (plt.figure object, plt.figure axis object)
    """
    fig_kws = util.none_to_empty_dict(fig_kws)
    mg_kws = util.none_to_empty_dict(mg_kws)
    patch_kws = util.none_to_empty_dict(patch_kws)

    fig = plt.figure(**fig_kws)
    if 'seaborn' in sys.modules:
        sns.set_context('poster', rc={'lines.linewidth': 2})

    bigax_kws = dict(xlabel='arcsec', ylabel='arcsec',
                     title_kws=dict(fontsize=20), mg_kws=mg_kws)
    bigAxes = make_big_axes(fig, **bigax_kws)

    n_ax = len(all_panel_kws)
    for i, panel_kws in enumerate(all_panel_kws):
        dx = 0.31 * i
        dy = 0.45
        if i >= (n_ax / 2):
            dx = 0.31 * (i - n_ax / 2)
            dy = 0
        left, bottom = (0.08+dx, 0.1+dy)
        if 'seaborn' in sys.modules:
            sns.set_context('talk', rc={'lines.linewidth': 2})
            sns.set_style(rc={'axes.facecolor': '#A8A8A8'})

        ax = fig.add_axes([left, bottom, 0.2, 0.3333])
        panel_kws['cb_kws']['axloc'] = [left + 0.21, bottom, 0.01, 0.3333]

        if 'seaborn' not in sys.modules:
            ax.set_axis_bgcolor('#A8A8A8')
            ax.grid(False, which='both', axis='both')

        fig, ax = plot_map(fig=fig, ax=ax, patch_kws=patch_kws,
                           fig_kws=fig_kws, **panel_kws)

    return fig, ax

def make_plots(columns, values, errors, spaxel_size=0.5, dapdata=None,
               val_no_measure=0, snr_thresh=1, mg_kws=None, titles=None,
               cblabels=None, cmaps=None, make_single=True, make_multi=True,
               make_binnum=False, savefig_single=True, savefig_multi=True,
               savefig_binnum=False, overwrite=False):
    """Make single panel plots and multi-panel plot for set of measurements.

    Args:
       columns (list): Columns of values and errors DataFrames to plot.
       values: Either string that references an attribute from dapdata or an
           array of values.
       errors: Either string that references an attribute from dapdata or an
           array of values.
       spaxel_size (float): Spaxel size in arcsec. Default is 0.5.
       dapdata: dap.DAP object. Default is None.
       val_no_measure (float): Value that corresponds to no measurement.
           Default is 0.
       snr_thresh (float): Signal-to-noise threshold for displaying a bin on a
           map. Default is 1.
       mg_kws (dict): Keyword args with identifying information about the
           galaxy and analysis run. Default is None.
       titles (list): Plot title for each map. Default is None.
       cblabels (list): Colorbar labels. Default is None.
       cmaps (list): Colormaps. Default is None.
       make_single (bool): Make single panel plots. Default is True.
       make_multi (bool): Make multi-panel plot. Default is True.
       make_binnum (bool): Make single panel bin number plot. Default is
           False.
       savefig_single (bool): Save single panel plots. Default is True.
       savefig_multi (bool): Save multi-panel plot. Default is True.
       savefig_binnum (bool): Save single panel bin number plots. Default is
           False.
       overwrite (bool): Overwrite plot if it exists. Default is False.
    """
    # Adjust arguments
    if isinstance(values, str):
        pname_base = copy.deepcopy(values)
        values = util.string_slice_multiindex_df(dapdata, values)

    if isinstance(errors, str):
        errors = util.string_slice_multiindex_df(dapdata, errors)

    mg_kws = util.none_to_empty_dict(mg_kws)
    cmaps = set_cmaps(cmaps, len(columns))

    # Set common plot elements
    xpos = dapdata.drps.xpos.values
    ypos = dapdata.drps.ypos.values
    binid = reorient(dapdata.drps.binid.values)
    delta = spaxel_size / 2.
    extent = set_extent(xpos=xpos, ypos=ypos, delta=delta)
    ax_kws, patch_kws = set_map_background_color(spaxel_size=spaxel_size)

    # Create images
    ims = []
    xys = []
    for col in columns:
        im, xy = make_image(val=values[col].values, err=errors[col].values,
                            xpos=xpos, ypos=ypos, binid=binid, delta=delta,
                            val_no_measure=val_no_measure,
                            snr_thresh=snr_thresh)
        ims.append(im)
        xys.append(xy)

    # Make plots
    all_panel_kws = []
    for i, (im, xy, col, cmap) in enumerate(zip(ims, xys, columns, cmaps)):
        tt, iw, cb = set_map_par(cmap=cmap, title=titles[col],
                                 cblabel=cblabels[col])
        sp_kws = dict(xy_nomeasure=xy, ax_kws=ax_kws, title_kws=tt,
                      fig_kws=dict(figsize=(10, 8)), patch_kws=patch_kws,
                      imshow_kws=iw, cb_kws=cb)
        # Make single panel maps
        if make_single:
            ig = plot_map(im, extent, **sp_kws)
            if savefig_single:
                pname = '_'.join([pname_base, col])
                util.saveplot(name=pname, path_data=dapdata.path_data,
                              plottype='maps', mg_kws=mg_kws, mkdir=True,
                              overwrite=overwrite)

        # Make single panel maps with bin numbers
        if make_binnum:
            binnum_kws = dict(val=values[col].values)
            ig = plot_map(im, extent, dapdata=dapdata, binnum_kws=binnum_kws,
                          **sp_kws)
            if savefig_binnum:
                pname = '_'.join([pname_base, col, 'binnum'])
                util.saveplot(name=pname, path_data=dapdata.path_data,
                              plottype='maps', mg_kws=mg_kws, ext='pdf',
                              mkdir=True, overwrite=overwrite)

        # create dictionaries for multi-panel maps
        t_kws, i_kws, c_kws = set_map_par(cmap=cmap, title=titles[col],
                                          cblabel=cblabels[col],
                                          titlefontsize=20, cbfontsize=16)
        all_panel_kws.append(dict(image=im, extent=extent, xy_nomeasure=xy,
                                  ax_kws=ax_kws, title_kws=t_kws,
                                  imshow_kws=i_kws, cb_kws=c_kws))

    # Make multi-panel maps
    if make_multi:
        ig = plot_multi_map(all_panel_kws=all_panel_kws, patch_kws=patch_kws,
                            mg_kws=mg_kws, fig_kws=dict(figsize=(20, 12)))
        if savefig_multi:
            util.saveplot(name=pname_base, path_data=dapdata.path_data,
                          plottype='maps', mg_kws=mg_kws, mkdir=True,
                          overwrite=overwrite)





# def plot_multi_radial_gradients(map_order,
#                                 args,
#                                 ylabel=r'Flux [10$^{-17}$ erg/s/cm$^2$]',
#                                 n_ax=6,
#                                 leg_kwargs=dict(handlelength=2, loc=1),
#                                 figsize=(20, 12)):
#     """
#     Plot multiple radial gradients at once.
#     """
#     fig = plt.figure(figsize=figsize)
#     if seaborn_installed:
#         sns.set_context('poster', rc={'lines.linewidth': 2})
# 
#     flux_units = 'Flux'
#     if self.dap_mode == 'CUBE':
#         flux_units += ' / spaxel'
#     elif self.dap_mode == 'RSS':
#         flux_units += ' / fiber'
# 
#     bigAxes = fig.add_axes([0.04, 0.05, 0.9, 0.88], frameon=False)
#     bigAxes.set_xticks([])
#     bigAxes.set_yticks([])
#     bigAxes.set_xlabel('R [arcsec]', fontsize=20)
#     bigAxes.set_ylabel('%s [10$^{-17}$ erg/s/cm$^2$]' % flux_units, fontsize=20)
#     bigAxes.set_title(
#         'pid-ifu %s     manga-id %s     %s     %s' % (
#         self.manga_pid, self.manga_id, self.dap_mode, self.exec_num),
#         fontsize=20)
#     
#     bin_edges = np.concatenate((self.binxrl,
#                                np.atleast_1d(self.binyru[-1])))
# 
#     for i, k in enumerate(map_order):
#         dx = 0.31 * i
#         dy = 0.45
#         if i >= (n_ax / 2):
#             dx = 0.31 * (i - n_ax / 2)
#             dy = 0
#         left, bottom = (0.08+dx, 0.1+dy)
#         if seaborn_installed:
#             sns.set_context('poster', rc={'lines.linewidth': 2})
#     
#         ax = fig.add_axes([left, bottom, 0.23, 0.33333])
#         axtitle = args[k]['kwargs']['title_text'].split(' (')[0]
#         ax.set_title(axtitle)
#         
#         if not seaborn_installed:
#             ax.set_axis_bgcolor('#EAEAF2')
#             ax.grid(False, which='both', axis='both')
#     
#         d = args[map_order[i]]
#         p = []
#         lab = []
#         if not np.isnan(d['val']).all():
#             self.plot_radial_gradients(k, args, c_ind=[2, 0],
#                                        leglabels=['F. Belfiore', 'E. Wang'],
#                                        fig=fig, ax=ax)
# 
# 
# 
# def plot_gradients(dapdata, gradname, args,
#                    c_ind=[0], leglabels=None,
#                    leg_kwargs=dict(handlelength=2, loc=1), fig=None, ax=None,
#                    figsize=(10, 8)):
#     """
#     Plot radial gradients.
#     """
#     c = set_spec_line_prop(lw=lw)
#     spec = make_spec_df(dapdata=dapdata, bin=bin, fits_to_plot=fits_to_plot,
#                         rest_frame=rest_frame)
#     if ax is None:
#         fig = plt.figure(figsize=figsize)
#         ax = fig.add_axes([0.17, 0.11, 2/3., 5/6.])
#         ax.set_xlabel('R [arcsec]', fontsize=28)
#         ax.set_ylabel(set_flux_units(dapdata), fontsize=28)
#         for tick in ax.xaxis.get_major_ticks():
#             tick.label.set_fontsize(20)
#         for tick in ax.yaxis.get_major_ticks():
#             tick.label.set_fontsize(20)
#         axtitle = args[gradname]['kwargs']['title_text'].split(' (')[0]
#         ax.set_title(axtitle, fontsize=28)
# 
#     if 'seaborn' not in sys.modules:
#         ax.set_axis_bgcolor('#A8A8A8')
#         ax.grid(False, which='both', axis='both')
#         c = ['b', 'r', 'c']
# 
#     d = args[gradname]
#     p = []
#     lab = []
#     if not np.isnan(d['val']).all():
#         for kk, kkerr, j, author in zip(['val2', 'val'], ['val2_err', 'val_err'],
#                                         c_ind, leglabels):
#             p.append(ax.hlines(args[gradname][kk], self.binxrl, self.binyru,
#                      color=c[j]))
#             #ytmp = np.concatenate((args[gradname][kk],
#             #                      np.atleast_1d(args[gradname][kk][-1])))
#             #p.append(ax.step(bin_edges, ytmp, c=c[j],
#             #         where='post')[0])
#             #p.append(ax.plot(self.binr, args[gradname][kk], c=c[j], zorder=8)[0])
#             ax.plot(self.binr, args[gradname][kk], c=c[j], zorder=8, lw=0.5)
#             ax.scatter(self.binr, args[gradname][kk], facecolor=c[j],
#                        edgecolor='None', s=60, zorder=9)
#             #ax.errorbar(self.binr, args[gradname][kk], yerr=args[gradname][kkerr],
#             #            ecolor=c[j], elinewidth=1, marker='None', ls='None')
#             label = args[gradname]['kwargs']['title_text'].split(' (')[0]
#             lab.append(author)
# 
#     leg = plt.legend(p, lab, **leg_kwargs)
#     plt.setp(leg.get_texts(), fontsize=24)
#     ax.set_xlim(left=0)
#     ax.set_ylim(bottom=0)
#         
#     if not np.isnan(d['val']).all():
#         for kk, kkerr, j in zip(['val2', 'val'], ['val2_err', 'val_err'], [2, 0]):
#             ax.errorbar(self.binr, args[gradname][kk], yerr=args[gradname][kkerr],
#                         ecolor=c[j], elinewidth=1, marker='None', ls='None')






def make_spec_df(dapdata, bin, fits_to_plot, rest_frame):
    """Create DataFrame with spectrum and fits.

    Args:
        dapdata: dap.DAP object.
        bin (int): Bin number. Default is 0.
        rest_frame (bool): If True, show spectrum in rest frame, otherwise show
            spectrum in observed frame. Default is True.
        fits_to_plot (tuple): Extension names of stellar continuum fits or
            stellar continuum plus emission line fits to plot. Default is
            ('smod', 'fullfit_ew', 'fullfit_fb').

    Returns:
        DataFrame
    """
    cols = ['wave', 'flux', 'ivar'] + list(fits_to_plot)
    data = {}
    for col in cols:
        attr = copy.deepcopy(col)
        if rest_frame:
            attr += '_rest'
        if attr is 'wave':
            data[col] = dapdata.__dict__[attr]
        else:
            data[col] = dapdata.__dict__[attr][bin]

    spec = pd.DataFrame(data, columns=cols)

    for col in cols:
        if 'fullfit' in col:
            spec['resid_' + col.split('_')[1]] = spec.flux - spec[col]

    return spec

def set_spec_lims(spec, xlim, ylim):
    """Set axis limits for spectrum.

    Args:
        spec (DataFrame): Spectral data and fits.
        xlim (list): Limits of x-axis.
        ylim (list): Limits of y-axes for spectrum and residuals plots.

    Returns:
        tuple: (limits of x-axis, limits of y-axis, indices where wavelength is
                within x-axis limits)
    """
    if xlim is None:
        xlim = [3600, 9300]
    
    ind = (spec.wave >= xlim[0]) & (spec.wave <= xlim[1])

    if ylim is None:
        p50 = np.percentile(spec.flux[ind], 50)
        ylim_resid_max = 0.1 + p50 * 0.2
        ylim = [[0., p50 * 3.], [-ylim_resid_max, ylim_resid_max]]

    return xlim, ylim, ind

def set_spec_line_prop(lw):
    """Set line properties for spectrum.

    Args:
        lw (int): Linewidth.

    Returns:
        list: Color palette.
    """
    if 'seaborn' in sys.modules:
        c = sns.color_palette('bright', 5)
        sns.set_context('poster', rc={'lines.linewidth': lw})
    else:
        c = ['b', 'lime', 'r', 'DarkOrchid', 'gold']
    return c

def make_stfit_masks(dapdata, bin):
    """Determine regions that were masked in the stellar continuum fit.

    Args:
        dapdata: dap.DAP object.
        bin (int): Bin number.

    Returns:
        tuple: (indices of wavelength for mask edges, )
    """
    # Find mask edges
    ind_split = np.where(np.diff(dapdata.smsk[bin]) != 0)[0]

    # mask values within each region (all the same value)
    smsk_vals = np.array_split(dapdata.smsk[bin], ind_split + 1)
    # mask value for each region
    smsk_val = np.array([item[0] for item in smsk_vals], dtype=int)

    # include starting pixel
    if 0 not in ind_split:
        ind_split = np.append([0], ind_split)
    # include ending pixel
    if dapdata.smsk.shape[1] not in ind_split:
        ind_split = np.append(ind_split, [dapdata.smsk.shape[1] - 1])

    return ind_split, smsk_val

def show_stfit_masks(ax, spec, ind_split, smsk_val, ylim):
    """Display stellar continuum fitting masks.

    Args:
        ax: plt.figure axis object.
        spec (DataFrame): Spectral data and fits.
        ind_split (array): Indices of wavelength mask edges.
        smsk_val (array): Mask value in each region.
        ylim (list): Minimum and maximum values of y-axis.

    Returns:
        plt.figure axis object
    """
    mskkws = dict(facecolor='gray', alpha=0.25)
    lglam = np.log10(spec.wave)
    half_dlglam = (lglam[1] - lglam[0]) / 2.
    for m in range(len(smsk_val)):
        # ind_split give the pixel to the left of the masked region
        # lglam gives the pixel center, hence the half pixel width shift
        x = np.array([10**(lglam[ind_split[m] + 1] - half_dlglam),
                      10**(lglam[ind_split[m+1]] + half_dlglam)])
        if smsk_val[m]:
            ax.fill_between(x, ylim[0], ylim[1], **mskkws)
    return ax

def make_spec_title(dapdata, bin):
    """Make a spectrum title to identify galaxy, analysis run, and bin number.

    Args:
        dapdata: dap.DAP object.
        bin (int): Bin number.

    Returns:
        str: Spectrum title.
    """
    d = dapdata.__dict__
    return '     '.join(('pid-ifu {plateifu}', 'manga-id {mangaid}',
                         '{bintype}-{niter}', 'bin {bin}')).format(bin=bin, **d)

def set_flux_units(dapdata):
    """Determine flux units for spectrum bases on analysis mode.
    Args:
        dapdata: dap.DAP object.

    Returns:
        string: Flux units for spectrum.
    """
    indiv_units = np.array(dapdata.flux_units.split(' ')[1].split('/'))
    indiv_units[indiv_units == 'cm^2'] = 'cm$^2$'
    indiv_units[indiv_units == 'Ang'] = '$\AA$'
    units = '/'.join(indiv_units)
    flux_units = 'Flux / {0} [10$^{{-17}}$ {1}]'.format(indiv_units[-1], units)
    return flux_units

def plot_spectra(dapdata, bins=[0],
                 fits_to_plot=('smod', 'fullfit_fb', 'fullfit_ew'),
                 rest_frame=True, xlim=None, ylim=None, stfit_masks=False, lw=1,
                 figsize=(20, 12), mg_kws=None, savefig=True, overwrite=False):
    for bin in bins:
        ig = plot_spectrum(dapdata, bin=bin, fits_to_plot=fits_to_plot,
                           rest_frame=rest_frame, xlim=xlim, ylim=ylim,
                           stfit_masks=stfit_masks, lw=lw, figsize=figsize,
                           mg_kws=mg_kws, savefig=savefig, overwrite=overwrite)

def plot_spectrum(dapdata, bin=0,
                  fits_to_plot=('smod', 'fullfit_fb', 'fullfit_ew'),
                  rest_frame=True, xlim=None, ylim=None, stfit_masks=False,
                  lw=1, figsize=(20, 12), mg_kws=None, savefig=True,
                  overwrite=False):
    """Plot spectrum and residuals.

    Args:
        dapdata: dap.DAP object.
        bin (int): Bin number. Default is 0.
        rest_frame (bool): If True, show spectrum in rest frame, otherwise show
            spectrum in observed frame. Default is True.
        fits_to_plot (tuple): Extension names of stellar continuum fits or
            stellar continuum plus emission line fits to plot. Default is
            ('smod', 'fullfit_ew', 'fullfit_fb').
        xlim (tuple): Limits of x-axis. Default is None.
        ylim (tuple): Limits of y-axes for spectrum and residuals plots. Default
            is None.
        stfit_masks (bool): Show masks used in stellar continuum fitting.
            Default is False.
        figsize (tuple): Figure size in inches. Default is (20, 12).
        mg_kws (dict): Keyword args with identifying information about the
           galaxy and analysis run. Default is None.
        savefig (bool): Save plot. Default is True.
        overwrite (bool): Overwrite plot if it exists. Default is False.

    Returns:
        plt.figure object
    """
    colors = set_spec_line_prop(lw=lw)
    spec = make_spec_df(dapdata=dapdata, bin=bin, fits_to_plot=fits_to_plot,
                        rest_frame=rest_frame)
    xlim, ylim, ind = set_spec_lims(spec, xlim, ylim)

    if stfit_masks:
        ind_split, smsk_val = make_stfit_masks(dapdata, bin)

    fig = plt.figure(figsize=figsize)

    panels = ['spectrum', 'residuals']
    axpar = pd.DataFrame(np.array([[0.32, 0.63], [0.1, 0.2]]), index=panels,
                         columns=['bottom', 'height'])
    fit_ids = [k.split('fullfit_')[1] for k in spec.columns if 'fullfit' in k]
    for i, panel in enumerate(panels):
        ax = fig.add_axes([0.1, axpar.bottom[panel], 0.85, axpar.height[panel]])

        if xlim is not None:
            ax.set_xlim(xlim)
        if ylim is not None:
            ax.set_ylim(ylim[i])

        p = []
        labels = []
        if panel is 'spectrum':
            ax.set_title(make_spec_title(dapdata, bin))
            plt.setp(ax.get_xticklabels(), visible=False)
            ax.set_ylabel(set_flux_units(dapdata))

            # spectrum
            p.append(ax.plot(spec.wave[ind], spec.flux[ind],
                     color='#808080')[0])
            labels.append('galaxy')

            # stellar continuum fit
            p.append(ax.plot(spec.wave[spec.smod > 0.],
                             spec.smod[spec.smod > 0.], color=colors[1])[0])
            labels.append('stellar continuum fit')

            emfits = [spec['fullfit_' + fit_id] for fit_id in fit_ids]
            em_labels = [fit_id.upper() + ' emline + cont fit'
                         for fit_id in fit_ids]

        elif panel is 'residuals':
            ax.set_xlabel(r'$\lambda_{\rm rest} [\AA]$')
            ax.set_ylabel(r'$\Delta$')
            ax.plot([xlim[0], xlim[1]], [0, 0], color='gray', alpha=0.75, 
                    lw=2, zorder=1)
            emfits = [spec['resid_' + fit_id] for fit_id in fit_ids]

        # emission line + stellar continuum fits
        for j, (emfit, em_label) in enumerate(zip(emfits, em_labels)):
            p.append(ax.plot(spec.wave[spec.smod > 0.], emfit[spec.smod > 0.],
                     color=colors[2-2*j])[0])
            labels.append(em_label)

        if panel is 'spectrum':
            plt.legend(p, labels, loc=2)

        if stfit_masks:
            show_stfit_masks(ax, spec, ind_split, smsk_val, ylim[i])

    if savefig:
        mg_kws['bin'] = bin
        util.saveplot(name='spec', path_data=dapdata.path_data, plottype='spec',
                      mg_kws=mg_kws, mkdir=True,overwrite=overwrite)

    return fig








def plot_emlines(dapdata, bin=0, mg_kws=None):
    """
    """
    win_cen = np.array([3727., 4861., 4985., 6565., 6565., 6723.])

    plot_emline_multi(dapdata=dapdata, bin=bin, mg_kws=mg_kws)
    for i in range(6):
        nii = False
        if i == 3:
            nii = True
        plot_emline(dapdata=dapdata, bin=bin, mg_kws=mg_kws, nii=nii,
                    win_cen=win_cen[i])




def plot_emline_multi(dapdata, bin=0, mg_kws=None,
                      kwargs={'alpha':0.75, 'lw':2}, figsize=(20, 12)):
    """Plot multiple panel zoom-ins of spectra near strong emission lines.

    Args:
        bin (int): bin number
        kwargs (dict): keyword args for ax.plot
        figsize (tuple): figure width and height in inches
    
    """
    win_cen = np.array([3727., 4861., 4985., 6565., 6565., 6723.])
    
    fig = plt.figure(figsize=figsize)
    bigax_kws = dict(xlabel=r'$\lambda \, [\AA]$',
                     ylabel=set_flux_units(dapdata), mg_kws=mg_kws,
                     title_kws=dict(label=make_spec_title(dapdata, bin),
                                    fontsize=20))
    bigAxes = make_big_axes(fig, axloc=(.04, 0.06, 0.9, 0.9), **bigax_kws)
    fig.subplots_adjust(wspace=0.2, hspace=0.15, left=0.1, bottom=0.1,
                         right=0.95, top=0.95)
    for i in range(6):
        ax = fig.add_subplot(2, 3, i+1)
        nii = False
        if i == 3:
            nii = True
        plot_emline(dapdata=dapdata, fig=fig, ax=ax, bin=bin, xlim=None,
                    ylim=None, win_cen=win_cen[i], nii=nii,
                    kwargs=kwargs, figsize=figsize)





def plot_emline(dapdata, fig=None, ax=None, bin=0, xlim=None, ylim=None,
                win_cen=None, nii=False, mg_kws=None,
                kwargs={'alpha':0.75, 'lw':2}, figsize=(10, 8)):
    """Plot data and model spectra near strong emission lines.

    Args:
        fig: figure object
        ax: axis object
        bin (int): bin number
        xlim (list): minimum and maximum x-axis values
        ylim (list): minimum and maximum y-axis values
        nii (bool): If True, maximum y value is determined by the
            [NII]6548,6583 lines (not Halpha)
        kwargs (dict): keyword args for ax.plot
        figsize (tuple): figure width and height in inches
    
    """
    colors = set_spec_line_prop(lw=kwargs['lw'])

    if ax is None:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_axes([0.17, 0.15, 0.72, 0.75])
        ax.set_xlabel(r'$\lambda \, [\AA]$', fontsize=24)
        ax.set_ylabel(set_flux_units(dapdata), fontsize=24)
        x_offs = [-9, -3, -9, -9, -12, -3, -5, -12, -5]
        indiv_ax = True
    else:
        x_offs = [-12, -4, -12, -12, -20, -4, -12, -20, -6]
        indiv_ax = False

    wave = dapdata.wave_rest[bin]
    flux = dapdata.flux_rest[bin]
    ivar = dapdata.ivar_rest[bin]
    noise = ivar**-0.5
    stmodel = dapdata.smod_rest[bin]
    fullfit_ew = dapdata.fullfit_ew_rest[bin]
    fullfit_fb = dapdata.fullfit_fb_rest[bin]

    names = ['[OII]3727', r'H$\beta$', '[OIII]4959', '[OIII]5007',
    '[NII]6548', r'H$\alpha$', '[NII]6583', '[SII]6716', '[SII]6731']

    if xlim is None:
        xmin = win_cen - 50.
        xmax = win_cen + 50.
    else:
        xmin, xmax = xlim
    if ylim is None:
        ind_w = np.where((wave > xmin) & (wave < xmax))
        ymin_tmp = flux[ind_w] - noise[ind_w]
        ymin_tmp = np.append(ymin_tmp, 1e6)
        ymax_tmp = flux[ind_w] + noise[ind_w]
        ymax_tmp = np.append(ymax_tmp, 1e-10)
        ymin = np.min(ymin_tmp[np.isfinite(ymin_tmp)])
        ymax = np.max(ymax_tmp[np.isfinite(ymax_tmp)])
        if nii:
            ind_w_nii = np.where((wave > 6575.) & (wave < 6591.))
            ymax_tmp = flux[ind_w_nii] + noise[ind_w_nii]
            ymax_tmp = np.append(ymax_tmp, 1e-10)
            ymax = np.max(ymax_tmp[np.isfinite(ymax_tmp)])
        dy = ymax - ymin
        ymin = ymin - (dy * 0.2)
        ymax = ymax + (dy * 0.2)
    else:
        ymin, ymax = ylim
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

    ind = np.where((wave > xmin-5.) & (wave < xmax+5.))[0]
    #ax.plot(wave[ind], flux[ind], 'gray', drawstyle='steps-mid', zorder=1)
    ax.scatter(wave[ind], flux[ind], c='gray', edgecolor='None', zorder=1)
    ax.errorbar(wave[ind], flux[ind], yerr=noise[ind], ecolor='gray',
                fmt='none', zorder=1)
    pst = ax.plot(wave[ind], stmodel[ind], color=colors[1], **kwargs)[0]
    pFB = ax.plot(wave[ind], fullfit_fb[ind], color=colors[2], **kwargs)[0]
    pEW = ax.plot(wave[ind], fullfit_ew[ind], color=colors[0], **kwargs)[0]

    for name, w, x_off in zip(names, dapdata.elopar.restwave, x_offs):
        if (w > xmin) and (w < xmax):
            if (not nii) or (name is not r'H$\alpha$'):
                ax.text(w + x_off, ymax-dy*0.08, name, color='k', fontsize=20)
                ax.plot([w, w], [ymax-dy*0.16, ymax-dy*0.11], color='k')

    leg = plt.legend([pst, pFB, pEW],
                     ['Stellar Cont.', 'Belfiore', 'Wang'],
                     borderaxespad=0.1, handlelength=1.25,
                     handletextpad=0.5, ncol=3, loc=8)
    if indiv_ax:
        plt.setp(leg.get_texts(), fontsize=24)
    else:
        plt.setp(leg.get_texts(), fontsize=16)



