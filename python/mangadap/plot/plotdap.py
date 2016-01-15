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
from matplotlib.colors import LogNorm

from astropy.stats import sigma_clip

from mangadap.plot import util

import warnings
try:
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        import seaborn as sns
except ImportError:
    print('Seaborn could not be imported. Continuing...')

def _set_extent(xpos, ypos, delta):
    """Set extent of map."""
    return np.array([-xpos.max() - delta, -xpos.min() + delta,
                     ypos.min() - delta, ypos.max() + delta])

def _reorient(x):
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
               snr_thresh, log_colorbar):
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
        log_colorbar (bool): If True, use log colorbar.

    Returns:
        tuple: (masked array of image,
                tuple of (x, y) coordinates of bins with no measurement)
    """
    # create a masked array of the data
    im = val[binid]
    im[binid == -1] = np.nan

    im_err = err[binid]
    im_err[binid == -1] = np.nan

    no_measure = _make_mask_no_measurement(im, im_err, val_no_measure,
                                           snr_thresh, log_colorbar)
    no_data = _make_mask_no_data(im, no_measure)
    image = np.ma.array(im, mask=no_data)

    # spaxels with data but no measurement
    xpos_re = _reorient(xpos)
    ypos_re = _reorient(ypos)
    xy_nomeasure = (-(xpos_re[no_measure]+delta), ypos_re[no_measure]-delta)
    return image, xy_nomeasure

def _make_mask_no_measurement(data, err, val_no_measure, snr_thresh,
                              log_colorbar):
    """Mask invalid measurements within a data array.

    Args:
        data (array): Data.
        err (array): Error.
        val_no_measure (float): Value in data array that corresponds to no
            measurement.
        snr_thresh (float): Signal-to-noise threshold for keeping a valid
            measurement.
        log_colorbar (bool): If True, use log colorbar.

    Returns:
        array: Boolean array for mask (i.e., True corresponds to value to be
            masked out).
    """
    no_measure = (data == val_no_measure)
    if np.all([[i in [None, np.nan] for i in x] for x in err]):
        err = None
    if err is not None:
        no_measure[err == 0.] = True
        no_measure[(np.abs(data / err) < snr_thresh)] = True
        # NEW for LogNorm
        if log_colorbar:
            no_measure[data <= 0.] = True

    return no_measure

def _make_mask_no_data(data, mask_no_measurement):
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

def _log_colorbar_ticks(cbrange):
    """Set ticks and ticklabels for a log normalized colorbar.

    Args:
        cbrange (list): Colorbar range.

    Returns:
        array
    """
    subs = [1., 2., 3., 6.]
    bottom = np.floor(np.log10(cbrange[0]))
    top = np.ceil(np.log10(cbrange[1]))
    decs = np.arange(bottom, top+1)
    tmp = np.array([sub * 10.**dec for dec in decs for sub in subs])
    return tmp[np.logical_and((tmp >= cbrange[0]), (tmp <= cbrange[1]))]
    
def _log_tick_format(value):
    """Format tick labels for log axis.

    If value between ___, return as ___:
       [0, 999], int
       [0.1, 0.99], 1 digit float
       otherwise: exponential notation

    Args:
        value (float): Input value.

    Returns:
        str
    """
    exp = np.floor(np.log10(value))
    base = value / 10**exp
    if exp in [0, 1, 2]:   
        return '{0:d}'.format(int(value))
    elif exp == -1:
        return '{0:.1f}'.format(value)
    else:
        return '{0:d}e{1:d}'.format(int(base), int(exp))

def _set_vmin_vmax(d, cbrange):
    """Set minimum and maximum values of the color map."""
    if 'vmin' not in d.keys():
        d['vmin'] = cbrange[0]
    if 'vmax' not in d.keys():
        d['vmax'] = cbrange[1]
    return d

def _cbrange_sigclip(image, sigma):
    """Sigma clip colorbar range.

    Args:
        image (masked array): Image.
        sigma (float): Sigma to clip.

    Returns:
        list: Colorbar range.
    """
    try:
        imclip = sigma_clip(image.data[~image.mask], sigma=sigma)
    except TypeError:
        imclip = sigma_clip(image.data[~image.mask], sig=sigma)

    try:
        cbrange = [imclip.min(), imclip.max()]
    except ValueError:
        cbrange = [image.min(), image.max()]

    return cbrange

def _cbrange_percentile_clip(image, lower, upper):
    """Clip colorbar range according to percentiles.

    Args:
        image (masked array): Image.
        lower (float): Lower percentile boundary.
        upper (float): Upper percentile boundary.

    Returns:
        list: Colorbar range.
    """
    cblow = np.percentile(image.data[~image.mask], lower)
    cbup = np.percentile(image.data[~image.mask], upper)
    return [cblow, cbup]

def _cbrange_user_defined(cbrange, cbrange_user):
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

def _set_cbrange(image, column, cb_kws):
    """Set colorbar range.

    Args:
        image (masked array): Image.
        column (str): Column name.
        cbrange (list): Colorbar kwargs.

    Returns:
        list: Colorbar range.
    """
    if cb_kws.get('sigclip') is not None:
        cbr = _cbrange_sigclip(image, cb_kws['sigclip'])
    elif cb_kws.get('percentile_clip', None) is not None:
        try:
            cbr = _cbrange_percentile_clip(image, *cb_kws['percentile_clip'])
        except IndexError:
            cbr = [0.1, 1]
    else:
        cbr = [image.min(), image.max()]
    
    if cb_kws.get('cbrange') is not None:
        cbr = _cbrange_user_defined(cbr, cb_kws['cbrange'])

    if cb_kws['symmetric']:
        cb_max = np.max(np.abs(cbr))
        cbr = [-cb_max, cb_max]

    cbr, cb_kws['ticks'] = _set_cbticks(cbr, cb_kws)

    if cb_kws.get('log_colorbar', False):
        try:
            im_min = np.min(image[image > 0.])
        except ValueError:
            im_min = 0.1
        if im_min is np.ma.masked:
            im_min = 0.1
        cbr[0] = np.max((cbr[0], im_min))

    cb_kws['cbrange'] = cbr

    return cb_kws

def _set_cbticks(cbrange, cb_kws):
    """Set colorbar ticks.

    Adjust colorbar range if using a discrete colorbar so that the ticks fall in
    the middle of each level.

    Args:
        cbrange (list): Colorbar range.
        cb_kws (dict): Keyword args to set and draw colorbar.

    Return:
        tuple: colorbar range, colorbar tick numbers
    """
    if cb_kws.get('log_colorbar'):
        ticks = _log_colorbar_ticks(cbrange)
    else:
        try:
            ticks = MaxNLocator(cb_kws['n_ticks']).tick_values(*cbrange)
        except AttributeError:
            print('AttributeError: MaxNLocator instance has no attribute'
                  ' "tick_values" ')

    # if discrete colorbar, offset upper and lower cbrange so that ticks are in
    # the center of each level
    if cb_kws['n_levels'] is not None:
        offset = (ticks[1] - ticks[0]) / 2.
        cbrange = [ticks[0] - offset, ticks[-1] + offset]
        if cb_kws.get('tick_everyother', False):
            ticks = ticks[::2]

    return cbrange, ticks

def draw_colorbar(fig, mappable, column=None, axloc=None, cbrange=None,
                  ticks=None, label_kws=None, tick_params_kws=None,
                  log_colorbar=False, **extras):
    """Make colorbar.

    Args:
        fig: plt.figure object.
        mappable: Plotting element to map to colorbar.
        column (str): Column name. Default is None.
        axloc (list): Specify (left, bottom, width, height) of colorbar axis.
            Defaults to None.
        cbrange (list): Colorbar min and max.
        ticks (list): Ticks on colorbar.
        label_kws (dict): Keyword args to set colorbar label. Default is None.
        tick_params_kws (dict): Keyword args to set colorbar tick parameters.
            Default is None.

    Returns:
        tuple: (plt.figure object, plt.figure axis object)
    """
    label_kws = util.none_to_empty_dict(label_kws)
    tick_params_kws = util.none_to_empty_dict(tick_params_kws)

    cax = (fig.add_axes(axloc) if axloc is not None else None)
    try:
        cb = fig.colorbar(mappable, cax, ticks=ticks)
    except ValueError:
        cb = None
    else:
        cb.ax.tick_params(**tick_params_kws)
        if label_kws.get('label') is not None:
            cb.set_label(**label_kws)
        if log_colorbar:
            cb.set_ticklabels([_log_tick_format(tick) for tick in ticks])

    return fig, cb

def _set_cmap(cm_name, n_levels=None):
    """Set the colormaps.

    Args:
        cm_name (str): Name of colormap.
        n_levels (int): Number of discrete levels of colormap. If None, then
            produce continuous colormap. Default is None.

    Returns:
        colormap
    """
    cmap = _string_to_cmap(cm_name)

    if n_levels is not None:
        cmap = util.cmap_discretize(cmap, n_levels)

    return cmap

def _string_to_cmap(cm_name):
    """Return colormap given name.

    Args:
        cm_name (str): Name of colormap.

    Returns:
        colormap
    """
    if cm_name is None:
        cm_name = 'linear_Lab'

    if 'linear_Lab' in cm_name:
        try:
            cmap, cmap_r = util.linear_Lab()
        except IOError:
            cmap = cm.Blues_r
        else:
            if '_r' in cm_name:
                cmap = cmap_r
    else:
        try:
            cmap = cm.__dict__[cm_name]
        except KeyError:
            cmap = cm.Blues_r
    return cmap

def _set_map_background_color(spaxel_size, color='#A8A8A8'):
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


def _set_title_kws(column, titles, fontsize=28):
    """Set title keyword args.

    Args:
        column (str): Column name.
        titles (str): Plot titles.
        fontsize (int): Title font size. Default is 28.

    Returns:
        dict
    """
    return dict(fontsize=fontsize, label=titles[column])

def _set_imshow_kws(column, cb_kws):
    """Set imshow keyword args.

    Args:
        column (str): Column name.
        cb_kws (dict): Color bar keyword args.
        
    Returns:
        dict
    """
    imshow_kws = dict(cmap=cb_kws['cmap'])
    if cb_kws.get('log_colorbar', False):
        imshow_kws['norm'] = LogNorm()

    return imshow_kws

def _set_cb_kws(column, cb_kws_master):
    """Set colorbar keyword args.

    Args:
        column (str): Column name.
        cb_kws_master (dict): Color bar keyword args.

    Returns:
        dict
    """
    cb_kws_default = dict(axloc=[0.82, 0.1, 0.02, 5/6.], cbrange=None,
                          symmetric=False, cmap=None, n_levels=None,
                          label_kws=dict(size=24),
                          tick_params_kws=dict(labelsize=24))

    # Load master kwargs
    cb_kws = {}
    for k, v in cb_kws_master.items():
        try:
            cb_kws[k] = v[column]
        except (TypeError, IndexError):
            cb_kws[k] = v

    # Load default kwargs
    for k, v in cb_kws_default.items():
        if k not in cb_kws:
            cb_kws[k] = v

    if 'label' in cb_kws:
        cb_kws['label_kws']['label'] = cb_kws.pop('label')

    cb_kws['cmap'] = _set_cmap(cb_kws['cmap'], cb_kws['n_levels'])

    return cb_kws

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

def _make_map_title(analysis_id):
    """Make a map title to identify galaxy and analysis run.

    Args:
        analysis_id (dict): Identifying info about galaxy and analysis run.

    Returns:
        str: Map title.
    """
    return '     '.join(('pid-ifu {plate}-{ifudesign}', 'manga-id {mangaid}',
                         '{bintype}-{niter}')).format(**analysis_id)

def _make_big_axes(fig, axloc=(0.04, 0.05, 0.9, 0.88), xlabel=None, ylabel=None,
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
            _make_map_title. Default is None.

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
        title_kws['label'] = _make_map_title(mg_kws)
    if 'label' in title_kws:
        bigAxes.set_title(**title_kws)
    return bigAxes

def _multi_panel_loc(panelnum, n_panels):
    """
    Args:
        panelnum: Panel number (0-indexed).
        n_panels: Number of panels.

    Returns:
        tuple: Left and bottom locations of panel.
    """
    dx = 0.31 * panelnum
    dy = 0.45
    if panelnum >= (n_panels / 2):
        dx = 0.31 * (panelnum - n_panels / 2)
        dy = 0
    left, bottom = (0.08 + dx, 0.1 + dy)
    return (left, bottom)

def _show_bin_num(dapdata, val, ax, imshow_kws, fontsize=6):
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
        fontsize_tmp = _set_bin_num_fontsize(fontsize, i, nb)
        color = _set_bin_num_color(v, imshow_kws)
        ax.text(-x, y, str(i), fontsize=fontsize_tmp, color=color,
                horizontalalignment='center', verticalalignment='center',
                zorder=10)
    return ax

def _set_bin_num_color(value, imshow_kws):
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

def _set_bin_num_fontsize(fontsize, number, nbin):
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

def plot_map(image, extent, xy_nomeasure=None, column=None, fig=None, ax=None,
             dapdata=None, fig_kws=None, ax_kws=None, title_kws=None,
             patch_kws=None, imshow_kws=None, cb_kws=None, binnum_kws=None,
             bindot_args=()):
    """Make single panel map or one panel of multi-panel map plot.

    Args:
        image (masked array): Image to display.
        extent (array): Minimum and maximum x- and y-values.
        xy_nomeasure (tuple): x- and y-coordinates of spaxels without
            measurements.
        column (str): Column name. Default is None.
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
        binnum_kws (dict): Keyword args to pass to _show_bin_num. Default is
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

    cb_kws = _set_cbrange(image, column, cb_kws)
    imshow_kws = _set_vmin_vmax(imshow_kws, cb_kws['cbrange'])

    # Plot regions with no measurement as hatched
    if xy_nomeasure is not None:
        for xh, yh in zip(*xy_nomeasure):
            ax.add_patch(mpl.patches.Rectangle((xh, yh), **patch_kws))

    p = ax.imshow(image, interpolation='none', extent=extent, **imshow_kws)

    fig, cb = draw_colorbar(fig, p, column, **cb_kws)

    if binnum_kws:
        ax = _show_bin_num(dapdata=dapdata, ax=ax, imshow_kws=imshow_kws,
                           **binnum_kws)

    if bindot_args:
        ax.plot(*bindot_args, color='k', marker='.', markersize=3, ls='None',
                zorder=10)

    if 'seaborn' in sys.modules:
        sns.set_style(rc={'axes.facecolor': '#EAEAF2'})

    return fig, ax

def plot_map_multi(all_panel_kws, fig_kws=None, patch_kws=None, mg_kws=None):
    """Make multi-panel map plot.

    Args:
        all_panel_kws (dict): Collection of keyword args for each panel
        fig_kws (dict): Keyword args to pass to plt.figure. Default is None.
        patch_kws (dict): Keyword args to pass to ax.add_patch. Default is
            None.
        mg_kws (dict): MaNGA analysis ID keyword args to pass to
            _make_map_title. Default is None.

    Returns:
        tuple: plt.figure object
    """
    fig_kws = util.none_to_empty_dict(fig_kws)
    mg_kws = util.none_to_empty_dict(mg_kws)
    patch_kws = util.none_to_empty_dict(patch_kws)

    fig = plt.figure(**fig_kws)
    if 'seaborn' in sys.modules:
        sns.set_context('poster', rc={'lines.linewidth': 2})

    bigax_kws = dict(xlabel='arcsec', ylabel='arcsec',
                     title_kws=dict(fontsize=20), mg_kws=mg_kws)
    bigAxes = _make_big_axes(fig, **bigax_kws)

    n_ax = len(all_panel_kws)
    for i, panel_kws in enumerate(all_panel_kws):
        left, bottom = _multi_panel_loc(panelnum=i, n_panels=n_ax)
        if 'seaborn' in sys.modules:
            sns.set_context('talk', rc={'lines.linewidth': 2})
            sns.set_style(rc={'axes.facecolor': '#A8A8A8'})

        ax = fig.add_axes([left, bottom, 0.2, 0.3333])
        panel_kws['cb_kws']['axloc'] = [left + 0.21, bottom, 0.01, 0.3333]

        if 'seaborn' not in sys.modules:
            ax.set_axis_bgcolor('#A8A8A8')
            ax.grid(False, which='both', axis='both')

        ig, ax = plot_map(fig=fig, ax=ax, patch_kws=patch_kws, fig_kws=fig_kws,
                          **panel_kws)

    return fig

def plot_maps(columns, values, errors, spaxel_size=0.5, dapdata=None,
              val_no_measure=0, snr_thresh=1, mg_kws=None, titles=None,
              cb_kws_master=None, main=True, make_single=True, make_multi=True,
              make_binnum=None, savefig_single=True, savefig_multi=True,
              savefig_binnum=None, savedir=None, overwrite=False):
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
        cb_kws_master (dict): Color bar keyword args. Default is None.
        main (bool): True is running as script. False is running interactively.
            Default is True.
        make_single (bool): Make single panel plots. Default is True.
        make_multi (bool): Make multi-panel plot. Default is True.
        make_binnum (list): Make single panel bin number plot. Default is
            None.
        savefig_single (bool): Save single panel plots. Default is True.
        savefig_multi (bool): Save multi-panel plot. Default is True.
        savefig_binnum (list): Save single panel bin number plots. Default is
            None.
        savedir (str): Directory to save plots in. If None, then default to data
            directory. Default is None.
        overwrite (bool): Overwrite plot if it exists. Default is False.
    """
    # Adjust arguments
    if isinstance(values, str):
        pname_base = copy.deepcopy(values)
        values = util.string_slice_multiindex_df(dapdata, values)

    if isinstance(errors, str):
        errors = util.string_slice_multiindex_df(dapdata, errors)

    mg_kws = util.none_to_empty_dict(mg_kws)
    cb_kws_master = util.none_to_empty_dict(cb_kws_master)
    if make_binnum is None:
        make_binnum = [False for _ in columns]
    if savefig_binnum is None:
        savefig_binnum = [False for _ in columns]

    # Set common plot elements
    xpos = dapdata.drps.xpos.values
    ypos = dapdata.drps.ypos.values
    binid = _reorient(dapdata.drps.binid.values)
    delta = spaxel_size / 2.
    extent = _set_extent(xpos=xpos, ypos=ypos, delta=delta)
    ax_kws, patch_kws = _set_map_background_color(spaxel_size=spaxel_size)

    # Create images
    ims = []
    xys = []
    for col in columns:
        log_cbar = cb_kws_master.get('log_colorbar', {col: False})[col]
        im, xy = make_image(val=values[col].values, err=errors[col].values,
                            xpos=xpos, ypos=ypos, binid=binid, delta=delta,
                            val_no_measure=val_no_measure,
                            snr_thresh=snr_thresh, log_colorbar=log_cbar)
        ims.append(im)
        xys.append(xy)

    # Make plots
    all_panel_kws = []
    for i, (im, xy, col) in enumerate(zip(ims, xys, columns)):
        tt = _set_title_kws(column=col, titles=titles)
        cb = _set_cb_kws(column=col, cb_kws_master=cb_kws_master)
        iw = _set_imshow_kws(column=col, cb_kws=cb)
        sp_kws = dict(xy_nomeasure=xy, column=col, ax_kws=ax_kws, title_kws=tt,
                      fig_kws=dict(figsize=(10, 8)), patch_kws=patch_kws,
                      imshow_kws=iw, cb_kws=cb)
        # Make single panel maps
        if make_single:
            fig, ax = plot_map(im, extent, **sp_kws)

            if savefig_single:
                pname = '_'.join([pname_base, col])
                util.saveplot(name=pname, path_data=dapdata.path_data,
                              category='maps', mg_kws=mg_kws, main=main,
                              mkdir=True, savedir=savedir, overwrite=overwrite)
            if main:
                plt.close(fig)

        # Make single panel maps with bin numbers
        if make_binnum[i]:
            binnum_kws = dict(val=values[col].values)
            fig, ax = plot_map(im, extent, dapdata=dapdata,
                               binnum_kws=binnum_kws, **sp_kws)
            if savefig_binnum[i]:
                pname = '_'.join([pname_base, col, 'binnum'])
                util.saveplot(name=pname, path_data=dapdata.path_data,
                              category='maps', mg_kws=mg_kws, ext='pdf',
                              main=main, mkdir=True, savedir=savedir,
                              overwrite=overwrite)
            if main:
                plt.close(fig)

        # create dictionaries for multi-panel maps
        t_kws = _set_title_kws(column=col, titles=titles, fontsize=20)
        c_kws = _set_cb_kws(column=col, cb_kws_master=cb_kws_master)
        i_kws = _set_imshow_kws(column=col, cb_kws=c_kws)
        all_panel_kws.append(dict(image=im, extent=extent, xy_nomeasure=xy,
                                  column=col, ax_kws=ax_kws, title_kws=t_kws,
                                  imshow_kws=i_kws, cb_kws=c_kws))

    # Make multi-panel maps
    if make_multi:
        fig = plot_map_multi(all_panel_kws=all_panel_kws, patch_kws=patch_kws,
                             mg_kws=mg_kws, fig_kws=dict(figsize=(20, 12)))
        if savefig_multi:
            util.saveplot(name=pname_base, path_data=dapdata.path_data,
                          category='maps', mg_kws=mg_kws, main=main, mkdir=True,
                          savedir=savedir, overwrite=overwrite)
        if main:
            plt.close(fig)


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

def _set_spec_lims(spec, xlim, ylim):
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

def _set_spec_line_prop(lw):
    """Set line properties for spectrum.

    Args:
        lw (int): Linewidth.

    Returns:
        list: Color palette.
    """
    if 'seaborn' in sys.modules:
        c_in = sns.color_palette('bright', 5)
        c = [c_in[i] for i in (1, 2, 0, 3, 4)]
        sns.set_context('poster', rc={'lines.linewidth': lw})
    else:
        c = ['b', 'lime', 'r', 'DarkOrchid', 'gold']
    return c

def _make_stfit_masks(dapdata, bin):
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

def _show_stfit_masks(ax, spec, ind_split, smsk_val, ylim):
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

def _make_spec_title(dapdata, bin):
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

def _set_flux_units(dapdata):
    """Determine flux units for spectrum bases on analysis mode.
    Args:
        dapdata: dap.DAP object.

    Returns:
        string: Flux units for spectrum.
    """
    indiv_units = np.array(dapdata.flux_units.split(' ')[1].split('/'))
    indiv_units = indiv_units.astype('<U32')
    indiv_units[indiv_units == 'cm^2'] = 'cm$^2$'
    indiv_units[indiv_units == 'Ang'] = '$\AA$'
    indiv_units[indiv_units == 'arcsec^2'] = 'arcsec$^2$'
    units = '/'.join(indiv_units)
    flux_units = 'Flux / {0} [10$^{{-17}}$ {1}]'.format(indiv_units[-1], units)
    return flux_units


def convert_bins_from_string_to_list(dapdata, bins):
    """If bins is specified by a string, convert it to an array.

    Args:
        dapdata: dap.DAP object.
        bins: If bins = 'all' or bins = *N*total (e.g., 100total) in config
            file instead of a list of integers, then create an array of
            integers.

    Returns:
        array (or list of intgers)
    """
    try:
        if bins.lower() == 'all':
            bins = np.arange(len(dapdata.bins))
        elif 'total' in bins:
            n_spec = int(bins.split('total')[0])
            if n_spec >= len(dapdata.bins):
                bins = np.arange(len(dapdata.bins))
            else:
                bins = np.linspace(0, len(dapdata.bins) - 1, n_spec, dtype=int)
    except AttributeError:
        pass
    return bins

def plot_spectra(dapdata, bins=(0),
                 fits_to_plot=('smod', 'fullfit_fb', 'fullfit_ew'),
                 rest_frame=True, xlim=None, ylim=None, stfit_masks=False, lw=1,
                 figsize=(20, 12), mg_kws=None, main=True, savefig=True,
                 overwrite=False):
    """Plot multiple spectra.

    Show all bins by setting bins = all in config file.
    Show N bins by setting bins = *N*total (e.g., 100total) in config file.
    """
    bins = convert_bins_from_string_to_list(dapdata, bins)
    for bin in bins:
        fig = plot_spectrum(dapdata, bin=bin, fits_to_plot=fits_to_plot,
                            rest_frame=rest_frame, xlim=xlim, ylim=ylim,
                            stfit_masks=stfit_masks, lw=lw, figsize=figsize,
                            mg_kws=mg_kws, main=main, savefig=savefig,
                            overwrite=overwrite)
        if main:
            plt.close(fig)

def plot_spectrum(dapdata, bin=0,
                  fits_to_plot=('smod', 'fullfit_fb', 'fullfit_ew'),
                  rest_frame=True, xlim=None, ylim=None, stfit_masks=False,
                  lw=1, figsize=(20, 12), mg_kws=None, main=True, savefig=True,
                  savedir=None, overwrite=False):
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
        main (bool): True is running as script. False is running interactively.
            Default is True.
        savefig (bool): Save plot. Default is True.
        savedir (str): Directory to save plots in. If None, then default to data
            directory. Default is None.
        overwrite (bool): Overwrite plot if it exists. Default is False.

    Returns:
        plt.figure object
    """
    colors = _set_spec_line_prop(lw=lw)
    spec = make_spec_df(dapdata=dapdata, bin=bin, fits_to_plot=fits_to_plot,
                        rest_frame=rest_frame)
    xlim, ylim, ind = _set_spec_lims(spec, xlim, ylim)

    if stfit_masks:
        ind_split, smsk_val = _make_stfit_masks(dapdata, bin)

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
            ax.set_title(_make_spec_title(dapdata, bin))
            plt.setp(ax.get_xticklabels(), visible=False)
            ax.set_ylabel(_set_flux_units(dapdata))

            # spectrum
            p.append(ax.plot(spec.wave[ind], spec.flux[ind],
                             color='#808080')[0])
            labels.append('galaxy')

            # stellar continuum fit
            p.append(ax.plot(spec.wave[spec.smod > 0.],
                             spec.smod[spec.smod > 0.], color=colors[0])[0])
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
                             color=colors[j+1])[0])
            labels.append(em_label)

        if panel is 'spectrum':
            plt.legend(p, labels, loc=2)

        if stfit_masks:
            _show_stfit_masks(ax, spec, ind_split, smsk_val, ylim[i])

    if savefig:
        mg_kws['bin'] = bin
        util.saveplot(name='spec', path_data=dapdata.path_data,
                      category='spectra', mg_kws=mg_kws, main=main, mkdir=True,
                      savedir=savedir, overwrite=overwrite)

    return fig


def plot_emline_spectra(dapdata, bins=(0), pnames=None, win_cen=None,
                        mg_kws=None, lw=2, main=True, make_single=True,
                        make_multi=True, savefig_single=True,
                        savefig_multi=True, savedir=None, overwrite=False):
    """Make single panel and multi-panel emission line spectra plots.

    Args:
        dapdata: dap.DAP object.
        bins: Bin numbers to plot. Default is (0).
        pnames (tuple): Plot names. Default is None.
        win_cen (tuple): Central wavelength of each panel. Default is None.
        mg_kws (dict): Keyword args with identifying information about the
            galaxy and analysis run. Default is None.
        lw (int): Linewidth. Default is 2.
        main (bool): True is running as script. False is running interactively.
            Default is True.
        make_single (bool): Make single panel plots. Default is True.
        make_multi (bool): Make multi-panel plot. Default is True.
        savefig_single (bool): Save single panel plots. Default is True.
        savefig_multi (bool): Save multi-panel plot. Default is True.
        savedir (str): Directory to save plots in. If None, then default to data
            directory. Default is None.
        overwrite (bool): Overwrite plot if it exists. Default is False.
    """
    bins = convert_bins_from_string_to_list(dapdata, bins)
    for bin in bins:
        mg_kws['bin'] = bin
        if make_multi:
            fig = plot_emline_spectra_multi(dapdata=dapdata, bin=bin,
                                            win_cen=win_cen, pnames=pnames,
                                            mg_kws=mg_kws, lw=lw)
            if savefig_multi:
                util.saveplot(name='multi', path_data=dapdata.path_data,
                              category='emline_spectra', mg_kws=mg_kws,
                              main=main, mkdir=True, savedir=savedir,
                              overwrite=overwrite)
            if main:
                plt.close(fig)

        if make_single:
            for wc, pname in zip(win_cen, pnames):
                nii = True if pname == 'nii' else False
                fig = plot_emline_spectrum(dapdata=dapdata, bin=bin, win_cen=wc,
                                           nii=nii, lw=lw)
                if savefig_single:
                    util.saveplot(name=pname, path_data=dapdata.path_data,
                                  category='emline_spectra', mg_kws=mg_kws,
                                  main=main, mkdir=True, savedir=savedir,
                                  overwrite=overwrite)
                if main:
                    plt.close(fig)

def plot_emline_spectra_multi(dapdata, bin=0, win_cen=None, pnames=None,
                              mg_kws=None, lw=2, figsize=(20, 12)):
    """Plot multiple panel zoom-ins of spectra near strong emission lines.

    Args:
        dapdata: dap.DAP object.
        bin (int): Bin number. Default is 0.
        win_cen (tuple): Central wavelength of each panel. Default is None.
        pnames (tuple): Plot names. Default is None.
        mg_kws (dict): Keyword args with identifying information about the
            galaxy and analysis run. Default is None.
        lw (int): Linewidth. Default is 2.
        figsize (tuple): Figure width and height in inches. Default is (20, 12).

    Returns:
        plt.figure axis object
    """
    fig = plt.figure(figsize=figsize)
    bigax_kws = dict(xlabel=r'$\lambda \, [\AA]$',
                     ylabel=_set_flux_units(dapdata), mg_kws=mg_kws,
                     title_kws=dict(label=_make_spec_title(dapdata, bin),
                                    fontsize=20))
    bigAxes = _make_big_axes(fig, axloc=(0.04, 0.06, 0.9, 0.9), **bigax_kws)
    fig.subplots_adjust(wspace=0.2, hspace=0.15, left=0.1, bottom=0.1,
                        right=0.95, top=0.95)
    for i, (wc, pname) in enumerate(zip(win_cen, pnames)):
        ax = fig.add_subplot(2, 3, i+1)
        nii = True if pname == 'nii' else False
        ig = plot_emline_spectrum(dapdata=dapdata, fig=fig, ax=ax, bin=bin,
                                  win_cen=wc, nii=nii, lw=lw)

    return fig

def plot_emline_spectrum(dapdata, fig=None, ax=None, bin=0, win_cen=None,
                         xlim=None, ylim=None, nii=False, lw=2,
                         figsize=(10, 8)):
    """Plot data and model spectra near strong emission lines.

    Args:
        dapdata: dap.DAP object.
        fig: plt.figure object. Default is None.
        ax: axis object. Default is None.
        bin (int): Bin number. Default is 0.
        win_cen (tuple): Central wavelength of each panel. Default is None.
        xlim (list): Minimum and maximum x-axis values. Default is None.
        ylim (list): Minimum and maximum y-axis values. Default is None.
        nii (bool): If True, maximum y value is determined by the[NII]6548,6583
            lines (not Halpha). Default is False.
        lw (int): Linewidth. Default is 2.
        figsize (tuple): Figure width and height in inches. Default is (10, 8).

    Returns:
        plt.figure axis object
    """
    colors = _set_spec_line_prop(lw=lw)

    if ax is None:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_axes([0.17, 0.15, 0.72, 0.75])
        ax.set_xlabel(r'$\lambda \, [\AA]$', fontsize=24)
        ax.set_ylabel(_set_flux_units(dapdata), fontsize=24)
        # x_offs = [-9, -3, -9, -9, -12, -3, -5, -12, -5]
        x_offs = [-24, 10, -8, -3, -9, -9, -9, -9, -12, -3, -5, -12, -5]
        indiv_ax = True
    else:
        # x_offs = [-12, -4, -12, -12, -20, -4, -12, -20, -6]
        x_offs = [-38, 10, -16, -4, -12, -12, -12, -12, -20, -4, -12, -20, -6]
        indiv_ax = False

    wave = dapdata.wave_rest[bin]
    flux = dapdata.flux_rest[bin]
    ivar = dapdata.ivar_rest[bin]
    noise = ivar**-0.5
    stmodel = dapdata.smod_rest[bin]
    fullfit_ew = dapdata.fullfit_ew_rest[bin]
    fullfit_fb = dapdata.fullfit_fb_rest[bin]

    #linenames = ['[OII]3727', r'H$\beta$', '[OIII]4959', '[OIII]5007',
    #'[NII]6548', r'H$\alpha$', '[NII]6583', '[SII]6716', '[SII]6731']

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
    p_kws = dict(alpha=0.75, lw=lw)
    pst = ax.plot(wave[ind], stmodel[ind], color=colors[0], **p_kws)[0]
    pFB = ax.plot(wave[ind], fullfit_fb[ind], color=colors[1], **p_kws)[0]
    pEW = ax.plot(wave[ind], fullfit_ew[ind], color=colors[2], **p_kws)[0]

    for name, w, x_off in zip(dapdata.elopar.elname_tex,
                              dapdata.elopar.restwave, x_offs):
        if ((w > xmin) and (w < xmax) and
            (name not in ['[OII]3727', '[OII]3729'])):
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

    return fig

def plot_gradient_multi(dapdata, values, errors, columns, leg_kws=None,
                        titles=None, labels=None, mg_kws=None,
                        figsize=(20, 12)):
    """
    Plot multiple radial gradients at once.

    Args:
        mg_kws (dict): MaNGA analysis ID keyword args to pass to
            _make_map_title. Default is None.

    Returns:
        tuple: (plt.figure object, plt.figure axis object)
    """
    mg_kws = util.none_to_empty_dict(mg_kws)

    fig = plt.figure(figsize=figsize)
    if 'seaborn' in sys.modules:
        sns.set_context('poster', rc={'lines.linewidth': 2})

    bigax_kws = dict(xlabel='R [arcsec]', ylabel=_set_flux_units(dapdata),
                     title_kws=dict(fontsize=20), mg_kws=mg_kws)
    bigAxes = _make_big_axes(fig, **bigax_kws)

    for i, column in enumerate(columns):
        left, bottom = _multi_panel_loc(panelnum=i, n_panels=len(columns))
        ax = fig.add_axes([left, bottom, 0.23, 0.33333])
        ig = plot_gradient(dapdata, values, errors, column, fig=fig, ax=ax,
                           leg_kws=leg_kws, title=titles[column], labels=labels)

    return fig



def plot_gradient(dapdata, values, errors, column, fig=None, ax=None,
                  leg_kws=None, title=None, labels=None, figsize=(10, 8)):
    """
    Plot radial gradients.
    """
    if leg_kws is None:
        leg_kws = dict(handlelength=2, loc=1)

    colors = _set_spec_line_prop(lw=2)

    if ax is None:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_axes([0.17, 0.11, 2/3., 5/6.])
        ax.set_xlabel('R [arcsec]', fontsize=28)
        ax.set_ylabel(_set_flux_units(dapdata), fontsize=28)
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(20)
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(20)
        title_kws = dict(fontsize=28)
    else:
        title_kws = dict(fontsize=20)

    ax.set_title(title, **title_kws)

    if 'seaborn' not in sys.modules:
        ax.set_axis_bgcolor('#A8A8A8')
        colors = ['b', 'r', 'c']

    p = []
    lab = []
    for value, error, color, label in zip(values, errors, colors[1:3], labels):
        binr = dapdata.bins.binr.values
        val = dapdata.__dict__[value][column]
        err = dapdata.__dict__[error][column]
        p.append(ax.hlines(val, dapdata.bins.binxrl, dapdata.bins.binyru,
                           color=color))
        ax.plot(binr, val, color=color, zorder=8, lw=0.5)
        ax.scatter(binr, val, facecolor=color, edgecolor='None', s=60, zorder=9)
        # Following line was commented out in MPL-3 run
        ax.errorbar(binr, val, yerr=err, ecolor=color, elinewidth=1,
                    marker='None', ls='None')
        lab.append(label)

    leg = plt.legend(p, lab, **leg_kws)
    plt.setp(leg.get_texts(), fontsize=24)
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)

    # This section was used to plot errors in MPL-3 run
    # if not np.isnan(d['val']).all():
    #     for kk, kkerr, j in zip(['val2', 'val'], ['val2_err', 'val_err'], [2, 0]):
    #         ax.errorbar(self.binr, args[gradname][kk], yerr=args[gradname][kkerr],
    #                     ecolor=c[j], elinewidth=1, marker='None', ls='None')

    return fig



def plot_gradients(dapdata, values, errors, columns, plotname=None, mg_kws=None,
                   leg_kws=None, titles=None, labels=None, figsize=(10, 8),
                   main=True, make_single=True, make_multi=True,
                   savefig_single=True, savefig_multi=True, savedir=None,
                   overwrite=False):
    """wrapper around plot_gradient and plot_gradient_multi

    Args:
        savedir (str): Directory to save plots in. If None, then default to data
            directory. Default is None.
    """
    if make_multi:
        fig = plot_gradient_multi(dapdata, values, errors, columns,
                                  leg_kws=leg_kws, titles=titles, labels=labels,
                                  mg_kws=mg_kws)
        if savefig_multi:
            util.saveplot(name=plotname, path_data=dapdata.path_data,
                          category='gradients', mg_kws=mg_kws, main=main,
                          mkdir=True, savedir=savedir, overwrite=overwrite)
        if main:
            plt.close(fig)

    for column in columns:
        if make_single:
            fig = plot_gradient(dapdata, values, errors, column,
                                leg_kws=leg_kws, title=titles[column],
                                labels=labels, figsize=figsize)
            if savefig_single:
                pname = '_'.join([plotname, column])
                util.saveplot(name=pname, path_data=dapdata.path_data,
                              category='gradients', mg_kws=mg_kws, main=main,
                              mkdir=True, savedir=savedir, overwrite=overwrite)
            if main:
                plt.close(fig)

def make_plots(plottype, dapdata, mg_kws, plot_kws):
    """general wrapper for all plotting functions"""
    print()
    print(plottype)
    ptype = '_'.join(plottype.split('_')[:-1])
    category = plottype.split('_')[-1]
    if category == 'spectra':
        if ptype == 'full':
            plot_spectra(dapdata=dapdata, mg_kws=mg_kws, **plot_kws)
        elif ptype == 'emline':
            plot_emline_spectra(dapdata=dapdata, mg_kws=mg_kws, **plot_kws)
    elif (category == 'gradients') and (mg_kws['bintype'] == 'RADIAL'):
        plot_gradients(dapdata=dapdata, mg_kws=mg_kws, plotname=ptype,
                       **plot_kws)
    elif (category == 'maps') and (mg_kws['bintype'] in ['NONE', 'STON']):
        plot_maps(dapdata=dapdata, mg_kws=mg_kws, **plot_kws)
