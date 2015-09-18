"""Methods for plotting DAP output.
"""

from __future__ import division, print_function, absolute_import

import sys
import copy

import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.ticker import MaxNLocator

from astropy.stats import sigma_clip

import util

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
    if err is not None:
        no_measure[(err == 0.)] = True
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
    """

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
        print('AttributeError: MaxNLocator instance has no attribute' +
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

def set_panel_par():
    """Set default plot parameters."""
    ax_kws = dict(facecolor='#A8A8A8')
    imshow_kws = dict(cmap=cm.Blues_r) # cm.RdBu
    return ax_kws, imshow_kws

def set_single_panel_par():
    """Set default parameters for a single panel plot.

    Returns:
        tuple: (figure keyword args, axes keyword args, title keyword args,
                imshow keyword args, colorbar keyword args)
    """
    ax_kws, imshow_kws = set_panel_par()
    fig_kws = dict(figsize=(10, 8))
    title_kws = dict(fontsize=28)
    cb_kws = dict(axloc=[0.82, 0.1, 0.02, 5/6.],
                  cbrange=None, sigclip=3, symmetric=False,
                  label_kws=dict(label=None, size=20),
                  tick_params_kws=dict(labelsize=20))
    return fig_kws, ax_kws, title_kws, imshow_kws, cb_kws

def set_multi_panel_par():
    """Set default parameters for a multi panel plot.

    Returns:
        tuple: (figure keyword args, axes keyword args, title keyword args,
                imshow keyword args, colorbar keyword args)
    """
    ax_kws, imshow_kws = set_panel_par()
    fig_kws = dict(figsize=(20, 12))
    title_kws = dict(fontsize=20)
    cb_kws = dict(cbrange=None, sigclip=3, symmetric=False,
                  label_kws=dict(label=None, size=16),
                  tick_params_kws=dict(labelsize=16))
    return fig_kws, ax_kws, title_kws, imshow_kws, cb_kws

def ax_setup(fig=None, ax=None, fig_kws=None, facecolor='#EAEAF2'):
    """Basic axes setup.

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
    mg_kws = util.none_to_empty_dict(title_kws)

    # make axis without a frame or x and y ticks
    bigAxes = fig.add_axes(axloc, frameon=False)
    bigAxes.set_xticks([])
    bigAxes.set_yticks([])

    if xlabel:
        bigAxes.set_xlabel(xlabel, fontsize=labelsize)
    if ylabel:
        bigAxes.set_ylabel(ylabel, fontsize=labelsize)

    # set title
    if mg_kws:
        title_kws['label'] = make_map_title(mg_kws)
    if 'label' in title_kws:
        bigAxes.set_title(**title_kws)
    return bigAxes

def show_bin_num(binxrl, binyru, nbin, val, ax, imshow_kws, fontsize=6):
    """Display bin number on map.

    Args:
        binxrl (array):
        binyru (array):
        nbin (array):
        val (array):
        ax :
        image:
        fontsize (int): Nominal font size. Defaults to 6.

    Returns:
        plt.figure axis object
    """
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

def plot_map(image, extent, xy_nomeasure=None, fig=None, ax=None,
             fig_kws=None, ax_kws=None, title_kws=None, patch_kws=None,
             imshow_kws=None, cb_kws=None, binnum_kws=None, bindot_args=()):
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
    for item in (fig_kws, ax_kws, title_kws, patch_kws, imshow_kws, cb_kws,
                 binnum_kws):
        item = util.none_to_empty_dict(item)

    fig, ax = ax_setup(fig, ax, fig_kws=fig_kws, **ax_kws)

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
        ax = show_bin_num(ax=ax, **binnum_kws)

    if bindot_args:
        ax.plot(*bindot_args, color='k', marker='.', markersize=3, ls='None',
                zorder=10)

    return fig, ax

def plot_multi_map(all_panel_kws, patch_kws=None, fig_kws=None, mg_kws=None):
    """Make multi-panel map plot.

    Args:
        all_panel_kws (dict): Collection of keyword args for each panel
        patch_kws (dict): Keyword args to pass to ax.add_patch. Default is
            None.
        fig_kws (dict): Keyword args to pass to plt.figure. Default is None.
        mg_kws (dict): MaNGA analysis ID keyword args to pass to
            make_map_title. Default is None.

    Returns:
        tuple: (plt.figure object, plt.figure axis object)
    """
    for item in (patch_kws, fig_kws, mg_kws):
        item = util.none_to_empty_dict(item)

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

def make_plots(columns, values, errors, extent, xpos, ypos, binid, binxrl,
               binyru, nbin, spaxel_size, delta, dapdata=None,
               val_no_measure=0, snr_thresh=1, mg_kws=None, patch_kws=None,
               titles=None, cblabels=None, show_binnum=False,
               show_bindot=False):
    """Make single panel plots and multi-panel plot for set of measurements.

    Add options to save to file.

    Args:
       columns (list): Columns of values and errors DataFrames to plot.
       values: Either string that references an attribute from dapdata or an
           array of values.
       errors: Either string that references an attribute from dapdata or an
           array of values.
       extent (array): Extent (xmin, xmax, ymin, ymax) of map in arcsec.
       xpos (array): x-coordinates of bins.
       ypos (array): y-coordinates of bins.
       binid (array): Bin ID numbers.
       binxrl (array): Luminosity-weighted on-sky x-coordinates of the binned
           spectra in arcsec.
       binyru (array): Luminosity-weighted on-sky y-coordinates of the binned
           spectra in arcsec.
       nbin (array): Number of DRP-produced spectra in each bin.
       spaxel_size (float): Spaxel size in arcsec.
       delta (float): Half of the spaxel size in arcsec.
       dapdata: dap.DAP object. Defaults to None.
       val_no_measure (float): Value that corresponds to no measurement.
           Defaults to 0.
       snr_thresh (float): Signal-to-noise threshold for displaying a bin on a
           map. Defaults to 1.
       mg_kws (dict): Keyword args with identifying information about the
           galaxy and analysis run. Default is None.
       patch_kws (dict): Keyword args for drawing hatched regions. Default is
           None.
       titles (list): Plot title for each map. Defaults to None.
       cblabels (list): Colorbar labels. Defaults to None.
       show_binnum (bool): Show bin numbers. Defaults to False.
       show_bindot (bool): Show bin dots. Defaults to False.
    """
    for item in (values, errors):
        if isinstance(item, str):
            item = getattr(dapdata, item)

    for item in (patch_kws, mg_kws):
        item = util.none_to_empty_dict(item)

    images = []
    xy_nomeasures = []
    for col in columns:
        im, xy = make_image(val=values[col].values, err=errors[col].values,
                            xpos=xpos, ypos=ypos, binid=binid, delta=delta,
                            val_no_measure=val_no_measure,
                            snr_thresh=snr_thresh)
        images.append(im)
        xy_nomeasures.append(xy)

    fig_kws, ax_kws, title_kws, imshow_kws, cb_kws = set_multi_panel_par()
    all_panel_kws = []
    for i, (im, xy, col) in enumerate(zip(images, xy_nomeasures, columns)):
        plot_title = titles[col]
        cblabel = cblabels[col]
        binnum_kws = {}
        bindot_args = ()
        if show_binnum:
            binnum_kws = dict(binxrl=binxrl, binyru=binyru, nbin=nbin,
                              val=values[col].values, spaxel_size=spaxel_size,
                              imshow_kws=imshow_kws, fontsize=6)
        if show_bindot:
            bindot_args = (-binxrl, binyru)

        # DEBUGGING: only plot one map
        if i == 0:
            # plot single panel maps
            fg, axs, tt, iw, cb = set_single_panel_par()
            tt['label'] = plot_title
            cb['label_kws']['label'] = cblabel
            ig = plot_map(im, extent, xy_nomeasure=xy, fig_kws=fg, ax_kws=axs,
                          title_kws=tt, patch_kws=patch_kws, imshow_kws=iw,
                          cb_kws=cb, binnum_kws=binnum_kws,
                          bindot_args=bindot_args)

        # create dictionaries for multi-panel maps
        kwdicts = (ax_kws, title_kws, imshow_kws, cb_kws)
        # ADD: allow for different color maps in multi-panel plot
        a_kws, t_kws, i_kws, c_kws = [copy.deepcopy(it) for it in kwdicts]
        t_kws['label'] = plot_title
        c_kws['label_kws']['label'] = cblabel
        all_panel_kws.append(dict(image=im, extent=extent, xy_nomeasure=xy,
                                  ax_kws=a_kws, title_kws=t_kws,
                                  imshow_kws=i_kws, cb_kws=c_kws))

    # plot multi-panel maps
    ig = plot_multi_map(all_panel_kws=all_panel_kws, fig_kws=fig_kws,
                        patch_kws=patch_kws, mg_kws=mg_kws)

