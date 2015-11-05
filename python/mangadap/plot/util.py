# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""Utility functions for DAP plotting."""

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import os
from os.path import join
import sys
import copy
import re

import numpy as np
import pandas as pd
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

from sdss_access.path import base_path

def fitsrec_to_dataframe(recarr, forceswap=False):
    """Convert astropy FITS_rec to pandas DataFrame.

    Args:
        recarr (astropy.io.fits.FITS_rec): FITS record array class.
        forceswap (bool): If True, do a byteswap regardless of python verion.
            Default is False.

    Returns:
        DataFrame
    """
    cols = recarr.columns.names
    dtmp = {}
    for col in cols:
        dtmp[col] = swap_byte(recarr[col])
    return pd.DataFrame(dtmp, columns=cols)

def make_df(data, columns):
    """Try to put data into a DataFrame."""
    try:
        return pd.DataFrame(data, columns=columns)
    except ValueError:
        return data

def read_line_names(dapf, ltype='emission'):
    """Read emission line or spectral index names.

    Args:
        dapf (dapfile): dapfile instance.
        ltype (str): 'emission' or 'specind'

    Returns:
        list

    """
    par_ext = dict(emission='ELOPAR', specind='SIPAR')
    name_ext = dict(emission='ELNAME', specind='SINAME')
    par = dapf.read_hdu_data(par_ext[ltype])
    names = list(par[name_ext[ltype]])
    if ltype == 'emission':
        names = [name.replace('-', '') for name in names]
    return names

def read_vals(dapf, hdu, ext, columns):
    """Read measurements into DataFrame.

    Args:
        dapf (dapfile): dapfile instance.
        hdu (str): HDU name.
        ext (str): Extension name.
        columns (list): Column names.

    Returns:
        DataFrame
    """
    recarr = dapf.read_hdu_data(hdu)
    df = pd.DataFrame(swap_byte(recarr[ext]), columns=columns)
    return df

def swap_byte(arr):
    """Swap byte order from big-endian (FITS) to little-endian (pandas).
    """
    if arr.dtype.byteorder == '>':
        return arr.byteswap().newbyteorder()
    else:
        return arr

def swap_byte_df(arr, columns=None):
    """Swap byte order and convert from array to DataFrame.

    Args:
        arr (array): Array read in from FITS files.
        columns (list): Column names.

    Returns:
        DataFrame
    """
    arr_swapped = swap_byte(arr)
    return pd.DataFrame(arr_swapped, columns=columns)

def lin_comb(df, columns, coeffs):
    """Do a linear combination of columns in a DataFrame.

    Args:
        df (DataFrame): DataFrame.
        columns (list): Column names.
        coeffs (list): Coefficients for linear combination.

    Returns:
        DataFrame

    """
    return (df[columns] * coeffs).sum(axis=1)

def lin_comb_err(df, columns, coeffs):
    """Error propogation for a linear combination of columns in a DataFrame.

    Args:
        df (DataFrame): DataFrame.
        columns (list): Column names.
        coeffs (list): Coefficients for linear combination.
        combination (str): Name of combined column.

    Returns:
        DataFrame

    """
    return np.sqrt((df[columns]**2. * coeffs**2.).sum(axis=1))

def remove_hyphen(names):
    """Remove hyphens from list of strings."""
    return [name.replace('-', '').strip() for name in names]

def lowercase_colnames(df):
    """Convert column names of a DataFrame to lowercase."""
    df.columns = [item.lower() for item in df.columns]
    return df

def texify_elnames(names):
    """Convert emission line names to pretty TeX format."""
    texnames = []
    for name in names:
        if name.lower() in ('hb4862', 'hbeta', 'hb'):
            name = r'H$\beta$'
        elif name.lower() in ('ha6564', 'halpha', 'ha'):
            name = r'H$\alpha$'
        else:
            linewave = re.split('(\d+)', name)
            line, wave = list(filter(None, linewave))
            name = '[' + line + ']' + wave
        texnames.append(name)
    return texnames

def none_to_empty_dict(x):
    """If a variable is None, return an empty dictionary."""
    if x is None:
        x = {}
    return x

def output_path(name, path_data, category, mg_kws, ext='png', mkdir=False):
    """Make plot output path and file name.

    Args:
        name (str): Plot name.
        path_data (str): Path to parent directory of *plots/*.
        category (str): Type of plot ('map,' 'spectra,' 'emline,' or
            'gradients').
        mg_kws (dict): MaNGA galaxy and analysis information.
        ext (str): File extension.
        mkdir (bool): Make directory if it does not exist. Default is False.

    Returns:
        str: Plot output path and file name.
    """
    stem = 'manga-{plateifu}-LOG{mode}_BIN-{bintype}-{niter}'.format(**mg_kws)
    if category == 'maps':
        filename = stem + '_{0}.{1}'.format(name, ext)
    elif category == 'gradients':
        filename = stem + '_{0}.png'.format(name)
    elif category == 'spectra':
        filename = stem + '_spec-{bin:0>4}.png'.format(**mg_kws)
    elif category == 'emline_spectra':
        filename = stem + '_spec-{bin:0>4}_emline_{0}.png'.format(name, **mg_kws)
        category = 'spectra'
    path_category = join(path_data, 'plots', category)
    fullpath = join(path_category, filename)
    if mkdir:
        if not os.path.isdir(path_category):
            os.makedirs(path_category)
            print('\nCreated directory: {}\n'.format(path_category))
    return fullpath

def saveplot(name, path_data, category, mg_kws, ext='png', dpi=200, main=True,
             mkdir=False, overwrite=False):
    """Save a figure.

    Args:
        name (str): Plot name.
        path_data (str): Path to parent directory of *plots/*.
        category (str): Type of plot ('map', 'spectra', or 'gradients').
        mg_kws (dict): MaNGA galaxy and analysis information.
        ext (str): File extension.
        dpi (int): If file is png, specify dots-per-inch. Default is 200.
        main (bool): True is running as script. False is running interactively.
            Default is True.
        mkdir (bool): Make directory if it does not exist. Default is False.
        overwrite (bool): Overwrite plot if it exists. Default is False.
    """
    path = output_path(name=name, path_data=path_data, category=category,
                       mg_kws=mg_kws, ext=ext, mkdir=mkdir)
    if overwrite or not os.path.isfile(path):
        kws = {}
        if ext == 'png':
            kws['dpi'] = dpi
        plt.savefig(path, **kws)
        filename = path.split('/')[-1]
        if not main:
            filename = '\n' + filename + '\n'
        print(filename)

def reverse_cmap(cdict):
    cdict_r = {}
    for k, v in cdict.items():
        data = []
        for it in v:
            data.append((1 - it[0], it[1], it[2]))
        cdict_r[k] = sorted(data)
    return cdict_r

def linear_Lab():
    """Make linear Lab color map.

    Returns:
        tuple: colormap and reversed colormap
    """
    LinL_file = join(os.environ['MANGADAP_DIR'], 'python', 'mangadap',
                     'plot', 'Linear_L_0-1.csv')
    LinL = np.loadtxt(LinL_file, delimiter=',')

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

    LinearL_r = reverse_cmap(LinearL)

    cmap = LinearSegmentedColormap('linearL', LinearL)
    cmap_r = LinearSegmentedColormap('linearL_r', LinearL_r)

    return (cmap, cmap_r)

def get_cmap_rgb(cmap, n_colors=256):
    """Return RGB values of a colormap.

    Args:
        cmap: Colormap.
        n_colors: Number of color tuples in colormap. Default is 256.

    Returns:
        array
    """
    rgb = np.zeros((n_colors, 3))
    for i in range(n_colors):
        rgb[i] = cmap(i)[:3]
    return rgb

def output_cmap_rgb(cmap, path=None, n_colors=256):
    """Print RGB values of a colormap to a file.

    Args:
        cmap: Colormap.
        path: Path to generate output file.
        n_colors: Number of color tuples in colormap. Default is 256.
    """
    rgb = get_cmap_rgb(cmap, n_colors)
    if path is None:
        home = os.path.expanduser('~')
        path = join(home, 'Downloads')
    filename = join(path, '{}.txt'.format(cmap.name))
    header = '{:22} {:24} {:22}'.format('Red', 'Green', 'Blue')
    np.savetxt(filename, rgb, header=header)
    print('Wrote: {}'.format(filename))

def read_drpall(paths_cfg):
    """Read DRPall file.

    Args:
        paths_cfg (str): Path to sdss_paths.ini config file.

    Returns:
        astropy fitsrec: DRPall table.
    """
    bp = base_path(paths_cfg)
    drpall_file = bp.full('drpall')
    fin = fits.open(drpall_file)
    drpall = fin[1].data
    fin.close()
    print('Read {}'.format(drpall_file))
    return drpall

def read_file_list(file_list):
    """Read file list.

    Args:
        file_list (str): Full path to file with list of FITS files to plot.

    Returns:
        dict: FITS file specifications parsed from file name.
    """
    files = np.genfromtxt(file_list, dtype='str')
    files = np.atleast_1d(files)
    f_kws = []
    for item in files:
        stem_file = item.strip('.fits')
        ig, plate, ifudesign, mode_in, bintype, niter = stem_file.split('-')
        mode = mode_in.split('_')[0].strip('LOG')
        f_kws.append(dict(plate=plate, ifudesign=ifudesign, mode=mode,
                          bintype=bintype, niter=niter))
    return f_kws

def make_data_path(paths_cfg, file_kws):
    """Make path to data files.

    Args:
        paths_cfg (str): Full path and file name for sdss_paths.ini file.
        file_kws (dict): Parameters that specify DAP FITS file.

    Returns:
        str: Path to data.
    """
    bp = base_path(paths_cfg)
    return bp.dir('dap', **file_kws)

def make_config_path(filename):
    """Make path to config files.

    Args:
        filename (str): Plot types config file name. If it does not include
            the full path, assume that the path points to the the config
            directory in MANGADAP.

    Returns:
        str: Config file directory.
    """
    if os.path.isfile(filename):
        cfg_dir = os.path.dirname(filename)
    else:
        cfg_dir = join(os.getenv('MANGADAP_DIR'), 'python', 'mangadap', 'plot',
                       'config')
    return cfg_dir

def fitsrec_to_multiindex_df(rec, cols1, cols2):
    """Convert a FITS recarray into a MultiIndex DataFrame.

    Args:
        rec (FITS recarray)
        cols1 (list): First level column names.
        cols2 (list): Second level column names.

    Returns:
        DataFrame
    """
    dt = np.concatenate([swap_byte(rec[c]).T for c in cols1]).T
    cols_out = pd.MultiIndex.from_product([cols1, cols2])
    return pd.DataFrame(dt, columns=cols_out)

def arr_to_multiindex_df(arr, cols1, cols2):
    """Convert a 3D array into a MultiIndex DataFrame.

    Args:
        arr (array): 3D array.
        cols1 (list): First level column names.
        cols2 (list): Second level column names.

    Returns:
        DataFrame
    """
    data = np.concatenate([arr[i] for i in range(len(cols1))]).T
    cols_out = pd.MultiIndex.from_product([cols1, cols2])
    return pd.DataFrame(data, columns=cols_out)

def string_slice_multiindex_df(dapdata, colnames):
    """Use dot notation to access a multi-indexed DataFrame.

    Args:
        dapdata: dap.DAP class instance.
        colnames (str): Column names with each level separated by a period.

    Returns:
        DataFrame
    """
    out = dapdata
    for level in colnames.split('.'):
        out = getattr(out, level)
    return out

def deredshift_velocities(redshift, vel, velerr):
    """Shift velocities to systemic frame.

    Args:
        redshift (float)
        vel (array): Velocities.
        velerr (array): Velocity errors.

    Returns:
        tuple: (rest-frame velocities, rest-frame velocity errors)
    """
    v_light = 299792.458
    vel_rest = vel - (redshift * v_light)
    # WHY IS IT SETTING ERRORS TO ZERO?
    velerr_rest = velerr * 0. + 1e-6
    mask = ((np.abs(vel_rest) > 250.) | (velerr > 250.) | (velerr == 0.))
    velerr_rest[mask] = 1e8
    return vel_rest, velerr_rest

def make_mg_kws(dapdata, file_kws):
    """Create a MaNGA kws dictionary.

    Args:
        dapdata: dap.DAP class instance.
        file_kws (dict): File description.

    Returns:
        dict
    """
    mg_kws = copy.deepcopy(file_kws)
    mg_kws['mangaid'] = dapdata.mangaid
    mg_kws['plateifu'] = '{plate}-{ifudesign}'.format(**mg_kws)
    return mg_kws

def check_h3_h4(bintype, binsn, plottypes_list):
    """If h3 and h4 were measured, use appropriate plottypes config file.

    Args:
        bintype (str): Bintype.
        binsn (float): Minimum S/N for Voronoi binning.
        plottypes_list (str): Plot types config file originally specified.

    Returns:
        str
    """
    if ((bintype == 'STON') and (binsn == 30) and
        (plottypes_list == 'dapqa_plottypes.ini')):
        new_list = 'dapqa_STON30_plottypes.ini'
        print('Switched to {} to make h3 and h4 maps.'.format(new_list))
        return new_list
    else:
        return plottypes_list
