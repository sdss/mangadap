from __future__ import division, print_function, absolute_import

import os
from os.path import join

import numpy as np
import pandas as pd
from astropy.io import fits
from matplotlib.colors import LinearSegmentedColormap

from sdss.files import base_path

def fitsrec_to_dataframe(recarr):
    """Convert astropy FITS_rec to pandas DataFrame.

    Args:
        recarr (astropy.io.fits.FITS_rec): FITS record array class.

    Returns:
        DataFrame
    """
    cols = recarr.columns.names
    dtmp = {}
    for col in cols:
        dtmp[col] = recarr[col].byteswap().newbyteorder()
    return pd.DataFrame(dtmp, columns=cols)

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
    df = pd.DataFrame(recarr[ext].byteswap().newbyteorder(), columns=columns)
    return df

def swap_byte(arr):
    """Swap byte order from big-endian (FITS) to little-endian (pandas)."""
    return arr.byteswap().newbyteorder()


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

def none_to_empty_dict(x):
    """If a variable is None, return an empty dictionary."""
    if x is None:
        x = {}
    return x

def output_path(name, path_data, plottype, mg_kws, mkdir=False):
    """Make plot output path and file name.

    Args:
        name (str): Plot name.
        path_data (str): Path to parent directory of *plots/*.
        plottype (str): Type of plot ('map', 'spectra', or 'gradients').
        mg_kws (dict): MaNGA galaxy and analysis information.
        mkdir (bool): Make directory if it does not exist. Default is False.

    Returns:
        str: Plot output path and file name.
    """
    filename = ('manga-{plate}-{ifudesign}-LOG{mode}_BIN-{bintype}-{niter}'
                '_{0}.png'.format(name, **mg_kws))
    path_plottype = join(path_data, 'plots', plottype)
    fullpath = join(path_plottype, filename)
    if mkdir:
        if not os.path.isdir(path_plottype):
            os.makedirs(path_plottype)
            print('\nCreated directory: {}\n'.format(path_plottype))
    return fullpath

def reverse_cmap(x):
    """Reverse colormap."""
    def out(y):
        return x(1. - y)
    return out

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

    LinearL_r = {}
    for k in LinearL:
        LinearL_r[k] = reverse_cmap(LinearL[k])

    cmap = LinearSegmentedColormap('linearL', LinearL)
    cmap_r = LinearSegmentedColormap('linearL_r', LinearL_r)

    return (cmap, cmap_r)

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

def make_data_path(file_kws):
    """Make path to data files.

    Args:
        file_kws (dict): Parameters that specify DAP FITS file.

    Returns:
        str: Path to data.
    """
    cfg_dir = join(os.getenv('MANGADAP_DIR'), 'python', 'mangadap', 'plot',
                   'config')
    paths_cfg = join(cfg_dir, 'sdss_paths.ini')
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
    """Convert a FITS recarray into a MultiIndex DataFrame."""
    dt = np.concatenate([rec[c].byteswap().newbyteorder().T for c in cols1]).T
    cols_out = pd.MultiIndex.from_product([cols1, cols2])
    return pd.DataFrame(dt, columns=cols_out)

def arr_to_multiindex_df(arr, cols1, cols2):
    """Convert a 3D array into a MultiIndex DataFrame."""
    data = np.concatenate([arr[i] for i in range(len(cols1))]).T
    cols_out = pd.MultiIndex.from_product([cols1, cols2])
    return pd.DataFrame(data, columns=cols_out)
