from __future__ import division, print_function, absolute_import

import os
from os.path import join

import numpy as np
import pandas as pd

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

def output_path(name, path_data, plottype, mg_kws):
    """Make plot output path and file name.

    Args:
        name (str): Plot name.
        path_data (str): Path to parent directory of *plots/*.
        plottype (str): Type of plot ('map', 'spectra', or 'gradients').
        mg_kws (dict): MaNGA galaxy and analysis information.

    Returns:
        str: Plot output path and file name.
    """
    filename = ('manga-{plate}-{ifudesign}-LOG{mode}_BIN-{bintype}-{niter}'
                '_{0}.png'.format(name, **mg_kws))
    fullpath = join(path_data, 'plots', plottype, filename)
    return fullpath



# TODO read drp3qual from drpall file
# drp class?

import os
from os.path import join
from astropy.io import fits
from sdss.files import base_path
home = os.path.expanduser('~')
bp = base_path(join(home, 'Dropbox', 'data', 'sdss_paths.ini'))
drpall_file = bp.full('drpall')
fin = fits.open(drpall_file)
tbl = fin[1]
plate = tbl.data['plate'].astype(str)
ifudesign = tbl.data['ifudsgn']
plateifu = 
drp3qual = tbl.data['drp3qual']

    def read_drpall(self):
        """
        Read in drpall FITS file.
        """
        try:
            fin = fits.open(self.drpall_file)
        except IOError:
            raise Exception('Cannot find drpall file: %s' % self.drpall_file)
        else:
            tbl = fin[1]
            plate = tbl.data['plate'].astype(str)
            if 'ifudsgn' in tbl.data.names:
                ifudesign = tbl.data['ifudsgn']
            elif 'ifudesign' in tbl.data.names:
                ifudesign = tbl.data['ifudesign']
            else:
                raise KeyError
            manga_id = tbl.data['mangaid']
            nsa_redshift = tbl.data['nsa_redshift']
            nsa_vdisp = tbl.data['nsa_vdisp']

            pid, ifu = self.manga_pid.split('-')

            ind_tbl = np.where((plate == pid) & (ifudesign == ifu))[0]
            if self.manga_id == 'n/a':
                self.manga_id = manga_id[ind_tbl][0]
            self.nsa_redshift = nsa_redshift[ind_tbl][0]
            self.nsa_vdisp = nsa_vdisp[ind_tbl][0]

