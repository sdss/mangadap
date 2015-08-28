from __future__ import division, print_function, absolute_import

import numpy as np
import pandas as pd

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


def read_hdu_par(dapf, ext, ltype='emission'):
    """Read emission line or spectral index names.
    
    Args:
        dapf (dapfile): dapfile instance.
        ext (str): Extension to read.
        ltype (str): 'emission' or 'specind'

    Returns:
        array

    """
    par_ext = dict(emission='ELOPAR', specind='SIPAR')
    par = dapf.read_hdu_data(par_ext[ltype])
    return par[ext]

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

def make_mask_no_measurement(data, err=None, val_no_measure=0.,
                             snr_thresh=1.):
    """Mask invalid measurements within a data array.

    Args:
        data (array): Data.
        err (array): Error. Defaults to None.
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

def set_extent(xpos, ypos, delta):
    """Set extent of map.
    """
    return np.array([-xpos.max() - delta, -xpos.min() + delta,
                    ypos.min() - delta, ypos.max() + delta])


def convert_specind_units_to_cblabel(units):
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
        raise('Unknown spectral index units.')
        cblabel = None
    return cblabel


def linear_combination(df, columns, coeffs):
    """Do a linear combination of columns in a DataFrame.

    Args:
        df (DataFrame): DataFrame.
        columns (list): Column names.
        coeffs (list): Coefficients for linear combination.

    Returns:
        DataFrame

    """
    return (df[columns] * coeffs).sum(axis=1)


def linear_combination_err(df, columns, coeffs):
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


def get_spec_ind_units(spec_ind, sinames, units):
    """Determine units of a spectral index measurement.

    Args:
        spec_ind (str): Spectral index name.
        sinames (str): Names of all spectral indices.
        units (array): Units of spectral index measurements.

    Returns:
        str: Units of spectral index measurement.
    """
    ind = np.where(sinames == spec_ind)[0][0]
    return units[ind]