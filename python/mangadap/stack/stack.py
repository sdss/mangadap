from __future__ import (division, print_function, absolute_import,
                        unicode_literals)
import numpy as np
# import pandas as pd

from mangadap.stack import select


def combine(val, err=None, bins=None, weights=None, scheme=None):
    """TEST
    """
    if scheme == 'median':
        return median(val, bins)
    elif scheme == 'mean':
        weights = None
    elif scheme == 'ivar_wtmean':
        # CHECK THIS
        weights = err.loc[bins]**-2. / np.sum(err.loc[bins]**-2.)
        # return ivar_wtmean(val, err, bins)
    return weighted_average(val, err, bins, weights)

def weighted_average(val, err=None, bins=None, weights=None):
    """combine with weights (mean = uniform weighting)"""
    pass






def mean(val):
    """Compute mean.

    Args:
        val (Series): Data.

    Returns:
        float
    """
    return np.mean(val)

def median(val, bins):
    """Compute median.

    Args:
        val (Series): Data.
        bins (array): Bins to combine.

    Returns:
        float
    """
    ind_notnan = select.get_notnan(val)
    b = select.join_logical_and([bins, ind_notnan])
    return np.median(val[b])

def ivar_wtmean(val, err, bins):
    """Compute inverse variance weighted mean.

    Args:
        val (Series): Data.
        err (Series): Errors.
        bins (array): Bins to combine.

    Returns:
        float
    """
    # Convert integer index array to boolean index array in config reading step
    # bins_bool = select.int_to_bool_index(bins, val.shape)
    ind_notnan = select.get_notnan(val)
    b = select.join_logical_and([bins, ind_notnan])
    return val.loc[b].dot(err.loc[b]**-2.) / np.sum(err.loc[b]**-2.)


