from __future__ import (division, print_function, absolute_import,
                        unicode_literals)
import numpy as np
# import pandas as pd


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



# Don't do nan selection or bin selection at this stage


def mean(val):
    """Compute mean.

    Args:
        val (Series): Data.

    Returns:
        float
    """
    return np.mean(val)

def median(val):
    """Compute median.

    Args:
        val (Series): Data.

    Returns:
        float
    """
    return np.median(val)

def ivar_wtmean(val, err, bins):
    """Compute inverse variance weighted mean.

    Args:
        val (Series): Data.
        err (Series): Errors.
        bins (array): Bins to combine.

    Returns:
        float
    """
    return val.dot(err**-2.) / np.sum(err**-2.)


