from __future__ import (division, print_function, absolute_import,
                        unicode_literals)
import numpy as np
# import pandas as pd

from mangadap.stack import select

def ivar_wtmean(val, err, bins):
    """Compute inverse variance weighted mean.

    Args:
        val (Series): Data.
        err (Series): Errors.
        bins (array): Bins to combine.

    Returns:
        float
    """
    bins_bool = select.int_to_bool_index(bins, val.shape)
    notnan = select.notnan(val)
    b = select.join_logical_and([bins_bool, notnan])
    return val.loc[b].dot(err.loc[b]**-2.) / np.sum(err.loc[b]**-2.)
