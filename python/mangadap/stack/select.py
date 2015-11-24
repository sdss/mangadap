from __future__ import (division, print_function, absolute_import,
                        unicode_literals)
import numpy as np
# import pandas as pd

def int_to_bool_index(ind_int, arr_shape):
    """Convert integer index array to a boolean index array.

    Args:
        ind: integer index or integer index array.
        arr_shape: integer length or tuple of shape to match.

    Returns:
        array: boolean index array.
    """
    ind_bool = np.zeros(arr_shape, dtype=bool)
    if len(ind_int) > 0:
        ind_bool[ind_int] = True
    return ind_bool

def join_conditions(ind_bool_list, operator):
    """Join conditions using logical operators
    
    Args:
        ind_bool_list (list): list of boolean index arrays.
        operator (str): logical operator to do.

    Returns:
       array: boolean index array.
    """
    if len(ind_bool_list) > 0:
        for i, it in enumerate(ind_bool_list):
            if i == 0:
                ind = it
            else:
                if operator is 'and':
                    ind = np.logical_and(ind, it)
                elif operator is 'or':
                    ind = np.logical_or(ind, it)
                elif operator is 'xor':
                    ind = np.logical_xor(ind, it)
                else:
                    raise ValueError('Must select a valid logical operator.')
        return ind
    else:
        return []

def join_logical_and(ind_bool_list):
    """Join conditions using logical AND.

    Args:
        ind_bool_list (list): list of boolean index arrays.

    Returns:
       array: boolean index array.
    """
    return join_conditions(ind_bool_list, operator='and')

def join_logical_or(ind_bool_list):
    """Join conditions using logical OR.

    Args:
        ind_bool_list (list): list of boolean index arrays.

    Returns:
       array: boolean index array.
    """
    return join_conditions(ind_bool_list, operator='or')

def join_logical_xor(ind_bool_list):
    """Join conditions using logical XOR.

    Args:
        ind_bool_list (list): list of boolean index arrays.

    Returns:
       array: boolean index array.
    """
    return join_conditions(ind_bool_list, operator='xor')

def notnan(arr, nanvals=None):
    """Find not NaN values.

    Args:
        nanvals: list of numbers or single number that corresponds to NaN.
        Default is None.

    Returns:
        array: boolean index array.
    """
    ind_bools = [~np.isnan(arr)]
    if nanvals is not None:
        try:
            for nanval in nanvals:
                ind_bools.append(arr != nanval)
        except TypeError:
            ind_bools.append(arr != nanvals)
    return join_logical_and(ind_bools)

    