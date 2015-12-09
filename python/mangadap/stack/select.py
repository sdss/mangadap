from __future__ import (division, print_function, absolute_import,
                        unicode_literals)
import operator as op
import numpy as np
import pandas as pd


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





#---- UNTESTED -----

def set_value_type(value, type='float'):
    """Convert value from string to another type.

    Args:
        value (str): Input value.
        type (str): Data type to convert value to. Default is 'float'.

    Returns:
        Value with new type.
    """
    if type == 'float':
        return float(value)
    elif type == 'int':
        return int(value)
    elif type == 'str':
        return str(value)
    else:
        print('No value type specified. Returning input value.')
        return value

def apply_selection_condition(cfg_in, data_refs):
    """Apply selection condition.

    Args:
        cfg_in (list): Strings the specify how to do selection.
        data_refs (dict): Mapping between name of data source and actual data
            source.

    Returns:
        boolean array
    """
    data_obj_name, column_in, operator, value_in, value_type = cfg_in
    value = set_value_type(value_in, value_type)
    column = column_in.split('.')
    data_obj = data_refs[data_obj_name]
    for i in range(len(column)-1):
        data_obj = data_obj[column[i]]
    ind_bool = op.__dict__[operator](data_obj[column[-1]], value)
    if isinstance(ind_bool, pd.Series):
        return ind_bool.values
    else:
        return ind_bool

def do_selection(raw_cfg_in, data_refs):
    """Do selection by applying multiple conditions.

    Args:
        raw_cfg_in (list): List of lists of string input from config file or
            Marvin.
        data_refs (dict): Mapping between name of data source and actual data
            source.

    Returns:
        boolean array
    """
    conditions = []
    for item in raw_cfg_in:
        conditions.append(apply_selection_condition(item, data_refs))
    # DROP NANs?
    return join_logical_and(conditions)

    