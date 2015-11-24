from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import unittest
from numpy.testing import assert_array_equal
from pandas.util.testing import assert_frame_equal

import numpy as np
import pandas as pd

from mangadap.stack import select

class SelectTestCase(unittest.TestCase):

    def test_int_to_bool_index_empty_list(self):
        int_ind = np.array([])
        arr_shape = 3
        desired = np.array([False, False, False])
        actual = select.int_to_bool_index(int_ind, arr_shape)
        assert_array_equal(actual, desired)

    def test_int_to_bool_index_partial_list(self):
        int_ind = np.arange(2)
        arr_shape = 3
        desired = np.array([True, True, False])
        actual = select.int_to_bool_index(int_ind, arr_shape)
        assert_array_equal(actual, desired)
    
    def test_int_to_bool_index_complete_list(self):
        int_ind = np.arange(3)
        arr_shape = 3
        desired = np.array([True, True, True])
        actual = select.int_to_bool_index(int_ind, arr_shape)
        assert_array_equal(actual, desired)

    def test_join_logical_and_empty_list(self):
        desired = []
        actual = select.join_logical_and(desired)
        assert_array_equal(actual, desired)
    
    def test_join_logical_and_single_condition(self):
        desired = np.array([True, False])
        actual = select.join_logical_and([desired])
        assert_array_equal(actual, desired)

    def test_join_logical_and_mixed_conditions(self):
        ind_bool = [np.array([True, True]), np.array([True, False])]
        desired = np.array([True, False])
        actual = select.join_logical_and(ind_bool)
        assert_array_equal(actual, desired)

    def test_join_logical_and_all_true(self):
        ind_bool = [np.array([True, True]), np.array([True, True])]
        desired = np.array([True, True])
        actual = select.join_logical_and(ind_bool)
        assert_array_equal(actual, desired)

    def test_join_logical_and_all_false(self):
        ind_bool = [np.array([False, False]), np.array([False, False])]
        desired = np.array([False, False])
        actual = select.join_logical_and(ind_bool)
        assert_array_equal(actual, desired)

    def test_join_logical_or_empty_list(self):
        desired = []
        actual = select.join_logical_or(desired)
        assert_array_equal(actual, desired)
    
    def test_join_logical_or_single_condition(self):
        desired = np.array([True, False])
        actual = select.join_logical_or([desired])
        assert_array_equal(actual, desired)

    def test_join_logical_or_mixed_conditions(self):
        ind_bool = [np.array([True, True]), np.array([True, False])]
        desired = np.array([True, True])
        actual = select.join_logical_or(ind_bool)
        assert_array_equal(actual, desired)

    def test_join_logical_or_all_true(self):
        ind_bool = [np.array([True, True]), np.array([True, True])]
        desired = np.array([True, True])
        actual = select.join_logical_or(ind_bool)
        assert_array_equal(actual, desired)

    def test_join_logical_or_all_false(self):
        ind_bool = [np.array([False, False]), np.array([False, False])]
        desired = np.array([False, False])
        actual = select.join_logical_or(ind_bool)
        assert_array_equal(actual, desired)

    def test_join_logical_xor_empty_list(self):
        desired = []
        actual = select.join_logical_xor(desired)
        assert_array_equal(actual, desired)
    
    def test_join_logical_xor_single_condition(self):
        desired = np.array([True, False])
        actual = select.join_logical_xor([desired])
        assert_array_equal(actual, desired)

    def test_join_logical_xor_mixed_conditions(self):
        ind_bool = [np.array([True, True]), np.array([True, False])]
        desired = np.array([False, True])
        actual = select.join_logical_xor(ind_bool)
        assert_array_equal(actual, desired)

    def test_join_logical_xor_all_true(self):
        ind_bool = [np.array([True, True]), np.array([True, True])]
        desired = np.array([False, False])
        actual = select.join_logical_xor(ind_bool)
        assert_array_equal(actual, desired)

    def test_join_logical_xor_all_false(self):
        ind_bool = [np.array([False, False]), np.array([False, False])]
        desired = np.array([False, False])
        actual = select.join_logical_and(ind_bool)
        assert_array_equal(actual, desired)

    def test_notnan_np_nan_only(self):
        data = np.array([np.nan, 1, 2])
        desired = np.array([False, True, True])
        actual = select.notnan(data)
        assert_array_equal(actual, desired)

    def test_notnan_special_nan_value(self):
        data = np.array([-9999, 1, 2])
        desired = np.array([False, True, True])
        actual = select.notnan(data, nanvals=-9999)
        assert_array_equal(actual, desired)

    def test_notnan_special_nan_values(self):
        data = np.array([-99, -9999, 2])
        desired = np.array([False, False, True])
        actual = select.notnan(data, nanvals=[-99, -9999])
        assert_array_equal(actual, desired)

    def test_notnan_mixed_nan_value(self):
        data = np.array([np.nan, -9999, 2])
        desired = np.array([False, False, True])
        actual = select.notnan(data, nanvals=-9999)
        assert_array_equal(actual, desired)

    def test_notnan_mixed_nan_values(self):
        data = np.array([np.nan, -99, -9999, 2])
        desired = np.array([False, False, False, True])
        actual = select.notnan(data, nanvals=[-99, -9999])
        assert_array_equal(actual, desired)

if __name__ == '__main__':
    unittest.main()
