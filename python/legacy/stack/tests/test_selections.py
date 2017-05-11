from __future__ import (division, print_function, absolute_import,
                        unicode_literals)
import numpy as np
import pandas as pd

import unittest
from numpy.testing import assert_array_equal
from pandas.util.testing import assert_series_equal
# from pandas.util.testing import assert_frame_equal

from mangadap.stack import selections

class SelectTestCase(unittest.TestCase):

    def test_int_to_bool_index_empty_list(self):
        int_ind = np.array([])
        arr_shape = 3
        desired = np.array([False, False, False])
        actual = selections.int_to_bool_index(int_ind, arr_shape)
        assert_array_equal(actual, desired)

    def test_int_to_bool_index_partial_list(self):
        int_ind = np.arange(2)
        arr_shape = 3
        desired = np.array([True, True, False])
        actual = selections.int_to_bool_index(int_ind, arr_shape)
        assert_array_equal(actual, desired)
    
    def test_int_to_bool_index_complete_list(self):
        int_ind = np.arange(3)
        arr_shape = 3
        desired = np.array([True, True, True])
        actual = selections.int_to_bool_index(int_ind, arr_shape)
        assert_array_equal(actual, desired)

    def test_join_logical_and_empty_list(self):
        desired = []
        actual = selections.join_logical_and(desired)
        assert_array_equal(actual, desired)
    
    def test_join_logical_and_single_condition(self):
        desired = np.array([True, False])
        actual = selections.join_logical_and([desired])
        assert_array_equal(actual, desired)

    def test_join_logical_and_mixed_conditions(self):
        ind_bool = [np.array([True, True]), np.array([True, False])]
        desired = np.array([True, False])
        actual = selections.join_logical_and(ind_bool)
        assert_array_equal(actual, desired)

    def test_join_logical_and_all_true(self):
        ind_bool = [np.array([True, True]), np.array([True, True])]
        desired = np.array([True, True])
        actual = selections.join_logical_and(ind_bool)
        assert_array_equal(actual, desired)

    def test_join_logical_and_all_false(self):
        ind_bool = [np.array([False, False]), np.array([False, False])]
        desired = np.array([False, False])
        actual = selections.join_logical_and(ind_bool)
        assert_array_equal(actual, desired)

    def test_join_logical_or_empty_list(self):
        desired = []
        actual = selections.join_logical_or(desired)
        assert_array_equal(actual, desired)
    
    def test_join_logical_or_single_condition(self):
        desired = np.array([True, False])
        actual = selections.join_logical_or([desired])
        assert_array_equal(actual, desired)

    def test_join_logical_or_mixed_conditions(self):
        ind_bool = [np.array([True, True]), np.array([True, False])]
        desired = np.array([True, True])
        actual = selections.join_logical_or(ind_bool)
        assert_array_equal(actual, desired)

    def test_join_logical_or_all_true(self):
        ind_bool = [np.array([True, True]), np.array([True, True])]
        desired = np.array([True, True])
        actual = selections.join_logical_or(ind_bool)
        assert_array_equal(actual, desired)

    def test_join_logical_or_all_false(self):
        ind_bool = [np.array([False, False]), np.array([False, False])]
        desired = np.array([False, False])
        actual = selections.join_logical_or(ind_bool)
        assert_array_equal(actual, desired)

    def test_get_notnan_np_nan_only(self):
        data = np.array([np.nan, 1, 2])
        desired = np.array([False, True, True])
        actual = selections.get_notnan(data)
        assert_array_equal(actual, desired)

    def test_get_notnan_special_nan_value(self):
        data = np.array([-9999, 1, 2])
        desired = np.array([False, True, True])
        actual = selections.get_notnan(data, nanvals=-9999)
        assert_array_equal(actual, desired)

    def test_get_notnan_special_nan_values(self):
        data = np.array([-99, -9999, 2])
        desired = np.array([False, False, True])
        actual = selections.get_notnan(data, nanvals=[-99, -9999])
        assert_array_equal(actual, desired)

    def test_get_notnan_mixed_nan_value(self):
        data = np.array([np.nan, -9999, 2])
        desired = np.array([False, False, True])
        actual = selections.get_notnan(data, nanvals=-9999)
        assert_array_equal(actual, desired)

    def test_get_notnan_mixed_nan_values(self):
        data = np.array([np.nan, -99, -9999, 2])
        desired = np.array([False, False, False, True])
        actual = selections.get_notnan(data, nanvals=[-99, -9999])
        assert_array_equal(actual, desired)

    def test_cfg_to_notnan(self):
        d = dict(key1=dict(key2=np.arange(5)))
        data_refs = dict(data=d)
        cfg_in = ['data', 'key1.key2', '3']
        desired = np.array([True, True, True, False, True])
        actual = selections.cfg_to_notnan(cfg_in, data_refs)
        assert_array_equal(actual, desired)

    def test_set_value_type_none(self):
        desired = 'foo'
        actual = selections.set_value_type('foo', value_type=None)
        self.assertEqual(actual, desired)

    def test_set_value_type_empty_string(self):
        desired = 'foo'
        actual = selections.set_value_type('foo', value_type='')
        self.assertEqual(actual, desired)

    def test_set_value_type_float(self):
        desired = 7.5
        actual = selections.set_value_type('7.5', value_type='float')
        self.assertEqual(actual, desired)

    def test_set_value_type_int(self):
        desired = 7
        actual = selections.set_value_type('7', value_type='int')
        self.assertEqual(actual, desired)

    def test_set_value_type_str(self):
        desired = 'foo'
        actual = selections.set_value_type('foo', value_type='str')
        self.assertEqual(actual, desired)

    def test_get_multilevel_attribute(self):
        d = dict(key1=dict(key2=np.arange(5)))
        desired = np.arange(5)
        actual = selections.get_multilevel_attribute(keys=['key1', 'key2'],
                                                     data=d)
        assert_array_equal(actual, desired)

    def test_apply_selection_condition_df(self):
        cfg_in = ['df', 'column1', 'gt', '2.5', 'float']
        df = pd.DataFrame(dict(column1=np.arange(5), column2=np.arange(5)*2))
        data_refs = dict(df=df)
        desired = np.array([False, False, False, True, True])
        actual = selections.apply_selection_condition(cfg_in, data_refs)
        assert_array_equal(actual, desired)

    def test_apply_selection_condition_dict(self):
        cfg_in = ['dict', 'column1', 'gt', '2.5', 'float']
        d = dict(column1=np.arange(5), column2=np.arange(5)*2)
        data_refs = dict(dict=d)
        desired = np.array([False, False, False, True, True])
        actual = selections.apply_selection_condition(cfg_in, data_refs)
        assert_array_equal(actual, desired)

    def test_apply_selection_condition_recarr(self):
        cfg_in = ['recarray', 'x', 'gt', '2.5', 'float']
        ra = np.rec.array([(1, 5), (6, 7)], dtype=[('x', 'int'), ('y', 'int')])
        data_refs = dict(recarray=ra)
        desired = np.array([False, True])
        actual = selections.apply_selection_condition(cfg_in, data_refs)
        assert_array_equal(actual, desired)

    def test_do_selection(self):
        data = dict(x=np.arange(5), y=np.arange(5))
        data_refs = dict(d=data)
        inp = [['d', 'x', 'gt', '1.5', 'float'],
               ['d', 'y', 'mod', '2', 'int']]
        desired = np.array([False, False, False, True, False])
        actual = selections.do_selection(inp, data_refs)
        assert_array_equal(actual, desired)

if __name__ == '__main__':
    unittest.main()
