from __future__ import (division, print_function, absolute_import,
                        unicode_literals)
import numpy as np
import pandas as pd

import unittest
from numpy.testing import assert_array_equal
from pandas.util.testing import assert_frame_equal

from mangadap.stack import stack

class SelectTestCase(unittest.TestCase):


    def test_mean_no_nan(self):
        arr = np.arange(4)
        ind = np.array([False, True, True, False])
        desired = 1.5
        actual = stack.mean(arr, ind)
        self.assertEqual(actual, desired)

    def test_mean_drop_nan(self):
        arr = np.array([0, np.nan, 2, 3])
        ind = np.array([False, True, True, True])
        desired = 2.5
        actual = stack.mean(arr, ind)
        self.assertEqual(actual, desired)

    def test_median_no_nan(self):
        arr = np.arange(4)
        ind = np.array([False, True, True, True])
        desired = 2
        actual = stack.median(arr, ind)
        self.assertEqual(actual, desired)

    def test_median_drop_nan(self):
        arr = np.array([0, np.nan, 2, 3])
        ind = np.array([True, True, True, True])
        desired = 2
        actual = stack.median(arr, ind)
        self.assertEqual(actual, desired)

    def test_ivar_wtmean_no_nan(self):
        arr = np.arange(1, 6)
        df = pd.DataFrame(dict(x=arr, y=arr*2))
        err = pd.DataFrame(dict(x=arr/10., y=arr/100.))
        ind = np.array([True, True, False, False, False])
        desired = 1.2
        actual = stack.ivar_wtmean(df.x, err.x, ind)
        self.assertEqual(actual, desired)

    def test_ivar_wtmean_drop_nan(self):
        arr = np.arange(1, 6, dtype=float)
        arr[2] = np.nan
        df = pd.DataFrame(dict(x=arr, y=arr*2))
        err = pd.DataFrame(dict(x=arr/10., y=arr/100.))
        ind = np.array([True, True, True, False, False])
        desired = 1.2
        actual = stack.ivar_wtmean(df.x, err.x, ind)
        self.assertEqual(actual, desired)


if __name__ == '__main__':
    unittest.main()
