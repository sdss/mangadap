from __future__ import division, print_function, absolute_import

import unittest
import numpy as np

import matplotlib.cm as cm

from mangadap.plot import plotdap
from mangadap.plot import util

class PlotdapTestCase(unittest.TestCase):

    def test_set_cmaps_all_given_as_matplotlib_cm(self):
        desired = [cm.RdBu, cm.jet]
        actual = plotdap.set_cmaps(cmaps=['RdBu', 'jet'], n_plots=2)
        self.assertListEqual(actual, desired)

    def test_set_cmaps_all_given_as_linear_Lab(self):
        linear_Lab, linear_Lab_r = util.linear_Lab()
        cm_desired = [linear_Lab, linear_Lab]
        desired = [it._segmentdata for it in cm_desired]
        cm_actual = plotdap.set_cmaps(cmaps=['linear_Lab', 'linear_Lab'],
                                      n_plots=2)
        actual = [it._segmentdata for it in cm_actual]
        self.assertListEqual(actual, desired)

    def test_set_cmaps_some_given_as_linear_Lab(self):
        linear_Lab, linear_Lab_r = util.linear_Lab()
        cm_desired = [cm.RdBu, linear_Lab]
        desired = [it._segmentdata for it in cm_desired]
        cm_actual = plotdap.set_cmaps(cmaps=['RdBu', 'linear_Lab'], n_plots=2)
        actual = [it._segmentdata for it in cm_actual]
        self.assertListEqual(actual, desired)

    def test_set_cmaps_all_given_as_linear_Lab_reversed(self):
        linear_Lab, linear_Lab_r = util.linear_Lab()
        cm_desired = [linear_Lab_r, linear_Lab_r]
        desired = [it._segmentdata for it in cm_desired]
        cm_actual = plotdap.set_cmaps(cmaps=['linear_Lab_r', 'linear_Lab_r'],
                                      n_plots=2)
        actual = [it._segmentdata for it in cm_actual]
        self.assertListEqual(actual, desired)

    def test_set_cmaps_some_given_as_linear_Lab_reversed(self):
        linear_Lab, linear_Lab_r = util.linear_Lab()
        cm_desired = [cm.RdBu, linear_Lab_r]
        desired = [it._segmentdata for it in cm_desired]
        cm_actual = plotdap.set_cmaps(cmaps=['RdBu', 'linear_Lab_r'], n_plots=2)
        actual = [it._segmentdata for it in cm_actual]
        self.assertListEqual(actual, desired)

    def test_set_cmaps_one_given_as_matplotlib_cm(self):
        desired = [cm.RdBu, cm.RdBu]
        actual = plotdap.set_cmaps(cmaps=['RdBu'], n_plots=2)
        self.assertListEqual(actual, desired)

    def test_set_cmaps_one_given_as_linear_Lab(self):
        linear_Lab, linear_Lab_r = util.linear_Lab()
        cm_desired = [linear_Lab, linear_Lab]
        desired = [it._segmentdata for it in cm_desired]
        cm_actual = plotdap.set_cmaps(cmaps=['linear_Lab'], n_plots=2)
        actual = [it._segmentdata for it in cm_actual]
        self.assertListEqual(actual, desired)

    def test_set_cmaps_one_given_as_linear_Lab_reversed(self):
        linear_Lab, linear_Lab_r = util.linear_Lab()
        cm_desired = [linear_Lab_r, linear_Lab_r]
        desired = [it._segmentdata for it in cm_desired]
        cm_actual = plotdap.set_cmaps(cmaps=['linear_Lab_r'], n_plots=2)
        actual = [it._segmentdata for it in cm_actual]
        self.assertListEqual(actual, desired)

    def test_set_cmaps_None(self):
        linear_Lab, linear_Lab_r = util.linear_Lab()
        cm_desired = [linear_Lab, linear_Lab]
        desired = [it._segmentdata for it in cm_desired]
        cm_actual = plotdap.set_cmaps(cmaps=None, n_plots=2)
        actual = [it._segmentdata for it in cm_actual]
        self.assertListEqual(actual, desired)


if __name__ == '__main__':
    unittest.main()


