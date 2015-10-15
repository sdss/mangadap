from __future__ import division, print_function, absolute_import

import unittest

from mangadap.plot import cfg_io

class CfgIOTestCase(unittest.TestCase):

    def test_string_to_float(self):
        d = dict(a='7.5', b='1')
        desired = dict(a=7.5, b=1)
        actual = cfg_io.string_to_float(d)
        self.assertDictEqual(actual, desired)

    def test_convert_to_number_int(self):
        desired = 7
        actual = cfg_io.convert_to_number('7')
        self.assertEqual(actual, desired)

    def test_convert_to_number_float(self):
        desired = 7.1
        actual = cfg_io.convert_to_number('7.1')
        self.assertEqual(actual, desired)

    def test_convert_to_number_int_overrides_float(self):
        desired = 7
        actual = cfg_io.convert_to_number('7.0')
        self.assertEqual(actual, desired)

    def test_convert_to_number_list_int(self):
        desired = [7, 8]
        actual = cfg_io.convert_to_number_list(['7', '8'])
        self.assertEqual(actual, desired)

    def test_convert_to_number_list_float(self):
        desired = [7.1, 8.3]
        actual = cfg_io.convert_to_number_list(['7.1', '8.3'])
        self.assertEqual(actual, desired)

    def test_convert_to_number_list_float_overrides_int(self):
        desired = [7.1, 8.0]
        actual = cfg_io.convert_to_number_list(['7.1', '8.0'])
        self.assertEqual(actual, desired)

    def test_convert_dtype_number_list(self):
        desired = [7.1, 8.0]
        actual = cfg_io.convert_dtype(['7.1', '8.0'])
        self.assertEqual(actual, desired)

if __name__ == '__main__':
    unittest.main()