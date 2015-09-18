from __future__ import division, print_function, absolute_import

import unittest

import cfg_io

class CfgIOTestCase(unittest.TestCase):

    def test_string_to_float(self):
        d = dict(a='7.5', b='1')
        desired = dict(a=7.5, b=1)
        actual = cfg_io.string_to_float(d)
        self.assertDictEqual(actual, desired)

if __name__ == '__main__':
    unittest.main()