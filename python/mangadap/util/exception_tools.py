# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""

Provides a set of tools that the MaNGA DAP uses to handle/report raised
exceptions.

*License*:
    Copyright (c) 2015, SDSS-IV/MaNGA Pipeline Group
        Licensed under BSD 3-clause license - see LICENSE.rst

*Source location*:
    $MANGADAP_DIR/python/mangadap/util/exception_tools.py

*Imports and python version compliance*:
    ::

        from __future__ import division
        from __future__ import print_function
        from __future__ import absolute_import
        from __future__ import unicode_literals

        import sys
        if sys.version > '3':
            long = int

        import inspect
        from os import environ

*Revision history*:
    | **2015**: Original implementation by K. Westfall (KBW)
    | **20 May 2015**: (KBW) Documentation and Sphinx tests
    | **03 Feb 2016**: (KBW) Added :func:`check_environment_variable`
"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import sys
if sys.version > '3':
    long = int

import inspect
from os import environ

def print_frame(prefix):
    """
    Print the frame stack.

    Args:
        prefix (str): Prefix (typically the exception type) for the
            printed statement
    """
    caller_frame = inspect.stack()[1]
    frame = caller_frame[0]
    info = inspect.getframeinfo(frame)
    print('{0}: FILE={1}, FUNCTION={2}, LINE={3}\n'.format(prefix, info.filename, info.function,
          info.lineno))


def check_environment_variable(name):
    """
    Check for the existence of an environment variable.

    Args:
        name (str): Name of a required environmental variable

    Raises:
        EnvironmentError: Raised if *name* is not defined.
    """
    if name not in environ:
        raise EnvironmentError('{0} not defined in current environment!')





