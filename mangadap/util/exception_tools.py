# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""

Provides a set of tools that the MaNGA DAP uses to handle/report raised
exceptions.

----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import inspect
import os

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
        name (:obj:`str`):
            Name of a required environmental variable

    Raises:
        EnvironmentError: Raised if *name* is not defined.
    """
    if name not in os.environ:
        raise EnvironmentError('{0} not defined in current environment!')





