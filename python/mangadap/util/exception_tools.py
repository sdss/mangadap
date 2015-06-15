"""

Provides a set of tools that the MaNGA DAP uses to handle/report raised
exceptions.

*Source location*:
    $MANGADAP_DIR/python/mangadap/util/exception_tools.py

*Python2/3 compliance*::

    from __future__ import division
    from __future__ import print_function
    from __future__ import absolute_import
    
    import sys
    if sys.version > '3':
        long = int

*Imports*::

    import inspect

*Revision history*:
    | **2015**: Original implementation by K. Westfall (KBW)
    | **20 May 2015**: (KBW) Documentation and Sphinx tests
"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import sys
if sys.version > '3':
    long = int

import inspect

__author__ = 'Kyle Westfall'

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


