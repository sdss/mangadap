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
    """

    caller_frame = inspect.stack()[1]

    frame = caller_frame[0]
    info = inspect.getframeinfo(frame)
    print('{0}: FILE={1}, FUNCTION={2}, LINE={3}\n'.format(prefix, info.filename, info.function,
          info.lineno))

