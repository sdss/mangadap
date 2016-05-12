# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Logging routines.

*License*:
    Copyright (c) 2015, SDSS-IV/MaNGA Pipeline Group
        Licensed under BSD 3-clause license - see LICENSE.rst

*Source location*:
    $MANGADAP_DIR/python/mangadap/util/log.py

*Imports and python version compliance*:
    ::

        from __future__ import division
        from __future__ import print_function
        from __future__ import absolute_import
        from __future__ import unicode_literals

        import sys
        if sys.version > '3':
            long = int

        import logging

*Revision history*:
    | **21 Mar 2016**: Original implementation by K. Westfall (KBW)

"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import sys
if sys.version > '3':
    long = int

import logging
import warnings

__author__ = 'Kyle Westfall'

def init_DAP_logging(log, simple_warnings=True):
    """
    .. todo::
        - Use a file with the logging configuration.  See:
        https://docs.python.org/3.5/howto/logging.html#handler-basic
    """

    # Remove existing root logging
    if len(logging.root.handlers) > 0:
        for handler in logging.root.handlers:
            logging.root.removeHandler(handler)

    # Init global logging
    logging.basicConfig(level=logging.INFO,
                        format='%(levelname)-8s ::  %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')

    # Set warnings format
    warnings.simplefilter('always')
    warnings.filterwarnings('ignore', 'unclosed file')
    logging.captureWarnings(True)
    if simple_warnings:
        warnings.formatwarning = short_warning

    # Add file handler if wanted
    if log is not None:
        logfile = logging.FileHandler(log, mode='w')
        logfile.setFormatter(
                logging.Formatter('%(asctime)s %(name)-10s %(levelname)-8s ::  %(message)s'))
        logfile.setLevel(logging.DEBUG)
        logging.getLogger('').addHandler(logfile)


def module_logging(name, verbose):
    if verbose == 0:
        return None
    loggers = []
    for i in range(verbose):
        loggers += [ logging.getLogger('{0} {1}'.format(name,i+1)) ]
    return loggers
    
    
def log_output(loggers, v, lvl, message):
    if loggers is None or len(loggers) < v:
        return
    loggers[v-1].log(lvl,message)


def short_warning(message, category, filename, lineno, file=None, line=None):
    return ' %s: %s' % (category.__name__, message)
    

