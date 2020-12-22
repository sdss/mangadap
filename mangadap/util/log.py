# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Logging routines.

----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import logging
import warnings

from astropy.wcs import FITSFixedWarning

def init_DAP_logging(log, simple_warnings=True, append=False, keep_fits_warnings=False):
    """
    Initialize the logging preferences for the DAP.

    Args:
        log (:obj:`str`):
            File with log output.
        simple_warnings (:obj:`bool`, optional):
            Shorten warning messages.
        append (:obj:`bool`, optional):
            Append new log messages to existing log file; if False, file
            is overwritten.
        keep_fits_warning (:obj:`bool`, optional):
            Flag to not ignore FITSFixedWarning messages.

    .. todo::
        Use a file with the logging configuration.  See:
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
    if not keep_fits_warnings:
        warnings.simplefilter('ignore', FITSFixedWarning)
    warnings.filterwarnings('ignore', 'unclosed file')
    logging.captureWarnings(True)
    if simple_warnings:
        warnings.formatwarning = short_warning
#    import numpy
#    warnings.simplefilter('error', numpy.VisibleDeprecationWarning)
#    warnings.simplefilter('error', SyntaxWarning)

    # Add file handler if wanted
    if log is not None:
        logfile = logging.FileHandler(log, mode=('a' if append else 'w'))
        logfile.setFormatter(
                logging.Formatter('%(asctime)s %(name)-10s %(levelname)-8s ::  %(message)s'))
        logfile.setLevel(logging.DEBUG)
        logging.getLogger('').addHandler(logfile)


def module_logging(name, verbose):
    """
    Return a number of `logging.Logger`_ objects, one for each verbosity
    level.
    """
    if verbose == 0:
        return None
    loggers = []
    for i in range(verbose):
        loggers += [ logging.getLogger('{0} {1}'.format(name,i+1)) ]
    return loggers
    
    
def log_output(loggers, v, lvl, message):
    """
    Write message to the log.

    Args:
        loggers (`logging.Logger`_):  Objects collecting the log for a
            given run of the DAP.
        v (int): Verbosity level of message.
        lvl (`logging.level`_):  Logging level
        message (str): Message to log.
    """
    if loggers is None:
        print(message)
        return

    if len(loggers) < v:
        return

    loggers[v-1].log(lvl,message)


def short_warning(message, category, filename, lineno, file=None, line=None):
    """
    Return the format for a short warning message.
    """
    return ' %s: %s' % (category.__name__, message)
    

