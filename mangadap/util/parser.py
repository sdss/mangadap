# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""

Provides a set of parsing utility functions.

.. todo::
    - Add function that will parse the default MaNGA fits file name (see
      :func:`mangadap.config.defaults.manga_fits_root`).

----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
from pathlib import Path
from os import environ
import warnings
from IPython import embed
from configparser import ConfigParser, ExtendedInterpolation

def arginp_to_list(inp, evaluate=False, quiet=True):
    """
    Separate a list of comma-separated values in the input string to
    a :obj:`list`.

    Args:
        inp (:obj:`str`, :obj:`list`):
            Input string with a list of comma-separated values
        evaluate (:obj:`bool`, optional):
            Attempt to evaluate the elements in the list using
            `eval`_.
        quiet (:obj:`bool`, optional):
            Suppress terminal output.

    Returns:
        :obj:`list`: The list of the comma-separated values
    """
    out = inp

    # Simply return a None value
    if out is None:
        return out

    # If the input value is a string, convert it to a list of strings
    if isinstance(out, str):
        out = out.replace("[", "")
        out = out.replace("]", "")
        out = out.replace(" ", "")
        out = out.split(',')

    # If the input is still not a list, make it one
    if not isinstance(out, list):
        out = [out]

    # Evaluate each string, if requested
    if evaluate:
        n = len(out)
        for i in range(0,n):
            try:
                tmp = eval(out[i])
            except (NameError, TypeError) as e:
                if not quiet:
                    warnings.warn(f'Could not evaluate value {out[i]}. Skipping.')
            else:
                out[i] = tmp

    return out


# TODO: Replace with ', '.join([str(f) for f in flist])
def list_to_csl_string(flist):
    """
    Convert a list to a comma-separated string; i.e. perform the inverse
    of :func:`arginp_to_list`.

    Args:
        flist (:obj:`list`):
            List to convert to a comma-separated string

    Returns:
        :obj:`str`: String with the values in the input list
        converted to strings separated by commas
    """
    n = len(flist)
    out=str(flist[0])
    for i in range(1,n):
        out += ', '+str(flist[i])
    return out


class DefaultConfig:
    """
    A wrapper for the ConfigParser class that handles None values and
    provides some convenience functions.
    """
    def __init__(self, f=None, interpolate=False):
        # TODO: May not need extended interpolation anymore...
        self.cnfg = ConfigParser(environ, allow_no_value=True,
                                 interpolation=ExtendedInterpolation()) \
                        if interpolate else ConfigParser(allow_no_value=True)
        if f is None:
            return
        self.read(f)

    def __getitem__(self, key):
        return self.cnfg['default'][key]

    def __iter__(self):
        return self.cnfg.options('default').__iter__()

    def read(self, f):
        _f = Path(f).resolve()
        if not _f.exists():
            raise FileNotFoundError(f'Configuration file not found: {_f}')
        self.cnfg.read(f)

    def keyword_specified(self, key):
        return key in self and self[key] is not None

    def all_required(self, keys):
        for k in keys:
            if k not in self or self[k] is None:
                return False
        return True

    def get(self, key, default=None):
        return default if key not in self or self.cnfg['default'][key] is None \
                        else self.cnfg['default'][key]

    def getint(self, key, default=None):
        return default if key not in self or self.cnfg['default'][key] is None \
                        else self.cnfg['default'].getint(key)

    def getbool(self, key, default=None):
        return default if key not in self or self.cnfg['default'][key] is None \
                        else self.cnfg['default'].getboolean(key)

    def getfloat(self, key, default=None):
        return default if key not in self or self.cnfg['default'][key] is None \
                        else self.cnfg['default'].getfloat(key)

    def getlist(self, key, evaluate=False, default=None):
        if key not in self or self.cnfg['default'][key] is None:
            return default
        return [ eval(e.strip()) if evaluate else e.strip() \
                    for e in self.cnfg['default'][key].split(',') ]



