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
import os
from configparser import ConfigParser, ExtendedInterpolation

from .exception_tools import print_frame

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
            except (NameError, TypeError):
                if not quiet:
                    print_frame('NameError/TypeError')
                    print('Could not evaluate value {0}. Skipping.'.format(out[i]))
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


def parse_drp_file_name(name):
    """
    Parse the name of a DRP file to provide the plate, ifudesign, and
    mode.

    Args:
        name (str): Name of the DRP file.

    Returns:
        int, int, str: The plate, ifudesign, and mode ('RSS' or 'CUBE')
        of the DRP file pulled from the file name.

    Raises:
        TypeError: Raised if *name* is not a string.
        ValueError: Raised if if the file name does not look like a DRP
            file because it does not include 'manga-', '-LOG', or
            '.fits.gz'.
    """

    if (type(name) != str):
        raise TypeError("File name must be a string.")
    if (str.find(name, 'manga-') == -1 or str.find(name, '-LOG') == -1 or 
        str.find(name, '.fits.gz') == -1):
        raise ValueError("String does not look like a DRP fits file name.")

    plate_start = str.find(name, '-')+1
    plate_end = str.find(name,'-',plate_start)
    plate = long(name[plate_start:plate_end])

    ifudesign_end = str.find(name,'-',plate_end+1)
    ifudesign = long(name[plate_end+1:ifudesign_end])

    mode_start = str.find(name,'LOG')+3
    mode = name[mode_start:str.find(name,'.fits.gz')]

    return plate, ifudesign, mode


def parse_dap_file_name(name):
    """
    Parse the name of a DAP file and return the plate, ifudesign, mode,
    binning type, and iteration number.

    Args:
        name (str): Name of the DAP file.

    Returns:
        int, int, str, str, int : The plate, ifudesign, mode ('RSS' or
        'CUBE'), bin type, and iteration number of the DAP file, pulled
        from the name of the file.

    Raises:
        TypeError: Raised if *name* is not a string.
        ValueError: Raised if if the file name does not look like a DRP
            file because it does not include 'manga-', '-BIN', or
            '.fits'.
    """
    if (type(name) != str):
        raise TypeError("File name must be a string.")
    if (str.find(name, 'manga-') == -1 or str.find(name, '-LOG') == -1
        or str.find(name, 'BIN-') == -1 or str.find(name, '.fits') == -1):
        raise ValueError("String does not look like a DAP fits file name.")

    plate_start = str.find(name, '-')+1
    plate_end = str.find(name,'-',plate_start)
    plate = int(name[plate_start:plate_end])

    ifudesign_end = str.find(name,'-',plate_end+1)
    ifudesign = int(name[plate_end+1:ifudesign_end])

    mode_start = str.find(name,'LOG')+3
    mode_end = str.find(name,'_BIN-')
    mode = name[mode_start:mode_end]

    bintype_end = str.find(name,'-',mode_end+5)
    bintype = name[mode_end+5:bintype_end]

    niter_end = str.find(name,'.fits',bintype_end)
    niter = int(name[bintype_end+1:niter_end])

    return plate, ifudesign, mode, bintype, niter

    
class DefaultConfig:
    """
    A wrapper for the ConfigParser class that handles None values and
    provides some convenience functions.
    """
    def __init__(self, f=None, interpolate=False):
        self.cnfg = ConfigParser(os.environ, allow_no_value=True,
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
        if not os.path.isfile(f):
            raise FileNotFoundError('Configuration file not found: {0}'.format(f))
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



