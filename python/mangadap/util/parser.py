"""

Provides a set of parsing utility functions.

*Source location*:
    $MANGADAP_DIR/python/mangadap/util/parser.py

*Python2/3 compliance*::

    from __future__ import division
    from __future__ import print_function
    from __future__ import absolute_import
    
    import sys
    if sys.version > '3':
        long = int

*Imports*::

    from mangadap.util.exception_tools import print_frame

*Revision history*:
    | **2015**: Original implementation by K. Westfall (KBW)
    | **20 May 2015**: (KBW) Documentation and Sphinx tests
    | **04 Jun 2015**: (KBW) Added :func:`parse_drp_file_name`,
        :func:`parse_dap_file_name`
"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import sys
if sys.version > '3':
    long = int

from mangadap.util.exception_tools import print_frame

def arginp_to_list(inp, evaluate=False, quiet=True):
    """
    Separate a list of comma-separated values in the input string to a
    *list* object.

    Args:
        inp (str or list): Input string with a list of comma-separated values
        evaluate (bool): (Optional) Attempt to evaluate the elements in
            the list using :func:`eval`.
        quiet (bool): (Optional) Suppress terminal output

    Returns:
        list : The list of the comma-separated values
    """
    out = inp

#   print('INP type: {0}'.format(type(inp)))

    # Simply return a None value
    if out is None:
        return out

    # If the input value is a string, convert it to a list of strings
    if type(out) == str:
        out = out.replace("[", "")
        out = out.replace("]", "")
        out = out.replace(" ", "")
        out = out.split(',')

    # If the input is still not a list, make it one
    if type(out) != list:
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

#   print('OUT type: {0}'.format(type(out)))
    return out


#def funcarg_to_list(*arg):
##   print(len(arg))
#    return arg[1:]


def list_to_csl_string(flist):
    """
    Convert a list to a comma-separated string; i.e. perform the inverse
    of :func:`arginp_to_list`.

    Args:
        flist (list): List to convert to a comma-separated string

    Returns:
        str : String with the values in the input list converted to
            string, using :func:`str`, separated by commas
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
        int, int, str : The plate, ifudesign, and mode ('RSS' or 'CUBE')
        of the DRP file pulled from the file name.

    Raises:
        Exception: Raised if *name* is not a string or if the file name
            does not look like a DRP file because it does not include
            'manga-', '-LOG', or '.fits.gz'.

    """

    if (type(name) != str):
        raise Exception("File name must be a string.")
    if (str.find(name, 'manga-') == -1 or str.find(name, '-LOG') == -1 or 
        str.find(name, '.fits.gz') == -1):
        raise Exception("String does not look like a DRP fits file name.")

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
            'CUBE'), bin type, and iteration number of the DAP file,
            pulled from the name of the file.

    Raises:
        Exception: Raised if *name* is not a string or if the file name
            does not look like a DAP file because it does not include
            'manga-', '-LOG', 'BIN-', or '.fits'.

    """
    if (type(name) != str):
        raise Exception("File name must be a string.")
    if (str.find(name, 'manga-') == -1 or str.find(name, '-LOG') == -1
        or str.find(name, 'BIN-') == -1 or str.find(name, '.fits') == -1):
        raise Exception("String does not look like a DAP fits file name.")

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

    

