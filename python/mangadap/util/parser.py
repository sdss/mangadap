# List of provided functions:
#
#       arginp_to_list - Convert a input argument to a list of values

# Force Python 3 behavior
from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

# Check long/int definition
import sys
if sys.version > '3':
    long = int

# Additional imports
from mangadap.util.exception_tools import print_frame


def arginp_to_list(inp, evaluate=False, quiet=True):
    """
    Convert an input argument to a list of separated values.

    If evaluate is True, evaluate each element in the list before
    returning it.
    """

    out = inp

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

    return out


def funcarg_to_list(*arg):
    print(len(arg))
    return arg[1:]


