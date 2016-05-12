# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Provides a set of utility function that should **only be used for the
survey-level execution of the DAP**.

*License*:
    Copyright (c) 2015, SDSS-IV/MaNGA Pipeline Group
        Licensed under BSD 3-clause license - see LICENSE.rst

*Source location*:
    $MANGADAP_DIR/python/mangadap/survey/util.py

*Imports and python version compliance*:
    ::

        from __future__ import division
        from __future__ import print_function
        from __future__ import absolute_import
        from __future__ import unicode_literals

        import sys
        if sys.version > '3':
            long = int

        import subprocess
        from os import environ
        from ..util.exception_tools import print_frame

*Revision history*:
    | **?? ??? 2014**: Original implementation by K. Westfall (KBW).
        Functions pulled (unedited) from file provided by J. Brownstein.
    | **27 Aug 2015**: (KBW) Sphinx documentation.
"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import sys
if sys.version > '3':
    long = int

import subprocess
from os import environ
from ..util.exception_tools import print_frame

def product_version(product='mangadap', simple=False):
    """
    Get currently loaded SDSS product version by parsing the output of
    the system call to {$product}_version.

    Args:
        product (str): (Optional) The name of the product
        simple (bool): (Optional) If True, only the first element of the
            reported version is returned.

    Returns:
        str : Version identifier
    """
    try:
        version = subprocess.check_output('%s_version' % product, shell=True)
        if type(version) is bytes:
            version = version.decode('utf-8')
        version = version.split(' ')[0].rstrip('\n') if simple else version.rstrip('\n')
    except Exception as e:
        print_frame('Exception')
        print(e)
        version = None

    return version


def module_version(product='mangadap'):
    """
    Return the version of the specified product among the currently
    loaded modules.

    Args:
        product (str): (Optional) The name of the product

    Returns:
        str : Version identifier
    """
    try:
        modules = environ['LOADEDMODULES']
    except:
        print_frame('Exception')
        modules = None
        return None
    # TODO: Re-raise the exception?
  
    # Parse the loaded version(s) of product
    versions = [module.split('/')[1] for module in modules.split(':')
                    if module.split('/')[0]==product]

    # If there is more than one version or no versions return None
    if len(versions) != 1:
        if len(versions) > 1:
            print('Multiple versions found for module {0}'.format(product))
        else:
            print('Module {0} is not loaded'.format(product))
        return None

    # Return the version
    return versions[0]





