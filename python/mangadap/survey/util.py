# List of provided functions:
#
#       product_version - system call to $product_version, parses output
#       module_version  - determines version of loaded module

# Force Python 3 behavior
from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

# Check long/int definition
import sys
if sys.version > '3':
    long = int

# Additional imports
import subprocess
from os import environ
from mangadap.util.exception_tools import print_frame

def product_version(product='mangadap', simple=False):
    """
    Gets the version for the SDSS-III or SDSS-IV product.  The
    default product is mangadap.
    """

    # Expects to find an executable called {$product}_version that
    # reports the SDSS-III/SDSS-IV product version.  If simple=True,
    # only the first element of the reported version is set to
    # 'version'
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
    Return the module version for the specified product.  The
    default product is mangadap.
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





