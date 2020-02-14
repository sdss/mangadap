
__version__ = '2.5.2'
__license__ = 'BSD3'
__author__ = 'Kyle B. Westfall'
__maintainer__ = 'Kyle B. Westfall'
__email__ = 'westfall@ucolick.org'
__copyright__ = '(c) 2014-2019, SDSS-IV/MaNGA Pipeline Group'
__credits__ = [ 'Kyle B. Westfall', 'Brett H. Andrews', 'Michele Cappellari',
                'Xihan Ji' ]
# 'Matthew Bershady', 'Francesco Belfiore', 'Kevin Bundy', 'David Law',
# 'Renbin Yan', , 'Daniel Goddard', 'Cheng Li', 'Taniya Parikh', 'Oliver
# Steele', 'Adam Schaefer', 'Shravan Shetty', 'Daniel Thomas', 'Christy
# Tremonti', 'David Wilkinson', 'Zheng Zheng'

import os

def short_warning(message, category, filename, lineno, file=None, line=None):
    """
    Return the format for a short warning message.
    """
    return ' %s: %s\n' % (category.__name__, message)

import warnings
warnings.formatwarning = short_warning

_MANGADRP_VER = 'v2_7_1'
_MANGACORE_VER = 'v1_8_0'

# WARNING: This has to be kept in sync with setup.py
# TODO: Put these in a config file
def default_paths():
    return { 'MANGADRP_VER': _MANGADRP_VER,
             'MANGA_SPECTRO_REDUX': os.path.join(os.environ['HOME'], 'MaNGA', 'redux'),
             'MANGADAP_VER': __version__,
             'MANGA_SPECTRO_ANALYSIS': os.path.join(os.environ['HOME'], 'MaNGA', 'analysis'),
             'MANGACORE_VER': _MANGACORE_VER,
             'MANGACORE_DIR': os.path.join(os.environ['HOME'], 'MaNGA', 'core', _MANGACORE_VER)
           }

def check_environment():
    ev = default_paths()
    for k in ev.keys():
        if k not in os.environ and 'MANGACORE' not in k:
            warnings.warn('{0} environmental variable undefined.  Using: {1}'.format(k,ev[k]))
            os.environ[k] = ev[k]

check_environment()

def dap_source_dir():
    """Return the root path to the DAP source directory."""
    dirlist = os.path.dirname(os.path.abspath(__file__)).split('/')[:-2]
    return os.path.join(os.sep, *dirlist) if dirlist[0] == '' else os.path.join(*dirlist)

os.environ['MANGADAP_DIR'] = dap_source_dir()

