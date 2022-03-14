
__version__ = '3.1.3dev'
__license__ = 'BSD3'
__author__ = 'Kyle B. Westfall'
__maintainer__ = 'Kyle B. Westfall'
__email__ = 'westfall@ucolick.org'
__copyright__ = '(c) 2014-2020, SDSS-IV/MaNGA Pipeline Group'
__credits__ = ['Kyle B. Westfall', 'Brett H. Andrews', 'Jorge Barrera-Ballesteros',
               'Francesco Belfiore', 'Matthew A. Bershady', 'Kevin Bundy', 'Joel R. Brownstein',
               'Michele Cappellari', 'Brian Cherinka', 'Lodovico Coccato', 'Niv Drory',
               'Cheng Du', 'Eric Emsellem', 'Daniel Goddard', 'Xihan Ji', 'David R. Law',
               'Niu Li', 'Claudia Maraston', 'Karen Masters', 'Héctor Javier Ibarra Medel',
               'Taniya Parikh', 'Sebastián F. Sánchez', 'José R. Sánchez-Gallego',
               'Adam Schaefer', 'Shravan Shetty', 'Daniel Thomas', 'Christy A. Tremonti',
               'Anne-Marie Weijmans', 'Renbin Yan', 'Meng Yang', 'Zheng Zheng', 'Shuang Zhou']

def short_warning(message, category, filename, lineno, file=None, line=None):
    """
    Return the format for a short warning message.
    """
    return ' %s: %s\n' % (category.__name__, message)

import warnings
warnings.formatwarning = short_warning

_MANGADRP_VER = 'v3_1_1'
_MANGACORE_VER = 'v1_9_1'

# WARNING: This has to be kept in sync with setup.py
# TODO: Put these in a config file
import os
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
    dirlist = os.path.dirname(os.path.abspath(__file__)).split('/')[:-1]
    return os.path.join(os.sep, *dirlist) if dirlist[0] == '' else os.path.join(*dirlist)

os.environ['MANGADAP_DIR'] = dap_source_dir()

