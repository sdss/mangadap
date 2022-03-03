import os
from ..version import version

def get_environ():

    # Default versions
    _MANGADRP_VER = 'v3_1_1'
    _MANGACORE_VER = 'v1_9_1'

    # First set the default environment
    environ = {'MANGADRP_VER': _MANGADRP_VER,
               'MANGA_SPECTRO_REDUX': os.path.join(os.getcwd(), 'MaNGA', 'redux'),
               'MANGADAP_VER': version,
               'MANGA_SPECTRO_ANALYSIS': os.path.join(os.getcwd(), 'MaNGA', 'analysis'),
               'MANGACORE_VER': _MANGACORE_VER,
               'MANGACORE_DIR': os.path.join(os.getcwd(), 'MaNGA', 'core', _MANGACORE_VER)}

    # Replace with variables from working environment if they exist
    for key in environ.keys():
        if key in os.environ:
            environ[key] = os.environ[key]

    return environ

manga_environ = get_environ()

#from .manga import MaNGAAnalysisPlan

