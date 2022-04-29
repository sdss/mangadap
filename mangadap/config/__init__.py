
def get_manga_environ():

    import os
    from ..version import version

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

    environ['sdss-access'] = os.environ['SDSS_ACCESS_DIR'].split('/')[-1] \
                                if 'SDSS_ACCESS_DIR' in os.environ else None
    environ['idlutils'] = os.environ['IDLUTILS_DIR'].split('/')[-1] \
                            if 'IDLUTILS_DIR' in os.environ else None

    return environ

manga_environ = get_manga_environ()

def get_python_versions():

    versions = {}
    import sys
    versions['python'] = '.'.join([str(v) for v in sys.version_info[:3]])
    import numpy
    versions['numpy'] = numpy.__version__
    import scipy
    versions['scipy'] = scipy.__version__
    import matplotlib
    versions['matplotlib'] = matplotlib.__version__
    import astropy
    versions['astropy'] = astropy.__version__
    import pydl
    versions['pydl'] = pydl.__version__
    import ppxf
    versions['ppxf'] = ppxf.__version__
    import vorbin
    versions['vorbin'] = vorbin.__version__
    from ..version import version
    versions['mangadap'] = version
    return versions

python_versions = get_python_versions()

