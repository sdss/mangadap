import os

import pytest

from mangadap.config import defaults

def data_file(filename):
    return os.path.join(defaults.dap_data_root(), 'tests', filename)

def cube_test_file():
    plt = 7815
    ifu = 3702
    return os.path.join(defaults.drp_directory_path(plt),
                        '{0}.fits.gz'.format(defaults.manga_fits_root(plt, ifu, 'LOGCUBE')))
    
cube_test_only = pytest.mark.skipif(not os.path.isfile(cube_test_file()),
                                    reason='requires example cube file')

