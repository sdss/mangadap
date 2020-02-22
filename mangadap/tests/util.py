import os
import warnings

import pytest

from mangadap.config import defaults

def test_data_file(filename=None):
    root = os.path.join(defaults.dap_data_root(), 'tests')
    return root if filename is None else os.path.join(root, filename)

def remote_data_file(filename=None):
    root = os.path.join(defaults.dap_data_root(), 'remote')
    return root if filename is None else os.path.join(root, filename)

def remote_data_files():
    return ['manga-7815-3702-LINCUBE.fits.gz', 'manga-7815-3702-LINRSS.fits.gz',
            'manga-7815-3702-LOGCUBE.fits.gz', 'manga-7815-3702-LOGRSS.fits.gz']

remote_available = all([os.path.isfile(remote_data_file(f)) for f in remote_data_files()])

requires_remote = pytest.mark.skipif(not remote_available, reason='Remote data files are missing.')

if not remote_available:
    warnings.warn('Remote data not available.  Some tests will be skipped.')

