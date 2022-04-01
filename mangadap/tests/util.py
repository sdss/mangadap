import os
import warnings

import pytest

from mangadap.config import defaults

def data_test_file(filename=None):
    root = os.path.join(defaults.dap_data_root(), 'tests')
    return root if filename is None else os.path.join(root, filename)

def remote_data_file(filename=None):
    root = os.path.join(defaults.dap_data_root(), 'remote')
    return root if filename is None else os.path.join(root, filename)

def remote_data_files():
    return ['manga-7815-3702-LINCUBE.fits.gz', 'manga-7815-3702-LINRSS.fits.gz',
            'manga-7815-3702-LOGCUBE.fits.gz', 'manga-7815-3702-LOGRSS.fits.gz']
#    return ['manga-7815-6101-LINCUBE.fits.gz', 'manga-7815-6101-LINRSS.fits.gz',
#            'manga-7815-6101-LOGCUBE.fits.gz', 'manga-7815-6101-LOGRSS.fits.gz']

drp_test_version = 'v3_1_1'
dap_test_version = '3.1.0'
remote_available = all([os.path.isfile(remote_data_file(f)) for f in remote_data_files()])
drpcomplete_available = os.path.isfile(remote_data_file('drpcomplete_{0}.fits'.format(
                                       drp_test_version)))
drpall_available = os.path.isfile(remote_data_file('drpall-{0}.fits'.format(drp_test_version)))

requires_remote = pytest.mark.skipif(not remote_available, reason='Remote data files are missing.')
requires_drpcomplete = pytest.mark.skipif(not drpcomplete_available,
                                          reason='Remote DRPComplete is missing.')
requires_drpall = pytest.mark.skipif(not drpall_available, reason='Remote DRPall is missing.')

if not remote_available:
    warnings.warn('Remote data not available.  Some tests will be skipped.')
if not drpcomplete_available:
    warnings.warn('Remote DRPComplete not available.  Some tests will be skipped.')
if not drpall_available:
    warnings.warn('Remote DRPall not available.  Some tests will be skipped.')

