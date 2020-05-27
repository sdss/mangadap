
import pytest

import os
import mangadap

# TODO: I don't think this is ever used...
def pytest_report_header(config):
    drp_root = os.path.join(os.environ['MANGA_SPECTRO_REDUX'], os.environ['MANGADRP_VER'])
    if not os.path.isdir(drp_root):
        return ['Root path with DRP files does not exist.  Some tests not executed.',
                'DRP path: {0}'.format(drp_root)]
    return 'Tests will use DRP files in {0}'.format(drp_root)

