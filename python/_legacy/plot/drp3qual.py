"""Read plateifu and drp3qual from DRPall file for DRP QA webpage."""
from __future__ import division, print_function, absolute_import

import os
from os.path import join

import util

home = os.path.expanduser('~')
path_cfg = join(os.getenv('MANGADAP_DIR'), 'python', 'mangadap', 'plot',
                'config', 'sdss_paths.ini')
drpall = util.read_drpall(paths_cfg)
plateifu = drpall['plateifu']
drp3qual = drpall['drp3qual']
