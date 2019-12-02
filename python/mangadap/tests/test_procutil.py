
import pytest
import os

import warnings
warnings.simplefilter("ignore", UserWarning)
warnings.simplefilter("ignore", RuntimeWarning)

from mangadap.proc import util

def test_database_key():
    assert util.get_database_key('junk') == 'JUNK'
    assert util.get_database_key('test.par') == 'TEST'
    assert util.get_database_key('/path/to/test.par') == 'TEST' 

def test_select_database():
    d = os.path.join(os.environ['MANGADAP_DIR'], 'data', 'emission_lines')
    assert os.path.split(util.select_database('ELPMPL9', d))[1] == 'elpmpl9.par'

