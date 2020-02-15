
import pytest
import os

from IPython import embed

from mangadap.par.absorptionindexdb import AbsorptionIndexDB

import warnings
warnings.simplefilter("ignore", UserWarning)
warnings.simplefilter("ignore", RuntimeWarning)

def test_read():
    dbs = AbsorptionIndexDB.available_databases()
    assert len(dbs) > 0, 'No absorption-line-index databases available'
    for key in dbs.keys():
        absdb = AbsorptionIndexDB.from_key(key)

def test_lick():
    absdb = AbsorptionIndexDB.from_key('LICKINDX')
    assert len(absdb) == 21, 'Incorrect number of Lick indices'
    assert 'Hb' in absdb['name'], 'Does not contain Hb in list'
